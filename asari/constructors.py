'''
Classes of MassGrid and CompositeMap.
'''

import os
import csv
from functools import lru_cache

import pandas as pd
import numpy as np
from scipy import interpolate
from scipy.ndimage import maximum_filter1d
from mass2chem.search import find_mzdiff_pairs_from_masstracks

from .mass_functions import (flatten_tuplelist, 
                            landmark_guided_mapping, 
                            calculate_selectivity)
from .chromatograms import (nn_cluster_by_mz_seeds,
                            rt_lowess_calibration_debug,
                            rt_lowess_calibration, 
                            remap_intensity_track,
                            __hacked_lowess__)
from .peaks import (quick_detect_unique_elution_peak, 
                    batch_deep_detect_elution_peaks, 
                    get_gaussian_peakarea_on_intensity_list)
from .samples import SimpleSample
from .utils import bulk_process


class MassGrid:
    '''
    MassGrid is the concept for m/z correspondence in asari.
    This shares similarity to FeatureMap in OpenMS, but the correspondence 
    in asari takes adavantage of high m/z resolution first before feature detection.
    '''
    def __init__(self, cmap, experiment):
        '''
        Initiating MassGrid by linking to CompositeMap and ext_Experiment instances.

        Parameters
        ----------
        cmap : CompositeMap intance, optional, default: None
            CompositeMap instance.
        experiment : ext_Experiment instance, optional, default: None
            ext_Experiment instance.

        
        Notes
        -----
        cmap and experiment are currently optional, but cmap appears to be required.
        should refactor this if possible.
        '''

        self.experiment = experiment
        self.CMAP = cmap
        self.reference_sample_instance = self.CMAP.reference_sample_instance
        self.max_ref_rtime = self.CMAP.max_ref_rtime
        self.list_sample_names = self.CMAP.list_sample_names
        self._number_of_samples_ = self.CMAP._number_of_samples_
        
    def build_grid_sample_wise(self):
        '''
        Align one sample at a time to reference m/z grid, based on their anchor m/z tracks.
        One of the two methods to build the grid.
        This is better for reliable assembly of small number of samples.

        See also
        --------
        build_grid_by_centroiding
        '''
        self._initiate_mass_grid()
        sample_ids = self.experiment.valid_sample_ids
        sample_ids.pop(self.experiment.reference_sample_id)
        for sid in sample_ids:
            SM = SimpleSample(self.experiment.sample_registry[sid],
                experiment=self.experiment, 
                database_mode=self.experiment.database_mode, 
                mode=self.experiment.mode)
            self.add_sample(SM)

    def build_grid_by_centroiding(self):
        '''
        Assemble mass grid by grouping m/z values to centroids.
        Each centroid can have no more than one mass track per sample.
        One of the two methods to build the grid.
        This is more efficient for large number of samples.

        See also
        --------
        build_grid_sample_wise
        '''
        all = []
        for ii in range(self._number_of_samples_):
            sid = self.experiment.valid_sample_ids[ii]
            for jj in self.experiment.sample_registry[sid]['track_mzs']:
                all.append(
                    (jj[0], jj[1], ii)      # m/z, masstrack_id, sample_index
                )
        all.sort()

        all_bins = self.bin_track_mzs(all, self.experiment.reference_sample_id)

        self.MassGrid = pd.DataFrame(
            np.full((len(all_bins), self._number_of_samples_), None),
            columns=self.CMAP.list_sample_names,
        )
        for jj in range(len(all_bins)):
            for (mz, track_id, sample_ii) in all_bins[jj][1]:
                self.MassGrid.iloc[jj, sample_ii] = track_id
        
        mz_list = [x[0] for x in all_bins]
        self.MassGrid.insert(0, "mz", mz_list)          # add extra col for mz

        # determine anchor and landmark mass tracks
        list_mass_tracks = []
        for ii in range(len(mz_list)):
            list_mass_tracks.append({'id_number': ii, 'mz': mz_list[ii],})
        self.anchor_mz_pairs = find_mzdiff_pairs_from_masstracks(
            list_mass_tracks, mz_tolerance_ppm=self.experiment.parameters['mz_tolerance_ppm']
            )
        self._mz_landmarks_ = flatten_tuplelist(self.anchor_mz_pairs)

        # make sample instances
        self.reference_sample_instance.rt_cal_dict = \
              self.reference_sample_instance.reverse_rt_cal_dict = {}
        self.experiment.all_samples.append(self.reference_sample_instance)
        for sid in self.experiment.valid_sample_ids:
            if sid != self.experiment.reference_sample_id:
                SM = SimpleSample(self.experiment.sample_registry[sid],
                    experiment=self.experiment, database_mode=self.experiment.database_mode, 
                    mode=self.experiment.mode
                    )
                self.experiment.all_samples.append(SM)


    def _initiate_mass_grid(self):
        '''
        Initiate self.MassGrid as pandas DataFrame.
        The reference sample is used to populate first m/z column.
        This sets 1st instance into self.experiment.all_samples.

        Updates
        -------
        self._mz_landmarks_ : 
            landmark m/z values that match to 13C/12C pattern
        self.MassGrid : 
            DataFrame with reference sample as first entry
        self.experiment.all_samples : 
            adding 1st sample (reference)
        '''
        reference_sample = self.reference_sample_instance
        reference_sample.rt_cal_dict = reference_sample.reverse_rt_cal_dict = {}
        ref_list_mass_tracks = reference_sample.list_mass_tracks

        self._mz_landmarks_ = reference_sample._mz_landmarks_
        reference_mzlist = [ x['mz'] for x in ref_list_mass_tracks ]
        # setting up DataFrame for MassGrid
        # not forcing dtype on DataFrame, to avoid unreported errors; convert to int when using MassGrid
        self.MassGrid = pd.DataFrame(
            np.full((len(reference_mzlist), 1+self._number_of_samples_), None),
            columns=['mz'] + self.list_sample_names,
        )
        # Add ref mz as a column to MassGrid; ref mzlist will be dynamic updated in MassGrid["mz"]
        self.MassGrid['mz'] = reference_mzlist
        self.MassGrid[ reference_sample.name ] = [ x['id_number'] for x in ref_list_mass_tracks ]
        self.experiment.all_samples.append(reference_sample)


    def add_sample(self, sample):
        '''
        This adds a sample to MassGrid, including the m/z alignment of the sample against the 
        existing reference m/z values in the MassGrid.

        Parameters
        ----------
        sample : SimpleSample instance
            instance of SimpleSample class.
        database_cursor : cursor object
            Not used now.

        Updates
        -------
        self._mz_landmarks_ :   
            landmark m/z values that match to 13C/12C and Na/H patterns
        self.MassGrid : 
            DataFrame with reference sample as first entry
        self.experiment.all_samples : 
            adding this sample 
        '''
        print("Adding sample to MassGrid,", sample.name)
        mzlist = [x[0] for x in sample.track_mzs]
        new_reference_mzlist, new_reference_map2, updated_REF_landmarks, _r = \
            landmark_guided_mapping(
                list(self.MassGrid['mz']), 
                self._mz_landmarks_, 
                mzlist, 
                sample._mz_landmarks_,
                std_ppm = self.experiment.parameters['mz_tolerance_ppm'],
                correction_tolerance_ppm = self.experiment.parameters['correction_tolerance_ppm']
                )

        NewGrid = pd.DataFrame(
            np.full((len(new_reference_mzlist), 1+self._number_of_samples_), None),
            columns=['mz'] + self.list_sample_names,
        )
        NewGrid[ :self.MassGrid.shape[0]] = self.MassGrid
        NewGrid['mz'] = new_reference_mzlist
        NewGrid[ sample.name ] = new_reference_map2
        self.MassGrid = NewGrid
        self._mz_landmarks_ = updated_REF_landmarks
        sample.mz_calibration_ratio = _r            # not used now
        
        self.experiment.all_samples.append(sample)


    def bin_track_mzs(self, tl, reference_id=None):
        '''
        Bin all track m/z values into centroids via clustering, to be used to build massGrid.

        Parameters
        ----------
        tl : list[tuple]
            sorted list of all track m/z values in experiment, [(m/z, track_id, sample_id), ...]
        reference_id: str?
            the sample_id of reference sample. Not used now.
        
        Returns
        -------
        list of bins: 
            [ (mean_mz, [(), (), ...]), (mean_mz, [(), (), ...]), ... ]

            
        Note
        ----
            Because the range of each bin cannot be larger than mz_tolerance, 
            and mass tracks in each sample cannot overlap within mz_tolerance,
            multiple entries from the same sample in same bin will not happen.
            Similar to nearest neighbor (NN) clustering used in initial mass track construction.
        '''
        def __get_bin__(bin_data_tuples):
            return (np.median([x[0] for x in bin_data_tuples]), bin_data_tuples)

        tol_ = 0.000001 * self.experiment.parameters['mz_tolerance_ppm']
        bins_of_bins = []
        tmp = [tl[0]]
        for ii in range(1, len(tl)):
            _delta = tl[ii][0] - tl[ii-1][0]
            # bin adjacent tuples if they are within ppm tolerance
            if _delta < tol_ * tl[ii-1][0]:
                tmp.append(tl[ii])
            else:
                bins_of_bins.append(tmp)
                tmp = [tl[ii]]
        bins_of_bins.append(tmp)
        good_bins = []
        for bin_data_tuples in bins_of_bins:
            mz_range = bin_data_tuples[-1][0] - bin_data_tuples[0][0]
            mz_tolerance = bin_data_tuples[0][0] * tol_
            # important: double tol_ range here as mean_mz falls in single tol_
            if mz_range < mz_tolerance * 2:
                good_bins.append( __get_bin__(bin_data_tuples) )

            else:
                good_bins += [__get_bin__(C) for C in nn_cluster_by_mz_seeds(
                    bin_data_tuples, mz_tolerance)]

        return good_bins


    def join(self, M2):
        '''
        Placeholder. Future option to join with another MassGrid via a common reference.

        Parameters
        ----------
        M2: MassGrid instance
            the mass grid to be merged with this MassGrid
        '''
        pass


class CompositeMap:
    '''
    Each experiment is summarized into a CompositeMap (CMAP), as a master feature map.
    The use of CompositeMap also facilitates data visualization and exploration.
    Related concepts:

    i) MassGrid: a matrix for recording correspondence of mass tracks to each sample 

    ii) FeatureList: list of feature definitions, i.e. elution peaks defined on composite mass tracks.
    
    iii) FeatureTable: a matrix for feature intensities per sample.
    '''
    def __init__(self, experiment):
        '''
        Composite map of mass tracks and features, with pointers to individual samples.

        Parameters
        ----------
        experiment: ext_experiment instance
            the object representing the experiment from which to build the composite map
        '''
        self.experiment = experiment
        self._number_of_samples_ = experiment.number_of_samples
        self.list_sample_names = [experiment.sample_registry[ii]['name'] for ii in experiment.valid_sample_ids]

        # designated reference sample; all RT is aligned to this sample
        self.reference_sample_instance = self.reference_sample = \
            self.get_reference_sample_instance(experiment.reference_sample_id)
        self.rt_length = self.experiment.number_scans
        self.dict_scan_rtime = self.get_reference_rtimes(self.rt_length)
        self.max_ref_rtime = self.dict_scan_rtime[self.rt_length-1]

        self.MassGrid = None                        # will be pandas DataFrame, = MassGrid.MassGrid
        self.FeatureTable = None
        self.FeatureList = []

        self._mz_landmarks_ = []                    # m/z landmarks as index numbers
        self.good_reference_landmark_peaks = []     # used for RT alignment and m/z calibration to DB
        # self.reference_mzdict = {}
        self.composite_mass_tracks = {}             # following MassGrid indices

    def get_reference_sample_instance(self, reference_sample_id):
        '''
        Wraps the reference_sample into a SimpleSample instance, so that
        it have same behaivors as other samples.

        Parameters
        ----------
        reference_sample_id: any valid sample_id
            this is used to retrieve the sample from the experiment's sample_registry

        Returns
        -------
        instance of SimpleSample class for the reference_sample.
        '''
        SM = SimpleSample(self.experiment.sample_registry[reference_sample_id],
                experiment=self.experiment, database_mode=self.experiment.database_mode, 
                mode=self.experiment.mode,
                is_reference=True)
        SM.list_mass_tracks = SM.get_masstracks_and_anchors()
        return SM

    def get_reference_rtimes(self, rt_length):
        '''
        Extrapolate retention time on self.reference_sample_instance to max scan number in the experiment.
        This will be used to calculate retention time in the end, as intermediary steps use scan numbers.

        Parameters
        ----------
        rt_length: int
            this represents the total number of scans

        Returns
        -------
        dictionary of scan number to retetion time in the reference_sample.
        '''
        X, Y = self.reference_sample.rt_numbers, self.reference_sample.list_retention_time
        interf = interpolate.interp1d(X, Y, fill_value="extrapolate")
        newX = range(rt_length)
        newY = interf(newX)
        return dict(zip(newX, newY))

    def construct_mass_grid(self):
        '''
        Constructing MassGrid for the whole experiment. 
        If the sample number is no more than a predefined parameter ('project_sample_number_small', default 10), 
        this is considered a small study and a pairwise alignment is performed.
        See `MassGrid.build_grid_sample_wise`, `MassGrid.add_sample`.
        Else, for a larger study, the mass alignment is performed by the same NN clustering method 
        that is used in initial mass track construction. 
        See `MassGrid.build_grid_by_centroiding`, `MassGrid.bin_track_mzs`.
        
        Updates
        -------
        self._mz_landmarks_ : 
            landmark m/z values that match to 13C/12C and Na/H patterns
        self.MassGrid : 
            DataFrame with reference sample as first entry. Use sample name as column identifiers.
        
            
        Note
        ----
            Number of samples dictate workflow: 
            build_grid_by_centroiding is fast, but build_grid_sample_wise is used for small studies 
            to compensate limited size for statistical distribution.
            All mass tracks are included at this stage, regardless if peaks are detected, because
            peak detection will be an improved process on the composite tracks.
        '''
        print("Constructing MassGrid, ...")
        MG = MassGrid( self, self.experiment )
        if self._number_of_samples_ <= self.experiment.parameters['project_sample_number_small']:
            print("Building Grid Sample Wise...")
            MG.build_grid_sample_wise()
        else:
            print("Building Grid by Centroiding...")
            MG.build_grid_by_centroiding()
           
        self.MassGrid = MG.MassGrid
        self._mz_landmarks_ = MG._mz_landmarks_

    def mock_rentention_alignment(self):
        '''
        Create empty mapping dictionaries if the RT alignment fails, e.g. for blank or exogenous samples.
        '''
        for sample in self.experiment.all_samples[1:]:      # first sample is reference
            sample.rt_cal_dict, sample.reverse_rt_cal_dict = {}, {}

    def perform_index_alignment(self):
        mzDict = dict(self.MassGrid['mz'])
        mzlist = list(self.MassGrid.index)                          # this gets indices as keys, per mass track
        basetrack = np.zeros(self.rt_length, dtype=np.int64)        # self.rt_length defines max rt number
        _comp_dict = {}
        for k in mzlist: 
            _comp_dict[k] = basetrack.copy()

        index_samples = []
        for ii, SM in enumerate(self.experiment.all_samples):
            if ii not in self.experiment.mapping:
                index_samples.append(SM)

        self.good_reference_landmark_peaks = self.set_RT_reference(self.experiment.parameters['cal_min_peak_height'])

        # align index standards
        master_index_sample = index_samples[0]
        for index_sample in index_samples:
            index_mass_tracks = SimpleSample.get_mass_tracks_for_sample(index_sample)
            candidate_landmarks = [self.MassGrid[index_sample.name].values[p['ref_id_num']] for p in self.good_reference_landmark_peaks]
            good_landmark_peaks, selected_reference_landmark_peaks = [], []
            for jj in range(len(self.good_reference_landmark_peaks)):
                ii = candidate_landmarks[jj]
                if not pd.isna(ii):
                    ii = int(ii)
                    this_mass_track = index_mass_tracks[ii]
                    Upeak = quick_detect_unique_elution_peak(this_mass_track['intensity'], 
                                                            min_peak_height=self.experiment.parameters['cal_min_peak_height'], 
                                                            min_fwhm=3, 
                                                            min_prominence_threshold_ratio=0.2)
                    if Upeak:
                        scan_no_delta = Upeak['apex'] - self.good_reference_landmark_peaks[jj]['apex']
                        if abs(scan_no_delta) < np.inf:
                            Upeak.update({'ref_id_num': ii})
                            good_landmark_peaks.append(Upeak)
                            selected_reference_landmark_peaks.append(self.good_reference_landmark_peaks[jj])
            _NN = len(good_landmark_peaks)
            print("\tgood_landmarks: ", index_sample.name, _NN)
            
            from .chromatograms import clean_rt_calibration_points

            sample_rt_bound = max(index_sample.list_scan_numbers)
            rt_rightend_ = 1.1 * sample_rt_bound
            xx, yy = [-0.1 * sample_rt_bound,]*3, [-0.1 * sample_rt_bound,]*3
            rt_cal = clean_rt_calibration_points(
                [(x[0]['apex'], x[1]['apex']) for x in zip(good_landmark_peaks, selected_reference_landmark_peaks)]
            )
            xx += [L[0] for L in rt_cal] + [rt_rightend_]*3
            yy += [L[1] for L in rt_cal] + [rt_rightend_]*3
            # scale frac parameter like a sigmoid of number of data points when len(rt_cal) is in (50,150).
            FRAC = 0.6 - 0.004*(len(rt_cal)-50)
            FRAC = max(0.2, min(FRAC, 0.6))    # bound frac in (0.2, 0.6)

            lowess_predicted = __hacked_lowess__(yy, xx, frac=FRAC, it=3, xvals=index_sample.list_scan_numbers)
            interf = interpolate.interp1d(lowess_predicted, index_sample.list_scan_numbers, fill_value="extrapolate", bounds_error=False)
            ref_interpolated = interf( master_index_sample.list_scan_numbers )
            lowess_predicted = [int(round(ii)) for ii in lowess_predicted]

            rt_cal_dict = dict( 
                [(x,y) for x,y in zip(index_sample.list_scan_numbers, lowess_predicted) if x!=y and 0<=y<=max(master_index_sample.list_scan_numbers)] )

            ref_interpolated = [int(round(ii)) for ii in ref_interpolated]
            reverse_rt_cal_dict = dict(
                [(x,y) for x,y in zip(master_index_sample.list_scan_numbers, ref_interpolated) if x!=y and 0<=y<=sample_rt_bound] )
            
            index_sample.rt_cal_dict = rt_cal_dict
            index_sample.reverse_rt_cal_dict = reverse_rt_cal_dict

        for ii, SM in enumerate(self.experiment.all_samples):
            list_mass_tracks = SimpleSample.get_mass_tracks_for_sample(SM)
            print("Aligning: ", SM.name)
            if ii in self.experiment.mapping:
                index_sample = self.experiment.all_samples[self.experiment.mapping[ii]]

                # convert study sample to retention index
                list_retention_index = self.experiment.RI_models[self.experiment.mapping[ii]](SM.list_retention_time)

                # convert retention index to scans in RI sample
                list_reference_scans = self.experiment.reverse_RI_models[self.experiment.mapping[ii]](list_retention_index)

                rt_cal_dict = {}
                reverse_rt_cal_dict = {}
                for jj, ref_scan in enumerate(list_reference_scans):
                    rt_cal_dict[jj] = index_sample.rt_cal_dict.get(int(ref_scan), max(index_sample.list_scan_numbers))
                    reverse_rt_cal_dict[jj] = index_sample.reverse_rt_cal_dict.get(int(ref_scan), max(index_sample.list_scan_numbers))

                SM.rt_cal_dict = rt_cal_dict
                SM.reverse_rt_cal_dict = reverse_rt_cal_dict

                if not self.experiment.parameters['drop_unaligned_samples'] or SM.is_rt_aligned:
                    for k in mzlist:
                        ref_index = self.MassGrid[SM.name][k]
                        if not pd.isna(ref_index): # ref_index can be NA 
                            _comp_dict[k] += remap_intensity_track(list_mass_tracks[int(ref_index)]['intensity'],  basetrack.copy(), SM.rt_cal_dict)
        result = {k: {'id_number': k, 'mz': mzDict[k], 'intensity': v} for k,v in _comp_dict.items()}
        self.composite_mass_tracks = result

    def build_composite_tracks_GC(self):
        self.perform_index_alignment()

    def START(self):
        # Spanning Tree Alignment of Retention Time (START)

        # An alternative to traditional alignment based on a reference sample, 
        # rather, each sample may have a chain of reference samples back to the 
        # master reference sample. This requires more alignments per sample possibly,
        # but allows similar samples to align to one another before attempting to align
        # across sample types. For instance, each blank can align with other blanks, 
        # and the blank most like a biological sample, will be used to align the blanks
        # to study samples. 

        # First we need to estimate the 'goodness' of each possible alignment. This needs 
        # to be fast / simple enough to evaluate for all samples. Here we will use the number of 
        # shared anchor peaks to start.
        import seaborn as sns
        import matplotlib.pyplot as plt
        import networkx as nx
        import numpy as np
        
        CAL_MIN_PEAK_HEIGHT = self.experiment.parameters['cal_min_peak_height']
        MIN_PEAK_NUM = self.experiment.parameters['peak_number_rt_calibration']
        NUM_ITERATIONS = self.experiment.parameters['num_lowess_iterations']
        MIN_C_SELECTIVITY = 0.99

        mg = self.MassGrid.copy()
        selectivities = calculate_selectivity(mg['mz'], self.experiment.parameters['mz_tolerance_ppm'])
        mgd_inv = {i: x for i, x in enumerate(self.MassGrid.to_dict(orient='records'))}
        mgd = {}
        for index, mz_row in mgd_inv.items():
            for k, v in mz_row.items():
                if not pd.isna(v):
                    mgd[(k, v)] = index


        # find all candidate peaks for alignment
        reference_peaks_per_sample = {}
        for sample in self.experiment.all_samples:
            reference_peaks_per_sample[sample.name] = []
            mass_tracks = SimpleSample.get_mass_tracks_for_sample(sample)
            for index, mass_track in enumerate(mass_tracks):
                mapped_index = mgd.get((sample.name, index), None)
                if mapped_index is not None:
                    if selectivities[mapped_index] > MIN_C_SELECTIVITY:
                        Upeak = quick_detect_unique_elution_peak(mass_track['intensity'], 
                                    min_peak_height=CAL_MIN_PEAK_HEIGHT, 
                                    min_fwhm=3, 
                                    min_prominence_threshold_ratio=0.2)
                        if Upeak:
                            Upeak['mz'] = round(mass_track['mz'], 3)
                            Upeak['index'] = mapped_index
                            reference_peaks_per_sample[sample.name].append(Upeak)                    

        def __similarity(s1, s2):
            # similarity must give back a metric, the larger the metric, the better the alignment
            # similarity must be order invariant, i.e., F(s1, s2) = F(s2, s1)
            # spanning tree wants costs not similarity, so -F(s1, s2) is the cost of 
            # aligning s1, s2.

            s1_indices = set([x['index'] for x in reference_peaks_per_sample[s1.name]])
            s2_indices = set([x['index'] for x in reference_peaks_per_sample[s2.name]])
            return len(s1_indices.intersection(s2_indices))/len(s1_indices.union(s2_indices))
        
        def __cost(s1, s2):
            return 1-__similarity(s1, s2)
        
        def __pairwise_cost(samples):
            # cost is the negative of similarity
            return np.array([[__cost(s1, s2) for s2 in samples] for s1 in samples], dtype=np.float16)

        def __pairwise_similarity(samples):
            return np.array([[__similarity(s1, s2) for s2 in samples] for s1 in samples], dtype=np.float16)

        def __distance_to_graph(dmatrix):
            # convert to networkx graph
            # use networkx to find the minimum spanning tree
            # return graph
            # note that __similarity is positive, but spanning
            # tree assumes that the edge weights are cost.
            return nx.from_numpy_array(dmatrix)

        def __find_graph_root(dgraph):
            # the root is the node that is the closest to all
            # other nodes. It is the most central node. We can
            # remove a "layer" of leaf nodes until we are left 
            # with the root of the graph. 
            root = nx.center(dgraph)[0]
            a, b = __align_pair(self.experiment.all_samples[root], self.experiment.all_samples[root])
            self.experiment.all_samples[root].rt_cal_dict = a
            self.experiment.all_samples[root].reverse_rt_cal_dict = b
            return nx.center(dgraph)[0]

        def __pairwise_traverse(distance_graph, root, target):
            for path in nx.shortest_simple_paths(distance_graph, source=root, target=target):
                return path # its a spanning tree, there is only one path

        @lru_cache(maxsize=128)
        def __align_pair(s1, s2):
            print("Aligning: ", s1.name, " to ", s2.name)
            peaks_mzs_s1 = set([x['mz'] for x in reference_peaks_per_sample[s1.name]])
            peaks_mzs_s2 = set([x['mz'] for x in reference_peaks_per_sample[s2.name]])
            shared_peak_mzs = peaks_mzs_s1.intersection(peaks_mzs_s2)
            print("\tPeaks - Shared / Sample 1 / Sample 2: ", len(shared_peak_mzs), " / ", len(peaks_mzs_s1), " / ", len(peaks_mzs_s2))
            reference_pairs = {}
            for peak in reference_peaks_per_sample[s1.name]:
                if peak['mz'] in shared_peak_mzs:
                    if peak['mz'] not in reference_pairs:
                        reference_pairs[peak['mz']] = [None, None]
                    reference_pairs[peak['mz']][0] = peak['apex']
            for peak in reference_peaks_per_sample[s2.name]:
                if peak['mz'] in shared_peak_mzs:
                    reference_pairs[peak['mz']][1] = peak['apex']

            X, Y = [], []
            for mz, (apex1, apex2) in reference_pairs.items():
                X.append(apex1)
                Y.append(apex2)

            from .chromatograms import clean_rt_calibration_points
            reference_rt_numbers = s1.rt_numbers
            sample_rt_numbers = s2.rt_numbers
            reference_rt_bound = max(s1.rt_numbers)
            sample_rt_bound = max(s2.rt_numbers)
            rt_rightend_ = 1.1 * sample_rt_bound
            xx, yy = [-0.1 * sample_rt_bound,]*3, [-0.1 * sample_rt_bound,]*3
            rt_cal = clean_rt_calibration_points([(x[0], x[1]) for x in reference_pairs.values()])
            xx += [L[0] for L in rt_cal] + [rt_rightend_]*3
            yy += [L[1] for L in rt_cal] + [rt_rightend_]*3

            from statsmodels.nonparametric.smoothers_lowess import lowess
            lowess_predicted = __hacked_lowess__(yy, xx, frac=0.5, it=3, xvals=sample_rt_numbers)            
            # scale frac parameter like a sigmoid of number of data points when len(rt_cal) is in (50,150).

            interf = interpolate.interp1d(lowess_predicted, sample_rt_numbers, fill_value="extrapolate")
            ref_interpolated = interf( reference_rt_numbers )
            lowess_predicted = [ii for ii in lowess_predicted]

            rt_cal_dict = dict([(x,y) for x,y in zip(sample_rt_numbers, lowess_predicted) if x!=y and 0<=y<=reference_rt_bound] )

            ref_interpolated = [ii for ii in ref_interpolated]
            reverse_rt_cal_dict = dict([(x,y) for x,y in zip(reference_rt_numbers, ref_interpolated) if x!=y and 0<=y<=sample_rt_bound])
                
            return rt_cal_dict, reverse_rt_cal_dict

        alignment_cache = {}
        def __align(path):
            # walk each path, start with root as s1 and neighbors as
            # various s2s. Then continue by letting those s2 be s1s, 
            # and their neighbor nodes various s2s. 
            #
            # dictionary mapping scan nos in s1 to s2
            calibrated_domain = self.experiment.all_samples[path[0]].rt_numbers
            for i in range(len(path)-1):
                s1, s2 = self.experiment.all_samples[path[i]], self.experiment.all_samples[path[i+1]]
                calibrated_domain = [__align_pair(s1, s2)[0].get(x, x) for x in calibrated_domain]
            calibrated_domain = [int(round(x)) for x in calibrated_domain]
            rt_cal_dict = dict(zip(self.experiment.all_samples[path[0]].rt_numbers, calibrated_domain))
            reverse_rt_cal_dict = dict(zip(calibrated_domain, self.experiment.all_samples[path[0]].rt_numbers))
            s2.rt_cal_dict, s2.reverse_rt_cal_dict, s2.is_rt_aligned = rt_cal_dict, reverse_rt_cal_dict, True
            return rt_cal_dict, reverse_rt_cal_dict

        # With this similarity metric, build the similarity matrix between all samples.
        # Note that similarity is simple the inverse of distance, so simply flip sign to get
        # a distance matrix.

        # build similarity matrix
        # sample_similarity = np.zeros((len(self.experiment.all_samples), len(self.experiment.all_samples)), dtype=np.float64)

        D = __pairwise_cost(self.experiment.all_samples)
        G = __distance_to_graph(D)
        T = nx.minimum_spanning_tree(G)
        root = __find_graph_root(T)
        for node in G.nodes:
            if node != root:
                path = __pairwise_traverse(T, root, node)
                __align(path)


        mzDict = dict(self.MassGrid['mz'])
        mzlist = list(self.MassGrid.index)                          # this gets indices as keys, per mass track
        basetrack = np.zeros(self.rt_length, dtype=np.int64)        # self.rt_length defines max rt number
        
        _comp_dict = {}
        for k in mzlist: 
            _comp_dict[k] = basetrack.copy()

        for SM in self.experiment.all_samples:
            SM.is_rt_aligned = True
            list_mass_tracks = SimpleSample.get_mass_tracks_for_sample(SM)
            for k in mzlist:
                ref_index = self.MassGrid[SM.name][k]
                if not pd.isna(ref_index): # ref_index can be NA 
                    _comp_dict[k] += remap_intensity_track( 
                        list_mass_tracks[int(ref_index)]['intensity'],  
                        basetrack.copy(), SM.rt_cal_dict 
                        )
        result = {}
        for k,v in _comp_dict.items():
            result[k] = {
                'id_number': k, 
                'mz': mzDict[k], 
                'intensity': v
                }
        self.composite_mass_tracks = result

    def build_composite_tracks(self):
        '''
        Perform RT calibration then make composite tracks.

        Updates
        -------
        self.good_reference_landmark_peaks : 
            [{'ref_id_num': 99, 'apex': 211, 'height': 999999}, ...]
        self.composite_mass_tracks : 
            list of composite mass tracks in this experiment.
        sample.rt_cal_dict and sample.reverse_rt_cal_dict for all samples.

        
        Note
        ----
        See calibrate_sample_RT for details in RT alignment. 
        '''
        print("\nBuilding composite mass tracks and calibrating retention time ...\n")

        cal_min_peak_height = self.experiment.parameters['cal_min_peak_height']
        MIN_PEAK_NUM = self.experiment.parameters['peak_number_rt_calibration']
        NUM_ITERATIONS = self.experiment.parameters['num_lowess_iterations']
        if self.experiment.parameters['max_retention_shift'] is None:
            MAX_RETENTION_SHIFT = np.inf
        else:
            MAX_RETENTION_SHIFT = self.experiment.parameters['max_retention_shift']

        self.good_reference_landmark_peaks = self.set_RT_reference(cal_min_peak_height)
        
        mzDict = dict(self.MassGrid['mz'])
        mzlist = list(self.MassGrid.index)                          # this gets indices as keys, per mass track
        basetrack = np.zeros(self.rt_length, dtype=np.int64)        # self.rt_length defines max rt number
        _comp_dict = {}
        for k in mzlist: 
            _comp_dict[k] = basetrack.copy()

        # add to export mz and rtime of good reference landmarks
        if self.experiment.parameters['debug_rtime_align']:
            self.export_reference_sample()

        batches = [[]]
        for SM in self.experiment.all_samples:
            if len(batches[-1]) == self.experiment.parameters['multicores']:
                batches.append([])
            batches[-1].append(SM)
        assert batches[-1], "Empty batch"

        for batch in batches:
            if self.experiment.parameters['database_mode'] == "memory":
                # this is a bug - to fix
                list_of_list_mass_tracks = bulk_process(SimpleSample.get_mass_tracks_for_sample, 
                                                        batch, 
                                                        dask_ip=False,
                                                        jobs_per_worker=self.experiment.paramters['multicores'])
            else:
                list_of_list_mass_tracks = bulk_process(SimpleSample.get_mass_tracks_for_sample, 
                                                        batch, 
                                                        dask_ip=self.experiment.parameters['dask_ip'],
                                                        jobs_per_worker=self.experiment.parameters['multicores'])
            for (SM, list_mass_tracks) in zip(batch, list_of_list_mass_tracks):
                print("   ", SM.name)
                if SM.is_reference:
                    print("\t\tgood_reference_landmark_peaks: ", len(self.good_reference_landmark_peaks))
                else:
                    if self.experiment.parameters['rt_align_on']:
                        if self.experiment.parameters['debug_rtime_align']:
                            cal_func = rt_lowess_calibration_debug
                        else:
                            cal_func = rt_lowess_calibration

                        self.calibrate_sample_RT(SM, list_mass_tracks, 
                                            calibration_fuction=cal_func,
                                            cal_min_peak_height=cal_min_peak_height, 
                                            MIN_PEAK_NUM=MIN_PEAK_NUM,
                                            MAX_RETENTION_SHIFT=MAX_RETENTION_SHIFT,
                                            NUM_ITERATIONS=NUM_ITERATIONS)

                # option to skip sample if not aligned
                if not self.experiment.parameters['drop_unaligned_samples'] or SM.is_rt_aligned:
                    for k in mzlist:
                        ref_index = self.MassGrid[SM.name][k]
                        if not pd.isna(ref_index): # ref_index can be NA 
                            _comp_dict[k] += remap_intensity_track( 
                                list_mass_tracks[int(ref_index)]['intensity'],  
                                basetrack.copy(), SM.rt_cal_dict 
                                )

        result = {}
        for k,v in _comp_dict.items():
            result[k] = { 'id_number': k, 'mz': mzDict[k], 'intensity': v }

        self.composite_mass_tracks = result

    def calibrate_sample_RT_by_standards(self, sample):
        '''
        Placeholder, to add RT calibration based on spike-in compound standards.

        Parameters
        ----------
        sample: 
            this will either be a SimpleSample object for the sample containing 
            the spike-in standards.
        '''
        pass


    def calibrate_sample_RT(self, 
                            sample, 
                            list_mass_tracks,
                            calibration_fuction=rt_lowess_calibration, 
                            cal_min_peak_height=100000,
                            MIN_PEAK_NUM=15,
                            MAX_RETENTION_SHIFT=np.inf,
                            NUM_ITERATIONS=3):
        '''
        Calibrate/align retention time per sample.

        Parameters
        ----------
        sample : SimpleSample instance
            instance of SimpleSample class
        list_mass_tracks : list
            list of mass tracks in sample. 
            This may not be kept in memeory with the sample instance, thus require retrieval.
        calibration_fuction : function, optional, default: rt_lowess_calibration
            RT calibration fuction to use, default to rt_lowess_calibration.
        cal_min_peak_height : float, optional, default: 100000
            minimal height required for a peak to be used for calibration.
            Only high-quality peaks unique in each mass track are used for calibration.
        MIN_PEAK_NUM : int, optional, default: 15
            minimal number of peaks required for calibration. Abort if not met.

        Updates
        -------
        sample.rt_cal_dict :   
            dictionary converting scan number in sample_rt_numbers to 
            calibrated integer values in self.reference_sample.
            Range matched. Only changed numbers are kept for efficiency.
        sample.reverse_rt_cal_dict : 
            dictionary from ref RT scan numbers to sample RT scan numbers. 
            Range matched. Only changed numbers are kept for efficiency.
        sample.rt_landmarks : 
            list of apex scan numbers for the peaks used in RT calibration.
            
        Note
        ----
            This is based on a set of unambiguous peaks: quich peak detection on anchor mass trakcs, 
            and peaks that are unique to each track are used for RT alignment.
            Only numbers different btw two samples are kept in the dictionaries for computing efficiency.
            When calibration_fuction fails, e.g. inf on lowess_predicted,
            it is assumed that this sample is not amendable to computational alignment,
            and the sample will be attached later without adjusting retention time.
            It will be good to have good_landmark_peaks to cover RT range evenly in the future.
            Using user-supplied internal standards will be an important option.
        '''

        candidate_landmarks = [self.MassGrid[sample.name].values[
                                p['ref_id_num']] for p in 
                                self.good_reference_landmark_peaks] # contains NaN
        good_landmark_peaks, selected_reference_landmark_peaks = [], []
        for jj in range(len(self.good_reference_landmark_peaks)):
            ii = candidate_landmarks[jj]
            if not pd.isna(ii):
                ii = int(ii)
                this_mass_track = list_mass_tracks[ii]
                Upeak = quick_detect_unique_elution_peak(this_mass_track['intensity'], 
                            min_peak_height=cal_min_peak_height, 
                            min_fwhm=3, min_prominence_threshold_ratio=0.2)
                
                if Upeak:
                    scan_no_delta = Upeak['apex'] - self.good_reference_landmark_peaks[jj]['apex']
                    if abs(scan_no_delta) < MAX_RETENTION_SHIFT:
                        Upeak.update({'ref_id_num': ii})
                        good_landmark_peaks.append(Upeak)
                        selected_reference_landmark_peaks.append(self.good_reference_landmark_peaks[jj])

        _NN = len(good_landmark_peaks)
        print("\tgood_landmark_peaks: ", _NN)

        sample.rt_landmarks = [p['apex'] for p in good_landmark_peaks]
        # only do RT calibration if MIN_PEAK_NUM is met.
        if _NN >  MIN_PEAK_NUM:
            try:
                sample.rt_cal_dict, sample.reverse_rt_cal_dict = calibration_fuction( 
                                        good_landmark_peaks, selected_reference_landmark_peaks, 
                                        sample.rt_numbers, self.reference_sample.rt_numbers, NUM_ITERATIONS, sample.name,
                                        self.experiment.parameters['outdir'])
                sample.is_rt_aligned = True
            except OverflowError:
                pass
            
        if not sample.is_rt_aligned:
                sample.rt_cal_dict, sample.reverse_rt_cal_dict =  {}, {}
                print("    ~warning~ Failure in retention time alignment (%d); %s." 
                                            %( _NN, sample.name))

    def set_RT_reference(self, cal_peak_intensity_threshold=100000):
        '''
        Start with the referecne samples, usually set for a sample of most landmark mass tracks.
        Do a quick peak detection for good peaks; use high selectivity m/z to avoid ambiguity 
        in peak definitions.

        Parameters
        ----------
        cal_peak_intensity_threshold: float, optional, default: 100000
            a peak must have an intensity above this value to be used as an RT_reference

        Returns
        ------- 
        good_reference_landmark_peaks: [{'ref_id_num': 99, 'apex': 211, 'height': 999999}, ...]

        Note
        ----
        Some members in good_reference_landmark_peaks may have the same RT apex.
        But the redundant numbers should be handled by rt_lowess_calibration, in which .frac is
        more important for stability.
        '''
        selectivities = calculate_selectivity(self.MassGrid['mz'][self._mz_landmarks_], self.experiment.parameters['mz_tolerance_ppm'])
        good_reference_landmark_peaks = []
        ref_list_mass_tracks = self.reference_sample.list_mass_tracks

        for ii, mz_landmark in enumerate(self._mz_landmarks_):
            if selectivities[ii] > 0.99:
                ref_ii = self.MassGrid[self.reference_sample.name][mz_landmark]
                if ref_ii and not pd.isna(ref_ii):
                    this_mass_track = ref_list_mass_tracks[int(ref_ii)]
                    Upeak = quick_detect_unique_elution_peak(this_mass_track['intensity'], 
                                min_peak_height=cal_peak_intensity_threshold, 
                                min_fwhm=3, min_prominence_threshold_ratio=0.2)
                    if Upeak:
                        Upeak.update({'ref_id_num': self._mz_landmarks_[ii]}) # as in MassGrid index
                        good_reference_landmark_peaks.append(Upeak)

        self.reference_sample.rt_landmarks = [p['apex'] for p in good_reference_landmark_peaks]
        
        return good_reference_landmark_peaks


    def global_peak_detection(self):
        '''
        Detects elution peaks on composite mass tracks, resulting to a list of features.
        Using peaks.batch_deep_detect_elution_peaks for parallel processing.

        Updates
        -------
        self.FeatureList :
            a list of JSON peaks
        self.FeatureTable : 
            a pandas dataframe for features across all samples.

            
        Note
        ----
            Because the composite mass tracks ar summarized on all samples, 
            the resulting elution peaks are really features at the experiment level.
            Peak area and height are cumulated from all samples, 
            not average because some peaks are in only few samples.
        '''
        print("\nPeak detection on %d composite mass tracks, ...\n" %len(self.composite_mass_tracks))

        self.FeatureList = batch_deep_detect_elution_peaks(
            self.composite_mass_tracks.values(), 
            self.experiment.number_scans, 
            self.experiment.parameters
        )
        for ii, peak in enumerate(self.FeatureList):
            peak['id_number'] = 'F'+str(ii)
            # convert scan numbers to rtime
            try:
                peak['rtime'] = self.dict_scan_rtime[peak['apex']]
            except KeyError:
                peak['rtime'] = self.max_ref_rtime                # imputed value set at max rtime
                print("Feature rtime out of bound - ", peak['id_number'], peak['apex'])
            try:
                peak['rtime_left_base'], peak['rtime_right_base'] = self.dict_scan_rtime[peak['left_base']], self.dict_scan_rtime[peak['right_base']]
            except KeyError:
                print("Feature rtime out of bound on", peak['id_number'], (peak['apex'], peak['left_base'], peak['right_base']))

        self.generate_feature_table()


    def get_peak_area_sum(self, track_intensity, left_base, right_base):
        '''
        Option to calculate peak area by sum of the intensity values on the track 
        within the peak boundaries.

        Parameters
        ----------
        track_intensity : np.array[dtype=INTENSITY_DATA_TYPE]
            np.array, i.e. mass_track['intensity']
        left_base : int
            index for peak left base
        right_base : int 
            index for peak right base

        Returns
        ------- 
        Integer of peak area value
        '''
        if isinstance(left_base, float):
            left_base = int(left_base)
        if isinstance(right_base, float):
            right_base = int(right_base)
        return track_intensity[left_base: right_base+1].sum()
    
    
    def get_peak_area_auc(self, track_intensity, left_base, right_base):
        '''
        Option to calculate peak area as area under the curve.
        This is approximated by a maximum filter to cover potential gaps.

        Parameters
        ----------
        track_intensity : np.array[dtype=INTENSITY_DATA_TYPE]
            np.array, i.e. mass_track['intensity']
        left_base : int
            index for peak left base
        right_base : int
            index for peak right base

        Returns
        ------- 
        Integer of peak area value
        
        '''
        return int(maximum_filter1d(track_intensity[left_base: right_base+1], size=2, mode='constant').sum())


    def get_peak_area_gaussian(self, track_intensity, left_base, right_base):
        '''
        Option to calculate peak area by fitting the data to a gaussian model.
        This is 

        Parameters
        ----------
        track_intensity : np.array[dtype=INTENSITY_DATA_TYPE]
            np.array, i.e. mass_track['intensity']
        left_base : int
            index for peak left base
        right_base : int
            index for peak right base

        Returns
        ------- 
        peak area, Integer value as gaussian integral.
        '''
        return int(get_gaussian_peakarea_on_intensity_list(track_intensity, left_base, right_base))


    def generate_feature_table(self):
        '''
        Initiate and populate self.FeatureTable, each sample per column in dataframe.
        '''
        peak_area_methods = {
            'auc': self.get_peak_area_auc,
            'sum': self.get_peak_area_sum,
            'gauss': self.get_peak_area_gaussian,
        }
        peak_area_function = peak_area_methods[self.experiment.parameters['peak_area']]

        FeatureTable = pd.DataFrame(self.FeatureList)
        for SM in self.experiment.all_samples:
            if not self.experiment.parameters['drop_unaligned_samples'] or SM.is_rt_aligned:
                FeatureTable[SM.name] = self.extract_features_per_sample(SM, peak_area_function)
        print("\nFeature Table: ", FeatureTable.shape)
        self.FeatureTable = FeatureTable


    def extract_features_per_sample(self, sample, peak_area_function):
        '''
        Extract and return peak area values in a sample, 
        based on the start and end positions defined in self.FeatureList.
        A peak area could be 0 if no real peak is present for a feature in this sample.

        Parameters
        ----------
        sample : SimpleSample instance
            instance of SimpleSample class.
        peak_area_function : function
            function to be used for peak area calculation

        Returns
        ------- 
        A list of peak area values, for all features in a sample.
        '''
        fList = []
        list_mass_tracks = sample.get_masstracks_and_anchors()
        for peak in self.FeatureList:
            track_number = self.MassGrid[sample.name][peak['parent_masstrack_id']]
            if not pd.isna(track_number):
                mass_track = list_mass_tracks[int(track_number)]
                left_base = sample.reverse_rt_cal_dict.get(peak['left_base'], peak['left_base'])
                right_base = sample.reverse_rt_cal_dict.get(peak['right_base'], peak['right_base'])
                peak_area = peak_area_function(mass_track['intensity'], left_base, right_base)
            else:
                peak_area = 0
            fList.append( peak_area )
        return fList
    
    def export_reference_sample(self):
        """Write mz and retention time of "good" ions to csv in reference sample

        Results
        -------
        mz,rtime
        84.04437446594238,196.3507106869998
        85.04770363867283,197.100775215
        90.05493021011353,160.75314731200018
        100.11204060912132,18.757312656
        101.11540949344635,19.138889808
        104.9922667145729,147.4066373920002
        105.99559181928635,147.7856911519998
        112.09949165582657,255.0619356640002
        114.06613251566887,74.11716273600001
        ......

        The file name would be reference sample name + _mz_rtime_landmarks under export dir
         
        """
        # extendable. could add height and other params
        mz_landmarks = [self.MassGrid['mz'].values[p['ref_id_num']] for p in self.good_reference_landmark_peaks] 
        rtime_landmarks = [self.dict_scan_rtime[p['apex']] for p in self.good_reference_landmark_peaks]
        reference_sample_name = self.reference_sample.name

        # example: batch14_MT_20210808_087_mz_rtime_landmarks.csv
        reference_path = os.path.join(self.experiment.parameters['outdir'], 'export', reference_sample_name + '_mz_rtime_landmarks.csv')
        with open(reference_path, 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(["mz", "rtime"])  #  headers
            writer.writerows(zip(mz_landmarks, rtime_landmarks)) 
