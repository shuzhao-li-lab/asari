"""
Classes of MassGrid and CompositeMap.
"""

import os
import csv
from functools import partial

import pandas as pd
import numpy as np
from scipy import interpolate
from scipy.ndimage import maximum_filter1d
from mass2chem.search import find_mzdiff_pairs_from_masstracks

from .mass_functions import (flatten_tuplelist,
                            landmark_guided_mapping,
                            calculate_selectivity)
from .chromatograms import (INTENSITY_DATA_TYPE,
                            nn_cluster_by_mz_seeds,
                            rt_lowess_calibration_debug,
                            rt_lowess_calibration,
                            remap_intensity_track,
                            __hacked_lowess__)
from .peaks import (quick_detect_unique_elution_peak,
                    batch_deep_detect_elution_peaks,
                    get_gaussian_peakarea_on_intensity_list)
from .samples import SimpleSample

from .utils import bulk_process
from .tools.ms2 import rt_cluster_msms


def nearest_peak_intensity(intensity, rtimes, rt, window=5):
    """
    Height of the MS1 peak nearest to rt on a mass track (precursor intensity of an MS/MS spectrum).

    intensity : intensity array of the mass track, co-indexed with rtimes (sample scans).
    rtimes : retention times of MS1 scans in seconds, ascending.
    rt : retention time of the MS/MS spectrum, in seconds.
    window : seconds on each side of rt to search.

    Returns the intensity of the local maximum closest to rt within the window
    (the window maximum if no local maximum inside it); 0 if no signal.
    """
    lo = int(np.searchsorted(rtimes, rt - window, side='left'))
    hi = min(int(np.searchsorted(rtimes, rt + window, side='right')), len(intensity))
    if hi <= lo:
        return 0
    seg = np.asarray(intensity[lo:hi], dtype=float)
    if seg.max() <= 0:
        return 0
    if len(seg) < 3:
        return float(seg.max())
    inner = np.where((seg[1:-1] > 0) & (seg[1:-1] >= seg[:-2]) & (seg[1:-1] >= seg[2:]))[0] + 1
    if not len(inner):
        return float(seg.max())
    best = inner[np.argmin(np.abs(np.asarray(rtimes[lo:hi])[inner] - rt))]
    return float(seg[best])


def _cluster_track_ms2_worker(item, similarity_function, rt_gap, mz_tolerance, distance_threshold,
                              method='greedy', similarity_threshold=0.7, score='entropy'):
    key, track = item
    if track.get('ms2_spectra'):
        track['ms2_spectra'] = rt_cluster_msms(
            track['ms2_spectra'], similarity_function, rt_gap, mz_tolerance, distance_threshold,
            method=method, similarity_threshold=similarity_threshold, score=score
        )
    return key, track


class MassGrid:
    """
    MassGrid is the concept for m/z correspondence in asari.
    This shares similarity to FeatureMap in OpenMS, but the correspondence
    in asari takes adavantage of high m/z resolution first before feature detection.
    """
    def __init__(self, cmap, experiment):
        """
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
        """

        self.experiment = experiment
        self.CMAP = cmap
        self.reference_sample_instance = self.CMAP.reference_sample_instance
        self.max_ref_rtime = self.CMAP.max_ref_rtime
        self.list_sample_names = self.CMAP.list_sample_names
        self._number_of_samples_ = self.CMAP._number_of_samples_

    def build_grid_sample_wise(self):
        """
        Align one sample at a time to reference m/z grid, based on their anchor m/z tracks.
        One of the two methods to build the grid.
        This is better for reliable assembly of small number of samples.

        See also
        --------
        build_grid_by_centroiding
        """
        self._initiate_mass_grid()
        # exclude reference sample by value (sample ids are not list positions when some
        # samples failed extraction); a new list, not modifying experiment.valid_sample_ids
        sample_ids = [sid for sid in self.experiment.valid_sample_ids
                      if sid != self.experiment.reference_sample_id]
        for sid in sample_ids:
            SM = SimpleSample(self.experiment.sample_registry[sid],
                experiment=self.experiment,
                database_mode=self.experiment.database_mode,
                mode=self.experiment.mode)
            self.add_sample(SM)

    def build_grid_by_centroiding(self):
        """
        Assemble mass grid by grouping m/z values to centroids.
        Each centroid can have no more than one mass track per sample.
        One of the two methods to build the grid.
        This is more efficient for large number of samples.

        See also
        --------
        build_grid_sample_wise
        """
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
        """
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
        """
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
        """
        This adds a sample to MassGrid, including the m/z alignment of the sample against the
        existing reference m/z values in the MassGrid.

        Parameters
        ----------
        sample : SimpleSample instance
            instance of SimpleSample class.

        Updates
        -------
        self._mz_landmarks_ :
            landmark m/z values that match to 13C/12C and Na/H patterns
        self.MassGrid :
            DataFrame with reference sample as first entry
        self.experiment.all_samples :
            adding this sample
        """
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
        NewGrid[ :self.MassGrid.shape[0]] = self.MassGrid   # replicate last grid
        NewGrid['mz'] = new_reference_mzlist                # update ref m/z list
        NewGrid[ sample.name ] = new_reference_map2         # update col of this sample
        self.MassGrid = NewGrid
        self._mz_landmarks_ = updated_REF_landmarks
        sample.mz_calibration_ratio = _r                    # not used now

        self.experiment.all_samples.append(sample)


    def bin_track_mzs(self, tl, reference_id=None):
        """
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
        """
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
        """
        Placeholder. Future option to join with another MassGrid via a common reference.

        Parameters
        ----------
        M2: MassGrid instance
            the mass grid to be merged with this MassGrid
        """


        pass


class CompositeMap:
    """
    Each experiment is summarized into a CompositeMap (CMAP), as a master feature map.
    The use of CompositeMap also facilitates data visualization and exploration.
    Related concepts:

    i) MassGrid: a matrix for recording correspondence of mass tracks to each sample

    ii) FeatureList: list of feature definitions, i.e. elution peaks defined on composite mass tracks.

    iii) FeatureTable: a matrix for feature intensities per sample.
    """
    def __init__(self, experiment):
        """
        Composite map of mass tracks and features, with pointers to individual samples.

        Parameters
        ----------
        experiment: ext_experiment instance
            the object representing the experiment from which to build the composite map
        """
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
        """
        Wraps the reference_sample into a SimpleSample instance, so that
        it have same behaivors as other samples.

        Parameters
        ----------
        reference_sample_id: any valid sample_id
            this is used to retrieve the sample from the experiment's sample_registry

        Returns
        -------
        instance of SimpleSample class for the reference_sample.
        """
        SM = SimpleSample(self.experiment.sample_registry[reference_sample_id],
                experiment=self.experiment, database_mode=self.experiment.database_mode,
                mode=self.experiment.mode,
                is_reference=True)
        SM.list_mass_tracks = SM.get_masstracks_and_anchors()
        return SM

    def get_reference_rtimes(self, rt_length):
        """
        Extrapolate retention time on self.reference_sample_instance to max scan number in the experiment.
        This will be used to calculate retention time in the end, as intermediary steps use scan numbers.
        This step of standardization is needed because
        scan numbers may not be consistent between samples, or continuous within samples.

        Parameters
        ----------
        rt_length: int
            this represents the total number of scans

        Returns
        -------
        dictionary of scan number to retetion time in the reference_sample.
        """
        X, Y = self.reference_sample.rt_numbers, self.reference_sample.list_retention_time
        interf = interpolate.interp1d(X, Y, fill_value="extrapolate")
        newX = range(rt_length)
        newY = interf(newX)
        return dict(zip(newX, newY))

    def construct_mass_grid(self):
        """
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
        """
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
        """
        Create empty mapping dictionaries if the RT alignment fails, e.g. for blank or exogenous samples.
        """
        for sample in self.experiment.all_samples[1:]:      # first sample is reference
            sample.rt_cal_dict, sample.reverse_rt_cal_dict = {}, {}

    def build_composite_tracks(self):
        """
        Perform RT calibration then make composite tracks.
        MS2 data are retrieved and attched composite to mass tracks here,
        otherwise requiring loading data likely from disk.

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
        """
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
        # LCMSMS precursor_features: MS/MS precursors per composite track, (sample, ref scan, MS1 intensity)
        collect_precursors = self.experiment.parameters.get('precursor_features', False)
        precursor_rt_window = self.experiment.parameters.get('precursor_rt_window', 5)
        self.ms2_precursors = {}
        basetrack = np.zeros(self.rt_length, dtype=INTENSITY_DATA_TYPE)        # self.rt_length defines max rt number
        _comp_dict, _comp_ms2 = {}, {}
        for k in mzlist:
            _comp_dict[k] = basetrack.copy()
            _comp_ms2[k] = []

        # add to export mz and rtime of good reference landmarks
        if self.experiment.parameters['debug_rtime_align']:
            self.export_reference_sample()
        #
        #  do not use bulk_process, too much slower than loading each pickle, which is limited by disk io anyway
        #
        for SM in self.experiment.all_samples:
            print("   ", SM.name)
            list_mass_tracks = SM.get_masstracks_and_anchors()

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
                        _this_track = list_mass_tracks[int(ref_index)]
                        _comp_dict[k] += remap_intensity_track(
                            _this_track['intensity'],
                            basetrack.copy(), SM.rt_cal_dict
                            )
                        for spec in _this_track['ms2_spectra']:
                            spec['sample'] = SM.name
                        _comp_ms2[k] += _this_track['ms2_spectra']
                        if collect_precursors and _this_track['ms2_spectra']:
                            self.ms2_precursors.setdefault(k, []).extend(
                                self._get_precursors(SM, _this_track, precursor_rt_window))

        result = {}
        for k,v in _comp_dict.items():
            result[k] = { 'id_number': k, 'mz': mzDict[k],
                         'intensity': v,                # MS1 intensity vector
                         'ms2_spectra': _comp_ms2[k],   # MS2 spectra list
                          }

        self.composite_mass_tracks = result

    @staticmethod
    def _get_precursors(sample, track, window):
        """
        Precursors of MS/MS spectra on a sample mass track, as (sample name, reference scan, MS1 intensity).
        Reference scan is the RT-aligned scan number (sample.rt_cal_dict), MS1 intensity by nearest_peak_intensity.
        """
        rtimes = np.asarray(sample.list_retention_time)
        result = []
        for spec in track['ms2_spectra']:
            scan = int(np.argmin(np.abs(rtimes - spec['rtime'])))
            ref_scan = sample.rt_cal_dict.get(scan, scan) if sample.rt_cal_dict else scan
            result.append((sample.name, int(round(ref_scan)),
                           nearest_peak_intensity(track['intensity'], rtimes, spec['rtime'], window)))
        return result

    def get_precursor_feature_table(self):
        """
        LCMSMS feature table from MS/MS precursors, without elution peak detection (parameters['precursor_features']).

        On each composite mass track, MS/MS precursors (from all samples, RT-aligned) are grouped by RT,
        a new group when the gap exceeds parameters['ms2_rt_gap'] (default 5 seconds); each group is a feature.
        Feature intensity per sample is the height of the MS1 peak nearest to the MS/MS precursor
        (nearest_peak_intensity, within parameters['precursor_rt_window'] seconds):
        the max over the sample's own MS/MS spectra in the group; for samples without MS/MS in the group,
        at the feature RT mapped to the sample (reverse RT calibration).
        Feature rtime is the median aligned RT of precursors, rtime_left_base and rtime_right_base the range.
        peak_area is the sum of sample intensities; cSelectivity, goodness_fitting and snr are not computed (0).
        Features without MS1 signal in any sample are dropped. Feature IDs follow track m/z and RT.

        Updates self.FeatureList, self.FeatureTable.
        """
        rt_gap = self.experiment.parameters.get('ms2_rt_gap', 5)
        window = self.experiment.parameters.get('precursor_rt_window', 5)
        max_scan = self.rt_length - 1
        ref_rt = lambda scan: self.dict_scan_rtime.get(min(max(int(scan), 0), max_scan), self.max_ref_rtime)
        self.FeatureList, group_intensities = [], []
        for k in sorted(self.ms2_precursors, key=lambda x: self.composite_mass_tracks[x]['mz']):
            precursors = sorted(self.ms2_precursors[k], key=lambda x: x[1])
            groups, current = [], [precursors[0]]
            for p in precursors[1:]:
                if ref_rt(p[1]) - ref_rt(current[-1][1]) > rt_gap:
                    groups.append(current)
                    current = [p]
                else:
                    current.append(p)
            groups.append(current)
            for G in groups:
                scans = [p[1] for p in G]
                apex = int(np.median(scans))
                intensities = {}
                for sample_name, _, x in G:
                    intensities[sample_name] = max(intensities.get(sample_name, 0), x)
                self.FeatureList.append({
                    'id_number': 'F' + str(len(self.FeatureList)),
                    'parent_masstrack_id': k,
                    'mz': self.composite_mass_tracks[k]['mz'],
                    'apex': apex, 'left_base': min(scans), 'right_base': max(scans),
                    'rtime': ref_rt(apex), 'rtime_left_base': ref_rt(min(scans)), 'rtime_right_base': ref_rt(max(scans)),
                    'peak_area': 0, 'cSelectivity': 0, 'goodness_fitting': 0, 'snr': 0,
                    'number_ms2_spectra': len(G),
                })
                group_intensities.append(intensities)
        print("\n%d precursor features from %d MS/MS spectra on %d mass tracks." %(
            len(self.FeatureList), sum(len(v) for v in self.ms2_precursors.values()), len(self.ms2_precursors)))

        FeatureTable = pd.DataFrame(self.FeatureList)
        total = np.zeros(len(self.FeatureList))
        for SM in self.experiment.all_samples:
            if self.experiment.parameters['drop_unaligned_samples'] and not SM.is_rt_aligned:
                continue
            list_mass_tracks = SM.get_masstracks_and_anchors()
            rtimes = np.asarray(SM.list_retention_time)
            fList = []
            for F, intensities in zip(self.FeatureList, group_intensities):
                if SM.name in intensities:
                    fList.append(intensities[SM.name])
                    continue
                track_number = self.MassGrid[SM.name][F['parent_masstrack_id']]
                if pd.isna(track_number):
                    fList.append(0)
                    continue
                scan = SM.reverse_rt_cal_dict.get(F['apex'], F['apex']) if SM.reverse_rt_cal_dict else F['apex']
                scan = min(max(int(round(scan)), 0), len(rtimes) - 1)
                fList.append(nearest_peak_intensity(
                    list_mass_tracks[int(track_number)]['intensity'], rtimes, rtimes[scan], window))
            FeatureTable[SM.name] = fList
            total += np.array(fList, dtype=float)
        FeatureTable['peak_area'] = total.round()
        keep = total > 0            # drop features without MS1 signal in any sample
        if (~keep).sum():
            print("Dropped %d precursor features without MS1 signal." %(~keep).sum())
        self.FeatureTable = FeatureTable[keep].reset_index(drop=True)
        self.FeatureList = [F for F, k in zip(self.FeatureList, keep) if k]

    def cluster_ms2_spectra(self,
                            similarity_function=None,
                            rt_gap=5,
                            mz_tolerance=0.01,
                            distance_threshold=0.8,
                            method='greedy',
                            similarity_threshold=0.7,
                            score='entropy'):
        """
        Cluster all MS/MS spectra on each mass track,
        replace track['ms2_spectra'] by cluster results (extended representative spectrum).

        similarity_function : default None uses built-in score.
        rt_gap : used to divide clusters on retention time.
        mz_tolerance : m/z tolerance used in MS/MS similarity calculation.
        distance_threshold : for method 'hcl'.
        method : 'greedy' (consensus) or 'hcl'.
        similarity_threshold, score : for method 'greedy'.
        """
        print("\nClustering MS/MS spectra per mass track (%s) ...\n" %method)
        worker = partial(_cluster_track_ms2_worker,
                         similarity_function=similarity_function,
                         rt_gap=rt_gap,
                         mz_tolerance=mz_tolerance,
                         distance_threshold=distance_threshold,
                         method=method,
                         similarity_threshold=similarity_threshold,
                         score=score)
        results = bulk_process(worker, list(self.composite_mass_tracks.items()))
        self.composite_mass_tracks = dict(results)


    def calibrate_sample_RT_by_standards(self, sample):
        """
        Placeholder, to add RT calibration based on spike-in compound standards.

        Parameters
        ----------
        sample:
            this will either be a SimpleSample object for the sample containing
            the spike-in standards.
        """
        pass


    def calibrate_sample_RT(self,
                            sample,
                            list_mass_tracks,
                            calibration_fuction=rt_lowess_calibration,
                            cal_min_peak_height=100000,
                            MIN_PEAK_NUM=15,
                            MAX_RETENTION_SHIFT=np.inf,
                            NUM_ITERATIONS=3):
        """
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
        """

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
        """
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
        """
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
        """
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
        """
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
        """
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
        """
        if isinstance(left_base, float):
            left_base = int(left_base)
        if isinstance(right_base, float):
            right_base = int(right_base)
        return track_intensity[left_base: right_base+1].sum()


    def get_peak_area_auc(self, track_intensity, left_base, right_base):
        """
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

        """
        return int(maximum_filter1d(track_intensity[left_base: right_base+1], size=2, mode='constant').sum())


    def get_peak_area_gaussian(self, track_intensity, left_base, right_base):
        """
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
        """
        return int(get_gaussian_peakarea_on_intensity_list(track_intensity, left_base, right_base))


    def get_DIMS_feature_table(self):
        """
        Mock DIMS feature as max intensity per mass track.
        """

        self.FeatureList = []
        for ii in self.MassGrid.index:
            self.FeatureList.append({
                'id_number': 'F'+str(ii),
                'parent_masstrack_id': ii,
                'mz': self.MassGrid['mz'][ii],
                'intensity': 0,
                'apex': 0,
                'left_base': 0,
                'right_base': 0,
                'rtime': 0,
                'rtime_left_base': 0,
                'rtime_right_base': 0,
                'peak_area': 0,
                'cSelectivity': 0,
                'goodness_fitting': 0,
                'snr': 0,
            })

        FeatureTable = pd.DataFrame(self.FeatureList)
        for sample in self.experiment.all_samples:
            ## if not self.experiment.parameters['drop_unaligned_samples'] or SM.is_rt_aligned:
            fList = []
            list_mass_tracks = sample.get_masstracks_and_anchors()
            for peak in self.FeatureList:
                track_number = self.MassGrid[sample.name][peak['parent_masstrack_id']]
                if not pd.isna(track_number):
                    mass_track = list_mass_tracks[int(track_number)]
                    # left_base = None
                    # right_base = None
                    peak_area = max(mass_track['intensity'])
                else:
                    peak_area = 0
                fList.append( peak_area )
            FeatureTable[sample.name] = fList
        # print("\nFeature Table: ", FeatureTable.shape)
        self.FeatureTable = FeatureTable


    def generate_feature_table(self):
        """
        Initiate and populate self.FeatureTable, each sample per column in dataframe.
        """
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
        # print("\nFeature Table: ", FeatureTable.shape)
        self.FeatureTable = FeatureTable


    def extract_features_per_sample(self, sample, peak_area_function):
        """
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
        """
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
        reference_path = os.path.join(self.experiment.parameters['outdir'], reference_sample_name + '_mz_rtime_landmarks.csv')
        with open(reference_path, 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(["mz", "rtime"])  #  headers
            writer.writerows(zip(mz_landmarks, rtime_landmarks))
