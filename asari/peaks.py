'''
Mass traces (i.e. XIC, EIC or chromatogram) and peaks.
'''

import numpy as np
from scipy.signal import find_peaks 
from scipy.optimize import curve_fit 

from metDataModel.core import MassTrace, Peak

from .utils import *


# peak detection is in this class
class ext_MassTrace(MassTrace):
    '''
    Extending metDataModel.core.MassTrace
    Peak detection using scipy.signal.find_peaks, a local maxima method with prominence control.
    Not keeping Peaks in this class; Peaks are attached to Sample. 
    To-do: further refine peak detection in future, 
    e.g. extra iteration of peak detection in remaining region; better deconvolution of overlapping peaks.
    '''
    def __init2__(self, mz, RT, INTS):
        # np.array or list -
        self.mz, self.list_retention_time, self.list_intensity = [mz, RT, INTS]
        # self.cal_list_retention_time = []             # not using now
        self.formula_mass = None
        self.raw_mzs = []                               # multiple values possible when merging traces
        self.__mass_corrected_by_asari__ = False        # True if updated by calibration
        self.features_assigned = False                  # tracking if assigned to Experiment features
        self.sample_name = ''

    def detect_peaks(self, min_intensity_threshold, min_fwhm, min_prominence_threshold, prominence_window, gaussian_shape, snr):
        list_peaks = []
        peaks, properties = find_peaks(self.list_intensity, height=min_intensity_threshold, width=min_fwhm, 
                                                        prominence=min_prominence_threshold, wlen=prominence_window) 
        # This is noise estimate, but doesn't have to correspond to final list of peaks
        _noise_level_ = self.__get_noise_level__(peaks, properties)
        if peaks.size > 1:
            # Rerun, raising min_prominence_threshold. Not relying on peak shape for small peaks, 
            # because chromatography is often not good so that small peaks can't be separated from noise.
            _tenth_height = 0.1* max(self.list_intensity)
            min_prominence_threshold = max(min_prominence_threshold, _tenth_height)
            peaks, properties = find_peaks(self.list_intensity, height=min_intensity_threshold, width=min_fwhm, 
                                                        prominence=min_prominence_threshold, wlen=prominence_window)
        
        for ii in range(peaks.size):
            if properties['peak_heights'][ii] > snr*_noise_level_:
                P = ext_Peak()
                # scipy.signal.find_peaks works on 1-D array. RT coordinates need to be added back.
                P.__init2__(parent_mass_trace=self, mz = self.mz, apex = peaks[ii],
                                peak_height = properties['peak_heights'][ii], left_base = properties['left_bases'][ii],
                                right_base = properties['right_bases'][ii])
                P.evaluate_peak_model()
                if P.goodness_fitting > gaussian_shape:
                    list_peaks.append(P)

        return list_peaks

    def extract_targeted_peak(self, rt_range):
        pass

    def __get_noise_level__(self, peaks, properties):
        peak_data_points = []
        for ii in range(peaks.size):
            peak_data_points += range(properties['left_bases'][ii], properties['right_bases'][ii]+1)
        noise_data_points = [ii for ii in range(len(self.list_intensity)) if ii not in peak_data_points]
        if noise_data_points:
            return np.median([self.list_intensity[ii] for ii in noise_data_points])
        else:
            return 0

# peak evaluation is in this class
class ext_Peak(Peak):
    '''
    Extending metDataModel.core.Peak, including pointer to parent MassTrace.
    Not storing RT or intensity arrays, but looking up in MassTrace if needed.
    '''
    def __init2__(self, parent_mass_trace, mz, apex, peak_height, left_base, right_base):
        [ self.parent_mass_trace, 
        self.mz, self.apex, self.peak_height, self.left_base, self.right_base
                ] = [ parent_mass_trace, mz, apex, peak_height, left_base, right_base ]


        # need to convert back to seconds in asari xics

        self.left_rtime = float(self.parent_mass_trace.list_retention_time[left_base])
        self.right_rtime = float(self.parent_mass_trace.list_retention_time[right_base])

    def evaluate_peak_model(self):
        '''
        Use Gaussian models to fit peaks, R^2 as goodness of fitting.
        Peak area is defined by Gaussian model, the integral of a Gaussian function being a * c *sqrt(2*pi).
        Good: Left to right base may not be full represenation of the peak. The Gaussian function will propose a good boundary.
        Less good: peaks are not always in Gaussian shape. But we are comparing same features across samples, 
        same bias in peak shape is applied to all samples.
        '''
        xx = self.parent_mass_trace.list_retention_time[self.left_base: self.right_base+1]
        # set initial parameters
        a, mu, sigma =  self.peak_height, \
                        self.parent_mass_trace.list_retention_time[self.apex], \
                        np.std(xx)
        try:
            popt, pcov = curve_fit(gaussian_function__, 
                            xx, self.parent_mass_trace.list_intensity[self.left_base: self.right_base+1],
                            p0=[a, mu, sigma])
            self.gaussian_parameters = popt
            self.rtime = popt[1]
            self.peak_area = popt[0]*abs(popt[2])*2.506628274631        # abs because curve_fit may return negative sigma
            self.goodness_fitting = goodness_fitting__(
                            self.parent_mass_trace.list_intensity[self.left_base: self.right_base+1], 
                            gaussian_function__(xx, *popt))

        # failure to fit
        except RuntimeError:
            self.rtime = mu
            self.peak_area = a*sigma*2.506628274631
            self.goodness_fitting = 0

    def extend_model_range(self):
        # extend model xrange, as the initial peak definition may not be complete, for visualization
        _extended = self.right_base - self.left_base
        # negative index does not work here, thus max to 0
        self.rt_extended = self.parent_mass_trace.list_retention_time[
                                            max(0, self.apex-_extended): self.apex+_extended]
        self.y_fitted_extended = gaussian_function__(self.rt_extended, *self.gaussian_parameters)
