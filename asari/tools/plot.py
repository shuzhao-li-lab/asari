'''

In [211]: plt.plot(Y, '.')                                                                                            
Out[211]: [<matplotlib.lines.Line2D at 0x7f30b8c33d30>]

In [212]: plt.plot(peaks, Y[peaks], 'x')                                                                              
Out[212]: [<matplotlib.lines.Line2D at 0x7f30b8c29bb0>]

In [213]: plt.vlines(x=peaks, ymin=Y[peaks] - properties["prominences"], ymax = Y[peaks], color = "C1")               
Out[213]: <matplotlib.collections.LineCollection at 0x7f30b8fb8df0>

In [214]: plt.hlines(y=properties["width_heights"], xmin=properties["left_ips"],xmax=properties["right_ips"], color = 
     ...: "C1")                                                                                                       
Out[214]: <matplotlib.collections.LineCollection at 0x7f30b8f0bb80>

In [215]: plt.show()                                                                                                  

color='green', marker='o', linestyle='dashed',
...      linewidth=2, markersize=12

'''

from matplotlib import pyplot as plt

def plot_peaks_masstrace(sample, mzstr, outfile='masstrace_plot.pdf'):
    '''
    To inspect how peak models fit the raw data. Mass traces and Peaks are indexed by mzstr in each Sample.
    A mass trace may correspond more than one peaks.

    Input
    -----
    Sample instance with detected Peak instances.
    mzstr has to be precise string that is used as a key in dictionaries dict_masstraces and dict_peaks.
    '''
    plt.figure()
    for mass_trace in sample.dict_masstraces[mzstr]:
        plt.plot(mass_trace.list_retention_time, mass_trace.list_intensity, marker='o', linewidth=0, markersize=1)
    
    dict_peaks = sample.create_peak_dict()
    for P in dict_peaks[mzstr]:
        P.extend_model_range()
        plt.plot(P.rt_extended, P.y_fitted_extended, color='red', alpha=0.5, linewidth=0.6)
    plt.title("mass trace " + str(round(mass_trace.mz, 6)))
    plt.savefig(outfile)
    plt.close()


def plot_peaks():

    '''
    
        # extend model xrange, as the initial peak definition may not be complete
        _extended = self.right_base - self.left_base
        self.rt_extended = self.parent_mass_trace.list_retention_time[self.apex-_extended: self.apex+_extended]
        self.y_fitted_extended = __gaussian_function__(self.rt_extended, *popt)
    '''

    pass

def plot_sample_rt_calibration(sample, outfile='rt_calibration.pdf'):
    plt.figure()
    plt.plot(*sample.__rt_calibration__data__, marker='o', linewidth=1, markersize=1)
    plt.title("rt_calibration")
    plt.savefig(sample.name+outfile)
    plt.close()
