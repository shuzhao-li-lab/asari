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

import numpy as np
from matplotlib import pyplot as plt
import pymzml

def get_plot_region_from_file(infile, min_scan_number, max_scan_number, min_mz, max_mz, ms_level=1):
    '''
    input
    -----
    infile: mzML file as input

    return 
    ------
    list, [(scan_number, mz, intensity value), ...]
    '''
    alldata = []
    exp = pymzml.run.Reader(infile)
    ii = 0   # scan_number starts with 0
    for spec in exp:
        if min_scan_number < ii < max_scan_number:
            if spec.ms_level == ms_level:
                _NN = spec.mz.size
                for jj in range(_NN):
                    if min_mz < spec.mz[jj] < max_mz:
                        alldata.append((ii, spec.mz[jj], int(spec.i[jj])))
                
        ii += 1
    return alldata


def plot_scatter_map_region(datapoints, figsize=(8,10), cmap=plt.cm.coolwarm, colorbar_orientation='horizontal'):
    '''
    Make a color map of LC-MS data, each dot a data point
    '''
    datapoints = np.array(datapoints).T
    normZ = np.log2(datapoints[2]+1)
    normZ = normZ / normZ.max()
    fig, ax = plt.subplots(figsize=figsize)
    im = ax.scatter(datapoints[0], datapoints[1], marker='.', c=normZ, cmap=cmap)
    fig.colorbar(im, ax=ax, orientation=colorbar_orientation, shrink=0.3, ticks=[0], label='relative intensity')
    #im.set_clim(0.0, 1.0)  # set the color limits 


def double_scatter_map_region(datapoints, figsize=(8,10), cmap=plt.cm.coolwarm, colorbar_orientation='horizontal'):
    datapoints = np.array(datapoints).T
    normZ = np.log2(datapoints[2]+1)
    normZ = normZ / normZ.max()
    fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=figsize)
    im = ax1.scatter(datapoints[0], datapoints[1], marker='.', c=normZ, cmap=cmap)
    ax1.set(ylabel='m/z')
    im2 = ax2.scatter(datapoints[0], datapoints[2], marker='.', c=normZ, cmap=cmap)
    ax2.set(ylabel='intensity')
    fig.colorbar(im2, ax=ax2, orientation=colorbar_orientation, shrink=0.3, ticks=[0], label='')
    

def with_line_scatter_map_region(datapoints, figsize=(8,10), cmap=plt.cm.coolwarm):
    datapoints = np.array(datapoints).T
    normZ = np.log2(datapoints[2]+1)
    normZ = normZ / normZ.max()
    fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=figsize)
    im = ax1.scatter(datapoints[0], datapoints[1], marker='.', c=normZ, cmap=cmap)
    ax1.set(ylabel='m/z')
    im2 = ax2.plot(datapoints[0], datapoints[2], '*-')
    ax2.set(ylabel='intensity')




# -----------------------------------------------------------------------------
# in progress

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




