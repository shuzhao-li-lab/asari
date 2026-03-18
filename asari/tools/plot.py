'''
Misc plot functions.
'''

import numpy as np
from matplotlib import pyplot as plt
import pymzml

# -----------------------------------------------------------------------------
# Plot mass tracks and raw data points
# -----------------------------------------------------------------------------

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


def plot_masstrack(track, color='m', start=100, end=400, yticks=[0, 5e7, 1e8]):
    plt.figure(figsize=(8, 4))
    X = range(start, end)
    Y = track['intensity'][start: end]
    plt.plot(X, Y, color=color, marker='o', markersize=4, linestyle='dashed')
    plt.yticks(yticks)

# -----------------------------------------------------------------------------
# Mirror plot 
# -----------------------------------------------------------------------------

def mirror_plot(
    peaks, peaks2, 
    figsize=(8,4),
    label1="GCMS features",
    label2="Spectrum in lib",
    normalize=True,
    match_tol=None
):
    """
    Plot MS/MS mirror plot, modified from ChatGPT.

    Parameters
    ----------
    peaks : array-like
        m/z and intensity for spectrum 1 (top)
    peaks2 : array-like
        m/z and intensity for spectrum 2 (bottom)
    label1, label2 : str
        Spectrum labels
    normalize : bool
        Normalize intensities to max = 1
    match_tol : float or None
        If set, draw vertical lines connecting matched peaks within tolerance
    """
    mz1 = np.asarray([x[0] for x in peaks])
    int1 = np.asarray([x[1] for x in peaks])
    mz2 = np.asarray([x[0] for x in peaks2])
    int2 = np.asarray([x[1] for x in peaks2])
    if normalize:
        if int1.max() > 0:
            int1 = int1 / int1.max()
        if int2.max() > 0:
            int2 = int2 / int2.max()

    fig, ax = plt.subplots(figsize=figsize)
    # top spectrum
    ax.vlines(mz1, 0, int1, color="tab:blue", linewidth=1)
    # bottom spectrum (inverted)
    ax.vlines(mz2, 0, -int2, color="tab:red", linewidth=1)
    # optional peak matching
    if match_tol is not None:
        for m1, i1 in zip(mz1, int1):
            diff = np.abs(mz2 - m1)
            idx = np.where(diff <= match_tol)[0]
            for j in idx:
                ax.plot(
                    [m1, mz2[j]],
                    [i1, -int2[j]],
                    color="gray",
                    linewidth=0.5,
                    alpha=0.5
                )
    ax.axhline(0, color="black", linewidth=1)
    ax.set_xlabel("m/z")
    ax.set_ylabel("Normalized Intensity")
    ax.set_title("GCMS Mirror Plot")
    ax.text(0.01, 0.95, label1, transform=ax.transAxes,
            verticalalignment="top", color="tab:blue")
    ax.text(0.01, 0.05, label2, transform=ax.transAxes,
            verticalalignment="bottom", color="tab:red")
    plt.tight_layout()
    plt.show()

# -----------------------------------------------------------------------------
# Selectivity plots for m, c, d-selectivities
# -----------------------------------------------------------------------------

def plot_mSelectivity(mzList, selectivities, figsize=(10,3), save_pdf=True, 
                            outfile="Figure_mSelectivity"):
    '''
    mSelectivity plot, showing m/z as pink vertical lines, and selectivity as blue curve.
    Used same way for dSelectivity.
    '''
    plt.figure(figsize=figsize)
    plt.plot(mzList, selectivities, '-o', markersize=7)
    plt.vlines(mzList, 0, 1, color='r',linewidth=0.2)
    plt.xlabel("m/z")
    plt.ylabel("Selectivity")
    plt.title("m/z selectivity (zoom in)")
    if save_pdf:
        plt.savefig(outfile + ".pdf")


def plot_cSelectivity(mass_track, list_ranges=[(0, 800), (100, 150)], 
                            save_pdf=True, outfile="Figure_cSelectivity"):
    '''
    Plot chromatographic selectivity (cSelectivity).
    list_ranges: list of tuples. 1st tuple indicates the full data points, 
    each following tuples a highlighted region.
    '''
    a, b = list_ranges[0]
    X = range(a, b)
    Y = mass_track['intensity'][X]
    plt.figure(figsize=(10, 4))
    plt.plot(X, Y, color='k', marker='o', markersize=3, linestyle='None')
    for pair in list_ranges[1:]:
        a, b = pair
        X2 = range(a, b)
        plt.fill_between(X2, mass_track['intensity'][X2], color='R0', alpha=0.2,)
    if save_pdf:
        plt.savefig(outfile + ".pdf")


# -----------------------------------------------------------------------------
# in progress. Do not use.
# -----------------------------------------------------------------------------


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
    visualizing peaks:

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
