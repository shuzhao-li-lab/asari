'''

In [7]: [X.number_of_peaks for X in SM1.list_MassTraces[10: 30]]                                                      
Out[7]: [1, 2, 1, 1, 1, 1, 1, 3, 1, 1, 2, 1, 1, 1, 2, 2, 1, 2, 1, 3]

In [21]: for ii in range( 7):  
    ...:     print(SM1.list_MassTraces[ii].mz, SM1.list_MassTraces[ii].number_of_peaks) 
    ...:     for P in SM1.list_MassTraces[ii].list_peaks:  
    ...:         print(P.goodness_fitting, P.peak_height, P.rtime) 

'''

import matplotlib.pyplot as plt

from .algorithms import *
from .plot import plot_peaks_masstrace

f = '/home/shuzhao/projects/pipelineMJ/preprocessing/test_R01batch9/MT_20210803_115_chrom.mzML'

SM1 = Sample(f)
SM1.detect_peaks()

plot_peaks_masstrace(SM1.list_MassTraces[4], 'masstrace4.pdf')

