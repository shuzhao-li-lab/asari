'''
Do peaks export per sample in a directory:
shuzhao@red5:~/projects/asari$ python3 -m asari.main /home/shuzhao/projects/pipelineMJ/preprocessing/test_R01batch9/


mydir = '/home/shuzhao/projects/pipelineMJ/preprocessing/test_R01batch9/'   

In [7]: [X.number_of_peaks for X in SM1.list_MassTraces[10: 30]]                                                      
Out[7]: [1, 2, 1, 1, 1, 1, 1, 3, 1, 1, 2, 1, 1, 1, 2, 2, 1, 2, 1, 3]

In [21]: for ii in range( 7):  
    ...:     print(SM1.list_MassTraces[ii].mz, SM1.list_MassTraces[ii].number_of_peaks) 
    ...:     for P in SM1.list_MassTraces[ii].list_peaks:  
    ...:         print(P.goodness_fitting, P.peak_height, P.rtime) 


def __rt_overlap__(RT1, RT2):   # over 90% of one RT in the other
    S1, S2 = set([str(x)[:6] for x in RT1]), set([str(x)[:6] for x in RT2])
    if len(S1.intersection(S2)) > 0.9 * min(RT1.size, RT2.size):
        return True
    else:
        return False


In [8]: len(SM1.list_MassTraces)                                                                
Out[8]: 3283

In [9]: good[0]                                                                                 
Out[9]: 
['C4H7N2_83.060374',
 83.06037446677,
 'C4H7N2',
 0.9965129342851303,
 [('C4H6N2_82.053098', 'M+H[1+]')]]

In [10]: good[1000]                                                                             
Out[10]: 
['C10H9ClNNaO3_250.024119',
 250.02411946677,
 'C10H9ClNNaO3',
 None,
 [('C10H9NO3_191.058243', 'M+NaCl[1+]')]]
 
In [11]: SM1.__mass_accuracy__                                                                  
Out[11]: -0.3369786631057757

In [12]: SM1.__mass_stdev__                                                                     
Out[12]: 3.766422552049111

'''

import matplotlib.pyplot as plt

from .algorithms import *
from .plot import plot_peaks_masstrace

f = '/home/shuzhao/projects/pipelineMJ/preprocessing/test_R01batch9/MT_20210803_115_chrom.mzML'

SM1 = Sample(f)
SM1.process_step_1()
SM1.export_peaklist('MT_20210803_115.peaklist')

plot_peaks_masstrace(SM1.list_MassTraces[4], 'masstrace4.pdf')

SM1._match_mass_formula_(DB_1)

good = [M.formula_mass for M in SM1.list_MassTraces if M.formula_mass]


