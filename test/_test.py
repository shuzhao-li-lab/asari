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

# ad hoc, will repalce with better lists for calibration, for each of pos and neg.
calibration_mass_dict_pos = {'C4H9N3O2_131.069477': 132.07675346677001, 'C6H11NO2_129.078979': 130.08625546677, 'C5H9NO4_147.053158': 148.06043446677, 'C6H10O6_178.047738': 179.05501446677002, 'C9H11NO3_181.073893': 182.08116946677, 'C9H11NO2_165.078979': 166.08625546677, 'C9H8O3_164.047344': 165.05462046677002, 'C23H45NO4_399.334859': 400.34213546677, 'C2H7NO3S_125.014664': 126.02194046676999, 'C5H7NO3_129.042593': 130.04986946677002, 'C5H4N4O3_168.02834': 169.03561646677, 'C5H8O5_148.037173': 149.04444946677, 'C5H12O5_152.068473': 153.07574946677002, 'C9H8O2_148.052429': 149.05970546677, 'C5H10N2O3_146.069142': 147.07641846677, 'C17H33NO4_315.240959': 316.24823546677, 'C5H11NO2S_149.051049': 150.05832546677001, 'C5H10O2_102.06808': 103.07535646676999, 'C15H29NO4_287.209658': 288.21693446677, 'C8H17NO2_159.125929': 160.13320546677002, 'H2O4S_97.967379': 98.97465546676999, 'C6H13NO5_179.079373': 180.08664946677, 'C16H22O4_278.151809': 279.15908546677, 'C20H28O2_300.20893': 301.21620646677, 'C7H8N4O2_180.064726': 181.07200246677002, 'C7H6O2_122.036779': 123.04405546676999, 'C14H22N2O3_266.163043': 267.17031946677, 'C19H16O4_308.104859': 309.11213546677, 'C21H39NO4_369.287909': 370.29518546677, 'C8H6O4_166.026609': 167.03388546677002, 'C18H35NO_281.271865': 282.27914146677, 'C19H37NO4_343.272259': 344.27953546677, 'C26H52NO7P_521.34814': 522.35541646677, 'C7H8N2O2_152.058578': 153.06585446677002, 'C9H17NO2_171.125929': 172.13320546677002, 'C25H47NO4_425.350509': 426.35778546677, 'C8H10O3_154.062994': 155.07027046677, 'C25H45NO4_423.334859': 424.34213546677, 'C22H46NO7P_467.301189': 468.30846546677003, 'C23H48NO7P_481.316839': 482.32411546677, 'C24H48NO7P_493.316839': 494.32411546677, 'C26H54NO7P_523.36379': 524.37106646677, 'C26H48NO7P_517.316839': 518.32411546677, 'C28H52NO7P_545.34814': 546.35541646677, 'C28H50NO7P_543.332489': 544.33976546677, 'C24H50NO6P_479.337575': 480.34485146677, 'C25H52NO7P_509.34814': 510.35541646677, 'C21H38O4_354.27701': 355.28428646677, 'C17H31NO4_313.225308': 314.23258446677, 'C15H27NO4_285.194008': 286.20128446677, 'C19H35NO4_341.256609': 342.26388546677003, 'C22H37NO3_363.277344': 364.28462046677004, 'C19H22ClN5O_371.151288': 372.15856446677003, 'C24H30O6_414.204239': 415.21151546677, 'C21H25N_291.1987': 292.20597646677, 'C15H23NO2_249.172879': 250.18015546677, 
    'C17H27NO4_309.194008': 310.20128446677, 'C21H25ClN2O3_388.15537': 389.16264646677, 'C16H18O10_370.089997': 371.09727346677}


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


