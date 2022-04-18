'''
dated.


✗ pip3 install asari-metabolomics

➜ time python3 extract_masstracks.py T04
Determination of memory status is not supported on this 
 platform, measuring for memoryleaks will never fail
Processing MG_20211022_011.mzML, found 4579 mass tracks.
Processing MG_20211022_007.mzML, found 4338 mass tracks.
Processing MG_20211022_015.mzML, found 5150 mass tracks.
Processing MG_20211022_005.mzML, found 4897 mass tracks.
Processing MG_20211022_013.mzML, found 4558 mass tracks.
Processing MG_20211022_009.mzML, found 4951 mass tracks.
python3 extract_masstracks.py T04  25.07s user 0.67s system 106% cpu 24.153 total

'''

import os
import numpy as np
from asari.samples import SimpleSample


def read_project_dir(directory, file_pattern='.mzML'):
    '''
    This reads centroided LC-MS files.
    '''
    return [os.path.join(directory, f) for f in os.listdir(directory) if file_pattern in f]


def process_single_file(infile, outdir=''):
    '''
    Write out
    id_number, mz, number_scans, max_intensity, median_intensity
    '''
    SS = SimpleSample(input_file = infile)
    SS.get_mass_tracks_(mz_tolerance_ppm=5, min_intensity=100, min_timepoints=5)
    s = '\t'.join(['id_number', 'mz', 'number_scans', 'max_intensity', 'median_intensity']) + '\n'
    for track in SS.list_mass_tracks:
        LL = np.array(track['intensity'])
        s += '\t'.join(
            [str(x) for x in [track['id_number'], track['mz'], 
                        LL.size, LL.max(), int(np.median(LL[LL>0]))]]
        ) + '\n'
    
    outfile = os.path.join(outdir, os.path.basename(infile).replace(".mzML", "_massTrack.tsv"))
    with open(outfile, 'w') as O:
        O.write(s)

#
# -----------------------------------------------------------------------------
#
if __name__ == '__main__':
    import sys
    directory = sys.argv[1]
    myfiles = read_project_dir(directory)
    for f in myfiles:
        process_single_file(f, outdir=directory)    # export to directory

