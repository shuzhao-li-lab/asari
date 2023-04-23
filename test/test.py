'''
Test functions.
Example mzML data can be found from https://github.com/shuzhao-li/data.


>>> indir = 'SZ22_Dataset'
>>>
>>> list_input_files = read_project_dir(indir)
Working on ~~ SZ22_Dataset ~~

'''

from asari.main import *

print(PARAMETERS)


# get mass tracks from a single mzML file
def extract_single(infile):
    """
    extract_massTracks_ function returns 
    {
        'rt_numbers': rt_numbers,
        'rt_times': rt_times,
        'tracks': updated_tracks,
    }
    """
    exp = pymzml.run.Reader(infile)
    xdict = extract_massTracks_(exp,)
    print(
        xdict['rt_numbers'],
        xdict['rt_times'],
        xdict['tracks'][55]
    )

