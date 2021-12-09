'''
LC-MS metabolomics data pre-processing

Example use
-----------
python3 -m asari.main neg /Users/shuzhao/li.projects/asari/T03

'''

import sys
from .lcms_experiment import *

def main(directory):
    print("\n\n~~~~~~~ Hello from Asari! ~~~~~~~~~\n")
    process_project(
            read_project_dir(directory), {}, PARAMETERS, directory   #setting output_dir as input dir
    )

#
# -----------------------------------------------------------------------------
#
if __name__ == '__main__':
    PARAMETERS['mode'] = sys.argv[1]
    directory = sys.argv[2]
    main(directory)
