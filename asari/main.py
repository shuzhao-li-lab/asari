'''
LC-MS metabolomics data pre-processing

- only for high resolution data
- performance conscious
- reproducible, trackable from features to XICs
- reference centril
- formula mass centric
- local maxima peak detection with prominence control

# will move to doc/
Typical chromatogram (XIC) extraction from mzML raw files:

> for file in ../mzML_IS_HILICposRPneg_05072021/*.mzML
>   do FeatureFinderMetabo -in $file -out ${file/.mzML/.featureXML} -out_chrom ${file/.mzML/_chrom.mzML} \
>   -algorithm:common:chrom_fwhm 1 -algorithm:mtd:mass_error_ppm 2 -algorithm:mtd:min_trace_length 2
> done

Parameters above: 1 sec, 2 ppm, 2 seconds
Ref: https://abibuilder.informatik.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/TOPP_FeatureFinderMetabo.html


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
