import asari
import sys


from asari.default_parameters import PARAMETERS
from asari.workflow import read_project_dir, process_project
from asari.main import update_peak_detection_params

project = read_project_dir(sys.argv[1])
PARAMETERS = update_peak_detection_params(PARAMETERS, None)
process_project(project, PARAMETERS)