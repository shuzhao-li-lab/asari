import os
from matchms import Spectrum
from matchms.exporting import save_as_msp
from matchms.importing import load_spectra

def combine_and_convert(dir_path, output_file):
    # Initialize an empty list to store all spectra
    combined_spectra = []

    # Iterate over each file in the directory
    for filename in os.listdir(dir_path):
        if filename.lower().endswith(".msp"):  # Assuming all files are .msp format
            for spectrum in load_spectra(os.path.join(dir_path, filename)):
                combined_spectra.append(spectrum)

    # Save as mzML format
    save_as_msp(combined_spectra, output_file)
    print(combined_spectra)

import sys
combine_and_convert(sys.argv[1], sys.argv[2])
