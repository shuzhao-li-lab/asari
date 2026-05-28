# utility code that doesn't fit into a specific class

import multiprocessing as mp
import os
import time
import hashlib
import zipfile
from io import BytesIO
from importlib import resources as pkg_resources

import json
import requests

import numpy as np
import pymzml
import tqdm


class NpEncoder(json.JSONEncoder):
    '''
    To handle numpy data types in JSON, using
    Jie Young's solution at StackOverflow:
    https://stackoverflow.com/questions/50916422/python-typeerror-object-of-type-int64-is-not-json-serializable

    Will change to use json_tricks
    '''
    def default(self, obj):
        '''
        This function converts obj into something that can be serialized by JSON, largely for 
        handling numpy datatypes that, despite being largely equivalent to their pure python
        equivalents, cannot be serialized. 

        Parameters
        ----------
        obj: np.integer or np.floating or np.ndarray or other serializable object instance
            for the numpy objects above, they are cast to their python 'equivalents', i.e., 
            np.integer -> int, np.floating -> float, np.ndarray -> list, else, the object
            is converted to its default serialization representation.
        '''
        
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NpEncoder, self).default(obj)

def bulk_process(command, arguments, jobs_per_worker=False):
    '''
    Multiprocessing for pooling tasks. 
    
    Not considering distributed computing (e.g. dask) because it's easier to parallelize projects than computing. 
    '''
    n_workers = jobs_per_worker if (jobs_per_worker and jobs_per_worker != 'auto') else mp.cpu_count()
    with mp.Pool(n_workers) as client:
        return client.starmap(command, [(arg,) for arg in arguments])
    
    
def download_and_unzip_to_pkg_resources(url, package, subdir="data"):
    """Downloads a ZIP archive from a URL and extracts it into a package's resource directory."""
    print("HERE")
    # Get the package directory
    package_dir = os.path.dirname(pkg_resources.files(package))
    
    # Target extraction directory within the package
    extract_to = os.path.join(package_dir, subdir)
    os.makedirs(extract_to, exist_ok=True)  # Ensure the directory exists

    # Download the ZIP file
    response = requests.get(url, stream=True)
    response.raise_for_status()

    # Extract the ZIP archive to the target directory
    with zipfile.ZipFile(BytesIO(response.content)) as zip_ref:
        zip_ref.extractall(extract_to)

    print(f"Extracted to: {extract_to}")

def download_and_unzip(url, extract_to):
    """Downloads a ZIP archive from a URL and extracts it to the specified directory."""
    response = requests.get(url, stream=True)
    response.raise_for_status()  # Raise an error for bad responses

    with zipfile.ZipFile(BytesIO(response.content)) as zip_ref:
        zip_ref.extractall(extract_to)

    print(f"Extracted to: {extract_to}")

def validate_mzml_file(file):
    try:
        with pymzml.run.Reader(file) as reader:
            for spec in reader:
                pass
    except:
        return False
    return True

def build_boolean_dict():
    return {
        'T': True, 
        'F': False, 
        1: True, 
        0: False, 
        'True': True, 
        'False': False, 
        'TRUE': True, 
        'FALSE': False, 
        'true': True, 
        'false': False
    }

def sizeof_fmt(num, suffix="B"):
    for unit in ("", "Ki", "Mi", "Gi", "Ti", "Pi", "Ei", "Zi"):
        if abs(num) < 1024.0:
            return f"{num:3.1f}{unit}{suffix}"
        num /= 1024.0
    return f"{num:.1f}Yi{suffix}"

def checksum_file(file, chunksize=16384):
    assert os.path.isfile(file)
    hash_md5 = hashlib.md5()
    with open(file, 'rb') as f:
        for chunk in iter(lambda: f.read(chunksize), b''):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

def wait_with_pbar(wait=5):
    for _ in tqdm.tqdm(range(wait), total=wait, desc="waiting..."):
        time.sleep(1)

def get_ionization_mode_mzml(mzml_file, limit=50):
    ion_modes = set()
    i = 0
    with pymzml.run.Reader(mzml_file.path) as reader:
        for spec in reader:
            i += 1
            assert 'positive scan' in spec or 'negative scan' in spec
            if spec['positive scan']:
                ion_modes.add("pos")
            else:
                ion_modes.add("neg")
            if len(ion_modes) > 1:
                return "mixed"
            if i > limit:
                break
    return list(ion_modes)[0]
