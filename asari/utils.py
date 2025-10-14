# Asari utils will be a catchall for reused code that doesn't fit into a specific class
# this should be importable in other code in related tools. This should be cleaner than
# having a dedicated util file in each tool and honeslty, if you have the other tools 
# installed, you probalby have asari too. 

import multiprocessing as mp
import os
import time
import hashlib
import pymzml
import tqdm

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

def _process_with_dask(command, arguments, scheduler_ip):
    """Processes arguments in parallel using a Dask cluster."""
    try:
        from dask.distributed import Client, progress
    except ImportError:
        raise ImportError(
            "Dask must be installed to use this feature. Run 'pip install dask distributed'."
        )

    print(f"Connecting to Dask scheduler at: tcp://{scheduler_ip}:8786")
    with Client(f"tcp://{scheduler_ip}:8786") as client:
        # client.map is the most straightforward and idiomatic approach.
        # It submits all tasks at once, and Dask's scheduler manages the queue
        # efficiently. Most importantly, it returns results in the same order
        # as the input arguments, eliminating the need for manual sorting.
        futures = client.map(command, arguments)
        
        # Use Dask's built-in progress bar for a rich console display.
        print("All tasks submitted to the cluster. Waiting for completion...")
        progress(futures, notebook=False)
        
        # client.gather collects results once all futures are complete.
        results = client.gather(futures)
        print("Processing complete.")

    return results

def _process_with_multiprocessing(command, arguments, num_workers=None):
    """Processes arguments in parallel using a local multiprocessing pool."""
    # Default to the number of CPU cores if num_workers is not specified.
    pool_size = num_workers if isinstance(num_workers, int) and num_workers > 0 else mp.cpu_count()
    print(f"Starting a local multiprocessing pool with {pool_size} workers.")

    with mp.Pool(pool_size) as pool:
        # pool.imap is memory-efficient and returns an ordered iterator.
        # We wrap it with tqdm for a clean, real-time progress bar.
        pbar = tqdm.tqdm(pool.imap(command, arguments), total=len(arguments), desc="Processing")
        # Collect the results from the progress bar iterator.
        results = list(pbar)
        
    return results

def _process_with_gen_multiprocessing(command, arguments, num_workers=None):
    """Processes arguments in parallel using a local multiprocessing pool."""
    # Default to the number of CPU cores if num_workers is not specified.
    pool_size = num_workers if isinstance(num_workers, int) and num_workers > 0 else mp.cpu_count()
    print(f"Starting a local multiprocessing pool with {pool_size} workers.")

    with mp.Pool(pool_size) as pool:
        # pool.imap is memory-efficient and returns an ordered iterator.
        # We wrap it with tqdm for a clean, real-time progress bar.
        pbar = tqdm.tqdm(pool.imap(command, arguments), total=len(arguments), desc="Processing")
        # Collect the results from the progress bar iterator.
        yield from pbar

def bulk_process(command, arguments, dask_ip=None, num_workers=None, streaming=False):
    """
    Executes a command in parallel for each item in a list of arguments.

    This function acts as a dispatcher, choosing between a distributed Dask cluster
    (if `dask_ip` is provided) or a local `multiprocessing` pool.

    Args:
        command (callable): The function to execute for each argument.
        arguments (list): A list of arguments to be passed to the command.
        dask_ip (str, optional): The IP address of the Dask scheduler. If provided,
                                 Dask will be used for processing. Defaults to None.
        num_workers (int, optional): The number of worker processes for local 
                                     multiprocessing. Defaults to the number of
                                     CPU cores. This parameter is ignored if 
                                     `dask_ip` is set.

    Returns:
        list: A list containing the results in the same order as the input arguments.
    """
    if not arguments:
        raise ValueError("The 'arguments' list cannot be empty.")

    if dask_ip:
        # Use Dask for distributed processing.
        return _process_with_dask(command, arguments, dask_ip)
    elif streaming:
        raise NotImplementedError()
        return _process_with_gen_multiprocessing(command, arguments, num_workers)
    else:
        # Use local multiprocessing.
        return _process_with_multiprocessing(command, arguments, num_workers)
    