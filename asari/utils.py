# Asari utils will be a catchall for reused code that doesn't fit into a specific class
# this should be importable in other code in related tools. This should be cleaner than
# having a dedicated util file in each tool and honeslty, if you have the other tools 
# installed, you probalby have asari too. 

import multiprocessing as mp
import tqdm
import os
import time
import hashlib
import pymzml
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def plot_correlations(feature_table_path):
    ft = pd.read_csv(feature_table_path, sep='\t')
    ft["sum_intensity_max"] = np.log2(ft.iloc[:, 11:]+1).max(axis=1)
    ft["sum_intensity_mean"] = np.log2(ft.iloc[:, 11:]+1).mean(axis=1)
    ft["sum_intensity_median"] = np.log2(ft.iloc[:, 11:]+1).median(axis=1)
    ft["sum_intensity_min"] = np.log2(ft.iloc[:, 11:]+1).min(axis=1)
    ft["peak_area"] = np.log2(ft["peak_area"]+1)

    fig, axs = plt.subplots(2, 2, figsize=(15, 10))

    # Max
    print("max")
    axs[0, 0].scatter(ft["sum_intensity_max"], ft["peak_area"])
    sns.regplot(x=ft["sum_intensity_max"], y=ft["peak_area"], scatter=False, color='blue', robust=True, ax=axs[0, 0])
    r2_value_max = round(np.corrcoef(ft['sum_intensity_max'], ft['peak_area'])[0, 1]**2, 2)
    axs[0, 0].set_title(f"max vs peak_area (R2 = {r2_value_max})")
    axs[0, 0].set_xlabel("max")
    axs[0, 0].set_ylabel("peak_area")

    # Mean
    print("mean")
    axs[0, 1].scatter(ft["sum_intensity_mean"], ft["peak_area"])
    sns.regplot(x=ft["sum_intensity_mean"], y=ft["peak_area"], scatter=False, color='blue', robust=True, ax=axs[0, 1])
    r2_value_mean = round(np.corrcoef(ft['sum_intensity_mean'], ft['peak_area'])[0, 1]**2, 2)
    axs[0, 1].set_title(f"mean vs peak_area (R2 = {r2_value_mean})")
    axs[0, 1].set_xlabel("mean")
    axs[0, 1].set_ylabel("peak_area")

    # Median
    print("median")
    axs[1, 0].scatter(ft["sum_intensity_median"], ft["peak_area"])
    sns.regplot(x=ft["sum_intensity_median"], y=ft["peak_area"], scatter=False, color='blue', robust=True, ax=axs[1, 0])
    r2_value_median = round(np.corrcoef(ft['sum_intensity_median'], ft['peak_area'])[0, 1]**2, 2)
    axs[1, 0].set_title(f"median vs peak_area (R2 = {r2_value_median})")
    axs[1, 0].set_xlabel("median")
    axs[1, 0].set_ylabel("peak_area")

    # Min
    print("min")
    axs[1, 1].scatter(ft["sum_intensity_min"], ft["peak_area"])
    sns.regplot(x=ft["sum_intensity_min"], y=ft["peak_area"], scatter=False, color='blue', robust=True, ax=axs[1, 1])
    r2_value_min = round(np.corrcoef(ft['sum_intensity_min'], ft['peak_area'])[0, 1]**2, 2)
    axs[1, 1].set_title(f"min vs peak_area (R2 = {r2_value_min})")
    axs[1, 1].set_xlabel("min")
    axs[1, 1].set_ylabel("peak_area")

    # Adjust layout
    plt.tight_layout()
    plt.show()

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

def bulk_process(command, arguments, dask_ip=None, jobs_per_worker=False, job_multiplier=1):
    DASK_ENABLE = False # to turn on Dask, set to True and install manually. EXPERIMENTAL!!!
    if arguments:
        if dask_ip:
            try:
                from dask.distributed import Client, as_completed, progress
            except:
                raise ImportError("Dask must be installed to use dask_ip=True")
            client = Client(f"tcp://{dask_ip}:8786")
            #client.scatter(arguments)
            #client.scatter(command)
            results = []
            if jobs_per_worker == 'auto':
                total_thread = sum(worker_info['nthreads'] for worker_info in client.scheduler_info()['workers'].values())
                total_workers = len(client.scheduler_info()['workers'])
                jobs_per_worker = total_thread // total_workers * job_multiplier
            if jobs_per_worker:
                pbar = tqdm.tqdm(total=len(arguments), desc="processing...")

                # Submit initial batch of jobs
                futures = {i: client.submit(command, arg) for i, arg in enumerate(arguments[:len(client.scheduler_info()['workers']) * int(jobs_per_worker)])}
                remaining_args = arguments[len(futures):]

                # Process jobs as they complete and submit new ones
                for i in range(len(arguments)):
                    completed_future = as_completed(futures.values()).next()
                    completed_index = list(futures.keys())[list(futures.values()).index(completed_future)]
                    pbar.update(1)
                    results.append((completed_index, completed_future.result()))
                    del futures[completed_index]
                    if remaining_args:
                        new_index = len(arguments) - len(remaining_args)
                        new_future = client.submit(command, remaining_args.pop(0))
                        futures[new_index] = new_future

                # Gather any remaining results
                if futures:
                    for future in futures.values():
                        pbar.update(1)
                        results.append((list(futures.keys())[list(futures.values()).index(future)], future.result()))

                # Sort results by original order
                results.sort(key=lambda x: x[0])
                return [result for _, result in results]
            else:
                # Default behavior: Submit all jobs at once
                futures = client.map(command, arguments)
                
                # Optionally show progress
                progress(futures)
                
                # Gather results once all tasks have completed
                return [x for x in client.gather(futures)]
        else:
            if jobs_per_worker and jobs_per_worker != 'auto':
                with mp.Pool(jobs_per_worker) as client:
                    pbar = tqdm.tqdm(client.imap(command, arguments), total=len(arguments), desc="processing...")
                    return [x for x in pbar]
            with mp.Pool(mp.cpu_count()) as client:
                pbar = tqdm.tqdm(client.imap(command, arguments), total=len(arguments), desc="processing...")
                return [x for x in pbar]
    else:
        raise Exception("No Arguments Provided")