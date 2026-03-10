import os
import subprocess
import multiprocessing as mp

from zipfile import ZipFile

import requests

from .utils import bulk_process


class mzMLconverter:
    versions = {
        "1.4.5": "https://github.com/compomics/ThermoRawFileParser/releases/download/v1.4.5/ThermoRawFileParser1.4.5.zip"
    }
    
    installed = {
    }

    def __init__(self, version="1.4.5", dask_ip=False, multicores=None):
        self.version = version
        self.dask_ip = dask_ip
        self.multicores = multicores if multicores is not None else mp.cpu_count()
        self.install_converter(version)
        self.command_template = self.__conversion_command_template()

    def __conversion_command_template(self):
        if self.version in mzMLconverter.installed:
            converter_path = os.path.join(mzMLconverter.installed[self.version], "ThermoRawFileParser.exe")
        else:
            converter_path = os.path.join(mzMLconverter.install_converter(self.version), "ThermoRawFileParser.exe")

        def __determine_executable():
            engines = [["mono"], []]
            for engine in engines:
                try:
                    subprocess.run(engine + [converter_path], capture_output=True)
                    return engine
                except:
                    pass
            raise Exception("Could not determine the executable to run the converter!")

        engine = __determine_executable()
        if engine is not None:
            command_template = engine + [converter_path, "-i", "INPUT", "-o", "OUTPUT", "-f", "2"]
        else:
            raise Exception("Could not determine the executable to run the converter")
        return command_template

    @staticmethod
    def install_converter(version="1.4.5"):
        if mzMLconverter.installed.get(version, None) is None:
            converter_release_url = mzMLconverter.versions[version]
            local_zip_path = os.path.join(os.path.dirname(__file__), f"ThermoRawFileParser{version}.zip")
            extract_path = os.path.join(os.path.dirname(__file__), f"ThermoRawFileParser{version}")

            if os.path.exists(extract_path):
                # Already downloaded and extracted
                mzMLconverter.installed[version] = os.path.abspath(extract_path)
                return os.path.abspath(extract_path)
            else:
                # Download the zip file
                try:
                    response = requests.get(converter_release_url)
                    with open(local_zip_path, 'wb') as f:
                        f.write(response.content)

                    # Extract the zip file
                    with ZipFile(local_zip_path, 'r') as zip_ref:
                        zip_ref.extractall(extract_path)

                    # Clean up the zip file
                    os.remove(local_zip_path)
                    if os.path.exists(extract_path):
                        mzMLconverter.installed[version] = os.path.abspath(extract_path)
                        return os.path.abspath(extract_path)
                    else:
                        raise Exception("Failed to install the converter")
                except Exception as e:
                    if os.path.exists(local_zip_path):
                        os.remove(local_zip_path)
        else:
            return mzMLconverter.installed[version] 
        
    @staticmethod
    def uninstall_converter(version="1.4.5"):
        if mzMLconverter.installed.get(version, False):
            os.remove(mzMLconverter.installed[version])
            mzMLconverter.installed.pop(version)
    
    def bulk_convert(self, raw_files):
        conversion_commands = [" ".join([x.replace("INPUT", file) for x in self.command_template]) for file in raw_files]
        conversion_commands = [x.replace("OUTPUT", os.path.dirname(file)) for (x, file) in zip(conversion_commands, raw_files)]
        results = bulk_process(os.system, conversion_commands, self.dask_ip)
        for (converted, to_convert) in zip(results, raw_files):
            if converted == 0:
                print(f"Successfully converted {to_convert}")
            else:
                print(f"Failed to convert {to_convert}")