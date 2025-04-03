
import requests
from io import BytesIO
import zipfile
import shutil
import os
import sys
import os
import subprocess
def install_asari_branch(branch=None):
    asari_patch = "https://github.com/shuzhao-li-lab/MatchMS_Asari/archive/refs/heads/main.zip"
    matchms_url = "https://github.com/matchms/matchms/archive/1e809be6ba8cdd75a5fc941858ec57d5956d8e65.zip"
    sparsetack_url = "https://github.com/matchms/sparsestack/archive/559ea0474c79cedb7f914809554ac39ba1e54fc1.zip"

    def download_and_unzip(url):
        response = requests.get(url)
        
        if response.status_code != 200:
            raise Exception('Failed to download file')
            
        zipped_source = zipfile.ZipFile(BytesIO(response.content))
        # Get the current working directory before extracting
        cwd = os.getcwd() 
        zipped_source.extractall()
        return cwd

    def copy_file(src, dest):
        try:
            shutil.copy2(src, dest)
            print("File copied successfully.")
        except FileNotFoundError:
            print("Source file not found.")
        except IsADirectoryError:
            print("Destination is a directory and it must be a file path.")
        except PermissionError:
            print("Permission denied.")
        except Exception as e:
            print(f"Error occurred while copying file. Error message: {str(e)}")
        return os.path.abspath(dest)

    patch_dir = download_and_unzip(asari_patch) + "/MatchMS_Asari-main/"
    matchms_dir = download_and_unzip(matchms_url) + "/matchms-1e809be6ba8cdd75a5fc941858ec57d5956d8e65/"
    sparsetack_dir = download_and_unzip(sparsetack_url) + "/sparsestack-559ea0474c79cedb7f914809554ac39ba1e54fc1/"

    # print the directories:

    print("patch_dir: f{patch_dir}")
    print("matchms_dir: f{matchms_dir}")
    print("sparsetack_dir: f{sparsetack_dir}")

    pathed_matchms_path = copy_file(os.path.join(patch_dir, "matchms_pyproject.toml"), os.path.join(matchms_dir, "pyproject.toml"))
    pathed_sparsetack_path = copy_file(os.path.join(patch_dir, "sparsestack_pyproject.toml"), os.path.join(sparsetack_dir, "pyproject.toml"))

    def install_asari(branch=None):
        if branch:
            try:
                env = os.environ.copy()
                env['PIP_REQUIRE_VIRTUALENV'] = 'false'
                    
                git_url = f"https://github.com/shuzhao-li-lab/asari"
                install_command = ["pip3", "install", "-e"]
                
                clone_result = subprocess.run(["git", "clone", git_url, branch], env=env)
                if clone_result.returncode == 0:
                    subprocess.run(install_command + [branch])
            except Exception as e:
                print("Error occurred while installing Asari with specified branch.")
                print('Error details - {}'.format(str(e))) 
        else:
            os.system("pip3 install git+https://github.com/shuzhao-li-lab/asari")

    def install_module(path):
        abs_path = os.path.abspath(path)
        print(f"Installing module at path: {abs_path}")
        try:
            # Change the current working directory to the module's directory
            # Run pip install command
            os.system(f"pip install {abs_path}")
        except Exception as e:
            print(f"failed to install! {path}")

    install_module(sparsetack_dir)
    install_module(matchms_dir)
    install_asari()


if __name__ == '__main__':
    # url = "https://github.com/shuzhao-li-lab/MatchMS_Asari"
    #if len(sys.argv) > 1:
    #    install_asari_branch(sys.argv[1])
    #else:
    install_asari_branch()

