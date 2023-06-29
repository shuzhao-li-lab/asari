__version__ = "1.11.6"

import sys
import importlib
import pkgutil

def import_submodules(package_name):
    '''
    Import all submodules of a module, recursively

    Parameters
    ----------
    package_name : str
        the name of a package

    Outputs
    -------
    dict[types.ModuleType]
    '''
    package = sys.modules[package_name]
    return {
        name: importlib.import_module(package_name + '.' + name)
        for loader, name, is_pkg in pkgutil.walk_packages(package.__path__)
    }
__all__ = import_submodules(__name__).keys()