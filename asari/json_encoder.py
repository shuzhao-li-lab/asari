import json
import numpy as np

class NpEncoder(json.JSONEncoder):
    '''
    To handle numpy data types in JSON, using
    Jie Young's solution at StackOverflow:
    https://stackoverflow.com/questions/50916422/python-typeerror-object-of-type-int64-is-not-json-serializable
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
