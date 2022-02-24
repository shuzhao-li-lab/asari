'''
Jie Young's solution at StackOverflow:
https://stackoverflow.com/questions/50916422/python-typeerror-object-of-type-int64-is-not-json-serializable
'''

import json
import numpy as np

class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NpEncoder, self).default(obj)
