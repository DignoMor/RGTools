
import numpy as np
import json

def str2bool(bool_str):
    '''
    Convert String to boolean value
    '''

    if not bool_str:
        return False

    if bool_str.upper() == "FALSE":
        return False

    if bool_str.upper() == "None":
        return False

    return True

def str2none(str_val):
    '''
    Convert string to None
    '''

    if str_val.upper() == "NONE":
        return None

    return str_val

class NumpyEncoder(json.JSONEncoder):
    '''
    Helper class to encode data with numpy arrays for json serialization.
    '''
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, (np.integer, np.floating)):
            return float(obj)  # Convert NumPy scalars to Python float
        return super().default(obj)
