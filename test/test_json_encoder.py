import json
import numpy as np
import unittest
from asari.json_encoder import NpEncoder

class TestNpEncoder(unittest.TestCase):
    def test_encode_np_integer(self):
        data = {'value': np.int64(42)}
        encoded = json.dumps(data, cls=NpEncoder)
        self.assertEqual(encoded, '{"value": 42}')

    def test_encode_np_floating(self):
        data = {'value': np.float64(3.14)}
        encoded = json.dumps(data, cls=NpEncoder)
        self.assertEqual(encoded, '{"value": 3.14}')

    def test_encode_np_array(self):
        data = {'values': np.array([1, 2, 3])}
        encoded = json.dumps(data, cls=NpEncoder)
        self.assertEqual(encoded, '{"values": [1, 2, 3]}')

    def test_encode_mixed_data(self):
        data = {
            'int': np.int64(42),
            'float': np.float64(3.14),
            'array': np.array([1, 2, 3]),
            'string': 'test'
        }
        encoded = json.dumps(data, cls=NpEncoder)
        self.assertEqual(encoded, '{"int": 42, "float": 3.14, "array": [1, 2, 3], "string": "test"}')

    def test_encode_nested_data(self):
        data = {
            'nested': {
                'int': np.int64(42),
                'float': np.float64(3.14),
                'array': np.array([1, 2, 3])
            }
        }
        encoded = json.dumps(data, cls=NpEncoder)
        self.assertEqual(encoded, '{"nested": {"int": 42, "float": 3.14, "array": [1, 2, 3]}}')

if __name__ == '__main__':
    unittest.main()