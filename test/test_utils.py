import unittest
import os
import tempfile
import time
import hashlib
from asari.utils import sizeof_fmt, checksum_file, wait_with_pbar, get_ionization_mode_mzml, bulk_process

class TestUtils(unittest.TestCase):

    def test_sizeof_fmt(self):
        self.assertEqual(sizeof_fmt(1024), "1.0KiB")
        self.assertEqual(sizeof_fmt(1048576), "1.0MiB")
        self.assertEqual(sizeof_fmt(1073741824), "1.0GiB")

    def test_checksum_file(self):
        with tempfile.NamedTemporaryFile(delete=False) as tmp_file:
            tmp_file.write(b"test content")
            tmp_file_path = tmp_file.name
        self.assertEqual(checksum_file(tmp_file_path), hashlib.md5(b"test content").hexdigest())
        os.remove(tmp_file_path)

    def test_wait_with_pbar(self):
        start_time = time.time()
        wait_with_pbar(wait=2)
        end_time = time.time()
        self.assertTrue(1.9 <= end_time - start_time <= 2.1)

    def test_get_ionization_mode_mzml(self):
        # todo - this needs to be implemented
        pass

    @staticmethod
    def dummy_command(x):
        return x * 2

    def test_bulk_process(self):
        arguments = [1, 2, 3, 4, 5]
        results = bulk_process(self.dummy_command, arguments)
        self.assertEqual(results, [2, 4, 6, 8, 10])

if __name__ == '__main__':
    unittest.main()