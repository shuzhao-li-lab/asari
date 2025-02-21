import unittest
import os
import tempfile
import time
import hashlib
from asari.utils import sizeof_fmt, checksum_file, wait_with_pbar, get_ionization_mode_mzml, bulk_process
from collections import namedtuple
import requests
import zipfile
import io
import os

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

        os.makedirs("./Datasets", exist_ok=True)

        datasets = [
            "https://github.com/shuzhao-li-lab/data/raw/main/data/MT02Dataset.zip",
        ]

        for dataset in datasets:
            r = requests.get(dataset)
            if dataset.endswith(".zip"):
                z = zipfile.ZipFile(io.BytesIO(r.content))
                z.extractall("./Datasets/")
            else:
                with open("./Datasets/" + os.path.basename(dataset), 'bw+') as out_fh:
                    out_fh.write(r.content)
        
        mzml_file = namedtuple("mzml_file", ["path"])
        to_test = mzml_file(path="./Datasets/MT02Dataset/batch11_MT_20210805_005.mzML")
        mode = get_ionization_mode_mzml(to_test)
        self.assertEqual(mode, "pos")

    @staticmethod
    def dummy_command(x):
        return x * 2

    def test_bulk_process(self):
        arguments = [1, 2, 3, 4, 5]
        results = bulk_process(self.dummy_command, arguments)
        self.assertEqual(results, [2, 4, 6, 8, 10])

if __name__ == '__main__':
    unittest.main()