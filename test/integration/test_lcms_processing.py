import os
import requests
import zipfile
import io
import pandas as pd
import subprocess
import unittest
import json

os.makedirs("./Datasets", exist_ok=True)

class TestLCMSProcessing(unittest.TestCase):
    def __init__(self, methodName = "runTest"):
        super().__init__(methodName)

    def setUpClass():
        instance = TestLCMSProcessing()
        #instance.download_and_extract_datasets()
        #instance.run_asari()
        #instance.test_outputs()

    @staticmethod
    def download_and_extract_datasets():
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
        
    def run_asari(self):
        command = [
            "asari", "process",
            "-i", "./Datasets/MT02Dataset/",
            "-o", os.path.dirname(__file__) + "/asari/",
            "-m", "pos",
            "-j", "MT02_results",
            "--keep_intermediates", "true",
            "--database_mode", "ondisk"
        ]
        results = subprocess.run(command, check=True)
        self.assertEqual(results.returncode, 0)

    def test_outputs(self):
        asari_subdir = None
        for x in os.listdir(os.path.dirname(__file__)):
            if x.startswith("asari_MT02_results"):
                asari_subdir = os.path.join(os.path.dirname(__file__), x)


        preferred_table = asari_subdir + "/preferred_Feature_table.tsv"
        self.assertTrue(os.path.exists(preferred_table), "Preferred feature table not found")
        
        ft = pd.read_csv(preferred_table, sep="\t")

        self.assertIn('cSelectivity', ft.columns, "cSelectivity column missing in feature table")
        self.assertIn('goodness_fitting', ft.columns, "goodness_fitting column missing in feature table")
        self.assertIn('snr', ft.columns, "snr column missing in feature table")
        self.assertIn('peak_area', ft.columns, "peak_area column missing in feature table")
        self.assertIn('rtime', ft.columns, "rtime column missing in feature table")
        self.assertIn('mz', ft.columns, "mz column missing in feature table")
        self.assertIn('id_number', ft.columns, "feature_id column missing in feature table")

        self.assertGreater(ft.shape[1] - 11, 0, "No samples found in feature table")
        self.assertGreater(ft.shape[0] - 1, 0, "No features found in feature table")
        self.assertEqual(ft.shape[1], 19)
        self.assertEqual(ft.shape[0], 5384)

        annotation_table = asari_subdir + "/Feature_annotation.tsv"
        self.assertTrue(os.path.exists(annotation_table), "Annotation table not found")
        
        annots = pd.read_csv(annotation_table, sep="\t")
        self.assertGreater(annots.shape[0], 0, "No annotations found in annotation table")

        empCpds = asari_subdir + "/Annotated_empiricalCompounds.json"
        self.assertTrue(os.path.exists(empCpds), "Annotated empirical compounds not found")
        with open(empCpds) as f:
            self.assertIsInstance(json.load(f), dict, "Annotated empirical compounds not found")

    def run_asari_memory(self):
        command = [
            "asari", "process",
            "-i", "./Datasets/MT02Dataset/",
            "-o", os.path.dirname(__file__) + "/asari/",
            "-m", "pos",
            "-j", "MT02_results",
            "--keep_intermediates", "true",
            "--database_mode", "memory"
        ]
        results = subprocess.run(command, check=True)
        self.assertEqual(results.returncode, 0)

    def test_outputs_memory(self):
        asari_subdir = None
        for x in os.listdir(os.path.dirname(__file__)):
            if x.startswith("asari_MT02_results"):
                asari_subdir = os.path.join(os.path.dirname(__file__), x)


        preferred_table = asari_subdir + "/preferred_Feature_table.tsv"
        self.assertTrue(os.path.exists(preferred_table), "Preferred feature table not found")
        
        ft = pd.read_csv(preferred_table, sep="\t")

        self.assertIn('cSelectivity', ft.columns, "cSelectivity column missing in feature table")
        self.assertIn('goodness_fitting', ft.columns, "goodness_fitting column missing in feature table")
        self.assertIn('snr', ft.columns, "snr column missing in feature table")
        self.assertIn('peak_area', ft.columns, "peak_area column missing in feature table")
        self.assertIn('rtime', ft.columns, "rtime column missing in feature table")
        self.assertIn('mz', ft.columns, "mz column missing in feature table")
        self.assertIn('id_number', ft.columns, "feature_id column missing in feature table")

        self.assertGreater(ft.shape[1] - 11, 0, "No samples found in feature table")
        self.assertGreater(ft.shape[0] - 1, 0, "No features found in feature table")
        self.assertEqual(ft.shape[1], 19)
        self.assertEqual(ft.shape[0], 5384)

        annotation_table = asari_subdir + "/Feature_annotation.tsv"
        self.assertTrue(os.path.exists(annotation_table), "Annotation table not found")
        
        annots = pd.read_csv(annotation_table, sep="\t")
        self.assertGreater(annots.shape[0], 0, "No annotations found in annotation table")

        empCpds = asari_subdir + "/Annotated_empiricalCompounds.json"
        self.assertTrue(os.path.exists(empCpds), "Annotated empirical compounds not found")
        with open(empCpds) as f:
            self.assertIsInstance(json.load(f), dict, "Annotated empirical compounds not found")

    def run_asari_large_project(self):
        command = [
            "asari", "process",
            "-i", "./Datasets/MT02Dataset/",
            "-o", os.path.dirname(__file__) + "/asari/",
            "-m", "pos",
            "-j", "MT02_results",
            "--keep_intermediates", "true",
            "--database_mode", "memory"
            "--project_sample_number_small", 0
        ]
        results = subprocess.run(command, check=True)
        self.assertEqual(results.returncode, 0)

    def test_outputs_large_project(self):
        asari_subdir = None
        for x in os.listdir(os.path.dirname(__file__)):
            if x.startswith("asari_MT02_results"):
                asari_subdir = os.path.join(os.path.dirname(__file__), x)


        preferred_table = asari_subdir + "/preferred_Feature_table.tsv"
        self.assertTrue(os.path.exists(preferred_table), "Preferred feature table not found")
        
        ft = pd.read_csv(preferred_table, sep="\t")

        self.assertIn('cSelectivity', ft.columns, "cSelectivity column missing in feature table")
        self.assertIn('goodness_fitting', ft.columns, "goodness_fitting column missing in feature table")
        self.assertIn('snr', ft.columns, "snr column missing in feature table")
        self.assertIn('peak_area', ft.columns, "peak_area column missing in feature table")
        self.assertIn('rtime', ft.columns, "rtime column missing in feature table")
        self.assertIn('mz', ft.columns, "mz column missing in feature table")
        self.assertIn('id_number', ft.columns, "feature_id column missing in feature table")

        self.assertGreater(ft.shape[1] - 11, 0, "No samples found in feature table")
        self.assertGreater(ft.shape[0] - 1, 0, "No features found in feature table")
        self.assertEqual(ft.shape[1], 19)
        self.assertEqual(ft.shape[0], 5307)

        annotation_table = asari_subdir + "/Feature_annotation.tsv"
        self.assertTrue(os.path.exists(annotation_table), "Annotation table not found")
        
        annots = pd.read_csv(annotation_table, sep="\t")
        self.assertGreater(annots.shape[0], 0, "No annotations found in annotation table")

        empCpds = asari_subdir + "/Annotated_empiricalCompounds.json"
        self.assertTrue(os.path.exists(empCpds), "Annotated empirical compounds not found")
        with open(empCpds) as f:
            self.assertIsInstance(json.load(f), dict, "Annotated empirical compounds not found")

if __name__ == "__main__":
    unittest.main()