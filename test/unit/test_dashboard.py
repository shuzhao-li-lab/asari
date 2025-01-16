import unittest
import os
import pickle
import json
import pandas as pd
from unittest.mock import patch, mock_open
from asari.dashboard import read_project, epd_convert, track_to_peaks, find_track_by_mz, find_a_good_peak, prepare_rt_alignment

class TestDashboard(unittest.TestCase):

    @patch("builtins.open", new_callable=mock_open, read_data='{"project_name": "test_project", "number_of_samples": 10}')
    @patch("os.path.join", return_value="dummy_path")
    @patch("pickle.load", return_value={"dummy_key": "dummy_value"})
    @patch("pandas.read_csv", return_value=pd.DataFrame({"mz": [100, 200], "peak_area": [1000, 2000], "rtime": [5, 10], "snr": [10, 20], "goodness_fitting": [0.95, 0.85], "cSelectivity": [0.9, 0.8]}))
    def test_read_project(self, mock_read_csv, mock_pickle_load, mock_path_join, mock_open):
        project_desc, cmap, epd, Ftable = read_project("dummy_dir")
        self.assertEqual(project_desc["project_name"], "test_project")
        self.assertEqual(cmap["dummy_key"], "dummy_value")
        self.assertEqual(epd["dummy_key"], "dummy_value")
        self.assertTrue(isinstance(Ftable, pd.DataFrame))

    def test_epd_convert(self):
        epd_dict = {
            "compound1": {
                "MS1_pseudo_Spectra": [{"id_number": "P1"}, {"id_number": "P2"}]
            }
        }
        peakDict, epdDict = epd_convert(epd_dict)
        self.assertEqual(peakDict["P1"]["id_number"], "P1")
        self.assertEqual(epdDict["compound1"]["MS1_pseudo_Spectra"], ["P1", "P2"])

    def test_track_to_peaks(self):
        peakDict = {
            "P1": {"parent_masstrack_id": "T1", "id_number": "P1"},
            "P2": {"parent_masstrack_id": "T1", "id_number": "P2"},
            "P3": {"parent_masstrack_id": "T2", "id_number": "P3"}
        }
        track2peaks = track_to_peaks(peakDict)
        self.assertEqual(track2peaks["T1"], ["P1", "P2"])
        self.assertEqual(track2peaks["T2"], ["P3"])

    def test_find_track_by_mz(self):
        cmap = {
            "list_mass_tracks": {
                "T1": {"mz": 100, "id_number": "T1"},
                "T2": {"mz": 200, "id_number": "T2"}
            }
        }
        rt_list = [0, 1, 2]
        track_id = find_track_by_mz(cmap, rt_list, 150)
        self.assertEqual(track_id, "T1")

    def test_find_a_good_peak(self):
        peakDict = {
            "P1": {"goodness_fitting": 0.95, "cSelectivity": 0.95},
            "P2": {"goodness_fitting": 0.85, "cSelectivity": 0.85}
        }
        good_peak = find_a_good_peak(peakDict)
        self.assertEqual(good_peak["goodness_fitting"], 0.95)
        self.assertEqual(good_peak["cSelectivity"], 0.95)

    def test_prepare_rt_alignment(self):
        cmap = {
            "rt_records": [
                {"reverse_rt_cal_dict": {0: 0, 1: 1, 2: 2}, "name": "sample1"},
                {"reverse_rt_cal_dict": {0: 0, 1: 1, 2: 2}, "name": "sample2"}
            ],
            "rt_length": 3,
            "dict_scan_rtime": {0: 0, 1: 1, 2: 2}
        }
        df = prepare_rt_alignment(cmap)
        self.assertTrue(isinstance(df, pd.DataFrame))
        self.assertEqual(df.shape, (3, 2))

if __name__ == "__main__":
    unittest.main()