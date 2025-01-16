import unittest
import pandas as pd
from asari.qc import asari_qc_plot, get_dataframe_from_file
import os

class TestAsariQC(unittest.TestCase):

    def setUp(self):
        # Create a sample dataframe for testing
        self.data = pd.DataFrame({
            'snr': [10, 20, 30],
            'peak_area': [100, 200, 300],
            'goodness_fitting': [0.9, 0.8, 0.85],
            'cSelectivity': [0.1, 0.2, 0.3]
        })
        self.outfile = "test_qc_plot.pdf"

    def tearDown(self):
        # Remove the generated file after test
        if os.path.exists(self.outfile):
            os.remove(self.outfile)

    def test_asari_qc_plot(self):
        # Test the asari_qc_plot function
        asari_qc_plot(self.data, outfile=self.outfile)
        self.assertTrue(os.path.exists(self.outfile))

    def test_get_dataframe_from_file(self):
        # Create a temporary file for testing
        test_file = "test_feature_table.tsv"
        self.data.to_csv(test_file, sep='\t')

        # Test the get_dataframe_from_file function
        df = get_dataframe_from_file(test_file)
        pd.testing.assert_frame_equal(df, self.data)

        # Clean up
        os.remove(test_file)

if __name__ == '__main__':
    unittest.main()