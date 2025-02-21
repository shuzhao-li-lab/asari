import unittest

from asari.default_parameters import PARAMETERS

class TestDefaultParameters(unittest.TestCase):

    def test_project_name(self):
        self.assertEqual(PARAMETERS['project_name'], 'asari_project')

    def test_outdir(self):
        self.assertEqual(PARAMETERS['outdir'], 'output')

    def test_pickle(self):
        self.assertFalse(PARAMETERS['pickle'])

    def test_database_mode(self):
        self.assertEqual(PARAMETERS['database_mode'], 'auto')

    def test_multicores(self):
        self.assertEqual(PARAMETERS['multicores'], 4)

    def test_mode(self):
        self.assertEqual(PARAMETERS['mode'], 'pos')

    def test_mass_range(self):
        self.assertEqual(PARAMETERS['mass_range'], (50, 2000))

    def test_mz_tolerance_ppm(self):
        self.assertEqual(PARAMETERS['mz_tolerance_ppm'], 5)

    def test_min_timepoints(self):
        self.assertEqual(PARAMETERS['min_timepoints'], 6)

    def test_signal_noise_ratio(self):
        self.assertEqual(PARAMETERS['signal_noise_ratio'], 2)

    def test_min_intensity_threshold(self):
        self.assertEqual(PARAMETERS['min_intensity_threshold'], 1000)

    def test_min_peak_height(self):
        self.assertEqual(PARAMETERS['min_peak_height'], 100000)

    def test_min_peak_ratio(self):
        self.assertEqual(PARAMETERS['min_peak_ratio'], 0.001)

    def test_wlen(self):
        self.assertEqual(PARAMETERS['wlen'], 25)

    def test_autoheight(self):
        self.assertFalse(PARAMETERS['autoheight'])

    def test_gaussian_shape(self):
        self.assertEqual(PARAMETERS['gaussian_shape'], 0.5)

    def test_peak_area(self):
        self.assertEqual(PARAMETERS['peak_area'], 'sum')

    def test_reference(self):
        self.assertIsNone(PARAMETERS['reference'])

    def test_rt_align_method(self):
        self.assertEqual(PARAMETERS['rt_align_method'], 'lowess')

    def test_rt_align_on(self):
        self.assertTrue(PARAMETERS['rt_align_on'])

    def test_drop_unaligned_samples(self):
        self.assertFalse(PARAMETERS['drop_unaligned_samples'])

    def test_rtime_tolerance(self):
        self.assertEqual(PARAMETERS['rtime_tolerance'], 50)

    def test_cal_min_peak_height(self):
        self.assertEqual(PARAMETERS['cal_min_peak_height'], 100000)

    def test_peak_number_rt_calibration(self):
        self.assertEqual(PARAMETERS['peak_number_rt_calibration'], 20)

    def test_max_retention_shift(self):
        self.assertIsNone(PARAMETERS['max_retention_shift'])

    def test_num_lowess_iterations(self):
        self.assertEqual(PARAMETERS['num_lowess_iterations'], 3)

    def test_project_sample_number_small(self):
        self.assertEqual(PARAMETERS['project_sample_number_small'], 10)

    def test_output_feature_table(self):
        self.assertEqual(PARAMETERS['output_feature_table'], 'Feature_table.tsv')

    def test_mass_grid_mapping(self):
        self.assertEqual(PARAMETERS['mass_grid_mapping'], "_mass_grid_mapping.csv")

    def test_annotation_filename(self):
        self.assertEqual(PARAMETERS['annotation_filename'], "Annotation_table.tsv")

    def test_json_empricalCompounds(self):
        self.assertEqual(PARAMETERS['json_empricalCompounds'], "_empCpd_json.json")

    def test_check_isotope_ratio(self):
        self.assertFalse(PARAMETERS['check_isotope_ratio'])

if __name__ == '__main__':
    unittest.main()