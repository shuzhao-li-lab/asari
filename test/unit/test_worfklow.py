import unittest
import os
import shutil
import tempfile
from unittest.mock import patch, MagicMock

from asari.workflow import (
    process_project,
    read_project_dir,
    register_samples,
    create_export_folders,
    remove_intermediate_pickles,
    batch_EIC_from_samples_,
    single_sample_EICs_,
    process_xics,
    get_mz_list
)

class TestWorkflow(unittest.TestCase):

    def setUp(self):
        self.test_dir = tempfile.mkdtemp()
        self.parameters = {
            'database_mode': 'auto',
            'project_sample_number_small': 5,
            'outdir': self.test_dir,
            'project_name': 'test_project',
            'mz_tolerance_ppm': 10,
            'min_intensity_threshold': 1000,
            'min_timepoints': 5,
            'min_peak_height': 500,
            'intensity_multiplier': 1.0,
            'multicores': 2,
            'mode': 'positive',
            'pickle': False,
            'anno': True,
            'tmp_pickle_dir': os.path.join(self.test_dir, 'tmp_pickle_dir')
        }
        self.list_input_files = [os.path.join(self.test_dir, f'test_file_{i}.mzML') for i in range(3)]
        for file in self.list_input_files:
            with open(file, 'w') as f:
                f.write("dummy content")

    def tearDown(self):
        shutil.rmtree(self.test_dir)

    def test_register_samples(self):
        sample_registry = register_samples(self.list_input_files)
        self.assertEqual(len(sample_registry), 3)
        for i, file in enumerate(self.list_input_files):
            self.assertEqual(sample_registry[i]['input_file'], file)

    def test_create_export_folders(self):
        time_stamp = '123456'
        create_export_folders(self.parameters, time_stamp)
        self.assertTrue(os.path.exists(self.parameters['outdir']))
        self.assertTrue(os.path.exists(self.parameters['tmp_pickle_dir']))
        self.assertTrue(os.path.exists(os.path.join(self.parameters['outdir'], 'export')))

    def test_remove_intermediate_pickles(self):
        os.makedirs(self.parameters['tmp_pickle_dir'])
        with open(os.path.join(self.parameters['tmp_pickle_dir'], 'test.pickle'), 'w') as f:
            f.write("dummy content")
        remove_intermediate_pickles(self.parameters)
        self.assertFalse(os.path.exists(self.parameters['tmp_pickle_dir']))

    @patch('asari.workflow.pymzml.run.Reader')
    @patch('asari.workflow.extract_massTracks_')
    @patch('asari.workflow.find_mzdiff_pairs_from_masstracks')
    def test_single_sample_EICs_(self, mock_find_mzdiff_pairs, mock_extract_massTracks, mock_pymzml_reader):
        mock_pymzml_reader.return_value = MagicMock()
        mock_extract_massTracks.return_value = {
            'rt_numbers': [1, 2, 3],
            'tracks': [(100, [10, 20, 30]), (200, [15, 25, 35])]
        }
        mock_find_mzdiff_pairs.return_value = [(100, 200)]

        shared_dict = {}
        single_sample_EICs_(0, self.list_input_files[0], 'positive', 'memory', 10, 1000, 5, 500, 1.5, 'outfile', shared_dict)
        self.assertIn(0, shared_dict)
        self.assertEqual(shared_dict[0][0], 'failed')

    @patch('asari.workflow.batch_EIC_from_samples_')
    @patch('asari.workflow.ext_Experiment')
    def test_process_project(self, mock_ext_Experiment, mock_batch_EIC_from_samples):
        mock_batch_EIC_from_samples.return_value = {
            0: ('passed', 'passed', 'outfile', 3, [1, 2, 3], [0.1, 0.2, 0.3], [(100, 0)], 1, [(100, 200)], {})
        }
        mock_ext_Experiment.return_value = MagicMock()
        response = process_project(self.list_input_files, self.parameters)
        mock_ext_Experiment.assert_called_once()

    def test_read_project_dir(self):
        files = read_project_dir(self.test_dir, file_pattern='.mzML')
        self.assertEqual(len(files), 3)
        for file in files:
            self.assertTrue(file.endswith('.mzML'))

    def test_get_mz_list(self):
        mz_list_file = os.path.join(self.test_dir, 'mz_list.txt')
        with open(mz_list_file, 'w') as f:
            f.write("m/z\n100.0\n200.0\n300.0\n")
        mz_list = get_mz_list(mz_list_file)
        self.assertEqual(mz_list, [100.0, 200.0, 300.0])

if __name__ == '__main__':
    unittest.main()