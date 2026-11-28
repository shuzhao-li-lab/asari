"""asari deconvolute passes the khipu workflow (regular or tracing) to khipu (khipu >= 3)."""
import json
import os
import tempfile
import unittest
from unittest.mock import patch

from asari.default_parameters import PARAMETERS
from asari.deconvolute import khipu_deconvolute_featuretable
from asari.main import build_parser, update_params_from_CLI

HERE = os.path.dirname(os.path.abspath(__file__))
FEATURE_TABLE = os.path.join(HERE, 'HighOnly_HILICpos_preferred_Feature_table.tsv')


def _count_13C(isotope):
    name = isotope.split(',')[0].strip()
    if name == 'M0':
        return 0
    if name == '13C/12C':
        return 1
    if name.startswith('13C/12C*'):
        return int(name.split('*')[1])
    return None


def _deconvolute(**params):
    parameters = dict(PARAMETERS, mode='pos', **params)
    with tempfile.TemporaryDirectory() as tmp, patch('sys.stdout'):
        out = os.path.join(tmp, 'khipus.json')
        khipu_deconvolute_featuretable(FEATURE_TABLE, parameters, out)
        with open(out) as f:
            return json.load(f)


class TestKhipuWorkflow(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.regular = _deconvolute(khipu_workflow='regular')
        cls.tracing = _deconvolute(khipu_workflow='tracing')

    def test_regular_isotopologues_decrease(self):
        """khipu regular rule: within each adduct, M0 > M1 > M2 ... (intensity from peak_area)."""
        checked = 0
        for E in self.regular:
            series = {}
            for p in E['MS1_pseudo_Spectra']:
                n = _count_13C(p.get('isotope', ''))
                mod = p.get('modification', '')
                if n is not None and all('charged' in x for x in mod.split(',')[1:]):
                    series.setdefault(mod, {})[n] = p['peak_area']
            for by_n in series.values():
                for n in by_n:
                    if n > 0:
                        self.assertIn(n - 1, by_n)
                        self.assertLess(by_n[n], by_n[n - 1])
                        checked += 1
        self.assertGreater(checked, 10)

    def test_tracing_searches_more_13C(self):
        max_n = lambda empCpds: max(_count_13C(p['isotope']) or 0 for E in empCpds
                                    for p in E['MS1_pseudo_Spectra'] if 'isotope' in p)
        self.assertLessEqual(max_n(self.regular), 3)
        self.assertNotEqual(self.regular, self.tracing)

    def test_ratio_check_ignored_for_tracing(self):
        self.assertEqual(_deconvolute(khipu_workflow='tracing', check_isotope_ratio=True), self.tracing)
        self.assertNotEqual(_deconvolute(khipu_workflow='regular', check_isotope_ratio=True), self.regular)

    def test_unknown_workflow(self):
        with self.assertRaises(ValueError):
            _deconvolute(khipu_workflow='bogus')


class TestKhipuWorkflowCLI(unittest.TestCase):

    def _params(self, *argv):
        with patch('sys.argv', ['asari', 'deconvolute', '-i', FEATURE_TABLE, *argv]):
            args = build_parser()
        parameters = dict(PARAMETERS)
        update_params_from_CLI(parameters, args)
        return parameters

    def test_default_regular(self):
        self.assertEqual(self._params()['khipu_workflow'], 'regular')

    def test_flag_and_synonyms(self):
        self.assertEqual(self._params('--khipu_workflow', 'tracing')['khipu_workflow'], 'tracing')
        self.assertEqual(self._params('--khipu_workflow', 'Labeled')['khipu_workflow'], 'tracing')

    def test_unknown_flag_value(self):
        with self.assertRaises(AssertionError):
            self._params('--khipu_workflow', 'bogus')


if __name__ == '__main__':
    unittest.main()
