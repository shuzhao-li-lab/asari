"""Since khipu 3, a feature found by khipu's extended search can be in several empCpds.
Each feature belongs to one empCpd: where it is a core member, else the first that has it."""
import unittest

from asari.annotate import combine_preannotation
from asari.join import empCpd_targets
from asari.utils import is_extended_modification


def _f(fid, isotope, modification):
    return {'id': fid, 'isotope': isotope, 'modification': modification,
            'ion_relation': f'{isotope},{modification}'}


LIB = {'glc': {'neutral_mass': 180.06339, 'ion': 'M+H+', 'isotope': 'M0'},
       'b': {'neutral_mass': 200.0, 'ion': 'M+H+', 'isotope': 'M0'},
       'x': {'neutral_mass': 145.05, 'ion': 'M+H+', 'isotope': 'M0'}}


def _empcpds(x_in_b='M+H+,NH3'):
    A = {'interim_id': 'kpA', 'neutral_formula_mass': 180.06339, 'MS1_pseudo_Spectra': [
        _f('F1', 'M0', 'M+H+'), _f('F2', '13C/12C', 'M+H+'), _f('X', 'M0', 'M+H+,H2O')]}
    B = {'interim_id': 'kpB', 'neutral_formula_mass': 200.0, 'MS1_pseudo_Spectra': [
        _f('G1', 'M0', 'M+H+'), _f('X', 'M0', x_in_b)]}
    return [A, B]


class TestSharedFeatures(unittest.TestCase):

    def test_is_extended_modification(self):
        for m, expected in [('M+H+', False), ('Na/H, 2x charged', False), ('M+H+,H2O', True),
                            ('M+H+, 2x charged,H2O, 2x charged', True)]:
            self.assertEqual(is_extended_modification(m), expected, m)

    def test_extended_feature_belongs_to_first_empcpd(self):
        combined, feature_empCpd = combine_preannotation(
            {'F1': ['glc'], 'G1': ['b']}, _empcpds(), LIB, mode='pos', mz_tolerance_ppm=5)
        self.assertEqual(feature_empCpd['X']['interim_id'], 'kpA')
        self.assertEqual(feature_empCpd['X']['modification'], 'M+H+,H2O')
        self.assertEqual(combined['X'], ['glc'])            # from kpA only, not kpB
        self.assertEqual(combined['F2'], ['glc'])
        self.assertEqual(combined['G1'], ['b'])

    def test_core_membership_wins_over_order(self):
        # X is a core member of kpB (listed second) and extended in kpA
        combined, feature_empCpd = combine_preannotation(
            {'F1': ['glc'], 'G1': ['b']}, _empcpds(x_in_b='Na/H'), LIB, mode='pos', mz_tolerance_ppm=5)
        self.assertEqual(feature_empCpd['X']['interim_id'], 'kpB')
        self.assertEqual(combined['X'], ['b'])

    def test_extended_feature_keeps_direct_match(self):
        combined, _ = combine_preannotation(
            {'F1': ['glc'], 'G1': ['b'], 'X': ['x']}, _empcpds(), LIB, mode='pos', mz_tolerance_ppm=5)
        self.assertEqual(combined['X'], ['x'])
        self.assertEqual(combined['F1'], ['glc'])

    def test_empcpd_targets(self):
        _, feature_empCpd = combine_preannotation({}, _empcpds(), LIB, mode='pos', mz_tolerance_ppm=5)
        members = {E['interim_id']: [F['id'] for F in E['MS1_pseudo_Spectra']] for E in _empcpds()}
        self.assertEqual(empCpd_targets('F1', 'kpA', members, feature_empCpd), ['F1', 'F2', 'X'])
        self.assertEqual(empCpd_targets('G1', 'kpB', members, feature_empCpd), ['G1'])   # X is kpA's
        self.assertEqual(empCpd_targets('S1', '', members, feature_empCpd), ['S1'])


if __name__ == '__main__':
    unittest.main()
