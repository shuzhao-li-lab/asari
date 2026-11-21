import unittest

from asari.deconvolute import khipu_feature_adaptor
from asari.default_parameters import isotope_search_patterns, adduct_search_patterns_pos, extended_adducts
from khipu.extended import peaklist_to_khipu_list


def _features(m1_area):
    # M0 (M+H+ of a 180.0634 compound) and its 13C M+1, same RT
    return [
        {'id': 'F1', 'id_number': 'F1', 'mz': 181.0707, 'rtime': 100.0, 'parent_masstrack_id': 1, 'peak_area': 1e6},
        {'id': 'F2', 'id_number': 'F2', 'mz': 182.0741, 'rtime': 100.0, 'parent_masstrack_id': 2, 'peak_area': m1_area},
    ]


def _grouped(features, check):
    khipus, _ = peaklist_to_khipu_list(
        khipu_feature_adaptor(features), isotope_search_patterns=isotope_search_patterns,
        adduct_search_patterns=adduct_search_patterns_pos, extended_adducts=extended_adducts,
        mz_tolerance_ppm=5, rt_tolerance=2, check_isotope_ratio=check, mode='pos',
        charges=[1], has_parent_masstrack=True)
    return any({'F1', 'F2'} <= set(k.nodes_to_use) for k in khipus)


class TestCheckIsotopeRatio(unittest.TestCase):
    def test_adaptor(self):
        f = khipu_feature_adaptor([{'peak_area': 5.0}, {'peak_area': 5.0, 'height': 3.0}])
        self.assertEqual((f[0]['representative_intensity'], f[0]['height']), (5.0, 5.0))
        self.assertEqual(f[1]['height'], 3.0)          # existing height kept

    def test_ratio_within_limit(self):
        self.assertTrue(_grouped(_features(1e5), check=True))      # M+1/M0 = 0.1

    def test_ratio_out_of_limit(self):
        self.assertFalse(_grouped(_features(2e6), check=True))     # M+1/M0 = 2
        self.assertTrue(_grouped(_features(2e6), check=False))


if __name__ == '__main__':
    unittest.main()
