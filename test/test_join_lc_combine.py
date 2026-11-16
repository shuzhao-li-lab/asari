"""Tests for LC annotation_combine helpers in asari.join."""
import unittest
from asari.join import map_ms2_to_features, compound_keys, ms2_entry_fields


class TestLCAnnotationCombine(unittest.TestCase):
    features = [
        {'id': 'F1', 'mz': 116.0707, 'rtime': 147.0, 'rtime_left_base': 140.0, 'rtime_right_base': 155.0, 'parent_masstrack_id': 5},
        {'id': 'F2', 'mz': 116.0707, 'rtime': 60.0, 'rtime_left_base': 55.0, 'rtime_right_base': 65.0, 'parent_masstrack_id': 5},
        {'id': 'F3', 'mz': 200.1000, 'rtime': 30.0, 'rtime_left_base': 28.0, 'rtime_right_base': 33.0, 'parent_masstrack_id': 9},
    ]

    def test_map_by_track_and_rtime(self):
        spectra = [
            {'id': 'a', 'precursor_mz': 116.0707, 'rtime': 150.0, 'parent_masstrack_id': 5},   # in F1 range
            {'id': 'b', 'precursor_mz': 116.0707, 'rtime': 70.0, 'parent_masstrack_id': 5},    # near F2, within 10 s
            {'id': 'c', 'precursor_mz': 116.0707, 'rtime': 100.0, 'parent_masstrack_id': 5},   # no feature
            {'id': 'd', 'precursor_mz': 200.1004, 'rtime': 31.0},                                # by m/z, 2 ppm
            {'id': 'e', 'precursor_mz': 200.1100, 'rtime': 31.0},                                # m/z too far
        ]
        m = map_ms2_to_features(spectra, self.features, mz_tolerance_ppm=5, rt_tolerance=10)
        self.assertEqual(m, {'a': 'F1', 'b': 'F2', 'd': 'F3'})

    def test_compound_keys(self):
        a = compound_keys('L-Proline', 'ONIBWKKTOPOVIA-BYPYZUCNSA-N')
        b = compound_keys('PROLINE', 'ONIBWKKTOPOVIA-UHFFFAOYSA-N')
        c = compound_keys('Proline', '')
        self.assertTrue(a & b)            # same skeleton InChIKey
        self.assertTrue(b & c)            # same name
        self.assertFalse(a & compound_keys('Choline', 'OEYIOHPDSNJKLS-UHFFFAOYSA-N'))
        self.assertEqual(compound_keys('', 'None'), set())

    def test_ms2_entry_fields(self):
        f = ms2_entry_fields({'NAME': 'Choline; CE30; OEYIOHPDSNJKLS-UHFFFAOYSA-N', 'PRECURSORTYPE': '[M]+',
                              'INCHIKEY': 'OEYIOHPDSNJKLS-UHFFFAOYSA-N', 'precursor_mz': 104.107})
        self.assertEqual((f['name'], f['ion'], f['lib_mz']), ('Choline', '[M]+', 104.107))


if __name__ == '__main__':
    unittest.main()
