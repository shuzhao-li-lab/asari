"""
Tests for asari.tools.msp_parser: MSP and MGF parsing, field-name standardization and RI values.
"""
import os
import tempfile
import unittest

from asari.tools.msp_parser import (parse_msp_to_listdict, parse_mgf, parse_ri_value,
                                    standardize_entry, _parse_mgf_charge)


def write_temp(text, suffix, newline=None):
    """Write text to a temporary file; return its path."""
    fd, path = tempfile.mkstemp(suffix=suffix)
    with os.fdopen(fd, 'w', newline=newline) as f:
        f.write(text)
    return path


MGF_TEXT = """# file comment
BEGIN IONS
TITLE=spec1 scan=12
PEPMASS=301.1412 15000
CHARGE=2-
RTINSECONDS=125.5
NAME=Compound A
RETENTION_INDEX=1105.2
101.05 450
120.07,30
140.1\t10\t"y1"
END IONS

begin ions
PEPMASS=abc
CHARGE=1+
50 1
END IONS

BEGIN IONS
PEPMASS=200.1
CHARGE=-1
; semicolon comment
60 5
BEGIN IONS
PEPMASS=400.2
RTINMINUTES=3.5
70 6
END IONS
BEGIN IONS
PEPMASS=500.5
80 8
"""

MSP_TEXT = ('NAME:Alanine 2TMS\r\n'
            'retention index: SemiStdNP=1105.2\r\n'
            'MW: 233\r\n'
            'Comments: "note: has colon" x\r\n'
            'Synon: first\r\n'
            'Synon: second\r\n'
            'Num.Peaks: 4\r\n'
            '73 999; 116 500;\r\n'
            '147\t200 "frag: a"\r\n'
            '(190 50)\r\n'
            '  \r\n'
            '# comment line\r\n'
            'Compound Name: B\r\n'
            'Kovats_RI: 1500\r\n'
            'NUM PEAKS: 2\r\n'
            '50,10 60,20\r\n'
            'Name: C\r\n'
            'RT: 3.5\r\n'
            'Num Peaks: 1\r\n'
            '100 1\r\n')


class TestParseMGF(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.path = write_temp(MGF_TEXT, '.mgf')
        cls.spectra = parse_mgf(cls.path, warn=False)

    @classmethod
    def tearDownClass(cls):
        os.remove(cls.path)

    def test_count_includes_unterminated(self):
        # 3rd spectrum has no END IONS before next BEGIN IONS; last one ends at EOF
        self.assertEqual(len(self.spectra), 5)

    def test_peaks_formats(self):
        # whitespace, comma and tab with annotation
        self.assertEqual(self.spectra[0]['peaks'], [(101.05, 450.0), (120.07, 30.0), (140.1, 10.0)])
        self.assertEqual(self.spectra[0]['mzs'], [101.05, 120.07, 140.1])
        self.assertEqual(self.spectra[0]['intensities'], [450.0, 30.0, 10.0])

    def test_precursor(self):
        self.assertEqual(self.spectra[0]['precursor_mz'], 301.1412)
        self.assertEqual(self.spectra[0]['precursor_intensity'], 15000.0)
        self.assertEqual(self.spectra[0]['params']['PEPMASS'], (301.1412, 15000.0))

    def test_nonnumeric_pepmass(self):
        self.assertIsNone(self.spectra[1]['precursor_mz'])
        self.assertEqual(self.spectra[1]['peaks'], [(50.0, 1.0)])

    def test_signed_charge(self):
        self.assertEqual([s['charge'] for s in self.spectra], [-2, 1, -1, None, None])

    def test_rtime_seconds(self):
        self.assertEqual(self.spectra[0]['rtime'], 125.5)
        self.assertEqual(self.spectra[3]['rtime'], 210.0)        # RTINMINUTES 3.5
        self.assertIsNone(self.spectra[4]['rtime'])

    def test_comments_skipped(self):
        self.assertEqual(self.spectra[2]['peaks'], [(60.0, 5.0)])

    def test_standardized_fields(self):
        self.assertEqual(self.spectra[0]['Name'], 'Compound A')
        self.assertEqual(self.spectra[0]['RI'], 1105.2)
        self.assertEqual(self.spectra[0]['params']['TITLE'], 'spec1 scan=12')

    def test_no_standardize(self):
        spectra = parse_mgf(self.path, standardize=False, warn=False)
        self.assertNotIn('Name', spectra[0])
        self.assertNotIn('RI', spectra[0])


class TestParseMGFCharge(unittest.TestCase):

    def test_charges(self):
        cases = {'2+': 2, '2-': -2, '-1': -1, '+3': 3, '1': 1,
                 '2+ and 3+': 2, '2+,3+': 2, 'abc': None, '': None}
        for value, expected in cases.items():
            self.assertEqual(_parse_mgf_charge(value), expected, value)


class TestParseMSP(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.path = write_temp(MSP_TEXT, '.msp', newline='')
        cls.entries = parse_msp_to_listdict(cls.path)

    @classmethod
    def tearDownClass(cls):
        os.remove(cls.path)

    def test_count_without_blank_separator(self):
        # records B and C are not separated by a blank line
        self.assertEqual([e['Name'] for e in self.entries], ['Alanine 2TMS', 'B', 'C'])

    def test_peak_formats(self):
        self.assertEqual(self.entries[0]['peaks'],
                         [(73.0, 999.0), (116.0, 500.0), (147.0, 200.0), (190.0, 50.0)])
        self.assertEqual(self.entries[1]['peaks'], [(50.0, 10.0), (60.0, 20.0)])

    def test_ri_variants(self):
        self.assertEqual(self.entries[0]['RI'], 1105.2)          # 'retention index: SemiStdNP=...'
        self.assertEqual(self.entries[1]['RI'], 1500.0)          # 'Kovats_RI'
        self.assertNotIn('RI', self.entries[2])

    def test_retention_time_kept_separate(self):
        self.assertEqual(self.entries[2]['RETENTIONTIME'], '3.5')

    def test_value_with_colon(self):
        self.assertEqual(self.entries[0]['Comments'], '"note: has colon" x')

    def test_repeated_key_keeps_first(self):
        self.assertEqual(self.entries[0]['Synon'], 'first')

    def test_canonical_keys(self):
        self.assertEqual(self.entries[0]['ExactMass'], '233')    # from MW
        self.assertEqual(self.entries[0]['Num Peaks'], '4')      # from Num.Peaks
        self.assertEqual(self.entries[1]['Num Peaks'], '2')      # from NUM PEAKS


class TestRIValue(unittest.TestCase):

    def test_parse_ri_value(self):
        cases = [('1817.1', 1817.1), ('SemiStdNP=3296.3', 3296.3), ('1817.1 (lit)', 1817.1),
                 ('1305/2', 1305.0), (1500, 1500.0), ('', None), (None, None), ('abc', None),
                 (float('nan'), None)]
        for value, expected in cases:
            self.assertEqual(parse_ri_value(value), expected, value)

    def test_standardize_entry(self):
        e = standardize_entry({'Retention_index': '1200.5', 'compound name': 'X', 'RI': None})
        self.assertIsNone(e['RI'])           # existing canonical key is not overwritten
        self.assertEqual(e['Name'], 'X')
        e = standardize_entry({'Retention_index': '1200.5'})
        self.assertEqual(e['RI'], 1200.5)
        self.assertEqual(e['Retention_index'], '1200.5')      # original key kept


if __name__ == '__main__':
    unittest.main()
