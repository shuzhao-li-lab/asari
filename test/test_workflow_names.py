"""Tests for workflow name normalization and synonyms (asari.main.normalize_workflow)."""
import unittest
from asari.main import normalize_workflow


class TestWorkflowNames(unittest.TestCase):
    def test_synonyms_and_case(self):
        for name, expected in [('LC', 'LC'), ('LCMS', 'LC'), ('lc-ms', 'LC'), ('GCMS', 'GC'), ('GC-MS', 'GC'),
                               ('LCMSMS', 'LCMSMS'), ('LC-MS/MS', 'LCMSMS'), ('ms/ms', 'MSMS'), ('dims', 'DIMS')]:
            self.assertEqual(normalize_workflow(name), expected)

    def test_unknown(self):
        self.assertIsNone(normalize_workflow('Bogus'))


if __name__ == '__main__':
    unittest.main()
