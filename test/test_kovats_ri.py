import os
import tempfile
import unittest
import numpy as np

from asari.gc_functions import read_fit_KovatsIndex_rtime, rtime_to_ri, report_ri_extrapolation


class TestKovatsRI(unittest.TestCase):
    def _model(self, rows):
        with tempfile.TemporaryDirectory() as d:
            f = os.path.join(d, 'kovats.tsv')
            with open(f, 'w') as O:
                O.write("Kovats_index\tMeasured_time(minutes)\n")
                for k, t in rows:
                    O.write(f"{k}\t{t}\n")
            return read_fit_KovatsIndex_rtime(f)

    def test_linear_between_alkanes(self):
        m = self._model([(1000, 2.0), (1100, 3.0), (1200, 5.0)])
        self.assertAlmostEqual(rtime_to_ri(2.0 * 60, m), 1000)
        self.assertAlmostEqual(rtime_to_ri(2.5 * 60, m), 1050)
        self.assertAlmostEqual(rtime_to_ri(4.0 * 60, m), 1150)
        self.assertAlmostEqual(rtime_to_ri(5.0 * 60, m), 1200)

    def test_extrapolation_capped(self):
        # end slopes: 100 RI/min below 1000, 50 RI/min above 1200
        m = self._model([(1000, 2.0), (1100, 3.0), (1200, 5.0)])
        self.assertAlmostEqual(rtime_to_ri(1.5 * 60, m), 950)
        self.assertAlmostEqual(rtime_to_ri(1.0 * 60, m), 900)
        self.assertAlmostEqual(rtime_to_ri(0.2 * 60, m), 900)     # beyond cap, held at bound
        self.assertAlmostEqual(rtime_to_ri(6.0 * 60, m), 1250)
        self.assertAlmostEqual(rtime_to_ri(7.0 * 60, m), 1300)
        self.assertAlmostEqual(rtime_to_ri(9.0 * 60, m), 1300)

    def test_batches_median(self):
        m = self._model([(1000, 2.0), (1000, 2.2), (1000, 2.1), (1100, 3.0), (1100, 3.2), (1200, 5.0)])
        np.testing.assert_allclose(m[1:-1], [[2.1, 1000], [3.1, 1100], [5.0, 1200]])

    def test_non_monotonic_dropped(self):
        m = self._model([(1000, 2.0), (1100, 3.0), (1200, 2.5), (1300, 4.0)])
        np.testing.assert_allclose(m[1:-1, 1], [1000, 1100, 1300])

    def test_report(self):
        import io, contextlib
        m = self._model([(1000, 2.0), (1100, 3.0), (1200, 5.0)])
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            report_ri_extrapolation([30, 150, 240, 360, 600], m)
        out = buf.getvalue()
        self.assertIn("1 features elute before", out)
        self.assertIn("2 after the last", out)
        self.assertIn("2 features beyond that set to the bound", out)
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            report_ri_extrapolation([150, 240], m)
        self.assertEqual(buf.getvalue(), "")

    def test_no_model(self):
        self.assertIsNone(rtime_to_ri(100, None))


if __name__ == '__main__':
    unittest.main()
