"""Tests for LCMSMS precursor features (constructors.nearest_peak_intensity)."""
import unittest
import numpy as np
from asari.constructors import nearest_peak_intensity


class TestNearestPeakIntensity(unittest.TestCase):
    rtimes = np.arange(0, 20, 1.0)                  # 1 s per scan
    intensity = np.zeros(20)
    intensity[3:8] = [10, 50, 100, 50, 10]          # peak at 5 s
    intensity[11:16] = [20, 80, 300, 80, 20]        # bigger peak at 13 s

    def test_nearest_not_largest(self):
        self.assertEqual(nearest_peak_intensity(self.intensity, self.rtimes, 7.0, window=10), 100)
        self.assertEqual(nearest_peak_intensity(self.intensity, self.rtimes, 10.5, window=10), 300)

    def test_window_and_no_signal(self):
        self.assertEqual(nearest_peak_intensity(self.intensity, self.rtimes, 7.0, window=1), 50)   # rising edge only
        self.assertEqual(nearest_peak_intensity(self.intensity, self.rtimes, 18.5, window=1), 0)
        self.assertEqual(nearest_peak_intensity(self.intensity, self.rtimes, 50.0, window=5), 0)


if __name__ == '__main__':
    unittest.main()
