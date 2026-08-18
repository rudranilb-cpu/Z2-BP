import unittest

import numpy as np
import pandas as pd

from analysis.common import fit_velocity, nearest_steps


class AnalysisTests(unittest.TestCase):
    def test_velocity_fit(self):
        time = np.linspace(0.0, 1.0, 11)
        frame = pd.DataFrame({"time": time, "front_quantile": 0.5 + 1.7 * time})
        fit = fit_velocity(frame, tmin=0.1, tmax=0.9)
        self.assertAlmostEqual(fit["velocity"], 1.7, places=12)
        self.assertAlmostEqual(fit["intercept"], 0.5, places=12)
        self.assertAlmostEqual(fit["r_squared"], 1.0, places=12)

    def test_nearest_steps(self):
        frame = pd.DataFrame({"step": [0, 1, 2, 3], "time": [0.0, 0.1, 0.2, 0.3]})
        self.assertEqual(nearest_steps(frame, [0.02, 0.19, 0.31]), [0, 2, 3])


if __name__ == "__main__":
    unittest.main()
