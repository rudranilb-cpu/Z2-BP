import math
import unittest

import numpy as np

from validation.exact_reference import (
    ExactConfig,
    continuous_evolve,
    evolve,
    initial_states,
    link_fluxes,
    links,
    phase_aligned_error,
    total_flux,
    with_dt,
)


class ExactReferenceTests(unittest.TestCase):
    def test_link_counts_and_boundary_convention(self):
        open_cfg = ExactConfig(nx=8, ny=6)
        full_cfg = ExactConfig(nx=8, ny=6, boundary="fixed_exterior")
        self.assertEqual(len(links(open_cfg)), 82)
        self.assertEqual(len(links(full_cfg)), 110)

    def test_bulk_single_flip_has_delta_flux_four(self):
        cfg = ExactConfig(nx=3, ny=3, defect=(2, 2))
        scv, loop = initial_states(cfg)
        self.assertAlmostEqual(total_flux(scv, cfg), 0.0, places=14)
        self.assertAlmostEqual(total_flux(loop, cfg) - total_flux(scv, cfg), 4.0, places=14)

    def test_first_lie_step_analytic_delta(self):
        cfg = ExactConfig(nx=3, ny=3, K=3.0, Gamma=1.0, dt=0.05, defect=(2, 2))
        scv, loop = initial_states(cfg)
        scv = evolve(scv, cfg, 1)
        loop = evolve(loop, cfg, 1)
        expected = 4.0 * math.cos(2.0 * cfg.K * cfg.dt) ** 2
        self.assertAlmostEqual(total_flux(loop, cfg) - total_flux(scv, cfg), expected, places=12)
        self.assertTrue(np.all(link_fluxes(scv, cfg) >= -1e-13))
        self.assertTrue(np.all(link_fluxes(loop, cfg) <= 1 + 1e-13))

    def test_trotter_order_differs_after_two_steps(self):
        base = ExactConfig(nx=3, ny=3, K=1.7, Gamma=0.9, dt=0.08, defect=(2, 2))
        _, loop = initial_states(base)
        xzz = evolve(loop, with_dt(base, base.dt, "lie_x_zz"), 2)
        zzx = evolve(loop, with_dt(base, base.dt, "lie_zz_x"), 2)
        self.assertGreater(phase_aligned_error(xzz, zzx), 1e-4)

    def test_lie_and_strang_convergence_orders(self):
        base = ExactConfig(nx=2, ny=2, K=1.3, Gamma=0.7, defect=(1, 1))
        _, loop = initial_states(base)
        final_time = 0.4
        reference = continuous_evolve(loop, base, final_time)

        lie_errors = []
        strang_errors = []
        for dt in (0.1, 0.05, 0.025):
            steps = round(final_time / dt)
            lie_errors.append(phase_aligned_error(reference, evolve(loop, with_dt(base, dt, "lie_x_zz"), steps)))
            strang_errors.append(phase_aligned_error(reference, evolve(loop, with_dt(base, dt, "strang"), steps)))

        self.assertGreater(lie_errors[0] / lie_errors[1], 1.7)
        self.assertGreater(lie_errors[1] / lie_errors[2], 1.7)
        self.assertGreater(strang_errors[0] / strang_errors[1], 3.2)
        self.assertGreater(strang_errors[1] / strang_errors[2], 3.2)
        self.assertLess(strang_errors[-1], lie_errors[-1])


if __name__ == "__main__":
    unittest.main()
