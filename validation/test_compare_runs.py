import unittest

from validation.compare_runs import metadata_mismatches


def metadata():
    return {
        "schema": {"version": "1.0.0"},
        "simulation": {
            "nx": 4, "ny": 4, "K": 3.0, "Gamma": 1.0, "dt": 0.05,
            "steps": 20, "trotter": "lie_x_zz", "boundary": "open_dual",
            "backend": "bp", "defect_x": 2, "defect_y": 2,
        },
        "peps": {"maxdim": 16, "cutoff": 1e-10, "normalize_tensors": True},
        "bp": {"tolerance": 1e-8, "maxiter": 50},
        "observables": {"measure_every": 1, "front_quantile": 0.9},
        "conventions": {"differential_sign": "loop_minus_scv"},
    }


class MetadataComparisonTests(unittest.TestCase):
    def test_only_dimension_may_change_in_dimension_channel(self):
        reference, candidate = metadata(), metadata()
        candidate["peps"]["maxdim"] = 24
        self.assertEqual(metadata_mismatches(reference, candidate, "bond_dimension"), [])
        candidate["simulation"]["K"] = 4.0
        self.assertTrue(metadata_mismatches(reference, candidate, "bond_dimension"))

    def test_bp_controls_are_isolated(self):
        reference, candidate = metadata(), metadata()
        candidate["bp"].update(tolerance=1e-10, maxiter=100)
        self.assertEqual(metadata_mismatches(reference, candidate, "bp_environment"), [])
        candidate["peps"]["maxdim"] = 24
        self.assertTrue(metadata_mismatches(reference, candidate, "bp_environment"))

    def test_trotter_requires_common_final_time(self):
        reference, candidate = metadata(), metadata()
        candidate["simulation"].update(dt=0.025, steps=40, trotter="strang")
        self.assertEqual(metadata_mismatches(reference, candidate, "trotter"), [])
        candidate["simulation"]["steps"] = 39
        self.assertTrue(any("tmax" in item for item in metadata_mismatches(reference, candidate, "trotter")))

    def test_finite_size_is_a_separate_channel(self):
        reference, candidate = metadata(), metadata()
        candidate["simulation"].update(nx=6, ny=6, defect_x=3, defect_y=3)
        candidate["conventions"].update(internal_bonds=60, measured_links=60,
                                        raw_pauli_observables=120, boundary_distance_sites=2)
        reference["conventions"].update(internal_bonds=24, measured_links=24,
                                        raw_pauli_observables=60, boundary_distance_sites=1)
        self.assertEqual(metadata_mismatches(reference, candidate, "finite_size"), [])
        candidate["simulation"]["K"] = 4.0
        self.assertTrue(metadata_mismatches(reference, candidate, "finite_size"))


if __name__ == "__main__":
    unittest.main()
