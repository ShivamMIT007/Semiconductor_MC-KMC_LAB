import unittest
import numpy as np
from semimc.models.implantation import simulate as implant
from semimc.models.dopant_diffusion import simulate as diffusion, build_diamond_lattice
from semimc.models.thin_film import simulate as film
from semimc.models.grain_growth import simulate as grain
from semimc.models.alloy import simulate_ordering, simulate_phase_separation
from semimc.models.defect_annealing import simulate as defects
from semimc.models.plasma_etch import simulate as etch
from semimc.models.oxidation import simulate as oxidation
from semimc.models.photolithography import simulate as lithography


class TestSimulationModels(unittest.TestCase):
    def test_ion_implantation_histogram_and_seed(self):
        p = {"n_ions": 800, "bins": 32, "seed": 17}
        a = implant(p); b = implant(p)
        self.assertTrue(np.allclose(a.after, b.after))
        self.assertGreater(a.metrics["sampled_mean_range_nm"], 0)
        self.assertLessEqual(a.metrics["fraction_stopped_within_target"], 1)

    def test_diamond_lattice_four_neighbours_and_reciprocal(self):
        coords, neighbours, deltas, pos, L = build_diamond_lattice(2, 0.5431)
        self.assertEqual(len(coords), 64)
        self.assertTrue(np.all((neighbours >= 0) & (neighbours < len(coords))))
        for i in range(len(coords)):
            self.assertEqual(len(set(neighbours[i].tolist())), 4)
            for j in neighbours[i]:
                self.assertIn(i, neighbours[j])

    def test_dopant_diffusion_kmc_finite(self):
        r = diffusion({"n_cells": 2, "events_per_run": 300, "replicates": 2, "seed": 8})
        self.assertTrue(np.isfinite(r.metrics["mean_tracer_diffusivity_m2_s"]))
        self.assertEqual(len(r.data), 100)
        self.assertTrue(np.all(np.diff(r.data["mean_physical_time_s"]) >= 0))

    def test_thin_film_kmc_time_and_height(self):
        r = film({"size": 8, "events": 250, "seed": 8})
        self.assertGreaterEqual(r.metrics["final_mean_thickness_ML"], 0)
        self.assertTrue(np.all(np.diff(r.data["time_s"]) >= 0))

    def test_grain_growth_potts(self):
        r = grain({"size": 16, "initial_grains": 6, "mcs": 3, "seed": 9})
        self.assertEqual(r.before.shape, (16, 16))
        self.assertGreaterEqual(r.metrics["final_grain_boundary_fraction"], 0)
        self.assertLessEqual(r.metrics["final_grain_boundary_fraction"], 1)

    def test_alloy_composition_conserved(self):
        for fn in (simulate_ordering, simulate_phase_separation):
            r = fn({"size": 12, "mcs": 3, "fraction_A": 0.5, "seed": 10})
            self.assertAlmostEqual(r.metrics["initial_A_fraction"], r.metrics["final_A_fraction"])

    def test_defect_population_does_not_increase(self):
        r = defects({"size": 12, "initial_pairs": 10, "events": 250, "seed": 11})
        self.assertLessEqual(r.metrics["final_surviving_defects"], r.metrics["initial_defect_objects"])

    def test_plasma_etch_output_finite(self):
        r = etch({"columns": 32, "ions": 400, "seed": 12})
        self.assertTrue(np.all(np.isfinite(r.after)))
        self.assertEqual(r.after.shape, (1, 32))

    def test_oxidation_monotonic(self):
        r = oxidation({"points": 30, "duration_s": 1000})
        self.assertTrue(np.all(np.diff(r.data["oxide_thickness_nm"]) >= 0))

    def test_photolithography_dimensions(self):
        r = lithography({"photons": 2000, "pixels_x": 24, "pixels_z": 12, "seed": 13})
        self.assertEqual(r.after.shape, (12, 24))
        self.assertGreaterEqual(r.metrics["fraction_absorbed_in_resist"], 0)
        self.assertLessEqual(r.metrics["fraction_absorbed_in_resist"], 1)


if __name__ == "__main__":
    unittest.main()
