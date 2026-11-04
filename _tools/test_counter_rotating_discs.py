"""The junction, density and stream interpretation of the reflected vacuum fields."""
import math
import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent / "derivations"))
import counter_rotating_discs as discs


class CounterRotatingDiscs(unittest.TestCase):
    def test_reflection_preserves_induced_geometry_and_reverses_normal_derivative(self):
        for kind in ("curzon", "rod", "kerr"):
            for rho in (0.2, 1.0, 5.0):
                metric = discs.seed_matrix(kind, rho, 0.0)
                upper = discs.seed_matrix(kind, rho, 1e-7)
                lower = discs.seed_matrix(kind, rho, 1e-7)  # z -> -z on the lower half
                np.testing.assert_array_equal(upper[:3, :3], lower[:3, :3])
                # Their one-sided slopes with respect to the same oriented coordinate z.
                derivative = discs.seed_matrix(kind, rho, 1j*1e-25).imag/1e-25
                np.testing.assert_allclose((upper-metric)/1e-7, derivative, rtol=3e-6, atol=1e-8)
                np.testing.assert_allclose((lower-metric)/(-1e-7), -derivative, rtol=3e-6, atol=1e-8)

    def test_israel_jump_equals_static_density_pressure_and_zero_radial_stress(self):
        for b in (1.0, 3.0, 8.0):
            for rho in (0.02, 0.5, 2.0, 10.0):
                metric, S = discs.seam_stress("curzon", rho, b=b)
                energy, pressure, _, _ = discs.static_matter(rho, b=b)
                self.assertAlmostEqual(S[0, 0]/(-metric[0, 0]), energy, places=13)
                self.assertAlmostEqual(S[2, 2]/metric[2, 2], pressure, places=13)
                self.assertAlmostEqual(S[1, 1], 0.0, places=13)

    def test_stream_speeds_follow_the_circular_geodesic_equation(self):
        for rho in np.geomspace(0.02, 30, 30):
            metric = discs.seed_matrix("curzon", rho, 0.0)
            dr = discs.seed_matrix("curzon", rho+1j*1e-25, 0.0).imag/1e-25
            omega2 = -dr[0, 0]/dr[2, 2]
            speed2 = omega2*metric[2, 2]/(-metric[0, 0])
            energy, pressure, speed, proper = discs.static_matter(rho)
            self.assertAlmostEqual(speed2, pressure/energy, places=13)
            self.assertAlmostEqual(speed*speed, speed2, places=13)
            self.assertGreaterEqual(speed2, 0)
            self.assertLess(speed2, 1)
            self.assertGreater(proper, 0)

    def test_two_streams_reconstruct_the_full_surface_tensor(self):
        for rho in (0.1, 1.0, 4.0, 20.0):
            metric, covariant = discs.seam_stress("curzon", rho)
            inv = np.linalg.inv(metric[:3, :3])
            target = inv @ covariant @ inv
            energy, pressure, speed, proper = discs.static_matter(rho)
            gamma = 1/math.sqrt(1-speed*speed)
            velocities = [np.array([gamma/math.sqrt(-metric[0, 0]), 0,
                                    sign*gamma*speed/math.sqrt(metric[2, 2])]) for sign in (-1, 1)]
            for u in velocities:
                self.assertAlmostEqual(float(u @ metric[:3, :3] @ u), -1, places=13)
            reconstructed = sum(proper*np.outer(u, u) for u in velocities)
            np.testing.assert_allclose(reconstructed, target, atol=1e-14, rtol=1e-12)
            self.assertAlmostEqual(2*proper, energy-pressure, places=14)
            self.assertAlmostEqual(2*proper*gamma*gamma, energy, places=14)
            self.assertAlmostEqual(2*proper*gamma*gamma*speed*speed, pressure, places=14)

    def test_rod_sources_have_positive_geodesic_streams_in_the_drawn_regime(self):
        for delta in (0.7, 1.0, 1.5):
            for rho in (0.1, 1.0, 5.0):
                g, S = discs.seam_stress("rod", rho, a=1, delta=delta, b=4)
                derivative = discs.seed_matrix("rod", rho+1j*1e-25, 0.0, a=1, delta=delta, b=4).imag/1e-25
                speed2 = -derivative[0, 0]/derivative[2, 2]*g[2, 2]/(-g[0, 0])
                energy, pressure = S[0, 0]/(-g[0, 0]), S[2, 2]/g[2, 2]
                self.assertAlmostEqual(speed2, pressure/energy, places=12)
                self.assertTrue(energy > pressure >= 0)
                self.assertLess(speed2, 1)

    def test_kerr_geodesic_streams_require_and_have_positive_unequal_densities(self):
        for a, b in ((0.8, 3.0), (1.2, 5.0)):
            for rho in (0.1, 1, 4, 15):
                g, Scov = discs.seam_stress("kerr", rho, a=a, b=b)
                inv = np.linalg.inv(g[:3, :3])
                S = inv @ Scov @ inv
                dr = discs.seed_matrix("kerr", rho+1j*1e-25, 0, a=a, b=b).imag/1e-25
                roots = np.sort(np.roots([dr[2, 2], 2*dr[0, 2], dr[0, 0]]))
                weights = np.linalg.solve(np.array([[1, 1], roots]), S[0, [0, 2]])
                reconstruction = np.zeros((3, 3))
                densities = []
                for omega, weight in zip(roots, weights):
                    norm = g[0, 0]+2*g[0, 2]*omega+g[2, 2]*omega*omega
                    self.assertLess(norm, 0)
                    density = -weight*norm
                    self.assertGreater(density, 0)
                    densities.append(density)
                    u = np.array([1, 0, omega])/math.sqrt(-norm)
                    reconstruction += density*np.outer(u, u)
                np.testing.assert_allclose(reconstruction, S, rtol=2e-10, atol=1e-13)
                self.assertNotAlmostEqual(densities[0], densities[1], delta=1e-10)

    def test_kerr_central_density_matches_the_1993_expression(self):
        for a, b in ((0.0, 2.0), (0.8, 3.0), (1.2, 5.0)):
            g, S = discs.seam_stress("kerr", 1e-4, a=a, b=b)
            eigenvalues = np.linalg.eigvals(np.linalg.solve(g[:3, :3], S))
            energy = -min(eigenvalues)
            expected = ((b+1)**2-a*a)/(2*math.pi*math.sqrt(b*b+a*a-1)*((b+1)**2+a*a)**1.5)
            self.assertAlmostEqual(energy, expected, delta=1e-9)


if __name__ == "__main__":
    unittest.main()
