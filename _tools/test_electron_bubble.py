"""Junction and causal-domain checks for the retained rotating geometric bubble.

The bulk tensor identities and source matrices are checked by the chart generator
and verify_metrics; these tests concern the interface and the global restrictions.
"""
import math
import unittest
from pathlib import Path
import sys

import numpy as np
from scipy.integrate import quad

sys.path.insert(0, str(Path(__file__).parent / "derivations"))


def metric(r, theta, m, q, a, exterior):
    """Equation (108), evaluated independently of the printed chart expressions."""
    s2 = np.sin(theta)**2
    sigma, b = r*r+a*a*np.cos(theta)**2, r*r+a*a
    f2 = 2*m*r-q*q if exterior else 0
    g = np.diag(np.array([-1, sigma/(b-f2), sigma, b*s2], dtype=np.result_type(r, float)))
    v = np.array([1, 0, 0, -a*s2])
    return g+f2/sigma*np.outer(v, v)


class Junction(unittest.TestCase):
    def test_continuity_and_israel_stress(self):
        for m, q, a in ((1, math.sqrt(2), 2), (2, -3, 4), (1, 2, 0)):
            r0, b = q*q/(2*m), (q*q/(2*m))**2+a*a
            for theta in (.21, .73, 1.3, 1.57, 2.1, 2.9):
                with self.subTest(m=m, q=q, a=a, theta=theta):
                    sigma = r0*r0+a*a*math.cos(theta)**2
                    outside, inside = [metric(r0, theta, m, q, a, side) for side in (True, False)]
                    np.testing.assert_allclose(outside, inside, atol=1e-13)
                    tangent = np.ix_([0, 2, 3], [0, 2, 3])
                    h = inside[tangent]
                    h_inv = np.linalg.inv(h)
                    # Complex differentiation of each independent side, outward normals
                    # both toward increasing r. A changed sign fails the tension check.
                    jumps = []
                    for side in (True, False):
                        dg = np.imag(metric(r0+1e-25j, theta, m, q, a, side))/1e-25
                        jumps.append(.5*math.sqrt(b/sigma)*dg[tangent])
                    jump = jumps[0]-jumps[1]
                    trace = np.sum(h_inv*jump)
                    surface = -(jump-h*trace)/(8*math.pi)
                    u = np.array([b, 0, a])/math.sqrt(b*sigma)
                    lower = h@u
                    tension = m/(8*math.pi*math.sqrt(b*sigma))
                    np.testing.assert_allclose(surface, -tension*(h+np.outer(lower, lower)), atol=3e-14)
                    self.assertAlmostEqual(u@h@u, -1)
                    self.assertAlmostEqual(u@surface@u, 0)
                    eigenvalues = np.sort(np.linalg.eigvals(h_inv@surface).real)
                    np.testing.assert_allclose(eigenvalues, [-tension, -tension, 0], atol=2e-14)
                    # Independent thin-limit p=-m delta(r-r0)/(8 pi Sigma0),
                    # converted to proper distance d ell=sqrt(Sigma0/B0) dr.
                    self.assertAlmostEqual(m/(8*math.pi*sigma)*math.sqrt(sigma/b), tension)
                    speed = a*math.sin(theta)/math.sqrt(b)
                    self.assertLess(abs(speed), 1)

    def test_retained_domain_has_global_time(self):
        for m, q, a in ((1, 1.01, 10), (2, -3, 4), (1, 2, 0)):
            r0 = q*q/(2*m)
            for r in r0*np.geomspace(1, 100, 15):
                for theta in (.03, .6, 1.4, 2.8):
                    g = metric(r, theta, m, q, a, True)
                    self.assertLess(np.linalg.inv(g)[0, 0], 0)
                    self.assertGreater(g[3, 3], 0)
                    self.assertGreater((r-m)**2+a*a+q*q-m*m, 0)
            self.assertGreater(r0, 0)

    def test_homogeneous_maxwell_jump_and_electric_charge(self):
        for m, q, a in ((1, 2, 3), (2, -3, 4)):
            r0, b = q*q/(2*m), (q*q/(2*m))**2+a*a
            for theta in (.31, .8, 1.2, 2.3):
                z = theta+1e-25j
                sigma_z = r0*r0+a*a*np.cos(z)**2
                derivative = np.imag(np.array([-q*r0/sigma_z, q*r0*a*np.sin(z)**2/sigma_z]))/1e-25
                sigma = r0*r0+a*a*math.cos(theta)**2
                coefficient = 2*q*r0*a*math.sin(theta)*math.cos(theta)/sigma**2
                np.testing.assert_allclose(derivative, coefficient*np.array([-a, b]), rtol=1e-13)
                self.assertGreater(np.linalg.norm(derivative), 0)
                self.assertAlmostEqual(derivative@np.array([b, a]), 0)
            # Gauss charge on a t=constant section, including the Lorentz factor
            # of the material relative to its flat-interior inertial observer.
            integral = quad(lambda mu: q*b*(r0*r0-a*a*mu*mu)/(2*(r0*r0+a*a*mu*mu)**2), -1, 1)[0]
            self.assertAlmostEqual(integral, q)


if __name__ == "__main__":
    unittest.main()
