"""Thin surface stresses and the flat massless exterior of the relativistic discs."""
import json
import math
from pathlib import Path
import sys
import unittest

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT/"_tools"/"derivations"))
import lynden_bell_pineault as lbp
import null_rays as nr

DATA = ROOT/"MFS"/"assets"/"data"


class PublishedDiscs(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        lbp.register_dimensions()
        cls.metric = json.loads((DATA/"metrics"/f"{lbp.METRIC}.json").read_text())

    def test_the_surface_is_separate_from_every_smooth_chart(self):
        self.assertEqual({c["id"] for c in self.metric["coordinates"]}, set(lbp.SYSTEMS))
        self.assertNotIn("vacuum", self.metric["tags"])
        self.assertIn("equatorial sheet", self.metric["convention"])
        for chart in self.metric["coordinates"]:
            domain = chart["domains"][2]
            self.assertIn("(0,", domain)
        diagrams = json.loads((DATA/"diagrams"/f"{lbp.METRIC}.json").read_text())
        self.assertEqual(set(diagrams["systems"]), set(lbp.SYSTEMS))

    def test_the_massless_vacuum_chart_is_a_rindler_region(self):
        _, entry, reader = nr.load(lbp.METRIC, "photons")
        g = nr.published_matrix(reader, entry, "metric_components")
        t, r, theta, phi = [reader.symbol[c] for c in entry["coords"]]
        A = reader.parameters["A"]
        rho, zeta = 2*sp.sqrt(A*r)*sp.sin(theta/2), 2*sp.sqrt(A*r)*sp.cos(theta/2)
        # Coordinates in Minkowski space, with x^0=ct in the published component array.
        image = sp.Matrix([zeta*sp.sinh(t/(2*A)), rho*sp.cos(phi), rho*sp.sin(phi), zeta*sp.cosh(t/(2*A))])
        J = image.jacobian([t, r, theta, phi])
        pulled = J.T*sp.diag(-1, 1, 1, 1)*J
        for radius, angle, azimuth in [(0.7, 0.4, 0.2), (2.1, 1.1, 1.4), (1.3, 0.8, -0.3)]:
            point = {t: 0.6, r: radius, theta: angle, phi: azimuth, A: 1.7}
            self.assertLess(max(abs(float(v)) for v in (pulled-g).subs(point)), 1e-12)
        self.assertEqual(entry["riemann"]["variants"]["llll"]["nonzero"], [])
        self.assertEqual(entry["kretschmann"], "K = 0")

    def surface_stress(self, system):
        """Israel's reflection jump, up to its common positive dimensional multiplier."""
        _, entry, reader = nr.load(lbp.METRIC, system)
        g = nr.published_matrix(reader, entry, "metric_components")
        theta = reader.symbol["\\theta"]
        point = {reader.symbol["r"]: 1.3, reader.parameters["A"]: 1, theta: math.pi/2}
        if system == "rotating":
            solution = lbp.rotating_solution()
            point[reader.parameters["n"]] = solution["n"]
            jets = lbp.rotating_jets(reader, math.pi/2, solution)
        else:
            point[reader.parameters["n"]] = 0.2
            jets = {}
        h = np.array(g.extract([0, 1, 3], [0, 1, 3]).xreplace(jets).subs(point), dtype=float)
        derivative = g.diff(theta).extract([0, 1, 3], [0, 1, 3])
        normal = math.sqrt(float(g[2, 2].xreplace(jets).subs(point)))
        K = np.array(derivative.xreplace(jets).subs(point), dtype=float)/(2*normal)
        stress = K-np.trace(np.linalg.solve(h, K))*h
        lapse = math.sqrt(-h[0, 0]+h[0, 2]**2/h[2, 2])
        omega = -h[0, 2]/h[2, 2]
        e0 = np.array([1, 0, omega])/lapse
        ephi = np.array([0, 0, 1])/math.sqrt(h[2, 2])
        return stress, float(e0@stress@e0), float(e0@stress@ephi), float(ephi@stress@ephi)

    def test_counter_rotating_streams_have_azimuthal_stress_and_no_radial_pressure(self):
        stress, energy, flux, pressure = self.surface_stress("spherical")
        self.assertGreater(energy, 0)
        self.assertAlmostEqual(stress[1, 1], 0, delta=1e-12)
        self.assertAlmostEqual(flux, 0, delta=1e-12)
        self.assertAlmostEqual(pressure/energy, 0.25, delta=1e-12)

    def test_the_rotating_surface_is_one_cold_stream_at_the_declared_speed(self):
        stress, energy, flux, pressure = self.surface_stress("rotating")
        self.assertGreater(energy, 0)
        self.assertAlmostEqual(stress[1, 1], 0, delta=2e-10)
        self.assertAlmostEqual(-flux/energy, 0.2, delta=2e-9)
        self.assertAlmostEqual(pressure/energy, 0.2**2, delta=2e-9)
        self.assertLess(abs(energy*pressure-flux*flux), 2e-10)


if __name__ == "__main__":
    unittest.main()
