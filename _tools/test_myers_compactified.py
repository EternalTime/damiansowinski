"""Global identification and geometric limits of the compactified black hole."""
import json
import math
from pathlib import Path
import unittest

import sympy as sp

ROOT = Path(__file__).resolve().parent.parent


def potential(r, z, mu, length):
    q = 2*math.pi/length
    return 1 + math.pi*mu/(length*r)*math.sinh(q*r)/(math.cosh(q*r)-math.cos(q*z))


class CompactifiedBlackHole(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        page = json.loads((ROOT / 'MFS/assets/data/metrics/myers_compactified.json').read_text())
        cls.charts = {chart['id']: chart for chart in page['coordinates']}

    def test_circle_identification_and_sourced_charts(self):
        self.assertEqual(set(self.charts), {'cylindrical', 'spherical'})
        self.assertIn(r'z \in [-L/2,L/2)', self.charts['cylindrical']['domains'])
        self.assertIn(r'-L/2 < \rho\cos\chi < L/2', self.charts['spherical']['domains'])
        for r, z in ((0.3, 0.2), (1.2, -0.7)):
            self.assertAlmostEqual(potential(r, z, 1, 4), potential(r, z + 4, 1, 4), places=13)

    def test_potential_is_periodic_array_of_extremal_holes(self):
        r, z, mu, length = 0.8, 0.3, 1.0, 4.0
        images = sum(mu/(r*r + (z + n*length)**2) for n in range(-10000, 10001))
        self.assertAlmostEqual(potential(r, z, mu, length), 1 + images, delta=1e-4)

    def test_spherical_spatial_metric_is_cylindrical_pullback(self):
        rho, chi = sp.symbols('rho chi', real=True)
        r, z = rho*sp.sin(chi), rho*sp.cos(chi)
        jac = sp.Matrix([[sp.diff(r, rho), sp.diff(r, chi)],
                         [sp.diff(z, rho), sp.diff(z, chi)]])
        self.assertEqual(sp.simplify(jac.T*jac), sp.diag(1, rho**2))
        self.assertEqual(sp.simplify(r*r), rho**2*sp.sin(chi)**2)

    def test_circle_length_is_written_the_same_way_in_every_view(self):
        data = ROOT / 'MFS/assets/data'
        diagrams = json.loads((data / 'diagrams/myers_compactified.json').read_text())
        views = [view for group in diagrams['systems'].values() for view in group]
        for folder in ('embedding', 'curvature'):
            views.extend(json.loads((data / folder / 'myers_compactified.json').read_text())['views'])
        self.assertEqual(len(views), 5)
        for view in views:
            self.assertIn(r'$L = 4\sqrt{\mu}$', view['settings'])
            if view.get('input'):
                self.assertIn(r'$L = 4\sqrt{\mu}$', view['input'])

    def test_horizon_area_limit(self):
        mu, length = 1.7, 4.0
        for chi in (0.3, 1.2, 2.6):
            rho = 1e-4
            u = potential(rho*math.sin(chi), rho*math.cos(chi), mu, length)
            radius_squared = u*rho*rho
            self.assertAlmostEqual(radius_squared, mu, delta=1e-7)
            self.assertAlmostEqual(2*math.pi**2*radius_squared**1.5,
                                   2*math.pi**2*mu**1.5, delta=1e-5)


if __name__ == '__main__':
    unittest.main()
