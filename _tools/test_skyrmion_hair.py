"""Boundary conditions and field equations of the spherical skyrmion black hole."""
import json
import math
import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / '_tools' / 'derivations'))
import skyrmion_hair as hole  # noqa: E402


class SkyrmionExterior(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.f_h, cls.solution, cls.logn_far = hole.solve()

    def test_horizon_and_exterior_baryon_number(self):
        f1, m1 = hole.horizon_data(self.f_h)
        self.assertTrue(0 < self.f_h < math.pi)
        self.assertLess(2*m1, 1)
        self.assertLess(f1, 0)
        self.assertEqual(hole.baryon_number(0), 0)
        self.assertAlmostEqual(hole.baryon_number(math.pi), 1, places=15)
        self.assertTrue(0 < hole.baryon_number(self.f_h) < 1)
        self.assertEqual(float(hole.jet(hole.HORIZON)['C']), 0)

    def test_boundary_conditions_and_numerical_convergence(self):
        y = self.solution.y[:, -1]
        self.assertLess(abs(hole.FAR*y[1] + 2*y[0]), 1e-8)
        self.assertLess(abs(y[0]), 0.01)
        self.assertAlmostEqual(float(hole.jet(hole.FAR)['N']), 1, places=12)
        for far, epsilon in ((30.0, 1e-5), (120.0, 1e-7)):
            f_h, result, _ = hole.solve(far, epsilon)
            self.assertLess(abs(f_h - self.f_h), 3e-7)
            self.assertLess(abs(result.y[2, -1] - y[2]), 5e-6)
        x = np.geomspace(0.101, hole.FAR, 400)
        f, p, mu, _ = self.solution.sol(x)
        self.assertTrue(np.all(f > 0) and np.all(p < 0))
        self.assertTrue(np.all(1 - 2*mu/x > 0))
        self.assertTrue(np.all(np.diff(mu) > 0))

    def test_published_metric_and_numeric_jets(self):
        from verify_metrics import Reader
        path = ROOT / 'MFS/assets/data/metrics/skyrmion_hair.json'
        chart = json.loads(path.read_text())['coordinates'][0]
        self.assertEqual(chart['id'], 'areal')
        reader = Reader(chart['coords'], [v['symbol'] for v in chart['parameters']], ())
        g = {tuple(v['indices']): reader(v['value']) for v in chart['metric_components']}
        r, n, c = reader.symbol['r'], reader.parameters['N'], reader.parameters['C']
        point = {r: 2, n: 0.8, c: 0.7}
        self.assertAlmostEqual(float(g[('t', 't')].subs(point)), -0.8**2*0.7)
        self.assertAlmostEqual(float(g[('r', 'r')].subs(point)), 1/0.7)
        x = np.array([0.15, 0.3, 0.8, 2.0, 4.0])
        h = 1e-5*x
        left, middle, right = hole.jet(x-h), hole.jet(x), hole.jet(x+h)
        for name in ('N', 'C'):
            for base, derivative in ((name, 'd'+name), ('d'+name, 'd2'+name)):
                numerical = (right[base] - left[base])/(2*h)
                scale = np.maximum(1e-4, np.abs(middle[derivative]))
                self.assertLess(float(np.max(np.abs(numerical-middle[derivative])/scale)), 2e-5)

    def test_field_equation_in_divergence_form(self):
        # This checks the flux form independently of the second-order equation used to shoot.
        for x in (0.2, 0.4, 0.8, 1.5, 3.0):
            h = 1e-5*x
            def flux(v):
                f, p, mu, logn = self.solution.sol(v)
                c = 1-2*mu/v
                return math.exp(logn)*(v*v + 8*math.sin(f)**2)*c*p
            f, p, mu, logn = self.solution.sol(x)
            c = 1-2*mu/x
            rhs = math.exp(logn)*math.sin(2*f)*(1+4*c*p*p+4*math.sin(f)**2/x**2)
            lhs = (flux(x+h)-flux(x-h))/(2*h)
            self.assertAlmostEqual(lhs, rhs, delta=2e-6*max(1, abs(rhs)))


if __name__ == '__main__':
    unittest.main()
