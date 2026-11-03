"""Source limits of the rotating deformed mass in the displayed chart."""
import json
import sys
import unittest
from pathlib import Path

import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / '_tools' / 'derivations'))
import quevedo_mashhoon as qm


class RotatingDeformedMass(unittest.TestCase):
    def setUp(self):
        mp.mp.dps = 40

    def close(self, actual, expected):
        self.assertLess(abs(actual - expected), mp.mpf('1e-30') * max(1, abs(expected)))

    def test_kerr_limit(self):
        for spin in (mp.mpf('.7'), mp.mpf('-.4')):
            mass, x, y = mp.mpf(2), mp.mpf('2.7'), mp.mpf('.3')
            sigma = mp.sqrt(mass**2-spin**2)
            f, w, gamma = qm.potentials(x, y, mass, spin, 0, 1)
            n = sigma**2*x**2 + spin**2*y**2-mass**2
            self.close(f, n / ((mass+sigma*x)**2+spin**2*y**2))
            self.close(w, 2*spin*mass*(mass+sigma*x)*(1-y*y)/n)
            self.close(gamma, n/(sigma**2*(x*x-y*y)))

    def test_zipoy_voorhees_limit(self):
        mass, x, y = mp.mpf('2'), mp.mpf('2.7'), mp.mpf('.3')
        for delta in (mp.mpf('.8'), mp.mpf('1.3')):
            f, w, gamma = qm.potentials(x, y, mass, 0, 0, delta)
            self.close(f, ((x-1)/(x+1))**delta)
            self.close(w, 0)
            self.close(gamma, ((x*x-1)/(x*x-y*y))**(delta**2))
            # These are the moments stated on the live Zipoy-Voorhees page.
            q_zv = delta - 1
            self.close(mass*delta, mass*(1+q_zv))
            self.close(-mass**3*delta*(delta**2-1)/3,
                       -mass**3*q_zv*(1+q_zv)*(2+q_zv)/3)

    def test_chart_domains_and_parameter_meaning(self):
        metric = json.loads((ROOT/'MFS/assets/data/metrics'/f'{qm.METRIC}.json').read_text())
        self.assertEqual({c['id'] for c in metric['coordinates']}, set(qm.SYSTEMS))
        for chart in metric['coordinates']:
            meanings = {p['symbol']: p['description'] for p in chart['parameters']}
            self.assertIn('M>0', meanings['M'])
            self.assertIn('|a|<M', meanings['a'])
            self.assertIn('M_0=M+', meanings['M'])
            self.assertIn('2q\\sigma^3/15', meanings['q'])

    def test_figure_derivatives_and_coordinate_map(self):
        x, y = mp.mpf('2.7'), mp.mpf('.3')
        spin = mp.mpf(2) / 5
        sigma = mp.sqrt(1 - spin**2)
        r, theta = 1 + sigma*x, mp.acos(y)
        prolate = qm.figure_values(x, y, 'prolate_spheroidal', 40)
        spherical = qm.figure_values(r, theta, 'spherical', 40)
        for i, name in enumerate(('f', 'omega', 'h')):
            fn = lambda u, v: qm.potentials(u, v, mp.mpf(1), spin, mp.mpf(2)/5, mp.mpf(1))[i]
            self.close(prolate[name]['1'], mp.diff(lambda u: fn(u, y), x))
            self.close(prolate[name]['12'], mp.diff(fn, (x, y), (1, 1)))
            self.close(spherical[name][''], prolate[name][''])
            self.close(spherical[name]['1'], prolate[name]['1'] / sigma)
            self.close(spherical[name]['2'], -mp.sin(theta) * prolate[name]['2'])


if __name__ == '__main__':
    unittest.main()
