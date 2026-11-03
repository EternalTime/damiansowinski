"""The two deficits and the regular axis of the cylindrical wave geometry."""

import json
import sys
import unittest
from pathlib import Path

import mpmath as mp
import sympy as sp

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / '_tools' / 'derivations'))
import verify_metrics as vm


class StringCylindricalWaves(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.metric = json.loads((ROOT / 'MFS/assets/data/metrics/string_cylindrical_waves.json').read_text())
        cls.chart = next(c for c in cls.metric['coordinates'] if c['id'] == 'prolate')
        cls.reader = vm.Reader(cls.chart['coords'], [p['symbol'] for p in cls.chart['parameters']], ())
        cls.g = vm.metric_from_line_element(cls.reader, cls.chart['line_element'], cls.chart['coords'])

    def circumference_ratio(self, rho, p, alpha):
        # At t=0, tau=0 and sigma=sqrt(1+rho^2).  The radial metric in rho is
        # g_sigma_sigma * (rho/sigma)^2, so this tends to circumference/(2*pi*proper radius).
        with mp.workdps(70):
            rho, p, alpha = (mp.mpf(str(v)) for v in (rho, p, alpha))
            sigma = mp.sqrt(1+rho*rho)
            R = self.reader
            at = {R.symbol[r'\tau']: 0, R.symbol[r'\sigma']: sp.Float(str(sigma), 70),
                  R.parameters['p']: sp.Float(str(p), 70), R.parameters['alpha']: sp.Float(str(alpha), 70),
                  R.parameters['ell']: 1}
            grr = mp.mpf(str(self.g[1,1].subs(at).evalf(60))) * (rho/sigma)**2
            gpp = mp.mpf(str(self.g[2,2].subs(at).evalf(60)))
            return mp.sqrt(gpp/grr)/rho

    def test_axis_deficit_and_regular_axis(self):
        with mp.workdps(60):
            for p, alpha in [('1', '2'), ('2', '1.25'), ('2', '.5')]:
                ratio = self.circumference_ratio('1e-8', p, alpha)
                expected = 1/(mp.mpf(alpha)*mp.mpf(p))
                self.assertLess(abs(ratio-expected), mp.mpf('1e-13'))
                deficit = 2*mp.pi*(1-ratio)
                self.assertLess(abs(deficit-2*mp.pi*(1-expected)), mp.mpf('1e-12'))
            # alpha*p=1 makes the ratio one and the axis deficit zero.
            self.assertLess(abs(1-self.circumference_ratio('1e-8', '2', '.5')), mp.mpf('1e-13'))

    def test_far_deficit_at_fixed_time(self):
        with mp.workdps(60):
            for p, alpha in [('1', '2'), ('2', '1.25')]:
                ratio = self.circumference_ratio('1e8', p, alpha)
                expected = 1/(mp.mpf(alpha)*mp.sqrt(1+mp.mpf(p)**2))
                self.assertLess(abs(ratio-expected), mp.mpf('1e-13'))
                self.assertLess(abs(2*mp.pi*(1-ratio)-2*mp.pi*(1-expected)), mp.mpf('1e-12'))

    def test_both_charts_state_exterior_domain(self):
        charts = {c['id']: c for c in self.metric['coordinates']}
        self.assertEqual(set(charts), {'prolate', 'hyperbolic'})
        self.assertIn(r'\sigma\in(1,\infty)', charts['prolate']['domains'])
        self.assertIn(r'\theta\in(0,\infty)', charts['hyperbolic']['domains'])


if __name__ == '__main__':
    unittest.main()
