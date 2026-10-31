"""The homogeneous stability ranges and parameter-dependent global causal domains.

The chart generator tests the matter equations and coordinate pullbacks.
These checks cover the phase-plane stability claims and improper null integrals.
"""
import json
import math
from pathlib import Path
import sys
import unittest

import numpy as np
from scipy.integrate import quad
import sympy as sp

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / '_tools' / 'derivations'))
import verify_metrics as vm


class ExponentialScalar(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.page = json.loads((ROOT/'MFS/assets/data/metrics/exponential_scalar_universes.json').read_text())
        cls.charts = {x['id']: x for x in cls.page['coordinates']}
        # Copeland, Liddle and Wands (1998), equations (6)-(7).
        x, y, lam, gamma = sp.symbols('x y lambda gamma', real=True)
        common = 2*x*x+gamma*(1-x*x-y*y)
        flow = sp.Matrix([-3*x+sp.sqrt(sp.Rational(3, 2))*lam*y*y+sp.Rational(3, 2)*x*common,
                          -sp.sqrt(sp.Rational(3, 2))*lam*x*y+sp.Rational(3, 2)*y*common])
        cls.flow = staticmethod(sp.lambdify((x, y, lam, gamma), flow, 'numpy'))
        cls.jacobian = staticmethod(sp.lambdify((x, y, lam, gamma), flow.jacobian((x, y)), 'numpy'))

    def test_scalar_attraction_depends_on_the_fluid(self):
        for lam2, gamma, stable in ((1, 1, True), (4, 1, False), (5, 1.9, True)):
            lam = math.sqrt(lam2)
            point = (lam/math.sqrt(6), math.sqrt(1-lam2/6), lam, gamma)
            np.testing.assert_allclose(self.flow(*point), 0, atol=1e-12)
            eigenvalues = np.linalg.eigvals(self.jacobian(*point))
            np.testing.assert_allclose(sorted(eigenvalues), sorted(((lam2-6)/2, lam2-3*gamma)), atol=1e-12)
            self.assertEqual(bool(np.max(eigenvalues.real) < 0), stable)
        # At the meeting of the branches one eigenvalue is zero.
        point = (math.sqrt(3/6), math.sqrt(1-3/6), math.sqrt(3), 1)
        self.assertLess(min(abs(np.linalg.eigvals(self.jacobian(*point)))), 1e-12)

    def test_positive_fluid_scaling_is_an_attractor(self):
        # Small gamma also tests the node regime without assuming 9 gamma - 2 > 0.
        for gamma, lam2 in ((1, 9), (4/3, 12), (0.1, 1), (0.5, 2)):
            lam = math.sqrt(lam2)
            x = math.sqrt(1.5)*gamma/lam
            y = math.sqrt(1.5*gamma*(2-gamma))/lam
            np.testing.assert_allclose(self.flow(x, y, lam, gamma), 0, atol=1e-12)
            self.assertGreater(1-x*x-y*y, 0)
            self.assertTrue(np.all(np.linalg.eigvals(self.jacobian(x, y, lam, gamma)).real < 0))

    def test_cosmic_metric_has_the_three_distinct_null_integrals(self):
        chart = self.charts['cosmic_cartesian']
        reader = vm.Reader(chart['coords'], [x['symbol'] for x in chart['parameters']], ())
        metric = vm.metric_from_line_element(reader, chart['line_element'], chart['coords'])
        time = reader.symbol['t']
        # dt/dx from the published null condition at a_0=t_0=1.
        a = sp.sqrt(metric[1, 1]/(-metric[0, 0]))
        for exponent in (2, 2/3):
            inv_a = sp.lambdify(time, (1/a).subs({reader.parameters['p']: exponent,
                                                reader.parameters['a_0']: 1, reader.parameters['t_0']: 1}), 'numpy')
            if exponent > 1:
                event, _ = quad(inv_a, 1, np.inf)
                self.assertAlmostEqual(event, 1/(exponent-1), places=10)
                self.assertGreater(quad(lambda z: inv_a(math.exp(z))*math.exp(z), math.log(1e-5), 0)[0], 1000)
            else:
                particle, _ = quad(inv_a, 0, 1)
                self.assertAlmostEqual(particle, 1/(1-exponent), places=10)
                self.assertGreater(quad(lambda z: inv_a(math.exp(z))*math.exp(z), 0, math.log(1e12))[0], 1000)
        # Coasting accumulates equal conformal lengths on every logarithmic interval.
        for lo in (1e-6, 1, 1e6):
            self.assertAlmostEqual(quad(lambda t: 1/t, lo, 10*lo)[0], math.log(10), places=12)
        domains = self.charts['conformal']['domains']
        self.assertIn(r'\eta \in (-\infty,0)\;\text{for}\;p>1', domains)
        self.assertIn(r'\eta \in (0,\infty)\;\text{for}\;p<1', domains)
        self.assertIn(r'\eta \in (-\infty,\infty)', self.charts['coasting']['domains'])


    def test_dust_display_covers_the_entire_positive_time_domain(self):
        """Exact cubic null curves bound every segment of the rounded drawing.

        The maximum is at a segment endpoint or where the cubic and segment
        slopes agree. No positive-time interval is excluded from this check.
        """
        view = json.loads((ROOT/'MFS/assets/data/diagrams/exponential_scalar_universes.json').read_text())['systems']['cosmic_spherical'][0]
        radius, time = view['box'][1], view['box'][3]
        self.assertEqual(view['box'][::2], [0, 0])

        def clip(a, b):
            lo, hi = 0., 1.
            for i in (0, 1):
                delta = b[i]-a[i]
                if delta == 0:
                    if not 0 <= a[i] <= 1:
                        return None
                else:
                    low, high = sorted(((0-a[i])/delta, (1-a[i])/delta))
                    lo, hi = max(lo, low), min(hi, high)
            if hi < lo:
                return None
            return ([a[i]+lo*(b[i]-a[i]) for i in (0,1)],
                    [a[i]+hi*(b[i]-a[i]) for i in (0,1)])

        vertical, past, covered, rays = 0., 0., 0, 0
        for family, sign in (('P', 1), ('M', -1)):
            for line in view['rays'][family]:
                segments = [pair for a,b in zip(line,line[1:]) if (pair := clip(a,b))]
                if not segments:
                    continue
                anchors = [p for p in line if 0 <= p[0] <= 1 and 0 < p[1] <= 1]
                anchor = max(anchors or [p for pair in segments for p in pair], key=lambda p:p[1])
                constant = 3*math.cbrt(time*anchor[1])+sign*radius*anchor[0]
                bang = constant/(sign*radius)
                rays += 1
                for a,b in segments:
                    if max(a[1],b[1]) > 0 and min(a[1],b[1]) < 0.02/time:
                        covered += 1
                    dx = b[0]-a[0]
                    if dx == 0:
                        candidates = [a,b]
                    else:
                        slope = (b[1]-a[1])/dx
                        lo, hi = sorted((a[0],b[0]))
                        xs = [lo,hi]
                        squared = -9*time*slope/(sign*radius)
                        if squared >= 0:
                            for eta in (math.sqrt(squared), -math.sqrt(squared)):
                                x = (constant-eta)/(sign*radius)
                                if lo < x < hi:
                                    xs.append(x)
                        if lo < bang < hi:
                            xs.append(bang)
                        candidates = [(x,a[1]+slope*(x-a[0])) for x in xs]
                    for x,y in candidates:
                        eta = constant-sign*radius*x
                        if eta >= 0:
                            vertical = max(vertical,abs(y-eta**3/(27*time)))
                        else:
                            past = max(past, math.hypot(x-bang,y))
        # 0.0006 thinning, four-decimal rounding and the maximal drawn slope
        # (below 1.43) permit 0.0012 in the unit window, below 0.624 pixels
        # even if both axes are conservatively taken as 520 pixels.
        self.assertGreater(rays, 0)
        self.assertGreater(covered, 0)
        self.assertLessEqual(vertical, 0.0012)
        self.assertLessEqual(past, 0.0012)
        print(f'Dust display: {rays} complete rays; {covered} segments meet 0<t<=0.02 t0; '
              f'vertical bound {vertical:.12g}; past-endpoint bound {past:.12g}; '
              f'pixel bound {520*max(vertical,past):.6g}.')


if __name__ == '__main__':
    unittest.main()
