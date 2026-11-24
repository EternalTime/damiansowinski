"""The homogeneous stability ranges and parameter-dependent global causal domains.

The chart generator tests the matter equations and coordinate pullbacks.
These checks cover the phase-plane stability claims and improper null integrals.
The embedding diagram's rings of free particles and horizons are held to the members drawn,
from the numbers in its file.
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

    def test_embedding_is_a_movie_of_a_ring_of_free_particles_for_each_member(self):
        """Each member's ring of free particles at rest in the chart is the circle of proper radius
        R t^p on the flat plane z = 0, R = 1, 3 and 1, from t = t_0/2 to 3t_0/2 a frame every t_0/40.
        The event horizon at p = 2, of proper radius t, and the particle horizon at p = 2/3, of
        proper radius 3t, each meet the ring at t = t_0, and the coasting member has neither."""
        embedding = json.loads((ROOT/'MFS/assets/data/embedding/exponential_scalar_universes.json').read_text())
        views = {view['id']: view for view in embedding['views']}
        self.assertEqual(list(views), ['inflation', 'tracking', 'coasting'])
        members = {'inflation': (2, 1, 'event'), 'tracking': (2/3, 3, 'particle'), 'coasting': (1, 1, None)}
        for vid, (p, ring, horizon) in members.items():
            view = views[vid]
            self.assertEqual([s['time'] for s in view['surfaces']], [0.5, 1, 1.5], vid)
            self.assertEqual(view['movie']['variable'], '$t$', vid)
            frames = view['movie']['frames']
            values = [frame['value'] for frame in frames]
            self.assertEqual(len(frames), 41, vid)
            for earlier, later in zip(values, values[1:]):
                self.assertAlmostEqual(later - earlier, 0.025, delta=1e-6, msg=vid)
            radii, horizons = [], []
            for frame in frames:
                t = frame['value']
                (curve,) = frame['curves']
                self.assertEqual(curve['class'], 'particles', vid)
                for X, Y, Z in curve['points']:
                    self.assertAlmostEqual(math.hypot(X, Y), ring*t**p, delta=2e-6, msg=f'{vid} at {t}')
                    self.assertEqual(Z, 0, vid)
                radii.append(math.hypot(*curve['points'][0][:2]))
                found = [r for r in frame['rings'] if r['class'] == 'horizon']
                if horizon is None:
                    self.assertEqual(found, [], f'{vid} at {t}')
                    continue
                (circle,) = found
                # The comoving radius 1/t of the event horizon, or 3t^(1/3) of the particle horizon, at the proper radius t/|p - 1|.
                self.assertAlmostEqual(circle['x'], t**(1-p)/abs(p-1), delta=1e-9, msg=f'{vid} at {t}')
                self.assertAlmostEqual(circle['rho'], t/abs(p-1), delta=1e-6, msg=f'{vid} at {t}')
                horizons.append(circle['rho'])
            # The frames stand a steady step of t apart: steps that lengthen are an expansion that speeds up.
            steps = [later-earlier for earlier, later in zip(radii, radii[1:])]
            changes = [later-earlier for earlier, later in zip(steps, steps[1:])]
            self.assertTrue(all(step > 0 for step in steps), vid)
            if vid == 'inflation':
                self.assertTrue(all(change > 0 for change in changes), vid)
            elif vid == 'tracking':
                self.assertTrue(all(change < 0 for change in changes), vid)
            else:
                self.assertLess(max(abs(change) for change in changes), 4e-6, vid)
            legend = {cls: text for _, cls, text in view['figure']['legend']}
            if horizon is None:
                self.assertNotIn('horizon', legend, vid)
                continue
            self.assertIn(f'the {horizon} horizon', legend['horizon'], vid)
            # The ring meets its horizon at t = t_0, from inside the event horizon and from outside the particle horizon.
            k = values.index(1)
            self.assertAlmostEqual(radii[k], horizons[k], delta=2e-6, msg=vid)
            inside = vid == 'inflation'
            self.assertTrue(all((r < h) == inside for r, h in zip(radii[:k], horizons[:k])), vid)
            self.assertTrue(all((r > h) == inside for r, h in zip(radii[k+1:], horizons[k+1:])), vid)


if __name__ == '__main__':
    unittest.main()
