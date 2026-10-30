"""Bekenstein's universes with no singularity, held to their stored charts and drawings.

Pin the field equations of dust coupled to a conformally coupled scalar field and radiation in
each chart, the amplitude relations at concrete constants, the condition for no singularity with
the free wave of the field in it and its reduction to Novello and Perez Bergliaffa's summary,
the sign of the energy density at the bounce, and every limit the page states.
Run only this module with the collection's Python environment.
"""
import importlib.util
import json
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / 'MFS' / 'assets' / 'data'
DERIVATIONS = ROOT / '_tools' / 'derivations'
HAS_SYMPY = importlib.util.find_spec('sympy') is not None
sys.path.insert(0, str(DERIVATIONS))
NAME = 'bekenstein_universes'
SYSTEMS = ('closed', 'flat', 'open')
CURVATURE = {'closed': 1, 'flat': 0, 'open': -1}


def data(kind, name=NAME):
    return json.loads((DATA / kind / (name + '.json')).read_text())


def closed_amplitude2(ell, f, N, B, P, Q):
    """A^2 of the closed universe, from the Friedmann constraint, G = c = 1."""
    return ell ** 2 + 4 * math.pi / 3 * (2 * B + P ** 2 + Q ** 2 - f ** 2 * N ** 2)


def open_amplitude2(ell, f, N, B, Dp, Dm):
    """A^2 of the open universe with the bounce at eta = 0, a = A cosh(eta) - l."""
    return ell ** 2 - 4 * math.pi / 3 * (2 * B + f ** 2 * N ** 2 - 4 * Dp * Dm)


def flat_bounce(N, mu, f, B, d, e):
    """a_b of the flat universe: 2 N mu a_b = -(2B + d^2 + 2 N f e)."""
    return -(2 * B + d ** 2 + 2 * N * f * e) / (2 * N * mu)


@unittest.skipUnless(HAS_SYMPY, 'sympy is not installed')
class FieldEquations(unittest.TestCase):
    """Each published chart, with constants chosen in its range of no singularity, satisfies
    G = 8 pi (T_field + T_dust + T_radiation) and the field equation at points, the field's stress
    tensor being the conformal one and the dust's mass mu + f psi."""

    @classmethod
    def setUpClass(cls):
        import sympy as sp
        import verify_metrics as vm
        cls.sp = sp
        cls.charts = {}
        for system in SYSTEMS:
            entry = next(c for c in data('metrics')['coordinates'] if c['id'] == system)
            reader = vm.Reader(entry['coords'], [p['symbol'] for p in entry['parameters']], ())
            g = sp.zeros(4)
            for item in entry['metric_components']:
                i, j = (entry['coords'].index(c) for c in item['indices'])
                g[i, j] = reader(item['value'])
            G = sp.zeros(4)
            for item in entry['einstein_tensor']['variants']['ll']['nonzero']:
                i, j = (entry['coords'].index(c) for c in item['indices'])
                G[i, j] = reader(item['value'])
            X = [reader.symbol[c] for c in entry['coords']]
            cls.charts[system] = (reader, X, g, G, entry)

    def sources(self, system, f, B, P, Q):
        """The Einstein tensor minus 8 pi T, the field equation and the trace identity, with
        N = l and mu = 3/4 pi, so that 4 pi N mu/3 = l."""
        sp = self.sp
        reader, X, g, G, _ = self.charts[system]
        eta = X[0]
        ell = reader.parameters['ell']
        N, mu = ell, sp.Rational(3, 4) / sp.pi
        F = {'closed': -f * N + P * sp.cos(eta) + Q * sp.sin(eta),
             'flat': -f * N * eta ** 2 / 2 + P * eta + Q,
             'open': f * N + P * sp.exp(eta) + Q * sp.exp(-eta)}[system]
        a = sp.sqrt(g[1, 1])
        psi = F / a
        ginv = sp.diag(*[1 / g[i, i] for i in range(4)])
        gamma = [[[sum(ginv[i, l] * (sp.diff(g[l, j], X[k]) + sp.diff(g[l, k], X[j]) - sp.diff(g[j, k], X[l]))
                       for l in range(4)) / 2 for k in range(4)] for j in range(4)] for i in range(4)]

        def hessian(s):
            d = [sp.diff(s, x) for x in X]
            return sp.Matrix(4, 4, lambda i, j: sp.diff(d[j], X[i]) - sum(gamma[l][i][j] * d[l] for l in range(4)))
        dpsi = [sp.diff(psi, x) for x in X]
        grad2 = sum(ginv[i, i] * dpsi[i] ** 2 for i in range(4))
        H2, H1 = hessian(psi ** 2), hessian(psi)
        box2 = sum(ginv[i, i] * H2[i, i] for i in range(4))
        box = sum(ginv[i, i] * H1[i, i] for i in range(4))
        R = sum(ginv[i, j] * (G[i, j]) for i in range(4) for j in range(4)) * -1
        n = N / a ** 3
        u = [-a, 0, 0, 0]
        rad = B / a ** 4
        T = sp.Matrix(4, 4, lambda i, j: dpsi[i] * dpsi[j] - g[i, j] * grad2 / 2
                      + (g[i, j] * box2 - H2[i, j] + G[i, j] * psi ** 2) / 6
                      + n * (mu + f * psi) * u[i] * u[j] + 4 * rad * u[i] * u[j] / 3 + rad * g[i, j] / 3)
        return G - 8 * sp.pi * T, (box - R * psi / 6 - f * n) * a ** 3, R - 8 * sp.pi * n * mu

    def check(self, system, f, B, P, Q, amplitude):
        sp = self.sp
        reader, X, *_ = self.charts[system]
        other = reader.parameters['a_b' if system == 'flat' else 'A']
        einstein, wave, trace = self.sources(system, f, B, P, Q)
        for time, angle in ((sp.Rational(3, 10), sp.Rational(7, 10)), (sp.Rational(-11, 10), sp.Rational(13, 10))):
            point = {reader.parameters['ell']: 1, other: amplitude, X[0]: time, X[1]: angle,
                     X[2]: sp.Rational(9, 10), X[3]: sp.Rational(1, 5)}
            gap = einstein.subs(point).evalf(30)
            self.assertLess(max(abs(v) for v in gap), 1e-20, f'{system}: G = 8 pi T at eta = {time}')
            self.assertLess(abs(wave.subs(point).evalf(30)), 1e-20, f'{system}: the field equation')
            self.assertLess(abs(trace.subs(point).evalf(30)), 1e-20, f'{system}: R = 8 pi N mu/a^3')

    def test_closed(self):
        sp = self.sp
        f, B, P, Q = sp.Rational(1, 2), sp.Rational(1, 20), sp.Rational(1, 5), sp.Rational(1, 10)
        A = sp.sqrt(1 + sp.Rational(4, 3) * sp.pi * (2 * B + P ** 2 + Q ** 2 - f ** 2))
        self.assertLess(float(A), 1)
        self.check('closed', f, B, P, Q, A)

    def test_flat(self):
        sp = self.sp
        f, B, d, e = sp.Rational(-9, 10), sp.Rational(1, 20), sp.Rational(1, 5), sp.Rational(1, 10)
        a_b = -(2 * B + d ** 2 + 2 * f * e) / (2 * sp.Rational(3, 4) / sp.pi)
        self.assertGreater(float(a_b), 0)
        self.check('flat', f, B, d, e, a_b)

    def test_open(self):
        sp = self.sp
        f, B, Dp, Dm = sp.Rational(2, 5), sp.Rational(1, 20), sp.Rational(1, 2), sp.Rational(2, 5)
        A = sp.sqrt(1 - sp.Rational(4, 3) * sp.pi * (2 * B + f ** 2 - 4 * Dp * Dm))
        self.assertGreater(float(A), 1)
        self.check('open', f, B, Dp, Dm, A)

    def test_a_wrong_amplitude_breaks_the_equations(self):
        sp = self.sp
        f, B, P, Q = sp.Rational(1, 2), sp.Rational(1, 20), sp.Rational(1, 5), sp.Rational(1, 10)
        A = sp.sqrt(1 + sp.Rational(4, 3) * sp.pi * (2 * B + f ** 2))
        with self.assertRaises(AssertionError):
            self.check('closed', f, B, P, Q, A)

    def test_the_ricci_scalar_is_six_l_over_a_cubed(self):
        sp = self.sp
        for system in SYSTEMS:
            reader, X, g, _, entry = self.charts[system]
            a = sp.sqrt(g[1, 1])
            R = reader(entry['ricci_scalar'].partition('=')[2])
            other = reader.parameters['a_b' if system == 'flat' else 'A']
            amplitude = sp.Integer(3) if system == 'open' else sp.Rational(1, 2)
            point = {reader.parameters['ell']: sp.Rational(7, 5), other: amplitude, X[0]: sp.Rational(2, 3),
                     X[1]: sp.Rational(1, 3), X[2]: 1}
            self.assertLess(abs((R - 6 * reader.parameters['ell'] / a ** 3).subs(point).evalf(30)), 1e-25, system)


class Conditions(unittest.TestCase):
    """The amplitude relations and the conditions for no singularity at concrete constants, the
    free wave of the field included, and their limits."""

    def test_the_closed_condition_holds_the_free_wave(self):
        ell, N = 1.0, 1.0
        # f^2 N^2 = 0.25 beats 2B = 0.1 alone but not 2B + P^2 + Q^2 = 0.35: the radius reaches zero.
        f, B, P, Q = 0.5, 0.05, 0.4, 0.3
        self.assertGreater(f ** 2 * N ** 2, 2 * B)
        self.assertLess(f ** 2 * N ** 2, 2 * B + P ** 2 + Q ** 2)
        self.assertGreater(math.sqrt(closed_amplitude2(ell, f, N, B, P, Q)), ell)
        # A weaker wave lets the interaction win, and the radius stays above zero.
        P, Q = 0.2, 0.1
        self.assertGreater(f ** 2 * N ** 2, 2 * B + P ** 2 + Q ** 2)
        self.assertLess(math.sqrt(closed_amplitude2(ell, f, N, B, P, Q)), ell)

    def test_the_closed_condition_is_exactly_the_sign_of_the_amplitude(self):
        for f, B, P, Q in ((0.3, 0.02, 0.1, 0.05), (0.55, 0.1, 0.3, 0.2), (0.62, 0.05, 0.0, 0.5), (0.45, 0.0, 0.2, 0.3)):
            A2 = closed_amplitude2(1.0, f, 1.0, B, P, Q)
            self.assertEqual(A2 < 1.0, f ** 2 > 2 * B + P ** 2 + Q ** 2, (f, B, P, Q))

    def test_with_no_free_wave_the_condition_is_the_reviews(self):
        for f, B in ((0.5, 0.1), (0.3, 0.05), (0.2, 0.05), (0.4, 0.08)):
            A2 = closed_amplitude2(1.0, f, 1.0, B, 0.0, 0.0)
            self.assertEqual(A2 < 1.0, f ** 2 > 2 * B, (f, B))

    def test_with_no_coupling_radiation_or_wave_the_closed_universe_is_friedmanns_dust(self):
        ell = 1.3
        A = math.sqrt(closed_amplitude2(ell, 0.0, 1.0, 0.0, 0.0, 0.0))
        self.assertAlmostEqual(A, ell, places=14)
        for eta in (0.4, 1.7, 3.0):
            self.assertAlmostEqual(ell - A * math.cos(eta), ell * (1 - math.cos(eta)), places=14)
        self.assertAlmostEqual(ell - A * math.cos(0.0), 0.0, places=14)

    def test_a_static_universe_at_zero_amplitude(self):
        ell = 1.0
        # A = 0 when f^2 N^2 - 2B - P^2 - Q^2 = 3 l^2/4 pi.
        f = math.sqrt(3 / (4 * math.pi) + 0.1 + 0.04)
        self.assertAlmostEqual(closed_amplitude2(ell, f, 1.0, 0.05, 0.2, 0.0), 0.0, places=14)

    def test_the_flat_condition(self):
        N, mu = 1.0, 3 / (4 * math.pi)
        self.assertGreater(flat_bounce(N, mu, -0.9, 0.05, 0.2, 0.1), 0)
        # -2Nfe = 0.18 against 2B + d^2 = 0.1 + 0.25: the bounce radius is negative.
        self.assertLess(flat_bounce(N, mu, -0.9, 0.05, 0.5, 0.1), 0)
        self.assertEqual(flat_bounce(N, mu, 0.0, 0.0, 0.0, 0.0), 0.0)

    def test_the_open_condition(self):
        ell = 1.0
        self.assertGreater(math.sqrt(open_amplitude2(ell, 0.4, 1.0, 0.05, 0.5, 0.4)), ell)
        # 4 D+ D- = 0.08 against 2B + f^2 N^2 = 0.26: the radius reaches zero.
        self.assertLess(math.sqrt(open_amplitude2(ell, 0.4, 1.0, 0.05, 0.2, 0.1)), ell)
        self.assertAlmostEqual(math.sqrt(open_amplitude2(ell, 0.0, 1.0, 0.0, 0.0, 0.0)), ell, places=14)

    def test_the_energy_density_at_the_bounce(self):
        # 3(a'^2 + k a^2)/8 pi a^4 with a' = 0: positive closed, zero flat, negative open.
        a = 0.5
        signs = {k: math.copysign(1, 3 * k * a * a / (8 * math.pi * a ** 4)) if k else 0 for k in (1, 0, -1)}
        self.assertEqual(signs, {1: 1.0, 0: 0, -1: -1.0})


class Drawings(unittest.TestCase):
    """The drawings hold what their captions say."""

    def test_the_movies_are_spheres_and_hyperboloids_of_radius_a(self):
        views = {v['id']: v for v in data('embedding')['views']}
        shapes = {'closed': (lambda e: 1 - 0.5 * math.cos(e), math.sin, lambda x, a: -a * math.cos(x)),
                  'open': (lambda e: 1.5 * math.cosh(e) - 1, math.sinh, lambda x, a: a * (math.cosh(x) - 1))}
        for vid, (scale, around, height) in shapes.items():
            frames = views[vid]['movie']['frames']
            self.assertGreater(len(frames), 30, vid)
            labels = [f['label'] for f in frames]
            self.assertEqual(len(set(labels)), len(labels), vid)
            steps = [b['value'] - a['value'] for a, b in zip(frames, frames[1:])]
            self.assertLess(max(steps) - min(steps), 1e-5, f'{vid}: frames evenly spaced in eta')
            for frame in frames:
                a = scale(frame['value'])
                for x, rho, z in frame['pieces'][0]['points']:
                    self.assertAlmostEqual(rho, a * around(x), places=5)
                    self.assertAlmostEqual(z, height(x, a), places=5)
            bounce = min(frames, key=lambda f: abs(f['value']))
            self.assertEqual(bounce['value'], 0.0)
            self.assertAlmostEqual(scale(0.0), 0.5, places=14)

    def test_the_conformal_diagram_has_each_universe(self):
        self.assertEqual([v['id'] for v in data('conformal')['views']], list(SYSTEMS))

    def test_the_curvature_figure_draws_each_universe(self):
        views = {v['id']: v for v in data('curvature')['views']}
        self.assertEqual(set(views), set(SYSTEMS))
        for system, view in views.items():
            self.assertEqual({c['panel'] for c in view['curves']}, {'K', 'R'}, system)


if __name__ == '__main__':
    unittest.main()
