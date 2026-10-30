"""Kerr radiating and Kerr in a universe, held to the stored charts.

The radiating metric comes from Vaidya and Patel's own restatement of 1989, so its claims are
pinned here on the stored mathematics: pure radiation along one null direction, Kerr's metric when
the radiation parameter vanishes, flat space when the mass does whatever the rotation parameter
does, and the fall of the mass and the angular momentum. Vaidya's Kerr metric in Einstein's
static universe is pinned to his perfect fluid, its two limits and its horizon, and his Kerr
metric in an expanding universe to its two equal pressures and its horizon condition.
Run only this module with the collection's Python environment.
"""
import importlib.util
import itertools
import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / 'MFS' / 'assets' / 'data'
DERIVATIONS = ROOT / '_tools' / 'derivations'
HAS_SYMPY = importlib.util.find_spec('sympy') is not None
sys.path.insert(0, str(DERIVATIONS))
NAME = 'vaidya_patel_kerr'
RADIATING = ('vaidya_patel', 'retarded')
EINSTEIN = ('einstein_vaidya', 'einstein_boyer_lindquist')
EXPANDING = ('expanding_vaidya',)


def data(kind, name=NAME):
    return json.loads((DATA / kind / (name + '.json')).read_text())


class Stored:
    """One stored chart read through the checker's Reader, its held names written out."""

    def __init__(self, system):
        import sympy as sp
        import verify_metrics as vm
        self.sp = sp
        self.entry = next(c for c in data('metrics')['coordinates'] if c['id'] == system)
        key = (NAME, system)
        self.reader = vm.Reader(self.entry['coords'], [p['symbol'] for p in self.entry['parameters']], (),
                                held=vm.HELD.get(key, ()), rates=vm.RATES.get(key))
        self.held = key in vm.RATES
        self.x = [self.reader.symbol[c] for c in self.entry['coords']]
        self.P = self.reader.parameters

    def value(self, text):
        sp = self.sp
        v = self.reader(text).subs(self.reader.c, 1)
        if self.held:
            v = self.reader.by_rates(v.doit())
        for _ in range(3):
            v = v.subs(self.reader.held)
        return v.doit()

    def matrix(self, items):
        out = self.sp.zeros(4)
        for item in items:
            i, j = (self.entry['coords'].index(c) for c in item['indices'])
            out[i, j] = self.value(item['value'])
        return out

    def g(self):
        return self.matrix(self.entry['metric_components'])

    def ginv(self):
        return self.matrix(self.entry['inverse_metric_components'])

    def ricci(self):
        return self.matrix(self.entry['ricci_tensor']['variants']['ll']['nonzero'])

    def einstein(self):
        return self.matrix(self.entry['einstein_tensor']['variants']['ll']['nonzero'])

    def scalar(self, field):
        return self.value(self.entry[field].split('=', 1)[1])

    def riemann(self):
        out = {}
        for item in self.entry['riemann']['variants']['llll']['nonzero']:
            index = tuple(self.entry['coords'].index(c) for c in item['indices'])
            out[index] = self.value(item['value'])
        return out


def number(sp, value, point):
    return complex(sp.N(sp.sympify(value).subs(point), 30))


@unittest.skipUnless(HAS_SYMPY, 'sympy is not installed')
class Radiating(unittest.TestCase):
    """The claims that rest on Vaidya and Patel's restatement rather than on the text of 1973."""

    @classmethod
    def setUpClass(cls):
        import sympy as sp
        cls.sp = sp
        cls.charts = {s: Stored(s) for s in RADIATING}

    def parameters(self, chart, a='3/10', k='7/10', C='1/2'):
        sp = self.sp
        return {chart.P['a']: sp.Rational(a), chart.P['k']: sp.Rational(k), chart.P['C']: sp.Rational(C)}

    def points(self, chart):
        sp = self.sp
        if chart.entry['id'] == 'vaidya_patel':
            # (u, alpha, beta, T), each at r = T - u - (a/2) sin^2 alpha > 0
            rows = ((0, '1/5', 0, '3/5'), (1, '3/2', '2/5', 3), (3, '9/10', 1, 9))
        else:
            rows = ((0, '3/5', '1/5', 0), (1, '37/20', '3/2', '2/5'), (3, '59/10', '9/10', 1))
        return [dict(zip(chart.x, (sp.Rational(v) for v in row))) for row in rows]

    def null_form(self, chart, point, numbers):
        sp = self.sp
        A = sp.sqrt(numbers[chart.P['k']] ** 2 + numbers[chart.P['a']] * point[chart.x[0]])
        if chart.entry['id'] == 'vaidya_patel':
            alpha = point[chart.x[1]]
            return sp.Matrix([1, 0, A * sp.sin(alpha) ** 2, 0])
        alpha = point[chart.x[2]]
        return sp.Matrix([1, 0, 0, A * sp.sin(alpha) ** 2])

    def test_the_ricci_tensor_is_pure_radiation_along_one_null_direction(self):
        sp = self.sp
        for system, chart in self.charts.items():
            ricci, ginv = chart.ricci(), chart.ginv()
            numbers = self.parameters(chart)
            for point in self.points(chart):
                at = {**point, **numbers}
                theta = self.null_form(chart, point, numbers)
                Gi = ginv.subs(at)
                self.assertAlmostEqual(number(sp, (theta.T * Gi * theta)[0], {}), 0, places=20, msg=system)
                R = ricci.subs(at)
                sigma = sp.N(R[0, 0], 30)
                self.assertGreater(sigma, 0, system)
                for i, j in itertools.product(range(4), repeat=2):
                    residual = abs(sp.N(R[i, j] - sigma * theta[i] * theta[j], 30))
                    self.assertLess(residual, 1e-20, f'{system} R[{i}{j}]')

    def test_with_no_mass_the_spacetime_is_flat_however_the_rotation_parameter_runs(self):
        sp = self.sp
        for system, chart in self.charts.items():
            riemann = chart.riemann()
            for a in ('3/10', '-1/4', 2):
                numbers = self.parameters(chart, a=a, C=0)
                for point in self.points(chart)[:2]:
                    if numbers[chart.P['k']] ** 2 + numbers[chart.P['a']] * point[chart.x[0]] <= 0:
                        continue
                    for index, value in riemann.items():
                        self.assertAlmostEqual(number(sp, value, {**point, **numbers}), 0, places=20,
                                               msg=f'{system} a={a} {index}')

    def test_with_no_radiation_parameter_it_is_kerrs_metric(self):
        """At a = 0 the Ricci tensor vanishes and the Kretschmann scalar is Kerr's, with spin k and
        mass C/k^3, at Kerr's r and polar angle."""
        sp = self.sp
        chart = self.charts['retarded']
        u, r, alpha, beta = chart.x
        numbers = self.parameters(chart, a=0)
        k = numbers[chart.P['k']]
        m = numbers[chart.P['C']] / k ** 3
        ricci = chart.ricci()
        K = chart.scalar('kretschmann')
        S = r ** 2 + k ** 2 * sp.cos(alpha) ** 2
        kerr = 48 * m ** 2 * (r ** 2 - k ** 2 * sp.cos(alpha) ** 2) * (S ** 2 - 16 * r ** 2 * k ** 2 * sp.cos(alpha) ** 2) / S ** 6
        for point in self.points(chart):
            at = {**point, **numbers}
            for value in ricci:
                self.assertAlmostEqual(number(sp, value, at), 0, places=20)
            self.assertAlmostEqual(number(sp, K - kerr, at), 0, places=18)
        # In the barred chart too, at a = 0, where q = u - T is minus Kerr's r.
        bar = self.charts['vaidya_patel']
        for value in bar.ricci():
            for point in self.points(bar):
                self.assertAlmostEqual(number(sp, value, {**point, **self.parameters(bar, a=0)}), 0, places=20)

    def test_the_mass_and_the_angular_momentum_fall_together_as_the_rotation_parameter_grows(self):
        sp = self.sp
        chart = self.charts['retarded']
        u = chart.x[0]
        a, k, C = (sp.Symbol(n, positive=True) for n in ('a', 'k', 'C'))
        A = chart.value('A').subs({chart.P['a']: a, chart.P['k']: k})
        m = chart.value('m').subs({chart.P['a']: a, chart.P['k']: k, chart.P['C']: C})
        self.assertEqual(sp.simplify(A ** 2 - (k ** 2 + a * u)), 0)
        self.assertEqual(sp.simplify(m * A ** 3 - C), 0)
        for value in (sp.diff(m, u), sp.diff(m * A, u)):
            for point in ({u: 0, a: sp.Rational(3, 10), k: 1, C: 1}, {u: 5, a: 2, k: sp.Rational(1, 3), C: 7}):
                self.assertLess(number(sp, value, point).real, 0)
        # The slopes the checker declares are those of the definitions.
        self.assertEqual(sp.simplify(sp.diff(A, u) - a / (2 * A)), 0)
        self.assertEqual(sp.simplify(sp.diff(m, u) + 3 * a * m / (2 * A ** 2)), 0)


@unittest.skipUnless(HAS_SYMPY, 'sympy is not installed')
class EinsteinUniverse(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        import sympy as sp
        cls.sp = sp
        cls.charts = {s: Stored(s) for s in EINSTEIN}

    def at(self, chart, triple, m='1/2', a='2/5', R=3):
        """(u, alpha, beta, t) of the Vaidya chart, read into each chart with r = t - u."""
        sp = self.sp
        u, alpha, beta, t = (sp.Rational(str(v)) for v in triple)
        values = (u, alpha, beta, t) if chart.entry['id'] == 'einstein_vaidya' else (t, t - u, alpha, beta)
        return {**dict(zip(chart.x, values)),
                chart.P['m']: sp.Rational(m), chart.P['a']: sp.Rational(a), chart.P['R']: sp.Rational(R)}

    def test_each_chart_holds_vaidyas_perfect_fluid(self):
        """Vaidya's p. 515: 8 pi p = -(1 - 2 m mu)/R^2 + Lambda and 8 pi rho = 3(1 - 2 m mu)/R^2 - Lambda,
        so in this signature G^i_k has -(1 - 2 m mu)/R^2 three times and -3(1 - 2 m mu)/R^2 once, with
        the single eigenvalue's eigenvector timelike."""
        import mpmath
        sp = self.sp
        for system, chart in self.charts.items():
            einstein, ginv, g = chart.einstein(), chart.ginv(), chart.g()
            for triple in ((0.3, 0.9, 0.2, 2.1), (-1.0, 1.7, 1.0, 1.4), (0.5, 1.2, 2.0, 4.5)):
                at = self.at(chart, triple)
                r = sp.Rational(str(triple[3])) - sp.Rational(str(triple[0]))
                alpha = sp.Rational(str(triple[1]))
                m, a, R = sp.Rational(1, 2), sp.Rational(2, 5), 3
                Sigma = (R ** 2 - a ** 2) * sp.sin(r / R) ** 2 + a ** 2 * sp.cos(alpha) ** 2
                shift = float(sp.N((1 - 2 * m * R * sp.sin(r / R) * sp.cos(r / R) / Sigma) / R ** 2, 30))
                mixed = mpmath.matrix(sp.N(ginv.subs(at) * einstein.subs(at), 30).tolist())
                values, vectors = mpmath.eig(mixed)
                values = [complex(v).real for v in values]
                single = min(range(4), key=lambda i: abs(values[i] + 3 * shift))
                self.assertAlmostEqual(values[single], -3 * shift, places=12, msg=system)
                for i in range(4):
                    if i != single:
                        self.assertAlmostEqual(values[i], -shift, places=12, msg=system)
                v = vectors[:, single]
                G = mpmath.matrix(sp.N(g.subs(at), 30).tolist())
                self.assertLess(float(mpmath.re((v.T * G * v)[0])), 0, system)

    def test_with_no_mass_it_is_einsteins_static_universe(self):
        """At m = 0 the Ricci scalar is 6/R^2 and the Kretschmann scalar 12/R^4 everywhere, the
        values of Einstein's universe of radius R, a product of time with a three sphere."""
        sp = self.sp
        for system, chart in self.charts.items():
            R = chart.scalar('ricci_scalar')
            K = chart.scalar('kretschmann')
            for triple in ((0.3, 0.9, 0.2, 2.1), (-1.0, 1.7, 1.0, 1.4)):
                at = self.at(chart, triple, m=0)
                self.assertAlmostEqual(number(sp, R, at), 6 / 9, places=20, msg=system)
                self.assertAlmostEqual(number(sp, K, at), 12 / 81, places=20, msg=system)

    def test_as_the_universe_grows_without_bound_it_is_kerrs_black_hole(self):
        """mu tends to Kerr's r/(r^2 + a^2 cos^2 alpha), and the Kretschmann scalar to Kerr's."""
        sp = self.sp
        chart = self.charts['einstein_boyer_lindquist']
        t, r, alpha, beta = chart.x
        mu = chart.value('\\mu')
        a = chart.P['a']
        kerr_mu = r / (r ** 2 + a ** 2 * sp.cos(alpha) ** 2)
        self.assertEqual(sp.simplify(sp.limit(mu.subs({chart.P['R']: 1 / sp.Symbol('e', positive=True)}),
                                              sp.Symbol('e', positive=True), 0) - kerr_mu), 0)
        K = chart.scalar('kretschmann')
        m = chart.P['m']
        S = r ** 2 + a ** 2 * sp.cos(alpha) ** 2
        kerr = 48 * m ** 2 * (r ** 2 - a ** 2 * sp.cos(alpha) ** 2) * (S ** 2 - 16 * r ** 2 * a ** 2 * sp.cos(alpha) ** 2) / S ** 6
        point = {t: 0, r: sp.Rational(5, 2), alpha: sp.Rational(7, 10), beta: 0, m: sp.Rational(1, 2), a: sp.Rational(2, 5)}
        gaps = [abs(number(sp, K - kerr, {**point, chart.P['R']: R})) for R in (10 ** 3, 10 ** 5)]
        self.assertLess(gaps[1], gaps[0] / 1000)
        self.assertLess(gaps[1], 1e-8)

    def test_the_horizons_are_vaidyas_and_curvature_pulls_them_in_by_a_third_of_r_squared_over_R_squared(self):
        sp = self.sp
        chart = self.charts['einstein_boyer_lindquist']
        r = chart.x[1]
        m, a, R = chart.P['m'], chart.P['a'], chart.P['R']
        delta = chart.value('\\Delta')
        # Delta/cos^2(r/R) is Vaidya's R^2 tan^2(r/R) - 2mR tan(r/R) + a^2, his (18), at any values.
        vaidya = R ** 2 * sp.tan(r / R) ** 2 - 2 * m * R * sp.tan(r / R) + a ** 2
        for values in ({r: 1, m: 1, a: sp.Rational(3, 5), R: 4}, {r: sp.Rational(7, 3), m: 2, a: 1, R: 9}):
            self.assertLess(abs(sp.N((delta / sp.cos(r / R) ** 2 - vaidya).subs(values), 30)), 1e-25)
        # g^rr = Delta/Sigma vanishes on R tan(r/R) = m + sqrt(m^2 - a^2).
        ginv = chart.ginv()
        values = {m: 1, a: sp.Rational(3, 5), R: 40}
        horizon = 40 * sp.atan((1 + sp.sqrt(1 - sp.Rational(9, 25))) / 40)
        self.assertAlmostEqual(number(sp, ginv[1, 1], {**values, r: horizon, chart.x[2]: sp.Rational(1, 3)}), 0, places=20)
        # The fraction by which the horizon falls inside Kerr's tends to r^2/3R^2, Vaidya's p. 516,
        # with a correction of order r^2/R^2 to that.
        kerr = 1 + sp.sqrt(1 - sp.Rational(9, 25))
        ratios = []
        for big in (40, 400):
            inside = big * sp.atan(kerr / big)
            ratios.append(float(sp.N((kerr - inside) / kerr / (kerr ** 2 / (3 * big ** 2)), 30)))
        self.assertLess(abs(ratios[0] - 1), 1e-2)
        self.assertLess(abs(ratios[1] - 1), abs(ratios[0] - 1) / 50)


@unittest.skipUnless(HAS_SYMPY, 'sympy is not installed')
class ExpandingUniverse(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        import sympy as sp
        cls.sp = sp
        cls.charts = {s: Stored(s) for s in EXPANDING}

    def test_the_fluid_has_two_equal_pressures_and_a_third(self):
        import mpmath
        sp = self.sp
        for system, chart in self.charts.items():
            einstein, ginv, g = chart.einstein(), chart.ginv(), chart.g()
            t = chart.x[3] if system == 'expanding_vaidya' else chart.x[0]
            scale = {chart.P['F']: sp.Rational(3, 10) * t + sp.Rational(1, 10) * t ** 2}
            for triple in ((0.3, 0.9, 0.2, 2.1), (-1.0, 1.7, 1.0, 1.4)):
                u, alpha, beta, tv = (sp.Rational(str(v)) for v in triple)
                values = (u, alpha, beta, tv) if system == 'expanding_vaidya' else (tv, tv - u, alpha, beta)
                at = {**dict(zip(chart.x, values)), chart.P['m']: sp.Rational(1, 2),
                      chart.P['a']: sp.Rational(2, 5), chart.P['R']: 3}
                mixed = (ginv * einstein).subs(scale).doit().subs(at)
                ev, vectors = mpmath.eig(mpmath.matrix(sp.N(mixed, 30).tolist()))
                ev = sorted(complex(v).real for v in ev)
                equal = [abs(ev[i] - ev[i + 1]) < 1e-12 for i in range(3)]
                self.assertEqual(sum(equal), 1, f'{system} {ev}')

    def test_the_horizon_needs_a_squared_over_m_squared_below_e_to_the_minus_four_F(self):
        """Vaidya's (19), R^2 tan^2 - 2 m e^{-2F} R tan + a^2 = 0, has a real root exactly while
        a^2/m^2 < e^{-4F}: its discriminant is 4R^2(m^2 e^{-4F} - a^2)."""
        sp = self.sp
        m, a, R, F, x = sp.symbols('m a R F x', positive=True)
        quadratic = sp.Poly(R ** 2 * x ** 2 - 2 * m * sp.exp(-2 * F) * R * x + a ** 2, x)
        self.assertEqual(sp.simplify(quadratic.discriminant() - 4 * R ** 2 * (m ** 2 * sp.exp(-4 * F) - a ** 2)), 0)


@unittest.skipUnless(HAS_SYMPY, 'sympy is not installed')
class ChartsAgree(unittest.TestCase):
    """Each chart after the first of its spacetime is the first pulled back, on the stored metrics."""

    def test_pullbacks(self):
        import sympy as sp
        pairs = {'retarded': 'vaidya_patel', 'einstein_boyer_lindquist': 'einstein_vaidya'}
        for system, first in pairs.items():
            here, there = Stored(system), Stored(first)
            X = here.x
            if system == 'retarded':
                u, r, alpha, beta = X
                images = [u, alpha, beta, r + u + here.P['a'] * sp.sin(alpha) ** 2 / 2]
            else:
                # Vaidya's chart holds u and t through r = t - u alone, so (t - r, alpha, phi, t) is
                # an image of the Boyer-Lindquist event, with dt -> dt + f dr and d(beta) -> d(phi) + h dr.
                t, r, alpha, phi = X
                images = [t - r, alpha, phi, t]
            J = sp.Matrix(4, 4, lambda i, j: sp.diff(images[i], X[j]))
            if system == 'einstein_boyer_lindquist':
                m, a, R = here.P['m'], here.P['a'], here.P['R']
                delta = here.value('\\Delta')
                f = -2 * m * R * sp.sin(r / R) * sp.cos(r / R) / delta
                J[0, 1], J[3, 1], J[2, 1] = f - 1, f, a / delta
            # The Reader names a coordinate or a parameter by one symbol in every chart, and F(t) is
            # the same function of the same t in both expanding charts, so only the coordinates move.
            old = there.g().subs(dict(zip(there.x, images)), simultaneous=True)
            pulled = J.T * old * J
            numbers = ({here.P['a']: sp.Rational(3, 10), here.P['k']: sp.Rational(7, 10), here.P['C']: sp.Rational(1, 2)}
                       if system == 'retarded' else
                       {here.P['m']: sp.Rational(1, 2), here.P['a']: sp.Rational(2, 5), here.P['R']: 3})
            point = dict(zip(X, (sp.Rational(1, 2), sp.Rational(13, 10), sp.Rational(7, 10), sp.Rational(1, 5))))
            if system == 'retarded':
                point = dict(zip(X, (sp.Rational(1, 2), sp.Rational(13, 10), sp.Rational(7, 10), sp.Rational(1, 5))))
            gap = pulled - here.g()
            for value in gap:
                self.assertAlmostEqual(number(sp, value, {**point, **numbers}), 0, places=20, msg=system)


if __name__ == '__main__':
    unittest.main()
