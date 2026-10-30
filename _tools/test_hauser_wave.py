"""Hauser's twisting type N wave, held to its stored chart and its declared drawings.

Pin the vacuum on Hauser's equation at points, Novoseller's solution of that equation, the twist
function Sigma = 2P^2 of the metric functions, the type N condition on the ray direction, the
homothety, the vanishing of both stored curvature scalars on the equation, and the declarations
that it has no conformal, embedding or curvature figure.
Run only this module with the collection's Python environment.
"""
import importlib.util
import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / 'MFS' / 'assets' / 'data'
DERIVATIONS = ROOT / '_tools' / 'derivations'
HAS_SYMPY = importlib.util.find_spec('sympy') is not None
sys.path.insert(0, str(DERIVATIONS))
sys.path.insert(0, str(ROOT / '_tools'))

# Points of rational coordinates, each with rational f and f' there; f'' and higher come from the equation.
POINTS = (({'t': (3, 10), 'r': (1, 2), 'x': (7, 10), 'y': (1, 5)}, (1, 1), (-1, 4)),
          ({'t': (-2, 1), 'r': (-3, 2), 'x': (6, 5), 'y': (0, 1)}, (4, 5), (1, 3)),
          ({'t': (5, 1), 'r': (1, 10), 'x': (2, 5), 'y': (1, 1)}, (3, 2), (-2, 3)))


def metric():
    return json.loads((DATA / 'metrics' / 'hauser_wave.json').read_text())


@unittest.skipUnless(HAS_SYMPY, 'sympy is not installed')
class StoredChart(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import sympy as sp
        import verify_metrics as vm
        cls.sp = sp
        cls.entry = metric()['coordinates'][0]
        cls.reader = vm.Reader(cls.entry['coords'], [p['symbol'] for p in cls.entry['parameters']], ())
        cls.symbols = [cls.reader.symbol[c] for c in cls.entry['coords']]
        cls.f = cls.reader.parameters['f']

    def read(self, text):
        return self.reader.surface(self.reader(text)).doit()

    def matrix(self, items, rank=2):
        sp = self.sp
        out = sp.MutableDenseNDimArray.zeros(*([4] * rank))
        for item in items:
            index = tuple(self.entry['coords'].index(c) for c in item['indices'])
            out[index] = self.read(item['value'])
        return out

    def on_shell(self, value, at, f0, f1, coefficient=3):
        """value at the point with f'' to f^(5) from 16(1 + t^2)f'' + coefficient f = 0."""
        sp = self.sp
        t = self.symbols[0]
        F0, F1, T = sp.symbols('F0 F1 T')
        second = -coefficient * F0 / (16 * (1 + T ** 2))
        derivatives = [F0, F1, second]
        for _ in range(3):
            last = derivatives[-1]
            derivatives.append(sp.diff(last, T) + sp.diff(last, F0) * F1 + sp.diff(last, F1) * second)
        known = {F0: f0, F1: f1, T: at[t]}
        value = sp.sympify(value)
        for k in range(5, 0, -1):
            value = value.subs(sp.Derivative(self.f, (t, k)), derivatives[k].subs(known))
        return value.subs(self.f, f0).subs(at)

    def points(self):
        sp = self.sp
        for coords, f0, f1 in POINTS:
            at = {self.reader.symbol[k]: sp.Rational(*v) for k, v in coords.items()}
            at[self.reader.parameters['L']] = 1
            yield at, sp.Rational(*f0), sp.Rational(*f1)

    def test_the_stored_ricci_tensor_vanishes_exactly_on_hauser_equation(self):
        ricci = self.matrix(self.entry['ricci_tensor']['variants']['ll']['nonzero'])
        for at, f0, f1 in self.points():
            values = [self.on_shell(ricci[i, j], at, f0, f1) for i in range(4) for j in range(4)]
            self.assertTrue(all(v == 0 for v in values), at)
            wrong = [self.on_shell(ricci[i, j], at, f0, f1, coefficient=4) for i in range(4) for j in range(4)]
            self.assertTrue(any(v != 0 for v in wrong), at)

    def test_both_stored_scalars_vanish_on_hauser_equation(self):
        # The reason the spacetime has no curvature figure: K and R are zero on the equation.
        for field in ('ricci_scalar', 'kretschmann'):
            value = self.read(self.entry[field].partition('=')[2])
            for at, f0, f1 in self.points():
                self.assertEqual(self.on_shell(value, at, f0, f1), 0, (field, at))

    def test_the_rays_are_a_fourfold_principal_null_direction(self):
        # d_r is null, and C^a_{bcd} with d = r vanishes on the equation, while C itself does not.
        g = self.matrix(self.entry['metric_components'])
        weyl = self.matrix(self.entry['weyl_tensor']['variants']['ulll']['nonzero'], rank=4)
        r = 1
        self.assertEqual(g[r, r], 0)
        for at, f0, f1 in self.points():
            self.assertTrue(all(self.on_shell(weyl[a, b, c, r], at, f0, f1) == 0
                                for a in range(4) for b in range(4) for c in range(4)), at)
            self.assertTrue(any(self.on_shell(weyl[a, b, c, d], at, f0, f1) != 0
                                for a in range(4) for b in range(4) for c in range(4) for d in range(4)), at)

    def test_sigma_is_twice_p_squared(self):
        # 2i Sigma = P^2(dbar L - d Lbar) with d = d_zeta - L d_u, L = 2i(zeta + zetabar) and
        # P = (zeta + zetabar)^{3/2} f(u/(zeta + zetabar)^2), so Sigma = 2P^2, which makes the twist
        # Sigma/(r^2 + Sigma^2) the page's 2P^2/(r^2 + 4P^4).
        sp = self.sp
        u, x, y = sp.symbols('u x y', real=True)
        F = sp.Function('f', real=True)
        P = (2 * x) ** sp.Rational(3, 2) * F(u / (2 * x) ** 2)
        L, Lbar = 4 * sp.I * x, -4 * sp.I * x

        def d(h):
            return (sp.diff(h, x) - sp.I * sp.diff(h, y)) / 2 - L * sp.diff(h, u)

        def dbar(h):
            return (sp.diff(h, x) + sp.I * sp.diff(h, y)) / 2 - Lbar * sp.diff(h, u)
        sigma = P ** 2 * (dbar(L) - d(Lbar)) / (2 * sp.I)
        self.assertEqual(sp.expand(sigma - 2 * P ** 2), 0)
        # In the chart, the transverse factor 2(r^2 + Sigma^2)/P^2 is the published (r^2 + 256x^6f^4)/(4x^3f^2).
        g = self.matrix(self.entry['metric_components'])
        t, r, xx, yy = self.symbols
        f = self.f
        self.assertEqual(sp.expand(g[3, 3] - self.reader.parameters['L'] ** 2 * (-256 * xx ** 3 * f * (3 * f - 4 * t * sp.Derivative(f, t))
                                              - 128 * xx ** 2 * self.read('H')
                                              + (r ** 2 + 256 * xx ** 6 * f ** 4) / (4 * xx ** 3 * f ** 2))).doit(), 0)

    def test_x_y_and_r_scaled_together_multiply_the_metric_by_five(self):
        # x d_x + y d_y + 3r d_r is a homothety, L g = 5g, Hauser's homothetic parameter 5/2.
        sp = self.sp
        g = self.matrix(self.entry['metric_components'])
        xi = [0, 3 * self.symbols[1], self.symbols[2], self.symbols[3]]
        for at, f0, f1 in self.points():
            for i in range(4):
                for j in range(i, 4):
                    lie = (sum(xi[k] * sp.diff(g[i, j], self.symbols[k]) for k in range(4))
                           + sum(g[k, j] * sp.diff(xi[k], self.symbols[i]) + g[i, k] * sp.diff(xi[k], self.symbols[j])
                                 for k in range(4)))
                    self.assertEqual(self.on_shell(lie - 5 * g[i, j], at, f0, f1), 0, (i, j, at))

    def test_y_is_a_killing_translation(self):
        g = self.matrix(self.entry['metric_components'])
        self.assertTrue(all(not self.sp.sympify(g[i, j]).has(self.symbols[3]) for i in range(4) for j in range(4)))


class HauserFunction(unittest.TestCase):
    def test_novoseller_solutions_solve_hauser_equation(self):
        import mpmath as mp
        mp.mp.dps = 30
        even = lambda t: mp.hyp2f1(-mp.mpf(1) / 8, -mp.mpf(3) / 8, mp.mpf(1) / 2, -t ** 2)
        odd = lambda t: t * mp.hyp2f1(mp.mpf(1) / 8, mp.mpf(3) / 8, mp.mpf(3) / 2, -t ** 2)
        for f in (even, odd):
            for t in (mp.mpf('0.3'), mp.mpf('2.5'), mp.mpf('-7')):
                self.assertLess(abs(16 * (1 + t ** 2) * mp.diff(f, t, 2) + 3 * f(t)), mp.mpf('1e-20'))
        self.assertEqual(even(0), 1)
        self.assertAlmostEqual(float(mp.findroot(even, 5)), 5.46728770730835, places=12)

    def test_the_declared_function_of_the_drawings_is_the_even_solution(self):
        import mpmath as mp
        import numpy as np
        import null_rays as nr
        import sympy as sp
        t = sp.Symbol('t')
        f, df = nr.DECLARED_FUNCTIONS['hauser_f'], nr.DECLARED_FUNCTIONS['hauser_df']
        self.assertEqual(sp.diff(f(t), t), df(t))
        self.assertEqual(sp.simplify(16 * (1 + t ** 2) * sp.diff(f(t), t, 2) + 3 * f(t)), 0)
        even = lambda s: mp.hyp2f1(-mp.mpf(1) / 8, -mp.mpf(3) / 8, mp.mpf(1) / 2, -s ** 2)
        for s in (0.0, 0.3, 2.5, 5.0):
            self.assertAlmostEqual(float(sp.lambdify(t, f(t), 'numpy')(np.float64(s))), float(even(mp.mpf(s))), places=12)
            self.assertAlmostEqual(float(sp.lambdify(t, df(t), 'numpy')(np.float64(s))), float(mp.diff(even, mp.mpf(s))),
                                   places=12)


class Drawings(unittest.TestCase):
    def test_no_conformal_embedding_or_curvature_file_and_each_declared(self):
        import build_mfs_data as build
        import conformal
        import embedding
        self.assertIn('hauser_wave', conformal.NOT_DRAWN)
        self.assertIn('hauser_wave', embedding.NOT_DRAWN)
        self.assertIn('hauser_wave', build.WITHOUT_EMBEDDING)
        for kind in ('conformal', 'embedding', 'curvature'):
            self.assertFalse((DATA / kind / 'hauser_wave.json').exists(), kind)
        note = (DERIVATIONS / 'hauser_wave.md').read_text()
        for step in ('No conformal diagram', 'No embedding diagram', 'No curvature figure'):
            self.assertIn(step, note)

    def test_the_spacetime_diagram_stays_inside_the_first_zero_of_f(self):
        self.assertTrue((DATA / 'diagrams' / 'hauser_wave.json').exists())
        import null_rays as nr
        row = next(d for d in nr.DIAGRAMS if d.metric == 'hauser_wave')
        xmin, xmax, ymin, ymax = row.box
        # X = r and Y = t + r, so t = Y - X runs over [ymin - xmax, ymax - xmin].
        self.assertLess(max(abs(ymin - xmax), abs(ymax - xmin)), 5.46728770730835)
        self.assertEqual(row.fixed, {'x': '1/2', 'y': '0'})


if __name__ == '__main__':
    unittest.main()
