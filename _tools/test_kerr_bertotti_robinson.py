"""Podolsky and Ovcharenko's Kerr black hole in Bertotti and Robinson's uniform magnetic field:
python3 -m unittest discover -s _tools

What its texts and drawings state, held on the published files. The first class needs nothing but
the files. The second needs sympy and mpmath and is skipped where they are absent: it reads each
published chart through the checker's Reader, writes its held names out, and holds the metric to
their line element, the Einstein-Maxwell equations with their potential, the electric charge the
spin parameter's description states, the horizons, ergosurface and orbit the captions quote, the
Bertotti-Robinson curvature at m = 0, and the radii and distances of the embedded equator.
"""
import importlib.util
import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
HAS_NUMBERS = all(importlib.util.find_spec(name) is not None for name in ("sympy", "mpmath"))
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))
ME = "kerr_bertotti_robinson"
# The values every drawing of the spinning hole takes, in units of m, and of the hole with no spin.
SPIN = {"m": 1, "a": "4/5", "B": "1/4", "C": "60025/61374"}
STILL = {"m": 1, "B": "1/4", "C": "16/17"}


def load(folder):
    return json.loads((DATA / folder / f"{ME}.json").read_text(encoding="utf-8"))


class Files(unittest.TestCase):
    def test_two_charts_and_what_is_drawn(self):
        metric = load("metrics")
        self.assertEqual([c["id"] for c in metric["coordinates"]], ["boyer_lindquist", "static"])
        diagrams = load("diagrams")["systems"]
        self.assertEqual({s: [v["id"] for v in views] for s, views in diagrams.items()},
                         {"boyer_lindquist": ["radial", "equator"], "static": ["radial"]})
        self.assertEqual([v["id"] for v in load("embedding")["views"]], ["equator"])
        self.assertFalse((DATA / "conformal" / f"{ME}.json").exists())

    def test_the_conicity_is_a_constant_of_each_chart(self):
        """C is no defined name, so no value spells it out; its description gives its value."""
        for chart in load("metrics")["coordinates"]:
            symbols = [p["symbol"] for p in chart["parameters"]]
            self.assertIn("C", symbols)
            described = next(p["description"] for p in chart["parameters"] if p["symbol"] == "C")
            self.assertIn("$C = ", described)

    def test_related_both_ways(self):
        metric = load("metrics")
        kinds = {r["id"]: r["kind"] for r in metric["related"]}
        self.assertEqual(kinds["kerr"], "special_case")
        self.assertEqual(kinds["bertotti_robinson"], "special_case")
        for other in kinds:
            theirs = json.loads((DATA / "metrics" / f"{other}.json").read_text(encoding="utf-8"))
            self.assertIn(ME, [r["id"] for r in theirs["related"]], other)


@unittest.skipUnless(HAS_NUMBERS, "needs sympy and mpmath")
class Physics(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import mpmath
        import sympy as sp
        import verify_metrics as vm
        cls.sp, cls.mp, cls.vm = sp, mpmath, vm
        mpmath.mp.dps = 30
        cls.charts = {c["id"]: c for c in load("metrics")["coordinates"]}

    def chart(self, system):
        """The published chart's Reader and its metric with every held name written out, in the
        chart x^0 = ct with c = 1."""
        sp, vm = self.sp, self.vm
        entry = self.charts[system]
        reader = vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], (),
                           held=vm.HELD[(ME, system)])
        g = vm.metric_from_line_element(reader, entry["line_element"], entry["coords"]).subs(reader.c, 1)
        return reader, self.written(reader, g), entry

    def written(self, reader, value):
        for _ in range(4):
            value = value.subs(reader.held).doit()
        return value

    def numbers(self, reader, values):
        return {reader.parameters[k]: self.sp.sympify(v) for k, v in values.items()}

    def potential(self, reader, a):
        """Their potential at gamma = 0, the complex one whose real part doubled is the field's."""
        sp = self.sp
        r, th = reader.symbol["r"], reader.symbol["\\theta"]
        B, C = reader.parameters["B"], reader.parameters["C"]
        Omega = self.written(reader, reader.parameters["Omega"])
        z = r + sp.I * a * sp.cos(th)
        first = [a, 0, 0, -(r ** 2 + a ** 2) * C]
        second = [1, 0, 0, -a * sp.sin(th) ** 2 * C]
        out = [(sp.diff(Omega, r) * first[i] / z + sp.I * sp.diff(Omega, th) / sp.sin(th) * second[i] / z) / (2 * B)
               for i in range(4)]
        out[3] += (Omega - 1) * C / (2 * B)
        return out

    def test_the_spinning_chart_is_their_line_element(self):
        sp = self.sp
        reader, g, _ = self.chart("boyer_lindquist")
        r, th = reader.symbol["r"], reader.symbol["\\theta"]
        a, C = reader.parameters["a"], reader.parameters["C"]
        Sigma, P, Q, Omega = (self.written(reader, reader.parameters[n]) for n in ("Sigma", "P", "Q", "Omega"))
        first = sp.Matrix([1, 0, 0, -a * sp.sin(th) ** 2 * C])
        second = sp.Matrix([a, 0, 0, -(r ** 2 + a ** 2) * C])
        theirs = (-Q / Sigma * first * first.T + P * sp.sin(th) ** 2 / Sigma * second * second.T) / Omega ** 2
        theirs[1, 1] += Sigma / (Q * Omega ** 2)
        theirs[2, 2] += Sigma / (P * Omega ** 2)
        at = {r: sp.Rational(23, 10), th: sp.Rational(7, 10), **self.numbers(reader, SPIN)}
        for i in range(4):
            for j in range(4):
                self.assertLess(abs(sp.N((g[i, j] - theirs[i, j]).subs(at), 30)), 1e-25)

    def einstein_maxwell(self, system, values, point):
        """max |G - 2 T| over max |G|, and |R| over the same, at a point, in 30 digits."""
        sp = self.sp
        reader, g, _ = self.chart(system)
        x = [reader.symbol[n] for n in self.charts[system]["coords"]]
        at = {x[1]: point[0], x[2]: point[1], **self.numbers(reader, values)}
        a = reader.parameters["a"] if "a" in values else sp.Integer(0)

        def number(e):
            return sp.N(e.subs(at), 30)
        G = g.applyfunc(number)
        d1 = [[[number(sp.diff(g[i, j], x[k])) for k in range(4)] for j in range(4)] for i in range(4)]
        d2 = [[[[number(sp.diff(g[i, j], x[k], x[l])) for l in range(4)] for k in range(4)] for j in range(4)]
              for i in range(4)]
        Gi = G.inv()
        dGi = [-Gi * sp.Matrix(4, 4, lambda i, j: d1[i][j][k]) * Gi for k in range(4)]
        low = lambda l, j, k: d1[l][j][k] + d1[l][k][j] - d1[j][k][l]  # noqa: E731
        gamma = [[[sum(Gi[i, l] * low(l, j, k) for l in range(4)) / 2 for k in range(4)] for j in range(4)]
                 for i in range(4)]

        def dgamma(i, j, k, s):
            return sum(dGi[s][i, l] * low(l, j, k) + Gi[i, l] * (d2[l][j][k][s] + d2[l][k][j][s] - d2[j][k][l][s])
                       for l in range(4)) / 2
        ricci = sp.Matrix(4, 4, lambda j, k: sum(
            dgamma(i, j, k, i) - dgamma(i, j, i, k)
            + sum(gamma[i][i][l] * gamma[l][j][k] - gamma[i][k][l] * gamma[l][j][i] for l in range(4))
            for i in range(4)))
        scalar = sum(Gi[i, j] * ricci[i, j] for i in range(4) for j in range(4))
        A = self.potential(reader, a)
        F = sp.Matrix(4, 4, lambda i, j: 2 * sp.re(number(sp.diff(A[j], x[i]) - sp.diff(A[i], x[j]))))
        square = sum((Gi * F * Gi)[i, j] * F[i, j] for i in range(4) for j in range(4))
        stress = sp.Matrix(4, 4, lambda i, j: sum(F[i, k] * (F * Gi)[j, k] for k in range(4)) - G[i, j] * square / 4)
        einstein = ricci - G * scalar / 2
        size = max(abs(e) for e in einstein)
        return max(abs(e) for e in einstein - 2 * stress) / size, abs(scalar) / size

    def test_both_charts_solve_the_einstein_maxwell_equations(self):
        sp = self.sp
        for system, values in (("boyer_lindquist", SPIN), ("static", STILL)):
            for point in ((sp.Rational(5, 2), sp.Rational(11, 10)), (sp.Rational(9, 2), sp.Rational(2, 5))):
                residual, trace = self.einstein_maxwell(system, values, point)
                self.assertLess(residual, 1e-24, (system, point))
                self.assertLess(trace, 1e-24, (system, point))

    def test_the_charge_the_spin_parameter_states(self):
        """The flux of the dual field through a sphere round the hole is the charge
        C m a B sqrt(I_2)/I_1 at every radius, with no magnetic charge: Ovcharenko and Podolsky's
        e_s times C, in this m. At a = 0 there is no charge."""
        sp, mp = self.sp, self.mp
        reader, g, _ = self.chart("boyer_lindquist")
        r, th = reader.symbol["r"], reader.symbol["\\theta"]
        at = self.numbers(reader, SPIN)
        a = reader.parameters["a"]
        A = self.potential(reader, a)
        x = [reader.symbol[n] for n in self.charts["boyer_lindquist"]["coords"]]
        gn = sp.lambdify((r, th), g.subs(at), "mpmath")
        F = [[sp.lambdify((r, th), (sp.diff(A[j], x[i]) - sp.diff(A[i], x[j])).subs(at), "mpmath") for j in range(4)]
             for i in range(4)]

        def flux(radius, dual):
            def integrand(angle):
                G = gn(radius, angle)
                f = lambda i, j: 2 * mp.re(F[i][j](radius, angle))  # noqa: E731
                if not dual:
                    return f(2, 3)
                D = G[0, 0] * G[3, 3] - G[0, 3] ** 2
                return mp.sqrt(G[2, 2] / G[1, 1]) * (G[3, 3] * f(0, 1) - G[0, 3] * f(3, 1)) / mp.sqrt(-D)
            # The angle phi runs over 2 pi, so the integral over it is 2 pi, and Gauss's 4 pi divides.
            return mp.quad(integrand, [mp.mpf("1e-12"), mp.pi / 2, mp.pi - mp.mpf("1e-12")]) / 2
        m, B, C = (sp.sympify(SPIN[k]) for k in ("m", "B", "C"))
        a = sp.sympify(SPIN["a"])
        I1, I2 = 1 - B ** 2 * a ** 2 / 2, 1 - B ** 2 * a ** 2
        stated = float(C * m * a * B * sp.sqrt(I2) / I1)
        for radius in (mp.mpf("2.5"), mp.mpf(7)):
            self.assertAlmostEqual(float(flux(radius, True)), stated, places=12)
            self.assertAlmostEqual(float(flux(radius, False)), 0.0, places=12)

    def test_the_numbers_the_captions_quote(self):
        """r_+ = 1.684 m, r_- = 0.405 m and g_tt = 0 at 2.021 m on the equator of the spinning hole,
        where Q = a^2; r_h = 32m/15 and the innermost stable circular orbit at 3 r_h = 32m/5 for the
        hole with no spin, found from the published equator as the least angular momentum of a
        circular orbit."""
        sp, mp = self.sp, self.mp
        reader, g, _ = self.chart("boyer_lindquist")
        r, th = reader.symbol["r"], reader.symbol["\\theta"]
        at = self.numbers(reader, SPIN)
        grr = sp.lambdify(r, (1 / g[1, 1]).subs(at).subs(th, sp.Rational(7, 10)), "mpmath")
        self.assertAlmostEqual(float(mp.findroot(grr, 1.7)), 1.6844809, places=6)
        self.assertAlmostEqual(float(mp.findroot(grr, 0.4)), 0.4052570, places=6)
        gtt = sp.lambdify(r, g[0, 0].subs(at).subs(th, sp.pi / 2), "mpmath")
        self.assertAlmostEqual(float(mp.findroot(gtt, 2.0)), 2.0210454, places=6)
        Q = self.written(reader, reader.parameters["Q"]).subs(at)
        self.assertAlmostEqual(float(Q.subs(r, sp.Float("2.0210454", 30))), 0.64, places=5)

        reader, g, _ = self.chart("static")
        r, th = reader.symbol["r"], reader.symbol["\\theta"]
        at = {**self.numbers(reader, STILL), th: sp.pi / 2}
        f = -g[0, 0].subs(at)
        self.assertEqual(sp.simplify(self.written(reader, reader.parameters["f"]).subs(at).subs(r, sp.Rational(32, 15))), 0)
        # A circular geodesic of the static equator has L^2 = g_phph^2 f'/(g_phph' f - g_phph f') per unit
        # mass squared, with f = -g_tt; the innermost stable one is where L^2 is least.
        h = g[3, 3].subs(at)
        L2 = h ** 2 * sp.diff(f, r) / (sp.diff(h, r) * f - h * sp.diff(f, r))
        least = mp.findroot(sp.lambdify(r, sp.diff(L2, r), "mpmath"), 6.5)
        self.assertAlmostEqual(float(least), 32 / 5, places=9)

    def test_the_massless_spinning_chart_is_bertotti_and_robinson(self):
        """At m = 0 the published Kretschmann scalar of the spinning chart is 8/e^4 with
        e = 1/(B sqrt(1 - B^2 a^2)), the radius the related entry states, at every point."""
        sp = self.sp
        reader, g, entry = self.chart("boyer_lindquist")
        r, th = reader.symbol["r"], reader.symbol["\\theta"]
        K = self.written(reader, reader(entry["kretschmann"].removeprefix("K = ")))
        values = {**SPIN, "m": 0, "C": "1/(1 - 1/25)"}
        at = self.numbers(reader, values)
        a, B = sp.Rational(4, 5), sp.Rational(1, 4)
        for point in ((sp.Rational(3, 2), sp.Rational(1, 3)), (sp.Rational(7, 3), sp.Rational(6, 5))):
            here = sp.N(K.subs(at).subs({r: point[0], th: point[1], reader.c: 1}), 30)
            self.assertLess(abs(here - 8 * B ** 4 * (1 - B ** 2 * a ** 2) ** 2), 1e-24)

    def test_the_embedded_equator(self):
        """The throat's radius 1.86 m, the rim's 3.84 m, and r = infinity within 7.56 m of the throat."""
        sp, mp = self.sp, self.mp
        reader, g, _ = self.chart("boyer_lindquist")
        r, th = reader.symbol["r"], reader.symbol["\\theta"]
        at = {**self.numbers(reader, SPIN), th: sp.pi / 2}
        rho = sp.lambdify(r, sp.sqrt(g[3, 3].subs(at)), "mpmath")
        length = sp.lambdify(r, sp.sqrt(g[1, 1].subs(at)), "mpmath")
        # The horizon, the larger root of the published g^rr, found to the working precision, so that
        # the integral of the distance starts on it and never below it.
        rp = mp.findroot(sp.lambdify(r, (1 / g[1, 1]).subs(at), "mpmath"), mp.mpf("1.7"))
        self.assertAlmostEqual(float(rho(rp)), 1.8607773, places=5)
        self.assertAlmostEqual(float(rho(mp.mpf(10) ** 8)), 3.8380256, places=5)
        self.assertAlmostEqual(float(mp.quad(length, [rp, 10, mp.inf])), 7.5602373, places=6)
        views = load("embedding")["views"]
        points = views[0]["surfaces"][0]["pieces"][0]["points"]
        self.assertAlmostEqual(min(p[1] for p in points), 1.8607773, places=4)
