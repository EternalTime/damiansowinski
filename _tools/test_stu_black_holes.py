"""The black holes of anti-de Sitter space with three charges:
.venv.noindex/bin/python -m unittest discover -s _tools

What the published files state, held to the articles. Each chart's published line element is read
through the checker's own Reader at points and held to its source: Behrndt, Cvetic and Sabra's (390)
with (110) and (170) for the three horizons, Myers and Tafjord's rho^2 = r^2 + q_3, and Son and
Starinets's (metric_u_3) with (identif). The History's numbers are held to the spherical chart: the
mass is BCS's with G = pi/4, the hole with no charge has its least temperature sqrt 2/(pi L) at
mu = 3L^2/4, and the brane's temperature is Son and Starinets's. Every drawing of the spherical and
continued charts is of the one member the embedding diagram draws. The tests that read the
mathematics need sympy and are skipped where it is absent.
"""
import importlib.util
import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
HAS_SYMPY = importlib.util.find_spec("sympy") is not None
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))
METRIC = json.loads((DATA / "metrics" / "stu_black_holes.json").read_text(encoding="utf-8"))
SYSTEMS = ("sphere", "flat", "hyperbolic", "continued", "brane")


def chart(system):
    return next(c for c in METRIC["coordinates"] if c["id"] == system)


def drawing(kind):
    return json.loads((DATA / kind / "stu_black_holes.json").read_text(encoding="utf-8"))


@unittest.skipUnless(HAS_SYMPY, "needs sympy")
class Charts(unittest.TestCase):
    """Each published chart at points against its source."""

    @classmethod
    def setUpClass(cls):
        import sympy as sp
        import stu_black_holes as stu
        import verify_metrics as vm
        cls.sp, cls.stu, cls.vm = sp, stu, vm

    def read(self, system):
        reader = self.stu.reader(system)
        c = chart(system)
        g = self.vm.metric_from_line_element(reader, c["line_element"].replace("c^2", ""), c["coords"])
        return reader, self.sp.Matrix(g).xreplace(reader.held).doit()

    def at(self, reader, system):
        sp, P, X = self.sp, reader.parameters, [reader.symbol[c] for c in chart(system)["coords"]]
        if system == "brane":
            point = {P["r_+"]: sp.Rational(6, 5), P["L"]: sp.Rational(3, 2), P["kappa_1"]: sp.Rational(1, 2),
                     P["kappa_2"]: sp.Rational(1, 5), P["kappa_3"]: sp.Rational(1, 10), X[4]: sp.Rational(2, 5)}
        else:
            point = {P["mu"]: sp.Rational(13, 10), P["q_1"]: sp.Rational(7, 10), P["q_2"]: sp.Rational(2, 5),
                     P["q_3"]: sp.Rational(1, 5), P["L"]: sp.Rational(3, 2), X[1]: sp.Rational(13, 10),
                     X[2]: sp.Rational(3, 5), X[3]: sp.Rational(11, 10), X[4]: sp.Rational(2, 5)}
        return X, point

    def source(self, reader, system, X):
        """The diagonal of the source's metric in the chart's own coordinates."""
        sp, P = self.sp, reader.parameters
        if system == "brane":
            u, rp, L = X[4], P["r_+"], P["L"]
            kappa = [P[f"kappa_{i}"] for i in (1, 2, 3)]
            H = sp.Mul(*[1 + k * u for k in kappa])
            f = H - u ** 2 * sp.Mul(*[1 + k for k in kappa])
            plane = rp ** 2 * H ** sp.Rational(1, 3) / (L ** 2 * u)
            return [-rp ** 2 * f / (L ** 2 * u * H ** sp.Rational(2, 3)), plane, plane, plane,
                    L ** 2 * H ** sp.Rational(1, 3) / (4 * u ** 2 * f)]
        mu, L = P["mu"], P["L"]
        q = [P[f"q_{i}"] for i in (1, 2, 3)]
        if system == "continued":
            rho = X[1]
            r2 = rho ** 2 - q[2]
            H = sp.Mul(*[1 + qi / r2 for qi in q])
            f = 1 - mu / r2 + r2 * H / L ** 2
            # ds^2 of BCS in r, with dr = rho drho/r and r^2 = rho^2 - q_3.
            sphere = H ** sp.Rational(1, 3) * r2
            return [-H ** sp.Rational(-2, 3) * f, H ** sp.Rational(1, 3) * rho ** 2 / (r2 * f), sphere,
                    sphere * sp.sin(X[2]) ** 2, sphere * sp.sin(X[2]) ** 2 * sp.sin(X[3]) ** 2]
        k = self.stu.CURVATURE[system]
        r = X[1]
        H = sp.Mul(*[1 + qi / r ** 2 for qi in q])
        f = k - mu / r ** 2 + r ** 2 * H / L ** 2
        s = H ** sp.Rational(1, 3) * r ** 2
        if system == "flat":
            return [-H ** sp.Rational(-2, 3) * f, H ** sp.Rational(1, 3) / f] + [s / L ** 2] * 3
        S = sp.sin(X[2]) if k == 1 else sp.sinh(X[2])
        return [-H ** sp.Rational(-2, 3) * f, H ** sp.Rational(1, 3) / f, s, s * S ** 2, s * S ** 2 * sp.sin(X[3]) ** 2]

    def test_each_line_element_is_its_source(self):
        for system in SYSTEMS:
            with self.subTest(system=system):
                reader, g = self.read(system)
                X, point = self.at(reader, system)
                for i, value in enumerate(self.source(reader, system, X)):
                    self.assertLess(abs(complex(self.sp.N((g[i, i] - value).subs(point), 30))), 1e-25)

    def test_the_history_states_the_published_metric_function_and_scalars(self):
        history = METRIC["history"]
        self.assertIn("$$f = k - \\frac{\\mu}{r^2} + \\frac{r^2}{L^2}H_1H_2H_3,$$", history)
        self.assertIn("$H_i = 1 + q_i/r^2$", history)
        self.assertIn("$X_i = (H_1H_2H_3)^{1/3}/H_i$", history)
        self.assertEqual(chart("sphere")["parameters"][-1]["symbol"],
                         "f = 1 - \\dfrac{\\mu}{r^2} + \\dfrac{r^2}{L^2}H_1\\,H_2\\,H_3")

    def test_the_hole_with_no_charge_is_coldest_at_three_quarters_of_l_squared(self):
        """f = 1 - mu/r^2 + r^2/L^2 with no charge; T = f'(r_+)/4 pi in units c = hbar = k_B = 1 is least,
        sqrt 2/(pi L), at r_+ = L/sqrt 2, where mu = r_+^2 + r_+^4/L^2 = 3L^2/4, as the History says."""
        sp = self.sp
        r, L = sp.symbols("r L", positive=True)
        mu = r ** 2 + r ** 4 / L ** 2                       # f(r_+) = 0
        T = sp.diff(1 - mu / sp.Symbol("x") ** 2 + sp.Symbol("x") ** 2 / L ** 2, sp.Symbol("x")).subs(sp.Symbol("x"), r) / (4 * sp.pi)
        least, = [s for s in sp.solve(sp.diff(T, r), r) if s.is_positive]
        self.assertEqual(sp.simplify(least - L / sp.sqrt(2)), 0)
        self.assertEqual(sp.simplify(T.subs(r, least) - sp.sqrt(2) / (sp.pi * L)), 0)
        self.assertEqual(sp.simplify(mu.subs(r, least) - 3 * L ** 2 / 4), 0)
        self.assertIn("$\\sqrt{2}/\\pi L$ in natural units, at $\\mu = 3L^2/4$", METRIC["history"])

    def test_the_mass_is_behrndt_cvetic_and_sabras_with_g_a_quarter_of_pi(self):
        sp = self.sp
        mu, q1, q2, q3, G, c = sp.symbols("mu q_1 q_2 q_3 G c", positive=True)
        published = sp.pi * c ** 2 * (3 * mu + 2 * q1 + 2 * q2 + 2 * q3) / (8 * G)
        self.assertEqual(sp.simplify(published.subs({G: sp.pi / 4, c: 1}) - (q1 + q2 + q3 + 3 * mu / 2)), 0)
        self.assertIn("M = \\pi c^2\\left(3\\mu + 2q_1 + 2q_2 + 2q_3\\right)/(8G)", chart("sphere")["convention"])
        self.assertIn("$M = q_1 + q_2 + q_3 + \\frac{3}{2}\\mu$", METRIC["history"])

    def test_the_brane_is_as_hot_as_son_and_starinets_say(self):
        """2 pi T = (d sqrt(-g_tt)/dl) on the horizon u = 1 of the published brane chart, in units
        hbar = k_B = c = 1, is Son and Starinets's T_0 (2 + sum kappa - kappa_1 kappa_2 kappa_3)/
        (2 sqrt(prod (1 + kappa))) with T_0 = r_+/(pi L^2)."""
        sp = self.sp
        reader, g = self.read("brane")
        X, point = self.at(reader, "brane")
        u = X[4]
        lapse, stretch = -g[0, 0], g[4, 4]
        kappa = [reader.parameters[f"kappa_{i}"] for i in (1, 2, 3)]
        rp, L = reader.parameters["r_+"], reader.parameters["L"]
        want = (rp / (sp.pi * L ** 2)) * (2 + sum(kappa) - sp.Mul(*kappa)) / (2 * sp.sqrt(sp.Mul(*[1 + k for k in kappa])))
        near = {**point, u: 1 - sp.Rational(1, 10 ** 15)}
        got = sp.N((-sp.diff(lapse, u) / (2 * sp.sqrt(lapse * stretch)) / (2 * sp.pi)).subs(near), 30)
        self.assertLess(abs(float(got - sp.N(want.subs(point), 30))), 1e-12)


class Drawings(unittest.TestCase):
    """Every drawing of the spherical and continued charts is of the member q_1 = 4 r_+^2, q_2 = 2 r_+^2,
    q_3 = r_+^2, L = 3 r_+, mu = 13 r_+^2/3, the one the embedding diagram draws and whose moment the
    spacetime and conformal diagrams of those charts show."""

    def test_every_drawing_of_the_spherical_and_continued_charts_is_of_one_member(self):
        settings = "$q_1 = 4r_+^2$, $q_2 = 2r_+^2$, $q_3 = r_+^2$, $L = 3r_+$, $\\mu = 13r_+^2/3$, $r_+ = 1$."
        self.assertEqual([v["settings"] for v in drawing("embedding")["views"]], [settings])
        self.assertEqual([(v["id"], v["settings"]) for v in drawing("conformal")["views"]],
                         [("sphere", settings), ("continued", settings)])
        self.assertEqual([(v["system"], v["settings"]) for v in drawing("curvature")["views"]],
                         [("continued", settings.replace(", $r_+ = 1$", ""))])
        diagrams = drawing("diagrams")["systems"]
        self.assertEqual(sorted(diagrams), sorted(SYSTEMS))
        for system in ("sphere", "continued"):
            with self.subTest(system=system):
                text = " ".join(v["settings"] for v in diagrams[system])
                for value in ("\\mu = \\frac{13}{3}", "q_1 = 4", "q_2 = 2", "q_3 = 1", "L = 3"):
                    self.assertIn(value, text)

    def test_the_embedding_reaches_the_level_circle(self):
        """The slice t = 0 lies level at r = 2.356 r_+, where (dR/dl)^2 = 1 for R = r (H_1H_2H_3)^{1/6}."""
        import math
        q, L, mu = (4.0, 2.0, 1.0), 3.0, 13.0 / 3.0

        def R(r):
            return r * ((1 + q[0] / r ** 2) * (1 + q[1] / r ** 2) * (1 + q[2] / r ** 2)) ** (1 / 6)

        def level(r, h=1e-6):
            H = (1 + q[0] / r ** 2) * (1 + q[1] / r ** 2) * (1 + q[2] / r ** 2)
            f = 1 - mu / r ** 2 + r ** 2 * H / L ** 2
            return ((R(r + h) - R(r - h)) / (2 * h)) ** 2 * f / H ** (1 / 3) - 1

        lo, hi = 2.0, 3.0
        for _ in range(60):
            mid = (lo + hi) / 2
            lo, hi = (mid, hi) if level(mid) < 0 else (lo, mid)
        self.assertAlmostEqual(lo, 2.356, places=3)
        self.assertIn("2.36", drawing("embedding")["views"][0]["figure"].__repr__())
        self.assertTrue(math.isfinite(lo))


if __name__ == "__main__":
    unittest.main()
