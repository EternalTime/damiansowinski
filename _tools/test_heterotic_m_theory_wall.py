"""The domain wall of heterotic M-theory:
.venv.noindex/bin/python -m unittest discover -s _tools

What the published files state, held to the articles. Each chart's published line element is read
through the checker's own Reader at points and held to its source: Lukas, Ovrut, Stelle and Waldram's
(sol), Lehners, McFadden and Turok's (domainwall) and (staticbulk) in their own alpha, c, mu and time,
and Lukas, Ovrut and Waldram's (burt6), (burt15) and (burt16). The static bulk's horizon is held to the
root of Lehners, McFadden and Turok's f, the History's exponents to Lukas, Ovrut and Waldram's numbers,
and the drawings to one member of the static walls and to the branch the History quotes. The tests that
read the mathematics need sympy and are skipped where it is absent.
"""
import importlib.util
import json
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
HAS_SYMPY = importlib.util.find_spec("sympy") is not None
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))
METRIC = json.loads((DATA / "metrics" / "heterotic_m_theory_wall.json").read_text(encoding="utf-8"))


def chart(system):
    return next(c for c in METRIC["coordinates"] if c["id"] == system)


def drawing(kind):
    return json.loads((DATA / kind / "heterotic_m_theory_wall.json").read_text(encoding="utf-8"))


@unittest.skipUnless(HAS_SYMPY, "needs sympy")
class Charts(unittest.TestCase):
    """Each published chart at points against its source, in the chart x^0 = ct."""

    @classmethod
    def setUpClass(cls):
        import sympy as sp
        import heterotic_m_theory_wall as hw
        import verify_metrics as vm
        cls.sp, cls.hw, cls.vm = sp, hw, vm

    def published(self, system):
        c = chart(system)
        key = (self.hw.METRIC, system)
        reader = self.vm.Reader(c["coords"], [p["symbol"] for p in c["parameters"]], (),
                                held=self.vm.HELD.get(key, ()), rates=self.vm.RATES.get(key))
        g = self.vm.metric_from_line_element(reader, c["line_element"].replace("c^2", ""), c["coords"])
        g = self.sp.Matrix(g).xreplace(reader.held).doit() if reader.held else self.sp.Matrix(g)
        return reader, g, [reader.symbol[x] for x in c["coords"]], reader.parameters

    def same(self, mine, theirs, at, what):
        for a in range(5):
            for b in range(5):
                gap = self.sp.N((mine[a, b] - theirs[a, b]).subs(at), 30)
                self.assertLess(abs(complex(gap)), 1e-20, f"{what}, slot {a}{b}")

    def test_static_is_their_solution(self):
        """LOSW (sol): a = a_0 H^{1/2}, b = b_0 H^2, H = (sqrt 2/3) alpha |y| + c_0, on both sides of y = 0."""
        sp = self.sp
        _, g, X, P = self.published("static")
        al, a0, b0, c0 = P["alpha"], P["a_0"], P["b_0"], P["c_0"]
        H = sp.sqrt(2) * al * sp.Abs(X[4]) / 3 + c0
        theirs = sp.diag(-a0 ** 2 * H, a0 ** 2 * H, a0 ** 2 * H, a0 ** 2 * H, b0 ** 2 * H ** 4)
        for y in (sp.Rational(2, 5), sp.Rational(-3, 4)):
            self.same(g, theirs, {al: sp.Rational(7, 10), a0: sp.Rational(13, 10), b0: sp.Rational(9, 10),
                                  c0: sp.Rational(11, 10), X[4]: y}, f"static at y = {y}")

    def test_conformal_is_their_domain_wall(self):
        """LMT (domainwall): h^{2/5}(eta + dy^2), h = 5 alpha_LMT (y + y_0) + c, alpha_LMT = alpha/3 sqrt 2."""
        sp = self.sp
        _, g, X, P = self.published("conformal")
        al, w0, h0 = P["alpha"], P["w_0"], P["h_0"]
        h = 5 * (al / (3 * sp.sqrt(2))) * (X[4] + w0) + h0
        theirs = sp.diag(-1, 1, 1, 1, 1) * h ** sp.Rational(2, 5)
        self.same(g, theirs, {al: sp.Rational(7, 10), w0: sp.Rational(1, 2), h0: sp.Rational(11, 10),
                              X[4]: sp.Rational(1, 5)}, "conformal")

    def test_static_bulk_is_their_static_bulk(self):
        """LMT (staticbulk): -f dt_L^2 + r^12 dr^2/f + r^2 dx^2, f = alpha_LMT^2 r^2 - mu r^4, with their time
        t_L = (c/alpha_LMT) t, here x^0/alpha_LMT, and mu = lambda^2 alpha_LMT^2, their lambda = sqrt(mu)/alpha."""
        sp = self.sp
        _, g, X, P = self.published("static_bulk")
        al, lam = P["alpha"], P["lambda"]
        aL = al / (3 * sp.sqrt(2))
        mu = lam ** 2 * aL ** 2
        r = X[4]
        f = aL ** 2 * r ** 2 - mu * r ** 4
        theirs = sp.diag(-f / aL ** 2, r ** 2, r ** 2, r ** 2, r ** 12 / f)
        self.same(g, theirs, {al: sp.Rational(7, 10), lam: sp.Rational(2, 5), r: sp.Rational(3, 5)}, "static bulk")

    def test_static_bulk_horizon_is_the_root_of_f(self):
        """f = alpha^2 r^2 - mu r^4 vanishes at r = alpha/sqrt(mu), which is 1/lambda: the published domain
        ends there, and at the reciprocal LMT print, sqrt(mu)/alpha = lambda, g_tt does not vanish."""
        sp = self.sp
        _, g, X, P = self.published("static_bulk")
        lam = P["lambda"]
        at = {P["alpha"]: sp.Rational(7, 10), lam: sp.Rational(2, 5)}
        self.assertEqual(sp.simplify(g[0, 0].subs({**at, X[4]: 1 / at[lam]})), 0)
        self.assertNotEqual(sp.simplify(g[0, 0].subs({**at, X[4]: at[lam]})), 0)
        self.assertIn("r \\in (0, 1/\\lambda)", chart("static_bulk")["domains"])
        self.assertIn("r = 1/\\lambda \\;\\text{(the cosmological horizon)}", chart("static_bulk")["domains"])

    def test_rolling_is_their_cosmological_solution(self):
        """LOW (burt6), (burt15): N = n a(y), a(t, y) = alpha(t) a(y), b(t, y) = beta(t) b(y), n = 1,
        alpha = A (t - t_0)^p, beta = B (t - t_0)^q, with LOSW's a(y) and b(y)."""
        sp = self.sp
        _, g, X, P = self.published("rolling")
        al, a0, b0, c0, A, B, t0, p, q = (P[n] for n in ("alpha", "a_0", "b_0", "c_0", "A", "B", "t_0", "p", "q"))
        H = sp.sqrt(2) * al * sp.Abs(X[4]) / 3 + c0
        since = X[0] - t0
        ay, by = a0 * H ** sp.Rational(1, 2), b0 * H ** 2
        alpha, beta = A * since ** p, B * since ** q
        theirs = sp.diag(-ay ** 2, alpha ** 2 * ay ** 2, alpha ** 2 * ay ** 2, alpha ** 2 * ay ** 2, beta ** 2 * by ** 2)
        for y in (sp.Rational(2, 5), sp.Rational(-3, 4)):
            self.same(g, theirs, {al: sp.Rational(7, 10), a0: sp.Rational(13, 10), b0: sp.Rational(9, 10),
                                  c0: sp.Rational(11, 10), A: sp.Rational(4, 5), B: sp.Rational(6, 5),
                                  t0: sp.Rational(1, 10), p: sp.Rational(3, 7), q: sp.Rational(-2, 5),
                                  X[0]: sp.Rational(9, 10), X[4]: y}, f"rolling at y = {y}")


class History(unittest.TestCase):
    """The History's numbers against Lukas, Ovrut and Waldram's (burt16) and the published charts."""

    def test_the_exponents(self):
        p = [3 / 11 * (1 - s * 4 / (3 * math.sqrt(3))) for s in (1, -1)]
        q = [2 / 11 * (1 + s * 2 * math.sqrt(3)) for s in (1, -1)]
        self.assertEqual([round(x, 2) for x in p + q], [0.06, 0.48, 0.81, -0.45])
        h = METRIC["history"]
        self.assertIn("p = \\frac{3}{11}\\left(1 \\mp \\frac{4}{3\\sqrt{3}}\\right)", h)
        self.assertIn("q = \\frac{2}{11}\\left(1 \\pm 2\\sqrt{3}\\right)", h)
        described = {x["symbol"]: x["description"] for x in chart("rolling")["parameters"]}
        self.assertIn("\\frac{3}{11}\\left(1 \\mp \\frac{4}{3\\sqrt{3}}\\right)", described["p"])
        self.assertIn("\\frac{2}{11}\\left(1 \\pm 2\\sqrt{3}\\right)", described["q"])

    def test_the_warp_factor(self):
        h = METRIC["history"]
        self.assertIn("H = \\frac{\\sqrt{2}}{3}\\alpha|y| + c_0", h)
        self.assertIn("V = b_0H^3", h)
        symbols = [x["symbol"] for x in chart("static")["parameters"]]
        self.assertIn("H = \\frac{\\sqrt{2}}{3}\\alpha|y| + c_0", symbols)
        self.assertIn("V = b_0 H^3", symbols)


class Drawings(unittest.TestCase):
    """The drawings of the static walls are of one member, and the moving walls of the quoted branch."""

    def test_one_member(self):
        import heterotic_m_theory_wall_drawings as hd
        k = math.sqrt(2) / 3
        h1 = 1 + k * hd.PIRHO
        # With b_0 = a_0^6 = 1 and c_0 = 1, h = H^{5/2}, so the conformal width of the interval is twice w_0.
        self.assertAlmostEqual(hd.W0, (h1 ** 2.5 - 1) / (5 * k), places=12)
        self.assertAlmostEqual(float(hd.W(hd.PIRHO)), 2 * hd.W0, places=12)
        self.assertEqual(hd.STATIC_PARAMS["rho"], "6/(5*pi)")
        self.assertAlmostEqual(math.pi * 6 / (5 * math.pi), hd.PIRHO, places=12)
        self.assertEqual(hd.CONFORMAL_PARAMS["h_0"], "1")
        if HAS_SYMPY:
            import sympy as sp
            self.assertAlmostEqual(float(sp.sympify(hd.CONFORMAL_PARAMS["w_0"])), hd.W0, places=12)

    def test_the_quoted_branch(self):
        """The branch on which the walls grow after t_0 while the interval shrinks, p > 0 and q < 0."""
        import heterotic_m_theory_wall_drawings as hd
        self.assertGreater(hd.P_DOWN, 0)
        self.assertLess(hd.Q_DOWN, 0)
        self.assertAlmostEqual(hd.P_DOWN, 3 / 11 * (1 + 4 / (3 * math.sqrt(3))), places=14)
        self.assertAlmostEqual(hd.Q_DOWN, 2 / 11 * (1 - 2 * math.sqrt(3)), places=14)

    def test_every_kind_is_drawn(self):
        for kind in ("diagrams", "conformal", "embedding", "curvature"):
            self.assertTrue((DATA / kind / "heterotic_m_theory_wall.json").exists(), kind)


if __name__ == "__main__":
    unittest.main()
