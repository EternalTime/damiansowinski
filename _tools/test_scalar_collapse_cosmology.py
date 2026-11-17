"""The scalar field universe of Husain, Martinez and Nunez and Fonarev's family:
.venv.noindex/bin/python -m unittest _tools/test_scalar_collapse_cosmology.py

What the captions and the drawings state, held on the published files. Every drawing takes b = 1 and
c = 1; the massless members are drawn at eta_0 = 2b and Fonarev's at gamma = 4/5, s = 9/7. The
apparent horizon is eta = 2 s r^2 f^(1 - gamma)/(2r - b(1 + gamma)), Husain, Martinez and Nunez's
(11) at s = 1/2, and the tortoise radius r_* is the integral of f^(-gamma) dr from r = b. The
generators check the rays, the surfaces and the scalars against the published metric; these tests
hold the files to the numbers and the shapes the captions put in words, and the published scalar
field and potential to the field equations the conventions state, which needs sympy.
"""
import importlib.util
import json
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
HAS_NUMERIC = all(importlib.util.find_spec(m) is not None for m in ("numpy", "scipy"))
HAS_SYMPY = importlib.util.find_spec("sympy") is not None
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))

ID = "scalar_collapse_cosmology"
ROOT3 = math.sqrt(3) / 2


def load(folder):
    return json.loads((DATA / folder / f"{ID}.json").read_text(encoding="utf-8"))


def horizon(r, gamma, s=0.5):
    f = (r - 1) / r
    return 2 * s * r * r * f ** (1 - gamma) / (2 * r - (1 + gamma))


def chart_point(view, at):
    X0, X1, Y0, Y1 = view["box"]
    return X0 + at[0] * (X1 - X0), Y0 + at[1] * (Y1 - Y0)


def flat(system, vid):
    return next(v for v in load("diagrams")["systems"][system] if v["id"] == vid)


def marked(view, kind):
    return [chart_point(view, at) for m in view["markers"] if m["kind"] == kind
            for line in m.get("lines", []) for at in line]


class Horizons(unittest.TestCase):
    """The apparent horizons of the spacetime diagrams and the numbers their captions give."""

    def turning(self, gamma, s=0.5):
        """The radii and times where eta_AH(r) turns, on a fine grid of r."""
        n = 400000
        rs = [1 + 4 * (k + 1) / n for k in range(n)]
        hs = [horizon(r, gamma, s) for r in rs]
        return [(rs[k], hs[k]) for k in range(1, n - 1) if (hs[k] - hs[k - 1]) * (hs[k + 1] - hs[k]) < 0]

    def test_at_root_three_halves_the_horizon_cuts_each_moment_three_times_between_1_67_and_3_58(self):
        """The caption's numbers: eta_AH rises to 3.58 b at r = 1.0125 b, within 0.013 b of r = b, falls
        to 1.67 b at r = 1.72 b and then grows without bound."""
        (r_max, h_max), (r_min, h_min) = self.turning(ROOT3)
        self.assertEqual(round(h_max, 2), 3.58)
        self.assertLess(r_max - 1, 0.013)
        self.assertEqual(round(h_min, 2), 1.67)
        self.assertEqual(round(r_min, 2), 1.72)
        self.assertGreater(horizon(1e6, ROOT3), 1e5)
        caption = flat("conformal_time", "positive")["caption"][0]
        for number in ("1.67\\,b", "3.58\\,b", "0.013\\,b"):
            self.assertIn(number, caption)

    def test_at_minus_root_three_halves_one_horizon_leaves_the_big_bang_and_grows_as_r_over_two(self):
        self.assertEqual(self.turning(-ROOT3), [])
        self.assertLess(horizon(1 + 1e-9, -ROOT3), 1e-12)
        self.assertAlmostEqual(horizon(1e8, -ROOT3) / 1e8, 0.5, places=6)
        self.assertIn("r \\approx 2\\eta", flat("conformal_time", "negative")["caption"][0])

    def test_every_drawn_apparent_horizon_is_the_closed_form(self):
        """Each point of the marker, read in the chart's own time and carried to eta, lies on the closed
        form, away from r = b where the curve stands nearly upright."""
        eta = {"conformal_time": lambda x: x, "comoving_time": lambda t: 2 * (t / (4 / 3)) ** (2 / 3),
               "fonarev_conformal": lambda x: x, "fonarev": lambda T: (9 / 7) * math.exp(7 * T / 9)}
        cases = (("conformal_time", "positive", ROOT3, 0.5), ("conformal_time", "negative", -ROOT3, 0.5),
                 ("comoving_time", "radial", ROOT3, 0.5), ("fonarev_conformal", "radial", 0.8, 9 / 7),
                 ("fonarev", "radial", 0.8, 9 / 7))
        for system, vid, gamma, s in cases:
            points = [(r, x) for r, x in marked(flat(system, vid), "apparent") if r > 1.02]
            self.assertGreater(len(points), 5, (system, vid))
            for r, x in points:
                self.assertAlmostEqual(eta[system](x), horizon(r, gamma, s), delta=0.01, msg=(system, vid, r))

    def test_the_positive_horizon_runs_from_beside_r_b_through_its_least_time(self):
        view = flat("conformal_time", "positive")
        points = marked(view, "apparent")
        start = min(points)
        self.assertLess(start[0] - 1, 0.02)
        self.assertAlmostEqual(start[1], 3.58, delta=0.02)
        low = min(points, key=lambda p: p[1])
        self.assertAlmostEqual(low[1], 1.67, delta=0.02)
        self.assertAlmostEqual(low[0], 1.72, delta=0.05)


class Singularities(unittest.TestCase):
    def test_the_big_bang_and_r_b_are_marked_singular_on_every_plane_that_reaches_them(self):
        for system, vid in (("conformal_time", "positive"), ("conformal_time", "negative"), ("comoving_time", "radial")):
            singular, = [m for m in flat(system, vid)["markers"] if m["kind"] == "singular"]
            self.assertEqual(sorted(singular["edges"]), ["bottom", "left"], (system, vid))
        for system, edges in (("fonarev_conformal", ["bottom"]), ("fonarev", [])):
            view = flat(system, "radial")
            singular, = [m for m in view["markers"] if m["kind"] == "singular"]
            self.assertEqual(singular["edges"], edges, system)
            points = marked(view, "singular")
            self.assertGreaterEqual(len(points), 2, system)
            for r, _ in points:
                self.assertAlmostEqual(r, 1.0, delta=1e-3, msg=system)

    def test_the_conformal_diagram_is_the_triangle_of_the_expanding_branch(self):
        """The big bang on T = 0 from i0 to the corner and r = b on X = 0 up to i+, one zigzag, and
        null infinity from i+ to i0."""
        for view in load("conformal")["views"]:
            zig, = [layer["points"] for layer in view["layers"] if layer["class"] == "singular"]
            for got, want in zip(zig, [[math.pi, 0], [0, 0], [0, math.pi]]):
                self.assertAlmostEqual(got[0], want[0], places=3)
                self.assertAlmostEqual(got[1], want[1], places=3)
            scri, = [layer["points"] for layer in view["layers"] if layer["class"] == "scri"]
            for X, T in scri:
                self.assertAlmostEqual(X + T, math.pi, places=3)
            self.assertIn("each point in the diagram a 2-sphere", view["caption"][0])


@unittest.skipUnless(HAS_NUMERIC, "needs numpy and scipy")
class Moments(unittest.TestCase):
    def setUp(self):
        import scalar_collapse_cosmology_drawings as drawings
        self.d = drawings

    def test_each_moment_is_eta_eta_0_on_every_drawing_of_its_member(self):
        """In the planes, the line eta = eta_0, or t = t_0, across the embedding's reach; in the conformal
        diagram, the same line carried by p, q = arctan((eta -+ r_*)/4b)."""
        reach = {}
        for view in load("embedding")["views"]:
            xs = [p[0] for piece in view["surfaces"][0]["pieces"] for p in piece["points"]]
            reach[view["id"]] = (min(xs), max(xs))
        for system, vid, level in (("conformal_time", "positive", 2.0), ("conformal_time", "negative", 2.0),
                                   ("comoving_time", "radial", 4 / 3)):
            view = flat(system, vid)
            mark, = view["slices"]
            self.assertEqual(mark["view"], "negative" if vid == "negative" else "positive")
            line, = mark["lines"]
            ends = [chart_point(view, at) for at in line]
            for r, x in ends:
                self.assertAlmostEqual(x, level, places=3)
            lo, hi = reach[mark["view"]]
            self.assertAlmostEqual(min(r for r, _ in ends), max(lo, 1.0), delta=1e-3)
            self.assertAlmostEqual(max(r for r, _ in ends), hi, delta=1e-3)
        for view in load("conformal")["views"]:
            if view["id"] not in ("positive", "negative", "comoving"):
                self.assertEqual(view.get("slices", []), [], view["id"])
                continue
            mark, = view["slices"]
            self.assertEqual(mark["view"], "negative" if view["id"] == "negative" else "positive")
            for X, T in mark["lines"][0]:
                # With l = 4b, tan p + tan q = eta/2 and tan q - tan p = r_*/2.
                p, q = (T - X) / 2, (T + X) / 2
                self.assertAlmostEqual(2 * (math.tan(p) + math.tan(q)), 2.0, delta=2e-3 * (1 + math.tan(q) ** 2))
                self.assertGreaterEqual(2 * (math.tan(q) - math.tan(p)), -1e-3)

    def test_the_level_circle_of_the_positive_member_is_r_1_005_b(self):
        """g_rr = (d rho/dr)^2 on r_e = (1 + gamma)^2 b/4 gamma, where the flat piece meets the piece in
        Minkowski space; the negative member is in Minkowski space throughout."""
        gamma = ROOT3
        level = (1 + gamma) ** 2 / (4 * gamma)
        self.assertEqual(round(level, 3), 1.005)
        views = {v["id"]: v for v in load("embedding")["views"]}
        near, far = views["positive"]["surfaces"][0]["pieces"]
        self.assertEqual(near.get("space"), "minkowski")
        self.assertNotIn("space", far)
        self.assertAlmostEqual(near["points"][-1][0], level, places=12)
        self.assertAlmostEqual(far["points"][0][0], level, places=12)
        self.assertIn("1.005\\,b", views["positive"]["caption"][0])
        whole, = views["negative"]["surfaces"][0]["pieces"]
        self.assertEqual(whole.get("space"), "minkowski")
        for view in views.values():
            for piece in view["surfaces"][0]["pieces"]:
                for r, rho, _ in piece["points"]:
                    self.assertAlmostEqual(rho, self.d.areal(r, -ROOT3 if view["id"] == "negative" else ROOT3),
                                           delta=1e-6)


@unittest.skipUnless(HAS_SYMPY, "needs sympy")
class Conventions(unittest.TestCase):
    """The scalar fields and the potentials the conventions print, held to the field equations they
    state, R_ab = (8 pi G/c^4)(d_a phi d_b phi + g_ab V), with the published Ricci tensor."""

    def setUp(self):
        import sympy
        import verify_metrics
        self.sp, self.vm = sympy, verify_metrics
        metric = json.loads((DATA / "metrics" / f"{ID}.json").read_text(encoding="utf-8"))
        self.charts = {c["id"]: c for c in metric["coordinates"]}

    def test_the_conventions_print_the_fields_the_ricci_tensor_holds(self):
        sp = self.sp
        for system, gammas in (("conformal_time", (ROOT3, -ROOT3)), ("fonarev_conformal", (0.8, -0.6))):
            entry = self.charts[system]
            symbols = [p["symbol"] for p in entry["parameters"]]
            reader = self.vm.Reader(entry["coords"], symbols, (), held=self.vm.HELD.get((ID, system), ()),
                                    rates=self.vm.RATES.get((ID, system)))
            eta, r = reader.symbol["\\eta"], reader.symbol["r"]
            P = reader.parameters
            ricci = {tuple(e["indices"]): reader(e["value"]) for e in entry["ricci_tensor"]["variants"]["ll"]["nonzero"]}
            g = {tuple(e["indices"]): reader(e["value"]) for e in entry["metric_components"]}
            f = 1 - P["b"] / r
            for gamma in gammas:
                if system == "conformal_time":
                    self.assertIn("\\sqrt{3}\\ln(\\eta/\\eta_0) + (\\gamma/\\sqrt{3})\\ln f", entry["convention"])
                    # phi in units of sqrt(c^4/16 pi G): (8 pi G/c^4) phi^2 = X^2/2.
                    X = sp.sqrt(3) * sp.log(eta / P["eta_0"]) + P["gamma"] * sp.log(f) / sp.sqrt(3)
                    half, V, on = sp.Rational(1, 2), 0, {}
                else:
                    self.assertIn("\\gamma\\ln a/\\sqrt{1-\\gamma^2}", entry["convention"])
                    self.assertIn("8\\pi GV/c^4 = s(2s-1)/(\\eta^2a^2f^{\\gamma})", entry["convention"])
                    s = (1 - P["gamma"] ** 2) / (2 * P["gamma"] ** 2 - 1)
                    a = (eta / P["eta_0"]) ** P["s"]
                    # phi in units of sqrt(c^4/4 pi G): (8 pi G/c^4) phi^2 = 2 X^2.
                    X = (P["gamma"] * sp.log(a) / sp.sqrt(1 - P["gamma"] ** 2)
                         + sp.sqrt(1 - P["gamma"] ** 2) * sp.log(f) / 2)
                    half, V, on = 2, P["s"] * (2 * P["s"] - 1) / (eta ** 2 * a ** 2 * f ** P["gamma"]), {P["s"]: s}
                at = {P["b"]: sp.Rational(1), P["eta_0"]: sp.Rational(3, 2), P["gamma"]: sp.nsimplify(gamma),
                      eta: sp.Rational(7, 5), r: sp.Rational(23, 10)}
                for (i, j), value in ricci.items():
                    want = half * sp.diff(X, reader.symbol[i]) * sp.diff(X, reader.symbol[j]) + V * g.get((i, j), 0)
                    got = (value - want).subs(reader.held).doit().subs(on).subs(reader.held).doit()
                    self.assertLess(abs(complex(sp.N(got.subs(on).subs(at), 30))), 1e-20, (system, gamma, i, j))


if __name__ == "__main__":
    unittest.main()
