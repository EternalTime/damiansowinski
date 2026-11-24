"""The Dvali-Gabadadze-Porrati braneworld:
.venv.noindex/bin/python -m unittest discover -s _tools

What the published files state, held to the articles. Each chart's published line element is read
through the checker's own Reader at points and held to its source: Deffayet's (bulkmet) on both
branches for a test scale factor, the same at a = r_c cosh(ct/r_c) and k = 1 for the global chart and at
a = e^{ct/r_c} and k = 0 for the planar one, and his flat chart (canon) with the hyperboloid (hyperbol)
for the inertial one. Each published chart passes its own checks again, among them the junction
condition with the curvature term of the wall, which is the Friedmann equation of either branch. The
History's numbers are held to the Friedmann equation, and the drawings to the member they draw. The tests
that read the mathematics need sympy and are skipped where it is absent.
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
METRIC = json.loads((DATA / "metrics" / "dgp_wall.json").read_text(encoding="utf-8"))
FRIEDMANN = "H^2 + kc^2/a^2 - (\\epsilon c/r_c)\\sqrt{H^2 + kc^2/a^2} = 8\\pi G\\rho/3"


def chart(system):
    return next(c for c in METRIC["coordinates"] if c["id"] == system)


@unittest.skipUnless(HAS_SYMPY, "needs sympy")
class Charts(unittest.TestCase):
    """Each published chart at points against its source, in the chart x^0 = ct."""

    @classmethod
    def setUpClass(cls):
        import sympy as sp
        import verify_metrics as vm
        cls.sp, cls.vm = sp, vm

    def published(self, system):
        c = chart(system)
        reader = self.vm.Reader(c["coords"], [p["symbol"] for p in c["parameters"]], ())
        g = self.vm.metric_from_line_element(reader, c["line_element"].replace("c^2dt^2", "dt^2").replace("c^2dT^2", "dT^2")
                                             .replace("{ct}", "{t}").replace("2ct", "2t"), c["coords"])
        g = self.sp.Matrix(g).xreplace(reader.held).doit() if reader.held else self.sp.Matrix(g)
        return reader, g, [reader.symbol[x] for x in c["coords"]], reader.parameters

    def same(self, mine, theirs, at, what):
        for a in range(5):
            for b in range(5):
                gap = self.sp.N((mine[a, b] - theirs[a, b]).subs(at), 30)
                self.assertLess(abs(complex(gap)), 1e-20, f"{what}, slot {a}{b}")

    def test_gaussian_normal_is_deffayets_bulk_on_both_branches(self):
        """Deffayet (bulkmet): a_D = a_0 + eps|y_D| sqrt(a_0'^2 + k), n = 1 + eps|y_D| a_0''/sqrt(a_0'^2 + k), the
        branch eps = +1 on y = |y_D| >= 0 and the branch eps = -1 on y = -|y_D| <= 0, for a test a(t)."""
        sp = self.sp
        reader, g, X, P = self.published("gaussian_normal")
        t, y, chi = X[0], X[1], X[2]
        a = reader.functions["a"][0]
        test = 2 + sp.cosh(t) / 3 + t ** 2 / 7
        g = g.subs(a, test).doit()
        for kv in (-1, 0, 1):
            for eps in (1, -1):
                distance = sp.Rational(3, 7)
                a0 = test
                S = sp.sqrt(sp.diff(a0, t) ** 2 + kv)
                aD = a0 + eps * distance * S
                n = 1 + eps * distance * sp.diff(a0, t, 2) / S
                theirs = sp.diag(-n ** 2, 1, aD ** 2 / (1 - kv * chi ** 2), aD ** 2 * chi ** 2,
                                 aD ** 2 * chi ** 2 * sp.sin(X[3]) ** 2)
                at = {P["k"]: kv, t: sp.Rational(3, 10), y: eps * distance, chi: sp.Rational(1, 3),
                      X[3]: sp.Rational(4, 7)}
                self.same(g, theirs, at, f"Gaussian normal at k = {kv}, eps = {eps}")

    def test_global_is_the_empty_wall_of_the_branch_plus_one(self):
        """Deffayet (bulkmet) at k = 1, eps = +1 and a_0 = r_c cosh(ct/r_c), with chi = sin psi, both sides."""
        sp = self.sp
        _, g, X, P = self.published("global")
        rc = P["r_c"]
        warp = (1 + sp.Abs(X[1]) / rc) ** 2
        theirs = sp.diag(-warp, 1, warp * rc ** 2 * sp.cosh(X[0] / rc) ** 2,
                         warp * rc ** 2 * sp.cosh(X[0] / rc) ** 2 * sp.sin(X[2]) ** 2,
                         warp * rc ** 2 * sp.cosh(X[0] / rc) ** 2 * sp.sin(X[2]) ** 2 * sp.sin(X[3]) ** 2)
        for side in (sp.Rational(2, 5), sp.Rational(-3, 4)):
            self.same(g, theirs, {rc: sp.Rational(7, 5), X[0]: sp.Rational(1, 3), X[1]: side, X[2]: sp.Rational(3, 5),
                                  X[3]: sp.Rational(4, 7)}, f"global at y = {side}")

    def test_planar_is_the_same_wall_in_flat_slices(self):
        sp = self.sp
        _, g, X, P = self.published("planar")
        rc = P["r_c"]
        warp = (1 + sp.Abs(X[4]) / rc) ** 2
        space = warp * sp.exp(2 * X[0] / rc)
        theirs = sp.diag(-warp, space, space, space, 1)
        for side in (sp.Rational(2, 5), sp.Rational(-3, 4)):
            self.same(g, theirs, {rc: sp.Rational(5, 8), X[0]: sp.Rational(1, 3), X[4]: side}, f"planar at y = {side}")

    def test_inertial_is_flat_space_outside_the_hyperboloid(self):
        sp = self.sp
        _, g, X, _ = self.published("inertial")
        theirs = sp.diag(-1, 1, X[1] ** 2, X[1] ** 2 * sp.sin(X[2]) ** 2, X[1] ** 2 * sp.sin(X[2]) ** 2 * sp.sin(X[3]) ** 2)
        self.same(g, theirs, {X[1]: sp.Rational(9, 4), X[2]: sp.Rational(3, 5), X[3]: sp.Rational(4, 7)}, "inertial")
        self.assertIn("R \\in [\\sqrt{r_c^2 + c^2T^2}, \\infty)", chart("inertial")["domains"])

    def test_every_published_chart_passes_its_checks(self):
        """The junction condition, flat space along the hyperboloid and the wall's delta, on the charts as published."""
        import chart_printer as cp
        import dgp_wall as dw
        for system in dw.CHARTS:
            spec = dw.chart_spec(system)
            published = chart(system)
            self.assertEqual(published["line_element"], spec["system"]["line_element"], system)
            self.assertEqual(published["domains"], spec["system"]["domains"], system)
            built = cp.Chart(spec["system"]["coords"], spec["system"]["parameters"], spec["chart_line_element"],
                             spec.get("printer"), spec.get("pretty"), spec.get("time"))
            dw.chart_check(built, system)


class History(unittest.TestCase):
    """The History's numbers against the Friedmann equation of the two branches."""

    def test_the_friedmann_equation(self):
        self.assertIn(FRIEDMANN, METRIC["history"])
        described = {x["symbol"]: x["description"] for x in chart("gaussian_normal")["parameters"]}
        self.assertIn(FRIEDMANN, described["a = a(t)"])
        self.assertIn("G_{ab} - 8\\pi GT_{ab}/c^4 = (K_{ab} - h_{ab}K)/r_c", METRIC["convention"])

    def test_the_rate_of_the_empty_wall(self):
        """With no matter the branch eps = +1 at k = 0 has H^2 = (c/r_c) H, so H = c/r_c, and the hyperboloid has radius r_c."""
        rc, c = 1.7, 2.3
        H = c / rc
        self.assertAlmostEqual(H ** 2 - (c / rc) * H, 0.0, places=12)
        self.assertIn("H = c/r_c", METRIC["history"])
        self.assertIn("-c^2T^2 + X_1^2 + X_2^2 + X_3^2 + X_4^2 = r_c^2", METRIC["history"])

    def test_the_crossover_length_of_a_flat_universe(self):
        """Deffayet, Dvali and Gabadadze: today H_0^2 - H_0 c/r_c = Omega_m H_0^2, so r_c = c/((1 - Omega_m) H_0)."""
        self.assertEqual(round(1 / (1 - 0.3), 1), 1.4)
        self.assertIn("r_c = 1.4\\,c/H_0", METRIC["history"])

    def test_the_tension_of_the_inside(self):
        """At a = r_c cosh(ct/r_c), k = 1, H^2 + c^2/a^2 = c^2/r_c^2 at every t, so the branch eps = -1 has
        8 pi G rho/3 = 2c^2/r_c^2 and the tension rho c^2 = 3c^4/(4 pi G r_c^2) the drawings declare."""
        rc, c, G = 1.3, 0.7, 0.9
        for t in (-1.0, 0.0, 0.6):
            a = rc * math.cosh(c * t / rc)
            H = (c / rc) * math.tanh(c * t / rc)
            root = math.sqrt(H ** 2 + c ** 2 / a ** 2)
            self.assertAlmostEqual(root, c / rc, places=12)
            rho = 3 * (root ** 2 + root * c / rc) / (8 * math.pi * G)
            self.assertAlmostEqual(rho * c ** 2, 3 * c ** 4 / (4 * math.pi * G * rc ** 2), places=12)
        import dgp_wall_drawings as dd
        self.assertIn("3c^4/(4\\pi Gr_c^2)", dd.GN_INPUT)


class Drawings(unittest.TestCase):
    """The drawings are of the empty wall of the branch eps = +1 and, on the Gaussian normal plane, of the
    inside of the same hyperboloid on the branch eps = -1."""

    def test_one_member(self):
        import dgp_wall_drawings as dd
        diagrams = json.loads((DATA / "diagrams" / "dgp_wall.json").read_text(encoding="utf-8"))
        gn = diagrams["systems"]["gaussian_normal"][0]
        self.assertIn("a = r_c\\cosh(ct/r_c)", gn["input"])
        self.assertEqual(dd.MOMENTS, (-1.0, -0.5, 0.0, 0.5, 1.0))
        self.assertGreaterEqual(dd.GN_BOX[0], -1.1)
        # The planar view ends where the embedding does, |y| = 1.5 r_c.
        self.assertEqual(dd.PLANAR_BOX[:2], (-dd.REACH, dd.REACH))

    def test_three_kinds_are_drawn_and_the_curvature_figure_is_not(self):
        for kind in ("diagrams", "conformal", "embedding"):
            self.assertTrue((DATA / kind / "dgp_wall.json").exists(), kind)
        self.assertFalse((DATA / "curvature" / "dgp_wall.json").exists())
        notes = (ROOT / "_tools" / "derivations" / "dgp_wall.md").read_text(encoding="utf-8")
        self.assertIn("- Curvature figure: not drawn.", notes)
        for c in METRIC["coordinates"]:
            off = c["kretschmann"].replace(" \\;\\text{for}\\; y \\neq 0", "")
            self.assertEqual(off, "K = 0", c["id"])


if __name__ == "__main__":
    unittest.main()
