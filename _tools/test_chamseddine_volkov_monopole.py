"""The Chamseddine-Volkov monopole of N = 4 gauged supergravity:
.venv.noindex/bin/python -m unittest discover -s _tools

What its chart, texts and drawings state, held on the published files. The published metric is read
through the checker's own Reader at points and held to Chamseddine and Volkov's eq. (13), and its
Einstein tensor to the field equations of their action with the fields of their eq. (14); the
Schwarzschild function N = rho^2/R^2 of their companion paper to its limits at the centre and far
away; the stored scalars to their values at the centre; the History to stating the line element
and the fields the chart has; the related spacetimes to answering; and the drawings to the monopole
they draw at a = 1. The tests that read the mathematics need sympy and are skipped where it is absent.
"""
import importlib.util
import json
import math
import re
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
HAS_SYMPY = importlib.util.find_spec("sympy") is not None
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))
ID = "chamseddine_volkov_monopole"
METRIC = json.loads((DATA / "metrics" / f"{ID}.json").read_text(encoding="utf-8"))
SYSTEM = "chamseddine_volkov"
CHART = next(c for c in METRIC["coordinates"] if c["id"] == SYSTEM)


def drawing(kind):
    return json.loads((DATA / kind / f"{ID}.json").read_text(encoding="utf-8"))


@unittest.skipUnless(HAS_SYMPY, "needs sympy")
class Chart(unittest.TestCase):
    """The chart at points against its sources."""

    @classmethod
    def setUpClass(cls):
        import sympy as sp
        import verify_metrics as vm
        import chamseddine_volkov_monopole as cv
        cls.sp, cls.vm, cls.cv = sp, vm, cv
        cls.reader = vm.Reader(CHART["coords"], [p["symbol"] for p in CHART["parameters"]], (),
                               held=vm.HELD[(ID, SYSTEM)])
        cls.t, cls.rho, cls.theta, cls.varphi = (cls.reader.symbol[c] for c in CHART["coords"])
        cls.a = cls.reader.parameters["a"]
        cls.points = [{cls.rho: sp.Rational(7, 10), cls.a: sp.Rational(13, 10), cls.theta: sp.Rational(1, 3)},
                      {cls.rho: sp.Rational(23, 10), cls.a: sp.Rational(1, 2), cls.theta: sp.Rational(6, 5)},
                      {cls.rho: 5, cls.a: 2, cls.theta: sp.Rational(5, 7)}]

    def written(self, text):
        """A published value with R written out as its definition."""
        return self.sp.sympify(self.reader(text)).xreplace(self.reader.held)

    def ours(self):
        sp = self.sp
        _, _, R = self.cv.fields(self.rho, self.a)
        return R, self.a ** 2 * sp.sinh(self.rho) / R * sp.diag(-1, 1, R ** 2, R ** 2 * sp.sin(self.theta) ** 2)

    def test_the_chart_is_their_eq_13(self):
        sp = self.sp
        g = self.vm.metric_from_line_element(self.reader, CHART["line_element"], CHART["coords"]).xreplace(self.reader.held)
        _, want = self.ours()
        for at in self.points:
            for i in range(4):
                with self.subTest(slot=i, at=at):
                    self.assertLess(abs(float(sp.N((g[i, i] - want[i, i]).subs(at), 40))), 1e-30)

    def test_the_published_einstein_tensor_is_twice_the_stress_of_their_action(self):
        sp = self.sp
        _, g = self.ours()
        twice = self.cv.stress(g, self.rho, self.theta, self.a)
        published = {tuple(c["indices"]): c["value"] for c in CHART["einstein_tensor"]["variants"]["ll"]["nonzero"]}
        self.assertEqual(sorted(published), sorted((c, c) for c in CHART["coords"]))
        for i, name in enumerate(CHART["coords"]):
            value = self.written(published[(name, name)])
            for at in self.points:
                with self.subTest(slot=name, at=at):
                    self.assertLess(abs(float(sp.N((value - twice[i]).subs(at), 40))), 1e-28)

    def test_the_stored_scalars_at_the_centre(self):
        """R^mu_mu = -2/a^2 and K = 20/(27 a^4) at the centre, read at rho = 1e-4 in 120 digits, since
        the stored values hold R unreduced and their terms cancel through more than sixty there."""
        import mpmath as mp
        sp = self.sp
        for field, want in (("ricci_scalar", -2), ("kretschmann", sp.Rational(20, 27))):
            f = sp.lambdify((self.rho, self.a), self.written(CHART[field].split("=", 1)[1]), "mpmath")
            for a in (1, 3):
                with self.subTest(field=field, a=a), mp.workdps(120):
                    value = f(mp.mpf("1e-4"), a) * a ** (2 if field == "ricci_scalar" else 4)
                    self.assertLess(abs(float(value - want)), 1e-7)


class Schwarzschild(unittest.TestCase):
    """The Schwarzschild function N = rho^2/R^2 of hep-th/9711181 eq. (5.9) against the areal radius
    r = a sqrt(R sinh rho), at a = 1: N = 1 + r^2/9 + O(r^4) near the centre and N/ln r -> 1 far away."""

    @staticmethod
    def R2(rho):
        import mpmath as mp
        return 2 * rho * mp.coth(rho) - rho ** 2 / mp.sinh(rho) ** 2 - 1

    def test_near_the_centre(self):
        import mpmath as mp
        mp.mp.dps = 60
        for rho in (mp.mpf("1e-2"), mp.mpf("1e-3")):
            N = rho ** 2 / self.R2(rho)
            r2 = mp.sqrt(self.R2(rho)) * mp.sinh(rho)
            with self.subTest(rho=rho):
                self.assertLess(abs(float((N - 1) / (r2 / 9) - 1)), 2 * float(rho) ** 2)

    def test_far_away_N_grows_as_ln_r(self):
        import mpmath as mp
        mp.mp.dps = 60
        ratios = []
        for rho in (mp.mpf(200), mp.mpf(2000), mp.mpf(20000)):
            N = rho ** 2 / self.R2(rho)
            ln_r = mp.log(mp.sqrt(self.R2(rho)) * mp.sinh(rho)) / 2
            ratios.append(float(N / ln_r))
        self.assertLess(abs(ratios[0] - 1), 0.02)
        self.assertTrue(abs(ratios[2] - 1) < abs(ratios[1] - 1) < abs(ratios[0] - 1))


class Texts(unittest.TestCase):
    """The History states the chart's mathematics, and every citation is held."""

    def test_history_states_the_line_element_and_the_fields(self):
        history = METRIC["history"]
        self.assertIn("ds^2 = a^2 (\\sinh\\rho / R)(-dt^2 + d\\rho^2 + R^2 d\\Omega^2)", history)
        self.assertIn("R^2 = 2\\rho\\coth\\rho - \\rho^2/\\sinh^2\\rho - 1", history)
        self.assertIn("w = \\rho/\\sinh\\rho", history)
        self.assertIn("r = a\\sqrt{R\\sinh\\rho}", history)
        self.assertIn("w = \\rho/\\sinh\\rho", CHART["convention"])
        self.assertIn("e^{2\\phi} = a^2\\sinh\\rho/(2R)", CHART["convention"])

    def test_references_follow_first_citation_and_are_in_the_bibliography(self):
        cited = []
        for group in re.findall(r"\[([a-z][^\]]*)\]", METRIC["history"]):
            for key in (k.strip() for k in group.split(",")):
                if key not in cited:
                    cited.append(key)
        self.assertEqual(METRIC["references"][:len(cited)], cited)
        bib = (ROOT / "assets" / "data" / "references.bib").read_text(encoding="utf-8")
        for key in METRIC["references"]:
            self.assertEqual(len(re.findall(r"^@\w+\{" + re.escape(key) + ",", bib, re.M)), 1, key)

    def test_every_related_spacetime_answers(self):
        for entry in METRIC["related"]:
            other = json.loads((DATA / "metrics" / f"{entry['id']}.json").read_text(encoding="utf-8"))
            answer = [e for e in other["related"] if e["id"] == ID]
            with self.subTest(other=entry["id"]):
                self.assertEqual(len(answer), 1)
                self.assertEqual(answer[0]["kind"], "family")
                for group in re.findall(r"\[([a-z][^\]]*)\]", answer[0]["text"]):
                    for key in (k.strip() for k in group.split(",")):
                        self.assertIn(key, other["references"])


class Drawings(unittest.TestCase):
    """The drawings of the monopole at a = 1."""

    def test_the_equator_is_drawn_in_minkowski_space_with_its_areal_radius(self):
        (view,) = drawing("embedding")["views"]
        self.assertEqual(view["id"], "equator")
        (piece,) = view["surfaces"][0]["pieces"]
        self.assertEqual(piece["space"], "minkowski")
        heights = [z for _, _, z in piece["points"]]
        self.assertEqual(heights, sorted(heights))
        self.assertEqual(piece["points"][0][:2], [0.0, 0.0])
        for rho, r, _ in piece["points"][1:]:
            R = math.sqrt(2 * rho / math.tanh(rho) - rho ** 2 / math.sinh(rho) ** 2 - 1) if rho > 1e-3 else rho
            self.assertAlmostEqual(r, math.sqrt(R * math.sinh(rho)), places=6)

    def test_curvature_figure_draws_only_the_classes_the_page_styles(self):
        (view,) = drawing("curvature")["views"]
        self.assertEqual([p["id"] for p in view["panels"]], ["K", "R"])
        self.assertLessEqual({c["class"] for c in view["curves"]}, {"K", "R"})
        self.assertEqual({c["params"]["a"] for c in view["curves"]}, {1.0})
        self.assertEqual(len(view["caption"]), 2)

    def test_the_conformal_diagram_is_minkowski_s_triangle_with_a_regular_centre(self):
        (view,) = drawing("conformal")["views"]
        classes = {layer["class"] for layer in view["layers"]}
        self.assertIn("centre", classes)
        self.assertNotIn("singular", classes)
        (centre,) = [layer for layer in view["layers"] if layer["class"] == "centre"]
        self.assertTrue(all(abs(x) < 1e-9 for x, _ in centre["points"]))

    def test_the_spacetime_diagram_marks_the_embedded_moment(self):
        (view,) = drawing("diagrams")["systems"][SYSTEM]
        self.assertEqual(view["box"], [0, 4, -2, 2])
        self.assertEqual([s["view"] for s in view["slices"]], ["equator"])


if __name__ == "__main__":
    unittest.main()
