"""Repulsons and massless black holes, Kallosh and Linde's charged solutions of the heterotic string:
.venv.noindex/bin/python -m unittest discover -s _tools

What its charts, texts and drawings state, held on the published files. Each chart's published
metric is read through the checker's own Reader at points and held to Kallosh and Linde's eq. (36)
or (43), and its Ricci scalar at M = 0 to Kallosh's closed form; the History to stating the line
element, the singular radius and the equilibrium radius the chart has; the related spacetimes to
answering; and the drawings to the massless member they draw: the embedding surfaces to their
radii, the curvature figure to the classes the page draws, and the diagrams to their slices.
The tests that read the mathematics need sympy and are skipped where it is absent.
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
METRIC = json.loads((DATA / "metrics" / "repulson.json").read_text(encoding="utf-8"))
SYSTEMS = ("isotropic", "string")


def chart(system):
    return next(c for c in METRIC["coordinates"] if c["id"] == system)


def drawing(kind):
    return json.loads((DATA / kind / "repulson.json").read_text(encoding="utf-8"))


@unittest.skipUnless(HAS_SYMPY, "needs sympy")
class Charts(unittest.TestCase):
    """Each chart at points against its source."""

    @classmethod
    def setUpClass(cls):
        import sympy as sp
        import verify_metrics as vm
        cls.sp, cls.vm = sp, vm

    def metric(self, system):
        """The chart's published line element read with Delta as a symbol of its own, which each
        point gives its value, 1 + 4M/r - 4q^2/r^2, so that its root is the positive one."""
        c = chart(system)
        reader = self.vm.Reader(c["coords"], ["M", "q", "\\Delta"], ())
        g = self.vm.metric_from_line_element(reader, c["line_element"].replace("c^2", ""), c["coords"])
        return c, reader, g

    def at(self, reader, values):
        """The point with M, q, r and theta given, Delta at its value there."""
        M, q, r = values["M"], values["q"], values["r"]
        point = {reader.parameters["M"]: M, reader.parameters["q"]: q, reader.symbol["r"]: r,
                 reader.symbol["\\theta"]: values["theta"], reader.parameters["Delta"]: 1 + 4 * M / r - 4 * q ** 2 / r ** 2}
        return point

    def test_each_chart_is_kallosh_and_linde_s(self):
        sp = self.sp
        for system in SYSTEMS:
            c, R, g = self.metric(system)
            t, r, th, ph = (R.symbol[name] for name in c["coords"])
            M, q = R.parameters["M"], R.parameters["q"]
            delta = 1 + 4 * M / r - 4 * q ** 2 / r ** 2
            power = sp.Rational(1, 2) if system == "isotropic" else 0
            want = sp.diag(-delta ** (power - 1), delta ** power, delta ** power * r ** 2,
                           delta ** power * r ** 2 * sp.sin(th) ** 2)
            for values in ({"M": sp.Rational(3, 10), "q": sp.Rational(7, 10), "r": sp.Rational(29, 10), "theta": sp.Rational(1, 3)},
                           {"M": 0, "q": 1, "r": sp.Rational(5, 2), "theta": sp.Rational(4, 5)}):
                at = self.at(R, values)
                for i in range(4):
                    with self.subTest(system=system, slot=i, at=values):
                        self.assertLess(abs(float(sp.N((g[i, i] - want[i, i]).subs(at), 30))), 1e-25)

    def test_ricci_scalar_of_the_massless_member_is_kallosh_s(self):
        sp = self.sp
        for system, closed in (("isotropic", lambda r: 4 * (2 + r ** 2) / (r * (r ** 2 - 4) ** sp.Rational(5, 2))),
                               ("string", lambda r: -8 * (8 + r ** 2) / (r ** 2 * (r ** 2 - 4) ** 2))):
            c, R, _ = self.metric(system)
            r = R.symbol["r"]
            published = R(c["ricci_scalar"].split("=", 1)[1])
            for x in (sp.Rational(5, 2), sp.Integer(3), sp.Rational(9, 2)):
                with self.subTest(system=system, r=x):
                    value = published.subs(self.at(R, {"M": 0, "q": 1, "r": x, "theta": sp.Rational(1, 3)}))
                    self.assertLess(abs(float(sp.N(value - closed(x), 30))), 1e-25)

    def test_delta_vanishes_on_the_singularity(self):
        for M, q in ((0.3, 0.7), (0.5, 0.2), (0.0, 1.0)):
            r0 = 2 * (math.sqrt(M * M + q * q) - M)
            self.assertAlmostEqual(1 + 4 * M / r0 - 4 * q * q / r0 ** 2, 0.0, places=12)


class Texts(unittest.TestCase):
    """The History states the chart's mathematics, and every citation is held."""

    def test_history_states_the_line_element_and_the_radii(self):
        history = METRIC["history"]
        self.assertIn("ds^2 = -\\Delta^{-1/2}dt^2 + \\Delta^{1/2}(dr^2 + r^2d\\Omega^2)", history)
        self.assertIn("\\Delta = 1 + \\frac{4M}{r} - \\frac{4q^2}{r^2}", history)
        self.assertIn("r_0 = 2(\\sqrt{M^2 + q^2} - M)", history)
        self.assertIn("r_c = 2q^2/M", history)
        self.assertIn("ds^2 = -\\dfrac{c^2dt^2}{\\sqrt{\\Delta}} + \\sqrt{\\Delta}", chart("isotropic")["line_element"])

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
            answer = [e for e in other["related"] if e["id"] == "repulson"]
            with self.subTest(other=entry["id"]):
                self.assertEqual(len(answer), 1)
                self.assertEqual(answer[0]["kind"], "family")
                for group in re.findall(r"\[([a-z][^\]]*)\]", answer[0]["text"]):
                    for key in (k.strip() for k in group.split(",")):
                        self.assertIn(key, other["references"])


class Drawings(unittest.TestCase):
    """The drawings of the massless member, M = 0 and q = 1."""

    def test_embedding_surfaces_have_their_radii(self):
        views = {v["id"]: v for v in drawing("embedding")["views"]}
        self.assertEqual(sorted(views), ["einstein", "string"])
        einstein = views["einstein"]["surfaces"][0]["pieces"][0]
        self.assertEqual(einstein["space"], "minkowski")
        heights = [z for _, _, z in einstein["points"]]
        self.assertEqual(heights, sorted(heights))
        for r, rho, _ in einstein["points"]:
            self.assertAlmostEqual(rho, math.sqrt(r) * (r * r - 4) ** 0.25, places=8)
        plane = views["string"]["surfaces"][0]["pieces"][0]
        for r, rho, z in plane["points"]:
            self.assertAlmostEqual(rho, r, places=8)
            self.assertEqual(z, 0.0)

    def test_curvature_figure_draws_only_the_classes_the_page_styles(self):
        view = drawing("curvature")["views"][0]
        self.assertEqual([p["id"] for p in view["panels"]], ["K", "R"])
        self.assertLessEqual({c["class"] for c in view["curves"]}, {"K", "R", "reference"})
        self.assertEqual({c["params"]["M"] for c in view["curves"]}, {0.0})
        self.assertEqual(len(view["caption"]), 2)

    def test_conformal_views_draw_the_singularity_as_a_timelike_line(self):
        views = drawing("conformal")["views"]
        self.assertEqual([v["id"] for v in views], list(SYSTEMS))
        for v in views:
            singular = [layer for layer in v["layers"] if layer["class"] == "singular"]
            with self.subTest(view=v["id"]):
                self.assertEqual([layer["kind"] for layer in singular], ["zig"])
                self.assertTrue(all(abs(x) < 1e-9 for x, _ in singular[0]["points"]))

    def test_each_spacetime_diagram_marks_its_own_moment(self):
        systems = drawing("diagrams")["systems"]
        for system, embedded in (("isotropic", "einstein"), ("string", "string")):
            with self.subTest(system=system):
                (view,) = systems[system]
                self.assertEqual(view["box"], [2, 6, -2, 2])
                self.assertEqual([s["view"] for s in view["slices"]], [embedded])


if __name__ == "__main__":
    unittest.main()
