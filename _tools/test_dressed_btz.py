"""The conformally dressed black hole of three dimensions, Martinez and Zanelli's, gr-qc/9604021:
.venv.noindex/bin/python -m unittest discover -s _tools

What its charts, texts and drawings state, held on the published files. Each chart's published
metric is read through the checker's own Reader and held at points to Martinez and Zanelli's (Fp),
its Ricci and Kretschmann scalars to their (R) and (KR), and at B = 0 to the published BTZ chart of
zero mass; the History to stating their formulas; the references to first citation and the
bibliography; the related spacetimes to answering; and the drawings at B = l/2 to the surfaces,
curves and moments they draw. The tests that read the mathematics need sympy and are skipped where
it is absent.
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
ME = "dressed_btz"
METRIC = json.loads((DATA / "metrics" / f"{ME}.json").read_text(encoding="utf-8"))
SYSTEMS = ("static", "eddington_finkelstein_ingoing")


def chart(system, metric=METRIC):
    return next(c for c in metric["coordinates"] if c["id"] == system)


def drawing(kind):
    return json.loads((DATA / kind / f"{ME}.json").read_text(encoding="utf-8"))


def cited(text):
    keys = []
    for group in re.findall(r"\[([a-z][^\]]*)\]", text):
        for key in (k.strip() for k in group.split(",")):
            if key not in keys:
                keys.append(key)
    return keys


@unittest.skipUnless(HAS_SYMPY, "needs sympy")
class Charts(unittest.TestCase):
    """Each chart at points against Martinez and Zanelli's statements."""

    @classmethod
    def setUpClass(cls):
        import sympy as sp
        import verify_metrics as vm
        cls.sp, cls.vm = sp, vm

    def metric(self, system, metric=METRIC):
        c = chart(system, metric)
        reader = self.vm.Reader(c["coords"], [p["symbol"] for p in c["parameters"]], ())
        g = self.vm.metric_from_line_element(reader, c["line_element"].replace("c^2", "").replace("2c\\,", "2\\,"),
                                             c["coords"])
        return c, reader, g

    POINTS = ((1.4, 0.3, 1.7), (1.0, 0.5, 0.3), (1.5, 2 / 7, 5.0))

    def test_the_lapse_is_martinez_and_zanelli_s(self):
        for system in SYSTEMS:
            c, R, g = self.metric(system)
            l, B, r = R.parameters["l"], R.parameters["B"], R.symbol["r"]
            for lv, bv, rv in self.POINTS:
                at = {l: lv, B: bv, r: rv}
                want = (rv * rv - 3 * bv * bv - 2 * bv ** 3 / rv) / lv ** 2
                with self.subTest(system=system, at=(lv, bv, rv)):
                    self.assertAlmostEqual(float(-g[0, 0].subs(at)), want, places=12)
                    if system == "static":
                        self.assertAlmostEqual(float(g[1, 1].subs(at)), 1 / want, places=10)
                    else:
                        self.assertEqual(float(g[0, 1].subs(at)), 1.0)

    def test_the_scalars_are_their_r_and_kr(self):
        sp = self.sp
        for system in SYSTEMS:
            c, R, _ = self.metric(system)
            l, B, r = R.parameters["l"], R.parameters["B"], R.symbol["r"]
            ricci = R(c["ricci_scalar"].split("=", 1)[1])
            kretschmann = R(c["kretschmann"].split("=", 1)[1])
            for lv, bv, rv in self.POINTS:
                at = {l: sp.nsimplify(lv), B: sp.nsimplify(bv), r: sp.nsimplify(rv)}
                with self.subTest(system=system, at=(lv, bv, rv)):
                    self.assertAlmostEqual(float(ricci.subs(at)), -6 / lv ** 2, places=12)
                    self.assertAlmostEqual(float(kretschmann.subs(at)),
                                           12 * (rv ** 6 + 2 * bv ** 6) / (lv ** 4 * rv ** 6), places=9)

    def test_at_b_zero_each_chart_is_the_btz_hole_of_zero_mass(self):
        btz = json.loads((DATA / "metrics" / "btz.json").read_text(encoding="utf-8"))
        for system, theirs in (("static", "stationary"), ("eddington_finkelstein_ingoing", "eddington_finkelstein_ingoing")):
            c, R, g = self.metric(system)
            o = chart(theirs, btz)
            other = self.vm.Reader(o["coords"], [p["symbol"] for p in o["parameters"]], ())
            values = {tuple(e["indices"]): e["value"] for e in o["metric_components"]}
            r = R.symbol["r"]
            for i, a in enumerate(c["coords"]):
                for j, b in enumerate(c["coords"]):
                    mine = g[i, j].subs(R.parameters["B"], 0).subs(R.parameters["l"], 1).subs(r, 2)
                    there = other(values.get((o["coords"][i], o["coords"][j]), "0")).subs(
                        {other.parameters["ell"]: 1, other.parameters["M"]: 0, other.parameters["J"]: 0,
                         other.symbol["r"]: 2})
                    with self.subTest(system=system, slot=(a, b)):
                        self.assertAlmostEqual(float(mine), float(there), places=12)


class Texts(unittest.TestCase):
    """The History states the source's formulas, and every citation is held."""

    def test_history_states_the_solution(self):
        history = METRIC["history"]
        for formula in ("g^{rr} = (r + B)^2(r - 2B)/(rl^2)", "\\Psi = \\sqrt{8B/(\\kappa(r + B))}",
                        "12(r^6 + 2B^6)/(l^4r^6)", "9r_+/(16\\pi l^2)", "r_+ = 2B", "\\pi r_+/3"):
            with self.subTest(formula=formula):
                self.assertIn(formula, history)
        self.assertEqual([len(p) > 0 for p in history.split("¶")], [True] * 6)

    def test_references_follow_first_citation_and_are_in_the_bibliography(self):
        keys = cited(METRIC["history"])
        self.assertEqual(METRIC["references"][:len(keys)], keys)
        for entry in METRIC["related"]:
            for key in cited(entry["text"]):
                self.assertIn(key, METRIC["references"])
        bib = (ROOT / "assets" / "data" / "references.bib").read_text(encoding="utf-8")
        for key in METRIC["references"]:
            self.assertEqual(len(re.findall(r"^@\w+\{" + re.escape(key) + ",", bib, re.M)), 1, key)

    def test_every_related_spacetime_answers(self):
        self.assertEqual(len(METRIC["related"]), 7)
        for entry in METRIC["related"]:
            other = json.loads((DATA / "metrics" / f"{entry['id']}.json").read_text(encoding="utf-8"))
            answer = [e for e in other["related"] if e["id"] == ME]
            with self.subTest(other=entry["id"]):
                self.assertEqual(len(answer), 1)
                self.assertEqual(answer[0]["kind"], "family")
                for key in cited(answer[0]["text"]):
                    self.assertIn(key, other["references"])


class Drawings(unittest.TestCase):
    """The drawings at l = 1 and B = l/2, where the horizon is r_+ = l."""

    LEVEL = 1.389228559129   # the largest root of 4r^3 - 7r - 1, where g_rr = 1

    def test_embedding_runs_through_the_throat_into_both_exteriors(self):
        (view,) = drawing("embedding")["views"]
        self.assertEqual(view["id"], "throat")
        pieces = {p["id"]: p for p in view["surfaces"][0]["pieces"]}
        self.assertEqual(sorted(pieces), ["exterior", "exterior_minkowski", "other_exterior", "other_exterior_minkowski"])
        for name, piece in pieces.items():
            with self.subTest(piece=name):
                self.assertEqual(piece.get("space") == "minkowski", name.endswith("minkowski"))
                for r, rho, _ in piece["points"]:
                    self.assertAlmostEqual(rho, r, places=6)
        self.assertAlmostEqual(pieces["exterior"]["points"][0][0], 1.0, places=12)
        self.assertAlmostEqual(pieces["exterior"]["points"][-1][0], self.LEVEL, places=9)
        self.assertAlmostEqual(4 * self.LEVEL ** 3 - 7 * self.LEVEL - 1, 0.0, places=9)

    def test_conformal_views_put_the_boundary_straight_and_bow_the_singularity(self):
        views = drawing("conformal")["views"]
        self.assertEqual([v["id"] for v in views], ["static", "ingoing"])
        middle = 2 * math.atan(math.sqrt(2 * math.exp(-1.5)))
        for v in views:
            with self.subTest(view=v["id"]):
                boundary = [layer for layer in v["layers"] if layer["class"] == "boundary"]
                self.assertTrue(all(abs(abs(x) - math.pi / 2) < 1e-3 for layer in boundary for x, _ in layer["points"]))
                top = [layer for layer in v["layers"] if layer["class"] == "singular"][0]["points"]
                self.assertEqual([layer["kind"] for layer in v["layers"] if layer["class"] == "singular"], ["zig", "zig"])
                centre = min(top, key=lambda xt: abs(xt[0]))
                self.assertAlmostEqual(abs(centre[1]), middle, places=2)
                self.assertEqual([s["view"] for s in v["slices"]], ["throat"])

    def test_curvature_figure_draws_k_alone(self):
        (view,) = drawing("curvature")["views"]
        self.assertEqual([p["id"] for p in view["panels"]], ["K"])
        self.assertEqual({c["class"] for c in view["curves"]}, {"K"})
        self.assertEqual(view["marks"], [{"class": "horizon", "x": 1.0}])
        self.assertIn("$-6/l^2$ at every point", view["caption"][0])

    def test_each_spacetime_diagram_marks_the_moment_embedded(self):
        systems = drawing("diagrams")["systems"]
        self.assertEqual(sorted(systems), sorted(SYSTEMS))
        for system, box in (("static", [0, 3, -1.5, 1.5]), ("eddington_finkelstein_ingoing", [0, 3, -3, 0])):
            with self.subTest(system=system):
                (view,) = systems[system]
                self.assertEqual(view["box"], box)
                self.assertEqual([s["view"] for s in view["slices"]], ["throat"])


if __name__ == "__main__":
    unittest.main()
