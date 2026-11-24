"""Krasnikov's singularity-free time machine that respects the weak energy condition:
.venv.noindex/bin/python -m unittest discover -s _tools

What its charts, texts and drawings state, held on the published files. The published metrics are
read through the checker's own Reader and held at points, with trial functions put in for the free
f, a and rho, to Krasnikov's eqs. (2), (5), (7), (9) and (10a); the published Einstein tensor of the
time machine to the energy density 1/R^2 and the condition K <= R^-2 of his eq. (6), and that of the
universe off the square to obeying the weak energy condition under his eq. (11) and to Friedmann's
closed universe where rho = cos chi; the History to stating the line elements the charts have; the
related spacetimes to answering; and the drawings to the example they declare. The tests that read
the mathematics need sympy and are skipped where it is absent.
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
ID = "krasnikov_time_machine"
METRIC = json.loads((DATA / "metrics" / f"{ID}.json").read_text(encoding="utf-8"))
CHARTS = {c["id"]: c for c in METRIC["coordinates"]}


def drawing(kind):
    return json.loads((DATA / kind / f"{ID}.json").read_text(encoding="utf-8"))


@unittest.skipUnless(HAS_SYMPY, "needs sympy")
class Charts(unittest.TestCase):
    """Both charts at points against Krasnikov's article."""

    @classmethod
    def setUpClass(cls):
        import sympy as sp
        import verify_metrics as vm
        cls.sp, cls.vm = sp, vm
        cls.readers = {s: vm.Reader(c["coords"], [p["symbol"] for p in c["parameters"]], ())
                       for s, c in CHARTS.items()}
        reader = cls.readers["time_machine"]
        cls.tau, cls.chi, cls.theta, cls.phi = (reader.symbol[c] for c in CHARTS["time_machine"]["coords"])

    def published(self, system, field):
        """A tensor's lowered components as sympy, read from the file."""
        chart, reader = CHARTS[system], self.readers[system]
        variants = chart[field]["variants"]
        entries = variants["ll"]["nonzero"]
        return {tuple(chart["coords"].index(i) for i in e["indices"]): self.sp.sympify(reader(e["value"]))
                for e in entries}

    def metric(self, system):
        chart, reader = CHARTS[system], self.readers[system]
        return self.vm.metric_from_line_element(reader, chart["line_element"], chart["coords"])

    def near(self, value, want, at):
        return abs(complex(self.sp.N((value - want).subs(at), 40))) < 1e-26

    def test_the_time_machine_is_his_eqs_2_5_and_7(self):
        sp = self.sp
        P = self.readers["time_machine"].parameters
        f, R = P["f"], P["R"]
        g = self.metric("time_machine")
        want = sp.diag(-f ** -2, f ** -2, R ** 2, R ** 2 * sp.sin(self.theta) ** 2)
        self.assertEqual(sp.simplify(g - want), sp.zeros(4, 4))

    def test_the_weak_energy_condition_is_his_eq_6(self):
        """The orthonormal Einstein tensor is diag(1/R^2, -1/R^2, -K, -K) with K the right side of
        eq. (6), and the Ricci scalar is 2K + 2/R^2, for a trial f with his four zeros."""
        sp = self.sp
        tau, chi, theta = self.tau, self.chi, self.theta
        P = self.readers["time_machine"].parameters
        f, R = P["f"], P["R"]
        h = sp.Rational(1, 5)
        trial = (sp.sin(sp.pi * (tau + h) / (2 * h)) ** 2 + (chi ** 2 - h ** 2) ** 2) / (12 * R)
        K = trial * (sp.diff(trial, chi, 2) - sp.diff(trial, tau, 2)) + sp.diff(trial, tau) ** 2 - sp.diff(trial, chi) ** 2
        G = self.published("time_machine", "einstein_tensor")
        self.assertEqual(sorted(G), [(i, i) for i in range(4)])
        g = [-trial ** -2, trial ** -2, R ** 2, R ** 2 * sp.sin(theta) ** 2]
        want = [1 / R ** 2, -1 / R ** 2, -K, -K]
        scalar = self.sp.sympify(self.readers["time_machine"](CHARTS["time_machine"]["ricci_scalar"].split("=", 1)[1]))
        for at in ({tau: sp.Rational(1, 7), chi: -sp.Rational(2, 9), R: 1, theta: sp.Rational(1, 3)},
                   {tau: -sp.Rational(3, 10), chi: sp.Rational(1, 4), R: sp.Rational(4, 5), theta: sp.Rational(6, 5)}):
            for i in range(4):
                with self.subTest(slot=i, at=at):
                    value = (G[(i, i)].subs(f, trial).doit() / abs(g[i]))
                    self.assertTrue(self.near(value, want[i], at))
            self.assertTrue(self.near(scalar.subs(f, trial).doit(), 2 * K + 2 / R ** 2, at))

    def test_the_universe_is_his_eq_9_and_on_the_square_eq_10a(self):
        sp = self.sp
        reader = self.readers["universe"]
        tau, chi, theta = (reader.symbol[c] for c in (r"\tau", r"\chi", r"\theta"))
        a, rho = reader.parameters["a"], reader.parameters["rho"]
        g = self.metric("universe")
        self.assertEqual(sp.simplify(g - a ** 2 * sp.diag(-1, 1, rho ** 2, rho ** 2 * sp.sin(theta) ** 2)), sp.zeros(4, 4))
        f, R = sp.Function("f")(tau, chi), sp.Symbol("R", positive=True)
        square = g.subs({a: 1 / f, rho: f * R})
        machine = sp.diag(-f ** -2, f ** -2, R ** 2, R ** 2 * sp.sin(theta) ** 2)
        self.assertEqual(sp.simplify(square - machine), sp.zeros(4, 4))

    def test_off_the_square_his_eq_11_gives_the_weak_energy_condition_and_friedmann_where_rho_is_cos_chi(self):
        sp = self.sp
        reader = self.readers["universe"]
        tau, chi, theta = (reader.symbol[c] for c in (r"\tau", r"\chi", r"\theta"))
        a, rho = reader.parameters["a"], reader.parameters["rho"]
        G = self.published("universe", "einstein_tensor")
        scale = 3 * sp.cos(tau / 2) ** 2 * (2 + sp.tanh(tau) / 3)
        Phi = sp.log(scale)
        d1, d2 = sp.diff(Phi, tau), sp.diff(Phi, tau, 2)
        for profile in (2 * (sp.cos(chi) + sp.Rational(1, 2)) / 3, sp.cos(chi)):
            r1, r2 = sp.diff(profile, chi), sp.diff(profile, chi, 2)
            g = [-scale ** 2, scale ** 2, scale ** 2 * profile ** 2, scale ** 2 * profile ** 2 * sp.sin(theta) ** 2]
            hat = {k: v.subs({a: scale, rho: profile}).doit() / sp.sqrt(abs(g[k[0]] * g[k[1]])) for k, v in G.items()}
            energy = 3 * d1 ** 2 + (1 - r1 ** 2) / profile ** 2 - 2 * r2 / profile
            along = 2 * d1 ** 2 - 2 * d2 - 2 * r2 / profile
            across = 2 * d1 ** 2 - 2 * d2 + (1 - r1 ** 2) / profile ** 2 - r2 / profile
            for at in ({tau: sp.Rational(1, 4), chi: sp.Rational(7, 10), theta: sp.Rational(1, 3)},
                       {tau: -sp.Rational(9, 10), chi: -sp.Rational(13, 10), theta: sp.Rational(6, 5)}):
                with self.subTest(profile=str(profile), at=at):
                    self.assertTrue(self.near(hat.get((0, 1), 0), 0, at))
                    self.assertTrue(self.near(scale ** 2 * hat[(0, 0)], energy, at))
                    self.assertTrue(self.near(scale ** 2 * (hat[(0, 0)] + hat[(1, 1)]), along, at))
                    self.assertTrue(self.near(scale ** 2 * (hat[(0, 0)] + hat[(2, 2)]), across, at))
                    for condition in (-d2, -r2 / profile, 1 - r1 ** 2, energy, along, across):
                        self.assertGreaterEqual(float(sp.N(condition.subs(at), 30)), 0)
                    if profile == sp.cos(chi):
                        self.assertTrue(self.near(hat[(1, 1)], hat[(2, 2)], at))
                        self.assertTrue(self.near(scale ** 2 * hat[(0, 0)], 3 * (d1 ** 2 + 1), at))


class Texts(unittest.TestCase):
    """The History states the charts' mathematics, and every citation is held."""

    def test_history_states_both_line_elements_and_the_condition(self):
        history = METRIC["history"]
        self.assertIn("ds^2 = f^{-2}(-d\\tau^2 + d\\chi^2)", history)
        self.assertIn("R^{-2} \\ge f(f_{,\\chi\\chi} - f_{,\\tau\\tau}) + f_{,\\tau}^2 - f_{,\\chi}^2", history)
        self.assertIn("ds^2 = a^2[-d\\tau^2 + d\\chi^2 + \\rho^2(d\\theta^2 + \\sin^2\\theta\\, d\\phi^2)]", history)
        self.assertIn("$a = 1/f$ and $\\rho = fR$", history)
        self.assertIn("f(\\tau, \\chi) = f(\\tau - 2h, \\chi)", history)
        self.assertIn("$\\rho = \\cos\\chi$", history)

    def test_no_text_of_its_own_says_convex_or_claims_the_dominant_energy_condition(self):
        own = [METRIC["history"], METRIC["description"], METRIC["convention"]]
        for chart in METRIC["coordinates"]:
            own += [chart["convention"]] + [p["description"] for p in chart["parameters"]]
        for text in own:
            self.assertNotIn("convex", text.lower())
            self.assertNotIn("dominant", text.lower())

    def test_references_follow_first_citation_and_are_in_the_bibliography(self):
        cited = []
        texts = [METRIC["history"]] + [e["text"] for e in METRIC["related"]]
        for text in texts:
            for group in re.findall(r"\[([a-z][^\]]*)\]", text):
                for key in (k.strip() for k in group.split(",")):
                    if key not in cited:
                        cited.append(key)
        self.assertEqual(METRIC["references"], cited)
        bib = (ROOT / "assets" / "data" / "references.bib").read_text(encoding="utf-8")
        for key in METRIC["references"]:
            self.assertEqual(len(re.findall(r"^@\w+\{" + re.escape(key) + ",", bib, re.M)), 1, key)

    def test_every_related_spacetime_answers_with_its_kind(self):
        answers = {"piece": "composite", "family": "family", "programme": "programme"}
        for entry in METRIC["related"]:
            other = json.loads((DATA / "metrics" / f"{entry['id']}.json").read_text(encoding="utf-8"))
            answer = [e for e in other["related"] if e["id"] == ID]
            with self.subTest(other=entry["id"]):
                self.assertEqual(len(answer), 1)
                self.assertEqual(answer[0]["kind"], answers[entry["kind"]])
                for group in re.findall(r"\[([a-z][^\]]*)\]", answer[0]["text"]):
                    for key in (k.strip() for k in group.split(",")):
                        self.assertIn(key, other["references"])


class Drawings(unittest.TestCase):
    """The drawings at m = 1/2, h = 1/5 and n = 1, with R = 4/5 and a = 1 at tau = 0."""

    def test_each_plane_marks_the_cuts_the_four_points_and_the_moment(self):
        systems = drawing("diagrams")["systems"]
        self.assertEqual(sorted(systems), ["time_machine", "universe"])
        for system, box in (("time_machine", [-0.5, 0.5, -0.5, 0.5]), ("universe", [-math.pi / 2, math.pi / 2, -1, 1])):
            (view,) = systems[system]
            with self.subTest(system=system):
                for want, got in zip(box, view["box"]):
                    self.assertAlmostEqual(want, got, places=9)
                kinds = [mk["kind"] for mk in view["markers"]]
                self.assertEqual(kinds.count("removed"), 4)
                self.assertEqual(kinds.count("surface"), 2)
                self.assertEqual([s["view"] for s in view["slices"]], ["moment"])
                self.assertIn("any", view["input"])

    def test_the_moment_is_a_cylinder_of_radius_r_between_two_caps_of_the_unit_sphere(self):
        (view,) = drawing("embedding")["views"]
        self.assertEqual(view["id"], "moment")
        (piece,) = view["surfaces"][0]["pieces"]
        heights = [z for _, _, z in piece["points"]]
        self.assertEqual(heights, sorted(heights))
        for chi, r, _ in piece["points"]:
            if abs(chi) <= 0.5:
                self.assertAlmostEqual(r, 0.8, places=6)
            elif abs(chi) >= 1:
                self.assertAlmostEqual(r, math.cos(chi), places=6)
            else:
                self.assertLessEqual(r, 0.8 + 1e-6)
                self.assertGreaterEqual(r, math.cos(1) - 1e-6)
        self.assertAlmostEqual(piece["points"][0][1], 0.0, places=6)
        self.assertAlmostEqual(piece["points"][-1][1], 0.0, places=6)

    def test_no_conformal_diagram_or_curvature_figure_is_published(self):
        self.assertFalse((DATA / "conformal" / f"{ID}.json").exists())
        self.assertFalse((DATA / "curvature" / f"{ID}.json").exists())


if __name__ == "__main__":
    unittest.main()
