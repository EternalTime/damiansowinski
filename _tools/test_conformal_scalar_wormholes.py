"""The wormholes of a conformally coupled scalar field, Barcelo and Visser's, gr-qc/9908029:
.venv.noindex/bin/python -m unittest discover -s _tools

What its charts, texts and drawings state, held on the published files. Each chart's published
metric, Christoffel symbols and Ricci tensor are read through the checker's own Reader and held at
points to Barcelo and Visser's field equations (2.2) and (2.3) with their field, the isotropic chart
at b = 0 to the published extremal hole of Majumdar and Papapetrou, and the chart in r to Omega^2 times
the published Fisher-Janis-Newman-Winicour metric and at w = 1/6 to the isotropic chart; the History to
stating their formulas; the references to first citation and the bibliography; the related
spacetimes to answering; and the drawings at mu = 1 and b = -1/4 to the surfaces, curves and moments
they draw. The tests that read the mathematics need sympy and are skipped where it is absent.
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
ME = "conformal_scalar_wormholes"
METRIC = json.loads((DATA / "metrics" / f"{ME}.json").read_text(encoding="utf-8"))
SYSTEMS = ("isotropic", "three_parameter")
HELD = {"three_parameter": ("F", "Omega", "Psi")}
RATES = {"three_parameter": {"F": {"r": "\\dfrac{4\\mu\\gamma F}{r\\left(r - 4\\mu\\right)}"},
                             "Omega": {"r": "\\dfrac{4\\mu\\sigma\\Psi}{r\\left(r - 4\\mu\\right)}"},
                             "Psi": {"r": "\\dfrac{4\\mu\\sigma\\Omega}{r\\left(r - 4\\mu\\right)}"}}}


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
    """Each chart at points against Barcelo and Visser's field equations and statements."""

    @classmethod
    def setUpClass(cls):
        import sympy as sp
        import verify_metrics as vm
        cls.sp, cls.vm = sp, vm

    def published(self, system, metric=METRIC):
        """The chart's reader and its published metric, Christoffel symbols and Ricci tensor as sympy,
        each with every held name written out."""
        sp = self.sp
        c = chart(system, metric)
        reader = self.vm.Reader(c["coords"], [p["symbol"] for p in c["parameters"]], (),
                                held=HELD.get(system, ()), rates=RATES.get(system))
        n = len(c["coords"])
        out = lambda text: reader(text).subs(reader.held).doit()  # noqa: E731
        index = {name: i for i, name in enumerate(c["coords"])}
        g = sp.zeros(n, n)
        for e in c["metric_components"]:
            g[index[e["indices"][0]], index[e["indices"][1]]] = out(e["value"])
        gamma = [[[0] * n for _ in range(n)] for _ in range(n)]
        for e in c["christoffel"]["variants"]["ull"]["nonzero"]:
            a, b, d = (index[i] for i in e["indices"])
            gamma[a][b][d] = out(e["value"])
        ricci = sp.zeros(n, n)
        for e in c["ricci_tensor"]["variants"]["ll"]["nonzero"]:
            ricci[index[e["indices"][0]], index[e["indices"][1]]] = out(e["value"])
        return c, reader, g, gamma, ricci

    def field_equations(self, system, phi, points):
        """Barcelo and Visser's (2.2) in every slot and (2.3), in units kappa = 1, at each point."""
        sp = self.sp
        c, reader, g, gamma, ricci = self.published(system)
        X = [reader.symbol[x] for x in c["coords"]]
        n = len(X)
        field = phi(reader)
        d = [sp.diff(field, x) for x in X]
        hessian = sp.Matrix(n, n, lambda a, b: sp.diff(field, X[a], X[b])
                            - sum(gamma[k][a][b] * d[k] for k in range(n)))
        ginv = g.inv()
        grad2 = sum(ginv[a, b] * d[a] * d[b] for a in range(n) for b in range(n))
        for at in points:
            values = {**at(reader), X[2]: sp.Rational(7, 10)}
            with self.subTest(system=system, at=values):
                for a in range(n):
                    for b in range(a, n):
                        source = sp.Rational(2, 3) * d[a] * d[b] - g[a, b] * grad2 / 6 - field * hessian[a, b] / 3
                        gap = ((1 - field ** 2 / 6) * ricci[a, b] - source).subs(values)
                        self.assertLess(abs(complex(sp.N(gap, 40))), 1e-25, (a, b))
                box = sum(ginv[a, b] * hessian[a, b] for a in range(n) for b in range(n)).subs(values)
                self.assertLess(abs(complex(sp.N(box, 40))), 1e-25)
                self.assertEqual(c["ricci_scalar"], "R = 0")

    def test_the_isotropic_chart_solves_barcelo_and_visser_s_equations_with_their_field(self):
        sp = self.sp

        def phi(R):
            mu, b, rho = R.parameters["mu"], R.parameters["b"], R.symbol["\\rho"]
            return sp.sqrt(6) * (b * rho - mu) / (rho - b * mu)

        def at(mu, b, rho):
            return lambda R: {R.parameters["mu"]: mu, R.parameters["b"]: b, R.symbol["\\rho"]: rho}
        self.field_equations("isotropic", phi, [at(1, -sp.Rational(1, 4), sp.Rational(3, 10)),
                                                at(sp.Rational(7, 10), -2, sp.Rational(17, 10))])

    def test_the_chart_in_r_solves_barcelo_and_visser_s_equations_with_their_field(self):
        sp = self.sp

        def phi(R):
            return sp.sqrt(6) * R.held[R.parameters["Psi"]] / R.held[R.parameters["Omega"]]

        def at(mu, b, w, r):
            return lambda R: {R.parameters["mu"]: mu, R.parameters["b"]: b, R.parameters["w"]: w, R.symbol["r"]: r}
        self.field_equations("three_parameter", phi, [at(1, -sp.Rational(1, 2), sp.Rational(1, 6), 5),
                                                      at(sp.Rational(3, 2), -2, 1, 9)])

    def test_the_throat_and_the_far_field_of_the_isotropic_chart(self):
        sp = self.sp
        c, R, g, _, _ = self.published("isotropic")
        mu, b, rho = R.parameters["mu"], R.parameters["b"], R.symbol["\\rho"]
        areal = sp.sqrt(g[2, 2])
        for bv in (-sp.Rational(1, 4), -1, -3):
            throat = sp.sqrt(-bv)
            with self.subTest(b=bv):
                self.assertEqual(sp.simplify(sp.diff(areal, rho).subs({mu: 1, b: bv, rho: throat})), 0)
                self.assertEqual(sp.simplify(areal.subs({mu: 1, b: bv, rho: throat}) - (1 + throat) ** 2), 0)
        self.assertEqual(sp.simplify(sp.limit((g[0, 0] + 1) * rho, rho, sp.oo) - 2 * mu * (1 + b)), 0)

    def test_at_b_zero_the_isotropic_chart_is_one_extremal_hole_of_majumdar_and_papapetrou(self):
        _, R, g, _, _ = self.published("isotropic")
        hole = json.loads((DATA / "metrics" / "majumdar_papapetrou.json").read_text(encoding="utf-8"))
        o = chart("isotropic", hole)
        other = self.vm.Reader(o["coords"], [p["symbol"] for p in o["parameters"]], ())
        values = {tuple(e["indices"]): e["value"] for e in o["metric_components"]}
        for i, a in enumerate(o["coords"]):
            mine = g[i, i].subs({R.parameters["b"]: 0, R.parameters["mu"]: 1, R.symbol["\\rho"]: 2,
                                 R.symbol["\\theta"]: 1})
            there = other(values.get((a, a), "0")).subs({other.parameters["m"]: 1, other.symbol["r"]: 2,
                                                         other.symbol["\\theta"]: 1})
            with self.subTest(slot=a):
                self.assertAlmostEqual(float(mine), float(there), places=12)

    def test_the_chart_in_r_is_omega_squared_times_fisher_janis_newman_winicour(self):
        sp = self.sp
        _, R, g, _, _ = self.published("three_parameter")
        seed = json.loads((DATA / "metrics" / "fisher_jnw.json").read_text(encoding="utf-8"))
        o = chart("spherical", seed)
        other = self.vm.Reader(o["coords"], [p["symbol"] for p in o["parameters"]], ())
        values = {tuple(e["indices"]): e["value"] for e in o["metric_components"]}
        P, r = R.parameters, R.symbol["r"]
        omega = R.held[P["Omega"]]
        for mu, b, w, rv in ((1, -sp.Rational(1, 2), sp.Rational(1, 6), 5), (sp.Rational(3, 2), -2, 1, 9)):
            gamma = (1 - 12 * w ** 2) / (1 + 12 * w ** 2)
            at = {P["mu"]: mu, P["b"]: b, P["w"]: w, r: rv, R.symbol["\\theta"]: 1}
            for i, a in enumerate(o["coords"]):
                there = other(values[(a, a)]).subs({other.parameters["b"]: 4 * mu, other.parameters["gamma"]: gamma,
                                                    other.symbol["r"]: rv, other.symbol["\\theta"]: 1})
                with self.subTest(at=(mu, b, w, rv), slot=a):
                    self.assertAlmostEqual(float(sp.N(g[i, i].subs(at), 30)),
                                           float(sp.N((omega ** 2).subs(at) * there, 30)), places=12)

    def test_at_w_one_sixth_the_chart_in_r_is_the_isotropic_chart(self):
        sp = self.sp
        _, R3, g3, _, _ = self.published("three_parameter")
        _, Ri, gi, _, _ = self.published("isotropic")
        for mu, b, rho in ((1, -sp.Rational(1, 4), sp.Rational(3, 2)), (sp.Rational(1, 2), -3, 2)):
            r = rho * (1 + mu / rho) ** 2
            dr = 1 - (mu / rho) ** 2
            here = {R3.parameters["mu"]: mu, R3.parameters["b"]: b, R3.parameters["w"]: sp.Rational(1, 6),
                    R3.symbol["r"]: r, R3.symbol["\\theta"]: 1}
            there = {Ri.parameters["mu"]: mu, Ri.parameters["b"]: b, Ri.symbol["\\rho"]: rho, Ri.symbol["\\theta"]: 1}
            for i, scale in enumerate((1, dr ** 2, 1, 1)):
                with self.subTest(at=(mu, b, rho), slot=i):
                    self.assertAlmostEqual(float(sp.N(g3[i, i].subs(here) * scale, 30)),
                                           float(sp.N(gi[i, i].subs(there), 30)), places=12)


class Texts(unittest.TestCase):
    """The History states the source's formulas, and every citation is held."""

    def test_history_states_the_solution(self):
        history = METRIC["history"]
        for formula in ("$\\tan(\\Delta/2) = (\\alpha_+ - \\alpha_-)/(\\alpha_+ + \\alpha_-)$",
                        "$\\bar r = (\\eta/2)|\\tan(\\Delta/2)|^{1/2}$", "$1/(8\\pi(\\kappa - \\phi_c^2/6))$",
                        "$\\phi_c^2 = 6\\kappa$", "$\\kappa = 1/(8\\pi G)$", "$r = 2\\eta$"):
            with self.subTest(formula=formula):
                self.assertIn(formula, history)
        self.assertEqual([len(p) > 0 for p in history.split("¶")], [True] * 6)

    def test_the_conventions_name_the_source_s_parameters(self):
        self.assertIn("$\\eta$", METRIC["convention"])
        self.assertIn("$b = \\tan(\\Delta/2)$", METRIC["convention"])
        self.assertIn("$\\bar r$", chart("isotropic")["convention"])
        self.assertIn("$w = \\tan(\\chi/2)/(2\\sqrt{3})$", chart("three_parameter")["convention"])

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
        kinds = {"fisher_jnw": "conformal", "majumdar_papapetrou": "conformal"}
        self.assertEqual(len(METRIC["related"]), 8)
        for entry in METRIC["related"]:
            other = json.loads((DATA / "metrics" / f"{entry['id']}.json").read_text(encoding="utf-8"))
            answer = [e for e in other["related"] if e["id"] == ME]
            with self.subTest(other=entry["id"]):
                self.assertEqual(entry["kind"], kinds.get(entry["id"], "family"))
                self.assertEqual(len(answer), 1)
                self.assertEqual(answer[0]["kind"], entry["kind"])
                for key in cited(answer[0]["text"]):
                    self.assertIn(key, other["references"])


class Drawings(unittest.TestCase):
    """The drawings at mu = 1 and b = -1/4, whose throat is rho = mu/2."""

    @staticmethod
    def areal(rho):
        return (rho + 0.25) * (rho + 1) / rho

    def test_embedding_is_one_piece_through_the_throat_the_same_on_both_sides(self):
        (view,) = drawing("embedding")["views"]
        self.assertEqual(view["id"], "wormhole")
        (piece,) = view["surfaces"][0]["pieces"]
        points = piece["points"]
        for x, rho, _ in points:
            self.assertAlmostEqual(rho, self.areal(x), places=5)
        self.assertAlmostEqual(points[0][0], 0.05, places=12)
        self.assertAlmostEqual(points[-1][0], 5.0, places=12)
        self.assertAlmostEqual(points[0][1], points[-1][1], places=5)
        throat = min(points, key=lambda p: p[1])
        self.assertAlmostEqual(throat[0], 0.5, places=2)
        self.assertAlmostEqual(throat[1], 2.25, places=4)
        self.assertAlmostEqual(throat[2] - points[0][2], points[-1][2] - throat[2], places=2)

    def test_conformal_view_is_the_full_diamond_with_the_throat_on_its_axis(self):
        (view,) = drawing("conformal")["views"]
        self.assertEqual(view["id"], "isotropic")
        (throat,) = [layer for layer in view["layers"] if layer["class"] == "throat"]
        self.assertEqual(throat["points"], [[0.0, -round(math.pi, 4)], [0.0, round(math.pi, 4)]])
        self.assertEqual(len([layer for layer in view["layers"] if layer["class"] == "infinity"]), 4)
        self.assertEqual([s["view"] for s in view["slices"]], ["wormhole"])

    def test_curvature_figure_plots_k_against_the_proper_distance_and_at_the_naked_member(self):
        views = drawing("curvature")["views"]
        self.assertEqual([v["id"] for v in views], ["wormhole", "naked"])
        for view in views:
            self.assertEqual([p["id"] for p in view["panels"]], ["K"])
            self.assertIn("The Ricci scalar is zero.", view["caption"][0])
        (curve,) = views[0]["curves"]
        rhos = curve["horizontal"]["values"]
        for (l, _), rho in zip(curve["points"], rhos):
            want = rho - 0.5 + 1.25 * math.log(2 * rho) - 0.25 * (1 / rho - 2)
            self.assertAlmostEqual(l, want, places=9)
        self.assertEqual(views[0]["marks"], [{"class": "front", "x": 0.0}])
        self.assertEqual({c["params"]["w"] for c in views[1]["curves"]}, {0.25})

    def test_spacetime_diagrams_mark_the_moment_embedded_on_the_wormhole_alone(self):
        systems = drawing("diagrams")["systems"]
        self.assertEqual(sorted(systems), sorted(SYSTEMS))
        (iso,) = systems["isotropic"]
        (naked,) = systems["three_parameter"]
        self.assertEqual(iso["box"], [0, 3, -1.5, 1.5])
        self.assertEqual(naked["box"], [4, 7, -1.5, 1.5])
        self.assertEqual([s["view"] for s in iso["slices"]], ["wormhole"])
        self.assertNotIn("slices", naked)
        self.assertEqual([m["kind"] for m in iso["markers"]], ["throat", "surface"])
        self.assertEqual([m["kind"] for m in naked["markers"]], ["grr", "singular"])


if __name__ == "__main__":
    unittest.main()
