"""Kasner's universe with a magnetic field, ds^2 = f^2(-c^2dt^2 + t^(2 p_1)dx^2 + t^(2 p_2)dy^2)
+ t^(2 p_3)dz^2/f^2 with f = 1 + b^2 t^(2 p_3): .venv.noindex/bin/python -m unittest discover -s _tools

What its texts and drawings state, held on the published files. The first class needs nothing
but the files: the exponents the drawings declare and the ones they end at, their square
windows, the ring of the embedding diagram, and the relations each side answers. The second
reads the published tensors through the checker's Reader and holds them to the stress of a
magnetic field along z, to Kasner's vacuum at b = 0, to a finite curvature at t = 0 in the
axisymmetric universe, and to Rosen's chart being that universe; it needs sympy and is skipped
where it is absent.
"""
import importlib.util
import json
import sys
import unittest
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
HAS_SYMPY = importlib.util.find_spec("sympy") is not None
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))

CHARTS = ["kasner_time", "rosen"]
P = (Fraction(-2, 7), Fraction(3, 7), Fraction(6, 7))


def load(folder, metric_id="kasner_magnetic"):
    return json.loads((DATA / folder / f"{metric_id}.json").read_text(encoding="utf-8"))


def late(p):
    """The exponents at late times, (24) of Kastor and Traschen: the field along the third axis."""
    k = 1 + 2 * p[2]
    return ((p[0] + 2 * p[2]) / k, (p[1] + 2 * p[2]) / k, -p[2] / k)


class Drawings(unittest.TestCase):
    def test_the_drawn_exponents_and_the_ones_they_end_at_lie_on_kasners_circle(self):
        for p in (P, late(P)):
            self.assertEqual(sum(p), 1)
            self.assertEqual(sum(k * k for k in p), 1)
        self.assertGreater(P[2], 0)
        self.assertEqual(late(P), (Fraction(10, 19), Fraction(15, 19), Fraction(-6, 19)))
        # The axisymmetric universe runs from flat space to the vacuum with one axis shrinking.
        self.assertEqual(late((Fraction(0), Fraction(0), Fraction(1))), (Fraction(2, 3), Fraction(2, 3), Fraction(-1, 3)))

    def test_every_chart_has_its_spacetime_diagrams_in_square_windows(self):
        systems = load("diagrams")["systems"]
        self.assertEqual(sorted(systems), sorted(CHARTS))
        self.assertEqual({k: [v["id"] for v in views] for k, views in systems.items()},
                         {"kasner_time": ["tx", "tz"], "rosen": ["etax", "etaz"]})
        for views in systems.values():
            for view in views:
                X0, X1, Y0, Y1 = view["box"]
                self.assertAlmostEqual((X1 - X0) / (Y1 - Y0), 1.0, places=9, msg=view["id"])

    def test_the_ring_is_longest_along_the_field_at_t_equal_to_1(self):
        """The ring at t is the ellipse of semi-axes f t^(-2/7) across the field and t^(6/7)/f along
        it, f = 1 + t^(12/7), in units of its radius in the chart."""
        views = {v["id"]: v for v in load("embedding")["views"]}
        self.assertEqual(sorted(views), ["ring", "tube"])
        surfaces = views["ring"]["surfaces"]
        self.assertEqual([s["time"] for s in surfaces], [0.25, 0.5, 1.0, 2.0])
        along = []
        for surface in surfaces:
            points = surface["curves"][0]["points"]
            a, b = max(abs(q[0]) for q in points), max(abs(q[1]) for q in points)
            t = surface["time"]
            f = 1 + t ** (12 / 7)
            self.assertAlmostEqual(a, f * t ** (-2 / 7), delta=4e-3)
            self.assertAlmostEqual(b, t ** (6 / 7) / f, delta=2e-3)
            along.append(b)
        self.assertEqual(max(along), along[2])
        self.assertAlmostEqual(along[2], 0.5, delta=2e-3)

    def test_it_has_no_conformal_diagram(self):
        self.assertFalse((DATA / "conformal" / "kasner_magnetic.json").exists())

    def test_each_relation_is_answered(self):
        answers = {"kasner": ("special_case", "generalisation"), "melvin": ("dual", "dual"),
                   "kasner_scalar": ("family", "family"), "bianchi": ("generalisation", "special_case")}
        mine = {r["id"]: r["kind"] for r in load("metrics")["related"]}
        self.assertEqual(mine, {k: v[0] for k, v in answers.items()})
        for other, (_, kind) in answers.items():
            back = [r for r in load("metrics", other)["related"] if r["id"] == "kasner_magnetic"]
            self.assertEqual([r["kind"] for r in back], [kind], other)


@unittest.skipUnless(HAS_SYMPY, "needs sympy: run under .venv.noindex/bin/python")
class Physics(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import sympy as sp
        import verify_metrics as vm
        cls.sp, cls.vm = sp, vm
        cls.charts = {c["id"]: c for c in load("metrics")["coordinates"]}

    def read(self, system, field, variant=None):
        """The chart's reader and one published tensor of rank two as a matrix, or a scalar."""
        sp = self.sp
        chart = self.charts[system]
        reader = self.vm.Reader(chart["coords"], [p["symbol"] for p in chart["parameters"]], ())
        if variant is None and isinstance(chart[field], str):
            return reader, reader(chart[field].partition("=")[2])
        entries = chart[field] if variant is None else chart[field]["variants"][variant]["nonzero"]
        out = sp.zeros(4)
        for e in entries:
            i, j = (chart["coords"].index(k) for k in e["indices"])
            out[i, j] = reader(e["value"])
        return reader, out

    def test_the_stress_is_that_of_a_magnetic_field_along_z(self):
        """Tension along the field and pressure across it: R^t_t = R^z_z = -R^x_x = -R^y_y, and no trace."""
        sp = self.sp
        for system in CHARTS:
            _, mixed = self.read(system, "ricci_tensor", "ul")
            self.assertNotEqual(mixed[0, 0], 0, system)
            self.assertEqual(sp.simplify(mixed[3, 3] - mixed[0, 0]), 0, system)
            for k in (1, 2):
                self.assertEqual(sp.simplify(mixed[k, k] + mixed[0, 0]), 0, system)
            self.assertEqual(self.charts[system]["ricci_scalar"], "R = 0")

    def test_the_field_the_parameter_states_is_the_one_in_the_curvature(self):
        """G B^2/c^4 = -R^t_t with B = (2 b p_3 c/sqrt(G)) t^(p_3 - 1)/f^2, as the parameter b is described."""
        sp = self.sp
        reader, mixed = self.read("kasner_time", "ricci_tensor", "ul")
        t, c, b, p3 = reader.symbol["t"], reader.c, reader.parameters["b"], reader.parameters["p_3"]
        time = sp.Symbol("T", positive=True)
        f = 1 + b ** 2 * time ** (2 * p3)
        B2 = (2 * b * p3) ** 2 * time ** (2 * p3 - 2) / f ** 4          # G B^2/c^2, with T = t
        have = mixed[0, 0].subs(t, time)
        self.assertEqual(sp.simplify(sp.powsimp(sp.powdenest(have + B2 / c ** 2, force=True), force=True)), 0)

    def test_with_no_field_the_metric_and_curvature_are_kasners(self):
        sp = self.sp
        reader, g = self.read("kasner_time", "metric_components")
        b = reader.parameters["b"]
        kasner = next(c for c in load("metrics", "kasner")["coordinates"] if c["id"] == "cartesian")
        theirs = self.vm.Reader(kasner["coords"], [p["symbol"] for p in kasner["parameters"]], ())
        same = {theirs.parameters[k]: reader.parameters[k] for k in ("p_1", "p_2", "p_3")}
        same[theirs.symbol["t"]] = reader.symbol["t"]
        for e in kasner["metric_components"]:
            i, j = (kasner["coords"].index(k) for k in e["indices"])
            self.assertEqual(sp.simplify(theirs(e["value"]).subs(same) - g[i, j].subs(b, 0)), 0, e["indices"])
        _, K = self.read("kasner_time", "kretschmann")
        p3, t, c = reader.parameters["p_3"], reader.symbol["t"], reader.c
        self.assertEqual(sp.simplify(K.subs(b, 0) - 16 * p3 ** 2 * (1 - p3) / (c * t) ** 4), 0)
        # On Kasner's circle 16 p_3^2 (1 - p_3) is -16 p_1 p_2 p_3, as Kasner's page has it.
        for p in (P, late(P)):
            self.assertEqual(16 * p[2] ** 2 * (1 - p[2]), -16 * p[0] * p[1] * p[2])

    def test_the_axisymmetric_universe_has_finite_curvature_at_its_first_instant(self):
        """At (0, 0, 1) the Kretschmann scalar is 64 b^4 (5 - 6 b^2 t^2 + 3 b^4 t^4)/f^8, 320 b^4 at t = 0,
        with t in the unit the powers are read in; at any other p_3 below 1 it diverges as t^-4."""
        sp = self.sp
        reader, K = self.read("kasner_time", "kretschmann")
        t, c, b, p3 = reader.symbol["t"], reader.c, reader.parameters["b"], reader.parameters["p_3"]
        time = sp.Symbol("T", positive=True)
        at_one = sp.simplify(K.subs(p3, 1).subs(t, time))
        self.assertEqual(sp.simplify(at_one - 64 * b ** 4 * (5 - 6 * b ** 2 * time ** 2 + 3 * b ** 4 * time ** 4)
                                     / (c ** 4 * (1 + b ** 2 * time ** 2) ** 8)), 0)
        general = K.subs(p3, sp.Rational(6, 7)).subs(t, time) * time ** 4
        self.assertEqual(sp.limit(general, time, 0), 16 * sp.Rational(36, 49) * sp.Rational(1, 7) / c ** 4)

    def test_rosens_chart_is_the_axisymmetric_universe(self):
        """tan(eta/2) = b t at (0, 0, 1): the Kretschmann scalars agree with ell = 2 c t_0/b, and the
        curvature at eta = 0 is 5120/ell^4, as the diagram's caption states."""
        sp = self.sp
        r1, K1 = self.read("kasner_time", "kretschmann")
        r2, K2 = self.read("rosen", "kretschmann")
        eta, ell = r2.symbol["\\eta"], r2.parameters["ell"]
        self.assertEqual(sp.simplify(K2.subs(eta, 0) - 5120 / ell ** 4), 0)
        b, time = sp.symbols("b T", positive=True)
        first = K1.subs({r1.parameters["p_3"]: 1, r1.parameters["b"]: b}).subs(r1.symbol["t"], time)
        # With t in the unit t_0, c t_0 = b ell/2.
        first = first.subs(r1.c, b * ell / 2)
        second = K2.subs(eta, 2 * sp.atan(b * time))
        self.assertEqual(sp.simplify(first - second), 0)


if __name__ == "__main__":
    unittest.main()
