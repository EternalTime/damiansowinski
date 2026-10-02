"""Kasner's universe with a scalar field, ds^2 = -c^2dt^2 + sum t^(2 p_i) dx_i^2 with sum p_i = 1
and sum p_i^2 = 1 - q^2: python3 -m unittest discover -s _tools

What its texts and drawings state, held on the published files. The first class needs nothing
but the files: the exponents the drawings declare, their square windows, the ring of the
embedding diagram, and the relations each side answers. The second reads the published tensors
through the checker's Reader and holds them, on the surface the exponents live on, to the stress
of a stiff fluid, to the wave equation, to Kasner's vacuum at q = 0, to Friedmann's universe at
q^2 = 2/3, and to the reduction of the vacuum of five dimensions; it needs sympy and is skipped
where it is absent.
"""
import importlib.util
import json
import math
import sys
import unittest
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
HAS_SYMPY = importlib.util.find_spec("sympy") is not None
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))

CHARTS = ["synchronous", "logarithmic", "kaluza_klein"]
P = (Fraction(2, 13), Fraction(4, 13), Fraction(7, 13))
Q = Fraction(10, 13)


def load(folder, metric_id="kasner_scalar"):
    return json.loads((DATA / folder / f"{metric_id}.json").read_text(encoding="utf-8"))


def five(p=P, q=Q):
    """The exponents of the vacuum of five dimensions that reduces to (p, q): Belinskii and
    Khalatnikov's fourth footnote."""
    root = math.sqrt(6)
    return [(root * float(k) - float(q)) / (root - float(q)) for k in p] + [2 * float(q) / (root - float(q))]


class Drawings(unittest.TestCase):
    def test_the_drawn_exponents_lie_on_the_surface_and_are_all_positive(self):
        self.assertEqual(sum(P), 1)
        self.assertEqual(sum(k * k for k in P), 1 - Q * Q)
        self.assertTrue(all(k > 0 for k in P))
        # Past q^2 = 1/2 no exponent can be negative, and q^2 stays below 2/3.
        self.assertGreater(Q * Q, Fraction(1, 2))
        self.assertLess(Q * Q, Fraction(2, 3))

    def test_the_exponents_of_five_dimensions_are_a_vacuum_and_reduce_to_the_drawn_ones(self):
        s = five()
        self.assertAlmostEqual(sum(s), 1, places=12)
        self.assertAlmostEqual(sum(k * k for k in s), 1, places=12)
        self.assertLess(s[0], 0)
        for k, p in zip(s, P):
            self.assertAlmostEqual((2 * k + s[3]) / (2 + s[3]), float(p), places=12)
        self.assertAlmostEqual(math.sqrt(6) * s[3] / (2 + s[3]), float(Q), places=12)
        for k, shown in zip(s, (-0.234, -0.009, 0.327, 0.916)):
            self.assertAlmostEqual(k, shown, delta=5e-4)

    def test_every_chart_has_its_spacetime_diagrams_in_square_windows(self):
        systems = load("diagrams")["systems"]
        self.assertEqual(sorted(systems), sorted(CHARTS))
        self.assertEqual({k: [v["id"] for v in views] for k, views in systems.items()},
                         {"synchronous": ["tx", "tz"], "logarithmic": ["taux"], "kaluza_klein": ["Tx", "Tw"]})
        for views in systems.values():
            for view in views:
                X0, X1, Y0, Y1 = view["box"]
                self.assertAlmostEqual((X1 - X0) / (Y1 - Y0), 1.0, places=9, msg=view["id"])

    def test_the_ring_shrinks_along_both_axes_toward_the_singularity(self):
        """The ring at t is the ellipse of semi-axes t^(2/13) and t^(7/13), in units of its radius at t = 1."""
        views = {v["id"]: v for v in load("embedding")["views"]}
        self.assertEqual(sorted(views), ["ring", "tube"])
        surfaces = views["ring"]["surfaces"]
        self.assertEqual([s["time"] for s in surfaces], [0.25, 0.5, 1.0, 2.0])
        reach = []
        for surface in surfaces:
            points = surface["curves"][0]["points"]
            a, b = max(abs(q[0]) for q in points), max(abs(q[1]) for q in points)
            t = surface["time"]
            self.assertAlmostEqual(a, t ** (2 / 13), delta=2e-3)
            self.assertAlmostEqual(b, t ** (7 / 13), delta=2e-3)
            reach.append((a, b))
        for (a0, b0), (a1, b1) in zip(reach, reach[1:]):
            self.assertLess(a0, a1)
            self.assertLess(b0, b1)

    def test_it_has_no_conformal_diagram(self):
        self.assertFalse((DATA / "conformal" / "kasner_scalar.json").exists())

    def test_each_relation_is_answered(self):
        answers = {"kasner": ("special_case", "generalisation"), "frw": ("special_case", "generalisation"),
                   "bianchi": ("generalisation", "special_case"), "mixmaster": ("family", "family"),
                   "fisher_jnw": ("family", "family")}
        mine = {r["id"]: r["kind"] for r in load("metrics")["related"]}
        self.assertEqual(mine, {k: v[0] for k, v in answers.items()})
        for other, (_, kind) in answers.items():
            back = [r for r in load("metrics", other)["related"] if r["id"] == "kasner_scalar"]
            self.assertEqual([r["kind"] for r in back], [kind], other)


@unittest.skipUnless(HAS_SYMPY, "needs sympy: run under /tmp/mfs-venv/bin/python")
class Physics(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import sympy as sp
        import verify_metrics as vm
        cls.sp, cls.vm = sp, vm
        cls.charts = {c["id"]: c for c in load("metrics")["coordinates"]}

    def read(self, system):
        """The chart's reader and its published metric, Ricci and Einstein tensors as matrices."""
        sp, vm = self.sp, self.vm
        chart = self.charts[system]
        reader = vm.Reader(chart["coords"], [p["symbol"] for p in chart["parameters"]], ())
        n = len(chart["coords"])

        def matrix(entries):
            out = sp.zeros(n)
            for e in entries:
                i, j = (chart["coords"].index(k) for k in e["indices"])
                out[i, j] = reader(e["value"])
            return out
        return reader, matrix(chart["metric_components"]), {
            name: matrix(chart[field]["variants"][variant]["nonzero"])
            for name, field, variant in (("ricci", "ricci_tensor", "ll"), ("mixed", "einstein_tensor", "ul"))}

    def scalar(self, system, field):
        chart = self.charts[system]
        reader = self.vm.Reader(chart["coords"], [p["symbol"] for p in chart["parameters"]], ())
        return reader, reader(chart[field].partition("=")[2])

    def test_the_stress_is_that_of_a_fluid_at_rest_with_pressure_equal_to_energy_density(self):
        """G^t_t = -8 pi G rho/c^4 and G^x_x = G^y_y = G^z_z = 8 pi G p/c^4 with p = rho."""
        for system in ("synchronous", "logarithmic"):
            _, _, tensors = self.read(system)
            G = tensors["mixed"]
            self.assertNotEqual(G[0, 0], 0, system)
            for k in (1, 2, 3):
                self.assertEqual(self.sp.simplify(G[k, k] + G[0, 0]), 0, system)

    def test_the_ricci_tensor_is_the_square_of_the_fields_gradient(self):
        """R_tt = (d_t phi)^2 with phi = q ln t in the chart x^0 = ct, and nothing else."""
        sp = self.sp
        reader, _, tensors = self.read("synchronous")
        t, q, c = reader.symbol["t"], reader.parameters["q"], reader.c
        ricci = tensors["ricci"]
        self.assertEqual(sp.simplify(ricci[0, 0] - q ** 2 / (c * t) ** 2), 0)
        self.assertEqual(sum(1 for e in ricci if e != 0), 1)
        reader, _, tensors = self.read("logarithmic")
        self.assertEqual(sp.simplify(tensors["ricci"][0, 0] - reader.parameters["q"] ** 2), 0)

    def test_the_energy_density_is_the_one_the_history_states(self):
        """rho c^2 = q^2 c^2/(16 pi G t^2): G_tt = q^2/(2 c^2 t^2) in the chart x^0 = ct."""
        sp = self.sp
        chart = self.charts["synchronous"]
        reader = self.vm.Reader(chart["coords"], [p["symbol"] for p in chart["parameters"]], ())
        (entry,) = [e for e in chart["einstein_tensor"]["variants"]["ll"]["nonzero"] if e["indices"] == ["t", "t"]]
        t, q, c = reader.symbol["t"], reader.parameters["q"], reader.c
        self.assertEqual(sp.simplify(reader(entry["value"]) - q ** 2 / (2 * c ** 2 * t ** 2)), 0)

    def test_with_no_field_the_curvature_is_that_of_kasners_vacuum(self):
        """At q = 0 the Kretschmann scalar is -16 p_1 p_2 p_3/(c t)^4, as Kasner's page has it."""
        sp = self.sp
        reader, K = self.scalar("synchronous", "kretschmann")
        P_ = reader.parameters
        t, c = reader.symbol["t"], reader.c
        vacuum = -16 * P_["p_1"] * P_["p_2"] * P_["p_3"] / (c * t) ** 4
        self.assertEqual(sp.simplify(K.subs(P_["q"], 0) - vacuum), 0)
        kasner = next(c_ for c_ in load("metrics", "kasner")["coordinates"] if c_["id"] == "cartesian")
        theirs = self.vm.Reader(kasner["coords"], [p["symbol"] for p in kasner["parameters"]], ())
        published = theirs(kasner["kretschmann"].partition("=")[2])
        same = {theirs.parameters[k]: P_[k] for k in ("p_1", "p_2", "p_3")}
        same[theirs.symbol["t"]] = t
        self.assertEqual(sp.simplify(published.subs(same) - vacuum), 0)

    def test_at_the_isotropic_point_the_weyl_tensor_vanishes(self):
        """p_1 = p_2 = p_3 = 1/3 and q^2 = 2/3: the flat Friedmann universe of a stiff fluid."""
        sp = self.sp
        for system in ("synchronous", "logarithmic"):
            chart = self.charts[system]
            reader = self.vm.Reader(chart["coords"], [p["symbol"] for p in chart["parameters"]], ())
            at = {reader.parameters[k]: sp.Rational(1, 3) for k in ("p_1", "p_2", "p_3")}
            at[reader.parameters["q"]] = sp.sqrt(sp.Rational(2, 3))
            entries = chart["weyl_tensor"]["variants"]["llll"]["nonzero"]
            self.assertTrue(entries)
            for e in entries:
                self.assertEqual(sp.simplify(reader(e["value"]).subs(at)), 0, (system, e["indices"]))

    def test_the_kretschmann_scalar_agrees_in_the_two_times(self):
        """t = t_0 e^(-tau): (c t)^4 is ell^4 e^(-4 tau)."""
        sp = self.sp
        r1, K1 = self.scalar("synchronous", "kretschmann")
        r2, K2 = self.scalar("logarithmic", "kretschmann")
        same = {r1.parameters[k]: r2.parameters[k] for k in ("p_1", "p_2", "p_3", "q")}
        same[r1.symbol["t"]] = r2.parameters["ell"] * sp.exp(-r2.symbol["\\tau"]) / r1.c
        self.assertEqual(sp.simplify(K1.subs(same) - K2), 0)

    def test_the_chart_of_five_dimensions_is_a_vacuum(self):
        chart = self.charts["kaluza_klein"]
        self.assertEqual(len(chart["coords"]), 5)
        for field in ("ricci_tensor", "einstein_tensor"):
            for variant in chart[field]["variants"].values():
                self.assertEqual(variant["nonzero"], [])
        self.assertEqual(chart["ricci_scalar"], "R = 0")
        # On a point of the surface the Kretschmann scalar is 4(sum s^2 (s - 1)^2 + sum s_a^2 s_b^2)/(c T)^4.
        sp = self.sp
        reader, K = self.scalar("kaluza_klein", "kretschmann")
        relations = self.vm.PARAMETER_RELATIONS[("kasner_scalar", "kaluza_klein")]
        at = {sp.Symbol("a"): sp.Rational(1, 3), sp.Symbol("b"): sp.Rational(-2, 5)}
        s = [sp.sympify(relations[k]).subs(at) for k in ("s_1", "s_2", "s_3", "s_5")]
        self.assertEqual(sum(s), 1)
        self.assertEqual(sum(k ** 2 for k in s), 1)
        want = 4 * (sum(k ** 2 * (k - 1) ** 2 for k in s)
                    + sum(s[i] ** 2 * s[j] ** 2 for i in range(4) for j in range(i + 1, 4)))
        values = dict(zip((reader.parameters[k] for k in ("s_1", "s_2", "s_3", "s_5")), s))
        T, c = reader.symbol["T"], reader.c
        self.assertEqual(sp.simplify(K.subs(values) * (c * T) ** 4 - want), 0)


if __name__ == "__main__":
    unittest.main()
