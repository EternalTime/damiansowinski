"""The A- and B-metrics of Ehlers and Kundt, the static vacuum fields of type D beside Schwarzschild's:
python3 -m unittest discover -s _tools

What the texts of the spacetime state, held on the published files: that BI is Schwarzschild's
metric with the time and the angle exchanged, that AII and AIII are the hyperbolic and the flat
topological black hole without the cosmological constant, that AIII with b made negative is Kasner's
universe, that BIII and AIII are Levi-Civita's cylinder at sigma = 1/4 and -1/2, that AII and BI
are flat at b = 0, that AIII's singularity repels, that Kruskal's radius solves
(r/b - 1) e^(r/b) = UV, and what the captions say of the drawings. The tests of the metric need
sympy and are skipped where it is absent; the tests of the drawings read the numbers in the files.
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


def read(kind, metric_id="ab_metrics"):
    return json.loads((DATA / kind / f"{metric_id}.json").read_text(encoding="utf-8"))


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class Family(unittest.TestCase):
    def setUp(self):
        import sympy
        import verify_metrics
        self.sp, self.vm = sympy, verify_metrics

    def chart(self, metric_id, system):
        """A published chart: its reader, its coordinates' symbols and its diagonal metric by coordinate."""
        entry = next(c for c in read("metrics", metric_id)["coordinates"] if c["id"] == system)
        held = self.vm.HELD.get((metric_id, system), ())
        reader = self.vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], (), held=held)
        g = {tuple(e["indices"]): reader(e["value"]) for e in entry["metric_components"]}
        return entry, reader, g

    def same(self, a, b):
        self.assertEqual(self.sp.simplify(a - b), 0, f"{a} is not {b}")

    def test_b_one_is_schwarzschild_with_the_time_and_the_angle_exchanged(self):
        _, there, schwarzschild = self.chart("schwarzschild", "spherical")
        _, here, b1 = self.chart("ab_metrics", "b1_static")
        swap = {there.symbol["r"]: here.symbol["r"], there.symbol["\\theta"]: here.symbol["\\theta"],
                there.parameters["r_s"]: here.parameters["b"]}
        at = lambda key: schwarzschild[key].subs(swap)
        # phi -> i tau and ct -> i z: g_tautau = -g_phiphi and g_zz = -g_tt.
        self.same(b1[("\\tau", "\\tau")], -at(("\\phi", "\\phi")))
        self.same(b1[("z", "z")], -at(("t", "t")))
        self.same(b1[("r", "r")], at(("r", "r")))
        self.same(b1[("\\theta", "\\theta")], at(("\\theta", "\\theta")))

    def test_a_two_and_a_three_are_topological_black_holes_without_the_cosmological_constant(self):
        sp = self.sp
        _, there, hole = self.chart("topological_black_hole", "hyperbolic")
        _, here, a2 = self.chart("ab_metrics", "a2_static")
        b, r = here.parameters["b"], here.symbol["r"]
        swap = {there.symbol["r"]: r, there.parameters["mu"]: -b, there.symbol["\\theta"]: here.symbol["\\chi"]}
        limit = lambda e: sp.limit(e.subs(swap), there.parameters["L"], sp.oo)
        for mine, theirs in ((("t", "t"), ("t", "t")), (("r", "r"), ("r", "r")),
                             (("\\phi", "\\phi"), ("\\phi", "\\phi"))):
            self.same(a2[mine].rewrite(sp.exp), limit(hole[theirs]).rewrite(sp.exp))
        _, there, hole = self.chart("topological_black_hole", "static")
        _, here, a3 = self.chart("ab_metrics", "a3")
        b, r = here.parameters["b"], here.symbol["r"]
        swap = {there.symbol["r"]: r, there.parameters["mu"]: -b, there.parameters["k"]: 0,
                there.symbol["\\rho"]: here.symbol["\\chi"]}
        limit = lambda e: sp.limit(e.subs(swap), there.parameters["L"], sp.oo)
        self.same(a3[("t", "t")], limit(hole[("t", "t")]))
        self.same(a3[("r", "r")], limit(hole[("r", "r")]))
        self.same(a3[("\\chi", "\\chi")], limit(hole[("\\rho", "\\rho")]))
        self.same(a3[("\\phi", "\\phi")], limit(hole[("\\phi", "\\phi")]))

    def test_a_three_with_b_negative_is_kasners_universe(self):
        sp = self.sp
        _, here, a3 = self.chart("ab_metrics", "a3")
        b, r = here.parameters["b"], here.symbol["r"]
        beta, T = sp.symbols("beta T", positive=True)
        flipped = {key: value.subs(b, -beta) for key, value in a3.items()}
        # The radius is now the time: (r/beta) dr^2 = dT^2 with r = (3T/2)^(2/3) beta^(1/3).
        radius = (3 * T / 2) ** sp.Rational(2, 3) * beta ** sp.Rational(1, 3)
        self.same((flipped[("r", "r")] * sp.diff(radius, T) ** 2).subs(r, radius), -1)
        for key, exponent in ((("t", "t"), sp.Rational(-1, 3)), (("\\chi", "\\chi"), sp.Rational(2, 3))):
            along = flipped[key].subs(r, radius)
            self.assertGreater(along.subs({T: 1, beta: 1}), 0)
            self.same(T * sp.diff(along, T) / along, 2 * exponent)

    def test_b_three_and_a_three_are_levi_civitas_cylinder(self):
        sp = self.sp
        _, there, cylinder = self.chart("levi_civita", "weyl")
        sigma, rho = there.parameters["sigma"], there.symbol["\\rho"]
        k = sp.Symbol("k", positive=True)
        for system, value, power, pairs in (
                ("b3", sp.Rational(1, 4), 2, ((("\\tau", "\\tau"), ("t", "t")), (("z", "z"), ("z", "z")))),
                ("a3", sp.Rational(-1, 2), sp.Rational(1, 2), ((("t", "t"), ("t", "t")), (("\\chi", "\\chi"), ("z", "z"))))):
            _, here, mine = self.chart("ab_metrics", system)
            r = here.symbol["r"]
            image = k * r ** power
            ratios = [mine[key] / cylinder[other].subs(sigma, value).subs(rho, image) for key, other in pairs]
            ratios.append(mine[("r", "r")] / (cylinder[("\\rho", "\\rho")].subs(sigma, value).subs(rho, image)
                                              * sp.diff(image, r) ** 2))
            for ratio in ratios:
                # The two agree up to a constant, which rescales the coordinate.
                numbers = [complex(sp.N(ratio.subs({r: x, k: sp.Rational(7, 5), here.parameters["b"]: sp.Rational(3, 2)}), 30))
                           for x in (sp.Rational(1, 2), 2, 5)]
                self.assertLess(max(abs(n - numbers[0]) for n in numbers), 1e-20)
                self.assertGreater(numbers[0].real, 0)

    def test_a_two_and_b_one_are_flat_without_b(self):
        sp = self.sp
        for system in ("a2_cartesian", "b1_cartesian"):
            entry, reader, g = self.chart("ab_metrics", system)
            b = reader.parameters["b"]
            for i, a in enumerate(entry["coords"]):
                for j, c in enumerate(entry["coords"]):
                    value = g.get((a, c), sp.Integer(0))
                    flat = (-1 if i == 0 else 1) if i == j else 0
                    self.same(sp.sympify(value).subs(b, 0), flat)
        for system in ("a3", "b2_static", "b3"):
            entry, reader, g = self.chart("ab_metrics", system)
            # No flat member: at b = 0 a component vanishes or diverges, or the signature is another.
            at = {symbol: sp.Rational(7, 10) for symbol in reader.symbol.values()}
            values = [sp.sympify(v).subs(reader.parameters["b"], 0).subs(at) for v in g.values()]
            broken = any(v == 0 or not v.is_finite for v in values)
            self.assertTrue(broken or sum(1 for v in values if v < 0) != 1, system)

    def test_every_chart_states_the_curvature_of_its_own_radius(self):
        sp = self.sp
        for chart in read("metrics")["coordinates"]:
            _, reader, _ = self.chart("ab_metrics", chart["id"])
            b = reader.parameters["b"]
            K = reader(chart["kretschmann"].partition("=")[2])
            if chart["id"] in ("b1_neck", "b2_neck"):
                rho = reader.symbol["\\rho"]
                radius = b / (1 - rho ** 2) if chart["id"] == "b1_neck" else b / (1 + rho ** 2)
            elif chart["id"] in ("a2_cone",):
                radius = reader.symbol["\\sigma"]
            elif chart["id"] == "a2_cartesian":
                radius = reader.parameters["sigma"]
            else:
                radius = reader.parameters["r"] if "r" in reader.parameters else reader.symbol["r"]
            self.same(K, 12 * b ** 2 / radius ** 6)
            self.assertEqual(chart["ricci_scalar"], "R = 0")
            self.assertEqual(chart["ricci_tensor"]["variants"]["ll"]["nonzero"], [])

    def test_a_particle_at_rest_in_a_three_moves_away_from_the_singularity(self):
        entry, reader, _ = self.chart("ab_metrics", "a3")
        gamma = next(e["value"] for e in entry["christoffel"]["variants"]["ull"]["nonzero"]
                     if e["indices"] == ["r", "t", "t"])
        value = reader(gamma).subs({reader.symbol["r"]: 2, reader.parameters["b"]: 1})
        # d^2r/dlambda^2 = -Gamma^r_tt (dt/dlambda)^2 is positive.
        self.assertLess(value, 0)

    def test_kruskals_radius_solves_its_equation(self):
        sp = self.sp
        entry = next(c for c in read("metrics")["coordinates"] if c["id"] == "a2_kruskal")
        reader = self.vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
        r, b = reader.parameters["r"], reader.parameters["b"]
        U, V = reader.symbol["U"], reader.symbol["V"]
        for u, v in ((0.3, 0.7), (-0.5, 0.9), (2.0, 1.5), (-0.2, -0.4)):
            at = {U: u, V: v, b: sp.Rational(3, 2)}
            radius = sp.N(r.subs(at), 30)
            self.assertAlmostEqual(float((radius / 1.5 - 1) * sp.exp(radius / 1.5)), u * v, places=12)

    def test_the_embedded_surface_lies_level_where_the_caption_says(self):
        _, reader, g = self.chart("ab_metrics", "b1_cone")
        sp = self.sp
        tau, r, b = reader.symbol["\\tau"], reader.symbol["r"], reader.parameters["b"]
        # dz/dr = sqrt(g_rr - (d rho/dr)^2) with rho = sqrt(g_phiphi) = r cosh(tau).
        slope = g[("r", "r")] - sp.diff(sp.sqrt(g[("\\phi", "\\phi")]), r) ** 2
        level = b * sp.cosh(tau) ** 2 / sp.sinh(tau) ** 2
        for t in (sp.Rational(2, 5), sp.Rational(4, 5)):
            at = {tau: t, b: sp.Rational(13, 10)}
            self.assertLess(abs(sp.N(slope.subs(r, level).subs(at), 30)), 1e-25)
            self.assertLess(sp.N(slope.subs(r, 2 * level).subs(at), 30), 0)
        self.assertLess(abs(sp.N((slope.subs(tau, 0) - b / (r - b)).subs({r: 3, b: sp.Rational(13, 10)}), 30)), 1e-25)


class Drawings(unittest.TestCase):
    def test_the_neck_is_flamms_paraboloid_at_the_moment_of_symmetry_and_widens_as_cosh(self):
        view = read("embedding")["views"][0]
        frames = view["movie"]["frames"]
        self.assertEqual([f["value"] for f in frames][0], -0.8)
        self.assertEqual(len(frames), 33)
        for frame in frames:
            tau = frame["value"]
            for piece in frame["pieces"]:
                points = piece["points"]
                # Each point is (r, rho, z): the circle of r has the radius r cosh(tau).
                for x, rho, z in points[:: max(1, len(points) // 40)]:
                    self.assertAlmostEqual(rho, x * math.cosh(tau), places=5)
                    if tau == 0:
                        self.assertAlmostEqual(z * z, 4 * (x - 1), places=4)
                self.assertAlmostEqual(min(p[1] for p in points), math.cosh(tau), places=6)
                if tau:
                    level = math.cosh(tau) ** 2 / math.sinh(tau) ** 2
                    self.assertLessEqual(max(p[0] for p in points), level + 1e-9)

    def test_the_moments_of_the_neck_are_marked_on_b_one_alone(self):
        systems = read("diagrams")["systems"]
        self.assertEqual(sum(len(views) for views in systems.values()), 16)
        for system, views in systems.items():
            for view in views:
                slices = view.get("slices", [])
                self.assertEqual(len(slices), 5 if system.startswith("b1") else 0, (system, view["id"]))
        for view in read("conformal")["views"]:
            self.assertEqual(len(view.get("slices", [])), 5 if view["id"].startswith("b1") else 0, view["id"])

    def test_a_twos_singularities_stand_at_the_sides_of_its_conformal_diagram(self):
        views = {v["id"]: v for v in read("conformal")["views"]}
        self.assertEqual(set(views), {"a2_static", "a2_cone", "a2_inertial", "a2_kruskal",
                                      "b1_cone", "b1_neck", "b1_static", "b1_inertial"})
        singular = [layer for layer in views["a2_kruskal"]["layers"] if layer.get("class") == "singular"]
        self.assertTrue(singular)
        for layer in singular:
            self.assertEqual({round(abs(x), 3) for x, _ in layer["points"]}, {round(math.pi / 2, 3)})
        # BI's plane has no singular line.
        self.assertFalse([layer for layer in views["b1_neck"]["layers"] if layer.get("class") == "singular"])


if __name__ == "__main__":
    unittest.main()
