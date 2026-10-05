"""Bonnor's magnetic dipole, the static field of a mass with a magnetic dipole moment:
.venv.noindex/bin/python -m unittest discover -s _tools

What its texts and drawings state, held on the published files. The first class needs nothing but
the files and holds the drawings, which draw the dipole at m = 1 and b = 2 sqrt 2, where
k = sqrt(m^2 + b^2) = 3 and the axis between the two black holes is r = 4m. The second needs sympy
and is skipped where it is absent: it reads the published metric components and holds them to the
mass 2m, to Emparan's excess of angle (1 + m^2/b^2)^2 on the axis between the holes and to a
regular axis beyond them, to the tension of the strings that would replace the strut, to Zipoy
and Voorhees's published metric at b = 0, and the Kretschmann scalar to a finite value at a hole.
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

EXCESS = 81 / 64        # (1 + m^2/b^2)^2 at b^2 = 8 m^2
EDGE = 4.0              # r = m + sqrt(m^2 + b^2), the axis between the holes, in units of m


def load(folder, metric="bonnor_magnetic_dipole"):
    return json.loads((DATA / folder / f"{metric}.json").read_text(encoding="utf-8"))


class Drawings(unittest.TestCase):
    def test_the_tip_of_the_embedded_plane_is_a_cone_of_81_64_of_a_turn(self):
        """At the strut a circle the metric distance l from the axis has circumference
        (1 + m^2/b^2)^2 2 pi l. The tip is drawn in Minkowski space, where the distance along the
        profile is sqrt(d rho^2 - d z^2)."""
        view = next(v for v in load("embedding")["views"] if v["id"] == "equator")
        tip, plane = view["surfaces"][0]["pieces"]
        self.assertEqual(tip.get("space"), "minkowski")
        self.assertNotIn("space", plane)
        self.assertEqual(tip["points"][0], [EDGE, 0.0, 0.0])
        _, rho, z = tip["points"][1]
        self.assertAlmostEqual(rho / math.sqrt(rho ** 2 - z ** 2), EXCESS, delta=2e-3)
        for a, b in zip(tip["points"][-1], plane["points"][0]):
            self.assertAlmostEqual(a, b, places=6)

    def test_the_circles_of_the_embedded_plane_are_r_root_Z_over_r_minus_2m(self):
        view = next(v for v in load("embedding")["views"] if v["id"] == "equator")
        for piece in view["surfaces"][0]["pieces"]:
            for r, rho, _ in piece["points"]:
                self.assertAlmostEqual(rho, r * math.sqrt(max((r - 4) * (r + 2), 0.0)) / (r - 2), delta=2e-4)

    def test_the_axis_between_the_holes_is_a_whole_diamond_between_two_horizons(self):
        """On the conformal diagram of the axis between the holes all four edges are horizons, and the
        moment embedded is the one event t = 0 midway, the centre of the diamond."""
        for system in ("spheroidal",):
            view = next(v for v in load("conformal")["views"] if v["id"] == f"{system}_strut")
            horizons = [layer for layer in view["layers"] if layer.get("class") == "horizon"]
            self.assertEqual(len(horizons), 4)
            self.assertFalse([layer for layer in view["layers"] if layer.get("class") in ("scri", "singular")])
            (mark,) = view["slices"]
            self.assertEqual(mark["lines"], [])
            (at,) = mark["points"]
            self.assertAlmostEqual(at[0], 0.0, places=9)
            self.assertAlmostEqual(at[1], 0.0, places=9)

    def test_the_axis_beyond_a_hole_ends_on_a_horizon_and_the_equator_on_the_strut(self):
        for system in ("spheroidal",):
            views = {v["id"]: v for v in load("conformal")["views"]}
            classes = lambda view: {layer.get("class") for layer in view["layers"]}  # noqa: E731
            self.assertIn("horizon", classes(views[f"{system}_axis"]))
            self.assertIn("scri", classes(views[f"{system}_axis"]))
            self.assertIn("centre", classes(views[f"{system}_equator"]))
            self.assertNotIn("horizon", classes(views[f"{system}_equator"]))

    def test_every_chart_draws_its_three_planes(self):
        systems = load("diagrams")["systems"]
        for system in ("spheroidal",):
            self.assertEqual([view["id"] for view in systems[system]], ["axis", "equator", "strut"])
            for view in systems[system]:
                self.assertGreater(len(view["cones"]), 20)
                self.assertFalse([m for m in view["markers"] if m["kind"] == "singular"])

    def test_the_rays_of_the_strut_take_longer_and_longer_toward_a_hole(self):
        """c dt/dtheta = 64 (1 + sin^2 theta)^2/(27 sin^3 theta) on the axis between the holes: the
        cones are narrowest in theta midway and close toward each end."""
        view = next(v for v in load("diagrams")["systems"]["spheroidal"] if v["id"] == "strut")
        X0, X1, Y0, Y1 = view["box"]

        for cone in view["cones"]:
            theta = X0 + cone["at"][0] * (X1 - X0)
            want = 27 * math.sin(theta) ** 3 / (64 * (1 + math.sin(theta) ** 2) ** 2)
            for edge in (cone["a"], cone["b"]):
                # The edge's direction in the unit square, as d theta over d(ct).
                got = abs(edge[0] * (X1 - X0)) / abs(edge[1] * (Y1 - Y0))
                self.assertAlmostEqual(got, want, delta=2e-3 + 0.02 * want)
        widths = sorted((abs(cone["at"][0] - 0.5), abs(cone["a"][0])) for cone in view["cones"])
        self.assertGreater(widths[0][1], 10 * widths[-1][1])


@unittest.skipUnless(HAS_SYMPY, "needs sympy")
class Physics(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import sympy as sp
        import verify_metrics as vm
        cls.sp, cls.vm = sp, vm
        cls.charts = {}
        for chart in load("metrics")["coordinates"]:
            reader = vm.Reader(chart["coords"], [p["symbol"] for p in chart["parameters"]], ())
            g = {tuple(e["indices"]): reader(e["value"]) for e in chart["metric_components"]}
            cls.charts[chart["id"]] = (chart, reader, g)

    def at(self, value, **numbers):
        _, reader, _ = self.charts["spheroidal"]
        names = {"r": reader.symbol["r"], "theta": reader.symbol["\\theta"],
                 "m": reader.parameters["m"], "b": reader.parameters["b"]}
        return float(value.subs({names[k]: self.sp.Rational(str(v)) for k, v in numbers.items()}))

    def test_the_mass_is_2m(self):
        """g_tt -> -(1 - 4m/r): the mass M has GM/c^2 = 2m."""
        _, reader, g = self.charts["spheroidal"]
        r = reader.symbol["r"]
        for m, b, theta in ((1, 2, "0.3"), (2, 1, "1.1"), (1, 5, "1.5")):
            far = self.at((1 + g[("t", "t")]) * r, r=10 ** 7, theta=theta, m=m, b=b)
            self.assertAlmostEqual(far, 4 * m, delta=1e-4)

    def test_the_axis_between_the_holes_has_emparans_excess_of_angle(self):
        """A small circle about r = r_+ has circumference (1 + m^2/b^2)^2 times 2 pi times its
        radius: d sqrt(g_phiphi)/(sqrt(g_rr) dr) -> (1 + m^2/b^2)^2 there, for any theta."""
        sp = self.sp
        _, reader, g = self.charts["spheroidal"]
        r, m, b = reader.symbol["r"], reader.parameters["m"], reader.parameters["b"]
        ratio = sp.diff(sp.sqrt(g[("\\phi", "\\phi")]), r) / sp.sqrt(g[("r", "r")])
        for mv, bv, theta in ((1, sp.sqrt(8), 1.0), (1, 2, 0.4), (3, 1, 2.0)):
            edge = mv + sp.sqrt(mv ** 2 + bv ** 2)
            near = ratio.subs({m: mv, b: bv, reader.symbol["\\theta"]: theta}).subs(r, edge + sp.Float("1e-12", 40))
            self.assertAlmostEqual(float(sp.N(near, 30)), float((1 + mv ** 2 / bv ** 2) ** 2), delta=1e-5)

    def test_the_axis_beyond_the_holes_is_regular(self):
        """d sqrt(g_phiphi)/(sqrt(g_thetatheta) dtheta) -> 1 at theta = 0: no string and no strut
        beyond the holes when phi runs over 2 pi."""
        sp = self.sp
        _, reader, g = self.charts["spheroidal"]
        th, m, b, r = reader.symbol["\\theta"], reader.parameters["m"], reader.parameters["b"], reader.symbol["r"]
        ratio = sp.diff(sp.sqrt(g[("\\phi", "\\phi")]), th) / sp.sqrt(g[("\\theta", "\\theta")])
        for mv, bv, rv in ((1, 2, 6), (2, 1, 9), (1, 4, 6)):
            near = ratio.subs({m: mv, b: bv, r: rv}).subs(th, sp.Float("1e-10", 40))
            self.assertAlmostEqual(float(sp.N(near, 30)), 1.0, delta=1e-8)

    def test_the_tension_of_the_strings_that_would_replace_the_strut(self):
        """With the period of phi chosen to smooth the axis between the holes, 2 pi/(1 + m^2/b^2)^2,
        the deficit beyond them is 8 pi T with T = (1 - b^4/(m^2 + b^2)^2)/4, which is m^2/(2 b^2)
        for b much larger than m."""
        for m, b in ((1.0, 2.0), (1.0, math.sqrt(8)), (3.0, 1.0)):
            period = 2 * math.pi / (1 + m ** 2 / b ** 2) ** 2
            tension = (2 * math.pi - period) / (8 * math.pi)
            self.assertAlmostEqual(tension, (1 - b ** 4 / (m ** 2 + b ** 2) ** 2) / 4, places=12)
        m, b = 1.0, 1000.0
        self.assertAlmostEqual((1 - b ** 4 / (m ** 2 + b ** 2) ** 2) / 4 / (m ** 2 / (2 * b ** 2)), 1.0, delta=1e-5)

    def test_without_the_dipole_it_is_zipoy_and_voorhees_at_q_1(self):
        """At b = 0 the published metric is the published spherical chart of Zipoy and Voorhees's
        deformed mass at q = 1, which is not Schwarzschild's."""
        sp, vm = self.sp, self.vm
        _, reader, g = self.charts["spheroidal"]
        zv = next(c for c in load("metrics", "zipoy_voorhees")["coordinates"] if c["id"] == "spherical")
        other = vm.Reader(zv["coords"], [p["symbol"] for p in zv["parameters"]], ())
        same = {other.symbol[n]: reader.symbol[n] for n in zv["coords"]}
        same.update({other.parameters["m"]: reader.parameters["m"], other.parameters["q"]: 1})
        for entry in zv["metric_components"]:
            theirs = other(entry["value"]).subs(same, simultaneous=True)
            ours = g[tuple(entry["indices"])].subs(reader.parameters["b"], 0)
            self.assertEqual(vm.norm(theirs - ours), 0)
        r, m = reader.symbol["r"], reader.parameters["m"]
        self.assertEqual(sp.simplify(g[("t", "t")].subs(reader.parameters["b"], 0) + (1 - 2 * m / r) ** 2), 0)

    def test_the_curvature_is_finite_at_a_hole(self):
        """The Kretschmann scalar tends to a finite value at a hole, along the axis between the
        holes and along the axis beyond it, though the two values differ: each hole is a whole
        horizon, which Bonnor's coordinates draw as one point."""
        sp = self.sp
        own, reader, _ = self.charts["spheroidal"]
        K = reader(own["kretschmann"].split("=", 1)[1])
        m, b, r, th = reader.parameters["m"], reader.parameters["b"], reader.symbol["r"], reader.symbol["\\theta"]
        at = K.subs({m: 1, b: sp.sqrt(8)})
        on_strut = sp.limit(sp.simplify(at.subs(r, 4)), th, 0)
        beyond = sp.limit(sp.simplify(at.subs(th, 0)), r, 4, "+")
        for value in (on_strut, beyond):
            self.assertTrue(value.is_finite)
            self.assertGreater(float(value), 0)
        self.assertNotAlmostEqual(float(on_strut), float(beyond), places=3)


if __name__ == "__main__":
    unittest.main()
