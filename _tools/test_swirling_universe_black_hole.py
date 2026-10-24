"""Schwarzschild's black hole in a swirling universe, Astorino, Martelli and Vigano's solution:
.venv.noindex/bin/python -m unittest discover -s _tools

What its texts and drawings state, held on the published files. The first class needs nothing
but the files. The second needs sympy and numpy and is skipped where they are absent: it reads
the published metric of the spherical chart and holds it to the horizon's area 16 pi m^2 and
temperature 1/(8 pi m) for every j, to a horizon at rest in the azimuth, to a regular axis, to no
closed timelike curve, to the vanishing angular momentum of Moreira, Herdeiro and Crispino's Komar
integral, to the horizon's equatorial circumference and Gaussian curvature, and to the widest
circle of the equator at r = 1/sqrt(j).
"""
import importlib.util
import json
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
HAS_NUMBERS = all(importlib.util.find_spec(name) is not None for name in ("sympy", "numpy"))
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))
ME = "swirling_universe_black_hole"


def load(folder):
    return json.loads((DATA / folder / f"{ME}.json").read_text(encoding="utf-8"))


class Drawings(unittest.TestCase):
    def test_each_chart_draws_its_plane_and_its_ergoregions(self):
        diagrams = load("diagrams")
        for system in ("spherical", "cylindrical"):
            self.assertEqual([v["id"] for v in diagrams["systems"][system]], ["radial"])
            self.assertEqual([f["id"] for f in diagrams["projections"][system]], ["ergo"])
        self.assertEqual({v["id"] for v in load("conformal")["views"]}, {"spherical", "cylindrical"})
        self.assertEqual([v["id"] for v in load("embedding")["views"]], ["universe", "hole"])
        self.assertEqual([v["id"] for v in load("curvature")["views"]], ["universe", "hole"])

    def test_the_ergoregions_are_two_and_clear_of_the_hole(self):
        """The figure of the hole fills one region above the equatorial plane and its mirror image
        below, each on both sides of the axis, none nearer the centre than r = 3.67 m, and covers the
        hole; the universe's are the same two with no hole."""
        for system, nearest in (("spherical", 3.66), ("cylindrical", 0.43)):
            figure, = load("diagrams")["projections"][system]
            fills = [layer for layer in figure["layers"] if layer["kind"] == "fill"]
            self.assertEqual([layer["class"] for layer in fills].count("double"), 4)
            self.assertEqual([layer["class"] for layer in fills].count("wedge"), 1 if system == "spherical" else 0)
            self.assertNotIn("turn", figure)
        figure, = load("diagrams")["projections"]["spherical"]
        caption = " ".join(figure["caption"])
        self.assertIn("$r = 3.67\\,m$", caption)

    def test_the_equator_widens_to_its_widest_circle_and_narrows(self):
        """The circle of r on the equator has radius r/sqrt(1 + j^2 r^4): 2m/sqrt(1 + 1/16) at the
        throat r = 2m, greatest, 2 sqrt 2 m, at r = 4m, and smaller beyond, out to r = 6m."""
        view = next(v for v in load("embedding")["views"] if v["id"] == "hole")
        points = view["surfaces"][0]["pieces"][0]["points"]
        self.assertAlmostEqual(points[0][0], 2.0, places=6)
        self.assertAlmostEqual(points[0][1], 2 / math.sqrt(1 + 1 / 16), places=5)
        radius, at = max((p[1], p[0]) for p in points)
        self.assertAlmostEqual(radius, 2 * math.sqrt(2), places=4)
        self.assertAlmostEqual(at, 4.0, delta=0.05)
        self.assertLess(points[-1][1], radius)
        self.assertAlmostEqual(points[-1][0], 6.0, places=6)


@unittest.skipUnless(HAS_NUMBERS, "needs sympy and numpy")
class Physics(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import sympy as sp
        import verify_metrics as vm
        cls.sp = sp
        metric = json.loads((DATA / "metrics" / f"{ME}.json").read_text(encoding="utf-8"))
        chart = next(c for c in metric["coordinates"] if c["id"] == "spherical")
        reader = vm.Reader(chart["coords"], [p["symbol"] for p in chart["parameters"]], ())
        cls.t, cls.r, cls.th, cls.ph = (reader.symbol[c] for c in chart["coords"])
        cls.m, cls.j = reader.parameters["m"], reader.parameters["j"]
        values = {tuple(e["indices"]): e["value"] for e in chart["metric_components"]}
        names = chart["coords"]
        cls.g = sp.Matrix(4, 4, lambda a, b: reader(values.get((names[a], names[b]), "0")))

    def test_the_horizon_keeps_schwarzschilds_area_for_every_j(self):
        sp = self.sp
        g = self.g.subs(self.r, 2 * self.m)
        area = sp.integrate(sp.simplify(sp.sqrt(g[2, 2] * g[3, 3])), (self.th, 0, sp.pi)) * 2 * sp.pi
        self.assertEqual(sp.simplify(area - 16 * sp.pi * self.m ** 2), 0)

    def test_the_horizon_is_at_rest_and_its_temperature_is_schwarzschilds(self):
        """The dragging rate -g_tphi/g_phiphi vanishes at r = 2m at every theta, so d_t is the null
        generator of the horizon, and g_tphi vanishes there as r - 2m does; the surface gravity is then
        (1/2) d_r(-g_tt)/sqrt(-g_tt g_rr) on the horizon, c^2/4m at every theta, so T = 1/(8 pi m) in
        units with G = c = 1."""
        sp = self.sp
        g = self.g
        drag = sp.simplify(-g[0, 3] / g[3, 3])
        self.assertEqual(sp.simplify(drag.subs(self.r, 2 * self.m)), 0)
        f = -g[0, 0]
        kappa = sp.limit(sp.diff(f, self.r) / (2 * sp.sqrt(f * g[1, 1])), self.r, 2 * self.m, "+")
        self.assertEqual(sp.simplify(kappa.subs(self.th, sp.pi / 3) - 1 / (4 * self.m)), 0)
        self.assertEqual(sp.simplify(kappa.subs(self.th, sp.pi / 7) - 1 / (4 * self.m)), 0)

    def test_the_axis_is_regular_and_no_circle_is_timelike(self):
        sp = self.sp
        g = self.g
        ratio = sp.simplify(g[3, 3] / g[2, 2])
        for pole in (0, sp.pi):
            self.assertEqual(sp.limit(ratio / sp.sin(self.th) ** 2, self.th, pole), 1)
        import numpy as np
        fn = sp.lambdify((self.r, self.th, self.m, self.j), g[3, 3], "numpy")
        R, T = np.meshgrid(np.linspace(0.05, 40, 200), np.linspace(0.01, math.pi - 0.01, 200))
        self.assertTrue(np.all(fn(R, T, 1.0, 1 / 16) > 0) and np.all(fn(R, T, 1.0, 4.0) > 0))

    def test_the_hemispheres_turn_oppositely_and_the_komar_angular_momentum_vanishes(self):
        """Moreira, Herdeiro and Crispino's star d psi^flat, 4 j r^4 sin^3 theta cos theta (r f' + f)/F^2,
        integrates to zero over a sphere, being odd about the equator."""
        sp = self.sp
        r, th, m, j = self.r, self.th, self.m, self.j
        drag = sp.simplify(-self.g[0, 3] / self.g[3, 3])
        self.assertEqual(sp.simplify(drag.subs(th, sp.pi - th) + drag), 0)
        f = 1 - 2 * m / r
        F = 1 + j ** 2 * r ** 4 * sp.sin(th) ** 4
        density = 4 * j * r ** 4 * sp.sin(th) ** 3 * sp.cos(th) * (r * sp.diff(f, r) + f) / F ** 2
        self.assertEqual(sp.simplify(density.subs(th, sp.pi - th) + density), 0)

    def test_the_horizon_narrows_at_its_equator_and_curves_negatively_past_jm2_one_quarter(self):
        """The equatorial circumference of the horizon is 4 pi m/sqrt(1 + 16 j^2 m^4), and its Gaussian
        curvature on the equator (1 - 16 j^2 m^4)/(4 m^2 (1 + 16 j^2 m^4)^2), Moreira, Herdeiro and
        Crispino's, negative once j m^2 > 1/4."""
        sp = self.sp
        m, j, th = self.m, self.j, self.th
        g = self.g.subs(self.r, 2 * m)
        # The circumference squared, (2 pi)^2 g_phiphi, against (4 pi m)^2/(1 + 16 j^2 m^4).
        self.assertEqual(sp.simplify(g[3, 3].subs(th, sp.pi / 2) - 4 * m ** 2 / (1 + 16 * j ** 2 * m ** 4)), 0)
        a, b = g[2, 2], g[3, 3]
        root = sp.sqrt(a * b)
        gauss = -sp.diff(sp.diff(b, th) / root, th) / (2 * root)
        at = sp.simplify(gauss.subs(th, sp.pi / 2))
        self.assertEqual(sp.simplify(at - (1 - 16 * j ** 2 * m ** 4) / (4 * m ** 2 * (1 + 16 * j ** 2 * m ** 4) ** 2)), 0)

    def test_the_widest_circle_of_the_equator_is_at_one_over_root_j(self):
        sp = self.sp
        radius = sp.sqrt(self.g[3, 3].subs(self.th, sp.pi / 2))
        rp = sp.symbols("rp", positive=True)
        jp = sp.symbols("jp", positive=True)
        stationary = sp.solve(sp.diff(radius.subs({self.r: rp, self.j: jp}), rp), rp)
        self.assertEqual(stationary, [1 / sp.sqrt(jp)])
