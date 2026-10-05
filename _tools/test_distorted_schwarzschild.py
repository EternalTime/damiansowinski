"""Schwarzschild's black hole in a quadrupole tidal field, the simplest of Geroch and Hartle's
distorted black holes: .venv.noindex/bin/python -m unittest discover -s _tools

What its texts and drawings state, held on the published files and on the two shapes every drawing
draws, the oblate q = 1/12 and the prolate q = -1/12 at m = 1. The first class needs nothing but
the files. The second needs numpy, scipy and sympy and is skipped where they are absent: it holds
the published U and V to Weyl's field equations, to Doroshkevich, Zel'dovich and Novikov's
Appendix IV and to Schwarzschild's metric at q = 0, the horizon to the surface gravity
e^{2q}/4m at every latitude, to the area 16 pi m^2 e^{-2q}, to the mass m and to the Gaussian
curvature the History and the captions quote, the published Kretschmann scalar to the values the
captions state on the horizon, and the tortoise coordinate the drawings use to the published
metric.
"""
import importlib.util
import json
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
HAS_NUMBERS = all(importlib.util.find_spec(name) is not None for name in ("sympy", "numpy", "scipy", "contourpy"))
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))
ME = "distorted_schwarzschild"
SHAPES = {"oblate": 1 / 12, "prolate": -1 / 12}
CHARTS = ("prolate_spheroidal", "spherical", "weyl")


def load(folder):
    return json.loads((DATA / folder / f"{ME}.json").read_text(encoding="utf-8"))


class Drawings(unittest.TestCase):
    def test_each_chart_draws_both_planes_of_both_shapes(self):
        systems = load("diagrams")["systems"]
        want = {f"{plane}_{shape}" for plane in ("axis", "equator") for shape in SHAPES}
        for system in CHARTS:
            self.assertEqual({v["id"] for v in systems[system]}, want)
        conformal = {v["id"] for v in load("conformal")["views"]}
        self.assertEqual(conformal, {f"{s}_{w}" for s in CHARTS for w in want})

    def test_only_the_oblate_axis_ends_on_a_timelike_edge(self):
        """On the axis at q > 0 the tortoise coordinate is finite at infinity, so that plane alone
        has a timelike edge where the others have null infinity; every plane has the horizon and
        the singularity."""
        for view in load("conformal")["views"]:
            classes = {layer.get("class") for layer in view["layers"]}
            timelike = view["id"].endswith("axis_oblate")
            self.assertEqual("boundary" in classes, timelike, view["id"])
            self.assertEqual("scri" in classes, not timelike, view["id"])
            self.assertTrue({"horizon", "singular"} <= classes, view["id"])

    def test_the_horizon_is_flattened_by_a_ring_and_stretched_by_masses_on_the_axis(self):
        """The equator of the horizon has radius 2m e^{q/2}, and a pole stands 1.225 m from the
        centre at q = 1/12 and 2.652 m at q = -1/12; the area is 16 pi m^2 e^{-2q}."""
        views = {v["id"]: v for v in load("embedding")["views"]}
        for shape, q in SHAPES.items():
            points = views[f"horizon_{shape}"]["surfaces"][0]["pieces"][0]["points"]
            equator = max(p[1] for p in points)
            half = (points[-1][2] - points[0][2]) / 2
            self.assertAlmostEqual(equator, 2 * math.exp(q / 2), places=5)
            self.assertAlmostEqual(half, {"oblate": 1.225, "prolate": 2.652}[shape], places=3)
            self.assertEqual(equator > half, shape == "oblate")
            area = sum(math.pi * (a[1] + b[1]) * math.hypot(b[1] - a[1], b[2] - a[2]) for a, b in zip(points, points[1:]))
            self.assertAlmostEqual(area / (16 * math.pi * math.exp(-2 * q)), 1.0, places=3)

    def test_the_equatorial_plane_meets_the_horizon_on_its_equator(self):
        """The throat of each equatorial plane is the horizon's equator, of radius 2m e^{q/2}, and
        the oblate plane lies level at r = 2.54 m, on a circle of radius 2.883 m, where it passes
        into Minkowski space."""
        views = {v["id"]: v for v in load("embedding")["views"]}
        for shape, q in SHAPES.items():
            pieces = views[shape]["surfaces"][0]["pieces"]
            self.assertAlmostEqual(pieces[0]["points"][0][0], 2.0)
            self.assertAlmostEqual(pieces[0]["points"][0][1], 2 * math.exp(q / 2), places=5)
        first, second = views["oblate"]["surfaces"][0]["pieces"][:2]
        self.assertAlmostEqual(first["points"][-1][0], 2.5388, places=3)
        self.assertAlmostEqual(first["points"][-1][1], 2.883, places=3)
        for a, b in zip(second["points"][0], first["points"][-1]):
            self.assertAlmostEqual(a, b, places=6)
        self.assertEqual(len(views["prolate"]["surfaces"][0]["pieces"]), 2)

    def test_the_spacetime_diagrams_are_square(self):
        for views in load("diagrams")["systems"].values():
            for view in views:
                X0, X1, Y0, Y1 = view["box"]
                self.assertAlmostEqual((X1 - X0) / (Y1 - Y0), 1.0)


@unittest.skipUnless(HAS_NUMBERS, "needs numpy, scipy and sympy")
class Physics(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import numpy as np
        import sympy as sp

        import null_rays as nr
        import verify_metrics as vm
        cls.np, cls.sp, cls.nr, cls.vm = np, sp, nr, vm
        cls.charts = {}
        for entry in load("metrics")["coordinates"]:
            reader = vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], (),
                               held=vm.HELD[(ME, entry["id"])])
            cls.charts[entry["id"]] = (entry, reader)
        _, r = cls.charts["prolate_spheroidal"]
        cls.x, cls.y = r.symbol["x"], r.symbol["y"]
        cls.m, cls.q = r.parameters["m"], r.parameters["q"]
        cls.U, cls.V = r.held[r.parameters["U"]], r.held[r.parameters["V"]]

    def test_the_published_functions_solve_weyls_equations(self):
        sp, x, y, U, V = self.sp, self.x, self.y, self.U, self.V
        self.assertEqual(sp.expand(sp.diff((x ** 2 - 1) * sp.diff(U, x), x) + sp.diff((1 - y ** 2) * sp.diff(U, y), y)), 0)
        # Weyl's quadrature for gamma = gamma_S + V with psi = psi_S + U.
        P = sp.log((x - 1) / (x + 1)) / 2 + U
        G = sp.log((x ** 2 - 1) / (x ** 2 - y ** 2)) / 2 + V
        px, py = sp.diff(P, x), sp.diff(P, y)
        along = (1 - y ** 2) / (x ** 2 - y ** 2) * (x * (x ** 2 - 1) * px ** 2 - x * (1 - y ** 2) * py ** 2
                                                      - 2 * y * (x ** 2 - 1) * px * py)
        across = (x ** 2 - 1) / (x ** 2 - y ** 2) * (y * (x ** 2 - 1) * px ** 2 - y * (1 - y ** 2) * py ** 2
                                                       + 2 * x * (1 - y ** 2) * px * py)
        self.assertEqual(sp.simplify(sp.diff(G, x) - along), 0)
        self.assertEqual(sp.simplify(sp.diff(G, y) - across), 0)
        self.assertEqual(sp.expand(V.subs(y, 1)), 0)
        self.assertEqual(sp.expand(V.subs(y, -1)), 0)

    def test_they_are_doroshkevich_zeldovich_and_novikovs(self):
        """Appendix IV: psi = (1/2) ln((l - 1)/(l + 1)) + (q/4)(3 l^2 - 1)(3 mu^2 - 1) and gamma =
        (1/2) ln((l^2 - 1)/(l^2 - mu^2)) - 3 q l (1 - mu^2) - (9/16) q^2 (l^2 - 1)(1 - mu^2)
        (9 mu^2 l^2 - l^2 - mu^2 + 1)."""
        sp, x, y, q = self.sp, self.x, self.y, self.q
        self.assertEqual(sp.expand(self.U - q * (3 * x ** 2 - 1) * (3 * y ** 2 - 1) / 4), 0)
        printed = (-3 * q * x * (1 - y ** 2)
                   - sp.Rational(9, 16) * q ** 2 * (x ** 2 - 1) * (1 - y ** 2) * (9 * y ** 2 * x ** 2 - x ** 2 - y ** 2 + 1))
        self.assertEqual(sp.expand(self.V - printed), 0)

    def test_without_the_tidal_field_it_is_schwarzschild(self):
        self.assertEqual(self.U.subs(self.q, 0), 0)
        self.assertEqual(self.V.subs(self.q, 0), 0)

    def test_the_spherical_and_weyl_charts_hold_the_same_functions(self):
        sp = self.sp
        _, r = self.charts["spherical"]
        rr, th, m, q = r.symbol["r"], r.symbol["\\theta"], r.parameters["m"], r.parameters["q"]
        at = {self.x: rr / m - 1, self.y: sp.cos(th), self.q: q, self.m: m}
        for name, mine in (("U", self.U), ("V", self.V)):
            self.assertEqual(sp.simplify(r.held[r.parameters[name]] - mine.subs(at, simultaneous=True)), 0)
        _, w = self.charts["weyl"]
        rho, z = w.symbol["\\rho"], w.symbol["z"]
        psi, gamma = w.held[w.parameters["psi"]], w.held[w.parameters["gamma"]]
        for x0, y0, q0 in ((sp.Rational(7, 3), sp.Rational(2, 5), sp.Rational(1, 12)),
                           (sp.Rational(3, 2), sp.Rational(-3, 4), sp.Rational(-1, 12))):
            here = {self.x: x0, self.y: y0, self.q: q0, self.m: 1}
            there = {rho: sp.sqrt((x0 ** 2 - 1) * (1 - y0 ** 2)), z: x0 * y0, w.parameters["m"]: 1, w.parameters["q"]: q0}
            want_psi = sp.log((x0 - 1) / (x0 + 1)) / 2 + self.U.subs(here)
            want_gamma = sp.log((x0 ** 2 - 1) / (x0 ** 2 - y0 ** 2)) / 2 + self.V.subs(here)
            self.assertAlmostEqual(float(sp.N(psi.subs(there), 30)), float(want_psi), places=12)
            self.assertAlmostEqual(float(sp.N(gamma.subs(there), 30)), float(want_gamma), places=12)

    def published(self, chart):
        """The published metric of a chart as a sympy matrix with U and V written out."""
        sp = self.sp
        entry, r = self.charts[chart]
        out = sp.zeros(4, 4)
        for component in entry["metric_components"]:
            i, j = (entry["coords"].index(n) for n in component["indices"])
            out[i, j] = out[j, i] = r(component["value"])
        written = {name: r.held[name] for name in r.held}
        return out.subs(written), r

    def test_the_surface_gravity_is_e_to_the_2q_over_4m_at_every_latitude(self):
        """From the published spherical metric: kappa^2 = lim (d_r g_tt)^2/(-4 g_tt g_rr) at r = 2m,
        which is e^{4q}/16m^2 whatever theta is, Frolov and Shoom's e^{2 u_0}/4m with u_0 = q, the
        value of U at the poles of the horizon. The exponent is +2q."""
        sp = self.sp
        g, r = self.published("spherical")
        rr, th, m, q = r.symbol["r"], r.symbol["\\theta"], r.parameters["m"], r.parameters["q"]
        squared = sp.simplify(sp.diff(g[0, 0], rr) ** 2 / (-4 * g[0, 0] * g[1, 1]))
        at_horizon = sp.simplify(squared.subs(rr, 2 * m))
        self.assertEqual(sp.simplify(at_horizon - sp.exp(4 * q) / (16 * m ** 2)), 0)
        self.assertFalse(at_horizon.has(th))
        U = r.held[r.parameters["U"]]
        self.assertEqual(sp.simplify(U.subs({rr: 2 * m, th: 0}) - q), 0)

    def test_the_horizon_has_the_area_the_mass_and_the_curvature_the_texts_state(self):
        """The area is 16 pi m^2 e^{-2q}, so kappa A/4 pi = m, and the Gaussian curvature is
        (1 - 12q) e^{2q}/4m^2 at a pole and (1 + 3q) e^{5q}/4m^2 on the equator: negative at the
        poles for q > 1/12, Frolov and Shoom's bound."""
        sp = self.sp
        g, r = self.published("spherical")
        rr, th, m, q = r.symbol["r"], r.symbol["\\theta"], r.parameters["m"], r.parameters["q"]
        E, G = sp.simplify(g[2, 2].subs(rr, 2 * m)), sp.simplify(g[3, 3].subs(rr, 2 * m))
        root_g = 2 * m * sp.exp(-q * (3 * sp.cos(th) ** 2 - 1) / 2) * sp.sin(th)
        self.assertEqual(sp.simplify(G - root_g ** 2), 0)
        element = sp.simplify(sp.sqrt(sp.simplify(E / (4 * m ** 2))) * 2 * m * root_g)
        area = 2 * sp.pi * sp.integrate(sp.simplify(element), (th, 0, sp.pi))
        self.assertEqual(sp.simplify(area - 16 * sp.pi * m ** 2 * sp.exp(-2 * q)), 0)
        kappa = sp.exp(2 * q) / (4 * m)
        self.assertEqual(sp.simplify(kappa * area / (4 * sp.pi) - m), 0)
        root_e = 2 * m * sp.exp(sp.simplify(sp.log(E / (4 * m ** 2)) / 2))
        curvature = sp.simplify(-sp.diff(sp.diff(root_g, th) / root_e, th) / (root_e * root_g))
        pole = sp.limit(curvature, th, 0)
        self.assertEqual(sp.simplify(pole - (1 - 12 * q) * sp.exp(2 * q) / (4 * m ** 2)), 0)
        equator = curvature.subs(th, sp.pi / 2)
        self.assertEqual(sp.simplify(equator - (1 + 3 * q) * sp.exp(5 * q) / (4 * m ** 2)), 0)

    def test_the_kretschmann_scalar_on_the_horizon_is_twelve_times_the_gaussian_curvature_squared(self):
        """The published scalar next to the horizon, on the axis and on the equator, against
        (3/4)(1 - 12q)^2 e^{4q}/m^4 and (3/4)(1 + 3q)^2 e^{10q}/m^4: zero and 2.70 at q = 1/12,
        3 e^{-1/3} and 0.183 at q = -1/12, the numbers the captions give."""
        import conformal as cf
        src = cf.Sources()
        stated = {("oblate", "axis"): 0.0, ("oblate", "equator"): 2.70,
                  ("prolate", "axis"): 3 * math.exp(-1 / 3), ("prolate", "equator"): 0.183}
        for shape, q in SHAPES.items():
            for plane, fixed, off in (("axis", {"theta": "0", "phi": "0"}, {"theta": "1e-30"}),
                                      ("equator", {"theta": "pi/2", "phi": "0"}, None)):
                pl = cf.Plane(src, ME, "spherical", ("t", "r"), fixed, {"m": 1, "q": "1/12" if q > 0 else "-1/12"},
                              kretschmann_fixed=off)
                want = (0.75 * (1 - 12 * q) ** 2 * math.exp(4 * q) if plane == "axis"
                        else 0.75 * (1 + 3 * q) ** 2 * math.exp(10 * q))
                got = float(pl.kretschmann(0, 2 + 1e-7))
                self.assertAlmostEqual(got, want, places=4)
                self.assertAlmostEqual(got, stated[(shape, plane)], places=2)

    def test_the_tortoise_coordinate_is_the_published_metrics(self):
        """dx_*/dx = sqrt(-g_xx/g_tt)/m from the published prolate spheroidal metric on each plane,
        the residue 2 e^{-2q} at the horizon on both, x_* = 0 at the singularity, and on the oblate
        axis the finite limit R = 0.3443 m far out."""
        sp, nr = self.sp, self.nr
        g, r = self.published("prolate_spheroidal")
        ratio = sp.sqrt(sp.simplify(-g[1, 1] / g[0, 0])) / self.m
        for shape, q in SHAPES.items():
            for plane, y in (("axis", 1), ("equator", 0)):
                star = nr._distorted_star(plane, q)
                for x in (-0.5, 0.4, 1.3, 2.5):
                    want = float(sp.Abs(ratio.subs({self.y: y, self.q: sp.nsimplify(q), self.x: sp.Float(x), self.m: 1})))
                    got = float(star(x + 1e-5) - star(x - 1e-5)) / 2e-5
                    self.assertAlmostEqual(abs(got) / want, 1.0, places=6)
                self.assertAlmostEqual(float(star(-1.0)), 0.0, places=9)
                near = float(star(1 + 1e-8) - star(1 + 1e-9)) / math.log(10)
                self.assertAlmostEqual(near, 2 * math.exp(-2 * q), places=6)
        self.assertAlmostEqual(float(nr._distorted_star("axis", 1 / 12)(500.0)), 0.3443, places=4)
        self.assertTrue(math.isinf(float(nr._distorted_star("equator", 1 / 12)(60.0))))


if __name__ == "__main__":
    unittest.main()
