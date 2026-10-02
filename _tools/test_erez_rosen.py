"""Erez and Rosen's quadrupole, Schwarzschild's mass with a quadrupole moment in Weyl's class:
python3 -m unittest discover -s _tools

What its texts and drawings state, held on the published files and on the two deformations every
drawing draws, the prolate q = 1 and the oblate q = -1/2 at m = 1. The first class needs nothing
but the files. The second needs numpy, scipy and sympy and is skipped where they are absent: it
holds the published psi and gamma to Weyl's field equations, to Schwarzschild's functions at
q = 0, to the mass m and the quadrupole moment (2/15) q m^3, the coefficient the translation of
Doroshkevich, Zel'dovich and Novikov's paper prints to failing those equations, the surface x = 1
to the exponents the captions state, and the floats the drawings are checked against to the
published functions.
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


def load(folder):
    return json.loads((DATA / folder / "erez_rosen.json").read_text(encoding="utf-8"))


class Drawings(unittest.TestCase):
    def test_each_chart_draws_both_planes_of_both_deformations(self):
        systems = load("diagrams")["systems"]
        want = {f"{plane}_{shape}" for plane in ("axis", "equator") for shape in ("prolate", "oblate")}
        for system in ("prolate_spheroidal", "spherical"):
            self.assertEqual({v["id"] for v in systems[system]}, want)
        conformal = {v["id"] for v in load("conformal")["views"]}
        self.assertEqual(conformal, {f"{s}_{w}" for s in ("prolate_spheroidal", "spherical") for w in want})

    def test_only_the_prolate_axis_has_no_singularity(self):
        """At q = 1 the Kretschmann scalar is finite on the axis, 3 e^6/(4 m^4) at x = 1, so that
        view alone is the whole diamond with the surface an edge of the chart."""
        for view in load("conformal")["views"]:
            singular = any(layer.get("class") == "singular" for layer in view["layers"])
            self.assertEqual(singular, not view["id"].endswith("axis_prolate"), view["id"])

    def test_the_prolate_equator_shrinks_and_the_oblate_one_has_a_neck(self):
        """The circle about the axis has radius r f^(q (3x^2 - 1)/8) e^(3qx/4): at q = 1 it grows
        with r all the way out, and at q = -1/2 it is least, 2.2867 m, at r = 2.1297 m."""
        views = {v["id"]: v for v in load("embedding")["views"]}
        prolate = [p[1] for p in views["prolate"]["surfaces"][0]["pieces"][0]["points"]]
        self.assertEqual(prolate, sorted(prolate))
        self.assertAlmostEqual(prolate[0], 1.371, places=2)
        oblate = views["oblate"]["surfaces"][0]["pieces"][0]["points"]
        radius, at = min((p[1], p[0]) for p in oblate)
        self.assertAlmostEqual(radius, 2.2867, places=3)
        self.assertAlmostEqual(at, 2.1297, places=2)
        self.assertGreater(oblate[0][1], radius)

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
        entry = next(c for c in load("metrics")["coordinates"] if c["id"] == "prolate_spheroidal")
        cls.reader = vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], (),
                               held=vm.HELD[("erez_rosen", "prolate_spheroidal")])
        r = cls.reader
        cls.x, cls.y = r.symbol["x"], r.symbol["y"]
        cls.m, cls.q = r.parameters["m"], r.parameters["q"]
        cls.psi, cls.gamma = r.held[r.parameters["psi"]], r.held[r.parameters["gamma"]]

    def quadrature(self, P):
        sp, x, y = self.sp, self.x, self.y
        px, py = sp.diff(P, x), sp.diff(P, y)
        along = (1 - y ** 2) / (x ** 2 - y ** 2) * (x * (x ** 2 - 1) * px ** 2 - x * (1 - y ** 2) * py ** 2
                                                      - 2 * y * (x ** 2 - 1) * px * py)
        across = (x ** 2 - 1) / (x ** 2 - y ** 2) * (y * (x ** 2 - 1) * px ** 2 - y * (1 - y ** 2) * py ** 2
                                                       + 2 * x * (1 - y ** 2) * px * py)
        return along, across

    def test_the_published_functions_solve_weyls_equations(self):
        sp, x, y, P, G = self.sp, self.x, self.y, self.psi, self.gamma
        laplace = sp.diff((x ** 2 - 1) * sp.diff(P, x), x) + sp.diff((1 - y ** 2) * sp.diff(P, y), y)
        self.assertEqual(self.vm.norm(laplace), 0)
        along, across = self.quadrature(P)
        self.assertEqual(self.vm.norm(sp.diff(G, x) - along), 0)
        self.assertEqual(self.vm.norm(sp.diff(G, y) - across), 0)
        self.assertEqual(self.vm.norm(G.subs(y, 1)), 0)

    def test_the_coefficient_the_translation_prints_fails_them(self):
        """The English translation of Doroshkevich, Zel'dovich and Novikov's Appendix I prints
        (1 + q + q^2)/2 before the first logarithm of gamma, where the published gamma has
        (1 + q)^2/2. The two differ by (q/2) ln((x^2 - 1)/(x^2 - y^2)), whose derivatives do not
        vanish, so with the published psi the printed gamma misses both quadratures."""
        sp, x, y, q = self.sp, self.x, self.y, self.q
        printed = self.gamma - q / 2 * sp.log((x ** 2 - 1) / (x ** 2 - y ** 2))
        self.assertEqual(sp.expand(printed.coeff(sp.log((x ** 2 - 1) / (x ** 2 - y ** 2))) - (1 + q + q ** 2) / 2), 0)
        along, across = self.quadrature(self.psi)
        point = {x: sp.Rational(17, 7), y: sp.Rational(3, 11), q: sp.Rational(5, 13)}
        for found, wanted in ((sp.diff(printed, x), along), (sp.diff(printed, y), across)):
            self.assertGreater(abs(sp.N((found - wanted).subs(point), 30)), 1e-3)

    def test_without_the_quadrupole_it_is_schwarzschild(self):
        sp, x, y, q = self.sp, self.x, self.y, self.q
        self.assertEqual(self.vm.norm(self.psi.subs(q, 0) - sp.log((x - 1) / (x + 1)) / 2), 0)
        self.assertEqual(self.vm.norm(self.gamma.subs(q, 0) - sp.log((x ** 2 - 1) / (x ** 2 - y ** 2)) / 2), 0)

    def test_the_mass_is_m_and_the_quadrupole_moment_two_fifteenths_of_q_m_cubed(self):
        """Geroch's moments of a static field are the coefficients of -tanh(psi) on the axis in
        powers of 1/z, and on the axis Weyl's z = m x."""
        sp, x, y, q = self.sp, self.x, self.y, self.q
        u = sp.Symbol("u", positive=True)
        xi = sp.series(-sp.tanh(self.psi.subs(y, 1).subs(x, 1 / u)), u, 0, 5).removeO()
        self.assertEqual(sp.expand(xi - u - sp.Rational(2, 15) * q * u ** 3), 0)

    def test_the_surface_behaves_as_a_rod_whose_density_varies_along_it(self):
        """Toward x = 1, 2 psi -> delta ln(x - 1) and 2 gamma -> delta^2 ln(x - 1) with
        delta = 1 + q P_2(y): 1 + q on the axis, where gamma vanishes, and 1 - q/2 on the equator."""
        sp, x, y, q = self.sp, self.x, self.y, self.q
        d = sp.Symbol("d", positive=True)
        for at, value in ((0, sp.Rational(1)), (0, sp.Rational(-1, 2)), (sp.Rational(1, 2), sp.Rational(1))):
            delta = 1 + value * (3 * at ** 2 - 1) / 2
            for function, power in ((self.psi, delta), (self.gamma, delta ** 2)):
                f = function.subs({q: value, y: at})
                slope = sp.limit((f.subs(x, 1 + d) / sp.log(d)), d, 0, "+")
                self.assertEqual(sp.nsimplify(2 * slope), power)

    def test_the_drawings_floats_are_the_published_functions(self):
        """_erez_rosen_numbers, which the tortoise coordinates of every drawing are integrated
        from, against the published psi and gamma at m = 1."""
        sp = self.sp
        for q in (1.0, -0.5):
            for x0, y0 in ((1.3, 0.0), (2.5, 0.0), (4.0, 1.0), (1.05, 0.6)):
                at = {self.x: x0, self.y: sp.Float(y0) if y0 != 1.0 else 1, self.q: sp.Rational(q), self.m: 1}
                psi, gamma = self.nr._erez_rosen_numbers(math.log(x0 - 1), y0, q)
                self.assertAlmostEqual(psi, float(self.psi.subs(at)), places=10)
                self.assertAlmostEqual(gamma, float(self.gamma.subs(at)), places=10)

    def test_the_tortoise_coordinates_have_the_slopes_the_captions_state(self):
        """dx_*/dx is e^{-2 psi} on the axis and x e^{gamma - 2 psi}/sqrt(x^2 - 1) in the plane, the
        prolate axis's x_* falls without bound toward x = 1, the other three vanish there, and in
        the oblate plane light from x = 1.1 takes 69.9 m/c to arrive."""
        nr = self.nr
        for q in (1.0, -0.5):
            for plane, y in (("axis", 1.0), ("equator", 0.0)):
                star = nr._erez_rosen_star(plane, q)
                for x in (1.3, 2.5, 4.0):
                    h = 1e-5
                    psi, gamma = nr._erez_rosen_numbers(math.log(x - 1), y, q)
                    slope = math.exp(-2 * psi) if plane == "axis" else math.exp(gamma - 2 * psi) * x / math.sqrt(x * x - 1)
                    self.assertAlmostEqual(float(star(x + h) - star(x - h)) / (2 * h) / slope, 1.0, places=6)
                if (plane, q) == ("axis", 1.0):
                    self.assertLess(float(star(1 + 1e-6)), -1e5)
                else:
                    self.assertEqual(float(star(1.0)), 0.0)
        self.assertAlmostEqual(float(nr._erez_rosen_star("equator", -0.5)(1.1)), 69.94, places=1)


if __name__ == "__main__":
    unittest.main()
