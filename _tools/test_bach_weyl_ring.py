"""Bach and Weyl's ring, the static vacuum field of a thin circular ring of matter in Weyl's class:
python3 -m unittest discover -s _tools

What its texts and drawings state, held on the published files and on the ring every drawing
draws, m = a/2 in units of the ring's radius a. The first class needs nothing but the files. The
second needs numpy, scipy and sympy and is skipped where they are absent: it holds the published
psi and gamma to Weyl's field equations by the derivatives of the two elliptic integrals, to
Curzon and Chazy's functions at a = 0, to the mass m, the three charts to one value at one event,
the ring to the behaviour the captions state on its two sides, and the floats the drawings are
checked against to the published functions.
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

CHARTS = ("weyl", "toroidal", "oblate_spheroidal")


def load(folder):
    return json.loads((DATA / folder / "bach_weyl_ring.json").read_text(encoding="utf-8"))


class Files(unittest.TestCase):
    def test_the_three_charts_are_published_in_order(self):
        self.assertEqual(tuple(c["id"] for c in load("metrics")["coordinates"]), CHARTS)

    def test_every_chart_is_a_vacuum(self):
        for chart in load("metrics")["coordinates"]:
            self.assertEqual(chart["ricci_scalar"], "R = 0", chart["id"])
            for variant in chart["ricci_tensor"]["variants"].values():
                for entry in variant["nonzero"]:
                    self.assertIn(entry["value"], ("0", "-0"), chart["id"])

    def test_each_chart_draws_the_axis_and_the_plane_on_both_sides_of_the_ring(self):
        systems = load("diagrams")["systems"]
        want = {"weyl": {"axis", "outside", "inside"}, "toroidal": {"axis", "outer", "inner"},
                "oblate_spheroidal": {"axis", "plane", "disc"}}
        self.assertEqual({s: {v["id"] for v in views} for s, views in systems.items()}, want)

    def test_the_spacetime_diagrams_are_near_square(self):
        for views in load("diagrams")["systems"].values():
            for view in views:
                X0, X1, Y0, Y1 = view["box"]
                self.assertGreaterEqual((X1 - X0) / (Y1 - Y0), 0.5)
                self.assertLessEqual((X1 - X0) / (Y1 - Y0), 2.0)

    def test_the_ring_is_singular_from_outside_and_at_infinity_from_inside(self):
        """Outside the ring light arrives in a finite time, so the ring is a timelike singular line
        of the triangle; inside it light never arrives, so the ring is the triangle's far edges."""
        views = {v["id"]: v for v in load("conformal")["views"]}
        outside = {"weyl_outside", "toroidal_outer", "oblate_plane"}
        self.assertEqual(len(views), 9)
        for name, view in views.items():
            singular = any(layer.get("class") == "singular" for layer in view["layers"])
            self.assertEqual(singular, name in outside, name)

    def test_the_embedded_plane_has_a_neck_outside_the_ring_and_none_inside(self):
        """The circle about the axis has radius rho e^(-psi): outside the ring it is least, 2.075 a,
        at rho = 1.220 a, and on the disc inside it grows all the way from the axis."""
        views = {v["id"]: v for v in load("embedding")["views"]}
        outside = views["outside"]["surfaces"][0]["pieces"][0]["points"]
        radius, at = min((p[1], p[0]) for p in outside)
        self.assertAlmostEqual(radius, 2.0751, places=3)
        self.assertAlmostEqual(at, 1.220, places=2)
        self.assertGreater(outside[0][1], radius)
        inside = [p[1] for p in views["inside"]["surfaces"][0]["pieces"][0]["points"]]
        self.assertEqual(inside, sorted(inside))


@unittest.skipUnless(HAS_NUMBERS, "needs numpy, scipy and sympy")
class Physics(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import numpy as np
        import sympy as sp

        import null_rays as nr
        import verify_metrics as vm
        cls.np, cls.sp, cls.nr, cls.vm = np, sp, nr, vm
        cls.readers = {}
        for entry in load("metrics")["coordinates"]:
            cls.readers[entry["id"]] = vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], (),
                                                 held=vm.HELD[("bach_weyl_ring", entry["id"])])
        r = cls.readers["weyl"]
        cls.rho, cls.z = r.symbol["\\rho"], r.symbol["z"]
        cls.m, cls.a = r.parameters["m"], r.parameters["a"]
        cls.psi, cls.gamma = (r.held[r.parameters[n]].subs(r.held) for n in ("psi", "gamma"))

    def test_the_elliptic_integrals_are_read_with_their_parameter(self):
        """K(1/2) = Gamma(1/4)^2/(4 sqrt(pi)), K(0) = E(0) = pi/2, and each derivative is written
        in the two integrals again."""
        sp, vm = self.sp, self.vm
        r = self.readers["weyl"]
        self.assertEqual(r("\\mathrm{K}\\left(\\kappa\\right)").func, vm.EllipticK)
        self.assertEqual(r("\\mathrm{E}\\left(\\kappa\\right)").func, vm.EllipticE)
        self.assertEqual(vm.EllipticK(0), sp.pi / 2)
        self.assertEqual(vm.EllipticE(0), sp.pi / 2)
        half = sp.Rational(1, 2)
        self.assertAlmostEqual(float(vm.EllipticK(half).evalf(30)), math.gamma(0.25) ** 2 / (4 * math.sqrt(math.pi)),
                               places=12)
        x = sp.Symbol("x", positive=True)
        for function in (vm.EllipticK, vm.EllipticE):
            slope = sp.diff(function(x), x).subs(x, half)
            step = sp.Rational(1, 10 ** 8)
            numeric = (function(half + step).evalf(40) - function(half - step).evalf(40)) / (2 * step)
            self.assertAlmostEqual(float(slope.evalf(30)), float(numeric), places=10)
        f = sp.lambdify(x, vm.EllipticK(x) + vm.EllipticE(x), "numpy")
        self.assertAlmostEqual(float(f(0.3)), float((vm.EllipticK(0.3) + vm.EllipticE(0.3)).evalf(20)), places=12)

    def test_the_published_functions_solve_weyls_equations(self):
        """Laplace's equation for psi and both quadratures for gamma, as identities in rho, z and
        the two elliptic integrals."""
        sp, rho, z, P, G = self.sp, self.rho, self.z, self.psi, self.gamma
        laplace = sp.diff(P, rho, 2) + sp.diff(P, rho) / rho + sp.diff(P, z, 2)
        self.assertEqual(self.vm.norm(laplace), 0)
        self.assertEqual(self.vm.norm(sp.diff(G, rho) - rho * (sp.diff(P, rho) ** 2 - sp.diff(P, z) ** 2)), 0)
        self.assertEqual(self.vm.norm(sp.diff(G, z) - 2 * rho * sp.diff(P, rho) * sp.diff(P, z)), 0)

    def test_a_misprinted_gamma_fails_them(self):
        """With the sign between its two terms turned, gamma misses the quadrature along rho."""
        sp, rho, z, P = self.sp, self.rho, self.z, self.psi
        vm, m, a = self.vm, self.m, self.a
        kappa = 4 * a * rho / ((rho + a) ** 2 + z ** 2)
        K, E = vm.EllipticK(kappa), vm.EllipticE(kappa)
        wrong = -m ** 2 / (4 * sp.pi ** 2 * a ** 2 * rho) * ((rho + a) * (E - K) ** 2
                                                           - (rho - a) * (E - (1 - kappa) * K) ** 2 / (1 - kappa))
        point = {rho: sp.Rational(17, 10), z: sp.Rational(3, 11), m: sp.Rational(1, 2), a: 1}
        miss = sp.diff(wrong, rho) - rho * (sp.diff(P, rho) ** 2 - sp.diff(P, z) ** 2)
        self.assertGreater(abs(sp.N(miss.subs(point), 30)), 1e-3)

    def test_a_ring_of_no_radius_is_the_curzon_chazy_particle(self):
        sp, rho, z, m, a = self.sp, self.rho, self.z, self.m, self.a
        R = sp.sqrt(rho ** 2 + z ** 2)
        self.assertEqual(sp.simplify(self.psi.subs(a, 0) + m / R), 0)
        point = {rho: sp.Rational(3, 2), z: sp.Rational(4, 5), m: sp.Rational(7, 5)}
        curzon = (-m ** 2 * rho ** 2 / (2 * R ** 4)).subs(point)
        self.assertAlmostEqual(float(sp.N(self.gamma.subs(point).subs(a, sp.Rational(1, 10 ** 10)), 40)), float(curzon),
                               places=12)

    def test_the_mass_is_m_and_the_axis_is_regular(self):
        sp, rho, z, m, a = self.sp, self.rho, self.z, self.m, self.a
        self.assertEqual(sp.simplify(self.psi.subs(rho, 0) + m / sp.sqrt(z ** 2 + a ** 2)), 0)
        point = {m: sp.Rational(1, 2), a: 1, z: sp.Rational(4, 5)}
        self.assertLess(abs(sp.N(self.gamma.subs(point).subs(rho, sp.Rational(1, 10 ** 12)), 40)), 1e-20)
        far = sp.Integer(10) ** 9
        self.assertAlmostEqual(float(sp.N((self.psi * far).subs(point).subs(rho, far), 30)), -0.5, places=8)

    def test_one_event_has_one_psi_and_one_gamma_in_all_three_charts(self):
        sp = self.sp
        m0, a0 = sp.Rational(1, 2), sp.Integer(1)
        for zeta, sigma in ((sp.Rational(3, 5), sp.Rational(7, 10)), (sp.Rational(2), sp.Rational(5, 2)),
                            (sp.Rational(1), sp.Rational(4))):
            d = sp.cosh(zeta) - sp.cos(sigma)
            rho0, z0 = sp.N(a0 * sp.sinh(zeta) / d, 40), sp.N(a0 * sp.sin(sigma) / d, 40)
            # rho = a sqrt((1 + xi^2)(1 - eta^2)) and z = a xi eta, solved for xi and eta.
            s = rho0 ** 2 + z0 ** 2 - 1
            xi2 = (s + sp.sqrt(s ** 2 + 4 * z0 ** 2)) / 2
            xi, eta = sp.sqrt(xi2), z0 / sp.sqrt(xi2)
            values = []
            for chart, at in (("weyl", {"\\rho": rho0, "z": z0}), ("toroidal", {"\\zeta": zeta, "\\sigma": sigma}),
                              ("oblate_spheroidal", {"\\xi": xi, "\\eta": eta})):
                r = self.readers[chart]
                point = {r.symbol[name]: value for name, value in at.items()}
                point.update({r.parameters["m"]: m0, r.parameters["a"]: a0})
                values.append([sp.N(r.held[r.parameters[n]].subs(r.held).subs(r.held).subs(point), 30)
                               for n in ("psi", "gamma")])
            for other in values[1:]:
                for one, two in zip(values[0], other):
                    self.assertAlmostEqual(float(one), float(two), places=20)

    def test_the_ring_is_near_from_outside_and_far_from_inside(self):
        """Toward the ring gamma -> -+(m^2/2 pi^2 a^2) e^zeta in its plane: it falls without bound
        outside the ring, sigma = 0, and rises without bound inside it, sigma = pi, while e^psi
        vanishes on both sides as (l_1/8a)^(m/pi a)."""
        nr = self.nr
        for zeta in (8.0, 10.0):
            out = nr._bach_weyl_numbers(1 / math.tanh(zeta / 2), 0.0)
            inn = nr._bach_weyl_numbers(math.tanh(zeta / 2), 0.0)
            lead = 0.25 / (2 * math.pi ** 2) * math.exp(zeta)
            self.assertAlmostEqual(out[1] / -lead, 1.0, delta=0.03)
            self.assertAlmostEqual(inn[1] / lead, 1.0, delta=0.03)
            for psi, rho in ((out[0], 1 / math.tanh(zeta / 2)), (inn[0], math.tanh(zeta / 2))):
                self.assertAlmostEqual(psi / (0.5 / math.pi * math.log(abs(rho - 1) / 8)), 1.0, delta=0.02)

    def test_the_drawings_floats_are_the_published_functions(self):
        sp = self.sp
        for rho0, z0 in ((0.3, 0.0), (0.9, 0.0), (1.5, 0.0), (2.0, 1.2), (0.6, -0.8)):
            at = {self.rho: sp.Float(rho0, 30), self.z: sp.Float(z0, 30), self.m: sp.Rational(1, 2), self.a: 1}
            psi, gamma = self.nr._bach_weyl_numbers(rho0, z0)
            self.assertAlmostEqual(psi, float(sp.N(self.psi.subs(at), 30)), places=12)
            self.assertAlmostEqual(gamma, float(sp.N(self.gamma.subs(at), 30)), places=12)
        self.assertEqual(self.nr._bach_weyl_numbers(0.0, 0.7)[1], 0.0)

    def test_the_tortoise_coordinates_have_the_slopes_and_the_times_the_captions_state(self):
        """dz_*/dz = e^(-2 psi) on the axis and d rho_*/d rho = e^(gamma - 2 psi) in the plane; light
        from the axis is at rho = 0.99 a at ct = 3.67 a and never at the ring, and light from
        rho = 2a outside reaches the ring at ct = 2.00 a."""
        nr = self.nr
        axis, plane = nr._bach_weyl_star("weyl_axis"), nr._bach_weyl_star("weyl_plane")
        h = 1e-5
        for z in (-2.0, 0.0, 1.5):
            slope = math.exp(-2 * nr._bach_weyl_numbers(0.0, z)[0])
            self.assertAlmostEqual(float(axis(z + h) - axis(z - h)) / (2 * h) / slope, 1.0, places=6)
        self.assertAlmostEqual(math.exp(-2 * nr._bach_weyl_numbers(0.0, 0.0)[0]), math.e, places=12)
        for rho in (0.4, 0.9, 1.3, 3.0):
            psi, gamma = nr._bach_weyl_numbers(rho, 0.0)
            self.assertAlmostEqual(float(plane(rho + h) - plane(rho - h)) / (2 * h) / math.exp(gamma - 2 * psi), 1.0,
                                   places=5)
        self.assertAlmostEqual(float(plane(0.99)), 3.67, places=2)
        self.assertGreater(float(plane(0.999)), 1e7)
        self.assertAlmostEqual(float(plane(2.0)), 2.00, places=2)
        # The same events through the other two charts.
        self.assertAlmostEqual(float(nr._bach_weyl_star("oblate_plane")(math.sqrt(3.0))), float(plane(2.0)), places=10)
        self.assertAlmostEqual(float(nr._bach_weyl_star("toroidal_outer")(2 * math.atanh(0.5))), float(plane(2.0)),
                               places=10)
        self.assertAlmostEqual(float(nr._bach_weyl_star("toroidal_inner")(2 * math.atanh(0.9))), float(plane(0.9)),
                               places=10)
        self.assertAlmostEqual(float(nr._bach_weyl_star("oblate_disc")(0.6)), float(plane(0.8)), places=10)

    def test_the_neck_and_the_photon_orbit_lie_outside_the_ring(self):
        """rho e^(-psi) is least at rho = 1.220 a, where rho d_rho psi = 1, and the circular photon
        orbit, 2 rho d_rho psi = 1, is at rho = 1.535 a."""
        from scipy.optimize import brentq
        nr = self.nr

        def slope(rho, h=1e-6):
            return rho * (nr._bach_weyl_numbers(rho + h, 0.0)[0] - nr._bach_weyl_numbers(rho - h, 0.0)[0]) / (2 * h)
        neck = brentq(lambda r: slope(r) - 1, 1.05, 2.0)
        self.assertAlmostEqual(neck, 1.2201, places=3)
        self.assertAlmostEqual(neck * math.exp(-nr._bach_weyl_numbers(neck, 0.0)[0]), 2.0751, places=3)
        self.assertAlmostEqual(brentq(lambda r: 2 * slope(r) - 1, 1.05, 3.0), 1.5350, places=3)


if __name__ == "__main__":
    unittest.main()
