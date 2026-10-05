"""The double Kerr solution, Kramer and Neugebauer's two Kerr black holes on one axis:
.venv.noindex/bin/python -m unittest discover -s _tools

What its texts and drawings state, held on the published files and on the pair every drawing
draws, Herdeiro and Rebelo's two holes of equal mass M and opposite angular momenta +-J at
J = M^2 and zeta = 4M, in units of M. The first class needs nothing but the files. The second
needs numpy and sympy and is skipped where they are absent: it holds the pair, in complex
arithmetic that no drawing uses, to Herdeiro and Rebelo's closed forms for the masses, the angular
momenta, the angular velocity and the surface gravity of the horizons and the force on the strut,
and Kramer and Neugebauer's potential to Ernst's equation.
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

SIGMA = math.sqrt(2 / 3)
POLES = (-2 - SIGMA, -2 + SIGMA, 2 - SIGMA, 2 + SIGMA)


def load(folder):
    return json.loads((DATA / folder / "double_kerr.json").read_text(encoding="utf-8"))


class Drawings(unittest.TestCase):
    def test_the_axis_marks_the_four_poles_and_draws_nothing_on_the_horizons(self):
        """f = 0 at |z| = 2 -+ sqrt(2/3), each marked; between two poles of one hole the axis is the
        horizon, hatched, and no ray or cone stands there."""
        view = next(v for v in load("diagrams")["systems"]["weyl"] if v["id"] == "axis")
        X0, X1 = view["box"][:2]
        marked = sorted(X0 + line[0][0] * (X1 - X0) for m in view["markers"] if m["kind"] == "grr" for line in m["lines"])
        self.assertEqual(len(marked), 4)
        for got, want in zip(marked, POLES):
            self.assertAlmostEqual(got, want, delta=2e-3)

        def on_a_rod(u):
            z = X0 + u * (X1 - X0)
            return POLES[0] + 0.02 < z < POLES[1] - 0.02 or POLES[2] + 0.02 < z < POLES[3] - 0.02
        self.assertTrue(view["hatch"])
        for family in view["rays"].values():
            for ray in family:
                self.assertFalse(any(on_a_rod(point[0]) for point in ray))
        self.assertFalse(any(on_a_rod(cone["at"][0]) for cone in view["cones"]))
        self.assertGreater(len(view["cones"]), 20)

    def test_the_tip_of_the_embedded_plane_is_a_cone_of_four_thirds_of_a_turn(self):
        """At the strut e^gamma = 3/4, so a circle the metric distance l from the axis has circumference
        (4/3) 2 pi l. The tip is drawn in Minkowski space, where the distance along the profile is
        sqrt(d rho^2 - d z^2)."""
        view = next(v for v in load("embedding")["views"] if v["id"] == "midplane")
        tip, plane = view["surfaces"][0]["pieces"]
        self.assertEqual(tip.get("space"), "minkowski")
        self.assertNotIn("space", plane)
        self.assertEqual(tip["points"][0], [0.0, 0.0, 0.0])
        _, rho, z = tip["points"][1]
        self.assertAlmostEqual(rho / math.sqrt(rho ** 2 - z ** 2), 4 / 3, delta=2e-3)
        # The two pieces meet in one circle.
        for a, b in zip(tip["points"][-1], plane["points"][0]):
            self.assertAlmostEqual(a, b, places=6)

    def test_the_circles_of_the_embedded_plane_are_rho_over_root_f(self):
        """On the axis midway between the holes f = 5/41, so the circles near the tip have radius
        rho sqrt(41/5); far out f -> 1 - 4m/r and the radius is rho/sqrt(f) > rho."""
        view = next(v for v in load("embedding")["views"] if v["id"] == "midplane")
        tip, plane = view["surfaces"][0]["pieces"]
        x, rho, _ = tip["points"][1]
        self.assertAlmostEqual(rho / x, math.sqrt(41 / 5), delta=2e-3)
        for x, rho, _ in plane["points"]:
            self.assertGreater(rho, x)

    def test_the_strut_is_a_whole_diamond_between_two_horizons(self):
        """On the conformal diagram of the axis between the holes all four edges are horizons and the
        moment embedded is the one event t = 0, z = 0, the centre of the diamond."""
        view = next(v for v in load("conformal")["views"] if v["id"] == "weyl_axis_between")
        edges = [layer for layer in view["layers"] if layer["kind"] == "line" and layer["class"] == "horizon"]
        self.assertEqual(len(edges), 4)
        self.assertFalse([layer for layer in view["layers"] if layer.get("class") == "scri"])
        (mark,) = view["slices"]
        self.assertEqual(mark["lines"], [])
        (at,) = mark["points"]
        self.assertAlmostEqual(at[0], 0.0, places=9)
        self.assertAlmostEqual(at[1], 0.0, places=9)


@unittest.skipUnless(HAS_NUMBERS, "numpy, scipy and sympy are needed")
class Pair(unittest.TestCase):
    """Herdeiro and Rebelo's (2.5) to (2.7) and (3.5) to (3.7) at M = J = 1 and zeta = 4, G = c = 1."""

    @classmethod
    def setUpClass(cls):
        import null_rays
        import numpy
        cls.nr, cls.np = null_rays, numpy
        cls.a = 2 * SIGMA                       # the length of each rod, their (3.6)

    def potential(self, z):
        """The imaginary part of the Ernst potential at a point of the axis, from f and gamma's
        companions: the complex arithmetic's own D_-/D_+."""
        ends, alpha = self.nr._double_kerr_floats()
        r = [abs(z - e) for e in ends]
        n = [alpha[j] * r[j] for j in range(4)]
        X = lambda i, k: (n[i] + n[k]) / (ends[i] - ends[k])  # noqa: E731
        det = lambda s: (X(0, 1) + s) * (X(2, 3) + s) - (X(0, 3) + s) * (X(2, 1) + s)  # noqa: E731
        return (det(-1) / det(1)).imag

    def test_the_rods_are_where_the_input_says(self):
        ends, _ = self.nr._double_kerr_floats()
        for got, want in zip(ends, POLES):
            self.assertAlmostEqual(got, want, places=12)

    def test_the_axis_is_an_axis_and_space_is_flat_far_away(self):
        """omega = 0 on all three stretches of the axis, gamma = 0 on the outer two, and
        f -> 1 - 2(2M)/r far away with omega f -> 0: no NUT charge and a total mass 2M."""
        np, numbers = self.np, self.nr._double_kerr_numbers
        z = np.array([-9.0, -4.0, -1.0, 0.0, 0.8, 3.5, 12.0])
        f, omega, e2g = numbers(0.0, z)
        self.assertLess(float(np.max(np.abs(omega))), 1e-12)
        self.assertTrue(np.all(f > 0))
        for i in (0, 1, 5, 6):
            self.assertAlmostEqual(float(e2g[i]), 1.0, places=12)
        f, omega, e2g = numbers(3.0e4, 4.0e4)
        self.assertAlmostEqual(float((1 - f) * 5.0e4), 4.0, places=3)
        self.assertAlmostEqual(float(e2g), 1.0, places=6)
        self.assertLess(abs(float(omega)), 1e-6)

    def test_the_strut_carries_the_force_of_two_schwarzschild_holes(self):
        """F = (e^(-gamma) - 1)/4 on the axis between the holes, their (2.5), is M^2/(zeta^2 - 4M^2) = 1/12."""
        e2g = self.nr._double_kerr_numbers(0.0, self.np.array([-0.9, 0.0, 0.5]))[2]
        for value in e2g:
            self.assertAlmostEqual((1 / math.sqrt(float(value)) - 1) / 4, 1 / 12, places=12)

    def test_each_hole_has_the_mass_the_spin_and_the_angular_velocity_stated(self):
        """On a rod omega = 1/Omega_H is constant; their (3.5) has Omega_H = (2M - a)/4J, opposite on
        the two holes. The Komar mass is -omega (Psi(top) - Psi(bottom))/4 and the angular momentum
        omega (M - a/2)/2, their (2.6) and (2.7): M = 1 and J = -+1."""
        np, numbers = self.np, self.nr._double_kerr_numbers
        want = 4 / (2 - self.a)
        for lo, hi, sign in ((POLES[0], POLES[1], 1), (POLES[2], POLES[3], -1)):
            z = np.linspace(lo + 0.05, hi - 0.05, 9)
            f, omega, _ = numbers(1e-9, z)
            self.assertTrue(np.all(f < 0))
            self.assertLess(float(np.max(np.abs(omega - sign * want))), 1e-6)
            mass = -sign * want * (self.potential(hi) - self.potential(lo)) / 4
            self.assertAlmostEqual(mass, 1.0, places=9)
            self.assertAlmostEqual(sign * want * (mass - self.a / 2) / 2, sign * 1.0, places=9)

    def test_the_poles_have_the_surface_gravity_of_the_temperature_stated(self):
        """Their (3.5): T = zeta M a (2M - a)/(16 pi J^2 (zeta - 2M)), so kappa = 2 pi T = a (2 - a)/4, and
        the tortoise coordinate of the axis runs off at a pole as ln|z - pole|/(2 kappa)."""
        kappa = self.a * (2 - self.a) / 4
        star = self.nr._double_kerr_star("axis")
        for pole, side in ((POLES[3], 1), (POLES[2], -1), (POLES[1], 1), (POLES[0], -1)):
            near, nearer = float(star(pole + side * 1e-5)), float(star(pole + side * 1e-7))
            self.assertAlmostEqual(side * (near - nearer) / math.log(100.0), 1 / (2 * kappa), places=3)

    def test_the_strings_the_drawings_use_are_the_pair(self):
        """The real expressions the spacetime diagrams are drawn from, and the closed form on the plane
        z = 0 the embedding is drawn from, against the complex arithmetic."""
        import sympy as sp
        np = self.np
        rho, z = sp.symbols("rho z", real=True)
        pair = {k: sp.lambdify((rho, z), sp.sympify(v, locals={"rho": rho, "z": z}), "numpy")
                for k, v in self.nr._double_kerr_pair().items()}
        plane = {k: sp.lambdify(rho, sp.sympify(v, locals={"rho": rho}), "numpy")
                 for k, v in self.nr._double_kerr_midplane().items()}
        for r, h in ((0.7, 0.0), (2.0, 3.0), (5.0, -1.5), (1.0, 6.0)):
            f, omega, e2g = (float(x) for x in self.nr._double_kerr_numbers(r, h))
            self.assertAlmostEqual(float(pair["f"](r, h)), f, places=11)
            self.assertAlmostEqual(float(pair["omega"](r, h)), omega, places=9)
            self.assertAlmostEqual(float(np.exp(2 * pair["gamma"](r, h))), e2g, places=11)
        for r in (0.0, 0.4, 1.0893994106890814, 3.0, 6.0):
            f, omega, e2g = (float(x) for x in self.nr._double_kerr_numbers(r, 0.0))
            self.assertAlmostEqual(float(plane["f"](r)), f, places=11)
            self.assertAlmostEqual(float(np.exp(2 * plane["gamma"](r))), e2g, places=11)
            self.assertLess(abs(omega), 1e-12)
        self.assertAlmostEqual(float(plane["f"](0.0)), 5 / 41, places=12)

    def test_the_potential_solves_ernsts_equation(self):
        """Kramer and Neugebauer's potential at the pair's rod ends and phases, in forty digits at four
        points off the axis: f Laplacian(E) = (grad E)^2."""
        import print_charts
        import sympy as sp
        x, y = sp.symbols("x y", positive=True)
        ends, phases = self.nr._double_kerr_data()
        minus, plus, _ = print_charts.double_kerr_potential(x, y, ends, [co + sp.I * si for co, si in phases])
        E = minus / plus
        lap = sp.diff(E, x, 2) + sp.diff(E, x) / x + sp.diff(E, y, 2)
        ernst = (E + sp.conjugate(E)) / 2 * lap - sp.diff(E, x) ** 2 - sp.diff(E, y) ** 2
        for at in ({x: sp.Rational(7, 10), y: sp.Rational(1, 3)}, {x: 2, y: 3}, {x: sp.Rational(1, 4), y: -2},
                   {x: 5, y: sp.Rational(-3, 2)}):
            self.assertLess(abs(sp.N(ernst.subs(at), 40)), sp.Float("1e-30") * (1 + abs(sp.N(lap.subs(at), 40))))


if __name__ == "__main__":
    unittest.main()
