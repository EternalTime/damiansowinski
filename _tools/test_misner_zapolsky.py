"""The star of infinite central density, held on the published file with the checker's own reader:
.venv.noindex/bin/python -m unittest discover -s _tools

The sphere of radiation is a perfect fluid with p = rho c^2/3, 8 pi G rho/c^2 = 3/(7 r^2) and the
mass 3r/14 inside r, Oppenheimer and Volkoff's (22), and it is the same sphere at every scale.
Tolman's solutions V and VI at n = 1/2 have it for their centre, lose their pressure at r_b and
meet Schwarzschild's published exterior there, with the masses of his (7.9) and (8.9). His
solution V with n free has p/rho = n/(2 - n) and is the sphere of radiation at n = 1/2. A star
of p = rho c^2/3 with a finite central density settles onto the sphere of radiation as
r^(-3/4) cos((sqrt(47)/4) ln r), and its 2Gm/(c^2 r) peaks at 0.493, Chavanis's numbers. The
charts need sympy and are skipped where it is absent; the star of finite central density is
integrated in plain Python.
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


def published(metric_id):
    return json.loads((DATA / "metrics" / f"{metric_id}.json").read_text(encoding="utf-8"))


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class TheCharts(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import sympy
        import verify_metrics
        cls.sp, cls.vm = sympy, verify_metrics
        cls.charts = {}
        for entry in published("misner_zapolsky")["coordinates"]:
            reader = cls.vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
            g = cls.vm.metric_from_line_element(reader, entry["line_element"], entry["coords"])
            cls.charts[entry["id"]] = (reader, g.subs(reader.c, 1))

    def fluid(self, chart):
        """(8 pi G rho/c^2, 8 pi G p/c^4 along r, the same across r) of a static sphere
        -A dt^2 + B dr^2 + r^2 dOmega^2, from its Einstein tensor in A and B."""
        sp = self.sp
        reader, g = self.charts[chart]
        r = reader.symbol["r"]
        A, B = -g[0, 0], g[1, 1]
        nu, lam = sp.log(A), sp.log(B)
        density = (sp.diff(lam, r) / r - 1 / r ** 2) / B + 1 / r ** 2
        radial = (sp.diff(nu, r) / r + 1 / r ** 2) / B - 1 / r ** 2
        across = (sp.diff(nu, r, 2) / 2 - sp.diff(lam, r) * sp.diff(nu, r) / 4 + sp.diff(nu, r) ** 2 / 4
                  + (sp.diff(nu, r) - sp.diff(lam, r)) / (2 * r)) / B
        return reader, r, density, radial, across

    def zero(self, expression):
        sp = self.sp
        positive = {s: sp.Symbol(s.name + "_positive", positive=True) for s in expression.free_symbols}
        self.assertEqual(sp.simplify(sp.powsimp(sp.powdenest(expression.subs(positive), force=True), force=True)), 0)

    def test_the_sphere_of_radiation_is_oppenheimer_and_volkoffs_exact_solution(self):
        reader, r, density, radial, across = self.fluid("areal")
        self.zero(density - 3 / (7 * r ** 2))
        self.zero(radial - density / 3)
        self.zero(across - radial)
        # u = 3r/14 in their (22): g^rr = 1 - 2u/r = 4/7.
        g = self.charts["areal"][1]
        self.zero(1 / g[1, 1] - (1 - 2 * self.sp.Rational(3, 14)))

    def test_the_sphere_of_radiation_is_the_same_at_every_scale(self):
        # r -> k r with t -> sqrt(k) t multiplies the line element by k^2.
        sp = self.sp
        reader, g = self.charts["areal"]
        r, k = reader.symbol["r"], sp.Symbol("k", positive=True)
        self.zero(g[0, 0].subs(r, k * r) * k - k ** 2 * g[0, 0])
        self.zero(g[1, 1].subs(r, k * r) * k ** 2 - k ** 2 * g[1, 1])
        self.zero(g[2, 2].subs(r, k * r) - k ** 2 * g[2, 2])

    def test_its_published_scalars(self):
        entry = next(c for c in published("misner_zapolsky")["coordinates"] if c["id"] == "areal")
        self.assertEqual(entry["ricci_scalar"], "R = 0")
        self.assertEqual(entry["kretschmann"], "K = \\dfrac{72}{49r^4}")

    def test_each_star_is_a_perfect_fluid_that_loses_its_pressure_at_its_surface(self):
        for chart in ("tolman_v", "tolman_vi"):
            reader, r, density, radial, across = self.fluid(chart)
            rb = reader.parameters["r_b"]
            self.zero(across - radial)
            self.zero(radial.subs(r, rb))
            # The density at the surface is still positive: Tolman's (7.7) and (8.7).
            self.assertGreater(float(density.subs({r: 1, rb: 1})), 0)

    def test_each_star_has_the_sphere_of_radiation_at_its_centre(self):
        sp = self.sp
        for chart in ("tolman_v", "tolman_vi"):
            reader, r, density, radial, across = self.fluid(chart)
            rb = reader.parameters["r_b"]
            at = {rb: 1}
            self.assertEqual(sp.limit((radial / density).subs(at), r, 0), sp.Rational(1, 3))
            self.assertEqual(sp.limit((density * r ** 2).subs(at), r, 0), sp.Rational(3, 7))
            g = self.charts[chart][1]
            self.assertEqual(sp.limit(g[1, 1].subs(at), r, 0), sp.Rational(7, 4))
            # Inside, the pressure never passes a third of the energy density.
            for x in (sp.Rational(1, 10), sp.Rational(1, 2), sp.Rational(9, 10)):
                self.assertLess(float((radial / density).subs(at).subs(r, x)), 1 / 3)

    def test_each_star_meets_schwarzschilds_published_exterior(self):
        sp = self.sp
        entry = next(c for c in published("schwarzschild")["coordinates"] if c["id"] == "spherical")
        other = self.vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
        vacuum = {tuple(e["indices"]): other(e["value"]) for e in entry["metric_components"]}
        # Tolman's (7.9), m = r_b/4, and (8.9), m = 3r_b/14, with r_s = 2m.
        for chart, mass in (("tolman_v", sp.Rational(1, 4)), ("tolman_vi", sp.Rational(3, 14))):
            reader, g = self.charts[chart]
            r, rb = reader.symbol["r"], reader.parameters["r_b"]
            there = {other.symbol["r"]: rb, other.parameters["r_s"]: 2 * mass * rb}
            self.zero(g[0, 0].subs(r, rb) - vacuum["t", "t"].subs(there))
            self.zero(g[1, 1].subs(r, rb) - vacuum["r", "r"].subs(there))

    def test_tolmans_equation_of_state_for_solution_vi(self):
        # His (8.5) with B/A = 1/(9 r_b): p = (rho/3)(1 - 9k)/(1 - k), k = (3/(56 pi rho))^(1/2) B/A,
        # and with 8 pi rho = 3/(7 r^2) that k is r/(9 r_b).
        reader, r, density, radial, across = self.fluid("tolman_vi")
        rb = reader.parameters["r_b"]
        k = r / (9 * rb)
        self.zero(radial - density / 3 * (1 - 9 * k) / (1 - k))

    def test_tolmans_exponent(self):
        sp = self.sp
        reader, r, density, radial, across = self.fluid("power_law")
        n = reader.parameters["n"]
        self.zero(across - radial)
        self.zero(radial / density - n / (2 - n))
        # n = 1 is the stiffest fluid, p = rho c^2, and n = 1/2 is the sphere of radiation.
        self.assertEqual(sp.simplify((radial / density).subs(n, 1)), 1)
        g, core = self.charts["power_law"][1], self.charts["areal"][1]
        same = {self.charts["areal"][0].symbol["r"]: r,
                self.charts["areal"][0].parameters["a"]: reader.parameters["a"]}
        for i in range(2):
            self.zero(g[i, i].subs(n, sp.Rational(1, 2)) - core[i, i].subs(same))

    def test_light_reaches_the_centre_in_a_finite_time_unless_the_fluid_is_stiffest(self):
        sp = self.sp
        reader, g = self.charts["power_law"]
        r, n, a = reader.symbol["r"], reader.parameters["n"], reader.parameters["a"]
        speed = sp.sqrt(g[1, 1] / (-g[0, 0]))            # d(ct)/dr along a radial ray
        affine = sp.sqrt(-g[0, 0] * g[1, 1])             # d(lambda)/dr, with the energy set to 1
        half = {n: sp.Rational(1, 2), a: 1}
        self.assertEqual(sp.simplify(sp.integrate(speed.subs(half), (r, 0, 3)) - sp.sqrt(21)), 0)   # sqrt(7 a r)
        stiff = {n: 1, a: 1}
        self.assertEqual(sp.integrate(speed.subs(stiff), (r, 0, 1)), sp.oo)
        self.assertEqual(sp.simplify(sp.integrate(affine.subs(stiff), (r, 0, 1)) - 1 / sp.sqrt(2)), 0)


class TheApproach(unittest.TestCase):
    """A star of p = rho c^2/3 with a finite central density, in u = Gm/(c^2 r) and
    v = 4 pi G rho r^2/c^2 against s = ln r: du/ds = v - u and dv/ds = 2v - 4v(u + v/3)/(1 - 2u),
    which is the equation of hydrostatic equilibrium. The sphere of radiation is the fixed point
    u = v = 3/14."""

    @staticmethod
    def flow(u, v):
        return v - u, 2 * v - 4 * v * (u + v / 3) / (1 - 2 * u)

    @classmethod
    def setUpClass(cls):
        step, s = 2e-4, math.log(1e-4)
        r = math.exp(s)
        u, v = 4 * math.pi / 3 * r * r, 4 * math.pi * r * r      # central density 1, G = c = 1
        cls.path = []
        while s < 9.5:
            k1 = cls.flow(u, v)
            k2 = cls.flow(u + step / 2 * k1[0], v + step / 2 * k1[1])
            k3 = cls.flow(u + step / 2 * k2[0], v + step / 2 * k2[1])
            k4 = cls.flow(u + step * k3[0], v + step * k3[1])
            u += step / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
            v += step / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
            s += step
            cls.path.append((s, u, v))
        cls.turns = [b for a, b, c in zip(cls.path, cls.path[1:], cls.path[2:])
                     if (b[1] - a[1]) * (c[1] - b[1]) < 0]

    def test_the_sphere_of_radiation_is_the_fixed_point(self):
        for value in self.flow(3 / 14, 3 / 14):
            self.assertAlmostEqual(value, 0, places=14)

    def test_small_departures_go_as_r_to_the_minus_three_quarters_times_a_cosine_of_root_47_over_4_ln_r(self):
        h, at = 1e-6, 3 / 14
        jacobian = [[(self.flow(at + h, at)[i] - self.flow(at - h, at)[i]) / (2 * h),
                     (self.flow(at, at + h)[i] - self.flow(at, at - h)[i]) / (2 * h)] for i in range(2)]
        trace = jacobian[0][0] + jacobian[1][1]
        determinant = jacobian[0][0] * jacobian[1][1] - jacobian[0][1] * jacobian[1][0]
        self.assertAlmostEqual(trace / 2, -3 / 4, places=8)
        self.assertAlmostEqual(math.sqrt(determinant - trace ** 2 / 4), math.sqrt(47) / 4, places=8)

    def test_the_mass_in_a_box_peaks_at_chavaniss_number(self):
        # 2Gm/(c^2 R) = 0.493 at alpha = sqrt(16 pi G rho_c/c^2) R = 4.7, where the density at the
        # centre is 22.4 times the density at the wall.
        s, u, v = self.turns[0]
        radius = math.exp(s)
        self.assertAlmostEqual(2 * u, 0.493, places=3)
        self.assertAlmostEqual(math.sqrt(16 * math.pi) * radius, 4.70, places=2)
        self.assertAlmostEqual(4 * math.pi * radius ** 2 / v, 22.4, places=1)
        self.assertTrue(all(2 * b[1] <= 2 * u for b in self.path))

    def test_the_mass_oscillates_about_three_sevenths_with_a_shrinking_swing(self):
        self.assertGreaterEqual(len(self.turns), 5)
        swings = [2 * u - 3 / 7 for s, u, v in self.turns]
        for first, second in zip(swings, swings[1:]):
            self.assertLess(first * second, 0)
            self.assertLess(abs(second), abs(first))
        # Far out, successive turns are half a period apart in ln r, 4 pi/sqrt(47), and each swing is
        # the last times -exp(-(3/4) 4 pi/sqrt(47)).
        gap = self.turns[4][0] - self.turns[3][0]
        self.assertAlmostEqual(gap, 4 * math.pi / math.sqrt(47), places=2)
        self.assertAlmostEqual(swings[4] / swings[3], -math.exp(-3 * math.pi / math.sqrt(47)), places=2)


if __name__ == "__main__":
    unittest.main()
