"""Horowitz and Myers's anti-de Sitter soliton:
.venv.noindex/bin/python -m unittest discover -s _tools

What its texts and drawings state, held on the published files. The numbers its captions quote,
the times light takes from the tip to the boundary, the period that makes the tip smooth, the
curvature at the tip and the energy of Horowitz and Myers's (3.15) and (3.16), are computed again
here, and the embedding diagram's sheet is held to the quadrature of its slope from the numbers
written and nothing else. Where sympy is installed, the published metric of every chart is held
to Einstein's equations with a negative cosmological constant, to the published black brane and
BTZ hole with their time and one direction exchanged, and to the other charts; those tests are
skipped where it is absent.
"""
import importlib.util
import json
import math
import sys
import types
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
HAS_SYMPY = importlib.util.find_spec("sympy") is not None
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))


def simpson(of, a, b, n=4000):
    h = (b - a) / n
    return h / 3 * (of(a) + of(b) + sum((4 if k % 2 else 2) * of(a + k * h) for k in range(1, n)))


def reach(n):
    """The time c t light takes from the tip to the boundary at r_0 = L = 1, the integral of
    dr/(r^2 sqrt(1 - r^-n)) from 1: with u = 1/r and u = 1 - s^2 it is the integral of
    2 s ds/sqrt(1 - (1 - s^2)^n) over 0 < s < 1, whose integrand is finite at both ends."""
    def of(s):
        return 2 / math.sqrt(n) if s == 0 else 2 * s / math.sqrt(1 - (1 - s * s) ** n)
    return simpson(of, 0.0, 1.0)


def lapse(rho):
    """sqrt(-g_tt) of the polar chart at r_0 = L = 1."""
    return math.cosh(1.5 * rho) ** (2 / 3)


def radius(rho):
    """sqrt(g_phiphi) of the polar chart at L = 1."""
    return 2 * math.sinh(1.5 * rho) / (3 * math.cosh(1.5 * rho) ** (1 / 3))


def slope(rho):
    """d radius/d rho = (2 cosh^2 + 1)/(3 cosh^(4/3)) of 3 rho/2."""
    c = math.cosh(1.5 * rho)
    return (2 * c * c + 1) / (3 * c ** (4 / 3))


class Numbers(unittest.TestCase):
    def test_light_reaches_the_boundary_in_a_beta_function_of_the_dimension(self):
        for n, quoted in ((3, "1.40"), (4, "1.31"), (2, "1.57")):
            exact = math.gamma(1 / n) * math.gamma(0.5) / (n * math.gamma(1 / n + 0.5))
            self.assertAlmostEqual(reach(n), exact, places=8)
            self.assertEqual(f"{exact:.2f}", quoted)
        self.assertAlmostEqual(reach(2), math.pi / 2, places=8)
        self.assertEqual(f"{2 * reach(3):.2f}", "2.80")
        self.assertEqual(f"{2 * reach(4):.2f}", "2.62")

    def test_light_round_the_boundary_circle_beats_light_through_the_tip(self):
        # Half the period 4 pi L^2/(3 r_0), against twice the time from the tip to the boundary.
        half_period = 2 * math.pi / 3
        self.assertEqual(f"{half_period:.2f}", "2.09")
        self.assertLess(half_period, 2 * reach(3))

    def test_the_period_makes_the_tip_the_origin_of_a_plane(self):
        # A circle of tau a small proper distance d from the tip has the circumference
        # beta (r/L) sqrt(f), and with beta = 4 pi L^2/(n r_0) that is 2 pi d.
        for n in (2, 3, 4):
            beta = 4 * math.pi / n
            for eps in (1e-4, 1e-5):
                r = 1 + eps
                f = 1 - r ** -n
                distance = simpson(lambda s: 2 / ((1 + s * s) * math.sqrt((1 - (1 + s * s) ** -n) / (s * s)))
                                   if s > 0 else 2 / math.sqrt(n), 0.0, math.sqrt(eps), 200)
                self.assertAlmostEqual(beta * r * math.sqrt(f) / (2 * math.pi * distance), 1.0, places=3)

    def test_the_polar_chart_is_horowitz_and_myerss_at_the_proper_distance(self):
        for rho in (0.1, 0.7, 2.0):
            r = math.cosh(1.5 * rho) ** (2 / 3)
            f = 1 - r ** -3
            self.assertAlmostEqual(lapse(rho), r, places=12)
            # (r/L) sqrt(f) d tau with tau = 2 L^2 phi/(3 r_0) is the circle of the polar chart.
            self.assertAlmostEqual(radius(rho), r * math.sqrt(f) * 2 / 3, places=12)
            # d rho = L dr/(r sqrt f).
            h = 1e-6
            dr = (math.cosh(1.5 * (rho + h)) ** (2 / 3) - math.cosh(1.5 * (rho - h)) ** (2 / 3)) / (2 * h)
            self.assertAlmostEqual(dr / (r * math.sqrt(f)), 1.0, places=7)

    def test_the_curvature_is_greatest_at_the_tip(self):
        # K = (12/L^4)(2 + r_0^6/r^6) in four dimensions and (8/L^4)(5 + 9 r_0^8/r^8) in five.
        self.assertEqual(12 * (2 + 1), 36)
        self.assertEqual(8 * (5 + 9), 112)
        self.assertEqual(12 * 2, 24)
        self.assertEqual(8 * 5, 40)

    def test_the_energy_of_horowitz_and_myers(self):
        # Their (3.15), E = -r_0^(p+1) beta V/(16 pi G l^(p+2)), with beta = 4 pi l^2/((p+1) r_0), is
        # their (3.16), E = -(V l^p/(16 pi G beta^p)) (4 pi/(p+1))^(p+1): negative, and as 1/beta^p.
        for p in (1, 2, 3, 5):
            for ell, r0 in ((1.0, 1.0), (1.3, 0.4)):
                beta = 4 * math.pi * ell ** 2 / ((p + 1) * r0)
                first = -r0 ** (p + 1) * beta / ell ** (p + 2)
                second = -ell ** p / beta ** p * (4 * math.pi / (p + 1)) ** (p + 1)
                self.assertAlmostEqual(first / second, 1.0, places=12)
                self.assertLess(first, 0)
        # Their (3.17) against (3.18): -pi^2 N^2/(8 beta^4) is 3/4 of -pi^2 N^2/(6 beta^4).
        self.assertAlmostEqual((1 / 8) / (1 / 6), 3 / 4, places=14)

    def test_a_particle_released_farther_out_takes_longer_over_each_swing(self):
        # A free particle released from rest at a has E = N(a), and with rho = a sin(theta) its
        # period in t is 4 times the integral of E a cos(theta)/(N sqrt(E^2 - N^2)).
        def period(a):
            E = lapse(a)

            def of(theta):
                if abs(theta - math.pi / 2) < 1e-12:
                    theta = math.pi / 2 - 1e-6
                N = lapse(a * math.sin(theta))
                return E * a * math.cos(theta) / (N * math.sqrt(E * E - N * N))
            return 4 * simpson(of, 0.0, math.pi / 2, 2000)
        small, half, one = period(0.01), period(0.5), period(1.0)
        # Near the tip N = 1 + 3 rho^2/4, so a small swing has the frequency sqrt(3/2).
        self.assertAlmostEqual(small, 2 * math.pi / math.sqrt(1.5), places=3)
        self.assertLess(small, half)
        self.assertLess(half, one)


class Surface(unittest.TestCase):
    def setUp(self):
        view, = json.loads((DATA / "embedding" / "ads_soliton.json").read_text(encoding="utf-8"))["views"]
        self.pieces = {p["id"]: p for p in view["surfaces"][0]["pieces"]}

    def test_the_sheet_stands_in_minkowski_space_at_the_quadrature_of_its_slope(self):
        sheet = self.pieces["sheet"]
        self.assertEqual(sheet["space"], "minkowski")
        points = sheet["points"]
        self.assertEqual(points[0], [0, 0, 0])
        self.assertEqual(points[-1][0], 2.0)
        for rho, R, Z in points[1::5] + [points[-1]]:
            self.assertLess(abs(R - radius(rho)), 2e-6)
            height = simpson(lambda x: math.sqrt(max(slope(x) ** 2 - 1, 0.0)), 0.0, rho, 400)
            self.assertLess(abs(Z - height), 2e-6, f"at {rho}")
        # Each chord is spacelike and as long as the proper distance between its circles, d rho.
        for (a, Ra, Za), (b, Rb, Zb) in zip(points, points[1:]):
            chord = math.sqrt((Rb - Ra) ** 2 - (Zb - Za) ** 2)
            self.assertLess(abs(chord - (b - a)), 2e-4 * (b - a), f"between {a} and {b}")

    def test_the_circles_outgrow_the_distance_everywhere_but_at_the_tip(self):
        self.assertAlmostEqual(slope(0.0), 1.0, places=14)
        for rho in (0.01, 0.3, 1.0, 2.0, 5.0):
            self.assertGreater(slope(rho), 1.0)

    def test_the_sheet_is_flat_at_the_tip_and_hyperbolic_far_out(self):
        # The curvature of a surface of revolution is -R''/R, here -tanh^2(3 rho/2)/L^2.
        h = 1e-4
        for rho in (0.05, 0.5, 1.5, 4.0):
            second = (radius(rho + h) - 2 * radius(rho) + radius(rho - h)) / h ** 2
            self.assertAlmostEqual(-second / radius(rho), -math.tanh(1.5 * rho) ** 2, places=5)
        self.assertLess(math.tanh(1.5 * 0.01) ** 2, 3e-4)
        self.assertAlmostEqual(math.tanh(1.5 * 8) ** 2, 1.0, places=9)

    def test_the_light_cone_drawn_is_the_one_the_sheet_nears(self):
        # R - Z tends to the integral of dR/drho - sqrt((dR/drho)^2 - 1), the depth of the apex.
        def of(x):
            return slope(x) - math.sqrt(max(slope(x) ** 2 - 1, 0.0))
        apex = simpson(of, 0.0, 30.0, 60000)
        cone = self.pieces["cone"]
        self.assertTrue(cone["reference"])
        self.assertLess(abs(cone["points"][0][2] + apex), 1e-5)
        self.assertIn(f"{apex:.2f}", cone["start"]["text"])
        for rho, R, Z in cone["points"]:
            self.assertLess(abs(R - Z - apex), 2e-6)
        # The sheet lies above the cone and closes on it.
        rho, R, Z = self.pieces["sheet"]["points"][-1]
        self.assertGreater(Z, R - apex)
        self.assertLess(Z - (R - apex), 0.2)


class Simplifier(unittest.TestCase):
    @unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
    def test_an_exponential_to_a_fractional_power_cancels_against_its_whole_powers(self):
        import sympy as sp
        import verify_metrics as vm
        x, L = sp.symbols("x L", positive=True)
        half = sp.cosh(3 * x / (2 * L))
        self.assertEqual(vm.norm(2 * half ** 2 - 1 - sp.cosh(3 * x / L)), 0)
        self.assertEqual(vm.norm(sp.sinh(3 * x / (2 * L)) ** 2 - half ** 2 + 1), 0)
        self.assertEqual(vm.norm(half ** sp.Rational(4, 3) / half ** sp.Rational(1, 3) - half), 0)
        self.assertNotEqual(vm.norm(2 * half ** 2 - sp.cosh(3 * x / L)), 0)
        # An expression with whole powers alone is handed back as it came.
        whole = sp.exp(2 * x / L) + sp.exp(-x / L)
        self.assertEqual(vm._whole_exponents(whole), (whole, {}))


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class FieldEquations(unittest.TestCase):
    def test_every_published_chart_is_the_soliton(self):
        import print_charts
        import verify_metrics as vm
        metric = json.loads((DATA / "metrics" / "ads_soliton.json").read_text(encoding="utf-8"))
        self.assertEqual([c["id"] for c in metric["coordinates"]], print_charts.SOLITON_CHARTS)
        for entry in metric["coordinates"]:
            reader = vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
            symbols = [reader.symbol[c] for c in entry["coords"]]
            n = len(symbols)
            published = {tuple(e["indices"]): reader(e["value"]) for e in entry["metric_components"]}
            g = vm.sp.Matrix(n, n, lambda i, j: published.get((entry["coords"][i], entry["coords"][j]), 0))
            chart = types.SimpleNamespace(geo=vm.Geometry(g, symbols, 600), reader=reader, symbols=symbols,
                                          coords_tex=entry["coords"])
            # Raises where the published metric misses R_ab = -((d - 1)/L^2) g_ab, the published
            # black brane or BTZ hole with its time and one direction exchanged, Horowitz and
            # Myers's chart under its map, or the plane at the tip of the polar chart.
            print_charts.ads_soliton_check(chart)


if __name__ == "__main__":
    unittest.main()
