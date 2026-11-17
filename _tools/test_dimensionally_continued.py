"""The dimensionally continued black holes of Banados, Teitelboim and Zanelli:
.venv.noindex/bin/python -m unittest _tools/test_dimensionally_continued.py

What the page's texts state about the family, held to the metric functions its History prints,
f = 1 - (2M/r)^(1/(n-1)) + r^2/l^2 in D = 2n dimensions and f = 1 - (M + 1)^(1/(n-1)) + r^2/l^2 in
D = 2n - 1, in numbers alone: the members at n = 2, the vacua, the gap of naked singularities, the
negative root that the 1994 article set aside, the temperatures and their limits, the least temperature
Muniain and Piriz found, the horizons the conventions state, and the values every drawing takes. The
History is held to printing the formulas these tests compute. The charts themselves, their field
equations and every drawing are checked by the generators that write them.
"""
import json
import math
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
METRIC = json.loads((ROOT / "MFS" / "assets" / "data" / "metrics" / "dimensionally_continued.json")
                    .read_text(encoding="utf-8"))
NS = (2, 3, 4, 5)
RADII = (0.05, 0.3, 0.9, 1.0, 1.7, 4.0, 25.0)


def even(r, M, n, ell=1.0, sign=1):
    """The metric function of D = 2n dimensions."""
    return 1 - sign * (2 * M / r) ** (1 / (n - 1)) + r * r / (ell * ell)


def odd(r, M, n, ell=1.0, sign=1):
    """The metric function of D = 2n - 1 dimensions; sign = -1 takes the negative root where 1/(n-1) is 1/2."""
    return 1 - sign * (M + 1) ** (1 / (n - 1)) + r * r / (ell * ell)


def slope(f, r, h=1e-6):
    return (f(r + h) - f(r - h)) / (2 * h)


def even_mass(rp, n, ell=1.0):
    """The mass whose even black hole has its horizon at rp: (2M/rp)^(1/(n-1)) = 1 + rp^2/l^2."""
    return rp * (1 + rp * rp / (ell * ell)) ** (n - 1) / 2


def odd_mass(rp, n, ell=1.0):
    return (1 + rp * rp / (ell * ell)) ** (n - 1) - 1


def even_temperature(rp, n, ell=1.0):
    """The History's T in units of hbar c/k_B."""
    return (1 + (2 * n - 1) * rp * rp / (ell * ell)) / (4 * math.pi * (n - 1) * rp)


def odd_temperature(rp, ell=1.0):
    return rp / (2 * math.pi * ell * ell)


class TheHistoryPrintsTheseFormulas(unittest.TestCase):
    def test_the_metric_functions_and_the_temperatures(self):
        h = METRIC["history"]
        for tex in ("f(r) = 1 - \\left(\\frac{2M}{r}\\right)^{1/(n-1)} + \\frac{r^2}{\\ell^2}",
                    "f(r) = 1 - (M + 1)^{1/(n-1)} + \\frac{r^2}{\\ell^2}",
                    "T = \\frac{\\hbar c}{4\\pi k_B}\\,\\frac{1 + (2n-1)\\,r_+^2/\\ell^2}{(n-1)\\,r_+}",
                    "T = \\hbar c\\,r_+/(2\\pi k_B\\ell^2)", "r_+ = \\ell/\\sqrt{D-1}"):
            self.assertIn(tex, h)


class Members(unittest.TestCase):
    def test_at_n_2_the_even_family_is_schwarzschild_anti_de_sitter_and_the_odd_one_the_black_hole_of_three(self):
        for r in RADII:
            for M in (0.3, 1.0, 2.5):
                self.assertAlmostEqual(even(r, M, 2, 1.3), 1 - 2 * M / r + r * r / 1.69, places=12)
                self.assertAlmostEqual(odd(r, M, 2, 1.3), -M + r * r / 1.69, places=12)

    def test_anti_de_sitter_space_is_M_minus_1_in_odd_dimensions_and_M_0_in_even_ones(self):
        for n in NS:
            for r in RADII:
                self.assertAlmostEqual(odd(r, -1.0, n), 1 + r * r, places=12)
                self.assertAlmostEqual(even(r, 0.0, n), 1 + r * r, places=12)

    def test_without_a_cosmological_constant_the_even_horizon_is_at_r_2M(self):
        for n in NS:
            for M in (0.5, 2.0):
                self.assertAlmostEqual(1 - (2 * M / (2 * M)) ** (1 / (n - 1)), 0.0, places=14)


class TheOddFamily(unittest.TestCase):
    def test_a_horizon_exists_for_M_at_least_0_and_none_in_the_gap_above_anti_de_sitter_space(self):
        rs = [k / 100 for k in range(1, 1000)]
        for n in NS:
            for M in (-0.99, -0.5, -0.01):
                self.assertTrue(all(odd(r, M, n) > 0 for r in rs), f"n = {n}, M = {M}")
            for M in (0.01, 0.5, 3.0):
                self.assertTrue(any(odd(r, M, n) < 0 for r in rs), f"n = {n}, M = {M}")

    def test_with_the_negative_root_of_odd_n_a_positive_mass_has_no_horizon(self):
        # n = 3: (M + 1)^(1/2) taken negative, f = 1 + sqrt(M + 1) + r^2 > 0.
        for M in (0.1, 3.0, 50.0):
            self.assertTrue(all(odd(r, M, 3, sign=-1) > 1 for r in RADII))


class Temperatures(unittest.TestCase):
    def test_the_historys_temperatures_are_f_prime_over_4_pi_at_the_horizon(self):
        for n in NS:
            for rp in (0.2, 0.7, 1.0, 2.3):
                M = even_mass(rp, n)
                self.assertAlmostEqual(even(rp, M, n), 0.0, places=12)
                self.assertAlmostEqual(slope(lambda r: even(r, M, n), rp) / (4 * math.pi),
                                       even_temperature(rp, n), places=7)
                M = odd_mass(rp, n)
                self.assertAlmostEqual(odd(rp, M, n), 0.0, places=12)
                self.assertAlmostEqual(slope(lambda r: odd(r, M, n), rp) / (4 * math.pi), odd_temperature(rp), places=7)

    def test_an_even_black_hole_grows_hotter_without_limit_as_it_shrinks_and_an_odd_one_cools_toward_zero(self):
        for n in NS:
            temps = [even_temperature(10.0 ** -k, n) for k in range(1, 8)]
            self.assertTrue(all(b > 9 * a for a, b in zip(temps, temps[1:])))
            self.assertLess(odd_temperature(1e-7), 1e-7)

    def test_the_even_temperature_is_least_at_r_plus_l_over_root_D_minus_1_and_the_odd_has_no_least(self):
        # The heat capacity dM/dT changes sign where dT/dr_+ does, M growing with r_+.
        for n in NS:
            least = 1 / math.sqrt(2 * n - 1)
            self.assertLess(even_temperature(least, n), even_temperature(least * 0.99, n))
            self.assertLess(even_temperature(least, n), even_temperature(least * 1.01, n))
            self.assertAlmostEqual(slope(lambda r: even_temperature(r, n), least), 0.0, places=8)
            self.assertGreater(even_mass(least * 1.01, n), even_mass(least, n))
            self.assertGreater(slope(odd_temperature, least), 0)


class TheChartsAndTheDrawings(unittest.TestCase):
    def test_the_horizons_the_conventions_state(self):
        for M in (0.2, 3.0, 15.0):
            for ell in (0.5, 1.0, 2.0):
                rp = ell * math.sqrt(math.sqrt(M + 1) - 1)
                self.assertAlmostEqual(odd(rp, M, 3, ell), 0.0, places=12)
        for rp in (0.3, 1.0, 2.0):
            for ell in (0.5, 1.0, 2.0):
                M = rp * (1 + rp * rp / (ell * ell)) ** 2 / 2
                self.assertAlmostEqual(even(rp, M, 3, ell), 0.0, places=12)

    def test_every_drawing_puts_the_horizon_at_l(self):
        self.assertAlmostEqual(odd(1.0, 3.0, 3), 0.0, places=14)
        self.assertAlmostEqual(even(1.0, 2.0, 3), 0.0, places=14)
        self.assertAlmostEqual(slope(lambda r: odd(r, 3.0, 3), 1.0) / 2, 1.0, places=8)
        self.assertAlmostEqual(slope(lambda r: even(r, 2.0, 3), 1.0) / 2, 1.5, places=8)

    def test_the_charts_are_of_five_and_six_dimensions_on_the_positive_root(self):
        charts = {c["id"]: c for c in METRIC["coordinates"]}
        self.assertEqual(sorted(len(c["coords"]) for c in charts.values()), [5, 5, 5, 6, 6])
        for system in ("spherical", "eddington_finkelstein_ingoing", "spherical_six",
                       "eddington_finkelstein_ingoing_six"):
            self.assertIn("positive square root", charts[system]["convention"])


if __name__ == "__main__":
    unittest.main()
