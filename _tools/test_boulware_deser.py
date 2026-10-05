"""Boulware and Deser's black hole of Einstein-Gauss-Bonnet gravity in five dimensions:
.venv.noindex/bin/python -m unittest discover -s _tools

What its texts and drawings state, held on the published files. The embedding diagram's two
surfaces are held to the quadrature of their slopes, from the numbers written and nothing else,
and the numbers the captions quote are computed again here. Where sympy is installed, the
published metric of every chart is held to the field equations of Einstein-Gauss-Bonnet gravity,
G_ab + (l^2/2) H_ab = 0, and to the horizon, the surface gravity and the centre the History
states; those tests are skipped where it is absent.
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

# The black hole as its diagrams draw it, in units of its horizon radius, and the other branch in units of l.
R0, ELL = 13 / 12, 5 / 12
PLUS_R0 = 1.0


def simpson(of, a, b, n=4000):
    h = (b - a) / n
    return h / 3 * (of(a) + of(b) + sum((4 if k % 2 else 2) * of(a + k * h) for k in range(1, n)))


def f_minus(r, r0=R0, ell=ELL):
    return 1 - 2 * r0 ** 2 / (r * r + math.sqrt(r ** 4 + 4 * ell ** 2 * r0 ** 2))


def f_plus(r, r0=PLUS_R0, ell=1.0):
    return 1 + (r * r + math.sqrt(r ** 4 + 4 * ell ** 2 * r0 ** 2)) / (2 * ell ** 2)


class Numbers(unittest.TestCase):
    def test_the_two_ways_of_writing_the_metric_function_agree(self):
        for r in (0.2, 0.9, 1.0, 1.7, 6.0):
            line = 1 + r * r / (2 * ELL ** 2) * (1 - math.sqrt(1 + 4 * ELL ** 2 * R0 ** 2 / r ** 4))
            self.assertAlmostEqual(line, f_minus(r), places=12)

    def test_the_horizon_the_centre_and_the_surface_gravity_of_the_hole_drawn(self):
        self.assertAlmostEqual(math.sqrt(R0 ** 2 - ELL ** 2), 1.0, places=14)
        self.assertAlmostEqual(f_minus(1.0), 0.0, places=14)
        self.assertAlmostEqual(f_minus(1e-9), 1 - R0 / ELL, places=9)
        self.assertAlmostEqual(1 - R0 / ELL, -8 / 5, places=14)
        h = 1e-6
        self.assertAlmostEqual((f_minus(1 + h) - f_minus(1 - h)) / (4 * h), 72 / 97, places=8)
        self.assertAlmostEqual(1 / (1 + 2 * ELL ** 2), 72 / 97, places=14)

    def test_the_hole_drawn_is_above_the_mass_below_which_small_perturbations_grow(self):
        # Beroiz, Dotti and Gleiser's range 3/2 < mu/alpha < 9/2 + 3 sqrt 2 is 1 < r_0/l < 1 + sqrt 2
        # with r_0^2 = 2 mu/3 and l^2 = alpha, and the specific heat is positive for r_h < sqrt 2 l.
        self.assertAlmostEqual(math.sqrt(2 / 3 * (4.5 + 3 * math.sqrt(2))), 1 + math.sqrt(2), places=12)
        self.assertGreater(R0 / ELL, 1 + math.sqrt(2))
        self.assertLess(math.sqrt(2 + 1), 1 + math.sqrt(2))

    def test_the_numbers_the_captions_quote(self):
        self.assertEqual(round(1 / f_minus(2.0), 2), 1.41)
        self.assertAlmostEqual(1 / (1 - 1 / 4), 4 / 3, places=14)
        self.assertAlmostEqual(2 / (1 - R0 / ELL), -5 / 4, places=12)
        self.assertAlmostEqual(f_plus(0.0), 2.0, places=14)
        # R, the tortoise coordinate of the other branch at infinity: in 1/s beyond s = 1.
        R = simpson(lambda s: 1 / f_plus(s), 0, 1) + simpson(
            lambda w: 1 / (f_plus(1 / w) * w * w) if w > 0 else 1.0, 0, 1)
        self.assertEqual(round(R, 4), 1.1981)


class Surfaces(unittest.TestCase):
    def setUp(self):
        self.views = {v["id"]: v for v in json.loads(
            (DATA / "embedding" / "boulware_deser.json").read_text(encoding="utf-8"))["views"]}

    def pieces(self, view):
        return {p["id"]: p for p in self.views[view]["surfaces"][0]["pieces"]}

    def test_the_holes_slice_stands_at_the_quadrature_of_its_slope(self):
        def in_q(q):
            # dz/dr = sqrt(1/f - 1) with r = 1 + q^2, where it is finite at the throat.
            r = 1 + q * q
            return 2 * math.sqrt((2 * R0 ** 2 - r * r + math.sqrt(r ** 4 + 4 * ELL ** 2 * R0 ** 2)) / (2 * (2 + q * q)))
        for sign, pid in ((1, "exterior"), (-1, "other_exterior")):
            points = self.pieces("hole")[pid]["points"]
            self.assertAlmostEqual(points[0][0], 1.0, places=12)
            for r, rho, z in points[::7] + [points[-1]]:
                self.assertLess(abs(rho - r), 2e-6)
                self.assertLess(abs(z - sign * simpson(in_q, 0.0, math.sqrt(r - 1), 400)), 2e-5, f"{pid} at {r}")
        self.assertEqual(round(self.pieces("hole")["exterior"]["points"][-1][2], 2), 2.77)
        self.assertEqual(round(math.acosh(6.0), 2), 2.48)

    def test_the_other_branch_is_a_sheet_of_minkowski_space_that_leaves_the_centre_along_a_cone(self):
        sheet = self.pieces("branch")["sheet"]
        self.assertEqual(sheet["space"], "minkowski")
        points = sheet["points"]
        self.assertEqual(points[0], [0, 0, 0])
        for r, rho, z in points[1::5] + [points[-1]]:
            self.assertLess(abs(rho - r), 2e-6)
            self.assertLess(abs(z - simpson(lambda x: math.sqrt(1 - 1 / f_plus(x)), 0.0, r, 400)), 2e-6, f"at {r}")
        r, _, z = points[1]
        self.assertLess(abs(z / r - math.sqrt(1 / 2)), 2e-3)
        # Each chord is spacelike and as long as the proper distance between its circles.
        for (a, ra, za), (b, rb, zb) in zip(points, points[1:]):
            chord = math.sqrt((rb - ra) ** 2 - (zb - za) ** 2)
            proper = simpson(lambda x: 1 / math.sqrt(f_plus(x)), a, b, 20)
            self.assertLess(abs(chord - proper), 2e-4 * proper, f"between {a} and {b}")


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class FieldEquations(unittest.TestCase):
    def test_every_published_chart_solves_the_field_equations_of_einstein_gauss_bonnet_gravity(self):
        import print_charts
        import verify_metrics as vm
        metric = json.loads((DATA / "metrics" / "boulware_deser.json").read_text(encoding="utf-8"))
        self.assertEqual([c["id"] for c in metric["coordinates"]], print_charts.BOULWARE_DESER_CHARTS)
        for entry in metric["coordinates"]:
            reader = vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
            symbols = [reader.symbol[c] for c in entry["coords"]]
            n = len(symbols)
            published = {tuple(e["indices"]): reader(e["value"]) for e in entry["metric_components"]}
            g = vm.sp.Matrix(n, n, lambda i, j: published.get((entry["coords"][i], entry["coords"][j]), 0))
            chart = types.SimpleNamespace(geo=vm.Geometry(g, symbols, 600), reader=reader, symbols=symbols,
                                          coords_tex=entry["coords"])
            # Raises where the published metric misses the field equations, the horizon
            # r_h^2 = r_0^2 - l^2, the surface gravity r_h/(r_h^2 + 2 l^2), f(0) = 1 - r_0/l,
            # Tangherlini's metric at l -> 0, or r^4 K -> 12 r_0^2/l^2 at the centre.
            print_charts.boulware_deser_check(chart, entry["id"])


if __name__ == "__main__":
    unittest.main()
