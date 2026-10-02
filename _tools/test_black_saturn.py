"""Black Saturn, a black ring in balance around a black hole in five dimensions:
python3 -m unittest discover -s _tools

Its charts leave four functions free, and its drawings declare Elvang and Figueras's on the plane
of the ring, where they are rational in z. These tests hold what the drawings and their captions
state to the published files: the Saturn drawn is in balance, the declared functions are the
limit on the plane of the functions the page defines, the frames are dragged at the ring's and
the hole's angular velocities on their horizons, and the embedding diagram's circles have the
radii its metric says. The tests of the declared functions need sympy and numpy and are skipped
where either is absent.
"""
import importlib.util
import json
import math
import sys
import unittest
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
HAS_TOOLS = all(importlib.util.find_spec(m) is not None for m in ("sympy", "numpy", "scipy", "contourpy"))
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))

K1, K2, K3 = Fraction(7, 8), Fraction(9, 16), Fraction(3, 7)


class Balance(unittest.TestCase):
    def test_the_saturn_drawn_is_in_balance_with_a_hole_of_no_spin(self):
        """Elvang and Figueras's balance condition at beta = 0:
        (kappa_1 - kappa_2)^2 = kappa_1 (1 - kappa_2)(1 - kappa_3)(kappa_1 - kappa_3)."""
        self.assertEqual((K1 - K2) ** 2, K1 * (1 - K2) * (1 - K3) * (K1 - K3))
        self.assertTrue(0 < K3 < K2 < K1 < 1)
        self.assertEqual(2 * K1 * K2 / K3, Fraction(147, 64))

    def test_the_horizons_turn_at_the_velocities_the_captions_state(self):
        """Their angular velocities at beta = 0, in units of c/L: sqrt(kappa_2 kappa_3/(2 kappa_1))
        for the hole and sqrt(kappa_1 kappa_3/(2 kappa_2)) for the ring."""
        self.assertAlmostEqual(math.sqrt(K2 * K3 / (2 * K1)), 3 * math.sqrt(3) / 14, places=14)
        self.assertAlmostEqual(math.sqrt(K1 * K3 / (2 * K2)), 1 / math.sqrt(3), places=14)


class Embedding(unittest.TestCase):
    def setUp(self):
        view = json.loads((DATA / "embedding" / "black_saturn.json").read_text(encoding="utf-8"))["views"][0]
        self.pieces = {p["id"]: p["points"] for p in view["surfaces"][0]["pieces"]}

    def test_each_circle_has_the_radius_of_its_circle_of_psi(self):
        for z, rho, _ in self.pieces["outside"]:
            q = 939 - 2527 * z + 2345 * z ** 2 - 784 * z ** 3
            self.assertAlmostEqual(rho, math.sqrt(2 * q / (49 * (1 - z) * (9 - 16 * z))), places=6)
        for z, rho, _ in self.pieces["between"]:
            self.assertAlmostEqual(rho, math.sqrt(14 * z * (1 - z) / (7 * z - 3)), places=6)

    def test_the_surface_runs_from_the_flat_plane_to_the_hole(self):
        outside, between = self.pieces["outside"], self.pieces["between"]
        self.assertEqual(outside[0][0], -8.0)
        self.assertAlmostEqual(outside[0][1], 4.354, delta=1e-3)        # sqrt(-2z) L = 4 L in the flat plane
        self.assertAlmostEqual(outside[-1][0], 0.242, delta=1e-3)       # where flat space stops carrying it
        self.assertAlmostEqual(between[0][0], 0.787, delta=1e-3)
        self.assertEqual(between[-1][0], 0.875)
        self.assertAlmostEqual(between[-1][1], 0.7, places=6)           # the hole's circle in the plane
        # The two circles where the surface stops stand at one height, and it lies level there.
        self.assertAlmostEqual(outside[-1][2], between[0][2], places=6)
        for points in (outside[-3:], between[:3]):
            (_, r0, h0), (_, r1, h1) = points[0], points[-1]
            self.assertLess(abs(h1 - h0), 0.1 * abs(r1 - r0))


@unittest.skipUnless(HAS_TOOLS, "sympy, numpy, scipy or contourpy is not installed")
class PlaneOfTheRing(unittest.TestCase):
    """The functions the drawings declare against the ones the page defines."""

    @classmethod
    def setUpClass(cls):
        import sympy
        import mpmath
        import null_rays
        import print_charts
        import verify_metrics
        cls.sp, cls.mp, cls.nr = sympy, mpmath, null_rays
        metric = json.loads((DATA / "metrics" / "black_saturn.json").read_text(encoding="utf-8"))
        entry = next(c for c in metric["coordinates"] if c["id"] == "weyl")
        reader = verify_metrics.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
        values = {"kappa_1": sympy.Rational(7, 8), "kappa_2": sympy.Rational(9, 16), "kappa_3": sympy.Rational(3, 7),
                  "beta": sympy.Integer(0)}
        cls.published, _ = print_charts.black_saturn_functions(reader, values, precision=60)

    def declared(self, functions, z):
        sp = self.sp
        return {n: sp.N(sp.sympify(text, locals={"z": sp.Symbol("z")}).subs(sp.Symbol("z"), sp.Rational(z)), 40)
                for n, text in functions.items()}

    def check(self, functions, places):
        mp = self.mp
        for z in places:
            declared = self.declared(functions, z)
            for name, value in declared.items():
                published = self.published[name](mp.mpf(10) ** -15, mp.mpf(self.sp.Rational(z).p) / self.sp.Rational(z).q)
                self.assertLess(abs(mp.mpmathify(value) - published), mp.mpf(10) ** -20, f"{name} at z = {z}")

    def test_outside_the_ring(self):
        self.check(self.nr.BS_OUTSIDE, ("-5", "-1/3", "1/5", "2/5"))

    def test_between_the_ring_and_the_hole(self):
        self.check(self.nr.BS_BETWEEN, ("3/5", "7/10", "17/20"))

    def test_the_frames_are_dragged_at_each_horizon_s_own_velocity(self):
        sp = self.sp
        z = sp.Symbol("z")
        at = lambda functions, value: sp.simplify(sp.sympify(functions["Omega"], locals={"z": z}).subs(z, value))
        self.assertEqual(at(self.nr.BS_OUTSIDE, sp.Rational(3, 7)), 1 / sp.sqrt(3))
        self.assertEqual(at(self.nr.BS_BETWEEN, sp.Rational(9, 16)), 1 / sp.sqrt(3))
        self.assertEqual(at(self.nr.BS_BETWEEN, sp.Rational(7, 8)), 3 * sp.sqrt(3) / 14)

    def test_the_gap_lies_inside_the_ergoregion_and_the_ergosurface_is_z_zero(self):
        sp = self.sp
        z = sp.Symbol("z")
        def gtt(functions):
            f = {n: sp.sympify(t, locals={"z": z}) for n, t in functions.items()}
            return sp.simplify(-sp.exp(2 * f["W"] - 2 * f["V"]) + f["Omega"] ** 2 * sp.exp(2 * f["V"]))
        outside, between = gtt(self.nr.BS_OUTSIDE), gtt(self.nr.BS_BETWEEN)
        self.assertEqual(sp.simplify(outside.subs(z, 0)), 0)
        self.assertLess(outside.subs(z, -1), 0)
        self.assertGreater(outside.subs(z, sp.Rational(1, 5)), 0)
        self.assertEqual(sp.simplify(between - (896 * z ** 2 - 1483 * z + 615) / (128 * (1 - z) * (7 * z - 3))), 0)
        for value in (sp.Rational(3, 5), sp.Rational(7, 10), sp.Rational(17, 20)):
            self.assertGreater(between.subs(z, value), 0)


if __name__ == "__main__":
    unittest.main()
