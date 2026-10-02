"""A root above the square in the checker's canonical form, and a function raised to a fractional
power in its reader, as Senovilla's cosh^(-2/3)(3a rho) needs both, and the physics his chart is
published with. Needs sympy, so it is skipped under a Python without it."""
import json
import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "derivations"))
try:
    import sympy as sp
    import verify_metrics as vm
except ImportError:
    sp = None


@unittest.skipIf(sp is None, "sympy is not installed")
class Roots(unittest.TestCase):
    def test_the_reader_takes_a_function_raised_to_a_negative_or_fractional_power(self):
        reader = vm.Reader(["t", "\\rho", "\\phi", "z"], ["a"], ())
        a, rho = reader.parameters["a"], reader.symbol["\\rho"]
        self.assertEqual(reader("\\cosh^{-2/3}(3a\\rho)"), sp.cosh(3 * a * rho) ** sp.Rational(-2, 3))
        self.assertEqual(reader("\\cosh^{-2}(3a\\rho)"), sp.cosh(3 * a * rho) ** -2)
        self.assertEqual(reader("\\cosh^2(3a\\rho)"), sp.cosh(3 * a * rho) ** 2)

    def test_a_cube_root_is_reduced_exactly(self):
        x = sp.Symbol("x", positive=True)
        w = sp.cosh(x) ** sp.Rational(1, 3)
        self.assertEqual(vm.norm(w ** 3 - sp.cosh(x)), 0)
        self.assertEqual(vm.norm(1 / w - w ** 2 / sp.cosh(x)), 0)
        self.assertEqual(vm.norm(sp.diff(w, x) - sp.sinh(x) / (3 * w ** 2)), 0)
        self.assertNotEqual(vm.norm(w ** 2 - w), 0)
        y = sp.Symbol("y", positive=True)
        self.assertEqual(vm.norm((x * y) ** sp.Rational(2, 3) - x ** sp.Rational(2, 3) * y ** sp.Rational(2, 3)), 0)

    def test_a_number_in_the_denominator_of_an_exponent_is_taken_out(self):
        # exp(x/(8(u + v))) squared is exp(x/(4(u + v))), however expand has spread the 8 or the 4
        # through the sum, as the conformal factor of Bonnor's dust cloud needs.
        x, u, v = sp.symbols("x u v", positive=True)
        base = sp.exp(x / (8 * (u + v) ** 2))
        self.assertEqual(vm.norm(base ** 2 - sp.exp(x / (4 * (u + v) ** 2))), 0)
        self.assertEqual(vm.norm(base * sp.exp(-x / (8 * u ** 2 + 16 * u * v + 8 * v ** 2)) - 1), 0)
        self.assertEqual(vm.norm(sp.diff(base, x) - base / (8 * u ** 2 + 16 * u * v + 8 * v ** 2)), 0)
        self.assertNotEqual(vm.norm(base - sp.exp(x / (4 * (u + v) ** 2))), 0)

    def test_a_square_root_is_reduced_as_before(self):
        x = sp.Symbol("x", positive=True)
        self.assertEqual(vm.norm(1 / (1 + sp.sqrt(x)) - (1 - sp.sqrt(x)) / (1 - x)), 0)


@unittest.skipIf(sp is None, "sympy is not installed")
class Senovilla(unittest.TestCase):
    """The published chart: radiation with a third of its density for pressure, a density greatest
    at the bounce on the axis, and the one circular light ray at cosh(3a rho) = 2."""

    def setUp(self):
        metric = json.loads((HERE.parent / "MFS/assets/data/metrics/senovilla.json").read_text(encoding="utf-8"))
        (self.chart,) = metric["coordinates"]
        self.reader = vm.Reader(self.chart["coords"], [p["symbol"] for p in self.chart["parameters"]], ())
        self.t, self.rho = self.reader.symbol["t"], self.reader.symbol["\\rho"]
        self.at = {self.reader.parameters["a"]: 1, self.reader.c: 1}

    def components(self, field, variant=None):
        block = self.chart[field]["variants"][variant]["nonzero"] if variant else self.chart[field]
        return {tuple(c["indices"]): self.reader(c["value"]).subs(self.at) for c in block}

    def test_the_fluid_is_radiation_and_its_density_is_greatest_at_the_bounce_on_the_axis(self):
        G = self.components("einstein_tensor", "ul")
        density = -G[("t", "t")]
        for x in ("\\rho", "\\phi", "z"):
            self.assertEqual(sp.simplify(G[(x, x)] - density / 3), 0)
        self.assertEqual(density.subs({self.t: 0, self.rho: 0}), 15)
        for t, rho in ((0.3, 0), (0, 0.2), (-1, 1.5), (2, 0.1)):
            value = float(density.subs({self.t: t, self.rho: rho}))
            self.assertTrue(0 < value < 15)

    def test_the_one_circular_light_ray_is_at_cosh_3a_rho_equal_to_2(self):
        g = self.components("metric_components")
        gamma = self.components("christoffel", "ull")
        turning = gamma[("\\rho", "t", "t")] - gamma[("\\rho", "\\phi", "\\phi")] * g[("t", "t")] / g[("\\phi", "\\phi")]
        ring = sp.acosh(2) / 3
        for t in (-1, 0, sp.Rational(1, 2)):
            self.assertAlmostEqual(float(turning.subs({self.t: t, self.rho: ring})), 0, places=12)
            self.assertLess(float(turning.subs({self.t: t, self.rho: ring / 2})), 0)
            self.assertGreater(float(turning.subs({self.t: t, self.rho: 2 * ring})), 0)


if __name__ == "__main__":
    unittest.main()
