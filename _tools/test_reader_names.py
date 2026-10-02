"""The checker's reader and the letter e: a power of e is the exponential, unless the system
declares e as a name of its own, as Ori's time machine does. Needs sympy, so it runs in the
environment the checker runs in and is skipped elsewhere."""
import importlib.util
import sys
import unittest
from pathlib import Path

HAS_SYMPY = importlib.util.find_spec("sympy") is not None
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class TheLetterE(unittest.TestCase):
    def setUp(self):
        import sympy
        import verify_metrics
        self.sp, self.vm = sympy, verify_metrics

    def test_a_power_of_e_is_the_exponential_where_no_name_is_declared(self):
        reader = self.vm.Reader(["t", "x"], ["a"], ())
        x = reader.symbol["x"]
        self.assertEqual(reader("e^{2x}"), self.sp.exp(2 * x))

    def test_a_power_of_a_declared_e_is_that_parameters(self):
        reader = self.vm.Reader(["t", "x"], ["a", "e"], ())
        e, a = reader.parameters["e"], reader.parameters["a"]
        self.assertEqual(self.sp.expand(reader("\\left(2e + a\\right)^2 - 4e^2")), self.sp.expand(4 * a * e + a ** 2))
        self.assertFalse(reader("e^2").has(self.sp.E))


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class ASubscriptedFunction(unittest.TestCase):
    """A run of \\partial takes the numeral subscript of a declared function with its name, as the
    1 of Wahlquist's h_1, and leaves a subscript that belongs to no declared function alone."""

    def setUp(self):
        import sympy
        import verify_metrics
        self.sp, self.vm = sympy, verify_metrics

    def test_a_partial_of_a_held_function_with_a_subscript(self):
        reader = self.vm.Reader(["t", "x"], ["h_1 = 1 + x^2", "h_2 = x^3"], (), held=("h_1", "h_2"))
        x, h1, h2 = reader.symbol["x"], reader.parameters["h_1"], reader.parameters["h_2"]
        self.assertEqual(reader("\\dfrac{\\partial_x h_1}{2h_2}"), self.sp.Derivative(h1, x) / (2 * h2))
        self.assertEqual(reader("\\partial_x^2 h_{2}"), self.sp.Derivative(h2, (x, 2)))
        self.assertEqual(reader("x\\,\\partial_x h_1").subs(reader.held).doit(), 2 * x ** 2)

    def test_a_function_without_a_subscript_keeps_what_follows_it(self):
        reader = self.vm.Reader(["t", "x"], ["m = m(x)", "x_0"], ())
        x, m, x0 = reader.symbol["x"], reader.parameters["m"], reader.parameters["x_0"]
        self.assertEqual(reader("\\partial_x m\\,x_0"), self.sp.Derivative(m, x) * x0)

    def test_the_inverse_sines(self):
        reader = self.vm.Reader(["t", "x"], ["k"], ())
        x, k = reader.symbol["x"], reader.parameters["k"]
        self.assertEqual(reader("\\arcsin(k x) + \\mathrm{arsinh}(k x)"), self.sp.asin(k * x) + self.sp.asinh(k * x))


if __name__ == "__main__":
    unittest.main()
