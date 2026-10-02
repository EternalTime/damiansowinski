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


if __name__ == "__main__":
    unittest.main()
