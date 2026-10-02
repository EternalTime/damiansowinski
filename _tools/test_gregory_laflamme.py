"""The Gregory-Laflamme instability of the black string, as _tools/derivations/gregory_laflamme.py
computes it: Gregory's perturbation equations against the linearised vacuum equations, and the
growth rates against Gregory and Laflamme's threshold. It needs sympy and scipy, so it runs under
/tmp/mfs-venv/bin/python and is skipped under a Python without them."""
import math
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))
try:
    import gregory_laflamme as gl
except ImportError:
    gl = None


@unittest.skipIf(gl is None, "needs sympy and scipy")
class GregoryLaflamme(unittest.TestCase):
    def test_the_perturbation_equations_solve_the_linearised_vacuum_equations(self):
        components = gl.linearised_ricci()
        self.assertEqual(len(components), 15)
        self.assertEqual({slot: value for slot, value in components.items() if value != 0}, {})

    def test_the_threshold_is_gregory_and_laflammes(self):
        self.assertEqual(round(gl.threshold(), 3), gl.THRESHOLD)
        # The threshold falls as the growth rate it is measured at rises.
        self.assertLess(gl.threshold(1e-2), gl.threshold(1e-4))

    def test_only_long_ripples_grow(self):
        # A wavelength of 10 r_s grows, at the rate the embedding diagram plays, and one of 7 r_s,
        # shorter than the critical 7.17 r_s, does not.
        self.assertEqual(round(gl.growth_rate(2 * math.pi / 10), 4), 0.0633)
        self.assertIsNone(gl.growth_rate(2 * math.pi / 7))
        self.assertGreater(gl.growth_rate(0.35), gl.growth_rate(0.6))
        self.assertGreater(gl.growth_rate(0.35), gl.growth_rate(0.15))


if __name__ == "__main__":
    unittest.main()
