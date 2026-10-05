"""Tests for Kerr-Taub-NUT's part of _tools/derivations: .venv.noindex/bin/python -m unittest discover -s _tools

The way print_charts.py writes its values, held to changing no value; its chart of Plebanski and
Demianski printed again byte for byte; and the tortoise coordinate its diagrams use on the regular
half of the axis. The tests need sympy and numpy, which the checker's own environment has, and are
skipped where they are absent.
"""
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path

DERIVATIONS = Path(__file__).resolve().parent / "derivations"
HAS_SYMPY = importlib.util.find_spec("sympy") is not None and importlib.util.find_spec("numpy") is not None


@unittest.skipUnless(HAS_SYMPY, "print_charts.py needs sympy")
class KerrTaubNut(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        sys.path.insert(0, str(DERIVATIONS))
        import print_charts
        import slices
        import sympy
        import verify_metrics
        cls.pc, cls.slices, cls.sp, cls.vm = print_charts, slices, sympy, verify_metrics

    @classmethod
    def tearDownClass(cls):
        sys.path.remove(str(DERIVATIONS))

    def forms(self, system):
        pc = self.pc
        reader = self.vm.Reader(["t", "r", "\\theta", "\\phi"], pc.kerr_taub_nut_parameters(system), ())
        return pc.KerrTaubNutForms(reader, north=system in pc.KTN_NORTH)

    def test_grouping_by_the_twist_changes_no_value(self):
        # chi = a sin^2 - 2 l cos with two strings and a sin^2 + 2 l (1 - cos) with one: a sum grouped
        # by its powers, with chi written out again, is the sum it came from.
        sp = self.sp
        for system in ("boyer_lindquist", "one_string"):
            forms = self.forms(system)
            r, a, l, m, C = forms.r, forms.a, forms.l, forms.m, forms.C
            twist = forms.twist
            for base in (r * (r ** 2 + a ** 2 + l ** 2) * twist ** 2 - 3 * m * a * twist * C + l * r ** 2,
                         twist ** 3 * (l + a * C) + (r - m) * twist * C + a * r):
                grouped = forms.by_twist(sp.expand(base))
                self.assertIsNotNone(grouped, system)
                written = grouped.subs({forms.chi: twist, forms.P: l + a * C})    # the names written out
                self.assertEqual(sp.expand(written - base), 0, system)
            # A sum whose coefficients would need a division by a is left to the other forms.
            self.assertIsNone(forms.by_twist(sp.expand((r - m) * (a * C ** 4 - 2 * l * C ** 3))))

    def test_the_twist_vanishes_on_the_regular_half_axis(self):
        sp = self.sp
        both, north = self.forms("boyer_lindquist"), self.forms("one_string")
        self.assertEqual(sp.expand(north.twist.subs(north.C, 1)), 0)
        self.assertEqual(sp.expand(north.twist.subs(north.C, -1) - 4 * north.l), 0)
        self.assertEqual(sp.expand(both.twist.subs(both.C, 1) + 2 * both.l), 0)
        self.assertEqual(sp.expand(north.twist - both.twist.subs(both.C, north.C) - 2 * north.l), 0)

    def test_the_tortoise_coordinate_of_the_regular_half_axis(self):
        # dr_*/dr = (r^2 + (a + l)^2)/Delta at m = 1, a = 1, l = 5/4, zero at r = 0.
        m, a, l = self.slices.KTN
        star = self.slices.ktn_rstar()
        self.assertAlmostEqual(float(star(0.0)), 0.0, places=14)
        for r in (-3.0, -1.0, 1.0, 3.0, 7.5):
            slope = float(star(r + 1e-6) - star(r - 1e-6)) / 2e-6
            self.assertAlmostEqual(slope, (r * r + (a + l) ** 2) / (r * r - 2 * m * r + a * a - l * l), places=6)

    def test_plebanski_and_demianskis_chart_prints_again_as_published(self):
        pc = self.pc
        published = pc.METRICS / "kerr_taub_nut.json"
        with tempfile.TemporaryDirectory() as folder:
            copy = Path(folder) / published.name
            copy.write_bytes(published.read_bytes())
            kept, pc.METRICS = pc.METRICS, Path(folder)
            try:
                pc.write(pc.kerr_taub_nut("plebanski"))
            finally:
                pc.METRICS = kept
            self.assertEqual(copy.read_bytes(), published.read_bytes())


if __name__ == "__main__":
    unittest.main()
