"""Tests for _tools/derivations/print_charts.py: python3 -m unittest discover -s _tools

A chart the script wrote is printed again and held to the published file byte for byte, so a
change to the printer that reorders a published value is caught here and not by a reader. The
tests need sympy, which the checker's own environment has, and are skipped where it is absent.
"""
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path

DERIVATIONS = Path(__file__).resolve().parent / "derivations"
HAS_SYMPY = importlib.util.find_spec("sympy") is not None


@unittest.skipUnless(HAS_SYMPY, "print_charts.py needs sympy")
class Reprint(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        sys.path.insert(0, str(DERIVATIONS))
        import print_charts
        cls.pc = print_charts

    @classmethod
    def tearDownClass(cls):
        sys.path.remove(str(DERIVATIONS))

    def reprinted(self, metric_id):
        """The file print_charts.py writes for one spacetime, printed into a copy of its own."""
        pc = self.pc
        published = pc.METRICS / f"{metric_id}.json"
        with tempfile.TemporaryDirectory() as folder:
            copy = Path(folder) / published.name
            copy.write_bytes(published.read_bytes())
            kept, pc.METRICS = pc.METRICS, Path(folder)
            try:
                builders = pc.CHARTS[metric_id]
                for build in builders if isinstance(builders, list) else [builders]:
                    specs = build()
                    for spec in specs if isinstance(specs, list) else [specs]:
                        pc.write(spec)
            finally:
                pc.METRICS = kept
            return published.read_bytes(), copy.read_bytes()

    def test_tov_prints_again_as_published(self):
        # Its collected bracket closes each product, r(r - 2m)((d_r Phi)^2 + d_r^2 Phi), which a
        # change to where the printer writes a named factor moved to the front on 30 September 2026.
        published, printed = self.reprinted("tov")
        self.assertEqual(printed, published)

    def test_a_chart_prints_the_same_twice(self):
        # The random point block() tells values apart at is seeded, so two prints agree.
        first = self.reprinted("einstein_static")[1]
        self.assertEqual(self.reprinted("einstein_static")[1], first)


if __name__ == "__main__":
    unittest.main()
