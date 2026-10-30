"""The global four-torus identifications in the two published five-dimensional charts.

The shared checker covers the local tensors and the generators cover the drawings.
These tests read the published coordinate ranges, which fix the global quotient.
"""
import json
from pathlib import Path
import sys
import unittest

import sympy as sp

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "_tools" / "derivations"))
import verify_metrics as vm


class CompactUniverse(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = ROOT / "MFS/assets/data/metrics/chodos_detweiler.json"
        cls.page = json.loads(path.read_text())
        cls.charts = {entry["id"]: entry for entry in cls.page["coordinates"]}

    def periods(self, chart):
        reader = vm.Reader(chart["coords"], [p["symbol"] for p in chart["parameters"]], ())
        periods = []
        for coordinate in chart["coords"][1:]:
            prefix = coordinate + " \\in [0, "
            domain, = [d for d in chart["domains"] if d.startswith(prefix)]
            self.assertTrue(domain.endswith(")"))
            periods.append(reader(domain[len(prefix):-1]))
        return reader, periods

    def test_both_charts_keep_four_compact_spatial_directions(self):
        self.assertEqual(set(self.charts), {"synchronous", "rescaled"})
        self.assertEqual(self.charts["synchronous"]["coords"], ["t", "x", "y", "z", "w"])
        self.assertEqual(self.charts["rescaled"]["coords"], ["t", "X", "Y", "Z", "W"])
        for chart in self.charts.values():
            self.assertIn("t \\in (0, \\infty)", chart["domains"])
            reader, periods = self.periods(chart)
            metric = vm.metric_from_line_element(reader, chart["line_element"], chart["coords"])
            for coordinate, period in zip(chart["coords"][1:], periods):
                self.assertNotEqual(period, 0)
                self.assertNotIn(reader.symbol[coordinate], metric.free_symbols)

    def test_rescaled_periods_follow_the_same_identifications(self):
        original, old_periods = self.periods(self.charts["synchronous"])
        self.assertEqual(old_periods, [original.parameters["L"]]*4)
        reader, periods = self.periods(self.charts["rescaled"])
        t0, tau, length = [reader.parameters[k] for k in ("t_0", "tau", "L")]
        for t0_value, tau_value, length_value in ((2, 5, 7), (7, 2, 11), (3, 7, 13)):
            at = {t0: sp.Rational(t0_value), tau: sp.Rational(tau_value), length: sp.Rational(length_value)}
            spatial_map = [sp.sqrt(at[tau]/at[t0])]*3 + [sp.sqrt(at[t0]/at[tau])]
            self.assertEqual([period.subs(at) for period in periods],
                             [at[length]*factor for factor in spatial_map])


if __name__ == "__main__":
    unittest.main()
