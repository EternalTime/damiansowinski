"""Scalar and electromagnetic matter on the Gowdy torus, with periodic initial data."""
import math
from pathlib import Path
import sys
import unittest

import numpy as np
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))
import chart_printer as cp
import null_rays as nr
import print_charts as pc
from test_build_mfs_data import checks


@checks("charach_malin")
class ToroidalMatter(unittest.TestCase):
    def test_each_chart_obeys_its_sourced_matter_equations(self):
        for system in pc.CHARACH_MALIN_CHARTS:
            with self.subTest(system=system):
                spec = pc.charach_malin(system)
                chart = cp.Chart(spec["system"]["coords"], spec["system"]["parameters"], spec["chart_line_element"])
                pc.charach_malin_check(chart, system)

    def test_logarithmic_metrics_are_areal_pullbacks(self):
        for matter in ("scalar", "electromagnetic"):
            specs = [pc.charach_malin(matter+"_"+time) for time in ("areal", "logarithmic")]
            charts = [cp.Chart(s["system"]["coords"], s["system"]["parameters"], s["chart_line_element"]) for s in specs]
            for tau in (-1.3, 0.2, 1.7):
                t = math.exp(-tau)
                matrices = []
                for chart, time in zip(charts, (t, tau)):
                    at = {chart.symbols[0]:time, chart.reader.parameters['L']:1.7,
                          chart.reader.parameters['P']:0.37, chart.reader.parameters['lambda']:-0.23}
                    matrices.append(np.array(chart.geo.g.subs(at), dtype=float))
                J = np.diag([-t, 1, 1, 1])
                np.testing.assert_allclose(J.T @ matrices[0] @ J, matrices[1], rtol=2e-14, atol=2e-14)
