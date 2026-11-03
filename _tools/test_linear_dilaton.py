"""Frame transformation and curvature of the linear dilaton cosmology."""

import json
from pathlib import Path
import sys
import unittest

import sympy as sp

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "_tools" / "derivations"))
import verify_metrics as vm


class LinearDilaton(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        page = json.loads((ROOT / "MFS/assets/data/metrics/linear_dilaton.json").read_text())
        cls.charts = {chart["id"]: chart for chart in page["coordinates"]}

    @staticmethod
    def read_chart(chart):
        reader = vm.Reader(chart["coords"], [p["symbol"] for p in chart["parameters"]], ())
        metric = vm.metric_from_line_element(reader, chart["line_element"], chart["coords"])
        return reader, sp.Matrix(metric)

    def test_ricci_scalar(self):
        chart = self.charts["cosmic"]
        reader, _ = self.read_chart(chart)
        published = reader(chart["ricci_scalar"].split("=", 1)[1])
        self.assertEqual(sp.simplify(published - 6 / reader.symbol[r"\tau"]**2), 0)

    def test_string_frame_factor_and_cosmic_time(self):
        cosmic_reader, cosmic = self.read_chart(self.charts["cosmic"])
        conformal_reader, conformal = self.read_chart(self.charts["conformal"])
        eta = conformal_reader.symbol[r"\eta"]
        q = conformal_reader.parameters["q"]
        tau = cosmic_reader.symbol[r"\tau"]
        cosmic = cosmic.subs({cosmic_reader.parameters["q"]: q, tau: sp.exp(q*eta)/q})
        jacobian = sp.diag(sp.exp(q*eta), 1, 1, 1)
        pulled = jacobian.T * cosmic * jacobian
        string_frame = sp.diag(-1, 1, 1, 1)
        expected = sp.exp(2*q*eta) * string_frame
        for entry in list(conformal - expected) + list(pulled - expected):
            self.assertEqual(sp.simplify(entry), 0)


if __name__ == "__main__":
    unittest.main()
