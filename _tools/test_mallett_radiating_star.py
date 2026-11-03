"""Flux sign and marginal sphere thresholds for the radiating de Sitter geometry."""
import json
from pathlib import Path
import sys
import unittest

import sympy as sp

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "_tools" / "derivations"))
import verify_metrics as vm


class RadiatingDeSitter(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        data = json.loads((ROOT / "MFS/assets/data/metrics/mallett_radiating_star.json").read_text())
        cls.charts = {c["id"]: c for c in data["coordinates"]}

    def test_published_flux_has_the_orientation_and_cosmological_term(self):
        for name, time, sign in (("eddington_finkelstein_outgoing", "u", -1),
                                 ("eddington_finkelstein_ingoing", "v", 1)):
            chart = self.charts[name]
            reader = vm.Reader(chart["coords"], [p["symbol"] for p in chart["parameters"]], ())
            ricci = {tuple(x["indices"]): reader(x["value"]) for x in chart["ricci_tensor"]["variants"]["ll"]["nonzero"]}
            metric = {tuple(x["indices"]): reader(x["value"]) for x in chart["metric_components"]}
            w, r = reader.symbol[time], reader.symbol["r"]
            mass, lam = reader.parameters["m"], reader.parameters["Lambda"]
            derivative = sp.Derivative(mass, w)
            for radius, m0, rate, cosm in ((2, sp.Rational(1, 2), sp.Rational(1, 20), sp.Rational(1, 10)),
                                           (3, sp.Rational(2, 3), -sp.Rational(1, 10), sp.Rational(1, 50))):
                at = {r: radius, w: 0, mass: m0, derivative: rate, lam: cosm,
                      reader.symbol[r"\theta"]: sp.Rational(1, 2), reader.symbol[r"\phi"]: 0}
                for i in chart["coords"]:
                    for j in chart["coords"]:
                        value = ricci.get((i, j), 0) - lam * metric.get((i, j), 0)
                        expected = sign * 2 * rate / radius**2 if i == j == time else 0
                        self.assertAlmostEqual(float(sp.N(sp.sympify(value).subs(at), 20)), float(expected), places=12)

    def test_subcritical_degenerate_and_supercritical_roots(self):
        # f has one positive maximum; its sign fixes the number of positive roots.
        for mass, cosm, expected in ((sp.Rational(1, 2), sp.Rational(1, 10), 1),
                                      (sp.Rational(1, 3), 1, 0),
                                      (sp.Rational(1, 2), 1, -1)):
            critical_r = (3 * mass / cosm) ** sp.Rational(1, 3)
            peak = 1 - 2 * mass / critical_r - cosm * critical_r**2 / 3
            self.assertEqual(int(sp.sign(peak)), expected)
            self.assertEqual(int(sp.sign(1 - 9 * cosm * mass**2)), expected)


if __name__ == "__main__":
    unittest.main()
