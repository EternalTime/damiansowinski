"""The dust exterior's source, past marginal tube, and excluded singular boundary."""
import json
import math
import sys
import unittest
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS/assets/data"
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))
import verify_metrics as vm


def load(kind):
    return json.loads((DATA / kind / "thakurta.json").read_text())


class DustExterior(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.entry = next(c for c in load("metrics")["coordinates"] if c["id"] == "conformal")
        cls.reader = vm.Reader(cls.entry["coords"], [p["symbol"] for p in cls.entry["parameters"]], ())
        cls.eta, cls.r = cls.reader.symbol[r"\eta"], cls.reader.symbol["r"]
        cls.dust = {cls.reader.parameters["a"]: cls.eta ** 2}
        cls.mass = {cls.reader.parameters["m"]: 1}

    def expression(self, value):
        return self.reader(value).subs(self.dust).doit().subs(self.mass)

    def components(self, field, variant=None):
        rows = self.entry[field] if variant is None else self.entry[field]["variants"][variant]["nonzero"]
        return {tuple(c["indices"]): self.expression(c["value"]) for c in rows}

    def number(self, expr, eta, r):
        return float(sp.sympify(expr).subs({self.eta: eta, self.r: r, self.reader.symbol[r"\theta"]: math.pi / 2}))

    def test_the_marginal_tube_is_past_and_the_near_boundary_region_is_anti_trapped(self):
        g = self.components("metric_components")
        area = g[(r"\theta", r"\theta")]
        def expansions(eta, r):
            radius = math.sqrt(self.number(area, eta, r))
            d0 = self.number(sp.diff(area, self.eta), eta, r) / (2 * radius)
            d1 = self.number(sp.diff(area, self.r), eta, r) / (2 * radius)
            k0 = 1 / math.sqrt(-2 * self.number(g[(r"\eta", r"\eta")], eta, r))
            k1 = 1 / math.sqrt(2 * self.number(g[("r", "r")], eta, r))
            return 2 * (k0 * d0 + k1 * d1) / radius, 2 * (k0 * d0 - k1 * d1) / radius
        for r in (2.5, 3, 4, 7, 12):
            outgoing, ingoing = expansions(2 * r * r / (r - 2), r)
            self.assertGreater(outgoing, 0)
            self.assertAlmostEqual(ingoing, 0, delta=1e-14)
        self.assertTrue(all(e > 0 for e in expansions(25, 2.05)))

    def test_zero_pressure_retains_inward_heat_flow_and_a_late_type_iv_region(self):
        g, G = self.components("metric_components"), self.components("einstein_tensor", "ll")
        for eta, r in ((3, 5), (100, 2.5)):
            rho = self.number(G[(r"\eta", r"\eta")] / -g[(r"\eta", r"\eta")], eta, r)
            pressure = self.number(G.get(("r", "r"), 0) / g[("r", "r")], eta, r)
            heat = -self.number(G[(r"\eta", "r")] / self.eta ** 4, eta, r)
            self.assertAlmostEqual(pressure, 0, delta=1e-15)
            self.assertLess(heat, 0)
            if eta == 100:
                self.assertLess((rho + pressure) ** 2 - 4 * heat ** 2, 0)

    def test_one_spatial_sheet_ends_outside_the_singularity(self):
        view, = load("embedding")["views"]
        surface, = view["surfaces"]
        piece, = surface["pieces"]
        self.assertGreater(min(p[0] for p in piece["points"]), 2)
        self.assertNotIn("movie", view)
        zigs = [l for v in load("conformal")["views"] for l in v["layers"] if l["kind"] == "zig"]
        self.assertEqual(len(zigs), 2)
        self.assertFalse(any(l["class"] == "horizon" for v in load("conformal")["views"] for l in v["layers"]))


if __name__ == "__main__":
    unittest.main()
