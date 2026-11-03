"""York's scalar radiation cavity, checked against its mass, wall and drawings."""
import json
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
sys.path.insert(0, str(ROOT / "_tools" / "derivations"))


class YorkEquilibrium(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import sympy as sp
        import verify_metrics as vm
        cls.sp = sp
        cls.entry = json.loads((DATA / "metrics" / "york_equilibrium.json").read_text())
        cls.chart = cls.entry["coordinates"][0]
        cls.reader = vm.Reader(cls.chart["coords"], [p["symbol"] for p in cls.chart["parameters"]], (),
                               kept=vm.ORDERS[("york_equilibrium", "areal")],
                               held=vm.HELD[("york_equilibrium", "areal")])

    def test_the_scalar_mass_and_energy_change_sign(self):
        sp, reader = self.sp, self.reader
        r, p = reader.symbol["r"], reader.parameters
        mu = reader.held[p["mu"]].subs(p["M"], 1)
        # The scalar contribution to mass is negative near the horizon and
        # positive farther out. Its slope has the sign of the energy density.
        near, far = sp.Rational(41, 20), sp.Integer(10)
        self.assertLess(float(mu.subs(r, near).evalf(30)), 0)
        self.assertGreater(float(mu.subs(r, far).evalf(30)), 0)
        self.assertLess(float(sp.diff(mu, r).subs(r, near).evalf(30)), 0)
        self.assertGreater(float(sp.diff(mu, r).subs(r, far).evalf(30)), 0)

    def test_matching_and_vacuum_limit(self):
        sp, reader = self.sp, self.reader
        p, r = reader.parameters, reader.symbol["r"]
        comps = {tuple(c["indices"]): reader(c["value"]).subs(reader.held).doit()
                 for c in self.chart["metric_components"]}
        M, R, eps = p["M"], p["R"], p["epsilon"]
        for radius in (3, 8, 40):
            at = {M: 1, R: 50, r: radius, eps: 0}
            self.assertEqual(sp.simplify(comps[("t", "t")].subs(at) + 1 - 2/sp.Integer(radius)), 0)
            self.assertEqual(sp.simplify(comps[("r", "r")].subs(at) - 1/(1 - 2/sp.Integer(radius))), 0)
        mu_wall = reader.held[p["mu"]].subs(r, R)
        expected = -(1 - 2*M*(1 + eps*mu_wall)/R)
        at = {M: 1, R: 50, eps: sp.Rational(1, 10000)}
        error = (comps[("t", "t")].subs(r, R) - expected).subs(at).evalf(40)
        self.assertLess(abs(float(error)), 1e-25)

    def test_drawings_use_the_cavity_and_have_bounded_shapes(self):
        embed = json.loads((DATA / "embedding" / "york_equilibrium.json").read_text())
        diagram = json.loads((DATA / "diagrams" / "york_equilibrium.json").read_text())
        view = diagram["systems"]["areal"][0]
        box = view["box"]
        self.assertTrue(0.5 <= (box[1] - box[0]) / (box[3] - box[2]) <= 2)
        piece = embed["views"][0]["surfaces"][0]["pieces"][0]
        self.assertAlmostEqual(piece["points"][0][0], 2)
        self.assertAlmostEqual(piece["points"][-1][0], 50)
        self.assertTrue(math.isfinite(piece["points"][-1][2]))


if __name__ == "__main__":
    unittest.main()
