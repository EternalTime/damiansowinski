"""York's scalar radiation cavity, checked against its wall and Page's stress."""
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

    def test_horizon_constants_and_page_energy_density(self):
        sp, reader = self.sp, self.reader
        r, p = reader.symbol["r"], reader.parameters
        M, epsilon = p["M"], p["epsilon"]
        mu, rho = (reader.held[p[n]] for n in ("mu", "rho"))
        self.assertEqual(sp.simplify(mu.subs(r, 2*M)), 0)
        self.assertEqual(sp.simplify(rho.subs(r, 2*M)), 0)
        # Taylor and Breen (2021), Eq. 42 at zero cosmological constant.
        f = 1 - 2*M/r
        kappa = 1/(4*M)
        energy = epsilon*M**2*(66*M**4 - 72*M**3*r + 20*M**2*r**2 - 2*kappa**4*r**8) / (960*sp.pi**2*r**8*f**2)
        residual = sp.diff(mu, r) + 4*sp.pi*r**2*energy/(epsilon*M)
        for radius in (3, 7, 20, 40):
            value = residual.subs({M: 1, r: radius, epsilon: sp.Rational(1, 10000)}).evalf(40)
            self.assertLess(abs(float(value)), 1e-25)

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
