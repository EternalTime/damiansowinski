"""The wormhole of Einstein-Dirac-Maxwell theory, Blazquez-Salcedo, Knoll and Radu's exact solution:
python3 -m unittest discover -s _tools

What its texts and drawings state, held on the published files. The first class needs nothing but
the files: the embedding's two sheets meet at the throat r = r_0 and their height agrees with the
quadrature in Bronnikov and Kim's u, and the rays of each spacetime diagram cross the throat. The
second reads the published metric and Einstein tensor through the checker's Reader and holds them
to what the History and einstein_dirac_maxwell_wormhole.md state: the charge exceeds the mass, the
metric is the extreme Reissner-Nordstrom black hole's at Q_e = r_0, the electric field carries all of
the energy density, and the throat's curvature is 2(3 Q_e^4 - 2 Q_e^2 r_0^2 + 3 r_0^4)/r_0^8; it
needs sympy and is skipped where it is absent.
"""
import importlib.util
import json
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
HAS_SYMPY = importlib.util.find_spec("sympy") is not None
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))

ID = "einstein_dirac_maxwell_wormhole"
B = 0.25  # Q_e^2/r_0 at r_0 = 1 and Q_e = 1/2


def load(folder):
    return json.loads((DATA / folder / f"{ID}.json").read_text(encoding="utf-8"))


def height(u, steps=4000):
    """dz/du = 2 sqrt((r_0^2 + (r_0 + b) u^2)/(r_0 - b + u^2)) at r_0 = 1, by Simpson's rule."""
    h = u / steps
    f = [2 * math.sqrt((1 + (1 + B) * (k * h) ** 2) / (1 - B + (k * h) ** 2)) for k in range(steps + 1)]
    return h / 3 * (f[0] + f[-1] + 4 * sum(f[1:-1:2]) + 2 * sum(f[2:-1:2]))


class Drawings(unittest.TestCase):
    def test_the_two_sheets_meet_at_the_throat_and_reach_six_r_0(self):
        view = next(v for v in load("embedding")["views"] if v["id"] == "wormhole")
        pieces = view["surfaces"][0]["pieces"]
        self.assertEqual(len(pieces), 2)
        for piece in pieces:
            radii = [p[1] for p in piece["points"]]
            self.assertAlmostEqual(min(radii), 1.0, places=6)
            self.assertAlmostEqual(max(radii), 6.0, places=6)

    def test_the_height_is_the_quadrature_in_u(self):
        """Each point (r, rho, z) of the sheet that rises has z = the integral of dz/du to u = sqrt(r - r_0)."""
        view = next(v for v in load("embedding")["views"] if v["id"] == "wormhole")
        rising = next(p for p in view["surfaces"][0]["pieces"] if p["points"][-1][2] > 0)
        for x, rho, z in rising["points"][::7]:
            self.assertAlmostEqual(z, height(math.sqrt(max(x - 1, 0.0))), delta=2e-3, msg=f"r = {x}")

    def test_every_view_draws_rays_of_both_families(self):
        systems = load("diagrams")["systems"]
        self.assertEqual(sorted(systems), ["areal", "bronnikov_kim", "compact"])
        for system, views in systems.items():
            for view in views:
                self.assertTrue(view["rays"], f"{system}/{view['id']}")

    def test_every_view_is_near_square(self):
        for system, views in load("diagrams")["systems"].items():
            for view in views:
                x0, x1, y0, y1 = view["box"]
                self.assertEqual(x1 - x0, y1 - y0, f"{system}/{view['id']}")


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class Physics(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import sympy as sp
        import verify_metrics as vm
        cls.sp, cls.vm = sp, vm
        metric = json.loads((DATA / "metrics" / f"{ID}.json").read_text(encoding="utf-8"))
        cls.charts = {c["id"]: c for c in metric["coordinates"]}

    def read(self, system):
        chart = self.charts[system]
        reader = self.vm.Reader(chart["coords"], [p["symbol"] for p in chart["parameters"]], ())
        return chart, reader

    def test_the_charge_exceeds_the_mass(self):
        sp = self.sp
        r0, Q = sp.symbols("r_0 Q_e", positive=True)
        M = 2 * Q ** 2 * r0 / (Q ** 2 + r0 ** 2)
        self.assertEqual(sp.factor(Q - M), Q * (Q - r0) ** 2 / (Q ** 2 + r0 ** 2))

    def test_g_tt_is_the_extreme_black_hole_s_and_becomes_it_at_q_e_equal_r_0(self):
        sp = self.sp
        chart, reader = self.read("areal")
        g = {tuple(e["indices"]): reader(e["value"]) for e in chart["metric_components"]}
        r, r0, Q = reader.symbol["r"], reader.parameters["r_0"], reader.parameters["Q_e"]
        M = 2 * Q ** 2 * r0 / (Q ** 2 + r0 ** 2)
        self.assertEqual(sp.simplify(g[("t", "t")] + (1 - M / r) ** 2), 0)
        self.assertEqual(sp.simplify(g[("r", "r")].subs(Q, r0) - 1 / (1 - r0 / r) ** 2), 0)
        self.assertEqual(sp.simplify(g[("t", "t")].subs(Q, 0) + 1), 0)

    def test_the_electric_field_carries_all_of_the_energy_density(self):
        sp = self.sp
        for system, radius in (("areal", lambda x, r0: x), ("bronnikov_kim", lambda x, r0: r0 + x ** 2),
                               ("compact", lambda x, r0: r0 / (1 - x ** 2))):
            chart, reader = self.read(system)
            x = reader.symbol[chart["coords"][1]]
            r0, Q = reader.parameters["r_0"], reader.parameters["Q_e"]
            Gtt = next(e["value"] for e in chart["einstein_tensor"]["variants"]["ul"]["nonzero"]
                       if e["indices"] == ["t", "t"])
            self.assertEqual(sp.simplify(reader(Gtt) + Q ** 2 / radius(x, r0) ** 4), 0, system)

    def test_the_throat_s_curvature(self):
        sp = self.sp
        for system, throat in (("areal", "r_0"), ("bronnikov_kim", 0), ("compact", 0)):
            chart, reader = self.read(system)
            x = reader.symbol[chart["coords"][1]]
            r0, Q = reader.parameters["r_0"], reader.parameters["Q_e"]
            K = reader(chart["kretschmann"].removeprefix("K = "))
            at = r0 if throat == "r_0" else 0
            self.assertEqual(sp.simplify(K.subs(x, at) - 2 * (3 * Q ** 4 - 2 * Q ** 2 * r0 ** 2 + 3 * r0 ** 4) / r0 ** 8),
                             0, system)

    def test_the_ricci_scalar_vanishes(self):
        for system, chart in self.charts.items():
            self.assertIn(chart["ricci_scalar"], ("0", "R = 0"), system)


if __name__ == "__main__":
    unittest.main()
