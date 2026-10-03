"""Interstellar's wormhole, the Dneg wormhole of James, von Tunzelmann, Franklin and Thorne:
python3 -m unittest discover -s _tools

What its texts and drawings state, held on the published files. The first class needs nothing but
the files: the embedding is a vertical tube of radius rho between the mouths and flares out beyond
them, turning from vertical to 45 degrees across the lensing width W = 1.42953 M, and every plane of
t and l is flat, its rays at 45 degrees. The second reads the published metric through the checker's
Reader and holds it to what the History states: on the cylinder the wormhole is the product of a flat
plane and a sphere of radius rho, Plebanski and Hacyan's, and the flare meets it at the mouth with
r = rho and dr/dl = 0, climbs as (2/pi) arctan x, and tends to dr/dl = 1 far away; it needs sympy and
is skipped where it is absent.
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

ID = "interstellar_wormhole"
A, M = 1.0, 0.5  # the drawings' a and M in units of rho, the wormhole of the paper's figure 3


def load(folder):
    return json.loads((DATA / folder / f"{ID}.json").read_text(encoding="utf-8"))


def radius(ell):
    x = 2 * (abs(ell) - A) / (math.pi * M)
    return 1.0 if abs(ell) <= A else 1 + M * (x * math.atan(x) - 0.5 * math.log1p(x * x))


class Drawings(unittest.TestCase):
    def setUp(self):
        view = next(v for v in load("embedding")["views"] if v["id"] == "wormhole")
        (self.piece,) = view["surfaces"][0]["pieces"]
        self.points = self.piece["points"]

    def test_the_embedding_is_a_tube_of_radius_rho_between_the_mouths(self):
        tube = [(x, rho, z) for x, rho, z in self.points if abs(x) <= A]
        self.assertGreater(len(tube), 2)
        for x, rho, z in tube:
            self.assertAlmostEqual(rho, 1.0, places=5, msg=f"l = {x}")
            self.assertAlmostEqual(z, x, places=5, msg=f"l = {x}")

    def test_every_circle_has_the_radius_r_of_l(self):
        for x, rho, _ in self.points:
            self.assertAlmostEqual(rho, radius(x), places=5, msg=f"l = {x}")

    def test_the_lensing_width_is_1_42953_M(self):
        """Where the profile turns through 45 degrees, rho - 1 is W = 1.42953 M, read off the published
        points: the slope dz/drho of each chord, at its midpoint, interpolated to where it is 1."""
        outside = [(rho, z) for x, rho, z in self.points if x >= A]
        mids = [((r0 + r1) / 2, (z1 - z0) / (r1 - r0)) for (r0, z0), (r1, z1) in zip(outside, outside[1:]) if r1 > r0]
        for (ra, sa), (rb, sb) in zip(mids, mids[1:]):
            if sa >= 1 > sb:
                width = ra + (rb - ra) * (sa - 1) / (sa - sb) - 1
                self.assertAlmostEqual(width, 1.42953 * M, delta=0.003)
                return
        self.fail("the profile never turns through 45 degrees")

    def test_every_plane_of_t_and_l_is_flat_and_drawn_square(self):
        systems = load("diagrams")["systems"]
        self.assertEqual(sorted(systems), ["cylinder", "flare", "proper_distance"])
        for system, views in systems.items():
            for view in views:
                self.assertTrue(view["rays"], f"{system}/{view['id']}")
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

    def reader(self, chart):
        c = self.charts[chart]
        return self.vm.Reader(c["coords"], [p["symbol"] for p in c["parameters"]])

    def metric(self, chart, i):
        c = self.charts[chart]
        values = {tuple(e["indices"]): e["value"] for e in c["metric_components"]}
        return self.reader(chart)(values[(c["coords"][i], c["coords"][i])])

    def test_the_cylinder_is_a_flat_plane_times_a_sphere_of_radius_rho(self):
        """Plebanski and Hacyan's flat plane times a sphere: the same line element with b = rho."""
        c = self.charts["cylinder"]
        R = self.reader("cylinder")
        rho = R.parameters["rho"]
        theta = R.symbol["\\theta"]
        expected = [-1, 1, rho ** 2, rho ** 2 * self.sp.sin(theta) ** 2]
        for i, want in enumerate(expected):
            self.assertEqual(self.sp.simplify(self.metric("cylinder", i) - want), 0, c["coords"][i])
        self.assertEqual(c["kretschmann"], "K = \\dfrac{4}{\\rho^4}")

    def test_the_flare_meets_the_cylinder_smoothly_and_opens_into_flat_space(self):
        sp = self.sp
        R = self.reader("flare")
        ell, rho, a, m = R.symbol["\\ell"], R.parameters["rho"], R.parameters["a"], R.parameters["M"]
        r = R.parameters["r"]
        self.assertEqual(sp.simplify(r.subs(ell, a) - rho), 0)
        slope = sp.diff(r, ell)
        self.assertEqual(sp.simplify(slope.subs(ell, a)), 0)
        x = 2 * (ell - a) / (sp.pi * m)
        self.assertEqual(sp.simplify(slope - 2 * sp.atan(x) / sp.pi), 0)
        positive = {m: sp.Symbol("m", positive=True), a: sp.Symbol("a", positive=True)}
        self.assertEqual(sp.limit(slope.subs(positive), ell, sp.oo), 1)
        g = self.metric("flare", 2)
        self.assertEqual(sp.simplify(g - r ** 2), 0)

    def test_the_clocks_keep_the_same_time_in_every_chart(self):
        for chart in self.charts:
            self.assertEqual(self.metric(chart, 0), -1, chart)


if __name__ == "__main__":
    unittest.main()
