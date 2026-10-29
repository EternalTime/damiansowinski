"""Kastor and Traschen's black holes in de Sitter space, from the numbers in their files alone.

The readers of the published files and MFS_SPACETIME are those of test_build_mfs_data.py."""
import json
import math
import unittest

import build_mfs_data as build
from test_build_mfs_data import checks


@checks("kastor_traschen")
class KastorTraschen(unittest.TestCase):
    """Kastor and Traschen's published diagrams against their metric, from the numbers in the files:
    two holes of m = 1 at z = +-2 with H = -3/32, where U = H tau + 2/sqrt(x^2 + 4) on the midplane,
    and one hole with H = -3/16 at H tau = 1/2."""

    H = -3 / 32

    @classmethod
    def setUpClass(cls):
        cls.diagrams = json.loads((build.DIAGRAMS_DIR / "kastor_traschen.json").read_text(encoding="utf-8"))["systems"]
        cls.views = {v["id"]: v for v in json.loads(
            (build.EMBEDDING_DIR / "kastor_traschen.json").read_text(encoding="utf-8"))["views"]}

    def midplane(self):
        view = next(v for v in self.diagrams["cartesian"] if v["id"] == "tx")
        X0, X1, Y0, Y1 = view["box"]
        chart = lambda line: [(X0 + u * (X1 - X0), Y0 + w * (Y1 - Y0)) for u, w in line]
        return view, chart

    def U(self, tau, x):
        return self.H * tau + 2 / math.sqrt(x * x + 4)

    def test_the_marked_rays_of_the_midplane_are_null_and_leave_the_axis_when_the_horizons_join(self):
        view, chart = self.midplane()
        lines = [chart(line) for m in view["markers"] if m["kind"] == "event" for line in m["lines"]]
        self.assertEqual(len(lines), 2)
        for line in lines:
            x, tau = line[0]
            self.assertAlmostEqual(x, 0.0, delta=5e-3)
            self.assertAlmostEqual(tau, -6.1995, delta=5e-3)
            # d(c tau)/dx = +-U^2 along the ray, by the midpoint rule over each chord.
            for (xa, ta), (xb, tb) in zip(line, line[1:]):
                want = self.U((ta + tb) / 2, (xa + xb) / 2) ** 2 * abs(xb - xa)
                self.assertAlmostEqual(tb - ta, want, delta=0.02 * want + 5e-3)
            # Far out the ray keeps H tau x = 2/3 more and more nearly, the horizon of the merged hole.
            x, tau = line[-1]
            self.assertAlmostEqual(self.H * tau * abs(x), 2 / 3, delta=0.02)

    def test_the_singular_curve_of_the_midplane_is_where_U_vanishes_and_no_ray_passes_it(self):
        view, chart = self.midplane()
        curve = [p for m in view["markers"] if m["kind"] == "singular" for line in m["lines"] for p in chart(line)]
        self.assertGreater(len(curve), 20)
        for x, tau in curve:
            self.assertAlmostEqual(self.U(tau, x), 0.0, delta=2e-3)
        for family in "PM":
            for line in view["rays"][family]:
                self.assertTrue(all(self.U(tau, x) > -2e-3 for x, tau in chart(line)))

    def test_the_midplane_of_two_holes_is_embedded_with_the_circumference_radius_rho_U(self):
        frames = self.views["two_holes"]["movie"]["frames"]
        self.assertAlmostEqual(frames[0]["value"], -8.0)
        self.assertAlmostEqual(frames[-1]["value"], -1.0)
        for frame in frames:
            tau, points = frame["value"], frame["pieces"][0]["points"]
            for (xa, ra, za), (xb, rb, zb) in zip(points, points[1:]):
                self.assertAlmostEqual(rb, xb * self.U(tau, xb), delta=1e-5)
                # A chord of the profile is as long as the metric distance U d rho between its circles.
                proper = self.U(tau, (xa + xb) / 2) * (xb - xa)
                self.assertAlmostEqual(math.hypot(rb - ra, zb - za), proper, delta=2e-3 * proper + 1e-6)
            # Every circle shrinks as the universe contracts.
        widest = [frame["pieces"][0]["points"][-1][1] for frame in frames]
        self.assertEqual(widest, sorted(widest, reverse=True))

    def test_one_hole_is_embedded_as_the_extremal_throat_over_its_areal_radius(self):
        points = self.views["one_hole"]["surfaces"][0]["pieces"][0]["points"]

        def throat(r):
            w = math.sqrt(r + 1)
            return 2 * w + math.log((w - 1) / (w + 1))
        for r, rho, z in points:
            self.assertAlmostEqual(rho, r / 2 + 1, delta=1e-6)
            self.assertAlmostEqual(z, throat(r) - throat(points[0][0]), delta=2e-4)


if __name__ == "__main__":
    unittest.main()
