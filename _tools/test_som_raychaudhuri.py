"""Som and Raychaudhuri's rotating dust, from the numbers in its files alone.

The readers of the published files and MFS_SPACETIME are those of test_build_mfs_data.py."""
import json
import math
import unittest

import build_mfs_data as build
from test_build_mfs_data import checks


@checks("som_raychaudhuri")
class SomRaychaudhuri(unittest.TestCase):
    """Som and Raychaudhuri's published diagrams against their metric, from the numbers in the
    files, in units of r_c = c/Omega: on the cylinder of t and phi at radius r the null curves are
    c dt = r (1 - r) dphi and c dt = -r (1 + r) dphi, and the circle of constant t and r is null at
    r = 1."""

    @classmethod
    def setUpClass(cls):
        data = json.loads((build.DIAGRAMS_DIR / "som_raychaudhuri.json").read_text(encoding="utf-8"))
        cls.views = {(system, v["id"]): v for system, views in data["systems"].items() for v in views}
        cls.figure = data["projections"]["cylindrical"][0]

    def slopes(self, view):
        """d(ct)/d(r phi) of every ray from end to end, by family; a ray is straight, and its points
        are rounded to a ten thousandth of the drawing, so a ray shorter than a tenth is passed over."""
        X0, X1, Y0, Y1 = view["box"]
        found = {}
        for family, lines in view["rays"].items():
            for line in lines:
                (ua, wa), (ub, wb) = line[0], line[-1]
                if abs(ub - ua) > 0.1:
                    found.setdefault(family, []).append((wb - wa) * (Y1 - Y0) / ((ub - ua) * (X1 - X0)))
        return found

    def test_the_rays_of_each_cylinder_are_the_null_lines_of_the_metric(self):
        for view, r in (("inside", 0.5), ("beyond", 1.5)):
            found = self.slopes(self.views["cylindrical", view])
            self.assertEqual(len(found), 2, view)
            want = {1 - r, -(1 + r)}
            for family, slopes in found.items():
                nearest = min(want, key=lambda k: abs(k - slopes[0]))
                want.discard(nearest)
                for k in slopes:
                    self.assertAlmostEqual(k, nearest, delta=2e-3, msg=f"{view} {family}")
                # -(c dt + r^2 dphi)^2 + r^2 dphi^2 = 0 with c dt = k r dphi.
                self.assertAlmostEqual(-(nearest * r + r * r) ** 2 + r * r, 0.0, places=12)
            self.assertEqual(want, set(), view)

    def test_beyond_the_null_circle_both_families_run_down_in_t_toward_plus_phi(self):
        inside = self.slopes(self.views["cylindrical", "inside"])
        beyond = self.slopes(self.views["cylindrical", "beyond"])
        self.assertEqual(sorted(k[0] > 0 for k in inside.values()), [False, True])
        self.assertTrue(all(k[0] < 0 for k in beyond.values()))

    def test_the_rays_of_the_plane_of_t_and_x_run_at_45_degrees(self):
        for slopes in self.slopes(self.views["cartesian", "tx"]).values():
            for k in slopes:
                self.assertAlmostEqual(abs(k), 1.0, delta=2e-3)

    def test_the_figure_draws_the_null_circle_at_r_c_and_every_cone_null(self):
        turn = self.figure["turn"]
        critical = [line["points"] for line in turn["lines"] if line["class"] == "critical"]
        self.assertEqual(len(critical), 1)
        for X, Y, T in critical[0]:
            self.assertAlmostEqual(math.hypot(X, Y), 1.0, delta=1e-5)
            self.assertAlmostEqual(T, 0.0, delta=1e-9)
        # A generator (dX, dY, dT) at (X, Y) is null: with r^2 dphi = X dY - Y dX,
        # -(dT + X dY - Y dX)^2 + dX^2 + dY^2 = 0.
        self.assertEqual(len(turn["cones"]), 13)
        for cone in turn["cones"]:
            X, Y, T = cone["apex"]
            for x, y, t in cone["rim"]:
                dX, dY, dT = x - X, y - Y, t - T
                size = dX * dX + dY * dY + dT * dT
                self.assertAlmostEqual((-(dT + X * dY - Y * dX) ** 2 + dX * dX + dY * dY) / size, 0.0, delta=2e-4)


if __name__ == "__main__":
    unittest.main()
