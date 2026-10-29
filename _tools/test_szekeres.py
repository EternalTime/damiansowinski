"""Szekeres's cloud: the marked ray on each half of its axis, from the numbers in the file.

The readers of the published files and MFS_SPACETIME are those of test_build_mfs_data.py."""
import math
import unittest

from test_build_mfs_data import checks, diagram_file


@checks("szekeres")
class SzekeresAxis(unittest.TestCase):
    """The marked ray on each half of the axis of Szekeres's cloud, from the numbers in the file:
    inside the cloud it obeys c dt = (dR/dr +- R S'/S) dr, the upper sign toward theta = 0, with
    the declared R^(3/2) = r^(3/2) - (3/2) sqrt(2M) ct, 2M = r^3(5 - 3r^2)/4 and S'/S = 2r(1 - r^2),
    and outside it stays on R = r_s = 1/2, the horizon of the exterior."""

    @staticmethod
    def slope(r, t, sign):
        root = math.sqrt(5 - 3 * r * r)
        u = 1 - 3 * root * t / 4
        R, dR = r * u ** (2 / 3), u ** (2 / 3) + 3 * r * r * t / (2 * u ** (1 / 3) * root)
        return dR + sign * R * 2 * r * (1 - r * r)

    def ray(self, sign, r_to, steps=2000):
        """ct at r_to along the ray that crosses the surface r = 1 where R = 1/2, by Runge and Kutta."""
        r, t = 1.0, (1 - 0.5 ** 1.5) / (1.5 * math.sqrt(0.5))
        h = (r_to - r) / steps
        for _ in range(steps):
            k1 = self.slope(r, t, sign)
            k2 = self.slope(r + h / 2, t + h * k1 / 2, sign)
            k3 = self.slope(r + h / 2, t + h * k2 / 2, sign)
            k4 = self.slope(r + h, t + h * k3, sign)
            r, t = r + h, t + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6
        return t

    def test_the_marked_ray_is_null_in_the_cloud_and_the_horizon_outside(self):
        views = {view["id"]: view for view in diagram_file("szekeres")["systems"]["axisymmetric"]}
        self.assertEqual(set(views), {"north", "south"})
        for half, sign in (("north", 1), ("south", -1)):
            view = views[half]
            x0, x1, y0, y1 = view["box"]
            line, = next(m for m in view["markers"] if m["kind"] == "event")["lines"]
            points = [(x0 + (x1 - x0) * x, y0 + (y1 - y0) * y) for x, y in line]
            inside = [(r, t) for r, t in points if r < 1]
            outside = [(r, t) for r, t in points if r >= 1]
            self.assertGreater(len(inside), 3, half)
            self.assertGreater(len(outside), 1, half)
            for r, t in inside:
                self.assertAlmostEqual(t, self.ray(sign, r), delta=3e-4, msg=f"{half}, r = {r}")
            for r, t in outside:
                self.assertAlmostEqual((r ** 1.5 - 1.5 * math.sqrt(0.5) * t) ** (2 / 3), 0.5, delta=3e-4, msg=f"{half}, r = {r}")
            # The times the caption states: the ray leaves the centre at -0.81 toward theta = 0 and at -0.12 toward pi.
            self.assertEqual(f"{self.ray(sign, 1e-9):.2f}", {"north": "-0.81", "south": "-0.12"}[half])
            self.assertIn({"north": "$ct = -0.81\\,r_b$", "south": "$ct = -0.12\\,r_b$"}[half], " ".join(view["caption"]))


if __name__ == "__main__":
    unittest.main()
