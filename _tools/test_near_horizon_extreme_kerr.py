"""The throat of the extreme Kerr black hole, from the numbers in its files alone.

The readers of the published files and MFS_SPACETIME are those of test_build_mfs_data.py."""
import math
import unittest

from test_build_mfs_data import bisect, checks, diagram_file, embedding_file


@checks("near_horizon_extreme_kerr")
class ExtremeKerrThroat(unittest.TestCase):
    """The throat of the extreme Kerr black hole, from the numbers in its files alone, at r_0 = 1.
    On the equator of Bardeen and Horowitz's global chart the metric is

        ds^2 = (1/2)(-(1 + y^2) dtau^2 + dy^2/(1 + y^2)) + 2 (dphi + y dtau)^2,

    and the sphere of theta and phi at any event has g_thetatheta = (1 + cos^2 theta)/2 and
    g_phiphi = 2 sin^2 theta/(1 + cos^2 theta), the horizon of the extreme hole."""

    NAME = "near_horizon_extreme_kerr"

    def setUp(self):
        self.embedding = {v["id"]: v for v in embedding_file(self.NAME)["views"]}
        self.diagrams = diagram_file(self.NAME)

    @staticmethod
    def meridian(a, b, n=400):
        """The proper distance along a meridian of the horizon from theta = a to b, by Simpson's rule."""
        f = [math.sqrt((1 + math.cos(a + (b - a) * k / n) ** 2) / 2) for k in range(n + 1)]
        return (b - a) / n / 3 * (f[0] + f[-1] + 4 * sum(f[1:-1:2]) + 2 * sum(f[2:-1:2]))

    def test_the_equator_at_one_moment_is_a_cylinder_of_radius_two_GM(self):
        (piece,) = self.embedding["throat"]["surfaces"][0]["pieces"]
        self.assertEqual((piece["system"], piece["coordinate"]), ("global", "y"))
        for y, rho, z in piece["points"]:
            self.assertAlmostEqual(rho, math.sqrt(2), delta=1e-6)
            self.assertAlmostEqual(z, math.asinh(y) / math.sqrt(2), delta=1e-6)
        self.assertEqual((piece["points"][0][0], piece["points"][-1][0]), (-2.25, 2.25))

    def test_the_horizon_stands_in_flat_space_about_its_equator_and_in_minkowski_space_toward_its_poles(self):
        north, belt, south = self.embedding["horizon"]["surfaces"][0]["pieces"]
        self.assertEqual([p.get("space") for p in (north, belt, south)], ["minkowski", None, "minkowski"])
        # The surface lies level where cos^3 + cos^2 + 3 cos = 1, at 72.81 degrees and its mirror image.
        level = math.acos(bisect(lambda c: c ** 3 + c * c + 3 * c - 1, 0, 1))
        self.assertAlmostEqual(math.degrees(level), 72.8066, delta=1e-4)
        self.assertAlmostEqual(north["points"][-1][0], level, delta=1e-12)
        self.assertAlmostEqual(belt["points"][0][0], level, delta=1e-12)
        self.assertAlmostEqual(belt["points"][-1][0], math.pi - level, delta=1e-12)
        self.assertAlmostEqual(south["points"][0][0], math.pi - level, delta=1e-12)
        for piece in (north, belt, south):
            sign = -1 if piece.get("space") == "minkowski" else 1
            points = piece["points"]
            for theta, rho, _ in points:
                self.assertAlmostEqual(rho, math.sqrt(2) * math.sin(theta) / math.sqrt(1 + math.cos(theta) ** 2),
                                       delta=1e-6)
            # Every chord against the metric distance along the meridian between its ends: in Minkowski
            # space a chord is sqrt(drho^2 - dZ^2), and near the pole of a cap, where the surface bends
            # most, it falls short of its arc by under a part in a thousand.
            for (a, ra, za), (b, rb, zb) in zip(points, points[1:]):
                chord = math.sqrt((rb - ra) ** 2 + sign * (zb - za) ** 2)
                self.assertAlmostEqual(chord / self.meridian(a, b, 20), 1, delta=1e-3, msg=f"{piece['id']} at {a}")
        # The two caps and the belt meet at one height, and the whole meridian is 2.70 r_0 long.
        self.assertAlmostEqual(north["points"][-1][2], belt["points"][0][2], delta=1e-7)
        self.assertAlmostEqual(belt["points"][-1][2], south["points"][0][2], delta=1e-7)
        self.assertAlmostEqual(self.meridian(0, math.pi), 2.7013, delta=1e-4)

    def test_the_rays_of_no_angular_momentum_cross_the_throat_in_pi(self):
        (view,) = self.diagrams["systems"]["global"]
        X0, X1, Y0, Y1 = view["box"]
        self.assertEqual(view["families"], ["moving left", "moving right"])
        for family, sign in (("P", 1), ("M", -1)):
            self.assertGreater(len(view["rays"][family]), 10)
            for ray in view["rays"][family]:
                kept = [u[1] * (Y1 - Y0) + Y0 + sign * math.atan(u[0] * (X1 - X0) + X0) for u in ray]
                self.assertLess(max(kept) - min(kept), 2e-3, family)
        # d/dtau is null on the equator where 3 y^2 = 1, which the plane marks on both sides.
        (marker,) = [m for m in view["markers"] if m["kind"] == "gtt"]
        marked = sorted(line[0][0] * (X1 - X0) + X0 for line in marker["lines"])
        for at, y in zip(marked, (-1 / math.sqrt(3), 1 / math.sqrt(3))):
            self.assertAlmostEqual(at, y, delta=1e-3)

    def test_the_cones_tip_over_beyond_the_circles_where_d_tau_is_null(self):
        (figure,) = self.diagrams["projections"]["global"]
        edge = 1 / math.sqrt(3)
        seen = set()
        for cone in figure["turn"]["cones"]:
            ax, ay, at = cone["apex"]
            radius = math.hypot(ax, ay)
            # A circle of constant y stands at the radius 2 + arsinh(y)/sqrt 2.
            y = math.sinh(math.sqrt(2) * (radius - 2))
            self.assertEqual(at, 0)
            turn = [(-ay * (x - ax) + ax * (v - ay)) / radius for x, v, _ in cone["rim"]]
            nearest = min((-1.5, -edge, 0.0, edge, 1.5), key=lambda c: abs(c - y))
            self.assertAlmostEqual(y, nearest, delta=1e-5)
            seen.add(nearest)
            if nearest == 0:
                self.assertLess(min(turn), -0.1)
                self.assertGreater(max(turn), 0.1)
                self.assertAlmostEqual(min(turn), -max(turn), delta=1e-6)
            elif abs(nearest) == 1.5:
                # Every future null direction turns toward -phi where y is positive and toward +phi
                # where it is negative.
                self.assertTrue(all(-math.copysign(1, y) * t > 0.02 for t in turn), y)
            else:
                # One edge of the cone stands vertical, and no generator turns back past it: of the
                # 96 drawn, the one nearest that edge turns by under a thousandth of its length.
                furthest = max(turn) if y > 0 else -min(turn)
                self.assertLessEqual(furthest, 1e-9)
                self.assertGreater(furthest, -3e-4)
        self.assertEqual(len(seen), 5)
        self.assertEqual(len(figure["turn"]["cones"]), 20)


if __name__ == "__main__":
    unittest.main()
