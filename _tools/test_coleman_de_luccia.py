"""Coleman and De Luccia's bubble, from the numbers in its files alone.

The readers of the published files and MFS_SPACETIME are those of test_build_mfs_data.py."""
import math
import unittest
from fractions import Fraction

from test_build_mfs_data import checks, conformal_file, diagram_file, embedding_file


@checks("coleman_de_luccia")
class ColemanDeLucciaBubble(unittest.TestCase):
    """Coleman and De Luccia's bubble, from the numbers in its files alone, in units of the curvature
    radius l. With their rho_0 = l a decay into flat space has its wall at rho_bar = 4/5 and a decay
    of flat space at 4/3, where rho' = sqrt(1 - Lambda rho^2/3) drops by rho_0 rho_bar/(2 l^2), which
    is Israel's condition for their S_1 = epsilon rho_0/3."""

    NAME = "coleman_de_luccia"

    def setUp(self):
        self.embedding = {v["id"]: v for v in embedding_file(self.NAME)["views"]}
        self.diagrams = diagram_file(self.NAME)["systems"]
        self.conformal = {v["id"]: v for v in conformal_file(self.NAME)["views"]}

    def view(self, system, view):
        return next(v for v in self.diagrams[system] if v["id"] == view)

    def test_the_junction_condition_gives_coleman_and_de_luccias_bubble_radii(self):
        for rho_0 in (Fraction(1, 3), Fraction(1), Fraction(3, 2)):
            jump = rho_0 / 2
            into_flat = rho_0 / (1 + rho_0 ** 2 / 4)            # their (3.15)
            out_of_flat = rho_0 / (1 - rho_0 ** 2 / 4)          # their (3.18)
            # 1 - sqrt(1 - rho^2) = jump rho and sqrt(1 + rho^2) - 1 = jump rho, squared with the root alone.
            self.assertEqual((1 - jump * into_flat) ** 2, 1 - into_flat ** 2)
            self.assertEqual((1 + jump * out_of_flat) ** 2, 1 + out_of_flat ** 2)
        self.assertEqual(Fraction(1) / (1 + Fraction(1, 4)), Fraction(4, 5))
        self.assertEqual(Fraction(1) / (1 - Fraction(1, 4)), Fraction(4, 3))
        # No bubble of the decay into flat space is wider than the de Sitter radius, and the decay of flat
        # space has none once rho_0 reaches 2 l.
        self.assertLessEqual(max(r / (1 + r * r / 4) for r in (Fraction(k, 10) for k in range(1, 80))), 1)
        with self.assertRaises(ZeroDivisionError):
            Fraction(2) / (1 - Fraction(2) ** 2 / 4)

    def test_each_moment_inside_the_bubble_is_a_hyperboloid_of_radius_sin_tau(self):
        view = self.embedding["hyperboloids"]
        times = [s["time"] for s in view["surfaces"]]
        self.assertEqual(len(times), 5)
        for k, (surface, t) in enumerate(zip(view["surfaces"], times), start=1):
            self.assertAlmostEqual(t, math.pi * k / 6, delta=1e-6)
            (piece,) = surface["pieces"]
            self.assertEqual((piece["system"], piece["coordinate"], piece.get("space")), ("open", "\\chi", "minkowski"))
            a = math.sin(math.pi * k / 6)
            for chi, rho, z in piece["points"]:
                self.assertAlmostEqual(rho, a * math.sinh(chi), delta=1e-6)
                self.assertAlmostEqual(z, a * math.cosh(chi), delta=1e-6)
                self.assertAlmostEqual(z * z - rho * rho, a * a, delta=1e-5)
        # The universe is largest at the middle moment and the moments either side of it are alike.
        radii = [s["pieces"][0]["points"][0][2] for s in view["surfaces"]]
        self.assertAlmostEqual(radii[2], 1.0, delta=1e-9)
        self.assertAlmostEqual(radii[0], radii[4], delta=1e-9)
        self.assertAlmostEqual(radii[1], radii[3], delta=1e-9)

    def test_the_rays_of_the_wall_chart_keep_psi_plus_and_minus_the_conformal_distance(self):
        def distance(xi, into_flat):
            if into_flat:
                return math.log(xi) if xi < 0.8 else math.log(math.tan((xi - 0.8 + math.asin(0.8)) / 2)) + math.log(1.6)
            return (math.log(math.tanh(xi / 2)) if xi < math.log(3)
                    else math.log(xi - math.log(3) + 4 / 3) + math.log(3 / 8))
        for name, into_flat in (("into_flat", True), ("out_of_flat", False)):
            view = self.view("wall", name)
            X0, X1, Y0, Y1 = view["box"]
            crossed = 0
            for family, sign in (("P", 1), ("M", -1)):
                for ray in view["rays"][family]:
                    points = [(Y0 + v * (Y1 - Y0), X0 + u * (X1 - X0)) for u, v in ray]
                    points = [(psi, xi) for psi, xi in points if 0.2 < xi < 2.9]
                    kept = [psi + s * distance(xi, into_flat) for s in (1, -1) for psi, xi in points]
                    half = len(points)
                    spread = min(max(kept[:half]) - min(kept[:half]), max(kept[half:]) - min(kept[half:])) if half > 1 else 0
                    self.assertLess(spread, 4e-3, (name, family))
                    wall = 0.8 if into_flat else math.log(3)
                    crossed += half > 1 and min(xi for _, xi in points) < wall < max(xi for _, xi in points)
            self.assertGreater(crossed, 10)

    def test_the_static_charts_end_at_the_wall(self):
        """Each static view hatches what lies across the wall, and the edge of the hatching inside the
        box is r_w(t) to the grid the hatching is found on, 1/160 of the box."""
        walls = {("static_inside", "into_flat"): lambda t: math.sqrt(0.64 + t * t),
                 ("static_outside", "out_of_flat"): lambda t: math.sqrt(16 / 9 + t * t),
                 ("static_outside", "into_flat"): lambda t: math.sqrt(0.64 + 0.36 * math.tanh(t) ** 2),
                 ("static_inside", "out_of_flat"): lambda t: math.sqrt(16 / 9 + 25 / 9 * math.tan(t) ** 2)}
        for (system, name), wall in walls.items():
            view = self.view(system, name)
            X0, X1, Y0, Y1 = view["box"]
            (ring,) = view["hatch"]
            on_wall = 0
            for u, v in ring:
                if min(u, v, 1 - u, 1 - v) < 1e-3:
                    continue        # the box's own edge
                r, t = X0 + u * (X1 - X0), Y0 + v * (Y1 - Y0)
                if system == "static_inside" and name == "out_of_flat" and t > 0.9:
                    continue        # the wall runs steeply out there, on its way to infinity at pi/2
                self.assertAlmostEqual(r, wall(t), delta=1.5 * (X1 - X0) / 160 + 1.5 * (Y1 - Y0) / 160, msg=(system, name, t))
                on_wall += 1
            self.assertGreater(on_wall, 10, (system, name))
            # The hatched side is the far side of the wall: the corner r = 0, t = 0 for a chart that
            # begins at the wall, the corner of greatest r for one that ends there.
            corner = [0.0, 0.0] if system == "static_outside" else [1.0, 0.0]
            self.assertIn(corner, [[round(u, 6), round(v, 6)] for u, v in ring], (system, name))

    def test_the_wall_of_each_conformal_view_runs_from_rest_to_where_the_light_cone_meets_infinity(self):
        for name in ("into_flat", "out_of_flat"):
            view = self.conformal[name]
            (wall,) = [layer["points"] for layer in view["layers"] if layer["class"] == "surface"]
            self.assertAlmostEqual(wall[0][0], 2 * math.atan(0.5), delta=1e-4)
            self.assertAlmostEqual(wall[0][1], 0.0, delta=1e-4)
            self.assertAlmostEqual(wall[-1][0], math.pi / 2, delta=2e-3)
            self.assertAlmostEqual(wall[-1][1], math.pi / 2, delta=2e-3)
            # The wall is timelike: it rises faster than it moves across, and it stays outside the light cone T = X.
            for (x0, t0), (x1, t1) in zip(wall, wall[1:]):
                self.assertGreater(t1 - t0, abs(x1 - x0) - 1e-4)
                self.assertGreaterEqual(x1 + 1e-4, t1)
            # On it tan p tan q is constant, the hyperboloid rho = rho_bar: -1/4 for both bubbles.
            for x, t in wall[::25]:
                p, q = (t - x) / 2, (t + x) / 2
                if q < math.pi / 2 - 1e-3:
                    self.assertAlmostEqual(math.tan(p) * math.tan(q), -0.25, delta=5e-3)


if __name__ == "__main__":
    unittest.main()
