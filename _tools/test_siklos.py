"""Siklos's waves in anti-de Sitter space, from the numbers in their files alone.

The readers of the published files and MFS_SPACETIME are those of test_build_mfs_data.py."""
import math
import unittest

from test_build_mfs_data import checks, conformal_file, diagram_file, embedding_file


@checks("siklos")
class SiklosWaves(unittest.TestCase):
    """The published diagrams of Siklos's waves, held to the metric from the numbers written and
    nothing else, in units of L."""

    def test_the_wave_front_is_the_hyperbolic_plane_and_its_curves_are_horocycles(self):
        view, = embedding_file("siklos")["views"]
        surface, = view["surfaces"]
        sheet = surface["pieces"][0]
        for r, rho, z in sheet["points"]:
            # The disc's circle of coordinate radius r is a proper distance s = 2 artanh(r/2) from its
            # centre, and on the hyperboloid (Z + 1)^2 - rho^2 = 1 it is rho = sinh s, Z = cosh s - 1.
            s = 2 * math.atanh(r / 2)
            self.assertAlmostEqual(rho, math.sinh(s), delta=1e-6)
            self.assertAlmostEqual(z, math.cosh(s) - 1, delta=1e-6)
        self.assertAlmostEqual(2 * math.atanh(sheet["points"][-1][0] / 2), 2.0, delta=1e-9)
        self.assertEqual(len(surface["curves"]), 3)
        for curve, x in zip(surface["curves"], (0.5, 1.0, 2.0)):
            points = curve["points"]
            for X, Y, Z in points:
                self.assertAlmostEqual((Z + 1) ** 2 - X * X - Y * Y, 1.0, delta=1e-6)
                # Siklos's x on the hyperboloid: with the ideal point x = infinity at (-1, 0) of the
                # disc, L/x = Z + 1 + X, the horocycle's defining height along that null direction.
                self.assertAlmostEqual(1 / (Z + 1 + X), x, delta=1e-6 * (1 + x) ** 2)
            # A chord of a horocycle, measured in Minkowski space, is its arc dy/x.
            chords = [math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 - (a[2] - b[2]) ** 2)
                      for a, b in zip(points, points[1:])]
            self.assertLess(max(chords) - min(chords), 1e-6)

    def test_every_ray_of_kaigorodovs_planes_keeps_its_null_coordinate(self):
        """With a Killing direction divided out the rays keep v -+ (2/5) x^(5/2), in each chart's
        own coordinates, and the Poincare chart's keep ct -+ z_*."""
        def zstar(z, n=400):
            h = z / n
            f = [math.sqrt(1 + (i * h) ** 3 / 2) for i in range(n + 1)]
            return h / 3 * (f[0] + f[-1] + 4 * sum(f[1:-1:2]) + 2 * sum(f[2:-1:2]))
        star = lambda x: 0.4 * x ** 2.5
        kept = {"kaigorodov": lambda v, x: (v, star(x)),
                "kaigorodov_horospheric": lambda v, rho: (v, star(math.exp(-rho))),
                "kaigorodov_homogeneous": lambda U, Z: (-U * math.exp(5 * Z), star(math.exp(2 * Z))),
                "kaigorodov_kundt": lambda V, x: (V * x * x / math.sqrt(2), star(x)),
                "kaigorodov_poincare": lambda t, z: (t, zstar(z))}
        systems = diagram_file("siklos")["systems"]
        checked = 0
        for system, of in kept.items():
            view, = systems[system]
            X0, X1, Y0, Y1 = view["box"]
            (a, b), (c, d) = view["to_display"]
            for family, rays in view["rays"].items():
                for ray in rays:
                    values = []
                    for px, py in ray:
                        if not (0.02 < px < 0.98 and 0.02 < py < 0.98):
                            continue
                        X, Y = X0 + px * (X1 - X0), Y0 + py * (Y1 - Y0)
                        det = a * d - b * c
                        x0, x1 = (d * X - b * Y) / det, (a * Y - c * X) / det
                        time, depth = of(x0, x1)
                        values.append((time - depth, time + depth))
                    if len(values) < 2:
                        continue
                    spreads = [max(v[k] for v in values) - min(v[k] for v in values) for k in (0, 1)]
                    # The published points are rounded to four decimals of the box.
                    scale = max(1.0, max(abs(v[k]) for v in values for k in (0, 1)))
                    self.assertLess(min(spreads), 2e-2 * scale, f"{system} {family}")
                    checked += 1
        self.assertGreater(checked, 40)

    def test_kaigorodovs_conformal_diagram_is_a_triangle_between_a_boundary_and_a_singularity(self):
        views = conformal_file("siklos")["views"]
        self.assertEqual([v["id"] for v in views], ["kaigorodov", "kaigorodov_poincare", "kaigorodov_horospheric",
                                                     "kaigorodov_homogeneous", "kaigorodov_kundt"])
        for view in views:
            self.assertIn("restriction", view)
            lines = {layer["class"]: layer for layer in view["layers"] if layer["kind"] == "line"}
            # The conformal boundary is the timelike line X = 0.
            self.assertTrue(all(abs(p[0]) < 1e-9 for p in lines["boundary"]["points"]))
            # The singularity is the pair of null edges from i^- and i^+ to the corner (pi, 0).
            edges = [layer["points"] for layer in view["layers"] if layer["class"] == "singular"]
            self.assertEqual(len(edges), 2)
            for edge in edges:
                for X, T in edge:
                    self.assertAlmostEqual(X + abs(T), math.pi, delta=1e-3)


if __name__ == "__main__":
    unittest.main()
