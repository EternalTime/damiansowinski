"""The travelling wave on a string: the ring of free particles on its cone, from the numbers in the file.

The readers of the published files and MFS_SPACETIME are those of test_build_mfs_data.py."""
import math
import unittest

from test_build_mfs_data import checks, embedding_file


@checks("string_wave")
class StringWaveRing(unittest.TestCase):
    """The ring of free particles on the cone of the travelling wave on a string, from the numbers
    in the file: at b = 1/2 and ell = 1 the isotropic chart's geodesics across the string are

        x'' = -rho A'' + (x (x'^2 - y'^2) + 2 y x' y')/(2 rho^2),
        y'' =            (y (y'^2 - x'^2) + 2 x x' y')/(2 rho^2),

    with a prime for d/du, rho^2 = x^2 + y^2 and the declared pulse A = exp(-4u^2)/2, and the
    particle at (x, y) stands on the cone at r = 2 sqrt(rho), which the surface of half angle 30
    degrees draws at the distance r/2 from its axis and the height sqrt(3) r/2."""

    START, STEP = -4.0, 0.001

    @staticmethod
    def rate(u, w):
        x, y, vx, vy = w
        rho2 = x * x + y * y
        pulse = (32 * u * u - 4) * math.exp(-4 * u * u)
        return (vx, vy,
                -math.sqrt(rho2) * pulse + (x * (vx * vx - vy * vy) + 2 * y * vx * vy) / (2 * rho2),
                (y * (vy * vy - vx * vx) + 2 * x * vx * vy) / (2 * rho2))

    def run_to(self, angle, times):
        """The particle that starts at rest at the angle given on rho = 1, at each of the times, by
        Runge and Kutta."""
        w, u, out = (math.cos(angle), math.sin(angle), 0.0, 0.0), self.START, []
        for target in times:
            steps = round((target - u) / self.STEP)
            h = (target - u) / steps
            for _ in range(steps):
                k1 = self.rate(u, w)
                k2 = self.rate(u + h / 2, tuple(a + h * k / 2 for a, k in zip(w, k1)))
                k3 = self.rate(u + h / 2, tuple(a + h * k / 2 for a, k in zip(w, k2)))
                k4 = self.rate(u + h, tuple(a + h * k for a, k in zip(w, k3)))
                w = tuple(a + h * (p + 2 * q + 2 * r + s) / 6 for a, p, q, r, s in zip(w, k1, k2, k3, k4))
                u += h
            out.append(w)
        return out

    def setUp(self):
        self.view, = embedding_file("string_wave")["views"]
        self.frames = {frame["value"]: frame for frame in self.view["movie"]["frames"]}

    def test_every_marked_particle_is_where_its_geodesic_puts_it_on_the_cone(self):
        times = [-1.5, 0.0, 1.0, 2.5]
        self.assertEqual([surface["time"] for surface in self.view["surfaces"]], times)
        for k in range(12):
            for u, (x, y, _, _) in zip(times, self.run_to(math.radians(30 * k), times)):
                rho = math.hypot(x, y)
                r = 2 * math.sqrt(rho)
                want = (r / 2 * x / rho, r / 2 * y / rho, math.sqrt(3) * r / 2)
                for got, expected in zip(self.frames[u]["dots"][k]["at"], want):
                    self.assertAlmostEqual(got, expected, delta=2e-5, msg=f"particle {k} at u = {u}")

    def test_the_ring_stays_on_the_cone_in_every_frame(self):
        for u, frame in self.frames.items():
            curve, = frame["curves"]
            for X, Y, Z in curve["points"]:
                self.assertAlmostEqual(Z, math.sqrt(3) * math.hypot(X, Y), delta=2e-6, msg=f"u = {u}")

    def test_the_caption_states_what_the_ring_does(self):
        def radii(u):
            return [2 * math.hypot(X, Y) for X, Y, _ in self.frames[u]["curves"][0]["points"]]
        caption = " ".join(self.view["caption"])
        self.assertIn("reaching from $r = 1.5\\,\\ell$ to $2.5\\,\\ell$ on the crest", caption)
        self.assertAlmostEqual(min(radii(0.0)), 1.5, delta=1e-3)
        self.assertAlmostEqual(max(radii(0.0)), 2.5, delta=1e-3)
        self.assertIn("left falling toward the string", caption)
        self.assertAlmostEqual(min(radii(-1.5)), 2.0, delta=1e-3)
        self.assertLess(max(radii(1.0)), 1.75)
        self.assertLess(max(radii(2.5)), max(radii(1.0)) - 0.3)


if __name__ == "__main__":
    unittest.main()
