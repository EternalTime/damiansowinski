"""Misner's and Brill and Lindquist's initial data for two black holes, one moment of a spacetime:
.venv.noindex/bin/python -m unittest discover -s _tools

The physics its texts and drawings state, held on the published files and on the finder of
minimal surfaces the drawings use. Brill and Lindquist's psi is harmonic and reads Newton's
potential energy off its three sheets; Misner's sum solves his equation on the cylinder of mu and
the 2-sphere, is his method of images in the flat coordinates, and is carried onto itself by the
reflection in each throat; a hole alone has the sphere of radius alpha for its throat; each
throat's area is stationary and lies below 16 pi times the square of the mass read beyond it; and
a surface surrounds two equal holes below the separation 1.532 in units of 2 alpha. The embedding
diagram's three views are then held to their closed forms from the numbers written. The finder
needs numpy and scipy, and its tests are skipped where they are absent.
"""
import json
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))
try:
    import numpy as np
    import two_holes
except ImportError:
    two_holes = None

CRITICAL = 1.5324      # the half separation, in alpha, below which a surface surrounds both equal holes
INSIDE = 1.4126        # and below which the sphere through both holes lies inside that surface


def brill_lindquist(x, y, z, a, alpha=(1.0, 1.0)):
    return 1 + alpha[0] / math.sqrt(x * x + y * y + (z - a) ** 2) + alpha[1] / math.sqrt(x * x + y * y + (z + a) ** 2)


def misner(mu, eta, mu0, terms=60):
    return sum((math.cosh(mu + 2 * n * mu0) - math.cos(eta)) ** -0.5 for n in range(-terms, terms + 1))


def Psi(mu, a):
    """Brill and Lindquist's conformal factor on the sphere eta = pi/2, alpha = 1."""
    return math.cosh(mu) ** -0.5 + math.sqrt(2) * math.cosh(mu / 2) / a


class TheConstraint(unittest.TestCase):
    def test_brill_and_lindquists_psi_is_harmonic(self):
        h = 1e-3
        for x, y, z in ((0.3, -0.2, 0.4), (1.5, 0.7, -2.0), (-0.6, 0.1, 1.9)):
            def f(dx=0.0, dy=0.0, dz=0.0):
                return brill_lindquist(x + dx, y + dy, z + dz, 1.25, (0.8, 1.3))
            laplacian = (f(h) + f(-h) + f(0, h) + f(0, -h) + f(0, 0, h) + f(0, 0, -h) - 6 * f()) / h ** 2
            self.assertLess(abs(laplacian), 1e-5)

    def test_misners_sum_solves_his_equation_and_repeats(self):
        mu0, h = 1.2, 1e-3
        for mu, eta in ((0.3, 1.1), (-0.8, 2.3), (1.0, 0.6)):
            def f(dm=0.0, de=0.0):
                return misner(mu + dm, eta + de, mu0)
            residual = ((f(h) + f(-h) - 2 * f()) / h ** 2 + (f(0, h) + f(0, -h) - 2 * f()) / h ** 2
                        + (f(0, h) - f(0, -h)) / (2 * h) / math.tan(eta) - f() / 4)
            self.assertLess(abs(residual), 1e-5)
            # The sum comes back to itself under mu -> mu + 2 mu_0, which closes the wormhole, and under
            # the reflection mu -> 2 mu_0 - mu in the throat mu = mu_0, which joins the two sheets.
            self.assertLess(abs(f(2 * mu0) - f()), 1e-12)
            self.assertLess(abs(misner(2 * mu0 - mu, eta, mu0) - f()), 1e-12)

    def test_misners_sum_is_his_method_of_images(self):
        # psi = Psi sqrt(cosh mu - cos eta) is 1 plus poles of strength a/sinh(n mu_0) at z = +-a coth(n mu_0).
        a, mu0 = 1.0, 1.2
        for mu, eta in ((0.3, 1.1), (-0.8, 2.3), (0.9, 0.4)):
            D = math.cosh(mu) - math.cos(eta)
            rho, z = a * math.sin(eta) / D, a * math.sinh(mu) / D
            images = 1 + sum(a / math.sinh(n * mu0) * (1 / math.hypot(rho, z - a / math.tanh(n * mu0))
                                                         + 1 / math.hypot(rho, z + a / math.tanh(n * mu0)))
                             for n in range(1, 60))
            self.assertLess(abs(misner(mu, eta, mu0) * math.sqrt(D) - images), 1e-10)

    def test_brill_and_lindquists_psi_in_misners_coordinates(self):
        a, alpha = 1.3, (0.8, 1.1)
        for mu, eta in ((0.3, 1.1), (-0.8, 2.3), (2.0, 0.4)):
            D = math.cosh(mu) - math.cos(eta)
            rho, z = a * math.sin(eta) / D, a * math.sinh(mu) / D
            here = D ** -0.5 + (alpha[0] * math.exp(mu / 2) + alpha[1] * math.exp(-mu / 2)) / (math.sqrt(2) * a)
            self.assertLess(abs(brill_lindquist(rho, 0.0, z, a, alpha) / math.sqrt(D) - here), 1e-12)


class TheThreeSheets(unittest.TestCase):
    def test_each_sheet_has_brill_and_lindquists_mass(self):
        a, alpha = 1.5, (0.8, 1.3)
        total, first, second = 2 * (alpha[0] + alpha[1]), 2 * alpha[0] * (1 + alpha[1] / (2 * a)), \
            2 * alpha[1] * (1 + alpha[0] / (2 * a))
        # Far out on the shared sheet psi = 1 + M/2r.
        r = 1e6
        self.assertLess(abs(2 * r * (brill_lindquist(r, 0, 0, a, alpha) - 1) - total), 1e-5)
        # Beyond the hole at z = a the radius r' = alpha_1^2/r_1 is flat far away, where the metric is
        # (alpha_1 psi/r')^4 times flat space in r', and that factor is 1 + M_1/2r'.
        for hole, sign, mass in ((0, 1, first), (1, -1, second)):
            far = 1e6
            r1 = alpha[hole] ** 2 / far
            factor = alpha[hole] / far * brill_lindquist(r1, 0, sign * a, a, alpha)
            self.assertLess(abs(2 * far * (factor - 1) - mass) / mass, 1e-5)

    def test_the_masses_differ_by_newtons_potential_energy(self):
        for a in (5.0, 50.0, 500.0):
            alpha = (0.8, 1.3)
            first, second = 2 * alpha[0] * (1 + alpha[1] / (2 * a)), 2 * alpha[1] * (1 + alpha[0] / (2 * a))
            binding = 2 * (alpha[0] + alpha[1]) - first - second
            self.assertLess(binding, 0)
            newton = -first * second / (2 * a)
            self.assertLess(abs(binding / newton - 1), 3 / a)


@unittest.skipIf(two_holes is None, "needs numpy and scipy")
class TheThroats(unittest.TestCase):
    def test_a_hole_alone_has_the_sphere_of_radius_alpha_for_its_throat(self):
        top = two_holes.throat(3.0, (1.0, 0.0))
        self.assertLess(abs(top - 4.0), 1e-8)
        path, end = two_holes.curve(3.0, top, "throat", (1.0, 0.0))
        for s in np.linspace(0.05, end - 0.05, 9):
            rho, z, _ = path(s)
            self.assertLess(abs(math.hypot(rho, z - 3.0) - 1.0), 1e-7)
        # Its area is 16 pi (2 alpha)^2, that of the horizon of the mass 2 alpha.
        self.assertLess(abs(two_holes.area(3.0, top, "throat", (1.0, 0.0)) / (16 * math.pi * 4) - 1), 1e-6)

    def scaled_area(self, a, top, factor):
        """The area of the throat's curve stretched by `factor` about its hole."""
        path, end = two_holes.curve(a, top, "throat")
        s = np.linspace(0, end, 4001)
        rho, z, _ = path(s)
        rho, z = factor * rho, a + factor * (z - a)
        f = np.array([two_holes.psi(a, r, zz)[0] ** 4 * r for r, zz in zip(rho, z)]) * factor
        h = s[1] - s[0]
        return 2 * math.pi * h / 3 * (f[0] + f[-1] + 4 * f[1:-1:2].sum() + 2 * f[2:-1:2].sum())

    def test_each_throat_is_a_minimal_surface_below_the_penrose_bound(self):
        for a in (0.5, 1.0, 2.0):
            top = two_holes.throat(a)
            self.assertLess(abs(two_holes.throat_miss(a, top)), 1e-7)
            area = self.scaled_area(a, top, 1.0)
            # Stationary: stretched or shrunk by a thousandth, the area changes at second order only.
            for factor in (0.999, 1.001):
                self.assertLess(abs(self.scaled_area(a, top, factor) / area - 1), 2e-5)
            self.assertLess(abs((self.scaled_area(a, top, 1.001) - self.scaled_area(a, top, 0.999)) / area), 2e-6)
            # Its area lies below 16 pi M_1^2, with M_1 the mass read beyond it, and nears it far apart.
            bound = 16 * math.pi * two_holes.masses(a)[1] ** 2
            self.assertLess(area, bound)
        far = two_holes.throat(8.0)
        self.assertLess(abs(two_holes.area(8.0, far, "throat") / (16 * math.pi * two_holes.masses(8.0)[1] ** 2) - 1), 1e-5)

    def test_a_surface_surrounds_two_equal_holes_below_the_separation_1_532(self):
        self.assertLess(abs(two_holes.critical() - CRITICAL), 2e-4)
        self.assertIsNone(two_holes.common(1.54))
        for a in (0.5, 1.0, 1.5):
            top = two_holes.common(a)
            self.assertLess(abs(two_holes.common_miss(a, top)), 1e-8)
            # It lies outside both throats and below the Penrose bound of the total mass, 4 alpha.
            self.assertGreater(top, two_holes.throat(a))
            self.assertLess(two_holes.area(a, top, "common"), 16 * math.pi * 16)
        # With the holes on top of each other it is the throat of the mass 4 alpha, the sphere of radius 2 alpha.
        top = two_holes.common(0.01)
        self.assertLess(abs(top - 2.0), 1e-3)
        self.assertLess(abs(two_holes.area(0.01, top, "common") / (16 * math.pi * 16) - 1), 1e-6)

    def test_the_sphere_through_both_holes_lies_inside_it_below_1_41(self):
        self.assertIsNone(two_holes.common_crossing(1.40))
        self.assertIsNotNone(two_holes.common_crossing(1.42))
        self.assertIsNone(two_holes.common_crossing(1.54))


class Embedding(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.views = {v["id"]: v for v in json.loads(
            (DATA / "embedding" / "misner_brill_lindquist.json").read_text(encoding="utf-8"))["views"]}

    def test_the_views_are_these(self):
        self.assertEqual(list(self.views), ["through", "between", "one_hole"])
        self.assertEqual(self.views["one_hole"]["system"], "isotropic")
        for view in ("through", "between"):
            self.assertNotIn("system", self.views[view])
            values = [f["value"] for f in self.views[view]["movie"]["frames"]]
            self.assertEqual((values[0], values[-1], len(values)), (0.5, 2.0, 31))

    def test_the_sphere_through_both_holes_is_drawn_to_its_metric(self):
        for frame in self.views["through"]["movie"]["frames"]:
            a = frame["value"]
            lower, middle, upper = frame["pieces"]
            self.assertEqual([p["class"] for p in frame["pieces"]], ["sheet2", "sheet", "sheet2"])
            points = lower["points"] + middle["points"][1:] + upper["points"][1:]
            # The pieces join at the throats' circles, and the surface is its own mirror image in z = 0.
            # Each piece is rounded to below a part in ten million of its own extent.
            for here, there in ((lower["points"][-1], middle["points"][0]), (middle["points"][-1], upper["points"][0])):
                self.assertLess(max(abs(here[1] - there[1]), abs(here[2] - there[2])), 2e-6)
            self.assertLess(abs(points[0][2] + points[-1][2]), 1e-5)
            for mu, rho, z in points:
                self.assertLess(abs(rho - a * Psi(mu, a) ** 2), 2e-6, f"a = {a}: rho at mu = {mu}")
            # Each chord is as long as the metric says, a Psi^2 dmu along the meridian.
            for (m0, r0, z0), (m1, r1, z1) in zip(points, points[1:]):
                proper = (m1 - m0) * a * (Psi(m0, a) ** 2 + 4 * Psi((m0 + m1) / 2, a) ** 2 + Psi(m1, a) ** 2) / 6
                self.assertLess(abs(math.hypot(r1 - r0, z1 - z0) / proper - 1), 3e-4, f"a = {a}: the chord at mu = {m0}")
            # The far sheets flatten: at each rim the radius grows nearly as fast as the distance.
            (m0, r0, z0), (m1, r1, z1) = points[-2], points[-1]
            self.assertGreater((r1 - r0) / math.hypot(r1 - r0, z1 - z0), 0.3)
            # The circle mu = 0 bulges out between two narrower circles while a > 1/sqrt(2), and below that
            # it is the narrowest circle of the whole surface.
            k = next(i for i, p in enumerate(points) if p[0] == 0.0)
            beside = (points[k - 1][1], points[k + 1][1])
            if a > 1 / math.sqrt(2):
                self.assertGreater(points[k][1], max(beside))
            else:
                self.assertEqual(points[k][1], min(p[1] for p in points))

    def test_the_marked_circles_are_the_throats_and_the_horizon(self):
        for frame in self.views["through"]["movie"]["frames"]:
            a = frame["value"]
            rings = {}
            for ring in frame["rings"]:
                rings.setdefault(ring["class"], []).append(ring)
            self.assertEqual(len(rings["throat"]), 2)
            self.assertAlmostEqual(rings["throat"][0]["x"], -rings["throat"][1]["x"], delta=1e-12)
            self.assertEqual("horizon" in rings, INSIDE < a < CRITICAL, f"a = {a}")
            if two_holes is not None:
                rho, z = two_holes.crossing(a, two_holes.throat(a))
                self.assertLess(abs(math.atanh(z / a) - rings["throat"][1]["x"]), 1e-7)
                if "horizon" in rings:
                    rho, z = two_holes.common_crossing(a)
                    self.assertLess(abs(math.atanh(z / a) - rings["horizon"][1]["x"]), 1e-7)

    def test_the_midplane_is_drawn_to_its_metric_with_the_horizon_below_1_532(self):
        for frame in self.views["between"]["movie"]["frames"]:
            a = frame["value"]
            plane, = frame["pieces"]

            def psi(r):
                return 1 + 2 / math.sqrt(r * r + a * a)
            self.assertEqual(plane["points"][0], [0.0, 0.0, 0.0])
            for r, rho, z in plane["points"]:
                self.assertLess(abs(rho - r * psi(r) ** 2), 2e-6)
            for (x0, r0, z0), (x1, r1, z1) in zip(plane["points"], plane["points"][1:]):
                proper = (x1 - x0) * (psi(x0) ** 2 + 4 * psi((x0 + x1) / 2) ** 2 + psi(x1) ** 2) / 6
                self.assertLess(abs(math.hypot(r1 - r0, z1 - z0) / proper - 1), 3e-4)
            horizon = [ring for ring in frame["rings"] if ring["class"] == "horizon"]
            self.assertEqual(len(horizon), 1 if a < CRITICAL else 0, f"a = {a}")
            if horizon and two_holes is not None:
                self.assertLess(abs(two_holes.waist(a) - horizon[0]["x"]), 1e-7)

    def test_one_hole_is_flamms_paraboloid_on_both_sheets(self):
        whole, = self.views["one_hole"]["surfaces"][0]["pieces"]
        for r, rho, z in whole["points"]:
            R = r * (1 + 1 / (4 * r)) ** 2
            self.assertLess(abs(rho - R), 2e-6)
            self.assertLess(abs(z - math.copysign(2 * math.sqrt(max(R - 1, 0)), r - 0.25)), 2e-5)
        # The inversion r -> r_s^2/16r carries each marked circle onto another, and the throat onto itself.
        marked = sorted(ring["x"] for ring in self.views["one_hole"]["surfaces"][0]["rings"])
        for x in marked:
            self.assertTrue(any(abs(1 / (16 * x) - y) < 1e-12 for y in marked), x)
        throat, = [ring for ring in self.views["one_hole"]["surfaces"][0]["rings"] if ring["class"] == "throat"]
        self.assertEqual((throat["x"], throat["rho"], throat["z"]), (0.25, 1.0, 0.0))


if __name__ == "__main__":
    unittest.main()
