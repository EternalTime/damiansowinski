"""Bonnor and Vaidya's charged shell, from the numbers in its files alone.

The readers of the published files and MFS_SPACETIME are those of test_build_mfs_data.py."""
import json
import math
import unittest

import build_mfs_data as build
from test_build_mfs_data import BV_Q, BV_RM, BV_RP, BV_TURN, bv_crossing, checks


@checks("bonnor_vaidya")
class BonnorVaidyaShell(unittest.TestCase):
    """The published drawings of Bonnor and Vaidya's charged shell against its closed forms."""

    @classmethod
    def setUpClass(cls):
        root = build.METRICS_DIR.parent
        cls.diagrams = json.loads((root / "diagrams" / "bonnor_vaidya.json").read_text(encoding="utf-8"))
        cls.conformal = json.loads((root / "conformal" / "bonnor_vaidya.json").read_text(encoding="utf-8"))
        cls.embedding = json.loads((root / "embedding" / "bonnor_vaidya.json").read_text(encoding="utf-8"))

    def marker(self, system, kind):
        (view,) = self.diagrams["systems"][system]
        return view, [m for m in view["markers"] if m["kind"] == kind]

    def test_the_horizons_stand_at_the_roots_of_one_minus_2M_over_r_plus_q_squared_over_r_squared(self):
        for system in ("eddington_finkelstein_ingoing", "eddington_finkelstein_outgoing"):
            view, (grr,) = self.marker(system, "grr")
            X0, X1 = view["box"][:2]
            radii = {round(X0 + x * (X1 - X0), 2) for line in grr["lines"] for x, _ in line}
            self.assertEqual({min(radii), max(radii)}, {BV_RM, BV_RP}, system)
            self.assertAlmostEqual(BV_RP * BV_RM, BV_Q ** 2)
            self.assertAlmostEqual(BV_RP + BV_RM, 2)

    def test_the_mark_on_the_shell_is_where_its_energy_density_changes_sign(self):
        """The jump of G_vv across the shell is 2(M r - q^2/2)/r^3, zero at r = q^2/2M, which lies inside r_-."""
        for system, sign in (("eddington_finkelstein_ingoing", -1), ("eddington_finkelstein_outgoing", 1)):
            view, (mark,) = self.marker(system, "mark")
            X0, X1, Y0, Y1 = view["box"]
            (x, y), = mark["points"]
            r, other = X0 + x * (X1 - X0), Y0 + y * (Y1 - Y0)
            self.assertAlmostEqual(r, BV_TURN, delta=2e-4, msg=system)
            self.assertAlmostEqual(other, sign * BV_TURN, delta=4e-4, msg=system)       # on the shell, v = 0 or u = 0
            self.assertAlmostEqual(r - BV_Q ** 2 / 2, 0, delta=2e-4)
            self.assertLess(r, BV_RM)

    def test_the_event_horizon_leaves_the_centre_at_minus_twice_r_plus(self):
        view, (event,) = self.marker("eddington_finkelstein_ingoing", "event")
        X0, X1, Y0, Y1 = view["box"]
        ends = sorted((X0 + x * (X1 - X0), Y0 + y * (Y1 - Y0)) for line in event["lines"] for x, y in (line[0], line[-1]))
        self.assertAlmostEqual(ends[0][0], 0, delta=1e-3)
        self.assertAlmostEqual(ends[0][1], -2 * BV_RP, delta=2e-3)
        self.assertAlmostEqual(ends[-1][0], BV_RP, delta=1e-3)

    def test_the_homothetic_rays_that_keep_their_radius_are_the_roots_of_the_cubic(self):
        """2 mu R^3 - M R^2 + 2 M^2 R - M Q^2 = 0 at mu = 1/20, and R^2 - 2MR + Q^2 = 0 for the trapped spheres."""
        view, (grr,) = self.marker("homothetic", "grr")
        X0, X1 = view["box"][:2]
        roots = sorted({round(X0 + x * (X1 - X0), 3) for line in grr["lines"] for x, _ in line})
        self.assertEqual(len(roots), 3)
        for R in roots:
            self.assertAlmostEqual(0.1 * R ** 3 - R ** 2 + 2 * R - BV_Q ** 2, 0, delta=2e-2)
        _, (apparent,) = self.marker("homothetic", "apparent")
        self.assertEqual(sorted({round(X0 + x * (X1 - X0), 2) for line in apparent["lines"] for x, _ in line}), [BV_RM, BV_RP])
        _, (shell,) = self.marker("homothetic", "shell")
        self.assertEqual({round(X0 + x * (X1 - X0), 3) for line in shell["lines"] for x, _ in line}, {round(BV_Q ** 2, 3)})

    def test_each_embedded_moment_is_a_flat_disc_inside_the_shell_and_the_closed_form_outside(self):
        """z = 2s - 2q arctan(s/q) with s = sqrt(2r - q^2), from the shell at r = -(v - r) to the rim at
        z = 0, and the last moment has the shell on q^2/2M, where the surface lies level."""
        def height(r):
            s = math.sqrt(max(2 * r - BV_Q ** 2, 0.0))
            return 2 * s - 2 * BV_Q * math.atan(s / BV_Q)
        (view,) = self.embedding["views"]
        frames = view["movie"]["frames"]
        self.assertEqual(len(frames), 42)
        for frame in frames:
            T = float(frame["value"])
            inside, outside = ([[float(x) for x in point] for point in piece["points"]] for piece in frame["pieces"])
            self.assertAlmostEqual(inside[-1][0], -T, delta=1e-9)
            self.assertTrue(all(z == inside[0][2] and abs(rho - r) < 1e-6 for r, rho, z in inside))
            self.assertEqual(outside[-1][0], 4.0)
            for r, rho, z in outside:
                self.assertAlmostEqual(rho, r, delta=1e-6)
                self.assertAlmostEqual(z, height(r) - height(4.0), delta=2e-6, msg=f"v - r = {T} at r = {r}")
        self.assertAlmostEqual(float(frames[-1]["value"]), -BV_TURN, delta=1e-9)
        last = [[float(x) for x in point] for point in frames[-1]["pieces"][1]["points"]]
        self.assertLess((last[1][2] - last[0][2]) / (last[1][1] - last[0][1]), 0.05)

    def test_the_turned_shell_leaves_from_the_mark_and_reaches_the_next_null_infinity(self):
        """On the conformal diagram of Ori's reading the shell is two straight null lines that meet at the
        mark, the first rising to the left and the second to the right, and the centre X = 0 is regular
        from i^- to the last i^+: the singular line is not on it."""
        view = next(v for v in self.conformal["views"] if v["id"] == "bounce")
        first, second = [layer["points"] for layer in view["layers"] if layer["class"] == "surface"]
        (mark,) = [layer["at"] for layer in view["layers"] if layer["class"] == "mark"]
        self.assertEqual(first[-1], mark)
        self.assertEqual(second[0], mark)
        for (Xa, Ta), (Xb, Tb), way in ((first[0], first[-1], -1), (second[0], second[-1], 1)):
            self.assertAlmostEqual(Tb - Ta, way * (Xb - Xa), delta=2e-4)
            self.assertGreater(Tb, Ta)
        self.assertAlmostEqual(second[-1][0] + second[-1][1], 3 * math.pi, delta=2e-4)      # q = 3 pi/2
        p = (mark[1] - mark[0]) / 2 + math.pi / 2
        self.assertAlmostEqual(p, bv_crossing(-2 * BV_TURN), delta=2e-4)
        (centre,) = [layer["points"] for layer in view["layers"] if layer["class"] == "centre"]
        self.assertEqual([X for X, _ in centre], [0, 0])
        self.assertAlmostEqual(centre[-1][1] - centre[0][1], 5 * math.pi, delta=2e-4)
        (singular,) = [layer["points"] for layer in view["layers"] if layer["class"] == "singular"]
        self.assertTrue(all(abs(X - math.pi) < 2e-4 for X, _ in singular))


if __name__ == "__main__":
    unittest.main()
