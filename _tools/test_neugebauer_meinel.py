"""Neugebauer and Meinel's rigidly rotating disc of dust:
.venv.noindex/bin/python -m unittest discover -s _tools

What its texts and drawings state, held on the published files and on the solution itself. The
first class needs nothing but the files. The second needs numpy and scipy and is skipped where they
are absent: it holds the theta functions of nm_disc.py, which every drawing is read from, to what
they were never told, Ernst's equation off the disc, the conditions Neugebauer and Meinel imposed
on the disc, their closed form for V_0, the Maclaurin disc at small mu and the extreme Kerr metric
near mu_0, and it holds the numbers the captions quote to the disc of mu = 3.
"""
import importlib.util
import json
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
HAS_NUMBERS = all(importlib.util.find_spec(name) is not None for name in ("numpy", "scipy"))
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))

MU = 3.0
ERGO = (0.15007976, 1.83091417)         # where the ergosurface crosses the plane of the disc, in rho_0
LIGHT = 1.37769858                      # where a point at rest in the turning frame moves at the speed of light


def load(folder):
    return json.loads((DATA / folder / "neugebauer_meinel.json").read_text(encoding="utf-8"))


class Drawings(unittest.TestCase):
    def view(self, system, name):
        return next(v for v in load("diagrams")["systems"][system] if v["id"] == name)

    def test_the_plane_marks_the_rim_and_both_crossings_of_the_ergosurface(self):
        view = self.view("bardeen_wagoner", "plane")
        X0, X1 = view["box"][:2]
        marked = {m["kind"]: sorted(X0 + line[0][0] * (X1 - X0) for line in m["lines"]) for m in view["markers"]}
        kinds = [k for k in marked if any(abs(x - ERGO[0]) < 5e-3 for x in marked[k])]
        self.assertEqual(len(kinds), 1, marked)
        for got, want in zip(marked[kinds[0]], ERGO):
            self.assertAlmostEqual(got, want, delta=5e-3)
        self.assertTrue(any(abs(x - 1.0) < 1e-3 for xs in marked.values() for x in xs), "the rim is marked")
        self.assertEqual(view["quotient"] if "quotient" in view else "phi", "phi")

    def test_rays_cross_the_disc_and_the_ergoregion(self):
        """The rays of no angular momentum run through the dust and through the ergoregion, where
        Bardeen and Wagoner's lapse is real: each family has rays on both sides of the rim."""
        view = self.view("bardeen_wagoner", "plane")
        X0, X1 = view["box"][:2]
        for family in view["rays"].values():
            xs = [X0 + p[0] * (X1 - X0) for ray in family for p in ray]
            self.assertLess(min(xs), 0.2)
            self.assertGreater(max(xs), 3.0)
        self.assertFalse(view.get("hatch"))

    def test_the_turning_frame_is_drawn_out_to_the_speed_of_light(self):
        view = self.view("corotating", "plane")
        self.assertAlmostEqual(view["box"][1], LIGHT, delta=2e-4)
        self.assertLess(view["box"][1], LIGHT)

    def test_the_axis_runs_through_the_centre_of_the_disc(self):
        view = self.view("weyl", "axis")
        X0, X1 = view["box"][:2]
        self.assertEqual((X0, X1), (-3, 3))
        lines = [X0 + line[0][0] * (X1 - X0) for m in view["markers"] for line in m["lines"]]
        self.assertTrue(any(abs(x) < 1e-3 for x in lines))
        (mark,) = view["slices"]
        (at,) = mark["points"]
        self.assertAlmostEqual(at[0], 0.5, places=9)
        self.assertAlmostEqual(at[1], 0.5, places=9)

    def test_the_limit_marks_the_ergosurface_of_extreme_kerr(self):
        """r = m in Bardeen and Wagoner's radius, 2m in Boyer and Lindquist's."""
        view = self.view("black_hole_limit", "equator")
        X0, X1 = view["box"][:2]
        marked = [X0 + line[0][0] * (X1 - X0) for m in view["markers"] for line in m["lines"]]
        self.assertTrue(any(abs(x - 1.0) < 5e-3 for x in marked), marked)
        self.assertNotIn("slices", {k for k, v in view.items() if v})

    def test_the_embedded_plane_passes_between_minkowski_space_and_flat_space_three_times(self):
        view = next(v for v in load("embedding")["views"] if v["id"] == "plane")
        pieces = view["surfaces"][0]["pieces"]
        self.assertEqual([p["id"] for p in pieces], ["centre", "disc", "rim", "plane"])
        self.assertEqual([p.get("space", "flat") for p in pieces], ["minkowski", "flat", "minkowski", "flat"])
        for a, b in zip(pieces, pieces[1:]):
            for x, y in zip(a["points"][-1], b["points"][0]):
                self.assertAlmostEqual(x, y, places=6)
        self.assertEqual(pieces[0]["points"][0][:2], [0.0, 0.0])
        # The rim lies in the narrow band that stands in Minkowski space.
        self.assertLess(pieces[2]["points"][0][0], 1.0)
        self.assertGreater(pieces[2]["points"][-1][0], 1.0)
        self.assertLess(pieces[2]["points"][-1][0] - pieces[2]["points"][0][0], 0.05)

    def test_the_circles_narrow_beyond_the_rim_and_widen_again(self):
        """The beginning of the throat: 4.41 rho_0 at 0.92, 4.24 at the rim, 4.07 at 1.20, wider beyond."""
        view = next(v for v in load("embedding")["views"] if v["id"] == "plane")
        points = [p for piece in view["surfaces"][0]["pieces"] for p in piece["points"]]
        radius = lambda x: min(points, key=lambda p: abs(p[0] - x))[1]  # noqa: E731
        widest = max(p[1] for p in points if p[0] < 1.0)
        narrowest = min(p[1] for p in points if 1.0 < p[0] < 2.0)
        self.assertAlmostEqual(widest, 4.413, delta=2e-3)
        self.assertAlmostEqual(narrowest, 4.066, delta=2e-3)
        self.assertAlmostEqual(radius(1.0), 4.241, delta=2e-3)
        self.assertGreater(radius(6.0), 7.9)

    def test_the_axis_of_the_disc_has_no_horizon_and_the_limit_has_one(self):
        views = {v["id"]: v for v in load("conformal")["views"]}
        disc, limit = views["weyl_axis"], views["limit_axis"]
        classes = lambda v: [layer.get("class") for layer in v["layers"] if layer["kind"] == "line"]  # noqa: E731
        self.assertNotIn("horizon", classes(disc))
        self.assertEqual(classes(disc).count("scri"), 4)
        self.assertEqual(classes(limit).count("horizon"), 2)
        self.assertEqual(classes(limit).count("scri"), 2)
        (mark,) = disc["slices"]
        (at,) = mark["points"]
        self.assertAlmostEqual(at[0], 0.0, places=9)
        self.assertAlmostEqual(at[1], 0.0, places=9)


@unittest.skipUnless(HAS_NUMBERS, "numpy and scipy are needed")
class Solution(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import numpy
        import nm_disc
        cls.np, cls.nm = numpy, nm_disc

    def test_the_theta_functions_solve_ernsts_equation_in_every_arrangement_of_the_branch_points(self):
        """e^(2U) Laplacian(f) = (grad f)^2, by fourth order differences, with the cut of zeta -+ i rho to
        the right of the others, over and under the cut of X_2, and between its ends."""
        nm = self.nm
        h = 1e-2
        slope = lambda v: (v[0] - 8 * v[1] + 8 * v[2] - v[3]) / (12 * h)  # noqa: E731
        bend = lambda v, c: (-v[0] + 16 * v[1] - 30 * c + 16 * v[2] - v[3]) / (12 * h * h)  # noqa: E731
        for mu, rho, zeta in ((1.0, 2.0, 3.0), (3.0, 1.5, 0.8), (3.0, 0.3, 0.08), (0.5, 0.5, 1.0), (4.0, 1.2, 0.05)):
            f0 = nm.Point(mu, rho, zeta).ernst()
            R = [nm.Point(mu, rho + n * h, zeta).ernst() for n in (-2, -1, 1, 2)]
            Z = [nm.Point(mu, rho, zeta + n * h).ernst() for n in (-2, -1, 1, 2)]
            fr, fz = slope(R), slope(Z)
            miss = abs(f0.real * (bend(R, f0) + bend(Z, f0) + fr / rho) - fr ** 2 - fz ** 2)
            self.assertLess(miss, 2e-4 * (abs(fr) ** 2 + abs(fz) ** 2), (mu, rho, zeta))

    def test_the_two_quotients_for_the_dragging_sum_to_two(self):
        """P_+ + P_- = 2, an identity of the theta functions that holds only with theta*(c) where the
        review of 2003 prints theta*(0); and its formula for e^(2U) is the real part of f."""
        nm = self.nm
        for mu, rho, zeta in ((0.5, 2.0, 3.0), (1.0, 0.5, 1.0), (3.0, 1.5, 0.0), (3.0, 0.6, 0.0), (4.5, 0.4, 0.3)):
            p = nm.Point(mu, rho, zeta)
            plus, minus = p.P()
            self.assertAlmostEqual(plus + minus, 2.0, places=8)
            self.assertAlmostEqual(p.F(), p.ernst().real, places=8)
            self.assertLess(p.asymmetry, 1e-8)

    def test_the_dust_is_at_rest_in_the_turning_frame_at_one_potential(self):
        """The condition Neugebauer and Meinel imposed, f' = e^(2 V_0) on the disc, which the theta
        functions are never told: e^(2U') is e^(2 V_0) at every radius, and outside the disc it is not."""
        nm = self.nm
        for mu in (0.5, 1.0, 3.0, 4.5):
            V, Om = nm.e2V0(mu), nm.omega_disc(mu)
            for rho in (0.1, 0.4, 0.7, 0.95, 0.999):
                self.assertAlmostEqual(nm.Point(mu, rho, 0.0).turning(Om) / V, 1.0, places=7)
            self.assertLess(nm.Point(mu, 1.05, 0.0).turning(Om) / V, 0.999)

    def test_the_metric_on_the_disc_follows_from_the_two_conditions(self):
        """A, g_phiphi and g_rhorho on the disc from e^(2U') = e^(2 V_0) and a'_zeta = 0 alone, with
        only f on the disc taken from the theta functions, against the theta functions' own."""
        nm, np = self.nm, self.np
        for mu in (1.0, 3.0):
            r = np.array([0.2, 0.5, 0.8, 0.95])
            F, A, G, grr = nm.plane(mu, r)
            for got, want, places in zip((A, G * r ** 2, grr), nm.disc_conditions(mu, r), (6, 6, 5)):
                self.assertLess(float(np.max(np.abs(got / want - 1))), 10.0 ** -places)

    def test_the_surface_potential_is_the_closed_form_of_1994(self):
        """V_0 = -arsinh(mu + (1 + mu^2)/(p - 2 mu/3))/2 with Weierstrass's function p at I(mu), and
        the parameter relation f_0 conj(f_0) + 4 Omega^2 rho_0^2 = 1."""
        nm = self.nm
        for mu in (0.5, 1.0, 3.0):
            self.assertAlmostEqual(nm.e2V0(mu), math.exp(2 * nm.V0_closed(mu)), places=7)
            self.assertAlmostEqual(abs(nm.centre(mu)) ** 2 + 4 * nm.omega_disc(mu) ** 2, 1.0, places=6)

    def test_the_axis_is_regular_and_space_is_flat_far_away(self):
        nm = self.nm
        for mu in (1.0, 3.0):
            Om = nm.omega_disc(mu)
            near = nm.Point(mu, 0.01, 0.7)
            F, A, gpp = near.metric(Om)
            self.assertLess(abs(A), 1e-3)                                           # a = 0 on the axis
            self.assertAlmostEqual(near.kappa_over_F() * F / nm.kappa_far(mu), 1.0, places=3)     # k = 0
            far = nm.Point(mu, 3000.0, 4000.0)
            self.assertAlmostEqual(far.ernst().real, 1 - 2 * nm.far_field(mu)[0] / 5000.0, places=5)

    def test_a_small_mu_is_the_maclaurin_disc(self):
        """f = 1 + mu f_1 with f_1 = -(1/pi)((4/3) arccot(xi) + (xi - (xi^2 + 1/3) arccot(xi))(1 - 3 eta^2)),
        their (10) of 1995."""
        nm = self.nm
        mu, rho, z = 0.01, 1.5, 0.8
        A = rho * rho + z * z - 1
        xi = math.sqrt((A + math.sqrt(A * A + 4 * z * z)) / 2)
        eta = z / xi
        acot = math.atan2(1, xi)
        f1 = -((4 / 3) * acot + (xi - (xi * xi + 1 / 3) * acot) * (1 - 3 * eta * eta)) / math.pi
        self.assertAlmostEqual((nm.Point(mu, rho, z).ernst().real - 1) / mu, f1, delta=1e-3 * abs(f1))

    def test_near_mu_0_the_field_outside_is_extreme_kerrs(self):
        """f -> (2 Omega R - 1 - i cos(theta))/(2 Omega R + 1 - i cos(theta)) at fixed Omega R, Meinel's
        (33) of 2002, and 2 Omega M -> 1: nearer at mu = 4.62 than at mu = 4.5."""
        nm = self.nm
        theta = math.pi / 3
        miss = []
        for mu in (4.5, 4.62):
            Om = nm.omega_disc(mu)
            R = 1 / Om
            f = nm.Point(mu, R * math.sin(theta), R * math.cos(theta)).ernst()
            kerr = (2 * Om * R - 1 - 1j * math.cos(theta)) / (2 * Om * R + 1 - 1j * math.cos(theta))
            miss.append(abs(f - kerr))
        self.assertLess(miss[0], 1e-2)
        self.assertLess(miss[1], 0.5 * miss[0])
        M, _ = nm.far_field(4.5)
        self.assertAlmostEqual(2 * nm.omega_disc(4.5) * M, 1.0, delta=0.03)
        self.assertLess(nm.e2V0(4.62), 1e-5)

    def test_the_ergoregion_appears_at_the_mu_the_history_states(self):
        nm = self.nm
        self.assertEqual(nm.ergosurface(1.68), [])
        self.assertEqual(len(nm.ergosurface(1.70)), 2)
        with self.assertRaises(ValueError):
            nm.Point(4.63, 1.0, 1.0)

    def test_the_disc_drawn_has_the_numbers_its_captions_quote(self):
        nm, np = self.nm, self.np
        V, Om = nm.e2V0(MU), nm.omega_disc(MU)
        self.assertAlmostEqual(V, 0.0303, places=4)
        self.assertAlmostEqual(1 / math.sqrt(V) - 1, 4.74, places=2)
        self.assertAlmostEqual(math.sqrt(V), 0.174, places=3)
        self.assertAlmostEqual(Om, 0.2133, places=4)
        self.assertAlmostEqual(nm.far_field(MU)[0], 1.716, places=3)
        for got, want in zip(nm.ergosurface(MU), ERGO):
            self.assertAlmostEqual(got, want, places=6)
        self.assertAlmostEqual(nm.light_cylinder(MU), LIGHT, places=6)
        self.assertAlmostEqual(float(nm.profile(MU, "omega", 0.0)) / Om, 0.951, places=3)
        self.assertAlmostEqual(math.sqrt(ERGO[1] ** 2 - 1), 1.534, places=3)
        from scipy.integrate import quad
        self.assertAlmostEqual(quad(lambda z: 1 / float(nm.axis(MU, z)), 0, 1)[0], 16.3, delta=0.05)
        self.assertAlmostEqual(quad(lambda r: math.sqrt(float(nm.profile(MU, "grr", r))), 0, 1)[0], 4.85, delta=0.005)
        self.assertAlmostEqual(float(nm.profile(MU, "F", np.array(1.0))), 1 - MU / 2, places=7)

    def test_the_functions_the_drawings_declare_are_the_tables_and_their_slopes(self):
        """Each declared function against the theta functions at the point itself, and its declared
        derivative against a difference of the function."""
        import null_rays as nr
        nm, np = self.nm, self.np
        Om, far = nm.omega_disc(MU), nm.kappa_far(MU)
        f = nr.DECLARED_FUNCTIONS
        for rho in (0.3, 0.8, 1.2, 2.5):
            p = nm.Point(MU, rho, 0.0)
            F, A, gpp = p.metric(Om)
            self.assertAlmostEqual(float(f["nm_g"]._imp_(rho)), gpp / rho ** 2, places=6)
            self.assertAlmostEqual(float(f["nm_w"]._imp_(rho)), A / gpp, places=7)
            self.assertAlmostEqual(float(f["nm_r"]._imp_(rho)), p.kappa_over_F() / far, places=6)
            self.assertAlmostEqual(float(f["nm_e"]._imp_(rho)), p.turning(Om), places=7)
            for name in ("nm_g", "nm_w", "nm_r", "nm_e", "nm_s"):
                h = 1e-6
                want = (float(f[name]._imp_(rho + h)) - float(f[name]._imp_(rho - h))) / (2 * h)
                self.assertAlmostEqual(float(f[name + "1"]._imp_(rho)), want, delta=1e-5 * (1 + abs(want)))
        self.assertAlmostEqual(float(f["nm_x"]._imp_(0.0)), nm.e2V0(MU), places=7)
        self.assertAlmostEqual(float(nr.NM_OMEGA), Om, places=12)


if __name__ == "__main__":
    unittest.main()
