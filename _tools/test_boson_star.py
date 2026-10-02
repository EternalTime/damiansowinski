"""The boson star as _tools/derivations/boson_star.py solves it, held to the physics: the one
field equation its construction does not use, Kaup's limit, the first law dM = omega dN, the
Newtonian limit, the regular centre, and the same star in the isotropic chart. Then the files
its drawings wrote, from their numbers alone: every ray drawn is a null curve of the star's
metric, and the embedded surface climbs at sqrt(a^2 - 1). It needs sympy and scipy, so it runs
under /tmp/mfs-venv/bin/python and is skipped under a Python without them."""
import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))
try:
    import numpy as np
    import boson_star as bs
except ImportError:
    bs = None

DATA = ROOT / "MFS" / "assets" / "data"


def tortoise(star, r):
    """r* = int_0^r (a/alpha) dr, by Simpson's rule on the solver's own functions."""
    out = []
    for top in np.atleast_1d(r):
        x = np.linspace(0.0, float(top), 4001)
        v = star.values(x)
        f = v["a"][0] / v["alpha"][0]
        h = x[1] - x[0]
        out.append(h / 3 * (f[0] + f[-1] + 4 * f[1:-1:2].sum() + 2 * f[2:-1:2].sum()))
    return np.array(out)


@unittest.skipIf(bs is None, "needs sympy and scipy")
class BosonStar(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.star = bs.star()
        cls.iso = cls.star.isotropic()

    def test_the_star_solves_the_field_equation_its_construction_did_not_use(self):
        # G^theta_theta is the Bianchi identity's consequence of the three equations integrated.
        x = np.linspace(1e-3, 0.99 * self.star.edge, 2000)
        self.assertLess(float(np.max(np.abs(self.star.theta_theta(x)))), 1e-10)

    def test_it_is_the_star_the_drawings_state(self):
        self.assertEqual(self.star.stated(), {})
        self.assertEqual(round(self.star.radius(), 2), 7.86)
        compactness, where = self.star.compactness()
        self.assertEqual((round(compactness, 2), round(where, 1)), (0.24, 3.8))
        # No sphere is trapped: the star has no horizon.
        self.assertLess(compactness, 1)

    def test_the_declared_star_is_the_heaviest_kaups_limit(self):
        field, mass = bs.heaviest()
        self.assertEqual(round(mass, 3), bs.KAUP_MASS)
        self.assertLess(abs(field - bs.SIGMA_C), 0.002)
        self.assertGreater(self.star.M, bs.star(0.251).M)
        self.assertGreater(self.star.M, bs.star(0.291).M)

    def test_the_first_law(self):
        # Along the family of ground states dM = omega dN: adding one boson of the star's own
        # frequency is what a neighbouring star costs.
        middle, step = 0.2, 0.002
        up, down = bs.star(middle + step), bs.star(middle - step)
        self.assertAlmostEqual((up.M - down.M) / (up.N - down.N), bs.star(middle).omega, places=4)

    def test_the_stable_stars_are_bound(self):
        # M < m N on the branch below Kaup's limit: the bosons weigh less together than apart.
        for field in (0.1, 0.2, bs.SIGMA_C):
            star = bs.star(field)
            self.assertLess(star.M, star.N)
            self.assertLess(star.omega, 1)

    def test_the_newtonian_limit(self):
        # For a weak field the equations are Schroedinger and Poisson's, whose solutions scale:
        # sigma_c -> sigma_c/4 halves the mass, quarters 1 - omega and doubles the radius.
        strong, weak = bs.star(0.04), bs.star(0.01)
        self.assertAlmostEqual(weak.M / strong.M, 0.5, delta=0.04)
        self.assertAlmostEqual((1 - weak.omega) / (1 - strong.omega), 0.25, delta=0.01)
        self.assertAlmostEqual(weak.radius() / strong.radius(), 2.0, delta=0.05)

    def test_the_centre_is_regular(self):
        star = self.star
        v = star.values([0.0])
        self.assertEqual((v["a"][0][0], v["a"][1][0], v["alpha"][1][0]), (1.0, 0.0, 0.0))
        self.assertEqual(round(v["alpha"][0][0], 3), bs.STATED["alpha_c"][0])
        # The series about the centre and the integration agree where they meet.
        for name in star.funcs:
            below, above = star.values(star.core * (1 - 1e-12))[name], star.values(star.core * (1 + 1e-12))[name]
            for k in range(3):
                self.assertAlmostEqual(float(below[k][0]), float(above[k][0]), places=10)

    def test_each_derivative_is_the_derivative(self):
        h = 1e-4
        for solver, x in ((self.star, np.array([0.03, 0.5, 3.0, 7.0, 20.0, 39.0, 60.0])),
                          (self.iso, np.array([0.03, 0.5, 3.0, 7.0, 20.0, 37.0, 60.0]))):
            for name in solver.funcs:
                v, up, down = (solver.values(x + d)[name] for d in (0.0, h, -h))
                for k in (0, 1):
                    miss = np.max(np.abs((up[k] - down[k]) / (2 * h) - v[k + 1]))
                    self.assertLess(float(miss), 1e-7, f"{type(solver).__name__} {name}, derivative {k + 1}")
                # The table the drawings read follows the functions.
                for k, column in enumerate(solver.at(name, 0 * x, x)):
                    self.assertLess(float(np.max(np.abs(column - v[k]))), 1e-10)

    def test_far_away_the_metric_is_schwarzschilds(self):
        star, M = self.star, self.star.M
        self.assertLess(abs(star.left), 1e-8)
        x = np.array([star.edge * (1 - 1e-10), star.edge * (1 + 1e-10)])
        for name in star.funcs:
            inner, outer = star.values(x)[name][0]
            self.assertAlmostEqual(float(inner), float(outer), places=8)
        v = star.values([200.0])
        self.assertAlmostEqual(float(v["alpha"][0][0] * v["a"][0][0]), 1.0, places=12)
        self.assertAlmostEqual(float(v["alpha"][0][0] ** 2), 1 - 2 * M / 200.0, places=12)

    def test_the_isotropic_chart_is_the_same_star(self):
        iso, star = self.iso, self.star
        R = np.array([0.0, 0.01, 0.5, 2.0, 7.0, 20.0, 36.0, 80.0])
        v = iso.values(R)
        psi, dpsi = v["psi"][0], v["psi"][1]
        r = psi ** 2 * R
        # r = psi^2 R is the areal radius, and the two charts' radial functions agree there.
        self.assertLess(float(np.max(np.abs(iso.r_of(R) - r))), 1e-10)
        self.assertLess(float(np.max(np.abs(psi / (psi + 2 * R * dpsi) - star.values(r)["a"][0]))), 1e-10)
        self.assertLess(float(np.max(np.abs(v["alpha"][0] - star.values(r)["alpha"][0]))), 1e-12)
        # Far away psi = 1 + M/2R, so the mass read off the conformal factor is the star's.
        self.assertAlmostEqual(float(2 * 80.0 * (v["psi"][0][-1] - 1)), star.M, places=10)
        self.assertEqual(round(float(v["psi"][0][0]), 2), 1.17)


@unittest.skipIf(bs is None, "needs sympy and scipy")
class Drawings(unittest.TestCase):
    """The files the generators wrote, measured against the solver from their numbers alone."""

    @classmethod
    def setUpClass(cls):
        cls.star = bs.star()

    def rays(self, system, view_id):
        view = next(v for v in json.loads((DATA / "diagrams" / "boson_star.json").read_text())["systems"][system]
                    if v["id"] == view_id)
        x_lo, x_hi, y_lo, y_hi = view["box"]
        for family in ("P", "M"):
            for line in view["rays"][family]:
                points = np.array(line)
                yield x_lo + points[:, 0] * (x_hi - x_lo), y_lo + points[:, 1] * (y_hi - y_lo)

    def test_every_ray_of_the_areal_chart_is_null(self):
        # c dt = +-(a/alpha) dr, so along a ray |ct| moves by the tortoise coordinate's change.
        count = 0
        for r, ct in self.rays("areal", "radial"):
            keep = r > 0
            if keep.sum() < 2:
                continue
            rs = tortoise(self.star, r[keep])
            self.assertLess(float(np.max(np.abs(np.abs(np.diff(ct[keep])) - np.abs(np.diff(rs))))), 0.02)
            count += 1
        self.assertGreater(count, 40)

    def test_every_ray_of_the_isotropic_chart_is_null(self):
        iso = self.star.isotropic()
        for R, ct in self.rays("isotropic", "radial"):
            keep = R > 0
            if keep.sum() < 2:
                continue
            rs = tortoise(self.star, iso.r_of(R[keep]))
            self.assertLess(float(np.max(np.abs(np.abs(np.diff(ct[keep])) - np.abs(np.diff(rs))))), 0.02)

    def test_the_marked_sphere_holds_99_percent_of_the_mass(self):
        views = json.loads((DATA / "diagrams" / "boson_star.json").read_text())["systems"]
        iso = self.star.isotropic()
        for system, radius in (("areal", lambda x: x), ("isotropic", lambda x: float(iso.r_of(x)[0]))):
            view = views[system][0]
            mark = next(m for m in view["markers"] if m["kind"] == "surface")
            x = view["box"][0] + mark["lines"][0][0][0] * (view["box"][1] - view["box"][0])
            self.assertAlmostEqual(float(self.star.mass(radius(x))[0]) / self.star.M, 0.99, places=3)

    def test_the_embedded_surface_climbs_at_the_metrics_rate(self):
        surface = json.loads((DATA / "embedding" / "boson_star.json").read_text())["views"][0]["surfaces"][0]
        pieces = {p["id"]: np.array(p["points"]) for p in surface["pieces"]}
        for name in ("core", "tail"):
            r, rho, z = pieces[name].T
            self.assertLess(float(np.max(np.abs(rho - r))), 1e-6)
            for a, b, rise in zip(r, r[1:], np.diff(z)):
                x = np.linspace(a, b, 201)
                slope = np.sqrt(np.maximum(self.star.values(x)["a"][0] ** 2 - 1, 0))
                h = x[1] - x[0]
                want = h / 3 * (slope[0] + slope[-1] + 4 * slope[1:-1:2].sum() + 2 * slope[2:-1:2].sum())
                self.assertAlmostEqual(float(rise), float(want), places=5)
        # The two pieces meet on the sphere of 99% of the mass, and the vacuum paraboloid of the
        # same mass, drawn under them, starts at its throat r = 2M and meets them at the rim.
        self.assertTrue(np.allclose(pieces["core"][-1], pieces["tail"][0], atol=1e-6))
        self.assertAlmostEqual(float(pieces["core"][-1][0]), self.star.radius(), places=9)
        self.assertAlmostEqual(float(pieces["vacuum"][0][0]), 2 * self.star.M, places=9)
        self.assertTrue(np.allclose(pieces["vacuum"][-1], pieces["tail"][-1], atol=1e-6))
        # The star's centre stands above the throat it does not have.
        self.assertGreater(float(pieces["core"][0][2]), float(pieces["vacuum"][0][2]))

    def test_the_conformal_diagram_puts_the_star_in_minkowskis_triangle(self):
        for view in json.loads((DATA / "conformal" / "boson_star.json").read_text())["views"]:
            for layer in view["layers"]:
                if layer["kind"] == "line" and "points" in layer:
                    points = np.array(layer["points"])
                    X, T = points[:, 0], points[:, 1]
                    # Every curve lies in the triangle 0 <= X, |T| + X <= pi.
                    self.assertGreater(float(X.min()), -1e-9)
                    self.assertLess(float(np.max(np.abs(T) + X)), np.pi + 1e-3)


if __name__ == "__main__":
    unittest.main()
