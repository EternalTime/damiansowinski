"""Herdeiro and Radu's Kerr black hole with scalar hair, configuration IV, as
_tools/derivations/kerr_scalar_hair.py reads it from the authors' published data, held to the
physics: the field equations the paper prints, which print_charts.py holds to the published Einstein
tensor, the boundary conditions, the mass and angular momentum the authors state, the split of both
between the horizon and the field, and Smarr's relation; Kerr's metric in the same chart, the
paper's appendix A; and then the files its drawings wrote, from their numbers alone. It needs numpy
and scipy, so it runs under /tmp/mfs-venv/bin/python and is skipped under a Python without them."""
import gzip
import hashlib
import json
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))
try:
    import numpy as np
    import kerr_scalar_hair as ksh
except ImportError:
    ksh = None

DATA = ROOT / "MFS" / "assets" / "data"


@unittest.skipIf(ksh is None, "needs numpy and scipy")
class Solution(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.hair = ksh.Hair()

    def test_the_file_is_the_authors_as_published(self):
        with gzip.open(ksh.DATA, "rb") as handle:
            digest = hashlib.sha256(handle.read()).hexdigest()
        self.assertEqual(digest, "cdf37db6c7eec5fcbef74158128cba4bc44ded4fa041d4664b3be493c9703b73")
        X, theta, values = ksh.table()
        self.assertEqual((len(X), len(theta)), (251, 30))
        self.assertEqual((X[0], X[-1]), (0.0, 1.0))
        self.assertAlmostEqual(theta[-1], math.pi / 2, places=12)

    def test_it_solves_the_field_equations_the_paper_prints(self):
        for name, (residual, scale) in ksh.field_equations(self.hair).items():
            with self.subTest(name):
                self.assertLess(residual / scale, 1e-3)

    def test_the_boundary_conditions(self):
        v = self.hair.values
        # On the horizon W is the horizon's angular velocity, w = m Omega_H, and F_0 - F_1 is constant.
        self.assertLess(float(np.max(np.abs(v["W"][:, 0] - ksh.OMEGA / ksh.M_WIND))), 1e-9)
        self.assertLess(float(np.ptp(v["F_0"][:, 0] - v["F_1"][:, 0])), 1e-6)
        # On the axis the field vanishes and F_1 = F_2, so the axis has no conical singularity.
        self.assertLess(float(np.max(np.abs(v["phi"][0]))), 1e-12)
        self.assertLess(float(np.max(np.abs(v["F_1"][0] - v["F_2"][0]))), 1e-5)
        # Far away every function vanishes, the metric is Minkowski's.
        for name in ksh.NAMES:
            self.assertLess(float(np.max(np.abs(v[name][:, -1]))), 1e-8, name)

    def test_the_mass_and_angular_momentum_the_authors_state(self):
        c = self.hair.charges
        self.assertAlmostEqual(c["M"], ksh.PUBLISHED["M"], delta=1e-3)
        self.assertAlmostEqual(c["J"], ksh.PUBLISHED["J"], delta=2e-3)

    def test_the_horizon_and_the_field_share_them_as_the_authors_state(self):
        c, f = self.hair.charges, self.hair.field
        self.assertAlmostEqual(c["J"] - ksh.M_WIND * f["Q"], ksh.PUBLISHED["J_H"], delta=1e-3)
        self.assertAlmostEqual(c["M"] - f["M_Psi"], ksh.PUBLISHED["M_H"], delta=1e-3)

    def test_smarrs_relation(self):
        """M = 2 T_H S + 2 Omega_H (J - mQ) + M_Psi, with S = A_H/4: the horizon's share computed from
        the horizon alone equals the mass less the field's energy."""
        c, f, h = self.hair.charges, self.hair.field, self.hair.horizon_quantities
        horizon = 2 * h["T_H"] * h["A_H"] / 4 + 2 * h["Omega_H"] * (c["J"] - ksh.M_WIND * f["Q"])
        self.assertLess(abs(horizon + f["M_Psi"] - c["M"]), 1e-4)

    def test_each_derivative_is_the_derivative(self):
        for name in ("F_0", "F_1", "F_2", "W"):
            for theta in (0.0, 0.7, math.pi / 2):
                for r in (0.2, 1.0, 4.0):
                    h = 1e-4
                    slope = (self.hair(name, r + h, theta) - self.hair(name, r - h, theta)) / (2 * h)
                    self.assertAlmostEqual(float(self.hair(name, r, theta, dr=1)), float(slope), delta=1e-6)
                    bend = (self.hair(name, r + h, theta, dr=1) - self.hair(name, r - h, theta, dr=1)) / (2 * h)
                    self.assertAlmostEqual(float(self.hair(name, r, theta, dr=2)), float(bend), delta=1e-4 * (1 + abs(bend)))

    def test_the_tortoise_table_is_the_quadrature(self):
        r = np.array([0.1000001, 0.2, 1.0, 5.0, 30.0])
        self.assertLess(float(np.max(np.abs(self.hair.tortoise_fast(r) - self.hair.tortoise(r)))), 1e-6)

    def test_the_ergosurface_on_the_equator(self):
        r = self.hair.ergosurface()
        self.assertAlmostEqual(r, 0.36495, delta=1e-4)
        gtt = lambda x: float(self.hair.metric(np.array(x), math.pi / 2)[0])
        self.assertGreater(gtt(r * 0.9), 0)
        self.assertLess(gtt(r * 1.1), 0)


@unittest.skipIf(ksh is None, "needs numpy and scipy")
class KerrMember(unittest.TestCase):
    """The paper's (A.1) is Kerr's metric: Boyer and Lindquist's with R = r - c_t."""

    def test_the_functions_are_kerrs_metric(self):
        M, J = 0.9330359, 0.7403448
        rH, ct = ksh.kerr_for(M, J)
        a = J / M
        r = np.linspace(rH * 1.01, 20, 50)[:, None]
        theta = np.linspace(0.1, math.pi / 2, 7)[None, :]
        F0, F1, F2, W = ksh.kerr_hr(rH, ct, r, theta)
        N = 1 - rH / r
        s2 = np.sin(theta) ** 2
        gtt = -np.exp(2 * F0) * N + np.exp(2 * F2) * r * r * s2 * W * W
        gtp = -np.exp(2 * F2) * r * r * s2 * W
        R = r - ct
        Sigma = R * R + a * a * np.cos(theta) ** 2
        Delta = R * R - 2 * M * R + a * a
        self.assertLess(float(np.max(np.abs(gtt + (1 - 2 * M * R / Sigma)))), 1e-12)
        self.assertLess(float(np.max(np.abs(gtp + 2 * M * a * R * s2 / Sigma))), 1e-12)
        self.assertLess(float(np.max(np.abs(np.exp(2 * F1) / N - Sigma / Delta))), 1e-11)
        # The horizon r = r_H is Boyer and Lindquist's R_+, and W there is Kerr's Omega_H.
        Rp = M + math.sqrt(M * M - a * a)
        self.assertAlmostEqual(rH - ct, Rp, places=12)
        self.assertAlmostEqual(float(ksh.kerr_hr(rH, ct, np.array(rH), 0.3)[3]), a / (Rp * Rp + a * a), places=12)


@unittest.skipIf(ksh is None, "needs numpy and scipy")
class Drawings(unittest.TestCase):
    """The files the generators wrote, measured against the solution from their numbers alone."""

    @classmethod
    def setUpClass(cls):
        cls.hair = ksh.Hair()
        cls.diagrams = json.loads((DATA / "diagrams" / "kerr_scalar_hair.json").read_text())["systems"]

    def rays(self, system, view_id):
        view = next(v for v in self.diagrams[system] if v["id"] == view_id)
        x_lo, x_hi, y_lo, y_hi = view["box"]
        for family in ("P", "M"):
            for line in view["rays"][family]:
                points = np.array(line)
                yield x_lo + points[:, 0] * (x_hi - x_lo), y_lo + points[:, 1] * (y_hi - y_lo)

    def check_null(self, system, view_id, star, lo):
        """Every ray of the view moves in |ct| by the change of star(r), from just outside the horizon lo."""
        count = 0
        for r, ct in self.rays(system, view_id):
            keep = r > lo
            if keep.sum() < 2:
                continue
            rs = star(r[keep])
            self.assertLess(float(np.max(np.abs(np.abs(np.diff(ct[keep])) - np.abs(np.diff(rs))))), 0.02)
            count += 1
        self.assertGreater(count, 20)

    def test_every_ray_on_the_axis_is_null(self):
        # c dt = +-e^(F_1 - F_0) dr/N, so along a ray |ct| moves by the tortoise coordinate's change.
        # From 1.5 r_H out: nearer the horizon e^(F_1 - F_0)/N passes 60, and the tracing's own step in r
        # moves r_* by more than the tolerance.
        self.check_null("herdeiro_radu", "axis", self.hair.tortoise_fast, self.hair.rH * 1.5)

    def test_every_ray_on_the_equator_is_null(self):
        # With varphi divided out the plane's g_tt is -e^(2F_0) N, so c dt = +-e^(F_1 - F_0) dr/N there too,
        # with the functions of the equator, integrated between the points of each ray.
        hair = self.hair

        def rate(x):
            return np.exp(hair("F_1", x, math.pi / 2) - hair("F_0", x, math.pi / 2)) / (1 - hair.rH / x)
        count = 0
        for r, ct in self.rays("herdeiro_radu", "equator"):
            keep = r > hair.rH * 1.05
            r, ct = r[keep], ct[keep]
            if len(r) < 2:
                continue
            for a, b, step in zip(r, r[1:], np.abs(np.diff(ct))):
                x = np.linspace(a, b, 21)
                f = rate(x)
                want = abs((x[1] - x[0]) / 3 * (f[0] + f[-1] + 4 * f[1:-1:2].sum() + 2 * f[2:-1:2].sum()))
                self.assertLess(abs(step - want), 0.02)
            count += 1
        self.assertGreater(count, 20)

    def test_every_ray_of_kerrs_axis_is_null(self):
        rk, ct = 0.98173, -0.44217
        M, Rp, Rm = rk / 2 - ct, rk - ct, -ct

        def star(r):
            return (r - ct) + 2 * M * Rp / rk * np.log((r - rk) / rk) - 2 * M * Rm / rk * np.log(r / rk)
        self.check_null("kerr_member", "axis", star, rk * 1.0001)

    def test_the_embedded_surface(self):
        surface = json.loads((DATA / "embedding" / "kerr_scalar_hair.json").read_text())["views"][0]["surfaces"][0]
        pieces = {p["id"]: np.array(p["points"]) for p in surface["pieces"]}
        r, rho, z = pieces["exterior"].T
        # rho = e^(F_2) r on the equator, the throat 0.687/mu round, and dz/dr = sqrt(g_rr - rho'^2).
        self.assertLess(float(np.max(np.abs(rho - self.hair.circumference_radius(r)))), 1e-6)
        self.assertAlmostEqual(float(rho[0]), 0.6872, places=3)
        for a, b, rise in list(zip(r, r[1:], np.diff(z)))[5::7]:
            x = np.linspace(a, b, 201)
            rate = np.exp(2 * self.hair("F_1", x, math.pi / 2)) / (1 - self.hair.rH / x)
            slope = np.sqrt(np.maximum(rate - np.gradient(self.hair.circumference_radius(x), x) ** 2, 0))
            h = x[1] - x[0]
            want = h / 3 * (slope[0] + slope[-1] + 4 * slope[1:-1:2].sum() + 2 * slope[2:-1:2].sum())
            self.assertAlmostEqual(float(rise), float(want), delta=1e-4 * (1 + abs(want)))
        # The other exterior is the mirror, and Kerr's surface of the same mass and angular momentum
        # has its throat at 2GM/c^2 = 1.866/mu, wider than the hairy hole's.
        other = pieces["other_exterior"]
        self.assertLess(float(np.max(np.abs(other[:, 1] - self.hair.circumference_radius(other[:, 0])))), 1e-6)
        self.assertLess(float(np.max(np.abs(other[:, 2] + np.interp(other[:, 0], r, z)))), 1e-3)
        self.assertAlmostEqual(float(pieces["kerr"][0][1]), 1.866, places=2)

    def test_the_conformal_diagram_puts_the_axis_in_a_diamond(self):
        for view in json.loads((DATA / "conformal" / "kerr_scalar_hair.json").read_text())["views"]:
            for layer in view["layers"]:
                if layer["kind"] == "line" and "points" in layer:
                    points = np.array(layer["points"])
                    self.assertLess(float(np.max(np.abs(points[:, 0]) + np.abs(points[:, 1]))), np.pi + 1e-3)


if __name__ == "__main__":
    unittest.main()
