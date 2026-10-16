"""Kodama's kink as _tools/derivations/kodama_kink.py solves it, held to the physics: the one field
equation its construction does not use, Kodama's own figures, the values at the throat his
CBPF-NF-003/78 derives, the thin wall the kink becomes as f grows, its parity, its vacuum outside,
and the areal chart's own Einstein equation. Then the files its drawings wrote, from their numbers
alone: every ray drawn is a null curve of the kink's metric, the embedded surface climbs at
sqrt(1 - r'^2), and the conformal diagram lies in the diamond. It needs sympy and scipy, so it
runs under .venv.noindex/bin/python and is skipped under a Python without them."""
import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))
try:
    import numpy as np
    import sympy as sp
    import kodama_kink as kk
except ImportError:
    kk = None

DATA = ROOT / "MFS" / "assets" / "data"


def tortoise(kink, R):
    """r* = int_0^R e^{-h} dR, by Simpson's rule on the solver's own h, odd about the throat."""
    out = []
    for top in np.atleast_1d(R):
        x = np.linspace(0.0, abs(float(top)), 4001)
        f = np.exp(-kink.values(x)["h"][0])
        h = x[1] - x[0]
        out.append(np.sign(top) * h / 3 * (f[0] + f[-1] + 4 * f[1:-1:2].sum() + 2 * f[2:-1:2].sum()))
    return np.array(out)


@unittest.skipIf(kk is None, "needs sympy and scipy")
class Kink(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.kink = kk.kink()

    def test_the_kink_solves_the_field_equation_its_construction_did_not_use(self):
        # G^R_R is the Bianchi identity's consequence of the three equations integrated.
        R = np.linspace(0.0, 0.99 * self.kink.edge, 2000)
        self.assertLess(float(np.max(np.abs(self.kink.radial(R)))), 1e-10)

    def test_it_is_the_kink_the_drawings_state(self):
        self.assertEqual(self.kink.stated(), {})
        self.assertGreater(self.kink.m, 0)
        self.assertGreater(self.kink.lapse, 0)

    def test_kodamas_figures(self):
        # CBPF-NF-003/78, Fig. 1: the throats of f = 0.77, 1 and 1.25 at mu r_0 about 0.25, 0.78 and
        # 1.47; CBPF-NF-036/78, Fig. 1c: h about -0.36 at the throat of f = 1.25.
        for f, read in ((0.77, 0.25), (1.0, 0.78), (1.25, 1.47)):
            self.assertAlmostEqual(kk.kink(f).r0, read, delta=0.03)
        self.assertAlmostEqual(float(self.kink.values([0.0])["h"][0][0]), -0.36, delta=0.03)

    def test_the_values_at_the_throat(self):
        # Kodama's (25) to (28) with epsilon = -1: S'(0)^2 = (1/2 + f^2/r_0^2)/f^2 and
        # r''(0) = E(0)/2 = (2 f^2 + r_0^2)/(2 r_0 f^2), with r'(0) = h'(0) = S(0) = 0.
        k, f = self.kink, self.kink.f
        r, r1, r2, h, h1, h2, S, S1 = k.inside(0.0)
        self.assertAlmostEqual(float(S1[0]) ** 2, (0.5 + f ** 2 / k.r0 ** 2) / f ** 2, places=12)
        self.assertAlmostEqual(float(r2[0]), (2 * f ** 2 + k.r0 ** 2) / (2 * k.r0 * f ** 2), places=10)
        self.assertEqual((float(r1[0]), float(h1[0]), float(S[0])), (0.0, 0.0, 0.0))

    def test_the_heavy_kink_is_a_thin_wall_on_the_photon_sphere(self):
        # For large f the kink is a wall of width 1/mu on a throat of radius about f^2/mu. Israel's
        # junction for a wall with p = -sigma between two Schwarzschild exteriors puts it at a = 3m,
        # and its tension, the flat kink's 4 sqrt(2)/3 over f^2, then gives a/f^2 = sqrt(6)/2.
        heavy = kk.kink(20.0)
        self.assertAlmostEqual(heavy.r0 / heavy.m, 3.0, delta=0.005)
        self.assertAlmostEqual(heavy.r0 / 400.0, np.sqrt(6) / 2, delta=0.002)

    def test_the_kink_closes_toward_kodamas_threshold(self):
        # Below f about 0.645 Kodama found no kink: the throat shrinks toward zero as f falls.
        throats = [kk.kink(f).r0 for f in (0.7, 0.66, 0.645)]
        self.assertTrue(throats[0] > throats[1] > throats[2] > 0)
        self.assertLess(throats[2], 0.03)

    def test_parity(self):
        R = np.linspace(0.1, 20, 50)
        plus, minus = self.kink.values(R), self.kink.values(-R)
        for name in ("h", "r"):
            for order, sign in ((0, 1), (1, -1), (2, 1)):
                self.assertTrue(np.allclose(minus[name][order], sign * plus[name][order], atol=1e-12))
        self.assertTrue(np.allclose(self.kink.field(-R)[0], -self.kink.field(R)[0], atol=1e-12))

    def test_the_vacuum_outside(self):
        # The field reaches its vacuum 1/f, and from there on the metric is Schwarzschild's of the
        # kink's mass: r'^2 = 1 - 2m/r = e^{2h}, smoothly across the edge.
        k = self.kink
        self.assertLess(abs(k.left), 1e-8)
        R = np.array([k.edge * 0.999, k.edge * 1.001, 30.0, 100.0])
        v = k.values(R)
        g = 1 - 2 * k.m / v["r"][0]
        self.assertTrue(np.allclose(v["r"][1] ** 2, g, atol=1e-7))
        self.assertTrue(np.allclose(np.exp(2 * v["h"][0]), g, atol=1e-7))

    def test_the_areal_chart_solves_its_own_einstein_equation(self):
        # The published G^t_t of the areal chart, in eta and alpha read at the same sphere, is the
        # ghost field's energy, e^{-2 alpha} (dS/dr)^2 + V = (dS/dR)^2 + V.
        entry = next(c for c in json.loads((DATA / "metrics" / "kodama_kink.json").read_text())["coordinates"]
                     if c["id"] == "areal")
        sys.path.insert(0, str(ROOT / "_tools" / "derivations"))
        import verify_metrics as vm
        reader = vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
        r = reader.symbol["r"]
        alpha = reader.parameters["alpha"]
        A0, A1 = sp.symbols("A0 A1")
        ul = entry["einstein_tensor"]["variants"]["ul"]["nonzero"]
        tt = reader(next(c["value"] for c in ul if c["indices"] == ["t", "t"]))
        tt = sp.lambdify((r, A0, A1), tt.subs(sp.Derivative(alpha, r), A1).subs(alpha, A0), "numpy")
        areal = self.kink.areal()
        R = np.linspace(0.3, 10, 60)
        x = self.kink.values(R)["r"][0]
        a = areal.values(x)["alpha"]
        S, S1 = self.kink.field(R)
        f = self.kink.f
        energy = S1 ** 2 + (1 - f ** 2 * S ** 2) ** 2 / (2 * f ** 2)
        self.assertLess(float(np.max(np.abs(tt(x, a[0], a[1]) - energy))), 1e-6)


@unittest.skipIf(kk is None, "needs sympy and scipy")
class Drawings(unittest.TestCase):
    """The files the generators wrote, measured against the solver from their numbers alone."""

    @classmethod
    def setUpClass(cls):
        cls.kink = kk.kink()

    def rays(self, system):
        view = json.loads((DATA / "diagrams" / "kodama_kink.json").read_text())["systems"][system][0]
        x_lo, x_hi, y_lo, y_hi = view["box"]
        for family in ("P", "M"):
            for line in view["rays"][family]:
                points = np.array(line)
                yield x_lo + points[:, 0] * (x_hi - x_lo), y_lo + points[:, 1] * (y_hi - y_lo)

    def test_every_ray_of_the_proper_radial_chart_is_null(self):
        # c dt = +-e^{-h} dR, so along a ray |ct| moves by the tortoise coordinate's change.
        count = 0
        for R, ct in self.rays("proper_radial"):
            if len(R) < 2:
                continue
            rs = tortoise(self.kink, R)
            self.assertLess(float(np.max(np.abs(np.abs(np.diff(ct)) - np.abs(np.diff(rs))))), 0.02)
            count += 1
        self.assertGreater(count, 20)

    def test_every_ray_of_the_areal_chart_is_null(self):
        areal = self.kink.areal()
        count = 0
        for r, ct in self.rays("areal"):
            keep = r > self.kink.r0 * 1.001
            if keep.sum() < 2:
                continue
            rs = tortoise(self.kink, areal.R_of(r[keep]))
            self.assertLess(float(np.max(np.abs(np.abs(np.diff(ct[keep])) - np.abs(np.diff(rs))))), 0.02)
            count += 1
        self.assertGreater(count, 10)

    def test_the_embedded_surface_climbs_at_the_metrics_rate(self):
        surface = json.loads((DATA / "embedding" / "kodama_kink.json").read_text())["views"][0]["surfaces"][0]
        pieces = {p["id"]: np.array(p["points"]) for p in surface["pieces"]}
        R, rho, z = pieces["whole"].T
        self.assertTrue(np.allclose(rho, self.kink.values(R)["r"][0], atol=1e-6))
        for a, b, rise in zip(R, R[1:], np.diff(z)):
            x = np.linspace(a, b, 201)
            slope = np.sqrt(np.maximum(1 - self.kink.values(x)["r"][1] ** 2, 0))
            h = x[1] - x[0]
            want = h / 3 * (slope[0] + slope[-1] + 4 * slope[1:-1:2].sum() + 2 * slope[2:-1:2].sum())
            self.assertAlmostEqual(float(rise), float(want), places=5)
        # The throat is the narrowest circle, of radius r_0, and Flamm's paraboloid of the same
        # mass, drawn beside it, starts lower down at its own throat r = 2m and meets it at the rim.
        self.assertAlmostEqual(float(rho.min()), self.kink.r0, places=6)
        self.assertAlmostEqual(float(pieces["vacuum"][0][0]), 2 * self.kink.m, places=9)
        self.assertLess(2 * self.kink.m, self.kink.r0)
        self.assertAlmostEqual(float(pieces["vacuum"][-1][2]), float(z[-1]), places=6)

    def test_the_conformal_diagram_lies_in_the_diamond(self):
        for view in json.loads((DATA / "conformal" / "kodama_kink.json").read_text())["views"]:
            for layer in view["layers"]:
                if layer["kind"] == "line" and "points" in layer:
                    X, T = np.array(layer["points"]).T
                    self.assertLess(float(np.max(np.abs(T) + np.abs(X))), np.pi + 1e-3)


if __name__ == "__main__":
    unittest.main()
