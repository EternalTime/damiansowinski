"""Tilted universes and the whimper: python3 -m unittest discover -s _tools

What the texts and drawings of `tilted_universes` state, held on the published files: Farnsworth's
chart is dust of the density its convention gives, its singularity and its Cauchy horizon stand
where the drawings and captions put them, the flat model is flat and is the inertial chart carried
along the map of its convention, its lines of matter are the straight lines the inertial drawing
draws, and the conformal diagram's corners and edges are where its construction derives them. The
tests that read the published components need sympy and are skipped where it is absent.
"""
import importlib.util
import json
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
HAS_SYMPY = importlib.util.find_spec("sympy") is not None
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))

ETA_S, ETA_H = 1.96839702026073, 2.37295716191676     # X = 0 and X = C at C = W = 1
U_S, U_H = 0.770679156116326, 1.47248850095861


def X_of(eta):
    return math.sinh(eta / 2) ** 2 - 1 / math.tanh(eta / 2)


def u_of(eta):
    return (math.sinh(eta) - eta) / 2


def published(system):
    metric = json.loads((DATA / "metrics" / "tilted_universes.json").read_text(encoding="utf-8"))
    return next(c for c in metric["coordinates"] if c["id"] == system)


def drawn(kind):
    return json.loads((DATA / kind / "tilted_universes.json").read_text(encoding="utf-8"))


class Numbers(unittest.TestCase):
    """The numbers the drawings and captions state, from the solution itself."""

    def test_the_singularity_and_the_horizon(self):
        self.assertAlmostEqual(X_of(ETA_S), 0.0, places=12)
        self.assertAlmostEqual(X_of(ETA_H), 1.0, places=12)
        self.assertAlmostEqual(u_of(ETA_S), U_S, places=12)
        self.assertAlmostEqual(u_of(ETA_H), U_H, places=12)
        self.assertEqual((round(ETA_S, 2), round(ETA_H, 2), round(U_S, 2), round(U_H, 2)), (1.97, 2.37, 0.77, 1.47))
        # The time each world line of the dust spends under the horizon, in W/c.
        self.assertEqual(round(U_H - U_S, 2), 0.70)

    def test_the_scale_factor_at_the_singularity(self):
        """X = 0 is Y = coth(eta/2), which is Y^3 = Y + 1 for Y = sinh^2(eta/2)."""
        Y = math.sinh(ETA_S / 2) ** 2
        self.assertAlmostEqual(Y ** 3, Y + 1, places=10)

    def test_the_tilt_at_the_top_of_the_homogeneous_view(self):
        """The dust crosses the surfaces of homogeneity at c C/X: 0.20 c at u = 5 W."""
        lo, hi = ETA_H, 10.0
        for _ in range(80):
            mid = (lo + hi) / 2
            lo, hi = (mid, hi) if u_of(mid) < 5 else (lo, mid)
        self.assertEqual(round(1 / X_of(lo), 2), 0.20)
        self.assertIn("$0.20\\,c$ at the top edge", drawn("diagrams")["systems"]["homogeneous"][0]["caption"][1])

    def test_the_density_is_positive_only_after_the_singularity(self):
        self.assertLess(X_of(ETA_S - 0.1), 0)
        self.assertGreater(X_of(ETA_S + 0.1), 0)
        self.assertTrue(all(X_of(ETA_S + 0.01 * k) < X_of(ETA_S + 0.01 * (k + 1)) for k in range(500)))


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class Theory(unittest.TestCase):
    """The published components, read as the checker reads them, with c = 1."""

    def read(self, system):
        import sympy as sp
        import verify_metrics as vm
        chart = published(system)
        reader = vm.Reader(chart["coords"], [p["symbol"] for p in chart["parameters"]], ())

        def tensor(field, variant=None):
            block = chart[field] if variant is None else chart[field]["variants"][variant]["nonzero"]
            out = sp.zeros(4, 4)
            for entry in block:
                i, j = (chart["coords"].index(x) for x in entry["indices"])
                out[i, j] = reader(entry["value"]).subs(reader.c, 1)
            return out
        return chart, reader, tensor

    def test_farnsworths_chart_is_dust_at_rest(self):
        """G_ab = rho u_a u_b with u_a = (-Y, C, 0, 0) and rho = 3/(W X sinh^4(eta/2)), at a point."""
        import sympy as sp
        chart, reader, tensor = self.read("farnsworth")
        eta, W, C = reader.symbol["\\eta"], reader.parameters["W"], reader.parameters["C"]
        at = {eta: sp.Rational(27, 10), W: sp.Rational(3, 2), C: sp.Rational(1, 2), reader.symbol["r"]: sp.Rational(1, 3)}
        Y = W * sp.sinh(eta / 2) ** 2
        X = Y - C * sp.cosh(eta / 2) / sp.sinh(eta / 2)
        rho = 3 / (W * X * sp.sinh(eta / 2) ** 4)
        u = [-Y, C, 0, 0]
        G = tensor("einstein_tensor", "ll")
        for i in range(4):
            for j in range(4):
                self.assertAlmostEqual(float((G[i, j] - rho * u[i] * u[j]).subs(at)), 0.0, places=10, msg=f"{i}{j}")
        scalar = reader(chart["ricci_scalar"].partition("=")[2]).subs(reader.c, 1)
        self.assertAlmostEqual(float((scalar - rho).subs(at)), 0.0, places=10)

    def test_the_horizon_is_null_and_the_curvature_finite_on_it(self):
        import sympy as sp
        chart, reader, tensor = self.read("farnsworth")
        eta, W, C = reader.symbol["\\eta"], reader.parameters["W"], reader.parameters["C"]
        at = {eta: sp.Float(ETA_H, 30), W: 1, C: 1, reader.symbol["r"]: 0}
        # g_rr = X^2 - C^2 and g^{eta eta} vanish together: the surface of constant eta is null.
        self.assertAlmostEqual(float(tensor("metric_components")[1, 1].subs(at)), 0.0, places=9)
        self.assertAlmostEqual(float(tensor("inverse_metric_components")[0, 0].subs(at)), 0.0, places=9)
        K = reader(chart["kretschmann"].partition("=")[2]).subs(reader.c, 1)
        self.assertLess(abs(float(K.subs(at))), 10)
        near, nearer = (abs(float(K.subs({**at, eta: sp.Float(ETA_S + d, 30)}))) for d in (1e-4, 1e-5))
        self.assertGreater(nearer, 50 * near)

    def test_the_flat_model_is_flat_and_is_the_inertial_chart(self):
        import sympy as sp
        chart, reader, tensor = self.read("flat_model")
        self.assertEqual(chart["riemann"]["variants"]["llll"]["nonzero"], [])
        u, r, y, z = (reader.symbol[n] for n in chart["coords"])
        C = reader.parameters["C"]
        plus = (u + C) * sp.exp(-r)
        minus = (u - C) * sp.exp(r) + (u + C) * sp.exp(-r) * (y ** 2 + z ** 2)
        image = [(plus + minus) / 2, (plus - minus) / 2, plus * y, plus * z]
        J = sp.Matrix(4, 4, lambda i, j: sp.diff(image[i], [u, r, y, z][j]))
        pulled = (J.T * sp.diag(-1, 1, 1, 1) * J - tensor("metric_components")).applyfunc(sp.simplify)
        self.assertEqual(pulled, sp.zeros(4, 4))
        interval = sp.simplify(-image[0] ** 2 + image[1] ** 2 + image[2] ** 2 + image[3] ** 2)
        self.assertEqual(sp.simplify(interval - (C ** 2 - u ** 2)), 0)
        _, _, flat = self.read("inertial")
        self.assertEqual(flat("metric_components"), sp.diag(-1, 1, 1, 1))


class Drawings(unittest.TestCase):
    def view(self, system):
        return drawn("diagrams")["systems"][system][0]

    def marker(self, view, kind):
        return [m for m in view["markers"] if m["kind"] == kind]

    def test_the_horizon_is_drawn_where_X_is_C(self):
        for system, level in (("farnsworth", ETA_H), ("homogeneous", U_H), ("flat_model", 1.0)):
            view = self.view(system)
            lo, hi = view["box"][2], view["box"][3]
            line = self.marker(view, "event")[0]["lines"][0]
            for _, y in line:
                self.assertAlmostEqual(y, (level - lo) / (hi - lo), places=3, msg=system)

    def test_the_singularity_is_drawn_where_X_is_zero(self):
        view = self.view("farnsworth")
        lo, hi = view["box"][2], view["box"][3]
        for _, y in self.marker(view, "singular")[0]["lines"][0]:
            self.assertAlmostEqual(y, (ETA_S - lo) / (hi - lo), places=3)

    def test_every_view_is_near_square(self):
        for views in drawn("diagrams")["systems"].values():
            for view in views:
                x0, x1, y0, y1 = view["box"]
                self.assertEqual(x1 - x0, y1 - y0, view["id"])

    def test_the_lines_of_matter_are_straight_and_pass_the_origin_at_C(self):
        """Each is x = sech(r) - T tanh(r): its distance from the origin, taken in the metric, is C = 1."""
        view = self.view("inertial")
        x0, x1, y0, y1 = view["box"]
        worlds = self.marker(view, "world")
        self.assertEqual(len(worlds), 7)
        for mark in worlds:
            (a, b), (c, d) = mark["lines"][0][0], mark["lines"][0][-1]
            xa, Ta, xb, Tb = x0 + a * (x1 - x0), y0 + b * (y1 - y0), x0 + c * (x1 - x0), y0 + d * (y1 - y0)
            slope = (xb - xa) / (Tb - Ta)                 # dx/dT = -tanh(r), below the speed of light
            self.assertLess(abs(slope), 1)
            at_rest = xa - slope * Ta                     # x at T = 0, which is sech(r)
            self.assertAlmostEqual(at_rest ** 2, 1 - slope ** 2, places=2)
            for p, q in mark["lines"][0]:
                x, T = x0 + p * (x1 - x0), y0 + q * (y1 - y0)
                self.assertAlmostEqual(x, at_rest + slope * T, places=2)
                self.assertGreater(T + x, 0)
                self.assertLess(x * x - T * T, 1 + 2e-2)

    def test_the_surfaces_of_homogeneity_are_at_a_constant_interval(self):
        view = self.view("inertial")
        x0, x1, y0, y1 = view["box"]
        for mark in self.marker(view, "surface") + self.marker(view, "shell"):
            intervals = [(x0 + p * (x1 - x0)) ** 2 - (y0 + q * (y1 - y0)) ** 2 for p, q in mark["lines"][0]]
            self.assertLess(max(intervals) - min(intervals), 0.05)

    def test_the_conformal_diagram(self):
        half = math.pi / 2
        views = {v["id"]: v for v in drawn("conformal")["views"]}
        self.assertEqual(set(views), {"farnsworth", "homogeneous", "flat_model", "inertial"})
        corners = [[half, -half], [-half, half], [0, math.pi], [half, half]]
        for vid, view in views.items():
            region = next(layer for layer in view["layers"] if layer["class"] == "region")
            for got, want in zip(region["points"], corners, strict=True):
                self.assertAlmostEqual(got[0], want[0], places=3, msg=vid)
                self.assertAlmostEqual(got[1], want[1], places=3, msg=vid)
            horizon = next(layer for layer in view["layers"] if layer["class"] == "horizon")
            self.assertEqual(horizon["points"][0], [0.0, 0.0])
            # Every line drawn lies in the quadrilateral: left of X = pi/2, above the lower left edge
            # T = -X, and below the two upper edges.
            for layer in view["layers"]:
                if layer["kind"] != "line":
                    continue
                for X, T in layer["points"]:
                    self.assertLessEqual(X, half + 1e-3, vid)
                    self.assertGreaterEqual(T + X, -1e-3, vid)
                    self.assertLessEqual(T - X, math.pi + 1e-3, vid)
                    self.assertLessEqual(T + X, math.pi + 1e-3, vid)
        # The whimper and the singularity are drawn singular for the dust, and as plain edges of the flat model.
        for vid in ("farnsworth", "homogeneous"):
            self.assertEqual(sum(layer["kind"] == "zig" for layer in views[vid]["layers"]), 2)
        for vid in ("flat_model", "inertial"):
            self.assertEqual(sum(layer["kind"] == "zig" for layer in views[vid]["layers"]), 0)

    def test_the_embedded_surface_is_a_pseudosphere_of_the_radius_of_curvature(self):
        """At eta = 3 the circle at r has the radius Y e^(-r), from the rim where that is sqrt(X^2 - C^2)."""
        Y = math.sinh(1.5) ** 2
        ell = math.sqrt(X_of(3.0) ** 2 - 1)
        self.assertEqual((round(Y, 2), round(X_of(3.0), 2), round(ell, 2), round(math.log(Y / ell), 2)),
                         (4.53, 3.43, 3.28, 0.32))
        view = drawn("embedding")["views"][0]
        piece = view["surfaces"][0]["pieces"][0]
        r0, rho0, z0 = piece["points"][0]
        self.assertAlmostEqual(r0, math.log(Y / ell), places=3)
        self.assertAlmostEqual(rho0, ell, places=3)
        for r, rho, z in piece["points"]:
            self.assertAlmostEqual(rho, Y * math.exp(-r), places=3)
            s = math.sqrt(max(0.0, 1 - (rho / ell) ** 2))
            self.assertAlmostEqual(z - z0, ell * (math.atanh(s) - s), places=2)

    def test_the_null_coordinates_are_regular_on_the_horizon(self):
        """V of tilted_null_tables is continuous through eta_H and vanishes there, and -V U = 1 on X = 0."""
        if importlib.util.find_spec("numpy") is None or importlib.util.find_spec("contourpy") is None:
            self.skipTest("numpy and contourpy are not installed")
        import numpy as np
        import null_rays as nr
        t = nr.tilted_null_tables()
        self.assertAlmostEqual(t["eta_s"], ETA_S, places=12)

        def V(eta, r):
            return (eta - t["eta_H"]) / (t["eta_H"] - t["eta_s"]) * math.exp((r - np.interp(eta, t["etas"], t["G"])) / t["a"])
        self.assertEqual(V(t["eta_H"], 0.3), 0.0)
        self.assertLess(abs(V(t["eta_H"] + 1e-6, 0.3) + V(t["eta_H"] - 1e-6, 0.3)), 1e-9)
        U = math.exp(-(0.3 - 0.0) / t["a"])
        self.assertAlmostEqual(-V(t["eta_s"], 0.3) * U, 1.0, places=10)
        # G is smooth across the horizon: its slope changes by little from one side to the other.
        k = int(np.searchsorted(t["etas"], t["eta_H"]))
        slopes = np.diff(t["G"][k - 50:k + 50]) / np.diff(t["etas"][k - 50:k + 50])
        self.assertLess(float(np.ptp(slopes)), 1e-2 * abs(float(slopes.mean())) + 1e-2)


if __name__ == "__main__":
    unittest.main()
