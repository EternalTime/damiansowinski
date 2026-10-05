"""The maximally supersymmetric plane wave of type IIB supergravity:
.venv.noindex/bin/python -m unittest discover -s _tools

What its texts state, held on the published files: every chart has ten dimensions; the Brinkmann
chart's one Ricci component is R_uu = 8 mu^2, which Metsaev's five-form F_u1234 = F_u5678 = 2 mu
supplies in his normalization R_MN = F_MPQRS F_N^PQRS / 24; every chart is conformally flat, with
vanishing Ricci and Kretschmann scalars; the Rosen and conformally flat charts are the Brinkmann
chart carried along x_i = y_i sin(mu c u) and mu c U = tan(mu c u); Berenstein and Nastase's chain
of maps carries the wave onto the Einstein static universe with the conformal factor
1/|e^(i psi) - cos(alpha) e^(i beta)|^2, and not a quarter of it; and the drawings hold what their
captions say. The tests of the published mathematics need sympy and mpmath and are skipped where
they are absent; the tests of the drawings read the files alone.
"""
import importlib.util
import json
import math
import random
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
HAS_SYMPY = importlib.util.find_spec("sympy") is not None and importlib.util.find_spec("mpmath") is not None
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))
ID = "maximally_supersymmetric_plane_wave"
CHARTS = ["brinkmann", "rosen", "conformally_flat"]


def published():
    return json.loads((DATA / "metrics" / f"{ID}.json").read_text(encoding="utf-8"))


def drawing(kind):
    return json.loads((DATA / kind / f"{ID}.json").read_text(encoding="utf-8"))


def chart(system):
    return next(c for c in published()["coordinates"] if c["id"] == system)


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class Published(unittest.TestCase):
    def setUp(self):
        import sympy
        import verify_metrics
        self.sp, self.vm = sympy, verify_metrics

    def metric(self, system):
        """The published metric of a chart as a sympy matrix with c = 1, and its reader."""
        entry = chart(system)
        reader = self.vm.Reader(entry["coords"], ["\\mu"], ())
        n = len(entry["coords"])
        g = self.sp.zeros(n, n)
        for e in entry["metric_components"]:
            i, j = (entry["coords"].index(k) for k in e["indices"])
            g[i, j] = reader(e["value"]).subs(reader.c, 1)
        return g, reader, entry

    def test_every_chart_has_ten_dimensions_and_no_scalar_curvature(self):
        for system in CHARTS:
            entry = chart(system)
            self.assertEqual(len(entry["coords"]), 10, system)
            self.assertEqual(entry["ricci_scalar"], "R = 0", system)
            self.assertEqual(entry["kretschmann"], "K = 0", system)
            for variant in entry["weyl_tensor"]["variants"].values():
                self.assertEqual(variant["nonzero"], [], f"{system}: the Weyl tensor does not vanish")

    def test_metsaev_s_five_form_supplies_the_one_ricci_component(self):
        # F_u PQRS F_u^PQRS sums 4! orderings of each of the two blocks of transverse indices, raised
        # with the identity, so the source is 24 (4 mu^2 + 4 mu^2)/24 = 8 mu^2.
        entry = chart("brinkmann")
        reader = self.vm.Reader(entry["coords"], ["\\mu"], ())
        mu = reader.parameters["mu"]
        source = math.factorial(4) * 2 * (2 * mu) ** 2 / 24
        ricci = entry["ricci_tensor"]["variants"]["ll"]["nonzero"]
        self.assertEqual([e["indices"] for e in ricci], [["u", "u"]])
        self.assertEqual(self.sp.simplify(reader(ricci[0]["value"]) - source), 0)
        self.assertIn("2\\mu c\\,du\\wedge", entry["parameters"][0]["description"])

    def test_the_rosen_and_conformally_flat_charts_are_the_brinkmann_chart_carried_over(self):
        sp = self.sp
        gB, rB, _ = self.metric("brinkmann")
        for system in ("rosen", "conformally_flat"):
            g, r, entry = self.metric(system)
            x = [r.symbol[c] for c in entry["coords"]]
            mu = r.parameters["mu"]
            u, v, across = x[0], x[1], x[2:]
            if system == "rosen":
                image = [u, v + mu * sum(y ** 2 for y in across) * sp.sin(mu * u) * sp.cos(mu * u) / 2] \
                    + [y * sp.sin(mu * u) for y in across]
            else:
                image = [sp.atan(mu * u) / mu,
                         v - mu ** 2 * u * sum(X ** 2 for X in across) / (2 * (1 + mu ** 2 * u ** 2))] \
                    + [X / sp.sqrt(1 + mu ** 2 * u ** 2) for X in across]
            at = {rB.symbol[c]: image[k] for k, c in enumerate(chart("brinkmann")["coords"])}
            at[rB.parameters["mu"]] = mu
            J = sp.Matrix(10, 10, lambda i, j: sp.diff(image[i], x[j]))
            pulled = J.T * gB.subs(at, simultaneous=True) * J
            for i in range(10):
                for j in range(i, 10):
                    self.assertEqual(sp.simplify(pulled[i, j] - g[i, j]), 0, f"{system}: slot {i}{j}")

    def test_berenstein_and_nastase_s_einstein_static_universe_has_the_factor_one(self):
        """Their chain from the Einstein static universe, (psi, alpha, beta) and the 7-sphere, back to
        the Brinkmann chart at mu = 1, followed numerically: the pulled back metric is
        1/|e^(i psi) - cos(alpha) e^(i beta)|^2 times -dpsi^2 + dalpha^2 + cos^2(alpha) dbeta^2 +
        sin^2(alpha) dOmega_7^2, on the three coordinates and on the 7-sphere."""
        import mpmath as mp
        mp.mp.dps = 40

        def image(p):
            psi, alpha, beta = p
            zeta = mp.acos(-mp.cos(alpha) * mp.cos(beta))
            sz = mp.sin(zeta)
            theta = mp.atan2(mp.sin(alpha) / sz, -mp.cos(alpha) * mp.sin(beta) / sz)
            up, vp = mp.tan((psi + zeta) / 2), -mp.tan((psi - zeta) / 2)
            r, tau = (up + vp) / 2, (up - vp) / 2
            sigma, xr = r * mp.cos(theta), r * mp.sin(theta)
            U, Xm = sigma + tau, (sigma - tau) / 2
            xplus = mp.atan(U)
            rho = xr * mp.cos(xplus)
            # Their metric is 2 dx^+ dx^- - x^2 (dx^+)^2 + dx^2; the published one has u = x^+, v = -x^-.
            return [xplus, -(Xm + rho ** 2 * U / 2), rho]

        rng = random.Random(4)
        for _ in range(4):
            p = [mp.mpf(rng.uniform(-0.6, 0.6)), mp.mpf(rng.uniform(0.2, 1.3)), mp.mpf(rng.uniform(0.3, 2.8))]
            u, v, rho = image(p)
            J = mp.matrix(3, 3)
            for i in range(3):
                for j in range(3):
                    J[i, j] = mp.diff(lambda s, i=i, j=j: image([p[k] if k != j else s for k in range(3)])[i], p[j])
            g = mp.matrix([[-rho ** 2, -1, 0], [-1, 0, 0], [0, 0, 1]])
            pulled = J.T * g * J
            psi, alpha, beta = p
            factor = 1 / abs(mp.exp(1j * psi) - mp.cos(alpha) * mp.exp(1j * beta)) ** 2
            want = [-factor, factor, factor * mp.cos(alpha) ** 2]
            for i in range(3):
                for j in range(3):
                    self.assertLess(abs(pulled[i, j] - (want[i] if i == j else 0)), mp.mpf(10) ** -30)
            # The 7-sphere: x = rho times a unit vector, so its coefficient is rho^2.
            self.assertLess(abs(rho ** 2 - factor * mp.sin(alpha) ** 2), mp.mpf(10) ** -30)
            self.assertGreater(abs(factor / 4 - factor), 0.1 * factor)


class Drawings(unittest.TestCase):
    def test_the_rays_of_the_brinkmann_chart_meet_again_at_the_foci(self):
        view = drawing("diagrams")["projections"]["brinkmann"][0]
        above = [layer["points"] for layer in view["layers"] if layer["class"] == "above"]
        below = [layer["points"] for layer in view["layers"] if layer["class"] == "below"]
        self.assertEqual((len(above), len(below)), (6, 6))
        for ray in above:
            amplitude = max(abs(x) for x, _ in ray)
            for x, t in ray:
                self.assertAlmostEqual(abs(x), abs(amplitude * math.sin(t)), delta=2e-3)
            self.assertAlmostEqual(ray[-1][1], 2 * math.pi, delta=1e-3)
        for ray in below:
            a = ray[0][0]
            for x, t in ray:
                self.assertAlmostEqual(x, a * math.cos(t), delta=2e-3)
        foci = sorted(layer["at"][1] for layer in view["layers"] if layer["kind"] == "point")
        self.assertEqual([round(t, 3) for t in foci], [round(math.pi, 3), round(2 * math.pi, 3)])

    def test_the_rosen_rays_from_the_origin_are_upright_and_the_others_cotangents(self):
        view = drawing("diagrams")["projections"]["rosen"][0]
        for layer in view["layers"]:
            if layer["class"] == "above":
                self.assertEqual(len({x for x, _ in layer["points"]}), 1)
            if layer["class"] == "below":
                # Each crosses the axis once, at mu c u = pi/2, where the light runs across it.
                pts = layer["points"]
                crossings = [(a, b) for a, b in zip(pts, pts[1:]) if (a[0] > 0) != (b[0] > 0)]
                self.assertEqual(len(crossings), 1)
                (xa, ta), (xb, tb) = crossings[0]
                self.assertAlmostEqual(ta + (tb - ta) * xa / (xa - xb), math.pi / 2, delta=1e-2)
                for x, t in layer["points"]:
                    self.assertTrue(0 < t < math.pi)

    def test_light_in_the_conformally_flat_chart_runs_straight(self):
        view = drawing("diagrams")["projections"]["conformally_flat"][0]
        for layer in view["layers"]:
            if layer["class"] in ("above", "below"):
                self.assertEqual(len(layer["points"]), 2, "a ray of the conformally flat chart is not one segment")

    def test_the_conformal_boundary_is_one_null_line_and_the_bands_are_where_the_charts_reach(self):
        views = {v["id"]: v for v in drawing("conformal")["views"]}
        self.assertEqual(set(views), set(CHARTS))
        for system, v in views.items():
            scri = [layer["points"] for layer in v["layers"] if layer["class"] == "scri"]
            pi2 = round(2 * math.pi, 4)
            self.assertEqual(scri, [[[0, 0], [pi2, pi2]], [[round(math.pi, 4), -round(math.pi, 4)], [pi2, 0]]],
                             system)
            area = 0.0
            for layer in v["layers"]:
                if layer["class"] == "cover":
                    P = layer["points"]
                    area += abs(sum(P[i][0] * P[i - 1][1] - P[i - 1][0] * P[i][1] for i in range(len(P)))) / 2
            # The whole window, one turn by three half turns; each band, between two lines of constant u
            # pi apart in mu c u, runs from one turn of the boundary to the next, an area of 2 pi^2.
            want = {"brinkmann": 6 * math.pi ** 2, "rosen": 2 * math.pi ** 2, "conformally_flat": 2 * math.pi ** 2}[system]
            self.assertAlmostEqual(area, want, delta=1e-2, msg=system)

    def test_the_ring_closes_on_the_axis_at_a_quarter_period(self):
        views = {v["id"]: v for v in drawing("embedding")["views"]}
        self.assertEqual(set(views), {"tube", "ring"})
        caption = " ".join(views["ring"]["caption"])
        self.assertIn("\\cos(\\mu cu)/\\mu", caption)
        self.assertIn("\\mu cu = \\pi/2", caption)


if __name__ == "__main__":
    unittest.main()
