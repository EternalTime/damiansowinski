"""Sultana and Dyer's black hole, Schwarzschild's metric times the scale factor of a universe of dust:
python3 -m unittest discover -s _tools

What its texts and drawings state, held on the published files. The first class needs nothing
but the files: the event horizon and the two trapping horizons where every spacetime diagram
puts them, the big bang of the plane of Schwarzschild's time, the paraboloid each moment embeds
as, and the triangle of the conformal diagram. The second reads the published tensors through the
checker's Reader and holds them to the two fluids, to the sign of the dust's density and to a
horizon of finite curvature, as the History and sultana_dyer.md state; it needs sympy and is
skipped where it is absent.
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

ETA0 = 3.0      # eta_0 in r_s, as every drawing takes it


def load(folder):
    return json.loads((DATA / folder / "sultana_dyer.json").read_text(encoding="utf-8"))


def inner(eta):
    """The inner trapping horizon, r_s = 1."""
    return (math.sqrt(eta * eta + 12 * eta + 4) - eta - 2) / 4


def chart_point(view, at):
    """A point of a view's unit square in the drawn axes."""
    X0, X1, Y0, Y1 = view["box"]
    return X0 + at[0] * (X1 - X0), Y0 + at[1] * (Y1 - Y0)


class Drawings(unittest.TestCase):
    def view(self, system):
        return load("diagrams")["systems"][system][0]

    def marked(self, view, kind):
        return [chart_point(view, at) for m in view["markers"] if m["kind"] == kind
                for line in m.get("lines", []) for at in line]

    def test_the_event_horizon_is_marked_on_r_s_where_a_chart_crosses_it(self):
        for system in ("kerr_schild", "eddington_finkelstein_ingoing"):
            points = self.marked(self.view(system), "grr")
            self.assertGreater(len(points), 1, system)
            for r, _ in points:
                self.assertAlmostEqual(r, 1.0, delta=0.02, msg=system)

    def test_the_trapping_horizons_are_the_hubble_radius_and_the_curve_inside_r_s(self):
        """|grad R|^2 = 0 for R = eta^2 r/eta_0^2 on r = eta/2 and on
        r = (sqrt(eta^2 + 12 r_s eta + 4 r_s^2) - eta - 2 r_s)/4, in the Kerr-Schild plane and, with
        eta = v - r, in the plane of v and r; both branches are drawn."""
        for system, eta_of in (("kerr_schild", lambda y, r: y), ("eddington_finkelstein_ingoing", lambda y, r: y - r)):
            points = self.marked(self.view(system), "apparent")
            on = {"hubble": 0, "inner": 0}
            for r, y in points:
                eta = eta_of(y, r)
                miss = {"hubble": abs(r - eta / 2), "inner": abs(r - inner(eta))}
                branch = min(miss, key=miss.get)
                self.assertLess(miss[branch], 0.03, f"{system}: ({r}, {y})")
                on[branch] += 1
            self.assertGreater(on["hubble"], 1, system)
            self.assertGreater(on["inner"], 3, system)

    def test_the_big_bang_of_schwarzschilds_time_is_the_curve_of_no_conformal_time(self):
        """eta = ct + r_s ln(r/r_s - 1) = 0, drawn as a singular curve up to the top of the plane."""
        view = self.view("schwarzschild_time")
        points = self.marked(view, "singular")
        self.assertGreater(len(points), 10)
        for r, t in points:
            self.assertAlmostEqual(r, 1 + math.exp(-t), delta=0.03)
        self.assertAlmostEqual(max(t for _, t in points), view["box"][3], delta=0.05)

    def test_every_moment_is_a_paraboloid_magnified_by_the_scale_factor(self):
        """rho = a r and z = 2 a sqrt(r_s r) with a = eta^2/eta_0^2, so z^2 = 4 a r_s rho, with the
        event horizon and both trapping horizons among the circles marked."""
        view = load("embedding")["views"][0]
        frames = view["movie"]["frames"]
        self.assertGreater(len(frames), 20)
        self.assertEqual([s["time"] for s in view["surfaces"]], [3.0, 4.0, 5.0, 6.0])
        for frame in frames:
            eta = frame["value"]
            a = (eta / ETA0) ** 2
            piece, = frame["pieces"]
            for r, rho, z in piece["points"]:
                self.assertAlmostEqual(rho, a * r, delta=1e-5)
                self.assertAlmostEqual(z, 2 * a * math.sqrt(r), delta=1e-5)
            marked = sorted(ring["x"] for ring in frame["rings"] if ring["class"] == "horizon")
            want = sorted([1.0, inner(eta)] + ([eta / 2] if eta / 2 < 4 else []))
            self.assertEqual(len(marked), len(want))
            for got, radius in zip(marked, want):
                self.assertAlmostEqual(got, radius, places=9)

    def test_the_conformal_diagram_is_the_part_of_kruskals_after_the_big_bang(self):
        """With p = (T - X)/2 and q = (T + X)/2, tan q = V = e^((eta + r)/2) and
        tan p = U = (1 - r) e^((r - eta)/2) at r_s = 1. The big bang eta = 0 is U = (1 - 2 ln V) V,
        from the singularity at (0, pi/2) to i0 at (pi, 0); r = 0 lies on T = pi/2 as far as i+ at
        (pi/2, pi/2); and each moment of the embedding lies on its own eta."""
        for view in load("conformal")["views"]:
            zigs = [layer["points"] for layer in view["layers"] if layer["class"] == "singular"]
            self.assertEqual(len(zigs), 2, view["id"])
            bang = max(zigs, key=len)
            centre = min(zigs, key=len)
            for X, T in bang[1:-1]:
                p, q = (T - X) / 2, (T + X) / 2
                V = math.tan(q)
                self.assertAlmostEqual(math.tan(p), (1 - 2 * math.log(V)) * V, delta=2e-2 * max(1, V * V), msg=view["id"])
            ends = sorted([bang[0], bang[-1]])
            self.assertAlmostEqual(ends[0][0], 0, delta=1e-3)
            self.assertAlmostEqual(ends[0][1], math.pi / 2, delta=1e-3)
            self.assertAlmostEqual(ends[1][0], math.pi, delta=1e-3)
            self.assertAlmostEqual(ends[1][1], 0, delta=1e-3)
            for X, T in centre:
                self.assertAlmostEqual(T, math.pi / 2, places=4)
            self.assertAlmostEqual(max(X for X, _ in centre), math.pi / 2, places=4)
            self.assertEqual(len(view["slices"]), 4, view["id"])
            for mark in view["slices"]:
                eta = float(mark["label"].split("=")[1].split("\\")[0])
                for X, T in mark["lines"][0][1:-1]:
                    p, q = (T - X) / 2, (T + X) / 2
                    r = 2 * math.log(math.tan(q)) - eta
                    self.assertGreater(r, -1e-2)
                    self.assertAlmostEqual(math.tan(p), (1 - r) * math.exp((r - eta) / 2), delta=2e-3 * (1 + abs(r)))


@unittest.skipUnless(HAS_SYMPY, "needs sympy")
class Fluids(unittest.TestCase):
    def setUp(self):
        import sympy
        import verify_metrics
        self.sp, self.vm = sympy, verify_metrics
        metric = json.loads((DATA / "metrics" / "sultana_dyer.json").read_text(encoding="utf-8"))
        self.charts = {c["id"]: c for c in metric["coordinates"]}

    def read(self, system, field, variant=None):
        entry = self.charts[system]
        reader = self.vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
        block = entry[field] if variant is None else entry[field]["variants"][variant]["nonzero"]
        return reader, {tuple(e["indices"]): reader(e["value"]) for e in block}

    def test_the_matter_is_a_dust_and_a_null_dust(self):
        """The published G_ab of the Kerr-Schild chart is mu u_a u_b + tau k_a k_b for Saida, Harada and
        Maeda's two fluids at M = r_s/2, with u a unit timelike vector, k null and k.u = -1, and no
        stress on the sphere."""
        sp = self.sp
        reader, G = self.read("kerr_schild", "einstein_tensor", "ll")
        _, g = self.read("kerr_schild", "metric_components")
        eta, r = reader.symbol["\\eta"], reader.symbol["r"]
        rs, e0 = reader.parameters["r_s"], reader.parameters["eta_0"]
        D = r ** 2 + rs * (r - eta)
        root = sp.sqrt(D)
        u = [e0 ** 2 * (2 * r ** 2 + rs * (2 * r - eta)) / (2 * r * eta ** 2 * root),
             -e0 ** 2 * rs * (2 * r - eta) / (2 * r * eta ** 2 * root)]
        k = [e0 ** 2 * root / (r * eta ** 2), -e0 ** 2 * root / (r * eta ** 2)]
        mu = 12 * e0 ** 4 * D / (r ** 2 * eta ** 6)
        tau = e0 ** 4 * rs * (8 * r ** 2 + 3 * rs * (2 * r - eta)) / (r ** 2 * eta ** 5 * D)
        plane = ["\\eta", "r"]
        h = [[g.get((a, b), 0) for b in plane] for a in plane]

        def dot(x, y):
            return sum(h[i][j] * x[i] * y[j] for i in range(2) for j in range(2))
        self.assertEqual(sp.simplify(dot(u, u) + 1), 0)
        self.assertEqual(sp.simplify(dot(k, k)), 0)
        self.assertEqual(sp.simplify(dot(u, k) + 1), 0)
        ul = [sum(h[i][j] * u[j] for j in range(2)) for i in range(2)]
        kl = [sum(h[i][j] * k[j] for j in range(2)) for i in range(2)]
        for i, a in enumerate(plane):
            for j, b in enumerate(plane):
                self.assertEqual(sp.simplify(G.get((a, b), 0) - mu * ul[i] * ul[j] - tau * kl[i] * kl[j]), 0, (a, b))
        self.assertNotIn(("\\theta", "\\theta"), G)
        self.assertNotIn(("\\phi", "\\phi"), G)
        # The density of the dust changes sign on eta = r(r + r_s)/r_s, which is eta = 2 r_s on the horizon.
        self.assertEqual(sp.solve(sp.Eq(D, 0), eta), [r * (r + rs) / rs])
        self.assertEqual(sp.simplify((r * (r + rs) / rs).subs(r, rs) - 2 * rs), 0)

    def test_the_horizon_has_finite_curvature_and_the_big_bang_and_the_centre_do_not(self):
        sp = self.sp
        reader, _ = self.read("kerr_schild", "metric_components")
        entry = self.charts["kerr_schild"]
        eta, r = reader.symbol["\\eta"], reader.symbol["r"]
        rs, e0 = reader.parameters["r_s"], reader.parameters["eta_0"]
        R = reader(entry["ricci_scalar"].partition("=")[2])
        K = reader(entry["kretschmann"].partition("=")[2])
        self.assertEqual(sp.simplify(R.subs(r, rs) - 12 * e0 ** 4 * (2 * rs - eta) / (eta ** 6 * rs)), 0)
        at = {rs: 1, e0: 3}
        self.assertLess(abs(float(K.subs(at).subs({eta: 3, r: 1}))), 1e3)
        self.assertGreater(float(K.subs(at).subs({eta: sp.Rational(1, 1000), r: 2})), 1e20)
        self.assertGreater(float(K.subs(at).subs({eta: 3, r: sp.Rational(1, 1000)})), 1e12)

    def test_the_metric_is_schwarzschilds_times_the_scale_factor_squared(self):
        """Each chart's published metric divided by eta^4/eta_0^4 is Schwarzschild's in that chart:
        Kerr-Schild, static and ingoing Eddington-Finkelstein, with eta = v - r in the last."""
        sp = self.sp
        for system, time, want in (
                ("kerr_schild", "\\eta", lambda f: (-f, 1 - f, 1 + (1 - f))),
                ("schwarzschild_time", "t", lambda f: (-f, 0, 1 / f)),
                ("eddington_finkelstein_ingoing", "v", lambda f: (-f, 1, 0))):
            reader, g = self.read(system, "metric_components")
            r, rs, e0 = reader.symbol["r"], reader.parameters["r_s"], reader.parameters["eta_0"]
            eta = {"kerr_schild": reader.symbol.get("\\eta"), "schwarzschild_time": reader.parameters.get("eta"),
                   "eddington_finkelstein_ingoing": reader.symbol.get("v", 0) - r}[system]
            scale = eta ** 4 / e0 ** 4
            f = 1 - rs / r
            got = (g.get((time, time), 0), g.get((time, "r"), 0), g.get(("r", "r"), 0))
            for have, wanted in zip(got, want(f)):
                self.assertEqual(sp.simplify(have / scale - wanted), 0, system)
            self.assertEqual(sp.simplify(g[("\\theta", "\\theta")] / scale - r ** 2), 0, system)


if __name__ == "__main__":
    unittest.main()
