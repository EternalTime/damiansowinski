"""The Kaluza-Klein black holes, the charged holes that vacuum gravity in five dimensions shows
in four: python3 -m unittest discover -s _tools

What the texts state about them, held on the published files. The drawings are held to the
horizon, the ergosurface and the closed forms of the embedded surfaces from the numbers written.
The physics is taken from the published metrics, in units with G = c = 1 in four dimensions: the
mass, the charge, the temperature, the entropy and the potential of the electric hole, which
have to satisfy Smarr's relation M = 2TS + Phi Q, the speed of the string along the circle, the
reduction of each chart of five dimensions to four, and the period of the circle that removes
the Dirac string. Those tests need sympy and are skipped where it is absent.
"""
import importlib.util
import json
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
NAME = "kaluza_klein_black_hole"
HAS_SYMPY = importlib.util.find_spec("sympy") is not None
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))


def load(kind):
    return json.loads((DATA / kind / f"{NAME}.json").read_text(encoding="utf-8"))


class Drawings(unittest.TestCase):
    def marker(self, system, kind):
        """The radii at which a spacetime diagram marks lines of one kind, in units of r_s."""
        view = load("diagrams")["systems"][system][0]
        x0, x1 = view["box"][:2]
        return sorted(round(x0 + line[0][0] * (x1 - x0), 3)
                      for m in view["markers"] if m["kind"] == kind for line in m.get("lines", []))

    def test_the_horizons_and_the_ergosurface_stand_where_the_metric_puts_them(self):
        """One charge, q = 2 r_s or p = 2 r_s: the cones close at r_s, and in five dimensions
        g_tt vanishes at r = q, the ergosurface of the moving string. Equal charges, p = 3 r_s:
        the horizons are (p + r_s)/2 and (p - r_s)/2."""
        for system in ("electric", "eddington_finkelstein_ingoing", "magnetic", "einstein", "einstein_eddington_finkelstein"):
            self.assertEqual(self.marker(system, "grr"), [1.0], system)
        for system in ("electric", "eddington_finkelstein_ingoing"):
            self.assertEqual(self.marker(system, "gtt"), [2.0], system)
        self.assertEqual(self.marker("magnetic", "gtt"), [])
        self.assertEqual(self.marker("dyonic", "grr"), [1.0, 2.0])

    def test_every_plane_ends_on_its_singularity(self):
        for system, views in load("diagrams")["systems"].items():
            edges = [m["edges"] for m in views[0]["markers"] if m["kind"] == "singular"]
            self.assertEqual(edges, [["left"]], system)

    def test_the_embedded_surfaces_are_the_ones_known_in_closed_form(self):
        """At r_s = 1 with a = q - r_s = p - r_s = 1. The Einstein metric's equator has circles of
        radius (r^3 (r + a))^(1/4); the electric hole's circle of the fifth dimension, drawn at
        L = 2 pi r_s, has radius sqrt(1 + a/r), widest on the horizon; the magnetic hole's has
        radius 2 sqrt(2) sqrt(r/(r + a)), narrowest there. Each profile climbs at
        sqrt(g_rr - (d rho/dr)^2), and both sheets start on the horizon."""
        k = 2 * math.sqrt(2)
        cases = {
            "einstein": (lambda r: (r ** 3 * (r + 1)) ** 0.25,
                         lambda r: math.sqrt(r * (r + 1)) / (r - 1) - ((4 * r + 3) / (4 * (r * (r + 1) ** 3) ** 0.25)) ** 2),
            "electric": (lambda r: math.sqrt(1 + 1 / r), lambda r: r / (r - 1) - 1 / (4 * r ** 3 * (r + 1))),
            "magnetic": (lambda r: k * math.sqrt(r / (r + 1)), lambda r: (r + 1) / (r - 1) - 2 / (r * (r + 1) ** 3)),
        }
        throats = {"einstein": 2 ** 0.25, "electric": math.sqrt(2), "magnetic": 2.0}
        views = {view["id"]: view["surfaces"][0] for view in load("embedding")["views"]}
        self.assertEqual(sorted(views), sorted(cases))
        for vid, (rho_of, rise2) in cases.items():
            for sign, pid in ((1, "exterior"), (-1, "other_exterior")):
                points = next(p for p in views[vid]["pieces"] if p["id"] == pid)["points"]
                self.assertEqual(points[0][0], 1.0)
                self.assertEqual(points[0][2], 0.0)
                self.assertAlmostEqual(points[0][1], throats[vid], places=6)
                self.assertEqual(points[-1][0], 6.0)
                for r, rho, z in points:
                    self.assertLess(abs(rho - rho_of(r)), 2e-6, f"{vid} {pid}, rho at {r}")
                    self.assertGreaterEqual(sign * z, 0.0)
                for (r0, _, z0), (r1, _, z1) in zip(points, points[1:]):
                    if r0 < 1.2:
                        continue
                    mid = 0.5 * (r0 + r1)
                    # A chord's slope against the slope at its middle, to the curvature of the profile over it.
                    want = math.sqrt(rise2(mid))
                    self.assertLess(abs(sign * (z1 - z0) / (r1 - r0) - want), 5e-3 * want,
                                    f"{vid} {pid}, dz/dr at {mid}")
        # The electric hole's circle grows toward the hole and the magnetic hole's shrinks.
        electric = next(p for p in views["electric"]["pieces"] if p["id"] == "exterior")["points"]
        magnetic = next(p for p in views["magnetic"]["pieces"] if p["id"] == "exterior")["points"]
        self.assertTrue(all(b[1] < a[1] for a, b in zip(electric, electric[1:])))
        self.assertTrue(all(b[1] > a[1] for a, b in zip(magnetic, magnetic[1:])))


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class Physics(unittest.TestCase):
    def setUp(self):
        import sympy
        import verify_metrics
        self.sp, self.vm = sympy, verify_metrics
        self.metric = load("metrics")

    def chart(self, system):
        """The published metric of a chart as a sympy matrix, with its coordinates and parameters."""
        entry = next(c for c in self.metric["coordinates"] if c["id"] == system)
        reader = self.vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
        x = [reader.symbol[name] for name in entry["coords"]]
        g = self.sp.zeros(len(x), len(x))
        for e in entry["metric_components"]:
            i, j = (entry["coords"].index(name) for name in e["indices"])
            g[i, j] = reader(e["value"])
        return g, x, reader.parameters

    def zero(self, expression, values):
        self.assertLess(abs(complex(self.sp.sympify(expression).subs(values).evalf(30))), 1e-20)

    def points(self, *names):
        """Three sets of values with r_s < the charge length < r."""
        sp = self.sp
        for rs, charge, r in ((1, sp.Rational(7, 3), 5), (sp.Rational(3, 2), 2, sp.Rational(11, 2)),
                              (sp.Rational(1, 4), 6, sp.Rational(9, 5))):
            yield rs, charge, r

    def test_the_electric_hole_reduces_to_the_einstein_metric_and_obeys_smarrs_relation(self):
        """ds_5^2 = g_yy (dy + 2A)^2 + g_yy^(-1/2) g_4. From the published charts: the mass
        M = (q + r_s)/4 from the Einstein metric far away, the charge Q = sqrt(q(q - r_s))/2 from
        A_t = -Q/R, the entropy a quarter of the horizon's area, the temperature kappa/2 pi, and
        the potential -A_t on the horizon, which satisfy M = 2TS + Phi Q; and Q < 2M, the bound
        the hole approaches as r_s -> 0."""
        sp = self.sp
        g5, x5, p5 = self.chart("electric")
        g4, x4, p4 = self.chart("einstein")
        r, th = x5[1], x5[2]
        rs, q = p5["r_s"], p5["q"]
        same = {x4[1]: r, x4[2]: th, p4["r_s"]: rs, p4["q"]: q}
        scalar = g5[4, 4]
        A = [g5[a, 4] / (2 * scalar) for a in range(4)]
        for RS, Qv, R in self.points():
            at = {rs: RS, q: Qv, r: R, th: sp.Rational(4, 5)}
            for a in range(4):
                for b in range(a, 4):
                    reduced = sp.sqrt(scalar) * (g5[a, b] - 4 * scalar * A[a] * A[b])
                    self.zero(reduced - g4[a, b].subs(same), at)
        f, grr, area = -g4[0, 0].subs(same), g4[1, 1].subs(same), g4[2, 2].subs(same)
        radius = sp.sqrt(area)
        M = sp.limit(radius * (1 - f) / 2, r, sp.oo)
        Q = sp.limit(-radius * A[0], r, sp.oo)
        kappa = sp.limit(sp.diff(f, r) / (2 * sp.sqrt(f * grr)), r, rs, "+")
        S = sp.pi * area.subs(r, rs)
        Phi = -A[0].subs(r, rs)
        for RS, Qv, _ in self.points():
            at = {rs: RS, q: Qv}
            self.zero(M - (q + rs) / 4, at)
            self.zero(Q - sp.sqrt(q * (q - rs)) / 2, at)
            self.zero(kappa - 1 / (2 * sp.sqrt(q * rs)), at)
            self.zero(M - 2 * kappa / (2 * sp.pi) * S - Phi * Q, at)
            self.assertLess(float((Q / M).subs(at)), 2.0)
        self.zero(sp.limit(Q / M, rs, 0) - 2, {q: 3})

    def test_the_string_moves_at_tanh_alpha_and_its_horizon_is_a_killing_horizon(self):
        """The Killing vector d_t + v d_y with v = -g_ty/g_yy on the horizon is null there, and
        v = sqrt(1 - r_s/q), the speed of a boost with cosh^2(alpha) = q/r_s."""
        sp = self.sp
        g, x, p = self.chart("electric")
        r, rs, q = x[1], p["r_s"], p["q"]
        speed = sp.limit(-g[0, 4] / g[4, 4], r, rs)
        norm = g[0, 0] + 2 * speed * g[0, 4] + speed ** 2 * g[4, 4]
        for RS, Qv, _ in self.points():
            at = {rs: RS, q: Qv}
            self.zero(speed - sp.sqrt(1 - rs / q), at)
            self.zero(sp.limit(norm, r, rs), at)

    def test_the_magnetic_circle_has_the_period_that_removes_the_dirac_string(self):
        """The fibre is dy + 2P(1 - cos theta) d phi with 2P = g_phi_y/(g_yy (1 - cos theta)), and
        the string on theta = pi is removed by y -> y + 4P phi, a change of chart when y has the
        period 8 pi P = 4 pi sqrt(p(p - r_s)), the period the chart publishes. The hole has the
        Einstein metric of the electric hole with p in place of q."""
        sp = self.sp
        g, x, par = self.chart("magnetic")
        g4, x4, p4 = self.chart("einstein")
        r, th, rs, p = x[1], x[2], par["r_s"], par["p"]
        twoP = g[3, 4] / (g[4, 4] * (1 - sp.cos(th)))
        entry = next(c for c in self.metric["coordinates"] if c["id"] == "magnetic")
        self.assertIn("y \\in \\left[0, 4\\pi\\sqrt{p\\left(p - r_s\\right)}\\right)", entry["domains"])
        same = {x4[1]: r, x4[2]: th, p4["r_s"]: rs, p4["q"]: p}
        scalar = g[4, 4]
        fibre = [g[a, 4] / scalar for a in range(4)]
        for RS, Pv, R in self.points():
            at = {rs: RS, p: Pv, r: R, th: sp.Rational(4, 5)}
            self.zero(twoP ** 2 - p * (p - rs), at)
            for a in range(4):
                for b in range(a, 4):
                    reduced = sp.sqrt(scalar) * (g[a, b] - scalar * fibre[a] * fibre[b])
                    self.zero(reduced - g4[a, b].subs(same), at)

    def test_equal_charges_keep_the_circle_and_leave_reissner_and_nordstrom(self):
        """g_yy = 1, the electric and the magnetic potential have one charge, 2Q = 2P =
        sqrt((p^2 - r_s^2)/2), and the metric orthogonal to the circles is -f dt^2 + d rho^2/f +
        rho^2 d Omega^2 with f = (1 - rho_+/rho)(1 - rho_-/rho), rho_+- = (p +- r_s)/2."""
        sp = self.sp
        g, x, par = self.chart("dyonic")
        rho, th, rs, p = x[1], x[2], par["r_s"], par["p"]
        self.assertEqual(g[4, 4], 1)
        f = (1 - (p + rs) / (2 * rho)) * (1 - (p - rs) / (2 * rho))
        want = sp.diag(-f, 1 / f, rho ** 2, rho ** 2 * sp.sin(th) ** 2)
        for RS, Pv, R in self.points():
            at = {rs: RS, p: Pv, rho: R + Pv, th: sp.Rational(4, 5)}
            self.zero((rho * g[0, 4]) ** 2 - (p ** 2 - rs ** 2) / 2, at)
            self.zero((g[3, 4] / (1 - sp.cos(th))) ** 2 - (p ** 2 - rs ** 2) / 2, at)
            for a in range(4):
                for b in range(a, 4):
                    self.zero(g[a, b] - g[a, 4] * g[b, 4] - want[a, b], at)


if __name__ == "__main__":
    unittest.main()
