"""The moving mirror of Fulling and Davies, flat spacetime of two dimensions to the right of a
perfectly reflecting mirror:
.venv.noindex/bin/python -m unittest discover -s _tools

What its texts and drawings state, held on the published files. The first class needs nothing
but the files: the mirror's world line on each spacetime diagram, each marked ray meeting its
reflection on the mirror at the retarded time the captions state, the horizon, the mirror and
the last ray on the conformal diagrams, the flux drawn as a height, and the relations each side
answers. The second reads the published tensors through the checker's Reader and holds every
chart flat and each mirror's flux to its closed form; it needs sympy and is skipped where it is
absent.
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

CHARTS = ["inertial", "null", "mirror_rest", "thermal", "collapse", "rindler"]
ARRIVE = {"thermal": (-2, -1, -0.5, -0.25, -0.125), "collapse": (-2, -1, -0.5, -0.25, -0.125),
          "uniform": (0.5, 1, 2, 4)}


def load(folder, metric_id="moving_mirror"):
    return json.loads((DATA / folder / f"{metric_id}.json").read_text(encoding="utf-8"))


def lambert(x):
    """Lambert's function on its principal branch for x >= 0, by Newton's method on w e^w = x."""
    w = math.log1p(x) if x < 3 else math.log(x) - math.log(math.log(x))
    for _ in range(60):
        e = math.exp(w)
        w -= (w * e - x) / (e * (w + 1))
    return w


# Each mirror at kappa = 1: the retarded time u at which the ray that arrives along v leaves it.
LEAVES = {"thermal": lambda v: -math.log(-v), "collapse": lambda v: v - math.log(-v), "uniform": lambda v: -1 / v}
# Its position at the time t.
WORLD = {"thermal": lambda t: -t - lambert(math.exp(-2 * t)), "collapse": lambda t: -t - lambert(2 * math.exp(-2 * t)) / 2,
         "uniform": lambda t: math.sqrt(1 + t * t)}


class Drawings(unittest.TestCase):
    def views(self, system):
        return {view["id"]: view for view in load("diagrams")["systems"][system]}

    def chart_point(self, view, point):
        """A point of a view's unit square as the view's chart coordinates (x^0, r)."""
        X0, X1, Y0, Y1 = view["box"]
        X, Y = X0 + point[0] * (X1 - X0), Y0 + point[1] * (Y1 - Y0)
        (a, b), (c, d) = view["to_display"]
        det = a * d - b * c
        return (d * X - b * Y) / det, (a * Y - c * X) / det

    def null_pair(self, system, view, point):
        """The inertial u and v of a point of a view, at kappa = 1."""
        x0, r = self.chart_point(view, point)
        case = view["id"] if view["id"] in LEAVES else system
        if system == "inertial":
            return x0 - r, x0 + r
        if system == "null":
            return x0, r
        if system == "mirror_rest":
            return LEAVES[case](x0) if x0 < 0 else math.inf, r
        if system in ("thermal", "collapse"):
            U = x0 - r
            return LEAVES[system](U) if U < 0 else math.inf, x0 + r
        return -math.exp(r - x0), math.exp(r + x0)

    def test_every_chart_has_its_spacetime_diagrams_and_every_window_is_square(self):
        systems = load("diagrams")["systems"]
        self.assertEqual(sorted(systems), sorted(CHARTS))
        self.assertEqual({s: sorted(v["id"] for v in views) for s, views in systems.items()},
                         {"inertial": ["collapse", "thermal", "uniform"], "null": ["collapse", "thermal"],
                          "mirror_rest": ["collapse", "thermal"], "thermal": ["tx"], "collapse": ["tx"], "rindler": ["tx"]})
        for system, views in systems.items():
            for view in views:
                X0, X1, Y0, Y1 = view["box"]
                self.assertAlmostEqual((X1 - X0) / (Y1 - Y0), 1.0, places=9, msg=f"{system} {view['id']}")

    def test_the_world_line_drawn_is_the_declared_mirror(self):
        """x = z(t) on the inertial planes, v = p(u) on the null ones, and the left edge where the
        mirror is at rest."""
        for case, view in self.views("inertial").items():
            (mark,) = [m for m in view["markers"] if m["kind"] == "world"]
            (line,) = mark["lines"]
            self.assertGreater(len(line), 5, case)
            for point in line:
                t, x = self.chart_point(view, point)
                self.assertAlmostEqual(x, WORLD[case](t), delta=6e-3, msg=f"{case} at t = {t}")
        for case, view in self.views("null").items():
            (mark,) = [m for m in view["markers"] if m["kind"] == "world"]
            for point in mark["lines"][0]:
                u, v = self.chart_point(view, point)
                self.assertAlmostEqual((v - u) / 2, WORLD[case]((u + v) / 2), delta=6e-3, msg=case)
        for system in ("mirror_rest", "thermal", "collapse", "rindler"):
            for view in self.views(system).values():
                (mark,) = [m for m in view["markers"] if m["kind"] == "world"]
                self.assertEqual(mark["lines"], [[[0.0, 0.0], [0.0, 1.0]]], system)

    def test_each_marked_ray_meets_its_reflection_on_the_mirror(self):
        """A ray arrives along v, meets the mirror, and leaves along u = f(v): -ln(-v) for Carlitz
        and Willey's mirror, v - ln(-v) for the mirror of Good, Anderson and Evans, and -1/v for
        uniform acceleration."""
        systems = load("diagrams")["systems"]
        for system, views in systems.items():
            for view in views:
                case = view["id"] if view["id"] in LEAVES else {"rindler": "uniform"}.get(system, system)
                marks = [m for m in view["markers"] if m["kind"] == "shell"]
                self.assertEqual(len(marks), len(ARRIVE[case]), f"{system} {view['id']}")
                for mark, v in zip(marks, ARRIVE[case]):
                    arriving, leaving = mark["lines"]
                    # Both halves have an end on the mirror, the same point.
                    ends = [min(line, key=lambda p: abs(self.gap(system, case, view, p))) for line in (arriving, leaving)]
                    for a, b in zip(*ends):
                        self.assertAlmostEqual(a, b, delta=1e-3, msg=f"{system} {view['id']} v = {v}")
                    self.assertAlmostEqual(self.gap(system, case, view, ends[0]), 0.0, delta=4e-3)
                    # The arriving half keeps v and the leaving half keeps u = f(v).
                    far = [max(line, key=lambda p: abs(self.gap(system, case, view, p))) for line in (arriving, leaving)]
                    _, got_v = self.null_pair(system, view, far[0])
                    got_u, _ = self.null_pair(system, view, far[1])
                    self.assertAlmostEqual(got_v, v, delta=6e-3, msg=f"{system} {view['id']}")
                    self.assertAlmostEqual(got_u, LEAVES[case](v), delta=2e-2, msg=f"{system} {view['id']} v = {v}")

    def gap(self, system, case, view, point):
        """How far a point of a view lies from the mirror, zero on it."""
        x0, r = self.chart_point(view, point)
        if system == "inertial":
            return r - WORLD[case](x0)
        if system == "null":
            return (r - x0) / 2 - WORLD[case]((x0 + r) / 2)
        if system == "mirror_rest":
            return r - x0
        return r

    def test_the_reflections_leave_as_the_captions_say(self):
        """Carlitz and Willey's five rays leave ln 2 apart; the other mirror's leave 1.69, 1.19, 0.94
        and 0.82 apart, at u = -2.69, -1, 0.19, 1.14 and 1.95; the hyperbola's at -2, -1, -1/2, -1/4."""
        thermal = [LEAVES["thermal"](v) for v in ARRIVE["thermal"]]
        for a, b in zip(thermal, thermal[1:]):
            self.assertAlmostEqual(b - a, math.log(2), places=12)
        collapse = [LEAVES["collapse"](v) for v in ARRIVE["collapse"]]
        for got, want in zip(collapse, (-2.69, -1, 0.19, 1.14, 1.95)):
            self.assertAlmostEqual(got, want, delta=5e-3)
        for got, want in zip((b - a for a, b in zip(collapse, collapse[1:])), (1.69, 1.19, 0.94, 0.82)):
            self.assertAlmostEqual(got, want, delta=5e-3)
        self.assertEqual([LEAVES["uniform"](v) for v in ARRIVE["uniform"]], [-2, -1, -0.5, -0.25])
        caption = " ".join(self.views("inertial")["collapse"]["caption"])
        for number in ("$1.69$", "$1.19$", "$0.94$", "$0.82$"):
            self.assertIn(number, caption)
        # Carlitz and Willey's mirror turns at t = -1/2, x = -1/2.
        self.assertAlmostEqual(WORLD["thermal"](-0.5), -0.5, places=12)
        self.assertTrue(all(WORLD["thermal"](-0.5) >= WORLD["thermal"](-0.5 + d) for d in (-0.2, -0.01, 0.01, 0.2)))
        self.assertIn("$x = -1/2\\kappa$", self.views("inertial")["thermal"]["caption"][0])

    def test_the_horizon_is_the_ray_v_equal_to_zero_and_only_the_mirrors_that_recede_have_one(self):
        systems = load("diagrams")["systems"]
        for system, views in systems.items():
            for view in views:
                events = [m for m in view["markers"] if m["kind"] == "event"]
                receding = system != "rindler" and view["id"] != "uniform"
                self.assertEqual(len(events), 1 if receding else 0, f"{system} {view['id']}")
                for mark in events:
                    for point in mark["lines"][0]:
                        _, v = self.null_pair(system, view, point)
                        self.assertAlmostEqual(v, 0.0, delta=4e-3, msg=f"{system} {view['id']}")

    def test_what_lies_behind_a_mirror_is_hatched(self):
        """Every view hatches the far side of its mirror, and the charts that bring a receding mirror
        to rest hatch what lies beyond U = 0 as well."""
        for system, views in load("diagrams")["systems"].items():
            for view in views:
                self.assertTrue(view["hatch"], f"{system} {view['id']}")
        for system in ("thermal", "collapse"):
            (view,) = load("diagrams")["systems"][system]
            X0, X1, Y0, Y1 = view["box"]
            corner = [((1.0 - X0) / (X1 - X0)), 1.0]       # X = cT = 1 on the top edge
            self.assertTrue(any(min(math.dist(corner, p) for p in ring) < 0.02 for ring in view["hatch"]), system)

    def conformal(self):
        return {view["id"]: view for view in load("conformal")["views"]}

    def test_the_conformal_diagrams_draw_each_mirror_where_its_ray_tracing_function_puts_it(self):
        """With p = arctan u and q = arctan v, X = q - p and T = p + q: v = -e^(-u) and
        v = -W(e^(-u)) from i- to (-pi/2, pi/2), and the hyperbola uv = -1 on the line X = pi/2."""
        views = self.conformal()
        self.assertEqual(sorted(views), ["collapse", "inertial_collapse", "inertial_thermal", "inertial_uniform",
                                         "mirror_rest", "null", "rindler", "thermal"])
        trace = {"thermal": lambda u: -math.exp(-u), "collapse": lambda u: -lambert(math.exp(-u))}
        cases = {"inertial_thermal": "thermal", "thermal": "thermal", "inertial_collapse": "collapse", "null": "collapse",
                 "mirror_rest": "collapse", "collapse": "collapse"}
        for vid, case in cases.items():
            lines = [layer["points"] for layer in views[vid]["layers"] if layer["class"] == "world"]
            self.assertEqual(len(lines), 1, vid)
            line = lines[0]
            for X, T in line[1:-1]:
                p, q = (T - X) / 2, (T + X) / 2
                if abs(abs(p) - math.pi / 2) < 0.05 or abs(q) < 1e-3 or abs(q + math.pi / 2) < 0.05:
                    continue
                self.assertAlmostEqual(math.atan(trace[case](math.tan(p))), q, delta=3e-3, msg=vid)
            for got, want in zip(line[0], (0, -math.pi)):
                self.assertAlmostEqual(got, want, delta=2e-3, msg=vid)
            for got, want in zip(line[-1], (-math.pi / 2, math.pi / 2)):
                self.assertAlmostEqual(got, want, delta=2e-3, msg=vid)
            (horizon,) = [layer["points"] for layer in views[vid]["layers"] if layer["class"] == "horizon"]
            for X, T in horizon:
                self.assertAlmostEqual(T + X, 0.0, places=3, msg=vid)      # q = 0, the ray v = 0
        for vid in ("inertial_uniform", "rindler"):
            (line,) = [layer["points"] for layer in views[vid]["layers"] if layer["class"] == "world"]
            for X, T in line:
                self.assertAlmostEqual(X, math.pi / 2, places=3, msg=vid)
            self.assertFalse([layer for layer in views[vid]["layers"] if layer["class"] == "horizon"], vid)

    def test_every_line_of_a_conformal_diagram_lies_to_the_right_of_its_mirror(self):
        views = self.conformal()
        trace = {"thermal": lambda u: -math.exp(-u), "collapse": lambda u: -lambert(math.exp(-u))}
        cases = {"inertial_thermal": "thermal", "thermal": "thermal", "inertial_collapse": "collapse", "null": "collapse",
                 "mirror_rest": "collapse", "collapse": "collapse"}
        for vid, view in views.items():
            for layer in view["layers"]:
                if layer["kind"] != "line" or layer["class"] not in ("r", "t", "null"):
                    continue
                for X, T in layer["points"]:
                    p, q = (T - X) / 2, (T + X) / 2
                    if vid in cases:
                        if abs(p) > math.pi / 2 - 0.02:
                            continue
                        self.assertGreaterEqual(q, math.atan(trace[cases[vid]](math.tan(p))) - 4e-3, f"{vid} {layer['class']}")
                    else:
                        self.assertGreaterEqual(X, math.pi / 2 - 2e-3, f"{vid} {layer['class']}")

    def test_the_height_is_the_flux_of_the_mirror_that_imitates_a_collapse(self):
        """Twice (4W + 1)/(W + 1)^4 with W = W(e^(-u)) on the ray u = ct - x, over -3 <= ct <= 3
        across and 0.5 <= x <= 6.5 away: (4.2) of Good, Anderson and Evans, in units of the thermal flux."""
        (view,) = load("embedding")["views"]
        self.assertEqual(view["id"], "flux")
        self.assertIn("height", view)
        (surface,) = view["surfaces"]
        (piece,) = surface["pieces"]
        grid = piece["grid"]
        self.assertEqual(grid["frame"], "cartesian")
        self.assertEqual((float(grid["u"][0]), float(grid["u"][-1]), float(grid["v"][0]), float(grid["v"][-1])),
                         (-3.0, 3.0, -3.0, 3.0))
        worst, top, low = 0.0, 0.0, 9.0
        for X, row in zip(grid["u"], grid["z"]):
            for Y, z in zip(grid["v"], row):
                w = lambert(math.exp(-(float(X) - (float(Y) + 3.5))))
                want = 2 * (4 * w + 1) / (w + 1) ** 4
                worst = max(worst, abs(float(z) - want))
                top, low = max(top, float(z)), min(low, float(z))
        self.assertLess(worst, 1e-6)
        self.assertLess(top, 2.0)                # below the thermal flux everywhere
        self.assertGreater(top, 1.9)             # and within 3% of it on the last ray drawn
        self.assertGreaterEqual(low, 0.0)        # no ray carries negative energy
        self.assertLess(low, 0.02)
        # The plane lies to the right of the mirror at every moment of it.
        for k in range(61):
            self.assertLess(WORLD["collapse"](-3 + 0.1 * k), 0.5)

    def test_the_marked_ray_is_where_the_flux_is_sixteen_twenty_sevenths_of_its_last_value(self):
        (view,) = load("embedding")["views"]
        (curve,) = view["surfaces"][0]["curves"]
        self.assertEqual(curve["class"], "path")
        u = math.log(2) - 0.5
        w = lambert(math.exp(-u))
        self.assertAlmostEqual(w, 0.5, places=12)                       # W(e^(1/2)/2) = 1/2
        self.assertAlmostEqual((4 * w + 1) / (w + 1) ** 4, 16 / 27, places=12)
        for X, Y, Z in curve["points"]:
            self.assertAlmostEqual(X - (Y + 3.5), u, delta=4e-3)
            self.assertAlmostEqual(Z, 32 / 27, delta=1e-3)

    def test_each_relation_is_answered(self):
        mirror = load("metrics")
        answers = {"locally_same": "locally_same", "family": "family"}
        self.assertEqual([r["id"] for r in mirror["related"]], ["minkowski", "schwarzschild", "vaidya", "witten_black_hole"])
        for relation in mirror["related"]:
            other = load("metrics", relation["id"])
            (back,) = [r for r in other["related"] if r["id"] == "moving_mirror"]
            self.assertEqual(back["kind"], answers[relation["kind"]], relation["id"])


@unittest.skipUnless(HAS_SYMPY, "needs sympy")
class Mathematics(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import sympy as sp
        import verify_metrics as vm
        cls.sp, cls.vm = sp, vm
        cls.metric = load("metrics")
        cls.charts = {chart["id"]: chart for chart in cls.metric["coordinates"]}

    def reader(self, system):
        chart = self.charts[system]
        return self.vm.Reader(chart["coords"], [p["symbol"] for p in chart["parameters"]], ())

    def test_every_chart_is_flat(self):
        self.assertEqual(sorted(self.charts), sorted(CHARTS))
        for system, chart in self.charts.items():
            for tensor in ("riemann", "ricci_tensor", "einstein_tensor", "weyl_tensor"):
                for name, variant in chart[tensor]["variants"].items():
                    self.assertEqual(variant["nonzero"], [], f"{system} {tensor} {name}")
            self.assertEqual(chart["ricci_scalar"], "R = 0", system)
            self.assertEqual(chart["kretschmann"], "K = 0", system)

    def flux(self, factor, U):
        """24 pi T_uu from the factor C of -C dU dv, by Davies, Fulling and Unruh's formula."""
        sp = self.sp
        return sp.simplify(-2 * sp.sqrt(factor) * sp.diff(1 / sp.sqrt(factor), U, 2) / factor ** 2)

    def factor(self, system):
        """The published g_XX of a chart of T and X as a function of U = cT - X, at c = 1."""
        sp = self.sp
        R = self.reader(system)
        (value,) = [c["value"] for c in self.charts[system]["metric_components"] if c["indices"] == ["X", "X"]]
        U = sp.Symbol("U", negative=True)
        T, X = R.symbol["T"], R.symbol["X"]
        factor = sp.simplify(R(value).subs(R.c, 1).subs(X, T - U))
        self.assertNotIn(T, factor.free_symbols)
        return factor, U, R.parameters["kappa"]

    def test_carlitz_and_willeys_mirror_radiates_the_thermal_flux_at_every_time(self):
        factor, U, kappa = self.factor("thermal")
        self.assertEqual(self.sp.simplify(self.flux(factor, U) - kappa ** 2 / 2), 0)

    def test_the_other_mirrors_flux_is_good_anderson_and_evanss(self):
        """(4W + 1)/(W + 1)^4 of the thermal flux with W = -kappa U, which is W(e^(-kappa u)) on the
        ray u = U - ln(-kappa U)/kappa: zero far in the past and the thermal flux at U = 0."""
        sp = self.sp
        factor, U, kappa = self.factor("collapse")
        w = -kappa * U
        want = kappa ** 2 * (4 * w + 1) / (2 * (w + 1) ** 4)
        self.assertEqual(sp.simplify(self.flux(factor, U) - want), 0)
        self.assertEqual(sp.limit(want, U, 0, "-"), kappa ** 2 / 2)
        self.assertEqual(sp.limit(want, U, -sp.oo), 0)

    def test_the_chart_that_brings_a_mirror_to_rest_has_one_christoffel_symbol(self):
        chart = self.charts["mirror_rest"]
        self.assertEqual(chart["line_element"], "ds^2 = -f'\\,dU\\,dv")
        (only,) = chart["christoffel"]["variants"]["ull"]["nonzero"]
        self.assertEqual(only["indices"], ["U", "U", "U"])
        R = self.reader("mirror_rest")
        f = R.parameters["f"]
        U = R.symbol["U"]
        self.assertEqual(self.sp.simplify(R(only["value"]) - self.sp.diff(f, U, 2) / self.sp.diff(f, U)), 0)

    def test_the_hyperbola_radiates_nothing(self):
        """Its ray tracing function is p = -1/(kappa^2 u), so u = f(U) = -1/(kappa^2 U) and the
        factor of the chart that brings it to rest is 1/(kappa U)^2, whose flux vanishes."""
        sp = self.sp
        U, kappa = sp.Symbol("U", positive=True), sp.Symbol("kappa", positive=True)
        self.assertEqual(self.flux(sp.diff(-1 / (kappa ** 2 * U), U), U), 0)


if __name__ == "__main__":
    unittest.main()
