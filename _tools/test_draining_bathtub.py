"""The draining bathtub, Visser's acoustic metric of a fluid that swirls down a drain:
.venv.noindex/bin/python -m unittest discover -s _tools

What its texts and drawings state, held on the published files. The first class needs nothing
but the files: the horizon and the ergosurface where every spacetime diagram marks them, the
rays of each chart against the speed of sound through moving water, the cones of the figure in
three dimensions, the two moments of the embedding diagram on every drawing, and the conformal
diagram's edges. The second holds the relations to the answers of the other files. The third
reads the published tensors through the checker's Reader and holds them to Visser's line
element, to sound moving at c past the fluid, to the map between the charts and to the curvature
the captions quote; it needs sympy and is skipped where it is absent.

Every drawing is of the drain A = -1, B = sqrt 3 at c = 1, so the horizon |A|/c is the unit of
length, the ergosurface lies at 2 and the water moves at (-1/r, sqrt(3)/r).
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

CHARTS = ["laboratory", "kerr_like", "vortex_filament"]
SWIRL = math.sqrt(3)


def load(folder, metric_id="draining_bathtub"):
    return json.loads((DATA / folder / f"{metric_id}.json").read_text(encoding="utf-8"))


def flat_views():
    return {(system, view["id"]): view for system, views in load("diagrams")["systems"].items() for view in views}


def chart_point(view, point):
    """A point of a flat view's unit square as (time, r) of its chart."""
    X0, X1, Y0, Y1 = view["box"]
    return Y0 + point[1] * (Y1 - Y0), X0 + point[0] * (X1 - X0)


class Drawings(unittest.TestCase):
    def marked(self, view, kind):
        X0, X1 = view["box"][:2]
        return sorted(X0 + line[0][0] * (X1 - X0) for m in view["markers"] if m["kind"] == kind for line in m["lines"])

    def test_every_view_marks_the_horizon_at_one_and_the_ergosurface_at_two(self):
        views = flat_views()
        self.assertEqual(sorted(views), [("kerr_like", "exterior"), ("laboratory", "drain"), ("laboratory", "spring"),
                                         ("vortex_filament", "drain")])
        for key, view in views.items():
            self.assertEqual(len(self.marked(view, "grr")), 1, key)
            self.assertAlmostEqual(self.marked(view, "grr")[0], 1.0, delta=2e-3, msg=key)
            self.assertEqual(len(self.marked(view, "gtt")), 1, key)
            self.assertAlmostEqual(self.marked(view, "gtt")[0], 2.0, delta=2e-3, msg=key)

    def test_the_drain_and_the_spring_end_on_a_singular_centre_and_the_kerr_like_chart_on_its_horizon(self):
        for key, view in flat_views().items():
            singular = [m for m in view["markers"] if m["kind"] == "singular"]
            if key[0] == "kerr_like":
                self.assertEqual(singular, [], key)
                self.assertEqual(view["box"][0], 1)
            else:
                self.assertEqual([m["edges"] for m in singular], [["left"]], key)
                self.assertEqual(view["box"][0], 0)

    def test_the_drawn_windows_are_square(self):
        for key, view in flat_views().items():
            X0, X1, Y0, Y1 = view["box"]
            self.assertAlmostEqual((X1 - X0) / (Y1 - Y0), 1.0, places=9, msg=key)

    def conserved(self, key):
        """What each family of rays keeps, at |A| = c = 1: the laboratory's rays have
        dr/dt = -+1/r -+ 1 and the Kerr-like chart's c dT/dr = -+r^2/(r^2 - 1)."""
        if key == ("laboratory", "spring"):
            return (lambda t, r: t + r + math.log(abs(r - 1)), lambda t, r: t - r + math.log(1 + r))
        if key[0] == "kerr_like":
            star = lambda r: r + math.log((r - 1) / (r + 1)) / 2
            return (lambda T, r: T + star(r), lambda T, r: T - star(r))
        return (lambda t, r: t + r - math.log(1 + r), lambda t, r: t - r - math.log(abs(r - 1)))

    def test_every_ray_is_sound_moving_at_c_through_the_moving_water(self):
        for key, view in flat_views().items():
            forms = self.conserved(key)
            checked = 0
            for family in view["rays"].values():
                for ray in family:
                    points = [chart_point(view, p) for p in ray]
                    points = [(t, r) for t, r in points if r > 0.15 and abs(r - 1) > 0.1]
                    if len(points) < 2:
                        continue
                    drift = min(max(f(*p) for p in points) - min(f(*p) for p in points) for f in forms)
                    self.assertLess(drift, 0.02, key)
                    checked += 1
            self.assertGreater(checked, 20, key)

    def test_inside_the_horizon_every_cone_of_the_drain_points_in_and_every_cone_of_the_spring_out(self):
        views = flat_views()
        for key, sign in ((("laboratory", "drain"), -1), (("laboratory", "spring"), 1), (("vortex_filament", "drain"), -1)):
            view = views[key]
            inside = [c for c in view["cones"] if chart_point(view, c["at"])[1] < 0.95]
            self.assertGreater(len(inside), 6, key)
            for cone in inside:
                for edge in (cone["a"], cone["b"]):
                    self.assertGreater(sign * edge[0], 0, key)
                    self.assertGreater(edge[1], 0, key)

    def test_the_laboratory_view_draws_both_moments_and_the_kerr_like_view_both(self):
        views = flat_views()
        lab = {m["view"]: m for m in views[("laboratory", "drain")]["slices"]}
        self.assertEqual(sorted(lab), ["funnel", "plane"])
        for t, r in (chart_point(views[("laboratory", "drain")], p) for line in lab["plane"]["lines"] for p in line):
            self.assertAlmostEqual(t, 0.0, delta=1e-3)
        for t, r in (chart_point(views[("laboratory", "drain")], p) for line in lab["funnel"]["lines"] for p in line):
            if r > 1.02:
                self.assertAlmostEqual(t, math.log(r * r - 1) / 2, delta=5e-3)
        self.assertEqual(lab["funnel"]["label"], "$T = 0$")
        kerr = {m["view"]: m for m in views[("kerr_like", "exterior")]["slices"]}
        self.assertEqual(sorted(kerr), ["funnel", "plane"])
        for T, r in (chart_point(views[("kerr_like", "exterior")], p) for line in kerr["funnel"]["lines"] for p in line):
            self.assertAlmostEqual(T, 0.0, delta=1e-3)
        for T, r in (chart_point(views[("kerr_like", "exterior")], p) for line in kerr["plane"]["lines"] for p in line):
            if r > 1.02:
                self.assertAlmostEqual(T, -math.log(r * r - 1) / 2, delta=5e-3)
        for key in (("laboratory", "spring"), ("vortex_filament", "drain")):
            self.assertNotIn("slices", views[key], key)

    def figure(self):
        (figure,) = load("diagrams")["projections"]["laboratory"]
        return figure

    def test_every_generator_of_every_cone_of_the_figure_is_sound_at_c_past_the_water(self):
        """A generator from an event at (x, y) runs along (dx, dy, dt) with |dx/dt - v| = c for
        the water's velocity v = (A r^ + B theta^)/r there."""
        cones = self.figure()["turn"]["cones"]
        self.assertEqual(len(cones), 12)
        for cone in cones:
            x, y, t0 = cone["apex"]
            r = math.hypot(x, y)
            vx = (-x - SWIRL * y) / (r * r)
            vy = (-y + SWIRL * x) / (r * r)
            for X, Y, T in cone["rim"]:
                dt = T - t0
                self.assertGreater(dt, 0)
                self.assertAlmostEqual(math.hypot((X - x) / dt - vx, (Y - y) / dt - vy), 1.0, delta=2e-4)

    def test_the_figures_cones_stand_on_the_horizon_the_ergosurface_and_outside(self):
        by_radius = {}
        for cone in self.figure()["turn"]["cones"]:
            x, y, t0 = cone["apex"]
            r = math.hypot(x, y)
            out = [((X - x) * x + (Y - y) * y) / r / (T - t0) for X, Y, T in cone["rim"]]
            rest = [math.hypot(X - x, Y - y) / (T - t0) for X, Y, T in cone["rim"]]
            by_radius.setdefault(round(r, 3), []).append((max(out), min(out), min(rest)))
        self.assertEqual(sorted(by_radius), [1.0, 2.0, 3.2])
        for outward, _, rest in by_radius[1.0]:
            # On the horizon the outermost generator runs along the cylinder, and none leads out.
            self.assertAlmostEqual(outward, 0.0, delta=1e-3)
            self.assertGreater(rest, 0.5)
        for outward, inward, rest in by_radius[2.0]:
            # On the ergosurface one generator stands vertical: sound sent upstream stays put.
            self.assertAlmostEqual(rest, 0.0, delta=2e-3)
            self.assertGreater(outward, 0)
        for outward, inward, rest in by_radius[3.2]:
            self.assertGreater(outward, 0)
            self.assertLess(inward, 0)
            self.assertGreater(rest, 0.3)

    def test_the_figures_floor_is_the_laboratorys_moment_with_spirals_of_the_flow_on_it(self):
        figure = self.figure()
        self.assertEqual([(m["view"], m["label"]) for m in figure["slices"]], [("plane", "$t = 0$")])
        spirals = [line["points"] for line in figure["turn"]["lines"] if line["class"] == "floor"][:6]
        for spiral in spirals:
            radii = [math.hypot(x, y) for x, y, _ in spiral]
            angles = [math.atan2(y, x) for x, y, _ in spiral]
            self.assertTrue(all(abs(t) < 1e-9 for _, _, t in spiral))
            # dr/dtheta = A r/B: ln r falls by 1/sqrt 3 for each radian the water turns.
            for i in range(0, len(spiral) - 1, 37):
                if radii[i + 1] < 0.5:
                    continue        # six decimals of a point near the drain no longer fix its angle
                turned = (angles[i + 1] - angles[i] + math.pi) % (2 * math.pi) - math.pi
                self.assertAlmostEqual(math.log(radii[i + 1] / radii[i]), -turned / SWIRL, delta=1e-5)

    def test_the_embedding_is_a_flat_plane_in_the_laboratorys_time_and_a_catenoid_in_the_kerr_like_charts(self):
        views = {view["id"]: view for view in load("embedding")["views"]}
        self.assertEqual(sorted(views), ["funnel", "plane"])
        (plane,) = views["plane"]["surfaces"][0]["pieces"]
        self.assertEqual((plane["system"], plane["coordinate"]), ("laboratory", "r"))
        self.assertEqual((plane["points"][0][0], plane["points"][-1][0]), (0.0, 4.0))
        for r, rho, z in plane["points"]:
            self.assertAlmostEqual(rho, r, places=6)
            self.assertEqual(z, 0.0)
        (funnel,) = views["funnel"]["surfaces"][0]["pieces"]
        self.assertEqual((funnel["system"], funnel["coordinate"]), ("kerr_like", "r"))
        self.assertEqual((funnel["points"][0][0], funnel["points"][-1][0]), (1.0, 5.0))
        for r, rho, z in funnel["points"]:
            self.assertAlmostEqual(rho, r, places=6)
            self.assertAlmostEqual(z, math.acosh(r), places=5)

    def test_the_conformal_diagram_is_the_outside_and_one_hole_with_the_laboratorys_time_ending_on_an_edge(self):
        views = {view["id"]: view for view in load("conformal")["views"]}
        self.assertEqual(list(views), ["drain", "spring", "kerr_like"])
        half = math.pi / 2
        for vid, up in (("drain", 1), ("spring", -1), ("kerr_like", 1)):
            lines = {}
            for layer in views[vid]["layers"]:
                if layer["class"] in ("horizon", "singular", "chartedge"):
                    lines[layer["class"]] = layer["points"]
            for (X, T), want in zip(lines["horizon"], ((0, 0), (half, up * half))):
                self.assertAlmostEqual(X, want[0], places=3)
                self.assertAlmostEqual(T, want[1], places=3)
            for (X, T), want in zip(lines["chartedge"], ((half, -up * half), (-half, up * half))):
                self.assertAlmostEqual(X, want[0], places=3)
                self.assertAlmostEqual(T, want[1], places=3)
            self.assertTrue(all(abs(T - up * half) < 0.06 for _, T in lines["singular"]), vid)

    def test_the_moments_on_the_conformal_diagram_run_from_the_drain_and_from_the_horizon(self):
        views = {view["id"]: view for view in load("conformal")["views"]}
        self.assertNotIn("slices", views["spring"])
        for vid in ("drain", "kerr_like"):
            marks = {m["view"]: m for m in views[vid]["slices"]}
            self.assertEqual(sorted(marks), ["funnel", "plane"])
            first = marks["plane"]["lines"][0][0]
            self.assertAlmostEqual(first[0], 0.0, places=3)
            self.assertAlmostEqual(first[1], math.pi / 2, places=3)
            # The laboratory's t = 0 at r = 4, the embedding's edge: v = 4 - ln 5 and UV = -(3/5) e^8.
            v = 4 - math.log(5)
            p, q = math.atan(-0.6 * math.exp(8 - v)), math.atan(math.exp(v))
            last = marks["plane"]["lines"][-1][-1]
            self.assertAlmostEqual(last[0], q - p, places=3)
            self.assertAlmostEqual(last[1], q + p, places=3)
            for line in marks["funnel"]["lines"]:
                for X, T in line:
                    self.assertAlmostEqual(T, 0.0, places=3)
            # T = 0 from the horizon, the corner, out to r = 5, where r* = 5 + ln(2/3)/2.
            self.assertAlmostEqual(marks["funnel"]["lines"][0][0][0], 0.0, places=3)
            self.assertAlmostEqual(marks["funnel"]["lines"][-1][-1][0],
                                   2 * math.atan(math.exp(5 + math.log(2 / 3) / 2)), places=3)


class Relations(unittest.TestCase):
    def test_each_relation_is_answered_with_the_kind_that_goes_with_it(self):
        answers = {"family": "family", "special_case": "generalisation"}
        own = load("metrics")["related"]
        self.assertEqual([r["id"] for r in own],
                         ["kerr", "btz", "schwarzschild", "natario", "spinning_string", "minkowski",
                          "unruh_acoustic_hole"])
        for relation in own:
            other = load("metrics", relation["id"])
            (back,) = [r for r in other["related"] if r["id"] == "draining_bathtub"]
            self.assertEqual(back["kind"], answers[relation["kind"]], relation["id"])

    def test_the_history_cites_every_reference_it_lists_before_the_relations_own(self):
        metric = load("metrics")
        history = metric["history"]
        cited = [key for key in metric["references"] if f"[{key}" in history or f", {key}" in history]
        self.assertEqual(cited, metric["references"][:len(cited)])
        self.assertEqual(metric["references"][len(cited):], ["banados1992", "natario2002"])

    def test_the_gain_quoted_is_the_published_one(self):
        """The abstract in Nature Physics gives 14% +- 8%; the preprint's 20% is not quoted."""
        history = load("metrics")["history"]
        self.assertIn("$14\\% \\pm 8\\%$", history)
        self.assertNotIn("20\\%", history)
        self.assertNotIn("20 percent", history)


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class Published(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import sympy as sp
        import verify_metrics as vm
        cls.sp, cls.vm = sp, vm
        cls.metric = load("metrics")
        cls.charts = {c["id"]: c for c in cls.metric["coordinates"]}

    def reader(self, system):
        chart = self.charts[system]
        return self.vm.Reader(chart["coords"], [p["symbol"] for p in chart["parameters"]], ())

    def matrix(self, system, field="metric_components"):
        reader = self.reader(system)
        coords = self.charts[system]["coords"]
        g = self.sp.zeros(len(coords), len(coords))
        for entry in self.charts[system][field]:
            i, j = (coords.index(x) for x in entry["indices"])
            g[i, j] = reader(entry["value"])
        return reader, g

    def names(self, reader, system):
        r = reader.symbol["r"]
        return r, reader.parameters["A"], reader.parameters["B"], reader.c

    def test_the_charts_are_the_three_the_literature_uses(self):
        self.assertEqual(list(self.charts), CHARTS)
        self.assertEqual([c["coords"] for c in self.metric["coordinates"]],
                         [["t", "r", "\\theta"], ["T", "r", "\\phi"], ["t", "r", "\\theta", "z"]])

    def test_the_laboratory_metric_is_vissers(self):
        """-(c^2 - (A^2 + B^2)/r^2) dt^2 - 2(A/r) dr dt - 2B dtheta dt + dr^2 + r^2 dtheta^2, in x^0 = ct."""
        sp = self.sp
        for system in ("laboratory", "vortex_filament"):
            reader, g = self.matrix(system)
            r, A, B, c = self.names(reader, system)
            want = sp.zeros(*g.shape)
            want[0, 0] = -(1 - (A ** 2 + B ** 2) / (c ** 2 * r ** 2))
            want[0, 1] = want[1, 0] = -A / (c * r)
            want[0, 2] = want[2, 0] = -B / c
            want[1, 1], want[2, 2] = 1, r ** 2
            if system == "vortex_filament":
                want[3, 3] = 1
            self.assertEqual((g - want).applyfunc(sp.simplify), sp.zeros(*g.shape), system)

    def test_sound_moves_at_c_past_the_water_in_every_direction(self):
        sp = self.sp
        reader, g = self.matrix("laboratory")
        r, A, B, c = self.names(reader, "laboratory")
        alpha = sp.Symbol("alpha", real=True)
        ray = sp.Matrix([1, A / (c * r) + sp.cos(alpha), (B / (c * r) + sp.sin(alpha)) / r])
        self.assertEqual(sp.simplify((ray.T * g * ray)[0, 0]), 0)

    def test_the_horizon_and_the_ergosurface_are_where_visser_puts_them(self):
        sp = self.sp
        reader, g = self.matrix("laboratory")
        _, ginv = self.matrix("laboratory", "inverse_metric_components")
        r, A, B, c = self.names(reader, "laboratory")
        self.assertEqual(sp.simplify(g * ginv - sp.eye(3)), sp.zeros(3, 3))
        self.assertEqual(sp.simplify(ginv[1, 1] - (1 - A ** 2 / (c ** 2 * r ** 2))), 0)
        self.assertEqual(sp.simplify(g[0, 0].subs(r, sp.sqrt(A ** 2 + B ** 2) / c)), 0)
        self.assertEqual(sp.simplify(ginv[0, 0] + 1), 0)

    def test_the_surface_gravity_is_c_cubed_over_the_strength_of_the_drain(self):
        """Visser's g_H = (1/2) d(c^2 - v_r^2)/dr on the horizon, with v_r = A/r read from g^rr."""
        sp = self.sp
        reader, ginv = self.matrix("laboratory", "inverse_metric_components")
        r, A, B, c = self.names(reader, "laboratory")
        a = sp.Symbol("a", positive=True)
        gravity = (sp.diff(c ** 2 * ginv[1, 1], r) / 2).subs(A, -a).subs(r, a / c)
        self.assertEqual(sp.simplify(gravity - c ** 3 / a), 0)

    def test_the_kerr_like_chart_is_the_laboratory_chart_pulled_back(self):
        """c dt = c dT - A c r dr/(c^2r^2 - A^2) and dtheta = dphi - A B dr/(r (c^2r^2 - A^2))."""
        sp = self.sp
        lab_reader, lab = self.matrix("laboratory")
        reader, g = self.matrix("kerr_like")
        r, A, B, c = self.names(reader, "kerr_like")
        lr, lA, lB, lc = self.names(lab_reader, "laboratory")
        lab = lab.subs({lr: r, lA: A, lB: B, lc: c}, simultaneous=True)
        J = sp.Matrix([[1, -A * c * r / (c ** 2 * r ** 2 - A ** 2), 0], [0, 1, 0],
                       [0, -A * B / (r * (c ** 2 * r ** 2 - A ** 2)), 1]])
        self.assertEqual((J.T * lab * J - g).applyfunc(sp.simplify), sp.zeros(3, 3))
        self.assertEqual((g[0, 1], g[1, 2]), (0, 0))

    def test_the_curvature_depends_on_the_speed_of_the_water_alone(self):
        sp = self.sp
        for system in CHARTS:
            reader = self.reader(system)
            r, A, B, c = self.names(reader, system)
            chart = self.charts[system]
            scalar = reader(chart["ricci_scalar"].partition("=")[2])
            kretschmann = reader(chart["kretschmann"].partition("=")[2])
            speed2 = (A ** 2 + B ** 2) / r ** 2
            self.assertEqual(sp.simplify(scalar - 2 * speed2 / (c ** 2 * r ** 2)), 0, system)
            self.assertEqual(sp.simplify(kretschmann - 44 * speed2 ** 2 / (c ** 4 * r ** 4)), 0, system)

    def test_with_the_water_at_rest_every_curvature_component_vanishes(self):
        for system in CHARTS:
            reader = self.reader(system)
            A, B = reader.parameters["A"], reader.parameters["B"]
            for entry in self.charts[system]["riemann"]["variants"]["ulll"]["nonzero"]:
                self.assertEqual(self.sp.simplify(reader(entry["value"]).subs({A: 0, B: 0})), 0, system)

    def test_the_plane_of_the_drain_has_no_weyl_tensor_and_the_filament_has_one(self):
        for system, empty in (("laboratory", True), ("kerr_like", True), ("vortex_filament", False)):
            self.assertEqual(self.charts[system]["weyl_tensor"]["variants"]["ulll"]["nonzero"] == [], empty, system)


if __name__ == "__main__":
    unittest.main()
