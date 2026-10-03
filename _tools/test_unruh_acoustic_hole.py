"""Unruh's acoustic black hole, a fluid of constant density falling into a point at c r_0^2/r^2:
python3 -m unittest discover -s _tools

What its texts and drawings state, held on the published files. The first class needs nothing
but the files: the horizon where every spacetime diagram marks it, the rays of each chart against
the speed of sound through the falling fluid, the two moments of the embedding diagram, and the
conformal diagram's edges. The second holds the relations to the answers of the other files and
the History to the numbers of its sources. The third reads the published tensors through the
checker's Reader and holds them to Visser's line element, to sound moving at c past the fluid, to
the map between the charts, to Unruh's temperature and to the curvature the captions quote; it
needs sympy and is skipped where it is absent.

Every drawing is at r_0 = c = 1, so the horizon is the unit of length and the fluid falls at 1/r^2.
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

ID = "unruh_acoustic_hole"
CHARTS = ["laboratory", "unruh"]


def load(folder, metric_id=ID):
    return json.loads((DATA / folder / f"{metric_id}.json").read_text(encoding="utf-8"))


def flat_views():
    return {(system, view["id"]): view for system, views in load("diagrams")["systems"].items() for view in views}


def chart_point(view, point):
    """A point of a flat view's unit square as (time, r) of its chart."""
    X0, X1, Y0, Y1 = view["box"]
    return Y0 + point[1] * (Y1 - Y0), X0 + point[0] * (X1 - X0)


def rstar(r):
    """The tortoise coordinate at r_0 = 1, the integral of r^4/(r^4 - 1), vanishing at r = 0."""
    return r + math.log(abs((r - 1) / (r + 1))) / 4 - math.atan(r) / 2


def unruh_shift(r):
    """Unruh's tau minus the laboratory's t at r_0 = c = 1."""
    return -math.log(abs((r - 1) / (r + 1))) / 4 - math.atan(r) / 2


class Drawings(unittest.TestCase):
    def marked(self, view, kind):
        X0, X1 = view["box"][:2]
        return sorted(X0 + line[0][0] * (X1 - X0) for m in view["markers"] if m["kind"] == kind for line in m["lines"])

    def test_each_chart_has_one_view_with_the_horizon_at_one(self):
        views = flat_views()
        self.assertEqual(sorted(views), [("laboratory", "infall"), ("unruh", "exterior")])
        self.assertEqual(len(self.marked(views[("laboratory", "infall")], "grr")), 1)
        self.assertAlmostEqual(self.marked(views[("laboratory", "infall")], "grr")[0], 1.0, delta=2e-3)

    def test_the_laboratory_ends_on_a_singular_sink_and_unruhs_time_on_its_horizon(self):
        for key, view in flat_views().items():
            singular = [m for m in view["markers"] if m["kind"] == "singular"]
            if key[0] == "unruh":
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
        """What each family of rays keeps: the laboratory's rays have dr/dt = -1/r^2 -+ 1, and
        Unruh's time c dtau/dr = -+ r^4/(r^4 - 1)."""
        if key[0] == "unruh":
            return (lambda T, r: T + rstar(r), lambda T, r: T - rstar(r))
        return (lambda t, r: t + r - math.atan(r), lambda t, r: t - r - math.log(abs((r - 1) / (r + 1))) / 2)

    def test_every_ray_is_sound_moving_at_c_through_the_falling_fluid(self):
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

    def test_inside_the_horizon_every_cone_points_to_the_sink(self):
        view = flat_views()[("laboratory", "infall")]
        inside = [c for c in view["cones"] if chart_point(view, c["at"])[1] < 0.95]
        self.assertGreater(len(inside), 6)
        for cone in inside:
            for edge in (cone["a"], cone["b"]):
                self.assertLess(edge[0], 0)
                self.assertGreater(edge[1], 0)

    def test_both_views_draw_both_moments(self):
        views = flat_views()
        lab_view, un_view = views[("laboratory", "infall")], views[("unruh", "exterior")]
        lab = {m["view"]: m for m in lab_view["slices"]}
        self.assertEqual(sorted(lab), ["funnel", "space"])
        for t, r in (chart_point(lab_view, p) for line in lab["space"]["lines"] for p in line):
            self.assertAlmostEqual(t, 0.0, delta=1e-3)
        for t, r in (chart_point(lab_view, p) for line in lab["funnel"]["lines"] for p in line):
            if r > 1.02:
                self.assertAlmostEqual(t, -unruh_shift(r), delta=5e-3)
        un = {m["view"]: m for m in un_view["slices"]}
        self.assertEqual(sorted(un), ["funnel", "space"])
        for T, r in (chart_point(un_view, p) for line in un["funnel"]["lines"] for p in line):
            self.assertAlmostEqual(T, 0.0, delta=1e-3)
        for T, r in (chart_point(un_view, p) for line in un["space"]["lines"] for p in line):
            if r > 1.02:
                self.assertAlmostEqual(T, unruh_shift(r), delta=5e-3)

    def test_the_embedding_is_a_flat_plane_in_the_laboratory_and_an_elliptic_funnel_in_unruhs_time(self):
        views = {view["id"]: view for view in load("embedding")["views"]}
        self.assertEqual(sorted(views), ["funnel", "space"])
        (space,) = views["space"]["surfaces"][0]["pieces"]
        self.assertEqual((space["system"], space["coordinate"]), ("laboratory", "r"))
        self.assertEqual((space["points"][0][0], space["points"][-1][0]), (0.0, 4.0))
        for r, rho, z in space["points"]:
            self.assertAlmostEqual(rho, r, places=6)
            self.assertEqual(z, 0.0)
        (funnel,) = views["funnel"]["surfaces"][0]["pieces"]
        self.assertEqual((funnel["system"], funnel["coordinate"]), ("unruh", "r"))
        self.assertEqual((funnel["points"][0][0], funnel["points"][-1][0]), (1.0, 5.0))
        previous = -1.0
        for r, rho, z in funnel["points"]:
            self.assertAlmostEqual(rho, r, places=6)
            # z = int_1^r ds/sqrt(s^4 - 1), by Simpson's rule in s = 1 + u^2, which takes the root away.
            n, top = 2000, math.sqrt(r - 1)
            f = lambda u: 2 / math.sqrt((2 + u * u) * (1 + (1 + u * u) ** 2))
            h = top / n
            want = h / 3 * sum((1 if k in (0, n) else 4 if k % 2 else 2) * f(k * h) for k in range(n + 1))
            self.assertAlmostEqual(z, want, places=4)
            self.assertGreaterEqual(z, previous)
            previous = z
        # The funnel levels off below K(1/2)/sqrt 2, about 1.311.
        self.assertLess(previous, 1.3111)

    def test_the_conformal_diagram_is_the_outside_and_the_hole_with_the_laboratorys_time_ending_on_an_edge(self):
        views = {view["id"]: view for view in load("conformal")["views"]}
        self.assertEqual(list(views), ["laboratory", "unruh"])
        half = math.pi / 2
        for vid in views:
            lines = {}
            for layer in views[vid]["layers"]:
                if layer["class"] in ("horizon", "singular", "chartedge"):
                    lines[layer["class"]] = layer["points"]
            for (X, T), want in zip(lines["horizon"], ((0, 0), (half, half))):
                self.assertAlmostEqual(X, want[0], places=3)
                self.assertAlmostEqual(T, want[1], places=3)
            for (X, T), want in zip(lines["chartedge"], ((half, -half), (-half, half))):
                self.assertAlmostEqual(X, want[0], places=3)
                self.assertAlmostEqual(T, want[1], places=3)
            self.assertTrue(all(abs(T - half) < 0.06 for _, T in lines["singular"]), vid)

    def test_the_moments_on_the_conformal_diagram_run_from_the_sink_and_from_the_horizon(self):
        views = {view["id"]: view for view in load("conformal")["views"]}
        for vid in views:
            marks = {m["view"]: m for m in views[vid]["slices"]}
            self.assertEqual(sorted(marks), ["funnel", "space"])
            first = marks["space"]["lines"][0][0]
            self.assertAlmostEqual(first[0], 0.0, places=3)
            self.assertAlmostEqual(first[1], math.pi / 2, places=3)
            # The laboratory's t = 0 at r = 4, the embedding's edge: v = 4 - arctan 4 and
            # UV = -(3/5) e^(16 - 2 arctan 4).
            v = 4 - math.atan(4)
            p = math.atan(-0.6 * math.exp(16 - 2 * math.atan(4) - 2 * v))
            q = math.atan(math.exp(2 * v))
            last = marks["space"]["lines"][-1][-1]
            self.assertAlmostEqual(last[0], q - p, places=3)
            self.assertAlmostEqual(last[1], q + p, places=3)
            for line in marks["funnel"]["lines"]:
                for X, T in line:
                    self.assertAlmostEqual(T, 0.0, places=3)
            # tau = 0 from the horizon, the corner, out to r = 5.
            self.assertAlmostEqual(marks["funnel"]["lines"][0][0][0], 0.0, places=3)
            self.assertAlmostEqual(marks["funnel"]["lines"][-1][-1][0], 2 * math.atan(math.exp(2 * rstar(5))), places=3)


class Relations(unittest.TestCase):
    def test_each_relation_is_answered_with_the_kind_that_goes_with_it(self):
        answers = {"family": "family", "special_case": "generalisation"}
        own = load("metrics")["related"]
        self.assertEqual([r["id"] for r in own], ["draining_bathtub", "schwarzschild", "natario", "alcubierre", "minkowski"])
        for relation in own:
            other = load("metrics", relation["id"])
            (back,) = [r for r in other["related"] if r["id"] == ID]
            self.assertEqual(back["kind"], answers[relation["kind"]], relation["id"])

    def test_the_history_cites_every_reference_it_lists_before_the_relations_own(self):
        metric = load("metrics")
        history = metric["history"]
        cited = [key for key in metric["references"] if f"[{key}" in history or f", {key}" in history]
        self.assertEqual(cited, metric["references"][:len(cited)])
        self.assertEqual(metric["references"][len(cited):], ["natario2002", "alcubierre1994"])

    def test_the_numbers_quoted_are_unruhs_and_steinhauers(self):
        """Unruh's estimate, 3e-7 K at c = 300 m/s and R = 1 mm, against T = (hbar/2 pi k) c/R; and
        Steinhauer's flow, subsonic outside the step and supersonic inside it."""
        hbar, k = 1.054571817e-34, 1.380649e-23
        self.assertAlmostEqual(hbar * 300 / 1e-3 / (2 * math.pi * k) / 1e-7, 3.6, delta=0.1)
        history = load("metrics")["history"]
        self.assertIn("$3\\times10^{-7}$ K", history)
        for value in ("$0.24$", "$0.57$", "$1.02$", "$0.25$"):
            self.assertIn(value + " mm/s", history)
        self.assertLess(0.24, 0.57)
        self.assertGreater(1.02, 0.25)


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

    def names(self, reader):
        return reader.symbol["r"], reader.symbol["\\theta"], reader.parameters["r_0"], reader.c

    def test_the_charts_are_the_two_of_unruh_and_visser(self):
        self.assertEqual(list(self.charts), CHARTS)
        self.assertEqual([c["coords"] for c in self.metric["coordinates"]],
                         [["t", "r", "\\theta", "\\phi"], ["\\tau", "r", "\\theta", "\\phi"]])

    def test_the_laboratory_metric_is_vissers_canonical_hole(self):
        """-c^2dt^2 + (dr + c r_0^2 dt/r^2)^2 + r^2 dOmega^2, in x^0 = ct."""
        sp = self.sp
        reader, g = self.matrix("laboratory")
        r, th, r0, c = self.names(reader)
        want = sp.diag(-(1 - r0 ** 4 / r ** 4), 1, r ** 2, r ** 2 * sp.sin(th) ** 2)
        want[0, 1] = want[1, 0] = r0 ** 2 / r ** 2
        self.assertEqual((g - want).applyfunc(sp.simplify), sp.zeros(4, 4))

    def test_sound_moves_at_c_past_the_fluid_in_every_direction(self):
        sp = self.sp
        reader, g = self.matrix("laboratory")
        r, th, r0, c = self.names(reader)
        a, b = sp.symbols("alpha beta", real=True)
        ray = sp.Matrix([1, -r0 ** 2 / r ** 2 + sp.cos(a), sp.sin(a) * sp.cos(b) / r,
                         sp.sin(a) * sp.sin(b) / (r * sp.sin(th))])
        self.assertEqual(sp.simplify((ray.T * g * ray)[0, 0]), 0)

    def test_the_horizon_is_where_the_fluid_falls_at_c(self):
        sp = self.sp
        reader, g = self.matrix("laboratory")
        _, ginv = self.matrix("laboratory", "inverse_metric_components")
        r, th, r0, c = self.names(reader)
        self.assertEqual(sp.simplify(g * ginv - sp.eye(4)), sp.zeros(4, 4))
        self.assertEqual(sp.simplify(ginv[1, 1] - (1 - r0 ** 4 / r ** 4)), 0)
        self.assertEqual(sp.simplify(ginv[1, 1].subs(r, r0)), 0)
        self.assertEqual(sp.simplify(g[0, 0].subs(r, r0)), 0)

    def test_unruhs_temperature_is_hbar_c_over_pi_k_r0(self):
        """Unruh's T = (hbar/2 pi k) dv/dr on the horizon, with v = -c r_0^2/r^2 read from g^rr =
        1 - v^2/c^2, is hbar c/(pi k r_0), and Visser's surface gravity (1/2) d(c^2 - v^2)/dr is
        2 c^2/r_0."""
        sp = self.sp
        reader, ginv = self.matrix("laboratory", "inverse_metric_components")
        r, th, r0, c = self.names(reader)
        v = -c * sp.sqrt(sp.simplify(1 - ginv[1, 1]))
        self.assertEqual(sp.simplify(v + c * r0 ** 2 / r ** 2), 0)
        hbar, k = sp.symbols("hbar k", positive=True)
        T = (hbar / (2 * sp.pi * k) * sp.diff(v, r)).subs(r, r0)
        self.assertEqual(sp.simplify(T - hbar * c / (sp.pi * k * r0)), 0)
        self.assertEqual(sp.simplify((sp.diff(c ** 2 * ginv[1, 1], r) / 2).subs(r, r0) - 2 * c ** 2 / r0), 0)

    def test_unruhs_time_is_the_laboratory_chart_pulled_back(self):
        """c dtau = c dt - r_0^2 r^2 dr/(r^4 - r_0^4), Unruh's tau = t + int v dr/(c^2 - v^2)."""
        sp = self.sp
        lab_reader, lab = self.matrix("laboratory")
        reader, g = self.matrix("unruh")
        r, th, r0, c = self.names(reader)
        lr, lth, lr0, lc = self.names(lab_reader)
        lab = lab.subs({lr: r, lth: th, lr0: r0, lc: c}, simultaneous=True)
        J = sp.eye(4)
        J[0, 1] = r0 ** 2 * r ** 2 / (r ** 4 - r0 ** 4)
        self.assertEqual((J.T * lab * J - g).applyfunc(sp.simplify), sp.zeros(4, 4))
        v = -c * r0 ** 2 / r ** 2
        self.assertEqual(sp.simplify(J[0, 1] - (-v * c / (c ** 2 - v ** 2))), 0)
        self.assertEqual(g[0, 1], 0)

    def test_the_curvature_quoted_in_every_chart(self):
        sp = self.sp
        for system in CHARTS:
            reader = self.reader(system)
            r, th, r0, c = self.names(reader)
            chart = self.charts[system]
            scalar = reader(chart["ricci_scalar"].partition("=")[2])
            kretschmann = reader(chart["kretschmann"].partition("=")[2])
            self.assertEqual(sp.simplify(scalar - 6 * r0 ** 4 / r ** 6), 0, system)
            self.assertEqual(sp.simplify(kretschmann - 468 * r0 ** 8 / r ** 12), 0, system)

    def test_with_the_fluid_at_rest_every_curvature_component_vanishes(self):
        for system in CHARTS:
            reader = self.reader(system)
            r0 = reader.parameters["r_0"]
            for entry in self.charts[system]["riemann"]["variants"]["ulll"]["nonzero"]:
                self.assertEqual(self.sp.simplify(reader(entry["value"]).subs(r0, 0)), 0, system)


if __name__ == "__main__":
    unittest.main()
