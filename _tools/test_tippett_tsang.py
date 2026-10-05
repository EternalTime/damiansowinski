"""Tippett and Tsang's time machine, a bubble whose riders go round a closed timelike curve:
.venv.noindex/bin/python -m unittest discover -s _tools

What the texts and drawings of `tippett_tsang` state, held on the published files. The bubble the
diagrams declare is held to Tippett and Tsang's own numbers. The published components are held to
Minkowski's metric at h = 0 and to the published interior chart at h = 1, to the determinant
-(1 - 4h(1 - h)c^2t^2/(x^2 + c^2t^2)) and to vanishing on the plane of t and x where x = 0 and
h = 1/2, to the published Rindler chart in the polar chart at h = 1, to the plane y = z = 0 holding
its light rays, and to the acceleration c^2/xi of the closed timelike curves. The spacetime diagram
is held to cones that point up in t inside the right half of the ring and down inside the left, to
the four singular events, and to rays the wall turns round in time; the figure in three dimensions
to its closed curve and its cones; and the height over the plane to the angle of the published
g_tx and g_xx. The tests of the published components need sympy and are skipped where it is absent.
"""
import importlib.util
import json
import math
import sys
import unittest
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
HAS_SYMPY = importlib.util.find_spec("sympy") is not None
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))

# Tippett and Tsang's bubble in units of A, the radius of the circle its centre runs round.
A, R, ALPHA = 1.0, 0.7, 50 / 3
INNER, OUTER = math.sqrt(A ** 2 - R ** 2), math.sqrt(A ** 2 + R ** 2)


def top_hat(t, x, y=0.0, z=0.0):
    return (1 + math.tanh(ALPHA * (R ** 4 - y ** 4 - z ** 4 - (x * x + t * t - A * A) ** 2))) / 2


def published(system):
    metric = json.loads((DATA / "metrics" / "tippett_tsang.json").read_text(encoding="utf-8"))
    return next(c for c in metric["coordinates"] if c["id"] == system)


def flat_view(system, view):
    diagrams = json.loads((DATA / "diagrams" / "tippett_tsang.json").read_text(encoding="utf-8"))
    return next(v for v in diagrams["systems"][system] if v["id"] == view)


class DeclaredBubble(unittest.TestCase):
    def test_the_bubble_drawn_is_tippett_and_tsangs_own(self):
        """Their A = 100, R = 70 and alpha = 1/6000000, in units of A."""
        a, r, alpha = Fraction(100), Fraction(70), Fraction(1, 6000000)
        self.assertEqual(r / a, Fraction(7, 10))
        self.assertEqual(alpha * a ** 4, Fraction(50, 3))
        self.assertEqual((r / a) ** 4, Fraction(2401, 10000))
        source = (ROOT / "_tools" / "derivations" / "null_rays.py").read_text(encoding="utf-8")
        self.assertIn('_TT_BOX = "Rational(2401, 10000) - y**4 - z**4 - ({ring} - 1)**2"', source)
        self.assertIn('TT_H = "(1 + tanh(Rational(50, 3)*(" + _TT_BOX.format(ring="x**2 + t**2") + ")))/2"', source)

    def test_the_box_is_all_but_flat_at_its_centre_and_the_walls_stand_where_the_captions_say(self):
        self.assertGreater(top_hat(0.0, A), 0.999)
        for edge in (INNER, OUTER):
            self.assertAlmostEqual(top_hat(0.0, edge), 0.5, places=12)
            self.assertAlmostEqual(top_hat(edge, 0.0), 0.5, places=12)
        self.assertAlmostEqual(top_hat(0.0, A, R), 0.5, places=12)


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class Components(unittest.TestCase):
    """The published metrics, read as the checker reads them, in the chart x^0 = ct with c = 1."""

    @classmethod
    def setUpClass(cls):
        import sympy as sp
        import null_rays as nr
        cls.sp = sp
        cls.charts = {}
        for system in ("cartesian", "polar", "interior", "rindler"):
            _, entry, reader = nr.load("tippett_tsang", system)
            g = nr.published_matrix(reader, entry, "metric_components").subs(reader.c, 1)
            cls.charts[system] = (entry, reader, g)

    def test_with_no_bubble_the_metric_is_minkowskis(self):
        sp = self.sp
        _, reader, g = self.charts["cartesian"]
        flat = g.subs(reader.parameters["h"], 0).doit().applyfunc(sp.simplify)
        self.assertEqual(flat, sp.diag(-1, 1, 1, 1))

    def test_inside_the_bubble_it_is_the_interior_chart_and_rindlers(self):
        sp = self.sp
        _, reader, g = self.charts["cartesian"]
        _, inner, g_in = self.charts["interior"]
        same = dict(zip((inner.symbol[c] for c in ("t", "x", "y", "z")), (reader.symbol[c] for c in ("t", "x", "y", "z"))))
        difference = g.subs(reader.parameters["h"], 1).doit() - g_in.subs(same)
        self.assertEqual(difference.applyfunc(sp.simplify), sp.zeros(4, 4))
        _, polar, g_p = self.charts["polar"]
        _, rindler, g_r = self.charts["rindler"]
        same = dict(zip((rindler.symbol[c] for c in ("\\lambda", "\\xi")), (polar.symbol[c] for c in ("\\lambda", "\\xi"))))
        difference = g_p.subs(polar.parameters["h"], 1).doit() - g_r.subs(same)
        self.assertEqual(difference.applyfunc(sp.simplify), sp.zeros(4, 4))

    def test_the_metric_is_degenerate_exactly_where_the_wall_crosses_x_equal_zero_at_one_half(self):
        sp = self.sp
        _, reader, g = self.charts["cartesian"]
        t, x, h = reader.symbol["t"], reader.symbol["x"], reader.parameters["h"]
        self.assertEqual(sp.simplify(g.det() + 1 - 4 * h * (1 - h) * t ** 2 / (x ** 2 + t ** 2)), 0)
        # 1 - 4h(1 - h)s with 0 <= s <= 1 vanishes only at h = 1/2 and s = 1, which is x = 0.
        at = {h: sp.Rational(1, 2), x: 0}
        for i in range(2):
            for j in range(2):
                self.assertEqual(sp.simplify(g[i, j].subs(at)), 0, (i, j))

    def test_the_plane_of_t_and_x_holds_its_light_rays(self):
        """Every Christoffel symbol that could turn a ray of the plane y = z = 0 out of it carries a
        first derivative of h along y or z, which vanishes there for a bubble even in y and z."""
        sp = self.sp
        entry, reader, _ = self.charts["cartesian"]
        h = reader.parameters["h"]
        y, z = reader.symbol["y"], reader.symbol["z"]
        seen = 0
        for item in entry["christoffel"]["variants"]["ull"]["nonzero"]:
            up, a, b = item["indices"]
            if up in ("y", "z") and a in ("t", "x") and b in ("t", "x"):
                value = reader(item["value"])
                leaves = value.subs({sp.Derivative(h, y): 0, sp.Derivative(h, z): 0})
                self.assertEqual(sp.simplify(leaves), 0, item["indices"])
                seen += 1
        self.assertEqual(seen, 8)

    def test_a_rider_on_a_closed_timelike_curve_feels_one_over_xi(self):
        sp = self.sp
        entry, reader, g = self.charts["rindler"]
        xi = reader.symbol["\\xi"]
        gamma = {tuple(i["indices"]): reader(i["value"]) for i in entry["christoffel"]["variants"]["ull"]["nonzero"]}
        self.assertEqual(g[0, 0], -xi ** 2)
        # u = d/d lambda / xi is a unit vector, and its acceleration is Gamma^xi_lambda lambda / xi^2.
        self.assertEqual(sp.simplify(gamma[("\\xi", "\\lambda", "\\lambda")] / xi ** 2 - 1 / xi), 0)
        self.assertIn("\\lambda \\sim \\lambda + 2\\pi \\;\\text{(closed timelike curves)}", entry["domains"])


class SpacetimeDiagram(unittest.TestCase):
    """The plane of t and x through the bubble, as the published drawing holds it."""

    @classmethod
    def setUpClass(cls):
        cls.view = flat_view("cartesian", "tx")
        x0, x1, y0, y1 = cls.view["box"]
        cls.chart = lambda self, u: (x0 + (x1 - x0) * u[0], y0 + (y1 - y0) * u[1])

    def test_inside_the_ring_the_future_runs_counterclockwise_and_outside_it_upward(self):
        inside = outside = 0
        for cone in self.view["cones"]:
            x, t = self.chart(cone["at"])
            xi = math.hypot(x, t)
            axis = [cone["a"][0] + cone["b"][0], cone["a"][1] + cone["b"][1]]
            if INNER + 0.1 < xi < OUTER - 0.1:
                # Counterclockwise is the direction (-t, x) of the drawing's (x, t).
                self.assertGreater(-t * axis[0] + x * axis[1], 0, (x, t))
                inside += 1
            elif xi < INNER - 0.1 or xi > OUTER + 0.1:
                self.assertGreater(axis[1], 0, (x, t))
                self.assertAlmostEqual(axis[0], 0, places=3)
                outside += 1
        self.assertEqual(inside, 12)
        self.assertGreaterEqual(outside, 24)

    def test_on_the_left_the_cones_inside_point_down_beside_cones_outside_that_point_up(self):
        down = [c for c in self.view["cones"] if abs(self.chart(c["at"])[0] + 1.0) < 1e-9
                and abs(self.chart(c["at"])[1]) < 0.3]
        self.assertEqual(len(down), 2)
        for cone in down:
            self.assertLess(cone["a"][1], 0)
            self.assertLess(cone["b"][1], 0)
        up = [c for c in self.view["cones"] if abs(self.chart(c["at"])[0] - 1.0) < 1e-9
              and abs(self.chart(c["at"])[1]) < 0.3]
        for cone in up:
            self.assertGreater(cone["a"][1], 0)
            self.assertGreater(cone["b"][1], 0)

    def test_the_four_singular_events_stand_where_the_walls_cross_x_equal_zero(self):
        marked = sorted(self.chart(p)[1] for m in self.view["markers"] if m["kind"] == "removed" for p in m["points"])
        for got, want in zip(marked, (-OUTER, -INNER, INNER, OUTER), strict=True):
            self.assertAlmostEqual(got, want, places=3)
        for m in self.view["markers"]:
            if m["kind"] == "removed":
                self.assertAlmostEqual(self.chart(m["points"][0])[0], 0.0, places=9)

    def test_the_wall_turns_some_rays_round_in_time(self):
        """Outside the bubble the family named moving left keeps x + ct, so a ray of it that enters
        through the bottom edge leaves through the left edge or the top. One that leaves through the
        right edge, or comes back to the edge it entered by, has been turned round by the wall."""
        def edge(p):
            if p[1] < 0.01:
                return "bottom"
            if p[1] > 0.99:
                return "top"
            return "left" if p[0] < 0.01 else "right" if p[0] > 0.99 else None
        ends = [{edge(ray[0]), edge(ray[-1])} for ray in self.view["rays"]["P"]]
        self.assertIn({"bottom", "right"}, ends)
        self.assertIn({"left"}, ends)
        self.assertIn({"right"}, ends)
        self.assertIn({"bottom", "top"}, ends)

    def test_the_moment_of_the_height_is_marked_at_half_a(self):
        (mark,) = self.view["slices"]
        self.assertEqual(mark["label"], "$ct = A/2$")
        for p in mark["lines"][0]:
            self.assertAlmostEqual(self.chart(p)[1], 0.5, places=3)


class Figure(unittest.TestCase):
    """The slice z = 0 of t, x and y in three dimensions."""

    @classmethod
    def setUpClass(cls):
        diagrams = json.loads((DATA / "diagrams" / "tippett_tsang.json").read_text(encoding="utf-8"))
        (cls.figure,) = diagrams["projections"]["cartesian"]

    def test_the_closed_curve_is_the_circle_of_radius_a_in_the_plane_of_t_and_x(self):
        (curve,) = [line for line in self.figure["turn"]["lines"] if line["class"] == "ctc"]
        self.assertEqual(curve["points"][0], curve["points"][-1])
        for X, Y, T in curve["points"]:
            self.assertAlmostEqual(math.hypot(X, T), A, places=5)
            self.assertEqual(Y, 0)

    def test_cones_round_the_curve_lean_along_it_and_cones_outside_stand_upright(self):
        cones = self.figure["turn"]["cones"]
        self.assertEqual(len(cones), 16)
        on_curve = 0
        for cone in cones:
            X, Y, T = cone["apex"]
            n = len(cone["rim"])
            axis = [sum(p[i] for p in cone["rim"]) / n - cone["apex"][i] for i in range(3)]
            if abs(math.hypot(X, T) - A) < 1e-6:
                # d/d lambda at (x, ct) is (-ct, x): the axis points counterclockwise round the circle.
                self.assertGreater(-T * axis[0] + X * axis[2], 0)
                self.assertAlmostEqual(X * axis[0] + T * axis[2], 0, delta=0.01)
                on_curve += 1
            else:
                self.assertGreater(axis[2], 0)
                self.assertAlmostEqual(axis[0], 0, places=5)
            self.assertAlmostEqual(axis[1], 0, places=5)
        self.assertEqual(on_curve, 8)

    def test_the_wall_is_drawn_where_the_top_hat_is_one_half(self):
        walls = [line for line in self.figure["turn"]["lines"] if line["class"] == "edge"]
        self.assertEqual(len(walls), 16)
        for line in walls:
            for X, Y, T in line["points"]:
                self.assertAlmostEqual(top_hat(T, X, Y), 0.5, delta=2e-4)


class Height(unittest.TestCase):
    """The turn of the light cone as a height over the plane z = 0 at ct = A/2."""

    @classmethod
    def setUpClass(cls):
        embedding = json.loads((DATA / "embedding" / "tippett_tsang.json").read_text(encoding="utf-8"))
        (cls.view,) = embedding["views"]
        (cls.piece,) = cls.view["surfaces"][0]["pieces"]

    @staticmethod
    def turn(x, y, t=0.5):
        h = top_hat(t, x, y)
        r2 = x * x + t * t
        return 0.5 * math.atan2(2 * h * x * t / r2, 1 - 2 * h * t * t / r2)

    def nodes(self):
        grid = self.piece["grid"]
        for i, x in enumerate(grid["u"]):
            for j, y in enumerate(grid["v"]):
                yield x, y, grid["z"][i][j]

    def test_the_height_is_half_the_angle_of_g_tx_and_g_xx(self):
        count = 0
        for x, y, z in self.nodes():
            self.assertAlmostEqual(z, self.turn(x, y), places=5)
            count += 1
        self.assertGreater(count, 10000)

    def test_the_two_boxes_are_turned_by_thirty_degrees_opposite_ways(self):
        centre = math.sqrt(3) / 2
        self.assertAlmostEqual(self.turn(centre, 0.0), math.pi / 6, places=3)
        self.assertAlmostEqual(self.turn(-centre, 0.0), -math.pi / 6, places=3)
        self.assertEqual(self.turn(0.0, 0.0), 0.0)
        heights = [z for _, _, z in self.nodes()]
        self.assertAlmostEqual(max(heights), -min(heights), places=6)
        # Nearest the inner walls the angle lambda is nearest 90 degrees, so the turn is greatest there.
        self.assertGreater(max(heights), math.pi / 6)
        self.assertLess(max(heights), math.pi / 4)

    def test_the_moment_is_a_moment_of_space(self):
        """g_xx = 1 - 2h c^2t^2/(x^2 + c^2t^2) > 0 everywhere at ct = A/2, since c^2t^2 < (A^2 - R^2)/2."""
        self.assertLess(0.25, (A * A - R * R) / 2)
        for x, y, _ in self.nodes():
            self.assertGreater(1 - 2 * top_hat(0.5, x, y) * 0.25 / (x * x + 0.25), 0)


if __name__ == "__main__":
    unittest.main()
