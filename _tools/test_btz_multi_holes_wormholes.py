"""Many black holes and wormholes in three dimensions, anti-de Sitter space cut along totally
geodesic timelike surfaces and glued back together: .venv.noindex/bin/python -m unittest discover -s _tools

What its texts and drawings state, held on the published files and on the geometry of the
hyperbolic plane. The tent the drawings declare is the region between the four surfaces
X = +-alpha V and Y = +-alpha V of -U^2 - V^2 + X^2 + Y^2 = -l^2 at alpha = sqrt(2/3). The first
class holds the numbers every caption leans on to the group that does the gluing, as Lorentz
transformations of the plane's embedding: the wormhole's one horizon is four times the distance
between adjacent cuts, the three holes' horizons are that distance twice and its double, and a
fold leaves infinity at ct = pi l/6. The second holds the drawings to those numbers, and the
third the charts and the relations.
"""
import json
import math
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
ME = "btz_multi_holes_wormholes"
CHARTS = ["sausage", "stereographic", "free_fall", "exterior"]

ALPHA = math.sqrt(2.0 / 3.0)
B = math.atanh(ALPHA)                   # each cut's distance from the axis, in units of l
D = math.acosh(2.0)                     # the distance between adjacent cuts
T_P = math.pi / 6                       # the sausage time at which a fold leaves infinity
ETA = np.diag([-1.0, 1.0, 1.0])         # the metric of (V, X, Y)


def load(folder, metric_id=ME):
    return json.loads((DATA / folder / f"{metric_id}.json").read_text(encoding="utf-8"))


def cut(k):
    """The unit normal (n_V, n_X, n_Y) of the plane through the origin that holds cut k, the
    geodesic X cos(k pi/2) + Y sin(k pi/2) = alpha V of the hyperbolic plane."""
    turn = k * math.pi / 2
    return np.array([math.sinh(B), math.cosh(B) * math.cos(turn), math.cosh(B) * math.sin(turn)])


def mirror(n):
    """The reflection of the hyperbolic plane in the geodesic with unit normal n."""
    n = np.asarray(n, dtype=float)
    return np.eye(3) - 2 * np.outer(n, ETA @ n)


def boost(axis, rapidity):
    """The translation of the hyperbolic plane along the X axis (1) or the Y axis (2)."""
    L = np.eye(3)
    L[0, 0] = L[axis, axis] = math.cosh(rapidity)
    L[0, axis] = L[axis, 0] = math.sinh(rapidity)
    return L


def length(L):
    """The distance a hyperbolic element of SO(2, 1) moves the points of its own geodesic, the
    length of the closed geodesic it makes: its trace is 1 + 2 cosh of it."""
    return math.acosh((np.trace(L) - 1) / 2)


class Tent(unittest.TestCase):
    def test_the_cuts_stand_where_the_captions_say(self):
        self.assertAlmostEqual(math.sinh(B), math.sqrt(2), places=12)
        # Adjacent cuts never meet: the product of their normals is cosh of their distance.
        for k in range(4):
            self.assertAlmostEqual(abs(cut(k) @ ETA @ cut((k + 1) % 4)), math.cosh(D), places=12)
        # Aminneborg, Bengtsson, Brill, Holst and Peldan's (9): an opening exists for alpha > 1/sqrt 2.
        self.assertGreater(ALPHA, 1 / math.sqrt(2))
        self.assertLess(ALPHA, 1)

    def test_a_fold_leaves_infinity_at_a_sixth_of_pi(self):
        """Their (10): tan t_P = sqrt(2 alpha^2 - 1). On the plane phi = pi/4 the fold is
        2 rho/(1 + rho^2) = sqrt(2) alpha cos t, which reaches rho = 1 at t_P and rho = 0 at pi/2."""
        self.assertAlmostEqual(math.atan(math.sqrt(2 * ALPHA ** 2 - 1)), T_P, places=12)
        self.assertAlmostEqual(math.sqrt(2) * ALPHA * math.cos(T_P), 1.0, places=12)
        self.assertAlmostEqual(math.sqrt(2) * ALPHA * math.cos(math.pi / 2), 0.0, places=12)
        # The horizon is the ray that reaches infinity at t_P, a time pi/2 after it leaves the axis.
        self.assertAlmostEqual(T_P - math.pi / 2, -math.pi / 3, places=12)

    def test_opposite_cuts_glued_make_one_horizon_of_four_segments(self):
        """The wormhole of their section V: a and b translate along X and Y by twice the cuts'
        distance from the axis, and the horizon is the closed geodesic of a b a^-1 b^-1."""
        a, b = boost(1, 2 * B), boost(2, 2 * B)
        for k, glue in ((2, a), (3, b)):
            image = glue @ cut(k)           # the normal of the cut's image, up to its sign
            np.testing.assert_allclose(np.abs(image @ ETA @ cut(k - 2)), 1.0, atol=1e-9)
        word = a @ b @ np.linalg.inv(a) @ np.linalg.inv(b)
        self.assertAlmostEqual(length(word), 4 * D, places=9)
        mass = (length(word) / (2 * math.pi)) ** 2
        self.assertAlmostEqual(mass, (2 * math.acosh(2) / math.pi) ** 2, places=9)
        self.assertAlmostEqual(mass, 0.703, places=3)

    def test_adjacent_cuts_glued_make_three_horizons(self):
        """Brill's doubling: the region between the diagonal and two cuts, and its mirror image.
        The three holes' horizons are the closed geodesics of the products of two reflections."""
        diagonal = np.array([0.0, math.cos(3 * math.pi / 4), math.sin(3 * math.pi / 4)])
        first = mirror(diagonal) @ mirror(cut(1))
        second = mirror(diagonal) @ mirror(cut(2))
        third = mirror(cut(1)) @ mirror(cut(2))
        self.assertAlmostEqual(length(first), D, places=9)
        self.assertAlmostEqual(length(second), D, places=9)
        self.assertAlmostEqual(length(third), 2 * D, places=9)


class Drawings(unittest.TestCase):
    def test_every_chart_has_a_spacetime_diagram_within_the_allowed_shapes(self):
        systems = load("diagrams")["systems"]
        self.assertEqual(sorted(systems), sorted(CHARTS))
        for name, views in systems.items():
            self.assertEqual(len(views), 1, name)
            X0, X1, Y0, Y1 = views[0]["box"]
            shape = (X1 - X0) * (2 if views[0]["mirror"] else 1) / (Y1 - Y0)
            self.assertTrue(0.5 <= shape <= 2.0, f"{name}: {shape}")

    def test_the_sausage_view_draws_the_folds_and_the_horizon_where_they_are(self):
        view, = load("diagrams")["systems"]["sausage"]
        X0, X1, Y0, Y1 = view["box"]

        def chart(point):
            return X0 + point[0] * (X1 - X0), Y0 + point[1] * (Y1 - Y0)
        folds = [m for m in view["markers"] if m["kind"] == "shell"]
        self.assertEqual(len(folds), 2)
        for fold in folds:
            for point in fold["lines"][0]:
                rho, t = chart(point)
                # Inside the drawing each point of a fold solves 2 rho/(1 + rho^2) = sqrt(2) alpha cos t.
                self.assertAlmostEqual(2 * rho / (1 + rho ** 2), math.sqrt(2) * ALPHA * math.cos(t), delta=4e-3)
            times = sorted(abs(chart(p)[1]) for p in fold["lines"][0])
            self.assertAlmostEqual(times[0], T_P, delta=4e-3)
            self.assertAlmostEqual(times[-1], math.pi / 2, delta=1e-3)
        horizon, = [m for m in view["markers"] if m["kind"] == "event"]
        line, = horizon["lines"]
        rho, t = chart(line[0])
        self.assertAlmostEqual(rho, 0.0, delta=1e-3)
        self.assertAlmostEqual(t, -math.pi / 3, delta=1e-3)
        for point in line:
            rho, t = chart(point)
            self.assertAlmostEqual(t, -math.pi / 3 + 2 * math.atan(rho), delta=4e-3)
        self.assertAlmostEqual(chart(line[-1])[1], T_P, delta=6e-3)

    def test_the_stereographic_view_draws_the_fold_as_a_straight_line_to_infinity(self):
        view, = load("diagrams")["systems"]["stereographic"]
        X0, X1, Y0, Y1 = view["box"]
        fold, = [m for m in view["markers"] if m["kind"] == "shell"]
        ends = [(X0 + p[0] * (X1 - X0), Y0 + p[1] * (Y1 - Y0)) for p in (fold["lines"][0][0], fold["lines"][0][-1])]
        far = max(ends)
        self.assertAlmostEqual(far[0], 4.0, delta=5e-3)
        self.assertAlmostEqual(far[1], -2 * math.sqrt(3), delta=5e-3)
        # That end is on infinity, x^2 - c^2 tau^2 = 4 l^2.
        self.assertAlmostEqual(far[0] ** 2 - far[1] ** 2, 4.0, delta=0.05)
        horizon, = [m for m in view["markers"] if m["kind"] == "event"]
        start = horizon["lines"][0][0]
        self.assertAlmostEqual(Y0 + start[1] * (Y1 - Y0), -(4 + 2 * math.sqrt(3)), delta=5e-3)

    def test_the_embedding_marks_four_cuts_and_four_segments_of_horizon_on_the_hyperboloid(self):
        view, = load("embedding")["views"]
        surface, = view["surfaces"]
        curves = surface["curves"]
        self.assertEqual([c["class"] for c in curves], ["cut"] * 4 + ["horizon"] * 4)
        for k, curve in enumerate(curves):
            P = np.array(curve["points"])
            V = P[:, 2] + 1
            # On the hyperboloid V^2 - X^2 - Y^2 = 1, as far as the file's decimals go.
            np.testing.assert_allclose(V ** 2 - P[:, 0] ** 2 - P[:, 1] ** 2, 1.0, atol=2e-2)
            if curve["class"] == "cut":
                np.testing.assert_allclose(np.column_stack([V, P[:, 0], P[:, 1]]) @ ETA @ cut(k), 0.0, atol=5e-3)
            else:
                step = np.diff(np.column_stack([V, P[:, 0], P[:, 1]]), axis=0)
                chords = np.sqrt(np.clip(np.einsum("ij,jk,ik->i", step, ETA, step), 0, None))
                self.assertAlmostEqual(float(chords.sum()), D, delta=2e-2)

    def test_the_conformal_diagram_puts_the_horizon_and_the_openings_where_the_caption_says(self):
        view, = load("conformal")["views"]
        horizon = [L["points"] for L in view["layers"] if L["class"] == "horizon"]
        self.assertEqual(len(horizon), 2)
        for line in horizon:
            (x0, t0), (x1, t1) = line[0], line[-1]
            self.assertAlmostEqual(x0, 0.0, places=3)
            self.assertAlmostEqual(t0, -math.pi / 3, places=3)
            self.assertAlmostEqual(abs(x1), math.pi / 2, places=3)
            self.assertAlmostEqual(t1, T_P, places=3)
            # A null line: it rises as far as it runs across.
            self.assertAlmostEqual(abs(x1 - x0), t1 - t0, places=3)
        for line in [L["points"] for L in view["layers"] if L["class"] == "boundary"]:
            self.assertAlmostEqual(abs(line[0][0]), math.pi / 2, places=3)
            self.assertEqual(sorted(round(p[1], 3) for p in (line[0], line[-1])), [round(-T_P, 3), round(T_P, 3)])
        for line in [L["points"] for L in view["layers"] if L["class"] == "singular"]:
            for x, t in line:
                self.assertAlmostEqual(math.sin(abs(x)), math.sqrt(2) * ALPHA * math.cos(t), delta=2e-3)


class Charts(unittest.TestCase):
    def setUp(self):
        self.metric = load("metrics")

    def test_every_chart_states_the_curvature_of_anti_de_sitter_space(self):
        self.assertEqual([c["id"] for c in self.metric["coordinates"]], CHARTS)
        for chart in self.metric["coordinates"]:
            self.assertEqual(chart["ricci_scalar"], "R = -\\dfrac{6}{\\ell^2}", chart["id"])
            self.assertEqual(chart["kretschmann"], "K = \\dfrac{12}{\\ell^4}", chart["id"])
            for variant in chart["weyl_tensor"]["variants"].values():
                self.assertEqual(variant["nonzero"], [], chart["id"])
            mixed = chart["ricci_tensor"]["variants"]["ul"]["nonzero"]
            self.assertEqual({entry["value"] for entry in mixed}, {"-\\dfrac{2}{\\ell^2}"}, chart["id"])
            self.assertEqual(len(mixed), 3, chart["id"])

    def test_the_spacetime_is_tagged_with_no_symmetry_in_time(self):
        """The exterior chart is static, and the whole spacetime has no Killing vector, so the
        chart that covers one outside alone does not speak for it."""
        tags = set(self.metric["tags"])
        self.assertFalse(tags & {"static", "stationary"})
        self.assertTrue({"Einstein space", "conformally flat", "three-dimensional"} <= tags)

    def test_each_relation_is_answered(self):
        for relation in self.metric["related"]:
            other = load("metrics", relation["id"])
            back = [r for r in other["related"] if r["id"] == ME]
            self.assertEqual(len(back), 1, relation["id"])
        kinds = {r["id"]: r["kind"] for r in self.metric["related"]}
        self.assertEqual(kinds["btz"], "piece")
        self.assertEqual(kinds["anti_de_sitter"], "locally_same")


if __name__ == "__main__":
    unittest.main()
