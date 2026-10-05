"""Gott and Li's self-creating universe, de Sitter space with each event identified with its images
under a boost: .venv.noindex/bin/python -m unittest discover -s _tools

What its texts and drawings state, held on the published files, with nothing but the files and
arithmetic: the horizons the drawings mark, the windows and periods they are drawn over, the
cylinders of the embedding diagram, the copies the conformal diagram tints, the numbers the History
and the captions quote, and the relations each side answers.
"""
import json
import math
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
ID = "self_creating_universe"
CHARTS = ["static", "kantowski_sachs", "steady_state", "conformal"]
BETA = 2 * math.pi          # the drawings' beta in units of r_0, Gott and Li's self-consistent value
MOMENTS = (0.25, 0.5, 1.0, 1.5)


def load(folder, metric_id=ID):
    return json.loads((DATA / folder / f"{metric_id}.json").read_text(encoding="utf-8"))


def chart_point(view, unit):
    """A point of a view's unit square in the drawn axes."""
    X0, X1, Y0, Y1 = view["box"]
    return X0 + unit[0] * (X1 - X0), Y0 + unit[1] * (Y1 - Y0)


class Charts(unittest.TestCase):
    def test_the_four_charts_are_gott_and_lis_and_each_states_its_identification(self):
        charts = {c["id"]: c for c in load("metrics")["coordinates"]}
        self.assertEqual(list(charts), CHARTS)
        self.assertIn("t \\in [0, \\beta/c)", charts["static"]["domains"])
        self.assertIn("l \\in [0, \\beta)", charts["kantowski_sachs"]["domains"])
        self.assertTrue(any("\\sim" in d and "e^{-\\beta/r_0}x" in d for d in charts["steady_state"]["domains"]))
        self.assertTrue(any("\\sim" in d and "e^{-\\beta/r_0}\\eta" in d for d in charts["conformal"]["domains"]))
        for chart in charts.values():
            self.assertEqual(chart["ricci_scalar"], "R = \\dfrac{12}{r_0^2}")
            self.assertEqual(chart["kretschmann"], "K = \\dfrac{24}{r_0^4}")
            self.assertEqual(chart["weyl_tensor"]["variants"]["llll"]["nonzero"], [])

    def test_the_boost_is_the_same_map_of_the_hyperboloid_in_every_chart(self):
        """(V, W, R) of Gott and Li's (63), (73) and (74) at r_0 = 1, before and after each chart's stated
        identification: V and W boosted by the rapidity beta, the radius of the sphere unchanged."""
        beta = 0.7
        ch, sh = math.cosh(beta), math.sinh(beta)

        def static(t, r):
            s = math.sqrt(1 - r * r)
            return s * math.sinh(t), s * math.cosh(t), r

        def kantowski_sachs(tau, l):
            return math.sinh(tau) * math.cosh(l), math.sinh(tau) * math.sinh(l), math.cosh(tau)

        def steady(tau, x):
            e = math.exp(tau)
            return math.sinh(tau) + e * x * x / 2, math.cosh(tau) - e * x * x / 2, e * x

        def conformal(eta, rho):
            return steady(-math.log(-eta), rho)

        shrink = math.exp(-beta)
        pairs = [(static, (0.3, 0.6), (0.3 + beta, 0.6)), (kantowski_sachs, (0.8, -0.4), (0.8, -0.4 + beta)),
                 (steady, (0.2, 1.7), (0.2 + beta, shrink * 1.7)), (conformal, (-0.9, 0.5), (-0.9 * shrink, 0.5 * shrink))]
        for chart, before, after in pairs:
            V, W, R = chart(*before)
            V2, W2, R2 = chart(*after)
            self.assertAlmostEqual(W * W + R * R - V * V, 1.0, places=12, msg=chart.__name__)
            self.assertAlmostEqual(V2, V * ch + W * sh, places=12, msg=chart.__name__)
            self.assertAlmostEqual(W2, W * ch + V * sh, places=12, msg=chart.__name__)
            self.assertAlmostEqual(R2, R, places=12, msg=chart.__name__)


class Drawings(unittest.TestCase):
    def setUp(self):
        self.systems = load("diagrams")["systems"]

    def view(self, system):
        view, = self.systems[system]
        return view

    def test_every_chart_is_drawn_over_one_period_at_the_self_consistent_beta(self):
        self.assertEqual(list(self.systems), CHARTS)
        for system in CHARTS:
            self.assertIn("$\\beta = 2 \\pi$", self.view(system)["settings"], system)
        self.assertAlmostEqual(self.view("static")["box"][3] - self.view("static")["box"][2], BETA)
        self.assertAlmostEqual(self.view("kantowski_sachs")["box"][1] - self.view("kantowski_sachs")["box"][0], BETA)
        self.assertAlmostEqual(self.view("steady_state")["box"][3] - self.view("steady_state")["box"][2], BETA)

    def test_the_static_plane_marks_the_cauchy_horizon_at_the_de_sitter_radius(self):
        view = self.view("static")
        marker, = [m for m in view["markers"] if m["kind"] == "grr"]
        for line in marker["lines"]:
            for point in line:
                self.assertAlmostEqual(chart_point(view, point)[0], 1.0, places=3)

    def test_the_steady_state_plane_marks_the_rays_x_equal_to_plus_and_minus_exp_minus_tau(self):
        view = self.view("steady_state")
        marked = [m for m in view["markers"] if m["kind"] == "event"]
        self.assertEqual(len(marked), 2)
        signs = set()
        for marker in marked:
            line, = marker["lines"]
            for point in line:
                x, tau = chart_point(view, point)
                self.assertAlmostEqual(abs(x), math.exp(-tau), delta=0.002)
                signs.add(x > 0)
        self.assertEqual(signs, {True, False})

    def test_the_conformal_plane_marks_the_ray_rho_equal_to_minus_eta(self):
        view = self.view("conformal")
        marker, = [m for m in view["markers"] if m["kind"] == "event"]
        line, = marker["lines"]
        for point in line:
            rho, eta = chart_point(view, point)
            self.assertAlmostEqual(rho, -eta, places=3)
        # Future infinity, where the conformal factor diverges, is no horizon and is not marked as one.
        self.assertFalse([m for m in view["markers"] if m["kind"] == "grr"])

    def test_the_moments_of_the_inflating_region_lie_where_the_charts_put_them(self):
        """c tau = tau_k is a line of the Kantowski-Sachs plane, x = cosh(tau_k) exp(-tau) on the steady
        state plane and rho = -eta cosh(tau_k) on the conformal one; no moment meets the static patch."""
        self.assertEqual(self.view("static").get("slices", []), [])
        for system in CHARTS[1:]:
            self.assertEqual(len(self.view(system)["slices"]), len(MOMENTS), system)
        for mark, tau_k in zip(self.view("kantowski_sachs")["slices"], MOMENTS):
            for line in mark["lines"]:
                for point in line:
                    self.assertAlmostEqual(chart_point(self.view("kantowski_sachs"), point)[1], tau_k, places=3)
        view = self.view("steady_state")
        for mark, tau_k in zip(view["slices"], MOMENTS):
            for line in mark["lines"]:
                for point in line:
                    x, tau = chart_point(view, point)
                    self.assertAlmostEqual(x, math.cosh(tau_k) * math.exp(-tau), delta=0.002)
        view = self.view("conformal")
        for mark, tau_k in zip(view["slices"], MOMENTS):
            for line in mark["lines"]:
                for point in line:
                    rho, eta = chart_point(view, point)
                    self.assertAlmostEqual(rho, -eta * math.cosh(tau_k), delta=0.002)


class Cylinders(unittest.TestCase):
    def test_each_moment_is_a_cylinder_of_radius_sinh_and_length_pi_cosh(self):
        view, = load("embedding")["views"]
        self.assertEqual([s["time"] for s in view["surfaces"]], list(MOMENTS))
        for surface in view["surfaces"]:
            tube, = surface["pieces"]
            tau = surface["time"]
            thetas, radii, heights = zip(*tube["points"])
            self.assertAlmostEqual(thetas[0], 0.0)
            self.assertAlmostEqual(thetas[-1], math.pi)
            for rho in radii:
                self.assertAlmostEqual(rho, math.sinh(tau), places=6)
            self.assertAlmostEqual(heights[-1] - heights[0], math.pi * math.cosh(tau), places=5)
            # With l running once round beta = 2 pi r_0 the circle's circumference is beta sinh(c tau/r_0).
            self.assertAlmostEqual(2 * math.pi * radii[0], BETA * math.sinh(tau), places=5)

    def test_the_movie_runs_from_the_first_moment_to_the_last(self):
        view, = load("embedding")["views"]
        values = [frame["value"] for frame in view["movie"]["frames"]]
        self.assertEqual((values[0], values[-1]), (MOMENTS[0], MOMENTS[-1]))
        self.assertEqual(values, sorted(values))


class Copies(unittest.TestCase):
    def test_the_conformal_diagram_tints_one_copy_in_each_chart(self):
        views = {v["id"]: v for v in load("conformal")["views"]}
        self.assertEqual(list(views), CHARTS)
        half = math.pi / 2
        for vid, view in views.items():
            cover, = [layer for layer in view["layers"] if layer["class"] == "cover"]
            xs, ts = zip(*cover["points"])
            if vid == "static":
                # Inside the triangle about the observer, between the diagonals of the square.
                self.assertTrue(all(abs(t) <= half - x + 1e-3 for x, t in zip(xs, ts)))
            elif vid == "kantowski_sachs":
                # Inside the triangle to the future of both horizons.
                self.assertTrue(all(t >= abs(x - half) - 1e-3 for x, t in zip(xs, ts)))
            else:
                # Above the observer's past horizon, the line from the bottom left to the top right.
                self.assertTrue(all(t >= x - half - 1e-3 for x, t in zip(xs, ts)))
        self.assertEqual(views["static"].get("slices", []), [])
        for vid in CHARTS[1:]:
            self.assertEqual(len(views[vid]["slices"]), len(MOMENTS), vid)


class Quoted(unittest.TestCase):
    def test_the_time_after_which_no_bubble_meets_its_image_is_0_086_at_beta_2_pi(self):
        """Gott and Li's tau_0 = r_0 ln((e^(b/2) + 1)/(e^(b/2) - 1)) at b = 2 pi, which the History and the
        Kantowski-Sachs caption quote: a ray from tau_0 covers half the circle, -ln tanh(tau_0/2) = b/2."""
        e = math.exp(BETA / 2)
        tau0 = math.log((e + 1) / (e - 1))
        self.assertEqual(f"{tau0:.3f}", "0.086")
        self.assertAlmostEqual(-math.log(math.tanh(tau0 / 2)), BETA / 2, places=12)
        self.assertIn("0.086\\,r_0/c", load("metrics")["history"])
        caption = " ".join(load("diagrams")["systems"]["kantowski_sachs"][0]["caption"])
        self.assertIn("0.086\\,r_0", caption)

    def test_the_region_of_closed_timelike_curves_holds_four_thirds_pi_beta_r_0_cubed(self):
        """The static patch's sqrt(-g) is r^2 sin(theta), so one period holds beta times the volume of a
        ball of radius r_0, summed here by the midpoint rule."""
        n = 2000
        ball = sum(4 * math.pi * ((k + 0.5) / n) ** 2 / n for k in range(n))
        self.assertAlmostEqual(BETA * ball, 4 * math.pi * BETA / 3, places=5)
        self.assertIn("\\tfrac{4}{3}\\pi\\beta r_0^3", load("metrics")["history"])


class Relations(unittest.TestCase):
    def test_each_related_spacetime_answers(self):
        related = {r["id"]: r["kind"] for r in load("metrics")["related"]}
        self.assertEqual(related, {"de_sitter": "locally_same", "misner": "conformal", "elliptic_de_sitter": "family",
                                   "gott_time_machine": "programme", "kantowski_sachs": "generalisation",
                                   "coleman_de_luccia": "family", "malament_hogarth": "generalisation",
                                   "universe_from_nothing": "family"})
        answers = {"locally_same": "locally_same", "conformal": "conformal", "family": "family",
                   "programme": "programme", "generalisation": "special_case"}
        for other, kind in related.items():
            back = [r["kind"] for r in load("metrics", other)["related"] if r["id"] == ID]
            self.assertEqual(back, [answers[kind]], other)


if __name__ == "__main__":
    unittest.main()
