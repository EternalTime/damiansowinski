"""The universe from nothing, half a 4-sphere joined at its equator to the waist of de Sitter's
closed universe: .venv.noindex/bin/python -m unittest discover -s _tools

What its texts and drawings state, held on the published files, with nothing but the files and
arithmetic: the charts and the signature each states, the join at rest on a 3-sphere of radius l,
the windows the drawings cover, the bowl and the skirt of the embedding diagram, the half of de
Sitter's square the conformal diagram draws, the action the History quotes, and the relations each
side answers.
"""
import json
import math
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
ID = "universe_from_nothing"
CHARTS = ["closed", "four_sphere", "scale_factor", "conformal", "lapse"]
LORENTZIAN = ["closed", "scale_factor", "conformal", "lapse"]


def load(folder, metric_id=ID):
    return json.loads((DATA / folder / f"{metric_id}.json").read_text(encoding="utf-8"))


def chart_point(view, unit):
    """A point of a view's unit square in the drawn axes."""
    X0, X1, Y0, Y1 = view["box"]
    return X0 + unit[0] * (X1 - X0), Y0 + unit[1] * (Y1 - Y0)


# Each chart at l = 1 as (X_0, a) on the sphere or the hyperboloid, the time coordinate in units of l,
# and the sign the time takes in the five-dimensional metric: +1 for Euclid's, -1 for Minkowski's.
EMBEDDINGS = {
    "closed": (lambda t: (math.sinh(t), math.cosh(t)), -1),
    "four_sphere": (lambda tau: (math.sin(tau), math.cos(tau)), 1),
    "conformal": (lambda eta: (math.tan(eta), 1 / math.cos(eta)), -1),
    "lapse": (lambda t: (t, math.sqrt(1 + t * t)), -1),
    "scale_factor_below": (lambda a: (-math.sqrt(1 - a * a), a), 1),
    "scale_factor_above": (lambda a: (math.sqrt(a * a - 1), a), -1),
}
# g_00 of each chart at l = 1, from its published line element.
G00 = {
    "closed": lambda t: -1.0,
    "four_sphere": lambda tau: 1.0,
    "conformal": lambda eta: -1 / math.cos(eta) ** 2,
    "lapse": lambda t: -1 / (1 + t * t),
    "scale_factor_below": lambda a: 1 / (1 - a * a),
    "scale_factor_above": lambda a: 1 / (1 - a * a),
}


class Charts(unittest.TestCase):
    def test_the_five_charts_each_of_constant_curvature(self):
        charts = {c["id"]: c for c in load("metrics")["coordinates"]}
        self.assertEqual(list(charts), CHARTS)
        for chart in charts.values():
            self.assertEqual(chart["ricci_scalar"], "R = \\dfrac{12}{\\ell^2}")
            self.assertEqual(chart["kretschmann"], "K = \\dfrac{24}{\\ell^4}")
            self.assertEqual(chart["weyl_tensor"]["variants"]["llll"]["nonzero"], [])

    def test_each_chart_is_the_sphere_or_the_hyperboloid_of_radius_l(self):
        """The time direction pulled back from flat space of five dimensions, sign (dX_0)^2 + (da)^2, by
        central differences, against g_00; the 3-sphere's radius is a in every chart."""
        h = 1e-6
        for name, (embed, sign) in EMBEDDINGS.items():
            samples = {"scale_factor_below": (0.2, 0.5, 0.9), "scale_factor_above": (1.1, 2.0, 4.0),
                       "four_sphere": (-1.4, -0.7, -0.1), "conformal": (0.1, 0.7, 1.3)}.get(name, (0.1, 1.0, 2.5))
            for x in samples:
                X0, a = embed(x)
                self.assertAlmostEqual(sign * X0 * X0 + a * a, 1.0, places=12, msg=name)
                dX0 = (embed(x + h)[0] - embed(x - h)[0]) / (2 * h)
                da = (embed(x + h)[1] - embed(x - h)[1]) / (2 * h)
                g00 = sign * dX0 * dX0 + da * da
                self.assertAlmostEqual(g00 / G00[name](x), 1.0, places=6, msg=f"{name} at {x}")

    def test_the_four_sphere_is_riemannian_and_the_scale_factor_changes_signature_at_l(self):
        self.assertGreater(G00["four_sphere"](-0.5), 0)
        self.assertGreater(G00["scale_factor_below"](0.5), 0)
        self.assertLess(G00["scale_factor_above"](2.0), 0)
        self.assertIn("(+,+,+,+)", load("metrics")["convention"])

    def test_the_halves_meet_at_rest_on_a_three_sphere_of_radius_l(self):
        """At t = 0 and tau = 0 the radius is l and its rate of change zero from both sides, so the extrinsic
        curvature of the join vanishes and the junction needs no shell."""
        h = 1e-6
        for name in ("closed", "four_sphere", "conformal", "lapse"):
            embed = EMBEDDINGS[name][0]
            self.assertAlmostEqual(embed(0.0)[1], 1.0, places=12, msg=name)
            self.assertAlmostEqual((embed(h)[1] - embed(-h)[1]) / (2 * h), 0.0, places=6, msg=name)
        # The south pole of the 4-sphere, tau = -pi l/2, and of the scale factor chart, a = 0.
        self.assertAlmostEqual(EMBEDDINGS["four_sphere"][0](-math.pi / 2)[1], 0.0, places=12)
        self.assertAlmostEqual(EMBEDDINGS["scale_factor_below"][0](0.0)[0], EMBEDDINGS["four_sphere"][0](-math.pi / 2)[0])


class Drawings(unittest.TestCase):
    def setUp(self):
        self.systems = load("diagrams")["systems"]

    def view(self, system):
        view, = self.systems[system]
        return view

    def test_every_lorentzian_chart_is_drawn_from_the_waist_up(self):
        self.assertEqual(list(self.systems), LORENTZIAN)
        for system in LORENTZIAN:
            X0, X1, Y0, Y1 = self.view(system)["box"]
            self.assertEqual((X0, X1), (0, math.pi), system)
            self.assertEqual(Y0, 1 if system == "scale_factor" else 0, system)
        self.assertAlmostEqual(self.view("conformal")["box"][3], math.pi / 2)

    def test_a_ray_from_the_waist_covers_at_most_a_quarter_turn(self):
        """chi changes by arctan sinh(ct/l) after the waist, which tends to pi/2, as the captions state, and is
        the arcsec of the radius a = l cosh(ct/l), as the scale factor chart's rays keep."""
        for t in (0.5, 3.0, 20.0):
            self.assertLess(math.atan(math.sinh(t)), math.pi / 2)
            self.assertAlmostEqual(math.acos(1 / math.cosh(t)), math.atan(math.sinh(t)), places=12)
        self.assertAlmostEqual(math.atan(math.sinh(40.0)), math.pi / 2, places=12)
        caption = " ".join(self.view("closed")["caption"])
        self.assertIn("\\arctan\\sinh(ct/\\ell)", caption)


class Bowl(unittest.TestCase):
    def setUp(self):
        view, = load("embedding")["views"]
        surface, = view["surfaces"]
        self.bowl, self.skirt = surface["pieces"]

    def test_the_bowl_is_a_hemisphere_in_flat_space_and_the_skirt_a_hyperboloid_in_minkowski_space(self):
        self.assertNotIn("space", self.bowl)
        self.assertEqual(self.skirt["space"], "minkowski")
        for a, rho, z in self.bowl["points"]:
            self.assertAlmostEqual(rho, a, places=6)
            self.assertAlmostEqual(rho * rho + z * z, 1.0, places=5)
            self.assertLessEqual(a, 1.0)
        for a, rho, z in self.skirt["points"]:
            self.assertAlmostEqual(rho, a, places=6)
            self.assertAlmostEqual(rho * rho - z * z, 1.0, places=6)
            self.assertGreaterEqual(a, 1.0)

    def test_the_two_meet_on_the_circle_a_equal_l_and_the_skirt_is_timelike(self):
        self.assertEqual(self.bowl["points"][0][:3], [0.0, 0.0, -1.0])
        self.assertEqual(self.bowl["points"][-1][1:], self.skirt["points"][0][1:])
        self.assertEqual(self.bowl["points"][-1][1:], [1.0, 0.0])
        for (_, r0, z0), (_, r1, z1) in zip(self.skirt["points"], self.skirt["points"][1:]):
            self.assertGreater((z1 - z0) ** 2, (r1 - r0) ** 2)

    def test_the_proper_time_up_the_skirt_is_the_closed_slicings(self):
        """From a = l to 2.5 l the skirt's chords add up to l arcosh(2.5), the proper time of the closed
        slicing from its waist to a = 2.5 l."""
        pts = self.skirt["points"]
        total = sum(math.sqrt((z1 - z0) ** 2 - (r1 - r0) ** 2) for (_, r0, z0), (_, r1, z1) in zip(pts, pts[1:]))
        self.assertAlmostEqual(total, math.acosh(2.5), places=3)


class Square(unittest.TestCase):
    def test_the_conformal_diagram_is_the_upper_half_of_de_sitters_square(self):
        views = {v["id"]: v for v in load("conformal")["views"]}
        self.assertEqual(list(views), LORENTZIAN)
        for vid, view in views.items():
            cover, = [layer for layer in view["layers"] if layer["class"] == "cover"]
            for x, t in cover["points"]:
                self.assertTrue(-1e-4 <= t <= math.pi / 2 + 1e-4 and -1e-4 <= x <= math.pi + 1e-4, vid)


class Quoted(unittest.TestCase):
    def test_the_four_spheres_action_is_minus_three_over_eight_g_squared_rho(self):
        """With hbar = c = 1 the Euclidean action of the whole 4-sphere is -rho_v times its volume
        8 pi^2 a_0^4/3, a_0^2 = 3/(8 pi G rho_v), which is -3/(8 G^2 rho_v) for any G and rho_v."""
        for G, rho in ((1.0, 1.0), (0.3, 2.7), (5.0, 0.01)):
            a0 = math.sqrt(3 / (8 * math.pi * G * rho))
            action = -rho * 8 * math.pi ** 2 * a0 ** 4 / 3
            self.assertAlmostEqual(action / (-3 / (8 * G * G * rho)), 1.0, places=12)
        self.assertIn("S_E = -3/(8G^2\\rho_v)", load("metrics")["history"])


class Relations(unittest.TestCase):
    def test_each_related_spacetime_answers(self):
        related = {r["id"]: r["kind"] for r in load("metrics")["related"]}
        self.assertEqual(related, {"de_sitter": "piece", "gravitational_instantons": "family",
                                   "coleman_de_luccia": "family", "self_creating_universe": "family",
                                   "nariai": "family", "kantowski_sachs": "family", "frw": "family",
                                   "einstein_static": "conformal"})
        answers = {"piece": "composite", "family": "family", "conformal": "conformal"}
        for other, kind in related.items():
            back = [r["kind"] for r in load("metrics", other)["related"] if r["id"] == ID]
            self.assertEqual(back, [answers[kind]], other)


if __name__ == "__main__":
    unittest.main()
