"""The barrier-top action, thermal approximation, and Euclidean/continuation distinction.

Tensor identities and ambient pullbacks belong to the chart generator and metric checker.
These checks exercise the physical interpretation and the generated drawing boundaries.
"""
import json
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
sys.path.insert(0, str(ROOT / "_tools" / "derivations"))
import hawking_moss as hm
import verify_metrics as vm


def load(folder, mid=hm.METRIC):
    return json.loads((DATA / folder / f"{mid}.json").read_text())


def integral(function, lo, hi, steps=2000):
    h = (hi-lo)/steps
    return h/3*(function(lo)+function(hi)+sum((4 if i % 2 else 2)*function(lo+i*h)
                                             for i in range(1, steps)))


def value(chart, field, at, indices=None):
    reader = vm.Reader(chart["coords"], [p["symbol"] for p in chart["parameters"]], ())
    text = (next(c["value"] for c in chart[field] if c["indices"] == list(indices)) if indices else
            chart[field].partition("=")[2])
    expression = reader(text)
    subs = {reader.symbol[c]: at[reader._plain(c)] for c in chart["coords"] if reader._plain(c) in at}
    subs.update({reader.parameters[k]: v for k, v in at.items() if k in reader.parameters})
    return float(expression.subs(subs))


class Saddle(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.metric = load("metrics")
        cls.charts = {c["id"]: c for c in cls.metric["coordinates"]}

    def test_sphere_and_continuations_have_explicit_distinct_signatures(self):
        self.assertEqual(self.metric["signature"], "(+,+,+,+)")
        self.assertEqual(tuple(self.charts), hm.SYSTEMS)
        for system in hm.LORENTZIAN:
            self.assertIn("continuation", self.charts[system]["name"].lower())
        self.assertIn("(+,-,+,+)", self.charts["wall_continuation"]["convention"])

    def test_the_thermal_period_removes_the_polar_conical_defect(self):
        chart = self.charts["thermal"]
        self.assertIn(r"\tau \in [0, 2\pi\ell)", chart["domains"])
        self.assertIn(r"\tau \in [0, 2\pi\ell)", self.charts["thermal_angle"]["domains"])
        for ell in (0.8, 1.7, 3.2):
            for eps in (1e-3, 1e-5, 1e-7):
                r = ell*(1-eps)
                at = {"ell": ell, "r": r, "theta": math.pi/2}
                circumference = 2*math.pi*ell*math.sqrt(value(chart, "metric_components", at, (r"\tau", r"\tau")))
                distance = ell*math.acos(r/ell)
                self.assertAlmostEqual(circumference/distance/(2*math.pi), 1, delta=eps)

    def test_full_sphere_action_from_its_volume_and_stored_ricci_scalar(self):
        chart = self.charts["four_sphere"]
        for G, V in ((0.7, 0.2), (1.3, 0.04), (2.1, 0.5)):
            ell = math.sqrt(3/(8*math.pi*G*V))
            at = {"ell": ell, "xi": ell*math.pi/2, "chi": math.pi/2, "theta": math.pi/2}
            # The equatorial angular radius follows from the stored angular component.
            radius = math.sqrt(value(chart, "metric_components", at, (r"\phi", r"\phi")))
            volume = 2*math.pi**2*radius**3*integral(lambda xi: math.sin(xi/ell)**3, 0, math.pi*ell)
            R = value(chart, "ricci_scalar", at)
            action = (-R/(16*math.pi*G)+V)*volume
            self.assertAlmostEqual(action/(-3/(8*G*G*V)), 1, places=10)
            self.assertAlmostEqual(volume/(8*math.pi**2*ell**4/3), 1, places=10)

    def test_false_vacuum_subtraction_and_thermal_expansion(self):
        history = self.metric["history"]
        self.assertIn(r"B = 3(1/V_F - 1/V_{\mathrm{top}})/(8G^2)", history)
        self.assertIn("potential differences small enough", history)
        G, VF = 0.6, 0.2
        ellF = math.sqrt(3/(8*math.pi*G*VF))
        for fraction in (0.01, 0.001, 0.0001):
            delta = VF*fraction
            top = VF+delta
            B = -3/(8*G*G*top)+3/(8*G*G*VF)
            # Static lapse times the proper spatial volume density cancels exactly.
            energy = integral(lambda r: 4*math.pi*r*r*delta, 0, ellF)
            thermal = energy*2*math.pi*ellF
            self.assertGreater(B, 0)
            self.assertAlmostEqual(B/thermal, VF/top, places=10)
            proper_volume = 4*math.pi*ellF**3*integral(lambda eta: math.sin(eta)**2, 0, math.pi/2)
            self.assertGreater(abs(proper_volume*delta/energy-1), 0.5)


class DrawingDomains(unittest.TestCase):
    def test_only_constant_field_lorentzian_sections_have_light_cones(self):
        systems = load("diagrams")["systems"]
        self.assertEqual(tuple(systems), hm.LORENTZIAN)
        for views in systems.values():
            for view in views:
                self.assertIn("constant field continuation", " ".join(view["caption"]).lower())
        import hawking_moss_drawings as drawings
        self.assertEqual(tuple(drawings.NOT_DRAWN), hm.EUCLIDEAN)
        self.assertTrue(all("Positive definite" in text for text in drawings.NOT_DRAWN.values()))

    def test_embedding_is_a_whole_euclidean_sphere(self):
        view, = load("embedding")["views"]
        piece, = view["surfaces"][0]["pieces"]
        points = piece["points"]
        self.assertAlmostEqual(points[0][0], 0)
        self.assertAlmostEqual(points[-1][0], math.pi, places=5)
        self.assertAlmostEqual(points[0][2], 1, places=5)
        self.assertAlmostEqual(points[-1][2], -1, places=5)
        self.assertTrue(any(z > 0.9 for _, _, z in points))
        self.assertTrue(any(z < -0.9 for _, _, z in points))
        self.assertNotIn("space", piece)

    def test_conformal_view_is_only_the_named_static_patch(self):
        view, = load("conformal")["views"]
        self.assertEqual(view["system"], "static_continuation")
        self.assertIn("constant field continuation", " ".join(view["caption"]))
        cover, = [layer for layer in view["layers"] if layer["class"] == "cover"]
        for X, T in cover["points"]:
            self.assertGreaterEqual(X, -1e-6)
            self.assertLessEqual(X+abs(T), math.pi/2+1e-5)
        self.assertNotIn("singular", {layer["class"] for layer in view["layers"]})


if __name__ == "__main__":
    unittest.main()
