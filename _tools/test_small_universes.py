"""Ellis's small universes: .venv.noindex/bin/python -m unittest discover -s _tools

What the texts and drawings of `small_universes` state, held on the published files. The torus
charts are held to dust at the scale factors their drawings declare and to the laps of light, the
ages and the redshifts their captions state; the hyperbolic chart and the horn to one Einstein
tensor, to the open dust their drawings declare, and the horn to Sokolov and Starobinskii's map
from the spherical form of hyperbolic space. The drawings are held to their captions: the past
light cone that reaches the bang three cells away, the galaxy drawn at six places, the laps of the
ray on the conformal diagram, the cylinder of the torus and the pseudosphere of the horn. The
tests that read the published components need sympy and are skipped where it is absent.
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

A = 11 / 50                                  # the open dust's a = A (cosh(eta) - 1), in a_0
NOW = 6 / 5 - A * math.log(11)               # the time at which a = a_0 and a' = 6/5


def published(system):
    metric = json.loads((DATA / "metrics" / "small_universes.json").read_text(encoding="utf-8"))
    return next(c for c in metric["coordinates"] if c["id"] == system)


def drawn(kind):
    return json.loads((DATA / kind / "small_universes.json").read_text(encoding="utf-8"))


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class Theory(unittest.TestCase):
    """The published components, read as the checker reads them, in the chart x^0 = ct with c = 1."""

    def read(self, system):
        import sympy as sp
        import verify_metrics as vm
        chart = published(system)
        reader = vm.Reader(chart["coords"], [p["symbol"] for p in chart["parameters"]], ())

        def tensor(field, variant=None):
            block = chart[field] if variant is None else chart[field]["variants"][variant]["nonzero"]
            out = sp.zeros(4, 4)
            for entry in block:
                i, j = (chart["coords"].index(x) for x in entry["indices"])
                out[i, j] = reader(entry["value"]).subs(reader.c, 1)
            return out
        return chart, reader, tensor

    def with_scale(self, expr, reader, time, scale):
        return expr.subs(reader.parameters["a"], scale).doit()

    def test_the_torus_drawn_is_dust(self):
        import sympy as sp
        for system, time, scale in (("torus", "t", lambda t: (3 * t / 2) ** sp.Rational(2, 3)),
                                    ("torus_conformal", "\\eta", lambda e: e ** 2 / 4)):
            chart, reader, tensor = self.read(system)
            t = reader.symbol[time]
            G = tensor("einstein_tensor", "ul").applyfunc(lambda e: self.with_scale(e, reader, t, scale(t)))
            for i in range(1, 4):
                self.assertEqual(sp.simplify(G[i, i]), 0, system)
            # The density is positive and falls as 1/a^3.
            density = sp.simplify(-G[0, 0] * scale(t) ** 3)
            self.assertFalse(density.has(t), system)
            self.assertGreater(density, 0, system)

    def test_no_component_of_the_torus_depends_on_a_periodic_coordinate(self):
        for system in ("torus", "torus_conformal"):
            chart, reader, tensor = self.read(system)
            g = tensor("metric_components")
            for name in ("x", "y", "z"):
                self.assertFalse(g.has(reader.symbol[name]), system)

    def test_the_hyperbolic_chart_and_the_horn_have_one_einstein_tensor(self):
        import sympy as sp
        mixed = []
        for system in ("hyperbolic", "horn"):
            chart, reader, tensor = self.read(system)
            G = tensor("einstein_tensor", "ul")
            a, t = reader.parameters["a"], reader.symbol["t"]
            rate = sp.Derivative(a, t)
            self.assertEqual(sp.simplify(G[0, 0] + 3 * (rate ** 2 - 1) / a ** 2), 0, system)
            for i in range(1, 4):
                self.assertEqual(sp.simplify(G[i, i] + (2 * a * sp.Derivative(a, (t, 2)) + rate ** 2 - 1) / a ** 2), 0, system)
            mixed.append(G)

    def test_the_open_dust_drawn_solves_the_published_equations(self):
        import sympy as sp
        eta = sp.Symbol("eta", positive=True)
        scale = sp.Rational(11, 50) * (sp.cosh(eta) - 1)
        rate = sp.sinh(eta) / (sp.cosh(eta) - 1)                 # da/d(ct) = (da/d eta)/a
        second = sp.diff(rate, eta) / scale
        for system in ("hyperbolic", "horn"):
            chart, reader, tensor = self.read(system)
            a, t = reader.parameters["a"], reader.symbol["t"]
            G = tensor("einstein_tensor", "ul")[1, 1]
            on = G.subs(sp.Derivative(a, (t, 2)), second).subs(sp.Derivative(a, t), rate).subs(a, scale)
            self.assertEqual(sp.simplify(on.rewrite(sp.exp)), 0, system)
        at = {eta: sp.log(11)}
        self.assertEqual(sp.simplify(scale.subs(at)), 1)
        self.assertEqual(sp.simplify(rate.subs(at)), sp.Rational(6, 5))
        self.assertAlmostEqual(float((sp.Rational(11, 50) * (sp.sinh(eta) - eta)).subs(at)), NOW, places=12)

    def test_the_horn_is_hyperbolic_space_by_sokolov_and_starobinskiis_map(self):
        import sympy as sp
        import print_charts as pc
        chart, reader, tensor = self.read("horn")
        a = reader.parameters["a"]
        chi, theta, phi = sp.symbols("chi theta phi", positive=True)
        image = pc.small_universes_horn_map(chi, theta, phi)
        horn = (tensor("metric_components")[1:, 1:] / a ** 2).subs(
            dict(zip((reader.symbol[n] for n in ("x", "y", "z")), image)), simultaneous=True)
        J = sp.Matrix(3, 3, lambda i, j: sp.diff(image[i], [chi, theta, phi][j]))
        _, there, sphere = self.read("hyperbolic")
        names = {there.symbol["\\chi"]: chi, there.symbol["\\theta"]: theta, there.symbol["\\phi"]: phi}
        wanted = (sphere("metric_components")[1:, 1:] / there.parameters["a"] ** 2).subs(names)
        for value in sp.flatten(J.T * horn * J - wanted):
            self.assertEqual(sp.simplify(value.rewrite(sp.exp)), 0)


class Numbers(unittest.TestCase):
    """The figures the captions state."""

    def test_the_laps_of_light_on_the_torus(self):
        # eta = 2 L sqrt(a) and a = (3ct/2L)^(2/3): once round at ct = L/12, twice at 2L/3, three times at 9L/4.
        for laps, ct in ((0.5, 1 / 96), (1, 1 / 12), (1.5, 9 / 32), (2, 2 / 3), (3, 9 / 4)):
            self.assertAlmostEqual(2 * math.sqrt((1.5 * ct) ** (2 / 3)), laps, places=12)
        self.assertAlmostEqual((1.5 * 9 / 4) ** (2 / 3), 9 / 4, places=12)

    def test_the_ages_and_redshifts_of_the_images(self):
        # An image at the comoving distance d, seen from eta = 3L, shows the galaxy at eta = 3L - d.
        age = lambda d: (3 - d) ** 3 / 12              # noqa: E731
        redshift = lambda d: 9 / (3 - d) ** 2          # noqa: E731
        self.assertEqual(round(age(0.3), 2), 1.64)
        self.assertEqual(round(age(2.7), 3), 0.002)
        self.assertEqual(round(redshift(0.3), 2), 1.23)
        self.assertAlmostEqual(redshift(2.7), 100, places=9)
        # Six images of the galaxy and two of ourselves on each side lie inside the horizon.
        self.assertEqual(sum(abs(n + 0.3) < 3 for n in range(-10, 10)), 6)
        self.assertEqual(sum(0 < n < 3 for n in range(10)), 2)

    def test_the_open_dust(self):
        self.assertAlmostEqual(A * (math.cosh(math.log(11)) - 1), 1.0, places=12)
        self.assertAlmostEqual(math.sqrt(1 + 2 * A), 6 / 5, places=12)
        self.assertAlmostEqual(1 - 1 / 1.2 ** 2, 11 / 36, places=12)
        self.assertEqual(round(math.log(11), 2), 2.40)
        self.assertEqual(round(NOW, 4), 0.6725)

    def test_where_an_observer_on_the_horn_has_seen_once_round(self):
        # Two images a period apart at x are 2a arsinh(b e^(-x)/2) apart, the chord of a horocycle.
        x = math.log(math.pi * math.sqrt(11) / 5)
        self.assertAlmostEqual(2 * math.asinh(2 * math.pi * math.exp(-x) / 2), math.log(11), places=12)
        self.assertEqual(round(x, 2), 0.73)

    def test_the_weeks_cell_lies_between_its_two_spheres(self):
        ball = lambda r: math.pi * (math.sinh(2 * r) - 2 * r)      # noqa: E731
        self.assertLess(ball(0.5192), 0.9427)
        self.assertLess(0.9427, ball(0.7525))
        self.assertLess(0.7525, math.log(11))


class Drawings(unittest.TestCase):
    """The drawings, held to what their captions say."""

    def view(self, system, view_id):
        return next(v for v in drawn("diagrams")["systems"][system] if v["id"] == view_id)

    def test_the_cells_are_one_period_wide_and_the_unrolled_views_six(self):
        for system in ("torus", "torus_conformal"):
            self.assertEqual(self.view(system, "cell")["box"][:2], [0, 1])
            self.assertEqual(self.view(system, "images")["box"][:2], [-3, 3])

    def test_the_unrolled_torus_draws_the_galaxy_six_times_and_the_faces_seven(self):
        for system in ("torus", "torus_conformal"):
            view = self.view(system, "images")
            kinds = [m["kind"] for m in view["markers"]]
            self.assertEqual(kinds.count("world"), 6, system)
            self.assertEqual(kinds.count("surface"), 7, system)
            places = sorted(m["lines"][0][0][0] for m in view["markers"] if m["kind"] == "world")
            for n, X in zip(range(-3, 3), places):
                self.assertAlmostEqual(X, (n + 0.3 + 3) / 6, places=3)

    def test_our_past_light_cone_reaches_the_bang_three_cells_away(self):
        for system, top in (("torus", 3.0), ("torus_conformal", 3.5)):
            view = self.view(system, "images")
            here = next(m for m in view["markers"] if m["kind"] == "mark")["points"][0]
            self.assertAlmostEqual(here[0], 0.5, places=3)
            self.assertAlmostEqual(here[1], (9 / 4 if system == "torus" else 3) / top, places=3)
            cone = next(m for m in view["markers"] if m["kind"] == "past")["lines"]
            ends = sorted(line[-1][0] for line in cone)
            self.assertEqual(len(cone), 2)
            for line in cone:
                self.assertLess(line[-1][1], 0.002, system)
            self.assertLess(ends[0], 0.005, system)
            self.assertGreater(ends[1], 0.995, system)

    def test_the_bang_is_marked_singular_on_every_view(self):
        for system, views in drawn("diagrams")["systems"].items():
            for view in views:
                singular = [m for m in view["markers"] if m["kind"] == "singular"]
                self.assertEqual([m["edges"] for m in singular], [["bottom"]], f"{system}/{view['id']}")

    def test_the_hyperbolic_view_marks_the_weeks_cell_and_now(self):
        view = self.view("hyperbolic", "radial")
        lines = {m["kind"]: m["lines"][0][0][0] for m in view["markers"] if m["kind"] in ("surface", "shell")}
        self.assertAlmostEqual(lines["surface"], 0.5192 / 1.5, places=3)
        self.assertAlmostEqual(lines["shell"], 0.7525 / 1.5, places=3)
        here = next(m for m in view["markers"] if m["kind"] == "mark")["points"][0]
        self.assertAlmostEqual(here[1], NOW / 1.5, places=3)
        reference = next(m for m in view["markers"] if m["kind"] == "reference")
        self.assertAlmostEqual(reference["y"], NOW / 1.5, places=3)

    def test_the_horn_marks_where_light_has_been_once_round(self):
        view = self.view("horn", "along")
        line = next(m for m in view["markers"] if m["kind"] == "surface")["lines"][0]
        self.assertAlmostEqual(line[0][0], math.log(math.pi * math.sqrt(11) / 5) / 3, places=3)

    def test_the_ray_on_the_conformal_diagram_laps_the_universe(self):
        for view in drawn("conformal")["views"][:2]:
            laps = [layer["points"] for layer in view["layers"] if layer["class"] == "null"]
            self.assertEqual(len(laps), 3)
            for k, ((X0, T0), (X1, T1)) in enumerate(laps):
                # Each lap starts on x = 0 at eta = k L, where T = 2 arctan(k), and runs at 45 degrees
                # to x = L at eta = (k + 1) L.
                self.assertAlmostEqual(X0, 0.0, places=3)
                self.assertAlmostEqual(T0, 2 * math.atan(k), places=3)
                self.assertAlmostEqual(X1 - X0, T1 - T0, places=3)
                self.assertAlmostEqual(X1, math.atan(k + 2) - math.atan(k), places=3)

    def test_the_torus_is_embedded_as_a_cylinder_whose_circle_is_its_period(self):
        view = next(v for v in drawn("embedding")["views"] if v["id"] == "torus")
        for surface, laps in zip(view["surfaces"], (0.5, 1, 1.5)):
            a = (1.5 * surface["time"]) ** (2 / 3)
            self.assertAlmostEqual(2 * math.sqrt(a), laps, places=4)
            points = surface["pieces"][0]["points"]
            for y, rho, z in points:
                # The file rounds the moment's time to six places, so the scale factor is known to five.
                self.assertAlmostEqual(rho, a / (2 * math.pi), places=5)
                self.assertAlmostEqual(z, a * (y - 0.5), places=5)
            self.assertEqual((points[0][0], points[-1][0]), (0.0, 1.0))

    def test_the_horn_is_embedded_as_the_pseudosphere(self):
        view = next(v for v in drawn("embedding")["views"] if v["id"] == "horn")
        points = view["surfaces"][0]["pieces"][0]["points"]
        self.assertEqual((points[0][0], points[-1][0]), (0.0, 3.0))
        for x, rho, z in points:
            s = math.sqrt(max(0.0, 1 - math.exp(-2 * x)))
            self.assertAlmostEqual(rho, math.exp(-x), places=6)
            self.assertAlmostEqual(z, math.atanh(s) - s, places=5)


if __name__ == "__main__":
    unittest.main()
