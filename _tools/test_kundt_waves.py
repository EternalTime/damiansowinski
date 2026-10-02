"""Kundt's waves, plane-fronted gravitational waves whose rays turn from front to front:
python3 -m unittest discover -s _tools

What its texts and drawings state, held on the published files. The first class needs nothing
but the files: the light cones of every plane drawn, one edge along the rays and the other
leaning by the Riccati slope its caption gives, the fronts of the figure touching the light cone
with parallel null rays on each, and the wave fronts on their sphere and their hyperboloid. The
second reads the published tensors through the checker's Reader and holds them to the field
equation, to flat space without the wave, to the scalars that vanish while the curvature does
not, and to the declared derivatives the Kerr-Schild chart is checked by; it needs sympy and is
skipped where it is absent.
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

PULSE = lambda u: math.exp(-4 * u * u)  # noqa: E731

# dv/du of the null curves that cross the rays on each plane of u and the coordinate along them,
# at ell = 1: v^2 + G/x in Podolsky and Belan's chart, w^2/2x^2 + 2xG in Kundt's canonical form
# with G = (x^2 - y^2) e^{-4u^2}, and (kappa v^2 + p h/q)/2 with a cosmological constant, where
# h = e^{-4u^2}(xi^2 - eta^2)(2 + p)/3p at xi = 1, eta = 0.
SLOPES = {
    ("kundt", "front"): lambda u, w: w * w / 2 + 2 * (1 - 1.21) * PULSE(u),
    ("podolsky_belan", "near"): lambda u, v: v * v + PULSE(u),
    ("podolsky_belan", "far"): lambda u, v: v * v + 2 * PULSE(u),
    ("simplest_wave", "front"): lambda u, v: v * v + 1,
    ("ozsvath_robinson_rozga", "de_sitter"): lambda u, v: (v * v + (5 / 3) * PULSE(u) * 13 / 15) / 2,
    ("ozsvath_robinson_rozga", "anti_de_sitter"): lambda u, v: (v * v + (3 / 4) * PULSE(u) * 11 / 9) / 2,
}


def load(folder):
    return json.loads((DATA / folder / "kundt_waves.json").read_text(encoding="utf-8"))


class Drawings(unittest.TestCase):
    def test_every_view_is_drawn(self):
        diagrams = load("diagrams")
        drawn = {(system, view["id"]) for system, views in diagrams["systems"].items() for view in views}
        self.assertEqual(drawn, set(SLOPES) | {("simplest_wave", "depth")})
        self.assertEqual({(system, view["id"]) for system, views in diagrams["projections"].items() for view in views},
                         {("kerr_schild", "fronts")})

    def test_one_edge_of_every_cone_runs_along_the_rays(self):
        """On a plane of u and the coordinate along the rays the rays keep u, and the other null
        direction is the Riccati slope, drawn with (v - u)/2 across and (u + v)/2 up."""
        for system, views in load("diagrams")["systems"].items():
            for view in views:
                if (system, view["id"]) not in SLOPES:
                    continue
                self.assertEqual(view["to_display"], [[-0.5, 0.5], [0.5, 0.5]])
                self.assertEqual(view["box"], [-2, 2, -2, 2])
                slope = SLOPES[(system, view["id"])]
                for cone in view["cones"]:
                    X, Y = (-2 + 4 * c for c in cone["at"])
                    u, v = Y - X, X + Y
                    k = slope(u, v)
                    across, up = (k - 1) / 2, (1 + k) / 2
                    norm = math.hypot(across, up)
                    self.assertAlmostEqual(cone["b"][0], math.sqrt(0.5), places=4)
                    self.assertAlmostEqual(cone["b"][1], math.sqrt(0.5), places=4)
                    self.assertAlmostEqual(cone["a"][0], across / norm, places=3, msg=f"{system}/{view['id']}")
                    self.assertAlmostEqual(cone["a"][1], up / norm, places=3, msg=f"{system}/{view['id']}")

    def test_the_cones_of_the_simplest_wave_narrow_toward_the_envelope(self):
        """With d/du divided out the plane of v and x is dx^2 - x^2 dv^2/(v^2 + x), so a cone's
        edges are dx = +-x dv/sqrt(v^2 + x), and it closes at x = 0."""
        view, = [v for v in load("diagrams")["systems"]["simplest_wave"] if v["id"] == "depth"]
        self.assertEqual(view["box"], [0, 3, -2, 2])
        widths = {}
        for cone in view["cones"]:
            x, v = 3 * cone["at"][0], -2 + 4 * cone["at"][1]
            across, up = x / math.sqrt(v * v + x) / 3, 1 / 4
            norm = math.hypot(across, up)
            for edge, sign in ((cone["a"], -1), (cone["b"], 1)):
                self.assertAlmostEqual(edge[0], sign * across / norm, places=3)
                self.assertAlmostEqual(edge[1], up / norm, places=3)
            widths.setdefault(cone["at"][1], []).append((x, across))
        for row in widths.values():
            ordered = [w for _, w in sorted(row)]
            self.assertEqual(ordered, sorted(ordered))
        self.assertEqual(len(view["hatch"]), 1)

    def test_the_embedded_fronts_are_marked_on_their_own_planes_alone(self):
        marked = {(system, view["id"]): [s["view"] for s in view.get("slices", [])]
                  for system, views in load("diagrams")["systems"].items() for view in views}
        self.assertEqual({k: v for k, v in marked.items() if v},
                         {("ozsvath_robinson_rozga", "de_sitter"): ["sphere"],
                          ("ozsvath_robinson_rozga", "anti_de_sitter"): ["hyperbolic"]})


class Fronts(unittest.TestCase):
    """The figure of the Kerr-Schild chart, from its pieces in the drawing's (X, Z, T)."""

    def setUp(self):
        figure, = load("diagrams")["projections"]["kerr_schild"]
        self.lines = {}
        for line in figure["turn"]["lines"]:
            self.lines.setdefault(line["class"], []).append(line["points"])
        self.figure = figure

    def test_it_turns_and_names_every_class_it_draws(self):
        self.assertEqual(self.figure["turn"]["cones"], [])
        self.assertEqual({entry[1] for entry in self.figure["legend"]}, {"critical", "world", "above"})

    def test_the_envelope_is_the_light_cone(self):
        """Its circles and the lines where the fronts touch it: X^2 + Z^2 = c^2T^2."""
        self.assertEqual(len(self.lines["critical"]), 2 + 6)
        for line in self.lines["critical"]:
            for X, Z, T in line:
                self.assertAlmostEqual(X * X + Z * Z, T * T, places=5)

    def test_every_front_is_a_half_line_touching_the_circle_of_its_time(self):
        """It starts on the circle of radius cT, runs at right angles to the radius there, and
        the distance along it is x = sqrt(X^2 + Z^2 - c^2T^2)."""
        self.assertEqual(len(self.lines["world"]), 12)
        for (X0, Z0, T0), (X1, Z1, T1) in self.lines["world"]:
            self.assertEqual(T0, T1)
            self.assertIn(T0, (1.0, 2.0))
            self.assertAlmostEqual(math.hypot(X0, Z0), T0, places=5)
            self.assertAlmostEqual((X1 - X0) * X0 + (Z1 - Z0) * Z0, 0.0, places=5)
            self.assertAlmostEqual(math.sqrt(X1 * X1 + Z1 * Z1 - T1 * T1), math.hypot(X1 - X0, Z1 - Z0), places=5)

    def test_the_rays_are_null_and_parallel_on_each_front(self):
        """Two rays to a front, each along (cT, X, Z) = (1, sin, cos) of its front's angle, the
        direction of the null line where the front touches the cone."""
        edges = [line for line in self.lines["critical"] if len(line) == 2]
        directions = set()
        for (X0, Z0, T0), (X1, Z1, T1) in edges:
            directions.add((round((X1 - X0) / (T1 - T0), 5), round((Z1 - Z0) / (T1 - T0), 5)))
        self.assertEqual(len(directions), 6)
        self.assertEqual(len(self.lines["above"]), 12)
        seen = {}
        for (X0, Z0, T0), (X1, Z1, T1) in self.lines["above"]:
            dX, dZ, dT = X1 - X0, Z1 - Z0, T1 - T0
            self.assertAlmostEqual(dX * dX + dZ * dZ, dT * dT, places=5)
            direction = (round(dX / dT, 5), round(dZ / dT, 5))
            self.assertIn(direction, directions)
            seen[direction] = seen.get(direction, 0) + 1
            # x is the same at both ends: a ray keeps its distance from the envelope.
            self.assertAlmostEqual(X0 * X0 + Z0 * Z0 - T0 * T0, X1 * X1 + Z1 * Z1 - T1 * T1, places=5)
        self.assertEqual(set(seen.values()), {2})


class WaveFronts(unittest.TestCase):
    def profile(self, view_id, piece_id):
        view, = [v for v in load("embedding")["views"] if v["id"] == view_id]
        piece, = [p for p in view["surfaces"][0]["pieces"] if p["id"] == piece_id]
        return [[float(c) for c in point] for point in piece["points"]]

    def test_in_de_sitter_space_the_front_is_a_hemisphere(self):
        """The circle r of the chart stands at rho = sin(theta), z = -cos(theta) with
        theta = 2 arctan(r/2), and the front ends on the equator, the envelope q = 0 at r = 2."""
        points = self.profile("sphere", "front")
        for r, rho, z in points:
            theta = 2 * math.atan(r / 2)
            self.assertAlmostEqual(rho, math.sin(theta), places=4)
            self.assertAlmostEqual(z, -math.cos(theta), places=4)
        self.assertEqual(points[0][:2], [0.0, 0.0])
        self.assertAlmostEqual(points[-1][0], 2.0, places=6)
        self.assertAlmostEqual(points[-1][2], 0.0, places=4)

    def test_in_anti_de_sitter_space_the_surface_is_a_sheet_of_a_hyperboloid(self):
        for r, rho, z in self.profile("hyperbolic", "sheet"):
            self.assertAlmostEqual((z + 1) ** 2 - rho ** 2, 1.0, places=3)
            self.assertAlmostEqual(rho, r / (1 - r * r / 4), places=3)

    def test_the_envelope_and_the_lines_of_constant_q_over_p_lie_on_one_side(self):
        """The envelope is the section X = 0 of the sheet, and the lines q/p = 1/2, 1 and 2 are
        the sections X = 1/2, 1 and 2, all on the half xi > 0 that is the front."""
        view, = [v for v in load("embedding")["views"] if v["id"] == "hyperbolic"]
        sections = []
        for curve in view["surfaces"][0]["curves"]:
            xs = {round(float(p[0]), 6) for p in curve["points"]}
            self.assertEqual(len(xs), 1)
            for X, Y, Z in ([float(c) for c in p] for p in curve["points"]):
                self.assertAlmostEqual((Z + 1) ** 2 - X * X - Y * Y, 1.0, places=3)
            sections.append((curve["class"], xs.pop()))
        self.assertEqual(sections, [("horizon", 0.0), ("flow", 0.5), ("flow", 1.0), ("flow", 2.0)])


@unittest.skipUnless(HAS_SYMPY, "needs sympy: /tmp/mfs-venv/bin/python -m unittest discover -s _tools")
class Physics(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import sympy as sp
        import verify_metrics as vm
        cls.sp, cls.vm = sp, vm
        cls.entry = json.loads((DATA / "metrics" / "kundt_waves.json").read_text(encoding="utf-8"))

    def chart(self, system_id):
        entry = next(c for c in self.entry["coordinates"] if c["id"] == system_id)
        reader = self.vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
        return entry, reader, [reader.symbol[c] for c in entry["coords"]]

    def components(self, entry, reader, field, variant=None):
        block = entry[field] if variant is None else entry[field]["variants"][variant]["nonzero"]
        return {tuple(e["indices"]): reader(e["value"]) for e in block}

    def test_the_wave_is_a_vacuum_exactly_where_its_profile_is_harmonic(self):
        sp = self.sp
        for system_id in ("kundt", "podolsky_belan"):
            entry, reader, (u, _, x, y) = self.chart(system_id)
            ricci = self.components(entry, reader, "ricci_tensor", "ll")
            G = reader.parameters["G"]
            self.assertEqual(set(ricci), {("u", "u")})
            self.assertEqual(sp.simplify(ricci[("u", "u")] + 2 * x * (sp.diff(G, x, 2) + sp.diff(G, y, 2))), 0)
        for system_id in ("simplest_wave", "kerr_schild"):
            entry, reader, _ = self.chart(system_id)
            self.assertEqual(self.components(entry, reader, "ricci_tensor", "ll"), {})

    def test_every_curvature_scalar_published_vanishes_and_the_curvature_does_not(self):
        """K = R = 0 in each chart without a cosmological constant, while the simplest wave has
        R_xuxu = -R_yuyu = -4x/l, which a falling observer feels."""
        for system_id in ("kundt", "podolsky_belan", "simplest_wave", "kerr_schild"):
            entry, _, _ = self.chart(system_id)
            self.assertEqual((entry["ricci_scalar"], entry["kretschmann"]), ("R = 0", "K = 0"))
        entry, reader, (u, v, x, y) = self.chart("simplest_wave")
        riemann = self.components(entry, reader, "riemann", "llll")
        ell = reader.parameters["ell"]
        self.assertEqual(self.sp.simplify(riemann[("x", "u", "x", "u")] + 4 * x / ell), 0)
        self.assertEqual(self.sp.simplify(riemann[("y", "u", "y", "u")] - 4 * x / ell), 0)

    def test_without_the_wave_the_metric_is_flat_space_outside_the_cone(self):
        """X = x(1 + 2uv), Z = x(v - u(1 + uv)), cT = x(v + u(1 + uv)) and Y = y pull Minkowski's
        metric back onto Podolsky and Belan's at G = 0, with X^2 + Z^2 - c^2T^2 = x^2."""
        sp = self.sp
        entry, reader, (u, v, x, y) = self.chart("podolsky_belan")
        published = self.components(entry, reader, "metric_components")
        names = entry["coords"]
        g = sp.Matrix(4, 4, lambda i, j: published.get((names[i], names[j]), 0)).subs(reader.parameters["G"], 0)
        image = [x * (v + u * (1 + u * v)), x * (1 + 2 * u * v), y, x * (v - u * (1 + u * v))]
        J = sp.Matrix(4, 4, lambda i, j: sp.diff(image[i], (u, v, x, y)[j]))
        self.assertEqual((J.T * sp.diag(-1, 1, 1, 1) * J - g).applyfunc(sp.simplify), sp.zeros(4, 4))
        self.assertEqual(sp.simplify(image[1] ** 2 + image[3] ** 2 - image[0] ** 2 - x ** 2), 0)
        # A front u = tan(alpha/2) is the null plane cT = X sin(alpha) + Z cos(alpha).
        front = (1 + u ** 2) * image[0] - 2 * u * image[1] - (1 - u ** 2) * image[3]
        self.assertEqual(sp.simplify(front), 0)

    def test_the_family_holds_the_scalars_of_its_cosmological_constant(self):
        entry, _, _ = self.chart("ozsvath_robinson_rozga")
        self.assertEqual(entry["ricci_scalar"], "R = 4\\Lambda")
        self.assertEqual(entry["kretschmann"], "K = \\dfrac{8\\Lambda^2}{3}")

    def test_the_declared_derivatives_of_the_kerr_schild_names_are_their_definitions(self):
        """The reader accepts the rates the checker declares, and their one-form for u is the
        null covector of the Kerr-Schild form over 2x."""
        sp, vm = self.sp, self.vm
        key = ("kundt_waves", "kerr_schild")
        entry = next(c for c in self.entry["coordinates"] if c["id"] == "kerr_schild")
        reader = vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]],
                           vm.time_coordinates(vm.DIMENSIONS[key], entry["coords"]), held=vm.HELD[key],
                           rates=vm.RATES[key])
        x, u = reader.parameters["x"], reader.parameters["u"]
        T, X, Z = (reader.symbol[n] for n in ("T", "X", "Z"))
        k = {T: reader.c * (1 + u ** 2), X: -2 * u, Z: -(1 - u ** 2)}
        for symbol, component in k.items():
            self.assertEqual(vm.norm(reader.by_rates(sp.diff(u, symbol)) - component / (2 * x)), 0)
        # A second derivative is written by the first: d^2x/dX^2 = (1 - (dx/dX)^2)/x.
        second = reader.by_rates(sp.diff(x, X, 2))
        self.assertFalse(second.atoms(sp.Derivative))
        self.assertEqual(vm.norm(second - (1 - reader.by_rates(sp.diff(x, X)) ** 2) / x), 0)


@unittest.skipUnless(HAS_SYMPY, "needs sympy: /tmp/mfs-venv/bin/python -m unittest discover -s _tools")
class Rates(unittest.TestCase):
    """The checker's declared derivatives of a held name, on a name small enough to read."""

    def setUp(self):
        import sympy as sp
        import verify_metrics as vm
        self.sp, self.vm = sp, vm
        self.parameters = ["r = \\sqrt{x^2 + y^2}"]

    def test_a_true_rate_is_taken_and_writes_every_derivative(self):
        sp, vm = self.sp, self.vm
        reader = vm.Reader(["x", "y"], self.parameters, (), held=("r",),
                           rates={"r": {"x": "\\dfrac{x}{r}", "y": "\\dfrac{y}{r}"}})
        r = reader.parameters["r"]
        x, y = reader.symbol["x"], reader.symbol["y"]
        laplacian = reader.surface(sp.diff(r, x, 2) + sp.diff(r, y, 2))
        self.assertFalse(laplacian.atoms(sp.Derivative))
        self.assertEqual(vm.norm(laplacian.subs(r, sp.sqrt(x ** 2 + y ** 2)) - 1 / sp.sqrt(x ** 2 + y ** 2)), 0)

    def test_a_false_rate_is_refused(self):
        with self.assertRaises(self.vm.LatexError):
            self.vm.Reader(["x", "y"], self.parameters, (), held=("r",),
                           rates={"r": {"x": "\\dfrac{y}{r}", "y": "\\dfrac{y}{r}"}})

    def test_rates_for_a_name_that_is_not_held_are_refused(self):
        with self.assertRaises(self.vm.LatexError):
            self.vm.Reader(["x", "y"], self.parameters, (), rates={"r": {"x": "\\dfrac{x}{r}", "y": "\\dfrac{y}{r}"}})


if __name__ == "__main__":
    unittest.main()
