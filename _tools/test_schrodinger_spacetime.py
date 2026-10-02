"""Schrodinger spacetime, anti-de Sitter space with one term in dt^2 added along a null direction:
python3 -m unittest discover -s _tools

What its texts and drawings state, held on the published files. The first class needs nothing
but the files: the light cones of every plane of the time and the null coordinate, one edge along
the null coordinate and the other leaning back by half the coefficient of -dt^2, the rays of the
global chart swinging in the trap, and the plane of x and r on its hyperboloid. The second reads
the published tensors through the checker's Reader and holds them to the matter the History
names, Son's massive vector field, which is Balasubramanian and McGreevy's dust, and to the
causal facts the History states; it needs sympy and is skipped where it is absent.
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

# The coefficient h of -c^2 dt^2 in the bracket on each plane drawn, at L = beta = 1 and
# omega = c/beta: beta^2/r^2, its image under rho = L^2/r, beta^2/R^2 + omega^2 R^2/c^2 on the
# axis of the trap, and (beta/r)^(2z - 2) at r = beta/sqrt(2).
PLANES = {
    ("poincare", "near"): 4.0, ("poincare", "middle"): 1.0, ("poincare", "far"): 0.25,
    ("inverse_radius", "near"): 4.0, ("inverse_radius", "middle"): 1.0, ("inverse_radius", "far"): 0.25,
    ("global", "near"): 4.25, ("global", "middle"): 2.0, ("global", "far"): 4.25,
    ("dynamical_exponent", "one"): 1.0, ("dynamical_exponent", "three_halves"): math.sqrt(2),
    ("dynamical_exponent", "three"): 4.0,
    ("poincare_5d", "middle"): 1.0, ("poincare_6d", "middle"): 1.0,
}


def load(folder):
    return json.loads((DATA / folder / "schrodinger_spacetime.json").read_text(encoding="utf-8"))


class Drawings(unittest.TestCase):
    def test_every_plane_is_drawn(self):
        drawn = {(system, view["id"]) for system, views in load("diagrams")["systems"].items() for view in views}
        self.assertEqual(drawn, set(PLANES))

    def test_one_edge_of_every_cone_runs_along_the_null_coordinate_at_one_time(self):
        """On a plane of the time and the null coordinate the metric is -(L^2/r^2)(2 dt dxi + h dt^2),
        so d/dxi is null and future directed, and the other null direction is dxi = -(h/2) dt."""
        for system, views in load("diagrams")["systems"].items():
            for view in views:
                h = PLANES[(system, view["id"])]
                self.assertEqual(view["to_display"], [[0.0, 1.0], [1.0, 0.0]])
                self.assertEqual(view["box"], [-2, 2, -2, 2])
                norm = math.hypot(h / 2, 1.0)
                for cone in view["cones"]:
                    self.assertAlmostEqual(cone["b"][0], 1.0, places=4)
                    self.assertAlmostEqual(cone["b"][1], 0.0, places=4)
                    self.assertAlmostEqual(cone["a"][0], -(h / 2) / norm, places=3)
                    self.assertAlmostEqual(cone["a"][1], 1 / norm, places=3)

    def test_the_cone_is_wider_the_nearer_the_boundary(self):
        """The angle between the two edges grows as r falls, toward the half plane after t."""
        by_view = {view["id"]: view["cones"][0] for view in load("diagrams")["systems"]["poincare"]}
        angle = {name: math.atan2(c["a"][1], c["a"][0]) for name, c in by_view.items()}
        self.assertLess(angle["far"], angle["middle"])
        self.assertLess(angle["middle"], angle["near"])
        self.assertLess(angle["near"], math.pi)

    def test_every_ray_keeps_its_time_or_its_sum(self):
        """A ray moving right keeps t, and a curve moving left keeps xi + (h/2) ct."""
        for system, views in load("diagrams")["systems"].items():
            for view in views:
                h = PLANES[(system, view["id"])]
                X0, X1, Y0, Y1 = view["box"]
                kept = {"P": lambda xi, t: xi + h * t / 2, "M": lambda xi, t: t}
                for family, rays in view["rays"].items():
                    for ray in rays:
                        values = [kept[family](X0 + x * (X1 - X0), Y0 + y * (Y1 - Y0)) for x, y in ray]
                        self.assertLess(max(values) - min(values), 2e-3, f"{system}/{view['id']}/{family}")

    def test_light_swings_in_the_trap_and_never_reaches_the_boundary(self):
        """A null geodesic of the surface X = 0 launched at T = 0 from the depth R_0 with no velocity
        along R has R^2 = R_0^2 cos^2 T + sin^2 T/R_0^2 at beta = 1 and omega = c/beta: it swings
        between R_0 and 1/R_0 once in every pi of T, and the one launched at R_0 = 1 keeps its depth."""
        figure, = load("diagrams")["projections"]["global"]
        self.assertEqual(figure["id"], "trap")
        self.assertNotIn("turn", figure)
        rays = [layer for layer in figure["layers"] if layer["class"] in ("above", "below")]
        self.assertEqual(len(rays), 6)
        launched = []
        for ray in rays:
            points = ray["points"]
            if ray["class"] == "below":
                # A straight line, which the file writes by its two ends.
                R0 = points[0][1]
            else:
                T0, R0 = min(points, key=lambda p: abs(p[0]))
                self.assertAlmostEqual(T0, 0.0, delta=0.02)
            launched.append(round(R0, 2))
            self.assertAlmostEqual(points[0][0], -math.pi, places=3)
            self.assertAlmostEqual(points[-1][0], math.pi, places=3)
            for T, R in points:
                self.assertAlmostEqual(R * R, R0 * R0 * math.cos(T) ** 2 + math.sin(T) ** 2 / (R0 * R0), delta=4e-3)
                self.assertGreaterEqual(R, min(R0, 1 / R0) - 1e-3)
            if ray["class"] == "below":
                self.assertTrue(all(abs(R - 1.0) < 1e-4 for _, R in points))
        self.assertEqual(sorted(launched), [0.3, 0.35, 0.45, 0.6, 0.8, 1.0])
        edges = sorted(layer["points"][0][0] for layer in figure["layers"] if layer["class"] == "edge")
        self.assertAlmostEqual(edges[0], -math.pi / 2, places=3)
        self.assertAlmostEqual(edges[1], math.pi / 2, places=3)

    def test_the_plane_of_x_and_r_is_a_sheet_of_a_hyperboloid(self):
        """(dx^2 + dr^2)/r^2 swept about x = 0, r = L: the circle a proper distance s = ln(r/L) out
        along x = 0 has the radius sinh(s) and stands at the height cosh(s) - 1, and each line of
        constant r drawn on it lies on that sheet."""
        view, = load("embedding")["views"]
        self.assertNotIn("system", view)
        surface, = view["surfaces"]
        sheet = next(p for p in surface["pieces"] if p["id"] == "sheet")
        self.assertEqual((sheet["system"], sheet["coordinate"], sheet["space"]), ("poincare", "r", "minkowski"))
        for r, rho, z in sheet["points"]:
            s = math.log(r)
            self.assertAlmostEqual(rho, math.sinh(s), places=5)
            self.assertAlmostEqual(z, math.cosh(s) - 1, places=5)
        self.assertAlmostEqual(math.log(sheet["points"][-1][0]), 2.0, places=9)
        curves = surface["curves"]
        self.assertEqual(len(curves), 3)
        depths = []
        for curve in curves:
            for X, Y, Z in curve["points"]:
                self.assertAlmostEqual((Z + 1) ** 2 - X * X - Y * Y, 1.0, delta=2e-5)
            # On the hyperboloid T = Z + 1, r = 1/(T - X) with X along the half line x = 0 toward larger r.
            depth = [1 / (Z + 1 - X) for X, Y, Z in curve["points"]]
            self.assertLess(max(depth) - min(depth), 2e-5 * max(depth))
            depths.append(round(depth[0], 3))
        self.assertEqual(sorted(depths), [0.5, 1.0, 2.0])


@unittest.skipUnless(HAS_SYMPY, "needs sympy: run with /tmp/mfs-venv/bin/python")
class FieldEquations(unittest.TestCase):
    """The published tensors against the matter and the causal facts the History states."""

    @classmethod
    def setUpClass(cls):
        import sympy as sp
        import verify_metrics as vm
        cls.sp, cls.vm = sp, vm
        cls.charts = {}
        for entry in load("metrics")["coordinates"]:
            reader = vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
            cls.charts[entry["id"]] = (entry, reader)

    def matrix(self, system, field, variant=None):
        entry, reader = self.charts[system]
        names = entry["coords"]
        items = entry[field] if variant is None else entry[field]["variants"][variant]["nonzero"]
        there = {tuple(e["indices"]): e["value"] for e in items}
        n = len(names)
        return self.sp.Matrix(n, n, lambda i, j: reader(there.get((names[i], names[j]), "0")))

    def test_the_source_is_dust_along_the_null_direction_beside_a_negative_cosmological_constant(self):
        """G_ab + Lambda g_ab with Lambda = -(d + 1)(d + 2)/2L^2 vanishes in every slot but the one
        of the time twice, where it is (z - 1)(2z + d) h/r^2: 4 + d for Schrodinger spacetime, and
        Balasubramanian and McGreevy's 2z^2 + z - 3 for d = 3."""
        sp, vm = self.sp, self.vm
        for system, (entry, reader) in self.charts.items():
            P = reader.parameters
            n = len(entry["coords"])
            d = n - 3
            L, beta = P["L"], P["beta"]
            radial = reader.symbol[entry["coords"][-1]]
            r = L ** 2 / radial if system == "inverse_radius" else radial
            z = P["z"] if "z" in P else sp.Integer(2)
            g = self.matrix(system, "metric_components")
            G = self.matrix(system, "einstein_tensor", "ll")
            Lam = -sp.Rational((d + 1) * (d + 2), 2) / L ** 2
            for a in range(n):
                for b in range(n):
                    want = (z - 1) * (2 * z + d) * (beta / r) ** (2 * z - 2) / r ** 2 if a == b == 0 else 0
                    self.assertEqual(vm.norm(G[a, b] + Lam * g[a, b] - want), 0, f"{system} {a}{b}")
        self.assertEqual(sp.expand((z - 1) * (2 * z + 3) - (2 * z ** 2 + z - 3)), 0)

    def test_no_chart_has_a_time_whose_gradient_is_timelike(self):
        """g^tt = 0 in every chart: the surfaces of constant time are null, d/dxi lies in them and
        is null, and the Killing vector of the time is timelike everywhere."""
        for system, (entry, reader) in self.charts.items():
            g = self.matrix(system, "metric_components")
            inverse = self.matrix(system, "inverse_metric_components")
            self.assertEqual(inverse[0, 0], 0, system)
            self.assertEqual(g[1, 1], 0, system)
            values = {s: self.sp.Rational(3, 2) for s in g[0, 0].free_symbols}
            self.assertLess(g[0, 0].subs(values), 0, system)

    def test_every_curvature_scalar_published_is_a_constant(self):
        """Homogeneous: R = -(d + 2)(d + 3)/L^2 and K = 2(d + 2)(d + 3)/L^4, anti-de Sitter's values."""
        for system, (entry, reader) in self.charts.items():
            n = len(entry["coords"])
            L = reader.parameters["L"]
            R = reader(entry["ricci_scalar"].partition("=")[2])
            K = reader(entry["kretschmann"].partition("=")[2])
            self.assertEqual(self.sp.simplify(R + n * (n - 1) / L ** 2), 0, system)
            self.assertEqual(self.sp.simplify(K - 2 * n * (n - 1) / L ** 4), 0, system)

    def test_light_launched_along_the_leaning_edge_is_turned_away_from_the_boundary(self):
        """The captions' accelerations, from the published Christoffel symbols: on the null direction
        dxi = -(h/2) dt the radial equation leaves r'' = (beta^2/r^3) t'^2 in the Poincare chart and
        R'' = (beta^2/R^3 - omega^2 R/c^2) T'^2 on the axis of the trap, which vanishes at
        R^2 = beta c/omega."""
        sp = self.sp
        for system, want in (("poincare", lambda P, r, c: P["beta"] ** 2 / r ** 3),
                             ("global", lambda P, r, c: P["beta"] ** 2 / r ** 3 - P["omega"] ** 2 * r / c ** 2)):
            entry, reader = self.charts[system]
            names = entry["coords"]
            gamma = {tuple(names.index(i) for i in e["indices"]): reader(e["value"])
                     for e in entry["christoffel"]["variants"]["ull"]["nonzero"]}
            g = self.matrix(system, "metric_components")
            k = [1, -g[0, 0] / (2 * g[0, 1]), 0, 0]
            self.assertEqual(sp.simplify(sum(g[a, b] * k[a] * k[b] for a in range(4) for b in range(4))), 0)
            r = reader.symbol[names[3]]
            push = -sum(gamma.get((3, a, b), 0) * k[a] * k[b] for a in range(4) for b in range(4))
            on_axis = {reader.symbol[names[2]]: 0}
            self.assertEqual(sp.simplify((push - want(reader.parameters, r, reader.c)).subs(on_axis)), 0, system)
            sideways = -sum(gamma.get((2, a, b), 0) * k[a] * k[b] for a in range(4) for b in range(4))
            self.assertEqual(sp.simplify(sideways.subs(on_axis)), 0, system)


if __name__ == "__main__":
    unittest.main()
