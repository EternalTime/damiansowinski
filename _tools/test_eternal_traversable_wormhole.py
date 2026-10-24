"""Maldacena and Qi's eternal traversable wormhole in two dimensions:
.venv.noindex/bin/python -m unittest discover -s _tools

What its texts state, held on the published files. On the checker's own Geometry: every chart has
R = -2/L^2 and no Einstein or Weyl tensor; the dilaton each convention states solves the equation the
metric imposes on it, D_a D_b phi - g_ab D^2 phi + g_ab phi/L^2 = -T_ab, with the traceless
T_tt = T_sigmasigma = -phi_h of the global chart; it is least at sigma = pi/2 and grows as
pi phi_h/2 sigma toward the boundary; a light ray crosses the strip in the global time pi; and the
global chart carried into the Poincare and thermal charts is each of them, metric and dilaton both.
On the drawings: the rays keep their null coordinates, the rays through the middle run from boundary
to boundary, the conformal strip has its boundaries on X = +-pi/2 and its middle on X = 0, and the
Euclidean cylinder has the radius L/sin(sigma) from sigma = pi/4 to 3 pi/4. The first tests need
sympy and are skipped where it is absent; the rest read the files alone.
"""
import importlib.util
import json
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
ID = "eternal_traversable_wormhole"
SYSTEMS = ("global", "poincare", "thermal")
HAS_SYMPY = importlib.util.find_spec("sympy") is not None
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))


def published(name, folder="metrics"):
    return json.loads((DATA / folder / f"{name}.json").read_text(encoding="utf-8"))


def chart_entry(system):
    return next(c for c in published(ID)["coordinates"] if c["id"] == system)


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class Charts(unittest.TestCase):
    def setUp(self):
        import sympy
        import verify_metrics
        self.sp, self.vm = sympy, verify_metrics

    def geometry(self, system):
        entry = chart_entry(system)
        declared = self.vm.DIMENSIONS[(ID, system)]
        reader = self.vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]],
                                self.vm.time_coordinates(declared, entry["coords"]))
        g = self.vm.metric_from_line_element(reader, entry["line_element"], entry["coords"])
        xs = [reader.symbol[name] for name in entry["coords"]]
        return reader, xs, self.vm.Geometry(g, xs, 60)

    def dilaton(self, system, reader, xs):
        """The dilaton each chart's convention states, phi_h (1 + x arctan x) with x = cot(sigma),
        (1 + t_P^2 - z^2)/2z and sinh(rho) cosh(t_R)."""
        sp, phi_h = self.sp, reader.parameters["phi_h"]
        a, b = xs
        if system == "global":
            return phi_h * (1 + (sp.pi / 2 - b) * sp.cot(b))
        x = {"poincare": (1 + a ** 2 - b ** 2) / (2 * b), "thermal": sp.sinh(b) * sp.cosh(a)}[system]
        return phi_h * (1 + x * sp.atan(x))

    def source(self, phi, xs, geo, L):
        sp = self.sp
        gamma, g, ginv = geo.christoffel_ull(), geo.g, geo.ginv
        d = [sp.diff(phi, x) for x in xs]
        hessian = sp.Matrix(2, 2, lambda i, j: sp.diff(phi, xs[i], xs[j]) - sum(gamma[k][i][j] * d[k] for k in range(2)))
        box = sum(ginv[i, j] * hessian[i, j] for i in range(2) for j in range(2))
        return hessian - g * box + g * phi / L ** 2

    def test_every_chart_has_the_curvature_of_anti_de_sitter_space_and_nothing_else(self):
        for system in SYSTEMS:
            reader, _, geo = self.geometry(system)
            L = reader.parameters["L"]
            self.assertEqual(self.vm.norm(geo.ricci_scalar() + 2 / L ** 2), 0, system)
            self.assertEqual(self.vm.norm(geo.kretschmann() - 4 / L ** 4), 0, system)
            for tensor, rank in ((geo.weyl_llll(), 4), (geo.einstein_ll(), 2)):
                for index in self.vm._indices(2, rank):
                    self.assertEqual(self.vm._at(tensor, index), 0, f"{system} {index}")

    def test_the_conventions_state_the_dilaton_the_tests_use(self):
        self.assertIn("\\phi = \\phi_h\\left(1 + (\\pi/2 - \\sigma)\\cot\\sigma\\right)", chart_entry("global")["convention"])
        self.assertIn("T_{tt} = T_{\\sigma\\sigma} = -\\phi_h", chart_entry("global")["convention"])
        self.assertIn("x = (1 + t_P^2 - z^2)/2z", chart_entry("poincare")["convention"])
        self.assertIn("x = \\sinh\\rho\\cosh t_R", chart_entry("thermal")["convention"])
        self.assertIn("= -T_{\\mu\\nu}", published(ID)["convention"])

    def test_the_global_dilaton_is_sourced_by_a_constant_negative_energy(self):
        sp = self.sp
        reader, xs, geo = self.geometry("global")
        L, phi_h = reader.parameters["L"], reader.parameters["phi_h"]
        E = self.source(self.dilaton("global", reader, xs), xs, geo, L)
        for i in range(2):
            for j in range(2):
                self.assertEqual(sp.simplify(E[i, j] - (phi_h if i == j else 0)), 0, (i, j))

    def test_every_stated_dilaton_has_a_traceless_source(self):
        sp = self.sp
        for system in SYSTEMS:
            reader, xs, geo = self.geometry(system)
            L, phi_h = reader.parameters["L"], reader.parameters["phi_h"]
            E = self.source(self.dilaton(system, reader, xs), xs, geo, L)
            trace = sum(geo.ginv[i, j] * E[i, j] for i in range(2) for j in range(2))
            for a0, b0 in ((0.3, 0.7), (-1.1, 1.9)):
                at = {xs[0]: a0, xs[1]: b0, L: 1.3, phi_h: 2.0}
                self.assertAlmostEqual(float(sp.N(trace.subs(at), 30)), 0.0, places=12, msg=system)

    def test_the_dilaton_is_least_in_the_middle_and_grows_toward_the_boundary(self):
        sp = self.sp
        reader, xs, _ = self.geometry("global")
        sigma, phi_h = xs[1], reader.parameters["phi_h"]
        phi = self.dilaton("global", reader, xs)
        self.assertEqual(sp.simplify(sp.diff(phi, sigma).subs(sigma, sp.pi / 2)), 0)
        self.assertEqual(sp.simplify(phi.subs(sigma, sp.pi / 2) - phi_h), 0)
        self.assertEqual(sp.limit(sigma * phi, sigma, 0, "+"), sp.pi * phi_h / 2)

    def test_a_light_ray_crosses_the_strip_in_the_global_time_pi(self):
        sp = self.sp
        reader, xs, geo = self.geometry("global")
        speed = sp.sqrt(sp.simplify(geo.g[1, 1] / -geo.g[0, 0]))
        self.assertEqual(sp.integrate(1 / speed, (xs[1], 0, sp.pi)), sp.pi)

    def test_the_global_chart_carried_into_the_others_is_each_of_them(self):
        # The Poincare chart's t_P +- z = tan((t +- sigma)/2) and the thermal chart's
        # tanh(t_R) = sin(t)/cos(sigma), cosh(rho) = cos(t)/sin(sigma), at events of both, to twelve places.
        sp = self.sp
        T, S = sp.symbols("T S", positive=True)
        maps = {
            "poincare": ((sp.tan((T + S) / 2) + sp.tan((T - S) / 2)) / 2,
                         (sp.tan((T + S) / 2) - sp.tan((T - S) / 2)) / 2),
            "thermal": (sp.atanh(sp.sin(T) / sp.cos(S)), sp.acosh(sp.cos(T) / sp.sin(S))),
        }
        _, gxs, _ = self.geometry("global")
        for system, image in maps.items():
            reader, xs, geo = self.geometry(system)
            L = reader.parameters["L"]
            here = dict(zip(xs, image))
            J = sp.Matrix([[sp.diff(f, s) for s in (T, S)] for f in image])
            pulled = J.T * geo.g.subs(here, simultaneous=True) * J
            phi = self.dilaton(system, reader, xs).subs(here, simultaneous=True)
            want = self.dilaton("global", reader, gxs).subs({gxs[0]: T, gxs[1]: S})
            for t0, s0 in ((0.1, 0.6), (-0.3, 1.0), (0.2, 0.4)):
                at = {T: t0, S: s0, L: 1.3, reader.parameters["phi_h"]: 2.0}
                for i in range(2):
                    for j in range(2):
                        expect = (-1 if i == 0 else 1) * 1.3 ** 2 / math.sin(s0) ** 2 if i == j else 0.0
                        self.assertAlmostEqual(float(sp.N(pulled[i, j].subs(at), 30)), expect, places=12,
                                               msg=f"{system} {(i, j)}")
                self.assertAlmostEqual(float(sp.N((phi - want).subs(at), 30)), 0.0, places=12, msg=system)


class TheDrawings(unittest.TestCase):
    def view(self, system):
        (view,) = published(ID, "diagrams")["systems"][system]
        return view

    def chart(self, view, u):
        """A point of the unit square in the drawn axes, X across and Y up."""
        X0, X1, Y0, Y1 = view["box"]
        return X0 + u[0] * (X1 - X0), Y0 + u[1] * (Y1 - Y0)

    def test_every_ray_keeps_its_null_coordinate(self):
        # The global and Poincare charts are conformally flat, so a ray keeps t +- sigma or t_P +- z; in
        # the thermal chart it keeps t_R -+ ln tanh(rho/2).
        keeps = {"global": lambda x, y: (y + x, y - x), "poincare": lambda x, y: (y + x, y - x),
                 "thermal": lambda x, y: (y + math.log(math.tanh(x / 2)), y - math.log(math.tanh(x / 2)))}
        for system, keep in keeps.items():
            view = self.view(system)
            for name, which in (("P", 0), ("M", 1)):
                for line in view["rays"][name]:
                    held = [keep(*self.chart(view, u))[which] for u in line if self.chart(view, u)[0] > 0.05]
                    self.assertLess(max(held) - min(held), 5e-3, f"{system} {name}")

    def test_the_rays_through_the_middle_run_from_boundary_to_boundary(self):
        view = self.view("global")
        (marked,) = [m for m in view["markers"] if m["kind"] == "event"]
        ends = sorted(round(self.chart(view, u)[0], 2) for line in marked["lines"] for u in (line[0], line[-1]))
        self.assertEqual(ends, [0.0, 0.0, round(math.pi, 2), round(math.pi, 2)])
        for line in marked["lines"]:
            for u in line:
                sigma, t = self.chart(view, u)
                self.assertAlmostEqual(abs(t), abs(sigma - math.pi / 2), delta=5e-3)

    def test_the_conformal_strip_has_its_boundaries_and_its_middle_where_the_text_puts_them(self):
        for view in published(ID, "conformal")["views"]:
            lines = {cls: [l["points"] for l in view["layers"] if l["class"] == cls]
                     for cls in ("boundary", "throat", "event")}
            self.assertEqual(len(lines["boundary"]), 2, view["id"])
            for points in lines["boundary"]:
                for X, _ in points:
                    self.assertAlmostEqual(abs(X), math.pi / 2, places=3)
            for points in lines["throat"]:
                for X, _ in points:
                    self.assertAlmostEqual(X, 0.0, places=3)
            for points in lines["event"]:
                for X, T in points:
                    self.assertAlmostEqual(abs(X), abs(T), places=3)

    def test_the_euclidean_cylinder_has_the_radius_l_over_sin_sigma_out_to_its_reach(self):
        (view,) = published(ID, "embedding")["views"]
        (piece,) = view["surfaces"][0]["pieces"]
        self.assertEqual((piece["system"], piece["coordinate"]), ("global", "\\sigma"))
        points = piece["points"]
        self.assertAlmostEqual(points[0][0], math.pi / 4, places=9)
        self.assertAlmostEqual(points[-1][0], 3 * math.pi / 4, places=9)
        for sigma, radius, z in points:
            self.assertAlmostEqual(radius, 1 / math.sin(sigma), places=6)
        self.assertAlmostEqual(min(p[1] for p in points), 1.0, places=6)
        self.assertAlmostEqual(points[0][2], -points[-1][2], places=6)


if __name__ == "__main__":
    unittest.main()
