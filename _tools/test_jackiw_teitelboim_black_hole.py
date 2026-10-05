"""The black hole of Jackiw and Teitelboim's gravity in two dimensions:
.venv.noindex/bin/python -m unittest discover -s _tools

What its texts state, held on the published files. On the checker's own Geometry: every chart has
R = -2/L^2, the field equation of the dilaton, and no Einstein or Weyl tensor; the dilaton each
convention states solves the equation the metric imposes on it and carries the one mass
phi^2 - L^2 (D phi)^2 = phi_r^2 r_h^2/L^4; the surface gravity is r_h/L^2, the temperature the
parameter's description states; the static chart pulled back through each map is the chart it is
carried into; and the plane of t and r of the BTZ black hole without rotation is the static chart
at L = l and r_h^2 = M l^2, the reduction of Achucarro and Ortiz. On the drawings: nothing is drawn
beyond the lines UV = +-1 of the Kruskal plane or the hyperbola where the dilaton vanishes on the
Poincare plane, the rays keep their null coordinates, the conformal square has its boundaries on
X = +-pi/2 and the lines where the dilaton vanishes on T = +-pi/2, and the Euclidean disc is the
sheet Z = L cosh(rho/L) - L of a hyperboloid. The first tests need sympy and are skipped where it
is absent; the rest read the files alone.
"""
import importlib.util
import json
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
ID = "jackiw_teitelboim_black_hole"
HAS_SYMPY = importlib.util.find_spec("sympy") is not None
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))


def published(name, folder="metrics"):
    return json.loads((DATA / folder / f"{name}.json").read_text(encoding="utf-8"))


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class Charts(unittest.TestCase):
    def setUp(self):
        import sympy
        import verify_metrics
        self.sp, self.vm = sympy, verify_metrics

    def geometry(self, system, metric_id=ID):
        metric = published(metric_id)
        entry = next(c for c in metric["coordinates"] if c["id"] == system)
        declared = self.vm.DIMENSIONS[(metric_id, system)]
        reader = self.vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]],
                                self.vm.time_coordinates(declared, entry["coords"]))
        # The line element carries c^2 dt^2; with c = 1 the time coordinate is x^0 itself, as the
        # conventions take it.
        g = self.vm.metric_from_line_element(reader, entry["line_element"], entry["coords"]).subs(reader.c, 1)
        return reader, self.vm.Geometry(g, [reader.symbol[name] for name in entry["coords"]], 60)

    def dilaton(self, system, reader):
        """The dilaton each chart's convention states, with x^0 = ct or cT, divided by phi_r."""
        sp, P, s = self.sp, reader.parameters, reader.symbol
        L, rh = P["L"], P["r_h"]
        a, b = (s[name] for name in next(c["coords"] for c in published(ID)["coordinates"] if c["id"] == system))
        return {
            "static": b / L ** 2,
            "proper_distance": rh * sp.cosh(b / L) / L ** 2,
            "kruskal": rh * (1 - a * b) / (L ** 2 * (1 + a * b)),
            "global": rh * sp.cos(a) / (L ** 2 * sp.sin(b)),
            "poincare": (1 - rh ** 2 * (a ** 2 - b ** 2) / (4 * L ** 4)) / b,
        }[system]

    def test_every_chart_has_the_curvature_the_dilaton_imposes_and_nothing_else(self):
        for system in ("static", "proper_distance", "kruskal", "global", "poincare"):
            reader, geo = self.geometry(system)
            L = reader.parameters["L"]
            self.assertEqual(self.vm.norm(geo.ricci_scalar() + 2 / L ** 2), 0, system)
            self.assertEqual(self.vm.norm(geo.kretschmann() - 4 / L ** 4), 0, system)
            for tensor, rank in ((geo.weyl_llll(), 4), (geo.einstein_ll(), 2)):
                for index in self.vm._indices(2, rank):
                    self.assertEqual(self.vm._at(tensor, index), 0, f"{system} {index}")

    def test_each_stated_dilaton_solves_its_equation_and_is_the_one_black_hole(self):
        sp = self.sp
        for system in ("static", "proper_distance", "kruskal", "global", "poincare"):
            reader, geo = self.geometry(system)
            L, rh = reader.parameters["L"], reader.parameters["r_h"]
            xs = [reader.symbol[name] for name in next(c["coords"] for c in published(ID)["coordinates"]
                                                       if c["id"] == system)]
            phi = self.dilaton(system, reader)
            gamma, g, ginv = geo.christoffel_ull(), geo.g, geo.ginv
            d = [sp.diff(phi, x) for x in xs]
            hessian = [[sp.diff(phi, xs[i], xs[j]) - sum(gamma[k][i][j] * d[k] for k in range(2))
                        for j in range(2)] for i in range(2)]
            box = sum(ginv[i, j] * hessian[i][j] for i in range(2) for j in range(2))
            for i in range(2):
                for j in range(2):
                    gap = hessian[i][j] - g[i, j] * box + g[i, j] * phi / L ** 2
                    self.assertEqual(sp.simplify(gap.rewrite(sp.exp)), 0, f"{system} {(i, j)}")
            mass = phi ** 2 - L ** 2 * sum(ginv[i, j] * d[i] * d[j] for i in range(2) for j in range(2))
            self.assertEqual(sp.simplify((mass - rh ** 2 / L ** 4).rewrite(sp.exp)), 0, system)

    def test_the_surface_gravity_is_the_temperature_the_description_states(self):
        # f = (r^2 - r_h^2)/L^2 vanishes at r_h, and f'(r_h)/2 = r_h/L^2, so the temperature is
        # hbar c r_h/(2 pi k_B L^2).
        reader, geo = self.geometry("static")
        r, L, rh = reader.symbol["r"], reader.parameters["L"], reader.parameters["r_h"]
        f = -geo.g[0, 0]
        self.assertEqual(self.sp.simplify(f.subs(r, rh)), 0)
        self.assertEqual(self.sp.simplify(self.sp.diff(f, r).subs(r, rh) / 2 - rh / L ** 2), 0)
        (static,) = [c for c in published(ID)["coordinates"] if c["id"] == "static"]
        (described,) = [p["description"] for p in static["parameters"] if p["symbol"] == "r_h"]
        self.assertIn("\\hbar c\\,r_h/2\\pi k_BL^2", described)

    def test_the_static_chart_carried_into_each_other_chart_is_that_chart(self):
        # Through the hyperboloid -X_0^2 - X_1^2 + X_2^2 = -L^2: the static chart's X_0 = L r/r_h,
        # X_1 = L sqrt(r^2 - r_h^2) sinh(k x^0)/r_h and X_2 = L sqrt(r^2 - r_h^2) cosh(k x^0)/r_h,
        # with k = r_h/L^2; evaluated at three events of the exterior, to twelve places.
        sp = self.sp
        events = ((0.3, 2.5, 1.3, 1.1), (-1.2, 1.4, 0.8, 0.6), (2.0, 7.0, 2.0, 3.0))
        for system in ("proper_distance", "kruskal", "global", "poincare"):
            reader, geo = self.geometry(system)
            L, rh = reader.parameters["L"], reader.parameters["r_h"]
            xs = [reader.symbol[name] for name in next(c["coords"] for c in published(ID)["coordinates"]
                                                       if c["id"] == system)]
            T, R = sp.symbols("T R", positive=True)
            k = rh / L ** 2
            root = sp.sqrt((R - rh) / (R + rh))
            X0, X1 = L * R / rh, L * sp.sqrt(R ** 2 - rh ** 2) * sp.sinh(k * T) / rh
            X2 = L * sp.sqrt(R ** 2 - rh ** 2) * sp.cosh(k * T) / rh
            image = {
                "proper_distance": (T, L * sp.acosh(R / rh)),
                "kruskal": (-sp.exp(-k * T) * root, sp.exp(k * T) * root),
                "global": (sp.atan2(X1, X0), sp.pi - sp.acot(X2 / L)),
                "poincare": (2 * L ** 2 * X1 / (rh * (X0 + X2)), 2 * L ** 3 / (rh * (X0 + X2))),
            }[system]
            J = sp.Matrix([[sp.diff(f, s) for s in (T, R)] for f in image])
            pulled = J.T * geo.g.subs(dict(zip(xs, image)), simultaneous=True) * J
            for t0, r0, l0, h0 in events:
                at = {T: t0, R: r0, L: l0, rh: h0}
                want = [[-((r0 ** 2 - h0 ** 2) / l0 ** 2), 0], [0, l0 ** 2 / (r0 ** 2 - h0 ** 2)]]
                for i in range(2):
                    for j in range(2):
                        got = sp.N(pulled[i, j].subs(at), 30)
                        self.assertAlmostEqual(float(got), want[i][j], places=12, msg=f"{system} {(i, j)} at t = {t0}, r = {r0}")

    def test_btz_without_rotation_reduces_to_the_static_chart(self):
        # Achucarro and Ortiz: the plane of t and r of BTZ at J = 0 is the static chart with L = l and
        # r_h^2 = M l^2, and the circle of phi has radius r, the dilaton phi_r r/L^2 at phi_r = L^2.
        sp = self.sp
        btz_reader, btz = self.geometry("stationary", "btz")
        reader, jt = self.geometry("static")
        ell, M, J = (btz_reader.parameters[name] for name in ("ell", "M", "J"))
        L, rh = reader.parameters["L"], reader.parameters["r_h"]
        r_btz, r_jt = btz_reader.symbol["r"], reader.symbol["r"]
        for i in range(2):
            for j in range(2):
                plane = btz.g[i, j].subs(J, 0)
                reduced = jt.g[i, j].subs({L: ell, rh: sp.sqrt(M) * ell, r_jt: r_btz})
                self.assertEqual(sp.simplify(plane - reduced), 0, (i, j))
        self.assertEqual(sp.simplify(btz.g[2, 2].subs(J, 0) - r_btz ** 2), 0)


class TheDrawings(unittest.TestCase):
    def view(self, system):
        (view,) = published(ID, "diagrams")["systems"][system]
        return view

    def chart(self, view, u):
        """A point of the unit square in the drawn axes, X across and Y up."""
        X0, X1, Y0, Y1 = view["box"]
        return X0 + u[0] * (X1 - X0), Y0 + u[1] * (Y1 - Y0)

    def drawn(self, view):
        points = [u for family in view["rays"].values() for line in family for u in line]
        return points + [cone["at"] for cone in view["cones"]]

    def test_nothing_on_the_kruskal_plane_lies_beyond_the_boundaries_or_where_the_dilaton_vanishes(self):
        view = self.view("kruskal")
        points = self.drawn(view)
        self.assertGreater(len(points), 30)
        for u in points:
            X, Y = self.chart(view, u)
            U, V = Y - X, Y + X
            self.assertLess(abs(U * V), 1 + 2e-2, u)

    def test_nothing_on_the_poincare_plane_lies_beyond_the_hyperbola_where_the_dilaton_vanishes(self):
        # At L = r_h = 1 the dilaton vanishes on c^2 T^2 - z^2 = 4.
        view = self.view("poincare")
        for u in self.drawn(view):
            z, cT = self.chart(view, u)
            self.assertLess(cT * cT - z * z, 4 + 5e-2, u)

    def test_every_ray_keeps_its_null_coordinate(self):
        # The Kruskal plane keeps U or V; the global and Poincare charts are conformally flat, so a
        # ray keeps tau +- sigma or cT +- z.
        for system, keep in (("kruskal", lambda X, Y: (Y + X, Y - X)), ("global", lambda X, Y: (Y + X, Y - X)),
                             ("poincare", lambda X, Y: (Y + X, Y - X))):
            view = self.view(system)
            for name, which in (("P", 0), ("M", 1)):
                for line in view["rays"][name]:
                    held = [keep(*self.chart(view, u))[which] for u in line]
                    self.assertLess(max(held) - min(held), 5e-3, f"{system} {name}")

    def test_the_static_view_marks_the_horizon_at_r_h(self):
        view = self.view("static")
        (horizon,) = [m for m in view["markers"] if m["kind"] == "grr"]
        for line in horizon["lines"]:
            for u in line:
                self.assertAlmostEqual(self.chart(view, u)[0], 1.0, delta=1e-3)

    def test_the_conformal_square_has_its_boundaries_and_its_ends_where_the_text_puts_them(self):
        for view in published(ID, "conformal")["views"]:
            lines = {cls: [l["points"] for l in view["layers"] if l["class"] == cls]
                     for cls in ("boundary", "singular", "horizon")}
            self.assertEqual(len(lines["boundary"]), 2, view["id"])
            for points in lines["boundary"]:
                for X, _ in points:
                    self.assertAlmostEqual(abs(X), math.pi / 2, places=3)
            self.assertEqual(len(lines["singular"]), 2, view["id"])
            for points in lines["singular"]:
                for _, T in points:
                    self.assertAlmostEqual(abs(T), math.pi / 2, places=3)
            for points in lines["horizon"]:
                for X, T in points:
                    self.assertAlmostEqual(abs(X), abs(T), places=3)

    def test_the_euclidean_disc_is_a_sheet_of_a_hyperboloid(self):
        # rho_circle = L sinh(rho/L) and Z = L cosh(rho/L) - L at L = 1, in Minkowski space, with the
        # light cone Z = sinh(rho) - 1 beside it.
        (view,) = published(ID, "embedding")["views"]
        pieces = {p["id"]: p for p in view["surfaces"][0]["pieces"]}
        disc, cone = pieces["disc"], pieces["cone"]
        self.assertEqual((disc["system"], disc["coordinate"], disc["space"]), ("proper_distance", "\\rho", "minkowski"))
        for rho, radius, z in disc["points"]:
            self.assertAlmostEqual(radius, math.sinh(rho), places=6)
            self.assertAlmostEqual(z, math.cosh(rho) - 1, places=5)
        self.assertEqual(disc["points"][0], [0, 0, 0])
        for rho, radius, z in cone["points"]:
            self.assertAlmostEqual(z, radius - 1, places=6)


if __name__ == "__main__":
    unittest.main()
