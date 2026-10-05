"""The three-brane and its throat, the extreme three-brane of type IIB supergravity in ten dimensions:
.venv.noindex/bin/python -m unittest discover -s _tools

What its texts state, held on the published files: the Ricci scalar vanishes and the mixed Ricci
tensor is -4 L^8/(rho^4 + L^4)^(5/2) on the brane's four dimensions and the radial one and the
opposite on the sphere, which is Duff and Lu's field equation R_MN = F_MPQRS F_N^PQRS/96 for the
five-form of the potential A_txyz = 1/H; the areal chart and Gibbons, Horowitz and Townsend's chart
are the isotropic chart carried along rho^4 = r^4 - L^4 and w = rho/r, and the second is even in w;
the inversion rho -> L^2/rho multiplies the metric by L^2/rho^2; the throat is the limit of the
isotropic chart toward rho = 0, has no Weyl tensor, the Kretschmann scalar 80/L^4 that the
three-brane has on its horizon, and the Ricci tensor of anti-de Sitter space of radius L times a
sphere of radius L; and the drawings hold what their captions say. The tests of the published
mathematics need sympy and are skipped where it is absent; the tests of the drawings read the files
alone.
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
CHARTS = ["isotropic", "areal", "horizon", "throat", "throat_proper"]


def published(metric_id="three_brane_throat"):
    return json.loads((DATA / "metrics" / f"{metric_id}.json").read_text(encoding="utf-8"))


def drawing(kind, metric_id="three_brane_throat"):
    return json.loads((DATA / kind / f"{metric_id}.json").read_text(encoding="utf-8"))


def tortoise(rho, steps=4000):
    """rho_* at L = 1 with the constant that makes it odd in rho: -1/rho plus the integral from 0 of
    (sqrt(1 + s^4) - 1)/s^2, by Simpson's rule."""
    def f(s):
        return s * s / (math.sqrt(1 + s ** 4) + 1)
    h = rho / steps
    total = f(0.0) + f(rho) + sum((4 if k % 2 else 2) * f(k * h) for k in range(1, steps))
    return -1 / rho + total * h / 3


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class Published(unittest.TestCase):
    def setUp(self):
        import sympy
        import verify_metrics
        self.sp, self.vm = sympy, verify_metrics

    def read(self, system):
        """The reader of a published chart, a function from a field's name to the diagonal of its
        published components with c = 1, and the chart."""
        entry = next(c for c in published()["coordinates"] if c["id"] == system)
        reader = self.vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
        names = entry["coords"]

        def diagonal(field, variant=None):
            block = entry[field] if variant is None else entry[field]["variants"][variant]["nonzero"]
            out = [self.sp.Integer(0)] * len(names)
            for e in block:
                i, j = (names.index(k) for k in e["indices"])
                self.assertEqual(i, j, f"{system}: {field} has a component off the diagonal")
                out[i] = reader(e["value"]).subs(reader.c, 1)
            return out
        return reader, diagonal, entry

    def number(self, value, at):
        return float(self.sp.N(self.sp.sympify(value).subs(at)))

    def at_point(self, reader, radial, L=1.3):
        """Numbers for the radial coordinate, L and the angles of the sphere."""
        names = ["\\alpha", "\\beta", "\\psi", "\\theta", "\\phi"]
        at = {reader.symbol[n]: v for n, v in zip(names, (0.7, 1.1, 0.9, 1.3, 0.4))}
        at[reader.parameters["L"]] = L
        at[reader.symbol[next(n for n in reader.symbol if n not in names + ["t", "x", "y", "z"])]] = radial
        return at

    def test_every_chart_has_ten_dimensions_and_no_ricci_scalar(self):
        for system in CHARTS:
            reader, diagonal, entry = self.read(system)
            self.assertEqual(len(entry["coords"]), 10, system)
            self.assertEqual(entry["ricci_scalar"], "R = 0", system)
            self.assertEqual(len(entry["geodesics"]), 10, system)

    def test_the_mixed_ricci_tensor_is_the_stress_of_a_self_dual_five_form(self):
        # -k on t, x, y, z and rho and +k on the sphere, with k = 4 L^8/(rho^4 + L^4)^(5/2).
        reader, diagonal, entry = self.read("isotropic")
        rho, L = 0.8, 1.3
        at = self.at_point(reader, rho, L)
        k = 4 * L ** 8 / (rho ** 4 + L ** 4) ** 2.5
        mixed = [self.number(v, at) for v in diagonal("ricci_tensor", "ul")]
        for i, value in enumerate(mixed):
            self.assertAlmostEqual(value, k if i >= 5 else -k, places=10)

    def test_the_isotropic_chart_solves_duff_and_lu_s_field_equation(self):
        # R_MN = F_MPQRS F_N^PQRS/96 with F_txyz(rho) = d(1/H)/d(rho): on the brane's block
        # R_MM = F^2 times the inverse metric of the block's other four, over 4, and on the sphere the
        # same with the dual, whose component is 4 L^4 times the unit sphere's volume.
        sp = self.sp
        reader, diagonal, entry = self.read("isotropic")
        rho, L = 0.8, 1.3
        at = self.at_point(reader, rho, L)
        g = [self.number(v, at) for v in diagonal("metric_components")]
        ricci = [self.number(v, at) for v in diagonal("ricci_tensor", "ll")]
        H = 1 + L ** 4 / rho ** 4
        electric = 4 * L ** 4 / (rho ** 5 * H ** 2)
        sines = [math.sin(a) for a in (0.7, 1.1, 0.9, 1.3)]
        magnetic = 4 * L ** 4 * sines[0] ** 4 * sines[1] ** 3 * sines[2] ** 2 * sines[3]
        for i in range(10):
            block = range(5) if i < 5 else range(5, 10)
            strength = electric if i < 5 else magnetic
            others = math.prod(1 / g[a] for a in block if a != i)
            self.assertAlmostEqual(ricci[i] / (strength ** 2 * others / 4), 1.0, places=9, msg=entry["coords"][i])
        # The dual is the same form: sqrt(-g) F^txyz(rho) is the sphere's component.
        root = math.sqrt(-math.prod(g))
        self.assertAlmostEqual(root * electric * math.prod(1 / g[a] for a in range(5)) / -magnetic, 1.0, places=9)
        self.assertEqual(sp.simplify(sp.sympify(reader(entry["ricci_scalar"].partition("=")[2]))), 0)

    def test_the_areal_and_the_horizon_charts_are_the_isotropic_chart_carried_along(self):
        # rho^4 = r^4 - L^4, Gibbons, Horowitz and Townsend's (4.1), and w = rho/r.
        iso_reader, iso, _ = self.read("isotropic")
        L, rho = 1.3, 0.8
        r = (rho ** 4 + L ** 4) ** 0.25
        w = rho / r
        there = [self.number(v, self.at_point(iso_reader, rho, L)) for v in iso("metric_components")]
        for system, x, slope in (("areal", r, r ** 3 / rho ** 3), ("horizon", w, L / (1 - w ** 4) ** 1.25)):
            reader, diagonal, entry = self.read(system)
            here = [self.number(v, self.at_point(reader, x, L)) for v in diagonal("metric_components")]
            for i in range(10):
                wanted = there[i] * (slope ** 2 if i == 4 else 1)
                self.assertAlmostEqual(here[i] / wanted, 1.0, places=9, msg=f"{system} {entry['coords'][i]}")

    def test_the_horizon_chart_is_even_in_w(self):
        # Gibbons, Horowitz and Townsend's reflection (2.22): the region behind the horizon is the exterior.
        reader, diagonal, entry = self.read("horizon")
        for field in ("metric_components", "kretschmann"):
            values = diagonal(field) if field == "metric_components" else [reader(entry[field].partition("=")[2])]
            for value in values:
                a = self.number(value, self.at_point(reader, 0.6))
                b = self.number(value, self.at_point(reader, -0.6))
                self.assertAlmostEqual(a, b, places=10, msg=field)

    def test_the_inversion_of_the_radius_is_a_conformal_isometry(self):
        # Gibbons and Townsend's (7) on the isotropic chart: rho -> L^2/rho multiplies ds^2 by L^2/rho^2.
        reader, diagonal, entry = self.read("isotropic")
        L, rho = 1.3, 0.8
        here = [self.number(v, self.at_point(reader, rho, L)) for v in diagonal("metric_components")]
        there = [self.number(v, self.at_point(reader, L ** 2 / rho, L)) for v in diagonal("metric_components")]
        there[4] *= (L ** 2 / rho ** 2) ** 2
        for a, b in zip(there, here):
            self.assertAlmostEqual(a / (b * L ** 2 / rho ** 2), 1.0, places=10)

    def test_the_throat_is_the_limit_of_the_three_brane_toward_its_horizon(self):
        # rho = lambda r with t, x, y and z divided by lambda, as lambda -> 0: Maldacena's limit.
        sp = self.sp
        iso_reader, iso, _ = self.read("isotropic")
        thr_reader, thr, _ = self.read("throat")
        rho, L = iso_reader.symbol["\\rho"], iso_reader.parameters["L"]
        r, scale = sp.Symbol("r", positive=True), sp.Symbol("lambda", positive=True)
        positive = sp.Symbol("L", positive=True)
        g = [sp.sympify(v).subs(L, positive).subs(rho, scale * r) for v in iso("metric_components")]
        limits = [sp.limit(g[1] / scale ** 2, scale, 0), sp.limit(g[4] * scale ** 2, scale, 0), sp.limit(g[5], scale, 0)]
        t = [sp.sympify(v).subs({thr_reader.symbol["r"]: r, thr_reader.parameters["L"]: positive})
             for v in thr("metric_components")]
        for a, b in zip(limits, (t[1], t[4], t[5])):
            self.assertEqual(sp.simplify(a - b), 0)

    def test_the_throat_is_anti_de_sitter_space_times_a_sphere_of_one_radius(self):
        sp = self.sp
        for system in ("throat", "throat_proper"):
            reader, diagonal, entry = self.read(system)
            L = reader.parameters["L"]
            self.assertEqual(entry["weyl_tensor"]["variants"]["llll"]["nonzero"], [], system)
            mixed = diagonal("ricci_tensor", "ul")
            self.assertEqual([sp.simplify(v * L ** 2) for v in mixed], [-4] * 5 + [4] * 5, system)
            self.assertEqual(sp.simplify(reader(entry["kretschmann"].partition("=")[2]) - 80 / L ** 4), 0, system)
        # The proper distance chart is the throat at r = L e^(sigma/L).
        reader, diagonal, _ = self.read("throat")
        proper_reader, proper, _ = self.read("throat_proper")
        sigma, L = 0.4, 1.3
        here = [self.number(v, self.at_point(proper_reader, sigma, L)) for v in proper("metric_components")]
        there = [self.number(v, self.at_point(reader, L * math.exp(sigma / L), L)) for v in diagonal("metric_components")]
        there[4] *= math.exp(sigma / L) ** 2
        for a, b in zip(here, there):
            self.assertAlmostEqual(a / b, 1.0, places=10)

    def test_the_curvature_on_the_horizon_is_the_throat_s(self):
        sp = self.sp
        reader, diagonal, entry = self.read("isotropic")
        rho, L = reader.symbol["\\rho"], reader.parameters["L"]
        K = reader(entry["kretschmann"].partition("=")[2])
        self.assertEqual(sp.simplify(K.subs(rho, 0) - 80 / L ** 4), 0)
        self.assertEqual(sp.limit(K * rho ** 12, rho, sp.oo), 960 * L ** 8)


class Drawings(unittest.TestCase):
    def test_every_chart_has_a_spacetime_diagram_and_a_conformal_diagram(self):
        diagrams = drawing("diagrams")
        self.assertEqual(list(diagrams["systems"]), CHARTS)
        for system in CHARTS:
            self.assertEqual([v["id"] for v in diagrams["systems"][system]], ["radial"], system)
        self.assertEqual([v["id"] for v in drawing("conformal")["views"]], CHARTS)

    def test_the_moment_is_marked_at_the_radius_each_chart_gives_it(self):
        # The embedding reads the isotropic radius from L/20 to 3L: the areal radius is
        # (rho^4 + L^4)^(1/4) and w is rho over it; the throat's cylinder runs from sigma = -L to L,
        # which is r = L/e to eL.
        ends = {"isotropic": (0.05, 3.0), "areal": ((0.05 ** 4 + 1) ** 0.25, 82 ** 0.25),
                "horizon": (0.05 / (0.05 ** 4 + 1) ** 0.25, 3 / 82 ** 0.25),
                "throat": (math.exp(-1), math.e), "throat_proper": (-1.0, 1.0)}
        for system, (lo, hi) in ends.items():
            view = drawing("diagrams")["systems"][system][0]
            X0, X1, Y0, Y1 = view["box"]
            (mark,) = view["slices"]
            self.assertEqual(mark["view"], "throat" if system.startswith("throat") else "brane", system)
            (line,) = mark["lines"]
            self.assertAlmostEqual(X0 + line[0][0] * (X1 - X0), lo, places=3, msg=system)
            self.assertAlmostEqual(X0 + line[-1][0] * (X1 - X0), hi, places=3, msg=system)
            self.assertTrue(all(abs(Y0 + u[1] * (Y1 - Y0)) < 1e-9 for u in line), system)

    def test_the_embedded_surface_is_a_plane_far_away_and_a_cylinder_of_radius_L_down_the_throat(self):
        views = {v["id"]: v for v in drawing("embedding")["views"]}
        self.assertEqual(sorted(views), ["brane", "throat"])
        points = views["brane"]["surfaces"][0]["pieces"][0]["points"]
        for rho, radius, height in points:
            self.assertAlmostEqual(radius, (rho ** 4 + 1) ** 0.25, places=5)
        by_rho = {round(p[0], 4): p for p in points}
        # Down the throat the circles close on L and the height falls as L ln(rho).
        self.assertLess(by_rho[0.05][1] - 1, 2e-6)
        self.assertAlmostEqual(by_rho[0.1][2] - by_rho[0.05][2], math.log(2), places=3)
        # Far away the surface flattens into a plane: its slope falls as sqrt(2) L^2/rho^2.
        (r1, h1), (r2, h2), (r3, h3) = by_rho[1.0][1:], by_rho[2.0][1:], by_rho[3.0][1:]
        self.assertLess((h3 - h2) / (r3 - r2), 0.25)
        self.assertLess((h3 - h2) / (r3 - r2), 0.4 * (h2 - h1) / (r2 - r1))
        for sigma, radius, height in views["throat"]["surfaces"][0]["pieces"][0]["points"]:
            self.assertAlmostEqual(radius, 1.0, places=6)
            self.assertAlmostEqual(height, sigma, places=6)

    def test_the_conformal_diagram_is_a_tower_of_exteriors_with_no_singularity(self):
        pi = round(math.pi, 4)
        for view in drawing("conformal")["views"][:3]:
            classes = {layer["class"] for layer in view["layers"]}
            self.assertNotIn("singular", classes, view["id"])
            self.assertEqual(len([layer for layer in view["layers"] if layer["class"] == "region"]), 5, view["id"])
            cover = next(layer for layer in view["layers"] if layer["class"] == "cover")
            self.assertEqual(sorted(map(tuple, cover["points"])), sorted([(0, -pi), (-pi, 0), (0, pi), (pi, 0)]))
            # A second infinity stands on the left, at X = -2 pi, behind the horizon.
            self.assertIn([round(-2 * math.pi, 4), pi], [layer["at"] for layer in view["layers"] if layer["kind"] == "point"])
            # The moment t = 0 from rho = L/20 to 3L: T = 0 and X = 2 arctan(rho_*).
            line = view["slices"][0]["lines"][0]
            self.assertTrue(all(abs(T) < 1e-6 for _, T in line), view["id"])
            self.assertAlmostEqual(line[0][0], 2 * math.atan(tortoise(0.05)), places=3)
            self.assertAlmostEqual(line[-1][0], 2 * math.atan(tortoise(3.0)), places=3)

    def test_the_tortoise_coordinate_is_odd_and_runs_to_minus_infinity_on_the_horizon(self):
        # d(rho_*)/d(rho) = sqrt(1 + 1/rho^4), and rho_*(1/rho) + rho_*(rho) is one constant, the
        # inversion's, 2 rho_*(1).
        for rho in (0.3, 0.7, 2.0):
            h = 1e-4
            slope = (tortoise(rho + h) - tortoise(rho - h)) / (2 * h)
            self.assertAlmostEqual(slope, math.sqrt(1 + rho ** -4), places=5)
            self.assertAlmostEqual(tortoise(rho) + tortoise(1 / rho), 2 * tortoise(1.0), places=6)
        self.assertLess(tortoise(1e-3), -999)

    def test_the_throat_s_views_are_the_poincare_wedge(self):
        half, pi = round(math.pi / 2, 4), round(math.pi, 4)
        for view in drawing("conformal")["views"][3:]:
            cover = next(layer for layer in view["layers"] if layer["class"] == "cover")
            self.assertEqual(cover["points"], [[-half, 0], [half, -pi], [half, pi]], view["id"])
            self.assertEqual(view["slices"][0]["view"], "throat")


if __name__ == "__main__":
    unittest.main()
