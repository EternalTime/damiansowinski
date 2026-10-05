"""Plebański and Hacyan's products and the anti-Nariai universe:
.venv.noindex/bin/python -m unittest discover -s _tools

What its texts state, held on the published files: every chart of a product has the Einstein
tensor of two surfaces of constant curvature K1 and K2, with the cosmological constant
(K1 + K2)/2 its parameter states and a field of energy density (K2 - K1)/2 in units of
c^4/8 pi G, which is |Lambda| in Plebański and Hacyan's two and nothing in anti-Nariai; the
exceptional chart has the curvature of the product with a flat plane whatever f and g are, and
no coordinate its metric is free of; the ultracold charged black hole of de Sitter space has
its triple horizon at the radius of the sphere and the field's energy density there; the
hyperbolic black hole of anti-de Sitter space at its least mass has a double horizon whose
throat is anti-Nariai's static chart; and the drawings hold what their captions say. The tests
of the published mathematics need sympy and are skipped where it is absent; the tests of the
drawings read the files alone.
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

# Each product: its charts, the letter of its radius, and K1, K2 and Lambda in units of 1/radius^2.
PRODUCTS = {
    "a flat plane times a sphere": (("sphere", "sphere_rindler"), "b", 0, 1, (1, 2)),
    "anti-de Sitter space times a flat plane": (("plane", "plane_null", "plane_static"), "a", -1, 0, (-1, 2)),
    "anti-Nariai": (("anti_nariai", "anti_nariai_static"), "a", -1, -1, (-1, 1)),
}


def published(metric_id="plebanski_hacyan"):
    return json.loads((DATA / "metrics" / f"{metric_id}.json").read_text(encoding="utf-8"))


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class Published(unittest.TestCase):
    def setUp(self):
        import sympy
        import verify_metrics
        self.sp, self.vm = sympy, verify_metrics

    def read(self, system, metric_id="plebanski_hacyan"):
        """The reader of a published chart, its coordinates, and a function from a field's name to
        the matrix of its published components, with c = 1."""
        entry = next(c for c in published(metric_id)["coordinates"] if c["id"] == system)
        reader = self.vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
        names = entry["coords"]
        x = [reader.symbol[n] for n in names]

        def matrix(field, variant=None):
            block = entry[field] if variant is None else entry[field]["variants"][variant]["nonzero"]
            M = self.sp.zeros(len(x), len(x))
            for e in block:
                i, j = (names.index(k) for k in e["indices"])
                M[i, j] = reader(e["value"]).subs(reader.c, 1)
            return M
        return reader, x, matrix, entry

    def scalar(self, reader, text):
        return reader(text.partition("=")[2]).subs(reader.c, 1)

    def test_every_chart_has_the_einstein_tensor_of_its_two_curvatures(self):
        # G^mu_nu = -K2 on the Lorentzian surface and -K1 on the other, so that G + Lambda delta is
        # (K2 - K1)/2 times diag(-1, -1, 1, 1): a uniform field along the first surface, or none.
        sp = self.sp
        for name, (charts, letter, k1, k2, _) in PRODUCTS.items():
            for system in charts:
                reader, _, matrix, _ = self.read(system)
                radius = reader.parameters[letter]
                G = matrix("einstein_tensor", "ul")
                wanted = sp.diag(-k2, -k2, -k1, -k1) / radius ** 2
                self.assertEqual(sp.simplify(G - wanted), sp.zeros(4, 4), f"{name}, {system}")

    def test_the_cosmological_constant_and_the_field_are_what_the_parameters_state(self):
        sp = self.sp
        for name, (charts, letter, k1, k2, (num, den)) in PRODUCTS.items():
            lam = sp.Rational(num, den)
            self.assertEqual(sp.Rational(k1 + k2, 2), lam, name)
            energy = sp.Rational(k2 - k1, 2)
            self.assertEqual(energy, 0 if name == "anti-Nariai" else abs(lam), name)
            sign = "" if num > 0 else "-"
            stated = f"\\Lambda = {sign}1/{'2' if den == 2 else ''}{letter}^2"
            for system in charts:
                reader, _, _, entry = self.read(system)
                radius = reader.parameters[letter]
                self.assertIn(stated, entry["parameters"][0]["description"], f"{name}, {system}")
                # R = 4 Lambda and K = 4(K1^2 + K2^2), the scalars of a product of two surfaces.
                self.assertEqual(sp.simplify(self.scalar(reader, entry["ricci_scalar"]) - 4 * lam / radius ** 2), 0)
                self.assertEqual(sp.simplify(self.scalar(reader, entry["kretschmann"])
                                             - 4 * (k1 ** 2 + k2 ** 2) / radius ** 4), 0)

    def test_the_exceptional_chart_has_the_planes_curvature_and_no_symmetry_of_its_coordinates(self):
        sp = self.sp
        from sympy.core.function import AppliedUndef
        reader, x, matrix, entry = self.read("exceptional")
        plane_reader, _, plane, plane_entry = self.read("plane")
        same = {plane_reader.parameters["a"]: reader.parameters["a"],
                **dict(zip([plane_reader.symbol[n] for n in plane_entry["coords"]], x))}
        for field, variant in (("einstein_tensor", "ul"), ("ricci_tensor", "ul")):
            mine = matrix(field, variant)
            self.assertFalse(mine.atoms(AppliedUndef), f"{field} holds f or g")
            self.assertEqual(sp.simplify(mine - plane(field, variant).subs(same)), sp.zeros(4, 4), field)
        self.assertEqual(entry["kretschmann"], plane_entry["kretschmann"])
        self.assertEqual(entry["ricci_scalar"], plane_entry["ricci_scalar"])
        # g_uu depends on every coordinate once f and g are functions of u: no coordinate is cyclic.
        guu = matrix("metric_components")[0, 0]
        for coordinate in x:
            self.assertNotEqual(sp.diff(guu, coordinate), 0, str(coordinate))
        # With f = g = 0 it is Plebański and Hacyan's chart of the product.
        f, g = (reader.parameters[k] for k in ("f", "g"))
        flat = matrix("metric_components").subs({f: 0, g: 0}).doit()
        self.assertEqual(sp.simplify(flat - plane("metric_components").subs(same)), sp.zeros(4, 4))

    def test_the_ultracold_hole_of_de_sitter_space_has_the_spheres_radius_and_field(self):
        # 9 r_s^2/4 = 8 r_q^2 = 2/Lambda: f, f' and f'' vanish at b = 1/sqrt(2 Lambda), and the Coulomb
        # field's energy density there, r_q^2/b^4 in units of c^4/8 pi G, is Lambda.
        sp = self.sp
        reader, x, matrix, entry = self.read("static", "reissner_nordstrom_de_sitter")
        r = reader.symbol["r"]
        lam = sp.Symbol("lam", positive=True)
        at = {reader.parameters["Lambda"]: lam, reader.parameters["r_s"]: sp.sqrt(8 / (9 * lam)),
              reader.parameters["r_q"]: 1 / (2 * sp.sqrt(lam))}
        f = -matrix("metric_components")[0, 0].subs(at)
        b = 1 / sp.sqrt(2 * lam)
        for order in range(3):
            self.assertEqual(sp.simplify(sp.diff(f, r, order).subs(r, b)), 0, f"derivative {order}")
        self.assertNotEqual(sp.simplify(sp.diff(f, r, 3).subs(r, b)), 0)
        self.assertEqual(sp.simplify(at[reader.parameters["r_q"]] ** 2 / b ** 4 - lam), 0)

    def test_the_least_hyperbolic_hole_of_anti_de_sitter_space_has_anti_nariais_throat(self):
        # mu = -2L/(3 sqrt 3): a double horizon at r_h = L/sqrt(3), where f = (r - r_h)^2/r_h^2 to
        # second order, the g_tt of anti-Nariai's static chart at a = r_h far from its horizon, and
        # the hyperbolic plane there has the radius r_h, with Lambda = -3/L^2 = -1/a^2.
        sp = self.sp
        reader, x, matrix, entry = self.read("hyperbolic", "topological_black_hole")
        r, L = reader.symbol["r"], sp.Symbol("L", positive=True)
        at = {reader.parameters["L"]: L, reader.parameters["mu"]: -2 * L / (3 * sp.sqrt(3))}
        g = matrix("metric_components").subs(at)
        f = -g[0, 0]
        rh = L / sp.sqrt(3)
        self.assertEqual(sp.simplify(f.subs(r, rh)), 0)
        self.assertEqual(sp.simplify(sp.diff(f, r).subs(r, rh)), 0)
        self.assertEqual(sp.simplify(sp.diff(f, r, 2).subs(r, rh) / 2 - 1 / rh ** 2), 0)
        self.assertEqual(sp.simplify(g[2, 2].subs(r, rh) - rh ** 2), 0)
        mine, _, static, _ = self.read("anti_nariai_static")
        a, rr = mine.parameters["a"], mine.symbol["r"]
        self.assertEqual(sp.simplify(sp.diff(-static("metric_components")[0, 0], rr, 2) / 2 - 1 / a ** 2), 0)
        self.assertEqual(sp.simplify(-3 / L ** 2 + 1 / rh ** 2), 0)


class Drawings(unittest.TestCase):
    def test_every_product_chart_has_a_spacetime_diagram_and_a_conformal_diagram(self):
        charts = [system for charts, *_ in PRODUCTS.values() for system in charts]
        diagrams = json.loads((DATA / "diagrams" / "plebanski_hacyan.json").read_text(encoding="utf-8"))
        conformal = json.loads((DATA / "conformal" / "plebanski_hacyan.json").read_text(encoding="utf-8"))
        self.assertEqual(sorted(diagrams["systems"]), sorted(charts))
        self.assertEqual(sorted(v["id"] for v in conformal["views"]), sorted(charts))

    def test_the_equator_is_a_cylinder_and_the_hyperbolic_plane_a_hyperboloid(self):
        views = {v["id"]: v for v in json.loads(
            (DATA / "embedding" / "plebanski_hacyan.json").read_text(encoding="utf-8"))["views"]}
        self.assertEqual(sorted(views), ["equator", "hyperbolic_plane", "sphere"])
        tube = views["equator"]["surfaces"][0]["pieces"][0]
        for z, rho, height in tube["points"]:
            self.assertAlmostEqual(rho, 1.0, places=6)
            self.assertAlmostEqual(height, z, places=6)
        ball = views["sphere"]["surfaces"][0]["pieces"][0]
        for theta, rho, height in ball["points"]:
            self.assertAlmostEqual(rho ** 2 + (height - 1) ** 2, 1.0, places=5)
        sheet = next(p for p in views["hyperbolic_plane"]["surfaces"][0]["pieces"] if p["id"] == "sheet")
        self.assertEqual(sheet["space"], "minkowski")
        for theta, rho, height in sheet["points"]:
            # (Z + a)^2 - rho^2 = a^2, with rho = a sinh(theta), at a = 1.
            self.assertAlmostEqual((height + 1) ** 2 - rho ** 2, 1.0, places=6)
            self.assertAlmostEqual(rho, math.sinh(theta), places=6)

    def test_the_static_charts_cover_the_wedge_between_a_crossing_of_horizons_and_the_boundary(self):
        conformal = json.loads((DATA / "conformal" / "plebanski_hacyan.json").read_text(encoding="utf-8"))
        half = round(math.pi / 2, 4)
        for view in conformal["views"]:
            if view["id"] not in ("plane_static", "anti_nariai", "anti_nariai_static"):
                continue
            cover = next(layer for layer in view["layers"] if layer["kind"] == "fill" and layer["class"] == "cover")
            self.assertEqual([[round(x, 4), round(t, 4)] for x, t in cover["points"]],
                             [[0, 0], [half, -half], [half, half]], view["id"])


if __name__ == "__main__":
    unittest.main()
