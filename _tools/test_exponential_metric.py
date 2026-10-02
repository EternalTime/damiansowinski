"""The exponential metric of Papapetrou and Yilmaz, g_tt = -1/g_rr = -e^(-2m/r):
python3 -m unittest discover -s _tools

What its texts and drawings state, held on the published files. The first class needs nothing
but the files: the throat where every spacetime diagram marks it, the circles and the level
circle of the embedding diagram, the diamond of the conformal diagrams with its singular left
edges, and the relations each side answers. The second reads the published tensors through the
checker's Reader and holds them to the scalar field of negative energy, to the throat, to the
orbits the History quotes and to the affine parameter that reaches r = 0; it needs sympy and is
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

CHARTS = ["isotropic", "cartesian", "areal", "harmonic"]


def load(folder, metric_id="exponential_metric"):
    return json.loads((DATA / folder / f"{metric_id}.json").read_text(encoding="utf-8"))


def areal(r):
    """The areal radius of the sphere of isotropic radius r, at m = 1."""
    return r * math.exp(1 / r)


class Drawings(unittest.TestCase):
    def marked(self, view, kind):
        X0, X1 = view["box"][:2]
        return sorted(X0 + line[0][0] * (X1 - X0) for m in view["markers"] if m["kind"] == kind for line in m["lines"])

    def test_every_spacetime_diagram_marks_the_throat_where_the_spheres_are_smallest(self):
        """r = m, x = +-m, R = e m and u = 1/m, at m = 1."""
        systems = load("diagrams")["systems"]
        self.assertEqual(sorted(systems), sorted(CHARTS))
        want = {"isotropic": [1.0], "cartesian": [-1.0, 1.0], "areal": [math.e], "harmonic": [1.0]}
        for system, throats in want.items():
            (view,) = systems[system]
            got = self.marked(view, "throat")
            self.assertEqual(len(got), len(throats), system)
            for a, b in zip(got, throats):
                self.assertAlmostEqual(a, b, delta=2e-3, msg=system)

    def test_the_drawn_windows_are_square(self):
        for system, (view,) in load("diagrams")["systems"].items():
            X0, X1, Y0, Y1 = view["box"]
            self.assertAlmostEqual((X1 - X0) / (Y1 - Y0), 1.0, places=9, msg=system)

    def test_the_faint_lines_are_the_spheres_of_one_size_on_either_side_of_the_throat(self):
        systems = load("diagrams")["systems"]
        for system, to_r in (("isotropic", lambda x: x), ("harmonic", lambda x: 1 / x)):
            (view,) = systems[system]
            X0, X1 = view["box"][:2]
            for level in view["areal"]:
                at = [to_r(X0 + line[0][0] * (X1 - X0)) for line in level["lines"]]
                self.assertEqual(len(at), 2, system)
                self.assertLess(min(at), 1.0)
                self.assertGreater(max(at), 1.0)
                for r in at:
                    self.assertAlmostEqual(areal(r), level["level"], delta=2e-2, msg=system)

    def test_the_moment_drawn_on_each_plane_is_the_embedding_diagrams_reach(self):
        """t = 0 from r = m/3 to 6m: the whole of it on the isotropic and harmonic planes, both
        halves on the Cartesian line, and the near side from the throat out in the areal radius."""
        systems = load("diagrams")["systems"]
        want = {"isotropic": [(1 / 3, 6.0)], "cartesian": [(-6.0, -1 / 3), (1 / 3, 6.0)],
                "areal": [(math.e, areal(6.0))], "harmonic": [(1 / 6, 3.0)]}
        for system, spans in want.items():
            (view,) = systems[system]
            X0, X1, Y0, Y1 = view["box"]
            (mark,) = view["slices"]
            self.assertEqual(len(mark["lines"]), len(spans), system)
            for line, (lo, hi) in zip(mark["lines"], spans):
                for point in line:
                    self.assertAlmostEqual(Y0 + point[1] * (Y1 - Y0), 0.0, places=6, msg=system)
                xs = sorted(X0 + point[0] * (X1 - X0) for point in line)
                self.assertAlmostEqual(xs[0], max(lo, X0), delta=2e-3, msg=system)
                self.assertAlmostEqual(xs[-1], min(hi, X1), delta=2e-3, msg=system)

    def test_the_embedding_has_the_circles_of_the_metric_and_lies_level_on_half_the_throats_radius(self):
        """rho = r e^(m/r): least on r = m, where it is e m, and the two parts, in Minkowski space
        inside r = m/2 and in flat space outside, meet level on the circle of radius e^2 m/2."""
        (view,) = load("embedding")["views"]
        pieces = {p["id"]: p for p in view["surfaces"][0]["pieces"]}
        self.assertEqual(sorted(pieces), ["deep", "outer"])
        self.assertEqual(pieces["deep"].get("space"), "minkowski")
        self.assertIsNone(pieces["outer"].get("space"))
        for piece in pieces.values():
            for r, rho, _ in piece["points"]:
                self.assertAlmostEqual(rho, areal(r), places=6)
        outer, deep = pieces["outer"]["points"], pieces["deep"]["points"]
        self.assertAlmostEqual(min(p[1] for p in outer), math.e, places=6)
        self.assertEqual(outer[0][0], 0.5)
        self.assertEqual(deep[-1][0], 0.5)
        for point in (outer[0], deep[-1]):
            self.assertAlmostEqual(point[1], math.e ** 2 / 2, places=6)
            self.assertAlmostEqual(point[2], 0.0, places=6)
        # Level there: the first chord of each part rises far less than it runs.
        self.assertLess(abs(outer[1][2] - outer[0][2]), 0.2 * abs(outer[1][1] - outer[0][1]))
        self.assertLess(abs(deep[-2][2] - deep[-1][2]), 0.2 * abs(deep[-2][1] - deep[-1][1]))
        # The height rises with r throughout, through the throat, where rho turns back.
        heights = [p[2] for p in deep] + [p[2] for p in outer[1:]]
        self.assertEqual(heights, sorted(heights))

    def test_the_embedding_in_minkowski_space_is_spacelike_and_nears_the_light_cone(self):
        (view,) = load("embedding")["views"]
        deep = next(p for p in view["surfaces"][0]["pieces"] if p["id"] == "deep")["points"]
        slopes = []
        for a, b in zip(deep, deep[1:]):
            drho, dz = abs(b[1] - a[1]), abs(b[2] - a[2])
            self.assertGreater(drho, dz)
            slopes.append(dz / drho)
        # dZ/d rho = sqrt((m/r)(m/r - 2))/(m/r - 1), which is sqrt(3)/2 at r = m/3 and falls to 0 at m/2.
        self.assertAlmostEqual(slopes[0], math.sqrt(3) / 2, delta=0.02)
        self.assertLess(slopes[-1], 0.2)

    def test_the_conformal_diagrams_are_a_diamond_with_infinity_on_the_right_and_a_singular_left(self):
        views = {v["id"]: v for v in load("conformal")["views"]}
        self.assertEqual(sorted(views), ["areal", "harmonic", "isotropic"])
        for vid, view in views.items():
            lines = {}
            for layer in view["layers"]:
                if layer["kind"] in ("line", "zig"):
                    lines.setdefault(layer["class"], []).append(layer)
            for cls, side in (("scri", 1), ("singular", -1)):
                self.assertEqual(len(lines[cls]), 2, vid)
                for layer in lines[cls]:
                    (x0, t0), (x1, t1) = layer["points"][0], layer["points"][-1]
                    self.assertAlmostEqual(abs(x1 - x0), abs(t1 - t0), places=3, msg=vid)
                    self.assertTrue(all(side * p[0] >= -1e-9 for p in layer["points"]), vid)
            self.assertTrue(all(layer["kind"] == "zig" for layer in lines["singular"]), vid)
            (throat,) = lines["throat"]
            self.assertTrue(all(abs(p[0]) < 1e-9 for p in throat["points"]), vid)
            cover = next(layer for layer in view["layers"] if layer["kind"] == "fill" and layer["class"] == "cover")
            leftmost = min(p[0] for p in cover["points"])
            # The areal radius covers the near side alone, the right half.
            self.assertAlmostEqual(leftmost, 0.0 if vid == "areal" else -math.pi, places=3, msg=vid)

    def test_the_moment_on_the_conformal_diagrams_crosses_the_throat(self):
        for view in load("conformal")["views"]:
            (mark,) = view["slices"]
            xs = [p[0] for line in mark["lines"] for p in line]
            self.assertLess(min(xs), -1.0, view["id"])
            self.assertGreater(max(xs), 0.5, view["id"])
            for line in mark["lines"]:
                for p in line:
                    self.assertAlmostEqual(p[1], 0.0, places=3)


class Relations(unittest.TestCase):
    def test_each_relation_is_answered_with_the_kind_that_goes_with_it(self):
        answers = {"family": "family", "generalisation": "special_case", "limit_source": "limit"}
        own = load("metrics")["related"]
        self.assertEqual([r["id"] for r in own],
                         ["schwarzschild", "ellis_bronnikov", "fisher_jnw", "morris_thorne", "ppn_metric"])
        for relation in own:
            other = load("metrics", relation["id"])
            (back,) = [r for r in other["related"] if r["id"] == "exponential_metric"]
            self.assertEqual(back["kind"], answers[relation["kind"]], relation["id"])

    def test_both_lines_of_the_limit_say_so(self):
        own = next(r for r in load("metrics")["related"] if r["id"] == "fisher_jnw")
        back = next(r for r in load("metrics", "fisher_jnw")["related"] if r["id"] == "exponential_metric")
        for text in (own["text"], back["text"]):
            self.assertIn("\\gamma \\to \\infty", text)
            self.assertIn("limit", text)


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class Published(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import sympy as sp
        import verify_metrics as vm
        cls.sp, cls.vm = sp, vm
        cls.metric = load("metrics")
        cls.charts = {c["id"]: c for c in cls.metric["coordinates"]}

    def reader(self, system):
        chart = self.charts[system]
        held = self.vm.HELD.get(("exponential_metric", system), ())
        return self.vm.Reader(chart["coords"], [p["symbol"] for p in chart["parameters"]], (), held=held)

    def components(self, system, field, variant=None):
        reader = self.reader(system)
        block = self.charts[system][field]
        entries = block["variants"][variant]["nonzero"] if variant else block
        return reader, {tuple(e["indices"]): reader(e["value"]) for e in entries}

    def test_the_charts_are_the_four_the_literature_uses(self):
        self.assertEqual(list(self.charts), CHARTS)

    def test_the_isotropic_ricci_tensor_is_that_of_a_scalar_field_of_negative_energy(self):
        """R_ab = -2 d_a phi d_b phi with phi = m/r, the one component R_rr = -2m^2/r^4, which
        grows without bound toward r = 0 while the Kretschmann scalar goes to zero."""
        sp = self.sp
        reader, ricci = self.components("isotropic", "ricci_tensor", "ll")
        r, m = reader.symbol["r"], reader.parameters["m"]
        self.assertEqual(list(ricci), [("r", "r")])
        field = m / r
        self.assertEqual(sp.simplify(ricci[("r", "r")] + 2 * sp.diff(field, r) ** 2), 0)
        K = reader(self.charts["isotropic"]["kretschmann"].partition("=")[2])
        self.assertEqual(sp.limit(K.subs(m, 1), r, 0, "+"), 0)
        self.assertEqual(sp.limit(ricci[("r", "r")].subs(m, 1), r, 0, "+"), -sp.oo)

    def test_the_energy_density_is_negative_everywhere(self):
        """G^t_t = -8 pi G rho/c^4 is positive, m^2 e^(-2m/r)/r^4."""
        sp = self.sp
        reader, einstein = self.components("isotropic", "einstein_tensor", "ul")
        r, m = reader.symbol["r"], reader.parameters["m"]
        self.assertEqual(sp.simplify(einstein[("t", "t")] - m ** 2 * sp.exp(-2 * m / r) / r ** 4), 0)

    def test_the_radius_is_an_affine_parameter_along_a_radial_ray(self):
        """g_tt g_rr = -1, so e^(-2m/r) d(ct)/d lambda is constant and dr/d lambda is too."""
        sp = self.sp
        reader, g = self.components("isotropic", "metric_components")
        self.assertEqual(sp.simplify(g[("t", "t")] * g[("r", "r")]), -1)

    def test_the_throat_is_at_r_equal_m_and_has_the_areal_radius_e_m(self):
        sp = self.sp
        reader, g = self.components("isotropic", "metric_components")
        r, m = reader.symbol["r"], reader.parameters["m"]
        area = g[("\\theta", "\\theta")]
        self.assertEqual(sp.simplify(sp.diff(area, r).subs(r, m)), 0)
        self.assertGreater(sp.simplify(sp.diff(area, r, 2).subs(r, m)), 0)
        self.assertEqual(sp.simplify(area.subs(r, m) - sp.E ** 2 * m ** 2), 0)

    def test_light_circles_at_2m_and_the_last_stable_orbit_is_at_three_plus_root_five_m(self):
        """The effective potentials of Boonserm, Ngampitipan, Simpson and Visser's (7.7) and (7.8),
        from the published metric: V_0 = -g_tt L^2/g_phiphi peaks at r = 2m, and the angular
        momentum of a circular orbit of a particle has its least value at r = (3 + sqrt 5) m."""
        sp = self.sp
        reader, g = self.components("isotropic", "metric_components")
        r, m = reader.symbol["r"], reader.parameters["m"]
        theta = reader.symbol["\\theta"]
        lapse = -g[("t", "t")]
        around = g[("\\phi", "\\phi")].subs(theta, sp.pi / 2)
        light = lapse / around
        self.assertEqual(sp.simplify(sp.diff(light, r).subs(r, 2 * m)), 0)
        L = sp.Symbol("L", positive=True)
        V = lapse * (1 + L ** 2 / around)
        L2 = sp.solve(sp.diff(V, r), L ** 2)[0]
        isco = (3 + sp.sqrt(5)) * m
        self.assertEqual(sp.simplify(sp.diff(L2, r).subs(r, isco)), 0)
        self.assertEqual(sp.simplify(L2 - r ** 2 * m * sp.exp(2 * m / r) / (r - 2 * m)), 0)

    def test_far_away_the_metric_is_the_post_newtonian_one_at_general_relativitys_values(self):
        """-g_tt = 1 - 2m/r + 2 beta m^2/r^2 and g_rr = 1 + 2 gamma m/r with beta = gamma = 1."""
        sp = self.sp
        reader, g = self.components("isotropic", "metric_components")
        r, m = reader.symbol["r"], reader.parameters["m"]
        e = sp.Symbol("e", positive=True)
        lapse = sp.series((-g[("t", "t")]).subs(m, e * r), e, 0, 3).removeO()
        space = sp.series(g[("r", "r")].subs(m, e * r), e, 0, 2).removeO()
        self.assertEqual(sp.expand(lapse - (1 - 2 * e + 2 * e ** 2)), 0)
        self.assertEqual(sp.expand(space - (1 + 2 * e)), 0)

    def test_it_is_the_limit_of_fisher_janis_newman_and_winicours_metric(self):
        """(1 - 2m/(gamma r))^gamma -> e^(-2m/r) as gamma -> infinity, in every component of the
        spherical chart of fisher_jnw with b = 2m/gamma."""
        sp = self.sp
        other = load("metrics", "fisher_jnw")
        chart = next(c for c in other["coordinates"] if c["id"] == "spherical")
        theirs = self.vm.Reader(chart["coords"], [p["symbol"] for p in chart["parameters"]], ())
        reader, mine = self.components("isotropic", "metric_components")
        r, m = reader.symbol["r"], reader.parameters["m"]
        gamma = theirs.parameters["gamma"]
        names = {theirs.symbol["r"]: r, theirs.symbol["\\theta"]: reader.symbol["\\theta"],
                 theirs.parameters["b"]: 2 * m / gamma}
        for entry in chart["metric_components"]:
            value = theirs(entry["value"]).subs(names)
            value = value.subs({r: 3, m: 1, reader.symbol["\\theta"]: 1})
            want = mine[tuple(entry["indices"])].subs({r: 3, m: 1, reader.symbol["\\theta"]: 1})
            self.assertAlmostEqual(float(sp.limit(value, gamma, sp.oo)), float(want), places=12,
                                   msg=str(entry["indices"]))

    def test_the_areal_charts_radius_solves_its_definition(self):
        """r = -m/W(-m/R) has r e^(m/r) = R, on the near side at R = 4m and on the far side
        on the branch W_-1."""
        import mpmath
        for branch, side in ((0, "near"), (-1, "far")):
            r = -1 / mpmath.lambertw(-mpmath.mpf(1) / 4, branch)
            self.assertAlmostEqual(float(r * mpmath.exp(1 / r)), 4.0, places=12, msg=side)
            self.assertTrue(r > 1 if side == "near" else r < 1)

    def test_the_harmonic_coordinate_is_the_reciprocal_of_the_isotropic_radius(self):
        sp = self.sp
        reader, g = self.components("harmonic", "metric_components")
        iso_reader, iso = self.components("isotropic", "metric_components")
        u, m = reader.symbol["u"], reader.parameters["m"]
        r, m0 = iso_reader.symbol["r"], iso_reader.parameters["m"]
        to = {r: 1 / u, m0: m, iso_reader.symbol["\\theta"]: reader.symbol["\\theta"]}
        self.assertEqual(sp.simplify(g[("t", "t")] - iso[("t", "t")].subs(to)), 0)
        self.assertEqual(sp.simplify(g[("u", "u")] - iso[("r", "r")].subs(to) / u ** 4), 0)
        self.assertEqual(sp.simplify(g[("\\theta", "\\theta")] - iso[("\\theta", "\\theta")].subs(to)), 0)


if __name__ == "__main__":
    unittest.main()
