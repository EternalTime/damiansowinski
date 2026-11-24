"""The black holes of constant curvature: .venv.noindex/bin/python -m unittest discover -s _tools

What the History, the charts and the drawings state, held on the published files. The first class holds
every chart to anti-de Sitter space at every point, the second the horizon and the singularity the History
describes, the third Holst and Peldan's trapped tori, and the fourth the drawings and the figure
constant_curvature_black_hole.md gives the reason for not drawing.
"""
import json
import sys
import unittest
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
ME = "constant_curvature_black_hole"
sys.path.insert(0, str(ROOT / "_tools" / "derivations"))
import verify_metrics as vm  # noqa: E402

FOUR = ("kruskal", "static", "exterior", "interior", "null")
FIVE = ("kruskal_five", "static_five", "exterior_five")


def load(folder, name=ME):
    return json.loads((DATA / folder / f"{name}.json").read_text(encoding="utf-8"))


def chart(system):
    return next(c for c in load("metrics")["coordinates"] if c["id"] == system)


def reader(system):
    c = chart(system)
    return vm.Reader(c["coords"], [p["symbol"] for p in c["parameters"]], ())


def matrix(system):
    """A published chart's metric as a sympy matrix, with its reader."""
    c = chart(system)
    rd = reader(system)
    values = {tuple(e["indices"]): e["value"] for e in c["metric_components"]}
    n = len(c["coords"])
    return rd, sp.Matrix(n, n, lambda i, k: rd(values.get((c["coords"][i], c["coords"][k]), "0")))


class Curvature(unittest.TestCase):
    def test_every_chart_has_the_scalars_of_anti_de_sitter_space(self):
        """R = -n(n - 1)/l^2 and K = 2n(n - 1)/l^4 in n dimensions, in every chart."""
        for system in FOUR + FIVE:
            n = 5 if system in FIVE else 4
            rd = reader(system)
            l = rd.parameters["l"]
            R = rd(chart(system)["ricci_scalar"].split("=", 1)[1])
            K = rd(chart(system)["kretschmann"].split("=", 1)[1])
            self.assertEqual(sp.simplify(R + n * (n - 1) / l ** 2), 0, system)
            self.assertEqual(sp.simplify(K - 2 * n * (n - 1) / l ** 4), 0, system)

    def test_every_chart_is_conformally_flat(self):
        for system in FOUR + FIVE:
            self.assertEqual(chart(system)["weyl_tensor"]["variants"]["llll"]["nonzero"], [], system)


class Horizon(unittest.TestCase):
    def test_the_norm_of_the_boost_is_r_squared_and_r_plus_squared_on_the_cone(self):
        """In Banados's Kruskal chart xi = d_phi has the squared norm r^2, which is r_+^2 on the cone
        y^2 = 0 and vanishes on the hyperboloid y^2 = -1, in four dimensions and in five."""
        for system, n in (("kruskal", 4), ("kruskal_five", 5)):
            rd, g = matrix(system)
            ys = [rd.symbol[f"y_{i}"] for i in range(n - 1)]
            rp = rd.parameters["r_+"]
            norm = g[n - 1, n - 1]
            cone = {ys[0]: sp.sqrt(sum(y ** 2 for y in ys[1:]))}
            self.assertEqual(sp.simplify(norm.subs(cone) - rp ** 2), 0, system)
            hyperboloid = {ys[0]: sp.sqrt(1 + sum(y ** 2 for y in ys[1:]))}
            self.assertEqual(sp.simplify(norm.subs(hyperboloid)), 0, system)

    def test_the_horizon_grows_along_one_cycle_at_constant_length_along_the_other(self):
        """On Holst and Peldan's event horizon u = 0 the circle of psi keeps the radius l, so its length is
        2 pi r_+, and the circle of phi has the radius l tan v, growing along the horizon."""
        rd, g = matrix("null")
        u, v = rd.symbol["u"], rd.symbol["v"]
        l = rd.parameters["l"]
        coords = chart("null")["coords"]
        psi, phi = coords.index("\\psi"), coords.index("\\varphi")
        self.assertEqual(sp.simplify(g[psi, psi].subs(u, 0) - l ** 2), 0)
        self.assertEqual(sp.simplify(g[phi, phi].subs(u, 0) - l ** 2 * sp.tan(v) ** 2), 0)

    def test_the_exterior_chart_is_a_spinless_hole_of_three_dimensions_on_every_slice_of_constant_angle(self):
        """At fixed theta the exterior chart is -(r^2 - r_+^2) dt^2 + l^2 dr^2/(r^2 - r_+^2) + r^2 dphi^2, the
        spinless hole of three dimensions at M = r_+^2/l^2 and radius l, whose c t is l t here."""
        rd, g = matrix("exterior")
        r, l, rp = rd.symbol["r"], rd.parameters["l"], rd.parameters["r_+"]
        coords = chart("exterior")["coords"]
        t, rr, ph = coords.index("t"), coords.index("r"), coords.index("\\varphi")
        self.assertEqual(sp.simplify(g[t, t] + (r ** 2 - rp ** 2)), 0)
        self.assertEqual(sp.simplify(g[rr, rr] - l ** 2 / (r ** 2 - rp ** 2)), 0)
        self.assertEqual(sp.simplify(g[ph, ph] - r ** 2), 0)


class TrappedTori(unittest.TestCase):
    def test_the_tori_of_zero_expansion_lie_on_holst_and_peldans_curve(self):
        """The tori of constant u and v have the area element sqrt(g_psipsi g_phiphi) in the published null
        chart, and their expansion along the future null normal d_v vanishes where Holst and Peldan's eq.
        (0ExpansionEquation) puts it, cosh 2t = (4r^4 - 3r^2 + 1)/((1 - 4r^2)(1 - r^2)), for every r < 1/2."""
        rd, g = matrix("null")
        u, v = rd.symbol["u"], rd.symbol["v"]
        coords = chart("null")["coords"]
        psi, phi = coords.index("\\psi"), coords.index("\\varphi")
        expansion = sp.diff(sp.log(g[psi, psi] * g[phi, phi]), v)
        for r in (sp.Rational(1, 5), sp.Rational(3, 10), sp.Rational(2, 5), sp.Rational(9, 20)):
            c2t = (4 * r ** 4 - 3 * r ** 2 + 1) / ((1 - 4 * r ** 2) * (1 - r ** 2))
            t = sp.acosh(c2t) / 2
            s = sp.sqrt((1 - r) / (1 + r))
            at = {u: sp.atan(sp.exp(-t) * s), v: sp.atan(sp.exp(t) * s), rd.parameters["l"]: 1}
            self.assertLess(abs(sp.N(expansion.subs(at), 40)), 1e-30, r)

    def test_no_torus_of_zero_expansion_reaches_past_r_one_half(self):
        """Their right side is negative beyond r = 1/2, where no cosh is, and runs to infinity as r reaches 1/2
        from below: the apparent horizon of the slices of constant r is r = 1/2."""
        for r in (0.51, 0.6, 0.8, 0.95):
            self.assertLess((4 * r ** 4 - 3 * r ** 2 + 1) / ((1 - 4 * r ** 2) * (1 - r ** 2)), 0)
        r = 0.5 - 1e-6
        self.assertGreater((4 * r ** 4 - 3 * r ** 2 + 1) / ((1 - 4 * r ** 2) * (1 - r ** 2)), 1e4)

    def test_the_apparent_horizon_is_drawn_where_the_domains_put_it(self):
        domains = chart("interior")["domains"] + chart("null")["domains"]
        self.assertIn("r = 1/2 \\;\\text{(the apparent horizon)}", domains)
        self.assertIn("\\tan u\\tan v = 1/3 \\;\\text{(the apparent horizon)}", domains)
        # tan u tan v = (1 - r)/(1 + r), which is 1/3 at r = 1/2.
        self.assertAlmostEqual((1 - 0.5) / (1 + 0.5), 1 / 3)


class Drawings(unittest.TestCase):
    def test_every_chart_draws_its_plane(self):
        drawn = {(system, v["id"]) for system, views in load("diagrams")["systems"].items() for v in views}
        planes = {system: ("plane" if system in ("kruskal", "null", "kruskal_five") else "radial")
                  for system in FOUR + FIVE}
        self.assertEqual(drawn, set(planes.items()))

    def test_the_moment_of_time_symmetry_is_marked_on_the_charts_of_four_dimensions_outside(self):
        marked = {system for system, views in load("diagrams")["systems"].items() for v in views if v.get("slices")}
        self.assertEqual(marked, {"kruskal", "static", "exterior"})
        self.assertTrue(load("conformal")["views"][0]["slices"])

    def test_every_kind_of_drawing_but_the_curvature_figure_is_published(self):
        """No curvature figure: every chart is anti-de Sitter space at every point, so both scalars are
        constants, as constant_curvature_black_hole.md says."""
        self.assertEqual([v["id"] for v in load("embedding")["views"]], ["throat"])
        self.assertEqual([v["id"] for v in load("conformal")["views"]], ["whole"])
        self.assertFalse((DATA / "curvature" / f"{ME}.json").exists())


if __name__ == "__main__":
    unittest.main()
