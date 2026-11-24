"""The flows of the renormalisation group of holography:
.venv.noindex/bin/python -m unittest discover -s _tools

What the published files state, held to the articles. Each chart's published line element is read
through the checker's own Reader at points and held to its source, GPPZ99's (m38) at C_1 = 0, BFS01's
Sec. 3.2 in u, or BFS01's Sec. 3.1 and FGPW00's (SUSYCases) for the disc; its published Einstein tensor
to the supersymmetric flow, G^5_5 = 6A'^2 with A' = -(2/3L)W and G^t_t - G^5_5 = 3A'' = -2 Phi'^2 with
Phi' = W'/L, BFS01's (1storder); the History to the warp factor the proper distance chart prints; and
every drawing of the disc to the one member ell = 3L/2 it is drawn at. The tests that read the
mathematics need sympy and are skipped where it is absent.
"""
import importlib.util
import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
HAS_SYMPY = importlib.util.find_spec("sympy") is not None
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))
METRIC = json.loads((DATA / "metrics" / "holographic_rg_flow.json").read_text(encoding="utf-8"))
SYSTEMS = ("proper_distance", "bfs_u", "coulomb_disc")


def chart(system):
    return next(c for c in METRIC["coordinates"] if c["id"] == system)


def drawing(kind):
    return json.loads((DATA / kind / "holographic_rg_flow.json").read_text(encoding="utf-8"))


@unittest.skipUnless(HAS_SYMPY, "needs sympy")
class Charts(unittest.TestCase):
    """Each published chart at points against its source."""

    @classmethod
    def setUpClass(cls):
        import sympy as sp
        import verify_metrics as vm
        cls.sp, cls.vm = sp, vm

    def read(self, system):
        c = chart(system)
        reader = self.vm.Reader(c["coords"], [p["symbol"] for p in c["parameters"]], ())
        g = self.vm.metric_from_line_element(reader, c["line_element"].replace("c^2", ""), c["coords"])
        return c, reader, g

    def number(self, expr, at):
        return complex(self.sp.N(self.sp.sympify(expr).subs(at), 40))

    def source(self, reader, system):
        """The warp factor e^{2A}, d(fifth)/dr, the scalar Phi and the superpotential W(Phi) of each
        chart's source, and points to hold them at."""
        sp, P = self.sp, reader.parameters
        L = P["L"]
        if system == "proper_distance":
            r = reader.symbol["r"]
            y = r / L
            m = sp.sqrt(3) / 2 * sp.log((1 + sp.exp(-y)) / (1 - sp.exp(-y)))
            W = lambda f: -sp.Rational(3, 4) * (sp.cosh(2 * f / sp.sqrt(3)) + 1)
            points = [{L: sp.Rational(3, 2), r: sp.Rational(7, 10)}, {L: 2, r: 3}]
            return r, sp.exp(y) * 2 * sp.sinh(y), sp.Integer(1), m, W, points
        if system == "bfs_u":
            u = reader.symbol["u"]
            root = sp.sqrt(1 - u)
            phi = sp.sqrt(3) / 2 * sp.log((1 + root) / (1 - root))
            W = lambda f: -sp.Rational(3, 4) * (sp.cosh(2 * f / sp.sqrt(3)) + 1)
            points = [{L: sp.Rational(3, 2), u: sp.Rational(1, 3)}, {L: 2, u: sp.Rational(19, 20)}]
            return u, u / (1 - u), 2 * (1 - u) / L, phi, W, points
        v, ell = reader.symbol["v"], P["ell"]
        W = lambda f: -sp.exp(-2 * f / sp.sqrt(6)) - sp.exp(4 * f / sp.sqrt(6)) / 2
        points = [{L: sp.Rational(3, 2), ell: 2, v: sp.Rational(1, 3)}, {L: 1, ell: sp.Rational(1, 2), v: sp.Rational(4, 5)}]
        return (v, ell ** 2 * v ** sp.Rational(2, 3) / (L ** 2 * (1 - v)), 2 * v ** sp.Rational(2, 3) * (1 - v) / L,
                sp.log(v) / sp.sqrt(6), W, points)

    def test_each_line_element_is_its_sources(self):
        for system in SYSTEMS:
            with self.subTest(system=system):
                c, reader, g = self.read(system)
                x, warp, rate, _, _, points = self.source(reader, system)
                for at in points:
                    for i in range(4):
                        self.assertLess(abs(self.number(g[i, i] - (-1 if i == 0 else 1) * warp, at)), 1e-30)
                    self.assertLess(abs(self.number(g[4, 4] - 1 / rate ** 2, at)), 1e-30)

    def test_the_published_einstein_tensor_is_the_supersymmetric_flow(self):
        """G^5_5 = 6A'^2 and G^t_t - G^5_5 = 3A'' for ds^2 = e^{2A} eta + dr^2, with BFS01's (1storder),
        A' = -(2/3L) W and Phi' = W'/L, so G^5_5 = (8/3) W^2/L^2 and G^t_t - G^5_5 = -2 W'^2/L^2."""
        sp = self.sp
        f = sp.Symbol("f")
        for system in SYSTEMS:
            with self.subTest(system=system):
                c, reader, g = self.read(system)
                x, _, _, phi, W, points = self.source(reader, system)
                L = reader.parameters["L"]
                fifth = c["coords"][-1]
                mixed = {tuple(e["indices"]): reader(e["value"])
                         for e in c["einstein_tensor"]["variants"]["ul"]["nonzero"]}
                w, dw = W(f).subs(f, phi), sp.diff(W(f), f).subs(f, phi)
                for at in points:
                    G55, Gtt = mixed[(fifth, fifth)], mixed[("t", "t")]
                    self.assertLess(abs(self.number(G55 - sp.Rational(8, 3) * w ** 2 / L ** 2, at)), 1e-25)
                    self.assertLess(abs(self.number(Gtt - G55 + 2 * dw ** 2 / L ** 2, at)), 1e-25)
                    for along in ("x_1", "x_2", "x_3"):
                        self.assertLess(abs(self.number(mixed[(along, along)] - Gtt, at)), 1e-25)


class History(unittest.TestCase):
    def test_it_states_the_warp_factor_the_proper_distance_chart_prints(self):
        history = METRIC["history"]
        self.assertIn("$e^{2A} = e^{2r/L} - 1$", history)
        self.assertIn("$ds^2 = e^{2A(r)}\\eta_{\\mu\\nu}dx^\\mu dx^\\nu + dr^2$", history)
        self.assertIn("\\left(e^{2r/L} - 1\\right)", chart("proper_distance")["line_element"])
        self.assertEqual(METRIC["references"][:3], ["maldacena1998", "girardello1998", "freedman1999"])


class Drawings(unittest.TestCase):
    """Every drawing of the disc is of the member ell = 3L/2, the one its embedding diagram draws
    and whose moment its other drawings mark, and every drawing of the confining flow is at L = 1."""

    def test_every_drawing_of_one_flow_is_of_one_member(self):
        embedding = {v["id"]: v for v in drawing("embedding")["views"]}
        self.assertEqual(sorted(embedding), ["disc", "flow"])
        self.assertIn("$\\ell = 3L/2$", embedding["disc"]["settings"])
        self.assertNotIn("\\ell", embedding["flow"]["settings"])
        conformal = {v["id"]: v for v in drawing("conformal")["views"]}
        self.assertEqual(conformal["coulomb_disc"]["settings"], "$L = 1$, $\\ell = 3L/2$.")
        self.assertEqual(conformal["proper_distance"]["settings"], "$L = 1$.")
        self.assertEqual([s["view"] for s in conformal["coulomb_disc"]["slices"]], ["disc"])
        self.assertEqual([s["view"] for s in conformal["proper_distance"]["slices"]], ["flow"])
        figure = {v["system"]: v["settings"] for v in drawing("curvature")["views"]}
        self.assertEqual(figure, {"proper_distance": "$L = 1$.", "coulomb_disc": "$L = 1$, $\\ell = 3L/2$."})
        diagrams = drawing("diagrams")["systems"]
        for system, member in (("proper_distance", "$L = 1$"), ("bfs_u", "$L = 1$"),
                               ("coulomb_disc", "$L = 1$, $\\ell = \\frac{3}{2}$")):
            views = diagrams[system] if isinstance(diagrams[system], list) else [diagrams[system]]
            self.assertEqual([v["settings"].split(", $x_1")[0] for v in views], [member])

    def test_the_disc_is_conformal_to_half_of_minkowski_space(self):
        """z = (L^2/ell) artanh(sqrt(1 - v)) runs from 0 at the boundary to infinity at the singularity,
        and dz/dv = -(L^2/2 ell)/(v sqrt(1 - v)) = -e^{-A} dr/dv at the published e^{2A} and g_vv."""
        import math
        import holographic_rg_flow_drawings as d
        for v in (0.05, 0.3, 0.7, 0.95):
            h = 1e-6
            slope = (d.disc_z(v + h) - d.disc_z(v - h)) / (2 * h)
            warp = d.ELL ** 2 * v ** (2 / 3) / (1 - v)
            dr_dv = 1 / (2 * v ** (2 / 3) * (1 - v))
            self.assertAlmostEqual(float(slope), -dr_dv / math.sqrt(warp), places=6)
        self.assertAlmostEqual(float(d.disc_z(1.0)), 0.0, places=12)
        self.assertGreater(float(d.disc_z(1e-300)), 200 * d.LAMBDA)


if __name__ == "__main__":
    unittest.main()
