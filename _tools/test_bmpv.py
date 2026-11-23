"""The spinning black hole of Breckenridge, Myers, Peet and Vafa: .venv.noindex/bin/python -m unittest discover -s _tools

What the History, the charts and the drawings state, held on the published files. The first class holds
the curvature the History and the captions quote, the second the horizon the History describes, the
third the drawings and the plane bmpv.md gives the reason for not drawing.
"""
import json
import sys
import unittest
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
ME = "bmpv"
sys.path.insert(0, str(ROOT / "_tools" / "derivations"))
import verify_metrics as vm  # noqa: E402


def load(folder, name=ME):
    return json.loads((DATA / folder / f"{name}.json").read_text(encoding="utf-8"))


def chart(system):
    return next(c for c in load("metrics")["coordinates"] if c["id"] == system)


def reader(system):
    c = chart(system)
    return vm.Reader(c["coords"], [p["symbol"] for p in c["parameters"]], ())


def scalar(system, field):
    """A published scalar of a chart, read as sympy, with its radius and the two parameters."""
    rd = reader(system)
    text = chart(system)[field].split("=", 1)[1]
    return rd(text), rd.symbol["r"], rd.parameters["mu"], rd.parameters["j"]


def matrix(system):
    """A published chart's metric as a sympy matrix, with its reader."""
    c = chart(system)
    rd = reader(system)
    values = {tuple(e["indices"]): e["value"] for e in c["metric_components"]}
    n = len(c["coords"])
    return rd, sp.Matrix(n, n, lambda i, k: rd(values.get((c["coords"][i], c["coords"][k]), "0")))


EVENTS = [(sp.Rational(13, 10), sp.Rational(7, 10), sp.Rational(17, 10)),
          (sp.Rational(1, 1), sp.Rational(9, 5), sp.Rational(5, 2)),
          (sp.Rational(1, 1), sp.Rational(1, 2), sp.Rational(31, 10))]


class Curvature(unittest.TestCase):
    def test_the_ricci_scalar_is_gauntlett_myers_and_townsends(self):
        """R = 2(2 sin^2 beta - lambda^2)/lambda^8 in units mu = 1, their eq. (singu), which is
        2(2j^2 - mu^2 r^2)/r^8 in the radius of the Hopf, Euler and ingoing charts and the same with
        r^2 + mu for r^2 in the isotropic one."""
        for system in ("bmpv", "euler", "ingoing", "isotropic"):
            R, r, mu, j = scalar(system, "ricci_scalar")
            rr = r ** 2 + mu if system == "isotropic" else r ** 2
            want = 2 * (2 * j ** 2 - mu ** 2 * rr) / rr ** 4
            for m0, j0, r0 in EVENTS:
                self.assertEqual(sp.nsimplify(sp.N((R - want).subs({mu: m0, j: j0, r: r0}), 40), [], tolerance=1e-30), 0,
                                 system)

    def test_the_kretschmann_scalar_is_one_function_of_the_radius(self):
        K0, r0_, mu0, j0 = scalar("euler", "kretschmann")
        for system in ("bmpv", "ingoing", "isotropic"):
            K, r, mu, j = scalar(system, "kretschmann")
            for m, jj, x in EVENTS:
                here = sp.sqrt(x ** 2 - m) if system == "isotropic" else x
                a = sp.N(K.subs({mu: m, j: jj, r: here}), 40)
                b = sp.N(K0.subs({mu0: m, j0: jj, r0_: x}), 40)
                self.assertLess(abs(a - b), 1e-25, system)

    def test_the_curvature_figure_quotes_the_horizon_values(self):
        K, r, mu, j = scalar("ingoing", "kretschmann")
        R, *_ = scalar("ingoing", "ricci_scalar")
        at = {mu: 1, j: sp.Rational(1, 2), r: 1}
        self.assertEqual(sp.nsimplify(K.subs(at)), sp.Rational(101, 2))
        self.assertEqual(sp.nsimplify(R.subs(at)), -1)


class Horizon(unittest.TestCase):
    def test_the_horizon_area_is_two_pi_squared_root_of_mu_cubed_less_j_squared(self):
        """The section v constant, r = sqrt(mu) of the ingoing chart: its angular block has the root of its
        determinant sqrt(r^6 - j^2) sin(theta)/8, and over theta, phi and psitilde the area
        2 pi^2 sqrt(mu^3 - j^2), as Gauntlett, Myers and Townsend find."""
        rd, g = matrix("ingoing")
        r, th = rd.symbol["r"], rd.symbol["\\theta"]
        mu, j = rd.parameters["mu"], rd.parameters["j"]
        block = g[2:, 2:]
        for m0, j0 in ((sp.Rational(13, 10), sp.Rational(7, 10)), (1, sp.Rational(1, 2))):
            at = {mu: m0, j: j0, r: sp.sqrt(m0), th: sp.Rational(2, 5)}
            density = sp.sqrt(sp.N(block.subs(at).det(), 40))
            area = density / sp.sin(sp.Rational(2, 5)) * 2 * (2 * sp.pi) * (4 * sp.pi)
            self.assertLess(abs(sp.N(area - 2 * sp.pi ** 2 * sp.sqrt(m0 ** 3 - j0 ** 2), 40)), 1e-25)

    def test_the_horizon_does_not_rotate(self):
        """On r = sqrt(mu) of the ingoing chart the Killing vector of v is null and orthogonal to the Hopf
        fibre, so it is the generator of the horizon with no turning added."""
        rd, g = matrix("ingoing")
        r, mu = rd.symbol["r"], rd.parameters["mu"]
        coords = chart("ingoing")["coords"]
        v, psi = coords.index("v"), coords.index("\\tilde\\psi")
        self.assertEqual(sp.simplify(g[v, v].subs(r, sp.sqrt(mu))), 0)
        self.assertEqual(sp.simplify(g[v, psi].subs(r, sp.sqrt(mu))), 0)

    def test_the_hopf_fibre_is_timelike_inside_the_velocity_of_light_surface(self):
        rd, g = matrix("euler")
        r, mu, j = rd.symbol["r"], rd.parameters["mu"], rd.parameters["j"]
        gamma = chart("euler")["coords"].index("\\gamma")
        fibre = sp.simplify(g[gamma, gamma] - (r ** 6 - j ** 2) / (4 * r ** 4))
        self.assertEqual(fibre, 0)


class Drawings(unittest.TestCase):
    def test_each_chart_but_the_hopf_one_draws_its_plane(self):
        drawn = {(system, v["id"]) for system, views in load("diagrams")["systems"].items() for v in views}
        self.assertEqual(drawn, {("euler", "radial"), ("isotropic", "radial"), ("ingoing", "radial")})

    def test_every_kind_of_drawing_is_published(self):
        self.assertEqual([v["id"] for v in load("embedding")["views"]], ["fibre"])
        self.assertEqual([v["id"] for v in load("conformal")["views"]], ["exterior"])
        self.assertEqual([v["id"] for v in load("curvature")["views"]], ["radial"])


if __name__ == "__main__":
    unittest.main()
