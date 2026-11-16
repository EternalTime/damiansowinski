"""Dust collapse in anti-de Sitter space:
.venv.noindex/bin/python -m unittest _tools/test_ads_dust_collapse.py

What its charts, its History and its drawings state, held on the published files. Each cloud's published
Einstein tensor is held to dust with its cosmological constant where the scale factor obeys the equation the
page states, each exterior's published metric function to the mass the History gives and to the energy of the
edge that the junction needs, the disc's edge to the closed forms its drawings declare, the numbers its captions
state to the Friedmann equation, and the interior of Lemos's article of 1998 to being no dust. The tests need
sympy and are skipped where it is absent.
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
ID = "ads_dust_collapse"
METRIC = json.loads((DATA / "metrics" / f"{ID}.json").read_text(encoding="utf-8"))
CHARTS = ["disc_interior", "disc_exterior", "disc_exterior_ingoing", "hyperbolic_interior", "hyperbolic_slab",
          "hyperbolic_exterior", "flat_interior", "flat_exterior"]
SOURCES = {"disc_interior": "mann1993", "disc_exterior": "mann1993", "disc_exterior_ingoing": "carlip1995",
           "hyperbolic_interior": "smith1997", "hyperbolic_slab": "mena2008toroidal",
           "hyperbolic_exterior": "smith1997", "flat_interior": "mena2008toroidal",
           "flat_exterior": "mena2008toroidal"}


def chart(system):
    return next(c for c in METRIC["coordinates"] if c["id"] == system)


class Page(unittest.TestCase):
    def test_the_eight_charts_stand_in_order_and_each_cites_its_source(self):
        self.assertEqual([c["id"] for c in METRIC["coordinates"]], CHARTS)
        for system, key in SOURCES.items():
            self.assertIn(f"[{key}]", chart(system)["convention"], system)
            self.assertIn(key, METRIC["references"])

    def test_no_chart_is_stationary_or_empty_as_a_whole(self):
        for tag in ("static", "stationary", "vacuum", "Einstein space", "spherically symmetric"):
            self.assertNotIn(tag, METRIC["tags"])
        for tag in ("three-dimensional", "four-dimensional"):
            self.assertIn(tag, METRIC["tags"])


@unittest.skipUnless(HAS_SYMPY, "needs sympy")
class Dust(unittest.TestCase):
    """The published Einstein tensors of the clouds, with one index raised, at points."""

    def mixed(self, system, values, scale, velocity, pull):
        import sympy as sp
        import verify_metrics as vm
        c = chart(system)
        R = vm.Reader(c["coords"], [p["symbol"] for p in c["parameters"]], ())
        a, t = R.parameters["a"], R.symbol["t"]
        at = [(sp.Derivative(a, (t, 2)), pull), (sp.Derivative(a, t), velocity), (a, scale), (R.c, 1)]
        at += [(R.parameters[k], v) for k, v in values.items() if k in R.parameters]
        at += [(R.symbol[name], sp.Rational(3, 5 + i)) for i, name in enumerate(c["coords"])]
        out = {}
        for entry in c["einstein_tensor"]["variants"]["ul"]["nonzero"]:
            out[tuple(entry["indices"])] = complex(sp.N(R(entry["value"]).subs(at), 30))
        return c["coords"], out

    def test_the_disc_is_dust_with_the_constant_of_three_dimensions(self):
        import sympy as sp
        ell, r0 = sp.Rational(7, 5), sp.Rational(9, 7)
        for k in (sp.Rational(-3, 11), sp.Rational(2, 9)):
            for a in (sp.Rational(9, 10), sp.Rational(1, 4)):
                coords, G = self.mixed("disc_interior", {"k": k, "r_0": r0, "ell": ell}, a,
                                       -sp.sqrt(1 - a * a) / ell, -a / ell ** 2)
                density = (k + 1 / ell ** 2) / a ** 2
                self.assertEqual(set(G), {(x, x) for x in coords})
                self.assertAlmostEqual(G[("t", "t")], complex(1 / ell ** 2 - density), 20)
                for x in coords[1:]:
                    self.assertAlmostEqual(G[(x, x)], complex(1 / ell ** 2), 20)

    def test_each_cloud_is_dust_with_the_constant_of_four_dimensions(self):
        import sympy as sp
        L, b = sp.Rational(7, 5), sp.Rational(9, 4)
        excess = 1 / L ** 2 - 1 / b ** 2
        for system in ("hyperbolic_interior", "hyperbolic_slab", "flat_interior"):
            for a in (sp.Rational(9, 10), sp.Rational(1, 4)):
                speed = sp.sqrt(1 / b ** 2 + excess / a - a * a / L ** 2)
                coords, G = self.mixed(system, {"b": b, "L": L, "r_0": 3, "chi_0": 1}, a, -speed,
                                       -excess / (2 * a * a) - a / L ** 2)
                self.assertAlmostEqual(G[("t", "t")], complex(3 / L ** 2 - 3 * excess / a ** 3), 20, system)
                for x in coords[1:]:
                    self.assertAlmostEqual(G[(x, x)], complex(3 / L ** 2), 20, system)

    def test_the_interior_of_the_article_of_1998_is_no_dust(self):
        # ds^2 = -dt^2 + a^2 (dl^2 + l^2 (dtheta^2 + dphi^2)) with flat theta and phi: G^l_l - G^theta_theta = 1/(l a)^2.
        import sympy as sp
        import verify_metrics as vm
        t, l, theta, phi = sp.symbols("t l theta phi", positive=True)
        a = sp.Function("a")(t)
        geo = vm.Geometry(sp.diag(-1, a ** 2, a ** 2 * l ** 2, a ** 2 * l ** 2), [t, l, theta, phi], 10 ** 6)
        G = geo.einstein_ll()
        difference = sp.simplify(G[1][1] / a ** 2 - G[2][2] / (a ** 2 * l ** 2))
        self.assertEqual(sp.simplify(difference - 1 / (l ** 2 * a ** 2)), 0)


class Junction(unittest.TestCase):
    """The masses the History gives, as the relations among the parameters that the drawings use."""

    def test_the_disc_makes_the_mass_of_mann_and_ross_and_a_horizon_inside_its_edge(self):
        import ads_dust_collapse_drawings as d
        # M = 8 pi G_3 rho_0 r_0^2/c^2 - 1 with k = 8 pi G_3 rho_0/c^2 - 1/ell^2, at ell = 1.
        self.assertAlmostEqual(d.MASS, (d.K + 1) * d.R0 ** 2 - 1, 15)
        self.assertAlmostEqual(d.ENERGY ** 2, d.R0 ** 2 - d.MASS, 15)
        self.assertLess(math.sqrt(d.MASS), d.R0)
        # The edge R = r_0 cos t obeys (dR/dt)^2 + R^2 - M = E^2 at every moment.
        for t in (0.0, 0.4, 1.0, 1.5):
            R = d.R0 * math.cos(t)
            self.assertAlmostEqual((d.R0 * math.sin(t)) ** 2 + R * R - d.MASS, d.ENERGY ** 2, 14)
        self.assertAlmostEqual(d.R0 * math.cos(d.T_HORIZON), math.sqrt(d.MASS), 15)

    def test_the_clouds_make_the_masses_of_smith_and_mann_and_of_mena_natario_and_tod(self):
        import ads_dust_collapse_drawings as d
        for b, r0, mu, horizon, shift in ((d.HYP_B, d.HYP_R0, d.HYP_MU, d.HYP_RH, -1.0),
                                          (d.FLAT_B, d.FLAT_R0, d.FLAT_MU, d.FLAT_RH, 0.0)):
            self.assertAlmostEqual(mu, (1 - 1 / b ** 2) * r0 ** 3, 14)
            self.assertAlmostEqual(shift - mu / horizon + horizon ** 2, 0.0, 14)
            energy = r0 ** 2 / b ** 2 + shift
            for a in (0.9, 0.6, 0.3):
                R = r0 * a
                self.assertAlmostEqual(r0 ** 2 * d.friedmann(a, b) + shift - mu / R + R * R, energy, 12)
            self.assertAlmostEqual(d.friedmann(1.0, b), 0.0, 15)
            self.assertGreater(b, 1.0)          # 8 pi G rho_0/c^2 < 3/L^2
        self.assertAlmostEqual(d.HYP_R0, d.HYP_B * math.cosh(d.HYP_CHI / d.HYP_B), 14)
        self.assertAlmostEqual(d.FLAT_R0, d.FLAT_B * math.exp(d.FLAT_CHI / d.FLAT_B), 14)


class Drawings(unittest.TestCase):
    """The closed forms and the numbers the drawings and their captions state."""

    def test_the_edge_of_the_disc_in_the_static_and_the_advanced_time(self):
        import ads_dust_collapse_drawings as d
        f = lambda R: R * R - d.MASS
        step = 1e-6
        for t in (0.2, 0.6, 0.9):
            # tan t = sqrt(3) tanh(T/2), and dT/dt = E/f on the edge.
            static = lambda x: 2 * math.atanh(math.tan(x) / math.sqrt(3))
            slope = (static(t + step) - static(t - step)) / (2 * step)
            self.assertAlmostEqual(slope, d.ENERGY / f(d.R0 * math.cos(t)), 6)
        for t in (0.2, 0.9, 1.3, 1.5):
            # tan(t/2) = sqrt(3) (e^(v/2) - 1)/(3 - e^(v/2)), and dv/dt = 1/(E + r_0 sin t), regular at the horizon.
            advanced = lambda x: 2 * math.log(3 * (math.sqrt(3) * math.tan(x / 2) + 1) / (math.sqrt(3) * math.tan(x / 2) + 3))
            slope = (advanced(t + step) - advanced(t - step)) / (2 * step)
            self.assertAlmostEqual(slope, 1 / (d.ENERGY + d.R0 * math.sin(t)), 6)
            s = math.sqrt(3) * (math.exp(advanced(t) / 2) - 1) / (3 - math.exp(advanced(t) / 2))
            self.assertAlmostEqual((1 - s * s) / (1 + s * s), math.cos(t), 12)
        self.assertAlmostEqual(advanced(math.pi / 2), math.log(3), 12)
        self.assertAlmostEqual(advanced(math.pi / 3), math.log(9 / 4), 12)

    def test_the_numbers_the_captions_state(self):
        import ads_dust_collapse_drawings as d
        # The event horizon leaves the centre of the disc at eta = eta_h - chi_0.
        eta = math.atanh(math.sin(d.T_HORIZON)) - d.CHI_DISC
        self.assertAlmostEqual(2 * math.atan(math.tanh(eta / 2)), 0.27, 2)
        self.assertAlmostEqual(d.HYP_CRUNCH, 1.34, 2)
        self.assertAlmostEqual(d.HYP_CROSS, 0.68, 2)
        self.assertAlmostEqual(d.FLAT_CRUNCH, 1.42, 2)
        self.assertAlmostEqual(d.FLAT_CROSS, 1.00, 2)
        self.assertAlmostEqual(d.LEVEL, math.sqrt(5) / 2, 15)
        self.assertAlmostEqual(2 * math.pi * d.HORN, 8 * math.sqrt(2) * math.pi / 21, 14)

    def test_every_chart_has_a_spacetime_diagram(self):
        drawn = json.loads((DATA / "diagrams" / f"{ID}.json").read_text(encoding="utf-8"))["systems"]
        self.assertEqual(sorted(drawn), sorted(CHARTS))
        for system, views in drawn.items():
            for view in views:
                x0, x1, y0, y1 = view["box"]
                width = (2 * x1 if view.get("mirror") else x1 - x0)
                self.assertTrue(0.5 <= width / (y1 - y0) <= 2.0, f"{system}/{view['id']}")


if __name__ == "__main__":
    unittest.main()
