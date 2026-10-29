"""The black hole in a universe filled by a scalar field, Sen and Banerjee's metric of 2001 that
Cardoso and Natario found again in 2024: .venv.noindex/bin/python -m unittest discover -s _tools

What its texts and drawings state, held on the published files at m = 1 and b = 1/3, the values
every drawing of the hole takes, and b = 1 for the universe with no black hole. The first class
needs nothing but the files. The second needs numpy, scipy and sympy and is skipped where they are
absent: it holds the published spherical chart to the Einstein-scalar equations, the published
scalars to the values the captions and the History state, the tortoise coordinates the drawings
use to the published metric, and Weyl's chart to the spherical one.
"""
import importlib.util
import json
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
HAS_NUMBERS = all(importlib.util.find_spec(name) is not None for name in ("sympy", "numpy", "scipy", "contourpy"))
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))
ME = "scalar_field_black_hole"
M, B = 1.0, 1 / 3


def load(folder):
    return json.loads((DATA / folder / f"{ME}.json").read_text(encoding="utf-8"))


def chart(system):
    return next(c for c in load("metrics")["coordinates"] if c["id"] == system)


class Drawings(unittest.TestCase):
    def test_each_chart_draws_its_planes(self):
        systems = load("diagrams")["systems"]
        self.assertEqual({v["id"] for v in systems["spherical"]}, {"equator", "axis"})
        self.assertEqual({v["id"] for v in systems["weyl"]}, {"equator", "axis"})
        self.assertEqual({v["id"] for v in systems["scalar_universe"]}, {"radial"})
        conformal = {v["id"] for v in load("conformal")["views"]}
        self.assertEqual(conformal, {"spherical_equator", "spherical_axis", "weyl_equator", "weyl_axis",
                                     "scalar_universe"})

    def test_the_equator_and_the_universe_end_on_a_timelike_singularity(self):
        """On the equator and in the universe the tortoise coordinate is finite at infinity, so the
        far edge is a singular curve; on the axis the metric is Schwarzschild's and the far edges
        are null infinity."""
        for view in load("conformal")["views"]:
            classes = {layer.get("class") for layer in view["layers"]}
            self.assertIn("singular", classes, view["id"])
            self.assertEqual("scri" in classes, view["id"].endswith("axis"), view["id"])
            if view["id"] != "scalar_universe":
                self.assertIn("horizon", classes, view["id"])

    def test_each_curvature_view_says_what_it_shows_in_at_most_200_characters(self):
        for view in load("curvature")["views"]:
            text = (" ".join(view["caption"]) + view["clip"] + " ".join(item[2] for item in view["legend"])
                    + view["settings"])
            self.assertLessEqual(len(text), 200, view["id"])
            self.assertTrue(view["label"] and view["settings"] and view["clip"], view["id"])
        scalars = {view["map"]["scalar"]: view for view in load("curvature")["views"]}
        self.assertEqual(set(scalars), {"kretschmann", "ricci_scalar"})
        self.assertTrue(scalars["kretschmann"]["caption"][0].startswith("The Kretschmann scalar $K$ of"))
        self.assertTrue(scalars["ricci_scalar"]["caption"][0].startswith("The Ricci scalar $R$ of"))

    def test_the_history_credits_sen_and_banerjee_first(self):
        """The metric is the third set of Sen and Banerjee, whose arXiv v4 and journal issue are of
        2001, and Cardoso and Natario found it again in 2024: both are cited, in that order."""
        entry = load("metrics")
        self.assertLessEqual({"sen2001", "cardoso2025", "herdeiro2025"}, set(entry["references"]))
        history = entry["history"]
        self.assertLess(history.index("Somasri Sen"), history.index("Vitor Cardoso"))
        self.assertIn("2001", entry["description"])


@unittest.skipUnless(HAS_NUMBERS, "needs sympy, numpy and scipy")
class Physics(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import sympy as sp

        import chart_printer as cp
        import verify_metrics as vm
        cls.sp, cls.cp, cls.vm = sp, cp, vm

    def published(self, system):
        """The published chart's geometry, built from its line element with c = 1."""
        c = chart(system)
        line = c["line_element"].replace("c^2", "")
        params = [p["symbol"] for p in c["parameters"]]
        return self.cp.Chart(c["coords"], params, line)

    def test_the_spherical_chart_solves_the_einstein_scalar_equations(self):
        """R_ab = 2 d_a phi d_b phi with phi = b(r - m) cos(theta), and at b = 0 Schwarzschild's
        metric with r_s = 2m."""
        sp, vm = self.sp, self.vm
        ch = self.published("spherical")
        t, r, th, ph = ch.symbols
        m, b = ch.reader.parameters["m"], ch.reader.parameters["b"]
        field = b * (r - m) * sp.cos(th)
        ricci = ch.geo.ricci_ll()
        for i in range(4):
            for j in range(i, 4):
                self.assertEqual(vm.norm(vm._at(ricci, (i, j)) - 2 * sp.diff(field, ch.symbols[i])
                                         * sp.diff(field, ch.symbols[j])), 0, (i, j))
        self.assertEqual(vm.norm(ch.geo.g[1, 1].subs(b, 0) - 1 / (1 - 2 * m / r)), 0)

    def test_the_published_scalars_are_the_ones_the_texts_quote(self):
        """R = 2b^2(r^2 - 2mr + m^2 sin^2 theta)/(r^2 F); K = 48m^2/r^6 at b = 0; on the horizon
        K = (3/4)(1 - 4 beta/3 + beta^2)/m^4, 0.648/m^4 on the equator at b = 1/3; on the axis
        K tends to 12 b^4; the universe's K is 12 b^4 e^{2b^2 rho^2}."""
        sp = self.sp
        c = chart("spherical")
        read = self.vm.Reader(c["coords"], [p["symbol"] for p in c["parameters"]], ())
        r, th = read.symbol["r"], read.symbol["\\theta"]
        m, b = read.parameters["m"], read.parameters["b"]
        K = read(c["kretschmann"].split("=", 1)[1])
        R = read(c["ricci_scalar"].split("=", 1)[1])
        F = sp.exp(-b ** 2 * (r ** 2 - 2 * m * r) * sp.sin(th) ** 2)
        self.assertEqual(sp.simplify(R - 2 * b ** 2 * (r ** 2 - 2 * m * r + m ** 2 * sp.sin(th) ** 2) / (r ** 2 * F)), 0)
        self.assertEqual(sp.simplify(K.subs(b, 0) - 48 * m ** 2 / r ** 6), 0)
        at = {m: 1, b: sp.Rational(1, 3)}
        equator = float(K.subs(at).subs({r: 2, th: sp.pi / 2}))
        self.assertAlmostEqual(equator, 0.75 * (1 - 4 / 27 + 1 / 81), places=12)
        self.assertEqual(f"{equator:.3f}", "0.648")
        far = [float(K.subs(at).subs({r: x, th: 0})) for x in (1e3, 1e5)]
        self.assertAlmostEqual(far[1] / (12 / 81), 1, places=4)
        u = chart("scalar_universe")
        ru = self.vm.Reader(u["coords"], [p["symbol"] for p in u["parameters"]], ())
        rho, bu = ru.symbol["\\rho"], ru.parameters["b"]
        self.assertEqual(sp.simplify(ru(u["kretschmann"].split("=", 1)[1]) - 12 * bu ** 4 * sp.exp(2 * bu ** 2 * rho ** 2)), 0)

    def test_light_reaches_infinity_in_the_times_the_captions_quote(self):
        """On the equator a ray from r = 3m reaches r = infinity within 3.85 m/c, and in the universe
        a ray from the axis within sqrt(pi/2)/(bc); the drawings' tortoise coordinates agree with
        the published metric's dr_*/dr = sqrt(g_rr/-g_tt)."""
        import numpy as np
        from scipy.integrate import quad

        import null_rays as nr
        slope = lambda s: math.exp(-B ** 2 * (s * s - 2 * s) / 2) * s / (s - 2)  # noqa: E731
        self.assertEqual(f"{quad(slope, 3, np.inf)[0]:.2f}", "3.85")
        for x in (0.5, 1.5, 2.6, 4.0, 7.0):
            d = 1e-6
            numeric = (nr._sfbh_equator(np.array([x + d]))[0] - nr._sfbh_equator(np.array([x - d]))[0]) / (2 * d)
            self.assertAlmostEqual(numeric / slope(x), 1, places=6)
        for p in (0.3, 1.0, 2.5):
            x = math.sqrt(1 + p * p)
            rate = math.exp(-B ** 2 * p * p / 2) * (x + 1) ** 2 / (x * p)
            d = 1e-6
            numeric = (nr._sfbh_weyl_equator(np.array([p + d]))[0] - nr._sfbh_weyl_equator(np.array([p - d]))[0]) / (2 * d)
            self.assertAlmostEqual(numeric / rate, 1, places=6)
        self.assertAlmostEqual(float(nr._sfbh_universe(np.array([50.0]))[0]), math.sqrt(math.pi / 2), places=12)

    def test_weyls_chart_is_the_spherical_chart_carried(self):
        """With rho = sqrt(r^2 - 2mr) sin(theta) and z = (r - m) cos(theta), Weyl's g_tt, conformal
        factor and g_phiphi are the spherical chart's, at points in forty digits."""
        sp = self.sp
        w = chart("weyl")
        read = self.vm.Reader(w["coords"], [p["symbol"] for p in w["parameters"]], ())
        rho, z = read.symbol["\\rho"], read.symbol["z"]
        m, b = read.parameters["m"], read.parameters["b"]
        psi, gam = read.parameters["psi"], read.parameters["gamma"]
        for rv, tv in ((3, sp.Rational(1, 2)), (sp.Rational(9, 4), 1), (6, sp.Rational(3, 2))):
            at = {m: 1, b: sp.Rational(1, 3), rho: sp.sqrt(rv ** 2 - 2 * rv) * sp.sin(tv), z: (rv - 1) * sp.cos(tv)}
            f = 1 - sp.Rational(2) / rv
            F = sp.exp(-sp.Rational(1, 9) * (rv ** 2 - 2 * rv) * sp.sin(tv) ** 2)
            self.assertLess(abs(sp.N((sp.exp(2 * psi) - f).subs(at), 40)), 1e-30)
            # d rho^2 + dz^2 = (r^2 - 2mr + m^2 sin^2 theta)(dr^2/(r^2 - 2mr) + dtheta^2)
            conformal = sp.exp(2 * gam - 2 * psi) * (rv ** 2 - 2 * rv + sp.sin(tv) ** 2) - rv ** 2 * F
            self.assertLess(abs(sp.N(conformal.subs(at), 40)), 1e-30)


if __name__ == "__main__":
    unittest.main()
