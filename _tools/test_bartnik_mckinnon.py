"""Bartnik and McKinnon's solitons, held on the numbers bartnik_mckinnon.py solves and on the
published file with the checker's own reader: python3 -m unittest discover -s _tools

No soliton is known in closed form, so the page's charts leave their functions free and every
drawing rests on a numerical solution. The solution is held here to what the literature tabulates
of it, to the field equations it was built from and to the one it was not built from, the
published G^theta_theta, to a vanishing Ricci scalar, and, in each of the three charts whose
radial coordinate is a function of the areal radius, to the same energy density on the same
sphere, which holds only if the declared functions and both of their derivatives are right.
Needs sympy, numpy and scipy, and is skipped where they are absent.
"""
import importlib.util
import json
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
HAS_ALL = all(importlib.util.find_spec(name) is not None for name in ("sympy", "numpy", "scipy", "contourpy"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))

# M. S. Volkov and D. V. Gal'tsov, Phys. Rep. 319, 1 (1999), Table 1: b_n, M_n, a_n and sigma_n.
TABLE = {1: (0.4537, 0.8286, 0.8933, 0.1264), 2: (0.6517, 0.9713, 8.8638, 0.0208), 3: (0.6970, 0.9953, 58.929, 0.0033)}
# P. Breitenlohner, P. Forgács and D. Maison, Commun. Math. Phys. 163, 141 (1994), Table 1: b_n and M_n.
PRECISE = {1: (0.45371627277, 0.82864698216), 2: (0.65172552552, 0.97134549426), 3: (0.69704005033, 0.99531647219)}


@unittest.skipUnless(HAS_ALL, "sympy, numpy, scipy and contourpy are not all installed")
class TheSeries(unittest.TestCase):
    def test_the_series_at_the_centre_and_at_infinity_solve_the_field_equations_to_their_order(self):
        import sympy as sp
        import bartnik_mckinnon as bm
        r, b, M, a = sp.symbols("r b M a", positive=True)

        def residuals(w, m):
            N = 1 - 2 * m / r
            return (sp.expand(r ** 2 * N * sp.diff(w, r, 2) + (2 * m - (1 - w ** 2) ** 2 / r) * sp.diff(w, r) + w * (1 - w ** 2)),
                    sp.expand(sp.diff(m, r) - N * sp.diff(w, r) ** 2 - (1 - w ** 2) ** 2 / (2 * r ** 2)))
        w, dw, m = bm.centre(b, r)
        self.assertEqual(sp.expand(sp.diff(w, r) - dw), 0)
        for residual, order in zip(residuals(w, m), (8, 8)):
            # The series holds w through r^6 and m through r^7, so each equation is met through r^6.
            self.assertTrue(all(sp.expand(residual).coeff(r, k) == 0 for k in range(order - 1)), residual)
        x = sp.Symbol("x", positive=True)
        w, dw, m = bm.far(M, a, r)
        self.assertEqual(sp.simplify(sp.diff(w, r) - dw), 0)
        for residual, order in zip(residuals(w, m), (5, 8)):
            # w through 1/r^5 and m through 1/r^7 meet the two equations through 1/r^5 and 1/r^8.
            series = sp.expand(residual.subs(r, 1 / x))
            self.assertTrue(all(sp.simplify(series.coeff(x, k)) == 0 for k in range(order + 1)), series.coeff(x, order))


@unittest.skipUnless(HAS_ALL, "sympy, numpy, scipy and contourpy are not all installed")
class TheSolitons(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import numpy
        import bartnik_mckinnon
        cls.np, cls.bm = numpy, bartnik_mckinnon
        cls.solitons = {n: bartnik_mckinnon.soliton(n) for n in (1, 2, 3)}

    def test_the_parameters_are_the_ones_the_literature_tabulates(self):
        for n, s in self.solitons.items():
            b, M, a, sigma = TABLE[n]
            self.assertAlmostEqual(s.b, b, delta=6e-5)
            self.assertAlmostEqual(s.M, M, delta=6e-5)
            self.assertAlmostEqual(s.a / a, 1, delta=1e-4)
            self.assertAlmostEqual(math.exp(-s.delta0), sigma, delta=1e-4)
            b, M = PRECISE[n]
            self.assertAlmostEqual(s.b, b, delta=2e-10)
            self.assertAlmostEqual(s.M, M, delta=2e-10)

    def test_the_two_halves_of_each_solution_meet(self):
        for s in self.solitons.values():
            self.assertLess(s.gap, 1e-9)

    def test_the_amplitude_has_n_zeros_stays_in_its_strip_and_ends_in_a_vacuum(self):
        np = self.np
        r = np.geomspace(1e-4, 1e7, 400001)
        for n, s in self.solitons.items():
            w = s.exact_state(r)[0]
            self.assertEqual(int(np.sum(np.sign(w[:-1]) != np.sign(w[1:]))), n)
            self.assertLessEqual(float(np.max(np.abs(w))), 1.0)
            self.assertAlmostEqual(float(w[0]), 1.0, places=6)
            self.assertAlmostEqual(float(w[-1]), (-1.0) ** n, delta=1e-4)

    def test_the_mass_and_the_lapse_rise_from_the_centre_and_no_horizon_forms(self):
        np = self.np
        r = np.geomspace(1e-4, 1e7, 400001)
        for n, s in self.solitons.items():
            j = s.exact_jet(r)
            self.assertTrue(bool(np.all(j["dm"] >= 0)))
            self.assertTrue(bool(np.all(j["ddelta"] <= 0)))
            self.assertGreater(float(np.min(j["N"])), 0)
            self.assertAlmostEqual(float(j["m"][-1]), s.M, places=9)
            self.assertAlmostEqual(float(j["delta"][-1]), 0, places=12)
        # The neck deepens toward the extreme Reissner-Nordström throat at r = ell as zeros are added.
        least = [float(np.min(self.solitons[n].exact_jet(r)["N"])) for n in (1, 2, 3)]
        self.assertAlmostEqual(least[0], 0.2424, places=4)
        self.assertTrue(least[0] > 10 * least[2] and least[1] > least[2] > 0)

    def test_the_solved_functions_keep_the_field_equations(self):
        np = self.np
        for n, s in self.solitons.items():
            r = np.geomspace(0.02, 50, 300)
            for residual in s.residual(r):
                # Central differences over r/10^4, whose own error is the square of that step.
                self.assertLess(float(np.max(np.abs(residual))), 1e-4 if n < 3 else 1e-3)

    def test_the_tables_are_the_integrations(self):
        np = self.np
        r = np.geomspace(1e-8, 1e9, 20001)
        for s in self.solitons.values():
            exact, table = s.exact_jet(r), s.jet(r)
            for key in exact:
                self.assertTrue(bool(np.all(np.abs(table[key] - exact[key]) <= 1e-9 + 1e-7 * np.abs(exact[key]))), key)

    def test_each_derivative_is_the_slope_of_the_function_before_it(self):
        np = self.np
        s = self.solitons[1]
        chains = (("m", "dm"), ("dm", "d2m"), ("delta", "ddelta"), ("ddelta", "d2delta"), ("N", "dN"), ("dN", "d2N"))
        r = np.array([0.05, 0.4, 1.0, 1.6, 2.3, 4.0, 30.0, 500.0])
        h = 1e-5 * r
        for value, slope in chains:
            difference = (s.exact_jet(r + h)[value] - s.exact_jet(r - h)[value]) / (2 * h)
            scale = np.maximum(np.abs(s.exact_jet(r)[slope]), 1e-3)
            self.assertLess(float(np.max(np.abs(difference - s.exact_jet(r)[slope]) / scale)), 1e-5, slope)
        # Inside R_NEAR every function is the series at the centre, differentiated term by term, and
        # it meets the integration there.
        edge = self.bm.R_NEAR
        inside, outside = s.exact_jet(edge * (1 - 1e-12)), s.exact_jet(edge * (1 + 1e-12))
        for key in inside:
            self.assertAlmostEqual(float(inside[key][0]), float(outside[key][0]), delta=1e-9 * max(1.0, abs(float(outside[key][0]))), msg=key)

    def test_the_two_other_radial_coordinates_invert_and_have_their_limits(self):
        np = self.np
        s = self.solitons[1]
        r = np.geomspace(1e-7, 1e8, 3001)
        self.assertLess(float(np.max(np.abs(s.from_isotropic(s.isotropic(r)) / r - 1))), 1e-12)
        self.assertLess(float(np.max(np.abs(s.from_tortoise(s.tortoise(r)) / r - 1))), 1e-12)
        # Far away the isotropic radius is Schwarzschild's, r = rho (1 + M/2 rho)^2.
        rho = float(s.isotropic(1e6)[0])
        self.assertAlmostEqual(rho * (1 + s.M / (2 * rho)) ** 2 / 1e6, 1, places=12)
        # At the centre x = r/sigma(0), and light climbs 1/sigma(0) of ct for each unit of r.
        self.assertAlmostEqual(float(s.tortoise(1e-5)[0]) / 1e-5, math.exp(s.delta0), places=8)
        # d rho/d r = rho/(r sqrt(N)) and d x/d r = 1/(sigma N), by differences.
        r = np.array([0.01, 0.3, 1.5, 4.0, 60.0])
        h = 1e-6 * r
        j = s.jet(r)
        slope = (s.isotropic(r + h) - s.isotropic(r - h)) / (2 * h)
        self.assertLess(float(np.max(np.abs(slope * r * np.sqrt(j["N"]) / s.isotropic(r) - 1))), 1e-7)
        slope = (s.tortoise(r + h) - s.tortoise(r - h)) / (2 * h)
        self.assertLess(float(np.max(np.abs(slope * np.exp(-j["delta"]) * j["N"] - 1))), 1e-7)


@unittest.skipUnless(HAS_ALL, "sympy, numpy, scipy and contourpy are not all installed")
class ThePublishedCharts(unittest.TestCase):
    """The published Einstein tensor and Ricci scalar of each chart, with the functions the
    drawings declare, on the soliton with one zero."""

    @classmethod
    def setUpClass(cls):
        import numpy
        import sympy
        import null_rays
        import verify_metrics
        cls.np, cls.sp, cls.nr, cls.vm = numpy, sympy, null_rays, verify_metrics
        cls.soliton = null_rays.bm_soliton.soliton(1)
        cls.metric = json.loads((DATA / "metrics" / "bartnik_mckinnon.json").read_text(encoding="utf-8"))

    def published(self, system, functions):
        """(the radial coordinate, {name: numpy function of it}) for the mixed Einstein tensor's
        diagonal and the Ricci scalar of one chart, with its free functions declared."""
        sp, nr = self.sp, self.nr
        entry = next(c for c in self.metric["coordinates"] if c["id"] == system)
        reader = self.vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
        x = reader.symbol[entry["coords"][1]]
        declared = {reader.parameters[name]: nr._as_lambda(reader, name, text)(x) for name, text in functions.items()}
        at = {reader.c: 1, reader.symbol["\\theta"]: sp.pi / 2}
        if "ell" in reader.parameters:
            at[reader.parameters["ell"]] = 1

        def number(text):
            value = reader(text)
            for fn, rep in declared.items():
                value = value.subs(fn, rep).doit()
            return sp.lambdify(x, value.subs(at), "numpy")
        mixed = {tuple(c["indices"]): c["value"] for c in entry["einstein_tensor"]["variants"]["ul"]["nonzero"]}
        out = {name: number(mixed[(index, index)]) for name, index in
               (("t", "t"), ("radial", entry["coords"][1]), ("theta", "\\theta"))}
        out["R"] = number(entry["ricci_scalar"].partition("=")[2])
        return out

    def stress(self, r):
        """2 T^t_t, 2 T^r_r and 2 T^theta_theta of the Yang-Mills field at areal radii r, ell = 1."""
        j = self.soliton.exact_jet(r)
        kinetic = j["N"] * j["dw"] ** 2
        V = (1 - j["w"] ** 2) ** 2 / (2 * r ** 2)
        return -2 * (kinetic + V) / r ** 2, 2 * (kinetic - V) / r ** 2, 2 * V / r ** 2

    def check(self, system, functions, coordinate_of):
        np = self.np
        r = np.array([0.05, 0.3, 0.8, 1.2, 1.6, 2.5, 5.0, 12.0, 40.0])
        x = coordinate_of(r)
        published = self.published(system, functions)
        scale = np.abs(self.stress(r)[0])
        for name, expected in zip(("t", "radial", "theta"), self.stress(r)):
            found = published[name](x)
            self.assertLess(float(np.max(np.abs(found - expected) / scale)), 1e-6, f"{system}: G^{name}_{name}")
        self.assertLess(float(np.max(np.abs(published["R"](x)) / scale)), 1e-6, f"{system}: the Ricci scalar")

    def test_the_areal_chart_carries_the_yang_mills_stress_and_no_ricci_scalar(self):
        self.check("areal", self.nr.BM_AREAL, lambda r: r)

    def test_the_isotropic_chart_carries_the_same_stress_on_the_same_spheres(self):
        self.check("isotropic", self.nr.BM_ISOTROPIC, self.soliton.isotropic)

    def test_the_tortoise_chart_carries_the_same_stress_on_the_same_spheres(self):
        self.check("tortoise", self.nr.BM_TORTOISE, self.soliton.tortoise)

    def test_the_chart_of_the_flow_carries_the_same_stress_on_the_same_spheres(self):
        self.check("flow", self.nr.BM_FLOW, lambda r: self.np.log(self.soliton.isotropic(r)))

    def test_the_drawings_mark_the_zero_of_the_amplitude(self):
        zero = self.soliton.facts()["nodes"][0]
        self.assertAlmostEqual(zero, self.nr.BM_NODE, places=4)
        self.assertAlmostEqual(float(self.soliton.w(self.nr.BM_NODE)[0]), 0, places=4)

    def test_the_slice_tests_know_where_the_embedding_ends_in_each_radial_coordinate(self):
        import test_build_mfs_data
        reach = test_build_mfs_data.BARTNIK_MCKINNON_REACH
        rho, xi = float(self.soliton.isotropic(8.0)[0]), float(self.soliton.tortoise(8.0)[0])
        self.assertAlmostEqual(reach["isotropic"], rho, places=7)
        self.assertAlmostEqual(reach["tortoise"], xi, places=7)
        self.assertAlmostEqual(reach["flow"], math.log(rho), places=7)

    def test_the_embedding_diagrams_climb_as_the_mass_function_says(self):
        np = self.np
        views = json.loads((DATA / "embedding" / "bartnik_mckinnon.json").read_text(encoding="utf-8"))["views"]
        self.assertEqual([v["id"] for v in views], ["n1", "n2", "n3"])
        for n, view in zip((1, 2, 3), views):
            piece = next(p for p in view["surfaces"][0]["pieces"] if not p.get("reference"))
            points = np.array(piece["points"], dtype=float)
            r, z = points[:, 0], points[:, -1]
            soliton = self.nr.bm_soliton.soliton(n)
            # The height is the integral of sqrt(2m/(r - 2m)) from the centre.
            grid = np.linspace(0, r[-1], 400001)
            N = soliton.jet(np.maximum(grid, 1e-12))["N"]
            slope = np.sqrt(np.maximum(1 - N, 0) / N)
            height = np.concatenate([[0], np.cumsum(0.5 * (slope[1:] + slope[:-1]) * np.diff(grid))])
            self.assertLess(float(np.max(np.abs(z - np.interp(r, grid, height)))), 2e-3 * float(z[-1]), n)


if __name__ == "__main__":
    unittest.main()
