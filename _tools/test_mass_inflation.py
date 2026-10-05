"""Mass inflation, the inside of a charged black hole crossed by two streams of radiation:
.venv.noindex/bin/python -m unittest discover -s _tools

What its texts state, held on the published files and on the numbers its diagrams are drawn
from. The published Einstein tensors are read back through the checker's Reader and held to the
field equations of Poisson and Israel's letter (Phys. Rev. Lett. 63, 1663, 1989): the charged
Vaidya chart to the Maxwell field of the charge plus one stream of null dust, the double null
chart to the two equations for the mass function and to its wave equation, whose source is the
product of the two fluxes, and Brady and Smith's chart to their two radial equations. Ori's
shell, as ori_shell.py integrates it, is held to the relation dm_2/f_2 = dm_1/f_1 along the
shell, to the law m_2 ~ v^(-p) exp(kappa v), to the shell's radius near the Cauchy horizon as
Bonanno, Droz, Israel and Morsink give it, and to an outgoing ray behind the shell arriving on
the Cauchy horizon at a radius above zero. The first class needs sympy and the second numpy and
scipy; each is skipped where they are absent.
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
HAS_SCIPY = importlib.util.find_spec("scipy") is not None and importlib.util.find_spec("numpy") is not None
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))


@unittest.skipUnless(HAS_SYMPY, "sympy is not installed")
class FieldEquations(unittest.TestCase):
    def setUp(self):
        import sympy
        import verify_metrics
        self.sp, self.vm = sympy, verify_metrics
        self.metric = json.loads((DATA / "metrics" / "mass_inflation.json").read_text(encoding="utf-8"))

    def published(self, system, field, variant):
        """A published tensor of one chart as a dictionary from index names to sympy values."""
        entry = next(c for c in self.metric["coordinates"] if c["id"] == system)
        reader = self.vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
        values = {tuple(e["indices"]): reader(e["value"]) for e in entry[field]["variants"][variant]["nonzero"]}
        return reader, values

    def zero(self, expression):
        self.assertEqual(self.sp.simplify(self.sp.sympify(expression).doit()), 0)

    def test_the_charged_vaidya_chart_is_the_charges_field_and_one_stream_of_null_dust(self):
        reader, G = self.published("ingoing", "einstein_tensor", "ul")
        r, v = reader.symbol["r"], reader.symbol["v"]
        m, rq = reader.parameters["m"], reader.parameters["r_q"]
        want = {("v", "v"): -rq ** 2 / r ** 4, ("r", "r"): -rq ** 2 / r ** 4, ("\\theta", "\\theta"): rq ** 2 / r ** 4,
                ("\\phi", "\\phi"): rq ** 2 / r ** 4, ("r", "v"): 2 * self.sp.diff(m, v) / r ** 2}
        self.assertEqual(set(G), set(want))
        for slot, value in want.items():
            self.zero(G[slot] - value)

    def test_the_mass_function_obeys_a_wave_equation_whose_source_is_the_product_of_the_fluxes(self):
        sp = self.sp
        reader, G = self.published("double_null", "einstein_tensor", "ll")
        u, v = reader.symbol["u"], reader.symbol["v"]
        r, sigma = reader.parameters["r"], reader.parameters["sigma"]
        q = sp.Symbol("q", positive=True)
        e = sp.exp(-2 * sigma)
        # Poisson and Israel's (1), with g^uv = -e^(-2 sigma).
        m = r / 2 * (1 + q ** 2 / r ** 2 + 2 * e * sp.diff(r, u) * sp.diff(r, v))
        Guu, Gvv, Guv, Gthth = G[("u", "u")], G[("v", "v")], G[("u", "v")], G[("\\theta", "\\theta")]
        Euv, Ethth = q ** 2 / (e * r ** 4), q ** 2 / r ** 2
        # The first of their (2): with G_vv = 2 L_in/r^2 and G_uv = E_uv, d_v m = -e^(-2 sigma) L_in d_u r.
        self.zero(sp.diff(m, v) + r ** 2 / 2 * e * (sp.diff(r, u) * Gvv - sp.diff(r, v) * (Guv - Euv)))
        self.zero(sp.diff(m, u) + r ** 2 / 2 * e * (sp.diff(r, v) * Guu - sp.diff(r, u) * (Guv - Euv)))
        # Where the two equations with no flux in them hold, d_u d_v m = (r^3/4) e^(-2 sigma) G_uu G_vv,
        # which is e^(-2 sigma) L_in L_out/r: no outflux, no inflation.
        D = sp.Derivative
        ruv, = sp.solve(sp.Eq((Guv - Euv).doit(), 0), D(r, u, v))
        suv, = sp.solve(sp.Eq((Gthth - Ethth).doit().subs(D(r, u, v), ruv), 0), D(sigma, u, v))
        on_shell = {D(r, u, v): ruv, D(sigma, u, v): suv}
        wave = (sp.diff(m, u, v) - r ** 3 / 4 * e * Guu * Gvv).doit()
        wave = wave.subs({D(r, (u, 2), v): sp.diff(ruv, u), D(r, u, (v, 2)): sp.diff(ruv, v)}).doit()
        self.zero(wave.subs(on_shell).doit().subs(on_shell))

    def test_brady_and_smiths_chart_carries_their_two_radial_equations(self):
        sp = self.sp
        reader, G = self.published("advanced", "einstein_tensor", "ll")
        _, mixed = self.published("advanced", "einstein_tensor", "ul")
        r = reader.symbol["r"]
        g, h = reader.parameters["g"], reader.parameters["h"]
        # Their (5), (ln g)_,r = (r/2) G_rr, and their (6), (r h)_,r = g (1 - r_q^2/r^2), which is
        # G^v_v + G^r_r = -2 r_q^2/r^4, the Maxwell field's, since a massless scalar field adds nothing to it.
        self.zero(G[("r", "r")] - 2 * sp.diff(g, r) / (r * g))
        self.zero(mixed[("v", "v")] + mixed[("r", "r")] - 2 * (sp.diff(r * h, r) - g) / (r ** 2 * g))


@unittest.skipUnless(HAS_SCIPY, "numpy and scipy are not installed")
class OrisShell(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import numpy
        import ori_shell
        cls.np, cls.ori, cls.S = numpy, ori_shell, ori_shell.shell()

    def test_the_hole_left_behind_is_reissner_nordstroms_own_diagram(self):
        # r_q = 0.96 m_0 is r_q = 0.48 r_s, the charge Reissner-Nordstrom's diagrams are drawn at.
        self.assertAlmostEqual(self.ori.R_PLUS, 1.28, places=12)
        self.assertAlmostEqual(self.ori.R_MINUS, 0.72, places=12)
        self.assertAlmostEqual(self.ori.KAPPA, 0.28 / 0.72 ** 2, places=12)
        self.assertAlmostEqual(float(self.ori.mass_before(1e6)), 1.0, places=12)

    def test_the_two_masses_along_the_shell_obey_dm2_over_f2_equals_dm1_over_f1(self):
        self.assertLess(self.S.consistent(), 1e-4)

    def test_before_the_tail_both_sides_are_static_and_differ_by_the_shells_mass(self):
        w = self.np.linspace(0.0, 9.9, 34)
        self.assertLess(float(self.np.max(self.np.abs(self.S.mass_behind_at(w) - 1.0))), 1e-9)
        self.assertLess(float(self.np.max(self.np.abs(self.ori.mass_before(w) - 0.98))), 1e-15)

    def test_behind_the_shell_the_mass_grows_as_oris_law(self):
        ori, S = self.ori, self.S
        # d ln m_2/dw = kappa - p/w, less a correction that falls off as 1/w^2.
        misses = []
        for w in (60.0, 70.0, 80.0):
            rate = float(S._by_v1["lnm2"](w, 1))
            misses.append((ori.KAPPA - ori.P / w) - rate)
        self.assertTrue(all(0 < b < a for a, b in zip(misses, misses[1:])), misses)
        self.assertLess(misses[-1], 1e-2)
        self.assertGreater(float(S.mass_behind_at(80.0)), 900)
        # In its own advanced time the mass passes 2 m_0 within 1e-10 m_0 of the Cauchy horizon.
        self.assertLess(float(S.mass_behind(-1e-9)), 2.0)
        self.assertGreater(float(S.mass_behind(-1e-11)), 2.0)

    def test_the_shell_closes_on_the_inner_horizon_as_the_tail_dies(self):
        # Bonanno, Droz, Israel and Morsink's (8): r_shell - r_- = dm (v_0/v)^(p-1) (1 + (p-1)/(kappa v))/(r_- kappa),
        # with corrections of order (p/(kappa v))^2, held here to twice that.
        ori, S = self.ori, self.S
        for v in (60.0, 80.0):
            tail = ori.DM * (ori.V0 / v) ** (ori.P - 1)
            want = tail * (1 + (ori.P - 1) / (ori.KAPPA * v)) / (ori.R_MINUS * ori.KAPPA)
            got = float(S.radius(v)) - ori.R_MINUS
            self.assertLess(abs(got / want - 1), 2 * (ori.P / (ori.KAPPA * v)) ** 2, f"at v = {v}")

    def test_an_outgoing_ray_behind_the_shell_arrives_on_the_cauchy_horizon_at_a_radius_above_zero(self):
        ori, S = self.ori, self.S
        start = float(S.v2_at(ori.V0))
        self.assertAlmostEqual(start, -3.7455, places=3)
        ends = [ori.outgoing_behind(start, 0.9, -10.0 ** -k) for k in (20, 40)]
        self.assertAlmostEqual(ends[1], 0.6715, places=3)
        self.assertLess(abs(ends[0] - ends[1]), 1e-6)

    def test_the_event_horizon_stands_outside_the_apparent_horizon_while_the_tail_falls_in(self):
        ori = self.ori
        horizon = float(ori.event_horizon(ori.V0))
        self.assertAlmostEqual(horizon, 1.2668, places=3)
        self.assertGreater(horizon, float(ori.horizons(ori.mass_before(ori.V0))[1]))
        self.assertAlmostEqual(float(ori.event_horizon(300.0)), ori.R_PLUS, places=9)


if __name__ == "__main__":
    unittest.main()
