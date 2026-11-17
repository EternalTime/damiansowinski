"""Chan, Horne and Mann's charged dilaton black holes:
.venv.noindex/bin/python -m unittest _tools/test_dilaton_unusual_asymptotics.py

What the page's History, conventions, related texts and captions state, held on the published files. The
first chart is held to the metric the History prints, with gamma^2 = r_0^{2/(1+a^2)}, and to the History's
horizon, temperature and fields; light reaches r -> infinity at a finite t for a^2 < 1 alone; the
curvature falls to zero far away; the hole has one horizon whatever its charge; the limits the related
texts name are the published charts' limits; the horizons the captions of the Liouville families state are
the zeros of the published g^rr at the drawn values; and each chart's own check fails on a metric that is
not the solution. The field equations themselves are checked once, by print_charts.py as each chart is
printed, and the tensors against the line element by verify_metrics.py, so neither is repeated here. The
tests need sympy and are skipped where it is absent.
"""
import importlib.util
import json
import math
import re
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
HAS_SYMPY = importlib.util.find_spec("sympy") is not None
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))
ID = "dilaton_unusual_asymptotics"
METRIC = json.loads((DATA / "metrics" / f"{ID}.json").read_text(encoding="utf-8"))
COUPLINGS = ("1/3", "1/2", "1", "7/5", "3")


def chart(system):
    return next(c for c in METRIC["coordinates"] if c["id"] == system)


def history():
    return METRIC["history"].replace("¶", " ")


@unittest.skipUnless(HAS_SYMPY, "needs sympy")
class Published(unittest.TestCase):
    """The published charts read back, with every defined name written out."""

    @classmethod
    def setUpClass(cls):
        import sympy as sp
        import verify_metrics as vm
        cls.sp, cls.vm = sp, vm
        cls.readers, cls.metrics, cls.inverses = {}, {}, {}
        for c in METRIC["coordinates"]:
            R = vm.Reader(c["coords"], [p["symbol"] for p in c["parameters"]], ())
            cls.readers[c["id"]] = R
            cls.metrics[c["id"]] = cls.matrix(R, c, "metric_components")
            cls.inverses[c["id"]] = cls.matrix(R, c, "inverse_metric_components")

    @classmethod
    def matrix(cls, R, c, field):
        out = cls.sp.zeros(4, 4)
        for entry in c[field]:
            i, j = (c["coords"].index(name) for name in entry["indices"])
            out[i, j] = R(entry["value"])
        return out

    def symbols(self, system):
        R = self.readers[system]
        return [R.symbol[name] for name in chart(system)["coords"]], R.parameters

    def scalar(self, system, field):
        return self.readers[system](chart(system)[field].split("=", 1)[1])

    def same(self, a, b, at):
        self.assertLess(abs(complex(self.sp.N((a - b).subs(at), 30))), 1e-20, (a, b, at))


class FirstChart(Published):

    def test_the_history_prints_the_first_chart_with_gamma_squared_r_0_to_the_2_over_1_plus_a_squared(self):
        sp = self.sp
        (t, r, th, ph), P = self.symbols("chan_horne_mann")
        a, r0, rh = P["a"], P["r_0"], P["r_h"]
        self.assertIn("U = \\frac{r^{2/(1+a^2)}}{\\gamma^2}\\left(1 - \\frac{r_h}{r}\\right)", history())
        self.assertIn("R = \\gamma r^N", history())
        self.assertIn("N = a^2/(1 + a^2)", history())
        gamma = r0 ** (1 / (1 + a ** 2))
        U = r ** (2 / (1 + a ** 2)) / gamma ** 2 * (1 - rh / r)
        R2 = (gamma * r ** (a ** 2 / (1 + a ** 2))) ** 2
        g = self.metrics["chan_horne_mann"]
        for value in COUPLINGS:
            at = {a: sp.Rational(value), r0: sp.Rational(5, 4), rh: sp.Rational(2, 3), r: sp.Rational(9, 4),
                  th: sp.Rational(1, 3)}
            self.same(g[0, 0], -U, at)
            self.same(g[1, 1], 1 / U, at)
            self.same(g[2, 2], R2, at)

    def test_the_horizon_is_r_h_and_the_mass_is_a_squared_r_h_over_2_1_plus_a_squared(self):
        sp = self.sp
        (t, r, th, ph), P = self.symbols("chan_horne_mann")
        self.assertIn("r_h = \\frac{2(1+a^2)}{a^2}\\,M", history())
        grr = self.inverses["chan_horne_mann"][1, 1]
        for value in COUPLINGS:
            at = {P["a"]: sp.Rational(value), P["r_0"]: 2, P["r_h"]: sp.Rational(3, 2)}
            self.same(grr.subs(r, P["r_h"]), 0, at)
            # One horizon for every charge: g^rr keeps its sign on either side of r_h, in 0 < r.
            for x in (sp.Rational(1, 100), sp.Rational(1, 2), sp.Rational(149, 100)):
                self.assertLess(float(grr.subs({**at, r: x})), 0)
            for x in (sp.Rational(151, 100), 3, 100):
                self.assertGreater(float(grr.subs({**at, r: x})), 0)
        self.assertIn("the hole has no inner horizon", history())

    def test_no_value_of_the_metric_holds_the_charge(self):
        free = set().union(*(self.sp.sympify(v).free_symbols for v in self.metrics["chan_horne_mann"]))
        self.assertNotIn("Q", {s.name for s in free})
        self.assertIn("its charge may take any value but zero", history())

    def test_the_temperature_is_the_historys_and_at_a_equal_one_the_same_for_every_mass(self):
        sp = self.sp
        (t, r, th, ph), P = self.symbols("chan_horne_mann")
        a, r0, rh = P["a"], P["r_0"], P["r_h"]
        self.assertIn("T = r_h^{(1-a^2)/(1+a^2)}/(4\\pi\\gamma^2)", history())
        U = -self.metrics["chan_horne_mann"][0, 0]
        T = sp.diff(U, r).subs(r, rh) / (4 * sp.pi)
        gamma2 = r0 ** (2 / (1 + a ** 2))
        for value in COUPLINGS:
            at = {a: sp.Rational(value), r0: sp.Rational(7, 4), rh: sp.Rational(5, 3)}
            self.same(T, rh ** ((1 - a ** 2) / (1 + a ** 2)) / (4 * sp.pi * gamma2), at)
        at_one = T.subs(a, 1)
        self.assertEqual(sp.simplify(sp.diff(at_one, rh)), 0)

    def test_light_reaches_infinity_at_a_finite_t_for_a_squared_below_one_alone(self):
        sp = self.sp
        (t, r, th, ph), P = self.symbols("chan_horne_mann")
        U = -self.metrics["chan_horne_mann"][0, 0]
        for value, finite in (("1/3", True), ("1/2", True), ("1", False), ("2", False)):
            # dt/dr = 1/U along a radial light ray; far away 1/U falls as r^{-2/(1+a^2)}.
            at = {P["a"]: sp.Rational(value), P["r_0"]: 1, P["r_h"]: 1}
            power = sp.limit(sp.log(U.subs(at)) / sp.log(r), r, sp.oo)
            self.assertEqual(power, 2 / (1 + sp.Rational(value) ** 2))
            self.assertEqual(bool(power > 1), finite, value)
        self.assertIn("for $a^2 < 1$ a light ray reaches infinity in a finite time $t$", history())

    def test_the_curvature_falls_to_zero_far_away_and_the_dilaton_grows_as_log_r(self):
        sp = self.sp
        (t, r, th, ph), P = self.symbols("chan_horne_mann")
        for field in ("ricci_scalar", "kretschmann"):
            value = self.scalar("chan_horne_mann", field)
            for c in COUPLINGS:
                at = {P["a"]: sp.Rational(c), P["r_0"]: 1, P["r_h"]: 1}
                self.assertEqual(sp.limit(value.subs(at), r, sp.oo), 0, (field, c))
        # The convention's e^{2 a phi} = r^2 (r_0/r)^{2/(1+a^2)}/((1 + a^2) Q^2) is R^2/((1 + a^2) Q^2), the
        # History's, so phi = (a/(1 + a^2)) ln r plus a constant.
        a, r0, Q = P["a"], P["r_0"], sp.Symbol("Q", positive=True)
        convention = r ** 2 * (r0 / r) ** (2 / (1 + a ** 2)) / ((1 + a ** 2) * Q ** 2)
        R2 = self.metrics["chan_horne_mann"][2, 2]
        for c in COUPLINGS:
            at = {a: sp.Rational(c), r0: sp.Rational(4, 3), Q: sp.Rational(2, 7), r: sp.Rational(11, 5)}
            self.same(convention, R2 / ((1 + a ** 2) * Q ** 2), at)
            phi = sp.log(convention) / (2 * a)
            self.same(sp.diff(phi, r), a / ((1 + a ** 2) * r), at)
        self.assertIn("e^{2a\\phi} = r^2(r_0/r)^{2/(1 + a^2)}/((1 + a^2)Q^2)", chart("chan_horne_mann")["convention"])


class Limits(Published):
    """The limits and members the related texts and the History name, on the published charts."""

    def test_at_a_to_infinity_the_first_chart_is_schwarzschilds(self):
        sp = self.sp
        (t, r, th, ph), P = self.symbols("chan_horne_mann")
        g = self.metrics["chan_horne_mann"]
        at = {P["r_0"]: sp.Rational(3, 2), P["r_h"]: 2}
        for i, expected in ((0, -(1 - 2 / r)), (1, 1 / (1 - 2 / r)), (2, r ** 2)):
            self.assertEqual(sp.simplify(sp.limit(g[i, i].subs(at), P["a"], sp.oo) - expected), 0)

    def test_at_a_equal_zero_the_first_chart_is_bertotti_and_robinsons_product_of_radius_r_0(self):
        sp = self.sp
        (t, r, th, ph), P = self.symbols("chan_horne_mann")
        at = {P["a"]: 0, P["r_0"]: sp.Rational(3, 2), P["r_h"]: 2}
        self.assertEqual(sp.simplify(self.metrics["chan_horne_mann"][2, 2].subs(at) - sp.Rational(9, 4)), 0)
        self.assertEqual(sp.simplify(self.scalar("chan_horne_mann", "ricci_scalar").subs(at)), 0)
        self.assertEqual(sp.simplify(self.scalar("chan_horne_mann", "kretschmann").subs(at) - 8 / sp.Rational(3, 2) ** 4), 0)

    def test_at_a_equal_zero_the_family_b_equal_one_over_a_is_reissner_and_nordstroms(self):
        sp = self.sp
        (t, r, th, ph), P = self.symbols("liouville_inverse")
        g = self.metrics["liouville_inverse"]
        at = {P["a"]: 0, P["r_0"]: 1, P["M"]: sp.Rational(3, 2), P["q"]: 1}
        self.assertEqual(sp.simplify(g[0, 0].subs(at) + (1 - 3 / r + 1 / r ** 2)), 0)
        self.assertEqual(sp.simplify(g[2, 2].subs(at) - r ** 2), 0)

    def test_the_family_b_equal_minus_a_at_sigma_one_and_minus_one_over_a_at_lambda_zero_are_the_first_chart(self):
        sp = self.sp
        (t, r, th, ph), P = self.symbols("chan_horne_mann")
        for system, name in (("liouville_minus_a", "sigma"), ("liouville_minus_inverse", "lambda")):
            (t2, r2, th2, ph2), Q = self.symbols(system)
            for c in COUPLINGS[:3]:
                at = {Q["a"]: sp.Rational(c), Q["r_0"]: 2, Q["r_h"]: 1, Q[name]: 1 if name == "sigma" else 0}
                first = {P["a"]: sp.Rational(c), P["r_0"]: 2, P["r_h"]: 1}
                for i in range(4):
                    self.same(self.metrics[system][i, i].subs({r2: r, th2: th}),
                              self.metrics["chan_horne_mann"][i, i], {**at, **first, r: sp.Rational(7, 3),
                                                                      th: sp.Rational(1, 2)})

    def test_the_massless_chart_is_flat_on_its_plane_with_spheres_of_area_pi_xi_squared_minus_eta_squared(self):
        sp = self.sp
        (eta, xi, th, ph), _ = self.symbols("massless")
        g = self.metrics["massless"]
        self.assertEqual(g[:2, :2], sp.diag(-1, 1))
        self.assertEqual(sp.simplify(4 * sp.pi * g[2, 2] - sp.pi * (xi ** 2 - eta ** 2)), 0)


class Captions(Published):
    """The horizons the captions of the Liouville families state are the zeros of the published g^rr."""

    def zeros(self, system, params):
        sp = self.sp
        (t, r, th, ph), P = self.symbols(system)
        grr = sp.lambdify(r, self.inverses[system][1, 1].subs({P[k]: sp.sympify(v) for k, v in params.items()}))
        xs = [0.05 + 0.001 * i for i in range(8000)]
        values = [float(grr(x)) for x in xs]
        return [round(x, 2) for x, y, z in zip(xs, values, values[1:]) if (y < 0) != (z < 0)]

    def test_b_equal_one_over_a_has_its_horizons_at_r_0_over_2_and_3_r_0_over_2(self):
        import dilaton_unusual_asymptotics_drawings as d
        self.assertEqual(self.zeros("liouville_inverse", d.INVERSE_PARAMS), [0.5, 1.5])

    def test_b_equal_minus_one_over_a_has_the_event_and_cosmological_horizons_its_caption_states(self):
        import dilaton_unusual_asymptotics_drawings as d
        found = self.zeros("liouville_minus_inverse", d.MINUS_INVERSE_PARAMS)
        self.assertEqual(found, [round(d.EVENT, 2), round(d.COSMOLOGICAL, 2)])
        caption = json.loads((DATA / "diagrams" / f"{ID}.json").read_text(encoding="utf-8"))
        text = next(v for v in caption["systems"]["liouville_minus_inverse"])["caption"][0]
        self.assertEqual([float(x) for x in re.findall(r"r = (\d\.\d\d)\\,r_h", text)], found)


@unittest.skipUnless(HAS_SYMPY, "needs sympy")
class Checks(unittest.TestCase):
    """Each chart's own check fails on a metric that is not the solution."""

    def test_each_charts_check_fails_on_a_metric_that_is_not_the_solution(self):
        import chart_printer as cp
        import dilaton_unusual_asymptotics as d
        broken = {"chan_horne_mann": ("\\dfrac{r^2}{h}", "\\dfrac{2r^2}{h}"),
                  "areal": ("\\dfrac{4r^2\\,dr^2}", "\\dfrac{3r^2\\,dr^2}"),
                  "string_magnetic": ("+ 2Q^2\\left(d", "+ 3Q^2\\left(d"),
                  "liouville_minus_inverse": ("\\dfrac{r^2}{h}\\left(d", "\\dfrac{2r^2}{h}\\left(d"),
                  "massless": ("\\dfrac{\\xi^2 - \\eta^2}{4}", "\\dfrac{\\xi^2 - \\eta^2}{3}")}
        for system, (old, new) in broken.items():
            spec = d.spec(system)
            line = spec["chart_line_element"]
            self.assertIn(old, line, system)
            chart_ = cp.Chart(spec["system"]["coords"], spec["system"]["parameters"], line.replace(old, new),
                              spec["printer"], None, None, None, spec.get("reduce"), None,
                              d.HELD.get(system, ()), False)
            with self.assertRaises(AssertionError, msg=system):
                spec["check"](chart_)


if __name__ == "__main__":
    unittest.main()
