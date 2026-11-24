"""The Hawking-Turok instanton: its five charts, the field equations they hold, its texts and drawings.
.venv.noindex/bin/python -m unittest discover -s _tools

Held on the published files: each chart's line element and domains to its source; the published
Einstein tensors, read through the checker's own Reader at points, to the field equations of a
minimally coupled scalar, Hawking and Turok's (field) for the instanton and its continuations,
Kirklin, Turok, and Wiseman's field and potential for the exponential instanton, and Garriga's massless
field for Vilenkin's; the published Ricci scalars to 2/(3 s^2) at a proper distance s from each
singularity, Turok and Hawking's coefficient; the History to the closed form the exponential chart
prints; the related entries to their answers; and the drawings to the closed forms they draw and to
the reasons the others are not drawn. The tests that read the mathematics need sympy and mpmath and are
skipped where they are absent.
"""
import importlib.util
import json
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
HAS_SYMPY = importlib.util.find_spec("sympy") is not None and importlib.util.find_spec("mpmath") is not None
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))
METRIC = json.loads((DATA / "metrics" / "hawking_turok_instanton.json").read_text(encoding="utf-8"))
SYSTEMS = ("instanton", "bubble_exterior", "open_universe", "exponential", "vilenkin")


def chart(system):
    return next(c for c in METRIC["coordinates"] if c["id"] == system)


class Charts(unittest.TestCase):
    def test_the_five_charts_are_the_sourced_ones_in_order(self):
        self.assertEqual([c["id"] for c in METRIC["coordinates"]], list(SYSTEMS))
        self.assertEqual(chart("instanton")["line_element"],
                         "ds^2 = d\\sigma^2 + b^2\\left(d\\psi^2 + \\sin^2\\psi\\,d\\theta^2"
                         " + \\sin^2\\psi\\sin^2\\theta\\,d\\phi^2\\right)")
        self.assertEqual(chart("bubble_exterior")["line_element"],
                         "ds^2 = d\\sigma^2 + b^2\\left(-d\\tau^2 + \\cosh^2\\tau\\,d\\theta^2"
                         " + \\cosh^2\\tau\\sin^2\\theta\\,d\\phi^2\\right)")
        self.assertEqual(chart("open_universe")["line_element"],
                         "ds^2 = -c^2dt^2 + a^2\\left(d\\chi^2 + \\sinh^2\\chi\\,d\\theta^2"
                         " + \\sinh^2\\chi\\sin^2\\theta\\,d\\phi^2\\right)")
        self.assertTrue(chart("exponential")["line_element"].startswith("ds^2 = \\ell^2\\cos\\chi\\left(d\\chi^2"))
        self.assertTrue(chart("vilenkin")["line_element"].startswith(
            "ds^2 = \\dfrac{T}{\\sqrt{A^2 + T^2}}\\left(dT^2 + \\left(A^2 + T^2\\right)"))

    def test_the_free_functions_and_the_ends_of_each_chart(self):
        self.assertEqual([p["symbol"] for p in chart("instanton")["parameters"]], ["b = b(\\sigma)", "\\sigma_f"])
        self.assertEqual([p["symbol"] for p in chart("open_universe")["parameters"]], ["a = a(t)"])
        self.assertIn("\\sigma = \\sigma_f \\;\\text{(the singularity)}", chart("instanton")["domains"])
        self.assertIn("\\chi = \\pi/2 \\;\\text{(the singularity)}", chart("exponential")["domains"])
        self.assertIn("T = 0 \\;\\text{(the singularity)}", chart("vilenkin")["domains"])
        self.assertIn("$a(t) = -ib(ict)$", chart("open_universe")["parameters"][0]["description"])
        self.assertIn("$\\ell^2 = 3c^4/(4\\pi GV_0)$", chart("exponential")["parameters"][0]["description"])

    def test_every_chart_is_conformally_flat(self):
        for system in SYSTEMS:
            self.assertEqual(chart(system)["weyl_tensor"]["variants"]["ulll"]["nonzero"], [], system)

    def test_the_scalars_of_the_closed_forms(self):
        self.assertEqual(chart("instanton")["ricci_scalar"], "R = \\dfrac{6 - 6\\left(b'\\right)^2 - 6b\\,b''}{b^2}")
        self.assertEqual(chart("exponential")["ricci_scalar"], "R = \\dfrac{45\\cos^2\\chi + 3}{2\\ell^2\\cos^3\\chi}")
        self.assertEqual(chart("vilenkin")["ricci_scalar"], "R = \\dfrac{3A^4}{2T^3\\left(A^2 + T^2\\right)^{3/2}}")
        self.assertEqual(chart("vilenkin")["kretschmann"], "K = \\dfrac{15A^8}{4T^6\\left(A^2 + T^2\\right)^3}")

    def test_the_convention_says_which_chart_is_riemannian_and_which_lorentzian(self):
        self.assertIn("The instanton and the two instantons in closed form are Riemannian", METRIC["convention"])
        self.assertIn("the region outside the bubble wall and the open universe are Lorentzian", METRIC["convention"])
        self.assertIn("with $t$ having dimensions of time", chart("open_universe")["convention"])


@unittest.skipUnless(HAS_SYMPY, "needs sympy and mpmath")
class FieldEquations(unittest.TestCase):
    """The published Einstein tensors at points, against a minimally coupled scalar field,
    G_ab = kappa (d_a phi d_b phi - g_ab ((grad phi)^2/2 + V)), kappa = 8 pi G, with G = 0.37."""

    @classmethod
    def setUpClass(cls):
        import mpmath
        import sympy as sp
        import verify_metrics as vm
        mpmath.mp.dps = 40
        cls.sp, cls.mpmath, cls.vm = sp, mpmath, vm
        cls.G = sp.Rational(37, 100)
        cls.kappa = 8 * sp.pi * cls.G

    def reader(self, system):
        c = chart(system)
        return c, self.vm.Reader(c["coords"], [p["symbol"] for p in c["parameters"]], ())

    def tensor(self, c, field, variant):
        return {tuple(m["indices"]): m["value"] for m in c[field]["variants"][variant]["nonzero"]}

    def mixed(self, system, point, carried=None):
        """The published G^a_b at the point, every value a complex number."""
        sp = self.sp
        c, R = self.reader(system)
        carried = carried or {}

        def at(text):
            return complex(sp.N(R(text).xreplace(carried).xreplace(point), 40))
        return c, R, {k: at(v) for k, v in self.tensor(c, "einstein_tensor", "ul").items()}

    def jets(self, R, function, variable, numbers):
        sp = self.sp
        f = R.parameters[function]
        x = R.symbol[variable]
        return {sp.Derivative(f, (x, 2)): numbers[2], sp.Derivative(f, x): numbers[1], f: numbers[0]}

    def test_hawking_and_turoks_field_equations_hold_on_the_instanton_and_outside_the_wall(self):
        sp, k = self.sp, self.kappa
        for system, radial in (("instanton", "\\sigma"), ("bubble_exterior", "\\sigma")):
            for b0, dphi, V in ((sp.Rational(7, 10), sp.Rational(3, 5), sp.Rational(1, 4)),
                                (sp.Rational(1, 3), sp.Rational(9, 5), sp.Rational(2, 3))):
                b1 = sp.sqrt(1 + k * b0 ** 2 * (dphi ** 2 / 2 - V) / 3)
                b2 = -k * b0 * (dphi ** 2 + V) / 3
                c, R = self.reader(system)
                point = {R.symbol[name]: sp.Rational(4 + i, 5) for i, name in enumerate(c["coords"])}
                _, _, einstein = self.mixed(system, point, self.jets(R, "b", radial, (b0, b1, b2)))
                along = complex(sp.N(k * (dphi ** 2 / 2 - V)))
                across = complex(sp.N(k * (-dphi ** 2 / 2 - V)))
                for (a, b), value in einstein.items():
                    self.assertEqual(a, b, f"{system}: G^{a}_{b}")
                    want = along if a == radial else across
                    self.assertLess(abs(value - want), 1e-12, f"{system}: G^{a}_{a} at b = {b0}")

    def test_the_open_friedmann_equations_hold_in_the_open_universe(self):
        sp, k = self.sp, self.kappa
        for a0, dphi, V in ((sp.Rational(7, 10), sp.Rational(3, 5), sp.Rational(1, 4)),
                            (sp.Rational(13, 10), sp.Rational(1, 2), sp.Rational(1, 10))):
            a1 = sp.sqrt(1 + k * a0 ** 2 * (dphi ** 2 / 2 + V) / 3)
            a2 = -k * a0 * (dphi ** 2 - V) / 3
            c, R = self.reader("open_universe")
            point = {R.symbol[name]: sp.Rational(4 + i, 5) for i, name in enumerate(c["coords"])}
            _, _, einstein = self.mixed("open_universe", point, self.jets(R, "a", "t", (a0, a1, a2)))
            for (a, b), value in einstein.items():
                want = -k * (dphi ** 2 / 2 + V) if a == "t" else k * (dphi ** 2 / 2 - V)
                self.assertLess(abs(value - complex(sp.N(want))), 1e-12, f"G^{a}_{a} at a = {a0}")

    def test_the_exponential_instanton_holds_kirklin_turok_and_wisemans_field(self):
        sp, k = self.sp, self.kappa
        c, R = self.reader("exponential")
        chi = R.symbol["\\chi"]
        B = sp.Rational(17, 10)                      # (2 pi/L) sqrt(6/kappa Lambda), which shifts the field
        lam = 6 * B / k                              # Lambda_4 at l = 1, so that l^2 = 6B/(kappa Lambda_4)
        field = -sp.sqrt(3 / (2 * k)) * sp.log(B * sp.cos(chi))
        potential = lam * sp.exp(sp.sqrt(2 * k / 3) * field)
        for x in (sp.Rational(3, 10), sp.Rational(6, 5)):
            point = {chi: x, R.symbol["\\psi"]: sp.Rational(4, 5), R.symbol["\\theta"]: sp.Rational(6, 5),
                     R.parameters["ell"]: 1}
            _, _, einstein = self.mixed("exponential", point)
            dphi2 = complex(sp.N((sp.diff(field, chi) ** 2 / sp.cos(chi)).subs(chi, x)))   # g^chichi = 1/(l^2 cos chi)
            V = complex(sp.N(potential.subs(chi, x)))
            for (a, b), value in einstein.items():
                want = k * (dphi2 / 2 - V) if a == "\\chi" else k * (-dphi2 / 2 - V)
                self.assertLess(abs(value - complex(sp.N(want))), 1e-12 * abs(complex(sp.N(want))), f"G^{a}_{a} at chi = {x}")
        # l^2 = 6/(kappa V_0): the potential at the pole is 6/kappa at l = 1.
        self.assertLess(abs(complex(sp.N(potential.subs(chi, 0) - 6 / k))), 1e-30)

    def test_vilenkins_instanton_holds_garrigas_massless_field(self):
        sp, k = self.sp, self.kappa
        c, R = self.reader("vilenkin")
        T, A = R.symbol["T"], R.parameters["A"]
        root_kappa = sp.sqrt(12 * sp.pi * self.G)
        field = -3 / (4 * root_kappa) * sp.log(T ** 2 / (A ** 2 + T ** 2))
        for x in (sp.Rational(3, 10), sp.Rational(5, 2)):
            point = {T: x, R.symbol["\\psi"]: sp.Rational(4, 5), R.symbol["\\theta"]: sp.Rational(6, 5), A: 1}
            _, _, einstein = self.mixed("vilenkin", point)
            gTT = (T / sp.sqrt(A ** 2 + T ** 2)).subs({A: 1, T: x})
            dphi2 = complex(sp.N(sp.diff(field, T).subs({A: 1, T: x}) ** 2 / gTT))
            for (a, b), value in einstein.items():
                want = complex(sp.N(k * dphi2 / 2 if a == "T" else -k * dphi2 / 2))
                self.assertLess(abs(value - want), 1e-12 * abs(want), f"G^{a}_{a} at T = {x}")

    def test_the_ricci_scalar_is_2_over_3s2_at_each_singularity(self):
        sp, mpmath = self.sp, self.mpmath
        c, R = self.reader("exponential")
        scalar = R(c["ricci_scalar"].split("=", 1)[1])
        gap = mpmath.mpf("1e-12")
        s = mpmath.quad(lambda u: mpmath.sqrt(mpmath.cos(u)), [mpmath.pi / 2 - gap, mpmath.pi / 2])
        value = sp.N(scalar.xreplace({R.symbol["\\chi"]: sp.pi / 2 - sp.Float(gap, 40), R.parameters["ell"]: 1}), 40)
        self.assertLess(abs(float(value) * float(s) ** 2 - 2 / 3), 1e-9)
        c, R = self.reader("vilenkin")
        scalar = R(c["ricci_scalar"].split("=", 1)[1])
        T = mpmath.mpf("1e-12")
        s = mpmath.quad(lambda u: mpmath.sqrt(u / mpmath.sqrt(1 + u ** 2)), [0, T])
        value = sp.N(scalar.xreplace({R.symbol["T"]: sp.Float(T, 40), R.parameters["A"]: 1}), 40)
        self.assertLess(abs(float(value) * float(s) ** 2 - 2 / 3), 1e-9)


class Texts(unittest.TestCase):
    def test_the_history_prints_the_closed_form_the_exponential_chart_holds(self):
        history = METRIC["history"]
        self.assertIn("$ds^2 \\propto \\cos\\chi\\,(d\\chi^2 + \\sin^2\\chi\\,d\\Omega_3^2)$", history)
        self.assertIn("$\\tfrac{2}{3}(\\sigma_f - \\sigma)^{-2}$", history)
        self.assertIn("$(\\sigma_f - \\sigma)^{1/3}$", history)
        self.assertIn("$a(t) = -ib(it)$", history)

    def test_every_related_spacetime_answers(self):
        kinds = {"special_case": "generalisation", "piece": "composite", "family": "family"}
        for entry in METRIC["related"]:
            other = json.loads((DATA / "metrics" / f"{entry['id']}.json").read_text(encoding="utf-8"))
            answer = [e for e in other["related"] if e["id"] == METRIC["id"]]
            self.assertEqual(len(answer), 1, entry["id"])
            self.assertEqual(answer[0]["kind"], kinds[entry["kind"]], entry["id"])


class Drawings(unittest.TestCase):
    def test_no_spacetime_or_conformal_diagram_is_written_and_each_chart_says_why(self):
        self.assertFalse((DATA / "diagrams" / "hawking_turok_instanton.json").exists())
        self.assertFalse((DATA / "conformal" / "hawking_turok_instanton.json").exists())
        import hawking_turok_instanton_drawings as drawings
        self.assertEqual(set(drawings.NOT_DRAWN), set(SYSTEMS))
        for system in drawings.EUCLIDEAN:
            self.assertIn("Positive definite", drawings.NOT_DRAWN[system])
        for system in ("bubble_exterior", "open_universe"):
            self.assertIn("left free", drawings.NOT_DRAWN[system])

    def test_the_exponential_instantons_surface_runs_from_its_pole_to_cos_chi_one_third(self):
        embedding = json.loads((DATA / "embedding" / "hawking_turok_instanton.json").read_text(encoding="utf-8"))
        view = next(v for v in embedding["views"] if v["id"] == "pole")
        piece = view["surfaces"][0]["pieces"][0]
        self.assertNotIn("space", piece)
        first, last = piece["points"][0], piece["points"][-1]
        self.assertEqual(first[:2], [0.0, 0.0])
        self.assertAlmostEqual(last[0], math.acos(1 / 3), places=12)
        for x, rho, _ in piece["points"]:
            self.assertAlmostEqual(rho, math.sqrt(math.cos(x)) * math.sin(x), places=7)
        # The surface lies level where it stops: the last chord runs inward with no rise to speak of.
        (x0, r0, z0), (x1, r1, z1) = piece["points"][-2], piece["points"][-1]
        self.assertLess(r1, r0)
        self.assertLess(abs(z1 - z0), abs(r1 - r0))

    def test_vilenkins_surface_lies_in_minkowski_space(self):
        embedding = json.loads((DATA / "embedding" / "hawking_turok_instanton.json").read_text(encoding="utf-8"))
        view = next(v for v in embedding["views"] if v["id"] == "vilenkin")
        piece = view["surfaces"][0]["pieces"][0]
        self.assertEqual(piece["space"], "minkowski")
        for T, rho, _ in piece["points"]:
            self.assertAlmostEqual(rho, math.sqrt(T) * (1 + T * T) ** 0.25, places=9)
        # Every chord is spacelike in dX^2 + dY^2 - dZ^2: the circles grow faster than the distance out.
        points = piece["points"]
        self.assertTrue(all((b[1] - a[1]) ** 2 > (b[2] - a[2]) ** 2 for a, b in zip(points, points[1:])))

    def test_the_curvature_figure_starts_at_the_values_of_the_regular_pole(self):
        curvature = json.loads((DATA / "curvature" / "hawking_turok_instanton.json").read_text(encoding="utf-8"))
        view = next(v for v in curvature["views"] if v["id"] == "exponential")
        first = {c["panel"]: c["points"][0] for c in view["curves"]}
        self.assertEqual(first["K"], [0.0, 96.0])
        self.assertEqual(first["R"], [0.0, 24.0])
        self.assertEqual([v["id"] for v in curvature["views"]], ["exponential", "vilenkin"])


if __name__ == "__main__":
    unittest.main()
