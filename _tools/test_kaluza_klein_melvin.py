"""The twisted space of the Kaluza-Klein Melvin universe: one flat geometry in two charts, its
reduction to a flux tube in four dimensions, and the period of its field.

Dowker, Gauntlett, Gibbons and Horowitz, hep-th/9507143, (2.1)-(2.6); Dowker, Gauntlett,
Giddings and Horowitz, hep-th/9312172, (1.5) and (5.6).
"""

import json
from pathlib import Path
import sys
import unittest

import sympy as sp

ROOT = Path(__file__).resolve().parent.parent
METRIC_FILE = ROOT / "MFS/assets/data/metrics/kaluza_klein_melvin.json"
sys.path.insert(0, str(ROOT / "_tools/derivations"))
import verify_metrics as vm  # noqa: E402

POINTS = [(sp.Rational(3, 2), sp.Rational(1, 4)), (sp.Rational(5, 3), sp.Rational(-2, 7)),
          (sp.Rational(1, 9), sp.Rational(7, 5))]


class KaluzaKleinMelvin(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        page = json.loads(METRIC_FILE.read_text(encoding="utf-8"))
        cls.charts = {c["id"]: c for c in page["coordinates"]}
        cls.metrics, cls.readers = {}, {}
        for system, chart in cls.charts.items():
            reader = vm.Reader(chart["coords"], [p["symbol"] for p in chart["parameters"]], ())
            cls.readers[system] = reader
            cls.metrics[system] = sp.Matrix(vm.metric_from_line_element(reader, chart["line_element"],
                                                                        chart["coords"]))

    def symbols(self, system):
        reader = self.readers[system]
        return [reader.symbol[c] for c in self.charts[system]["coords"]], reader.parameters["B"]

    def test_the_adapted_chart_is_the_minkowski_chart_pulled_back(self):
        # phi = varphi + B y, (2.1) with the coordinate of (2.3).
        (t, z, rho, angle, y), B = self.symbols("adapted")
        (_, _, rho_m, phi, y_m), B_m = self.symbols("minkowski")
        jacobian = sp.eye(5)
        jacobian[3, 4] = B
        pulled = jacobian.T * self.metrics["minkowski"].subs({rho_m: rho, B_m: B}) * jacobian
        for x, b in POINTS:
            difference = (pulled - self.metrics["adapted"]).subs({rho: x, B: b})
            self.assertEqual(difference, sp.zeros(5))

    def test_every_published_curvature_vanishes(self):
        for system, chart in self.charts.items():
            self.assertEqual(chart["kretschmann"], "K = 0", system)
            self.assertEqual(chart["ricci_scalar"], "R = 0", system)
            for field in ("riemann", "ricci_tensor", "einstein_tensor", "weyl_tensor"):
                for variant in chart[field]["variants"].values():
                    self.assertEqual(variant["nonzero"], [], f"{system} {field}")

    def test_the_fields_of_four_dimensions_rebuild_the_adapted_chart(self):
        # (2.4) with (2.6): e^{-4 phi/sqrt3} = Lambda, A = B rho^2/2 Lambda d varphi, and the
        # Einstein metric Lambda^{1/2}(-dt^2 + dz^2 + d rho^2) + Lambda^{-1/2} rho^2 d varphi^2, in the
        # chart x^0 = ct, so c = 1 in the published g_tt.
        (t, z, rho, angle, y), B = self.symbols("adapted")
        lam = 1 + B**2 * rho**2
        potential = B * rho**2 / (2 * lam)
        einstein = sp.diag(-sp.sqrt(lam), sp.sqrt(lam), sp.sqrt(lam), rho**2 / sp.sqrt(lam), 0)
        rebuilt = einstein / sp.sqrt(lam)
        one_form = sp.Matrix([0, 0, 0, 2 * potential, 1])
        rebuilt += lam * one_form * one_form.T
        for x, b in POINTS:
            at = {rho: x, B: b, self.readers["adapted"].c: 1}
            self.assertEqual(sp.simplify((rebuilt - self.metrics["adapted"]).subs(at)), sp.zeros(5))

    def test_the_square_of_the_field_falls_off_as_lambda_to_the_fourth(self):
        # F_{rho varphi} = d/d rho (B rho^2/2 Lambda), raised with the Einstein metric of (2.6):
        # F_{mu nu}F^{mu nu} = 2B^2/Lambda^4, which is 2B^2 on the axis, so B is the field there.
        rho, B = sp.symbols("rho B", positive=True)
        lam = 1 + B**2 * rho**2
        field = sp.diff(B * rho**2 / (2 * lam), rho)
        square = 2 * field**2 / sp.sqrt(lam) / (rho**2 / sp.sqrt(lam))
        for x, b in POINTS:
            at = {rho: x, B: abs(b)}
            self.assertEqual(sp.nsimplify(square.subs(at)), (2 * B**2 / lam**4).subs(at))
        self.assertEqual(sp.limit(square, rho, 0), 2 * B**2)

    def test_the_fifth_circle_grows_in_proportion_to_rho(self):
        # Its length is 2 pi R sqrt(g_yy) = 2 pi R Lambda^{1/2}, R on the axis and B rho R far away.
        (t, z, rho, angle, y), B = self.symbols("adapted")
        g_yy = self.metrics["adapted"][4, 4]
        self.assertEqual(sp.expand(g_yy - (1 + B**2 * rho**2)), 0)
        self.assertEqual(g_yy.subs(rho, 0), 1)
        self.assertEqual(sp.limit(sp.sqrt(g_yy) / rho, rho, sp.oo).subs(B, sp.Rational(1, 2)), sp.Rational(1, 2))

    def test_fields_differing_by_one_over_r_give_one_identification(self):
        # (2.2): (phi, y) ~ (phi + 2 pi n_1 R B + 2 pi n_2, y + 2 pi n_1 R), so B and B + 1/R
        # generate the same lattice of identifications.
        R, B = sp.symbols("R B", positive=True)
        n1, n2 = sp.symbols("n_1 n_2", integer=True)
        shifted = 2 * sp.pi * n1 * R * (B + 1 / R) + 2 * sp.pi * n2
        self.assertEqual(sp.expand(shifted - (2 * sp.pi * n1 * R * B + 2 * sp.pi * (n1 + n2))), 0)


if __name__ == "__main__":
    unittest.main()
