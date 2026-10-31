"""The non-Abelian magnetic source and the closed radiation solution.

Gibbons and Steif (2.1)-(2.12), specialized to Hosotani's f=1/2.
The gauge curvature is derived from Pauli matrices and the Maurer-Cartan equation.
The source tests are separate from the general tensor verification.
"""
import json
import math
import sys
import unittest
from pathlib import Path

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "_tools" / "derivations"))
import verify_metrics as vm


def charts():
    return json.loads((ROOT / "MFS/assets/data/metrics/hosotani_universe.json").read_text())["coordinates"]


class MagneticSource(unittest.TestCase):
    def test_curvature_includes_the_su2_commutator(self):
        generators = [sp.Matrix([[0, 1], [1, 0]]) / 2,
                      sp.Matrix([[0, -sp.I], [sp.I, 0]]) / 2,
                      sp.diag(1, -1) / 2]
        f, e, a = sp.symbols("f e a", positive=True)
        connection = [-f * T / e for T in generators]
        for i, j in ((0, 1), (1, 2), (2, 0)):
            # d omega^a = epsilon_abc omega^b wedge omega^c / 2.
            derivative = sum((-f / e * sp.LeviCivita(k, i, j) * generators[k]
                              for k in range(3)), sp.zeros(2))
            commutator = -sp.I * e * (connection[i] * connection[j] - connection[j] * connection[i])
            curvature = derivative + commutator
            expected = sum((f * (f - 1) / e * sp.LeviCivita(k, i, j) * generators[k]
                            for k in range(3)), sp.zeros(2))
            self.assertEqual((curvature - expected).applyfunc(sp.expand), sp.zeros(2))
            self.assertNotEqual(curvature.subs(f, sp.Rational(1, 2)), derivative.subs(f, sp.Rational(1, 2)))
            # theta^i=a omega^i/2 supplies the factor 4 in an orthonormal frame.
            frame = 4 * curvature.subs(f, sp.Rational(1, 2)) / a**2
            self.assertEqual(frame, -generators[(i + 2) % 3] / (e * a**2))
            self.assertEqual(curvature.subs(f, 0), sp.zeros(2))
            self.assertEqual(curvature.subs(f, 1), sp.zeros(2))

    def test_magnetic_stress_is_isotropic_and_fixes_the_radius(self):
        for radius, coupling, grav in ((0.3, 0.7, 1.1), (1.2, 1.3, 0.8), (2.7, 0.4, 1.6)):
            magnetic = -np.eye(3) / (coupling * radius**2)
            rho = np.sum(magnetic**2) / 2
            pressure = np.eye(3) * rho - magnetic.T @ magnetic
            self.assertAlmostEqual(rho, 3 / (2 * coupling**2 * radius**4))
            np.testing.assert_allclose(pressure, np.eye(3) * rho / 3, atol=1e-12)
            ell2 = 4 * math.pi * grav / coupling**2
            self.assertAlmostEqual(8 * math.pi * grav * rho / 3, ell2 / radius**4)

    def test_fixed_gauge_equilibrium_and_its_homogeneous_instability(self):
        f = sp.Symbol("f")
        potential = 2 * (f**2 - f)**2
        self.assertEqual(potential.subs(f, sp.Rational(1, 2)), sp.Rational(1, 8))
        self.assertEqual(sp.diff(potential, f).subs(f, sp.Rational(1, 2)), 0)
        self.assertEqual(-sp.diff(potential, f, 2).subs(f, sp.Rational(1, 2)), 2)


class ClosedRadiation(unittest.TestCase):
    def test_friedmann_equations_and_both_phases(self):
        for ell in (0.7, 1.0, 2.3):
            for eta in (0.2, 0.9, math.pi / 2, 2.0, 2.9):
                a = ell * math.sin(eta)
                adot = math.cos(eta) / math.sin(eta)
                addot = -ell**2 / a**3
                self.assertAlmostEqual(adot**2 + 1, ell**2 / a**2, places=10)
                self.assertAlmostEqual(addot / a, -ell**2 / a**4, places=10)
                self.assertAlmostEqual((ell * math.cos(eta))**2 + a**2, ell**2)
                self.assertAlmostEqual(a**2, (ell * (1 - math.cos(eta))) * (ell * (1 + math.cos(eta))))
                if eta != math.pi / 2:
                    self.assertEqual(adot > 0, eta < math.pi / 2)
                else:
                    self.assertAlmostEqual(adot, 0)

    def test_stored_einstein_tensor_has_the_fixed_magnetic_source(self):
        for entry in charts():
            reader = vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
            time = reader.symbol[vm.Reader._plain(entry["coords"][0])]
            ell = reader.parameters["ell"]
            stored = entry["einstein_tensor"]["variants"]["ul"]["nonzero"]
            self.assertEqual(len(stored), 4)
            for eta, lv in ((0.6, 0.8), (math.pi / 2, 1.0), (2.5, 1.7)):
                at = {ell: lv, time: lv * (1 - math.cos(eta)) if entry["id"].startswith("cosmic") else eta}
                a = lv * math.sin(eta)
                for component in stored:
                    expected = (-3 if component["indices"][0] == entry["coords"][0] else 1) * lv**2 / a**4
                    self.assertAlmostEqual(float(reader(component["value"]).subs(at)), expected, delta=1e-9)

    def test_coordinate_map_preserves_the_full_round_sphere(self):
        cosmic, conformal = charts()
        self.assertEqual(cosmic["coords"][1:], conformal["coords"][1:])
        self.assertIn(r"\chi \in [0, \pi]", cosmic["domains"])
        def metric_at(entry, time, ell, chi):
            reader = vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
            at = {reader.parameters["ell"]: ell,
                  **{reader.symbol[vm.Reader._plain(c)]: value
                     for c, value in zip(entry["coords"], (time, chi, 1.1, 0.4))}}
            metric = np.zeros((4, 4))
            for component in entry["metric_components"]:
                i, j = [entry["coords"].index(c) for c in component["indices"]]
                metric[i, j] = float(reader(component["value"]).subs(at))
            return metric
        for eta, chi in ((0.4, 0.2), (1.3, 1.8), (2.8, 2.6)):
            ell = 1.4
            t = ell * (1 - math.cos(eta))
            a = ell * math.sin(eta)
            J = np.diag([a, 1, 1, 1])
            np.testing.assert_allclose(J.T @ metric_at(cosmic, t, ell, chi) @ J,
                                       metric_at(conformal, eta, ell, chi), atol=1e-12)
        for entry in (cosmic, conformal):
            self.assertEqual(entry["ricci_scalar"], "0")


if __name__ == "__main__":
    unittest.main()
