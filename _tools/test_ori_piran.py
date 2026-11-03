"""Self similar collapse of a perfect fluid: the facts its page states, held to its numerical solution."""
import json
import math
import sys
import unittest
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "_tools" / "derivations"))
import ori_piran as op  # noqa: E402

K = op.K
METRIC = json.loads((ROOT / "MFS/assets/data/metrics/ori_piran.json").read_text(encoding="utf-8"))


class FlatFriedmann(unittest.TestCase):
    """The member the page gives in closed form, D = 2/(3(1+k)^2)."""

    def test_its_D_is_the_flat_universe_with_pressure_k_rho(self):
        # a = (-t)^(2/(3(1+k))) gives H = 2/(3(1+k)(-t)) and rho = 3H^2/(8 pi), so D = 4 pi rho t^2.
        for k in (Fraction(1, 100), Fraction(1, 3), Fraction(1, 2)):
            hubble_t = Fraction(2, 1) / (3 * (1 + k))
            rho_t2_times_8pi = 3 * hubble_t ** 2
            self.assertEqual(rho_t2_times_8pi / 2, Fraction(2, 1) / (3 * (1 + k) ** 2))

    def test_its_centre_series_solves_the_equations_exactly(self):
        D = 2 / (3 * (1 + K) ** 2)
        for z in (1e-8, 1e-3, 0.0123, 0.5, 30.0):
            Y = np.array(op.centre(D, z))
            M, S, eta, _, y, V2, num = op.parts(math.log(z), Y, D)
            self.assertAlmostEqual(y, 1 / 3, places=12)
            self.assertLess(abs(num), 1e-12 * max(1.0, 4 * y * V2))
            self.assertLess(abs(op.constraint(math.log(z), Y, D)), 1e-12)
            self.assertAlmostEqual(math.exp(2 * K / (1 + K) * math.log(2 * D * z * z / eta)), 1.0, places=12)


class SoundAndDensity(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.D, cls.ys, cls.xs, cls.pieces = op.solve()

    def test_the_sonic_point_moves_at_the_speed_of_sound(self):
        # For P = k rho the speed of sound is sqrt(dP/drho) = sqrt(k); there V_z^2 = k, the
        # regularity condition holds and the constraint with it.
        x, Y = op.sonic_state(self.ys, self.D)
        M, S, eta, z, y, V2, num = op.parts(x, Y, self.D)
        self.assertAlmostEqual(V2, K, places=14)
        self.assertLess(abs(num), 1e-14)
        self.assertLess(abs(op.constraint(x, Y, self.D)), 1e-13)
        j = op.jet(math.exp(self.xs))
        vz2 = math.exp(self.xs) ** 2 * math.exp(float(j["omega"]) - float(j["sigma"]))
        self.assertAlmostEqual(vz2, K, places=7)

    def test_the_two_densities_agree_at_the_centre(self):
        # rho_c = D/(4 pi t^2) and rho = eta/(8 pi r^2) agree as r -> 0, since eta -> 2 D z^2.
        # The next term of the series is of order z^(2(1+3k)/(3(1+k))), about 1e-6 at z = 1e-10.
        for z, places in ((1e-12, 9), (1e-10, 5)):
            lnS, lneta = op.logs(math.log(z))[1:]
            rho_r2_8pi = math.exp(lneta)
            self.assertAlmostEqual(rho_r2_8pi / (2 * self.D * z * z), 1.0, places=places)
        self.assertAlmostEqual(math.exp(float(op.jet(1e-12)["sigma"])), 1.0, places=8)

    def test_the_solution_keeps_its_constraint_and_its_matching(self):
        for piece in self.pieces:
            for x, Y in zip(piece.t[::25], piece.y.T[::25]):
                self.assertLess(abs(op.constraint(x, Y, self.D)), 1e-9)
        inner, middle, _ = self.pieces
        self.assertLess(float(np.max(np.abs(inner.y[1:, -1] - middle.y[1:, -1]))), 1e-9)
        left, right = op.logs(self.xs - 1e-4), op.logs(self.xs + 1e-4)
        self.assertLess(abs((right[2] - left[2]) / 2e-4 - op.rates(self.xs + 1e-3, op.logs(self.xs + 1e-3), self.D)[2]), 2e-2)

    def test_harada_and_maeda_value_of_D(self):
        # Harada and Maeda give D of about 1.439 at k = 0.01; this integration finds 1.43771.
        self.assertLess(abs(self.D / 1.439 - 1), 2e-3)
        self.assertTrue(0 < self.ys < 1)

    def test_no_trapped_surface_and_an_embeddable_slice_before_the_singularity(self):
        for z in np.geomspace(1e-8, 1e6, 60):
            M, S, eta, *_ = op.parts(math.log(z), op.logs(math.log(z)), self.D)
            self.assertLess(M / S, 1.0)
            j = op.jet(z)
            gamma2 = math.exp(-float(j["omega"])) * (float(j["S"]) + z * float(j["dS"])) ** 2
            self.assertLess(gamma2, 1.0)

    def test_jets_are_derivatives_of_the_functions(self):
        # Beyond z of a few the functions change by less than a difference of doubles can see.
        z = np.array([1e-6, 1e-3, 0.005, 0.02, 0.3, 3.0])
        h = 1e-5 * z
        left, middle, right = op.jet(z - h), op.jet(z), op.jet(z + h)
        for name in ("sigma", "omega", "S"):
            for base, derivative in ((name, "d" + name), ("d" + name, "d2" + name)):
                numerical = (right[base] - left[base]) / (2 * h)
                scale = np.maximum(1e-12, np.abs(middle[derivative]))
                self.assertLess(float(np.max(np.abs(numerical - middle[derivative]) / scale)), 1e-4)


class PublishedEinsteinTensor(unittest.TestCase):
    """The comoving chart's published G^mu_nu, evaluated on the numerical solution, equals the
    stresses of the fluid at rest in the chart, 8 pi (rho + P) u^mu u_nu + 8 pi P delta^mu_nu."""

    def test_on_the_solution(self):
        import curvature as cv
        import null_rays as nr

        class Sources:
            def note(self, *args):
                pass
        rho = {}
        for z in (0.002, 0.005, 0.03, 0.3, 2.0):
            lneta = op.logs(math.log(z))[2]
            rho[z] = math.exp(lneta) / (8 * math.pi * z * z)
        for (a, b), want in ((("t", "t"), lambda p: -8 * math.pi * p), (("r", "r"), lambda p: 8 * math.pi * K * p),
                             (("\\theta", "\\theta"), lambda p: 8 * math.pi * K * p), (("t", "r"), lambda p: 0.0)):
            value = cv.Held(cv.Published(Sources(), "ori_piran", "comoving", "einstein_tensor", ["r", "t", "theta"], [],
                                         indices=[a, b], functions=nr.OP_COMOVING), {"t": -1, "theta": 1.1})
            for z, p in rho.items():
                got = float(value(np.array([z]))[0])
                self.assertAlmostEqual(got, want(p), delta=2e-5 * abs(8 * math.pi * p),
                                       msg=f"G^{a}_{b} at z = {z}")
                if a == "r":
                    # The same comparison refuses a fluid whose pressure is twice as large.
                    self.assertGreater(abs(got - 2 * want(p)), 1e-3 * abs(8 * math.pi * p))


class Drawings(unittest.TestCase):
    def test_the_marked_lines_are_the_sonic_point_and_the_null_orbit(self):
        import null_rays as nr
        for z, square in ((nr.OP_SONIC, K), (nr.OP_NULL, 1.0)):
            j = op.jet(z)
            # The rows hold each z to six decimals.
            self.assertAlmostEqual(z * z * math.exp(float(j["omega"]) - float(j["sigma"])), square, delta=1e-4 * square)
        diagrams = json.loads((ROOT / "MFS/assets/data/diagrams/ori_piran.json").read_text(encoding="utf-8"))
        kinds = {sid: sorted(m["kind"] for m in views[0]["markers"]) for sid, views in diagrams["systems"].items()}
        self.assertEqual(kinds["comoving"], ["mark", "past"])
        self.assertIn("mark", kinds["zooming"])
        self.assertIn("past", kinds["zooming"])

    def test_no_conformal_diagram_and_the_page_says_which_quantities_are_numerical(self):
        import conformal
        self.assertIn("ori_piran", conformal.NOT_DRAWN)
        self.assertFalse((ROOT / "MFS/assets/data/conformal/ori_piran.json").exists())
        described = " ".join(p["description"] for c in METRIC["coordinates"] for p in c["parameters"])
        self.assertIn("known only numerically", described)
        self.assertIn("D = 2/(3(1+k)^2)", described)


if __name__ == "__main__":
    unittest.main()
