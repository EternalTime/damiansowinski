"""Circumference bound, particle stress, and Israel jumps from the stored cylindrical metrics."""

import json
import random
import sys
import unittest
from pathlib import Path

import mpmath as mp
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "_tools/derivations"))
import verify_metrics as vm
import rotation_halts as model


class RotationHalts(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mp.mp.dps = 80
        cls.metric = json.loads((ROOT / "MFS/assets/data/metrics/rotation_halts.json").read_text())

    def test_bound_follows_from_positive_lapses_and_a_regular_axis(self):
        rng = random.Random(1992)
        admitted = 0
        for _ in range(400):
            alpha = mp.mpf(str(rng.uniform(0.1, 3)))
            lam = mp.mpf(str(rng.uniform(0.001, 0.15)))
            R = mp.mpf(str(rng.uniform(0.01, 8)))
            psi = mp.mpf(str(rng.uniform(-1, 1)))
            gamma_in = mp.mpf(str(rng.uniform(0, 1)))
            V = mp.mpf(str(rng.uniform(-1, 1)))
            Xm = mp.sqrt(mp.exp(2*psi-2*gamma_in) + V*V)
            Xp = Xm - 4*lam*mp.sqrt(R*R + alpha*alpha*mp.exp(2*psi))/R
            if Xp <= abs(V):
                continue
            admitted += 1
            self.assertLess(Xm-Xp, mp.exp(psi))
            self.assertGreater(R, 4*alpha*lam)
            C, Lambda = 2*mp.pi*R*mp.exp(-psi), lam*mp.exp(-psi)
            self.assertGreater(C, 8*mp.pi*alpha*Lambda)
        self.assertGreater(admitted, 30)
        for alpha in map(mp.mpf, ["0.1", "1", "3"]):
            for lam in map(mp.mpf, ["0.01", "0.1"]):
                for factor in map(mp.mpf, ["0.1", "0.9", "1"]):
                    R = 4*alpha*lam*factor
                    # Even the largest interior lapse (gamma-=0) is insufficient.
                    Xp = 1 - 4*lam*mp.sqrt(R*R+alpha*alpha)/R
                    self.assertLess(Xp, 0)

    def test_flat_interior_is_a_restriction_and_not_a_generic_vacuum(self):
        entry = next(c for c in self.metric["coordinates"] if c["id"] == "interior")
        rd = vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
        psi, gamma = rd.parameters["psi"], rd.parameters["gamma"]
        t, rho = rd.symbol["t"], rd.symbol["\\rho"]
        tensor = entry["riemann"]["variants"]["ulll"]["nonzero"]
        self.assertTrue(tensor)
        flat = {psi: sp.S.Zero, gamma: sp.S.Zero}
        wave = {psi: t/7, gamma: rho*rho/98}
        maximum = 0
        for component in tensor:
            value = rd(component["value"])
            zero = value.subs(flat).doit().subs({t: 1, rho: 2})
            self.assertEqual(zero, 0)
            curved = value.subs(wave).doit().subs({t: 1, rho: 2}).evalf(60)
            maximum = max(maximum, abs(complex(curved)))
        self.assertGreater(maximum, 0.001)


if __name__ == "__main__":
    unittest.main()
