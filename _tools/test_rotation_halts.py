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
        cls.faces = {side: cls.face(side) for side in ("interior", "exterior")}

    @classmethod
    def face(cls, side):
        """Read the metric and connection, keeping independently specified boundary jets."""
        entry = next(c for c in cls.metric["coordinates"] if c["id"] == side)
        rd = vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
        t, rho = rd.symbol["t"], rd.symbol["\\rho"]
        psi, gamma = rd.parameters["psi"], rd.parameters["gamma"]
        g0, pr, pt, gr, gt = sp.symbols("g0 pr pt gr gt")
        jets = {sp.diff(psi, rho): pr, sp.diff(psi, t): pt,
                sp.diff(gamma, rho): gr, sp.diff(gamma, t): gt}

        def compile_value(text):
            value = rd(text).xreplace(jets).subs({psi: 0, gamma: g0})
            return sp.lambdify((rho, g0, pr, pt, gr, gt), value, "mpmath")

        metric = {tuple(v["indices"]): compile_value(v["value"])
                  for v in entry["metric_components"]}
        connection = {tuple(v["indices"]): compile_value(v["value"])
                      for v in entry["christoffel"]["variants"]["ull"]["nonzero"]}
        return metric, connection

    def face_curvature(self, side, R, V, X, A, Xdot, jet):
        """Pull back the stored metric and connection with the same outward normal."""
        g, connection = self.faces[side]
        args = (R,) + jet
        names = ("t", "\\rho", "\\phi", "z")
        U, n = (X, V), (V, X)
        covnormal = [sum(g.get((names[i], names[j]), lambda *a: 0)(*args)*n[j]
                        for j in range(2)) for i in range(2)]

        def G(a, b, c):
            fn = connection.get((a, b, c))
            return fn(*args) if fn else mp.mpf(0)

        norm = sum(g.get((names[i], names[j]), lambda *a: 0)(*args)*U[i]*U[j]
                   for i in range(2) for j in range(2))
        self.assertLess(abs(norm+1), mp.mpf("1e-65"))
        acc = [q + sum(G(names[a], names[b], names[c])*U[b]*U[c]
                       for b in range(2) for c in range(2))
               for a, q in enumerate((Xdot, A))]
        Ktt = -sum(covnormal[a]*acc[a] for a in range(2))
        Kpp = -sum(covnormal[a]*G(names[a], "\\phi", "\\phi") for a in range(2))/R**2
        Kzz = -sum(covnormal[a]*G(names[a], "z", "z") for a in range(2))
        return (Ktt, Kpp, Kzz)

    def matching_case(self, R, V, alpha, lam):
        R, V, alpha, lam = map(mp.mpf, (R, V, alpha, lam))
        s = mp.sqrt(R*R+alpha*alpha)
        Xm = mp.sqrt(1+V*V)
        Xp = Xm - 4*lam*s/R
        self.assertGreater(Xp, abs(V))
        D = Xp*Xp-V*V
        P = -2*lam/s
        pr, pt = Xp*P/D, -V*P/D
        gr, gt = R*(pr*pr+pt*pt), 2*R*pr*pt
        A = ((1+V*V)*alpha*alpha/(R*(R*R+alpha*alpha))
             - lam*Xm*(R*R+2*alpha*alpha)**2/(R*R*s**3))
        Xmdot = V*A/Xm
        Xpdot = Xmdot + 4*lam*alpha*alpha*V/(s*R*R)
        inner = self.face_curvature("interior", R, V, Xm, A, Xmdot, (0, 0, 0, 0, 0))
        outer = self.face_curvature("exterior", R, V, Xp, A, Xpdot,
                                    (-mp.log(D)/2, pr, pt, gr, gt))
        # The stress is independently obtained by summing the equal dust streams.
        u = alpha/R
        sigma = lam*mp.sqrt(1+u*u)/(2*mp.pi*R)
        pressure = lam*u*u/(2*mp.pi*R*mp.sqrt(1+u*u))
        trace = -sigma+pressure
        rhs = [-8*mp.pi*(sij-eta*trace/2)
               for sij, eta in zip((sigma, pressure, 0), (-1, 1, 1))]
        residual = max(abs(kp-km-source) for kp, km, source in zip(outer, inner, rhs))
        self.assertLess(residual, mp.mpf("1e-65"))
        self.assertLess(abs(Xp*pr+V*pt-P), mp.mpf("1e-65"))
        self.assertLess(abs(V*pr+Xp*pt), mp.mpf("1e-65"))
        return residual, A

    def test_complete_israel_jump_from_published_christoffel_symbols(self):
        maximum = mp.mpf(0)
        for R, V, alpha, lam in [("2", "0", "1", ".08"), ("2.1", ".1", "1", ".08"),
                                  ("2.1", "-.1", "1", ".08"), ("3", ".2", ".5", ".03"),
                                  ("3", "-.2", ".5", ".03")]:
            residual, _ = self.matching_case(R, V, alpha, lam)
            maximum = max(maximum, residual)
        print("maximum full Israel-jump residual:", mp.nstr(maximum, 5))

    def test_equilibrium_and_the_turning_event(self):
        for R, alpha in [(mp.mpf(2), mp.mpf(1)), (mp.mpf(3), mp.mpf("0.5"))]:
            lam = alpha*alpha*R*mp.sqrt(R*R+alpha*alpha)/(R*R+2*alpha*alpha)**2
            _, A = self.matching_case(R, 0, alpha, lam)
            self.assertLess(abs(A), mp.mpf("1e-65"))
        _, A = self.matching_case("2", "0", "1", ".08")
        self.assertGreater(A, 0)
        # With no rotation the same restricted shell accelerates inward.
        self.assertLess(model.acceleration(2, 0, 0, .08), 0)

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
