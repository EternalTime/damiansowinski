"""The semiclassical equation, its regular trace evolution and its numerical realization.

The chart printer checks the FLRW Ricci normalization.
The collection checker checks every stored tensor and drawing, so those are not repeated here.
"""
import json
import sys
import unittest
from pathlib import Path

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "_tools" / "derivations"))
import starobinsky_universe as su


class SemiclassicalEquation(unittest.TestCase):
    def test_covariant_00_normalization(self):
        h, hp, hpp, K = sp.symbols("h hp hpp K")
        Q, D = h*h+K, hp+h*h
        Rp = 6*(hpp+4*h*hp-2*h*K)
        bracket = 2*h*hpp+6*h*h*hp-hp*hp-2*K*h*h+K*K
        self.assertEqual(sp.expand(Q*Q-D*D+h*Rp/3-bracket), 0)

    def test_trace_preserves_constraint_for_every_spatial_curvature(self):
        h, R, S, K, m = sp.symbols("h R S K m")
        Q = h*h+K
        C = m*m*(Q-Q*Q)+Q*R/3-R*R/36+h*S/3
        hp = R/6-2*h*h-K
        Sp = -3*h*S-m*m*R+12*m*m*Q*(R/6-Q)
        derivative = sp.diff(C,h)*hp+sp.diff(C,R)*S+sp.diff(C,S)*Sp-2*h*K*sp.diff(C,K)
        self.assertEqual(sp.expand(derivative+4*h*C), 0)
        # A changed anomaly coefficient must fail this identity.
        wrong = derivative+sp.diff(C,S)*m*m*Q*(R/6-Q)+4*h*C
        self.assertNotEqual(sp.expand(wrong), 0)

    def test_founding_phase_plane_equation(self):
        z, y, yx, M, H = sp.symbols("z y yx M H", positive=True)
        h, hp, hpp = H*z*z, 2*H*H*z*y, 2*H**3*y*y+6*H**3*z**3*y*yx
        full = M*M*(h*h-h**4/H**2)+2*h*hpp+6*h*h*hp-hp*hp
        target = 12*H**4*z**5*(y*yx+y-M*M/(12*H*H)*(z**3-1/z))
        self.assertEqual(sp.expand(full-target), 0)

    def test_departure_exponent_and_second_order_initialization(self):
        q, u = sp.symbols("q u")
        b = q*(4*q+3)/(6*(q+1))
        h, f = 1+u, q*u+b*u*u
        residual = (q*q+3*q)*(h*h-h**4)+2*h*sp.diff(f,u)*f+6*h*h*f-f*f
        for order in (0,1,2):
            self.assertEqual(sp.factor(sp.expand(residual).coeff(u,order)),0)
        self.assertGreater(su.EXPONENT,0)
        self.assertLess(su.EXPONENT,2)


class NumericalBranch(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.sol = su.trajectory()

    def test_constraint_and_initialization_convergence(self):
        samples = np.linspace(-6,160,2001)
        N,h,R,S,eta = self.sol.state(samples)
        residual = float(np.max(np.abs(su.constraint(h,R,S))))
        self.assertLess(residual,2e-10)
        refined = su.Trajectory(displacement=su.DISPLACEMENT/2,tolerance=5e-13)
        fine = refined.state(samples)
        # Normalize the time and a at the same h=1/sqrt(2) crossing.
        gap = float(np.max(np.abs(fine[:4]-np.array([N,h,R,S]))))
        self.assertLess(gap,2e-8)
        print(f"Starobinsky trajectory: constraint {residual:.3e}; refined N,h,R,S gap {gap:.3e}")

    def test_conformal_quadrature_inverse_and_oscillations(self):
        samples = np.linspace(-6,150,157)
        eta = self.sol.eta(samples)
        self.assertTrue(np.all(np.diff(eta)>0))
        gap = float(np.max(np.abs(self.sol.time(eta)-samples)))
        self.assertLess(gap,2e-6)
        step = 1e-4
        for t in (-3,0,5,30,120):
            slope = (self.sol.eta(t+step)-self.sol.eta(t-step))/(2*step)
            self.assertLess(abs(slope*self.sol.a(t)-1),2e-7)
        R = self.sol.state(np.linspace(20,160,2001))[2]
        self.assertGreater(np.count_nonzero(R[:-1]*R[1:]<0),10)
        print(f"Starobinsky conformal inverse: worst time gap {gap:.3e}")


if __name__ == "__main__":
    unittest.main()
