"""Momentum tensors and numerical spin-throat constraints of Bowen-York initial data.

The chart printer checks generic intrinsic tensors and coordinate maps.
The embedding generator checks distances on the drawn surface.
These checks cover the remaining extrinsic data and elliptic boundary problem.
"""
import math
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "_tools" / "derivations"))
from bowen_york_data import solution


def tensor(x, p, s, image=0):
    """Brandt-Bruegmann (5), plus either sign of Cook's a=1 image term (72)."""
    r = np.linalg.norm(x)
    n = np.asarray(x)/r
    pn = np.dot(p,n)
    cross = np.cross(s,n)
    nn = np.outer(n,n)
    momentum = 1.5/r**2*(np.outer(p,n)+np.outer(n,p)-(np.eye(3)-nn)*pn)
    spin = 3/r**3*(np.outer(cross,n)+np.outer(n,cross))
    mirror = image*1.5/r**4*(np.outer(p,n)+np.outer(n,p)+(np.eye(3)-5*nn)*pn)
    return momentum+spin+mirror


class MomentumConstraint(unittest.TestCase):
    def test_trace_and_divergence_including_both_image_signs(self):
        p, s = np.array([0.7,-0.3,0.5]), np.array([-0.2,0.9,0.4])
        for x in (np.array([1.3,0.6,0.7]), np.array([-2.1,1.2,0.3])):
            for image in (-1,0,1):
                A = tensor(x,p,s,image)
                self.assertLess(abs(np.trace(A)), 2e-15)
                h = 1e-4
                # Fourth-order finite differences differentiate the sourced Cartesian tensor.
                divergence = np.zeros(3)
                for j in range(3):
                    e = np.eye(3)[j]*h
                    derivative = (-tensor(x+2*e,p,s,image)+8*tensor(x+e,p,s,image)
                                  -8*tensor(x-e,p,s,image)+tensor(x-2*e,p,s,image))/(12*h)
                    divergence += derivative[:,j]
                self.assertLess(np.max(np.abs(divergence)), 1e-10)

    def test_spin_contraction_and_spherical_component(self):
        for r,t,phi in ((1.2,0.8,0.4),(3.1,1.4,2.3),(1.0,math.pi/2,0.1)):
            n=np.array([math.sin(t)*math.cos(phi),math.sin(t)*math.sin(phi),math.cos(t)])
            dphi=r*np.array([-math.sin(t)*math.sin(phi),math.sin(t)*math.cos(phi),0])
            A=tensor(r*n,np.zeros(3),np.array([0,0,1]))
            self.assertAlmostEqual(np.sum(A*A),18*math.sin(t)**2/r**6,places=13)
            self.assertAlmostEqual(n@A@dphi,3*math.sin(t)**2/r**2,places=13)

    def test_surface_integrals_identify_single_hole_momenta(self):
        p,s=np.array([0.7,-0.3,0.5]),np.array([-0.2,0.9,0.4])
        mu,w=np.polynomial.legendre.leggauss(8)
        linear,angular=np.zeros(3),np.zeros(3)
        r=3.0
        for m,weight in zip(mu,w):
            for phi in np.linspace(0,2*math.pi,32,endpoint=False):
                n=np.array([math.sqrt(1-m*m)*math.cos(phi),math.sqrt(1-m*m)*math.sin(phi),m])
                flux=tensor(r*n,p,s)@n*r*r*weight*(2*math.pi/32)/(8*math.pi)
                linear+=flux
                angular+=np.cross(r*n,flux)
        np.testing.assert_allclose(linear,p,atol=2e-14)
        np.testing.assert_allclose(angular,s,atol=2e-14)


class SpinThroatConstraint(unittest.TestCase):
    def test_off_grid_hamiltonian_boundary_conditions_and_resolution(self):
        w=solution()
        finer=solution(nodes=40,modes=10)
        r,mu=np.meshgrid(np.geomspace(1,100,29),np.linspace(-1,1,23))
        self.assertTrue(np.all(w.field(r,mu)>0))
        self.assertLess(np.max(np.abs(w.residual(r,mu))),2e-9)
        self.assertLess(np.max(np.abs(w.field(r,mu)-finer.field(r,mu))),2e-11)
        angles=np.linspace(-1,1,17)
        self.assertLess(np.max(np.abs(w.field(1,angles,1)+w.field(1,angles)/2)),2e-11)
        self.assertLess(np.max(np.abs(w.field(1e9,angles)-1)),2e-9)

    def test_zero_spin_and_inversion_extension(self):
        w=solution(spin=0)
        r=np.geomspace(1,40,30)
        np.testing.assert_allclose(w.field(r,0.3),1+1/r,atol=2e-14)
        w=solution()
        # Check the Hamiltonian source on the other end using the inverted numerical solution.
        # Finite differences are independent of the exterior spectral derivatives.
        for r in (0.2,0.6,0.9):
            for mu in (-0.8,0,0.6):
                def f(rv,m):
                    return float(w.field(1/rv,m))/rv
                h=2e-4
                def first(axis):
                    values=[f(r+k*h,mu) if axis==0 else f(r,mu+k*h) for k in (-2,-1,1,2)]
                    return (values[0]-8*values[1]+8*values[2]-values[3])/(12*h)
                def second(axis):
                    values=[f(r+k*h,mu) if axis==0 else f(r,mu+k*h) for k in (-2,-1,0,1,2)]
                    return (-values[0]+16*values[1]-30*values[2]+16*values[3]-values[4])/(12*h*h)
                lap=second(0)+2*first(0)/r+((1-mu*mu)*second(1)-2*mu*first(1))/r**2
                source=2.25*(1-mu*mu)/r**6/f(r,mu)**7
                self.assertLess(abs(lap+source),2e-6)


if __name__ == '__main__':
    unittest.main()
