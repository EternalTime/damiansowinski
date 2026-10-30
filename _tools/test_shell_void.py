"""Connection-based junctions and accretion for the cosmological dust shell.

Both metrics and their connections come from the stored smooth tensors.
No bulk scalar is evaluated on the distributional shell.
"""
import json
import math
import sys
import unittest
from pathlib import Path

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / '_tools' / 'derivations'))
import shell_void as shell
import verify_metrics as vm


class Junction(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        metric = json.loads((ROOT / 'MFS/assets/data/metrics/shell_void.json').read_text())
        cls.sides = {}
        for chart in metric['coordinates']:
            coords = chart['coords']
            rd = vm.Reader(coords, [p['symbol'] for p in chart['parameters']], ())
            t, r = [rd.symbol[c] for c in coords[:2]]
            at = {rd.c: 1, rd.symbol[r'\theta']: sp.pi/2}
            if chart['id'] == 'exterior':
                at[rd.parameters['a']] = t ** sp.Rational(2,3)
                at[rd.parameters['k']] = 0
            def expr(s):
                return sp.sympify(rd(s)).subs(at).doit()
            g = {tuple(e['indices']): expr(e['value']) for e in chart['metric_components']}
            connection = {tuple(e['indices']): expr(e['value']) for e in chart['christoffel']['variants']['ull']['nonzero']}
            h = sp.Matrix(2,2,lambda i,j:g.get((coords[i],coords[j]),0))
            G = sp.Array([[[connection.get((coords[a],coords[b],coords[c]),0) for c in range(2)]
                           for b in range(2)] for a in range(2)])
            angular = [connection.get((coords[a],r'\theta',r'\theta'),0) for a in range(2)]
            cls.sides[chart['id']] = (sp.lambdify((t,r),h,'numpy'),sp.lambdify((t,r),G.tolist(),'numpy'),
                                     sp.lambdify((t,r),angular,'numpy'),sp.lambdify((t,r),g[(r'\theta',r'\theta')],'numpy'))

    @staticmethod
    def velocity(side,t,R,v):
        H,rho,gamma,B,sigma = shell.quantities(t,R,v)
        return np.array([gamma, gamma*v/shell.scale(t)]) if side=='exterior' else np.array([gamma*B,gamma*(H*R+v)])

    def curvature(self,side,t,y,perturb=0):
        R,v,T=y
        H,rho,gamma,B,sigma=shell.quantities(t,R,v)
        rate=shell.rhs(t,y)
        rate[1]+=perturb
        # Complex differentiation of the velocity along the shell, d/dtau = gamma d/dt.
        eps=1e-24
        u=self.velocity(side,t,R,v)
        du=gamma*np.imag(self.velocity(side,t+1j*eps,R+1j*eps*rate[0],v+1j*eps*rate[1]))/eps
        time,radius=(t,R/shell.scale(t)) if side=='exterior' else (T,R)
        h,G,angular,area=(f(time,radius) for f in self.sides[side])
        h,G=np.asarray(h,float),np.asarray(G,float)
        n=np.sqrt(-np.linalg.det(h))*np.array([-u[1],u[0]])
        accel=du+np.einsum('abc,b,c->a',G,u,u)
        return -np.dot(n,angular)/area,np.dot(n,accel),u@h@u,n@np.linalg.solve(h,n),area

    def residuals(self,t,perturb=0):
        y=shell.state(t)
        R,v,T=y
        H,rho,gamma,B,sigma=shell.quantities(t,R,v)
        inside=self.curvature('interior',t,y,perturb)
        outside=self.curvature('exterior',t,y,perturb)
        return inside,outside,4*math.pi*sigma

    def test_jumps_and_induced_metric_at_unfitted_points(self):
        for t in (1.0,1.183,1.729,2.417,3.083,4.0):
            with self.subTest(t=t):
                inside,outside,stress=self.residuals(t)
                R,v,T=shell.state(t)
                self.assertGreater(stress,0)
                self.assertLess(abs(v),1)
                for side in (inside,outside):
                    self.assertAlmostEqual(side[2],-1,places=11)
                    self.assertAlmostEqual(side[3],1,places=11)
                    self.assertAlmostEqual(side[4],R*R,places=11)
                self.assertAlmostEqual(outside[0]-inside[0],-stress,places=10)
                self.assertAlmostEqual(outside[1]-inside[1],stress,places=10)

    def test_surface_mass_flux_from_the_angular_constraint(self):
        for t in (1.0,1.381,2.139,3.617,4.0):
            R,v,T=shell.state(t)
            dR,dv,dT=shell.rhs(t,[R,v,T])
            eps=1e-24
            def mass(t,R,v):
                return 4*np.pi*R*R*shell.quantities(t,R,v)[4]
            derivative=np.imag(mass(t+1j*eps,R+1j*eps*dR,v+1j*eps*dv))/eps
            H,rho,gamma,B,sigma=shell.quantities(t,R,v)
            self.assertAlmostEqual(derivative,4*np.pi*R*R*rho*gamma*v,places=11)

    def test_changed_acceleration_is_detected(self):
        inside,outside,stress=self.residuals(2.417,perturb=.03)
        self.assertGreater(abs(outside[1]-inside[1]-stress),1e-5)

    def test_shell_density_is_not_set_to_removed_mass(self):
        R,v,T=shell.state(3.083)
        H,rho,gamma,B,sigma=shell.quantities(3.083,R,v)
        self.assertGreater(abs(4*np.pi*R*R*sigma-4*np.pi*rho*R**3/3),1e-7)
        self.assertGreater(v,0)


if __name__=='__main__':
    unittest.main()
