"""Source domain, regular axis, and central rotation of the bounded 1984 fluid.

Tensor generation already checks the source stress and the complete coordinate pullback.
These checks cover the endpoints and the physical claims those point checks do not establish.
"""
import json
import math
from pathlib import Path
import sys
import unittest

import mpmath as mp
import sympy as sp

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT/'_tools'/'derivations'))
import verify_metrics as vm


class KramerFluid(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = json.loads((ROOT/'MFS/assets/data/metrics/kramer_fluid.json').read_text())
        cls.chart = next(c for c in cls.data['coordinates'] if c['id']=='radial')
        cls.reader = vm.Reader(cls.chart['coords'], [p['symbol'] for p in cls.chart['parameters']], (),
                               held=vm.HELD[('kramer_fluid','radial')])
        cls.g = vm.metric_from_line_element(cls.reader, cls.chart['line_element'], cls.chart['coords'])
        cls.g = cls.g.applyfunc(lambda e: e.subs(cls.reader.held).doit())

    def test_strict_parameter_endpoints_and_pressure_boundary(self):
        with mp.workdps(50):
            eta_star = mp.findroot(lambda a: (a+1)*(2*mp.exp(a-1)-mp.exp(2*a-2))-2, (1.2,1.5))
            b_star = (eta_star+1)*mp.exp(eta_star)
            self.assertAlmostEqual(float(b_star), 8.603, places=3)
            # The upper endpoint collapses the equatorial orbit and is excluded by every chart.
            def orbit(b, eta, xi):
                return 4*(eta-1-b*mp.cos(xi)*(mp.exp(-eta)+mp.exp(eta-2)-2/mp.e))
            self.assertLess(abs(orbit(b_star,eta_star,0)),mp.mpf('1e-45'))
            self.assertLess(orbit(b_star+mp.mpf('.01'),eta_star,0),0)
            for chart in self.data['coordinates']:
                self.assertIn(r'2e<b<b_\star',''.join(chart['domains']))
            for b in (mp.mpf('5.5'),mp.mpf('6.1'),mp.mpf('8.5')):
                edge = mp.findroot(lambda a:(a+1)*mp.exp(a)-b,(1,1.5))
                for f in (mp.mpf('.01'),mp.mpf('.5'),mp.mpf('.99')):
                    eta = 1+(edge-1)*f
                    bound = mp.acos((eta+1)*mp.exp(eta)/b)
                    for fraction in (mp.mpf('-.99'),mp.mpf('0'),mp.mpf('.99')):
                        xi = bound*fraction
                        F = b*mp.cos(xi)*mp.exp(-eta)
                        p, density = (1+eta-F)/2,(1-3*eta+3*F)/2
                        self.assertLess(p,0)
                        self.assertGreater(density,abs(p))
                        self.assertLess(abs(density+3*p-2),mp.mpf('1e-45'))
                        self.assertGreater(orbit(b,eta,xi),0)
                    self.assertLess(abs(b*mp.cos(bound)*mp.exp(-eta)-1-eta),mp.mpf('1e-45'))
            for b in (mp.e, 2*mp.e):
                self.assertLessEqual(b/mp.e-2, 0)  # Even the central maximum F has no strict interior.

    def test_elementary_flatness_and_regular_twist_on_the_axis(self):
        R = self.reader
        w, xi = R.symbol['w'], R.symbol[r'\xi']
        m,b = R.parameters['m'],R.parameters['b']
        for latitude in ('-.3','0','.3'):
            at = {w:sp.Float('1e-12',70),xi:sp.Float(latitude,70),m:sp.Rational(1,2),b:sp.Rational(61,10)}
            gpp = self.g[3,3].subs(at).evalf(60)
            # proper radial distance=2 ell w, so circumference/radius ->2 pi.
            self.assertLess(abs(gpp/(4*at[w]**2)-1),sp.Float('1e-20'))
            twist = self.g[0,3].subs(at).evalf(60)/at[w]**2
            expected = 2*(1+at[b]*sp.cos(at[xi])/sp.E)
            self.assertLess(abs(twist-expected),sp.Float('1e-20'))
            self.assertEqual(self.g[1,1].subs(at),4)

    def test_vorticity_at_the_centre_and_no_static_endpoint(self):
        R = self.reader
        w,xi = R.symbol['w'],R.symbol[r'\xi']
        m,b = R.parameters['m'],R.parameters['b']
        # u is the unit stationary Killing direction. Its projected exterior derivative is
        # twice the vorticity, independently of the source's expression for its magnitude.
        u = self.g[:,0]/sp.sqrt(-self.g[0,0])
        du = sp.Matrix(4,4,lambda i,j:(sp.diff(u[j],R.symbol[self.chart['coords'][i]])
                                      -sp.diff(u[i],R.symbol[self.chart['coords'][j]]))/2
        for bv in (sp.Rational(11,2),sp.Rational(61,10),sp.Rational(17,2)):
            at = {w:sp.Float('1e-9',65),xi:0,m:sp.Rational(1,2),b:bv}
            g = self.g.subs(at).evalf(55)
            cov = u.subs(at).evalf(55)
            up = sp.Matrix([1/sp.sqrt(-g[0,0]),0,0,0])
            project = sp.eye(4)+cov*up.T
            omega = project*du.subs(at).evalf(55)*project.T
            raised = g.inv()*omega*g.inv()
            square = sum(omega[i,j]*raised[i,j] for i in range(4) for j in range(4))/2
            expected = (bv+sp.E)**2*at[m]/(2*bv*sp.E)
            self.assertLess(abs(square-expected).evalf(40),sp.Float('1e-14'))
            self.assertGreater(expected.evalf(),0)


if __name__=='__main__':
    unittest.main()
