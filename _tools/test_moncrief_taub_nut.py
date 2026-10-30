"""Analytic horizon regularity and the inhomogeneous Bessel drawing member.

The printer checks the full two-amplitude vacuum equations on generic jets and both pullbacks.
These checks cover the distinct horizon limit and the explicit member used for the drawings.
"""
import json
import sys
import unittest
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT/'_tools'/'derivations'))
import moncrief_taub_nut as moncrief
import verify_metrics as vm


class AnalyticHorizon(unittest.TestCase):
    def test_bessel_member_solves_wave_equation_and_both_constraints(self):
        with mp.workdps(40):
            for t,theta in [(mp.mpf('0.3'),mp.mpf('0.7')),(mp.mpf('1.4'),mp.mpf('2.1'))]:
                b=lambda t,theta: moncrief.drawing_functions(t,theta)[0]
                f=lambda t,theta: moncrief.drawing_functions(t,theta)[1]
                bt=mp.diff(lambda v:b(v,theta),t)
                by=mp.diff(lambda v:b(t,v),theta)
                wave=mp.diff(lambda v:b(v,theta),t,2)+bt/t-mp.diff(lambda v:b(t,v),theta,2)
                first=mp.diff(lambda v:f(v,theta),t)-t*(bt*bt+by*by)
                second=mp.diff(lambda v:f(t,v),theta)-2*t*bt*by
                self.assertLess(float(max(abs(wave),abs(first),abs(second))),1e-35)
                self.assertNotEqual(float(by),0)
            self.assertEqual(moncrief.drawing_functions(0,mp.mpf('0.7'))[1],0)

    def test_both_published_charts_have_a_regular_nondegenerate_horizon(self):
        metric=json.loads((ROOT/'MFS/assets/data/metrics/moncrief_taub_nut.json').read_text())
        for entry in metric['coordinates'][1:]:
            R=vm.Reader(entry['coords'],[p['symbol'] for p in entry['parameters']],())
            T=R.symbol['T'];B,F,Q=(R.parameters[n] for n in ('B','F','Q'))
            # F=log(1+4TH/kappa^2)/2 writes out the removable quotient before T=0.
            H=sp.Rational(17,100)
            vals={B:sp.Rational(13,100),Q:sp.Rational(21,100),R.parameters['L']:sp.Rational(6,5),
                  R.parameters['kappa']:sp.Rational(9,10)}
            matrix=sp.zeros(4)
            for c in entry['metric_components']:
                i,j=(entry['coords'].index(x) for x in c['indices'])
                value=R(c['value']).subs(F,sp.log(1+4*T*H/R.parameters['kappa']**2)/2).subs(vals)
                # Cancel the removable quotient before assigning its analytic value at T=0.
                value=sp.cancel(sp.expand_power_exp(value))
                matrix[i,j]=value
                matrix[j,i]=value
            horizon=np.array(matrix.subs(T,0),dtype=float)
            determinant=np.linalg.det(horizon)
            expected=-float(vals[R.parameters['L']]**8*vals[R.parameters['kappa']]**4*sp.exp(4*vals[B]))/4
            self.assertAlmostEqual(determinant/expected,1,places=12)
            self.assertEqual(np.count_nonzero(np.linalg.eigvalsh(horizon)<0),1)
            self.assertGreater(float(matrix[2,2].subs(T,sp.Rational(1,10))),0)
            self.assertLess(float(matrix[2,2].subs(T,-sp.Rational(1,10))),0)

    def test_regular_member_coefficients_continue_to_negative_T(self):
        # The entire Bessel series fixes the extension at negative T and the horizon.
        self.assertEqual(float(moncrief.j0_T(0)),1)
        self.assertEqual(float(moncrief.h0_T(0)),0)
        for T in [-1,-0.01,0.01,1]:
            with mp.workdps(40):
                z=mp.mpf(T)
                a,d=mp.hyp0f1(1,-z/4),mp.hyp0f1(2,-z/4)/2
                f=(z*(a*a+z*d*d)/2-z*a*d)/16
                exact=mp.expm1(2*f)/(4*z)
            self.assertAlmostEqual(float(moncrief.h0_T(T)),float(exact),places=14)


if __name__=='__main__':
    unittest.main()
