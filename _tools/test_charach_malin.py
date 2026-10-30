"""Scalar and electromagnetic matter on the Gowdy torus, with periodic initial data."""
import math
from pathlib import Path
import sys
import unittest

import numpy as np
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))
import chart_printer as cp
import null_rays as nr
import print_charts as pc
from test_build_mfs_data import checks


@checks("charach_malin")
class ToroidalMatter(unittest.TestCase):
    def test_each_chart_obeys_its_sourced_matter_equations(self):
        for system in pc.CHARACH_MALIN_CHARTS:
            with self.subTest(system=system):
                spec = pc.charach_malin(system)
                chart = cp.Chart(spec["system"]["coords"], spec["system"]["parameters"], spec["chart_line_element"])
                pc.charach_malin_check(chart, system)

    def test_logarithmic_metrics_are_areal_pullbacks(self):
        for matter in ("scalar", "electromagnetic"):
            specs = [pc.charach_malin(matter+"_"+time) for time in ("areal", "logarithmic")]
            charts = [cp.Chart(s["system"]["coords"], s["system"]["parameters"], s["chart_line_element"]) for s in specs]
            for tau in (-1.3, 0.2, 1.7):
                t = math.exp(-tau)
                matrices = []
                for chart, time in zip(charts, (t, tau)):
                    at = {chart.symbols[0]:time, chart.reader.parameters['L']:1.7,
                          chart.reader.parameters['P']:0.37, chart.reader.parameters['lambda']:-0.23}
                    matrices.append(np.array(chart.geo.g.subs(at), dtype=float))
                J = np.diag([-t, 1, 1, 1])
                np.testing.assert_allclose(J.T @ matrices[0] @ J, matrices[1], rtol=2e-14, atol=2e-14)

    @staticmethod
    def profiles(matter):
        t, theta = sp.symbols('t theta', positive=True)
        fields = {name:sp.sympify(text, locals={'t':t, 'theta':theta})
                  for name,text in nr._charach_malin_wave(matter).items()}
        return t, theta, fields

    def test_the_drawn_standing_waves_satisfy_both_local_constraints(self):
        for matter in ('scalar','electromagnetic'):
            t, h, fields = self.profiles(matter)
            P, lam = fields['P'], fields['lambda']
            f = fields['Phi' if matter=='scalar' else 'chi']
            pt, ph, ft, fh = [sp.diff(a,b) for a,b in ((P,t),(P,h),(f,t),(f,h))]
            if matter=='scalar':
                residuals = [sp.diff(P,t,2)+pt/t-sp.diff(P,h,2),
                             sp.diff(f,t,2)+ft/t-sp.diff(f,h,2),
                             sp.diff(lam,t)-t*(pt**2+ph**2+2*ft**2+2*fh**2),
                             sp.diff(lam,h)-2*t*(pt*ph+2*ft*fh)]
            else:
                residuals = [sp.diff(P,t,2)+pt/t-sp.diff(P,h,2)-2*sp.exp(P)/t*(ft**2-fh**2),
                             sp.diff(f,t,2)-sp.diff(f,h,2)+pt*ft-ph*fh,
                             sp.diff(lam,t)-t*(pt**2+ph**2)-4*sp.exp(P)*(ft**2+fh**2),
                             sp.diff(lam,h)-2*t*pt*ph-8*sp.exp(P)*ft*fh]
            for i, residual in enumerate(residuals):
                value = sp.lambdify((t,h),residual,['scipy','numpy'])
                for T in (0.2,0.7,1.0,2.3,5.0):
                    error=np.max(np.abs(value(T,np.linspace(0,2*math.pi,63))))
                    self.assertLess(float(error),2e-11,(matter,i,T))

    def test_toroidal_identifications_and_zero_total_momentum(self):
        from scipy.integrate import quad
        for matter in ('scalar','electromagnetic'):
            t,h,fields=self.profiles(matter)
            for T in (0.2,1.0,5.0):
                for name,field in fields.items():
                    value=sp.lambdify((t,h),field,['scipy','numpy'])
                    for angle in (0.1,1.3,4.7):
                        self.assertAlmostEqual(value(T,angle),value(T,angle+2*math.pi),delta=2e-13,
                                               msg=(matter,name,T))
                momentum=sp.lambdify((t,h),sp.diff(fields['lambda'],h),['scipy','numpy'])
                integral=quad(lambda angle:momentum(T,angle),0,2*math.pi,epsabs=1e-12)[0]
                self.assertLess(abs(integral),1e-11,(matter,T))
            field=sp.lambdify((t,h),fields['Phi' if matter=='scalar' else 'chi'],['scipy','numpy'])
            self.assertGreater(abs(field(1,0)-field(1,math.pi)),0.1)


if __name__ == '__main__':
    unittest.main()
