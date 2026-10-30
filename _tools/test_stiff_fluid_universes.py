"""Independent published-data checks of the canonical scalar and stiff-fluid family.

Every chart is held to Einstein's equation with eight pi explicit, its timelike fluid,
the scalar dictionary, the density and its boundary, and the isotropic Friedmann member.
"""
import importlib.util
import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / '_tools' / 'derivations'))
HAS_SYMPY = importlib.util.find_spec('sympy') is not None


@unittest.skipUnless(HAS_SYMPY, 'sympy is not installed')
class StiffFluidUniverses(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import sympy as sp
        import verify_metrics as vm
        cls.sp = sp
        cls.metric = json.loads((ROOT / 'MFS/assets/data/metrics/stiff_fluid_universes.json').read_text())
        cls.charts = []
        for entry in cls.metric['coordinates']:
            r = vm.Reader(entry['coords'], [p['symbol'] for p in entry['parameters']], ())
            X = [r.symbol[c] for c in entry['coords']]
            a, b, ell = (r.parameters[c] for c in ('a', 'b', 'ell'))
            if entry['id'] == 'characteristic':
                t, x = (X[0]+X[1])/2, (X[1]-X[0])/2
            else:
                t = X[0]
                x = X[1] if entry['id'] == 'areal' else X[1]-b*t**2/(2*a)
            sigma = a*sp.log(t)+b*x
            phi = sigma/sp.sqrt(8*sp.pi)
            def matrix(rows, rank=2):
                out = sp.zeros(4)
                for row in rows:
                    i,j = (entry['coords'].index(c) for c in row['indices'])
                    out[i,j] = r(row['value'])
                return out
            g = matrix(entry['metric_components'])
            gi = matrix(entry['inverse_metric_components'])
            G = matrix(entry['einstein_tensor']['variants']['ll']['nonzero'])
            Ric = matrix(entry['ricci_tensor']['variants']['ll']['nonzero'])
            cls.charts.append((entry,r,X,a,b,ell,t,x,sigma,phi,g,gi,G,Ric))

    def close(self, mine, theirs, at, message=''):
        self.assertLess(abs(self.sp.N((mine-theirs).subs(at),35)),1e-26,message)

    def points(self, chart, av=1, bv=1, tv=None, xv=None, lv=1):
        entry,r,X,a,b,ell,t,x,*_ = chart
        sp=self.sp
        tv=sp.Rational(2,5) if tv is None else tv
        xv=sp.Rational(-1,4) if xv is None else xv
        at={a:av,b:bv,ell:lv,X[2]:sp.Rational(1,7),X[3]:sp.Rational(-2,7),r.c:1}
        if entry['id']=='characteristic':
            at.update({X[0]:tv-xv,X[1]:tv+xv})
        else:
            at.update({X[0]:tv,X[1]:xv if entry['id']=='areal' else xv+bv*tv**2/(2*av)})
        return at

    def test_every_chart_has_eight_pi_scalar_and_stiff_fluid_stress(self):
        sp=self.sp
        for chart in self.charts:
            entry,r,X,a,b,ell,t,x,sigma,phi,g,gi,G,Ric=chart
            d=sp.Matrix([sp.diff(phi,c) for c in X])
            for av,bv,tv,xv,lv in ((1,1,sp.Rational(2,5),sp.Rational(-1,4),1),
                                   (2,sp.Rational(1,2),1,sp.Rational(1,3),2),
                                   (sp.Rational(3,2),0,2,1,sp.Rational(3,2)),
                                   (1,1,sp.Rational(2,5),2,1)):
                at=self.points(chart,av,bv,tv,xv,lv)
                gn=g.subs(at).evalf(40);gin=gi.subs(at).evalf(40);dn=d.subs(at).evalf(40)
                norm=(dn.T*gin*dn)[0]
                self.assertLess(norm,0,entry['id'])
                rho=-norm/2
                # Future velocity is minus the raised gradient because a>0.
                U=-gin*dn/sp.sqrt(-norm)
                Ul=gn*U
                self.close((U.T*gn*U)[0],-1,{},entry['id'])
                for i in range(4):
                    for j in range(4):
                        canonical=dn[i]*dn[j]-gn[i,j]*norm/2
                        fluid=2*rho*Ul[i]*Ul[j]+rho*gn[i,j]
                        self.close(G[i,j],8*sp.pi*canonical,at,(entry['id'],i,j,'scalar'))
                        self.close(G[i,j],8*sp.pi*fluid,at,(entry['id'],i,j,'fluid'))
                        self.close(Ric[i,j],8*sp.pi*dn[i]*dn[j],at,(entry['id'],i,j,'Ricci'))

    def test_review_dictionary_and_scalar_wave_equation_in_every_chart(self):
        sp=self.sp
        self.assertIn(r'\sigma=\sqrt{8\pi}\Phi',self.metric['convention'])
        self.assertIn(r'G_{\mu\nu}=8\pi T_{\mu\nu}',self.metric['convention'])
        for chart in self.charts:
            entry,r,X,a,b,ell,t,x,sigma,phi,g,gi,G,Ric=chart
            d=sp.Matrix([sp.diff(phi,c) for c in X])
            ds=sp.Matrix([sp.diff(sigma,c) for c in X])
            # The determinant is known by the coordinate Jacobians, avoiding a symbolic determinant.
            F=t**(a*a-sp.Rational(1,2))*sp.exp(b*b*t*t/2+2*a*b*x)
            volume=ell**4*F*t/(2 if entry['id']=='characteristic' else 1)
            wave=sum(sp.diff(volume*sum(gi[i,j]*d[j] for j in range(4)),X[i]) for i in range(4))/volume
            for at in (self.points(chart),self.points(chart,2,sp.Rational(1,2),1,sp.Rational(1,3),2)):
                self.close(wave,0,at,entry['id'])
                for i in range(4):
                    for j in range(4):
                        self.close(ds[i]*ds[j],8*sp.pi*d[i]*d[j],at,entry['id'])

    def test_density_and_null_gradient_boundary(self):
        sp=self.sp
        for chart in self.charts:
            entry,r,X,a,b,ell,t,x,sigma,phi,g,gi,G,Ric=chart
            d=sp.Matrix([sp.diff(phi,c) for c in X])
            norm=(d.T*gi*d)[0]
            F=t**(a*a-sp.Rational(1,2))*sp.exp(b*b*t*t/2+2*a*b*x)
            expected=(a*a/t**2-b*b)/(16*sp.pi*ell**2*F)
            for av,bv,lv in ((1,1,1),(2,sp.Rational(3,2),2)):
                for factor in (sp.Rational(3,4),1,sp.Rational(5,4)):
                    at=self.points(chart,av,bv,factor*av/bv,sp.Rational(1,5),lv)
                    self.close(-norm/2,expected,at,entry['id'])
                    value=sp.N(expected.subs(at),35)
                    if factor<1:self.assertGreater(value,0)
                    elif factor==1:self.close(value,0,{},entry['id'])
                    else:self.assertLess(value,0)
                at=self.points(chart,av,bv,av/bv,0,lv)
                self.assertTrue(all(v.is_finite for v in G.subs(at).evalf()))
                self.assertGreater(max(abs(v) for v in G.subs(at).evalf()),0)

    def test_coordinate_changes_are_isometries_and_comoving_labels_follow_fluid(self):
        sp=self.sp
        base=next(c for c in self.charts if c[0]['id']=='areal')
        for chart in self.charts:
            entry,r,X,a,b,ell,t,x,sigma,phi,g,gi,G,Ric=chart
            old=sp.Matrix([t,x,X[2],X[3]])
            J=old.jacobian(X)
            F=t**(a*a-sp.Rational(1,2))*sp.exp(b*b*t*t/2+2*a*b*x)
            pull=J.T*sp.diag(-ell**2*F,ell**2*F,ell**2*t,ell**2*t)*J
            at=self.points(chart)
            for i in range(4):
                for j in range(4):self.close(g[i,j],pull[i,j],at,entry['id'])
            if entry['id']=='comoving':
                grad=gi*sp.Matrix([sp.diff(phi,c) for c in X])
                self.close(grad[1],0,at,'fluid labels must stay constant')

    def test_isotropic_member_has_three_one_third_powers(self):
        sp=self.sp
        # Differentiate the published areal spatial metric at b=0, using d tau=ell sqrt(F)dt.
        chart=next(c for c in self.charts if c[0]['id']=='areal')
        entry,r,X,a,b,ell,t,x,sigma,phi,g,gi,G,Ric=chart
        tau=ell*t**(a*a/2+sp.Rational(3,4))/(a*a/2+sp.Rational(3,4))
        powers=[]
        for i in range(1,4):
            spatial=g[i,i].subs(b,0)
            exponent=(sp.diff(spatial,t)/(2*spatial))/(sp.diff(tau,t)/tau)
            at=self.points(chart,sp.sqrt(sp.Rational(3,2)),0,2,1,2)
            self.close(exponent,sp.Rational(1,3),at,'isotropic proper-time power')
            powers.append(exponent)
        for av in (sp.Rational(1,2),1,2):
            self.close(sum(powers),1,self.points(chart,av,0,2,1,2),'Kasner sum')


if __name__=='__main__':
    unittest.main()
