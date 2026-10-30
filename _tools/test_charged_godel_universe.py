"""The charged rotating universes, held to the stored charts and drawings.

Pin the Einstein-Maxwell source, physical parameter relations, the two published
limits, circular chronology threshold and Euclidean embedding obstruction.
Run only this module with the collection's Python environment.
"""
import ast
import importlib.util
import json
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / 'MFS' / 'assets' / 'data'
DERIVATIONS = ROOT / '_tools' / 'derivations'
HAS_SYMPY = importlib.util.find_spec('sympy') is not None
sys.path.insert(0, str(DERIVATIONS))
SYSTEMS = ('hyperbolic', 'linear', 'trigonometric')


def data(kind, name='charged_godel_universe'):
    return json.loads((DATA / kind / (name + '.json')).read_text())


def chart(name, system):
    import sympy as sp
    import verify_metrics as vm
    entry = next(c for c in data('metrics', name)['coordinates'] if c['id'] == system)
    reader = vm.Reader(entry['coords'], [p['symbol'] for p in entry['parameters']], ())

    def matrix(items):
        out = sp.zeros(4)
        for item in items:
            i, j = (entry['coords'].index(c) for c in item['indices'])
            out[i, j] = reader(item['value']).subs(reader.c, 1)
        return out
    return entry, reader, matrix


@unittest.skipUnless(HAS_SYMPY, 'sympy is not installed')
class PublishedPhysics(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import sympy as sp
        cls.sp = sp
        cls.charts = {}
        for system in SYSTEMS:
            entry, reader, matrix = chart('charged_godel_universe', system)
            cls.charts[system] = (reader, matrix(entry['metric_components']),
                                  matrix(entry['inverse_metric_components']),
                                  matrix(entry['einstein_tensor']['variants']['ll']['nonzero']))

    def close_matrix(self, actual, expected, point, label):
        difference = (actual - expected).subs(point).evalf(35)
        self.assertLess(max(abs(complex(v)) for v in difference), 1e-10, label)

    def parameters(self, system, reader, w, q, r):
        sp = self.sp
        values = {reader.parameters['omega']: sp.Rational(w), reader.symbol['r']: sp.Rational(r)}
        if system != 'linear':
            values[reader.parameters['m' if system == 'hyperbolic' else 'k']] = sp.Rational(q)
        return values

    def test_each_stored_chart_satisfies_einstein_maxwell_with_force_free_charged_dust(self):
        sp = self.sp
        for system, (reader, g, gi, einstein) in self.charts.items():
            r = reader.symbol['r']
            w = reader.parameters['omega']
            q = reader.parameters.get('m', reader.parameters.get('k', sp.S.Zero))
            m2 = q*q if system == 'hyperbolic' else -q*q if system == 'trigonometric' else sp.S.Zero
            D = sp.sinh(q*r)/q if system == 'hyperbolic' else sp.sin(q*r)/q if system == 'trigonometric' else r
            B = sp.sqrt(w*w-m2/2)
            field = sp.zeros(4)
            field[1, 2], field[2, 1] = B*D, -B*D
            raised = gi*field*gi
            coordinates = [reader.symbol[c] for c in ('t', 'r', '\\phi', 'z')]
            current = sp.Matrix([sum(sp.diff(D*raised[i, j], coordinates[j]) for j in range(4))/D for i in range(4)])
            for wv, qv, rv in ((1, 1, '3/10'), ('3/2', '1/2', '4/5'), (2, '3/2', '11/10')):
                at = self.parameters(system, reader, wv, qv, rv)
                gn, inv, fn = (a.subs(at).evalf(35) for a in (g, gi, field))
                square = sum(fn[i,j]*(inv*fn*inv)[i,j] for i in range(4) for j in range(4))
                electromagnetic = 2*(fn*inv*fn.T-gn*square/4)  # 8 pi T_EM in Gaussian units.
                covelocity = -gn[:, 0]
                dust = (2*w*w).subs(at)*covelocity*covelocity.T
                self.close_matrix(einstein.subs(at)-m2.subs(at)*gn/2,
                                  dust+electromagnetic, {}, system)
                self.close_matrix(current, sp.Matrix([2*w*B,0,0,0]), at, system+' current')
                self.close_matrix(field*current, sp.zeros(4,1), at, system+' Lorentz force')
                self.assertAlmostEqual(float((g[2,2]+g[0,2]**2-D**2).subs(at)), 0, places=10)

    def test_parameter_relations_are_inferred_from_the_stored_einstein_tensor(self):
        sp = self.sp
        for system, B2, cosmological in (('hyperbolic', 0.5, -0.5), ('linear', 1, 0), ('trigonometric', 1.5, 0.5)):
            reader, g, _, einstein = self.charts[system]
            at = self.parameters(system, reader, 1, 1, '2/5')
            gn, En = g.subs(at).evalf(35), einstein.subs(at).evalf(35)
            H, D = -gn[0,2], sp.sqrt(gn[2,2]+gn[0,2]**2)
            inverse_frame = sp.Matrix([[1,0,-H/D,0],[0,1,0,0],[0,0,1/D,0],[0,0,0,1]])
            frame = inverse_frame.T*En*inverse_frame
            inferred_lambda = -(frame[1,1]+frame[3,3])/2
            inferred_B2 = (frame[1,1]-frame[3,3])/2
            density = (frame[0,0]-inferred_lambda-inferred_B2)/(8*sp.pi)
            self.assertAlmostEqual(float(inferred_B2), B2, places=12)
            self.assertAlmostEqual(float(inferred_lambda), cosmological, places=12)
            self.assertAlmostEqual(float(density), 1/(4*math.pi), places=12)

    def test_zero_magnetic_field_is_the_stored_godel_metric(self):
        sp = self.sp
        reader, g, _, einstein = self.charts['hyperbolic']
        old, R, matrix = chart('godel', 'cylindrical')
        theirs = matrix(old['metric_components'])
        jacobian = sp.diag(1/sp.sqrt(2),1/sp.sqrt(2),1,1/sp.sqrt(2))
        for rv in ('1/5','3/5','6/5'):
            at = self.parameters('hyperbolic', reader, 1, 1, rv)
            at[reader.parameters['m']] = sp.sqrt(2)  # B^2=1-m^2/2=0.
            other = {R.parameters['omega']:1, R.symbol['r']:sp.Rational(rv)/sp.sqrt(2)}
            self.close_matrix(g.subs(at), jacobian.T*theirs.subs(other)*jacobian, {}, 'Godel endpoint')
            self.assertAlmostEqual(float((-reader.parameters['m']**2/2).subs(at)), -1, places=12)
            self.assertAlmostEqual(float(einstein[0,0].subs(at)), 1, places=12)

    def test_m_zero_limit_has_som_raychaudhuri_cylindrical_functions_and_metric(self):
        sp = self.sp
        R, g, _, _ = self.charts['hyperbolic']
        r, w, m = R.symbol['r'], R.parameters['omega'], R.parameters['m']
        H = -g[0,2]
        D2 = g[2,2]+H*H
        self.assertEqual(sp.limit(H,m,0), w*r*r)
        self.assertEqual(sp.limit(D2,m,0), r*r)
        entry, SR, matrix = chart('som_raychaudhuri','cylindrical')
        theirs = matrix(entry['metric_components'])
        for rv in ('1/5','3/5','6/5'):
            at = {r:sp.Rational(rv),w:1,m:sp.Rational(1,10**6)}
            other = {SR.symbol['r']:sp.Rational(rv),SR.parameters['Omega']:1}
            self.close_matrix(g.subs(at),theirs.subs(other),{},'Som-Raychaudhuri limit')
            self.assertAlmostEqual(float(sp.sqrt(D2.subs(at)).evalf(35)),float(rv.split('/')[0])/float(rv.split('/')[1]),places=10)

    def test_the_first_null_radius_is_a_zero_of_the_stored_circular_tangent(self):
        for system, (R,g,_,_) in self.charts.items():
            for w,q in ((1,1),(1.5,0.5),(1.2,1.4)):
                rc = 2*math.atanh(q/(2*w))/q if system=='hyperbolic' else 2*math.atan(q/(2*w))/q if system=='trigonometric' else 1/w
                at = self.parameters(system,R,str(w),str(q),str(rc))
                self.assertAlmostEqual(float(g[2,2].subs(at).evalf(35)),0,places=10)
                self.assertGreater(float(g[2,2].subs({**at,R.symbol['r']:.99*rc})),0)
                self.assertLess(float(g[2,2].subs({**at,R.symbol['r']:1.01*rc})),0)
        self.assertAlmostEqual(2*math.atanh(.5),math.log(3),places=14)

    def test_the_hyperbolic_embedding_ends_at_its_euclidean_obstruction(self):
        sp = self.sp
        R,g,_,_ = self.charts['hyperbolic']
        view = next(v for v in data('embedding')['views'] if v['id']=='hyperbolic')
        stop = view['surfaces'][0]['pieces'][0]['points'][-1][0]
        self.assertAlmostEqual(stop,0.9516430642375219,places=12)
        self.assertLess(stop,math.log(3))
        r = R.symbol['r']
        gpp = g[2,2].subs({R.parameters['omega']:1,R.parameters['m']:1})
        defect = 1-sp.diff(gpp,r)**2/(4*gpp)
        self.assertAlmostEqual(float(defect.subs(r,stop).evalf(35)),0,places=10)
        self.assertGreater(float(defect.subs(r,stop-.001)),0)
        self.assertLess(float(defect.subs(r,stop+.001)),0)
        for v in data('embedding')['views']:
            self.assertTrue(all(math.isfinite(x) for s in v['surfaces'] for p in s['pieces'] for point in p['points'] for x in point),v['id'])


class Explanations(unittest.TestCase):
    def test_conformal_omission_has_its_geometric_reason_in_the_derivation(self):
        tree = ast.parse((DERIVATIONS/'conformal.py').read_text())
        declaration = next(n.value for n in tree.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='NOT_DRAWN' for t in n.targets))
        self.assertIn('charged_godel_universe',ast.literal_eval(declaration))
        self.assertIn('a global time-oriented radial conformal square does not encode its closed angular timelike curves',(DERIVATIONS/'charged_godel_universe.md').read_text())

    def test_slice_exceptions_explain_the_geometric_carve_outs(self):
        exceptions = json.loads((DERIVATIONS/'exceptions'/'charged_godel_universe.json').read_text())
        self.assertIn('beyond the first null circle',exceptions['hidden'][0]['reason'])
        self.assertIn('different curvature parameters',exceptions['hidden_views'][0]['reason'])
        for system in SYSTEMS:
            self.assertIn('charged_godel_universe/'+system+'/beyond',exceptions['hidden'][0]['drawings'])
            self.assertEqual(set(exceptions['hidden_views'][0]['drawings']['charged_godel_universe/'+system+'/inside']),set(SYSTEMS)-{system})


if __name__ == '__main__':
    unittest.main()
