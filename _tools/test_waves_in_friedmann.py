"""Independent checks of the published gravitational waves in Friedmann universes.

Every chart is held to Einstein's equation with a stiff fluid at rest, p = rho, the density its
convention states, the backgrounds to a vanishing Weyl tensor, and the plane wave to the flat
universe on its front z = eta, where the profile and its slope vanish.
"""
import importlib.util
import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / '_tools' / 'derivations'))
HAS_SYMPY = importlib.util.find_spec('sympy') is not None

DENSITY = {
    'plane_wave': '3*exp(S)/(64*pi*gamma*eta**3)',
    'flat_planar': '3/(64*pi*gamma*eta**3)', 'flat_cylindrical': '3/(64*pi*gamma*eta**3)',
    'open_planar': '3/(8*pi*gamma*sinh(2*eta)**3)', 'open_cylindrical': '3/(8*pi*gamma*sinh(2*eta)**3)',
    'closed_toroidal': '3/(8*pi*gamma*sin(2*eta)**3)',
}
POINTS = ((0.7, 0.2, 1.3), (1.1, -0.5, 0.6), (0.45, 0.3, 2.0))


@unittest.skipUnless(HAS_SYMPY, 'sympy is not installed')
class WavesInFriedmann(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import sympy as sp
        import verify_metrics as vm
        cls.sp = sp
        cls.metric = json.loads((ROOT / 'MFS/assets/data/metrics/waves_in_friedmann.json').read_text())
        cls.charts = {}
        for entry in cls.metric['coordinates']:
            r = vm.Reader(entry['coords'], [p['symbol'] for p in entry['parameters']], ())

            def matrix(rows, entry=entry, r=r):
                out = sp.zeros(4)
                for row in rows:
                    i, j = (entry['coords'].index(c) for c in row['indices'])
                    out[i, j] = r(row['value'])
                return out
            cls.charts[entry['id']] = (entry, r, matrix(entry['metric_components']),
                                       matrix(entry['einstein_tensor']['variants']['ll']['nonzero']),
                                       r(entry['ricci_scalar'].split('=', 1)[1]),
                                       entry['weyl_tensor']['variants']['llll']['nonzero'])

    def at(self, cid, eta, x, amplitude):
        entry, r, *_ = self.charts[cid]
        X = [r.symbol[c] for c in entry['coords']]
        point = {X[0]: eta, X[1]: x if cid != 'closed_toroidal' else 0.3 + x / 4, X[2]: 0.25, X[3]: -0.5,
                 r.parameters['gamma']: 1.7, r.c: 1}
        if 'A' in r.parameters:
            point[r.parameters['A']] = amplitude
        if cid == 'closed_toroidal':
            point[X[0]] = eta / 2
        return point

    def density(self, cid, point):
        sp = self.sp
        entry, r, *_ = self.charts[cid]
        names = {'eta': r.symbol['\\eta'], 'gamma': r.parameters['gamma'], 'pi': sp.pi}
        if cid == 'plane_wave':
            names['S'] = r('S')
        return sp.sympify(DENSITY[cid], locals=names).subs(point)

    def test_stiff_fluid_at_rest(self):
        """G_ab = 16 pi rho u_a u_b + 8 pi rho g_ab with u along eta, in every chart."""
        sp = self.sp
        for cid, (entry, r, g, G, R, _) in self.charts.items():
            for eta, x, amplitude in POINTS:
                point = self.at(cid, eta, x, amplitude)
                rho = self.density(cid, point)
                gn = g.subs(point)
                u = sp.zeros(4, 1)
                u[0] = -sp.sqrt(-gn[0, 0])
                expected = 16 * sp.pi * rho * u * u.T + 8 * sp.pi * rho * gn
                gap = max(abs(sp.N(v, 30)) for v in (G.subs(point) - expected))
                self.assertLess(gap, 1e-12, f'{cid} at {point}')
                self.assertLess(abs(sp.N(R.subs(point) + 16 * sp.pi * rho, 30)), 1e-12, f'{cid}: R = -16 pi rho')

    def test_backgrounds_conformally_flat_and_wave_not(self):
        for cid, (*_, weyl) in self.charts.items():
            if cid == 'plane_wave':
                self.assertTrue(weyl, 'the plane wave has a Weyl tensor')
            else:
                self.assertEqual(weyl, [], f'{cid} is conformally flat')

    def test_wave_meets_the_flat_universe_on_its_front(self):
        """On z = eta the profile H_2 and its slope vanish, so V = S = 0 and the plane wave's
        metric is the flat planar chart's."""
        sp = self.sp
        entry, r, g, *_ = self.charts['plane_wave']
        _, rf, gf, *_ = self.charts['flat_planar']
        eta, z = r.symbol['\\eta'], r.symbol['z']
        for value in (sp.Rational(1, 2), 1, 3):
            near = {eta: sp.Float(value, 50), z: sp.Float(value, 50) * (1 - sp.Float('1e-24', 50)), r.parameters['A']: 2,
                    r.parameters['gamma']: 1, r.c: 1}
            flat = {rf.symbol['\\eta']: value, rf.parameters['gamma']: 1, rf.c: 1}
            for i in range(4):
                gap = sp.N(g[i, i].subs(near) - gf[i, i].subs(flat), 40)
                self.assertLess(abs(gap), 1e-12, f'g_{i}{i} on the front at eta = {value}')

    def test_profile_normalised(self):
        """H_2(-1) = 1, as Alekseev and Griffiths normalise it."""
        sp = self.sp
        entry, r, *_ = self.charts['plane_wave']
        eta, z = r.symbol['\\eta'], r.symbol['z']
        H = r('H')
        self.assertLess(abs(sp.N(H.subs({eta: 1, z: sp.Float('-1', 50) + sp.Float('1e-30', 50)}), 40) - 1), 1e-12)


if __name__ == '__main__':
    unittest.main()
