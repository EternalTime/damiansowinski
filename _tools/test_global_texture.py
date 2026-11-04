"""Physical claims of the collapsing texture, derived from its published metric."""
import json
import sys
import unittest
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / '_tools' / 'derivations'))
import verify_metrics as vm


class GlobalTexture(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.metric = json.loads((ROOT / 'MFS/assets/data/metrics/global_texture.json').read_text())
        cls.chart = cls.metric['coordinates'][0]
        cls.reader = vm.Reader(cls.chart['coords'], [p['symbol'] for p in cls.chart['parameters']], ())
        cls.g = vm.metric_from_line_element(cls.reader, cls.chart['line_element'], cls.chart['coords'])
        cls.t, cls.r, cls.theta, _ = [cls.reader.symbol[c] for c in cls.chart['coords']]
        cls.epsilon = cls.reader.parameters['epsilon']

    def test_missing_solid_angle_from_sphere_area(self):
        # g_theta_theta and g_phi_phi/sin(theta)^2 fix the sphere area.
        area_radius_squared = self.g[2, 2]
        self.assertEqual(sp.simplify(self.g[3, 3] - area_radius_squared*sp.sin(self.theta)**2), 0)
        missing = 4*sp.pi*(1-area_radius_squared/self.r**2)
        x = sp.Symbol('x', positive=True)
        result = sp.simplify(missing.subs(self.t, -x*self.r))
        self.assertEqual(sp.simplify(result-8*sp.pi*self.epsilon*(1-x*sp.atan(1/x))), 0)
        self.assertEqual(sp.limit(result, x, sp.oo), 0)
        self.assertEqual(sp.limit(result, x, 0, dir='+'), 8*sp.pi*self.epsilon)

    def test_deficit_grows_with_time_at_fixed_proper_radius(self):
        missing = 4*sp.pi*(1-self.g[2, 2]/self.r**2)
        growth = sp.diff(missing, self.t)
        x = sp.Symbol('x', positive=True)
        expected = 8*sp.pi*self.epsilon/self.r*(sp.atan(1/x)-x/(1+x**2))
        self.assertEqual(sp.simplify(growth.subs(self.t, -x*self.r)-expected), 0)
        # atan(z)-z/(1+z^2) starts at zero and has a positive derivative for z>0.
        z = sp.Symbol('z', positive=True)
        positive_factor = sp.atan(z)-z/(1+z**2)
        self.assertEqual(positive_factor.subs(z, 0), 0)
        self.assertEqual(sp.simplify(sp.diff(positive_factor, z)-2*z**2/(1+z**2)**2), 0)
        for time in [-100, -10, -1, -.01]:
            self.assertGreater(float(growth.subs({self.t:time,self.r:2,self.epsilon:.05})), 0)

    def test_source_angular_factors_agree_to_first_order(self):
        # Durrer et al. (0.46), squared only through first order, against g_theta_theta/r^2.
        omega = 1+self.epsilon*(self.t/self.r*sp.atan(self.r/self.t)-1)
        linear_square = 1+2*(omega-1)
        self.assertEqual(sp.simplify(self.g[2, 2]/self.r**2-linear_square), 0)
        self.assertEqual(sp.simplify(omega**2-linear_square-(omega-1)**2), 0)

    def test_radial_geodesics_are_exactly_straight(self):
        # Constant metric components with no t-r cross term give zero radial connections.
        self.assertEqual(self.g[0, 0], -1)
        self.assertEqual(self.g[1, 1], 1)
        self.assertEqual(self.g[0, 1], 0)
        coords = [self.reader.symbol[c] for c in self.chart['coords']]
        for a in (0, 1):
            for b in (0, 1):
                for c in (0, 1):
                    connection = (1 if a==1 else -1)*(
                        sp.diff(self.g[a,c],coords[b])+sp.diff(self.g[a,b],coords[c])
                        -sp.diff(self.g[b,c],coords[a]))/2
                    self.assertEqual(connection, 0)
        self.assertFalse(any(set(e['indices']) <= {'t','r'}
                             for e in self.chart['christoffel']['variants']['ull']['nonzero']))

    def test_dust_density_follows_the_spatial_volume(self):
        # g_rr=1 and the sphere metric fix sqrt(det h) = r^2 H sin(theta).
        volume_squared = self.g[1,1]*self.g[2,2]*self.g[3,3]
        H = self.g[2,2]/self.r**2
        self.assertEqual(sp.simplify(volume_squared-(self.r**2*H*sp.sin(self.theta))**2), 0)
        density_ratio = 1/H
        for radius in [.1,1,4]:
            earlier = float(density_ratio.subs({self.t:-10,self.r:radius,self.epsilon:.05}))
            later = float(density_ratio.subs({self.t:-.1,self.r:radius,self.epsilon:.05}))
            self.assertGreater(later, earlier)


if __name__ == '__main__':
    unittest.main()
