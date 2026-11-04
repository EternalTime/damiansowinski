"""Physical claims of the averaged string, beyond the shared tensor checks."""
import json
from pathlib import Path
import sys
import unittest

import sympy as sp

sys.path.insert(0, str(Path(__file__).parent / "derivations"))
import verify_metrics as vm


class AveragedString(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        data = json.loads((Path(__file__).parents[1] / "MFS/assets/data/metrics/wiggly_string.json").read_text())
        cls.system = data["coordinates"][0]
        cls.reader = vm.Reader(cls.system["coords"], [p["symbol"] for p in cls.system["parameters"]], ())
        cls.g = vm.metric_from_line_element(cls.reader, cls.system["line_element"], cls.system["coords"])
        cls.t, cls.r, cls.phi, cls.z = [cls.reader.symbol[n] for n in cls.system["coords"]]
        cls.U, cls.T, cls.r0 = [cls.reader.parameters[n] for n in ("U", "T", "r_0")]
        cls.geo = vm.Geometry(cls.g, [cls.t, cls.r, cls.phi, cls.z], 30)

    def first_order(self, expression):
        e = sp.Symbol("epsilon")
        return sp.series(expression.subs({self.U:e*self.U, self.T:e*self.T}), e, 0, 2).removeO().subs(e, 1).expand()

    def test_rest_acceleration_and_straight_string(self):
        # At release dr/dt=dphi/dt=dz/dt=0, so d²r/dt²=-Gamma^r_tt.
        acceleration = -self.geo.christoffel_ull()[1][0][0]
        expected = -2*(self.U-self.T)/self.r
        self.assertEqual(sp.cancel(self.first_order(acceleration)-expected), 0)
        self.assertEqual(sp.cancel(acceleration.subs(self.T, self.U)), 0)
        potential = -(self.g[0, 0]+1)/2
        self.assertEqual(sp.cancel(-sp.diff(potential,self.r)-expected), 0)

    def test_equation_of_state_makes_the_force_inward(self):
        mu, alpha = sp.symbols("mu alpha", positive=True)
        energy, tension = alpha*mu, mu/alpha
        self.assertEqual(sp.cancel(energy*tension-mu**2), 0)
        acceleration = self.first_order(-self.geo.christoffel_ull()[1][0][0])
        reduced = sp.factor(acceleration.subs({self.U:energy, self.T:tension}))
        self.assertEqual(sp.cancel(reduced+2*mu*(alpha-1/alpha)/self.r), 0)
        for wiggliness in (sp.Rational(11,10), 2, 3):
            self.assertLess(reduced.subs({mu:sp.Rational(1,1000),alpha:wiggliness,self.r:2}), 0)
        self.assertEqual(reduced.subs(alpha,1), 0)

    def test_levi_civita_parameter_map_by_radial_pullback(self):
        # Camps-Emparan (2.1), (2.6), (2.7), in units r0=1.
        e = sp.Symbol("epsilon")
        r, U, T = self.r, self.U, self.T
        m, C = 2*e*(U-T), 1+4*e*T
        rho = r*(1+4*e*T*(1-sp.log(r)))
        weyl = [-sp.exp(2*m*sp.log(rho)),
                sp.exp(2*m*(m-1)*sp.log(rho))*sp.diff(rho,r)**2,
                sp.exp(2*(1-m)*sp.log(rho))/C**2,
                sp.exp(2*m*(m-1)*sp.log(rho))]
        for i, value in enumerate(weyl):
            linear = sp.series(value,e,0,2).removeO().expand()
            original = self.g[i,i].subs(self.r0,1).subs({U:e*U,T:e*T})
            self.assertEqual(sp.simplify(linear-original),0)
        sigma = m/2
        self.assertEqual(sp.expand(sigma-e*(U-T)),0)

    def test_exterior_ricci_begins_beyond_first_order(self):
        ricci = sp.Matrix(self.geo.ricci_ll())
        for value in ricci:
            self.assertEqual(sp.cancel(self.first_order(value)),0)
        # It is not an exact vacuum completion of the first-order line element.
        at = {self.r:2,self.r0:1,self.U:sp.Rational(1,100),self.T:sp.Rational(1,400)}
        self.assertTrue(any(abs(complex(sp.N(value.subs(at),20)))>1e-12 for value in ricci))


if __name__ == "__main__":
    unittest.main()
