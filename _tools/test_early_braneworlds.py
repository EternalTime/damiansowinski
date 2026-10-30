"""Localization and source-coordinate relations beyond the tensor generator's checks."""
import json
import math
import sys
import unittest
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.special import gamma, jv, jvp

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "_tools" / "derivations"))
import verify_metrics as vm


@lru_cache(None)
def chart(system, metric="early_braneworlds"):
    d = json.loads((ROOT / "MFS/assets/data/metrics" / (metric+".json")).read_text())
    c = next(c for c in d["coordinates"] if c["id"] == system)
    reader = vm.Reader(c["coords"], [p["symbol"] for p in c["parameters"]], ())
    values = {tuple(e["indices"]): reader(e["value"]) for e in c["metric_components"]}
    g = vm.sp.Matrix(len(c["coords"]), len(c["coords"]),
                     lambda i,j: values.get((c["coords"][i], c["coords"][j]), 0))
    return c, reader, g


def matrix(system, coordinates, parameters, metric="early_braneworlds"):
    c, reader, g = chart(system, metric)
    at = {reader.symbol[k]: v for k,v in zip(c["coords"], coordinates)}
    at.update({reader.parameters[k]: v for k,v in parameters.items()})
    return np.array(g.subs(at).evalf(30), dtype=float)


class Localization(unittest.TestCase):
    def test_visser_massive_geodesic_turns_at_its_energy_bound(self):
        c, reader, _ = chart("visser_electric")
        entry = next(e for e in c["christoffel"]["variants"]["ull"]["nonzero"]
                     if e["indices"] == ["w", "t", "t"])
        connection = reader(entry["value"])
        w, L = reader.symbol["w"], reader.parameters["L"]
        acceleration = vm.sp.lambdify(w, connection.subs(L,1), "numpy")
        w0, v0 = 0.2, 1.0
        energy = math.cosh(w0)*math.sqrt(1+v0*v0)
        def rhs(tau, state):
            return [state[1], -acceleration(state[0])*energy**2/math.cosh(state[0])**4]
        def turn(tau, state):
            return state[1]
        turn.terminal, turn.direction = True, -1
        solution = solve_ivp(rhs, (0,20), [w0,v0], events=turn, rtol=1e-10, atol=1e-12)
        self.assertEqual(len(solution.t_events[0]), 1)
        self.assertAlmostEqual(solution.y_events[0][0,0], math.acosh(energy), places=9)
        self.assertLess(energy**2/math.cosh(math.acosh(energy)+0.1)**2-1, 0)

    def test_scalar_wall_fermion_zero_mode_is_normalizable_for_one_chirality(self):
        # Rubakov (2001), (2)-(5): phi=v*tanh(k*w), f=sech(k*w)^(h*v/k).
        for exponent in (0.2,0.7,2.0):
            density = lambda w: math.exp(-2*exponent*(abs(w)+math.log1p(math.exp(-2*abs(w)))-math.log(2)))
            norm = 2*quad(density,0,np.inf,epsabs=1e-10)[0]
            self.assertAlmostEqual(norm,math.sqrt(math.pi)*gamma(exponent)/gamma(exponent+0.5),places=9)
            for w in (0.1,0.7,2):
                f = 1/math.cosh(w)**exponent
                derivative = -exponent*math.tanh(w)*f
                self.assertAlmostEqual(derivative+exponent*math.tanh(w)*f,0)
            self.assertGreater(1/density(20),1/density(10))

    def test_magnetic_flux_and_axisymmetric_scalar_continuum(self):
        # Source (3.7b), (5.2), (6.1c): flux pi*B0*a^2 and the Bessel radial equation.
        a, B0 = 1.3, 0.8
        flux = 2*math.pi*quad(lambda r: B0*r/(1+r*r/a**2)**2,0,np.inf)[0]
        self.assertAlmostEqual(flux,math.pi*B0*a*a,places=10)
        for r in (0.2,0.7,2,7):
            g = matrix("magnetic_flat", [0,0,0,0,r,0], {"a":a})
            volume = math.sqrt(-np.linalg.det(g))
            # Both kinetic and mass measures reduce to rho in the source's flat membrane.
            self.assertAlmostEqual(volume/g[4,4],r,places=10)
            self.assertAlmostEqual(volume/g[1,1],r,places=10)
            for mass in (1e-6,0.01,0.4,1.7):
                residual = -mass**2*jvp(0,mass*r,2)-mass*jvp(0,mass*r)/r-mass**2*jv(0,mass*r)
                self.assertAlmostEqual(residual,0,places=10)
        # The constant scalar has unbounded norm in the radial measure rho*d rho.
        self.assertEqual(quad(lambda r:r,0,20)[0]/quad(lambda r:r,0,10)[0],4)


class CoordinateRelations(unittest.TestCase):
    def test_flat_membrane_in_radial_warp_coordinates(self):
        for a in (0.7,1.3):
            params = {"C":3*a*a/2,"mu":4/(9*a*a),"Q":2/(3*a),"Lambda":0}
            for rho in (0.2,0.8,2.4):
                H = 1+rho*rho/a**2
                beta = H**(1/3)
                J = np.eye(6)
                J[4,4] = 2*rho/(3*a*a)*H**(-2/3)
                pulled = J.T@matrix("magnetic_general",[0,0,0,0,beta,0],params)@J
                direct = matrix("magnetic_flat",[0,0,0,0,rho,0],{"a":a})
                np.testing.assert_allclose(pulled,direct,rtol=2e-13,atol=2e-13)

    def test_squires_regular_axis_is_the_five_dimensional_ads_soliton(self):
        for L in (0.8,1.4):
            for w in (0.2,0.8,1.7):
                r = 2*L*math.sqrt(math.cosh(w/L))
                J = np.zeros((5,5))
                J[0,0], J[1,3], J[2,4], J[3,1], J[4,2] = 1, math.sinh(w/L)/math.sqrt(math.cosh(w/L)), L, 1, 1
                old = matrix("five_dimensional",[0,r,0,0,0],{"L":2*L,"r_0":2*L},"ads_soliton")
                direct = matrix("squires_axis",[0,0,0,w,0],{"L":L})
                np.testing.assert_allclose(J.T@old@J,direct,rtol=2e-13,atol=2e-13)


if __name__ == "__main__":
    unittest.main()
