#!/usr/bin/env python3
"""Numerical exterior of a spherical Einstein-Skyrme black hole.

Shiiki and Sawado, arXiv:gr-qc/0501025, equations (21)-(35), use
x = a F_pi r, mu = a F_pi G m and alpha = 4 pi G F_pi^2.
Here C = 1 - 2 mu/x and a prime means d/dx. The massless-pion field is
integrated outward from a regular, nondegenerate horizon; the shooting
condition x f' + 2 f = 0 at a distant boundary selects the decaying tail.
"""
import math
from functools import lru_cache

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

ALPHA = 0.01
HORIZON = 0.1
FAR = 60.0
EPSILON = 1e-6


def horizon_data(f_h, alpha=ALPHA, x_h=HORIZON):
    """Horizon derivatives from equations (33) and (34)."""
    s = math.sin(f_h)
    mu1 = alpha / 4 * (s**2 + 2 * s**4 / x_h**2)
    f1 = ((x_h**2 + 4 * s**2) * math.sin(2 * f_h)
          / (x_h * (x_h**2 + 8 * s**2) * (1 - 2 * mu1)))
    return f1, mu1


def equations(x, y, alpha=ALPHA):
    """f, f', mu and log(N); divergence form of the static field equation."""
    f, p, mu, _ = y
    s = math.sin(f)
    u = x*x + 8*s*s
    c = 1 - 2*mu/x
    mu_p = alpha/8 * (c*u*p*p + 2*s*s + 4*s**4/x**2)
    c_p = 2*mu/x**2 - 2*mu_p/x
    logn_p = alpha*u*p*p/(4*x)
    u_p = 2*x + 8*math.sin(2*f)*p
    p_p = (math.sin(2*f)*(1 + 4*c*p*p + 4*s*s/x**2)
           - (logn_p*c*u + c_p*u + c*u_p)*p)/(c*u)
    return [p, p_p, mu_p, logn_p]


def integrate(f_h, far=FAR, epsilon=EPSILON, rtol=2e-10, atol=1e-12):
    f1, mu1 = horizon_data(f_h)
    if 2*mu1 >= 1:
        raise ValueError('The selected horizon is not nondegenerate and outer')
    y0 = [f_h + epsilon*f1, f1, HORIZON/2 + epsilon*mu1, 0.0]
    result = solve_ivp(equations, (HORIZON + epsilon, far), y0,
                       method='DOP853', rtol=rtol, atol=atol, dense_output=True)
    if not result.success:
        raise RuntimeError(result.message)
    return result


@lru_cache(maxsize=4)
def solve(far=FAR, epsilon=EPSILON):
    """Return the large-f_h branch and a normalized lapse at the far boundary."""
    def residual(f_h):
        state = integrate(f_h, far, epsilon).y[:, -1]
        return far*state[1] + 2*state[0]

    f_h = brentq(residual, 3.02, 3.06, xtol=1e-13)
    result = integrate(f_h, far, epsilon)
    return f_h, result, float(result.y[3, -1])


def baryon_number(f_h):
    return (2*f_h - math.sin(2*f_h))/(2*math.pi)


def jet(x):
    """N, C and their first two derivatives in the dimensionless radius."""
    x = np.asarray(x, dtype=float)
    shape = x.shape
    values = np.atleast_1d(x).ravel()
    result = {key: np.full(values.shape, np.nan) for key in ('N', 'dN', 'd2N', 'C', 'dC', 'd2C')}
    valid = (values >= HORIZON) & (values <= FAR)
    if not np.any(valid):
        return {key: value.reshape(shape) for key, value in result.items()}
    xx = np.maximum(values[valid], HORIZON + EPSILON)
    _, sol, logn_far = solve()
    f, p, mu, logn = sol.sol(xx)
    sinf, cosf = np.sin(f), np.cos(f)
    u = xx**2 + 8*sinf**2
    c = 1 - 2*mu/xx
    mu_p = ALPHA/8*(c*u*p**2 + 2*sinf**2 + 4*sinf**4/xx**2)
    c_p = 2*mu/xx**2 - 2*mu_p/xx
    q = ALPHA*u*p**2/(4*xx)
    u_p = 2*xx + 8*np.sin(2*f)*p
    p_p = (np.sin(2*f)*(1 + 4*c*p**2 + 4*sinf**2/xx**2)
           - (q*c*u + c_p*u + c*u_p)*p)/(c*u)
    mu_pp = ALPHA/8*(c_p*u*p**2 + c*u_p*p**2 + 2*c*u*p*p_p
                       + 2*np.sin(2*f)*p + 16*sinf**3*cosf*p/xx**2 - 8*sinf**4/xx**3)
    c_pp = -2*mu_pp/xx + 4*mu_p/xx**2 - 4*mu/xx**3
    q_p = ALPHA/4*((u_p/xx - u/xx**2)*p**2 + 2*u*p*p_p/xx)
    n = np.exp(logn - logn_far)
    data = {'N': n, 'dN': n*q, 'd2N': n*(q*q + q_p),
            'C': c, 'dC': c_p, 'd2C': c_pp}
    # The horizon itself is the exact zero; the nearby lapse uses its limiting value.
    data['C'] = np.where(values[valid] == HORIZON, 0.0, data['C'])
    for key in result:
        result[key][valid] = data[key]
    return {key: value.reshape(shape) for key, value in result.items()}


def declared_functions():
    """SymPy functions with the solved first and second derivatives for drawing."""
    import sympy as sp

    def make(name, key, next_function=None):
        def value(x):
            return jet(x)[key]

        def fdiff(self, argindex=1):
            if argindex != 1 or next_function is None:
                raise sp.ArgumentIndexError(self, argindex)
            return next_function(self.args[0])

        return type(name, (sp.Function,), {'nargs': 1, 'is_real': True,
                                           '_imp_': staticmethod(value), 'fdiff': fdiff,
                                           '_eval_evalf': lambda self, prec: sp.Float(
                                               float(self._imp_(float(self.args[0]))), prec)})

    result = {}
    for key in ('N', 'C'):
        second = make('skyrme_d2' + key, 'd2' + key)
        first = make('skyrme_d' + key, 'd' + key, second)
        base = make('skyrme_' + key, key, first)
        result.update({fn.__name__: fn for fn in (base, first, second)})
    return result


if __name__ == '__main__':
    f_h, solution, logn_far = solve()
    print(f'f_h={f_h:.12g}  mu_far={solution.y[2,-1]:.12g}  '
          f'B_out={baryon_number(f_h):.12g}  logN_far={logn_far:.12g}')
