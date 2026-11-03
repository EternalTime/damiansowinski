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


if __name__ == '__main__':
    f_h, solution, logn_far = solve()
    print(f'f_h={f_h:.12g}  mu_far={solution.y[2,-1]:.12g}  '
          f'B_out={baryon_number(f_h):.12g}  logN_far={logn_far:.12g}')
