#!/usr/bin/env python3
"""The quantum Oppenheimer-Snyder black hole as every diagram draws it, and the path of its dust.

Outside the dust the metric is Lewandowski, Ma, Yang and Zhang's, f = 1 - 2m/r + alpha m^2/r^4,
with m = GM/c^2 and alpha an area. Every diagram takes m = 1 and alpha = 5/4, so that

    r_b = (alpha/2)^(1/3) = 0.8550,   r_- = 1.1272,   r_+ = 1.7774,   kappa_-/kappa_+ = 3.34,

a hole with 1.16 times the least mass 4 sqrt(alpha)/3 sqrt(3) that has a horizon. Inside, the dust
is the flat Friedmann ball of loop quantum cosmology, a^3 = 1 + 9 c^2 tau^2/alpha, whose surface
chi_0 = r_b has the areal radius

    R(tau)^3 = r_b^3 + 9 m c^2 tau^2/2,

which turns round at tau = 0. It is the radial geodesic of the exterior that falls from rest at
infinity, E = f dt/dtau = 1 and (dR/d(c tau))^2 = 1 - f(R) = N(R)^2, so along it

    on the way in,  dR/dtau = -N,   dv/dtau = 1/(1 + N),    the Painleve-Gullstrand T = tau;
    on the way out, dR/dtau = +N,   dv/dtau = (1 + N)/f,    dT/dtau = (1 + N^2)/f,

with v = ct + r_* and T = ct + h(r), dr_*/dr = 1/f and dh/dr = N/f, both zero at r_b, and the
static time t zero at the turn, so the surface is at v = T = 0 there. The way out reaches r_- at a
finite proper time TAU_M, where v and T run off to infinity: it leaves both ingoing charts there,
through the horizon of the white hole. The outgoing chart is the time reverse of the ingoing one,
u = -v with tau -> -tau, and R is even in tau, so its surface is the ingoing chart's reflected,
radius_at_retarded(u) = radius_at_advanced(-u).

The tortoise coordinate is one partial fraction per zero of f, two real and two complex,
r_* = r + sum_i A_i ln(r - r_i) - r_*(r_b) with A_i = r_i^4/P'(r_i), P = r^4 - 2r^3 + alpha, the
complex pair taken as twice the real part of one term.

_tools/test_quantum_oppenheimer_snyder.py holds these forms to the published metrics and to the
junction, and quantum_oppenheimer_snyder.md beside this file is the derivation.
"""
import math

import numpy as np

ALPHA = 1.25                                # alpha in units of m^2, with m = 1
RB = (ALPHA / 2) ** (1 / 3)                 # where the surface turns round, f(r_b) = 1
ROOTS = np.roots([1.0, -2.0, 0.0, 0.0, ALPHA])
REAL = sorted(float(z.real) for z in ROOTS if abs(z.imag) < 1e-12)
RM, RP = REAL                               # the inner and the outer horizon
COMPLEX = next(complex(z) for z in ROOTS if z.imag > 1e-12)
M_LEAST = 4 * math.sqrt(ALPHA) / (3 * math.sqrt(3))     # the least mass with a horizon, in m


def f(r):
    r = np.asarray(r, dtype=float)
    return 1 - 2 / r + ALPHA / r ** 4


def speed(r):
    """N = sqrt(1 - f) = sqrt(2m/r - alpha m^2/r^4), the speed at which the surface falls, real for r >= r_b."""
    r = np.asarray(r, dtype=float)
    return np.sqrt(np.maximum(2 / r - ALPHA / r ** 4, 0.0))


def kappa(r):
    """|f'(r)|/2, the surface gravity of a horizon at r."""
    return abs(2 / r ** 2 - 4 * ALPHA / r ** 5) / 2


def _residue(z):
    return z ** 4 / (4 * z ** 3 - 6 * z ** 2)


def _tortoise_raw(r):
    r = np.asarray(r, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        out = r + sum(_residue(z) * np.log(np.abs(r - z)) for z in REAL)
        return out + 2 * np.real(_residue(COMPLEX) * np.log(r - COMPLEX))


_TORTOISE_RB = float(_tortoise_raw(RB))


def tortoise(r):
    """r_* with dr_*/dr = 1/f and r_* = 0 at r_b."""
    return _tortoise_raw(r) - _TORTOISE_RB


def radius(tau):
    """R(tau), the areal radius of the surface at the dust's proper time c tau, in units of m."""
    tau = np.asarray(tau, dtype=float)
    return np.cbrt(RB ** 3 + 4.5 * tau ** 2)


def scale_factor(tau):
    """a(tau) = (1 + 9 c^2 tau^2/alpha)^(1/3), so that R = a r_b."""
    return np.cbrt(1 + 9 * np.asarray(tau, dtype=float) ** 2 / ALPHA)


# The proper time at which the way out reaches r_-: R(TAU_M) = r_-.
TAU_M = math.sqrt((RM ** 3 - RB ** 3) / 4.5)
TAU_P = math.sqrt((RP ** 3 - RB ** 3) / 4.5)        # and r_+, on the way in at -TAU_P


# Near TAU_M, f(R) = f'(r_-) N(r_-) (tau - TAU_M) with N(r_-) = 1, so both rates on the way out,
# (1 + N)/f and (1 + N^2)/f, have the simple pole 2/(f'(r_-)(tau - TAU_M)) = POLE/(TAU_M - tau).
POLE = 1 / kappa(RM)


def _in_rate(t):
    return 1 / (1 + speed(radius(t)))


def _out_v_rate(t):
    R = radius(t)
    return (1 + speed(R)) / f(R)


def _out_T_rate(t):
    R = radius(t)
    return (1 + speed(R) ** 2) / f(R)


def _cumulative(rate, grid):
    """The integral of rate from grid[0] along grid, by the trapezium rule on a grid fine enough
    that its error is below 1e-9 in the units of m."""
    y = rate(grid)
    return np.concatenate([[0.0], np.cumsum((y[1:] + y[:-1]) * np.diff(grid) / 2)])


# The way in, from far in the past to the turn, and the way out, in s = TAU_M - tau with the pole's
# part POLE ln(TAU_M/s) taken in closed form and the smooth rest summed.
_TAU_IN = -np.concatenate([np.linspace(0, 1, 40001)[1:], np.geomspace(1, 1e4, 40001)[1:]])[::-1]
_TAU_IN = np.concatenate([_TAU_IN, [0.0]])
_V_IN = _cumulative(_in_rate, _TAU_IN)
_V_IN -= _V_IN[-1]
_SS = np.concatenate([TAU_M * np.geomspace(1e-12, 1e-3, 20001)[:-1], TAU_M * np.linspace(1e-3, 1, 60001)])[::-1]
_TAU_OUT = TAU_M - _SS


def _way_out(rate):
    smooth = _cumulative(lambda t: rate(t) - POLE / (TAU_M - t), _TAU_OUT)
    return smooth + POLE * np.log(TAU_M / _SS)


_V_OUT, _T_OUT = _way_out(_out_v_rate), _way_out(_out_T_rate)
_TAU_V, _V = np.concatenate([_TAU_IN, _TAU_OUT[1:]]), np.concatenate([_V_IN, _V_OUT[1:]])


def advanced(tau):
    """v of the surface at its proper time tau, for tau < TAU_M."""
    return np.interp(np.asarray(tau, dtype=float), _TAU_V, _V)


def slice_time(tau):
    """The Painleve-Gullstrand T of the surface, tau on the way in."""
    tau = np.asarray(tau, dtype=float)
    return np.where(tau <= 0, tau, np.interp(tau, _TAU_OUT, _T_OUT))


def _invert(times, values, x):
    """The time at which a table reaches x; beyond its last value the surface is at r_- to within
    the table's last step, 1e-12 TAU_M."""
    return np.interp(np.asarray(x, dtype=float), values, times)


def radius_at_advanced(v):
    """The radius of the surface at the advanced time v of the ingoing chart."""
    v = np.asarray(v, dtype=float)
    far = v < _V[0]
    out = radius(_invert(_TAU_V, _V, v))
    # Far in the past the surface falls almost freely, v = tau + O(tau^(1/3)) and so R = R(v) to that order.
    return np.where(far, radius(v), out)


def radius_at_retarded(u):
    """The radius of the surface at the retarded time u of the outgoing chart, the time reverse."""
    return radius_at_advanced(-np.asarray(u, dtype=float))


def radius_at_slice_time(T):
    """The radius of the surface at the Painleve-Gullstrand time T of the ingoing chart."""
    T = np.asarray(T, dtype=float)
    return np.where(T <= 0, radius(T), radius(_invert(_TAU_OUT, _T_OUT, np.maximum(T, 0.0))))


def lead(r):
    """v - T on a sphere of radius r, the integral of 1/(1 + N) from r_b, which is smooth across both
    horizons: the Painleve-Gullstrand T = ct + h(r) with h = r_* - lead."""
    from scipy.integrate import quad
    out = [quad(lambda x: 1 / (1 + float(speed(x))), RB, float(x), epsabs=1e-13, epsrel=1e-12, limit=200)[0]
           for x in np.atleast_1d(np.asarray(r, dtype=float))]
    return np.array(out).reshape(np.shape(r))


def conformal_time(tau):
    """eta = int_0^tau c dtau/a, the conformal time of the dust, in closed form:
    tau 2F1(1/3, 1/2; 3/2; -9 tau^2/alpha), odd in tau and zero at the bounce."""
    from scipy.special import hyp2f1
    tau = np.asarray(tau, dtype=float)
    return tau * hyp2f1(1 / 3, 0.5, 1.5, -9 * tau ** 2 / ALPHA)
