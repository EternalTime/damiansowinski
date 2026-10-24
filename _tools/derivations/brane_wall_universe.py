#!/usr/bin/env python3
"""The universe on a wall moving past a black hole in five dimensions: the wall's motion, and the one
wall every diagram draws.

The wall is r = a(tau) in the black hole of anti-de Sitter space in five dimensions,
f = k + r^2/l^2 - mu/r^2, with the mirror symmetry of Randall and Sundrum's wall and each side kept
where r falls away from the wall. With rho the energy density of the matter on it, sigma =
3c^4/(4 pi G_5 l) its tension and a dot for d/d(c tau) along it, the junction condition is Kraus's (26),

    sqrt(f(a) + adot^2) = (a/l)(1 + rho/sigma),

so adot^2 = a^2 (1 + rho/sigma)^2/l^2 - f(a), and the static time runs as f(a) dt/d(c tau) =
(a/l)(1 + rho/sigma), the energy of the wall along the static Killing vector.

The wall drawn has no matter beside its tension, rho = 0, a flat space, k = 0, and mu = l^2, in units
of l: f = r^2 - 1/r^2, with the horizon at r_h = 1 and the surface gravity f'(1)/2 = 2. Then

    adot^2 = 1/a^2,  a^2 = 2 c tau,    the dark radiation of the black hole alone,
    dt/da = a^4/(a^4 - 1),             t(a) = a + r_*(a) - arctan(a),
    r_*(r) = (1/4) ln|(r - 1)/(r + 1)| + (1/2) arctan(r),   dr_*/dr = 1/f,

so on the wall the null times are u = t - r_* = a - arctan(a) and v = t + r_* = a + 2 r_*(a) - arctan(a).
Outside the horizon v runs from -infinity as a -> 1, where the wall comes out through the past
horizon at u = 1 - pi/4, and both run to infinity as a does: the wall reaches the conformal boundary
only at future timelike infinity. In the Gaussian normal chart Q = 1, and Mukohyama,
Shiromizu and Maeda's phi and psi are

    Phi = e^(-2w) + sinh^2(w)/(4 t^2),  Psi = e^(-2w) - sinh^2(w)/(4 t^2),

with t the wall's proper time and w the proper distance from it, so the lapse Psi/sqrt(Phi) vanishes
on w_*(t) = ln(1 + 4t)/2, where the normals from neighbouring moments meet.

_tools/test_brane_wall_universe.py holds these forms to the published metrics and to the junction
conditions, and brane_wall_universe.md beside this file is the derivation.
"""
import numpy as np
from scipy.optimize import brentq


# ---------------------------------------------------------------- any wall

def lapse_squared(r, ell, mu, k):
    """f = k + r^2/l^2 - mu/r^2, the black hole of anti-de Sitter space in five dimensions."""
    r = np.asarray(r, dtype=float)
    return k + r * r / ell ** 2 - mu / (r * r)


def speed_squared(a, ell, mu, k, rho_over_sigma=0.0):
    """(da/d(c tau))^2 = a^2 (1 + rho/sigma)^2/l^2 - f(a), the junction condition."""
    a = np.asarray(a, dtype=float)
    return a * a * (1 + rho_over_sigma) ** 2 / ell ** 2 - lapse_squared(a, ell, mu, k)


def static_rate(a, ell, mu, k, rho_over_sigma=0.0):
    """dt/d(c tau) on the wall, (a/l)(1 + rho/sigma)/f(a), outside the horizon."""
    a = np.asarray(a, dtype=float)
    return a * (1 + rho_over_sigma) / ell / lapse_squared(a, ell, mu, k)


# ---------------------------------------------------------------- the wall drawn, k = 0, rho = 0, mu = l^2 = 1

HORIZON = 1.0


def tortoise(r):
    """r_*, with dr_*/dr = 1/f = r^2/(r^4 - 1), zero where r -> 0 and pi/4 at the boundary."""
    r = np.asarray(r, dtype=float)
    return 0.25 * np.log(np.abs((r - 1) / (r + 1))) + 0.5 * np.arctan(r)


def static_time(a):
    """The static time at which the wall has the radius a > 1, t = a + r_*(a) - arctan(a)."""
    a = np.asarray(a, dtype=float)
    return a + tortoise(a) - np.arctan(a)


def radius_at_static_time(t):
    """The wall's radius at the static time t, the root a > 1 of static_time(a) = t."""
    def one(t_):
        # static_time rises from -infinity at a = 1 to infinity, so the root is bracketed.
        hi = max(2.0, t_ + 4.0)
        return brentq(lambda a: static_time(a) - t_, 1 + 1e-15, hi, xtol=1e-14, rtol=1e-14)
    t = np.asarray(t, dtype=float)
    return np.vectorize(one, otypes=[float])(t) if t.shape else one(float(t))


def proper_time(a):
    """c tau on the wall, a^2/2, zero at a = 0."""
    a = np.asarray(a, dtype=float)
    return a * a / 2


def scale_factor(tau):
    """a = sqrt(2 c tau), the wall's areal radius at its own proper time."""
    return np.sqrt(2 * np.asarray(tau, dtype=float))


def null_times(a):
    """(u, v) = (t - r_*, t + r_*) on the wall outside the horizon."""
    a = np.asarray(a, dtype=float)
    return a - np.arctan(a), a + 2 * tortoise(a) - np.arctan(a)


def phi(t, w):
    """Mukohyama, Shiromizu and Maeda's phi on the wall drawn, A^2 = a^2 phi."""
    t, w = np.asarray(t, dtype=float), np.asarray(w, dtype=float)
    return np.exp(-2 * w) + np.sinh(w) ** 2 / (4 * t * t)


def psi(t, w):
    """Their psi, the lapse times sqrt(phi)."""
    t, w = np.asarray(t, dtype=float), np.asarray(w, dtype=float)
    return np.exp(-2 * w) - np.sinh(w) ** 2 / (4 * t * t)


def focal_distance(t):
    """w_*(t) = ln(1 + 4t)/2, where psi and the lapse vanish."""
    return 0.5 * np.log1p(4 * np.asarray(t, dtype=float))
