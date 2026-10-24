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


# ---------------------------------------------------------------- one moment of the wall's own time

KAPPA = 2.0                                     # the surface gravity f'(1)/2 of the drawn hole


def kruskal_product(r):
    """Kruskal's UV = (1 - r) e^(2 arctan r)/(1 + r) at r, with U = -e^(-2u) and V = e^(2v) outside."""
    r = np.asarray(r, dtype=float)
    return (1 - r) / (1 + r) * np.exp(2 * np.arctan(r))


class Moment:
    """The surface of constant t of the Gaussian normal chart, the wall's own moment at its proper time
    t, followed out from the wall along the normal geodesics, which keep E = da/d(c tau) = 1/a along the
    static Killing vector. Mukohyama, Shiromizu and Maeda's (7) at k = 0 and l = mu = 1 is

        2 r^2 + E^2 = sqrt(4 + E^4) cosh(2 (w_0 - w)),   cosh(2 w_0) = (2 a^2 + E^2)/sqrt(4 + E^4),

    so r falls from the wall to r_min^2 = (sqrt(4 + E^4) - E^2)/2 < 1 at w = w_0, inside the black hole,
    and grows again into the other exterior. With dr/dw = -+sqrt(f + E^2) before and after w_0 and
    dt/dw = E/f, the null times change as dv/dw = -1/(E + sqrt(f + E^2)) before w_0 and
    du/dw = -1/(E + sqrt(f + E^2)) after it, both smooth through the horizons, so V = e^(2v) is
    followed from the wall to w_0 and U from w_0 on, the other from UV = kruskal_product(r)."""

    def __init__(self, t):
        self.t = float(t)
        self.a = float(scale_factor(self.t))
        self.E = 1 / self.a
        root = np.sqrt(4 + self.E ** 4)
        self.root = root
        self.w0 = 0.5 * np.arccosh((2 * self.a ** 2 + self.E ** 2) / root)
        self.r_min = np.sqrt((root - self.E ** 2) / 2)
        u_w, v_w = null_times(self.a) if self.a > 1 else (None, None)
        self.v_wall = v_w

    def radius(self, w):
        w = np.asarray(w, dtype=float)
        return np.sqrt((self.root * np.cosh(2 * (self.w0 - w)) - self.E ** 2) / 2)

    def rate(self, w):
        """1/(E + sqrt(f + E^2)), with f + E^2 = (dr/dw)^2."""
        r = self.radius(w)
        return 1 / (self.E + np.sqrt(np.maximum(lapse_squared(r, 1, 1, 0) + self.E ** 2, 0)))

    def kruskal(self, w):
        """(U, V) at the proper distances w from the wall, for a wall outside the horizon."""
        from scipy.integrate import quad
        w = np.atleast_1d(np.asarray(w, dtype=float))
        V0 = np.exp(KAPPA * self.v_wall)
        V_min = V0 * np.exp(-KAPPA * quad(self.rate, 0, self.w0, epsabs=1e-13, epsrel=1e-13)[0])
        U_min = float(kruskal_product(self.r_min)) / V_min
        U, V = np.empty_like(w), np.empty_like(w)
        for i, x in enumerate(w):
            if x <= self.w0:
                V[i] = V0 * np.exp(-KAPPA * quad(self.rate, 0, x, epsabs=1e-13, epsrel=1e-13)[0])
                U[i] = float(kruskal_product(self.radius(x))) / V[i]
            else:
                U[i] = U_min * np.exp(KAPPA * quad(self.rate, self.w0, x, epsabs=1e-13, epsrel=1e-13)[0])
                V[i] = float(kruskal_product(self.radius(x))) / U[i]
        return U, V

    def static(self, w):
        """(t, r) of the static chart along the moment, t = (ln|V| - ln|U|)/2k in each cell."""
        U, V = self.kruskal(w)
        return (np.log(np.abs(V)) - np.log(np.abs(U))) / (2 * KAPPA), self.radius(w)
