#!/usr/bin/env python3
"""Israel's collapsing shell of dust: the shell's motion, and the one shell every diagram draws.

A thin spherical shell of dust of rest mass m, flat space inside it and Schwarzschild's vacuum of
mass M outside. With mu = G m/c^2 and r_s = 2 G M/c^2, both lengths, and a dot for d/d(c tau)
along the shell, Israel's junction condition for dust is

    sqrt(1 + Rdot^2) - sqrt(1 - r_s/R + Rdot^2) = mu/R,

which is one first integral, r_s/2 = mu sqrt(1 + Rdot^2) - mu^2/(2R). So with E = r_s/(2 mu), the
energy of the shell per unit rest mass,

    gamma = dT/dtau         = sqrt(1 + Rdot^2)         = E + mu/(2R),      the flat time inside,
    beta  = (1 - r_s/R) dt/dtau = +-sqrt(1 - r_s/R + Rdot^2) = E - mu/(2R),  Schwarzschild's time outside,
    Rdot^2 = gamma^2 - 1 = beta^2 - (1 - r_s/R).

gamma, beta and speed below are those three for any r_s and mu. beta changes sign at R = mu/(2E),
which lies inside r_s for every shell.

The shell drawn is the one that falls from rest at infinity, E = 1, mu = r_s/2, in units of r_s.
In s = sqrt(1 + 8R) every time along it is a polynomial, or a polynomial and one logarithm:

    R     = (s^2 - 1)/8
    tau   = -(s - 1)^2 (s + 2)/24                      proper time, zero where the shell reaches R = 0
    T     = -(s^3 + 3s - 4)/24                          the flat time inside, zero there too
    T + R = -(s - 1)^3/24,   T - R = -((s + 1)^3 - 8)/24
    v     = -(s^3/3 - s^2 + 5s - 15 - 16 ln((s + 3)/6))/8   the advanced time outside, zero on r_s

since dR/dtau = -2s/(s^2 - 1), gamma = (s^2 + 1)/(s^2 - 1), beta = (s^2 - 3)/(s^2 - 1),
1 - 1/R = (s^2 - 9)/(s^2 - 1) and dv/dtau = (beta + Rdot)/(1 - 1/R) = (s + 1)/(s + 3). The shell
crosses r_s at s = 3, where tau = -5/6, T = -4/3 and v = 0, and reaches R = 0 at s = 1, where
v = V_END. The event horizon is the outgoing ray T - r = -7/3 inside the shell and r = r_s outside
it, so it leaves the centre at T = -7/3, before the shell arrives.

_tools/test_israel_shell.py holds these forms to the published metrics and to the junction
conditions, and israel_shell.md beside this file is the derivation.
"""
import math

import numpy as np

HORIZON_S = 3.0                                 # s where the shell crosses r_s
HORIZON_U = -7.0 / 3                            # T - r of the event horizon inside the shell
V_END = (32.0 / 3 - 16 * math.log(1.5)) / 8     # the advanced time at which the shell reaches R = 0


# ---------------------------------------------------------------- any shell

def energy(r_s, mu):
    """E = r_s/(2 mu), the mass seen from outside over the rest mass of the dust."""
    return r_s / (2 * mu)


def gamma(R, r_s, mu):
    """dT/dtau on the shell, sqrt(1 + Rdot^2), for T the flat time inside."""
    return energy(r_s, mu) + mu / (2 * R)


def beta(R, r_s, mu):
    """(1 - r_s/R) dt/dtau on the shell, for t Schwarzschild's time outside."""
    return energy(r_s, mu) - mu / (2 * R)


def speed(R, r_s, mu):
    """dR/d(c tau) of a shell falling in."""
    return -np.sqrt(gamma(R, r_s, mu) ** 2 - 1)


# ---------------------------------------------------------------- the shell drawn, r_s = 1 and mu = 1/2

def s_of_radius(R):
    return np.sqrt(1 + 8 * np.asarray(R, dtype=float))


def radius(s):
    s = np.asarray(s, dtype=float)
    return (s * s - 1) / 8


def proper_time(s):
    s = np.asarray(s, dtype=float)
    return -(s - 1) ** 2 * (s + 2) / 24


def inner_time(s):
    s = np.asarray(s, dtype=float)
    return -(s ** 3 + 3 * s - 4) / 24


def inner_advanced(s):
    """T + R on the shell."""
    s = np.asarray(s, dtype=float)
    return -(s - 1) ** 3 / 24


def inner_retarded(s):
    """T - R on the shell."""
    s = np.asarray(s, dtype=float)
    return -((s + 1) ** 3 - 8) / 24


def advanced(s):
    """Schwarzschild's advanced time v = ct + r + ln|r - 1| on the shell, zero where it crosses r_s."""
    s = np.asarray(s, dtype=float)
    return -(s ** 3 / 3 - s * s + 5 * s - 15 - 16 * np.log((s + 3) / 6)) / 8


def tortoise(r):
    r = np.asarray(r, dtype=float)
    with np.errstate(divide="ignore"):
        return r + np.log(np.abs(r - 1))


def outer_time(s):
    """Schwarzschild's time ct on the shell, for s > 3; it grows without limit as s falls to 3."""
    return advanced(s) - tortoise(radius(s))


def kruskal(s):
    """Kruskal's U = -e^(-u/2) and V = e^(v/2) on the shell, with u = v - 2r_*, so UV = (1 - R)e^R;
    both are finite and smooth through r_s."""
    R, v = radius(s), advanced(s)
    return (1 - R) * np.exp(R - v / 2), np.exp(v / 2)


def s_of_inner_time(T):
    """The shell at the flat time T <= 0: the root of s^3 + 3s = 4 - 24T, by Cardano's formula."""
    T = np.minimum(np.asarray(T, dtype=float), 0.0)
    return 2 * np.sinh(np.arcsinh(2 - 12 * T) / 3)


def s_of_inner_advanced(w):
    """The shell where the ingoing ray T + r = w <= 0 of the interior crosses it."""
    w = np.minimum(np.asarray(w, dtype=float), 0.0)
    return 1 + np.cbrt(-24 * w)


def s_of_inner_retarded(w):
    """The shell where the outgoing ray T - r = w <= 0 of the interior crosses it."""
    w = np.minimum(np.asarray(w, dtype=float), 0.0)
    return np.cbrt(8 - 24 * w) - 1


def _invert(f, df, target, lo, hi, guess):
    """The root s of f(s) = target for a function falling from lo to hi, by Newton's rule kept
    inside the bracket, for arrays of targets."""
    target = np.asarray(target, dtype=float)
    a, b = np.full(target.shape, float(lo)), np.full(target.shape, float(hi))
    s = np.clip(np.asarray(guess, dtype=float) * np.ones(target.shape), a, b)
    for _ in range(200):
        miss = f(s) - target
        a, b = np.where(miss > 0, s, a), np.where(miss > 0, b, s)      # f falls, so a miss above is left of the root
        with np.errstate(all="ignore"):
            step = s - miss / df(s)
        s = np.where((step > a) & (step < b) & np.isfinite(step), step, (a + b) / 2)
        if np.all(np.abs(b - a) < 1e-15 * np.maximum(1.0, np.abs(s))):
            break
    return s


def _d_advanced(s):
    return -(s + 1) ** 2 * (s - 1) / (8 * (s + 3))


def s_of_advanced(v):
    """The shell at Schwarzschild's advanced time v <= V_END."""
    v = np.minimum(np.asarray(v, dtype=float), V_END)
    top = np.maximum(4.0, 2 * np.cbrt(24 * np.abs(v) + 64))
    return _invert(advanced, _d_advanced, v, 1.0, float(np.max(top)), np.cbrt(24 * np.abs(v) + 27))


def s_of_slice(w):
    """The shell on the slice v - r = w of the ingoing chart, w <= V_END."""
    w = np.minimum(np.asarray(w, dtype=float), V_END)
    top = np.maximum(4.0, 2 * np.cbrt(24 * np.abs(w) + 64))
    return _invert(lambda s: advanced(s) - radius(s), lambda s: _d_advanced(s) - s / 4, w, 1.0, float(np.max(top)),
                   np.cbrt(24 * np.abs(w) + 27))


def s_of_outer_time(t):
    """The shell at Schwarzschild's time ct, outside r_s: the root above s = 3 of outer_time(s) = ct."""
    t = np.asarray(t, dtype=float)

    def f(x):
        # In x = ln(s - 3) the time falls from +infinity at x = -infinity, and the root is well conditioned.
        return outer_time(3 + np.exp(x))

    def df(x):
        s = 3 + np.exp(x)
        R = radius(s)
        return (_d_advanced(s) - (s / 4) * R / (R - 1)) * np.exp(x)
    hi = np.log(np.maximum(4.0, 2 * np.cbrt(24 * np.abs(t) + 64)))
    x = _invert(f, df, t, -700.0, float(np.max(hi)), 0.0)
    return 3 + np.exp(x)


def radius_at_inner_time(T):
    """R(T), the shell in the flat chart inside it; zero from T = 0 on."""
    T = np.asarray(T, dtype=float)
    return np.where(T < 0, radius(s_of_inner_time(T)), 0.0)


def radius_at_advanced(v):
    """R(v), the shell in the ingoing chart; zero from v = V_END on."""
    v = np.asarray(v, dtype=float)
    return np.where(v < V_END, radius(s_of_advanced(v)), 0.0)


def radius_at_outer_time(t):
    """R(t), the shell in Schwarzschild's chart, which it leaves only as ct grows without limit."""
    return radius(s_of_outer_time(t))


def radius_on_slice(w):
    """The shell's radius on the slice v - r = w; zero from w = V_END on."""
    w = np.asarray(w, dtype=float)
    return np.where(w < V_END, radius(s_of_slice(w)), 0.0)


# ---------------------------------------------------------------- the conformal diagram

# Inside the shell p, q = Phi(T -+ r) with Phi(w) = arctan(w/CONFORMAL_LENGTH), Minkowski's own
# compactification, which puts the centre on the straight line p = q. The length is 7/3, so that
# the event horizon T - r = -7/3 is the line p = -pi/4. Outside, each ray keeps the p or the q it
# has where it crosses the shell, which makes both continuous there: an outgoing ray left the shell
# at exit_s, an ingoing ray with v < V_END enters it at s_of_advanced(v), and an ingoing ray with
# v >= V_END, which ends on the singularity at Kruskal's U = e^(-v/2), takes q = -p of the outgoing
# ray that ends there with it, so the singularity is the straight line p + q = 0.
CONFORMAL_LENGTH = 7.0 / 3


def phi(w):
    return np.arctan(np.asarray(w, dtype=float) / CONFORMAL_LENGTH)


def _retarded(s):
    """Schwarzschild's retarded time u = v - 2r_* on the shell outside r_s, s > 3."""
    return advanced(s) - 2 * tortoise(radius(s))


def exit_s(v, r):
    """Where the outgoing ray through the event (v, r) outside the shell left the shell. Outside r_s
    the ray keeps u = v - 2r_*, and the shell had that u at one s > 3; inside r_s it keeps
    Kruskal's U = (1 - r) e^(r - v/2) > 0, which the shell had at one s < 3; on r_s it is s = 3."""
    v, r = np.broadcast_arrays(np.asarray(v, dtype=float), np.asarray(r, dtype=float))
    out = np.full(v.shape, HORIZON_S)
    far = r > 1
    if far.any():
        u = v[far] - 2 * tortoise(r[far])

        def f(x):
            return _retarded(3 + np.exp(x))

        def df(x):
            s = 3 + np.exp(x)
            R = radius(s)
            return (_d_advanced(s) - (s / 2) * R / (R - 1)) * np.exp(x)
        hi = np.log(np.maximum(4.0, 2 * np.cbrt(24 * np.abs(u) + 64)))
        out[far] = 3 + np.exp(_invert(f, df, u, -700.0, float(np.max(hi)), 0.0))
    near = r < 1
    if near.any():
        U = (1 - r[near]) * np.exp(r[near] - v[near] / 2)

        def g(s):
            return kruskal(s)[0]

        def dg(s):
            R = radius(s)
            return -np.exp(R - advanced(s) / 2) * (R * s / 4 + (1 - R) * _d_advanced(s) / 2)
        out[near] = _invert(g, dg, np.minimum(U, kruskal(1.0)[0]), 1.0, HORIZON_S, 2.0)
    return out


def conformal_inside(T, r):
    """(p, q) of an event of the flat interior."""
    T, r = np.asarray(T, dtype=float), np.asarray(r, dtype=float)
    return phi(T - r), phi(T + r)


def conformal_outside(v, r):
    """(p, q) of an event outside the shell, in the ingoing chart's v and r."""
    v, r = np.broadcast_arrays(np.asarray(v, dtype=float), np.asarray(r, dtype=float))
    p = phi(inner_retarded(exit_s(v, r)))
    early = v < V_END
    q = np.empty(v.shape)
    if early.any():
        q[early] = phi(inner_advanced(s_of_advanced(v[early])))
    if not early.all():
        # The outgoing ray that reaches r = 0 at the advanced time v has U = e^(-v/2): the event (v, 0).
        q[~early] = -phi(inner_retarded(exit_s(v[~early], np.zeros(int((~early).sum())))))
    return p, q
