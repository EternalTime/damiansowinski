#!/usr/bin/env python3
"""The charged shell of dust: the shell's motion, and the shells every diagram draws.

A thin spherical shell of dust of rest mass m and charge Q, flat space inside it and
Reissner-Nordstrom's field of mass M and charge Q outside. With mu = G m/c^2, r_s = 2 G M/c^2 and
r_q^2 = G Q^2/(4 pi eps_0 c^4), all lengths, f = 1 - r_s/R + r_q^2/R^2 and a dot for d/d(c tau)
along the shell, Israel's junction condition for dust is

    sqrt(1 + Rdot^2) - sqrt(f + Rdot^2) = mu/R,

which is one first integral, r_s/2 = mu sqrt(1 + Rdot^2) + (r_q^2 - mu^2)/(2R), Kuchar's law of
conservation of energy, M = m sqrt(1 + Rdot^2) + (Q^2 - m^2)/(2R) in G = c = 1. So

    gamma = dT/dtau   = sqrt(1 + Rdot^2)    = r_s/(2 mu) - (r_q^2 - mu^2)/(2 mu R),  the flat time inside,
    beta  = f dt/dtau = +-sqrt(f + Rdot^2)  = r_s/(2 mu) - (r_q^2 + mu^2)/(2 mu R),  the static time outside,
    Rdot^2 = gamma^2 - 1 = beta^2 - f.

gamma, beta and speed_squared below are those for any r_s, r_q and mu, and turn is the radius
(r_q^2 - mu^2)/(r_s - 2 mu) at which gamma = 1 and the shell is at rest.

The shell drawn falls in from infinity and turns round inside the inner horizon: r_s = 1,
r_q = 12/25 and mu = 1/5, so that r_+ = 16/25, r_- = 9/25, E = r_s/(2 mu) = 5/2,
k = (r_q^2 - mu^2)/(2 mu) = 119/250, l = (r_q^2 + mu^2)/(2 mu) = 169/250 and a = E^2 - 1 = 21/4.
Its motion inside is that of a charge in a repulsive Coulomb field in flat space, a hyperbola's
Kepler problem, and in the parameter eta, zero at the turn,

    R   = k (E + cosh eta)/a                   least at eta = 0, R = k/(E - 1) = 119/375
    tau = k (E eta + sinh eta)/a^(3/2)         proper time, zero at the turn
    T   = k (eta + E sinh eta)/a^(3/2)         the flat time inside, zero there too
    Rdot = sqrt(a) sinh eta/(E + cosh eta),    gamma = (1 + E cosh eta)/(E + cosh eta).

Outside, the advanced time v = t + r_* is regular along the way in and the retarded time
u = t - r_* along the way out, with dv/dtau = 1/(beta - Rdot), du/dtau = 1/(beta + Rdot) and the
static time t zero at the turn, so u(eta) = -v(-eta). v is one quadrature over eta <= 0, held in
a table, and v = u + 2 r_* for eta > 0, which runs off to infinity where the shell comes back to
r_-, at eta = ETA_M. The tortoise coordinate is the Tower's of conformal.py,
r_* = r + A_+ ln|r/r_+ - 1| + A_- ln|r/r_- - 1| with A_+- = +-r_+-^2/(r_+ - r_-), zero at r = 0.

_tools/test_charged_shell.py holds these forms to the published metrics and to the junction
conditions, and charged_shell.md beside this file is the derivation.
"""
import math

import numpy as np

# ---------------------------------------------------------------- any shell

def gamma(R, r_s, r_q, mu):
    """dT/dtau on the shell, sqrt(1 + Rdot^2), for T the flat time inside."""
    return r_s / (2 * mu) - (r_q ** 2 - mu ** 2) / (2 * mu * R)


def beta(R, r_s, r_q, mu):
    """f dt/dtau on the shell, for t the static time outside."""
    return r_s / (2 * mu) - (r_q ** 2 + mu ** 2) / (2 * mu * R)


def speed_squared(R, r_s, r_q, mu):
    """(dR/d(c tau))^2 of the shell."""
    return gamma(R, r_s, r_q, mu) ** 2 - 1


def turn(r_s, r_q, mu):
    """The radius at which the shell is at rest, gamma = 1."""
    return (r_q ** 2 - mu ** 2) / (r_s - 2 * mu)


def mass_at_rest(R, r_q, mu):
    """r_s of the shell at rest at the areal radius R: M = m + (Q^2 - m^2)/(2R)."""
    return 2 * mu + (r_q ** 2 - mu ** 2) / R


def isotropic_radius(R, r_s, r_q):
    """The isotropic radius of the sphere of areal radius R outside the outer horizon, the larger
    root of R = eps + r_s/2 + (r_s^2 - 4 r_q^2)/(16 eps)."""
    half = (R - r_s / 2) / 2
    return half + np.sqrt(half ** 2 - (r_s ** 2 - 4 * r_q ** 2) / 16)


def areal_radius(eps, r_s, r_q):
    """The areal radius of the sphere of isotropic radius eps."""
    return eps + r_s / 2 + (r_s ** 2 - 4 * r_q ** 2) / (16 * eps)


def adm_mass(eps, r_q, mu):
    """r_s of Arnowitt, Deser and Misner's shell at rest at the isotropic radius eps, their (7.1):
    M = -eps + sqrt(eps^2 + 2 m_0 eps + Q^2) in G = c = 1, with m_0 the rest mass."""
    return 2 * (-eps + np.sqrt(eps ** 2 + 2 * mu * eps + r_q ** 2))


# ---------------------------------------------------------------- the shell drawn: r_s = 1, r_q = 12/25, mu = 1/5

RQ = 12 / 25
MU = 1 / 5
RP, RM = 16 / 25, 9 / 25                    # the horizons r_+ and r_-
E = 1 / (2 * MU)                            # the mass seen from outside over the rest mass of the dust
K = (RQ ** 2 - MU ** 2) / (2 * MU)
ELL = (RQ ** 2 + MU ** 2) / (2 * MU)
A = E * E - 1
ROOT_A = math.sqrt(A)
R_TURN = K / (E - 1)                        # 119/375, where the shell turns round
ETA_P = math.acosh(A * RP / K - E)          # the shell is on r_+ at eta = -+ETA_P
ETA_M = math.acosh(A * RM / K - E)          # and on r_- at eta = -+ETA_M
A_PLUS, A_MINUS = RP ** 2 / (RP - RM), -RM ** 2 / (RP - RM)


def tortoise(r):
    """r_* = r + A_+ ln|r/r_+ - 1| + A_- ln|r/r_- - 1|, zero at r = 0."""
    r = np.asarray(r, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        return r + A_PLUS * np.log(np.abs(r / RP - 1)) + A_MINUS * np.log(np.abs(r / RM - 1))


def radius(eta):
    return K * (E + np.cosh(np.asarray(eta, dtype=float))) / A


def proper_time(eta):
    eta = np.asarray(eta, dtype=float)
    return K * (E * eta + np.sinh(eta)) / A ** 1.5


def inner_time(eta):
    eta = np.asarray(eta, dtype=float)
    return K * (eta + E * np.sinh(eta)) / A ** 1.5


def inner_advanced(eta):
    """T + R on the shell."""
    return inner_time(eta) + radius(eta)


def inner_retarded(eta):
    """T - R on the shell."""
    return inner_time(eta) - radius(eta)


def speed(eta):
    """dR/d(c tau) on the shell."""
    eta = np.asarray(eta, dtype=float)
    return ROOT_A * np.sinh(eta) / (E + np.cosh(eta))


def eta_of_radius(R, way=-1):
    """The shell at the radius R >= R_TURN, on the way in (way = -1) or out (way = +1)."""
    return way * np.arccosh(np.maximum(A * np.asarray(R, dtype=float) / K - E, 1.0))


def _rate(eta):
    """dv/d(eta) on the way in: (dtau/d eta)/(beta - Rdot), with dtau/d eta = R/sqrt(a)."""
    R = radius(eta)
    return R / ROOT_A / (E - ELL / R - speed(eta))


# v over eta <= 0 as a table: a panel every STEP of eta from -SPAN, each by Gauss and Legendre's
# rule of eight points, and from the last node below eta by the same rule.
STEP, SPAN = 1 / 64, 40.0
_NODES, _WEIGHTS = np.polynomial.legendre.leggauss(8)
_GRID = -SPAN + STEP * np.arange(int(round(SPAN / STEP)) + 1)


def _panel(lo, hi):
    """The integral of _rate over [lo, hi], for arrays of ends."""
    mid, half = (hi + lo) / 2, (hi - lo) / 2
    return half * np.sum(_WEIGHTS * _rate(mid[..., None] + half[..., None] * _NODES), axis=-1)


# Summed from the turn outward, so that the small values near the turn are not differences of the
# large ones far out; t = 0 at the turn, so v(0) = r_*(R_TURN).
_SUMS = float(tortoise(R_TURN)) - np.concatenate([np.cumsum(_panel(_GRID[:-1], _GRID[1:])[::-1])[::-1], [0.0]])


def _advanced_in(eta):
    """v on the shell for -SPAN <= eta <= 0."""
    eta = np.clip(np.asarray(eta, dtype=float), -SPAN, 0.0)
    i = np.minimum(np.floor((eta + SPAN) / STEP).astype(int), _GRID.size - 2)
    return _SUMS[i] + _panel(_GRID[i], eta)


def advanced(eta):
    """The advanced time v = t + r_* on the shell, for eta < ETA_M, where it runs off to infinity."""
    eta = np.asarray(eta, dtype=float)
    with np.errstate(invalid="ignore"):
        return np.where(eta <= 0, _advanced_in(eta), -_advanced_in(-eta) + 2 * tortoise(radius(eta)))


def retarded(eta):
    """The retarded time u = t - r_* on the shell, for eta > -ETA_M: the way in turned over in time."""
    return -advanced(-np.asarray(eta, dtype=float))


def outer_time(eta):
    """The static time t on the shell, zero at the turn; it runs off to infinity at each horizon."""
    eta = np.asarray(eta, dtype=float)
    return np.where(eta <= 0, 1.0, -1.0) * (_advanced_in(-np.abs(eta)) - tortoise(radius(eta)))


def _bisect(f, target, lo, hi, rising=True):
    """The eta in (lo, hi) at which the monotonic f is target, by bisection, for arrays of targets;
    a target f does not reach is given the end it lies beyond."""
    target = np.asarray(target, dtype=float)
    a, b = np.full(target.shape, float(lo)), np.full(target.shape, float(hi))
    for _ in range(90):
        mid = (a + b) / 2
        with np.errstate(all="ignore"):
            below = (f(mid) < target) == rising
        a, b = np.where(below, mid, a), np.where(below, b, mid)
    return (a + b) / 2


_TABLES = {}


def _invert(f, target, lo, hi, rising=True):
    """The same root, found faster: from a table of f over (lo, hi), crowded toward both ends, by
    the secant rule between the two entries that hold the target and then inside that bracket, and
    by bisection wherever the secant rule has not closed on the root."""
    target = np.asarray(target, dtype=float)
    key = (f, lo, hi)
    if key not in _TABLES:
        x = (1 - np.cos(np.linspace(0, np.pi, 4001))) / 2
        eta = lo + (hi - lo) * x
        with np.errstate(all="ignore"):
            values = f(eta)
        good = np.isfinite(values)
        _TABLES[key] = (eta[good], values[good] if rising else -values[good])
    eta, values = _TABLES[key]
    want = target if rising else -target
    i = np.clip(np.searchsorted(values, want), 1, eta.size - 1)
    a, b = eta[i - 1], eta[i]
    fa, fb = values[i - 1] - want, values[i] - want
    inside = (fa <= 0) & (fb >= 0)
    x = np.where(inside, a, (a + b) / 2)
    for _ in range(60):
        with np.errstate(all="ignore"):
            x = np.where(inside & (fb != fa), a - fa * (b - a) / (fb - fa), (a + b) / 2)
            # The Illinois rule: a secant step, and the end that stays is halved so that both ends close in.
            fx = (f(x) if rising else -f(x)) - want
        left = fx < 0
        fa, fb = np.where(left, fx, np.where(inside, fa / 2, fa)), np.where(left, np.where(inside, fb / 2, fb), fx)
        a, b = np.where(left, x, a), np.where(left, b, x)
        if np.all(~inside | (np.abs(b - a) < 4e-16 * np.maximum(1.0, np.abs(x))) | (fx == 0)):
            break
    done = inside & ((np.abs(b - a) < 1e-14 * np.maximum(1.0, np.abs(x))) | (fx == 0))
    if not np.all(done):
        x = np.where(done, x, _bisect(f, target, lo, hi, rising))
    return x


TINY = 1e-13
FAR = 38.0


def eta_of_inner_time(T):
    return _invert(inner_time, T, -FAR, FAR)


def eta_of_inner_advanced(w):
    """The shell where the ingoing ray T + r = w of the interior meets it."""
    return _invert(inner_advanced, w, -FAR, FAR)


def eta_of_inner_retarded(w):
    """The shell where the outgoing ray T - r = w of the interior meets it."""
    return _invert(inner_retarded, w, -FAR, FAR)


def eta_of_advanced(v):
    """The shell at the advanced time v, before it comes back to r_-."""
    return _invert(advanced, v, -FAR, ETA_M - TINY)


def eta_of_retarded(u):
    """The shell at the retarded time u, after it has crossed r_- on the way in."""
    return -eta_of_advanced(-np.asarray(u, dtype=float))


def slice_time(eta):
    """v - R on the shell, which labels the slices of constant v - r of the ingoing chart."""
    return advanced(eta) - radius(eta)


def eta_of_slice(w):
    """The shell on the slice v - r = w of the ingoing chart."""
    return _invert(slice_time, w, -FAR, ETA_M - TINY)


def eta_of_outer_time(t, region="I"):
    """The shell at the static time t: outside r_+ on the way in (I), inside r_- (III), or outside
    r_+ on the way out (I')."""
    if region == "I":
        return _invert(outer_time, t, -FAR, -ETA_P - TINY)
    if region == "III":
        return _invert(outer_time, t, -ETA_M + TINY, ETA_M - TINY)
    return _invert(outer_time, t, ETA_P + TINY, FAR)


def radius_at_inner_time(T):
    """R(T), the shell in the flat chart inside it."""
    return radius(eta_of_inner_time(T))


def radius_at_advanced(v):
    """R(v), the shell in the ingoing chart."""
    return radius(eta_of_advanced(v))


def radius_at_retarded(u):
    """R(u), the shell in the outgoing chart."""
    return radius(eta_of_retarded(u))


def radius_at_outer_time(t):
    """R(t) outside r_+, the shell on its way in."""
    return radius(eta_of_outer_time(t, "I"))


def radius_at_static_time_inside(t):
    """R(t) inside r_-, where the shell turns round at t = 0."""
    return radius(eta_of_outer_time(t, "III"))


def radius_at_outer_time_leaving(t):
    """R(t) outside r_+ of the next exterior, the shell on its way out."""
    return radius(eta_of_outer_time(t, "I'"))


def radius_on_slice(w):
    """The shell's radius on the slice v - r = w."""
    return radius(eta_of_slice(w))


# ---------------------------------------------------------------- the conformal diagram

# Inside the shell p, q = Phi(T -+ r) with Phi(w) = arctan(w/CONFORMAL_LENGTH), Minkowski's own
# compactification, which puts the centre on the straight line p = q. Outside, each ray keeps the
# p or the q it has where it crosses the shell, which makes both continuous there: every outgoing
# ray outside the shell left it once, and every ingoing ray meets it once. The length is T + R of
# the shell where it comes back to r_-, so that the two branches of the inner horizon it crosses
# are the lines p = -pi/4 and q = pi/4.
CONFORMAL_LENGTH = float(inner_advanced(ETA_M))


def phi(w):
    return np.arctan(np.asarray(w, dtype=float) / CONFORMAL_LENGTH)


def _other(eta):
    """v - 2r_* on the shell, constant along the other family of rays of each region of the
    ingoing chart: the retarded time u outside r_+ and inside r_-."""
    return advanced(eta) - 2 * tortoise(radius(eta))


def exit_eta(v, r):
    """Where the outgoing ray through the event (v, r) of the ingoing chart left the shell. The
    ray keeps v - 2r_*, which the shell has once on each stretch of its way: before r_+ for an
    event outside r_+, between the horizons for one between them, and after r_- for one inside r_-."""
    v, r = np.broadcast_arrays(np.asarray(v, dtype=float), np.asarray(r, dtype=float))
    other = v - 2 * tortoise(r)
    out = np.empty(v.shape)
    far, near = r > RP, r < RM
    mid = ~far & ~near
    if far.any():
        out[far] = _invert(_other, other[far], -FAR, -ETA_P - TINY)
    if mid.any():
        on = (r[mid] == RP) | (r[mid] == RM)
        found = _invert(_other, np.where(on, 0.0, other[mid]), -ETA_P + TINY, -ETA_M - TINY, rising=False)
        out[mid] = np.where(r[mid] == RP, -ETA_P, np.where(r[mid] == RM, -ETA_M, found))
    if near.any():
        out[near] = eta_of_retarded(other[near])
    return out


def conformal_inside(T, r):
    """(p, q) of an event of the flat interior."""
    T, r = np.asarray(T, dtype=float), np.asarray(r, dtype=float)
    return phi(T - r), phi(T + r)


def conformal_ingoing(v, r):
    """(p, q) of an event outside the shell in the ingoing chart's v and r."""
    v, r = np.broadcast_arrays(np.asarray(v, dtype=float), np.asarray(r, dtype=float))
    return phi(inner_retarded(exit_eta(v, r))), phi(inner_advanced(eta_of_advanced(v)))


def conformal_outgoing(u, r):
    """(p, q) of an event outside the shell in the outgoing chart's u and r: the ingoing chart
    turned over in time about the turn."""
    p, q = conformal_ingoing(-np.asarray(u, dtype=float), r)
    return -q, -p


def _between(other):
    """p of the outgoing ray that left the shell between the horizons with v - 2r_* = other."""
    return phi(inner_retarded(_invert(_other, other, -ETA_P + TINY, -ETA_M - TINY, rising=False)))


def conformal_beyond(t, r):
    """(p, q) of an event inside r_- on the far side of the Cauchy horizon, in the static time t
    of that region, zero on the moment of time symmetry and growing to the future: its outgoing ray
    left the shell between the horizons, where it has v - 2r_* = -t - r_*, and its ingoing ray is the
    mirror image of such a ray."""
    t, r = np.asarray(t, dtype=float), np.asarray(r, dtype=float)
    rstar = tortoise(r)
    return _between(-t - rstar), -_between(t - rstar)


# ---------------------------------------------------------------- the point charge: the shell at rest

# Arnowitt, Deser and Misner's shell: charge r_q = 1, the unit of every length there, at rest at
# the isotropic radius eps. With no rest mass, r_s = 2(sqrt(eps^2 + 1) - eps), the shell's areal
# radius is R = eps + r_s/4 = r_q^2/r_s, where g_rr = 1 on both sides, and r_s rises to 2 r_q, the
# mass of the extreme Reissner-Nordstrom field, as eps falls to zero.

def point_charge(eps, mu=0.0):
    """(r_s, R) of the shell of unit charge radius and rest mass mu at rest at the isotropic radius eps."""
    r_s = adm_mass(eps, 1.0, mu)
    return r_s, areal_radius(eps, r_s, 1.0)
