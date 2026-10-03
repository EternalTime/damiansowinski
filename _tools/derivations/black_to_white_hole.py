#!/usr/bin/env python3
"""Haggard and Rovelli's black hole fireworks: the one bounce every diagram draws.

Hal M. Haggard and Carlo Rovelli, Phys. Rev. D 92, 104020 (2015), arXiv:1407.0989. A thin shell
of light falls in from past null infinity. Inside it space is flat, region I; outside it, below the
quantum region, it is a piece of Kruskal's spacetime, region II; the quantum region III is bounded
by a spacelike line from the point E on the shell to the point Delta on the surface of time
symmetry; and the whole is glued to its time reverse along that surface outside Delta, so that the
shell comes back out of a white hole into the same flat space it fell from. In units of r_s = 2m,
with c = 1:

    region II    ds^2 = -(4/r) e^(-r) dU dV + r^2 dOmega^2,   (1 - r) e^r = UV,
                 their (3) and (4) with 32 m^3 = 4 r_s^3, the shell falling in along V = V0, and
                 the surface of time symmetry U + V = 0 outside Delta;
    region I     ds^2 = -du dv + ((v - u)/2)^2 dOmega^2, the shell along v = 0;
    the junction r = (v - u)/2 on v = 0, so an outgoing ray u of region I crosses the shell into
                 region II at U = (1 + u/2) e^(-u/2)/V0. Haggard and Rovelli's (28) to (30) print
                 the exponent with the other sign, e^(+u/4m); the radius r = -u/2 on the shell
                 makes (1 - r) e^r = (1 + u/2) e^(-u/2), the form used here, which is
                 black_to_white_hole.md's Step 2;
    E            the point of the shell with u = -2 eps, their choice, at r = eps;
    Delta        the point of U + V = 0 with r = 2m + delta, their delta = m/3, so r = 7/6;
    E to Delta   the spacelike geodesic between the two, their choice.

The drawings take V0 = 3/5 and eps = 1/2, the values for which a spacelike geodesic from Delta
reaches E without leaving the region between E's radius and Delta's. For V0 = 1/2 it reaches the
shell no further in than r = 0.497, and the smaller V0 is, the further out that radius lies; the
derivation is black_to_white_hole.md's Step 4. Haggard and Rovelli's own V0 is
exp(-k m/2 l_P), far smaller, and their eps is Planckian: the drawing is a picture of the
construction at sizes where its parts can be seen, and every view says so in its input.

The time reverse is (U, V) -> (-V, -U), which maps region II onto the flap after the bounce, where
the shell goes out along U = -V0. Every chart's time is set to zero on the surface of time
symmetry: Schwarzschild's t, the Painleve-Gullstrand time and Lemaitre's tau.

_tools/test_black_to_white_hole.py holds the geodesic, the junction and the maps to the published
metrics.
"""
import math
from functools import lru_cache

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq
from scipy.special import lambertw

V0 = 0.6                                         # the shell falls in along V = V0
EPS = 0.5                                        # E is the point u = -2 eps of the shell, at r = eps
R_DELTA = 7.0 / 6                                # Delta, at r = 2m + m/3
U_E = (1 - EPS) * math.exp(EPS) / V0             # Kruskal's U at E
V_DELTA = math.sqrt((R_DELTA - 1) * math.exp(R_DELTA))   # Delta is (-V_DELTA, V_DELTA)
U_HORIZON = -2.0                                 # the outgoing ray of region I that becomes r = r_s


def radius(U, V):
    """Kruskal's areal radius, (1 - r) e^r = UV, on the principal branch of Lambert's function."""
    return 1 + np.real(lambertw(-np.asarray(U, dtype=float) * np.asarray(V, dtype=float) / np.e))


def shell_U(u):
    """Kruskal's U of the outgoing ray u of the flat interior where it crosses the shell."""
    u = np.asarray(u, dtype=float)
    return (1 + u / 2) * np.exp(-u / 2) / V0


def interior_u(U):
    """The inverse of shell_U on u <= 0, where it rises: the outgoing ray of region I that is U."""
    U = np.asarray(U, dtype=float)
    out = np.full(U.shape, np.nan)
    for i, x in np.ndenumerate(U):
        if x <= shell_U(0.0):
            lo = -2.0
            while shell_U(lo) > x:
                lo *= 2
            out[i] = brentq(lambda u: shell_U(u) - x, lo, 0.0, xtol=1e-14)
    return out


def shell_crossing():
    """The radius at which the shell crosses U + V = 0: (1 - r) e^r = -V0^2."""
    return brentq(lambda r: (1 - r) * math.exp(r) + V0 ** 2, 1.0 + 1e-12, 3.0, xtol=1e-14)


def _rhs(s, y):
    """The radial geodesic equations of -(4/r) e^(-r) dU dV: U'' = -(d_U ln F) U'^2 and the same in V,
    with d ln F/dr = -1/r - 1, d_U r = -V e^(-r)/r and d_V r = -U e^(-r)/r."""
    U, V, dU, dV = y
    r = float(radius(U, V))
    k = (-1 / r - 1) * math.exp(-r) / r
    return [dU, dV, k * V * dU ** 2, k * U * dV ** 2]


def _shoot(angle):
    """The geodesic leaving Delta toward smaller V with dU = cos(angle), dV = -sin(angle), to the shell."""
    hit = lambda s, y: y[1] - V0  # noqa: E731
    hit.terminal = True
    return solve_ivp(_rhs, (0, 50), [-V_DELTA, V_DELTA, math.cos(angle), -math.sin(angle)], events=hit,
                     rtol=1e-12, atol=1e-13, dense_output=True)


@lru_cache(maxsize=1)
def geodesic():
    """The spacelike geodesic from Delta to E, as arrays (U, V, r) from Delta to E. Two geodesics
    leave Delta for E; the one taken keeps r between E's and Delta's, falling all the way, where the
    other runs in toward r = 0.012 and back out."""
    def miss(angle):
        sol = _shoot(angle)
        return sol.y_events[0][0][0] - U_E if sol.t_events[0].size else -U_E
    angles = np.geomspace(1e-5, 1.0, 300)
    misses = np.array([miss(a) for a in angles])
    roots = [brentq(miss, angles[i], angles[i + 1], xtol=1e-15)
             for i in range(len(angles) - 1) if np.sign(misses[i]) != np.sign(misses[i + 1])]
    best = None
    for angle in roots:
        sol = _shoot(angle)
        s = np.linspace(0, sol.t_events[0][0], 4001)
        U, V = sol.sol(s)[:2]
        r = radius(U, V)
        if np.all(np.diff(r) < 0):
            best = (U, V, r)
    if best is None:
        raise AssertionError("black_to_white_hole: no geodesic from Delta to E keeps r between theirs")
    U, V, r = best
    U[-1], V[-1], r[-1] = U_E, V0, EPS
    return U, V, r


@lru_cache(maxsize=1)
def _boundary_of_U():
    U, V, _ = geodesic()
    return PchipInterpolator(U, V)


@lru_cache(maxsize=1)
def _boundary_of_r():
    U, V, r = geodesic()
    return PchipInterpolator(r[::-1], U[::-1]), PchipInterpolator(r[::-1], V[::-1])


def sigma_V(U):
    """V on the future edge of region II above the outgoing ray U: the surface of time symmetry
    V = -U outside Delta, the geodesic to E, and the shell V = V0 beyond E, so that region II is
    V0 <= V < sigma_V(U)."""
    U = np.asarray(U, dtype=float)
    out = np.where(U <= -V_DELTA, -U, V0)
    mid = (U > -V_DELTA) & (U < U_E)
    if np.any(mid):
        out = np.where(mid, _boundary_of_U()(np.clip(U, -V_DELTA, U_E)), out)
    return out


def boundary_UV(r):
    """The point of the future edge of region II at the radius r, for r >= eps: on the geodesic
    for r <= 7/6 and on U + V = 0 beyond."""
    r = np.asarray(r, dtype=float)
    fU, fV = _boundary_of_r()
    rr = np.clip(r, EPS, R_DELTA)
    far = np.sqrt(np.maximum(r - 1, 0) * np.exp(r))
    return np.where(r <= R_DELTA, fU(rr), -far), np.where(r <= R_DELTA, fV(rr), far)


def region_two(U, V):
    """Positive in region II: outside the shell and below the quantum region and the surface of time symmetry."""
    return np.minimum(np.asarray(V, dtype=float) - V0, sigma_V(U) - np.asarray(V, dtype=float))


# -- the charts outside the shell, each put in Kruskal's U and V -------------------------------

def schwarzschild_UV(t, r):
    """Schwarzschild's t and r outside r_s: U = -sqrt(r - 1) e^((r - t)/2), V = sqrt(r - 1) e^((r + t)/2)."""
    t, r = np.asarray(t, dtype=float), np.asarray(r, dtype=float)
    a = np.sqrt(np.maximum(r - 1, 0)) * np.exp(r / 2)
    return -a * np.exp(-t / 2), a * np.exp(t / 2)


def painleve_UV(t, r):
    """The ingoing Painleve-Gullstrand time t = T + 2 sqrt(r) + ln|(sqrt(r) - 1)/(sqrt(r) + 1)|, with T
    Schwarzschild's, on either side of r_s: V = (sqrt(r) + 1) e^((t - 2 sqrt(r) + r)/2), U = (1 - r) e^r/V."""
    t, r = np.asarray(t, dtype=float), np.asarray(r, dtype=float)
    V = (np.sqrt(r) + 1) * np.exp((t - 2 * np.sqrt(r) + r) / 2)
    return (1 - r) * np.exp(r) / V, V


def painleve_time(U, V):
    """The ingoing Painleve-Gullstrand time of a point of region II."""
    r = radius(U, V)
    return 2 * np.log(np.asarray(V, dtype=float) / (np.sqrt(r) + 1)) + 2 * np.sqrt(r) - r


def lemaitre_UV(tau, rho):
    """Lemaitre's tau is the Painleve-Gullstrand time and r = (3 (rho - tau)/2)^(2/3)."""
    tau, rho = np.asarray(tau, dtype=float), np.asarray(rho, dtype=float)
    r = np.cbrt(1.5 * np.maximum(rho - tau, 0)) ** 2
    return painleve_UV(tau, r)


def shell_painleve_time(r):
    """The Painleve-Gullstrand time at which the shell, V = V0, passes the radius r."""
    r = np.asarray(r, dtype=float)
    return 2 * np.log(V0 / (np.sqrt(r) + 1)) + 2 * np.sqrt(r) - r


def shell_schwarzschild_time(r):
    """Schwarzschild's time at which the shell passes the radius r > r_s."""
    r = np.asarray(r, dtype=float)
    return 2 * math.log(V0) - r - np.log(r - 1)


def _inverse(f, x, lo, hi):
    """The r in (lo, hi) at which the falling f(r) is x, NaN where there is none."""
    out = np.full(np.shape(x), np.nan)
    for i, value in np.ndenumerate(np.asarray(x, dtype=float)):
        if f(hi) <= value <= f(lo):
            out[i] = brentq(lambda r: f(r) - value, lo, hi, xtol=1e-13)
    return out if np.ndim(x) else float(out)


def shell_radius_painleve(t):
    """The radius of the shell at the Painleve-Gullstrand time t."""
    return _inverse(lambda r: float(shell_painleve_time(r)), t, 1e-12, 1e6)


def shell_radius_schwarzschild(t):
    """The radius of the shell at Schwarzschild's time t."""
    return _inverse(lambda r: float(shell_schwarzschild_time(r)), t, 1 + 1e-15, 1e6)


def edge_painleve_time(r):
    """The Painleve-Gullstrand time of the future edge of region II at the radius r >= eps."""
    U, V = boundary_UV(r)
    return painleve_time(U, V)


def edge_schwarzschild_time(r):
    """Schwarzschild's time of the future edge of region II at the radius r > r_s: zero beyond Delta,
    and ln(-V/U) on the geodesic, which runs to t = infinity on r_s."""
    U, V = boundary_UV(r)
    return np.log(-V / U)


# -- the conformal diagram ---------------------------------------------------------------------

def _H(V):
    """-arctan of the interior ray that reaches the geodesic at V, for V0 <= V <= V_DELTA."""
    U, Vg, _ = geodesic()
    Ug = np.interp(V, Vg[::-1], U[::-1])
    return -np.arctan(interior_u(Ug))


def q_of_V(V):
    """The drawing's q of an ingoing ray V of region II. Beyond Delta's ray it mirrors p across the
    surface of time symmetry, q(V) = -p(-V), so that U + V = 0 is the line T = 0. Below it,
    q = H(V) w(V) with w = (V - V0)/(V_DELTA - V0), zero on the shell, where region I's q = arctan(v)
    is zero too, and below H, the value that would put the geodesic's point on V on T = 0, so that
    the geodesic lies below T = 0 and its mirror image above, with the quantum region between them.
    q has a corner on V = V_DELTA, the ray through Delta, which carries the corner of the gluing at Delta."""
    V = np.asarray(V, dtype=float)
    out = np.full(V.shape, np.nan)
    far = V >= V_DELTA
    out[far] = -np.arctan(interior_u(-V[far]))
    near = (V >= V0) & ~far
    if np.any(near):
        out[near] = _H(V[near]) * (V[near] - V0) / (V_DELTA - V0)
    return out


def p_of_U(U):
    """The drawing's p of an outgoing ray U of region II: arctan of the interior ray it continues."""
    return np.arctan(interior_u(U))


def conformal_two(U, V):
    """(p, q) of a point of region II."""
    return p_of_U(U), q_of_V(V)


def conformal_two_reversed(U, V):
    """(p, q) of a point (U, V) of the flap after the bounce: the time reverse of region II's point (-V, -U)."""
    p, q = conformal_two(-np.asarray(V, dtype=float), -np.asarray(U, dtype=float))
    return -q, -p


def conformal_one(u, v):
    """(p, q) of a point of the flat interior: Minkowski's own arctangents."""
    return np.arctan(np.asarray(u, dtype=float)), np.arctan(np.asarray(v, dtype=float))
