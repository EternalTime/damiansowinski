"""The minimal surfaces of Brill and Lindquist's initial data for two black holes at rest.

The slice is psi^4 times flat space, psi = 1 + alpha_1/r_1 + alpha_2/r_2, with the holes on the
axis at z = +a and z = -a of cylindrical coordinates (rho, phi, z), every length in one unit. A
surface of revolution about the axis is a curve (rho(s), z(s)) of the half plane, s its flat arc
length and theta the angle of its tangent, (cos theta, sin theta). Its area is the integral of
2 pi psi^4 rho ds, so it is minimal exactly when the curve is a geodesic of the half plane with
the metric (psi^4 rho)^2 (drho^2 + dz^2), which is

    dtheta/ds = n . grad ln(psi^4 rho),    n = (-sin theta, cos theta),

and on the axis, where the surface closes smoothly, dtheta/ds = 2 d_z psi/psi. At a moment of time
symmetry a minimal surface is a marginally trapped one, and the outermost is the apparent horizon.

A curve is shot from a point z_0 of the axis at right angles to it. It is the throat of the hole
at z = +a when it comes back to the axis at right angles between the holes, and for two equal
holes it is a surface around both when it crosses the midplane z = 0 at right angles. Brill and
Lindquist found the second kind for separations below 1.56 in units of 2 alpha (Phys. Rev. 131,
471, 1963), Bishop 1.53 (Gen. Rel. Grav. 14, 717, 1982), and Alcubierre and others 1.532
(Class. Quantum Grav. 17, 2159, 2000); critical() finds 1.5319.

It needs numpy and scipy and nothing else, so the tests run it as it stands.
"""
import math

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar

RTOL, ATOL = 1e-11, 1e-13
AXIS = 1e-7         # how near the axis a returning curve is stopped
FAR = 60.0          # a curve that strays this far has missed


def psi(a, rho, z, alpha=(1.0, 1.0)):
    """psi and its two derivatives at (rho, z), for holes of parameters alpha at z = +a and -a."""
    r1, r2 = math.hypot(rho, z - a), math.hypot(rho, z + a)
    return (1 + alpha[0] / r1 + alpha[1] / r2,
            -alpha[0] * rho / r1 ** 3 - alpha[1] * rho / r2 ** 3,
            -alpha[0] * (z - a) / r1 ** 3 - alpha[1] * (z + a) / r2 ** 3)


def turning(a, rho, z, theta, alpha=(1.0, 1.0)):
    """dtheta/ds of the minimal surface through (rho, z) with tangent angle theta."""
    p, pr, pz = psi(a, rho, z, alpha)
    if rho < 1e-9:
        return 2 * pz / p * (1 if math.cos(theta) > 0 else -1)
    return 4 * (-math.sin(theta) * pr + math.cos(theta) * pz) / p - math.sin(theta) / rho


def shoot(a, z0, events, alpha=(1.0, 1.0), dense=False):
    """The minimal surface that leaves the axis at z0, until the first of `events` fires."""
    def rhs(s, y):
        return [math.cos(y[2]), math.sin(y[2]), turning(a, y[0], y[1], y[2], alpha)]

    def lost(s, y):
        return FAR - math.hypot(y[0], y[1])
    lost.terminal = True
    for e in events:
        e.terminal = True
    return solve_ivp(rhs, [0, 4 * FAR], [0.0, z0, 0.0], events=[*events, lost], rtol=RTOL, atol=ATOL,
                     dense_output=dense, first_step=1e-6)


def throat_miss(a, z0, alpha=(1.0, 1.0)):
    """How far the curve from z0 misses closing on the axis. A curve that bends too much turns to
    point straight back at the axis while still off it, and its distance from the axis there is the
    miss, positive. A curve that bends too little is turned away by the axis before it points
    straight at it, and how far its tangent was from doing so is the miss, negative. Zero for a
    smooth closed surface."""
    def level(s, y):
        return y[2] + math.pi

    def away(s, y):
        return turning(a, y[0], y[1], y[2], alpha) if s > 1e-3 else -1.0
    away.direction = 1
    sol = shoot(a, z0, [level, away], alpha)
    if len(sol.t_events[0]):
        return float(sol.y_events[0][0][0])
    if len(sol.t_events[1]):
        return -(float(sol.y_events[1][0][2]) + math.pi)
    return float("nan")


def common_miss(a, z0, alpha=(1.0, 1.0)):
    """For two equal holes: the cosine of the tangent's angle where the curve from z0 crosses the
    midplane z = 0, which is zero for a surface that the reflection z -> -z closes."""
    def midplane(s, y):
        return y[1] if y[0] > 1e-6 else 1.0
    sol = shoot(a, z0, [midplane], alpha)
    if not len(sol.t_events[0]):
        return float("nan")
    return math.cos(float(sol.y_events[0][0][2]))


def _roots(miss, lo, hi, n):
    zs = np.linspace(lo, hi, n)
    v = [miss(z) for z in zs]
    out = []
    for i in range(n - 1):
        if math.isfinite(v[i]) and math.isfinite(v[i + 1]) and v[i] * v[i + 1] < 0:
            out.append(brentq(miss, zs[i], zs[i + 1], xtol=1e-13, rtol=1e-14))
    return out


def throat(a, alpha=(1.0, 1.0)):
    """The top z_0 of the throat of the hole at z = +a, the minimal surface around that hole
    alone: going up the axis from the hole, the first point whose curve stops bending too much.
    Alone, the hole's throat is the sphere of radius alpha about it."""
    def miss(z):
        return throat_miss(a, z, alpha)
    zs = np.linspace(a + 0.02 * alpha[0], a + 1.5 * alpha[0], 75)
    before = miss(zs[0])
    for lo, hi in zip(zs, zs[1:]):
        after = miss(hi)
        if before > 0 and after < 0:
            return brentq(miss, lo, hi, xtol=1e-13, rtol=1e-14)
        before = after
    raise AssertionError(f"no throat found around the hole at z = {a}")


def deepest(a):
    """Where on the axis above the holes common_miss is least, and its value there. The miss is
    positive far up the axis and negative only between the two surfaces that surround both holes,
    an inner and an outer one, so they exist exactly when the least value is negative."""
    zs = np.linspace(a + 0.3, a + 2.5, 45)
    v = np.array([common_miss(a, z) for z in zs])
    v[~np.isfinite(v)] = np.inf
    k = int(np.argmin(v))
    found = minimize_scalar(lambda z: common_miss(a, z), bounds=(zs[max(k - 1, 0)], zs[min(k + 1, len(zs) - 1)]),
                            method="bounded", options={"xatol": 1e-12})
    return float(found.x), float(found.fun)


def common(a):
    """The top z_0 of the outermost minimal surface around two equal holes of alpha = 1 at
    z = +-a, the apparent horizon, or None where there is none."""
    z, least = deepest(a)
    if least >= 0:
        return None
    return brentq(lambda z: common_miss(a, z), z, a + 2.5, xtol=1e-13, rtol=1e-14)


def curve(a, z0, kind, alpha=(1.0, 1.0)):
    """The curve from z0 as a dense solution with its end: kind "throat" runs to the axis, and
    "common" to the midplane."""
    if kind == "throat":
        def stop(s, y):
            return y[0] - 10 * AXIS if s > 1e-3 else 1.0
    else:
        def stop(s, y):
            return y[1] if y[0] > 1e-6 else 1.0
    sol = shoot(a, z0, [stop], alpha, dense=True)
    return sol.sol, float(sol.t_events[0][0])


def crossing(a, z0, alpha=(1.0, 1.0)):
    """Where the throat from z0 cuts the coordinate sphere rho^2 + z^2 = a^2 through both holes:
    (rho, z) of the circle."""
    path, end = curve(a, z0, "throat", alpha)
    s = brentq(lambda s: math.hypot(*path(s)[:2]) - a, 1e-9, end, xtol=1e-14)
    rho, z, _ = path(s)
    return float(rho), float(z)


def waist(a):
    """rho of the circle in which the apparent horizon of two equal holes cuts the midplane, or
    None where no surface surrounds both."""
    z0 = common(a)
    if z0 is None:
        return None
    path, end = curve(a, z0, "common")
    return float(path(end)[0])


def common_crossing(a):
    """Where the apparent horizon of two equal holes cuts the coordinate sphere rho^2 + z^2 = a^2
    above the midplane, (rho, z), or None where there is no such surface or the sphere lies inside
    it. Its waist is its nearest point to the centre, so it cuts the sphere exactly while the waist
    is narrower than a, which is for a between 1.4126 and the critical 1.5324."""
    z0 = common(a)
    if z0 is None:
        return None
    path, end = curve(a, z0, "common")
    if math.hypot(*path(end)[:2]) >= a:
        return None
    s = brentq(lambda s: math.hypot(*path(s)[:2]) - a, 1e-9, end, xtol=1e-14)
    rho, z, _ = path(s)
    return float(rho), float(z)


def area(a, z0, kind, alpha=(1.0, 1.0), n=4001):
    """The area of the minimal surface from z0, 2 pi times the integral of psi^4 rho ds, the
    whole of it: a common surface is doubled across the midplane."""
    path, end = curve(a, z0, kind, alpha)
    s = np.linspace(0, end, n)
    y = path(s)
    f = np.array([psi(a, r, z, alpha)[0] ** 4 * r for r, z in zip(y[0], y[1])])
    h = s[1] - s[0]
    total = h / 3 * (f[0] + f[-1] + 4 * f[1:-1:2].sum() + 2 * f[2:-1:2].sum())
    return 2 * math.pi * total * (2 if kind == "common" else 1)


def critical():
    """The half separation a, in units of alpha, below which a minimal surface surrounds two equal
    holes: 1.5319, a coordinate separation of 1.532 in units of 2 alpha."""
    return brentq(lambda a: deepest(a)[1], 1.45, 1.6, xtol=1e-9)


def masses(a, alpha=(1.0, 1.0)):
    """Brill and Lindquist's masses, as lengths GM/c^2: the total, read at the infinity of the
    sheet both holes share, and each hole's own, read at the infinity beyond its throat."""
    total = 2 * (alpha[0] + alpha[1])
    cross = 2 * alpha[0] * alpha[1] / (2 * a)
    return total, 2 * alpha[0] + cross, 2 * alpha[1] + cross
