#!/usr/bin/env python3
"""Ori's exact model of mass inflation, as numbers, in units of the final mass m_0, G = c = 1.

A charged black hole of charge q takes in a stream of null dust, the tail of its collapse, and
one thin shell of null dust runs outward inside it. On each side of the shell the metric is the
charged Vaidya metric of mass_inflation's ingoing chart (A. Ori, Phys. Rev. Lett. 67, 789, 1991,
as A. Bonanno, S. Droz, W. Israel and S. M. Morsink restate it, Proc. R. Soc. Lond. A 450, 553,
1995, their (3) to (9)),

    ds^2 = -f dv^2 + 2 dv dr + r^2 dOmega^2,      f = 1 - 2 m(v)/r + q^2/r^2,

with its own advanced time and its own mass function: v_1 and m_1 before the shell, v_2 and m_2
behind it. Before the shell the mass is the declared one, Price's tail,

    m_1(v_1) = m_0 - dm (v_0/v_1)^(p - 1),        p = 12,

from v_0 on, so that the influx dm_1/dv_1 falls off as v_1^(-p), and the constant m_0 - dm before
v_0, a Reissner-Nordstrom hole that the tail has not reached yet.

The shell is the outgoing light ray r = R, so on each side dR/dv_i = f_i/2, Bonanno's (5). Let
lambda be an affine parameter along it and z_i = R/(dv_i/dlambda). The geodesic equation of each
side, v'' = -(m_i/R^2 - q^2/R^3) v'^2, gives

    dz_i/dlambda = (1 - q^2/R^2)/2

on both sides, so z_2 - z_1 is a constant, Z, and R' = f_i R/(2 z_i) gives the mass behind,

    m_2 = m_1 - Z dR/dlambda,

the shell's own mass being -Z dR/dlambda. The influx is then continuous across the shell with
nothing more assumed: Raychaudhuri's equation on each side is R'' = -(dm_i/dlambda)/z_i, and the
two agree. Bonanno's (7), dm_2/f_2 = dm_1/f_1, follows, and consistent() checks it on the numbers.

Everything is integrated in v_1, which runs to infinity at the Cauchy horizon, with
x = R - r_-(v_1) in place of R so that f_1 = (R - r_+) x/R^2 keeps its digits as x -> 0:

    dx/dv_1 = f_1/2 - dr_-/dv_1,   dln z_1/dv_1 = (1 - q^2/R^2)/(2R),   dv_2/dv_1 = z_1/(z_1 + Z).

z_1 falls as exp(-kappa v_1), kappa the surface gravity of the inner horizon, so v_2 reaches the
Cauchy horizon at a finite value, which is taken as v_2 = 0: behind the shell the advanced time
runs from negative values up to zero, and m_2 grows as (-v_2)^(-1) ln(-1/v_2)^(-p) there, which
is m_2 ~ v_1^(-p) exp(kappa v_1), Ori's law. Its integral over v_2 converges, so an outgoing ray
behind the shell reaches the Cauchy horizon at a radius above zero.

    .venv.noindex/bin/python _tools/derivations/ori_shell.py
"""
import math

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

Q = 0.96            # the charge, so that the hole left behind is the one Reissner-Nordstrom's diagrams draw
P = 12              # Price's tail of a quadrupole: the influx falls as v^(-12)
V0 = 10.0           # the advanced time the tail is switched on at
V_START = 0.0       # where the table begins, in the static hole before the tail
DM = 0.02           # the mass still to fall in at V0
R0 = 1.0            # the shell's radius at V0, between the apparent horizons
SHELL = 0.02        # the shell's mass at V0
V_END = 85.0        # where the table ends: behind the shell the Cauchy horizon is 4e-18 m_0 of v_2 away

ROOT = math.sqrt(1 - Q * Q)
R_PLUS, R_MINUS = 1 + ROOT, 1 - ROOT        # the horizons of the hole left behind, 1.28 and 0.72
KAPPA = ROOT / R_MINUS ** 2                 # the surface gravity of its inner horizon, 0.540


def mass_before(v):
    """m_1 at the advanced time v before the shell: constant until V0, Price's tail after."""
    return 1 - DM * (V0 / np.maximum(np.asarray(v, dtype=float), V0)) ** (P - 1)


def influx_before(v):
    """dm_1/dv_1."""
    v = np.asarray(v, dtype=float)
    return np.where(v >= V0, DM * (P - 1) / V0 * (V0 / np.maximum(v, V0)) ** P, 0.0)


def horizons(m):
    """The two roots of r^2 - 2 m r + q^2, the apparent horizons at mass m."""
    m = np.asarray(m, dtype=float)
    root = np.sqrt(m * m - Q * Q)
    return m - root, m + root


def f(m, r):
    return 1 - 2 * m / r + Q * Q / r ** 2


class Split:
    """A cubic spline on each side of a knot where the function's derivative jumps, as every
    function of the shell does where the tail is switched on."""

    def __init__(self, x, y, knot):
        k = int(np.argmin(np.abs(x - knot)))
        self.knot = x[k]
        self.low, self.high = CubicSpline(x[:k + 1], y[:k + 1]), CubicSpline(x[k:], y[k:])

    def __call__(self, x, nu=0):
        x = np.asarray(x, dtype=float)
        return np.where(x < self.knot, self.low(x, nu), self.high(x, nu))


class Shell:
    """The shell from V_START to V_END, every quantity a function of v_1 on the shell."""

    def __init__(self):
        inner0, _ = horizons(mass_before(V0))
        self.Z = -SHELL / (f(mass_before(V0), R0) / 2)      # z_1 = R at V0, so lambda' = 1 there
        start = [R0 - inner0, math.log(R0), 0.0, 0.0]
        legs = [solve_ivp(self._rhs, [V0, end], start, method="DOP853", rtol=1e-13, atol=1e-40, max_step=0.05)
                for end in (V_START, V_END)]
        self.v = np.concatenate([legs[0].t[:0:-1], legs[1].t])
        x, lnz, v2, lam = np.concatenate([legs[0].y[:, :0:-1], legs[1].y], axis=1)
        m1 = mass_before(self.v)
        inner, outer = horizons(m1)
        self.x, self.R, self.z, self.lam = x, inner + x, np.exp(lnz), lam
        self.f1 = (self.R - outer) * x / self.R ** 2
        self.m2 = m1 - self.Z * self.f1 * self.R / (2 * self.z)
        self.f2 = f(self.m2, self.R)
        # v_2 measured from the Cauchy horizon: beyond the table z_1 = z_end exp(-kappa (v - V_END)),
        # which leaves z_end/(kappa Z) of v_2, and the rest is summed backward so that it keeps its digits.
        step = np.diff(self.v)
        rate = self.z / (self.z + self.Z)
        left = self.z[-1] / (KAPPA * self.Z)
        pieces = (rate[1:] + rate[:-1]) / 2 * step
        self.v2 = -(left + np.concatenate([np.cumsum(pieces[::-1])[::-1], [0.0]]))
        # The trapezium rule is good to the step squared; the dense solution refines it.
        fine = -(left + (v2[-1] - v2))
        self.v2 = np.where(np.abs(fine) > 1e-6, fine, self.v2)
        self._by_v1 = {name: Split(self.v, values, V0) for name, values in
                       (("R", self.R), ("lnz", lnz), ("lnm2", np.log(self.m2)), ("s", -np.log(-self.v2)),
                        ("lnx", np.log(x)), ("lam", lam))}
        s = -np.log(-self.v2)
        self._v1_of_s = Split(s, self.v, s[int(np.argmin(np.abs(self.v - V0)))])
        self.s_end = s[-1]

    def _rhs(self, v, y):
        x, lnz, _, _ = y
        m1 = float(mass_before(v))
        root = math.sqrt(m1 * m1 - Q * Q)
        inner, outer = m1 - root, m1 + root
        R, z = inner + x, math.exp(lnz)
        f1 = (R - outer) * x / R ** 2
        return [f1 / 2 + float(influx_before(v)) * inner / root, (1 - Q * Q / R ** 2) / (2 * R),
                z / (z + self.Z), z / R]

    def radius(self, v1):
        """The shell's radius where the ingoing ray v_1 crosses it."""
        return self._by_v1["R"](v1)

    def z1(self, v1):
        return np.exp(self._by_v1["lnz"](v1))

    def mass_behind_at(self, v1):
        """m_2 on the ingoing ray that crosses the shell at v_1."""
        return np.exp(self._by_v1["lnm2"](v1))

    def v2_at(self, v1):
        """The advanced time behind the shell of the ingoing ray v_1, negative, zero on the Cauchy horizon."""
        return -np.exp(-self._by_v1["s"](v1))

    def v1_at(self, v2):
        """The advanced time before the shell of the ingoing ray v_2 < 0. Beyond the table,
        v_2 = -z/(kappa Z) with z falling as exp(-kappa v_1)."""
        s = -np.log(-np.asarray(v2, dtype=float))
        inside = np.minimum(s, self.s_end)
        w = self._v1_of_s(inside)
        # Two steps of Newton's rule on the spline of s against v_1, so that this is its inverse to
        # rounding and the rays traced in v_1 are the rays of the metric written in v_2.
        for _ in range(2):
            w = np.clip(w - (self._by_v1["s"](w) - inside) / self._by_v1["s"](w, 1), V_START, V_END)
        return np.where(s <= self.s_end, w, V_END + (s - self.s_end) / KAPPA)

    def mass_behind(self, v2):
        """m_2 as a function of the advanced time behind the shell, v_2 < 0. Beyond the table it is
        continued by Ori's law, m_2 ~ v_1^(-p) exp(kappa v_1)."""
        v1 = self.v1_at(v2)
        inside = np.minimum(v1, V_END)
        return self.mass_behind_at(inside) * np.exp(KAPPA * (v1 - inside)) * (inside / v1) ** P

    def influx_behind(self, v2, h=1e-6):
        """dm_2/dv_2, by a centred difference in ln(-v_2)."""
        v2 = np.asarray(v2, dtype=float)
        up, down = self.mass_behind(v2 * math.exp(-h)), self.mass_behind(v2 * math.exp(h))
        return (up - down) / (v2 * (math.exp(-h) - math.exp(h)))

    def radius_behind(self, v2):
        """The shell's radius as a function of v_2 < 0."""
        return self.radius(np.minimum(self.v1_at(v2), V_END))

    def consistent(self):
        """Bonanno, Droz, Israel and Morsink's (7), dm_2/f_2 = dm_1/f_1 along the shell, on the
        table: the largest relative miss of (dm_2/dv_1) f_1 - (dm_1/dv_1) f_2."""
        v = self.v[(self.v > V0 + 0.5) & (self.v < V_END - 20)]
        dm2 = self._by_v1["lnm2"](v, 1) * self.mass_behind_at(v)
        m1 = mass_before(v)
        R = self.radius(v)
        f1, f2 = f(m1, R), f(self.mass_behind_at(v), R)
        return float(np.max(np.abs(dm2 * f1 - influx_before(v) * f2) / np.abs(influx_before(v) * f2)))


_SHELL = None


def shell():
    global _SHELL
    if _SHELL is None:
        _SHELL = Shell()
    return _SHELL


def outgoing_before(v_from, r_from, v_to):
    """The outgoing ray through (v_from, r_from) before the shell, carried to v_to: dr/dv = f_1/2."""
    sol = solve_ivp(lambda v, y: [f(float(mass_before(v)), y[0]) / 2], [v_from, v_to], [r_from],
                    method="DOP853", rtol=1e-11, atol=1e-13)
    return float(sol.y[0, -1])


def outgoing_behind(v2_from, r_from, v2_to):
    """The outgoing ray through (v2_from, r_from) behind the shell, carried to v2_to, integrated
    in s = -ln(-v_2), where dr/ds = -v_2 f_2/2 stays bounded up to the Cauchy horizon."""
    S = shell()

    def rhs(s, y):
        v2 = -math.exp(-s)
        return [-v2 * f(float(S.mass_behind(v2)), y[0]) / 2]
    sol = solve_ivp(rhs, [-math.log(-v2_from), -math.log(-v2_to)], [r_from], method="DOP853", rtol=1e-10, atol=1e-13)
    return float(sol.y[0, -1])


def rate_behind(w):
    """dv_2/dv_1 = z_1/(z_1 + Z) on the ingoing ray that crosses the shell at v_1 = w; beyond the
    table it falls as exp(-kappa w)."""
    S = shell()
    w = np.asarray(w, dtype=float)
    inside = np.minimum(w, V_END)
    # Inside the table, the derivative of the very spline v2_at reads, v_2 = -exp(-s).
    return (np.exp(-S._by_v1["s"](inside)) * S._by_v1["s"](inside, 1)) * np.exp(-KAPPA * np.maximum(w - V_END, 0.0))


def _behind(w, x, lnz):
    """From the shell's state on the ingoing ray w, its x = R - r_- and ln z_1: their rates of
    change in w, the mass behind the shell and dv_2/dw there, each a smooth function of the state,
    which a ray carries along with it so that nothing it integrates is read from a table."""
    m1 = mass_before(w)
    root = np.sqrt(m1 * m1 - Q * Q)
    inner, outer = m1 - root, m1 + root
    R, z = inner + x, np.exp(lnz)
    f1 = (R - outer) * x / R ** 2
    Z = shell().Z
    return (f1 / 2 + influx_before(w) * inner / root, (1 - Q * Q / R ** 2) / (2 * R),
            m1 - Z * f1 * R / (2 * z), z / (z + Z))


def _state(w):
    S = shell()
    return np.exp(S._by_v1["lnx"](w)), S._by_v1["lnz"](w)


def carried_behind(w, r, w_stop):
    """The radius at which the outgoing ray through each event (w, r) behind the shell crosses
    the ingoing ray w_stop, for rays that do not meet the centre on the way: dr/dw =
    (dv_2/dw) f_2/2, every ray over its own stretch, cut at V0."""
    w, r = np.broadcast_arrays(np.asarray(w, dtype=float), np.asarray(r, dtype=float))
    start = w.ravel()
    n = start.size
    x, lnz = _state(start)
    y = np.concatenate([x, lnz, r.ravel()])
    for a, b in _stretches(start, w_stop):
        span = b - a

        def rhs(s, y, a=a, span=span):
            x, lnz, rr = y[:n], y[n:2 * n], y[2 * n:]
            dx, dlnz, m2, rate = _behind(a + s * span, x, lnz)
            return np.concatenate([span * dx, span * dlnz, span * rate * f(m2, rr) / 2])
        y = solve_ivp(rhs, [0.0, 1.0], y, method="DOP853", rtol=1e-12, atol=1e-300).y[:, -1]
    return y[2 * n:].reshape(w.shape)


def born(w, r):
    """Where the outgoing ray through each event behind the shell began, the event given by the
    advanced time w its ingoing ray has before the shell, at most V_END, and its radius r: (True,
    r_a) for a ray that crosses the first ingoing ray of the table, w = V_START, at the radius r_a,
    and (False, w_s) for one that left the singularity r = 0 on the ingoing ray w_s. Behind the
    shell dr/dv_2 = f_2/2 and dv_2 = rate dw, which with d sigma = dw/(2 r^2) is

        dw/dsigma = 2 r^2,      dr/dsigma = rate (r^2 - 2 m_2 r + q^2),

    regular at r = 0, where a ray leaves with r^3 = (3 q^2/2)(v_2 - v_s), and bounded up to the
    Cauchy horizon, where rate m_2 falls as w^(-p). A first pass in sigma says which rays meet the
    centre; those that do not are then carried in w by carried_behind, which ends on the first
    ingoing ray exactly, and the others again in sigma to a tighter tolerance."""
    def rhs(_, y):
        ww, rr, x, lnz = y
        ww = max(ww, V_START)
        dx, dlnz, m2, rate = _behind(ww, x, lnz)
        step = -2 * rr * rr
        return [step, -float(rate) * (rr * rr - 2 * float(m2) * rr + Q * Q), step * float(dx), step * float(dlnz)]

    def first_ray(_, y):
        return y[0] - V_START

    def centre(_, y):
        return y[1]
    first_ray.terminal = centre.terminal = True
    first_ray.direction = centre.direction = -1

    def back(ww, rr, tol):
        x, lnz = _state(ww)
        return solve_ivp(rhs, [0.0, 1e9], [ww, rr, float(x), float(lnz)], method="DOP853", rtol=tol, atol=tol * 1e-2,
                         events=(first_ray, centre))
    w, r = np.broadcast_arrays(np.asarray(w, dtype=float), np.asarray(r, dtype=float))
    crossed, value = np.zeros(w.shape, dtype=bool), np.zeros(w.shape)
    for i in np.ndindex(w.shape):
        if w[i] <= V_START:
            crossed[i], value[i] = True, r[i]
            continue
        sol = back(w[i], r[i], 1e-9)
        if sol.t_events[0].size:
            crossed[i] = True
        elif sol.t_events[1].size:
            # The tighter pass decides for a ray that all but grazes the corner where the first
            # ingoing ray meets the centre.
            sol = back(w[i], r[i], 1e-13)
            if sol.t_events[1].size:
                value[i] = sol.y_events[1][0][0]
            else:
                crossed[i] = True
        else:
            raise AssertionError(f"the outgoing ray through w = {w[i]}, r = {r[i]} behind the shell was not traced back")
    later = crossed & (w > V_START)
    if later.any():
        value[later] = carried_behind(w[later], r[later], V_START)
    return crossed, value


_HORIZON = None


def event_horizon(v):
    """The radius of the event horizon at the advanced time v before the shell: the outgoing ray
    that ends on r_+ of the hole left behind, carried back from v = 400, where the mass is m_0 to
    one part in 10^19. Going back in v the rays inside and outside close on it, so the integration
    is stable."""
    global _HORIZON
    if _HORIZON is None:
        _HORIZON = solve_ivp(lambda v, y: [f(float(mass_before(v)), y[0]) / 2], [400.0, V_START], [R_PLUS],
                             method="DOP853", rtol=1e-13, atol=1e-15, dense_output=True).sol
    v = np.asarray(v, dtype=float)
    return np.where(v >= 400.0, R_PLUS, _HORIZON(np.clip(v, V_START, 400.0))[0])


def offset_before(v, r, v_ref):
    """How far outside the event horizon the outgoing ray through each event (v, r) before the
    shell crosses the ingoing ray v_ref, negative inside it. With d = r - r_H the rays obey

        dd/dv = (f(r_H + d) - f(r_H))/2 = d (m/(r_H r) - q^2 (r_H + r)/(2 r_H^2 r^2)),    r = r_H + d,

    which keeps the digits of d however small it is, as it is for every late event of the
    exterior and of the hole, whose rays left the horizon a hair to one side of it."""
    v, r = np.broadcast_arrays(np.asarray(v, dtype=float), np.asarray(r, dtype=float))
    start = v.ravel()
    n = start.size
    rh0 = event_horizon(start)
    y = np.concatenate([rh0, r.ravel() - rh0])
    for a, b in _stretches(start, v_ref):
        span = b - a

        def rhs(s, y, a=a, span=span):
            # The horizon's own radius is carried along each ray's stretch beside d, so that it is a
            # ray to the tolerance of this integration at every step.
            rh, d = y[:n], y[n:]
            m = mass_before(a + s * span)
            rr = rh + d
            return np.concatenate([span * f(m, rh) / 2,
                                   span * d * (m / (rh * rr) - Q * Q * (rh + rr) / (2 * rh ** 2 * rr ** 2))])
        y = solve_ivp(rhs, [0.0, 1.0], y, method="DOP853", rtol=1e-12, atol=1e-300).y[:, -1]
    return y[n:].reshape(v.shape)


def _stretches(start, v_ref):
    """Each ray's stretch from its event to v_ref, cut at V0, where the influx jumps: the two
    pairs (from, to), the first empty for an event that does not have V0 between it and v_ref."""
    mid = np.where((start - V0) * (v_ref - V0) < 0, V0, start)
    return (start, mid), (mid, np.full_like(start, v_ref))


def label_before(v, r, v_ref):
    """The radius at which the outgoing ray through each event (v, r) before the shell crosses the
    ingoing ray v_ref: every ray carried along dr/dv = f_1/2 at once, each over its own stretch."""
    v, r = np.broadcast_arrays(np.asarray(v, dtype=float), np.asarray(r, dtype=float))
    y = r.ravel().copy()
    for a, b in _stretches(v.ravel(), v_ref):
        span = b - a
        y = solve_ivp(lambda s, y, a=a, span=span: span * f(mass_before(a + s * span), y) / 2, [0.0, 1.0], y,
                      method="DOP853", rtol=1e-12, atol=1e-14).y[:, -1]
    return y.reshape(v.shape)


def label_behind(v2, r, v2_ref):
    """The same behind the shell, in s = -ln(-v_2)."""
    S = shell()
    v2, r = np.broadcast_arrays(np.asarray(v2, dtype=float), np.asarray(r, dtype=float))
    start = -np.log(-v2.ravel())
    span = -math.log(-v2_ref) - start

    def rhs(t, y):
        at = -np.exp(-(start + t * span))
        return -span * at * f(S.mass_behind(at), y) / 2
    sol = solve_ivp(rhs, [0.0, 1.0], r.ravel().copy(), method="DOP853", rtol=1e-11, atol=1e-13)
    return sol.y[:, -1].reshape(v2.shape)


if __name__ == "__main__":
    S = shell()
    print(f"Z = {S.Z:.6f}, kappa = {KAPPA:.6f}, r_- = {R_MINUS:.4f}, r_+ = {R_PLUS:.4f}")
    print(f"v_2 at V0 = {float(S.v2_at(V0)):.6f}")
    print(f"dm_2/f_2 = dm_1/f_1 to {S.consistent():.2e}")
    for v in (10, 12, 16, 20, 30, 40, 50, 60, 70, 80):
        print(f"v_1 = {v:3d}: R = {float(S.radius(v)):.12f}, z_1 = {float(S.z1(v)):.3e}, v_2 = {float(S.v2_at(v)):.3e}, "
              f"m_2 = {float(S.mass_behind_at(v)):.6g}, m_2 v^p e^(-kappa v) = "
              f"{float(S.mass_behind_at(v)) * v ** P * math.exp(-KAPPA * v):.4g}")
