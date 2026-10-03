"""Herdeiro and Radu's Kerr black hole with scalar hair, as every drawing of it is made: their
configuration IV, read from the data they published with their paper of 2015.

kerr_scalar_hair_IV.dat.gz beside this file is the authors' configuration-IV.dat, gzipped and
otherwise as published with C. Herdeiro and E. Radu, Class. Quantum Grav. 32, 144001 (2015),
arXiv:1501.04319, at http://gravitation.web.ua.pt/node/416 (the paper's reference [gravwebsite]);
its SHA-256 before gzipping is cdf37db6c7eec5fcbef74158128cba4bc44ded4fa041d4664b3be493c9703b73.
It tabulates F_1, F_2, F_0, phi and W of the paper's line element (2.5) on 251 points of
X = x/(1 + x), x = sqrt(r^2 - r_H^2), from 0 to 1, by 30 of theta from 0 to pi/2, in units with
G = c = 1 and the field's inverse length mu = 1. The field is scaled by their (3.21) so that the
8 pi G phi^2 of their equations is 2 phi^2: the tabulated phi is the dimensionless amplitude
sigma = sqrt(4 pi G) phi/c^2. The solution has m = 1, w = Omega_H = 0.82 and r_H = 0.1, and the
authors' page gives its mass 0.933, angular momentum 0.739, and horizon mass and angular momentum
0.234 and 0.114.

Hair evaluates the five functions and their first two derivatives along r and theta anywhere
outside the horizon by bicubic splines in X and theta, the angle reflected through the axis and
through the equator by each function's parity, and computes from them what the drawings and the
tests use: the mass and angular momentum from the fall of g_tt and g_tvarphi, the horizon's angular
velocity, temperature and area, the Noether charge and the field's energy outside the horizon, the
ergosurface, the tortoise coordinate along the axis and the embedding of the equatorial plane.
kerr_hr gives Kerr's metric functions in the same chart, the paper's (A.1), and kerr_for the two
constants of the Kerr black hole of a given mass and angular momentum. kerr_scalar_hair.md beside
this file is the derivation.
"""
import gzip
import math
from functools import cached_property
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.interpolate import RectBivariateSpline
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
DATA = HERE / "kerr_scalar_hair_IV.dat.gz"

R_H = 0.1          # the horizon radius of the chart, in units of 1/mu
OMEGA = 0.82       # the field's frequency w = m Omega_H, in units of mu c
M_WIND = 1         # the winding number m
# The authors' page: ADM mass and angular momentum, and their horizon parts.
PUBLISHED = {"M": 0.933, "J": 0.739, "M_H": 0.234, "J_H": 0.114}
NAMES = ("F_1", "F_2", "F_0", "phi", "W")
# Parity through the axis theta -> -theta; every function is even through the equator.
AXIS_PARITY = {"F_1": 1, "F_2": 1, "F_0": 1, "phi": -1, "W": 1}


def table():
    """The authors' grid: X (251), theta (30) and the five functions as arrays [theta, X]."""
    with gzip.open(DATA, "rt") as handle:
        rows = np.loadtxt(handle)
    X = np.unique(rows[:, 0])
    theta = np.unique(rows[:, 1])
    grid = rows.reshape(len(theta), len(X), rows.shape[1])
    if not (np.allclose(grid[0, :, 0], X) and np.allclose(grid[:, 0, 1], theta)):
        raise ValueError("the data is not ordered X fastest within theta")
    return X, theta, {name: grid[:, :, 2 + k] for k, name in enumerate(NAMES)}


class Hair:
    """Configuration IV: the five functions and their derivatives along r and theta."""

    def __init__(self):
        self.rH, self.w, self.m = R_H, OMEGA, M_WIND
        X, theta, values = table()
        self.X, self.theta, self.values = X, theta, values
        # theta reflected to [-pi/2, pi]: through the axis by each function's parity, through the
        # equator evenly, so that the splines see a smooth periodic function at both ends.
        below, above = -theta[:0:-1], np.pi - theta[-2::-1]
        angles = np.concatenate([below, theta, above])
        self.splines = {}
        for name, v in values.items():
            whole = np.concatenate([AXIS_PARITY[name] * v[:0:-1], v, v[-2::-1]], axis=0)
            self.splines[name] = RectBivariateSpline(angles, X, whole, kx=3, ky=3, s=0)

    def radius(self, X):
        x = np.asarray(X, dtype=float) / (1 - np.asarray(X, dtype=float))
        return np.sqrt(x * x + self.rH ** 2)

    def compact(self, r):
        x = np.sqrt(np.asarray(r, dtype=float) ** 2 - self.rH ** 2)
        return x / (1 + x)

    def __call__(self, name, r, theta, dr=0, dtheta=0):
        """The function `name` at (r, theta), or its derivative dr times along r and dtheta times along
        theta, dr + dtheta at most 2. The first derivative along r holds on the horizon too, by its limit
        there; the second only off it, where X'(r) is finite."""
        r = np.asarray(r, dtype=float)
        theta = np.asarray(theta, dtype=float)
        x = np.sqrt(r * r - self.rH ** 2)
        X = x / (1 + x)
        s = self.splines[name]
        if dr == 0:
            return s.ev(theta, X, dx=dtheta, dy=0)
        with np.errstate(all="ignore"):
            Xr = r / (x * (1 + x) ** 2)
        if dr == 1:
            # On the horizon the function is even in x, so its slope along X vanishes as X does and
            # dF/dr tends to r_H d^2F/dX^2 there.
            on = x < 1e-9
            return np.where(on, self.rH * s.ev(theta, X, dx=dtheta, dy=2), s.ev(theta, X, dx=dtheta, dy=1) * np.where(on, 1.0, Xr))
        Xrr = -self.rH ** 2 / (x ** 3 * (1 + x) ** 2) - 2 * r * r / (x * x * (1 + x) ** 3)
        return s.ev(theta, X, dy=2) * Xr ** 2 + s.ev(theta, X, dy=1) * Xrr

    def horizon(self, name, theta):
        """A function on the horizon, X = 0."""
        return self.splines[name].ev(np.asarray(theta, dtype=float), np.zeros_like(np.asarray(theta, dtype=float)))

    # -- the metric ------------------------------------------------------------------------

    def metric(self, r, theta):
        """g_tt, g_tvarphi, g_varphivarphi, g_rr and g_thetatheta at (r, theta), in units G = c = mu = 1."""
        F0, F1, F2, W = (self(n, r, theta) for n in ("F_0", "F_1", "F_2", "W"))
        N = 1 - self.rH / np.asarray(r, dtype=float)
        rs = np.asarray(r, dtype=float) * np.sin(theta)
        gpp = np.exp(2 * F2) * rs ** 2
        return (-np.exp(2 * F0) * N + gpp * W ** 2, -gpp * W, gpp, np.exp(2 * F1) / N, np.exp(2 * F1) * r ** 2)

    # -- charges at infinity ---------------------------------------------------------------

    @cached_property
    def charges(self):
        """M and J from F_0 = c_t/r + ... and W = c_phi/r^3 + ... on the equator, the paper's (B.3):
        2GM/c^2 = r_H - 2c_t and 2GJ/c^3 = c_phi, fitted over the points 0.95 < X < 0.999, where the
        field has fallen below 10^-4 of its greatest value and only powers of 1/r are left."""
        pick = (self.X > 0.95) & (self.X < 0.999)
        r = self.radius(self.X[pick])
        F0 = self.values["F_0"][-1, pick]
        W = self.values["W"][-1, pick]
        ct = np.polyfit(1 / r, F0, 4)[-2]
        cphi = np.polyfit(1 / r, W * r ** 3, 3)[-1]
        return {"M": (self.rH - 2 * ct) / 2, "J": cphi / 2, "c_t": ct}

    # -- the horizon -----------------------------------------------------------------------

    @cached_property
    def horizon_quantities(self):
        """Omega_H, T_H, A_H and the spread of F_0 - F_1 over the horizon, from the tabulated row X = 0:
        T_H = e^(F_0 - F_1)/(4 pi r_H) and A_H = 2 pi r_H^2 int sin(theta) e^(F_1 + F_2) dtheta over the
        sphere, the letter's formulas (with F^(0), the horizon values; the paper of 2015 prints F^(2))."""
        th = self.theta
        F0, F1, F2, W = (self.values[n][:, 0] for n in ("F_0", "F_1", "F_2", "W"))
        gap = F0 - F1
        area = 2 * (2 * np.pi * self.rH ** 2 * simpson(np.sin(th) * np.exp(F1 + F2), th))
        return {"Omega_H": float(W.mean()), "Omega_spread": float(np.ptp(W)), "T_H": float(np.exp(gap.mean()) / (4 * np.pi * self.rH)),
                "gap_spread": float(np.ptp(gap)), "A_H": float(area)}

    # -- the field outside the horizon -----------------------------------------------------

    def _volume(self, integrand):
        """The integral over r > r_H and the whole sphere of an expression in (r, theta, F_0, F_1, F_2,
        W, sigma, N), with the angle by Simpson's rule on the tabulated angles and r by the tabulated X."""
        X = self.X[1:-1]
        r = self.radius(X)
        x = X / (1 - X)
        drdX = (x / r) / (1 - X) ** 2
        out = []
        for k, th in enumerate(self.theta):
            v = {n: self.values[n][k, 1:-1] for n in NAMES}
            N = 1 - self.rH / r
            f = integrand(r, th, v["F_0"], v["F_1"], v["F_2"], v["W"], v["phi"], N) * drdX
            out.append(simpson(f, X))
        # Twice the half sphere, and 2 pi round the axis.
        return 2 * 2 * np.pi * simpson(np.array(out), self.theta)

    @cached_property
    def field(self):
        """The Noether charge Q and the field's energy outside the horizon M_Psi, the paper's (3.17) and
        (3.16) with 4 pi G phi^2 = sigma^2 (G = c = mu = 1), and the horizon's share by Smarr's
        relation M = 2 T_H S + 2 Omega_H (J - mQ) + M_Psi: M_H = M - M_Psi and J_H = J - mQ."""
        w, m = self.w, self.m
        # Their 4 pi int dr dtheta (...) phi^2 is int dr dtheta (...) sigma^2, the volume integral over 2 pi.
        Q = self._volume(lambda r, th, F0, F1, F2, W, s, N:
                         r ** 2 * np.sin(th) * np.exp(-F0 + 2 * F1 + F2) * m * (w - m * W) / N * s ** 2) / (2 * np.pi)
        Mpsi = -self._volume(lambda r, th, F0, F1, F2, W, s, N:
                             r ** 2 * np.sin(th) * np.exp(F0 + 2 * F1 + F2) * (1 - 2 * np.exp(-2 * F0) * w * (w - m * W) / N)
                             * s ** 2) / (2 * np.pi)
        return {"Q": Q, "M_Psi": Mpsi}

    # -- surfaces and drawings -------------------------------------------------------------

    def ergosurface(self, theta=np.pi / 2):
        """The outermost radius on the ray of angle theta at which g_tt vanishes, or None."""
        gtt = lambda r: float(self.metric(np.array(r), theta)[0])
        rs = np.geomspace(self.rH * 1.0005, 60.0, 4000)
        values = np.array([gtt(r) for r in rs])
        sign = np.where(np.diff(np.sign(values)) != 0)[0]
        if len(sign) == 0:
            return None
        k = sign[-1]
        return brentq(gtt, rs[k], rs[k + 1], xtol=1e-13)

    def lapse_ratio(self, r, theta=0.0):
        """e^(F_1 - F_0) at (r, theta), the factor in dr_*/dr = e^(F_1 - F_0)/N along the axis."""
        return np.exp(self("F_1", r, theta) - self("F_0", r, theta))

    @cached_property
    def axis_gap(self):
        """e^(F_1 - F_0) on the horizon at the axis, the coefficient of the logarithm in r_*."""
        return float(np.exp(self.horizon("F_1", 0.0) - self.horizon("F_0", 0.0)))

    def tortoise(self, r):
        """r_* on the axis, dr_*/dr = e^(F_1 - F_0) r/(r - r_H). In x = sqrt(r^2 - r_H^2) it is
        dr_*/dx = C(x) (r + r_H)/x with C = e^(F_1 - F_0), which is even in x, so
        r_* = 2 r_H C(0) ln(x/r_H) + int_0^x (C(s)(sqrt(s^2 + r_H^2) + r_H) - 2 r_H C(0)) ds/s, the
        logarithm carrying the whole divergence at the horizon and the quadrature smooth."""
        C0, rH = self.axis_gap, self.rH

        def C(x):
            return float(np.exp(self.splines["F_1"].ev(0.0, x / (1 + x)) - self.splines["F_0"].ev(0.0, x / (1 + x))))

        def smooth(x):
            return (C(x) * (math.sqrt(x * x + rH * rH) + rH) - 2 * rH * C0) / x

        def one(r):
            x = math.sqrt(r * r - rH * rH)
            return 2 * rH * C0 * math.log(x / rH) + quad(smooth, 0.0, x, limit=400, epsabs=1e-12, epsrel=1e-12)[0]
        return np.vectorize(one, otypes=[float])(np.asarray(r, dtype=float))

    @cached_property
    def _smooth_tortoise(self):
        """The regular part of tortoise(), int_0^x of its smooth integrand, tabulated on 1600 points of
        x from 10^-6 to 10^4 by one quadrature per interval, as a cubic spline in ln x."""
        from scipy.interpolate import CubicSpline
        C0, rH = self.axis_gap, self.rH
        xs = np.geomspace(1e-6, 1e4, 1600)

        def smooth(x):
            X = x / (1 + x)
            C = float(np.exp(self.splines["F_1"].ev(0.0, X) - self.splines["F_0"].ev(0.0, X)))
            return (C * (math.sqrt(x * x + rH * rH) + rH) - 2 * rH * C0) / x
        steps = [quad(smooth, 0.0, xs[0], epsabs=1e-14)[0]]
        steps += [quad(smooth, a, b, limit=200, epsabs=1e-13, epsrel=1e-12)[0] for a, b in zip(xs[:-1], xs[1:])]
        return CubicSpline(np.log(xs), np.cumsum(steps))

    def tortoise_fast(self, r):
        """tortoise() from the table, for the many points a conformal map is drawn and checked at."""
        r = np.asarray(r, dtype=float)
        x = np.sqrt(np.maximum(r * r - self.rH ** 2, 1e-300))
        smooth = np.where(x > 1e-6, self._smooth_tortoise(np.log(np.maximum(x, 1e-6))), 0.0)
        return 2 * self.rH * self.axis_gap * np.log(x / self.rH) + smooth

    def circumference_radius(self, r, theta=np.pi / 2):
        """e^(F_2) r sin(theta), the radius of the circle of constant t, r and theta."""
        return np.exp(self("F_2", r, theta)) * np.asarray(r, dtype=float) * np.sin(theta)


def field_equations(hair, rmin=0.15, rmax=30.0):
    """The residuals of the field equations; the axis and the horizon divide by zero and lie outside
    the points kept, so numpy is told not to warn of it."""
    with np.errstate(all="ignore"):
        return _field_equations(hair, rmin, rmax)


def _field_equations(hair, rmin, rmax):
    """The paper's (2.10) to (2.13) on the authors' grid: the Klein-Gordon equation, the equations for
    F_1, F_2, F_0 and W drawn from the combinations (2.11), and the two constraints, which print_charts.kerr_scalar_hair_check holds to the
    published Einstein tensor of the free chart, each times one factor. Derivatives along X and theta
    are fourth order central differences on the tabulated points, the angle continued through the axis
    and the equator by each function's parity, and turned into derivatives along r through
    X = x/(1 + x), x = sqrt(r^2 - r_H^2). Returns, for each equation, the largest residual over the
    points with rmin < r < rmax and theta two grid steps or more from the axis, and the largest
    single term there, its scale."""
    X, th, values = hair.X, hair.theta, hair.values
    hX, ht = X[1] - X[0], th[1] - th[0]
    T = np.concatenate([-th[:0:-1], th, np.pi - th[-2::-1]])
    d = {}
    for name, v in values.items():
        f = np.concatenate([AXIS_PARITY[name] * v[:0:-1], v, v[-2::-1]], axis=0)
        fX, fXX, fT, fTT = (np.full_like(f, np.nan) for _ in range(4))
        fX[:, 2:-2] = (f[:, :-4] - 8 * f[:, 1:-3] + 8 * f[:, 3:-1] - f[:, 4:]) / (12 * hX)
        fXX[:, 2:-2] = (-f[:, :-4] + 16 * f[:, 1:-3] - 30 * f[:, 2:-2] + 16 * f[:, 3:-1] - f[:, 4:]) / (12 * hX ** 2)
        fT[2:-2] = (f[:-4] - 8 * f[1:-3] + 8 * f[3:-1] - f[4:]) / (12 * ht)
        fTT[2:-2] = (-f[:-4] + 16 * f[1:-3] - 30 * f[2:-2] + 16 * f[3:-1] - f[4:]) / (12 * ht ** 2)
        fXT = np.full_like(f, np.nan)
        fXT[2:-2] = (fX[:-4] - 8 * fX[1:-3] + 8 * fX[3:-1] - fX[4:]) / (12 * ht)
        d[name] = (f, fX, fXX, fT, fTT, fXT)
    XX, TT = np.meshgrid(X, T)
    with np.errstate(all="ignore"):
        x = XX / (1 - XX)
        r = np.sqrt(x * x + hair.rH ** 2)
        Xr = r / (x * (1 + x) ** 2)
        Xrr = -hair.rH ** 2 / (x ** 3 * (1 + x) ** 2) - 2 * r * r / (x * x * (1 + x) ** 3)
    D = {}
    for name, (f, fX, fXX, fT, fTT, fXT) in d.items():
        D[name] = {"": f, "r": fX * Xr, "rr": fXX * Xr ** 2 + fX * Xrr, "t": fT, "tt": fTT, "rt": fXT * Xr}
    F0, F1, F2, W, s = (D[n] for n in ("F_0", "F_1", "F_2", "W", "phi"))
    w, m, mu = hair.w, hair.m, 1.0
    N = 1 - hair.rH / r
    Np = hair.rH / r ** 2
    sn, ct = np.sin(TT), np.cos(TT) / np.sin(TT)
    a = 1 / (r * r * N)
    grad = lambda f, g: f["r"] * g["r"] + a * f["t"] * g["t"]
    minus = lambda f, g: f["r"] * g["r"] - a * f["t"] * g["t"]
    lap = lambda f: f["rr"] + a * f["tt"]
    sum02 = {k: F0[k] + F2[k] for k in ("r", "t")}
    twist = np.exp(-2 * F0[""] + 2 * F2[""]) * r * r * sn * sn / N
    e2F1 = np.exp(2 * F1[""])
    beat = (w - m * W[""]) ** 2
    k2 = 2 * s[""] ** 2
    kg = (lap(s) + s["r"] * sum02["r"] + a * s["t"] * sum02["t"] + (1 + r * Np / (2 * N)) * 2 / r * s["r"] + ct * a * s["t"]
          - (np.exp(-2 * F2[""]) * m * m / (r * r * sn * sn) - np.exp(-2 * F0[""]) * beat / N + mu * mu) * e2F1 / N * s[""])
    eF1 = (lap(F1) - grad(F0, F2) - twist / 4 * grad(W, W) - F0["r"] / r - Np * F2["r"] / (2 * N)
           + (1 + r * Np / (2 * N)) * F1["r"] / r - ct * a * F0["t"]
           + 2 * (s["r"] ** 2 + a * s["t"] ** 2 + e2F1 / N ** 2 * (np.exp(-2 * F0[""]) * beat
                                                                - np.exp(-2 * F2[""]) * m * m * N / (r * r * sn * sn)) * s[""] ** 2))
    eF2 = (lap(F2) + grad(F2, F2) + grad(F0, F2) + twist / 2 * grad(W, W) + (F0["r"] + ct * F0["t"] / (r * N)) / r
           + (1 + r * Np / (3 * N)) * 3 * F2["r"] / r + 2 * ct * a * F2["t"]
           + k2 * e2F1 / N * (mu * mu + 2 * np.exp(-2 * F2[""]) * m * m / (r * r * sn * sn)))
    eF0 = (lap(F0) + grad(F0, F0) + grad(F0, F2) - twist / 2 * grad(W, W) + (1 + 3 * r * Np / (4 * N)) * 2 * F0["r"] / r
           + ct * a * F0["t"] + Np * F2["r"] / (2 * N) - k2 * e2F1 / N * (2 * np.exp(-2 * F0[""]) * beat / N - mu * mu))
    eW = (lap(W) + (3 * F2["r"] - F0["r"]) * W["r"] + a * (3 * F2["t"] - F0["t"]) * W["t"]
          + 4 / r * (W["r"] + 3 * ct * W["t"] / (4 * r * N)) + 4 * k2 * np.exp(2 * F1[""] - 2 * F2[""]) * m * (w - m * W[""]) / (r * r * sn * sn * N))
    c1 = (F0["rr"] - a * F0["tt"] + F2["rr"] - a * F2["tt"] + minus(F0, F0) - 2 * minus(F0, F1) - 2 * minus(F1, F2)
          - twist / 2 * minus(W, W) + minus(F2, F2) + (3 * r * Np / (2 * N) - 1) * F0["r"] / r
          + (1 + r * Np / (2 * N)) * (F2["r"] - 2 * F1["r"]) / r + 2 * ct * a * (F1["t"] - F2["t"]) + 4 * minus(s, s))
    c2 = (F0["rt"] + F2["rt"] + F0["r"] * F0["t"] + F2["r"] * F2["t"] - (F0["r"] * F1["t"] + F1["r"] * F0["t"])
          - (F1["r"] * F2["t"] + F2["r"] * F1["t"]) + (r * Np / (2 * N) - 1) * F0["t"] / r - (1 + r * Np / (2 * N)) * F1["t"] / r
          - ct * (F1["r"] - F2["r"]) - twist / 2 * W["r"] * W["t"] + 4 * s["r"] * s["t"])
    scale = {"Klein-Gordon": np.abs(s["rr"]), "F_1": np.abs(F1["rr"]), "F_2": np.abs(F2["rr"]), "F_0": np.abs(F0["rr"]),
             "W": np.abs(W["rr"]), "E^r_r - E^theta_theta": np.abs(F0["rr"]), "E^r_theta": np.abs(F0["rt"])}
    keep = (r > rmin) & (r < rmax) & (TT >= th[2]) & (TT <= np.pi / 2 + 1e-12)
    out = {}
    for name, value in (("Klein-Gordon", kg), ("F_1", eF1), ("F_2", eF2), ("F_0", eF0), ("W", eW),
                        ("E^r_r - E^theta_theta", c1), ("E^r_theta", c2)):
        out[name] = (float(np.nanmax(np.abs(value[keep]))), float(np.nanmax(scale[name][keep])))
    return out


def simpson(f, x):
    """Simpson's rule on the points x, which are equally spaced, by scipy's composite rule."""
    from scipy.integrate import simpson as integrate
    return integrate(f, x=x)


def kerr_hr(rH, ct, r, theta):
    """Kerr's F_0, F_1, F_2 and W in Herdeiro and Radu's chart, the paper's (A.1), c_t < 0."""
    e2F1 = (1 - ct / r) ** 2 + ct * (ct - rH) * np.cos(theta) ** 2 / r ** 2
    e2F2 = (((1 - ct / r) ** 2 + ct * (ct - rH) / r ** 2) ** 2 + ct * (rH - ct) * (1 - rH / r) * np.sin(theta) ** 2 / r ** 2) / e2F1
    W = np.sqrt(ct * (ct - rH)) * (rH - 2 * ct) * (1 - ct / r) / (r ** 3 * e2F1 * e2F2)
    return -np.log(e2F2) / 2, np.log(e2F1) / 2, np.log(e2F2) / 2, W


def kerr_for(M, J):
    """r_H and c_t of the Kerr black hole of mass M and angular momentum J (G = c = 1): the paper's
    (A.2) inverted, r_H = 2 sqrt(M^2 - a^2) and c_t = -(M - sqrt(M^2 - a^2)), a = J/M."""
    a = J / M
    root = math.sqrt(M * M - a * a)
    return 2 * root, -(M - root)
