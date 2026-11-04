#!/usr/bin/env python3
"""Mukhanov and Brandenberger's closed universe with limiting curvature, as numbers, in units
with c = H_0 = 1, so that times are in 1/H_0 and lengths in c/H_0.

The theory (V. Mukhanov and R. Brandenberger, Phys. Rev. Lett. 68, 1969, 1992; R. Brandenberger,
V. Mukhanov and A. Sornborger, Phys. Rev. D 48, 1629, 1993) has the action

    sqrt(-g) [R + phi_1 I_1 - phi_2 sqrt(I_2) - V_1(phi_1) - V_2(phi_2)],
    I_1 = R - sqrt(3) sqrt(I_2),        I_2 = 4 R_mu_nu R^mu_nu - R^2,

in the signature (-,+,+,+) with the standard contraction for the Ricci tensor. On the closed
Robertson-Walker metric ds^2 = -dt^2 + a^2 (dchi^2 + sin^2 chi dOmega^2), with H = (da/dt)/a and
k = 1, I_2 = 12 (dH/dt - k/a^2)^2, and the root is the signed one, sqrt(I_2) = sqrt(12)(dH/dt -
k/a^2); then I_1 = 12 (H^2 + k/a^2). limiting_curvature.md beside this file derives that sign
from the action, and _tools/test_limiting_curvature.py holds it. Varying phi_1, phi_2 and the
lapse gives the paper's (2.25) and (2.26) and the Hamiltonian equation

    H^2 + k/a^2 = V_1'(phi_1)/12,                       (2.25)
    dH/dt - k/a^2 = -V_2'(phi_2)/sqrt(12),               (2.26)
    -(V_1 + V_2)/2 + 3H^2 (1 - 2 phi_1) + 3 (k/a^2)(2 phi_1 + 1)
        = sqrt(3) H (dphi_2/dt + 3H phi_2 - k phi_2/(H a^2)),

the last with the coefficient 2 phi_1 + 1 where the paper prints 4 phi_1 + 1 in its (2.27) and
4 phi_1 - 1 in its (4.2). Differentiating (2.25) and using (2.26) gives its (4.1),

    dphi_1/dt = -4 sqrt(3) H V_2'(phi_2)/V_1''(phi_1).

The Hamiltonian equation divides by H, which vanishes at a bounce, so the integration uses the
equation of the scale factor in its place, the variation of a, which is second order in phi_2 and
divides by nothing:

    d^2phi_2/dt^2 = -4 sqrt(3) H dphi_1/dt - 6 H dphi_2/dt
                    + phi_1 (-6 sqrt(3) H^2 - 4 sqrt(3) dH/dt + 2 sqrt(3) k/a^2)
                    + phi_2 (-9 H^2 - 6 dH/dt + k/a^2)
                    + 3 sqrt(3) H^2 + 2 sqrt(3) dH/dt + sqrt(3) k/a^2 - sqrt(3)(V_1 + V_2)/2.

The Hamiltonian equation is its first integral, so it holds along the solution once it holds at
the start, and both constraints are kept here as checks, to 1e-10 over the whole run.

The potentials are the paper's (4.7) and (4.8),

    V_1 = 12 H_0^2 phi_1^2/(1 + phi_1) (1 - ln(1 + phi_1)/(1 + phi_1)),
    V_2 = 2 sqrt(3) H_0^2 phi_2^2/(1 + phi_2^2),

whose slopes keep H^2 + k/a^2 below H_0^2, V_1' rising from 0 toward 12 H_0^2 without reaching
it, and |dH/dt - k/a^2| at most 3 sqrt(3) H_0^2/8, the largest |V_2'|/sqrt(12).

The universe drawn is the time-symmetric one through a bounce at t = 0 with a = A_BOUNCE: there
H = 0, phi_1 solves V_1'(phi_1) = 12/a^2, phi_2 is the root of the Hamiltonian equation at H = 0,

    -(V_1 + V_2)/2 + 3 (2 phi_1 + 1)/a^2 + sqrt(3) phi_2/a^2 = 0,

and dphi_2/dt = 0. At H = 0 the Hamiltonian equation fixes phi_2 and leaves its rate free, since
the coefficient of the rate in it vanishes there, so the bounces form a family of two numbers, the
radius and that rate, and the rate zero is the one symmetric under t -> -t, which takes a, phi_1
and phi_2 into themselves and H and the rates into their negatives. The solution is integrated
from t = 0 to T_FAR and read at negative times through that symmetry. Past the bounce phi_2 runs
off toward -infinity and phi_1 grows, so dH/dt - k/a^2 = -V_2'(phi_2)/sqrt(12) falls toward zero
and the universe approaches de Sitter space, slowly: at a_b = 2 it is 3e-2 at t = 2, 3e-3 at t = 10
and 5e-4 at T_FAR = 60, while H climbs from 0.82 to 0.92 and the scale factor grows as e^(Ht). The
conformal time eta, d eta/dt = 1/a, rides along, zero at the bounce, and has stopped changing in its
tenth digit by t = 30. Its limit is eta_end = 0.8430 at a_b = 2, so the whole history spans 1.686
of conformal time, less than pi, and no light ray ever reaches the antipode of the event it left.

    .venv.noindex/bin/python _tools/derivations/limiting_curvature.py
"""
import math

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

A_BOUNCE = 2.0      # the radius of the 3-sphere at the bounce, in c/H_0
T_FAR = 60.0        # how far the integration runs, in 1/H_0, where a is about 1e22
TOL = dict(rtol=1e-12, atol=1e-14, method="DOP853", dense_output=True)
K = 1.0             # the closed universe
S3 = math.sqrt(3.0)


def V1(p):
    p = np.asarray(p, dtype=float)
    return 12 * p ** 2 / (1 + p) * (1 - np.log1p(p) / (1 + p))


def V1p(p):
    """dV_1/dphi_1 = 12 phi (phi^2 + 2 phi + 2 - 2 ln(1 + phi))/(1 + phi)^3, which rises from 0 to 12."""
    p = np.asarray(p, dtype=float)
    return 12 * p * (p ** 2 + 2 * p + 2 - 2 * np.log1p(p)) / (1 + p) ** 3


def V1pp(p):
    """d^2V_1/dphi_1^2 = 12 (phi^2 - 2 phi + 2 + (4 phi - 2) ln(1 + phi))/(1 + phi)^4, positive for phi >= 0."""
    p = np.asarray(p, dtype=float)
    return 12 * (p ** 2 - 2 * p + 2 + (4 * p - 2) * np.log1p(p)) / (1 + p) ** 4


def V2(p):
    p = np.asarray(p, dtype=float)
    return 2 * S3 * p ** 2 / (1 + p ** 2)


def V2p(p):
    """dV_2/dphi_2 = 4 sqrt(3) phi/(1 + phi^2)^2, largest in size, 9/4, at phi = +-1/sqrt(3)."""
    p = np.asarray(p, dtype=float)
    return 4 * S3 * p / (1 + p ** 2) ** 2


def V2pp(p):
    p = np.asarray(p, dtype=float)
    return 4 * S3 * (1 - 3 * p ** 2) / (1 + p ** 2) ** 3


def rate_of_H(a, f2):
    """(2.26): dH/dt = k/a^2 - V_2'(phi_2)/sqrt(12)."""
    return K / a ** 2 - V2p(f2) / math.sqrt(12)


def rates(t, y):
    """The derivatives of (a, H, phi_1, phi_2, dphi_2/dt, eta)."""
    a, H, f1, f2, q, _ = y
    Hd = rate_of_H(a, f2)
    d1 = -4 * S3 * H * V2p(f2) / V1pp(f1)
    p2dd = (-4 * S3 * d1 * H - 6 * H * q
            + f1 * (-6 * S3 * H ** 2 - 4 * S3 * Hd + 2 * S3 * K / a ** 2)
            + f2 * (-9 * H ** 2 - 6 * Hd + K / a ** 2)
            + 3 * S3 * H ** 2 + 2 * S3 * Hd + S3 * K / a ** 2 - S3 * (V1(f1) + V2(f2)) / 2)
    return [a * H, Hd, d1, q, p2dd, 1 / a]


def constraints(y):
    """(2.25) and the Hamiltonian equation, each as a residual that vanishes on a solution."""
    a, H, f1, f2, q = y[:5]
    first = H ** 2 + K / a ** 2 - V1p(f1) / 12
    hamiltonian = (-(V1(f1) + V2(f2)) / 2 + 3 * H ** 2 * (1 - 2 * f1) + 3 * (K / a ** 2) * (2 * f1 + 1)
                   - S3 * H * (q + 3 * H * f2) + S3 * K * f2 / a ** 2)
    return first, hamiltonian


def bounce_data(a_b):
    """(a, H, phi_1, phi_2, dphi_2/dt, eta) at the symmetric bounce of radius a_b > 1, with phi_2
    the root of the Hamiltonian equation at H = 0 nearest zero."""
    if a_b <= 1:
        raise ValueError("a bounce needs a > c/H_0, since H^2 + k/a^2 = V_1'/12 < H_0^2")
    f1 = brentq(lambda p: V1p(p) - 12 * K / a_b ** 2, 1e-15, 1e12, xtol=1e-15, rtol=1e-15)

    def at_rest(p):
        return -(V1(f1) + V2(p)) / 2 + 3 * (K / a_b ** 2) * (2 * f1 + 1) + S3 * K * p / a_b ** 2
    grid = np.linspace(-40, 40, 80001)
    values = at_rest(grid)
    crossings = grid[:-1][np.sign(values[:-1]) != np.sign(values[1:])]
    if not len(crossings):
        raise ValueError(f"no phi_2 satisfies the Hamiltonian equation at a bounce of radius {a_b}")
    near = crossings[np.argmin(np.abs(crossings))]
    f2 = brentq(at_rest, near - 1e-3, near + 2e-3, xtol=1e-15, rtol=1e-15)
    return np.array([a_b, 0.0, f1, f2, 0.0, 0.0])


class Bounce:
    """The time-symmetric closed universe through a bounce of radius a_b, as numbers."""

    def __init__(self, a_b=A_BOUNCE):
        self.a_b = a_b
        self.start = bounce_data(a_b)
        self.run = solve_ivp(rates, (0.0, T_FAR), self.start, **TOL)
        if self.run.status != 0:
            raise RuntimeError(f"the integration stopped: {self.run.message}")
        self.eta_end = float(self.run.y[5, -1])

    def state(self, t):
        """(a, H, phi_1, phi_2, dphi_2/dt, eta) at the times t, rows by quantity, through the
        symmetry t -> -t, which keeps a, phi_1 and phi_2 and reverses H, the rate of phi_2 and eta."""
        t = np.asarray(t, dtype=float)
        flat = np.abs(t.ravel())
        if np.any(flat > T_FAR):
            raise ValueError(f"a time beyond T_FAR = {T_FAR}")
        y = self.run.sol(flat)
        sign = np.sign(t.ravel())
        for row in (1, 4, 5):
            y[row] = y[row] * np.where(sign == 0, 1, sign)
        return y.reshape((6,) + t.shape)

    def jet(self, t):
        """a and its first three derivatives along t, as a dict of arrays: da = aH, d2a = a(dH/dt +
        H^2), d3a = da (dH/dt + H^2) + a (d^2H/dt^2 + 2H dH/dt), with d^2H/dt^2 = -2kH/a^2 -
        V_2''(phi_2) (dphi_2/dt)/sqrt(12) by (2.26)."""
        a, H, f1, f2, q, eta = self.state(t)
        Hd = rate_of_H(a, f2)
        Hdd = -2 * K * H / a ** 2 - V2pp(f2) * q / math.sqrt(12)
        da = a * H
        return {"a": a, "da": da, "d2a": a * (Hd + H ** 2), "d3a": da * (Hd + H ** 2) + a * (Hdd + 2 * H * Hd),
                "H": H, "dH": Hd, "phi1": f1, "phi2": f2, "eta": eta}

    def a(self, t):
        return self.state(t)[0]

    def eta(self, t):
        return self.state(t)[5]

    def time_at(self, eta):
        """The time t at the conformal time eta, |eta| < ETA_END, by bisection on the stored eta(t)."""
        eta = np.asarray(eta, dtype=float)
        out = np.empty(eta.size)
        for i, e in enumerate(eta.ravel()):
            if abs(e) >= self.eta_end:
                raise ValueError(f"eta = {e} lies beyond the end of the universe, {self.eta_end}")
            out[i] = math.copysign(brentq(lambda t: self.run.sol(t)[5] - abs(e), 0.0, T_FAR, xtol=1e-14), e)
        return out.reshape(eta.shape)

    def residual(self, t):
        """The largest size of the two constraints at the times t."""
        y = self.state(t)
        return float(np.max(np.abs(np.vstack(constraints(y)))))

    def facts(self):
        times = np.linspace(0, T_FAR, 4001)
        y = self.state(times)
        far = self.state(np.array([T_FAR]))
        return {
            "a_b": self.a_b, "phi1(0)": float(self.start[2]), "phi2(0)": float(self.start[3]),
            "residual": self.residual(times), "eta_end": self.eta_end,
            "H(T_FAR)": float(far[1][0]), "dH - k/a^2 at T_FAR": float(-V2p(far[3][0]) / math.sqrt(12)),
            "largest H^2 + k/a^2": float(np.max(y[1] ** 2 + K / y[0] ** 2)),
        }


_BOUNCES = {}


def bounce(a_b=A_BOUNCE):
    """The universe through a bounce of radius a_b, solved once."""
    if a_b not in _BOUNCES:
        _BOUNCES[a_b] = Bounce(a_b)
    return _BOUNCES[a_b]


if __name__ == "__main__":
    for key, value in bounce().facts().items():
        print(f"{key:28s} {value:.12g}")
