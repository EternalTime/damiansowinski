#!/usr/bin/env python3
"""The relativistic Larson-Penston solution of self similar fluid collapse, at k = 1/100.

Harada and Maeda, arXiv:gr-qc/0101064, section II: in x = ln z, with z = r/(-t),

    (ln M)' = k (1/y - 1)/(1 + k),   (ln S)' = -(1 - y)/(1 + k),
    (ln eta)' = 2 - (4 y V^2 - (1 + k)^2 S^-4 eta^-(1-k)/(1+k)) / (2 (V^2 - k)),

y = M/(eta S^3), V^2 = z^(2(1-k)/(1+k)) eta^(-2(1-k)/(1+k)) S^-4 / a and a = (2D)^(2k/(1+k)), the
gauge in which e^sigma = 1 at the centre, with the constraint (constraint) holding throughout.
A regular centre fixes the solution near z = 0 up to D; the sonic point, V^2 = k, is crossed
analytically along the direction Harada calls secondary. The solution is the one that matches:
the centre's series integrated outward, and the analytic expansion at a sonic point on the sonic
line integrated inward, agree in ln S and ln eta at z = 1/1000, which fixes D and y at the sonic
point by Newton's method. Outside the sonic point the analytic expansion is integrated outward.
The metric functions are then e^sigma = (2D z^2/eta)^(2k/(1+k)) and e^omega = eta^(-2/(1+k)) S^-4.
"""
import math
from functools import lru_cache

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from scipy.optimize import fsolve

K = 0.01
CENTRE = math.log(1e-12)
MATCH = math.log(1e-3)
FAR = math.log(1e8)
STEP = 1e-7
TOLERANCES = {"rtol": 1e-12, "atol": 1e-14, "method": "DOP853"}
GUESS = (1.4377, 0.4146)


def parts(x, Y, D, k=K):
    M, S, eta = math.exp(Y[0]), math.exp(Y[1]), math.exp(Y[2])
    z = math.exp(x)
    a = (2 * D) ** (2 * k / (1 + k))
    y = M / (eta * S ** 3)
    V2 = z ** (2 * (1 - k) / (1 + k)) * eta ** (-2 * (1 - k) / (1 + k)) * S ** -4 / a
    num = 4 * y * V2 - (1 + k) ** 2 * S ** -4 * eta ** (-(1 - k) / (1 + k))
    return M, S, eta, z, y, V2, num


def rates(x, Y, D, k=K):
    M, S, eta, z, y, V2, num = parts(x, Y, D, k)
    return [k / (1 + k) * (1 / y - 1), -(1 - y) / (1 + k), 2 - num / (2 * (V2 - k))]


def constraint(x, Y, D, k=K):
    M, S, eta, z, y, V2, num = parts(x, Y, D, k)
    return (1 - y) ** 2 * V2 - (k + y) ** 2 + (1 + k) ** 2 * eta ** (-2 / (1 + k)) * S ** -6 * (1 - M / S)


def centre(D, z, k=K):
    """Harada and Maeda's leading terms at a regular centre, ln M, ln S and ln eta."""
    ck = 3 * (1 + k) / (1 + 3 * k)
    return [math.log(ck / 3 * (2 * D) ** (k / (1 + k)) * z ** (2 * k / (1 + k))),
            math.log(ck ** (1 / 3) * (2 * D) ** (-1 / (3 * (1 + k))) * z ** (-2 / (3 * (1 + k)))),
            math.log(2 * D * z * z)]


def sonic_state(ys, D, k=K):
    """The point of the sonic line with y = ys: V^2 = k, the regularity condition and the
    constraint solved for z, M, S and eta."""
    p = 2 * (1 - k) / (1 + k)
    a = (2 * D) ** (2 * k / (1 + k))
    w = 4 * ys * k / (1 + k) ** 2
    q = (1 + k) ** 2 * w ** 1.5 / ((k + ys) ** 2 + (1 + k) ** 2 * ys * w - (1 - ys) ** 2 * k)
    eta = q ** (2 * (1 + k) / (1 + 3 * k))
    S = (w * eta ** (p / 2)) ** -0.25
    M = ys * eta * S ** 3
    z = math.sqrt(eta * (k * a / w) ** (2 / p))
    return math.log(z), np.array([math.log(M), math.log(S), math.log(eta)])


def sonic_slopes(x, Y, D, k=K, h=1e-7):
    """The two analytic values of (ln eta)' at a sonic point, smallest first, from l'Hopital's
    rule: num and V^2 - k both vanish, and their rates along the flow are linear in (ln eta)'."""
    M, S, eta, z, y, V2, num = parts(x, Y, D, k)
    dlM, dlS = k / (1 + k) * (1 / y - 1), -(1 - y) / (1 + k)

    def along(lam):
        *_, V2h, numh = parts(x + h, Y + h * np.array([dlM, dlS, lam]), D, k)
        return (numh - num) / h, (V2h - V2) / h
    n0, v0 = along(0.0)
    n1, v1 = along(1.0)
    A, B, C, E = n1 - n0, n0, v1 - v0, v0
    roots = np.roots([2 * C, 2 * E - 4 * C + A, -4 * E + B])
    return sorted(float(np.real(root)) for root in roots), (dlM, dlS)


def from_sonic(D, ys, end, k=K):
    x, Y = sonic_state(ys, D, k)
    (lam, _), (dlM, dlS) = sonic_slopes(x, Y, D, k)
    sign = 1 if end > x else -1
    start = Y + sign * STEP * np.array([dlM, dlS, lam])
    return solve_ivp(rates, (x + sign * STEP, end), start, args=(D, k), dense_output=True, **TOLERANCES)


def from_centre(D, end, k=K):
    return solve_ivp(rates, (CENTRE, end), centre(D, math.exp(CENTRE), k), args=(D, k),
                     dense_output=True, **TOLERANCES)


def mismatch(values, k=K):
    D, ys = values
    inner, outer = from_centre(D, MATCH, k), from_sonic(D, ys, MATCH, k)
    if inner.status != 0 or outer.status != 0:
        return [10.0, 10.0]
    return [outer.y[1, -1] - inner.y[1, -1], outer.y[2, -1] - inner.y[2, -1]]


@lru_cache(maxsize=2)
def solve(k=K):
    """D, y at the sonic point, ln z at the sonic point, and the three pieces of the solution."""
    D, ys = fsolve(mismatch, GUESS, args=(k,), xtol=1e-13)
    if max(abs(v) for v in mismatch((D, ys), k)) > 1e-9:
        raise RuntimeError("the centre and the sonic point do not match")
    xs, _ = sonic_state(ys, D, k)
    pieces = (from_centre(D, MATCH, k), from_sonic(D, ys, MATCH, k), from_sonic(D, ys, FAR, k))
    if any(piece.status != 0 for piece in pieces):
        raise RuntimeError("an integration of the solution did not finish")
    return float(D), float(ys), float(xs), pieces


def logs(x, k=K):
    """ln M, ln S and ln eta at ln z = x, from the piece that covers it."""
    D, _, xs, (inner, middle, outer) = solve(k)
    if x <= MATCH:
        return inner.sol(x)
    return middle.sol(x) if x < xs else outer.sol(x)


@lru_cache(maxsize=2)
def tables(k=K, n=9001):
    """ln S and ln eta and their exact rates along x = ln z on a grid, each a cubic spline, with
    the rates' decay past the last point: there y tends to 1 and both rates fall as a power of z."""
    D, _, xs, _ = solve(k)
    grid = np.linspace(CENTRE, FAR, n)
    states = np.array([logs(x, k) for x in grid])
    slopes = np.array([rates(x, Y, D, k) for x, Y in zip(grid, states)])
    splines = {name: (CubicSpline(grid, states[:, i]), CubicSpline(grid, slopes[:, i]))
               for name, i in (("lnS", 1), ("lneta", 2))}
    decay = {name: math.log(slopes[-2, i] / slopes[-1, i]) / (grid[-1] - grid[-2])
             for name, i in (("lnS", 1), ("lneta", 2))}
    return grid, D, splines, decay


def logs_and_rates(x, k=K):
    """ln S and ln eta with their first two derivatives along x = ln z. Below the grid the
    centre's limiting rates hold; beyond it each rate decays as e^(-q x), q from the last points."""
    grid, D, splines, decay = tables(k)
    out = {}
    for name, (value, rate) in splines.items():
        inside = np.clip(x, grid[0], grid[-1])
        v, v1, v2 = value(inside), rate(inside), rate(inside, 1)
        low, high = x < grid[0], x > grid[-1]
        q = decay[name]
        r_end = rate(grid[-1])
        beyond = r_end * np.exp(-q * (x - grid[-1]))
        v = np.where(low, value(grid[0]) + rate(grid[0]) * (x - grid[0]),
                     np.where(high, value(grid[-1]) + (r_end - beyond) / q, v))
        v1 = np.where(low, rate(grid[0]), np.where(high, beyond, v1))
        v2 = np.where(low, 0.0, np.where(high, -q * beyond, v2))
        out[name] = (v, v1, v2)
    return D, out


def jet(z, k=K):
    """sigma, omega and S with their first two derivatives along z itself, at the z given, from
    e^sigma = (2D z^2/eta)^(2k/(1+k)) and e^omega = eta^(-2/(1+k)) S^-4."""
    z = np.asarray(z, dtype=float)
    shape = z.shape
    flat = np.atleast_1d(z).ravel()
    x = np.log(np.where(flat > 0, flat, np.nan))
    D, out = logs_and_rates(x, k)
    lS, lS1, lS2 = out["lnS"]
    le, le1, le2 = out["lneta"]
    c = 2 * k / (1 + k)
    along_x = {"sigma": (c * (math.log(2 * D) + 2 * x - le), c * (2 - le1), -c * le2),
               "omega": (-2 / (1 + k) * le - 4 * lS, -2 / (1 + k) * le1 - 4 * lS1, -2 / (1 + k) * le2 - 4 * lS2)}
    result = {}
    for name, (v, v1, v2) in along_x.items():
        result[name], result["d" + name], result["d2" + name] = v, v1 / flat, (v2 - v1) / flat ** 2
    S = np.exp(lS)
    result["S"], result["dS"], result["d2S"] = S, S * lS1 / flat, S * (lS2 + lS1 * lS1 - lS1) / flat ** 2
    return {key: value.reshape(shape) for key, value in result.items()}


def marked_lines(k=K):
    """ln z of the sonic point, V_z^2 = k, and of the null orbit of the scaling symmetry,
    V_z^2 = 1, where V_z = -z e^((omega - sigma)/2)."""
    D, _, xs, _ = solve(k)
    from scipy.optimize import brentq

    def speed(x):
        j = jet(math.exp(x), k)
        return math.exp(x) ** 2 * math.exp(float(j["omega"]) - float(j["sigma"])) - 1
    return xs, brentq(speed, xs, math.log(10.0))


def declared_functions(k=K):
    """SymPy functions of z carrying sigma, omega and S with their first two derivatives."""
    import sympy as sp

    def make(name, key, next_function=None):
        def value(z):
            return jet(z, k)[key]

        def fdiff(self, argindex=1):
            if argindex != 1 or next_function is None:
                raise sp.ArgumentIndexError(self, argindex)
            return next_function(self.args[0])

        return type(name, (sp.Function,), {
            "nargs": 1, "is_real": True, "_imp_": staticmethod(value), "fdiff": fdiff,
            "_eval_evalf": lambda self, prec: sp.Float(float(self._imp_(float(self.args[0]))), prec)})

    result = {}
    for key in ("sigma", "omega", "S"):
        second = make("op_d2" + key, "d2" + key)
        first = make("op_d" + key, "d" + key, second)
        base = make("op_" + key, key, first)
        result.update({fn.__name__: fn for fn in (base, first, second)})
    return result


if __name__ == "__main__":
    D, ys, xs, pieces = solve()
    zs, z1 = (math.exp(v) for v in marked_lines())
    print(f"D = {D:.10f}, y at the sonic point = {ys:.10f}, sonic z = {zs:.8f}, null orbit z = {z1:.8f}")
