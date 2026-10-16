#!/usr/bin/env python3
"""Kodama's kink as numbers: a static wormhole held open by a ghost scalar field S with the
potential V = (mu/f)^2 (1 - f^2 S^2)^2/2, in units where c = 1 and the unit of length is 1/mu.

No closed form is known for the metric, so its two functions are integrated. Along the proper
radial distance R of T. Kodama, L. C. S. de Oliveira and F. C. Santos (CBPF-NF-036/78, their (2)),
ds^2 = -e^{2h} dt^2 + dR^2 + r^2 dOmega^2, the ghost field's stresses are

    G^t_t         = S'^2 + V,
    G^R_R         = V - S'^2,
    G^theta_theta = S'^2 + V,

Kodama's (3) to (5) of CBPF-NF-003/78 with epsilon = -1, and the field obeys Box S = S dV/d(S^2),

    S'' = -(h' + 2 r'/r) S' - mu^2 S (1 - f^2 S^2).

The published G^t_t and G^theta_theta of the spacetime's own proper radial chart, read from its
file with h and r as plain symbols, give r'' and h''; the published G^R_R, which the Bianchi
identity makes a consequence of the other three, is left for a check, Kink.radial, and at the
throat, where r' = 0, it fixes the slope of the field, S'(0)^2 = (1/2 + f^2/r_0^2)/f^2.

The solution is even in R for r and h and odd for S, so it is integrated from the throat outward,
from r(0) = r_0, r'(0) = 0, h'(0) = 0, S(0) = 0. The field reaches its vacuum S = 1/f far away only
for one radius of the throat, which Kink finds by bisection on r_0: too small and the field turns
back before it gets there, too large and it overshoots. The field of the last radius that turns
back is followed to its turning point, R = edge, where 1 - f S has fallen to about 10^-9, and from
there on the metric is Schwarzschild's of the mass inside, m = r (1 - r'^2)/2, written along the
proper distance, and h is shifted to meet it there, so that h vanishes far away.

The areal chart of Kodama's CBPF-NF-003/78, ds^2 = -e^{2 eta} dt^2 + e^{2 alpha} dr^2 + r^2 dOmega^2,
is the same kink on either side of the throat, with eta = h and e^{-alpha} = dr/dR; Areal hands
it over as functions of r.

    .venv.noindex/bin/python _tools/derivations/kodama_kink.py

prints the declared kink and a few others. null_rays.py, conformal.py and embedding.py draw with
these classes, and _tools/test_kodama_kink.py holds them to the physics.
"""
import json
import sys
from functools import lru_cache
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
import build_mfs_data as build  # noqa: E402
import verify_metrics as vm  # noqa: E402

METRIC, PROPER, AREAL = "kodama_kink", "proper_radial", "areal"
F = 1.25                 # the declared kink: the f of Kodama's figures 1 to 3 in CBPF-NF-036/78
# What the declared kink is stated to be, each to the digits the drawings print.
STATED = {"r0": (1.475, 3), "m": (0.390, 3), "lapse": (0.703, 3)}

# The declared kink in words, which every drawing of it prints.
INPUT = ("Kodama's kink at $f = 1.25$, in units where $c = \\mu = 1$, solved from this spacetime's own "
         "$G^t{}_t$ and $G^\\theta{}_\\theta$ and the field's wave equation: a throat of areal radius "
         "$r_0 = 1.475/\\mu$, a mass $M$ on either side with $GM/c^2 = 0.390/\\mu$, and clocks at the throat running at "
         "$e^{h} = 0.703$ of the rate far away.")
SETTINGS = "$\\mu = 1$ and $f = 1.25$, so that the unit of every length is $1/\\mu$."


def _load(metric_id, system_id):
    metric = json.loads((build.METRICS_DIR / f"{metric_id}.json").read_text(encoding="utf-8"))
    entry = next(e for e in metric["coordinates"] if e["id"] == system_id)
    reader = vm.Reader(entry["coords"], [p["symbol"] for p in entry.get("parameters", [])],
                       vm.time_coordinates(vm.DIMENSIONS[(metric_id, system_id)], entry["coords"]))
    return entry, reader


@lru_cache(maxsize=None)
def equations(metric_id=METRIC, system_id=PROPER):
    """r'' and h'' as functions of (R, r, r', h', S, S', f) from the published G^t_t and
    G^theta_theta, with mu = 1, and the published G^R_R less the field's stress along the radius
    as a function of (R, r, r', h', S, S', f)."""
    entry, reader = _load(metric_id, system_id)
    R = reader.symbol["R"]
    fh, fr = reader.parameters["h"], reader.parameters["r"]
    held = dict(zip(entry["coords"], (0, R, sp.pi / 2, 0)))
    angles = {reader.symbol[c]: v for c, v in held.items() if c not in ("t", "R")}
    X, X1, X2, H1, H2, S, S1, f = sp.symbols("r r1 r2 h1 h2 S S1 f")
    ul = entry["einstein_tensor"]["variants"]["ul"]["nonzero"]

    def published(index):
        e = reader(next(c["value"] for c in ul if c["indices"] == [index, index]))
        e = e.subs({sp.Derivative(fh, (R, 2)): H2, sp.Derivative(fr, (R, 2)): X2})
        e = e.subs({sp.Derivative(fh, R): H1, sp.Derivative(fr, R): X1})
        return e.subs({fr: X, reader.c: 1}).subs(angles)
    V = (1 - f ** 2 * S ** 2) ** 2 / (2 * f ** 2)
    d2r = sp.solve(published("t") - (S1 ** 2 + V), X2)[0]
    d2h = sp.solve(published("\\theta").subs(X2, d2r) - (S1 ** 2 + V), H2)[0]
    d2S = -(H1 + 2 * X1 / X) * S1 - S * (1 - f ** 2 * S ** 2)
    state = (R, X, X1, H1, S, S1, f)
    return {
        "second": sp.lambdify(state, [sp.simplify(d2r), sp.simplify(d2h), d2S], "numpy"),
        "radial": sp.lambdify(state, published("R") - (V - S1 ** 2), "numpy"),
    }


class Kink:
    """The kink with the parameter f. `funcs` are the declared functions of the proper radial
    chart, which `values` and `at` hand back with their first two derivatives along R, for R of
    either sign."""

    funcs = ("h", "r")
    TABLE = 40001

    def __init__(self, f=F, metric_id=METRIC, system_id=PROPER):
        self.f = float(f)
        self.eq = equations(metric_id, system_id)
        low, high = 1e-3, 4.0 * self.f ** 2
        fate_low, fate_high = self._fate(low), self._fate(high)
        if fate_low == fate_high:
            raise ArithmeticError(f"Kodama's kink: no throat between {low} and {high} at f = {f}")
        while high - low > 4e-16 * high:
            mid = 0.5 * (low + high)
            if self._fate(mid) == fate_low:
                low = mid
            else:
                high = mid
        # The last radius whose field turns back is followed to where it turns.
        turning = low if fate_low == "turns" else high
        self.r0 = 0.5 * (low + high)
        self.ivp = self._shoot(turning, dense=True)
        self.edge = float(self.ivp.t_events[1][0])
        r, r1, h, h1, S, S1 = self.ivp.sol(self.edge)
        self.left = float(1 - self.f * S)                    # what is left of the field at the edge
        self.m = float(0.5 * r * (1 - r1 ** 2))
        self.r_edge = float(r)
        self.shift = float(0.5 * np.log(1 - 2 * self.m / r) - h)
        self.lapse = float(np.exp(self.ivp.sol(0.0)[2] + self.shift))   # e^h at the throat

    def _start(self, r0):
        slope = np.sqrt((0.5 + self.f ** 2 / r0 ** 2)) / self.f
        return [r0, 0.0, 0.0, 0.0, 0.0, slope]

    def _rhs(self, R, y):
        r, r1, h, h1, S, S1 = y
        d2r, d2h, d2S = self.eq["second"](R, r, r1, h1, S, S1, self.f)
        return [r1, d2r, h1, d2h, S1, d2S]

    def _shoot(self, r0, dense=False):
        def overshoots(R, y):
            return self.f * y[4] - 1.0

        def turns(R, y):
            return y[5]
        overshoots.terminal = turns.terminal = True
        turns.direction = -1
        return solve_ivp(self._rhs, (0.0, 200.0), self._start(r0), events=(overshoots, turns), rtol=1e-12,
                         atol=1e-14, method="DOP853", dense_output=dense)

    def _fate(self, r0):
        sol = self._shoot(r0)
        if sol.t_events[0].size:
            return "overshoots"
        if sol.t_events[1].size:
            return "turns"
        return "neither"

    # ---------------------------------------------------------------- the solution

    def inside(self, R):
        """(r, r', r'', h, h', h'', S, S') at proper distances 0 <= R <= edge."""
        R = np.atleast_1d(np.asarray(R, dtype=float))
        r, r1, h, h1, S, S1 = self.ivp.sol(R)
        d2r, d2h, _ = self.eq["second"](R, r, r1, h1, S, S1, self.f)
        return r, r1, d2r, h + self.shift, h1, d2h, S, S1

    def _outside_r(self, R):
        """The areal radius at proper distances R >= edge, from Schwarzschild's
        R - edge = F(r) - F(r_edge), F(r) = sqrt(r (r - 2m)) + 2m ln(sqrt(r) + sqrt(r - 2m))."""
        m = self.m

        def F(r):
            return np.sqrt(r * (r - 2 * m)) + 2 * m * np.log(np.sqrt(r) + np.sqrt(r - 2 * m))
        target = F(self.r_edge) + (np.asarray(R, dtype=float) - self.edge)
        r = np.maximum(self.r_edge, np.asarray(R, dtype=float) - self.edge + self.r_edge)
        for _ in range(60):
            step = (F(r) - target) * np.sqrt(1 - 2 * m / r)
            r = np.maximum(r - step, 2 * m * (1 + 1e-12))
            if np.all(np.abs(step) <= 1e-13 * r):
                break
        return r

    def values(self, R):
        """(value, first, second derivative along R) of h and of r at proper distances R of
        either sign: r and h are even about the throat."""
        R = np.atleast_1d(np.asarray(R, dtype=float))
        a = np.abs(R)
        sign = np.where(R < 0, -1.0, 1.0)
        out = {name: [np.empty_like(R) for _ in range(3)] for name in self.funcs}
        near = a < self.edge
        if near.any():
            r, r1, r2, h, h1, h2, _, _ = self.inside(a[near])
            for store, v in zip(out["r"], (r, sign[near] * r1, r2)):
                store[near] = v
            for store, v in zip(out["h"], (h, sign[near] * h1, h2)):
                store[near] = v
        if (~near).any():
            r, m = self._outside_r(a[~near]), self.m
            g = 1 - 2 * m / r
            r1 = np.sqrt(g)
            r2 = m / r ** 2
            h1 = m / (r ** 2 * np.sqrt(g))
            h2 = -m * (2 * r - 3 * m) / (r ** 4 * g)
            for store, v in zip(out["r"], (r, sign[~near] * r1, r2)):
                store[~near] = v
            for store, v in zip(out["h"], (0.5 * np.log(g), sign[~near] * h1, h2)):
                store[~near] = v
        return out

    def field(self, R):
        """S and S' at proper distances R of either sign: S is odd about the throat, and 1/f from
        the edge on, where what is left is dropped."""
        R = np.atleast_1d(np.asarray(R, dtype=float))
        a = np.minimum(np.abs(R), self.edge)
        sign = np.where(R < 0, -1.0, 1.0)
        _, _, _, _, _, _, S, S1 = self.inside(a)
        inside = np.abs(R) < self.edge
        return (sign * np.where(inside, S, 1 / self.f), np.where(inside, S1, 0.0))

    def _table(self):
        if getattr(self, "_spline", None) is None:
            grid = np.linspace(0.0, self.edge, self.TABLE)
            v = self.values(np.minimum(grid, self.edge * (1 - 1e-15)))
            self._spline = CubicSpline(grid, np.column_stack([v[name][k] for name in self.funcs for k in range(3)]))
        return self._spline

    def at(self, name, x0, x):
        """A declared function and its first two derivatives at chart points (x^0, R), from a cubic
        spline inside the edge and closed forms outside it, as the drawings ask point by point."""
        shape = np.broadcast(np.asarray(x0), np.asarray(x)).shape
        flat = np.broadcast_to(np.asarray(x, dtype=float), shape).ravel()
        a, sign = np.abs(flat), np.where(flat < 0, -1.0, 1.0)
        column = 3 * self.funcs.index(name)
        inside = a < self.edge
        out = np.empty((3, flat.size))
        if inside.any():
            out[:, inside] = self._table()(a[inside])[:, column:column + 3].T
        if (~inside).any():
            out[:, ~inside] = np.array(self.values(a[~inside])[name])
        out[1] *= sign
        return [v.reshape(shape) for v in out]

    def radial(self, R):
        """The published G^R_R less the field's stress along the radius, along the solution,
        against the energy density at the throat: zero if the kink solves the one field equation
        its construction did not use."""
        r, r1, _, _, h1, _, S, S1 = self.inside(np.asarray(R, dtype=float))
        miss = self.eq["radial"](R, r, r1, h1, S, S1, self.f)
        return miss / (self.inside(0.0)[7][0] ** 2 + 1 / (2 * self.f ** 2))

    def stated(self):
        """What of STATED the kink is not, to the digits stated."""
        mine = {"r0": self.r0, "m": self.m, "lapse": self.lapse}
        return {k: mine[k] for k, (value, digits) in STATED.items() if round(mine[k], digits) != value}

    def areal(self):
        return Areal(self)

    def embedding_height(self, R):
        """The height z of the slice t = 0, theta = pi/2 embedded in flat space, as a function of R:
        dz/dR = sqrt(1 - r'^2), zero at the throat and odd about it."""
        R = np.atleast_1d(np.asarray(R, dtype=float))
        grid = np.linspace(0.0, np.max(np.abs(R)), 20001)
        slope = np.sqrt(np.clip(1 - self.values(grid)["r"][1] ** 2, 0.0, None))
        z = np.concatenate([[0.0], np.cumsum(0.5 * (slope[1:] + slope[:-1]) * np.diff(grid))])
        return np.sign(R) * np.interp(np.abs(R), grid, z)


class Areal:
    """The kink in Kodama's areal radius r >= r_0, on either side of the throat: eta and alpha,
    each with its first two derivatives along r. eta = h and e^{-alpha} = dr/dR at the R >= 0
    where r(R) = r, which is found from a table in rho = sqrt(r - r_0), since r - r_0 grows as R^2
    from the throat. alpha diverges at the throat, where dr/dR = 0."""

    funcs = ("eta", "alpha")

    def __init__(self, kink):
        self.kink = kink
        self.r0 = kink.r0
        grid = np.concatenate([np.linspace(0.0, kink.edge, 20001), kink.edge * np.geomspace(1.0001, 50.0, 4000)])
        r = kink.values(grid)["r"][0]
        self.rho = np.sqrt(np.maximum(r - self.r0, 0.0))
        self.R = grid

    def R_of(self, r):
        """The proper distance R >= 0 to the sphere of areal radius r, by Newton's method from the table."""
        r = np.atleast_1d(np.asarray(r, dtype=float))
        R = np.interp(np.sqrt(np.maximum(r - self.r0, 0.0)), self.rho, self.R)
        for _ in range(8):
            v = self.kink.values(R)["r"]
            step = np.where(v[1] > 1e-12, (v[0] - r) / np.where(v[1] > 1e-12, v[1], 1.0), 0.0)
            R = np.maximum(R - step, 0.0)
        return R

    def values(self, r):
        """(value, first, second derivative along r) of eta and of alpha at areal radii r > r_0."""
        r = np.atleast_1d(np.asarray(r, dtype=float))
        R = self.R_of(r)
        v = self.kink.values(R)
        (_, r1, r2), (h, h1, h2) = v["r"], v["h"]
        with np.errstate(divide="ignore", invalid="ignore"):
            return {"eta": [h, h1 / r1, (h2 - h1 * r2 / r1) / r1 ** 2],
                    "alpha": [-np.log(r1), -r2 / r1 ** 2, (2 * r2 ** 2 / r1 - self._r3(R)) / r1 ** 3]}

    def _r3(self, R):
        """r''' along R, by differencing r'' on a step small against the scale of the kink."""
        d = 1e-5
        return (self.kink.values(R + d)["r"][2] - self.kink.values(np.maximum(R - d, 0.0))["r"][2]) / (
            (R + d) - np.maximum(R - d, 0.0))

    def at(self, name, x0, x):
        shape = np.broadcast(np.asarray(x0), np.asarray(x)).shape
        flat = np.broadcast_to(np.asarray(x, dtype=float), shape).ravel()
        out = np.array(self.values(flat)[name])
        return [v.reshape(shape) for v in out]


@lru_cache(maxsize=None)
def kink(f=F):
    return Kink(f)


def threshold(low=0.6, high=0.7):
    """The least f with a kink, where the throat closes, by bisection on whether one exists."""
    def exists(f):
        try:
            return kink(round(f, 10)).r0 > 0
        except ArithmeticError:
            return False
    for _ in range(30):
        mid = 0.5 * (low + high)
        if exists(mid):
            high = mid
        else:
            low = mid
    return high


if __name__ == "__main__":
    k = kink()
    print(f"declared, f = {k.f}: r_0 = {k.r0:.6f}, m = {k.m:.6f}, e^h at the throat = {k.lapse:.6f}, "
          f"edge = {k.edge:.2f} with 1 - fS = {k.left:.1e}")
    R = np.linspace(0.01, 0.98 * k.edge, 400)
    print(f"G^R_R misses by {np.max(np.abs(k.radial(R))):.1e}")
    for f in (0.77, 1.0, 2.0, 5.0, 20.0):
        o = kink(f)
        print(f"f = {f}: r_0 = {o.r0:.5f}, r_0/f^2 = {o.r0 / f ** 2:.5f}, m = {o.m:.5f}, r_0/m = {o.r0 / o.m:.4f}")
