#!/usr/bin/env python3
"""A boson star as numbers: the ground state of a complex scalar field held together by its own
gravity, in units where G = c = 1 and the unit of length is 1/mu = hbar/(m c), m the mass of
the boson, so that a mass is counted in M_Pl^2/m.

No closed form is known for the metric, so its two functions are integrated. The field is
phi = phi_0(r) e^{i omega t} with sigma = sqrt(4 pi G) phi_0/c^2, and its stress tensor,
T_ab = (d_a phi* d_b phi + d_a phi d_b phi*)/2 - g_ab (g^cd d_c phi* d_d phi + mu^2 |phi|^2)/2, has

    8 pi T^t_t         = -(omega^2 sigma^2/alpha^2 + sigma'^2/a^2 + sigma^2),
    8 pi T^r_r         =   omega^2 sigma^2/alpha^2 + sigma'^2/a^2 - sigma^2,
    8 pi T^theta_theta =   omega^2 sigma^2/alpha^2 - sigma'^2/a^2 - sigma^2

in the areal chart ds^2 = -alpha^2 dt^2 + a^2 dr^2 + r^2 dOmega^2 (S. L. Liebling and
C. Palenzuela, Living Rev. Relativ. 26, 1, 2023, arXiv:1202.5809v5, their (33)). The published G^t_t and G^r_r of
the spacetime's own areal chart, read from its file with alpha and a as plain symbols, are
linear in a' and alpha' and give them; the wave equation gives

    sigma'' = -(2/r + alpha'/alpha - a'/a) sigma' - a^2 (omega^2/alpha^2 - 1) sigma;

and the published G^theta_theta, which the Bianchi identity makes a consequence of those three,
is left for a check, Star.theta_theta.

A solution regular at the centre has a(0) = 1, sigma'(0) = 0 and a chosen sigma(0), and decays
far away only for a discrete set of frequencies. Star finds the lowest, the state with no node,
by bisection on omega^2 with alpha(0) = 1: too low and the field turns back up before it reaches
zero, too high and it crosses zero. Within `core` of the centre the three functions are their
power series in r^2, centre_series, and each integration starts from the series there. The field
of the last frequency that turns up is followed to its turning point, r = edge, where sigma has
fallen by eight or nine orders of magnitude, and from there on the metric is Schwarzschild's of
the mass inside, which differs from the star's by the square of what is left of the field, and
alpha is scaled to meet it, which scales omega with it. The star has no surface: the field
falls as exp(-sqrt(1 - omega^2) r) for ever.

Isotropic carries the same star into the isotropic radius rho of
ds^2 = -alpha^2 dt^2 + psi^4 (d rho^2 + rho^2 dOmega^2), with r = psi^2 rho and
d ln(rho)/d ln(r) = a, fixed by Schwarzschild's rho = (r - M + sqrt(r^2 - 2 M r))/2 at the edge.

    .venv.noindex/bin/python _tools/derivations/boson_star.py

prints the declared star and the heaviest one. null_rays.py, conformal.py and embedding.py draw
with these classes, and _tools/test_boson_star.py holds them to the physics.
"""
import json
import sys
from functools import lru_cache
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq, minimize_scalar

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
import build_mfs_data as build  # noqa: E402
import verify_metrics as vm  # noqa: E402

METRIC, AREAL = "boson_star", "areal"
SIGMA_C = 0.271          # the declared star: the central field of the heaviest ground state
# What the declared star is stated to be, each to the digits the drawings print.
STATED = {"M": (0.633, 3), "omega": (0.853, 3), "alpha_c": (0.686, 3)}
KAUP_MASS = 0.633        # the largest mass of a ground state, in M_Pl^2/m


# The declared star in words, which every drawing of it prints.
INPUT = ("The ground state of a free complex scalar field with the central amplitude $\\sigma(0) = 0.271$, in "
         "units where $G = c = \\mu = 1$, solved from this spacetime's own $G^t{}_t$ and $G^r{}_r$ and the "
         "field's wave equation: the heaviest boson star, of mass $M = 0.633\\,M_P^2/m$ and frequency "
         "$\\omega = 0.853\\,\\mu c$, where $M_P = \\sqrt{\\hbar c/G}$ is the Planck mass.")
SETTINGS = "$\\mu = 1$, so that the unit of every length is $\\hbar/mc$, the reduced Compton wavelength of the boson."


def _load(metric_id, system_id):
    metric = json.loads((build.METRICS_DIR / f"{metric_id}.json").read_text(encoding="utf-8"))
    entry = next(e for e in metric["coordinates"] if e["id"] == system_id)
    reader = vm.Reader(entry["coords"], [p["symbol"] for p in entry.get("parameters", [])],
                       vm.time_coordinates(vm.DIMENSIONS[(metric_id, system_id)], entry["coords"]))
    return entry, reader


@lru_cache(maxsize=None)
def equations(metric_id=METRIC, system_id=AREAL):
    """a', alpha' and sigma'' as functions of (r, a, alpha, sigma, sigma', omega^2), then a'' and
    alpha'' by the chain rule, and the published G^theta_theta less the field's stress across the
    radius as a function of (r, a, a', alpha, alpha', alpha'', sigma, sigma', omega^2)."""
    entry, reader = _load(metric_id, system_id)
    r = reader.symbol["r"]
    fa, falpha = reader.parameters["a"], reader.parameters["alpha"]
    held = dict(zip(entry["coords"], (0, r, sp.pi / 2, 0)))
    angles = {reader.symbol[c]: v for c, v in held.items() if c not in ("t", "r")}
    A, A1, L, L1, L2, S, S1, W = sp.symbols("a a1 alpha alpha1 alpha2 sigma sigma1 omega2")
    ul = entry["einstein_tensor"]["variants"]["ul"]["nonzero"]

    def published(index):
        e = reader(next(c["value"] for c in ul if c["indices"] == [index, index]))
        e = e.subs(sp.Derivative(falpha, (r, 2)), L2).subs({sp.Derivative(falpha, r): L1, sp.Derivative(fa, r): A1})
        return e.subs({falpha: L, fa: A, reader.c: 1}).subs(angles)
    kinetic, gradient, mass = W * S ** 2 / L ** 2, S1 ** 2 / A ** 2, S ** 2
    solved = sp.solve([published("t") + kinetic + gradient + mass,
                       published("r") - (kinetic + gradient - mass)], [A1, L1], dict=True)[0]
    da, dalpha = sp.simplify(solved[A1]), sp.simplify(solved[L1])
    d2sigma = -(2 / r + dalpha / L - da / A) * S1 - A ** 2 * (W / L ** 2 - 1) * S

    def along(e):
        return sp.diff(e, r) + sp.diff(e, A) * da + sp.diff(e, L) * dalpha + sp.diff(e, S) * S1 + sp.diff(e, S1) * d2sigma
    state = (r, A, L, S, S1, W)
    miss = published("\\theta") - (kinetic - gradient - mass)
    return {
        "first": sp.lambdify(state, [da, dalpha, d2sigma], "numpy"),
        "second": sp.lambdify(state, [along(da), along(dalpha)], "numpy"),
        "theta_theta": sp.lambdify((r, A, A1, L, L1, L2, S, S1, W), miss, "numpy"),
    }


class _Tabled:
    """`at`, which the drawings call point by point along every ray: the functions and their first
    two derivatives from a cubic spline through `values` on a grid of TABLE points inside the edge,
    which follows them to parts in 10^10, and from `values` itself from the edge on, where they
    are closed forms."""

    TABLE = 40001

    def _table(self):
        if getattr(self, "_spline", None) is None:
            grid = np.linspace(0.0, self.edge, self.TABLE)
            v = self.values(np.minimum(grid, self.edge * (1 - 1e-15)))
            self._spline = CubicSpline(grid, np.column_stack([v[name][k] for name in self.funcs for k in range(3)]))
        return self._spline

    def at(self, name, x0, x):
        """A declared function and its first two derivatives at chart points (x^0, x)."""
        shape = np.broadcast(np.asarray(x0), np.asarray(x)).shape
        flat = np.broadcast_to(np.asarray(x, dtype=float), shape).ravel()
        column = 3 * self.funcs.index(name)
        inside = flat < self.edge
        out = np.empty((3, flat.size))
        if inside.any():
            out[:, inside] = self._table()(flat[inside])[:, column:column + 3].T
        if (~inside).any():
            out[:, ~inside] = np.array(self.values(flat[~inside])[name])
        return [v.reshape(shape) for v in out]


def _times(p, q):
    """The product of two power series, cut to the length of the first."""
    return np.convolve(p, q)[:len(p)]


def _over(p):
    """The reciprocal of a power series whose first coefficient is not zero."""
    out = np.zeros_like(p)
    out[0] = 1 / p[0]
    for k in range(1, len(p)):
        out[k] = -np.dot(p[1:k + 1], out[k - 1::-1][:k]) / p[0]
    return out


def centre_series(sigma_c, w2, terms=12):
    """a, alpha and sigma about the centre as power series in x = r^2, with alpha(0) = 1.

    With d/dr = 2r d/dx the three equations read

        4x A' = A (-(A^2 - 1) + x E),          E = (w2/L^2 + 1) A^2 S^2 + 4x S'^2,
        4x L' = L ( (A^2 - 1) + x P),          P = (w2/L^2 - 1) A^2 S^2 + 4x S'^2,
        6 S' + 4x S'' = -4x (L'/L - A'/A) S' - A^2 (w2/L^2 - 1) S,

    so the coefficient of x^k on each left hand side, (4k + 2) A_k once the -2 A_k of -(A^2 - 1)
    is brought over, 4k L_k and (k + 1)(4k + 6) S_(k+1), is fixed by the coefficients below it."""
    A, L, S = (np.zeros(terms) for _ in range(3))
    A[0], L[0], S[0] = 1.0, 1.0, sigma_c
    k = np.arange(terms)
    x = np.zeros(terms)
    x[1] = 1.0

    def slope(p):
        return np.append(p[1:] * k[1:], 0.0)
    for order in range(1, terms):
        A2, iL2, dS = _times(A, A), _over(_times(L, L)), slope(S)
        field, grad = _times(_times(A2, _times(S, S)), w2 * iL2), 4 * _times(x, _times(dS, dS))
        mass = _times(A2, _times(S, S))
        E, P = field + mass + grad, field - mass + grad
        curve = A2.copy()
        curve[0] -= 1.0
        rhs_a = _times(A, _times(x, E) - curve)
        A[order] = rhs_a[order] / (4 * order + 2)
        curve[order] += 2 * A[order]
        L[order] = _times(L, curve + _times(x, P))[order] / (4 * order)
        drift = _times(slope(L), _over(L)) - _times(slope(A), _over(A))
        w = w2 * iL2
        w[0] -= 1.0
        rhs_s = -4 * _times(x, _times(drift, dS)) - _times(_times(A2, w), S)
        S[order] = rhs_s[order - 1] / (order * (4 * order + 2))
    return A, L, S


class Star(_Tabled):
    """The ground state with central field sigma_c. `funcs` are the declared functions of the
    areal chart, which `values` and `at` hand back with their first two derivatives along r.

    Within `core` of the centre every function is its power series in r^2, centre_series, since
    there the equations divide differences that vanish as r^2 by r, and the integration starts
    from the series at `core`."""

    funcs = ("alpha", "a")
    core = 0.1

    def __init__(self, sigma_c, metric_id=METRIC, system_id=AREAL):
        self.sigma_c = float(sigma_c)
        self.eq = equations(metric_id, system_id)
        low, high = 0.0, 1.0
        while not self._shoot(high).t_events[0].size:
            low, high = high, 2 * high
        while high - low > 4e-16 * high:
            mid = 0.5 * (low + high)
            if self._shoot(mid).t_events[0].size:
                high = mid
            else:
                low = mid
        self.w2 = low
        self.series = centre_series(self.sigma_c, low)
        self.ivp = self._shoot(low, dense=True)
        if not self.ivp.t_events[1].size:
            raise ArithmeticError(f"boson star: the field of sigma_c = {sigma_c} never turned")
        self.edge = float(self.ivp.t_events[1][0])
        a, alpha, left, _, log_ratio, number = self.ivp.sol(self.edge)
        self.left = float(left)                       # the field at the edge
        self.M = float(0.5 * self.edge * (1 - 1 / a ** 2))
        self.scale = float(np.sqrt(1 - 2 * self.M / self.edge) / alpha)
        self.alpha_c = self.scale                     # the lapse at the centre, where it was 1
        self.omega = float(np.sqrt(low)) * self.scale
        # N in M_Pl^2/m^2: int (omega/alpha) a sigma^2 r^2 dr, with alpha(0) = 1 while integrating.
        self.N = float(np.sqrt(low) * number)
        rho_edge = 0.5 * (self.edge - self.M + np.sqrt(self.edge ** 2 - 2 * self.M * self.edge))
        self.log_shift = float(np.log(rho_edge / self.edge) - log_ratio)

    def _rhs(self, x, y, w2):
        a, alpha, sigma, slope = y[:4]
        da, dalpha, d2sigma = self.eq["first"](x, a, alpha, sigma, slope, w2)
        return [da, dalpha, slope, d2sigma, (a - 1) / x, a * sigma ** 2 * x ** 2 / alpha]

    @staticmethod
    def _centre(series, r):
        """From the series: a, alpha, sigma, sigma', ln(rho/r) less its value at the centre, the
        number inside r over omega, and then (a - 1)/r^2, a'/r, a'', alpha'/r and alpha''."""
        A, L, S = series
        r = np.asarray(r, dtype=float)
        x = r * r
        k = np.arange(len(A))

        def value(p):
            return np.polynomial.polynomial.polyval(x, p)

        def slope(p):
            return np.append(p[1:] * k[1:], 0.0)
        a, alpha, sigma = value(A), value(L), value(S)
        dA, dL = value(slope(A)), value(slope(L))
        log_ratio = value(np.append(0.0, A[1:] / (2 * k[1:])))
        density = _times(_times(A, _times(S, S)), _over(L))
        number = r ** 3 * value(density / (2 * k + 3))
        return (a, alpha, sigma, 2 * r * value(slope(S)), log_ratio, number,
                value(np.append(A[1:], 0.0)), 2 * dA, 2 * dA + 4 * x * value(slope(slope(A))),
                2 * dL, 2 * dL + 4 * x * value(slope(slope(L))))

    def _shoot(self, w2, dense=False):
        def crosses(x, y, w2):
            return y[2]

        def turns(x, y, w2):
            return y[3]
        crosses.terminal = turns.terminal = True
        turns.direction = 1
        start = [float(v) for v in self._centre(centre_series(self.sigma_c, w2), self.core)[:6]]
        return solve_ivp(self._rhs, (self.core, 400.0), start, args=(w2,), events=(crosses, turns), rtol=1e-12,
                         atol=1e-18, method="DOP853", dense_output=dense)

    def inside(self, x):
        """(a, alpha, sigma, sigma', ln(rho/r), (a - 1)/r^2, a'/r) at radii x inside the edge,
        alpha scaled to 1 far away."""
        x = np.atleast_1d(np.asarray(x, dtype=float))
        near = x < self.core
        out = np.empty((7, x.size))
        if near.any():
            c = self._centre(self.series, x[near])
            out[:, near] = np.array([c[0], c[1], c[2], c[3], c[4], c[6], c[7]])
        if (~near).any():
            xo = x[~near]
            a, alpha, sigma, slope, log_ratio, _ = self.ivp.sol(xo)
            da = self.eq["first"](xo, a, alpha, sigma, slope, self.w2)[0]
            out[:, ~near] = np.array([a, alpha, sigma, slope, log_ratio, (a - 1) / xo ** 2, da / xo])
        out[1] *= self.scale
        out[4] += self.log_shift
        return out

    def field(self, x):
        """sigma and its slope at radii x; zero from the edge on, where what is left is dropped."""
        x = np.atleast_1d(np.asarray(x, dtype=float))
        inner = self.inside(np.minimum(x, self.edge))
        return np.where(x < self.edge, inner[2], 0.0), np.where(x < self.edge, inner[3], 0.0)

    def mass(self, x):
        """The mass inside radius x, r (1 - 1/a^2)/2."""
        return 0.5 * np.asarray(x, dtype=float) * (1 - 1 / self.values(x)["a"][0] ** 2)

    def radius(self, fraction=0.99):
        """The radius holding that fraction of the mass, since the star has no surface."""
        return float(brentq(lambda x: self.mass(x)[0] - fraction * self.M, self.core, self.edge))

    def compactness(self):
        """The largest 2M(r)/r, and the radius where it stands."""
        best = minimize_scalar(lambda x: -2 * self.mass(x)[0] / x, bounds=(self.core, self.edge), method="bounded",
                               options={"xatol": 1e-10})
        return float(-best.fun), float(best.x)

    def values(self, x):
        """(value, first, second derivative) of alpha and of a at radii x."""
        x = np.atleast_1d(np.asarray(x, dtype=float))
        out = {name: [np.empty_like(x) for _ in range(3)] for name in self.funcs}
        near, far = x < self.core, x >= self.edge
        between = ~near & ~far
        if near.any():
            xn = x[near]
            c = self._centre(self.series, xn)
            for store, v in zip(out["alpha"], (c[1], xn * c[9], c[10])):
                store[near] = v * self.scale
            for store, v in zip(out["a"], (c[0], xn * c[7], c[8])):
                store[near] = v
        if between.any():
            xi = x[between]
            a, alpha, sigma, slope = self.inside(xi)[:4]
            w2 = self.omega ** 2
            da, dalpha, _ = self.eq["first"](xi, a, alpha, sigma, slope, w2)
            d2a, d2alpha = self.eq["second"](xi, a, alpha, sigma, slope, w2)
            for store, v in zip(out["alpha"], (alpha, dalpha, d2alpha)):
                store[between] = v
            for store, v in zip(out["a"], (a, da, d2a)):
                store[between] = v
        if far.any():
            xo, M = x[far], self.M
            f = 1 - 2 * M / xo
            for store, v in zip(out["alpha"], (np.sqrt(f), M / (xo ** 2 * np.sqrt(f)),
                                               -M * (2 * xo - 3 * M) / (xo ** 4 * f ** 1.5))):
                store[far] = v
            for store, v in zip(out["a"], (f ** -0.5, -M / (xo ** 2 * f ** 1.5),
                                           M * (2 * xo - M) / (xo ** 4 * f ** 2.5))):
                store[far] = v
        return out

    def theta_theta(self, x):
        """The published G^theta_theta less the field's stress across the radius along the
        solution, against the central energy density: zero if the star solves the one field
        equation its construction did not use."""
        x = np.asarray(x, dtype=float)
        v = self.values(x)
        sigma, slope = self.field(x)
        miss = self.eq["theta_theta"](x, v["a"][0], v["a"][1], v["alpha"][0], v["alpha"][1], v["alpha"][2],
                                      sigma, slope, self.omega ** 2)
        return miss / ((self.omega ** 2 / self.alpha_c ** 2 + 1) * self.sigma_c ** 2)

    def stated(self):
        """What of STATED the star is not, to the digits stated."""
        mine = {"M": self.M, "omega": self.omega, "alpha_c": self.alpha_c}
        return {k: mine[k] for k, (value, digits) in STATED.items() if round(mine[k], digits) != value}

    def isotropic(self):
        return Isotropic(self)


class Isotropic(_Tabled):
    """The star in the isotropic radius R: alpha and psi, each with its first two derivatives
    along R. Inside the edge R = r exp(L(r)) with L' = (a - 1)/r, inverted by Newton's method, and
    psi = exp(-L/2); outside it is Schwarzschild's isotropic chart, psi = 1 + M/(2R). The
    derivatives are written in (a - 1)/r^2 and a'/r, which stay finite at the centre:
    psi' = -g psi^5 R/(2a) with g = (a - 1)/r^2, since 1/a - 1 = -g r^2/a and r = psi^2 R."""

    funcs = ("alpha", "psi")

    def __init__(self, star):
        self.star = star
        self.M = star.M
        self.grid_r = np.concatenate([[0.0], np.geomspace(1e-6, star.edge, 4000)])
        self.grid_rho = self.rho_of(self.grid_r)
        self.edge = float(self.grid_rho[-1])

    def rho_of(self, r):
        """The isotropic radius of the areal radius r."""
        r = np.atleast_1d(np.asarray(r, dtype=float))
        M, edge = self.M, self.star.edge
        inner = r * np.exp(self.star.inside(np.minimum(r, edge))[4])
        with np.errstate(invalid="ignore"):
            outer = 0.5 * (r - M + np.sqrt(r * r - 2 * M * r))
        return np.where(r < edge, inner, outer)

    def r_of(self, rho):
        """The areal radius of the isotropic radius rho."""
        rho = np.atleast_1d(np.asarray(rho, dtype=float))
        M = self.M
        r = np.interp(rho, self.grid_rho, self.grid_r)
        for _ in range(6):
            inner = self.star.inside(np.minimum(r, self.star.edge))
            r = r - (r * np.exp(inner[4]) - rho) / (inner[0] * np.exp(inner[4]))
        with np.errstate(divide="ignore", invalid="ignore"):
            outer = rho * (1 + M / (2 * rho)) ** 2
        return np.where(rho < self.edge, r, outer)

    def values(self, rho):
        """(value, first, second derivative along R) of alpha and of psi at isotropic radii rho."""
        rho = np.atleast_1d(np.asarray(rho, dtype=float))
        out = {name: [np.empty_like(rho) for _ in range(3)] for name in self.funcs}
        inside = rho < self.edge
        if inside.any():
            p = rho[inside]
            r = self.r_of(p)
            a, _, _, _, log_ratio, g, da_r = self.star.inside(r)
            alpha, dalpha, d2alpha = self.star.values(r)["alpha"]
            psi = np.exp(-0.5 * log_ratio)
            rate = psi ** 2 / a                       # dr/dR
            # d(rate)/dR, with (1/a - 1)/R = -g psi^4 R/a and r a'/R = psi^4 R a'/r.
            drate = -rate * g * psi ** 4 * p / a - psi ** 6 * p * da_r / a ** 3
            dpsi = -g * psi ** 5 * p / (2 * a)
            d2psi = (-dpsi * g * psi ** 4 * p / (2 * a) + g * psi ** 5 / (2 * a)
                     - psi ** 5 * da_r / (2 * a ** 3))
            for store, value in zip(out["alpha"], (alpha, dalpha * rate, d2alpha * rate ** 2 + dalpha * drate)):
                store[inside] = value
            for store, value in zip(out["psi"], (psi, dpsi, d2psi)):
                store[inside] = value
        if (~inside).any():
            p, M = rho[~inside], self.M
            u = M / (2 * p)
            for store, value in zip(out["alpha"], ((1 - u) / (1 + u), 2 * u / (p * (1 + u) ** 2),
                                                   -4 * u / (p * p * (1 + u) ** 3))):
                store[~inside] = value
            for store, value in zip(out["psi"], (1 + u, -u / p, 2 * u / (p * p))):
                store[~inside] = value
        return out


@lru_cache(maxsize=None)
def star(sigma_c=SIGMA_C):
    return Star(sigma_c)


def heaviest(step=0.02):
    """The central field whose ground state is heaviest and its mass, Kaup's limit, from the
    parabola through the declared star and one on either side of it."""
    fields = np.array([SIGMA_C - step, SIGMA_C, SIGMA_C + step])
    c2, c1, c0 = np.polyfit(fields, [star(round(float(f), 6)).M for f in fields], 2)
    top = -c1 / (2 * c2)
    return float(top), float(c0 + c1 * top + c2 * top * top)


if __name__ == "__main__":
    s = star()
    print(f"declared, sigma_c = {s.sigma_c}: M = {s.M:.6f}, omega = {s.omega:.6f}, alpha_c = {s.alpha_c:.6f}, "
          f"N = {s.N:.6f}, R_99 = {s.radius():.4f}, edge = {s.edge:.2f} with sigma = {s.left:.1e}, "
          f"largest 2m/r = {s.compactness()[0]:.4f} at r = {s.compactness()[1]:.3f}")
    x = np.linspace(0.05, 0.98 * s.edge, 400)
    print(f"G^theta_theta misses by {np.max(np.abs(s.theta_theta(x))):.1e}")
    print("heaviest: sigma_c = {:.4f}, M = {:.6f}".format(*heaviest()))
