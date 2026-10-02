#!/usr/bin/env python3
"""The Gregory-Laflamme instability of the black string in five dimensions, in units of r_s.

A perturbation h_ab = exp(Omega t + i k z) H_ab(r) of Schwarzschild's black hole times a line,
spherically symmetric and with no components along z, solves the linearised vacuum equations
when Gregory's first order system holds (R. Gregory, "The Gregory-Laflamme instability",
arXiv:1107.5821, eqs. 1.21 to 1.23): with V = 1 - 1/r, H = H_tr, H_+- = H_tt/V +- V H_rr and
K = r^2 H_-/2 on the sphere,

    H_+ = H_- (2 r^2 Omega^2 + r^2 k^2 V - (1 - V^2)/2)/(V D) - (r H/Omega)(4 Omega^2 + k^2 (1 - 3V))/D,
    H'  = Omega (H_+ + H_-)/(2V) - (1 + V) H/(r V),
    H_-' = k^2 H/Omega + H_+/r + (1 - 5V) H_-/(2 r V),          D = r^2 k^2 + 1 - V.

linearised_ricci() checks that in sympy: every component of the Ricci tensor of g + epsilon h
vanishes at first order once the three equations are put in, with no gauge condition assumed.

A mode is unstable when it is regular on the future horizon, H ~ (Omega - 1/2)(r - 1)^(Omega - 1)
and H_- ~ (k^2/Omega + 2)(r - 1)^Omega, and decays as exp(-sqrt(Omega^2 + k^2) r) far away.
mismatch() integrates the regular solution out from the horizon and the decaying one in from
far away, in H/Omega and H_-, which leaves Omega^2 alone in the equations, and returns the
Wronskian of the two where they meet; growth_rate(k) is its zero. Gregory and Laflamme's
threshold is k r_s = 0.876, and nothing grows at shorter wavelengths.

    /tmp/mfs-venv/bin/python _tools/derivations/gregory_laflamme.py
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar

THRESHOLD = 0.876        # k r_s where the growth rate vanishes, Gregory and Laflamme's value


def rhs(r, y, Om, k):
    Hh, Hm = y
    V = 1 - 1 / r
    D = r * r * k * k + 1 - V
    Hp = Hm * (2 * r * r * Om * Om + r * r * k * k * V - (1 - V * V) / 2) / (V * D) - r * Hh * (4 * Om * Om + k * k * (1 - 3 * V)) / D
    return [(Hp + Hm) / (2 * V) - (1 + V) * Hh / (r * V), k * k * Hh + Hp / r + (1 - 5 * V) * Hm / (2 * r * V)]


def mismatch(Om, k, meet=2.5, x0=1e-6):
    """The Wronskian at r = meet of the solution regular on the future horizon and the one
    that decays far away, each normalised; it vanishes for an unstable mode."""
    far = 1 + 30 / np.hypot(Om, k)
    start = [x0 ** (Om - 1) * (Om - 0.5), x0 ** Om * (k * k + 2 * Om)]
    left = solve_ivp(rhs, [1 + x0, meet], start, args=(Om, k), rtol=1e-10, atol=0, method="LSODA")
    # Far away the system has constant coefficients, and the decaying solution is the
    # eigenvector of the negative eigenvalue.
    J = np.array([rhs(far, [1, 0], Om, k), rhs(far, [0, 1], Om, k)]).T
    w, v = np.linalg.eig(J)
    right = solve_ivp(rhs, [far, meet], np.real(v[:, np.argmin(np.real(w))]), args=(Om, k), rtol=1e-10, atol=0, method="LSODA")
    a, b = left.y[:, -1], right.y[:, -1]
    return (a[0] * b[1] - a[1] * b[0]) / (np.hypot(*a) * np.hypot(*b))


def growth_rate(k, lo=2e-3, hi=0.15, n=40):
    """Omega r_s/c of the unstable mode of wavenumber k r_s, or None where no mode grows."""
    grid = np.linspace(lo, hi, n)
    values = [mismatch(Om, k) for Om in grid]
    for i in range(n - 1):
        if values[i] * values[i + 1] < 0:
            return brentq(mismatch, grid[i], grid[i + 1], args=(k,), xtol=1e-12)
    return None


def threshold(Om=1e-5):
    """The wavenumber at which a mode of the small growth rate Om is found, which tends to the
    threshold as Om tends to zero."""
    return brentq(lambda k: mismatch(Om, k), 0.80, 0.90, xtol=1e-12)


def fastest():
    """(k r_s, Omega r_s/c) of the mode that grows fastest."""
    best = minimize_scalar(lambda k: -growth_rate(k, 0.07, 0.11, 8), bounds=(0.28, 0.48), method="bounded",
                           options={"xatol": 1e-7})
    return best.x, -best.fun


def linearised_ricci():
    """Every first order component of the Ricci tensor of the perturbed string, with Gregory's
    three equations put in, simplified: a dictionary that holds only zeros when they solve the
    linearised vacuum equations."""
    import sympy as sp
    t, r, th, ph, z, Om, k = sp.symbols("t r theta phi z Omega k", positive=True)
    eps = sp.Symbol("epsilon")
    V = 1 - 1 / r
    X = [t, r, th, ph, z]
    Hm, H, Hp = sp.Function("Hm")(r), sp.Function("H")(r), sp.Function("Hp")(r)
    E = sp.exp(Om * t + sp.I * k * z)
    g0 = sp.diag(-V, 1 / V, r ** 2, r ** 2 * sp.sin(th) ** 2, 1)
    h = sp.zeros(5)
    h[0, 0], h[1, 1] = V * (Hp + Hm) / 2 * E, (Hp - Hm) / (2 * V) * E
    h[0, 1] = h[1, 0] = H * E
    h[2, 2], h[3, 3] = r ** 2 * Hm / 2 * E, r ** 2 * Hm / 2 * sp.sin(th) ** 2 * E
    g, ginv0 = g0 + eps * h, g0.inv()
    ginv = ginv0 - eps * ginv0 * h * ginv0

    def first(e):
        e = sp.expand(e)
        return e.coeff(eps, 0) + eps * e.coeff(eps, 1)
    G = [[[first(sum(ginv[a, d] * (sp.diff(g[d, b], X[c]) + sp.diff(g[d, c], X[b]) - sp.diff(g[b, c], X[d]))
                     for d in range(5)) / 2) for c in range(5)] for b in range(5)] for a in range(5)]
    D = r ** 2 * k ** 2 + 1 - V
    Hp_is = Hm * (2 * r ** 2 * Om ** 2 + r ** 2 * k ** 2 * V - (1 - V ** 2) / 2) / (V * D) - (r * H / Om) * (4 * Om ** 2 + k ** 2 * (1 - 3 * V)) / D
    dH = Om * (Hp_is + Hm) / (2 * V) - (1 + V) * H / (r * V)
    dHm = k ** 2 * H / Om + Hp_is / r + (1 - 5 * V) * Hm / (2 * r * V)
    out = {}
    for b in range(5):
        for c in range(b, 5):
            e = 0
            for a in range(5):
                e += sp.diff(G[a][b][c], X[a]) - sp.diff(G[a][b][a], X[c])
                for d in range(5):
                    e += G[a][a][d] * G[d][b][c] - G[a][c][d] * G[d][b][a]
            e = sp.expand(e).coeff(eps, 1).subs(Hp, Hp_is).doit()
            for _ in range(3):
                e = e.subs(sp.Derivative(H, (r, 2)), sp.diff(dH, r)).subs(sp.Derivative(Hm, (r, 2)), sp.diff(dHm, r))
                e = e.subs(sp.Derivative(H, r), dH).subs(sp.Derivative(Hm, r), dHm)
            out[(str(X[b]), str(X[c]))] = sp.simplify(e / E)
    return out


if __name__ == "__main__":
    wrong = {slot: value for slot, value in linearised_ricci().items() if value != 0}
    print("Gregory's system solves the linearised vacuum equations:", "FAILED " + str(wrong) if wrong else "ok")
    k, Om = fastest()
    print(f"the fastest mode: k r_s = {k:.4f}, wavelength {2 * np.pi / k:.2f} r_s, Omega r_s/c = {Om:.4f}")
    print(f"the threshold: k r_s = {threshold():.4f}, wavelength {2 * np.pi / threshold():.2f} r_s")
    print(f"on a circle of length 10 r_s: Omega r_s/c = {growth_rate(2 * np.pi / 10):.4f}")
    raise SystemExit(1 if wrong else 0)
