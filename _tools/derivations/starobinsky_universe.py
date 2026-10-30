"""The homogeneous semiclassical geometry of Starobinsky (1980), Eqs. (4), (7), (12).

All times in the numerical solution are H_0 tau and all lengths are H_0^-1.
The numerical initial displacement is a linear approximation at the de Sitter saddle.
The evolution retains the full homogeneous equation, including the anomaly term.
"""
from functools import lru_cache

import numpy as np
import sympy as sp

import chart_printer as cp
import verify_metrics as vm

ID = "starobinsky_universe"
CHART_IDS = ("cosmic", "conformal", "closed_de_sitter", "flat_de_sitter", "open_de_sitter")
ANGLES = r"\left(d\theta^2 + \sin^2\theta\,d\phi^2\right)"


def chart_spec(system):
    """Two evolving FLRW charts and the three exact special solutions of Eq. (7)."""
    evolving = system in ("cosmic", "conformal")
    first = r"\eta" if system == "conformal" else r"\tau"
    radial = r"\chi" if system in ("closed_de_sitter", "open_de_sitter") else "r"
    coords = [first, radial, r"\theta", r"\phi"]
    if evolving:
        parameters = [f"a = a({first})", "k", "H_0", "M"]
        space = r"\dfrac{dr^2}{1-kr^2} + r^2" + ANGLES
        line = (r"ds^2 = a^2\left(-d\eta^2 + " + space + r"\right)" if system == "conformal"
                else r"ds^2 = -d\tau^2 + a^2\left(" + space + r"\right)")
        domains = [first + r"\in I", r"a>0", r"r\in[0,\infty)\;\text{for}\;k\le0",
                   r"r\in[0,1/\sqrt{k})\;\text{for}\;k>0"]
    else:
        parameters = ["H_0"]
        if system == "flat_de_sitter":
            space = r"dr^2 + r^2" + ANGLES
            factor = r"e^{2H_0\tau}"
        else:
            trig = r"\sin" if system == "closed_de_sitter" else r"\sinh"
            scale = r"\cosh" if system == "closed_de_sitter" else r"\sinh"
            space = r"d\chi^2 + " + trig + r"^2\chi" + ANGLES
            factor = scale + r"^2(H_0\tau)/H_0^2"
        line = r"ds^2 = -d\tau^2 + " + factor + r"\left(" + space + r"\right)"
        domains = [r"\tau\in(0,\infty)" if system == "open_de_sitter" else r"\tau\in(-\infty,\infty)",
                   (r"\chi\in[0,\pi]" if system == "closed_de_sitter" else radial + r"\in[0,\infty)")]
    domains += [r"\theta\in[0,\pi]", r"\phi\in[0,2\pi)"]
    probe = vm.Reader(coords, parameters, ())
    x = probe.symbol[first]
    pretty = cp.hyperbolic(probe.parameters["H_0"] * x) if system in ("closed_de_sitter", "open_de_sitter") else sp.factor
    return {"metric_id": ID,
            "system": {"id": system, "name": {"cosmic": "Cosmic Time", "conformal": "Conformal Time",
                       "closed_de_sitter": "Closed de Sitter", "flat_de_sitter": "Flat de Sitter",
                       "open_de_sitter": "Open de Sitter"}[system], "coords": coords, "domains": domains,
                       "parameters": parameters, "line_element": line},
            "chart_line_element": line, "printer": {"primed": ["a"]}, "pretty": pretty,
            "check": lambda chart: chart_check(chart, system)}


def chart_check(chart, system):
    """Check the independent FLRW Ricci formula at regular points, with no symbolic determinant."""
    x, r, theta, phi = chart.symbols
    p = chart.reader.parameters
    g = chart.geo.g
    ric = chart.geo.ricci_ll()
    evolving = system in ("cosmic", "conformal")
    for time, radius, curvature in ((sp.Rational(2, 5), sp.Rational(1, 3), -1),
                                    (sp.Rational(7, 5), sp.Rational(1, 4), 1)):
        point = {x: time, r: radius, theta: sp.Rational(4, 5), phi: sp.Rational(1, 5),
                 p["H_0"]: sp.Rational(3, 2), chart.reader.c: 1}
        if evolving:
            a = p["a"]
            values = {a: 2, sp.diff(a, x): sp.Rational(3, 7), sp.diff(a, x, 2): sp.Rational(-1, 9)}
            point.update({p["k"]: curvature, p["M"]: sp.Rational(1, 2)})
            A, B, C = (values[a], values[sp.diff(a, x)], values[sp.diff(a, x, 2)])
            if system == "cosmic":
                temporal, spatial = -3*C/A, C/A+2*(B*B+curvature)/A**2
            else:
                temporal, spatial = -3*(C/A-B*B/A**2), (C/A+B*B/A**2+2*curvature)/A**2
        else:
            values = {}
            temporal, spatial = -3*p["H_0"]**2, 3*p["H_0"]**2
        for i in range(4):
            for j in range(4):
                expected = temporal if i == j == 0 else spatial*g[i,j] if i == j else 0
                delta = (sp.sympify(vm._at(ric, (i,j))) - expected).subs(values).subs(point)
                if abs(complex(sp.N(delta, 30))) > 1e-24:
                    raise AssertionError(f"{ID}/{system}: FLRW Ricci component {i}{j} at {time}")


MASS = 0.5
DISPLACEMENT = 1e-4
EXPONENT = np.sqrt(2.25 + MASS**2) - 1.5


def trace_rhs(_, state, mass=MASS):
    """Regular trace equation, a'/a=h, h'=R/6-2h^2, R'=S, for k=0 and H_0=1."""
    N, h, R, S, eta = state
    Q = h*h
    return [h, R/6-2*Q, S, -3*h*S-mass**2*R+12*mass**2*Q*(R/6-Q), np.exp(-N)]


def constraint(h, R, S, mass=MASS, curvature=0.0):
    """M^2 times Eq. (4), using R=6(h'+2h^2+k/a^2), with H_0=1."""
    Q = h*h+curvature
    return mass**2*(Q-Q*Q)+Q*R/3-R*R/36+h*S/3


class Trajectory:
    """Full trace evolution with the exact constraint at the finite initial point.

    Time and a are normalized at h=1/sqrt(2): tau=0 and a=1 there.
    The second-order unstable-manifold expansion is the only initialization approximation.
    The sample extends 160 H_0^-1 beyond that crossing, through scalaron oscillations.
    """

    def __init__(self, displacement=DISPLACEMENT, tolerance=2e-12):
        from scipy.integrate import solve_ivp
        self.displacement = displacement
        q = EXPONENT
        b = q*(4*q+3)/(6*(q+1))
        h = 1-displacement
        hp = -q*displacement+b*displacement**2
        R = 6*(hp+2*h*h)
        S = (R*R/12-h*h*R-3*MASS**2*(h*h-h**4))/h
        initial = [0.0, h, R, S, 0.0]

        def crossing(_, state):
            return state[1]-1/np.sqrt(2)
        crossing.terminal = True
        first = solve_ivp(trace_rhs, (0, 250), initial, method="DOP853", events=crossing,
                          dense_output=True, rtol=tolerance, atol=tolerance*0.01, max_step=0.25)
        if not first.success or not len(first.t_events[0]):
            raise RuntimeError("Starobinsky trajectory did not reach h=1/sqrt(2)")
        self.shift = first.t_events[0][0]
        middle = first.y_events[0][0]
        self.N0, self.eta0 = middle[0], middle[4]
        self.scale = np.exp(self.N0)
        # Anchor eta at the crossing: subtracting two early-time conformal integrals
        # would lose all late-time digits after the long near-de Sitter expansion.
        self.past_eta = solve_ivp(lambda t, _: [np.exp(self.N0-first.sol(t)[0])],
                                  (self.shift, 0), [0.0], method="DOP853", dense_output=True,
                                  rtol=tolerance, atol=tolerance*0.01, max_step=0.2)
        normalized = [0.0, middle[1], middle[2], middle[3], 0.0]
        last = solve_ivp(trace_rhs, (0, 160), normalized, method="DOP853",
                         dense_output=True, rtol=tolerance, atol=tolerance*0.01, max_step=0.2)
        if not last.success:
            raise RuntimeError(last.message)
        self.first, self.last = first, last
        self.lo, self.hi = -self.shift, 160.0

    def state(self, time):
        time = np.asarray(time, dtype=float)
        if np.any((time < self.lo) | (time > self.hi)):
            raise ValueError(f"trajectory time outside [{self.lo}, {self.hi}]")
        native = time+self.shift
        left = self.first.sol(np.minimum(native, self.shift))
        left[0] -= self.N0
        left[4] = self.past_eta.sol(np.minimum(native, self.shift))[0]
        right = self.last.sol(np.maximum(time, 0))
        state = np.where(native <= self.shift, left, right)
        return state

    def a(self, time, order=0):
        N, h, R, _, _ = self.state(time)
        a = np.exp(N)
        return a if order == 0 else a*h if order == 1 else a*(R/6-h*h)

    def eta(self, time):
        """Conformal time from quadrature of 1/a in the same evolution."""
        return self.state(time)[4]

    def time(self, eta):
        from scipy.optimize import brentq
        values = np.asarray(eta, dtype=float)
        result = [brentq(lambda t: float(self.eta(t))-e, self.lo, self.hi, xtol=1e-11)
                  for e in values.ravel()]
        return np.asarray(result).reshape(values.shape)

    def conformal_a(self, eta, order=0):
        time = self.time(eta)
        a, h, R = self.a(time), self.state(time)[1], self.state(time)[2]
        return a if order == 0 else a*a*h if order == 1 else a**3*R/6


@lru_cache(maxsize=1)
def trajectory():
    return Trajectory()
