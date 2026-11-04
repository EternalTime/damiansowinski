"""A local flat-interior bounce and a separate exact equilibrium of the cylindrical shell."""

from functools import lru_cache
import math

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from rotation_halts import METRIC, acceleration

BOUNCE = {"alpha": 1, "lambda": "2/25"}
EQUILIBRIUM = {"alpha": 1, "lambda": "sqrt(5)/18"}
STATIC = {"psi": "-log(rho/2)/2", "gamma": "2*log(3/2)+log(rho/2)/4", "R": "2"}


@lru_cache(maxsize=1)
def trajectory():
    def rhs(tau, state):
        R, V, _ = state
        return V, acceleration(R, V, 1.0, .08), math.sqrt(1+V*V)
    sol = solve_ivp(rhs, (0, 5), (2, 0, 0), method="DOP853", rtol=2e-13, atol=2e-14,
                    dense_output=True, max_step=.04)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.sol


def proper_values(tau):
    t = np.asarray(tau, dtype=float)
    values = trajectory()(np.abs(t))
    return values[0], values[1]*np.sign(t)


def radius_at_time(time):
    def scalar(T):
        tau = brentq(lambda s: float(trajectory()(s)[2])-abs(T), 0, 5, xtol=1e-13)
        return float(trajectory()(tau)[0])
    return np.vectorize(scalar, otypes=[float])(np.asarray(time, dtype=float))


class rh_R(sp.Function):
    nargs = 1
    is_real = True
    _imp_ = staticmethod(lambda t: proper_values(t)[0])

    def fdiff(self, argindex=1):
        return rh_V(self.args[0])


class rh_V(sp.Function):
    nargs = 1
    is_real = True
    _imp_ = staticmethod(lambda t: proper_values(t)[1])

    def fdiff(self, argindex=1):
        return rh_A(self.args[0])


class rh_A(sp.Function):
    nargs = 1
    is_real = True
    _imp_ = staticmethod(lambda t: acceleration(*proper_values(t), 1.0, .08))


class rh_RT(sp.Function):
    nargs = 1
    is_real = True
    _imp_ = staticmethod(radius_at_time)


def register_null(nr):
    nr.DECLARED_FUNCTIONS.update({c.__name__: c for c in (rh_R, rh_V, rh_A, rh_RT)})
    specs = [
        nr.Diagram(METRIC, "interior", "bounce", "The flat interior", ("t", "\\rho"),
                   (-.1, 2.5, -2, 2), "$\\rho/\\alpha$", "$t/\\alpha$", BOUNCE,
                   {"phi": 0, "z": 0}, functions={"psi": "0", "gamma": "0", "R": "rh_RT(t)"},
                   where="rh_RT(t)-rho", curves=(("shell", "rh_RT(t)", "the shell"),),
                   points=(("mark", (0, 2), "the turning point"),),
                   input="$R(0)=2\\alpha$, $dR/d\\tau(0)=0$, $\\psi_-=\\gamma_-=0$."),
        nr.Diagram(METRIC, "exterior", "equilibrium", "The equilibrium exterior", ("t", "\\rho"),
                   (1.7, 5.3, -2, 2), "$\\rho/\\alpha$", "$t/\\alpha$", EQUILIBRIUM,
                   {"phi": 0, "z": 0}, functions=STATIC, where="rho-2",
                   curves=(("shell", "2", "the equilibrium shell"),),
                   input="$R_0=2\\alpha$, $\\psi=-\\frac12\\ln(\\rho/R_0)$."),
        nr.Diagram(METRIC, "shell", "axial", "Along the shell", ("\\tau", "z"),
                   (-2.3, 2.3, -2.3, 2.3), "$z/\\alpha$", "$\\tau/\\alpha$", BOUNCE,
                   {"phi": 0}, tau="tau", functions={"psi": "0", "R": "rh_R(tau)"},
                   families=("toward negative $z$", "toward positive $z$"),
                   input="$R(0)=2\\alpha$, $dR/d\\tau(0)=0$, $\\psi_\\Sigma=0$."),
    ]
    nr.DIAGRAMS.extend(specs)
    nr.CAPTIONS.update({
        (METRIC, "interior", "bounce"): ["Radial light rays inside the shell through its turning point ($t=0$), with the moving material boundary."],
        (METRIC, "exterior", "equilibrium"): ["Radial light rays outside the static matched shell ($R_0=2\\alpha$), with an independent exterior clock."],
        (METRIC, "shell", "axial"): ["Axial light rays on the material surface through the turn ($\\phi=0$), with the shell's proper time vertically."],
    })
    nr.CLOSED_FORMS.update({
        (METRIC, "interior", "bounce"): (lambda t, r: t+r, lambda t, r: t-r, None),
        (METRIC, "exterior", "equilibrium"): (lambda t, r: t+r, lambda t, r: t-r, None),
        (METRIC, "shell", "axial"): (lambda t, z: t+z, lambda t, z: t-z, None),
    })
