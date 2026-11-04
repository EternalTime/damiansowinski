"""A local flat-interior bounce and a separate exact equilibrium of the cylindrical shell."""

from functools import lru_cache
import math

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

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


@lru_cache(maxsize=1)
def radius_interpolant():
    tau = np.linspace(0, 5, 2001)
    R, V, T = trajectory()(tau)
    interpolation = CubicSpline(T, R)
    probe = np.linspace(.001, 4.999, 701)
    reference = trajectory()(probe)
    error = np.max(np.abs(interpolation(reference[2])-reference[0]))
    if error > 2e-12:
        raise AssertionError(f"interior-time interpolation misses by {error}")
    return interpolation


def radius_at_time(time):
    return radius_interpolant()(np.abs(np.asarray(time, dtype=float)))


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


def embedding(ck, src, em):
    """Three spatial slices read from the stored metrics, with their junction circles marked."""
    views = []
    definitions = [
        ("interior_turn", "The interior at the turn", "interior", "\\rho", {"t": 0, "z": 0},
         BOUNCE, {"psi": "0", "gamma": "0", "R": "2"}, 0, 2, 0,
         (("axis", "the regular axis"), ("edge", "the shell")), [(2, "r", None)], "$t=0$"),
        ("shell_turn", "The shell at the turn", "shell", "z", {"tau": 0},
         BOUNCE, {"psi": "0", "R": "2"}, -2, 2, -2,
         (("edge", "the shell continues along the axis"), ("edge", "the shell continues along the axis")),
         [(-2, "r", None), (0, "r", None), (2, "r", None)], "$\\tau=0$"),
        ("equilibrium", "The equilibrium exterior", "exterior", "\\rho", {"t": 0, "z": 0},
         EQUILIBRIUM, STATIC, 2, 5, 0,
         (("edge", "the equilibrium shell"), ("edge", "the exterior continues outward")),
         [(2, "r", None), (5, "r", None)], "$t=0$"),
    ]
    for vid, label, system, coordinate, fixed, params, functions, lo, hi, z0, ends, marks, moment in definitions:
        sl = em.Slice(src, METRIC, system, coordinate, "\\phi", fixed, params, functions=functions)
        size = 4 if vid != "equilibrium" else 12
        piece = em.Piece(vid, "sheet", sl, lo, hi, z0, ends=ends, marks=marks, size=size)
        ck.isometry(vid, piece)
        if vid == "interior_turn":
            ck.radius("flat interior circles have their canonical radius", piece, lambda r: r, size)
            ck.form("flat interior is a plane", piece, lambda r: np.zeros_like(r), size)
        elif vid == "shell_turn":
            ck.radius("shell circles have radius two", piece, lambda z: np.zeros_like(z)+2, size)
        else:
            ck.radius("equilibrium circumferential radius", piece, lambda r: r*np.sqrt(r/2), size)
        surface = em.Surface([piece])
        fig = em.figure_of([surface], {"sheet": "cover"}, size, camera=em.Camera(-90, 40))
        fig.legend("fill", "cover", "the spatial surface")
        fig.legend("line", "r", "the marked circles")
        setting = "$\\alpha=1$, $\\lambda=0.08$, $R=2\\alpha$." if vid != "equilibrium" else "$\\alpha=1$, $R_0=2\\alpha$."
        views.append(em.view(vid, label, "$\\alpha$", [surface], fig.done(), settings=setting))
    return views


def register_embedding(em):
    em.DRAWN[METRIC] = lambda ck, src: embedding(ck, src, em)
    em.CAPTIONS.update({
        (METRIC, "interior_turn"): ["The flat interior disc at the turning point ($t=0$, $z=0$), bounded by the shell."],
        (METRIC, "shell_turn"): ["The material cylinder at its minimum local radius ($\\tau=0$, $|z|\\leq2\\alpha$), with its central ring marked."],
        (METRIC, "equilibrium"): ["The static exterior slice ($t=0$, $2\\alpha\\leq\\rho\\leq5\\alpha$), with the material boundary at its inner circle."],
    })


def register_slices(sl):
    def marks(system):
        out = []
        moments = [sl.moments(METRIC, vid, label=label)[0] for vid, label in
                   [("interior_turn", "$t=0$"), ("shell_turn", "$\\tau=0$"), ("equilibrium", "$t=0$")]]
        for m in moments:
            vid = m.view["id"]
            if system == "exterior":
                if vid == "equilibrium":
                    out.append(sl.Mark(m, sl.along(0, *m.reach("exterior", "\\rho"))))
            elif vid == "interior_turn":
                out.append(sl.Mark(m, sl.along(0, *m.reach("interior", "\\rho")))
                           if system == "interior" else sl.Mark(m, points=[(0, 0)]))
            elif vid == "shell_turn":
                out.append(sl.Mark(m, sl.along(0, *m.reach("shell", "z")))
                           if system == "shell" else sl.Mark(m, points=[(0, 2)]))
        return out
    for system, view in [("interior", "bounce"), ("exterior", "equilibrium"), ("shell", "axial")]:
        sl.FLAT[(METRIC, system, view)] = lambda system=system: marks(system)
    sl.FLAT_METRICS.add(METRIC)


def curvature(ck, src, cv):
    """Read the smooth equilibrium bulk scalars; a squared shell delta is undefined."""
    inside = cv.Published(src, METRIC, "interior", "kretschmann", ["rho"], [],
                          functions={"psi": "0", "gamma": "0", "R": "2"})
    outside = cv.Published(src, METRIC, "exterior", "kretschmann", ["rho"], [], functions=STATIC)
    ricci = cv.Published(src, METRIC, "exterior", "ricci_scalar", ["rho"], [], functions=STATIC)
    ck("equilibrium bulk Ricci scalar is zero", ricci.expr == 0)
    ck("restricted flat interior Kretschmann scalar is zero", inside.expr == 0)
    rs = np.linspace(2.00001, 5, 200)
    values = outside(rs)
    ck("equilibrium exterior curvature is positive and decreases outward",
       np.all(values > 0) and np.all(np.diff(values) < 0))
    ck("equilibrium exterior radial power is minus seven",
       abs(float(outside(4)/outside(2))-2**-7) < 1e-14)
    axis = cv.axis("$\\rho/\\alpha$", 0, 5)
    scale = {"kind": "linear", "from": 0, "to": .05}
    aspect = 1.6
    cv.shape_checks(ck, METRIC, [aspect])
    curves = (cv.line_curves(ck, METRIC+" interior", "K", "K", inside, "\\rho", axis, scale,
                             cv.PLOT_W/aspect, {}, within=(0, 1.99999))
              + cv.line_curves(ck, METRIC+" exterior", "K", "K", outside, "\\rho", axis, scale,
                               cv.PLOT_W/aspect, {}, within=(2.00001, 5)))
    return [{"id": "equilibrium", "label": "The equilibrium bulk", "kind": "line", "x": axis,
             "panels": [{"id": "K", "name": "$K\\alpha^4$", "aspect": aspect,
                         "scale": scale, "ticks": cv.value_ticks(scale)}],
             "curves": curves, "marks": [{"class": "front", "x": 2}],
             "legend": [["line", "K", "the smooth bulk curvature"], ["line", "front", "the shell"]],
             "caption": ["The Kretschmann scalar $K$ in the smooth equilibrium bulk ($R_0=2\\alpha$), zero inside the shell and decreasing as $\\rho^{-7}$ outside."],
             "settings": "$\\alpha=1$, $R_0=2\\alpha$.", "clip": "Linear ($0\\leq K\\alpha^4\\leq0.05$)."}]


def register_curvature(cv):
    cv.DRAWN[METRIC] = lambda ck, src: curvature(ck, src, cv)
