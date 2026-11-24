"""The drawings of the black holes of anti-de Sitter space with three charges, each chart at one member,
in units of the radius r_+ of its event horizon.

The spherical and continued charts are drawn at q_1 = 4, q_2 = 2, q_3 = 1 and L = 3, with mu = 13/3, which
puts the event horizon at r_+ = 1: in x = r^2, r^4 f L^2 = (x - 1)(x^2 + 17x - 8), so the inner horizon is
r_-^2 = (sqrt(321) - 17)/2 = 0.4582. In Myers and Tafjord's rho^2 = r^2 + q_3 the horizons are rho_+ = sqrt 2
and rho_- = 1.2075, the radius r = 0 is the regular sphere rho = 1, and the singularity is rho = 0, where
D = q_3(q_3 + mu) > 0, so it is timelike. The flat chart and the brane chart are drawn at the same charges,
kappa_i = q_i/r_+^2, and L = 3, with mu = r_+^4 (1 + kappa_1)(1 + kappa_2)(1 + kappa_3)/L^2 = 10/3, whose inner
horizon is x^2 + 8x - 8 = 0, x = 2 sqrt 6 - 4, u = 1/x = 1.1124. The hyperbolic chart is drawn at L = 1 and
q = (1, 1/2, 1/4), mu = 11/4, which keeps mu >= q_i so that its gauge fields are real.

The slice t = 0 of the spherical chart is embedded in flat space from the throat r_+ out to the circle where
it lies level, (dR/dl)^2 = 1 with R = r (H_1 H_2 H_3)^{1/6} the radius of its circles, at r = 2.356. The
conformal diagram of the spherical and the continued charts is the tower of Reissner-Nordstrom-anti-de Sitter,
its tortoise coordinate r_*(rho) = int_0^rho s sqrt(P(s))/D(s) ds by quadrature, zero on the singularity and
finite on the conformal boundary. stu_black_holes.md beside this file records each choice.
"""
import functools
import math

import numpy as np

METRIC = "stu_black_holes"
Q = (4.0, 2.0, 1.0)
LL = 3.0
MU = 13.0 / 3.0
MU_FLAT = 10.0 / 3.0
RM = math.sqrt((math.sqrt(321.0) - 17.0) / 2.0)        # the inner horizon of the spherical chart drawn
RHO_P, RHO_M = math.sqrt(1.0 + Q[2]), math.sqrt(RM ** 2 + Q[2])
U_M = 1.0 / (2.0 * math.sqrt(6.0) - 4.0)                 # the inner horizon of the brane chart drawn
SPHERE_PARAMS = {"mu": "13/3", "q_1": "4", "q_2": "2", "q_3": "1", "L": "3"}
FLAT_PARAMS = {"mu": "10/3", "q_1": "4", "q_2": "2", "q_3": "1", "L": "3"}
HYPERBOLIC_PARAMS = {"mu": "11/4", "q_1": "1", "q_2": "1/2", "q_3": "1/4", "L": "1"}
BRANE_PARAMS = {"r_+": "1", "L": "3", "kappa_1": "4", "kappa_2": "2", "kappa_3": "1"}
ANGLES = {"chi": "pi/2", "theta": "pi/2", "phi": "0"}
PLANE = {"x": "0", "y": "0", "z": "0"}
SETTINGS = "$q_1 = 4r_+^2$, $q_2 = 2r_+^2$, $q_3 = r_+^2$, $L = 3r_+$, $\\mu = 13r_+^2/3$, $r_+ = 1$."
FLAT_SETTINGS = "$q_1 = 4r_+^2$, $q_2 = 2r_+^2$, $q_3 = r_+^2$, $L = 3r_+$, $\\mu = 10r_+^2/3$, $r_+ = 1$."
HYPERBOLIC_SETTINGS = "$q_1 = r_+^2$, $q_2 = r_+^2/2$, $q_3 = r_+^2/4$, $L = r_+$, $\\mu = 11r_+^2/4$, $r_+ = 1$."
BRANE_SETTINGS = "$\\kappa_1 = 4$, $\\kappa_2 = 2$, $\\kappa_3 = 1$, $L = 3r_+$, $r_+ = 1$."


def product(r, q=Q):
    """H_1 H_2 H_3 at r."""
    r = np.asarray(r, dtype=float)
    return (1 + q[0] / r ** 2) * (1 + q[1] / r ** 2) * (1 + q[2] / r ** 2)


def areal(r):
    """R = r (H_1 H_2 H_3)^{1/6}, the radius of the circles of the slice t = 0 of the spherical chart drawn."""
    r = np.asarray(r, dtype=float)
    return r * product(r) ** (1 / 6)


# ---------------------------------------------------------------- spacetime diagrams

def register_null(nr):
    nr.DIAGRAMS.extend((
        nr.Diagram(METRIC, "sphere", "radial", "$t$ and $r$", ("t", "r"), (0, 3, -1.5, 1.5), "$r/r_+$", "$ct/r_+$",
                   SPHERE_PARAMS, ANGLES, orient="ingoing", areal=True),
        nr.Diagram(METRIC, "flat", "radial", "$t$ and $r$", ("t", "r"), (0, 3, -1.5, 1.5), "$r/r_+$", "$ct/r_+$",
                   FLAT_PARAMS, PLANE, orient="ingoing"),
        nr.Diagram(METRIC, "hyperbolic", "radial", "$t$ and $r$", ("t", "r"), (0, 3, -1.5, 1.5), "$r/r_+$",
                   "$ct/r_+$", HYPERBOLIC_PARAMS, ANGLES, orient="ingoing"),
        nr.Diagram(METRIC, "continued", "radial", "$t$ and $\\rho$", ("t", "\\rho"), (0, 3, -1.5, 1.5), "$\\rho/r_+$",
                   "$ct/r_+$", SPHERE_PARAMS, ANGLES, orient="ingoing", areal=True),
        nr.Diagram(METRIC, "brane", "radial", "$t$ and $u$", ("t", "u"), (0, 2, -1, 1), "$u$", "$ct/r_+$",
                   BRANE_PARAMS, PLANE, orient="outgoing"),
    ))
    nr.CAPTIONS[(METRIC, "sphere", "radial")] = [
        "The plane of $t$ and $r$ ($\\chi = \\theta = \\pi/2$, $\\phi = 0$) of the black hole with three charges "
        "($q_1 = 4r_+^2$, $q_2 = 2r_+^2$, $q_3 = r_+^2$, $L = 3r_+$), each point in the diagram a 3-sphere, with the "
        "cones closing at the event horizon $r_+$ and at the inner horizon."]
    nr.CAPTIONS[(METRIC, "flat", "radial")] = [
        "The plane of $t$ and $r$ ($x = y = z = 0$) of the black hole with a flat horizon ($q_1 = 4r_+^2$, "
        "$q_2 = 2r_+^2$, $q_3 = r_+^2$, $L = 3r_+$), each point in the diagram a flat 3-space, with the cones closing "
        "at the event horizon $r_+$ and at the inner horizon."]
    nr.CAPTIONS[(METRIC, "hyperbolic", "radial")] = [
        "The plane of $t$ and $r$ ($\\chi = \\theta = \\pi/2$, $\\phi = 0$) of the black hole with a hyperbolic "
        "horizon ($q_1 = r_+^2$, $q_2 = r_+^2/2$, $q_3 = r_+^2/4$, $L = r_+$), each point in the diagram a "
        "hyperbolic 3-space, with the cones closing at the event horizon $r_+$ and at the inner horizon."]
    nr.CAPTIONS[(METRIC, "continued", "radial")] = [
        "The plane of $t$ and $\\rho$ ($\\chi = \\theta = \\pi/2$, $\\phi = 0$) of the black hole with three charges "
        "($q_1 = 4r_+^2$, $q_2 = 2r_+^2$, $q_3 = r_+^2$, $L = 3r_+$), each point in the diagram a 3-sphere, with the "
        "cones closing at the two horizons and at the singularity $\\rho = 0$."]
    nr.CAPTIONS[(METRIC, "brane", "radial")] = [
        "The plane of $t$ and $u = r_+^2/r^2$ ($x = y = z = 0$) of the black brane ($\\kappa_1 = 4$, $\\kappa_2 = 2$, "
        "$\\kappa_3 = 1$, $L = 3r_+$), each point in the diagram a flat 3-space, with the boundary at $u = 0$ and the "
        "cones closing at the horizons $u = 1$ and $u = 1.11$."]


def register_slices(sl):
    sl.FLAT[(METRIC, "sphere", "radial")] = lambda: sl.one(METRIC, lambda m: sl.along(0.0, *m.reach("sphere", "r")))
    sl.FLAT[(METRIC, "continued", "radial")] = lambda: sl.one(
        METRIC, lambda m: sl.along(0.0, *(math.sqrt(r * r + Q[2]) for r in m.reach("sphere", "r"))))
    sl.FLAT_METRICS.add(METRIC)
    sl.HIDDEN[(METRIC, "flat", "radial")] = "another black hole than the one embedded, with a flat horizon"
    sl.HIDDEN[(METRIC, "hyperbolic", "radial")] = "another black hole than the one embedded, with a hyperbolic horizon"
    sl.HIDDEN[(METRIC, "brane", "radial")] = "another black hole than the one embedded, with a flat horizon"


# ---------------------------------------------------------------- embedding diagram

def register_embedding(eb):
    eb.DRAWN[METRIC] = lambda ck, src: embedding(eb, ck, src)
    eb.CAPTIONS[(METRIC, "sphere")] = [
        "The equatorial plane ($\\chi = \\theta = \\pi/2$) at $t = 0$ of the black hole with three charges, each "
        "circle of radius $r(H_1H_2H_3)^{1/6}$, with a throat at the horizon $r = r_+$ joining two exteriors."]


def level_radius(sl):
    """The radius where g_rr = (dR/dr)^2 on the slice, the surface level, beyond the horizon: in eighty digits."""
    import mpmath
    import sympy as sp
    value = sp.lambdify(sl.x, sl.defect, "mpmath")
    with mpmath.workdps(90):
        root = mpmath.findroot(value, (mpmath.mpf(2), mpmath.mpf(3)), solver="anderson", tol=mpmath.mpf(10) ** -80)
        z = sp.Float(mpmath.nstr(root, 85), 85)
    sl.exact[float(z)] = z
    sl.known = {float(z): z}
    return float(z)


def embedding(eb, ck, src):
    """The slice t = 0 of the spherical chart drawn, from the horizon out to the level circle, each
    exterior joined at the bifurcation sphere r = r_+."""
    name = "black hole with three charges"
    sl = eb.Slice(src, METRIC, "sphere", "r", "\\phi", {"t": 0, **eb.EQUATOR, "chi": "pi/2"}, SPHERE_PARAMS)
    rp = sl.horizon_between(0.9, 1.1)
    rm = sl.horizon_between(0.6, 0.75)
    ck.add(f"{name}: the horizons are r_+ = 1 and r_-^2 = (sqrt 321 - 17)/2", abs(rp - 1.0) + abs(rm - RM), 1e-12)
    top = level_radius(sl)
    ck.add(f"{name}: the surface lies level at r = 2.356", abs(top - 2.356) * 1e-3, 1e-5)
    ck.stops(f"{name}, between the horizons", sl, np.linspace(rm, rp, 402)[1:-1])
    ck.stops(f"{name}, beyond the level circle in flat space", sl, np.linspace(top, 40, 402)[1:])
    size = 2 * top
    radii = (1.25, 1.5, 2.0)
    edge = "the circle $r = 2.36\\,r_+$, where the surface lies level"
    near = eb.Piece("exterior", "sheet", sl, rp, top, 0.0, 1,
                    (("throat", "the throat $r = r_+$, the bifurcation sphere, where the other exterior begins"),
                     ("edge", edge)),
                    [(rp, "horizon", "$r = r_+$")] + [(r, "r", None) for r in radii] + [(top, "r", None)], size)
    far = eb.Piece("other_exterior", "sheet2", sl, rp, top, 0.0, -1,
                   (("throat", "the throat $r = r_+$"), ("edge", edge)),
                   [(r, "r2", None) for r in radii] + [(top, "r2", None)], size)
    ck.join(f"{name}, the two exteriors at the throat", near, rp, far, rp)
    for p in (near, far):
        ck.isometry(f"{name}, {p.id}", p)
        ck.radius(f"{name}, {p.id}, rho = R", p, areal, size)
    surface = eb.Surface([near, far])
    fig = eb.figure_of([surface], {"sheet": "cover"}, size)
    eb.ring_label(fig, [0, 0, 0], *near.at(rp), "$r = r_+$", dx=14)
    eb.ring_label(fig, [0, 0, 0], *near.at(1.5), "$1.5\\,r_+$")
    eb.ring_label(fig, [0, 0, 0], *near.at(top), "$2.36\\,r_+$")
    fig.legend("fill", "cover", "the exterior $r > r_+$")
    fig.legend("line", "r", "$r$ constant")
    fig.legend("line", "r2", "the other exterior")
    fig.legend("line", "horizon", "the throat $r = r_+$")
    fig.legend("line", "meridian", "$\\phi$ constant, every $15°$")
    stops = ["Between the two horizons $r$ is a time, and a slice of constant $t$ is a timelike surface."]
    return [eb.view("sphere", "Three Charges", "$r_+$", [surface], fig.done(), settings=SETTINGS, stops=stops)]


# ---------------------------------------------------------------- conformal diagram

def _P(s):
    """P = rho^2 (rho^2 + q_1 - q_3)(rho^2 + q_2 - q_3) of the member drawn."""
    s = np.asarray(s, dtype=float)
    return s * s * (s * s + Q[0] - Q[2]) * (s * s + Q[1] - Q[2])


def _D(s):
    """D = (rho^2 - q_3)^2 - mu (rho^2 - q_3) + P/L^2 of the member drawn."""
    s = np.asarray(s, dtype=float)
    x = s * s - Q[2]
    return x * x - MU * x + _P(s) / LL ** 2


def _D_slope(s):
    s = np.asarray(s, dtype=float)
    p1, p2 = s * s + Q[0] - Q[2], s * s + Q[1] - Q[2]
    return 2 * s * (2 * s * s - 2 * Q[2] - MU) + 2 * s * (p1 * p2 + s * s * p1 + s * s * p2) / LL ** 2


def slope(s):
    """dr_*/d rho = rho sqrt(P)/D, from -g_tt = D/(rho^{4/3}(P_1 P_2)^{2/3}) and g_rho rho = rho^{8/3}(P_1 P_2)^{1/3}/D."""
    s = np.asarray(s, dtype=float)
    return s * np.sqrt(_P(s)) / _D(s)


RESIDUES = [h * math.sqrt(float(_P(h))) / float(_D_slope(h)) for h in (RHO_P, RHO_M)]
JOIN = 4.0                                   # beyond it r_* is integrated from the far side


def _rest(s):
    """rho sqrt(P)/D less its two poles, smooth through both horizons, where it takes its limit."""
    for h in (RHO_P, RHO_M):
        if abs(s - h) < 1e-6:
            return 0.5 * (_rest(h - 2e-6) + _rest(h + 2e-6))
    return float(slope(s)) - sum(A / (s - h) for A, h in zip(RESIDUES, (RHO_P, RHO_M)))


def _far(w):
    """The integrand of r_* beyond JOIN in w = JOIN/rho, finite at w = 0, where it is L^2/JOIN."""
    if w < 1e-12:
        return LL ** 2 / JOIN
    rho = JOIN / w
    return float(slope(rho)) * JOIN / (w * w)


@functools.lru_cache(maxsize=None)
def _tables():
    """The smooth part of r_* on [0, JOIN] and the far part on w in [0, 1], each integrated piece by piece
    in fifteen digits and interpolated by cubic splines."""
    from scipy.integrate import quad
    from scipy.interpolate import CubicSpline
    import warnings
    from scipy.integrate import IntegrationWarning
    tol = dict(epsabs=1e-12, epsrel=1e-12, limit=200)
    near = np.linspace(0.0, JOIN, 4001)
    # Beside a horizon the rest is the difference of two terms near 1/(rho - rho_h), so its last digits
    # are rounding and quad says so; the checks below hold r_* and its slope to what the drawing needs.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", IntegrationWarning)
        steps = [quad(_rest, a, b, **tol)[0] for a, b in zip(near[:-1], near[1:])]
    inner = CubicSpline(near, np.concatenate([[0.0], np.cumsum(steps)]))
    ws = np.linspace(0.0, 1.0, 2001)
    steps = [quad(_far, a, b, **tol)[0] for a, b in zip(ws[:-1], ws[1:])]
    beyond = np.concatenate([[0.0], np.cumsum(steps)])          # int_0^w
    outer = CubicSpline(ws, beyond[-1] - beyond)                 # int_w^1 = int_JOIN^rho
    return inner, outer, float(beyond[-1])


def _poles(rho):
    with np.errstate(divide="ignore"):
        return sum(A * np.log(np.abs(rho / h - 1)) for A, h in zip(RESIDUES, (RHO_P, RHO_M)))


def rstar(rho):
    """r_*(rho) = int_0^rho s sqrt(P)/D ds: its two logarithms at the horizons written out and the rest
    tabulated by quadrature."""
    inner, outer, _ = _tables()
    rho = np.asarray(rho, dtype=float)
    at_join = float(inner(JOIN) + _poles(np.array(JOIN)))
    with np.errstate(divide="ignore", invalid="ignore"):
        inside = inner(np.minimum(rho, JOIN)) + _poles(rho)
        outside = at_join + outer(np.clip(JOIN / np.maximum(rho, JOIN), 0.0, 1.0))
    return np.where(rho <= JOIN, inside, outside)


def far_rstar():
    """The limit of r_* on the conformal boundary, rho -> infinity."""
    inner, _, total = _tables()
    return float(inner(JOIN) + _poles(np.array(JOIN))) + total


def register_conformal(cf):
    cf.DRAWN[METRIC] = lambda ck, src: conformal(cf, ck, src)
    words = "each point in the diagram a 3-sphere, with a timelike singularity at $\\rho = 0$."
    cf.CAPTIONS[(METRIC, "sphere")] = ["The plane of $t$ and $r$ of the black hole with three charges, "
                                       + words[:-1] + ", where $\\rho^2 = r^2 + q_3$."]
    cf.CAPTIONS[(METRIC, "continued")] = ["The plane of $t$ and $\\rho$ of the black hole with three charges, " + words]


def conformal(cf, ck, src):
    HALF, PI = cf.HALF, cf.PI
    name = "black hole with three charges"

    class STUTower(cf.Tower):
        """The tower of the member drawn, in rho, with r_* by quadrature."""

        def __init__(self):
            self.rf = [RHO_P, RHO_M]
            self.Af = list(RESIDUES)
            self.kappa = [1 / (2 * abs(A)) for A in RESIDUES]
            self.kp = self.kappa[0]
            self.far = far_rstar()

        def rstar(self, r):
            return rstar(r)

        def f(self, r):
            return 1 / slope(r)

    T = STUTower()
    fixed = {"chi": "pi/2", **cf.EQUATOR}
    cont = cf.Plane(src, METRIC, "continued", ("t", "\\rho"), fixed, SPHERE_PARAMS)
    sph = cf.Plane(src, METRIC, "sphere", ("t", "r"), fixed, SPHERE_PARAMS)
    rho = np.array([0.05, 0.5, 1.0, 1.3, 1.6, 3.0, 30.0])
    ck.limit(f"{name}: D vanishes at rho_+ = sqrt 2 and rho_- = 1.2075, and slope() is the published "
             "sqrt(g_rho rho/-g_tt)", [float(_D(RHO_P)), float(_D(RHO_M))], [0.0, 0.0], 1e-12)
    published = np.sqrt(np.abs(cont.g_at(0, rho)[2] / cont.g_at(0, rho)[0])) if hasattr(cont, "g_at") else None
    if published is not None:
        ck.limit(f"{name}: dr_*/drho is sqrt(g_rho rho/-g_tt)", np.abs(slope(rho)), published, 1e-10)
    ck.limit(f"{name}: dr_*/drho is rho sqrt(P)/D", (T.rstar(rho + 1e-6) - T.rstar(rho - 1e-6)) / 2e-6, slope(rho), 1e-5)
    ck.limit(f"{name}: r_* vanishes at the singularity rho = 0", [float(T.rstar(0.0))], [0.0], 1e-12)
    ck.limit(f"{name}: r_* is continuous across the join", [float(T.rstar(JOIN - 1e-9))], [float(T.rstar(JOIN + 1e-9))], 1e-7)
    ck.limit(f"{name}: r_* tends to its boundary value", [float(T.rstar(1e7))], [T.far], 1e-5)
    cf.tower_checks(ck, f"{name}, continued", cont, T, 0.05, 4)
    for cell, lo, hi, future in (("I", 1.0 + 1e-3, 40.0, (1, 0)), ("II", RM + 1e-3, 1.0 - 1e-3, (0, -1))):
        ck.chart(f"{name}, the spherical chart in cell {cell}", sph,
                 lambda t, r, c=cell: T.pq(c, t, np.sqrt(np.asarray(r, dtype=float) ** 2 + Q[2])),
                 ck.uniform(-4, 4), ck.uniform(lo, hi), lambda t, r, d=future: d)
    p, q = T.pq("III", np.array([-2.0, 0, 2]), np.full(3, 1e-9))
    ck.limit(f"{name}: rho -> 0 lands on the vertical line X = pi/2", q - p, [HALF] * 3, 1e-6)
    K = cont.kretschmann
    ck.diverges(f"{name}: the Kretschmann scalar diverges at rho = 0", K(0, 1e-2), K(0, 1e-3))
    ck.finite(f"{name}: the Kretschmann scalar is finite at both horizons and at r = 0",
              K(np.zeros(3), np.array([RHO_M, RHO_P, 1.0])))

    D = cf.HyperbolicHoleDrawing(T, T.far)
    box = [-HALF - 0.85, HALF + 0.85, -PI - 0.1, 3 * PI + 0.1]
    times = [c / T.kp for c in (-1.6, -0.6, 0, 0.6, 1.6)]
    rI, rII, rIII = (1.5, 2.0, 3.0, 6.0), (1.25, 1.3, 1.37), (0.3, 0.6, 1.0)
    grids = {"I": (rI, times), "II": (rII, times), "IV": (rII, times), "III": (rIII, times)}
    moment, = cf.slices.moments(METRIC, "sphere")
    lo, hi = moment.reach("sphere", "r")
    line = cf.through_bifurcation(T, ("I'", "I"), math.sqrt(hi * hi + Q[2]), math.sqrt(lo * lo + Q[2]))
    views = []
    for vid, title, cover in (("sphere", "Spherical Horizon", "the exterior, which $t$ and $r$ cover"),
                              ("continued", "Past $r = 0$", "the exterior, which $t$ and $\\rho$ cover")):
        v = cf.View(vid, title, box, vid)
        D.draw(v, grids, cover=[("I", False)])
        v.label_xt([0, HALF + 0.1], "black hole", cls="region")
        v.label_xt([0, 1.5 * PI - 0.1], "white hole", cls="region")
        v.label_xt([0, -HALF], "white hole", cls="region")
        v.label_xt([0, 2.5 * PI], "black hole", cls="region")
        for sx in (1, -1):
            for base in (0, 2 * PI):
                v.label_xt([sx * 0.9, base + 0.3], "exterior", cls="region")
            v.label_xt([sx * HALF, PI + 0.7], "$\\rho = 0$", "l" if sx > 0 else "r", dx=8 * sx)
        v.label_xt([cf.Q4, cf.Q4], "$r_+$", "tl", "small", dx=5, dy=1)
        v.set(fade={"top": 0.9, "bottom": 0.9}, settings=SETTINGS)
        v.legend("cover", cover)
        v.legend("horizon", "the event horizon $r_+$ and the inner horizon")
        v.legend("singular", "the singularity $\\rho = 0$")
        v.legend("boundary", "the conformal boundary")
        v.legend("r", "$\\rho$ constant")
        v.legend("t", "$ct$ constant")
        v.slice(moment, [line])
        views.append(v)
    return views


# ---------------------------------------------------------------- curvature figure

def register_curvature(cv):
    cv.DRAWN[METRIC] = lambda ck, src: curvature(cv, ck, src)


def curvature(cv, ck, src):
    """K and the Ricci scalar against rho from the stored scalars of the continued chart at the member
    drawn: both finite at the horizons and at r = 0, where rho = r_+, both growing without bound toward the
    singularity rho = 0, and far out tending to anti-de Sitter space of five dimensions, K = 40/L^4 and
    R = -20/L^2."""
    params = {"mu": MU, "q_1": Q[0], "q_2": Q[1], "q_3": Q[2], "L": LL}
    names = list(params)
    Kp = cv.Published(src, METRIC, "continued", "kretschmann", ["rho"], names)
    Rp = cv.Published(src, METRIC, "continued", "ricci_scalar", ["rho"], names)
    lo, hi = 0.25, 6.0
    x_axis = cv.axis("$\\rho/r_+$", lo, hi)
    K_scale = {"kind": "log", "from": 1e-1, "to": 1e9}
    R_scale = {"kind": "symlog", "from": -1e2, "to": 1e4, "linear": 1}
    panels = [{"id": "K", "name": "$K\\,r_+^4$", "aspect": 2, "scale": K_scale, "ticks": cv.value_ticks(K_scale)},
              {"id": "R", "name": "$\\mathcal{R}\\,r_+^2$", "aspect": 2, "scale": R_scale,
               "ticks": cv.value_ticks(R_scale)}]
    cv.shape_checks(ck, METRIC, [p["aspect"] for p in panels])
    height = {p["id"]: cv.PLOT_W / p["aspect"] for p in panels}
    curves = (cv.line_curves(ck, METRIC, "K", "K", Kp, "rho", x_axis, K_scale, height["K"], params)
              + cv.line_curves(ck, METRIC, "R", "R", Rp, "rho", x_axis, R_scale, height["R"], params))
    vals = list(params.values())
    rs = np.linspace(lo, hi, cv.SAMPLES)
    name = "black hole with three charges"
    ck(f"{name}: K is positive", bool(np.all(Kp(rs, *vals) > 0)))
    ck(f"{name}: K and R are finite at the horizons and at r = 0",
       bool(np.isfinite(Kp(np.array([RHO_M, RHO_P, 1.0]), *vals)).all()
            and np.isfinite(Rp(np.array([RHO_M, RHO_P, 1.0]), *vals)).all()))
    ck(f"{name}: K grows without bound toward rho = 0",
       float(Kp(np.array([1e-3]), *vals)[0]) > 1e6 * float(Kp(np.array([1e-1]), *vals)[0]))
    far = np.array([1e4])
    ck(f"{name}: far out K tends to 40/L^4", abs(float(Kp(far, *vals)[0]) - 40 / LL ** 4) < 1e-6)
    ck(f"{name}: far out R tends to -20/L^2", abs(float(Rp(far, *vals)[0]) + 20 / LL ** 2) < 1e-6)
    legend = [["line", "K", "$K$"], ["line", "R", "$\\mathcal{R}$"], ["line", "horizon", "the horizons"]]
    return [{
        "id": "radial", "label": "Along $\\rho$", "kind": "line", "system": "continued", "x": x_axis,
        "panels": panels, "curves": curves,
        "marks": [{"class": "horizon", "x": RHO_M}, {"class": "horizon", "x": RHO_P}],
        "legend": legend,
        "caption": ["The Kretschmann scalar $K$ and the Ricci scalar $\\mathcal{R}$ of the black hole against $\\rho$. "
                    "$K$ is finite at the horizons, diverges toward the singularity $\\rho = 0$, and tends to $40/L^4$ "
                    "toward the boundary. "
                    "$\\mathcal{R}$ diverges too and tends to $-20/L^2$ there."],
        "ricci_sentences": [[0, 2]],
        "per_scalar": cv.scalar_clips("A logarithmic scale from $10^{-1}$ to $10^9$.",
                                      "A signed logarithmic scale from $-10^2$ to $10^4$, linear within $1$."),
        "settings": SETTINGS.replace(", $r_+ = 1$", ""),
        "clip": "$K$ logarithmic, $\\mathcal{R}$ signed logarithmic.",
    }]
