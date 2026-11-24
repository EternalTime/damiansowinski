"""The drawings of Horowitz and Ross's naked black holes, each chart at one member, in units of its horizon.

The dilaton chart is drawn at a = 1/2 and r_- = r_+/2, so the power of R is n = a^2/(1 + a^2) = 1/5 and
R = r(1 - r_-/r)^{1/5}. Its radial light rays obey c dt = +-(R^2/((r - r_+)(r - r_-))) dr, and the tortoise
coordinate

    r_* = A ln|(r - r_+)/(r_+ - r_-)| + int_{r_-}^r ((R^2/(s - r_-) - A)/(s - r_+)) ds,   A = R(r_+)^2/(r_+ - r_-),

has no closed form at a = 1/2: the integrand is smooth at r_+ and has the integrable power (s - r_-)^{2n - 1}
at r_-, so it is summed by quadrature, and r_* vanishes at the singularity r = r_-. Near the horizon
r_* = A ln|r - r_+| + finite, so the surface gravity is kappa = 1/(2A) and Kruskal's U = -+exp(-kappa u),
V = +-exp(kappa v) draw the plane as Schwarzschild's square with the singularity r_- on UV = 1.

The chart with two charges is drawn at r_- = 3r_+/5 and Sigma = 3r_+/10, where 1/f = (r^2 - Sigma^2)/((r - r_+)
(r - r_-)) is 1 plus two simple poles, so its conformal diagram is the tower of Reissner and Nordstrom with
r_* shifted by a constant to vanish at the singularity r = |Sigma|, which puts that singularity on the
vertical lines X = +-pi/2, and kappa_-/kappa_+ = (r_+^2 - Sigma^2)/(r_-^2 - Sigma^2) = 91/27.

The string chart is drawn at r_1 = r_p = r_0/2, where its light rays obey c dt = +-((r + r_1)/(r - r_0)) dr, so
r_* = r + (r_0 + r_1) ln|r/r_0 - 1|, which vanishes at the singularity r = 0, with kappa = 1/(2(r_0 + r_1)).

The embedding diagram draws the equatorial plane of the dilaton chart's slice t = 0, whose circles have the
radius R, out from the throat at the horizon r_+, where the second exterior begins. The curvature figure is
the stored Kretschmann and Ricci scalars of the dilaton chart against r, out from the singularity r_-.
naked_black_holes.md beside this file records each choice.
"""
import functools
import math

import numpy as np

METRIC = "naked_black_holes"
DIL, TWO, STR = "dilaton", "two_charges", "string_ns"
EQUATOR = {"theta": "pi/2", "phi": "0"}
DIL_PARAMS = {"r_+": "1", "r_-": "1/2", "a": "1/2"}
TWO_PARAMS = {"r_+": "1", "r_-": "3/5", "Sigma": "3/10"}
STR_PARAMS = {"r_0": "1", "r_1": "1/2", "r_p": "1/2"}
RP, RM, N = 1.0, 0.5, 0.2                 # the dilaton chart drawn: r_+, r_-, a^2/(1 + a^2) at a = 1/2
SIG = 0.3                                  # |Sigma| of the chart with two charges
R0, R1 = 1.0, 0.5                          # the string chart drawn: r_0 and r_1 = r_p
TOP = 6.0                                  # how far out the slice t = 0 is embedded
DIL_SETTINGS = "$a = 1/2$, $r_- = r_+/2$, $r_+ = 1$."
TWO_SETTINGS = "$r_- = 0.6\\,r_+$, $\\Sigma = 0.3\\,r_+$, $r_+ = 1$."
STR_SETTINGS = "$r_1 = r_p = r_0/2$, $r_0 = 1$."


def areal(r):
    """R = r(1 - r_-/r)^{1/5} of the dilaton chart drawn."""
    r = np.asarray(r, dtype=float)
    return r * (1 - RM / r) ** N


def _g(s):
    """R^2/(s - r_-) = s^{2 - 2n}(s - r_-)^{2n - 1}."""
    return s ** (2 - 2 * N) * (s - RM) ** (2 * N - 1)


A_DIL = _g(RP)                             # R(r_+)^2/(r_+ - r_-)
KAPPA_DIL = 1 / (2 * A_DIL)
_SLOPE_AT_HORIZON = A_DIL * ((2 - 2 * N) / RP + (2 * N - 1) / (RP - RM))


def _h(s):
    d = s - RP
    if abs(d) < 1e-7:
        return _SLOPE_AT_HORIZON
    return (_g(s) - A_DIL) / d


K_SUB = 1 / (2 * N)                        # s - r_- = w^k takes the power (s - r_-)^{2n - 1} out of the integrand
TOL = dict(epsabs=1e-14, epsrel=1e-13, limit=400)


def _hw(w):
    """h ds in the variable w of s - r_- = w^k: k(s^{2 - 2n} - A w^{k - 1})/(s - r_+), smooth at w = 0."""
    s = RM + w ** K_SUB
    if abs(s - RP) < 1e-7:
        return _SLOPE_AT_HORIZON * K_SUB * w ** (K_SUB - 1)
    return K_SUB * (s ** (2 - 2 * N) - A_DIL * w ** (K_SUB - 1)) / (s - RP)


def _near(r):
    """The integral of the smooth part of dr_*/dr from the singularity out to r <= r_+, in w."""
    from scipy.integrate import quad
    return quad(_hw, 0.0, (r - RM) ** (1 / K_SUB), **TOL)[0]


@functools.lru_cache(maxsize=1)
def _inside():
    """The same integral from the singularity to the horizon."""
    return _near(RP)


def _rstar_one(r):
    from scipy.integrate import quad
    if r <= RM:
        return 0.0
    log = A_DIL * math.log(abs((r - RP) / (RP - RM))) if r != RP else -math.inf
    if r <= RP:
        return log + _near(r)
    if r <= 10.0:
        return log + _inside() + quad(_h, RP, r, **TOL)[0]
    middle = quad(_h, RP, 10.0, **TOL)[0]
    # Far out h tends to 1, so the rest is r - 10 and the integral of h - 1 in ln s, with
    # R^2/(s - r_-) - s = s((1 - r_-/s)^{2n - 1} - 1) written without the cancellation.
    def rest(x):
        s = math.exp(x)
        return s * (s * math.expm1((2 * N - 1) * math.log1p(-RM / s)) - A_DIL + RP) / (s - RP)
    tail = quad(rest, math.log(10.0), math.log(r), **TOL)[0]
    return log + _inside() + middle + (r - 10.0) + tail


_rstar_vec = np.vectorize(_rstar_one, otypes=[float])


def rstar_dilaton(r):
    """The tortoise coordinate of the dilaton chart drawn, zero at the singularity r = r_-."""
    return _rstar_vec(np.asarray(r, dtype=float))


def rstar_string(r):
    """The tortoise coordinate of the string chart drawn, r + (r_0 + r_1) ln|r/r_0 - 1|, zero at r = 0."""
    r = np.asarray(r, dtype=float)
    with np.errstate(divide="ignore"):
        return r + (R0 + R1) * np.log(np.abs(r / R0 - 1))


# ---------------------------------------------------------------- spacetime diagrams

def register_null(nr):
    nr.DIAGRAMS.extend((
        nr.Diagram(METRIC, DIL, "radial", "$t$ and $r$", ("t", "r"), (0.5, 6.5, -3, 3), "$r/r_+$", "$ct/r_+$",
                   DIL_PARAMS, EQUATOR, orient="ingoing", areal=True),
        nr.Diagram(METRIC, TWO, "radial", "$t$ and $r$", ("t", "r"), (0.3, 6.3, -3, 3), "$r/r_+$", "$ct/r_+$",
                   TWO_PARAMS, EQUATOR, orient="ingoing", areal=True),
        nr.Diagram(METRIC, STR, "radial", "$t$ and $r$", ("t", "r"), (0, 6, -3, 3), "$r/r_0$", "$ct/r_0$",
                   STR_PARAMS, EQUATOR, orient="ingoing", areal=True),
    ))
    nr.CAPTIONS[(METRIC, DIL, "radial")] = [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) of the black hole of coupling $a = 1/2$ "
        "($r_- = r_+/2$), each point in the diagram a 2-sphere of area $4\\pi R^2$, with the cones closing at the "
        "horizon $r = r_+$ and the singularity $r = r_-$."]
    nr.CAPTIONS[(METRIC, TWO, "radial")] = [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) of the black hole with two charges "
        "($r_- = 0.6\\,r_+$, $\\Sigma = 0.3\\,r_+$), each point in the diagram a 2-sphere of area "
        "$4\\pi(r^2 - \\Sigma^2)$, with the cones closing at the horizons $r_+$ and $r_-$ and the singularity "
        "$r = |\\Sigma|$."]
    nr.CAPTIONS[(METRIC, STR, "radial")] = [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) in the string metric ($r_1 = r_p = r_0/2$), "
        "each point in the diagram a 2-sphere of area $4\\pi r^2$, with the cones closing at the horizon "
        "$r = r_0$ and the singularity $r = 0$."]


def register_slices(sl):
    sl.FLAT[(METRIC, DIL, "radial")] = lambda: sl.one(METRIC, lambda m: sl.along(0.0, *m.reach(DIL, "r")))
    sl.FLAT_METRICS.add(METRIC)
    sl.HIDDEN[(METRIC, TWO, "radial")] = "another black hole than the one embedded, of coupling $a = 1/2$"
    sl.HIDDEN[(METRIC, STR, "radial")] = ("the string metric, which measures the manifold another way than the "
                                          "Einstein metric embedded")


# ---------------------------------------------------------------- embedding diagram

def register_embedding(eb):
    eb.DRAWN[METRIC] = lambda ck, src: embedding(eb, ck, src)
    eb.CAPTIONS[(METRIC, DIL)] = [
        "The equatorial plane ($\\theta = \\pi/2$) at $t = 0$ of the black hole of coupling $a = 1/2$ "
        "($r_- = r_+/2$), each circle of radius $R$, with a throat at the horizon $r = r_+$ joining two exteriors."]


def embedding(eb, ck, src):
    """The slice t = 0 of the dilaton chart drawn, from the horizon out to r = 6 r_+, each exterior joined
    at the bifurcation sphere r = r_+, where the circles' radius R is smallest."""
    name = "naked black hole of coupling a = 1/2"
    sl = eb.Slice(src, METRIC, DIL, "r", "\\phi", {"t": 0, **eb.EQUATOR}, DIL_PARAMS)
    ck.add(f"{name}: the horizon is at r_+", abs(sl.horizons()[0] - RP), 1e-12)
    ck.stops(f"{name}, inside the horizon", sl, np.linspace(RM, RP, 402)[1:-1])
    size = 2 * TOP
    radii = (2.0, 3.0, 4.0, 5.0)
    near = eb.Piece("exterior", "sheet", sl, RP, TOP, 0.0, 1,
                    (("throat", "the throat $r = r_+$, the bifurcation sphere, where the other exterior begins"),
                     ("edge", "the surface runs on as $r \\to \\infty$, toward a plane")),
                    [(RP, "horizon", "$r = r_+$")] + [(r, "r", None) for r in radii] + [(TOP, "r", None)], size)
    far = eb.Piece("other_exterior", "sheet2", sl, RP, TOP, 0.0, -1,
                   (("throat", "the throat $r = r_+$"), ("edge", "the other exterior, as $r \\to \\infty$")),
                   [(r, "r2", None) for r in radii] + [(TOP, "r2", None)], size)
    ck.join(f"{name}, the two exteriors at the throat", near, RP, far, RP)
    for p in (near, far):
        ck.isometry(f"{name}, {p.id}", p)
        ck.radius(f"{name}, {p.id}, rho = R", p, areal, size)
    surface = eb.Surface([near, far])
    fig = eb.figure_of([surface], {"sheet": "cover"}, size)
    eb.ring_label(fig, [0, 0, 0], *near.at(RP), "$r = r_+$", dx=14)
    eb.ring_label(fig, [0, 0, 0], *near.at(3.0), "$3\\,r_+$")
    eb.ring_label(fig, [0, 0, 0], *near.at(TOP), "$6\\,r_+$")
    fig.legend("fill", "cover", "the exterior $r > r_+$")
    fig.legend("line", "r", "$r$ constant, $2$ to $6\\,r_+$")
    fig.legend("line", "r2", "the other exterior")
    fig.legend("line", "horizon", "the throat $r = r_+$")
    fig.legend("line", "meridian", "$\\phi$ constant, every $15°$")
    stops = ["Inside the horizon, $r_- < r < r_+$, $r$ is a time, and a slice of constant $t$ is a timelike "
             "surface."]
    return [eb.view(DIL, "Coupling $a = 1/2$", "$r_+$", [surface], fig.done(), settings=DIL_SETTINGS, stops=stops)]


# ---------------------------------------------------------------- conformal diagrams

def kruskal(cf, cell, t, rs, kappa):
    """(p, q) of points (t, r_*) of a cell: U = -+exp(-kappa u), V = +-exp(kappa v)."""
    t = np.asarray(t, dtype=float)
    u, v = t - rs, t + rs
    p, q = cf.atan_exp(-kappa * u), cf.atan_exp(kappa * v)
    return {"I": (-p, q), "II": (p, q), "I'": (p, -q), "IV": (-p, -q)}[cell]


def register_conformal(cf):
    cf.DRAWN[METRIC] = lambda ck, src: conformal(cf, ck, src)
    cf.CAPTIONS[(METRIC, DIL)] = [
        "The plane of $t$ and $r$ of the black hole of coupling $a = 1/2$ ($r_- = r_+/2$), each point in the "
        "diagram a 2-sphere, with Schwarzschild's causal structure and a spacelike singularity at $r = r_-$."]
    cf.CAPTIONS[(METRIC, TWO)] = [
        "The plane of $t$ and $r$ of the black hole with two charges ($r_- = 0.6\\,r_+$, $\\Sigma = 0.3\\,r_+$), "
        "each point in the diagram a 2-sphere, with a tower of horizons and a timelike singularity at "
        "$r = |\\Sigma|$."]
    cf.CAPTIONS[(METRIC, STR)] = [
        "The plane of $t$ and $r$ of the string metric ($r_1 = r_p = r_0/2$), each point in the diagram a "
        "2-sphere, with Schwarzschild's causal structure and a spacelike singularity at $r = 0$."]


def _square(cf, ck, src, vid, system, params, rstar, kappa, horizon, singular, names, settings, moment):
    """A plane with one horizon and a spacelike singularity where r_* = 0, as Kruskal's square."""
    HALF, PI = cf.HALF, cf.PI
    title, h_text, s_text, unit = names
    plane = cf.Plane(src, METRIC, system, ("t", "r"), cf.EQUATOR, params)
    far = 30.0
    for cell, rs, future, where in (("I", (horizon + 1e-3, far), (1, 0), "exterior"),
                                    ("II", (singular + 1e-2, horizon - 1e-3), (0, -1), "black hole"),
                                    ("I'", (horizon + 1e-3, far), (-1, 0), "other exterior"),
                                    ("IV", (singular + 1e-2, horizon - 1e-3), (0, 1), "white hole")):
        ck.chart(f"{title}, {where}", plane, lambda t, r, c=cell: kruskal(cf, c, t, rstar(r), kappa),
                 ck.uniform(-15, 15), ck.uniform(*rs), lambda t, r, d=future: d)
    p, q = kruskal(cf, "II", np.array([-3.0, 0.0, 3.0]), rstar(np.full(3, singular + 1e-12)), kappa)
    ck.limit(f"{title}: the singularity lands on T = pi/2", p + q, [HALF] * 3, 1e-5)
    p, q = kruskal(cf, "I", np.array([3.0]), rstar(np.array([horizon + 1e-12])), kappa)
    ck.limit(f"{title}: the horizon at fixed t lands on the bifurcation sphere", cf.point(p[0], q[0]), [0, 0], 1e-3)
    p, q = kruskal(cf, "I", np.array([0.0]), rstar(np.array([1e8])), kappa)
    ck.limit(f"{title}: r -> infinity at t = 0 lands on i0", cf.point(p[0], q[0]), [PI, 0], 1e-3)
    K = plane.kretschmann
    ck.diverges(f"{title}: the Kretschmann scalar diverges at the singularity", K(0, singular + 1e-4),
                K(0, singular + 1e-5))
    ck.finite(f"{title}: the Kretschmann scalar is finite at the horizon",
              K(np.zeros(3), np.array([horizon - 1e-3, horizon, horizon + 1e-3])))
    box = [-PI - 0.25, PI + 0.25, -HALF - 0.25, HALF + 0.25]
    outline = [[PI, 0], [HALF, HALF], [-HALF, HALF], [-PI, 0], [-HALF, -HALF], [HALF, -HALF]]
    exterior = [[0, 0], [HALF, -HALF], [PI, 0], [HALF, HALF]]
    v = cf.View(vid, title, box, system)
    v.fill("region", outline)
    v.fill("cover", exterior)
    t = cf.S_ALL
    for r in (1.25, 1.5, 2.0, 3.0, 6.0):
        v.curve("r", *kruskal(cf, "I", t, rstar(np.full_like(t, r * horizon)), kappa))
    rr = cf.spread(horizon, np.inf, 600, 14)
    for tt in (-4, -2, -1, 0, 1, 2, 4):
        v.curve("t", *kruskal(cf, "I", np.full_like(rr, tt * horizon), rstar(rr), kappa))
    v.line("horizon", [[[-HALF, -HALF], [HALF, HALF]], [[HALF, -HALF], [-HALF, HALF]]])
    v.line("singular", [[[-HALF, HALF], [HALF, HALF]], [[-HALF, -HALF], [HALF, -HALF]]], zig=True)
    v.line("scri", [[[PI, 0], [HALF, HALF]], [[PI, 0], [HALF, -HALF]],
                    [[-PI, 0], [-HALF, HALF]], [[-PI, 0], [-HALF, -HALF]]])
    v.label_xt([PI, 0], "$i^0$", "l", dx=6)
    v.label_xt([-PI, 0], "$i^0$", "r", dx=-6)
    v.label_xt([0, HALF], s_text, "b", dy=-8)
    v.label_xt([0, -HALF], s_text, "t", dy=8)
    v.label_xt([-HALF / 2, HALF / 2], h_text, "tr", "small", dx=-6, dy=2)
    v.legend("cover", f"the exterior $r > {unit}$")
    v.legend("r", "$r$ constant")
    v.legend("t", "$t$ constant")
    v.legend("scri", "null infinity, $r \\to \\infty$")
    v.legend("horizon", f"the horizon {h_text}")
    v.legend("singular", f"{s_text}, a curvature singularity")
    if moment is not None:
        lo, hi = moment.reach(system, "r")
        r_m = np.linspace(lo + 1e-9, hi, 400)
        left = kruskal(cf, "I'", 0 * r_m, rstar(r_m[::-1]), kappa)
        right = kruskal(cf, "I", 0 * r_m, rstar(r_m), kappa)
        v.slice(moment, [(np.concatenate([left[0], right[0]]), np.concatenate([left[1], right[1]]))])
    v.set(settings=settings)
    return v


def conformal(cf, ck, src):
    import sympy as sp

    views = []
    moment = cf.slices.moments(METRIC, DIL)[0]
    views.append(_square(cf, ck, src, DIL, DIL, DIL_PARAMS, rstar_dilaton, KAPPA_DIL, RP, RM,
                         ("Coupling $a = 1/2$", "$r = r_+$", "$r = r_-$", "r_+"), DIL_SETTINGS, moment))
    # Near the horizon the integrand of r_* is A/(r - r_+), so r_* - A ln|r - r_+| stays finite there.
    ck.limit("naked black hole a = 1/2: r_* - A ln|r - r_+| is finite at the horizon",
             [float(rstar_dilaton(RP + 1e-9) - A_DIL * math.log(1e-9))],
             [float(rstar_dilaton(RP + 1e-7) - A_DIL * math.log(1e-7))], 1e-5)
    views.append(_square(cf, ck, src, STR, STR, STR_PARAMS, rstar_string, 1 / (2 * (R0 + R1)), R0, 0.0,
                         ("String metric", "$r = r_0$", "$r = 0$", "r_0"), STR_SETTINGS, None))

    # The chart with two charges: the tower of Reissner and Nordstrom, r_* shifted to vanish at r = |Sigma|.
    pl = cf.Plane(src, METRIC, TWO, ("t", "r"), cf.EQUATOR, TWO_PARAMS)
    roots = sorted(sp.solve(sp.numer(sp.together(pl.gi[1, 1])), pl.x1), reverse=True)
    base = cf.Tower(-pl.g[0, 0], pl.x1, roots)

    class Shifted(type(base)):
        def __init__(self):
            self.__dict__.update(base.__dict__)
            self.shift = float(type(base).rstar(base, np.array(SIG)))

        def rstar(self, r):
            return type(base).rstar(self, r) - self.shift

    T = Shifted()
    rp, rm = T.rf
    name = "naked black hole with two charges"
    cf.tower_checks(ck, name, pl, T, SIG + 1e-3, 15)
    p, q = T.pq("III", np.array([-6.0, 0, 6]), np.full(3, SIG + 1e-12))
    ck.limit(f"{name}: r -> |Sigma| lands on the vertical line X = pi/2", q - p, [cf.HALF] * 3, 1e-9)
    ck.limit(f"{name}: the horizons are the roots of the published g^rr", [rp, rm], [1.0, 0.6], 1e-12)
    ck.diverges(f"{name}: the Kretschmann scalar diverges at r = |Sigma|",
                pl.kretschmann(0, SIG + 1e-3), pl.kretschmann(0, SIG + 1e-4))
    ck.finite(f"{name}: the Kretschmann scalar is finite at both horizons",
              pl.kretschmann(np.zeros(2), np.array([rp, rm])))
    D = cf.TowerDrawing(T, True, SIG)
    PI = cf.PI
    box = [-PI - 0.45, PI + 0.45, -PI - 0.1, 3 * PI + 0.1]
    times = [c / T.kp for c in (-1.6, -0.6, 0, 0.6, 1.6)]
    grids = {"I": ([1.5, 2.0, 3.0, 6.0], times), "II": ([0.7, 0.8, 0.9], times), "IV": ([0.7, 0.8, 0.9], times),
             "III": ([0.35, 0.45, 0.55], times)}
    v = cf.View(TWO, "Two charges", box, TWO)
    D.draw(v, grids, cover=[("I", False), ("II", False), ("III", False)])
    D.labels(v)
    for layer in v.layers:
        if layer.get("text") == "$r = 0$":
            layer["text"] = "$r = |\\Sigma|$"
    v.set(fade={"top": 0.9, "bottom": 0.9})
    v.legend("cover", "the regions the chart covers")
    v.legend("r", "$r$ constant")
    v.legend("t", "$t$ constant")
    v.legend("horizon", "the horizons $r_+$ and $r_-$")
    v.legend("singular", "the singularity $r = |\\Sigma|$")
    v.legend("scri", "null infinity $\\mathscr{I}^\\pm$")
    v.set(settings=TWO_SETTINGS)
    views.append(v)
    return views


# ---------------------------------------------------------------- curvature figure

def register_curvature(cv):
    cv.DRAWN[METRIC] = lambda ck, src: curvature(cv, ck, src)


def curvature(cv, ck, src):
    """K and the Ricci scalar against r from the stored scalars of the dilaton chart at a = 1/2 and
    r_- = r_+/2. The Ricci scalar is held to 2(grad phi)^2, the trace of Einstein's equations with the
    dilaton e^{-2 phi} = (1 - r_-/r)^{2a/(1 + a^2)} of HR98, the Maxwell stress being traceless: it is
    negative inside the horizon, zero on it and positive outside; K is positive, finite at the horizon
    and grows without bound toward the singularity r_-."""
    params = {"r_+": RP, "r_-": RM, "a": 0.5}
    names = list(params)
    Kp = cv.Published(src, METRIC, DIL, "kretschmann", ["r"], names)
    Rp = cv.Published(src, METRIC, DIL, "ricci_scalar", ["r"], names)
    lo = 0.55
    x_axis = cv.axis("$r/r_+$", lo, TOP)
    K_scale = {"kind": "log", "from": 1e-4, "to": 1e4}
    R_scale = {"kind": "symlog", "from": -1e2, "to": 1, "linear": 1e-3}
    panels = [{"id": "K", "name": "$K\\,r_+^4$", "aspect": 2, "scale": K_scale, "ticks": cv.value_ticks(K_scale)},
              {"id": "R", "name": "$\\mathcal{R}\\,r_+^2$", "aspect": 2, "scale": R_scale,
               "ticks": cv.value_ticks(R_scale)}]
    cv.shape_checks(ck, METRIC, [p["aspect"] for p in panels])
    height = {p["id"]: cv.PLOT_W / p["aspect"] for p in panels}
    curves = (cv.line_curves(ck, METRIC, "K", "K", Kp, "r", x_axis, K_scale, height["K"], params)
              + cv.line_curves(ck, METRIC, "R", "R", Rp, "r", x_axis, R_scale, height["R"], params))
    vals = list(params.values())
    rs = np.linspace(lo, TOP, cv.SAMPLES)
    a = 0.5
    name = "naked black hole at a = 1/2"
    dphi = -(a / (1 + a * a)) * RM / (rs * (rs - RM))
    lapse = (rs - RP) * (rs - RM) / areal(rs) ** 2
    ck(f"{name}: the Ricci scalar is 2(grad phi)^2 = 2 F phi'^2",
       bool(np.allclose(Rp(rs, *vals), 2 * lapse * dphi ** 2, rtol=1e-11)))
    ck(f"{name}: K is positive", bool(np.all(Kp(rs, *vals) > 0)))
    ck(f"{name}: the Ricci scalar vanishes at the horizon", abs(float(Rp.exact([RP] + vals))) < 1e-12)
    ck(f"{name}: K is finite at the horizon", bool(np.isfinite(Kp(np.array([RP]), *vals)).all()))
    ck(f"{name}: K grows without bound toward r_-",
       float(Kp(np.array([RM + 1e-6]), *vals)[0]) > 1e3 * float(Kp(np.array([RM + 1e-2]), *vals)[0]))
    inside = np.linspace(lo, RP - 1e-3, 50)
    ck(f"{name}: the Ricci scalar is negative inside the horizon", bool(np.all(Rp(inside, *vals) < 0)))
    legend = [["line", "K", "$K$"], ["line", "R", "$\\mathcal{R}$"], ["line", "horizon", "$r = r_+$"]]
    return [{
        "id": "radial", "label": "Along $r$", "kind": "line", "system": DIL, "x": x_axis,
        "panels": panels, "curves": curves, "marks": [{"class": "horizon", "x": RP}],
        "legend": legend,
        "caption": ["The Kretschmann scalar $K$ and the Ricci scalar $\\mathcal{R}$ of the black hole against $r$. "
                    "$K$ is finite at the horizon and diverges toward the singularity $r = r_-$. "
                    "$\\mathcal{R}$ is twice the squared gradient of the dilaton, zero at the horizon, negative "
                    "inside it."],
        "ricci_sentences": [[0, 2]],
        "per_scalar": cv.scalar_clips("A logarithmic scale from $10^{-4}$ to $10^4$.",
                                      "A signed logarithmic scale from $-10^2$ to $1$, linear within $10^{-3}$."),
        "settings": DIL_SETTINGS,
        "clip": "$K$ logarithmic, $\\mathcal{R}$ signed logarithmic.",
    }]
