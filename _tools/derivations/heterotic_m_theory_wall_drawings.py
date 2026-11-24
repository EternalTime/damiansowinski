"""The drawings of the domain wall of heterotic M-theory, in units of 1/alpha, the inverse of Lukas,
Ovrut, Stelle and Waldram's constant.

The static wall is drawn at alpha = a_0 = b_0 = c_0 = 1 with the walls pi rho = 1.2/alpha apart in y, so
H = (sqrt 2/3)|y| + 1 runs from 1 on the wall at y = 0 to H_1 = 1 + 2 sqrt(2)/5 = 1.566 on the
wall at y = +-pi rho. With b_0 = a_0^6 it is Lehners, McFadden and Turok's conformally flat chart
at h_0 = 1, h = H^{5/2}, through dw/dy = H^{3/2}, and the conformal chart is drawn at that member, its
walls at w = +-w_0 with 2 w_0 = W(pi rho), W(y) = sgn(y)(2/5k)(H^{5/2} - 1) and k = sqrt(2)/3, so that
w_0 = (3/(5 sqrt 2))(H_1^{5/2} - 1) = 0.877. On the plane of t and y the rays keep ct -+ W(y).

The static bulk is drawn at lambda = 1/2, its cosmological horizon at r = 1/lambda = 2, and its rays keep ct -+ r_*, r_* = 3 sqrt(2) int_0^r s^4 ds/(1 - lambda^2 s^2), zero on the
singularity and without bound at the horizon. The moving walls are drawn at the same constants as the
static wall, with A = B = 1 and t_0 = 0, on the branch p = (3/11)(1 + 4/(3 sqrt 3)), q = (2/11)(1 - 2 sqrt 3)
on which the walls grow and the interval shrinks after t_0; their rays keep eta -+ W(y), with the
conformal time eta = (ct)^{1-q}/(1 - q), zero on the singularity.

The slice t = 0 of the static chart is embedded with x_1 rolled up once round, each circle of radius
H^{1/2}. heterotic_m_theory_wall.md beside this file records each choice.
"""
import math

import numpy as np

METRIC = "heterotic_m_theory_wall"
K = math.sqrt(2) / 3                       # the slope of H at alpha = 1
PIRHO = 1.2                                # pi rho, in 1/alpha
H1 = 1 + K * PIRHO                         # H on the wall at y = +-pi rho
W0 = (H1 ** 2.5 - 1) / (5 * K)             # half the conformal width of the interval
LAMBDA = 0.5
P_DOWN = 3 / 11 * (1 + 4 / (3 * math.sqrt(3)))
Q_DOWN = 2 / 11 * (1 - 2 * math.sqrt(3))
FIXED = {"x_1": "0", "x_2": "0", "x_3": "0"}
STATIC_PARAMS = {"alpha": "1", "a_0": "1", "b_0": "1", "c_0": "1", "rho": "6/(5*pi)"}
CONFORMAL_PARAMS = {"alpha": "1", "w_0": "3/(5*sqrt(2))*((1 + 2*sqrt(2)/5)**(5/2) - 1)", "h_0": "1"}
BULK_PARAMS = {"alpha": "1", "lambda": "1/2"}
ROLLING_PARAMS = {**STATIC_PARAMS, "A": "1", "B": "1", "t_0": "0",
                  "p": "3/11*(1 + 4/(3*sqrt(3)))", "q": "2/11*(1 - 2*sqrt(3))"}
E6 = "the wall at $y = 0$"


def H(y):
    return K * np.abs(np.asarray(y, dtype=float)) + 1


def W(y):
    """The conformal distance from the wall at y = 0, int_0^y H^{3/2} dy, odd in y."""
    y = np.asarray(y, dtype=float)
    return np.sign(y) * (H(y) ** 2.5 - 1) / (2.5 * K)


def rstar(r):
    """r_* = 3 sqrt(2) int_0^r s^4/(1 - lambda^2 s^2) ds at lambda = 1/2:
    3 sqrt(2)(32 artanh(r/2) - 16 r - 4 r^3/3)."""
    r = np.asarray(r, dtype=float)
    return 3 * math.sqrt(2) * (32 * np.arctanh(r / 2) - 16 * r - 4 * r ** 3 / 3)


def eta(t):
    """The conformal time c int dt/b = (ct)^{1 - q}/(1 - q) of the moving walls drawn."""
    return np.asarray(t, dtype=float) ** (1 - Q_DOWN) / (1 - Q_DOWN)


# ---------------------------------------------------------------- spacetime diagrams

def register_null(nr):
    nr.DIAGRAMS.extend((
        nr.Diagram(METRIC, "static", "ty", "$t$ and $y$", ("t", "y"), (-PIRHO, PIRHO, -PIRHO, PIRHO),
                   "$\\alpha y$", "$\\alpha ct$", STATIC_PARAMS, FIXED, families=nr.SIDEWAYS, cones=(6, 7),
                   periodic=("y",), lines=(("shell", "r", "0", E6),)),
        nr.Diagram(METRIC, "conformal", "tw", "$t$ and $w$", ("t", "w"), (-W0, W0, -W0, W0),
                   "$\\alpha w$", "$\\alpha ct$", CONFORMAL_PARAMS, FIXED, families=nr.SIDEWAYS, cones=(6, 7),
                   rounded={"w_0": "0.877"}),
        nr.Diagram(METRIC, "static_bulk", "tr", "$t$ and $r$", ("t", "r"), (0, 2, -1.5, 1.5),
                   "$r$", "$\\alpha ct$", BULK_PARAMS, FIXED, cones=(6, 7)),
        nr.Diagram(METRIC, "rolling", "ty", "$t$ and $y$", ("t", "y"), (-PIRHO, PIRHO, 0, 2 * PIRHO),
                   "$\\alpha y$", "$\\alpha ct$", ROLLING_PARAMS, FIXED, families=nr.SIDEWAYS, cones=(6, 7),
                   periodic=("y",), lines=(("shell", "r", "0", E6),), rounded={"p": "0.483", "q": "-0.448"}),
    ))
    nr.CAPTIONS[(METRIC, "static", "ty")] = [
        "The plane of $t$ and $y$ ($x_1 = x_2 = x_3 = 0$), each ray a null geodesic, with one wall at $y = 0$ and "
        "the edges $y = \\pm\\pi\\rho$ one line, the other wall, toward which the cones narrow as $H^{-3/2}$."]
    nr.CAPTIONS[(METRIC, "conformal", "tw")] = [
        "The plane of $t$ and $w$ ($x_1 = x_2 = x_3 = 0$) of the static walls, conformal to the plane of $t$ and "
        "$w$ of Minkowski space, "
        "every ray a null geodesic at 45°, between the wall of negative tension at $w = -w_0$ and the wall of "
        "positive tension at $w = w_0$."]
    nr.CAPTIONS[(METRIC, "static_bulk", "tr")] = [
        "The plane of $t$ and $r$ ($x_1 = x_2 = x_3 = 0$) of the static bulk, each ray a null geodesic, from the "
        "timelike naked singularity $r = 0$ to the cosmological horizon $r = 1/\\lambda$, where the cones close."]
    nr.CAPTIONS[(METRIC, "rolling", "ty")] = [
        "The plane of $t$ and $y$ ($x_1 = x_2 = x_3 = 0$) of the moving walls on the branch that shrinks the "
        "interval, its cones opening as $1/b$ after the singularity $t = t_0$."]
    nr.CLOSED_FORMS[(METRIC, "static", "ty")] = (lambda t, y: t + W(y), lambda t, y: t - W(y), None)
    nr.CLOSED_FORMS[(METRIC, "conformal", "tw")] = (lambda t, w: t + w, lambda t, w: t - w, None)
    nr.CLOSED_FORMS[(METRIC, "static_bulk", "tr")] = (lambda t, r: t + rstar(r), lambda t, r: t - rstar(r),
                                                      lambda t, r: (r > 0.02) & (r < 1.98))
    nr.CLOSED_FORMS[(METRIC, "rolling", "ty")] = (lambda t, y: eta(t) + W(y), lambda t, y: eta(t) - W(y),
                                                  lambda t, y: t > 0.02)


def register_slices(sl):
    sl.FLAT[(METRIC, "static", "ty")] = lambda: sl.one(METRIC, lambda m: sl.along(0.0, *m.reach("static", "y")))
    sl.FLAT[(METRIC, "conformal", "tw")] = lambda: sl.one(METRIC, lambda m: sl.along(0.0, -W0, W0))
    sl.FLAT_METRICS.add(METRIC)
    sl.HIDDEN[(METRIC, "static_bulk", "tr")] = "another spacetime than the one embedded, the bulk in which the walls move"
    sl.HIDDEN[(METRIC, "rolling", "ty")] = "another spacetime than the one embedded, the walls moving"


# ---------------------------------------------------------------- embedding diagram

def register_embedding(eb):
    eb.DRAWN[METRIC] = lambda ck, src: embedding(eb, ck, src)
    eb.CAPTIONS[(METRIC, "static")] = [
        "The surface of $y$ and $x_1$ at $t = 0$ ($x_2 = x_3 = 0$), each circle of radius $a_0H^{1/2}$, "
        "narrowest on the wall at $y = 0$, where the two copies of the interval meet."]


def embedding(eb, ck, src):
    """The slice t = 0 of the static chart drawn, x_2 = x_3 = 0, with x_1 rolled up once round, so that
    each circle has radius H^{1/2}: the circles grow more slowly than the distance across,
    |d rho/dl| = k/(2 H^{5/2}) <= sqrt(2)/6, so the whole slice lies in flat space. The points y and -y are
    one point, so the two halves are two copies of the interval, joined at the wall at y = 0,
    where the slope of the profile changes sign."""
    name = "heterotic M-theory wall"
    fixed_at = {"t": 0, "x_2": 0, "x_3": 0}
    sl = eb.Slice(src, METRIC, "static", "y", "x_1", fixed_at, STATIC_PARAMS)
    size = 4.5
    edge = "the wall at $y = \\pm\\pi\\rho$, where this end and the other are one circle"
    near = eb.Piece("near", "sheet", sl, -PIRHO, 0.0, 0.0, 1,
                    (("edge", edge), ("crease", "the wall at $y = 0$, where the two copies meet")),
                    [(-PIRHO, "surface", None), (-0.6, "r", None), (0.0, "surface", None)], size)
    near.z = near.z - float(near.z[-1])
    far = eb.Piece("far", "sheet2", sl, 0.0, PIRHO, 0.0, 1,
                   (("crease", "the wall at $y = 0$"), ("edge", edge)),
                   [(0.6, "r2", None), (PIRHO, "surface", None)], size)
    ck.add(f"{name}, the two copies meet at the wall: one point",
           float(np.max(np.abs(np.array(near.at(0.0)) - far.at(0.0)))), eb.JOIN)
    for piece in (near, far):
        ck.isometry(f"{name}, {piece.id}", piece)
        ck.radius(f"{name}, {piece.id}: rho = H^(1/2)", piece, lambda y: H(y) ** 0.5, size)
    surface = eb.Surface([near, far])
    fig = eb.figure_of([surface], {"sheet": "cover"}, size)
    eb.ring_label(fig, [0, 0, 0], *near.at(0.0), "$y = 0$")
    eb.ring_label(fig, [0, 0, 0], *far.at(PIRHO), "$y = \\pi\\rho$")
    eb.ring_label(fig, [0, 0, 0], *near.at(-PIRHO), "$-\\pi\\rho$")
    fig.legend("fill", "cover", "the copy $y < 0$")
    fig.legend("line", "surface", "the walls")
    fig.legend("line", "r", "$y = -\\pi\\rho/2$")
    fig.legend("line", "r2", "$y = \\pi\\rho/2$")
    fig.legend("line", "meridian", "$x_1$ constant")
    return [eb.view("static", "Static Walls", "$1/\\alpha$", [surface], fig.done(),
                    settings="$\\alpha = a_0 = b_0 = c_0 = 1$, $\\pi\\rho = 1.2/\\alpha$, $x_1$ rolled up with "
                             "period $2\\pi/\\alpha$.")]


# ---------------------------------------------------------------- conformal diagrams

UNIT = "$a_0 = b_0 = c_0 = 1$ and $\\pi\\rho = 1.2/\\alpha$."


def register_conformal(cf):
    cf.DRAWN[METRIC] = lambda ck, src: conformal(cf, ck, src)
    cf.CAPTIONS[(METRIC, "static")] = [
        "The plane of $t$ and $y$ ($x_1 = x_2 = x_3 = 0$), each point a flat 3-space, a strip from $i^-$ to $i^+$ "
        "with one wall down its middle and the other on both edges."]
    cf.CAPTIONS[(METRIC, "conformal")] = [
        "The plane of $t$ and $w$ ($x_1 = x_2 = x_3 = 0$), each point a flat 3-space, one copy of the interval from "
        "$i^-$ to $i^+$ between its two walls."]
    cf.CAPTIONS[(METRIC, "static_bulk")] = [
        "The plane of $t$ and $r$ ($x_1 = x_2 = x_3 = 0$), each point a flat 3-space, with the naked singularity "
        "$r = 0$ on the left and the cosmological horizon on the two right edges."]
    cf.CAPTIONS[(METRIC, "rolling")] = [
        "The plane of $t$ and $y$ ($x_1 = x_2 = x_3 = 0$) of the moving walls, each point a flat 3-space, the upper "
        "half of the static strip with the singularity $t = t_0$ on its lower edge."]


def _strip(cf, lo, hi, S):
    """The polygon in X and T between the lines of constant x = lo and x = hi under Minkowski's maps,
    over the times S."""
    right = np.array(cf.xt(*cf.mink_pq(S, np.full_like(S, hi)))).T
    left = np.array(cf.xt(*cf.mink_pq(S, np.full_like(S, lo)))).T
    return [list(map(float, P)) for P in right] + [list(map(float, P)) for P in left[::-1]]


def conformal(cf, ck, src):
    """Four views, in units of 1/alpha, each by Minkowski's maps p, q = arctan(T -+ X) of a pair in which
    the plane is conformally flat. On the plane of t and y of the static walls the metric is
    a_0^2 H(-c^2dt^2 + dW^2), W(y) = int_0^y (b_0/a_0)H^{3/2} dy, so the drawing is the strip |W| <= W(pi rho)
    from i^- to i^+, its edges the one wall at y = +-pi rho and its axis the wall at y = 0
    field, and the conformal chart is the copy w = W - w_0 of one side. The moving walls have
    a_0^2 H b^2(-d eta^2 + dW^2), eta = c int dt/b = (ct)^{1-q}/(1 - q) at A = B = 1 and t_0 = 0, finite at
    the singularity, which is the line eta = 0, T = 0, so they fill the upper half of the static strip. The
    static bulk has r^2(1 - lambda^2 r^2)(-c^2dt^2 + dr_*^2), r_* = 3 sqrt(2) int_0^r s^4 ds/(1 - lambda^2 s^2)
    from 0 on the singularity to infinity at the cosmological horizon, so it fills the right half of the
    diamond; along a ray the affine parameter is int sqrt(-g_tt g_rr) dr = 3 sqrt(2) int r^6 dr, which is
    checked to stay finite, 384 sqrt(2)/7 from the singularity to the horizon, so its edges are drawn as
    the horizon and no infinity. Every map is checked against the published metric, the Kretschmann scalar
    to diverge on the two singularities and to stay finite on the walls and the horizon."""
    PI = cf.PI
    st = cf.Plane(src, METRIC, "static", ("t", "y"), FIXED, STATIC_PARAMS)
    co = cf.Plane(src, METRIC, "conformal", ("t", "w"), FIXED, CONFORMAL_PARAMS)
    bu = cf.Plane(src, METRIC, "static_bulk", ("t", "r"), FIXED, BULK_PARAMS)
    ro = cf.Plane(src, METRIC, "rolling", ("t", "y"), FIXED, ROLLING_PARAMS)

    def static(t, y):
        return cf.mink_pq(t, W(y))

    def conformal_flat(t, w):
        return cf.mink_pq(t, w)

    def bulk(t, r):
        return cf.mink_pq(t, rstar(r))

    def rolling(t, y):
        return cf.mink_pq(eta(t), W(y))

    future = lambda t, x: (1, 0)  # noqa: E731
    for side, sign in (("y < 0", -1), ("y > 0", 1)):
        ck.chart(f"heterotic M-theory wall, static, {side}", st, static, ck.uniform(-4, 4),
                 sign * ck.uniform(0.001, PIRHO - 0.001), future)
        ck.chart(f"heterotic M-theory wall, rolling, {side}", ro, rolling, ck.uniform(0.01, 4),
                 sign * ck.uniform(0.001, PIRHO - 0.001), future)
    ck.chart("heterotic M-theory wall, conformal", co, conformal_flat, ck.uniform(-4, 4),
             ck.uniform(-W0 + 0.001, W0 - 0.001), future)
    ck.chart("heterotic M-theory wall, static bulk", bu, bulk, ck.uniform(-4, 4), ck.uniform(0.01, 1.99), future)
    ck.limit("heterotic M-theory wall: the conformal width of the interval is 2 w_0", float(W(PIRHO)), 2 * W0, 1e-12)
    times = ck.uniform(-3, 3, 40)
    ck.finite("heterotic M-theory wall: the curvature is finite on the walls of the static chart",
              np.concatenate([st.kretschmann(times, np.full(40, y)) for y in (-PIRHO, -1e-6, 1e-6, PIRHO)]))
    ck.finite("heterotic M-theory wall: the curvature is finite on the walls of the conformal chart",
              np.concatenate([co.kretschmann(times, np.full(40, w)) for w in (-W0, W0)]))
    ck.diverges("heterotic M-theory wall: the static bulk is singular at r = 0",
                bu.kretschmann(times, np.full(40, 1e-1)), bu.kretschmann(times, np.full(40, 1e-2)))
    ck.finite("heterotic M-theory wall: the curvature is finite at the cosmological horizon",
              bu.kretschmann(times, np.full(40, 1.999)))
    r = np.linspace(1e-6, 2 - 1e-9, 400001)
    g = bu.metric(np.zeros_like(r), r)
    root = np.sqrt(-g[0] * g[2])
    affine = float(np.sum((root[1:] + root[:-1]) / 2 * np.diff(r)))
    ck.limit("heterotic M-theory wall: a ray reaches the horizon at the affine parameter 384 sqrt(2)/7",
             affine, 384 * math.sqrt(2) / 7, 1e-6)
    ys = ck.uniform(-PIRHO + 0.01, PIRHO - 0.01, 40)
    ck.diverges("heterotic M-theory wall: the moving walls are singular at t = t_0",
                ro.kretschmann(np.full(40, 1e-2), ys), ro.kretschmann(np.full(40, 1e-3), ys))
    ck.limit("heterotic M-theory wall: the singularity t = t_0 lies on T = 0",
             cf.xt(*rolling(np.full(40, 1e-12), ys))[1], np.zeros(40), 1e-6)

    S = cf.spread(-np.inf, np.inf, 600, 10)
    up = cf.spread(0.0, np.inf, 600, 10)
    edge = W(PIRHO)
    labels = (((0, PI), "$i^+$", "b", -6), ((0, -PI), "$i^-$", "t", 6))

    def ends(v, which=labels):
        for at, text, anchor, dy in which:
            v.layers.append({"kind": "point", "class": "infinity", "at": [round(at[0], 4), round(at[1], 4)]})
            v.label_xt(at, text, anchor, dy=dy)

    views = []
    v = cf.View("static", "Static Walls", [-2.6, 2.6, -PI - 0.25, PI + 0.25], "static")
    strip = _strip(cf, -edge, edge, S)
    v.fill("region", strip)
    v.fill("cover", strip)
    cf.grid(v, "r", lambda c, s: static(s, c), (-0.6, 0.6), S)
    for t0 in (-2, -1, 0, 1, 2):
        span = np.linspace(-PIRHO, PIRHO, 801)
        v.curve("t", *static(np.full_like(span, t0), span))
    for y in (-PIRHO, 0.0, PIRHO):
        v.curve("surface", *static(S, np.full_like(S, y)))
    ends(v)
    v.label_xt([0, 0.35], "$y = 0$", "l", "small", dx=5)
    v.legend("cover", "two copies of the interval")
    v.legend("surface", "the walls")
    v.legend("r", "$y = \\pm\\pi\\rho/2$")
    v.legend("t", "$ct$ constant")
    views.append(v)

    v = cf.View("conformal", "Conformally Flat", [-2.0, 2.0, -PI - 0.25, PI + 0.25], "conformal")
    strip = _strip(cf, -W0, W0, S)
    v.fill("region", strip)
    v.fill("cover", strip)
    cf.grid(v, "r", lambda c, s: conformal_flat(s, c), (-W0 / 2, 0.0, W0 / 2), S)
    for t0 in (-2, -1, 0, 1, 2):
        span = np.linspace(-W0, W0, 401)
        v.curve("t", *conformal_flat(np.full_like(span, t0), span))
    for w in (-W0, W0):
        v.curve("surface", *conformal_flat(S, np.full_like(S, w)))
    ends(v)
    v.legend("cover", "one copy of the interval")
    v.legend("surface", "the walls, $w = \\pm w_0$")
    v.legend("r", "$w$ constant")
    v.legend("t", "$ct$ constant")
    views.append(v)

    v = cf.View("static_bulk", "Static Bulk", [-0.3, PI + 0.45, -PI - 0.25, PI + 0.25], "static_bulk")
    triangle = [[0, -PI], [PI, 0], [0, PI]]
    v.fill("region", triangle)
    v.fill("cover", triangle)
    cf.grid(v, "r", lambda c, s: bulk(s, c), (0.5, 1.0, 1.5), S)
    for t0 in (-2, -1, 0, 1, 2):
        span = np.linspace(1e-4, 2 - 1e-9, 2001)
        v.curve("t", *bulk(np.full_like(span, t0), span))
    v.line("singular", [[[0, -PI], [0, PI]]])
    v.line("horizon", [[[0, PI], [PI, 0]], [[PI, 0], [0, -PI]]])
    v.label_xt([cf.HALF, cf.HALF], "$r = 1/\\lambda$", "bl", dx=5, dy=-3)
    v.legend("cover", "the static region")
    v.legend("singular", "$r = 0$")
    v.legend("horizon", "$r = 1/\\lambda$, at a finite affine parameter")
    v.legend("r", "$r$ constant")
    v.legend("t", "$ct$ constant")
    v.set(settings="$\\lambda = 1/2$.")
    views.append(v)

    v = cf.View("rolling", "Rolling Walls", [-2.6, 2.6, -0.25, PI + 0.25], "rolling")
    half = _strip(cf, -edge, edge, up)
    v.fill("region", half)
    v.fill("cover", half)
    cf.grid(v, "r", lambda c, s: rolling(s, c), (-0.6, 0.6), up)
    for t0 in (0.5, 1, 2):
        span = np.linspace(-PIRHO, PIRHO, 801)
        v.curve("t", *rolling(np.full_like(span, t0), span))
    for y in (-PIRHO, 0.0, PIRHO):
        v.curve("surface", *rolling(up, np.full_like(up, y)))
    span = np.linspace(-PIRHO, PIRHO, 801)
    v.curve("singular", *rolling(np.full_like(span, 1e-12), span))
    ends(v, labels[:1])
    v.legend("cover", "two copies of the interval")
    v.legend("singular", "$t = t_0$")
    v.legend("surface", "the walls")
    v.legend("r", "$y = \\pm\\pi\\rho/2$")
    v.legend("t", "$t$ constant")
    v.set(settings="$a_0 = b_0 = c_0 = A = B = 1$, $t_0 = 0$, $\\pi\\rho = 1.2/\\alpha$, $p \\approx 0.483$, $q \\approx -0.448$.")
    views.append(v)

    moment, = cf.slices.moments(METRIC, "static")
    lo, hi = moment.reach("static", "y")
    y = np.linspace(lo, hi, 601)
    views[0].set(settings=UNIT)
    views[0].slice(moment, [static(np.zeros_like(y), y)])
    w = np.linspace(-W0, W0, 401)
    views[1].set(settings="$h_0 = 1$ and $w_0 \\approx 0.877/\\alpha$.")
    views[1].slice(moment, [conformal_flat(np.zeros_like(w), w)])
    return views


# ---------------------------------------------------------------- curvature figure

def register_curvature(cv):
    cv.DRAWN[METRIC] = lambda ck, src: curvature(cv, ck, src)


def curvature(cv, ck, src):
    """K and the Ricci scalar across the interval from the stored scalars of the conformal chart at the
    member drawn, h_0 = 1 and w_0 = 0.877/alpha, which is the static chart at a_0 = b_0 = c_0 = 1 and
    pi rho = 1.2/alpha with h = H^{5/2}: K = 106 alpha^4/81 h^{24/5} and R = 14 alpha^2/9 h^{12/5}, both
    finite on the walls. The conformal chart's scalars hold no delta, so both read as the bulk's; the
    static chart's stored K is checked to agree on both walls."""
    params = {"alpha": 1.0, "w_0": W0, "h_0": 1.0}
    names = list(params)
    vals = list(params.values())
    Kp = cv.Published(src, METRIC, "conformal", "kretschmann", ["w"], names)
    Rp = cv.Published(src, METRIC, "conformal", "ricci_scalar", ["w"], names)
    Ks = cv.Published(src, METRIC, "static", "kretschmann", ["y"], ["alpha", "a_0", "b_0", "c_0"])
    x_axis = cv.axis("$\\alpha w$", -W0, W0)
    K_scale = {"kind": "log", "from": 1e-3, "to": 1e1}
    R_scale = {"kind": "log", "from": 1e-2, "to": 1e1}
    panels = [{"id": "K", "name": "$K/\\alpha^4$", "aspect": 2, "scale": K_scale, "ticks": cv.value_ticks(K_scale)},
              {"id": "R", "name": "$R/\\alpha^2$", "aspect": 2, "scale": R_scale,
               "ticks": cv.value_ticks(R_scale)}]
    cv.shape_checks(ck, METRIC, [p["aspect"] for p in panels])
    height = {p["id"]: cv.PLOT_W / p["aspect"] for p in panels}
    curves = (cv.line_curves(ck, METRIC, "K", "K", Kp, "w", x_axis, K_scale, height["K"], params)
              + cv.line_curves(ck, METRIC, "R", "R", Rp, "w", x_axis, R_scale, height["R"], params))
    ws = np.linspace(-W0, W0, cv.SAMPLES)
    name = "heterotic M-theory wall"
    K = Kp(ws, *vals)
    ck(f"{name}: K is positive and falls across the interval", bool(np.all(K > 0) and np.all(np.diff(K) < 0)))
    ck(f"{name}: K = 106/81 and R = 14/9 on the wall of negative tension, where h = 1",
       abs(float(Kp(np.array([-W0]), *vals)[0]) - 106 / 81) < 1e-12
       and abs(float(Rp(np.array([-W0]), *vals)[0]) - 14 / 9) < 1e-12)
    ck(f"{name}: the static chart's K agrees on both walls",
       abs(float(Ks(np.array([1e-12]), 1, 1, 1, 1)[0]) - float(Kp(np.array([-W0]), *vals)[0])) < 1e-9
       and abs(float(Ks(np.array([PIRHO]), 1, 1, 1, 1)[0]) - float(Kp(np.array([W0]), *vals)[0])) < 1e-9)
    legend = [["line", "K", "$K$"], ["line", "R", "$R$"]]
    return [{
        "id": "across", "label": "Across the Interval", "kind": "line", "system": "conformal", "x": x_axis,
        "panels": panels, "curves": curves, "marks": [], "legend": legend,
        "caption": ["The Kretschmann scalar $K$ and the Ricci scalar $R$ of the static walls across the "
                    "interval, from the wall of negative tension at $w = -w_0$ to the wall of positive tension at "
                    "$w = w_0$. $K$ falls as $h^{-24/5}$, finite on both walls. "
                    "$R$ falls as $h^{-12/5}$."],
        "ricci_sentences": [[0, 2]],
        "per_scalar": cv.scalar_clips("A logarithmic scale from $10^{-3}$ to $10$.",
                                      "A logarithmic scale from $10^{-2}$ to $10$."),
        "settings": "$h_0 = 1$ and $w_0 \\approx 0.877/\\alpha$.",
        "clip": "$K$ and $R$ logarithmic.",
    }]
