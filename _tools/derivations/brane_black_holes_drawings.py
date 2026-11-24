"""The drawings of the black holes on a wall of Emparan, Horowitz and Myers, in units of the radius of
anti-de Sitter space, ell = 1/A = 1.

Every drawing of the flat wall is of one hole, mu = 1/8, on the wall x = 0. There the plane of t and y
has the metric -(1 + 2 mu y) dt^2 + dy^2/(y^4 (1 + 2 mu y)), which in r = -ell/y and t-hat = ell t is
Schwarzschild's plane with r_s = 2 mu ell = 1/4, the first paper's \\bhbrane: the horizon y = -1/(2 mu) = -4,
r = r_s, and the wall's infinity y -> 0. The rays keep t -+ y_* with y_* = -1/y + 2 mu ln|(1 + 2 mu y)/y|,
which is r_*/ell, and the conformal diagram is Kruskal and Szekeres's of that plane, cells I and II, the
exterior and the black hole the chart covers. The moment t = 0 of the wall has the metric
dy^2/(y^4 (1 + 2 mu y)) + dvarphi^2/y^2, varphi of period Delta varphi = 2 pi s, s = 1/(x_2 (1 + 3 mu x_2))
with x_2 the positive root of 1 - x^2 - 2 mu x^3; swept as varphi = s phi it is a surface of revolution
with circles of radius s r, Aryal, Ford and Vilenkin's equator at b = s, which near the horizon rises
as Flamm's paraboloid does and far out opens into a cone of half angle arcsin s.

The two-wall chart is drawn at A = 1, m = 1, lambda = 1/2 and k = 1, where G = (1 - x)(1 + x + 2x^2) has
x_2 = 1 and H = -(y + 1/2)(2y^2 - 2y + 1) the one root y_0 = -1/2, on the axis x = x_2, the fixed set
of the rotation, between the horizon and the second wall y = 0. The rotating chart is drawn on a flat
wall, lambda = 0, k = -1, at mA = 1/8 and a = 3/25, on the wall x = 0 with varphi divided out, where
the metric orthogonal to the circles is (H dt^2 - dy^2/H)/(A^2 y^2) with H = -y^2 (1 + 2 mA y + a^2 y^2),
whose roots y = -25/4 and y = -100/9 are Kerr's r_+ = 4/25 and r_- = 9/100, and g_tt vanishes on the
ergosurface y = -1/(2 mA) = -4, Kerr's r = 2m.
"""
import math

import numpy as np
import sympy as sp

METRIC = "brane_black_holes"
MU = 1 / 8
RS = 2 * MU                      # r_s = 2 mu ell on the wall
HORIZON = -1 / (2 * MU)          # y = -4
FLAT = {"ell": 1, "mu": "1/8"}
TWO = {"A": 1, "m": 1, "lambda": "1/2", "k": 1}
ROT = {"A": 1, "m": "1/8", "lambda": 0, "k": -1, "a": "3/25"}
ROT_HORIZONS = (-25 / 4, -100 / 9)
_X = sp.Symbol("x", positive=True)
X2 = float(sp.nsolve(1 - _X ** 2 - 2 * sp.Rational(1, 8) * _X ** 3, _X, 0.9, prec=40))
FOLD = 1 / (X2 * (1 + 3 * MU * X2))          # Delta varphi / 2 pi
DEFICIT = 360 * (1 - FOLD)                    # the deficit angle, in degrees
REACH = 3.0                                  # the embedding's outer radius r, in units of ell
SETTINGS = f"$\\ell = 1$ and $\\mu = 1/8$, so that $\\Delta\\varphi = {360 * FOLD:.0f}°$."
TWO_SETTINGS = "$A = 1$, $m = 1$, $\\lambda = 1/2$, and $k = 1$, so that $x_2 = 1$ and $y_0 = -1/2$."
ROT_SETTINGS = "$A = 1$, $m = 1/8$, $\\lambda = 0$, $k = -1$, and $a = 3/25$."


def ystar(y, mu=MU):
    """y_* = -1/y + 2 mu ln|(1 + 2 mu y)/y|, with dy_*/dy = 1/(y^2 (1 + 2 mu y)), which is r_*/ell."""
    y = np.asarray(y, dtype=float)
    with np.errstate(divide="ignore"):
        return -1 / y + 2 * mu * np.log(np.abs((1 + 2 * mu * y) / y))


def ystar_axis(y):
    """The tortoise coordinate of the two-wall chart's axis, int dy/F with F = (y + 1/2)(2y^2 - 2y + 1):
    (2/5) ln|y + 1/2| - (1/5) ln(2y^2 - 2y + 1) + (4/5) arctan(2y - 1)."""
    y = np.asarray(y, dtype=float)
    with np.errstate(divide="ignore"):
        return 0.4 * np.log(np.abs(y + 0.5)) - 0.2 * np.log(2 * y * y - 2 * y + 1) + 0.8 * np.arctan(2 * y - 1)


def _rotating_ystar():
    """The tortoise coordinate of the rotating chart's wall, int dy/F with F = y^2 (1 + y/4 + 9 y^2/625)
    = (9/625) y^2 (y + 25/4)(y + 100/9), by partial fractions, each logarithm of an absolute value."""
    y = sp.Symbol("y")
    F = y ** 2 * (1 + y / 4 + sp.Rational(9, 625) * y ** 2)
    integral = sp.integrate(sp.apart(1 / F, y), y)
    assert sp.simplify(sp.diff(integral, y) - 1 / F) == 0
    return sp.lambdify(y, integral.replace(sp.log, lambda arg: sp.log(sp.Abs(arg))), "numpy")


ystar_rotating = _rotating_ystar()


# ---------------------------------------------------------------- spacetime diagrams

def register_null(nr):
    nr.DIAGRAMS.append(nr.Diagram(
        METRIC, "flat_wall", "wall", "the wall, $x = 0$", ("t", "y"), (-6, 0, -3, 3), "$y$", "$t$",
        FLAT, {"x": "0", "varphi": "0"}, orient="ingoing",
        # y -> 0 is the wall's infinity, r = -ell/y -> infinity, where g^yy = y^4 (1 + 2 mu y)/ell^2 vanishes
        # with the conformal factor's divergence, and no horizon.
        where="-y", where_is_infinity=True))
    nr.DIAGRAMS.append(nr.Diagram(
        METRIC, "two_walls", "axis", "the axis, $x = x_2$", ("t", "y"), (-1.5, 0, -0.75, 0.75), "$y$", "$t$",
        TWO, {"x": "1", "varphi": "0"}, orient="ingoing",
        lines=(("surface", "r", "0", "the second wall $y = 0$"),)))
    nr.DIAGRAMS.append(nr.Diagram(
        METRIC, "rotating", "wall", "the wall, $x = 0$", ("t", "y"), (-14, 0, -7, 7), "$y$", "$t$",
        ROT, {"x": "0"}, orient="ingoing", quotient="varphi", mark_gtt="the ergosurface",
        cone="future cone of no angular momentum", where="-y", where_is_infinity=True))
    nr.CAPTIONS[(METRIC, "flat_wall", "wall")] = [
        "The plane of $t$ and $y$ on the wall ($x = 0$, $\\varphi = 0$) of the hole with $\\mu = 1/8$, which in "
        "$r = -\\ell/y$ is Schwarzschild's plane with $r_s = 2\\mu\\ell$: the cones close at the horizon $y = -4$, "
        "narrow toward the wall's infinity $y \\to 0$, and inside the horizon point to the singularity "
        "$y \\to -\\infty$."]
    nr.CAPTIONS[(METRIC, "two_walls", "axis")] = [
        "The plane of $t$ and $y$ on the axis ($x = x_2 = 1$, $\\varphi = 0$), which the rotation about it holds "
        "fixed, so its null curves are null geodesics: the cones close at the horizon $y_0 = -1/2$, and light "
        "reaches the second wall $y = 0$ in a finite $t$."]
    nr.CAPTIONS[(METRIC, "rotating", "wall")] = [
        "The wall ($x = 0$) in $t$ and $y$ with $\\varphi$ divided out, the equatorial plane of Kerr's black hole "
        "of mass $m$ and rotation parameter $a/A$: the cones close at the horizons $y = -25/4$ and $-100/9$, "
        "with the ergosurface $y = -4$ dotted."]
    nr.CLOSED_FORMS[(METRIC, "flat_wall", "wall")] = (
        lambda t, y: t + ystar(y), lambda t, y: t - ystar(y),
        lambda t, y: (np.abs(y - HORIZON) > 0.02) & (y < -0.02))
    nr.CLOSED_FORMS[(METRIC, "two_walls", "axis")] = (
        lambda t, y: t + ystar_axis(y), lambda t, y: t - ystar_axis(y),
        lambda t, y: np.abs(y + 0.5) > 0.01)
    nr.CLOSED_FORMS[(METRIC, "rotating", "wall")] = (
        lambda t, y: t + ystar_rotating(y), lambda t, y: t - ystar_rotating(y),
        lambda t, y: (np.abs(y - ROT_HORIZONS[0]) > 0.05) & (np.abs(y - ROT_HORIZONS[1]) > 0.05) & (y < -0.05))


HIDDEN = {(METRIC, "two_walls", "axis"): "another member of the family than the hole on a flat wall whose moment "
                                         "t = 0 is embedded, between two walls curved like anti-de Sitter space",
          (METRIC, "rotating", "wall"): "another member of the family than the hole whose moment t = 0 is "
                                        "embedded, a spinning hole"}


def slice_table(slices):
    """The moment t = 0 of the embedding view, along y on the flat wall's plane of t and y."""
    def wall():
        return slices.one(METRIC, lambda m: slices.along(0.0, *m.reach("flat_wall", "y")), view_id="wall")
    return {(METRIC, "flat_wall", "wall"): wall}


# ---------------------------------------------------------------- embedding diagram

def register_embedding(em):
    em.DRAWN[METRIC] = lambda ck, src: embedding_views(ck, src, em)
    em.CAPTIONS[(METRIC, "wall")] = [
        f"The wall ($x = 0$) at $t = 0$ from the horizon $r = 2\\mu\\ell$ to $r = 3\\ell$, steep at the horizon "
        f"and opening far out into a cone of deficit angle $2\\pi - \\Delta\\varphi = {DEFICIT:.0f}°$."]


def embedding_views(ck, src, em):
    """The moment t = 0 of the flat wall's chart on the wall x = 0, swept as varphi = s phi with
    s = Delta varphi/2 pi, from the horizon y = -4 to y = -1/3, r = 3 ell. Its circles have the radius
    s ell/|y| = s r and dz/dy = sqrt(g_yy - (drho/dy)^2), real everywhere outside the horizon, since in r
    it is sqrt(1/(1 - r_s/r) - s^2) with s < 1. It is checked an isometry of the published metric, its
    circles s r, and its height the closed form of Aryal, Ford and Vilenkin's equator at b = s with
    r_s = 1/4, z = r_s (w sqrt(1 + k^2 w^2) + arsinh(k w)/k) with w = sqrt(r/r_s - 1) and k = sqrt(1 - s^2),
    which tends to a cone of slope sqrt(1 - s^2)/s far out."""
    name = "black hole on a flat wall"
    sl = em.Slice(src, METRIC, "flat_wall", "y", "\\varphi", {"t": 0, "x": 0}, FLAT,
                  swept={"varphi": f"{FOLD!r}*varphi"})
    lo, hi = HORIZON, -1 / REACH
    size = 2 * FOLD * REACH
    k = math.sqrt(1 - FOLD * FOLD)
    radii = (0.5, 1.0, 1.5, 2.0, 2.5)
    marks = ([(lo, "horizon", "$r = 2\\mu\\ell$")] + [(-1 / r, "r", None) for r in radii]
             + [(hi, "r", "$r = 3\\ell$")])
    wall = em.Piece("wall", "sheet", sl, lo, hi, 0.0, 1,
                    (("stops", "the horizon $y = -4$, where the slice meets the black hole"),
                     ("edge", "the wall runs on, opening into its cone, to its infinity $y \\to 0$")),
                    marks, size)
    ck.isometry(f"{name}, the wall at t = 0", wall)
    ck.radius(f"{name}, the circles of radius s r = s ell/|y|", wall, lambda y: FOLD / np.abs(y), size)

    def closed(y):
        w = np.sqrt(np.maximum(-1 / (np.asarray(y, dtype=float) * RS) - 1, 0))
        return RS * (w * np.sqrt(1 + k * k * w * w) + np.arcsinh(k * w) / k)
    ck.form(f"{name}, the height of Aryal, Ford and Vilenkin's equator at b = s", wall, closed, size)
    far = sl.slope(-1e-6, "-")
    ck.add(f"{name}: far out the wall is a cone, dz/drho = sqrt(1 - s^2)/s",
           abs(far[1] / far[0] - k / FOLD), 1e-4)
    ck.add(f"{name}: the period of varphi is 2 pi/(x_2 (1 + 3 mu x_2)), x_2 a root of 1 - x^2 - 2 mu x^3",
           abs(1 - X2 ** 2 - 2 * MU * X2 ** 3), 1e-15)
    surface = em.Surface([wall])
    fig = em.figure_of([surface], {"sheet": "cover"}, size)
    em.ring_label(fig, [0, 0, 0], *wall.at(lo), "$r = 2\\mu\\ell$", dx=14)
    em.ring_label(fig, [0, 0, 0], *wall.at(-1.0), "$\\ell$")
    em.ring_label(fig, [0, 0, 0], *wall.at(hi), "$3\\ell$")
    fig.legend("fill", "cover", "the wall at $t = 0$")
    fig.legend("line", "r", "$r$ constant, every $\\ell/2$")
    fig.legend("line", "horizon", "the horizon $r = 2\\mu\\ell$")
    fig.legend("line", "meridian", "$\\varphi$ constant, every $\\Delta\\varphi/24$")
    return [em.view("wall", "The wall", "$\\ell$", [surface], fig.done(), settings=SETTINGS)]


# ---------------------------------------------------------------- conformal diagram

def register_conformal(cf):
    cf.DRAWN[METRIC] = lambda ck, src: conformal_views(ck, src, cf)
    cf.CAPTIONS[(METRIC, "flat_wall")] = [
        "The wall's plane of $t$ and $y$ ($x = 0$), Schwarzschild's in $r = -\\ell/y$, in Kruskal and Szekeres's "
        "coordinates: the exterior and the black hole, whose singularity lies on the wall."]


def conformal_views(ck, src, cf):
    """The flat wall's plane of t and y at x = 0 in cells I and II of Schwarzschild's tower with
    r_s = 2 mu ell = 1/4, through r = -ell/y and t-hat = ell t, checked against the published metric at
    random points of both cells; the singularity y -> -infinity lands on T = pi/2, the wall's infinity at
    t = 0 on i0, and the horizon at fixed t on the bifurcation circle, and the published Kretschmann
    scalar diverges at y -> -infinity and is finite on the horizon."""
    name = "black hole on a flat wall"
    pl = cf.Plane(src, METRIC, "flat_wall", ("t", "y"), {"x": "0", "varphi": "0"}, FLAT)
    r = sp.Symbol("r", positive=True)
    T = cf.Tower(1 - sp.Rational(1, 4) / r, r, [sp.Rational(1, 4)])

    def outside(t, y):
        return T.pq("I", t, -1 / np.asarray(y, dtype=float))

    def inside(t, y):
        return T.pq("II", t, -1 / np.asarray(y, dtype=float))
    rng = np.random.default_rng(1999)
    # The surface gravity is 1/(2 r_s) = 2/ell, so the arctangents flatten within a few ell of the horizon:
    # the samples keep to |t| < 4 and r between 0.01 ell and 10 ell.
    ck.chart(f"{name}, exterior", pl, outside, rng.uniform(-4, 4, 4000), -rng.uniform(0.1, 3.999, 4000),
             lambda t, y: (1, 0))
    ck.chart(f"{name}, black hole", pl, inside, rng.uniform(-4, 4, 4000), -10.0 ** rng.uniform(0.603, 2, 4000),
             lambda t, y: (0, -1))
    p, q = inside(np.array([-5.0, 0, 5]), np.full(3, -1e9))
    ck.limit(f"{name}: y -> -infinity in the black hole lands on T = pi/2", p + q, [cf.HALF] * 3)
    p, q = outside(np.array([0.0]), np.array([-1e-9]))
    ck.limit(f"{name}: y -> 0 at t = 0 lands on i0, (X, T) = (pi, 0)", cf.point(p[0], q[0]), [cf.PI, 0], 1e-3)
    p, q = outside(np.array([0.5]), np.array([HORIZON + 1e-12]))
    ck.limit(f"{name}: y -> -4 at fixed t lands on the bifurcation circle", cf.point(p[0], q[0]), [0, 0], 1e-4)
    K = pl.kretschmann
    ck.diverges(f"{name}: the Kretschmann scalar diverges at y -> -infinity", K(0, -1e2), K(0, -1e3))
    ck.finite(f"{name}: the Kretschmann scalar is finite on the horizon y = -4",
              K(np.zeros(3), np.array([HORIZON - 1e-3, HORIZON, HORIZON + 1e-3])))

    half, pi, q4 = cf.HALF, cf.PI, cf.Q4
    box = [-half - 0.3, pi + 0.3, -half - 0.3, half + 0.3]
    region = [[0, 0], [half, -half], [pi, 0], [half, half], [-half, half]]
    v = cf.View("flat_wall", "Flat Wall", box, "flat_wall")
    v.fill("region", region)
    v.fill("cover", region)
    t = cf.spread(-np.inf, np.inf, 500, 9)
    for r_out in (0.3, 0.5, 1.0, 2.0):
        v.curve("r", *outside(t, np.full_like(t, -1 / r_out)))
    for r_in in (0.05, 0.125, 0.2):
        v.curve("r", *inside(t, np.full_like(t, -1 / r_in)))
    rr = RS + (cf.spread(0, np.inf, 500, 14))
    for tt in (-4, -2, -1, 0, 1, 2, 4):
        v.curve("t", *outside(np.full_like(rr, tt), -1 / rr))
    v.line("scri", [[[pi, 0], [half, half]], [[pi, 0], [half, -half]]])
    v.line("horizon", [[[0, 0], [half, half]], [[0, 0], [half, -half]], [[0, 0], [-half, half]]])
    v.line("singular", [[[-half, half], [half, half]]], zig=True)
    for at in ((pi, 0), (half, half), (half, -half)):
        v.layers.append({"kind": "point", "class": "infinity", "at": cf.rounded(at)})
    v.label_xt([pi, 0], "$i^0$", "l", dx=6)
    v.label_xt([half, half], "$i^+$", "b", dy=-6)
    v.label_xt([half, -half], "$i^-$", "t", dy=6)
    v.label_xt([3 * q4, q4], "$\\mathscr{I}^+$", "bl", dx=4, dy=-4)
    v.label_xt([3 * q4, -q4], "$\\mathscr{I}^-$", "tl", dx=4, dy=4)
    v.label_xt([0, half], "$y \\to -\\infty$", "b", dy=-8)
    v.label_xt([half, -0.95], "exterior", cls="region")
    v.label_xt([0, 1.15], "black hole", cls="region")
    cf.label_on(v, outside(0.0, -1 / 0.5), "$r = \\ell/2$")
    v.legend("cover", "the region $t$ and $y$ cover")
    v.legend("r", "$y$ constant")
    v.legend("t", "$t$ constant")
    v.legend("horizon", "the horizon $y = -4$")
    v.legend("singular", "the singularity $y \\to -\\infty$")
    v.legend("scri", "the wall's null infinity, $y \\to 0$")
    v.set(settings=SETTINGS)
    moment = cf.slices.moments(METRIC)[0]
    lo, hi = moment.reach("flat_wall", "y")
    ys = np.linspace(lo, hi, 200)
    v.slice(moment, [outside(0 * ys, ys)])
    return [v]


# ---------------------------------------------------------------- curvature figure

def register_curvature(cv):
    cv.DRAWN[METRIC] = lambda ck, src: curvature_views(ck, src, cv)


def curvature_views(ck, src, cv):
    """The stored Kretschmann scalar of the flat wall's chart, read from the metric file and never
    recomputed, on the wall x = 0 against r = -ell/y, from r = 0.1 ell through the horizon r_s = ell/4 to
    r = 4 ell, beside Schwarzschild's published scalar at the same r_s: on the wall K = 24/ell^4 +
    12 r_s^2/r^6, Schwarzschild's plus the constant of anti-de Sitter space, which is checked at every
    sampled point, with K falling outward to within a part in 10^3 of 24/ell^4 at r = 4 ell. The stored
    Ricci scalar is -12/ell^2 everywhere."""
    lo, hi = 0.1, 4.0
    x_axis = cv.axis("$r/\\ell$", lo, hi)
    scale = {"kind": "log", "from": 1e-3, "to": 1e6}
    panels = [{"id": "K", "name": "$K\\,\\ell^4$", "aspect": 2, "scale": scale, "ticks": cv.value_ticks(scale)}]
    cv.shape_checks(ck, METRIC, [p["aspect"] for p in panels])
    height = cv.PLOT_W / panels[0]["aspect"]
    params = {"ell": 1.0, "mu": MU}
    stored = cv.nr_published(src, METRIC, "flat_wall", "kretschmann", ["y", "x"], ["ell", "mu"])
    K = cv.Held(stored, {"x": 0.0})
    Ric = cv.Held(cv.nr_published(src, METRIC, "flat_wall", "ricci_scalar", ["y", "x"], ["ell", "mu"]), {"x": 0.0})
    for y in (-10.0, HORIZON, -1.0, -0.25):
        ck(f"{METRIC} flat_wall R = -12/ell^2 at y = {y}", abs(float(Ric.exact([y, 1.0, MU])) + 12) < 1e-30)
    rs = np.geomspace(lo, hi, cv.SAMPLES)
    ys = -1 / rs
    values = K(ys, 1.0, MU)
    want = 24 + 12 * RS ** 2 / rs ** 6
    ck(f"{METRIC} flat_wall K on the wall is 24/ell^4 + 12 r_s^2/r^6", float(np.max(np.abs(values / want - 1))) < 1e-12)
    ck(f"{METRIC} flat_wall K lies within its scale", bool(np.all(values > scale["from"]) and np.all(values < scale["to"])))
    ck(f"{METRIC} flat_wall K falls outward", bool(np.all(np.diff(values) < 0)))
    ck(f"{METRIC} flat_wall K at r = 4 ell is within a part in 10^3 of 24/ell^4", abs(values[-1] / 24 - 1) < 1e-3)
    p0, p1 = cv.span(scale)
    X = cv.PLOT_W * (rs - lo) / (hi - lo)
    Y = height * (1 - (cv.position(scale, values) - p0) / (p1 - p0))
    keep = cv.thinned(X, Y, cv.TOLERANCE)
    entry = cv.curve_entry("K", "K", K, "y", params, [[float(rs[k]), float(values[k])] for k in keep])
    entry["horizontal"] = {"expression": "-1/y", "values": [float(ys[k]) for k in keep]}
    cv.agrees(ck, f"{METRIC} flat_wall K", K, [[float(ys[k]), 1.0, MU] for k in keep], [float(values[k]) for k in keep])
    K_s = cv.nr_published(src, "schwarzschild", "spherical", "kretschmann", ["r"], ["r_s"])
    reference = cv.line_curves(ck, f"{METRIC} Schwarzschild", "K", "reference", K_s, "r", x_axis, scale, height,
                               {"r_s": RS})
    schwarzschild = K_s(rs, RS)
    ck(f"{METRIC} flat_wall K is Schwarzschild's plus 24/ell^4", float(np.max(np.abs(values - schwarzschild - 24))) < 1e-6)
    marks = [{"class": "horizon", "x": RS}]
    ginv = cv.Held(cv.nr_published(src, METRIC, "flat_wall", "inverse_metric_components", ["y", "x"], ["ell", "mu"],
                                   ("y", "y")), {"x": 0.0})
    ck(f"{METRIC} flat_wall g^yy vanishes on the wall at the horizon y = -4", abs(float(ginv.exact([HORIZON, 1.0, MU]))) < 1e-30)
    return [{
        "id": "wall", "label": "The wall", "kind": "line", "system": "flat_wall", "x": x_axis,
        "panels": panels, "curves": reference + [entry], "marks": marks,
        "legend": [["line", "K", "$K$"], ["line", "reference", "Schwarzschild's black hole of the same $r_s$"],
                   ["line", "horizon", "the horizon $r = r_s$"]],
        "caption": ["The Kretschmann scalar $K$ on the wall ($x = 0$) against $r = -\\ell/y$, with "
                    "$r_s = 2\\mu\\ell$: $K = 24/\\ell^4 + 12r_s^2/r^6$, Schwarzschild's plus the $24/\\ell^4$ of "
                    "anti-de Sitter space, which it approaches far out. The Ricci scalar is $-12/\\ell^2$."],
        "settings": "$\\ell = 1$ and $\\mu = 1/8$.", "clip": "$K\\,\\ell^4$ logarithmic, from $10^{-3}$ to $10^6$."}]
