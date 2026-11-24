"""The drawings of the wall of Dvali, Gabadadze and Porrati, in units of the crossover length r_c.

The empty wall of the branch epsilon = +1 is the hyperboloid R^2 - c^2T^2 = r_c^2 of flat space of five
dimensions with two copies of its outside glued along it. On the plane of t and y of its global and
planar charts, the line psi = theta = pi/2 and x_i = 0, the point at (t, y) is cT = (r_c + |y|) sinh(ct/r_c)
and R = (r_c + |y|) cosh(ct/r_c), the metric is -(1 + |y|/r_c)^2 c^2dt^2 + dy^2 and the rays keep
ct/r_c -+ sgn(y) ln(1 + |y|/r_c). The global chart's plane is the planar chart's, so only the planar one
is drawn as a spacetime diagram.

Deffayet's Gaussian normal chart is drawn at a = r_c cosh(ct/r_c) and k = 1 on the whole line of y: the
side y > 0 is the global chart, kept by the empty wall of the branch epsilon = +1, and the side y < 0 is
the inside of the same hyperboloid out to the light cone of its centre, kept by a wall of tension
3c^4/(4 pi G r_c^2) on the branch epsilon = -1, with N = 1 + y/r_c and A = (r_c + y) cosh(ct/r_c), both
zero at the centre y = -r_c. Its rays keep ct/r_c -+ ln(1 + y/r_c).

The embedding diagram draws the closed moments of the global chart, the equator psi = theta = pi/2 of
each, from ct = -r_c to r_c, as domain_wall's are drawn: its circles grow at |drho/dy| = cosh(ct/r_c), so
each moment lies in three dimensional Minkowski space, where each side is a cone. dgp_wall.md beside
this file records each choice.
"""
import math

import numpy as np

METRIC = "dgp_wall"
REACH = 1.5                                # how far from the wall each moment is drawn, in r_c
MOMENTS = (-1.0, -0.5, 0.0, 0.5, 1.0)       # ct/r_c of the moments the other drawings mark
GN_BOX = (-1.1, 1.4, -1.25, 1.25)
PLANAR_BOX = (-REACH, REACH, -REACH, REACH)
INERTIAL_BOX = (0, 2.5, -2.5, 2.5)
GN_FIXED = {"chi": "1/2", "theta": "pi/2", "phi": "0"}
PLANAR_FIXED = {"x_1": "0", "x_2": "0", "x_3": "0"}
POLE = {"psi": "pi/2", "theta": "pi/2", "phi": "0"}
GN_INPUT = ("$a = r_c\\cosh(ct/r_c)$ and $k = 1$: on the side $y > 0$ the empty wall of the branch $\\epsilon = +1$, "
            "and on the side $y < 0$ a wall of tension $3c^4/(4\\pi Gr_c^2)$ on the branch $\\epsilon = -1$.")
WALL = "the wall, $y = 0$"


def warp(y):
    return 1 + np.abs(np.asarray(y, dtype=float))


# ---------------------------------------------------------------- spacetime diagrams

def register_null(nr):
    nr.DIAGRAMS.extend((
        nr.Diagram(METRIC, "gaussian_normal", "ty", "$t$ and $y$", ("t", "y"), GN_BOX, "$y/r_c$", "$ct/r_c$",
                   {"k": "1"}, GN_FIXED, families=nr.SIDEWAYS, cones=(6, 7), functions={"a": "cosh(t)"},
                   where="1 + y", input=GN_INPUT,
                   lines=(("shell", "r", "0", WALL), ("mark", "r", "-1", "$y = -r_c$, where $N = A = 0$"))),
        nr.Diagram(METRIC, "planar", "ty", "$t$ and $y$", ("t", "y"), PLANAR_BOX, "$y/r_c$", "$ct/r_c$",
                   {"r_c": "1"}, PLANAR_FIXED, families=nr.SIDEWAYS, cones=(6, 7), lines=(("shell", "r", "0", WALL),)),
        nr.Diagram(METRIC, "inertial", "through", "through the centre", ("T", "R"), INERTIAL_BOX, "$x/r_c$",
                   "$cT/r_c$", {"r_c": "1"}, POLE, mirror=True, tau="T", families=nr.SIDEWAYS, cones=(4, 8),
                   areal=True, curves=(("shell", "sqrt(1 + T**2)", "the wall, $R^2 - c^2T^2 = r_c^2$"),)),
    ))
    nr.CAPTIONS[(METRIC, "gaussian_normal", "ty")] = [
        "The plane of $t$ and $y$ ($\\chi = 1/2$), each ray a null geodesic, with $dy/d(ct) = \\pm(1 + y/r_c)$, "
        "the cones closing on the side $y < 0$ toward the centre, $y = -r_c$."]
    nr.CAPTIONS[(METRIC, "planar", "ty")] = [
        "The plane of $t$ and $y$ ($x_1 = x_2 = x_3 = 0$), each ray a null geodesic, with "
        "$dy/d(ct) = \\pm(1 + |y|/r_c)$, the cones opening away from the wall on both sides."]
    nr.CAPTIONS[(METRIC, "inertial", "through")] = [
        "The plane of $T$ and $R$ along a line through the centre ($\\psi = \\theta = \\pi/2$) outside the wall, the "
        "hyperbola $R^2 - c^2T^2 = r_c^2$, every ray a straight 45° line. The wall is nearest the centre at $T = 0$, "
        "where $R = r_c$, and across it lies the same region of a second copy."]
    nr.CLOSED_FORMS[(METRIC, "gaussian_normal", "ty")] = (lambda t, y: t + np.log1p(y), lambda t, y: t - np.log1p(y),
                                                         lambda t, y: y > -0.95)
    nr.CLOSED_FORMS[(METRIC, "planar", "ty")] = (lambda t, y: t + np.sign(y) * np.log(warp(y)),
                                                lambda t, y: t - np.sign(y) * np.log(warp(y)), None)
    nr.CLOSED_FORMS[(METRIC, "inertial", "through")] = (lambda T, R: T + R, lambda T, R: T - R, None)


def _inertial_moment(m, edge=INERTIAL_BOX[1], N=400):
    """A moment ct/r_c of the global chart in the inertial chart, on both sides of the centre:
    cT = |R| tanh(ct/r_c), from the wall at |R| = cosh(ct/r_c) out to the embedding's reach,
    (1 + 1.5) cosh(ct/r_c), or to the edge of the drawing."""
    lo = math.cosh(m.time)
    hi = min((1 + m.reach("global", "y")[1]) * lo, edge)
    R = np.linspace(lo, hi, N)
    return [np.column_stack([R * math.tanh(m.time), -R[::-1]]), np.column_stack([R * math.tanh(m.time), R])]


def register_slices(sl):
    sl.FLAT[(METRIC, "planar", "ty")] = lambda: sl.one(
        METRIC, lambda m: [[(m.time, m.reach("global", "y")[0]), (m.time, m.reach("global", "y")[1])]])
    sl.FLAT[(METRIC, "gaussian_normal", "ty")] = lambda: sl.one(
        METRIC, lambda m: [[(m.time, 0.0), (m.time, min(m.reach("global", "y")[1], GN_BOX[1]))]])
    sl.FLAT[(METRIC, "inertial", "through")] = lambda: sl.one(METRIC, _inertial_moment)
    sl.FLAT_METRICS.add(METRIC)


# ---------------------------------------------------------------- embedding diagram

def register_embedding(eb):
    eb.DRAWN[METRIC] = lambda ck, src: embedding(eb, ck, src)
    eb.CAPTIONS[(METRIC, "moments")] = [
        "The equator ($\\psi = \\theta = \\pi/2$) of the global chart from $ct = -r_c$ to $r_c$ in three dimensional "
        "Minkowski space, its circles of radius $(r_c + |y|)\\cosh(ct/r_c)$, each side a cone, $Z = y\\sinh(ct/r_c)$, "
        "meeting the other at the wall."]


def embedding(eb, ck, src):
    """The equator of a moment of the global chart, g_yy = 1 and g_phiphi = (r_c + |y|)^2 cosh^2(ct/r_c),
    at r_c = 1, from |y| = 1.5 r_c on one side through the wall to the same distance on the other. Its
    circles grow at |drho/dy| = cosh(ct/r_c), faster than the distance out to them at every moment but
    ct = 0, which is checked, so it is drawn in three dimensional Minkowski space, dX^2 + dY^2 - dZ^2,
    where it climbs at dZ/dy = sinh(ct/r_c): rho = (1 + |y|) cosh(ct/r_c) and Z = y sinh(ct/r_c), two cones
    joined at their narrowest circle, the wall, where the profile turns. Each side is its own moment in
    the inertial chart of its copy, Z = cT and rho = R, the cone cT = R tanh(ct/r_c) from the wall out;
    the side y < 0 is the mirror image of the side y > 0 in the wall, an isometry of Minkowski space. At
    ct = 0 both cones are the flat plane outside the circle of radius r_c, and the moment is that plane
    taken twice, joined at the circle. Drawn at five moments from ct = -r_c to r_c and played as a movie
    with a frame every 0.05 r_c of ct."""
    size = 8.0
    rings = (0.5, 1.0)

    def moment(t):
        sl = eb.Slice(src, METRIC, "global", "y", "\\phi", {"t": repr(t), **POLE}, {"r_c": 1}, space="minkowski")
        sense = -1 if t < 0 else 1
        where = f"DGP wall, ct = {t:+.2f} r_c"
        edge = "the copy cut at $|y| = 3r_c/2$"
        near = eb.Piece("near", "sheet", sl, -REACH, 0.0, 0.0, sense,
                        (("edge", edge), ("crease", "the wall, where the two copies meet")),
                        [(-x, "r", None) for x in reversed(rings)] + [(0.0, "surface", None)], size)
        far = eb.Piece("far", "sheet2", sl, 0.0, REACH, float(near.z[-1]), sense,
                       (("crease", "the wall"), ("edge", edge)), [(x, "r2", None) for x in rings], size)
        shift = -float(near.z[-1])
        near.z, far.z = near.z + shift, far.z + shift
        ck.add(f"{where}, the two copies meet at the wall: one point",
               float(np.max(np.abs(np.array(near.at(0.0)) - far.at(0.0)))), eb.JOIN)
        flat = eb.Slice(src, METRIC, "global", "y", "\\phi", {"t": repr(t), **POLE}, {"r_c": 1})
        ys = np.concatenate([np.linspace(-REACH, -1e-3, 200), np.linspace(1e-3, REACH, 200)])
        if t == 0:
            ck.plane(f"{where}, the equator is a plane", flat, ys)
        else:
            ck.stops(f"{where}, the equator in flat space", flat, ys)
        for piece in (near, far):
            ck.isometry(f"{where}, {piece.id}", piece)
            ck.radius(f"{where}, {piece.id}: rho = (r_c + |y|) cosh(ct/r_c)", piece,
                      lambda y: (1 + np.abs(y)) * math.cosh(t), size)
            ck.form(f"{where}, {piece.id}: Z = y sinh(ct/r_c)", piece, lambda y: y * math.sinh(t), size)
        return eb.Surface([near, far], label=f"$ct/r_c = {t:.2f}$", time=t)

    values, keys = eb.movie_values(list(MOMENTS), 0.05)
    frames = [moment(t) for t in values]
    surfaces = [frames[i] for i in keys]
    fig = eb.movie_figure(frames, {"sheet": "cover"}, size, meridians=12)
    fig.legend("fill", "cover", "the copy $y < 0$")
    fig.legend("line", "r", "$y = -r_c/2$ and $-r_c$")
    fig.legend("line", "r2", "$y = r_c/2$ and $r_c$")
    fig.legend("line", "surface", "the wall")
    fig.legend("line", "meridian", "$\\phi$ constant")
    return [eb.view("moments", "The closed moments", "$r_c$", surfaces, fig.done(),
                    movie=eb.movie(frames, "$ct/r_c$", values),
                    settings="$r_c = 1$, lengths measured with $dX^2 + dY^2 - dZ^2$.",
                    stops=["At every moment but $ct = 0$ the circles grow faster than the distance out to them, "
                           "so the slice lies in Minkowski space."])]


# ---------------------------------------------------------------- conformal diagrams

def register_conformal(cf):
    cf.DRAWN[METRIC] = lambda ck, src: conformal(cf, ck, src)
    cf.CAPTIONS[(METRIC, "gaussian_normal")] = [
        "The plane of $t$ and $y$ ($\\chi = 1/2$) at $a = r_c\\cosh(ct/r_c)$ and $k = 1$, each point a 3-sphere, one "
        "copy of flat space with both sides of the wall."]
    cf.CAPTIONS[(METRIC, "global")] = [
        "The plane of $t$ and $y$ ($\\psi = \\theta = \\pi/2$), each point a 3-sphere, two copies of the outside of "
        "the wall glued along it, from $\\mathscr{I}^-$ to $\\mathscr{I}^+$."]
    cf.CAPTIONS[(METRIC, "planar")] = [
        "The plane of $t$ and $y$ ($x_1 = x_2 = x_3 = 0$), each point a flat 3-space, two copies of the outside of "
        "the wall glued along it, from $\\mathscr{I}^-$ to $\\mathscr{I}^+$."]
    cf.CAPTIONS[(METRIC, "inertial")] = [
        "The plane of $T$ and $R$ ($\\psi = \\theta = \\pi/2$), each point a 3-sphere, the outside of the wall, the "
        "hyperbola $R^2 - c^2T^2 = r_c^2$, and its mirror image, the second copy."]


def conformal(cf, ck, src):
    """Four views, in units of r_c, each by Minkowski's maps p = arctan(cT - R) and q = arctan(cT + R) of
    flat space's plane of T and R, which send the wall, tan p tan q = -1, to the line X = pi/2 from the
    middle of I- at T = -pi/2 to the middle of I+ at T = pi/2. The side y > 0 is the outside of the wall,
    the part X >= pi/2 of Minkowski's triangle, a triangle with its corner i^0 at (pi, 0), and the side
    y < 0 of the global, planar and inertial charts is its mirror image in the wall, X -> pi - X, which is
    p = arctan(cT + R) - pi/2, q = arctan(cT - R) + pi/2; the two copies make the diamond of the plane's
    two edges, each meeting the other along the wall at one T. On the plane of t and y the point at (t, y)
    is cT = (1 + |y|) sinh(ct), R = (1 + |y|) cosh(ct). The Gaussian normal chart at a = cosh(ct), k = 1 is
    one copy of flat space, cT = (1 + y) sinh(ct) and R = (1 + y) cosh(ct): its side y > 0 is the outside,
    and its side -1 < y < 0 the inside of the wall where R > c|T|, between the wall and the light cone of
    the centre, which t reaches only as it runs to either end, its centre y = -1 the point R = T = 0."""
    PI, HALF, Q4 = cf.PI, cf.HALF, cf.Q4

    def right(T, R):
        return cf.mink_pq(T, R)

    def left(T, R):
        T, R = np.asarray(T, dtype=float), np.asarray(R, dtype=float)
        return np.arctan(T + R) - HALF, np.arctan(T - R) + HALF

    def sides(T, R, y):
        pr, qr = right(T, R)
        pl, ql = left(T, R)
        y = np.asarray(y, dtype=float)
        return np.where(y > 0, pr, pl), np.where(y > 0, qr, ql)

    def planar(t, y):
        t, y = np.asarray(t, dtype=float), np.asarray(y, dtype=float)
        return sides(warp(y) * np.sinh(t), warp(y) * np.cosh(t), y)

    def gaussian_normal(t, y):
        t, y = np.asarray(t, dtype=float), np.asarray(y, dtype=float)
        return right((1 + y) * np.sinh(t), (1 + y) * np.cosh(t))

    gn = cf.Plane(src, METRIC, "gaussian_normal", ("t", "y"), GN_FIXED, {"k": 1}, functions={"a": "cosh(t)"})
    pl = cf.Plane(src, METRIC, "planar", ("t", "y"), PLANAR_FIXED, {"r_c": 1})
    gl = cf.Plane(src, METRIC, "global", ("t", "y"), POLE, {"r_c": 1})
    ip = cf.Plane(src, METRIC, "inertial", ("T", "R"), POLE, {"r_c": 1})
    future = lambda t, x: (1, 0)  # noqa: E731
    for name, plane in (("planar", pl), ("global", gl)):
        for side, sign in (("<", -1), (">", 1)):
            ck.chart(f"DGP wall {name}, y {side} 0", plane, planar, ck.uniform(-4, 4), sign * ck.uniform(0.001, 6),
                     future)
    ck.chart("DGP wall Gaussian normal", gn, gaussian_normal, ck.uniform(-4, 4), ck.uniform(-0.999, 6), future)
    T = ck.uniform(-10, 10)
    ck.chart("DGP wall inertial", ip, right, T, np.sqrt(1 + T ** 2) * ck.uniform(1.001, 6), future)
    ck.chart("DGP wall inertial, mirrored", ip, left, T, np.sqrt(1 + T ** 2) * ck.uniform(1.001, 6), future)
    T = np.linspace(-20, 20, 81)
    wall = np.sqrt(1 + T ** 2)
    ck.limit("DGP wall: the wall R^2 - c^2T^2 = r_c^2 is the line X = pi/2 on both sides",
             np.concatenate([cf.xt(*right(T, wall))[0], cf.xt(*left(T, wall))[0]]), np.full(162, HALF), 1e-12)
    t = np.linspace(-4, 4, 81)
    ck.limit("DGP wall: both copies reach the wall at one T for every t",
             np.concatenate(cf.xt(*planar(t, np.full(81, -1e-15)))), np.concatenate(cf.xt(*planar(t, np.full(81, 1e-15)))),
             1e-12)
    ck.limit("DGP wall: the wall ends at the middle of I+ and I-, (X, T) = (pi/2, +-pi/2)",
             cf.xt(*right(np.array([1e9, -1e9]), np.sqrt(1 + 1e18)))[1], [HALF, -HALF], 1e-8)
    X, Tc = cf.xt(*gaussian_normal(t, np.full(81, -1 + 1e-12)))
    ck.limit("DGP wall: the Gaussian normal chart's y = -r_c is the centre at T = 0", np.abs(X) + np.abs(Tc),
             np.zeros(81), 1e-9)
    T = ck.uniform(-5, 5, 50)
    ck.finite("DGP wall: the inertial chart is flat outside the wall", ip.kretschmann(T, np.sqrt(1 + T ** 2) + 1))

    diamond = [[0, 0], [HALF, HALF], [PI, 0], [HALF, -HALF]]
    outside = [[HALF, -HALF], [HALF, HALF], [PI, 0]]
    inside = [[0, 0], [HALF, HALF], [HALF, -HALF]]
    triangle = [[0, -PI], [0, PI], [PI, 0]]
    S = cf.spread(-np.inf, np.inf, 600, 10)

    def infinity(v, points):
        for at, text, anchor, dy, dx in points:
            v.layers.append({"kind": "point", "class": "infinity", "at": [round(at[0], 4), round(at[1], 4)]})
            v.label_xt(at, text, anchor, dx=dx, dy=dy)

    def two_copies(v, cover):
        v.fill("region", diamond)
        v.fill("cover", diamond)
        v.line("scri", [[[0, 0], [HALF, HALF]], [[HALF, HALF], [PI, 0]], [[0, 0], [HALF, -HALF]], [[HALF, -HALF], [PI, 0]]])
        v.line("surface", [[[HALF, -HALF], [HALF, HALF]]])
        infinity(v, (((0, 0), "$i^0$", "r", 0, -6), ((PI, 0), "$i^0$", "l", 0, 6)))
        v.label_xt([Q4, Q4], "$\\mathscr{I}^+$", "br", dx=-4, dy=-3)
        v.label_xt([3 * Q4, Q4], "$\\mathscr{I}^+$", "bl", dx=4, dy=-3)
        v.label_xt([Q4, -Q4], "$\\mathscr{I}^-$", "tr", dx=-4, dy=3)
        v.label_xt([3 * Q4, -Q4], "$\\mathscr{I}^-$", "tl", dx=4, dy=3)
        v.legend("cover", cover)
        v.legend("surface", "the wall")
        v.legend("scri", "null infinity")

    views = []
    box = [-0.45, PI + 0.45, -HALF - 0.45, HALF + 0.45]
    for vid, label, cover in (("global", "Global", "two copies of the outside"),
                              ("planar", "Planar", "two copies of the outside")):
        v = cf.View(vid, label, box, vid)
        two_copies(v, cover)
        cf.grid(v, "r", lambda c, s: planar(s, c), (-1.0, -0.5, 0.5, 1.0), S)
        span = np.linspace(-40, 40, 1601)
        for t0 in (-1, 0, 1):
            for half in (span[span < 0], span[span > 0]):
                v.curve("t", *planar(np.full_like(half, t0), half))
        v.legend("r", "$y = \\pm r_c/2$ and $\\pm r_c$")
        v.legend("t", "$ct/r_c = -1$, $0$, and $1$")
        views.append(v)

    v = cf.View("inertial", "Inertial", box, "inertial")
    two_copies(v, "the outside and its mirror image")
    for R0 in (1.5, 2.0, 3.0):
        T = np.sinh(np.linspace(-8, 8, 1200))
        out = np.where(R0 >= np.sqrt(1 + T ** 2), T, np.nan)
        v.curve("r", *right(out, np.full_like(T, R0)))
        v.curve("r", *left(out, np.full_like(T, R0)))
    for T0 in (-1, 0, 1):
        R = np.sqrt(1 + T0 ** 2) * np.exp(np.linspace(0, 10, 600))
        v.curve("t", *right(np.full_like(R, T0), R))
        v.curve("t", *left(np.full_like(R, T0), R))
    v.legend("r", "$R = 3r_c/2$, $2r_c$, and $3r_c$")
    v.legend("t", "$cT/r_c = -1$, $0$, and $1$")
    inertial = v

    v = cf.View("gaussian_normal", "Deffayet, Gaussian Normal", [-0.45, PI + 0.45, -PI - 0.25, PI + 0.25], "gaussian_normal")
    v.fill("region", triangle)
    v.fill("cover", outside)
    v.fill("cover2", inside)
    v.line("centre", [[[0, -PI], [0, PI]]])
    v.line("scri", [[[0, PI], [PI, 0]], [[PI, 0], [0, -PI]]])
    v.line("horizon", [[[0, 0], [HALF, HALF]], [[0, 0], [HALF, -HALF]]])
    v.line("surface", [[[HALF, -HALF], [HALF, HALF]]])
    infinity(v, (((0, PI), "$i^+$", "b", -6, 0), ((0, -PI), "$i^-$", "t", 6, 0), ((PI, 0), "$i^0$", "l", 0, 6)))
    cf.grid(v, "r", lambda c, s: gaussian_normal(s, c), (-0.5, 0.5, 1.0), S)
    for t0 in (-1, 0, 1):
        span = np.linspace(-0.999999, 40, 2001)
        v.curve("t", *gaussian_normal(np.full_like(span, t0), span))
    v.legend("cover", "$y > 0$, kept on the branch $\\epsilon = +1$")
    v.legend("cover2", "$y < 0$, kept on the branch $\\epsilon = -1$")
    v.legend("surface", "the wall")
    v.legend("horizon", "the centre's light cone")
    v.legend("centre", "the centre")
    v.legend("scri", "null infinity")
    v.legend("r", "$y = -r_c/2$, $r_c/2$, and $r_c$")
    v.legend("t", "$ct/r_c = -1$, $0$, and $1$")
    views = [v] + views + [inertial]

    # Each moment of the embedding, ct/r_c = t_k, out to |y| = 1.5 r_c on either side of the wall.
    y = np.concatenate([np.linspace(-REACH, -1e-12, 400), np.linspace(1e-12, REACH, 400)])
    yp = np.linspace(1e-12, REACH, 400)
    for v in views:
        v.set(settings="$r_c$, the unit of every length and of $ct$.")
        for m in cf.slices.moments(METRIC):
            if v.d["id"] == "gaussian_normal":
                v.slice(m, [gaussian_normal(np.full_like(yp, m.time), yp)])
            elif v.d["id"] == "inertial":
                R = (1 + yp) * math.cosh(m.time)
                T = R * math.tanh(m.time)
                v.slice(m, [right(T, R), left(T, R)])
            else:
                v.slice(m, [planar(np.full_like(y, m.time), y)])
    return views
