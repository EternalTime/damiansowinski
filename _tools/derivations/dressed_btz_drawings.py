"""The drawings of Martinez and Zanelli's conformally dressed black hole, at l = 1, the unit of every
length, and B = l/2, so that the horizon is r_+ = 2B = l.

F = (r + B)^2(r - 2B)/(r l^2) has a simple root at r = 2B and a double one at r = -B, and

    r_* = int dr/F = l^2 (2/(9B)) ln|(r - 2B)/(r + B)| - l^2/(3(r + B)),

which vanishes as r -> infinity and is finite, -(l^2/B)(1/3 - (2/9) ln 2), at r = 0. The surface
gravity is kappa = F'(2B)/2 = 9B/(4 l^2), and 2 kappa (2 l^2/(9B)) = 1, so Kruskal's
U = -exp(-kappa u), V = exp(kappa v), u, v = ct -+ r_*, give

    UV = ((2B - r)/(r + B)) exp(-3B/(2(r + B))),

analytic through the horizon: -1 on the conformal boundary r -> infinity, 0 on the horizon and
2 exp(-3/2) = 0.446 at the singularity r = 0. With p = arctan U and q = arctan V the boundary
tan p tan q = -1 is the pair of vertical lines X = +-pi/2 of the BTZ hole's square, and the
singularity tan p tan q = 2 exp(-3/2) is a curve from corner to corner of that square that bows
toward its centre, reaching T = 2 arctan(sqrt(2) exp(-3/4)) = 1.178 on the axis, below the pi/2
of a straight line. The ingoing chart's v is fixed by v - t -> 0 as r -> infinity, as its
convention says, so the moment t = 0 is cv = r_* on its plane.

The embedding diagram draws the moment t = 0 through both exteriors, as the BTZ hole's is drawn:
in flat space from the bifurcation circle r_+ out to the circle where F = 1, and in three
dimensional Minkowski space beyond it, where 1/F < 1. The curvature figure is the stored
Kretschmann scalar against r; the Ricci scalar is -6/l^2 at every point.
"""
import math

import numpy as np

METRIC = "dressed_btz"
STATIC, INGOING = "static", "eddington_finkelstein_ingoing"
L, B = 1.0, 0.5
PARAMS = {"l": 1, "B": "1/2"}
FIXED = {"theta": "0"}
RP = 2 * B
KAPPA = 9 * B / (4 * L ** 2)
TOP = 3.0
SETTINGS = "$l = 1$, the unit of every length, and $B = l/2$, so that $r_+ = l$."


def F(r):
    r = np.asarray(r, dtype=float)
    return (r + B) ** 2 * (r - 2 * B) / (r * L ** 2)


def rstar(r):
    """The tortoise coordinate, zero as r -> infinity."""
    r = np.asarray(r, dtype=float)
    with np.errstate(divide="ignore"):
        return L ** 2 * (2 / (9 * B)) * np.log(np.abs((r - 2 * B) / (r + B))) - L ** 2 / (3 * (r + B))


def UV(r):
    """Kruskal's UV, one formula for every r > 0."""
    r = np.asarray(r, dtype=float)
    return (2 * B - r) / (r + B) * np.exp(-3 * B / (2 * (r + B)))


def bisect(f, lo, hi, n=200):
    a, b = lo, hi
    for _ in range(n):
        m = 0.5 * (a + b)
        if (f(a) < 0) == (f(m) < 0):
            a = m
        else:
            b = m
    return 0.5 * (a + b)


LEVEL = bisect(lambda r: float(F(r)) - 1.0, RP + 1e-9, TOP)   # F = 1, where the slice lies level


# ---------------------------------------------------------------- spacetime diagrams

def register_null(nr):
    nr.DIAGRAMS.append(nr.Diagram(
        METRIC, STATIC, "radial", "$t$ and $r$", ("t", "r"), (0, 3, -1.5, 1.5), "$r/l$", "$ct/l$",
        PARAMS, FIXED, orient="ingoing"))
    nr.DIAGRAMS.append(nr.Diagram(
        METRIC, INGOING, "radial", "$v$ and $r$", ("v", "r"), (0, 3, -3, 0), "$r/l$", "$(cv - r)/l$",
        PARAMS, FIXED, to_display=nr.FINKELSTEIN_IN, orient="ingoing"))
    nr.CAPTIONS[(METRIC, STATIC, "radial")] = [
        "The plane of $t$ and $r$ ($\\theta = 0$) of the dressed black hole ($B = l/2$), the same at every "
        "$\\theta$, each point in the diagram a circle. The cones close at the horizon $r_+ = 2B$, and a ray "
        "reaches $r \\to \\infty$ in a finite time. Inside $r_+$ the future is taken from the ingoing chart, and "
        "every cone points to the singularity $r = 0$."]
    nr.CAPTIONS[(METRIC, INGOING, "radial")] = [
        "The plane of $v$ and $r$ ($\\theta = 0$) of the dressed black hole ($B = l/2$), drawn with $cv - r$ "
        "upward so that the ingoing rays run at 45°. The outgoing rays stand vertical at the horizon "
        "$r_+ = 2B$, which the cones cross smoothly; inside it every future directed ray ends at $r = 0$."]


def slice_table(slices):
    """The moment t = 0 of the embedding view on each chart's spacetime diagram: along r in the
    static chart, and as cv = r_* in the ingoing chart, crowding toward the horizon."""
    def static():
        return slices.one(METRIC, lambda m: slices.along(0.0, *m.reach(STATIC, "r")), view_id="throat")

    def ingoing():
        def lines(m):
            lo, hi = m.reach(STATIC, "r")
            r = slices.near(lo, hi)
            return [np.column_stack([rstar(r), r])]
        return slices.one(METRIC, lines, view_id="throat")
    return {(METRIC, STATIC, "radial"): static, (METRIC, INGOING, "radial"): ingoing}


# ---------------------------------------------------------------- embedding diagram

def register_embedding(em):
    em.DRAWN[METRIC] = lambda ck, src: embedding_views(ck, src, em)
    em.CAPTIONS[(METRIC, "throat")] = [
        f"The moment $t = 0$ of the dressed hole through both exteriors, joined at the bifurcation "
        f"circle $r = r_+$, in Euclidean space out to $r = {LEVEL:.2f}\\,l$ and in Minkowski space beyond."]


def embedding_views(ck, src, em):
    """Both exteriors of the moment t = 0, as the BTZ hole's is drawn, each a piece in flat space from
    the throat to the level circle F = 1 and a piece in Minkowski space from there to r = 3l."""
    name = "dressed BTZ"
    sl = em.Slice(src, METRIC, STATIC, "r", "\\theta", {"t": 0}, PARAMS)
    msl = em.Slice(src, METRIC, STATIC, "r", "\\theta", {"t": 0}, PARAMS, space="minkowski")
    # The level circle F = 1 is the largest root of (r + 1/2)^2 (r - 1) = r, 4r^3 - 7r - 1 = 0, named
    # exactly so that the joins are taken at the root itself.
    import sympy as sp
    x = sp.Symbol("x")
    exact = sp.CRootOf(4 * x ** 3 - 7 * x - 1, 2)
    ck.add(f"{name}: the level circle is the largest root of 4r^3 - 7r - 1", abs(float(exact) - LEVEL), 1e-12)
    sl.known = msl.known = {LEVEL: exact}
    rp = sl.horizons()[0]
    size = 2 * TOP
    ck.add(f"{name}: the horizon is at r_+ = 2B = l", abs(rp - RP), 1e-12)
    ck.add(f"{name}: the surface lies level where F = 1, at r = {LEVEL:.6f} l",
           abs(float(sl.defect_at(np.array([LEVEL]))[0])), 1e-9)
    ck.stops(f"{name}, beyond F = 1 in flat space", sl, np.linspace(LEVEL, 40, 402)[1:])
    inward = np.linspace(rp, LEVEL, 402)[1:-1]
    ck.add(f"{name}: between r_+ and F = 1 no surface in Minkowski space carries the slice, (drho/dr)^2 - g_rr < 0",
           float(max(0.0, np.max(-msl.defect_at(inward)))), 0.0)
    if not np.all(msl.defect_at(inward) > 0):
        ck.items[-1]["ok"] = False

    join = "where $g_{rr} = 1$ the surface lies level, in Euclidean space nearer the horizon and in Minkowski space beyond"
    edge = "the sheet runs on toward a light cone, to $r \\to \\infty$"
    near = em.Piece("exterior", "sheet", sl, rp, LEVEL, 0.0, 1,
                    (("throat", "the throat $r = r_+$, the bifurcation circle, where the other exterior begins"),
                     ("join", join)),
                    [(rp, "horizon", "$r = r_+$"), (1.2, "r", None), (LEVEL, "space", None)], size)
    out = em.Piece("exterior_minkowski", "sheet", msl, LEVEL, TOP, near.at(LEVEL)[1], 1, (("join", join), ("edge", edge)),
                   [(2.0, "r", None), (TOP, "r", None)], size)
    far = em.Piece("other_exterior", "sheet2", sl, rp, LEVEL, 0.0, -1, (("throat", "the throat $r = r_+$"), ("join", join)),
                   [(1.2, "r2", None), (LEVEL, "space", None)], size)
    far_out = em.Piece("other_exterior_minkowski", "sheet2", msl, LEVEL, TOP, far.at(LEVEL)[1], -1,
                       (("join", join), ("edge", edge)), [(2.0, "r2", None), (TOP, "r2", None)], size)
    for p in (near, out, far, far_out):
        space = "in Minkowski space" if p.sl.lorentz else "in flat space"
        ck.isometry(f"{name}, {p.id} {space}", p)
        ck.radius(f"{name}, {p.id}, rho = r {space}", p, lambda r: r, size)
    ck.join(f"{name}, the two sheets at the throat", near, rp, far, rp)
    for a, b in ((near, out), (far, far_out)):
        ck.join(f"{name}, {a.id} in flat space and {b.id} in Minkowski space where F = 1", a, LEVEL, b, LEVEL)
        pa, pb = a.data()["points"][-1], b.data()["points"][0]
        ck.add(f"{name}, {a.id} and {b.id} as written: one point where F = 1",
               max(abs(pa[0] - pb[0]), abs(pa[1] - pb[1]), abs(pa[2] - pb[2])), 10.0 ** -min(a.decimals, b.decimals))

    surface = em.Surface([near, out, far, far_out])
    fig = em.figure_of([surface], {"sheet": "cover"}, size)
    em.ring_label(fig, [0, 0, 0], rp, 0.0, "$r = r_+$", dx=14)
    em.ring_label(fig, [0, 0, 0], *out.at(TOP), "$3l$")
    fig.legend("fill", "cover", "the exterior $t$ and $r$ cover")
    fig.legend("line", "r", "$r$ constant")
    fig.legend("line", "r2", "$r$ constant, the other exterior")
    fig.legend("line", "horizon", "the throat $r = r_+$")
    fig.legend("line", "space", "$g_{rr} = 1$")
    fig.legend("line", "meridian", "$\\theta$ constant, every $15°$")
    return [em.view("throat", "Through the throat", "$l$", [surface], fig.done(),
                    settings="$l = 1$ and $B = l/2$, so $r_+ = l$; beyond $g_{rr} = 1$, lengths in $dX^2 + dY^2 - dZ^2$.")]


# ---------------------------------------------------------------- conformal diagram

S_SINGULAR = 2 * math.exp(-1.5)          # UV at r = 0


def pq(cf, cell, t, r):
    """(p, q) of points (t, r) of a cell, t the static time of its region: U = -+exp(-kappa u),
    V = +-exp(kappa v), u, v = ct -+ r_*, p = arctan U, q = arctan V, the cells of a Tower."""
    t = np.asarray(t, dtype=float)
    rs = rstar(r)
    u, v = t - rs, t + rs

    def G(w):
        return cf.atan_exp(-KAPPA * np.asarray(w, dtype=float))
    return {"I": (-G(u), G(-v)), "II": (G(u), G(-v)), "IV": (-G(u), -G(-v)), "I'": (G(u), -G(-v))}[cell]


def ingoing(cf, w, r):
    """(p, q) of points (cv, r) of the ingoing chart: V = exp(kappa v), U = UV(r)/V."""
    w, r = np.asarray(w, dtype=float), np.asarray(r, dtype=float)
    return np.arctan(UV(r) * np.exp(-KAPPA * w)), cf.atan_exp(KAPPA * w)


def register_conformal(cf):
    cf.DRAWN[METRIC] = lambda ck, src: conformal_views(ck, src, cf)
    end = ("each point a circle, with a timelike boundary $r \\to \\infty$ and a spacelike singularity $r = 0$ "
           "bowed toward the centre.")
    cf.CAPTIONS[(METRIC, "static")] = [
        "The dressed black hole ($B = l/2$) maximally extended in Kruskal's coordinates, " + end]
    cf.CAPTIONS[(METRIC, "ingoing")] = [
        "The dressed black hole ($B = l/2$) in its ingoing chart, one exterior and the black hole, " + end]


def conformal_views(ck, src, cf):
    """The static and ingoing charts on Kruskal's extension of the module docstring: the conformal
    boundary on the straight lines X = +-pi/2 and the singularity on tan p tan q = 2 exp(-3/2)."""
    import sympy as sp
    HALF = cf.HALF
    name = "dressed BTZ"
    # r_* and UV against the published metric: dr_*/dr = 1/F = g^rr of the static chart, and
    # exp(2 kappa r_*) = |UV|.
    r = sp.Symbol("r", positive=True)
    Bs, ls = sp.Rational(1, 2), sp.Integer(1)
    rs_sym = ls ** 2 * sp.Rational(2, 9) / Bs * sp.log((r - 2 * Bs) / (r + Bs)) - ls ** 2 / (3 * (r + Bs))
    st = cf.Plane(src, METRIC, STATIC, ("t", "r"), FIXED, PARAMS)
    gi = st.gi[1, 1].subs(st.x1, r)
    ck.limit(f"{name}: dr_*/dr is the published g^rr inverted", [float(sp.simplify(sp.diff(rs_sym, r) * gi - 1))], [0.0], 1e-15)
    rr = np.array([0.05, 0.3, 0.7, 1.4, 3.0, 20.0])
    ck.limit(f"{name}: exp(2 kappa r_*) is |UV|", np.exp(2 * KAPPA * rstar(rr)), np.abs(UV(rr)), 1e-12)
    ck.limit(f"{name}: kappa = F'(r_+)/2 = 9B/(4l^2)", [KAPPA], [float(sp.diff(gi, r).subs(r, 1) / 2)], 1e-12)
    ck.limit(f"{name}: UV is -1 as r -> infinity and 2 exp(-3/2) at r = 0", UV(np.array([1e12, 0.0])),
             [-1.0, S_SINGULAR], 1e-9)
    for cell, ranges, future, where in (("I", (1.001, 40), (1, 0), "exterior"), ("II", (0.01, 0.999), (0, -1), "black hole"),
                                        ("IV", (0.01, 0.999), (0, 1), "white hole"),
                                        ("I'", (1.001, 40), (-1, 0), "other exterior")):
        ck.chart(f"{name} static, {where}", st, lambda t, r, c=cell: pq(cf, c, t, r),
                 ck.uniform(-4, 4), ck.uniform(*ranges), lambda t, r, d=future: d)
    ein = cf.Plane(src, METRIC, INGOING, ("v", "r"), FIXED, PARAMS)
    ck.chart(f"{name} ingoing Eddington-Finkelstein", ein, lambda w, r: ingoing(cf, w, r),
             ck.uniform(-4, 4), ck.uniform(0.01, 40), lambda w, r: (1, -60))
    ck.limit(f"{name}: the ingoing and static coordinates put one event at one point",
             ingoing(cf, 0.4 + rstar(3.0), 3.0), pq(cf, "I", 0.4, 3.0), 1e-12)
    p, q = pq(cf, "I", np.array([-0.5, 0.0, 0.5]), np.full(3, 1e9))
    ck.limit(f"{name}: r -> infinity lands on the vertical line X = pi/2", q - p, [HALF] * 3, 1e-6)
    p, q = pq(cf, "II", np.array([-1.0, 0.0, 1.0]), np.full(3, 1e-9))
    ck.limit(f"{name}: r -> 0 lands on tan p tan q = 2 exp(-3/2)", np.tan(p) * np.tan(q), [S_SINGULAR] * 3, 1e-6)
    mid = 2 * math.atan(math.sqrt(S_SINGULAR))
    ck.limit(f"{name}: the singularity crosses the axis at T = 2 arctan sqrt(2 exp(-3/2)) = 1.178", (p + q)[1], mid, 1e-6)
    p, q = pq(cf, "I", np.array([0.5]), np.array([RP * (1 + 1e-12)]))
    ck.limit(f"{name}: r -> r_+ at fixed t lands on the bifurcation circle", cf.point(p[0], q[0]), [0, 0], 1e-4)
    K = st.kretschmann
    ck.diverges(f"{name}: the Kretschmann scalar diverges at r = 0", K(0, 1e-2), K(0, 1e-3))
    ck.finite(f"{name}: the Kretschmann scalar is finite at r = r_+", K(np.zeros(3), np.array([0.999, 1, 1.001])))

    box = [-HALF - 0.3, HALF + 0.3, -HALF - 0.25, HALF + 0.25]
    a = np.linspace(1e-6, HALF - 1e-6, 400)
    top = [cf.point(x, math.atan(S_SINGULAR / math.tan(x))) for x in a]      # from the left corner to the right
    bottom = [[X, -T] for X, T in top]
    right = [[HALF, -HALF], [HALF, HALF]]
    whole = [[-HALF, -HALF]] + bottom + [[HALF, -HALF], [HALF, HALF]] + top[::-1] + [[-HALF, HALF]]
    R_OUT, R_IN, TS = (1.2, 1.5, 2.0, 4.0), (0.3, 0.6, 0.9), (-0.8, -0.4, -0.2, 0, 0.2, 0.4, 0.8)

    def edges(v):
        v.line("boundary", [right, [[-HALF, -HALF], [-HALF, HALF]]])
        v.line("horizon", [[[-HALF, -HALF], [HALF, HALF]], [[HALF, -HALF], [-HALF, HALF]]])
        v.line("singular", [top, bottom], zig=True)
        v.label_xt([0, mid], "$r = 0$", "b", dy=-8)
        v.label_xt([0, -mid], "$r = 0$", "t", dy=8)
        v.label_xt([HALF, 0.9], "$r \\to \\infty$", "l", "small", dx=6)
        v.label_xt([-HALF, 0.9], "$r \\to \\infty$", "r", "small", dx=-6)
        v.label_xt([-cf.Q4, cf.Q4], "$r_+$", "tr", "small", dx=-6, dy=2)
        v.legend("horizon", "the horizon $r_+$")
        v.legend("boundary", "$r \\to \\infty$")
        v.legend("singular", "the singularity $r = 0$")

    views = []
    t = cf.spread(-np.inf, np.inf, 500, 6)
    v = cf.View("static", "Static", box, STATIC)
    v.fill("region", whole)
    v.fill("cover", [[0, 0], [HALF, -HALF], [HALF, HALF]])
    for rv in R_OUT:
        v.curve("r", *pq(cf, "I", t, np.full_like(t, rv)))
    rr = cf.spread(RP, np.inf, 500, 14)
    for tt in TS:
        v.curve("t", *pq(cf, "I", np.full_like(rr, tt), rr))
    edges(v)
    v.legend("cover", "covered by $t$ and $r > r_+$")
    v.legend("r", "$r$ constant")
    v.legend("t", "$ct$ constant")
    views.append(v)

    v = cf.View("ingoing", "Ingoing Eddington-Finkelstein", box, INGOING)
    v.fill("region", whole)
    v.fill("cover", [[0, 0], [HALF, -HALF], [HALF, HALF]] + top[::-1] + [[-HALF, HALF]])
    for rv in R_OUT + R_IN:
        v.curve("r", *ingoing(cf, t, np.full_like(t, rv)))
    rr = cf.spread(0, np.inf, 600, 14)
    for w in (-0.4, 0, 0.4, 0.8, 1.2):
        v.curve("null", *ingoing(cf, np.full_like(rr, w), rr))
    edges(v)
    v.legend("cover", "covered by $v$ and $r > 0$")
    v.legend("r", "$r$ constant")
    v.legend("null", "$v$ constant, an ingoing ray")
    views.append(v)
    moment, = cf.slices.moments(METRIC, "throat")
    lo, hi = moment.reach(STATIC, "r")
    rm = np.linspace(hi, lo, 300)
    a_, b_ = pq(cf, "I'", 0 * rm, rm), pq(cf, "I", 0 * rm, rm[::-1])
    line = (np.concatenate([a_[0], b_[0]]), np.concatenate([a_[1], b_[1]]))
    for view in views:
        view.set(settings="$l = 1$ and $B = l/2$, so $r_+ = l$.")
        view.slice(moment, [line])
    return views


# ---------------------------------------------------------------- curvature figure

def register_curvature(cv):
    cv.DRAWN[METRIC] = lambda ck, src: curvature_views(ck, src, cv)


def curvature_views(ck, src, cv):
    """The stored Kretschmann scalar of the static chart at B = l/2 against r, from r = l/4 to 3l,
    read from the metric file and never recomputed; the Ricci scalar, -6/l^2 at every point, is
    stated in the caption. Checks hold K to Martinez and Zanelli's 12(r^6 + 2B^6)/(l^4 r^6), to
    falling outward toward 12/l^4, and the stored R to -6/l^2."""
    lo, hi = 0.25, TOP
    xaxis = cv.axis("$r/l$", lo, hi)
    scale = {"kind": "log", "from": 10, "to": 1e4}
    panels = [{"id": "K", "name": "$K\\,l^4$", "aspect": 2, "scale": scale, "ticks": cv.value_ticks(scale)}]
    cv.shape_checks(ck, METRIC, [p["aspect"] for p in panels])
    K = cv.nr_published(src, METRIC, STATIC, "kretschmann", ["r"], list(PARAMS))
    R = cv.nr_published(src, METRIC, STATIC, "ricci_scalar", ["r"], list(PARAMS))
    for at in (0.3, 1.0, 2.5):
        want = 12 * (at ** 6 + 2 * B ** 6) / (L ** 4 * at ** 6)
        ck(f"{METRIC} K at r = {at}l is Martinez and Zanelli's (KR)", abs(float(K.exact([at, 1, 0.5])) - want) < 1e-9 * want)
        ck(f"{METRIC} R at r = {at}l is -6/l^2", abs(float(R.exact([at, 1, 0.5])) + 6) < 1e-12)
    xs = np.linspace(lo, hi, 400)
    values = K(xs, 1, 0.5)
    ck(f"{METRIC} K lies within its scale and falls outward",
       bool(np.all(values > scale["from"]) and np.all(values < scale["to"]) and np.all(np.diff(values) < 0)))
    ck(f"{METRIC} K tends to 12/l^4 far out", abs(float(K(np.array([40.0]), 1, 0.5)[0]) - 12) < 1e-6)
    curves = cv.line_curves(ck, METRIC, "K", "K", K, "r", xaxis, scale, cv.PLOT_W / 2, {"l": 1.0, "B": 0.5})
    marks = [{"class": "horizon", "x": RP}]
    return [{"id": "hole", "label": "Dressed black hole", "kind": "line", "system": STATIC, "x": xaxis,
             "panels": panels, "curves": curves, "marks": marks,
             "legend": [["line", "K", "$K$"], ["line", "horizon", "the horizon $r_+ = 2B$"]],
             "caption": ["The Kretschmann scalar $K$ of the dressed black hole ($B = l/2$) against $r$, finite on "
                         "the horizon, growing without bound toward the singularity $r = 0$ and falling to $12/l^4$ "
                         "far away. The Ricci scalar is $-6/l^2$ at every point."],
             "settings": SETTINGS, "clip": "$K\\,l^4$ logarithmic, from $r = l/4$."}]
