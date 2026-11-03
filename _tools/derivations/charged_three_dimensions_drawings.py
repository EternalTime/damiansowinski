"""Drawings of the point charge in three dimensions, at kQ^2 = pi with r_c the unit of length.

There f = (kQ^2/2pi) ln(r_c/r) = ln(r_c/r)/2, written a ln(r_c/r) with a = 1/2, and the horizon r = r_c has
the surface gravity kappa = |f'(r_c)|/2 = a/(2 r_c) = 1/4. The tortoise coordinate that vanishes at r = 0,

    r_* = E_1(ln(r_c/r))/a        inside the horizon,      r_* = -Ei(ln(r/r_c))/a        outside it,

has dr_*/dr = 1/f on each side, tends to +infinity at r_c from both, stays finite at r = 0 and runs to
-infinity as r -> infinity, as ln(r)/r -> 0 and Ei(y) ~ e^y/y. With u = t - r_* and v = t + r_*, Kruskal's
coordinates are U = +-exp(kappa u) and V = +-exp(-kappa v), so that UV = F(r) with

    F(r) = -exp(-E_1(ln(r_c/r)))   inside,      F(r) = exp(Ei(ln(r/r_c)))   outside,      F(r_c) = 0,

which is e^gamma ln(r/r_c) to first order on both sides of the horizon, gamma Euler's constant, so the map is
smooth across it. The cells, with p = arctan U and q = arctan V, T = p + q up and X = q - p across:

    I    the static region on the left,    U > 0, V < 0, future along +t;
    II   the expanding region above,      U > 0, V > 0, future toward larger r;
    III  the static region on the right,   U < 0, V > 0, future along -t;
    IV   the contracting region below,    U < 0, V < 0, future toward smaller r.

r = 0 is UV = -1, the vertical lines X = -+pi/2, and r -> infinity is U or V -> infinity, the null lines that
close the hexagon above and below. The outgoing chart has U = exp(kappa u), V = F(r) exp(-kappa u) and covers
I and II; the ingoing chart has V = -exp(-kappa v), U = -F(r) exp(kappa v) and covers I and IV.

The slice t = 0 has g_rr = 1/f = 2/ln(r_c/r), at least 1 from r = e^{-2} r_c, where g_rr = 1 and the surface
lies level, out to the horizon, where it is the widest circle, so it is embedded in flat space between those
radii on each of the two static regions, joined at r_c into one closed surface.
"""
import math

import numpy as np
from scipy import special

METRIC = "charged_three_dimensions"
STATIC, OUT, IN = "static", "eddington_finkelstein_outgoing", "eddington_finkelstein_ingoing"
PARAMS = {"k": "pi", "Q": "1", "r_c": "1"}
A = 0.5                         # kQ^2/2pi
KAPPA = A / 2                   # the surface gravity of the horizon, in units of 1/r_c
LEVEL = math.exp(-1 / A)        # g_rr = 1 at r = r_c e^{-2pi/kQ^2} = e^{-2} r_c
SETTINGS = "$kQ^2 = \\pi$, and $r_c = 1$, the unit of every length."
FIXED = {"phi": "0"}


def rstar(r):
    """The tortoise coordinate that vanishes at r = 0, in units of r_c."""
    r = np.asarray(r, dtype=float)
    out = np.full(r.shape, np.nan)
    inside, outside = (r > 0) & (r < 1), r > 1
    out[inside] = special.exp1(np.log(1 / r[inside])) / A
    out[outside] = -special.expi(np.log(r[outside])) / A
    out[r == 0] = 0.0
    return out


def product(r):
    """UV = F(r), negative in the static regions and positive beyond the horizon."""
    r = np.asarray(r, dtype=float)
    out = np.zeros(r.shape)
    inside, outside = (r > 0) & (r < 1), r > 1
    out[inside] = -np.exp(-special.exp1(np.log(1 / r[inside])))
    out[outside] = np.exp(special.expi(np.log(r[outside])))
    out[r == 0] = -1.0
    return out


SIGNS = {"I": (1, -1), "II": (1, 1), "III": (-1, 1), "IV": (-1, -1)}


def kruskal(cf, cell, t, r):
    """(p, q) of the points (t, r) of a cell, t the static coordinate there."""
    t, rs = np.asarray(t, dtype=float), rstar(r)
    su, sv = SIGNS[cell]
    return su * cf.atan_exp(KAPPA * (t - rs)), sv * cf.atan_exp(-KAPPA * (t + rs))


def outgoing(cf, u, r):
    u = np.asarray(u, dtype=float)
    return cf.atan_exp(KAPPA * u), np.arctan(product(r) * np.exp(-KAPPA * u))


def ingoing(cf, v, r):
    v = np.asarray(v, dtype=float)
    return np.arctan(-product(r) * np.exp(KAPPA * v)), -cf.atan_exp(-KAPPA * v)


# -- Spacetime diagrams ----------------------------------------------------------------------------

def register_null(nr):
    rows = (
        nr.Diagram(METRIC, STATIC, "radial", "$t$ and $r$", ("t", "r"), (0, 2, -1, 1), "$r/r_c$", "$ct/r_c$",
                   PARAMS, FIXED, orient="outgoing"),
        nr.Diagram(METRIC, OUT, "radial", "against $u + r$", ("u", "r"), (0, 2, -1, 1), "$r/r_c$",
                   "$(cu + r)/r_c$", PARAMS, FIXED, to_display=nr.FINKELSTEIN_OUT, orient="outgoing"),
        nr.Diagram(METRIC, IN, "radial", "against $v - r$", ("v", "r"), (0, 2, -1, 1), "$r/r_c$",
                   "$(cv - r)/r_c$", PARAMS, FIXED, to_display=nr.FINKELSTEIN_IN, orient="ingoing"),
    )
    nr.DIAGRAMS.extend(rows)
    nr.CAPTIONS[(METRIC, STATIC, "radial")] = [
        "The plane of $t$ and $r$ ($\\phi = 0$), the same at every angle, with the charge at $r = 0$. "
        "The cones close at the horizon $r = r_c$ and open beyond it toward larger $r$, which is a time there."]
    nr.CAPTIONS[(METRIC, OUT, "radial")] = [
        "The plane of $u$ and $r$ ($\\phi = 0$) against $u + r$, the outgoing rays at 45°, crossing the horizon "
        "$r = r_c$ into the expanding region, where every cone opens toward larger $r$."]
    nr.CAPTIONS[(METRIC, IN, "radial")] = [
        "The plane of $v$ and $r$ ($\\phi = 0$) against $v - r$, the ingoing rays at 45°, crossing the horizon "
        "$r = r_c$ out of the contracting region, where every cone opens toward smaller $r$."]

    def away(w, r):
        return np.abs(np.asarray(r, dtype=float) - 1) > 0.05

    nr.CLOSED_FORMS[(METRIC, STATIC, "radial")] = (
        lambda t, r: t + rstar(r), lambda t, r: t - rstar(r), away)
    nr.CLOSED_FORMS[(METRIC, OUT, "radial")] = (
        lambda u, r: u + 2 * rstar(r), lambda u, r: u, away)
    nr.CLOSED_FORMS[(METRIC, IN, "radial")] = (
        lambda v, r: v, lambda v, r: v - 2 * rstar(r), away)


# -- The slice t = 0 on the spacetime diagrams ---------------------------------------------------------

def register_slices(sl):
    def drawn(sign):
        m = sl.moments(METRIC)[0]
        lo, hi = m.reach(STATIC, "r")
        r = sl.near(lo, hi)
        return [sl.Mark(m, [np.column_stack([sign * rstar(r), r])])]

    sl.FLAT[(METRIC, STATIC, "radial")] = lambda: sl.one(
        METRIC, lambda moment: sl.along(0.0, *moment.reach(STATIC, "r")))
    sl.FLAT[(METRIC, OUT, "radial")] = lambda: drawn(-1)
    sl.FLAT[(METRIC, IN, "radial")] = lambda: drawn(1)
    sl.FLAT_METRICS.add(METRIC)


# -- Embedding diagram ----------------------------------------------------------------------------------

def register_embedding(eb):
    eb.DRAWN[METRIC] = lambda ck, src: embedding(eb, ck, src)
    eb.CAPTIONS[(METRIC, "football")] = [
        "The slice $t = 0$ across both static regions, joined at the horizon $r = r_c$ into a closed surface, "
        "with $Q$ beyond one end and $-Q$ beyond the other."]


def embedding(eb, ck, src):
    """The two static regions at t = 0, each from r = e^{-2} r_c, where the surface lies level, up to the
    widest circle r_c, where the slice crosses the horizon into the other static region."""
    sl = eb.Slice(src, METRIC, STATIC, "r", "\\phi", {"t": 0, **FIXED}, PARAMS)
    ck.add("Point charge: g_rr = 1 at r = e^{-2} r_c", abs(float(sl.defect_at(np.array([LEVEL]))[0])), 1e-9)
    ck.stops("Point charge, nearer the charge than e^{-2} r_c", sl, np.linspace(0, LEVEL, 402)[1:-1])
    ck.stops("Point charge, beyond the horizon", sl, np.linspace(1, 3, 402)[1:-1])
    probe = eb.Piece("probe", "sheet", sl, LEVEL, 1.0, 0.0, 1)
    size = max(2.0, 2 * float(probe.z[-1] - probe.z[0]))
    lo = eb.Piece("static", "sheet", sl, LEVEL, 1.0, 0.0, 1,
                  (("stops", "at $r = e^{-2}r_c$ the surface lies level, and nearer the charge the circles grow "
                             "faster than the distance out to them"),
                   ("join", "the horizon $r = r_c$, the widest circle, where the slice runs on into the other "
                            "static region")),
                  [(LEVEL, "chartedge", None), (0.5, "r", None), (1.0, "horizon", "$r = r_c$")], size)
    lo.z = lo.z - lo.z[-1]
    hi = eb.Piece("other_static", "sheet2", sl, LEVEL, 1.0, -lo.z[0], -1,
                  (("stops", "at $r = e^{-2}r_c$"), ("join", "the horizon $r = r_c$")),
                  [(LEVEL, "chartedge", None), (0.5, "r2", None)], size)
    for p in (lo, hi):
        ck.isometry(f"Point charge, {p.id}", p)
    ck.join("Point charge, the two static regions at r_c", lo, 1.0, hi, 1.0)
    surface = eb.Surface([lo, hi])
    fig = eb.figure_of([surface], {"sheet": "cover"}, size, eb.Camera(-90, 22))
    eb.ring_label(fig, [0, 0, 0], 1.0, 0.0, "$r = r_c$", dx=10)
    eb.ring_label(fig, [0, 0, 0], *lo.at(LEVEL), "$e^{-2}r_c$", side=-1)
    fig.legend("fill", "cover", "the static region that $t$ and $r$ cover")
    fig.legend("line", "r", "$r = 0.5\\,r_c$")
    fig.legend("line", "r2", "the same circle on the other side")
    fig.legend("line", "horizon", "the horizon $r = r_c$, the widest circle")
    fig.legend("line", "chartedge", "$r = e^{-2}r_c$, where $g_{rr} = 1$")
    fig.legend("line", "meridian", "$\\phi$ constant, every $15°$")
    stops = ["Nearer the charge than $r = e^{-2}r_c$, $g_{rr} < 1$: the circles grow faster than the distance "
             "out to them, and no surface in flat space carries that part of the slice.",
             "Beyond the horizon, $r > r_c$, $g_{rr} < 0$: $r$ is a time there, and a slice of constant $t$ is "
             "not a moment of space."]
    return [eb.view("football", "The slice $t = 0$", "$r_c$", [surface], fig.done(), settings=SETTINGS,
                    stops=stops)]


# -- Conformal diagram ----------------------------------------------------------------------------------

HEXAGON = [[-math.pi / 2, -math.pi / 2], [-math.pi / 2, math.pi / 2], [0.0, math.pi],
           [math.pi / 2, math.pi / 2], [math.pi / 2, -math.pi / 2], [0.0, -math.pi]]
CELLS = {"I": [[-math.pi / 2, -math.pi / 2], [-math.pi / 2, math.pi / 2], [0.0, 0.0]],
         "II": [[-math.pi / 2, math.pi / 2], [0.0, math.pi], [math.pi / 2, math.pi / 2], [0.0, 0.0]],
         "III": [[math.pi / 2, math.pi / 2], [math.pi / 2, -math.pi / 2], [0.0, 0.0]],
         "IV": [[math.pi / 2, -math.pi / 2], [0.0, -math.pi], [-math.pi / 2, -math.pi / 2], [0.0, 0.0]]}
INSIDE, BEYOND, TIMES, NULLS = (0.2, 0.5, 0.8), (1.5, 3.0, 10.0), (-8, -4, 0, 4, 8), (-8, -4, 0, 4, 8)


def register_conformal(cf):
    cf.DRAWN[METRIC] = lambda ck, src: conformal(cf, ck, src)
    cf.CAPTIONS[(METRIC, "static")] = [
        "The whole spacetime ($\\phi = 0$) in the static chart, each point a circle: every line of constant $ct$ "
        "leaves the crossing of the horizons, and beyond $r_c$ each line of constant $r$ is a moment of space."]
    cf.CAPTIONS[(METRIC, "outgoing")] = [
        "The whole spacetime ($\\phi = 0$) in the outgoing chart, each point a circle: each outgoing ray of "
        "constant $u$ leaves the charge, crosses the future horizon and reaches null infinity."]
    cf.CAPTIONS[(METRIC, "ingoing")] = [
        "The whole spacetime ($\\phi = 0$) in the ingoing chart, each point a circle: each ingoing ray of "
        "constant $v$ comes from null infinity, crosses the past horizon and ends at the charge."]


def conformal(cf, ck, src):
    """The hexagon of the module docstring, one view for each chart, each tinting what its chart covers."""
    st = cf.Plane(src, METRIC, STATIC, ("t", "r"), FIXED, PARAMS)
    out = cf.Plane(src, METRIC, OUT, ("u", "r"), FIXED, PARAMS)
    inn = cf.Plane(src, METRIC, IN, ("v", "r"), FIXED, PARAMS)
    for cell, rs, future, name in (("I", (0.01, 0.99), (1, 0), "static region"),
                                   ("II", (1.01, 30), (0, 1), "expanding region"),
                                   ("III", (0.01, 0.99), (-1, 0), "other static region"),
                                   ("IV", (1.01, 30), (0, -1), "contracting region")):
        ck.chart(f"point charge, {name}", st, lambda t, r, c=cell: kruskal(cf, c, t, r),
                 ck.uniform(-30, 30), ck.uniform(*rs), lambda t, r, f=future: f)
    ck.chart("point charge, outgoing Eddington-Finkelstein", out, lambda u, r: outgoing(cf, u, r),
             ck.uniform(-30, 30), ck.uniform(0.01, 30), lambda u, r: (1, 60))
    ck.chart("point charge, ingoing Eddington-Finkelstein", inn, lambda v, r: ingoing(cf, v, r),
             ck.uniform(-30, 30), ck.uniform(0.01, 30), lambda v, r: (1, -60))

    rr = np.array([0.05, 0.3, 0.7, 0.95, 1.05, 2.0, 5.0])
    for cell, sel in (("I", rr < 1), ("II", rr > 1), ("III", rr < 1), ("IV", rr > 1)):
        p, q = kruskal(cf, cell, 0.7 + 0 * rr[sel], rr[sel])
        ck.limit(f"point charge {cell}: tan p tan q is UV = F(r)", np.tan(p) * np.tan(q), product(rr[sel]), 1e-9)
    for cell, side in (("I", -1), ("III", 1)):
        p, q = kruskal(cf, cell, np.array([-6.0, 0.0, 6.0]), np.full(3, 1e-300))
        ck.limit(f"point charge {cell}: r -> 0 lands on X = {side} pi/2", q - p, [side * math.pi / 2] * 3, 1e-9)
    p, q = outgoing(cf, np.array([-3.0, 0.0, 3.0]), np.full(3, 1e12))
    ck.limit("point charge: r -> infinity along an outgoing ray lands on q = pi/2", q, [math.pi / 2] * 3, 1e-6)
    p, q = ingoing(cf, np.array([-3.0, 0.0, 3.0]), np.full(3, 1e12))
    ck.limit("point charge: r -> infinity along an ingoing ray lands on p = -pi/2", p, [-math.pi / 2] * 3, 1e-6)
    p, q = kruskal(cf, "I", np.array([3.0]), np.array([1 - 1e-12]))
    ck.limit("point charge: r -> r_c at fixed t lands on the bifurcation circle", cf.point(p[0], q[0]), [0, 0], 1e-3)
    ck.limit("point charge: the outgoing and static coordinates put one event at one point",
             outgoing(cf, 2.0 - rstar(0.4), 0.4), kruskal(cf, "I", 2.0, 0.4), 1e-12)
    ck.limit("point charge: the ingoing and static coordinates put one event at one point",
             ingoing(cf, 2.0 + rstar(0.4), 0.4), kruskal(cf, "I", 2.0, 0.4), 1e-12)
    K = st.kretschmann
    ck.diverges("point charge: the Kretschmann scalar diverges at r = 0", K(0, 1e-2), K(0, 1e-3))
    ck.finite("point charge: the Kretschmann scalar is finite at r = r_c",
              K(np.zeros(3), np.array([0.999, 1.0, 1.001])))

    moment = cf.slices.moments(METRIC)[0]
    lo, hi = moment.reach(STATIC, "r")
    r_m = np.linspace(lo, hi - 1e-9, 400)
    left, right = kruskal(cf, "I", 0 * r_m, r_m), kruskal(cf, "III", 0 * r_m, r_m[::-1])
    through = (np.concatenate([left[0], right[0]]), np.concatenate([left[1], right[1]]))

    box = [-math.pi / 2 - 0.9, math.pi / 2 + 0.9, -math.pi - 0.25, math.pi + 0.25]
    t_all = cf.S_ALL
    inside = cf.spread(0, 1, 600, 16)
    beyond = cf.spread(1, np.inf, 600, 14)
    radii = cf.spread(0, np.inf, 800, 14)

    def edges(v):
        v.line("singular", [[[-math.pi / 2, -math.pi / 2], [-math.pi / 2, math.pi / 2]],
                            [[math.pi / 2, -math.pi / 2], [math.pi / 2, math.pi / 2]]], zig=True)
        v.line("horizon", [[[-math.pi / 2, -math.pi / 2], [math.pi / 2, math.pi / 2]],
                           [[math.pi / 2, -math.pi / 2], [-math.pi / 2, math.pi / 2]]])
        v.line("scri", [[[-math.pi / 2, math.pi / 2], [0.0, math.pi], [math.pi / 2, math.pi / 2]],
                        [[-math.pi / 2, -math.pi / 2], [0.0, -math.pi], [math.pi / 2, -math.pi / 2]]])
        v.label_xt([-math.pi / 2, 0.9], "$Q$", "r", dx=-8)
        v.label_xt([math.pi / 2, 0.9], "$-Q$", "l", dx=8)
        v.label_xt([-math.pi / 2, -0.9], "$r = 0$", "r", "small", dx=-8)
        v.label_xt([-math.pi / 4, math.pi / 4], "$r = r_c$", "br", "small", dx=-4, dy=-2)
        v.label_xt([-math.pi / 4, 3 * math.pi / 4], "$r \\to \\infty$", "br", "small", dx=-4, dy=-3)
        v.label_xt([-1.0, 0.0], "static", cls="region")
        v.label_xt([1.0, 0.0], "static", cls="region")
        v.label_xt([0.0, 1.7], "expanding", cls="region")
        v.label_xt([0.0, -1.7], "contracting", cls="region")
        v.legend("horizon", "the horizon $r = r_c$")
        v.legend("singular", "$r = 0$, a curvature singularity")
        v.legend("scri", "null infinity, $r \\to \\infty$")

    def static_grid(v, cells):
        for cell in cells:
            for r in (INSIDE if cell in ("I", "III") else BEYOND):
                v.curve("r", *kruskal(cf, cell, t_all, np.full_like(t_all, r)))
            span = inside if cell in ("I", "III") else beyond
            for t in TIMES:
                v.curve("t", *kruskal(cf, cell, np.full_like(span, t), span))

    views = []
    v = cf.View("static", "Static", box, STATIC)
    v.fill("region", HEXAGON)
    v.fill("cover", CELLS["I"])
    v.fill("cover", CELLS["II"])
    static_grid(v, ("I", "II"))
    edges(v)
    v.legend("cover", "what $t$ and $r$ cover")
    v.legend("r", "$r$ constant")
    v.legend("t", "$ct$ constant")
    v.slice(moment, [through])
    v.set(settings=SETTINGS)
    views.append(v)

    for vid, label, system, cells, fmap, w_name, name in (
            ("outgoing", "Outgoing Eddington-Finkelstein", OUT, ("I", "II"), outgoing, "u",
             "$u$ constant, an outgoing light ray"),
            ("ingoing", "Ingoing Eddington-Finkelstein", IN, ("I", "IV"), ingoing, "v",
             "$v$ constant, an ingoing light ray")):
        v = cf.View(vid, label, box, system)
        v.fill("region", HEXAGON)
        for cell in cells:
            v.fill("cover", CELLS[cell])
        for r in INSIDE + BEYOND:
            v.curve("r", *fmap(cf, t_all, np.full_like(t_all, r)))
        for w in NULLS:
            v.curve("null", *fmap(cf, np.full_like(radii, w), radii))
        edges(v)
        v.legend("cover", f"what ${w_name}$ and $r$ cover")
        v.legend("r", "$r$ constant")
        v.legend("null", name)
        v.slice(moment, [through])
        v.set(settings=SETTINGS)
        views.append(v)
    return views


# -- Curvature figure -----------------------------------------------------------------------------------

def register_curvature(cv):
    cv.DRAWN[METRIC] = lambda ck, src: curvature(cv, ck, src)


def curvature(cv, ck, src):
    """K and R against r from the stored scalars of the static chart, which every chart shares: K = 3R^2,
    both falling as powers of r, smooth through the horizon r = r_c, which is marked."""
    params = {"k": math.pi, "Q": 1.0, "r_c": 1.0}
    names = list(params)
    K = cv.Published(src, METRIC, STATIC, "kretschmann", ["r"], names)
    R = cv.Published(src, METRIC, STATIC, "ricci_scalar", ["r"], names)
    x_axis = cv.axis("$r/r_c$", 0.1, 3)
    K_scale = {"kind": "log", "from": 1e-3, "to": 1e4}
    R_scale = {"kind": "log", "from": 1e-2, "to": 1e2}
    panels = [{"id": "K", "name": "$K\\,r_c^4$", "aspect": 2, "scale": K_scale, "ticks": cv.value_ticks(K_scale)},
              {"id": "R", "name": "$R\\,r_c^2$", "aspect": 2, "scale": R_scale, "ticks": cv.value_ticks(R_scale)}]
    cv.shape_checks(ck, METRIC, [p["aspect"] for p in panels])
    height = {p["id"]: cv.PLOT_W / p["aspect"] for p in panels}
    curves = (cv.line_curves(ck, METRIC, "K", "K", K, "r", x_axis, K_scale, height["K"], params)
              + cv.line_curves(ck, METRIC, "R", "R", R, "r", x_axis, R_scale, height["R"], params))
    vals = list(params.values())
    rs = np.linspace(0.1, 3, cv.SAMPLES)
    ck("Point charge: K = 3 R^2", bool(np.allclose(K(rs, *vals), 3 * R(rs, *vals) ** 2, rtol=1e-12)))
    ck("Point charge: R = kQ^2/(2 pi r^2)", bool(np.allclose(R(rs, *vals), 1 / (2 * rs ** 2), rtol=1e-12)))
    ck("Point charge: K and R fall outward", bool(np.all(np.diff(K(rs, *vals)) < 0) and np.all(np.diff(R(rs, *vals)) < 0)))
    for radius in (0.25, 1.0, 2.0):
        got = float(K.exact([radius] + vals))
        ck(f"Point charge: K = 3 (kQ^2)^2/(4 pi^2 r^4) at r = {radius}", abs(got / (3 / (4 * radius ** 4)) - 1) < 1e-12)
    legend = [["line", "K", "$K$"], ["line", "R", "$R$"], ["line", "horizon", "the horizon $r = r_c$"]]
    return [{
        "id": "radial", "label": "Along $r$", "kind": "line", "x": x_axis,
        "panels": panels, "curves": curves, "marks": [{"class": "horizon", "x": 1.0}], "legend": legend,
        "caption": ["The Kretschmann scalar $K$ and the Ricci scalar $R$ of the point charge against $r$ "
                    "($kQ^2 = \\pi$), the same at every $t$ and $\\phi$ and smooth through the horizon. "
                    "$K = 3R^2$, growing as $r^{-4}$ toward the charge. $R$ grows as $r^{-2}$ toward the charge."],
        "ricci_sentences": [[0, 2]],
        "per_scalar": cv.scalar_clips("A logarithmic scale from $10^{-3}$ to $10^4$.",
                                      "A logarithmic scale from $10^{-2}$ to $10^2$."),
        "settings": SETTINGS,
        "clip": "Both on logarithmic scales; every value drawn lies within them.",
    }]
