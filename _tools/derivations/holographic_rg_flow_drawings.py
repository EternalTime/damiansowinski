"""The drawings of the flows of the renormalisation group of holography, in units L = 1 of the radius of
the anti-de Sitter space each approaches at its boundary, every x_i but the one rolled up held at 0.

The confining flow, e^{2A} = e^{2r/L} - 1, has radial light rays c dt = +-e^{-A} dr, and
z = int_r^infinity e^{-A} dr' = L arcsin(e^{-r/L}) is finite at both ends: the plane of t and r is
e^{2A}(-c^2 dt^2 + dz^2), conformal to the strip 0 < z < pi L/2 as it stands, with the boundary r -> infinity
at z = 0 and the singularity r = 0 at z = pi L/2, both timelike. It is drawn as X = pi/2 - z/L =
arccos(e^{-r/L}) and T = ct/L, so p, q = (T -+ X)/2, the singularity on the left at X = 0 and the
boundary on the right at X = pi/2. The chart of u = 1 - e^{-2r/L} covers the same plane and is drawn on
its spacetime diagram alone.

The disc of branes, e^{2A} = (ell^2/L^2) v^{2/3}/(1 - v), dv/dr = (2/L) v^{2/3}(1 - v), has
dz = e^{-A} dr = (L^2/2 ell) dv/(v sqrt(1 - v)), so z = (L^2/ell) artanh(sqrt(1 - v)) runs from 0 at the
boundary v = 1 to infinity at v = 0: the plane is conformal to a half of Minkowski's plane, as Freedman,
Gubser, Pilch, and Warner found for n = 2 [hep-th/9906194, after (zMetric)]. With w = -z and the
unit lambda = L^2/ell, p, q = arctan((ct -+ w)/lambda) draw it as the left half of Minkowski's diamond,
the boundary on X = 0. Radial light reaches v = 0 at a finite affine parameter, int e^{A} dr ~ int v^{-1/3} dv,
and at infinite t, so the singularity is the two null edges of the half diamond and their corner i^0.

The embedding diagrams draw the surface of r, or v, and x_1 at t = 0 with the strip 0 <= x_1 < 2 pi L rolled
up so that x_1 runs once round, as Randall and Sundrum's page draws its wall. Its circles have the radius
rho = e^A L, and in the proper distance d rho/dr = rho A'. For the confining flow
d rho/dr = e^{2r/L}/sqrt(e^{2r/L} - 1) >= 2, least at e^{2r/L} = 2, so every circle grows faster than the
distance out to it and the surface lies in Minkowski space, near its light cone at both ends. For the
disc d rho/dr = (2 ell/3L)(1 + v/2)/sqrt(1 - v), which is at least 1 everywhere exactly when
ell >= 3L/2, and every drawing of the disc takes ell = 3L/2, where it is 1 at the tip v = 0 and the
surface leaves the axis level. The curvature figure is the stored Kretschmann and Ricci scalars of the
proper distance chart against r and of the disc against v. holographic_rg_flow.md records each choice.
"""
import math

import numpy as np

METRIC = "holographic_rg_flow"
PD, BU, CD = "proper_distance", "bfs_u", "coulomb_disc"
WALL = {"x_1": "0", "x_2": "0", "x_3": "0"}
FLOW = {"L": "1"}
DISC = {"L": "1", "ell": "3/2"}
ELL = 1.5                                   # the disc drawn, ell/L
LAMBDA = 1 / ELL                            # L^2/ell, the unit of z on the disc
FLOW_TOP = 1.5                              # how far out the confining flow is embedded, r/L
FLOW_LO = 0.05                              # how near the singularity it starts, r/L, where it nears the light cone
DISC_TOP = 0.9                              # how far out the disc is embedded, v
FLOW_SETTINGS = "$L = 1$."
DISC_SETTINGS = "$L = 1$, $\\ell = 3L/2$."


def flow_rho(r):
    """The radius e^A L of the circles of the confining flow at L = 1."""
    return np.sqrt(np.expm1(2 * np.asarray(r, dtype=float)))


def disc_rho(v):
    """The radius e^A L of the circles of the disc at L = 1, ell = 3/2."""
    v = np.asarray(v, dtype=float)
    return ELL * np.cbrt(v) / np.sqrt(1 - v)


def flow_x(r):
    """X = arccos(e^{-r/L}) of the confining flow's strip."""
    return np.arccos(np.exp(-np.asarray(r, dtype=float)))


def flow_pq(t, r):
    X, T = flow_x(r), np.asarray(t, dtype=float)
    return (T - X) / 2, (T + X) / 2


def disc_z(v):
    """z = (L^2/ell) artanh(sqrt(1 - v)), zero at the boundary v = 1, written as
    (L^2/ell)(ln(1 + sqrt(1 - v)) - ln(v)/2), since (1 + s)/(1 - s) = (1 + s)^2/v for s = sqrt(1 - v),
    which keeps every digit where 1 - v rounds to 1."""
    v = np.asarray(v, dtype=float)
    return LAMBDA * (np.log1p(np.sqrt(1 - v)) - 0.5 * np.log(v))


def disc_pq(t, v):
    w = -disc_z(v)
    t = np.asarray(t, dtype=float)
    return np.arctan((t - w) / LAMBDA), np.arctan((t + w) / LAMBDA)


# ---------------------------------------------------------------- spacetime diagrams

def register_null(nr):
    nr.DIAGRAMS.extend((
        nr.Diagram(METRIC, PD, "tr", "$t$ and $r$", ("t", "r"), (0, 2, -1, 1), "$r/L$", "$ct/L$",
                   FLOW, WALL, families=nr.SIDEWAYS, cones=(6, 7)),
        nr.Diagram(METRIC, BU, "tu", "$t$ and $u$", ("t", "u"), (0, 1, -0.8, 0.8), "$u$", "$ct/L$",
                   FLOW, WALL, families=nr.SIDEWAYS, cones=(6, 7)),
        nr.Diagram(METRIC, CD, "tv", "$t$ and $v$", ("t", "v"), (0, 1, -0.8, 0.8), "$v$", "$ct/L$",
                   DISC, WALL, families=nr.SIDEWAYS, cones=(6, 7), singular_zero="v"),
    ))
    nr.CAPTIONS[(METRIC, PD, "tr")] = [
        "The plane of $t$ and $r$ ($x_1 = x_2 = x_3 = 0$) of the confining flow, each point in the diagram a "
        "single event, with the cones closing toward the singularity $r = 0$ and opening toward the boundary as "
        "$\\sqrt{e^{2r/L} - 1}$."]
    nr.CAPTIONS[(METRIC, BU, "tu")] = [
        "The plane of $t$ and $u$ ($x_1 = x_2 = x_3 = 0$) of the confining flow, each point in the diagram a "
        "single event, with light rays $du/d(ct) = \\pm 2\\sqrt{u(1 - u)}/L$ and the cones closing at the "
        "singularity $u = 0$ and at the boundary $u = 1$."]
    nr.CAPTIONS[(METRIC, CD, "tv")] = [
        "The plane of $t$ and $v$ ($x_1 = x_2 = x_3 = 0$) of the disc of branes ($\\ell = 3L/2$), each point in "
        "the diagram a single event, with light rays $dv/d(ct) = \\pm 2\\ell v\\sqrt{1 - v}/L^2$ and the cones "
        "closing toward the singularity $v = 0$, which light reaches only as $t \\to \\pm\\infty$."]


def register_slices(sl):
    sl.FLAT[(METRIC, PD, "tr")] = lambda: sl.one(
        METRIC, lambda m: sl.along(0.0, *m.reach(PD, "r")), view_id="flow")
    sl.FLAT[(METRIC, BU, "tu")] = lambda: sl.one(
        METRIC, lambda m: sl.along(0.0, *(-math.expm1(-2 * r) for r in m.reach(PD, "r"))), view_id="flow")
    sl.FLAT[(METRIC, CD, "tv")] = lambda: sl.one(
        METRIC, lambda m: sl.along(0.0, *m.reach(CD, "v")), view_id="disc")
    sl.FLAT_METRICS.add(METRIC)


# ---------------------------------------------------------------- embedding diagrams

def register_embedding(eb):
    eb.DRAWN[METRIC] = lambda ck, src: embedding(eb, ck, src)
    eb.CAPTIONS[(METRIC, "flow")] = [
        "The surface of $r$ and $x_1$ at $t = 0$ of the confining flow in Minkowski space, from $r = L/20$ near "
        "the singularity out toward the boundary, close to a light cone throughout."]
    eb.CAPTIONS[(METRIC, "disc")] = [
        "The surface of $v$ and $x_1$ at $t = 0$ of the disc of branes ($\\ell = 3L/2$) in Minkowski space, level "
        "at the singularity $v = 0$ on its axis and rising toward the boundary."]


def embedding(eb, ck, src):
    """The surface of r and x_1 of the confining flow and of v and x_1 of the disc, each at t = 0,
    x_2 = x_3 = 0, the strip 0 <= x_1 < 2 pi L rolled up, in Minkowski space."""
    fixed_at = {"t": 0, "x_2": 0, "x_3": 0}
    rolled = "$L = 1$, the strip $0 \\le x_1 < 2\\pi L$ rolled up, every length measured with $dX^2 + dY^2 - dZ^2$."
    stop = ("Every circle grows faster than the distance out to it, so no surface of revolution in Euclidean space "
            "has this geometry; one in Minkowski space does.")
    views = []

    name = "confining flow"
    sl = eb.Slice(src, METRIC, PD, "r", "x_1", fixed_at, {"L": 1}, space="minkowski")
    ck.stops(f"{name}, the surface of r and x_1 in flat space", sl, np.linspace(1e-3, 4, 400))
    rs = np.linspace(1e-3, 6, 600)
    slope = np.exp(2 * rs) / np.sqrt(np.expm1(2 * rs))
    least = math.log(2) / 2
    ck.add(f"{name}: d rho/dr = e^(2r)/sqrt(e^(2r) - 1) is at least 2, and 2 at e^(2r) = 2",
           float(np.sum(slope < 2)) + abs(math.exp(2 * least) / math.sqrt(math.expm1(2 * least)) - 2), 1e-12)
    size = 2 * float(flow_rho(FLOW_TOP))
    rings = (0.25, 0.5, 1.0, FLOW_TOP)
    sheet = eb.Piece("sheet", "sheet", sl, FLOW_LO, FLOW_TOP, 0.0, 1,
                     (("edge", "the surface runs on to the singularity $r = 0$, nearing the light cone"),
                      ("edge", "the surface runs on toward the boundary, $r \\to \\infty$")),
                     [(r, "r", None) for r in rings], size, digits=1e-12,
                     knots=tuple(float(k) for k in np.linspace(FLOW_LO, FLOW_TOP, 721)[1:-1]))
    ck.isometry(f"{name}, the sheet", sheet)
    ck.radius(f"{name}: rho = L sqrt(e^(2r/L) - 1)", sheet, flow_rho, size)
    ck.add(f"{name}: from r = L/20 to r = 3L/2 is 29L/20", abs(sl.proper(FLOW_LO, FLOW_TOP) - (FLOW_TOP - FLOW_LO)), 1e-9)
    surface = eb.Surface([sheet])
    fig = eb.figure_of([surface], {"sheet": "cover"}, size)
    eb.ring_label(fig, [0, 0, 0], *sheet.at(0.5), "$r = L/2$")
    eb.ring_label(fig, [0, 0, 0], *sheet.at(FLOW_TOP), "$3L/2$")
    fig.legend("fill", "cover", "the surface of $r$ and $x_1$")
    fig.legend("line", "r", "$r$ constant, at $L/4$, $L/2$, $L$, and $3L/2$")
    fig.legend("line", "meridian", "$x_1$ constant, every $\\pi L/12$")
    views.append(eb.view("flow", "Confining flow", "$L$", [surface], fig.done(), settings=rolled, stops=[stop]))

    name = "disc of branes"
    sd = eb.Slice(src, METRIC, CD, "v", "x_1", fixed_at, {"L": 1, "ell": ELL}, space="minkowski")
    ck.stops(f"{name}, the surface of v and x_1 in flat space", sd, np.linspace(1e-3, 0.99, 400))
    vs = np.linspace(0.0, 0.999, 600)
    ck.add(f"{name}: d rho/dr = (1 + v/2)/sqrt(1 - v) at ell = 3L/2, 1 at the tip and greater beyond",
           abs(float((1 + vs[0] / 2) / math.sqrt(1 - vs[0])) - 1)
           + float(np.sum(np.diff((1 + vs / 2) / np.sqrt(1 - vs)) <= 0)), 1e-12)
    size = 2 * float(disc_rho(DISC_TOP))
    rings = (0.25, 0.5, 0.75, DISC_TOP)
    sheet = eb.Piece("sheet", "sheet", sd, 0.0, DISC_TOP, 0.0, 1,
                     (("axis", "the singularity $v = 0$, where the circles shrink to a point"),
                      ("edge", "the surface runs on toward the boundary, $v \\to 1$")),
                     [(v, "r", None) for v in rings], size)
    ck.isometry(f"{name}, the sheet", sheet)
    ck.radius(f"{name}: rho = ell v^(1/3)/sqrt(1 - v)", sheet, disc_rho, size)
    surface = eb.Surface([sheet])
    fig = eb.figure_of([surface], {"sheet": "cover"}, size)
    eb.ring_label(fig, [0, 0, 0], *sheet.at(0.5), "$v = 1/2$")
    eb.ring_label(fig, [0, 0, 0], *sheet.at(DISC_TOP), "$0.9$")
    fig.legend("fill", "cover", "the surface of $v$ and $x_1$")
    fig.legend("line", "r", "$v$ constant, at $1/4$, $1/2$, $3/4$, and $0.9$")
    fig.legend("line", "meridian", "$x_1$ constant, every $\\pi L/12$")
    views.append(eb.view("disc", "Disc of branes", "$L$", [surface], fig.done(),
                         settings=rolled.replace("$L = 1$,", "$L = 1$, $\\ell = 3L/2$,"), stops=[stop]))
    return views


# ---------------------------------------------------------------- conformal diagrams

def register_conformal(cf):
    cf.DRAWN[METRIC] = lambda ck, src: conformal(cf, ck, src)
    cf.CAPTIONS[(METRIC, PD)] = [
        "The plane of $t$ and $r$ of the confining flow, each point in the diagram a single event, a strip "
        "between the timelike boundary and the timelike singularity $r = 0$, which light crosses in the time "
        "$\\pi L/2c$."]
    cf.CAPTIONS[(METRIC, CD)] = [
        "The plane of $t$ and $v$ of the disc of branes, each point in the diagram a single event, "
        "half of Minkowski's plane, with the timelike boundary on its right and the singularity $v = 0$ on its "
        "two null edges."]


def conformal(cf, ck, src):
    HALF, PI = cf.HALF, cf.PI
    views = []

    # The confining flow, the strip 0 < X < pi/2.
    plane = cf.Plane(src, METRIC, PD, ("t", "r"), WALL, {"L": 1})
    name = "confining flow"
    ck.chart(name, plane, flow_pq, ck.uniform(-6, 6), ck.uniform(1e-3, 12), lambda t, r: (1, 0))
    p, q = flow_pq(np.array([-1.0, 0.0, 1.0]), np.full(3, 40.0))
    ck.limit(f"{name}: r -> infinity lands on the boundary X = pi/2", q - p, [HALF] * 3, 1e-9)
    p, q = flow_pq(np.array([-1.0, 0.0, 1.0]), np.full(3, 1e-14))
    ck.limit(f"{name}: r -> 0 lands on the singularity X = 0", q - p, [0.0] * 3, 1e-6)
    ck.diverges(f"{name}: the Kretschmann scalar diverges at r = 0", plane.kretschmann(0, 1e-4),
                plane.kretschmann(0, 1e-5))
    ck.finite(f"{name}: the Kretschmann scalar is finite toward the boundary",
              plane.kretschmann(np.zeros(3), np.array([1.0, 5.0, 20.0])))
    T0, T1 = -HALF - 0.1, HALF + 0.1
    v = cf.View(PD, "Confining flow", [-0.3, HALF + 0.3, T0, T1], PD)
    v.fill("region", [[0, T0], [HALF, T0], [HALF, T1], [0, T1]])
    v.fill("cover", [[0, T0], [HALF, T0], [HALF, T1], [0, T1]])
    for r in (0.25, 0.5, 1.0, 2.0):
        X = float(flow_x(r))
        v.line("r", [[[X, T0], [X, T1]]])
    for k in (-1.0, -0.5, 0.0, 0.5, 1.0):
        v.line("t", [[[0, k], [HALF, k]]])
    v.line("null", [[[HALF, -HALF], [0, 0]], [[0, 0], [HALF, HALF]]])
    v.line("boundary", [[[HALF, T0], [HALF, T1]]])
    v.line("singular", [[[0, T0], [0, T1]]], zig=True)
    v.label_xt([float(flow_x(1.0)), T0], "$r = L$", "t", "coord", dy=4)
    v.set(fade={"top": 0.7, "bottom": 0.7})
    v.legend("cover", "the plane of $t$ and $r$")
    v.legend("r", "$r$ constant, at $L/4$, $L/2$, $L$, and $2L$")
    v.legend("t", "$ct$ constant, every $L/2$")
    v.legend("null", "a radial light ray")
    v.legend("boundary", "the conformal boundary $r \\to \\infty$")
    v.legend("singular", "the singularity $r = 0$")
    moment = cf.slices.moments(METRIC, "flow")[0]
    lo, hi = moment.reach(PD, "r")
    r = np.linspace(lo, hi, 400)
    v.slice(moment, [flow_pq(0 * r, r)])
    v.set(settings=FLOW_SETTINGS)
    views.append(v)

    # The disc, the left half of Minkowski's diamond.
    plane = cf.Plane(src, METRIC, CD, ("t", "v"), WALL, {"L": 1, "ell": ELL})
    name = "disc of branes"
    ck.chart(name, plane, disc_pq, ck.uniform(-6, 6), ck.uniform(1e-3, 1 - 1e-3), lambda t, v: (1, 0))
    p, q = disc_pq(np.array([-1.0, 0.0, 1.0]), np.full(3, 1 - 1e-14))
    ck.limit(f"{name}: v -> 1 lands on the boundary X = 0", q - p, [0.0] * 3, 1e-6)
    # z grows only as -(L^2/2 ell) ln v, so at v = 10^-300 the point is still 2 lambda/z = 0.006 from i0.
    p, q = disc_pq(np.zeros(1), np.full(1, 1e-300))
    ck.limit(f"{name}: v -> 0 at fixed t lands on i0", cf.point(p[0], q[0]), [-PI, 0], 1e-2)
    # Light toward v = 0 keeps ct + w = 0, z = ct, where v = 4e^{-2z/lambda} to within a part in 1/v, and runs
    # to ct - w -> infinity, the edge p = pi/2, which it nears as pi/2 - lambda/2ct.
    t = np.array([100.0, 200.0])
    p, q = disc_pq(t, 4 * np.exp(-2 * t / LAMBDA))
    ck.limit(f"{name}: light toward v = 0 reaches the null edge p = pi/2", p, HALF - LAMBDA / (2 * t), 1e-6)
    ck.limit(f"{name}: light toward v = 0 keeps q = 0", q, [0.0, 0.0], 1e-9)
    ck.diverges(f"{name}: the Kretschmann scalar diverges at v = 0", plane.kretschmann(0, 1e-6),
                plane.kretschmann(0, 1e-9))
    ck.finite(f"{name}: the Kretschmann scalar is finite toward the boundary",
              plane.kretschmann(np.zeros(3), np.array([0.5, 0.9, 0.999])))
    box = [-PI - 0.3, 0.3, -PI - 0.1, PI + 0.1]
    v = cf.View(CD, "Disc of branes", box, CD)
    v.fill("region", [[0, -PI], [0, PI], [-PI, 0]])
    v.fill("cover", [[0, -PI], [0, PI], [-PI, 0]])
    tt = cf.S_ALL
    for vv in (0.25, 0.5, 0.75, 0.95):
        v.curve("r", *disc_pq(tt, np.full_like(tt, vv)))
    vv = cf.spread(0.0, 1.0, 600, 14)
    for k in (-2.0, -1.0, 0.0, 1.0, 2.0):
        v.curve("t", *disc_pq(np.full_like(vv, k * LAMBDA), vv))
    v.line("boundary", [[[0, -PI], [0, PI]]])
    v.line("singular", [[[0, PI], [-PI, 0]], [[-PI, 0], [0, -PI]]], zig=True)
    v.label_xt([-PI, 0], "$i^0$", "r", dx=-6)
    v.legend("cover", "the plane of $t$ and $v$")
    v.legend("r", "$v$ constant, at $1/4$, $1/2$, $3/4$, and $0.95$")
    v.legend("t", "$ct$ constant, every $L^2/\\ell$")
    v.legend("boundary", "the conformal boundary $v \\to 1$")
    v.legend("singular", "the singularity $v = 0$")
    moment = cf.slices.moments(METRIC, "disc")[0]
    lo, hi = moment.reach(CD, "v")
    vs = np.linspace(max(lo, 1e-12), hi, 400)
    v.slice(moment, [disc_pq(0 * vs, vs)])
    v.set(settings=DISC_SETTINGS)
    views.append(v)
    return views


# ---------------------------------------------------------------- curvature figure

def register_curvature(cv):
    cv.DRAWN[METRIC] = lambda ck, src: curvature(cv, ck, src)


def _view(cv, ck, src, system, coordinate, label, x_name, lo, hi, K_scale, R_scale, mark, caption, settings,
          params, checks):
    names = list(params)
    Kp = cv.Published(src, METRIC, system, "kretschmann", [coordinate], names)
    Rp = cv.Published(src, METRIC, system, "ricci_scalar", [coordinate], names)
    x_axis = cv.axis(x_name, lo, hi)
    panels = [{"id": "K", "name": "$K\\,L^4$", "aspect": 2, "scale": K_scale, "ticks": cv.value_ticks(K_scale)},
              {"id": "R", "name": "$\\mathcal{R}\\,L^2$", "aspect": 2, "scale": R_scale,
               "ticks": cv.value_ticks(R_scale)}]
    cv.shape_checks(ck, METRIC, [p["aspect"] for p in panels])
    height = {p["id"]: cv.PLOT_W / p["aspect"] for p in panels}
    curves = (cv.line_curves(ck, METRIC, "K", "K", Kp, coordinate, x_axis, K_scale, height["K"], params)
              + cv.line_curves(ck, METRIC, "R", "R", Rp, coordinate, x_axis, R_scale, height["R"], params))
    checks(Kp, Rp, list(params.values()))
    return {
        "id": system, "label": label, "kind": "line", "system": system, "x": x_axis,
        "panels": panels, "curves": curves, "marks": mark,
        "legend": [["line", "K", "$K$"], ["line", "R", "$\\mathcal{R}$"]],
        "caption": caption, "ricci_sentences": [[0, 2]],
        "per_scalar": cv.scalar_clips("A logarithmic scale from $10$ to $10^8$.",
                                      "A signed logarithmic scale from $-10^4$ to $-10$."),
        "settings": settings,
        "clip": "$K$ logarithmic, $\\mathcal{R}$ signed logarithmic.",
    }


def curvature(cv, ck, src):
    """K and the Ricci scalar from the stored scalars: of the confining flow against r, of the disc
    against v. Each tends to anti-de Sitter's, K = 40/L^4 and R = -20/L^2, at the boundary, each grows
    without bound toward the singularity, K is positive and R negative throughout."""
    K_scale = {"kind": "log", "from": 10, "to": 1e8}
    R_scale = {"kind": "symlog", "from": -1e4, "to": -10, "linear": 1}

    def held(name, x_far, x_near, x_all):
        def checks(Kp, Rp, vals):
            ck(f"{name}: K tends to 40/L^4 at the boundary",
               abs(float(Kp(np.array([x_far]), *vals)[0]) - 40) < 1e-6)
            ck(f"{name}: R tends to -20/L^2 at the boundary",
               abs(float(Rp(np.array([x_far]), *vals)[0]) + 20) < 1e-6)
            ck(f"{name}: K is positive and R negative", bool(np.all(Kp(x_all, *vals) > 0) and np.all(Rp(x_all, *vals) < 0)))
            ck(f"{name}: K grows without bound toward the singularity",
               float(Kp(np.array([x_near]), *vals)[0]) > 1e6)
        return checks

    rs = np.linspace(0.05, 3.0, 400)
    vs = np.linspace(1e-4, 0.99, 400)
    return [
        _view(cv, ck, src, PD, "r", "Confining flow", "$r/L$", 0.05, 3.0, K_scale, R_scale, [],
              ["The Kretschmann scalar $K$ and the Ricci scalar $\\mathcal{R}$ of the confining flow against $r$, "
               "diverging toward the singularity $r = 0$. $K$ tends to anti-de Sitter's $40/L^4$ at the "
               "boundary. $\\mathcal{R}$ is negative and tends to anti-de Sitter's $-20/L^2$ at the boundary."],
              FLOW_SETTINGS, {"L": 1.0}, held("confining flow", 30.0, 1e-3, rs)),
        _view(cv, ck, src, CD, "v", "Disc of branes", "$v$", 0.0, 1.0, K_scale, R_scale, [],
              ["The Kretschmann scalar $K$ and the Ricci scalar $\\mathcal{R}$ of the disc of branes against $v$, "
               "diverging toward the singularity $v = 0$. $K$ tends to anti-de Sitter's $40/L^4$ at the "
               "boundary. $\\mathcal{R}$ is negative and tends to anti-de Sitter's $-20/L^2$ at the boundary."],
              DISC_SETTINGS, {"L": 1.0, "ell": ELL}, held("disc of branes", 1 - 1e-12, 1e-8, vs)),
    ]
