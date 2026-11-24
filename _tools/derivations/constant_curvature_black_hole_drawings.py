"""The drawings of the black holes of constant curvature, anti-de Sitter space with its points identified
along a boost.

Every drawing takes l = 1, the unit of every length, and r_+ = l. constant_curvature_black_hole.md
beside this file is the derivation.

The spacetime diagrams. The Kruskal charts are conformally flat in (y_0, y_1, ..., y_{n-2}) with the factor
4l^2/(1 - y^2)^2, and the plane y_2 = 0 (and y_3 = 0 in five dimensions) at phi = 0 is the fixed set of two
reflections, so it is totally geodesic and light runs on it at 45 degrees. On it the event horizon
y_0^2 = y_1^2 + y_2^2 is the pair of rays y_0 = +-y_1, the singularity y_0^2 = 1 + y_1^2, where the circle of
phi has radius r = 0, and infinity y_1^2 - y_0^2 = 1. The static chart at theta = pi/2 and the exterior chart
at theta = 0 both draw the wedge y_1 > |y_0| of that plane, y_0 = f sinh(r_+ t/l), y_1 = f cosh(r_+ t/l),
f = sqrt((r - r_+)/(r + r_+)), with the rays dt = +-l dr/(r^2 - r_+^2). Holst and Peldan's interior is drawn
in its plane of r and t with -r up, since the future of the black hole runs toward the singularity r = 0,
and its null chart in (v - u)/2 across and (u + v)/2 up.

The embedding diagram. The moment t = 0 of the exterior chart through theta = 0 and theta = pi is the line
y_0 = 0, y_2 = 0 of the Kruskal chart times the circle of phi, totally geodesic as the fixed set of the
reflections t -> -t and y_2 -> -y_2. Its metric is l^2 dr^2/(r^2 - r_+^2) + r^2 dphi^2, the moment of the
spinless hole of three dimensions at M = r_+^2/l^2, so it is drawn as that hole's is: in flat space from the
throat r = r_+, where the moment meets the horizon in one circle, to r = sqrt(r_+^2 + l^2), where it lies
level, and in Minkowski space of three dimensions beyond, where its circles grow faster than the distance
out to them. The two sheets are the halves theta = 0 and theta = pi of the one exterior, which the angle
theta joins round the axis r = r_+.

The conformal diagram. In the Kruskal chart of four dimensions, with y_1 = rho cos(theta), y_2 = rho sin(theta),
the half plane of y_0 and rho >= 0 has the metric 4l^2(-dy_0^2 + drho^2)/(1 - rho^2 + y_0^2)^2, each point a
torus of theta and phi, a circle on the axis rho = 0. p = arctan(y_0 - rho) and q = arctan(y_0 + rho) bring it
to the rectangle 0 <= X < pi/2, |T| < pi/2, X = q - p and T = p + q: tan p tan q = y_0^2 - rho^2 puts infinity
on X = pi/2 and the singularities on T = +-pi/2, and the event horizon is the ray p = 0, X = T, from the
middle of the axis to the corner where infinity meets the future singularity. In five dimensions the same
rectangle holds, each point a 2-sphere times a circle.
"""
import math

import numpy as np

METRIC = "constant_curvature_black_hole"
PARAMS = {"l": 1, "r_+": 1}
SETTINGS = "$l = 1$, the unit of every length, and $r_+ = l$."
TOP = 3.0            # the embedding is drawn out to r = 3l
LEVEL = math.sqrt(2.0)


def f(r):
    """Banados, Gomberoff and Martinez's f = sqrt((r - r_+)/(r + r_+)) at r_+ = 1, the y_1 of the moment."""
    r = np.asarray(r, dtype=float)
    return np.sqrt((r - 1) / (r + 1))


# ---------------------------------------------------------------- the spacetime diagrams

KRUSKAL_WHERE = "(1 + y_0**2 - y_1**2)*(1 - y_0**2 + y_1**2)"
HORIZONS = (("event", {"x0": "0", "r": "0"}, "both", "the event horizon, $y_0 = \\pm y_1$"),)


def register_null(nr):
    nr.DIAGRAMS.extend([
        nr.Diagram(METRIC, "kruskal", "plane", "$y_0$ and $y_1$", ("y_0", "y_1"), (-1.5, 1.5, -1.5, 1.5),
                   "$y_1$", "$y_0$", PARAMS, {"y_2": "0", "varphi": "0"}, families=nr.SIDEWAYS, tau="y_0",
                   where=KRUSKAL_WHERE, marked=HORIZONS),
        nr.Diagram(METRIC, "static", "radial", "$t$ and $r$", ("t", "r"), (0, 3, -1.5, 1.5),
                   "$r/l$", "$t$", PARAMS, {"theta": "pi/2", "varphi": "0"}, orient="ingoing"),
        nr.Diagram(METRIC, "exterior", "radial", "$t$ and $r$", ("t", "r"), (0, 3, -1.5, 1.5),
                   "$r/l$", "$t$", PARAMS, {"theta": "0", "varphi": "0"}, orient="ingoing"),
        nr.Diagram(METRIC, "interior", "radial", "$r$ and $t$", ("r", "t"), (0, 1.5, -1, 0),
                   "$t$", "$-r$", PARAMS, {"psi": "0", "varphi": "0"}, to_display=((0, 1), (-1, 0)), tau="-r",
                   families=nr.SIDEWAYS,
                   lines=(("apparent", "x0", "1/2", "the apparent horizon $r = 1/2$"),)),
        nr.Diagram(METRIC, "null", "plane", "$u$ and $v$", ("u", "v"), (0, math.pi / 4, 0, math.pi / 4),
                   "$(v - u)/2$", "$(u + v)/2$", PARAMS, {"psi": "0", "varphi": "0"}, to_display=nr.NULL_TO_TR,
                   tau="u + v", families=nr.SIDEWAYS, where="u*(pi/2 - u - v)",
                   curves=(("apparent", "atan(1/(3*tan(u)))", "the apparent horizon $\\tan u\\tan v = 1/3$"),)),
        nr.Diagram(METRIC, "kruskal_five", "plane", "$y_0$ and $y_1$", ("y_0", "y_1"), (-1.5, 1.5, -1.5, 1.5),
                   "$y_1$", "$y_0$", PARAMS, {"y_2": "0", "y_3": "0", "varphi": "0"}, families=nr.SIDEWAYS,
                   tau="y_0", where=KRUSKAL_WHERE, marked=HORIZONS),
        nr.Diagram(METRIC, "static_five", "radial", "$t$ and $r$", ("t", "r"), (0, 3, -1.5, 1.5),
                   "$r/l$", "$t$", PARAMS, {"theta": "0", "chi": "0", "varphi": "0"}, orient="ingoing"),
        nr.Diagram(METRIC, "exterior_five", "radial", "$t$ and $r$", ("t", "r"), (0, 3, -1.5, 1.5),
                   "$r/l$", "$t$", PARAMS, {"theta": "0", "chi": "0", "varphi": "0"}, orient="ingoing"),
    ])
    kruskal = ("light runs at 45°. The event horizon is the pair of rays $y_0 = \\pm y_1$, the singularity "
               "the hyperbola $y_0^2 - y_1^2 = 1$, where $r = 0$, and infinity the "
               "hyperbola $y_1^2 - y_0^2 = 1$.")
    static = "The rays follow $dt = \\pm l\\,dr/(r^2 - r_+^2)$, and the cones close at the horizon $r = r_+$."
    nr.CAPTIONS.update({
        (METRIC, "kruskal", "plane"): [
            "The plane of $y_0$ and $y_1$ ($y_2 = 0$, $\\varphi = 0$, $r_+ = l$), on which " + kruskal],
        (METRIC, "static", "radial"): [
            "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\varphi = 0$, $r_+ = l$). " + static],
        (METRIC, "exterior", "radial"): [
            "The plane of $t$ and $r$ ($\\theta = 0$, $\\varphi = 0$, $r_+ = l$), also the static chart's plane at "
            "$\\theta = \\pi/2$. " + static],
        (METRIC, "interior", "radial"): [
            "The plane of $r$ and $t$ ($\\psi = 0$, $\\varphi = 0$) between the event horizon $r = 1$ and the "
            "singularity $r = 0$, where the circle of $\\psi$ is null. The rays follow "
            "$dt = \\pm dr/(1 - r^2)$, and the apparent horizon of the slices of constant $r$ lies at $r = 1/2$."],
        (METRIC, "null", "plane"): [
            "The plane of $u$ and $v$ ($\\psi = 0$, $\\varphi = 0$) of the interior, on which the lines of "
            "constant $u$ and of constant $v$ are light rays. The event horizon is $u = 0$, the singularity "
            "$u + v = \\pi/2$, and the apparent horizon $\\tan u\\tan v = 1/3$."],
        (METRIC, "kruskal_five", "plane"): [
            "The plane of $y_0$ and $y_1$ ($y_2 = y_3 = 0$, $\\varphi = 0$, $r_+ = l$), on which " + kruskal],
        (METRIC, "static_five", "radial"): [
            "The plane of $t$ and $r$ ($\\theta = 0$, $\\chi = 0$, $\\varphi = 0$, $r_+ = l$). " + static],
        (METRIC, "exterior_five", "radial"): [
            "The plane of $t$ and $r$ ($\\theta = 0$, $\\chi = 0$, $\\varphi = 0$, $r_+ = l$). " + static],
    })


# ---------------------------------------------------------------- slices: the moment t = 0

FIVE = "The hole of five dimensions is another spacetime than the hole of four, whose moment is embedded."
INSIDE = ("The interior charts cover the black hole between the event horizon and the singularity, which the "
          "moment t = 0 meets only in the circle r = r_+ at their edge.")
HIDDEN_FIVE = [(METRIC, "kruskal_five", "plane"), (METRIC, "static_five", "radial"), (METRIC, "exterior_five", "radial")]
HIDDEN_INSIDE = [(METRIC, "interior", "radial"), (METRIC, "null", "plane")]


def slice_marks(slices, system):
    """The moment t = 0 of the exterior chart: along t = 0 from the throat out on the static and exterior
    charts, and on the Kruskal chart along y_0 = 0 through both halves, y_1 = +-f(r)."""
    m = slices.moments(METRIC, "throat")[0]
    lo, hi = m.reach("exterior", "r")
    if system == "kruskal":
        a = float(f(hi))
        return [slices.Mark(m, [[(0.0, -a), (0.0, a)]])]
    return [slices.Mark(m, slices.along(0.0, lo, hi))]


def register_slices(slices):
    slices.FLAT.update({(METRIC, system, view): (lambda system=system: slice_marks(slices, system))
                        for system, view in (("kruskal", "plane"), ("static", "radial"), ("exterior", "radial"))})
    slices.FLAT_METRICS.add(METRIC)
    slices.HIDDEN.update({key: FIVE for key in HIDDEN_FIVE})
    slices.HIDDEN.update({key: INSIDE for key in HIDDEN_INSIDE})


# ---------------------------------------------------------------- the embedding diagram

def embedding_views(ck, src, em):
    name = "constant curvature black hole"
    sl = em.Slice(src, METRIC, "exterior", "r", "\\varphi", {"t": 0, "theta": 0}, PARAMS)
    msl = em.Slice(src, METRIC, "exterior", "r", "\\varphi", {"t": 0, "theta": 0}, PARAMS, space="minkowski")
    rp = sl.horizons()[0]
    level, top = LEVEL, TOP
    size = 2 * top
    ck.add(f"{name}: the horizon of the exterior chart is at r_+ = l", abs(rp - 1.0), 1e-12)
    ck.add(f"{name}: the surface lies level where l^2/(r^2 - r_+^2) = 1, at r = sqrt(r_+^2 + l^2)",
           abs(float(sl.defect_at(np.array([level]))[0])), 1e-12)
    ck.stops(f"{name}, beyond r = sqrt 2 l in flat space", sl, np.linspace(level, 40, 402)[1:])
    inward = np.linspace(rp, level, 402)[1:-1]
    ck.add(f"{name}: between r_+ and sqrt 2 l no surface in Minkowski space carries the slice, (drho/dr)^2 - g_rr < 0",
           float(max(0.0, np.max(-msl.defect_at(inward)))), 0.0)
    if not np.all(msl.defect_at(inward) > 0):
        ck.items[-1]["ok"] = False
    join = ("at $r = \\sqrt{2}\\,l$ the surface lies level, in Euclidean space nearer the throat and in Minkowski "
            "space beyond")
    edge = "the sheet runs on toward a light cone, to $r \\to \\infty$"
    near = em.Piece("half", "sheet", sl, rp, level, 0.0, 1,
                    (("throat", "the throat $r = r_+$, where the half $\\theta = \\pi$ begins"), ("join", join)),
                    [(rp, "horizon", "$r = r_+$"), (1.2, "r", None), (level, "space", None)], size)
    out = em.Piece("half_minkowski", "sheet", msl, level, top, near.at(level)[1], 1, (("join", join), ("edge", edge)),
                   [(2.0, "r", None), (top, "r", None)], size)
    far = em.Piece("other_half", "sheet2", sl, rp, level, 0.0, -1, (("throat", "the throat $r = r_+$"), ("join", join)),
                   [(1.2, "r2", None), (level, "space", None)], size)
    far_out = em.Piece("other_half_minkowski", "sheet2", msl, level, top, far.at(level)[1], -1,
                       (("join", join), ("edge", edge)), [(2.0, "r2", None), (top, "r2", None)], size)
    for p in (near, out, far, far_out):
        space = "in Minkowski space" if p.sl.lorentz else "in flat space"
        ck.isometry(f"{name}, {p.id} {space}", p)
        ck.radius(f"{name}, {p.id}, rho = r {space}", p, lambda r: r, size)
    ck.join(f"{name}, the two halves at the throat", near, rp, far, rp)
    for a, b in ((near, out), (far, far_out)):
        ck.join(f"{name}, {a.id} in flat space and {b.id} in Minkowski space at r = sqrt 2 l", a, level, b, level)
        pa, pb = a.data()["points"][-1], b.data()["points"][0]
        ck.add(f"{name}, {a.id} and {b.id} as written: one point at r = sqrt 2 l",
               max(abs(pa[0] - pb[0]), abs(pa[1] - pb[1]), abs(pa[2] - pb[2])), 10.0 ** -min(a.decimals, b.decimals))
    surface = em.Surface([near, out, far, far_out])
    fig = em.figure_of([surface], {"sheet": "cover"}, size)
    em.ring_label(fig, [0, 0, 0], rp, 0.0, "$r = r_+$", dx=14)
    em.ring_label(fig, [0, 0, 0], *near.at(level), "$\\sqrt{2}\\,l$")
    em.ring_label(fig, [0, 0, 0], *out.at(top), "$3l$")
    fig.legend("fill", "cover", "the half $\\theta = 0$")
    fig.legend("line", "r", "$r$ at $1.2\\,l$, $2l$ and $3l$")
    fig.legend("line", "r2", "the same at $\\theta = \\pi$")
    fig.legend("line", "horizon", "the throat $r = r_+$")
    fig.legend("line", "space", "$r = \\sqrt{2}\\,l$, Euclidean space inside and Minkowski space beyond")
    fig.legend("line", "meridian", "$\\varphi$ every $15°$")
    return [em.view("throat", "Through the throat", "$l$", [surface], fig.done(),
                    settings="$l = 1$ and $r_+ = l$, and $dX^2 + dY^2 - dZ^2$ in Minkowski space.")]


def embedding_captions():
    return {(METRIC, "throat"): [
        "The plane of $r$ and $\\varphi$ through $\\theta = 0$ and $\\pi$ at $t = 0$, its throat the "
        "horizon's circle $r = r_+$ and its two sheets halves of one exterior."]}


# ---------------------------------------------------------------- the conformal diagram

def conformal_views(ck, src, cf):
    name = "constant curvature black hole, the Kruskal half plane"
    plane = cf.Plane(src, METRIC, "kruskal", ("y_0", "y_1"), {"y_2": "0", "varphi": "0"}, PARAMS)

    def E(y0, rho):
        y0, rho = np.asarray(y0, dtype=float), np.asarray(rho, dtype=float)
        return np.arctan(y0 - rho), np.arctan(y0 + rho)

    rng = np.random.default_rng(7)
    # The exterior and both interiors, inside -1 < rho^2 - y_0^2 < 1.
    rho = rng.uniform(0.0, 0.95, 400)
    y0 = rng.uniform(-1, 1, 400) * np.sqrt(np.maximum(rho ** 2 + 0.95 - 0.0, 0.0))
    keep = (rho ** 2 - y0 ** 2 > -0.98) & (rho ** 2 - y0 ** 2 < 0.98)
    ck.chart(name, plane, E, y0[keep], rho[keep], lambda a, b: (1, 0))
    p, q = E(np.array([-0.5, 0.0, 0.7]), np.zeros(3))
    ck.limit(f"{name}: the axis rho = 0 is X = 0", q - p, np.zeros(3), 1e-12)
    s = np.array([0.0, 2.0, -3.0])
    p, q = E(s, np.sqrt(1 + s ** 2) * (1 - 1e-12))
    ck.limit(f"{name}: infinity rho^2 - y_0^2 = 1 lands on X = pi/2", q - p, [cf.HALF] * 3, 1e-6)
    p, q = E(np.sqrt(1 + s ** 2) * (1 - 1e-12), np.abs(s))
    ck.limit(f"{name}: the future singularity y_0^2 - rho^2 = 1 lands on T = pi/2", p + q, [cf.HALF] * 3, 1e-6)
    p, q = E(-np.sqrt(1 + s ** 2) * (1 - 1e-12), np.abs(s))
    ck.limit(f"{name}: the past singularity lands on T = -pi/2", p + q, [-cf.HALF] * 3, 1e-6)
    p, q = E(np.array([0.3, 0.6, 2.0]), np.array([0.3, 0.6, 2.0]))
    ck.limit(f"{name}: the event horizon y_0 = rho is the ray p = 0", p, np.zeros(3), 1e-12)
    ck.finite(f"{name}: the Kretschmann scalar is finite on the singularity, where the circle of phi shrinks",
              plane.kretschmann(np.array([1 - 1e-9, 1.2, 2.0]), np.array([0.0, math.sqrt(1.44 - 1) * (1 + 1e-9),
                                                                          math.sqrt(3) * (1 + 1e-9)])))
    H = cf.HALF
    box = [-0.3, H + 0.45, -H - 0.3, H + 0.3]
    v = cf.View("whole", "The whole spacetime", box)
    rectangle = [[0, -H], [H, -H], [H, H], [0, H]]
    v.fill("region", rectangle)
    v.fill("cover", rectangle)
    t = cf.spread(-np.inf, np.inf, 500, 10)
    s_up = cf.spread(0.0, np.inf, 500, 10)
    # rho^2 - y_0^2 = c is the surface r = r_+(1 + c)/(1 - c): r = 2r_+ and 4r_+ outside, at c = 1/3 and 3/5,
    # and inside, at c = -1/3, r = r_+/2, the apparent horizon of Holst and Peldan's tori of constant r.
    for c in (1 / 3, 3 / 5):
        v.curve("r", *E(math.sqrt(c) * np.sinh(t), math.sqrt(c) * np.cosh(t)))
    a = math.sqrt(1 / 3)
    v.curve("apparent", *E(a * np.cosh(s_up), a * np.sinh(s_up)))
    ck.limit(f"{name}: rho^2 - y_0^2 = -1/3, 1/3 and 3/5 are r = r_+/2, 2r_+ and 4r_+",
             [(1 + c) / (1 - c) for c in (-1 / 3, 1 / 3, 3 / 5)], [0.5, 2.0, 4.0], 1e-12)
    v.line("horizon", [[[0, 0], [H, H]], [[0, 0], [H, -H]]])
    v.line("boundary", [[[H, -H], [H, H]]])
    v.line("singular", [[[0, H], [H, H]], [[0, -H], [H, -H]]], zig=True)
    v.line("centre", [[[0, -H], [0, H]]])
    v.label_xt([H / 2, H], "$r = 0$", "b", dy=-8)
    v.label_xt([H / 2, -H], "$r = 0$", "t", dy=8)
    v.label_xt([H, 0], "$r \\to \\infty$", "l", "small", dx=6)
    v.label_xt([H * 0.62, 0], "exterior", cls="region")
    v.label_xt([H * 0.28, H * 0.72], "black hole", cls="region")
    v.label_xt([H * 0.28, -H * 0.72], "white hole", cls="region")
    v.legend("horizon", "the event horizon $y_0 = \\rho$")
    v.legend("boundary", "infinity")
    v.legend("singular", "$r = 0$, where the curvature is finite")
    v.legend("centre", "the axis $\\rho = 0$")
    v.legend("r", "$r$ at $2r_+$ and $4r_+$")
    v.legend("apparent", "$r = r_+/2$, the apparent horizon")
    moment = cf.slices.moments(METRIC, "throat")[0]
    lo, hi = moment.reach("exterior", "r")
    rr = np.linspace(0.0, float(f(hi)), 400)
    v.slice(moment, [E(np.zeros_like(rr), rr)])
    v.set(settings="$l = 1$ and $r_+ = l$.",
          restriction="The half plane of $y_0$ and $\\rho = \\sqrt{y_1^2 + y_2^2}$.")
    return [v]


def conformal_captions():
    return {(METRIC, "whole"): [
        "The whole spacetime, $p = \\arctan(y_0 - \\rho)$ and $q = \\arctan(y_0 + \\rho)$, each point in the diagram "
        "a torus of $\\theta$ and $\\varphi$ in four dimensions and a 2-sphere times a circle in five."]}
