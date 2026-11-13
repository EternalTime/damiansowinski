"""Drawings of Horne and Horowitz's charged black string, at M = 1, Q^2 = M^2/2 and k = 52/23.

There the event horizon is r = M = 1 and the inner horizon r = Q^2/M = 1/2. With s = sqrt(1 - Q^2/Mr) and
a = sqrt(M^2 - Q^2)/M, the tortoise coordinate of the plane of t and r, dr_*/dr = sqrt(k/8)/(r f sqrt(h)), is

    r_* = sqrt(k/8) [ln(|s - a|/(s + a))/a + ln((1 + s)/(1 - s))],      r > Q^2/M,

zero at the inner horizon, -infinity at the event horizon and +infinity as r -> infinity. Its log at r = M
has the coefficient sqrt(k/8)/a = 1/2 kappa, kappa = sqrt(2(M^2 - Q^2)/k)/M, and

    F(r) = ((s - a)/(s + a)) ((1 + s)/(1 - s))^a = sign(r - M) exp(2 kappa r_*)

is smooth through r = M. With u = t - r_*, v = t + r_*, U = -+exp(-kappa u), V = +-exp(kappa v) and
p = arctan U, q = arctan V, the plane of t and r is Kruskal and Szekeres's square, UV = -F(r), with the
inner horizon UV = 1, T = +-pi/2, across its top and bottom: there h = 0 and inside it the plane of t and r
is spacelike, so the square ends there.

With sigma = sqrt(M/r - 1) and b = sqrt(M^2 - Q^2)/Q, the plane of x and r inside the event horizon has

    x_* = sqrt(k/8) [pi - 2 arctan(sigma) - ln(|b + sigma|/|b - sigma|)/b],      0 < r < M,

with dx_*/dr = sqrt(k/8)/(r h sqrt(-f)), zero at the singularity r = 0, sqrt(k/8) pi at r = M and -infinity
at the inner horizon, whose surface gravity for d/dx is kappa_- = b/(2 sqrt(k/8)). With w = x -+ x_*,

    W(r) = exp(b(pi - 2 arctan(sigma))) (b - sigma)/(b + sigma) = sign(Q^2/M - r) exp(2 kappa_- x_*)

is smooth through r = Q^2/M, and the cells of that plane are the middle of Reissner and Nordstrom's tower:
r = 0 at UV = -1, the vertical lines X = +-pi/2, and r = M at UV = e^(b pi), where the plane of x and r
turns spacelike.

The slice t = 0 of the black string is drawn with x a circle of length 2 pi M, and the slice t = 0 beyond
extremality, at |Q| = 2M, with x of the period that closes its circles smoothly at rho = 0, Horne and
Horowitz's cigar.
"""
import math

import numpy as np

METRIC = "three_dim_black_string"
STATIC, WZW, IN = "static", "wzw", "eddington_finkelstein_ingoing"
EXTREMAL, EXTREMAL_IN, SUPER = "extremal", "extremal_interior", "superextremal"
M, Q2, K = 1.0, 0.5, 52 / 23
Q = math.sqrt(Q2)
PARAMS = {"M": "1", "Q": "sqrt(2)/2", "k": "52/23"}
EXTREMAL_PARAMS = {"M": "1", "k": "52/23"}
SUPER_PARAMS = {"M": "1", "Q": "2", "k": "52/23"}
QS = 2.0                                        # the charge beyond extremality
A = math.sqrt(M * M - Q2) / M
B = math.sqrt(M * M - Q2) / Q
SQ = math.sqrt(K / 8)
KAPPA = math.sqrt(2 * (M * M - Q2) / K) / M     # d/dt at r = M
KAPPA_IN = B / (2 * SQ)                         # d/dx at r = Q^2/M
INNER = Q2 / M
TOP = 6.0                                       # how far out the slice t = 0 is embedded, in M
SCALE = QS * math.sqrt(K / (2 * (QS * QS - M * M)))   # the period of x over 2 pi beyond extremality
SCALE_TEXT = "2*sqrt(26/69)"                         # the same, exactly, at Q = 2M and k = 52/23
assert abs(SCALE - 2 * math.sqrt(26 / 69)) < 1e-15
SETTINGS = "$M = 1$, $Q = M/\\sqrt{2}$, and $k = 52/23$."
SUPER_SETTINGS = "$M = 1$, $Q = 2M$, and $k = 52/23$."
EXTREMAL_SETTINGS = "$M = 1$, $|Q| = M$, and $k = 52/23$."


def s_of(r):
    return np.sqrt(1 - Q2 / (M * np.asarray(r, dtype=float)))


def rstar(r):
    """The tortoise coordinate of the plane of t and r, zero at the inner horizon."""
    s = s_of(r)
    with np.errstate(divide="ignore"):
        return SQ * (np.log(np.abs(s - A) / (s + A)) / A + np.log((1 + s) / (1 - s)))


def signed(r):
    """F(r) = sign(r - M) exp(2 kappa r_*), smooth through the event horizon."""
    s = s_of(r)
    return (s - A) / (s + A) * ((1 + s) / (1 - s)) ** A


def xstar(r):
    """The tortoise coordinate of the plane of x and r inside the event horizon, zero at r = 0."""
    sigma = np.sqrt(M / np.asarray(r, dtype=float) - 1)
    with np.errstate(divide="ignore"):
        return SQ * (math.pi - 2 * np.arctan(sigma) - np.log(np.abs(B + sigma) / np.abs(B - sigma)) / B)


def inner_signed(r):
    """W(r) = sign(Q^2/M - r) exp(2 kappa_- x_*), smooth through the inner horizon."""
    sigma = np.sqrt(M / np.asarray(r, dtype=float) - 1)
    return np.exp(B * (math.pi - 2 * np.arctan(sigma))) * (B - sigma) / (B + sigma)


# -- Spacetime diagrams ----------------------------------------------------------------------------

# Kruskal's U + V in the group's u and v at the drawing's values, U = v sqrt(G) and V = u sqrt(G) with G of
# group() below: u + v alone is no time function far out when Q is not zero, and U + V is one everywhere.
GROUP_R = "(1 - u*v/2)"
GROUP_S = f"sqrt(1 - 1/(2*{GROUP_R}))"
GROUP_TIME = (f"(u + v)*sqrt(1/(4*{GROUP_R}*({GROUP_S} + sqrt(2)/2)**2)"
              f"*((1 + {GROUP_S})/(1 - {GROUP_S}))**(sqrt(2)/2))")


def register_null(nr):
    rows = (
        nr.Diagram(METRIC, STATIC, "radial", "$t$ and $r$", ("t", "r"), (0, 2, -1, 1), "$r/M$", "$t/M$",
                   PARAMS, {"x": "0"}, orient="ingoing"),
        nr.Diagram(METRIC, STATIC, "along", "$x$ and $r$", ("x", "r"), (0, 1, -0.5, 0.5), "$r/M$", "$x/M$",
                   PARAMS, {"t": "0"}, orient="ingoing"),
        nr.Diagram(METRIC, WZW, "plane", "$u$ and $v$", ("u", "v"), (-2, 2, -2, 2), "$(v - u)/2$", "$(u + v)/2$",
                   PARAMS, {"x": "0"}, to_display=nr.NULL_TO_TR, tau=GROUP_TIME, families=nr.SIDEWAYS,
                   where="1 - u*v"),
        nr.Diagram(METRIC, IN, "chart", "against $v$", ("v", "r"), (0.5, 2.5, -1, 1), "$r/M$", "$v/M$",
                   PARAMS, {"x": "0"}, orient="ingoing"),
        nr.Diagram(METRIC, EXTREMAL, "radial", "$t$ and $\\rho$", ("t", "\\rho"), (0, 2, -1, 1),
                   "$\\rho/\\sqrt{M}$", "$t/M$", EXTREMAL_PARAMS, {"x": "0"}),
        nr.Diagram(METRIC, EXTREMAL_IN, "along", "$x$ and $\\rho$", ("x", "\\rho"), (0, 1, -0.5, 0.5),
                   "$\\rho/\\sqrt{M}$", "$x/M$", EXTREMAL_PARAMS, {"t": "0"}, tau="x"),
        nr.Diagram(METRIC, SUPER, "radial", "$t$ and $\\rho$", ("t", "\\rho"), (0, 2, -1, 1),
                   "$\\rho/\\sqrt{M}$", "$t/M$", SUPER_PARAMS, {"x": "0"}),
    )
    nr.DIAGRAMS.extend(rows)
    nr.CAPTIONS[(METRIC, STATIC, "radial")] = [
        "The plane of $t$ and $r$ ($x = 0$, $Q^2 = M^2/2$), its cones closing at the event horizon $r = M$ and "
        "tipping inward between the horizons; inside the inner horizon $r = M/2$ the plane is spacelike."]
    nr.CAPTIONS[(METRIC, STATIC, "along")] = [
        "The plane of $x$ and $r$ inside the event horizon ($t = 0$, $Q^2 = M^2/2$), where $r$ is the time down to "
        "the inner horizon $r = M/2$ and $x$ the time inside it, out to the timelike singularity $r = 0$."]
    nr.CAPTIONS[(METRIC, WZW, "plane")] = [
        "The plane of the group's $u$ and $v$ ($x = 0$, $Q^2 = M^2/2$), across both exteriors and through the "
        "event horizon $uv = 0$ up to the inner horizon $uv = 1$."]
    nr.CAPTIONS[(METRIC, IN, "chart")] = [
        "The plane of $v$ and $r$ ($x = 0$, $Q^2 = M^2/2$), the ingoing rays crossing the event horizon $r = M$ on "
        "toward the inner horizon $r = M/2$."]
    nr.CAPTIONS[(METRIC, EXTREMAL, "radial")] = [
        "The plane of $t$ and $\\rho$ at $|Q| = M$ ($x = 0$), the cones closing toward the event horizon "
        "$\\rho = 0$, an infinite distance down the throat."]
    nr.CAPTIONS[(METRIC, EXTREMAL_IN, "along")] = [
        "The plane of $x$ and $\\rho$ inside the horizon at $|Q| = M$ ($t = 0$), $x$ the time throughout, from the "
        "horizon $\\rho = 0$ out to the singularity $\\rho = \\sqrt{M}$."]
    nr.CAPTIONS[(METRIC, SUPER, "radial")] = [
        "The plane of $t$ and $\\rho$ at $|Q| = 2M$ ($x = 0$), static throughout, with no horizon and the smooth "
        "centre $\\rho = 0$."]

    def own(f):
        return lambda w, r: f(np.asarray(w, dtype=float), np.asarray(r, dtype=float))

    outside = own(lambda t, r: np.abs(r - M) > 0.05)
    nr.CLOSED_FORMS[(METRIC, STATIC, "radial")] = (
        own(lambda t, r: t + rstar(r)), own(lambda t, r: t - rstar(r)),
        own(lambda t, r: (np.abs(r - M) > 0.05) & (r > INNER + 0.02)))
    nr.CLOSED_FORMS[(METRIC, STATIC, "along")] = (
        own(lambda x, r: x - xstar(r)), own(lambda x, r: x + xstar(r)),
        own(lambda x, r: (np.abs(r - INNER) > 0.03) & (r < M - 0.02) & (r > 0.02)))
    nr.CLOSED_FORMS[(METRIC, IN, "chart")] = (
        own(lambda v, r: v), own(lambda v, r: v - 2 * rstar(r)),
        own(lambda v, r: (np.abs(r - M) > 0.05) & (r > INNER + 0.02)))

    def extremal_star(rho):
        # dt/drho = sqrt(k/2) sqrt(rho^2 + M)/rho^2 on the cone.
        return math.sqrt(K / 2) * (-np.sqrt(rho * rho + M) / rho + np.arcsinh(rho / math.sqrt(M)))

    def extremal_inside_star(rho):
        # dx/drho = sqrt(k/2) sqrt(M - rho^2)/rho^2 on the cone.
        return math.sqrt(K / 2) * (-np.sqrt(M - rho * rho) / rho - np.arcsin(rho / math.sqrt(M)))

    nr.CLOSED_FORMS[(METRIC, EXTREMAL, "radial")] = (
        own(lambda t, rho: t + extremal_star(rho)), own(lambda t, rho: t - extremal_star(rho)),
        own(lambda t, rho: rho > 0.05))
    nr.CLOSED_FORMS[(METRIC, EXTREMAL_IN, "along")] = (
        own(lambda x, rho: x + extremal_inside_star(rho)), own(lambda x, rho: x - extremal_inside_star(rho)),
        own(lambda x, rho: (rho > 0.05) & (rho < 0.98)))


# -- The slices on the spacetime diagrams --------------------------------------------------------------

def register_slices(sl):
    def string():
        return sl.moments(METRIC, "string")[0]

    def static():
        m = string()
        return [sl.Mark(m, sl.along(0.0, *m.reach(STATIC, "r")))]

    def wzw():
        # t = 0 is u = -v, and r = M + (M^2 - Q^2)u^2/M on it.
        m = string()
        _, hi = m.reach(STATIC, "r")
        u = math.sqrt(M * (hi - M) / (M * M - Q2))
        return [sl.Mark(m, [[(-u, u), (u, -u)]])]

    def ingoing():
        m = string()
        lo, hi = m.reach(STATIC, "r")
        r = sl.near(lo, hi)
        return [sl.Mark(m, [np.column_stack([rstar(r), r])])]

    def cigar():
        m = sl.moments(METRIC, "cigar")[0]
        return [sl.Mark(m, sl.along(0.0, *m.reach(SUPER, "\\rho")))]

    sl.FLAT[(METRIC, STATIC, "radial")] = static
    sl.FLAT[(METRIC, WZW, "plane")] = wzw
    sl.FLAT[(METRIC, IN, "chart")] = ingoing
    sl.FLAT[(METRIC, SUPER, "radial")] = cigar
    for key in ((METRIC, STATIC, "along"), (METRIC, EXTREMAL, "radial"), (METRIC, EXTREMAL_IN, "along")):
        sl.FLAT[key] = lambda: []
    sl.FLAT_METRICS.add(METRIC)


# -- Embedding diagram ----------------------------------------------------------------------------------

def register_embedding(eb):
    eb.DRAWN[METRIC] = lambda ck, src: embedding(eb, ck, src)
    eb.CAPTIONS[(METRIC, "string")] = [
        "The slice $t = 0$ of the black string ($Q^2 = M^2/2$), $x$ a circle of length $2\\pi M$: two exteriors "
        "joined at the event horizon, the narrowest circle."]
    eb.CAPTIONS[(METRIC, "cigar")] = [
        "The slice $t = 0$ beyond extremality ($|Q| = 2M$), $x$ a circle that closes smoothly at $\\rho = 0$: an "
        "infinite cigar."]


def embedding(eb, ck, src):
    """The black string's slice t = 0 from the event horizon out to r = 6M, each exterior, joined at r = M,
    where the circle of x is narrowest, rho = M sqrt(1 - Q^2/Mr) with x of period 2 pi M; and the slice t = 0
    at |Q| = 2M, where x has the period 2 pi SCALE that makes rho = 0 a smooth centre, rho_circle =
    SCALE rho sqrt(M/(Q^2 + M rho^2))."""
    sl = eb.Slice(src, METRIC, STATIC, "r", "x", {"t": 0}, PARAMS)
    ck.stops("Charged black string, between the horizons", sl, np.linspace(INNER, M, 402)[1:-1])
    size = 2.6
    marks = [(2.0, "r", None), (4.0, "r", None), (TOP, "r", "$r = 6M$")]
    out = eb.Piece("exterior", "sheet", sl, M, TOP, 0.0, 1,
                   (("throat", "the event horizon $r = M$, the narrowest circle, where the slice runs on into the "
                               "other exterior"),
                    ("edge", "the surface runs on as $r \\to \\infty$ toward a cylinder of radius $M$")),
                   [(M, "horizon", "$r = M$")] + marks, size)
    other = eb.Piece("other_exterior", "sheet2", sl, M, TOP, 0.0, -1,
                     (("throat", "the event horizon $r = M$"), ("edge", "the other exterior, as $r \\to \\infty$")),
                     [(r, "r2", None) for r, _, _ in marks], size)
    for p in (out, other):
        ck.isometry(f"Charged black string, {p.id}", p)
        ck.radius(f"Charged black string, {p.id}, rho = M sqrt(1 - Q^2/Mr)", p,
                  lambda r: M * np.sqrt(1 - Q2 / (M * r)), size)
    ck.join("Charged black string, the two exteriors at r = M", out, M, other, M)
    string = eb.Surface([out, other])
    fig = eb.figure_of([string], {"sheet": "cover", "sheet2": "cover"}, size, eb.Camera(-90, 18))
    eb.ring_label(fig, [0, 0, 0], *out.at(M), "$r = M$", dx=10)
    eb.ring_label(fig, [0, 0, 0], *out.at(TOP), "$6M$")
    fig.legend("fill", "cover", "the slice $t = 0$ outside the event horizon")
    fig.legend("line", "horizon", "the event horizon $r = M$")
    fig.legend("line", "r", "$r$ constant, at $2M$, $4M$, and $6M$")
    fig.legend("line", "r2", "the same circles in the other exterior")
    fig.legend("line", "meridian", "$x$ constant, every $\\pi M/12$")
    stops = ["Inside the event horizon $r$ is a time, so a slice of constant $t$ is a timelike surface."]
    views = [eb.view("string", "The black string", "$M$", [string], fig.done(),
                     settings=SETTINGS, stops=stops)]

    sc = eb.Slice(src, METRIC, SUPER, "\\rho", "x", {"t": 0}, SUPER_PARAMS, swept={"x": f"{SCALE_TEXT}*x"})
    top = 6.0
    size = 2 * SCALE + 0.6
    cigar = eb.Piece("cigar", "sheet", sc, 0.0, top, 0.0, 1,
                     (("axis", "the centre $\\rho = 0$, where the circle of $x$ closes to a point"),
                      ("edge", "the surface runs on as $\\rho \\to \\infty$ toward a cylinder")),
                     [(rho, "r", None) for rho in (1.0, 2.0, 4.0)] + [(top, "r", "$\\rho = 6\\sqrt{M}$")], size)
    ck.isometry("Charged string beyond extremality, the cigar", cigar)
    ck.radius("Charged string beyond extremality, the cigar, rho = SCALE rho sqrt(M/(Q^2 + M rho^2))", cigar,
              lambda rho: SCALE * rho * np.sqrt(M / (QS * QS + M * rho * rho)), size)
    tip = sc.slope(0.0, "+")
    ck.add("Charged string beyond extremality: the cigar closes smoothly at rho = 0, drho/ds = 1",
           float(np.hypot(tip[0] - 1, tip[1])), 1e-9)
    ck.add("Charged string beyond extremality: far away the circle's radius is the period over 2 pi",
           abs(float(sc.rho_at(1e9)) - SCALE), 1e-6)
    surface = eb.Surface([cigar])
    fig = eb.figure_of([surface], {"sheet": "cover"}, size)
    eb.ring_label(fig, [0, 0, 0], *cigar.at(2.0), "$\\rho = 2\\sqrt{M}$", dx=10)
    eb.ring_label(fig, [0, 0, 0], *cigar.at(top), "$6\\sqrt{M}$")
    fig.legend("fill", "cover", "the slice $t = 0$, which $\\rho$ and $x$ cover")
    fig.legend("line", "r", "$\\rho$ constant, at $1$, $2$, $4$, and $6\\,\\sqrt{M}$")
    fig.legend("line", "meridian", "$x$ constant, every $1/24$ of its period")
    views.append(eb.view("cigar", "Beyond extremality", "$M$", [surface], fig.done(),
                         settings=SUPER_SETTINGS[:-1] + ", with $x$ of period $\\pi|Q|\\sqrt{2k/(Q^2 - M^2)}$."))
    return views


# -- Conformal diagram ----------------------------------------------------------------------------------

HALF, PI = math.pi / 2, math.pi
OUTLINE = [[-PI, 0.0], [-HALF, HALF], [HALF, HALF], [PI, 0.0], [HALF, -HALF], [-HALF, -HALF]]
CELLS = {"I": [[0.0, 0.0], [HALF, HALF], [PI, 0.0], [HALF, -HALF]],
         "II": [[0.0, 0.0], [HALF, HALF], [-HALF, HALF]],
         "III": [[0.0, 0.0], [-HALF, HALF], [-PI, 0.0], [-HALF, -HALF]],
         "IV": [[0.0, 0.0], [HALF, -HALF], [-HALF, -HALF]]}


def square(cf, cell, t, r):
    """(p, q) of the points (t, r) of a cell of the plane of t and r, t the static coordinate there."""
    t, rs = np.asarray(t, dtype=float), rstar(r)
    u, v = t - rs, t + rs
    p, q = cf.atan_exp(-KAPPA * u), cf.atan_exp(KAPPA * v)
    return {"I": (-p, q), "II": (p, q), "III": (p, -q), "IV": (-p, -q)}[cell]


def group(u, v):
    """(p, q) of the group's (u, v): U = v sqrt(G), V = u sqrt(G), G = -F(r)/uv, which is smooth at r = M."""
    u, v = np.asarray(u, dtype=float), np.asarray(v, dtype=float)
    r = M - (M * M - Q2) * u * v / M
    s = s_of(r)
    G = Q2 / (M * M * r * (s + A) ** 2) * ((1 + s) / (1 - s)) ** A * (M * M - Q2) / M
    return np.arctan(v * np.sqrt(G)), np.arctan(u * np.sqrt(G))


def finkelstein(cf, v, r):
    """(p, q) of the ingoing chart's (v, r): V = exp(kappa v) and U = -exp(-kappa v) F(r)."""
    v = np.asarray(v, dtype=float)
    return np.arctan(-np.exp(-KAPPA * v) * signed(r)), cf.atan_exp(KAPPA * v)


def middle(cf, cell, x, r):
    """(p, q) of the points (x, r) of a cell of the plane of x and r inside the event horizon."""
    x, xs = np.asarray(x, dtype=float), xstar(r)
    p, q = cf.atan_exp(-KAPPA_IN * (x - xs)), cf.atan_exp(KAPPA_IN * (x + xs))
    return {"right": (-p, q), "left": (p, -q), "below": (-p, -q), "above": (p, q)}[cell]


def register_conformal(cf):
    cf.DRAWN[METRIC] = lambda ck, src: conformal(cf, ck, src)
    cf.CAPTIONS[(METRIC, "static")] = [
        "The plane of $t$ and $r$ outside the inner horizon ($Q^2 = M^2/2$), each point a line along $x$: two "
        "exteriors through the event horizon, ending at $r = Q^2/M$, where $t$ and $r$ turn spacelike."]
    cf.CAPTIONS[(METRIC, "interior")] = [
        "The plane of $x$ and $r$ inside the event horizon ($Q^2 = M^2/2$), each point a line along $t$, with the "
        "timelike singularity $r = 0$ beyond the inner horizon."]
    cf.CAPTIONS[(METRIC, "wzw")] = [
        "The plane of $t$ and $r$ outside the inner horizon ($Q^2 = M^2/2$), each point a line along $x$, all of "
        "it covered by the group's $u$ and $v$, $uv < 1$."]
    cf.CAPTIONS[(METRIC, "ingoing")] = [
        "The plane of $t$ and $r$ outside the inner horizon ($Q^2 = M^2/2$), each point a line along $x$, with the "
        "ingoing rays of constant $v$ crossing the event horizon."]


def conformal(cf, ck, src):
    """Two planes of the string, neither of which carries its causal structure alone: the plane of t and r,
    Lorentzian outside the inner horizon, drawn as Kruskal's square from the closed form r_* of the module
    docstring under every chart that covers it, and the plane of x and r, Lorentzian inside the event
    horizon, drawn from x_* as the middle of Reissner and Nordstrom's tower under the static chart."""
    st = cf.Plane(src, METRIC, STATIC, ("t", "r"), {"x": "0"}, PARAMS)
    gr = cf.Plane(src, METRIC, WZW, ("u", "v"), {"x": "0"}, PARAMS)
    inn = cf.Plane(src, METRIC, IN, ("v", "r"), {"x": "0"}, PARAMS)
    xr = cf.Plane(src, METRIC, STATIC, ("x", "r"), {"t": "0"}, PARAMS)
    gap = 0.01
    for cell, rs, future, name in (("I", (M + gap, 30), (1, 0), "exterior"),
                                   ("II", (INNER + gap, M - gap), (0, -1), "black hole"),
                                   ("III", (M + gap, 30), (-1, 0), "other exterior"),
                                   ("IV", (INNER + gap, M - gap), (0, 1), "white hole")):
        ck.chart(f"charged string, t and r, {name}", st, lambda t, r, c=cell: square(cf, c, t, r),
                 ck.uniform(-20, 20), ck.uniform(*rs), lambda t, r, f=future: f)
    u, v = ck.uniform(-3, 3), ck.uniform(-3, 3)
    keep = u * v < 0.98
    ck.chart("charged string, the group's u and v", gr, group, u[keep], v[keep], lambda u, v: (1, 1))
    ck.chart("charged string, ingoing Eddington-Finkelstein", inn, lambda v, r: finkelstein(cf, v, r),
             ck.uniform(-20, 20), ck.uniform(INNER + gap, 30), lambda v, r: (1, -1))
    for cell, rs, future, name in (("right", (gap, INNER - gap), (1, 0), "inside the inner horizon"),
                                   ("left", (gap, INNER - gap), (-1, 0), "the other side"),
                                   ("below", (INNER + gap, M - gap), (0, -1), "the black hole"),
                                   ("above", (INNER + gap, M - gap), (0, 1), "the region above")):
        ck.chart(f"charged string, x and r, {name}", xr, lambda x, r, c=cell: middle(cf, c, x, r),
                 ck.uniform(-20, 20), ck.uniform(*rs), lambda x, r, f=future: f)

    rr = np.array([0.55, 0.7, 0.95, 1.05, 2.0, 5.0])
    for cell, sel in (("I", rr > M), ("II", rr < M), ("III", rr > M), ("IV", rr < M)):
        p, q = square(cf, cell, 0.4 + 0 * rr[sel], rr[sel])
        ck.limit(f"charged string {cell}: tan p tan q is UV = -F(r)", np.tan(p) * np.tan(q), -signed(rr[sel]), 1e-9)
    p, q = square(cf, "II", np.array([-3.0, 0.0, 3.0]), np.full(3, INNER))
    ck.limit("charged string: the inner horizon r = Q^2/M lands on T = pi/2", p + q, [HALF] * 3, 1e-12)
    p, q = square(cf, "I", np.array([-3.0, 0.0, 3.0]), np.full(3, 1e12))
    # r_* grows as sqrt(k/8) ln r without bound, so a point of fixed t runs to i0, short of it by
    # 2 arctan(exp(-kappa r_*)), which falls as r^(-a/2): about 1.3e-4 at r = 1e12.
    ck.limit("charged string: r_* grows as sqrt(k/8) ln r", (rstar(1e8) - rstar(1e5)) / math.log(1e3), SQ, 1e-6)
    far, ts = rstar(np.full(3, 1e12)), np.array([-3.0, 0.0, 3.0])
    ck.limit("charged string: r -> infinity at fixed t runs to i0", PI - (q - p),
             np.arctan(np.exp(KAPPA * (ts - far))) + np.arctan(np.exp(-KAPPA * (ts + far))), 1e-12)
    p, q = square(cf, "I", np.array([2.0]), np.array([M + 1e-12]))
    ck.limit("charged string: r -> M at fixed t lands on the bifurcation line", cf.point(p[0], q[0]), [0, 0], 1e-3)
    uu, vv = np.array([0.7, -1.2, 0.4]), np.array([-0.9, 0.5, 0.6])
    r_uv = M - (M * M - Q2) * uu * vv / M
    p, q = group(uu, vv)
    ck.limit("charged string: the group's chart has UV = -F(r)", np.tan(p) * np.tan(q), -signed(r_uv), 1e-9)
    t_uv = np.log(np.abs(uu / vv)) / (2 * KAPPA)
    p2, q2 = square(cf, "I", t_uv[:1], r_uv[:1])
    ck.limit("charged string: the group's and the static coordinates put one event at one point",
             np.array(group(uu[:1], vv[:1])).ravel(), np.array([p2[0], q2[0]]), 1e-9)
    ck.limit("charged string: the ingoing and static coordinates put one event at one point",
             np.array(finkelstein(cf, 1.5 + rstar(2.0), 2.0)).ravel(), np.array(square(cf, "I", 1.5, 2.0)).ravel(),
             1e-12)
    ck.limit("charged string: the ingoing chart is smooth through r = M",
             np.array(finkelstein(cf, 0.3, M)).ravel(), [0.0, float(cf.atan_exp(KAPPA * 0.3))], 1e-12)
    rm = np.array([0.1, 0.3, 0.45, 0.55, 0.8, 0.95])
    for cell, sel in (("right", rm < INNER), ("left", rm < INNER), ("below", rm > INNER), ("above", rm > INNER)):
        p, q = middle(cf, cell, 0.2 + 0 * rm[sel], rm[sel])
        ck.limit(f"charged string {cell}: tan p tan q is W(r)", np.tan(p) * np.tan(q), inner_signed(rm[sel]), 1e-9)
    for cell, side in (("right", 1), ("left", -1)):
        p, q = middle(cf, cell, np.array([-1.0, 0.0, 1.0]), np.full(3, 1e-14))
        ck.limit(f"charged string {cell}: r -> 0 lands on X = {side} pi/2", q - p, [side * HALF] * 3, 1e-6)
    K_tr, K_xr = st.kretschmann, xr.kretschmann
    ck.finite("charged string: the Kretschmann scalar is finite at r = Q^2/M",
              K_tr(np.zeros(3), np.array([INNER - 1e-3, INNER, INNER + 1e-3])))
    ck.finite("charged string: the Kretschmann scalar is finite at r = M",
              K_xr(np.zeros(3), np.array([M - 1e-3, M, M + 1e-3])))
    ck.diverges("charged string: the Kretschmann scalar diverges at r = 0", K_xr(0, 1e-3), K_xr(0, 1e-4))

    moment = cf.slices.moments(METRIC, "string")[0]
    lo, hi = moment.reach(STATIC, "r")
    r_m = np.linspace(lo + 1e-9, hi, 400)
    left, right = square(cf, "III", 0 * r_m, r_m[::-1]), square(cf, "I", 0 * r_m, r_m)
    through = (np.concatenate([left[0], right[0]]), np.concatenate([left[1], right[1]]))

    box = [-PI - 0.25, PI + 0.25, -HALF - 0.3, HALF + 0.3]
    t_all = cf.S_ALL
    between = np.linspace(INNER + 1e-6, M - 1e-9, 600)
    beyond = cf.spread(M, np.inf, 600, 14)
    radii_in, radii_out = (0.6, 0.75, 0.9), (1.5, 3.0, 10.0)
    times = (-6, -3, 0, 3, 6)

    def edges(v):
        v.line("horizon", [[[-HALF, -HALF], [HALF, HALF]], [[HALF, -HALF], [-HALF, HALF]]])
        v.line("chartedge", [[[-HALF, HALF], [HALF, HALF]], [[-HALF, -HALF], [HALF, -HALF]]])
        v.line("scri", [[[HALF, HALF], [PI, 0.0], [HALF, -HALF]], [[-HALF, HALF], [-PI, 0.0], [-HALF, -HALF]]])
        v.label_xt([0.0, HALF], "$r = Q^2/M$", "b", "small", dy=-4)
        v.label_xt([HALF / 2, HALF / 2], "$r = M$", "br", "small", dx=-4, dy=-2)
        v.label_xt([3 * HALF / 2, HALF / 2], "$r \\to \\infty$", "bl", "small", dx=4, dy=-3)
        v.legend("horizon", "the event horizon $r = M$")
        v.legend("chartedge", "$r = Q^2/M$, past which $t$ and $r$ are both spacelike")
        v.legend("scri", "null infinity, $r \\to \\infty$")

    def static_grid(v, cells):
        for cell in cells:
            inside = cell in ("II", "IV")
            for r in (radii_in if inside else radii_out):
                v.curve("r", *square(cf, cell, t_all, np.full_like(t_all, r)))
            span = between if inside else beyond
            for t in times:
                v.curve("t", *square(cf, cell, np.full_like(span, t), span))

    views = []
    for vid, label, system, cells, cover in (("static", "$t$ and $r$", STATIC, ("I", "II", "III", "IV"), ("I", "II")),
                                             ("wzw", "Wess-Zumino-Witten", WZW, ("I", "II", "III", "IV"),
                                              ("I", "II", "III", "IV")),
                                             ("ingoing", "Ingoing Eddington-Finkelstein", IN, ("I", "II", "III", "IV"),
                                              ("I", "II"))):
        v = cf.View(vid, label, box, system)
        v.fill("region", OUTLINE)
        for cell in cover:
            v.fill("cover", CELLS[cell])
        if vid == "ingoing":
            radii = cf.spread(INNER, np.inf, 800, 14)
            for r in radii_in + radii_out:
                v.curve("r", *finkelstein(cf, t_all, np.full_like(t_all, r)))
            for w in times:
                v.curve("null", *finkelstein(cf, np.full_like(radii, w), radii))
            v.legend("null", "$v$ constant, an ingoing light ray")
        elif vid == "wzw":
            span = np.linspace(-6, 6, 801)
            for c in (-1.0, -0.5, 0.5, 1.0):
                for which in ("u", "v"):
                    a = np.full_like(span, c)
                    w = span[a * span < 0.999]
                    a = a[: len(w)]
                    v.curve("t" if which == "u" else "r", *(group(a, w) if which == "u" else group(w, a)))
            v.legend("t", "$u$ constant")
            v.legend("r", "$v$ constant")
        else:
            static_grid(v, cells)
            v.legend("r", "$r$ constant")
            v.legend("t", "$t$ constant")
        edges(v)
        v.legend("cover", {"static": "what $t$ and $r$ cover", "wzw": "what $u$ and $v$ cover",
                           "ingoing": "what $v$ and $r$ cover"}[vid])
        v.slice(moment, [through])
        v.set(settings=SETTINGS)
        views.append(v)

    # The plane of x and r: r = 0 on X = +-pi/2, r = M on tan p tan q = e^(b pi), above and below.
    edge_p = np.linspace(-HALF + 1e-6, -1e-6, 400)
    edge_q = np.arctan(math.exp(B * PI) / np.tan(edge_p))
    bottom = [[float(q - p), float(q + p)] for p, q in zip(edge_p, edge_q)]
    top = [[-X, -T] for X, T in bottom]
    height = max(-T for _, T in bottom)
    ck.limit("charged string: r = M in the plane of x and r is W = e^(b pi)", np.tan(edge_p) * np.tan(edge_q),
             np.full(edge_p.size, math.exp(B * PI)), 1e-6)
    outline = [[HALF, -HALF], [HALF, HALF]] + top + [[-HALF, HALF], [-HALF, -HALF]] + bottom
    v = cf.View("interior", "$x$ and $r$", [-HALF - 0.6, HALF + 0.6, -height - 0.25, height + 0.25], STATIC)
    v.fill("region", outline)
    v.fill("cover", [[0.0, 0.0], [HALF, HALF], [HALF, -HALF]])
    v.fill("cover", [[0.0, 0.0], [HALF, -HALF]] + bottom[::-1][1:-1] + [[-HALF, -HALF]])
    xs = cf.S_ALL
    inside_r = np.linspace(1e-6, INNER - 1e-9, 600)
    between_r = np.linspace(INNER + 1e-9, M - 1e-6, 600)
    for cell in ("right", "left", "below", "above"):
        span = inside_r if cell in ("right", "left") else between_r
        for r in ((0.15, 0.3, 0.42) if cell in ("right", "left") else (0.6, 0.75, 0.9)):
            v.curve("r", *middle(cf, cell, xs, np.full_like(xs, r)))
        for x in (-3, -1.5, 0, 1.5, 3):
            v.curve("t", *middle(cf, cell, np.full_like(span, x), span))
    v.line("singular", [[[-HALF, -HALF], [-HALF, HALF]], [[HALF, -HALF], [HALF, HALF]]], zig=True)
    v.line("horizon", [[[-HALF, -HALF], [HALF, HALF]], [[HALF, -HALF], [-HALF, HALF]]])
    v.line("chartedge", [bottom, top])
    v.label_xt([HALF, 0.0], "$r = 0$", "l", "small", dx=8)
    v.label_xt([HALF / 2, HALF / 2], "$r = Q^2/M$", "br", "small", dx=-4, dy=-2)
    v.label_xt([0.0, -height], "$r = M$", "t", "small", dy=6)
    v.legend("cover", "what $x$ and $r$ cover")
    v.legend("r", "$r$ constant")
    v.legend("t", "$x$ constant")
    v.legend("horizon", "the inner horizon $r = Q^2/M$")
    v.legend("singular", "$r = 0$, a curvature singularity")
    v.legend("chartedge", "$r = M$, past which $x$ and $r$ are both spacelike")
    v.set(settings=SETTINGS)
    views.insert(1, v)
    return views


# -- Curvature figure -----------------------------------------------------------------------------------

def register_curvature(cv):
    cv.DRAWN[METRIC] = lambda ck, src: curvature(cv, ck, src)


def curvature(cv, ck, src):
    """K and R against r from the stored scalars of the static chart: K positive everywhere, finite through
    both horizons and growing as r^-4 toward r = 0; R = 4(2M^2 r - 7MQ^2 + 2Q^2 r)/(Mkr^2), which changes sign
    at r = 7MQ^2/(2M^2 + 2Q^2) = 7M/6 and grows as -r^-2 toward r = 0."""
    params = {"M": M, "Q": Q, "k": K}
    names = list(params)
    Kp = cv.Published(src, METRIC, STATIC, "kretschmann", ["r"], names)
    Rp = cv.Published(src, METRIC, STATIC, "ricci_scalar", ["r"], names)
    x_axis = cv.axis("$r/M$", 0.1, 3)
    K_scale = {"kind": "log", "from": 0.1, "to": 1e6}
    R_scale = {"kind": "symlog", "from": -1e3, "to": 10, "linear": 0.1}
    panels = [{"id": "K", "name": "$K\\,M^4$", "aspect": 2, "scale": K_scale, "ticks": cv.value_ticks(K_scale)},
              {"id": "R", "name": "$R\\,M^2$", "aspect": 2, "scale": R_scale, "ticks": cv.value_ticks(R_scale)}]
    cv.shape_checks(ck, METRIC, [p["aspect"] for p in panels])
    height = {p["id"]: cv.PLOT_W / p["aspect"] for p in panels}
    curves = (cv.line_curves(ck, METRIC, "K", "K", Kp, "r", x_axis, K_scale, height["K"], params)
              + cv.line_curves(ck, METRIC, "R", "R", Rp, "r", x_axis, R_scale, height["R"], params))
    vals = list(params.values())
    rs = np.linspace(0.1, 3, cv.SAMPLES)
    R_hand = 4 * (2 * M * M * rs - 7 * M * Q2 + 2 * Q2 * rs) / (M * K * rs ** 2)
    K_hand = 16 * (4 * M ** 4 * rs ** 2 - 12 * M ** 3 * Q2 * rs + 19 * M * M * Q2 ** 2 - 12 * M * Q2 ** 2 * rs
                   + 4 * Q2 ** 2 * rs ** 2) / (M * M * K * K * rs ** 4)
    ck("Charged black string: R = 4(2M^2 r - 7MQ^2 + 2Q^2 r)/(Mkr^2)",
       bool(np.allclose(Rp(rs, *vals), R_hand, rtol=1e-12)))
    ck("Charged black string: K as Horne and Horowitz's curvature gives it", bool(np.allclose(Kp(rs, *vals), K_hand, rtol=1e-12)))
    ck("Charged black string: K is positive", bool(np.all(Kp(rs, *vals) > 0)))
    ck("Charged black string: R vanishes at r = 7M/6", abs(float(Rp.exact([7 / 6] + vals))) < 1e-12)
    ck("Charged black string: K is finite at both horizons",
       bool(np.all(np.isfinite(Kp(np.array([INNER, M]), *vals)))))
    legend = [["line", "K", "$K$"], ["line", "R", "$R$"],
              ["line", "horizon", "the horizons $r = M/2$ and $r = M$"]]
    return [{
        "id": "radial", "label": "Along $r$", "kind": "line", "x": x_axis,
        "panels": panels, "curves": curves, "marks": [{"class": "horizon", "x": INNER}, {"class": "horizon", "x": M}],
        "legend": legend,
        "caption": ["The Kretschmann scalar $K$ and the Ricci scalar $R$ of the string against $r$ ($Q^2 = M^2/2$). "
                    "$K$ is finite at both horizons and diverges at $r = 0$. "
                    "$R$ changes sign at $r = 7M/6$ and diverges at $r = 0$."],
        "ricci_sentences": [[0, 2]],
        "per_scalar": cv.scalar_clips("A logarithmic scale from $10^{-1}$ to $10^6$.",
                                      "A signed logarithmic scale from $-10^3$ to $10$, linear within $10^{-1}$."),
        "settings": SETTINGS,
        "clip": "$K$ logarithmic from $10^{-1}$ to $10^6$; $R$ signed logarithmic from $-10^3$ to $10$, "
                "linear within $10^{-1}$.",
    }]
