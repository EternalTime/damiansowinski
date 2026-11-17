"""Drawings of Chan, Horne and Mann's charged dilaton black holes, the hole embedded at a = 1, r_0 = r_h.

At the coupling a = 1 of string theory and r_0 = r_h = 1, the unit of every length, the first chart has
h = r and U = r - 1, so the tortoise coordinate dr_*/dr = 1/U is

    r_* = ln|r - 1|,

which runs to +infinity with r: light reaches no r -> infinity at finite t, and the plane of t and r is
drawn as Kruskal and Szekeres's square with kappa = U'(1)/2 = 1/2, UV = -(r - 1) outside and UV = 1 - r
inside, so that the singularity r = 0 is UV = 1, as Schwarzschild's is. At a^2 = 1/3, h = r^{3/2} and
U = r^{1/2}(r - 1), so

    r_* = ln|(sqrt(r) - 1)/(sqrt(r) + 1)|,

which tends to 0 as r -> infinity: a ray reaches r -> infinity at a finite t, Chan, Horne and Mann's
section 3.2, and in Kruskal's U = -exp(-u/2), V = exp(v/2) the boundary r -> infinity is the hyperbola
UV = -1, a vertical line in the conformal diagram, and the singularity r = 0 is UV = 1, as before.

The embedding diagram draws the equatorial plane of the slice t = 0 at a = 1: its metric is
dr^2/(r - 1) + r dphi^2 with r_0 = 1, so the circle at r has the radius rho = sqrt(r), and
(dz/drho)^2 = 4r/(r - 1) - 1, which runs to 3 as r -> infinity: far away the surface is a cone, the
slice's deficit of angle, and it has a throat at the bifurcation sphere r = 1, where it runs on into a
second exterior.
"""
import math

import numpy as np

METRIC = "dilaton_unusual_asymptotics"
CHM, AREAL, ELECTRIC, MAGNETIC = "chan_horne_mann", "areal", "string_electric", "string_magnetic"
INVERSE, MINUS_A, MINUS_INVERSE, MASSLESS = ("liouville_inverse", "liouville_minus_a", "liouville_minus_inverse",
                                            "massless")
EQUATOR = {"theta": "pi/2", "phi": "0"}
ONE = {"a": "1", "r_0": "1", "r_h": "1"}                      # the hole embedded
THIRD = {"a": "sqrt(3)/3", "r_0": "1", "r_h": "1"}             # a^2 = 1/3, infinity timelike
STRING = {"r_0": "1", "r_h": "1", "Q": "sqrt(2)/2"}            # horizons at r = 1 in both string metrics
INVERSE_PARAMS = {"a": "1/2", "r_0": "1", "M": "4/3", "q": "1"}   # horizons r = 1/2 and 3/2
MINUS_A_PARAMS = {"a": "1/2", "r_0": "1", "r_h": "1", "sigma": "2"}
MINUS_INVERSE_PARAMS = {"a": "2", "r_0": "1", "r_h": "1", "lambda": "-1/10"}
TOP = 6.0                                                     # how far out the slice t = 0 is embedded


def rstar_one(r):
    """The tortoise coordinate at a = 1, r_0 = r_h = 1."""
    return np.log(np.abs(np.asarray(r, dtype=float) - 1))


def rstar_third(r):
    """The tortoise coordinate at a^2 = 1/3, r_0 = r_h = 1, zero at r -> infinity."""
    s = np.sqrt(np.asarray(r, dtype=float))
    return np.log(np.abs((s - 1) / (s + 1)))


def minus_inverse_u(r):
    """U at b = -1/a, a = 2, r_0 = r_h = 1 and lambda = -1/10."""
    r = np.asarray(r, dtype=float)
    h = r ** 0.4
    return h * (1 - 1 / r) - 0.1 * r * r / h


def bisect(f, lo, hi, n=200):
    a, b = lo, hi
    for _ in range(n):
        m = (a + b) / 2
        if (f(a) < 0) == (f(m) < 0):
            a = m
        else:
            b = m
    return (a + b) / 2


# The event horizon and the cosmological horizon at b = -1/a, a = 2, lambda = -1/10, and the radius
# between them where U is largest, at which the cones turn from ingoing to outgoing.
EVENT = bisect(lambda r: float(minus_inverse_u(r)), 1.0, 2.0)
COSMOLOGICAL = bisect(lambda r: float(minus_inverse_u(r)), 3.0, 10.0)
STATIC = bisect(lambda r: float(minus_inverse_u(r + 1e-7) - minus_inverse_u(r - 1e-7)), EVENT + 0.1, COSMOLOGICAL - 0.1)


def register_null(nr):
    rows = (
        nr.Diagram(METRIC, CHM, "radial", "$t$ and $r$, $a = 1$", ("t", "r"), (0, 6, -3, 3), "$r/r_h$", "$t/r_h$",
                   ONE, EQUATOR, orient="ingoing", areal=True),
        nr.Diagram(METRIC, CHM, "third", "$t$ and $r$, $a^2 = 1/3$", ("t", "r"), (0, 6, -3, 3), "$r/r_h$",
                   "$t/r_h$", THIRD, EQUATOR, orient="ingoing", areal=True),
        nr.Diagram(METRIC, AREAL, "radial", "$t$ and $r$", ("t", "r"), (0, 6, -3, 3), "$r/r_h$", "$t/r_h$",
                   {"r_0": "1", "r_h": "1"}, EQUATOR, orient="ingoing", areal=True),
        nr.Diagram(METRIC, ELECTRIC, "radial", "$t$ and $r$", ("t", "r"), (0, 6, -3, 3), "$r/r_h$", "$t/r_h$",
                   STRING, EQUATOR, orient="ingoing", areal=True),
        nr.Diagram(METRIC, MAGNETIC, "radial", "$t$ and $r$", ("t", "r"), (0, 6, -3, 3), "$r/r_h$", "$t/r_h$",
                   STRING, EQUATOR, orient="ingoing"),
        nr.Diagram(METRIC, INVERSE, "radial", "$t$ and $r$", ("t", "r"), (0, 3, -1.5, 1.5), "$r/r_0$", "$t/r_0$",
                   INVERSE_PARAMS, EQUATOR, orient="ingoing", areal=True),
        nr.Diagram(METRIC, MINUS_A, "radial", "$t$ and $r$", ("t", "r"), (0, 6, -3, 3), "$r/r_h$", "$t/r_h$",
                   MINUS_A_PARAMS, EQUATOR, orient="ingoing", areal=True),
        nr.Diagram(METRIC, MINUS_INVERSE, "radial", "$t$ and $r$", ("t", "r"), (0, 7, -3.5, 3.5), "$r/r_h$",
                   "$t/r_h$", MINUS_INVERSE_PARAMS, EQUATOR, orient="split", split=STATIC, areal=True),
        nr.Diagram(METRIC, MASSLESS, "plane", "$\\eta$ and $\\xi$", ("\\eta", "\\xi"), (0, 4, -2, 2), "$\\xi/r_0$",
                   "$\\eta/r_0$", {}, EQUATOR, tau="eta", where="xi**2 - eta**2"),
    )
    nr.DIAGRAMS.extend(rows)
    nr.CAPTIONS[(METRIC, CHM, "radial")] = [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) of the linear dilaton black hole ($a = 1$, "
        "$r_0 = r_h$), each point a 2-sphere of area $4\\pi r_0r$. The cones close at the horizon $r = r_h$ and "
        "point to the singularity $r = 0$ inside it."]
    nr.CAPTIONS[(METRIC, CHM, "third")] = [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) at $a^2 = 1/3$ ($r_0 = r_h$), where "
        "$dt/dr = 1/U$ falls as $r^{-3/2}$ and an outgoing ray reaches $r \\to \\infty$ at a finite $t$."]
    nr.CAPTIONS[(METRIC, AREAL, "radial")] = [
        "The plane of $t$ and the areal radius $r$ ($\\theta = \\pi/2$, $\\phi = 0$) of the linear dilaton black "
        "hole ($r_0 = r_h$), with the horizon at $r = \\sqrt{r_0r_h}$ and the singularity at $r = 0$."]
    nr.CAPTIONS[(METRIC, ELECTRIC, "radial")] = [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) of the string metric of the electric hole "
        "($a = 1$, $r_0 = r_h = \\sqrt{2}\\,Q$), conformal to the Einstein metric, with the same null rays and the "
        "horizon at $r = r_h$."]
    nr.CAPTIONS[(METRIC, MAGNETIC, "radial")] = [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) of the string metric of the magnetic hole "
        "($a = 1$, $r_0 = r_h = \\sqrt{2}\\,Q$), Witten's two dimensional black hole, each point a sphere of fixed "
        "area $8\\pi Q^2$."]
    nr.CAPTIONS[(METRIC, INVERSE, "radial")] = [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) with the potential $b = 1/a$ ($a = 1/2$, "
        "$M = 4r_0/3$, $q = r_0$), with an event horizon at $r = 3r_0/2$ and an inner horizon at $r = r_0/2$, "
        "between which $t$ is a space coordinate."]
    nr.CAPTIONS[(METRIC, MINUS_A, "radial")] = [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) with the potential $b = -a$ ($a = 1/2$, "
        "$r_0 = r_h$, $\\sigma = 2$), with one horizon at $r = r_h$ and the singularity $r = 0$ inside it."]
    nr.CAPTIONS[(METRIC, MINUS_INVERSE, "radial")] = [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) with the potential $b = -1/a$ ($a = 2$, "
        f"$r_0 = r_h$, $\\lambda = -1/(10r_h^2)$), with the event horizon at $r = {EVENT:.2f}\\,r_h$ and a "
        f"cosmological horizon at $r = {COSMOLOGICAL:.2f}\\,r_h$."]
    nr.CAPTIONS[(METRIC, MASSLESS, "plane")] = [
        "The plane of $\\eta$ and $\\xi$ ($\\theta = \\pi/2$, $\\phi = 0$) of the hole of no mass at $a = 1$, flat, "
        "each point a 2-sphere of area $\\pi(\\xi^2 - \\eta^2)$, which closes on the null singularity "
        "$\\xi = |\\eta|$."]


# -- The slices on the spacetime diagrams ------------------------------------------------------------------

def register_slices(sl):
    def first():
        return sl.one(METRIC, lambda m: sl.along(0.0, *m.reach(CHM, "r")))

    def areal():
        # The areal radius is sqrt(r_0 r) for the radius r of the first chart, r_0 = 1.
        return sl.one(METRIC, lambda m: sl.along(0.0, *(math.sqrt(r) for r in m.reach(CHM, "r"))))

    sl.FLAT[(METRIC, CHM, "radial")] = first
    sl.FLAT[(METRIC, AREAL, "radial")] = areal
    sl.FLAT_METRICS.add(METRIC)
    other = "another member of the family than the hole embedded at $a = 1$ and $r_0 = r_h$"
    for system, view in ((CHM, "third"), (INVERSE, "radial"), (MINUS_A, "radial"), (MINUS_INVERSE, "radial"),
                         (MASSLESS, "plane")):
        sl.HIDDEN[(METRIC, system, view)] = other
    for system in (ELECTRIC, MAGNETIC):
        sl.HIDDEN[(METRIC, system, "radial")] = (
            "the string metric, which measures the manifold another way than the Einstein metric embedded")


# -- Embedding diagram -------------------------------------------------------------------------------------

def register_embedding(eb):
    eb.DRAWN[METRIC] = lambda ck, src: embedding(eb, ck, src)
    eb.CAPTIONS[(METRIC, "einstein")] = [
        "The equatorial plane ($\\theta = \\pi/2$) of the linear dilaton black hole ($a = 1$, $r_0 = r_h$) at "
        "$t = 0$, with a throat at the horizon $r = r_h$ and a cone far away ($dz/d\\rho \\to \\sqrt{3}$)."]


def embedding(eb, ck, src):
    """The slice t = 0 of the hole at a = 1 and r_0 = r_h = 1, from the horizon out to r = 6, each exterior,
    joined at the bifurcation sphere r = 1, where rho = sqrt(r) is smallest; far away (dz/drho)^2 =
    (3r + 1)/(r - 1) runs to 3, the slope of a cone."""
    name = "Chan-Horne-Mann dilaton black hole"
    sl = eb.Slice(src, METRIC, CHM, "r", "\\phi", {"t": 0, **eb.EQUATOR}, ONE)
    ck.add(f"{name}: the horizon is at r_h", abs(sl.horizons()[0] - 1.0), 1e-12)
    ck.stops(f"{name}, inside the horizon", sl, np.linspace(0.0, 1.0, 402)[1:-1])
    size = 2 * TOP
    radii = (2.0, 3.0, 4.0, 5.0)
    near = eb.Piece("exterior", "sheet", sl, 1.0, TOP, 0.0, 1,
                    (("throat", "the throat $r = r_h$, the bifurcation sphere, where the other exterior begins"),
                     ("edge", "the surface runs on as $r \\to \\infty$, toward a cone")),
                    [(1.0, "horizon", "$r = r_h$")] + [(r, "r", None) for r in radii] + [(TOP, "r", None)], size)
    far = eb.Piece("other_exterior", "sheet2", sl, 1.0, TOP, 0.0, -1,
                   (("throat", "the throat $r = r_h$"), ("edge", "the other exterior, as $r \\to \\infty$")),
                   [(r, "r2", None) for r in radii] + [(TOP, "r2", None)], size)
    ck.join(f"{name}, the two exteriors at the throat", near, 1.0, far, 1.0)
    for p in (near, far):
        ck.isometry(f"{name}, {p.id}", p)
        ck.radius(f"{name}, {p.id}, rho = sqrt(r_0 r)", p, lambda r: np.sqrt(r), size)
    # The slope of the profile, dz/drho = sqrt((3r + 1)/(r - 1)), at points of the exterior.
    for r in (2.0, 4.0, TOP):
        i = int(np.argmin(np.abs(near.x - r)))
        lo, hi = max(i - 1, 0), min(i + 1, len(near.x) - 1)
        slope = (near.z[hi] - near.z[lo]) / (near.rho[hi] - near.rho[lo])
        ck.add(f"{name}: dz/drho = sqrt((3r + 1)/(r - 1)) near r = {r:g}",
               abs(slope - math.sqrt((3 * near.x[i] + 1) / (near.x[i] - 1))) / math.sqrt(3), 2e-2)
    surface = eb.Surface([near, far])
    fig = eb.figure_of([surface], {"sheet": "cover"}, size)
    eb.ring_label(fig, [0, 0, 0], *near.at(1.0), "$r = r_h$", dx=14)
    eb.ring_label(fig, [0, 0, 0], *near.at(3.0), "$3\\,r_h$")
    eb.ring_label(fig, [0, 0, 0], *near.at(TOP), "$6\\,r_h$")
    fig.legend("fill", "cover", "the exterior $r > r_h$")
    fig.legend("line", "r", "$r$ constant, $2$ to $6\\,r_h$")
    fig.legend("line", "r2", "the same radii on the other exterior")
    fig.legend("line", "horizon", "the throat $r = r_h$, on the horizon")
    fig.legend("line", "meridian", "$\\phi$ constant, every $15°$")
    stops = ["Inside the horizon, $0 < r < r_h$, $r$ is a time, and a slice of constant $t$ is not a moment of "
             "space."]
    return [eb.view("einstein", "Einstein metric", "$r_h$", [surface], fig.done(),
                    settings="$a = 1$ and $r_0 = r_h = 1$, the unit of every length.", stops=stops)]


# -- Conformal diagram -------------------------------------------------------------------------------------

def kruskal(cf, cell, t, rs):
    """(p, q) of points (t, r_*) of a cell, kappa = 1/2: U = -+exp(-u/2), V = +-exp(v/2)."""
    t = np.asarray(t, dtype=float)
    u, v = t - rs, t + rs
    p, q = cf.atan_exp(-u / 2), cf.atan_exp(v / 2)
    return {"I": (-p, q), "II": (p, q), "I'": (p, -q), "IV": (-p, -q)}[cell]


def register_conformal(cf):
    cf.DRAWN[METRIC] = lambda ck, src: conformal(cf, ck, src)
    cf.CAPTIONS[(METRIC, "one")] = [
        "The plane of $t$ and $r$ of the linear dilaton black hole ($a = 1$, $r_0 = r_h$), each point in the diagram "
        "a 2-sphere, with Schwarzschild's causal structure and a spacelike singularity at $r = 0$."]
    cf.CAPTIONS[(METRIC, "third")] = [
        "The plane of $t$ and $r$ at $a^2 = 1/3$ ($r_0 = r_h$), each point in the diagram a 2-sphere, with a "
        "spacelike singularity at $r = 0$ and a timelike line $r \\to \\infty$ that light reaches at a finite $t$."]


def conformal(cf, ck, src):
    """The plane of t and r of the first chart at a = 1 and at a^2 = 1/3, r_0 = r_h = 1, from the closed
    forms r_* of the module docstring, as Kruskal's square with kappa = 1/2 in each."""
    HALF, PI = cf.HALF, cf.PI
    views = []
    for vid, params, rstar, label in (("one", ONE, rstar_one, "$a = 1$"), ("third", THIRD, rstar_third, "$a^2 = 1/3$")):
        name = f"Chan-Horne-Mann {label}"
        plane = cf.Plane(src, METRIC, CHM, ("t", "r"), cf.EQUATOR, params)
        far = 30.0
        for cell, rs, future, where in (("I", (1.001, far), (1, 0), "exterior"), ("II", (0.01, 0.999), (0, -1), "black hole"),
                                        ("I'", (1.001, far), (-1, 0), "other exterior"),
                                        ("IV", (0.01, 0.999), (0, 1), "white hole")):
            ck.chart(f"{name}, {where}", plane, lambda t, r, c=cell, f=rstar: kruskal(cf, c, t, f(r)),
                     ck.uniform(-15, 15), ck.uniform(*rs), lambda t, r, d=future: d)
        rr = np.array([0.2, 0.5, 0.8, 1.5, 3.0, 9.0])
        inside, outside = rr < 1, rr > 1
        p, q = kruskal(cf, "I", 0.4 + 0 * rr[outside], rstar(rr[outside]))
        ck.limit(f"{name}: tan p tan q is UV = -exp(r_*) outside", np.tan(p) * np.tan(q),
                 -np.exp(rstar(rr[outside])), 1e-9)
        p, q = kruskal(cf, "II", 0.4 + 0 * rr[inside], rstar(rr[inside]))
        ck.limit(f"{name}: tan p tan q is UV = exp(r_*) inside", np.tan(p) * np.tan(q), np.exp(rstar(rr[inside])), 1e-9)
        p, q = kruskal(cf, "II", np.array([-3.0, 0.0, 3.0]), rstar(np.full(3, 1e-14)))
        ck.limit(f"{name}: r -> 0 lands on T = pi/2", p + q, [HALF] * 3, 1e-6)
        p, q = kruskal(cf, "I", np.array([3.0]), rstar(np.array([1 + 1e-12])))
        ck.limit(f"{name}: r -> r_h at fixed t lands on the bifurcation sphere", cf.point(p[0], q[0]), [0, 0], 1e-3)
        K = plane.kretschmann
        ck.diverges(f"{name}: the Kretschmann scalar diverges at r = 0", K(0, 1e-4), K(0, 1e-5))
        ck.finite(f"{name}: the Kretschmann scalar is finite at r = r_h", K(np.zeros(3), np.array([0.999, 1, 1.001])))
        moment = cf.slices.moments(METRIC, "einstein")[0] if vid == "one" else None
        if vid == "one":
            # r_* = ln(r - 1) grows without bound: r -> infinity at fixed t runs to i0.
            p, q = kruskal(cf, "I", np.array([0.0]), rstar(np.array([1e8])))
            ck.limit(f"{name}: r -> infinity at t = 0 lands on i0", cf.point(p[0], q[0]), [PI, 0], 1e-3)
            box = [-PI - 0.25, PI + 0.25, -HALF - 0.25, HALF + 0.25]
            outline = [[PI, 0], [HALF, HALF], [-HALF, HALF], [-PI, 0], [-HALF, -HALF], [HALF, -HALF]]
            exterior = [[0, 0], [HALF, -HALF], [PI, 0], [HALF, HALF]]
        else:
            # r_* tends to 0: r -> infinity at every t lands on X = pi/2, where UV = -1.
            ts = np.array([-3.0, 0.0, 3.0])
            p, q = kruskal(cf, "I", ts, rstar(np.full(3, 1e12)))
            ck.limit(f"{name}: r -> infinity lands on the timelike line X = pi/2", q - p, [HALF] * 3, 1e-5)
            ck.limit(f"{name}: r_* tends to 0 as r -> infinity", np.array([float(rstar(1e12))]), [0.0], 1e-5)
            box = [-HALF - 0.3, HALF + 0.3, -HALF - 0.3, HALF + 0.3]
            outline = [[HALF, HALF], [-HALF, HALF], [-HALF, -HALF], [HALF, -HALF]]
            exterior = [[0, 0], [HALF, -HALF], [HALF, HALF]]
        v = cf.View(vid, label, box, CHM)
        v.fill("region", outline)
        v.fill("cover", exterior)
        t = cf.S_ALL
        for r in (1.25, 1.5, 2.0, 3.0, 6.0):
            v.curve("r", *kruskal(cf, "I", t, rstar(np.full_like(t, r))))
        rr = cf.spread(1, np.inf, 600, 14)
        for tt in (-4, -2, -1, 0, 1, 2, 4):
            v.curve("t", *kruskal(cf, "I", np.full_like(rr, tt), rstar(rr)))
        v.line("horizon", [[[-HALF, -HALF], [HALF, HALF]], [[HALF, -HALF], [-HALF, HALF]]])
        v.line("singular", [[[-HALF, HALF], [HALF, HALF]], [[-HALF, -HALF], [HALF, -HALF]]], zig=True)
        if vid == "one":
            v.line("scri", [[[PI, 0], [HALF, HALF]], [[PI, 0], [HALF, -HALF]],
                            [[-PI, 0], [-HALF, HALF]], [[-PI, 0], [-HALF, -HALF]]])
            v.label_xt([PI, 0], "$i^0$", "l", dx=6)
            v.label_xt([-PI, 0], "$i^0$", "r", dx=-6)
            v.legend("scri", "null infinity, $r \\to \\infty$")
        else:
            v.line("boundary", [[[HALF, -HALF], [HALF, HALF]], [[-HALF, -HALF], [-HALF, HALF]]])
            v.label_xt([HALF, 0], "$r \\to \\infty$", "l", "small", dx=6)
            v.legend("boundary", "$r \\to \\infty$, a timelike line")
        v.label_xt([0, HALF], "$r = 0$", "b", dy=-8)
        v.label_xt([0, -HALF], "$r = 0$", "t", dy=8)
        v.label_xt([-HALF / 2, HALF / 2], "$r = r_h$", "tr", "small", dx=-6, dy=2)
        v.legend("cover", "the exterior $r > r_h$")
        v.legend("r", "$r$ constant")
        v.legend("t", "$t$ constant")
        v.legend("horizon", "the horizon $r = r_h$")
        v.legend("singular", "$r = 0$, a curvature singularity")
        if moment is not None:
            lo, hi = moment.reach(CHM, "r")
            r_m = np.linspace(lo + 1e-9, hi, 400)
            left, right = kruskal(cf, "I'", 0 * r_m, rstar(r_m[::-1])), kruskal(cf, "I", 0 * r_m, rstar(r_m))
            v.slice(moment, [(np.concatenate([left[0], right[0]]), np.concatenate([left[1], right[1]]))])
        v.set(settings=f"{label} and $r_0 = r_h = 1$, the unit of every length.")
        views.append(v)
    return views


# -- Curvature figure --------------------------------------------------------------------------------------

def register_curvature(cv):
    cv.DRAWN[METRIC] = lambda ck, src: curvature(cv, ck, src)


def curvature(cv, ck, src):
    """K and R against r from the stored scalars of the first chart at a = 1 and r_0 = r_h = 1, where
    h = r: K = (11r^2 + 6r + 3)/(4r^4), positive, finite at the horizon, growing as 3/(4r^4) toward r = 0
    and falling as 11/(4r^2) far away; R = (r - 1)/(2r^2), negative inside the horizon, zero on it, and
    falling as 1/(2r) far away. Both fall more slowly than Schwarzschild's K = 12 r_s^2/r^6."""
    params = {"a": 1.0, "r_0": 1.0, "r_h": 1.0}
    names = list(params)
    Kp = cv.Published(src, METRIC, CHM, "kretschmann", ["r"], names)
    Rp = cv.Published(src, METRIC, CHM, "ricci_scalar", ["r"], names)
    x_axis = cv.axis("$r/r_h$", 0.1, 6)
    K_scale = {"kind": "log", "from": 1e-2, "to": 1e4}
    R_scale = {"kind": "symlog", "from": -1e2, "to": 1, "linear": 1e-2}
    panels = [{"id": "K", "name": "$K\\,r_h^4$", "aspect": 2, "scale": K_scale, "ticks": cv.value_ticks(K_scale)},
              {"id": "R", "name": "$R\\,r_h^2$", "aspect": 2, "scale": R_scale, "ticks": cv.value_ticks(R_scale)}]
    cv.shape_checks(ck, METRIC, [p["aspect"] for p in panels])
    height = {p["id"]: cv.PLOT_W / p["aspect"] for p in panels}
    curves = (cv.line_curves(ck, METRIC, "K", "K", Kp, "r", x_axis, K_scale, height["K"], params)
              + cv.line_curves(ck, METRIC, "R", "R", Rp, "r", x_axis, R_scale, height["R"], params))
    vals = list(params.values())
    rs = np.linspace(0.1, 6, cv.SAMPLES)
    name = "Chan-Horne-Mann dilaton black hole at a = 1"
    ck(f"{name}: K = (11r^2 + 6r r_h + 3r_h^2)/(4r_0^2 r^4)",
       bool(np.allclose(Kp(rs, *vals), (11 * rs ** 2 + 6 * rs + 3) / (4 * rs ** 4), rtol=1e-12)))
    ck(f"{name}: R = (r - r_h)/(2r_0 r^2)", bool(np.allclose(Rp(rs, *vals), (rs - 1) / (2 * rs ** 2), rtol=1e-12)))
    ck(f"{name}: K is positive", bool(np.all(Kp(rs, *vals) > 0)))
    ck(f"{name}: R vanishes at the horizon", abs(float(Rp.exact([1.0] + vals))) < 1e-12)
    ck(f"{name}: K is finite at the horizon", bool(np.isfinite(Kp(np.array([1.0]), *vals)).all()))
    far = np.array([1e6])
    ck(f"{name}: K falls as 11/(4r^2)", abs(float(Kp(far, *vals)[0]) * 4e12 / 11 - 1) < 1e-5)
    ck(f"{name}: R falls as 1/(2r)", abs(float(Rp(far, *vals)[0]) * 2e6 - 1) < 1e-5)
    legend = [["line", "K", "$K$"], ["line", "R", "$R$"], ["line", "horizon", "$r = r_h$"]]
    return [{
        "id": "radial", "label": "Along $r$", "kind": "line", "x": x_axis,
        "panels": panels, "curves": curves, "marks": [{"class": "horizon", "x": 1.0}],
        "legend": legend,
        "caption": ["The Kretschmann scalar $K$ and the Ricci scalar $R$ of the linear dilaton black hole against $r$ "
                    "($a = 1$, $r_0 = r_h$). "
                    "$K$ is finite at the horizon, diverges at $r = 0$, and falls as $11/(4r_0^2r^2)$. "
                    "$R$ vanishes at the horizon and falls as $1/(2r_0r)$."],
        "ricci_sentences": [[0, 2]],
        "per_scalar": cv.scalar_clips("A logarithmic scale from $10^{-2}$ to $10^4$.",
                                      "A signed logarithmic scale from $-10^2$ to $1$, linear within $10^{-2}$."),
        "settings": "$a = 1$, $r_0 = r_h = 1$.",
        "clip": "$K$ logarithmic, $10^{-2}$ to $10^4$; $R$ signed logarithmic, $-10^2$ to $1$, linear within $10^{-2}$.",
    }]
