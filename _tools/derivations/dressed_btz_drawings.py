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
