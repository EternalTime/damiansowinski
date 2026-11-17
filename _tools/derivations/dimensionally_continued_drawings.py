"""Drawings of the dimensionally continued black holes, each with its horizon at r_+ = l, the unit of length.

Five dimensions, M = 3: sqrt(M + 1) = 2, so f = 1 - 2 + r^2 = r^2 - 1, the metric function of the black
hole of three dimensions without rotation at its M = 1. The surface gravity is k = f'(1)/2 = 1 and

    r_* = (1/2) ln|(r - 1)/(r + 1)|

vanishes at r = 0 and as r -> infinity, so Kruskal's U = -exp(-k u), V = exp(k v), u, v = ct -+ r_*, have
UV = (1 - r)/(1 + r): 1 at r = 0, 0 on the horizon and -1 on the conformal boundary, and with p = arctan U
and q = arctan V the spacetime is a square, r = 0 on T = +-pi/2 and the boundary on X = +-pi/2. In the
conformal chart r = sec(rho), so r_* = ln tan(rho/2).

Six dimensions, M = 2: f = 1 - 2/sqrt(r) + r^2, which vanishes at r = 1 with f'(1) = 3, k = 3/2. With
x = sqrt(r), dr/f = 2x^2 dx/(x^5 + x - 2), a rational function with simple poles at the five roots x_i of
x^5 + x - 2 = (x - 1)(x^4 + x^3 + x^2 + x + 2), so

    r_* = Re sum_i A_i ln(1 - sqrt(r)/x_i),      A_i = 2 x_i^2/(5 x_i^4 + 1),

which vanishes at r = 0 and tends to the finite R = Re sum_i A_i ln(-1/x_i) as r -> infinity, since the
A_i sum to zero. The residue at x = 1 is 1/3 = 1/2k, so UV = (1 - sqrt(r)) exp(2k Re sum' A_i ln(1 -
sqrt(r)/x_i)), the sum over the other four roots, is analytic through the horizon: 1 at r = 0 and -B,
B = exp(2kR), on the boundary, which lies on tan p tan q = -B as Schwarzschild-anti-de Sitter's does.

The slice t = 0 on the plane of r and phi has g_rr = 1/f, at least 1 from the horizon out to the circle
where f = 1, r = sqrt 2 in five dimensions and r = 2^(2/5) in six, and less than 1 beyond it: it is drawn in
flat space from the throat to that circle and in three dimensional Minkowski space beyond, through both
exteriors, as the black hole of three dimensions and Schwarzschild-anti-de Sitter are.
"""
import math

import numpy as np

METRIC = "dimensionally_continued"
STATIC5, IN5, CONF5 = "spherical", "eddington_finkelstein_ingoing", "conformal"
STATIC6, IN6 = "spherical_six", "eddington_finkelstein_ingoing_six"
P5, P6 = {"M": 3, "ell": 1}, {"M": 2, "ell": 1}
FIXED5 = {"psi": "pi/2", "theta": "pi/2", "phi": "0"}
FIXED6 = {"chi": "pi/2", **FIXED5}
ANGLES5 = "$\\psi = \\theta = \\pi/2$, $\\phi = 0$"
ANGLES6 = "$\\chi = \\psi = \\theta = \\pi/2$, $\\phi = 0$"
SETTINGS5 = "$M = 3$ and $\\ell = 1$, so that $r_+ = \\ell$."
SETTINGS6 = "$M = 2\\ell$ and $\\ell = 1$, so that $r_+ = \\ell$."
KAPPA5, KAPPA6 = 1.0, 1.5
LEVEL5, LEVEL6 = math.sqrt(2.0), 2.0 ** 0.4
TOP = 3.0

ROOTS6 = np.roots([1, 0, 0, 0, 1, -2])
RESIDUES6 = 2 * ROOTS6 ** 2 / (5 * ROOTS6 ** 4 + 1)
HORIZON6 = int(np.argmin(np.abs(ROOTS6 - 1)))
FAR6 = float(np.sum(RESIDUES6 * np.log(-1 / ROOTS6)).real)
B6 = math.exp(2 * KAPPA6 * FAR6)


def rstar5(r):
    """The tortoise coordinate of five dimensions at M = 3, zero at r = 0 and as r -> infinity."""
    r = np.asarray(r, dtype=float)
    with np.errstate(divide="ignore"):
        return 0.5 * np.log(np.abs((r - 1) / (r + 1)))


def rstar6(r):
    """The tortoise coordinate of six dimensions at M = 2, zero at r = 0."""
    x = np.sqrt(np.asarray(r, dtype=float)).astype(complex)
    with np.errstate(divide="ignore"):
        return sum(A * np.log(1 - x / xi) for A, xi in zip(RESIDUES6, ROOTS6)).real


def product6(r):
    """Kruskal's UV in six dimensions, one formula for every r > 0."""
    x = np.sqrt(np.asarray(r, dtype=float)).astype(complex)
    rest = sum(A * np.log(1 - x / xi) for i, (A, xi) in enumerate(zip(RESIDUES6, ROOTS6)) if i != HORIZON6).real
    return (1 - x.real) * np.exp(2 * KAPPA6 * rest)


def rho_of(r):
    """The conformal chart's rho at the areal radius r, r = r_+ sec(rho) with r_+ = l = 1."""
    return np.arccos(1 / np.asarray(r, dtype=float))


def rstar_conformal(rho):
    """The tortoise coordinate of five dimensions in the conformal chart, ln tan(rho/2), which is
    r_*(sec rho) and vanishes as rho -> pi/2."""
    with np.errstate(divide="ignore"):
        return np.log(np.tan(np.asarray(rho, dtype=float) / 2))


class Hole:
    """The plane of t and r of one black hole, drawn as Kruskal and Szekeres draw Schwarzschild's, in the
    cells of conformal.Tower: G(u) = arctan exp(-k u) with u, v = ct -+ r_*, and

        I    exterior, right       p = -G(u),   q = G(-v)
        I'   exterior, left        p = G(u),    q = -G(-v)
        II   black hole            p = G(u),    q = G(-v)
        IV   white hole            p = -G(u),   q = -G(-v)

    so that tan p tan q is Kruskal's UV, 1 at r = 0 and -B on the conformal boundary."""

    def __init__(self, cf, rstar, product, kappa, far):
        self.cf, self.rstar, self.UV, self.kp, self.far = cf, rstar, product, kappa, far
        self.B = math.exp(2 * kappa * far)

    def G(self, u):
        return self.cf.atan_exp(-self.kp * np.asarray(u, dtype=float))

    def pq(self, cell, t, r):
        rs = self.rstar(r)
        t = np.asarray(t, dtype=float)
        u, v = t - rs, t + rs
        return {"I": (-self.G(u), self.G(-v)), "I'": (self.G(u), -self.G(-v)),
                "II": (self.G(u), self.G(-v)), "IV": (-self.G(u), -self.G(-v))}[cell]

    def ingoing(self, w, r):
        """(p, q) of the ingoing chart's (v, r): V = exp(k v) and U = UV(r)/V, one formula for every r > 0."""
        w = np.asarray(w, dtype=float)
        return np.arctan(self.UV(r) * np.exp(-self.kp * w)), self.cf.atan_exp(self.kp * w)

    def boundary(self, n=400):
        """r -> infinity in cell I: tan p tan q = -B, from the corner (X, T) = (pi/2, -pi/2) to (pi/2, pi/2)."""
        HALF = math.pi / 2
        p = -HALF / (1 + np.exp(np.linspace(-14, 14, n)))
        return p, np.arctan(-self.B / np.tan(p))


def product5(r):
    """Kruskal's UV in five dimensions at M = 3, (1 - r)/(1 + r)."""
    r = np.asarray(r, dtype=float)
    return (1 - r) / (1 + r)


# -- Conformal diagrams ----------------------------------------------------------------------------------

def register_conformal(cf):
    cf.DRAWN[METRIC] = lambda ck, src: conformal(cf, ck, src)
    for vid, d, sphere in (("static", "five", "3-sphere"), ("static_six", "six", "4-sphere")):
        mass = "$M = 3$" if d == "five" else "$M = 2\\ell$"
        cf.CAPTIONS[(METRIC, vid)] = [
            f"The black hole of {d} dimensions ({mass}), maximally extended, each point in the diagram a {sphere}, "
            "with the spacelike singularity $r = 0$ inside the horizon."]
    for vid, d, sphere in (("ingoing", "five", "3-sphere"), ("ingoing_six", "six", "4-sphere")):
        mass = "$M = 3$" if d == "five" else "$M = 2\\ell$"
        cf.CAPTIONS[(METRIC, vid)] = [
            f"The black hole of {d} dimensions ({mass}), each point in the diagram a {sphere}, with the ingoing "
            "rays of constant $v$ crossing the horizon to the singularity $r = 0$."]
    cf.CAPTIONS[(METRIC, "conformal")] = [
        "The black hole of five dimensions ($M = 3$), each point in the diagram a 3-sphere, with $t$ and $\\rho$ "
        "covering one exterior from the horizon $\\rho = 0$ to the boundary $\\rho = \\pi/2$."]


def conformal(cf, ck, src):
    """Each black hole maximally extended in the cells of Hole, from the closed forms of the module
    docstring, checked under every chart that covers it. In five dimensions B = 1 and the diagram is the
    square of the black hole of three dimensions, its conformal boundary on X = +-pi/2; in six B =
    exp(2kR) = 4.79 and the boundary lies on tan p tan q = -B, bowed outward to X = 2 arctan sqrt(B) =
    2.28 on T = 0, as Schwarzschild-anti-de Sitter's is."""
    HALF, Q4 = math.pi / 2, math.pi / 4
    views = []
    for d, static_id, in_id, H, params, fixed, settings, sphere in (
            ("five", STATIC5, IN5, Hole(cf, rstar5, product5, KAPPA5, 0.0), P5, FIXED5, SETTINGS5, "3-sphere"),
            ("six", STATIC6, IN6, Hole(cf, rstar6, product6, KAPPA6, FAR6), P6, FIXED6, SETTINGS6, "4-sphere")):
        name = f"Dimensionally continued, {d}"
        st = cf.Plane(src, METRIC, static_id, ("t", "r"), fixed, params)
        ein = cf.Plane(src, METRIC, in_id, ("v", "r"), fixed, params)
        for cell, rs, future, what in (("I", (1.001, 40), (1, 0), "exterior"),
                                       ("II", (0.01, 0.999), (0, -1), "black hole"),
                                       ("IV", (0.01, 0.999), (0, 1), "white hole"),
                                       ("I'", (1.001, 40), (-1, 0), "other exterior")):
            ck.chart(f"{name}, static, {what}", st, lambda t, r, c=cell: H.pq(c, t, r),
                     ck.uniform(-4, 4), ck.uniform(*rs), lambda t, r, f=future: f)
        f_of = (lambda r: 1 - 2.0 + r * r) if d == "five" else (lambda r: 1 - 2 / np.sqrt(r) + r * r)
        ck.chart(f"{name}, ingoing Eddington-Finkelstein", ein, H.ingoing, ck.uniform(-4, 4), ck.uniform(0.01, 40),
                 lambda w, r: (np.ones_like(r), -np.maximum(1, np.abs(f_of(r)))))
        rr = np.linspace(0.05, 5, 50)
        for cell, sel in (("I", rr > 1.001), ("II", rr < 0.999)):
            p, q = H.pq(cell, 0.3 + 0 * rr[sel], rr[sel])
            ck.limit(f"{name}, {cell}: tan p tan q is Kruskal's UV", np.tan(p) * np.tan(q), H.UV(rr[sel]), 1e-8)
        ck.limit(f"{name}: r_* vanishes at r = 0", H.rstar(np.array([1e-12])), [0.0], 1e-6)
        ck.limit(f"{name}: r_* tends to R as r -> infinity", H.rstar(np.array([1e9])), [H.far], 1e-4)
        ck.limit(f"{name}: UV is 1 at r = 0 and 0 on the horizon", H.UV(np.array([1e-12, 1.0])), [1.0, 0.0], 1e-5)
        p, q = H.pq("II", np.array([-1.0, 0, 1]), np.full(3, 1e-12))
        ck.limit(f"{name}: r -> 0 in the black hole lands on T = pi/2", p + q, [HALF] * 3, 1e-5)
        p, q = H.pq("IV", np.array([-1.0, 0, 1]), np.full(3, 1e-12))
        ck.limit(f"{name}: r -> 0 in the white hole lands on T = -pi/2", p + q, [-HALF] * 3, 1e-5)
        p, q = H.pq("I", np.array([-0.5, 0, 0.5]), np.full(3, 1e12))
        ck.limit(f"{name}: r -> infinity lands on the boundary tan p tan q = -B", np.tan(p) * np.tan(q),
                 [-H.B] * 3, 1e-4)
        p, q = H.pq("I", np.array([0.5]), np.array([1 + 1e-12]))
        ck.limit(f"{name}: r -> r_+ at fixed t lands on the bifurcation sphere", cf.point(p[0], q[0]), [0, 0], 1e-4)
        ck.limit(f"{name}: the ingoing and static coordinates put one event at one point",
                 np.array(H.ingoing(0.4 + H.rstar(3.0), 3.0)).ravel(), np.array(H.pq("I", 0.4, 3.0)).ravel(), 1e-12)
        ck.limit(f"{name}: the ingoing and static coordinates put one event at one point in the black hole",
                 np.array(H.ingoing(0.4 + H.rstar(0.5), 0.5)).ravel(), np.array(H.pq("II", 0.4, 0.5)).ravel(),
                 1e-12)
        K = st.kretschmann
        ck.diverges(f"{name}: the Kretschmann scalar diverges at r = 0", K(0, 1e-2), K(0, 1e-3))
        ck.finite(f"{name}: the Kretschmann scalar is finite at r = r_+", K(np.zeros(3), np.array([0.999, 1, 1.001])))
        if d == "five":
            ck.limit(f"{name}: B = 1, so the boundary is the vertical line X = pi/2", H.B, 1.0, 1e-15)
            cv = cf.Plane(src, METRIC, CONF5, ("t", "\\rho"), fixed, params)
            ck.chart(f"{name}, conformal", cv, lambda t, rho: H.pq("I", t, 1 / np.cos(rho)),
                     ck.uniform(-4, 4), ck.uniform(0.01, HALF - 0.01), lambda t, rho: (1, 0))
            ck.limit(f"{name}: ln tan(rho/2) is r_* at r = sec(rho)", rstar_conformal(np.array([0.3, 0.9, 1.4])),
                     rstar5(1 / np.cos(np.array([0.3, 0.9, 1.4]))), 1e-12)
        else:
            ck.limit(f"{name}: R = 0.522 l and B = exp(2kR) = 4.79", [H.far, H.B], [0.5223, 4.792], 1e-3)
            ck.limit(f"{name}: the residue of r_* at the horizon is 1/2k", RESIDUES6[HORIZON6].real, 1 / (2 * KAPPA6),
                     1e-12)

        reach = 2 * math.atan(math.sqrt(H.B))
        box = [-reach - 0.3, reach + 0.3, -HALF - 0.25, HALF + 0.25]
        if d == "five":
            right = [[HALF, -HALF], [HALF, HALF]]
            named = [HALF, 1.0]
        else:
            bp, bq = H.boundary()
            right = [cf.point(a, b) for a, b in zip(bp, bq)]
            named = min(right, key=lambda at: abs(at[1] - 1.2))
        left = [[-x, t_] for x, t_ in right[::-1]]
        whole = right + left
        t_all = cf.S_ALL
        r_out, r_in, times = (1.2, 1.5, 2.0, 4.0), (0.3, 0.6, 0.9), (-1.0, -0.5, 0.0, 0.5, 1.0)
        beyond = cf.spread(1.0, np.inf, 500, 14)
        inside = np.linspace(1e-6, 1 - 1e-9, 500)
        side = 1.05 if d == "five" else 1.6

        def edges(v):
            if d == "five":
                v.line("boundary", [right, [[-HALF, -HALF], [-HALF, HALF]]])
            else:
                v.curve("boundary", bp, bq)
                v.curve("boundary", -bq, -bp)
            v.line("horizon", [[[-HALF, -HALF], [HALF, HALF]], [[HALF, -HALF], [-HALF, HALF]]])
            v.line("singular", [[[-HALF, HALF], [HALF, HALF]], [[-HALF, -HALF], [HALF, -HALF]]], zig=True)
            v.label_xt([0, HALF], "$r = 0$", "b", dy=-8)
            v.label_xt([0, -HALF], "$r = 0$", "t", dy=8)
            v.label_xt(named, "$r \\to \\infty$", "l", "small", dx=6)
            v.label_xt([-named[0], named[1]], "$r \\to \\infty$", "r", "small", dx=-6)
            v.label_xt([-Q4, Q4], "$r_+$", "tr", "small", dx=-6, dy=2)
            v.label_xt([side, -0.45], "exterior", cls="region")
            v.label_xt([-side, -0.45], "exterior", cls="region")
            v.label_xt([0, 1.15], "black hole", cls="region")
            v.label_xt([0, -1.15], "white hole", cls="region")
            v.legend("horizon", "the horizon $r = r_+$")
            v.legend("boundary", "the conformal boundary $r \\to \\infty$, timelike")
            v.legend("singular", "$r = 0$, a curvature singularity")

        moment = cf.slices.moments(METRIC, d)[0]
        lo, hi = moment.reach(static_id, "r")
        through = cf.through_bifurcation(H, ("I'", "I"), hi, lo)
        suffix = "" if d == "five" else "_six"

        v = cf.View("static" + suffix, "Static", box, static_id)
        v.fill("region", whole)
        v.fill("cover", [[0, 0]] + right + [[-HALF, HALF]])
        for r in r_out:
            v.curve("r", *H.pq("I", t_all, np.full_like(t_all, r)))
        for r in r_in:
            v.curve("r", *H.pq("II", t_all, np.full_like(t_all, r)))
        for tt in times:
            v.curve("t", *H.pq("I", np.full_like(beyond, tt), beyond))
            v.curve("t", *H.pq("II", np.full_like(inside, tt), inside))
        edges(v)
        v.legend("cover", "what $t$ and $r$ cover")
        v.legend("r", "$r$ constant")
        v.legend("t", "$t$ constant")
        views.append(v)

        v = cf.View("ingoing" + suffix, "Ingoing Eddington-Finkelstein", box, in_id)
        v.fill("region", whole)
        v.fill("cover", [[0, 0]] + right + [[-HALF, HALF]])
        for r in r_out + r_in:
            v.curve("r", *H.ingoing(t_all, np.full_like(t_all, r)))
        radii = cf.spread(0, np.inf, 600, 14)
        for w in (-1.0, -0.5, 0.0, 0.5, 1.0):
            v.curve("null", *H.ingoing(np.full_like(radii, w), radii))
        edges(v)
        v.legend("cover", "what $v$ and $r$ cover")
        v.legend("r", "$r$ constant")
        v.legend("null", "$v$ constant, an ingoing light ray")
        views.append(v)

        if d == "five":
            v = cf.View("conformal", "Conformal", box, CONF5)
            v.fill("region", whole)
            v.fill("cover", [[0, 0]] + right)
            rhos = np.linspace(1e-6, HALF - 1e-9, 500)
            for rho in (math.pi / 8, math.pi / 4, 3 * math.pi / 8):
                v.curve("r", *H.pq("I", t_all, np.full_like(t_all, 1 / math.cos(rho))))
            for tt in times:
                v.curve("t", *H.pq("I", np.full_like(rhos, tt), 1 / np.cos(rhos)))
            edges(v)
            v.legend("cover", "what $t$ and $\\rho$ cover")
            v.legend("r", "$\\rho$ constant, every $\\pi/8$")
            v.legend("t", "$t$ constant")
            views.append(v)

        for view in views[-3 if d == "five" else -2:]:
            view.slice(moment, [through])
            view.set(settings=settings)
    return views


# -- Spacetime diagrams ----------------------------------------------------------------------------

def register_null(nr):
    """One plane of the time and the radius for each chart: the static charts through the horizon to
    r = 0, oriented by the ingoing rays as the ingoing charts are, the ingoing charts against v - r, and
    the conformal chart of five dimensions outside the horizon, against its own rho."""
    rows = (
        nr.Diagram(METRIC, STATIC5, "radial", "$t$ and $r$", ("t", "r"), (0, 3, -1.5, 1.5), "$r/\\ell$",
                   "$ct/\\ell$", P5, FIXED5, orient="ingoing", areal=True),
        nr.Diagram(METRIC, IN5, "finkelstein", "against $v - r$", ("v", "r"), (0, 3, -3, 0), "$r/\\ell$",
                   "$(v - r)/\\ell$", P5, FIXED5, to_display=nr.FINKELSTEIN_IN, orient="ingoing", areal=True),
        nr.Diagram(METRIC, CONF5, "radial", "$t$ and $\\rho$", ("t", "\\rho"), (0, 1.5, -0.75, 0.75), "$\\rho$",
                   "$ct/\\ell$", P5, FIXED5, areal=True),
        nr.Diagram(METRIC, STATIC6, "radial", "$t$ and $r$", ("t", "r"), (0, 3, -1.5, 1.5), "$r/\\ell$",
                   "$ct/\\ell$", P6, FIXED6, orient="ingoing", areal=True),
        nr.Diagram(METRIC, IN6, "finkelstein", "against $v - r$", ("v", "r"), (0, 3, -3, 0), "$r/\\ell$",
                   "$(v - r)/\\ell$", P6, FIXED6, to_display=nr.FINKELSTEIN_IN, orient="ingoing", areal=True),
    )
    nr.DIAGRAMS.extend(rows)
    nr.CAPTIONS[(METRIC, STATIC5, "radial")] = [
        f"The plane of $t$ and $r$ ({ANGLES5}) in five dimensions, its cones closing at the horizon $r_+ = \\ell$ "
        "and pointing to the singularity $r = 0$ inside it."]
    nr.CAPTIONS[(METRIC, IN5, "finkelstein")] = [
        f"The plane of $v$ and $r$ ({ANGLES5}) in five dimensions against $v - r$, the ingoing rays crossing the "
        "horizon $r_+ = \\ell$ on to the singularity $r = 0$."]
    nr.CAPTIONS[(METRIC, CONF5, "radial")] = [
        f"The plane of $t$ and $\\rho$ ({ANGLES5}) in five dimensions outside the horizon, its cones closing toward "
        "the horizon $\\rho = 0$, with $\\rho$ in radians."]
    nr.CAPTIONS[(METRIC, STATIC6, "radial")] = [
        f"The plane of $t$ and $r$ ({ANGLES6}) in six dimensions, its cones closing at the horizon $r_+ = \\ell$ "
        "and pointing to the singularity $r = 0$ inside it."]
    nr.CAPTIONS[(METRIC, IN6, "finkelstein")] = [
        f"The plane of $v$ and $r$ ({ANGLES6}) in six dimensions against $v - r$, the ingoing rays crossing the "
        "horizon $r_+ = \\ell$ on to the singularity $r = 0$."]

    def own(f):
        return lambda w, r: f(np.asarray(w, dtype=float), np.asarray(r, dtype=float))

    away = own(lambda w, r: np.abs(r - 1) > 0.05)
    for static, ingoing, rstar in ((STATIC5, IN5, rstar5), (STATIC6, IN6, rstar6)):
        nr.CLOSED_FORMS[(METRIC, static, "radial")] = (
            own(lambda t, r, rs=rstar: t + rs(r)), own(lambda t, r, rs=rstar: t - rs(r)), away)
        nr.CLOSED_FORMS[(METRIC, ingoing, "finkelstein")] = (
            own(lambda v, r: v), own(lambda v, r, rs=rstar: v - 2 * rs(r)), away)
    nr.CLOSED_FORMS[(METRIC, CONF5, "radial")] = (
        own(lambda t, rho: t + rstar_conformal(rho)), own(lambda t, rho: t - rstar_conformal(rho)),
        own(lambda t, rho: rho > 0.02))


# -- The slices on the spacetime diagrams --------------------------------------------------------------

def register_slices(sl):
    """The moment t = 0 of each embedding view on the planes of its dimension: along r in the static
    chart, along rho = arcsec(r/r_+) in the conformal one, and as v = r_* in the ingoing one, crowding
    toward the horizon, where the curve runs off. Each chart covers one exterior, and the moment's other
    exterior lies over the same r."""
    def static(view, system):
        m = sl.moments(METRIC, view)[0]
        return [sl.Mark(m, sl.along(0.0, *m.reach(system, "r")))]

    def ingoing(view, system, rstar):
        m = sl.moments(METRIC, view)[0]
        lo, hi = m.reach(system, "r")
        r = sl.near(lo, hi)
        return [sl.Mark(m, [np.column_stack([rstar(r), r])])]

    def conformal():
        m = sl.moments(METRIC, "five")[0]
        lo, hi = m.reach(STATIC5, "r")
        return [sl.Mark(m, sl.along(0.0, float(rho_of(lo)), float(rho_of(hi))))]

    sl.FLAT[(METRIC, STATIC5, "radial")] = lambda: static("five", STATIC5)
    sl.FLAT[(METRIC, IN5, "finkelstein")] = lambda: ingoing("five", STATIC5, rstar5)
    sl.FLAT[(METRIC, CONF5, "radial")] = conformal
    sl.FLAT[(METRIC, STATIC6, "radial")] = lambda: static("six", STATIC6)
    sl.FLAT[(METRIC, IN6, "finkelstein")] = lambda: ingoing("six", STATIC6, rstar6)
    sl.FLAT_METRICS.add(METRIC)


# -- Embedding diagrams ----------------------------------------------------------------------------------

def register_embedding(eb):
    eb.DRAWN[METRIC] = lambda ck, src: embedding(eb, ck, src)
    eb.CAPTIONS[(METRIC, "five")] = [
        "The slice $t = 0$ ($\\psi = \\theta = \\pi/2$) in five dimensions through both exteriors, in flat space "
        "to $r = \\sqrt{2}\\,\\ell$ and in Minkowski space beyond."]
    eb.CAPTIONS[(METRIC, "six")] = [
        "The slice $t = 0$ ($\\chi = \\psi = \\theta = \\pi/2$) in six dimensions through both exteriors, in "
        "flat space to $r = 2^{2/5}\\ell$ and in Minkowski space beyond."]


def embedding(eb, ck, src):
    """The slice t = 0 of each black hole through its throat, in flat space out to the circle where the
    metric function is 1 and in three dimensional Minkowski space beyond it, to r = 3."""
    views = []
    for vid, label, system, params, held, level, level_tex, near_r, settings in (
            ("five", "Five dimensions", STATIC5, P5, {"psi": "pi/2", "theta": "pi/2"}, LEVEL5, "\\sqrt{2}\\,\\ell",
             1.2, SETTINGS5),
            ("six", "Six dimensions", STATIC6, P6, {"chi": "pi/2", "psi": "pi/2", "theta": "pi/2"}, LEVEL6,
             "2^{2/5}\\ell", 1.15, SETTINGS6)):
        name = f"Dimensionally continued, {vid}"
        fixed_at = {"t": 0, **held}
        sl = eb.Slice(src, METRIC, system, "r", "\\phi", fixed_at, params)
        msl = eb.Slice(src, METRIC, system, "r", "\\phi", fixed_at, params, space="minkowski")
        # In six dimensions 1/g_rr holds sqrt(1/r), which no polynomial in r carries, so its zero is found
        # in numbers and remembered, as a root horizons() finds is.
        rp = sl.horizons()[0] if vid == "five" else sl.horizon_between(0.5, 1.5)
        size = 2 * TOP
        ck.add(f"{name}: the horizon is at r_+ = l", abs(rp - 1.0), 1e-12)
        ck.add(f"{name}: the surface lies level where the metric function is 1",
               abs(float(sl.defect_at(np.array([level]))[0])), 1e-12)
        ck.stops(f"{name}, beyond the level circle in flat space", sl, np.linspace(level, 40, 402)[1:])
        inward = np.linspace(rp, level, 402)[1:-1]
        ck.add(f"{name}: between r_+ and the level circle no surface in Minkowski space carries the slice",
               float(max(0.0, np.max(-msl.defect_at(inward)))), 0.0)
        if not np.all(msl.defect_at(inward) > 0):
            ck.items[-1]["ok"] = False
        join = (f"at $r = {level_tex}$ the surface lies level, in flat space nearer the horizon and in Minkowski "
                "space beyond")
        edge = "the sheet runs on toward a light cone, to $r \\to \\infty$"
        near = eb.Piece("exterior", "sheet", sl, rp, level, 0.0, 1,
                        (("throat", "the throat $r = r_+$, where the other exterior begins"), ("join", join)),
                        [(rp, "horizon", "$r = r_+$"), (near_r, "r", None), (level, "space", None)], size)
        out = eb.Piece("exterior_minkowski", "sheet", msl, level, TOP, near.at(level)[1], 1,
                       (("join", join), ("edge", edge)), [(2.0, "r", None), (TOP, "r", None)], size)
        far = eb.Piece("other_exterior", "sheet2", sl, rp, level, 0.0, -1,
                       (("throat", "the throat $r = r_+$"), ("join", join)),
                       [(near_r, "r2", None), (level, "space", None)], size)
        far_out = eb.Piece("other_exterior_minkowski", "sheet2", msl, level, TOP, far.at(level)[1], -1,
                           (("join", join), ("edge", edge)), [(2.0, "r2", None), (TOP, "r2", None)], size)
        for p in (near, out, far, far_out):
            space = "in Minkowski space" if p.sl.lorentz else "in flat space"
            ck.isometry(f"{name}, {p.id} {space}", p)
            ck.radius(f"{name}, {p.id}, rho = r {space}", p, lambda r: r, size)
        ck.join(f"{name}, the two sheets at the throat", near, rp, far, rp)
        for a, b in ((near, out), (far, far_out)):
            ck.join(f"{name}, {a.id} in flat space and {b.id} in Minkowski space at the level circle", a, level, b, level)
            pa, pb = a.data()["points"][-1], b.data()["points"][0]
            ck.add(f"{name}, {a.id} and {b.id} as written: one point at the level circle",
                   max(abs(pa[0] - pb[0]), abs(pa[1] - pb[1]), abs(pa[2] - pb[2])), 10.0 ** -min(a.decimals, b.decimals))
        surface = eb.Surface([near, out, far, far_out])
        fig = eb.figure_of([surface], {"sheet": "cover"}, size)
        eb.ring_label(fig, [0, 0, 0], rp, 0.0, "$r = r_+$", dx=14)
        eb.ring_label(fig, [0, 0, 0], *near.at(level), f"${level_tex}$")
        eb.ring_label(fig, [0, 0, 0], *out.at(TOP), "$3\\ell$")
        fig.legend("fill", "cover", "the exterior $r > r_+$")
        fig.legend("line", "r", f"$r = {near_r}\\,\\ell$, $2\\ell$, $3\\ell$")
        fig.legend("line", "r2", "the same radii on the other exterior")
        fig.legend("line", "horizon", "the horizon $r = r_+$")
        fig.legend("line", "space", f"$r = {level_tex}$, where $g_{{rr}} = 1$")
        fig.legend("line", "meridian", "$\\phi$ constant, every $15°$")
        views.append(eb.view(vid, label, "$\\ell$", [surface], fig.done(), settings=settings))
    return views


# -- Curvature figure -----------------------------------------------------------------------------------

def register_curvature(cv):
    cv.DRAWN[METRIC] = lambda ck, src: curvature(cv, ck, src)


def curvature(cv, ck, src):
    """K and R against r from the stored scalars of each static chart, at the values every other drawing
    takes. In five dimensions, with s = sqrt(M + 1), K = 40/l^4 - 24s/l^2r^2 + 12s^2/r^4 is least, 28/l^4
    at every M, at r^2 = s l^2, r = sqrt(2) l at M = 3, and R = -20/l^2 + 6s/r^2 changes sign at
    r^2 = 3s l^2/10, r = sqrt(3/5) l. In six, with y = (l/r)^(5/2) at M = 2l, K = (60 - 70y + 425y^2/4)/l^4
    is least, 824/17l^4 = 48.5/l^4, at y = 28/85, r = (85/28)^(2/5) l = 1.56 l, and R = (-30 + 35y/2)/l^2
    changes sign at y = 12/7, r = (7/12)^(2/5) l = 0.806 l. Far out both tend to anti-de Sitter space's,
    2D(D - 1)/l^4 and -D(D - 1)/l^2, and toward r = 0 both diverge. Each statement is checked here in
    sixty digits, and K is positive at every r drawn: its quadratic in 1/r^2, or in y, has no real root."""
    views = []
    x_axis = cv.axis("$r/\\ell$", 0, 3)
    K_scale = {"kind": "log", "from": 10, "to": 1e6}
    R_scale = {"kind": "symlog", "from": -100, "to": 1e4, "linear": 1}
    clip = "$K$ logarithmic, $10$ to $10^6$; $R$ signed logarithmic, $-10^2$ to $10^4$, linear within $1$."
    for d, system, params, settings, least, least_at, zero_at, empty_K, empty_R, least_tex, zero_tex in (
            ("five", STATIC5, {"M": 3.0, "\\ell": 1.0}, SETTINGS5, 28, math.sqrt(2), math.sqrt(3 / 5), 40, -20,
             "$28/\\ell^4$, at $r = \\sqrt{2}\\ell$", "$r = \\sqrt{3/5}\\ell$"),
            ("six", STATIC6, {"M": 2.0, "\\ell": 1.0}, SETTINGS6, 824 / 17, (85 / 28) ** 0.4, (7 / 12) ** 0.4, 60,
             -30, "$48.5/\\ell^4$, at $r = 1.56\\,\\ell$", "$r = 0.806\\,\\ell$")):
        name = f"Dimensionally continued, {d}"
        names = list(params)
        Kp = cv.nr_published(src, METRIC, system, "kretschmann", ["r"], names)
        Rp = cv.nr_published(src, METRIC, system, "ricci_scalar", ["r"], names)
        grr = cv.nr_published(src, METRIC, system, "inverse_metric_components", ["r"], names, ("r", "r"))
        vals = list(params.values())
        panels = [{"id": "K", "name": "$K\\,\\ell^4$", "aspect": 2, "scale": K_scale, "ticks": cv.value_ticks(K_scale)},
                  {"id": "R", "name": "$R\\,\\ell^2$", "aspect": 2, "scale": R_scale, "ticks": cv.value_ticks(R_scale)}]
        cv.shape_checks(ck, f"{METRIC} {d}", [p["aspect"] for p in panels])
        height = {p["id"]: cv.PLOT_W / p["aspect"] for p in panels}
        curves = (cv.line_curves(ck, f"{METRIC} {d}", "K", "K", Kp, "r", x_axis, K_scale, height["K"], params)
                  + cv.line_curves(ck, f"{METRIC} {d}", "R", "R", Rp, "r", x_axis, R_scale, height["R"], params))

        rs = np.linspace(0.05, 3, 600)
        if d == "five":
            s = 2.0
            K_hand, R_hand = 40 - 24 * s / rs ** 2 + 12 * s * s / rs ** 4, -20 + 6 * s / rs ** 2
        else:
            y = rs ** -2.5
            K_hand, R_hand = 60 - 70 * y + 106.25 * y * y, -30 + 17.5 * y
        ck(f"{name}: K as the module docstring writes it", bool(np.allclose(Kp(rs, *vals), K_hand, rtol=1e-12)))
        ck(f"{name}: R as the module docstring writes it", bool(np.allclose(Rp(rs, *vals), R_hand, rtol=1e-12)))
        ck(f"{name}: K is positive at every r drawn", bool(np.all(Kp(rs, *vals) > 0)))
        ck(f"{name}: the horizon, where g^rr = 0, is r_+ = l", abs(float(grr.exact([1.0] + vals))) < 1e-50)
        ck(f"{name}: K is least, {least:.4g}/l^4, at r = {least_at:.4g} l",
           abs(float(Kp.exact([least_at] + vals)) - least) < 1e-12
           and float(Kp.exact([least_at * 0.999] + vals)) > least and float(Kp.exact([least_at * 1.001] + vals)) > least)
        ck(f"{name}: the figures the caption prints", (round(least, 1), round(least_at, 2), round(zero_at, 3))
           == ((48.5, 1.56, 0.806) if d == "six" else (28, round(math.sqrt(2), 2), round(math.sqrt(0.6), 3))))
        ck(f"{name}: R changes sign at r = {zero_at:.4g} l", abs(float(Rp.exact([zero_at] + vals))) < 1e-12
           and float(Rp.exact([zero_at * 0.99] + vals)) > 0 > float(Rp.exact([zero_at * 1.01] + vals)))
        ck(f"{name}: far out K and R are anti-de Sitter space's", abs(float(Kp.exact([1e8] + vals)) - empty_K) < 1e-9
           and abs(float(Rp.exact([1e8] + vals)) - empty_R) < 1e-9)
        ck(f"{name}: K and R diverge toward r = 0", float(Kp.exact([1e-6] + vals)) > 1e20
           and float(Rp.exact([1e-6] + vals)) > 1e10)
        legend = [["line", "K", "$K$"], ["line", "R", "$R$"], ["line", "horizon", "the horizon"]]
        views.append({
            "id": d, "label": f"{d.capitalize()} dimensions", "kind": "line", "x": x_axis,
            "panels": panels, "curves": curves, "marks": [{"class": "horizon", "x": 1.0}],
            "legend": legend,
            "caption": [f"The Kretschmann scalar $K$ and the Ricci scalar $R$ of the black hole of {d} dimensions "
                        f"against $r$. $K$ diverges at $r = 0$ and is least, {least_tex}. "
                        f"$R$ changes sign at {zero_tex} and tends to ${empty_R}/\\ell^2$."],
            "ricci_sentences": [[0, 2]],
            "per_scalar": cv.scalar_clips("A logarithmic scale from $10$ to $10^6$.",
                                          "A signed logarithmic scale from $-10^2$ to $10^4$, linear within $1$."),
            "settings": settings,
            "clip": clip,
        })
    return views
