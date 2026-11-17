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
