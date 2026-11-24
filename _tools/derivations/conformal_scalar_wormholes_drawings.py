"""The drawings of Barcelo and Visser's wormholes of a conformally coupled scalar field, in units of mu,
the half of their length eta, mu = 1.

Every drawing of the wormholes is of one, b = tan(Delta/2) = -1/4, drawn in the isotropic chart: its
throat is the sphere rho = mu sqrt(-b) = mu/2, of areal radius mu (1 + sqrt(-b))^2 = 9 mu/4, and the
sphere rho = mu, where the chart in r ends at r = 4 mu, lies on the side of the first flat region. On
the plane of t and rho the metric is (1 - b mu/rho)^2 (1 + mu/rho)^{-2} (-c^2 dt^2 + (1 + mu/rho)^4 drho^2),
so the rays keep ct -+ rho_* with rho_* = rho + 2 mu ln(rho/mu) - mu^2/rho, which runs over the whole
line, -infinity at rho -> 0 and +infinity at rho -> infinity, for every b: the conformal diagram is the
full diamond, p, q = arctan((ct -+ (rho_* - rho_*t))/mu), with the throat on its axis, rho_*t the
throat's rho_*. The slice t = 0 of the space is symmetric about the throat for every b < 0, since
rho -> mu^2 (-b)/rho keeps the areal radius rho + mu (1 - b) - b mu^2/rho and the proper distance
along it, and the proper distance from the throat is
l = rho - rho_t + mu (1 - b) ln(rho/rho_t) + b mu^2 (1/rho - 1/rho_t), against which the curvature
figure plots the stored Kretschmann scalar.

The chart in r is drawn at another member of the class, w = 1/4 and b = -1/4, whose angle chi has
tan(chi/2) = sqrt 3/2, chi = 81.8 degrees, gamma = 1/7 and sigma = 2/7: its spheres shrink to a point at
r = 4 mu, a naked curvature singularity, as the areal radius r f^{(1 - gamma)/2} Omega goes as
f^{(1 - gamma)/2 - sigma} = f^{1/7} there.
"""
import math

import numpy as np

METRIC = "conformal_scalar_wormholes"
MU, B = 1.0, -0.25
THROAT = MU * math.sqrt(-B)
THROAT_AREAL = MU * (1 + math.sqrt(-B)) ** 2
W_NAKED = 0.25
ISO = {"mu": 1, "b": "-1/4"}
THREE = {"mu": 1, "b": "-1/4", "w": "1/4"}
SETTINGS = "$\\mu = 1$, the unit of every length, and $b = -1/4$."
NAKED = "$\\mu = 1$, $b = -1/4$, and $w = 1/4$, where $\\chi = 81.8°$, $\\gamma = 1/7$, and $\\sigma = 2/7$."
# The reach of the embedding, two spheres each the other's image under rho -> mu^2 (-b)/rho.
NEAR, FAR = 0.05, 5.0


def areal(rho):
    """The areal radius of the drawn wormhole, (rho - b mu)(rho + mu)/rho."""
    rho = np.asarray(rho, dtype=float)
    return (rho - B * MU) * (rho + MU) / rho


def rstar(rho):
    """The tortoise coordinate of the plane of t and rho, rho + 2 mu ln(rho/mu) - mu^2/rho, less its value
    at the throat."""
    def bare(x):
        x = np.asarray(x, dtype=float)
        return x + 2 * MU * np.log(x / MU) - MU ** 2 / x
    return bare(rho) - bare(THROAT)


def proper(rho):
    """The proper distance from the throat along a radius of the slice t = 0, negative on the side of the
    second flat region."""
    rho = np.asarray(rho, dtype=float)
    return (rho - THROAT) + MU * (1 - B) * np.log(rho / THROAT) + B * MU ** 2 * (1 / rho - 1 / THROAT)


# ---------------------------------------------------------------- spacetime diagrams

def register_null(nr):
    nr.DIAGRAMS.append(nr.Diagram(
        METRIC, "isotropic", "radial", "$t$ and $\\rho$", ("t", "\\rho"), (0, 3, -1.5, 1.5), "$\\rho/\\mu$",
        "$ct/\\mu$", ISO, nr.EQUATOR, families=nr.SIDEWAYS, areal=True, areal_contours=(2.5, 3.0, 4.0),
        lines=(("surface", "r", "1", "the sphere $\\rho = \\mu$, where the chart in $r$ ends"),),
        # rho -> 0 is the second flat region, where the conformal factor (1 - b mu/rho)^2 diverges and g^rho rho
        # vanishes with it, and no horizon.
        where="rho", where_is_infinity=True))
    nr.DIAGRAMS.append(nr.Diagram(
        METRIC, "three_parameter", "radial", "$t$ and $r$", ("t", "r"), (4, 7, -1.5, 1.5), "$r/\\mu$", "$ct/\\mu$",
        THREE, nr.EQUATOR, areal=True))
    nr.CAPTIONS[(METRIC, "isotropic", "radial")] = [
        "The plane of $t$ and $\\rho$ ($\\theta = \\pi/2$, $\\phi = 0$) of the wormhole with $b = -1/4$, through "
        "the throat $\\rho = \\mu/2$, where the spheres are smallest, to the second flat region $\\rho \\to 0$, an "
        "infinite distance away, where the cones narrow."]
    nr.CAPTIONS[(METRIC, "three_parameter", "radial")] = [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) of a member of the class with a naked singularity at "
        "$r = 4\\mu$, where the spheres shrink to a point and the cones close."]
    nr.CLOSED_FORMS[(METRIC, "isotropic", "radial")] = (
        lambda t, r: t + rstar(r), lambda t, r: t - rstar(r), None)


HIDDEN = {(METRIC, "three_parameter", "radial"):
          "another member of the class than the wormhole embedded, with a naked singularity at r = 4 mu"}


def slice_table(slices):
    """The moment t = 0 of the embedding view, along rho on the isotropic chart's plane."""
    def isotropic():
        return slices.one(METRIC, lambda m: slices.along(0.0, *m.reach("isotropic", "\\rho")), view_id="wormhole")
    return {(METRIC, "isotropic", "radial"): isotropic}


# ---------------------------------------------------------------- embedding diagram

def register_embedding(em):
    em.DRAWN[METRIC] = lambda ck, src: embedding_views(ck, src, em)
    em.CAPTIONS[(METRIC, "wormhole")] = [
        "The equatorial plane ($\\theta = \\pi/2$) of the wormhole at $t = 0$, from one flat region through "
        "the throat $\\rho = \\mu/2$ to the other, the same on both sides."]


def embedding_views(ck, src, em):
    """The slice t = 0 of the isotropic chart on the equator, one piece through the throat, the surface
    rising with rho: its circles have the areal radius (rho - b mu)(rho + mu)/rho and dz/drho is
    sqrt(g_rho rho - (d areal/d rho)^2), real everywhere, since g_rho rho - (d areal/d rho)^2 =
    (mu (1 - b)/rho - 2 b mu^2/rho^2)(2 + mu (1 - b)/rho) for b < 0. The reach runs between two spheres
    each the other's image under rho -> mu^2 (-b)/rho, which the checks hold to equal circumferences and
    to equal heights above and below the throat."""
    name = "Barcelo-Visser wormhole"
    sl = em.Slice(src, METRIC, "isotropic", "\\rho", "\\phi", {"t": 0, **em.EQUATOR}, ISO)
    size = 2 * float(areal(FAR))
    marks = ([(NEAR, "r", None), (0.1, "r", None), (0.25, "r", None), (THROAT, "throat", "$\\rho = \\mu/2$"),
              (1.0, "r", None), (2.5, "r", None), (FAR, "r", None)])
    whole = em.Piece("whole", "sheet", sl, NEAR, FAR, 0.0, 1,
                     (("edge", "the side $\\rho < \\mu/2$ runs on, flattening, to $\\rho \\to 0$"),
                      ("edge", "the side $\\rho > \\mu/2$ runs on, flattening, to $\\rho \\to \\infty$")),
                     marks, size)
    ck.isometry(f"{name}, through the throat", whole)
    ck.radius(f"{name}, the areal radius (rho - b mu)(rho + mu)/rho", whole, areal, size)
    zt = whole.at(THROAT)[1]
    for rho in (NEAR, 0.1, 0.25):
        dual = THROAT ** 2 / rho
        (r1, z1), (r2, z2) = whole.at(rho), whole.at(dual)
        ck.add(f"{name}: rho = {rho} and its image {dual} have one circumference", abs(r1 - r2), 1e-6 * size)
        ck.add(f"{name}: rho = {rho} and its image {dual} stand as far below and above the throat",
               abs((zt - z1) - (z2 - zt)), 2e-3 * size)
    surface = em.Surface([whole])
    fig = em.figure_of([surface], {"sheet": "cover"}, size)
    em.ring_label(fig, [0, 0, 0], *whole.at(THROAT), "$\\rho = \\mu/2$", dx=14)
    em.ring_label(fig, [0, 0, 0], *whole.at(1.0), "$\\mu$")
    em.ring_label(fig, [0, 0, 0], *whole.at(0.25), "$\\mu/4$")
    fig.legend("fill", "cover", "the slice $t = 0$, which $t$ and $\\rho$ cover")
    fig.legend("line", "r", "$\\rho$ constant")
    fig.legend("line", "throat", "the throat, of circumference $9\\pi\\mu/2$")
    fig.legend("line", "meridian", "$\\phi$ constant, every $15°$")
    return [em.view("wormhole", "The wormhole", "$\\mu$", [surface], fig.done(), settings=SETTINGS)]


# ---------------------------------------------------------------- conformal diagram

def register_conformal(cf):
    cf.DRAWN[METRIC] = lambda ck, src: conformal_views(ck, src, cf)
    cf.CAPTIONS[(METRIC, "isotropic")] = [
        "The wormhole in $t$ and $\\rho$, each point a sphere, with its two flat regions, $\\rho \\to \\infty$ on the "
        "right and $\\rho \\to 0$ on the left, joined at the throat, and no horizon."]


def conformal_views(ck, src, cf):
    """The isotropic chart's plane of t and rho in the full diamond, p, q = arctan((ct -+ (rho_* - rho_*t))/mu),
    checked against the published metric at random points of both sides; the curvature is finite on the
    throat and at the sphere rho = mu, and rho_* runs to -+infinity at both ends, which are the two
    spatial infinities."""
    name = "Barcelo-Visser wormhole"
    pl = cf.Plane(src, METRIC, "isotropic", ("t", "\\rho"), cf.EQUATOR, ISO)
    rng = np.random.default_rng(1999)

    def fmap(t, rho):
        return cf.mink_pq(t, rstar(rho))
    ck.chart(name, pl, fmap, rng.uniform(-20, 20, 4000), 10.0 ** rng.uniform(-3, 3, 4000), lambda t, r: (1, 0))
    ck.finite(f"{name}: the curvature is finite on the throat and on the sphere rho = mu",
              pl.kretschmann(rng.uniform(-5, 5, 50), np.concatenate([np.full(25, THROAT), np.full(25, 1.0)])))
    ck.limit(f"{name}: rho_* runs to -infinity at rho -> 0 and to infinity at rho -> infinity",
             [float(np.tanh(rstar(1e-9) / 1e3)), float(np.tanh(rstar(1e9) / 1e3))], [-1.0, 1.0], 1e-6)
    v = cf.View("isotropic", "Isotropic", [-cf.PI - 0.35, cf.PI + 0.35, -cf.PI - 0.25, cf.PI + 0.25], "isotropic")
    v.fill("region", cf.DIAMOND)
    v.fill("cover", cf.DIAMOND)
    cf.grid(v, "r", lambda r, t: fmap(t, r), (0.05, 0.1, 0.25, 1.0, 2.5, 5.0), cf.S_ALL)
    cf.grid(v, "t", fmap, (-4, -2, -1, 0, 1, 2, 4), cf.S_ALL)
    v.line("throat", [[[0, -cf.PI], [0, cf.PI]]])
    cf.diamond_edges(v)
    cf.label_on(v, fmap(0, 1.0), "$\\rho = \\mu$")
    cf.label_on(v, fmap(0, 0.25), "$\\mu/4$")
    v.label_xt([0, 0.3], "throat", "l", "small", dx=6)
    v.legend("cover", "the whole spacetime, which $t$ and $\\rho$ cover")
    v.legend("r", "$\\rho$ constant")
    v.legend("t", "$ct$ constant")
    v.legend("throat", "the throat $\\rho = \\mu/2$")
    v.set(settings=SETTINGS)
    wormhole = cf.slices.moments(METRIC)[0]
    rho = np.geomspace(*wormhole.reach("isotropic", "\\rho"), 200)
    v.slice(wormhole, [fmap(0 * rho, rho)])
    return [v]


# ---------------------------------------------------------------- curvature figure

def register_curvature(cv):
    cv.DRAWN[METRIC] = lambda ck, src: curvature_views(ck, src, cv)


def curvature_views(ck, src, cv):
    """The stored Kretschmann scalars, read from the metric file and never recomputed: the isotropic chart's
    of the drawn wormhole against the proper distance l from the throat, through the closed horizontal
    expression of l, finite everywhere, greatest near the throat and falling toward both flat regions; and
    the chart in r at the drawn member with a naked singularity against r, growing without bound toward
    r = 4 mu. Both stored Ricci scalars vanish."""
    import sympy as sp
    from scipy.optimize import brentq
    figures = []

    lo, hi = -7.5, 7.5
    x_axis = cv.axis("$l/\\mu$", lo, hi)
    scale = {"kind": "log", "from": 1e-4, "to": 1}
    panels = [{"id": "K", "name": "$K\\,\\mu^4$", "aspect": 2, "scale": scale, "ticks": cv.value_ticks(scale)}]
    cv.shape_checks(ck, METRIC, [p["aspect"] for p in panels])
    params = {"mu": MU, "b": B}
    K = cv.nr_published(src, METRIC, "isotropic", "kretschmann", ["rho"], ["mu", "b"])
    Ric = cv.nr_published(src, METRIC, "isotropic", "ricci_scalar", ["rho"], ["mu", "b"])
    for at in (NEAR, THROAT, 1.0, FAR):
        ck(f"{METRIC} isotropic R = 0 at rho = {at}", abs(float(Ric.exact([at, MU, B]))) < 1e-30)
    # The ends a hair inside the axis, so that rounding leaves no point outside it.
    ends = [brentq(lambda x: float(proper(x)) - edge, 1e-6, 1e3, xtol=1e-15) for edge in (lo + 1e-9, hi - 1e-9)]
    rhos = np.geomspace(ends[0], ends[1], cv.SAMPLES)
    values = K(rhos, MU, B)
    ls = proper(rhos)
    ck(f"{METRIC} isotropic K lies within its scale", bool(np.all(values > scale["from"]) and np.all(values < scale["to"])))
    top = int(np.argmax(values))
    ck(f"{METRIC} isotropic K is greatest within mu/10 of the throat", abs(float(ls[top])) < 0.1)
    ck(f"{METRIC} isotropic K falls from its greatest value toward both flat regions",
       bool(np.all(np.diff(values[:top + 1]) > 0) and np.all(np.diff(values[top:]) < 0)))
    p0, p1 = cv.span(scale)
    X = cv.PLOT_W * (ls - lo) / (hi - lo)
    Y = (cv.PLOT_W / 2) * (1 - (cv.position(scale, values) - p0) / (p1 - p0))
    keep = cv.thinned(X, Y, cv.TOLERANCE)
    entry = cv.curve_entry("K", "K", K, "rho", params, [[float(ls[k]) + 0.0, float(values[k]) + 0.0] for k in keep])
    entry["horizontal"] = {"expression": "rho - 1/2 + (5/4)*log(2*rho) - (1/4)*(1/rho - 2)",
                           "values": [float(rhos[k]) for k in keep]}
    rho = sp.Symbol("rho", positive=True)
    l_text = sp.lambdify(rho, sp.sympify(entry["horizontal"]["expression"]), "numpy")
    ck(f"{METRIC} isotropic the horizontal expression is the proper distance from the throat",
       float(np.max(np.abs(l_text(rhos) - ls))) < 1e-12)
    cv.agrees(ck, f"{METRIC} isotropic K", K, [[float(rhos[k]), MU, B] for k in keep], [float(values[k]) for k in keep])
    figures.append({
        "id": "wormhole", "label": "The wormhole", "kind": "line", "system": "isotropic", "x": x_axis,
        "panels": panels, "curves": [entry], "marks": [{"class": "front", "x": 0.0}],
        "legend": [["line", "K", "$K$"], ["line", "front", "the throat $l = 0$"]],
        "caption": ["The Kretschmann scalar $K$ of the wormhole with $b = -1/4$ against the proper distance $l$ from "
                    "the throat, with $l < 0$ toward the second flat region, finite everywhere and greatest near "
                    "the throat. The Ricci scalar is zero."],
        "settings": SETTINGS, "clip": "$K\\,\\mu^4$ logarithmic, none clipped."})

    lo, hi = 4.1, 8.0
    x_axis = cv.axis("$r/\\mu$", lo, hi)
    scale = {"kind": "log", "from": 1e-5, "to": 1e2}
    panels = [{"id": "K", "name": "$K\\,\\mu^4$", "aspect": 2, "scale": scale, "ticks": cv.value_ticks(scale)}]
    cv.shape_checks(ck, METRIC, [p["aspect"] for p in panels])
    names = ["mu", "b", "w"]
    K = cv.nr_published(src, METRIC, "three_parameter", "kretschmann", ["r"], names)
    Ric = cv.nr_published(src, METRIC, "three_parameter", "ricci_scalar", ["r"], names)
    for at in (4.1, 5.0, 8.0):
        ck(f"{METRIC} three_parameter R = 0 at r = {at}", abs(float(Ric.exact([at, MU, B, W_NAKED]))) < 1e-30)
    near, nearer = (float(K.exact([4 + d, MU, B, W_NAKED])) for d in (1e-3, 1e-4))
    ck(f"{METRIC} three_parameter K grows without bound toward r = 4 mu", near > 1e6 and nearer > 50 * near)
    rs = np.linspace(lo, hi, cv.SAMPLES)
    values = K(rs, MU, B, W_NAKED)
    ck(f"{METRIC} three_parameter K lies within its scale",
       bool(np.all(values > scale["from"]) and np.all(values < scale["to"])))
    ck(f"{METRIC} three_parameter K falls outward", bool(np.all(np.diff(values) < 0)))
    curves = cv.line_curves(ck, METRIC, "K", "K", K, "r", x_axis, scale, cv.PLOT_W / 2,
                            {"mu": MU, "b": B, "w": W_NAKED})
    figures.append({
        "id": "naked", "label": "A naked singularity", "kind": "line", "system": "three_parameter", "x": x_axis,
        "panels": panels, "curves": curves, "marks": [],
        "legend": [["line", "K", "$K$"]],
        "caption": ["The Kretschmann scalar $K$ of a member of the class with a naked singularity against $r$, growing without "
                    "bound toward the singularity $r = 4\\mu$. The Ricci scalar is zero."],
        "settings": NAKED, "clip": "$K\\,\\mu^4$ logarithmic, from $r = 4.1\\mu$."})
    return figures
