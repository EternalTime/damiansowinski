"""The drawings of the repulson, every one at its massless member, M = 0, in units of the charge q.

Then Delta = 1 - 4q^2/r^2 and the singularity is r_0 = 2q. In both charts the radial light rays obey
c dt = +-sqrt(Delta) dr, since the string metric is the Einstein metric times 1/sqrt(Delta), so
ct -+ r_* is constant along them with the tortoise coordinate
r_* = sqrt(r^2 - 4q^2) - 2q arccos(2q/r), which vanishes at r = 2q and grows as r far out. A ray
moving in reaches the singularity after the finite time r_*/c.

The conformal diagram is Minkowski's triangle by p, q' = arctan(ct -+ r_*), with the singularity
r = 2q where Minkowski's centre stands. The embedding diagram draws each chart's moment t = 0 on its
equator: the Einstein metric's in three dimensional Minkowski space, since every circle there grows
faster than the distance out to it, (d rho/dr)^2 - g_rr = 4(1 - 3q^2/r^2)/(r^2 Delta^{3/2}) > 0 for
r > 2q, and the string metric's as the flat plane outside the circle r = 2q, which is what
g_rr = 1 and g_phiphi = r^2 make it. The curvature figure is the two stored scalars of the Einstein
metric against r. repulson.md beside this file records the charts.
"""
import numpy as np

METRIC = "repulson"
SYSTEMS = ("isotropic", "string")
PARAMS = {"M": 0, "q": 1}
EQUATOR = {"theta": "pi/2", "phi": "0"}
R0, TOP = 2.0, 6.0
# The Einstein metric's moment is drawn from r = 2.05q: nearer the singularity the surface turns toward a
# light cone of Minkowski space, its metric length falls to about (r - 2q)/q of its drawn length, and its
# chords, crowded toward the inner end by KNOTS and written to eleven decimals of the piece's extent, so
# that rounding does not swamp the small interval of a nearly null chord, keep the proper distance within ALONG.
EMB_LO = 2.05
KNOTS = 200
VIEWS = {"isotropic": "einstein", "string": "string"}
SETTINGS = "$M = 0$, $q = 1$."


def rstar(r):
    """The tortoise coordinate at M = 0, q = 1: sqrt(r^2 - 4) - 2 arccos(2/r), zero at r = 2."""
    r = np.asarray(r, dtype=float)
    return np.sqrt(np.maximum(r * r - 4.0, 0.0)) - 2.0 * np.arccos(np.minimum(2.0 / r, 1.0))


# ---------------------------------------------------------------- spacetime diagrams

def register_null(nr):
    for system in SYSTEMS:
        nr.DIAGRAMS.append(nr.Diagram(
            METRIC, system, "radial", "Radial plane", ("t", "r"), (R0, TOP, -2, 2), "$r/q$", "$ct/q$",
            PARAMS, EQUATOR, cones=(5, 5)))
        nr.CLOSED_FORMS[(METRIC, system, "radial")] = (
            lambda t, r: t + rstar(r), lambda t, r: t - rstar(r), lambda t, r: r > R0 + 5e-4)
    nr.CAPTIONS[(METRIC, "isotropic", "radial")] = [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) of the massless repulson, its radial light "
        "rays reaching the singularity $r = 2q$ in a finite time."]
    nr.CAPTIONS[(METRIC, "string", "radial")] = [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) of the massless repulson's string metric, "
        "whose radial light rays are those of the Einstein metric."]


def slice_table(slices):
    """The moment t = 0 of each chart's embedding view, on that chart's spacetime diagram."""
    return {(METRIC, system, "radial"):
            (lambda system=system: slices.one(METRIC, lambda m: slices.along(0.0, *m.reach(system, "r")),
                                              view_id=VIEWS[system]))
            for system in SYSTEMS}


# ---------------------------------------------------------------- embedding diagram

def register_embedding(em):
    em.DRAWN[METRIC] = lambda ck, src: embedding_views(ck, src, em)
    em.CAPTIONS[(METRIC, "einstein")] = [
        "The equator of the massless repulson at $t = 0$ in its Einstein metric, drawn in Minkowski space from "
        "$r = 2.05q$, beside the singularity $r = 2q$, out to $6q$."]
    em.CAPTIONS[(METRIC, "string")] = [
        "The equator of the massless repulson at $t = 0$ in its string metric, the flat plane outside the "
        "circle $r = 2q$."]


def embedding_views(ck, src, em):
    size = 2 * TOP
    fixed = {"t": 0, "theta": "pi/2"}
    marks = (2.5, 3.0, 4.0, 5.0)

    # The Einstein metric: rho = r Delta^(1/4), and the profile climbs in Minkowski space.
    sl = em.Slice(src, METRIC, "isotropic", "r", "\\phi", fixed, PARAMS, space="minkowski")
    rr = np.linspace(R0 + 1e-3, 40.0, 400)
    defect = sl.defect_at(rr)
    expected = -4 * (1 - 3 / rr ** 2) / (rr ** 2 * (1 - 4 / rr ** 2) ** 1.5)
    ck.add("repulson: (drho/dr)^2 - g_rr is 4(1 - 3q^2/r^2)/(r^2 Delta^(3/2)), so every circle grows faster "
           "than the distance out to it", float(np.max(np.abs(defect - expected) / np.abs(expected))), 1e-9)
    ck.stops("repulson Einstein metric in flat space", em.Slice(src, METRIC, "isotropic", "r", "\\phi", fixed, PARAMS),
             rr)
    einstein = em.Piece("einstein", "sheet", sl, EMB_LO, TOP, 0.0, 1,
                        (("stops", "the surface turns toward a light cone of Minkowski space as it nears the "
                                   "singularity $r = 2q$"),
                         ("edge", "the surface runs on to $r \\to \\infty$")),
                        [(EMB_LO, "r", "$r = 2.05q$")] + [(r, "r", None) for r in marks] + [(TOP, "r", "$6q$")], size,
                        knots=EMB_LO + (TOP - EMB_LO) * np.geomspace(2e-4, 1, KNOTS), digits=1e-11)
    ck.isometry("repulson Einstein metric in Minkowski space", einstein)
    ck.radius("repulson Einstein metric, rho = r Delta^(1/4)", einstein,
              lambda r: np.sqrt(r) * (r * r - 4.0) ** 0.25, size)
    surface = em.Surface([einstein])
    fig = em.figure_of([surface], {"sheet": "cover"}, size)
    em.ring_label(fig, [0, 0, 0], *einstein.at(TOP), "$6q$")
    fig.legend("fill", "cover", "the equator at $t = 0$")
    fig.legend("line", "r", "$r$ constant, from $2.05q$ to $6q$")
    fig.legend("line", "meridian", "$\\phi$ constant, every $15°$")
    views = [em.view("einstein", "Einstein metric", "$q$", [surface], fig.done(),
                     settings=SETTINGS + " Every length along the surface is measured with $dX^2 + dY^2 - dZ^2$.")]

    # The string metric: g_rr = 1 and g_phiphi = r^2, the plane outside r = 2q.
    st = em.Slice(src, METRIC, "string", "r", "\\phi", fixed, PARAMS)
    ck.plane("repulson string metric: the equator is flat", st, np.linspace(R0 + 1e-3, TOP, 80))
    plane = em.Piece("string", "sheet", st, R0, TOP, 0.0, 1,
                     (("edge", "the singularity $r = 2q$, a circle of circumference $4\\pi q$"),
                      ("edge", "the plane runs on to $r \\to \\infty$")),
                     [(R0, "r", "$r = 2q$")] + [(r, "r", None) for r in marks] + [(TOP, "r", "$6q$")], size)
    ck.isometry("repulson string metric", plane)
    ck.radius("repulson string metric, rho = r", plane, lambda r: r, size)
    ck.form("repulson string metric, the height is zero", plane, lambda r: 0 * r, size)
    flat = em.Surface([plane])
    fig = em.figure_of([flat], {"sheet": "cover"}, size)
    fig.legend("fill", "cover", "the equator at $t = 0$")
    fig.legend("line", "r", "$r$ constant, from $2q$ to $6q$")
    fig.legend("line", "meridian", "$\\phi$ constant, every $15°$")
    views.append(em.view("string", "String metric", "$q$", [flat], fig.done(), settings=SETTINGS))
    return views


# ---------------------------------------------------------------- conformal diagram

def register_conformal(cf):
    cf.DRAWN[METRIC] = lambda ck, src: conformal_views(ck, src, cf)
    cf.CAPTIONS.update({(METRIC, system): [
        "The massless repulson" + ("" if system == "isotropic" else " in its string metric")
        + " on its plane of $t$ and $r$, each point in the diagram a 2-sphere, Minkowski's triangle with the "
        "singularity $r = 2q$ where the centre stands."] for system in SYSTEMS})


def conformal_views(ck, src, cf):
    PI, HALF = cf.PI, cf.HALF
    box = [-0.35, PI + 0.35, -PI - 0.25, PI + 0.25]

    def fmap(t, r):
        return cf.mink_pq(t, rstar(r))
    labels = {"isotropic": "Einstein Metric, Isotropic", "string": "String Metric"}
    views = []
    for system in SYSTEMS:
        plane = cf.Plane(src, METRIC, system, ("t", "r"), EQUATOR, PARAMS)
        ck.chart(f"repulson {system}", plane, fmap, ck.uniform(-20, 20), ck.uniform(R0 + 0.01, 30),
                 lambda t, r: (1, 0))
        ck.diverges(f"repulson {system}: the Kretschmann scalar diverges at r = 2q",
                    plane.kretschmann(0, R0 + 1e-5), plane.kretschmann(0, R0 + 1e-6))
        ck.limit(f"repulson {system}: r_* vanishes at the singularity", [float(rstar(R0))], [0.0], 1e-12)
        v = cf.View(system, labels[system], box, system)
        v.fill("region", cf.TRIANGLE)
        v.fill("cover", cf.TRIANGLE)
        cf.grid(v, "r", lambda r, t: fmap(t, r), (2.5, 3, 4, 6), cf.S_ALL)
        cf.grid(v, "t", fmap, (-4, -2, -1, 0, 1, 2, 4), R0 + cf.S_POS)
        v.line("singular", [[[0, -PI], [0, PI]]], zig=True)
        v.line("scri", [[[0, PI], [PI, 0]], [[PI, 0], [0, -PI]]])
        for at, text, anchor, dx, dy in (((PI, 0), "$i^0$", "l", 6, 0), ((0, PI), "$i^+$", "b", 0, -6),
                                         ((0, -PI), "$i^-$", "t", 0, 6)):
            v.layers.append({"kind": "point", "class": "infinity", "at": [round(at[0], 4), round(at[1], 4)]})
            v.label_xt(at, text, anchor, dx=dx, dy=dy)
        v.label_xt([HALF, HALF], "$\\mathscr{I}^+$", "bl", dx=5, dy=-3)
        v.label_xt([HALF, -HALF], "$\\mathscr{I}^-$", "tl", dx=5, dy=3)
        v.label_xt([0, 0.25], "$r = 2q$", "l", dx=8)
        v.legend("cover", "the whole spacetime")
        v.legend("r", "$r$ constant")
        v.legend("t", "$ct$ constant")
        v.legend("singular", "the singularity $r = 2q$")
        v.legend("scri", "null infinity $\\mathscr{I}^\\pm$")
        m = cf.slices.moments(METRIC, VIEWS[system])[0]
        lo, hi = m.reach(system, "r")
        rr = np.linspace(lo, hi, 60)
        v.slice(m, [fmap(0 * rr, rr)])
        v.set(settings=SETTINGS + " Light rays keep $ct \\pm r_*$, with $r_* = \\sqrt{r^2 - 4q^2} - 2q\\arccos(2q/r)$.")
        views.append(v)
    return views


# ---------------------------------------------------------------- curvature figure

def register_curvature(cv):
    cv.DRAWN[METRIC] = lambda ck, src: curvature_views(ck, src, cv)


def curvature_views(ck, src, cv):
    """The stored Kretschmann and Ricci scalars of the Einstein metric at M = 0 against r, from the
    singularity r = 2q out to 6q, read from the metric file and never recomputed. Checks hold the
    Ricci scalar to Kallosh's 4q^2(2q^2 + r^2)/(r(r^2 - 4q^2)^{5/2}) and both scalars to growing without
    bound toward r = 2q and falling outward."""
    import sympy as sp

    xaxis = cv.axis("$r/q$", R0, TOP)
    scales = {"K": {"kind": "log", "from": 1e-4, "to": 1e4}, "R": {"kind": "log", "from": 1e-3, "to": 1e3}}
    panels = [{"id": field, "name": "$" + field + ("\\,q^4$" if field == "K" else "\\,q^2$"),
               "aspect": 2, "scale": scales[field], "ticks": cv.value_ticks(scales[field])} for field in ("K", "R")]
    cv.shape_checks(ck, METRIC, [p["aspect"] for p in panels])
    K = cv.nr_published(src, METRIC, "isotropic", "kretschmann", ["r"], list(PARAMS))
    R = cv.nr_published(src, METRIC, "isotropic", "ricci_scalar", ["r"], list(PARAMS))
    r = R.coords[0]
    kallosh = 4 * (2 + r ** 2) / (r * (r ** 2 - 4) ** sp.Rational(5, 2))
    for at in (2.5, 3.0, 4.5):
        ck(f"{METRIC} R at r = {at}q is Kallosh's", abs(float(R.exact([at, 0, 1])) - float(kallosh.subs(r, at))) < 1e-12)
    xs = np.linspace(R0 + 0.02, TOP, 400)
    curves = []
    for name, fn in (("K", K), ("R", R)):
        values = fn(xs, 0, 1)
        ck(f"{METRIC} {name} is positive and falls outward", bool(np.all(values > 0) and np.all(np.diff(values) < 0)))
        ck(f"{METRIC} {name} grows without bound toward r = 2q",
           float(fn(np.array([R0 + 1e-6]), 0, 1)[0]) > 1e3 * float(fn(np.array([R0 + 1e-2]), 0, 1)[0]))
        curves += cv.line_curves(ck, METRIC, name, name, fn, "r", xaxis, scales[name], cv.PLOT_W / 2, PARAMS)
    k_clip = "$K\\,q^4$ logarithmic, cut near $r = 2q$."
    r_clip = "$R\\,q^2$ logarithmic, cut near $r = 2q$."
    return [{"id": "massless", "label": "Massless repulson", "kind": "line", "system": "isotropic", "x": xaxis,
             "panels": panels, "curves": curves, "marks": [],
             "legend": [["line", "K", "$K$"], ["line", "R", "$R$"]],
             "caption": ["The Kretschmann scalar $K$ and the Ricci scalar $R$ of the massless repulson's Einstein metric "
                         "against $r$, out from the singularity $r = 2q$. $K$ grows without bound toward the singularity.",
                         "$R$ is positive and grows without bound toward the singularity."],
             "ricci_sentences": [[1, 0]], "settings": SETTINGS, "clip": k_clip + " " + r_clip,
             "per_scalar": cv.scalar_clips(k_clip, r_clip)}]
