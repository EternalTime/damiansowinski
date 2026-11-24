"""The drawings of the Chamseddine-Volkov monopole, every one at a = 1, the unit of every length.

The plane of t and rho is a^2 sinh(rho)/R times the flat plane -dt^2 + d rho^2, so the radial light
rays keep t -+ rho, the spacetime diagram is the half plane rho >= 0 with every ray at 45 degrees, and
the conformal diagram is Minkowski's triangle by p, q = arctan(t -+ rho), with the regular centre
rho = 0 on its vertical edge, the authors' Fig. 1.

The equator at t = 0 has g_rho_rho = a^2 sinh(rho)/R and the circles of radius r = a sqrt(R sinh rho),
so (dr/d rho)^2 = N g_rho_rho with N = rho^2/R^2 > 1 off the centre: every circle grows faster than the
distance out to it, no surface of revolution in flat space holds the moment, and it is drawn in three
dimensional Minkowski space, level at the centre, where N = 1 + rho^2/9 + O(rho^4).

The curvature figure is the two stored scalars against rho, both finite at the centre,
K = 20/(27 a^4) and R = -2/a^2, and falling toward zero far away.
chamseddine_volkov_monopole.md beside this file records the chart.
"""
import numpy as np

METRIC = "chamseddine_volkov_monopole"
SYSTEM = "chamseddine_volkov"
PARAMS = {"a": 1}
EQUATOR = {"theta": "pi/2", "varphi": "0"}
VIEW = "equator"
TOP = 3.0           # the embedding's outer edge in rho, where r = 4.72a
SETTINGS = "$a = 1$, the unit of every length."


SERIES = 40         # terms of the series below, enough for rho < 1 to the last digit


def _numerator(rho):
    """R^2 sinh^2 rho = rho sinh(2 rho) - rho^2 - sinh^2 rho and its derivative along rho, from the
    series sum_{n >= 2} 2^(2n-1) (2n-1) rho^(2n)/(2n)!, whose terms are all positive, where rho < 1:
    near the centre the closed form subtracts numbers that agree to all but the last few digits,
    since R^2 = rho^2 - rho^4/9 + O(rho^6) while each of its terms is of order one."""
    rho = np.asarray(rho, dtype=float)
    f, df = np.zeros_like(rho), np.zeros_like(rho)
    term = rho ** 4                            # n = 2: 2^3 . 3 rho^4 / 4!
    for n in range(2, SERIES):
        f = f + term
        df = df + term * 2 * n / np.where(rho == 0, 1.0, rho)
        term = term * 4 * rho ** 2 / ((2 * n - 1) * (2 * n + 2))
    return f, df


def R_of(rho):
    """R at a = 1, from the series where rho < 1 and from the closed form beyond."""
    rho = np.asarray(rho, dtype=float)
    small = rho < 1.0
    safe = np.where(small, 1.0, rho)
    closed = np.sqrt(2 * safe / np.tanh(safe) - safe ** 2 / np.sinh(safe) ** 2 - 1)
    f, _ = _numerator(np.where(small, rho, 0.0))
    sinh = np.sinh(np.where(small, rho, 1.0))
    series = np.sqrt(f) / np.where(rho == 0, 1.0, sinh)
    return np.where(small, np.where(rho == 0, 0.0, series), closed)


def dR_of(rho):
    """dR/drho at a = 1, half the derivative of R^2 over R, which tends to 1 at the centre."""
    rho = np.asarray(rho, dtype=float)
    small = rho < 1.0
    R = R_of(rho)
    safe = np.where(small, 1.0, rho)
    closed = (np.cosh(safe) * np.sinh(safe) ** 2 - 2 * safe * np.sinh(safe) + safe ** 2 * np.cosh(safe)) / (
        np.where(small, 1.0, R) * np.sinh(safe) ** 3)
    x = np.where(small, rho, 0.5)
    f, df = _numerator(x)
    s, c = np.sinh(x), np.cosh(x)
    with np.errstate(divide="ignore", invalid="ignore"):
        series = (df * s - 2 * f * c) / (2 * np.sqrt(f) * s ** 2)
    return np.where(small, np.where(rho == 0, 1.0, series), closed)


def areal(rho):
    """The areal radius a sqrt(R sinh rho) at a = 1."""
    rho = np.asarray(rho, dtype=float)
    return np.sqrt(R_of(rho) * np.sinh(rho))


# ---------------------------------------------------------------- spacetime diagram

def register_null(nr):
    # The rays stop at the centre: continued to rho < 0 the metric's sinh(rho) changes sign and the
    # future cone turns round, so a ray traced past rho = 0 would come back along it.
    nr.DIAGRAMS.append(nr.Diagram(
        METRIC, SYSTEM, "radial", "Radial plane", ("t", "\\rho"), (0, 4, -2, 2), "$\\rho$", "$t$",
        PARAMS, EQUATOR, cones=(5, 5), inside=True))
    nr.CLOSED_FORMS[(METRIC, SYSTEM, "radial")] = (
        lambda t, r: t + r, lambda t, r: t - r, lambda t, r: r > 5e-4)
    nr.CAPTIONS[(METRIC, SYSTEM, "radial")] = [
        "The plane of $t$ and $\\rho$ ($\\theta = \\pi/2$, $\\varphi = 0$) of the monopole, $a^2\\sinh\\rho/R$ "
        "times the flat plane, its radial light rays at 45° with $t \\pm \\rho$ constant."]


def slice_table(slices):
    """The moment t = 0 of the embedding view, on the spacetime diagram."""
    return {(METRIC, SYSTEM, "radial"):
            (lambda: slices.one(METRIC, lambda m: slices.along(0.0, *m.reach(SYSTEM, "\\rho")), view_id=VIEW))}


# ---------------------------------------------------------------- embedding diagram

def register_embedding(em):
    em.DRAWN[METRIC] = lambda ck, src: embedding_views(ck, src, em)
    em.CAPTIONS[(METRIC, VIEW)] = [
        "The equator of the monopole at $t = 0$, drawn in Minkowski space from the regular centre out to "
        "$\\rho = 3$ ($r = 4.72a$)."]


def embedding_views(ck, src, em):
    size = 2 * float(areal(TOP))
    fixed = {"t": 0, "theta": "pi/2"}
    marks = (1.0, 2.0)

    # R is handed to the slice as the function R_of computes, the published definition summed without
    # cancellation near the centre, where the closed form loses its digits below rho = 1e-5.
    numeric = {"R": (R_of, dR_of)}
    sl = em.Slice(src, METRIC, SYSTEM, "\\rho", "\\varphi", fixed, PARAMS, numeric=numeric, space="minkowski")
    rr = np.linspace(0.05, 12.0, 400)
    R = R_of(rr)
    expected = np.sinh(rr) / R * (1 - rr ** 2 / R ** 2)
    ck.add("chamseddine_volkov_monopole: g_rhorho - (dr/drho)^2 is a^2 (sinh rho/R)(1 - rho^2/R^2), so every "
           "circle grows faster than the distance out to it",
           float(np.max(np.abs(sl.defect_at(rr) - expected) / np.abs(expected))), 1e-9)
    ck.stops("chamseddine_volkov_monopole in flat space",
             em.Slice(src, METRIC, SYSTEM, "\\rho", "\\varphi", fixed, PARAMS, numeric=numeric), rr)
    ball = em.Piece("equator", "sheet", sl, 0.0, TOP, 0.0, 1,
                    (("axis", "the regular centre $\\rho = 0$, where the surface is level"),
                     ("edge", "the surface runs on to $\\rho \\to \\infty$")),
                    [(r, "r", None) for r in marks] + [(TOP, "r", "$\\rho = 3$")], size)
    ck.isometry("chamseddine_volkov_monopole in Minkowski space", ball)
    ck.radius("chamseddine_volkov_monopole, r = a sqrt(R sinh rho)", ball, areal, size)
    surface = em.Surface([ball])
    fig = em.figure_of([surface], {"sheet": "cover"}, size)
    em.ring_label(fig, [0, 0, 0], *ball.at(TOP), "$\\rho = 3$")
    fig.legend("fill", "cover", "the equator at $t = 0$")
    fig.legend("line", "r", "$\\rho$ constant, at $1$, $2$, and $3$")
    fig.legend("line", "meridian", "$\\varphi$ constant, every $15°$")
    return [em.view(VIEW, "Equator", "$a$", [surface], fig.done(),
                    settings=SETTINGS + " Every length along the surface is measured with $dX^2 + dY^2 - dZ^2$.")]


# ---------------------------------------------------------------- conformal diagram

def published_scalar(nr, field):
    """The stored scalar at a = 1 as a function of rho, read through the checker's Reader with R
    written out and evaluated in 120 digits: near the centre its terms cancel through more than sixty,
    and a plane's own Kretschmann scalar, simplified symbolically first, takes over half an hour."""
    import mpmath
    import sympy as sp

    _, entry, reader = nr.load(METRIC, SYSTEM)
    expr = reader(nr.strip_lhs(entry[field])).subs(reader.held).doit().subs(reader.parameters["a"], 1)
    f = sp.lambdify(reader.symbol["\\rho"], expr, "mpmath")

    def at(rho):
        with mpmath.workdps(120):
            return float(f(mpmath.mpf(rho)))
    return at


def register_conformal(cf):
    cf.DRAWN[METRIC] = lambda ck, src: conformal_views(ck, src, cf)
    cf.CAPTIONS[(METRIC, SYSTEM)] = [
        "The monopole on its plane of $t$ and $\\rho$, each point in the diagram a 2-sphere, Minkowski's triangle "
        "with the regular centre $\\rho = 0$ on its vertical edge."]


def conformal_views(ck, src, cf):
    PI = cf.PI
    box = [-0.35, PI + 0.35, -PI - 0.25, PI + 0.25]
    fmap = cf.mink_pq
    plane = cf.Plane(src, METRIC, SYSTEM, ("t", "\\rho"), EQUATOR, PARAMS)
    ck.chart("Chamseddine-Volkov monopole", plane, fmap, ck.uniform(-20, 20), ck.uniform(0.01, 30),
             lambda t, r: (1, 0))
    centre = published_scalar(cf.nr, "kretschmann")
    near = [centre(x) for x in ("1e-2", "1e-3", "1e-4")]
    ck.finite("Chamseddine-Volkov monopole: the centre is regular", near)
    ck.limit("Chamseddine-Volkov monopole: K tends to 20/27 a^-4 at the centre", near[-1], 20 / 27, 1e-7)
    v = cf.View(SYSTEM, "Chamseddine-Volkov", box, SYSTEM)
    v.fill("region", cf.TRIANGLE)
    v.fill("cover", cf.TRIANGLE)
    cf.grid(v, "r", lambda r, t: fmap(t, r), (0.5, 1, 2, 4), cf.S_ALL)
    cf.grid(v, "t", fmap, (-4, -2, -1, 0, 1, 2, 4), cf.S_POS)
    cf.triangle_edges(v, centre="$\\rho = 0$")
    v.legend("cover", "the whole spacetime")
    v.legend("r", "$\\rho$ constant, at $1/2$, $1$, $2$, and $4$")
    v.legend("t", "$t$ constant")
    v.legend("centre", "$\\rho = 0$, a regular centre")
    m = cf.slices.moments(METRIC, VIEW)[0]
    lo, hi = m.reach(SYSTEM, "\\rho")
    rr = np.linspace(lo, hi, 60)
    v.slice(m, [fmap(0 * rr, rr)])
    v.set(settings=SETTINGS + " Light rays keep $t \\pm \\rho$.")
    return [v]


# ---------------------------------------------------------------- curvature figure

def register_curvature(cv):
    cv.DRAWN[METRIC] = lambda ck, src: curvature_views(ck, src, cv)


def in_eighty_digits(published):
    """The published scalar evaluated as it stands, in eighty digits at every sample. Its value holds R
    unreduced, as the checker compares it, so near the centre and far out its terms, each of order
    R^4 sinh^6(rho), cancel to a small remainder that doubles cannot hold: at rho = 1e-3, where the curves start,
    about thirty digits cancel. The sixty-digit check of curvature.py reads the same expression
    through Published.exact."""
    import mpmath
    import sympy as sp

    f = sp.lambdify(published.coords + published.params, published.expr, "mpmath")

    def at(*values):
        arrays = np.broadcast_arrays(*[np.asarray(v, dtype=float) for v in values])
        out = np.empty(arrays[0].shape)
        with mpmath.workdps(80):
            for index in np.ndindex(out.shape):
                try:
                    out[index] = float(f(*[mpmath.mpf(float(a[index])) for a in arrays]))
                except (ZeroDivisionError, ValueError):
                    out[index] = np.nan
        return out
    published.f = at
    return published


def curvature_views(ck, src, cv):
    """The stored Kretschmann and Ricci scalars against rho, from the regular centre out to rho = 8,
    read from the metric file and never recomputed. Checks hold both to their limits at the centre,
    K = 20/(27 a^4) and R = -2/a^2, the second the trace of the field equations there, and to falling
    toward zero far away, K from above and R from below."""
    top, gap = 8.0, 1e-3
    xaxis = cv.axis("$\\rho$", 0.0, top)
    scales = {"K": {"kind": "log", "from": 1e-6, "to": 1.0}, "R": {"kind": "linear", "from": -2.2, "to": 0.2}}
    names = {"K": "$K\\,a^4$", "R": "$R^\\mu{}_\\mu\\,a^2$"}
    panels = [{"id": field, "name": names[field], "aspect": 2, "scale": scales[field],
               "ticks": cv.value_ticks(scales[field])} for field in ("K", "R")]
    cv.shape_checks(ck, METRIC, [p["aspect"] for p in panels])
    K = in_eighty_digits(cv.nr_published(src, METRIC, SYSTEM, "kretschmann", ["rho"], list(PARAMS)))
    R = in_eighty_digits(cv.nr_published(src, METRIC, SYSTEM, "ricci_scalar", ["rho"], list(PARAMS)))
    # At the centre the terms cancel through more than sixty digits, so the limits are read at
    # rho = 1e-4 in 120, where K and R stand within 1e-8 of them.
    ck(f"{METRIC} K at the centre is 20/27 a^-4", abs(float(K.exact([1e-4, 1.0], 120)) - 20 / 27) < 1e-8)
    ck(f"{METRIC} R at the centre is -2/a^2", abs(float(R.exact([1e-4, 1.0], 120)) + 2) < 1e-7)
    xs = np.linspace(gap, top, 400)
    k, r = K(xs, 1.0), R(xs, 1.0)
    ck(f"{METRIC} K is positive and falls outward", bool(np.all(k > 0) and np.all(np.diff(k) < 0)))
    ck(f"{METRIC} R is negative and rises outward", bool(np.all(r < 0) and np.all(np.diff(r) > 0)))
    ck(f"{METRIC} both scalars tend to zero far away",
       float(K.exact([30.0, 1.0])) < 1e-20 and abs(float(R.exact([30.0, 1.0]))) < 1e-9)
    curves = []
    for name, fn in (("K", K), ("R", R)):
        curves += cv.line_curves(ck, METRIC, name, name, fn, "rho", xaxis, scales[name], cv.PLOT_W / 2, PARAMS,
                                 within=(gap, top))
    k_clip = "$K\\,a^4$ logarithmic, from $10^{-6}$ to $1$."
    r_clip = "$R^\\mu{}_\\mu\\,a^2$ linear, from $-2.2$ to $0.2$."
    return [{"id": "monopole", "label": "Monopole", "kind": "line", "system": SYSTEM, "x": xaxis,
             "panels": panels, "curves": curves, "marks": [],
             "legend": [["line", "K", "$K$"], ["line", "R", "$R^\\mu{}_\\mu$"]],
             "caption": ["The Kretschmann scalar $K$ and the Ricci scalar $R^\\mu{}_\\mu$ of the monopole against "
                         "$\\rho$, out from the regular centre. $K$ is $20/27a^4$ at the centre and falls toward zero "
                         "far away.",
                         "$R^\\mu{}_\\mu$ is $-2/a^2$ at the centre and rises toward zero far away."],
             "ricci_sentences": [[1, 0]], "settings": SETTINGS, "clip": k_clip + " " + r_clip,
             "per_scalar": cv.scalar_clips(k_clip, r_clip)}]
