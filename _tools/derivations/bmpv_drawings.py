"""The drawings of the spinning black hole of Breckenridge, Myers, Peet and Vafa.

Every drawing takes mu = 1, the unit of every area, and j = 1/2, a hole whose closed timelike curves lie
behind its horizon, j^2 < mu^3: the horizon is r = sqrt(mu) and the velocity of light surface of Gibbons
and Herdeiro, where the Hopf fibre is null, is r = j^(1/3) = 0.794 sqrt(mu), in the radius r of the Hopf,
Euler and ingoing charts. The isotropic radius is sqrt(r^2 - mu).

The plane of t and r is drawn with the Hopf fibre divided out, as Myers and Perry's hole with equal spins
is: the quotient metric is -Delta^2/(1 - j^2/r^6) c^2dt^2 + dr^2/Delta^2 with Delta = 1 - mu/r^2, whose
null curves are the projections of the null geodesics with no angular momentum, and which is positive
definite inside the velocity of light surface (Gibbons and Herdeiro, sec. 4). Its tortoise coordinate is
the shift of Gauntlett, Myers and Townsend's advanced time,

    r_*(r) = r - int_r^inf (sqrt(1 - j^2/x^6)/(1 - mu/x^2)^2 - 1) dx,   c dv = c dt + dr_*,

which runs to -infinity at the degenerate horizon, as the extreme Reissner-Nordstrom hole's does. The
exterior is drawn as Majumdar and Papapetrou's single extreme hole is, p = arctan(ct - r_*), q =
arctan(ct + r_*). bmpv.md beside this file is the derivation.
"""
import functools
import math

import numpy as np
from scipy import integrate

METRIC = "bmpv"
PARAMS = {"mu": 1, "j": "1/2"}
MU, J = 1.0, 0.5
VLS = J ** (1 / 3)
CONE = "future cone of no angular momentum"
SETTINGS = "$\\mu = 1$ and $j = \\mu^{3/2}/2$."


def slope(x):
    """dr_*/dr, the slope of the rays of no angular momentum on the quotient plane."""
    return math.sqrt(1 - J ** 2 / x ** 6) / (1 - MU / x ** 2) ** 2


SPLIT = 2.0


def _G0(x):
    """A primitive of g = x^4/(x^2 - 1)^2, the factor 1/Delta^2 at mu = 1."""
    return (-2 * x + (x * x - 1) * (4 * x + 3 * math.log(x - 1) - 3 * math.log(x + 1))) / (4 * (x * x - 1))


def _G1(x):
    """A primitive of (x - 1) g."""
    return ((x + 1) * (2 * x * x - 4 * x + math.log(x - 1) + 7 * math.log(x + 1)) + 2) / (4 * (x + 1))


FAR = 50.0


def _series(x):
    """int_x^inf (slope - 1) ds by its series far out at mu = 1, 2/s^2 + 3/s^4 + (4 - j^2/2)/s^6 + O(s^-8),
    whose next term leaves below 1e-12 beyond FAR."""
    return 2 / x + 1 / x ** 3 + (4 - J ** 2 / 2) / (5 * x ** 5)


@functools.lru_cache(maxsize=None)
def _far(x):
    """int_x^inf (slope - 1) ds for x >= SPLIT, where the integrand is regular: by quadrature up to FAR and
    by the series beyond it."""
    if x >= FAR:
        return _series(x)
    return integrate.quad(lambda s: slope(s) - 1, x, FAR, limit=400, epsabs=1e-13, epsrel=1e-12)[0] + _series(FAR)


def _tail(x):
    """int_x^inf (slope - 1) ds. Beyond SPLIT the integrand is regular. Below it slope = s g with
    s = sqrt(1 - j^2/s^6), and the double pole of g at the horizon is integrated in closed form:
    s g = (s - s1 - s1'(s - 1)) g + s1 g + s1'(s - 1) g, s1 = s(1), s1' = s'(1), whose first part is
    regular at s = 1."""
    if x >= SPLIT:
        return _far(x)
    far = _far(SPLIT)
    s1 = math.sqrt(1 - J ** 2)
    s1p = 3 * J ** 2 / s1

    def regular(s):
        return (math.sqrt(1 - J ** 2 / s ** 6) - s1 - s1p * (s - 1)) * s ** 4 / ((s - 1) * (s + 1)) ** 2
    near = (integrate.quad(regular, x, SPLIT, limit=400, epsabs=1e-12, epsrel=1e-12)[0]
            + s1 * (_G0(SPLIT) - _G0(x)) + s1p * (_G1(SPLIT) - _G1(x)) - (SPLIT - x))
    return near + far


def rstar(r):
    """r_*(r) = r - int_r^inf (slope - 1) dx, outside the horizon (mu = 1)."""
    r = np.atleast_1d(np.asarray(r, dtype=float))
    return np.array([x - _tail(x) for x in r])


def register_null(nr):
    # The Hopf chart draws no plane of its own: its angles hold the Hopf fibre as a coordinate direction
    # only on the axes theta = 0 and pi/2, where they are degenerate, and its plane of t and r with the
    # fibre divided out is the Euler chart's, the same radius and time, drawn there.
    nr.DIAGRAMS.extend([
        nr.Diagram(METRIC, "euler", "radial", "$t$ and $r$", ("t", "r"), (0.85, 3.25, -1.2, 1.2),
                   "$r/\\sqrt{\\mu}$", "$ct/\\sqrt{\\mu}$", PARAMS, {"alpha": "0", "beta": "pi/2"}, orient="ingoing",
                   quotient="gamma", cone=CONE),
        nr.Diagram(METRIC, "isotropic", "radial", "$t$ and $r$", ("t", "r"), (0, 3, -1.5, 1.5),
                   "$r/\\sqrt{\\mu}$", "$ct/\\sqrt{\\mu}$", PARAMS, {"theta": "pi/2", "phi": "0"}, orient="ingoing",
                   quotient="psi", cone=CONE),
        nr.Diagram(METRIC, "ingoing", "radial", "$v$ and $r$", ("v", "r"), (0.85, 3.25, -1.2, 1.2),
                   "$r/\\sqrt{\\mu}$", "$(cv - r)/\\sqrt{\\mu}$", PARAMS, {"theta": "pi/2", "phi": "0"},
                   to_display=nr.FINKELSTEIN_IN, tau="v - r", orient="ingoing", quotient="tildepsi", cone=CONE),
    ])
    shadow = "Its null curves are the projections of the null geodesics with no angular momentum."
    nr.CAPTIONS.update({
        (METRIC, "euler", "radial"): [
            "The plane of $t$ and $r$ ($\\beta = \\pi/2$, $\\alpha = 0$, $j = \\mu^{3/2}/2$) with the Hopf fibre "
            "$\\gamma$ divided out, each point a squashed 3-sphere. " + shadow + " The cones close at the "
            "horizon $r = \\sqrt{\\mu}$."],
        (METRIC, "isotropic", "radial"): [
            "The plane of $t$ and the isotropic $r$ ($\\theta = \\pi/2$, $\\phi = 0$, $j = \\mu^{3/2}/2$) with the "
            "Hopf fibre $\\psi$ divided out, each point a squashed 3-sphere. " + shadow + " The cones close at "
            "the horizon $r = 0$."],
        (METRIC, "ingoing", "radial"): [
            "The plane of $v$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$, $j = \\mu^{3/2}/2$) with $\\tilde\\psi$ "
            "divided out. " + shadow + " The ingoing rays cross the horizon $r = \\sqrt{\\mu}$."],
    })


# ---------------------------------------------------------------- slices: the moment t = 0

def slice_marks(slices, system):
    """The moment t = 0 of the embedding, read in the Hopf chart, on the spacetime diagram of each chart:
    t = 0 along the same r in the Hopf and Euler charts, along the isotropic sqrt(r^2 - mu), and on the
    ingoing chart along cv = r_*(r), its (cv - r, r) through the view's own display."""
    m = slices.moments(METRIC, "fibre")[0]
    lo, hi = m.reach("bmpv", "r")
    if system == "euler":
        return [slices.Mark(m, slices.along(0.0, lo, hi))]
    r = slices.near(lo, hi)
    if system == "isotropic":
        return [slices.Mark(m, [np.column_stack([np.zeros_like(r), np.sqrt(r ** 2 - MU)])])]
    return [slices.Mark(m, [np.column_stack([rstar(r), r])])]


# ---------------------------------------------------------------- the embedding diagram
#
# The moment t = 0 of the Hopf chart on the axis theta = 0, phi = 0: the plane of r and psi, whose metric
# is dr^2/Delta^2 + (r^6 - j^2)/r^4 dpsi^2, so each circle of psi, the Hopf fibre there, has the radius
# sqrt(r^6 - j^2)/r^2. Toward the degenerate horizon g_rr grows as r^4/(r^2 - mu)^2 and the throat runs on
# without end, its circles closing on sqrt(mu^3 - j^2)/mu = sqrt(3)/2 sqrt(mu).

EMBED_LO, EMBED_TOP = 1.01, 4.0


def circle_radius(r):
    r = np.asarray(r, dtype=float)
    return np.sqrt(r ** 6 - J ** 2) / r ** 2


def embedding_views(ck, src, em):
    sl = em.Slice(src, METRIC, "bmpv", "r", "\\psi", {"t": 0, "theta": 0, "phi": 0}, PARAMS)
    name = "BMPV black hole, the Hopf fibre"
    lo, top = EMBED_LO, EMBED_TOP
    size = 2 * float(circle_radius(top))
    end = math.sqrt(MU ** 3 - J ** 2) / MU
    radii = (1.1, 1.5, 2.0, 3.0)
    throat = em.Piece("throat", "sheet", sl, lo, top, 0.0, 1,
                      (("edge", "the throat, running on toward the horizon $r = \\sqrt{\\mu}$"),
                       ("edge", "the surface runs on to $r \\to \\infty$")),
                      [(lo, "r", None)] + [(r, "r", None) for r in radii] + [(top, "r", None)], size)
    ck.isometry(name, throat)
    ck.radius(f"{name}, the circle radius sqrt(r^6 - j^2)/r^2", throat, circle_radius, size)
    ck.add(f"{name}: the circles close on sqrt(mu^3 - j^2)/mu", abs(float(circle_radius(1 + 1e-12)) - end), 1e-9)
    surface = em.Surface([throat])
    fig = em.figure_of([surface], {"sheet": "cover"}, size)
    em.ring_label(fig, [0, 0, 0], *throat.at(lo), "$r = 1.01\\sqrt{\\mu}$")
    em.ring_label(fig, [0, 0, 0], *throat.at(top), "$4\\sqrt{\\mu}$")
    fig.legend("fill", "cover", "one moment of $t$")
    fig.legend("line", "r", "$r$ constant, from $1.01$ to $4\\sqrt{\\mu}$")
    fig.legend("line", "meridian", "$\\psi$ constant, every $15°$")
    return [em.view("fibre", "The Hopf fibre", "$\\sqrt{\\mu}$", [surface], fig.done(),
                    settings=SETTINGS)]


def embedding_captions():
    return {(METRIC, "fibre"): [
        "The plane of $r$ and $\\psi$ on the axis $\\theta = 0$ at one moment of $t$, each circle the Hopf fibre, of "
        "radius $\\sqrt{r^6 - j^2}/r^2$. An infinitely long throat, drawn down to $r = 1.01\\sqrt{\\mu}$."]}


# ---------------------------------------------------------------- the conformal diagram

def conformal_views(ck, src, cf):
    plane = cf.Plane(src, METRIC, "euler", ("t", "r"), {"alpha": "0", "beta": "pi/2"}, PARAMS, quotient="gamma")
    name = "BMPV black hole, exterior"
    box = [-cf.PI - 0.25, cf.PI + 0.25, -cf.PI - 0.25, cf.PI + 0.25]
    t = cf.spread(-np.inf, np.inf, 500, 10)

    def E(tt, x):
        rs = rstar(x)
        tt = np.asarray(tt, dtype=float)
        return np.arctan(tt - rs), np.arctan(tt + rs)

    ck.chart(f"{name}", plane, E, ck.uniform(-20, 20, 400), ck.uniform(1.02, 40, 400), lambda tt, x: (1, 0))
    g00, _, g11 = plane.metric(np.zeros(4), np.array([1.1, 1.5, 2.0, 5.0]))[:3]
    ck.limit(f"{name}: the quotient metric is -Delta^2/(1 - j^2/r^6) dt^2 + dr^2/Delta^2",
             [float(g00[k] / -((1 - 1 / x ** 2) ** 2 / (1 - J ** 2 / x ** 6))) for k, x in enumerate((1.1, 1.5, 2.0, 5.0))]
             + [float(g11[k] * (1 - 1 / x ** 2) ** 2) for k, x in enumerate((1.1, 1.5, 2.0, 5.0))], [1.0] * 8, 1e-10)
    pp, qq = E(np.array([-6.0, 0.0, 6.0]), np.full(3, 1 + 1e-6))
    ck.limit(f"{name}: r -> sqrt(mu) lands on the horizon p = pi/2", pp, [cf.HALF] * 3, 1e-4)
    pp, qq = E(np.array([0.0]), np.array([1e6]))
    ck.limit(f"{name}: r -> infinity at t = 0 lands on i0, (X, T) = (pi, 0)", cf.point(pp[0], qq[0]), [cf.PI, 0], 1e-3)
    ck.finite(f"{name}: the Kretschmann scalar is finite toward the horizon",
              plane.kretschmann(np.zeros(3), np.array([1 + 1e-4, 1 + 1e-3, 1 + 1e-2])))
    v = cf.View("exterior", "The exterior", box, "euler")
    diamond = [[0, -cf.PI], [-cf.PI, 0], [0, cf.PI], [cf.PI, 0]]
    v.fill("region", diamond)
    v.fill("cover", diamond)
    r_out = (1.1, 1.5, 2.0, 3.0)
    for x in r_out:
        v.curve("r", *E(t, np.full_like(t, x)))
    xs = 1 + cf.spread(0.0, np.inf, 600, 16)
    for tt in (-3.0, -1.0, 0.0, 1.0, 3.0):
        v.curve("t", *E(np.full_like(xs, tt), xs))
    v.line("horizon", [[[-cf.PI, 0], [0, cf.PI]], [[0, -cf.PI], [-cf.PI, 0]]])
    v.line("scri", [[[0, cf.PI], [cf.PI, 0]], [[cf.PI, 0], [0, -cf.PI]]])
    for at, text, anchor, dx in (((cf.PI, 0), "$i^0$", "l", 6), ((0, cf.PI), "$i^+$", "l", 6),
                                 ((0, -cf.PI), "$i^-$", "l", 6)):
        v.layers.append({"kind": "point", "class": "infinity", "at": [round(at[0], 4), round(at[1], 4)]})
        v.label_xt(list(at), text, anchor, dx=dx)
    v.label_xt([3 * cf.Q4 * 2 / 3 + 0.8, cf.HALF], "$\\mathscr{I}^+$", "bl", dx=4, dy=-3)
    v.label_xt([3 * cf.Q4 * 2 / 3 + 0.8, -cf.HALF], "$\\mathscr{I}^-$", "tl", dx=4, dy=3)
    v.label_xt([-cf.HALF, cf.HALF], "$r = \\sqrt{\\mu}$", "br", "small", dx=-5, dy=-1)
    v.label_xt([cf.HALF, 0.45], "exterior", cls="region")
    v.legend("cover", "the exterior $r > \\sqrt{\\mu}$")
    v.legend("r", "$r$ constant, at $1.1$, $1.5$, $2$, and $3\\sqrt{\\mu}$")
    v.legend("t", "$t$ constant")
    v.legend("horizon", "the horizon $r = \\sqrt{\\mu}$")
    v.legend("scri", "null infinity $\\mathscr{I}^\\pm$")
    moment = cf.slices.moments(METRIC, "fibre")[0]
    lo, hi = moment.reach("bmpv", "r")
    rr = np.geomspace(lo, hi, 400)
    v.slice(moment, [E(np.zeros_like(rr), rr)])
    v.set(settings=SETTINGS,
          restriction="The plane $\\beta = \\pi/2$, $\\alpha = 0$ with $\\gamma$ divided out.")
    return [v]


def conformal_captions():
    return {(METRIC, "exterior"): [
        "The exterior with the Hopf fibre divided out, $p = \\arctan(ct - r_*)$ and $q = \\arctan(ct + r_*)$, where "
        "$dr_*/dr = \\sqrt{1 - j^2/r^6}/(1 - \\mu/r^2)^2$. The degenerate horizon lies on the two left edges."]}


# ---------------------------------------------------------------- the curvature figure
#
# The published Kretschmann and Ricci scalars of the ingoing chart against r, over the range the chart
# covers at j = 1/2, from the velocity of light surface out: both are finite at the horizon,
# K = 50.5/mu^2 and R = -1/mu there, and R = 2(2j^2 - mu^2 r^2)/r^8 is negative over the whole range.

def curvature_views(ck, src, cv):
    params = {"\\mu": MU, "j": J}
    values = [MU, J]
    K = cv.nr_published(src, METRIC, "ingoing", "kretschmann", ["r"], ["\\mu", "j"])
    R = cv.nr_published(src, METRIC, "ingoing", "ricci_scalar", ["r"], ["\\mu", "j"])
    x_axis = cv.axis("$r/\\sqrt{\\mu}$", 0.8, 4.0)
    K_scale = {"kind": "log", "from": 1e-3, "to": 1e3}
    R_scale = {"kind": "symlog", "from": -10.0, "to": 0.0, "linear": 1e-4}
    panels = [{"id": "K", "name": "$K\\,\\mu^2$", "aspect": 2, "scale": K_scale, "ticks": cv.value_ticks(K_scale)},
              {"id": "R", "name": "$R\\,\\mu$", "aspect": 2, "scale": R_scale, "ticks": cv.value_ticks(R_scale)}]
    height = {p["id"]: cv.PLOT_W / p["aspect"] for p in panels}
    cv.shape_checks(ck, METRIC, [p["aspect"] for p in panels])
    curves = (cv.line_curves(ck, METRIC, "K", "K", K, "r", x_axis, K_scale, height["K"], params)
              + cv.line_curves(ck, METRIC, "R", "R", R, "r", x_axis, R_scale, height["R"], params))
    ck(f"{METRIC} K at the horizon is 50.5/mu^2", abs(float(K.exact([1.0] + values)) - 50.5) < 1e-10)
    ck(f"{METRIC} R at the horizon is -1/mu", abs(float(R.exact([1.0] + values)) + 1.0) < 1e-12)
    rs = np.linspace(0.8, 4.0, cv.SAMPLES)
    ck(f"{METRIC} R is 2(2j^2 - mu^2 r^2)/r^8 and negative on the range drawn",
       np.allclose(R(rs, *values), 2 * (2 * J ** 2 - rs ** 2) / rs ** 8, rtol=1e-10, atol=0)
       and bool(np.all(R(rs, *values) < 0)))
    ck(f"{METRIC} K is positive and falls from the velocity of light surface outward",
       bool(np.all(K(rs, *values) > 0)) and bool(np.all(np.diff(K(rs, *values)) < 0)))
    ck(f"{METRIC} both curves stay inside their panels",
       not any(cv.beyond(K, x_axis, K_scale, values)) and not any(cv.beyond(R, x_axis, R_scale, values)))
    return [{
        "id": "radial", "label": "against $r$", "kind": "line", "x": x_axis, "panels": panels,
        "curves": curves, "marks": [{"class": "horizon", "x": 1.0}],
        "legend": [["line", "K", "$K$"], ["line", "R", "$R$"], ["line", "horizon", "the horizon"]],
        "caption": [
            "The Kretschmann scalar $K$ and the Ricci scalar $R$ of the ingoing chart against $r$ "
            "($j = \\mu^{3/2}/2$) outside the velocity of light surface. The Kretschmann scalar "
            "falls outward, $50.5/\\mu^2$ at the horizon.",
            "The Ricci scalar is $R = 2(2j^2 - \\mu^2r^2)/r^8$, $-1/\\mu$ at the horizon."],
        "ricci_sentences": [[1, 0]],
        "per_scalar": cv.scalar_clips("Logarithmic, $10^{-3}$ to $10^3$.",
                                      "Signed logarithmic, linear within $10^{-4}$ of zero, $-10$ to $0$."),
        "settings": SETTINGS,
        "clip": "Logarithmic $K$ and signed logarithmic $R$.",
    }]
