"""The drawings of the scalar field universe of Husain, Martinez and Nunez and of Fonarev's family.

Every drawing takes b = 1, the unit of every length, and c = 1. The two massless members,
gamma = +sqrt(3)/2 and gamma = -sqrt(3)/2, are drawn in the expanding branch with eta_0 = 2b, where
the scale factor a = (eta/eta_0)^(1/2) is one, so the comoving chart has t_0 = 2 eta_0/3c = 4b/3c.
Fonarev's family is drawn at gamma = 4/5, where s = (1 - gamma^2)/(2 gamma^2 - 1) = 9/7 and
n = gamma^2/(1 - gamma^2) = 16/9, with k = c/b, so that T = 0 is eta = eta_0 = s c/k = 9b/7.

On the plane of the time and r every chart is a factor times f^gamma (-d eta^2 + dr_*^2), with
f = 1 - b/r and the tortoise radius

    r_* = b f^(1 - gamma) 2F1(2, 1 - gamma; 2 - gamma; f)/(1 - gamma),

the integral of f^(-gamma) dr from r = b, finite at r = b for every |gamma| < 1 and without bound
as r grows. So the radial light rays keep eta + r_* or eta - r_*, and p, q = arctan((eta -+ r_*)/l)
bring the expanding branch eta > 0 into the triangle of the half of Minkowski's plane with t > 0:
the big bang eta = 0 on T = 0, the singularity r = b on X = 0, and null infinity from i^+ to i^0.
"""
import math

import numpy as np
from scipy.special import hyp2f1


METRIC = "scalar_collapse_cosmology"
ROOT = math.sqrt(3) / 2
ETA0 = 2.0
T0 = 2 * ETA0 / 3
# Fonarev's family as the drawings take it.
FON_GAMMA, FON_S, FON_N, FON_K = 0.8, 9 / 7, 16 / 9, 1.0
FON_ETA0 = FON_S / FON_K
SIGNS = {"positive": ROOT, "negative": -ROOT}
GAMMA_TEXT = {"positive": "\\sqrt{3}/2", "negative": "-\\sqrt{3}/2"}


def params(system, sign="positive"):
    """The drawing's parameter values in the chart's plain names, as strings sympy reads."""
    gamma = "sqrt(3)/2" if sign == "positive" else "-sqrt(3)/2"
    if system == "conformal_time":
        return {"b": "1", "eta_0": "2", "gamma": gamma}
    if system == "comoving_time":
        return {"b": "1", "t_0": "4/3", "gamma": gamma}
    if system == "fonarev_conformal":
        return {"b": "1", "eta_0": "9/7", "gamma": "4/5", "s": "9/7"}
    return {"b": "1", "k": "1", "gamma": "4/5", "n": "16/9"}


def rstar(r, gamma):
    """The tortoise radius from r = b = 1, the integral of f^(-gamma) dr."""
    r = np.asarray(r, dtype=float)
    f = (r - 1) / r
    return f ** (1 - gamma) * hyp2f1(2.0, 1 - gamma, 2 - gamma, f) / (1 - gamma)


def eta_of(system, x0):
    """The conformal time of the chart's x^0, in units of b."""
    x0 = np.asarray(x0, dtype=float)
    if system == "comoving_time":
        return ETA0 * np.cbrt(x0 / T0) ** 2
    if system == "fonarev":
        return FON_ETA0 * np.exp(FON_K * x0 / FON_S)
    return x0


def horizon(r, gamma, s=0.5):
    """The apparent horizon, eta = 2 s r^2 f^(1 - gamma)/(2r - b(1 + gamma)): Husain, Martinez and
    Nunez's (11) at s = 1/2, where the gradient of the areal radius a r f^((1 - gamma)/2) is null."""
    r = np.asarray(r, dtype=float)
    f = (r - 1) / r
    return 2 * s * r ** 2 * f ** (1 - gamma) / (2 * r - (1 + gamma))


# (system, view, sign, box, ylabel): each chart's plane of its time and r.
PLANES = (
    ("conformal_time", "positive", "positive", (1, 5, 0, 4), "$\\eta/b$"),
    ("conformal_time", "negative", "negative", (1, 5, 0, 4), "$\\eta/b$"),
    ("comoving_time", "radial", "positive", (1, 5, 0, 4), "$ct/b$"),
    ("fonarev_conformal", "radial", None, (1, 7, 0, 6), "$\\eta/b$"),
    ("fonarev", "radial", None, (1, 5, -2, 2), "$cT/b$"),
)
CLOCK = {"conformal_time": "\\eta", "comoving_time": "t", "fonarev_conformal": "\\eta", "fonarev": "T"}
TAU = {"conformal_time": "eta", "comoving_time": "t", "fonarev_conformal": "eta", "fonarev": "T"}
LABEL = {("conformal_time", "positive"): "$\\gamma = \\sqrt{3}/2$",
         ("conformal_time", "negative"): "$\\gamma = -\\sqrt{3}/2$",
         ("comoving_time", "radial"): "$t$ and $r$",
         ("fonarev_conformal", "radial"): "$\\eta$ and $r$",
         ("fonarev", "radial"): "$T$ and $r$"}


def gamma_of(system, sign):
    return FON_GAMMA if system in ("fonarev_conformal", "fonarev") else SIGNS[sign]


def null_rows(nr):
    """The rows, the captions and the closed forms of the null rays of every plane."""
    rows, captions, forms = [], {}, {}
    for system, vid, sign, box, ylabel in PLANES:
        clock = CLOCK[system]
        # Beneath the big bang the chart runs on into the time reverse, where its signature turns over and
        # the two families of rays trade orientation, so each plane that reaches eta = 0 declares the time
        # positive in `where` and every ray ends there. At gamma = 4/5 the Kretschmann scalar passes 10^8
        # beside r = b too slowly for the test at 10^-5 of the drawing, so Fonarev's rows declare r = b,
        # which the generator checks in sixty digits.
        extra = {"conformal_time": {"where": "eta"}, "comoving_time": {"where": "t"},
                 "fonarev_conformal": {"where": "eta", "singular_zero": "r - b"},
                 "fonarev": {"singular_zero": "r - b"}}[system]
        rows.append(nr.Diagram(METRIC, system, vid, LABEL[(system, vid)], (clock, "r"), box, "$r/b$", ylabel,
                               params(system, sign or "positive"), nr.EQUATOR, tau=TAU[system], areal=True, **extra))
        captions[(METRIC, system, vid)] = [CAPTIONS[(system, vid)]]
        gamma = gamma_of(system, sign or "positive")
        forms[(METRIC, system, vid)] = (
            (lambda x0, r, s=system, g=gamma: eta_of(s, x0) + rstar(r, g)),
            (lambda x0, r, s=system, g=gamma: eta_of(s, x0) - rstar(r, g)),
            # Along r = b the rays stand nearly upright where r_* is steep, and beside t = 0 the conformal
            # time eta = eta_0 (t/t_0)^(2/3) is steep in t, so that a step of the tracing moves either
            # by more than the check allows; the forms are held beyond r = 1.01 b and t = 0.1 b/c.
            (lambda x0, r: (r > 1.01) & (x0 > 0.1)) if system == "comoving_time" else (lambda x0, r: r > 1.01))
    return rows, captions, forms


CAPTIONS = {
    ("conformal_time", "positive"):
        "The plane of $\\eta$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) at $\\gamma = \\sqrt{3}/2$, with the radial "
        "null geodesics of Fisher's metric and the apparent horizon dotted. Between $\\eta = 1.67\\,b$ and "
        "$3.58\\,b$ the horizon cuts each moment three times, once within $0.013\\,b$ of the singularity $r = b$.",
    ("conformal_time", "negative"):
        "The plane of $\\eta$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) at $\\gamma = -\\sqrt{3}/2$, with the radial "
        "null geodesics of Fisher's metric and one apparent horizon, dotted, which leaves the big bang at $r = b$ and "
        "grows as $r \\approx 2\\eta$.",
    ("comoving_time", "radial"):
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) at $\\gamma = \\sqrt{3}/2$, in the proper time of "
        "comoving observers far from the centre, $\\eta = \\eta_0(t/t_0)^{2/3}$, with the apparent horizon dotted.",
    ("fonarev_conformal", "radial"):
        "The plane of $\\eta$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) of Fonarev's solution at $\\gamma = 4/5$, "
        "where $a = (\\eta/\\eta_0)^{9/7}$, with the radial null geodesics of Fisher's metric at that exponent and "
        "the apparent horizon dotted.",
    ("fonarev", "radial"):
        "The plane of $T$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) of Fonarev's solution at $\\gamma = 4/5$, where "
        "$\\eta = \\eta_0e^{7kT/9}$, with the apparent horizon dotted; its big bang is $T \\to -\\infty$.",
}


def slice_lines(system, vid):
    """The moment eta = eta_0 of the embedding view of the plane's own exponent, as a line of the
    plane's time across the embedding's reach in r."""
    def lines(m):
        lo, hi = m.reach("conformal_time", "r")
        clock = T0 if system == "comoving_time" else ETA0
        return [[(clock, lo), (clock, hi)]]
    return lines


def flat_marks(slices):
    """The flat views that mark a moment, each the moment of its own exponent."""
    label = {"comoving_time": "$t = t_0$", "conformal_time": "$\\eta = \\eta_0$"}
    out = {}
    for system, vid, sign, _, _ in PLANES:
        if system in ("conformal_time", "comoving_time"):
            out[(METRIC, system, vid)] = (lambda system=system, vid=vid, sign=sign: slices.one(
                METRIC, slice_lines(system, vid), label=label[system], view_id=sign))
    return out


HIDDEN_REASON = ("Fonarev's solution at gamma = 4/5, with an exponential potential, another member of the family "
                 "than the two massless members whose moments are embedded")


def level_radius(gamma):
    """Where g_rr = (d rho/dr)^2 on a moment: rho = a r f^((1 - gamma)/2) gives
    g_rr - (d rho/dr)^2 = a^2 f^(-1 - gamma)(gamma/r - (1 + gamma)^2/4r^2), zero at r_e = (1 + gamma)^2 b/4 gamma,
    which lies outside r = b only for gamma > 0."""
    return (1 + gamma) ** 2 / (4 * gamma)


def defect(r, gamma):
    """g_rr - (d rho/dr)^2 at a = 1 and b = 1, in the closed form above, with f written (r - 1)/r."""
    r = np.asarray(r, dtype=float)
    return ((r - 1) / r) ** (-1 - gamma) * (gamma / r - (1 + gamma) ** 2 / (4 * r * r))


def areal(r, gamma):
    r = np.asarray(r, dtype=float)
    return r * ((r - 1) / r) ** ((1 - gamma) / 2)


def height(lo, hi, gamma, space):
    """The rise of the profile from lo to hi, a quadrature of sqrt(+-defect) in y = ln(r - b), which
    spreads the decades next to the singularity evenly; it never reads the published metric."""
    from scipy.integrate import quad
    sign = 1 if space == "flat" else -1

    def rate(y):
        r = 1 + math.exp(y)
        return math.sqrt(max(sign * float(defect(r, gamma)), 0.0)) * math.exp(y)
    return quad(rate, math.log(lo - 1), math.log(hi - 1), epsabs=1e-13, epsrel=1e-12, limit=400)[0]


# Where each profile stops short of the singular point: toward it the surface closes on the light cone,
# a chord's interval being the small difference of two nearly equal lengths, about 15 (r - b)^(1/2) of
# its extent in rho at gamma = sqrt(3)/2 and (r - b)^(1/2) at gamma = -sqrt(3)/2, which fourteen decimals
# hold to the checks' part in 5000 down to these.
NEAR_TIP, WHOLE_TIP = 1e-7, 6e-4


def embedding_views(ck, src, em):
    """The equatorial plane of the moment eta = eta_0, a = 1, of each massless member, read in the
    conformal chart from the singularity out to r = 4b.

    At gamma = sqrt(3)/2 the plane lies level on r_e = (1 + gamma)^2 b/4 gamma = 1.00518 b. Beyond it the
    surface is drawn in flat space, climbing from level as Fisher's does; inside it the circles grow
    faster than the distance out to them, which is checked, and the surface is drawn in three
    dimensional Minkowski space, falling from level toward the singular point along the light cone.
    Its circles close as (r - b)^0.067, so the profile stops at r = b + 10^-14 b, a circle of radius
    0.115 b, the last that doubles hold. At gamma = -sqrt(3)/2 the circles grow faster than the
    distance out to them everywhere, which is checked, and the whole plane is drawn in Minkowski
    space, from the singular point, which it leaves along the light cone, out to r = 4b."""
    import sympy as sp
    views = []
    fixed = {"eta": "2", "theta": "pi/2"}
    rewrite = lambda e, x: e.subs(1 - 1 / x, (x - 1) / x)
    top = 4.0
    for vid in ("positive", "negative"):
        gamma = SIGNS[vid]
        P = params("conformal_time", vid)
        name = f"Husain-Martinez-Nunez, gamma = {'' if gamma > 0 else '-'}sqrt(3)/2"
        msl = em.Slice(src, METRIC, "conformal_time", "r", "\\phi", fixed, P, space="minkowski", rewrite=rewrite)
        sl = em.Slice(src, METRIC, "conformal_time", "r", "\\phi", fixed, P, rewrite=rewrite)
        size = 2 * float(areal(top, gamma))
        rings = [(1.5, "r", None), (2.0, "r", None), (3.0, "r", None), (top, "r", None)]
        if gamma > 0:
            level, tip = level_radius(gamma), 1 + NEAR_TIP
            exact = (1 + sp.sqrt(3) / 2) ** 2 / (2 * sp.sqrt(3))
            sl.known = msl.known = {level: exact}
            ck.add(f"{name}: the surface lies level where g_rr = (d rho/dr)^2, at r = (1 + gamma)^2 b/4 gamma, "
                   "as a part of g_rr", abs(float(sl.defect_at(np.array([level]))[0]) / float(sl.gxx_at(level))), 1e-12)
            ck.stops(f"{name}, inside r_e in flat space", sl, 1 + np.geomspace(1e-14, level - 1, 400)[:-1])
            outward = np.linspace(level, 40, 400)[1:]
            ck.add(f"{name}: beyond r_e no surface in Minkowski space carries the slice",
                   float(max(0.0, np.max(-msl.defect_at(outward)))), 0.0)
            if not np.all(msl.defect_at(outward) > 0):
                ck.items[-1]["ok"] = False
            join = ("at $r = 1.005\\,b$ the surface lies level, in Minkowski space nearer the singularity and in flat "
                    "space beyond")
            far = em.Piece("far", "sheet", sl, level, top, 0.0, 1,
                           (("join", join), ("edge", "the surface runs on, to $r \\to \\infty$")),
                           [(level, "space", None)] + rings, size, digits=1e-12)
            # A two hundredth of rho at a time toward the tip, where rho grows as (r - b)^0.067.
            knots = [1 + (level - 1) * 1.08 ** -k for k in range(1, 400) if (level - 1) * 1.08 ** -k > NEAR_TIP]
            drop = msl.rise(tip, level)
            near = em.Piece("near", "sheet", msl, tip, level, -drop, 1,
                            (("edge", "the surface runs on into the axis along a light cone, to the singularity "
                                      "$r = b$"), ("join", join)),
                            [(level, "space", None)], size, digits=1e-14, knots=knots)
            pieces = [near, far]
            ck.join(f"{name}, near in Minkowski space and far in flat space at r_e", near, level, far, level)
            ck.form(f"{name}, far, z = the integral of sqrt(g_rr - (d rho/dr)^2)", far,
                    lambda xs: np.array([height(level, x, gamma, "flat") if x > level else 0.0 for x in xs]), size)
            ck.form(f"{name}, near, Z = -the integral of sqrt((d rho/dr)^2 - g_rr)", near,
                    lambda xs: np.array([-height(x, level, gamma, "minkowski") if x < level else 0.0 for x in xs]),
                    size)
        else:
            tip = 1 + WHOLE_TIP
            ck.add(f"{name}: no point of the plane has g_rr >= (d rho/dr)^2, so all of it is drawn in Minkowski space",
                   float(max(0.0, np.max(sl.defect_at(1 + np.geomspace(1e-12, 39, 600))))), 0.0)
            # A fiftieth of rho at a time toward the tip, where rho grows as (r - b)^0.933.
            knots = [1 + 1.02 ** -k for k in range(1, 800) if 1.02 ** -k > WHOLE_TIP]
            whole = em.Piece("whole", "sheet", msl, tip, top, 0.0, 1,
                             (("edge", "the surface leaves the singular point $r = b$ along a light cone"),
                              ("edge", "the surface runs on, to $r \\to \\infty$")),
                             rings, size, digits=1e-14, knots=knots)
            pieces = [whole]
            ck.form(f"{name}, Z = the integral of sqrt((d rho/dr)^2 - g_rr)", whole,
                    lambda xs: np.array([height(tip, x, gamma, "minkowski") if x > tip else 0.0 for x in xs]), size)
        for p in pieces:
            space = "in Minkowski space" if p.sl.lorentz else "in flat space"
            ck.isometry(f"{name}, {p.id} {space}", p)
            ck.radius(f"{name}, {p.id}, rho = a r f^((1 - gamma)/2) {space}", p, lambda r: areal(r, gamma), size)
        surface = em.Surface(pieces)
        fig = em.figure_of([surface], {"sheet": "cover"}, size)
        em.ring_label(fig, [0, 0, 0], *pieces[-1].at(2.0), "$2b$")
        em.ring_label(fig, [0, 0, 0], *pieces[-1].at(top), "$4b$")
        fig.legend("fill", "cover", "the plane $\\eta = \\eta_0$")
        fig.legend("line", "r", "$r = 1.5\\,b$, $2b$, $3b$, $4b$")
        if gamma > 0:
            fig.legend("line", "space", "the level circle $r = 1.005\\,b$")
        fig.legend("line", "meridian", "$\\phi$ every $15°$")
        settings = f"$\\gamma = {GAMMA_TEXT[vid]}$, $b = 1$, $\\eta_0 = 2b$, $a = 1$; "
        settings += ("inside $r = 1.005\\,b$ lengths are $dX^2 + dY^2 - dZ^2$." if gamma > 0
                     else "lengths are $dX^2 + dY^2 - dZ^2$.")
        views.append(em.view(vid, f"$\\gamma = {GAMMA_TEXT[vid]}$", "$b$", [surface], fig.done(), settings=settings))
    return views


EMBEDDING_CAPTIONS = {
    "positive": "The equatorial plane at $\\eta = \\eta_0$ ($\\gamma = \\sqrt{3}/2$), in flat space outside "
                "$r = 1.005\\,b$ and in Minkowski space inside, where it closes on a light cone.",
    "negative": "The equatorial plane at $\\eta = \\eta_0$ ($\\gamma = -\\sqrt{3}/2$), in Minkowski space from the "
                "singularity $r = b$ out.",
}


# (system, view, exponent's sign): one conformal view for each chart, both massless members in the
# conformal chart.
CONFORMAL = (("conformal_time", "positive", "positive"), ("conformal_time", "negative", "negative"),
             ("comoving_time", "comoving", "positive"), ("fonarev_conformal", "fonarev_conformal", None),
             ("fonarev", "fonarev", None))
ELL = 4.0
# Each chart's own time, at the moments its coordinate lines are drawn, and where its big bang is.
TIMES = {"conformal_time": (0.5, 1, 2, 4, 8), "comoving_time": (0.5, 1, 2, 4, 8),
         "fonarev_conformal": (0.5, 1, 2, 4, 8), "fonarev": (-2, -1, 0, 1, 2)}
BANG = {"conformal_time": 0.0, "comoving_time": 0.0, "fonarev_conformal": 0.0, "fonarev": -400.0}
MIDDLE = {"conformal_time": ETA0, "comoving_time": T0, "fonarev_conformal": FON_ETA0, "fonarev": 0.0}
TIME_TEXT = {"conformal_time": "$\\eta = b/2$, $b$, $2b$, $4b$, $8b$", "comoving_time": "$ct = b/2$, $b$, $2b$, $4b$, $8b$",
             "fonarev_conformal": "$\\eta = b/2$, $b$, $2b$, $4b$, $8b$", "fonarev": "$cT = 0$, $\\pm b$, $\\pm 2b$"}
CONFORMAL_LABEL = {"positive": "$\\gamma = \\sqrt{3}/2$", "negative": "$\\gamma = -\\sqrt{3}/2$",
                   "comoving": "Comoving time", "fonarev_conformal": "Fonarev, $\\eta$", "fonarev": "Fonarev, $T$"}


def conformal_views(ck, src, cf):
    """Each chart's plane carried by p, q = arctan((eta -+ r_*)/l), l = 4b, onto the triangle
    T = p + q > 0, X = q - p > 0 of Minkowski's compactification, each point a 2-sphere: the big bang
    eta = 0 is its base T = 0 and r = b its side X = 0, both curvature singularities where the
    Kretschmann scalar diverges, which is checked, and null infinity runs from i^+ at (0, pi) to i^0
    at (pi, 0). The apparent horizon is Husain, Martinez and Nunez's (11) for the massless members
    and its generalisation to a = (eta/eta_0)^s for Fonarev's, both of which print_charts.py checks,
    carried by the same map."""
    views = []
    pi, half = cf.PI, cf.HALF
    for system, vid, sign in CONFORMAL:
        gamma = gamma_of(system, sign or "positive")
        s = FON_S if system.startswith("fonarev") else 0.5
        plane = cf.Plane(src, METRIC, system, (CLOCK[system], "r"), cf.EQUATOR, params(system, sign or "positive"))

        def fmap(x0, r, system=system, gamma=gamma):
            return cf.mink_pq(eta_of(system, x0), rstar(r, gamma), ELL)

        name = f"Husain-Martinez-Nunez {system}" + (f", gamma = {'' if gamma > 0 else '-'}sqrt(3)/2" if sign else "")
        x0s = ck.uniform(-6, 4) if system == "fonarev" else ck.uniform(0.02, 20)
        ck.chart(name, plane, fmap, x0s, 1 + np.exp(ck.uniform(-6, 3.5)), lambda x0, r: (1, 0))
        K = plane.kretschmann
        mid = np.full(3, MIDDLE[system])
        ck.diverges(f"{name}: the Kretschmann scalar diverges on r = b",
                    K(mid, np.full(3, 1 + 1e-4)), K(mid, np.full(3, 1 + 1e-5)))
        rr = np.array([1.5, 2.0, 4.0])
        near, nearer = ((-6.0, -6.0 - 9 * math.log(10) / 7) if system == "fonarev" else
                        (0.01, 0.001) if system != "comoving_time" else (0.01, 0.01 * 10 ** -1.5))
        ck.diverges(f"{name}: the Kretschmann scalar diverges on the big bang",
                    K(np.full(3, near), rr), K(np.full(3, nearer), rr))
        ck.limit(f"{name}: the big bang meets r = b at (X, T) = (0, 0)",
                 cf.point(*fmap(BANG[system], 1.0)), [0, 0], 1e-9)
        ck.limit(f"{name}: the big bang runs to i0 at (pi, 0)", cf.point(*fmap(BANG[system], 1e9)), [pi, 0], 1e-6)
        ck.limit(f"{name}: r = b is the side X = 0", [cf.point(*fmap(x, 1.0))[0] for x in TIMES[system]],
                 [0] * 5, 1e-12)
        later = {"fonarev": 60.0, "comoving_time": 1e12}.get(system, 1e9)
        ck.limit(f"{name}: at a fixed r the end of time is i+ at (0, pi)", cf.point(*fmap(later, 3.0)), [0, pi], 1e-6)

        v = cf.View(vid, CONFORMAL_LABEL[vid], [-0.35, pi + 0.35, -0.3, pi + 0.3], system)
        tri = [[0, 0], [pi, 0], [0, pi]]
        v.fill("region", tri)
        v.fill("cover", tri)
        out = cf.spread(1, np.inf, 600, 14)
        for x0 in TIMES[system]:
            v.curve("t", *fmap(np.full_like(out, x0), out))
        own = cf.S_ALL if system == "fonarev" else cf.spread(0, np.inf, 600, 12)
        for r in (1.5, 2.0, 4.0):
            v.curve("r", *fmap(own, np.full_like(own, r)))
        rh = cf.spread(1, np.inf, 1500, 16)
        v.curve("apparent", *cf.mink_pq(horizon(rh, gamma, s), rstar(rh, gamma), ELL))
        v.line("singular", [[[pi, 0], [0, 0], [0, pi]]], zig=True)
        v.line("scri", [[[0, pi], [pi, 0]]])
        for at, text, anchor, dx, dy in (((0, pi), "$i^+$", "b", 0, -6), ((pi, 0), "$i^0$", "l", 6, 0)):
            v.layers.append({"kind": "point", "class": "infinity", "at": [round(at[0], 4), round(at[1], 4)]})
            v.label_xt(at, text, anchor, dx=dx, dy=dy)
        v.label_xt([half, half], "$\\mathscr{I}^+$", "bl", dx=5, dy=-3)
        v.label_xt([0, half], "$r = b$", "r", dx=-6)
        v.legend("cover", "the expanding branch")
        v.legend("t", TIME_TEXT[system])
        v.legend("r", "$r = 1.5\\,b$, $2b$, $4b$")
        v.legend("apparent", "the apparent horizon")
        v.legend("singular", "the big bang and $r = b$")
        v.legend("scri", "future null infinity $\\mathscr{I}^+$")
        if sign:
            moment = cf.slices.moments(METRIC, view_id=sign)[0]
            lo, hi = moment.reach("conformal_time", "r")
            rr = lo + (hi - lo) * np.linspace(0, 1, 600) ** 3
            v.slice(moment, [cf.mink_pq(np.full_like(rr, ETA0), rstar(rr, gamma), ELL)],
                    label="$t = t_0$" if system == "comoving_time" else "$\\eta = \\eta_0$")
        tail = "$p, q = \\arctan((\\eta \\mp r_*)/\\ell)$ with $r_* = \\int_b^r f^{-\\gamma}dr$."
        if system == "fonarev":
            settings = "$b = 1$, $k = c/b$, $\\ell = 4b$; " + tail[:-1] + " and $\\eta = \\eta_0e^{7kT/9}$."
        elif system == "fonarev_conformal":
            settings = "$b = 1$, $\\eta_0 = 9b/7$, $\\ell = 4b$; " + tail
        else:
            settings = "$b = 1$, $\\eta_0 = 2b$, $\\ell = 4b$; " + tail
        v.set(settings=settings)
        views.append(v)
    return views


CONFORMAL_CAPTIONS = {
    "positive": "The expanding universe of Husain, Martinez, and Núñez ($\\gamma = \\sqrt{3}/2$), each point in the "
                "diagram a 2-sphere.",
    "negative": "The expanding universe of Husain, Martinez, and Núñez ($\\gamma = -\\sqrt{3}/2$), each point in the "
                "diagram a 2-sphere.",
    "comoving": "The universe of Husain, Martinez, and Núñez in comoving time ($\\gamma = \\sqrt{3}/2$), "
                "each point in the diagram a 2-sphere.",
    "fonarev_conformal": "Fonarev's expanding universe with an exponential potential ($\\gamma = 4/5$), each point in "
                         "the diagram a 2-sphere.",
    "fonarev": "Fonarev's expanding universe in his time $T$ ($\\gamma = 4/5$), each point in the diagram a 2-sphere.",
}


def curvature_views(ck, src, cv):
    """The Kretschmann scalar K and the Ricci scalar R of the moment eta = eta_0, the embedding
    diagram's, against r from the singularity r = b to 6b, one line view for each massless member,
    read from the conformal chart with eta held.

    What the captions state is checked here: both scalars grow without bound toward r = b; far out
    they tend to the flat Friedmann universe's K = 15 eta_0^2/4 eta^6 and R = -3 eta_0/2 eta^3, which
    at eta = eta_0 are 15/4 eta_0^4 and -3/2 eta_0^2, read in sixty digits at r = 10^8 b and at b = 0;
    R changes sign once, at the radius each caption states; at gamma = sqrt(3)/2 K falls at every r
    drawn, and at gamma = -sqrt(3)/2 it is least at the radius its caption states, K = 0.153/b^4 there."""
    x_axis = cv.axis("$r/b$", 1, 6)
    K_scale = {"kind": "log", "from": 0.1, "to": 1e4}
    R_scale = {"kind": "symlog", "from": -10, "to": 1e3, "linear": 0.1}
    panels = [{"id": "K", "name": "$K\\,b^4$", "aspect": 2, "scale": K_scale, "ticks": cv.value_ticks(K_scale)},
              {"id": "R", "name": "$R\\,b^2$", "aspect": 2, "scale": R_scale, "ticks": cv.value_ticks(R_scale)}]
    cv.shape_checks(ck, METRIC, [p["aspect"] for p in panels])
    height = {p["id"]: cv.PLOT_W / p["aspect"] for p in panels}
    names = (["r", "eta"], ["b", "eta_0", "gamma"])
    K = cv.Held(cv.Published(src, METRIC, "conformal_time", "kretschmann", *names), {"eta": ETA0})
    R = cv.Held(cv.Published(src, METRIC, "conformal_time", "ricci_scalar", *names), {"eta": ETA0})
    within = (1 + 1e-4, 6.0)
    rs = np.linspace(1 + 1e-4, 6, cv.SAMPLES)
    views = []
    for vid in ("positive", "negative"):
        gamma = SIGNS[vid]
        params = {"b": 1.0, "eta_0": ETA0, "gamma": gamma}
        values = list(params.values())
        name = f"{METRIC} at gamma = {'' if gamma > 0 else '-'}sqrt(3)/2"
        curves = (cv.line_curves(ck, name, "K", "K", K, "r", x_axis, K_scale, height["K"], params, within)
                  + cv.line_curves(ck, name, "R", "R", R, "r", x_axis, R_scale, height["R"], params, within))
        ck(f"{name}: K and R grow without bound toward r = b",
           float(K.exact([1 + 1e-12] + values)) > 1e12 and float(R.exact([1 + 1e-12] + values)) > 1e6)
        far_K, far_R = 15 / (4 * ETA0 ** 4), -3 / (2 * ETA0 ** 2)
        for at, b in ((1e8, 1.0), (3.0, 0.0)):
            ck(f"{name}: K tends to the Friedmann universe's 15/4 eta_0^4 (r = {at:g}, b = {b:g})",
               abs(float(K.exact([at, b, ETA0, gamma])) / far_K - 1) < 1e-6)
            ck(f"{name}: R tends to the Friedmann universe's -3/2 eta_0^2 (r = {at:g}, b = {b:g})",
               abs(float(R.exact([at, b, ETA0, gamma])) / far_R - 1) < 1e-6)
        sign = cv.zeros(ck, f"{name} R", R, *within, values)
        want = 1.014 if gamma > 0 else 1.72
        ck(f"{name}: R changes sign once, at r = {want} b", len(sign) == 1 and round(sign[0], 3 if gamma > 0 else 2) == want,
           f"{sign}")
        if gamma > 0:
            ck(f"{name}: K falls at every r drawn", bool(np.all(np.diff(K(rs, *values)) < 0)))
            text = "$K$ diverges at $r = b$ and falls toward the Friedmann universe's $15/4\\eta_0^4$."
        else:
            fine = np.linspace(2, 6, 400001)
            least = int(np.argmin(K(fine, *values)))
            r_least, K_least = float(fine[least]), float(K(fine[least], *values))
            ck(f"{name}: K is least, 0.153/b^4, at r = 4.01 b", round(r_least, 2) == 4.01 and round(K_least, 3) == 0.153,
               f"r = {r_least:.4f}, K = {K_least:.5f}")
            text = ("$K$ diverges at $r = b$, is least at $r \\approx 4.01\\,b$ and tends to the Friedmann universe's "
                    "$15/4\\eta_0^4$.")
        for label, fn, scale in (("K", K, K_scale), ("R", R, R_scale)):
            up, down = cv.beyond(fn, x_axis, scale, values, within)
            ck(f"{name}: {label} leaves its panel only beside r = b, upward",
               len(up) == 1 and up[0][0] < 1.001 and up[0][1] < 1.5 and not down, f"{up} {down}")
        views.append({
            "id": vid, "label": f"$\\gamma = {GAMMA_TEXT[vid]}$", "kind": "line", "x": x_axis, "panels": panels,
            "curves": curves, "marks": [], "legend": [["line", "K", "$K$"], ["line", "R", "$R$"]],
            "caption": [f"The Kretschmann scalar $K$ and the Ricci scalar $R$ of the moment $\\eta = \\eta_0$ against "
                        f"$r$ ($\\gamma = {GAMMA_TEXT[vid]}$). " + text + " $R$ diverges at $r = b$, changes sign at "
                        f"$r \\approx {want}\\,b$ and tends to $-3/2\\eta_0^2$."],
            "ricci_sentences": [[0, 2]],
            "per_scalar": cv.scalar_clips("A logarithmic scale from $0.1$ to $10^4$, cut beside $r = b$.",
                                          "A signed logarithmic scale from $-10$ to $10^3$, linear within $0.1$, cut "
                                          "beside $r = b$."),
            "settings": "$b = 1$, $\\eta_0 = 2b$.",
            "clip": "Logarithmic scales, signed for $R$.",
        })
    return views
