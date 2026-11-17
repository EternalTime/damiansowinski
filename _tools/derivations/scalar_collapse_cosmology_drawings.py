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
