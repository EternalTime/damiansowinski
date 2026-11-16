"""Drawings of dust collapsing in anti-de Sitter space.

The disc of three dimensions, in units of ell: k = 1/4, r_0 = 1, so that M = (k + 1) r_0^2 - 1 = 1/4, the
horizon is R = 1/2 and the energy of the edge is E = sqrt(1 - k r_0^2) = sqrt(3)/2. The scale factor is
a = cos t, the edge R_s = cos t, which crosses the horizon at t = pi/3 and reaches zero at t = pi/2. The
disc's conformal time is eta = artanh(sin t), which has no end, and its conformal radius chi, with
r = 2 sin(chi/2), runs to chi_0 = pi/3. Outside, with kappa = sqrt(M) = 1/2 and
R_* = ln((R - 1/2)/(R + 1/2)), which vanishes as R -> infinity, the edge has

    tan t = sqrt(3) tanh(T/2)                         in the static time, T = 0 at rest,
    tan(t/2) = sqrt(3) (e^(v/2) - 1)/(3 - e^(v/2))    in the advanced time, v = 0 at rest,

and reaches R = 0 at v = ln 3.

The hyperbolic cloud, in units of L: b = 8/7 and r_0 = 2, so that mu = (1 - 1/b^2) r_0^3 = 15/8, the
horizon is the root R_h = 3/2 of R^3 - R - mu, and chi_0 = (8/7) arccosh(7/4). The flat cloud:
b = sqrt(8/7) and R_0 = 2, so that mu = (1 - 1/b^2) R_0^3 = 1, the horizon is R = 1, and
chi_0 = b ln(R_0/b). Their scale factors are solved from each chart's own G^i_i = 3/L^2, which is dust with
the cosmological constant -3/L^2, from rest at a = 1.
"""
import math

import numpy as np
from scipy.integrate import quad

METRIC = "ads_dust_collapse"
DISC, DISC_OUT, DISC_IN = "disc_interior", "disc_exterior", "disc_exterior_ingoing"
HYP, SLAB, HYP_OUT = "hyperbolic_interior", "hyperbolic_slab", "hyperbolic_exterior"
FLAT, FLAT_OUT = "flat_interior", "flat_exterior"

# The disc, ell = 1.
K, R0, MASS = 0.25, 1.0, 0.25
KAPPA = math.sqrt(MASS)
ENERGY = math.sqrt(1 - K * R0 * R0)
T_HORIZON, T_CRUNCH = math.pi / 3, math.pi / 2
CHI_DISC = 2 * math.asin(0.5)
DISC_P = {"k": "1/4", "r_0": 1, "ell": 1}
DISC_OUT_P = {"ell": 1, "M": "1/4"}
DISC_SETTINGS = "$k = 1/(4\\ell^2)$ and $r_0 = \\ell$, so that $M = 1/4$."

# The clouds of four dimensions, L = 1.
HYP_B, HYP_R0, HYP_MU, HYP_RH = 8 / 7, 2.0, 15 / 8, 1.5
HYP_CHI = HYP_B * math.acosh(HYP_R0 / HYP_B)
FLAT_B, FLAT_R0, FLAT_MU, FLAT_RH = math.sqrt(8 / 7), 2.0, 1.0, 1.0
FLAT_CHI = FLAT_B * math.log(FLAT_R0 / FLAT_B)
HYP_P = {"b": "8/7", "r_0": 2, "L": 1}
SLAB_P = {"b": "8/7", "chi_0": "8*acosh(Rational(7, 4))/7", "L": 1}
HYP_OUT_P = {"mu": "15/8", "L": 1}
FLAT_P = {"b": "sqrt(Rational(8, 7))", "chi_0": "sqrt(Rational(8, 7))*log(sqrt(Rational(7, 2)))", "L": 1}
FLAT_OUT_P = {"mu": 1, "L": 1, "b": "sqrt(Rational(8, 7))"}
PLANES = {"theta": "1", "phi": "0"}
FLATS = {"x": "0", "y": "0"}


def friedmann(a, b):
    """The square of da/d(ct) of a cloud of four dimensions released from rest at a = 1, L = 1."""
    return 1 / b ** 2 + (1 - 1 / b ** 2) / a - a ** 2


def time_to(a, b):
    """The proper time from the release to the scale factor a, with a = 1 - s^2 taking the root out of the
    integrand at the release."""
    return quad(lambda s: 2 * s / math.sqrt(friedmann(1 - s * s, b)), 0, math.sqrt(1 - a), limit=200)[0]


HYP_CRUNCH, FLAT_CRUNCH = time_to(0.0, HYP_B), time_to(0.0, FLAT_B)
HYP_CROSS, FLAT_CROSS = time_to(HYP_RH / HYP_R0, HYP_B), time_to(FLAT_RH / FLAT_R0, FLAT_B)


def face_at_static_time(b, r0, mu, hyperbolic):
    """The areal radius of a cloud's face at the exterior's static time T, T = 0 at rest, L = 1. With
    a = 1 - s^2 the face obeys ds/dT = sqrt(F(a)) f(r_0 a)/(2 s E), F the Friedmann function, f the exterior's
    metric function and E the energy of the face, which is regular at the release, where sqrt(F)/s tends to
    sqrt(-F'(1)), and stalls on the horizon, where f vanishes."""
    from scipy.integrate import solve_ivp
    energy = math.sqrt(r0 * r0 / (b * b) - 1) if hyperbolic else r0 / b
    slope = math.sqrt((1 - 1 / b ** 2) + 2)        # sqrt(-F'(1)), F'(1) = -(1 - 1/b^2) - 2

    def f(R):
        return (-1.0 if hyperbolic else 0.0) - mu / R + R * R

    def rhs(_, y):
        s = y[0]
        root = slope if s < 1e-7 else math.sqrt(max(friedmann(1 - s * s, b), 0.0)) / s
        return [root * f(r0 * (1 - s * s)) / (2 * energy)]
    sol = solve_ivp(rhs, (0.0, 40.0), [0.0], rtol=1e-12, atol=1e-14, dense_output=True, method="DOP853",
                    max_step=0.01)

    def radius(T):
        T = np.asarray(T, dtype=float)
        s = sol.sol(np.clip(T.ravel(), 0.0, 40.0))[0].reshape(T.shape)
        return np.where(T >= 0, r0 * (1 - s * s), np.nan)
    return radius


def rstar(R):
    """The disc's exterior tortoise coordinate, zero as R -> infinity."""
    R = np.asarray(R, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.log(np.abs((R - KAPPA) / (R + KAPPA))) / (2 * KAPPA)


# The edge of the disc in the advanced time, v = 0 at rest: R_s = (1 - s^2)/(1 + s^2), s = tan(t/2).
_HALF = "(sqrt(3)*(exp(v/2) - 1)/(3 - exp(v/2)))"
EDGE_IN = f"(1 - {_HALF}**2)/(1 + {_HALF}**2)"
# What lies outside the edge and after the moment of rest, T = v - R_*(R) + R_*(r_0) >= 0 outside the horizon.
AFTER_REST = (f"Min(R - {EDGE_IN}, Piecewise((v - log((2*R - 1)/(2*R + 1)) - log(3), R > Rational(1, 2)), "
              "(1, True)))")
DUST = {"funcs": ["a"], "rates": [0.0], "start": [1.0], "origin": "reference", "sources": ["3"]}


# -- Spacetime diagrams ----------------------------------------------------------------------------

def register_null(nr):
    import sympy as sp

    class adc_face_hyperbolic(sp.Function):
        """The areal radius of the hyperbolic cloud's face at the static time T."""
        nargs = 1
        is_real = True
        _imp_ = staticmethod(face_at_static_time(HYP_B, HYP_R0, HYP_MU, True))

    class adc_face_flat(sp.Function):
        """The areal radius of the flat cloud's face at the static time T."""
        nargs = 1
        is_real = True
        _imp_ = staticmethod(face_at_static_time(FLAT_B, FLAT_R0, FLAT_MU, False))
    nr.DECLARED_FUNCTIONS.update({"adc_face_hyperbolic": adc_face_hyperbolic, "adc_face_flat": adc_face_flat})
    fallen = ("The face of the cloud of $b = {b}$ and ${r} = 2L$, at rest at $R = 2L$ at $T = 0$, on the radial "
              "geodesic it follows, with its areal radius ${r}a$ solved from the cloud's own Friedmann equation.")
    edge = "the edge of the disc, falling freely from rest"
    face = "the face of the cloud, falling freely from rest"
    released = ("The {what} released from rest at $R = {at}$ at $T = 0$, the {which} of the comoving "
                "coordinates.")
    solved = ("Dust released from rest at $t = 0$ with $a = 1$, $a(t)$ solved from this spacetime's own "
              "$G^{{{c}}}{{}}_{{{c}}} = 3/L^2$, with $b = {b}$ and {rest}.")
    rows = (
        nr.Diagram(METRIC, DISC, "through", "through the centre", ("t", "r"), (0, 1, 0, T_CRUNCH), "$r/\\ell$",
                   "$ct/\\ell$", DISC_P, {"theta": "0"}, mirror=True, families=nr.SIDEWAYS, cones=(4, 8),
                   areal=True, functions={"a": "cos(t)"},
                   lines=(("surface", "r", "1", "the edge of the disc, $r = r_0$"),),
                   marked=(("event", {"x0": "pi/3", "r": "1"}, 1, "the event horizon"),),
                   input="The scale factor $a = \\cos(ct/\\ell)$ of a disc released from rest, with "
                         + DISC_SETTINGS),
        nr.Diagram(METRIC, DISC_OUT, "radial", "$T$ and $R$", ("T", "R"), (0, 3, 0, 4), "$R/\\ell$", "$cT/\\ell$",
                   DISC_OUT_P, {"theta": "0"}, orient="ingoing", surface="1", surface_legend=edge,
                   input=released.format(what="edge", at="\\ell", which="disc of $k = 1/(4\\ell^2)$ and $r_0 = \\ell$")),
        nr.Diagram(METRIC, DISC_IN, "finkelstein", "against $v - R$", ("v", "R"), (0, 2.5, -2, 1.5), "$R/\\ell$",
                   "$(v - R)/\\ell$", DISC_OUT_P, {"theta": "0"}, to_display=nr.FINKELSTEIN_IN, orient="ingoing",
                   where=AFTER_REST, curves=(("shell", EDGE_IN, "the edge of the disc"),),
                   input="The edge of the disc of $k = 1/(4\\ell^2)$ and $r_0 = \\ell$, "
                         "$R_s = \\ell\\,(1 - s^2)/(1 + s^2)$ with "
                         "$s = \\sqrt{3}\\,(e^{v/2\\ell} - 1)/(3 - e^{v/2\\ell})$ and $v = 0$ at the moment of rest."),
        nr.Diagram(METRIC, HYP, "radial", "$t$ and $r$", ("t", "r"), (1, 2, 0, HYP_CRUNCH), "$r/L$", "$ct/L$",
                   HYP_P, PLANES, areal=True, dust={**DUST, "eqs": [["r", "r"]], "params": {"b": 8 / 7}},
                   lines=(("surface", "r", "2", "the face of the cloud, $r = r_0$"),),
                   input=solved.format(c="r", b="8L/7", rest="$r_0 = 2L$")),
        nr.Diagram(METRIC, SLAB, "through", "through the smallest surface", ("t", "\\chi"),
                   (-HYP_CHI, HYP_CHI, 0, HYP_CRUNCH), "$\\chi/L$", "$ct/L$", SLAB_P, PLANES, families=nr.SIDEWAYS,
                   areal=True, dust={**DUST, "eqs": [["\\chi", "\\chi"]], "params": {"b": 8 / 7}},
                   rounded={"chi_0": "1.324"},
                   input=solved.format(c="\\chi", b="8L/7", rest="$r_0 = 2L$ at the faces $\\chi = \\pm\\chi_0$")),
        nr.Diagram(METRIC, HYP_OUT, "radial", "$T$ and $R$", ("T", "R"), (0, 4, 0, 4), "$R/L$", "$cT/L$",
                   HYP_OUT_P, PLANES, orient="ingoing", areal=True, where="R - adc_face_hyperbolic(T)",
                   curves=(("shell", "adc_face_hyperbolic(T)", face),), input=fallen.format(b="8L/7", r="r_0")),
        nr.Diagram(METRIC, FLAT, "radial", "$t$ and $\\chi$", ("t", "\\chi"), (FLAT_CHI - 2.4, FLAT_CHI, 0, FLAT_CRUNCH),
                   "$\\chi/L$", "$ct/L$", FLAT_P, FLATS, families=nr.SIDEWAYS,
                   dust={**DUST, "eqs": [["\\chi", "\\chi"]], "params": {"b": FLAT_B}},
                   rounded={"b": "1.069", "chi_0": "0.670"},
                   input=solved.format(c="\\chi", b="\\sqrt{8/7}\\,L", rest="$R_0 = 2L$ at the face $\\chi = \\chi_0$")),
        nr.Diagram(METRIC, FLAT_OUT, "radial", "$T$ and $R$", ("T", "R"), (0, 4, 0, 4), "$R/L$", "$cT/L$",
                   FLAT_OUT_P, FLATS, orient="ingoing", where="R - adc_face_flat(T)",
                   curves=(("shell", "adc_face_flat(T)", face),), rounded={"b": "1.069"},
                   input=fallen.format(b="\\sqrt{8/7}\\,L", r="R_0")),
    )
    nr.DIAGRAMS.extend(rows)
    nr.CAPTIONS[(METRIC, DISC, "through")] = [
        "The line through the centre of the disc ($\\theta = 0$ on the right, $\\theta = \\pi$ on the left), from "
        "the release to the crunch at $ct = \\pi\\ell/2$. The event horizon leaves the centre at "
        "$ct = 0.27\\,\\ell$ and meets the edge at $ct = \\pi\\ell/3$."]
    nr.CAPTIONS[(METRIC, DISC_OUT, "radial")] = [
        "The plane of $T$ and $R$ ($\\theta = 0$) outside the disc. The edge falls from rest at $R = \\ell$ and "
        "nears the horizon $R = \\ell/2$ only as $T \\to \\infty$, and the outgoing rays reach $R \\to \\infty$ in a "
        "finite $T$."]
    nr.CAPTIONS[(METRIC, DISC_IN, "finkelstein")] = [
        "The plane of $v$ and $R$ ($\\theta = 0$) outside the disc after its moment of rest, drawn against "
        "$v - R$. The edge crosses the horizon $R = \\ell/2$ at $v = 0.72\\,\\ell$ and reaches $R = 0$ at "
        "$v = \\ell\\ln 3$, and the ingoing rays are the lines of constant $v$."]
    nr.CAPTIONS[(METRIC, HYP, "radial")] = [
        "The plane of $t$ and $r$ ($\\theta = 1$, $\\phi = 0$) from the smallest surface $r = b$ to a face "
        f"$r = r_0$ of the hyperbolic cloud, from the release to the crunch at $ct = {HYP_CRUNCH:.2f}\\,L$. The "
        f"face crosses the horizon $R = 3L/2$ at $ct = {HYP_CROSS:.2f}\\,L$."]
    nr.CAPTIONS[(METRIC, SLAB, "through")] = [
        "The plane of $t$ and $\\chi$ ($\\theta = 1$, $\\phi = 0$) across the hyperbolic cloud, from one face to "
        f"the other through the smallest surface $\\chi = 0$, up to the crunch at $ct = {HYP_CRUNCH:.2f}\\,L$. Each "
        "face meets a black hole of its own."]
    nr.CAPTIONS[(METRIC, HYP_OUT, "radial")] = [
        "The plane of $T$ and $R$ ($\\theta = 1$, $\\phi = 0$) outside a face of the hyperbolic cloud. The face "
        "falls from rest at $R = 2L$ and nears the horizon $R = 3L/2$ only as $T \\to \\infty$."]
    nr.CAPTIONS[(METRIC, FLAT, "radial")] = [
        "The plane of $t$ and $\\chi$ ($x = 0$, $y = 0$) in the flat cloud, from the face $\\chi = \\chi_0$ toward "
        f"the cusp at $\\chi \\to -\\infty$, up to the crunch at $ct = {FLAT_CRUNCH:.2f}\\,L$. The face crosses "
        f"the horizon $R = L$ at $ct = {FLAT_CROSS:.2f}\\,L$."]
    nr.CAPTIONS[(METRIC, FLAT_OUT, "radial")] = [
        "The plane of $T$ and $R$ ($x = 0$, $y = 0$) outside the face of the flat cloud. The face falls from rest "
        "at $R = 2L$ and nears the horizon $R = L$ only as $T \\to \\infty$."]

    def own(f):
        return lambda w, r: f(np.asarray(w, dtype=float), np.asarray(r, dtype=float))

    def conformal_time(t):
        return np.arctanh(np.sin(t))

    def conformal_radius(r):
        return 2 * np.arcsin(np.clip(r / 2, -1, 1))

    # The disc: light keeps eta -+ chi, with eta = artanh(sin t) and r = 2 sin(chi/2).
    nr.CLOSED_FORMS[(METRIC, DISC, "through")] = (
        own(lambda t, r: conformal_time(t) + conformal_radius(r)),
        own(lambda t, r: conformal_time(t) - conformal_radius(r)),
        own(lambda t, r: (t < T_CRUNCH - 0.05) & (np.abs(r) < 0.98)))
    nr.CLOSED_FORMS[(METRIC, DISC_OUT, "radial")] = (
        own(lambda T, R: T + rstar(R)), own(lambda T, R: T - rstar(R)), own(lambda T, R: R > KAPPA + 0.05))
    nr.CLOSED_FORMS[(METRIC, DISC_IN, "finkelstein")] = (
        own(lambda v, R: v), own(lambda v, R: v - 2 * rstar(R)), own(lambda v, R: np.abs(R - KAPPA) > 0.05))


# -- Embedding diagrams ---------------------------------------------------------------------------------

LEVEL = math.sqrt(1 + MASS)                 # where the disc's exterior has g_RR = 1, sqrt(5)/2
DISC_TOP = 2.5
# The flat cloud with x closed into a circle of length 2 pi HORN: the circle at chi has the radius
# HORN e^(chi/b) and outside the radius HORN R/b. HORN is chosen so that the exterior lies level, with
# (d rho/dR)^2 = g_RR, at R = 4: f(4) = 63/4 = b^2/HORN^2, so HORN = 4 sqrt(2)/21, and the face's circle has
# the radius 4 sqrt(7)/21, under half of b.
FLAT_LEVEL = 4
HORN = 4 * math.sqrt(2) / 21
HORN_TEXT = "4*sqrt(2)/21"
HORN_LENGTH = 3.0                           # how far in from the face the horn is drawn, in b
FLAT_TOP = 5.5
assert abs(FLAT_LEVEL ** 2 - FLAT_MU / FLAT_LEVEL - (FLAT_B / HORN) ** 2) < 1e-13


def register_embedding(eb):
    eb.DRAWN[METRIC] = lambda ck, src: embedding(eb, ck, src)
    eb.CAPTIONS[(METRIC, "rest")] = [
        "The disc at its moment of rest ($k = 1/(4\\ell^2)$, $r_0 = \\ell$). The dust is a cap of a sphere of "
        "radius $2\\ell$. The exterior lies in flat space out to $R = \\sqrt{5}\\,\\ell/2$ and in Minkowski "
        "space ($dX^2 + dY^2 - dZ^2$) beyond."]
    eb.CAPTIONS[(METRIC, "horn")] = [
        "The surface $y = 0$ of the flat cloud at rest, with $x$ closed into a circle. The dust is a horn narrowing "
        "toward its cusp at $\\chi \\to -\\infty$. The exterior lies in flat space out to $R = 4L$ and in Minkowski "
        "space beyond."]


def embedding(eb, ck, src):
    """The disc at rest: the published interior at t = 0, a = 1, the cap of a sphere of radius 1/sqrt(k) = 2
    out to r_0 = 1, and the published exterior at T = 0 from R = 1, in flat space to R = sqrt(1 + M), where
    g_RR = 1 and the surface lies level, and in Minkowski space beyond. The two meet with one tangent, since
    g_rr = 1/(1 - k r_0^2) inside and 1/(r_0^2 - M) outside are one number. The flat cloud at rest, with x a
    circle of length 2 pi HORN: inside, the circles have the radius HORN e^(chi/b), a pseudosphere of radius b,
    and outside HORN R/b, in flat space to R = 4, where f = b^2/HORN^2, and in Minkowski space beyond."""
    inner = eb.Slice(src, METRIC, DISC, "r", "\\theta", {"t": 0}, DISC_P, {"a": "cos(t)"})
    outer = eb.Slice(src, METRIC, DISC_OUT, "R", "\\theta", {"T": 0}, DISC_OUT_P)
    outer_m = eb.Slice(src, METRIC, DISC_OUT, "R", "\\theta", {"T": 0}, DISC_OUT_P, space="minkowski")
    size = 2 * DISC_TOP
    ck.add("Anti-de Sitter dust collapse: the disc's exterior lies level where R^2 - M = 1",
           abs(float(outer.defect_at(np.array([LEVEL]))[0])), 1e-12)
    ck.stops("Anti-de Sitter dust collapse, the disc's exterior beyond sqrt(1 + M) in flat space", outer,
             np.linspace(LEVEL, 40, 402)[1:])
    join = ("at $R = \\sqrt{1 + M}\\,\\ell$ the surface lies level, in flat space nearer the disc and in Minkowski "
            "space beyond")
    dust = eb.Piece("dust", "star", inner, 0.0, R0, 0.0, 1,
                    (("axis", "the centre $r = 0$, where the cap is smooth"), ("join", "the edge of the disc, $r = r_0$")),
                    [(R0 / 3, "r", None), (2 * R0 / 3, "r", None), (R0, "surface", "$r = r_0$")], size)
    near = eb.Piece("exterior", "sheet", outer, R0, LEVEL, dust.at(R0)[1], 1,
                    (("join", "the edge of the disc, $R = r_0$"), ("join", join)), [(LEVEL, "space", None)], size)
    far = eb.Piece("exterior_minkowski", "sheet", outer_m, LEVEL, DISC_TOP, near.at(LEVEL)[1], 1,
                   (("join", join), ("edge", "the sheet runs on toward a light cone, to $R \\to \\infty$")),
                   [(1.5, "r", None), (2.0, "r", None), (DISC_TOP, "r", None)], size)
    for p in (dust, near, far):
        space = "in Minkowski space" if p.sl.lorentz else "in flat space"
        ck.isometry(f"Anti-de Sitter dust collapse, the disc, {p.id} {space}", p)
        ck.radius(f"Anti-de Sitter dust collapse, the disc, {p.id}, rho = R {space}", p, lambda r: r, size)
    ck.form("Anti-de Sitter dust collapse: the dust is a cap of a sphere of radius 1/sqrt(k)", dust,
            lambda r: (1 - np.sqrt(1 - K * r * r)) / math.sqrt(K), size)
    ck.join("Anti-de Sitter dust collapse: the dust meets its exterior with one tangent", dust, R0, near, R0)
    ck.join("Anti-de Sitter dust collapse: the disc's exterior in flat and in Minkowski space", near, LEVEL, far, LEVEL)
    disc = eb.Surface([dust, near, far])
    fig = eb.figure_of([disc], {"star": "star", "sheet": "cover"}, size)
    eb.ring_label(fig, [0, 0, 0], *dust.at(R0), "$r = r_0$", dx=10)
    eb.ring_label(fig, [0, 0, 0], *far.at(DISC_TOP), "$R = 2.5\\,\\ell$")
    fig.legend("fill", "star", "the dust")
    fig.legend("fill", "cover", "the exterior")
    fig.legend("line", "r", "$r$ and $R$ constant")
    fig.legend("line", "surface", "the edge $r = r_0$")
    fig.legend("line", "space", "$R = \\sqrt{5}\\,\\ell/2$")
    fig.legend("line", "meridian", "$\\theta$ constant")
    views = [eb.view("rest", "The disc at rest", "$\\ell$", [disc], fig.done(),
                     settings="$\\ell = 1$, $M = 1/4$, $t = 0$, and $T = 0$.")]

    level = FLAT_LEVEL
    swept = {"x": f"{HORN_TEXT}*x"}
    cloud = eb.Slice(src, METRIC, FLAT, "\\chi", "x", {"t": 0, "y": 0}, FLAT_P, {"a": "1"}, swept=swept)
    out = eb.Slice(src, METRIC, FLAT_OUT, "R", "x", {"T": 0, "y": 0}, FLAT_OUT_P, swept=swept)
    out_m = eb.Slice(src, METRIC, FLAT_OUT, "R", "x", {"T": 0, "y": 0}, FLAT_OUT_P, swept=swept, space="minkowski")
    tip = FLAT_CHI - HORN_LENGTH * FLAT_B
    size = 2 * HORN * FLAT_TOP / FLAT_B + 1.0
    ck.add("Anti-de Sitter dust collapse: the flat cloud's exterior lies level where f = b^2/HORN^2",
           abs(float(out.defect_at(np.array([float(level)]))[0])), 1e-12)
    ck.stops("Anti-de Sitter dust collapse, the flat cloud's exterior beyond that radius in flat space", out,
             np.linspace(level, 40, 402)[1:])
    join = "the surface lies level here, in flat space nearer the cloud and in Minkowski space beyond"
    marks = [(FLAT_CHI - k * FLAT_B, "r", None) for k in (2, 1)] + [(FLAT_CHI, "surface", "$\\chi = \\chi_0$")]
    horn = eb.Piece("dust", "star", cloud, tip, FLAT_CHI, 0.0, 1,
                    (("edge", "the horn runs on, ever narrower, to its cusp at $\\chi \\to -\\infty$"),
                     ("join", "the face of the cloud, $\\chi = \\chi_0$")), marks, size)
    near = eb.Piece("exterior", "sheet", out, FLAT_R0, level, horn.at(FLAT_CHI)[1], 1,
                    (("join", "the face of the cloud, $R = R_0$"), ("join", join)), [(3.0, "r", None), (level, "space", None)],
                    size)
    far = eb.Piece("exterior_minkowski", "sheet", out_m, level, FLAT_TOP, near.at(level)[1], 1,
                   (("join", join), ("edge", "the sheet runs on toward a light cone, to $R \\to \\infty$")),
                   [(FLAT_TOP, "r", None)], size)
    for p in (horn, near, far):
        space = "in Minkowski space" if p.sl.lorentz else "in flat space"
        ck.isometry(f"Anti-de Sitter dust collapse, the flat cloud, {p.id} {space}", p)
    ck.radius("Anti-de Sitter dust collapse: the horn's circles have the radius HORN e^(chi/b)", horn,
              lambda chi: HORN * np.exp(chi / FLAT_B), size)
    for p in (near, far):
        ck.radius(f"Anti-de Sitter dust collapse: outside the flat cloud, {p.id}, rho = HORN R/b", p,
                  lambda R: HORN * R / FLAT_B, size)
    ck.join("Anti-de Sitter dust collapse: the horn meets its exterior with one tangent", horn, FLAT_CHI, near, FLAT_R0)
    ck.join("Anti-de Sitter dust collapse: the flat cloud's exterior in flat and in Minkowski space", near, level, far, level)
    # A pseudosphere of radius b: along the horn the circle's radius grows at d rho/d chi = rho/b.
    chi = np.linspace(tip, FLAT_CHI, 50)
    ck.add("Anti-de Sitter dust collapse: the horn is a pseudosphere of radius b, d rho/d chi = rho/b",
           float(np.max(np.abs(np.gradient(cloud.rho_at(chi), chi)[1:-1] * FLAT_B / cloud.rho_at(chi)[1:-1] - 1))), 1e-3)
    flat = eb.Surface([horn, near, far])
    fig = eb.figure_of([flat], {"star": "star", "sheet": "cover"}, size)
    eb.ring_label(fig, [0, 0, 0], *horn.at(FLAT_CHI), "$\\chi = \\chi_0$", dx=10)
    eb.ring_label(fig, [0, 0, 0], *far.at(FLAT_TOP), "$R = 5.5\\,L$")
    fig.legend("fill", "star", "the dust")
    fig.legend("fill", "cover", "the exterior")
    fig.legend("line", "r", "$\\chi$ and $R$ constant")
    fig.legend("line", "surface", "the face $\\chi = \\chi_0$")
    fig.legend("line", "space", "$R = 4L$")
    fig.legend("line", "meridian", "$x$ constant")
    views.append(eb.view("horn", "The flat cloud at rest", "$L$", [flat], fig.done(),
                         settings="$b = \\sqrt{8/7}\\,L$, $R_0 = 2L$, and $x$ of period "
                                  "$8\\sqrt{2}\\,\\pi L/21$."))
    return views


# -- The moments on the spacetime diagrams ----------------------------------------------------------------

def register_slices(sl):
    def rest():
        return sl.moments(METRIC, "rest")[0]

    def horn():
        return sl.moments(METRIC, "horn")[0]

    def disc():
        m = rest()
        return [sl.Mark(m, sl.across(0.0, *m.reach(DISC, "r")))]

    def disc_static():
        m = rest()
        return [sl.Mark(m, sl.along(0.0, *m.reach(DISC_OUT, "R")))]

    def disc_ingoing():
        # T = 0 is v = R_*(R) - R_*(r_0), with v = 0 on the edge.
        m = rest()
        R = np.linspace(*m.reach(DISC_OUT, "R"), 401)
        return [sl.Mark(m, [np.column_stack([rstar(R) - rstar(R0), R])])]

    def cloud():
        m = horn()
        return [sl.Mark(m, sl.along(0.0, *m.reach(FLAT, "\\chi")))]

    def cloud_static():
        m = horn()
        return [sl.Mark(m, sl.along(0.0, *m.reach(FLAT_OUT, "R")))]

    sl.FLAT[(METRIC, DISC, "through")] = disc
    sl.FLAT[(METRIC, DISC_OUT, "radial")] = disc_static
    sl.FLAT[(METRIC, DISC_IN, "finkelstein")] = disc_ingoing
    sl.FLAT[(METRIC, FLAT, "radial")] = cloud
    sl.FLAT[(METRIC, FLAT_OUT, "radial")] = cloud_static
    for key in ((METRIC, HYP, "radial"), (METRIC, SLAB, "through"), (METRIC, HYP_OUT, "radial")):
        sl.FLAT[key] = lambda: []
    sl.FLAT_METRICS.add(METRIC)
