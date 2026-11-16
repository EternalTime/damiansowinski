"""Drawings of dust collapsing in anti-de Sitter space.

The disc of three dimensions, in units of ell: k = 1/4, r_0 = 1, so that M = (k + 1) r_0^2 - 1 = 1/4, the
horizon is R = 1/2 and the energy of the edge is E = sqrt(1 - k r_0^2) = sqrt(3)/2. The scale factor is
a = cos t, the edge R_s = cos t, which crosses the horizon at t = pi/3 and reaches zero at t = pi/2. The
disc's conformal time is eta = artanh(sin t), which has no end, and its conformal radius chi, with
r = 2 sin(chi/2), runs to chi_0 = pi/3. Outside, with kappa = sqrt(M) = 1/2 and
R_* = ln((R - 1/2)/(R + 1/2)), which vanishes as R -> infinity, the edge has

    tan t = sqrt(3) tanh(T/2)                         in the static time, T = 0 at rest,
    tan(t/2) = sqrt(3) (e^(v/2) - 1)/(3 - e^(v/2))    in the advanced time, v = 0 at rest,

crosses the horizon at v = ln(9/4) and reaches R = 0 at v = ln 3.

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
        "nears the horizon $R = \\ell/2$ only as $T \\to \\infty$. The outgoing rays reach $R \\to \\infty$ in a "
        "finite $T$."]
    nr.CAPTIONS[(METRIC, DISC_IN, "finkelstein")] = [
        "The plane of $v$ and $R$ ($\\theta = 0$) outside the disc after its moment of rest, drawn against "
        "$v - R$. The edge crosses the horizon $R = \\ell/2$ at $v = \\ell\\ln(9/4)$ and reaches $R = 0$ at "
        "$v = \\ell\\ln 3$. The ingoing rays are the lines of constant $v$."]
    nr.CAPTIONS[(METRIC, HYP, "radial")] = [
        "The plane of $t$ and $r$ ($\\theta = 1$, $\\phi = 0$) from the smallest surface $r = b$ to a face "
        f"$r = r_0$ of the hyperbolic cloud, from the release to the crunch at $ct = {HYP_CRUNCH:.2f}\\,L$. The "
        f"face crosses the horizon $R = 3L/2$ at $ct = {HYP_CROSS:.2f}\\,L$."]
    nr.CAPTIONS[(METRIC, SLAB, "through")] = [
        "The plane of $t$ and $\\chi$ ($\\theta = 1$, $\\phi = 0$) across the hyperbolic cloud, from one face to "
        f"the other through the smallest surface $\\chi = 0$, up to the crunch at $ct = {HYP_CRUNCH:.2f}\\,L$. Each "
        "face is joined to a black hole of its own."]
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


# -- Conformal diagram ----------------------------------------------------------------------------------
#
# The disc, by the construction of Oppenheimer and Snyder's diagram. The dust is conformal to the strip of
# its conformal time eta = artanh(sin t), which has no end, and chi, 0 <= chi <= chi_0, and is drawn with
# p = arctan(eta - chi), q = arctan(eta + chi): its centre is the line X = 0, its moment of rest the line
# T = 0 and its crunch the one point p = q = pi/2. Outside, Kruskal's U = -e^(-kappa u), V = e^(kappa v), with
# kappa = sqrt(M) and UV = (kappa - R)/(kappa + R), are drawn as P(U) and Q(V), fixed by three conditions:
#
#   1. on the edge the two sides agree. With s = tan(t/2) = tanh(eta/2) the edge has
#      V_s = C (r_0 - kappa + E s)/(r_0 + kappa + E s) and U_s = -C (r_0 - kappa - E s)/(r_0 + kappa - E s),
#      C = sqrt((r_0 + kappa)/(r_0 - kappa)), so P(U_s) = arctan(eta - chi_0) and Q(V_s) = arctan(eta + chi_0);
#   2. the moment of rest is the line T = 0, P(U) = -Q(-U), for the rays that start outside the disc;
#   3. R = 0, where UV = 1, is the line T = pi, Q(V) = pi - P(1/V), for the rays that arrive after the crunch.
#
# U_s vanishes at s = (r_0 - kappa)/E, which is t = pi/3, the horizon, and U_s V_s = 1 at s = 1, the crunch.

EDGE_C = math.sqrt((R0 + KAPPA) / (R0 - KAPPA))
U_REST, V_REST = -1 / EDGE_C, 1 / EDGE_C


def edge_kruskal(s):
    """Kruskal's U and V on the edge at s = tan(t/2)."""
    s = np.asarray(s, dtype=float)
    return (-EDGE_C * (R0 - KAPPA - ENERGY * s) / (R0 + KAPPA - ENERGY * s),
            EDGE_C * (R0 - KAPPA + ENERGY * s) / (R0 + KAPPA + ENERGY * s))


def _eta(s):
    with np.errstate(divide="ignore"):
        return 2 * np.arctanh(np.clip(s, 0.0, 1.0))


def disc_Q(V):
    """Q(V) outside the disc: on the rays that reach the edge, and by condition 3 after the crunch."""
    V = np.asarray(V, dtype=float)
    w = np.clip(V, V_REST, 1.0) / EDGE_C
    on_edge = np.arctan(_eta((w * (R0 + KAPPA) - (R0 - KAPPA)) / (ENERGY * (1 - w))) + CHI_DISC)
    return np.where(V < 1, on_edge, math.pi - _disc_P_edge(1 / np.maximum(V, 1.0)))


def _disc_P_edge(U):
    w = -np.clip(U, U_REST, 1.0) / EDGE_C
    return np.arctan(_eta((R0 - KAPPA - w * (R0 + KAPPA)) / (ENERGY * (1 - w))) - CHI_DISC)


def disc_P(U):
    """P(U) outside the disc: on the rays that leave the edge, and by condition 2 for those that start
    outside it at the moment of rest."""
    U = np.asarray(U, dtype=float)
    return np.where(U >= U_REST, _disc_P_edge(U), -disc_Q(-np.minimum(U, U_REST)))


def register_conformal(cf):
    cf.DRAWN[METRIC] = lambda ck, src: conformal(cf, ck, src)
    cf.CAPTIONS[(METRIC, "disc")] = [
        "The disc and its exterior from the moment of rest ($M = 1/4$), each point in the diagram a circle. The edge "
        "crosses the event horizon and reaches $R = 0$ at the crunch."]


def conformal(cf, ck, src):
    PI = math.pi
    inner = cf.Plane(src, METRIC, DISC, ("t", "r"), {"theta": "0"}, DISC_P, functions={"a": "cos(t)"})
    ext = cf.Plane(src, METRIC, DISC_OUT, ("T", "R"), {"theta": "0"}, DISC_OUT_P)
    ingoing = cf.Plane(src, METRIC, DISC_IN, ("v", "R"), {"theta": "0"}, DISC_OUT_P)

    def inside(t, r):
        eta, chi = np.arctanh(np.sin(t)), 2 * np.arcsin(np.asarray(r, dtype=float) / 2)
        return np.arctan(eta - chi), np.arctan(eta + chi)

    def outside(T, R):
        root = np.sqrt((R - KAPPA) / (R + KAPPA))
        return disc_P(-np.exp(-KAPPA * T) * root), disc_Q(np.exp(KAPPA * T) * root)

    def advanced(v, R):
        V = np.exp(KAPPA * v) / EDGE_C
        return disc_P((KAPPA - R) / ((KAPPA + R) * V)), disc_Q(V)

    def inner_region(T, R):
        """The exterior inside its horizon, R < kappa, on its own static time T."""
        root = np.sqrt((KAPPA - R) / (KAPPA + R))
        return disc_P(np.exp(-KAPPA * T) * root), disc_Q(np.exp(KAPPA * T) * root)
    name = "Anti-de Sitter dust collapse"
    ck.chart(f"{name}, the disc", inner, inside, ck.uniform(0.01, 1.5), ck.uniform(0.001, R0), lambda t, r: (1, 0))
    ck.chart(f"{name}, the disc's exterior", ext, outside, ck.uniform(0.05, 8), ck.uniform(1.0, 30),
             lambda T, R: (1, 0))
    ck.chart(f"{name}, the disc's exterior on the advanced time", ingoing, advanced, ck.uniform(1.2, 4),
             ck.uniform(0.02, 3), lambda v, R: (1, -1))
    s = np.linspace(0, 1, 2001)
    U, V = edge_kruskal(s)
    radius = R0 * (1 - s * s) / (1 + s * s)
    ck.limit(f"{name}: on the edge UV = (kappa - R)/(kappa + R)", U * V, (KAPPA - radius) / (KAPPA + radius), 1e-12)
    ck.limit(f"{name}: the edge is at rest on T = 0, U = -V", [U[0] + V[0]], [0.0], 1e-12)
    ck.limit(f"{name}: the edge crosses the horizon at t = pi/3", edge_kruskal(math.tan(PI / 6))[0], 0.0, 1e-12)
    ck.limit(f"{name}: the edge reaches R = 0 on UV = 1", [U[-1] * V[-1]], [1.0], 1e-12)
    ck.limit(f"{name}: the two sides agree on the edge",
             np.concatenate([disc_P(U[:-1]) - np.arctan(_eta(s[:-1]) - CHI_DISC),
                             disc_Q(V[:-1]) - np.arctan(_eta(s[:-1]) + CHI_DISC)]), 0.0, 1e-9)
    ck.limit(f"{name}: the moment of rest is one line outside, P(-V) + Q(V) = 0",
             [float(disc_P(-x) + disc_Q(x)) for x in (0.6, 0.8, 0.95)], [0, 0, 0], 1e-12)
    ck.limit(f"{name}: R = 0 is one line, P(1/V) + Q(V) = pi",
             [float(disc_P(1 / x) + disc_Q(x)) for x in (1.5, 10.0, 1e4)], [PI] * 3, 1e-12)
    # Marginally trapped circles: the radius a r has a null gradient where r^2 (sin^2 t + k) = 1.
    t = np.linspace(T_HORIZON, T_CRUNCH - 0.02, 60)
    trapped = 1 / np.sqrt(np.sin(t) ** 2 + K)
    _, _, _, h00, h01, h11 = inner.metric(t, trapped)
    dR = [-np.sin(t) * trapped, np.cos(t)]
    ck.limit(f"{name}: the radius a r has a null gradient on r^2 (sin^2 t + k) = 1",
             (h00 * dR[0] ** 2 + 2 * h01 * dR[0] * dR[1] + h11 * dR[1] ** 2) / (dR[1] ** 2 * np.abs(h11)), 0, 1e-9)
    ck.limit(f"{name}: the trapped circles begin on the edge at the horizon", [trapped[0]], [R0], 1e-12)
    ck.diverges(f"{name}: the Kretschmann scalar diverges at the crunch of the disc",
                inner.kretschmann(T_CRUNCH - 1e-2, 0.3), inner.kretschmann(T_CRUNCH - 1e-3, 0.3))
    ck.finite(f"{name}: outside, the Kretschmann scalar is finite at R = 0",
              ext.kretschmann(np.zeros(3), np.array([1e-3, 1e-2, 0.1])))

    p0 = float(disc_P(0.0))                        # the event horizon, U = 0
    corner = [PI - 2 * p0, PI]                     # where the boundary meets R = 0
    eta = np.concatenate([np.linspace(0, 6, 1200), np.geomspace(6, 1e6, 400)[1:]])
    edge_xt = [[float(np.arctan(e + CHI_DISC) - np.arctan(e - CHI_DISC)),
                float(np.arctan(e + CHI_DISC) + np.arctan(e - CHI_DISC))] for e in eta] + [[0.0, PI]]
    big = np.geomspace(1.0, 1e9, 1600)
    far_p, far_q = disc_P(-1 / big), disc_Q(big)   # the conformal boundary, UV = -1, from T = 0 up
    far_xt = [[float(q - p), float(p + q)] for p, q in zip(far_p, far_q)] + [corner]
    v = cf.View("disc", "The disc and its exterior", [-0.35, PI + 0.35, -0.25, PI + 0.3])
    dust = [[0.0, 0.0]] + edge_xt
    outer = edge_xt + [corner] + far_xt[::-1]
    v.fill("region", outer)
    v.fill("region", dust)
    v.fill("star", dust)
    v.fill("cover", outer)
    chi = np.linspace(0, CHI_DISC, 200)
    for time in (0.4, 0.8, 1.1, 1.3, 1.45):
        e = math.atanh(math.sin(time))
        v.curve("t2", np.arctan(e - chi), np.arctan(e + chi))
    for r in (R0 / 3, 2 * R0 / 3):
        c = 2 * math.asin(r / 2)
        v.curve("r2", np.arctan(eta - c), np.arctan(eta + c))

    def ext_curve(cls, p, q):
        # Outside the edge, q >= arctan(tan p + 2 chi_0) on the rays that leave it, and after the moment of rest.
        with np.errstate(invalid="ignore"):
            keep = (p + q > -1e-9) & ((p < -math.atan(CHI_DISC)) | (q > np.arctan(np.tan(p) + 2 * CHI_DISC) - 1e-9))
        v.curve(cls, np.where(keep, p, np.nan), np.where(keep, q, np.nan))
    T = cf.spread(-np.inf, np.inf, 2000, 11)
    for R in (0.75, 1.5, 3.0):
        ext_curve("r", *outside(T, np.full_like(T, R)))
    ext_curve("r", *inner_region(T, np.full_like(T, 0.25)))
    radii = KAPPA + np.geomspace(1e-9, 1e6, 2000)
    for time in (0.5, 1.0, 2.0, 4.0):
        ext_curve("t", *outside(np.full_like(radii, time), radii))
    v.line("event", [[[0.0, 2 * p0], corner]])
    trapped_eta, trapped_chi = np.arctanh(np.sin(t)), 2 * np.arcsin(trapped / 2)
    v.curve("apparent", np.arctan(trapped_eta - trapped_chi), np.arctan(trapped_eta + trapped_chi))
    v.line("surface", [edge_xt])
    v.line("centre", [[[0.0, 0.0], [0.0, PI]]])
    v.line("singular", [[[0.0, PI], corner]], zig=True)
    v.line("boundary", [far_xt])
    v.line("chartedge", [[[0.0, 0.0], [PI, 0.0]]])
    v.label_xt([corner[0] / 2, PI], "$R = 0$", "b", dy=-8)
    v.label_xt([2.9, 1.4], "$R \\to \\infty$", "l", "small", dx=6)
    v.label_xt([PI / 2, 0], "$t = 0$ inside, $T = 0$ outside", "t", "coord", dy=7)
    v.label_xt([0.6, 0.35], "dust", cls="region")
    v.label_xt([2.2, 0.9], "exterior", cls="region")
    v.legend("star", "the dust")
    v.legend("cover", "the exterior")
    v.legend("surface", "the edge")
    v.legend("t2", "$t$ constant")
    v.legend("r2", "$r$ constant")
    v.legend("r", "$R$ constant")
    v.legend("t", "$T$ constant")
    v.legend("event", "the event horizon")
    v.legend("apparent", "marginally trapped circles")
    v.legend("singular", "$R = 0$")
    v.legend("boundary", "the boundary")
    v.legend("centre", "the centre")
    v.legend("chartedge", "the moment of rest")
    v.set(settings="$k = 1/(4\\ell^2)$, $r_0 = \\ell$, and $M = 1/4$.")
    from slices import moments
    for m in moments(METRIC, "rest"):
        c = np.linspace(0, CHI_DISC, 101)
        lo, hi = m.reach(DISC_OUT, "R")
        root = np.sqrt((np.linspace(lo, hi, 401) - KAPPA) / (np.linspace(lo, hi, 401) + KAPPA))
        v.slice(m, [(np.arctan(-c), np.arctan(c)), (disc_P(-root), disc_Q(root))])
    return [v]


# -- Curvature figure ---------------------------------------------------------------------------------

def register_curvature(cv):
    cv.DRAWN[METRIC] = lambda ck, src: curvature(cv, ck, src)


def curvature(cv, ck, src):
    """K and R of the disc against its proper time from the stored scalars of the comoving chart at
    a = cos t, ell = 1: K = 4(2 a^4 + (1 + k - a^2)^2)/a^4 and R = 2(1 + k - 3 a^2)/a^2, the same at every r,
    with the exterior's stored scalars, 12 and -6, checked beside them."""
    params = {"k": K}
    names = list(params)
    scale_factor = {"a": "cos(t)"}
    Kp = cv.Published(src, METRIC, DISC, "kretschmann", ["t"], names, functions=scale_factor)
    Rp = cv.Published(src, METRIC, DISC, "ricci_scalar", ["t"], names, functions=scale_factor)
    Ko = cv.Published(src, METRIC, DISC_OUT, "kretschmann", ["R"], ["ell"])
    Ro = cv.Published(src, METRIC, DISC_OUT, "ricci_scalar", ["R"], ["ell"])
    x_axis = cv.axis("$ct/\\ell$", 0, 1.5)
    K_scale = {"kind": "log", "from": 1, "to": 1e6}
    R_scale = {"kind": "symlog", "from": -10, "to": 1e3, "linear": 1}
    panels = [{"id": "K", "name": "$K\\,\\ell^4$", "aspect": 2, "scale": K_scale, "ticks": cv.value_ticks(K_scale)},
              {"id": "R", "name": "$R\\,\\ell^2$", "aspect": 2, "scale": R_scale, "ticks": cv.value_ticks(R_scale)}]
    cv.shape_checks(ck, METRIC, [p["aspect"] for p in panels])
    height = {p["id"]: cv.PLOT_W / p["aspect"] for p in panels}
    curves = (cv.line_curves(ck, METRIC, "K", "K", Kp, "t", x_axis, K_scale, height["K"], params)
              + cv.line_curves(ck, METRIC, "R", "R", Rp, "t", x_axis, R_scale, height["R"], params))
    ts = np.linspace(0, 1.5, cv.SAMPLES)
    a = np.cos(ts)
    ck("Anti-de Sitter dust collapse: K = 4(2 a^4 + (1 + k - a^2)^2)/a^4 at a = cos t",
       bool(np.allclose(Kp(ts, K), 4 * (2 * a ** 4 + (1 + K - a * a) ** 2) / a ** 4, rtol=1e-12)))
    ck("Anti-de Sitter dust collapse: R = 2(1 + k - 3 a^2)/a^2 at a = cos t",
       bool(np.allclose(Rp(ts, K), 2 * (1 + K - 3 * a * a) / (a * a), rtol=1e-12, atol=1e-12)))
    ck("Anti-de Sitter dust collapse: K = 33/4 and R = -7/2 at the release",
       abs(float(Kp.exact([0, K])) - 8.25) < 1e-12 and abs(float(Rp.exact([0, K])) + 3.5) < 1e-12)
    ck("Anti-de Sitter dust collapse: R vanishes where cos^2 t = (1 + k)/3",
       abs(float(Rp.exact([math.acos(math.sqrt((1 + K) / 3)), K]))) < 1e-12)
    ck("Anti-de Sitter dust collapse: outside the disc K = 12 and R = -6",
       abs(float(Ko.exact([2.0, 1.0])) - 12) < 1e-12 and abs(float(Ro.exact([2.0, 1.0])) + 6) < 1e-12)
    legend = [["line", "K", "$K$"], ["line", "R", "$R$"],
              ["line", "horizon", "the edge crosses the horizon"]]
    return [{
        "id": "collapse", "label": "The collapse", "kind": "line", "x": x_axis,
        "panels": panels, "curves": curves, "marks": [{"class": "horizon", "x": T_HORIZON}],
        "legend": legend,
        "caption": ["The Kretschmann scalar $K$ and the Ricci scalar $R$ of the dust against its proper time. "
                    "$K$ rises from $33/(4\\ell^4)$ and diverges at the crunch. Outside the disc $K = 12/\\ell^4$. "
                    "$R$ rises from $-7/(2\\ell^2)$ through zero at $ct = 0.87\\,\\ell$. "
                    "Outside the disc $R = -6/\\ell^2$."],
        "ricci_sentences": [[0, 3], [0, 4]],
        "per_scalar": cv.scalar_clips("A logarithmic scale from $1$ to $10^6$.",
                                      "A signed logarithmic scale from $-10$ to $10^3$, linear within $1$."),
        "settings": "$k = 1/(4\\ell^2)$.",
        "clip": "$K$ logarithmic; $R$ signed logarithmic, linear within $1$.",
    }]
