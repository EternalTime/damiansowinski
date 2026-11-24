"""The drawings of the Hawking-Turok instanton.

The instanton and the two instantons in closed form are positive definite, so those charts have no
nonzero null tangent and no causal boundary: no spacetime diagram. The two Lorentzian charts, outside
the bubble wall and the open universe, leave b(sigma) and a(t) free, and no source prints in closed
form a b that reaches the singularity: the one closed form, a field at rest, b = ell sin(sigma/ell), is
de Sitter space, whose planes the Hawking-Moss and Coleman-De Luccia pages draw, and it ends on a
second regular pole where the instanton ends on a timelike singularity, which light outside the wall
reaches at a finite tau since int d sigma/b converges there for b ~ (sigma_f - sigma)^{1/3}. So
neither plane is drawn, and for the same reason no conformal diagram.

The embedding diagram draws the surface of chi and phi at psi = theta = pi/2 of the exponential
instanton at ell = 1. Its circles have the radius rho = ell sqrt(cos chi) sin chi, and in the proper
distance s, ds = ell sqrt(cos chi) d chi, d rho/ds = cos chi - sin^2 chi/(2 cos chi), which is 1 at the
regular pole, 0 at the widest circle, tan chi = sqrt 2, and -1 at cos chi = 1/3, the root of
3c^2 + 2c - 1 = 0; there the surface lies level, and nearer the singularity the circles shrink faster
than the distance between them, so the surface runs from the pole to cos chi = 1/3. It draws as well
the surface of T and phi of Vilenkin's instanton at A = 1, whose circles have the radius
rho = T^{1/2}(A^2 + T^2)^{1/4} with d rho/ds = (A^2 + 2T^2)/(2T sqrt(A^2 + T^2)) > 1 at every T, since
(A^2 + 2T^2)^2 - 4T^2(A^2 + T^2) = A^4: every circle grows faster than the distance out to it, and the
surface lies in Minkowski space, on its light cone at the singularity and level toward flat space.

The curvature figure is the stored Kretschmann and Ricci scalars of the exponential instanton against
chi and of Vilenkin's instanton against T. hawking_turok_instanton.md records each choice.
"""
import math

import numpy as np
import sympy as sp

import hawking_turok_instanton as hti

METRIC = hti.METRIC
EUCLIDEAN = ("instanton", "exponential", "vilenkin")
POSITIVE = "Positive definite: no nonzero null tangent or causal boundary."
FREE = ("The function {f} is left free, and no source prints in closed form one that reaches the singularity; "
        "the field at rest, de Sitter space, is drawn on the {page} page.")
NOT_DRAWN = {"instanton": POSITIVE, "exponential": POSITIVE, "vilenkin": POSITIVE,
             "bubble_exterior": FREE.format(f="$b(\\sigma)$", page="Hawking-Moss"),
             "open_universe": FREE.format(f="$a(t)$", page="Coleman-De Luccia")}
STOP = math.acos(1 / 3)                     # cos chi = 1/3, where the exponential instanton's surface lies level
WIDEST = math.atan(math.sqrt(2))            # tan chi = sqrt 2, its widest circle
VILENKIN_TOP = 2.0                          # how far out Vilenkin's surface is drawn, T/A
VILENKIN_LO = 0.05                          # how near the singularity it starts, T/A, where it nears the light cone
POLE_SETTINGS = "$\\ell = 1$, the unit of every length."
VILENKIN_SETTINGS = "$A = 1$, the unit of every length."


def pole_rho(chi):
    """The radius ell sqrt(cos chi) sin chi of the circles of the exponential instanton at ell = 1."""
    chi = np.asarray(chi, dtype=float)
    return np.sqrt(np.cos(chi)) * np.sin(chi)


def vilenkin_rho(T):
    """The radius T^{1/2}(A^2 + T^2)^{1/4} of the circles of Vilenkin's instanton at A = 1."""
    T = np.asarray(T, dtype=float)
    return np.sqrt(T) * (1 + T ** 2) ** 0.25


# ---------------------------------------------------------------- spacetime and conformal diagrams

def register_null(nr):
    nr.NOT_DRAWN.update({(METRIC, system): reason for system, reason in NOT_DRAWN.items()})


def register_conformal(cf):
    cf.NOT_DRAWN.add(METRIC)


# ---------------------------------------------------------------- embedding diagrams

def register_embedding(eb):
    eb.DRAWN[METRIC] = lambda ck, src: embedding(eb, ck, src)
    eb.CAPTIONS[(METRIC, "pole")] = [
        "The surface of $\\chi$ and $\\phi$ ($\\psi = \\theta = \\pi/2$) of the exponential instanton, from its "
        "regular pole out to $\\cos\\chi = 1/3$, where it lies level."]
    eb.CAPTIONS[(METRIC, "vilenkin")] = [
        "The surface of $T$ and $\\phi$ ($\\psi = \\theta = \\pi/2$) of Vilenkin's instanton in Minkowski space, "
        "nearing the light cone toward the singularity $T = 0$ and level toward flat space."]


def embedding(eb, ck, src):
    """The surface of chi and phi of the exponential instanton in Euclidean space, and of T and phi of
    Vilenkin's instanton in Minkowski space, each at psi = theta = pi/2."""
    fixed_at = {"psi": "pi/2", "theta": "pi/2"}
    views = []

    name = "exponential instanton"
    sl = eb.Slice(src, METRIC, "exponential", "\\chi", "\\phi", fixed_at, {"ell": 1})
    ck.add(f"{name}: d rho/ds = cos chi - sin^2 chi/(2 cos chi) is -1 at cos chi = 1/3",
           abs(math.cos(STOP) - math.sin(STOP) ** 2 / (2 * math.cos(STOP)) + 1), 1e-12)
    ck.add(f"{name}: d rho/ds is 0 at the widest circle, tan chi = sqrt 2",
           abs(math.cos(WIDEST) - math.sin(WIDEST) ** 2 / (2 * math.cos(WIDEST))), 1e-12)
    ck.stops(f"{name}, the surface of chi and phi past cos chi = 1/3", sl,
             np.linspace(STOP + 1e-6, math.pi / 2 - 1e-6, 400))
    size = 2 * float(pole_rho(WIDEST))
    rings = (math.pi / 8, math.pi / 4, 3 * math.pi / 8)
    cap = eb.Piece("cap", "sheet", sl, 0.0, STOP, 0.0, -1,
                   (("axis", "the regular pole $\\chi = 0$"),
                    ("stops", "at $\\cos\\chi = 1/3$ the surface lies level, and nearer the singularity the circles "
                              "shrink faster than the distance between them")),
                   [(x, "r", None) for x in rings] + [(STOP, "chartedge", None)], size)
    ck.isometry(f"{name}, the cap", cap)
    ck.radius(f"{name}: rho = l sqrt(cos chi) sin chi", cap, pole_rho, size)
    surface = eb.Surface([cap])
    fig = eb.figure_of([surface], {"sheet": "cover"}, size)
    eb.ring_label(fig, [0, 0, 0], *cap.at(math.pi / 4), "$\\chi = \\pi/4$")
    fig.legend("fill", "cover", "the surface of $\\chi$ and $\\phi$")
    fig.legend("line", "r", "$\\chi$ constant, at $\\pi/8$, $\\pi/4$, and $3\\pi/8$")
    fig.legend("line", "chartedge", "$\\cos\\chi = 1/3$, where the surface lies level and the drawing stops")
    fig.legend("line", "meridian", "$\\phi$ constant, every $15°$")
    views.append(eb.view("pole", "Exponential potential", "$\\ell$", [surface], fig.done(), settings=POLE_SETTINGS,
                         stops=["Between $\\cos\\chi = 1/3$ and the singularity the circles shrink faster than the "
                                "distance between them, so no surface of revolution in Euclidean space has the geometry "
                                "of that part of the slice."]))

    name = "Vilenkin's instanton"
    sv = eb.Slice(src, METRIC, "vilenkin", "T", "\\phi", fixed_at, {"A": 1}, space="minkowski")
    flat = eb.Slice(src, METRIC, "vilenkin", "T", "\\phi", fixed_at, {"A": 1})
    ck.stops(f"{name}, the surface of T and phi in flat space", flat, np.linspace(1e-4, 20, 400))
    A, T = sp.symbols("A T", positive=True)
    ck.exact(f"{name}: (A^2 + 2T^2)^2 - 4T^2(A^2 + T^2) = A^4, so d rho/ds > 1",
             sp.expand((A ** 2 + 2 * T ** 2) ** 2 - 4 * T ** 2 * (A ** 2 + T ** 2) - A ** 4) == 0)
    size = 2 * float(vilenkin_rho(VILENKIN_TOP))
    rings = (0.25, 0.5, 1.0, 1.5, VILENKIN_TOP)
    sheet = eb.Piece("sheet", "sheet", sv, VILENKIN_LO, VILENKIN_TOP, 0.0, 1,
                     (("edge", "the surface runs on to the singularity $T = 0$, nearing the light cone"),
                      ("edge", "the surface runs on toward flat space, $T \\to \\infty$")),
                     [(x, "r", None) for x in rings], size, digits=1e-12,
                     knots=tuple(float(k) for k in np.linspace(VILENKIN_LO, VILENKIN_TOP, 1441)[1:-1]))
    ck.isometry(f"{name}, the sheet", sheet)
    ck.radius(f"{name}: rho = T^(1/2)(A^2 + T^2)^(1/4)", sheet, vilenkin_rho, size)
    surface = eb.Surface([sheet])
    fig = eb.figure_of([surface], {"sheet": "cover"}, size)
    eb.ring_label(fig, [0, 0, 0], *sheet.at(1.0), "$T = A$")
    eb.ring_label(fig, [0, 0, 0], *sheet.at(VILENKIN_TOP), "$2A$")
    fig.legend("fill", "cover", "the surface of $T$ and $\\phi$")
    fig.legend("line", "r", "$T$ constant, at $A/4$, $A/2$, $A$, $3A/2$, and $2A$")
    fig.legend("line", "meridian", "$\\phi$ constant, every $15°$")
    views.append(eb.view("vilenkin", "Vilenkin's instanton", "$A$", [surface], fig.done(),
                         settings=VILENKIN_SETTINGS + " Every length measured with $dX^2 + dY^2 - dZ^2$.",
                         stops=["Every circle grows faster than the distance out to it, so no surface of revolution "
                                "in Euclidean space has this geometry; one in Minkowski space does."]))
    return views


# ---------------------------------------------------------------- curvature figure

def register_curvature(cv):
    cv.DRAWN[METRIC] = lambda ck, src: curvature(cv, ck, src)


def _view(cv, ck, src, system, coordinate, label, x_name, lo, hi, unit, K_scale, R_scale, caption, settings,
          params, checks, clips):
    names = list(params)
    Kp = cv.Published(src, METRIC, system, "kretschmann", [coordinate], names)
    Rp = cv.Published(src, METRIC, system, "ricci_scalar", [coordinate], names)
    x_axis = cv.axis(x_name, lo, hi)
    panels = [{"id": "K", "name": f"$K\\,{unit}^4$", "aspect": 2, "scale": K_scale, "ticks": cv.value_ticks(K_scale)},
              {"id": "R", "name": f"$R\\,{unit}^2$", "aspect": 2, "scale": R_scale, "ticks": cv.value_ticks(R_scale)}]
    cv.shape_checks(ck, METRIC, [p["aspect"] for p in panels])
    height = {p["id"]: cv.PLOT_W / p["aspect"] for p in panels}
    curves = (cv.line_curves(ck, METRIC, "K", "K", Kp, coordinate, x_axis, K_scale, height["K"], params)
              + cv.line_curves(ck, METRIC, "R", "R", Rp, coordinate, x_axis, R_scale, height["R"], params))
    checks(Kp, Rp, list(params.values()))
    return {
        "id": system, "label": label, "kind": "line", "system": system, "x": x_axis,
        "panels": panels, "curves": curves, "marks": [],
        "legend": [["line", "K", "$K$"], ["line", "R", "$R$"]],
        "caption": caption, "ricci_sentences": [[0, 2]],
        "per_scalar": cv.scalar_clips(*clips),
        "settings": settings,
        "clip": clips[0] + " " + clips[1],
    }


def curvature(cv, ck, src):
    """K and R from the stored scalars: of the exponential instanton against chi from its regular pole
    to the singularity chi = pi/2, and of Vilenkin's instanton against T from the singularity T = 0 out
    toward flat space. Both scalars are positive throughout, finite at the regular pole, where
    K = 96/l^4 and R = 24/l^2, and growing without bound toward each singularity, where R s^2 tends to
    2/3 at the proper distance s from it, as Turok and Hawking found; Vilenkin's fall off as 15A^8/4T^12
    and 3A^4/2T^6 far out."""
    def exponential_checks(Kp, Rp, vals):
        x = np.linspace(0.0, 1.55, 400)
        ck(f"{METRIC} exponential: K = 96/l^4 at the regular pole", abs(float(Kp(np.array([0.0]), *vals)[0]) - 96) < 1e-9)
        ck(f"{METRIC} exponential: R = 24/l^2 at the regular pole", abs(float(Rp(np.array([0.0]), *vals)[0]) - 24) < 1e-9)
        ck(f"{METRIC} exponential: K and R are positive and grow toward the singularity",
           bool(np.all(Kp(x, *vals) > 0) and np.all(Rp(x, *vals) > 0)
                and np.all(np.diff(Kp(x, *vals)) > 0) and np.all(np.diff(Rp(x, *vals)) > 0)))
        gap = 1e-4
        s = 2 / 3 * gap ** 1.5                  # the proper distance to chi = pi/2 at l = 1, to leading order
        ck(f"{METRIC} exponential: R s^2 tends to 2/3 at the singularity",
           abs(float(Rp(np.array([math.pi / 2 - gap]), *vals)[0]) * s ** 2 - 2 / 3) < 1e-3)

    def vilenkin_checks(Kp, Rp, vals):
        x = np.linspace(0.01, 4.0, 400)
        ck(f"{METRIC} vilenkin: K and R are positive and fall outward",
           bool(np.all(Kp(x, *vals) > 0) and np.all(Rp(x, *vals) > 0)
                and np.all(np.diff(Kp(x, *vals)) < 0) and np.all(np.diff(Rp(x, *vals)) < 0)))
        far = np.array([1e3])
        ck(f"{METRIC} vilenkin: K falls off as 15A^8/4T^12", abs(float(Kp(far, *vals)[0]) * 4e36 / 15 - 1) < 1e-5)
        ck(f"{METRIC} vilenkin: R falls off as 3A^4/2T^6", abs(float(Rp(far, *vals)[0]) * 2e18 / 3 - 1) < 1e-5)
        T = 1e-6
        s = 2 / 3 * T ** 1.5                    # int_0^T sqrt(T/sqrt(1 + T^2)) dT at A = 1, to leading order
        ck(f"{METRIC} vilenkin: R s^2 tends to 2/3 at the singularity",
           abs(float(Rp(np.array([T]), *vals)[0]) * s ** 2 - 2 / 3) < 1e-6)

    return [
        _view(cv, ck, src, "exponential", "chi", "Exponential potential", "$\\chi$", 0.0, math.pi / 2, "\\ell",
              {"kind": "log", "from": 10, "to": 1e8}, {"kind": "log", "from": 10, "to": 1e5},
              ["The Kretschmann scalar $K$ and the Ricci scalar $R$ of the exponential instanton against $\\chi$, "
               "up to the singularity $\\chi = \\pi/2$. $K$ is $96/\\ell^4$ at the regular pole. $R$ is $24/\\ell^2$ "
               "there and grows as $2/(3s^2)$ at a proper distance $s$ from the singularity."],
              POLE_SETTINGS, {"ell": 1.0}, exponential_checks,
              ("A logarithmic scale from $10$ to $10^8$.", "A logarithmic scale from $10$ to $10^5$.")),
        _view(cv, ck, src, "vilenkin", "T", "Vilenkin's instanton", "$T/A$", 0.0, 4.0, "A",
              {"kind": "log", "from": 1e-8, "to": 1e8}, {"kind": "log", "from": 1e-5, "to": 1e5},
              ["The Kretschmann scalar $K$ and the Ricci scalar $R$ of Vilenkin's instanton against $T$, out from "
               "the singularity $T = 0$. $K$ falls off as $15A^8/(4T^{12})$. $R$ grows as $2/(3s^2)$ at a proper "
               "distance $s$ from the singularity and falls off as $3A^4/(2T^6)$."],
              VILENKIN_SETTINGS, {"A": 1.0}, vilenkin_checks,
              ("A logarithmic scale from $10^{-8}$ to $10^8$.", "A logarithmic scale from $10^{-5}$ to $10^5$.")),
    ]
