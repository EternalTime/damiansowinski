"""The drawings of Krasnikov's time machine, at m = 1/2, h = 1/5 and n = 1, with R = 4/5 and the scale
factor a = 1 at tau = 0, the unit of every length.

The plane of tau and chi is conformally flat in both charts, f^-2 and a^2 times the flat plane, so its
light rays keep tau -+ chi for every f and every a: each spacetime diagram is drawn with a declared
function and checked to give the null directions of the flat plane. The time machine's f is
(sin^2(5 pi (tau + 1/5)/2) + (chi^2 - 1/25)^2)/12R, which vanishes in the closed square only at
(+-h, +-h), repeats near the upper cut its values near the lower one, eq. (4), since sin^2 has the
period 2h, and keeps the right side of eq. (6) below 0.92/R^2, so that R^-2 bounds it.

The moment tau = 0 of the universe on the great sphere theta = pi/2 has a = 1 and rho = R on the
square, where f = 1, the plateau eq. (10a) and the join with eq. (10b) ask for when f is constant near
the square's edge; rho = cos chi beyond |chi| = n, and between the two the polynomial of degree five
that matches rho, rho' and rho'' at both ends, which bends down, rho'' <= 0, and keeps |rho'| <= sin 1,
so Krasnikov's eq. (11) holds and the moment lies in flat space: a cylinder of radius R round the
equator between two caps of the unit sphere.

Neither plane is given a conformal diagram: tau and chi are conformal coordinates of finite range on
both, so each spacetime diagram is the conformal diagram of its plane. No curvature figure is drawn,
since the scalars are those of f, a and rho, which Krasnikov fixes only by inequalities.
krasnikov_time_machine.md beside this file records each decision.
"""
import math

import numpy as np

METRIC = "krasnikov_time_machine"
M, H, N, RADIUS = 0.5, 0.2, 1.0, 0.8
EQUATOR = {"theta": "pi/2", "phi": "0"}
VIEW = "moment"
F_TIME_MACHINE = "(sin(5*pi*(tau + 1/5)/2)**2 + (chi**2 - 1/25)**2)/(12*R)"
A_UNIVERSE = "exp(sin(tau + 2*chi)/2)"
SETTINGS = "$m = 1/2$, $h = 1/5$, and $n = 1$."


def _quintic():
    """The coefficients, lowest power first, of the polynomial on [m, n] equal to R at m with no
    slope or curvature there, and to cos chi at n with its slope and second derivative."""
    rows, values = [], [RADIUS, 0.0, 0.0, math.cos(N), -math.sin(N), -math.cos(N)]
    for x in (M, N):
        rows.append([x ** i for i in range(6)])
        rows.append([i * x ** (i - 1) if i else 0.0 for i in range(6)])
        rows.append([i * (i - 1) * x ** (i - 2) if i > 1 else 0.0 for i in range(6)])
    return np.linalg.solve(np.array(rows), np.array(values))


QUINTIC = _quintic()


def rho_hat(chi, derivative=0):
    """rho at tau = 0 and its derivatives along chi: R on |chi| <= m, the quintic out to n and
    cos chi beyond, even in chi."""
    chi = np.asarray(chi, dtype=float)
    x = np.abs(chi)
    poly = np.polynomial.polynomial
    join = poly.polyval(x, poly.polyder(QUINTIC, derivative) if derivative else QUINTIC)
    outer = {0: np.cos(x), 1: -np.sin(x), 2: -np.cos(x)}[derivative]
    inner = RADIUS if derivative == 0 else 0.0
    value = np.where(x <= M, inner, np.where(x < N, join, outer))
    return value * (np.sign(chi) if derivative == 1 else 1.0)


# ---------------------------------------------------------------- spacetime diagrams

def _removed():
    return tuple(("removed", (repr(s * H), repr(c * H)), "the removed points $(\\pm h, \\pm h)$")
                 for s in (1, -1) for c in (1, -1))


CUT = "the cuts, each upper bank glued to the other's lower bank"


def register_null(nr):
    cuts = tuple(("surface", "x0", repr(s * H), CUT, (repr(-H), repr(H))) for s in (1, -1))
    nr.DIAGRAMS.append(nr.Diagram(
        METRIC, "time_machine", "square", "$\\tau$ and $\\chi$", ("\\tau", "\\chi"), (-M, M, -M, M),
        "$\\chi$", "$\\tau$", {"R": "1", "h": "1/5", "m": "1/2"}, EQUATOR, tau="tau", families=nr.SIDEWAYS,
        functions={"f": F_TIME_MACHINE}, any_factor="f", lines=cuts, points=_removed(),
        input="Drawn with $f = (\\sin^2(\\pi(\\tau + h)/2h) + (\\chi^2 - h^2)^2)/12R$; the rays are the same for any $f$."))
    nr.CLOSED_FORMS[(METRIC, "time_machine", "square")] = (lambda t, x: t + x, lambda t, x: t - x, lambda t, x: np.isfinite(t))
    nr.CAPTIONS[(METRIC, "time_machine", "square")] = [
        "The plane of $\\tau$ and $\\chi$ ($\\theta = \\pi/2$, $\\phi = 0$) of the time machine, $f^{-2}$ times the "
        "flat plane, with light rays at 45°. The causality violating set is the hexagon "
        "$|\\tau| < h$, $|\\chi| < 2h - |\\tau|$."]
    edge = "the square $|\\tau|, |\\chi| < m$"
    square = tuple(("shell", which, repr(s * M), edge, (repr(-M), repr(M))) for which in ("x0", "r") for s in (1, -1))
    friedmann = tuple(("world", "r", repr(s * N), "$|\\chi| = n$, the edge of Friedmann's universe")
                      for s in (1, -1))
    cuts = tuple(("surface", "x0", repr(s * H), CUT, (repr(-H), repr(H))) for s in (1, -1))
    nr.DIAGRAMS.append(nr.Diagram(
        METRIC, "universe", "plane", "$\\tau$ and $\\chi$", ("\\tau", "\\chi"), (-math.pi / 2, math.pi / 2, -1, 1),
        "$\\chi$", "$\\tau$", {"h": "1/5", "m": "1/2", "n": "1"}, EQUATOR, tau="tau", families=nr.SIDEWAYS,
        functions={"a": A_UNIVERSE}, any_factor="a", kretschmann=False, lines=square + cuts + friedmann,
        points=_removed(),
        input="Drawn with $a = e^{\\sin(\\tau + 2\\chi)/2}$; the rays are the same for any $a$."))
    nr.CLOSED_FORMS[(METRIC, "universe", "plane")] = (lambda t, x: t + x, lambda t, x: t - x, lambda t, x: np.isfinite(t))
    nr.CAPTIONS[(METRIC, "universe", "plane")] = [
        "The plane of $\\tau$ and $\\chi$ ($\\theta = \\pi/2$, $\\phi = 0$) of the universe, $a^2$ times the flat "
        "plane, with light rays at 45°."]


def slice_table(slices):
    """The moment tau = 0 of the embedding view, on both planes."""
    moment = "$\\tau = 0$"
    return {(METRIC, "time_machine", "square"): (lambda: slices.one(
                METRIC, lambda m: slices.along(0.0, -M, M), label=moment, view_id=VIEW)),
            (METRIC, "universe", "plane"): (lambda: slices.one(
                METRIC, lambda m: slices.along(0.0, *m.reach("universe", "\\chi")), label=moment, view_id=VIEW))}


# ---------------------------------------------------------------- embedding diagram

def register_embedding(em):
    em.DRAWN[METRIC] = lambda ck, src: embedding_views(ck, src, em)
    em.CAPTIONS[(METRIC, VIEW)] = [
        "The moment $\\tau = 0$ of the universe ($\\theta = \\pi/2$), a cylinder of radius $R$ round its equator "
        "between two caps of a sphere of radius $a$."]


def embedding_views(ck, src, em):
    size = 2.0
    numeric = {"rho": (lambda x: rho_hat(x), lambda x: rho_hat(x, 1))}
    sl = em.Slice(src, METRIC, "universe", "\\chi", "\\phi", {"tau": 0, "theta": "pi/2"},
                  {"h": "1/5", "m": "1/2", "n": "1"}, functions={"a": "1"}, numeric=numeric)
    xs = np.linspace(-math.pi / 2, math.pi / 2, 4001)
    ck.add("krasnikov_time_machine: rho bends down at tau = 0, rho'' <= 0, eq. (11)",
           float(max(0.0, np.max(rho_hat(xs, 2)))), 1e-12)
    ck.add("krasnikov_time_machine: |rho'| <= 1 at tau = 0, eq. (11)", float(max(0.0, np.max(np.abs(rho_hat(xs, 1))) - 1)),
           1e-12)
    for x in (M, N):
        for k in (0, 1, 2):
            jump = abs(float(rho_hat(x - 1e-9, k)) - float(rho_hat(x + 1e-9, k)))
            ck.add(f"krasnikov_time_machine: rho and its derivative {k} join at |chi| = {x}", jump, 1e-6)
    sheet = em.Piece("moment", "sheet", sl, -math.pi / 2, math.pi / 2, 0.0, 1,
                     (("axis", "the pole $\\chi = -\\pi/2$"), ("axis", "the pole $\\chi = \\pi/2$")),
                     [(-N, "r", None), (-M, "surface", None), (M, "surface", "$\\chi = m$"), (N, "r", "$\\chi = n$")],
                     size, knots=tuple(np.linspace(-math.pi / 2, math.pi / 2, 41)[1:-1]))
    ck.isometry("krasnikov_time_machine, the moment tau = 0", sheet)
    ck.radius("krasnikov_time_machine, circles of radius a rho", sheet, lambda x: rho_hat(x), size)
    surface = em.Surface([sheet])
    fig = em.figure_of([surface], {"sheet": "cover"}, size)
    em.ring_label(fig, [0, 0, 0], *sheet.at(M), "$\\chi = m$")
    em.ring_label(fig, [0, 0, 0], *sheet.at(N), "$\\chi = n$")
    fig.legend("fill", "cover", "the moment $\\tau = 0$")
    fig.legend("line", "surface", "$|\\chi| = m$, the time machine's edge")
    fig.legend("line", "r", "$|\\chi| = n$")
    fig.legend("line", "meridian", "$\\phi$ constant, every $15°$")
    return [em.view(VIEW, "The moment $\\tau = 0$", "$a$", [surface], fig.done(),
                    settings="$a = 1$ at $\\tau = 0$, the unit of length, $m = 1/2$, $n = 1$, and "
                             "$\\rho = R = 4/5$ on $|\\chi| \\le m$, joined to $\\cos\\chi$ by a polynomial of "
                             "degree five.")]


# ---------------------------------------------------------------- conformal diagram

def register_conformal(cf):
    # tau and chi are conformal coordinates of finite range on both planes, so each spacetime
    # diagram is already the conformal diagram of its plane; krasnikov_time_machine.md says so.
    cf.NOT_DRAWN.add(METRIC)
