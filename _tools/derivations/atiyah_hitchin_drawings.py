"""The drawings of the Atiyah-Hitchin space.

The metric is positive definite, so neither chart has a nonzero null tangent or a causal boundary:
no spacetime diagram and no conformal diagram. The embedding diagrams are the two geodesic surfaces
of Atiyah and Hitchin's book that Bachmaier, Dvali, Seitz and Valbuena-Bermudez (2025) describe, the
smoothed cone and the trumpet, drawn in Gibbons and Manton's radius with a, b and c evaluated from
their solution at each r. The curvature figure is the Kretschmann scalar of the modulus chart.

Each circle on these surfaces closes after an angle pi once the two identifications are made, so
each is drawn swept by twice its angle, and its radius is half the coefficient of the form."""
import math

import mpmath
import numpy as np

import atiyah_hitchin as ah

METRIC = ah.METRIC
REASON = "Positive definite: no nonzero null tangent or causal boundary."
NOT_DRAWN = {system: REASON for system in ah.SYSTEMS}
TOP = 6.0               # the trumpet is drawn out to r = 6n
CONE_TOP = 10.0         # and the cone, which leaves the core almost flat, to r = 10n
SETTINGS = "$n = 1$, the unit of every length."


def modulus_of(r):
    """The modulus k of Gibbons and Manton's radius at n = 1, the root of 2K(k^2) = r."""
    if r <= math.pi:
        return mpmath.mpf(0)
    return mpmath.findroot(lambda q: 2 * mpmath.ellipk(q) - r, (mpmath.mpf("1e-30"), 1 - mpmath.mpf("1e-30")),
                           solver="anderson") ** 0.5


def solution(r):
    """a, b and c at the radius r, n = 1, from the products bc = -4KE, ca = -4K(E - (1 - k^2)K) and
    ab = 4K(K - E), with their first derivatives along r from the equations of self duality."""
    mpmath.mp.dps = max(mpmath.mp.dps, 40)
    r = mpmath.mpf(r)
    if r <= math.pi:
        a, b, c = mpmath.mpf(0), mpmath.pi, -mpmath.pi
        return (a, b, c), (mpmath.mpf(2), mpmath.mpf("0.5"), mpmath.mpf("0.5"))
    m = modulus_of(r) ** 2
    K, E = mpmath.ellipk(m), mpmath.ellipe(m)
    ab, bc, ca = 4 * K * (K - E), -4 * K * E, -4 * K * (E - (1 - m) * K)
    a, b, c = mpmath.sqrt(ab * ca / bc), mpmath.sqrt(ab * bc / ca), -mpmath.sqrt(bc * ca / ab)
    return (a, b, c), ah.rates(a, b, c, r)


def numeric(index):
    """The function a, b or c of r and its derivative, as numpy functions for a slice."""
    def f(r):
        return np.array([float(solution(x)[0][index]) for x in np.atleast_1d(r)]).reshape(np.shape(r))

    def df(r):
        return np.array([float(solution(x)[1][index]) for x in np.atleast_1d(r)]).reshape(np.shape(r))
    return f, df


FUNCTIONS = {name: numeric(i) for i, name in enumerate("abc")}


def register_null(nr):
    nr.NOT_DRAWN.update({(METRIC, system): reason for system, reason in NOT_DRAWN.items()})


def register_conformal(cf):
    cf.NOT_DRAWN.add(METRIC)


def register_embedding(eb):
    eb.DRAWN[METRIC] = lambda ck, src: embedding(eb, ck, src)
    eb.CAPTIONS[(METRIC, "cone")] = [
        "The surface of $r$ and $\\phi$ ($\\theta = \\pi/2$, $\\psi = 0$), where the two monopoles stay in one plane: "
        "smooth at the core, far away a cone of half angle $30°$."]
    eb.CAPTIONS[(METRIC, "trumpet")] = [
        "The surface through the core ($\\phi = 0$), $\\psi$ turning below it and $\\theta$ above: "
        "a tube of radius $n$ at one end, a cone of half angle $30°$ at the other."]


def embedding(eb, ck, src):
    size = 3.5
    marks = (4.0, 5.0)
    cone_slice = eb.Slice(src, METRIC, "radius", "r", "\\phi", {"theta": "pi/2", "psi": 0}, {"n": 1},
                          numeric=FUNCTIONS, swept={"phi": "phi/2"})
    cone = eb.Piece("cone", "sheet", cone_slice, math.pi, CONE_TOP, 0.0, 1,
                    (("axis", "the core $r = \\pi n$, where the circle of $\\phi$ closes to a point"),
                     ("edge", "the surface runs on to $r \\to \\infty$, a cone of half angle $30°$")),
                    [(r, "r", None) for r in (4.0, 6.0, 8.0)] + [(CONE_TOP, "r", "$r = 10n$")], size)
    ck.isometry("Atiyah-Hitchin, the cone", cone)
    ck.radius("Atiyah-Hitchin, the cone, rho = a/2", cone, lambda r: FUNCTIONS["a"][0](r) / 2, size)
    far = solution(30)[0]
    ck.add("Atiyah-Hitchin: at r = 30n, a and b are r sqrt(1 - 2n/r) and c is -2n/sqrt(1 - 2n/r), Taub-NUT's "
           "with a negative mass, to the exponentially small corrections",
           max(abs(float(far[0] / (30 * math.sqrt(1 - 2 / 30))) - 1), abs(float(far[1] / (30 * math.sqrt(1 - 2 / 30))) - 1),
               abs(float(-far[2] * math.sqrt(1 - 2 / 30) / 2) - 1)), 1e-9)
    cone_surface = eb.Surface([cone])
    fig = eb.figure_of([cone_surface], {"sheet": "cover"}, size)
    eb.ring_label(fig, [0, 0, 0], *cone.at(CONE_TOP), "$10n$")
    fig.legend("fill", "cover", "the plane of the two monopoles, $\\theta = \\pi/2$, $\\psi = 0$")
    fig.legend("line", "r", "$r$ constant, at $4$, $6$, $8$ and $10\\,n$")
    fig.legend("line", "meridian", "$\\phi$ constant")
    views = [eb.view("cone", "The plane of the pair", "$n$", [cone_surface], fig.done(),
                     settings=SETTINGS + " Each circle is drawn swept by $2\\phi$.")]

    tube_slice = eb.Slice(src, METRIC, "radius", "r", "\\psi", {"theta": "pi/2", "phi": 0}, {"n": 1},
                          numeric=FUNCTIONS, swept={"psi": "psi/2"})
    flare_slice = eb.Slice(src, METRIC, "radius", "r", "\\theta", {"phi": 0, "psi": 0}, {"n": 1},
                           numeric=FUNCTIONS, swept={"theta": "theta/2"})
    tube = eb.Piece("tube", "sheet", tube_slice, math.pi, TOP, 0.0, -1,
                    (("throat", "the core $r = \\pi n$, where the two halves meet"),
                     ("edge", "the surface runs on to $r \\to \\infty$, a tube of radius $n$")),
                    [(math.pi, "throat", "$r = \\pi n$")] + [(r, "r", None) for r in marks] + [(TOP, "r", None)], size)
    flare = eb.Piece("flare", "sheet2", flare_slice, math.pi, TOP, 0.0, 1,
                     (("throat", "the core $r = \\pi n$"),
                      ("edge", "the surface runs on to $r \\to \\infty$, a cone of half angle $30°$")),
                     [(r, "r2", None) for r in marks] + [(TOP, "r2", "$r = 6n$")], size)
    for piece, index in ((tube, 2), (flare, 1)):
        ck.isometry(f"Atiyah-Hitchin, the trumpet, {piece.id}", piece)
        ck.radius(f"Atiyah-Hitchin, the trumpet, {piece.id}, rho = |{'cb'[index - 1]}|/2", piece,
                  lambda r, i=index: np.abs(FUNCTIONS["abc"[i]][0](r)) / 2, size)
    ck.join("Atiyah-Hitchin, the two halves of the trumpet at the core", tube, math.pi, flare, math.pi)
    trumpet_surface = eb.Surface([tube, flare])
    fig = eb.figure_of([trumpet_surface], {"sheet": "cover", "sheet2": "cover"}, size)
    eb.ring_label(fig, [0, 0, 0], *tube.at(math.pi), "$r = \\pi n$", dx=14)
    eb.ring_label(fig, [0, 0, 0], *flare.at(TOP), "$6n$")
    fig.legend("fill", "cover", "the surface through the core")
    fig.legend("line", "r", "$r = 4$, $5$, $6\\,n$ on the half of $\\psi$")
    fig.legend("line", "r2", "the same on the half of $\\theta$")
    fig.legend("line", "throat", "the core $r = \\pi n$")
    fig.legend("line", "meridian", "the angle constant")
    views.append(eb.view("trumpet", "The trumpet", "$n$", [trumpet_surface], fig.done(),
                         settings=SETTINGS + " Each circle is drawn swept by twice its angle."))
    return views


def register_curvature(cv):
    cv.DRAWN[METRIC] = lambda ck, src: curvature(cv, ck, src)


# The modulus is sampled from k = 0.02, where r - pi n = 3 x 10^-4 n, to k = 0.9995, where r = 9.7n.
K_FROM, K_TO = 0.02, 0.9995
K_LINE = {"kind": "linear", "from": 0.0, "to": 0.6}


class Precise:
    """A stored scalar of the modulus chart evaluated in sixty digits wherever it is sampled.

    The stored Kretschmann scalar is a rational function of k, K and E whose numerator and
    denominator both vanish to high order at the core, where K - E and E - (1 - k^2)K are of
    order k^2: in doubles it loses every digit below k = 0.6, and in sixty digits it holds to a
    part in 10^19 at k = 0.02. The value drawn is the stored one; only its arithmetic is longer."""

    def __init__(self, published, digits=60):
        import sympy as sp
        self.published, self.digits = published, digits
        self.metric, self.system, self.field = published.metric, published.system, published.field
        self.functions, self.expr = published.functions, published.expr
        self.f = sp.lambdify(published.coords + published.params, published.expr, "mpmath")

    def __call__(self, x, *params):
        with mpmath.workdps(self.digits):
            values = [float(self.f(mpmath.mpf(float(v)), *[mpmath.mpf(float(p)) for p in params]))
                      for v in np.atleast_1d(x)]
        return np.array(values).reshape(np.shape(x))

    def exact(self, values, digits=60):
        return self.published.exact(values, digits)


def curvature(cv, ck, src):
    """The Kretschmann scalar of the modulus chart against k at n = 1, from k = 0.02, effectively
    the core, out to k = 0.9995, r = 9.7n. The stored Ricci scalar of the modulus chart vanishes."""
    params = {"n": 1.0}
    stored = cv.Published(src, METRIC, "modulus", "kretschmann", ["k"], ["n"])
    ck(f"{METRIC} the modulus chart's Ricci scalar vanishes",
       cv.Published(src, METRIC, "modulus", "ricci_scalar", ["k"], ["n"]).expr == 0)
    K = Precise(stored)
    core, far = float(K(K_FROM, 1.0)), float(K(K_TO, 1.0))
    ck(f"{METRIC} K is largest near the core and falls away from it",
       core > float(K(0.5, 1.0)) > float(K(0.9, 1.0)) > far > 0, f"{core:.6f} at k = {K_FROM}, {far:.3e} at k = {K_TO}")
    near = float(K(0.03, 1.0))
    ck(f"{METRIC} K is flat at the core, the same at k = 0.02 and 0.03 to a part in 10^3",
       abs(near / core - 1) < 1e-3, f"{near:.8f} against {core:.8f}")
    x_axis = cv.axis("$k$", 0.0, 1.0)
    panels = [{"id": "K", "name": "$K\\,n^4$", "aspect": 2, "scale": K_LINE, "ticks": cv.value_ticks(K_LINE)}]
    cv.shape_checks(ck, f"{METRIC} line", [p["aspect"] for p in panels])
    curves = cv.line_curves(ck, f"{METRIC} modulus", "K", "K", K, "k", x_axis, K_LINE, cv.PLOT_W / 2, params,
                            within=(K_FROM, K_TO))
    up, down = cv.beyond(K, x_axis, K_LINE, [1.0], within=(K_FROM, K_TO))
    ck(f"{METRIC} the curve stays inside its panel", not up and not down)
    view = {
        "id": "modulus", "label": "the modulus $k$", "kind": "line", "x": x_axis, "panels": panels,
        "curves": curves, "marks": [],
        "legend": [["line", "K", "$K$"]],
        "caption": ["The Kretschmann scalar from the core $k = 0$, where it is largest, out to $k = 0.9995$, "
                    "Gibbons and Manton's radius $r = 9.7n$; the Ricci scalar vanishes."],
        "settings": "$n = 1$, elliptic modulus chart.",
        "clip": "Linear scale $0$ to $0.6$, none cut.",
    }
    return [view]
