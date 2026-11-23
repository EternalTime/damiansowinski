"""The drawings of the generalised Vaidya spacetime.

Every drawing is of one sourced mass function, in units of a length L, with the null time scaled so
that the mass settles within the drawing.

The ingoing chart draws Husain's collapse of radiation, P = rho/3, his eq. (21) at k = 1/3 with C = 0,
2m/r = A(1 + tanh)/r + B(1 + tanh)/r^(2/3), at A = L/2 and B = L^(2/3)/8 and with the tanh of 2cv/L:

    m = (1 + tanh(2cv/L)) (L/4 + L^(2/3) r^(1/3)/16).

It is flat space as v -> -infinity and a black hole as v -> +infinity, whose apparent horizon is
Husain's cubic (r - 2A)^3 = (2B)^3 r, r = 1.27079 L. Husain's chart draws his hairy black hole, his
eq. (20) at k = 3/4 with A = L and B = 1/2: 2f = A(1 + tanh) and 4(q/r)^(3/2) = (1 - B tanh)/r^(3/2),

    f = L(1 + tanh(2cv/L))/2,   q = L((1 - tanh(2cv/L)/2)/4)^(2/3),

a naked singularity at early times, with 1 + (3/2)(L/r)^(3/2) in the lapse, and a black hole with hair
at late times, 1 - 2L/r + (L/r)^(3/2)/2, whose horizons are r = 1.60538 L and 0.06690 L. The outgoing
chart draws Glass and Krisch's diffusing strings, their eq. (18a) with rho_0 = 0,

    m = m_0 + 2 pi k_1 (r^2 - 2 D c u),   m_0 = L/2, k_1 = 1/(80 pi L), D = L/2:

Schwarzschild's mass L/2 with strings whose density k_1/r falls with the radius and which carry mass
inward as the light streams out. The free mass function leaves the global completion undetermined, so
no conformal diagram is drawn, as for Mallett's radiating star. generalised_vaidya.md beside this file
is the derivation.
"""
import numpy as np

METRIC = "generalised_vaidya"

COLLAPSE = "(1 + tanh(2*v))*(1/4 + r**(1/3)/16)"
HAIR = {"f": "(1 + tanh(2*v))/2", "q": "((1 - tanh(2*v)/2)/4)**(2/3)"}
STRINGS = "1/2 + (r**2 - u)/40"
COLLAPSE_HORIZON = 1.27078909342472
HAIR_HORIZONS = (1.60537794047960, 0.0669008687454887)
UNIT = "with $L$ the unit of every length"
COLLAPSE_INPUT = ("Husain's collapse of radiation, $P = \\rho/3$: $m = (1 + \\tanh(2cv/L))(L/4 + L^{2/3}r^{1/3}/16)$, "
                  + UNIT + ".")
HAIR_INPUT = ("Husain's hairy black hole at $k = 3/4$: $f = L(1 + \\tanh(2cv/L))/2$ and "
              "$q = L\\left((1 - \\tanh(2cv/L)/2)/4\\right)^{2/3}$, " + UNIT + ".")
STRINGS_INPUT = ("Glass and Krisch's diffusing strings: $m = m_0 + 2\\pi k_1(r^2 - 2Dcu)$ with $m_0 = L/2$, "
                 "$k_1 = 1/80\\pi L$, and $D = L/2$, " + UNIT + ".")


def register_null(nr):
    # Near the naked singularity of Husain's early times g^rr = 1 - 2m/r exceeds 2, where cv - r is no time,
    # so that view orients its cones by the ingoing rays, which every ingoing chart has future directed.
    nr.DIAGRAMS.extend([
        nr.Diagram(METRIC, "ingoing", "collapse", "radiation collapsing", ("v", "r"), (0, 4, -3, 1.5),
                   "$r/L$", "$(cv - r)/L$", {}, nr.EQUATOR, to_display=nr.FINKELSTEIN_IN, tau="v - r",
                   areal=True, functions={"m": COLLAPSE}, singular_runs=True,
                   marked=(("event", {"x0": "2.6", "r": repr(COLLAPSE_HORIZON)}, 1, "the event horizon"),),
                   input=COLLAPSE_INPUT),
        nr.Diagram(METRIC, "husain", "hair", "a black hole growing hair", ("v", "r"), (0, 4, -3, 1.5),
                   "$r/L$", "$(cv - r)/L$", {"k": "3/4"}, nr.EQUATOR, to_display=nr.FINKELSTEIN_IN, tau="v - r",
                   orient="ingoing", areal=True, functions=dict(HAIR), singular_runs=True,
                   marked=(("event", {"x0": "2.6", "r": repr(HAIR_HORIZONS[0])}, 1, "the event horizon"),),
                   input=HAIR_INPUT),
        nr.Diagram(METRIC, "outgoing", "strings", "diffusing strings", ("u", "r"), (0, 4, -1.5, 3),
                   "$r/L$", "$(cu + r)/L$", {}, nr.EQUATOR, to_display=nr.FINKELSTEIN_OUT, tau="u + r",
                   areal=True, functions={"m": STRINGS}, singular_runs=True, input=STRINGS_INPUT),
    ])
    nr.CAPTIONS.update({
        (METRIC, "ingoing", "collapse"): [
            "The plane of $v$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), with $cv - r$ vertical. Radiation "
            "collapses from nearly flat space to a black hole whose horizon settles at $r = 1.27L$."],
        (METRIC, "husain", "hair"): [
            "The plane of $v$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), with $cv - r$ vertical. A naked "
            "singularity at early times, and a black hole with hair once the fluid has fallen in, with horizons "
            "at $r = 1.61L$ and $0.067L$."],
        (METRIC, "outgoing", "strings"): [
            "The plane of $u$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), with $cu + r$ vertical. Light "
            "streams out of a central mass $m_0 = L/2$ through an atmosphere of strings."],
    })


# ---------------------------------------------------------------- the embedding diagram
#
# Husain's collapse of radiation, the ingoing chart with COLLAPSE. A slice of constant v is null, so
# the moments are slices of constant v - r = T, as Vaidya's are, spacelike wherever m > -r/2: with
# v = T + r the published metric pulls back to (1 + 2m/r) dr^2 + r^2 dphi^2, so dz/dr = sqrt(2m/r).
# At every finite v the centre holds a mass, m(v, 0) = (1 + tanh(2v))/4, and the surface closes there
# in a spike. The marginally trapped sphere on each moment is the root of r = 2m(T + r, r). The moments
# start at cv - r = -L: before it the central mass is below a hundredth of L and its spike narrower than
# the drawing resolves.

EMBED_TOP = 4.0
EMBED_TIMES = [-1.0 + k / 16 for k in range(41)]
EMBED_SHOWN = (0, 16, 28, 40)


def collapse_mass(v, r):
    return (1 + np.tanh(2 * v)) * (0.25 + np.cbrt(r) / 16)


def trapped_on(T):
    """The marginally trapped sphere on the moment v - r = T, r = 2m(T + r, r), by bisection."""
    lo, hi = 1e-12, EMBED_TOP
    f = lambda r: r - 2 * collapse_mass(T + r, r)  # noqa: E731
    if f(lo) >= 0 or f(hi) <= 0:
        return None
    for _ in range(200):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if f(mid) < 0 else (lo, mid)
    return (lo + hi) / 2


def embedding_views(ck, src, em):
    top, size = EMBED_TOP, 2 * EMBED_TOP

    def moment(T):
        where = f"generalised_vaidya, v - r = {T:g}"
        sl = em.Slice(src, METRIC, "ingoing", "r", "\\phi", {"theta": "pi/2"}, {}, {"m": COLLAPSE},
                      along={"v": f"{T!r} + r"})
        rings = [(x, "r", None) for x in (2.0, 3.0, top)]
        trapped = trapped_on(T)
        marks = ([(trapped, "horizon", None)] if trapped is not None else []) + rings
        whole = em.Piece("whole", "sheet", sl, 0.0, top, 0.0, 1,
                         (("apex", "the singularity $r = 0$, where the surface closes in a spike"),
                          ("edge", "the surface runs on to $r \\to \\infty$")), marks, size)
        if trapped is not None:
            ck.add(f"{where}, the marginally trapped sphere is r = 2m", abs(trapped - 2 * collapse_mass(T + trapped, trapped)),
                   1e-12)
        # The rim of the drawing, r = 4L, stands at z = 0 at every moment.
        whole.z = whole.z - whole.z[-1]
        ck.isometry(f"{where}, whole", whole)
        ck.radius(f"{where}, areal radius", whole, lambda r: r, size)
        return em.Surface([whole], label=f"$cv - r = {T:g}\\,L$", time=T)

    frames = [moment(T) for T in EMBED_TIMES]
    fig = em.movie_figure(frames, {"sheet": "cover"}, size)
    fig.legend("fill", "cover", "the slice of constant $cv - r$")
    fig.legend("line", "r", "$r$ constant, at $2L$, $3L$ and $4L$")
    fig.legend("line", "horizon", "the marginally trapped sphere, $r = 2m$")
    fig.legend("line", "meridian", "$\\phi$ constant, every $30°$")
    return [em.view("collapse", "Radiation collapsing", "$L$", [frames[k] for k in EMBED_SHOWN], fig.done(),
                    movie=em.movie(frames, "$cv - r$", EMBED_TIMES),
                    settings="$L = 1$; each moment is a slice of constant $cv - r$.", input=COLLAPSE_INPUT)]


def embedding_captions():
    return {(METRIC, "collapse"): [
        "Husain's radiation collapsing to a black hole, at moments of constant $cv - r$. Each moment closes in "
        "a spike at $r = 0$, and the marginally trapped sphere grows to $r = 1.27L$."]}


def slice_marks(slices):
    """Each moment on the spacetime diagram of the ingoing chart, v = T + r in its (cv, r), out to the
    embedding's reach."""
    return slices.one(METRIC, lambda m: [[(m.time + r, r) for r in m.reach("ingoing", "r")]])


HIDDEN_REASON = ("another spacetime than the collapse of radiation embedded: Husain's hairy black hole in his "
                 "own chart, and Glass and Krisch's strings in the outgoing chart")


# ---------------------------------------------------------------- the curvature figure
#
# Husain's collapse of radiation in the ingoing chart, the published Kretschmann and Ricci scalars over
# the plane of r and v, with the marginally trapped spheres, the zero of the published g^rr. The Ricci
# scalar is (4 d_r m + 2r d_r^2 m)/r^2 = (1 + tanh(2v)) r^(-8/3)/18 here, the trace of the fluid, which
# vanishes for the null dust alone.

CURVATURE_SCALES = {"kretschmann": {"kind": "log", "from": 1e-4, "to": 1e4},
                    "ricci_scalar": {"kind": "log", "from": 1e-4, "to": 1e2}}
CURVATURE_CELLS = 100


def curvature_views(ck, src, cv):
    n = CURVATURE_CELLS
    functions = {"m": COLLAPSE}
    x_axis, y_axis = cv.axis("$r/L$", 0.0, 3.0), cv.axis("$cv/L$", -1.5, 1.5)
    cv.shape_checks(ck, METRIC, [(x_axis["to"] - x_axis["from"]) / (y_axis["to"] - y_axis["from"])])
    grr = cv.nr_published(src, METRIC, "ingoing", "inverse_metric_components", ["r", "v"], [], ("r", "r"), functions)
    track = np.linspace(y_axis["from"], y_axis["to"], 1201)
    curve = []
    for t in track:
        zeros = cv._g_zeros(lambda r, at: grr(r, at), t, x_axis["to"])
        if len(zeros) == 1:
            curve.append([zeros[0], float(t)])
    ck(f"{METRIC} the marginally trapped sphere is r = 2m at every moment drawn",
       len(curve) > 1000 and all(abs(p[0] - 2 * collapse_mass(p[1], p[0])) < 1e-6 for p in curve))
    Xp = cv.PLOT_W * np.array([p[0] for p in curve]) / x_axis["to"]
    Yp = cv.PLOT_W * (1 - (np.array([p[1] for p in curve]) - y_axis["from"]) / (y_axis["to"] - y_axis["from"]))
    curve = [curve[k] for k in cv.thinned(Xp, Yp, cv.TOLERANCE)]
    out = []
    for field, label in (("kretschmann", "$K\\,L^4$"), ("ricci_scalar", "$R\\,L^2$")):
        scale = CURVATURE_SCALES[field]
        published = cv.nr_published(src, METRIC, "ingoing", field, ["r", "v"], [], functions=functions)
        xs, ys, cells = cv.cells(published, x_axis, y_axis, n, n, [])
        X, Y = np.meshgrid(xs, ys)
        ck(f"{METRIC} {field} is finite at every cell", np.all(np.isfinite(cells)))
        if field == "ricci_scalar":
            ck(f"{METRIC} R is (1 + tanh(2v)) r^(-8/3)/18 at every cell",
               np.allclose(cells, (1 + np.tanh(2 * Y)) * X ** (-8 / 3) / 18, rtol=1e-10, atol=0))
        rng = np.random.default_rng(1996)
        picks = sorted({(int(i), int(j)) for i, j in zip(rng.integers(0, n, 40), rng.integers(0, n, 40))}
                       | {(0, n - 1), (0, n // 2), (n - 1, n - 1)})
        cv.agrees(ck, f"{METRIC} {field} map", published, [[xs[i], ys[j]] for i, j in picks],
                  [cells[j, i] for i, j in picks])
        word = "Kretschmann scalar $K$" if field == "kretschmann" else "Ricci scalar $R$"
        out.append({
            "id": field, "label": "$K$" if field == "kretschmann" else "$R$", "system": "ingoing", "kind": "map",
            "x": x_axis, "y": y_axis,
            "colour": {"name": label, "scale": scale, "ticks": cv.value_ticks(scale)},
            "map": {"metric": METRIC, "system": "ingoing", "scalar": field, "x": "r", "y": "v", "params": {},
                    "functions": functions, "nx": n, "ny": n, "cells": [[cv.rounded(c) for c in row] for row in cells]},
            "marks": [{"class": "horizon", "points": curve}],
            "legend": [["line", "horizon", "the marginally trapped spheres, $r = 2m$"]],
            "caption": [f"The {word} of Husain's collapse of radiation over the plane of $v$ and $r$."],
            "settings": "$L = 1$.", "input": COLLAPSE_INPUT,
            "clip": ("Logarithmic, $10^{-4}$ to $10^4$." if field == "kretschmann" else "Logarithmic, $10^{-4}$ to $10^2$."),
        })
    return out
