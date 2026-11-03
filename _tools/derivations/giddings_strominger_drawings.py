"""The drawings of the Giddings-Strominger wormhole.

The metric is positive definite, so no chart has a nonzero null tangent or a causal boundary: no
spacetime diagram and no conformal diagram. The embedding diagram is the surface of r and phi at
chi = theta = pi/2 of the areal chart at a = 1, both halves through the throat r = a, a funnel with
dz/dr = a^2/sqrt(r^4 - a^4) whose two sheets level off at the height
sqrt(pi) Gamma(1/4) a/(4 Gamma(3/4)) = 1.3110 a, z = z_inf - (a^2/r) 2F1(1/2, 1/4; 5/4; a^4/r^4). The
curvature figure is the two stored scalars of the conformal chart, the one chart that covers both
halves, against eta."""
import mpmath
import numpy as np
import sympy as sp

import giddings_strominger as gs

METRIC = gs.METRIC
REASON = "Positive definite: no nonzero null tangent or causal boundary."
NOT_DRAWN = {system: REASON for system in gs.SYSTEMS}
TOP = 4.0               # both sheets are drawn out to r = 4a
RADII = (1.5, 2.0, 3.0)


def register_null(nr):
    nr.NOT_DRAWN.update({(METRIC, system): reason for system, reason in NOT_DRAWN.items()})


def register_conformal(cf):
    cf.NOT_DRAWN.add(METRIC)


def register_embedding(eb):
    eb.DRAWN[METRIC] = lambda ck, src: embedding(eb, ck, src)
    eb.CAPTIONS[(METRIC, "funnel")] = [
        "The surface of $r$ and $\\phi$ through the throat $r = a$, both halves out to $r = 4a$: a funnel whose "
        "sheets level off far out at the height $1.311a$ above and below the throat."]


def funnel(r):
    """The height of the funnel above the throat at a = 1, z_inf - 2F1(1/2, 1/4; 5/4; 1/r^4)/r."""
    from scipy.special import hyp2f1
    r = np.asarray(r, dtype=float)
    return float(gs.height()) - hyp2f1(0.5, 0.25, 1.25, np.minimum(r ** -4.0, 1.0)) / r


def embedding(eb, ck, src):
    with mpmath.workdps(30):
        z_inf = gs.height()
        # r = 1 + u^2 takes the inverse square root at the throat out of the integrand.
        quadrature = mpmath.quad(lambda u: 2 / mpmath.sqrt((2 + u ** 2) * ((1 + u ** 2) ** 2 + 1)), [0, 1, mpmath.inf])
        ck.add("Giddings-Strominger: the height at infinity is sqrt(pi) Gamma(1/4)/(4 Gamma(3/4)) = 1.3110287771 a",
               abs(float(z_inf) - 1.3110287771), 1e-10)
        ck.add("Giddings-Strominger: the height at infinity is the quadrature of a^2/sqrt(r^4 - a^4)",
               abs(float(quadrature - z_inf)), 1e-25)
    ck.add("Giddings-Strominger: the closed form of z vanishes at the throat", abs(float(funnel(1.0))), 1e-12)
    size = 2 * TOP
    sl = eb.Slice(src, METRIC, "areal", "r", "\\phi", {"chi": "pi/2", "theta": "pi/2"}, {"a": 1})
    throat = sl.horizons()[0]
    ck.add("Giddings-Strominger: g_rr diverges at the throat r = a", abs(throat - 1), 1e-12)
    end = "the surface runs on to $r \\to \\infty$, where it levels off"
    near = eb.Piece("near", "sheet", sl, throat, TOP, 0.0, 1,
                    (("throat", "the throat $r = a$, where the two halves meet"), ("edge", end)),
                    [(throat, "throat", "$r = a$")] + [(r, "r", None) for r in RADII] + [(TOP, "r", "$r = 4a$")],
                    size)
    far = eb.Piece("far", "sheet2", sl, throat, TOP, 0.0, -1,
                   (("throat", "the throat $r = a$"), ("edge", end)),
                   [(r, "r2", None) for r in RADII] + [(TOP, "r2", None)], size)
    for p in (near, far):
        ck.isometry(f"Giddings-Strominger, {p.id}", p)
        ck.form(f"Giddings-Strominger, {p.id}, the closed form of z", p, lambda r, s=p.sense: s * funnel(r), size)
        ck.radius(f"Giddings-Strominger, {p.id}, rho = r", p, lambda r: r, size)
    ck.join("Giddings-Strominger, the two halves at the throat", near, throat, far, throat)
    surface = eb.Surface([near, far])
    fig = eb.figure_of([surface], {"sheet": "cover", "sheet2": "cover"}, size)
    eb.ring_label(fig, [0, 0, 0], throat, 0.0, "$r = a$", dx=14)
    eb.ring_label(fig, [0, 0, 0], *near.at(TOP), "$4a$")
    fig.legend("fill", "cover", "the surface $\\chi = \\theta = \\pi/2$")
    fig.legend("line", "r", "$r$ constant, at $1.5$, $2$, $3$ and $4\\,a$")
    fig.legend("line", "r2", "the same radii on the other half")
    fig.legend("line", "throat", "the throat $r = a$")
    fig.legend("line", "meridian", "$\\phi$ constant, every $15°$")
    return [eb.view("funnel", "Through the throat", "$a$", [surface], fig.done(),
                    settings="$a = 1$, the unit of every length.")]


def register_curvature(cv):
    cv.DRAWN[METRIC] = lambda ck, src: curvature(cv, ck, src)


ETA = 1.5               # |eta| <= 1.5 reaches the areal radius a sqrt(cosh 3) = 3.17a on either side


def curvature(cv, ck, src):
    """The two stored scalars of the conformal chart against eta at a = 1, from eta = -1.5 through
    the throat eta = 0 to 1.5, r = 3.17a on either side; they depend on eta alone, so a line shows
    their whole variation. K is logarithmic from 10^-5 to 10^2 and R signed logarithmic from -10 to
    zero, linear within 10^-3; both curves lie wholly inside their panels.

    Checks hold every sample to the stored scalar in sixty digits, the conformal factor to
    r^2 = a^2 cosh 2 eta, both scalars to their closed forms 60/(a^4 cosh^6 2eta) and
    -6/(a^2 cosh^3 2eta), which are 60a^8/r^12 and -6a^4/r^6, their even symmetry, the extrema
    60/a^4 and -6/a^2 at the throat, the monotone fall of K and rise of R outward, and the unit
    scaling."""
    x_axis = cv.axis("$\\eta$", -ETA, ETA)
    K_scale = {"kind": "log", "from": 1e-5, "to": 1e2}
    R_scale = {"kind": "symlog", "from": -10, "to": 0, "linear": 1e-3}
    panels = [{"id": "K", "name": "$K\\,a^4$", "aspect": 2, "scale": K_scale, "ticks": cv.value_ticks(K_scale)},
              {"id": "R", "name": "$R\\,a^2$", "aspect": 2, "scale": R_scale, "ticks": cv.value_ticks(R_scale)}]
    cv.shape_checks(ck, METRIC, [p["aspect"] for p in panels])
    K = cv.nr_published(src, METRIC, "conformal", "kretschmann", ["\\eta"], ["a"])
    R = cv.nr_published(src, METRIC, "conformal", "ricci_scalar", ["\\eta"], ["a"])
    sphere = cv.nr_published(src, METRIC, "conformal", "metric_components", ["\\eta"], ["a"], ("\\chi", "\\chi"))
    eta, a = K.coords[0], K.params[0]
    ck(f"{METRIC} the 3-sphere at eta has the areal radius a sqrt(cosh 2eta)",
       sp.simplify(sphere.expr - a ** 2 * sp.cosh(2 * eta)) == 0)
    ck(f"{METRIC} K is 60/(a^4 cosh^6 2eta), 60a^8/r^12",
       sp.simplify(K.expr - 60 / (a ** 4 * sp.cosh(2 * eta) ** 6)) == 0)
    ck(f"{METRIC} R is -6/(a^2 cosh^3 2eta), -6a^4/r^6",
       sp.simplify(R.expr + 6 / (a ** 2 * sp.cosh(2 * eta) ** 3)) == 0)
    ck(f"{METRIC} K at the throat is finite and greatest, 60/a^4", float(K.exact([0, 1])) == 60)
    ck(f"{METRIC} R at the throat is least, -6/a^2", float(R.exact([0, 1])) == -6)
    xs = np.linspace(-ETA, ETA, cv.SAMPLES)
    curves = []
    for name, fn, scale, power in (("K", K, K_scale, 4), ("R", R, R_scale, 2)):
        cv.agrees(ck, f"{METRIC} all {cv.SAMPLES} {name} samples", fn, [[x, 1] for x in xs], fn(xs, 1))
        ck(f"{METRIC} {name} is the same on both halves", sp.simplify(fn.expr.subs(eta, -eta) - fn.expr) == 0)
        ck(f"{METRIC} {name} has the declared length scaling",
           all(np.allclose(fn(xs, unit) * unit ** power, fn(xs, 1), rtol=1e-12, atol=0) for unit in (0.5, 2)))
        outward = np.diff(fn(xs[xs >= 0], 1))
        ck(f"{METRIC} {name} {'falls' if name == 'K' else 'rises'} outward on either half",
           np.all(outward < 0) if name == "K" else np.all(outward > 0))
        curves += cv.line_curves(ck, METRIC, name, name, fn, "\\eta", x_axis, scale, cv.PLOT_W / 2, {"a": 1})
        ck(f"{METRIC} {name} stays inside its panel", not any(cv.beyond(fn, x_axis, scale, [1])))
    k_clip = "We draw $K\\,a^4$ on a logarithmic scale from $10^{-5}$ to $10^{2}$; every value drawn lies within it."
    r_clip = ("We draw $R\\,a^2$ on a signed logarithmic scale, linear within $10^{-3}$ of zero, from $-10$ to $0$; "
              "every value drawn lies within it.")
    return [{
        "id": "throat", "label": "both halves of the throat", "kind": "line", "x": x_axis,
        "panels": panels, "curves": curves, "marks": [],
        "legend": [["line", "K", "$K$, the Kretschmann scalar"], ["line", "R", "$R$, the Ricci scalar"]],
        "caption": [
            "The Kretschmann scalar $K$ and the Ricci scalar $R$ of the Giddings-Strominger wormhole against the "
            "conformal coordinate $\\eta$ ($a = 1$), across both halves from the throat $\\eta = 0$ to the areal "
            "radius $r = 3.17a$. $K$ is greatest at the throat, $60/a^4$.",
            "$R$ is negative everywhere and least at the throat, $-6/a^2$.",
        ],
        "ricci_sentences": [[1, 0]], "per_scalar": cv.scalar_clips(k_clip, r_clip),
        "settings": "$a = 1$, the unit of every length.",
        "clip": k_clip + " " + r_clip,
    }]
