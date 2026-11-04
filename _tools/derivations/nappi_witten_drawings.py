"""Draw the cosmological region and its round moment of maximum volume."""

import math

METRIC = "nappi_witten"
MOMENT = r"$\psi=\pi/4$"
PARAMS = {"L": 1, "b": 1}


def diagram_specs(nr):
    """The cosmological null planes and a periodic-time plane in the extension."""
    half = math.pi/2
    rows = [
        nr.Diagram(METRIC, "cosmological", "polar", "Cosmological region", (r"\psi", "s"),
                   (-0.12, half+0.12, -0.12, half+0.12), "$s$", r"$\psi$", PARAMS,
                   {"rho": "0", "lambda": "0"}, tau="psi", families=nr.SIDEWAYS),
        nr.Diagram(METRIC, "bang_regular", "polar", "Cosmological region", ("x", "s"),
                   (-0.12, half+0.12, -0.10, 3.35), "$s$", "$x$", PARAMS,
                   {"rho": "0"}, tau="x", families=nr.SIDEWAYS, quotient="w", where="x",
                   cone="future cone of zero angular momentum"),
        nr.Diagram(METRIC, "invariant", "polar", "Cosmological region", ("r", "s"),
                   (-0.12, half+0.12, -0.12, 1.12), "$s$", "$r$", PARAMS,
                   {"rho": "0", "lambda": "0"}, tau="r", families=nr.SIDEWAYS,
                   where="r*(1-r)"),
        nr.Diagram(METRIC, "whisker", "periodic", "Periodic time", (r"\lambda", r"\chi"),
                   (-0.15, 6.45, -0.15, 2*math.pi+0.15), r"$\chi$", r"$\lambda$", PARAMS,
                   {"rho": "0", "s": "pi/4"}, tau="lambda_", families=nr.SIDEWAYS,
                   periodic=(r"\lambda",)),
    ]
    return rows


def diagram_captions():
    return {
        (METRIC, "cosmological", "polar"):
            [r"Null geodesics across the closed universe ($\rho=\lambda=0$, $b=L=1$)."],
        (METRIC, "bang_regular", "polar"):
            [r"The cosmological region in $x=\tan^2\psi$, with zero momentum around $w$ ($\rho=0$, $b=L=1$)."],
        (METRIC, "invariant", "polar"):
            [r"The expanding and recollapsing region in $r=\sin^2\psi$ ($\rho=\lambda=0$, $b=L=1$)."],
        (METRIC, "whisker", "periodic"):
            [r"Null curves at fixed $s$ in the static extension, with periodic time $\lambda$ ($s=\pi/4$, $\rho=0$, $b=L=1$)."],
    }


def null_forms():
    import numpy as np
    return {
        (METRIC, "cosmological", "polar"):
            (lambda t, s: t+s, lambda t, s: t-s, None),
        (METRIC, "bang_regular", "polar"):
            (lambda x, s: np.arctan(np.sqrt(x))+s, lambda x, s: np.arctan(np.sqrt(x))-s, None),
        (METRIC, "invariant", "polar"):
            (lambda r, s: np.arcsin(np.sqrt(r))+s, lambda r, s: np.arcsin(np.sqrt(r))-s, None),
        (METRIC, "whisker", "periodic"):
            (lambda t, x: t+np.log(np.tanh(x/2)), lambda t, x: t-np.log(np.tanh(x/2)), None),
    }


def embedding_views(ck, src, em):
    """A hemisphere cut from the round spatial sphere at maximum volume."""
    import numpy as np
    sl = em.Slice(src, METRIC, "cosmological", "s", r"\lambda",
                  {"psi": math.pi/4, "rho": 0}, PARAMS)
    piece = em.Piece("hemisphere", "sheet", sl, 0, math.pi/2, -1, 1,
                     (("axis", "$s=0$"), ("edge", r"$s=\pi/2$")),
                     [(s, "r", None) for s in (math.pi/6, math.pi/3)], 2.5)
    ck.isometry("Nappi-Witten round section", piece)
    ck.radius("Nappi-Witten circle radius", piece, np.sin, 2.5)
    ck.form("Nappi-Witten hemisphere height", piece, lambda s: -np.cos(s), 2.5)
    surface = em.Surface([piece], label=MOMENT, time=math.pi/4)
    fig = em.figure_of([surface], {"sheet": "cover"}, 2.5)
    fig.legend("fill", "cover", "Round spatial section")
    fig.legend("line", "r", "$s$ constant")
    fig.legend("line", "meridian", r"$\lambda$ constant")
    data = fig.done()
    data["box"] = [-1.4, 1.4, -1.4, 1.4]
    return [em.view("turnaround", "Maximum volume", "$L$", [surface], data,
                    settings=r"$b=L=1$, $\psi=\pi/4$, $\rho=0$.")]


def embedding_captions():
    return {(METRIC, "turnaround"):
            [r"A round hemispherical section at maximum spatial volume ($\psi=\pi/4$, $b=1$)."]}


def slice_marks(slices, system):
    clock = {"cosmological": math.pi/4, "bang_regular": 1, "invariant": 0.5}[system]
    def lines(moment):
        lo, hi = moment.reach("cosmological", "s")
        return slices.along(clock, lo, hi)
    return slices.one(METRIC, lines, label=MOMENT, view_id="turnaround")


def conformal_views(ck, src, cf):
    """The finite cosmological causal plane, with both angular circles suppressed."""
    import numpy as np
    half = math.pi/2
    plane = cf.Plane(src, METRIC, "cosmological", (r"\psi", "s"),
                     {"rho": "0", "lambda": "0"}, PARAMS)
    null_map = lambda t, s: (t-s, t+s)
    ck.chart("Nappi-Witten cosmological causal plane", plane, null_map,
             ck.uniform(0.01, half-0.01), ck.uniform(0.01, half-0.01), lambda t, s: (1, 0))
    view = cf.View("cosmological", "Closed cosmological region",
                   [-0.18, half+0.18, -0.18, half+0.18], "cosmological")
    rectangle = [[0, 0], [half, 0], [half, half], [0, half]]
    view.fill("region", rectangle)
    view.fill("cover", rectangle)
    view.line("centre", [[[0, 0], [0, half]], [[half, 0], [half, half]]])
    view.line("boundary", [[[0, 0], [half, 0]], [[0, half], [half, half]]])
    for intercept in np.linspace(-half, half, 9)[1:-1]:
        s = np.linspace(max(0, -intercept), min(half, half-intercept), 80)
        view.curve("null", *null_map(s+intercept, s))
    for intercept in np.linspace(0, math.pi, 9)[1:-1]:
        s = np.linspace(max(0, intercept-half), min(half, intercept), 80)
        view.curve("null", *null_map(intercept-s, s))
    view.point("singular", (0, 0))
    view.point("singular", (0, math.pi))
    view.legend("null", "Light rays")
    view.legend("centre", "Circular axes")
    view.legend("boundary", "Big bang and big crunch")
    view.legend("singular", "Curvature singularities")
    moment = cf.slices.moments(METRIC, "turnaround", label=MOMENT)[0]
    lo, hi = moment.reach("cosmological", "s")
    s = np.linspace(lo, hi, 80)
    view.slice(moment, [null_map(np.full_like(s, math.pi/4), s)])
    view.set(settings=r"$b=L=1$, $p=\psi-s$, $q=\psi+s$.",
             restriction="Each interior point a 2-torus.")
    return [view]


def conformal_captions():
    return {(METRIC, "cosmological"):
            [r"The closed cosmological region between its big bang and big crunch, with both circular directions suppressed."]}


def curvature_views(ck, src, cv):
    """The stored four-dimensional scalars on the moment of maximum volume."""
    import numpy as np
    xaxis = cv.axis("$s$", 0, math.pi/2)
    scales = {"K": {"kind": "linear", "from": -135, "to": 80},
              "R": {"kind": "linear", "from": -0.2, "to": 14.5}}
    panels = [{"id": letter, "name": "$"+letter+(r"L^4$" if letter == "K" else r"L^2$"),
               "aspect": 2, "scale": scales[letter], "ticks": cv.value_ticks(scales[letter])}
              for letter in ("K", "R")]
    cv.shape_checks(ck, METRIC, [p["aspect"] for p in panels])
    curves = []
    for field, letter in (("kretschmann", "K"), ("ricci_scalar", "R")):
        stored = cv.Held(cv.Published(src, METRIC, "cosmological", field,
                                      ["s", "psi"], ["L", "b"]), {"psi": math.pi/4})
        curves += cv.line_curves(ck, METRIC, letter, letter, stored, "s", xaxis,
                                 scales[letter], cv.PLOT_W/2, PARAMS)
        s = np.linspace(0, math.pi/2, 101)
        c = np.cos(2*s)
        expected = 14*c*c if letter == "R" else 76*c**4+128*c*c-128
        ck(METRIC+" "+letter+" at maximum volume", np.max(np.abs(stored(s, 1, 1)-expected)) < 1e-9)
    return [{"id": "turnaround", "label": "At maximum volume", "system": "cosmological",
             "kind": "line", "x": xaxis, "panels": panels, "curves": curves, "marks": [],
             "legend": [["line", "K", "$K$"], ["line", "R", "$R$"]],
             "caption": [r"Spacetime curvature across the 3-sphere at maximum volume ($\psi=\pi/4$, $b=1$)."],
             "ricci_sentences": [], "settings": r"$L=1$.", "clip": "Linear vertical scales.",
             "per_scalar": cv.scalar_clips("Linear scale.", "Linear scale.")}]
