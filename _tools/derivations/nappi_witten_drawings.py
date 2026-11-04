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
                   {"rho": "0", "lambda_": "0"}, tau="psi", families=nr.SIDEWAYS),
        nr.Diagram(METRIC, "bang_regular", "polar", "Cosmological region", ("x", "s"),
                   (-0.12, half+0.12, -0.10, 3.35), "$s$", "$x$", PARAMS,
                   {"rho": "0"}, tau="x", families=nr.SIDEWAYS, quotient="w", where="x",
                   cone="future cone of zero angular momentum"),
        nr.Diagram(METRIC, "invariant", "polar", "Cosmological region", ("r", "s"),
                   (-0.12, half+0.12, -0.12, 1.12), "$s$", "$r$", PARAMS,
                   {"rho": "0", "lambda_": "0"}, tau="r", families=nr.SIDEWAYS,
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
            [r"Null curves in the static extension, with periodic time $\lambda$ ($s=\pi/4$, $\rho=0$, $b=L=1$)."],
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
