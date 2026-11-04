"""Isometric slices and stored curvature of the reflected Curzon field."""
import numpy as np

import counter_rotating_discs as discs

METRIC = discs.METRIC


def embedding_views(ck, src, em):
    params = {"m": 1, "b": 3}
    plane = em.Slice(src, METRIC, "static", r"\rho", r"\phi", {"t": 0, "z": 0}, params,
                     functions=discs.curzon_functions())
    size = 16
    sheet = em.Piece("disc", "sheet", plane, 0, 6, 0, 1,
                     (("axis", r"$\rho=0$"), ("edge", "The disc extends to infinity.")),
                     [(r, "r", None) for r in (1, 2, 3, 4, 5, 6)], size)
    ck.isometry("counter rotating disc: intrinsic matter plane", sheet)
    surface = em.Surface([sheet])
    fig = em.figure_of([surface], {"sheet": "cover"}, size)
    fig.legend("fill", "cover", r"The matter disc at $t=0$")
    fig.legend("line", "r", r"$\rho/m=1,2,3,4,5,6$")
    fig.legend("line", "meridian", "Radial lines")
    disc = em.view("disc", "The matter disc", "$m$", [surface], fig.done(),
                   settings=r"Curzon seed, $m=1$, $b=3m$.")

    upper = em.Slice(src, METRIC, "static", "z", r"\phi", {"t": 0, "rho": 2}, params,
                     functions=discs.curzon_functions(1))
    lower = em.Slice(src, METRIC, "static", "z", r"\phi", {"t": 0, "rho": 2}, params,
                     functions=discs.curzon_functions(-1))
    hi = em.Piece("upper", "sheet", upper, 0, 2, 0, 1,
                  (("join", r"The seam at $z=0$."), ("edge", "The upper vacuum region extends onward.")),
                  [(0, "surface", None), (1, "r", None), (2, "r", None)], 8)
    lo = em.Piece("lower", "sheet", lower, -2, 0, -lower.rise(-2, 0), 1,
                  (("edge", "The lower vacuum region extends onward."), ("join", r"The seam at $z=0$.")),
                  [(-2, "r", None), (-1, "r", None)], 8)
    ck.isometry("counter rotating disc: upper cylinder", hi)
    ck.isometry("counter rotating disc: lower cylinder", lo)
    ck.add("counter rotating disc: both cylinder halves meet at the seam",
           abs(hi.at(0)[0]-lo.at(0)[0])+abs(hi.at(0)[1]-lo.at(0)[1]), 1e-10)
    surface = em.Surface([lo, hi])
    fig = em.figure_of([surface], {"sheet": "cover"}, 8)
    fig.legend("fill", "cover", "Two reflected vacuum halves")
    fig.legend("line", "surface", "Their circle on the matter seam")
    cylinder = em.view("cylinder", "The joined vacuum halves", "$m$", [surface], fig.done(),
                       settings=r"Curzon seed, $m=1$, $b=3m$, $\rho=2m$.")
    return [disc, cylinder]


def embedding_captions():
    return {
        (METRIC, "disc"): [r"Intrinsic spatial geometry of the matter disc in Euclidean space ($t=z=0$, Curzon seed)."],
        (METRIC, "cylinder"): [r"A cylinder through the junction, with the reflected vacuum halves meeting at the matter circle ($t=0$, $\rho=2m$)."],
    }


def diagram_specs(nr):
    reflected = {k: v.replace("(b+(1)*z)", "(b+Abs(z))") for k, v in discs.curzon_functions().items()}
    rows = [nr.Diagram(METRIC, "static", "disc", "Light in the matter plane", ("t", r"\rho"),
                       (0, 6, -3, 3), r"$\rho/m$", "$t/m$", {"m": 1, "b": 3}, {"z": "0", "phi": "0"},
                       functions=reflected, kretschmann=False, cones=(6, 6), step=0.004),
            nr.Diagram(METRIC, "static", "axis", "Light across the disc centre", ("t", "z"),
                       (-3, 3, -3, 3), "$z/m$", "$t/m$", {"m": 1, "b": 3}, {"rho": "0", "phi": "0"},
                       functions=reflected, kretschmann=False, families=nr.SIDEWAYS, cones=(6, 6), step=0.004,
                       lines=(("shell", "r", "0", "The matter seam"),))]
    # On the axis, Kerr's f is this rational function, g=1 and A=0.
    # All transverse first derivatives vanish there, so this restriction also gives the null geodesics.
    Z = "(b+Abs(z))"
    functions = {"f": f"({Z}**2+a**2-m**2)/(({Z}+m)**2+a**2)", "g": "1", "A": "0"}
    rows.append(nr.Diagram(METRIC, "stationary", "axis", "Light on the Kerr disc axis", ("t", "z"),
                           (-3, 3, -3, 3), "$z/m$", "$t/m$", {"m": 1, "a": 0.8, "b": 3},
                           {"rho": "0", "phi": "0"}, functions=functions, kretschmann=False,
                           families=nr.SIDEWAYS, cones=(6, 6), step=0.004,
                           input=r"Kerr axis restriction: $g=1$, $A=0$.",
                           lines=(("shell", "r", "0", "The Kerr matter seam"),)))
    return rows


def diagram_captions():
    return {
        (METRIC, "static", "disc"): [r"Radial light in the reflected Curzon disc ($z=0$, $b=3m$)."],
        (METRIC, "static", "axis"): [r"Light crossing the centre of the Curzon disc, where the two vacuum halves meet ($\rho=0$)."],
        (METRIC, "stationary", "axis"): [r"Light along the Kerr disc axis, crossing the matter seam at $z=0$ ($a=0.8m$, $b=3m$)."],
    }


def null_forms():
    from scipy.integrate import quad

    def primitive(system, plane, coordinate):
        def speed(x):
            if system == "stationary":
                Z = 3+abs(x)
                return ((Z+1)**2+0.64)/(Z*Z+0.64-1)
            if plane == "axis":
                return np.exp(2/(3+abs(x)))
            r = np.hypot(x, 3)
            return np.exp(2/r-x*x/(2*r**4))
        return np.vectorize(lambda x: quad(speed, 0, float(x), epsabs=1e-11, epsrel=1e-11)[0])(coordinate)

    return {(METRIC, system, plane):
            (lambda t, x, s=system, p=plane: t+primitive(s, p, x),
             lambda t, x, s=system, p=plane: t-primitive(s, p, x), None)
            for system, plane in (("static", "disc"), ("static", "axis"), ("stationary", "axis"))}


def slice_marks(sl, plane):
    disc = sl.moments(METRIC, "disc", label="$t=0$")[0]
    if plane == "axis":
        return [sl.Mark(disc, points=[(0, 0)])]
    cylinder = sl.moments(METRIC, "cylinder", label=r"$t=0$, $\rho=2m$")[0]
    return [sl.Mark(disc, sl.along(0, *disc.reach("static", r"\rho"))),
            sl.Mark(cylinder, points=[(0, 2)])]


def curvature_views(ck, src, cv):
    """Read the stored K; substitute the seed before fixing the vacuum plane z=m."""
    params = {"m": 1, "b": 3}
    K = cv.Published(src, METRIC, "static", "kretschmann", ["rho", "z"], list(params),
                     functions=discs.curzon_functions())
    held = cv.Held(K, {"z": 1})
    axis = cv.axis(r"$\rho/m$", 0.05, 6)
    scale = {"kind": "log", "from": 1e-5, "to": 0.01}
    panels = [{"id": "K", "name": "$K m^4$", "aspect": 1, "scale": scale, "ticks": cv.value_ticks(scale)}]
    cv.shape_checks(ck, METRIC, [1])
    curves = cv.line_curves(ck, METRIC, "K", "K", held, "rho", axis, scale, cv.PLOT_W, params)
    points = np.linspace(0.05, 6, 61)
    ck("counter rotating discs: finite positive vacuum curvature above the junction",
       np.all(np.isfinite(held(points, 1, 3))) and np.all(held(points, 1, 3)>0))
    return [{"id": "vacuum", "label": "Curvature above the disc", "system": "static", "kind": "line",
             "x": axis, "panels": panels, "curves": curves, "marks": [], "legend": [["line", "K", "$K$"]],
             "caption": [r"The Kretschmann scalar $K$ above the matter seam ($z=m$, Curzon seed)."],
             "settings": r"$m=1$, $b=3m$.", "clip": "Logarithmic vertical scale."}]
