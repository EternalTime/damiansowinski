"""Eleven dimensional wave: the two sourced charts and their static drawings.

Sources and normalization are recorded in kowalski_glikman.md.
"""

ID = "kowalski_glikman"
CHARTS = ("brinkmann", "rosen")

import itertools
import math

import numpy as np
import sympy as sp


def chart_spec(system):
    """Figueroa-O'Farrill and Papadopoulos (12), (14); Blau et al. (14)-(17)."""
    import chart_printer as cp
    import verify_metrics as vm

    rosen = system == "rosen"
    null, letter = ("w", "y") if rosen else ("v", "x")
    coords = ["u", null] + [f"{letter}_{i}" for i in range(1, 10)]
    blocks = [" + ".join(f"{letter}_{i}^2" for i in block) for block in (range(1, 4), range(4, 10))]
    differentials = [" + ".join(f"d{letter}_{i}^2" for i in block) for block in (range(1, 4), range(4, 10))]
    if rosen:
        line = (r"ds^2 = -2{c}\,du\,dw + \sin^2\left(\mu {c}u/3\right)\left(" + differentials[0]
                + r"\right) + \sin^2\left(\mu {c}u/6\right)\left(" + differentials[1] + r"\right)")
    else:
        line = (r"ds^2 = -\mu^2\left(\frac{" + blocks[0] + r"}{9} + \frac{" + blocks[1]
                + r"}{36}\right){c2}du^2 - 2{c}\,du\,dv + " + " + ".join(differentials))
    domains = [c + r" \in (-\infty,\infty)" for c in coords]
    if rosen:
        domains[0] = r"u \in (0,3\pi/(\mu c))"
    reader = vm.Reader(coords, [r"\mu"], ())
    return {"metric_id": ID, "system": {"id": system, "name": system.title(), "coords": coords,
            "domains": domains, "parameters": [r"\mu"],
            "line_element": line.replace("{c2}", "c^2").replace("{c}", "c")},
            "chart_line_element": line.replace("{c2}", "").replace("{c}", ""),
            "printer": {"lead": [reader.parameters["mu"], *reader.symbol.values()], "flip": False},
            "check": lambda chart: check_chart(chart, system), **({"time": "u"} if rosen else {})}


def check_chart(chart, system):
    """Source curvature and the Rosen pullback at three regular points, without symbolic reduction."""
    import verify_metrics as vm

    mu = chart.reader.parameters["mu"]
    x = chart.symbols
    points = [{**{s: sp.Rational((i + 2) * n, 13) for i, s in enumerate(x)},
               x[0]: sp.Rational(n + 1, 3), mu: sp.Rational(n + 2, 4)} for n in (1, 2, 3)]

    def near(value, label):
        value = sp.sympify(value)
        if value == 0:
            return
        for p in points:
            if abs(complex(value.subs(p).evalf(45))) > 1e-30:
                raise AssertionError(f"{ID}/{system}: {label}")

    if system == "brinkmann":
        # Source (5): all independent R_uiuj, including their vanishing off-diagonal entries.
        r = chart.geo.riemann_llll()
        C = chart.geo.weyl_llll()
        for a in range(11):
            for b in range(a + 1, 11):
                for c in range(11):
                    for d in range(c + 1, 11):
                        tidal = mu**2 / (9 if b < 5 else 36) if a == c == 0 and b == d and b >= 2 else 0
                        weyl = mu**2 / (18 if b < 5 else -36) if a == c == 0 and b == d and b >= 2 else 0
                        near(vm._at(r, (a, b, c, d)) - tidal, f"R_{a}{b}{c}{d}")
                        near(vm._at(C, (a, b, c, d)) - weyl, f"C_{a}{b}{c}{d}")
        ricci = chart.geo.ricci_ll()
        for a in range(11):
            for b in range(11):
                near(vm._at(ricci, (a, b)) - (mu**2/2 if a == b == 0 else 0), f"Ricci {a},{b}")
        return

    # The map is evaluated before the matrices are multiplied: no symbolic determinant or simplify.
    u, w, *ys = x
    omega = [mu / (3 if i < 3 else 6) for i in range(9)]
    scale = [sp.sin(o*u) for o in omega]
    image = sp.Matrix([u, w + sum(a*sp.diff(a, u)*y*y/2 for a, y in zip(scale, ys))]
                      + [a*y for a, y in zip(scale, ys)])
    J = image.jacobian(x)
    H = sum(o*o*X*X for o, X in zip(omega, image[2:]))
    gB = sp.diag(0, 0, *([1]*9))
    gB[0, 0], gB[0, 1], gB[1, 0] = -H, -1, -1
    for point in points:
        j = J.subs(point).evalf(45)
        pulled = j.T * gB.subs(point).evalf(45) * j
        target = chart.geo.g.subs(point).evalf(45)
        if max(abs(v) for v in pulled - target) > 1e-30:
            raise AssertionError(f"{ID}/{system}: Brinkmann metric pulled back")


def embedding_views(ck, src, e):
    """The induced transverse plane, pulled back to polar coordinates and embedded isometrically."""
    fixed = {"u": 3, "v": 0, **{f"x_{i}": 0 for i in range(2, 10) if i != 4}}
    sl = e.Slice(src, ID, "brinkmann", "x_1", None, fixed, {"mu": 1}, turn="x_4")
    piece = e.Piece("plane", "sheet", sl, 0, 3, size=6,
                    ends=(("centre", None), ("edge", r"the transverse plane extending to infinity")),
                    marks=[(r, "r", None) for r in (1, 2, 3)])
    ck.isometry("Eleven dimensional wave: transverse plane", piece)
    ck.form("Eleven dimensional wave: Z=0", piece, lambda r: np.zeros_like(r), 6)
    surface = e.Surface([piece], label=r"$\mu cu=3$", time=3)
    fig = e.figure_of([surface], {"sheet": "cover"}, 6)
    fig.legend("fill", "cover", "the transverse Euclidean plane")
    fig.legend("line", "r", r"$\mu\rho=1,2,3$")
    fig.legend("line", "meridian", r"$\varphi$ constant")
    return [e.view("plane", "Transverse plane", r"$1/\mu$", [surface], fig.done(),
                   settings=r"$\mu=1$, $\mu cu=3$, $v=0$; the remaining transverse coordinates zero.")]


EMBEDDING_CAPTION = [
    r"The transverse plane $(x_1,x_4)=(\rho\cos\varphi,\rho\sin\varphi)$ at fixed $u,v$, "
    r"with Euclidean embedding $(X,Y,Z)=(x_1,x_4,0)$ and induced metric $d\rho^2+\rho^2d\varphi^2$."
]


def ray_figure(spec, pj):
    """Full null geodesics in u,v,x1,x4, projected to the two transverse axes and cu.

    The numerical ODE uses only the published Christoffel symbols and metric.
    Its transverse positions are compared to independent oscillator solutions.
    """
    from types import SimpleNamespace
    from scipy.integrate import solve_ivp
    import null_rays as nr
    import slices

    _, entry, reader = nr.load(ID, spec.system)
    rosen = spec.system == "rosen"
    letter, null = ("y", "w") if rosen else ("x", "v")
    kept = ["u", null, f"{letter}_1", f"{letter}_4"]
    names = entry["coords"]
    symbols = [reader.symbol[c] for c in kept]
    fixed = {reader.c: 1, reader.parameters["mu"]: 1,
             **{s: 0 for c, s in reader.symbol.items() if c not in kept}}
    g = sp.zeros(4)
    for component in entry["metric_components"]:
        if all(c in kept for c in component["indices"]):
            a, b = [kept.index(c) for c in component["indices"]]
            g[a, b] = reader(component["value"]).subs(fixed)
    metric = sp.lambdify(symbols, g, "numpy")
    gamma = {}
    for component in entry["christoffel"]["variants"]["ull"]["nonzero"]:
        a, b, c = component["indices"]
        if b not in kept or c not in kept:
            continue
        value = reader(component["value"]).subs(fixed)
        if a not in kept:
            if value != 0:
                raise AssertionError("A transverse geodesic leaves its four coordinate submanifold")
            continue
        gamma[tuple(kept.index(d) for d in (a, b, c))] = sp.lambdify(symbols, value, "numpy")

    def rhs(_, state):
        acceleration = np.zeros(4)
        for (a, b, c), fun in gamma.items():
            acceleration[a] -= float(fun(*state[:4])) * state[4+b]*state[4+c]
        return np.r_[state[4:], acceleration]

    fig = pj.Figure(spec.view, spec.label, pj.Camera(-55, 25))
    lo, hi = (1.2, 3*math.pi-1.2) if rosen else (0.0, 6*math.pi)
    mid = 1.5*math.pi if rosen else 0.0
    frequencies = np.array([1/3, 1/6])
    for family in ("above", "below"):
        for angle in np.linspace(0, 2*math.pi, 8, endpoint=False):
            A = (6 if family == "above" else 2)*np.array([math.cos(angle), math.sin(angle)])
            def position(t):
                t = np.asarray(t)
                phase = t[..., None]*frequencies
                base = A*(np.sin(phase) if family == "above" else np.cos(phase))
                return base/np.sin(phase) if rosen else base
            if rosen:
                xy = A if family == "above" else A/np.tan(mid*frequencies)
                derivative = np.zeros(2) if family == "above" else -A*frequencies/np.sin(mid*frequencies)**2
                hidden = 0 if family == "above" else -float(np.sum(A*A*frequencies/np.tan(mid*frequencies)))/2
            else:
                xy = np.zeros(2) if family == "above" else A
                derivative = A*frequencies if family == "above" else np.zeros(2)
                hidden = 0
            point = np.r_[mid, hidden, xy]
            G = np.asarray(metric(*point), dtype=float)
            tangent = np.r_[1.0, 0.0, derivative]
            tangent[1] = -float(tangent@G@tangent)/(2*G[0, 1])
            pieces = []
            for endpoint in ([lo, hi] if rosen else [hi]):
                times = np.linspace(mid, endpoint, 501)
                sol = solve_ivp(rhs, (mid, endpoint), np.r_[point, tangent], t_eval=times,
                                rtol=2e-11, atol=2e-12, method="DOP853", max_step=0.08)
                if not sol.success:
                    raise AssertionError(sol.message)
                if np.max(np.abs(sol.y[2:4].T-position(times))) > 2e-8:
                    raise AssertionError("Transverse ray misses its sourced oscillator solution")
                for p, k in zip(sol.y[:4].T[::25], sol.y[4:].T[::25]):
                    G = np.asarray(metric(*p), dtype=float)
                    if abs(k@G@k) > 1e-8*max(1, np.linalg.norm(k)**2):
                        raise AssertionError("Projected ray misses the full null constraint")
                P = np.column_stack([sol.y[2], sol.y[3], sol.y[0]])
                pieces.append(P[::-1] if endpoint == lo else P)
            fig.line(family, np.vstack(pieces))
    fig.line("axis", np.array([[0, 0, lo], [0, 0, hi]]))
    # The embedded plane is u=3,v=0; its projection is the disc at that value of u.
    for m in slices.moments(ID):
        if lo < m.time < hi:
            theta = np.linspace(0, 2*math.pi, 181)
            a = np.array([math.sin(m.time/3), math.sin(m.time/6)]) if rosen else np.ones(2)
            ring = np.column_stack([3*np.cos(theta)/a[0], 3*np.sin(theta)/a[1], np.full(len(theta), m.time)])
            fig.slice(m, fills=[[ring]])
    fig.legend("line", "above", "rays from the Brinkmann origin")
    fig.legend("line", "below", "rays initially parallel to the axis")
    fig.legend("line", "axis", "the wave axis")
    fig.label(np.array([6, 0, mid]), f"$\\mu {letter}_1$", "l", cls="small")
    fig.label(np.array([0, 6, mid]), f"$\\mu {letter}_4$", "l", cls="small")
    fig.label(np.array([0, 0, hi]), r"$\mu cu$", "b", cls="small")
    return fig.done(pad=0.08), SimpleNamespace(entry=entry)


RAY_CAPTIONS = {
    "brinkmann": [r"Null geodesics in $(cu,x_1,x_4)$, with $v$ suppressed: transverse oscillations "
                  r"at $\mu/3$ and $\mu/6$, all other transverse coordinates zero."],
    "rosen": [r"Null geodesics in $(cu,y_1,y_4)$, with $w$ suppressed, inside "
              r"$0<\mu cu<3\pi$; the remaining transverse coordinates zero."]}
