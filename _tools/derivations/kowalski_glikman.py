"""Eleven dimensional wave: the two sourced charts and their static drawings.

Sources and normalization are recorded in kowalski_glikman.md.
"""

ID = "kowalski_glikman"
CHARTS = ("brinkmann", "rosen")

import math

import numpy as np
import sympy as sp


def chart_spec(system):
    """Figueroa-O'Farrill and Papadopoulos (12), (14); Blau et al. (14)-(17)."""
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


# The embedding diagram, in units of 1/mu: the wave fronts of the plane of x_1, one of the three
# directions of the flux, and x_4, one of the other six, from u = 0 to mu c u = 6 pi, a frame
# every pi/10 of mu c u, with the moments where the ring is a circle or a segment named.
TOP = 1.5                           # the proper radius of each wave front's disc
TENTHS = 60                         # mu c u = 6 pi, in tenths of pi
KEYS = (0, 15, 30, 45, 60)          # mu c u = 0, 3 pi/2, 3 pi, 9 pi/2 and 6 pi, in tenths of pi
LIFT = 1 / 6                        # the tube's height for each 1/mu of cu
REST = {f"x_{i}": 0 for i in (2, 3, 5, 6, 7, 8, 9)}


def tenths(k):
    """mu c u = k pi/10, the time of the frame k, to the nine decimals every frame is drawn at."""
    return round(k * math.pi / 10, 9)


def pi_tex(k):
    """k pi/10 as TeX, in lowest terms."""
    from fractions import Fraction

    f = Fraction(k, 10)
    if f == 0:
        return "0"
    return ("" if f.numerator == 1 else str(f.numerator)) + r"\pi" + ("" if f.denominator == 1 else f"/{f.denominator}")


def embedding_views(ck, src, e):
    """The wave front at sixty-one values of u of the Brinkmann chart, mu c u = 0 to 6 pi in units of
    1/mu, the plane of x_1 and x_4 with v and every other direction across the wave at zero: a surface
    of constant u has the metric dx_1^2 + dx_4^2 whatever v is on it, flat. The ring is 360 free
    particles at rest at u = 0 on x_1^2 + x_4^2 = 1/mu^2, each run with the published Christoffel
    symbols, u being affine on every geodesic since no Gamma^u is published: x_i'' = -Gamma^x_i_uu.
    The wave pulls them toward the axis at two rates, x_1'' = -mu^2 x_1/9 and x_4'' = -mu^2 x_4/36, so
    the ring is the ellipse of semi-axes |cos(mu c u/3)|/mu along x_1 and |cos(mu c u/6)|/mu along x_4:
    a segment of the x_4 axis of half length 1/(sqrt2 mu) at mu c u = 3 pi/2 and 9 pi/2, a segment of
    the x_1 axis of half length 1/mu at 3 pi, and back at rest on its first circle at 6 pi, every
    particle at its mirror image across the x_1 axis. Drawn first as the ring's world tube, each
    ellipse at the height of its u, and beside it as a movie of the flat wave fronts."""
    from scipy.integrate import solve_ivp
    from scipy.optimize import brentq

    name = "Kowalski-Glikman wave"
    gamma, R = e.published_christoffel(src, ID, "brinkmann")
    ck.exact(f"{name}: no published Gamma^u, so u is an affine parameter", not any(ix[0] == "u" for ix in gamma))
    across = [f"x_{i}" for i in range(1, 10)]
    ck.exact(f"{name}: the only published Gamma^x_i are Gamma^x_i_uu",
             {ix for ix in gamma if ix[0] in across} == {(c, "u", "u") for c in across})
    names = {R._plain(n): s for n, s in R.symbol.items()}
    mu = R.parameters["mu"]
    held = {names[c]: 0 for c in REST}
    ck.exact(f"{name}: a particle at rest with every x_i but x_1 and x_4 zero stays in their plane",
             all(gamma[(c, "u", "u")].subs(held) == 0 for c in REST))
    args = (names["x_1"], names["x_4"])
    g1 = sp.lambdify(args, gamma[("x_1", "u", "u")].subs(held).subs(mu, 1), "numpy")
    g4 = sp.lambdify(args, gamma[("x_4", "u", "u")].subs(held).subs(mu, 1), "numpy")
    alpha = np.linspace(0, 2 * math.pi, 361)[:-1]
    n = len(alpha)

    def rhs(u, w):
        x, y = w[:n], w[n:2 * n]
        return np.concatenate([w[2 * n:3 * n], w[3 * n:], -g1(x, y) * np.ones(n), -g4(x, y) * np.ones(n)])
    start = np.concatenate([np.cos(alpha), np.sin(alpha), np.zeros(2 * n)])
    end = tenths(TENTHS)
    run = solve_ivp(rhs, (0, end + 0.2), start, rtol=1e-12, atol=1e-13, dense_output=True, max_step=0.01)

    def form(u):
        """The ring at u in closed form, each particle where it set out times cos(mu c u/3) along x_1
        and cos(mu c u/6) along x_4, which the run is held to and never drawn from."""
        return np.column_stack([math.cos(u / 3) * np.cos(alpha), math.cos(u / 6) * np.sin(alpha)])
    us = np.linspace(0, end, 241)
    ck.add(f"{name}: every particle moves as cos(mu c u/3) along x_1 and cos(mu c u/6) along x_4",
           max(float(np.max(np.abs(np.column_stack([run.sol(u)[:n], run.sol(u)[n:2 * n]]) - form(u)))) for u in us), 1e-8)
    along_1 = lambda u: float(run.sol(u)[0])        # the particle that sets out on the x_1 axis
    along_4 = lambda u: float(run.sol(u)[n + 90])   # and the one on the x_4 axis
    for lo, hi, want, axis_, f in ((4.0, 5.5, 1.5, "x_4", along_1), (13.5, 15.0, 4.5, "x_4", along_1),
                                   (8.5, 10.5, 3.0, "x_1", along_4)):
        ck.add(f"{name}: the ring closes onto the {axis_} axis at mu c u = {want:g} pi",
               abs(brentq(f, lo, hi, xtol=1e-13) - want * math.pi), 1e-8)
    for k, half, axis_ in ((15, 1 / math.sqrt(2), "x_4"), (30, 1.0, "x_1"), (45, 1 / math.sqrt(2), "x_4")):
        w = run.sol(tenths(k))
        on, off = (w[n:2 * n], w[:n]) if axis_ == "x_4" else (w[:n], w[n:2 * n])
        ck.add(f"{name}, mu c u = {pi_tex(k)}: every particle is on the {axis_} axis, within {half:.6f}/mu of the centre",
               float(max(np.max(np.abs(off)), abs(np.max(np.abs(on)) - half))), 1e-8)
    w = run.sol(end)
    ck.add(f"{name}: the ring is back at rest on its first circle at mu c u = 6 pi, every particle at its mirror "
           "image across the x_1 axis",
           float(np.max(np.abs(w[:n] - np.cos(alpha))) + np.max(np.abs(w[n:2 * n] + np.sin(alpha)))
                 + np.max(np.abs(w[2 * n:]))), 1e-7)

    def ring_at(u):
        """The ring at u as the curve (a cos alpha, b sin alpha), a and b its semi-axes, through the
        particles run: past a focus along x_1 the particle that set out at alpha stands at pi - alpha
        of the ellipse, and past one along x_4 at -alpha, so each is put at its angle there."""
        w = run.sol(u)
        sa, sb = (1.0 if w[0] >= 0 else -1.0), (1.0 if w[n + 90] >= 0 else -1.0)
        at = np.round(np.degrees(np.arctan2(sb * np.sin(alpha), sa * np.cos(alpha)))).astype(int) % n
        if len(set(at.tolist())) != n:
            raise AssertionError(f"{name}: the particles at mu c u = {u} are not one to each angle")
        P = np.empty((n, 2))
        P[at] = np.column_stack([w[:n], w[n:2 * n]])
        return P, abs(float(w[0])), abs(float(w[n + 90]))

    frames = []
    for k in range(TENTHS + 1):
        u = tenths(k)
        sl = e.FlatPlane(src, ID, "brinkmann", "x_1", "x_4", {"u": repr(u), "v": 0, **REST}, params={"mu": 1})
        where = f"{name}, mu c u = {pi_tex(k)}"
        plane = e.disc(sl, "plane", TOP, "the centre of the ring", "the wave front runs on, flat, to infinity", [], 2 * TOP)
        ck.isometry(where, plane)
        P, a, b = ring_at(u)
        ck.add(f"{where}: the ring is the ellipse of semi-axes |cos(mu c u/3)| and |cos(mu c u/6)|",
               float(np.max(np.abs(P - np.column_stack([abs(math.cos(u / 3)) * np.cos(alpha),
                                                        abs(math.cos(u / 6)) * np.sin(alpha)])))), 1e-8)
        ring, dots = e.particles(ck, where, sl, plane, np.append(P[:, 0], P[0, 0]), np.append(P[:, 1], P[0, 1]))
        frames.append(e.Surface([plane], label=f"$\\mu cu = {pi_tex(k)}$", time=u, curves=[ring], dots=dots))
    surfaces = [frames[k] for k in KEYS]

    def rows(u):
        w = run.sol(np.atleast_1d(u))
        return np.abs(w[0]), np.abs(w[n + 90])

    def on_row(u, a, b):
        # Every particle of the ring, run with the published Christoffel symbols, on the row, each at
        # its angle of the ellipse.
        P, _, _ = ring_at(u)
        return float(np.max(np.abs(P - np.column_stack([a * np.cos(alpha), b * np.sin(alpha)]))))
    tube = e.stack(ck, name, surfaces, rows, LIFT, on_row, 2 * TOP,
                   "the world tube runs on before $u = 0$ and after $\\mu cu = 6\\pi$, repeating every $6\\pi$ of "
                   "$\\mu cu$")
    tube_fig = e.stack_figure(tube, 2 * TOP, "$\\mu cu$", [
        ("fill", "cover", "the ring at each $u$"),
        ("line", "particles", "the ring at five moments"),
        ("line", "worldline", "twelve particles' world lines"),
        ("line", "axis", "the axis of $u$")])
    fig, played = e.ring_movie(frames, 2 * TOP, "$\\mu cu$")
    fig.legend("fill", "cover", "the wave front, flat")
    fig.legend("line", "particles", "the ring, twelve particles marked")
    fig.legend("line", "meridian", "lines from the centre every $30°$")
    settings = "$\\mu = 1$, every $x_i$ but $x_1$ and $x_4$ zero."
    return [e.view("tube", "The ring's world tube", "$1/\\mu$", [tube], tube_fig, settings=settings,
                   height="$\\mu cu$, a height of $1/(6\\mu)$ for each $1/\\mu$ of $cu$"),
            e.view("ring", "A ring of particles", "$1/\\mu$", surfaces, fig.done(), movie=played, settings=settings)]


EMBEDDING_CAPTIONS = {
    (ID, "tube"): [
        r"The world tube of a ring of free particles ($0 \le \mu cu \le 6\pi$), each ellipse at the height of its "
        r"$u$, of semiaxes $|\cos(\mu cu/3)|/\mu$ along $x_1$ and $|\cos(\mu cu/6)|/\mu$ along $x_4$, closing onto a "
        r"segment at $\mu cu = 3\pi/2$, $3\pi$, and $9\pi/2$."],
    (ID, "ring"): [
        r"A ring of free particles at rest on $x_1^2 + x_4^2 = 1/\mu^2$ at $u = 0$, an ellipse of semiaxes "
        r"$|\cos(\mu cu/3)|/\mu$ along $x_1$, a direction of the flux, and $|\cos(\mu cu/6)|/\mu$ along "
        r"$x_4$, back on its circle at $\mu cu = 6\pi$."],
}


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
    # Each wave front of the embedding diagram is the plane of constant u, drawn out to the proper
    # radius TOP from the axis, in the ring's view and in the tube's alike: the disc of that radius at
    # the height of its u, and in the Rosen chart, where x_i = sin(omega_i cu) y_i, the ellipse of
    # semi-axes TOP/sin(mu cu/3) and TOP/sin(mu cu/6), for the fronts within its window, a front's time
    # being written to six decimals.
    for m in slices.moments(ID):
        if lo - 1e-6 <= m.time <= hi + 1e-6:
            theta = np.linspace(0, 2*math.pi, 181)
            a = np.array([math.sin(m.time/3), math.sin(m.time/6)]) if rosen else np.ones(2)
            ring = np.column_stack([TOP*np.cos(theta)/a[0], TOP*np.sin(theta)/a[1], np.full(len(theta), m.time)])
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


def curvature_views(ck, src, cv):
    """Only the two stored scalar contractions; no new curvature is calculated."""
    horizontal = cv.axis(r"$\mu cu$", 0.0, 6*math.pi)
    panels, curves = [], []
    for field, pid, power in (("kretschmann", "K", 4), ("ricci_scalar", "R", 2)):
        published = cv.nr_published(src, ID, "brinkmann", field, ["u"], ["mu"])
        scale = {"kind": "linear", "from": -1.0, "to": 1.0}
        panels.append({"id": pid, "name": f"${pid}/\\mu^{power}$", "aspect": 2.0,
                       "scale": scale, "ticks": cv.value_ticks(scale)})
        curves.extend(cv.line_curves(ck, ID, pid, pid, published, "u", horizontal, scale,
                                     cv.PLOT_W/2, {"mu": 1.0}))
    cv.shape_checks(ck, ID, [2.0, 2.0])
    return [{"id": "null_scalars", "label": "Null curvature contractions", "kind": "line", "x": horizontal,
             "panels": panels, "curves": curves, "marks": [],
             "legend": [["line", "K", "Kretschmann scalar"], ["line", "R", "Ricci scalar"]],
             "caption": [r"The Kretschmann scalar $K$ and the Ricci scalar $R$ of the eleven dimensional wave.",
                         r"$K=0$, with tidal eigenvalues $\mu^2/9$ and $\mu^2/36$.",
                         r"The vanishing Ricci scalar, $R=0$."],
             "ricci_sentences": [[2, 0]],
             "settings": r"$\mu=1$; length unit $1/\mu$.",
             "per_scalar": cv.scalar_clips("Linear scale; no clipping.", "Linear scale; no clipping."),
             "clip": "Linear scales; both curves within their panels."}]
