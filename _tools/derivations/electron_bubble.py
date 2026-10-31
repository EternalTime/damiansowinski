"""The rotating geometric shell, from gr-qc/0109085v2, (94) and (108).

The shell is outside every smooth chart's open domain.
The electromagnetic qualification and the Israel tensor are derived in
electron_bubble.md; these bulk tensors do not contain either surface layer.
"""

METRIC = "electron_bubble"
CHARTS = ("exterior", "interior", "exterior_kerr", "interior_kerr", "cartesian")
PARAMETERS = ["m", "q", "a", "r_0 = q^2/(2m)"]
SIGMA = r"\Sigma = r^2 + a^2\cos^2\theta"
B = "B = r^2 + a^2"
DELTA = r"\Delta = r^2 + a^2 - 2mr + q^2"
DIMENSIONS = {
    (METRIC, system): {
        **({"T": "T", "x": "L", "y": "L", "z": "L"} if system == "cartesian" else
           {("T" if system.endswith("kerr") else "t"): "T", "r": "L", r"\theta": "1",
            (r"\psi" if system.endswith("kerr") else r"\phi"): "1",
            r"\Sigma": "L**2", "B": "L**2",
            **({r"\Delta": "L**2"} if system == "exterior" else {})}),
        "m": "L", "q": "L", "a": "L", "r_0": "L",
    } for system in CHARTS
}


def line_element(system, clock):
    """The source's f=mr-q²/2 or f=0, with x^0=c times the printed time."""
    if system == "cartesian":
        return f"ds^2 = -{clock}dT^2 + dx^2 + dy^2 + dz^2"
    kerr = system.endswith("kerr")
    exterior = system.startswith("exterior")
    time, angle = ("T", r"\psi") if kerr else ("t", r"\phi")
    c = "c\\," if clock else ""
    metric = f"ds^2 = -{clock}d{time}^2"
    radial = "dr^2" if kerr else (r"\dfrac{\Sigma}{\Delta}dr^2" if exterior else r"\dfrac{\Sigma}{B}dr^2")
    metric += f" + {radial}" + r" + \Sigma d\theta^2 + B\sin^2\theta\,d" + angle + "^2"
    if kerr:
        metric += r" - 2a\sin^2\theta\,dr\,d\psi"
    if exterior:
        metric += r" + \dfrac{2mr - q^2}{\Sigma}\left(" + c + f"d{time}"
        metric += " + dr" if kerr else ""
        metric += r" - a\sin^2\theta\,d" + angle + r"\right)^2"
    return metric


def chart_spec(system):
    """Five source-defined charts and a numerical comparison with its matrix."""
    import verify_metrics as vm

    cartesian = system == "cartesian"
    kerr = system.endswith("kerr")
    exterior = system.startswith("exterior")
    time, angle = ("T", r"\psi") if kerr else ("t", r"\phi")
    coords = ["T", "x", "y", "z"] if cartesian else [time, "r", r"\theta", angle]
    parameters = PARAMETERS + ([] if cartesian else [SIGMA, B] + ([DELTA] if system == "exterior" else []))
    names = {"exterior": "Boyer-Lindquist Exterior", "interior": "Flat Oblate Interior",
             "exterior_kerr": "Kerr Angular Exterior", "interior_kerr": "Flat Kerr Angular Interior",
             "cartesian": "Flat Cartesian Interior"}
    if cartesian:
        domains = [f"{v} \\in (-\\infty, \\infty)" for v in coords]
        domains += [r"\dfrac{x^2 + y^2}{r_0^2 + a^2} + \dfrac{z^2}{r_0^2} < 1"]
    else:
        domains = [time + r" \in (-\infty, \infty)",
                   r"r \in (r_0, \infty)" if exterior else r"r \in (0, r_0)",
                   r"\theta \in [0, \pi]", angle + r" \in [0, 2\pi)"]
    probe = vm.Reader(coords, parameters, ())
    return {
        "metric_id": METRIC,
        "system": {"id": system, "name": names[system], "coords": coords,
                   "parameters": parameters, "domains": domains,
                   "line_element": line_element(system, "c^2")},
        "chart_line_element": line_element(system, ""),
        "printer": {"lead": list(probe.symbol.values()), "flip": False},
        "ricci_scalar": "0",
        "check": lambda chart: chart_check(chart, system),
    }


def chart_check(chart, system):
    """Source matrices at points, including both charge signs and a static control.

    No symbolic determinant or tensor recomputation is needed here.
    The angular-chart Jacobian is the Appendix C transformation.
    """
    import numpy as np

    reader = chart.reader
    exterior = system.startswith("exterior")
    for m, q, a, theta in ((1, 2, 3, .7), (2, -3, 1, 1.1), (1, 1.5, 0, .9)):
        r0 = q*q/(2*m)
        r = r0 * (1.4 if exterior else .6)
        at = {reader.parameters[k]: v for k, v in (("m", m), ("q", q), ("a", a))}
        if system == "cartesian":
            want = np.diag([-1., 1., 1., 1.])
        else:
            at.update({reader.symbol["r"]: r, reader.symbol[r"\theta"]: theta})
            sine = np.sin(theta)**2
            sigma, b = r*r + a*a*np.cos(theta)**2, r*r + a*a
            f2 = 2*m*r - q*q if exterior else 0
            want = np.diag([-1., sigma/(b-f2), sigma, b*sine])
            v = np.array([1., 0., 0., -a*sine])
            want += f2/sigma * np.outer(v, v)
            if system.endswith("kerr"):
                jac = np.eye(4)
                jac[0, 1], jac[3, 1] = -f2/(b-f2), -a/(b-f2)
                want = jac.T @ want @ jac
        got = np.array(chart.geo.g.subs(at).evalf(), dtype=float)
        if not np.allclose(got, want, atol=2e-12, rtol=2e-12):
            raise AssertionError(f"{METRIC}/{system}: source matrix or Appendix C pullback disagrees")


def register_charts(pc):
    pc.CHARTS[METRIC] = [lambda system=system: chart_spec(system) for system in CHARTS]


def diagram_specs(nr):
    """Axial null geodesics in each retained smooth coordinate region."""
    rows = []
    for system in CHARTS:
        cartesian = system == "cartesian"
        kerr = system.endswith("kerr")
        exterior = system.startswith("exterior")
        time, angle = ("T", "psi") if kerr else ("t", "phi")
        fixed = {"x": "0", "y": "0"} if cartesian else {"theta": "0", angle: "0"}
        box = (-1.4, 1.4, -1.4, 1.4) if cartesian else ((0, 4, -2, 2) if exterior else (0, 1.4, -.7, .7))
        coordinate = "z" if cartesian else "r"
        clock = "T" if cartesian else time
        where = "1-z**2" if cartesian else None
        rows.append(nr.Diagram(METRIC, system, "axis", "The rotation axis", (clock, coordinate), box,
                               f"${coordinate}/m$", f"$c{clock}/m$", {"m": 1, "q": "sqrt(2)", "a": 2},
                               fixed, tau=clock, where=where,
                               lines=(("shell", "r", "1", "the shell"),) if not cartesian else
                                     (("shell", "r", "1", "the shell"), ("shell", "r", "-1", None))))
    return rows


def diagram_captions():
    out = {}
    for system in CHARTS:
        time = "T" if system.endswith("kerr") or system == "cartesian" else "t"
        plane = r"$x=y=0$" if system == "cartesian" else r"$\theta=0$"
        region = "flat interior" if system.startswith("interior") or system == "cartesian" else "retained exterior"
        out[(METRIC, system, "axis")] = [
            f"The rotation axis ({plane}) in the {region}, with null geodesics in $c{time}$ and "
            + ("$z$" if system == "cartesian" else "$r$") + r"; the shell at $r_0=m$."
        ]
    return out


def null_forms():
    """Axial ray integrals in the exact parameter choice of the drawings."""
    import numpy as np
    forms = {}
    for system in CHARTS:
        if system == "exterior":
            plus = lambda t, r: t+r+np.log(((r-1)**2+5)/5)
            minus = lambda t, r: t-r-np.log(((r-1)**2+5)/5)
        elif system == "exterior_kerr":
            plus = lambda T, r: T+r
            minus = lambda T, r: T-r-2*np.log(((r-1)**2+5)/5)
        else:
            plus = lambda t, r: t+r
            minus = lambda t, r: t-r
        forms[(METRIC, system, "axis")] = (plus, minus, None)
    return forms


def axis_star(z):
    """Optical axial distance, joined continuously to flat z at |z|=r0=1."""
    import numpy as np
    z = np.asarray(z, dtype=float)
    radius = np.abs(z)
    return np.sign(z) * np.where(radius <= 1, radius, radius + np.log(((radius-1)**2+5)/5))


def embedding_views(ck, src, em):
    """The shell at one common time, embedded by rho=sqrt(5)sin(theta), z=cos(theta).

    Its induced spatial metric is the continuous limit from the flat oblate interior.
    Isometry is checked against those stored metric components, including both poles.
    """
    import math
    import numpy as np
    sl = em.Slice(src, METRIC, "interior", r"\theta", r"\phi", {"t": 0, "r": 1},
                  {"m": 1, "q": "sqrt(2)", "a": 2})
    piece = em.Piece("shell", "shell", sl, 0, math.pi, 1, -1,
                     (("axis", "the north pole"), ("axis", "the south pole")),
                     [(math.pi/2, "r", None)], 5)
    ck.isometry("electron bubble, the induced shell geometry", piece)
    ck.radius("electron bubble, the oblate equatorial radius", piece, lambda theta: np.sqrt(5)*np.sin(theta), 5)
    ck.form("electron bubble, the polar height", piece, np.cos, 5)
    surface = em.Surface([piece])
    fig = em.figure_of([surface], {"shell": "cover"}, 5)
    fig.legend("fill", "cover", "the shell")
    fig.legend("line", "r", "the equator")
    return [em.view("shell", "The oblate shell", "$m$", [surface], fig.done(),
                    settings=r"$m=1$, $q^2=2m^2$, $a=2m$.")]


def flat_slices(spec, sl):
    """The embedded shell's poles at the common instant t=T=0."""
    return [sl.Mark(m, points=[(0, -1), (0, 1)] if spec.system == "cartesian" else [(0, 1)])
            for m in sl.moments(METRIC)]


def conformal_views(ck, src, cf):
    """The totally geodesic joined axis, compactified by arctan(t +/- z_star).

    z_star=z inside; outside its derivative is B/Delta on each half-axis.
    The complete joined axis is a diamond with regular interior and two timelike
    shell crossings, no Kerr ring or second radial sheet.
    """
    import numpy as np
    params = {"m": 1, "q": "sqrt(2)", "a": 2}

    def mapped(t, z):
        return np.arctan((t-axis_star(z))/2), np.arctan((t+axis_star(z))/2)

    exterior = cf.Plane(src, METRIC, "exterior", ("t", "r"), {"theta": "0", "phi": "0"}, params)
    ck.chart("electron bubble exterior axis", exterior, mapped, ck.uniform(-8, 8), ck.uniform(1.01, 8),
             lambda t, r: (1, 0))
    inside = cf.Plane(src, METRIC, "cartesian", ("T", "z"), {"x": "0", "y": "0"}, params)
    ck.chart("electron bubble flat axis", inside, mapped, ck.uniform(-8, 8), ck.uniform(-.99, .99),
             lambda t, z: (1, 0))
    v = cf.View("axis", "The joined axis", [-cf.PI-.3, cf.PI+.3, -cf.PI-.3, cf.PI+.3])
    v.fill("region", cf.DIAMOND)
    v.fill("cover", cf.DIAMOND)
    cf.grid(v, "r", lambda z, t: mapped(t, z), (-4, -2, 0, 2, 4), cf.S_ALL)
    cf.grid(v, "surface", lambda z, t: mapped(t, z), (-1, 1), cf.S_ALL)
    cf.grid(v, "t", lambda t, s: (np.arctan((t-s)/2), np.arctan((t+s)/2)), (-4, 0, 4), cf.S_ALL)
    v.line("scri", [[[0, cf.PI], [cf.PI, 0], [0, -cf.PI], [-cf.PI, 0], [0, cf.PI]]])
    v.legend("cover", "the joined axis")
    v.legend("surface", "the shell at both poles")
    v.legend("scri", "null infinity")
    v.slice(cf.slices.moments(METRIC)[0], points=[mapped(0, -1), mapped(0, 1)])
    v.set(restriction="The rotation axis only, each point a single event.",
          settings=r"$m=1$, $q^2=2m^2$, $a=2m$.")
    return [v]


def curvature_views(ck, src, cv):
    """Stored exterior bulk K along the axis and equator; no shell scalar."""
    import numpy as np
    x = cv.axis("$r/m$", 1, 6)
    scale = {"kind": "symlog", "from": -100, "to": 100, "linear": .001}
    published = cv.nr_published(src, METRIC, "exterior", "kretschmann", ["r", "theta"], ["m", "q", "a"])
    curves = []
    for cls, theta in (("axis", 0), ("equator", np.pi/2)):
        held = cv.Held(published, {"theta": theta})
        curves += cv.line_curves(ck, METRIC, "K", cls, held, "r", x, scale, cv.PLOT_W,
                                 {"m": 1, "q": np.sqrt(2), "a": 2})
    return [{"id": "exterior", "label": "Smooth exterior", "kind": "line", "x": x,
             "panels": [{"id": "K", "name": "$Km^4$", "aspect": 1, "scale": scale,
                         "ticks": cv.value_ticks(scale)}], "curves": curves, "marks": [],
             "legend": [["line", "axis", "the rotation axis"], ["line", "equator", "the equator"]],
             "caption": [r"The smooth exterior Kretschmann scalar $K$ ($\theta=0$ and $\theta=\pi/2$), from its outer limit at the shell."],
             "settings": r"$m=1$, $q^2=2m^2$, $a=2m$.",
             "clip": "Signed logarithmic scale, linear near zero; both curves within the panel."}]
