"""Gravitational waves in Friedmann universes filled with a stiff fluid, in units G = c = 1.

Bicak and Griffiths (1996), arXiv gr-qc/9610041, (1.1), (2.4), (2.5), (3.15), (3.16), (4.7),
(5.2) and (5.4)-(5.8), give the metric form and the five stiff backgrounds; their (7.1)-(7.4)
give the plane wave in the spatially flat universe. Alekseev and Griffiths (1995), Phys. Rev.
D 52, 4497, (2.12)-(2.13) and (4.1)-(4.4), give the profile H_2 and the closed form of M for
one term. Every line element here is theirs multiplied by -1, for the signature (-,+,+,+).
See waves_in_friedmann.md beside this file for the derivation and the choices.
"""

METRIC = "waves_in_friedmann"
CHARTS = ("plane_wave", "flat_planar", "flat_cylindrical", "open_planar", "open_cylindrical",
          "closed_toroidal")
NAMES = {"plane_wave": "Plane Wave, Flat Universe", "flat_planar": "Flat Universe, Planar",
         "flat_cylindrical": "Flat Universe, Cylindrical", "open_planar": "Open Universe, Planar",
         "open_cylindrical": "Open Universe, Cylindrical", "closed_toroidal": "Closed Universe, Toroidal"}
COORDS = {"plane_wave": ["\\eta", "z", "x", "y"], "flat_planar": ["\\eta", "z", "x", "y"],
          "flat_cylindrical": ["\\eta", "r", "\\phi", "z"], "open_planar": ["\\eta", "\\mu", "x", "y"],
          "open_cylindrical": ["\\eta", "r", "\\phi", "y"],
          "closed_toroidal": ["\\eta", "\\zeta", "\\delta", "\\sigma"]}

# The profile H = H_2(z/eta) of Alekseev and Griffiths's (4.3), normalised to H_2(-1) = 1, its
# slope K = dH_2/d(z/eta), the wave V and the function S of e^{-M} = 2 eta e^{-S}.
H_DEF = r"H = \dfrac{\left(\eta^2+2z^2\right)\arccos(z/\eta)-3z\sqrt{\eta^2-z^2}}{3\pi\eta^2}"
K_DEF = r"K = \dfrac{4\left(z\arccos(z/\eta)-\sqrt{\eta^2-z^2}\right)}{3\pi\eta}"
V_DEF = r"V = A\eta^2H"
S_DEF = r"S = -\dfrac{A^2\eta^2}{8}\left(4\eta^2H^2+\left(\eta^2-z^2\right)K^2\right)"
WAVE_PARAMETERS = ["\\gamma", "A", H_DEF, K_DEF, V_DEF, S_DEF]

DIMENSIONS = {(METRIC, s): {**{x: "1" for x in COORDS[s]}, "\\gamma": "L**2",
                            **({"A": "1", "H": "1", "K": "1", "V": "1", "S": "1"} if s == "plane_wave" else {})}
              for s in CHARTS}
HELD = {(METRIC, "plane_wave"): ("H", "K", "V", "S")}
# First derivatives in the names: H' = K and, from the Euler-Poisson-Darboux equation of
# Alekseev and Griffiths's (4.2) with n = 2, (1-w^2)K' + 3wK - 4H = 0 for w = z/eta; then
# V_eta = A(2 eta H - zK), V_z = A eta K and Bicak and Griffiths's (6.7) for S.
RATES = {(METRIC, "plane_wave"): {
    "H": {"\\eta": r"-zK/\eta^2", "z": r"K/\eta"},
    "K": {"\\eta": r"-\dfrac{z\left(4\eta H-3zK\right)}{\eta\left(\eta^2-z^2\right)}",
          "z": r"\dfrac{4\eta H-3zK}{\eta^2-z^2}"},
    "V": {"\\eta": r"A\left(2\eta H-zK\right)", "z": r"A\eta K"},
    "S": {"\\eta": r"-\dfrac{A^2\eta}{2}\left(\left(2\eta H-zK\right)^2+\eta^2K^2\right)",
          "z": r"-A^2\eta^2K\left(2\eta H-zK\right)"},
}}

LINES = {
    "plane_wave": r"ds^2 = 2\gamma\eta\left(e^{-S}\left(-d\eta^2+dz^2\right)+e^{V}dx^2+e^{-V}dy^2\right)",
    "flat_planar": r"ds^2 = 2\gamma\eta\left(-d\eta^2+dz^2+dx^2+dy^2\right)",
    "flat_cylindrical": r"ds^2 = 2\gamma\eta\left(-d\eta^2+dr^2+r^2d\phi^2+dz^2\right)",
    "open_planar": r"ds^2 = \gamma\sinh(2\eta)\left(-d\eta^2+d\mu^2+e^{-2\mu}\left(dx^2+dy^2\right)\right)",
    "open_cylindrical": r"ds^2 = \gamma\sinh(2\eta)\left(-d\eta^2+dr^2+\sinh^2r\,d\phi^2+\cosh^2r\,dy^2\right)",
    "closed_toroidal": r"ds^2 = \gamma\sin(2\eta)\left(-d\eta^2+d\zeta^2+\sin^2\zeta\,d\delta^2+\cos^2\zeta\,d\sigma^2\right)",
}
REALS = r" \in (-\infty,\infty)"
DOMAINS = {
    "plane_wave": [r"\eta \in (0,\infty)", r"z \in (-\eta,\eta)", "x"+REALS, "y"+REALS,
                   r"z=\eta\;\text{(wave front)}", r"\eta\to0\;\text{(big bang)}"],
    "flat_planar": [r"\eta \in (0,\infty)", "z"+REALS, "x"+REALS, "y"+REALS,
                    r"\eta\to0\;\text{(big bang)}"],
    "flat_cylindrical": [r"\eta \in (0,\infty)", r"r \in [0,\infty)", r"\phi \in [0,2\pi)", "z"+REALS,
                         r"\eta\to0\;\text{(big bang)}"],
    "open_planar": [r"\eta \in (0,\infty)", r"\mu"+REALS, "x"+REALS, "y"+REALS,
                    r"\eta\to0\;\text{(big bang)}"],
    "open_cylindrical": [r"\eta \in (0,\infty)", r"r \in [0,\infty)", r"\phi \in [0,2\pi)", "y"+REALS,
                         r"\eta\to0\;\text{(big bang)}"],
    "closed_toroidal": [r"\eta \in (0,\pi/2)", r"\zeta \in [0,\pi/2]", r"\delta \in [0,2\pi)",
                        r"\sigma \in [0,2\pi)", r"\eta\to0\;\text{(big bang)}", r"\eta\to\pi/2\;\text{(big crunch)}"],
}
# The scale factor squared, a^2 = gamma f(eta), of Bicak and Griffiths's (5.2), and the fluid
# potential of their (5.11)-(5.13), with Ricci = (1/2) d(sigma) d(sigma): their (5.9) reads
# 16 pi rho for what the Friedmann normalisation of their (5.1) and BG94 (3.1)-(3.4) makes 32 pi rho.
CURVATURE = {"plane_wave": "flat", "flat_planar": "flat", "flat_cylindrical": "flat",
             "open_planar": "open", "open_cylindrical": "open", "closed_toroidal": "closed"}


def chart_spec(system):
    """One chart of Bicak and Griffiths's stiff universes, or the plane wave behind its front."""
    import verify_metrics as vm

    key = (METRIC, system)
    params = WAVE_PARAMETERS if system == "plane_wave" else ["\\gamma"]
    probe = vm.Reader(COORDS[system], params, (), held=HELD.get(key, ()), rates=RATES.get(key))
    spec = {"metric_id": METRIC,
            "system": {"id": system, "name": NAMES[system], "coords": COORDS[system],
                       "domains": DOMAINS[system], "parameters": params, "line_element": LINES[system]},
            "chart_line_element": LINES[system], "printer": {},
            "check": lambda chart: chart_check(chart, system)}
    if system == "plane_wave":
        spec["reduce"] = lambda value: vm.norm(probe.by_rates(value))
    elif CURVATURE[system] == "open":
        import chart_printer as cp
        spec["pretty"] = cp.hyperbolic(probe.symbol["\\eta"])
        spec["ricci_scalar"] = r"-\dfrac{6}{\gamma\sinh^3(2\eta)}"
        spec["kretschmann"] = r"\dfrac{60}{\gamma^2\sinh^6(2\eta)}"
    elif CURVATURE[system] == "closed":
        spec["pretty"] = trigonometric(probe.symbol["\\eta"])
        spec["ricci_scalar"] = r"-\dfrac{6}{\gamma\sin^3(2\eta)}"
        spec["kretschmann"] = r"\dfrac{60}{\gamma^2\sin^6(2\eta)}"
    return spec


def trigonometric(x):
    """The closed universe's values, which the checker hands back in sin x and cos x, rewritten
    as chart_printer.hyperbolic rewrites sinh and cosh: reduced by cos^2 = 1 - sin^2, the
    denominator cleared of cos by its conjugate, factored, and 1 - sin^2 written cos^2."""
    import sympy as sp

    S, C = sp.symbols("_S _C", positive=True)
    relation = [C ** 2 - 1 + S ** 2]

    def reduce(p):
        return sp.expand(sp.reduced(sp.expand(p), relation, C, S)[1])

    def pretty(value):
        v = sp.sympify(value)
        v = v.subs({sp.tan(x): sp.sin(x)/sp.cos(x), sp.cot(x): sp.cos(x)/sp.sin(x),
                    sp.sec(x): 1/sp.cos(x), sp.csc(x): 1/sp.sin(x)}).subs({sp.sin(x): S, sp.cos(x): C})
        num, den = (reduce(q) for q in sp.fraction(sp.together(v)))
        d0, d1 = sp.Poly(den, C).coeff_monomial(1), sp.Poly(den, C).coeff_monomial(C)
        if d1 != 0:
            num, den = reduce(num * (d0 - d1 * C)), sp.expand(d0 ** 2 - d1 ** 2 * (1 - S ** 2))
        out = sp.factor(sp.cancel(sp.factor(num) / sp.factor(den)))
        # (S - 1)^k (S + 1)^k, which factor makes of 1 - S^2, is read as (-C^2)^k and cancelled again.
        powers = sp.Mul.make_args(out)
        pair = {}
        for f in powers:
            base, k = (f.base, f.exp) if f.is_Pow else (f, 1)
            if sp.expand(base - (S - 1)) == 0 or sp.expand(base - (S + 1)) == 0:
                pair[sp.expand(base)] = k
        if len(pair) == 2 and len(set(pair.values())) == 1:
            k = next(iter(pair.values()))
            rest = [f for f in powers if sp.expand(f.base if f.is_Pow else f) not in pair]
            out = sp.factor(sp.Mul(*rest) * (-C ** 2) ** k)
        return out.subs({S: sp.sin(x), C: sp.cos(x)})
    return pretty


def chart_check(chart, system):
    """Every Ricci slot against the stiff fluid, and the held slopes, at points; no simplify."""
    import sympy as sp

    R = chart.reader
    X = chart.symbols
    eta = X[0]
    kind = CURVATURE[system]
    sigma = sp.sqrt(3)*(sp.log(eta) if kind == "flat" else sp.log(sp.tanh(eta)) if kind == "open"
                        else sp.log(sp.tan(eta)))
    gradient = [sp.diff(sigma, c) for c in X]
    ricci = chart.geo.ricci_ll()
    gamma = R.parameters["gamma"]
    points = ((sp.Rational(7, 10), sp.Rational(1, 5)), (sp.Rational(11, 10), sp.Rational(-1, 2)),
              (sp.Rational(1, 2), sp.Rational(1, 3)))
    for k, (ev, sv) in enumerate(points):
        at = {eta: ev, X[1]: sv if system != "closed_toroidal" else sp.Rational(1, 2)+sv/3, X[2]: sp.Rational(1, 3),
              X[3]: sp.Rational(1, 4), gamma: (1, sp.Rational(5, 2), sp.Rational(7, 10))[k]}
        if system == "plane_wave":
            at[R.parameters["A"]] = (sp.Rational(13, 10), sp.Rational(3, 5), 2)[k]

        def number(value):
            return sp.N(sp.sympify(value).subs(R.held).doit().subs(at), 40)
        for i in range(4):
            for j in range(4):
                gap = number(ricci[i][j]) - number(gradient[i]*gradient[j]/2)
                if abs(gap) > 1e-25:
                    raise AssertionError(f"{METRIC}/{system}: Ricci slot {i},{j} misses the stiff fluid by {gap}")
        for name, rates in RATES.get((METRIC, system), {}).items():
            for coord, rate in rates.items():
                gap = number(R(rate)) - number(sp.diff(R.held[R.parameters[name]], R.symbol[coord]))
                if abs(gap) > 1e-25:
                    raise AssertionError(f"{METRIC}/{system}: held slope of {name} along {coord} is off by {gap}")


# ------------------------------------------------------------------ drawings
# Every chart's plane of eta and its spatial coordinate carries a conformal factor times
# -d eta^2 + dx^2, so its light rays are the lines eta -+ x constant, and the wave fronts of
# Bicak and Griffiths's (7.1), (8.1), (9.1), (10.1) and (11.x) are the rays x = eta.
DRAW = {"\\gamma": 1, "A": 2}
PLANES = (
    ("plane_wave", "wave", "Behind the front", ("\\eta", "z"), (-2.0, 2.0, 0.0, 3.0), {"x": "0", "y": "0"},
     "eta**2 - z**2", ("front", "eta", "wave front $z=\\eta$", "signed")),
    ("flat_planar", "flat", "Flat universe", ("\\eta", "z"), (-2.0, 2.0, 0.0, 3.0), {"x": "0", "y": "0"},
     None, ("front", "eta", "wave front $z=\\eta$", "signed")),
    ("flat_cylindrical", "flat_cylinder", "Cylinders, flat", ("\\eta", "r"), (0.0, 3.0, 0.0, 3.0),
     {"phi": "0", "z": "0"}, None, ("front", "eta", "cylindrical front $r=\\eta$")),
    ("open_planar", "open", "Open universe", ("\\eta", "\\mu"), (-2.0, 2.0, 0.0, 3.0), {"x": "0", "y": "0"},
     None, ("front", "eta", "wave front $\\mu=\\eta$", "signed")),
    ("open_cylindrical", "open_cylinder", "Cylinders, open", ("\\eta", "r"), (0.0, 3.0, 0.0, 3.0),
     {"phi": "0", "y": "0"}, None, ("front", "eta", "cylindrical front $r=\\eta$")),
    ("closed_toroidal", "closed", "Closed universe", ("\\eta", "\\zeta"), (0.0, 1.5707963267948966, 0.0, 1.5707963267948966),
     {"delta": "0", "sigma": "0"}, None, ("front", "eta", "toroidal front $\\zeta=\\eta$")),
)
CAPTION = {
    "wave": r"The plane of $\eta$ and $z$ ($x=y=0$) behind the wave front $z=\eta$, light rays at $45°$.",
    "flat": r"The plane of $\eta$ and $z$ ($x=y=0$) of the flat stiff universe ahead of the front $z=\eta$.",
    "flat_cylinder": r"The plane of $\eta$ and $r$ ($\phi=z=0$) of the flat stiff universe, with the cylindrical front $r=\eta$.",
    "open": r"The plane of $\eta$ and $\mu$ ($x=y=0$) of the open stiff universe, with the plane front $\mu=\eta$.",
    "open_cylinder": r"The plane of $\eta$ and $r$ ($\phi=y=0$) of the open stiff universe, with the cylindrical front $r=\eta$.",
    "closed": r"The plane of $\eta$ and $\zeta$ ($\delta=\sigma=0$) of the closed stiff universe, from the axis $\zeta=0$ to the axis $\zeta=\pi/2$, with the toroidal front $\zeta=\eta$.",
}


def diagram_specs(Diagram, sideways):
    """One view of light rays and future cones per chart, with its wave front."""
    rows = []
    for system, vid, label, plane, box, fixed, where, front in PLANES:
        cylinder = vid.endswith("cylinder")
        params = {"gamma": 1, **({"A": 2} if system == "plane_wave" else {})}
        xname = plane[1]
        extra = {"where": where} if where else {}
        rows.append(Diagram(METRIC, system, vid, label, plane, box, f"${xname}$", "$\\eta$", params, fixed,
                            tau="eta", families=sideways if not cylinder else ("ingoing", "outgoing"),
                            curves=(front,), cones=(5, 5), **extra))
    return rows


def diagram_captions():
    return {(METRIC, system, vid): [CAPTION[vid]] for system, vid, *_ in PLANES}


def null_forms():
    """The two null families of every plane, eta + x and eta - x."""
    return {(METRIC, system, vid): ((lambda e, x: e + x), (lambda e, x: e - x), None)
            for system, vid, *_ in PLANES}


MOMENT = 0.7853981633974483   # eta = pi/4, the moment of greatest expansion of the closed universe


def slice_marks(slices, system, vid):
    """The embedding's moment eta = pi/4 on the closed universe's plane, from axis to axis."""
    def line(moment):
        lo, hi = moment.reach("closed_toroidal", "\\zeta")
        return slices.along(MOMENT, lo, hi)
    return slices.one(METRIC, line, label="$\\eta=\\pi/4$", view_id="hemisphere")


def conformal_views(ck, src, cf):
    """Each plane is a conformal factor times Minkowski's half plane eta > 0, by p, q = arctan(eta -+ x),
    the closed universe's square eta, zeta in (0, pi/2) under p, q = (eta -+ zeta)/2 as it stands."""
    import math
    import numpy as np

    pi = math.pi
    views = []
    for system, vid, label, plane_names, _, fixed, _, _ in PLANES:
        params = {"gamma": 1, **({"A": 2} if system == "plane_wave" else {})}
        plane = cf.Plane(src, METRIC, system, plane_names, fixed, params)
        closed = system == "closed_toroidal"
        cylinder = vid.endswith("cylinder")
        if closed:
            def pq(t, x):
                t, x = np.asarray(t, float), np.asarray(x, float)
                return (t - x) / 2, (t + x) / 2
            ck.chart(f"{METRIC}/{system}", plane, pq, ck.uniform(0.05, 1.5), ck.uniform(0.05, 1.5), lambda t, x: (1, 0))
            ck.diverges(f"{system}: K diverges toward the big bang", plane.kretschmann(0.01, 0.7), plane.kretschmann(0.001, 0.7))
            ck.diverges(f"{system}: K diverges toward the big crunch", plane.kretschmann(pi/2 - 0.01, 0.7),
                        plane.kretschmann(pi/2 - 0.001, 0.7))
            box = [-0.25, pi/2 + 0.25, -0.25, pi/2 + 0.25]
        else:
            pq = cf.mink_pq
            if system == "plane_wave":
                ck.chart(f"{METRIC}/{system}", plane, pq, ck.uniform(0.6, 5), ck.uniform(-0.5, 0.5), lambda t, x: (1, 0))
                ck.diverges(f"{system}: K diverges toward the big bang", plane.kretschmann(0.01, 0.0), plane.kretschmann(0.001, 0.0))
            else:
                lo = 0.01 if cylinder else -5
                ck.chart(f"{METRIC}/{system}", plane, pq, ck.uniform(0.03, 5), ck.uniform(lo, 5), lambda t, x: (1, 0))
                ck.diverges(f"{system}: K diverges toward the big bang", plane.kretschmann(0.01, 1.0), plane.kretschmann(0.001, 1.0))
                ck.limit(f"{system}: K tends to zero in the future", [plane.kretschmann(1e6, 1.0)], [0], 1e-12)
            box = [-0.25, pi + 0.25, -0.25, pi + 0.25] if cylinder else [-pi - 0.25, pi + 0.25, -0.25, pi + 0.25]
        v = cf.View(vid, label, box, system)
        if closed:
            square = [[0, 0], [pi/2, 0], [pi/2, pi/2], [0, pi/2]]
            v.fill("region", square)
            v.fill("cover", square)
            v.line("singular", [[[0, 0], [pi/2, 0]], [[0, pi/2], [pi/2, pi/2]]], zig=True)
            v.line("centre", [[[0, 0], [0, pi/2]], [[pi/2, 0], [pi/2, pi/2]]])
            v.line("null", [[[0, 0], [pi/2, pi/2]]])
            for t in (pi/8, pi/4, 3*pi/8):
                v.line("t", [[[0, t], [pi/2, t]]])
            v.legend("t", "$\\eta$ constant")
            v.legend("singular", "big bang $\\eta=0$ and big crunch $\\eta=\\pi/2$")
            v.legend("centre", "the axes $\\zeta=0$ and $\\zeta=\\pi/2$")
            v.legend("null", "toroidal front $\\zeta=\\eta$")
            moment = cf.slices.moments(METRIC, "hemisphere", label="$\\eta=\\pi/4$")[0]
            reach = np.linspace(*moment.reach("closed_toroidal", "\\zeta"), 81)
            v.slice(moment, [pq(np.full_like(reach, MOMENT), reach)])
            v.set(settings="$\\gamma=1$; drawn as it stands, $T=\\eta$, $X=\\zeta$.")
        else:
            whole = [[0, 0], [pi, 0], [0, pi]] if cylinder else [[-pi, 0], [pi, 0], [0, pi]]
            v.fill("region", whole)
            v.fill("cover", whole)
            v.line("singular", [[[0 if cylinder else -pi, 0], [pi, 0]]], zig=True)
            v.line("scri", [[[pi, 0], [0, pi]]] if cylinder else [[[-pi, 0], [0, pi]], [[pi, 0], [0, pi]]])
            if cylinder:
                v.line("centre", [[[0, 0], [0, pi]]])
                v.legend("centre", "the axis $r=0$")
            if system == "plane_wave":
                wave = [[0, 0], [pi/2, pi/2], [0, pi], [-pi/2, pi/2]]
                v.fill("cover2", wave)
                v.legend("cover2", "the chart, behind the front and ahead of $z=-\\eta$")
            v.line("null", [[[0, 0], [pi/2, pi/2]]])
            rr = cf.spread(0, np.inf, 500, 10)
            xs = rr if cylinder else np.concatenate([-rr[::-1], rr])
            for t in (0.25, 0.5, 1, 2, 4):
                v.curve("t", *pq(np.full_like(xs, t), xs))
            v.legend("t", "$\\eta$ constant")
            v.legend("singular", "big bang $\\eta=0$")
            v.legend("scri", "future null infinity")
            v.legend("null", "the wave front, $" + plane_names[1] + "=\\eta$")
            v.set(settings=("$\\gamma=1$, $A=2$; " if system == "plane_wave" else "$\\gamma=1$; ")
                  + "$p,q=\\arctan(\\eta\\mp " + plane_names[1] + ")$, $T=p+q$, $X=q-p$.")
        views.append(v)
    return views


CONFORMAL_CAPTION = {
    "wave": r"The flat stiff universe in $\eta$ and $z$, each point in the diagram a plane, the wave behind its front shaded.",
    "flat": r"The flat stiff universe in $\eta$ and $z$, each point in the diagram a plane.",
    "flat_cylinder": r"The flat stiff universe in $\eta$ and $r$, each point in the diagram a cylinder.",
    "open": r"The open stiff universe in $\eta$ and $\mu$, each point in the diagram a plane.",
    "open_cylinder": r"The open stiff universe in $\eta$ and $r$, each point in the diagram a cylinder.",
    "closed": r"The closed stiff universe in $\eta$ and $\zeta$, each point in the diagram a torus.",
}


def conformal_captions():
    return {(METRIC, vid): [CONFORMAL_CAPTION[vid]] for _, vid, *_ in PLANES}


def embedding_views(ck, src, em):
    """The closed universe at eta = pi/4, cut at sigma = 0: g_zetazeta = a^2 and g_deltadelta =
    a^2 sin^2 zeta with a^2 = gamma sin 2 eta, the hemisphere rho = a sin zeta, Z = a cos zeta."""
    import math
    import numpy as np

    sl = em.Slice(src, METRIC, "closed_toroidal", "\\zeta", "\\delta", {"eta": "pi/4", "sigma": "0"}, {"gamma": 1})
    size = 1.5
    cap = em.Piece("sphere", "sheet", sl, 0.0, math.pi / 2, -1.0, 1,
                   (("axis", "the axis $\\zeta=0$"), ("edge", "the axis $\\zeta=\\pi/2$")),
                   [(math.pi / 8, "r", None), (math.pi / 4, "r", None), (3 * math.pi / 8, "r", None)], size)
    ck.isometry("closed universe hemisphere", cap)
    ck.radius("closed universe: rho = a sin zeta", cap, np.sin, size)
    ck.form("closed universe: Z = -a cos zeta", cap, lambda z: -np.cos(z), size)
    surface = em.Surface([cap])
    fig = em.figure_of([surface], {"sheet": "cover"}, size)
    fig.legend("fill", "cover", "$\\eta=\\pi/4$, $\\sigma=0$")
    fig.legend("line", "r", "tori $\\zeta=\\pi/8$, $\\pi/4$, $3\\pi/8$")
    fig.legend("line", "meridian", "$\\delta$ constant")
    data = fig.done()
    return [em.view("hemisphere", "Closed universe", "$\\sqrt{\\gamma}$", [surface], data, settings="$\\gamma=1$.")]


def embedding_captions():
    return {(METRIC, "hemisphere"): [r"A hemisphere of the closed stiff universe at $\eta=\pi/4$ ($\sigma=0$), radius $a=\sqrt{\gamma}$, from the axis $\zeta=0$ at its pole to the axis $\zeta=\pi/2$ at its rim."]}


def curvature_views(ck, src, cv):
    """The Kretschmann and Ricci scalars across the wave region at eta = 1, against the flat
    universe's ahead of the front, R = -16 pi rho with rho = p = 3e^S/(64 pi gamma eta^3)."""
    xaxis = cv.axis("$z$", -1, 1)
    scales = {"K": {"kind": "linear", "from": -0.6, "to": 1.2}, "R": {"kind": "linear", "from": -0.8, "to": 0.0}}
    panels = [{"id": f, "name": "$\\gamma^2K$" if f == "K" else "$\\gamma R$", "aspect": 2, "scale": scales[f],
               "ticks": cv.value_ticks(scales[f])} for f in ("K", "R")]
    cv.shape_checks(ck, METRIC + "/wave", [p["aspect"] for p in panels])
    curves = []
    for field, scalar in (("K", "kretschmann"), ("R", "ricci_scalar")):
        for system, cls, params in (("plane_wave", "wave", {"gamma": 1, "A": 2}), ("flat_planar", "ahead", {"gamma": 1})):
            stored = cv.Held(cv.nr_published(src, METRIC, system, scalar, ["z", "eta"], list(params)), {"eta": 1})
            curves += cv.line_curves(ck, METRIC + "/wave", field, cls, stored, "z", xaxis, scales[field],
                                     cv.PLOT_W / 2, params)
    return [{"id": "wave", "label": "Across the wave", "kind": "line", "x": xaxis, "panels": panels,
             "curves": curves, "marks": [],
             "legend": [["line", "wave", "behind the front"], ["line", "ahead", "the flat universe ahead"]],
             "caption": [r"The Kretschmann scalar $K$ and the Ricci scalar $R$ of the wave region at $\eta=1$, from $z=-\eta$ to the front $z=\eta$. "
                         r"On the front $K$ meets the flat universe's $15/16$, and behind it $K$ falls below zero. "
                         r"Behind the front $R=-3e^{S}/4$ lies nearer zero, where the fluid is less dense."],
             "ricci_sentences": [[0, 2]], "settings": "$\\gamma=1$, $A=2$, $\\eta=1$.",
             "clip": "Linear vertical scales.", "per_scalar": cv.scalar_clips("Linear scale.", "Linear scale.")}]
