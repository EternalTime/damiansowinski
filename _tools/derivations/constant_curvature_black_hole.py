"""The charts of the black holes of constant curvature, anti-de Sitter space with its points identified
along a boost. constant_curvature_black_hole.md beside this file records each chart's source.

Four dimensions:

    kruskal       Banados, Gomberoff and Martinez, CQG 15, 3575 (1998), hep-th/9805087, eqs. (K), (r-y)
                  and (ds/krus): ds^2 = l^2(r + r_+)^2/r_+^2 (-dy_0^2 + dy_1^2 + dy_2^2) + r^2 dphi^2 with
                  r = r_+(1 + y^2)/(1 - y^2), y^2 = -y_0^2 + y_1^2 + y_2^2 between -1 and 1;
    static        their eqs. (schw1), (ds/sch) and (sph2), the static patch y_0 = f sin(theta) sinh(r_+ t/l),
                  y_1 = f sin(theta) cosh(r_+ t/l), y_2 = f cos(theta), f = sqrt((r - r_+)/(r + r_+));
    exterior      their eq. (ds2/sch), the whole exterior y_0 = f sinh(r_+ t/l), y_1 = f cos(theta)
                  cosh(r_+ t/l), y_2 = f sin(theta) cosh(r_+ t/l), with the square on the cosh that Cai's
                  eq. (2.15) and Anand, Campos Delgado and Samart's eq. (2.5) print;
    interior      Holst and Peldan, CQG 14, 3433 (1997), gr-qc/9705067, eqs. (SchwarzschildCoord) and
                  (SchwarzschildMetric), with their unit radius restored as l and their period 2 pi alpha
                  of psi written 2 pi r_+/l;
    null          their eqs. (uvCoordinates) and (uvMetric), tan u = e^-t sqrt((1 - r)/(1 + r)) and
                  tan v = e^t sqrt((1 - r)/(1 + r)) in the interior chart.

Five dimensions:

    kruskal_five  Banados, PRD 57, 1068 (1998), gr-qc/9703040, eqs. (y) and (ds/krus) at n = 5;
    static_five   his eqs. (cc), (ds/sch) and (sph), y_0 = f cos(theta) sinh(r_+ t/l), y_1 = f cos(theta)
                  cosh(r_+ t/l), y_2 = f sin(theta) sin(chi), y_3 = f sin(theta) cos(chi);
    exterior_five Cai, PLB 544, 176 (2002), hep-th/0206223, eqs. (2.13) to (2.15).

Every chart is checked at points: its Riemann tensor is -(g_ac g_bd - g_ad g_bc)/l^2, and it is the
pullback of its source's map, Banados's embedding in flat space of two times for the Kruskal charts,
the Kruskal chart for the static and exterior ones, Holst and Peldan's embedding for the interior, and
the interior chart for the null one.
"""
import sympy as sp

METRIC = "constant_curvature_black_hole"
CHARTS = ("kruskal", "static", "exterior", "interior", "null", "kruskal_five", "static_five", "exterior_five")
NAMES = {"kruskal": "Kruskal", "static": "Static", "exterior": "Exterior",
         "interior": "Holst-Peld\u00e1n Interior", "null": "Holst-Peld\u00e1n Null",
         "kruskal_five": "Kruskal, Five Dimensions", "static_five": "Static, Five Dimensions",
         "exterior_five": "Exterior, Five Dimensions"}

REALS = " \\in (-\\infty, \\infty)"
TURN = " \\in [0, 2\\pi)"
PHI = "\\varphi"


def _y2(n):
    """-y_0^2 + y_1^2 + ... + y_{n-2}^2, written out."""
    return "-y_0^2 + " + " + ".join(f"y_{i}^2" for i in range(1, n - 1))


def _radius(n):
    """Banados's r = r_+(1 + y^2)/(1 - y^2) in the Kruskal coordinates of n dimensions."""
    plus = " + ".join(f"y_{i}^2" for i in range(1, n - 1))
    minus = " - ".join(f"y_{i}^2" for i in range(1, n - 1))
    return f"r = r_+\\dfrac{{1 - y_0^2 + {plus}}}{{1 + y_0^2 - {minus}}}"


def _kruskal_line(n):
    flat = " + ".join(f"dy_{i}^2" for i in range(1, n - 1))
    return ("ds^2 = \\dfrac{l^2\\left(r + r_+\\right)^2}{r_+^2}\\left(-dy_0^2 + " + flat + "\\right) + r^2d\\varphi^2")


STATIC_PART = "\\left(r^2 - r_+^2\\right)"
RADIAL = " + \\dfrac{l^2\\,dr^2}{r^2 - r_+^2} + r^2d\\varphi^2"
COSH = "\\cosh^2\\left(\\dfrac{r_+t}{l}\\right)"

LINES = {
    "kruskal": _kruskal_line(4),
    "static": "ds^2 = " + STATIC_PART + "\\left(-\\sin^2\\theta\\,dt^2 + \\dfrac{l^2}{r_+^2}d\\theta^2\\right)" + RADIAL,
    "exterior": "ds^2 = " + STATIC_PART + "\\left(-dt^2 + \\dfrac{l^2}{r_+^2}" + COSH + "\\,d\\theta^2\\right)" + RADIAL,
    "interior": ("ds^2 = l^2\\left(-\\dfrac{dr^2}{1 - r^2} + \\left(1 - r^2\\right)dt^2 + r^2d\\psi^2"
                 " + \\left(1 - r^2\\right)\\sinh^2t\\,d\\varphi^2\\right)"),
    # Holst and Peldan's eq. (uvMetric) with alpha = u + v and beta = v - u.
    "null": ("ds^2 = \\dfrac{l^2}{\\cos^2\\beta}\\left(-4\\,du\\,dv + \\cos^2\\alpha\\,d\\psi^2"
             " + \\sin^2\\beta\\,d\\varphi^2\\right)"),
    "kruskal_five": _kruskal_line(5),
    "static_five": ("ds^2 = " + STATIC_PART + "\\left(-\\cos^2\\theta\\,dt^2 + \\dfrac{l^2}{r_+^2}\\left(d\\theta^2"
                    " + \\sin^2\\theta\\,d\\chi^2\\right)\\right)" + RADIAL),
    "exterior_five": ("ds^2 = " + STATIC_PART + "\\left(-dt^2 + \\dfrac{l^2}{r_+^2}" + COSH + "\\left(d\\theta^2"
                      " + \\sin^2\\theta\\,d\\chi^2\\right)\\right)" + RADIAL),
}

COORDS = {
    "kruskal": ["y_0", "y_1", "y_2", PHI],
    "static": ["t", "r", "\\theta", PHI],
    "exterior": ["t", "r", "\\theta", PHI],
    "interior": ["r", "t", "\\psi", PHI],
    "null": ["u", "v", "\\psi", PHI],
    "kruskal_five": ["y_0", "y_1", "y_2", "y_3", PHI],
    "static_five": ["t", "r", "\\theta", "\\chi", PHI],
    "exterior_five": ["t", "r", "\\theta", "\\chi", PHI],
}

PARAMETERS = {system: ["l", "r_+"] for system in CHARTS}
PARAMETERS["kruskal"] = ["l", "r_+", _radius(4)]
PARAMETERS["kruskal_five"] = ["l", "r_+", _radius(5)]
PARAMETERS["null"] = ["l", "r_+", "\\alpha = u + v", "\\beta = v - u"]


def _kruskal_domains(n):
    ys = [f"y_{i}" + REALS for i in range(n - 1)]
    rest = " + ".join(f"y_{i}^2" for i in range(1, n - 1))
    return ys + [PHI + TURN, f"-1 < {_y2(n)} < 1",
                 f"y_0^2 = {rest} \\;\\text{{(the event horizon)}}",
                 f"y_0^2 = 1 + {rest} \\;\\text{{(the singularity)}}"]


PSI = "\\psi \\in [0, 2\\pi r_+/l)"
DOMAINS = {
    "kruskal": _kruskal_domains(4),
    "static": ["t" + REALS, "r \\in (r_+, \\infty)", "\\theta \\in (0, \\pi)", PHI + TURN],
    "exterior": ["t" + REALS, "r \\in (r_+, \\infty)", "\\theta" + TURN, PHI + TURN],
    "interior": ["r \\in (0, 1)", "t \\in [0, \\infty)", PSI, PHI + TURN,
                 "r = 0 \\;\\text{(the singularity)}", "r = 1 \\;\\text{(the event horizon)}",
                 "r = 1/2 \\;\\text{(the apparent horizon)}"],
    "null": ["u \\in (0, \\pi/4)", "v \\in [u, \\pi/2 - u)", PSI, PHI + TURN,
             "u = 0 \\;\\text{(the event horizon)}", "u + v = \\pi/2 \\;\\text{(the singularity)}",
             "\\tan u\\tan v = 1/3 \\;\\text{(the apparent horizon)}"],
    "kruskal_five": _kruskal_domains(5),
    "static_five": ["t" + REALS, "r \\in (r_+, \\infty)", "\\theta \\in [0, \\pi/2)", "\\chi" + TURN, PHI + TURN],
    "exterior_five": ["t" + REALS, "r \\in (r_+, \\infty)", "\\theta \\in [0, \\pi]", "\\chi" + TURN, PHI + TURN],
}


def _dimensions(system):
    dims = {"l": "L", "r_+": "L"}
    for name in COORDS[system]:
        dims[name] = "L" if (name == "r" and system not in ("interior",)) else "1"
    if system.startswith("kruskal"):
        dims["r"] = "L"
    if system == "null":
        dims.update({"\\alpha": "1", "\\beta": "1"})
    return dims


DIMENSIONS = {(METRIC, system): _dimensions(system) for system in CHARTS}
# The radius r of the Kruskal charts is written out wherever it stands: held with its rates, the
# curvature is constant only through r's relation with y^2, which the rates leave out.
# The null chart holds alpha = u + v and beta = v - u, a pair of coordinates with no relation between
# them, so that its values are printed in the cosines and sines of the two.
HELD = {(METRIC, "null"): ("alpha", "beta")}
RATES = {(METRIC, "null"): {"alpha": {"u": "1", "v": "1"}, "beta": {"u": "-1", "v": "1"}}}

# The Kretschmann scalar of anti-de Sitter space of n dimensions, 2n(n - 1)/l^4.
KRETSCHMANN = {system: ("\\dfrac{40}{l^4}" if system.endswith("_five") else "\\dfrac{24}{l^4}") for system in CHARTS}


def register(module):
    module.CHARTS[METRIC] = [lambda system=system: chart_spec(system) for system in CHARTS]


def trigonometric(value):
    """A value of the static chart of five dimensions in secants and tangents of theta, where the
    factored form writes 1 - sin^2 theta as a product of two brackets."""
    return sp.factor(sp.trigsimp(sp.factor(sp.sympify(value))))


def folded(x):
    """A factored value with each pair (sin a + 1)^k (sin a - 1)^k written (-1)^k cos^2k a, and each
    pair (cos a + 1)^k (cos a - 1)^k written (-1)^k sin^2k a."""
    factors = {}
    rest = sp.Integer(1)
    for f in sp.Mul.make_args(x):
        base, k = (f.base, f.exp) if f.is_Pow else (f, sp.Integer(1))
        factors[base] = factors.get(base, 0) + k
    out = sp.Integer(1)
    done = set()
    for base, k in factors.items():
        if base in done:
            continue
        for trig, other in ((sp.sin, sp.cos), (sp.cos, sp.sin)):
            if isinstance(base, sp.Add) and len(base.args) == 2:
                terms = dict(sp.Add.make_args(base).__iter__() and [(t.is_Number, t) for t in base.args])
                number, call = terms.get(True), terms.get(False)
                if number in (1, -1) and isinstance(call, trig):
                    partner = call - number
                    if factors.get(partner) == k:
                        out *= (-1) ** k * other(call.args[0]) ** (2 * k)
                        done |= {base, partner}
                        break
        else:
            out *= base ** k
            done.add(base)
            continue
    return out * rest


def held_angles(value):
    """A value of the null chart in the cosines and sines of alpha and beta, each printed as a name
    the reader reads back as the held function it stands for."""
    from sympy.core.function import AppliedUndef
    x = sp.sympify(value)
    x = x.xreplace({f: sp.Symbol(f.func.__name__, real=True) for f in x.atoms(AppliedUndef)})
    for _ in range(2):
        x = sp.expand_trig(x)
    # Every even power of a sine written in cosines, sin^2 = 1 - cos^2, before the value is factored.
    alpha, beta = sp.Symbol("alpha", real=True), sp.Symbol("beta", real=True)
    sa, ca, sb, cb = sp.symbols("_sa _ca _sb _cb")
    names = {sp.sin(alpha): sa, sp.cos(alpha): ca, sp.sin(beta): sb, sp.cos(beta): cb}
    x = sp.together(x.xreplace(names))
    relation = [sa ** 2 + ca ** 2 - 1, sb ** 2 + cb ** 2 - 1]
    num, den = (sp.expand(sp.reduced(sp.expand(part), relation, sa, sb, ca, cb)[1]) for part in sp.fraction(x))
    x = sp.factor(sp.cancel(num / den)).xreplace({v: k for k, v in names.items()})
    return folded(x)


def chart_spec(system):
    """One chart: its line element, which holds no c since no coordinate has dimensions of time,
    the Kretschmann scalar of anti-de Sitter space, and the check of the chart at points."""
    import chart_printer as cp
    import verify_metrics as vm
    held = HELD.get((METRIC, system), ())
    rates = RATES.get((METRIC, system), {})
    probe = vm.Reader(COORDS[system], PARAMETERS[system], (), held=held, rates=rates)
    pretty = {}
    if system in ("exterior", "exterior_five"):
        # cosh(r_+ t/l), which the checker's Geometry hands back as exponentials.
        P = probe.parameters
        pretty = {"pretty": cp.hyperbolic(P["r_+"] * probe.symbol["t"] / P["l"])}
    elif system == "interior":
        pretty = {"pretty": cp.hyperbolic(probe.symbol["t"])}
    elif system == "static_five":
        pretty = {"pretty": trigonometric}
    elif system == "null":
        pretty = {"pretty": held_angles}
    return {
        **pretty,
        "metric_id": METRIC,
        "system": {"id": system, "name": NAMES[system], "coords": COORDS[system], "domains": DOMAINS[system],
                   "parameters": PARAMETERS[system], "line_element": LINES[system]},
        "chart_line_element": LINES[system],
        "printer": {},
        "reduce": probe.surface if rates else None,
        "kretschmann": KRETSCHMANN[system],
        "check": lambda chart: check_chart(chart, system),
    }


# Parameters (l, r_+) and events of every check, inside the domain of each chart.
L_VALUE, RP_VALUE = sp.Rational(7, 5), sp.Rational(3, 5)
EVENTS = {
    "kruskal": [(sp.Rational(1, 5), sp.Rational(2, 5), sp.Rational(1, 3), sp.Rational(1, 2)),
                (sp.Rational(3, 5), sp.Rational(1, 5), sp.Rational(1, 4), sp.Integer(2)),
                (sp.Rational(1, 10), sp.Rational(7, 10), sp.Rational(1, 2), sp.Integer(3))],
    "static": [(sp.Rational(2, 5), sp.Integer(1), sp.Rational(1, 3), sp.Rational(1, 2)),
               (sp.Integer(-1), sp.Rational(9, 4), sp.Integer(2), sp.Integer(1)),
               (sp.Rational(3, 2), sp.Rational(4, 5), sp.Rational(5, 4), sp.Integer(2))],
    "exterior": [(sp.Rational(2, 5), sp.Integer(1), sp.Rational(1, 3), sp.Rational(1, 2)),
                 (sp.Integer(-1), sp.Rational(9, 4), sp.Integer(4), sp.Integer(1)),
                 (sp.Rational(3, 2), sp.Rational(4, 5), sp.Rational(5, 4), sp.Integer(2))],
    "interior": [(sp.Rational(1, 3), sp.Rational(1, 2), sp.Rational(1, 5), sp.Integer(1)),
                 (sp.Rational(7, 10), sp.Integer(2), sp.Rational(1, 2), sp.Integer(3)),
                 (sp.Rational(1, 10), sp.Rational(1, 5), sp.Integer(1), sp.Integer(2))],
    "null": [(sp.Rational(1, 5), sp.Rational(3, 5), sp.Rational(1, 5), sp.Integer(1)),
             (sp.Rational(1, 10), sp.Integer(1), sp.Rational(1, 2), sp.Integer(3)),
             (sp.Rational(1, 2), sp.Rational(7, 10), sp.Integer(1), sp.Integer(2))],
    "kruskal_five": [(sp.Rational(1, 5), sp.Rational(2, 5), sp.Rational(1, 3), sp.Rational(1, 4), sp.Rational(1, 2)),
                     (sp.Rational(3, 5), sp.Rational(1, 5), sp.Rational(1, 4), sp.Rational(1, 6), sp.Integer(2)),
                     (sp.Rational(1, 10), sp.Rational(1, 2), sp.Rational(1, 2), sp.Rational(1, 3), sp.Integer(3))],
    "static_five": [(sp.Rational(2, 5), sp.Integer(1), sp.Rational(1, 3), sp.Rational(1, 2), sp.Integer(1)),
                    (sp.Integer(-1), sp.Rational(9, 4), sp.Rational(6, 5), sp.Integer(4), sp.Integer(2)),
                    (sp.Rational(3, 2), sp.Rational(4, 5), sp.Rational(1, 5), sp.Integer(2), sp.Integer(3))],
    "exterior_five": [(sp.Rational(2, 5), sp.Integer(1), sp.Rational(1, 3), sp.Rational(1, 2), sp.Integer(1)),
                      (sp.Integer(-1), sp.Rational(9, 4), sp.Rational(5, 2), sp.Integer(4), sp.Integer(2)),
                      (sp.Rational(3, 2), sp.Rational(4, 5), sp.Rational(1, 5), sp.Integer(2), sp.Integer(3))],
}


def kruskal_metric(n, ys, phi, l, rp):
    """Banados's eq. (ds/krus) of n dimensions as a matrix in (y_0, ..., y_{n-2}, phi), written out."""
    y2 = -ys[0] ** 2 + sum(y ** 2 for y in ys[1:])
    r = rp * (1 + y2) / (1 - y2)
    eta = [-1] + [1] * (n - 2)
    factor = l ** 2 * (r + rp) ** 2 / rp ** 2
    return sp.diag(*[factor * e for e in eta], r ** 2)


def pulled(matrix, functions, coords):
    """The metric `matrix`, in the variables `functions` are values of, pulled back to `coords`."""
    jacobian = sp.Matrix(len(functions), len(coords), lambda a, b: sp.diff(functions[a], coords[b]))
    return jacobian.T * matrix * jacobian


def check_chart(chart, system):
    """The chart at points: Riemann of constant curvature -1/l^2, and the pullback of its source's map."""
    import verify_metrics as vm
    geo, g = chart.geo, chart.geo.g
    X = list(chart.symbols)
    n = len(X)
    P = chart.reader.parameters
    l, rp = P["l"], P["r_+"]
    held = chart.reader.held or {}
    name = f"{METRIC}/{system}"

    def written(e):
        e = sp.sympify(e)
        return e.subs(held).doit() if held else e

    def at(e, event):
        values = {l: L_VALUE, rp: RP_VALUE, **dict(zip(X, event))}
        return complex(sp.N(written(e).subs(values), 30))

    def agree(a, b, what):
        for event in EVENTS[system]:
            if abs(at(a, event) - at(b, event)) > 1e-20:
                raise AssertionError(f"{name}: {what} fails at {event}")

    # Constant curvature: R_abcd = -(g_ac g_bd - g_ad g_bc)/l^2.
    riemann = geo.riemann_llll()
    for a in range(n):
        for b in range(a + 1, n):
            for c in range(n):
                for d in range(c + 1, n):
                    expected = -(g[a, c] * g[b, d] - g[a, d] * g[b, c]) / l ** 2
                    agree(vm._at(riemann, (a, b, c, d)), expected,
                          f"constant curvature in R_{{{chart.coords_tex[a]}{chart.coords_tex[b]}"
                          f"{chart.coords_tex[c]}{chart.coords_tex[d]}}}")

    target = None
    if system.startswith("kruskal"):
        # Banados's eq. (y): x_alpha = 2l y_alpha/(1 - y^2), x_{n-1} = (l r/r_+) sinh(r_+ phi/l),
        # x_n = (l r/r_+) cosh(r_+ phi/l), on -x_0^2 + x_1^2 + ... + x_{n-1}^2 - x_n^2 = -l^2.
        ys, phi = X[:-1], X[-1]
        y2 = -ys[0] ** 2 + sum(y ** 2 for y in ys[1:])
        r = rp * (1 + y2) / (1 - y2)
        x = [2 * l * y / (1 - y2) for y in ys] + [l * r / rp * sp.sinh(rp * phi / l), l * r / rp * sp.cosh(rp * phi / l)]
        ambient = sp.diag(-1, *([1] * (n - 1)), -1)
        agree(sum(ambient[i, i] * x[i] ** 2 for i in range(n + 1)), -l ** 2, "Banados's hyperboloid (ads/n)")
        target = pulled(ambient, x, X)
    elif system in ("static", "exterior", "static_five", "exterior_five"):
        t, r = X[0], X[1]
        f = sp.sqrt((r - rp) / (r + rp))
        s, ch = sp.sinh(rp * t / l), sp.cosh(rp * t / l)
        theta = X[2]
        if system == "static":
            ys = [f * sp.sin(theta) * s, f * sp.sin(theta) * ch, f * sp.cos(theta)]
        elif system == "exterior":
            ys = [f * s, f * sp.cos(theta) * ch, f * sp.sin(theta) * ch]
        elif system == "static_five":
            chi = X[3]
            ys = [f * sp.cos(theta) * s, f * sp.cos(theta) * ch, f * sp.sin(theta) * sp.sin(chi),
                  f * sp.sin(theta) * sp.cos(chi)]
        else:
            chi = X[3]
            ys = [f * s, f * sp.cos(theta) * ch, f * sp.sin(theta) * sp.cos(chi) * ch,
                  f * sp.sin(theta) * sp.sin(chi) * ch]
        Y = sp.symbols(f"Y0:{n - 1}", real=True)
        Phi = sp.Symbol("Phi", real=True)
        kruskal = kruskal_metric(n, list(Y), Phi, l, rp)
        jac = sp.Matrix(n, n, lambda a, b: sp.diff((ys + [X[-1]])[a], X[b]))
        target = jac.T * kruskal.subs(dict(zip(Y, ys))) * jac
    elif system == "interior":
        # Holst and Peldan's eq. (SchwarzschildCoord) at radius l, on X^2 + Y^2 + Z^2 - U^2 - V^2 = -l^2.
        r, t, psi, phi = X
        w = sp.sqrt(1 - r ** 2)
        x = [-l * w * sp.sinh(t) * sp.cos(phi), l * r * sp.sinh(psi), -l * w * sp.sinh(t) * sp.sin(phi),
             l * r * sp.cosh(psi), -l * w * sp.cosh(t)]
        ambient = sp.diag(1, 1, 1, -1, -1)
        agree(sum(ambient[i, i] * x[i] ** 2 for i in range(5)), -l ** 2, "Holst and Peldan's hyperboloid")
        target = pulled(ambient, x, X)
    elif system == "null":
        # Their eq. (uvCoordinates) solved for r and t: tan u tan v = (1 - r)/(1 + r), e^2t = tan v/tan u.
        u, v, psi, phi = X
        q = sp.tan(u) * sp.tan(v)
        r = (1 - q) / (1 + q)
        t = sp.log(sp.tan(v) / sp.tan(u)) / 2
        R, T, S, F = sp.symbols("R T S F", real=True)
        interior = sp.diag(-l ** 2 / (1 - R ** 2), l ** 2 * (1 - R ** 2), l ** 2 * R ** 2,
                           l ** 2 * (1 - R ** 2) * sp.sinh(T) ** 2)
        jac = sp.Matrix(4, 4, lambda a, b: sp.diff([r, t, psi, phi][a], X[b]))
        target = jac.T * interior.subs({R: r, T: t}) * jac
    for a in range(n):
        for b in range(a, n):
            agree(g[a, b], target[a, b], f"the pullback in slot {chart.coords_tex[a]}{chart.coords_tex[b]}")
