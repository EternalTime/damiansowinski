"""The domain wall of heterotic M-theory: the two flat walls at the ends of the interval of
Horava and Witten's theory reduced on a Calabi-Yau space to gauged N = 1 supergravity in five
dimensions, with the volume V of the Calabi-Yau space as its scalar. Lukas, Ovrut, Stelle and
Waldram, Phys. Rev. D 59, 086001 (1999), hep-th/9803235.

Every chart is in the signature (-,+,+,+,+), with alpha the inverse length of their constant
and phi = ln V:

    static       Lukas, Ovrut, Stelle and Waldram's solution, their eq. (sol),
                 ds^2 = a_0^2 H eta + b_0^2 H^4 dy^2, H = (sqrt 2/3) alpha |y| + c_0, V = b_0 H^3,
                 on the circle -pi rho <= y <= pi rho with y and -y one point;
    conformal    Lehners, McFadden and Turok's conformally flat form, their eq. (domainwall),
                 hep-th/0611259, ds^2 = h^{2/5}(eta + dw^2), V = h^{6/5}, h linear on -w_0 <= w <= w_0,
                 with their alpha = alpha/(3 sqrt 2) written in this alpha and their c written h_0;
    static_bulk  Chamblin and Reall's static bulk, hep-th/9903225, in Lehners, McFadden and Turok's
                 form, their eq. (staticbulk), -f dt^2 + r^12 dr^2/f + r^2 dx^2 with f = alpha^2 r^2 - mu r^4
                 in their alpha and e^phi = r^6; their time is written (c/alpha)t in their alpha, so
                 that t is a time, and mu is written lambda^2 times their alpha^2, their own
                 lambda = sqrt(mu)/alpha;
    rolling      Lukas, Ovrut and Waldram's moving walls, Phys. Rev. D 60, 086001 (1999),
                 hep-th/9806022, their eqs. (burt6), (burt15) and (burt16) in the gauge n = 1:
                 ds^2 = a_0^2 H(-c^2dt^2 + A^2(t - t_0)^{2p} dx^2) + B^2(t - t_0)^{2q} b_0^2 H^4 dy^2,
                 V = B(t - t_0)^q b_0 H^3, on the branch t > t_0.

The field equations of the bulk are those of the action of Lehners, McFadden and Turok, their
(Action5d), R - (1/2)(d phi)^2 - 6 alpha_LMT^2 e^{-2 phi}, which is Lukas, Ovrut, Stelle and
Waldram's (actparts) with the other fields zero and their curvature sign turned to this one:
    R_ab = (1/2) d_a phi d_b phi + (alpha^2/9) e^{-2 phi} g_ab,   box phi = -(2 alpha^2/3) e^{-2 phi},
and the walls are the boundary terms of both, Lehners, McFadden and Turok's (eom1) and (eom2):
the delta of y in G_mu nu is sqrt 2 alpha V^{-1} g_mu nu delta(y)/sqrt(g_yy), and in box phi it is
2 sqrt 2 alpha V^{-1} delta(y)/sqrt(g_yy). heterotic_m_theory_wall.md beside this file derives
every chart and records every drawing decision.
"""
import sympy as sp

METRIC = "heterotic_m_theory_wall"
CHARTS = ("static", "conformal", "static_bulk", "rolling")
NAMES = {"static": "Lukas, Ovrut, Stelle, and Waldram, Static Walls",
         "conformal": "Lehners, McFadden, and Turok, Conformally Flat",
         "static_bulk": "Chamblin and Reall, Static Bulk",
         "rolling": "Lukas, Ovrut, and Waldram, Rolling Walls"}

ALONG = ["x_1", "x_2", "x_3"]
FLAT = "dx_1^2 + dx_2^2 + dx_3^2"
H_DEF = "H = \\frac{\\sqrt{2}}{3}\\alpha|y| + c_0"
H_DEF_W = "h = \\frac{5\\sqrt{2}}{6}\\alpha\\left(w + w_0\\right) + h_0"
COORDS = {"static": ["t"] + ALONG + ["y"], "conformal": ["t"] + ALONG + ["w"],
          "static_bulk": ["t"] + ALONG + ["r"], "rolling": ["t"] + ALONG + ["y"]}
PARAMETERS = {
    "static": ["\\alpha", "a_0", "b_0", "c_0", "\\rho", H_DEF, "V = b_0 H^3"],
    "conformal": ["\\alpha", "w_0", "h_0", H_DEF_W, "V = h^{6/5}"],
    "static_bulk": ["\\alpha", "\\lambda", "V = r^6"],
    "rolling": ["\\alpha", "a_0", "b_0", "c_0", "\\rho", "A", "B", "t_0", "p", "q", H_DEF,
                "a = A\\left(t - t_0\\right)^p", "b = B\\left(t - t_0\\right)^q", "V = b_0\\,b\\,H^3"],
}
# The rolling chart holds its two scale factors as names, so every value is rational in a, b, the
# time since t_0, H and the parameters. Its tensors are built in the chart x^0 = ct, where the
# names are written in x^0/c, beside the file's own definitions in t, as Newman and Tamburino's are.
HELD = {(METRIC, "rolling"): ("a", "b")}
RATES = {(METRIC, "rolling"): {"a": {"t": "\\dfrac{p\\,a}{t - t_0}"}, "b": {"t": "\\dfrac{q\\,b}{t - t_0}"}}}
CHART_PARAMETERS = [p.replace("\\left(t - t_0\\right)", "\\left(\\dfrac{t}{c} - t_0\\right)") for p in PARAMETERS["rolling"]]
CHART_RATES = {"a": {"t": "\\dfrac{p\\,a}{t - c\\,t_0}"}, "b": {"t": "\\dfrac{q\\,b}{t - c\\,t_0}"}}


def line(system, c2):
    """The line element of a chart, with c2 standing before dt^2."""
    if system == "static":
        return "ds^2 = a_0^2H\\left(-" + c2 + "dt^2 + " + FLAT + "\\right) + b_0^2H^4dy^2"
    if system == "conformal":
        return "ds^2 = h^{2/5}\\left(-" + c2 + "dt^2 + " + FLAT + " + dw^2\\right)"
    if system == "static_bulk":
        return ("ds^2 = -r^2\\left(1 - \\lambda^2r^2\\right)" + c2 + "dt^2 + \\dfrac{18r^{10}}{\\alpha^2\\left(1 - "
                "\\lambda^2r^2\\right)}dr^2 + r^2\\left(" + FLAT + "\\right)")
    return "ds^2 = a_0^2H\\left(-" + c2 + "dt^2 + a^2\\left(" + FLAT + "\\right)\\right) + b^2b_0^2H^4dy^2"


REALS = " \\in (-\\infty, \\infty)"
INTERVAL = ["y \\in [-\\pi\\rho, \\pi\\rho]",
            "y = 0 \\;\\text{(a wall)}",
            "y = \\pm\\pi\\rho \\;\\text{(the other wall)}",
            "H > 0"]
DOMAINS = {
    "static": ["t" + REALS] + [x + REALS for x in ALONG] + INTERVAL,
    "conformal": ["t" + REALS] + [x + REALS for x in ALONG]
    + ["w \\in [-w_0, w_0]", "w = -w_0 \\;\\text{(the wall of negative tension)}",
       "w = w_0 \\;\\text{(the wall of positive tension)}", "h > 0"],
    "static_bulk": ["t" + REALS] + [x + REALS for x in ALONG]
    + ["r \\in (0, 1/\\lambda)", "\\lambda \\ge 0", "r = 0 \\;\\text{(the timelike naked singularity)}",
       "r = 1/\\lambda \\;\\text{(the cosmological horizon)}"],
    "rolling": ["t \\in (t_0, \\infty)"] + [x + REALS for x in ALONG] + INTERVAL
    + ["t = t_0 \\;\\text{(the curvature singularity)}"],
}
_DIMS = {"t": "T", "x_1": "L", "x_2": "L", "x_3": "L", "y": "L", "w": "L", "r": "1", "\\alpha": "1/L",
         "a_0": "1", "b_0": "1", "c_0": "1", "\\rho": "L", "H": "1", "V": "1", "w_0": "L", "h_0": "1", "h": "1",
         "\\lambda": "1", "A": "1", "B": "1", "t_0": "T", "p": "1", "q": "1", "a": "1", "b": "1"}
DIMENSIONS = {(METRIC, s): {name: _DIMS[name] for name in COORDS[s] + [p.split(" = ")[0] for p in PARAMETERS[s]]}
              for s in CHARTS}
DIGITS = 30
TOLERANCE = sp.Float("1e-22")


def register(module):
    module.CHARTS[METRIC] = [lambda system=system: chart_spec(system) for system in CHARTS]


def _in_h(probe, across, name, slope, at_zero):
    """A `then` for a value whose coordinate across the interval stands only inside a function
    linear in it, which is printed as that function's name: the coordinate is written through the
    name, and the value factored. `across` is the coordinate or the placeholder of its absolute
    value, `slope` and `at_zero` the function's slope and its value where `across` vanishes."""
    named = sp.Symbol("_" + name, positive=True)

    def then(value):
        return sp.factor(sp.sympify(value).subs(across, (named - at_zero) / slope))
    return then, named


def _held_apart(value, change):
    """`change` made to a value with the held names and their derivatives set aside, so that they
    keep their arguments, as chart_printer.Chart.named_time keeps them."""
    from sympy.core.function import AppliedUndef
    value = sp.sympify(value)
    held = {}
    for kind in (sp.Derivative, AppliedUndef):
        for atom in sorted(value.atoms(kind), key=sp.default_sort_key):
            held[atom] = sp.Dummy()
        value = value.xreplace(held)
    value = change(value)
    back = {dummy: atom for atom, dummy in held.items()}
    for _ in range(2):
        value = value.xreplace(back)
    return value


def _named_time(value, time, c):
    """A value in the chart x^0 = ct with its bare time written as c times the time."""
    return _held_apart(value, lambda v: v.subs(time, c * time))


def chart_spec(system):
    """One of the four charts, its values printed in H or h where the chart has one."""
    import chart_printer as cp
    import verify_metrics as vm

    coords = COORDS[system]
    probe = vm.Reader(coords, PARAMETERS[system], ())
    P = probe.parameters
    alpha = P["alpha"]
    spec = {"metric_id": METRIC,
            "system": {"id": system, "name": NAMES[system], "coords": coords, "domains": DOMAINS[system],
                       "parameters": PARAMETERS[system], "line_element": line(system, "c^2")},
            "chart_line_element": line(system, ""),
            "check": lambda chart: chart_check(chart, system)}
    if system in ("static", "rolling"):
        y = probe.symbol["y"]
        absolute = sp.Symbol("_abs" + y.name, positive=True)
        then, named = _in_h(probe, absolute, "H", sp.sqrt(2) * alpha / 3, P["c_0"])
        kinked, overrides = cp.kink(y, then)
        delta = sp.DiracDelta(y)

        def pretty(value):
            # The bulk and the wall are printed apart: the part without the delta as the kink
            # prints it, in H, and the delta with its coefficient taken on the wall, where H = c_0.
            value = sp.sympify(value)
            if not value.has(delta):
                return kinked(value)
            regular = value.subs(delta, 0)
            on_wall = sp.simplify(sp.simplify((value - regular) / delta).subs(y, 0))
            if on_wall.has(delta) or on_wall.has(sp.sign) or on_wall.has(y):
                raise ValueError(f"{value} is not a function plus a delta on the wall")
            return kinked(regular) + kinked(on_wall * delta)
        placeholder = next(p for p in overrides if p.name.startswith("_delta"))
        overrides[named] = "H"
        lead = [alpha, P["a_0"], P["b_0"], named]
        if system == "rolling":
            held = vm.Reader(coords, CHART_PARAMETERS, (), held=HELD[(METRIC, system)], rates=CHART_RATES)
            lead = [alpha, P["a_0"], P["b_0"], named]
            spec["time"] = "t"
            spec["chart_parameters"] = CHART_PARAMETERS
            spec["reduce"] = lambda value: vm.norm(held.by_rates(value))
        spec["pretty"] = spec["bracketed"] = pretty
        if system == "rolling":
            # The time stands only as the time since t_0, printed as one bracket. The printer hands
            # a value to `bracketed` in the chart x^0 = ct, where `pretty` has it with x^0 written
            # as c times the time, so the bracket writes the time first.
            time, since = probe.symbol["t"], sp.Symbol("_since", positive=True)
            overrides[since] = "\\left(t - t_0\\right)"
            kinked_in_h = pretty

            def pretty(value):
                return kinked_in_h(_held_apart(value, lambda v: v.subs(time, since + P["t_0"])))
            spec["pretty"] = pretty
            spec["bracketed"] = lambda value: pretty(_named_time(value, time, probe.c))
        spec["printer"] = {"lead": lead, "last": [placeholder], "overrides": overrides}
        spec["kretschmann_where"] = "y"
    elif system == "conformal":
        w = probe.symbol["w"]
        then, named = _in_h(probe, w, "h", 5 * sp.sqrt(2) * alpha / 6, 5 * sp.sqrt(2) * alpha * P["w_0"] / 6 + P["h_0"])
        spec["pretty"] = then
        spec["printer"] = {"lead": [alpha, named], "overrides": {named: "h"}}
    else:
        # Every value holds 1 - lambda^2 r^2 as the line element writes it, a factor bracketed whole.
        r, lam = probe.symbol["r"], P["lambda"]
        F = sp.Symbol("_F", positive=True)

        def fold(side):
            # Each pair (lambda r - 1)(lambda r + 1) of a factored side is -(1 - lambda^2 r^2).
            powers = {}
            for f in sp.Mul.make_args(side):
                base, k = f.as_base_exp()
                powers[base] = powers.get(base, 0) + k
            minus, plus = lam * r - 1, lam * r + 1
            pairs = min(powers.get(minus, 0), powers.get(plus, 0))
            out = (-F) ** pairs
            for base, k in powers.items():
                out *= base ** (k - pairs if base in (minus, plus) else k)
            return out

        def pretty(value):
            numerator, denominator = sp.fraction(sp.factor(value))
            return fold(numerator) / fold(denominator)
        spec["pretty"] = pretty
        spec["printer"] = {"lead": [alpha, lam, r], "named": {F: "1 - \\lambda^2r^2"}}
    return spec


# -- the checks ------------------------------------------------------------------------------------
#
# Each chart is checked at events in thirty digits, independently of the printed tensors: the field
# equations of the bulk on both sides of the wall at y = 0 where the chart crosses it, the delta on
# that wall against the boundary terms, and each chart against the static chart it reduces to or
# is pulled back from, with b_0 = a_0^6, the normalisation of V in Lehners, McFadden and Turok's
# charts, where h = a_0^5 H^{5/2} and r^2 = a_0^2 H.

VALUES = [{"alpha": sp.Rational(7, 10), "a_0": sp.Rational(13, 10), "b_0": sp.Rational(9, 10),
           "c_0": sp.Rational(11, 10), "rho": sp.Rational(3, 2), "A": sp.Rational(4, 5), "B": sp.Rational(6, 5),
           "t_0": sp.Rational(1, 10), "w_0": sp.Rational(1, 2), "h_0": sp.Rational(11, 10),
           "lambda": sp.Rational(2, 5)},
          {"alpha": sp.Rational(-1, 4), "a_0": sp.Rational(7, 10), "b_0": sp.Rational(3, 2),
           "c_0": sp.Rational(17, 10), "rho": sp.Rational(1, 2), "A": sp.Rational(3, 2), "B": sp.Rational(2, 3),
           "t_0": sp.Rational(-1, 3), "w_0": sp.Rational(3, 10), "h_0": sp.Rational(7, 10),
           "lambda": sp.Rational(9, 10)}]
ACROSS = [sp.Rational(2, 5), sp.Rational(-3, 4)]      # y on both sides of the wall at y = 0, |y| < pi rho
INSIDE = [sp.Rational(1, 5), sp.Rational(-1, 10)]     # w, |w| < w_0
RADII = [sp.Rational(3, 5), sp.Rational(1, 1)]        # r < 1/lambda for both cases
LIGHT = sp.Rational(3, 2)                              # c, so that a factor of c that is lost is caught
BRANCHES = [(sp.Rational(3, 11) * (1 - s * 4 / (3 * sp.sqrt(3))), sp.Rational(2, 11) * (1 + s * 2 * sp.sqrt(3)))
            for s in (1, -1)]


def _zero(value, what):
    if abs(complex(sp.N(value, DIGITS))) > TOLERANCE:
        raise AssertionError(f"{METRIC}: {what} fails, by {complex(sp.N(value, DIGITS))}")


def _explicit(chart, value):
    """A value with the held names of the rolling chart written out."""
    value = sp.sympify(value)
    return value.xreplace(chart.reader.held).doit() if chart.reader.held else value


def _volume(chart, system):
    """V in the chart's own symbols, the time standing as x^0 = ct."""
    X, P = chart.symbols, chart.reader.parameters
    if system == "static":
        return P["b_0"] * P["H"] ** 3
    if system == "conformal":
        return P["h"] ** sp.Rational(6, 5)
    if system == "static_bulk":
        return X[4] ** 6
    return _explicit(chart, P["b_0"] * P["b"] * P["H"] ** 3)


def _events(chart, system, value):
    X, P = chart.symbols, chart.reader.parameters
    at = {P[name]: v for name, v in value.items() if name in P and P[name].is_Symbol}
    at[chart.reader.c] = LIGHT
    at.update({X[1]: sp.Rational(1, 7), X[2]: sp.Rational(-2, 9), X[3]: sp.Rational(5, 11)})
    at[X[0]] = LIGHT * (value["t_0"] + sp.Rational(7, 10)) if system == "rolling" else sp.Rational(1, 3)
    if system in ("static", "rolling"):
        return [{**at, X[4]: v * value["rho"]} for v in ACROSS]
    if system == "conformal":
        return [{**at, X[4]: v * value["w_0"]} for v in INSIDE]
    return [{**at, X[4]: v} for v in RADII]


def _field_equations(chart, system, at, what=""):
    """The Einstein and scalar equations of the bulk at an event, a delta off its support being zero."""
    X = chart.symbols
    alpha = chart.reader.parameters["alpha"]
    g = sp.Matrix(chart.geo.g).applyfunc(lambda e: _explicit(chart, e))
    ginv = sp.Matrix(chart.geo.ginv).applyfunc(lambda e: _explicit(chart, e))
    ricci = sp.Matrix(chart.geo.ricci_ll()).applyfunc(lambda e: _explicit(chart, e))
    phi = sp.log(_volume(chart, system))
    grad = [sp.diff(phi, x) for x in X]
    potential = alpha ** 2 * sp.exp(-2 * phi)
    for a in range(5):
        for b in range(a, 5):
            _zero((ricci[a, b] - grad[a] * grad[b] / 2 - potential * g[a, b] / 9).subs(at),
                  f"{system}{what}: the Einstein equation in slot {chart.coords_tex[a]}{chart.coords_tex[b]} at {at}")
    root = sp.sqrt(-g.det())
    box = sum(sp.diff(root * ginv[a, a] * grad[a], X[a]) for a in range(5)) / root
    _zero((box + 2 * potential / 3).subs(at), f"{system}{what}: the equation of the scalar at {at}")


def _wall(chart, system, at):
    """The delta on the wall at y = 0 against the boundary terms of the action, taken just off it."""
    X = chart.symbols
    y = X[4]
    alpha = chart.reader.parameters["alpha"]
    D = sp.Symbol("_D")
    g = sp.Matrix(chart.geo.g).applyfunc(lambda e: _explicit(chart, e))
    ginv = sp.Matrix(chart.geo.ginv).applyfunc(lambda e: _explicit(chart, e))
    einstein = chart.geo.einstein_ll()
    import verify_metrics as vm
    V = _volume(chart, system)
    phi = sp.log(V)
    root = sp.sqrt(-g.det())
    box = sum(sp.diff(root * ginv[a, a] * sp.diff(phi, X[a]), X[a]) for a in range(5)) / root
    near = {**at, y: sp.Rational(1, 10 ** 40)}
    for a in range(5):
        value = _explicit(chart, vm._at(einstein, (a, a))).xreplace({sp.DiracDelta(y): D})
        on = sp.diff(value, D).subs(D, 0)
        want = sp.sqrt(2) * alpha * g[a, a] / (V * sp.sqrt(g[4, 4])) if a < 4 else 0
        _zero((on - want).subs(near), f"{system}: the wall's delta in G_{chart.coords_tex[a]}{chart.coords_tex[a]}")
    on = sp.diff(sp.expand(box).xreplace({sp.DiracDelta(y): D}), D).subs(D, 0)
    _zero((on - 2 * sp.sqrt(2) * alpha / (V * sp.sqrt(g[4, 4]))).subs(near), f"{system}: the wall's delta in box phi")


def _static_chart():
    import chart_printer as cp
    spec = chart_spec("static")
    return cp.Chart(spec["system"]["coords"], spec["system"]["parameters"], spec["chart_line_element"])


SIDE = [sp.Rational(1, 5), sp.Rational(1, 2)]          # y on the side y > 0, with H > 0 for both cases


def _against_static(chart, system, image, names, value, what):
    """This chart's metric pulled back along `image`, a function of the static chart's symbols
    giving this chart's coordinate across the interval on the side y > 0, against the static
    chart's metric with `names` putting its parameters in this chart's."""
    source = _static_chart()
    Y = source.symbols
    X = chart.symbols
    across = image(Y[4])
    J = sp.diag(1, 1, 1, 1, sp.diff(across, Y[4]))
    pulled = J.T * sp.Matrix(chart.geo.g).subs(dict(zip(X, list(Y[:4]) + [across])), simultaneous=True) * J
    theirs = sp.Matrix(source.geo.g).subs(names, simultaneous=True)
    at = {chart.reader.parameters[name]: v for name, v in value.items()
          if name in chart.reader.parameters and chart.reader.parameters[name].is_Symbol}
    at.update({Y[0]: sp.Rational(1, 3), Y[1]: sp.Rational(1, 7), Y[2]: sp.Rational(-2, 9), Y[3]: sp.Rational(5, 11)})
    for side in SIDE:
        e = {**at, Y[4]: side}
        for a in range(5):
            for b in range(a, 5):
                _zero((pulled[a, b] - theirs[a, b]).subs(e), f"{system}: {what}, slot {a}{b}")


def chart_check(chart, system):
    """Every chart against the field equations of the bulk, the kinked charts against the
    boundary terms on the wall at y = 0, and each other chart against the static one.

    static       the bulk on both sides of y = 0, and the delta on the wall;
    rolling      the bulk at both of Lukas, Ovrut and Waldram's branches (p, q), and the delta on
                 the wall at every time, which holds because the interval and V scale alike; at
                 p = 1/3, away from both branches, the bulk fails, so the check can fail;
    conformal    the bulk, and the static chart with b_0 = a_0^6 pulled back along
                 w = 2(a_0^5 H^{5/2} - h_0)/(5k) - w_0, k = sqrt(2) alpha/3, where h = a_0^5 H^{5/2};
    static_bulk  the bulk, and at lambda = 0 the static chart with b_0 = a_0^6 pulled back along
                 r = a_0 H^{1/2}, Lehners, McFadden and Turok's recovery of the static wall at mu = 0."""
    X, P = chart.symbols, chart.reader.parameters
    a0 = sp.Rational(6, 5)
    for value in VALUES:
        events = _events(chart, system, value)
        if system == "rolling":
            for p, q in BRANCHES:
                for e in events:
                    e = {**e, P["p"]: p, P["q"]: q}
                    _field_equations(chart, system, e, f" at p = {p}")
                    if e[X[4]] > 0:
                        _wall(chart, system, e)
            try:
                _field_equations(chart, system, {**events[0], P["p"]: sp.Rational(1, 3), P["q"]: BRANCHES[0][1]})
            except AssertionError:
                pass
            else:
                raise AssertionError(f"{METRIC}: the rolling chart solves the bulk away from its branches")
            continue
        for e in events:
            _field_equations(chart, system, e)
            if system == "static" and e[X[4]] > 0:
                _wall(chart, system, e)
        if system == "static":
            continue
        Q = _static_chart().reader.parameters
        alpha = P["alpha"]
        k = sp.sqrt(2) * alpha / 3
        names = {Q["alpha"]: alpha, Q["a_0"]: a0, Q["b_0"]: a0 ** 6}
        if system == "conformal":
            c0 = value["c_0"]
            names[Q["c_0"]] = c0
            _against_static(chart, system,
                            lambda s: 2 * (a0 ** 5 * (k * s + c0) ** sp.Rational(5, 2) - P["h_0"]) / (5 * k) - P["w_0"],
                            names, value, "the static chart with b_0 = a_0^6")
        else:
            c0 = value["c_0"]
            names[Q["c_0"]] = c0
            _against_static(chart, system,
                            lambda s: a0 * sp.sqrt(k * s + c0), names, {**value, "lambda": 0},
                            "the static chart at lambda = 0")
