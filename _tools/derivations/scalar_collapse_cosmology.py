"""The four charts of the scalar field universe of Husain, Martinez and Nunez and of Fonarev's family.

Each is a factor of time on Fisher's static bracket, with f = 1 - b/r:

conformal_time     ds^2 = (eta/eta_0) [-f^gamma d eta^2 + dr^2/f^gamma + f^(1-gamma) r^2 dOmega^2],
                   equations (3) and (4) of Husain, Martinez and Nunez (1994), with their a t + b
                   written eta/eta_0, their 2c written b and their alpha written gamma;
comoving_time      ds^2 = -f^gamma c^2 dt^2 + a^2 [dr^2/f^gamma + f^(1-gamma) r^2 dOmega^2] with
                   a = (t/t_0)^(1/3), the solution in comoving time of Faraoni (2013), the equation
                   after equation (3.161) in the arXiv text, with c dt = a d eta;
fonarev_conformal  ds^2 = a^2 [-f^gamma d eta^2 + ...] with a = (eta/eta_0)^s, Fonarev's family as
                   Kastor and Traschen (2017) write it in their (12) and Maeda (2007) in equations (6)
                   to (11), whose 2/(lambda^2 - 2) is s;
fonarev            ds^2 = -e^(2nkT) f^gamma c^2 dT^2 + e^(2kT) [...], equation (2.18) of Fonarev
                   (1995) with Fonarev's a written k/c, 2m written b, delta written gamma and
                   4 alpha^2 written n.

Every value is printed with gamma, s and n free, so it holds for any of them. The field equations
hold where gamma^2 = (s + 1)/(2s + 1) and n = gamma^2/(1 - gamma^2), and with no potential at
gamma = +-sqrt(3)/2, where s = 1/2 and n = 3; chart_check tests each of these at points, and each
chart after the first against the chart before it pulled back. scalar_collapse_cosmology.md beside
this file records each chart's source.
"""
import sympy as sp

METRIC = "scalar_collapse_cosmology"
CHARTS = ("conformal_time", "comoving_time", "fonarev_conformal", "fonarev")
ANGLES = ["\\theta \\in [0, \\pi]", "\\phi \\in [0, 2\\pi)"]
SPHERE = "\\left(d\\theta^2 + \\sin^2\\theta\\,d\\phi^2\\right)"
F = "f = 1 - \\dfrac{b}{r}"
SPACE = "\\dfrac{dr^2}{f^{\\gamma}} + f^{1-\\gamma}r^2" + SPHERE
NAMES = {"conformal_time": "Conformal Time", "comoving_time": "Comoving Time",
         "fonarev_conformal": "Fonarev Conformal Time", "fonarev": "Fonarev Time"}
CLOCK = {"conformal_time": "\\eta", "comoving_time": "t", "fonarev_conformal": "\\eta", "fonarev": "T"}
COORDS = {system: [CLOCK[system], "r", "\\theta", "\\phi"] for system in CHARTS}
SCALE = {"comoving_time": "a = \\left(\\dfrac{t}{t_0}\\right)^{1/3}",
         "fonarev_conformal": "a = \\left(\\dfrac{\\eta}{\\eta_0}\\right)^s"}
PARAMETERS = {
    "conformal_time": ["b", "\\eta_0", "\\gamma", F],
    "comoving_time": ["b", "t_0", "\\gamma", F, SCALE["comoving_time"]],
    "fonarev_conformal": ["b", "\\eta_0", "\\gamma", "s", F, SCALE["fonarev_conformal"]],
    "fonarev": ["b", "k", "\\gamma", "n", F],
}
# The charts whose time is a time, which the chart x^0 multiplies by c; the conformal time is a length.
TIMES = {"comoving_time": "t", "fonarev": "T"}

DIMENSIONS = {
    (METRIC, "conformal_time"): {"\\eta": "L", "r": "L", "\\theta": "1", "\\phi": "1",
                                 "b": "L", "\\eta_0": "L", "\\gamma": "1", "f": "1"},
    (METRIC, "comoving_time"): {"t": "T", "r": "L", "\\theta": "1", "\\phi": "1",
                                "b": "L", "t_0": "T", "\\gamma": "1", "f": "1", "a": "1"},
    (METRIC, "fonarev_conformal"): {"\\eta": "L", "r": "L", "\\theta": "1", "\\phi": "1",
                                    "b": "L", "\\eta_0": "L", "\\gamma": "1", "s": "1", "f": "1", "a": "1"},
    (METRIC, "fonarev"): {"T": "T", "r": "L", "\\theta": "1", "\\phi": "1",
                          "b": "L", "k": "1/T", "\\gamma": "1", "n": "1", "f": "1"},
}
# A scale factor that is a power of the time is held as a function with its slope declared, so
# that every value is a rational function of it and no root or symbolic power of the time is taken.
HELD = {(METRIC, "comoving_time"): ("a",), (METRIC, "fonarev_conformal"): ("a",)}
RATES = {(METRIC, "comoving_time"): {"a": {"t": "a/(3t)"}},
         (METRIC, "fonarev_conformal"): {"a": {"\\eta": "sa/\\eta"}}}

# The exponent of the massless field, and the points every check is made at: a value of the factor
# of time's argument, a distance from r = b, and the two lengths.
ROOT = sp.sqrt(3) / 2
POINTS = ((sp.Rational(7, 5), sp.Rational(23, 10), sp.Rational(9, 10), sp.Rational(6, 5)),
          (sp.Rational(1, 3), sp.Rational(31, 10), sp.Rational(13, 10), sp.Rational(3, 4)),
          (sp.Rational(29, 7), sp.Rational(11, 10), sp.Rational(1, 2), sp.Rational(17, 9)))
LIGHT = sp.Rational(13, 10)


def register(pc):
    """Hand the printer its charts; pc is print_charts itself."""
    pc.CHARTS[METRIC] = [lambda system=system: chart_spec(pc, system) for system in CHARTS]


def line_element(system, c2="c^2"):
    """The line element as the file prints it, or with c2 empty as the chart x^0 has it."""
    if system == "conformal_time":
        return "ds^2 = \\dfrac{\\eta}{\\eta_0}\\left(-f^{\\gamma}d\\eta^2 + " + SPACE + "\\right)"
    if system == "comoving_time":
        return f"ds^2 = -f^{{\\gamma}}{c2}dt^2 + a^2\\left(" + SPACE + "\\right)"
    if system == "fonarev_conformal":
        return "ds^2 = a^2\\left(-f^{\\gamma}d\\eta^2 + " + SPACE + "\\right)"
    over = "" if c2 else "/c"
    return f"ds^2 = -e^{{2nkT{over}}}f^{{\\gamma}}{c2}dT^2 + e^{{2kT{over}}}\\left(" + SPACE + "\\right)"


def domains(system):
    centre = ["r \\in (b, \\infty)"] + ANGLES
    edge = "r = b \\;\\text{(curvature singularity)}"
    if system == "comoving_time":
        return (["t \\in (0, \\infty) \\;\\text{for}\\; t_0 > 0", "t \\in (-\\infty, 0) \\;\\text{for}\\; t_0 < 0"] + centre
                + ["t = 0 \\;\\text{(curvature singularity, the big bang for}\\; t_0 > 0\\text{)}", edge])
    if system == "fonarev":
        return (["T \\in (-\\infty, \\infty)"] + centre
                + ["kT \\to -\\infty \\;\\text{(curvature singularity, the big bang for}\\; k > 0\\text{)}", edge])
    return (["\\eta \\in (0, \\infty) \\;\\text{for}\\; \\eta_0 > 0", "\\eta \\in (-\\infty, 0) \\;\\text{for}\\; \\eta_0 < 0"]
            + centre + ["\\eta = 0 \\;\\text{(curvature singularity, the big bang for}\\; \\eta_0 > 0\\text{)}", edge])


def chart_parameters(system):
    """The parameters as the chart x^0 = ct reads them: a scale factor named in the time holds c t_0."""
    return [p.replace("{t_0}", "{ct_0}") for p in PARAMETERS[system]]


def by_hand(system):
    """The metric components the line element writes its own way, each read back against the chart's."""
    clock = CLOCK[system]
    lapse, scale = {"conformal_time": ("\\dfrac{\\eta}{\\eta_0}", "\\dfrac{\\eta}{\\eta_0}"),
                    "comoving_time": ("", "a^2"),
                    "fonarev_conformal": ("a^2", "a^2"),
                    "fonarev": ("e^{2nkT}", "e^{2kT}")}[system]
    return {"metric_components": {(clock, clock): "-" + lapse + "f^{\\gamma}",
                                  ("\\theta", "\\theta"): scale + "r^2f^{1-\\gamma}",
                                  ("\\phi", "\\phi"): scale + "r^2f^{1-\\gamma}\\sin^2\\theta"}}


def chart_spec(pc, system):
    vm = pc.vm
    state = {}
    coords, parameters = COORDS[system], chart_parameters(system)
    key = (METRIC, system)
    probe = vm.Reader(coords, parameters, (), held=HELD.get(key, ()), rates=RATES.get(key))
    r, b, gamma = probe.symbol["r"], probe.parameters["b"], probe.parameters["gamma"]
    clock = probe.symbol[coords[0]]
    whole = sp.Symbol("HMNsum")
    named = {whole: "2r - b\\left(1 + \\gamma\\right)"}
    wholes = [(whole, 2 * r - (1 + gamma) * b)]
    powers = pc.named_powers({"f": {r - b: 1, r: -1}}, state)
    table = {}
    # The two sums an exponential of Fonarev's time holds, each written as one factor of its exponent.
    rates = []
    if "n" in probe.parameters:
        n = probe.parameters["n"]
        rates = [(sp.Symbol("HMNless"), 1 - n), (sp.Symbol("HMNmore"), 1 + n)]
        named.update({rates[0][0]: "1 - n", rates[1][0]: "1 + n"})

    def exponential(value):
        """The exponentials of a product as one, its exponent factored around 1 - n or 1 + n."""
        value = sp.powsimp(value)
        if not isinstance(value, sp.exp):
            return value
        out = sp.Integer(1)
        for factor in sp.Mul.make_args(sp.factor(value.args[0])):
            for placeholder, poly in rates:
                if sp.expand(factor - poly) == 0:
                    factor = placeholder
                elif sp.expand(factor + poly) == 0:
                    factor = -placeholder
            out *= factor
        return sp.exp(out)

    def collected(base):
        """A polynomial in r as a named sum grouped by its powers of r, as Fisher's chart writes
        its own, leading with a term that has no minus sign; and the sign taken out of it."""
        printer = state["printer"]
        terms = list(pc.cp.collect_by(base, [r], printer).terms)
        sign = -1 if all(printer.leads_negative(pair) for pair in terms) else 1
        terms = [(sign * c, rest) for c, rest in terms]
        for i, pair in enumerate(terms):
            if not printer.leads_negative(pair):
                terms = [pair] + terms[:i] + terms[i + 1:]
                break
        text = printer.sum_text(pc.cp.Sum(terms))
        placeholder = table.setdefault(text, sp.Symbol(f"HMNpoly{len(table)}"))
        printer.named[placeholder] = text
        return sign, placeholder

    def product(value):
        """One product of powers, written around f^gamma, with the sum the horizon holds written whole."""
        # An exponential of the time is no power of f, so it is set aside and put back as one factor.
        factors = sp.Mul.make_args(value)
        out = exponential(sp.Mul(*[factor for factor in factors if factor.has(sp.exp)]))
        value = sp.Mul(*[factor for factor in factors if not factor.has(sp.exp)])
        for factor in sp.Mul.make_args(powers(value)):
            base, e = (factor.base, factor.exp) if factor.is_Pow else (factor, sp.Integer(1))
            if base.is_Add and e.is_Integer:
                content, primitive = base.as_content_primitive()
                for placeholder, poly in wholes:
                    if sp.expand(primitive - poly) == 0:
                        base, out = placeholder, out * content ** e
                    elif sp.expand(primitive + poly) == 0:
                        base, out = placeholder, out * (-content) ** e
                if base.is_Add and base.has(r) and len(base.args) > 2:
                    sign, base = collected(base)
                    out *= sign ** e
            out *= base ** e
        return out

    def pretty(value):
        # The static bracket gives each value one power of f^gamma and the factor of time another,
        # so a value is written as a sum, one term for each power of f^gamma it holds.
        return sp.Add(*[product(part) for part in by_power(value, r, b, gamma)])

    def collect(poly, printer):
        # A value that is a sum of parts leads with a part that is written with no minus sign;
        # a polynomial inside a part keeps the printer's own order.
        if not any(str(x).startswith("NAMEDPOWER") for x in poly.free_symbols) and not poly.has(sp.exp):
            return printer.sum_of(poly)
        terms = list(printer.sum_of(poly).terms)
        for i, pair in enumerate(terms):
            if not printer.leads_negative(pair):
                terms = [pair] + terms[:i] + terms[i + 1:]
                break
        return pc.cp.Sum(terms)

    def check(chart):
        state["printer"], state["chart"] = chart.printer, chart
        chart_check(pc, chart, system)

    constants = [probe.parameters[name] for name in ("n", "s", "k", "eta_0", "t_0") if name in probe.parameters]
    scale = [probe.functions[name][0] for name in ("a",) if name in probe.functions]
    spec = {
        "metric_id": METRIC,
        "system": {"id": system, "name": NAMES[system], "coords": coords, "domains": domains(system),
                   "parameters": PARAMETERS[system], "line_element": line_element(system)},
        "chart_parameters": parameters,
        "chart_line_element": line_element(system, ""),
        "printer": {"lead": [r, clock, b, gamma], "factors": [rate[0] for rate in rates] + [b] + constants + [gamma] + scale + [probe.c, clock, r],
                    "rising": [gamma], "flip": False, "named": named, "collect": collect},
        "pretty": pretty,
        # A sum that stands in a bracket is handed over with the chart's x^0 standing for the time.
        "bracketed": lambda value: pretty(state["chart"].named_time(value)),
        "components": by_hand(system),
        "check": check,
    }
    if system in TIMES:
        spec["time"] = TIMES[system]
    if key in RATES:
        spec["reduce"] = probe.surface
    return spec


def by_power(value, r, b, gamma):
    """The parts of a value, one for each power of f^gamma = ((r - b)/r)^gamma it holds, each a
    rational function times that power, written (r - b)^(k gamma) r^(-k gamma)."""
    Y, G = sp.Symbol("HMNpowerY", positive=True), sp.Symbol("HMNpowerG", positive=True)

    def named(power):
        base, exponent = power.base, power.exp
        while base.is_Pow:
            # (1/r)^gamma is r^(-gamma).
            base, exponent = base.base, exponent * base.exp
        whole, part = sp.expand(exponent).as_coeff_Add()
        k = sp.cancel(part / gamma)
        if not k.is_Rational:
            raise AssertionError(f"{power} is no power of f^gamma")
        if sp.expand(base - r + b) == 0:
            return (r - b) ** whole * (Y * G) ** k
        if base == r:
            return r ** whole * G ** k
        if base == -1:
            # The reader's phase of a base written with the other sign: the ratio is positive.
            return sp.Integer(-1) ** whole
        raise AssertionError(f"{power} is a power of neither r - b nor r")

    value = sp.expand_power_base(sp.powdenest(sp.together(sp.sympify(value)), force=True), force=True)
    value = value.replace(lambda e: e.is_Pow and not e.exp.is_Number and e.exp.has(gamma), named)
    numerator, denominator = sp.fraction(sp.together(value))
    below = sp.factor(denominator).as_powers_dict()
    rest = sp.Mul(*[base ** e for base, e in below.items() if base not in (Y, G)])
    if rest.has(Y, G):
        raise AssertionError(f"{value} has a power of f^gamma inside a sum below the line")
    parts = {}
    for term in sp.Add.make_args(sp.expand(numerator)):
        held = term.as_powers_dict()
        if held.get(G, 0) != below.get(G, 0):
            raise AssertionError(f"{term} holds a power of r that is no power of f^gamma")
        k = held.get(Y, 0) - below.get(Y, 0)
        parts[k] = parts.get(k, 0) + term / (Y ** held.get(Y, 0) * G ** held.get(G, 0))
    return [sp.factor(part / rest) * (r - b) ** (k * gamma) * r ** (-k * gamma)
            for k, part in sorted(parts.items(), key=lambda item: -item[0]) if sp.expand(part) != 0]


def solution(system, X, P, c):
    """What the chart is on Fonarev's solutions, in the chart x^0: the scaled field phi, with
    R_ab = 2 d_a phi d_b phi + W g_ab, the potential term W, and the relation that puts the
    chart's parameters on the solutions, each as a function of the exponent gamma."""
    clock, r = X[:2]
    gamma = P["gamma"]
    f = 1 - P["b"] / r
    width = sp.sqrt(1 - gamma ** 2)
    s = (1 - gamma ** 2) / (2 * gamma ** 2 - 1)
    n = gamma ** 2 / (1 - gamma ** 2)
    if system == "conformal_time":
        # The massless member alone, in the article's own form: either sign of the field for either exponent.
        phi = sp.sqrt(3) * sp.log(clock / P["eta_0"]) / 2 + gamma * sp.log(f) / (2 * sp.sqrt(3))
        return phi, sp.Integer(0), {}
    if system == "comoving_time":
        phi = sp.log(clock / (c * P["t_0"])) / sp.sqrt(3) + gamma * sp.log(f) / (2 * sp.sqrt(3))
        return phi, sp.Integer(0), {}
    if system == "fonarev_conformal":
        phi = gamma * s * sp.log(clock / P["eta_0"]) / width + width * sp.log(f) / 2
        W = s * (2 * s - 1) * f ** (-gamma) * (clock / P["eta_0"]) ** (-2 * s) / clock ** 2
        return phi, W, {P["s"]: s}
    rate = P["k"] / c
    phi = gamma * rate * clock / width + width * sp.log(f) / 2
    W = (3 - n) * rate ** 2 * f ** (-gamma) * sp.exp(-2 * n * rate * clock)
    return phi, W, {P["n"]: n}


def earlier(system, X, P, c):
    """The chart before this one that it is pulled back from: that chart's own time as a function
    of this chart's x^0, and that chart's parameters in this one's."""
    clock = X[0]
    if system == "comoving_time":
        # c dt = a d eta with a = (eta/eta_0)^(1/2) = (t/t_0)^(1/3) and t_0 = 2 eta_0/3c.
        unit = 3 * c * P["t_0"] / 2
        return "conformal_time", unit * (clock / (c * P["t_0"])) ** sp.Rational(2, 3), {"eta_0": unit}
    if system == "fonarev_conformal":
        # The family holds the first chart at s = 1/2.
        return "conformal_time", clock, {"eta_0": P["eta_0"]}
    # e^(nkT) c dT = e^(kT) d eta, so eta = eta_0 e^(kT/s) with s = 1/(n - 1) and eta_0 = sc/k.
    s = 1 / (P["n"] - 1)
    unit = s * c / P["k"]
    return "fonarev_conformal", unit * sp.exp(P["k"] * clock / (c * s)), {"eta_0": unit, "s": s}


def chart_check(pc, chart, system):
    """At three points, in forty digits: on the solutions the Ricci tensor is 2 d phi d phi + W g
    and the field solves its wave equation with that potential, for both signs of the massless
    exponent and for two more exponents where the chart holds the family; off the solutions the
    first fails, so the check can fail; the first chart has the Ricci scalar Faraoni states and a
    null gradient of the areal radius on the apparent horizon of the article's (11); and each
    chart after the first is the chart before it pulled back. No value is simplified."""
    vm, cp = pc.vm, pc.cp
    X, P, c = chart.symbols, chart.reader.parameters, chart.reader.c
    clock, r = X[:2]
    b, gamma = P["b"], P["gamma"]
    held = dict(chart.reader.held)
    family = system in ("fonarev_conformal", "fonarev")
    exponents = (ROOT, -ROOT) + ((sp.Rational(4, 5), -sp.Rational(3, 5)) if family else ())
    unit = [P[name] for name in ("eta_0", "t_0", "k") if name in P][0]
    g, ginv = chart.geo.g, chart.geo.ginv
    ricci = chart.geo.ricci_ll()
    phi, W, on = solution(system, X, P, c)
    f = 1 - b / r

    def number(value, at):
        value = sp.sympify(value).subs(held).doit().subs(on).subs(held).doit()
        return complex(sp.N(value.subs(at), 40))

    def near(values, at, label):
        worst = max((abs(number(v, at)) for v in values), default=0)
        if worst > 1e-25:
            raise AssertionError(f"{METRIC}/{system}: {label} fails by {worst:.3g}")

    def at_point(point, exponent):
        tv, rv, bv, uv = point
        return {clock: tv * (LIGHT if system in TIMES else 1), r: rv + bv, b: bv, unit: uv, gamma: exponent,
                X[2]: sp.Rational(4, 5), X[3]: sp.Rational(1, 3), c: LIGHT}

    root = sp.sqrt(-g.det())
    residual = [vm._at(ricci, (i, j)) - 2 * sp.diff(phi, X[i]) * sp.diff(phi, X[j]) - W * g[i, j]
                for i, j in vm._indices(4, 2)]
    wave = (sum(sp.diff(root * ginv[i, j] * sp.diff(phi, X[j]), X[i]) for i in range(4) for j in range(4)) / root
            - sp.diff(W, X[1]) / (2 * sp.diff(phi, X[1])))
    for exponent in exponents:
        for point in POINTS:
            at = at_point(point, exponent)
            near(residual, at, f"R_ab = 2 d_a phi d_b phi + W g_ab at gamma = {exponent}")
            # W depends on the place through the field alone, so its slope along the field is
            # its slope along r over the field's; the wave equation is box phi = W'(phi)/2.
            near([wave], at, f"the wave equation at gamma = {exponent}")
    stray = at_point(POINTS[0], sp.Rational(7, 10) if family else sp.Rational(4, 5))
    if family:
        stray[[P[name] for name in ("s", "n") if name in P][0]] = sp.Rational(5, 3)
        on = {}
    if max(abs(number(v, stray)) for v in residual) < 1e-6:
        raise AssertionError(f"{METRIC}/{system}: the field equations hold off the solutions, where they should fail")
    phi, W, on = solution(system, X, P, c)

    if system == "conformal_time":
        eta0 = P["eta_0"]
        scalar = chart.geo.ricci_scalar()
        area = sp.sqrt(g[2, 2])
        gradient = sum(ginv[i, j] * sp.diff(area, X[i]) * sp.diff(area, X[j]) for i in range(2) for j in range(2))
        horizon = r ** 2 * f ** (1 - gamma) / (2 * (r - b * (1 + gamma) / 2))
        for exponent in (ROOT, -ROOT):
            for point in POINTS:
                at = at_point(point, exponent)
                near([scalar - (gamma ** 2 * b ** 2 * f ** (gamma - 2) * eta0 / (6 * r ** 4 * clock)
                                - 3 * eta0 / (2 * clock ** 3 * f ** gamma))], at, "Faraoni's Ricci scalar")
                on_horizon = dict(at)
                on_horizon[clock] = horizon.subs({r: at[r], b: at[b], gamma: exponent})
                near([gradient], on_horizon, "the apparent horizon")
        return

    source, time, renamed = earlier(system, X, P, c)
    own = cp.Chart(COORDS[source], chart_parameters(source), line_element(source, ""), None, None, None, None, None,
                   None, HELD.get((METRIC, source), ()))
    names = dict(zip(own.symbols[1:], X[1:]))
    names.update({own.reader.parameters[name]: P[name] for name in ("b", "gamma")})
    names.update({own.reader.parameters[name]: value for name, value in renamed.items()})
    names[own.reader.c] = c
    J = sp.diag(sp.diff(time, clock), 1, 1, 1)
    before = own.geo.g.subs(own.reader.held).doit().subs(names, simultaneous=True).subs(own.symbols[0], time)
    pulled = J.T * before * J
    # The first chart is the family's massless member alone, and Fonarev's chart is the family's for any exponent.
    for exponent in (exponents if system == "fonarev" else (ROOT, -ROOT)):
        for point in POINTS:
            at = at_point(point, exponent)
            near(list(g - pulled), at, f"the {source} chart pulled back at gamma = {exponent}")
