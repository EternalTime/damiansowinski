"""The two charts of the pre-big bang dilaton vacuum.

synchronous  ds^2 = -dt^2 + sum_i (-t)^{2 beta_i} dx_i^2 on t < 0, with the dilaton
             phi = (sum_i beta_i - 1) ln(-t) and sum_i beta_i^2 = 1, equation (3.1) of
             Gasperini and Veneziano (2003) and the four dimensional case of Veneziano's
             (1991) equation (4.9);
conformal    the two members with beta_1 = beta_2 = beta_3 = beta = -+1/sqrt3 in the
             conformal time a d eta = dt of Veneziano's (2.3), whose scale factor is the
             power |eta|^k of Gasperini and Veneziano's (2.94), k = beta/(1 - beta).

Every value of the synchronous chart assumes sum_i beta_i^2 = 1: check_synchronous holds
the Ricci tensor and the two scalars stated there to the tensors of free exponents on the
checker's parametrisation of that sphere, and the stated Ricci tensor is minus the second
covariant derivative of the dilaton, which is the first field equation. The conformal
chart's values hold for every k. pre_big_bang.md beside this file is the derivation.
"""
import itertools

import sympy as sp

METRIC = "pre_big_bang"
CHARTS = ("synchronous", "conformal")
SPACE = ["x \\in (-\\infty, \\infty)", "y \\in (-\\infty, \\infty)", "z \\in (-\\infty, \\infty)"]
BETAS = ("\\beta_1", "\\beta_2", "\\beta_3")

# Units with c = 1: the times and the comoving coordinates are lengths, and a power whose
# exponent holds a parameter is read with its base in a fixed unit.
DIMENSIONS = {
    (METRIC, "synchronous"): {"t": "L", "x": "L", "y": "L", "z": "L", **{b: "1" for b in BETAS}},
    (METRIC, "conformal"): {"\\eta": "L", "x": "L", "y": "L", "z": "L", "k": "1", "a": "1"},
}

# The sphere sum_i beta_i^2 = 1, projected stereographically from (0, 0, 1): every point of
# it but that one is reached, and that one is the flat member (0, 0, 1).
PARAMETER_RELATIONS = {
    (METRIC, "synchronous"): {
        "beta_1": "2*u/(1 + u**2 + v**2)",
        "beta_2": "2*v/(1 + u**2 + v**2)",
        "beta_3": "(u**2 + v**2 - 1)/(1 + u**2 + v**2)",
    },
}

# The two isotropic members: k = beta/(1 - beta) at beta = -1/sqrt3 and beta = 1/sqrt3, the two
# roots of 2k^2 - 2k - 1, which is 3 beta^2 = 1 written in k.
MEMBERS = {"expanding": -1 / (1 + sp.sqrt(3)), "contracting": 1 / (sp.sqrt(3) - 1)}
SCALE = "a = \\left(-\\eta/(1 + k)\\right)^k"
SPACETIME = "\\left(-d\\eta^2 + dx^2 + dy^2 + dz^2\\right)"


def register(pc):
    """Hand the printer its two charts; pc is print_charts itself."""
    pc.CHARTS[METRIC] = [lambda system=system: chart_spec(pc, system) for system in CHARTS]


def line_element(system):
    if system == "synchronous":
        return ("ds^2 = -dt^2 + (-t)^{2\\beta_1}dx^2 + (-t)^{2\\beta_2}dy^2"
                " + (-t)^{2\\beta_3}dz^2")
    return "ds^2 = a^2" + SPACETIME


def chart_line_element(system):
    """The line element with every named quantity written out, as the tensors are built from it."""
    if system == "synchronous":
        return line_element(system)
    return "ds^2 = \\left(-\\dfrac{\\eta}{1 + k}\\right)^{2k}" + SPACETIME


def chart_spec(pc, system):
    vm = pc.vm
    state = {}
    if system == "synchronous":
        coords, parameters = ["t", "x", "y", "z"], list(BETAS)
        probe = vm.Reader(coords, parameters, ())
        time = probe.symbol["t"]
        names = [probe.parameters[f"beta_{i}"] for i in (1, 2, 3)]

        def check(chart):
            state["printer"] = chart.printer
            check_synchronous(pc, chart)
        return {
            "metric_id": METRIC,
            "system": {"id": system, "name": "Synchronous", "coords": coords,
                       "domains": ["t \\in (-\\infty, 0)"] + SPACE,
                       "parameters": parameters, "line_element": line_element(system)},
            "chart_line_element": chart_line_element(system),
            "printer": {"lead": names + [time], "factors": names + [time]},
            "pretty": minus_powers(time, "t", state),
            "check": check,
            "components": {
                "metric_components": {(k, k): f"(-t)^{{2\\beta_{i}}}" for i, k in enumerate("xyz", 1)},
                "inverse_metric_components": {(k, k): f"(-t)^{{-2\\beta_{i}}}"
                                              for i, k in enumerate("xyz", 1)}},
            "kretschmann": KRETSCHMANN,
            # The printer leads a bare sum with its first positive term; the mixed R^t_t reads
            # better as minus R_tt.
            "rewrite": [("\\dfrac{1 - \\beta_3 - \\beta_2 - \\beta_1}{t^2}",
                         "-\\dfrac{\\beta_1 + \\beta_2 + \\beta_3 - 1}{t^2}")],
        }
    coords, parameters = ["\\eta", "x", "y", "z"], ["k", SCALE]
    probe = vm.Reader(coords, parameters, ())
    eta, k = probe.symbol["\\eta"], probe.parameters["k"]
    a = sp.Symbol("a", positive=True)

    def check(chart):
        state["printer"] = chart.printer
        check_conformal(pc, chart)
    return {
        "metric_id": METRIC,
        "system": {"id": system, "name": "Conformal Time", "coords": coords,
                   "domains": ["\\eta \\in (-\\infty, 0)"] + SPACE,
                   "parameters": parameters, "line_element": line_element(system)},
        "chart_line_element": chart_line_element(system),
        "printer": {"lead": [k, eta, a]},
        "pretty": scale_powers(eta, k, a),
        "check": check,
    }


# K = 4 (sum_i R^{ti}_{ti}^2 + sum_{i<j} R^{ij}_{ij}^2), with R^{ti}_{ti} = beta_i(beta_i - 1)/t^2
# and R^{ij}_{ij} = beta_i beta_j/t^2: a sum of squares, which holds for free exponents.
KRETSCHMANN = ("\\dfrac{4\\left(\\beta_1^2\\left(\\beta_1 - 1\\right)^2 + \\beta_2^2\\left(\\beta_2 - 1\\right)^2"
               " + \\beta_3^2\\left(\\beta_3 - 1\\right)^2 + \\beta_1^2\\,\\beta_2^2 + \\beta_1^2\\,\\beta_3^2"
               " + \\beta_2^2\\,\\beta_3^2\\right)}{t^4}")


def minus_powers(time, time_tex, state):
    """A pretty printer for a value that is one monomial in -t raised to an exponent holding the
    exponents, as (-t)^{2beta_1 - 1}. The checker writes (-t)^e as (-1)^e t^e, so the parts of e
    that hold a parameter are gathered back into one power of -t, with any integer power of t
    taken into it at its sign, and the power is handed to the printer as a placeholder whose text
    is -t with its whole exponent, above the line, or below it where every term of the exponent
    is negative. A value with no parameter in its exponent is factored and printed in t."""
    table = {}

    def pretty(value):
        value = sp.together(sp.sympify(value))
        value = sp.powsimp(sp.expand_power_base(sp.powdenest(value, force=True), force=True), combine="exp")
        of_time, of_sign, rest = sp.Integer(0), sp.Integer(0), sp.Integer(1)
        for f in sp.Mul.make_args(value):
            base, k = f.as_base_exp()
            if base == time:
                of_time += k
            elif base == -1 and not k.is_Number:
                of_sign += k
            else:
                rest *= f
        if rest.has(time):
            raise AssertionError(f"{value} is not one monomial in -{time}")
        of_time, of_sign = sp.expand(of_time), sp.expand(of_sign)
        whole, held = of_time.as_coeff_Add()
        sign_whole, sign_held = of_sign.as_coeff_Add()
        if sp.expand(held - sign_held) != 0:
            raise AssertionError(f"{value} is not a power of -{time}: the signs and the powers differ")
        # (-1)^(sign_whole) t^whole (-t)^held, with t^whole = (-1)^whole (-t)^whole.
        out = sp.factor(rest * sp.Integer(-1) ** (sign_whole + whole))
        exponent = whole + held
        if held == 0:
            return sp.factor(rest * sp.Integer(-1) ** sign_whole * time ** whole)
        printer = state["printer"]
        negative = all(term.as_coeff_Mul()[0].is_negative for term in sp.Add.make_args(exponent))
        shown, side = (-exponent, -1) if negative else (exponent, 1)
        text = "\\left(-" + time_tex + "\\right)^{" + printer.positive_first(shown) + "}"
        placeholder = table.setdefault(text, sp.Symbol(f"MINUSPOWER{len(table)}", positive=True))
        printer.overrides[placeholder] = text
        return out * placeholder ** side

    return pretty


def scale_powers(eta, k, a):
    """A pretty printer for the conformal chart, whose every value is a power of the scale factor
    a = (-eta/(1 + k))^k times a rational function of eta and k. The checker writes that power as
    (-1)^(nk) eta^(nk) (1 + k)^(-nk), so the three are gathered back into a^n, and the rest is
    factored."""
    scale = 1 + k

    def pretty(value):
        value = sp.together(sp.sympify(value))
        value = sp.powsimp(sp.expand_power_base(sp.powdenest(value, force=True), force=True), combine="exp")
        powers = {"eta": sp.Integer(0), "sign": sp.Integer(0), "scale": sp.Integer(0)}
        rest = sp.Integer(1)
        for f in sp.Mul.make_args(value):
            base, e = f.as_base_exp()
            while base.is_Pow:
                # (1/(1 + k))^(2k) is (1 + k)^(-2k).
                base, e = base.base, e * base.exp
            if base == eta:
                powers["eta"] += e
            elif base == -1 and not e.is_Number:
                powers["sign"] += e
            elif sp.expand(base - scale) == 0 and e.has(k):
                powers["scale"] += e
            else:
                rest *= f
        if rest.has(eta):
            raise AssertionError(f"{value} is not a power of the scale factor times a rational function")
        split = {name: sp.expand(e).as_coeff_Add() for name, e in powers.items()}
        held = split["eta"][1]
        if sp.expand(split["sign"][1] - held) != 0 or sp.expand(split["scale"][1] + held) != 0:
            raise AssertionError(f"{value} holds powers of -eta and 1 + k that are no power of a")
        n = sp.cancel(held / k)
        if not n.is_Integer:
            raise AssertionError(f"{value} holds a power of the scale factor that is not whole")
        out = (rest * sp.Integer(-1) ** split["sign"][0] * eta ** split["eta"][0]
               * scale ** split["scale"][0])
        return sp.factor(out) * a ** n

    return pretty


def on_sphere(pc, betas):
    """print_charts's OnSurface with each polynomial in the exponents that the Weyl tensor holds
    written in its shortest form modulo sum_i beta_i^2 = 1, over every order in which the
    exponents can be eliminated, so that a component holds eight terms where it would hold ten."""
    sphere = sum(b ** 2 for b in betas) - 1

    def shortest(polynomial):
        forms = {sp.reduced(polynomial, [sphere], *order, order="lex")[1]
                 for order in itertools.permutations(betas)} | {polynomial}
        return min(forms, key=lambda f: (sp.count_ops(f), sp.default_sort_key(f)))

    def reduced(value):
        out = sp.Integer(1)
        for f in sp.Mul.make_args(sp.factor(value)):
            base, e = f.as_base_exp()
            if base.free_symbols and base.free_symbols <= set(betas) and sp.Poly(base, *betas).total_degree() > 1:
                f = shortest(sp.expand(base)) ** e
            out *= f
        return sp.factor(out)

    class Sphere(pc.OnSurface):
        def weyl_llll(self):
            weyl = super().weyl_llll()
            return [[[[reduced(c) for c in row] for row in plane] for plane in block] for block in weyl]
    return Sphere


def check_synchronous(pc, chart):
    """The two field equations of the string frame on the sphere sum_i beta_i^2 = 1, with
    phi = (sum_i beta_i - 1) ln(-t): R_mu nu + nabla_mu nabla_nu phi = 0, Veneziano's (3.7) and
    Gasperini and Veneziano's (2.10), and R + 2 box phi - (nabla phi)^2 = 0. The Ricci tensor
    stated is then minus the second covariant derivative of the dilaton, the Ricci scalar
    (sum_i beta_i - 1)^2/t^2, and the Kretschmann scalar KRETSCHMANN, each held to the tensor of
    free exponents on the sphere; the stated curvature stands in for the chart's, so that it is
    what is printed."""
    vm = pc.vm
    geo, reader = chart.geo, chart.reader
    n = geo.n
    t = chart.symbols[0]
    betas = [reader.parameters[f"beta_{i}"] for i in (1, 2, 3)]
    surface = {reader.parameters[name]: sp.sympify(value)
               for name, value in PARAMETER_RELATIONS[(METRIC, "synchronous")].items()}

    def on(value):
        return vm.norm(sp.sympify(value).subs(surface))

    total = sum(betas)
    phi = (total - 1) * sp.log(-t)
    gamma = geo.christoffel_ull()
    gradient = [sp.diff(phi, x) for x in chart.symbols]
    hessian = [[vm.norm(sp.diff(phi, chart.symbols[a], chart.symbols[b])
                        - sum(gamma[c][a][b] * gradient[c] for c in range(n)))
                for b in range(n)] for a in range(n)]
    ricci = geo.ricci_ll()
    for a in range(n):
        for b in range(n):
            if on(ricci[a][b] + hessian[a][b]) != 0:
                raise AssertionError(f"{METRIC}: R + nabla nabla phi is not zero in slot {a}{b}")
    box = sum(geo.ginv[a, b] * hessian[a][b] for a in range(n) for b in range(n))
    square = sum(geo.ginv[a, b] * gradient[a] * gradient[b] for a in range(n) for b in range(n))
    if on(geo.ricci_scalar() + 2 * box - square) != 0:
        raise AssertionError(f"{METRIC}: the dilaton equation does not hold on the sphere")
    # Off the sphere the tt equation fails by -(sum beta_i^2 - 1)/t^2 and the dilaton equation
    # by (sum beta_i^2 - 1)/t^2, and the space equations hold for free exponents.
    off = (sum(b ** 2 for b in betas) - 1) / t ** 2
    if vm.norm(ricci[0][0] + hessian[0][0] + off) != 0 or vm.norm(geo.ricci_scalar() + 2 * box - square - off) != 0:
        raise AssertionError(f"{METRIC}: the field equations off the sphere are misstated")
    stated_ricci = [[-hessian[a][b] for b in range(n)] for a in range(n)]
    scalar = (total - 1) ** 2 / t ** 2
    kretschmann = chart.reader(KRETSCHMANN)
    if on(geo.ricci_scalar() - scalar) != 0:
        raise AssertionError(f"{METRIC}: the Ricci scalar is misstated")
    if vm.norm(geo.kretschmann() - kretschmann) != 0:
        raise AssertionError(f"{METRIC}: the Kretschmann scalar is misstated for free exponents")
    stated = on_sphere(pc, betas)(geo, stated_ricci, scalar, kretschmann)
    for name, rank in (("ricci_ll", 2), ("einstein_ll", 2), ("weyl_llll", 4)):
        free, there = getattr(geo, name)(), getattr(stated, name)()
        for index in vm._indices(n, rank):
            if on(vm._at(free, index) - vm._at(there, index)) != 0:
                raise AssertionError(f"{METRIC}: {name} is misstated in slot {index}")
    chart.geo = stated


def check_conformal(pc, chart):
    """The conformal chart is the synchronous one at beta_i = k/(1 + k) pulled back through
    -t = (-eta/(1 + k))^(1 + k), for every k, and with phi = (2k - 1) ln(-eta/(1 + k)) both field
    equations hold where 2k^2 - 2k - 1 = 0, at the two members and nowhere else."""
    vm = pc.vm
    geo = chart.geo
    n = geo.n
    eta = chart.symbols[0]
    k = chart.reader.parameters["k"]
    quadratic = 2 * k ** 2 - 2 * k - 1
    if any(sp.simplify(quadratic.subs(k, value)) != 0 for value in MEMBERS.values()):
        raise AssertionError(f"{METRIC}/conformal: the members are not the roots of 2k^2 - 2k - 1")
    if [sp.nsimplify(sp.simplify(v / (1 + v))) for v in MEMBERS.values()] != [-1 / sp.sqrt(3), 1 / sp.sqrt(3)]:
        raise AssertionError(f"{METRIC}/conformal: the members are not beta = -1/sqrt3 and 1/sqrt3")
    a = (-eta / (1 + k)) ** k
    minus_t = (-eta / (1 + k)) ** (1 + k)
    pulled = [geo.g[0, 0] + a ** 2, geo.g[1, 1] - a ** 2, sp.diff(-minus_t, eta) - a,
              sp.powdenest(minus_t ** (2 * k / (1 + k)), force=True) - a ** 2]
    if any(vm.norm(value) != 0 for value in pulled):
        raise AssertionError(f"{METRIC}/conformal: the chart is not the synchronous one in conformal time")
    phi = (2 * k - 1) * sp.log(-eta / (1 + k))
    gamma = geo.christoffel_ull()
    gradient = [sp.diff(phi, s) for s in chart.symbols]
    hessian = [[sp.diff(phi, chart.symbols[i], chart.symbols[j])
                - sum(gamma[c][i][j] * gradient[c] for c in range(n)) for j in range(n)] for i in range(n)]
    ricci = geo.ricci_ll()
    box = sum(geo.ginv[i, j] * hessian[i][j] for i in range(n) for j in range(n))
    square = sum(geo.ginv[i, j] * gradient[i] * gradient[j] for i in range(n) for j in range(n))
    residues = [ricci[i][j] + hessian[i][j] for i in range(n) for j in range(n)]
    residues.append(geo.ricci_scalar() + 2 * box - square)
    for residue in residues:
        value = sp.factor(vm.norm(residue))
        if value == 0:
            continue
        top = sp.fraction(value)[0]
        if not any(sp.expand(f.as_base_exp()[0] - quadratic) == 0 or sp.expand(f.as_base_exp()[0] + quadratic) == 0
                   for f in sp.Mul.make_args(top)):
            raise AssertionError(f"{METRIC}/conformal: a field equation fails at the members: {value}")
