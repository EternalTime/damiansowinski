"""The charged dilaton black holes of Chan, Horne and Mann (Nucl. Phys. B 447, 441, 1995), whose
region far away is neither flat nor de Sitter nor anti-de Sitter space.

The action is their (1.1), S = int d^4x sqrt(-g) (R - 2(grad phi)^2 - V(phi) - e^{-2 a phi} F^2), in
units G = c = 1 and signature (-,+,+,+), and the field equations are their (2.1):

    R_ab = 2 d_a phi d_b phi + (1/2) g_ab V + 2 e^{-2 a phi} (F_ac F_b^c - (1/4) g_ab F^2),
    d_a (sqrt(-g) e^{-2 a phi} F^ab) = 0,
    box phi = (1/4) dV/dphi - (a/2) e^{-2 a phi} F^2.

The charts and their sources:

chan_horne_mann: their (2.2) and (3.7), ds^2 = -U dt^2 + dr^2/U + R^2 dOmega^2 with
    U = r^{2/(1+a^2)} (1 - r_h/r)/gamma^2 and R = gamma r^{a^2/(1+a^2)}, r_h = 2(1+a^2)M/a^2 their
    (3.8), written with Clement and Leygnac's length r_0 of their (2.8) (Phys. Rev. D 70, 084018,
    2004), gamma^2 = r_0^{2/(1+a^2)}. The constant gamma carries the dimension of a length to the
    power 1/(1+a^2), and the checker reads a power whose exponent holds a parameter with its base in
    a fixed unit, so the chart is written in the pure number h = (r/r_0)^{2/(1+a^2)}, with
    U = h(1 - r_h/r) and R^2 = r^2/h. h is held as a function of r with the slope 2h/((1+a^2)r),
    vm.HELD and vm.RATES, so every value is a rational function of r, h, a and r_h; printed with
    the power written out, sympy's factoring of exponents in 1/(1+a^2) did not finish in ten minutes.

The electric field is F_tr = Q e^{2 a phi}/R^2, their (2.3), with e^{2 a phi} = R^2/((1+a^2)Q^2) from
their (3.7c), so F_tr = 1/((1+a^2)Q); the magnetic hole has F_theta phi = Q sin(theta) and the
dilaton reversed, phi -> -phi, as their section 2 says. Every check is made at points, in forty
digits, independent of the printed tensors: the metric against the source's own form, the three
field equations for the electric and the magnetic hole, and the Ricci scalar against their (3.9).
dilaton_unusual_asymptotics.md beside this file records each chart's source.
"""
import sympy as sp

METRIC = "dilaton_unusual_asymptotics"
SYSTEMS = ("chan_horne_mann",)
NAMES = {"chan_horne_mann": "Chan-Horne-Mann"}
COORDS = {"chan_horne_mann": ["t", "r", "\\theta", "\\phi"]}

H = "h = \\left(\\dfrac{r}{r_0}\\right)^{2/(1 + a^2)}"
SLOPE = "\\dfrac{2h}{\\left(1 + a^2\\right)r}"
PARAMETERS = {"chan_horne_mann": ["a", "r_0", "r_h", H]}
HELD = {"chan_horne_mann": ("h",)}
SPHERE = "\\left(d\\theta^2 + \\sin^2\\theta\\,d\\phi^2\\right)"
LAPSE = "\\left(1 - \\dfrac{r_h}{r}\\right)"
LINES = {
    "chan_horne_mann": ("ds^2 = -h" + LAPSE + "dt^2 + \\dfrac{dr^2}{h" + LAPSE + "} + \\dfrac{r^2}{h}" + SPHERE),
}
LINE = "(-\\infty, \\infty)"
ANGLES = ["\\theta \\in [0, \\pi]", "\\phi \\in [0, 2\\pi)"]
DOMAINS = {
    "chan_horne_mann": ["t \\in " + LINE, "r \\in (0, \\infty)"] + ANGLES + [
        "a \\neq 0,\\quad r_0 > 0,\\quad r_h > 0",
        "r = r_h \\;\\text{(the event horizon)}", "r = 0 \\;\\text{(the singularity)}"],
}
TOLERANCE = 1e-30


def register_dimensions(vm):
    """With G = c = 1 every coordinate but the angles is a length, as are r_0 and r_h, and the
    coupling a is a pure number."""
    dims = {"t": "L", "r": "L", "\\theta": "1", "\\phi": "1", "a": "1", "r_0": "L", "r_h": "L", "h": "1"}
    for system in SYSTEMS:
        names = [p.split(" = ")[0] for p in PARAMETERS[system]]
        vm.DIMENSIONS[(METRIC, system)] = {c: dims[c] for c in COORDS[system] + names}
        if system in HELD:
            vm.HELD[(METRIC, system)] = HELD[system]
            vm.RATES[(METRIC, system)] = {"h": {"r": SLOPE}}


def spec(system):
    import verify_metrics as vm
    probe = vm.Reader(COORDS[system], PARAMETERS[system], (), held=HELD.get(system, ()))
    P, X = probe.parameters, probe.symbol
    a, r, h = P["a"], X["r"], P.get("h")

    def reduce(value):
        # The slope of h is 2h/((1 + a^2) r), so no derivative of h is left.
        value = sp.sympify(value)
        while value.has(sp.Derivative):
            for d in sorted(value.atoms(sp.Derivative), key=lambda d: -d.derivative_count):
                if d.expr != h:
                    raise AssertionError(f"{METRIC}: {d} is no derivative of h")
                value = value.xreplace({d: sp.diff(2 * h / ((1 + a ** 2) * r), r, d.derivative_count - 1)})
            value = value.doit()
        return vm.norm(value)

    out = {
        "metric_id": METRIC,
        "system": {"id": system, "name": NAMES[system], "coords": COORDS[system], "domains": DOMAINS[system],
                   "parameters": PARAMETERS[system], "line_element": LINES[system]},
        "chart_line_element": LINES[system],
        "printer": {"lead": [a, r, P["r_h"], h] if h is not None else [a, r, P["r_h"]], "rising": [a]},
        "check": lambda chart: check(chart, system),
    }
    if system == "chan_horne_mann":
        out["reduce"] = reduce
        out["components"] = {
            "metric_components": {("t", "t"): "-h" + LAPSE, ("r", "r"): "\\dfrac{1}{h" + LAPSE + "}"},
            "inverse_metric_components": {("t", "t"): "-\\dfrac{1}{h" + LAPSE + "}", ("r", "r"): "h" + LAPSE}}
        out["ricci_scalar"] = "\\dfrac{2a^2h\\left(r - r_h\\right)}{\\left(1 + a^2\\right)^2r^3}"
    return out


def register(module):
    module.CHARTS[METRIC] = [lambda system=system: spec(system) for system in SYSTEMS]


# -- Checks at points, independent of the printed tensors ---------------------------------------

def small(values, where, what):
    worst = max([abs(complex(sp.N(x, 40))) for x in values] + [0])
    if worst > TOLERANCE:
        raise AssertionError(f"{METRIC}/{where}: {what} ({worst:.3g})")


def field_equations(g, ricci, symbols, at, phi, F, a, V=0, dV=0):
    """The residuals of the three field equations at the point `at`, with the metric g, its Ricci
    tensor, the dilaton phi and the Maxwell field F held as expressions in the chart's symbols."""
    gi = g.inv()
    gv, giv = g.subs(at), gi.subs(at)
    dphi = [sp.diff(phi, x) for x in symbols]
    Fup = gi * F * gi
    F2 = sum(F[i, j] * Fup[i, j] for i in range(4) for j in range(4))
    out = []
    for i in range(4):
        for j in range(i, 4):
            FF = sum(F[i, k] * F[j, l] * giv[k, l] for k in range(4) for l in range(4))
            rhs = (2 * dphi[i] * dphi[j] + gv[i, j] * V / 2
                   + 2 * sp.exp(-2 * a * phi) * (FF - gv[i, j] * F2 / 4))
            out.append((ricci[i][j] - rhs).subs(at))
    root = sp.sqrt(-g.det())
    box = sum(sp.diff(root * gi[i, j] * dphi[j], symbols[i]) for i in range(4) for j in range(4)) / root
    out.append((box - dV / 4 + a * sp.exp(-2 * a * phi) * F2 / 2).subs(at))
    for b in range(4):
        out.append(sum(sp.diff(root * sp.exp(-2 * a * phi) * Fup[i, b], symbols[i]) for i in range(4)).subs(at))
    return out


POINTS = [
    # (a, r_0, r_h, Q), then (r, theta)
    ((sp.Rational(1, 2), sp.Rational(6, 5), sp.Rational(3, 2), sp.Rational(4, 7)), (sp.Rational(7, 3), sp.Rational(9, 10))),
    ((sp.Integer(1), sp.Rational(5, 4), sp.Rational(2, 3), sp.Rational(3, 5)), (sp.Rational(19, 4), sp.Rational(1, 3))),
    ((sp.Rational(7, 5), sp.Rational(2, 3), sp.Rational(5, 2), sp.Rational(9, 7)), (sp.Rational(11, 3), sp.Rational(6, 5))),
    ((sp.Integer(3), sp.Rational(3, 2), sp.Rational(1, 4), sp.Rational(2, 9)), (sp.Rational(5, 7), sp.Rational(2, 7))),
]


def check(chart, system):
    s = chart.symbols
    t, r, th, ph = s
    P = chart.reader.parameters
    a, r0, rh = P["a"], P["r_0"], P["r_h"]
    # The held name written out, so that the field equations below differentiate the metric itself.
    written = {P["h"]: (r / r0) ** (2 / (1 + a ** 2))}
    g = chart.geo.g.subs(written)
    ricci = [[chart.geo.ricci_ll()[i][j].subs(written) for j in range(4)] for i in range(4)]
    scalar = chart.geo.ricci_scalar().subs(written)
    Q = sp.Symbol("Q", positive=True)
    # Chan, Horne and Mann's (3.7) with gamma^2 = r_0^{2/(1+a^2)}, and M from r_h by their (3.8).
    gamma = r0 ** (1 / (1 + a ** 2))
    U = r ** (2 / (1 + a ** 2)) / gamma ** 2 * (1 - rh / r)
    R = gamma * r ** (a ** 2 / (1 + a ** 2))
    source = sp.diag(-U, 1 / U, R ** 2, R ** 2 * sp.sin(th) ** 2)
    phi = sp.log(R ** 2 / ((1 + a ** 2) * Q ** 2)) / (2 * a)
    electric = sp.zeros(4, 4)
    electric[0, 1], electric[1, 0] = Q * sp.exp(2 * a * phi) / R ** 2, -Q * sp.exp(2 * a * phi) / R ** 2
    magnetic = sp.zeros(4, 4)
    magnetic[2, 3], magnetic[3, 2] = Q * sp.sin(th), -Q * sp.sin(th)
    M = a ** 2 * rh / (2 * (1 + a ** 2))
    scalar_39 = 2 * a ** 2 * r ** (-2 * a ** 2 / (1 + a ** 2)) / ((1 + a ** 2) ** 2 * gamma ** 2) * (
        1 - 2 * (1 + a ** 2) * M / (a ** 2 * r))
    for (av, r0v, rhv, Qv), (rv, thv) in POINTS:
        values = {a: av, r0: r0v, rh: rhv, Q: Qv}
        at = {**values, r: rv, th: thv, t: sp.Rational(1, 3), ph: sp.Rational(2, 5)}
        where = f"{system} at (r, theta) = {(rv, thv)}, (a, r_0, r_h, Q) = {(av, r0v, rhv, Qv)}"
        small(list((g - source).subs(at)), where, "the metric is not Chan, Horne and Mann's (3.7)")
        small([(Q * sp.exp(2 * a * phi) / R ** 2 - 1 / ((1 + a ** 2) * Q)).subs(at)], where,
              "F_tr is not 1/((1 + a^2) Q)")
        # The coupling and the constants are set before anything is differentiated.
        gv = g.subs(values)
        rv_ = [[ricci[i][j].subs(values) for j in range(4)] for i in range(4)]
        small(field_equations(gv, rv_, s, at, phi.subs(values), electric.subs(values), av), where,
              "the electric hole does not solve the field equations")
        small(field_equations(gv, rv_, s, at, -phi.subs(values), magnetic.subs(values), av), where,
              "the magnetic hole does not solve the field equations")
        small([(scalar - scalar_39).subs(at)], where, "the Ricci scalar is not their (3.9)")
