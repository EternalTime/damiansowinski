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
areal: their (3.2) at a = 1, the areal radius sqrt(r_0 r) of the first chart taken as the
    coordinate, with gamma^2 = r_0 and 4M = r_h.
string_electric: their (3.5), the string metric e^{2 phi} g of the electric hole at a = 1, in
    their rescaled radius. They print the rescaling as r -> sqrt(2) Q r/gamma; the metric they
    print follows from r -> sqrt(2) Q r/gamma^2, which the check below uses, and only that one
    gives their (3.5) from their (3.1).
string_magnetic: their (3.6), the string metric e^{-2 phi} g of the magnetic hole at a = 1, in the
    radius of the first chart.
massless: Clement, Gal'tsov and Leygnac's (2.10) and (2.11) (Phys. Rev. D 67, 024012, 2003), the
    hole of no mass at a = 1 with their signature reversed. Their r_0 = sqrt(2) Qbar is this r_0
    and their alpha^{-2} = 4 sqrt(2) Qbar is 4 r_0, so xi = 2 sqrt(r_0 r) cosh(t/2r_0) and
    eta = 2 sqrt(r_0 r) sinh(t/2r_0) for the t and r of the first chart at r_h = 0.

The electric field is F_tr = Q e^{2 a phi}/R^2 in the first chart, their (2.3), with
e^{2 a phi} = R^2/((1+a^2)Q^2) from their (3.7c), so F_tr = 1/((1+a^2)Q); the magnetic hole has
F_theta phi = Q sin(theta) and the dilaton reversed, phi -> -phi, as their section 2 says. Every
check is made at points, in forty digits, independent of the printed tensors: each chart against
its source's own form, each chart after the first against the first carried along its map, the
three field equations for the electric and the magnetic hole in each chart of the Einstein metric,
each string metric as e^{+-2 phi} times the Einstein metric, and the Ricci scalar of the first
chart against their (3.9). dilaton_unusual_asymptotics.md beside this file records each chart's source.
"""
import sympy as sp

METRIC = "dilaton_unusual_asymptotics"
SYSTEMS = ("chan_horne_mann", "areal", "string_electric", "string_magnetic", "massless")
NAMES = {"chan_horne_mann": "Chan-Horne-Mann", "areal": "Areal Radius, $a = 1$",
         "string_electric": "String Metric, Electric Charge", "string_magnetic": "String Metric, Magnetic Charge",
         "massless": "No Mass, $a = 1$"}
COORDS = {system: ["t", "r", "\\theta", "\\phi"] for system in SYSTEMS}
COORDS["massless"] = ["\\eta", "\\xi", "\\theta", "\\phi"]

H = "h = \\left(\\dfrac{r}{r_0}\\right)^{2/(1 + a^2)}"
SLOPE = "\\dfrac{2h}{\\left(1 + a^2\\right)r}"
PARAMETERS = {"chan_horne_mann": ["a", "r_0", "r_h", H], "areal": ["r_0", "r_h"],
              "string_electric": ["r_0", "r_h", "Q"], "string_magnetic": ["r_0", "r_h", "Q"], "massless": []}
HELD = {"chan_horne_mann": ("h",)}
SPHERE = "\\left(d\\theta^2 + \\sin^2\\theta\\,d\\phi^2\\right)"
LAPSE = "\\left(1 - \\dfrac{r_h}{r}\\right)"
STRING_LAPSE = "\\left(1 - \\dfrac{r_0\\,r_h}{\\sqrt{2}\\,Qr}\\right)"
LINES = {
    "chan_horne_mann": "ds^2 = -h" + LAPSE + "dt^2 + \\dfrac{dr^2}{h" + LAPSE + "} + \\dfrac{r^2}{h}" + SPHERE,
    "areal": ("ds^2 = -\\dfrac{r^2 - r_0\\,r_h}{r_0^2}dt^2 + \\dfrac{4r^2\\,dr^2}{r^2 - r_0\\,r_h} + r^2" + SPHERE),
    "string_electric": ("ds^2 = -\\dfrac{r^2}{r_0^2}" + STRING_LAPSE + "dt^2 + \\dfrac{dr^2}{1 - \\dfrac{r_0\\,r_h}{\\sqrt{2}\\,Qr}}"
                        " + r^2" + SPHERE),
    "string_magnetic": ("ds^2 = -\\dfrac{2Q^2}{r_0^2}" + LAPSE + "dt^2 + \\dfrac{2Q^2\\,dr^2}{r^2" + LAPSE + "} + 2Q^2"
                        + SPHERE),
    "massless": "ds^2 = -d\\eta^2 + d\\xi^2 + \\dfrac{\\xi^2 - \\eta^2}{4}" + SPHERE,
}
LINE = "(-\\infty, \\infty)"
ANGLES = ["\\theta \\in [0, \\pi]", "\\phi \\in [0, 2\\pi)"]
SINGULARITY = "r = 0 \\;\\text{(the singularity)}"
DOMAINS = {
    "chan_horne_mann": ["t \\in " + LINE, "r \\in (0, \\infty)"] + ANGLES + [
        "a \\neq 0,\\quad r_0 > 0,\\quad r_h > 0", "r = r_h \\;\\text{(the event horizon)}", SINGULARITY],
    "areal": ["t \\in " + LINE, "r \\in (0, \\infty)"] + ANGLES + [
        "r_0 > 0,\\quad r_h > 0", "r = \\sqrt{r_0\\,r_h} \\;\\text{(the event horizon)}", SINGULARITY],
    "string_electric": ["t \\in " + LINE, "r \\in (0, \\infty)"] + ANGLES + [
        "r_0 > 0,\\quad r_h > 0,\\quad Q > 0", "r = r_0\\,r_h/(\\sqrt{2}\\,Q) \\;\\text{(the event horizon)}",
        SINGULARITY],
    "string_magnetic": ["t \\in " + LINE, "r \\in (0, \\infty)"] + ANGLES + [
        "r_0 > 0,\\quad r_h > 0,\\quad Q > 0", "r = r_h \\;\\text{(the event horizon)}", SINGULARITY],
    "massless": ["\\eta \\in " + LINE, "\\xi \\in (|\\eta|, \\infty)"] + ANGLES + [
        "\\xi = |\\eta| \\;\\text{(the null singularity)}"],
}
TOLERANCE = 1e-30


def register_dimensions(vm):
    """With G = c = 1 every coordinate but the angles is a length, as are r_0, r_h and the charge Q,
    and the coupling a and the name h are pure numbers."""
    dims = {"t": "L", "r": "L", "\\eta": "L", "\\xi": "L", "\\theta": "1", "\\phi": "1", "a": "1", "r_0": "L",
            "r_h": "L", "Q": "L", "h": "1"}
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
    out = {
        "metric_id": METRIC,
        "system": {"id": system, "name": NAMES[system], "coords": COORDS[system], "domains": DOMAINS[system],
                   "parameters": PARAMETERS[system], "line_element": LINES[system]},
        "chart_line_element": LINES[system],
        "printer": {"lead": [X[COORDS[system][1]]] + [P[name] for name in ("r_0", "r_h", "Q") if name in P]},
        "check": lambda chart: check(chart, system),
    }
    if system == "chan_horne_mann":
        a, r, h = P["a"], X["r"], P["h"]

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

        out["printer"] = {"lead": [a, r, P["r_h"], h], "rising": [a]}
        out["reduce"] = reduce
        out["components"] = {
            "metric_components": {("t", "t"): "-h" + LAPSE, ("r", "r"): "\\dfrac{1}{h" + LAPSE + "}"},
            "inverse_metric_components": {("t", "t"): "-\\dfrac{1}{h" + LAPSE + "}", ("r", "r"): "h" + LAPSE}}
        out["ricci_scalar"] = "\\dfrac{2a^2h\\left(r - r_h\\right)}{\\left(1 + a^2\\right)^2r^3}"
    elif system == "areal":
        out["components"] = {
            "metric_components": {("t", "t"): "-\\dfrac{r^2 - r_0\\,r_h}{r_0^2}"},
            "inverse_metric_components": {("t", "t"): "-\\dfrac{r_0^2}{r^2 - r_0\\,r_h}"}}
    elif system == "massless":
        # The radius of the spheres, (xi^2 - eta^2)/4, is kept whole.
        whole = "\\dfrac{\\xi^2 - \\eta^2}{4}"
        out["components"] = {
            "metric_components": {("\\theta", "\\theta"): whole, ("\\phi", "\\phi"): whole + "\\sin^2\\theta"},
            "inverse_metric_components": {("\\theta", "\\theta"): "\\dfrac{4}{\\xi^2 - \\eta^2}",
                                          ("\\phi", "\\phi"): "\\dfrac{4}{\\left(\\xi^2 - \\eta^2\\right)\\sin^2\\theta}"}}
        out["ricci_scalar"] = "\\dfrac{2}{\\xi^2 - \\eta^2}"
        out["kretschmann"] = "\\dfrac{44}{\\left(\\xi^2 - \\eta^2\\right)^2}"
    elif system == "string_magnetic":
        out["components"] = {
            "metric_components": {("t", "t"): "-\\dfrac{2Q^2}{r_0^2}" + LAPSE,
                                  ("r", "r"): "\\dfrac{2Q^2}{r^2" + LAPSE + "}"},
            "inverse_metric_components": {("t", "t"): "-\\dfrac{r_0^2}{2Q^2" + LAPSE + "}",
                                          ("r", "r"): "\\dfrac{r^2}{2Q^2}" + LAPSE}}
    elif system == "string_electric":
        out["components"] = {
            "metric_components": {("t", "t"): "-\\dfrac{r^2}{r_0^2}" + STRING_LAPSE,
                                  ("r", "r"): "\\dfrac{1}{1 - \\dfrac{r_0\\,r_h}{\\sqrt{2}\\,Qr}}"},
            "inverse_metric_components": {("t", "t"): "-\\dfrac{r_0^2}{r^2" + STRING_LAPSE + "}",
                                          ("r", "r"): "1 - \\dfrac{r_0\\,r_h}{\\sqrt{2}\\,Qr}"}}
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


def two_form(value, i, j):
    F = sp.zeros(4, 4)
    F[i, j], F[j, i] = value, -value
    return F


def first_chart(a, r0, rh, t, r, th):
    """Chan, Horne and Mann's (3.7) in the first chart, with gamma^2 = r_0^{2/(1+a^2)}: the metric,
    the areal radius R, and the dilaton of the electric hole of charge Q."""
    gamma = r0 ** (1 / (1 + a ** 2))
    U = r ** (2 / (1 + a ** 2)) / gamma ** 2 * (1 - rh / r)
    R = gamma * r ** (a ** 2 / (1 + a ** 2))
    return sp.diag(-U, 1 / U, R ** 2, R ** 2 * sp.sin(th) ** 2), R


# (a, r_0, r_h, Q) and (first coordinate, second coordinate, theta) at each point.
POINTS = {
    "chan_horne_mann": [
        ((sp.Rational(1, 2), sp.Rational(6, 5), sp.Rational(3, 2), sp.Rational(4, 7)), (sp.Rational(1, 3), sp.Rational(7, 3), sp.Rational(9, 10))),
        ((sp.Integer(1), sp.Rational(5, 4), sp.Rational(2, 3), sp.Rational(3, 5)), (sp.Rational(1, 3), sp.Rational(19, 4), sp.Rational(1, 3))),
        ((sp.Rational(7, 5), sp.Rational(2, 3), sp.Rational(5, 2), sp.Rational(9, 7)), (sp.Rational(1, 3), sp.Rational(11, 3), sp.Rational(6, 5))),
        ((sp.Integer(3), sp.Rational(3, 2), sp.Rational(1, 4), sp.Rational(2, 9)), (sp.Rational(1, 3), sp.Rational(5, 7), sp.Rational(2, 7))),
    ],
    "a1": [
        ((sp.Integer(1), sp.Rational(5, 4), sp.Rational(2, 3), sp.Rational(3, 5)), (sp.Rational(2, 7), sp.Rational(19, 4), sp.Rational(1, 3))),
        ((sp.Integer(1), sp.Rational(2, 3), sp.Rational(7, 2), sp.Rational(9, 7)), (sp.Rational(-3, 5), sp.Rational(1, 2), sp.Rational(6, 5))),
        ((sp.Integer(1), sp.Rational(9, 4), sp.Rational(1, 5), sp.Rational(2, 9)), (sp.Rational(5, 3), sp.Rational(13, 2), sp.Rational(2, 7))),
    ],
}


def check(chart, system):
    s = chart.symbols
    x0, x1, th, ph = s
    P = chart.reader.parameters
    g = chart.geo.g
    if system == "chan_horne_mann":
        # The held name written out, so that the field equations differentiate the metric itself.
        a, r0, rh = P["a"], P["r_0"], P["r_h"]
        g = g.subs({P["h"]: (x1 / r0) ** (2 / (1 + a ** 2))})
    else:
        a, r0, rh = sp.Integer(1), P.get("r_0", sp.Symbol("r_0", positive=True)), P.get("r_h", sp.Symbol("r_h", positive=True))
    Q = P.get("Q", sp.Symbol("Q", positive=True))
    held = {P["h"]: (x1 / r0) ** (2 / (1 + a ** 2))} if "h" in P else {}
    ricci = [[chart.geo.ricci_ll()[i][j].subs(held) for j in range(4)] for i in range(4)]
    scalar = chart.geo.ricci_scalar().subs(held)
    # The first chart's coordinates as functions of this chart's, with the Jacobian of (t, r).
    t1, r1 = sp.symbols("t1 r1", positive=True)
    if system == "areal":
        image = (x0, x1 ** 2 / r0)
    elif system == "string_electric":
        image = (x0, sp.sqrt(2) * Q * x1 / r0)
    elif system == "massless":
        image = (2 * r0 * sp.atanh(x0 / x1), (x1 ** 2 - x0 ** 2) / (4 * r0))
    else:
        image = (x0, x1)
    J = sp.Matrix(4, 4, lambda i, j: (sp.diff(image[i], s[j]) if i < 2 else (1 if i == j else 0)))
    first, R1 = first_chart(a, r0, rh, t1, r1, th)
    pulled = J.T * first.subs({t1: image[0], r1: image[1]}) * J
    R = R1.subs(r1, image[1])
    # The electric hole: e^{2 a phi} = R^2/((1 + a^2) Q^2), F = F_tr dt ^ dr with F_tr = 1/((1 + a^2) Q),
    # carried along the map; the magnetic hole: F = Q sin(theta) dtheta ^ dphi with the dilaton reversed.
    phi = sp.log(R ** 2 / ((1 + a ** 2) * Q ** 2)) / (2 * a)
    electric = (J[:2, :2].det() / ((1 + a ** 2) * Q)) * two_form(1, 0, 1)
    magnetic = two_form(Q * sp.sin(th), 2, 3)
    if system == "chan_horne_mann":
        M = a ** 2 * rh / (2 * (1 + a ** 2))
        gamma = r0 ** (1 / (1 + a ** 2))
        source = first
        source = source.subs({t1: x0, r1: x1})
        scalar_39 = 2 * a ** 2 * x1 ** (-2 * a ** 2 / (1 + a ** 2)) / ((1 + a ** 2) ** 2 * gamma ** 2) * (
            1 - 2 * (1 + a ** 2) * M / (a ** 2 * x1))
    elif system == "areal":
        # Their (3.2) with gamma^2 = r_0 and M = r_h/4.
        gamma2, M = r0, rh / 4
        source = sp.diag(-(x1 ** 2 - 4 * gamma2 * M) / gamma2 ** 2, 4 * x1 ** 2 / (x1 ** 2 - 4 * gamma2 * M),
                         x1 ** 2, x1 ** 2 * sp.sin(th) ** 2)
    elif system == "string_electric":
        # Their (3.5) with gamma^2 = r_0 and M = r_h/4.
        gamma2, M = r0, rh / 4
        f = 1 - 2 * sp.sqrt(2) * gamma2 * M / (Q * x1)
        source = sp.diag(-x1 ** 2 / gamma2 ** 2 * f, 1 / f, x1 ** 2, x1 ** 2 * sp.sin(th) ** 2)
    elif system == "string_magnetic":
        # Their (3.6) with gamma^2 = r_0 and M = r_h/4.
        gamma2, M = r0, rh / 4
        f = 1 - 4 * M / x1
        source = sp.diag(-2 * Q ** 2 / gamma2 ** 2 * f, 2 * Q ** 2 / (x1 ** 2 * f), 2 * Q ** 2,
                         2 * Q ** 2 * sp.sin(th) ** 2)
    else:
        # Clement, Gal'tsov and Leygnac's (2.11), d eta^2 - d xi^2 - (xi^2 - eta^2) dOmega^2/4, reversed.
        source = sp.diag(-1, 1, (x1 ** 2 - x0 ** 2) / 4, (x1 ** 2 - x0 ** 2) / 4 * sp.sin(th) ** 2)
    if system == "string_electric":
        einstein, sign = g * sp.exp(-2 * phi), 1
    elif system == "string_magnetic":
        einstein, sign = g * sp.exp(2 * phi), -1
    else:
        einstein, sign = g, 0
    points = POINTS["chan_horne_mann" if system == "chan_horne_mann" else "a1"]
    for (av, r0v, rhv, Qv), (u, v, thv) in points:
        values = {a: av, r0: r0v, rh: (0 if system == "massless" else rhv), Q: Qv} if system == "chan_horne_mann" \
            else {r0: r0v, rh: (0 if system == "massless" else rhv), Q: Qv}
        if system == "massless":
            u, v = v * sp.Rational(1, 3), v
        at = {**values, x0: u, x1: v, th: thv, ph: sp.Rational(2, 5)}
        where = f"{system} at {(u, v, thv)}, (a, r_0, r_h, Q) = {(av, r0v, values[rh], Qv)}"
        small(list((g - source).subs(at)), where, "the metric is not its source's")
        small(list((einstein - pulled).subs(at)), where,
              "the Einstein metric is not the first chart carried along its map" if sign == 0 else
              "the string metric is not e^(+-2 phi) times the Einstein metric")
        if sign == 0:
            gv = g.subs(values)
            rc = [[ricci[i][j].subs(values) for j in range(4)] for i in range(4)]
            av_ = values.get(a, a)
            small(field_equations(gv, rc, s, at, phi.subs(values), electric.subs(values), av_), where,
                  "the electric hole does not solve the field equations")
            small(field_equations(gv, rc, s, at, -phi.subs(values), magnetic.subs(values), av_), where,
                  "the magnetic hole does not solve the field equations")
        if system == "chan_horne_mann":
            small([(scalar - scalar_39).subs(at)], where, "the Ricci scalar is not their (3.9)")
            small([(Q * sp.exp(2 * a * phi) / R ** 2 - 1 / ((1 + a ** 2) * Q)).subs(at)], where,
                  "F_tr is not 1/((1 + a^2) Q)")
