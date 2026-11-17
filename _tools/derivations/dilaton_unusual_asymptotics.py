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
liouville_inverse, liouville_minus_a and liouville_minus_inverse: their three families with the
    potential V = 2 Lambda e^{2 b phi}, their (4.1), at b = 1/a, -a and -1/a, their (4.2), (4.6) and
    (4.9). Each is written as the first chart is, with r_0 for gamma and a pure number for the power
    of r/r_0: (4.2) has N = 1/(1+a^2), so gamma^2 = r_0^{2a^2/(1+a^2)}, k = (r/r_0)^{2a^2/(1+a^2)},
    U_1 = (1+a^2) k (1/(1-a^2) - 2M/r + q^2/r^2) and R^2 = r^2/k with q^2 = Q^2 e^{2a phi_0}/gamma^2;
    (4.6) and (4.9) have the first chart's N and h, U_2 = sigma h (1 - r_h/r) with
    sigma = (1+a^2)(2Q^2 e^{2a phi_0}/gamma^2 - 1)/(1-a^2) and r_h their (4.7), and
    U_3 = h (1 - r_h/r) + lambda r^2/h with lambda = (1+a^2)^{(2a^2+1)/a^2} (Q/r_0)^{2/a^2} Lambda
    /((1-3a^2) a^2), their r^{-2(1-a^2)/(1+a^2)} term written through h. Each is checked against
    their own form and against the three field equations with the potential, at points where the
    dilaton's constant phi_0 and the charge Q are chosen and Lambda follows by their (4.2d), (4.6d)
    or is read from lambda.
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
SYSTEMS = ("chan_horne_mann", "areal", "string_electric", "string_magnetic", "liouville_inverse", "liouville_minus_a",
           "liouville_minus_inverse", "massless")
LIOUVILLE = ("liouville_inverse", "liouville_minus_a", "liouville_minus_inverse")
NAMES = {"chan_horne_mann": "Chan-Horne-Mann", "areal": "Areal Radius, $a = 1$",
         "string_electric": "String Metric, Electric Charge", "string_magnetic": "String Metric, Magnetic Charge",
         "liouville_inverse": "Liouville Potential, $b = 1/a$", "liouville_minus_a": "Liouville Potential, $b = -a$",
         "liouville_minus_inverse": "Liouville Potential, $b = -1/a$", "massless": "No Mass, $a = 1$"}
COORDS = {system: ["t", "r", "\\theta", "\\phi"] for system in SYSTEMS}
COORDS["massless"] = ["\\eta", "\\xi", "\\theta", "\\phi"]

H = "h = \\left(\\dfrac{r}{r_0}\\right)^{2/(1 + a^2)}"
SLOPE = "\\dfrac{2h}{\\left(1 + a^2\\right)r}"
K = "k = \\left(\\dfrac{r}{r_0}\\right)^{2a^2/(1 + a^2)}"
K_SLOPE = "\\dfrac{2a^2k}{\\left(1 + a^2\\right)r}"
PARAMETERS = {"chan_horne_mann": ["a", "r_0", "r_h", H], "areal": ["r_0", "r_h"],
              "string_electric": ["r_0", "r_h", "Q"], "string_magnetic": ["r_0", "r_h", "Q"],
              "liouville_inverse": ["a", "r_0", "M", "q", K], "liouville_minus_a": ["a", "r_0", "r_h", "\\sigma", H],
              "liouville_minus_inverse": ["a", "r_0", "r_h", "\\lambda", H], "massless": []}
HELD = {"chan_horne_mann": ("h",), "liouville_inverse": ("k",), "liouville_minus_a": ("h",),
        "liouville_minus_inverse": ("h",)}
RATES = {"chan_horne_mann": {"h": {"r": SLOPE}}, "liouville_inverse": {"k": {"r": K_SLOPE}},
         "liouville_minus_a": {"h": {"r": SLOPE}}, "liouville_minus_inverse": {"h": {"r": SLOPE}}}
# The Ricci scalar of each Liouville family, grouped by the powers of r; print_charts checks each.
LIOUVILLE_SCALAR = {
    "liouville_inverse": ("-\\dfrac{2a^2k\\left(\\left(1 + 2a^2\\right)r^2 + 2\\left(1 - a^2\\right)Mr - \\left(1 - a^2\\right)q^2\\right)}"
                          "{\\left(1 - a^2\\right)\\left(1 + a^2\\right)r^4}"),
    "liouville_minus_a": ("\\dfrac{2h\\left(\\left(1 + a^2\\right)^2r - \\left(1 + a^2 + a^4\\right)\\sigma r - a^2\\sigma r_h\\right)}"
                          "{\\left(1 + a^2\\right)^2r^3}"),
    "liouville_minus_inverse": ("\\dfrac{2a^2\\left(3\\left(1 - 2a^2\\right)\\lambda r^3 + h^2\\left(r - r_h\\right)\\right)}"
                                "{\\left(1 + a^2\\right)^2r^3h}"),
}
# The factor of each Liouville family's g_tt, Chan, Horne and Mann's U.
LIOUVILLE_U = {
    "liouville_inverse": "\\left(1 + a^2\\right)k\\left(\\dfrac{1}{1 - a^2} - \\dfrac{2M}{r} + \\dfrac{q^2}{r^2}\\right)",
    "liouville_minus_a": "\\sigma h\\left(1 - \\dfrac{r_h}{r}\\right)",
    "liouville_minus_inverse": "\\left(h\\left(1 - \\dfrac{r_h}{r}\\right) + \\dfrac{\\lambda r^2}{h}\\right)",
}
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
    **{system: ("ds^2 = -" + LIOUVILLE_U[system] + "dt^2 + \\dfrac{dr^2}{" + LIOUVILLE_U[system] + "} + \\dfrac{r^2}{"
                + ("k" if system == "liouville_inverse" else "h") + "}" + SPHERE) for system in LIOUVILLE},
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
    "liouville_inverse": ["t \\in " + LINE, "r \\in (0, \\infty)"] + ANGLES + [
        "a \\neq 0,\\quad a^2 \\neq 1,\\quad r_0 > 0",
        "\\dfrac{1}{1 - a^2} - \\dfrac{2M}{r} + \\dfrac{q^2}{r^2} = 0 \\;\\text{(the horizons)}", SINGULARITY],
    "liouville_minus_a": ["t \\in " + LINE, "r \\in (0, \\infty)"] + ANGLES + [
        "a \\neq 0,\\quad a^2 \\neq 1,\\quad r_0 > 0,\\quad r_h > 0,\\quad \\sigma > 0",
        "r = r_h \\;\\text{(the event horizon)}", SINGULARITY],
    "liouville_minus_inverse": ["t \\in " + LINE, "r \\in (0, \\infty)"] + ANGLES + [
        "a \\neq 0,\\quad a^2 \\neq 1/3,\\quad r_0 > 0",
        "h\\left(1 - \\dfrac{r_h}{r}\\right) + \\dfrac{\\lambda r^2}{h} = 0 \\;\\text{(the horizons)}", SINGULARITY],
    "massless": ["\\eta \\in " + LINE, "\\xi \\in (|\\eta|, \\infty)"] + ANGLES + [
        "\\xi = |\\eta| \\;\\text{(the null singularity)}"],
}
TOLERANCE = 1e-30


def register_dimensions(vm):
    """With G = c = 1 every coordinate but the angles is a length, as are r_0, r_h and the charge Q,
    and the coupling a and the name h are pure numbers."""
    dims = {"t": "L", "r": "L", "\\eta": "L", "\\xi": "L", "\\theta": "1", "\\phi": "1", "a": "1", "r_0": "L",
            "r_h": "L", "Q": "L", "h": "1", "k": "1", "M": "L", "q": "L", "\\sigma": "1", "\\lambda": "L**-2"}
    for system in SYSTEMS:
        names = [p.split(" = ")[0] for p in PARAMETERS[system]]
        vm.DIMENSIONS[(METRIC, system)] = {c: dims[c] for c in COORDS[system] + names}
        if system in HELD:
            vm.HELD[(METRIC, system)] = HELD[system]
            vm.RATES[(METRIC, system)] = RATES[system]


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
    if system in HELD:
        a, r = P["a"], X["r"]
        name = HELD[system][0]
        n = P[name]
        # The slope of h is 2h/((1 + a^2) r), and of k 2a^2 k/((1 + a^2) r), so no derivative is left.
        slope = 2 * n / ((1 + a ** 2) * r) * (a ** 2 if name == "k" else 1)

        def reduce(value):
            value = sp.sympify(value)
            while value.has(sp.Derivative):
                for d in sorted(value.atoms(sp.Derivative), key=lambda d: -d.derivative_count):
                    if d.expr != n:
                        raise AssertionError(f"{METRIC}: {d} is no derivative of {name}")
                    value = value.xreplace({d: sp.diff(slope, r, d.derivative_count - 1)})
                value = value.doit()
            return vm.norm(value)

        out["printer"] = {"lead": [a, r] + [P[x] for x in ("r_h", "M", "q", "sigma", "lambda") if x in P] + [n],
                          "rising": [a]}
        out["reduce"] = reduce
    if system in LIOUVILLE:
        U = LIOUVILLE_U[system]
        out["ricci_scalar"] = LIOUVILLE_SCALAR[system]
        out["components"] = {
            "metric_components": {("t", "t"): "-" + U, ("r", "r"): "\\dfrac{1}{" + U + "}"},
            "inverse_metric_components": {("t", "t"): "-\\dfrac{1}{" + U + "}", ("r", "r"): U}}
    if system == "chan_horne_mann":
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


# (a, r_0, the chart's third and fourth parameter, and the constant chosen at the point: the
# dilaton's phi_0 for b = 1/a and b = -a, the charge Q for b = -1/a), then (r, theta).
LIOUVILLE_POINTS = {
    "liouville_inverse": [  # (a, r_0, M, q, phi_0)
        ((sp.Rational(1, 2), sp.Rational(6, 5), sp.Rational(3, 10), sp.Rational(1, 7), sp.Rational(1, 5)), (sp.Rational(7, 3), sp.Rational(9, 10))),
        ((sp.Rational(3, 2), sp.Rational(2, 3), sp.Rational(1, 4), sp.Rational(2, 5), sp.Rational(-1, 3)), (sp.Rational(5, 2), sp.Rational(1, 3))),
        ((sp.Rational(2, 5), sp.Rational(7, 4), sp.Rational(5, 9), sp.Rational(1, 3), sp.Rational(3, 7)), (sp.Rational(9, 4), sp.Rational(6, 5))),
    ],
    "liouville_minus_a": [  # (a, r_0, r_h, sigma, phi_0)
        ((sp.Rational(1, 2), sp.Rational(6, 5), sp.Rational(3, 2), sp.Rational(5, 4), sp.Rational(1, 5)), (sp.Rational(7, 3), sp.Rational(9, 10))),
        ((sp.Rational(3, 2), sp.Rational(2, 3), sp.Rational(1, 4), sp.Rational(1, 2), sp.Rational(-1, 3)), (sp.Rational(5, 2), sp.Rational(1, 3))),
        ((sp.Integer(2), sp.Rational(7, 4), sp.Rational(5, 9), sp.Rational(2, 3), sp.Rational(3, 7)), (sp.Rational(9, 4), sp.Rational(6, 5))),
    ],
    "liouville_minus_inverse": [  # (a, r_0, r_h, lambda, Q)
        ((sp.Rational(1, 2), sp.Rational(6, 5), sp.Rational(3, 2), sp.Rational(-1, 3), sp.Rational(4, 7)), (sp.Rational(7, 3), sp.Rational(9, 10))),
        ((sp.Integer(1), sp.Rational(5, 4), sp.Rational(2, 3), sp.Rational(1, 9), sp.Rational(3, 5)), (sp.Rational(19, 4), sp.Rational(1, 3))),
        ((sp.Rational(3, 2), sp.Rational(2, 3), sp.Rational(1, 4), sp.Rational(2, 5), sp.Rational(9, 7)), (sp.Rational(5, 2), sp.Rational(6, 5))),
    ],
}


def check_liouville(chart, system):
    """Each family against Chan, Horne and Mann's own form, (4.2), (4.6) or (4.9) with the constants
    of the point, and against the field equations with V = 2 Lambda e^{2 b phi} and the electric
    field F_tr = Q e^{2 a phi}/R^2 of their (2.3)."""
    s = chart.symbols
    t, r, th, ph = s
    P = chart.reader.parameters
    a, r0 = P["a"], P["r_0"]
    name = HELD[system][0]
    power = 2 * a ** 2 / (1 + a ** 2) if name == "k" else 2 / (1 + a ** 2)
    written = {P[name]: (r / r0) ** power}
    g = chart.geo.g.subs(written)
    ricci = [[chart.geo.ricci_ll()[i][j].subs(written) for j in range(4)] for i in range(4)]
    # The reader keys a Greek name without its backslash.
    third, fourth = [P[x.lstrip("\\")] for x in PARAMETERS[system][2:4]]
    for (av, r0v, cv, dv, ev), (rv, thv) in LIOUVILLE_POINTS[system]:
        if system == "liouville_inverse":
            M, q, phi0 = cv, dv, ev
            gamma = r0v ** (av ** 2 / (1 + av ** 2))
            Q = q * gamma * sp.exp(-av * phi0)
            Lam = -av ** 2 * sp.exp(-2 * phi0 / av) / (gamma ** 2 * (1 - av ** 2))
            b = 1 / av
            U = r ** (2 * av ** 2 / (1 + av ** 2)) * ((1 + av ** 2) / ((1 - av ** 2) * gamma ** 2)
                                                      - 2 * (1 + av ** 2) * M / (gamma ** 2 * r)
                                                      + Q ** 2 * (1 + av ** 2) * sp.exp(2 * av * phi0) / (gamma ** 4 * r ** 2))
            R = gamma * r ** (1 / (1 + av ** 2))
            phi = phi0 - av / (1 + av ** 2) * sp.log(r)
        elif system == "liouville_minus_a":
            rh, sigma, phi0 = cv, dv, ev
            gamma = r0v ** (1 / (1 + av ** 2))
            square = (1 + sigma * (1 - av ** 2) / (1 + av ** 2)) / 2      # Q^2 e^{2 a phi_0}/gamma^2
            Q = gamma * sp.exp(-av * phi0) * sp.sqrt(square)
            M = av ** 2 * sigma * rh / (2 * (1 + av ** 2))
            Lam = (sp.exp(2 * av * phi0) / ((1 - av ** 2) * gamma ** 2)
                   - (1 + av ** 2) * Q ** 2 * sp.exp(4 * av * phi0) / ((1 - av ** 2) * gamma ** 4))
            b = -av
            U = r ** (2 / (1 + av ** 2)) * ((1 + av ** 2) / ((1 - av ** 2) * gamma ** 2)
                                            * (-1 + 2 * Q ** 2 * sp.exp(2 * av * phi0) / gamma ** 2)
                                            - 2 * (1 + av ** 2) * M / (av ** 2 * gamma ** 2 * r))
            R = gamma * r ** (av ** 2 / (1 + av ** 2))
            phi = phi0 + av / (1 + av ** 2) * sp.log(r)
            # Their (4.7): the horizon.
            small([2 * M * (1 - av ** 2) / (av ** 2 * (-1 + 2 * Q ** 2 * sp.exp(2 * av * phi0) / gamma ** 2)) - rh],
                  system, "r_h is not their (4.7)")
        else:
            rh, lam, Q = cv, dv, ev
            gamma = r0v ** (1 / (1 + av ** 2))
            Lam = lam * (1 - 3 * av ** 2) * av ** 2 / ((1 + av ** 2) ** ((2 * av ** 2 + 1) / av ** 2)
                                                       * (Q / r0v) ** (2 / av ** 2))
            M = av ** 2 * rh / (2 * (1 + av ** 2))
            phi0 = -sp.log(Q ** 2 * (1 + av ** 2) / gamma ** 2) / (2 * av)
            b = -1 / av
            U = r ** (2 / (1 + av ** 2)) * (1 / gamma ** 2 - 2 * (1 + av ** 2) * M / (av ** 2 * gamma ** 2 * r)
                                            + (1 + av ** 2) ** ((2 * av ** 2 + 1) / av ** 2) * Lam * Q ** (2 / av ** 2)
                                            / ((1 - 3 * av ** 2) * av ** 2 * gamma ** (2 / av ** 2))
                                            * r ** (-2 * (1 - av ** 2) / (1 + av ** 2)))
            R = gamma * r ** (av ** 2 / (1 + av ** 2))
            phi = phi0 + av / (1 + av ** 2) * sp.log(r)
        values = {a: av, r0: r0v, third: cv, fourth: dv}
        at = {**values, r: rv, th: thv, t: sp.Rational(1, 3), ph: sp.Rational(2, 5)}
        where = f"{system} at (r, theta) = {(rv, thv)}, parameters {(av, r0v, cv, dv)}, constant {ev}"
        source = sp.diag(-U, 1 / U, R ** 2, R ** 2 * sp.sin(th) ** 2)
        gv = g.subs(values)
        small(list((gv - source).subs(at)), where, "the metric is not Chan, Horne and Mann's")
        V = 2 * Lam * sp.exp(2 * b * phi)
        F = two_form(Q * sp.exp(2 * av * phi) / R ** 2, 0, 1)
        rc = [[ricci[i][j].subs(values) for j in range(4)] for i in range(4)]
        small(field_equations(gv, rc, s, at, phi, F, av, V, 2 * b * V), where,
              "the field equations with the potential do not hold")


def check(chart, system):
    if system in LIOUVILLE:
        return check_liouville(chart, system)
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
