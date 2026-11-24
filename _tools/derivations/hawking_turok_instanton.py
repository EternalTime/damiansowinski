"""The instanton of Hawking and Turok, its two Lorentzian continuations, and two singular instantons
of the same kind in closed form.

instanton        Hawking and Turok, Phys. Lett. B 425, 25 (1998), hep-th/9802030, eq. (emetric):
                 ds^2 = d sigma^2 + b^2(sigma)(d psi^2 + sin^2 psi d Omega_2^2), Riemannian, with b free,
                 b ~ sigma at the regular pole and b ~ (sigma_f - sigma)^{1/3} at the singular one, and
                 the field equations (field), phi'' + 3(b'/b) phi' = V_phi and b'' = -(8 pi G/3) b(phi'^2 + V).
bubble_exterior  the same, eq. (twometric), continued with psi = pi/2 + i tau: region II, outside the
                 bubble wall, Lorentzian, with tau a dimensionless time.
open_universe    the same, eq. (onemetric), continued with sigma = i t and tau = i pi/2 + chi: region I,
                 the open universe with a(t) = -i b(i t), Lorentzian; here sigma = ict and x^0 = ct.
exponential      Kirklin, Turok, and Wiseman, Phys. Rev. D 63, 083509 (2001), hep-th/0005062, Sec. VII,
                 eq. (solchi): ds^2 = (2 pi/L)(6/kappa_4 Lambda_4)^{3/2} cos chi (d chi^2 + sin^2 chi
                 d Omega_3^2), with phi = -sqrt(3/2 kappa_4) log((2 pi/L) sqrt(6/kappa_4 Lambda_4) cos chi)
                 in the potential Lambda_4 exp(sqrt(2 kappa_4/3) phi), kappa_4 = 8 pi G; the prefactor is
                 written ell^2, and ell^2 = 6/(kappa_4 V_0) for the potential V_0 at chi = 0. The
                 potential is Garriga's lambda e^{2 kappa phi/3} with kappa = (12 pi G)^{1/2},
                 hep-th/9804106, eq. (expo).
vilenkin         Chiba, Phys. Lett. B 442, 59 (1998), gr-qc/9810002, Sec. II.A: ds^2 = T (mu + T^2)^{-1/2}
                 (dT^2 + (mu + T^2) d Omega_3^2) with phi ~ -ln(T^2/(mu + T^2)); Garriga, hep-th/9804106,
                 eqs. (instanton), (sol) and (solphi) with mu = A^2: b = T^{1/2}(A^2 + T^2)^{1/4},
                 phi = -(3/4 kappa) ln(T^2/(A^2 + T^2)), C = 3A^2/2 kappa; Vilenkin's massless field,
                 phi' = -C/b^3 and b'^2 - 4 pi G C^2/3 b^4 = 1, hep-th/9803084.

Each chart is checked at points, in forty digits, with no symbolic simplification: the instanton and its
continuations against the field equations of Hawking and Turok and against one another, and the de
Sitter space of a field at rest; the two closed forms against Einstein's equations with their scalar
fields, their sources' own relations, and the Ricci scalar 2/(3 s^2) of Turok and Hawking at a proper
distance s from the singularity, hep-th/9803156. hawking_turok_instanton.md beside this file records
each chart's source and every drawing decision.
"""
import mpmath as mp
import sympy as sp

METRIC = "hawking_turok_instanton"
CHARTS = ("instanton", "bubble_exterior", "open_universe", "exponential", "vilenkin")
NAMES = {"instanton": "Instanton", "bubble_exterior": "Outside the Bubble Wall", "open_universe": "Open Universe",
         "exponential": "Exponential Potential", "vilenkin": "Vilenkin's Instanton"}
ANGLES = ["\\theta \\in [0, \\pi]", "\\phi \\in [0, 2\\pi)"]
S2 = "d\\theta^2 + \\sin^2\\theta\\,d\\phi^2"
S3 = "d\\psi^2 + \\sin^2\\psi\\,d\\theta^2 + \\sin^2\\psi\\sin^2\\theta\\,d\\phi^2"
COORDS = {"instanton": ["\\sigma", "\\psi", "\\theta", "\\phi"],
          "bubble_exterior": ["\\tau", "\\sigma", "\\theta", "\\phi"],
          "open_universe": ["t", "\\chi", "\\theta", "\\phi"],
          "exponential": ["\\chi", "\\psi", "\\theta", "\\phi"],
          "vilenkin": ["T", "\\psi", "\\theta", "\\phi"]}
PARAMETERS = {"instanton": ["b = b(\\sigma)", "\\sigma_f"], "bubble_exterior": ["b = b(\\sigma)", "\\sigma_f"],
              "open_universe": ["a = a(t)"], "exponential": ["\\ell"], "vilenkin": ["A"]}
DOMAINS = {
    "instanton": ["\\sigma \\in [0, \\sigma_f)", "\\psi \\in [0, \\pi]"] + ANGLES
    + ["\\sigma = 0 \\;\\text{(the regular pole)}", "\\sigma = \\sigma_f \\;\\text{(the singularity)}"],
    "bubble_exterior": ["\\tau \\in (-\\infty, \\infty)", "\\sigma \\in (0, \\sigma_f)"] + ANGLES
    + ["\\sigma = 0 \\;\\text{(the bubble wall)}", "\\sigma = \\sigma_f \\;\\text{(the singularity)}"],
    "open_universe": ["t \\in (0, \\infty)", "\\chi \\in [0, \\infty)"] + ANGLES
    + ["t = 0 \\;\\text{(the bubble wall)}"],
    "exponential": ["\\chi \\in [0, \\pi/2)", "\\psi \\in [0, \\pi]"] + ANGLES
    + ["\\chi = 0 \\;\\text{(the regular pole)}", "\\chi = \\pi/2 \\;\\text{(the singularity)}"],
    "vilenkin": ["T \\in (0, \\infty)", "\\psi \\in [0, \\pi]"] + ANGLES
    + ["T = 0 \\;\\text{(the singularity)}", "T \\to \\infty \\;\\text{(flat space)}"],
}
LINES = {
    "instanton": "ds^2 = d\\sigma^2 + b^2\\left(" + S3 + "\\right)",
    "bubble_exterior": ("ds^2 = d\\sigma^2 + b^2\\left(-d\\tau^2 + \\cosh^2\\tau\\,d\\theta^2"
                        " + \\cosh^2\\tau\\sin^2\\theta\\,d\\phi^2\\right)"),
    "open_universe": ("ds^2 = -%sdt^2 + a^2\\left(d\\chi^2 + \\sinh^2\\chi\\,d\\theta^2"
                      " + \\sinh^2\\chi\\sin^2\\theta\\,d\\phi^2\\right)"),
    "exponential": ("ds^2 = \\ell^2\\cos\\chi\\left(d\\chi^2 + \\sin^2\\chi\\,d\\psi^2"
                    " + \\sin^2\\chi\\sin^2\\psi\\,d\\theta^2 + \\sin^2\\chi\\sin^2\\psi\\sin^2\\theta\\,d\\phi^2\\right)"),
    "vilenkin": ("ds^2 = \\dfrac{T}{\\sqrt{A^2 + T^2}}\\left(dT^2 + \\left(A^2 + T^2\\right)\\left("
                 + S3 + "\\right)\\right)"),
}
_DIMS = {"\\sigma": "L", "\\psi": "1", "\\theta": "1", "\\phi": "1", "\\tau": "1", "t": "T", "\\chi": "1", "T": "L",
         "b": "L", "\\sigma_f": "L", "a": "L", "\\ell": "L", "A": "L"}
DIMENSIONS = {(METRIC, s): {name.split(" = ")[0]: _DIMS[name.split(" = ")[0]] for name in COORDS[s] + PARAMETERS[s]}
              for s in CHARTS}
DIGITS = 40
TOLERANCE = sp.Float("1e-28")


def register(module):
    module.CHARTS[METRIC] = [lambda system=system: chart_spec(system) for system in CHARTS]


def line(system, c2=""):
    return LINES[system] % c2 if "%s" in LINES[system] else LINES[system]


def chart_spec(system):
    """One of the five charts, its values printed with the free function or the parameter first."""
    import verify_metrics as vm

    probe = vm.Reader(COORDS[system], PARAMETERS[system], ())
    P = probe.parameters
    spec = {"metric_id": METRIC,
            "system": {"id": system, "name": NAMES[system], "coords": COORDS[system], "domains": DOMAINS[system],
                       "parameters": PARAMETERS[system], "line_element": line(system, "c^2")},
            "chart_line_element": line(system),
            "printer": {},
            "check": lambda chart: chart_check(chart, system)}
    if system in ("instanton", "bubble_exterior"):
        spec["printer"] = {"lead": [P["b"]], "primed": ["b"]}
        if system == "bubble_exterior":
            spec["pretty"] = held_hyperbolic(probe.symbol["\\tau"])
    elif system == "open_universe":
        spec["printer"] = {"lead": [P["a"]], "primed": ["a"]}
        spec["pretty"] = held_hyperbolic(probe.symbol["\\chi"])
    elif system == "exponential":
        spec["printer"] = {"lead": [P["ell"]]}
        spec["pretty"] = in_cosine(probe.symbol["\\chi"])
    else:
        spec["printer"] = {"lead": [P["A"], probe.symbol["T"]]}
    return spec


def held_hyperbolic(angle):
    """chart_printer.hyperbolic in `angle`, with the free function and its derivatives held as symbols
    while the hyperbolic functions are reduced, as coleman_de_luccia_pretty in print_charts.py does."""
    import chart_printer as cp
    hyperbolic = cp.hyperbolic(angle)

    def pretty(value):
        value = sp.sympify(value)
        held = {d: sp.Dummy(f"d{i}") for i, d in enumerate(sorted(
            value.atoms(sp.Derivative) | value.atoms(sp.core.function.AppliedUndef), key=sp.default_sort_key))}
        back = {v: k for k, v in held.items()}
        return hyperbolic(value.xreplace(held)).xreplace(back)
    return pretty


def in_cosine(angle):
    """A `pretty` that writes every even power of sin(angle) in cos(angle), cancels the fraction and
    keeps the shortest factoring of each side, all in the cosine or all in the sine, as szekeres_polar
    in print_charts.py does; the conformal factor of the exponential chart is a power of the cosine."""
    c, s = sp.Symbol("_c"), sp.Symbol("_s")
    held = {sp.cos(angle): c, sp.sin(angle): s}
    back = {c: sp.cos(angle), s: sp.sin(angle)}
    circle = s ** 2 + c ** 2 - 1

    def cosine(value):
        return sp.expand(value).xreplace(held).replace(
            lambda e: e.is_Pow and e.base == s and e.exp.is_Integer and e.exp > 1,
            lambda e: (1 - c ** 2) ** (e.exp // 2) * s ** (e.exp % 2))

    def side(polynomial):
        polynomial = sp.expand(polynomial)
        forms = [polynomial, sp.reduced(polynomial, [circle], s, c)[1], sp.reduced(polynomial, [circle], c, s)[1]]
        return min((sp.factor(form) for form in forms), key=sp.count_ops).xreplace(back)

    def pretty(value):
        numerator, denominator = sp.fraction(sp.together(sp.sympify(value)))
        numerator, denominator = sp.fraction(sp.cancel(cosine(numerator) / cosine(denominator)))
        return side(numerator) / side(denominator)
    return pretty


# -- the checks ------------------------------------------------------------------------------------

def _numeric(expr, at):
    return sp.N(sp.sympify(expr).xreplace(at), DIGITS)


def _zero(value, what, tolerance=TOLERANCE):
    if abs(complex(value)) > tolerance:
        raise AssertionError(f"{METRIC}: {what} fails, by {value}")


def _chart(system):
    import chart_printer as cp
    spec = chart_spec(system)
    return cp.Chart(spec["system"]["coords"], spec["system"]["parameters"], spec["chart_line_element"])


def _profile(value, function, profile, variable):
    """A value with the free function replaced by a profile of its variable, derivatives done."""
    return sp.sympify(value).subs(function, profile.subs(sp.Symbol("x", real=True), variable)).doit()


def _jets(value, function, variable, numbers):
    """A value with the free function and its first two derivatives replaced by three numbers."""
    f0, f1, f2 = numbers
    value = sp.sympify(value)
    value = value.subs(sp.Derivative(function, (variable, 2)), f2).subs(sp.Derivative(function, variable), f1)
    return value.subs(function, f0)


def _stress(chart, at, gradient, potential, kappa, what, function=None, variable=None, numbers=None):
    """Einstein's equations G_ab = kappa (d_a phi d_b phi - g_ab ((grad phi)^2/2 + V)) at the point,
    for the gradient d_a phi and the potential V there, of the Euclidean action
    R/2kappa - (grad phi)^2/2 - V, or its Lorentzian form, for kappa = 8 pi G."""
    import verify_metrics as vm

    n = len(chart.symbols)
    g, ginv = chart.geo.g, chart.geo.ginv
    einstein = chart.geo.einstein_ll()

    def here(value):
        if function is not None:
            value = _jets(value, function, variable, numbers)
        return _numeric(value, at)

    square = sum(here(ginv[a, b]) * gradient[a] * gradient[b] for a in range(n) for b in range(n))
    for a in range(n):
        for b in range(a, n):
            T = gradient[a] * gradient[b] - here(g[a, b]) * (square / 2 + potential)
            _zero(here(vm._at(einstein, (a, b))) - kappa * T,
                  f"{what}: Einstein's equation in slot {chart.coords_tex[a]}{chart.coords_tex[b]}")


def _de_sitter(chart, function, profile, variable, ell, at, what):
    """With the radius of a field at rest, R_ab = (3/ell^2) g_ab: the round 4-sphere, or de Sitter space."""
    ricci = sp.Matrix(chart.geo.ricci_ll())
    for i in range(4):
        for j in range(4):
            value = _profile(ricci[i, j] - 3 * chart.geo.g[i, j] / ell ** 2, function, profile, variable)
            _zero(_numeric(value, at), f"{what}: R_ab = 3 g_ab/l^2 in slot {i}{j}")


def chart_check(chart, system):
    {"instanton": _check_instanton, "bubble_exterior": _check_bubble_exterior, "open_universe": _check_open_universe,
     "exponential": _check_exponential, "vilenkin": _check_vilenkin}[system](chart)


KAPPA = sp.Rational(8, 1) * sp.pi * sp.Rational(37, 100)        # 8 pi G, with G = 0.37 in units c = 1
X = sp.Symbol("x", real=True)


def _ht_jets(b0, dphi, V, kappa):
    """HT98 (field) with the constraint it integrates: b'^2 = 1 + (kappa/3) b^2 (phi'^2/2 - V) and
    b'' = -(kappa/3) b (phi'^2 + V), the instanton's radius and its first two derivatives at a point."""
    b1 = sp.sqrt(1 + kappa * b0 ** 2 * (dphi ** 2 / 2 - V) / 3)
    b2 = -kappa * b0 * (dphi ** 2 + V) / 3
    return b0, b1, b2


JETS = [(sp.Rational(7, 10), sp.Rational(3, 5), sp.Rational(1, 4)), (sp.Rational(13, 10), sp.Rational(1, 2), sp.Rational(1, 10)),
        (sp.Rational(1, 3), sp.Rational(9, 5), sp.Rational(2, 3))]


def _check_instanton(chart):
    sigma, psi, theta, phi = chart.symbols
    b = chart.reader.parameters["b"]
    ell = sp.Rational(17, 10)
    # HT98: the O(5) instanton of a field at rest, b = H^{-1} sin(H sigma), is the round 4-sphere.
    for point in ((sp.Rational(1, 2), sp.Rational(7, 10), sp.Rational(11, 10)), (2, sp.Rational(5, 4), sp.Rational(3, 10))):
        at = {sigma: point[0], psi: point[1], theta: point[2]}
        _de_sitter(chart, b, ell * sp.sin(X / ell), sigma, ell, at, "HT98, the O(5) instanton")
    # HT98 (field): with b, b' and b'' of its field equations, G_ab = kappa T_ab of a scalar phi(sigma).
    for b0, dphi, V in JETS:
        at = {sigma: sp.Rational(9, 10), psi: sp.Rational(4, 5), theta: sp.Rational(6, 5)}
        _stress(chart, at, [dphi, 0, 0, 0], V, KAPPA, "HT98 (field)", b, sigma, _ht_jets(b0, dphi, V, KAPPA))
    # HT98: b ~ sigma at the regular pole, and TH98's R = 2/(3 (sigma_f - sigma)^2) where b ~ (sigma_f - sigma)^{1/3}.
    sigma_f, k = sp.Rational(5, 2), sp.Rational(6, 5)
    scalar = chart.geo.ricci_scalar()
    for gap in (sp.Rational(1, 10 ** 12), sp.Rational(1, 10 ** 15)):
        value = _profile(scalar, b, k * (sigma_f - X) ** sp.Rational(1, 3), sigma)
        _zero(_numeric(value * gap ** 2, {sigma: sigma_f - gap, psi: 1, theta: 1}) - sp.Rational(2, 3),
              "TH98, R (sigma_f - sigma)^2 = 2/3 at the singularity", sp.Float("1e-12"))
    value = _profile(scalar * sigma ** 2, b, X + X ** 3 / 5, sigma)
    _zero(_numeric(value, {sigma: sp.Rational(1, 10 ** 15), psi: 1, theta: 1}),
          "HT98, the regular pole, R sigma^2 -> 0 with b = sigma + O(sigma^3)", sp.Float("1e-12"))


def _check_bubble_exterior(chart):
    tau, sigma, theta, phi = chart.symbols
    b = chart.reader.parameters["b"]
    source = _chart("instanton")
    s_sigma, s_psi, s_theta, s_phi = source.symbols
    s_b = source.reader.parameters["b"]
    # HT98 (twometric): psi = pi/2 + i tau, d psi = i d tau, from (emetric), with the profile b = sin x + x^3/7.
    profile = sp.sin(X) + X ** 3 / 7
    g = source.geo.g.subs(s_b, profile.subs(X, s_sigma)).doit().subs(
        {s_sigma: sigma, s_psi: sp.pi / 2 + sp.I * tau, s_theta: theta, s_phi: phi}, simultaneous=True)
    jacobian, order = [1, sp.I, 1, 1], [1, 0, 2, 3]
    for point in ((sp.Rational(3, 10), sp.Rational(7, 10), sp.Rational(11, 10)), (-sp.Rational(6, 5), 2, sp.Rational(3, 10))):
        at = {tau: point[0], sigma: point[1], theta: point[2]}
        for i in range(4):
            for j in range(4):
                value = g[order[i], order[j]] * jacobian[order[i]] * jacobian[order[j]]
                mine = _profile(chart.geo.g[i, j], b, profile, sigma)
                _zero(_numeric(value - mine, at), f"HT98 (twometric) is (emetric) at psi = pi/2 + i tau, slot {i}{j}")
        _de_sitter(chart, b, sp.Rational(17, 10) * sp.sin(X / sp.Rational(17, 10)), sigma, sp.Rational(17, 10), at,
                   "HT98, de Sitter space outside the wall")
    for b0, dphi, V in JETS:
        at = {tau: sp.Rational(2, 5), sigma: sp.Rational(9, 10), theta: sp.Rational(6, 5)}
        _stress(chart, at, [0, dphi, 0, 0], V, KAPPA, "HT98 (field) outside the wall", b, sigma,
                _ht_jets(b0, dphi, V, KAPPA))


def _check_open_universe(chart):
    x0, chi, theta, phi = chart.symbols          # x0 = ct
    a = chart.reader.parameters["a"]
    source = _chart("bubble_exterior")
    s_tau, s_sigma, s_theta, s_phi = source.symbols
    s_b = source.reader.parameters["b"]
    # HT98 (onemetric): sigma = i x0 and tau = i pi/2 + chi, so d sigma = i dx0 and d tau = d chi, with
    # b(i x0) = i a(x0), a(t) = -i b(i t); the profile b = sin x + x^3/7 is a = sinh x - x^3/7.
    for point in ((sp.Rational(3, 10), sp.Rational(7, 10), sp.Rational(11, 10)), (sp.Rational(9, 5), 2, sp.Rational(3, 10))):
        at = {x0: point[0], chi: point[1], theta: point[2]}
        g = source.geo.g.subs(s_b, sp.sin(s_sigma) + s_sigma ** 3 / 7).doit().subs(
            {s_sigma: sp.I * x0, s_tau: sp.I * sp.pi / 2 + chi, s_theta: theta, s_phi: phi}, simultaneous=True)
        jacobian, order = [1, sp.I, 1, 1], [1, 0, 2, 3]
        for i in range(4):
            for j in range(4):
                value = g[order[i], order[j]] * jacobian[order[i]] * jacobian[order[j]]
                mine = _profile(chart.geo.g[i, j], a, sp.sinh(X) - X ** 3 / 7, x0)
                _zero(_numeric(value - mine, at), f"HT98 (onemetric) is (twometric) continued, slot {i}{j}")
        ell = sp.Rational(17, 10)
        _de_sitter(chart, a, ell * sp.sinh(X / ell), x0, ell, at, "HT98, the open de Sitter universe")
    # The Friedmann equations of an open universe with a scalar phi(t): a'^2 = 1 + (kappa/3) a^2 (phi'^2/2 + V)
    # and a'' = -(kappa/3) a (phi'^2 - V), HT98 (field) continued.
    for a0, dphi, V in JETS:
        a1 = sp.sqrt(1 + KAPPA * a0 ** 2 * (dphi ** 2 / 2 + V) / 3)
        a2 = -KAPPA * a0 * (dphi ** 2 - V) / 3
        at = {x0: sp.Rational(2, 5), chi: sp.Rational(9, 10), theta: sp.Rational(6, 5)}
        _stress(chart, at, [dphi, 0, 0, 0], V, KAPPA, "HT98 (field) in the open universe", a, x0, (a0, a1, a2))


def _near_singularity(chart, radius_of, proper_from, coordinate, at_singularity, toward, what):
    """Near the singularity R s^2 -> 2/3 for the proper distance s from it, TH98, and b^3/s tends to a
    constant, b ~ s^{1/3}, HT98; s is integrated in forty digits."""
    scalar = chart.geo.ricci_scalar()
    ratios = []
    for gap in (sp.Rational(1, 10 ** 10), sp.Rational(1, 10 ** 13)):
        x = at_singularity - toward * gap
        s = proper_from(x)
        R = _numeric(scalar, {coordinate: x, **chart_point(chart)})
        _zero(R * s ** 2 - sp.Rational(2, 3), f"{what}: TH98, R s^2 = 2/3 near the singularity", sp.Float("1e-6"))
        ratios.append(_numeric(radius_of(x), {}) ** 3 / s)
    _zero(ratios[0] / ratios[1] - 1, f"{what}: HT98, b ~ s^(1/3) near the singularity", sp.Float("1e-6"))


def chart_point(chart):
    return {chart.symbols[1]: sp.Rational(4, 5), chart.symbols[2]: sp.Rational(6, 5)}


def _check_exponential(chart):
    chi, psi, theta, phi = chart.symbols
    ell = chart.reader.parameters["ell"]
    mp.mp.dps = DIGITS
    for kappa, lam, L in ((KAPPA, sp.Rational(3, 7), sp.Rational(5, 2)), (sp.Rational(2, 1), sp.Rational(5, 3), sp.Rational(3, 1))):
        # KTW00 (solchi): the prefactor and the field, in the potential Lambda_4 exp(sqrt(2 kappa_4/3) phi).
        B = 2 * sp.pi / L * sp.sqrt(6 / (kappa * lam))
        size = 2 * sp.pi / L * (6 / (kappa * lam)) ** sp.Rational(3, 2)
        field = -sp.sqrt(3 / (2 * kappa)) * sp.log(B * sp.cos(chi))
        potential = lam * sp.exp(sp.sqrt(2 * kappa / 3) * field)
        _zero(_numeric(6 / (kappa * potential.subs(chi, 0)) - size, {}), "KTW00 (solchi), l^2 = 6/(kappa V_0)")
        for point in ((sp.Rational(3, 10), sp.Rational(7, 10), sp.Rational(11, 10)), (sp.Rational(3, 2), 2, sp.Rational(3, 10))):
            at = {chi: point[0], psi: point[1], theta: point[2], ell: sp.sqrt(size)}
            gradient = [_numeric(sp.diff(field, x), at) for x in chart.symbols]
            _stress(chart, at, gradient, _numeric(potential, at), kappa, "KTW00 (solchi), Einstein's equations")
            # The wave equation box phi = dV/dphi.
            g, ginv = chart.geo.g, chart.geo.ginv
            root = sp.sqrt(g[0, 0] * g[1, 1] * g[2, 2] * g[3, 3])       # the metric is diagonal
            box = sum(sp.diff(root * ginv[a, b] * sp.diff(field, chart.symbols[b]), chart.symbols[a])
                      for a in range(4) for b in range(4)) / root
            slope = sp.sqrt(2 * kappa / 3) * potential
            _zero(_numeric(box - slope, at), "KTW00 (solchi), the wave equation")
    # The regular pole, where b = l sqrt(cos chi) sin chi grows as the proper distance, and the singularity.
    size = sp.Rational(9, 4)
    radius = lambda x: size * sp.sqrt(sp.cos(x)) * sp.sin(x)
    slope = sp.diff(radius(chi), chi) / (size * sp.sqrt(sp.cos(chi)))
    _zero(_numeric(slope, {chi: sp.Rational(1, 10 ** 12)}) - 1, "KTW00, db/d sigma = 1 at the regular pole",
          sp.Float("1e-20"))
    proper = lambda x: mp.mpf(str(size)) * mp.quad(lambda u: mp.sqrt(mp.cos(u)), [mp.mpf(str(x)), mp.pi / 2])
    _near_singularity(_with(chart, {ell: size}), radius, lambda x: sp.Float(proper(sp.N(x, DIGITS)), DIGITS), chi,
                      sp.pi / 2, 1, "KTW00 (solchi)")


class _with:
    """A chart with its parameters set, for the scalars near a singularity."""

    def __init__(self, chart, values):
        self.symbols = chart.symbols
        self._values = values
        self._chart = chart

    @property
    def geo(self):
        outer = self

        class Geo:
            def ricci_scalar(self):
                return outer._chart.geo.ricci_scalar().xreplace(outer._values)
        return Geo()


def _check_vilenkin(chart):
    T, psi, theta, phi = chart.symbols
    A = chart.reader.parameters["A"]
    mp.mp.dps = DIGITS
    G = sp.Rational(37, 100)
    kappa = 8 * sp.pi * G
    root_kappa = sp.sqrt(12 * sp.pi * G)                        # Garriga's kappa
    for size in (sp.Rational(3, 2), sp.Rational(2, 5)):
        field = -3 / (4 * root_kappa) * sp.log(T ** 2 / (A ** 2 + T ** 2))
        radius = sp.sqrt(T) * (A ** 2 + T ** 2) ** sp.Rational(1, 4)
        charge = 3 * A ** 2 / (2 * root_kappa)
        rate = 1 / sp.sqrt(chart.geo.g[0, 0])                   # d/d sigma = (dT/d sigma) d/dT
        for point in ((sp.Rational(3, 10), sp.Rational(7, 10), sp.Rational(11, 10)), (sp.Rational(5, 2), 2, sp.Rational(3, 10))):
            at = {T: point[0], psi: point[1], theta: point[2], A: size}
            gradient = [_numeric(sp.diff(field, x), at) for x in chart.symbols]
            _stress(chart, at, gradient, 0, kappa, "G98b (solphi), Einstein's equations of a massless field")
            # The radius of the 3-spheres, G98b (sol), and Vilenkin's phi' = -C/b^3 and b'^2 - 4 pi G C^2/3 b^4 = 1.
            _zero(_numeric(chart.geo.g[1, 1] - radius ** 2, at), "G98b (sol), b^2 = g_psi psi")
            _zero(_numeric(rate * sp.diff(field, T) + charge / radius ** 3, at), "V98, phi' = -C/b^3")
            _zero(_numeric((rate * sp.diff(radius, T)) ** 2 - 4 * sp.pi * G * charge ** 2 / (3 * radius ** 4) - 1, at),
                  "V98, b'^2 - 4 pi G C^2/3 b^4 = 1")
            # G98b (instanton) and (ansatz): the four-dimensional metric is e^{-2 kappa phi/3} times the
            # five-dimensional one, whose circle has g_yy = T^2/(A^2 + T^2) = e^{-4 kappa phi/3}.
            weight = sp.exp(-2 * root_kappa * field / 3)
            _zero(_numeric(chart.geo.g[0, 0] - weight, at), "G98b (ansatz), g_TT")
            _zero(_numeric(chart.geo.g[1, 1] - weight * (A ** 2 + T ** 2), at), "G98b (ansatz), g_psi psi")
            _zero(_numeric(sp.exp(-4 * root_kappa * field / 3) - T ** 2 / (A ** 2 + T ** 2), at), "G98b (ansatz), g_yy")
    # Garriga's five-dimensional instanton is Ricci flat, Chiba's Euclidean Schwarzschild solution.
    import verify_metrics as vm
    Y = sp.Symbol("Y", real=True)
    five = sp.diag(1, A ** 2 + T ** 2, (A ** 2 + T ** 2) * sp.sin(psi) ** 2,
                   (A ** 2 + T ** 2) * sp.sin(psi) ** 2 * sp.sin(theta) ** 2, T ** 2 / (A ** 2 + T ** 2))
    ricci = sp.Matrix(vm.Geometry(five, [T, psi, theta, phi, Y], 120).ricci_ll())
    for i in range(5):
        for j in range(5):
            _zero(_numeric(ricci[i, j], {T: sp.Rational(7, 10), psi: sp.Rational(4, 5), theta: sp.Rational(6, 5),
                                         A: sp.Rational(3, 2)}), f"G98b (instanton) is Ricci flat, slot {i}{j}")
    # Flat far away, and TH98's 2/(3 s^2) near the singularity, with s = int_0^T sqrt(g_TT) dT.
    size = sp.Rational(3, 2)
    kretschmann = chart.geo.kretschmann()
    _zero(_numeric(kretschmann, {T: 10 ** 8, psi: 1, theta: 1, A: size}), "C98, flat at large T", sp.Float("1e-20"))
    proper = lambda x: mp.quad(lambda u: mp.sqrt(u / mp.sqrt(mp.mpf(str(size)) ** 2 + u ** 2)), [0, mp.mpf(str(x))])
    radius_of = lambda x: sp.sqrt(x) * (size ** 2 + x ** 2) ** sp.Rational(1, 4)
    _near_singularity(_with(chart, {A: size}), radius_of, lambda x: sp.Float(proper(sp.N(x, DIGITS)), DIGITS), T,
                      sp.Integer(0), -1, "C98")
