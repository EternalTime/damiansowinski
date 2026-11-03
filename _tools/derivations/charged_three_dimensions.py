"""The static field of a point charge in three dimensions, in three charts.

The static chart is the line element of Gott, Simon and Alpert's abstract (1986), as Cataldo and
Salgado (1996) print it and as Carlip's charged black hole (1995) reduces to it with no cosmological
constant,

    ds^2 = -(kQ^2/2pi) ln(r_c/r) c^2dt^2 + 2pi dr^2/(kQ^2 ln(r_c/r)) + r^2 dphi^2,

with Einstein's equations G = kT, the stress tensor T = (F F - g F^2/4)/2pi of the charge's field
E = Q/r, and r_c a constant of integration. The outgoing chart is Podolsky and Papajcik's
Robinson-Trautman form (2022) with alpha = 0 and Lambda = 0, -f c^2du^2 - 2c du dr + r^2 dphi^2, and
the ingoing chart Mao and Zhao's (2024) with omega = 0 and no cosmological constant,
-f c^2dv^2 + 2c dv dr + r^2 dphi^2, f = (kQ^2/2pi) ln(r_c/r), the first reached from the static chart
by c du = c dt - dr/f and the second by c dv = c dt + dr/f. charged_three_dimensions.md beside this
file records the sources and the map of each author's constants onto k, Q and r_c.

Every check is made at points: each chart's metric against the source form, each null chart
against the static one pulled back, and each chart against Einstein's equations with the field of
the charge, G_ab - k T_ab = 0, from the chart's own Einstein tensor.
"""
import sympy as sp

METRIC = "charged_three_dimensions"
SYSTEMS = ("static", "eddington_finkelstein_outgoing", "eddington_finkelstein_ingoing")
PARAMETERS = ["k", "Q", "r_c"]
LN = "\\ln\\left(\\dfrac{r_c}{r}\\right)"
F_TEXT = f"\\dfrac{{kQ^2}}{{2\\pi}}{LN}"
PHI = ["\\phi \\in [0, 2\\pi)"]
TIME = {"static": "t", "eddington_finkelstein_outgoing": "u", "eddington_finkelstein_ingoing": "v"}
NAMES = {"static": "Static", "eddington_finkelstein_outgoing": "Outgoing Eddington-Finkelstein",
         "eddington_finkelstein_ingoing": "Ingoing Eddington-Finkelstein"}
# The sign of the cross term c dw dr, -1 for the retarded time and +1 for the advanced one.
CROSS = {"eddington_finkelstein_outgoing": -1, "eddington_finkelstein_ingoing": 1}
# Sample points (k Q^2, r_c, r), inside the horizon, outside it and on either side close to it.
POINTS = ((sp.Rational(1, 3), 2, sp.Rational(1, 2)), (sp.Rational(2, 5), 3, 5),
          (sp.Rational(3, 7), 1, sp.Rational(99, 100)), (sp.Rational(1, 2), sp.Rational(5, 2), 7))


def line_element(system, c="c"):
    """The line element as the page prints it, with c; c = "" gives the chart's own, in x^0 = cw."""
    w = TIME[system]
    if system == "static":
        return (f"ds^2 = -{F_TEXT}{c}{'^2' if c else ''}dt^2 + \\dfrac{{2\\pi\\,dr^2}}{{kQ^2{LN}}}"
                " + r^2d\\phi^2")
    sign = "-" if CROSS[system] < 0 else "+"
    return (f"ds^2 = -{F_TEXT}{c}{'^2' if c else ''}d{w}^2 {sign} 2{c + chr(92) + ',' if c else ''}d{w}\\,dr"
            " + r^2d\\phi^2")


def domains(system):
    w = TIME[system]
    return [f"{w} \\in (-\\infty, \\infty)", "r \\in (0, \\infty)"] + PHI + [
        "k > 0,\\quad r_c > 0", "r = r_c \\;\\text{(the horizon)}"]


def register_dimensions(vm):
    """t, u and v are times; k carries T^2/(M L^2) and Q carries sqrt(M) L/T, so k Q^2 is a number."""
    for system in SYSTEMS:
        vm.DIMENSIONS[(METRIC, system)] = {TIME[system]: "T", "r": "L", "\\phi": "1",
                                           "k": "T**2/(M*L**2)", "Q": "M**(1/2)*L/T", "r_c": "L"}


def spec(system):
    """Every value is a rational function of r, k, Q and ln(r_c/r): pretty writes the logarithm as the
    placeholder L before factoring, and the printer writes L back as the line element does."""
    import verify_metrics as vm
    coords = [TIME[system], "r", "\\phi"]
    probe = vm.Reader(coords, PARAMETERS, ())
    r, rc = probe.symbol["r"], probe.parameters["r_c"]
    k, Q = probe.parameters["k"], probe.parameters["Q"]
    L = sp.Symbol("_L", positive=True)

    def pretty(value):
        v = sp.sympify(value).subs(sp.log(rc / r), L)
        if v.atoms(sp.log):
            raise ValueError(f"{METRIC}/{system}: a logarithm other than ln(r_c/r) in {value}")
        return sp.factor(v)

    line = line_element(system)
    return {
        "metric_id": METRIC,
        "system": {"id": system, "name": NAMES[system], "coords": coords,
                   "domains": domains(system), "parameters": PARAMETERS, "line_element": line},
        "chart_line_element": line_element(system, c=""),
        "printer": {"lead": [r, L, Q, k, rc], "factors": [k, Q, sp.pi, r, L, rc], "overrides": {sp.pi: "\\pi", L: LN}},
        "pretty": pretty,
        "check": lambda chart: check(chart, system),
    }


def register(module):
    module.CHARTS[METRIC] = [lambda system=system: spec(system) for system in SYSTEMS]


# -- Checks at points, independent of the printed tensors ---------------------------------------

def source_metric(system, kq2, rc, r):
    """The source form at a point, written out by hand in each chart's coordinates (w, r, phi)."""
    f = kq2 / (2 * sp.pi) * sp.log(sp.nsimplify(rc) / r)
    if system == "static":
        return sp.Matrix([[-f, 0, 0], [0, 1 / f, 0], [0, 0, r ** 2]])
    s = CROSS[system]
    return sp.Matrix([[-f, s, 0], [s, 0, 0], [0, 0, r ** 2]])


def pulled_back(system, kq2, rc, r):
    """The static metric pulled back through c dw = c dt + s dr/f, s the sign of the cross term:
    dt = dw - s dr/f, so the Jacobian d(t, r, phi)/d(w, r, phi) has the entry -s/f."""
    s = CROSS[system]
    f = kq2 / (2 * sp.pi) * sp.log(sp.nsimplify(rc) / r)
    jacobian = sp.Matrix([[1, -s / f, 0], [0, 1, 0], [0, 0, 1]])
    return jacobian.T * source_metric("static", kq2, rc, r) * jacobian


def maxwell(chart, at):
    """The two form of the charge's field, F = (Q/r) dr ^ c dw, F_rw = Q/r, in every chart: it is the
    two form Podolsky and Papajcik print for their retarded time, and since c dw differs from c dt by
    a multiple of dr, the same two form in the static and the advanced chart."""
    r = chart.symbols[1]
    Q = chart.reader.parameters["Q"]
    F = sp.zeros(3)
    F[1, 0] = Q / r
    F[0, 1] = -Q / r
    return F.subs(at)


def check(chart, system):
    w, r, phi = chart.symbols
    p = chart.reader.parameters
    k, Q, rc = p["k"], p["Q"], p["r_c"]
    einstein = chart.geo.einstein_ll()
    for kq2, rcv, rv in POINTS:
        kv, Qv = sp.Rational(1, 2), sp.sqrt(2 * kq2)
        at = {w: sp.Rational(3, 10), r: rv, phi: 1, k: kv, Q: Qv, rc: rcv}
        g = chart.geo.g.subs(at)
        expected = source_metric(system, kq2, rcv, rv)
        if max(abs(complex(x)) for x in (g - expected).evalf(30)) > 1e-25:
            raise AssertionError(f"{METRIC}/{system}: the metric differs from its source at r = {rv}")
        if system != "static":
            if max(abs(complex(x)) for x in (expected - pulled_back(system, kq2, rcv, rv)).evalf(30)) > 1e-25:
                raise AssertionError(f"{METRIC}/{system}: not the static chart pulled back at r = {rv}")
        gi = g.inv()
        F = maxwell(chart, at)
        F2 = sum(F[a, b] * F[c, d] * gi[a, c] * gi[b, d] for a in range(3) for b in range(3)
                 for c in range(3) for d in range(3))
        T = sp.Matrix(3, 3, lambda a, b: (sum(F[a, c] * F[b, d] * gi[c, d] for c in range(3) for d in range(3))
                                          - g[a, b] * F2 / 4) / (2 * sp.pi))
        residual = (sp.Matrix(3, 3, lambda a, b: einstein[a][b]).subs(at) - kv * T).evalf(30)
        if max(abs(complex(x)) for x in residual) > 1e-25:
            raise AssertionError(f"{METRIC}/{system}: G - kT does not vanish at r = {rv}")
        # Maxwell's equations: the divergence of sqrt(-g) F^ab vanishes, from the chart's metric.
        Fu = chart.geo.ginv * maxwell(chart, {}) * chart.geo.ginv
        root = sp.sqrt(-chart.geo.g.det())
        for a in range(3):
            div = sum(sp.diff(root * Fu[a, b], chart.symbols[b]) for b in range(3))
            if abs(complex(div.subs(at).evalf(30))) > 1e-25:
                raise AssertionError(f"{METRIC}/{system}: Maxwell's equations fail at r = {rv}")
