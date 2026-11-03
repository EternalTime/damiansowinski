"""The wormhole of Giddings and Strominger (Nucl. Phys. B 306, 890 (1988)): a Riemannian space of four
dimensions with no time, so no coordinate is scaled by c, held open by the field strength H of an
axion, a 3-form threading every 3-sphere with the same flux.

The areal chart is Hebecker, Mikhail and Soler's eq. `GSwormhole` (Front. Astron. Space Sci. 5, 35
(2018), arXiv:1807.00824), ds^2 = (1 + C/r^4)^{-1} dr^2 + r^2 dOmega_3^2 with C = -a^4, one half of the
wormhole, r >= a. The proper distance chart is Giddings and Strominger's own ansatz (2.7),
ds^2 = dt^2 + a^2(t) dOmega_3^2 (their t is rho here and their a(t) is r), with the radius left a free
function of rho and no component assuming their equation (2.12), (dr/drho)^2 = 1 - a^4/r^4; Witten
writes the same, `yogo` and `winto` (arXiv:2601.01587). The conformal chart is Loges, Shiu and Sudhir's
(J. High Energy Phys. 2022(08), 064, arXiv:2203.01956), the Euclidean solution a(r) = sqrt(n cosh 2r)
in conformal time, here a^2 cosh(2 eta)(d eta^2 + dOmega_3^2), which covers both halves.

Every check compares at points in forty digits, no simplify: the areal chart against the curvature
scalar -6a^4/r^6, the mixed Einstein tensor (-3, 3, 3, 3) a^4/r^6 and the mixed Ricci tensor
(-6, 0, 0, 0) a^4/r^6, whose eigenvalues are non-positive with one negative; the proper distance chart
against the same values once r, r' and r'' are carried by the equation of the solution, and against the
areal chart pulled back; the conformal chart against the areal chart pulled back along
r^2 = a^2 cosh 2 eta. giddings_strominger.md beside this file records each chart's source."""
import mpmath
import sympy as sp

METRIC = "giddings_strominger"
SYSTEMS = ("areal", "proper", "conformal")
NAMES = {"areal": "Areal Radius", "proper": "Proper Distance (Giddings-Strominger)",
         "conformal": "Conformal (Loges-Shiu-Sudhir)"}
SPHERE = "\\left(d\\chi^2 + \\sin^2\\chi\\,d\\theta^2 + \\sin^2\\chi\\sin^2\\theta\\,d\\phi^2\\right)"
ANGLES = ["\\chi \\in [0, \\pi]", "\\theta \\in [0, \\pi]", "\\phi \\in [0, 2\\pi)"]
FIRST = {"areal": "r", "proper": "\\rho", "conformal": "\\eta"}
PARAMETERS = {"areal": ["a"], "proper": ["a", "r = r(\\rho)"], "conformal": ["a"]}
POINTS = [("1.37", "0.7", "0.4", "0.3"), ("2.21", "1.9", "2.2", "1.1"), ("4.05", "1.2", "0.9", "2.6")]


def register_dimensions(vm):
    """The coordinates and parameters of each chart: a and r are lengths, and so is the proper
    distance rho; the conformal eta and the three angles are pure numbers."""
    angles = {"\\chi": "1", "\\theta": "1", "\\phi": "1"}
    vm.DIMENSIONS[(METRIC, "areal")] = {"r": "L", **angles, "a": "L"}
    vm.DIMENSIONS[(METRIC, "proper")] = {"\\rho": "L", **angles, "a": "L", "r": "L"}
    vm.DIMENSIONS[(METRIC, "conformal")] = {"\\eta": "1", **angles, "a": "L"}


def line_element(system):
    if system == "areal":
        return f"ds^2 = \\dfrac{{dr^2}}{{1 - \\dfrac{{a^4}}{{r^4}}}} + r^2{SPHERE}"
    if system == "proper":
        return f"ds^2 = d\\rho^2 + r^2{SPHERE}"
    return f"ds^2 = a^2\\cosh(2\\eta)\\left(d\\eta^2 + d\\chi^2 + \\sin^2\\chi\\,d\\theta^2 + \\sin^2\\chi\\sin^2\\theta\\,d\\phi^2\\right)"


def domains(system):
    if system == "areal":
        return ["r \\in [a, \\infty)"] + ANGLES + ["r = a \\;\\text{(the throat, where the two halves meet)}"]
    if system == "proper":
        return ["\\rho \\in (-\\infty, \\infty)"] + ANGLES + ["\\rho = 0 \\;\\text{(the throat)}"]
    return ["\\eta \\in (-\\infty, \\infty)"] + ANGLES + ["\\eta = 0 \\;\\text{(the throat)}"]


def spec(system):
    import chart_printer as cp
    import verify_metrics as vm
    coords = [FIRST[system], "\\chi", "\\theta", "\\phi"]
    probe = vm.Reader(coords, PARAMETERS[system], ())
    x, chi, th = probe.symbol[coords[0]], probe.symbol["\\chi"], probe.symbol["\\theta"]
    a = probe.parameters["a"]
    line = line_element(system)
    if system == "areal":
        # Every value is a rational function of r whose sums are powers of r^4 - a^4, named whole.
        quartic = sp.Symbol("GSquartic", positive=True)
        printer = {"lead": [a, x], "flip": False, "named": {quartic: "r^4 - a^4"}}
        pretty = quartic_pretty(x, a, quartic)
    elif system == "proper":
        printer = {"lead": [probe.parameters["r"], a, x], "flip": False}
        pretty = None
    else:
        printer = {"lead": [a], "flip": False, "arguments": {2 * x: "2\\eta"}}
        pretty = cp.hyperbolic(2 * x)
    out = {
        "metric_id": METRIC,
        "system": {"id": system, "name": NAMES[system], "coords": coords, "domains": domains(system),
                   "parameters": PARAMETERS[system], "line_element": line},
        "chart_line_element": line,
        "printer": printer,
        "check": lambda chart: check(chart, system),
    }
    if pretty is not None:
        out["pretty"] = pretty
    return out


def quartic_pretty(r, a, quartic):
    """A `pretty` that writes each value as a power of r^4 - a^4, named, times the factored rest."""
    divisor = r ** 4 - a ** 4

    def split(poly):
        poly, power = sp.expand(poly), 0
        while poly != 0 and sp.rem(poly, divisor, r) == 0:
            poly, power = sp.quo(poly, divisor, r), power + 1
        return sp.factor(poly) * quartic ** power

    def pretty(value):
        num, den = sp.fraction(sp.cancel(sp.together(sp.sympify(value))))
        return split(num) / split(den)
    return pretty


# -- The solution, evaluated independently of the printed chart ---------------------------------

def areal_values(r, a=1):
    """The curvature scalar, the mixed Einstein tensor and the mixed Ricci tensor at the areal
    radius r, written out from the solution: R = -6a^4/r^6, G = (-3, 3, 3, 3)a^4/r^6 and
    Ricci = (-6, 0, 0, 0)a^4/r^6."""
    q = mpmath.mpf(a) ** 4 / mpmath.mpf(r) ** 6
    return -6 * q, [-3 * q, 3 * q, 3 * q, 3 * q], [-6 * q, 0, 0, 0]


def radius_at(rho, a=1):
    """The areal radius at the proper distance rho from the throat, the root of
    rho = int_a^r dr'/sqrt(1 - a^4/r'^4), Giddings and Strominger's (2.12) integrated."""
    a, rho = mpmath.mpf(a), mpmath.mpf(rho)
    distance = lambda r: mpmath.quad(lambda s: 1 / mpmath.sqrt(1 - a ** 4 / s ** 4), [a, r])  # noqa: E731
    return mpmath.findroot(lambda r: distance(r) - abs(rho), abs(rho) + a)


def height(a=1):
    """The rise of the embedded funnel from the throat to infinity, sqrt(pi) Gamma(1/4) a / (4 Gamma(3/4))."""
    return mpmath.sqrt(mpmath.pi) * mpmath.gamma(mpmath.mpf(1) / 4) * a / (4 * mpmath.gamma(mpmath.mpf(3) / 4))


def check(chart, system):
    """Point checks only, at forty digits."""
    mpmath.mp.dps = 40
    R = chart.reader
    symbols = chart.symbols
    a = R.parameters["a"]
    geo = chart.geo
    ginv = geo.ginv
    ricci = sp.Matrix(geo.ricci_ll())
    scalar = geo.ricci_scalar()
    einstein = sp.Matrix(geo.einstein_ll()) if hasattr(geo, "einstein_ll") else ricci - scalar * geo.g / 2
    mixed_g = ginv * einstein
    mixed_r = ginv * ricci
    for xv, cv, tv, pv in POINTS:
        at = {symbols[1]: sp.Rational(cv), symbols[2]: sp.Rational(tv), symbols[3]: sp.Rational(pv), a: 1}
        if system == "areal":
            r = mpmath.mpf(xv)
            at[symbols[0]] = sp.Float(r, 40)
            value = lambda v: complex(sp.N(sp.sympify(v).xreplace(at), 40))  # noqa: E731
        elif system == "proper":
            r = mpmath.mpf(xv)
            fn = R.parameters["r"]
            rho = symbols[0]
            first, second = mpmath.sqrt(1 - 1 / r ** 4), 2 / r ** 5
            carried = {sp.Derivative(fn, (rho, 2)): sp.Float(second, 40), sp.Derivative(fn, rho): sp.Float(first, 40),
                       fn: sp.Float(r, 40)}
            at[rho] = sp.Float(radius_at_inverse(r), 40)
            value = lambda v: complex(sp.N(sp.sympify(v).xreplace(carried).xreplace(at), 40))  # noqa: E731
            # Pulled back along (dr/drho)^2 = 1 - a^4/r^4, g_rr dr^2 of the areal chart is drho^2.
            if abs(complex(first ** 2 / (1 - 1 / r ** 4)) - value(geo.g[0, 0])) > 1e-30:
                raise AssertionError(f"{METRIC} proper: g_rho rho at r = {r} is not the areal chart's pulled back")
            if abs(value(geo.g[1, 1]) - complex(r ** 2)) > 1e-30:
                raise AssertionError(f"{METRIC} proper: g_chi chi at r = {r} is not r^2")
        else:
            eta = mpmath.acosh(mpmath.mpf(xv)) / 2
            r = mpmath.sqrt(mpmath.cosh(2 * eta))
            at[symbols[0]] = sp.Float(eta, 40)
            value = lambda v: complex(sp.N(sp.sympify(v).xreplace(at), 40))  # noqa: E731
            # Pulled back along r^2 = a^2 cosh 2eta: g_rr (dr/deta)^2 and r^2 of the areal chart.
            dr = mpmath.sinh(2 * eta) / r
            if abs(complex(dr ** 2 / (1 - 1 / r ** 4)) - value(geo.g[0, 0])) > 1e-30:
                raise AssertionError(f"{METRIC} conformal: g_eta eta at eta = {eta} is not the areal chart's pulled back")
            if abs(value(geo.g[1, 1]) - complex(r ** 2)) > 1e-30:
                raise AssertionError(f"{METRIC} conformal: g_chi chi at eta = {eta} is not r^2")
        want_r, want_g, want_ricci = areal_values(r)
        if abs(value(scalar) - complex(want_r)) > 1e-25:
            raise AssertionError(f"{METRIC} {system}: the Ricci scalar at r = {r} is not -6a^4/r^6")
        for i in range(4):
            for j in range(4):
                g_want = want_g[i] if i == j else 0
                ricci_want = want_ricci[i] if i == j else 0
                if abs(value(mixed_g[i, j]) - complex(g_want)) > 1e-25:
                    raise AssertionError(f"{METRIC} {system}: G^{i}_{j} at r = {r} is not the axion's stress")
                if abs(value(mixed_r[i, j]) - complex(ricci_want)) > 1e-25:
                    raise AssertionError(f"{METRIC} {system}: R^{i}_{j} at r = {r} is not (-6, 0, 0, 0)a^4/r^6")


def radius_at_inverse(r, a=1):
    """The proper distance from the throat at the areal radius r, int_a^r dr'/sqrt(1 - a^4/r'^4)."""
    a = mpmath.mpf(a)
    return mpmath.quad(lambda s: 1 / mpmath.sqrt(1 - a ** 4 / s ** 4), [a, mpmath.mpf(r)])


def register(pc):
    pc.CHARTS[METRIC] = [lambda system=system: spec(system) for system in SYSTEMS]
