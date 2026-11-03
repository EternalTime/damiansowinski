"""The Atiyah-Hitchin space, the moduli space of two slowly moving magnetic monopoles: a Riemannian
space of four dimensions with no time, so no coordinate is scaled by c.

ds^2 = f^2 dr^2 + a^2 sigma_1^2 + b^2 sigma_2^2 + c^2 sigma_3^2, with the left invariant forms
sigma_1 = sin psi d theta - cos psi sin theta d phi, sigma_2 = -cos psi d theta - sin psi sin theta
d phi and sigma_3 = d psi + cos theta d phi of Hanany and Pioline (2000), section 3.1, and the
three equations of self duality 2bc a' = f((b - c)^2 - a^2) and their cyclic permutations
(Ghezelbash and Mann 2004, `conditions`; Tong 2005, the subsection "The Moduli Space").

The modulus chart is Gibbons and Manton's solution of 1986 in the modulus k = sin(beta/2) of
their angle beta, f = -b/r, r = 2nK (Weinberg and Yi 2007, the subsection on two identical monopoles;
Ghezelbash and Mann 2004, `abc`, `GAMMA`, `gama`): bc = -4n^2 K E, ca = -4n^2 K (E - (1 - k^2) K) and
ab = 4n^2 K (K - E), with K and E the complete elliptic integrals at the parameter k^2, held as
functions of k with their rates, so that every value is a rational function of k, K and E and
the Ricci tensor is exactly zero. The radius chart is Gibbons and Manton's radius r with a, b
and c free functions of it and no component assuming the equations of self duality, which its
check holds at points. atiyah_hitchin.md beside this file records each chart's source and why
Hanany and Pioline's chart in t is not printed."""
import mpmath
import sympy as sp

METRIC = "atiyah_hitchin"
SYSTEMS = ("modulus", "radius")

# The squares of the three forms, sigma_2 written with its sign taken out.
FORMS = ("\\left(\\sin\\psi\\,d\\theta - \\cos\\psi\\sin\\theta\\,d\\phi\\right)^2",
         "\\left(\\cos\\psi\\,d\\theta + \\sin\\psi\\sin\\theta\\,d\\phi\\right)^2",
         "\\left(d\\psi + \\cos\\theta\\,d\\phi\\right)^2")
ANGLES = ["\\theta \\in [0, \\pi]", "\\phi \\in [0, 2\\pi)", "\\psi \\in [0, \\pi)"]
# 1 - k^2, the square of the complementary modulus.
KP2 = "\\left(1 - k^2\\right)"
A2 = f"\\dfrac{{4n^2K\\left(K - E\\right)\\left(E - {KP2}K\\right)}}{{E}}"
B2 = f"\\dfrac{{4n^2KE\\left(K - E\\right)}}{{E - {KP2}K}}"
C2 = f"\\dfrac{{4n^2KE\\left(E - {KP2}K\\right)}}{{K - E}}"
GKK = f"\\dfrac{{4n^2E\\left(K - E\\right)\\left(E - {KP2}K\\right)}}{{k^2{KP2}^2K}}"
MODULUS_PARAMETERS = ["n", "K = \\mathrm{K}\\left(k^2\\right)", "E = \\mathrm{E}\\left(k^2\\right)"]
RADIUS_PARAMETERS = ["n", "a = a(r)", "b = b(r)", "c = c(r)"]
# dK/dk and dE/dk, from dK/dm = (E - (1 - m)K)/(2m(1 - m)) and dE/dm = (E - K)/(2m) at m = k^2.
MODULUS_RATES = {"K": {"k": f"\\dfrac{{E - {KP2}K}}{{k{KP2}}}"}, "E": {"k": "\\dfrac{E - K}{k}"}}


def parameters(system):
    return MODULUS_PARAMETERS if system == "modulus" else RADIUS_PARAMETERS


def register_dimensions(vm):
    """The coordinates and parameters of each chart, its held names and their rates."""
    vm.DIMENSIONS[(METRIC, "modulus")] = {"k": "1", "\\theta": "1", "\\phi": "1", "\\psi": "1",
                                          "n": "L", "K": "1", "E": "1"}
    vm.DIMENSIONS[(METRIC, "radius")] = {"r": "L", "\\theta": "1", "\\phi": "1", "\\psi": "1",
                                         "n": "L", "a": "L", "b": "L", "c": "L"}
    vm.HELD[(METRIC, "modulus")] = ("K", "E")
    vm.RATES[(METRIC, "modulus")] = MODULUS_RATES


def line_element(system):
    if system == "modulus":
        return (f"ds^2 = {GKK}dk^2 + {A2}{FORMS[0]} + {B2}{FORMS[1]} + {C2}{FORMS[2]}")
    return (f"ds^2 = \\dfrac{{b^2}}{{r^2}}dr^2 + a^2{FORMS[0]} + b^2{FORMS[1]} + c^2{FORMS[2]}")


def domains(system):
    if system == "modulus":
        first, core, point = "k \\in [0, 1)", "k = 0", "k"
    else:
        first, core, point = "r \\in [\\pi n, \\infty)", "r = \\pi n", "r"
    return [first] + ANGLES + [
        f"({point}, \\theta, \\phi, \\psi) \\sim ({point}, \\pi - \\theta, \\phi + \\pi, -\\psi)"
        " \\;\\text{(the exchange of the two monopoles)}",
        f"{core} \\;\\text{{(the core, a real projective plane)}}"]


def spec(system):
    import chart_printer as cp
    import verify_metrics as vm
    coords = ["k" if system == "modulus" else "r", "\\theta", "\\phi", "\\psi"]
    held = vm.HELD.get((METRIC, system), ())
    rates = vm.RATES.get((METRIC, system))
    probe = vm.Reader(coords, parameters(system), (), held=held, rates=rates)
    x, th, psi = probe.symbol[coords[0]], probe.symbol["\\theta"], probe.symbol["\\psi"]
    trig = (sp.sin(th), sp.cos(th)), (sp.sin(psi), sp.cos(psi))

    def collect(poly, printer):
        return cp.collect_nested(poly, trig[0], trig[1], printer)

    line = line_element(system)
    if system == "modulus":
        K, E = probe.parameters["K"], probe.parameters["E"]
        printer = {"lead": [K, E, x, probe.parameters["n"], *trig[1], *trig[0]], "collect": collect}
        extra = {"reduce": lambda value: vm.norm(probe.by_rates(value))}
    else:
        abc = [probe.parameters[s] for s in "abc"]
        printer = {"primed": ["a", "b", "c"], "lead": abc + [x, *trig[1], *trig[0]], "collect": collect}
        extra = {}
    return {
        "metric_id": METRIC,
        "system": {"id": system,
                   "name": "Elliptic Modulus" if system == "modulus" else "Gibbons-Manton Radius",
                   "coords": coords, "domains": domains(system), "parameters": parameters(system),
                   "line_element": line},
        "chart_line_element": line,
        "printer": printer,
        "pretty": lambda v: cp.merge_squares(sp.factor(v)),
        "check": lambda chart: check(chart, system),
        **extra,
    }


# -- The solution, evaluated independently of the printed chart ---------------------------------

def products(k, n=1):
    """Weinberg and Yi's three products at the modulus k = sin(beta/2), from their r dr/dbeta:
    ab = -sin(beta) r r_beta + (1 - cos beta) r^2/2, bc = -sin(beta) r r_beta - (1 + cos beta) r^2/2,
    ca = -sin(beta) r r_beta, with r = 2nK(sin(beta/2)) differentiated numerically."""
    k = mpmath.mpf(k)
    beta = 2 * mpmath.asin(k)
    radius = lambda b: 2 * n * mpmath.ellipk(mpmath.sin(b / 2) ** 2)  # noqa: E731
    r, slope = radius(beta), mpmath.diff(radius, beta)
    s, c = mpmath.sin(beta), mpmath.cos(beta)
    return (-s * r * slope + (1 - c) * r ** 2 / 2, -s * r * slope - (1 + c) * r ** 2 / 2, -s * r * slope)


def lengths(k, n=1):
    """a > 0, b > 0 and c < 0 from the three products, as Hanany and Pioline give the signs."""
    ab, bc, ca = products(k, n)
    return mpmath.sqrt(ab * ca / bc), mpmath.sqrt(ab * bc / ca), -mpmath.sqrt(bc * ca / ab)


def rates(a, b, c, r):
    """The first derivatives along r that self duality with f = -b/r gives a, b and c."""
    return ((a ** 2 - (b - c) ** 2) / (2 * c * r),
            b * (b ** 2 - (c - a) ** 2) / (2 * c * a * r),
            b * (c ** 2 - (a - b) ** 2) / (2 * a * b * r))


def check(chart, system):
    """Point checks only, at forty digits: the modulus chart against Weinberg and Yi's products and
    against a vanishing Ricci tensor; the radius chart against a vanishing Ricci tensor once its
    functions are carried along by the equations of self duality, and against the modulus chart
    pulled back along r = 2nK."""
    mpmath.mp.dps = 40
    R = chart.reader
    symbols = chart.symbols
    points = [("0.31", "0.7", "0.4", "0.3"), ("0.62", "1.9", "2.2", "1.1"), ("0.93", "1.2", "5.0", "2.6")]
    if system == "modulus":
        K, E, n = R.parameters["K"], R.parameters["E"], R.parameters["n"]
        ricci = [R.by_rates(v) for v in sp.Matrix(chart.geo.ricci_ll())]
        for kv, tv, pv, sv in points:
            kv = mpmath.mpf(kv)
            names = {K: sp.Float(mpmath.ellipk(kv ** 2), 40), E: sp.Float(mpmath.ellipe(kv ** 2), 40)}
            at = {symbols[0]: sp.Float(kv, 40), symbols[1]: sp.Rational(tv), symbols[2]: sp.Rational(pv),
                  symbols[3]: sp.Rational(sv), n: 1}
            value = lambda v: complex(sp.N(sp.sympify(v).xreplace(names).xreplace(at), 40))  # noqa: E731
            a, b, c = lengths(kv)
            radius = lambda q: 2 * mpmath.ellipk(q ** 2)  # noqa: E731
            r, slope = radius(kv), mpmath.diff(radius, kv)
            if abs(value(chart.geo.g[0, 0]) - complex((b * slope / r) ** 2)) > 1e-30:
                raise AssertionError(f"atiyah_hitchin modulus: g_kk at k = {kv} is not (b dr/dk / r)^2")
            if abs(value(chart.geo.g[3, 3]) - complex(c ** 2)) > 1e-30:
                raise AssertionError(f"atiyah_hitchin modulus: g_psipsi at k = {kv} is not c^2")
            worst = max(abs(value(v)) for v in ricci)
            if worst > 1e-25:
                raise AssertionError(f"atiyah_hitchin modulus: the Ricci tensor at k = {kv} is {worst}")
        return
    a, b, c = (R.parameters[s] for s in "abc")
    n, r = R.parameters["n"], symbols[0]
    ricci = sp.Matrix(chart.geo.ricci_ll())
    for kv, tv, pv, sv in points:
        kv = mpmath.mpf(kv)
        av, bv, cv = lengths(kv)
        rv = 2 * mpmath.ellipk(kv ** 2)
        first = rates(av, bv, cv, rv)
        # The second derivatives, the first carried along once more.
        step = mpmath.mpf("1e-12")
        ahead = rates(*(q + step * dq for q, dq in zip((av, bv, cv), first)), rv + step)
        behind = rates(*(q - step * dq for q, dq in zip((av, bv, cv), first)), rv - step)
        second = [(p - m) / (2 * step) for p, m in zip(ahead, behind)]
        seconds = {sp.Derivative(f, (r, 2)): sp.Float(d2, 40) for f, d2 in zip((a, b, c), second)}
        firsts = {sp.Derivative(f, r): sp.Float(d1, 40) for f, d1 in zip((a, b, c), first)}
        functions = {f: sp.Float(v, 40) for f, v in zip((a, b, c), (av, bv, cv))}
        at = {r: sp.Float(rv, 40), symbols[1]: sp.Rational(tv), symbols[2]: sp.Rational(pv),
              symbols[3]: sp.Rational(sv), n: 1}
        value = lambda v: complex(sp.N(sp.sympify(v).xreplace(seconds).xreplace(firsts)  # noqa: E731
                                       .xreplace(functions).xreplace(at), 40))
        worst = max(abs(value(v)) for v in ricci)
        if worst > 1e-15:
            raise AssertionError(f"atiyah_hitchin radius: the Ricci tensor at r = {rv} is {worst}")
        grr = value(chart.geo.g[0, 0])
        # Pulled back along r = 2nK(k), g_rr dr^2 is g_kk dk^2 of the modulus chart.
        radius = lambda q: 2 * mpmath.ellipk(q ** 2)  # noqa: E731
        modulus_gkk = (bv * mpmath.diff(radius, kv) / rv) ** 2
        if abs(grr * complex(mpmath.diff(radius, kv)) ** 2 - complex(modulus_gkk)) > 1e-25:
            raise AssertionError(f"atiyah_hitchin radius: g_rr at r = {rv} is not the modulus chart's")


def register(pc):
    pc.CHARTS[METRIC] = [lambda system=system: spec(system) for system in SYSTEMS]
