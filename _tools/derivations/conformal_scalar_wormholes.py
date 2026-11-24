"""The charts of the wormholes of a massless conformally coupled scalar field, Barcelo and Visser,
gr-qc/9908029, Phys. Lett. B 466, 127 (1999). conformal_scalar_wormholes.md beside this file records
each chart's source and every drawing decision.

Barcelo and Visser write the field equations as kappa G = T for the field phi_c with the stress tensor
of their (2.1), traceless, so R = 0 and their (2.2) and (2.3),

    R_ab = (kappa - phi_c^2/6)^{-1} ((2/3) d_a phi_c d_b phi_c - (1/6) g_ab (d phi_c)^2 - (1/3) phi_c nabla_a nabla_b phi_c),
    box phi_c = 0,

with kappa = c^4/(8 pi G). Their solutions are the Janis-Newman-Winicour-Wyman metrics (2.6) times
Omega^2, Omega(r) = alpha_+ f^{sin chi/2 sqrt 3} + alpha_- f^{-sin chi/2 sqrt 3} with f = 1 - 2 eta/r,
their (2.17), and tan(Delta/2) = (alpha_+ - alpha_-)/(alpha_+ + alpha_-).

    isotropic   the branch chi = pi/3 in the isotropic radius rbar of their (3.21), r = rbar (1 + eta/2 rbar)^2,
                their (3.22) with alpha_+ + alpha_- = 1, written rho here, and in Ray's parameters,
                CQG 42, 017002 (2025), section 4 at D = 4, mu = eta/2 and b = tan(Delta/2):
                ds^2 = (1 - b mu/rho)^2 (-(1 + mu/rho)^{-2} c^2 dt^2 + (1 + mu/rho)^2 (d rho^2 + rho^2 dOmega^2)),
                with the field of their (3.23), phi_c = sqrt(6 kappa) (b rho - mu)/(rho - b mu). For b < 0 the
                chart runs from one flat region at rho -> infinity through the throat rho = mu sqrt(-b) to a
                second at rho -> 0.
    three_parameter
                the whole class in their radius r, (2.6) times Omega^2 with (2.17), alpha_+ + alpha_- = 1,
                eta = 2 mu, gamma = cos chi and sigma = sin chi/2 sqrt 3, for chi in [0, pi):
                ds^2 = Omega^2 (-F c^2 dt^2 + dr^2/F + (f/F) r^2 dOmega^2),  f = 1 - 4 mu/r,  F = f^gamma,
                Omega = ((1 + b) f^sigma + (1 - b) f^-sigma)/2, with the field of their (2.14),
                phi_c = sqrt(6 kappa) Psi/Omega, Psi = ((1 + b) f^sigma - (1 - b) f^-sigma)/2.
                F, Omega and Psi are held as names with declared slopes, so every value is a rational
                function of them, of r and of the parameters, and no exponent holding gamma appears in a
                value. gamma and sigma lie on the conic gamma^2 + 12 sigma^2 = 1 and are defined through its
                rational parametrisation by w = tan(chi/2)/2 sqrt 3, gamma = (1 - 12 w^2)/(1 + 12 w^2) and
                sigma = 2w/(1 + 12 w^2), so that the conic holds in every value with no radical: the checker's
                canonical form writes the root of a radicand with a negative leading term as i times the root
                of its negative, which evaluates to minus the root, and sqrt(1 - 12 sigma^2) came back as -gamma.
                `three_parameter_pretty` writes each value back in gamma and sigma, by w = sigma/(1 + gamma).
"""
import json

import sympy as sp

METRIC = "conformal_scalar_wormholes"
CHARTS = ("isotropic", "three_parameter")
NAMES = {"isotropic": "Isotropic, $\\chi = \\pi/3$", "three_parameter": "All Three Parameters"}
COORDS = {"isotropic": ["t", "\\rho", "\\theta", "\\phi"], "three_parameter": ["t", "r", "\\theta", "\\phi"]}
GAMMA = "\\gamma = \\dfrac{1 - 12w^2}{1 + 12w^2}"
SIGMA = "\\sigma = \\dfrac{2w}{1 + 12w^2}"
F_DEF = "f = 1 - \\dfrac{4\\mu}{r}"
POWER = "F = f^{\\gamma}"
OMEGA = "\\Omega = \\dfrac{(1 + b)f^{\\sigma} + (1 - b)f^{-\\sigma}}{2}"
PSI = "\\Psi = \\dfrac{(1 + b)f^{\\sigma} - (1 - b)f^{-\\sigma}}{2}"
PARAMETERS = {"isotropic": ["\\mu", "b"],
              "three_parameter": ["\\mu", "b", "w", GAMMA, SIGMA, F_DEF, POWER, OMEGA, PSI]}
SLOPE = "r\\left(r - 4\\mu\\right)"
RATE_TEXTS = {"three_parameter": {"F": {"r": "\\dfrac{4\\mu\\gamma F}{" + SLOPE + "}"},
                                  "Omega": {"r": "\\dfrac{4\\mu\\sigma\\Psi}{" + SLOPE + "}"},
                                  "Psi": {"r": "\\dfrac{4\\mu\\sigma\\Omega}{" + SLOPE + "}"}}}
SPHERE = "\\left(d\\theta^2 + \\sin^2\\theta\\,d\\phi^2\\right)"
ANGLES = ["\\theta \\in [0, \\pi]", "\\phi \\in [0, 2\\pi)"]
DOMAINS = {
    "isotropic": ["t \\in (-\\infty, \\infty)", "\\rho \\in (0, \\infty)"] + ANGLES
    + ["b < 0", "\\rho = \\mu\\sqrt{-b} \\;\\text{(the throat)}",
       "\\rho \\to 0 \\;\\text{(the second flat region)}"],
    "three_parameter": ["t \\in (-\\infty, \\infty)", "r \\in (4\\mu, \\infty)"] + ANGLES
    + ["w \\ge 0,\\quad b < 1",
       "r = 4\\mu \\;\\text{(the horizon, for } w = 0\\text{)}",
       "r = 4\\mu \\;\\text{(a curvature singularity, for } w \\neq 0, 1/6\\text{)}"],
}
DIMENSIONS = {
    (METRIC, "isotropic"): {"t": "T", "\\rho": "L", "\\theta": "1", "\\phi": "1", "\\mu": "L", "b": "1"},
    (METRIC, "three_parameter"): {"t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "\\mu": "L", "b": "1", "w": "1",
                                  "\\gamma": "1", "\\sigma": "1", "f": "1", "F": "1", "\\Omega": "1", "\\Psi": "1"},
}
HELD = {(METRIC, "three_parameter"): ("F", "Omega", "Psi")}
RATES = {(METRIC, s): rates for s, rates in RATE_TEXTS.items()}
DIGITS = 40
TOLERANCE = sp.Float("1e-28")


def register(module):
    module.CHARTS[METRIC] = [lambda system=system: chart_spec(system) for system in CHARTS]


def line(system, c2):
    """The line element of a chart, with c2 standing before dt^2."""
    if system == "three_parameter":
        return ("ds^2 = \\Omega^2\\left(-F\\," + c2 + "dt^2 + \\dfrac{dr^2}{F} + \\dfrac{f}{F}\\,r^2" + SPHERE
                + "\\right)")
    return ("ds^2 = \\left(1 - \\dfrac{b\\mu}{\\rho}\\right)^2\\left(-\\left(1 + \\dfrac{\\mu}{\\rho}\\right)^{-2}" + c2
            + "dt^2 + \\left(1 + \\dfrac{\\mu}{\\rho}\\right)^2\\left(d\\rho^2 + \\rho^2" + SPHERE + "\\right)\\right)")


def reader(system):
    import verify_metrics as vm
    key = (METRIC, system)
    return vm.Reader(COORDS[system], PARAMETERS[system], (), held=HELD.get(key, ()), rates=RATES.get(key))


def chart_spec(system):
    """One chart: its line element with c, the same in x^0 = ct, and the check at points."""
    import verify_metrics as vm
    probe = reader(system)
    if system == "three_parameter":
        P = probe.parameters
        r = probe.symbol["r"]
        gamma, sigma = sp.Symbol("gamma", real=True), sp.Symbol("sigma", real=True)
        pretty = three_parameter_pretty(P["w"], gamma, sigma)
        return {
            "metric_id": METRIC,
            "system": {"id": system, "name": NAMES[system], "coords": COORDS[system], "domains": DOMAINS[system],
                       "parameters": PARAMETERS[system], "line_element": line(system, "c^2")},
            "chart_line_element": line(system, ""),
            "printer": {"lead": [P["Omega"], P["Psi"], P["F"], r, P["mu"], P["b"], gamma, sigma],
                        "factors": [r, P["Omega"], P["F"]], "flip": False},
            "reduce": lambda value: vm.norm(probe.by_rates(value)),
            "pretty": pretty,
            "check": lambda chart: check_three(chart),
        }
    rho, mu, b = probe.symbol["\\rho"], probe.parameters["mu"], probe.parameters["b"]
    return {
        "metric_id": METRIC,
        "system": {"id": system, "name": NAMES[system], "coords": COORDS[system], "domains": DOMAINS[system],
                   "parameters": PARAMETERS[system], "line_element": line(system, "c^2")},
        "chart_line_element": line(system, ""),
        "printer": {"lead": [rho, mu, b], "factors": [mu, b, rho], "flip": False},
        "check": lambda chart: check_chart(chart, system),
    }


def published_metric(vm, metric_id, system_id):
    """A live page's chart as (reader, matrix of its metric components, its coordinates)."""
    data = json.loads((vm.METRICS_DIR / f"{metric_id}.json").read_text(encoding="utf-8"))
    chart = next(c for c in data["coordinates"] if c["id"] == system_id)
    reader = vm.Reader(chart["coords"], [p["symbol"] for p in chart["parameters"]], ())
    values = {tuple(e["indices"]): e["value"] for e in chart["metric_components"]}
    n = len(chart["coords"])
    matrix = sp.Matrix(n, n, lambda i, j: reader(values.get((chart["coords"][i], chart["coords"][j]), "0")))
    return reader, matrix, chart["coords"]


# Points (mu, b, rho) at which the isotropic chart is checked in forty digits: inside the throat, on the
# side of the second flat region, at the throat's own radius and far out, for three wormholes.
POINTS = [(sp.Integer(1), -sp.Rational(1, 2), sp.Rational(3, 10)),
          (sp.Rational(7, 10), -sp.Integer(2), sp.Rational(17, 10)),
          (sp.Rational(3, 2), -sp.Rational(1, 3), sp.Integer(5))]


def check_chart(chart, system):
    """The isotropic chart at points against Barcelo and Visser's (2.2) and (2.3) with their field
    (3.23), against R = 0, the throat at rho = mu sqrt(-b) with the areal radius mu (1 + sqrt(-b))^2
    there, the far field of Schwarzschild's metric with r_s = 2 mu (1 + b), their (3.26); against the
    published isotropic chart of majumdar_papapetrou at b = 0 with m = mu, and, on rho > mu, against
    Omega^2 times the published isotropic chart of fisher_jnw at gamma = 1/2 and Fisher's b = 4 mu,
    with Omega = ((1 + b) h^{1/2} + (1 - b) h^{-1/2})/2 and h = (rho - mu)/(rho + mu), their (2.17)
    with sin(pi/3)/2 sqrt 3 = 1/4 and h^2 = 1 - 2 eta/r."""
    import verify_metrics as vm
    geo, g, ginv = chart.geo, chart.geo.g, chart.geo.ginv
    X = list(chart.symbols)
    n = len(X)
    P = chart.reader.parameters
    mu, b = P["mu"], P["b"]
    rho = X[1]
    name = f"{METRIC}/{system}"
    gamma = geo.christoffel_ull()

    def hessian(f):
        return sp.Matrix(n, n, lambda a, c: sp.diff(f, X[a], X[c])
                         - sum(gamma[k][a][c] * sp.diff(f, X[k]) for k in range(n)))

    def vanishes(e, what, points=POINTS):
        for mm, bb, rr in points:
            value = sp.N(sp.sympify(e).subs({mu: mm, b: bb, rho: rr, X[2]: sp.Rational(7, 10)}), DIGITS)
            if abs(complex(value)) > TOLERANCE:
                raise AssertionError(f"{name}: {what} fails at mu = {mm}, b = {bb}, rho = {rr}: {value}")

    # Their field (3.23) in units kappa = 1, so phi_c^2 = 6 is the sphere where kappa - phi_c^2/6 vanishes.
    phi = sp.sqrt(6) * (b * rho - mu) / (rho - b * mu)
    d = [sp.diff(phi, x) for x in X]
    grad2 = sum(ginv[a, c] * d[a] * d[c] for a in range(n) for c in range(n))
    h = hessian(phi)
    ricci = sp.Matrix(n, n, lambda a, c: vm._at(geo.ricci_ll(), (a, c)))
    for a in range(n):
        for c in range(a, n):
            source = sp.Rational(2, 3) * d[a] * d[c] - g[a, c] * grad2 / 6 - phi * h[a, c] / 3
            vanishes((1 - phi ** 2 / 6) * ricci[a, c] - source,
                     f"Barcelo and Visser's (2.2) in slot {chart.coords_tex[a]}{chart.coords_tex[c]}")
    vanishes(sum(ginv[a, c] * h[a, c] for a in range(n) for c in range(n)), "the field equation (2.3)")
    vanishes(geo.ricci_scalar(), "R = 0")

    # The areal radius has its least value at the throat.
    areal = sp.sqrt(g[2, 2])
    throats = [(mm, bb, mm * sp.sqrt(-bb)) for mm, bb, _ in POINTS]
    vanishes(sp.diff(areal, rho), "the throat at rho = mu sqrt(-b)", throats)
    for mm, bb, rr in throats:
        if not (sp.diff(areal, rho, 2).subs({mu: mm, b: bb, rho: rr}) > 0):
            raise AssertionError(f"{name}: the areal radius is not least at the throat for b = {bb}")
        miss = sp.N((areal - mu * (1 + sp.sqrt(-b)) ** 2).subs({mu: mm, b: bb, rho: rr}), DIGITS)
        if abs(miss) > TOLERANCE:
            raise AssertionError(f"{name}: the areal radius at the throat is not mu (1 + sqrt(-b))^2 for b = {bb}")
    far = sp.limit((g[0, 0] + 1) * rho, rho, sp.oo)
    if sp.simplify(far - 2 * mu * (1 + b)) != 0:
        raise AssertionError(f"{name}: far out g_tt is not -1 + 2 mu (1 + b)/rho, but has {far}")

    # At b = 0, the published single extremal hole of Majumdar and Papapetrou with m = mu.
    reader, hole, coords = published_metric(vm, "majumdar_papapetrou", "isotropic")
    swap = dict(zip((reader.symbol[x] for x in coords), X))
    swap[reader.parameters["m"]] = mu
    extremal = hole.subs(swap)
    for a in range(n):
        for c in range(n):
            if vm.norm(extremal[a, c] - g[a, c].subs(b, 0)) != 0:
                raise AssertionError(f"{name}: at b = 0 the chart misses the published majumdar_papapetrou/"
                                     f"isotropic in slot {chart.coords_tex[a]}{chart.coords_tex[c]}")

    # On rho > mu, Omega^2 times the published Fisher-Janis-Newman-Winicour chart at gamma = 1/2.
    reader, seed, coords = published_metric(vm, "fisher_jnw", "isotropic")
    swap = dict(zip((reader.symbol[x] for x in coords), X))
    swap.update({reader.parameters["b"]: 4 * mu, reader.parameters["gamma"]: sp.Rational(1, 2)})
    seed = seed.subs(swap, simultaneous=True)
    hh = (rho - mu) / (rho + mu)
    omega = ((1 + b) * sp.sqrt(hh) + (1 - b) / sp.sqrt(hh)) / 2
    outside = [(mm, bb, mm + rr) for mm, bb, rr in POINTS]
    for a in range(n):
        for c in range(n):
            vanishes(omega ** 2 * seed[a, c] - g[a, c],
                     f"Omega^2 times fisher_jnw/isotropic in slot {chart.coords_tex[a]}{chart.coords_tex[c]}", outside)


# Points (mu, b, w, r) of the class at which the chart in r is checked: a wormhole, w = 1/6, and members
# with a naked singularity on each side of it, one with gamma < 0, each outside r = 4 mu.
CLASS_POINTS = [(sp.Integer(1), -sp.Rational(1, 2), sp.Rational(1, 6), sp.Integer(5)),
                (sp.Rational(7, 10), sp.Rational(3, 10), sp.Rational(1, 10), sp.Integer(4)),
                (sp.Rational(3, 2), -sp.Integer(2), sp.Integer(1), sp.Integer(9)),
                (sp.Rational(1, 2), sp.Rational(1, 2), sp.Rational(1, 4), sp.Rational(21, 10))]


def three_parameter_pretty(w, gamma, sigma):
    """A `pretty` for the chart in r: each value comes from the canonical form as a rational function of w,
    which is written in gamma and sigma by w = sigma/(1 + gamma), its numerator and denominator reduced by
    12 sigma^2 = 1 - gamma^2 to the first power of sigma, and the quotient factored."""
    conic = 12 * sigma ** 2 + gamma ** 2 - 1

    def reduce(p):
        return sp.expand(sp.rem(sp.expand(p), conic, sigma)) if p.has(sigma) else sp.expand(p)

    def pretty(value):
        x = sp.sympify(value)
        if not x.has(w):
            return sp.factor(x)
        num, den = sp.fraction(sp.together(x.subs(w, sigma / (1 + gamma))))
        num, den = reduce(num), reduce(den)
        # A denominator left with sigma is cleared by its conjugate, as (d0 + d1 sigma)(d0 - d1 sigma) is
        # d0^2 - d1^2 sigma^2, which the conic writes without sigma.
        d1 = sp.Poly(den, sigma).coeff_monomial(sigma)
        if d1 != 0:
            d0 = sp.Poly(den, sigma).coeff_monomial(1)
            num, den = reduce(num * (d0 - d1 * sigma)), reduce(d0 ** 2 - d1 ** 2 * sigma ** 2)
        out = sp.factor(sp.cancel(num / den))
        if out.has(w):
            raise ValueError(f"{METRIC}: {value} is not written in gamma and sigma")
        return out
    return pretty


def check_three(chart):
    """The chart in r at points (mu, b, w, r): Barcelo and Visser's (2.2) and (2.3) with their field (2.14),
    phi_c = sqrt(6 kappa) Psi/Omega; R = 0; the metric as Omega^2 times the published spherical chart
    of fisher_jnw with Fisher's b = 4 mu and the same gamma; the far field of Schwarzschild's metric
    with r_s = 4 mu (gamma + 2 b sigma), twice their G_eff M = eta (cos chi + tan(Delta/2) sin chi/sqrt 3)
    of footnote 5; and at gamma = 1/2 the isotropic chart of this file along r = rho (1 + mu/rho)^2,
    their (3.21), on rho > mu. Every held name is written by its declared slope and then out."""
    import verify_metrics as vm
    probe = reader("three_parameter")
    geo, g, ginv = chart.geo, chart.geo.g, chart.geo.ginv
    X = list(chart.symbols)
    n = len(X)
    P = probe.parameters
    mu, b, w, r = P["mu"], P["b"], P["w"], X[1]
    name = f"{METRIC}/three_parameter"
    gamma = geo.christoffel_ull()

    def at(e, point):
        mm, bb, gg, rr = point
        e = probe.by_rates(sp.sympify(e)).subs(probe.held).doit()
        return sp.N(e.subs({mu: mm, b: bb, w: gg}).subs({r: rr, X[2]: sp.Rational(7, 10)}), DIGITS)

    def vanishes(e, what, points=CLASS_POINTS):
        for point in points:
            value = at(e, point)
            if abs(complex(value)) > TOLERANCE:
                raise AssertionError(f"{name}: {what} fails at (mu, b, w, r) = {point}: {value}")

    def hessian(f):
        return sp.Matrix(n, n, lambda a, c: sp.diff(f, X[a], X[c])
                         - sum(gamma[k][a][c] * sp.diff(f, X[k]) for k in range(n)))

    phi = sp.sqrt(6) * P["Psi"] / P["Omega"]
    d = [sp.diff(phi, x) for x in X]
    grad2 = sum(ginv[a, c] * d[a] * d[c] for a in range(n) for c in range(n))
    h = hessian(phi)
    ricci = sp.Matrix(n, n, lambda a, c: vm._at(geo.ricci_ll(), (a, c)))
    for a in range(n):
        for c in range(a, n):
            source = sp.Rational(2, 3) * d[a] * d[c] - g[a, c] * grad2 / 6 - phi * h[a, c] / 3
            vanishes((1 - phi ** 2 / 6) * ricci[a, c] - source,
                     f"Barcelo and Visser's (2.2) in slot {chart.coords_tex[a]}{chart.coords_tex[c]}")
    vanishes(sum(ginv[a, c] * h[a, c] for a in range(n) for c in range(n)), "the field equation (2.3)")
    vanishes(geo.ricci_scalar(), "R = 0")

    # Omega^2 times the published spherical chart of Fisher, Janis, Newman and Winicour.
    seed_reader, seed, coords = published_metric(vm, "fisher_jnw", "spherical")
    swap = dict(zip((seed_reader.symbol[x] for x in coords), X))
    swap.update({seed_reader.parameters["b"]: 4 * mu, seed_reader.parameters["gamma"]: P["gamma"]})
    seed = seed.subs(swap, simultaneous=True)
    for a in range(n):
        for c in range(n):
            vanishes(P["Omega"] ** 2 * seed[a, c] - g[a, c],
                     f"Omega^2 times fisher_jnw/spherical in slot {chart.coords_tex[a]}{chart.coords_tex[c]}")

    # Far out, g_tt = -1 + r_s/r + O(1/r^2) with r_s = 4 mu (gamma + 2 b sigma).
    for mm, bb, gg, _ in CLASS_POINTS:
        far = sp.Integer(10) ** 20
        value = at((g[0, 0] + 1) * r, (mm, bb, gg, far))
        want = sp.N(4 * mm * ((1 - 12 * gg ** 2) + 4 * bb * gg) / (1 + 12 * gg ** 2), DIGITS)
        if abs(value - want) > sp.Float("1e-15"):
            raise AssertionError(f"{name}: far out r (g_tt + 1) is {value}, not 4 mu (gamma + 2 b sigma) = {want}")

    # At gamma = 1/2, the isotropic chart along r = rho (1 + mu/rho)^2, on rho > mu.
    iso = reader("isotropic")
    iso_g = vm.metric_from_line_element(iso, line("isotropic", ""), COORDS["isotropic"])
    rho = iso.symbol["\\rho"]
    swap = {iso.parameters["mu"]: mu, iso.parameters["b"]: b, iso.symbol["t"]: X[0], iso.symbol["\\theta"]: X[2],
            iso.symbol["\\phi"]: X[3]}
    iso_g = iso_g.subs(swap, simultaneous=True)
    radius = rho * (1 + mu / rho) ** 2
    J = sp.diag(1, sp.diff(radius, rho), 1, 1)
    for mm, bb, _, rr in CLASS_POINTS:
        point = {mu: mm, b: bb, w: sp.Rational(1, 6), rho: mm + rr / 3, X[2]: sp.Rational(7, 10)}
        mine = probe.by_rates(g).subs(probe.held).doit()
        pulled = J.T * mine.subs(r, radius) * J
        for a in range(n):
            for c in range(n):
                gap = sp.N((pulled[a, c] - iso_g[a, c]).subs(point), DIGITS)
                if abs(complex(gap)) > TOLERANCE:
                    raise AssertionError(f"{name}: at gamma = 1/2 the chart misses the isotropic chart along "
                                         f"r = rho (1 + mu/rho)^2 in slot {chart.coords_tex[a]}{chart.coords_tex[c]}")
