"""The charts of Schein and Aichelburg's wormhole between two charged shells, gr-qc/9606069,
Phys. Rev. Lett. 77, 4130 (1996), and Aichelburg and Schein, Acta Phys. Pol. B 29, 1025 (1998).
charged_shell_wormhole.md beside this file records each chart's source.

Outside the two spheres of radius R centred at z = +-d_1 the metric is Majumdar and Papapetrou's,
their (1),

    ds^2 = -c^2 dT^2/V^2 + V^2 (dx^2 + dy^2 + dz^2),

with V harmonic, equal to V_0 = 1 + m_1/R on both spheres and to 1 at infinity: the image series
of their (6), V = 1 + sum_n m_n (1/|x - d_n z| + 1/|x + d_n z|), d_n = a coth(n mu_0) and
m_n = (-1)^{n+1} m_1 sinh(mu_0)/sinh(n mu_0), where R = a/sinh(mu_0) and d_1 = a coth(mu_0). The
founders write the focal length a as c; it is a here, since c is the speed of light. In their
bispherical coordinates, coth(mu) = (x^2 + y^2 + z^2 + a^2)/(2az),
cot(eta) = (x^2 + y^2 + z^2 - a^2)/(2a sqrt(x^2 + y^2)), the metric is their (bisph) and V the series
after it; in polar coordinates about z = -d_1 it is their (sphcoor). Between the shells the metric is
Reissner and Nordstrom's, their (RN) with f = 1 - 2m/r + e^2/r^2, the charge e written q here, since
e is Euler's number on the page; each shell sits at the areal radius R V_0 = R + m_1.
"""
import json

import sympy as sp

METRIC = "charged_shell_wormhole"
CHARTS = ("cartesian", "bispherical", "polar", "interior")

DIMENSIONS = {
    (METRIC, "cartesian"): {"T": "T", "x": "L", "y": "L", "z": "L", "V": "1", "R": "L", "d_1": "L", "m_1": "L"},
    (METRIC, "bispherical"): {"T": "T", "\\mu": "1", "\\eta": "1", "\\phi": "1", "V": "1", "W": "1", "a": "L",
                              "\\mu_0": "1", "m_1": "L"},
    (METRIC, "polar"): {"T": "T", "\\rho": "L", "\\theta": "1", "\\phi": "1", "V": "1", "R": "L", "d_1": "L",
                        "m_1": "L"},
    (METRIC, "interior"): {"t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "m": "L", "q": "L", "R": "L",
                           "m_1": "L"},
}

NAMES = {"cartesian": "Cartesian Exterior", "bispherical": "Bispherical Exterior", "polar": "Polar Exterior",
         "interior": "Reissner-Nordström Interior"}
COORDS = {"cartesian": ["T", "x", "y", "z"], "bispherical": ["T", "\\mu", "\\eta", "\\phi"],
          "polar": ["T", "\\rho", "\\theta", "\\phi"], "interior": ["t", "r", "\\theta", "\\phi"]}
PARAMETERS = {"cartesian": ["V = V(x,y,z)", "R", "d_1", "m_1"],
              "bispherical": ["V = V(\\mu,\\eta)", "W = \\cosh\\mu - \\cos\\eta", "a", "\\mu_0", "m_1"],
              "polar": ["V = V(\\rho,\\theta)", "R", "d_1", "m_1"],
              "interior": ["m", "q", "R", "m_1"]}
SPHERE = "\\left(d\\theta^2 + \\sin^2\\theta\\,d\\phi^2\\right)"
ANGLES = ["\\theta \\in [0, \\pi]", "\\phi \\in [0, 2\\pi)"]
DOMAINS = {
    "cartesian": ["T \\in (-\\infty, \\infty)", "x \\in (-\\infty, \\infty)", "y \\in (-\\infty, \\infty)",
                  "z \\in (-\\infty, \\infty)",
                  "x^2 + y^2 + (z \\mp d_1)^2 \\geq R^2 \\;\\text{(outside both spheres)}",
                  "x^2 + y^2 + (z \\mp d_1)^2 = R^2 \\;\\text{(the shells)}"],
    "bispherical": ["T \\in (-\\infty, \\infty)", "\\mu \\in [-\\mu_0, \\mu_0]", "\\eta \\in [0, \\pi]",
                    "\\phi \\in [0, 2\\pi)", "\\mu = \\pm\\mu_0 \\;\\text{(the shells)}",
                    "(\\mu, \\eta) = (0, 0) \\;\\text{(spatial infinity)}"],
    "polar": ["T \\in (-\\infty, \\infty)", "\\rho \\in [R, \\infty)"] + ANGLES
             + ["\\rho^2 - 4d_1\\rho\\cos\\theta + 4d_1^2 \\geq R^2 \\;\\text{(outside the other sphere)}",
                "\\rho = R \\;\\text{(one shell)}"],
    "interior": ["t \\in (-\\infty, \\infty)", "r \\in (0, R + m_1]"] + ANGLES
                + ["r = R + m_1 \\;\\text{(the shells)}",
                   "r = m \\pm \\sqrt{m^2 - q^2} \\;\\text{(the horizons, for}\\; |q| < m\\text{)}"],
}
F = "\\left(1 - \\dfrac{2m}{r} + \\dfrac{q^2}{r^2}\\right)"
F_BARE = "1 - \\dfrac{2m}{r} + \\dfrac{q^2}{r^2}"
# The bispherical chart holds W = cosh(mu) - cos(eta), the factor of the founders' (bisph), as a name,
# and bispherical_pretty writes each value back in it.
HELD = {(METRIC, "bispherical"): ("W",)}
# W has no declared rates, so every comparison, the checker's among them, writes it out.
RATES = {}
KRETSCHMANN = "\\dfrac{8\\left(6m^2r^2 - 12mq^2r + 7q^4\\right)}{r^8}"


def line(system, c2):
    """The line element of a chart, c2 standing before the square of the time."""
    if system == "cartesian":
        return "ds^2 = -\\dfrac{" + c2 + "dT^2}{V^2} + V^2\\left(dx^2 + dy^2 + dz^2\\right)"
    if system == "bispherical":
        return ("ds^2 = -\\dfrac{" + c2 + "dT^2}{V^2} + \\dfrac{V^2a^2}{W^2}"
                "\\left(d\\mu^2 + d\\eta^2 + \\sin^2\\eta\\,d\\phi^2\\right)")
    if system == "polar":
        return "ds^2 = -\\dfrac{" + c2 + "dT^2}{V^2} + V^2\\left(d\\rho^2 + \\rho^2" + SPHERE + "\\right)"
    return "ds^2 = -" + F + c2 + "dt^2 + \\dfrac{dr^2}{" + F_BARE + "} + r^2" + SPHERE


def chart_spec(system):
    """One chart: its line element with c, the same in x^0 = cT or ct, and the check at points."""
    import verify_metrics as vm
    import chart_printer as cp
    coords, parameters = COORDS[system], PARAMETERS[system]
    key = (METRIC, system)
    probe = vm.Reader(coords, parameters, (), held=HELD.get(key, ()), rates=RATES.get(key))
    spec = {
        "metric_id": METRIC,
        "system": {"id": system, "name": NAMES[system], "coords": coords, "domains": DOMAINS[system],
                   "parameters": parameters, "line_element": line(system, "c^2")},
        "chart_line_element": line(system, ""),
        "check": lambda chart: check_chart(chart, system),
    }
    if system == "interior":
        r, m, q = probe.symbol["r"], probe.parameters["m"], probe.parameters["q"]
        spec["printer"] = {"lead": [r, m, q], "factors": [m, q, r], "flip": False}
        spec["components"] = {"metric_components": {("t", "t"): "-" + F, ("r", "r"): F + "^{-1}"},
                              "inverse_metric_components": {("t", "t"): "-" + F + "^{-1}", ("r", "r"): F_BARE}}
        spec["kretschmann"] = KRETSCHMANN
    elif system == "bispherical":
        V, W = probe.parameters["V"], probe.parameters["W"]
        spec["printer"] = {"lead": [V, W], "flip": False}
        mu, eta = probe.symbol["\\mu"], probe.symbol["\\eta"]
        # Every value is reduced with W written out, and the printed form writes it back.
        spec["reduce"] = lambda value: vm.norm(probe.surface(sp.sympify(value)))
        spec["pretty"] = bispherical_pretty(cp, probe.symbol["\\mu"], probe.symbol["\\eta"], W)
    else:
        V = probe.parameters["V"]
        spec["printer"] = {"lead": [V], "flip": False,
                           "collect": lambda poly, printer: cp.collect_by(poly, [V], printer)}
    return spec


def bispherical_pretty(cp, mu, eta, W):
    """A `pretty` for the bispherical chart: the printer's hyperbolic rewriting of the exponentials
    the checker hands back, then cosh(mu) written W + cos(eta) and sinh(mu)^2 reduced by
    sinh(mu)^2 = (W + cos(eta))^2 - 1, so that the curvature of the coordinates alone cancels, and
    the value factored."""
    hyperbolic = cp.hyperbolic(mu)
    S, ce, sn = sp.symbols("_S _ce _sn")
    relation = [S ** 2 - (W + ce) ** 2 + 1, sn ** 2 + ce ** 2 - 1]

    def reduce(p):
        return sp.expand(sp.reduced(sp.expand(p), relation, S, sn, ce, W)[1])

    def pretty(value):
        x = sp.sympify(hyperbolic(value))
        x = x.xreplace({sp.cosh(mu): W + ce, sp.sinh(mu): S, sp.cos(eta): ce, sp.sin(eta): sn,
                        sp.csc(eta): 1 / sn, sp.cot(eta): ce / sn})
        num, den = (reduce(p) for p in sp.fraction(sp.together(x)))
        return sp.factor(sp.cancel(num / den)).subs({S: sp.sinh(mu), ce: sp.cos(eta), sn: sp.sin(eta)})
    return pretty


def published_metric(vm, metric_id, system_id):
    """A live page's chart as (reader, matrix of its metric components, its coordinates)."""
    data = json.loads((vm.METRICS_DIR / f"{metric_id}.json").read_text(encoding="utf-8"))
    chart = next(c for c in data["coordinates"] if c["id"] == system_id)
    reader = vm.Reader(chart["coords"], [p["symbol"] for p in chart["parameters"]], ())
    values = {tuple(e["indices"]): e["value"] for e in chart["metric_components"]}
    n = len(chart["coords"])
    matrix = sp.Matrix(n, n, lambda i, j: reader(values.get((chart["coords"][i], chart["coords"][j]), "0")))
    return reader, matrix, chart["coords"]


def images(a, mu0, m1, terms):
    """The image charges and their heights, (m_n, d_n) for n = 1 .. terms, the founders' (4) and (5)
    solved by Lindquist's parameters."""
    return [((-1) ** (n + 1) * m1 * sp.sinh(mu0) / sp.sinh(n * mu0), a * sp.coth(n * mu0)) for n in range(1, terms + 1)]


def cartesian_potential(x, y, z, a, mu0, m1, terms):
    """The founders' (6), V = 1 + sum_n m_n (1/|x - d_n z| + 1/|x + d_n z|), to the given number of images."""
    return 1 + sum(mn * (1 / sp.sqrt(x ** 2 + y ** 2 + (z - dn) ** 2) + 1 / sp.sqrt(x ** 2 + y ** 2 + (z + dn) ** 2))
                   for mn, dn in images(a, mu0, m1, terms))


def bispherical_potential(mu, eta, mu0, m1, R, terms):
    """The founders' V in bispherical coordinates, the series after their (bisph), n from -terms to terms."""
    root = sp.sqrt(sp.cosh(mu) - sp.cos(eta))
    return 1 + m1 / R * (1 + root * sum((-1) ** (n + 1) / sp.sqrt(sp.cosh(mu + 2 * n * mu0) - sp.cos(eta))
                                       for n in range(-terms, terms + 1)))


def to_cartesian(mu, eta, phi, a):
    """The point (x, y, z) with bispherical coordinates (mu, eta, phi) and focal length a."""
    den = sp.cosh(mu) - sp.cos(eta)
    return a * sp.sin(eta) * sp.cos(phi) / den, a * sp.sin(eta) * sp.sin(phi) / den, a * sp.sinh(mu) / den


# Points of the parameters and of the coordinates at which the charts are checked in forty digits.
SHAPE = {"a": sp.Rational(3, 2), "mu0": sp.Rational(6, 5), "m1": sp.Rational(1, 5)}
INTERIOR = [(sp.Rational(1, 3), sp.Rational(1, 4), sp.Rational(17, 10)),
            (sp.Rational(1, 2), sp.Rational(1, 2), sp.Rational(9, 4)),
            (sp.Rational(2, 5), sp.Rational(1, 10), sp.Rational(3, 10))]


def check_chart(chart, system):
    """Each chart at points.

    The Cartesian chart is the published Cartesian chart of majumdar_papapetrou with U written V and
    t written T, slot by slot, and with V the image series truncated, which is harmonic, its Einstein
    tensor equals 2(F_ma F_n^a - g_mn F^2/4) for A = c dT/V in units G = 4 pi epsilon_0 = 1, at points
    outside the spheres. The bispherical and polar charts are the Cartesian chart pulled back, at
    points, with V carried along, and the bispherical series is the Cartesian series at those points.
    The interior is the published spherical chart of rn_metric at r_s = 2m and r_q = q, slot by slot,
    with the Einstein tensor of a point charge's field, G^a_b = (q^2/r^4) diag(-1, -1, 1, 1), at points."""
    import verify_metrics as vm
    geo, g = chart.geo, chart.geo.g
    X = list(chart.symbols)
    P = chart.reader.parameters
    name = f"{METRIC}/{system}"

    def small(value, what):
        if abs(complex(sp.N(value, 40))) > 1e-25:
            raise AssertionError(f"{name}: {what}")

    if system == "interior":
        m, q, r = P["m"], P["q"], X[1]
        reader, hole, coords = published_metric(vm, "rn_metric", "spherical")
        swap = dict(zip((reader.symbol[x] for x in coords), X))
        swap.update({reader.parameters["r_s"]: 2 * m, reader.parameters["r_q"]: q})
        rn = hole.subs(swap)
        for i in range(4):
            for j in range(4):
                if vm.norm(rn[i, j] - g[i, j]) != 0:
                    raise AssertionError(f"{name}: the chart misses the published rn_metric/spherical at "
                                         f"r_s = 2m, r_q = q in slot {chart.coords_tex[i]}{chart.coords_tex[j]}")
        mixed = geo.ginv * sp.Matrix(geo.einstein_ll())
        for mm, qq, rr in INTERIOR:
            at = {m: mm, q: qq, r: rr, X[2]: sp.Rational(7, 10)}
            for i in range(4):
                for j in range(4):
                    want = 0 if i != j else (-1 if i < 2 else 1) * qq ** 2 / rr ** 4
                    small(mixed[i, j].subs(at) - want,
                          f"G^a_b is not a point charge's field in slot {chart.coords_tex[i]}{chart.coords_tex[j]}")
        return

    a, mu0, m1 = SHAPE["a"], SHAPE["mu0"], SHAPE["m1"]
    R, d1 = a / sp.sinh(mu0), a * sp.coth(mu0)
    V = P["V"]
    if system == "cartesian":
        reader, mp, coords = published_metric(vm, "majumdar_papapetrou", "cartesian")
        swap = dict(zip((reader.symbol[x] for x in coords), X))
        swap[reader.parameters["U"]] = V
        mp = mp.subs(swap)
        for i in range(4):
            for j in range(4):
                if vm.norm(mp[i, j] - g[i, j]) != 0:
                    raise AssertionError(f"{name}: the chart misses the published majumdar_papapetrou/cartesian "
                                         f"in slot {chart.coords_tex[i]}{chart.coords_tex[j]}")
        # The truncated image series is harmonic; V and its derivatives are evaluated at each point
        # and put into the Einstein tensor the chart printed, which holds V free.
        concrete = cartesian_potential(X[1], X[2], X[3], a, mu0, m1, 3)
        einstein = sp.Matrix(geo.einstein_ll())
        derivatives = sorted(einstein.atoms(sp.Derivative), key=lambda d: -len(d.variables))
        for point in [(sp.Rational(3, 10), sp.Rational(-1, 5), sp.Rational(1, 7)),
                      (sp.Rational(2, 1), sp.Rational(1, 2), sp.Rational(5, 2))]:
            at = dict(zip(X[1:], point))
            values = {d: sp.N(sp.diff(concrete, *d.variables).subs(at), 50) for d in derivatives}
            v = sp.N(concrete.subs(at), 50)
            grad = [sp.N(sp.diff(concrete, x).subs(at), 50) for x in X[1:]]
            gp = sp.diag(-1 / v ** 2, v ** 2, v ** 2, v ** 2)
            gip = gp.inv()
            # A = c dT/V, so F_i0 = d_i(1/V) = -d_iV/V^2 for each spatial i.
            Fp = sp.zeros(4, 4)
            for k in range(3):
                Fp[k + 1, 0] = -grad[k] / v ** 2
                Fp[0, k + 1] = grad[k] / v ** 2
            Fu = gip * Fp * gip
            F2 = sum(Fp[i, j] * Fu[i, j] for i in range(4) for j in range(4))
            Ep = einstein.subs(values).subs(V, v).subs(at)
            for i in range(4):
                for j in range(i, 4):
                    stress = 2 * (sum(Fp[i, k] * Fp[j, l] * gip[k, l] for k in range(4) for l in range(4))
                                  - gp[i, j] * F2 / 4)
                    small(Ep[i, j] - stress,
                          f"G is not the Maxwell stress of A = c dT/V in slot {chart.coords_tex[i]}{chart.coords_tex[j]}")
        return

    # The bispherical and polar charts: the Cartesian metric pulled back at points.
    if system == "bispherical":
        g = g.subs(P["W"], sp.cosh(X[1]) - sp.cos(X[2]))
        here = {P["a"]: a, P["mu_0"]: mu0, P["m_1"]: m1}
        image = [X[0], *to_cartesian(X[1], X[2], X[3], a)]
        points = [(sp.Rational(1, 3), sp.Rational(7, 5), sp.Rational(1, 2)),
                  (sp.Rational(-9, 10), sp.Rational(2, 5), sp.Rational(3, 1))]
    else:
        here = {P["R"]: R, P["d_1"]: d1, P["m_1"]: m1}
        image = [X[0], X[1] * sp.sin(X[2]) * sp.cos(X[3]), X[1] * sp.sin(X[2]) * sp.sin(X[3]),
                 -d1 + X[1] * sp.cos(X[2])]
        points = [(sp.Rational(13, 10) * R, sp.Rational(9, 10), sp.Rational(1, 3)),
                  (sp.Rational(5, 2), sp.Rational(5, 2), sp.Rational(-2, 3))]
    J = sp.Matrix(4, 4, lambda i, j: sp.diff(image[i], X[j]))
    flat = sp.diag(-1, 1, 1, 1)
    for point in points:
        at = dict(zip(X[1:], point))
        Jp = J.subs(at)
        pulled = Jp.T * flat * Jp
        for i in range(4):
            for j in range(4):
                # The chart's metric is -c^2dT^2/V^2 and V^2 times the flat spatial metric: compare it
                # with V = 1, the flat metric pulled back, in each slot.
                small(g[i, j].subs(V, 1).subs(here).subs(at) - pulled[i, j],
                      f"the flat metric pulled back misses slot {chart.coords_tex[i]}{chart.coords_tex[j]}")
        if system == "bispherical":
            xyz = [image[k].subs(at) for k in (1, 2, 3)]
            cart = cartesian_potential(*xyz, a, mu0, m1, 60)
            bisp = bispherical_potential(at[X[1]], at[X[2]], mu0, m1, R, 60)
            small(cart - bisp, f"the bispherical series misses the Cartesian series at {point}")
        for i in range(4):
            for j in range(4):
                if i != j and g[i, j] != 0:
                    raise AssertionError(f"{name}: slot {chart.coords_tex[i]}{chart.coords_tex[j]} is not zero")
    if vm.norm(g[0, 0] + 1 / V ** 2) != 0:
        raise AssertionError(f"{name}: g_00 is not -1/V^2")
