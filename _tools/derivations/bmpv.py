"""The charts of the spinning black hole of Breckenridge, Myers, Peet and Vafa, Phys. Lett. B 391, 93
(1997), hep-th/9602065. bmpv.md beside this file records each chart's source.

Every chart is in the signature (-,+,+,+,+) and holds the two parameters of Gauntlett, Myers and
Townsend, CQG 16, 1 (1999), hep-th/9810204: the mass parameter mu, an area, and the spin parameter j, a
length cubed. Breckenridge, Myers, Peet and Vafa and Gibbons and Herdeiro write omega = j/mu, a length:

    bmpv       Breckenridge, Myers, Peet and Vafa's extremal metric of their section 3, after their shift
               r^2 -> r^2 + mu, with the Hopf angles of the round 3-sphere,
               d theta^2 + sin^2 theta d phi^2 + cos^2 theta d psi^2;
    euler      Gibbons and Herdeiro's eq. (2.1), CQG 16, 3619 (1999), hep-th/9906098, the same radius with
               the Euler angles (alpha, beta, gamma) and f = mu omega/(2(r^2 - mu));
    isotropic  Gauntlett, Myers and Townsend's eq. (fivemet), the harmonic H = 1 + mu/r^2 on flat space of
               four dimensions with its Euler angles (theta, phi, psi);
    ingoing    Gauntlett, Myers and Townsend's eq. (nicemetric), regular through the horizon, written with
               the radius of the first two charts, r = sqrt(mu) lambda, and the advanced time v, their
               sqrt(mu) u, where r^6 > j^2.

The isotropic radius is sqrt(r^2 - mu) in the radius of the other charts, and the Euler angles of the
isotropic and ingoing charts are Gibbons and Herdeiro's (beta, alpha, gamma).
"""
import itertools

import sympy as sp

METRIC = "bmpv"
CHARTS = ("bmpv", "euler", "isotropic", "ingoing")
NAMES = {"bmpv": "Breckenridge-Myers-Peet-Vafa Radius, Hopf Angles",
         "euler": "Breckenridge-Myers-Peet-Vafa Radius, Euler Angles",
         "isotropic": "Isotropic, Euler Angles",
         "ingoing": "Ingoing, Euler Angles"}

DIMENSIONS = {
    (METRIC, "bmpv"): {"t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "\\psi": "1", "\\mu": "L^2", "j": "L^3"},
    (METRIC, "euler"): {"t": "T", "r": "L", "\\alpha": "1", "\\beta": "1", "\\gamma": "1", "\\mu": "L^2", "j": "L^3"},
    (METRIC, "isotropic"): {"t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "\\psi": "1", "\\mu": "L^2", "j": "L^3",
                            "H": "1"},
    (METRIC, "ingoing"): {"v": "T", "r": "L", "\\theta": "1", "\\phi": "1", "\\tilde\\psi": "1", "\\mu": "L^2",
                          "j": "L^3"},
}
HELD = {}
RATES = {}

REALS = " \\in (-\\infty, \\infty)"
TURN = " \\in [0, 2\\pi)"
HORIZON = "r = \\sqrt{\\mu} \\;\\text{(the horizon, for } j^2 < \\mu^3\\text{)}"
VLS = "r^6 = j^2 \\;\\text{(the velocity of light surface)}"
DELTA = "\\left(1 - \\dfrac{\\mu}{r^2}\\right)"
EULER = "d\\alpha^2 + d\\beta^2 + d\\gamma^2 + 2\\cos\\beta\\,d\\alpha\\,d\\gamma"
SIGMA = "\\left(d\\psi + \\cos\\theta\\,d\\phi\\right)"
SIGMA_TILDE = "\\left(d\\tilde\\psi + \\cos\\theta\\,d\\phi\\right)"
H_DEF = "H = 1 + \\dfrac{\\mu}{r^2}"


def _bmpv(c):
    ct = "c\\,dt" if c else "dt"
    return ("ds^2 = -" + DELTA + "^2\\left(" + ct + " - \\dfrac{j\\sin^2\\theta}{r^2 - \\mu}d\\phi"
            " + \\dfrac{j\\cos^2\\theta}{r^2 - \\mu}d\\psi\\right)^2 + " + DELTA + "^{-2}dr^2"
            " + r^2\\left(d\\theta^2 + \\sin^2\\theta\\,d\\phi^2 + \\cos^2\\theta\\,d\\psi^2\\right)")


def _euler(c):
    ct = "c\\,dt" if c else "dt"
    return ("ds^2 = -" + DELTA + "^2\\left(" + ct + " + \\dfrac{j}{2\\left(r^2 - \\mu\\right)}\\left(d\\gamma"
            " + \\cos\\beta\\,d\\alpha\\right)\\right)^2 + " + DELTA + "^{-2}dr^2 + \\dfrac{r^2}{4}\\left(" + EULER
            + "\\right)")


def _isotropic(c):
    ct = "c\\,dt" if c else "dt"
    return ("ds^2 = -\\dfrac{1}{H^2}\\left(" + ct + " + \\dfrac{j}{2r^2}" + SIGMA + "\\right)^2 + H\\left(dr^2"
            " + \\dfrac{r^2}{4}\\left(d\\theta^2 + \\sin^2\\theta\\,d\\phi^2 + " + SIGMA + "^2\\right)\\right)")


def _ingoing(c):
    cv = "c\\,dv" if c else "dv"
    c2 = "c^2" if c else ""
    return ("ds^2 = -" + DELTA + "^2" + c2 + "dv^2 + \\dfrac{2r^3}{\\sqrt{r^6 - j^2}}" + cv + "\\,dr - " + DELTA
            + "\\dfrac{j}{r^2}" + cv + "\\," + SIGMA_TILDE + " + \\dfrac{r^2}{4}\\left(d\\theta^2"
            " + \\sin^2\\theta\\,d\\phi^2\\right) + \\dfrac{r^6 - j^2}{4r^4}" + SIGMA_TILDE + "^2")


SYSTEMS = {
    "bmpv": {"coords": ["t", "r", "\\theta", "\\phi", "\\psi"], "parameters": ["\\mu", "j"],
             "domains": ["t" + REALS, "r \\in (\\sqrt{\\mu}, \\infty)", "\\theta \\in [0, \\pi/2]",
                         "\\phi" + TURN, "\\psi" + TURN, HORIZON, VLS],
             "line": _bmpv},
    "euler": {"coords": ["t", "r", "\\alpha", "\\beta", "\\gamma"], "parameters": ["\\mu", "j"],
              "domains": ["t" + REALS, "r \\in (\\sqrt{\\mu}, \\infty)", "\\alpha" + TURN,
                          "\\beta \\in [0, \\pi]", "\\gamma \\in [0, 4\\pi)", HORIZON, VLS],
              "line": _euler},
    "isotropic": {"coords": ["t", "r", "\\theta", "\\phi", "\\psi"], "parameters": ["\\mu", "j", H_DEF],
                  "domains": ["t" + REALS, "r \\in (0, \\infty)", "\\theta \\in [0, \\pi]", "\\phi" + TURN,
                              "\\psi \\in [0, 4\\pi)",
                              "r = 0 \\;\\text{(the horizon, for } j^2 < \\mu^3\\text{)}"],
                  "line": _isotropic},
    "ingoing": {"coords": ["v", "r", "\\theta", "\\phi", "\\tilde\\psi"], "parameters": ["\\mu", "j"],
                "domains": ["v" + REALS, "r^6 > j^2", "\\theta \\in [0, \\pi]", "\\phi" + TURN,
                            "\\tilde\\psi \\in [0, 4\\pi)", HORIZON],
                "line": _ingoing},
}

PRINT = None
# The sign of sqrt(-g) nabla_a F^ab against (1/2 sqrt 3) eps^{bcdef} F_cd F_ef, with eps^{01234} = 1 in the
# chart's order of coordinates, as each check last found it.
ORIENTATION = {}


def register(module):
    global PRINT
    PRINT = module
    module.CHARTS[METRIC] = [lambda system=system: chart_spec(system) for system in CHARTS]


def chart_spec(system):
    """One chart: its line element, the same in the chart x^0 = ct, the printer's order of names, and
    the check of the chart at points."""
    import verify_metrics as vm
    chart = SYSTEMS[system]
    probe = vm.Reader(chart["coords"], chart["parameters"], ())
    r = probe.symbol["r"]
    mu, j = probe.parameters["mu"], probe.parameters["j"]
    return {
        "metric_id": METRIC,
        "system": {"id": system, "name": NAMES[system], "coords": chart["coords"], "domains": chart["domains"],
                   "parameters": chart["parameters"], "line_element": chart["line"](True)},
        "chart_line_element": chart["line"](False),
        "printer": {"lead": [r, mu, j], "factors": [mu, j, r], "flip": False},
        "check": lambda built: check_chart(built, system),
    }


# Parameters and events of every check: one hole with its closed timelike curves behind the horizon,
# j^2 < mu^3, and one spun past the bound, j^2 > mu^3. Each event lies off the axes of its angles, with
# r^2 > mu in the charts of Breckenridge, Myers, Peet and Vafa's radius and r^6 > j^2 in the ingoing one.
CASES = [(sp.Rational(13, 10), sp.Rational(7, 10)), (sp.Rational(1, 1), sp.Rational(9, 5))]
RADII = [sp.Rational(17, 10), sp.Rational(5, 2), sp.Rational(31, 10)]
ANGLE_POINTS = [(sp.Rational(3, 5), sp.Rational(7, 5), sp.Rational(2, 5)),
                (sp.Rational(11, 10), sp.Rational(1, 3), sp.Rational(5, 2)),
                (sp.Rational(2, 7), sp.Rational(9, 4), sp.Rational(6, 5))]


def potential(system, s, mu, j):
    """Gauntlett, Myers and Townsend's A = (sqrt 3/2)(e^0 - dt), e^0 = H^{-1}(dt + a), in the chart's own
    coordinates, x^0 = ct; in Breckenridge, Myers, Peet and Vafa's radius H^{-1} = 1 - mu/r^2 and
    a = j sigma_3/(2(r^2 - mu))."""
    k = sp.sqrt(3) / 2
    if system == "isotropic":
        t, r, th, ph, ps = s
        Hinv = 1 / (1 + mu / r ** 2)
        a = [0, 0, 0, j * sp.cos(th) / (2 * r ** 2), j / (2 * r ** 2)]
    elif system == "euler":
        t, r, al, be, ga = s
        Hinv = 1 - mu / r ** 2
        a = [0, 0, j * sp.cos(be) / (2 * (r ** 2 - mu)), 0, j / (2 * (r ** 2 - mu))]
    else:
        raise ValueError(system)
    e0 = [Hinv * (1 + a[0])] + [Hinv * x for x in a[1:]]
    return [k * (e0[0] - 1)] + [k * x for x in e0[1:]]


def check_chart(chart, system):
    """Each chart at points against what its sources state.

    The two charts of Breckenridge, Myers, Peet and Vafa's radius with the Euler angles, and of the
    isotropic radius, solve the bosonic equations of minimal supergravity of Gauntlett, Myers and
    Townsend's Lagrangian, L = R - F^2 - (2/3 sqrt 3) eps A F F, with their potential: the Einstein
    equation R_ab = 2 F_ac F_b^c - (1/3) g_ab F^2, and Maxwell's equation with the Chern-Simons term,
    sqrt(-g) nabla_a F^ab = (1/2 sqrt 3) eps^{bcdef} F_cd F_ef up to the sign the orientation fixes, the same
    sign at every event. Their Ricci scalar is Gauntlett, Myers and Townsend's eq. (singu),
    R = 2(2 sin^2 beta - lambda^2)/lambda^8 in units mu = 1, which is 2(2j^2 - mu^2 R^2)/R^8 in the radius
    R of Breckenridge, Myers, Peet and Vafa. The Hopf chart is the Euler one pulled back along
    beta = 2 theta, alpha = psi + phi, gamma = psi - phi, and the ingoing chart is the Euler one pulled
    back along Gauntlett, Myers and Townsend's c dt = c dv - sqrt(1 - j^2/r^6) dr/(1 - mu/r^2)^2 and
    d gamma = d psitilde - 2j dr/((1 - mu/r^2) r sqrt(r^6 - j^2)), with (alpha, beta) their (phi, theta)."""
    import verify_metrics as vm
    name = f"{METRIC}/{system}"
    s = list(chart.symbols)
    g = chart.geo.g
    P = chart.reader.parameters
    mu, j = P["mu"], P["j"]

    def vanishes(value, what, at):
        number = complex(sp.N(sp.sympify(value).subs(at), 30))
        if abs(number) > 1e-18:
            raise AssertionError(f"{name}: {what} at {at}, off by {number}")

    def events():
        for mu0, j0 in CASES:
            for r0 in RADII:
                radius = sp.sqrt(r0 ** 2 - mu0) if system == "isotropic" else r0
                for angles in ANGLE_POINTS:
                    yield {mu: mu0, j: j0, s[0]: sp.Rational(1, 3), s[1]: radius,
                           **dict(zip(s[2:], angles))}

    if system in ("bmpv", "ingoing"):
        source = PRINT_CHART("euler")
        t, r, al, be, ga = source.symbols
        smu, sj = source.reader.parameters["mu"], source.reader.parameters["j"]
        if system == "bmpv":
            _, rr, th, ph, ps = s
            image = [s[0], rr, ps + ph, 2 * th, ps - ph]
            J = sp.Matrix(5, 5, lambda a, b: sp.diff(image[a], s[b]))
        else:
            v, rr, th, ph, pst = s
            image = [v, rr, ph, th, pst]
            J = sp.Matrix(5, 5, lambda a, b: sp.diff(image[a], s[b]))
            J[0, 1] = -sp.sqrt(1 - j ** 2 / rr ** 6) / (1 - mu / rr ** 2) ** 2
            J[4, 1] = -2 * j / ((1 - mu / rr ** 2) * rr * sp.sqrt(rr ** 6 - j ** 2))
        at = dict(zip(source.symbols, image))
        at.update({smu: mu, sj: j})
        pulled = J.T * source.geo.g.subs(at, simultaneous=True) * J
        for e in events():
            for a in range(5):
                for b in range(a, 5):
                    vanishes(pulled[a, b] - g[a, b],
                             f"the Euler chart pulled back misses slot {chart.coords_tex[a]}{chart.coords_tex[b]}", e)
        return

    A = potential(system, s, mu, j)
    F = sp.Matrix(5, 5, lambda a, b: sp.diff(A[b], s[a]) - sp.diff(A[a], s[b]))
    ginv = chart.geo.ginv
    Fup = ginv * F * ginv
    ricci = chart.geo.ricci_ll()
    gamma = chart.geo.christoffel_ull()
    divergence = [sum(sp.diff(Fup[m, x], s[m]) for m in range(5))
                  + sum(vm._at(gamma, (m, m, k)) * Fup[k, x] for m in range(5) for k in range(5))
                  for x in range(5)]
    scalar = chart.geo.ricci_scalar()
    sign = None
    for e in events():
        gn = g.subs(e).evalf(30)
        Fn = F.subs(e).evalf(30)
        Fupn = Fup.subs(e).evalf(30)
        ginvn = ginv.subs(e).evalf(30)
        F2 = sum(Fn[a, b] * Fupn[a, b] for a in range(5) for b in range(5))
        FF = Fn * ginvn * Fn.T
        for a in range(5):
            for b in range(a, 5):
                value = vm._at(ricci, (a, b))
                vanishes(sp.sympify(value).subs(e).evalf(30) - 2 * FF[a, b] + gn[a, b] * F2 / 3,
                         f"the Einstein equation of minimal supergravity fails in slot {a}{b}", {})
        root = sp.sqrt(-gn.det())
        for x in range(5):
            others = [i for i in range(5) if i != x]
            source_term = 0
            for perm in itertools.permutations(others):
                order = [x, *perm]
                sgn = (-1) ** sum(1 for a in range(5) for b in range(a + 1, 5) if order[a] > order[b])
                source_term += sgn * Fn[perm[0], perm[1]] * Fn[perm[2], perm[3]]
            lhs = root * sp.sympify(divergence[x]).subs(e).evalf(30)
            if abs(complex(source_term)) < 1e-20:
                vanishes(lhs, f"Maxwell's equation fails where the Chern-Simons source vanishes, x^{x}", {})
                continue
            ratio = complex(lhs / source_term) * 2 * float(sp.sqrt(3))
            here = 1 if ratio.real > 0 else -1
            sign = here if sign is None else sign
            if abs(ratio - sign) > 1e-15:
                raise AssertionError(f"{name}: Maxwell's equation with the Chern-Simons term fails in x^{x} "
                                     f"at {e}: sqrt(-g) div F / (eps F F) is {ratio / (2 * float(sp.sqrt(3)))}")
        radius = e[s[1]] if system == "euler" else sp.sqrt(e[s[1]] ** 2 + e[mu])
        want = 2 * (2 * e[j] ** 2 - e[mu] ** 2 * radius ** 2) / radius ** 8
        value = sp.sympify(scalar)
        vanishes(value.subs(e).evalf(30) - want, "the Ricci scalar is not Gauntlett, Myers and Townsend's", {})
    if sign is None:
        raise AssertionError(f"{name}: the Chern-Simons source vanished at every event, so nothing was compared")
    ORIENTATION[system] = sign


def PRINT_CHART(system):
    """A chart built by the printer's own Chart, for the pullbacks."""
    import chart_printer as cp
    spec = chart_spec(system)
    return cp.Chart(spec["system"]["coords"], spec["system"]["parameters"], spec["chart_line_element"])
