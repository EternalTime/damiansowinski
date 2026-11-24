"""Chamseddine and Volkov's monopole of N = 4 gauged supergravity, hep-th/9707176, Phys. Rev. Lett.
79, 3343 (1997), in their units 4 pi G = c = 1 with the gauge coupling one, their signature
(+,-,-,-) turned to (-,+,+,+).

The line element is their eq. (13), TeX label 10,
ds^2 = a^2 (sinh rho / R)(-dt^2 + d rho^2 + R^2 dOmega^2), with R^2 = 2 rho coth rho - rho^2/sinh^2 rho - 1
their eq. (12), TeX label R, and the fields are their eq. (14), TeX label 11: the Yang-Mills amplitude
w = rho / sinh rho and the dilaton e^{2 phi} = a^2 sinh rho / (2R). R is held as a function of rho by its
slope. chamseddine_volkov_monopole.md beside this file is the derivation, and records why the
Schwarzschild coordinates of their companion paper, hep-th/9711181 eq. (5.9), are no chart of the page.
"""

METRIC = "chamseddine_volkov_monopole"
CHARTS = ("chamseddine_volkov",)
COORDS = ["t", r"\rho", r"\theta", r"\varphi"]
R_DEFINITION = r"R = \sqrt{2\rho\coth\rho - \dfrac{\rho^2}{\sinh^2\rho} - 1}"
DIMENSIONS = {(METRIC, s): {"t": "1", r"\rho": "1", r"\theta": "1", r"\varphi": "1", "a": "L", "R": "1"}
              for s in CHARTS}
HELD = {(METRIC, s): ("R",) for s in CHARTS}
RATES = {(METRIC, s): {"R": {r"\rho": r"\dfrac{\cosh\rho\,\sinh^2\rho - 2\rho\sinh\rho + \rho^2\cosh\rho}{R\sinh^3\rho}"}}
         for s in CHARTS}
NAMES = {"chamseddine_volkov": "Chamseddine-Volkov"}
DOMAINS = [r"t \in (-\infty, \infty)", r"\rho \in [0, \infty)", r"\theta \in [0, \pi]", r"\varphi \in [0, 2\pi)"]
LINE = (r"ds^2 = \dfrac{a^2\sinh\rho}{R}\left(-dt^2 + d\rho^2 + R^2\left(d\theta^2"
        r" + \sin^2\theta\,d\varphi^2\right)\right)")


def chart_spec(system):
    """The one chart, its values printed in R, rho and the hyperbolic functions of rho."""
    import chart_printer as cp
    import verify_metrics as vm

    key = (METRIC, system)
    parameters = ["a", R_DEFINITION]
    probe = vm.Reader(COORDS, parameters, (), held=HELD[key], rates=RATES[key])
    rho, R = probe.symbol[r"\rho"], probe.parameters["R"]
    return {"metric_id": METRIC,
            "system": {"id": system, "name": NAMES[system], "coords": COORDS, "domains": DOMAINS,
                       "parameters": parameters, "line_element": LINE},
            "chart_line_element": LINE,
            "printer": {"lead": [R, rho], "flip": False},
            "pretty": cp.hyperbolic(rho),
            "reduce": lambda value: vm.norm(probe.by_rates(value)),
            "check": lambda chart: chart_check(chart, system)}


def fields(rho, a):
    """The Yang-Mills amplitude, the dilaton and R of their eqs. (12) and (14), written out."""
    import sympy as sp

    R = sp.sqrt(2 * rho * sp.cosh(rho) / sp.sinh(rho) - rho ** 2 / sp.sinh(rho) ** 2 - 1)
    return rho / sp.sinh(rho), sp.log(a ** 2 * sp.sinh(rho) / (2 * R)) / 2, R


def stress(g, rho, theta, a):
    """2 T_{mu nu} of their action, eq. (1), in the signature (-,+,+,+): the dilaton phi with the
    potential -e^{-2 phi}/8, so the cosmological term -e^{-2 phi}/4, and the purely magnetic SU(2) field
    with F_{rho theta} = w', F_{rho varphi} = w' sin(theta) and F_{theta varphi} = (1 - w^2) sin(theta)
    in its three directions of the gauge group, with e^{2 phi} before F^2. -R/4 in the action is
    4 pi G = 1, so G_{mu nu} = 2 T_{mu nu}."""
    import sympy as sp

    w, phi, _ = fields(rho, a)
    wp, dphi = sp.diff(w, rho), sp.diff(phi, rho)
    grr, gthth = g[1, 1], g[2, 2]
    squares = [0, 2 * wp ** 2 / gthth, wp ** 2 / grr + (1 - w ** 2) ** 2 / gthth,
               (wp ** 2 / grr + (1 - w ** 2) ** 2 / gthth) * sp.sin(theta) ** 2]
    F2 = 2 * (2 * wp ** 2 / (grr * gthth) + (1 - w ** 2) ** 2 / gthth ** 2)
    potential = -sp.exp(-2 * phi) / 8
    lagrangian = dphi ** 2 / grr / 2 + potential + sp.exp(2 * phi) * F2 / 4
    return [2 * ((dphi ** 2 if i == 1 else 0) + sp.exp(2 * phi) * squares[i] - g[i, i] * lagrangian)
            for i in range(4)]


def chart_check(chart, system):
    """The chart at points against what its sources state.

    The metric is their eq. (13) with the signature turned; it and the fields of eq. (14) solve the
    field equations of their action, eq. (1), G_{mu nu} = 2 T_{mu nu}, slot by slot. The sphere at rho
    has the areal radius r = a sqrt(R sinh rho) and the Schwarzschild function N = g^{rr} = rho^2/R^2,
    their companion paper's eq. (5.9), hep-th/9711181."""
    import sympy as sp

    X = chart.symbols
    t, rho, theta, varphi = X
    a = chart.reader.parameters["a"]
    held = chart.reader.held
    g = sp.Matrix(chart.geo.g).xreplace(held)
    G = sp.Matrix(chart.geo.einstein_ll()).xreplace(held)
    w, phi, R = fields(rho, a)
    theirs = a ** 2 * sp.sinh(rho) / R * sp.diag(-1, 1, R ** 2, R ** 2 * sp.sin(theta) ** 2)
    twice_stress = stress(theirs, rho, theta, a)
    r = a * sp.sqrt(R * sp.sinh(rho))
    N = sp.diff(r, rho) ** 2 / theirs[1, 1]
    points = [{rho: sp.Rational(7, 10), a: sp.Rational(13, 10), theta: sp.Rational(1, 3)},
              {rho: sp.Rational(23, 10), a: sp.Rational(1, 2), theta: sp.Rational(6, 5)},
              {rho: 4, a: 3, theta: sp.Rational(5, 7)}]

    def vanishes(e, what):
        for at in points:
            if abs(complex(sp.N(sp.sympify(e).subs(at), 40))) > 1e-28:
                raise AssertionError(f"{METRIC}/{system}: {what} fails at {at}")

    for i in range(4):
        for j in range(4):
            vanishes(g[i, j] - theirs[i, j], f"their eq. (13), slot {i}{j}")
            if i != j:
                vanishes(G[i, j], f"the Einstein tensor off the diagonal, slot {i}{j}")
        vanishes(G[i, i] - twice_stress[i], f"G = 2T for their eq. (1), slot {i}{i}")
    vanishes(g[2, 2] - r ** 2, "the areal radius a sqrt(R sinh rho)")
    vanishes(N - rho ** 2 / R ** 2, "their eq. (5.9), N = rho^2/R^2")
