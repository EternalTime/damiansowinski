"""Krasnikov's singularity-free time machine that respects the weak energy condition,
Class. Quantum Grav. 15, 997 (1998), gr-qc/9702031, in his signature (-,+,+,+) and with no unit fixed.

time_machine  L_N = Q_N x S, his eqs. (2), (5) and (7): ds^2 = f^-2 (-dtau^2 + dchi^2) + R^2 dOmega^2 on
              the square |tau|, |chi| < m, eq. (1), with f a smooth bounded function that vanishes only
              at (+-h, +-h), eq. (3), repeats near the upper cut its values near the lower one, eq. (4),
              and R bounded by eq. (6); the four points are removed and the cuts tau = +-h, |chi| < h
              glued, the upper bank of each to the lower bank of the other.
universe      N on R x S^3, his eq. (9): ds^2 = a^2 (-dtau^2 + dchi^2 + rho^2 dOmega^2), with
              -pi/2 <= chi <= pi/2, a = 1/f and rho = fR on the square, eq. (10a), and a = a(tau),
              rho = rho(chi) off it, eq. (10b), with rho = cos chi where |chi| > n.

f, a and rho are free functions of tau and chi, so every value is printed in them and their
derivatives. krasnikov_time_machine.md beside this file is the derivation.
"""

METRIC = "krasnikov_time_machine"
CHARTS = ("time_machine", "universe")
COORDS = [r"\tau", r"\chi", r"\theta", r"\phi"]
DIMENSIONS = {
    (METRIC, "time_machine"): {r"\tau": "1", r"\chi": "1", r"\theta": "1", r"\phi": "1", "f": "1/L", "R": "L",
                               "h": "1", "m": "1"},
    (METRIC, "universe"): {r"\tau": "1", r"\chi": "1", r"\theta": "1", r"\phi": "1", "a": "L", r"\rho": "1",
                           "h": "1", "m": "1", "n": "1"},
}
NAMES = {"time_machine": "Time Machine", "universe": "Almost Friedmann Universe"}
PARAMETERS = {"time_machine": [r"f = f(\tau,\chi)", "R", "h", "m"],
              "universe": [r"a = a(\tau,\chi)", r"\rho = \rho(\tau,\chi)", "h", "m", "n"]}
CUTS = (r"\tau = \pm h,\; |\chi| < h \;\text{(the two cuts, the upper bank of each glued to the lower bank"
        r" of the other)}")
DOMAINS = {
    "time_machine": [r"\tau \in (-m, m)", r"\chi \in (-m, m)", r"\theta \in [0, \pi]", r"\phi \in [0, 2\pi)",
                     r"(\tau, \chi) \neq (\pm h, \pm h)", CUTS,
                     r"|\tau| < h,\; |\chi| < 2h - |\tau| \;\text{(the causality violating set)}"],
    "universe": [r"\tau \in (-\infty, \infty)", r"\chi \in [-\pi/2, \pi/2]", r"\theta \in [0, \pi]",
                 r"\phi \in [0, 2\pi)",
                 r"a > 0 \;\text{(between the cosmological singularities, where } a = 0\text{)}",
                 r"|\tau|, |\chi| < m \;\text{(the time machine)}",
                 r"(\tau, \chi) \neq (\pm h, \pm h)", CUTS],
}
SPHERE = r"\left(d\theta^2 + \sin^2\theta\,d\phi^2\right)"
LINES = {"time_machine": r"ds^2 = f^{-2}\left(-d\tau^2 + d\chi^2\right) + R^2" + SPHERE,
         "universe": r"ds^2 = a^2\left(-d\tau^2 + d\chi^2 + \rho^2" + SPHERE + r"\right)"}


def chart_spec(system):
    """One of the two charts, its values printed in its functions and their derivatives."""
    import sympy as sp
    import verify_metrics as vm

    probe = vm.Reader(COORDS, PARAMETERS[system], ())
    tau, chi = probe.symbol[r"\tau"], probe.symbol[r"\chi"]
    D = sp.Derivative
    functions = [probe.parameters["f"]] if system == "time_machine" else [probe.parameters["a"],
                                                                          probe.parameters["rho"]]
    lead = list(functions)
    for F in functions:
        lead += [D(F, tau), D(F, chi)]
    for F in functions:
        lead += [D(F, (tau, 2)), D(F, tau, chi), D(F, (chi, 2))]
    if system == "time_machine":
        lead.append(probe.parameters["R"])
    return {"metric_id": METRIC,
            "system": {"id": system, "name": NAMES[system], "coords": COORDS, "domains": DOMAINS[system],
                       "parameters": PARAMETERS[system], "line_element": LINES[system]},
            "chart_line_element": LINES[system],
            "printer": {"lead": lead, "factors": lead},
            "check": lambda chart: chart_check(chart, system)}


def wec_bound(f, tau, chi):
    """The right side of Krasnikov's eq. (6) at a point, f (f_chichi - f_tautau) + f_tau^2 - f_chi^2."""
    import sympy as sp

    return (f * (sp.diff(f, chi, 2) - sp.diff(f, tau, 2)) + sp.diff(f, tau) ** 2 - sp.diff(f, chi) ** 2)


def orthonormal(G, g):
    """The Einstein tensor of a diagonal metric in its orthonormal frame, G_ab / sqrt(|g_aa g_bb|)."""
    import sympy as sp

    return sp.Matrix(4, 4, lambda i, j: G[i, j] / sp.sqrt(abs(g[i, i] * g[j, j])))


def chart_check(chart, system):
    """The chart at points against what Krasnikov states.

    time_machine: the metric is his eqs. (2), (5) and (7) slot by slot; in the orthonormal frame the
    Einstein tensor is diag(1/R^2, -1/R^2, -K, -K) with K the right side of his eq. (6), so the
    energy density is 1/R^2, its sum with the pressure along chi vanishes and its sum with either
    pressure on the sphere is 1/R^2 - K, and the weak energy condition holds exactly where
    K <= R^-2. The Ricci scalar is 2K + 2/R^2.

    universe: on the square, a = 1/f and rho = fR give the time machine's metric, eq. (10a). Off it,
    with a(tau) and rho(chi), eq. (10b), the energy density and its sums with the pressures are,
    times a^2, 3 Phi'^2 + (1 - rho'^2)/rho^2 - 2 rho''/rho, 2 Phi'^2 - 2 Phi'' - 2 rho''/rho and
    2 Phi'^2 - 2 Phi'' + (1 - rho'^2)/rho^2 - rho''/rho with Phi = ln a, each non-negative under his
    eq. (11), Phi'' <= 0, rho''/rho <= 0 and 1 - rho'^2 >= 0; where rho = cos chi the three
    pressures agree and the universe is Friedmann's closed one, with a^2 times the energy density
    3 (Phi'^2 + 1)."""
    import sympy as sp

    X = chart.symbols
    tau, chi, theta, phi = X
    P = chart.reader.parameters
    g = sp.Matrix(chart.geo.g)
    G = sp.Matrix(chart.geo.einstein_ll())
    where = f"{METRIC}/{system}"
    angles = [{theta: sp.Rational(1, 3)}, {theta: sp.Rational(6, 5)}, {theta: sp.Rational(5, 2)}]

    def vanishes(e, what, points):
        for at in points:
            if abs(complex(sp.N(sp.sympify(e).subs(at), 40))) > 1e-28:
                raise AssertionError(f"{where}: {what} fails at {at}")

    def nonnegative(e, what, points):
        for at in points:
            if sp.N(sp.sympify(e).subs(at), 40) < -1e-28:
                raise AssertionError(f"{where}: {what} is negative at {at}")

    if system == "time_machine":
        f, R = P["f"], P["R"]
        theirs = sp.diag(-f ** -2, f ** -2, R ** 2, R ** 2 * sp.sin(theta) ** 2)
        # Two smooth functions held in for f, neither one Krasnikov's, since the identities hold
        # for any f: one with his zeros at (+-1/3, +-1/3), and one with none.
        h = sp.Rational(1, 3)
        trials = [((tau ** 2 - h ** 2) ** 2 + (chi ** 2 - h ** 2) ** 2) * (2 + sp.sin(tau + 2 * chi)),
                  sp.exp(tau * chi / 3) + sp.cos(tau - chi) / 5 + 2]
        points = [{tau: sp.Rational(1, 7), chi: -sp.Rational(2, 9), R: sp.Rational(3, 2)},
                  {tau: -sp.Rational(5, 11), chi: sp.Rational(3, 8), R: sp.Rational(7, 5)},
                  {tau: sp.Rational(2, 5), chi: sp.Rational(1, 2), R: sp.Rational(1, 3)}]
        for trial in trials:
            on = lambda e: sp.sympify(e).subs(f, trial).doit()  # noqa: E731
            K = wec_bound(trial, tau, chi)
            gt = theirs.subs(f, trial)
            Gt = orthonormal(G.applyfunc(on), gt)
            stated = sp.diag(1 / R ** 2, -1 / R ** 2, -K, -K)
            at = [{**p, **q} for p in points for q in angles]
            for i in range(4):
                for j in range(4):
                    vanishes(on(g[i, j]) - gt[i, j], f"his eqs. (2), (5) and (7), slot {i}{j}", at)
                    vanishes(Gt[i, j] - stated[i, j], f"the orthonormal Einstein tensor, slot {i}{j}", at)
            vanishes(on(chart.geo.ricci_scalar()) - 2 * K - 2 / R ** 2, "the Ricci scalar 2K + 2/R^2", at)
            vanishes((Gt[0, 0] + Gt[2, 2]) - (1 / R ** 2 - K), "the energy density plus the pressure on the sphere",
                     at)
        return

    a, rho = P["a"], P["rho"]
    theirs = a ** 2 * sp.diag(-1, 1, rho ** 2, rho ** 2 * sp.sin(theta) ** 2)
    for i in range(4):
        for j in range(4):
            vanishes(g[i, j] - theirs[i, j], f"his eq. (9), slot {i}{j}",
                     [{tau: 0, chi: 0, theta: sp.Rational(1, 3)}])
    # Eq. (10a): on the square, a = 1/f and rho = fR give the time machine's metric, at any f.
    f = sp.Function("f")(tau, chi)
    R = sp.Symbol("R", positive=True)
    trial = sp.exp(tau * chi / 3) + sp.cos(tau - chi) / 5 + 2
    machine = sp.diag(-f ** -2, f ** -2, R ** 2, R ** 2 * sp.sin(theta) ** 2)
    square = [{tau: sp.Rational(1, 7), chi: -sp.Rational(2, 9), R: sp.Rational(3, 2), theta: sp.Rational(1, 3)},
              {tau: -sp.Rational(5, 11), chi: sp.Rational(3, 8), R: sp.Rational(7, 5), theta: sp.Rational(6, 5)}]
    for i in range(4):
        vanishes(g[i, i].subs({a: 1 / f, rho: f * R}).subs(f, trial) - machine[i, i].subs(f, trial),
                 f"his eq. (10a), slot {i}{i}", square)
    # Eq. (10b) and (11): a of tau alone with ln a bending down, rho of chi alone, positive, bending
    # down, and with |rho'| <= 1; then where rho = cos chi, Friedmann's closed universe.
    outer = [{tau: sp.Rational(1, 4), chi: sp.Rational(7, 10), theta: sp.Rational(1, 3)},
             {tau: -sp.Rational(9, 10), chi: -sp.Rational(13, 10), theta: sp.Rational(6, 5)},
             {tau: sp.Rational(3, 2), chi: sp.Rational(1, 5), theta: sp.Rational(5, 2)}]
    scale = 3 * sp.cos(tau / 2) ** 2 * (2 + sp.tanh(tau) / 3)
    for profile in (2 * (sp.cos(chi) + sp.Rational(1, 2)) / 3, sp.cos(chi)):
        on = lambda e: sp.sympify(e).subs({a: scale, rho: profile}).doit()  # noqa: E731
        gt = theirs.subs({a: scale, rho: profile})
        Gt = orthonormal(G.applyfunc(on), gt)
        Phi = sp.log(scale)
        d1, d2 = sp.diff(Phi, tau), sp.diff(Phi, tau, 2)
        r1, r2 = sp.diff(profile, chi), sp.diff(profile, chi, 2)
        nonnegative(-d2, "his Phi'' <= 0, eq. (11)", outer)
        nonnegative(-r2 / profile, "his rho''/rho <= 0, eq. (11)", outer)
        nonnegative(1 - r1 ** 2, "his 1 - rho'^2 >= 0, eq. (11)", outer)
        energy = 3 * d1 ** 2 + (1 - r1 ** 2) / profile ** 2 - 2 * r2 / profile
        along = 2 * d1 ** 2 - 2 * d2 - 2 * r2 / profile
        across = 2 * d1 ** 2 - 2 * d2 + (1 - r1 ** 2) / profile ** 2 - r2 / profile
        for i in range(4):
            for j in range(4):
                if i != j:
                    vanishes(Gt[i, j], f"the Einstein tensor off the diagonal, slot {i}{j}", outer)
        vanishes(scale ** 2 * Gt[0, 0] - energy, "the energy density", outer)
        vanishes(scale ** 2 * (Gt[0, 0] + Gt[1, 1]) - along, "the energy density plus the pressure along chi", outer)
        for k in (2, 3):
            vanishes(scale ** 2 * (Gt[0, 0] + Gt[k, k]) - across,
                     f"the energy density plus the pressure on the sphere, slot {k}{k}", outer)
        for name, value in (("the energy density", energy), ("its sum with the pressure along chi", along),
                            ("its sum with the pressure on the sphere", across)):
            nonnegative(value, f"the weak energy condition, {name},", outer)
        if profile == sp.cos(chi):
            vanishes(Gt[1, 1] - Gt[2, 2], "Friedmann's equal pressures", outer)
            vanishes(scale ** 2 * Gt[0, 0] - 3 * (d1 ** 2 + 1), "Friedmann's energy density", outer)
