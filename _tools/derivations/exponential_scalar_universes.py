"""Exact flat exponential scalar cosmologies, in units c = 1.

Copeland, Liddle and Wands (1998), (1)-(4), Table I and (A7)-(A8), and
Heard and Wands (2002), (10)-(11) and (22)-(24), fix the two matter branches.
Abbott and Wise (1984), (2.1)-(2.3), fix the inflationary conformal map.
The continuation to deceleration and the logarithmic coasting map follow
by integrating dt/a. See exponential_scalar_universes.md for the domains.
"""

METRIC = "exponential_scalar_universes"
CHARTS = ("cosmic_cartesian", "cosmic_spherical", "conformal", "coasting")
SCALES = {
    "cosmic_cartesian": r"a = a_0\left(t/t_0\right)^p",
    "cosmic_spherical": r"a = a_0\left(t/t_0\right)^p",
    "conformal": r"a = a_0\left((1-p)a_0\eta/t_0\right)^{p/(1-p)}",
    "coasting": r"a = a_0e^{a_0\eta/t_0}",
}
COORDS = {s: (["t", "r", r"\theta", r"\phi"] if s == "cosmic_spherical" else
              ["t" if s.startswith("cosmic") else r"\eta", "x", "y", "z"]) for s in CHARTS}
DIMENSIONS = {(METRIC, s): {**{x: ("1" if x in (r"\theta", r"\phi") else "L") for x in COORDS[s]},
                            "p": "1", "a_0": "1", "t_0": "L", "a": "1"} for s in CHARTS}
HELD = {(METRIC, s): ("a",) for s in CHARTS}
RATES = {(METRIC, s): {"a": {COORDS[s][0]: ("pa/t" if s.startswith("cosmic") else
                    r"pa/((1-p)\eta)" if s == "conformal" else "a_0 a/t_0")}} for s in CHARTS}


def chart_spec(system):
    """Four source-defined charts with a positive scale factor held by its slope."""
    import verify_metrics as vm

    cosmic = system.startswith("cosmic")
    params = (["p"] if system != "coasting" else []) + ["a_0", "t_0", SCALES[system]]
    if system == "cosmic_spherical":
        line = r"ds^2 = -dt^2+a^2\left(dr^2+r^2d\theta^2+r^2\sin^2\theta\,d\phi^2\right)"
        space = [r"r \in [0,\infty)", r"\theta \in [0,\pi]", r"\phi \in [0,2\pi)"]
    else:
        line = (r"ds^2 = -dt^2+a^2\left(dx^2+dy^2+dz^2\right)" if cosmic else
                r"ds^2 = a^2\left(-d\eta^2+dx^2+dy^2+dz^2\right)")
        space = [x+r" \in (-\infty,\infty)" for x in ("x", "y", "z")]
    domain = ([r"t \in (0,\infty)"] if cosmic else
              [r"\eta \in (-\infty,0)\;\text{for}\;p>1",
               r"\eta \in (0,\infty)\;\text{for}\;p<1"] if system == "conformal" else
              [r"\eta \in (-\infty,\infty)"])
    key = (METRIC, system)
    probe = vm.Reader(COORDS[system], params, (), held=HELD[key], rates=RATES[key])
    return {"metric_id": METRIC,
            "system": {"id": system, "name": {"cosmic_cartesian": "Cosmic Cartesian",
                        "cosmic_spherical": "Cosmic Spherical", "conformal": "Conformal Time",
                        "coasting": "Coasting Conformal Time"}[system],
                       "coords": COORDS[system], "parameters": params, "domains": domain+space,
                       "line_element": line},
            "chart_line_element": line, "printer": {},
            "reduce": lambda value: vm.norm(probe.by_rates(value)),
            "check": lambda chart: chart_check(chart, system)}


def chart_check(chart, system):
    """Point checks of the coordinate maps and the two exact matter branches."""
    import sympy as sp

    R = chart.reader
    P = R.parameters
    X = chart.symbols
    clock = X[0]
    a0, t0 = P['a_0'], P['t_0']
    p = P.get('p', sp.Integer(1))
    cosmic = system.startswith('cosmic')
    t = (clock if cosmic else t0*sp.exp(a0*clock/t0) if system == 'coasting' else
         t0*((1-p)*a0*clock/t0)**(1/(1-p)))
    a = a0*(t/t0)**p
    lapse = sp.diff(t, clock)
    if system == 'cosmic_spherical':
        r, theta, phi = X[1:]
        image = [t, r*sp.sin(theta)*sp.cos(phi), r*sp.sin(theta)*sp.sin(phi), r*sp.cos(theta)]
    else:
        image = [t, *X[1:]]
    J = sp.Matrix(image).jacobian(X)
    pulled = J.T*sp.diag(-1, a*a, a*a, a*a)*J
    G = sp.Matrix(chart.geo.einstein_ll())
    # Scalar-only and positive-fluid examples span inflation, deceleration and coasting.
    for pv, tv, a0v, t0v in ((2, 0.7, 1.2, 0.9), (2/3, 1.3, 0.8, 1.4),
                               (1, 2.1, 1.1, 0.7), (0.4, 0.9, 1.4, 1.2)):
        if system == 'coasting':
            pv = 1.0
        if system == 'conformal' and pv == 1:
            continue
        cv = (tv if cosmic else t0v/a0v*__import__('math').log(tv/t0v) if system == 'coasting' else
              t0v*(tv/t0v)**(1-pv)/(a0v*(1-pv)))
        at = {clock: cv, X[1]: 0.8, X[2]: 0.9, X[3]: 0.4, a0: a0v, t0: t0v, R.c: 1}
        if 'p' in P:
            at[p] = pv
        def number(v):
            return complex(sp.N(sp.sympify(v).subs(R.held).doit().subs(at), 35))
        def near(values, label):
            if max((abs(number(v)) for v in values), default=0) > 2e-8:
                raise AssertionError(f'{METRIC}/{system}: {label}')
        near(list(chart.geo.g-pulled), 'cosmic Cartesian pullback')
        for branch, lam2 in (('scalar', 2/pv), ('fluid', 9/pv)):
            kappa, lam = 1.3, lam2**0.5
            gamma = 2/(3*pv)
            rho = 0 if branch == 'scalar' else 3*pv*(pv-2/lam2)/(kappa*kappa*t*t)
            V = 2*(3*pv-1)/(kappa*kappa*lam2*t*t)
            gradient = sp.Matrix([2*sp.diff(t, x)/(lam*kappa*t) for x in X])
            velocity = sp.Matrix([-lapse, 0, 0, 0])
            # (grad phi)^2 = -dot(phi)^2 in these comoving charts.
            scalar = gradient*gradient.T-chart.geo.g*(-2/(kappa*kappa*lam2*t*t)+V)
            fluid = gamma*rho*velocity*velocity.T+(gamma-1)*rho*chart.geo.g
            near(list(G-kappa*kappa*(scalar+fluid)), branch+' Einstein source')
            near([-2/(lam*kappa*t*t)+6*pv/(lam*kappa*t*t)-lam*kappa*V],
                 branch+' scalar equation')
            near([-2*rho/t+3*gamma*pv*rho/t], branch+' fluid conservation')
        for name, rates in RATES[(METRIC, system)].items():
            for coord, rate in rates.items():
                near([R(rate)-sp.diff(R.held[P[name]], R.symbol[coord])], 'scale factor derivative')
