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


MATTER_TEXT = {"inflation": r"Scalar inflation ($\lambda=1$)",
               "tracking": r"Dust scaling ($\gamma=1$, $\lambda=3$, $\Omega_\phi=1/3$)",
               "coasting": r"Coasting scalar ($\lambda^2=2$)"}

EXAMPLES = (("inflation", "Scalar inflation", 2.0),
            ("tracking", "Dust scaling", 2/3),
            ("coasting", "Coasting scalar", 1.0))


def eta_of(t, p):
    """Conformal time in drawing units a_0 = t_0 = 1."""
    import numpy as np
    t = np.asarray(t, float)
    return np.log(t) if p == 1 else t**(1-p)/(1-p)


def embedding_views(ck, src, em):
    """Three individual moments t=t_0, each an isometric Euclidean plane."""
    import numpy as np

    views = []
    for vid, label, p in EXAMPLES:
        sl = em.Slice(src, METRIC, "cosmic_spherical", "r", r"\phi", {"t": 1, **em.EQUATOR},
                      {"p": p, "a_0": 1, "t_0": 1})
        size = 5.0
        disc = em.Piece("plane", "sheet", sl, 0, 2, 0, 1,
                        (("axis", "$r=0$"), ("edge", r"$r\to\infty$")),
                        [(0.5, "r", None), (1, "r", None), (1.5, "r", None), (2, "r", None)], size)
        ck.plane(label, sl, np.linspace(0.001, 2, 80))
        ck.isometry(label, disc)
        ck.radius(label+": rho = a r", disc, lambda r: r, size)
        ck.form(label+": Z = 0", disc, lambda r: 0*r, size)
        surface = em.Surface([disc])
        fig = em.figure_of([surface], {"sheet": "cover"}, size, meridians=12)
        fig.legend("fill", "cover", r"$t=t_0$, $\theta=\pi/2$")
        fig.legend("line", "r", "$r/t_0=1/2,1,3/2,2$")
        fig.legend("line", "meridian", r"$\phi$ constant")
        data = fig.done()
        data['box'] = [-2.5, 2.5, -2.5, 2.5]
        views.append(em.view(vid, label, "$t_0$", [surface], data,
                            settings='$p='+('2/3' if vid == 'tracking' else f'{p:g}')+'$, $a_0=t_0=1$.'))
    return views


def embedding_captions():
    """Single-moment equatorial planes."""
    return {(METRIC, vid): [r"A flat equatorial plane at $t=t_0$, with $\rho=ar$ and $Z=0$; the disc cut at $r=2t_0$."]
            for vid, _, _ in EXAMPLES}


PLANES = (("cosmic_cartesian", "inflation", 2.0),
          ("cosmic_spherical", "tracking", 2/3),
          ("conformal", "inflation", 2.0),
          ("conformal", "tracking", 2/3),
          ("coasting", "coasting", 1.0))


def diagram_specs(Diagram, sideways):
    """Null geodesics in every chart, with horizons from their time integrals."""
    rows = []
    for system, vid, p in PLANES:
        cosmic = system.startswith('cosmic')
        spherical = system == 'cosmic_spherical'
        clock, radial = COORDS[system][:2]
        params = {'a_0': 1, 't_0': 1, **({'p': p} if system != 'coasting' else {})}
        box = ((0, 2.4, 0, 4.8) if spherical else (-2.4, 2.4, 0, 4.8) if cosmic or p < 1 else
               (-2.4, 2.4, -2.4, 2.4) if p == 1 else (-2.4, 2.4, -4.8, 0))
        fixed = {'theta': 'pi/2', 'phi': '0'} if spherical else {'y': '0', 'z': '0'}
        curve = (() if p == 1 else
                 (('event', '1/t' if cosmic else '-eta', 'event horizon of the observer at $x=0$'),) if p > 1 else
                 (('past', '3*t**(1/3)' if cosmic else 'eta', 'particle horizon of the observer at $r=0$'),))
        rows.append(Diagram(METRIC, system, vid, dict((x[0], x[1]) for x in EXAMPLES)[vid],
                            (clock, radial), box, f'${radial}/t_0$', f'${clock}/t_0$', params, fixed,
                            tau='t' if cosmic else 'eta', mirror=spherical,
                            families=sideways if not spherical else ('ingoing', 'outgoing'),
                            curves=curve, cones=(5, 5), step=0.00025 if spherical else 0.0025))
    return rows


def diagram_captions():
    """Compact descriptions of the null planes."""
    return {(METRIC, system, vid): [MATTER_TEXT[vid]+
        (r", radial null geodesics through $r=0$." if system == 'cosmic_spherical' else
         r", null geodesics on $y=z=0$.")]
        for system, vid, _ in PLANES}


def null_forms():
    """First integrals of the null slopes, checked against traced rays."""
    return {(METRIC, system, vid):
            ((lambda t, r, p=p: eta_of(t, p)+r), (lambda t, r, p=p: eta_of(t, p)-r),
             (lambda t, r: t > 0.02) if system == 'cosmic_spherical' else None)
            if system.startswith('cosmic') else ((lambda e, r: e+r), (lambda e, r: e-r), None)
            for system, vid, p in PLANES}


def slice_marks(slices, system, vid):
    """The embedding moment t=t_0 in the null chart at the same parameters."""
    p = dict((x[0], x[2]) for x in EXAMPLES)[vid]
    clock = 1 if system.startswith('cosmic') else float(eta_of(1, p))
    def line(moment):
        hi = moment.reach('cosmic_spherical', 'r')[1]
        return slices.along(clock, 0 if system == 'cosmic_spherical' else -hi, hi)
    return slices.one(METRIC, line, label="$t=t_0$", view_id=vid)


def conformal_views(ck, src, cf):
    """Compactified radial domains from eta(t), with curvature limits on singular edges."""
    import math
    import numpy as np

    views = []
    for system, vid, p in PLANES:
        cosmic = system.startswith('cosmic')
        spherical = system == 'cosmic_spherical'
        clock, radial = COORDS[system][:2]
        params = {'a_0': 1, 't_0': 1, **({'p': p} if system != 'coasting' else {})}
        fixed = cf.EQUATOR if spherical else {'y': '0', 'z': '0'}
        plane = cf.Plane(src, METRIC, system, (clock, radial), fixed, params)
        def pq(t, r, p=p, cosmic=cosmic):
            return cf.mink_pq(eta_of(t, p) if cosmic else t, r)
        ts = ck.uniform(0.03, 6)
        clocks = ts if cosmic else eta_of(ts, p)
        ck.chart(f'{METRIC}/{system}/{vid}', plane, pq, clocks, ck.uniform(0.01, 6), lambda t, r: (1, 0))
        c1, c2 = (0.01, 0.001) if cosmic else (float(eta_of(t, p)) for t in (0.01, 0.001))
        ck.diverges(f'{system}/{vid}: K diverges toward t=0', plane.kretschmann(c1, 1), plane.kretschmann(c2, 1))
        cfuture = 1e6 if cosmic else float(eta_of(1e6, p))
        ck.limit(f'{system}/{vid}: K tends to zero in the future', [plane.kretschmann(cfuture, 1)], [0], 1e-12)
        pi = math.pi
        lower, upper = ((-pi, 0) if p > 1 else (0, pi) if p < 1 else (-pi, pi))
        box = [-0.25, pi+0.25, lower-0.25, upper+0.25]
        v = cf.View(system+'_'+vid, dict((x[0], x[1]) for x in EXAMPLES)[vid], box, system)
        tri = [[0, -pi], [pi, 0], [0, 0]] if p > 1 else [[0, 0], [pi, 0], [0, pi]] if p < 1 else cf.TRIANGLE
        v.fill('region', tri)
        v.fill('cover', tri)
        past = [[0, -pi], [pi, 0]] if p >= 1 else [[0, 0], [pi, 0]]
        future = [[0, 0], [pi, 0]] if p > 1 else [[0, pi], [pi, 0]]
        v.line('singular', [past], zig=True)
        v.line('scri', [future])
        v.line('centre', [[[0, lower], [0, upper]]])
        rr = cf.spread(0, np.inf, 500, 10)
        for tv in (0.25, 0.5, 1, 2, 4):
            e = float(eta_of(tv, p))
            v.curve('t', *cf.mink_pq(np.full_like(rr, e), rr))
        es = (-rr if p > 1 else rr if p < 1 else np.linspace(-30, 30, 1000))
        for r in (0.5, 1, 2):
            v.curve('r', *cf.mink_pq(es, np.full_like(es, r)))
        if p != 1:
            # The null ray u=0 for a particle horizon, v=0 for an event horizon.
            v.line('horizon', [[[0, 0], [pi/2, -pi/2 if p > 1 else pi/2]]])
            v.legend('horizon', 'event horizon at $V=0$' if p > 1 else 'particle horizon at $U=0$')
        v.legend('r', '$r$ constant')
        v.legend('t', '$t$ constant')
        v.legend('singular', '$t=0$')
        v.legend('scri', r'$\eta=0$, future infinity' if p > 1 else 'future null infinity')
        v.legend('centre', '$r=0$')
        moment = cf.slices.moments(METRIC, vid, label="$t=t_0$")[0]
        reach = np.linspace(*moment.reach('cosmic_spherical', 'r'), 81)
        v.slice(moment, [cf.mink_pq(np.full_like(reach, float(eta_of(1, p))), reach)])
        v.set(settings='$p='+('2/3' if vid == 'tracking' else f'{p:g}')+'$, $a_0=t_0=1$, '+r'$U,V=\arctan(\eta\mp r)$, $T=U+V$, $X=V-U$.')
        views.append(v)
    return views


def conformal_captions():
    """Each point represents a sphere, with the initial and future limits explicit."""
    return {(METRIC, system+'_'+vid): [MATTER_TEXT[vid]+r", each point a 2-sphere."]
            for system, vid, _ in PLANES}


def curvature_views(ck, src, cv):
    """Line figures of stored curvature scalars at the three drawing parameters."""
    views = []
    for vid, label, p in EXAMPLES:
        xaxis = cv.axis('$t/t_0$', 0.25, 4)
        scales = {'K': {'kind': 'log', 'from': 0.01, 'to': 1e5},
                  'R': {'kind': 'log', 'from': 0.01, 'to': 1e4}}
        panels = [{'id': field, 'name': '$'+field+('t_0^4$' if field == 'K' else 't_0^2$'),
                   'aspect': 2, 'scale': scales[field], 'ticks': cv.value_ticks(scales[field])}
                  for field in ('K', 'R')]
        cv.shape_checks(ck, METRIC+'/'+vid, [x['aspect'] for x in panels])
        params = {'p': p, 'a_0': 1, 't_0': 1}
        curves = []
        for field, scalar in (('K', 'kretschmann'), ('R', 'ricci_scalar')):
            stored = cv.nr_published(src, METRIC, 'cosmic_cartesian', scalar, ['t'], list(params))
            curves += cv.line_curves(ck, METRIC+'/'+vid, field, field, stored, 't', xaxis,
                                     scales[field], cv.PLOT_W/2, params)
        ptext = {'inflation': '2', 'tracking': '2/3', 'coasting': '1'}[vid]
        views.append({'id': vid, 'label': label, 'kind': 'line', 'x': xaxis, 'panels': panels,
                      'curves': curves, 'marks': [],
                      'legend': [['line', 'K', '$K$'], ['line', 'R', '$R$']],
                      'caption': [r'The Kretschmann scalar $K$ and the Ricci scalar $R$ of power law expansion, with $K\propto t^{-4}$ and $R\propto t^{-2}$.'],
                      'ricci_sentences': [],
                      'settings': f'$p={ptext}$, $a_0=t_0=1$.', 'clip': 'Logarithmic vertical scales.',
                      'per_scalar': cv.scalar_clips('Logarithmic scale.', 'Logarithmic scale.')})
    return views
