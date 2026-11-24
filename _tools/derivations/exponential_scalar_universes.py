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


# The embedding's movies, one a member, over the same span of cosmic time, a frame every STEP of it,
# with the key moments the other drawings mark. Each ring of free particles at rest in the chart has
# the comoving radius RING, so that it meets the member's horizon at t = t_0, and each moment is the
# disc of proper radius TOP about the ring's centre.
MOVIE_KEYS = (("1/2", r"$t=0.5\,t_0$"), ("1", "$t=t_0$"), ("3/2", r"$t=1.5\,t_0$"))
STEP = 0.025
RING = {"inflation": 1, "tracking": 3, "coasting": 1}
TOP = {"inflation": 2.5, "tracking": 5.0, "coasting": 1.75}
HORIZON = {"inflation": "event", "tracking": "particle"}
PTEXT = {"inflation": "2", "tracking": "2/3", "coasting": "1"}


def horizon_radius(src, p):
    """The comoving radius of the horizon of the observer at r = 0 as a function of t, read from the
    published metric of the cosmic Cartesian chart at a_0 = t_0 = 1: the integral of dt/sqrt(g_xx)
    from t to infinity, the event horizon, for p > 1, and from 0 to t, the particle horizon, for p < 1.
    Returns the sympy expression in a positive t, or None where both integrals diverge, at p = 1."""
    import sympy as sp
    import null_rays as nr

    src.note(METRIC, "cosmic_cartesian", ["metric_components"])
    _, entry, R = nr.load(METRIC, "cosmic_cartesian")
    t, s = sp.Symbol("t", positive=True), sp.Symbol("s", positive=True)
    subs = {R.c: 1, R.parameters["p"]: sp.nsimplify(p), R.parameters["a_0"]: 1, R.parameters["t_0"]: 1}
    g = nr.published_matrix(R, entry, "metric_components").subs(R.held).doit().subs(subs)
    rate = sp.powsimp(sp.simplify(1 / sp.sqrt(g[1, 1].subs(R.symbol["t"], s))), force=True)
    future = sp.integrate(rate, (s, t, sp.oo))
    past = sp.integrate(rate, (s, 0, t))
    if p > 1 and future.is_finite and not past.is_finite:
        return sp.simplify(future)
    if p < 1 and past.is_finite and not future.is_finite:
        return sp.simplify(past)
    if p == 1 and not future.is_finite and not past.is_finite:
        return None
    raise AssertionError(f"{METRIC}: the horizon integrals of p = {p} are {future} and {past}")


def embedding_views(ck, src, em):
    """For each member, the plane z = 0 of the cosmic Cartesian chart at the moments t = t_0/2 to
    3t_0/2, a frame every t_0/40, at a_0 = t_0 = 1. Each moment is flat, a^2 (dx^2 + dy^2) with the
    chart's own a = t^p, and a ring of free particles at rest in the chart on x^2 + y^2 = RING^2,
    which the published Christoffel symbols keep at rest, is the circle of proper radius RING t^p:
    growing faster and faster at p = 2, more and more slowly at p = 2/3 and steadily at p = 1. Beside
    it, where the member has one, is the horizon of the observer at the centre, read from the
    published metric: at p = 2 the event horizon, of proper radius t/(p - 1) = t, and at p = 2/3 the
    particle horizon, of proper radius t/(1 - p) = 3t, each of which the ring meets at t = t_0, the
    ring of comoving radius 1 inside the event horizon before t_0 and outside it after, the ring of
    comoving radius 3 outside the particle horizon before t_0 and inside it after. At p = 1 both
    integrals of dt/a diverge and no horizon is drawn."""
    import math
    import numpy as np
    import sympy as sp
    import null_rays as nr

    em.at_rest(ck, src, METRIC, "cosmic_cartesian")
    _, _, reader = nr.load(METRIC, "cosmic_cartesian")
    (scale,) = (str(rule) for rule in reader.held.values())
    keys = [float(sp.Rational(k)) for k, _ in MOVIE_KEYS]
    times, at = em.movie_values(keys, STEP)
    alpha = np.linspace(0, 2 * math.pi, 361)
    views = []
    for vid, label, p in EXAMPLES:
        name, R, top = f"{label} ($p={PTEXT[vid]}$)", RING[vid], TOP[vid]
        size = 2 * top
        params = {"p": PTEXT[vid], "a_0": 1, "t_0": 1}
        closed = horizon_radius(src, p)
        r_h = None if closed is None else sp.lambdify(sp.Symbol("t", positive=True), closed, "math")
        if vid == "inflation":
            ck.exact(f"{name}: the event horizon's comoving radius is 1/t", sp.simplify(closed - 1 / sp.Symbol("t", positive=True)) == 0)
        elif vid == "tracking":
            ck.exact(f"{name}: the particle horizon's comoving radius is 3t^(1/3)",
                     sp.simplify(closed - 3 * sp.cbrt(sp.Symbol("t", positive=True))) == 0)
        else:
            ck.exact(f"{name}: both integrals of dt/a diverge, so the observer at the centre has no horizon", closed is None)
        frames, radii, horizons = [], [], []
        for i, t in enumerate(times):
            key = at.index(i) if i in at else None
            fixed = {"t": MOVIE_KEYS[key][0] if key is not None else repr(round(t, 9)), "z": 0}
            moment = MOVIE_KEYS[key][1] if key is not None else f"$t={t:g}\\,t_0$"
            where = f"{name}, t = {t:g}"
            sl = em.FlatPlane(src, METRIC, "cosmic_cartesian", "x", "y", fixed, params, {"a": scale})
            ck.add(f"{where}: the plane is a^2(dx^2 + dy^2) with a = t^p", float(np.max(np.abs(sl.scale - t ** p))), 1e-12)
            marks = [] if r_h is None else [(r_h(t), "horizon", None)]
            plane = em.disc(sl, "plane", top, "the centre of the ring", "the plane continues, flat, to infinity", marks, size)
            ck.isometry(where, plane)
            ring, dots = em.particles(ck, where, sl, plane, R * np.cos(alpha), R * np.sin(alpha))
            radius = np.hypot(ring.points[:, 0], ring.points[:, 1])
            ck.add(f"{where}: the ring is the circle of proper radius {R} t^p", float(np.max(np.abs(radius - R * t ** p))), 1e-12)
            radii.append(float(radius[0]))
            if r_h is not None:
                rho, _ = plane.at(r_h(t))
                ck.add(f"{where}: the {HORIZON[vid]} horizon's proper radius is t/|p - 1|", abs(float(rho) - t / abs(p - 1)), 1e-12)
                ck.exact(f"{where}: the horizon lies on the disc", t / abs(p - 1) < top)
                horizons.append(float(rho))
            frames.append(em.Surface([plane], label=moment, time=t, curves=[ring], dots=dots))
        steps = np.diff(radii)
        growth = np.diff(steps)
        if vid == "inflation":
            ck.exact(f"{name}: the ring grows faster and faster", bool(np.all(steps > 0) and np.all(growth > 0)))
        elif vid == "tracking":
            ck.exact(f"{name}: the ring grows more and more slowly", bool(np.all(steps > 0) and np.all(growth < 0)))
        else:
            ck.add(f"{name}: the ring grows steadily", float(np.max(np.abs(growth))), 1e-12)
        if horizons:
            k = at[1]
            ck.add(f"{name}: the ring meets the {HORIZON[vid]} horizon at t = t_0", abs(radii[k] - horizons[k]), 1e-12)
            inside_before = vid == "inflation"
            ck.exact(f"{name}: the ring lies {'inside' if inside_before else 'outside'} the horizon before t_0 and "
                     f"{'outside' if inside_before else 'inside'} it after",
                     all((r < h) == inside_before for r, h in zip(radii[:k], horizons[:k]))
                     and all((r > h) == inside_before for r, h in zip(radii[k + 1:], horizons[k + 1:])))
        fig, played = em.ring_movie(frames, size, "$t$")
        fig.legend("fill", "cover", "the plane $z=0$ at each moment, flat")
        fig.legend("line", "particles", "free particles at rest in the chart, twelve marked")
        if vid in HORIZON:
            fig.legend("line", "horizon", f"the {HORIZON[vid]} horizon of the observer at the centre")
        fig.legend("line", "meridian", "radii, every $30°$")
        views.append(em.view(vid, label, "$t_0$", [frames[i] for i in at], fig.done(), movie=played,
                             settings=f"$p={PTEXT[vid]}$, $a_0=t_0=1$."))
    return views


def embedding_captions():
    """A ring of free particles on the plane z = 0 of each member, with the horizon where it has one."""
    span = r"($t_0/2\le t\le3t_0/2$)"
    return {
        (METRIC, "inflation"): [
            rf"The plane $z=0$ through a ring of free particles {span}, flat at every moment, the ring growing as $t^2$ "
            r"and crossing the event horizon, of proper radius $t$, at $t=t_0$."],
        (METRIC, "tracking"): [
            rf"The plane $z=0$ through a ring of free particles {span}, flat at every moment, the ring growing as "
            r"$t^{2/3}$ and crossed by the particle horizon, of proper radius $3t$, at $t=t_0$."],
        (METRIC, "coasting"): [
            rf"The plane $z=0$ through a ring of free particles {span}, flat at every moment, the ring growing as $t$."],
    }


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
                            (clock, radial), box, '$x/t_0$' if spherical else f'${radial}/t_0$', f'${clock}/t_0$', params, fixed,
                            tau='t' if cosmic else 'eta', mirror=spherical,
                            families=sideways if not spherical else ('ingoing', 'outgoing'),
                            curves=curve, cones=(5, 5), step=0.00025 if spherical else 0.0025))
    return rows


def diagram_captions():
    """Compact descriptions of the null planes."""
    return {(METRIC, system, vid): [MATTER_TEXT[vid]+
        (r", radial null geodesics along the signed diameter $x=\pm r$." if system == 'cosmic_spherical' else
         r", null geodesics on $y=z=0$.")]
        for system, vid, _ in PLANES}


def null_forms():
    """First integrals of the null slopes, checked against traced rays."""
    return {(METRIC, system, vid):
            ((lambda t, r, p=p: eta_of(t, p)+r), (lambda t, r, p=p: eta_of(t, p)-r), None)
            if system.startswith('cosmic') else ((lambda e, r: e+r), (lambda e, r: e-r), None)
            for system, vid, p in PLANES}


def slice_marks(slices, system, vid):
    """Each moment t of the member's movie in the drawing's chart at the same parameters: the moment
    is the whole plane z = 0, so its line of constant time runs across the drawing, t in the cosmic
    charts and eta(t) in the conformal ones."""
    p = dict((x[0], x[2]) for x in EXAMPLES)[vid]
    lo = 0.0 if system == 'cosmic_spherical' else -slices.BIG
    def line(moment):
        clock = moment.time if system.startswith('cosmic') else float(eta_of(moment.time, p))
        return slices.along(clock, lo, slices.BIG)
    return slices.one(METRIC, line, view_id=vid)


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
        # Each moment of the embedding's movie is the whole plane z = 0, the curve of its eta(t).
        for moment in cf.slices.moments(METRIC, vid):
            v.slice(moment, [cf.mink_pq(np.full_like(rr, float(eta_of(moment.time, p))), rr)])
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


def dust_trace(plot, seed, family, direction):
    """Integrate the published null slope in s=t^(1/3), regular at t=0.

    For the dust example the metric makes ds/dr constant. Its value is read
    from the published null direction, independently of the closed integrals
    used by the ray checker. The cubic maps back to cosmic time only after
    integration, so no positive-time interval needs to be excluded.
    """
    import math
    import numpy as np

    spec = plot.c.spec
    if spec.params != {'a_0': 1, 't_0': 1, 'p': 2/3} or spec.box != (0, 2.4, 0, 4.8):
        raise ValueError('the regular dust integrator requires the declared dust drawing units')
    time, radius = map(float, plot.to_chart(plot.from_unit(seed)))
    if time <= 0:
        return np.array([seed], float)
    root = np.cbrt(time)
    times = np.array([time, 0.1, 1.0, 4.8])
    radii = np.full_like(times, radius)
    tangent = plot.c.null_dirs(times, radii)[family]
    rates = tangent[:, 0]/(3*np.cbrt(times)**2*tangent[:, 1])
    rate = float(rates[0])
    if not np.all(np.isfinite(rates)) or np.max(np.abs(rates-rate)) > 1e-12:
        raise AssertionError('the published dust metric does not have constant ds/dr')
    # The two limits in s, intersected with the drawn radius interval.
    rlimits = sorted((radius-root/rate, radius+(np.cbrt(4.8)-root)/rate))
    low, high = max(0.0, rlimits[0]), min(2.4, rlimits[1])
    end = high if direction*tangent[0, 1] > 0 else low
    # A bound on normalized dt/dr fixes the maximum distance between samples.
    slope = 3*abs(rate)*np.cbrt(4.8)**2*2.4/4.8
    count = max(2, math.ceil(abs(end-radius)/2.4*math.sqrt(1+slope*slope)/spec.step)+1)
    r = np.linspace(radius, end, count)
    s = np.maximum(root+rate*(r-radius), 0.0)
    t = s*s*s
    return (np.column_stack((r,t))-plot.lo)/plot.span
