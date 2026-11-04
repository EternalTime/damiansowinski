"""The first-order collapsing texture, Nötzold (10) and Durrer et al. (0.46).

The displayed tensors belong exactly to the truncated metric.
The texture stress matches its Einstein tensor only through first order.
"""
import sympy as sp

ANGULAR = r"H = 1 - 2\epsilon + \dfrac{2\epsilon t}{r}\arctan\left(\dfrac{r}{t}\right)"
LINE = r"ds^2=-dt^2+dr^2+r^2H\left(d\theta^2+\sin^2\theta\,d\phi^2\right)"


def chart_spec():
    import verify_metrics as vm
    coords = ["t", "r", r"\theta", r"\phi"]
    probe = vm.Reader(coords, [r"\epsilon", ANGULAR], (), held=("H",))
    H = probe.parameters["H"]
    r = probe.symbol["r"]
    return {
        "metric_id": "global_texture",
        "system": {"id": "collapse", "name": "Before Collapse", "coords": coords,
                   "domains": [r"t \in (-\infty, 0)", r"r \in (0, \infty)",
                               r"\theta \in [0, \pi]", r"\phi \in [0, 2\pi)"],
                   "parameters": [r"\epsilon", ANGULAR], "line_element": LINE},
        "chart_line_element": LINE,
        "printer": {"lead": [H, r], "flip": False},
        "check": chart_check,
    }


def chart_check(chart):
    """Compare the linear Einstein coefficient to the flat-space sigma-model stress.

    All comparisons are numerical at rational points, without symbolic determinants
    or simplification of a check expression.
    """
    t, r, theta, _ = chart.symbols
    epsilon = chart.reader.parameters["epsilon"]
    explicit = chart.reader.held
    g0 = sp.diag(-1, 1, r**2, r**2 * sp.sin(theta)**2)
    chi = sp.pi + 2 * sp.atan(t/r)
    dt, dr = sp.diff(chi, t), sp.diff(chi, r)
    inner = sp.zeros(4)
    inner[0, 0], inner[0, 1], inner[1, 0], inner[1, 1] = dt**2, dt*dr, dt*dr, dr**2
    inner[2, 2] = sp.sin(chi)**2
    inner[3, 3] = sp.sin(chi)**2 * sp.sin(theta)**2
    trace = -dt**2 + dr**2 + 2*sp.sin(chi)**2/r**2
    stress = inner - g0*trace/2
    G = chart.geo.einstein_ll()
    points = [(-3, sp.Rational(2, 3)), (-1, 2), (sp.Rational(-1, 4), 3)]
    for a in range(4):
        for b in range(a, 4):
            exact = G[a][b].subs(explicit).doit()
            leading = sp.diff(exact, epsilon).subs(epsilon, 0)
            for time, radius in points:
                at = {t: time, r: radius, theta: sp.Rational(4, 5)}
                if abs(sp.N((leading-stress[a, b]).subs(at), 40)) > sp.Float('1e-28'):
                    raise AssertionError(f"global_texture: first-order Einstein equation fails in {a}{b}")
                if abs(sp.N(exact.subs(epsilon, 0).subs(at), 40)) > sp.Float('1e-28'):
                    raise AssertionError("global_texture: flat limit is not vacuum")
    gamma = chart.geo.christoffel_ull()
    for a in (0, 1):
        for b in (0, 1):
            for c in (0, 1):
                if gamma[a][b][c] != 0:
                    raise AssertionError("global_texture: radial geodesics are not straight")


METRIC = 'global_texture'
PARAMS = {'epsilon': '1/20'}


def embedding_views(ck, src, em):
    """The t=-1 equator, isometrically embedded using the stored metric."""
    sl = em.Slice(src, METRIC, 'collapse', 'r', r'\phi', {'t': -1, **em.EQUATOR}, PARAMS)
    size = 8.0
    piece = em.Piece('equator', 'sheet', sl, 0, 4, 0, 1,
                     (('axis', r'$r=0$'), ('edge', r'$r=4\ell$')),
                     [(x, 'r', None) for x in (1, 2, 3, 4)], size)
    ck.isometry('global texture: pre-collapse equator', piece)
    surface = em.Surface([piece], label=r'$t=-\ell$', time=-1)
    fig = em.figure_of([surface], {'sheet': 'cover'}, size, em.Camera(-70, 30))
    fig.legend('fill', 'cover', 'Equatorial space')
    fig.legend('line', 'r', 'Constant proper radius')
    fig.legend('line', 'meridian', 'Radial lines')
    data = fig.done()
    data['box'] = [-4.5, 4.5, -4.5, 4.5]
    return [em.view('collapse', 'Before collapse', r'$\ell$', [surface], data,
                    settings=r'$t=-\ell$, $\epsilon=1/20$; any length $\ell$ as unit.')]


def embedding_captions():
    return {(METRIC, 'collapse'): [r'Equatorial space with circumference below $2\pi r$ ($t=-\ell$).']}


def slice_marks(slices):
    return slices.one(METRIC, lambda m: slices.along(-1, *m.reach('collapse', 'r')),
                      view_id='collapse', label=r'$t=-\ell$')


def diagram_specs(nr):
    return [nr.Diagram(METRIC, 'collapse', 'radial', 'Before collapse', ('t', 'r'),
                       (0, 4, -4, -.05), r'$r/\ell$', r'$t/\ell$', PARAMS, nr.EQUATOR,
                       curves=(('world', '-t', r'$r=-t$, the self similar light cone'),),
                       input=r'$\epsilon=1/20$; any length $\ell$ as unit.')]


def diagram_captions():
    return {(METRIC, 'collapse', 'radial'): [r'Radial light rays before collapse, with the embedded moment marked at $t=-\ell$.']}


def null_forms():
    return {(METRIC, 'collapse', 'radial'): (lambda t, r: t+r, lambda t, r: t-r, None)}


def conformal_views(ck, src, cf):
    """Compactify the pre-collapse radial region, without a post-collapse choice."""
    import numpy as np
    plane = cf.Plane(src, METRIC, 'collapse', ('t', 'r'), cf.EQUATOR, PARAMS)
    ck.chart('global texture before collapse', plane, cf.mink_pq,
             ck.uniform(-8, -.05), ck.uniform(.05, 8), lambda t, r: (1, 0))
    pi = np.pi
    v = cf.View('collapse', 'Before collapse', [-.32, pi+.32, -pi-.32, .32], 'collapse')
    triangle = [[0, -pi], [0, 0], [pi, 0]]
    v.fill('region', triangle)
    v.fill('cover', triangle)
    v.line('centre', [[[0, -pi], [0, 0]]])
    v.line('chartedge', [[[0, 0], [pi, 0]]])
    v.line('scri', [[[0, -pi], [pi, 0]]])
    times = -np.geomspace(1e-5, 1e5, 600)
    for radius in (.5, 1, 2, 4):
        v.curve('r', *cf.mink_pq(times, np.full_like(times, radius)))
    radii = np.r_[0, np.geomspace(1e-5, 1e5, 600)]
    for time in (-.5, -1, -2, -4):
        v.curve('t', *cf.mink_pq(np.full_like(radii, time), radii))
    v.curve('null', *cf.mink_pq(-radii, radii))
    v.label_xt([0, -pi], r'$i^-$', 't', dy=8)
    v.label_xt([pi, 0], r'$i^0$', 'l', dx=8)
    v.label_xt([pi/2, 0], r'$t=0$', 'b', dy=-7)
    v.label_xt([0, -pi/2], r'$r=0$', 'r', dx=-7)
    v.label_xt([pi/2, -pi/2], r'$\mathscr I^-$', 'tl', dx=7, dy=5)
    v.legend('r', r'$r$ constant')
    v.legend('t', r'$t$ constant')
    v.legend('null', r'$r=-t$')
    v.legend('chartedge', r'$t=0$')
    v.legend('scri', 'Past null infinity')
    moment = cf.slices.moments(METRIC, 'collapse')[0]
    xs = np.linspace(*moment.reach('collapse', 'r'), 81)
    v.slice(moment, [cf.mink_pq(-np.ones_like(xs), xs)])
    v.set(settings=r'Any length $\ell$ as unit; $p=\arctan((t-r)/\ell)$, $q=\arctan((t+r)/\ell)$.',
          restriction='Each point represents a symmetry sphere.')
    return [v]


def conformal_captions():
    return {(METRIC, 'collapse'): [r'Radial causal structure before collapse, with the light cone $r=-t$ ending at the collapse.']}


def curvature_views(ck, src, cv):
    """Read both exact scalars of the truncated metric along the embedded moment."""
    xaxis = cv.axis(r'$r/\ell$', .05, 4)
    scales = {'K': {'kind': 'log', 'from': 1e-4, 'to': 1},
              'R': {'kind': 'symlog', 'from': -.02, 'to': 1, 'linear': .01}}
    panels = [{'id': name, 'name': '$'+name+(r'\ell^4$' if name=='K' else r'\ell^2$'),
               'aspect': 1.6, 'scale': scales[name], 'ticks': cv.value_ticks(scales[name])}
              for name in ('K', 'R')]
    cv.shape_checks(ck, METRIC, [p['aspect'] for p in panels])
    curves = []
    for name, field in [('K', 'kretschmann'), ('R', 'ricci_scalar')]:
        stored = cv.nr_published(src, METRIC, 'collapse', field, ['r', 't'], ['epsilon'])
        held = cv.Held(stored, {'t': -1})
        curves += cv.line_curves(ck, METRIC, name, name, held, 'r', xaxis, scales[name],
                                 cv.PLOT_W/1.6, {'epsilon': .05})
    return [{'id': 'collapse', 'label': 'Before collapse', 'kind': 'line', 'x': xaxis,
             'panels': panels, 'curves': curves, 'marks': [],
             'legend': [['line', 'K', '$K$, Kretschmann scalar'], ['line', 'R', '$R$, Ricci scalar']],
             'caption': [r'The Kretschmann scalar $K$ and the Ricci scalar $R$ of the truncated metric along the embedded moment ($t=-\ell$).'],
             'ricci_sentences': [], 'settings': r'$t=-\ell$, $\epsilon=1/20$; any length $\ell$ as unit.',
             'clip': 'Logarithmic $K$; signed logarithmic $R$.',
             'per_scalar': cv.scalar_clips('Logarithmic scale.', 'Signed logarithmic scale.')}]
