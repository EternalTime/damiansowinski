"""Drawings of the six pole constant modulus limit of the modular string geometry."""
import numpy as np

METRIC = 'stringy_cosmic_strings'
VIEW = 'six_string_cone'
SETTING = r'Six coincident strings in the constant shape limit ($\ell=1$).'


def register_embedding(eb):
    eb.DRAWN[METRIC] = lambda ck, src: embedding(eb, ck, src)
    eb.CAPTIONS[(METRIC, VIEW)] = [
        r'The transverse cone at $t=w=0$ in the constant shape limit of six coincident strings.']


def embedding(eb, ck, src):
    """The annular slice of H=ell/r: its circumference is half that of a flat plane."""
    sl = eb.Slice(src, METRIC, 'polar', 'r', r'\phi', {'t': '0', 'w': '0'},
                  functions={'H': '1/r'})
    size = 2 * float(sl.rho_at(6))
    piece = eb.Piece('cone', 'sheet', sl, .25, 6, 0., 1,
                     (('edge', 'the coincident core lies at $r=0$'),
                      ('edge', 'the cone continues outward')),
                     [(.25, 'chartedge', None), (1, 'r', None), (4, 'r', None), (6, 'chartedge', None)], size)
    ck.isometry('Six string cone', piece)
    surface = eb.Surface([piece])
    fig = eb.figure_of([surface], {'sheet': 'cover'}, size)
    fig.legend('fill', 'cover', 'the transverse cone')
    fig.legend('line', 'chartedge', r'$r=\ell/4$ and $r=6\ell$')
    fig.legend('line', 'r', '$r$ constant')
    fig.legend('line', 'meridian', r'$\phi$ constant')
    return [eb.view(VIEW, 'Six string cone', r'$\ell$', [surface], fig.done(), settings=SETTING,
                    input=r'The constant shape limit has $H=\ell/r$.',
                    stops=[r'The slice begins at $r=\ell/4$, outside the coincident core at $r=0$.'])]


def register_slices(sl):
    for system in ('cartesian', 'polar'):
        sl.FLAT[(METRIC, system, 'radial')] = lambda: sl.one(
            METRIC, lambda moment: sl.along(0., *moment.reach('polar', 'r')),
            label=r'$t=0$', view_id=VIEW)
    sl.FLAT_METRICS.add(METRIC)


def register_null(nr):
    for system, coordinate, fixed, function in (
        ('cartesian', 'x', {'w': '0', 'y': '0'}, '1/sqrt(x**2+y**2)'),
        ('polar', 'r', {'w': '0', 'phi': '0'}, '1/r'),
    ):
        nr.DIAGRAMS.append(nr.Diagram(
            METRIC, system, 'radial', f'${coordinate}$ and $t$', ('t', coordinate),
            (.25, 6, -3, 3), f'${coordinate}/\\ell$', r'$ct/\ell$', {}, fixed,
            functions={'H': function}, input=r'Six coincident strings with constant torus shape ($\ell=1$).'))
        fixed_label = r'y=w=0' if system == 'cartesian' else r'\phi=w=0'
        nr.CAPTIONS[(METRIC, system, 'radial')] = [
            f'The plane of $t$ and ${coordinate}$ (${fixed_label}$) in the six string cone.']
        nr.CLOSED_FORMS[(METRIC, system, 'radial')] = (
            lambda t, r: t+2*np.sqrt(r), lambda t, r: t-2*np.sqrt(r), lambda t, r: r > 0)


def register_conformal(cf):
    cf.DRAWN[METRIC] = lambda ck, src: conformal(cf, ck, src)
    cf.CAPTIONS[(METRIC, VIEW)] = [
        r'The radial half plane brought to a finite diagram, with the coincident core and infinity at its boundary ($\phi=w=0$).']


def conformal(cf, ck, src):
    """The cone's radial half plane, with its coincident core at its edge."""
    def pq(t, r):
        return cf.mink_pq(t, 2*np.sqrt(np.asarray(r, dtype=float)))
    plane = cf.Plane(src, METRIC, 'polar', ('t', 'r'), {'phi': '0', 'w': '0'},
                     functions={'H': '1/r'})
    ck.chart('Six string cone radial half plane', plane, pq,
             ck.uniform(-20, 20), ck.uniform(.01, 20), lambda t, r: (1, 0))
    v = cf.View(VIEW, 'Six string cone', [-.5, cf.PI+.35, -cf.PI-.25, cf.PI+.25], 'polar')
    v.fill('region', cf.TRIANGLE)
    v.fill('cover', cf.TRIANGLE)
    cf.grid(v, 'r', lambda r, t: pq(t, r), (.25, 1, 4), cf.S_ALL)
    cf.grid(v, 't', pq, (-4, -2, -1, 0, 1, 2, 4), cf.S_POS+1e-3)
    cf.triangle_edges(v, centre=r'$r=0$', centre_class='surface')
    v.legend('cover', 'the radial half plane')
    v.legend('r', '$r$ constant')
    v.legend('t', '$t$ constant')
    v.legend('surface', 'the coincident core, $r=0$')
    moment = cf.slices.moments(METRIC, VIEW)[0]
    r = np.linspace(*moment.reach('polar', 'r'), 801)
    v.slice(moment, [pq(0*r, r)])
    v.set(restriction=r'The radial half plane of the constant shape limit.', settings=SETTING)
    return [v]


def register_curvature(cv):
    cv.DRAWN[METRIC] = lambda ck, src: curvature(cv, ck, src)


def curvature(cv, ck, src):
    """The stored scalars on the local solution tau=z/ell and H=y/ell."""
    functions = {'H': 'y'}
    axis = cv.axis(r'$y/\ell$', .5, 4.)
    scales = {'K': {'kind': 'log', 'from': 1e-4, 'to': 100.},
              'R': {'kind': 'log', 'from': 1e-2, 'to': 10.}}
    panels, curves = [], []
    for key, field, label in (('K', 'kretschmann', r'$K\ell^4$'),
                              ('R', 'ricci_scalar', r'$R\ell^2$')):
        scalar = cv.Published(src, METRIC, 'cartesian', field, ['y'], [], functions=functions)
        panels.append({'id': key, 'name': label, 'aspect': 1., 'scale': scales[key],
                       'ticks': cv.value_ticks(scales[key])})
        curves += cv.line_curves(ck, METRIC, key, key, scalar, 'y', axis, scales[key], cv.PLOT_W, {})
    cv.shape_checks(ck, METRIC, [1., 1.])
    return [{'id': 'local', 'label': 'Local solution', 'kind': 'line', 'x': axis,
             'panels': panels, 'curves': curves, 'marks': [],
             'legend': [['line', 'K', r'Kretschmann scalar $K$'], ['line', 'R', r'Ricci scalar $R$']],
             'caption': [r'The Kretschmann scalar $K$ and the Ricci scalar $R$ of the local solution ($\tau=z/\ell$, $H=y/\ell$).'],
             'ricci_sentences': [],
             'per_scalar': cv.scalar_clips('The logarithmic scale runs from $10^{-4}$ to $10^2$.',
                                           'The logarithmic scale runs from $10^{-2}$ to $10$.'),
             'settings': r'The unit of length is $\ell=1$.',
             'clip': 'Both panels use logarithmic scales.'}]
