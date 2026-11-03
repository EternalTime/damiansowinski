"""Drawings of the cylindrical wave geometry at p=1, alpha=2 and ell=1."""

import numpy as np

METRIC = 'string_cylindrical_waves'
PARAMS = {'p': '1', 'alpha': '2', 'ell': '1'}


def register_null(nr):
    charts = [
        ('prolate', r'\tau', r'\sigma', (.8, 4.2, -2., 2.), 'tau'),
        ('hyperbolic', r'\psi', r'\theta', (-.15, 2.85, -1.5, 1.5), 'psi'),
    ]
    for system, time, radius, box, time_name in charts:
        nr.DIAGRAMS.append(nr.Diagram(METRIC, system, 'radial', f'${time}$ and ${radius}$',
                                     (time, radius), box, f'${radius}$', f'${time}$',
                                     PARAMS, {'phi': '0', 'z': '0'}, tau=time_name))
        nr.CAPTIONS[(METRIC, system, 'radial')] = [
            f'The plane of ${time}$ and ${radius}$ ($\\phi=0$, $z=0$), with radial null geodesics.']
    nr.CLOSED_FORMS[(METRIC, 'prolate', 'radial')] = (
        lambda t, r: np.arcsinh(t)+np.arccosh(np.maximum(r, 1.)),
        lambda t, r: np.arcsinh(t)-np.arccosh(np.maximum(r, 1.)), lambda t, r: r > 1.001)
    nr.CLOSED_FORMS[(METRIC, 'hyperbolic', 'radial')] = (
        lambda t, r: t+r, lambda t, r: t-r, lambda t, r: r > .001)


def register_conformal(cf):
    cf.DRAWN[METRIC] = lambda ck, src: conformal(cf, ck, src)
    for system in ('prolate', 'hyperbolic'):
        cf.CAPTIONS[(METRIC, system)] = [
            r'The radial half plane ($\phi=0$, $z=0$), with null infinity and the string axis.']


def conformal(cf, ck, src):
    def prolate(tau, sigma):
        tau, sigma = np.asarray(tau), np.asarray(sigma)
        return cf.mink_pq(tau*sigma, np.sqrt((1+tau*tau)*(sigma*sigma-1)))

    def hyperbolic(psi, theta):
        return cf.mink_pq(np.sinh(psi)*np.cosh(theta), np.cosh(psi)*np.sinh(theta))

    views = []
    for system, coords, mapping, radial in [
            ('prolate', (r'\tau', r'\sigma'), prolate, (1.15, 1.5, 2., 3.)),
            ('hyperbolic', (r'\psi', r'\theta'), hyperbolic, (.25, .6, 1., 2.))]:
        plane = cf.Plane(src, METRIC, system, coords, {'phi': '0', 'z': '0'}, PARAMS)
        x = ck.uniform(-2, 2)
        y = ck.uniform(1.01, 3) if system == 'prolate' else ck.uniform(.01, 2)
        ck.chart(f'string waves: {system} radial plane', plane, mapping, x, y, lambda a, b: (1, 0))
        v = cf.View(system, 'Radial half plane', [-.35, cf.PI+.35, -cf.PI-.25, cf.PI+.25], system)
        v.fill('region', cf.TRIANGLE)
        v.fill('cover', cf.TRIANGLE)
        cf.grid(v, 'r', lambda r, t: mapping(t, r), radial, cf.S_ALL)
        if system == 'prolate':
            cf.grid(v, 't', lambda t, r: mapping(t, r+1), (-2, -1, 0, 1, 2), cf.S_POS)
        else:
            cf.grid(v, 't', mapping, (-2, -1, 0, 1, 2), cf.S_POS)
        cf.triangle_edges(v, centre=f'${coords[1]}={1 if system == "prolate" else 0}$', centre_class='surface')
        v.legend('cover', 'the radial half plane')
        v.legend('r', f'${coords[1]}$ constant')
        v.legend('t', f'${coords[0]}$ constant')
        v.legend('surface', 'the string axis')
        v.set(settings=r'$p=1$, $\alpha=2$, $\ell=1$.',
              restriction=r'The totally geodesic half plane of fixed $\phi$ and $z$.')
        for moment in cf.slices.moments(METRIC):
            radius = np.linspace(*moment.reach('hyperbolic', r'\theta'), 100)
            if system == 'prolate':
                radius = np.cosh(radius)
            v.slice(moment, [mapping(np.zeros_like(radius), radius)])
        views.append(v)
    return views


def register_embedding(eb):
    eb.DRAWN[METRIC] = lambda ck, src: embedding(eb, ck, src)
    eb.CAPTIONS[(METRIC, 'moment')] = [
        r'The moment $\psi=0$ on the plane $z=0$, from the string outward.']


def embedding(eb, ck, src):
    sl = eb.Slice(src, METRIC, 'hyperbolic', r'\theta', r'\phi', {'psi': 0, 'z': 0}, PARAMS)
    original_rho, original_defect = sl._rho, sl._defect
    sl._rho = lambda x: np.where(np.asarray(x) == 0, 0., original_rho(x))
    sl._defect = lambda x: np.where(np.asarray(x) == 0, 9., original_defect(x))
    points = np.linspace(.001, 2., 500)
    ck.add('String waves: Euclidean embedding across the slice',
           max(0., float(np.max(-sl.defect_at(points)))), 1e-10)
    piece = eb.Piece('exterior', 'sheet', sl, 0., 2., 0., 1,
                     (('axis', 'the string axis'), ('edge', r'$\theta=2$')),
                     [(1., 'r', None), (2., 'r', None)], 12.)
    ck.isometry('String waves: moment psi=0', piece)
    surface = eb.Surface([piece])
    fig = eb.figure_of([surface], {'sheet': 'cover'}, 12., eb.Camera(-90, 30))
    fig.legend('fill', 'cover', 'the moment of the wave geometry')
    fig.legend('line', 'r', r'$\theta$ constant')
    fig.legend('line', 'meridian', r'$\phi$ constant')
    return [eb.view('moment', 'Constant time', r'$\ell$', [surface], fig.done(),
                    settings=r'$p=1$, $\alpha=2$, $\ell=1$.')]


def register_slices(sl):
    sl.FLAT[(METRIC, 'hyperbolic', 'radial')] = lambda: sl.one(
        METRIC, lambda m: sl.along(0., *m.reach('hyperbolic', r'\theta')), label=r'$\psi=0$')
    sl.FLAT[(METRIC, 'prolate', 'radial')] = lambda: sl.one(
        METRIC, lambda m: sl.along(0., *np.cosh(m.reach('hyperbolic', r'\theta'))), label=r'$\tau=0$')
    sl.FLAT_METRICS.add(METRIC)


def register_curvature(cv):
    cv.DRAWN[METRIC] = lambda ck, src: curvature(cv, ck, src)


def curvature(cv, ck, src):
    """Read the stored Kretschmann scalar on the prolate coordinate plane."""
    x, y = cv.axis(r'$\sigma$', 1., 4.), cv.axis(r'$\tau$', -1.5, 1.5)
    nx = ny = 96
    stored = cv.Published(src, METRIC, 'prolate', 'kretschmann', ['sigma', 'tau'], ['p', 'alpha', 'ell'])
    _, _, values = cv.cells(stored, x, y, nx, ny, [1., 2., 1.])
    ck('String waves: the stored curvature is finite on the displayed plane', np.all(np.isfinite(values)))
    low, high = float(np.min(values)), float(np.max(values))
    width = max(high-low, .01)
    scale = {'kind': 'linear', 'from': low-.05*width, 'to': high+.05*width}
    return [{'id': 'K', 'label': '$K$', 'kind': 'map', 'x': x, 'y': y,
             'colour': {'name': r'$K\ell^4$', 'scale': scale, 'ticks': cv.value_ticks(scale)},
             'map': {'metric': METRIC, 'system': 'prolate', 'scalar': 'kretschmann',
                     'x': r'\sigma', 'y': r'\tau', 'params': {'p': 1., 'alpha': 2., 'ell': 1.},
                     'nx': nx, 'ny': ny, 'cells': [[cv.rounded(v) for v in row] for row in values]},
             'marks': [], 'legend': [],
             'caption': [r'The Kretschmann scalar $K$ in $\sigma$ and $\tau$.'],
             'settings': r'$p=1$, $\alpha=2$, $\ell=1$.', 'clip': 'Linear colour scale.'}]
