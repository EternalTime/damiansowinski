"""Xanthopoulos's vacuum solution in Gurtug and Sakalli, equations (1)-(3), (55)."""

import mpmath as mp
import sympy as sp

METRIC = 'string_cylindrical_waves'
SYSTEMS = ('prolate', 'cylindrical', 'hyperbolic')


def parameters(system):
    names = [r'\ell', r'\alpha', 'p', r'q = \sqrt{1+p^2}']
    if system == 'cylindrical':
        names += [r'D = \sqrt{(1+\rho^2-t^2)^2+4t^2}',
                  r'\sigma = \sqrt{(D+1+\rho^2-t^2)/2}', r'\tau = t/\sigma']
    elif system == 'hyperbolic':
        names += [r'\tau = \sinh\psi', r'\sigma = \cosh\theta']
    names += [r'\Delta = 1+\tau^2', r'\delta = \sigma^2-1',
              r'X = (1-p\tau)^2+q^2\sigma^2', r'Y = X-2(1-p\tau)',
              r'\omega = 2q(1-p\tau)\delta/(pY)']
    return names


def register_dimensions(vm):
    for system in SYSTEMS:
        key = (METRIC, system)
        coords = {'prolate': (r'\tau', r'\sigma'), 'cylindrical': ('t', r'\rho'),
                  'hyperbolic': (r'\psi', r'\theta')}[system]
        if system != 'prolate':
            vm.HELD[key] = ('X', 'Y', 'omega')
        dims = {name: '1' for name in coords + (r'\phi', 'z')}
        dims.update({r'\ell': 'L', r'\alpha': '1', 'p': '1', 'q': '1',
                     r'\Delta': '1', r'\delta': '1', 'X': '1', 'Y': '1', r'\omega': '1'})
        if system != 'prolate':
            dims.update({r'\tau': '1', r'\sigma': '1'})
        if system == 'cylindrical':
            dims['D'] = '1'
        vm.DIMENSIONS[key] = dims


def spec(system):
    coords = {'prolate': [r'\tau', r'\sigma', r'\phi', 'z'],
              'cylindrical': ['t', r'\rho', r'\phi', 'z'],
              'hyperbolic': [r'\psi', r'\theta', r'\phi', 'z']}[system]
    plane = {'prolate': r'\alpha^2 X(-d\tau^2/\Delta+d\sigma^2/\delta)',
             'cylindrical': r'\alpha^2X(-dt^2+d\rho^2)/(\tau^2+\sigma^2)',
             'hyperbolic': r'\alpha^2X(-d\psi^2+d\theta^2)'}[system]
    radius = r'\rho^2' if system == 'cylindrical' else r'\Delta\delta'
    line = (r'ds^2=\ell^2\left(' + plane + '+' + radius + r' X\,d\phi^2/Y'
            r'+(Y/X)(dz-\omega\,d\phi)^2\right)')
    domain = {'prolate': [r'\tau\in(-\infty,\infty)', r'\sigma\in(1,\infty)'],
              'cylindrical': [r't\in(-\infty,\infty)', r'\rho\in(0,\infty)'],
              'hyperbolic': [r'\psi\in(-\infty,\infty)', r'\theta\in(0,\infty)']}[system]
    domain += [r'\phi\in[0,2\pi)', r'z\in(-\infty,\infty)', r'\ell>0,\quad\alpha>0,\quad p>0']
    names = {'prolate': 'Prolate Coordinates', 'cylindrical': 'Cylindrical',
             'hyperbolic': 'Hyperbolic Coordinates'}
    return {'metric_id': METRIC, 'system': {'id': system, 'name': names[system], 'coords': coords,
            'domains': domain, 'parameters': parameters(system), 'line_element': line},
            'chart_line_element': line, 'printer': {},
            'check': lambda chart: check(chart, system)}


def source_metric(system, first, second, p, alpha, ell):
    """The source metric in each chart, evaluated without a symbolic determinant."""
    q = mp.sqrt(1+p*p)
    if system == 'prolate':
        tau, sigma = first, second
    elif system == 'hyperbolic':
        tau, sigma = mp.sinh(first), mp.cosh(second)
    else:
        disc = mp.sqrt((1+second**2-first**2)**2+4*first**2)
        sigma = mp.sqrt((disc+1+second**2-first**2)/2)
        tau = first/sigma
    delta, capital_delta = sigma*sigma-1, 1+tau*tau
    X = (1-p*tau)**2+q*q*sigma*sigma
    Y = p*p*capital_delta+q*q*delta
    omega = 2*q*delta*(1-p*tau)/(p*Y)
    plane = alpha*alpha*X
    if system == 'prolate':
        tpart, rpart = -plane/capital_delta, plane/delta
    elif system == 'cylindrical':
        tpart, rpart = -plane/(tau*tau+sigma*sigma), plane/(tau*tau+sigma*sigma)
    else:
        tpart, rpart = -plane, plane
    orbit = capital_delta*delta*X/Y
    return ell**2*mp.matrix([[tpart, 0, 0, 0], [0, rpart, 0, 0],
                            [0, 0, orbit+Y*omega*omega/X, -Y*omega/X],
                            [0, 0, -Y*omega/X, Y/X]])


def check(chart, system):
    """Hold the printed metric and the vacuum equations to independent point checks."""
    reader = chart.reader
    g = chart.geo.g
    mp.mp.dps = 45
    for first, second, p, alpha, ell in [(0.7, 1.6, 1, 1, 1), (-0.4, 2.1, 0.8, 1.3, 2)]:
        first, second, p, alpha, ell = map(lambda v: mp.mpf(str(v)), (first, second, p, alpha, ell))
        expected = source_metric(system, first, second, p, alpha, ell)
        at = {chart.symbols[0]: sp.Float(str(first), 45), chart.symbols[1]: sp.Float(str(second), 45),
              chart.symbols[2]: 0, chart.symbols[3]: 0, reader.parameters['p']: sp.Float(str(p), 45),
              reader.parameters['alpha']: sp.Float(str(alpha), 45), reader.parameters['ell']: sp.Float(str(ell), 45)}
        for i in range(4):
            for j in range(4):
                value = g[i, j].subs(reader.held).doit().subs(at).evalf(40)
                if abs(mp.mpf(str(value))-expected[i,j]) > mp.mpf('1e-30'):
                    raise AssertionError(f'{system}: source metric differs at {i},{j}')
        vacuum_at_point(system, first, second, p, alpha, ell)


def vacuum_at_point(system, first, second, p, alpha, ell):
    """Ricci from metric jets at a point, with no symbolic tensor reduction."""
    def metric(x, y):
        return source_metric(system, x, y, p, alpha, ell)

    g = metric(first, second)
    inv = g**-1
    dg = [mp.matrix(4) for _ in range(4)]
    ddg = [[mp.matrix(4) for _ in range(4)] for _ in range(4)]
    for i, j in ((0,0), (1,1), (2,2), (2,3), (3,3)):
        f = lambda x, y: metric(x, y)[i,j]
        for k in range(2):
            orders = tuple(int(n == k) for n in range(2))
            dg[k][i,j] = dg[k][j,i] = mp.diff(f, (first, second), orders)
            for l in range(2):
                orders2 = tuple(int(n == k)+int(n == l) for n in range(2))
                ddg[k][l][i,j] = ddg[k][l][j,i] = mp.diff(f, (first, second), orders2)
    dinv = [-inv*dg[k]*inv for k in range(4)]
    gamma = [[[mp.mpf('0') for _ in range(4)] for _ in range(4)] for _ in range(4)]
    dgamma = [[[[mp.mpf('0') for _ in range(4)] for _ in range(4)] for _ in range(4)] for _ in range(4)]
    for a in range(4):
        for b in range(4):
            for c in range(4):
                for d in range(4):
                    h = dg[b][d,c]+dg[c][d,b]-dg[d][b,c]
                    gamma[a][b][c] += inv[a,d]*h/2
                    for k in range(2):
                        dh = ddg[k][b][d,c]+ddg[k][c][d,b]-ddg[k][d][b,c]
                        dgamma[k][a][b][c] += (dinv[k][a,d]*h+inv[a,d]*dh)/2
    for b in range(4):
        for d in range(4):
            value = sum(dgamma[a][a][b][d]-dgamma[d][a][b][a] for a in range(4))
            value += sum(gamma[a][a][e]*gamma[e][b][d]-gamma[a][d][e]*gamma[e][b][a]
                         for a in range(4) for e in range(4))
            if abs(value) > mp.mpf('1e-30'):
                raise AssertionError(f'{system}: Ricci {b},{d} = {value}')


def register(pc):
    pc.CHARTS[METRIC] = [lambda system=system: spec(system) for system in SYSTEMS]
