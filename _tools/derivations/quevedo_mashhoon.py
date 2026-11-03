"""The four-parameter vacuum field in prolate and spheroidal coordinates.

The functions are those of Kurmanov et al., arXiv:2503.18028, equations (1)-(11).
The regular form of alpha extends the expression in that paper to zero rotation.
"""
import mpmath as mp

METRIC = 'quevedo_mashhoon'
SYSTEMS = ('prolate_spheroidal',)


def potentials(x, y, mass, spin, deformation, delta, backend=None):
    """Return f, omega, exp(2 gamma) on the real exterior chart."""
    backend = backend or {"sqrt": mp.sqrt, "log": mp.log, "exp": mp.exp}
    sqrt, log, exp = backend["sqrt"], backend["log"], backend["exp"]
    sigma = sqrt(mass * mass - spin * spin)
    alpha = -spin / (mass + sigma)
    p = (3 * y * y - 1) / 2
    Q0 = log((x + 1) / (x - 1)) / 2
    Q1 = x * Q0 - 1
    Q2 = (3 * x * x - 1) * Q0 / 2 - 3 * x / 2
    Q2p = 3 * x * Q0 - (3 * x * x - 1) / (2 * (x * x - 1)) - 3 / 2
    def d(sign):
        return (log((x + sign * y) ** 2 / (x * x - 1)) / 2
                + 3 / 2 * (1 - y * y - sign * x * y)
                + 3 / 4 * (x * (1 - y * y) - sign * y * (x * x - 1))
                * log((x - 1) / (x + 1)))
    lp = alpha * (x*x-1)**(1-delta) * (x+y)**(2*delta-2) * exp(2*deformation*delta*d(1))
    lm = alpha * (x*x-1)**(1-delta) * (x-y)**(2*delta-2) * exp(2*deformation*delta*d(-1))
    def ab(sign):
        power = (x + sign)**(delta - 1)
        return (power * (x*(1-lp*lm) + sign*(1+lp*lm)),
                power * (y*(lp+lm) - sign*(lp-lm)))
    ap, bp = ab(1)
    am, bm = ab(-1)
    R, L = ap*am + bp*bm, ap*ap + bp*bp
    H = (x+1)**(delta-1) * (x*(1-y*y)*(lp+lm)*ap + y*(x*x-1)*(1-lp*lm)*bp)
    chi = ((1+deformation)**2 * log((x*x-1)/(x*x-y*y))/2
           + 2*deformation*(1-p)*Q1
           + deformation**2*(1-p)*((1+p)*(Q1*Q1-Q2*Q2)
               + (x*x-1)*(2*Q2*Q2-3*x*Q1*Q2+3*Q0*Q2-Q2p)/2))
    e = exp(2*deformation*delta*p*Q2)
    return R/L/e, -2*spin - 2*sigma*H/R*e, ((1+mass/sigma)**2/4 * R/(x*x-1)**delta
                                                    * exp(2*delta*delta*chi))


def chart_spec(system):
    """Print derivatives of the source functions without expanding their definitions."""
    import sympy as sp
    import verify_metrics as vm
    import json
    from pathlib import Path

    data = json.loads((Path(__file__).resolve().parents[2] / 'MFS/assets/data/metrics'
                       / f'{METRIC}.json').read_text())
    saved = next(c for c in data['coordinates'] if c['id'] == system)
    coords = saved['coords']
    args = '(x,y)' if system == 'prolate_spheroidal' else '(r,\\theta)'
    parameters = ['\\sigma', 'M', 'a', 'f = f' + args,
                  '\\omega = \\omega' + args, 'h = h' + args]
    probe = vm.Reader(coords, parameters, ())
    u, v = probe.symbol[coords[1]], probe.symbol[coords[2]]
    funcs = [probe.parameters[n] for n in ('f', 'omega', 'h')]
    lead = list(funcs)
    for fn in funcs:
        lead.extend([sp.Derivative(fn, u), sp.Derivative(fn, v),
                     sp.Derivative(fn, (u, 2)), sp.Derivative(fn, u, v),
                     sp.Derivative(fn, (v, 2))])
    lead.extend([u, v, probe.parameters['sigma'], probe.parameters['M'], probe.parameters['a']])
    spec = {'metric_id': METRIC,
            'system': {'id': system, 'name': saved['name'], 'coords': coords,
                       'domains': saved['domains'], 'parameters': [p['symbol'] for p in saved['parameters']],
                       'line_element': saved['line_element']},
            'chart_parameters': parameters, 'chart_line_element': saved['line_element'],
            'printer': {'lead': lead, 'factors': lead},
            'pretty': lambda value: sp.powsimp(sp.factor(sp.sympify(value)), combine='exp')}

    def check(chart):
        check_chart(chart, system)
        value, text = compact_kretschmann(chart)
        # The contraction is exact, kept as products of six-index-pair blocks
        # so that no massive expanded numerator has to be simplified.
        chart.geo.kretschmann = lambda: value
        spec['kretschmann'] = text

    spec['check'] = check
    spec['scalar_check'] = check_scalar_at_points
    return spec


def compact_kretschmann(chart):
    """The Riemann contraction as 4 tr((G R)^2) on antisymmetric pairs."""
    import sympy as sp
    from itertools import combinations
    pairs = list(combinations(range(4), 2))
    lower, inverse = chart.geo.riemann_llll(), chart.geo.ginv
    bivector = [[inverse[a, c] * inverse[b, d] - inverse[a, d] * inverse[b, c]
                 for c, d in pairs] for a, b in pairs]
    components = [[lower[a][b][c][d] for c, d in pairs] for a, b in pairs]
    blocks, rendered = {}, {}
    for i in range(6):
        for j in range(6):
            terms = [(bivector[i][k], components[k][j]) for k in range(6)
                     if bivector[i][k] != 0 and components[k][j] != 0]
            if terms:
                blocks[i, j] = sp.Add(*(a*b for a, b in terms), evaluate=False)
                rendered[i, j] = ' + '.join('\\left(' + chart.printer(a) + '\\right)'
                                           '\\left(' + chart.printer(b) + '\\right)' for a, b in terms)
    squares = [(blocks[i, i], rendered[i, i]) for i in range(6) if (i, i) in blocks]
    cross = [(blocks[i, j], blocks[j, i], rendered[i, j], rendered[j, i])
             for i in range(6) for j in range(i + 1, 6) if (i, j) in blocks and (j, i) in blocks]
    value = 4 * (sp.Add(*(a*a for a, _ in squares), evaluate=False)
                 + 2*sp.Add(*(a*b for a, b, _, _ in cross), evaluate=False))
    square_text = ' + '.join('\\left(' + text + '\\right)^2' for _, text in squares) or '0'
    cross_text = ' + '.join('\\left(' + a + '\\right)\\left(' + b + '\\right)'
                            for _, _, a, b in cross) or '0'
    return value, '4\\left(' + square_text + ' + 2\\left(' + cross_text + '\\right)\\right)'


def check_scalar_at_points(chart, field, difference):
    """Compare the compact printed contraction at two numerical points."""
    import random
    import sympy as sp
    from sympy.core.function import AppliedUndef
    if field != 'kretschmann':
        raise AssertionError(field)
    atoms = difference.atoms(sp.Derivative) | difference.atoms(AppliedUndef) | difference.free_symbols
    rng = random.Random(1985)
    for _ in range(2):
        values = {atom: sp.Rational(rng.randint(20, 80), 17)
                  for atom in sorted(atoms, key=sp.default_sort_key)}
        values[chart.symbols[1]] = sp.Rational(13, 5)
        values[chart.symbols[2]] = sp.Rational(2, 5)
        error = complex(difference.xreplace(values).evalf(40))
        if abs(error) > 1e-20:
            raise AssertionError(f'{METRIC}: compact Kretschmann contraction differs by {error}')


def check_chart(chart, system):
    """Compare the displayed metric to the sourced form at nontrivial points."""
    import sympy as sp
    with mp.workdps(40):
        for mass, spin, deformation, delta, x, y in [
                (2, mp.mpf('.7'), mp.mpf('.4'), mp.mpf('1.2'), mp.mpf('2.7'), mp.mpf('.3')),
                (1, mp.mpf('-.4'), mp.mpf('-.3'), mp.mpf('.8'), mp.mpf('3.2'), mp.mpf('-.2'))]:
            sigma = mp.sqrt(mass*mass-spin*spin)
            f, omega, eg = potentials(x, y, mass, spin, deformation, delta)
            if system == 'prolate_spheroidal':
                point = (x, y)
                xx = sigma*sigma * eg*(x*x-y*y)/(f*(x*x-1))
                yy = sigma*sigma * eg*(x*x-y*y)/(f*(1-y*y))
            else:
                r, theta = mass+sigma*x, mp.acos(y)
                point = (r, theta)
                xx = eg*(x*x-y*y)/(f*(x*x-1))
                yy = sigma*sigma*eg*(x*x-y*y)/f
            expected = [[-f,0,0,f*omega], [0,xx,0,0], [0,0,yy,0],
                        [f*omega,0,0,sigma*sigma*(x*x-1)*(1-y*y)/f-f*omega*omega]]
            values = {chart.symbols[1]: sp.Float(str(point[0]), 40), chart.symbols[2]: sp.Float(str(point[1]), 40)}
            for name, val in [('sigma',sigma), ('M',mass), ('a',spin), ('f',f), ('omega',omega), ('h',eg)]:
                values[chart.reader.parameters[name]] = sp.Float(str(val), 40)
            for i in range(4):
                for j in range(4):
                    got = mp.mpf(str(sp.N(chart.geo.g[i,j].subs(values, simultaneous=True), 35)))
                    if abs(got-expected[i][j]) > mp.mpf('1e-28') * max(1,abs(expected[i][j])):
                        raise AssertionError(f'{METRIC}/{system}: metric component {i},{j} disagrees at {point}')


# One common source of the functions and their derivatives for the figures.
# The fixed parameters are dimensionless in units of the mass parameter.
FIGURE_PARAMETERS = {'M': 1, 'a': '2/5', 'q': '2/5', 'delta': 1}


def _figure_values(u, v, system, digits=None):
    """Evaluate both derivatives of the three displayed functions as second order jets."""
    import numpy as np
    import black_saturn as bs
    if digits is None:
        u, v = np.asarray(u, dtype=float), np.asarray(v, dtype=float)
        one, zero = 1.0, 0.0
        mass, spin, deform, delta = 1.0, .4, .4, 1.0
        backend = bs.NUMPY
    else:
        one, zero = mp.mpf(1), mp.mpf(0)
        u, v = mp.mpf(u), mp.mpf(v)
        mass, spin, deform, delta = one, mp.mpf(2)/5, mp.mpf(2)/5, one
        backend = bs.MPMATH
    first, second = bs.Jet(u, one, zero), bs.Jet(v, zero, one)
    if system == 'spherical':
        sigma = backend['sqrt'](mass*mass-spin*spin)
        x, y = (first-mass)/sigma, backend['cos'](second)
    elif system == 'prolate_spheroidal':
        x, y = first, second
    else:
        raise ValueError(system)
    result = potentials(x, y, mass, spin, deform, delta, backend)
    return {name: dict(zip(bs.ORDERS, value.parts())) for name, value in zip(('f', 'omega', 'h'), result)}


_LAST = {}


def figure_values(u, v, system, digits=None):
    """Reuse the jet at a point while its six derivative components are read."""
    import numpy as np
    if digits is None:
        u, v = np.asarray(u, dtype=float), np.asarray(v, dtype=float)
        key = (system, digits, u.shape, u.tobytes(), v.shape, v.tobytes())
    else:
        key = (system, digits, repr(u), repr(v))
    if key not in _LAST:
        if digits is None:
            result = _figure_values(u, v, system, digits)
        else:
            with mp.workdps(digits):
                result = _figure_values(u, v, system, digits)
        if len(_LAST) >= 4:
            _LAST.clear()
        _LAST[key] = result
    return _LAST[key]
