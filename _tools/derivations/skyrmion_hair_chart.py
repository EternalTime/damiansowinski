"""Areal exterior of the spherical Einstein-Skyrme black hole."""
import sympy as sp

import verify_metrics as vm


def chart_spec():
    coords = ['t', 'r', '\\theta', '\\phi']
    parameters = ['N = N(r)', 'C = C(r)', 'f = f(r)', 'm = m(r)',
                  'r_h', 'G', 'M', '\\alpha', 'F_\\pi', 'a']
    probe = vm.Reader(coords, parameters, ())
    r = probe.symbol['r']
    n, c = probe.parameters['N'], probe.parameters['C']
    dn, dc = sp.Derivative(n, r), sp.Derivative(c, r)
    return {
        'metric_id': 'skyrmion_hair',
        'system': {
            'id': 'areal', 'name': 'Areal Radius', 'coords': coords,
            'domains': ['t \\in (-\\infty, \\infty)', 'r \\in (r_h, \\infty)',
                        '\\theta \\in [0, \\pi]', '\\phi \\in [0, 2\\pi)'],
            'parameters': parameters,
            'line_element': ('ds^2 = -N^2 C\\,dt^2 + \\dfrac{dr^2}{C}'
                             ' + r^2\\left(d\\theta^2 + \\sin^2\\theta\\,d\\phi^2\\right)'),
        },
        'chart_line_element': ('ds^2 = -N^2 C\\,dt^2 + \\dfrac{dr^2}{C}'
                               ' + r^2\\left(d\\theta^2 + \\sin^2\\theta\\,d\\phi^2\\right)'),
        'printer': {'lead': [sp.Derivative(n, (r, 2)), sp.Derivative(c, (r, 2)), dn, dc, r, n, c],
                    'factors': [r, n, c, dn, dc]},
        'check': check,
    }


def check(chart):
    """Check the radial Einstein equations at independent numerical points."""
    r = chart.reader.symbol['r']
    n, c = chart.reader.parameters['N'], chart.reader.parameters['C']
    dn, dc = sp.Derivative(n, r), sp.Derivative(c, r)
    mixed = chart.geo.ginv * sp.Matrix(chart.geo.einstein_ll())
    expected = ((r*dc + c - 1)/r**2,
                (r*dc + c - 1)/r**2 + 2*c*dn/(r*n))
    for point in ({r: sp.Rational(2), n: sp.Rational(4, 5), c: sp.Rational(2, 3),
                   dn: sp.Rational(1, 17), dc: sp.Rational(1, 9)},
                  {r: sp.Rational(5, 2), n: sp.Rational(9, 10), c: sp.Rational(3, 4),
                   dn: sp.Rational(1, 27), dc: sp.Rational(1, 7)}):
        for i in range(2):
            if abs(complex((mixed[i, i] - expected[i]).xreplace(point).evalf(30))) > 1e-20:
                raise AssertionError('skyrmion_hair/areal: radial Einstein equation')


def tortoise_spec():
    """The tortoise radius in Shiiki and Sawado's equation (54)."""
    coords = ['t', '\\xi', '\\theta', '\\phi']
    parameters = ['F = F(\\xi)', 'r = r(\\xi)']
    probe = vm.Reader(coords, parameters, ())
    x = probe.symbol['\\xi']
    f, r = probe.parameters['F'], probe.parameters['r']
    first = [sp.Derivative(v, x) for v in (f, r)]
    lead = [sp.Derivative(v, (x, 2)) for v in (f, r)] + first + [f, r]
    line = ('ds^2 = F\\left(-dt^2 + d\\xi^2\\right)'
            ' + r^2\\left(d\\theta^2 + \\sin^2\\theta\\,d\\phi^2\\right)')
    return {
        'metric_id': 'skyrmion_hair',
        'system': {'id': 'tortoise', 'name': 'Tortoise', 'coords': coords,
                   'domains': ['t \\in (-\\infty, \\infty)', '\\xi \\in (-\\infty, \\infty)',
                               '\\theta \\in [0, \\pi]', '\\phi \\in [0, 2\\pi)'],
                   'parameters': parameters, 'line_element': line},
        'chart_line_element': line,
        'printer': {'lead': lead, 'factors': [f, r] + first},
        'check': check_tortoise,
    }


def check_tortoise(chart):
    """The radial plane is conformally flat and the sphere has area 4 pi r^2."""
    f, r = chart.reader.parameters['F'], chart.reader.parameters['r']
    expected = sp.diag(-f, f, r**2, r**2*sp.sin(chart.reader.symbol['\\theta'])**2)
    if chart.geo.g != expected:
        raise AssertionError('skyrmion_hair/tortoise: metric transformation')
