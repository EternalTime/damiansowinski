"""Charts of the modular string geometry and pointwise field equation checks."""
import sympy as sp

METRIC = 'stringy_cosmic_strings'
SYSTEMS = ('cartesian', 'polar')
DIMENSIONS = {
    (METRIC, 'cartesian'): {'t': 'T', 'w': 'L', 'x': 'L', 'y': 'L', 'H': '1'},
    (METRIC, 'polar'): {'t': 'T', 'w': 'L', 'r': 'L', '\\phi': '1', 'H': '1'},
}


def chart_spec(system):
    import verify_metrics as vm
    import chart_printer as cp
    if system == 'cartesian':
        coords = ['t', 'w', 'x', 'y']
        params = ['H = H(x,y)']
        plane = 'H\\left(dx^2 + dy^2\\right)'
        domains = ['x \\in (-\\infty, \\infty)', 'y \\in (-\\infty, \\infty)']
        name = 'Cartesian transverse coordinates'
    else:
        coords = ['t', 'w', 'r', '\\phi']
        params = ['H = H(r,\\phi)']
        plane = 'H\\left(dr^2 + r^2d\\phi^2\\right)'
        domains = ['r \\in (0, \\infty)', '\\phi \\in [0, 2\\pi)']
        name = 'Polar transverse coordinates'
    line = 'ds^2 = -c^2dt^2 + dw^2 + ' + plane
    probe = vm.Reader(coords, params, ())
    H = probe.parameters['H']
    return {'metric_id': METRIC,
            'system': {'id': system, 'name': name, 'coords': coords,
                       'domains': ['t \\in (-\\infty, \\infty)', 'w \\in (-\\infty, \\infty)'] + domains,
                       'parameters': params, 'line_element': line},
            'chart_line_element': line.replace('c^2dt^2', 'dt^2'),
            'printer': {'lead': [H], 'collect': lambda poly, printer: cp.collect_by(poly, [H], printer)},
            'check': lambda chart: check(chart, system)}


def check(chart, system):
    """The local tau=z solution obeys R_ab=(d_a tau d_b tau* + d_b tau d_a tau*)/(4 Im(tau)^2)."""
    t, w, a, b = chart.symbols
    H = chart.reader.parameters['H']
    if system == 'cartesian':
        sample_h = b
        points = ((sp.Rational(2, 3), sp.Rational(3, 2)), (sp.Rational(-1, 2), sp.Rational(7, 3)))
        expected = lambda av, bv: sp.diag(0, 0, 1 / (2 * bv**2), 1 / (2 * bv**2))
    else:
        sample_h = a * sp.sin(b)
        points = ((sp.Rational(2, 1), sp.pi / 4), (sp.Rational(3, 2), sp.pi / 3))
        def expected(av, bv):
            factor = 1 / (2 * sp.sin(bv)**2)
            return sp.diag(0, 0, factor / av**2, factor)
    ricci = chart.geo.ricci_ll()
    for av, bv in points:
        actual = sp.Matrix(4, 4, lambda i, j: ricci[i][j].subs(H, sample_h).doit().subs({a: av, b: bv}))
        error = (actual - expected(av, bv)).evalf(30)
        if max(abs(complex(v)) for v in error) > 1e-24:
            raise AssertionError(f'{METRIC}/{system}: Ricci field equation failed at {(av, bv)}')
