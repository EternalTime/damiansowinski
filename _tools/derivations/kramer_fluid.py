"""The bounded 1984 rotating fluid, from Chinea and Pareja (1999), (1)-(4)."""
import sympy as sp

METRIC = 'kramer_fluid'
SYSTEMS = ('published', 'radial')


def parameters(system):
    radial = system == 'radial'
    return ['m', 'b'] + ([r'\eta = 1+w^2'] if radial else []) + [
        r'Q = e^{w^2}' if radial else r'Q = e^{\eta-1}',
        r'F = b\cos\xi\,e^{-1-w^2}' if radial else r'F = b\cos\xi\,e^{-\eta}', r'\eta_e']


def register_dimensions(vm):
    for system in SYSTEMS:
        key = (METRIC, system)
        x = 'w' if system == 'radial' else r'\eta'
        vm.HELD[key] = ('Q', 'F')
        vm.RATES[key] = {'Q': {x: '2wQ' if system == 'radial' else 'Q'},
                         'F': {x: '-2wF' if system == 'radial' else '-F', r'\xi': r'-F\tan\xi'}}
        vm.DIMENSIONS[key] = {r'\tau': '1', x: '1', r'\xi': '1', r'\phi': '1',
                              'm': 'L^-2', 'b': '1', 'Q': '1', 'F': '1', r'\eta_e': '1', r'\eta': '1'}


def spec(system):
    import verify_metrics as vm
    radial = system == 'radial'
    x = 'w' if radial else r'\eta'
    coords = [r'\tau', x, r'\xi', r'\phi']
    radial_term = r'4dw^2' if radial else r'\dfrac{d\eta^2}{\eta-1}'
    line = (r'ds^2=\dfrac{1}{2m}\left((\eta-1-F)d\tau^2'
            r'+4(\eta-1-F+FQ)d\tau\,d\phi'
            r'+4(\eta-1-F(Q-1)^2)d\phi^2+' + radial_term + r'+\dfrac{d\xi^2}{F}\right)')
    probe = vm.Reader(coords, parameters(system), (), held=('Q', 'F'), rates=vm.RATES[(METRIC, system)])
    return {'metric_id': METRIC, 'system': {'id': system,
            'name': 'Square Root Radial Coordinate' if radial else 'Published Radial Coordinate',
            'coords': coords, 'parameters': parameters(system), 'line_element': line,
            'domains': [r'\tau\in(-\infty,\infty)',
                        r'w\in[0,\sqrt{\eta_e-1}]' if radial else r'\eta\in[1,\eta_e]',
                        r'|\xi|\leq\arccos\left(\dfrac{(\eta+1)e^{\eta}}{b}\right)',
                        r'\phi\in[0,2\pi)', r'm>0,\quad 2e<b<b_\star',
                        r'\eta=1\;\text{(regular rotation axis, where the polar chart fails)}']},
            'chart_line_element': line,
            'printer': {'lead': [probe.parameters['m'], probe.parameters['F'], probe.parameters['Q'], probe.symbol[x]]},
            'reduce': lambda value: vm.norm(probe.by_rates(value)),
            'check': lambda chart: check(chart, system)}


def check(chart, system):
    """Point checks against the source metric and isotropic fluid, without symbolic checks."""
    R = chart.reader
    m, b, Q, F = (R.parameters[k] for k in ('m', 'b', 'Q', 'F'))
    x, xi = chart.symbols[1:3]
    eta = 1+x*x if system == 'radial' else x
    G = sp.Matrix(chart.geo.einstein_ll())
    for bv, ev, xiv, mv in [(6.1, 1.03, .09, .5), (7.2, 1.12, -.18, .8), (8.2, 1.22, .13, .3)]:
        ev, xiv, bv, mv = [sp.Rational(str(q)) for q in (ev, xiv, bv, mv)]
        h = ev-1
        fv = bv*sp.cos(xiv)*sp.exp(-ev)
        at = {m: mv, b: bv, x: sp.sqrt(h) if system == 'radial' else ev, xi: xiv,
              Q: sp.exp(h), F: fv, chart.symbols[0]: 0, chart.symbols[3]: 0}
        def matrix(value):
            return sp.Matrix(value).subs(at).evalf(40)
        def near(left, right, label):
            gap = matrix(left)-matrix(right)
            if max(abs(complex(v)) for v in gap) > 1e-24:
                raise AssertionError(f'Kramer {system}: {label}: {gap}')
        # Sourced equation (1), independently expanded in its exponentials.
        B = 2*(h-bv*sp.cos(xiv)*(sp.exp(-ev)-sp.exp(-1)))
        C = 4*(h-bv*sp.cos(xiv)*(sp.exp(-ev)+sp.exp(ev-2)-2*sp.exp(-1)))
        original = sp.Matrix([[h-fv, 0, 0, B], [0, 4 if system == 'radial' else 1/h, 0, 0],
                              [0, 0, 1/fv, 0], [B, 0, 0, C]])/(2*mv)
        g = matrix(chart.geo.g)
        near(chart.geo.g, original, 'source metric and radial pullback')
        if not (g[0, 0] < 0 and g[3, 3] > 0):
            raise AssertionError('Kramer: stationary or azimuthal orbit has wrong sign')
        near([[g[0,0]*g[3,3]-g[0,3]**2]], [[-h*fv*sp.exp(2*h)/mv**2]], 'orbit determinant')
        # Equations (2)-(4) use exactly the standard Einstein coupling kappa0=8 pi G/c^4.
        kp = mv*(1+ev-fv)/2
        krho = mv*(1-3*ev+3*fv)/2
        velocity = g[:, 0]/sp.sqrt(-g[0, 0])
        near(G, (krho+kp)*velocity*velocity.T+kp*g, 'perfect fluid and source coupling')


def register(pc):
    pc.CHARTS[METRIC] = [lambda system=system: spec(system) for system in SYSTEMS]
