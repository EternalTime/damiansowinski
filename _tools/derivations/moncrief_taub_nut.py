"""Toroidal analytic horizon family, derived from Chrusciel-Lake (2004), (1.4)-(1.5), (4.2).

The two wave amplitudes remain free in printed tensors.
The vacuum equations and analytic horizon conditions select the spacetime family.
The drawing member and its extension are recorded in moncrief_taub_nut.md.
"""
import random

import mpmath as mp
import sympy as sp

import chart_printer as cp
import verify_metrics as vm

SYSTEMS = ('areal', 'extension_plus', 'extension_minus')
AREAL = (r'ds^2=L^2\left(\kappa^2e^{2b+2f}(-dt^2+d\theta^2)'
         r'+t^2e^{2b}(d\psi+q\,d\phi)^2+e^{-2b}d\phi^2\right)')


def spec(system):
    """Areal chart and the two logarithmic continuations, with all coordinates dimensionless."""
    if system == 'areal':
        coords = ['t', r'\theta', r'\psi', r'\phi']
        params = ['L', r'\kappa'] + [f'{f} = {f}(t,\\theta)' for f in ('b', 'f', 'q')]
        line, name = AREAL, 'Toroidal Areal Time'
        domains = [r't \in (0,\infty)']
    else:
        coords = ['T', r'\theta', r'\chi', r'\phi']
        params = ['L', r'\kappa'] + [f'{f} = {f}(T,\\theta)' for f in ('B', 'F', 'Q')]
        sign = '-' if system == 'extension_plus' else '+'
        line = (r'ds^2=L^2\left(e^{2B}\left(T(d\chi+Q\,d\phi)^2'
                + sign + r'\kappa\,dT(d\chi+Q\,d\phi)'
                r'-\dfrac{\kappa^2(e^{2F}-1)}{4T}dT^2\right)'
                r'+\kappa^2e^{2B+2F}d\theta^2+e^{-2B}d\phi^2\right)')
        name = 'Horizon Extension (' + ('+' if system == 'extension_plus' else '-') + ')'
        domains = [r'T \in (-\delta,\delta)']
    domains += [c + r' \in [0,2\pi)' for c in coords[1:]]
    if system != 'areal':
        domains += [r'\delta>0\;\text{within the analytic continuation}']
    return {'metric_id': 'moncrief_taub_nut',
            'system': {'id': system, 'name': name, 'coords': coords, 'domains': domains,
                       'parameters': params, 'line_element': line},
            'chart_line_element': line, 'printer': {},
            'pretty': lambda v: sp.powsimp(sp.factor(sp.sympify(v)), combine='exp'),
            'kretschmann_text': scalar_text,
            'check': lambda chart: check(chart, system)}


def scalar_text(chart):
    """Collect the free-function scalar without factoring its expanded numerator."""
    top, bottom = sp.fraction(sp.together(chart.geo.kretschmann()))
    return chart.printer(sp.powsimp(sp.expand(top) / sp.factor(bottom), combine='exp'))


def equations(x, y, b, f, q, squared=False):
    """Wave equations and both constraints, in areal t or in T=t^2."""
    D = sp.diff
    bx, by, qx, qy = D(b, x), D(b, y), D(q, x), D(q, y)
    e = sp.exp(4*b)
    if squared:
        waves = {D(b, x, 2): (D(b, y, 2)-4*bx+2*x*x*e*qx*qx-x*e*qy*qy/2)/(4*x),
                 D(q, x, 2): (D(q, y, 2)-8*qx-16*x*bx*qx+4*by*qy)/(4*x)}
        first = {x: 2*x*bx*bx+by*by/2+x*e*(4*x*qx*qx+qy*qy)/8,
                 y: 4*x*bx*by+x*x*e*qx*qy}
    else:
        waves = {D(b, x, 2): D(b, y, 2)-bx/x+x*x*e*(qx*qx-qy*qy)/2,
                 D(q, x, 2): D(q, y, 2)-3*qx/x-4*(bx*qx-by*qy)}
        first = {x: x*(bx*bx+by*by)+x**3*e*(qx*qx+qy*qy)/4,
                 y: 2*x*bx*by+x**3*e*qx*qy/2}
    second = {D(f, x, 2): D(first[x], x), D(f, y, 2): D(first[y], y),
              D(f, x, y): D(first[x], y)}
    second = {k: v.subs(waves) for k, v in second.items()}
    return waves | {D(f, k): v for k, v in first.items()} | second


def numeric(expr, values):
    return float(sp.sympify(expr).subs(values, simultaneous=True).evalf(25))


def check(chart, system):
    """Evaluate Ricci on generic vacuum jets, including nonzero q derivatives; check both pullbacks."""
    x, y = chart.symbols[:2]
    names = ('b', 'f', 'q') if system == 'areal' else ('B', 'F', 'Q')
    b, f, q = (chart.reader.parameters[n] for n in names)
    rules = equations(x, y, b, f, q, system != 'areal')
    ricci = [sp.sympify(v).subs(rules, simultaneous=True) for row in chart.geo.ricci_ll() for v in row]
    jets = set().union(*(v.atoms(sp.Derivative) for v in ricci))
    rng = random.Random(1981312)
    for time in (0.7, 1.3, 2.1):
        values = {x: time, y: 0.4, b: 0.13, f: 0.07, q: 0.21,
                  chart.reader.parameters['L']: 1.2, chart.reader.parameters['kappa']: 0.9}
        values.update({j: rng.uniform(0.03, 0.25) for j in sorted(jets, key=str)})
        if max(abs(numeric(v, values)) for v in ricci) > 2e-10:
            raise AssertionError('Moncrief family: vacuum equations leave a Ricci component')
    if system == 'areal':
        return
    s = 1 if system == 'extension_plus' else -1
    source = cp.Chart(['t', r'\theta', r'\psi', r'\phi'],
                      ['L', r'\kappa'] + [f'{v} = {v}(t,\\theta)' for v in ('b', 'f', 'q')], AREAL)
    for time in (0.3, 1.2, 2.0):
        k, length, bv, fv, qv = 0.9, 1.2, 0.13, 0.07, 0.21
        there = {source.symbols[0]: time, source.reader.parameters['L']: length,
                 source.reader.parameters['kappa']: k}
        there.update({source.reader.parameters[n]: v for n, v in zip(('b','f','q'), (bv,fv,qv))})
        here = {x: time*time, chart.reader.parameters['L']: length, chart.reader.parameters['kappa']: k,
                b: bv, f: fv, q: qv}
        J = sp.eye(4)
        J[0, 0], J[2, 0] = 1/(2*time), -s*k/(2*time*time)
        delta = J.T * source.geo.g.subs(there, simultaneous=True) * J - chart.geo.g.subs(here, simultaneous=True)
        if max(abs(float(v)) for v in delta) > 1e-10:
            raise AssertionError('Moncrief family: logarithmic pullback differs')


def drawing_functions(t, theta):
    """Inhomogeneous analytic polarised solution, epsilon=1/4, L=kappa=1."""
    epsilon = mp.mpf(1)/4
    j0, j1 = mp.besselj(0, t), mp.besselj(1, t)
    b = epsilon*j0*mp.cos(theta)
    f = epsilon**2*(t*t*(j0*j0+j1*j1)/2-t*j0*j1*mp.cos(theta)**2)
    return b, f
