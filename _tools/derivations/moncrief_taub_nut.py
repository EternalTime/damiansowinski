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
        params = ['L', r'\kappa', r'\delta'] + [f'{f} = {f}(T,\\theta)' for f in ('B', 'F', 'Q')]
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


# The extension drawings use the totally geodesic section theta=phi=0.
# F/T is evaluated by its entire power series so the coefficient H stays regular at T=0.
INPUT = r'Polarised Bessel member ($\epsilon=1/4$, $q=0$, $L=\kappa=1$).'


def areal_functions():
    return {'b': 'besselj(0,t)*cos(theta)/4', 'q': '0',
            'f': '(t**2*(besselj(0,t)**2+besselj(1,t)**2)/2-t*besselj(0,t)*besselj(1,t)*cos(theta)**2)/16'}


def extension_functions():
    """Declared functions on theta=0; F's transverse first derivative also vanishes there."""
    return {'B': 'mon_j0(T)*cos(theta)/4', 'Q': '0', 'F': 'log(1+4*T*mon_h0(T))/2'}


def f_over_T(T):
    import math
    import numpy as np
    # D(T)=J1(sqrt(T))/sqrt(T); F_T=T*D(T)^2/32 on theta=0.
    d = [(-1)**n / (2*4**n*math.factorial(n)*math.factorial(n+1)) for n in range(40)]
    square = np.convolve(d, d)[:40]
    coeff = [c/(32*(n+2)) for n, c in enumerate(square)]
    return np.asarray(T)*np.polynomial.polynomial.polyval(T, coeff)


def j0_T(T, order=0):
    import numpy as np
    from scipy.special import hyp0f1
    import math
    return (-0.25)**order/math.factorial(order)*hyp0f1(1+order, -np.asarray(T)/4)


def h0_T(T, order=0):
    import numpy as np
    from scipy.special import exprel
    if order:
        def derivative(x):
            with mp.workdps(40):
                def value(z):
                    if not z:
                        return mp.mpf(0)
                    a, d = mp.hyp0f1(1,-z/4), mp.hyp0f1(2,-z/4)/2
                    f = (z*(a*a+z*d*d)/2-z*a*d)/16
                    return mp.expm1(2*f)/(4*z)
                return float(mp.diff(value, mp.mpf(float(x)), order))
        return np.vectorize(derivative, otypes=[float])(T)
    ft = f_over_T(T)
    return ft*exprel(2*np.asarray(T)*ft)/2


def conformal_w(T):
    """Integral of (exp(F)-1)/(2T), regular on the horizon section."""
    import numpy as np
    from scipy.special import exprel
    nodes, weights = np.polynomial.legendre.leggauss(24)
    ts = np.asarray(T)[..., None]*(nodes+1)/2
    ft = f_over_T(ts)
    values = ft*exprel(ts*ft)/2
    return np.asarray(T)*np.sum(values*weights, axis=-1)/2


def null_rows(nr):
    params, fixed = {'L': 1, 'kappa': 1}, {'theta': 0, 'phi': 0}
    rows = [nr.Diagram('moncrief_taub_nut', 'areal', 'compact', '$t$ and $\\psi$',
                       ('t', r'\psi'), (0, 2*nr.math.pi, 0, 4), '$\\psi$', '$t$', params, fixed,
                       tau='t', families=nr.SIDEWAYS, functions=areal_functions(), input=INPUT,
                       periodic=(r'\psi',), kretschmann=False, where='t', inside=True)]
    for system in SYSTEMS[1:]:
        rows.append(nr.Diagram('moncrief_taub_nut', system, 'compact', '$T$ and $\\chi$',
                               ('T', r'\chi'), (0, 2*nr.math.pi, -2, 2), '$\\chi$', '$T$', params | {'delta': 3}, fixed,
                               orient='outgoing' if system.endswith('plus') else 'ingoing',
                               families=nr.SIDEWAYS, functions=extension_functions(), input=INPUT,
                               
                               periodic=(r'\chi',), kretschmann=False, mark_g00=False,
                               lines=(('horizon','x0','0', '$T=0$, compact Cauchy horizon'),)))
    return rows


def null_captions():
    captions = {('moncrief_taub_nut','areal','compact'): [
        r'The section $\theta=\phi=0$, with the edges $\psi=0,2\pi$ identified and the compact direction shrinking toward the Cauchy horizon $t=0$.']}
    for system in SYSTEMS[1:]:
        captions[('moncrief_taub_nut',system,'compact')] = [
            r'The regular section $\theta=\phi=0$, with the edges $\chi=0,2\pi$ identified; spacelike circles at $T>0$, the compact null horizon at $T=0$, and closed timelike circles at $T<0$.']
    return captions


def embedding_view(em, ck, src):
    """The isometric surface of t=1, phi=0, cut at theta=0 and 2pi with ends identified."""
    import math
    import numpy as np
    sl = em.Slice(src, 'moncrief_taub_nut', 'areal', r'\theta', r'\psi',
                  {'t': 1, 'phi': 0}, {'L': 1, 'kappa': 1}, areal_functions())
    size, top = 8.0, 2*math.pi
    ck.add('Moncrief: positive embedding defect on the entire cut torus',
           max(0.0, -float(np.min(sl.defect_at(np.linspace(0, top, 721))))), 0.0)
    tube = em.Piece('torus', 'sheet', sl, 0.0, top, 0.0, 1,
                    (('edge', r'$\theta=0$, identified with the other end'),
                     ('edge', r'$\theta=2\pi$, identified with the other end')),
                    [(math.pi, 'r', None)], size)
    ck.isometry('Moncrief cut torus', tube)
    ck.radius('Moncrief: rho=e^b at t=1', tube,
              lambda theta: np.exp(float(mp.besselj(0,1))*np.cos(theta)/4), size)
    ck.add('Moncrief: identified ends have the same radius', abs(tube.rho[0]-tube.rho[-1]), em.FORM)
    surface = em.Surface([tube])
    fig = em.figure_of([surface], {'sheet': 'cover'}, size, meridians=12)
    fig.legend('fill', 'cover', r'$t=1$, $\phi=0$, with the two ends identified')
    fig.legend('line', 'r', r'$\theta=\pi$')
    fig.legend('line', 'meridian', r'$\psi$ constant')
    return [em.view('cut_torus', 'Toroidal spatial slice', '$L$', [surface], fig.done(),
                    settings=r'$L=\kappa=1$, $t=1$, $\phi=0$.', input=INPUT)]


def conformal_view(cf, ck, src):
    """Null coordinates on the regular theta=phi=0 section, with its compact identification.

    w'= (exp(F)-1)/(2T), U=exp(chi+w), V=T exp(-chi+w) for the plus extension.
    For the minus extension U=T exp(chi+w), V=exp(-chi+w).
    Thus ds^2=-exp(2B-2w)dU dV, and chi+2pi is the boost (U,V)->(e^2pi U,e^-2pi V).
    The boundaries T=+-1 are finite section cuts, not an assumed global continuation.
    """
    import numpy as np
    views = []
    period = 2*np.pi
    for system in SYSTEMS:
        areal = system == 'areal'
        s = 1 if system != 'extension_minus' else -1
        def regular(T, chi, s=s):
            T, chi = np.asarray(T), np.asarray(chi)
            w = conformal_w(T)
            U, V = ((np.exp(chi+w), T*np.exp(-chi+w)) if s == 1
                    else (T*np.exp(chi+w), np.exp(-chi+w)))
            return np.arctan(U), np.arctan(V)
        if areal:
            # Rebind regular in the function default so the loop cannot change it later.
            def fmap(t, psi, regular=regular):
                return regular(np.asarray(t)**2, np.asarray(psi)+np.log(t))
            plane = cf.Plane(src, 'moncrief_taub_nut', system, ('t', r'\psi'),
                             {'theta': 0, 'phi': 0}, {'L': 1, 'kappa': 1}, areal_functions(), induced=True)
            ck.chart('Moncrief areal section', plane, fmap, ck.uniform(0.04,1,80),
                     ck.uniform(0,period,80), lambda t, psi: (1,0))
            x = np.linspace(0.001,1,401)
            across = np.linspace(0,period,401)
            boundary = [*zip(*fmap(x,0*x)), *zip(*fmap(1+0*across,across)),
                        *zip(*fmap(x[::-1],0*x+period)), *zip(*fmap(0*across+0.001,across[::-1]))]
        else:
            fmap = regular
            plane = cf.Plane(src, 'moncrief_taub_nut', system, ('T', r'\chi'),
                             {'theta': 0, 'phi': 0}, {'L': 1, 'kappa': 1}, extension_functions(), induced=True)
            ck.chart('Moncrief '+system+' section', plane, fmap, ck.uniform(-1,1,80),
                     ck.uniform(0,period,80), lambda T, chi, s=s: (1,s/(1+np.abs(T))))
            x = np.linspace(-1,1,401)
            across = np.linspace(0,period,401)
            boundary = [*zip(*fmap(x,0*x)), *zip(*fmap(1+0*across,across)),
                        *zip(*fmap(x[::-1],0*x+period)), *zip(*fmap(0*across-1,across[::-1]))]
        v = cf.View(system, plane.entry['name'], (-1.9,1.9,-1.9,1.9), system)
        v.fill('region', [cf.point(*pq) for pq in boundary])
        for edge in (0,period):
            p,q=fmap(x,0*x+edge)
            v.curve('r',p,q)
        for edge in ((0.25,0.5,1) if areal else (-1,-0.5,0,0.5,1)):
            p,q=fmap(0*across+edge,across)
            v.curve('horizon' if edge == 0 else 't',p,q)
        v.legend('r', r'identified angular edges')
        v.legend('t', r'finite cuts at constant time')
        if not areal:
            v.legend('horizon', r'$T=0$, compact Cauchy horizon')
        v.set(settings=INPUT, restriction=r'The section $\theta=\phi=0$, with its periodic angular edges identified.')
        v.set(caption=[r'The local horizon section in conformal null coordinates, with finite time cuts.'])
        for moment in cf.slices.moments('moncrief_taub_nut', label='$t=1$'):
            p,q=fmap(0*across+1,across)
            v.slice(moment, [(p,q)])
        views.append(v)
    return views


def curvature_view(cv, ck, src):
    """The stored Kretschmann scalar of the Bessel member on a finite areal-time rectangle."""
    import math
    import numpy as np
    wave = areal_functions()
    K = cv.nr_published(src, 'moncrief_taub_nut', 'areal', 'kretschmann',
                        ['theta','t'], ['L','kappa'], functions=wave)
    xaxis, yaxis = cv.axis('$\\theta$', 0,2*math.pi), cv.axis('$t$',0.2,4)
    nx, ny = 96,64
    xs,ys,values=cv.cells(K,xaxis,yaxis,nx,ny,[1,1])
    ck('Moncrief stored K finite in every cell', np.all(np.isfinite(values)))
    bound=10**math.ceil(math.log10(max(1e-6,float(np.max(np.abs(values))))))
    scale={'kind':'linear','from':-bound,'to':bound}
    cv.shape_checks(ck,'moncrief_taub_nut',[(2*math.pi)/3.8])
    return [{'id':'bessel','label':'Toroidal Bessel member','kind':'map','x':xaxis,'y':yaxis,
             'colour':{'name':'$K L^4$','scale':scale,'ticks':cv.value_ticks(scale)},
             'map':{'metric':'moncrief_taub_nut','system':'areal','scalar':'kretschmann',
                    'x':r'\theta','y':'t','params':{'L':1,'kappa':1},'functions':wave,
                    'nx':nx,'ny':ny,'cells':[[cv.rounded(v) for v in row] for row in values]},
             'marks':[],'legend':[], 'settings':INPUT,
             'clip':'Colour range covering the displayed values.',
             'caption':[r'The Kretschmann scalar $K$ of the polarised Bessel member over one spatial period, the same at every $\psi$ and $\phi$.']}]
