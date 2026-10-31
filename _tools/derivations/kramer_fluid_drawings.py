"""Interior-only drawings of the bounded rotating fluid at m=1/2, b=6.1."""
import math

import numpy as np
from scipy.optimize import brentq

METRIC = 'kramer_fluid'
B = 6.1
ETA_EDGE = brentq(lambda eta: (eta+1)*math.exp(eta)-B, 1, 2, xtol=1e-14)
W_EDGE = math.sqrt(ETA_EDGE-1)
SETTINGS = r'$m=1/2$, $b=6.1$, $\ell=1$.'


def register_null(nr):
    for system, coordinate, lo, hi in [('published', r'\eta', 1, ETA_EDGE), ('radial', 'w', 0, W_EDGE)]:
        pad = .107*(hi-lo)
        span = hi-lo+2*pad
        row = nr.Diagram(METRIC, system, 'equator', 'Equatorial projection', (r'\tau', coordinate),
                         (lo-pad, hi+pad, -span/2, span/2), f'${coordinate}$', r'$\tau$',
                         {'m': '.5', 'b': '6.1', 'eta_e': repr(ETA_EDGE)}, {'xi': '0'},
                         quotient='phi', tau='tau', where=f'({nr.vm.Reader._plain(coordinate)}-({lo}))*(({hi})-{nr.vm.Reader._plain(coordinate)})',
                         cone='future cone of the quotient',
                         lines=(('surface', 'r', repr(lo), 'regular axis'), ('surface', 'r', repr(hi), '$p=0$')))
        nr.DIAGRAMS.append(row)
        nr.CAPTIONS[(METRIC, system, 'equator')] = [
            f'The equatorial projection in $\\tau$ and ${coordinate}$ ($\\xi=0$), '
            r'with $\phi$ divided out; null geodesics of zero angular momentum within the fluid.']


def register_embedding(eb):
    eb.DRAWN[METRIC] = lambda ck, src: embedding(eb, ck, src)
    eb.CAPTIONS[(METRIC, 'disk')] = [
        r'The equatorial disk ($\tau=0$, $\xi=0$), from the regular axis to the surface of zero pressure.']


def embedding(eb, ck, src):
    def stable(value, w):
        # exp(-u)(exp(u)-1)^2=4 sinh(u/2)^2 avoids subtracting equal constants at the axis.
        a, u = eb.sp.Rational(61,10)/eb.sp.E, w*w
        gpp = 4*(u-4*a*eb.sp.sinh(u/2)**2)
        defect = 16*a*(u*eb.sp.sinh(u)-eb.sp.sinh(u/2)**2-a*u*eb.sp.sinh(u)**2)/(u-4*a*eb.sp.sinh(u/2)**2)
        new = gpp if not stable.calls else defect
        stable.calls += 1
        for q in ('0.000001', '0.01', '0.1', str(W_EDGE)):
            point = eb.sp.Float(q,65)
            gap = abs((value-new).subs(w,point).evalf(55))
            ck.add('Kramer: stable induced metric expression', float(gap), 1e-40)
        return new
    stable.calls = 0
    sl = eb.Slice(src, METRIC, 'radial', 'w', r'\phi', {'tau': 0, 'xi': 0},
                  {'m': '1/2', 'b': '61/10'}, rewrite=stable)
    # The removable central limits are part of the regular polar axis.
    original_rho, original_gpp, original_defect = sl._rho, sl._gpp, sl._defect
    sl._rho = lambda w: np.where(np.asarray(w)==0, 0., original_rho(w))
    sl._gpp = lambda w: np.where(np.asarray(w)==0, 0., original_gpp(w))
    sl._defect = lambda w: np.where(np.asarray(w)==0, 0., original_defect(w))
    # The source's induced constant-time metric, equation (15), not its fluid rest-space quotient.
    sample = np.linspace(1e-4, W_EDGE, 601)
    ck.add('Kramer: Euclidean embedding condition on the whole equatorial disk',
           max(0., float(np.max(-sl.defect_at(sample)))), 1e-12)
    size = 1.5
    disk = eb.Piece('disk', 'star', sl, 0., W_EDGE, 0., 1,
                    (('axis', 'regular axis'), ('edge', '$p=0$')),
                    [(W_EDGE/2, 'r', None), (W_EDGE, 'surface', '$p=0$')], size)
    ck.isometry('Kramer: induced equatorial disk', disk)
    surface = eb.Surface([disk])
    fig = eb.figure_of([surface], {'star': 'star'}, size, eb.Camera(-90, 35))
    fig.legend('fill', 'star', 'fluid disk')
    fig.legend('line', 'surface', '$p=0$')
    fig.legend('line', 'r', '$w$ constant')
    fig.legend('line', 'meridian', r'$\phi$ constant')
    return [eb.view('disk', 'Equatorial disk', r'$\ell$', [surface], fig.done(), settings=SETTINGS)]


def register_slices(sl):
    sl.FLAT[(METRIC, 'radial', 'equator')] = lambda: sl.one(
        METRIC, lambda m: sl.along(0., *m.reach('radial', 'w')))
    sl.FLAT[(METRIC, 'published', 'equator')] = lambda: sl.one(
        METRIC, lambda m: sl.along(0., *[1+q*q for q in m.reach('radial', 'w')]))
    sl.FLAT_METRICS.add(METRIC)


def register_curvature(cv):
    cv.DRAWN[METRIC] = lambda ck, src: curvature(cv, ck, src)


def curvature(cv, ck, src):
    """Read the stored scalars on the actual fluid body, with a plain coordinate box."""
    x_axis, y_axis = cv.axis('$w$', -.08, .52), cv.axis(r'$\xi$', -.55, .55)
    nx, ny = 108, 198
    views = []
    for field, letter, name in [('kretschmann', 'K', 'Kretschmann'), ('ricci_scalar', 'R', 'Ricci')]:
        stored = cv.Published(src, METRIC, 'radial', field, ['w', 'xi'], ['m', 'b'])
        xs, ys, values = cv.cells(stored, x_axis, y_axis, nx, ny, [.5, B])
        W, XI = np.meshgrid(xs, ys)
        inside = (W >= 0) & (B*np.cos(XI)*np.exp(-1-W*W) >= 2+W*W)
        kept = values[inside]
        ck(f'Kramer {letter}: finite stored values throughout the fluid', np.all(np.isfinite(kept)))
        low, high = float(np.min(kept)), float(np.max(kept))
        width = max(high-low, .1)
        scale = {'kind': 'linear', 'from': low-.1*width, 'to': high+.1*width}
        rows = [[cv.rounded(v) if ok else None for v, ok in zip(row, mask)] for row, mask in zip(values, inside)]
        views.append({'id': letter, 'label': f'${letter}$', 'kind': 'map', 'x': x_axis, 'y': y_axis,
                      'colour': {'name': f'${letter}\\ell^{4 if letter=="K" else 2}$', 'scale': scale,
                                 'ticks': cv.value_ticks(scale)},
                      'map': {'metric': METRIC, 'system': 'radial', 'scalar': field, 'x': 'w', 'y': r'\xi',
                              'params': {'m': .5, 'b': B},
                              'domain': '(w >= 0) & (b*cos(xi)*exp(-1-w**2) >= 2+w**2)',
                              'nx': nx, 'ny': ny, 'cells': rows},
                      'marks': [], 'legend': [],
                      'caption': [f'The {name} scalar ${letter}$ across the bounded fluid in $w$ and $\\xi$.'],
                      'settings': SETTINGS, 'clip': 'Linear colour scale over the fluid interior.'})
    return views
