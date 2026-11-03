"""Drawings of the current-carrying string's cylindrical exterior."""
import numpy as np

METRIC = 'superconducting_string'
SYSTEM = 'cylindrical'
PARAMS = {'m': '1/5', 'k': '7/10', 'W': '4/5', 'ell': '1', 'r_c': '1'}
SETTINGS = r'$m=1/5$, $k=7/10$, $W=4/5$, $\ell=r_c=1$.'


def register_null(nr):
    row = nr.Diagram(METRIC, SYSTEM, 'radial', r'$t$ and $r$', ('t', 'r'),
                     (0, 6, -3, 3), r'$r/\ell$', r'$t/\ell$', PARAMS,
                     {'phi': '0', 'z': '0'})
    nr.DIAGRAMS.append(row)
    nr.CAPTIONS[(METRIC, SYSTEM, 'radial')] = [
        r'The plane of $t$ and $r$ ($\phi=z=0$), with radial null geodesics '
        r'$t\mp r=\mathrm{constant}$. The exterior begins beyond $r_c$.']
    nr.CLOSED_FORMS[(METRIC, SYSTEM, 'radial')] = (
        lambda t, r: t+r, lambda t, r: t-r, lambda t, r: r > 1)


def register_embedding(eb):
    eb.DRAWN[METRIC] = lambda ck, src: embedding(eb, ck, src)
    eb.CAPTIONS[(METRIC, 'exterior')] = [
        r'The plane $z=0$ at $t=0$, outside the inner boundary $r=r_c$.']


def embedding(eb, ck, src):
    """Embed the exterior annulus where its circles fit in flat space."""
    sl = eb.Slice(src, METRIC, SYSTEM, 'r', r'\phi', {'t': 0, 'z': 0}, PARAMS)
    sample = np.linspace(1, 6, 501)
    ck.add('String: the displayed annulus embeds in flat space',
           max(0., float(np.max(-sl.defect_at(sample)))), 1e-12)
    size = 2 * float(sl.rho_at(6))
    piece = eb.Piece('exterior', 'sheet', sl, 1, 6, 0., 1,
                     (('edge', 'the inner limit $r=r_c$'), ('edge', 'the exterior continues outward')),
                     [(1, 'chartedge', None), (2, 'r', None), (4, 'r', None), (6, 'r', None)], size)
    ck.isometry('String: exterior annulus', piece)
    surface = eb.Surface([piece])
    fig = eb.figure_of([surface], {'sheet': 'cover'}, size)
    fig.legend('fill', 'cover', 'the exterior slice')
    fig.legend('line', 'chartedge', '$r=r_c$')
    fig.legend('line', 'r', '$r$ constant')
    fig.legend('line', 'meridian', '$\\phi$ constant')
    return [eb.view('exterior', 'The plane $z=0$', r'$\ell$', [surface], fig.done(), settings=SETTINGS)]


def register_slices(sl):
    sl.FLAT[(METRIC, SYSTEM, 'radial')] = lambda: sl.one(
        METRIC, lambda moment: sl.along(0., *moment.reach(SYSTEM, 'r')), label=r'$t=0$')
    sl.FLAT_METRICS.add(METRIC)


def register_conformal(cf):
    cf.DRAWN[METRIC] = lambda ck, src: conformal(cf, ck, src)
    cf.CAPTIONS[(METRIC, 'exterior')] = [
        r'The plane of $t$ and $r$ ($\phi=z=0$), compactified by '
        r'$p,q=\arctan((t\mp(r-r_c))/\ell)$.']


def conformal(cf, ck, src):
    """Compactify the totally geodesic exterior half plane at r_c=ell=1."""
    def pq(t, r):
        return cf.mink_pq(t, np.asarray(r, dtype=float)-1)
    plane = cf.Plane(src, METRIC, SYSTEM, ('t', 'r'), {'phi': '0', 'z': '0'}, PARAMS)
    ck.chart('String: exterior radial half plane', plane, pq,
             ck.uniform(-20, 20), ck.uniform(1.01, 20), lambda t, r: (1, 0))
    v = cf.View('exterior', 'Cylindrical exterior', [-.5, cf.PI+.35, -cf.PI-.25, cf.PI+.25], SYSTEM)
    v.fill('region', cf.TRIANGLE)
    v.fill('cover', cf.TRIANGLE)
    cf.grid(v, 'r', lambda r, t: pq(t, r), (2, 3, 5), cf.S_ALL)
    cf.grid(v, 't', pq, (-4, -2, -1, 0, 1, 2, 4), cf.S_POS+1)
    cf.triangle_edges(v, centre_class='surface')
    v.legend('cover', 'the exterior half plane')
    v.legend('r', '$r$ constant')
    v.legend('t', '$t$ constant')
    v.legend('surface', 'the inner limit $r=r_c$')
    moment = cf.slices.moments(METRIC)[0]
    r = np.linspace(*moment.reach(SYSTEM, 'r'), 2)
    v.slice(moment, [pq(0*r, r)])
    v.set(restriction=r'The totally geodesic half plane of fixed $\phi$ and $z$ in the exterior.',
          settings=SETTINGS)
    return [v]
