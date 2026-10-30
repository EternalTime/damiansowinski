"""Matched static stars at GM/(c^2 r_b)=1/5, in units r_b=1."""
import math

import numpy as np

import spheroidal_stars as stars

METRIC = stars.METRIC
SETTINGS = "$r_b=1$, $GM/(c^2r_b)=1/5$."


class StableScalar:
    """Preserve the removable central cancellations in the stored scalar."""
    def __init__(self, published):
        self.published = published

    def __getattr__(self, name):
        return getattr(self.published, name)

    def __call__(self, *values):
        arrays = np.broadcast_arrays(*values)
        out = self.published(*values).copy()
        near = np.abs(arrays[0]) < 0.01
        for index in np.ndindex(out.shape):
            if near[index]:
                out[index] = float(self.published.exact([a[index] for a in arrays]))
        return out


def register_null(nr):
    for system in stars.SYSTEMS:
        row = nr.Diagram(METRIC, system, "radial", "$t$ and $r$", ("t", "r"),
                         (0, 1, -0.5, 0.5), "$r/r_b$", "$ct/r_b$",
                         stars.drawing_parameters(system), nr.EQUATOR, areal=True)
        nr.DIAGRAMS.append(row)
        nr.CAPTIONS[(METRIC, system, "radial")] = [
            "The radial plane of $t$ and $r$ ($\\theta=\\pi/2$, $\\phi=0$), from the regular centre "
            "to the surface of zero pressure."]


def embedding(ck, src):
    """Read both sides of each matched equatorial slice and check their common tangent."""
    import embedding as eb
    views = []
    for system, name in zip(stars.SYSTEMS, stars.NAMES):
        inner = eb.Slice(src, METRIC, system, "r", "\\phi", {"t": 0, **eb.EQUATOR},
                         stars.drawing_parameters(system))
        outer = eb.Slice(src, "schwarzschild", "spherical", "r", "\\phi", {"t": 0, **eb.EQUATOR},
                         {"r_s": "2/5"})
        size, top = 6.0, 3.0
        star = eb.Piece("star", "star", inner, 0.0, 1.0, 0.0, 1,
                        (("axis", "the regular centre"), ("join", "$r=r_b$")),
                        [(0.5, "r", None), (1.0, "surface", "$r=r_b$")], size)
        ext = eb.Piece("exterior", "sheet", outer, 1.0, top, star.at(1.0)[1], 1,
                       (("join", "$r=r_b$"), ("edge", "$r\\to\\infty$")),
                       [(2.0, "r", None), (top, "r", None)], size)
        surface = eb.Surface([star, ext])
        ck.isometry(f"{name}: interior", star)
        ck.isometry(f"{name}: exterior", ext)
        ck.join(f"{name}: common surface tangent", star, 1.0, ext, 1.0)
        fig = eb.figure_of([surface], {"star": "star", "sheet": "cover"}, size, eb.Camera(-90, 32))
        eb.ring_label(fig, [0, 0, 0], *star.at(1.0), "$r_b$")
        fig.legend("fill", "star", "fluid interior")
        fig.legend("fill", "cover", "Schwarzschild exterior")
        fig.legend("line", "surface", "$r=r_b$")
        fig.legend("line", "r", "$r/r_b=0.5,2,3$")
        fig.legend("line", "meridian", "$\\phi$ constant")
        views.append(eb.view(system, name, "$r_b$", [surface], fig.done(), settings=SETTINGS))
    return views


def register_embedding(eb):
    eb.DRAWN[METRIC] = embedding
    for system in stars.SYSTEMS:
        eb.CAPTIONS[(METRIC, system)] = [
            "The equatorial plane ($t=0$, $\\theta=\\pi/2$), with the fluid bowl meeting "
            "Schwarzschild's paraboloid at one tangent plane."]


def register_slices(sl):
    for system in stars.SYSTEMS:
        sl.FLAT[(METRIC, system, "radial")] = lambda s=system: sl.one(
            METRIC, lambda m, s=s: sl.along(0.0, *m.reach(s, "r")), view_id=s)
    sl.FLAT_METRICS.add(METRIC)


def conformal(ck, src):
    """r*=integral dr/(y sqrt(Z)) inside, continued into the vacuum exterior."""
    import conformal as cd
    from scipy.integrate import cumulative_simpson
    class StellarPlane(cd.Plane):
        def prep(self, expr, limit=True):
            # Numeric comparisons need no global simplification of these radical lapses.
            e = cd.sp.sympify(expr).subs(self.reader.held).doit()
            return e.subs(self.subs).subs(self.fixed)
    views = []
    for system, name in zip(stars.SYSTEMS, stars.NAMES):
        inner = StellarPlane(src, METRIC, system, ("t", "r"), cd.EQUATOR, stars.drawing_parameters(system))
        outer = cd.Plane(src, "schwarzschild", "spherical", ("t", "r"), cd.EQUATOR, {"r_s": "2/5"})
        ck.limit(f"{name}: continuous lapse", [float(inner.g[0, 0].subs(inner.x1, 1))],
                 [float(outer.g[0, 0].subs(outer.x1, 1))], 1e-12)
        ck.limit(f"{name}: continuous radial metric", [float(inner.g[1, 1].subs(inner.x1, 1))],
                 [float(outer.g[1, 1].subs(outer.x1, 1))], 1e-12)
        speed = cd.sp.lambdify(inner.x1, cd.sp.sqrt(inner.g[1, 1]/(-inner.g[0, 0])), "numpy")
        radii = np.linspace(0, 1, 20001)
        optical = cumulative_simpson(speed(radii), x=radii, initial=0)
        def rstar(r):
            rr = np.asarray(r, dtype=float)
            # Clamp only the unused exterior expression inside the star.
            safe = np.maximum(rr, 1.0)
            exterior = optical[-1]+safe-1+0.4*np.log((safe-0.4)/0.6)
            return np.where(rr <= 1, np.interp(rr, radii, optical), exterior)
        def pq(t, r):
            return cd.mink_pq(t, rstar(r), 2.0)
        ck.chart(f"{name}: fluid", inner, pq, ck.uniform(-8, 8), ck.uniform(0.01, 0.99), lambda t, r: (1, 0))
        ck.chart(f"{name}: vacuum", outer, pq, ck.uniform(-8, 8), ck.uniform(1.01, 20), lambda t, r: (1, 0))
        import curvature as cv
        params = stars.drawing_parameters(system)
        scalar = cv.nr_published(src, METRIC, system, "kretschmann", ["r"], list(params))
        ck.finite(f"{name}: regular centre", float(scalar.exact([1e-5]+[float(v) for v in params.values()])))
        ck.limit(f"{name}: spatial infinity", cd.point(*pq(0, 1e10)), [math.pi, 0], 1e-4)
        v = cd.View(system, name, [-0.35, math.pi+0.35, -math.pi-0.25, math.pi+0.25], system)
        v.fill("region", cd.TRIANGLE)
        times = cd.S_ALL*2.0
        ps, qs = pq(times, np.ones_like(times))
        v.fill("cover", [[0, -math.pi]]+[cd.point(p, q) for p, q in zip(ps, qs)]+[[0, math.pi]])
        for r in (0.5, 2.0, 4.0):
            v.curve("r", *pq(times, np.full_like(times, r)))
        rr = np.concatenate([np.linspace(0, 1, 60)[:-1], 1+np.geomspace(1e-6, 1e5, 250)])
        for t in (-8, -4, 0, 4, 8):
            v.curve("t", *pq(np.full_like(rr, t), rr))
        v.curve("surface", ps, qs)
        cd.triangle_edges(v)
        cd.label_on(v, pq(0, 1), "$r=r_b$")
        v.legend("cover", "fluid interior")
        v.legend("r", "$r$ constant")
        v.legend("t", "$t$ constant")
        v.legend("surface", "$r=r_b$")
        v.legend("centre", "regular centre")
        v.set(settings=SETTINGS)
        moment = cd.slices.moments(METRIC, view_id=system)[0]
        lo, hi = moment.reach(coordinate="r")
        rr = np.linspace(lo, hi, 200)
        v.slice(moment, [pq(0*rr, rr)])
        views.append(v)
    return views


def register_conformal(cd):
    cd.DRAWN[METRIC] = conformal
    for system in stars.SYSTEMS:
        cd.CAPTIONS[(METRIC, system)] = [
            "The matched star and its vacuum exterior, with a timelike surface and a regular centre; "
            "each point a 2-sphere."]


def curvature(ck, src):
    """Only read the stored scalars, with no curvature calculation or extra matter audit."""
    import curvature as cv
    views = []
    for system, name in zip(stars.SYSTEMS, stars.NAMES):
        params = {k: float(v) for k, v in stars.drawing_parameters(system).items()}
        xaxis = cv.axis("$r/r_b$", 0, 1)
        curves, panels = [], []
        for field, letter, power in (("kretschmann", "K", 4), ("ricci_scalar", "R", 2)):
            fn = StableScalar(cv.nr_published(src, METRIC, system, field, ["r"], list(params)))
            # A finite centre is approached at 1e-5; exact decimal evaluation avoids cancellation.
            samples = [float(fn.exact([r]+list(params.values()))) for r in (1e-5, 0.1, 0.3, 0.6, 1.0)]
            low, high = min(0.0, min(samples)), max(samples)
            scale = {"kind": "linear", "from": low-0.1*max(abs(low), high), "to": 1.2*high}
            panels.append({"id": letter, "name": f"${letter}\\,r_b^{power}$", "aspect": 2,
                           "scale": scale, "ticks": cv.value_ticks(scale)})
            curves += cv.line_curves(ck, name, letter, letter, fn, "r", xaxis, scale, cv.PLOT_W/2,
                                     params, (1e-5, 1.0))
        cv.shape_checks(ck, name, [2, 2])
        views.append({"id": system, "label": name, "system": system, "kind": "line", "x": xaxis,
                      "panels": panels, "curves": curves, "marks": [],
                      "legend": [["line", "K", "$K$"], ["line", "R", "$R$"]],
                      "caption": [f"The Kretschmann scalar $K$ and the Ricci scalar $R$ of {name}, from the centre to the fluid surface."],
                      "ricci_sentences": [], "settings": SETTINGS,
                      "clip": "Linear scales, every value within its panel.",
                      "per_scalar": cv.scalar_clips("Linear, none clipped.", "Linear, none clipped.")})
    return views


def register_curvature(cv):
    cv.DRAWN[METRIC] = curvature
