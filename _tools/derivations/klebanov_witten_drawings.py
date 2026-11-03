"""Induced radial planes and a cylindrical spatial surface of the conifold throat."""
import math
import numpy as np

METRIC = "klebanov_witten"
FIXED = {"x": "0", "y": "0", "z": "0", "theta": "pi/2", "phi": "0",
         "alpha": "pi/2", "beta": "0", "psi": "0"}
LO, HI = 0.5, 2.5


def register_null(nr):
    nr.DIAGRAMS.append(nr.Diagram(
        METRIC, "poincare", "radial", "Radial plane", ("t", "r"), (0, 3, -1.5, 1.5),
        "$r/L$", "$ct/L$", {"L": 1}, FIXED, step=0.002, cones=(5, 5)))
    nr.CAPTIONS[(METRIC, "poincare", "radial")] = [
        r"Radial null geodesics ($x=y=z=\phi=\beta=\psi=0$, $\theta=\alpha=\pi/2$), with $dr/d(ct)=\pm r^2/L^2$."]
    nr.CLOSED_FORMS[(METRIC, "poincare", "radial")] = (
        lambda t, r: t-1/r, lambda t, r: t+1/r, lambda t, r: r>0.05)


def embedding(ck, src, eb):
    fixed = {k: v for k, v in FIXED.items() if k != "phi"}
    sl = eb.Slice(src, METRIC, "poincare", "r", r"\phi", {"t": 0, **fixed}, {"L": 1})
    size = 2.0
    tube = eb.Piece("cylinder", "sheet", sl, LO, HI, math.log(LO), 1,
                    (("edge", r"Continues toward $r=0$."), ("edge", r"Continues toward $r\to\infty$.")),
                    [(r, "r", None) for r in (LO, 1, 1.5, 2, HI)], size)
    ck.isometry("Klebanov-Witten: cylindrical spatial surface", tube)
    ck.radius("Klebanov-Witten: circle radius L/sqrt(6)", tube, lambda r: np.ones_like(r)/math.sqrt(6), size)
    ck.form("Klebanov-Witten: height L log(r/L)", tube, np.log, size)
    surface = eb.Surface([tube])
    fig = eb.figure_of([surface], {"sheet": "cover"}, size)
    eb.ring_label(fig, [0, 0, 0], *tube.at(LO), "$r=L/2$")
    eb.ring_label(fig, [0, 0, 0], *tube.at(HI), "$r=5L/2$")
    fig.legend("fill", "cover", r"$t=0$, $\theta=\alpha=\pi/2$")
    fig.legend("line", "r", "$r$ constant")
    fig.legend("line", "meridian", r"$\phi$ constant")
    return [eb.view("throat", "A circle down the throat", "$L$", [surface], fig.done(), settings="$L=1$.")]


def register_embedding(eb):
    eb.DRAWN[METRIC] = lambda ck, src: embedding(ck, src, eb)
    eb.CAPTIONS[(METRIC, "throat")] = [
        r"A spatial cylinder of radius $L/\sqrt6$ ($x=y=z=\beta=\psi=0$), with height $L\ln(r/L)$."]


def register_slices(sl):
    sl.FLAT[(METRIC, "poincare", "radial")] = lambda: sl.one(
        METRIC, lambda m: sl.along(0, *m.reach("poincare", "r")), view_id="throat")
    sl.FLAT_METRICS.add(METRIC)


def conformal(ck, src, cf):
    pl = cf.Plane(src, METRIC, "poincare", ("t", "r"), FIXED, {"L": 1})
    pq = lambda t, r: cf.poincare_pq(t, 1/np.asarray(r, dtype=float))
    ck.chart("Klebanov-Witten: radial AdS2 slice", pl, pq, ck.uniform(-5, 5), ck.uniform(0.05, 5), lambda t, r: (1, 0))
    ck.limit("Klebanov-Witten: finite throat curvature", pl.kretschmann(np.zeros(3), np.array([0.001, 1, 1000])), [176]*3, 1e-10)
    pi, half = math.pi, math.pi/2
    v = cf.View("poincare", "Radial AdS slice", [-half-0.4, half+0.4, -pi-0.4, pi+0.4], "poincare")
    v.fill("cover", [[-half, 0], [half, -pi], [half, pi]])
    v.line("outline", [[[half, -pi], [half, pi]]])
    v.line("chartedge", [[[-half, 0], [half, -pi]], [[-half, 0], [half, pi]]])
    for r in (LO, 1, 2, HI):
        ts = np.linspace(-1000, 1000, 2001)
        v.curve("r", *pq(ts, np.full_like(ts, r)))
    for t in (-2, -1, 0, 1, 2):
        rr = np.geomspace(0.001, 1000, 501)
        v.curve("meridian", *pq(np.full_like(rr, t), rr))
    v.legend("cover", "Poincaré patch")
    v.legend("outline", "Timelike boundary")
    v.legend("chartedge", "Poincaré horizon")
    v.legend("r", "$r$ constant")
    v.legend("meridian", "$t$ constant")
    moment = cf.slices.moments(METRIC, "throat")[0]
    rr = np.linspace(*moment.reach("poincare", "r"), 101)
    v.slice(moment, [pq(np.zeros_like(rr), rr)])
    v.set(settings="$L=1$.")
    return [v]


def register_conformal(cf):
    cf.DRAWN[METRIC] = lambda ck, src: conformal(ck, src, cf)
    cf.CAPTIONS[(METRIC, "poincare")] = [
        r"The radial anti-de Sitter slice ($x=y=z=\phi=\beta=\psi=0$, $\theta=\alpha=\pi/2$), compactified through its null coordinates."]
