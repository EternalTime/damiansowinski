#!/usr/bin/env python3
"""Where each embedding diagram is cut from, on its spacetime's other diagrams.

Every embedding diagram draws part of one moment of its spacetime, or of several moments in
turn, and embedding_slices.md beside this file works out, spacetime by spacetime, how that
moment meets every spacetime diagram and conformal diagram drawn of it. This module is that
table in code. It reads the embedding files themselves, so a moment is drawn where the
embedding stands and over the part of it the embedding reaches, and it says what the moment
is in the chart of each drawing.

    Moment      one surface of one embedding view: its moment's label, its time, its reach
                along each coordinate the embedding read it in, and the stamp over them all
                that build_mfs_data.py checks every drawn slice against;
    flat(spec)  the moments one flat view of null_rays.py shows, each as polylines, points
                or regions in the chart coordinates (x^0, r) of that view's plane, which the
                view carries into its unit square with the map it draws everything else with;
    HIDDEN      the flat views and figures on which a moment of their spacetime lies but is
                not drawn, each with the reason in the physics.

projections.py and conformal.py draw their slices inside the functions that draw each
figure and each conformal view, with the maps those functions already use, from the
Moments and reaches read here.

A drawing on which no moment of its spacetime's embedding lies carries no slice: the flat
universe of FRW and the marginally bound cloud of Tolman-Bondi are other spacetimes than
the ones embedded, and Vaidya's outgoing chart draws the exploding shell, the time reverse
of the one embedded.

The transformations a moment is carried through, where the drawing's chart is not the
embedding's, are each checked here by pulling the published metric of one chart back onto
the other, at points of the drawn plane: `python3 slices.py` runs those checks and prints
them, and `null_rays.py --slices` rewrites only the slices of every diagram file.
"""

import functools
import json
import math
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import build_mfs_data as build  # noqa: E402

BIG = 1e6                       # a line of a homogeneous moment runs across any box
N = 801                         # points along a curved moment before it is thinned


@functools.lru_cache(maxsize=None)
def embedding(metric_id):
    path = build.EMBEDDING_DIR / f"{metric_id}.json"
    return json.loads(path.read_text(encoding="utf-8"))


class Moment:
    """One surface of one embedding view, or one ring of a stack of moments, and the moment of the
    spacetime it is cut from."""

    def __init__(self, metric_id, view, index, label=None, curve=None):
        self.metric, self.view, self.index, self.curve = metric_id, view, index, curve
        self.surface = view["surfaces"][index]
        if curve is not None:
            # A stack names each moment on the ring it marks at the height of its time.
            ring = self.surface["curves"][curve]
            self.time, self.label = ring["time"], ring["label"]
        else:
            self.time = self.surface.get("time")
            # A sequence names each moment on its surface; a single moment is drawn at t = 0 of
            # the chart the embedding read it in, save where the view names another.
            self.label = self.surface.get("label") or label or "$t = 0$"
        self.version = build.embedding_moment_version(view, index, curve)

    def reach(self, system=None, coordinate=None, reference=False):
        """The least and greatest value of the coordinate over the pieces of the surface read
        in `system`, reference pieces left out unless asked for, as (lo, hi)."""
        lo, hi = math.inf, -math.inf
        for piece in self.surface["pieces"]:
            if piece.get("reference") and not reference:
                continue
            if system and piece.get("system") != system:
                continue
            if coordinate and piece.get("coordinate") != coordinate:
                continue
            if "grid" in piece:
                raise ValueError(f"{self.metric}: a grid piece has no reach along a coordinate")
            x = [p[0] for p in piece["points"]]
            lo, hi = min(lo, *x), max(hi, *x)
        if lo > hi:
            raise ValueError(f"{self.metric}/{self.view['id']}: no piece read in {system} along {coordinate}")
        return lo, hi

    def grid(self):
        """The one grid piece of a height over a plane."""
        pieces = [p for p in self.surface["pieces"] if "grid" in p]
        if len(pieces) != 1:
            raise ValueError(f"{self.metric}/{self.view['id']}: not one grid piece")
        return pieces[0]["grid"]

    def json(self):
        out = {"view": self.view["id"], "surface": self.index}
        if self.curve is not None:
            out["curve"] = self.curve
        return dict(out, label=self.label, version=self.version)


def moments(metric_id, view_id=None, label=None):
    """Every moment of a spacetime's embedding views, in the order of views and surfaces, a
    stack's in the order of the rings it marks at their times."""
    out = []
    for view in embedding(metric_id)["views"]:
        if view_id is None or view["id"] == view_id:
            for i, surface in enumerate(view["surfaces"]):
                rings = [k for k, c in enumerate(surface.get("curves", [])) if "time" in c]
                out += [Moment(metric_id, view, i, label, k) for k in rings] or [Moment(metric_id, view, i, label)]
    if not out:
        raise ValueError(f"{metric_id}: no embedding view {view_id!r}")
    return out


class Mark:
    """A moment on one drawing: polylines, points and regions, each region a list of rings
    filled by the even odd rule, in whatever coordinates the drawing hands it."""

    def __init__(self, moment, lines=(), points=(), fills=(), label=None):
        self.moment = moment
        self.lines = [np.asarray(line, dtype=float) for line in lines]
        self.points = [np.asarray(p, dtype=float) for p in points]
        self.fills = [[np.asarray(ring, dtype=float) for ring in rings] for rings in fills]
        self.label = label or moment.label


# ---------------------------------------------------------------- the moments in the flat views

def across(t, lo, hi):
    """A moment of constant time t on a plane through the centre of a surface of revolution,
    whose drawn coordinate runs through the axis from -hi to hi: both profiles, and a gap
    about the axis where the embedding stops short of it."""
    if lo <= 0:
        return [[(t, -hi), (t, hi)]]
    return [[(t, -hi), (t, -lo)], [(t, lo), (t, hi)]]


def along(t, lo, hi):
    """A moment of constant time t on a plane of a radius, from lo to hi."""
    return [[(t, lo), (t, hi)]]


def near(lo, hi, n=N, crowd=1e-9):
    """Values on (lo, hi] crowding toward lo, where a curve through a horizon runs off."""
    return lo + (hi - lo) * np.geomspace(crowd, 1.0, n)


def gravastar_x(r):
    """The tortoise coordinate inside the gravastar the diagrams draw, L = 2 and C = 64/195."""
    return 2 / math.sqrt(64 / 195) * math.atanh(r / 2)


def schwarzschild_t(sign):
    """Schwarzschild's t = 0 in an Eddington-Finkelstein chart, r_s = 1: v = r + ln(r - 1)
    in the ingoing chart and u = -r - ln(r - 1) in the outgoing one, outside r_s."""
    m = moments("schwarzschild")[0]
    lo, hi = m.reach("spherical", "r")
    r = near(lo, hi)
    return [Mark(m, [np.column_stack([sign * (r + np.log(r - 1)), r])])]


def kottler_rstar(r):
    """Kottler's tortoise coordinate at r_s = 1 and Lambda = 1/5, as the Eddington-Finkelstein
    charts fix it, r_* = sum_i ln|1 - r/r_i|/f'(r_i) over the three roots of 3r - 3 - r^3/5,
    which vanishes at r = 0."""
    roots = np.roots([-0.2, 0, 3, -3]).real
    return sum(np.log(np.abs(1 - r / ri)) / (1 / ri ** 2 - 0.4 * ri / 3) for ri in roots)


def kottler_t(sign):
    """Kottler's static t = 0 in an Eddington-Finkelstein chart: v = r_* in the ingoing chart and
    u = -r_* in the outgoing one, between the horizons, as far as the embedding reaches."""
    m = moments("schwarzschild_de_sitter")[0]
    lo, hi = m.reach("static", "r")
    r = np.concatenate([near(lo, 0.5 * (lo + hi)), near(hi, 0.5 * (lo + hi))[::-1]])
    return [Mark(m, [np.column_stack([sign * kottler_rstar(r), r])])]


KISELEV_ROOTS = (4 - 2 * math.sqrt(2), 4 + 2 * math.sqrt(2))


def kiselev_rstar(r):
    """Kiselev's tortoise coordinate at w = -2/3, r_s = 1 and r_q = 8, as the Eddington-Finkelstein
    charts fix it: 1/f = -8r/((r - a)(r - b)) with a, b = 4 -+ 2 sqrt 2, so
    r_* = (4 sqrt 2 - 4) ln|1 - r/a| - (4 sqrt 2 + 4) ln|1 - r/b|, which vanishes at r = 0."""
    a, b = KISELEV_ROOTS
    return ((4 * math.sqrt(2) - 4) * np.log(np.abs(1 - r / a))
            - (4 * math.sqrt(2) + 4) * np.log(np.abs(1 - r / b)))


def kiselev_t(sign):
    """Kiselev's static t = 0 in an Eddington-Finkelstein chart: v = r_* in the ingoing chart and
    u = -r_* in the outgoing one, between the horizons, as far as the embedding reaches."""
    m = moments("kiselev", "black_hole")[0]
    lo, hi = m.reach("linear", "r")
    r = np.concatenate([near(lo, 0.5 * (lo + hi)), near(hi, 0.5 * (lo + hi))[::-1]])
    return [Mark(m, [np.column_stack([sign * kiselev_rstar(r), r])])]


def kiselev_free(flat):
    """The static t = 0 of Kiselev's matter alone in his two charts of it. It is eta = 0, the
    whole half line of chi, since r = r_q(1 - e^(-2 chi)) reaches the horizon only as chi grows
    without bound, and in the conformally flat chart, tau = e^eta cosh(chi) and
    rho = e^eta sinh(chi), the hyperbola tau = sqrt(1 + rho^2)."""
    m = moments("kiselev", "free")[0]
    if not flat:
        return [Mark(m, along(0.0, 0.0, BIG), label="$\\eta = 0$")]
    rho = np.linspace(0.0, BIG ** 0.25, 4 * N + 1)
    return [Mark(m, [np.column_stack([np.sqrt(1 + rho * rho), rho])], label="static $t = 0$")]


def monopole_rstar(r):
    """The tortoise coordinate of Letelier's black hole at Delta = 0.19 and r_s = 1, as the
    Eddington-Finkelstein charts fix it, r_* = r/(1 - Delta) + ln|(1 - Delta)r - 1|/(1 - Delta)^2,
    which vanishes at r = 0."""
    return r / 0.81 + np.log(np.abs(0.81 * r - 1)) / 0.81 ** 2


def monopole_t(sign):
    """Letelier's black hole's static t = 0 in an Eddington-Finkelstein chart: v = r_* in the
    ingoing chart and u = -r_* in the outgoing one, outside r_h, as far as the embedding reaches.
    The monopole's cone is another spacetime, with no mass at its centre, and is not drawn here."""
    m = moments("global_monopole", "black_hole")[0]
    lo, hi = m.reach("static", "r")
    r = near(lo, hi)
    return [Mark(m, [np.column_stack([sign * monopole_rstar(r), r])])]

def string_hole(sign=0):
    """The threaded black hole's two moments, r_s = 1. On the static planes, the equator's t = 0
    over the r it reaches and the horizon's bifurcation sphere, the point t = 0, r = r_s. In an
    Eddington-Finkelstein chart the equator's t = 0 is v = r + ln(r - 1) or u = -r - ln(r - 1),
    Schwarzschild's, and the bifurcation sphere lies at v -> -infinity or u -> +infinity, off the chart."""
    equator = moments("string_black_hole", "equator")[0]
    lo, hi = equator.reach("static", "r")
    if not sign:
        horizon = moments("string_black_hole", "horizon", label="$t = 0$, $r = r_s$")[0]
        return [Mark(equator, along(0.0, lo, hi)), Mark(horizon, points=[(0.0, 1.0)])]
    r = near(lo, hi)
    return [Mark(equator, [np.column_stack([sign * (r + np.log(r - 1)), r])])]


def teo_proper(r):
    """Teo's proper radial distance from the throat at b_0 = 1, his eq. (28)."""
    return math.sqrt(r * (r - 1)) + math.log(math.sqrt(r) + math.sqrt(r - 1))


def teo(system, view):
    """Teo's wormhole's two moments. The equatorial plane's t = 0, embedded at a = 1, lies on the
    equator's planes, over the r it reaches on either side of the throat. The throat r = b_0 at
    t = 0, embedded at a = 1/4, meets the axis at its poles, one point of the axis's planes, drawn
    at that spin. Neither lies on the other's drawings: the equatorial plane does not meet the axis."""
    if view == "axis":
        at, label = ((0.0, 1.0), "$t = 0$, $r = b_0$") if system == "spherical" else ((0.0, 0.0), "$t = 0$, $l = 0$")
        return [Mark(moments("teo_wormhole", "throat", label=label)[0], points=[at])]
    equator = moments("teo_wormhole", "equator")[0]
    lo, hi = equator.reach("spherical", "r")
    if system == "spherical":
        return [Mark(equator, along(0.0, lo, hi))]
    return [Mark(equator, along(0.0, -teo_proper(hi), teo_proper(hi)))]


def dilaton_t(sign):
    """The dilaton black hole's static t = 0 in an Eddington-Finkelstein chart, r_s = 1: its plane
    of t and r is Schwarzschild's, so v = r + ln(r - 1) in the ingoing chart and u = -r - ln(r - 1)
    in the outgoing one, outside r_s, as far as the Einstein metric's embedding reaches. The string
    metrics' moments are the same events measured by other metrics, and are marked on their own charts."""
    m = moments("dilaton_black_hole", "einstein")[0]
    lo, hi = m.reach("static", "r")
    r = near(lo, hi)
    return [Mark(m, [np.column_stack([sign * (r + np.log(r - 1)), r])])]


def kkbh_smooth(r, a=1.0):
    """The part of the Kaluza-Klein black holes' tortoise coordinate that is smooth through the
    horizon, at r_s = 1 with a = q - r_s or p - r_s: r_* less sqrt(1 + a) ln|r - 1|, where
    r_* = s + (1 + a/2) ln((2s + 2r + a)/a) - sqrt(1 + a) ln(((a + 2) r + a + 2 sqrt(1 + a) s)/(a |r - 1|))
    with s = sqrt(r(r + a)) is the integral of sqrt(r(r + a))/(r - 1) from the singularity r = 0."""
    r = np.maximum(np.asarray(r, dtype=float), 0.0)
    s = np.sqrt(r * (r + a))
    k = math.sqrt(1 + a)
    return s + (1 + a / 2) * np.log((2 * s + 2 * r + a) / a) - k * np.log(((a + 2) * r + a + 2 * k * s) / a)


def kkbh_rstar(r, a=1.0):
    """The tortoise coordinate of a Kaluza-Klein black hole of one charge at r_s = 1, zero at r = 0:
    dr_*/dr = sqrt(r(r + a))/(r - 1), the rays with no momentum along the circle and the rays of
    the Einstein metric of four dimensions alike."""
    with np.errstate(divide="ignore"):
        return kkbh_smooth(r, a) + math.sqrt(1 + a) * np.log(np.abs(np.asarray(r, dtype=float) - 1))


def kkbh_t(view_id):
    """The static t = 0 of a Kaluza-Klein black hole of one charge, q = 2 r_s, in its ingoing
    charts, outside r_s, as far as the embedding reaches: v = sqrt(q/r_s) (r + r_s ln(r/r_s - 1)) in
    the chart of five dimensions, whose v is built on Schwarzschild's tortoise coordinate, and
    v = r_* in the Einstein metric's chart."""
    system = {"electric": "electric", "einstein": "einstein"}[view_id]
    m = moments("kaluza_klein_black_hole", view_id)[0]
    lo, hi = m.reach(system, "r")
    r = near(lo, hi)
    v = math.sqrt(2) * (r + np.log(r - 1)) if view_id == "electric" else kkbh_rstar(r)
    return [Mark(m, [np.column_stack([v, r])])]


def tangherlini_t(sign):
    """Tangherlini's static t = 0 in five dimensions in an Eddington-Finkelstein chart, r_h = 1:
    v = r_* in the ingoing chart and u = -r_* in the outgoing one, r_* = r + ln((r - 1)/(r + 1))/2,
    outside r_h, as far as the embedding reaches."""
    m = moments("tangherlini", "five")[0]
    lo, hi = m.reach("spherical", "r")
    r = near(lo, hi)
    return [Mark(m, [np.column_stack([sign * (r + 0.5 * np.log((r - 1) / (r + 1))), r])])]


# Boulware and Deser's black hole as every diagram draws it, in units of its horizon radius: the
# mass radius r_0 = 13/12 and the Gauss-Bonnet length l = 5/12, so that r_h^2 = r_0^2 - l^2 = 1; and
# the other branch in units of l, at r_0 = l. The black hole's ratio r_0/l = 13/5 is above
# 1 + sqrt 2, below which Beroiz, Dotti and Gleiser found the black hole unstable (Phys. Rev. D
# 76, 024012).
BD_R0, BD_ELL = 13 / 12, 5 / 12
BD_PLUS_R0 = 1.0


def _gauss(of, a, b, panels):
    """The integral of a smooth function from a to b by Gauss and Legendre's rule of twelve points
    on each of `panels` equal panels."""
    nodes, weights = np.polynomial.legendre.leggauss(12)
    edges = np.linspace(a, b, panels + 1)
    half, mid = 0.5 * np.diff(edges), 0.5 * (edges[:-1] + edges[1:])
    return float(np.sum(half[:, None] * weights * of(mid[:, None] + half[:, None] * nodes)))


def _bd_part(r, of, far):
    """The integral from 0 to r of `of`, which is smooth and falls as far/s^2: in s out to
    s = 1 and in 1/s beyond it, where the integrand of[1/w]/w^2 is smooth down to w = 0."""
    def one(x):
        if x <= 1.0:
            return _gauss(of, 0.0, x, 8) if x > 0 else 0.0
        inner = _gauss(of, 0.0, 1.0, 8)
        w = 0.0 if math.isinf(x) else 1.0 / x
        return inner + _gauss(lambda q: np.where(q > 0, of(1 / np.maximum(q, 1e-300)) / np.maximum(q, 1e-300) ** 2, far),
                              w, 1.0, 16)
    return np.array([one(float(x)) for x in np.atleast_1d(np.asarray(r, dtype=float))]).reshape(np.shape(r))


def boulware_deser_rstar(r):
    """The tortoise coordinate of Boulware and Deser's black hole at r_0 = 13/12 and l = 5/12, as the
    Eddington-Finkelstein charts fix it, vanishing at r = 0. With W = sqrt(r^4 + 4 l^2 r_0^2) and
    W_h = r_h^2 + 2 l^2 its value at the horizon r_h = 1,
    1/f = (r^2 + 2 l^2 + W)/(2 (r^2 - r_h^2)) = 1 + W_h/(r^2 - r_h^2) - (W - r^2 + 2 l^2)/(2 (W + W_h)),
    so r_* = r + (W_h/2 r_h) ln|(r - r_h)/(r + r_h)| - J(r)/2: Tangherlini's, with the surface
    gravity r_h/W_h in place of 1/r_h, less half the integral J of (W - s^2 + 2 l^2)/(W + W_h)
    from 0, which is smooth, falls as 2 l^2/s^2 and vanishes at l = 0."""
    l2, k = BD_ELL ** 2, 4 * BD_ELL ** 2 * BD_R0 ** 2
    wh = 1 + 2 * l2

    def of(s):
        w = np.sqrt(s ** 4 + k)
        # W - s^2 is written k/(W + s^2), which keeps its digits far out.
        return (k / (w + s * s) + 2 * l2) / (w + wh)
    r = np.asarray(r, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        pole = np.where(np.isinf(r), 0.0, 0.5 * wh * np.log(np.abs((r - 1) / (r + 1))))
    return r + pole - 0.5 * _bd_part(r, of, 2 * l2)


def boulware_deser_plus_rstar(r):
    """The tortoise coordinate of the other branch at l = 1 and r_0 = 1, vanishing at r = 0:
    the integral of 1/f_+ = 2 l^2/(r^2 + 2 l^2 + W), which falls as l^2/s^2, so r_* tends to a
    finite value R as r -> infinity."""
    k = 4 * BD_PLUS_R0 ** 2
    return _bd_part(r, lambda s: 2 / (s * s + 2 + np.sqrt(s ** 4 + k)), 1.0)


def boulware_deser_t(sign):
    """The static t = 0 of Boulware and Deser's black hole in an Eddington-Finkelstein chart, r_h = 1:
    v = r_* in the ingoing chart and u = -r_* in the outgoing one, outside r_h, as far as the
    embedding reaches."""
    m = moments("boulware_deser", "hole")[0]
    lo, hi = m.reach("spherical", "r")
    r = near(lo, hi)
    return [Mark(m, [np.column_stack([sign * boulware_deser_rstar(r), r])])]


def black_string_t(kerr_schild=False):
    """The black string's static t = 0 across the string, r_s = 1, where the plane of the time and r
    is Schwarzschild's: v = r + ln(r - 1) in the ingoing Eddington-Finkelstein chart, and
    cT = v - r = ln(r - 1) in the Kerr-Schild chart, outside r_s, as far as the embedding reaches.
    The rippled horizon is the perturbed string, another spacetime, and is marked on no drawing."""
    m = moments("black_string", "across")[0]
    lo, hi = m.reach("static", "r")
    r = near(lo, hi)
    return [Mark(m, [np.column_stack([np.log(r - 1) + (0 if kerr_schild else r), r])])]


def myers_perry_t(view_id):
    """Myers and Perry's Boyer-Lindquist t = 0 in the ingoing chart, one spin in five dimensions at
    mu = 1 and a = 3/5: v = r_*, r_* = r + (5/8) ln((r - 4/5)/(r + 4/5)), outside r+ = 4/5, the same on
    the plane transverse to the rotation and on the plane of rotation, as far as each embedding
    reaches."""
    m = moments("myers_perry", view_id)[0]
    lo, hi = m.reach("boyer_lindquist", "r")
    r = near(lo, hi)
    return [Mark(m, [np.column_stack([r + 0.625 * np.log((r - 0.8) / (r + 0.8)), r])])]


def black_saturn_t(view):
    """Black Saturn's moment t = 0 on the plane of its ring, as far as its embedding reaches: the piece
    outside the ring on Weyl's plane of t and z there, the piece about the hole on the plane between the
    ring and the hole, and the outside piece's part with z < 0 on the polar chart's theta = pi/2, where
    z = -r^2/2."""
    m = moments("black_saturn", "plane")[0]
    pieces = {p["id"]: [q[0] for q in p["points"]] for p in m.surface["pieces"]}
    if view == "far":
        return [Mark(m, along(0.0, 0.0, math.sqrt(-2 * min(pieces["outside"]))))]
    z = pieces["outside" if view == "outside" else "between"]
    return [Mark(m, along(0.0, min(z), max(z)))]


def novikov(R, tau):
    """A shell of dust released from rest at areal radius R at t = 0, r_s = 1, at its proper
    time tau: its areal radius r and its Schwarzschild t, from the cycloid
    r = (R/2)(1 + cos eta), tau = (R/2) sqrt(R) (eta + sin eta), and Misner, Thorne and
    Wheeler's (31.10), t = ln|(k + tan(eta/2))/(k - tan(eta/2))| + k (eta + (R/2)(eta + sin eta))
    with k = sqrt(R - 1). The logarithm diverges where the shell crosses r_s."""
    R = np.asarray(R, dtype=float)
    eta = np.empty_like(R)
    for i, Ri in enumerate(R):
        scale = 0.5 * Ri * math.sqrt(Ri)
        lo, hi = 0.0, math.pi
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if scale * (mid + math.sin(mid)) < tau else (lo, mid)
        eta[i] = 0.5 * (lo + hi)
    k = np.sqrt(R - 1)
    tan = np.tan(eta / 2)
    r = 0.5 * R * (1 + np.cos(eta))
    with np.errstate(divide="ignore"):
        t = np.log(np.abs((k + tan) / (k - tan))) + k * (eta + 0.5 * R * (eta + np.sin(eta)))
    return r, t, eta


def novikov_kruskal(R, tau):
    """The same shells in Kruskal's U and V, r_s = 1, regular through r_s:
    V = (k cos(eta/2) + sin(eta/2)) exp((r + k(eta + (R/2)(eta + sin eta)))/2) and
    U V = (1 - r) e^r."""
    r, _, eta = novikov(R, tau)
    R = np.asarray(R, dtype=float)
    k = np.sqrt(R - 1)
    V = (k * np.cos(eta / 2) + np.sin(eta / 2)) * np.exp((r + k * (eta + 0.5 * R * (eta + np.sin(eta)))) / 2)
    return (1 - r) * np.exp(r) / V, V, r


OS_R0 = 2.0                                  # the release, in r_s, as the conformal diagram declares
OS_AM = OS_R0 / math.sin(math.pi / 4)        # a_m = R_0 / sin chi_0 = 2 sqrt 2 r_s


def os_exterior():
    """Novikov's slice outside the star: the shells from the surface out as far as the
    embedding reaches, each at the dust's proper time, drawn where r > r_s, where
    Schwarzschild's chart ends."""
    out = []
    for m in moments("oppenheimer_snyder"):
        lo, hi = m.reach("comoving_synchronous", "r")
        R = near(lo, hi, 2001, 1e-6)
        r, t, _ = novikov(R, m.time)
        keep = (r > 1) & np.isfinite(t)
        out.append(Mark(m, [np.column_stack([t[keep], r[keep]])]))
    return out


def kerr_above(metric_id):
    """The equator of the moment seen from above with t left out: the whole plane outside
    the horizon, as far as the embedding reaches, as a region of (phi, r)."""
    m = moments(metric_id)[0]
    lo, hi = m.reach("boyer_lindquist", "r")
    phi = np.linspace(0, 2 * np.pi, 721)
    return [Mark(m, fills=[[np.column_stack([phi, np.full_like(phi, hi)]),
                            np.column_stack([phi, np.full_like(phi, lo)])]])]


def levi_civita_r(rho, sigma=0.25):
    """The proper distance from Levi-Civita's axis, the radius of the Kasner form:
    r = rho^Sigma/Sigma with Sigma = 4 sigma^2 - 2 sigma + 1."""
    Sigma = 4 * sigma * sigma - 2 * sigma + 1
    return rho ** Sigma / Sigma


def vdb_radius(l):
    """The comoving radius, in units of R, of the sphere at proper distance l from the middle of
    the neck of Van Den Broeck's declared pocket: d ln(rho) = dl/r(l), with rho = l outside the
    neck, the inverse of null_rays._vdb_distance."""
    neck = math.exp(-math.pi) / 2
    a = math.atanh(1 / math.sqrt(7))
    if l < -2.5:
        return (l + 4) / 1.5 * neck / 3 * math.exp(-4 * a / math.sqrt(7))
    if l < -1.5:
        return neck / 3 * math.exp(2 / math.sqrt(7) * (math.atanh(2 * (l + 2) / math.sqrt(7)) - a))
    if l < -0.5:
        return -math.exp(-math.pi) / (4 * l)
    if l < 0.5:
        return math.exp(2 * math.atan(2 * l) - math.pi / 2) / 2
    return l


MCV_H0 = 1 / math.sqrt(15)      # McVittie's H_0 r_s/c, as its diagrams and its embedding declare it


def mcvittie_areal(t, r):
    """The areal radius R = ar(1 + r_s/(4ar))^2 of McVittie's comoving r at ct, in r_s, with
    a = sinh^(2/3)(3 H_0 t/2)."""
    x = math.sinh(1.5 * MCV_H0 * t) ** (2 / 3) * r
    return x * (1 + 1 / (4 * x)) ** 2


def _mcvittie_areal(m):
    """A moment of McVittie's cosmic time in the areal chart, which keeps that time: R from the
    throat r_s out to the areal radius of the comoving r the embedding reaches."""
    lo, hi = m.reach("isotropic", "r")
    return along(m.time, mcvittie_areal(m.time, lo), mcvittie_areal(m.time, hi))


def _ds_isotropic(R):
    """The isotropic radius outside the throat of Damour and Solodukhin's sphere of areal radius R, at
    r_s = 1: R = r(1 + 1/4r)^2."""
    return (R - 0.5 + math.sqrt(R * (R - 1))) / 2


SV_A = {"bounce": 0.5, "null": 1.0, "wormhole": 2.0}      # a in r_s of each of Simpson and Visser's geometries
SV_MOMENT = {"bounce": "outside", "null": "null", "wormhole": "wormhole"}


def sv_rstar(r, a):
    """The tortoise coordinate of Simpson and Visser's black bounce at r_s = 1, as null_rays._sv_rstar."""
    r = np.asarray(r, float)
    rho = np.sqrt(r * r + a * a)
    out = r + np.arcsinh(r / a)
    if a < 1:
        h = math.sqrt(1 - a * a)
        return out + np.log(np.abs((r - h) * (r - h * rho) / ((r + h) * (r + h * rho)))) / (2 * h)
    if a == 1:
        return out - (1 + rho) / r
    k = math.sqrt(a * a - 1)
    return out + (np.arctan(r / k) + np.arctan(r / (k * rho))) / k


def sv_reach(case):
    """The moment t = 0 of one of Simpson and Visser's geometries and the r it reaches: the black
    bounce's exterior is read in the areal radius rho, where r = sqrt(rho^2 - a^2)."""
    m = moments("simpson_visser", SV_MOMENT[case])[0]
    if case == "bounce":
        return m, tuple(math.sqrt(rho * rho - 0.25) for rho in m.reach("areal", "\\rho"))
    return m, m.reach("spherical", "r")


def sv_inside_r(m):
    """The r of a moment of the black bounce between its horizons, at a = r_s/2: its cylinder's radius
    is sqrt(r^2 + a^2), and r is positive before the proper time of r = 0, the third moment's."""
    radius = m.surface["pieces"][0]["points"][0][1]
    turn = moments("simpson_visser", "inside")[2].time
    return math.copysign(math.sqrt(max(radius * radius - 0.25, 0.0)), turn - m.time) if radius > 0.5 else 0.0


def simpson_visser(chart, case):
    """The moments of Simpson and Visser's geometry `case` on its plane in `chart`, as (x^0, r). The
    moment t = 0 is level in the two charts of t, from r to rho = sqrt(r^2 + a^2) in the areal chart,
    and v = r_* or u = -r_* in the Eddington-Finkelstein charts. The black bounce's moments of
    constant r between the horizons are each the stretch |ct| <= r_s the embedding reaches, upright
    at r about t = 0 and about v = r_* in the ingoing chart, and at rho in the areal chart, which is
    drawn on the side r > 0 and so holds the moments down to r = 0."""
    a = SV_A[case]
    m, (lo, hi) = sv_reach(case)
    sign = {"eddington_finkelstein_ingoing": 1, "eddington_finkelstein_outgoing": -1}.get(chart, 0)
    if sign:
        r = near(lo, hi) if case == "bounce" else np.linspace(lo, hi, N)
        marks = [Mark(m, [np.column_stack([sign * sv_rstar(r, a), r])])]
    elif chart == "areal":
        # A moment through the throat reaches it from both sides, each half the line from rho = a out.
        nearest = 0.0 if lo < 0 < hi else min(abs(lo), abs(hi))
        marks = [Mark(m, along(0.0, math.sqrt(nearest ** 2 + a * a), math.sqrt(max(lo * lo, hi * hi) + a * a)))]
    else:
        marks = [Mark(m, along(0.0, lo, hi))]
    if case != "bounce" or sign < 0:
        # The outgoing chart's region between the horizons is the white hole, where the black hole's
        # moments do not lie.
        return marks
    for inside in moments("simpson_visser", "inside"):
        r = sv_inside_r(inside)
        if chart == "areal" and r < 0:
            # The areal chart is drawn on the side r > 0, up to the sphere of least area.
            continue
        t0, t1 = inside.reach("spherical", "t")
        x = math.sqrt(r * r + a * a) if chart == "areal" else r
        shift = sign * float(sv_rstar(r, a)) if sign else 0.0
        marks.append(Mark(inside, [[(t0 + shift, x), (t1 + shift, x)]]))
    return marks


# Bonnor and Vaidya's charged shell, in units of its mass M: its charge, the value of Reissner-Nordstrom's
# drawings, r_q = 0.48 r_s, and the roots r_+ and r_- of 1 - 2M/r + q^2/r^2.
BV_Q = 0.96
BV_RP, BV_RM = 1.28, 0.72


ROBERTS_P = {"disperses": 0.9, "threshold": 1.0, "collapses": 2.0}     # p of each of the Roberts solution's outcomes


def roberts_null(p, t, rho):
    """The null coordinates u and v of the point (t, rho) of Roberts's diagonal chart, at l = c = 1:
    u = sqrt(1 + p)(t - rho) and v = (t + rho)/sqrt(1 + p)."""
    s = math.sqrt(1 + p)
    return s * (t - rho), (t + rho) / s


def roberts_point(chart, p, t, rho):
    """The point (t, rho) of Roberts's diagonal chart as the (x^0, r) of another chart's plane: u and
    v; v and r = ((1 + p) v - u)/2; v and the areal radius sqrt(rho (rho - p t)); and the scaling
    chart's tau = -ln(-u/2) and x = ln(1 - 2v/u)/2."""
    rho = np.asarray(rho, dtype=float)
    u, v = roberts_null(p, t, rho)
    if chart == "double_null":
        return u, v
    if chart == "advanced":
        return v, ((1 + p) * v - u) / 2
    if chart == "areal":
        return v, np.sqrt(np.maximum(rho * (rho - p * t), 0.0))
    if chart == "diagonal":
        return np.full_like(rho, t), rho
    with np.errstate(divide="ignore", invalid="ignore"):
        return -np.log(-u / 2), np.log(1 - 2 * v / u) / 2


def roberts(chart, case):
    """Each moment of Roberts's time t that the embedding of `case` marks, over the stretch of rho it
    reaches in the diagonal chart, on the plane of `chart`."""
    p = ROBERTS_P[case]
    marks = []
    for m in moments("roberts", case):
        lo, hi = m.reach("diagonal", "\\rho")
        rho = near(lo, hi)
        marks.append(Mark(m, [np.column_stack(roberts_point(chart, p, m.time, rho))]))
    return marks


def _fjnw_reach(m, of_r):
    """The embedding of Fisher, Janis, Newman and Winicour's equator is read in the harmonic chart at
    k = 1/2, b = 1, where e^(-u) = 1 - b/r: the radii r it reaches, least first, each carried to a
    chart's own radial coordinate by of_r."""
    lo, hi = m.reach("harmonic", "u")
    return [of_r(1 / (1 - math.exp(-u))) for u in (hi, lo)]


def witten(chart):
    """The moment t = 0 of Witten's black hole outside the horizon, the meridian theta = 0 of his
    cigar, as far as the embedding reaches in his proper distance r, at lambda = m = 1: level in
    each chart of t, from the horizon out, with e^(2x) = w = cosh^2 r and e^sigma = sinh r; the
    curve v = sigma(x) or u = -sigma(x), sigma = ln(e^(2x) - 1)/2, in an Eddington-Finkelstein chart;
    and on the Kruskal plane the line U + V = 0 through the bifurcation point, V = -U = sinh r on
    one side and its mirror on the other, the meridian theta = pi."""
    m = moments("witten_black_hole")[0]
    lo, hi = m.reach("witten", "r")
    if chart == "witten":
        return [Mark(m, along(0.0, lo, hi))]
    if chart == "schwarzschild_gauge":
        return [Mark(m, along(0.0, math.log(math.cosh(lo)), math.log(math.cosh(hi))))]
    if chart == "dilaton":
        return [Mark(m, along(0.0, math.cosh(lo) ** 2, math.cosh(hi) ** 2))]
    if chart == "conformal":
        return [Mark(m, along(0.0, -BIG, math.log(math.sinh(hi))))]
    if chart == "kruskal":
        return [Mark(m, [[(math.sinh(hi), -math.sinh(hi)), (-math.sinh(hi), math.sinh(hi))]])]
    sign = 1 if chart == "eddington_finkelstein_ingoing" else -1
    x = near(math.log(math.cosh(lo)), math.log(math.cosh(hi)))
    return [Mark(m, [np.column_stack([sign * np.log(np.expm1(2 * x)) / 2, x])])]


def _hiscock(view):
    """Hiscock's evaporating hole: each moment is the slice v - r = T over the pieces read in the
    ingoing chart and the slice u + r = null_rays.hiscock_outer_time(T) over those read in the
    outgoing chart, whose part beyond the last ray, u > 8, is the flat space after the hole in its
    double null chart, u = T - r and v = T + r."""
    import null_rays as nr

    def lines(m):
        try:
            lo, hi = m.reach("outgoing" if view == "after" else view, "r")
        except ValueError:
            return []
        if view == "ingoing":
            return [[(m.time + r, r) for r in np.linspace(lo, hi, 200)]]
        outer = nr.hiscock_outer_time(m.time)
        if view == "outgoing":
            return [[(outer - r, r) for r in np.linspace(lo, hi, 200)]]
        hi = min(hi, (m.time - nr.HISCOCK_V0) / 2)
        return [[(m.time - r, m.time + r) for r in np.linspace(lo, hi, 50)]] if hi > lo else []
    return [Mark(m, found) for m in moments("hiscock") for found in [lines(m)] if found]


def _sultana_dyer_t(m):
    """A moment of Sultana and Dyer's conformal time on the plane of Schwarzschild's t and r, r_s = 1:
    ct = eta - ln(r - 1), outside r_s, as far as the embedding reaches."""
    r = near(1.0, m.reach("kerr_schild", "r")[1])
    return [np.column_stack([m.time - np.log(r - 1), r])]


def _boson_star_isotropic(m):
    """The boson star's moment t = 0 in the isotropic radius: the embedding's reach along the areal
    radius r carried to R, where r = psi^2 R, by the solver that drew it."""
    import boson_star
    lo, hi = m.reach("areal", "r")
    return along(0.0, *(float(x) for x in boson_star.star().isotropic().rho_of([lo, hi])))
def _bm_radial(which, r):
    """The isotropic radius or the tortoise coordinate of the sphere of areal radius r on Bartnik
    and McKinnon's soliton with one zero, in units of ell."""
    import bartnik_mckinnon
    return float(getattr(bartnik_mckinnon.soliton(1), which)(r)[0])


def one(metric_id, lines_of, label=None, view_id=None):
    """Each moment of a spacetime as the lines lines_of(moment) returns."""
    return [Mark(m, lines_of(m), label=label) for m in moments(metric_id, view_id)]


# Kerr-de Sitter as its diagrams draw it, (r_s, a, Lambda): Lambda > 0 in units of r_s and
# Lambda < 0 in units of l = sqrt(-3/Lambda).
KDS = {"de_sitter": (1.0, 0.45, 0.2), "anti_de_sitter": (2.0, 0.5, -3.0)}


def primitive(numerator, denominator):
    """The primitive of the rational function N/D that vanishes at r = 0, for a D with simple
    roots r_i, none of them zero, and a N of lower degree: Re sum_i A_i ln(1 - r/r_i) with
    A_i = N(r_i)/D'(r_i), the principal logarithm, whose real part is ln|1 - r/r_i| at a real
    root. Both polynomials are numpy's, highest power first."""
    roots = np.roots(denominator)
    slope = np.polyder(denominator)
    residues = [np.polyval(numerator, z) / np.polyval(slope, z) for z in roots]

    def value(r):
        r = np.asarray(r, dtype=float)
        with np.errstate(divide="ignore"):
            return sum(A * np.log((1 - r / z).astype(complex)) for A, z in zip(residues, roots)).real
    return value


def kds_delta(sign):
    """Delta_r = (r^2 + a^2)(1 - Lambda r^2/3) - r_s r as a polynomial in r."""
    rs, a, L = KDS[sign]
    return [-L / 3, 0.0, 1 - L * a * a / 3, -rs, a * a]


def kds_rstar(sign):
    """Kerr-de Sitter's r_*, dr_*/dr = (r^2 + a^2)/Delta_r with r_* = 0 at r = 0, as the Kerr
    charts fix it: v = ct + r_* and u = ct - r_*."""
    return primitive([1.0, 0.0, KDS[sign][1] ** 2], kds_delta(sign))


# Kerr-Taub-NUT as its diagrams draw it, (m, a, l), in units of m.
KTN = (1.0, 1.0, 1.25)


def ktn_rstar():
    """Kerr-Taub-NUT's r_* on the half of the axis that is regular, dr_*/dr = (r^2 + (a + l)^2)/Delta
    with Delta = r^2 - 2mr + a^2 - l^2 and r_* = 0 at r = 0, as the Kerr charts fix it:
    v = c t_N + r_* and u = c t_N - r_*. The integrand is 1 + (2mr + 2l(a + l))/Delta."""
    m, a, l = KTN
    rest = primitive([2 * m, 2 * l * (a + l)], [1.0, -2 * m, a * a - l * l])
    return lambda r: np.asarray(r, dtype=float) + rest(r)


def kds_schild(sign):
    """c(tau - t) of Kerr-de Sitter's Kerr-Schild chart, d/dr of it r_s r/((1 - Lambda r^2/3) Delta_r),
    zero at r = 0."""
    rs, a, L = KDS[sign]
    return primitive([rs, 0.0], np.polymul([-L / 3, 0.0, 1.0], kds_delta(sign)))


def _kds(sign, shift=None, scale=1):
    """The moment t = 0 of Kerr-de Sitter on a plane of a time and r: the line t = 0 as far as
    the embedding reaches in Carter's chart, or in a chart whose time is ct + scale shift(r) the
    curve scale shift(r), crowding toward each horizon it runs off at."""
    m, = moments("kerr_de_sitter", sign)
    lo, hi = m.reach("boyer_lindquist", "r")
    if shift is None:
        return [Mark(m, along(0.0, lo, hi))]
    s = np.linspace(-30, 30, N)
    r = lo + (hi - lo) / (1 + np.exp(-s)) if sign == "de_sitter" else near(lo, hi)
    return [Mark(m, [np.column_stack([scale * shift(sign)(r), r])])]


def _kds_above(sign):
    """The equator of the moment from above, the time left out: the whole plane between the
    radii the embedding reaches, every angle, since each chart's angle differs from Carter's by
    a function of r or of t alone."""
    m, = moments("kerr_de_sitter", sign)
    lo, hi = m.reach("boyer_lindquist", "r")
    phi = np.linspace(0, 2 * np.pi, 721)
    return [Mark(m, fills=[[np.column_stack([phi, np.full_like(phi, hi)]),
                            np.column_stack([phi, np.full_like(phi, lo)])]])]


def flat(spec):
    """The moments the flat view `spec` of null_rays.py shows, in its chart's (x^0, r)."""
    key = (spec.metric, spec.system, spec.view)
    if key in HIDDEN or spec.metric not in FLAT_METRICS:
        return []
    return FLAT[key]()


def _godel_cylinder(r):
    """A cylinder of Godel's or van Stockum's at radius r inside the embedding's reach: the
    whole line t = 0, every phi."""
    def lines(m):
        if not r <= m.reach(None, "r")[1]:
            raise ValueError("the cylinder lies beyond the embedding")
        return [[(0.0, -math.pi), (0.0, math.pi)]]
    return lines


def _bonnor_circle(rho):
    """A circle of the plane z = 0 of Bonnor's dust cloud through time, at a rho within the
    embedding's reach: the whole line t = 0, every phi."""
    def lines(m):
        lo, hi = m.reach("cylindrical", "\\rho")
        if not lo <= rho <= hi:
            raise ValueError("the circle lies outside the embedding")
        return [[(0.0, -math.pi), (0.0, math.pi)]]
    return lines


def _spinning_cylinder(R, line=None):
    """A cylinder about the spinning string at the circumference radius R = sqrt(b^2 r^2 - a^2),
    which must lie within the embedding's reach: the whole line t = 0, every phi, or the line
    given."""
    def lines(m):
        lo, hi = m.reach("circumference_radius", "R")
        if not lo <= R <= hi:
            raise ValueError("the cylinder lies outside the embedding")
        return line or [[(0.0, -math.pi), (0.0, math.pi)]]
    return lines


def _lewis_cylinder(r):
    """A cylinder about Lewis's axis at the radius r, which must lie within the embedding's reach:
    the whole line t = 0, every phi."""
    def lines(m):
        lo, hi = m.reach("canonical", "r")
        if not lo <= r <= hi:
            raise ValueError("the cylinder lies outside the embedding")
        return [[(0.0, -math.pi), (0.0, math.pi)]]
    return lines


def _os_interior():
    return one("oppenheimer_snyder", lambda m: [[(m.time / OS_AM, lo) for lo in m.reach("interior_comoving", "\\chi")]])


def _krasnikov():
    """The tube's moment ct = 5 along the axis, x = X + 2, the frame's X running from the
    path's start at x = 0: the whole plane of the height, -1 <= x <= 5."""
    m = moments("krasnikov", label="$ct = 5$")[0]
    path = next(c for c in m.surface["curves"] if c["class"] == "path")
    offset = -path["points"][0][0]
    u = m.grid()["u"]
    return [Mark(m, [[(KRASNIKOV_T, u[0] + offset), (KRASNIKOV_T, u[-1] + offset)]])]


KRASNIKOV_T = 5.0               # ct of the height, in rho_0, as the embedding's settings state


def _rn(view_id):
    return one("rn_metric", lambda m: along(0.0, *m.reach("spherical", "r")), view_id=view_id)


def _br(chart):
    """Bertotti-Robinson's moment t = 0: the equator along r from b/e to eb, and the sphere at
    the event r = b. The Poincare chart has the same t and x = b^2/r, which carries the one
    line element onto the other, so the same stretch is x from b/e to eb."""
    equator = moments("bertotti_robinson", "equator")[0]
    sphere = moments("bertotti_robinson", "sphere", label="$t = 0$, $r = b$")[0]
    lo, hi = equator.reach("static", "r")
    if chart == "poincare":
        lo, hi = 1 / hi, 1 / lo
        sphere.label = "$t = 0$, $x = b$"
    return [Mark(equator, along(0.0, lo, hi)), Mark(sphere, points=[(0.0, 1.0)])]


def nhek_radius(y):
    """Bardeen and Horowitz's Poincare radius, in r_0, of the event at y on the moment tau = 0 of
    their global chart, where t = 0 too: r = sqrt(1 + y^2) + y."""
    return math.sqrt(1 + y * y) + y


def _nhek(chart):
    """The throat of extreme Kerr at the moment tau = 0 of the global chart: the equator along y
    as far as the cylinder reaches, and the horizon's sphere at the event y = 0. Bardeen and
    Horowitz's map puts tau = 0 on t = 0 of their Poincare chart with r = r_0 (sqrt(1 + y^2) + y)
    and the same phi, and the inverse radius is x = r_0^2/r, so the same stretch is a line of
    t = 0 on each of those planes and the sphere stands at r = r_0 and at x = r_0."""
    throat = moments("near_horizon_extreme_kerr", "throat", label="$\\tau = 0$")[0]
    sphere = moments("near_horizon_extreme_kerr", "horizon", label="$\\tau = 0$, $y = 0$")[0]
    lo, hi = throat.reach("global", "y")
    at = 0.0
    if chart != "global":
        lo, hi, at = nhek_radius(lo), nhek_radius(hi), 1.0
        throat.label = "$t = 0$"
        sphere.label = "$t = 0$, $r = r_0$"
    if chart == "inverse_radius":
        lo, hi = 1 / hi, 1 / lo
        sphere.label = "$t = 0$, $x = r_0$"
    return [Mark(throat, along(0.0, lo, hi)), Mark(sphere, points=[(0.0, at)])]


def _ds_flat():
    """de Sitter's static moment t = 0 in the flat slicing, H = 1: with
    X_0 = sinh t_f + rho^2 e^(t_f)/2 of the embedding space it is X_0 = 0, which is
    t_f = -ln(1 + rho^2)/2, and rho = |x| on the plane y = z = 0. The whole of the observer's
    hemisphere lies on it, its r = rho e^(t_f) reaching 1 only as rho -> infinity."""
    m = moments("de_sitter", label="static $t = 0$")[0]
    x = np.linspace(-BIG ** 0.25, BIG ** 0.25, 4 * N + 1)
    return [Mark(m, [np.column_stack([-0.5 * np.log1p(x * x), x])])]


def _es_areal(m):
    """The Einstein static universe's moment in its areal chart, R = 1: r = sin chi over the
    near hemisphere the embedding reaches, chi from 0 to pi/2."""
    lo, hi = m.reach("hyperspherical", "\\chi")
    return along(0.0, math.sin(lo), math.sin(min(hi, math.pi / 2)))


def _milne_inertial(m):
    """A moment ct of the Milne universe in its inertial chart: T = t cosh chi and R = ct sinh chi,
    the hyperbola cT = sqrt(c^2t^2 + R^2) through the centre, out to the chi the embedding reaches."""
    hi = m.time * math.sinh(m.reach("comoving_hyperbolic", "\\chi")[1])
    R = np.linspace(-hi, hi, N)
    return [np.column_stack([np.sqrt(m.time ** 2 + R * R), R])]


def _cdl_static(m):
    """A moment c tau of the open universe inside Coleman and De Luccia's anti-de Sitter bubble in
    the static chart inside the wall, l = 1: r = sin(c tau) sinh(chi) and
    ct = atan2(sqrt(sin^2(c tau) + r^2), cos(c tau)), out to the chi the embedding reaches."""
    hi = math.sin(m.time) * math.sinh(m.reach("open", "\\chi")[1])
    r = np.linspace(0.0, hi, N)
    return [np.column_stack([np.arctan2(np.sqrt(math.sin(m.time) ** 2 + r * r), math.cos(m.time)), r])]


def _wall_inertial(m):
    """A moment kct of the domain wall's global chart in the inertial chart of the side z < 0, k = 1:
    cT = (1 - |z|) sinh(kct) and R = (1 - |z|) cosh(kct), the line cT = |R| tanh(kct) through the
    centre out to the wall at R = cosh(kct) on either side."""
    hi = math.cosh(m.time)
    R = np.linspace(-hi, hi, N)
    return [np.column_stack([np.abs(R) * math.tanh(m.time), R])]


def ks_vacuum_T(m):
    """The areal radius T, in units of r_s, of a vacuum moment of the Kantowski-Sachs embedding,
    which is the radius of its cylinder; its time is the proper time since the horizon."""
    return m.surface["pieces"][0]["points"][0][1]


def _kantowski_sachs(chart):
    """The dust universe's moments eta_k, each every r: level in the dust chart, and in the comoving
    chart at ct = pi/2 + eta + sin(eta) cos(eta), counted from the first singularity at b_0 = 1. The
    vacuum moments, each every r at its T, lie on the plane inside Schwarzschild's horizon alone."""
    if chart == "schwarzschild_interior":
        return one("kantowski_sachs", lambda m: across(ks_vacuum_T(m), 0.0, BIG), view_id="vacuum")
    time = {"dust": lambda e: e, "comoving": lambda e: math.pi / 2 + e + math.sin(e) * math.cos(e)}[chart]
    return one("kantowski_sachs", lambda m: across(time(m.time), 0.0, BIG), view_id="dust")


def _ads_poincare():
    """Anti-de Sitter's static moment t = 0 is the Poincare moment t = 0; on the plane y = 0,
    z = L its static radius is r^2 = x^2 + x^4/4L^2, so the embedding's reach r <= 4L is
    |x| <= sqrt(2(sqrt(17) - 1)) L."""
    m = moments("anti_de_sitter")[0]
    hi = m.reach("static_global", "r")[1]
    x = math.sqrt(2 * (math.sqrt(1 + hi * hi) - 1))
    return [Mark(m, across(0.0, 0.0, x))]


def btz_rstar(r):
    """The BTZ hole's r_* = (1/2) ln|(r - 1)/(r + 1)| at M = 1, l = 1, vanishing as r -> infinity,
    which fixes the Eddington-Finkelstein charts' v = ct + r_* and u = ct - r_*."""
    r = np.asarray(r, dtype=float)
    return 0.5 * np.log(np.abs((r - 1) / (r + 1)))


def _btz(sign=0):
    """The moment t = 0 of the BTZ hole without rotation, from the throat r_+ = 1 out: along r
    in its stationary chart (sign 0), and in its ingoing (1) or outgoing (-1) chart as v = r_* or
    u = -r_*, crowding toward the horizon, where the curve runs off. Each chart covers one
    exterior, and the moment's other exterior lies over the same r."""
    m, = moments("btz")
    lo, hi = m.reach("stationary", "r")
    if not sign:
        return [Mark(m, along(0.0, lo, hi))]
    r = near(lo, hi)
    return [Mark(m, [np.column_stack([sign * btz_rstar(r), r])])]


def sads_rstar(r):
    """Schwarzschild-anti-de Sitter's tortoise coordinate at r_s = 2 and L = 1, where 1/f =
    r/((r - 1)(r^2 + r + 2)), as the Eddington-Finkelstein charts fix it, vanishing at r = 0:
    r_* = (1/4) ln|1 - r| - (1/8) ln((r^2 + r + 2)/2) + (5/(4 sqrt 7))(arctan((2r + 1)/sqrt 7) -
    arctan(1/sqrt 7))."""
    r = np.asarray(r, dtype=float)
    w = math.sqrt(7.0)
    return (0.25 * np.log(np.abs(1 - r)) - 0.125 * np.log((r * r + r + 2) / 2)
            + 5 / (4 * w) * (np.arctan((2 * r + 1) / w) - math.atan(1 / w)))


def _sads(sign=0):
    """The moment t = 0 of the Schwarzschild-anti-de Sitter hole, from the throat r_h = 1 out: along
    r in its static chart (sign 0), and in its ingoing (1) or outgoing (-1) chart as v = r_* or
    u = -r_*, crowding toward the horizon, where the curve runs off. Each chart covers one
    exterior, and the moment's other exterior lies over the same r."""
    m, = moments("schwarzschild_ads")
    lo, hi = m.reach("static", "r")
    if not sign:
        return [Mark(m, along(0.0, lo, hi))]
    r = near(lo, hi)
    return [Mark(m, [np.column_stack([sign * sads_rstar(r), r])])]


RNDS_ROOTS = (2.0, 2 / 3, (2 * math.sqrt(7) - 4) / 3, -(2 * math.sqrt(7) + 4) / 3)
RNDS_H = 3 / 8                  # H r_s/c of the lukewarm hole every diagram draws
RNDS_LABEL = "$\\tau = 1/H$"      # the name of that moment on every drawing
RNDS_TAU = 1 / RNDS_H           # the moment of the cosmological chart that is embedded, H tau = 1


def rnds_rstar(r):
    """The lukewarm hole's tortoise coordinate at r_s = 1, r_q = 1/2 and Lambda = 27/64, as the
    Eddington-Finkelstein charts fix it, r_* = sum_i ln|1 - r/r_i|/f'(r_i) over the four roots of
    9r^4 - 64r^2 + 64r - 16, with 1/f'(r_i) = -64 r_i^2/(9 prod_j (r_i - r_j)), which vanishes at r = 0."""
    r = np.asarray(r, dtype=float)
    return sum(np.log(np.abs(1 - r / a)) * (-64 * a * a / (9 * math.prod(a - b for b in RNDS_ROOTS if b != a)))
               for a in RNDS_ROOTS)


def rnds_static_t(tau, r):
    """The static time cT of the event of the cosmological chart at tau with the areal radius
    r = H tau rho + 1/2: cT = ln|H tau|/H + F(r) up to a constant, with F' = H r^2/((r - 1/2) f), a
    sum of logarithms c_a ln|r - a| over a = 1/2 and the four roots of f,
    c_a = -(8/3) a^4/prod_b (a - b). c at 1/2 is -1/H, so cT is continuous through tau = 0, where
    r = 1/2, and its one constant is taken so that cT vanishes at H tau = 1 on the static radius
    r = 1.298, where f is greatest, the root of 9r^4 - 32r + 16."""
    r = np.asarray(r, dtype=float)
    poles = (0.5,) + RNDS_ROOTS

    def F(x):
        return sum(np.log(np.abs(x - a)) * (-(8 / 3) * a ** 4 / math.prod(a - b for b in poles if b != a)) for a in poles)
    return np.log(np.abs(RNDS_H * tau)) / RNDS_H + F(r) - F(1.2979848366419)


# The static time inside r_- is fixed only up to a constant of its own, and the moment embedded
# there is the one through H tau = -1/2 on r = 0.35, which lies on the drawn cosmological plane.
RNDS_INSIDE_T = float(rnds_static_t(-0.5 / RNDS_H, 0.35))


# Kastor and Traschen's holes as every diagram draws them, falling together, H < 0, in units of
# the mass parameter m of one hole. One hole alone at H = -3c/(16m) is the lukewarm hole of the
# Reissner-Nordstrom-de Sitter diagrams, r_s = 2m and 4m|H|/c = 3/4, so its horizons are twice
# RNDS_ROOTS. Two holes of mass parameter m each, on the axis at z = +-2m, at H = -3c/(32m), half
# that rate, merge into the same lukewarm hole of mass parameter 2m, whose horizons are at the
# areal radii 8m/3 and 8m, the values H tau r + 2m takes far from the pair.
KT_H_ONE = -3 / 16
KT_H_TWO = -3 / 32
KT_ONE_ROOTS = tuple(2 * a for a in RNDS_ROOTS)


def kt_potential(plane, x):
    """V of the two holes on the axis through them and on the plane midway between them."""
    x = np.asarray(x, dtype=float)
    if plane == "axis":
        return 1 / np.abs(x - 2) + 1 / np.abs(x + 2)
    return 2 / np.sqrt(x * x + 4)


@functools.lru_cache(maxsize=None)
def kt_last_ray(plane):
    """c tau along the last ray of a plane of the two holes to reach infinity, as a function of
    the distance x > 2 along the axis or x >= 0 across the midplane: the event horizon there.

    An outgoing ray has d(c tau)/dx = U^2, U = H tau + V. With R = H tau x and y = ln x that is
    dR/dy = R + H(R + xV)^2, and far from the pair, where xV -> 2, the right hand side vanishes at
    R = 2/3 and R = 6, the horizons of the merged hole less its mass parameter 2: a ray above
    R = 2/3 runs on to R = 6, the cosmological horizon, a ray below it falls to U = 0, and the
    one ray that tends to R = 2/3 divides them, as Brill, Horowitz, Kastor and Traschen argue
    from their (4.6). The slope of the right hand side at R = 2/3 is +1/2, so that ray is found
    by integrating inward from x = 1e7, where every ray closes on it."""
    from scipy.integrate import solve_ivp
    H = KT_H_TWO
    lo = 2 + 1e-4 if plane == "axis" else 1e-9

    def slope(y, R):
        x = math.exp(y)
        return [R[0] + H * (R[0] + x * float(kt_potential(plane, x))) ** 2]
    ray = solve_ivp(slope, [math.log(1e7), math.log(lo)], [2 / 3], rtol=1e-12, atol=1e-14, dense_output=True)

    def tau(x):
        x = np.maximum(np.asarray(x, dtype=float), lo)
        return ray.sol(np.log(x))[0] / (H * x)
    return tau


def kt_merger():
    """c tau at which the event horizon first reaches the midplane, at the midpoint of the pair,
    where the last ray of the midplane leaves the axis: -6.1995 m."""
    return float(kt_last_ray("midplane")(0.0))


def kt_comoving_t(tau):
    """The comoving chart's ct at the cosmological time tau < 0 of two holes: H tau = e^{Ht}."""
    return np.log(KT_H_TWO * np.asarray(tau, dtype=float)) / KT_H_TWO


def kt_static_t(tau, r):
    """The static time cT of the event of one hole's isotropic chart at tau and r, m = 1 and
    H = -3/16, with the areal radius R = H tau r + 1: cT = ln|H tau|/H + F(R), F' = H R^2/((R - 1) f)
    and f = (1 - 1/R)^2 - 9R^2/256, Brill, Horowitz, Kastor and Traschen's (2.4). F is a sum of
    logarithms c_a ln|R - a| over a = 1 and the four roots of f, c_a = (16/3) a^4/prod_b (a - b)."""
    R = KT_H_ONE * np.asarray(tau, dtype=float) * np.asarray(r, dtype=float) + 1
    poles = (1.0,) + KT_ONE_ROOTS
    F = sum(np.log(np.abs(R - a)) * ((16 / 3) * a ** 4 / math.prod(a - b for b in poles if b != a)) for a in poles)
    return np.log(np.abs(KT_H_ONE * tau)) / KT_H_ONE + F


def kt_rstar(R):
    """One hole's tortoise coordinate at m = 1 and H = -3/16, sum_i ln|1 - R/R_i|/f'(R_i) over the
    four roots, with 1/f'(R_i) = -256 R_i^2/(9 prod_j (R_i - R_j)), which vanishes at R = 0."""
    R = np.asarray(R, dtype=float)
    return sum(np.log(np.abs(1 - R / a)) * (-256 * a * a / (9 * math.prod(a - b for b in KT_ONE_ROOTS if b != a)))
               for a in KT_ONE_ROOTS)


def _rnds(chart):
    """The three moments of the lukewarm hole on each of its planes: the static t = 0 between r_+
    and r_c and inside r_-, and the moment H tau = 1 of the cosmological chart, on which the
    areal radius is rho + 1/2. In an Eddington-Finkelstein chart the static moments are v = r_*
    and u = -r_*, and the cosmological one is v = cT + r_* or u = cT - r_* with cT of
    rnds_static_t; on the cosmological plane the static t = 0 is |H tau| = exp(-H F(r)),
    rho = (r - 1/2)/(H tau), and inside r_-, where tau < 0, it is the moment cT = RNDS_INSIDE_T of that time. The cosmological moment runs through the
    white hole between r_s/2 and r_+, which the outgoing chart alone covers among the others: the
    static plane reads that range of r as the black hole and the ingoing chart covers the black
    hole, so they draw the moment from r_+ out, and the ingoing chart no further than r_c."""
    between, = moments("reissner_nordstrom_de_sitter", "between")
    inside, = moments("reissner_nordstrom_de_sitter", "inside")
    cosmic, = moments("reissner_nordstrom_de_sitter", "cosmological", label=RNDS_LABEL)
    rc, rp, rm, _ = RNDS_ROOTS
    b_lo, b_hi = between.reach("static", "r")
    i_lo, i_hi = inside.reach("static", "r")
    c_lo, c_hi = (x + 0.5 for x in cosmic.reach("cosmological", "\\rho"))

    def span(lo, hi, open_lo=True, open_hi=True):
        """Values from lo to hi, crowding toward each end the curve runs off at."""
        mid = 0.5 * (lo + hi)
        left = near(lo, mid) if open_lo else np.linspace(lo, mid, N)
        right = near(hi, mid)[::-1] if open_hi else np.linspace(mid, hi, N)
        return np.concatenate([left, right[1:]])

    r_b, r_i = span(b_lo, b_hi), span(i_lo, i_hi, open_lo=False)
    r_static = span(rp, rc)
    r_beyond = span(rc, c_hi, open_hi=False)
    r_white = span(c_lo, rp, open_lo=False)
    if chart == "static":
        return [Mark(between, along(0.0, b_lo, b_hi)), Mark(inside, along(0.0, i_lo, i_hi)),
                Mark(cosmic, [np.column_stack([rnds_static_t(RNDS_TAU, r), r]) for r in (r_static, r_beyond)])]
    if chart == "cosmological":
        def static_moment(r, sign, at):
            tau = sign * np.exp(RNDS_H * (at - rnds_static_t(1 / RNDS_H, r))) / RNDS_H
            return np.column_stack([tau, (r - 0.5) / (RNDS_H * tau)])
        return [Mark(between, [static_moment(r_b, 1, 0.0)]), Mark(inside, [static_moment(r_i, -1, RNDS_INSIDE_T)]),
                Mark(cosmic, along(RNDS_TAU, c_lo - 0.5, c_hi - 0.5))]
    sign = 1 if chart == "ingoing" else -1
    # In the outgoing chart the moment crosses r_+ and r_c as one curve.
    pieces = (r_static,) if sign == 1 else (np.concatenate([r_white, r_static, r_beyond]),)
    return [Mark(between, [np.column_stack([sign * rnds_rstar(r_b), r_b])]),
            Mark(inside, [np.column_stack([sign * rnds_rstar(r_i), r_i])]),
            Mark(cosmic, [np.column_stack([rnds_static_t(RNDS_TAU, r) + sign * rnds_rstar(r), r]) for r in pieces])]


BARDEEN_G = 1 / 3                         # g in r_s, as every diagram of Bardeen's black hole takes it
BARDEEN_HORIZONS = (0.3009628011521132, 0.7754191789541744)     # r_- and r_+, the positive zeros of f


def bardeen_f(r):
    """Bardeen's f at r_s = 1 and g = 1/3."""
    return 1 - r * r / (r * r + BARDEEN_G ** 2) ** 1.5


@functools.lru_cache(maxsize=None)
def _bardeen_smooth():
    """What is left of 1/f once its two poles are taken out, and its integral from 0 to each
    multiple of a tenth of r_s up to 8 r_s and on from there to each half step of ln r. Within a
    twentieth of r_s of a horizon it is worked in forty digits, since 1/f and its pole cancel
    there to the last digits of a float."""
    import mpmath
    g2 = mpmath.mpf(1) / 9
    with mpmath.workdps(40):
        def f(x):
            return 1 - x * x / (x * x + g2) ** mpmath.mpf("1.5")
        roots = [mpmath.findroot(f, ri) for ri in BARDEEN_HORIZONS]
        exact = [(ri, (ri * ri + g2) ** mpmath.mpf("2.5") / (ri * (ri * ri - 2 * g2))) for ri in roots]
    poles = [(float(ri), float(a)) for ri, a in exact]

    def smooth(x):
        if min(abs(x - ri) for ri, _ in poles) > 0.05:
            return 1 / bardeen_f(x) - sum(a / (x - ri) for ri, a in poles)
        with mpmath.workdps(40):
            x = mpmath.mpf(x)
            return float(1 / f(x) - sum(a / (x - ri) for ri, a in exact))

    nodes, weights = np.polynomial.legendre.leggauss(12)

    def panel(a, b, of=smooth):
        half, mid = 0.5 * (b - a), 0.5 * (a + b)
        return half * sum(w * of(mid + half * n) for n, w in zip(nodes, weights))

    def in_log(u):
        return smooth(math.exp(u)) * math.exp(u)
    edges = [0.1 * k for k in range(81)]
    sums = np.concatenate([[0.0], np.cumsum([panel(a, b) for a, b in zip(edges, edges[1:])])])
    logs = [math.log(8.0) + 0.5 * k for k in range(61)]
    tails = sums[-1] + np.concatenate([[0.0], np.cumsum([panel(a, b, in_log) for a, b in zip(logs, logs[1:])])])

    def integral(x):
        """The integral of the smooth part from 0 to x."""
        if x <= 8.0:
            k = min(int(x / 0.1), 80)
            return float(sums[k]) + (panel(0.1 * k, x) if x > 0.1 * k else 0.0)
        u = math.log(x)
        k = min(int((u - logs[0]) / 0.5), 60)
        return float(tails[k]) + (panel(logs[k], u, in_log) if u > logs[k] else 0.0)
    return integral, poles


def bardeen_rstar(r):
    """Bardeen's tortoise coordinate at r_s = 1 and g = 1/3, dr_*/dr = 1/f, as the
    Eddington-Finkelstein charts fix it, vanishing at r = 0. 1/f has a simple pole at each
    horizon, of residue 1/f'(r_i), so r_* is the sum over both of ln|1 - r/r_i|/f'(r_i) and the
    integral from 0 of what is left of 1/f once those poles are taken out, which is smooth and is
    summed by Gauss and Legendre's rule, on panels a tenth of r_s wide out to 8 r_s and half a
    unit of ln r wide beyond."""
    integral, poles = _bardeen_smooth()

    def one(x):
        if math.isinf(x):
            return math.inf
        with np.errstate(divide="ignore"):
            return integral(x) + sum(a * float(np.log(abs(1 - x / ri))) for ri, a in poles)
    return np.array([one(float(x)) for x in np.atleast_1d(np.asarray(r, dtype=float))]).reshape(np.shape(r))


def _bardeen(view, sign=0):
    """The moment t = 0 of Bardeen's black hole, outside r_+ or inside r_-: along r in its static
    chart (sign 0), and in its ingoing (1) or outgoing (-1) chart as v = r_* or u = -r_*, which
    runs off toward the horizon the view ends on. Each chart covers one side of the moment, and
    the other side lies over the same r."""
    m, = moments("bardeen", view)
    lo, hi = m.reach("static", "r")
    if not sign:
        return [Mark(m, along(0.0, lo, hi))]
    r = near(lo, hi) if view == "outside" else np.concatenate([[0.0], (hi - (hi - lo) * np.geomspace(1.0, 1e-9, N))[1:]])
    return [Mark(m, [np.column_stack([sign * bardeen_rstar(r), r])])]


def _bardeen_both(sign=0):
    return _bardeen("outside", sign) + _bardeen("inside", sign)


def hayward_rstar(r):
    """Hayward's tortoise coordinate at m = 1 and ell = 12/(7 sqrt 7), where 1/F = 1 + 2r^2/((r - 6/7)
    (r - 12/7)(r + 4/7)), as the Eddington-Finkelstein charts fix it, vanishing at r = 0:
    r_* = r - (6/5) ln|1 - 7r/6| + 3 ln|1 - 7r/12| + (1/5) ln(1 + 7r/4)."""
    r = np.asarray(r, dtype=float)
    return r - 1.2 * np.log(np.abs(1 - 7 * r / 6)) + 3 * np.log(np.abs(1 - 7 * r / 12)) + 0.2 * np.log(1 + 7 * r / 4)


def _hayward(sign=0):
    """The moment t = 0 of Hayward's black hole, outside r_+ = 12/7 and inside r_- = 6/7: along r in
    its static chart (sign 0), and in its ingoing (1) or outgoing (-1) chart as v = r_* or
    u = -r_*, each part crowding toward its horizon, where the curve runs off. The static t of the
    region inside r_- is the one the same r_* gives there."""
    out = []
    for view_id in ("outside", "inside"):
        m, = moments("hayward", view_id)
        lo, hi = m.reach("static", "r")
        if not sign:
            out.append(Mark(m, along(0.0, lo, hi)))
            continue
        r = near(lo, hi) if view_id == "outside" else (lo + hi) - near(lo, hi)[::-1]
        out.append(Mark(m, [np.column_stack([sign * hayward_rstar(r), r])]))
    return out


# The charged black hole in anti-de Sitter space as every diagram draws it, L = 1, r_s = 27/8 and
# r_q^2 = 11/8: r^2 f = (r - 1)(r - 1/2)(r^2 + 3r/2 + 11/4), with the roots below and the residues
# 1/f'(r_i) of 1/f at them, which sum to zero.
RNADS_ROOTS = (1.0, 0.5, complex(-0.75, math.sqrt(35) / 4), complex(-0.75, -math.sqrt(35) / 4))
RNADS_RESIDUES = (8 / 21, -2 / 15, complex(-13 / 105, -math.sqrt(35) / 35), complex(-13 / 105, math.sqrt(35) / 35))


def rnads_rstar(r):
    """The tortoise coordinate of the charged black hole in anti-de Sitter space at L = 1, r_s = 27/8
    and r_q^2 = 11/8, as the Eddington-Finkelstein charts fix it, vanishing at r = 0:
    r_* = Re sum_i ln(1 - r/r_i)/f'(r_i) over the four roots of r^2 f, which is
    (8/21) ln|1 - r| - (2/15) ln|1 - 2r| plus the complex pair's part, and tends to 0.4052 as
    r -> infinity."""
    r = np.asarray(r, dtype=float)
    with np.errstate(divide="ignore"):
        return sum(A * np.log((1 - r / z).astype(complex)) for A, z in zip(RNADS_RESIDUES, RNADS_ROOTS)).real


def _rnads(sign=0):
    """The moment t = 0 of the charged black hole in anti-de Sitter space, outside r_+ = 1 and inside
    r_- = 1/2: along r in its static chart (sign 0), and in its ingoing (1) or outgoing (-1) chart
    as v = r_* or u = -r_*, each part crowding toward its horizon, where the curve runs off. The
    static t of the region inside r_- is the one the same r_* gives there."""
    out = []
    for view_id in ("outside", "inside"):
        m, = moments("reissner_nordstrom_ads", view_id)
        lo, hi = m.reach("static", "r")
        if not sign:
            out.append(Mark(m, along(0.0, lo, hi)))
            continue
        r = near(lo, hi) if view_id == "outside" else (lo + hi) - near(lo, hi)[::-1]
        out.append(Mark(m, [np.column_stack([sign * rnads_rstar(r), r])]))
    return out


def tbh_rstar(r):
    """The flat topological black hole's tortoise coordinate at mu = 1 and L = 1, where 1/f =
    r/((r - 1)(r^2 + r + 1)), as the Eddington-Finkelstein charts fix it, vanishing at r = 0:
    r_* = (1/3) ln|1 - r| - (1/6) ln(r^2 + r + 1) + (arctan((2r + 1)/sqrt 3) - pi/6)/sqrt 3."""
    r = np.asarray(r, dtype=float)
    w = math.sqrt(3.0)
    return (np.log(np.abs(1 - r)) / 3 - np.log(r * r + r + 1) / 6
            + (np.arctan((2 * r + 1) / w) - math.pi / 6) / w)


def _tbh(sign=0, brane=False):
    """The black string's moment t = 0, the flat hole's, from the throat r_h = 1 out: along r in
    the static chart and in Lemos's (sign 0), in the ingoing (1) or outgoing (-1) chart as
    v = r_* or u = -r_*, crowding toward the horizon, where the curve runs off, and on the
    brane's plane along z = L^2/r, from the horizon z_h = 1 toward the boundary."""
    m, = moments("topological_black_hole", "string")
    lo, hi = m.reach("black_string", "r")
    if brane:
        return [Mark(m, along(0.0, 1 / hi, 1 / lo))]
    if not sign:
        return [Mark(m, along(0.0, lo, hi))]
    r = near(lo, hi)
    return [Mark(m, [np.column_stack([sign * tbh_rstar(r), r])])]


def _tbh_horizon():
    """The hyperbolic hole's horizon at one moment, mu = 0: the bifurcation surface, the point
    t = 0, r = r_h = 1 of the static planes. In an Eddington-Finkelstein chart it lies at
    v -> -infinity or u -> +infinity, off the chart."""
    m, = moments("topological_black_hole", "horizon", label="$t = 0$, $r = r_h$")
    return [Mark(m, points=[(0.0, 1.0)])]


def _siklos(depth=None):
    """Siklos's wave front u = 0, v = 0, read on the disc out to the circle a proper distance 2L
    from its centre. A plane of u and v at one place on the front meets it at the event (0, 0).
    Kaigorodov's planes with a Killing direction divided out and y = 0 meet it along the line of
    zero v, or t, or U, or V, over the diameter eta = 0 of the disc, where Siklos's
    x = L(2L - xi)/(2L + xi), and `depth` writes that x in the plane's own coordinate."""
    m, = moments("siklos", "front", label="$u = v = 0$")
    if depth is None:
        return [Mark(m, points=[(0.0, 0.0)])]
    _, top = m.reach("ozsvath_robinson_rozga", "\\xi")
    ends = sorted(depth(x) for x in ((2 - top) / (2 + top), (2 + top) / (2 - top)))
    return [Mark(m, along(0.0, *ends))]


def _kundt(view_id):
    """A wave front u = 0, v = 0 of the family with a cosmological constant, the hemisphere in
    de Sitter space or the half plane in anti-de Sitter space. The plane of u and v at one place
    on the front, drawn for the same member of the family, meets it at the event (0, 0)."""
    m, = moments("kundt_waves", view_id, label="$u = v = 0$")
    return [Mark(m, points=[(0.0, 0.0)])]


def soliton_radius(rho):
    """Horowitz and Myers's radius at the proper distance rho from the soliton's tip, at
    r_0 = L = 1: r = cosh^(2/3)(3 rho/2)."""
    return math.cosh(1.5 * rho) ** (2 / 3)


def _soliton(chart):
    """The soliton's moment t = 0, the surface of rho and phi the embedding reads in the polar
    chart, from the tip out: along rho on the polar planes, along r = r_0 cosh^(2/3)(3 rho/2L) on
    Horowitz and Myers's plane, and along z = L^2/r on the Poincare plane, where the tip is
    z_0 = 1 and the boundary lies toward z = 0."""
    m, = moments("ads_soliton")
    lo, hi = m.reach("polar", "\\rho")
    if chart == "polar":
        return [Mark(m, along(0.0, lo, hi))]
    if chart == "horowitz_myers":
        return [Mark(m, along(0.0, soliton_radius(lo), soliton_radius(hi)))]
    return [Mark(m, along(0.0, 1 / soliton_radius(hi), 1 / soliton_radius(lo)))]


def _schrodinger():
    """Schrodinger spacetime's surface t = 0, xi = 0, the plane of x and r, which is T = 0, V = 0
    of the global chart and the same hyperbolic plane for every dynamical exponent and in every
    dimension. A plane of the time and the null coordinate at one depth meets it at the event (0, 0)."""
    m, = moments("schrodinger_spacetime", "plane", label="$t = 0$, $\\xi = 0$")
    return [Mark(m, points=[(0.0, 0.0)])]


def _c_metric(y):
    """The C-metric's two moments on a plane of its axis: the equator's t = 0, which meets the
    axis along t = 0 over the same r as it reaches on the equator, and the black hole horizon,
    the bifurcation sphere at t = 0 and r = 2m, whose poles lie on the axis. On the Hong-Teo
    plane the same is tau = 0 over y = 1/(alpha r), alpha = 1/6."""
    equator = moments("c_metric", "equator")[0]
    horizon = moments("c_metric", "horizon", label="$t = 0$, $r = 2m$")[0]
    lo, hi = equator.reach("spherical", "r")
    if y:
        lo, hi = 6 / hi, 6 / lo
        horizon.label = "$\\tau = 0$, $y = 3$"
    return [Mark(equator, along(0.0, lo, hi)), Mark(horizon, points=[(0.0, 3.0 if y else 2.0)])]


def nariai_static_t(tau, r):
    """The Nariai universe's global moment ct = tau, Lambda = 1, in the static chart: with the static
    patch 0 < chi < pi, r = -cosh(tau) cos(chi) and sinh(t) = sinh(tau)/sqrt(1 - r^2), both from the
    embedding -Z0^2 + Z1^2 + Z2^2 = 1 of the de Sitter factor, Z0 = sinh(tau) = sqrt(1 - r^2) sinh(t)."""
    r = np.asarray(r, dtype=float)
    return np.arcsinh(math.sinh(tau) / np.sqrt(1 - r * r))


def _nariai(chart):
    """The Nariai universe's moments: each global moment of the movie, the whole circle of chi,
    which on the static plane is the curve nariai_static_t across the patch from horizon to
    horizon, and the sphere at the event t = 0, r = 0, which is chi = pi/2 of the global chart."""
    out = []
    for m in moments("nariai", "universe"):
        if chart == "global":
            out.append(Mark(m, [[(m.time, 0.0), (m.time, 2 * math.pi)]]))
        else:
            r = np.tanh(np.linspace(-12, 12, N))
            out.append(Mark(m, [np.column_stack([nariai_static_t(m.time, r), r])]))
    sphere = moments("nariai", "sphere", label="$t = 0$, $r = 0$")[0]
    if chart == "global":
        sphere.label = "$t = 0$, $\\chi = \\pi/2$"
        return out + [Mark(sphere, points=[(0.0, math.pi / 2)])]
    return out + [Mark(sphere, points=[(0.0, 0.0)])]


def _ori_hyperbola(m):
    """Ori's moment t = t_k < 0 on the Brinkmann plane through the central circle, uv = -2 t_k with
    u and v both negative, as (u, v), out to where either is 50 in size."""
    u = -np.geomspace(-2 * m.time / 50, 50, N)
    return [np.column_stack([u, -2 * m.time / u])]


def _mp_axis():
    """The midplane z = 0 between Majumdar and Papapetrou's two holes, at t = 0, meets the axis
    through them at one event, t = 0 and z = 0."""
    m = moments("majumdar_papapetrou", "two_holes", label="$t = 0$, $z = 0$")[0]
    return [Mark(m, points=[(0.0, 0.0)])]


def _nm_centre(label):
    """The plane of Neugebauer and Meinel's disc at t = 0 meets the axis at one event, the centre of
    the disc, the centre of the embedded surface."""
    m = moments("neugebauer_meinel", "plane", label=label)[0]
    return [Mark(m, points=[(0.0, 0.0)])]


def _double_kerr_axis():
    """The plane z = 0 midway between Kramer and Neugebauer's two holes, at t = 0, meets the axis
    through them at one event, t = 0 and z = 0, the tip of the embedded cone."""
    m = moments("double_kerr", "midplane", label="$t = 0$, $z = 0$")[0]
    return [Mark(m, points=[(0.0, 0.0)])]


def _morgan_morgan_centre():
    """The plane z = 0 of the first Morgan-Morgan disc, at t = 0, meets the axis at one event, the
    centre of the disc: z = 0 in Weyl's chart and xi = 0 in the oblate spheroidal one."""
    m = moments("morgan_morgan", "plane", label="$t = 0$, the centre of the disc")[0]
    return [Mark(m, points=[(0.0, 0.0)])]


def _bonnor_dipole_strut():
    """The equatorial plane of Bonnor's dipole at t = 0 meets the axis between the two black holes
    at one event, t = 0 and theta = pi/2, the tip of the embedded cone."""
    m = moments("bonnor_magnetic_dipole", "equator", label="$t = 0$, the equatorial plane")[0]
    return [Mark(m, points=[(0.0, math.pi / 2)])]


KT_LABEL = "$c\\tau = -8m/3$"     # the moment of one of Kastor and Traschen's holes that is embedded, H tau = 1/2
KT_TAU = -8 / 3


def _kt(chart, plane):
    """Kastor and Traschen's two moments on their flat views. One hole's moment c tau = -8m/3 lies
    on its own plane of tau and r. Each moment of the midplane between two holes is a line of
    constant tau across the midplane, at ct = ln(H tau)/H in the comoving chart, and meets the
    axis through the holes at one event, z = 0."""
    if chart == "isotropic":
        return one("kastor_traschen", lambda m: along(KT_TAU, *m.reach("isotropic", "r")), label=KT_LABEL,
                   view_id="one_hole")
    marks = []
    for m in moments("kastor_traschen", "two_holes"):
        t = float(kt_comoving_t(m.time)) if chart == "comoving" else m.time
        lo, hi = m.reach("cylindrical", "\\rho")
        if plane == "axis":
            marks.append(Mark(m, points=[(t, 0.0)]))
        else:
            marks.append(Mark(m, across(t, lo, hi) if plane == "across" else along(t, lo, hi)))
    return marks


def _rt_fronts():
    """Robinson and Trautman's fronts: the fronts of one retarded time differ only in size, and
    the embedding draws each with its own r as the unit, so a moment is the whole outgoing ray
    u = u_k of the plane, from r = 0 out."""
    return [Mark(m, [[(m.time, 0.0), (m.time, 100.0)]]) for m in moments("robinson_trautman", "fronts")]


def string_wave_V(u, X):
    """V of the moving string chart on the surface v = 0 of the isotropic chart, at the line X, for
    the pulse A = exp(-4u^2)/2: 2(X - A)A' + int_0^u A'^2, where A'^2 = 16u^2 exp(-8u^2) and its
    integral from 0 is sqrt(pi/8) erf(2 sqrt(2) u)/2 - u exp(-8u^2)."""
    A, dA = math.exp(-4 * u * u) / 2, -4 * u * math.exp(-4 * u * u)
    integral = math.sqrt(math.pi / 8) * math.erf(2 * math.sqrt(2) * u) / 2 - u * math.exp(-8 * u * u)
    return 2 * (X - A) * dA + integral


FLAT = {
    **{("kerr_de_sitter", system, view + suffix): (lambda sign=sign: _kds(sign))
       for sign, suffix in (("de_sitter", ""), ("anti_de_sitter", "_ads"))
       for system, view in (("boyer_lindquist", "axis"), ("boyer_lindquist", "principal"))},
    **{("kerr_de_sitter", system, "above" + suffix): (lambda sign=sign: _kds_above(sign))
       for sign, suffix in (("de_sitter", ""), ("anti_de_sitter", "_ads"))
       for system in ("boyer_lindquist", "nonrotating") if (system, suffix) != ("boyer_lindquist", "_ads")},
    **{("kerr_de_sitter", system, "axis" + suffix): (lambda sign=sign, scale=scale: _kds(sign, kds_rstar, scale))
       for sign, suffix in (("de_sitter", ""), ("anti_de_sitter", "_ads"))
       for system, scale in (("kerr_ingoing", 1), ("kerr_outgoing", -1))},
    **{("kerr_de_sitter", "kerr_schild", "axis" + suffix): (lambda sign=sign: _kds(sign, kds_schild))
       for sign, suffix in (("de_sitter", ""), ("anti_de_sitter", "_ads"))},
    ("kerr_taub_nut", "boyer_lindquist", "principal"): lambda: one("kerr_taub_nut", lambda m: along(0.0, *m.reach("boyer_lindquist", "r"))),
    ("kerr_taub_nut", "boyer_lindquist", "above"): lambda: kerr_above("kerr_taub_nut"),
    ("robinson_trautman", "axisymmetric", "axis"): _rt_fronts,
    ("robinson_trautman", "axisymmetric", "equator"): _rt_fronts,
    ("btz", "stationary", "static"): lambda: _btz(),
    ("btz", "eddington_finkelstein_ingoing", "static"): lambda: _btz(1),
    ("btz", "eddington_finkelstein_outgoing", "static"): lambda: _btz(-1),
    ("reissner_nordstrom_de_sitter", "static", "radial"): lambda: _rnds("static"),
    ("reissner_nordstrom_de_sitter", "eddington_finkelstein_ingoing", "finkelstein"): lambda: _rnds("ingoing"),
    ("reissner_nordstrom_de_sitter", "eddington_finkelstein_ingoing", "chart"): lambda: _rnds("ingoing"),
    ("reissner_nordstrom_de_sitter", "eddington_finkelstein_outgoing", "finkelstein"): lambda: _rnds("outgoing"),
    ("reissner_nordstrom_de_sitter", "eddington_finkelstein_outgoing", "chart"): lambda: _rnds("outgoing"),
    ("reissner_nordstrom_de_sitter", "cosmological", "plane"): lambda: _rnds("cosmological"),
    ("bardeen", "static", "radial"): lambda: _bardeen_both(),
    ("bardeen", "eddington_finkelstein_ingoing", "finkelstein"): lambda: _bardeen_both(1),
    ("bardeen", "eddington_finkelstein_ingoing", "chart"): lambda: _bardeen_both(1),
    ("bardeen", "eddington_finkelstein_outgoing", "finkelstein"): lambda: _bardeen_both(-1),
    ("bardeen", "eddington_finkelstein_outgoing", "chart"): lambda: _bardeen_both(-1),
    ("hayward", "static", "radial"): lambda: _hayward(),
    ("hayward", "eddington_finkelstein_ingoing", "finkelstein"): lambda: _hayward(1),
    ("hayward", "eddington_finkelstein_ingoing", "chart"): lambda: _hayward(1),
    ("hayward", "eddington_finkelstein_outgoing", "finkelstein"): lambda: _hayward(-1),
    ("hayward", "eddington_finkelstein_outgoing", "chart"): lambda: _hayward(-1),
    # v - r = T, every r the embedding reaches.
    ("hiscock", "ingoing", "history"): lambda: _hiscock("ingoing"),
    ("hiscock", "outgoing", "history"): lambda: _hiscock("outgoing"),
    ("hiscock", "flat", "after"): lambda: _hiscock("after"),
    ("hayward", "evaporating", "history"): lambda: one(
        "hayward", lambda m: [[(m.time + r, r) for r in m.reach("evaporating", "r")]], view_id="history"),
    ("reissner_nordstrom_ads", "static", "radial"): lambda: _rnads(),
    ("reissner_nordstrom_ads", "eddington_finkelstein_ingoing", "finkelstein"): lambda: _rnads(1),
    ("reissner_nordstrom_ads", "eddington_finkelstein_ingoing", "chart"): lambda: _rnads(1),
    ("reissner_nordstrom_ads", "eddington_finkelstein_outgoing", "finkelstein"): lambda: _rnads(-1),
    ("reissner_nordstrom_ads", "eddington_finkelstein_outgoing", "chart"): lambda: _rnads(-1),
    ("schwarzschild_ads", "static", "radial"): lambda: _sads(),
    ("schwarzschild_ads", "eddington_finkelstein_ingoing", "finkelstein"): lambda: _sads(1),
    ("schwarzschild_ads", "eddington_finkelstein_ingoing", "chart"): lambda: _sads(1),
    ("schwarzschild_ads", "eddington_finkelstein_outgoing", "finkelstein"): lambda: _sads(-1),
    ("schwarzschild_ads", "eddington_finkelstein_outgoing", "chart"): lambda: _sads(-1),
    ("ads_soliton", "horowitz_myers", "radial"): lambda: _soliton("horowitz_myers"),
    ("ads_soliton", "poincare", "tz"): lambda: _soliton("poincare"),
    ("ads_soliton", "polar", "radial"): lambda: _soliton("polar"),
    ("ads_soliton", "polar", "through"): lambda: _soliton("polar"),
    ("topological_black_hole", "static", "flat"): lambda: _tbh(),
    ("topological_black_hole", "black_string", "radial"): lambda: _tbh(),
    ("topological_black_hole", "eddington_finkelstein_ingoing", "flat"): lambda: _tbh(1),
    ("topological_black_hole", "eddington_finkelstein_outgoing", "flat"): lambda: _tbh(-1),
    ("topological_black_hole", "brane", "tz"): lambda: _tbh(brane=True),
    ("topological_black_hole", "static", "massless"): lambda: _tbh_horizon(),
    ("topological_black_hole", "hyperbolic", "massless"): lambda: _tbh_horizon(),
    ("schwarzschild", "spherical", "radial"): lambda: one("schwarzschild", lambda m: along(0.0, *m.reach("spherical", "r"))),
    ("schwarzschild", "eddington_finkelstein_ingoing", "finkelstein"): lambda: schwarzschild_t(1),
    ("schwarzschild", "eddington_finkelstein_ingoing", "chart"): lambda: schwarzschild_t(1),
    ("schwarzschild", "eddington_finkelstein_outgoing", "finkelstein"): lambda: schwarzschild_t(-1),
    ("schwarzschild", "eddington_finkelstein_outgoing", "chart"): lambda: schwarzschild_t(-1),
    ("ellis_bronnikov", "spherical", "radial"): lambda: one("ellis_bronnikov", lambda m: along(0.0, *m.reach("spherical", "r"))),
    ("morris_thorne", "spherical", "radial"): lambda: one("morris_thorne", lambda m: along(0.0, *m.reach("spherical", "r"))),
    # The gravastar's moment t = 0: inside the shell from the centre to it, in the areal radius and in
    # the tortoise coordinate x = (L/sqrt(C)) artanh(r/L) at L = 2 and C = 64/195, and outside it from the
    # shell as far as the embedding reaches.
    ("gravastar", "interior", "radial"): lambda: one("gravastar", lambda m: along(0.0, *m.reach("interior", "r"))),
    ("gravastar", "interior", "through"): lambda: one("gravastar", lambda m: along(0.0, *m.reach("interior", "r"))),
    ("gravastar", "interior_tortoise", "radial"): lambda: one(
        "gravastar", lambda m: along(0.0, *(gravastar_x(r) for r in m.reach("interior", "r")))),
    ("gravastar", "exterior", "radial"): lambda: one("gravastar", lambda m: along(0.0, *m.reach("exterior", "r"))),
    ("thin_shell_wormhole", "spherical", "radial"): lambda: one(
        "thin_shell_wormhole", lambda m: along(0.0, *m.reach("spherical", "r"))),
    # Both sides are read in the areal chart from the throat r = a out, and l = +-(r - a).
    ("thin_shell_wormhole", "throat", "radial"): lambda: one(
        "thin_shell_wormhole", lambda m: along(0.0, m.reach("spherical", "r")[0] - m.reach("spherical", "r")[1],
                                               m.reach("spherical", "r")[1] - m.reach("spherical", "r")[0])),
    **{("teo_wormhole", system, view): (lambda system=system, view=view: teo(system, view))
       for system in ("spherical", "proper_radial") for view in ("axis", "equator")},
    # Damour and Solodukhin's wormhole at r_s = 1 and lambda = 1/5: one side in its own chart and in the
    # rescaled time, from the throat r = r_s out, and both sides through the throat, where
    # cosh rho = 52 r - 51, the isotropic radius runs from r_s^2/16r to r, and u = +-sqrt(r - r_s).
    ("damour_solodukhin", "spherical", "radial"): lambda: one(
        "damour_solodukhin", lambda m: along(0.0, *m.reach("spherical", "r"))),
    ("damour_solodukhin", "rescaled", "radial"): lambda: one(
        "damour_solodukhin", lambda m: along(0.0, *m.reach("spherical", "r"))),
    ("damour_solodukhin", "throat", "radial"): lambda: one(
        "damour_solodukhin", lambda m: along(0.0, -math.acosh(52 * m.reach("spherical", "r")[1] - 51),
                                             math.acosh(52 * m.reach("spherical", "r")[1] - 51))),
    ("damour_solodukhin", "isotropic", "radial"): lambda: one(
        "damour_solodukhin", lambda m: along(0.0, 1 / (16 * _ds_isotropic(m.reach("spherical", "r")[1])),
                                             _ds_isotropic(m.reach("spherical", "r")[1]))),
    ("damour_solodukhin", "einstein_rosen", "radial"): lambda: one(
        "damour_solodukhin", lambda m: along(0.0, -math.sqrt(m.reach("spherical", "r")[1] - 1),
                                             math.sqrt(m.reach("spherical", "r")[1] - 1))),
    # Simpson and Visser's three geometries, each on its own plane in each of the four charts.
    **{("simpson_visser", chart, case): (lambda chart=chart, case=case: simpson_visser(chart, case))
       for chart in ("spherical", "areal", "eddington_finkelstein_ingoing", "eddington_finkelstein_outgoing")
       for case in SV_A},
    # Roberts's collapse: the moments of his time t of each outcome, on that outcome's plane in each chart.
    **{("roberts", chart, case): (lambda chart=chart, case=case: roberts(chart, case))
       for chart in ("double_null", "advanced", "areal", "diagonal", "scaling") for case in ROBERTS_P},
    # Fisher, Janis, Newman and Winicour's scalar field at gamma = 1/2 and b = 1: the moment t = 0 from
    # the singularity out, in Wyman's r, in R = r - 3b/4, in the isotropic radius, and in the harmonic
    # coordinate, which its plane draws in units of 1/k at k = 1, so at half the embedding's own u.
    ("fisher_jnw", "spherical", "radial"): lambda: one("fisher_jnw", lambda m: along(0.0, *_fjnw_reach(m, lambda r: r))),
    ("fisher_jnw", "jnw", "radial"): lambda: one("fisher_jnw", lambda m: along(0.0, *_fjnw_reach(m, lambda r: r - 0.75))),
    ("fisher_jnw", "isotropic", "radial"): lambda: one(
        "fisher_jnw", lambda m: along(0.0, *_fjnw_reach(m, lambda r: (r - 0.5 + math.sqrt(r * (r - 1))) / 2))),
    ("fisher_jnw", "harmonic", "radial"): lambda: one(
        "fisher_jnw", lambda m: along(0.0, *(u / 2 for u in m.reach("harmonic", "u")))),
    ("hartle_thorne", "hartle_thorne", "equator"): lambda: one(
        "hartle_thorne", lambda m: along(0.0, *m.reach("hartle_thorne", "r"))),
    # Witten's black hole in two dimensions: the moment t = 0 outside the horizon in each chart, and
    # on both sides of it on the Kruskal plane.
    **{("witten_black_hole", chart, view): (lambda chart=chart: witten(chart))
       for chart, view in (("witten", "radial"), ("schwarzschild_gauge", "radial"), ("dilaton", "radial"),
                           ("conformal", "radial"), ("kruskal", "plane"),
                           ("eddington_finkelstein_ingoing", "finkelstein"),
                           ("eddington_finkelstein_outgoing", "finkelstein"))},
    # Einstein and Rosen's bridges, three spacetimes of one page, each drawing marked with the moment of its
    # own: the neutral bridge at r_s = 1, read in Schwarzschild's r, where u = +-sqrt(r - r_s) and the
    # isotropic radius runs from r_s^2/16r to r; the charged bridge with a mass in its areal radius;
    # and the charged bridge with no mass in their u.
    ("einstein_rosen_bridge", "bridge", "radial"): lambda: one(
        "einstein_rosen_bridge", lambda m: along(0.0, -math.sqrt(m.reach("spherical", "r")[1] - 1),
                                                 math.sqrt(m.reach("spherical", "r")[1] - 1)), view_id="neutral"),
    ("einstein_rosen_bridge", "spherical", "radial"): lambda: one(
        "einstein_rosen_bridge", lambda m: along(0.0, *m.reach("spherical", "r")), view_id="neutral"),
    ("einstein_rosen_bridge", "isotropic", "radial"): lambda: one(
        "einstein_rosen_bridge", lambda m: along(0.0, 1 / (16 * _ds_isotropic(m.reach("spherical", "r")[1])),
                                                 _ds_isotropic(m.reach("spherical", "r")[1])), view_id="neutral"),
    ("einstein_rosen_bridge", "charged_spherical", "radial"): lambda: one(
        "einstein_rosen_bridge", lambda m: along(0.0, *m.reach("charged_spherical", "r")), view_id="charged_mass"),
    ("einstein_rosen_bridge", "charged_bridge", "radial"): lambda: one(
        "einstein_rosen_bridge", lambda m: along(0.0, *m.reach("charged_bridge", "u")), view_id="charged"),
    ("minkowski", "spherical", "radial"): lambda: one("minkowski", lambda m: along(0.0, *m.reach("spherical", "r"))),
    # t = (u + v)/2 and r = (v - u)/2, so the moment is u = -r, v = r.
    ("minkowski", "spherical_null", "radial"): lambda: one(
        "minkowski", lambda m: [[(-r, r) for r in m.reach("spherical", "r")]]),
    ("minkowski", "cartesian", "tx"): lambda: one("minkowski", lambda m: across(0.0, 0.0, m.reach("spherical", "r")[1])),
    # ct = X sinh(aT/c) vanishes in the wedge only at T = 0, where x = X.
    ("minkowski", "rindler", "tx"): lambda: one("minkowski", lambda m: along(0.0, 0.0, m.reach("spherical", "r")[1])),
    # Kiselev's black hole at w = -2/3 on the charts that hold it, and his matter alone, another
    # spacetime, on his two charts of it.
    ("kiselev", "static", "radial"): lambda: one(
        "kiselev", lambda m: along(0.0, *m.reach("linear", "r")), view_id="black_hole"),
    ("kiselev", "linear", "radial"): lambda: one(
        "kiselev", lambda m: along(0.0, *m.reach("linear", "r")), view_id="black_hole"),
    ("kiselev", "eddington_finkelstein_ingoing", "finkelstein"): lambda: kiselev_t(1),
    ("kiselev", "eddington_finkelstein_ingoing", "chart"): lambda: kiselev_t(1),
    ("kiselev", "eddington_finkelstein_outgoing", "finkelstein"): lambda: kiselev_t(-1),
    ("kiselev", "eddington_finkelstein_outgoing", "chart"): lambda: kiselev_t(-1),
    ("kiselev", "hyperbolic", "radial"): lambda: kiselev_free(False),
    ("kiselev", "conformally_flat", "radial"): lambda: kiselev_free(True),
    ("schwarzschild_de_sitter", "static", "radial"): lambda: one("schwarzschild_de_sitter", lambda m: along(0.0, *m.reach("static", "r"))),
    ("schwarzschild_de_sitter", "eddington_finkelstein_ingoing", "finkelstein"): lambda: kottler_t(1),
    ("schwarzschild_de_sitter", "eddington_finkelstein_ingoing", "chart"): lambda: kottler_t(1),
    ("schwarzschild_de_sitter", "eddington_finkelstein_outgoing", "finkelstein"): lambda: kottler_t(-1),
    ("schwarzschild_de_sitter", "eddington_finkelstein_outgoing", "chart"): lambda: kottler_t(-1),
    # Tangherlini's black hole in five dimensions on its three charts, and in six on its static chart;
    # each dimension is another spacetime, marked on its own drawings alone.
    ("tangherlini", "spherical", "radial"): lambda: one(
        "tangherlini", lambda m: along(0.0, *m.reach("spherical", "r")), view_id="five"),
    ("tangherlini", "eddington_finkelstein_ingoing", "finkelstein"): lambda: tangherlini_t(1),
    ("tangherlini", "eddington_finkelstein_ingoing", "chart"): lambda: tangherlini_t(1),
    ("tangherlini", "eddington_finkelstein_outgoing", "finkelstein"): lambda: tangherlini_t(-1),
    ("tangherlini", "eddington_finkelstein_outgoing", "chart"): lambda: tangherlini_t(-1),
    ("tangherlini", "spherical_six", "radial"): lambda: one(
        "tangherlini", lambda m: along(0.0, *m.reach("spherical_six", "r")), view_id="six"),
    ("boulware_deser", "spherical", "radial"): lambda: one(
        "boulware_deser", lambda m: along(0.0, *m.reach("spherical", "r")), view_id="hole"),
    ("boulware_deser", "eddington_finkelstein_ingoing", "finkelstein"): lambda: boulware_deser_t(1),
    ("boulware_deser", "eddington_finkelstein_ingoing", "chart"): lambda: boulware_deser_t(1),
    ("boulware_deser", "eddington_finkelstein_outgoing", "finkelstein"): lambda: boulware_deser_t(-1),
    ("boulware_deser", "eddington_finkelstein_outgoing", "chart"): lambda: boulware_deser_t(-1),
    ("boulware_deser", "spherical_plus", "radial"): lambda: one(
        "boulware_deser", lambda m: along(0.0, *m.reach("spherical_plus", "r")), view_id="branch"),
    # The black string across the string, in five dimensions on its three charts and in six on its
    # static chart; five and six are two spacetimes, each marked on its own charts alone.
    ("black_string", "static", "radial"): lambda: one(
        "black_string", lambda m: along(0.0, *m.reach("static", "r")), view_id="across"),
    ("black_string", "eddington_finkelstein_ingoing", "finkelstein"): lambda: black_string_t(),
    ("black_string", "eddington_finkelstein_ingoing", "chart"): lambda: black_string_t(),
    ("black_string", "kerr_schild", "radial"): lambda: black_string_t(kerr_schild=True),
    ("black_string", "static_six", "radial"): lambda: one(
        "black_string", lambda m: along(0.0, *m.reach("static_six", "r")), view_id="six"),
    # Myers and Perry's black hole: each plane of t and r marks the moment embedded on its own surface,
    # the plane transverse to the rotation or the plane of rotation in five dimensions, the Hopf fibre of
    # the hole with equal spins, and the transverse plane in six dimensions.
    ("myers_perry", "boyer_lindquist", "transverse"): lambda: one(
        "myers_perry", lambda m: along(0.0, *m.reach("boyer_lindquist", "r")), view_id="transverse"),
    ("myers_perry", "boyer_lindquist", "rotation"): lambda: one(
        "myers_perry", lambda m: along(0.0, *m.reach("boyer_lindquist", "r")), view_id="rotation"),
    ("myers_perry", "ingoing_kerr", "transverse"): lambda: myers_perry_t("transverse"),
    ("myers_perry", "ingoing_kerr", "rotation"): lambda: myers_perry_t("rotation"),
    ("myers_perry", "equal_spins", "radial"): lambda: one(
        "myers_perry", lambda m: along(0.0, *m.reach("equal_spins", "\\rho")), view_id="fibre"),
    ("myers_perry", "boyer_lindquist_six", "transverse"): lambda: one(
        "myers_perry", lambda m: along(0.0, *m.reach("boyer_lindquist_six", "r")), view_id="six"),
    # Black Saturn's plane of the ring at t = 0, on each plane of Weyl's chart and of the polar chart.
    ("black_saturn", "weyl", "outside"): lambda: black_saturn_t("outside"),
    ("black_saturn", "weyl", "between"): lambda: black_saturn_t("between"),
    ("black_saturn", "polar", "far"): lambda: black_saturn_t("far"),
    # The Kaluza-Klein monopole's cigar, the half axis theta = 0 at t = 0, on its plane of t and the radius
    # in each chart; the Taub-NUT radius is rho = r + 2m.
    ("kaluza_klein_monopole", "gross_perry", "radial"): lambda: one(
        "kaluza_klein_monopole", lambda m: along(0.0, *m.reach("gross_perry", "r"))),
    ("kaluza_klein_monopole", "gross_perry", "through"): lambda: one(
        "kaluza_klein_monopole", lambda m: along(0.0, *m.reach("gross_perry", "r"))),
    ("kaluza_klein_monopole", "hopf", "radial"): lambda: one(
        "kaluza_klein_monopole", lambda m: along(0.0, *m.reach("gross_perry", "r"))),
    ("kaluza_klein_monopole", "taub_nut", "radial"): lambda: one(
        "kaluza_klein_monopole", lambda m: along(0.0, *(r + 2 for r in m.reach("gross_perry", "r")))),
    **{("kaluza_klein_black_hole", system, "radial"): lambda system=system: one(
        "kaluza_klein_black_hole", lambda m: along(0.0, *m.reach(system, "r")), view_id=system)
       for system in ("electric", "magnetic", "einstein")},
    ("kaluza_klein_black_hole", "eddington_finkelstein_ingoing", "finkelstein"): lambda: kkbh_t("electric"),
    ("kaluza_klein_black_hole", "einstein_eddington_finkelstein", "finkelstein"): lambda: kkbh_t("einstein"),
    # Letelier's black hole, r_s = 1, on the static and Eddington-Finkelstein planes, and the monopole
    # with no mass at its centre on the Barriola-Vilenkin plane; each is another spacetime than the other.
    ("string_black_hole", "static", "radial"): lambda: string_hole(),
    ("string_black_hole", "wedge", "radial"): lambda: string_hole(),
    ("string_black_hole", "eddington_finkelstein_ingoing", "finkelstein"): lambda: string_hole(1),
    ("string_black_hole", "eddington_finkelstein_ingoing", "chart"): lambda: string_hole(1),
    ("string_black_hole", "eddington_finkelstein_outgoing", "finkelstein"): lambda: string_hole(-1),
    ("string_black_hole", "eddington_finkelstein_outgoing", "chart"): lambda: string_hole(-1),
    ("global_monopole", "static", "radial"): lambda: one(
        "global_monopole", lambda m: along(0.0, *m.reach("static", "r")), view_id="black_hole"),
    ("global_monopole", "conical", "radial"): lambda: one(
        "global_monopole", lambda m: along(0.0, *m.reach("conical", "r")), view_id="monopole"),
    ("global_monopole", "eddington_finkelstein_ingoing", "finkelstein"): lambda: monopole_t(1),
    ("global_monopole", "eddington_finkelstein_ingoing", "chart"): lambda: monopole_t(1),
    ("global_monopole", "eddington_finkelstein_outgoing", "finkelstein"): lambda: monopole_t(-1),
    ("global_monopole", "eddington_finkelstein_outgoing", "chart"): lambda: monopole_t(-1),
    # The dilaton black hole at r_d = r_s/2: the Einstein metric's moment on the static and
    # Eddington-Finkelstein planes, and each string metric's moment on its own plane.
    ("dilaton_black_hole", "static", "radial"): lambda: one(
        "dilaton_black_hole", lambda m: along(0.0, *m.reach("static", "r")), view_id="einstein"),
    ("dilaton_black_hole", "eddington_finkelstein_ingoing", "finkelstein"): lambda: dilaton_t(1),
    ("dilaton_black_hole", "eddington_finkelstein_ingoing", "chart"): lambda: dilaton_t(1),
    ("dilaton_black_hole", "eddington_finkelstein_outgoing", "finkelstein"): lambda: dilaton_t(-1),
    ("dilaton_black_hole", "eddington_finkelstein_outgoing", "chart"): lambda: dilaton_t(-1),
    ("dilaton_black_hole", "string_magnetic", "radial"): lambda: one(
        "dilaton_black_hole", lambda m: along(0.0, *m.reach("string_magnetic", "r")), view_id="string_magnetic"),
    ("dilaton_black_hole", "string_electric", "radial"): lambda: one(
        "dilaton_black_hole", lambda m: along(0.0, *m.reach("string_electric", "r")), view_id="string_electric"),
    # The moments ct = T of the Einstein-Rosen pulse, out to where the embedding reaches, and in the
    # null chart u = T - rho, v = T + rho along the same stretch.
    ("einstein_rosen_waves", "cylindrical", "radial"): lambda: one(
        "einstein_rosen_waves", lambda m: along(m.time, *m.reach("cylindrical", "\\rho"))),
    ("einstein_rosen_waves", "null", "radial"): lambda: one(
        "einstein_rosen_waves", lambda m: [[(m.time - r, m.time + r) for r in m.reach("cylindrical", "\\rho")]]),
    ("de_sitter", "static_spherical", "radial"): lambda: one("de_sitter", lambda m: along(0.0, *m.reach("static_spherical", "r"))),
    ("de_sitter", "static_spherical", "through"): lambda: one("de_sitter", lambda m: along(0.0, *m.reach("static_spherical", "r"))),
    ("de_sitter", "flat_slicing", "tx"): _ds_flat,
    # The moment t = 0 runs from the pole to the antipode, chi from 0 to pi, and the areal chart
    # carries its near hemisphere, r = R sin chi from 0 to R.
    ("einstein_static", "hyperspherical", "radial"): lambda: one("einstein_static", lambda m: along(0.0, *m.reach("hyperspherical", "\\chi"))),
    ("einstein_static", "hyperspherical", "through"): lambda: one("einstein_static", lambda m: along(0.0, *m.reach("hyperspherical", "\\chi"))),
    ("einstein_static", "static_areal", "radial"): lambda: one("einstein_static", _es_areal),
    # A moment ct of the Milne universe reaches chi from 0 to its edge: r = sinh chi, c tau = ln(ct)
    # in units of ct_0, and in the inertial chart the hyperbola c^2T^2 - R^2 = c^2t^2 across the centre.
    ("milne", "comoving_hyperbolic", "through"): lambda: one("milne", lambda m: across(m.time, 0.0, m.reach("comoving_hyperbolic", "\\chi")[1])),
    ("milne", "comoving_spherical", "radial"): lambda: one("milne", lambda m: along(m.time, 0.0, math.sinh(m.reach("comoving_hyperbolic", "\\chi")[1]))),
    ("milne", "logarithmic_time", "radial"): lambda: one("milne", lambda m: along(math.log(m.time), *m.reach("comoving_hyperbolic", "\\chi"))),
    ("milne", "inertial", "through"): lambda: one("milne", _milne_inertial),
    # A moment c tau of the open universe inside the anti-de Sitter bubble, l = 1, reaches chi from 0
    # to its edge: eta = ln tan(c tau/2) in the conformal time, and in the static chart inside the
    # wall the curve sqrt(1 + r^2) cos(ct) = cos(c tau) with r = sin(c tau) sinh(chi).
    ("coleman_de_luccia", "open", "negative"): lambda: one("coleman_de_luccia", lambda m: along(m.time, *m.reach("open", "\\chi"))),
    ("coleman_de_luccia", "open_conformal", "negative"): lambda: one(
        "coleman_de_luccia", lambda m: along(math.log(math.tan(m.time / 2)), *m.reach("open", "\\chi"))),
    ("coleman_de_luccia", "static_inside", "out_of_flat"): lambda: one("coleman_de_luccia", _cdl_static),
    # A moment kct of the domain wall's global chart meets the plane x = y = 0 of the planar chart
    # along t = const, every z, and in the inertial chart of the side z < 0 it is the cone
    # cT = R tanh(kct) from the centre at T = 0 out to the wall, across the centre on either side.
    # Randall and Sundrum's moment t = 0, out to the embedding's reach in the proper distance y: in the
    # conformally flat chart w = sgn(y)(e^{k|y|} - 1)/k, on one side in the Poincare chart z = e^{ky}/k, at
    # k = 1; between two walls phi and -phi are one point, so the moment is the whole line of phi.
    ("randall_sundrum", "proper_distance", "ty"): lambda: one(
        "randall_sundrum", lambda m: along(0.0, *m.reach("proper_distance", "y")), view_id="pseudosphere"),
    ("randall_sundrum", "conformal", "tw"): lambda: one(
        "randall_sundrum", lambda m: along(0.0, *(math.copysign(math.expm1(abs(y)), y) for y in m.reach("proper_distance", "y"))),
        view_id="pseudosphere"),
    ("randall_sundrum", "poincare", "tz"): lambda: one(
        "randall_sundrum", lambda m: along(0.0, 1.0, math.exp(m.reach("proper_distance", "y")[1])), view_id="pseudosphere"),
    ("randall_sundrum", "two_walls", "tphi"): lambda: one(
        "randall_sundrum", lambda m: along(0.0, -m.reach("two_walls", "\\phi")[1], m.reach("two_walls", "\\phi")[1]),
        view_id="two_walls"),
    # Lifshitz spacetime's moment t = 0, over the embedding's reach in the proper distance rho, at z = 2
    # and L = 1: r = e^rho, u = e^-rho, w = e^(-2 rho)/2 and s = e^(2 rho), the curve v = -w in the two
    # null charts, and on each plane of t and x the line t = 0 over the strip 0 <= x < 2 pi L.
    ("lifshitz_spacetime", "kachru_liu_mulligan", "tr"): lambda: one(
        "lifshitz_spacetime", lambda m: along(0.0, *(math.exp(y) for y in m.reach("proper_distance", "\\rho")))),
    ("lifshitz_spacetime", "poincare", "tu"): lambda: one(
        "lifshitz_spacetime", lambda m: along(0.0, *sorted(math.exp(-y) for y in m.reach("proper_distance", "\\rho")))),
    **{("lifshitz_spacetime", "poincare", view): lambda: one("lifshitz_spacetime", lambda m: along(0.0, 0.0, 2 * math.pi))
       for view in ("tx_half", "tx_one", "tx_two")},
    ("lifshitz_spacetime", "proper_distance", "trho"): lambda: one(
        "lifshitz_spacetime", lambda m: along(0.0, *m.reach("proper_distance", "\\rho"))),
    ("lifshitz_spacetime", "tortoise", "tw"): lambda: one(
        "lifshitz_spacetime", lambda m: along(0.0, *sorted(0.5 * math.exp(-2 * y) for y in m.reach("proper_distance", "\\rho")))),
    ("lifshitz_spacetime", "eddington_finkelstein", "vr"): lambda: one(
        "lifshitz_spacetime", lambda m: [np.column_stack([-0.5 * np.exp(-2 * rho), np.exp(rho)])
                                         for rho in [np.linspace(*m.reach("proper_distance", "\\rho"), N)]]),
    ("lifshitz_spacetime", "affine", "vs"): lambda: one(
        "lifshitz_spacetime", lambda m: [np.column_stack([-0.5 * np.exp(-2 * rho), np.exp(2 * rho)])
                                         for rho in [np.linspace(*m.reach("proper_distance", "\\rho"), N)]]),
    ("domain_wall", "planar", "tz"): lambda: one("domain_wall", lambda m: across(m.time, 0.0, 1.0)),
    ("domain_wall", "inertial", "through"): lambda: one("domain_wall", _wall_inertial),
    ("anti_de_sitter", "static_global", "radial"): lambda: one("anti_de_sitter", lambda m: along(0.0, *m.reach("static_global", "r"))),
    ("anti_de_sitter", "static_global", "through"): lambda: one("anti_de_sitter", lambda m: along(0.0, *m.reach("static_global", "r"))),
    ("anti_de_sitter", "poincare", "tx"): _ads_poincare,
    ("rn_metric", "spherical", "radial"): lambda: _rn("outside") + _rn("inside"),
    # One hole's equator lies on the plane of t and r about it; the two holes' midplane z = 0 is
    # the plane of t and x through its centre, the plane of t and rho out from the axis, and one
    # event of the axis.
    ("majumdar_papapetrou", "isotropic", "radial"): lambda: one(
        "majumdar_papapetrou", lambda m: along(0.0, *m.reach("isotropic", "r")), view_id="one_hole"),
    ("majumdar_papapetrou", "cartesian", "tx"): lambda: one(
        "majumdar_papapetrou", lambda m: across(0.0, *m.reach("cylindrical", "\\rho")), view_id="two_holes"),
    ("majumdar_papapetrou", "cylindrical", "radial"): lambda: one(
        "majumdar_papapetrou", lambda m: along(0.0, *m.reach("cylindrical", "\\rho")), view_id="two_holes"),
    ("majumdar_papapetrou", "cartesian", "tz"): _mp_axis,
    ("israel_wilson_perjes", "cylindrical", "midplane"): lambda: one(
        "israel_wilson_perjes", lambda m: along(0.0, *m.reach("cylindrical", "\\rho")), view_id="two_sources"),
    ("israel_wilson_perjes", "spheroidal", "axis"): lambda: one(
        "israel_wilson_perjes", lambda m: along(0.0, *m.reach("spheroidal", "r")), view_id="spinning"),
    ("israel_wilson_perjes", "spheroidal", "principal"): lambda: one(
        "israel_wilson_perjes", lambda m: along(0.0, *m.reach("spheroidal", "r")), view_id="spinning"),
    ("israel_wilson_perjes", "spherical", "radial"): lambda: one(
        "israel_wilson_perjes", lambda m: along(0.0, *m.reach("spherical", "r")), view_id="charged_nut"),
    ("kastor_traschen", "cartesian", "tz"): lambda: _kt("cartesian", "axis"),
    ("kastor_traschen", "cartesian", "tx"): lambda: _kt("cartesian", "across"),
    ("kastor_traschen", "cylindrical", "radial"): lambda: _kt("cylindrical", "along"),
    ("kastor_traschen", "isotropic", "radial"): lambda: _kt("isotropic", None),
    ("kastor_traschen", "comoving", "tz"): lambda: _kt("comoving", "axis"),
    ("kastor_traschen", "comoving", "tx"): lambda: _kt("comoving", "across"),
    ("taub_nut", "spherical", "radial"): lambda: one("taub_nut", lambda m: along(0.0, *m.reach("spherical", "r"))),
    ("near_horizon_extreme_kerr", "poincare", "equator"): lambda: _nhek("poincare"),
    ("near_horizon_extreme_kerr", "inverse_radius", "equator"): lambda: _nhek("inverse_radius"),
    ("near_horizon_extreme_kerr", "global", "equator"): lambda: _nhek("global"),
    ("bertotti_robinson", "static", "radial"): lambda: _br("static"),
    ("bertotti_robinson", "poincare", "tx"): lambda: _br("poincare"),
    ("nariai", "static", "patch"): lambda: _nariai("static"),
    ("nariai", "global", "circle"): lambda: _nariai("global"),
    ("interior_schwarzschild", "spherical", "radial"): lambda: one("interior_schwarzschild", lambda m: along(0.0, *m.reach("spherical", "r"))),
    ("interior_schwarzschild", "spherical", "through"): lambda: one("interior_schwarzschild", lambda m: along(0.0, *m.reach("spherical", "r"))),
    # On the axis the moment meets the plane of t and r off its equator, along the whole axis
    # outside r_+; on the equator it is the equator itself.
    ("kerr", "boyer_lindquist", "radial"): lambda: one("kerr", lambda m: along(0.0, *m.reach("boyer_lindquist", "r"))),
    ("kerr", "boyer_lindquist", "principal"): lambda: one("kerr", lambda m: along(0.0, *m.reach("boyer_lindquist", "r"))),
    ("kerr", "boyer_lindquist", "above"): lambda: kerr_above("kerr"),
    ("kerr_newman", "boyer_lindquist", "radial"): lambda: one("kerr_newman", lambda m: along(0.0, *m.reach("boyer_lindquist", "r"))),
    ("kerr_newman", "boyer_lindquist", "principal"): lambda: one("kerr_newman", lambda m: along(0.0, *m.reach("boyer_lindquist", "r"))),
    ("kerr_newman", "boyer_lindquist", "above"): lambda: kerr_above("kerr_newman"),
    # Homogeneous planes: every moment runs across the whole drawing.
    # A moment of the areal time t is the line tau = -ln t of the logarithmic chart.
    ("gowdy", "areal", "plane"): lambda: one("gowdy", lambda m: along(m.time, *m.reach("areal", "\\theta"))),
    ("gowdy", "logarithmic", "plane"): lambda: one(
        "gowdy", lambda m: along(-math.log(m.time), *m.reach("areal", "\\theta"))),
    ("kasner", "cartesian", "tx"): lambda: one("kasner", lambda m: across(m.time, 0.0, BIG)),
    ("kasner", "cartesian", "tz"): lambda: one("kasner", lambda m: across(m.time, 0.0, BIG)),
    ("bianchi", "type_i_cartesian", "tx"): lambda: one("bianchi", lambda m: across(m.time, 0.0, BIG)),
    ("kantowski_sachs", "comoving", "tr"): lambda: _kantowski_sachs("comoving"),
    ("kantowski_sachs", "dust", "etar"): lambda: _kantowski_sachs("dust"),
    ("kantowski_sachs", "schwarzschild_interior", "Tr"): lambda: _kantowski_sachs("schwarzschild_interior"),
    # Godel's t_x = 2t + sqrt(2)(2 arctan(e^(-2r) tan(phi/2)) - phi) and e^x = cosh 2r + cos(phi) sinh 2r
    # put the plane y = 0 at phi = 0 and pi, where t_x = 2t and x = +-2r.
    ("godel", "cartesian", "tx"): lambda: one("godel", lambda m: across(0.0, 0.0, 2 * m.reach("cylindrical", "r")[1])),
    ("godel", "cylindrical", "inside"): lambda: one("godel", _godel_cylinder(math.asinh(1.0) / 2)),
    ("stockum_dust", "cylindrical", "inside"): lambda: one("stockum_dust", _godel_cylinder(0.5)),
    # Som and Raychaudhuri's plane y = 0 is phi = 0 and pi with the same t, where x = +-r.
    ("som_raychaudhuri", "cartesian", "tx"): lambda: one(
        "som_raychaudhuri", lambda m: across(0.0, 0.0, m.reach("cylindrical", "r")[1])),
    ("som_raychaudhuri", "cylindrical", "inside"): lambda: one("som_raychaudhuri", _godel_cylinder(0.5)),
    # The spinning string's moment t = 0 outside the null circle: the whole line t = 0 of each
    # cylinder outside it, and in the helical chart c tau = a phi~/b, one turn of the helix.
    **{("spinning_string", system, "outside"): lambda R=R: one("spinning_string", _spinning_cylinder(R))
       for system, R in (("proper_radius", 0.9 * math.sqrt(1.25)), ("rescaled_radius", 0.9 * math.sqrt(1.25)),
                         ("circumference_radius", 0.9))},
    ("spinning_string", "helical", "outside"): lambda: one(
        "spinning_string", _spinning_cylinder(0.9 * math.sqrt(1.25), [[(0.0, 0.0), (2 * math.pi * 0.9, 2 * math.pi * 0.9)]])),
    ("alcubierre", "cartesian", "tx"): lambda: one("alcubierre", lambda m: across(0.0, 0.0, m.grid()["u"][-1])),
    ("natario", "cartesian_flow", "tx"): lambda: one("natario", lambda m: across(0.0, 0.0, m.grid()["u"][-1])),
    ("krasnikov", "cylindrical", "tx"): _krasnikov,
    # Van Den Broeck's pocket at t = 0, read in the proper distance l: the comoving radius of the
    # circle at l is vdb_radius(l), which is also |x| on the axis of the Cartesian chart at t = 0.
    ("van_den_broeck", "proper_radial", "radial"): lambda: one(
        "van_den_broeck", lambda m: along(0.0, *m.reach("proper_radial", "l"))),
    ("van_den_broeck", "pocket", "radial"): lambda: one(
        "van_den_broeck", lambda m: along(0.0, *(vdb_radius(l) for l in m.reach("proper_radial", "l")))),
    ("van_den_broeck", "cartesian", "tx"): lambda: one(
        "van_den_broeck", lambda m: across(0.0, 0.0, vdb_radius(m.reach("proper_radial", "l")[1]))),
    # Misner space's moments are the Milne chart's t = t_k, every chi, and in Misner's
    # coordinates T = -c^2t_k^2/4, every psi; checks() carries the one chart onto the other.
    ("misner", "misner", "plane"): lambda: one("misner", lambda m: across(-m.time * m.time / 4, 0.0, BIG)),
    ("misner", "milne", "plane"): lambda: one("misner", lambda m: across(m.time, 0.0, BIG)),
    # Gott's moments are Grant's Milne time tau = tau_k, every chi.
    ("gott_time_machine", "grant_milne", "plane"): lambda: one("gott_time_machine", lambda m: across(m.time, 0.0, BIG)),
    # Ori's moments are his t = t_k on the slice y = 0, every z. The central circle has T = t, and
    # the circle at x = 4 l has T = t + (a/2 - e) x^2 = t - 3/2 at a = 1/16 and e = 1/8; on the
    # Brinkmann plane through the central circle T = -uv/2, so the moment is the hyperbola
    # uv = -2 t_k with both negative. checks() carries each chart onto the vacuum core's.
    ("ori_time_machine", "foliation", "centre"): lambda: one("ori_time_machine", lambda m: across(m.time, 0.0, BIG)),
    ("ori_time_machine", "foliation", "off_centre"): lambda: one("ori_time_machine", lambda m: across(m.time, 0.0, BIG)),
    ("ori_time_machine", "vacuum_core", "centre"): lambda: one("ori_time_machine", lambda m: across(m.time, 0.0, BIG)),
    ("ori_time_machine", "vacuum_core", "off_centre"): lambda: one(
        "ori_time_machine", lambda m: across(m.time - 1.5, 0.0, BIG)),
    ("ori_time_machine", "brinkmann", "plane"): lambda: one("ori_time_machine", _ori_hyperbola),
    # The travelling wave on a string: each moment is the cone of constant u = u_k and v = 0, which
    # a plane of fixed place across the string meets at one event, (u_k, 0) in the null conical and
    # isotropic charts. In the moving string chart V = v + 2xA' + int_0^u A'^2 with x = X - A(u), so
    # the event is (u_k, 2(X - A)A' + int_0^(u_k) A'^2) for the declared pulse A = exp(-4u^2)/2.
    **{("string_wave", system, view): lambda: [Mark(m, points=[(m.time, 0.0)]) for m in moments("string_wave")]
       for system, view in (("null_conical", "toward"), ("null_conical", "away"), ("isotropic", "beside"))},
    **{("string_wave", "moving_string", view): lambda X=X: [
        Mark(m, points=[(m.time, string_wave_V(m.time, X))]) for m in moments("string_wave")]
       for view, X in (("behind", -0.25), ("ahead", 0.75))},
    # Siklos's wave front u = v = 0: an event on each plane of u and v, and the line of zero v, t, U
    # or V on each of Kaigorodov's planes, with x = e^(-rho) = e^(2Z) in units of L.
    **{("siklos", system, view): lambda: _siklos()
       for system, view in (("siklos", "near"), ("siklos", "far"), ("ozsvath_robinson_rozga", "centre"))},
    **{("siklos", system, "depth"): lambda: _siklos(lambda x: x)
       for system in ("kaigorodov", "kaigorodov_poincare", "kaigorodov_kundt")},
    ("siklos", "kaigorodov_horospheric", "depth"): lambda: _siklos(lambda x: -math.log(x)),
    ("siklos", "kaigorodov_homogeneous", "depth"): lambda: _siklos(lambda x: math.log(x) / 2),
    # The fronts of Kundt's kind with a cosmological constant: one event on the plane of u and v of
    # the same member of the family.
    ("kundt_waves", "ozsvath_robinson_rozga", "de_sitter"): lambda: _kundt("sphere"),
    ("kundt_waves", "ozsvath_robinson_rozga", "anti_de_sitter"): lambda: _kundt("hyperbolic"),
    # Schrodinger spacetime's plane of x and r: one event on each plane of the time and the null coordinate.
    **{("schrodinger_spacetime", system, view): lambda: _schrodinger()
       for system, views in (("poincare", ("near", "middle", "far")), ("inverse_radius", ("near", "middle", "far")),
                             ("global", ("near", "middle", "far")),
                             ("dynamical_exponent", ("one", "three_halves", "three")),
                             ("poincare_5d", ("middle",)), ("poincare_6d", ("middle",)))
       for view in views},
    # The wave front u = u_k, every v.
    ("pp_wave", "exact_plane_wave", "tz"): lambda: one("pp_wave", lambda m: [[(m.time, -BIG), (m.time, BIG)]]),
    **{("aichelburg_sexl", "null_cartesian", view): lambda: one("aichelburg_sexl", lambda m: [[(m.time, -BIG), (m.time, BIG)]])
       for view in ("half", "eighth", "thirtysecond")},
    # Bonnor's wave front u = u_k, every v, with sqrt(2) u = ct - z in the Cartesian chart.
    **{("light_beam", "cartesian", view): lambda: one(
        "light_beam", lambda m: [[(-BIG, -BIG - math.sqrt(2) * m.time), (BIG, BIG - math.sqrt(2) * m.time)]])
       for view in ("axis", "beside")},
    **{("light_beam", system, view): lambda: one("light_beam", lambda m: [[(m.time, -BIG), (m.time, BIG)]])
       for system, view in (("null_cylindrical_interior", "edge"), ("null_cylindrical_exterior", "twice"),
                            ("null_cylindrical_exterior", "four"))},
    # Each moment is the plane of x and y at one event, on sigma = 0 at tau = t: u = v = sin(t/2).
    ("khan_penrose", "double_null", "plane"): lambda: [
        Mark(m, points=[(math.sin(m.time / 2), math.sin(m.time / 2))]) for m in moments("khan_penrose")],
    ("khan_penrose", "cosmological", "plane"): lambda: [Mark(m, points=[(m.time, 0.0)]) for m in moments("khan_penrose")],
    # Each moment of Bell and Szekeres's ring is the plane of x and y at one event, on eta = 0 at
    # xi = t with a = b = 1: u = v = t/2; rho = 0 and chi = t - pi/2; U = V = -cos t/(sqrt 2 (1 + sin t));
    # and t = tan xi, r = sec xi on y = 0.
    ("bell_szekeres", "double_null", "plane"): lambda: [
        Mark(m, points=[(m.time / 2, m.time / 2)]) for m in moments("bell_szekeres")],
    ("bell_szekeres", "time_space", "plane"): lambda: [Mark(m, points=[(m.time, 0.0)]) for m in moments("bell_szekeres")],
    ("bell_szekeres", "global", "plane"): lambda: [
        Mark(m, points=[(m.time - math.pi / 2, 0.0)]) for m in moments("bell_szekeres")],
    ("bell_szekeres", "kruskal_szekeres", "plane"): lambda: [
        Mark(m, points=[(-math.cos(m.time) / (math.sqrt(2) * (1 + math.sin(m.time))),) * 2])
        for m in moments("bell_szekeres")],
    ("bell_szekeres", "bertotti_robinson", "plane"): lambda: [
        Mark(m, points=[(math.tan(m.time), 1 / math.cos(m.time))]) for m in moments("bell_szekeres")],
    # The moments act = T of Senovilla's universe, out to where the embedding reaches.
    ("senovilla", "cylindrical", "radial"): lambda: one(
        "senovilla", lambda m: along(m.time, *m.reach("cylindrical", "\\rho"))),
    # Melvin's plane z = 0 at t = 0 and Ernst's equator at t = 0, each on its own chart's plane.
    ("melvin", "cylindrical", "radial"): lambda: one(
        "melvin", lambda m: along(0.0, *m.reach("cylindrical", "\\rho")), view_id="universe"),
    ("melvin", "ernst", "radial"): lambda: one("melvin", lambda m: along(0.0, *m.reach("ernst", "r")), view_id="ernst"),
    # Levi-Civita's plane z = 0 at t = 0, read in Weyl's coordinates; in the Kasner form r is the
    # proper distance, r = rho^Sigma/Sigma with Sigma = 3/4 at sigma = 1/4.
    # The particle's cone, with Gott and Alpert's planet at its apex or without, is read in the conical
    # chart out from the apex, whose r is the wedge chart's, R/alpha of the circumference radius and
    # (ell/alpha)(rho/ell)^alpha of the isotropic radius, at alpha = 3/4 and ell = 1; the planet's own
    # plane carries the planet alone.
    **{("point_particle_2plus1", system, "radial"): (lambda of=of: one(
        "point_particle_2plus1", lambda m: along(0.0, *(of(r) for r in m.reach("conical", "r", reference=True)))))
       for system, of in (("conical", lambda r: r), ("wedge", lambda r: r), ("circumference", lambda r: 0.75 * r),
                          ("isotropic", lambda r: (0.75 * r) ** (4 / 3)))},
    ("point_particle_2plus1", "planet", "radial"): lambda: one(
        "point_particle_2plus1", lambda m: along(0.0, *m.reach("planet", "\\chi")), view_id="cone"),
    # Lewis's moment t = 0 outside the null circle, embedded in the canonical chart: its radial line
    # on that chart's plane, and the whole line t = 0 of each cylinder outside r = ell, in the
    # canonical chart and in Lewis's, whose time is the same and whose angle turns rigidly.
    ("lewis", "canonical", "radial"): lambda: one("lewis", lambda m: along(0.0, *m.reach("canonical", "r"))),
    **{("lewis", system, view): lambda r=r: one("lewis", _lewis_cylinder(r))
       for system, view, r in (("canonical", "outside", 1.5), ("lewis", "between", 2.0), ("lewis", "beyond", 6.0))},
    ("levi_civita", "weyl", "radial"): lambda: one("levi_civita", lambda m: along(0.0, *m.reach("weyl", "\\rho"))),
    ("levi_civita", "kasner", "radial"): lambda: one(
        "levi_civita", lambda m: along(0.0, *(levi_civita_r(x) for x in m.reach("weyl", "\\rho")))),
    # The plane z = 0 midway between the two holes at t = 0, and the one event where it meets the axis.
    ("double_kerr", "weyl", "midplane"): lambda: one("double_kerr", lambda m: along(0.0, *m.reach("weyl", "\\rho"))),
    ("double_kerr", "weyl", "axis"): _double_kerr_axis,
    # The plane of Neugebauer and Meinel's disc at t = 0: the line t = 0 as far as the embedding reaches,
    # in the turning frame as far as its chart does, the plane beyond the rim in the spheroidal xi, where
    # rho = sqrt(1 + xi^2), and the one event where the plane meets the axis, the centre of the disc.
    ("neugebauer_meinel", "bardeen_wagoner", "plane"): lambda: one(
        "neugebauer_meinel", lambda m: along(0.0, *m.reach("bardeen_wagoner", "\\rho"))),
    ("neugebauer_meinel", "corotating", "plane"): lambda: one(
        "neugebauer_meinel", lambda m: along(0.0, *m.reach("bardeen_wagoner", "\\rho"))),
    ("neugebauer_meinel", "spheroidal", "plane"): lambda: one(
        "neugebauer_meinel", lambda m: along(0.0, 0.0, math.sqrt(m.reach("bardeen_wagoner", "\\rho")[1] ** 2 - 1))),
    ("neugebauer_meinel", "weyl", "axis"): lambda: _nm_centre("$t = 0$, $z = 0$"),
    ("neugebauer_meinel", "spheroidal", "axis"): lambda: _nm_centre("$t = 0$, $\\xi = 0$"),
    # The plane z = 0 of the first Morgan-Morgan disc at t = 0: Weyl's rho from the axis out, the
    # oblate spheroidal chart's eta across the disc and its xi = sqrt(rho^2/a^2 - 1) outside the rim.
    ("morgan_morgan", "weyl", "plane"): lambda: one("morgan_morgan", lambda m: along(0.0, *m.reach("weyl", "\\rho"))),
    ("morgan_morgan", "weyl", "axis"): _morgan_morgan_centre,
    ("morgan_morgan", "oblate_spheroidal", "axis"): _morgan_morgan_centre,
    ("morgan_morgan", "oblate_spheroidal", "disc"): lambda: one("morgan_morgan", lambda m: along(0.0, 0.0, 1.0)),
    ("morgan_morgan", "oblate_spheroidal", "plane"): lambda: one(
        "morgan_morgan", lambda m: along(0.0, 0.0, math.sqrt(m.reach("weyl", "\\rho")[1] ** 2 - 1))),
    # The equatorial plane of Bonnor's dipole at t = 0, and the one event where it meets the axis between the holes.
    ("bonnor_magnetic_dipole", "spheroidal", "equator"): lambda: one(
        "bonnor_magnetic_dipole", lambda m: along(0.0, *m.reach("spheroidal", "r"))),
    ("bonnor_magnetic_dipole", "spheroidal", "strut"): _bonnor_dipole_strut,
    # Bonnor's dust cloud: the moment t = 0 of the plane z = 0 outside the null circle, along the
    # radius of each equatorial view and the whole line t = 0 of the circle at 3a/2.
    **{("bonnor_rotating_dust", system, "equator"): lambda: one(
        "bonnor_rotating_dust", lambda m: along(0.0, *m.reach("cylindrical", "\\rho")))
       for system in ("cylindrical", "spherical")},
    ("bonnor_rotating_dust", "cylindrical", "outside"): lambda: one("bonnor_rotating_dust", _bonnor_circle(1.5)),
    # The plane z = 0 at t = 0, where the spherical chart's r is Weyl's rho.
    ("curzon_chazy", "weyl", "equator"): lambda: one("curzon_chazy", lambda m: along(0.0, *m.reach("weyl", "\\rho"))),
    ("curzon_chazy", "spherical", "equator"): lambda: one("curzon_chazy", lambda m: along(0.0, *m.reach("weyl", "\\rho"))),
    ("mcvittie", "isotropic", "radial"): lambda: one("mcvittie", lambda m: along(m.time, *m.reach("isotropic", "r"))),
    ("mcvittie", "areal", "radial"): lambda: one("mcvittie", _mcvittie_areal),
    # Sultana and Dyer's moments of the conformal time eta: a line of the Kerr-Schild plane, the curve
    # ct = eta - r_s ln(r/r_s - 1) of Schwarzschild's time, which runs off the drawing toward the horizon,
    # and the line v = eta + r of the advanced time.
    ("sultana_dyer", "kerr_schild", "radial"): lambda: one(
        "sultana_dyer", lambda m: along(m.time, *m.reach("kerr_schild", "r"))),
    ("sultana_dyer", "schwarzschild_time", "radial"): lambda: one("sultana_dyer", _sultana_dyer_t),
    ("sultana_dyer", "eddington_finkelstein_ingoing", "chart"): lambda: one(
        "sultana_dyer", lambda m: [[(m.time + r, r) for r in m.reach("kerr_schild", "r")]]),
    # The equatorial plane at t = 0 for each deformation, where the prolate spheroidal x is r/m - 1.
    **{("zipoy_voorhees", "spherical", f"equator_{shape}"): lambda shape=shape: one(
        "zipoy_voorhees", lambda m: along(0.0, *m.reach("spherical", "r")), view_id=shape)
       for shape in ("oblate", "prolate")},
    **{("zipoy_voorhees", "prolate_spheroidal", f"equator_{shape}"): lambda shape=shape: one(
        "zipoy_voorhees", lambda m: along(0.0, *(r - 1 for r in m.reach("spherical", "r"))), view_id=shape)
       for shape in ("oblate", "prolate")},
    ("tolman_vii", "spherical", "radial"): lambda: one("tolman_vii", lambda m: along(0.0, *m.reach("spherical", "r"))),
    ("tolman_vii", "spherical", "through"): lambda: one("tolman_vii", lambda m: along(0.0, *m.reach("spherical", "r"))),
    ("tolman_vii", "tolman", "radial"): lambda: one("tolman_vii", lambda m: along(0.0, *m.reach("spherical", "r"))),
    ("tolman_vii", "tolman", "through"): lambda: one("tolman_vii", lambda m: along(0.0, *m.reach("spherical", "r"))),
    ("boson_star", "areal", "radial"): lambda: one("boson_star", lambda m: along(0.0, *m.reach("areal", "r"))),
    ("boson_star", "areal", "through"): lambda: one("boson_star", lambda m: along(0.0, *m.reach("areal", "r"))),
    ("boson_star", "isotropic", "radial"): lambda: one("boson_star", _boson_star_isotropic),
    ("tov", "spherical", "radial"): lambda: one("tov", lambda m: along(0.0, *m.reach("spherical", "r"))),
    ("tov", "spherical", "through"): lambda: one("tov", lambda m: along(0.0, *m.reach("spherical", "r"))),
    # Bartnik and McKinnon's soliton with one zero at t = 0, from its centre to the edge of its embedding
    # diagram, in the areal radius and in the three radial coordinates that are functions of it.
    ("bartnik_mckinnon", "areal", "radial"): lambda: one("bartnik_mckinnon", lambda m: along(0.0, *m.reach("areal", "r")), view_id="n1"),
    ("bartnik_mckinnon", "areal", "through"): lambda: one("bartnik_mckinnon", lambda m: along(0.0, *m.reach("areal", "r")), view_id="n1"),
    ("bartnik_mckinnon", "isotropic", "radial"): lambda: one(
        "bartnik_mckinnon", lambda m: along(0.0, 0.0, _bm_radial("isotropic", m.reach("areal", "r")[1])), view_id="n1"),
    ("bartnik_mckinnon", "tortoise", "radial"): lambda: one(
        "bartnik_mckinnon", lambda m: along(0.0, 0.0, _bm_radial("tortoise", m.reach("areal", "r")[1])), view_id="n1"),
    ("bartnik_mckinnon", "flow", "radial"): lambda: one(
        "bartnik_mckinnon", lambda m: along(0.0, -50.0, math.log(_bm_radial("isotropic", m.reach("areal", "r")[1]))), view_id="n1"),
    ("einstein_cluster", "areal", "radial"): lambda: one("einstein_cluster", lambda m: along(0.0, *m.reach("areal", "r")), view_id="core"),
    ("einstein_cluster", "areal", "through"): lambda: one("einstein_cluster", lambda m: along(0.0, *m.reach("areal", "r")), view_id="core"),
    ("einstein_cluster", "constant_speed", "radial"): lambda: one("einstein_cluster", lambda m: along(0.0, *m.reach("constant_speed", "r")), view_id="speed"),
    # The isotropic radius and the polar angle of the same moments, from the centre to the surface.
    ("einstein_cluster", "isotropic", "radial"): lambda: one("einstein_cluster", lambda m: along(0.0, 0.0, 3 / (6 - 2 * math.sqrt(6)) ** 2), view_id="speed"),
    ("einstein_cluster", "uniform", "radial"): lambda: one("einstein_cluster", lambda m: along(0.0, *m.reach("uniform", "r")), view_id="uniform"),
    ("einstein_cluster", "uniform", "through"): lambda: one("einstein_cluster", lambda m: along(0.0, *m.reach("uniform", "r")), view_id="uniform"),
    ("einstein_cluster", "hyperspherical", "radial"): lambda: one("einstein_cluster", lambda m: along(0.0, 0.0, math.asin(1 / math.sqrt(3))), view_id="uniform"),
    ("malament_hogarth", "cartesian", "tx"): lambda: one("malament_hogarth", lambda m: across(m.time, *m.reach("cartesian", "x"))),
    ("oppenheimer_snyder", "interior_comoving", "through"): _os_interior,
    ("oppenheimer_snyder", "exterior_schwarzschild", "radial"): os_exterior,
    # v - r = w, every r the embedding reaches.
    ("c_metric", "spherical", "inner"): lambda: _c_metric(False),
    ("c_metric", "spherical", "outer"): lambda: _c_metric(False),
    ("c_metric", "hong_teo", "inner"): lambda: _c_metric(True),
    # The tail falling into the charged hole: v - r = w, every r the embedding reaches.
    ("mass_inflation", "ingoing", "tail"): lambda: one(
        "mass_inflation", lambda m: [[(m.time + r, r) for r in m.reach("ingoing", "r")]]),
    ("vaidya", "eddington_finkelstein_ingoing", "shell"): lambda: one(
        "vaidya", lambda m: [[(m.time + r, r) for r in m.reach("eddington_finkelstein_ingoing", "r")]]),
    # The same slices of Bonnor and Vaidya's charged shell.
    ("bonnor_vaidya", "eddington_finkelstein_ingoing", "shell"): lambda: one(
        "bonnor_vaidya", lambda m: [[(m.time + r, r) for r in m.reach("eddington_finkelstein_ingoing", "r")]]),
    # Wahlquist's fluid at the moment t = 0 of its own rest frame: on the equatorial plane the line
    # t = 0 from the ring to the surface of zero pressure, on the disc from the ring to the axis, and
    # on Whittaker's sphere from the centre to its surface.
    ("wahlquist", "wahlquist", "equator"): lambda: one(
        "wahlquist", lambda m: along(0.0, *m.reach("wahlquist", "\\xi")), view_id="rotating"),
    ("wahlquist", "wahlquist", "disc"): lambda: one(
        "wahlquist", lambda m: along(0.0, *m.reach("wahlquist", "\\eta")), view_id="rotating"),
    ("wahlquist", "whittaker", "radial"): lambda: one(
        "wahlquist", lambda m: along(0.0, *m.reach("whittaker", "X")), view_id="static"),
    ("wahlquist", "whittaker", "through"): lambda: one(
        "wahlquist", lambda m: along(0.0, *m.reach("whittaker", "X")), view_id="static"),
}
FLAT_METRICS = {key[0] for key in FLAT}

# Where a moment of the spacetime lies on the drawing and is not drawn, and why.
HIDDEN = {
    ("kundt_waves", "kundt", "front"): "a wave with no cosmological constant, another spacetime than the waves in de Sitter and anti-de Sitter space whose fronts are embedded",
    ("kundt_waves", "podolsky_belan", "near"): "a wave with no cosmological constant, another spacetime than the waves in de Sitter and anti-de Sitter space whose fronts are embedded",
    ("kundt_waves", "podolsky_belan", "far"): "a wave with no cosmological constant, another spacetime than the waves in de Sitter and anti-de Sitter space whose fronts are embedded",
    ("kundt_waves", "simplest_wave", "front"): "a wave with no cosmological constant, another spacetime than the waves in de Sitter and anti-de Sitter space whose fronts are embedded",
    ("kundt_waves", "simplest_wave", "depth"): "a wave with no cosmological constant, another spacetime than the waves in de Sitter and anti-de Sitter space whose fronts are embedded",
    ("kundt_waves", "kerr_schild", "fronts"): "a wave with no cosmological constant, another spacetime than the waves in de Sitter and anti-de Sitter space whose fronts are embedded",
    ("siklos", "kaigorodov_stationary", "plane"): "the region x < 0 of Siklos's chart, another region than the one whose wave front is embedded",
    ("btz", "stationary", "rotating"): "the rotating hole, J = 4l/5, another spacetime than the hole without rotation whose moment is embedded",
    ("btz", "eddington_finkelstein_ingoing", "rotating"): "the rotating hole, J = 4l/5, another spacetime than the hole without rotation whose moment is embedded",
    ("btz", "eddington_finkelstein_outgoing", "rotating"): "the rotating hole, J = 4l/5, another spacetime than the hole without rotation whose moment is embedded",
    ("btz", "rotating"): "the rotating hole, J = 4l/5, another spacetime than the hole without rotation whose moment is embedded",
    ("ads_soliton", "five_dimensional", "radial"): "the soliton of five dimensions, another spacetime than the soliton of four whose moment is embedded",
    ("ads_soliton", "three_dimensional", "radial"): "the soliton of three dimensions, which is anti-de Sitter space, another spacetime than the soliton of four whose moment is embedded",
    **{("topological_black_hole", system, "negative"): "the hyperbolic hole of negative mass, another spacetime than the flat hole and the hyperbolic hole without mass whose moments are embedded"
       for system in ("static", "eddington_finkelstein_ingoing", "eddington_finkelstein_outgoing", "hyperbolic")},
    **{("topological_black_hole", f"{system}_negative"): "the hyperbolic hole of negative mass, another spacetime than the flat hole and the hyperbolic hole without mass whose moments are embedded"
       for system in ("static", "ingoing", "outgoing", "hyperbolic")},
    **{("topological_black_hole", f"eddington_finkelstein_{way}", "massless"): "the bifurcation surface of the hyperbolic hole without mass lies off both Eddington-Finkelstein charts, and the string's moment is the flat hole's"
       for way in ("ingoing", "outgoing")},
    **{("coleman_de_luccia", system, view): "a bubble with another vacuum inside than the anti-de Sitter space whose moments are embedded"
       for system, view in (("wall", "into_flat"), ("open", "zero"), ("open", "positive"), ("static_inside", "into_flat"),
                            ("static_outside", "into_flat"))},
    ("coleman_de_luccia", "into_flat"): "a bubble with another vacuum inside than the anti-de Sitter space whose moments are embedded",
    **{("coleman_de_luccia", system, "out_of_flat"): "the region outside the light cone of the bubble's centre, which no moment of the open universe inside meets"
       for system in ("wall", "static_outside")},
    **{("kerr_taub_nut", system, "axis"): "the regular half of the axis, which the embedded equatorial plane does not meet"
       for system in ("one_string", "kerr_ingoing", "kerr_outgoing")},
    **{("kerr_taub_nut", view): "the regular half of the axis, which the embedded equatorial plane does not meet"
       for view in ("axis", "ingoing", "outgoing")},
    ("kerr_taub_nut", "plebanski", "principal"): "the moment of constant t is a surface on which tau changes with sigma, the coordinate the drawing leaves out",
    ("frw", "comoving_spherical", "radial"): "the flat universe, k = 0, whose moments are planes; the moments embedded are the closed universe's",
    ("frw", "comoving_spherical", "through"): "the flat universe, k = 0, whose moments are planes; the moments embedded are the closed universe's",
    ("frw", "conformal_spherical", "radial"): "the flat universe, k = 0, whose moments are planes; the moments embedded are the closed universe's",
    ("godel", "cylindrical", "beyond"): "beyond r_c the circles are closed timelike curves and no surface of constant t is a moment of space; the embedding stops at sinh^2 r = 1/sqrt 2",
    ("stockum_dust", "cylindrical", "beyond"): "beyond r = R the circles are closed timelike curves; the embedding stops at r = 0.83 R",
    ("som_raychaudhuri", "cylindrical", "beyond"): "beyond r_c the circles are closed timelike curves; the embedding stops at r = sqrt(3) r_c/2",
    **{("lewis", system, "inside"): "inside r = ell the circles are closed timelike curves; the embedding begins at ell"
       for system in ("lewis", "canonical")},
    **{("lewis", system, view): "another member of Lewis's family than the cylinder of the Weyl class whose moment is embedded"
       for system, views in (("lewis_class", ("first", "second")), ("stockum_light", ("surface", "beyond")),
                             ("stockum_critical", ("surface", "beyond")), ("stockum_heavy", ("surface", "band")))
       for view in views},
    **{("spinning_string", system, "inside"): "inside r_c the circles are closed timelike curves; the embedding begins at r_c"
       for system in ("proper_radius", "rescaled_radius", "helical")},
    **{("point_particle_2plus1", "two_bodies", view): "two particles at rest, another spacetime than the one particle whose cone is embedded"
       for view in ("between", "beyond")},
    ("point_particle_2plus1", "two_bodies"): "two particles at rest, another spacetime than the one particle whose cone is embedded",
    ("point_particle_2plus1", "moving", "wedge"): "the particle in motion, whose moment of the frame's t is not the moment of its rest frame that is embedded",
    ("tolman_bondi", "comoving_synchronous", "collapse"): "the marginally bound cloud, E = 0, whose moments are planes; the cloud embedded is released from rest",
    ("vaidya", "eddington_finkelstein_outgoing", "shell"): "the exploding shell, the time reverse of the imploding shell embedded",
    ("bonnor_vaidya", "eddington_finkelstein_outgoing", "shell"): "the leaving shell, the time reverse of the falling shell embedded",
    ("bonnor_vaidya", "leaving"): "the leaving shell, the time reverse of the falling shell embedded",
    ("kaluza_klein_black_hole", "dyonic", "radial"): "the hole of equal charges, another member of the family than the holes of one charge embedded",
    ("kaluza_klein_black_hole", "equal"): "the hole of equal charges, another member of the family than the holes of one charge embedded",
    ("bonnor_vaidya", "homothetic", "scaling"): "the collapse of a mass and a charge that grow with the advanced time, another spacetime than the shell embedded",
    ("mass_inflation", "ingoing", "behind"): "behind Ori's shell, where the mass function is another one than the tail's, whose moments are embedded",
    ("mass_inflation", "shell"): "Ori's shell and the region behind it, another spacetime than the tail falling in alone, whose moments are embedded",
    ("hiscock", "ingoing", "shells"): "the simplest model, a hole made and removed by two shells, another spacetime than the one embedded",
    **{("bonnor_rotating_dust", system, "axis"): "the axis, which the embedded plane z = 0 meets only at the centre, inside where the embedding begins"
       for system in ("cylindrical", "spherical")},
    ("bonnor_rotating_dust", "cylindrical", "inside"): "inside rho = a the circles are closed timelike curves; the embedding begins at rho = a",
    **{("curzon_chazy", system, "axis"): "the axis, which the embedded plane z = 0 meets only at rho = 0, inside where the embedding stops"
       for system in ("weyl", "spherical")},
    **{("bonnor_rotating_dust", f"{system}_axis"): "the axis, which the embedded plane z = 0 meets only at the centre, inside where the embedding begins"
       for system in ("cylindrical", "spherical")},
    **{("curzon_chazy", f"{system}_axis"): "the axis, which the embedded plane z = 0 meets only at rho = 0, inside where the embedding stops"
       for system in ("weyl", "spherical")},
    ("bonnor_magnetic_dipole", "spheroidal", "axis"): "the axis beyond a hole, which the embedded equatorial plane does not meet",
    ("bonnor_magnetic_dipole", "spheroidal_axis"): "the axis beyond a hole, which the embedded equatorial plane does not meet",
    ("double_kerr", "weyl_axis_outside"): "the axis above the upper hole, which the embedded plane z = 0 does not meet",
    **{("neugebauer_meinel", "black_hole_limit", view): "the limit mu -> mu_0, the extreme Kerr metric; the moment embedded is the disc's at mu = 3"
       for view in ("axis", "equator")},
    ("neugebauer_meinel", "limit_axis"): "the limit mu -> mu_0, the extreme Kerr metric; the moment embedded is the disc's at mu = 3",
    **{("zipoy_voorhees", system, f"axis_{shape}"): "the axis, which the embedded equatorial plane does not meet"
       for system in ("spherical", "prolate_spheroidal") for shape in ("oblate", "prolate")},
    **{("zipoy_voorhees", f"{system}_axis_{shape}"): "the axis, which the embedded equatorial plane does not meet"
       for system in ("spherical", "prolate_spheroidal") for shape in ("oblate", "prolate")},
    **{("szekeres", "axisymmetric", half): "the axis of symmetry, which the embedded surface through the equators of the shells meets only at the centre r = 0"
       for half in ("north", "south")},
    **{("photon_rocket", system, half): "the axis of flight, which the embedded surface of the rays that leave the rocket sideways meets only at r = 0"
       for system in ("rectilinear", "robinson_trautman") for half in ("behind", "ahead")},
    **{("photon_rocket", half): "the axis of flight, which the embedded surface of the rays that leave the rocket sideways meets only at r = 0"
       for half in ("behind", "ahead")},
    ("frw", "flat"): "the flat universe's conformal diagram; the moments embedded are the closed universe's",
    ("misner", "rindler", "plane"): "the region T > 0 beyond the chronology horizon, which no moment of the contracting region meets",
    ("misner", "rindler"): "the region T > 0 beyond the chronology horizon, which no moment of the contracting region meets",
    ("gott_time_machine", "grant_rindler", "plane"): "the region of closed timelike curves beyond the chronology horizon, which no moment of Grant's Milne time meets",
    ("gott_time_machine", "grant_rindler"): "the region of closed timelike curves beyond the chronology horizon, which no moment of Grant's Milne time meets",
    ("gott_time_machine", "centre_of_momentum", "loop"): "the centre of momentum chart about the strings; the moments embedded are Grant's, away from the strings",
    ("bell_szekeres", "regular", "plane"): "a plane of constant X and Y with X^2 + Y^2 < 1, off eta = 0, where the embedded ring's centre lies on the rim X^2 + Y^2 = 1 of the regular chart",
    ("bell_szekeres", "regular"): "a plane of constant X and Y with X^2 + Y^2 < 1, off eta = 0, where the embedded ring's centre lies on the rim X^2 + Y^2 = 1 of the regular chart",
    ("gowdy", "sphere"): "the inside of Schwarzschild's horizon, a universe on S^2 x S^1; the moments embedded are the torus universe's",
    ("gowdy", "sphere", "plane"): "the inside of Schwarzschild's horizon, a universe on S^2 x S^1; the moments embedded are the torus universe's",
    **{("light_beam", "null_cartesian", view): "two beams side by side, another spacetime than the single beam whose wave fronts are embedded"
       for view in ("midway", "one")},
    ("light_beam", "midway"): "two beams side by side, another spacetime than the single beam whose wave fronts are embedded",
    ("light_beam", "cartesian", "lens"): "the plane y = 0 with t left out, which every wave front covers whole",
    ("lifshitz_spacetime", "poincare", "rays"): "the plane y = 0 with t left out, which every moment of the static spacetime covers whole",
    ("schrodinger_spacetime", "global", "trap"): "the surface X = 0 with V left out, where a line of constant T holds every V and the embedded plane only V = 0",
    ("wormhole_time_machine", "wormhole", "speeding"): "the axis of the acceleration, theta = 0, which the embedded plane theta = pi/2 meets nowhere",
    ("wormhole_time_machine", "wormhole", "slowing"): "the axis of the acceleration, theta = 0, which the embedded plane theta = pi/2 meets nowhere",
    ("wormhole_time_machine", "short_throat", "radial"): "the mouth of the throat of zero length; the moment embedded is the smooth wormhole's",
    ("wormhole_time_machine", "lorentz", "trip"): "the flat space outside the mouths, with the mouths drawn as world lines; the moment embedded runs through the throat",
    ("wormhole_time_machine", "lorentz"): "the flat space outside the mouths, with the mouths drawn as world lines; the moment embedded runs through the throat",
    ("frw", "open"): "the open universe's conformal diagram; the moments embedded are the closed universe's",
    ("myers_perry", "boyer_lindquist_six", "rotation"): "the plane of rotation in six dimensions, theta = pi/2, which the embedded transverse plane theta = 0 meets nowhere outside the horizon",
    ("black_saturn", "ring", "outside"): "the ring alone, with no hole inside it; the moment embedded is the Saturn's",
    ("black_saturn", "ring", "inside"): "the ring alone, with no hole inside it; the moment embedded is the Saturn's",
    ("near_horizon_extreme_kerr", "near_nhek", "equator"): "the patch ct > r_0^2/r of the Poincare chart, to the future of the ray that leaves the boundary at t = 0, which the moment tau = 0 embedded does not enter",
    ("near_horizon_extreme_kerr", "near_nhek"): "the patch ct > r_0^2/r of the Poincare chart, to the future of the ray that leaves the boundary at t = 0, which the moment tau = 0 embedded does not enter",
    **{("wahlquist", system, "equator"): "the plane with Mars's angle divided out, whose circles each run through every moment of Wahlquist's t"
       for system in ("mars", "mars_ingoing")},
    **{("hartle_thorne", system, "axis"): "the axis of rotation, which the embedded equatorial plane does not meet"
       for system in ("hartle_thorne", "painleve_gullstrand")},
    ("hartle_thorne", "painleve_gullstrand", "equator"): "the Painleve-Gullstrand line element, which agrees with Hartle and Thorne's to first order in the spin and no further; the moment embedded is one of Hartle and Thorne's t",
    ("hartle_thorne", "painleve_gullstrand"): "the Painleve-Gullstrand line element, which agrees with Hartle and Thorne's to first order in the spin and no further; the moment embedded is one of Hartle and Thorne's t",
}


# ---------------------------------------------------------------- geometry the drawings share

def clip_runs(P, lo=(0.0, 0.0), hi=(1.0, 1.0)):
    """A polyline cut to the box [lo, hi], as the runs of it inside, each cut on the box's edge."""
    P = np.asarray(P, dtype=float)
    lo, hi = np.asarray(lo, dtype=float), np.asarray(hi, dtype=float)
    runs, run = [], []
    for a, b in zip(P[:-1], P[1:]):
        seg = clip_segment(a, b, lo, hi)
        if seg is None:
            if len(run) > 1:
                runs.append(np.array(run))
            run = []
            continue
        if run and np.allclose(run[-1], seg[0], atol=1e-12, rtol=0):
            run.append(seg[1])
        else:
            if len(run) > 1:
                runs.append(np.array(run))
            run = [seg[0], seg[1]]
        if not np.allclose(seg[1], b, atol=1e-12, rtol=0):
            runs.append(np.array(run))
            run = []
    if len(run) > 1:
        runs.append(np.array(run))
    return [r for r in runs if np.linalg.norm(r[-1] - r[0]) > 1e-9 or len(r) > 2]


def clip_segment(a, b, lo, hi):
    """The part of the segment from a to b inside the box, or None."""
    s0, s1 = 0.0, 1.0
    d = b - a
    for k in range(2):
        if abs(d[k]) < 1e-15:
            if not lo[k] - 1e-12 <= a[k] <= hi[k] + 1e-12:
                return None
            continue
        t0, t1 = sorted(((lo[k] - a[k]) / d[k], (hi[k] - a[k]) / d[k]))
        s0, s1 = max(s0, t0), min(s1, t1)
    if s0 > s1 or (s0 == s1 and np.linalg.norm(d) > 0):
        return None
    A, B = a + s0 * d, a + s1 * d
    return np.array([np.clip(A, lo, hi), np.clip(B, lo, hi)])


def clip_ring(P, lo=(0.0, 0.0), hi=(1.0, 1.0)):
    """A closed polygon cut to the box, by Sutherland and Hodgman."""
    out = [np.asarray(p, dtype=float) for p in P]
    for axis, bound, below in ((0, hi[0], True), (0, lo[0], False), (1, hi[1], True), (1, lo[1], False)):
        def inside(p):
            return p[axis] <= bound if below else p[axis] >= bound
        new = []
        for i in range(len(out)):
            a, b = out[i - 1], out[i]
            if inside(b):
                if not inside(a):
                    new.append(a + (b - a) * (bound - a[axis]) / (b[axis] - a[axis]))
                new.append(b)
            elif inside(a):
                new.append(a + (b - a) * (bound - a[axis]) / (b[axis] - a[axis]))
        out = new
        if not out:
            return None
    return np.array(out)


# ---------------------------------------------------------------- where a slice's label stands

# A slice is drawn 2.6 wide, so its edge is 1.3 from its line, and its point is a dot of radius
# 4.5 with an outline 1.2 wide, so its edge is 5.1 from its centre, both in the units the page
# draws the drawing in: a spacetime diagram's plot 520 wide, a conformal diagram 628.
EDGE, DOT = 1.3, 5.1


# The parts of a label's mathematics that MathJax sets wider than a letter: a relation, with
# the room it leaves on either side, and a function's name, set upright letter by letter.
RELATION = r"=|<|>|\\to(?![a-zA-Z])|\\leq?(?![a-zA-Z])|\\geq?(?![a-zA-Z])|\\approx|\\neq?(?![a-zA-Z])|\\sim(?![a-zA-Z])|\\equiv|\\in(?![a-zA-Z])|\\rightarrow|\\mapsto"
FUNCTION = r"\\(sinh|cosh|tanh|sin|cos|tan|ln|log|exp|arctan|min|max)(?![a-zA-Z])"
# The letters MathJax sets about an em wide, where every other is near half of one.
WIDE = "mMW"


def label_size(text):
    """A label's box as the page sets it, in ems of its own size, no smaller than MathJax sets
    it: in mathematics 0.7 em a letter, 1.2 an m, M or W, 1.1 a relation with its room, 1.0 any
    other command, 0.5 a letter of a function's name and 0.3 a space in the source; prose at the 0.6 em of
    Source Code Pro; the padding of its ground, 0.3 em either side; and 1.5 em tall, or 1.75
    with a root or a superscript. Measured in Chrome on 30 September 2026 against every label
    of every drawing, 632 of them set at the page's caption size, it is never smaller, and on
    average 1.5 em wider: "$r \\to \\infty$" 3.9 em wide, estimated 4.0, and "$R/\\sqrt{2}$"
    1.66 em tall, estimated 1.75. An m is 0.96 em wide there, an M 1.07 and a W 1.04, where a
    digit is 0.55 and an r 0.5, so a label that is little but an m was wider than 0.7 em a letter
    made it until 1 October 2026: "$4m$" 2.11 em wide, estimated 2.0, and now 2.5.
    _layouts/mfs.html carries the same function as cdLabelSize()."""
    width, tall = 0.6, False
    for part in re.split(r"(\$[^$]*\$)", text):
        if part.startswith("$"):
            math = re.sub(r"\\[,;:!]", " ", part[1:-1])
            tall = tall or "^" in math or "\\sqrt" in math
            width += 0.5 * sum(len(f) for f in re.findall(FUNCTION, math))
            math = re.sub(FUNCTION, "", math)
            width += 1.1 * len(re.findall(RELATION, math))
            math = re.sub(RELATION, "", math)
            math = re.sub(r"\\(bar|hat|tilde|vec|dot|mathrm|text|left|right|mathscr|mathcal|mathfrak|operatorname)(?![a-zA-Z])",
                          "", math)
            width += 1.0 * len(re.findall(r"\\[a-zA-Z]+", math))
            math = re.sub(r"[{}^_]", "", re.sub(r"\\[a-zA-Z]+", "", math))
            width += 0.7 * len(math.replace(" ", "")) + 0.5 * sum(math.count(c) for c in WIDE) + 0.3 * math.count(" ")
        else:
            width += 0.6 * len(part)
    return width, 1.75 if tall else 1.5


def place(marks, px, box, size, others=()):
    """Where each mark's label stands, on the edge of its line at the line's right hand end or
    on the edge of its point, so that no label overlaps another label of the drawing, covers
    a line of a mark or leaves the box. Labels may touch, and nothing is added between them.
    A mark that is a region alone names its moment in the legend and gets none.

    marks   the Marks, their lines and points already in the drawing's own coordinates;
    px      the map from those coordinates to the units the page draws the drawing in, y down;
    box     (width, height) of the drawing in those units;
    size    the size of a label in those units;
    others  the boxes (x0, y0, x1, y1) of the drawing's other labels.

    The sides are tried in order: at the right hand end of the lines on the line's upper edge
    and over the rest of it, on its lower edge, past the end above the line and past it below;
    then the same at the left hand end; then the four corners of the point 85% of the way
    along from the left hand end; and for a point from above and to its right round to below
    and to its left, each against the dot's edge. The side that stays in the box, overlaps no
    label, covers none of its own line and none of another mark's wins, in that order of what
    matters most, the earliest of equals. Returns, for each mark, None or
    {"at": [x, y], "anchor": ..., "dx": ..., "dy": ...}: the label's corner `anchor` stands at
    `at` moved by dx and dy in the drawing's units, to the edge of the line or the dot."""
    W, H = box

    def traced(mark):
        out = []
        for line in mark.lines:
            P = np.array([px(q) for q in line])
            for a, b in zip(P[:-1], P[1:]):
                n = max(1, int(np.ceil(np.linalg.norm(b - a) / 2)))
                out += list(a + (b - a) * np.linspace(0, 1, n + 1)[:, None])
        return np.array(out) if out else np.zeros((0, 2))
    traces = [traced(mark) for mark in marks]
    taken = [tuple(o) for o in others]

    def covers(rect, trace):
        return bool(np.any((trace[:, 0] > rect[0]) & (trace[:, 0] < rect[2])
                           & (trace[:, 1] > rect[1]) & (trace[:, 1] < rect[3])))
    out = []
    for i, mark in enumerate(marks):
        if not (mark.lines or mark.points):
            out.append(None)
            continue
        w, h = (v * size for v in label_size(mark.label))
        corners = (("bl", 1, -1), ("tl", 1, 1), ("br", -1, -1), ("tr", -1, 1))
        if mark.lines:
            ends = [q for line in mark.lines for q in (line[0], line[-1])]
            right = max(ends, key=lambda q: (px(q)[0], -px(q)[1]))
            left = min(ends, key=lambda q: (px(q)[0], px(q)[1]))
            line = max(mark.lines, key=lambda L: px(L[-1])[0] if px(L[-1])[0] >= px(L[0])[0] else px(L[0])[0])
            P = np.array(line if px(line[-1])[0] >= px(line[0])[0] else line[::-1], dtype=float)
            arc = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(np.array([px(q) for q in P]), axis=0), axis=1))])
            along = np.array([np.interp(0.85 * arc[-1], arc, P[:, j]) for j in range(2)])
            e = EDGE
            tries = ([(right, "br", 0, -e), (right, "tr", 0, e), (right, "bl", e, -e), (right, "tl", e, e),
                      (left, "bl", 0, -e), (left, "tl", 0, e), (left, "br", -e, -e), (left, "tr", -e, e)]
                     + [(along, a, sx * e, sy * e) for a, sx, sy in corners])
        else:
            tries = [(mark.points[0], a, sx * DOT, sy * DOT) for a, sx, sy in corners]
        best = None
        for order, (at, anchor, dx, dy) in enumerate(tries):
            x, y = px(at)
            x0 = x + dx - (w if anchor[1] == "r" else 0)
            y0 = y + dy - (h if anchor[0] == "b" else 0)
            rect = (x0, y0, x0 + w, y0 + h)
            score = (not (rect[0] >= 0 and rect[1] >= 0 and rect[2] <= W and rect[3] <= H),
                     any(rect[0] < o[2] and o[0] < rect[2] and rect[1] < o[3] and o[1] < rect[3] for o in taken),
                     covers(rect, traces[i]),
                     any(covers(rect, t) for j, t in enumerate(traces) if j != i), order)
            if best is None or score < best[0]:
                best = (score, at, anchor, dx, dy, rect)
        _, at, anchor, dx, dy, rect = best
        taken.append(rect)
        out.append({"at": [round(float(v), 4) for v in at], "anchor": anchor, "dx": dx, "dy": dy})
    return out


# ---------------------------------------------------------------- the checks

def checks():
    """Every transformation a moment is carried through, checked by pulling one published
    metric back onto the other at points of the drawn plane, and the closed forms against
    what they stand for. Returns the failures."""
    import sympy as sp
    import null_rays as nr

    failures, rng = [], np.random.default_rng(3)

    def report(name, miss, tol):
        ok = miss < tol
        print(f"  {'ok  ' if ok else 'FAIL'} {name}: {miss:.1e}")
        if not ok:
            failures.append(name)

    def metric(metric_id, system, params):
        """The published metric of a system, c = 1 and its parameters set, with its coordinates."""
        _, entry, reader = nr.load(metric_id, system)
        g = nr.published_matrix(reader, entry, "metric_components")
        subs = {reader.c: 1}
        subs.update({reader.parameters[k]: nr.number(v) for k, v in params.items()})
        return g.subs(subs), [reader.symbol[c] for c in entry["coords"]]

    # Van Den Broeck: with rho = vdb_radius(l), the comoving chart's B^2(d rho^2 + rho^2 d phi^2)
    # is the proper distance chart's dl^2 + r(l)^2 d phi^2 for the declared pocket.
    g_c, (_, rho, *_) = metric("van_den_broeck", "pocket", {"R": 1})
    g_l, (_, lv, *_) = metric("van_den_broeck", "proper_radial", {"l_0": 4})
    _, _, reader_c = nr.load("van_den_broeck", "pocket")
    _, _, reader_l = nr.load("van_den_broeck", "proper_radial")
    B = sp.lambdify(rho, nr._as_lambda(reader_c, "B", nr._vdb_factor("r"))(rho), "numpy")
    areal = sp.lambdify(lv, nr._as_lambda(reader_l, "r", nr._VDB_R)(lv), "numpy")
    ls = rng.uniform(-3.99, 0.99, 400)
    radii = np.array([vdb_radius(v) for v in ls])
    step = 1e-6
    slope = np.array([(vdb_radius(v + step) - vdb_radius(v - step)) / (2 * step) for v in ls])
    with np.errstate(all="ignore"):
        factor = B(radii)
    miss = max(float(np.max(np.abs(factor * radii - areal(ls)))), float(np.max(np.abs(factor * slope - 1))),
               float(np.max(np.abs(nr._vdb_distance(radii) - ls))))
    report("Van Den Broeck: rho(l) pulls the comoving chart of the pocket back onto the proper distance one",
           miss, 1e-6)

    # Bertotti-Robinson: x = b^2/r and the same t carry the static plane onto the Poincare one.
    g_s, (t, r, *_) = metric("bertotti_robinson", "static", {"b": 1})
    g_p, (tp, x, *_) = metric("bertotti_robinson", "poincare", {"b": 1})
    J = sp.Matrix([[1, 0], [0, sp.diff(1 / r, r)]])
    pulled = J.T * g_p[:2, :2].subs({tp: t, x: 1 / r}) * J
    miss = max(abs(float((pulled - g_s[:2, :2]).subs(r, rv)[i, j])) for rv in rng.uniform(0.3, 3, 20)
               for i in range(2) for j in range(2))
    report("Bertotti-Robinson: x = b^2/r pulls the Poincare plane back onto the static one", miss, 1e-12)

    # Kiselev: the tortoise coordinate's slope is the published g_rr of the chart of the linear
    # term, it vanishes at the centre, and the chart that keeps w is that chart at w = -2/3.
    g_k, (tk, rk, *_) = metric("kiselev", "linear", {"r_s": 1, "r_q": 8})
    g_w, (tw, rw, *_) = metric("kiselev", "static", {"r_s": 1, "r_q": 8, "w": "-2/3"})
    grr = sp.lambdify(rk, g_k[1, 1], "numpy")
    radii = np.concatenate([rng.uniform(0.05, 1.1, 20), rng.uniform(1.25, 6.7, 20), rng.uniform(6.95, 9, 20)])
    step = 1e-6
    miss = max(float(np.max(np.abs((kiselev_rstar(radii + step) - kiselev_rstar(radii - step)) / (2 * step)
                                   / grr(radii) - 1))), abs(float(kiselev_rstar(0.0))),
               max(abs(float(sp.simplify(g_w[i, i].subs({tw: tk, rw: rk}) - g_k[i, i]).subs(rk, 3))) for i in range(2)))
    report("Kiselev: dr_*/dr is the published g_rr, r_* = 0 at the centre, and w = -2/3 is the linear term", miss, 1e-6)
    # Kiselev's matter alone: ct = 2 r_q eta and r = r_q(1 - e^(-2 chi)) carry the static plane onto the
    # hyperbolic one, and eta = ln(tau^2 - rho^2)/2, chi = artanh(rho/tau) the hyperbolic plane onto the
    # conformally flat one, where eta = 0 is tau^2 - rho^2 = 1.
    g_f, (tf, rf, *_) = metric("kiselev", "linear", {"r_s": 0, "r_q": 1})
    g_h, (eta, chi, *_) = metric("kiselev", "hyperbolic", {"r_q": 1})
    g_c, (tau, rho, *_) = metric("kiselev", "conformally_flat", {"r_q": 1})
    image = [2 * eta, 1 - sp.exp(-2 * chi)]
    J = sp.Matrix(2, 2, lambda i, j: sp.diff(image[i], [eta, chi][j]))
    pulled = J.T * g_f[:2, :2].subs({tf: image[0], rf: image[1]}, simultaneous=True) * J
    miss = max(abs(float((pulled - g_h[:2, :2]).subs(chi, v)[i, j])) for v in rng.uniform(0.1, 3, 20)
               for i in range(2) for j in range(2))
    image = [sp.log(tau ** 2 - rho ** 2) / 2, sp.atanh(rho / tau)]
    J = sp.Matrix(2, 2, lambda i, j: sp.diff(image[i], [tau, rho][j]))
    pulled = J.T * g_h[:2, :2].subs({eta: image[0], chi: image[1]}, simultaneous=True) * J
    miss = max(miss, max(abs(float((pulled - g_c[:2, :2]).subs({tau: a + b, rho: b})[i, j]))
                         for a, b in zip(rng.uniform(0.1, 3, 20), rng.uniform(0.1, 3, 20))
                         for i in range(2) for j in range(2)))
    report("Kiselev: his map to the hyperbolic chart and Fock's transformation pull the static plane of the matter alone back onto "
           "the hyperbolic and the conformally flat ones", miss, 1e-10)

    # Bardeen: the tortoise coordinate's slope is the published g_rr of the static chart, it
    # vanishes at the centre, and the horizons named here are the zeros of the published g^rr.
    g_b, (tb, rb, *_) = metric("bardeen", "static", {"r_s": 1, "g": "1/3"})
    grr = sp.lambdify(rb, g_b[1, 1], "numpy")
    pts = np.concatenate([rng.uniform(0.02, 0.28, 8), rng.uniform(0.33, 0.74, 8), rng.uniform(0.82, 6, 8)])
    h = 1e-5
    miss = max(abs((bardeen_rstar(v + h) - bardeen_rstar(v - h)) / (2 * h) / grr(v) - 1) for v in pts)
    report("Bardeen: dr_*/dr is the published g_rr", float(miss), 1e-7)
    report("Bardeen: r_* = 0 at the centre", abs(float(bardeen_rstar(0.0))), 1e-15)
    report("Bardeen: the horizons are the zeros of the published g^rr",
           max(abs(1 / grr(ri)) for ri in BARDEEN_HORIZONS), 1e-13)

    # de Sitter: the static chart of the flat slicing, r = rho e^t_f and
    # t_s = t_f - ln(1 - rho^2 e^(2 t_f))/2, pulls the static plane back onto the flat one, and the
    # moment t_s = 0 is t_f = -ln(1 + rho^2)/2.
    g_s, (ts, rs, th, _) = metric("de_sitter", "static_spherical", {"Lambda": 3})
    g_f, (tf, xf, *_) = metric("de_sitter", "flat_slicing", {"H": 1})
    T_s = tf - sp.log(1 - xf ** 2 * sp.exp(2 * tf)) / 2
    R_s = xf * sp.exp(tf)
    J = sp.Matrix([[sp.diff(T_s, tf), sp.diff(T_s, xf)], [sp.diff(R_s, tf), sp.diff(R_s, xf)]])
    pulled = J.T * g_s[:2, :2].subs({rs: R_s, th: sp.pi / 2}, simultaneous=True) * J
    pts = [(a, b) for a, b in zip(rng.uniform(-0.5, 0.2, 20), rng.uniform(0.1, 0.8, 20)) if b * math.exp(a) < 0.95]
    miss = max(abs(float((pulled - g_f[:2, :2]).subs({tf: a, xf: b})[i, j])) for a, b in pts for i in range(2) for j in range(2))
    report("de Sitter: the static chart pulls back onto the flat slicing's plane y = z = 0", miss, 1e-12)
    xs = np.linspace(-3, 3, 13)
    miss = max(abs(float(T_s.subs({tf: -0.5 * math.log1p(a * a), xf: a}))) for a in xs)
    report("de Sitter: t_f = -ln(1 + x^2)/2 is the static t = 0", miss, 1e-14)

    # Godel: the published transformation carries the Cartesian metric onto the cylindrical one
    # on the plane y = 0, where phi = 0 or pi, t_x = 2t and x = +-2r.
    g_c, (tx_, xx, yy, zz) = metric("godel", "cartesian", {"omega": 1})
    g_y, (ty, ry, py, zy) = metric("godel", "cylindrical", {"omega": 1})
    X = sp.log(sp.cosh(2 * ry) + sp.cos(py) * sp.sinh(2 * ry))
    Y = sp.sqrt(2) * sp.sin(py) * sp.sinh(2 * ry) / sp.exp(X)
    Tx = 2 * ty + sp.sqrt(2) * (2 * sp.atan(sp.exp(-2 * ry) * sp.tan(py / 2)) - py)
    new = [Tx, X, Y, 2 * zy]
    J = sp.Matrix([[sp.diff(f, v) for v in (ty, ry, py, zy)] for f in new])
    pulled = J.T * g_c.subs({tx_: Tx, xx: X, yy: Y, zz: 2 * zy}, simultaneous=True) * J
    miss = 0.0
    for a, b, c in zip(rng.uniform(-1, 1, 12), rng.uniform(0.05, 0.7, 12), rng.uniform(-2.5, 2.5, 12)):
        diff = (pulled - g_y).subs({ty: a, ry: b, py: c, zy: 0.3})
        miss = max(miss, max(abs(float(diff[i, j])) for i in range(4) for j in range(4)))
    report("Godel: the published transformation pulls the Cartesian metric onto the cylindrical", miss, 1e-10)
    miss = max(abs(float(Tx.subs({ty: 0, ry: b, py: c}))) + abs(float(X.subs({ry: b, py: c}) - (2 * b if c == 0 else -2 * b)))
               for b in (0.1, 0.4, 0.76) for c in (0, sp.pi - 1e-12))
    report("Godel: t = 0 on the plane y = 0 is t_x = 0 at x = +-2r", miss, 1e-9)

    # Anti-de Sitter: the Poincare plane y = 0, z = L at t = 0 has static radius
    # r^2 = x^2 + x^4/4L^2, read from the embedding space.
    x = np.linspace(-2.5, 2.5, 11)
    X1, X3 = x, x * x / 2
    report("anti-de Sitter: r^2 = x^2 + x^4/4 on the Poincare plane y = 0, z = L",
           float(np.max(np.abs(X1 ** 2 + X3 ** 2 - (x ** 2 + x ** 4 / 4)))), 1e-14)
    hi = moments("anti_de_sitter")[0].reach("static_global", "r")[1]
    xm = math.sqrt(2 * (math.sqrt(1 + hi * hi) - 1))
    report("anti-de Sitter: the reach r = 4L is |x| = sqrt(2(sqrt 17 - 1)) L", abs(xm ** 2 + xm ** 4 / 4 - hi * hi), 1e-12)

    # Kantowski-Sachs: ct = eta + sin(eta) cos(eta), a = 1 + eta tan(eta) and b = cos^2(eta) carry the
    # comoving chart onto the dust chart at b_0 = 1 and kappa = 0, so the moment eta_k is one t.
    _, entry, reader = nr.load("kantowski_sachs", "comoving")
    g_c = nr.published_matrix(reader, entry, "metric_components").subs(reader.c, 1)
    g_d, (ed, *_) = metric("kantowski_sachs", "dust", {"b_0": 1, "kappa": 0})
    lapse = sp.diff(ed + sp.sin(ed) * sp.cos(ed), ed)
    g_c = g_c.subs({reader.parameters["a"]: 1 + ed * sp.tan(ed), reader.parameters["b"]: sp.cos(ed) ** 2})
    pulled = sp.diag(lapse, 1, 1, 1) * g_c * sp.diag(lapse, 1, 1, 1)
    theta = reader.symbol["\\theta"]
    miss = max(abs(float((pulled - g_d).subs({ed: e, theta: 1.1})[i, i])) for e in rng.uniform(-1.4, 1.4, 20) for i in range(4))
    report("Kantowski-Sachs: the comoving chart with the dust's a, b and t pulls back onto the dust chart", miss, 1e-10)

    # Misner: T = -t^2/4 and psi = 2 chi - ln(t^2/4) carry the Milne plane of t and chi onto
    # Misner's plane of T and psi, both being the covering Minkowski plane's t - x = t e^(-chi) =
    # -2 e^(-psi/2) and t + x = t e^chi = 2T e^(psi/2), so the moment t = t_k is T = -t_k^2/4.
    g_m, (Tm, pm, *_) = metric("misner", "misner", {"psi_0": "4*pi"})
    g_n, (tn, cn, *_) = metric("misner", "milne", {"psi_0": "4*pi"})
    new = [-tn ** 2 / 4, 2 * cn - sp.log(tn ** 2 / 4)]
    J = sp.Matrix([[sp.diff(f, v) for v in (tn, cn)] for f in new])
    pulled = J.T * g_m[:2, :2].subs({Tm: new[0], pm: new[1]}, simultaneous=True) * J
    miss = max(abs(float((pulled - g_n[:2, :2]).subs({tn: a, cn: b})[i, j]))
               for a, b in zip(rng.uniform(-3, -0.1, 20), rng.uniform(-3, 3, 20)) for i in range(2) for j in range(2))
    report("Misner: T = -t^2/4, psi = 2 chi - ln(t^2/4) pulls Misner's plane back onto the Milne plane", miss, 1e-12)

    # Ori: T = t + a(x^2 - y^2)/2 - e(x^2 + y^2) carries the foliation onto the vacuum core at
    # f = a(x^2 - y^2)/2, so the moment t = t_k is T = t_k on the central circle and T = t_k - 3/2 at
    # x = 4, y = 0; and T = -uv/2, z = -2 ln(-u/2) carries the Brinkmann chart onto it, so the same
    # moment is uv = -2 t_k there.
    ori = dict(nr.ORI)
    g_c, (Tc, xc, yc, zc) = metric("ori_time_machine", "vacuum_core", {"L": ori["L"]})
    _, _, reader = nr.load("ori_time_machine", "vacuum_core")
    a_, e_ = sp.Rational(ori["a"]), sp.Rational(ori["e"])
    g_c = g_c.subs(reader.parameters["f"], a_ * (xc ** 2 - yc ** 2) / 2)
    g_f, (tf, xf, yf, zf) = metric("ori_time_machine", "foliation", ori)
    new = [tf + a_ * (xf ** 2 - yf ** 2) / 2 - e_ * (xf ** 2 + yf ** 2), xf, yf, zf]
    J = sp.Matrix([[sp.diff(f, v) for v in (tf, xf, yf, zf)] for f in new])
    pulled = J.T * g_c.subs(dict(zip((Tc, xc, yc, zc), new)), simultaneous=True) * J
    miss = max(abs(float((pulled - g_f).subs({tf: p[0], xf: p[1], yf: p[2], zf: p[3]})[i, j]))
               for p in rng.uniform(-3, 3, (20, 4)) for i in range(4) for j in range(4))
    report("Ori: T = t + a(x^2 - y^2)/2 - e(x^2 + y^2) pulls the vacuum core back onto the foliation", miss, 1e-12)
    report("Ori: at x = 4, y = 0 the moment t is T = t - 3/2",
           abs(float(new[0].subs({tf: 0, xf: 4, yf: 0})) + 1.5), 1e-12)
    g_b, (ub, vb, xb, yb) = metric("ori_time_machine", "brinkmann", {"a": ori["a"], "L": ori["L"]})
    new = [-ub * vb / 2, xb, yb, -2 * sp.log(-ub / 2)]
    J = sp.Matrix([[sp.diff(f, v) for v in (ub, vb, xb, yb)] for f in new])
    pulled = J.T * g_c.subs(dict(zip((Tc, xc, yc, zc), new)), simultaneous=True) * J
    miss = max(abs(float((pulled - g_b).subs({ub: -abs(p[0]) - 0.1, vb: p[1], xb: p[2], yb: p[3]})[i, j]))
               for p in rng.uniform(-3, 3, (20, 4)) for i in range(4) for j in range(4))
    report("Ori: T = -uv/2, z = -2 ln(-u/2) pulls the vacuum core back onto the Brinkmann chart", miss, 1e-11)

    # Novikov: each shell at its proper time is a radial geodesic from rest of the exterior,
    # which keeps its energy E = sqrt(1 - r_s/R) = (1 - r_s/r) dt/dtau, and the same shell in
    # Kruskal's U and V is the same event as its t and r.
    for R in (2.0, 3.0, 4.0):
        taus = np.linspace(0.05, 0.95, 9) * 0.5 * R * math.sqrt(R) * math.pi
        E, h, miss = math.sqrt(1 - 1 / R), 1e-6, 0.0
        for tau in taus:
            r, _, _ = novikov(np.array([R]), tau)
            if r[0] > 1.02:
                _, ta, _ = novikov(np.array([R]), tau + h)
                _, tb, _ = novikov(np.array([R]), tau - h)
                miss = max(miss, abs((ta[0] - tb[0]) / (2 * h) * (1 - 1 / r[0]) - E) / E)
        report(f"Novikov: the shell from R = {R:g} r_s keeps its energy along Schwarzschild's t", miss, 1e-6)
        r, t, _ = novikov(np.full(len(taus), R), 0.0)
        U, V, r = np.array([novikov_kruskal(np.array([R]), tau) for tau in taus])[:, :, 0].T
        _, t, _ = np.array([novikov(np.array([R]), tau) for tau in taus])[:, :, 0].T
        out = r > 1.02
        v = t[out] + r[out] + np.log(r[out] - 1)
        report(f"Novikov: V = exp(v/2) along the shell from R = {R:g} r_s outside r_s",
               float(np.max(np.abs(np.log(V[out]) - v / 2))), 1e-10)
        report(f"Novikov: U V = (1 - r) e^r along the shell from R = {R:g} r_s",
               float(np.max(np.abs(U * V - (1 - r) * np.exp(r)))), 1e-12)

    # Nariai: r = -cosh(tau) cos(chi) and sinh(t) = sinh(tau)/sqrt(1 - r^2) carry the global plane
    # into the static patch 0 < chi < pi, sin(chi) > |tanh(tau)|, pulling the static metric back onto
    # the global one, and nariai_static_t is the same map at fixed tau.
    g_s, (ts, rs, *_) = metric("nariai", "static", {"Lambda": 1})
    g_g, (tg, cg, *_) = metric("nariai", "global", {"Lambda": 1})
    R_s = -sp.cosh(tg) * sp.cos(cg)
    T_s = sp.asinh(sp.sinh(tg) / sp.sqrt(1 - R_s ** 2))
    J = sp.Matrix([[sp.diff(T_s, tg), sp.diff(T_s, cg)], [sp.diff(R_s, tg), sp.diff(R_s, cg)]])
    pulled = J.T * g_s[:2, :2].subs({rs: R_s}, simultaneous=True) * J
    pts = [(a, b) for a, b in zip(rng.uniform(-1.5, 1.5, 40), rng.uniform(0.05, math.pi - 0.05, 40))
           if math.sin(b) > abs(math.tanh(a)) + 0.05]
    miss = max(abs(float((pulled - g_g[:2, :2]).subs({tg: a, cg: b})[i, j])) for a, b in pts for i in range(2) for j in range(2))
    report("Nariai: the static chart pulls back onto the global plane over the patch 0 < chi < pi", miss, 1e-12)
    miss = max(abs(float(T_s.subs({tg: a, cg: b})) - float(nariai_static_t(a, float(R_s.subs({tg: a, cg: b})))))
               for a, b in pts)
    report("Nariai: nariai_static_t is the global moment in the static chart", miss, 1e-12)

    # McVittie: R = ar(1 + 1/(4ar))^2 and the same t carry the isotropic plane onto the areal one,
    # with H = (da/dt)/a and sqrt(1 - 1/R) = (4ar - 1)/(4ar + 1) outside the throat, at the scale
    # factor the diagrams declare; mcvittie_areal is that R.
    g_i, (ti, ri, *_) = metric("mcvittie", "isotropic", {"r_s": 1})
    g_a, (ta, Ra, *_) = metric("mcvittie", "areal", {"r_s": 1})
    _, _, reader_i = nr.load("mcvittie", "isotropic")
    _, _, reader_a = nr.load("mcvittie", "areal")
    scale = sp.sinh(sp.Rational(3, 2) * ti / sp.sqrt(15)) ** sp.Rational(2, 3)
    g_i = g_i.replace(reader_i.parameters["a"].func, sp.Lambda(ti, scale)).doit()
    hubble = sp.diff(scale, ti) / scale
    g_a = g_a.replace(reader_a.parameters["H"].func, sp.Lambda(ta, hubble.subs(ti, ta))).doit()
    R_of = scale * ri * (1 + 1 / (4 * scale * ri)) ** 2
    J = sp.Matrix([[1, 0], [sp.diff(R_of, ti), sp.diff(R_of, ri)]])
    pulled = J.T * g_a[:2, :2].subs({ta: ti, Ra: R_of}, simultaneous=True) * J
    pts = [(a, b) for a, b in zip(rng.uniform(0.5, 8, 40), rng.uniform(0.05, 4, 40))
           if 4 * math.sinh(1.5 * MCV_H0 * a) ** (2 / 3) * b > 1.05]
    miss = max(abs(complex((pulled - g_i[:2, :2]).subs({ti: a, ri: b})[i, j])) for a, b in pts for i in range(2) for j in range(2))
    report("McVittie: the areal chart pulls back onto the isotropic plane outside the throat", miss, 1e-10)
    miss = max(abs(float(R_of.subs({ti: a, ri: b})) - mcvittie_areal(a, b)) for a, b in pts)
    report("McVittie: mcvittie_areal is the areal radius of the comoving r", miss, 1e-12)
    # The lukewarm hole: r = H tau rho + r_s/2 and the static time of rnds_static_t pull the static
    # plane back onto the cosmological one, and r_* of rnds_rstar has dr_*/dr = 1/f.
    g_s, (ts, rs_, *_) = metric("reissner_nordstrom_de_sitter", "static", {"r_s": 1, "r_q": "1/2", "Lambda": "27/64"})
    g_c, (tc, xc, *_) = metric("reissner_nordstrom_de_sitter", "cosmological", {"r_s": 1, "H": "3/8"})
    f_s = sp.lambdify(rs_, -g_s[0, 0], "numpy")
    plane = sp.lambdify((tc, xc), g_c[:2, :2], "numpy")
    miss = 0.0
    for tau, rho in zip(rng.uniform(0.5, 4, 40), rng.uniform(0.2, 2.5, 40)):
        r = RNDS_H * tau * rho + 0.5
        if min(abs(r - a) for a in RNDS_ROOTS) < 0.05:
            continue
        h = 1e-6
        dT = [(rnds_static_t(tau + h, RNDS_H * (tau + h) * rho + 0.5) - rnds_static_t(tau - h, RNDS_H * (tau - h) * rho + 0.5)) / (2 * h),
              (rnds_static_t(tau, RNDS_H * tau * (rho + h) + 0.5) - rnds_static_t(tau, RNDS_H * tau * (rho - h) + 0.5)) / (2 * h)]
        J = np.array([dT, [RNDS_H * rho, RNDS_H * tau]])
        pulled = J.T @ np.diag([-f_s(r), 1 / f_s(r)]) @ J
        there = np.array(plane(tau, rho), dtype=float)
        miss = max(miss, float(np.max(np.abs(pulled - there) / (1 + np.abs(there)))))
    report("Reissner-Nordstrom-de Sitter: rnds_static_t pulls the static plane back onto the cosmological one", miss, 1e-6)
    rr = rng.uniform(0.05, 2.5, 40)
    rr = rr[np.min(np.abs(rr[:, None] - np.array(RNDS_ROOTS)[None, :]), axis=1) > 0.05]
    slope = (rnds_rstar(rr + 1e-6) - rnds_rstar(rr - 1e-6)) / 2e-6
    report("Reissner-Nordstrom-de Sitter: rnds_rstar has dr_*/dr = 1/f and vanishes at r = 0",
           float(np.max(np.abs(slope * f_s(rr) - 1))) + abs(float(rnds_rstar(0.0))), 1e-6)

    # Kastor and Traschen's one hole: R = H tau r + m and the static time of kt_static_t pull the
    # lukewarm hole's static plane, at r_s = 2m, back onto the isotropic chart's, and r_* of kt_rstar
    # has dr_*/dR = 1/f. Two holes: H tau = e^{Ht} pulls the Cartesian chart's time back onto the
    # comoving chart's, where g_tt = -(d tau/dt)^2/U^2 and g_xx = U^2 with a Omega for U.
    g_s, (ts, rs_, *_) = metric("reissner_nordstrom_de_sitter", "static", {"r_s": 2, "r_q": 1, "Lambda": "27/256"})
    g_k, (tk, rk, *_) = metric("kastor_traschen", "isotropic", {"m": 1, "H": "-3/16"})
    f_s = sp.lambdify(rs_, -g_s[0, 0], "numpy")
    plane = sp.lambdify((tk, rk), g_k[:2, :2], "numpy")
    miss = 0.0
    for tau, r in zip(rng.uniform(-4, 4, 60), rng.uniform(0.2, 8, 60)):
        R = KT_H_ONE * tau * r + 1
        if R < 0.05 or abs(tau) < 0.05 or min(abs(R - a) for a in (1.0,) + KT_ONE_ROOTS) < 0.05:
            continue
        h = 1e-6
        dT = [(kt_static_t(tau + h, r) - kt_static_t(tau - h, r)) / (2 * h),
              (kt_static_t(tau, r + h) - kt_static_t(tau, r - h)) / (2 * h)]
        J = np.array([dT, [KT_H_ONE * r, KT_H_ONE * tau]])
        pulled = J.T @ np.diag([-f_s(R), 1 / f_s(R)]) @ J
        there = np.array(plane(tau, r), dtype=float)
        miss = max(miss, float(np.max(np.abs(pulled - there) / (1 + np.abs(there)))))
    report("Kastor-Traschen: kt_static_t pulls the lukewarm hole's static plane back onto one hole's isotropic plane", miss, 1e-6)
    RR = rng.uniform(0.05, 5, 40)
    RR = RR[np.min(np.abs(RR[:, None] - np.array(KT_ONE_ROOTS)[None, :]), axis=1) > 0.05]
    slope = (kt_rstar(RR + 1e-6) - kt_rstar(RR - 1e-6)) / 2e-6
    report("Kastor-Traschen: kt_rstar has dr_*/dR = 1/f and vanishes at R = 0",
           float(np.max(np.abs(slope * f_s(RR) - 1))) + abs(float(kt_rstar(0.0))), 1e-6)
    _, entry_c, reader_c = nr.load("kastor_traschen", "cartesian")
    _, entry_m, reader_m = nr.load("kastor_traschen", "comoving")
    g_c = nr.published_matrix(reader_c, entry_c, "metric_components").subs(reader_c.held).doit()
    g_m = nr.published_matrix(reader_m, entry_m, "metric_components").subs(reader_m.held).doit()
    tau_c, t_m = reader_c.symbol["\\tau"], reader_m.symbol["t"]
    W = sp.Symbol("W", positive=True)
    at_c = {reader_c.c: 1, reader_c.parameters["H"]: KT_H_TWO, reader_c.parameters["V"]: W}
    at_m = {reader_m.c: 1, reader_m.parameters["H"]: KT_H_TWO, reader_m.parameters["V"]: W}
    image = sp.exp(KT_H_TWO * t_m) / KT_H_TWO
    pulled = [g_c[0, 0].subs(at_c).subs(tau_c, image) * sp.diff(image, t_m) ** 2, g_c[1, 1].subs(at_c).subs(tau_c, image)]
    own = [g_m[0, 0].subs(at_m), g_m[1, 1].subs(at_m)]
    miss = max(abs(float((a - b).subs({t_m: t, W: w}))) / (1 + abs(float(b.subs({t_m: t, W: w}))))
               for a, b in zip(pulled, own) for t, w in zip(rng.uniform(0, 24, 20), rng.uniform(0.1, 3, 20)))
    report("Kastor-Traschen: H tau = e^{Ht} pulls the Cartesian chart back onto the comoving one", miss, 1e-10)
    report("Kastor-Traschen: kt_comoving_t inverts H tau = e^{Ht}",
           float(np.max(np.abs(np.exp(KT_H_TWO * kt_comoving_t(np.array([-12.0, -3.0, -1.0]))) / KT_H_TWO
                               - np.array([-12.0, -3.0, -1.0])))), 1e-12)

    # The travelling wave on a string: x = X - A, y = Y - B and v = V - 2A'(X - A) - 2B'(Y - B) -
    # int_0^u (A'^2 + B'^2) carry the moving string chart onto the isotropic one, for the declared
    # pulse, and string_wave_V is the V of the surface v = 0 on a line of fixed X.
    pulse = {"A": "exp(-4*u**2)/2", "B": "0"}
    _, entry_i, R_i = nr.load("string_wave", "isotropic")
    _, entry_m, R_m = nr.load("string_wave", "moving_string")

    def declared(R, entry):
        g = nr.published_matrix(R, entry, "metric_components").applyfunc(R.surface)
        for name, rep in pulse.items():
            g = g.replace(R.parameters[name].func, nr._as_lambda(R, name, rep)).doit()
        return g.subs({R.parameters["b"]: sp.Rational(1, 2), R.parameters["ell"]: 1}), [R.symbol[c] for c in entry["coords"]]
    g_i, (ui, vi, xi, yi) = declared(R_i, entry_i)
    g_m, (um, Vm, Xm, Ym) = declared(R_m, entry_m)
    A = sp.exp(-4 * um ** 2) / 2
    s = sp.Symbol("s", real=True)
    image = [um, Vm - 2 * sp.diff(A, um) * (Xm - A) - sp.integrate(sp.diff(A, um).subs(um, s) ** 2, (s, 0, um)), Xm - A, Ym]
    J = sp.Matrix(4, 4, lambda i, j: sp.diff(image[i], (um, Vm, Xm, Ym)[j]))
    pulled = J.T * g_i.subs(dict(zip((ui, vi, xi, yi), image)), simultaneous=True) * J
    miss = 0.0
    for a, b, c, d in zip(rng.uniform(-2, 2, 30), rng.uniform(-2, 2, 30), rng.uniform(-2, 2, 30), rng.uniform(0.1, 2, 30)):
        at = {um: a, Vm: b, Xm: c, Ym: d}
        miss = max(miss, max(abs(complex((pulled - g_m).subs(at)[i, j])) for i in range(4) for j in range(4)))
    report("string wave: the isotropic chart pulls back onto the moving string chart", miss, 1e-10)
    miss = max(abs(float((Vm - image[1]).subs({um: a, Xm: X})) - string_wave_V(a, X))
               for a in (-1.5, -0.3, 0.0, 1.0, 2.5) for X in (-0.25, 0.75))
    report("string wave: string_wave_V is V on the surface v = 0", miss, 1e-12)

    # The throat of extreme Kerr: Bardeen and Horowitz's map from their global chart, r = r_0 (sqrt(1 +
    # y^2) cos tau + y), ct = r_0 sqrt(1 + y^2) sin tau/(sqrt(1 + y^2) cos tau + y) and phi moved by
    # ln((cos tau + y sin tau)/(1 + sqrt(1 + y^2) sin tau)), pulls the Poincare chart back onto the
    # global one in every slot, and on tau = 0 it is t = 0, r = r_0 nhek_radius(y) and the same phi.
    g_p, (tp, rp, thp, php) = metric("near_horizon_extreme_kerr", "poincare", {"r_0": 1})
    g_g, (tg, yg, thg, phg) = metric("near_horizon_extreme_kerr", "global", {"r_0": 1})
    g_x, (tx, xx, thx, phx) = metric("near_horizon_extreme_kerr", "inverse_radius", {"r_0": 1})
    root = sp.sqrt(1 + yg ** 2)
    image = [root * sp.sin(tg) / (root * sp.cos(tg) + yg), root * sp.cos(tg) + yg, thg,
             phg + sp.log((sp.cos(tg) + yg * sp.sin(tg)) / (1 + root * sp.sin(tg)))]
    J = sp.Matrix(4, 4, lambda i, j: sp.diff(image[i], (tg, yg, thg, phg)[j]))
    pulled = J.T * g_p.subs(dict(zip((tp, rp, thp, php), image)), simultaneous=True) * J
    miss = 0.0
    for a, b, c in zip(rng.uniform(-0.4, 0.4, 20), rng.uniform(-2, 2, 20), rng.uniform(0.2, 2.9, 20)):
        at = {tg: a, yg: b, thg: c}
        miss = max(miss, max(abs(complex((pulled - g_g).subs(at)[i, j])) for i in range(4) for j in range(4)))
    report("NHEK: Bardeen and Horowitz's map pulls the Poincare chart back onto the global one", miss, 1e-10)
    miss = max(max(abs(float(e.subs({tg: 0, yg: b, phg: 0.3}) - v)) for e, v in zip(image, (0.0, nhek_radius(b), None, 0.3))
                   if v is not None) for b in rng.uniform(-2.25, 2.25, 20))
    report("NHEK: the moment tau = 0 is t = 0 with r = r_0 (sqrt(1 + y^2) + y) and the same phi", miss, 1e-12)
    J = sp.diag(1, sp.diff(1 / xx, xx), 1, 1)
    pulled = J.T * g_p.subs({tp: tx, rp: 1 / xx, thp: thx, php: phx}, simultaneous=True) * J
    miss = max(abs(complex((pulled - g_x).subs({xx: a, thx: c})[i, j])) for a, c in zip(
        rng.uniform(0.2, 4, 20), rng.uniform(0.2, 2.9, 20)) for i in range(4) for j in range(4))
    report("NHEK: x = r_0^2/r pulls the Poincare chart back onto the inverse radius one", miss, 1e-10)
    return failures


if __name__ == "__main__":
    print("the transformations the slices are carried through:")
    sys.exit(1 if checks() else 0)
