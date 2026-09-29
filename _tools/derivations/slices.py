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
    """One surface of one embedding view, and the moment of the spacetime it is cut from."""

    def __init__(self, metric_id, view, index, label=None):
        self.metric, self.view, self.index = metric_id, view, index
        self.surface = view["surfaces"][index]
        self.time = self.surface.get("time")
        # A sequence names each moment on its surface; a single moment is drawn at t = 0 of
        # the chart the embedding read it in, save where the view names another.
        self.label = self.surface.get("label") or label or "$t = 0$"
        self.version = build.embedding_moment_version(view, index)

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
        return {"view": self.view["id"], "surface": self.index, "label": self.label, "version": self.version}


def moments(metric_id, view_id=None, label=None):
    """Every moment of a spacetime's embedding views, in the order of views and surfaces."""
    out = []
    for view in embedding(metric_id)["views"]:
        if view_id is None or view["id"] == view_id:
            out += [Moment(metric_id, view, i, label) for i in range(len(view["surfaces"]))]
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


def schwarzschild_t(sign):
    """Schwarzschild's t = 0 in an Eddington-Finkelstein chart, r_s = 1: v = r + ln(r - 1)
    in the ingoing chart and u = -r - ln(r - 1) in the outgoing one, outside r_s."""
    m = moments("schwarzschild")[0]
    lo, hi = m.reach("spherical", "r")
    r = near(lo, hi)
    return [Mark(m, [np.column_stack([sign * (r + np.log(r - 1)), r])])]


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


def one(metric_id, lines_of, label=None, view_id=None):
    """Each moment of a spacetime as the lines lines_of(moment) returns."""
    return [Mark(m, lines_of(m), label=label) for m in moments(metric_id, view_id)]


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


def _ds_flat():
    """de Sitter's static moment t = 0 in the flat slicing, H = 1: with
    X_0 = sinh t_f + rho^2 e^(t_f)/2 of the embedding space it is X_0 = 0, which is
    t_f = -ln(1 + rho^2)/2, and rho = |x| on the plane y = z = 0. The whole of the observer's
    hemisphere lies on it, its r = rho e^(t_f) reaching 1 only as rho -> infinity."""
    m = moments("de_sitter", label="static $t = 0$")[0]
    x = np.linspace(-BIG ** 0.25, BIG ** 0.25, 4 * N + 1)
    return [Mark(m, [np.column_stack([-0.5 * np.log1p(x * x), x])])]


def _ads_poincare():
    """Anti-de Sitter's static moment t = 0 is the Poincare moment t = 0; on the plane y = 0,
    z = L its static radius is r^2 = x^2 + x^4/4L^2, so the embedding's reach r <= 4L is
    |x| <= sqrt(2(sqrt(17) - 1)) L."""
    m = moments("anti_de_sitter")[0]
    hi = m.reach("static_global", "r")[1]
    x = math.sqrt(2 * (math.sqrt(1 + hi * hi) - 1))
    return [Mark(m, across(0.0, 0.0, x))]


FLAT = {
    ("schwarzschild", "spherical", "radial"): lambda: one("schwarzschild", lambda m: along(0.0, *m.reach("spherical", "r"))),
    ("schwarzschild", "eddington_finkelstein_ingoing", "finkelstein"): lambda: schwarzschild_t(1),
    ("schwarzschild", "eddington_finkelstein_ingoing", "chart"): lambda: schwarzschild_t(1),
    ("schwarzschild", "eddington_finkelstein_outgoing", "finkelstein"): lambda: schwarzschild_t(-1),
    ("schwarzschild", "eddington_finkelstein_outgoing", "chart"): lambda: schwarzschild_t(-1),
    ("ellis_bronnikov", "spherical", "radial"): lambda: one("ellis_bronnikov", lambda m: along(0.0, *m.reach("spherical", "r"))),
    ("morris_thorne", "spherical", "radial"): lambda: one("morris_thorne", lambda m: along(0.0, *m.reach("spherical", "r"))),
    ("minkowski", "spherical", "radial"): lambda: one("minkowski", lambda m: along(0.0, *m.reach("spherical", "r"))),
    # t = (u + v)/2 and r = (v - u)/2, so the moment is u = -r, v = r.
    ("minkowski", "spherical_null", "radial"): lambda: one(
        "minkowski", lambda m: [[(-r, r) for r in m.reach("spherical", "r")]]),
    ("minkowski", "cartesian", "tx"): lambda: one("minkowski", lambda m: across(0.0, 0.0, m.reach("spherical", "r")[1])),
    # ct = X sinh(aT/c) vanishes in the wedge only at T = 0, where x = X.
    ("minkowski", "rindler", "tx"): lambda: one("minkowski", lambda m: along(0.0, 0.0, m.reach("spherical", "r")[1])),
    ("de_sitter", "static_spherical", "radial"): lambda: one("de_sitter", lambda m: along(0.0, *m.reach("static_spherical", "r"))),
    ("de_sitter", "static_spherical", "through"): lambda: one("de_sitter", lambda m: along(0.0, *m.reach("static_spherical", "r"))),
    ("de_sitter", "flat_slicing", "tx"): _ds_flat,
    ("anti_de_sitter", "static_global", "radial"): lambda: one("anti_de_sitter", lambda m: along(0.0, *m.reach("static_global", "r"))),
    ("anti_de_sitter", "static_global", "through"): lambda: one("anti_de_sitter", lambda m: along(0.0, *m.reach("static_global", "r"))),
    ("anti_de_sitter", "poincare", "tx"): _ads_poincare,
    ("rn_metric", "spherical", "radial"): lambda: _rn("outside") + _rn("inside"),
    ("taub_nut", "spherical", "radial"): lambda: one("taub_nut", lambda m: along(0.0, *m.reach("spherical", "r"))),
    ("bertotti_robinson", "static", "radial"): lambda: _br("static"),
    ("bertotti_robinson", "poincare", "tx"): lambda: _br("poincare"),
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
    ("kasner", "cartesian", "tx"): lambda: one("kasner", lambda m: across(m.time, 0.0, BIG)),
    ("kasner", "cartesian", "tz"): lambda: one("kasner", lambda m: across(m.time, 0.0, BIG)),
    ("bianchi", "type_i_cartesian", "tx"): lambda: one("bianchi", lambda m: across(m.time, 0.0, BIG)),
    # Godel's t_x = 2t + sqrt(2)(2 arctan(e^(-2r) tan(phi/2)) - phi) and e^x = cosh 2r + cos(phi) sinh 2r
    # put the plane y = 0 at phi = 0 and pi, where t_x = 2t and x = +-2r.
    ("godel", "cartesian", "tx"): lambda: one("godel", lambda m: across(0.0, 0.0, 2 * m.reach("cylindrical", "r")[1])),
    ("godel", "cylindrical", "inside"): lambda: one("godel", _godel_cylinder(math.asinh(1.0) / 2)),
    ("stockum_dust", "cylindrical", "inside"): lambda: one("stockum_dust", _godel_cylinder(0.5)),
    ("alcubierre", "cartesian", "tx"): lambda: one("alcubierre", lambda m: across(0.0, 0.0, m.grid()["u"][-1])),
    ("natario", "cartesian_flow", "tx"): lambda: one("natario", lambda m: across(0.0, 0.0, m.reach("cartesian_flow", "y")[1])),
    ("krasnikov", "cylindrical", "tx"): _krasnikov,
    # The wave front u = u_k, every v.
    ("pp_wave", "exact_plane_wave", "tz"): lambda: one("pp_wave", lambda m: [[(m.time, -BIG), (m.time, BIG)]]),
    ("tov", "spherical", "radial"): lambda: one("tov", lambda m: along(0.0, *m.reach("spherical", "r"))),
    ("tov", "spherical", "through"): lambda: one("tov", lambda m: along(0.0, *m.reach("spherical", "r"))),
    ("malament_hogarth", "cartesian", "tx"): lambda: one("malament_hogarth", lambda m: across(m.time, *m.reach("cartesian", "x"))),
    ("oppenheimer_snyder", "interior_comoving", "through"): _os_interior,
    ("oppenheimer_snyder", "exterior_schwarzschild", "radial"): os_exterior,
    # v - r = w, every r the embedding reaches.
    ("vaidya", "eddington_finkelstein_ingoing", "shell"): lambda: one(
        "vaidya", lambda m: [[(m.time + r, r) for r in m.reach("eddington_finkelstein_ingoing", "r")]]),
}
FLAT_METRICS = {key[0] for key in FLAT}

# Where a moment of the spacetime lies on the drawing and is not drawn, and why.
HIDDEN = {
    ("frw", "comoving_spherical", "radial"): "the flat universe, k = 0, whose moments are planes; the moments embedded are the closed universe's",
    ("frw", "comoving_spherical", "through"): "the flat universe, k = 0, whose moments are planes; the moments embedded are the closed universe's",
    ("frw", "conformal_spherical", "radial"): "the flat universe, k = 0, whose moments are planes; the moments embedded are the closed universe's",
    ("godel", "cylindrical", "beyond"): "beyond r_c the circles are closed timelike curves and no surface of constant t is a moment of space; the embedding stops at sinh^2 r = 1/sqrt 2",
    ("stockum_dust", "cylindrical", "beyond"): "beyond r = R the circles are closed timelike curves; the embedding stops at r = 0.83 R",
    ("tolman_bondi", "comoving_synchronous", "collapse"): "the marginally bound cloud, E = 0, whose moments are planes; the cloud embedded is released from rest",
    ("vaidya", "eddington_finkelstein_outgoing", "shell"): "the exploding shell, the time reverse of the imploding shell embedded",
    ("frw", "flat"): "the flat universe's conformal diagram; the moments embedded are the closed universe's",
    ("frw", "open"): "the open universe's conformal diagram; the moments embedded are the closed universe's",
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

    # Bertotti-Robinson: x = b^2/r and the same t carry the static plane onto the Poincare one.
    g_s, (t, r, *_) = metric("bertotti_robinson", "static", {"b": 1})
    g_p, (tp, x, *_) = metric("bertotti_robinson", "poincare", {"b": 1})
    J = sp.Matrix([[1, 0], [0, sp.diff(1 / r, r)]])
    pulled = J.T * g_p[:2, :2].subs({tp: t, x: 1 / r}) * J
    miss = max(abs(float((pulled - g_s[:2, :2]).subs(r, rv)[i, j])) for rv in rng.uniform(0.3, 3, 20)
               for i in range(2) for j in range(2))
    report("Bertotti-Robinson: x = b^2/r pulls the Poincare plane back onto the static one", miss, 1e-12)

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
    return failures


if __name__ == "__main__":
    print("the transformations the slices are carried through:")
    sys.exit(1 if checks() else 0)
