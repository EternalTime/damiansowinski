"""The drawings of Schein and Aichelburg's wormhole between two charged shells, at R = 1, the unit of
every length: the spheres of radius R are centred at z = +-d_1 with d_1 = 2R, and each has the mass
parameter m_1 = R/2, so V_0 = 1 + m_1/R = 3/2 on both and the shells sit at the areal radius
R + m_1 = 3R/2 of the black hole between them, whose mass and charge are m = R/4 and q = 6R/25,
q/2m = 0.48 as the live Reissner-Nordstrom tower is drawn, with horizons r_+ = 8R/25 and r_- = 9R/50.

Lindquist's parameters are a = sqrt(d_1^2 - R^2) = sqrt(3) R and cosh(mu_0) = d_1/R = 2, so
sinh(mu_0) = sqrt(3); the images are m_n = (-1)^{n+1} (sqrt(3)/2) R / sinh(n mu_0) at
z = +-d_n, d_n = sqrt(3) R coth(n mu_0), and sinh(n mu_0) grows as (2 + sqrt(3))^n / 2, so twelve
images leave V within 1e-7 of the whole series everywhere outside the spheres. V is harmonic outside
the spheres, so by the maximum principle it lies between 1 and V_0 there, largest on the shells.

The spacetime diagrams draw the axis through both spheres and the midplane between them in the
Cartesian chart, the stretch of the axis between the shells in the bispherical chart (eta = pi), the
axis outward from one shell in the polar chart (theta = pi), and the black hole between the shells
in Reissner and Nordstrom's chart. The conformal diagram is the black hole's maximal extension with
the two shells, one in an exterior and one in an exterior of the next universe, in its causal future.
The embedding diagram is the midplane z = 0 at T = 0, a surface of revolution about the axis, and the
curvature figure is the stored Kretschmann scalar inside, against r, and on the midplane, against the
distance from the axis.
"""
import math

import numpy as np

METRIC = "charged_shell_wormhole"
R, D1, M1 = 1.0, 2.0, 0.5
M, Q = 0.25, 0.24
A = math.sqrt(D1 ** 2 - R ** 2)
MU0 = math.acosh(D1 / R)
IMAGES = 12
SHELL = R + M1
RP, RM = M + math.sqrt(M ** 2 - Q ** 2), M - math.sqrt(M ** 2 - Q ** 2)

PARAMS = {"cartesian": {"R": 1, "d_1": 2, "m_1": "1/2"},
          "bispherical": {"a": "sqrt(3)", "mu_0": "acosh(2)", "m_1": "1/2"},
          "polar": {"R": 1, "d_1": 2, "m_1": "1/2"},
          "interior": {"m": "1/4", "q": "6/25", "R": 1, "m_1": "1/2"}}


def image_terms(n):
    """(m_n, d_n) of the n-th image as text in twenty digits, R = 1; written exactly, as
    sqrt(3)/(2 sinh(n acosh 2)), sympy carries the radicals through every derivative and a slice
    of the midplane takes minutes."""
    import mpmath
    mpmath.mp.dps = 30
    mu0 = mpmath.acosh(2)
    mn = (-1) ** (n + 1) * mpmath.sqrt(3) / (2 * mpmath.sinh(n * mu0))
    dn = mpmath.sqrt(3) / mpmath.tanh(n * mu0)
    return (mpmath.nstr(mn, 20), mpmath.nstr(dn, 20))


def cartesian_text(rho2="x**2 + y**2", z="z"):
    """V as sympy text, the first IMAGES images of the founders' (6), with the squared distance from
    the axis and the height given as text."""
    terms = []
    for n in range(1, IMAGES + 1):
        mn, dn = image_terms(n)
        terms.append(f"({mn})*(1/sqrt({rho2} + ({z} - {dn})**2) + 1/sqrt({rho2} + ({z} + {dn})**2))")
    return "1 + " + " + ".join(terms)


SHIFT = {n: repr(2 * n * MU0) for n in range(-IMAGES, IMAGES + 1)}
V_CARTESIAN = cartesian_text()
V_POLAR = cartesian_text("(rho*sin(theta))**2", "(rho*cos(theta) - 2)")
V_BISPHERICAL = ("1 + (1/2)*(1 + sqrt(cosh(mu) - cos(eta))*("
                 + " + ".join(f"(-1)**{n + 1}/sqrt(cosh(mu + {SHIFT[n]}) - cos(eta))"
                              for n in range(-IMAGES, IMAGES + 1)) + "))")
INPUT = ("Two spheres of radius $R$ centred at $z = \\pm 2R$, each of mass parameter $m_1 = R/2$, with $V$ "
         "the series of the first twelve image charges.")
SETTINGS = "$R = 1$, the unit of every length, $d_1 = 2R$, and $m_1 = R/2$, so that $V = 3/2$ on both shells."
INSIDE = "$m = R/4$, $q = 6R/25$, $R = 1$, and $m_1 = R/2$, so that the shells sit at $r = 3R/2$."


def v_numeric(rho, z):
    """V at the distance rho from the axis and the height z, to twelve images, R = 1."""
    rho, z = np.asarray(rho, dtype=float), np.asarray(z, dtype=float)
    total = np.ones(np.broadcast(rho, z).shape)
    for n in range(1, IMAGES + 1):
        mn = (-1) ** (n + 1) * M1 * math.sinh(MU0) / math.sinh(n * MU0)
        dn = A / math.tanh(n * MU0)
        total = total + mn * (1 / np.sqrt(rho ** 2 + (z - dn) ** 2) + 1 / np.sqrt(rho ** 2 + (z + dn) ** 2))
    return total


# ---------------------------------------------------------------- spacetime diagrams

def register_null(nr):
    side = nr.SIDEWAYS
    shells = (("shell", "r", "1", "the shells"), ("shell", "r", "3", "the shells"),
              ("shell", "r", "-1", "the shells"), ("shell", "r", "-3", "the shells"))
    nr.DIAGRAMS.append(nr.Diagram(
        METRIC, "cartesian", "axis", "the axis", ("T", "z"), (-4, 4, -4, 4), "$z/R$", "$cT/R$",
        PARAMS["cartesian"], {"x": "0", "y": "0"}, tau="T", families=side, functions={"V": V_CARTESIAN}, input=INPUT,
        where="((z - 2)**2 - 1)*((z + 2)**2 - 1)", lines=shells))
    nr.DIAGRAMS.append(nr.Diagram(
        METRIC, "cartesian", "midplane", "the midplane", ("T", "x"), (-4, 4, -4, 4), "$x/R$", "$cT/R$",
        PARAMS["cartesian"], {"y": "0", "z": "0"}, tau="T", families=side, functions={"V": V_CARTESIAN}, input=INPUT))
    reach = A * MU0
    nr.DIAGRAMS.append(nr.Diagram(
        METRIC, "bispherical", "axis", "between the shells", ("T", "\\mu"), (-reach, reach, -reach, reach),
        "$a\\mu/R$", "$cT/R$", PARAMS["bispherical"], {"eta": "pi", "phi": "0"}, tau="T", families=side,
        to_display=((0, A), (1, 0)), functions={"V": V_BISPHERICAL}, input=INPUT,
        lines=(("shell", "r", "acosh(2)", "the shells $\\mu = \\pm\\mu_0$"),
               ("shell", "r", "-acosh(2)", "the shells $\\mu = \\pm\\mu_0$"))))
    nr.DIAGRAMS.append(nr.Diagram(
        METRIC, "polar", "axis", "away from the other sphere", ("T", "\\rho"), (1, 5, -2, 2), "$\\rho/R$", "$cT/R$",
        PARAMS["polar"], {"theta": "pi", "phi": "0"}, tau="T", functions={"V": V_POLAR}, input=INPUT,
        lines=(("shell", "r", "1", "the shell $\\rho = R$"),)))
    nr.DIAGRAMS.append(nr.Diagram(
        METRIC, "interior", "radial", "$t$ and $r$", ("t", "r"), (0, 1.5, -0.75, 0.75), "$r/R$", "$ct/R$",
        PARAMS["interior"], nr.EQUATOR, orient="ingoing", areal=True,
        lines=(("shell", "r", "3/2", "the shells, $r = R + m_1$"),)))
    nr.CAPTIONS[(METRIC, "cartesian", "axis")] = [
        "The plane of $T$ and $z$ ($x = y = 0$), the axis through the centres of both spheres, hatched inside "
        "them. The cones are narrowest at the shells, where the potential takes its largest value "
        "$V = 1 + m_1/R$, and open toward $z = 0$ and toward infinity, where $V \\to 1$."]
    nr.CAPTIONS[(METRIC, "cartesian", "midplane")] = [
        "The plane of $T$ and $x$ ($y = 0$, $z = 0$), the midplane between the spheres. The cones are "
        "narrowest at $x = 0$, nearest both shells, and open outward as $V$ falls toward 1."]
    nr.CAPTIONS[(METRIC, "bispherical", "axis")] = [
        "The plane of $T$ and $\\mu$ ($\\eta = \\pi$, $\\phi = 0$), the stretch of the axis between the "
        "shells $\\mu = \\pm\\mu_0$, drawn against $a\\mu$. The cones widen toward the shells, where the "
        "coordinate's scale $a/W$ is smallest."]
    nr.CAPTIONS[(METRIC, "polar", "axis")] = [
        "The plane of $T$ and $\\rho$ ($\\theta = \\pi$, $\\phi = 0$), the axis outward from one shell on "
        "the side away from the other. The cones are narrowest at the shell $\\rho = R$ and open outward "
        "as $V$ falls toward 1."]
    nr.CAPTIONS[(METRIC, "interior", "radial")] = [
        "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$) of the black hole between the shells "
        "($m = R/4$, $q = 6R/25$), out to the shells at $r = R + m_1$. The cones close "
        "at the horizons $r_\\pm$, and inside $r_+$ the future is taken from the ingoing rays."]


def slice_table(slices):
    """The moment T = 0 of the midplane z = 0, the embedding view's: across the midplane's plane,
    and the one event T = 0, z = 0 where it meets the axis, in the Cartesian and bispherical charts."""
    def midplane():
        return slices.one(METRIC, lambda m: slices.across(0.0, *m.reach("cartesian", "x")), view_id="midplane")

    def axis():
        m = slices.moments(METRIC, "midplane", label="$T = 0$, $z = 0$")[0]
        return [slices.Mark(m, points=[(0.0, 0.0)])]
    return {(METRIC, "cartesian", "midplane"): midplane, (METRIC, "cartesian", "axis"): axis,
            (METRIC, "bispherical", "axis"): axis}


# ---------------------------------------------------------------- embedding diagram

def register_embedding(em):
    em.DRAWN[METRIC] = lambda ck, src: embedding_views(ck, src, em)
    em.CAPTIONS[(METRIC, "midplane")] = [
        "The midplane $z = 0$ between the spheres at the moment $T = 0$, turned about the axis through them, "
        "each circle at the distance $s$ from the axis having circumference $2\\pi sV$."]


def embedding_views(ck, src, em):
    """The midplane z = 0 at T = 0: V^2 (ds^2 + s^2 dphi^2) with V = V(s, 0) on it, a surface of
    revolution with circumference radius s V and dZ/ds = sqrt(V^2 - (V + s V')^2), real since
    V + s V' lies between 0 and V there; flat at the axis, where V' = 0, and drawn out to s = 4R."""
    name = "Schein-Aichelburg wormhole"
    sl = em.Slice(src, METRIC, "cartesian", "x", None, {"T": 0, "z": 0}, PARAMS["cartesian"],
                  functions={"V": V_CARTESIAN}, turn="y")
    top = 4.0
    size = 2 * float(sl.rho_at(top))
    s = np.linspace(1e-3, top, 400)
    v = v_numeric(s, 0.0)
    dv = np.gradient(v, s)
    ck.add(f"{name}: V + s V' lies between 0 and V on the midplane", float(max(0.0, np.max(-(v + s * dv)))
                                                                           + max(0.0, np.max(v + s * dv - v))), 0.0)
    plane = em.Piece("midplane", "sheet", sl, 0.0, top, 0.0, 1,
                     (("axis", "the axis $s = 0$, midway between the spheres"),
                      ("edge", "the surface runs on to $s \\to \\infty$")),
                     [(r, "r", None) for r in (1.0, 2.0, 3.0)] + [(top, "r", None)], size)
    ck.isometry(f"{name}, the midplane", plane)
    ck.radius(f"{name}, the midplane, rho = s V", plane, lambda r: r * v_numeric(r, 0.0), size)
    sheet = em.Surface([plane])
    fig = em.figure_of([sheet], {"sheet": "cover"}, size)
    em.ring_label(fig, [0, 0, 0], *plane.at(2.0), "$s = 2R$")
    em.ring_label(fig, [0, 0, 0], *plane.at(top), "$4R$")
    fig.legend("fill", "cover", "the plane $z = 0$ at one moment, which $T$, $x$ and $y$ cover")
    fig.legend("line", "r", "$s$ constant, at $1$, $2$, $3$ and $4$ times $R$")
    fig.legend("line", "meridian", "the angle about the axis constant, every $15°$")
    return [em.view("midplane", "The midplane", "$R$", [sheet], fig.done(),
                    settings="$R = 1$, $d_1 = 2R$, $m_1 = R/2$, and $s = \\sqrt{x^2 + y^2}$.")]


# ---------------------------------------------------------------- conformal diagram

def register_conformal(cf):
    cf.DRAWN[METRIC] = lambda ck, src: conformal_views(ck, src, cf)
    cf.CAPTIONS[(METRIC, "tower")] = [
        "The black hole between the shells maximally extended, each point a sphere, with one shell in an "
        "exterior and the other in an exterior of the next universe, in its causal future."]


def conformal_views(ck, src, cf):
    """The interior chart on the tower of Reissner-Nordstrom, horizons at the roots of its published
    g^rr, with the two shells r = R + m_1 drawn in the exterior I and in the exterior I above the
    inner horizon, the next universe's, every point of which lies to the future of every point of
    the first: p and q both grow along future causal curves, and the cell above has the larger p and q."""
    import sympy as sp
    name = "Schein-Aichelburg interior"
    pl = cf.Plane(src, METRIC, "interior", ("t", "r"), cf.EQUATOR, PARAMS["interior"])
    roots = sorted(sp.solve(sp.numer(sp.together(pl.gi[1, 1])), pl.x1), reverse=True)
    T = cf.Tower(-pl.g[0, 0], pl.x1, roots)
    rp, rm = T.rf
    cf.tower_checks(ck, name, pl, T, 0.01, 15)
    ck.limit(f"{name}: the horizons are the roots of the published g^rr", [rp, rm], [RP, RM], 1e-12)
    ck.diverges(f"{name}: the Kretschmann scalar diverges at r = 0", pl.kretschmann(0, 1e-2), pl.kretschmann(0, 1e-3))
    # In units of r_-, since K r^8 is of order q^4 and the horizons lie well inside R.
    ck.finite(f"{name}: the Kretschmann scalar is finite at both horizons, K r_-^4",
              pl.kretschmann(np.zeros(2), np.array([rp, rm])) * rm ** 4)
    t = cf.spread(-np.inf, np.inf, 500, 10)
    low = T.pq("I", t, np.full_like(t, SHELL))
    high = cf.reflect(*T.pq("I", t, np.full_like(t, SHELL)), True)
    ck.limit(f"{name}: every point of the upper shell has larger p and q than every point of the lower",
             [float(np.min(high[0]) - np.max(low[0]) > 0), float(np.min(high[1]) - np.max(low[1]) > 0)], [1.0, 1.0], 0.5)

    D = cf.TowerDrawing(T, True, 0.0)
    PI = cf.PI
    box = [-PI - 0.45, PI + 0.45, -PI - 0.1, 3 * PI + 0.1]
    times = [c / T.kp for c in (-1.6, -0.6, 0, 0.6, 1.6)]
    rI = cf.nice_all(cf.even_radii(T, "I", 3, rp, SHELL), [rp])
    rII = cf.nice_all(cf.even_radii(T, "II", 3, rm, rp), [rm, rp])
    rIII = cf.nice_all(cf.even_radii(T, "III", 2, 0, rm, xmax=cf.HALF), [0, rm])
    grids = {"I": (rI, times), "II": (rII, times), "IV": (rII, times), "III": (rIII, times)}
    v = cf.View("tower", "Between the shells", box, "interior")
    D.draw(v, grids, cover=[("I", False), ("II", False), ("III", False)])
    D.labels(v)
    # Beyond each shell the exterior is cut away: the region between the shell, which runs from the
    # corner i^- = (-pi/2, 0) of cell I to its corner i^+ = (0, pi/2), and spatial infinity (-pi/2, pi/2).
    for up in (False, True):
        p, q = cf.reflect(*T.pq("I", t, np.full_like(t, SHELL)), up)
        corner = cf.reflect(np.array([-cf.HALF]), np.array([cf.HALF]), up)
        v.fill("removed", [cf.point(a, b) for a, b in zip(p, q)] + [cf.point(corner[0][0], corner[1][0])])
        D.curve(v, "surface", "I", t, np.full_like(t, SHELL), up=up)
    v.set(fade={"top": 0.9, "bottom": 0.9}, settings="$m = R/4$, $q = 6R/25$, and $R + m_1 = 3R/2$.")
    v.legend("cover", "covered by $t$ and $r > 0$")
    v.legend("surface", "the shells $r = R + m_1$")
    v.legend("removed", "cut away beyond each shell")
    v.legend("r", "$r$ constant")
    v.legend("t", "$t$ constant")
    v.legend("horizon", "the horizons $r_\\pm$")
    v.legend("singular", "the singularity $r = 0$")
    v.legend("scri", "null infinity $\\mathscr{I}^\\pm$")
    return [v]


# ---------------------------------------------------------------- curvature figure

def register_curvature(cv):
    cv.DRAWN[METRIC] = lambda ck, src: curvature_views(ck, src, cv)


def curvature_views(ck, src, cv):
    """The stored Kretschmann scalars, read from the metric file and never recomputed: the interior
    chart's against r from R/10 to the shells, held to Reissner and Nordstrom's
    8(6m^2r^2 - 12mq^2r + 7q^4)/r^8 at points; and the Cartesian chart's on the midplane z = 0
    against x, with V the twelve images, where the stored Ricci scalar vanishes."""
    figures = []
    lo, hi = 0.1, SHELL
    xaxis = cv.axis("$r/R$", lo, hi)
    scale = {"kind": "log", "from": 1e-1, "to": 1e7}
    panels = [{"id": "K", "name": "$K\\,R^4$", "aspect": 2, "scale": scale, "ticks": cv.value_ticks(scale)}]
    cv.shape_checks(ck, METRIC, [p["aspect"] for p in panels])
    inside = ["m", "q", "R", "m_1"]
    K = cv.nr_published(src, METRIC, "interior", "kretschmann", ["r"], inside)
    Ric = cv.nr_published(src, METRIC, "interior", "ricci_scalar", ["r"], inside)
    for at in (0.12, RM, RP, 1.0, SHELL):
        want = 8 * (6 * M ** 2 * at ** 2 - 12 * M * Q ** 2 * at + 7 * Q ** 4) / at ** 8
        ck(f"{METRIC} interior K at r = {at:.3g}R is Reissner and Nordstrom's",
           abs(float(K.exact([at, M, Q, R, M1])) - want) < 1e-9 * abs(want))
        ck(f"{METRIC} interior R = 0 at r = {at:.3g}R", abs(float(Ric.exact([at, M, Q, R, M1]))) < 1e-12)
    xs = np.linspace(lo, hi, 400)
    values = K(xs, M, Q, R, M1)
    ck(f"{METRIC} interior K lies within its scale", bool(np.all(values > scale["from"]) and np.all(values < scale["to"])))
    curves = cv.line_curves(ck, METRIC, "K", "K", K, "r", xaxis, scale, cv.PLOT_W / 2,
                            {"m": M, "q": Q, "R": R, "m_1": M1})
    marks = [{"class": "horizon", "x": RP}, {"class": "horizon", "x": RM}]
    figures.append({
        "id": "interior", "label": "Between the shells", "kind": "line", "system": "interior", "x": xaxis,
        "panels": panels, "curves": curves, "marks": marks,
        "legend": [["line", "K", "$K$"], ["line", "horizon", "the horizons $r_\\pm$"]],
        "caption": ["The Kretschmann scalar $K$ of the black hole between the shells ($m = R/4$, $q = 6R/25$) "
                    "against $r$, finite at both horizons and growing without bound toward the singularity "
                    "$r = 0$. The Ricci scalar is zero."],
        "settings": INSIDE, "clip": "$K\\,R^4$ logarithmic, from $r = R/10$."})

    lo, hi = 0.0, 4.0
    xaxis = cv.axis("$x/R$", lo, hi)
    scale = {"kind": "log", "from": 1e-5, "to": 1}
    panels = [{"id": "K", "name": "$K\\,R^4$", "aspect": 2, "scale": scale, "ticks": cv.value_ticks(scale)}]
    cv.shape_checks(ck, METRIC, [p["aspect"] for p in panels])
    outside = ["R", "d_1", "m_1"]
    K = cv.Held(cv.nr_published(src, METRIC, "cartesian", "kretschmann", ["x", "y", "z"], outside,
                                functions={"V": V_CARTESIAN}), {"y": 0, "z": 0})
    Ric = cv.Held(cv.nr_published(src, METRIC, "cartesian", "ricci_scalar", ["x", "y", "z"], outside,
                                  functions={"V": V_CARTESIAN}), {"y": 0, "z": 0})
    for at in (0.0, 1.0, 3.0):
        ck(f"{METRIC} midplane R = 0 at x = {at}R, V harmonic", abs(float(Ric.exact([at, R, D1, M1]))) < 1e-9)
    xs = np.linspace(lo, hi, 200)
    values = K(xs, R, D1, M1)
    ck(f"{METRIC} midplane K lies within its scale", bool(np.all(values > scale["from"]) and np.all(values < scale["to"])))
    curves = cv.line_curves(ck, METRIC, "K", "K", K, "x", xaxis, scale,
                            cv.PLOT_W / 2, {"R": R, "d_1": D1, "m_1": M1})
    figures.append({
        "id": "midplane", "label": "The midplane", "kind": "line", "system": "cartesian", "x": xaxis,
        "panels": panels, "curves": curves, "marks": [],
        "legend": [["line", "K", "$K$"]],
        "caption": ["The Kretschmann scalar $K$ on the midplane $z = 0$ between the spheres against $x$ "
                    "($y = 0$), the same at every $T$. The Ricci scalar vanishes, since $V$ is harmonic."],
        "settings": SETTINGS, "clip": "$K\\,R^4$ logarithmic."})
    return figures
