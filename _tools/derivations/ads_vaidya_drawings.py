"""The drawings of a shell of light falling in anti-de Sitter space.

Every drawing is of a thin shell about v = 0, in units of the anti-de Sitter radius, with the final
horizon at r = l in every dimension: m runs from -1 to 1 in three dimensions, 2Gm/c^2 from 0 to 2l
in four, m from 0 to l with a flat front, and m from 0 to 2l^2 in five. Before the shell the
spacetime is anti-de Sitter space and after it the black hole. The spacetime diagrams draw the
shell with the thickness v_0 = l/100c of Balasubramanian and coauthors' thin shell, their equation
(2), m = (M/2)(1 + tanh(v/v_0)), and in three dimensions Hubeny, Rangamani and Takayanagi's profile
from -1 to m_0, which at m_0 = 1 is tanh(v/v_0): the tracer follows a ray through a smooth shell,
and across a step it loses the rays that cross deep inside the horizon. The conformal diagram is
the limit of no thickness (Anous, Hartman, Rovai and Sonner's equation (2.4)), where both sides
have null coordinates in closed form.

The ingoing charts are drawn against v itself. In anti-de Sitter space the outgoing rays are
r = l tan((cv - cv_0)/2l), which reach the boundary in a finite v, so no line v - a r with a > 0 is
a time over the whole plane; v - (l/c) arctan(r/l) is one everywhere, and orients the cones. The
Poincare charts are drawn against cv + z, which is the Poincare time before the shell and a time
everywhere, so every ray there runs at 45 degrees before the shell.

The event horizon is the outgoing ray that ends on r = l. In the round charts it leaves the centre
at cv = -2 l arctan(1) = -pi l/2 for a shell of no thickness, Hubeny, Liu and Rangamani's t_h =
t_s - 2 arctan(r_+), and at -1.571 l as drawn; behind a flat front it is r = 2l^2/(2l - cv), which
comes up from r = 0 in the far past. ads_vaidya.md beside this file is the derivation."""
import math

import numpy as np

METRIC = "ads_vaidya"

# Each chart's mass function, a thin shell of thickness v_0 = l/100c about v = 0, Balasubramanian and
# coauthors' equation (2) and, in three dimensions, Hubeny, Rangamani and Takayanagi's profile from
# -1 to m_0 at m_0 = 1: as sympy reads it, and the same in words.
V0 = "100*v"
STEP = {
    "ingoing_three": (f"tanh({V0})", "$m(v) = \\tanh(v/v_0)$"),
    "poincare_three": (f"(1 + tanh({V0}))/2", "$m(v) = (1 + \\tanh(v/v_0))/2\\ell^2$"),
    "ingoing_four": (f"(1 + tanh({V0}))/2", "$m(v) = \\tfrac{1}{2}M(1 + \\tanh(v/v_0))$, $2GM/c^2 = 2\\ell$,"),
    "flat_front": (f"(1 + tanh({V0}))/2", "$m(v) = \\tfrac{1}{2}\\ell(1 + \\tanh(v/v_0))$"),
    "poincare_four": (f"(1 + tanh({V0}))/2", "$m(v) = (1 + \\tanh(v/v_0))/2\\ell^3$"),
    "ingoing_five": (f"1 + tanh({V0})", "$m(v) = \\ell^2(1 + \\tanh(v/v_0))$"),
}
FIXED = {
    "ingoing_three": {"theta": "0"}, "poincare_three": {"x": "0"},
    "ingoing_four": {"theta": "pi/2", "phi": "0"}, "flat_front": {"x": "0", "y": "0"},
    "poincare_four": {"x": "0", "y": "0"}, "ingoing_five": {"chi": "pi/2", "theta": "pi/2", "phi": "0"},
}
HELD = {
    "ingoing_three": "$\\theta = 0$", "poincare_three": "$x = 0$", "ingoing_four": "$\\theta = \\pi/2$, $\\phi = 0$",
    "flat_front": "$x = y = 0$", "poincare_four": "$x = y = 0$",
    "ingoing_five": "$\\chi = \\theta = \\pi/2$, $\\phi = 0$",
}
ROUND = ("ingoing_three", "ingoing_four", "ingoing_five")
POINCARE = ("poincare_three", "poincare_four")
SHELL = "the shell of light, $v = 0$"


def shell_input(system):
    """The declared shell of a chart, as its spacetime diagram and its curvature map both draw it."""
    shell = "A thin shell of light of mass $M$" if system == "ingoing_four" else "A thin shell of light"
    return shell + " falling in along $v = 0$, with " + STEP[system][1] + " and $cv_0 = \\ell/100$."


def register_null(nr):
    for system in STEP:
        function = STEP[system][0]
        common = dict(params={"ell": 1, **({"G": 1} if system == "ingoing_four" else {})}, fixed=FIXED[system],
                      functions={"m": function}, lines=(("shell", "x0", "0", SHELL),),
                      # The horizon is seeded on r = l well after the shell, where the mass is the final one
                      # to rounding, and traced both ways.
                      marked=(("event", {"x0": "2/5", "r": "1"}, 0 if system in POINCARE else 1, "the event horizon"),),
                      input=shell_input(system))
        if system in POINCARE:
            # X = z, Y = cv + z, a time everywhere; the future inside the hole is toward larger z.
            row = nr.Diagram(METRIC, system, "shell", "a shell of light", ("v", "z"), (0, 3, -1.5, 1.5),
                             "$z/\\ell$", "$(cv + z)/\\ell$", to_display=nr.FINKELSTEIN_OUT, tau="v + z",
                             families=("toward the boundary", "into the interior"), **common)
        else:
            # The sphere's radius is r in every round chart, so the marginally trapped spheres are
            # where g^rr vanishes, which the view marks already.
            tau = "v + 1/r" if system == "flat_front" else "v - atan(r)"
            row = nr.Diagram(METRIC, system, "shell", "a shell of light", ("v", "r"), (0, 3, -2, 1),
                             "$r/\\ell$", "$cv/\\ell$", tau=tau, singular_runs=system != "ingoing_three", **common)
        nr.DIAGRAMS.append(row)
    for system in ROUND:
        nr.CAPTIONS[(METRIC, system, "shell")] = [
            f"The plane of $v$ and $r$ ({HELD[system]}), with the ingoing rays the horizontal lines of constant $v$. "
            "The event horizon leaves the centre at $cv = -1.57\\ell$, before the shell arrives, and lies on "
            "$r = \\ell$ after it."]
    nr.CAPTIONS[(METRIC, "flat_front", "shell")] = [
        f"The plane of $v$ and $r$ ({HELD['flat_front']}), with the ingoing rays the horizontal lines of constant "
        "$v$. The event horizon is $r = 2\\ell^2/(2\\ell - cv)$ before the shell and $r = \\ell$ after it."]
    for system in POINCARE:
        nr.CAPTIONS[(METRIC, system, "shell")] = [
            f"The plane of $v$ and $z$ ({HELD[system]}), with the boundary at $z = 0$. Every ray runs at 45° "
            "before the shell. The event horizon is $z = \\ell - cv/2$ before the shell and $z = \\ell$ after it."]
    # Before the shell every chart is anti-de Sitter space, where the rays have closed forms: the family
    # along the shell keeps v, and the other keeps arctan(r/l) - cv/2l, l/r + cv/2l or z + cv/2.
    before = lambda v, r: v < -0.1  # noqa: E731
    nr.CLOSED_FORMS.update({
        **{(METRIC, system, "shell"): (lambda v, r: v, lambda v, r: np.arctan(r) - v / 2, before) for system in ROUND},
        # Near r = 0 the flat front's m/r is as large as r^2 even in the shell's far tail, so its rays keep
        # 1/r + cv/2l only from r = l/20 out.
        (METRIC, "flat_front", "shell"): (lambda v, r: v, lambda v, r: 1 / r + v / 2,
                                          lambda v, r: (v < -0.1) & (r > 0.05)),
        **{(METRIC, system, "shell"): (lambda v, z: z + v / 2, lambda v, z: v, before) for system in POINCARE},
    })


# ---------------------------------------------------------------- the conformal diagrams
#
# Each conformal view is the shell of no thickness, Anous, Hartman, Rovai and Sonner's equation (2.4),
# at l = 1 with the final horizon r_h = l, so that both sides have null coordinates in closed form.
#
# Round front. Before the shell f = 1 + r^2, r_* = arctan r, and the plane of v and r is anti-de
# Sitter space's global strip, p, q = (v - 2 arctan r)/2, v/2, which puts the centre on X = 0, the
# boundary on X = pi/2 and the shell v = 0 on the null line q = 0, from (X, T) = (pi/2, -pi/2) to the
# centre at (0, 0). After it the hole has Kruskal's V = exp(k v) and U = UV(r)/V, UV(r) falling from
# UV(0) at r = 0 through 0 on the horizon to -B on the boundary. Each outgoing ray keeps the p it has
# where it crosses the shell, at the radius r_s with UV(r_s) = U, so p = P(U) = -arctan r_s(U), and
# q = Q(V) = P(-B/V) + pi/2 puts the boundary UV = -B on X = pi/2 after the shell too. The event
# horizon U = 0 is p = -pi/4, the line T - X = -pi/2 in every dimension, from the centre at
# cv = -pi l/2 (Hubeny, Liu and Rangamani's t_h = t_s - 2 arctan r_+) through the shell at r = l to
# the corner (pi/2, 0) where the boundary meets r = 0. In three dimensions UV = (1 - r)/(1 + r),
# P(U) = arctan U - pi/4 and Q(V) = arctan V - pi/4, the hole is the square of Banados, Henneaux,
# Teitelboim and Zanelli moved down by pi/2, and r = 0 is the straight line T = 0.
#
# Flat front. Before the shell f = r^2, r_* = -1/r = -z, and the plane is the Poincare patch,
# p, q = arctan(v + 2z) - pi/4, arctan v + pi/4, with the boundary z = 0 on X = pi/2 and the two
# Poincare horizons v -> -infinity and z -> infinity, where the chart ends, on the null lines through
# (-pi/2, 0). After it each outgoing ray keeps its p, P(U) = arctan(2 z_s(U)) - pi/4, and
# Q(V) = P(-B/V) + pi/2 again puts the boundary on X = pi/2.

ROUND_VIEWS = {
    "ingoing_three": dict(dims=3, fixed=FIXED["ingoing_three"], params={"ell": 1}, before="-1", after="1",
                          label="Three dimensions"),
    "ingoing_four": dict(dims=4, fixed=FIXED["ingoing_four"], params={"ell": 1, "G": 1}, before="0", after="1",
                         label="Four dimensions"),
    "ingoing_five": dict(dims=5, fixed=FIXED["ingoing_five"], params={"ell": 1}, before="0", after="2",
                         label="Five dimensions"),
}
PLANAR_VIEWS = {
    "poincare_three": dict(dims=3, fixed=FIXED["poincare_three"], radius="z", label="Poincaré, three dimensions"),
    "flat_front": dict(dims=4, fixed=FIXED["flat_front"], radius="r", label="Flat front"),
    "poincare_four": dict(dims=4, fixed=FIXED["poincare_four"], radius="z", label="Poincaré, four dimensions"),
}


def btz_uv(r):
    """Kruskal's UV of the hole of three dimensions, f = r^2 - 1, k = 1, with r_* -> 0 far out."""
    r = np.asarray(r, dtype=float)
    return (1 - r) / (1 + r)


def five_uv(r):
    """UV of the hole of five dimensions, f = (r^2 - 1)(r^2 + 2)/r^2, k = 3, with r_* -> 0 far out:
    r_* = ln|(r - 1)/(r + 1)|/6 + (sqrt 2/3)(arctan(r/sqrt 2) - pi/2)."""
    r = np.asarray(r, dtype=float)
    return (1 - r) / (1 + r) * np.exp(2 * np.sqrt(2) * (np.arctan(r / np.sqrt(2)) - np.pi / 2))


def brane_uv(r):
    """UV of the black brane of four dimensions, f = (r^3 - 1)/r, k = 3/2, with r_* -> 0 far out:
    r_* = ln|r - 1|/3 - ln(r^2 + r + 1)/6 + (arctan((2r + 1)/sqrt 3) - pi/2)/sqrt 3."""
    r = np.asarray(r, dtype=float)
    return (1 - r) / np.sqrt(r * r + r + 1) * np.exp(np.sqrt(3) * (np.arctan((2 * r + 1) / np.sqrt(3)) - np.pi / 2))


def shell_angle(uv, U, n=64):
    """arctan r_s for the radius r_s at which UV(r_s) = U, UV falling with r, by bisection on (0, pi/2)."""
    U = np.asarray(U, dtype=float)
    a, b = np.zeros(U.shape), np.full(U.shape, np.pi / 2)
    for _ in range(n):
        m = 0.5 * (a + b)
        beyond = uv(np.tan(m)) > U
        a, b = np.where(beyond, m, a), np.where(beyond, b, m)
    return 0.5 * (a + b)


class Shell:
    """The maps of one view: inside(v, x) before the shell and outside(v, x) after it, x the chart's r
    for a round front and z = l^2/r for a flat one."""

    def __init__(self, uv, kappa, B, planar):
        self.uv, self.k, self.B, self.planar = uv, kappa, B, planar

    def P(self, U):
        s = shell_angle(self.uv, U)
        return np.arctan(2 / np.tan(s)) - np.pi / 4 if self.planar else -s

    def Q(self, V):
        return self.P(-self.B / np.asarray(V, dtype=float)) + np.pi / 2

    def inside(self, v, x):
        v, x = np.asarray(v, dtype=float), np.asarray(x, dtype=float)
        if self.planar:
            return np.arctan(v + 2 * x) - np.pi / 4, np.arctan(v) + np.pi / 4
        return (v - 2 * np.arctan(x)) / 2, v / 2

    def outside(self, v, x):
        v, x = np.asarray(v, dtype=float), np.asarray(x, dtype=float)
        V = np.exp(self.k * v)
        r = 1 / x if self.planar else x
        return self.P(self.uv(r) / V), self.Q(V)


def round_moment_lines(sh, T0, lo, hi):
    """The moment cv - l arctan(r/l) = T0 l of the round chart of three dimensions, from lo to hi in
    r: inside the shell, r < tan(-T0), it is the strip's line T = T0, and outside it each point is
    placed by the hole's map at v = T0 + arctan r."""
    R = np.tan(-T0) if T0 < 0 else 0.0
    lines = []
    if R > lo:
        r = np.linspace(lo, min(R, hi), 200)
        lines.append(sh.inside(T0 + np.arctan(r), r))
    if R < hi:
        r = max(lo, R) + (hi - max(lo, R)) * np.linspace(0, 1, 400) ** 2
        lines.append(sh.outside(T0 + np.arctan(r), r))
    return lines


def conformal_views(ck, src, cf):
    """One view for each chart, each the whole of the plane its spacetime diagram draws, with the
    shell's two sides joined by the maps of Shell; the docstring of the block above is the derivation."""
    HALF, Q4 = np.pi / 2, np.pi / 4
    views = []
    for system, spec in ROUND_VIEWS.items():
        plane = ("v", "r")
        inner = cf.Plane(src, METRIC, system, plane, spec["fixed"], spec["params"], functions={"m": spec["before"]})
        hole = cf.Plane(src, METRIC, system, plane, spec["fixed"], spec["params"], functions={"m": spec["after"]})
        name = f"ads_vaidya {system}"
        # A future directed timelike vector at every r: (1, -a) has the length -f - 2a, and -f/2 is
        # below 1/r^2 in every round chart, inside the horizon too.
        future = lambda w, r: (1, -60 - 1 / np.asarray(r) ** 2)  # noqa: E731
        if spec["dims"] == 4:
            tower = cf.AdSHoleTower(hole.gi[1, 1], hole.x1)
            ck.limit(f"{name}: the horizon is the positive root of the published g^rr, at r_h = l", tower.rf, [1.0], 1e-12)
            sh = Shell(tower.UV, tower.kp, tower.B, False)
        else:
            sh = Shell(btz_uv if spec["dims"] == 3 else five_uv, 1.0 if spec["dims"] == 3 else 3.0, 1.0, False)
            grr = hole.lambdify([hole.gi[1, 1]])
            r = ck.uniform(0.02, 20, 400)
            r = r[np.abs(r - 1) > 1e-3]
            dr = 1e-6 * r
            slope = (np.log(np.abs(sh.uv(r + dr))) - np.log(np.abs(sh.uv(r - dr)))) / (2 * dr)
            f = hole.call(grr, np.zeros_like(r), r)[0]
            ck.limit(f"{name}: d ln|UV|/dr = 2k/f, with f the published g^rr", np.max(np.abs(slope * f - 2 * sh.k)), 0, 1e-6)
        ck.chart(f"{name}, anti-de Sitter space before the shell", inner, sh.inside,
                 ck.uniform(-15, -0.01), ck.uniform(0.01, 40), future)
        ck.chart(f"{name}, the hole after the shell", hole, sh.outside,
                 ck.uniform(0.01, 8), ck.uniform(0.01, 40), future)
        r = ck.uniform(0.01, 40, 300)
        ck.limit(f"{name}: the two sides put the shell v = 0 at one place",
                 np.concatenate(sh.inside(0 * r, r)), np.concatenate(sh.outside(0 * r, r)), 1e-12)
        p_in, q_in = sh.inside(np.array([-9.0, -2.0, -0.5]), np.full(3, 1e12))
        p_out, q_out = sh.outside(np.array([0.5, 2.0, 6.0]), np.full(3, 1e12))
        ck.limit(f"{name}: r -> infinity lands on the boundary X = pi/2 on both sides",
                 np.concatenate([q_in - p_in, q_out - p_out]), HALF, 1e-9)
        ck.limit(f"{name}: the event horizon U = 0 meets the shell at r_h = l and reaches the centre at cv = -pi l/2",
                 [float(sh.P(0.0)), *sh.inside(-HALF, 0.0)], [-Q4, -Q4, -Q4], 1e-12)
        late = np.array([0.3, 1.0, 3.0])
        ck.limit(f"{name}: the boundary meets r = 0 at (X, T) = (pi/2, 0), where the horizon ends",
                 cf.point(sh.P(-sh.B / 1e12), sh.Q(1e12)) + cf.point(sh.P(sh.uv(0.0) / 1e12), sh.Q(1e12)),
                 [HALF, 0, HALF, 0], 1e-9)
        ck.finite(f"{name}: r = 0 before the shell is a regular centre",
                  inner.kretschmann(ck.uniform(-9, -0.1, 20), np.full(20, 1e-6)))
        if spec["dims"] == 3:
            p, q = sh.outside(late, np.full(3, 1e-12))
            ck.limit(f"{name}: r -> 0 after the shell lands on the straight line T = 0", p + q, 0, 1e-9)
            ck.finite(f"{name}: the curvature is finite at r = 0 after the shell, a singularity of the causal structure",
                      hole.kretschmann(late, np.full(3, 1e-9)))
        else:
            ck.diverges(f"{name}: the Kretschmann scalar diverges at r = 0 after the shell",
                        hole.kretschmann(1, 1e-2), hole.kretschmann(1, 1e-3))

        T0 = -2.45
        view = cf.View(system, spec["label"], [-0.42, HALF + 0.42, T0, 0.35], system)
        V = 1 + cf.spread(0, np.inf, 500, 12)
        sp_, sq_ = sh.P(sh.uv(0.0) / V), sh.Q(V)
        top = [cf.point(a, b) for a, b in zip(sp_, sq_)]
        top = [[0.0, 0.0]] + top + [[HALF, 0.0]]
        inside = [[0, T0], [HALF, T0], [HALF, -HALF], [0, 0]]
        outside = [[0, 0], [HALF, -HALF]] + top[::-1]
        view.fill("region", inside)
        view.fill("region", outside)
        view.fill("cover", outside)
        view.fill("cover2", inside)
        for r in (0.25, 0.5, 1.0, 2.0):
            X = float(np.arctan(r))
            view.line("r", [[[X, T0], [X, -X]]])
        w = cf.spread(0, np.inf, 500, 10)
        for r in (0.5, 0.8, 1.25, 2.0, 4.0):
            view.curve("r", *sh.outside(w, np.full_like(w, r)))
        view.segment("surface", (-HALF, 0), (0, 0))
        view.segment("event", (-Q4, -Q4), (-Q4, Q4))
        view.line("centre", [[[0, T0], [0, 0]]])
        view.line("boundary", [[[HALF, T0], [HALF, 0]]])
        view.line("singular", [top], zig=True)
        view.label_xt([0, -1.0], "$r = 0$", "r", "small", dx=-6)
        view.label_xt([Q4, 0.0], "$r = 0$", "b", "small", dy=-8)
        view.label_xt([HALF, -1.0], "$r \\to \\infty$", "l", "small", dx=6)
        point = {3: "a circle", 4: "a 2-sphere", 5: "a 3-sphere"}[spec["dims"]]
        view.legend("cover", "the black hole")
        view.legend("cover2", "anti-de Sitter space")
        view.legend("r", "$r$ constant")
        view.legend("surface", "the shell $v = 0$")
        view.legend("event", "the event horizon")
        view.legend("centre", "the centre $r = 0$")
        view.legend("boundary", "the boundary")
        view.legend("singular", "$r = 0$, " + ("a singularity of the causal structure" if spec["dims"] == 3
                                               else "where the curvature diverges"))
        view.set(fade={"top": 0, "bottom": 0.7}, settings=CONFORMAL_SETTINGS[system])
        if system == "ingoing_three":
            for m in cf.slices.moments(METRIC):
                lo, hi = m.reach(system, "r")
                view.slice(m, round_moment_lines(sh, m.time, lo, hi))
        views.append(view)

    for system, spec in PLANAR_VIEWS.items():
        x = spec["radius"]
        plane = ("v", x)
        inner = cf.Plane(src, METRIC, system, plane, spec["fixed"], {"ell": 1}, functions={"m": "0"})
        hole = cf.Plane(src, METRIC, system, plane, spec["fixed"], {"ell": 1}, functions={"m": "1"})
        name = f"ads_vaidya {system}"
        three = spec["dims"] == 3
        sh = Shell(btz_uv if three else brane_uv, 1.0 if three else 1.5, 1.0, True)
        to_z = (lambda c: c) if x == "z" else (lambda c: 1 / np.asarray(c, dtype=float))
        # (1, a) in v and z has the length (l/z)^2(-(1 - mz^d) - 2a), so a = 60 + z^3 is timelike at every z,
        # and in v and r, where (1, -a) has the length -f - 2a with -f/2 below 1/r^2, so is a = 60 + 1/r^2.
        future = ((lambda w, c: (1, 60 + np.asarray(c) ** 3)) if x == "z"
                  else (lambda w, c: (1, -60 - 1 / np.asarray(c) ** 2)))
        z = ck.uniform(0.02, 20, 400)
        z = z[np.abs(z - 1) > 1e-3]
        dz = 1e-6 * z
        slope = (np.log(np.abs(sh.uv(1 / (z + dz)))) - np.log(np.abs(sh.uv(1 / (z - dz))))) / (2 * dz)
        f = (1 / z ** 2 - 1) if three else (1 / z ** 2 - z)
        ck.limit(f"{name}: d ln|UV|/dz = -2k/(f z^2), f = r^2 - 1 or r^2 - l/r", np.max(np.abs(slope * f * z * z + 2 * sh.k)),
                 0, 1e-6)
        ck.chart(f"{name}, the Poincaré patch before the shell", inner, lambda w, c: sh.inside(w, to_z(c)),
                 ck.uniform(-15, -0.01), ck.uniform(0.02, 40), future)
        ck.chart(f"{name}, the {'hole' if three else 'black brane'} after the shell", hole,
                 lambda w, c: sh.outside(w, to_z(c)), ck.uniform(0.01, 8), ck.uniform(0.02, 40), future)
        c = ck.uniform(0.02, 40, 300)
        ck.limit(f"{name}: the two sides put the shell v = 0 at one place",
                 np.concatenate(sh.inside(0 * c, to_z(c))), np.concatenate(sh.outside(0 * c, to_z(c))), 1e-12)
        p_in, q_in = sh.inside(np.array([-9.0, -2.0, -0.5]), np.full(3, 1e-12))
        p_out, q_out = sh.outside(np.array([0.5, 2.0, 6.0]), np.full(3, 1e-12))
        ck.limit(f"{name}: the boundary z -> 0 lands on X = pi/2 on both sides",
                 np.concatenate([q_in - p_in, q_out - p_out]), HALF, 1e-9)
        horizon = float(sh.P(0.0))
        ck.limit(f"{name}: the event horizon U = 0 is z = l - cv/2 before the shell",
                 sh.inside(np.array([-6.0, -1.0]), 1 - np.array([-6.0, -1.0]) / 2)[0], horizon, 1e-12)
        ck.finite(f"{name}: the Poincaré horizon z -> infinity before the shell is regular",
                  inner.kretschmann(ck.uniform(-9, -0.1, 20), np.full(20, 1e6) if x == "z" else np.full(20, 1e-6)))
        far = (lambda c: c) if x == "z" else (lambda c: 1 / c)
        if three:
            ck.finite(f"{name}: the curvature is finite at r = 0 after the shell",
                      hole.kretschmann(np.array([0.3, 1.0, 3.0]), np.full(3, far(1e9))))
        else:
            ck.diverges(f"{name}: the Kretschmann scalar diverges at r = 0 after the shell",
                        hole.kretschmann(1, far(1e2)), hole.kretschmann(1, far(1e3)))

        V = 1 + cf.spread(0, np.inf, 500, 12)
        sp_, sq_ = sh.P(sh.uv(0.0) / V), sh.Q(V)
        T_top = float(2 * horizon + HALF)
        top = [[0.0, HALF]] + [cf.point(a, b) for a, b in zip(sp_, sq_)] + [[HALF, T_top]]
        view = cf.View(system, spec["label"], [-HALF - 0.3, HALF + 0.42, -np.pi - 0.25, T_top + 0.3], system)
        inside = [[HALF, -np.pi], [-HALF, 0], [0, HALF], [HALF, 0]]
        outside = [[HALF, 0], [0, HALF]] + top[1:]
        view.fill("region", inside)
        view.fill("region", outside)
        view.fill("cover", outside)
        view.fill("cover2", inside)
        w = cf.spread(-np.inf, 0, 500, 10)
        for zz in (0.25, 0.5, 1.0, 2.0):
            view.curve("r", *sh.inside(w, np.full_like(w, zz)))
        w = cf.spread(0, np.inf, 500, 10)
        for zz in (0.25, 0.5, 0.8, 1.25, 2.0):
            view.curve("r", *sh.outside(w, np.full_like(w, zz)))
        view.segment("surface", (-Q4, Q4), (Q4, Q4))
        view.segment("event", (horizon, -Q4), (horizon, float(sh.Q(1e12))))
        view.line("chartedge", [[[HALF, -np.pi], [-HALF, 0]], [[-HALF, 0], [0, HALF]]])
        view.line("boundary", [[[HALF, -np.pi], [HALF, T_top]]])
        view.line("singular", [top], zig=True)
        view.label_xt([HALF, -1.4], "$z \\to 0$" if x == "z" else "$r \\to \\infty$", "l", "small", dx=6)
        word = "z" if x == "z" else "r"
        view.legend("cover", "the " + ("black hole" if three else "black brane"))
        view.legend("cover2", "the Poincaré patch")
        view.legend("r", f"${word}$ constant")
        view.legend("surface", "the shell $v = 0$")
        view.legend("event", "the event horizon")
        view.legend("chartedge", "the Poincaré horizons")
        view.legend("boundary", "the boundary")
        view.legend("singular", ("$r = 0$, " if x == "r" else "$z \\to \\infty$, ")
                    + ("a singularity of the causal structure" if three else "where the curvature diverges"))
        view.set(settings=CONFORMAL_SETTINGS[system])
        views.append(view)
    return views


CONFORMAL_SETTINGS = {
    "ingoing_three": "$\\ell = 1$, $m = -1$ before the shell and $m = 1$ after it.",
    "ingoing_four": "$\\ell = G = 1$, $m = 0$ before the shell and $2Gm/c^2 = 2\\ell$ after it.",
    "ingoing_five": "$\\ell = 1$, $m = 0$ before the shell and $m = 2\\ell^2$ after it.",
    "poincare_three": "$\\ell = 1$, $m = 0$ before the shell and $m = 1/\\ell^2$ after it.",
    "flat_front": "$\\ell = 1$, $m = 0$ before the shell and $m = \\ell$ after it.",
    "poincare_four": "$\\ell = 1$, $m = 0$ before the shell and $m = 1/\\ell^3$ after it.",
}


def conformal_captions():
    point = {"ingoing_three": "a circle", "ingoing_four": "a 2-sphere", "ingoing_five": "a 3-sphere",
             "poincare_three": "a line", "flat_front": "a plane", "poincare_four": "a plane"}
    space = {"ingoing_three": "anti-de Sitter space of three dimensions",
             "ingoing_four": "anti-de Sitter space of four dimensions",
             "ingoing_five": "anti-de Sitter space of five dimensions",
             "poincare_three": "the Poincaré patch of anti-de Sitter space of three dimensions",
             "flat_front": "the Poincaré patch of anti-de Sitter space of four dimensions",
             "poincare_four": "the Poincaré patch of anti-de Sitter space of four dimensions"}
    out = {}
    for system in (*ROUND_VIEWS, *PLANAR_VIEWS):
        first = f"A thin shell of light falling into {space[system]}, each point in the diagram {point[system]}."
        if system in ROUND_VIEWS:
            second = ("The event horizon leaves the centre at $cv = -\\pi\\ell/2$ and meets the shell at "
                      "$r = \\ell$.")
        else:
            second = ("The event horizon meets the shell at "
                      + ("$z = \\ell$." if PLANAR_VIEWS[system]["radius"] == "z" else "$r = \\ell$."))
        out[(METRIC, system)] = [first + " " + second]
    return out


# ---------------------------------------------------------------- the embedding diagram
#
# The round chart of three dimensions, the shell of no thickness at l = 1, m = -1 before it and 1
# after it. A slice of constant v is null, so the moments are the slices of constant
# cv - l arctan(r/l) = T l, which is the global time of anti-de Sitter space before the shell. With
# v = T + arctan r the published metric pulls back to (2(1 + r^2) - f) dr^2/(1 + r^2)^2 + r^2 dtheta^2,
# spacelike at every r, inside the horizon too. Inside the shell, r < R = tan(-T), f = 1 + r^2 and the
# moment is the hyperbolic plane dr^2/(1 + r^2) + r^2 dtheta^2, whose circles grow faster than the
# distance out to them, so it is drawn in three dimensional Minkowski space, the hyperboloid
# Z = sqrt(1 + r^2) - 1, as anti_de_sitter's static slice is. Outside it f = r^2 - 1 and
# g_rr = (r^2 + 3)/(1 + r^2)^2, which is 1 at r = l, above 1 inside and below 1 outside: the moment
# is drawn in flat space between the shell and r = l and in Minkowski space beyond, the two lying
# level at r = l, as btz's moment t = 0 does at sqrt 2. Once the shell has reached the centre, T >= 0,
# the moment ends at the singularity r = 0 in the tip of a cone, dz/dr = sqrt 2. The event horizon is
# the circle r = 1/R inside the shell while R > l, from the strip's p = -pi/4, and r = l outside it.

EMBED_TOP = 3.0
EMBED_TIMES = [-1.25 + k / 20 for k in range(36)]
EMBED_SHOWN = (0, 10, 20, 35)


def moment_label(T):
    return f"$cv - \\ell\\arctan(r/\\ell) = {T:g}\\,\\ell$"


def embedding_views(ck, src, em):
    top, size = EMBED_TOP, 2 * EMBED_TOP

    def slice_at(T, m, space):
        return em.Slice(src, METRIC, "ingoing_three", "r", "\\theta", {}, {"ell": 1}, {"m": m},
                        along={"v": f"{T!r} + atan(r)"}, space=space)
    shell_end = ("crease", "the shell of light, where the surface folds")
    level = ("join", "the circle $r = \\ell$, where the surface lies level")
    rim = ("edge", "the surface runs on toward a light cone of Minkowski space, to $r \\to \\infty$")

    def rings(lo, hi):
        return [(x, "r", None) for x in (2.0, top) if lo < x <= hi]

    def moment(T):
        where = f"ads_vaidya, cv - l arctan(r/l) = {T:g} l"
        R = math.tan(-T) if T < 0 else 0.0
        pieces = []
        if R > 0:
            end = min(R, top)
            marks = ([(1 / R, "horizon", None)] if R > 1 else []) + ([(R, "surface", None)] if R < top else [])
            inside = em.Piece("inside", "sheet", slice_at(T, "-1", "minkowski"), 0.0, end, 0.0, 1,
                              (("axis", "the centre $r = 0$"), shell_end if R < top else rim),
                              marks + rings(0.0, end), size)
            ck.form(f"{where}, inside the shell the hyperboloid Z = sqrt(l^2 + r^2) - l", inside,
                    lambda r: np.sqrt(1 + r * r) - 1, size)
            pieces.append(inside)
        if R < top:
            z0 = pieces[-1].at(R)[1] if pieces else 0.0
            if R < 1:
                ends = (shell_end if R > 0 else ("apex", "the singularity $r = 0$, where the surface ends in the tip of a cone"),
                        level)
                near = em.Piece("between", "sheet", slice_at(T, "1", "flat"), R, 1.0, z0, 1, ends,
                                [(1.0, "horizon", None)], size)
                pieces.append(near)
                z0 = near.at(1.0)[1]
            far = em.Piece("beyond", "sheet", slice_at(T, "1", "minkowski"), max(R, 1.0), top, z0, 1,
                           (level if R < 1 else shell_end, rim), rings(max(R, 1.0), top), size)
            pieces.append(far)
            if R < 1:
                ck.join(f"{where}, flat space and Minkowski space at r = l", near, 1.0, far, 1.0)
            if R > 0:
                ck.add(f"{where}, the two sides meet at the shell: one point",
                       float(np.max(np.abs(np.array(pieces[0].at(R)) - pieces[1].at(R)))), em.JOIN)
        if R == 0:
            ck.add(f"{where}, the tip of the cone, dz/dr = sqrt 2 at r = 0",
                   abs(float(pieces[0].sl.defect_at(np.array([0.0]))[0]) - 2.0), 1e-12)
        # The rim of the drawing, r = 3l, stands at z = 0 at every moment.
        shift = -pieces[-1].z[-1]
        for p in pieces:
            p.z = p.z + shift
            ck.isometry(f"{where}, {p.id}", p)
        return em.Surface(pieces, label=moment_label(T), time=T)

    frames = [moment(T) for T in EMBED_TIMES]
    fig = em.movie_figure(frames, {"sheet": "cover"}, size)
    fig.legend("fill", "cover", "the moment")
    fig.legend("line", "r", "$r = 2\\ell$ and $3\\ell$")
    fig.legend("line", "surface", "the shell")
    fig.legend("line", "horizon", "the event horizon")
    fig.legend("line", "meridian", "$\\theta$ constant")
    return [em.view("shell", "The falling shell", "$\\ell$", [frames[k] for k in EMBED_SHOWN], fig.done(),
                    movie=em.movie(frames, "$cv - \\ell\\arctan(r/\\ell)$", EMBED_TIMES),
                    settings="$\\ell = 1$, with lengths in Minkowski space inside the shell and beyond $r = \\ell$.")]


def embedding_captions():
    return {(METRIC, "shell"): [
        "A thin shell of light falling into anti-de Sitter space of three dimensions, at moments of constant "
        "$cv - \\ell\\arctan(r/\\ell)$. Each moment is the hyperbolic plane inside the shell and ends in the tip "
        "of a cone after the shell reaches the centre."]}


def slice_marks(slices):
    """Each moment on the spacetime diagram of the round chart of three dimensions, v = T + arctan r
    in its (cv, r), out to the embedding's reach."""
    def lines(m):
        lo, hi = m.reach("ingoing_three", "r")
        r = np.linspace(lo, hi, slices.N)
        return [list(zip(m.time + np.arctan(r), r))]
    return slices.one(METRIC, lines)


HIDDEN_REASON = ("another spacetime than the shell embedded, which falls into anti-de Sitter space of three "
                 "dimensions with a round front")


# ---------------------------------------------------------------- the curvature figure
#
# The charts of four and five dimensions, each a map of its published Kretschmann scalar over the
# plane of its radius and v, with the thin shell of the spacetime diagrams declared: anti-de Sitter's
# constant before the shell and the hole's growth toward r = 0 after it. In three dimensions K is
# 12/l^4 at every event, and no map is drawn of it. The published Ricci scalar is the constant
# -d(d + 1)/l^2, since the light is traceless, which is checked and not drawn.

CURVATURE_VIEWS = {
    "ingoing_four": dict(x="r", params={"G": 1.0, "ell": 1.0}, ads=24.0, hole=lambda r, m: 48 * m * m / r ** 6,
                         ricci=-12.0, mass=1.0, name="the hole of four dimensions"),
    "flat_front": dict(x="r", params={"ell": 1.0}, ads=24.0, hole=lambda r, m: 12 * m * m / r ** 6,
                       ricci=-12.0, mass=1.0, name="the black brane"),
    "poincare_four": dict(x="z", params={"ell": 1.0}, ads=24.0, hole=lambda z, m: 12 * m * m * z ** 6,
                          ricci=-12.0, mass=1.0, name="the black brane"),
    "ingoing_five": dict(x="r", params={"ell": 1.0}, ads=40.0, hole=lambda r, m: 72 * m * m / r ** 8,
                         ricci=-20.0, mass=2.0, name="the hole of five dimensions"),
}
CURVATURE_SCALE = {"kind": "log", "from": 10.0, "to": 1e4}
CURVATURE_CELLS = 100


def curvature_views(ck, src, cv):
    out = []
    n = CURVATURE_CELLS
    for system, spec in CURVATURE_VIEWS.items():
        x = spec["x"]
        name = f"ads_vaidya {system}"
        functions = {"m": STEP[system][0]}
        names = list(spec["params"])
        K = cv.nr_published(src, METRIC, system, "kretschmann", [x, "v"], names, functions=functions)
        R = cv.nr_published(src, METRIC, system, "ricci_scalar", [x, "v"], names, functions=functions)
        grr = cv.nr_published(src, METRIC, system, "inverse_metric_components", [x, "v"], names, (x, x), functions)
        values = list(spec["params"].values())
        x_axis = cv.axis(f"${x}/\\ell$", 0.0, 2.0)
        y_axis = cv.axis("$cv/\\ell$", -1.0, 1.0)
        cv.shape_checks(ck, name, [(x_axis["to"] - x_axis["from"]) / (y_axis["to"] - y_axis["from"])])
        xs, ys, cells = cv.cells(K, x_axis, y_axis, n, n, values)
        X, Y = np.meshgrid(xs, ys)
        m = spec["mass"] * (1 + np.tanh(100 * Y)) / 2
        ck(f"{name} K at every cell is anti-de Sitter's plus the hole's term at that moment's m",
           np.allclose(cells, spec["ads"] + spec["hole"](X, m), rtol=1e-12, atol=0))
        ck(f"{name} K is anti-de Sitter's {spec['ads']:g}/l^4 before the shell, v < -0.2 l/c",
           np.allclose(cells[Y < -0.2], spec["ads"], rtol=1e-12, atol=0))
        ck(f"{name} the published Ricci scalar is the constant {spec['ricci']:g}/l^2",
           np.allclose(R(X, Y, *values), spec["ricci"], rtol=1e-14, atol=0))
        ck(f"{name} every cell is finite and inside the bottom of the colours",
           np.all(np.isfinite(cells)) and np.all(cells >= CURVATURE_SCALE["from"]))
        rng = np.random.default_rng(1994)
        picks = sorted({(int(i), int(j)) for i, j in zip(rng.integers(0, n, 40), rng.integers(0, n, 40))}
                       | {(0, n - 1), (0, n // 2), (n - 1, n - 1)})
        cv.agrees(ck, f"{name} map", K, [[xs[i], ys[j], *values] for i, j in picks], [cells[j, i] for i, j in picks])
        # Where K passes the top of its colours after the shell, at the mass M.
        top = CURVATURE_SCALE["to"]
        if x == "r":
            power = 6 if spec["ads"] == 24.0 else 8
            coefficient = spec["hole"](1.0, spec["mass"])
            cut = (coefficient / (top - spec["ads"])) ** (1 / power)
        else:
            cut = ((top - spec["ads"]) / spec["hole"](1.0, spec["mass"])) ** (1 / 6)
        ck(f"{name} K reaches the top of its colours at the {x} stated, at v = l/2c",
           abs(float(K(cut, 0.5, *values)) - top) < 1e-6 * top, f"{cut:.6f}")
        if x == "z":
            ck(f"{name} K stays inside its colours over the whole map, which ends before z = {cut:.2f} l",
               cut > x_axis["to"] and np.all(cells <= top))
        # The apparent horizon, the zero of the published g^xx at each moment.
        track = np.linspace(y_axis["from"], y_axis["to"], 2001)
        curve = []
        for t in track:
            zeros = cv._g_zeros(lambda r, at: grr(r, at, *values), t, x_axis["to"])
            if len(zeros) == 1:
                curve.append([zeros[0], float(t)])
        late = [p for p in curve if p[1] > 0.2]
        ck(f"{name} the apparent horizon is {x} = l after the shell",
           len(late) > 100 and np.allclose([p[0] for p in late], 1.0, rtol=0, atol=1e-9))
        Xp = cv.PLOT_W * np.array([p[0] for p in curve]) / x_axis["to"]
        Yp = cv.PLOT_W * (1 - (np.array([p[1] for p in curve]) - y_axis["from"]) / (y_axis["to"] - y_axis["from"]))
        curve = [curve[k] for k in cv.thinned(Xp, Yp, cv.TOLERANCE)]
        word = "r" if x == "r" else "z"
        dims = "five" if spec["ads"] == 40.0 else "four"
        out.append({
            "id": system, "label": CURVATURE_LABEL[system], "system": system, "kind": "map", "x": x_axis, "y": y_axis,
            "colour": {"name": "$K\\,\\ell^4$", "scale": CURVATURE_SCALE, "ticks": cv.value_ticks(CURVATURE_SCALE)},
            "map": {"metric": METRIC, "system": system, "scalar": "kretschmann", "x": x, "y": "v",
                    "params": dict(spec["params"]), "functions": functions, "nx": n, "ny": n,
                    "cells": [[cv.rounded(c) for c in row] for row in cells]},
            "marks": [{"class": "horizon", "points": curve}],
            "legend": [["line", "horizon", f"the apparent horizon, ${word} = \\ell$ after the shell"]],
            "caption": [CURVATURE_CAPTION[system]],
            "settings": "$\\ell = G = 1$." if "G" in spec["params"] else "$\\ell = 1$.",
            "input": shell_input(system),
            "clip": ("Logarithmic, $10$ to $10^4$" + (f"; $K$ passes $10^4/\\ell^4$ inside $r = {cut:.2f}\\,\\ell$."
                                                       if x == "r" else ".")),
        })
    return out


CURVATURE_LABEL = {"ingoing_four": "Four dimensions", "flat_front": "Flat front",
                   "poincare_four": "Poincaré, four dimensions", "ingoing_five": "Five dimensions"}
CURVATURE_CAPTION = {
    "ingoing_four": ("The Kretschmann scalar $K$ of a thin shell of light falling into anti-de Sitter space of four "
                     "dimensions, over the plane of $v$ and $r$. It is $24/\\ell^4$ before the shell and "
                     "$24/\\ell^4 + 48G^2M^2/c^4r^6$ after it."),
    "flat_front": ("The Kretschmann scalar $K$ of a thin shell of light with a flat front, over the plane of $v$ and "
                   "$r$. It is $24/\\ell^4$ before the shell and $24/\\ell^4 + 12\\ell^2/r^6$ after it."),
    "poincare_four": ("The Kretschmann scalar $K$ of a thin shell of light with a flat front, over the plane of $v$ "
                      "and $z$. It is $24/\\ell^4$ before the shell and $24/\\ell^4 + 12z^6/\\ell^{10}$ after it."),
    "ingoing_five": ("The Kretschmann scalar $K$ of a thin shell of light falling into anti-de Sitter space of five "
                     "dimensions, over the plane of $v$ and $r$. It is $40/\\ell^4$ before the shell and "
                     "$40/\\ell^4 + 288\\ell^4/r^8$ after it."),
}
