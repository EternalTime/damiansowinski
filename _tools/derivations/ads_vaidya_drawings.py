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


def register_null(nr):
    for system in STEP:
        function, words = STEP[system]
        shell = "A thin shell of light of mass $M$" if system == "ingoing_four" else "A thin shell of light"
        common = dict(params={"ell": 1, **({"G": 1} if system == "ingoing_four" else {})}, fixed=FIXED[system],
                      functions={"m": function}, lines=(("shell", "x0", "0", SHELL),),
                      # The horizon is seeded on r = l well after the shell, where the mass is the final one
                      # to rounding, and traced both ways.
                      marked=(("event", {"x0": "2/5", "r": "1"}, 0 if system in POINCARE else 1, "the event horizon"),),
                      input=shell + " falling in along $v = 0$, with " + words + " and $cv_0 = \\ell/100$.")
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
        (METRIC, "flat_front", "shell"): (lambda v, r: v, lambda v, r: 1 / r + v / 2, before),
        **{(METRIC, system, "shell"): (lambda v, z: z + v / 2, lambda v, z: v, before) for system in POINCARE},
    })
