"""The drawings of the Koutras-McIntosh spacetimes.

The two planes of u and the ray coordinate draw Wils's metric, a = 1 and b = 0 in Koutras and
McIntosh's chart and g = h = 0 in Edgar and Ludwig's, under a pulse of radiation switched on and off,
f = e^(-4u^2)/(2 l) with l any length, a sandwich wave of the kind Griffiths and Podolsky's sec. 5
treats, with f >= 0 so that the energy density is positive. The Kerr-Schild chart draws the wave fronts
round their envelope with the figure of Kundt's waves, whose fronts and envelope are the same, Griffiths
and Podolsky's eqs. (8)-(10). koutras_mcintosh.md beside this file is the derivation.
"""

METRIC = "koutras_mcintosh"
PULSE = "exp(-4*u**2)/2"
INPUT = ("A pulse of radiation in Wils's metric, $f = e^{-4u^2}/2\\ell$ with $\\ell$ any length, "
         "positive, so that the energy density is positive.")
CREST = (("shell", "x0", "0", "the crest of the pulse, $u = 0$"),)


def register_null(nr):
    nr.DIAGRAMS.extend([
        nr.Diagram(METRIC, "koutras_mcintosh", "front", "$x = \\ell$, $y = 0$", ("u", "w"), (-2, 2, -2, 2),
                   "$(w/\\ell - u)/2$", "$(u + w/\\ell)/2$", {"a": "1", "b": "0"}, {"x": "1", "y": "0"},
                   to_display=nr.NULL_TO_TR, tau="w + 30*u", families=nr.SIDEWAYS,
                   functions={"f": PULSE}, lines=CREST, input=INPUT),
        nr.Diagram(METRIC, "edgar_ludwig", "front", "$x = \\ell$, $y = 0$", ("u", "v"), (-2, 2, -2, 2),
                   "$(v/\\ell^2 - u)/2$", "$(u + v/\\ell^2)/2$", {}, {"x": "1", "y": "0"},
                   to_display=nr.NULL_TO_TR, tau="v + 30*u", families=nr.SIDEWAYS,
                   functions={"f": PULSE, "g": "0", "h": "0"}, lines=CREST, input=INPUT),
    ])
    nr.CAPTIONS.update({
        (METRIC, "koutras_mcintosh", "front"): [
            "The plane of $u$ and $w$ ($x = \\ell$, $y = 0$) of Wils's metric ($a = 1$, $b = 0$) under a pulse "
            "of radiation, with $(w/\\ell - u)/2$ across. The rays moving right keep their $u$; a curve moving "
            "left obeys $dw/du = (w^2 - 2fx^3)/2x$."],
        (METRIC, "edgar_ludwig", "front"): [
            "The plane of $u$ and $v$ ($x = \\ell$, $y = 0$) of Wils's metric ($g = h = 0$) under the same pulse, "
            "with $(v/\\ell^2 - u)/2$ across. The rays moving right keep their $u$; a curve moving left obeys "
            "$dv/du = (v^2/x^2 - 2fx^3)/2$."],
    })


def register_projection(pj):
    pj.FIGURES.append(pj.Projection(METRIC, "kerr_schild", "fronts", "the wave fronts round their envelope",
                                    lambda spec: pj.kundt_fronts(spec), {"ell": 1}, {"Y": "0"}))
    pj.CAPTIONS[(METRIC, "kerr_schild", "fronts")] = [
        "The slice $Y = 0$ of Wils's metric ($f = 1/\\ell$), $T$ up: the envelope $X^2 + Z^2 = c^2T^2$ at "
        "$cT = \\ell$ and $2\\ell$, and six wave fronts touching it."]
