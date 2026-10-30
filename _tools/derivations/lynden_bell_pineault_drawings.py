"""Restricted light sections and conical spatial sections of the relativistic discs."""
import math

import numpy as np

import lynden_bell_pineault as lbp

METRIC = lbp.METRIC
ANGLE = math.pi/3
REACH = (0.3, 3.5)


def angular_values(theta=ANGLE):
    """Restrictions of the integrated angular functions to one fixed polar angle."""
    solution = lbp.rotating_solution()
    N, _, Q, _, Z = solution["solution"].sol(theta)
    return solution["n"], {"N": float(N), "P": math.log(math.sin(theta)), "Q": float(Q), "Z": float(Z)}


def register_null(nr):
    """Four declared sections, with closed null coordinates checked by null_rays.py."""
    nrot, angular = angular_values()
    fns = {"nu": "n*log((sqrt(R**2+z**2)+z)/(2*A))",
           "xi": "2*n**2*log((1+z/sqrt(R**2+z**2))/2)"}
    rows = [nr.Diagram(METRIC, "weyl", "axis", "The upper axis", ("t", "z"), (0.3, 3.5, -1.6, 1.6),
                       "$z/A$", "$ct/A$", {"n": "1/5", "A": 1}, {"R": "0", "phi": "0"},
                       functions=fns, kretschmann=False, input="The infinite counter-rotating member.")]
    captions = {"weyl": ["The axial plane of $t$ and $z$ ($R = 0$, $\\phi = 0$), null geodesics above the disc."]}
    forms = {"weyl": lambda z: np.asarray(z)**0.6/0.6}
    for system in ("spherical", "photons", "rotating"):
        params = {"A": 1}
        functions = {}
        note = None
        if system == "spherical":
            params["n"] = "1/5"
            factor = (3/4)**(2/25-2/5)
            forms[system] = lambda r, a=factor: a*np.asarray(r)**0.6/0.6
        elif system == "photons":
            forms[system] = lambda r: np.log(np.asarray(r))/math.sqrt(3/4)
        else:
            params["n"] = repr(nrot)
            functions = {k: repr(v) for k, v in angular.items()}
            factor = math.exp(angular["Z"]/2-angular["N"])/math.sqrt(1-angular["Q"]**2)
            forms[system] = lambda r, a=factor, n=nrot: a*np.asarray(r)**(1-2*n)/(1-2*n)
            note = "The cold solution at $V = 0.2c$, restricted to $\\theta = \\pi/3$."
        rows.append(nr.Diagram(METRIC, system, "radial", "A radial vacuum section", ("t", "r"),
                               (0.3, 3.5, -1.6, 1.6), "$r/A$", "$ct/A$", params,
                               {"theta": "pi/3", "phi": "0"}, functions=functions,
                               kretschmann=system != "rotating", input=note))
        captions[system] = ["The plane of $t$ and $r$ ($\\theta = \\pi/3$, $\\phi = 0$), radial null curves in the vacuum exterior."]
    for row in rows:
        place = (METRIC, row.system, row.view)
        nr.DIAGRAMS.append(row)
        nr.CAPTIONS[place] = captions[row.system]
        star = forms[row.system]
        nr.CLOSED_FORMS[place] = (lambda t, r, f=star: t+f(r), lambda t, r, f=star: t-f(r), None)


def embedding(ck, src):
    """Cones at fixed theta=pi/3, read from the induced metric, at one constant time.

    These are vacuum sections above the matter sheet, not surfaces made of disc matter.
    Their semiangle satisfies sin(alpha)=(1-n)sin(theta)exp(-Z/2) in the rotating chart.
    """
    import embedding as eb
    views = []
    nrot, angular = angular_values()
    for system in ("spherical", "photons", "rotating"):
        params = {"A": 1}
        fns = {}
        if system == "rotating":
            n = nrot
            params["n"] = repr(n)
            fns = {k: repr(v) for k, v in angular.items()}
            radial = math.exp((angular["Z"]-angular["N"])/2)
            radius = math.exp(angular["P"]-angular["N"]/2)
        else:
            n = 0.5 if system == "photons" else 0.2
            if system == "spherical":
                params["n"] = "1/5"
            h = 3/4
            radial = h**(2*n*n-n)
            radius = math.sin(ANGLE)*h**(-n)
        sl = eb.Slice(src, METRIC, system, "r", "\\phi", {"t": 0, "theta": "pi/3"}, params, functions=fns)
        lo, hi = REACH
        size = 2*radius*hi**(1-n)
        piece = eb.Piece("cone", "sheet", sl, lo, hi, 0.0, 1,
                         (("edge", "$r = 0.3A$"), ("edge", "$r = 3.5A$")),
                         [(1.0, "r", None), (2.0, "r", None), (hi, "r", None)], size)
        ck.isometry(f"{METRIC} {system}", piece)
        alpha = (1-n)*radius/radial
        ck.form(f"{METRIC} {system}: conical profile", piece,
                lambda r, a=alpha, c=radius, n=n: c*math.sqrt(1/a**2-1)*(np.asarray(r)**(1-n)-lo**(1-n)), size)
        surface = eb.Surface([piece])
        fig = eb.figure_of([surface], {"sheet": "cover"}, size, eb.Camera(-90, 28))
        fig.legend("fill", "cover", "vacuum cone")
        fig.legend("line", "r", "$r/A = 1,2,3.5$")
        fig.legend("line", "meridian", "$\\phi$ constant")
        views.append(eb.view(system, system.capitalize(), "$A$", [surface], fig.done(),
                             settings=("$V = 0.2c$, $A = 1$." if system == "rotating" else
                                       "$n = 1/2$, $A = 1$." if system == "photons" else "$n = 1/5$, $A = 1$."),
                             input=("The cold rotating solution at $\\theta = \\pi/3$." if fns else None)))
    return views


def register_embedding(eb):
    eb.DRAWN[METRIC] = embedding
    for system in ("spherical", "photons", "rotating"):
        eb.CAPTIONS[(METRIC, system)] = ["The spatial cone of $r$ and $\\phi$ ($t = 0$, $\\theta = \\pi/3$), in the vacuum above the matter sheet."]


def register_slices(sl):
    for system in ("spherical", "photons", "rotating"):
        sl.FLAT[(METRIC, system, "radial")] = lambda s=system: sl.one(
            METRIC, lambda m, s=s: sl.along(0.0, *m.reach(s, "r")), view_id=s)
    sl.FLAT_METRICS.add(METRIC)
    sl.HIDDEN[(METRIC, "weyl", "axis")] = "the embedded cones at theta=pi/3 do not meet the upper axis"
