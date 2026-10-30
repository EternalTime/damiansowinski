"""Radial light sections and the magnetic transverse surface.

All views are stills. A fixed-coordinate section is labelled as a section;
the conformal triangle describes the magnetic radial section alone.
"""
import math

import numpy as np

import early_braneworlds as source

METRIC = source.METRIC
MAGNETIC = "magnetic_flat"


def register_null(nr):
    """One source-specific spacetime section for each chart."""
    for system in source.SYSTEMS:
        d = source.definition(system)
        radial = ("u" if system == "vortex_cartesian" else "w" if system == "scalar_wall" or
                  system.startswith("visser") or system.startswith("squires") else
                  r"\beta" if system in ("magnetic_general", "magnetic_curved") else r"\rho")
        fixed = {nr.vm.Reader._plain(s): "0" for s in d["coords"] if s not in ("t", radial)}
        params, functions, note = {}, {}, None
        box = (-3, 3, -3, 3) if radial in ("u", "w") else (0, 3, -1.5, 1.5)
        if system.startswith("squires"):
            box = (0, 3, -1.5, 1.5)
        if any(p["symbol"] == "L" for p in d["parameters"]):
            params["L"] = 1
        if system == "vacuum_six":
            box = (0, 1.6, -1.1, 1.1)
        if system == "visser_lapse":
            functions["N"] = "cosh(w)"
            note = r"The electric lapse $N=\cosh(w/L)$, with $L=1$."
        if system == "squires_family":
            params.update({f"q_{i}": "3/4" if i == 3 else "-1/4" for i in range(4)})
        if system == MAGNETIC:
            params["a"] = 1
            box = (0, 5, -2.5, 2.5)
        if system in ("magnetic_general", "magnetic_curved"):
            params.update({"C": "3/2", "mu": "4/9", "Q": "2/3", "Lambda": "0"})
            box = (1.05, 3.05, -1.5, 1.5)
            if system == "magnetic_curved":
                params.update({"C": 1, "mu": 2, "Q": "1/2", "Lambda": "-1/10", "K": "1/4"})
                fixed.update({"r": "1/2", "theta": "pi/2"})
        length = "a" if system == MAGNETIC else "L" if "L" in params or system == "visser_lapse" else None
        xlabel = "$" + radial + ("/"+length if length else "") + "$"
        ylabel = "$ct" + ("/"+length if length else "") + "$"
        if radial == r"\beta":
            ylabel = r"$ct\,Q$"
        display = ((0, 1), ((0.5 if system == "magnetic_curved" else 2/3) if radial == r"\beta" else 1, 0))
        row = nr.Diagram(METRIC, system, "radial", "Radial section", ("t", radial), box,
                         xlabel, ylabel, params, fixed, functions=functions, input=note,
                         families=nr.SIDEWAYS, cones=(5, 6), to_display=display)
        key = (METRIC, system, "radial")
        nr.DIAGRAMS.append(row)
        held = ", ".join(s + " = " + (r"\pi/2" if value == "pi/2" else value)
                         for s in d["coords"] if (value := fixed.get(nr.vm.Reader._plain(s))) is not None)
        nr.CAPTIONS[key] = ["The section of $t$ and $"+radial+"$ ($"+held+"$), radial null curves."]


def register_slices(sl):
    """The magnetic embedding meets the magnetic sections at t=0."""
    sl.FLAT_METRICS.add(METRIC)
    sl.FLAT[(METRIC, MAGNETIC, "radial")] = lambda: sl.one(
        METRIC, lambda m: sl.along(0.0, *m.reach(MAGNETIC, r"\rho")), view_id="magnetic")
    sl.FLAT[(METRIC, "magnetic_general", "radial")] = lambda: sl.one(
        METRIC, lambda m: sl.along(0.0, 1.0, 26**(1/3)), view_id="magnetic")
    for system in source.SYSTEMS:
        if system not in (MAGNETIC, "magnetic_general"):
            sl.HIDDEN[(METRIC, system, "radial")] = "a different localization background from the flat magnetic membrane"


def embedding(ck, src):
    """The source's actual polar angle, with radius rho/H and radial length H^(1/3) d rho."""
    import embedding as eb
    sl = eb.Slice(src, METRIC, MAGNETIC, r"\rho", r"\phi", {"t": 0, "x": 0, "y": 0, "z": 0}, {"a": 1})
    top, size = 5.0, 12.0
    piece = eb.Piece("membrane", "sheet", sl, 0, top, 0, 1,
                     (("axis", r"$\rho=0$"), ("edge", r"$\rho=5a$")),
                     [(1, "surface", None), (3, "r", None), (5, "r", None)], size)
    ck.isometry("six dimensional magnetic membrane", piece)
    ck.radius("magnetic circles rho/H", piece, lambda r: r/(1+r*r), size)
    surface = eb.Surface([piece])
    fig = eb.figure_of([surface], {"sheet": "cover"}, size)
    fig.legend("fill", "cover", "magnetic transverse surface")
    fig.legend("line", "surface", r"$\rho=a$")
    fig.legend("line", "r", r"$\rho/a=3,5$")
    fig.legend("line", "meridian", r"$\phi$ constant")
    return [eb.view("magnetic", "Magnetic transverse surface", "$a$", [surface], fig.done(), settings="$a=1$.")]


def register_embedding(eb):
    eb.DRAWN[METRIC] = embedding
    eb.CAPTIONS[(METRIC, "magnetic")] = [
        r"The transverse surface of $\rho$ and $\phi$ ($t=x=y=z=0$), with its widest circle at $\rho=a$."]


def conformal(ck, src):
    """The complete magnetic radial section, conformal to the Minkowski half plane.

    ds^2=H^(2/3)(-dt^2+d rho^2) on this totally geodesic section.
    The conformal factor is positive and the radial null affine length diverges
    toward infinity, so the far null edges are infinities of this section.
    """
    import conformal as cd
    plane = cd.Plane(src, METRIC, MAGNETIC, ("t", r"\rho"),
                     {"x": "0", "y": "0", "z": "0", "phi": "0"}, {"a": 1})
    ck.chart("magnetic radial section", plane, cd.mink_pq,
             ck.uniform(-12, 12, 120), ck.uniform(0.01, 12, 120), lambda t, r: (1, 0))
    v = cd.View("magnetic", "Magnetic radial section",
                [-0.4, math.pi+0.4, -math.pi-0.3, math.pi+0.3], MAGNETIC)
    v.fill("region", cd.TRIANGLE)
    v.fill("cover", cd.TRIANGLE)
    cd.grid(v, "r", lambda r, t: cd.mink_pq(t, r), (0.5, 1, 2, 5), cd.S_ALL)
    cd.grid(v, "t", cd.mink_pq, (-4, -2, 0, 2, 4), cd.S_POS)
    cd.triangle_edges(v, r"$\rho=0$")
    v.legend("cover", "magnetic radial section")
    v.legend("r", r"$\rho$ constant")
    v.legend("t", "$t$ constant")
    v.legend("centre", r"$\rho=0$")
    moment = cd.slices.moments(METRIC, "magnetic")[0]
    r = np.linspace(*moment.reach(MAGNETIC, r"\rho"), 2)
    v.slice(moment, [cd.mink_pq(0*r, r)])
    v.set(restriction=r"The totally geodesic section $x=y=z=0$, $\phi=0$.",
          settings=r"$a=1$, $p,q=\arctan((ct\mp\rho)/a)$.")
    return [v]


def register_conformal(cd):
    cd.DRAWN[METRIC] = conformal
    cd.CAPTIONS[(METRIC, "magnetic")] = [
        r"The magnetic radial section, with its regular axis at $\rho=0$ and radial null infinity on the sloping edges."]


def curvature(ck, src):
    """One line view from the stored magnetic Kretschmann scalar alone."""
    import curvature as cv
    xaxis = cv.axis(r"$\rho/a$", 0, 5)
    scale = {"kind": "log", "from": 1e-7, "to": 1e3}
    value = cv.Published(src, METRIC, MAGNETIC, "kretschmann", ["rho"], ["a"])
    curves = cv.line_curves(ck, METRIC, "K", "K", value, "rho", xaxis, scale, cv.PLOT_W, {"a": 1})
    return [{"id": "magnetic", "label": "Magnetic curvature", "system": MAGNETIC, "kind": "line",
             "x": xaxis, "panels": [{"id": "K", "name": "$K a^4$", "aspect": 1,
                                      "scale": scale, "ticks": cv.value_ticks(scale)}],
             "curves": curves, "marks": [], "legend": [["line", "K", "$K$"]],
             "caption": [r"The Kretschmann scalar $K$ along the magnetic transverse radius $\rho$."],
             "settings": "$a=1$.", "clip": "Logarithmic, none clipped."}]


def register_curvature(cv):
    cv.DRAWN[METRIC] = curvature
