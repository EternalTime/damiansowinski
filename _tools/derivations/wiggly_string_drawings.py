"""The weak-field annulus of the averaged string, with one registered moment."""

METRIC = "wiggly_string"
SYSTEM = "cylindrical"
PARAMS = {"U": 0.01, "T": 0.0025, "r_0": 1}
LOW, HIGH = 0.5, 3.5
SETTINGS = r"$U=1/100$, $T=1/400$, $r_0=1$."


def register_null(nr):
    nr.DIAGRAMS.append(nr.Diagram(
        METRIC, SYSTEM, "radial", "Radial plane", ("t", "r"),
        (0.25, 3.75, -1.75, 1.75), "$r/r_0$", "$t/r_0$", PARAMS,
        {"phi": "0", "z": "0"}, cones=(5, 5)))
    nr.CAPTIONS[(METRIC, SYSTEM, "radial")] = [
        r"Radial null geodesics of the averaged string ($\phi=z=0$)."]


def embedding_views(ck, src, em):
    import numpy as np
    sl = em.Slice(src, METRIC, SYSTEM, "r", r"\phi", {"t": 0, "z": 0}, PARAMS)
    piece = em.Piece("annulus", "sheet", sl, LOW, HIGH, 0, 1,
                     (("edge", r"$r/r_0=1/2$"), ("edge", r"$r/r_0=7/2$")),
                     [(r, "r", None) for r in (LOW, 1, 2, 3, HIGH)], 7.0)
    ck.isometry("wiggly string: the transverse annulus", piece)
    ck.radius("wiggly string: circumferential radius", piece,
              lambda r: r*np.sqrt(1-4*(PARAMS["U"]+PARAMS["T"])*np.log(r)), 7.0)
    # The derived defect B-(d(r sqrt(B))/dr)^2=4(U+T)(1-(U+T)/B) is positive here.
    radii = np.linspace(LOW, HIGH, 101)
    B = 1-4*(PARAMS["U"]+PARAMS["T"])*np.log(radii)
    expected = 4*(PARAMS["U"]+PARAMS["T"])*(1-(PARAMS["U"]+PARAMS["T"])/B)
    ck.add("wiggly string: Euclidean embedding defect", float(np.max(np.abs(sl.defect_at(radii)-expected))), 1e-12)
    surface = em.Surface([piece])
    fig = em.figure_of([surface], {"sheet": "cover"}, 7.0)
    fig.legend("fill", "cover", r"$t=z=0$")
    fig.legend("line", "r", "$r$ constant")
    fig.legend("line", "meridian", r"$\phi$ constant")
    data = fig.done()
    data["box"] = [-4, 4, -4, 4]
    return [em.view("annulus", "Transverse annulus", "$r_0$", [surface], data,
                    settings=SETTINGS,
                    stops=[r"The annulus is restricted to $1/2\le r/r_0\le7/2$."])]


def register_embedding(em):
    em.DRAWN[METRIC] = lambda ck, src: embedding_views(ck, src, em)
    em.CAPTIONS[(METRIC, "annulus")] = [
        r"The transverse spatial annulus of the averaged string ($t=z=0$)."]


def slice_marks(slices):
    return slices.one(METRIC, lambda moment: slices.along(0, *moment.reach(SYSTEM, "r")),
                      label="$t=0$", view_id="annulus")


def curvature_views(ck, src, cv):
    """Read the stored exact scalars of the truncated exterior; recompute no tensor."""
    xaxis = cv.axis("$r/r_0$", LOW, HIGH)
    scales = {"K": {"kind": "log", "from": 1e-5, "to": 0.1},
              "R": {"kind": "log", "from": 1e-4, "to": 0.02}}
    panels = [{"id": field, "name": "$"+field+("r_0^4$" if field == "K" else "r_0^2$"),
               "aspect": 2, "scale": scales[field], "ticks": cv.value_ticks(scales[field])}
              for field in ("K", "R")]
    cv.shape_checks(ck, METRIC, [p["aspect"] for p in panels])
    curves = []
    for field, scalar in (("K", "kretschmann"), ("R", "ricci_scalar")):
        stored = cv.nr_published(src, METRIC, SYSTEM, scalar, ["r"], list(PARAMS))
        curves += cv.line_curves(ck, METRIC, field, field, stored, "r", xaxis,
                                 scales[field], cv.PLOT_W/2, PARAMS)
    return [{"id": "annulus", "label": "Radial curvature", "kind": "line", "system": SYSTEM,
             "x": xaxis, "panels": panels, "curves": curves, "marks": [],
             "legend": [["line", "K", "$K$"], ["line", "R", "$R$"]],
             "caption": [r"The Kretschmann scalar $K$ and the Ricci scalar $R$ of the averaged string across its transverse annulus ($1/2\le r/r_0\le7/2$)."],
             "ricci_sentences": [], "settings": SETTINGS, "clip": "Logarithmic vertical scales.",
             "per_scalar": cv.scalar_clips("Logarithmic scale.", "Logarithmic scale.")}]


def register_curvature(cv):
    cv.DRAWN[METRIC] = lambda ck, src: curvature_views(ck, src, cv)
