"""Drawings of the expanding dust exterior: m=eta_0=1, eta>0, r>2.

The boundary r=2 is singular; only one exterior spatial sheet is embedded.
"""
import math

import numpy as np


SETTINGS = r"$m=\eta_0=1$, $a=(\eta/\eta_0)^2$."


def null_specs(api):
    """One radial plane in each chart, with every ray restricted to the exterior."""
    D = api.Diagram
    common = dict(params={"m": 1}, fixed=api.EQUATOR, areal=True)
    data = [
        ("conformal", (r"\eta", "r"), (1, 25, 0, 24), "eta", {"a": "eta**2"}, "r-2", "r-2"),
        ("cosmic", ("t", "r"), (1, 9, 0, 8), "t", {"a": "(3*t)**Rational(2,3)"}, "r-2", "r-2"),
        ("areal", ("t", "R"), (0, 32, 0, 32), "t", {"a": "(3*t)**Rational(2,3)"},
         "R-2*(3*t)**Rational(2,3)", "R-2*(3*t)**Rational(2,3)"),
        ("isotropic", (r"\eta", r"\sigma"), (0, 24, 0, 24), "eta", {"a": "eta**2"}, "sigma-Rational(1,2)",
         "sigma-Rational(1,2)"),
        ("tortoise", (r"\eta", "x"), (-12, 12, 0, 24), "eta", {"a": "eta**2"}, None, None),
    ]
    rows, captions = [], {}
    for system, plane, box, tau, functions, where, singular in data:
        rows.append(D("thakurta", system, "radial", "Radial plane", plane, box,
                      rf"${plane[1]}/m$", rf"${plane[0]}/m$", tau=tau, functions=functions,
                      where=where, singular_zero=singular, inside=system != "areal", input=SETTINGS, **common))
        captions[("thakurta", system, "radial")] = [
            rf"The plane of ${plane[0]}$ and ${plane[1]}$ ($\theta=\pi/2$, $\phi=0$), in the expanding exterior."]
    rows.append(D("thakurta", "kruskal", "radial", "Null plane", ("U", "V"), (0, 8, 0, 8),
                  "$(V-U)/2$", "$(V+U)/2$", {"m": 1, "eta_0": 1}, api.EQUATOR,
                  to_display=((-0.5, 0.5), (0.5, 0.5)), tau="2*log(-V/U)",
                  where="V+U", areal=True, null_radius=True, other_trapped=True, input=SETTINGS))
    captions[("thakurta", "kruskal", "radial")] = [
        r"The plane of $U$ and $V$ ($\theta=\pi/2$, $\phi=0$), in the expanding exterior $U<0$, $V>-U$."]
    return rows, captions


def embedding(ck, src, api):
    """The exterior Flamm profile at eta=1, with a strict inner cutoff r=2.05."""
    sl = api.Slice(src, "thakurta", "conformal", "r", r"\phi", {"eta": 1, **api.EQUATOR},
                   {"m": 1}, functions={"a": "eta**2"})
    lo, hi, size = 2.05, 8.0, 16.0
    # The height is zero at the cutoff, not at the excluded singular sphere.
    piece = api.Piece("exterior", "sheet", sl, lo, hi, 0.0, 1,
                      (("edge", r"the cutoff $r=2.05m$"), ("edge", r"the exterior continuing to $r\to\infty$")),
                      [(r, "r", None) for r in (lo, 3, 4, 6, hi)], size)
    ck.isometry("Thakurta exterior at eta=1", piece)
    ck.form("Thakurta one exterior Flamm profile", piece,
            lambda r: 2 * np.sqrt(2) * (np.sqrt(r - 2) - np.sqrt(lo - 2)), size)
    surface = api.Surface([piece], label=r"$\eta=\eta_0$", time=1.0)
    fig = api.figure_of([surface], {"sheet": "cover"}, size)
    fig.legend("fill", "cover", r"the equator at $\eta=\eta_0$")
    fig.legend("line", "r", r"constant $r$")
    fig.legend("line", "meridian", r"constant $\phi$")
    return [api.view("exterior", "Exterior spatial slice", "$m$", [surface], fig.done(), settings=SETTINGS)]


def slice_lines(moment, system):
    """The eta=1 embedding slice, transformed into each chart's own coordinates."""
    lo, hi = moment.reach("conformal", "r")
    rr = np.linspace(lo, hi, 801)
    eta = moment.time
    if system in ("conformal", "cosmic", "areal"):
        time = eta if system == "conformal" else eta ** 3 / 3
        return [[(time, float(r * eta ** 2 if system == "areal" else r)) for r in rr]]
    if system == "isotropic":
        return [[(eta, float((r - 1 + np.sqrt(r * (r - 2))) / 2)) for r in rr]]
    x = rr + 2 * np.log(rr / 2 - 1)
    if system == "tortoise":
        return [[(eta, float(v)) for v in x]]
    return [[(-math.exp(-(eta - v) / 4), math.exp((eta + v) / 4)) for v in x]]


def conformal(ck, src, api):
    """The eta>0 exterior is the upper half of the conformal Schwarzschild exterior diamond.

    p=atan((eta-x)/4), q=atan((eta+x)/4), x=r+2 ln(r/2-1).
    Its left future null edge is the curvature singularity approached at fixed eta+x.
    The marginal tube eta=2r^2/(r-2) has both ends at future timelike infinity.
    """
    plane = api.Plane(src, "thakurta", "conformal", (r"\eta", "r"), api.EQUATOR,
                      {"m": 1}, functions={"a": "eta**2"})
    def compact(eta, r):
        x = np.asarray(r) + 2 * np.log(np.asarray(r) / 2 - 1)
        return np.arctan((np.asarray(eta) - x) / 4), np.arctan((np.asarray(eta) + x) / 4)
    eta, rr = ck.uniform(0.1, 40), ck.uniform(2.01, 25)
    ck.chart("Thakurta dust exterior", plane, compact, eta, rr, lambda eta, r: (1, 0))
    K = plane.kretschmann
    ck.diverges("Thakurta initial curvature singularity", K(np.full(3, 0.01), np.array([3, 4, 8])),
                K(np.full(3, 0.001), np.array([3, 4, 8])))
    def at_null(epsilon):
        r = 2 + epsilon
        eta = -(r + 2 * np.log(r / 2 - 1))
        return K(np.array([eta]), np.array([r]))
    ck.diverges("Thakurta future null curvature singularity", at_null(1e-10), at_null(1e-14))
    ck.limit("Thakurta initial boundary has T=0", sum(compact(np.zeros(3), np.array([2.01, 4, 20]))), [0] * 3, 1e-12)
    half = np.pi / 2
    ck.limit("Thakurta future null edge at p=pi/2", compact(1000, 2 + 1e-8)[0], half, 0.005)
    v = api.View("exterior", "Expanding dust exterior", [-np.pi - 0.4, np.pi + 0.4, -0.4, np.pi + 0.4])
    v.fill("cover", [[-np.pi, 0], [np.pi, 0], [0, np.pi]])
    v.segment("singular", (half, -half), (-half, half), zig=True)
    v.segment("singular", (half, -half), (half, half), zig=True)
    v.segment("scri", (-half, half), (half, half))
    r = 2 + np.geomspace(1e-7, 1e5, 1400)
    v.curve("apparent", *compact(2 * r ** 2 / (r - 2), r))
    for radius in (2.1, 3, 4, 8):
        e = api.spread(0, np.inf, 400, 15)
        v.curve("r", *compact(e, np.full_like(e, radius)))
    for moment in api.slices.moments("thakurta"):
        lo, hi = moment.reach("conformal", "r")
        r = np.linspace(lo, hi, 250)
        v.slice(moment, [compact(np.full_like(r, moment.time), r)])
    v.legend("cover", r"$\eta>0$, $r>2m$")
    v.legend("r", r"constant $r$")
    v.legend("apparent", r"past marginal tube")
    v.legend("singular", r"$\eta=0$ and the future null boundary $r=2m$")
    v.legend("scri", r"future null infinity")
    v.label_xt([0, np.pi], "$i^+$", "b", dy=-8)
    v.set(settings=SETTINGS)
    return [v]


def curvature(ck, src, api):
    """Two maps from the stored scalars on the dust exterior, with a signed scale for K."""
    x_axis, y_axis = api.axis("$r/m$", 2.02, 26.02), api.axis(r"$\eta/m$", 0.2, 24.2)
    api.shape_checks(ck, "Thakurta scalar maps", [1.0])
    views = []
    for field, symbol, scale in (
            ("kretschmann", "K", {"kind": "symlog", "from": -1e5, "to": 1e5, "linear": 1e-10}),
            ("ricci_scalar", "R", {"kind": "log", "from": 1e-8, "to": 1e6})):
        published = api.Published(src, "thakurta", "conformal", field, ["r", "eta"], ["m"], functions={"a": "eta**2"})
        nx = ny = 120
        xs, ys, values = api.cells(published, x_axis, y_axis, nx, ny, [1])
        picks = {(i, j) for i in (0, 1, 30, 60, 90, 119) for j in (0, 1, 30, 60, 90, 119)}
        api.agrees(ck, f"Thakurta stored {symbol} map", published,
                   [[xs[i], ys[j], 1] for i, j in picks], [values[j, i] for i, j in picks])
        rr = np.linspace(2.02, 26.02, 1500)
        ee = 2 * rr ** 2 / (rr - 2)
        marked = [[float(r), float(e)] for r, e in zip(rr, ee) if 0.2 <= e <= 24.2]
        views.append({"id": symbol.lower(), "label": rf"${symbol}$ in the exterior", "system": "conformal",
                      "kind": "map", "x": x_axis, "y": y_axis,
                      "colour": {"name": f"${symbol}$", "scale": scale, "ticks": api.value_ticks(scale)},
                      "map": {"metric": "thakurta", "system": "conformal", "scalar": field, "x": "r", "y": r"\eta",
                              "params": {"m": 1}, "functions": {"a": "eta**2"}, "nx": nx, "ny": ny,
                              "cells": [[api.rounded(v) for v in row] for row in values]},
                      "marks": [{"class": "horizon", "points": marked}],
                      "legend": [["line", "horizon", "past marginal tube"]],
                      "caption": [rf"The scalar ${symbol}$ over the expanding exterior ($r\ge2.02m$, $\eta>0$)."],
                      "settings": SETTINGS,
                      "clip": (r"Signed logarithmic colours, linear within $10^{-10}$; values beyond $\pm10^5$ clipped."
                               if symbol == "K" else r"Logarithmic colours from $10^{-8}$ to $10^6$; values beyond the range clipped.")})
    return views
