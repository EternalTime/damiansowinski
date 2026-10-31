#!/usr/bin/env python3
"""Compute, for every published chart, the facts its metric decides, which the tags rest on.

For each coordinate system in MFS/assets/data/metrics the script reads the line element as
verify_metrics.py does, and computes with the same Geometry:

- the dimension, the number of coordinates;
- whether the Ricci tensor vanishes, which is vacuum;
- whether the Ricci tensor is a constant multiple of the metric, R_{ab} = k g_{ab} with
  k not zero, which is an Einstein space, the vacuum of a cosmological constant
  Lambda = (n - 2) k / 2, and the sign of k where the parameters decide it;
- whether the chart is conformally flat: the Weyl tensor vanishes, or in three dimensions,
  where the Weyl tensor vanishes identically, the Cotton tensor does, and in two always;
- every coordinate the metric does not depend on, whose coordinate vector is therefore a
  Killing vector, with three things about it: whether it is timelike somewhere in the
  chart's domain; whether it is orthogonal to a family of hypersurfaces,
  xi_{[a} d_b xi_{c]} = 0; and whether it is a translation of the chart, which is to say
  that the coordinate runs over the whole real line and no other entry of the chart's
  domains names it. A timelike translation is what makes a chart stationary, and one that
  is hypersurface orthogonal as well makes it static. The inertial chart of the Milne
  universe has a timelike Killing vector in d/dT, but T runs from 0 and the chart ends on
  the cone R = cT, so a translation in T carries the chart off itself, and it is no
  translation of the chart;
- whether the chart shows a round sphere: angles a_1 ... a_k, k >= 2, entering only through
  F (da_1^2 + sin^2 a_1 da_2^2 + ...), with F and every other component free of them, the
  last angle running over [0, 2 pi) and every other over [0, pi], so that a sphere with a
  wedge missing or a piece of one does not count.

A fact is a statement about one chart. decided_tags() below turns the facts of a
spacetime's charts into its tags: REGIONS says which charts count and how they divide into
regions, and OVERRULED names the few tags the charts cannot decide.

The answers are written to metric_tags.json beside this file, each under a stamp of what it
was computed from: the coordinates, their domains, the line element and the parameters.
The test suite has no sympy, so it holds every published chart to a stamp in that file and
every spacetime's tags to the facts under it; a chart that is new or whose line element
changed fails there until this script has been run again. A chart whose stamp still matches
is not computed again, so a run after one new spacetime takes seconds.

    .venv.noindex/bin/python _tools/derivations/metric_tags.py
    .venv.noindex/bin/python _tools/derivations/metric_tags.py --system kerr/boyer_lindquist
    .venv.noindex/bin/python _tools/derivations/metric_tags.py --all     # every chart again
    .venv.noindex/bin/python _tools/derivations/metric_tags.py --check   # write nothing
"""
import argparse
import hashlib
import json
import random
import re
import sys
import time
import unicodedata
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
METRICS_DIR = ROOT / "MFS" / "assets" / "data" / "metrics"
FACTS = HERE / "metric_tags.json"

# How many points a Killing vector's norm is sampled at, and the seed they are drawn with,
# so that every run samples the same points.
SAMPLES = 400
SEED = 20261002


def stamp(entry):
    """What a chart's facts were computed from. Needs no sympy, so the tests share it."""
    source = {
        "coords": entry["coords"],
        "domains": entry.get("domains", []),
        "line_element": entry["line_element"],
        "parameters": [p["symbol"] for p in entry.get("parameters", [])],
    }
    text = json.dumps(source, sort_keys=True, ensure_ascii=False)
    return hashlib.sha256(text.encode("utf-8")).hexdigest()[:16]


def published_charts():
    """Every published chart as (metric id, entry), in the order of the files."""
    for path in sorted(METRICS_DIR.glob("*.json")):
        metric = json.loads(path.read_text(encoding="utf-8"))
        for entry in metric.get("coordinates", []):
            yield metric["id"], entry


def load_facts():
    if not FACTS.exists():
        return {}
    return json.loads(FACTS.read_text(encoding="utf-8"))


def write_facts(facts):
    ordered = {key: facts[key] for key in sorted(facts)}
    FACTS.write_text(json.dumps(ordered, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")


# ---------------------------------------------------------------------------------------
# From the facts of a spacetime's charts to its tags.
# ---------------------------------------------------------------------------------------

DIMENSION = {2: "two-dimensional", 3: "three-dimensional", 4: "four-dimensional",
             5: "five-dimensional", 6: "six-dimensional", 10: "ten-dimensional", 11: "eleven-dimensional"}

# The tags the facts decide. A spacetime carries each of them exactly when its charts say
# so, which the tests hold every metric file to, so none of them is written by hand.
OWNED = ("vacuum", "Einstein space", "conformally flat", "flat", "stationary", "static",
         "spherically symmetric") + tuple(DIMENSION.values())

# A tag that was merged into another, and the one it was merged into, so that the pair
# cannot come back as two filters that each show part of one set. None is a tag dropped
# outright: "non-vacuum" said only that "vacuum" is absent, and a search for "vacuum"
# found it.
RETIRED = {
    "vacuum solution": "vacuum",
    "non-vacuum": None,
    "cosmological": "cosmology",
    "dynamical": "time-dependent",
    "singularity": "curvature singularity",
    "gravitational waves": "gravitational wave",
    "gravitational radiation": "gravitational wave",
    "axial symmetry": "axisymmetric",
    "causality violation": "closed timelike curves",
    "Gregory-Laflamme": "Gregory-Laflamme instability",
    "collapse": "gravitational collapse",
    "interior": "stellar interior",
    "extra dimensions": "higher dimensions",
    "nonsingular": "regular black hole",
    "nonsingular black hole": "regular black hole",
    "acceleration": "accelerating",
    "plane-fronted wave": "pp-wave",
    "pure radiation": "null dust",
    "three-dimensional gravity": "2+1 dimensional gravity",
    # Merged when the whole collection was revisited on 3 October 2026.
    "alternative theory of gravity": "alternative theories of gravity",
    "gravitomagnetic": "gravitomagnetism",
    "Lense-Thirring": "frame dragging",
    "Painlevé-Gullstrand": "Painlevé-Gullstrand coordinates",
    "point masses": "point mass",
    "point particle": "point mass",
    "NUT charge": "NUT parameter",
    "nut": "nuts and bolts",
    "bolt": "nuts and bolts",
    "extremal black hole": "extremal",
    "extremal horizon": "extremal",
    "degenerate horizon": "extremal",
    "circular symmetry": "circularly symmetric",
    "plane symmetric": "planar symmetry",
    "Bekenstein-Hawking entropy": "black hole entropy",
    "binary black hole": "black hole collision",
    "colliding black holes": "black hole collision",
    "black hole merger": "black hole collision",
    "periastron advance": "perihelion precession",
    "tachyonic momentum": "tachyon",
    "accelerating black hole": "accelerating",
    "boost symmetry": "boost",
    "Hawking-Page transition": "phase transition",
    "membrane": "brane",
    "three-brane": "D3-brane",
    "positive mass theorem": "positive energy",
    "Schwarzschild interior": "black hole interior",
    "supersymmetry": "supersymmetric",
    "Kaup limit": "maximum mass",
    "magnetic dipole": "dipole",
    "cosmic dipole": "dipole",
    "Euclidean signature": "Riemannian",
    "instanton": "gravitational instanton",
    "dumb hole": "acoustic black hole",
    "homothetic": "self-similar",
    "ergosphere": "ergoregion",
    "charged": "electric charge",
    "regular": "no singularity",
    "conical": "angular deficit",
    "solid angle deficit": "angular deficit",
}

# A spacetime's charts are taken as charts of one region, each a chart of the whole
# spacetime, unless it is named here. A spacetime named here gives its regions, each a list
# of charts, and a chart it leaves out is not counted: it draws a part or a special case
# that does not speak for the spacetime. A statement about the curvature has to hold in
# every chart counted, and a symmetry has to show in some chart of every region.
REGIONS = {

    "electron_bubble": {"regions": [["exterior", "exterior_kerr"], ["interior", "interior_kerr", "cartesian"]],
                        "why": "the flat interior and electrovacuum exterior are distinct smooth regions"},
