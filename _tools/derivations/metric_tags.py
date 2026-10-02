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

    /tmp/mfs-venv/bin/python _tools/derivations/metric_tags.py
    /tmp/mfs-venv/bin/python _tools/derivations/metric_tags.py --system kerr/boyer_lindquist
    /tmp/mfs-venv/bin/python _tools/derivations/metric_tags.py --all     # every chart again
    /tmp/mfs-venv/bin/python _tools/derivations/metric_tags.py --check   # write nothing
"""
import argparse
import hashlib
import json
import random
import re
import sys
import time
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
             5: "five-dimensional", 6: "six-dimensional"}

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
}

# A spacetime's charts are taken as charts of one region, each a chart of the whole
# spacetime, unless it is named here. A spacetime named here gives its regions, each a list
# of charts, and a chart it leaves out is not counted: it draws a part or a special case
# that does not speak for the spacetime. A statement about the curvature has to hold in
# every chart counted, and a symmetry has to show in some chart of every region.
REGIONS = {
    "black_saturn": {
        "regions": [["weyl", "polar"]],
        "why": "ring is the black ring alone, the Saturn with its hole taken away",
    },
    "cosmic_string": {
        "regions": [["conical"]],
        "why": "interior_cap is Gott's model of the string's own core; the spacetime of "
               "the entry is the cone outside it",
    },
    "gott_time_machine": {
        "regions": [["centre_of_momentum", "grant_rindler", "grant_milne"]],
        "why": "string_rest covers the half space of one string in that string's rest "
               "frame, where it is static; the two halves share no rest frame",
    },
    "hartle_thorne": {
        "regions": [["hartle_thorne", "lense_thirring"]],
        "why": "those two are kept to an order and are a vacuum through it; "
               "painleve_gullstrand is exact for its line element as written, which is "
               "a vacuum only to first order in the spin",
    },
    "kaluza_klein_black_hole": {
        "regions": [["electric", "eddington_finkelstein_ingoing", "magnetic", "dyonic"]],
        "why": "einstein and einstein_eddington_finkelstein print the metric of four dimensions the "
               "holes reduce to, which carries a Maxwell field and a scalar field; the entry's "
               "spacetime is the vacuum of five dimensions",
    },
    "kastor_traschen": {
        "regions": [["cartesian", "comoving"]],
        "why": "isotropic is the single hole, which alone is spherically symmetric, and "
               "cylindrical the holes on one axis",
    },
    "majumdar_papapetrou": {
        "regions": [["cartesian"]],
        "why": "isotropic is the single hole, which alone is spherically symmetric, and "
               "cylindrical the holes on one axis",
    },
    "misner_brill_lindquist": {
        "regions": [["cartesian"]],
        "why": "isotropic is the single hole, which alone is spherically symmetric, cylindrical and "
               "bispherical the holes on one axis, and charged the slice with an electric field on it",
    },
    "oppenheimer_snyder": {
        "regions": [["interior_comoving"], ["exterior_schwarzschild"]],
        "why": "the dust ball and the vacuum outside it are two regions of one spacetime",
    },
    "ori_time_machine": {
        "regions": [["foliation", "brinkmann"]],
        "why": "vacuum_core leaves f free; the core is a vacuum where f is harmonic in x "
               "and y, as the other two charts take it",
    },
    "point_particle_2plus1": {
        "regions": [["conical", "wedge", "circumference", "isotropic", "two_bodies", "moving"]],
        "why": "planet is the inside of a body of finite size, a cap of a sphere; the "
               "entry's spacetime is the flat cone outside a particle",
    },
    "pp_wave": {
        "regions": [["exact_plane_wave"]],
        "why": "brinkmann leaves the profile H free; the wave is a vacuum where H is "
               "harmonic in x and y, as the plane wave's is",
    },
    "spinning_string": {
        "regions": [["proper_radius", "rescaled_radius", "circumference_radius", "helical"]],
        "why": "extended_source is the inside of a string of finite thickness",
    },
    "string_wave": {
        "regions": [["isotropic", "moving_string"]],
        "why": "null_conical leaves the profile F free; the wave is a vacuum where F is "
               "harmonic on the cone, as the other two charts take it",
    },
    "van_den_broeck": {
        "regions": [["cartesian"]],
        "why": "pocket and proper_radial cover the inside of the bubble alone, where "
               "f = 1 and nothing depends on the time",
    },
    "wormhole_time_machine": {
        "regions": [["lorentz"], ["wormhole"]],
        "why": "the flat space outside the mouths and the wormhole whose mouth is moved "
               "are two regions; short_throat is a mouth in its own rest frame",
    },
}

# Where the charts cannot decide a tag and the entry's history and sources do.
OVERRULED = {
    ("coleman_de_luccia", "stationary"): (
        False, "each vacuum is static in its own chart, and the wall between them, at r = r_w, "
               "accelerates outward through both, so the bubble as a whole has no time translation"),
    ("coleman_de_luccia", "static"): (
        False, "each vacuum is static in its own chart, and the wall between them, at r = r_w, "
               "accelerates outward through both, so the bubble as a whole has no time translation"),
    ("black_saturn", "vacuum"): (
        True, "V, Omega, W and nu are left free in Weyl's chart and the polar chart; black Saturn "
              "is the solution of the vacuum equations their parameters define, which "
              "print_charts.py holds to a vanishing Ricci tensor"),
    ("double_kerr", "vacuum"): (
        True, "f, omega and gamma are left free; the two Kerr black holes are a solution of "
              "the vacuum equations, Ernst's equation for f and omega and a quadrature for gamma"),
    ("morgan_morgan", "vacuum"): (
        True, "psi and gamma are left free in Weyl's chart; off the disc the field is the solution "
              "of the vacuum equations, Laplace's equation and a quadrature, that the oblate "
              "spheroidal chart writes out for the first disc, whose Ricci tensor vanishes"),
    ("einstein_rosen_waves", "vacuum"): (
        True, "psi and gamma are left free in both charts; the waves are the solutions of "
              "the vacuum equations, a wave equation for psi and a quadrature for gamma"),
    ("gowdy", "vacuum"): (
        True, "P, Q and lambda are left free in every chart; Gowdy's universes are the "
              "solutions of the vacuum equations, two wave equations and a quadrature"),
    ("mixmaster", "vacuum"): (
        True, "the three scale factors are left free; the Mixmaster universe is the "
              "vacuum of Bianchi type IX, whose equations fix them"),
    ("robinson_trautman", "vacuum"): (
        True, "P and H are left free in both charts; the spacetimes are the vacuum "
              "solutions, on which P obeys the Robinson-Trautman equation"),
    ("spinning_string", "static"): (
        False, "the time translation is orthogonal to the surfaces of constant "
               "t + a phi / c, which do not close up round the string, so the string is "
               "static in every patch and only stationary as a whole"),
    ("lewis", "static"): (
        False, "the canonical chart's time translation is orthogonal to the surfaces of constant "
               "t + 4j phi / (c(1 - 4 sigma)), which do not close up round the axis, so the Weyl "
               "class is static in every patch and only stationary as a whole, and the Lewis "
               "class is static in no patch"),
}


def facts_of(metric, facts):
    """The facts of one metric file's charts, by chart id."""
    return {entry["id"]: facts[f"{metric['id']}/{entry['id']}"]
            for entry in metric.get("coordinates", [])}


def _time_translation(chart, orthogonal=False):
    return any(k["timelike"] and k["translation"] and (k["hypersurface_orthogonal"] or not orthogonal)
               for k in chart["killing"])


def _round(chart):
    return chart["sphere"] is not None and len(chart["sphere"]) >= chart["dimension"] - 2


def decided_tags(metric_id, charts):
    """What the charts of one spacetime decide: (owned, implied).

    `charts` is its facts by chart id. `owned` gives every tag of OWNED as True or False,
    the tags the spacetime has to carry and the ones it must not. `implied` is the set of
    further tags it has to carry and that others may carry by hand: an Einstein space has
    a cosmological constant, and its sign where the parameters decide it.
    """
    regions = REGIONS.get(metric_id, {}).get("regions") or [list(charts)]
    counted = [charts[c] for region in regions for c in region]

    def everywhere(test):
        return all(test(chart) for chart in counted)

    def in_every_region(test):
        return all(any(test(charts[c]) for c in region) for region in regions)

    owned = {
        "vacuum": everywhere(lambda c: c["ricci_flat"]),
        "Einstein space": everywhere(lambda c: c["einstein"] is not None),
        "conformally flat": everywhere(lambda c: c["conformally_flat"]),
        "flat": everywhere(lambda c: c["ricci_flat"] and c["conformally_flat"]),
        "stationary": in_every_region(_time_translation),
        "static": in_every_region(lambda c: _time_translation(c, orthogonal=True)),
        "spherically symmetric": in_every_region(_round),
    }
    dimensions = {chart["dimension"] for chart in charts.values()}
    for n, tag in DIMENSION.items():
        owned[tag] = n in dimensions
    unnamed = dimensions - set(DIMENSION)
    if unnamed:
        raise ValueError(f"{metric_id} has a chart of {sorted(unnamed)} dimensions, which DIMENSION has no tag for")
    for (overruled, tag), (value, _) in OVERRULED.items():
        if overruled == metric_id:
            owned[tag] = value

    implied = set()
    if owned["Einstein space"]:
        implied.add("cosmological constant")
        signs = {chart["einstein"]["sign"] for chart in counted} - {None}
        if signs == {1}:
            implied.add("positive cosmological constant")
        if signs == {-1}:
            implied.add("negative cosmological constant")
    return owned, implied


def tag_faults(metric, facts):
    """What is wrong with one metric file's tags, as sentences; empty if nothing is."""
    tags = metric["tags"]
    faults = [f"carries {tag!r} twice" for tag in sorted(set(tags)) if tags.count(tag) > 1]
    for tag in tags:
        if tag in RETIRED:
            into = RETIRED[tag]
            faults.append(f"carries {tag!r}, which was " + (f"merged into {into!r}" if into else "dropped"))
    owned, implied = decided_tags(metric["id"], facts_of(metric, facts))
    for tag, holds in owned.items():
        if holds and tag not in tags:
            faults.append(f"lacks {tag!r}, which its charts decide it has")
        if not holds and tag in tags:
            faults.append(f"carries {tag!r}, which its charts decide it has not")
    for tag in sorted(implied - set(tags)):
        faults.append(f"lacks {tag!r}, which its charts imply")
    return faults



INTERVAL = re.compile(r"^(?P<open>[\[(])(?P<body>.*)(?P<close>[\])])$")
TOKEN = re.compile(r"\\[A-Za-z]+(?:_(?:\{[^{}]*\}|\\[A-Za-z]+|\w))?|[A-Za-z](?:_(?:\{[^{}]*\}|\\[A-Za-z]+|\w))?")


def _plain_domain(text):
    """A domain entry with its prose and its spacing commands taken out.

    A remark in parentheses, as `\\text{(the circles of constant } t, r, z \\text{ are closed
    timelike curves)}`, goes whole, with the mathematics inside it.
    """
    text = re.sub(r"\\text\{\(.*\)\}", " ", text)
    text = re.sub(r"\\text\{[^{}]*\}", " ", text)
    text = text.replace("\\left", "").replace("\\right", "")
    return re.sub(r"\\[,;!]", " ", text).strip()


def interval(entry, coord):
    """(lower, upper) as written if the entry is `coord \\in [lower, upper)`, else None.

    An entry that goes on past its interval, as FRW's `r \\in [0, \\infty) for k <= 0`, is
    one of several for its coordinate and is not read.
    """
    text = _plain_domain(entry)
    head = coord + " \\in"
    if not text.startswith(head):
        return None
    match = INTERVAL.match(text[len(head):].strip())
    if match is None:
        return None
    body, depth = match.group("body"), 0
    for at, character in enumerate(body):
        depth += character in "([{"
        depth -= character in ")]}"
        if character == "," and depth == 0:
            return body[:at].strip(), body[at + 1:].strip()
    return None


def names(entry):
    """The symbols a domain entry names: every letter and every command, with its subscript."""
    return set(TOKEN.findall(_plain_domain(entry)))


def is_translation(entry, coord):
    """Whether the coordinate runs over the whole line and no other domain entry names it."""
    domains = entry.get("domains", [])
    own = [d for d in domains if _plain_domain(d).startswith(coord + " \\in")]
    if len(own) != 1 or interval(own[0], coord) != ("-\\infty", "\\infty"):
        return False
    return not any(coord in names(d) for d in domains if d is not own[0])


def is_whole_sphere(entry, angles):
    """Whether the last angle runs over [0, 2 pi) and every other over [0, pi]."""
    domains = entry.get("domains", [])

    def runs(coord):
        found = [interval(d, coord) for d in domains if _plain_domain(d).startswith(coord + " \\in")]
        return found[0] if len(found) == 1 else None

    return (all(runs(a) == ("0", "\\pi") for a in angles[:-1])
            and runs(angles[-1]) == ("0", "2\\pi"))


def compute(metric_id, entry, seconds):
    """The facts of one chart. Raises verify_metrics.Timeout if a tensor does not finish."""
    import sympy as sp
    import verify_metrics as vm

    coords = entry["coords"]
    declared = vm.DIMENSIONS[(metric_id, entry["id"])]
    relations = vm.PARAMETER_RELATIONS.get((metric_id, entry["id"]), {})
    orders = vm.ORDERS.get((metric_id, entry["id"]))
    reader = vm.Reader(coords, [p["symbol"] for p in entry.get("parameters", [])],
                       vm.time_coordinates(declared, coords), relations, orders,
                       vm.HELD.get((metric_id, entry["id"]), ()))
    g = vm.metric_from_line_element(reader, entry["line_element"], coords)
    symbols = [reader.symbol[name] for name in coords]
    # A chart kept to an order, as Hartle and Thorne's is to the second in the spin, has its
    # tensors cut to that order as the checker cuts them, and vanishes where they do.
    geometry = vm.Geometry(g, symbols, seconds, reader if reader.order else None)
    n = len(coords)

    def zero(expression):
        return vm.norm(reader.surface(sp.sympify(expression))) == 0

    ricci = geometry.ricci_ll()
    ricci_flat = all(zero(ricci[a][b]) for a in range(n) for b in range(a, n))

    einstein = None
    if not ricci_flat:
        k = vm.norm(reader.surface(geometry.ricci_scalar())) / n
        constant = not (k.free_symbols & set(symbols)) and not k.atoms(sp.core.function.AppliedUndef) \
            and not k.atoms(sp.DiracDelta, sp.sign, sp.Abs)
        if constant and all(zero(ricci[a][b] - k * g[a, b]) for a in range(n) for b in range(a, n)):
            einstein = {"k": sp.sstr(sp.factor(k)), "sign": _sign(k)}

    if n >= 4:
        weyl = geometry.weyl_llll()
        conformally_flat = all(
            zero(weyl[a][b][c][d])
            for a in range(n) for b in range(a + 1, n) for c in range(n) for d in range(c + 1, n)
        )
    elif n == 3:
        conformally_flat = _cotton_vanishes(geometry, zero)
    else:
        # Every metric of two dimensions is conformally flat.
        conformally_flat = True

    killing = []
    for i, x in enumerate(symbols):
        if any(g[a, b].has(x) for a in range(n) for b in range(n)):
            continue
        timelike = _negative_somewhere(reader.surface(g[i, i]), reader, entry)
        xi = [g[a, i] for a in range(n)]
        curl = [[sp.diff(xi[b], symbols[a]) - sp.diff(xi[a], symbols[b]) for b in range(n)]
                for a in range(n)]
        orthogonal = all(
            zero(xi[a] * curl[b][c] + xi[b] * curl[c][a] + xi[c] * curl[a][b])
            for a in range(n) for b in range(a + 1, n) for c in range(b + 1, n)
        )
        killing.append({"coord": coords[i], "timelike": timelike,
                        "hypersurface_orthogonal": orthogonal,
                        "translation": is_translation(entry, coords[i])})

    sphere = _sphere(g, symbols, coords, zero, lambda angles: is_whole_sphere(entry, angles))

    return {
        "stamp": stamp(entry),
        "relations": relations,
        "orders": list(orders) if orders else None,
        "dimension": n,
        "ricci_flat": ricci_flat,
        "einstein": einstein,
        "conformally_flat": conformally_flat,
        "killing": killing,
        "sphere": sphere,
    }


def _sign(k):
    """+1 or -1 where k has one sign for every positive value of its parameters, or None.

    A parameter spelled Lambda is the cosmological constant itself, which the collection
    lets take either sign, so it is left free and a k that rests on it has no sign here.
    """
    import sympy as sp
    swap = {s: (sp.Symbol(s.name, real=True) if "Lambda" in s.name else sp.Symbol(s.name, positive=True))
            for s in k.free_symbols}
    value = sp.factor(k.subs(swap))
    if value.is_positive:
        return 1
    if value.is_negative:
        return -1
    return None


def _cotton_vanishes(geometry, zero):
    """In three dimensions: C_{abc} = D_c S_{ab} - D_b S_{ac}, S the Schouten tensor."""
    import sympy as sp
    n, g, x = geometry.n, geometry.g, geometry.coords
    gamma = geometry.christoffel_ull()
    ricci = geometry.ricci_ll()
    scalar = geometry.ricci_scalar()
    schouten = [[ricci[a][b] - scalar * g[a, b] / (2 * (n - 1)) for b in range(n)] for a in range(n)]

    def derivative(a, b, c):
        return sp.diff(schouten[a][b], x[c]) - sum(
            gamma[d][c][a] * schouten[d][b] + gamma[d][c][b] * schouten[a][d] for d in range(n))

    return all(zero(derivative(a, b, c) - derivative(a, c, b))
               for a in range(n) for b in range(n) for c in range(b + 1, n))


def _negative_somewhere(norm_squared, reader, entry):
    """Whether a Killing vector's norm is negative at one of the sampled points of the chart.

    The same points every run, drawn with a fixed seed. Each parameter and each free
    function takes a positive value between 1/50 and 50. Each coordinate is drawn inside
    the interval its domain gives it, where the reader can read the ends of that interval:
    between the two when both are finite, within a factor of 50 of the finite one when
    the other is infinite, and of either sign when both are. An end the reader cannot
    read, as r_+ is, is taken to be 0 below and infinite above. A point where the norm is
    not a real number is passed over.

    So a norm that is negative only where the parameters leave the range the entry means,
    or only beyond an end that could not be read, can read as timelike, and one negative
    only on a set the points miss as not. OVERRULED names the spacetimes where the answer
    is overruled, each with its reason.
    """
    import sympy as sp
    if norm_squared == 0:
        return False
    functions = sorted(norm_squared.atoms(sp.core.function.AppliedUndef), key=sp.sstr)
    held = {f: sp.Symbol(f"F{i}", positive=True) for i, f in enumerate(functions)}
    expression = norm_squared.subs(held)
    if expression.atoms(sp.Derivative, sp.DiracDelta):
        return False
    coordinate = {reader.symbol[name]: name for name in reader.coords}
    ends = {symbol: _ends(reader, entry, name) for symbol, name in coordinate.items()}
    wanted = set(expression.free_symbols)
    for lower, upper in ends.values():
        for end in (lower, upper):
            if end is not None and end is not sp.oo and end is not -sp.oo:
                wanted |= end.free_symbols
    wanted |= set(coordinate)
    ordered = sorted(wanted, key=lambda s: s.name)
    arguments = sorted(expression.free_symbols, key=lambda s: s.name)
    evaluate = sp.lambdify(arguments, expression, "mpmath")
    source = random.Random(SEED)

    def value_of(end, point):
        try:
            value = complex(end.subs(point).evalf())
        except (TypeError, ValueError):
            return None
        return value.real if abs(value.imag) < 1e-12 else None

    for _ in range(SAMPLES):
        point = {s: sp.Float(50 ** source.uniform(-1, 1)) for s in ordered if s not in coordinate}
        left = [s for s in ordered if s in coordinate]
        inside = True
        while left and inside:
            ready = [s for s in left
                     if all(e is None or e in (sp.oo, -sp.oo) or e.free_symbols <= set(point)
                            for e in ends[s])] or left[:1]
            for symbol in ready:
                left.remove(symbol)
                spread, share = 50 ** source.uniform(-1, 1), source.random()
                lower, upper = ends[symbol]
                low = None if lower in (None, -sp.oo) else value_of(lower, point) \
                    if lower.free_symbols <= set(point) else None
                high = None if upper in (None, sp.oo) else value_of(upper, point) \
                    if upper.free_symbols <= set(point) else None
                if lower is None or (lower is not -sp.oo and low is None):
                    low, lower = 0.0, sp.Integer(0)
                if lower is -sp.oo and high is None:
                    value = spread if share < 0.5 else -spread
                elif lower is -sp.oo:
                    value = high - spread
                elif high is None:
                    value = low + spread
                elif high > low:
                    value = low + share * (high - low)
                else:
                    inside = False
                    break
                point[symbol] = sp.Float(value)
        if not inside:
            continue
        try:
            value = complex(evaluate(*(float(point[s]) for s in arguments)))
        except (TypeError, ValueError, ZeroDivisionError, OverflowError):
            continue
        if value != value:
            continue
        if abs(value.imag) < 1e-12 * max(1.0, abs(value.real)) and value.real < 0:
            return True
    return False


def _ends(reader, entry, coord):
    """The two ends of a coordinate's interval as sympy expressions; None for one not read."""
    import sympy as sp
    import verify_metrics as vm
    found = [interval(d, coord) for d in entry.get("domains", [])
             if _plain_domain(d).startswith(coord + " \\in")]
    if len(found) != 1 or found[0] is None:
        return None, None

    def read(text):
        if text == "\\infty":
            return sp.oo
        if text == "-\\infty":
            return -sp.oo
        try:
            value = reader.surface(sp.sympify(reader(text)))
        except Exception:
            return None
        if value.atoms(sp.core.function.AppliedUndef) or value.has(sp.oo, sp.zoo, sp.nan):
            return None
        return value

    return read(found[0][0]), read(found[0][1])


def _sphere(g, symbols, coords, zero, whole):
    """The longest chain of angles the chart shows a whole round sphere in, as coordinate names."""
    import sympy as sp
    n = len(symbols)

    def linked(chain):
        """The angles of the chain have no cross terms, and each one's component is the
        one before it times the square of that one's sine."""
        for i in chain:
            if any(g[i, j] != 0 for j in range(n) if j != i):
                return False
        expected = g[chain[0], chain[0]]
        for previous, i in zip(chain, chain[1:]):
            expected = expected * sp.sin(symbols[previous]) ** 2
            if not zero(g[i, i] - expected):
                return False
        return True

    def alone(chain):
        """Nothing but the sphere itself depends on the angles of the chain."""
        angles = [symbols[i] for i in chain]
        rest = [i for i in range(n) if i not in chain]
        return not g[chain[0], chain[0]].has(*angles) and not any(
            g[a, b].has(*angles) for a in rest for b in rest)

    best = []

    def extend(chain):
        nonlocal best
        if len(chain) >= 2:
            if not linked(chain):
                return
            if len(chain) > len(best) and alone(chain) and whole([coords[i] for i in chain]):
                best = list(chain)
        for i in range(n):
            if i not in chain:
                extend(chain + [i])

    for start in range(n):
        extend([start])
    return [coords[i] for i in best] or None


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--system", action="append", default=[],
                        help="compute only <metric_id>/<system_id>, repeatable")
    parser.add_argument("--all", action="store_true", help="compute every chart again")
    parser.add_argument("--check", action="store_true",
                        help="compute, compare with metric_tags.json and write nothing")
    parser.add_argument("--budget", type=float, default=300,
                        help="seconds sympy may spend on one tensor")
    arguments = parser.parse_args()

    sys.path.insert(0, str(HERE))
    import verify_metrics as vm

    facts = load_facts()
    seen, differing, unfinished = set(), [], []
    for metric_id, entry in published_charts():
        name = f"{metric_id}/{entry['id']}"
        seen.add(name)
        if arguments.system and name not in arguments.system:
            continue
        relations = vm.PARAMETER_RELATIONS.get((metric_id, entry["id"]), {})
        current = facts.get(name)
        orders = vm.ORDERS.get((metric_id, entry["id"]))
        fresh = (current is not None and current["stamp"] == stamp(entry)
                 and current["relations"] == relations
                 and current.get("orders") == (list(orders) if orders else None))
        if fresh and not (arguments.all or arguments.system or arguments.check):
            continue
        started = time.time()
        try:
            computed = compute(metric_id, entry, arguments.budget)
        except vm.Timeout as error:
            unfinished.append(f"{name}: {error}")
            print(f"  UNFINISHED {name}: {error}", flush=True)
            continue
        print(f"  {name}  {time.time() - started:.1f}s", flush=True)
        if current != computed:
            differing.append(name)
        if not arguments.check:
            facts[name] = computed
            write_facts(facts)

    if not arguments.check and not arguments.system:
        gone = [name for name in facts if name not in seen]
        for name in gone:
            del facts[name]
        if gone:
            write_facts(facts)
    if differing:
        verb = "differ from" if arguments.check else "were written to"
        print(f"\n{len(differing)} charts {verb} {FACTS.name}:")
        for name in differing:
            print(f"  {name}")
    if unfinished:
        print(f"\n{len(unfinished)} UNFINISHED:", file=sys.stderr)
        return 1
    return 1 if arguments.check and differing else 0


if __name__ == "__main__":
    sys.exit(main())
