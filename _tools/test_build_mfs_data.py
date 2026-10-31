#!/usr/bin/env python3
"""Tests for the MFS data generator.

    .venv.noindex/bin/python -m unittest discover -s _tools

MFS_SPACETIME checks the files of one spacetime alone, or of several named with commas between
them, and passes over the tests of the page's own sources, which a new spacetime does not touch:

    MFS_SPACETIME=kerr .venv.noindex/bin/python -m unittest discover -s _tools -p test_build_mfs_data.py
"""

import ast
import collections
import contextlib
import copy
import io
import itertools
import json
import cmath
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
import unicodedata
import unittest
from pathlib import Path
from unittest import mock

import build_mfs_data as build
from voice import (CAPTION_VOICE, CONTRAST, COST, GIVES, MACHINERY, MACHINERY_OUTSIDE_A_HISTORY,
                   METRIC_GIVES, VOICE, refused_hyphens, without_mathematics_or_quotations)

# The dashes the prose is not allowed to carry: U+2010 to U+2015, which are the
# hyphen, the non breaking hyphen, the figure dash, the en dash, the em dash and
# the horizontal bar, together with U+2212, the minus sign. The ASCII hyphen at
# U+002D is deliberately absent, because it spells Reissner-Nordstrom, Kerr-Newman,
# Lanczos-van Stockum, anti-de Sitter and pp-wave, where it is part of the name.
# _tools/README.md carries the rule this stands for.
DASHES = "\u2010\u2011\u2012\u2013\u2014\u2015\u2212"

# The templates and pages whose words reach a reader, beside the generated files, and the
# data the site hands to agents.
TEMPLATES = ("_layouts/*.html", "_includes/*.html", "_includes/*.txt", "MFS/*.markdown", "*.markdown", "llms*.txt",
             "data/*.json")


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


prose = build.metric_prose
diagram_prose = build.diagram_prose
conformal_prose = build.conformal_prose
embedding_prose = build.embedding_prose
curvature_prose = build.curvature_prose


def load_folder(files):
    """Load metrics from a scratch folder holding `files`, a map of file name to metric."""
    with tempfile.TemporaryDirectory() as folder:
        for filename, metric in files.items():
            (Path(folder) / filename).write_text(json.dumps(metric), encoding="utf-8")
        with mock.patch.object(build, "METRICS_DIR", Path(folder)):
            return build.load_metrics()


def spacetimes_named(value):
    """The spacetimes a value of MFS_SPACETIME names, refusing a name with no metric file, which
    would check nothing and pass."""
    names = frozenset(name.strip() for name in value.split(",") if name.strip())
    unknown = names - {p.stem for p in build.METRICS_DIR.glob("*.json")}
    if unknown:
        raise ValueError(f"MFS_SPACETIME names {', '.join(sorted(unknown))}, which has no file in "
                         f"{build.METRICS_DIR.relative_to(build.ROOT)}")
    return names


# The spacetimes MFS_SPACETIME names, so that a worker adding one checks its own files alone. Empty,
# as when it is unset, it chooses every spacetime, and a run checks everything.
CHOSEN = spacetimes_named(os.environ.get("MFS_SPACETIME", ""))


def chosen(name):
    """Whether this run checks the spacetime `name`."""
    return not CHOSEN or name in CHOSEN


def checks(*names):
    """Pass over a test of the spacetimes `names` alone when MFS_SPACETIME chooses none of them."""
    return unittest.skipUnless(any(map(chosen, names)), f"checks {', '.join(names)} alone")


# The tests of the page's own sources, its layout, stylesheets and scripts, which a new spacetime
# does not touch, passed over when MFS_SPACETIME chooses spacetimes.
page_source = unittest.skipIf(CHOSEN, "holds the page's own sources, which a spacetime's files do not touch")


# Each spacetime's exceptions to the tests every spacetime is held to, every one with its reason, in
# _tools/derivations/exceptions/<id>.json, which _tools/README.md, "Exceptions", describes.
EXCEPTIONS_DIR = build.ROOT / "_tools" / "derivations" / "exceptions"
EXCEPTIONS = {p.stem: read(p) for p in sorted(EXCEPTIONS_DIR.glob("*.json"))}


def exceptions(field):
    """Each group of one field of every spacetime's exceptions, with its spacetime."""
    return [(name, group) for name, data in EXCEPTIONS.items() for group in data.get(field, [])]


def found(test, things, what):
    """Fail a test that found no `what` to check, or pass over it where the spacetimes MFS_SPACETIME
    chooses have none."""
    if not things:
        if CHOSEN:
            test.skipTest(f"{', '.join(sorted(CHOSEN))} has no {what}")
        test.fail(f"no {what} was found, so nothing was checked")


def counted(test, count, least, what):
    """Hold a sweep over every spacetime to having checked more than `least` of `what`, or one over
    the spacetimes MFS_SPACETIME chooses to having checked any, passing over it where they have none."""
    if CHOSEN:
        found(test, count, what)
    else:
        test.assertGreater(count, least, f"the {what} were not all checked")


def chosen_entries(table, spacetime=lambda key: key):
    """The entries of a table kept by spacetime, as `spacetime` reads one from a key, that this run
    checks."""
    return {key: value for key, value in table.items() if chosen(spacetime(key))}


def chosen_metrics():
    """The metric files of the spacetimes this run checks."""
    return [metric for metric in build.load_metrics() if chosen(metric["id"])]


def folder_files(folder, every=False):
    """The parsed files of one folder by spacetime, of the spacetimes this run checks, or of every
    one for a test that takes them as the input it spoils."""
    return {p.stem: read(p) for p in sorted(folder.glob("*.json")) if every or chosen(p.stem)}


def diagram_files(every=False):
    return folder_files(build.DIAGRAMS_DIR, every)


def conformal_files():
    return folder_files(build.CONFORMAL_DIR)


def embedding_files():
    return folder_files(build.EMBEDDING_DIR)


class NamedFiles(dict):
    """The parsed files of one folder by spacetime, each read the first time a test names it,
    whichever spacetimes MFS_SPACETIME chooses."""

    def __init__(self, folder):
        super().__init__()
        self.folder = folder

    def __missing__(self, name):
        self[name] = read(self.folder / f"{name}.json")
        return self[name]


def diagram_file(name):
    """One spacetime's diagram file, whichever spacetimes this run checks."""
    return read(build.DIAGRAMS_DIR / f"{name}.json")


def conformal_file(name):
    return read(build.CONFORMAL_DIR / f"{name}.json")


def embedding_file(name):
    return read(build.EMBEDDING_DIR / f"{name}.json")


def grid_nodes(piece):
    """The nodes (X, Y, Z) of a grid piece, row i of u and column j of v, a polar grid's X and Y
    being u cos v and u sin v, as _tools/README.md defines them."""
    grid = piece["grid"]
    if grid["frame"] == "ellipses":
        # Row i the ellipse (a_i cos v, b_i sin v) at the height z_i, sheared where the grid gives
        # sx and sy: (a_i cos v + sx_i sin v, sy_i cos v + b_i sin v).
        zero = [0.0] * len(grid["u"])
        return [[(a * math.cos(v) + sx * math.sin(v), sy * math.cos(v) + b * math.sin(v), z) for v in grid["v"]]
                for a, b, z, sx, sy in zip(grid["a"], grid["b"], grid["z"], grid.get("sx", zero), grid.get("sy", zero))]
    polar = grid["frame"] == "polar"
    return [[(u * math.cos(v), u * math.sin(v), z) if polar else (u, v, z)
             for v, z in zip(grid["v"], row)] for u, row in zip(grid["u"], grid["z"])]


def grid_height(piece, X, Y, reach=1e-6):
    """The height of a grid piece's triangles over the point (X, Y), each cell cut along its
    diagonal from (i, j) to (i + 1, j + 1), or None where none lies within `reach` of it, which
    allows the rounding of a point written on an edge of a thin triangle. The cell is found by u
    and v and checked with its neighbours, since a polar grid's cells have straight sides."""
    grid, P = piece["grid"], grid_nodes(piece)
    polar = grid["frame"] == "polar"
    m, n = len(grid["u"]), len(grid["v"])
    u, v = (math.hypot(X, Y), math.atan2(Y, X) % (2 * math.pi)) if polar else (X, Y)
    i = max(0, min(m - 2, sum(1 for a in grid["u"] if a <= u) - 1))
    j = sum(1 for b in grid["v"] if b <= v) - 1
    best = None
    for di in (0, -1, 1):
        for dj in (0, -1, 1):
            ii, jj = i + di, j + dj
            if not 0 <= ii < m - 1:
                continue
            if polar:
                jj %= n
            elif not 0 <= jj < n - 1:
                continue
            j1 = (jj + 1) % n
            a, b, c, d = P[ii][jj], P[ii + 1][jj], P[ii + 1][j1], P[ii][j1]
            for A, B, C in ((a, b, c), (a, c, d)):
                area = (B[0] - A[0]) * (C[1] - A[1]) - (C[0] - A[0]) * (B[1] - A[1])
                if abs(area) < 1e-15:
                    continue
                wa = ((B[0] - X) * (C[1] - Y) - (C[0] - X) * (B[1] - Y)) / area
                wb = ((C[0] - X) * (A[1] - Y) - (A[0] - X) * (C[1] - Y)) / area
                wc = 1 - wa - wb
                # How far outside the triangle the point lies, along the normal of each edge it
                # lies beyond: a weight below zero is that distance over the corner's height.
                edges = ((wa, B, C), (wb, C, A), (wc, A, B))
                off = max(-w * abs(area) / math.hypot(Q[0] - P[0], Q[1] - P[1]) for w, P, Q in edges)
                if off <= reach and (best is None or off < best[0]):
                    best = (off, wa * A[2] + wb * B[2] + wc * C[2])
    return best and best[1]


def curvature_files():
    return folder_files(build.CURVATURE_DIR)


BARDEEN_G2 = 1 / 9                # g^2 in r_s^2, as every diagram of Bardeen's black hole takes it


def bardeen_f(r):
    """Bardeen's f = 1 - r_s r^2/(r^2 + g^2)^(3/2) at r_s = 1 and g = 1/3."""
    return 1 - r * r / (r * r + BARDEEN_G2) ** 1.5


def bardeen_horizons():
    """The two positive zeros of f, by bisection on each side of its minimum at r = sqrt(2) g."""
    turn = math.sqrt(2 * BARDEEN_G2)
    roots = []
    for lo, hi in ((1e-6, turn), (turn, 2.0)):
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if bardeen_f(lo) * bardeen_f(mid) > 0 else (lo, mid)
        roots.append(0.5 * (lo + hi))
    return roots


QOS_ALPHA = 1.25                  # alpha in units of m^2, as every diagram of the quantum Oppenheimer-Snyder hole takes it


def qos_roots():
    """The zeros of r^4 - 2r^3 + alpha at alpha = 5/4, without numpy: the two real ones by bisection and
    the complex pair from the quadratic left over, r^2 + br + c with b = r_- + r_+ - 2, c = alpha/(r_- r_+)."""
    def P(r):
        return r ** 4 - 2 * r ** 3 + QOS_ALPHA
    real = []
    for lo, hi in ((1.0, 1.5), (1.5, 2.0)):
        for _ in range(200):
            mid = (lo + hi) / 2
            lo, hi = (mid, hi) if P(lo) * P(mid) > 0 else (lo, mid)
        real.append((lo + hi) / 2)
    b, c = real[0] + real[1] - 2, QOS_ALPHA / (real[0] * real[1])
    return real, complex(-b / 2, math.sqrt(4 * c - b * b) / 2)


def qos_rstar(r):
    """The tortoise coordinate at m = 1 and alpha = 5/4, dr_*/dr = 1/f, zero at r_b = (alpha/2)^(1/3)."""
    real, z = qos_roots()

    def raw(x):
        out = x + sum(w ** 4 / (4 * w ** 3 - 6 * w ** 2) * math.log(abs(x - w)) for w in real)
        return out + 2 * (z ** 4 / (4 * z ** 3 - 6 * z ** 2) * cmath.log(x - z)).real
    return raw(r) - raw((QOS_ALPHA / 2) ** (1 / 3))


def qos_lead(r, n=2000):
    """The integral of 1/(1 + N) from r_b to r, N = sqrt(2/r - alpha/r^4), by Simpson's rule in the
    variable s with r = r_b + s^2, which takes away the square root N has at r_b."""
    rb = (QOS_ALPHA / 2) ** (1 / 3)
    top = math.sqrt(max(r - rb, 0.0))
    if top == 0:
        return 0.0

    def g(s):
        x = rb + s * s
        return 2 * s / (1 + math.sqrt(max(2 / x - QOS_ALPHA / x ** 4, 0.0)))
    h = top / n
    return h / 3 * (g(0) + g(top) + sum((4 if k % 2 else 2) * g(k * h) for k in range(1, n)))


def eos_rstar(r):
    """The eye of the storm's tortoise coordinate on its axis at m = 1, a = 3/5 and ell = 1/5,
    dr_*/dr = (r^2 + a^2)/Delta with Delta = r^2 + a^2 - 2r e^(-ell/r), zero at r = 0, without numpy:
    the horizons found by bisection, the logarithm of each, and the integral of what is left by
    Gauss and Legendre's rule of five points on panels that end on the horizons."""
    a2, ell = 0.36, 0.2

    def delta(x):
        return x * x + a2 - 2 * x * math.exp(-ell / x)

    def root(lo, hi):
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if (delta(mid) > 0) == (delta(lo) > 0) else (lo, mid)
        return 0.5 * (lo + hi)
    poles = [(ri, (ri * ri + a2) / (2 * ri - 2 * math.exp(-ell / ri) * (1 + ell / ri)))
             for ri in (root(0.3, 0.9), root(1.0, 1.8))]
    nodes = (0.0, 0.5384693101056831, -0.5384693101056831, 0.906179845938664, -0.906179845938664)
    weights = (0.5688888888888889, 0.47862867049936647, 0.47862867049936647, 0.23692688505618908, 0.23692688505618908)

    def smooth(x):
        return (x * x + a2) / delta(x) - sum(c / (x - ri) for ri, c in poles)
    cuts = sorted({0.0, r} | {ri for ri, _ in poles if ri < r})
    total = 0.0
    for lo, hi in zip(cuts, cuts[1:]):
        panels = max(1, math.ceil((hi - lo) / 0.05))
        width = (hi - lo) / panels
        for k in range(panels):
            mid = lo + (k + 0.5) * width
            total += 0.5 * width * sum(w * smooth(mid + 0.5 * width * n) for n, w in zip(nodes, weights))
    return total + sum(c * math.log(abs(1 - r / ri)) for ri, c in poles)


def bardeen_rstar(r):
    """Bardeen's tortoise coordinate, dr_*/dr = 1/f and r_* = 0 at the centre, without numpy: the
    logarithm of each horizon, ln|1 - r/r_i|/f'(r_i), and the integral from 0 of what is left of
    1/f, which is smooth, by Gauss and Legendre's rule of five points on panels that end on the
    horizons, so that no point of the rule comes near a pole."""
    poles = [(ri, (ri * ri + BARDEEN_G2) ** 2.5 / (ri * (ri * ri - 2 * BARDEEN_G2))) for ri in bardeen_horizons()]
    nodes = (0.0, 0.5384693101056831, -0.5384693101056831, 0.906179845938664, -0.906179845938664)
    weights = (0.5688888888888889, 0.47862867049936647, 0.47862867049936647, 0.23692688505618908, 0.23692688505618908)

    def smooth(x):
        return 1 / bardeen_f(x) - sum(a / (x - ri) for ri, a in poles)
    cuts = sorted({0.0, r} | {ri for ri, _ in poles if ri < r})
    total = 0.0
    for lo, hi in zip(cuts, cuts[1:]):
        panels = max(1, math.ceil((hi - lo) / 0.05))
        width = (hi - lo) / panels
        for k in range(panels):
            mid = lo + (k + 0.5) * width
            total += 0.5 * width * sum(w * smooth(mid + 0.5 * width * n) for n, w in zip(nodes, weights))
    return total + sum(a * math.log(abs(1 - r / ri)) for ri, a in poles if r != ri)


# Born and Infeld's point charge as every diagram draws it, in units of r_0 at r_q = r_0/2: (r_s, r_q)
# of Hoffmann's particle, whose whole mass is the energy of its field, and of the black hole.
BORN_INFELD = {"particle": (math.gamma(0.25) ** 2 / (12 * math.sqrt(math.pi)), 0.5), "hole": (2.0, 0.5)}


def carlson_rf(x, y, z):
    """Carlson's symmetric elliptic integral R_F, by his duplication theorem."""
    for _ in range(60):
        lam = math.sqrt(x * y) + math.sqrt(y * z) + math.sqrt(z * x)
        x, y, z = (x + lam) / 4, (y + lam) / 4, (z + lam) / 4
        mean = (x + y + z) / 3
        if max(abs(x - mean), abs(y - mean), abs(z - mean)) < 1e-9 * mean:
            break
    X, Y, Z = 1 - x / mean, 1 - y / mean, 1 - z / mean
    e2, e3 = X * Y - Z * Z, X * Y * Z
    return (1 - e2 / 10 + e3 / 14 + e2 * e2 / 24 - 3 * e2 * e3 / 44) / math.sqrt(mean)


def elliptic_f(phi, m):
    """The incomplete elliptic integral of the first kind, F(phi | m), for an amplitude up to pi."""
    if phi > math.pi / 2:
        return 2 * carlson_rf(0.0, 1 - m, 1.0) - elliptic_f(math.pi - phi, m)
    s = math.sin(phi)
    return s * carlson_rf(1 - s * s, 1 - m * s * s, 1.0)


def born_infeld_f(r, case):
    """g^rr = 1 - 2m/r of Born and Infeld's point charge, with the mass function as the charts define it,
    m = r_s/2 + r_q^2 r/(3 (r^2 + W)) - (r_q^2/(3 r_0)) F(2 arctan(r_0/r) | 1/2), W = sqrt(r^4 + r_0^4)."""
    rs, rq = BORN_INFELD[case]
    W = math.sqrt(r ** 4 + 1)
    m = rs / 2 + rq ** 2 * r / (3 * (r * r + W)) - rq ** 2 / 3 * elliptic_f(2 * math.atan2(1.0, r), 0.5)
    return 1 - 2 * m / r


def born_infeld_horizon():
    """The black hole's one horizon, the zero of g^rr, by bisection."""
    lo, hi = 1.5, 2.2
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if born_infeld_f(lo, "hole") * born_infeld_f(mid, "hole") > 0 else (lo, mid)
    return 0.5 * (lo + hi)


def born_infeld_rstar(r, case):
    """The tortoise coordinate of Born and Infeld's point charge, dr_*/dr = 1/f and r_* = 0 at the
    centre, without numpy, as Bardeen's is summed above: the particle's 1/f is smooth, and the black
    hole's has one pole, at the horizon, of residue 1/f'(r_h) = r_h/(1 - 2 r_q^2/(r_h^2 + W))."""
    poles = []
    if case == "hole":
        rh, rq = born_infeld_horizon(), BORN_INFELD[case][1]
        poles = [(rh, rh / (1 - 2 * rq ** 2 / (rh * rh + math.sqrt(rh ** 4 + 1))))]
    nodes = (0.0, 0.5384693101056831, -0.5384693101056831, 0.906179845938664, -0.906179845938664)
    weights = (0.5688888888888889, 0.47862867049936647, 0.47862867049936647, 0.23692688505618908, 0.23692688505618908)

    def smooth(x):
        return 1 / born_infeld_f(x, case) - sum(a / (x - ri) for ri, a in poles)
    cuts = sorted({0.0, r} | {ri for ri, _ in poles if ri < r})
    total = 0.0
    for lo, hi in zip(cuts, cuts[1:]):
        panels = max(1, math.ceil((hi - lo) / 0.05))
        width = (hi - lo) / panels
        for k in range(panels):
            mid = lo + (k + 0.5) * width
            total += 0.5 * width * sum(w * smooth(mid + 0.5 * width * n) for n, w in zip(nodes, weights))
    return total + sum(a * math.log(abs(1 - r / ri)) for ri, a in poles if r != ri)


class PublishedFilesAreCurrent(unittest.TestCase):
    def test_check_mode_passes(self):
        self.assertEqual(build.main(["--check"]), 0)


class Index(unittest.TestCase):
    def setUp(self):
        self.metrics = chosen_metrics()
        self.index = [entry for entry in read(build.INDEX_FILE) if chosen(entry["id"])]

    def test_one_entry_per_file_and_no_more(self):
        self.assertEqual(
            [e["id"] for e in self.index],
            [m["id"] for m in self.metrics],
        )
        self.assertEqual(
            sorted(e["id"] for e in self.index),
            sorted(p.stem for p in build.METRICS_DIR.glob("*.json") if chosen(p.stem)),
        )

    def test_entry_carries_the_search_fields_and_a_stamp(self):
        diagrams, conformal, embedding = diagram_files(), conformal_files(), embedding_files()
        curvature = curvature_files()
        for entry, metric in zip(self.index, self.metrics):
            expected = ({"id", "name", "tags", "version"} | ({"diagrams"} if metric["id"] in diagrams else set())
                        | ({"conformal"} if metric["id"] in conformal else set())
                        | ({"embedding"} if metric["id"] in embedding else set())
                        | ({"curvature"} if metric["id"] in curvature else set()))
            self.assertEqual(set(entry), expected)
            self.assertEqual(entry["name"], metric["short_name"])
            self.assertEqual(entry["tags"], metric["tags"])
            self.assertEqual(entry["version"], build.content_version(metric))

    def test_a_spacetime_with_diagrams_carries_their_stamp(self):
        diagrams = diagram_files()
        found(self, diagrams, "diagram file")
        stamped = {e["id"]: e["diagrams"] for e in self.index if "diagrams" in e}
        self.assertEqual(set(stamped), set(diagrams))
        for metric_id, diagram in diagrams.items():
            self.assertEqual(stamped[metric_id], build.content_version(diagram))

    def test_a_spacetime_with_a_conformal_diagram_carries_its_stamp(self):
        conformal = conformal_files()
        found(self, conformal, "conformal diagram file")
        stamped = {e["id"]: e["conformal"] for e in self.index if "conformal" in e}
        self.assertEqual(set(stamped), set(conformal))
        for metric_id, data in conformal.items():
            self.assertEqual(stamped[metric_id], build.content_version(data))

    def test_a_spacetime_with_an_embedding_diagram_carries_its_stamp(self):
        embedding = embedding_files()
        found(self, embedding, "embedding diagram file")
        stamped = {e["id"]: e["embedding"] for e in self.index if "embedding" in e}
        self.assertEqual(set(stamped), set(embedding))
        for metric_id, data in embedding.items():
            self.assertEqual(stamped[metric_id], build.content_version(data))

    def test_a_spacetime_with_a_curvature_figure_carries_its_stamp(self):
        curvature = curvature_files()
        found(self, curvature, "curvature figure file")
        stamped = {e["id"]: e["curvature"] for e in self.index if "curvature" in e}
        self.assertEqual(set(stamped), set(curvature))
        for metric_id, data in curvature.items():
            self.assertEqual(stamped[metric_id], build.content_version(data))

    def test_a_metric_missing_its_short_name_is_refused(self):
        with self.assertRaises(build.DataError):
            self.load_one("x.json", {"id": "x", "name": "X", "tags": ["t"]})

    def test_a_metric_whose_id_is_not_its_filename_is_refused(self):
        with self.assertRaises(build.DataError):
            self.load_one("x.json", {"id": "y", "name": "X", "short_name": "X", "tags": ["t"]})

    def test_a_metric_with_a_sort_name_is_refused(self):
        with self.assertRaises(build.DataError):
            self.load_one("x.json", {"id": "x", "name": "X", "short_name": "X", "sort_name": "A", "tags": ["t"]})

    def test_a_cloud_sync_conflict_copy_is_passed_over(self):
        good = {"id": "x", "name": "X", "short_name": "X", "tags": ["t"]}
        stale = dict(good, name="An older X")
        loaded = self.load_folder({"x.json": good, "x 2.json": stale})
        self.assertEqual([m["id"] for m in loaded], ["x"])
        self.assertEqual(loaded[0]["name"], "X")

    def test_a_normally_named_file_with_a_wrong_id_still_fails(self):
        with self.assertRaises(build.DataError):
            self.load_folder({"x2.json": {"id": "y", "name": "X", "short_name": "X", "tags": ["t"]}})

    def load_one(self, filename, metric):
        return self.load_folder({filename: metric})

    def load_folder(self, files):
        return load_folder(files)


@page_source
class ListOrder(unittest.TestCase):
    """The search list reads in alphabetical order of the names it shows."""

    def test_the_published_list_is_in_order_of_its_displayed_names(self):
        names = [e["name"] for e in read(build.INDEX_FILE)]
        for before, after in zip(names, names[1:]):
            self.assertLessEqual(build.name_key(before), build.name_key(after), f"{before!r} is listed before {after!r}")

    def test_the_key_ignores_case_and_accents_and_nothing_else(self):
        self.assertEqual(build.name_key("Gödel"), build.name_key("Godel"))
        self.assertEqual(build.name_key("Natário"), build.name_key("Natario"))
        self.assertEqual(build.name_key("Reissner-Nordström"), build.name_key("reissner-nordstrom"))
        in_order = ["Anti-de Sitter", "Bianchi", "de Sitter", "Ellis-Bronnikov", "Gödel", "Kasner",
                    "Natário", "Oppenheimer-Snyder", "pp-wave", "Reissner-Nordström", "Vilenkin-Gott"]
        self.assertEqual(sorted(reversed(in_order), key=build.name_key), in_order)

    def test_a_new_spacetime_takes_its_place_by_name_not_by_id(self):
        folder = {
            "a.json": {"id": "a", "name": "Z", "short_name": "Zeta", "tags": ["t"]},
            "b.json": {"id": "b", "name": "E", "short_name": "Éta", "tags": ["t"]},
            "c.json": {"id": "c", "name": "A", "short_name": "alpha", "tags": ["t"]},
        }
        loaded = load_folder(folder)
        self.assertEqual([m["short_name"] for m in loaded], ["alpha", "Éta", "Zeta"])


class Stamps(unittest.TestCase):
    def setUp(self):
        self.metric = build.load_metrics()[0]

    def test_content_change_moves_the_stamp(self):
        changed = copy.deepcopy(self.metric)
        changed["history"] = changed.get("history", "") + "."
        self.assertNotEqual(build.content_version(self.metric), build.content_version(changed))

    def test_reordering_the_file_leaves_the_stamp_alone(self):
        reordered = dict(reversed(list(self.metric.items())))
        self.assertEqual(build.content_version(self.metric), build.content_version(reordered))

    def test_stamps_are_distinct(self):
        stamps = [build.content_version(m) for m in build.load_metrics()]
        self.assertEqual(len(set(stamps)), len(stamps))


class Prose(unittest.TestCase):
    def test_no_metric_on_disk_carries_a_dash_in_its_prose(self):
        fields = [(f"{m['id']}.json: {field}", value)
                  for m in chosen_metrics() for field, value in prose(m)]
        self.assert_no_dashes(fields)

    def test_no_diagram_on_disk_carries_a_dash_in_its_prose(self):
        fields = [(field, value) for name, diagram in diagram_files().items()
                  for field, value in diagram_prose(name, diagram)]
        self.assert_no_dashes(fields)

    def test_no_conformal_diagram_on_disk_carries_a_dash_in_its_prose(self):
        fields = [(field, value) for name, data in conformal_files().items()
                  for field, value in conformal_prose(name, data)]
        self.assert_no_dashes(fields)

    def test_no_conformal_diagram_spells_a_character_as_an_escape(self):
        for name, data in conformal_files().items():
            for field, value in conformal_prose(name, data):
                self.assertNotRegex(value, r"\\u[0-9a-fA-F]{4}", field)

    def test_no_embedding_diagram_on_disk_carries_a_dash_in_its_prose(self):
        fields = [(field, value) for name, data in embedding_files().items()
                  for field, value in embedding_prose(name, data)]
        self.assert_no_dashes(fields)

    def test_no_embedding_diagram_spells_a_character_as_an_escape(self):
        for name, data in embedding_files().items():
            for field, value in embedding_prose(name, data):
                self.assertNotRegex(value, r"\\u[0-9a-fA-F]{4}", field)

    def test_no_curvature_figure_on_disk_carries_a_dash_in_its_prose(self):
        fields = [(field, value) for name, data in curvature_files().items()
                  for field, value in curvature_prose(name, data)]
        self.assert_no_dashes(fields)

    def test_no_curvature_figure_spells_a_character_as_an_escape(self):
        for name, data in curvature_files().items():
            for field, value in curvature_prose(name, data):
                self.assertNotRegex(value, r"\\u[0-9a-fA-F]{4}", field)

    def test_no_metric_on_disk_spells_a_character_as_an_escape(self):
        # A \u escaped twice in the JSON reaches the page as the six characters of the
        # escape, as the î of Lemaître did in the Tolman-Bondi conventions, since neither
        # TeX nor the page reads it. The character itself is what the files carry.
        for metric in chosen_metrics():
            for field, value in prose(metric):
                self.assertNotRegex(value, r"\\u[0-9a-fA-F]{4}",
                                    f"{metric['id']}.json: {field} spells a character as an escape")

    def test_no_diagram_spells_a_character_as_an_escape(self):
        for name, diagram in diagram_files().items():
            for field, value in diagram_prose(name, diagram):
                self.assertNotRegex(value, r"\\u[0-9a-fA-F]{4}", field)

    def test_every_caption_names_its_plane_and_never_a_block(self):
        """A reader is told which slice is drawn, in words, and never handed the shorthand."""
        views = 0
        for name, diagram in diagram_files().items():
            for part in ("systems", "projections"):
                for system, drawn in diagram.get(part, {}).items():
                    for view in drawn:
                        views += 1
                        where = f"diagrams/{name}.json {system}/{view['id']}"
                        self.assertTrue(opens_with_a_noun_phrase(view["caption"][0]), where)
                        one = {"systems": {}, part: {system: [view]}}
                        for field, value in diagram_prose(name, one):
                            self.assertNotRegex(value, r"(?i)\bblocks?\b", field)
        found(self, views, "view of a spacetime diagram")

    def test_no_template_hands_a_noun_back_to_what_gives_it(self):
        """A reader is told the metric distance, never the distance the metric gives. The build
        holds the data's prose to it; the templates and pages are held here."""
        templates = sorted({p for pattern in TEMPLATES for p in build.ROOT.glob(pattern)})
        self.assertTrue(templates, "no template was read, so none was checked")
        fields = [(str(p.relative_to(build.ROOT)), " ".join(p.read_text(encoding="utf-8").split()))
                  for p in templates]
        for where, value in fields:
            self.assertNotRegex(value, GIVES, where)
            self.assertNotRegex(value, METRIC_GIVES, where)

    def test_the_noun_handed_back_is_caught_and_plain_speech_is_not(self):
        for text in ("every distance along it is the distance the metric gives.",
                     "the distances the metric gives", "against the length the metric gives the line",
                     "the multiplicities the symmetries give each of them", "the value the chart yields",
                     "the area the integral returns",
                     "along the radial geodesic the Christoffel symbols give,"):
            self.assertRegex(text, GIVES, text)
        for text in ("every distance along it is the metric distance.", "the condition he gives is",
                     "the Einstein tensor returns the density", "where $f = 1$, $ds^2 = -c^2dt^2$",
                     "for any $b_0$ the formula gives a metric",
                     "along its radial geodesic,", "reach the edge of the universe and return",
                     "The first against the second gives $16$", "the metric given above"):
            self.assertNotRegex(text, GIVES, text)
            self.assertNotRegex(text, METRIC_GIVES, text)
        for text in ("where $f = 1$ and the metric gives $ds^2 = -c^2dt^2$", "The metric gives the interval"):
            self.assertNotRegex(text, GIVES, text)
            self.assertRegex(text, METRIC_GIVES, text)

    def test_the_machinery_words_are_caught_and_the_physics_is_not(self):
        def caught(text, history=False):
            patterns = MACHINERY + (() if history else MACHINERY_OUTSIDE_A_HISTORY)
            return any(re.search(p, text) for p in patterns)
        for text in ("the statements of the entry are on it", "no published value leans on a symbol the entry has not declared",
                     "the two published blocks", "the second published chart", "with the same $\\Gamma$ printed above them",
                     "the mixed variant is the shortest", "the one component flow that `alcubierre` carries",
                     "geometries already in this collection", "the Rindler coordinates displayed here"):
            self.assertTrue(caught(text), text)
        self.assertTrue(caught("which is why no metric is displayed here", history=True))
        for text, history in (("the local building block of every generic singularity", False),
                              ("Einstein and Rosen published it in 1937", True),
                              ("the translation the Monthly Notices printed in 1931", True),
                              ("the chart covers the exterior", False), ("a uniform electromagnetic field", False)):
            self.assertFalse(caught(text, history), text)

    def assert_no_dashes(self, fields):
        found(self, fields, "prose")
        for where, value in fields:
            for character in value:
                if character in DASHES:
                    self.fail(
                        f"{where} carries "
                        f"U+{ord(character):04X} {unicodedata.name(character)}. "
                        "Write a full stop, a semicolon or a comma instead, "
                        "or rewrite the sentence so it does not want the break."
                    )


def opens_with_a_noun_phrase(caption):
    """Whether a caption opens with a noun phrase naming what is drawn, as "The plane of $t$ and
    $r$ ($\\theta = \\pi/2$, $\\phi = 0$)" does: it starts with a capital, never with "This", and
    its head, up to the first comma, parenthesis, stop or clause of its own, holds no "is" or "are"."""
    text = without_mathematics_or_quotations(caption)
    head = re.split(r"[,(.:;]|\b(?:which|where|when|that|as|whose|who|so|since)\b", text, maxsplit=1)[0]
    return (bool(re.match(r"[A-Z]", text)) and not text.startswith("This ")
            and not re.search(r"\b(?:is|are|was|were)\b", head))


class Voice(unittest.TestCase):
    """The rules of the captain's voice that a pattern can hold catch what they are for and leave
    plain physics alone; the build holds every history, convention and caption to them."""

    def test_the_caption_rules_catch_the_old_openings_and_leave_the_compact_ones(self):
        def caught(text):
            text = without_mathematics_or_quotations(text)
            return [rule for rule, pattern in CAPTION_VOICE if re.search(pattern, text)]
        for text in ("This is the whole of a ball of dust collapsing from rest at $R_0 = 2\\,r_s$, and each point of "
                     "the diagram stands for a sphere.",
                     "This is the plane of $t$ and $r$ at $\\theta = \\pi/2$ and $\\phi = 0$.",
                     "and each point of this plane stands for one such sphere.",
                     "the height stands for the tilt alone.", "so that the axes stand for $t$ and $z$"):
            self.assertTrue(caught(text), text)
        for text in ("A spherically symmetric distribution of dust collapsing from rest ($R_0 = 2\\,r_s$), each "
                     "point in the diagram a 2-sphere.",
                     "The plane of $t$ and $r$ ($\\theta = \\pi/2$, $\\phi = 0$), the same at every other angle "
                     "by spherical symmetry.",
                     "The whole Schwarzschild spacetime with the ingoing Eddington-Finkelstein coordinates $v$ "
                     "and $r$ on it.",
                     "Inside $r_s$ both edges of every future cone point to larger $r$: this is the white hole.",
                     "no surface in flat space carries the whole of it."):
            self.assertEqual(caught(text), [], text)
        self.assertFalse(opens_with_a_noun_phrase("This is the plane of $t$ and $x$ at $y = z = 0$."))
        self.assertFalse(opens_with_a_noun_phrase("The metric on this plane is $-dt^2 + dx^2$."))
        self.assertTrue(opens_with_a_noun_phrase("The plane of $t$ and $x$ ($y = z = 0$). The metric on it is flat."))
        self.assertTrue(opens_with_a_noun_phrase("The static chart along a line through the observer: $x = r$ on the "
                                                 "right is $\\phi = 0$."))

    def test_the_hyphen_rule_keeps_the_names_of_physics_and_refuses_english_compounds(self):
        for text in ("the no-boundary proposal of Hartle and Hawking", "half a 4-sphere", "for a 3-geometry",
                     "sixty e-foldings of inflation", "one e-folding", "the no-boundary proposal's sign",
                     "the 4-sphere's action", "Reissner-Nordström", "anti-de Sitter space",
                     "a pp-wave", "Taub-NUT", "a 3-sphere", "the Kerr-Schild chart [kerr1965, kerr-schild]",
                     "\"a future-directed\" quotation", "$x - y$"):
            self.assertEqual(refused_hyphens(text), [], text)
        for text, word in (("a future-directed curve", "future-directed"), ("a non-singular core", "non-singular"),
                           ("a time-like line", "time-like"), ("a well-known result", "well-known"),
                           ("half a four-sphere", "four-sphere"), ("a three-geometry", "three-geometry")):
            self.assertEqual(refused_hyphens(text), [word], text)

    def test_the_rules_catch_what_he_flagged_and_leave_plain_physics_alone(self):
        def caught(text):
            text = without_mathematics_or_quotations(text)
            return [rule for rule, pattern in VOICE if re.search(pattern, text)]
        for text in ("What comes out is a soliton you can build in pieces.",
                     "Three papers landed in the same window and were read together.", "each paper reports",
                     "That reads as a refutation,", "and that is a question with answers",
                     "so the spacetime does not stand on its own.", "every distance is the distance the metric gives",
                     "It is the surface the Morris-Thorne wormhole draws, and the same metric.",
                     "as the spacetime diagram declares.", "no conformal diagram is given because the tube is free",
                     "so nothing is drawn.", "so the Weyl tensor is not the Riemann tensor but its trace free part",
                     "every single time", "It is worth noting that $r = 0$ is regular.", "Crucially, $k$ is free.",
                     "as noted above", "we will take the chart $x^0 = ct$", "with $x$, $y$ and $z$ lengths",
                     "the mass actually inside it", "What carries the curvature is the Kretschmann scalar",
                     "a dash - trailing"):
            self.assertTrue(caught(text), text)
        for text in ("so anti-de Sitter space is not globally hyperbolic.", "every distance is the metric distance",
                     "the Einstein tensor returns the density", "with $x$, $y$, and $z$ lengths",
                     "the energy density measured by the observers who ride the slices",
                     "they differ from it in what the geometry does with the traveller",
                     "stationary black holes \"are all, every single one of them, described exactly\"",
                     "$R = -\\left(g_{\\mu\\rho}g_{\\nu\\sigma} - g_{\\mu\\sigma}g_{\\nu\\rho}\\right)/L^2$",
                     "the one parameter is the anti-de Sitter radius $L$", "whether it holds at every event"):
            self.assertEqual(caught(text), [], text)


class Contrast(unittest.TestCase):
    """No text a reader of the spacetimes page sees defines a thing by a contrast with what it is
    not, as the captain asked on 29 September 2026: no "X rather than Y", "instead of", "not a X
    but Y", "is X, not Y" or "less X, more Y". Each says what the thing is and does, plainly."""

    # A sentence where the phrase states physics plainly, by its place and the phrase, each with
    # its reason. None is needed while every text says what a thing is.
    ALLOWED = {}

    @staticmethod
    def layout_text():
        """The words of the page's layout that reach a reader: its markup's text and the strings
        its scripts set, with every comment left out."""
        page = (build.ROOT / "_layouts" / "mfs.html").read_text(encoding="utf-8")
        text = re.sub(r"<[^>]+>", " ", re.sub(r"<script\b.*?</script>|<style\b.*?</style>|<!--.*?-->", "", page, flags=re.S))
        strings = []
        for script in re.findall(r"<script\b[^>]*>(.*?)</script>", page, re.S):
            script = re.sub(r"(?m)(^|[^:\\])//.*$", r"\1", re.sub(r"/\*.*?\*/", "", script, flags=re.S))
            strings += re.findall(r"'((?:[^'\\\n]|\\.)*)'", script)
        return " ".join(text.split()) + " | " + " | ".join(strings)

    def texts(self):
        """The words of the page a reader sees; the build holds every text of the data to the rule."""
        return [("_layouts/mfs.html", self.layout_text()),
                ("MFS/index.markdown", (build.ROOT / "MFS" / "index.markdown").read_text(encoding="utf-8"))]

    def caught(self, text):
        return [m.group(0) for m in re.finditer(CONTRAST, without_mathematics_or_quotations(text))]

    def test_no_text_of_the_page_defines_a_thing_by_what_it_is_not(self):
        fields = self.texts()
        self.assertTrue(all(value.strip() for _, value in fields), "the page was not all read")
        for where, value in fields:
            for phrase in self.caught(value):
                with self.subTest(where, phrase=phrase):
                    self.assertIn((where, phrase), self.ALLOWED, f"{where} has {phrase!r}: {value[:160]!r}")

    def test_the_rule_catches_the_contrasts_and_leaves_plain_statements(self):
        for text in ("It reads like a book you can search rather than a calculator you drive.",
                     "the first warp geometries built from familiar matter rather than the forbidden kind",
                     "slides space sideways instead of compressing it.", "They began instead from what a traveler would insist on.",
                     "It is not a metric but a property of spacetimes.", "The solution is a laboratory, not a proposal.",
                     "Less calculation, more reading."):
            with self.subTest(text):
                self.assertTrue(self.caught(text))
        for text in ("Kornel Lanczos wrote down the rotating dust cylinder in 1924, not long after Einstein's field "
                     "equations had settled into their final form.",
                     "It was a clean result, not yet joined to any other solution.",
                     "the singularity theorems no longer apply", "each point in the diagram a single event.",
                     "The quote \"rather than\" is his own.", "where $a \\neq b$, not $a = b$"):
            with self.subTest(text):
                self.assertEqual(self.caught(text), [])


strings = build.strings


class NoCost(unittest.TestCase):
    """No text a reader of the spacetimes page sees says what physics requires in the words of
    money, as the captain asked on 29 September 2026: a shortcut that stays open requires exotic
    matter, and no field of any spacetime or diagram, nor the page, says cost, costs, costly or
    costing."""

    def texts(self):
        fields = Contrast().texts()
        # Every string of every file, labels, legends and mathematics among them, where the build
        # holds only the fields written in sentences to the rule.
        for folder in (build.METRICS_DIR, build.DIAGRAMS_DIR, build.CONFORMAL_DIR, build.EMBEDDING_DIR,
                       build.CURVATURE_DIR):
            for path in sorted(p for p in folder.glob("*.json") if chosen(p.stem)):
                fields += strings(read(path), path.name)
        return fields

    def test_no_text_a_reader_sees_says_cost(self):
        fields = self.texts()
        counted(self, len(fields), 10000, "texts")
        for where, value in fields:
            found = re.search(COST, value)
            self.assertIsNone(found, f"{where} has {found and found.group(0)!r}; state what is required, "
                                     "what grows or what was spent")

    def test_the_rule_catches_every_form_and_leaves_other_words_alone(self):
        for text in ("a shortcut that stays open costs exotic matter", "at the cost of seven years of work",
                     "what the trip would cost", "a costly bubble", "Costing the trip in fuel"):
            with self.subTest(text):
                self.assertRegex(text, COST)
        for text in ("$\\cos t$", "the coast of the Baltic", "a cosmic string", "Costa Rica", "accosted",
                     "a shortcut that stays open requires exotic matter"):
            with self.subTest(text):
                self.assertNotRegex(text, COST)


class HistoryShape(unittest.TestCase):
    """Every history has at least five paragraphs of three to six sentences each."""

    @staticmethod
    def metric(history, metric_id="x"):
        return {"id": metric_id, "name": "X", "short_name": "X", "tags": ["t"], "history": history}

    @staticmethod
    def paragraph(count, word="Physics"):
        return " ".join(f"{word} happened {n}." for n in range(count))

    def history(self, *counts):
        return "¶".join(self.paragraph(count) for count in counts)

    def test_sentences_are_counted_as_a_reader_counts_them(self):
        cases = {
            "J. Robert Oppenheimer met D. M. Chitre. They talked.": 2,
            'He called it "a simple example." Then he moved on.': 2,
            "The horizons sit at $$r_\\pm = \\frac{r_s}{2}.$$ These coalesce.": 2,
            "The value $0.378$ is small. So is $1.5\\ell$.": 2,
            "It is cited [kasner1921]. Relativists still read it.": 2,
            "Or would it? It would!": 2,
            "It was 1921. 1922 came next.": 2,
            "Cartan was there. Élie Cartan, that is.": 2,
            "A single sentence with no stop at the end": 1,
        }
        for text, count in cases.items():
            self.assertEqual(len(build.sentences(text)), count, text)

    def test_the_rule_holds_at_its_edges(self):
        self.assertIsNone(build.check_prose_shape([self.metric(self.history(3, 6, 3, 6, 3))]))

    def test_a_history_of_four_paragraphs_is_refused(self):
        with self.assertRaises(build.DataError) as raised:
            build.check_prose_shape([self.metric(self.history(4, 4, 4, 4))])
        self.assertIn("4 paragraphs", str(raised.exception))

    def test_a_paragraph_too_short_or_too_long_is_refused_by_its_place(self):
        for counts, place in (((4, 4, 2, 4, 4), 3), ((4, 4, 4, 4, 7), 5)):
            with self.assertRaises(build.DataError) as raised:
                build.check_prose_shape([self.metric(self.history(*counts))])
            self.assertIn(f"paragraph {place}", str(raised.exception))

    def test_the_longest_paragraph_may_be_at_most_twice_the_shortest(self):
        with mock.patch.object(build, "PARAGRAPH_SENTENCES", (1, 10)):
            self.assertIsNone(build.check_prose_shape([self.metric(self.history(2, 4, 3, 4, 2))]))
            with self.assertRaises(build.DataError) as raised:
                build.check_prose_shape([self.metric(self.history(2, 5, 3, 4, 2))])
        self.assertIn("more than 2 times", str(raised.exception))

    def test_a_table_is_not_a_paragraph(self):
        history = self.history(3, 3, 3, 3) + "¶TABLE:: a | b ;; 1 | 2¶" + self.paragraph(3)
        self.assertEqual(build.paragraph_shape(history), [3, 3, 3, 3, 3])
        self.assertIsNone(build.check_prose_shape([self.metric(history)]))

    def test_every_history_out_of_shape_is_named_at_once(self):
        with self.assertRaises(build.DataError) as raised:
            build.check_prose_shape([self.metric(self.history(3, 3), "one"),
                                       self.metric(self.history(4, 4, 4, 4, 4), "fine"),
                                       self.metric(self.history(9, 3, 3, 3, 3), "two")])
        self.assertIn("one.json", str(raised.exception))
        self.assertIn("two.json", str(raised.exception))
        self.assertNotIn("fine.json", str(raised.exception))

    def test_a_history_out_of_shape_leaves_both_published_files_alone(self):
        before = {path: path.read_text(encoding="utf-8") for path in (build.INDEX_FILE, build.REFERENCES_FILE)}
        broken = build.load_metrics()
        broken[0] = dict(broken[0], history=self.history(3, 3, 3))
        for argv in (["--check"], []):
            with mock.patch.object(build, "load_metrics", return_value=broken), \
                    contextlib.redirect_stderr(io.StringIO()):
                self.assertEqual(build.main(argv), 2)
        for path, text in before.items():
            self.assertEqual(path.read_text(encoding="utf-8"), text)


def page_rules():
    """Every rule of the spacetimes page's stylesheets, in order, as (selectors, declarations,
    print), print being whether the rule sits inside an @media print block."""
    page = (build.ROOT / "_layouts" / "mfs.html").read_text(encoding="utf-8")
    rules = []
    for sheet in re.findall(r"<style>(.*?)</style>", page, re.S):
        sheet = re.sub(r"/\*.*?\*/", "", sheet, flags=re.S)
        blocks, start = [], 0
        for i, char in enumerate(sheet):
            if char == "{":
                blocks.append(sheet[start:i].strip())
                start = i + 1
            elif char == "}":
                selector = blocks.pop()
                if not selector.startswith("@"):
                    declarations = dict(
                        (name.strip(), value.strip()) for name, _, value in
                        (d.partition(":") for d in sheet[start:i].split(";") if ":" in d))
                    rules.append(([s.strip() for s in selector.split(",")], declarations,
                                  any(b.startswith("@media print") for b in blocks)))
                start = i + 1
    return page, rules


@page_source
class WrittenAreas(unittest.TestCase):
    """Every written area of the spacetimes page is set in the font of the history's prose,
    as the captain asked on 28 September 2026: "every written area must use the same font as
    the History section." Bianchi's list of the nine types had come up in the page's serif
    among the history's monospaced prose, and the references, the placeholders, the loading
    line and the word standing in for a vanishing tensor still did."""

    HISTORY = "#mfs-content-panel .mfs-history p"
    # Each written area, by the rule that sets it. Headings, the choices, the list of names
    # and the words inside a drawing are not among them and keep their own fonts.
    WRITTEN = (
        HISTORY,
        "#mfs-content-panel .mfs-table",
        "#mfs-content-panel .mfs-convention p",
        "#mfs-content-panel .mfs-nr-caption p",
        "#mfs-content-panel .mfs-nr-legend",
        "#mfs-content-panel .mfs-nr-note",
        "#mfs-content-panel .mfs-cd-restriction",
        "#mfs-content-panel .mfs-reference",
        "#mfs-content-panel .mfs-zero",
        "#mfs-content-panel .mfs-placeholder",
        "#mfs-loading",
        "#mfs-coffee-text",
    )
    # What the page's script sets as prose, in the classes wireProse() watches, and the
    # written area that holds each.
    PROSE_HELD_BY = {
        ".mfs-history": "#mfs-content-panel .mfs-history p",
        ".mfs-convention": "#mfs-content-panel .mfs-convention p",
        ".mfs-nr-caption": "#mfs-content-panel .mfs-nr-caption p",
        ".mfs-nr-note": "#mfs-content-panel .mfs-nr-note",
        ".mfs-nr-itext": "#mfs-content-panel .mfs-nr-legend",
        ".mfs-cd-restriction": "#mfs-content-panel .mfs-cd-restriction",
    }

    @classmethod
    def setUpClass(cls):
        cls.page, cls.rules = page_rules()
        cls.tokens = {re.findall(r"[.#][\w-]+", selector)[-1] for selector in cls.WRITTEN} | {".mfs-history"}

    def resolved(self, selector, prop, printed=False):
        """What the last rule naming selector exactly gives prop, an !important one winning."""
        found = [(d[prop].endswith("!important"), n) for n, (selectors, d, p) in enumerate(self.rules)
                 if selector in selectors and prop in d and p == printed]
        if not found:
            return None
        return self.rules[max(found)[1]][1][prop].replace("!important", "").strip()

    def concerns_written_area(self, selector):
        return any(re.search(re.escape(token) + r"(?![\w-])", selector) for token in self.tokens)

    def test_the_history_is_source_code_pro_regular_and_upright(self):
        self.assertEqual(self.resolved(self.HISTORY, "font-family"), "'Source Code Pro', monospace")
        self.assertIn(self.resolved(self.HISTORY, "font-weight"), (None, "400", "normal"))
        self.assertIn(self.resolved(self.HISTORY, "font-style"), (None, "normal"))

    def test_every_written_area_is_set_in_the_history_font(self):
        family = self.resolved(self.HISTORY, "font-family")
        for selector in self.WRITTEN:
            with self.subTest(selector):
                self.assertEqual(self.resolved(selector, "font-family"), family)
                self.assertIn(self.resolved(selector, "font-weight"), (None, "400", "normal"))
                self.assertIn(self.resolved(selector, "font-style"), (None, "normal"))

    def test_every_written_area_is_set_at_the_history_size(self):
        # The note beside the coffee button stands in the list's column and keeps its size.
        self.assertEqual(self.resolved("#mfs-content-panel .mfs-history", "font-size"), "var(--mfs-prose)")
        self.assertEqual(self.resolved(self.HISTORY, "font-size"), "1em")
        for selector in set(self.WRITTEN) - {self.HISTORY, "#mfs-coffee-text"}:
            with self.subTest(selector):
                self.assertEqual(self.resolved(selector, "font-size"), "var(--mfs-prose)")

    def test_no_rule_gives_a_written_area_another_font(self):
        # A rule that reaches into a written area, such as the table's cells, the legend's
        # entries or the phone's sizes, may not name a family, a size or a style of its own; on
        # paper, where the copy is set whole in the print root's font at its 12pt, it may not
        # name a family or a size at all.
        family = self.resolved(self.HISTORY, "font-family")
        for selectors, declarations, printed in self.rules:
            for selector in filter(self.concerns_written_area, selectors):
                with self.subTest(selector, printed=printed):
                    named = declarations.get("font-family", "").replace("!important", "").strip()
                    self.assertIn(named, ("",) if printed else ("", family))
                    size = declarations.get("font-size", "").replace("!important", "").strip()
                    self.assertIn(size, ("",) if printed else ("", "var(--mfs-prose)", "1em"))
                    if selector.split()[-1] in self.tokens or selector.split()[-1] == "p":
                        self.assertNotEqual(declarations.get("font-style"), "italic")

    def test_no_written_area_carries_a_font_in_its_own_markup(self):
        for tag in re.findall(r"<[a-z]+\b[^>]*>", self.page):
            names = re.findall(r'(?:id|class)="([^"]*)"', tag)
            if any(("." + n in self.tokens or "#" + n in self.tokens) for ns in names for n in ns.split()):
                self.assertNotIn("font-family", tag, tag[:80])
        # The script sets no font of its own on anything; fitWords() only reads them.
        self.assertIsNone(re.search(r"\.style\.fontFamily\s*=|cssText[^\n]*font-family", self.page),
                          "the page's script sets a font")

    def test_every_prose_the_page_fits_is_a_written_area(self):
        watched = re.search(r"var PROSE = '([^']*)';", self.page)
        self.assertIsNotNone(watched, "PROSE no longer names the prose the page fits")
        self.assertRegex(self.page, r"root\.querySelectorAll\(PROSE\)\.forEach\(function\(block\) \{ _proseObserver",
                         "wireProse() no longer observes PROSE")
        for cls in (c.strip() for c in watched.group(1).split(",")):
            with self.subTest(cls):
                self.assertIn(cls, self.PROSE_HELD_BY, f"{cls} is prose; add its rule to WRITTEN")
                self.assertIn(self.PROSE_HELD_BY[cls], self.WRITTEN)

    def test_the_page_serves_the_history_font_itself(self):
        # The site's own declarations in custom.css are not loaded on this page, and without
        # its own a reader who had not installed Source Code Pro read the browser's monospace.
        self.assertIn("/assets/css/mfs.css", self.page)
        sheet = build.ROOT / "assets" / "css" / "mfs.css"
        faces = re.findall(r"@font-face\s*\{([^}]*)\}", sheet.read_text(encoding="utf-8"))
        styles = set()
        for face in faces:
            if "'Source Code Pro'" not in face:
                continue
            source = re.search(r"url\('([^']+)'\)", face).group(1)
            self.assertTrue((sheet.parent / source).resolve().is_file(), source)
            self.assertIn("font-weight: 100 900", face)
            styles.add(re.search(r"font-style:\s*(\w+)", face).group(1))
        self.assertEqual(styles, {"normal", "italic"})

    def test_the_print_copy_sets_its_prose_in_the_print_root_font(self):
        self.assertEqual(self.resolved("#mfs-print-root", "font-family", printed=True), "'EB Garamond', serif")


@page_source
class DividedWords(unittest.TestCase):
    """hyphenateWords() in _layouts/mfs.html divides a word of the prose wider than its line, and
    leaves whole all mathematics MathJax has yet to set, which MathJax finds only with both its
    delimiters in one piece of text: on a phone at three times the usual text, the signature
    (-,+,+,+) of the first spacetime opened, set as MathJax started, stood as TeX."""

    def setUp(self):
        page = (build.ROOT / "_layouts" / "mfs.html").read_text(encoding="utf-8")
        self.body = page[page.index("function hyphenateWords("):page.index("function buildMetricMatrix(")]
        found = re.search(r"math = /(.+?)/g,", self.body)
        self.assertIsNotNone(found, "hyphenateWords() looks for no mathematics")
        self.math = re.compile(found.group(1))
        # Every pair of delimiters MathJax is told to read, as the page's configuration spells them.
        self.pairs = [(a.replace("\\\\", "\\"), b.replace("\\\\", "\\"))
                      for line in re.findall(r"(?:inlineMath|displayMath): \[(.*)\]", page)
                      for a, b in re.findall(r"\['([^']*)', '([^']*)'\]", line)]
        self.assertEqual(len(self.pairs), 4)

    def test_mathematics_in_every_delimiter_is_found_whole(self):
        for opening, closing in self.pairs:
            tex = opening + "(-,+,+,+)" + closing
            with self.subTest(tex):
                spans = [m.span() for m in self.math.finditer("the signature " + tex + " holds")]
                self.assertEqual(spans, [(14, 14 + len(tex))])

    def test_no_word_inside_mathematics_is_divided(self):
        self.assertRegex(self.body, r"if \(tex\.some\(function\(t\) \{ return at < t\[1\] && to > t\[0\]; \}\)\) continue;")


@page_source
class FixedSizes(unittest.TestCase):
    """The spacetime's name, the word above it and every button keep one size in px at any
    text size, the size each had before the prose grew with the reader's text, and the prose's
    headings grow with it, as the captain asked on 29 September 2026: "When I said I wanted
    the title and subtitle fonts to increase along with the text, I meant the font in the
    actual prose, not in the buttons." Each fixed size is the desktop's and then the phone's."""

    FIXED = {
        "#mfs-content-panel .mfs-title": ("min(40px, var(--mfs-fit, 40px))", "min(26px, var(--mfs-fit, 26px))"),
        "#mfs-content-panel .mfs-label": ("min(21px, var(--mfs-fit, 21px))", "min(14px, var(--mfs-fit, 14px))"),
        "#mfs-content-panel #mfs-print-btn": ("15px", "12px"),
        "#mfs-content-panel .mfs-choice": ("min(15px, var(--mfs-fit, 15px))", "min(12px, var(--mfs-fit, 12px))"),
        ".mfs-result": ("min(18px, var(--mfs-fit, 18px))",),
        ".mfs-toc-entry": ("18px",),
        ".mfs-toc-back": ("18px",),
    }
    HEADINGS = "#mfs-content-panel .mfs-section-label"

    @classmethod
    def setUpClass(cls):
        cls.page, cls.rules = page_rules()

    def sizes(self, selector, prop="font-size"):
        return tuple(d[prop].replace("!important", "").strip() for selectors, d, printed in self.rules
                     if selector in selectors and prop in d and not printed)

    def test_the_name_and_every_button_keep_their_size_at_any_text_size(self):
        for selector, sizes in self.FIXED.items():
            with self.subTest(selector):
                self.assertEqual(self.sizes(selector), sizes)

    def test_nothing_on_the_screen_sizes_a_button_or_the_name_by_the_text(self):
        # A rule that reaches a button or the name, a hover or a phone's among them, may give
        # it no size or spacing that follows the reader's text.
        names = [re.findall(r"[.#][\w-]+", s)[-1] for s in self.FIXED]
        for selectors, declarations, printed in self.rules:
            for selector in selectors:
                # The element a selector sizes is named in its last compound, as in .mfs-choice:hover.
                target = re.findall(r"[.#][\w-]+", selector.split()[-1]) if selector.split() else []
                if printed or not set(target) & set(names):
                    continue
                for prop in ("font-size", "padding", "margin", "gap", "height", "width"):
                    with self.subTest(selector, prop=prop):
                        self.assertNotRegex(declarations.get(prop, ""), r"rem\b|--mfs-prose|--mfs-text")
        for selector in ("#mfs-content-panel .mfs-choices",):
            self.assertEqual(self.sizes(selector, "gap"), ("6px",))

    def test_the_reset_is_drawn_at_the_usual_text_size(self):
        self.assertTrue(re.search(r"drawing\.parentElement\.getBoundingClientRect\(\)\.width / cdFrame\(it\.figure\)\.vw",
                                  self.page), "the reset is not drawn at the frame's scale")
        self.assertEqual(self.sizes("#mfs-content-panel .mfs-turn-reset"), ("calc(13 * var(--em-u, 1px))",))

    def test_the_headings_of_the_prose_grow_with_the_text(self):
        sizes = self.sizes(self.HEADINGS)
        self.assertEqual(sizes, ("min(1.3125rem, var(--mfs-fit, 1.3125rem))", "min(0.875rem, var(--mfs-fit, 0.875rem))"))
        for selector in WrittenAreas.WRITTEN[1:]:
            for size in self.sizes(selector):
                with self.subTest(selector, size=size):
                    self.assertRegex(size, r"rem\b|--mfs-prose|^1em$")


class ReadableDrawings(unittest.TestCase):
    """Every word of every drawing is white and no smaller than the caption beside it, as the
    captain asked on 30 September 2026: "Make all diagram label text white, matching the
    spacetime and conformal diagrams" and "Set one readable minimum for every diagram kind:
    tick numbers and axis labels no smaller than the page's caption text at the reader's size".
    The page sets a spacetime diagram's words and a drawing's labels at --mfs-prose, the
    caption's size, and draws a conformal diagram, a figure or an embedding diagram no narrower
    than keeps a label at CD_LAB of its 628 units, which the generators compose them at.
    `node _tools/page_timing.mjs` holds the page as drawn to the same rule."""

    # The rules that set a drawing's words, and the elements a word of a drawing is.
    GROUNDS = (".mfs-nr-figure", ".mfs-cd-drawing")
    WORDS = (".mfs-nr-figure", ".nr-tick", ".nr-name", ".nr-ref", ".mfs-slice-label", ".mfs-nr-readout", ".cv-word",
             ".mfs-cd-labels", ".cd-at", ".cd-lab-small", ".cd-lab-coord", ".cd-lab-region", ".cd-lab-lab")

    @classmethod
    def setUpClass(cls):
        cls.page, cls.rules = page_rules()

    def screen(self):
        for selectors, declarations, printed in self.rules:
            if printed:
                continue
            for selector in selectors:
                last = selector.split()[-1] if selector.split() else ""
                if any(word in last for word in self.WORDS):
                    yield selector, declarations

    def test_every_word_of_a_drawing_is_white(self):
        coloured = [(selector, d["color"]) for selector, d in self.screen() if "color" in d]
        for base in (".mfs-nr-figure", ".mfs-cd-labels", ".mfs-slice-label", "#mfs-content-panel .mfs-nr-readout"):
            self.assertIn((base, "var(--text-bright)"), coloured)
        for selector, colour in coloured:
            with self.subTest(selector):
                self.assertEqual(colour, "var(--text-bright)")

    def test_every_word_of_a_drawing_is_set_at_the_caption_size(self):
        for ground in self.GROUNDS:
            sizes = [d["font-size"] for selectors, d, printed in self.rules
                     if ground in selectors and "font-size" in d and not printed]
            self.assertEqual(sizes, ["var(--mfs-prose)"], ground)
        for selector, declarations in self.screen():
            size = declarations.get("font-size")
            if size is None or selector in self.GROUNDS:
                continue
            with self.subTest(selector, size=size):
                self.assertRegex(size, r"^(inherit|1(\.\d+)?em)$")
        self.assertNotIn("100cqw * ", self.page)

    def test_the_label_size_a_drawing_is_composed_at_is_one_size_everywhere(self):
        page = int(re.search(r"var CD_LAB = (\d+);", self.page).group(1))
        derivations = build.ROOT / "_tools" / "derivations"
        conformal = (derivations / "conformal.py").read_text(encoding="utf-8")
        embedding = (derivations / "embedding.py").read_text(encoding="utf-8")
        projections = (derivations / "projections.py").read_text(encoding="utf-8")
        turn = (build.ROOT / "MFS" / "assets" / "turn.js").read_text(encoding="utf-8")
        found = {
            "conformal.py": re.search(r"^CD_LABEL_SIZE = \{([^}]*)\}", conformal, re.M).group(1),
            "embedding.py": re.search(r"^LAB = \{([^}]*)\}", embedding, re.M).group(1),
            "projections.py": "lab: " + re.search(r"^CD_W, CD_LAB = 560, (\d+)$", projections, re.M).group(1),
            "turn.js": re.search(r"var LAB = \{([^}]*)\}", turn).group(1),
        }
        for name, sizes in found.items():
            with self.subTest(name):
                self.assertEqual({int(v) for v in re.findall(r":\s*(\d+)", sizes)}, {page})

    def test_the_page_and_the_generators_give_every_label_one_box(self):
        # The margin of 1.5 em is kept from the box cdLabelSize() gives a label, and slices.py's
        # label_size() places and spaces labels by the same box, each character MathJax sets
        # wider than 0.7 em at its own width: "$4m$" stood 1.4 em from the edge of
        # Majumdar-Papapetrou's embedding diagram while an m was counted 0.7, and "$X = X_s$"
        # 1.45 em from the edge of Whittaker's sphere while an X was. Both are run on every label
        # of every drawing, the page's in Node. `node _tools/label_sizes.mjs` holds the box to
        # the one MathJax sets.
        if shutil.which("node") is None:
            self.skipTest("node is not installed")
        sys.path.insert(0, str(build.ROOT / "_tools" / "derivations"))
        import slices
        start = self.page.index("      var CD_RELATION = ")
        end = self.page.index("      window._mfsLabelSize = cdLabelSize;")
        data = build.ROOT / "MFS" / "assets" / "data"
        texts = set()
        for folder in ("conformal", "embedding", "diagrams"):
            for path in sorted((data / folder).glob("*.json")):
                def walk(node):
                    if isinstance(node, dict):
                        for key, value in node.items():
                            if key in ("text", "label") and isinstance(value, str):
                                texts.add(value)
                            walk(value)
                    elif isinstance(node, list):
                        for value in node:
                            walk(value)
                walk(json.loads(path.read_text(encoding="utf-8")))
        texts = sorted(texts)
        self.assertGreater(len(texts), 1000)
        script = (self.page[start:end] + "const texts = JSON.parse(require('fs').readFileSync(0, 'utf8'));"
                  "console.log(JSON.stringify(texts.map(cdLabelSize)));")
        run = subprocess.run(["node", "-e", script], input=json.dumps(texts), capture_output=True, text=True, check=True)
        for text, page in zip(texts, json.loads(run.stdout)):
            with self.subTest(text):
                self.assertEqual([round(v, 9) for v in page], [round(v, 9) for v in slices.label_size(text)])


class Anchors(unittest.TestCase):
    """Every label of every drawing stands by an anchor the page knows, CD_ANCHOR in
    _layouts/mfs.html and ANCHOR in MFS/assets/turn.js, which name the same nine. Until
    3 October 2026 two labels of Schrödinger's spacetime stood by "bc", which neither knows, and
    its global chart raised an error and was never drawn."""

    def test_every_label_stands_by_an_anchor_the_page_knows(self):
        page = (build.ROOT / "_layouts" / "mfs.html").read_text(encoding="utf-8")
        turn = (build.ROOT / "MFS" / "assets" / "turn.js").read_text(encoding="utf-8")
        known = set(re.findall(r"(\w+): \[", re.search(r"var CD_ANCHOR = \{(.*?)\};", page, re.S).group(1)))
        self.assertEqual(known, {"c", "l", "r", "t", "b", "tl", "tr", "bl", "br"})
        self.assertEqual(set(re.findall(r"(\w+): \[", re.search(r"var ANCHOR = \{(.*?)\};", turn, re.S).group(1))), known)

        def anchors(value, where):
            if isinstance(value, dict):
                if "anchor" in value:
                    yield where, value["anchor"]
                for key, item in value.items():
                    yield from anchors(item, f"{where}/{key}")
            elif isinstance(value, list):
                for i, item in enumerate(value):
                    yield from anchors(item, f"{where}[{i}]")

        data = build.ROOT / "MFS" / "assets" / "data"
        for folder in ("diagrams", "conformal", "embedding"):
            for path in sorted((data / folder).glob("*.json")):
                for where, anchor in anchors(json.loads(path.read_text(encoding="utf-8")), ""):
                    with self.subTest(f"{folder}/{path.name}{where}"):
                        self.assertIn(anchor, known)


@page_source
class NoGlow(unittest.TestCase):
    """Nothing on the spacetimes page glows, and a button chosen or pressed turns pink, as the
    captain asked on 29 September 2026: "Get rid of the glow on the pressed button in MFS.
    Only turn it pink. I already said to get rid of the glow earlier. No glow anywhere."
    `node _tools/page_timing.mjs` holds the page as drawn to the same rule."""

    # Each way a stylesheet, a script or a drawing can make light, and what it looks like.
    GLOW = (
        (r"(?i)(?<![\w-])(?:text|box)-shadow\s*:(?!\s*none\b)", "a shadow"),
        (r"\.style\.(?:textShadow|boxShadow)\s*=(?!\s*['\"]none['\"])", "a shadow set by script"),
        (r"(?i)drop-shadow\s*\(", "a drop shadow"),
        (r"(?i)(?<![\w-])filter\s*:[^;}\"]*blur\s*\(", "a blur"),
        (r"\.style\.(?:filter|webkitFilter)\s*=\s*[^;]*(?:blur|drop-shadow)", "a blur set by script"),
        (r"\bshadow(?:Blur|Color|OffsetX|OffsetY)\s*=", "a canvas shadow"),
        (r"globalCompositeOperation\s*=\s*['\"](?:lighter|screen|plus-lighter)", "light added on a canvas"),
        (r"(?i)<filter\b|\bfe(?:GaussianBlur|DropShadow|Morphology)\b", "an SVG filter"),
        (r"(?i)mix-blend-mode\s*:\s*(?:screen|lighten|plus-lighter|color-dodge)", "light blended in"),
        (r"(?i)(?:class(?:Name|List\.(?:add|toggle))?\s*[=(]\s*['\"][^'\"]*|[.#][\w-]*)(?:glow|neon|halo|bloom)",
         "a glow by name"),
    )

    @classmethod
    def setUpClass(cls):
        root = build.ROOT
        cls.sources = [root / "_layouts" / "mfs.html", root / "assets" / "css" / "mfs.css",
                       root / "assets" / "css" / "palette.css"]
        cls.sources += sorted((root / "MFS").glob("*.markdown")) + sorted((root / "MFS" / "assets").glob("*.js"))
        cls.page, cls.rules = page_rules()

    def glows(self, text):
        return [what for pattern, what in self.GLOW if re.search(pattern, text)]

    def test_no_source_of_the_page_makes_anything_glow(self):
        for path in self.sources:
            with self.subTest(str(path.relative_to(build.ROOT))):
                found = [(what, text.strip()[:80]) for text in path.read_text(encoding="utf-8").splitlines()
                         for what in self.glows(text)]
                self.assertEqual(found, [])

    def test_the_rules_catch_every_glow_the_page_carried(self):
        for text in ("box-shadow:0 0 2px var(--cyan),0 0 10px var(--cyan);",
                     "  box-shadow: 0 0 6px var(--cyan), 0 0 14px var(--cyan);",
                     "a:hover { text-shadow: 0 0 8px var(--cyan), 0 0 20px var(--cyan); }",
                     "span.style.textShadow = '0 0 6px ' + tealLight;",
                     "span.style.textShadow =",
                     "  filter: drop-shadow(1px 1px 0px rgba(0,0,0,0.9));",
                     "  filter: blur(6px) brightness(2);",
                     "ctx.shadowBlur  = glowStrength * 18;",
                     "ctx.shadowColor = 'rgba(' + tealRgb + ')';",
                     "ctx.globalCompositeOperation = 'lighter';",
                     '<filter id="g"><feGaussianBlur stdDeviation="3"/></filter>',
                     ".mfs-choice-glow { color: red; }",
                     "el.classList.add('neon');"):
            with self.subTest(text):
                self.assertTrue(self.glows(text))

    def test_the_rules_leave_plain_colour_and_the_print_reset_alone(self):
        for text in ("        text-shadow: none !important;",
                     "        box-shadow: none !important;",
                     "span.style.textShadow = 'none';",
                     "  backdrop-filter: blur(4px);",
                     ".mfs-choice.mfs-choice-on { color: var(--pink-light); border-color: var(--pink-light); }",
                     "var shadow = ring.filter(function(p) { return p.hidden; });",
                     "The Shadow of the Supermassive Black Hole"):
            with self.subTest(text):
                self.assertEqual(self.glows(text), [])

    def last(self, selector, prop):
        found = [d[prop] for selectors, d, printed in self.rules if selector in selectors and prop in d and not printed]
        return found[-1] if found else None

    def test_the_chosen_button_and_the_spacetime_shown_turn_pink(self):
        self.assertEqual(self.last("#mfs-content-panel .mfs-choice.mfs-choice-on", "color"), "var(--pink-light)")
        self.assertEqual(self.last("#mfs-content-panel .mfs-choice.mfs-choice-on", "border-color"), "var(--pink-light)")
        self.assertEqual(self.last(".mfs-result.mfs-result-active", "color"), "var(--pink-light)")
        # At once, so a press as short as a tap shows it whole.
        for selector in ("#mfs-content-panel .mfs-choice.mfs-choice-on", "#mfs-content-panel .mfs-choice:active",
                         "#mfs-content-panel #mfs-print-btn:active"):
            self.assertEqual(self.last(selector, "transition"), "none", selector)

    def test_every_button_turns_pink_while_it_is_pressed(self):
        # Each :active rule comes after the button's :hover, which a mouse holds while it presses.
        for button, has_border in (("#mfs-content-panel .mfs-choice", True), ("#mfs-content-panel #mfs-print-btn", True),
                                   ("#mfs-content-panel .mfs-turn-reset", True), (".mfs-result", False),
                                   (".mfs-toc-entry", False), (".mfs-toc-back", False)):
            with self.subTest(button):
                self.assertEqual(self.last(button + ":active", "color"), "var(--pink-light)")
                if has_border:
                    self.assertEqual(self.last(button + ":active", "border-color"), "var(--pink-light)")
                places = [n for n, (selectors, _, printed) in enumerate(self.rules) if not printed
                          for s in selectors if s in (button + ":hover", button + ":active")]
                pressed = max(n for n, (selectors, _, _) in enumerate(self.rules) if button + ":active" in selectors)
                self.assertEqual(max(places), pressed)


@page_source
class CyanNotTeal(unittest.TestCase):
    """The spacetimes page is drawn in VoidFlux's cyans, #29a3c3 and #37dfff, where the rest of
    the site keeps its teal, as the captain asked on 30 September 2026: mfs.css sets them for the
    page alone, and nothing the page loads names the teal, as a token, a hex or an rgb, but
    palette.css, which the whole site shares."""

    TEAL = r"(?i)--teal-|--bg-active-teal|#23bbad|#25d9c8|\b35,\s*187,\s*173\b|\b37,\s*217,\s*200\b"

    def test_the_page_names_no_teal(self):
        root = build.ROOT
        sources = [root / "_layouts" / "mfs.html", root / "assets" / "css" / "mfs.css"]
        sources += sorted((root / "MFS").glob("*.markdown")) + sorted((root / "MFS" / "assets").glob("*.js"))
        sources += sorted((root / "MFS" / "assets" / "data").rglob("*.json"))
        for path in sources:
            with self.subTest(str(path.relative_to(root))):
                self.assertIsNone(re.search(self.TEAL, path.read_text(encoding="utf-8")))

    def test_mfs_css_sets_the_cyans(self):
        sheet = (build.ROOT / "assets" / "css" / "mfs.css").read_text(encoding="utf-8")
        root = re.search(r":root\s*\{([^}]*)\}", sheet).group(1)
        tokens = dict(re.findall(r"(--[\w-]+)\s*:\s*([^;]+);", root))
        self.assertEqual(tokens["--mfs-cyan-dark"], "#29a3c3")
        self.assertEqual(tokens["--mfs-cyan-light"], "#37dfff")
        self.assertEqual(tokens["--mfs-cyan-dark-rgb"], "41, 163, 195")
        self.assertEqual(tokens["--mfs-cyan-light-rgb"], "55, 223, 255")


class ConventionShape(unittest.TestCase):
    """Every convention says only what a reader needs to read the mathematics: each chart has
    its own, the reader sees it with the text its spacetime shares across its charts as one
    paragraph of at most five sentences, and none of it says the words the captain named on
    29 September 2026."""

    @staticmethod
    def metric(shared="We keep factors of $c$ explicit.", charts=None, metric_id="x"):
        charts = {"one": "We use coordinates $(t, r)$."} if charts is None else charts
        coordinates = [{"id": chart_id, "convention": text} if text is not None else {"id": chart_id}
                       for chart_id, text in charts.items()]
        metric = dict(HistoryShape.metric("¶".join(HistoryShape.paragraph(3) for _ in range(5)), metric_id),
                      coordinates=coordinates)
        if shared is not None:
            metric["convention"] = shared
        return metric

    @staticmethod
    def sentences(count):
        return " ".join(f"Coordinate {n} is a length." for n in range(count))

    def test_the_reader_sees_the_charts_text_then_the_shared_text(self):
        metric = self.metric("We keep factors of $c$ explicit.",
                             {"a": "We use coordinates $(t, r)$.", "b": "We use coordinates $(u, r)$."})
        self.assertEqual(build.chart_conventions(metric),
                         [("a", "We use coordinates $(t, r)$. We keep factors of $c$ explicit."),
                          ("b", "We use coordinates $(u, r)$. We keep factors of $c$ explicit.")])

    def test_one_to_five_sentences_in_all_are_taken(self):
        for chart, shared in ((1, 0), (1, 1), (3, 2), (5, 0), (1, 4)):
            metric = self.metric(self.sentences(shared), {"one": self.sentences(chart)})
            self.assertIsNone(build.check_conventions([metric]))

    def test_six_sentences_in_all_are_refused_by_chart(self):
        metric = self.metric(self.sentences(2), {"short": self.sentences(3), "long": self.sentences(4)})
        with self.assertRaises(build.DataError) as raised:
            build.check_conventions([metric])
        self.assertIn("x.json: the convention of the chart 'long' has 6 sentences", str(raised.exception))
        self.assertNotIn("'short'", str(raised.exception))

    def test_the_words_the_captain_named_are_refused(self):
        for text in ("A thin, fast bubble costs so much.", "The cost of the throat is exotic matter.",
                     "It has no singularity, since nowhere does the geometry break.",
                     "This is the standing objection to the drive.",
                     "Its Weyl tensor is the Riemann tensor with the traces removed.",
                     "Thirty six slots of $C_{\\mu\\nu\\rho\\sigma}$ behave the same way."):
            with self.subTest(text):
                with self.assertRaises(build.DataError):
                    build.check_conventions([self.metric("We keep factors of $c$ explicit.", {"one": text})])
        for text in ("We use coordinates $(t, r)$, with $t$ carrying dimensions of time.",
                     "The Ricci tensor is the standard contraction, $R_{\\mu\\nu} = R^\\alpha{}_{\\mu\\alpha\\nu}$."):
            self.assertIsNone(build.check_conventions([self.metric(text)]))

    def test_a_chart_without_its_own_convention_is_refused(self):
        for text in (None, "", "  "):
            with self.assertRaises(build.DataError) as raised:
                build.check_conventions([self.metric(charts={"one": "We use $t$.", "bare": text})])
            self.assertIn("the chart 'bare' has no convention of its own", str(raised.exception))

    def test_a_spacetime_keeps_its_shared_convention_even_when_empty(self):
        self.assertIsNone(build.check_conventions([self.metric("")]))
        with self.assertRaises(build.DataError) as raised:
            build.check_conventions([self.metric(None)])
        self.assertIn("x.json has no convention shared by its charts", str(raised.exception))

    def test_a_convention_is_one_paragraph(self):
        with self.assertRaises(build.DataError) as raised:
            build.check_conventions([self.metric("We keep $c$.¶We keep $G$.")])
        self.assertIn("runs to more than one paragraph", str(raised.exception))

    def test_every_convention_out_of_rule_is_named_at_once(self):
        with self.assertRaises(build.DataError) as raised:
            build.check_conventions([self.metric(self.sentences(6), metric_id="one"), self.metric(metric_id="fine"),
                                     self.metric(charts={"a": "The standing objection."}, metric_id="two")])
        self.assertIn("one.json: the convention", str(raised.exception))
        self.assertIn("two.json: the convention", str(raised.exception))
        self.assertNotIn("fine.json", str(raised.exception))

    def test_a_convention_out_of_rule_leaves_both_published_files_alone(self):
        before = {path: path.read_text(encoding="utf-8") for path in (build.INDEX_FILE, build.REFERENCES_FILE)}
        broken = build.load_metrics()
        broken[0] = dict(broken[0], convention=self.sentences(8))
        for argv in (["--check"], []):
            with mock.patch.object(build, "load_metrics", return_value=broken), \
                    contextlib.redirect_stderr(io.StringIO()):
                self.assertEqual(build.main(argv), 2)
        for path, text in before.items():
            self.assertEqual(path.read_text(encoding="utf-8"), text)

    def test_every_break_of_a_history_stands_where_a_sentence_ends(self):
        # A break takes the place of the space between two sentences, so turning every break
        # back into a space gives the prose as it reads unbroken, and splitting at the breaks
        # gives each paragraph with nothing to trim. A paragraph that leads into a table may
        # end on a colon, since the table is not prose.
        for metric in build.load_metrics():
            text = metric.get("history")
            if not text:
                continue
            with self.subTest(f"{metric['id']}.json history"):
                self.assertNotRegex(text, r"\s¶|¶\s|^¶|¶$|¶¶")
                paragraphs = text.split("¶")
                for before, after in zip(paragraphs, paragraphs[1:]):
                    if "TABLE::" in (before[:7], after[:7]):
                        continue
                    self.assertEqual(len(build.sentences(f"{before} {after}")),
                                     len(build.sentences(before)) + len(build.sentences(after)),
                                     f"a break inside a sentence, before {after[:40]!r}")


class Diagrams(unittest.TestCase):
    """A diagram is drawn from what its metric publishes, and stops being published when that changes."""

    def setUp(self):
        self.metrics = build.load_metrics()
        self.diagrams = diagram_files()
        self.named = NamedFiles(build.DIAGRAMS_DIR)
        # The file the rule tests spoil, the first by name whichever spacetimes this run checks.
        self.first = min(p.stem for p in build.DIAGRAMS_DIR.glob("*.json"))
        found(self, self.diagrams, "diagram file")

    def assert_some_flat_view(self, views):
        """A test of the flat views ran on at least one, or the spacetimes this run checks draw none,
        as a spacetime whose spacetime diagrams are all projections does."""
        if not views and CHOSEN and not any(d["systems"] for d in self.diagrams.values()):
            self.skipTest("the spacetimes this run checks draw no flat view")
        self.assertTrue(views)

    def test_every_view_was_drawn_from_what_its_metric_publishes(self):
        by_id = {m["id"]: m for m in self.metrics}
        views = 0
        for metric_id, diagram in self.diagrams.items():
            systems = {s["id"]: s for s in by_id[metric_id]["coordinates"]}
            for system_id, drawn in diagram["systems"].items():
                for view in drawn:
                    views += 1
                    self.assertEqual(
                        build.diagram_source_version(systems[system_id], view["source"]["fields"]),
                        view["source"]["version"], f"{metric_id}/{system_id}/{view['id']}")
        self.assertEqual({name for name in build.load_diagrams(self.metrics) if chosen(name)}, set(self.diagrams))
        self.assert_some_flat_view(views)

    def test_a_changed_component_is_refused_and_leaves_the_files_alone(self):
        metric_id, system_id, changed = self.metric_with_a_changed_component()
        before = {path: path.read_text(encoding="utf-8") for path in (build.INDEX_FILE, build.REFERENCES_FILE)}
        for argv in (["--check"], []):
            stderr = io.StringIO()
            with mock.patch.object(build, "load_metrics", return_value=changed), \
                    contextlib.redirect_stderr(stderr):
                self.assertEqual(build.main(argv), 2)
            self.assertIn(f"diagrams/{metric_id}.json", stderr.getvalue())
            self.assertIn(system_id, stderr.getvalue())
        for path, text in before.items():
            self.assertEqual(path.read_text(encoding="utf-8"), text)

    def test_a_change_that_draws_nothing_leaves_the_diagrams_standing(self):
        metric_id = self.first
        changed = copy.deepcopy(self.metrics)
        for metric in changed:
            if metric["id"] == metric_id:
                metric["history"] = metric.get("history", "") + " More."
                for system in metric["coordinates"]:
                    for parameter in system.get("parameters", []):
                        parameter["description"] = parameter.get("description", "") + " Reworded."
        self.assertIn(metric_id, build.load_diagrams(changed))

    def test_the_principal_rays_are_tied_to_the_weyl_tensor_and_the_christoffel_symbols(self):
        """Kerr's principal null rays are found from its Weyl tensor and checked against its
        Christoffel symbols, so a change to either stops its diagrams, as a changed metric does."""
        for field in ("weyl_tensor", "christoffel"):
            changed = copy.deepcopy(self.metrics)
            kerr = next(m for m in changed if m["id"] == "kerr")
            system = next(s for s in kerr["coordinates"] if s["id"] == "boyer_lindquist")
            variant = system[field]["variants"][system[field]["default"]]
            variant["nonzero"][0]["value"] = "2" + variant["nonzero"][0]["value"]
            with self.assertRaises(build.DataError) as raised:
                build.load_diagrams(changed)
            self.assertIn("diagrams/kerr.json", str(raised.exception), field)
            self.assertIn("'principal'", str(raised.exception), field)

    def test_every_axis_is_labelled_as_tex_and_ticked_inside_its_box(self):
        """The page and the application both set these labels as TeX and draw no tick of their own."""
        views = 0
        for metric_id, diagram in self.diagrams.items():
            for system_id, drawn in diagram["systems"].items():
                for view in drawn:
                    views += 1
                    where = f"{metric_id}/{system_id}/{view['id']}"
                    references = [mk["label"] for mk in view["markers"] if mk["kind"] == "reference"]
                    for label in [view["xlabel"], view["ylabel"], *references]:
                        self.assertRegex(label, r"^\$[^$]+\$$", where)
                    x0, x1, y0, y1 = view["box"]
                    for axis, lo, hi in (("x", -x1 if view["mirror"] else x0, x1), ("y", y0, y1)):
                        at = [tick["at"] for tick in view["ticks"][axis]]
                        self.assertGreaterEqual(len(at), 2, f"{where} {axis}")
                        self.assertEqual(at, sorted(at), f"{where} {axis}")
                        self.assertTrue(lo <= at[0] and at[-1] <= hi, f"{where} {axis}")
                        for tick in view["ticks"][axis]:
                            self.assertRegex(tick["label"], r"^\$-?\d+(\.\d+)?\$$", where)
                            self.assertEqual(float(tick["label"][1:-1]), tick["at"], where)
        self.assert_some_flat_view(views)

    def test_a_diagram_of_a_system_its_metric_lacks_is_refused(self):
        metric = {"id": "x", "name": "X", "short_name": "X", "tags": ["t"], "coordinates": [{"id": "a"}]}
        diagram = {"metric": "x", "systems": {"b": []}}
        with self.assertRaises(build.DataError) as raised:
            self.load_diagram_folder({"x.json": diagram}, [metric])
        self.assertIn("'b'", str(raised.exception))

    def test_a_diagram_without_a_metric_is_refused(self):
        with self.assertRaises(build.DataError) as raised:
            self.load_diagram_folder({"ghost.json": {"metric": "ghost", "systems": {}}}, [])
        self.assertIn("ghost", str(raised.exception))

    def test_a_diagram_file_that_draws_nothing_is_refused(self):
        """A spacetime with no diagram has no file, rather than one the page would find empty."""
        metric = {"id": "x", "name": "X", "short_name": "X", "tags": ["t"], "coordinates": [{"id": "a"}]}
        for diagram in ({"metric": "x", "systems": {}}, {"metric": "x", "systems": {"a": []}, "projections": {}}):
            with self.assertRaises(build.DataError) as raised:
                self.load_diagram_folder({"x.json": diagram}, [metric])
            self.assertIn("draws nothing", str(raised.exception))

    def metric_with_a_changed_component(self):
        metric_id = self.first
        system_id = next(iter(self.named[metric_id]["systems"]))
        changed = copy.deepcopy(self.metrics)
        for metric in changed:
            if metric["id"] == metric_id:
                system = next(s for s in metric["coordinates"] if s["id"] == system_id)
                system["metric_components"][0]["value"] = "2" + system["metric_components"][0]["value"]
        return metric_id, system_id, changed

    def load_diagram_folder(self, files, metrics):
        with tempfile.TemporaryDirectory() as folder:
            for filename, diagram in files.items():
                (Path(folder) / filename).write_text(json.dumps(diagram), encoding="utf-8")
            with mock.patch.object(build, "DIAGRAMS_DIR", Path(folder)):
                return build.load_diagrams(metrics)


class Figures(unittest.TestCase):
    """A figure in three dimensions is drawn from what its metric publishes, inside its box,
    and names in its legend only what it draws."""

    def setUp(self):
        self.metrics = build.load_metrics()
        self.diagrams = diagram_files()
        self.figures = [(name, system, figure) for name, diagram in self.diagrams.items()
                        for system, figures in diagram.get("projections", {}).items() for figure in figures]
        found(self, self.figures, "figure in three dimensions")

    def test_every_figure_draws_inside_its_box_and_names_only_what_it_draws(self):
        for name, system, figure in self.figures:
            where = f"{name}/{system}/{figure['id']}"
            x0, x1, y0, y1 = figure["box"]
            self.assertTrue(x0 < x1 and y0 < y1, where)
            drawn = set()
            for layer in figure["layers"]:
                self.assertIn(layer["kind"], {"fill", "line", "point"}, where)
                drawn.add(layer["class"])
                for x, y in [layer["at"]] if layer["kind"] == "point" else layer["points"]:
                    self.assertTrue(x0 <= x <= x1 and y0 <= y <= y1, f"{where} {layer['class']} at {x}, {y}")
            for label in figure["labels"]:
                x, y = label["at"]
                self.assertTrue(x0 <= x <= x1 and y0 <= y <= y1, f"{where} label {label['text']}")
            for kind, cls, _ in figure["legend"]:
                self.assertIn(kind, {"fill", "line", "point", "cone"}, where)
                self.assertIn(cls, drawn, f"{where} legend {cls}")

    def test_no_label_of_a_figure_or_a_conformal_diagram_stands_over_another(self):
        # The generators keep every label clear of the others by slices.py's clear_labels() at
        # the size the page sets them at, which a label they have placed leaves where it is:
        # until 3 October 2026 the figures were not held to it, and the name $c\tau/r_0$ of the
        # wormhole time machine's proper times stood over the number 0 above it.
        sys.path.insert(0, str(build.ROOT / "_tools" / "derivations"))
        import slices
        drawings = [(f"{name}/{system}/{figure['id']}", figure) for name, system, figure in self.figures]
        drawings += [(f"conformal {name}/{view['id']}", view) for name, data in conformal_files().items()
                     for view in data["views"]]
        for where, drawing in drawings:
            x0, x1, y0, y1 = drawing["box"]
            scale = 560 / (x1 - x0)
            labels = copy.deepcopy(drawing["labels"])
            with self.subTest(where):
                slices.clear_labels(labels, lambda L, anchor, dx, dy: slices.label_box(
                    (L["at"][0] - x0) * scale + dx, (y1 - L["at"][1]) * scale + dy, L["text"], anchor, 21), where)
                self.assertEqual(labels, drawing["labels"])

    def test_every_class_a_figure_paints_is_styled_on_the_page_and_in_print(self):
        # A class the stylesheet does not know is painted as a black line on the dark page.
        page = (build.ROOT / "_layouts" / "mfs.html").read_text(encoding="utf-8")
        for name, system, figure in self.figures:
            for layer in figure["layers"]:
                selector = "pj-" + layer["class"] + ("-fill" if layer["kind"] == "fill" else "")
                self.assertRegex(page, r"\n    \." + re.escape(selector) + r"(?![\w-])[^{\n]*\{",
                                 f"{name}/{system}/{figure['id']} paints {selector}, which is not styled")
                self.assertRegex(page, r"\.mfs-print-body \." + re.escape(selector) + r"(?![\w-])",
                                 f"{name}/{system}/{figure['id']} paints {selector}, which is not styled in print")

    def test_every_text_is_tex_with_its_mathematics_closed(self):
        for name, diagram in self.diagrams.items():
            for field, value in diagram_prose(name, diagram):
                self.assertTrue(value.strip(), field)
                self.assertEqual(value.replace("\\$", "").count("$") % 2, 0, field)
                self.assertEqual(value.count("{"), value.count("}"), field)

    def test_a_changed_component_stops_a_figure(self):
        by_id = {m["id"]: m for m in self.metrics}
        for name, system, _ in self.figures:
            changed = copy.deepcopy(self.metrics)
            chart = next(s for s in next(m for m in changed if m["id"] == name)["coordinates"] if s["id"] == system)
            chart["metric_components"][0]["value"] = "2" + chart["metric_components"][0]["value"]
            with self.assertRaises(build.DataError) as raised:
                build.load_diagrams(changed)
            self.assertIn(f"diagrams/{name}.json", str(raised.exception))
            self.assertIn(system, str(raised.exception))
        self.assertTrue(by_id)


class ConformalDiagrams(unittest.TestCase):
    """A conformal diagram is of a whole spacetime, drawn from what its metrics publish, and
    stops being published when any of that changes."""

    SLICES = {"kerr", "kerr_newman", "cosmic_string"}

    def setUp(self):
        self.metrics = build.load_metrics()
        self.conformal = conformal_files()
        self.named = NamedFiles(build.CONFORMAL_DIR)
        found(self, self.conformal, "conformal diagram file")

    def test_every_file_was_drawn_from_what_its_metrics_publish(self):
        by_id = {m["id"]: m for m in self.metrics}
        for metric_id, data in self.conformal.items():
            self.assertTrue(data["source"], metric_id)
            for source in data["source"]:
                system = next(s for s in by_id[source["metric"]]["coordinates"] if s["id"] == source["system"])
                self.assertEqual(build.diagram_source_version(system, source["fields"]), source["version"],
                                 f"{metric_id} from {source['metric']}/{source['system']}")
        self.assertEqual({name for name in build.load_conformal(self.metrics) if chosen(name)}, set(self.conformal))

    def test_a_changed_component_is_refused_and_leaves_the_files_alone(self):
        # The interior Schwarzschild star is drawn with schwarzschild.json's exterior, so a
        # change there must stop it too, and not only schwarzschild's own diagram.
        changed = copy.deepcopy(self.metrics)
        for metric in changed:
            if metric["id"] == "schwarzschild":
                system = next(s for s in metric["coordinates"] if s["id"] == "spherical")
                system["metric_components"][0]["value"] = "2" + system["metric_components"][0]["value"]
        before = {path: path.read_text(encoding="utf-8") for path in (build.INDEX_FILE, build.REFERENCES_FILE)}
        for argv in (["--check"], []):
            stderr = io.StringIO()
            with mock.patch.object(build, "load_metrics", return_value=changed), \
                    mock.patch.object(build, "load_diagrams", return_value={}), \
                    contextlib.redirect_stderr(stderr):
                self.assertEqual(build.main(argv), 2)
            self.assertIn("conformal/", stderr.getvalue())
            self.assertIn("conformal.py --metric", stderr.getvalue())
        for path, text in before.items():
            self.assertEqual(path.read_text(encoding="utf-8"), text)
        with self.assertRaises(build.DataError) as raised:
            self.load_folder({"interior_schwarzschild.json": self.named["interior_schwarzschild"]}, changed)
        self.assertIn("schwarzschild.json", str(raised.exception))

    def test_r_equals_zero_and_the_region_inside_it_stand_apart(self):
        # A tower's r = 0 and the name of the region beside it, r < r_-, stand across one edge
        # of the drawing. Their boxes, 1.5 high and for r < r_- 3.8 wide in units of the 21 a
        # label is composed at, as cdLabelSize() in _layouts/mfs.html gives them, share no line
        # of the page or stand a label's size apart along it, the most a reader's text makes
        # them.
        size, seen = 21, 0
        for metric_id, data in self.conformal.items():
            for view in data["views"]:
                scale = 560 / (view["box"][1] - view["box"][0])
                for zero in (L for L in view["labels"] if L["text"] == "$r = 0$" and L["anchor"] in ("l", "r")):
                    for inner in (L for L in view["labels"] if L["text"] == "$r < r_-$" and L["anchor"] == "c"):
                        if zero["at"][0] * inner["at"][0] <= 0 or abs(zero["at"][0]) < abs(inner["at"][0]):
                            continue
                        seen += 1
                        lines = abs((zero["at"][1] - inner["at"][1]) * scale - (zero["dy"] - inner["dy"])) - 1.5 * size
                        along = (abs(zero["at"][0] * scale + zero["dx"])
                                 - abs(inner["at"][0] * scale + inner["dx"]) - 1.9 * size)
                        self.assertTrue(lines >= 0 or along >= size,
                                        f"{metric_id}/{view['id']}: $r = 0$ and $r < r_-$ are {along:.1f} apart on one line")
        found(self, seen, "tower with r = 0 outside a region r < r_-")

    def test_a_change_that_draws_nothing_leaves_the_diagrams_standing(self):
        changed = copy.deepcopy(self.metrics)
        for metric in changed:
            metric["history"] = metric.get("history", "") + " More."
            for system in metric.get("coordinates") or []:
                system["domains"] = (system.get("domains") or []) + ["x \\in \\mathbb{R}"]
                for parameter in system.get("parameters", []):
                    parameter["description"] = parameter.get("description", "") + " Reworded."
        self.assertEqual({name for name in build.load_conformal(changed) if chosen(name)}, set(self.conformal))

    def test_a_file_without_a_metric_or_a_source_is_refused(self):
        metric = {"id": "x", "name": "X", "short_name": "X", "tags": ["t"], "coordinates": [{"id": "a"}]}
        with self.assertRaises(build.DataError) as raised:
            self.load_folder({"ghost.json": {"metric": "ghost", "source": [], "views": []}}, [metric])
        self.assertIn("ghost", str(raised.exception))
        with self.assertRaises(build.DataError) as raised:
            self.load_folder({"x.json": {"metric": "x", "source": [], "views": []}}, [metric])
        self.assertIn("drawn from", str(raised.exception))
        source = [{"metric": "x", "system": "b", "fields": ["coords"], "version": "0"}]
        with self.assertRaises(build.DataError) as raised:
            self.load_folder({"x.json": {"metric": "x", "source": source, "views": []}}, [metric])
        self.assertIn("'b'", str(raised.exception))

    def test_a_file_that_draws_nothing_is_refused(self):
        """A spacetime with no conformal diagram has no file, rather than one the page would find empty."""
        kerr = self.named["kerr"]
        for data in (dict(kerr, views=[]), {key: value for key, value in kerr.items() if key != "views"}):
            with self.assertRaises(build.DataError) as raised:
                self.load_folder({"kerr.json": data}, self.metrics)
            self.assertIn("draws nothing", str(raised.exception))

    def test_every_caption_names_what_is_drawn(self):
        for metric_id, data in self.conformal.items():
            for view in data["views"]:
                self.assertTrue(view["caption"], f"{metric_id} {view['id']}")
                self.assertTrue(opens_with_a_noun_phrase(view["caption"][0]), f"{metric_id} {view['id']}")

    @checks(*SLICES)
    def test_the_spacetimes_drawn_on_a_slice_say_so_on_the_diagram(self):
        """A view of a surface that is not the whole spacetime carries its restriction, and the
        three spacetimes with no faithful picture of the whole carry it on every view."""
        slices = set(filter(chosen, self.SLICES))
        self.assertLessEqual(slices, set(self.conformal))
        for metric_id in slices:
            for view in self.conformal[metric_id]["views"]:
                self.assertTrue(view.get("restriction"), f"{metric_id} {view['id']}")
                self.assertRegex(view["restriction"], r"only", f"{metric_id} {view['id']}")

    def test_every_text_is_tex_with_its_mathematics_closed(self):
        for name, data in self.conformal.items():
            for field, value in conformal_prose(name, data):
                self.assertTrue(value.strip(), field)
                self.assertEqual(value.replace("\\$", "").count("$") % 2, 0, field)
                self.assertEqual(value.count("{"), value.count("}"), field)

    def test_every_view_draws_inside_its_box_and_names_only_what_it_draws(self):
        kinds = {"fill", "line", "zig", "point"}
        for name, data in self.conformal.items():
            ids = [view["id"] for view in data["views"]]
            self.assertEqual(len(ids), len(set(ids)), name)
            for view in data["views"]:
                where = f"{name} {view['id']}"
                x0, x1, t0, t1 = view["box"]
                self.assertTrue(x0 < x1 and t0 < t1, where)
                drawn = {}
                for layer in view["layers"]:
                    self.assertIn(layer["kind"], kinds, where)
                    drawn.setdefault(layer["class"], layer["kind"])
                    points = [layer["at"]] if layer["kind"] == "point" else layer["points"]
                    for x, t in points:
                        self.assertTrue(x0 <= x <= x1 and t0 <= t <= t1, f"{where} {layer['class']} at {x}, {t}")
                for label in view["labels"]:
                    x, t = label["at"]
                    self.assertTrue(x0 <= x <= x1 and t0 <= t <= t1, f"{where} label {label['text']}")
                self.assertIn("region", drawn, where)
                for kind, cls, _ in view["legend"]:
                    self.assertEqual(drawn.get(cls), kind, f"{where} legend {cls}")

    def load_folder(self, files, metrics):
        with tempfile.TemporaryDirectory() as folder:
            for filename, data in files.items():
                (Path(folder) / filename).write_text(json.dumps(data), encoding="utf-8")
            with mock.patch.object(build, "CONFORMAL_DIR", Path(folder)):
                return build.load_conformal(metrics)


class EmbeddingDiagrams(unittest.TestCase):
    """An embedding diagram is a slice of a spacetime drawn as a surface of revolution, from
    what its metrics publish, and stops being published when any of that changes. The file is
    the definition in _tools/README.md, which the application builds against."""

    ENDS = {"axis", "apex", "join", "crease", "throat", "edge", "stops"}

    def setUp(self):
        self.metrics = build.load_metrics()
        self.embedding = embedding_files()
        self.named = NamedFiles(build.EMBEDDING_DIR)
        found(self, self.embedding, "embedding diagram file")

    @checks("btz")
    def test_btz_is_one_surface_through_both_exteriors(self):
        """The BTZ moment t = 0 as the file holds it: each exterior a piece in flat space from the
        throat r_+ = l to sqrt(2) l and a piece in Minkowski space beyond, rho = r throughout, the
        two parts of each exterior meeting level in one circle marked `space`, and the two
        exteriors meeting at the throat."""
        views = self.named["btz"]["views"]
        self.assertEqual(len(views), 1)
        surface = views[0]["surfaces"][0]
        pieces = {p["id"]: p for p in surface["pieces"]}
        level = math.sqrt(2)
        for inner, outer in (("exterior", "exterior_minkowski"), ("other_exterior", "other_exterior_minkowski")):
            a, b = pieces[inner], pieces[outer]
            self.assertNotIn("space", a, inner)
            self.assertEqual(b.get("space"), "minkowski", outer)
            self.assertEqual((a["points"][0][0], a["points"][-1][0], b["points"][0][0]), (1.0, level, level))
            for x, rho, z in a["points"] + b["points"]:
                self.assertLess(abs(rho - x), 1e-6, f"{inner} rho at {x}")
            for i in (1, 2):
                self.assertLess(abs(a["points"][-1][i] - b["points"][0][i]), 1e-7, f"{inner} meets {outer}")
            # Level on both sides of the circle: the chords next to it climb far less than they run.
            for P, Q in ((a["points"][-2], a["points"][-1]), (b["points"][0], b["points"][1])):
                self.assertLess(abs(Q[2] - P[2]), 0.1 * abs(Q[1] - P[1]), f"{inner} and {outer} level at the join")
            self.assertEqual(a["end"]["kind"], "join")
            self.assertEqual(b["start"]["kind"], "join")
            joins = [r for r in surface["rings"] if r["class"] == "space" and r["piece"] == inner]
            self.assertEqual([r["x"] for r in joins], [level], inner)
        self.assertEqual(pieces["exterior"]["points"][0][1:], pieces["other_exterior"]["points"][0][1:])
        self.assertEqual(pieces["exterior"]["start"]["kind"], "throat")

    @checks("spinning_string")
    def test_the_spinning_string_leaves_the_axis_along_the_light_cone_and_ends_on_its_cone(self):
        """The spinning string's moment t = 0 as the file holds it, at b = 0.9 and a = 0.9: a piece
        in Minkowski space from the null circle, where it runs along the light cone, to the level
        circle R = a b/sqrt(1 - b^2), and a piece in flat space beyond, whose last chord climbs at
        sqrt(g_RR - 1), a little below the slope sqrt(1/b^2 - 1) of the cosmic string's cone; every chord of both is the
        proper distance sqrt(g_RR) dR of the published metric, g_RR = R^2/(b^2 (R^2 + a^2))."""
        a = b = 0.9
        level = a * b / math.sqrt(1 - b * b)
        surface, = self.named["spinning_string"]["views"][0]["surfaces"]
        near, far = surface["pieces"]
        self.assertEqual(near.get("space"), "minkowski")
        self.assertNotIn("space", far)
        self.assertLess(abs(near["points"][-1][0] - level), 1e-12)
        self.assertEqual(near["points"][-1][0], far["points"][0][0])
        self.assertLess(near["points"][0][1], 0.006)
        P, Q = near["points"][:2]
        self.assertLess(abs((Q[2] - P[2]) / (Q[1] - P[1]) - 1), 1e-3, "along the light cone at the null circle")
        P, Q = far["points"][-2:]
        cone = math.sqrt(1 / b ** 2 - 1)
        slope, mid = (Q[2] - P[2]) / (Q[1] - P[1]), (P[0] + Q[0]) / 2
        self.assertLess(abs(slope - math.sqrt(mid ** 2 / (b * b * (mid ** 2 + a * a)) - 1)), 1e-3)
        self.assertTrue(0.85 * cone < slope < cone, "nearing the cone's slope from below at 5 r_c")
        proper = lambda R: math.sqrt(R * R + a * a) / b      # the proper radius, whose differences are distances
        for piece, sign in ((near, -1), (far, 1)):
            for P, Q in zip(piece["points"], piece["points"][1:]):
                chord = math.sqrt((Q[1] - P[1]) ** 2 + sign * (Q[2] - P[2]) ** 2)
                want = proper(Q[0]) - proper(P[0])
                self.assertLess(abs(chord - want), 3e-4 * want, f"{piece['id']} at R = {P[0]}")
            for R, rho, _ in piece["points"]:
                self.assertLess(abs(rho - R), 1e-6, f"{piece['id']} rho at {R}")

    def test_every_file_was_drawn_from_what_its_metrics_publish(self):
        by_id = {m["id"]: m for m in self.metrics}
        for metric_id, data in self.embedding.items():
            self.assertTrue(data["source"], metric_id)
            for source in data["source"]:
                system = next(s for s in by_id[source["metric"]]["coordinates"] if s["id"] == source["system"])
                self.assertEqual(build.diagram_source_version(system, source["fields"]), source["version"],
                                 f"{metric_id} from {source['metric']}/{source['system']}")
        self.assertEqual({name for name in build.load_embedding(self.metrics) if chosen(name)}, set(self.embedding))

    def test_a_changed_component_is_refused_and_leaves_the_files_alone(self):
        # The star is drawn with schwarzschild.json's exterior, so a change there stops it too.
        changed = copy.deepcopy(self.metrics)
        for metric in changed:
            if metric["id"] == "schwarzschild":
                system = next(s for s in metric["coordinates"] if s["id"] == "spherical")
                system["metric_components"][1]["value"] = "2" + system["metric_components"][1]["value"]
        before = {path: path.read_text(encoding="utf-8") for path in (build.INDEX_FILE, build.REFERENCES_FILE)}
        for argv in (["--check"], []):
            stderr = io.StringIO()
            with mock.patch.object(build, "load_metrics", return_value=changed), \
                    mock.patch.object(build, "load_diagrams", return_value={}), \
                    mock.patch.object(build, "load_conformal", return_value={}), \
                    contextlib.redirect_stderr(stderr):
                self.assertEqual(build.main(argv), 2)
            self.assertIn("embedding/", stderr.getvalue())
            self.assertIn("embedding.py --metric", stderr.getvalue())
        for path, text in before.items():
            self.assertEqual(path.read_text(encoding="utf-8"), text)
        with self.assertRaises(build.DataError) as raised:
            self.load_folder({"interior_schwarzschild.json": self.named["interior_schwarzschild"]}, changed)
        self.assertIn("schwarzschild.json", str(raised.exception))

    def test_a_change_that_draws_nothing_leaves_the_diagrams_standing(self):
        changed = copy.deepcopy(self.metrics)
        for metric in changed:
            metric["history"] = metric.get("history", "") + " More."
            for system in metric.get("coordinates") or []:
                system["domains"] = (system.get("domains") or []) + ["x \\in \\mathbb{R}"]
                for parameter in system.get("parameters", []):
                    parameter["description"] = parameter.get("description", "") + " Reworded."
        self.assertEqual({name for name in build.load_embedding(changed) if chosen(name)}, set(self.embedding))

    def test_a_file_without_a_metric_or_a_source_or_a_surface_is_refused(self):
        metric = {"id": "x", "name": "X", "short_name": "X", "tags": ["t"], "coordinates": [{"id": "a"}]}
        with self.assertRaises(build.DataError) as raised:
            self.load_folder({"ghost.json": {"metric": "ghost", "source": [], "views": []}}, [metric])
        self.assertIn("ghost", str(raised.exception))
        with self.assertRaises(build.DataError) as raised:
            self.load_folder({"x.json": {"metric": "x", "source": [], "views": []}}, [metric])
        self.assertIn("drawn from", str(raised.exception))
        flamm = self.named["schwarzschild"]
        for data in (dict(flamm, views=[]), dict(flamm, views=[dict(flamm["views"][0], surfaces=[])])):
            with self.assertRaises(build.DataError) as raised:
                self.load_folder({"schwarzschild.json": data}, self.metrics)
            self.assertIn("draws nothing", str(raised.exception))

    def test_every_spacetime_draws_and_a_file_that_draws_nothing_says_why(self):
        # Every spacetime but those WITHOUT_EMBEDDING names has a surface drawn, and none says beside
        # its views that it draws nothing.
        self.assertEqual(set(self.embedding),
                         {m["id"] for m in self.metrics if chosen(m["id"])} - build.WITHOUT_EMBEDDING)
        for name, data in self.embedding.items():
            self.assertTrue(data["views"], name)
            self.assertNotIn("stops", data, f"{name} keeps what it does not draw in its views")
        # A file with no surface to draw says why in sentences, and a file with views may not.
        flamm = self.named["schwarzschild"]
        stated = {"metric": "schwarzschild", "source": flamm["source"], "views": [], "stops": ["Nothing to draw."]}
        self.assertEqual(set(self.load_folder({"schwarzschild.json": stated}, self.metrics)), {"schwarzschild"})
        with self.assertRaises(build.DataError) as raised:
            self.load_folder({"schwarzschild.json": dict(flamm, stops=["Nothing."])}, self.metrics)
        self.assertIn("says beside them", str(raised.exception))
        for data in (dict(stated, stops=[]), dict(stated, stops=[" "])):
            with self.assertRaises(build.DataError):
                self.load_folder({"schwarzschild.json": data}, self.metrics)

    def test_every_curve_and_point_marked_on_a_surface_lies_on_its_piece(self):
        # A curve or a point stands at a distance from the axis that the piece's profile reaches,
        # at the height the profile has there, between the two points of it that it falls between.
        marked = 0
        for name, data in self.embedding.items():
            for view in data["views"]:
                for number, surface in enumerate(view["surfaces"]):
                    pieces = {p["id"]: p for p in surface["pieces"]}
                    for mark in surface.get("curves", []) + surface.get("dots", []):
                        where = f"{name} {view['id']} surface {number} {mark['class']}"
                        points = mark["points"] if "points" in mark else [mark["at"]]
                        self.assertGreaterEqual(len(points), 1 if "at" in mark else 2, where)
                        self.assertIn(mark.get("closed", False), (False, True), where)
                        if pieces[mark["piece"]].get("grid", {}).get("frame") == "ellipses":
                            # A mark on a stack of ellipses is a node of the row at its height.
                            nodes = grid_nodes(pieces[mark["piece"]])
                            for X, Y, Z in points:
                                row = next(r for r in nodes if r[0][2] == Z)
                                self.assertLess(min(math.dist((X, Y, Z), q) for q in row), 1e-6, f"{where} at {X}, {Y}")
                            marked += 1
                            continue
                        if "grid" in pieces[mark["piece"]]:
                            for X, Y, Z in points:
                                z = grid_height(pieces[mark["piece"]], X, Y)
                                self.assertIsNotNone(z, f"{where} at {X}, {Y}")
                                self.assertLess(abs(Z - z), 1e-6, f"{where} at {X}, {Y}")
                            marked += 1
                            continue
                        profile = sorted((rho, z) for _, rho, z in pieces[mark["piece"]]["points"])
                        for X, Y, Z in points:
                            r = math.hypot(X, Y)
                            self.assertTrue(profile[0][0] - 1e-6 <= r <= profile[-1][0] + 1e-6, f"{where} at {r}")
                            k = next((i for i in range(1, len(profile)) if profile[i][0] >= r), len(profile) - 1)
                            (r0, z0), (r1, z1) = profile[k - 1], profile[k]
                            z = z0 if r1 == r0 else z0 + (z1 - z0) * (r - r0) / (r1 - r0)
                            self.assertLess(abs(Z - z), 1e-4, f"{where} at {r}")
                        marked += 1
        found(self, marked, "curve or point marked on a surface")

    def test_a_piece_that_doubles_back_or_crosses_the_axis_is_refused(self):
        for spoil, words in ((lambda points: points.reverse() or points.insert(1, points[0]), "one way"),
                             (lambda points: points[3].__setitem__(1, -0.5), "below the axis")):
            data = copy.deepcopy(self.named["schwarzschild"])
            spoil(data["views"][0]["surfaces"][0]["pieces"][0]["points"])
            with self.assertRaises(build.DataError) as raised:
                self.load_folder({"schwarzschild.json": data}, self.metrics)
            self.assertIn(words, str(raised.exception))

    def test_no_file_writes_a_negative_zero(self):
        # Rounding a small negative number leaves -0.0, which a client would print as "-0".
        # The match alone is reported, since the file is too long to print.
        for path in sorted(p for p in build.EMBEDDING_DIR.glob("*.json") if chosen(p.stem)):
            text = path.read_text(encoding="utf-8")
            found = re.search(r"(?<![\d.])-0(\.0*)?(?![\d.eE])", text)
            self.assertIsNone(found and text[max(found.start() - 40, 0):found.end() + 10], path.name)

    def test_every_surface_has_the_shape_its_definition_promises(self):
        for name, data in self.embedding.items():
            for view in data["views"]:
                where = f"{name} {view['id']}"
                self.assertTrue(view["unit"].startswith("$") and view["unit"].endswith("$"), where)
                sequence = len(view["surfaces"]) > 1
                times = []
                for surface in view["surfaces"]:
                    self.assertEqual("label" in surface, sequence, where)
                    self.assertEqual("time" in surface, sequence, where)
                    times.append(surface.get("time", 0))
                    ids = [piece["id"] for piece in surface["pieces"]]
                    self.assertEqual(len(ids), len(set(ids)), where)
                    for piece in surface["pieces"]:
                        if "grid" in piece:
                            grid = piece["grid"]
                            self.assertIn(grid["frame"], ("polar", "cartesian", "ellipses"), where)
                            self.assertEqual(len(grid["z"]), len(grid["u"]), where)
                            if grid["frame"] == "ellipses":
                                self.assertEqual((len(grid["a"]), len(grid["b"])), (len(grid["u"]), len(grid["u"])), where)
                            else:
                                self.assertTrue(all(len(row) == len(grid["v"]) for row in grid["z"]), where)
                            self.assertEqual(piece["edge"]["kind"], "edge", where)
                            self.assertIn("height", view, f"{where}: a grid says what its height is")
                            continue
                        self.assertIn(piece["start"]["kind"], self.ENDS, where)
                        self.assertIn(piece["end"]["kind"], self.ENDS, where)
                        for point in piece["points"]:
                            self.assertEqual(len(point), 3, where)
                        if piece["start"]["kind"] in ("axis", "apex"):
                            self.assertEqual(piece["points"][0][1], 0, f"{where} {piece['id']} starts on the axis")
                    for ring in surface["rings"]:
                        piece = next(p for p in surface["pieces"] if p["id"] == ring["piece"])
                        self.assertIn([ring["x"], ring["rho"], ring["z"]], piece["points"], f"{where} ring {ring}")
                self.assertEqual(times, sorted(times), where)

    def test_a_ring_is_classed_by_the_piece_it_lies_on(self):
        # The README gives r to a circle on a sheet or star piece and r2 to one on a sheet2
        # piece, which is how a client tells the two sides of a throat apart.
        allowed = {"sheet": {"r"}, "star": {"r"}, "sheet2": {"r2"}, "reference": {"reference"}}
        for name, data in self.embedding.items():
            for view in data["views"]:
                for surface in view["surfaces"]:
                    kinds = {p["id"]: p["class"] for p in surface["pieces"]}
                    for ring in surface["rings"]:
                        if ring["class"] in ("r", "r2", "reference"):
                            self.assertIn(ring["class"], allowed[kinds[ring["piece"]]],
                                          f"{name} {view['id']} ring {ring['class']} on {ring['piece']}")

    @checks("teo_wormhole")
    def test_teos_equator_is_flamms_paraboloid_and_his_throat_a_dumbbell(self):
        """Teo's eq. (27), z = +-2 sqrt(b_0 (r - b_0)) on the equator for every spin, and the
        throat's circles of radius b_0 (1 + cos^2 theta) sin theta at a = 1/4, the same on either
        side of its equator, from the numbers written and nothing else."""
        views = {view["id"]: view["surfaces"][0] for view in self.named["teo_wormhole"]["views"]}
        sides = {p["id"]: p["points"] for p in views["equator"]["pieces"]}
        for sign, pid in ((1, "near"), (-1, "far")):
            for r, rho, z in sides[pid]:
                self.assertLess(abs(rho - r), 2e-6, f"Teo's equator, rho at {r}")
                self.assertLess(abs(z - sign * 2 * math.sqrt(max(r - 1, 0))), 2e-6, f"Teo's equator, z at {r}")
        ergo = [ring for ring in views["equator"]["rings"] if ring["class"] == "ergo"]
        self.assertEqual([round(ring["rho"] ** 2, 6) for ring in ergo], [2.0, 2.0])
        throat = views["throat"]["pieces"][0]["points"]
        for theta, rho, z in throat:
            self.assertLess(abs(rho - (1 + math.cos(theta) ** 2) * math.sin(theta)), 2e-6, f"Teo's throat at {theta}")
        self.assertLess(abs(throat[0][2] + throat[-1][2]), 2e-6)
        self.assertGreater(max(rho for _, rho, _ in throat), 1.08)

    @checks("bardeen", "bertotti_robinson", "black_string", "cosmic_string", "de_sitter", "einstein_static",
            "ellis_bronnikov", "frw", "global_monopole", "godel", "gravastar", "hayward", "interior_schwarzschild",
            "israel_shell", "kaluza_klein_monopole", "kantowski_sachs", "kerr", "kerr_newman", "mass_inflation",
            "morris_thorne", "myers_perry", "nariai", "oppenheimer_snyder", "photon_rocket",
            "point_particle_2plus1", "rn_metric", "schwarzschild", "som_raychaudhuri", "stockum_dust",
            "string_black_hole", "tangherlini", "taub_nut", "tolman_bondi", "tov", "vaidya", "van_den_broeck")
    def test_the_surfaces_are_the_ones_known_in_closed_form(self):
        """Flamm's paraboloid, the interior Schwarzschild cap, the catenoid, the cone, Gott's cap
        and the sphere, from the numbers written and nothing else, to their rounding."""
        def piece(metric_id, piece_id, number=0, view=0):
            surface = self.named[metric_id]["views"][view]["surfaces"][number]
            points = next(p for p in surface["pieces"] if p["id"] == piece_id)["points"]
            return [tuple(point) for point in points]

        def near(a, b, where):
            self.assertLess(abs(a - b), 2e-6, where)
        for sign, pid in ((1, "exterior"), (-1, "other_exterior")):
            for r, rho, z in piece("schwarzschild", pid):
                near(rho, r, f"Flamm rho at {r}")
                near(z, sign * 2 * math.sqrt(r - 1), f"Flamm z at {r}")
            # Tangherlini's slice in five dimensions is the catenoid z = r_h arcosh(r/r_h).
            for r, rho, z in piece("tangherlini", pid):
                near(rho, r, f"Tangherlini's catenoid rho at {r}")
                near(z, sign * math.acosh(r), f"Tangherlini's catenoid z at {r}")
            # The black string across the string: Flamm's paraboloid in five dimensions, the catenoid in six.
            for r, rho, z in piece("black_string", pid):
                near(rho, r, f"the black string's paraboloid rho at {r}")
                near(z, sign * 2 * math.sqrt(r - 1), f"the black string's paraboloid z at {r}")
            for r, rho, z in piece("black_string", pid, view=2):
                near(rho, r, f"the black string's catenoid rho at {r}")
                near(z, sign * math.acosh(r), f"the black string's catenoid z at {r}")
        # The black string's rippling horizon on a circle of length 10 r_s: at each moment the circle at z
        # has the radius r_s (1 + a cos(2 pi z/10)), the ripple starts at a = 0.02 and grows as
        # exp(0.0633 v), Gregory and Laflamme's rate at that wavelength, and every chord of the profile
        # is as long as the step along z it spans, since the mode has no component along z.
        ripple = self.named["black_string"]["views"][1]
        self.assertEqual(ripple["id"], "ripple")
        for surface in ripple["surfaces"]:
            points = [tuple(point) for point in surface["pieces"][0]["points"]]
            a = 0.02 * math.exp(0.0633 * surface["time"])
            # The rate is quoted to four decimals, half a unit of the last over 42 r_s of v.
            self.assertLess(abs(max(rho for _, rho, _ in points) - 1 - a), 3e-3 * a, surface["label"])
            a = max(rho for _, rho, _ in points) - 1
            for (x, rho, z), (x2, rho2, z2) in zip(points, points[1:]):
                near(rho, 1 + a * math.cos(2 * math.pi * x / 10), f"the rippled horizon's rho at {x}, {surface['label']}")
                self.assertLess(abs(math.hypot(rho2 - rho, z2 - z) - (x2 - x)), 2e-4 * (x2 - x), f"the chord at {x}")
            self.assertEqual((points[0][0], points[-1][0]), (-5.0, 5.0))
            # Myers and Perry's hole with one spin in five dimensions, mu = 1 and a = 3/5: in the plane of
            # rotation the circle at r has the radius sqrt(r^2 + a^2 + mu a^2/r^2), which is mu/r_+ = 5/4 at
            # the throat, and in the transverse plane the slice is z = sqrt(mu) arcosh(r/r_+) with r_+ = 4/5.
            # With equal spins, a = 2/5, the Hopf fibre at rho has the radius sqrt(rho^4 + mu a^2)/rho.
            rotation = piece("myers_perry", pid, view=0)
            self.assertEqual(tuple(rotation[0]), (0.8, 1.25, 0.0))
            for r, rho, z in rotation:
                near(rho, math.sqrt(r * r + 0.36 + 0.36 / (r * r)), f"Myers-Perry's plane of rotation rho at {r}")
            for r, rho, z in piece("myers_perry", pid, view=1):
                near(rho, r, f"Myers-Perry's transverse plane rho at {r}")
                near(z, sign * math.acosh(max(r / 0.8, 1.0)), f"Myers-Perry's transverse plane z at {r}")
            fibre = piece("myers_perry", pid, view=2)
            self.assertAlmostEqual(fibre[0][1], 1.0, places=6)
            for r, rho, z in fibre:
                near(rho, math.sqrt(r ** 4 + 0.16) / r, f"Myers-Perry's Hopf fibre rho at {r}")
            for r, rho, z in piece("myers_perry", pid, view=3):
                near(rho, r, f"Myers-Perry's transverse plane in six dimensions rho at {r}")
        # The Kaluza-Klein monopole's cigar: the circle of the fifth dimension at r has radius
        # 8m sqrt(r/(r + 4m)), and the surface rises from the nut at dz/dr = 2 toward dz/dr = 1.
        cigar = piece("kaluza_klein_monopole", "cigar")
        self.assertEqual(list(cigar[0]), [0.0, 0.0, 0.0])
        for r, rho, z in cigar:
            near(rho, 8 * math.sqrt(r / (r + 4)), f"the cigar's rho at {r}")
        for (r0, _, z0), (r1, _, z1) in zip(cigar, cigar[1:]):
            slope = (z1 - z0) / (r1 - r0)
            mid = 0.5 * (r0 + r1)
            want = math.sqrt((mid ** 3 + 16 * mid ** 2 + 96 * mid + 256) / (mid + 4) ** 3)
            self.assertLess(abs(slope - want), 2e-3, f"the cigar's dz/dr at {mid}")
        cap = piece("interior_schwarzschild", "star")
        for r, rho, z in cap:
            near(rho, r, f"cap rho at {r}")
            near(z - cap[0][2], math.sqrt(27 / 8) - math.sqrt(27 / 8 - r * r), f"cap z at {r}")
        for r, rho, z in piece("interior_schwarzschild", "exterior"):
            near(z, 2 * math.sqrt(r - 1), f"star exterior z at {r}")
        # The gravastar at R = 5/4 and L = 2: a cap of the sphere of radius L inside the shell, Flamm's
        # paraboloid outside it, and one circle where they meet, which the cap reaches at dz/dr = 0.80.
        core, outside = piece("gravastar", "interior"), piece("gravastar", "exterior")
        for r, rho, z in core:
            near(rho, r, f"gravastar cap rho at {r}")
            near(z - core[0][2], 2 - math.sqrt(4 - r * r), f"gravastar cap z at {r}")
        for r, rho, z in outside:
            near(z, 2 * math.sqrt(r - 1), f"gravastar exterior z at {r}")
        self.assertEqual(core[-1][:1] + core[-1][1:], outside[0][:1] + outside[0][1:], "the gravastar's shell is one circle")
        self.assertEqual(core[-1][0], 1.25)
        # The neutron star's exterior is Flamm's paraboloid of its mass, whose throat, at 2M, is
        # where the vacuum drawn under the star begins.
        two_m = piece("tov", "vacuum")[0][0]
        for r, rho, z in piece("tov", "exterior"):
            near(rho, r, f"neutron star rho at {r}")
            near(z, 2 * math.sqrt(two_m * (r - two_m)), f"neutron star exterior z at {r}")
        for sign, pid in ((1, "near"), (-1, "far")):
            for r, rho, z in piece("morris_thorne", pid):
                near(z, sign * math.acosh(r), f"catenoid z at {r}")
        # Ellis-Bronnikov's proper r runs through the throat, one piece of the same catenoid.
        for r, rho, z in piece("ellis_bronnikov", "whole"):
            near(rho, math.sqrt(r * r + 1), f"Ellis-Bronnikov rho at {r}")
            near(z, math.asinh(r), f"Ellis-Bronnikov z at {r}")
        # Van Den Broeck's pocket in the proper distance l, in units of R: a floor of radius 3/2, a rim
        # that widens to 7/4, a lid in to 1/2, a neck that narrows to 1/4, and the plane outside, pi/2
        # above the floor; the rim and the neck each climb pi/4 as (w sqrt(1 - w^2) + arcsin w)/4 + pi/8.
        def arc(w):
            return (w * math.sqrt(1 - w * w) + math.asin(w)) / 4 + math.pi / 8
        pocket = piece("van_den_broeck", "pocket")
        self.assertEqual([pocket[0][0], pocket[-1][0]], [-4.0, 1.0])
        for l, rho, z in pocket:
            want = ((l + 4, 0.0) if l < -2.5 else (1.75 - (l + 2) ** 2, arc(2 * (l + 2))) if l < -1.5
                    else (-l, math.pi / 4) if l < -0.5 else (l * l + 0.25, math.pi / 4 + arc(2 * l)) if l < 0.5
                    else (l, math.pi / 2))
            near(rho, want[0], f"the pocket's rho at {l}")
            near(z, want[1], f"the pocket's z at {l}")
        self.assertAlmostEqual(min(rho for l, rho, z in pocket if l > -1), 0.25, places=6)
        self.assertAlmostEqual(max(rho for l, rho, z in pocket), 1.75, places=6)
        # Reissner-Nordstrom's circles have their areal radius on both views, and inside r- each
        # side starts level at r_q^2/r_s, 0.2304 r_s, and ends at r- = 0.36 r_s.
        for view, pids in ((0, ("exterior", "other_exterior")), (1, ("inside", "other_inside"))):
            for pid in pids:
                points = piece("rn_metric", pid, view=view)
                for r, rho, z in points:
                    near(rho, r, f"Reissner-Nordstrom rho at {r}")
                if view:
                    self.assertEqual((points[0][0], points[-1][0]), (0.2304, 0.36), pid)
        # Bardeen's circles have their areal radius on both views, at g = r_s/3. Outside, the throat is
        # the outer zero of f, and inside, each side runs from the centre to the inner zero. Each chord
        # climbs as dz/dr = sqrt(1/f - 1), and at the centre the cap has the curvature of the sphere
        # of radius sqrt(g^3/r_s), de Sitter's.
        for view, pids in ((0, ("exterior", "other_exterior")), (1, ("inside", "other_inside"))):
            for pid in pids:
                points = piece("bardeen", pid, view=view)
                for r, rho, z in points:
                    near(rho, r, f"Bardeen rho at {r}")
                horizon = points[-1][0] if view else points[0][0]
                self.assertLess(abs(bardeen_f(horizon)), 1e-12, f"Bardeen {pid} ends on a horizon")
                self.assertEqual(points[0][0], 0.0 if view else horizon, pid)
                for (r0, _, z0), (r1, _, z1) in zip(points, points[1:]):
                    if min(abs(r0 - horizon), abs(r1 - horizon)) > 0.02:
                        # Simpson's rule for the rise along the chord, which lies clear of the horizon.
                        slopes = [math.sqrt(max(1 / bardeen_f(r) - 1, 0.0)) for r in (r0, 0.5 * (r0 + r1), r1)]
                        rise = (r1 - r0) * (slopes[0] + 4 * slopes[1] + slopes[2]) / 6
                        self.assertLess(abs(abs(z1 - z0) - rise), 1e-3 * rise + 2e-7, f"Bardeen {pid} from {r0} to {r1}")
        cap = piece("bardeen", "inside", view=1)
        a = (1 / 3) ** 1.5
        for r, _, z in cap[1:]:
            if r < 0.02:
                self.assertLess(abs(2 * a * (z - cap[0][2]) / r ** 2 - 1), 2e-2, f"Bardeen's cap at {r}")
        # Hayward's static slice at ell = 12m/(7 sqrt 7): outside, from the throat r_+ = 12m/7, and inside,
        # a closed surface from the centre to r_- = 6m/7 and on to the other centre, each climbing at
        # dz/dr = sqrt(1/F - 1) and leaving the centre as the sphere of radius ell does.
        ell2 = 144 / 343

        def hayward_slope(r, mass=1.0):
            return math.sqrt(2 * mass * r * r / (r ** 3 - 2 * mass * r * r + 2 * mass * ell2))
        for view, pids, ends in ((0, ("exterior", "other_exterior"), (12 / 7, 6.0)),
                                 (1, ("inside", "other_inside"), (0.0, 6 / 7))):
            for sign, pid in zip((1, -1), pids):
                points = piece("hayward", pid, view=view)
                self.assertAlmostEqual(points[0][0], ends[0], places=12)
                self.assertAlmostEqual(points[-1][0], ends[1], places=12)
                for r, rho, z in points:
                    near(rho, r, f"Hayward rho at {r}")
                for (r0, _, z0), (r1, _, z1) in zip(points, points[1:]):
                    mid = 0.5 * (r0 + r1)
                    if min(abs(mid - 6 / 7), abs(mid - 12 / 7)) > 0.05 and r1 - r0 > 1e-4:
                        self.assertLess(abs(sign * (z1 - z0) / (r1 - r0) - hayward_slope(mid)),
                                        2e-2 * (1 + hayward_slope(mid)), f"Hayward's dz/dr at {mid}")
        cap = piece("hayward", "inside", view=1)
        for r, _, z in cap:
            if r <= 0.03:
                near(z - cap[0][2], r * r / (2 * math.sqrt(ell2)), f"Hayward's centre, the sphere of radius ell, at {r}")
        # The hole that forms and evaporates: flat before the radiation reaches the rim and after the
        # last of it has left, and between them the trapping horizons marked where r^3 - 2m r^2 +
        # 2m ell^2 vanishes at the mass of the advanced time v = T + r, two of them or none.
        def hayward_mass(v):
            a, b = min(max(v, 0.0), 2.0), min(max(v, 4.0), 8.0)
            return math.sin(math.pi * a / 4) ** 2 * math.cos(math.pi * (b - 4) / 8) ** 2
        frames = self.named["hayward"]["views"][2]["movie"]["frames"]
        self.assertEqual((frames[0]["value"], frames[-1]["value"]), (-6.0, 8.0))
        trapped = 0
        for frame in frames:
            T, points = frame["value"], frame["pieces"][0]["points"]
            self.assertEqual(points[-1][2], 0)
            if T <= -6 or T >= 8:
                self.assertTrue(all(z == 0 for _, _, z in points), f"Hayward at v - r = {T}")
            horizons = [ring["rho"] for ring in frame["rings"] if ring["class"] == "horizon"]
            self.assertIn(len(horizons), (0, 2), f"Hayward at v - r = {T}")
            trapped += bool(horizons)
            for r in horizons:
                mass = hayward_mass(T + r)
                self.assertGreater(mass, 3 * math.sqrt(3 * ell2) / 4 - 1e-9)
                self.assertLess(abs(r ** 3 - 2 * mass * r * r + 2 * mass * ell2), 1e-5, f"Hayward's horizon at v - r = {T}")
        self.assertGreater(trapped, 10)
        # The tail falling into the charged hole, r_q = 0.96 m_0: each slice of constant v - r climbs at
        # dz/dr = sqrt(2m/r - r_q^2/r^2) with the mass of each circle's own advanced time, from the
        # circle where that vanishes, and marks the two apparent horizons, the roots of
        # r^2 - 2m r + r_q^2, which move from 0.78 and 1.18 m_0 to 0.72 and 1.28 m_0.
        def tail_mass(v):
            return 1 - (10 / max(v, 10.0)) ** 11 / 50
        frames = self.named["mass_inflation"]["views"][0]["movie"]["frames"]
        self.assertEqual((frames[0]["value"], frames[-1]["value"]), (7.5, 16.5))
        for frame in frames:
            T, points = frame["value"], frame["pieces"][0]["points"]
            self.assertEqual(points[-1][2], 0)
            start = points[0][0]
            self.assertLess(abs(2 * tail_mass(T + start) * start - 0.9216), 3e-4, f"mass inflation at v - r = {T}")
            for (r0, _, z0), (r1, _, z1) in zip(points, points[1:]):
                mid = (r0 + r1) / 2
                if r1 - r0 > 1e-3:
                    slope = math.sqrt(2 * tail_mass(T + mid) / mid - 0.9216 / mid ** 2)
                    self.assertLess(abs((z1 - z0) / (r1 - r0) - slope), 2e-2 * (1 + slope), f"mass inflation dz/dr at {mid}")
            horizons = sorted(ring["x"] for ring in frame["rings"] if ring["class"] == "horizon")
            self.assertEqual(len(horizons), 2, f"mass inflation at v - r = {T}")
            for r in horizons:
                self.assertLess(abs(r * r - 2 * tail_mass(T + r) * r + 0.9216), 1e-9)
        first = sorted(r["x"] for r in frames[0]["rings"] if r["class"] == "horizon")
        last = sorted(r["x"] for r in frames[-1]["rings"] if r["class"] == "horizon")
        self.assertEqual([round(x, 2) for x in first + last], [0.78, 1.18, 0.72, 1.28])
        # On Kerr's equator the throat's circumference radius is 2GM/c^2 whatever the spin, and
        # Kerr-Newman's charge pulls it in to 2GM/c^2 - r_Q^2/r+, at a = 0.6 and r_Q = 0.5.
        self.assertAlmostEqual(piece("kerr", "exterior")[0][1], 2.0, places=6)
        r_plus = 1 + math.sqrt(0.39)
        self.assertAlmostEqual(piece("kerr_newman", "exterior")[0][0], r_plus, places=12)
        self.assertAlmostEqual(piece("kerr_newman", "exterior")[0][1], 2 - 0.25 / r_plus, places=6)
        # de Sitter's static slice is the sphere of radius l, a hemisphere to each static patch.
        for sign, pid in ((-1, "near"), (1, "far")):
            for r, rho, z in piece("de_sitter", pid):
                near(rho, r, f"de Sitter rho at {r}")
                near(z, sign * math.sqrt(max(1 - r * r, 0)), f"de Sitter z at {r}")
        # The Einstein static universe's moment is the round sphere of radius R, pole to antipode.
        for chi, rho, z in piece("einstein_static", "sphere"):
            near(rho, math.sin(chi), f"Einstein static rho at {chi}")
            near(z, -math.cos(chi), f"Einstein static z at {chi}")
        # Vaidya's slices of constant v - r: a flat disc inside the shell, and outside it Flamm's
        # paraboloid moved in by r_s, z = 2 sqrt(r), from the shell at z = 0 or from the axis.
        for number, surface in enumerate(self.named["vaidya"]["views"][0]["surfaces"]):
            ids = [p["id"] for p in surface["pieces"]]
            # The rim of the drawing, r = 4 r_s, stands at z = 0 at every moment.
            if "inside" in ids:
                inside = piece("vaidya", "inside", number)
                self.assertTrue(all(z == inside[0][2] for _, _, z in inside))
                outside = piece("vaidya", "outside", number)
                shell = outside[0][0]
                for r, rho, z in outside:
                    near(z - outside[0][2], 2 * (math.sqrt(r) - math.sqrt(shell)), f"Vaidya z at {r}")
                self.assertEqual(outside[-1][2], 0)
            else:
                whole = piece("vaidya", "whole", number)
                for r, rho, z in whole:
                    near(z, 2 * math.sqrt(r) - 4, f"Vaidya z at {r}")
        # Israel's shell of dust on the same slices: a flat disc inside the shell, whose radius is the
        # shell's on that slice, and z = 2 sqrt(r) outside it; the shell is one circle from both sides.
        for number, surface in enumerate(self.named["israel_shell"]["views"][0]["surfaces"]):
            ids = [p["id"] for p in surface["pieces"]]
            if "inside" in ids:
                inside = piece("israel_shell", "inside", number)
                self.assertTrue(all(z == inside[0][2] for _, _, z in inside))
                outside = piece("israel_shell", "outside", number)
                shell = outside[0][0]
                self.assertEqual(inside[-1], outside[0], "the shell is one circle")
                self.assertAlmostEqual(shell, israel_radius(israel_on_slice(surface["time"])), delta=1e-9)
                for r, rho, z in outside:
                    near(z - outside[0][2], 2 * (math.sqrt(r) - math.sqrt(shell)), f"Israel's shell z at {r}")
                self.assertEqual(outside[-1][2], 0)
            else:
                self.assertGreater(surface["time"], ISRAEL_END)
                for r, rho, z in piece("israel_shell", "whole", number):
                    near(z, 2 * math.sqrt(r) - 4, f"Israel's shell z at {r}")
        # The photon rocket's slices of constant cu + r at theta = pi/2: before the light of the burn
        # reaches them, z^2 = 8 m_0 r, and once it has passed the rim, z^2 = 8 m r with m = m_0 e^(-6/5),
        # the mass the burn leaves; the rim r = 12 m_0 stands at z = 0.
        rocket = self.named["photon_rocket"]["views"][0]["surfaces"]
        for number, mass in ((0, 1.0), (len(rocket) - 1, math.exp(-1.2))):
            for r, rho, z in piece("photon_rocket", "whole", number):
                near(rho, r, f"photon rocket rho at {r}")
                near(z, 2 * math.sqrt(2 * mass) * (math.sqrt(r) - math.sqrt(12)), f"photon rocket z at {r}")
        # Oppenheimer-Snyder's dust is a cap of a sphere of radius a out to chi0 = pi/4 at every
        # moment, and at the release, the first moment, the outside is Flamm's paraboloid.
        for number in range(len(self.named["oppenheimer_snyder"]["views"][0]["surfaces"])):
            dust = piece("oppenheimer_snyder", "dust", number)
            a = dust[-1][1] / math.sin(math.pi / 4)
            for chi, rho, z in dust:
                near(rho, a * math.sin(chi), f"Oppenheimer-Snyder cap rho at {chi}")
                near(z - dust[0][2], a * (1 - math.cos(chi)), f"Oppenheimer-Snyder cap z at {chi}")
        release = piece("oppenheimer_snyder", "exterior")
        for r, rho, z in release:
            near(z - release[0][2], 2 * math.sqrt(r - 1) - 2, f"Oppenheimer-Snyder release z at {r}")
        self.assertEqual(release[-1][2], 0, "the rim of the drawing stands at z = 0")
        # Tolman-Bondi's cloud at its release: every shell at its label, R = r, and outside it
        # Flamm's paraboloid of 2GM/c^2 = r_b/2 from the surface, drawn twice as tall.
        vertical = self.named["tolman_bondi"]["views"][0]["vertical"]
        cloud, outside = piece("tolman_bondi", "cloud"), piece("tolman_bondi", "exterior")
        for r, rho, z in cloud + outside:
            near(rho, r, f"Tolman-Bondi release rho at {r}")
        for r, rho, z in outside:
            near((z - outside[0][2]) / vertical, 2 * math.sqrt(0.5 * (r - 0.5)) - 1, f"Tolman-Bondi release z at {r}")
        # Bertotti-Robinson's equator is a cylinder of radius b, z = b ln r, and its sphere of radius b.
        for r, rho, z in piece("bertotti_robinson", "cylinder"):
            near(rho, 1.0, f"Bertotti-Robinson rho at {r}")
            near(z, math.log(r), f"Bertotti-Robinson z at {r}")
        for theta, rho, z in piece("bertotti_robinson", "sphere", view=1):
            near(rho, math.sin(theta), f"Bertotti-Robinson sphere rho at {theta}")
            near(z, 1 - math.cos(theta), f"Bertotti-Robinson sphere z at {theta}")
        # The Nariai universe's moments are cylinders of radius cosh(ct) and height 2 pi, a = 1, and
        # its sphere has radius a.
        for number, t in enumerate((-1.5, -0.75, 0.0, 0.75, 1.5)):
            for pid, z_of in (("near", lambda th: th), ("far", lambda th: 2 * math.pi - th)):
                for theta, rho, z in piece("nariai", pid, number):
                    near(rho, math.cosh(t), f"Nariai rho at ct = {t}, theta = {theta}")
                    near(z, z_of(theta), f"Nariai z at ct = {t}, theta = {theta}")
        # The Kantowski-Sachs moments are cylinders of radius b on which |r| <= 1 is 2a long: the dust
        # at a = 1 + eta tan(eta), b = cos^2(eta), and the vacuum at a = sqrt(1/T - 1), b = T.
        for number, eta in enumerate((0.0, 0.3, 0.6, 0.9, 1.2)):
            for r, rho, z in piece("kantowski_sachs", "tube", number):
                near(rho, math.cos(eta) ** 2, f"Kantowski-Sachs dust rho at eta = {eta}")
                near(z, (1 + eta * math.tan(eta)) * r, f"Kantowski-Sachs dust z at eta = {eta}, r = {r}")
        for number, T in enumerate((0.9, 0.7, 0.5, 0.3, 0.1)):
            for r, rho, z in piece("kantowski_sachs", "tube", number, view=1):
                near(rho, T, f"Kantowski-Sachs vacuum rho at T = {T}")
                near(z, math.sqrt(1 / T - 1) * r, f"Kantowski-Sachs vacuum z at T = {T}, r = {r}")
        for theta, rho, z in piece("nariai", "sphere", view=1):
            near(rho, math.sin(theta), f"Nariai sphere rho at {theta}")
            near(z, 1 - math.cos(theta), f"Nariai sphere z at {theta}")
        # The global monopole at Delta = 0.19: the Barriola-Vilenkin equator is the cone rho = 0.9 r,
        # z = sqrt(0.19) r, and Letelier's black hole at r_s = 1 rises from its throat at 100/81 as
        # z = (w sqrt(0.19 w^2 + 1) + arsinh(sqrt(0.19) w)/sqrt(0.19))/0.81^(3/2), w = sqrt(0.81 r - 1).
        for r, rho, z in piece("global_monopole", "cone"):
            near(rho, 0.9 * r, f"global monopole cone rho at {r}")
            near(z, math.sqrt(0.19) * r, f"global monopole cone z at {r}")
        for sign, pid in ((1, "exterior"), (-1, "other_exterior")):
            for r, rho, z in piece("global_monopole", pid, view=1):
                w = math.sqrt(max(0.81 * r - 1, 0))
                near(rho, r, f"Letelier rho at {r}")
                near(z, sign * (w * math.sqrt(0.19 * w * w + 1) + math.asinh(math.sqrt(0.19) * w) / math.sqrt(0.19))
                     / 0.81 ** 1.5, f"Letelier z at {r}")
        # The black hole on a cosmic string at b = 0.9 and r_s = 1: its equator rises from the throat as
        # z = w sqrt(1 + 0.19 w^2) + arsinh(sqrt(0.19) w)/sqrt(0.19), w = sqrt(r - 1), at rho = 0.9 r, and
        # its horizon is the spindle rho = 0.9 sin(theta), 0.9 r_s wide at its equator.
        for sign, pid in ((1, "exterior"), (-1, "other_exterior")):
            for r, rho, z in piece("string_black_hole", pid):
                w = math.sqrt(max(r - 1, 0))
                near(rho, 0.9 * r, f"threaded black hole rho at {r}")
                near(z, sign * (w * math.sqrt(1 + 0.19 * w * w) + math.asinh(math.sqrt(0.19) * w) / math.sqrt(0.19)),
                     f"threaded black hole z at {r}")
        spindle = piece("string_black_hole", "horizon", view=1)
        for theta, rho, z in spindle:
            near(rho, 0.9 * math.sin(theta), f"threaded black hole horizon rho at {theta}")
        self.assertLess(abs(sum(math.pi * (a[1] + b[1]) * math.hypot(b[1] - a[1], b[2] - a[2])
                                for a, b in zip(spindle, spindle[1:])) - 3.6 * math.pi), 2e-3)
        # Van Stockum's circles grow to r = R/sqrt(2) and shrink after, and the drawing stops at 0.83 R.
        dust = piece("stockum_dust", "dust")
        for r, rho, z in dust:
            near(rho, r * math.sqrt(1 - r * r), f"van Stockum rho at {r}")
        self.assertAlmostEqual(dust[-1][0], 0.8337, places=4)
        # Taub-NUT's equator has circles of radius sqrt(r^2 + l^2), from the horizon r+ = m + sqrt(m^2 + l^2).
        nut = piece("taub_nut", "exterior")
        self.assertAlmostEqual(nut[0][0], 1 + math.sqrt(1.25), places=12)
        for r, rho, z in nut:
            near(rho, math.sqrt(r * r + 0.25), f"Taub-NUT rho at {r}")
        # Godel's circles about one world line, which the drawing follows to sinh^2 r = 1/sqrt(2).
        dust = piece("godel", "dust")
        for r, rho, z in dust:
            s = math.sinh(r)
            near(rho, math.sqrt(2) * s * math.sqrt(max(1 - s * s, 0)), f"Godel rho at {r}")
        self.assertAlmostEqual(dust[-1][0], math.asinh(2 ** -0.25), places=12)
        # Som and Raychaudhuri's circles about one world line, of radius r sqrt(1 - r^2) in units of
        # c/Omega, which the drawing follows to r = sqrt(3)/2, at the height the quadrature of
        # dz/dr = r sqrt((3 - 4r^2)/(1 - r^2)) has in closed form.
        dust = piece("som_raychaudhuri", "dust")
        for r, rho, z in dust:
            near(rho, r * math.sqrt(1 - r * r), f"Som-Raychaudhuri rho at {r}")
            w = math.sqrt(max(3 - 4 * r * r, 0) / (4 - 4 * r * r))
            near(z, (math.sqrt(3) - math.sqrt((1 - r * r) * max(3 - 4 * r * r, 0))) / 2
                 - (math.atanh(math.sqrt(3) / 2) - math.atanh(w)) / 4, f"Som-Raychaudhuri z at {r}")
        self.assertAlmostEqual(dust[-1][0], math.sqrt(3) / 2, places=12)
        self.assertAlmostEqual(max(rho for _, rho, _ in dust), 0.5, places=6)
        for r, rho, z in piece("cosmic_string", "exterior", view=1):
            near(rho, 0.9 * r, f"cone rho at {r}")
            near(z, math.sqrt(0.19) * r, f"cone z at {r}")
        for r, rho, z in piece("cosmic_string", "cone"):
            near(rho, 0.9 * r, f"the ideal string's cone rho at {r}")
            near(z, math.sqrt(0.19) * r, f"the ideal string's cone z at {r}")
        core = piece("cosmic_string", "core", view=1)
        for chi, rho, z in core:
            near(rho, math.sin(chi), f"Gott rho at {chi}")
            near(z - core[0][2], 1 - math.cos(chi), f"Gott z at {chi}")
        # The point particle at alpha = 3/4: the cone outside Gott and Alpert's planet and down to
        # its apex, and the planet, a cap of a unit sphere out to cos(chi_0) = 3/4, where they meet.
        for pid in ("exterior", "apex"):
            for r, rho, z in piece("point_particle_2plus1", pid, view=1):
                near(rho, 0.75 * r, f"the particle's cone rho at {r}")
                near(z, math.sqrt(1 - 0.75 ** 2) * r, f"the particle's cone z at {r}")
        planet = piece("point_particle_2plus1", "core", view=1)
        for chi, rho, z in planet:
            near(rho, math.sin(chi), f"the planet's rho at {chi}")
            near(z - planet[0][2], 1 - math.cos(chi), f"the planet's z at {chi}")
        self.assertAlmostEqual(planet[-1][0], math.acos(0.75), places=12)
        outside = piece("point_particle_2plus1", "exterior", view=1)
        self.assertAlmostEqual(outside[0][0], math.tan(math.acos(0.75)), places=12)
        near(planet[-1][1], outside[0][1], "the planet's edge and the cone's are one circle")
        near(planet[-1][2], outside[0][2], "at one height")
        for number in range(len(self.named["frw"]["views"][0]["surfaces"])):
            hemisphere = piece("frw", "near", number)
            a = -hemisphere[0][2]
            for sign, pid in ((-1, "near"), (1, "far")):
                for r, rho, z in piece("frw", pid, number):
                    near(rho, a * r, f"sphere rho at {r}")
                    near(z, sign * a * math.sqrt(max(1 - r * r, 0)), f"sphere z at {r}")

    @checks("aichelburg_sexl", "anti_de_sitter", "bianchi", "kasner", "light_beam", "malament_hogarth", "milne",
            "minkowski", "mixmaster", "natario", "pp_wave")
    def test_the_eleven_drawn_last_are_the_surfaces_known_in_closed_form(self):
        """The plane, the hyperboloid, the tube of the Malament-Hogarth toy, Taub's round
        moment, and the flat planes with what is marked on them, from the numbers written."""
        def view(metric_id):
            return self.named[metric_id]["views"][0]

        def piece(surface, piece_id):
            return [tuple(point) for point in next(p for p in surface["pieces"] if p["id"] == piece_id)["points"]]

        def near(a, b, where, tol=2e-6):
            self.assertLess(abs(a - b), tol, where)
        # The flat plane: every point level, and rho the distance along the profile.
        for metric_id, piece_id in (("minkowski", "plane"),):
            for x, rho, z in piece(view(metric_id)["surfaces"][0], piece_id):
                near(rho, x, f"{metric_id} rho at {x}")
                self.assertEqual(z, 0, f"{metric_id} z at {x}")
        # Anti-de Sitter's static equator on the hyperboloid Z = sqrt(L^2 + r^2) - L in Minkowski
        # space, and the light cone it nears, Z = rho - L.
        ads = view("anti_de_sitter")
        for p in ads["surfaces"][0]["pieces"]:
            self.assertEqual(p.get("space"), "minkowski", f"anti-de Sitter {p['id']}")
        for r, rho, z in piece(ads["surfaces"][0], "sheet"):
            near(rho, r, f"anti-de Sitter rho at {r}")
            near(z, math.sqrt(1 + r * r) - 1, f"anti-de Sitter z at {r}")
        for r, rho, z in piece(ads["surfaces"][0], "cone"):
            near(z, rho - 1, f"the light cone at {r}")
        # The Milne universe's moments on the hyperboloids Z = ct cosh chi, rho = ct sinh chi, the
        # moment itself in its inertial chart, nested inside the light cone Z = rho of their apex.
        milne = view("milne")
        for surface in milne["surfaces"] + milne["movie"]["frames"]:
            for p in surface["pieces"]:
                self.assertEqual(p.get("space"), "minkowski", f"Milne {p['id']} at ct = {surface.get('time', surface.get('value'))}")
            t = surface.get("time", surface.get("value"))
            for chi, rho, z in piece(surface, "sheet"):
                near(rho, t * math.sinh(chi), f"Milne rho at ct = {t}, chi = {chi}")
                near(z, t * math.cosh(chi), f"Milne z at ct = {t}, chi = {chi}")
            for chi, rho, z in piece(surface, "cone"):
                near(z, rho, f"Milne's light cone at ct = {t}")
        # The Malament-Hogarth plane at ct = 0: the circle through s has radius s Omega, s plus
        # e^(1 - 1/(1 - s^2)), which closes in on 1 down the tube, and it is flat beyond s = 1.
        tube = view("malament_hogarth")["surfaces"][-1]
        for s, rho, z in piece(tube, "well"):
            near(rho, s + (math.exp(1 - 1 / (1 - s * s)) if s < 1 else 0.0), f"the tube at {s}", 1e-5)
        self.assertEqual(piece(tube, "well")[-1][0], 1.0, "the well ends on the edge of the unit ball")
        flat = piece(tube, "flat")
        for s, rho, z in flat:
            near(rho, s, f"beyond the unit ball at {s}")
            self.assertEqual(z, flat[0][2], f"beyond the unit ball at {s}")
        # Taub's round moment: the great sphere of radius 2a, rho = 2a sin(theta/2) and
        # z = -2a cos(theta/2) on the near hemisphere; at every moment both hemispheres end at
        # the equator theta = pi.
        taub = view("mixmaster")["surfaces"]
        hemisphere = piece(taub[2], "near")
        a = hemisphere[-1][1] / 2
        for theta, rho, z in hemisphere:
            near(rho, 2 * a * math.sin(theta / 2), f"Taub's round moment rho at {theta}")
            near(z, -2 * a * math.cos(theta / 2), f"Taub's round moment z at {theta}")
        for surface in taub:
            self.assertEqual(piece(surface, "near")[-1][0], math.pi)
            self.assertEqual(piece(surface, "far")[-1][0], math.pi)
        # Kasner's ring at t: the ellipse reaching t^(-2/7) along x and t^(6/7) along z. The dust's
        # and the wave's rings are ellipses on the axes, and the wave's last a segment on x.
        def rings(metric_id):
            return next(v for v in self.named[metric_id]["views"] if v["id"] == "ring")
        for surface in rings("kasner")["surfaces"]:
            t = surface["time"]
            for X, Y, Z in surface["curves"][0]["points"]:
                near((X / t ** (-2 / 7)) ** 2 + (Y / t ** (6 / 7)) ** 2, 1, f"Kasner's ring at t = {t}", 1e-5)
        for metric_id in ("bianchi", "pp_wave"):
            for surface in rings(metric_id)["surfaces"]:
                P = surface["curves"][0]["points"]
                A, B = max(abs(p[0]) for p in P), max(abs(p[1]) for p in P)
                for X, Y, Z in P:
                    near((X / A) ** 2 + ((Y / B) ** 2 if B else 1 - (X / A) ** 2), 1, f"{metric_id}'s ring", 1e-5)
        self.assertEqual(max(abs(p[1]) for p in rings("pp_wave")["surfaces"][-1]["curves"][0]["points"]), 0)
        # The shock's ring stays a circle about the axis, of radius 1 - u/2 behind the shock to the
        # width of the declared pulse.
        for surface in rings("aichelburg_sexl")["surfaces"]:
            P = surface["curves"][0]["points"]
            radii = [math.hypot(X, Y) for X, Y, _ in P]
            near(max(radii) - min(radii), 0, "the Aichelburg-Sexl ring is a circle", 1e-5)
            near(radii[0], 1 - max(surface["time"], 0) / 2, "the Aichelburg-Sexl ring's radius", 2e-3)
        # The ring about Bonnor's uniform beam, at rest on rho = 2R at u = 0, stays a circle. Outside the
        # beam rho'^2 = ln(2/rho)/2, so it stands at rho when u = 2 sqrt(2 pi) erf(sqrt(ln(2/rho))); inside,
        # rho'' = -rho/4 from the edge, reached at u_1 with the speed v_1 = sqrt(ln(2)/2), so its radius is
        # |cos(s/2) - 2 v_1 sin(s/2)| at s = u - u_1. It closes on the axis at u_f = u_1 + 2 arctan(1/(2 v_1))
        # and opens again as it closed, its radius at u the radius at 2 u_f - u.
        u_1, v_1 = 2 * math.sqrt(2 * math.pi) * math.erf(math.sqrt(math.log(2))), math.sqrt(math.log(2) / 2)
        u_f = u_1 + 2 * math.atan(1 / (2 * v_1))
        for surface in rings("light_beam")["surfaces"]:
            P = surface["curves"][0]["points"]
            radii = [math.hypot(X, Y) for X, Y, _ in P]
            near(max(radii) - min(radii), 0, "the ring about the beam of light is a circle", 1e-5)
            u = min(surface["time"], 2 * u_f - surface["time"])
            if u < u_1:
                near(2 * math.sqrt(2 * math.pi) * math.erf(math.sqrt(math.log(2 / radii[0]))), u,
                     "the ring falling toward the beam of light", 1e-5)
            else:
                near(radii[0], abs(math.cos((u - u_1) / 2) - 2 * v_1 * math.sin((u - u_1) / 2)),
                     "the ring inside the beam of light", 1e-5)

        # Natario's lines of flow, closed, each on one level of the stream function n(r_s) y^2,
        # y being the drawing's Y, and each on the triangles of the height.
        n = lambda r: (math.tanh(4 * (r + 1)) - math.tanh(4 * (r - 1))) / (4 * math.tanh(4))  # noqa: E731
        natario = view("natario")["surfaces"][0]
        flows = [c for c in natario["curves"] if c["class"] == "flow"]
        self.assertEqual(len(flows), 6)
        for curve in flows:
            self.assertTrue(curve["closed"])
            levels = [n(math.hypot(X, Y)) * Y * Y for X, Y, Z in curve["points"]]
            self.assertLess(max(levels) - min(levels), 1e-6 * max(levels), "a line of Natario's flow")
            for X, Y, Z in curve["points"][::7]:
                self.assertLess(abs(grid_height(natario["pieces"][0], X, Y) - Z), 1e-6,
                                f"a line of Natario's flow at {X}, {Y}")

    def test_every_caption_names_what_is_drawn(self):
        for metric_id, data in self.embedding.items():
            for view in data["views"]:
                self.assertTrue(view["caption"], f"{metric_id} {view['id']}")
                self.assertTrue(opens_with_a_noun_phrase(view["caption"][0]), f"{metric_id} {view['id']}")

    def test_every_text_is_tex_with_its_mathematics_closed(self):
        for name, data in self.embedding.items():
            for field, value in embedding_prose(name, data):
                self.assertTrue(value.strip(), field)
                self.assertEqual(value.replace("\\$", "").count("$") % 2, 0, field)
                self.assertEqual(value.count("{"), value.count("}"), field)

    def test_every_figure_draws_inside_its_box_and_names_only_what_it_draws(self):
        for name, data in self.embedding.items():
            for view in data["views"]:
                figure, where = view["figure"], f"{name} {view['id']}"
                x0, x1, y0, y1 = figure["box"]
                self.assertTrue(x0 < x1 and y0 < y1, where)
                drawn = {}
                for layer in figure["layers"]:
                    self.assertIn(layer["kind"], {"fill", "line", "point"}, where)
                    drawn.setdefault(layer["class"].removesuffix("-far"), layer["kind"])
                    rings = ([[layer["at"]]] if layer["kind"] == "point"
                             else [layer["points"]] + layer.get("holes", []))
                    for ring in rings:
                        for x, y in ring:
                            self.assertTrue(x0 <= x <= x1 and y0 <= y <= y1, f"{where} {layer['class']} at {x}, {y}")
                for label in figure["labels"]:
                    x, y = label["at"]
                    self.assertTrue(x0 <= x <= x1 and y0 <= y <= y1, f"{where} label {label['text']}")
                # A movie's legend names what any of its frames marks.
                for frame in view.get("movie", {}).get("frames", []):
                    for mark in frame["rings"] + frame.get("curves", []):
                        drawn.setdefault(mark["class"], "line")
                for kind, cls, _ in figure["legend"]:
                    self.assertEqual(drawn.get(cls), kind, f"{where} legend {cls}")

    def test_every_class_a_figure_paints_is_styled_on_the_page_and_in_print(self):
        page = (build.ROOT / "_layouts" / "mfs.html").read_text(encoding="utf-8")
        for name, data in self.embedding.items():
            for view in data["views"]:
                for layer in view["figure"]["layers"]:
                    selector = "em-" + layer["class"] + ("-fill" if layer["kind"] == "fill" else "")
                    self.assertRegex(page, r"\n    [^\n{]*\." + re.escape(selector) + r"(?![\w-])[^{\n]*\{",
                                     f"{name} {view['id']} paints {selector}, which is not styled")
                    self.assertRegex(page, r"\.mfs-print-body \." + re.escape(selector) + r"(?![\w-])",
                                     f"{name} {view['id']} paints {selector}, which is not styled in print")

    def grid(self, metric_id):
        view = self.named[metric_id]["views"][0]
        piece = view["surfaces"][0]["pieces"][0]
        self.assertIn("grid", piece, metric_id)
        return view, piece

    @checks("alcubierre")
    def test_alcubierre_draws_the_expansion_as_a_height(self):
        """theta = v_s (x - x_s)/r_s df/dr_s of the observers who ride the slices, with v_s = 2 and
        Alcubierre's profile at R = 1 and sigma = 4, drawn as the height theta R^2/4c over the plane
        of the path, the ship heading toward +X: negative ahead of the ship, where space contracts,
        and positive behind it, where it expands, as Alcubierre drew it in 1994."""
        view, piece = self.grid("alcubierre")
        grid = piece["grid"]
        self.assertEqual(grid["frame"], "polar")

        def theta(X, Y):
            r = math.hypot(X, Y)
            if r == 0:
                return 0.0
            df = 4 * (1 / math.cosh(4 * (r + 1)) ** 2 - 1 / math.cosh(4 * (r - 1)) ** 2) / (2 * math.tanh(4))
            return 2 * X / r * df
        nodes = [P for row in grid_nodes(piece) for P in row]
        # Every node stands at theta R^2/4c, to the rounding of the file.
        for X, Y, Z in nodes:
            self.assertLess(abs(Z - theta(X, Y) / 4), 1e-6, f"Alcubierre's height at {X}, {Y}")
        # Space contracts ahead of the ship and expands behind it, everywhere on the plane.
        for X, Y, Z in nodes:
            if abs(X) > 1e-9 and abs(Z) > 1e-6:
                self.assertEqual(math.copysign(1, Z), -math.copysign(1, X), f"the sign of theta at {X}, {Y}")
        # At chosen points, in units of c/R.
        u, v = grid["u"], grid["v"]
        for r, phi, want in ((1.0, math.pi, 4.002683), (1.0, 0.0, -4.002683), (0.5, math.pi, 0.282695),
                             (0.5, 0.0, -0.282695), (1.0, math.pi / 4, -2.830324), (2.0, 0.0, -0.005367),
                             (1.0, math.pi / 2, 0.0)):
            i = min(range(len(u)), key=lambda k: abs(u[k] - r))
            j = min(range(len(v)), key=lambda k: abs(v[k] - phi))
            self.assertAlmostEqual(u[i], r, places=12)
            self.assertAlmostEqual(v[j], phi, places=12)
            self.assertLess(abs(4 * grid["z"][i][j] - want), 1e-5, f"theta at r_s = {r}, phi = {phi}")
        # The greatest expansion, 4.0027 c/R one R behind the ship, is the highest node, and the
        # fastest contraction the lowest, one R ahead.
        top = max(nodes, key=lambda P: P[2])
        bottom = min(nodes, key=lambda P: P[2])
        self.assertAlmostEqual(top[2], 1.000671, places=5)
        self.assertAlmostEqual(bottom[2], -1.000671, places=5)
        self.assertLess(top[0], 0)
        self.assertGreater(bottom[0], 0)
        # The level lines lie at half the greatest expansion and contraction, the one ahead of
        # the ship and the other behind it; the circle v_s f = 1 at 1.0001676 R.
        curves = {c["class"]: c for c in view["surfaces"][0]["curves"]}
        for cls, sign in (("contract", -1), ("expand", 1)):
            for X, Y, Z in curves[cls]["points"]:
                self.assertAlmostEqual(Z, sign * 0.500335, places=5, msg=f"Alcubierre's {cls} at {X}, {Y}")
                self.assertEqual(math.copysign(1, X), -sign, f"Alcubierre's {cls} at {X}, {Y}")
                self.assertLess(abs(theta(X, Y) / 4 - Z), 0.015, f"Alcubierre's {cls} at {X}, {Y}")
        # The circle runs from node to node of the grid's circle there, on its 72 chords.
        for X, Y, Z in curves["wall"]["points"]:
            r = math.hypot(X, Y)
            self.assertTrue(1.0001676 * math.cos(math.pi / 72) - 1e-6 <= r <= 1.0001676 + 1e-6,
                            f"the circle v_s f = 1 at {X}, {Y}")
        self.assertIn("4c/R", view["height"])

    @checks("natario")
    def test_natario_draws_the_energy_density_as_a_height(self):
        """epsilon = -c^4 K_ij K^ij/16 pi G of the observers who ride the slices, with the zero
        expansion field of Natario's drive for Alcubierre's profile at R = 1 and sigma = 4, n = f/2
        and v_s = 2, drawn as the height G^tt R^2/16 = -K_ij K^ij R^2/32 over the plane of the
        path: zero at the ship and far outside, negative everywhere else, and deepest in the
        wall beside the ship. K_ij is taken here from the declared field by differences, with
        no sympy, and the height checked against it."""
        view, piece = self.grid("natario")
        grid = piece["grid"]
        self.assertEqual(grid["frame"], "polar")

        def n(r):
            return (math.tanh(4 * (r + 1)) - math.tanh(4 * (r - 1))) / (4 * math.tanh(4))

        def dn(r):
            return (1 / math.cosh(4 * (r + 1)) ** 2 - 1 / math.cosh(4 * (r - 1)) ** 2) / math.tanh(4)

        def field(x, y, z):
            rho = math.sqrt(x * x + y * y + z * z)
            return (2 * (2 * n(rho) + rho * dn(rho) - dn(rho) * x * x / rho),
                    -2 * dn(rho) * x * y / rho, -2 * dn(rho) * x * z / rho)

        def height(X, Y, h=1e-5):
            if math.hypot(X, Y) < 1e-3:
                return 0.0
            d = []
            for axis in range(3):
                e = [0.0, 0.0, 0.0]
                e[axis] = h
                a, b = field(X + e[0], Y + e[1], e[2]), field(X - e[0], Y - e[1], -e[2])
                d.append([(p - q) / (2 * h) for p, q in zip(a, b)])
            KK = sum(((d[i][j] + d[j][i]) / 2) ** 2 for i in range(3) for j in range(3))
            self.assertLess(abs(d[0][0] + d[1][1] + d[2][2]), 1e-6, f"the divergence at {X}, {Y}")
            return -KK / 32
        nodes = [P for row in grid_nodes(piece) for P in row]
        for X, Y, Z in nodes:
            self.assertLess(abs(Z - height(X, Y)), 1e-5, f"Natario's height at {X}, {Y}")
            self.assertLessEqual(Z, 0, f"Natario's height at {X}, {Y}")
        # The deepest point, -1.0614 R, beside the ship at r_s = 0.878 R across the path, is
        # below every node, and a node lies within a few thousandths of it.
        deepest = min(nodes, key=lambda P: P[2])
        self.assertGreater(deepest[2], -1.06139)
        self.assertLess(deepest[2], -1.05)
        self.assertLess(abs(deepest[0]), 1e-9)
        self.assertLess(abs(abs(deepest[1]) - 0.878), 0.13)
        u, v, z = grid["u"], grid["v"], grid["z"]
        self.assertEqual(z[0], [0.0] * len(v), "the ship's own place is level")
        self.assertLess(max(abs(h) for h in z[-1]), 1e-6, "the rim at 3R is level")
        on_path = min(z[u.index(1.0)][v.index(0.0)], z[u.index(1.0)][min(range(len(v)), key=lambda j: abs(v[j] - math.pi))])
        self.assertAlmostEqual(on_path, -12.016102 / 16, places=5)
        curves = {c["class"]: c for c in view["surfaces"][0]["curves"]}
        for X, Y, Z in curves["wall"]["points"]:
            r = math.hypot(X, Y)
            self.assertTrue(math.cos(math.pi / 72) - 1e-6 <= r <= 1 + 1e-6, f"the circle r_s = R at {X}, {Y}")
        path = curves["path"]["points"]
        self.assertEqual((path[0][0], path[-1][0]), (-3.0, 3.0))
        self.assertTrue(all(abs(Y) < 1e-12 for _, Y, _ in path))
        self.assertIn("K_{ij}K^{ij}", view["height"])

    def test_every_embedding_diagram_stands_in_relief_but_the_planes_of_flat_slices(self):
        """No drawing is a bare flat plane: every surface of revolution rises or falls by at least
        RELIEF of its width and every height over a plane by RELIEF of its extent, save the flat
        planes that carry what their spacetime does on them, Minkowski's and the rings of free
        particles. Lentz's class has flat slices for every potential and no soliton that can be
        computed, so it has no embedding diagram."""
        RELIEF = 0.05
        # The spacetimes whose embedding is the flat plane that carries what they do, and the single
        # moments that are flat planes, from the `flat` and `flat_moments` of their exceptions.
        flat = {name for name, data in EXCEPTIONS.items() if "flat" in data}
        flat_moments = {(name, view, k) for name, group in exceptions("flat_moments") for view, k in group["moments"]}
        self.assertNotIn("lentz", self.embedding)
        self.assertNotIn("embedding", next(m for m in read(build.INDEX_FILE) if m["id"] == "lentz"))
        for name, data in self.embedding.items():
            if name in flat:
                continue
            for view in data["views"]:
                for k, surface in enumerate(view["surfaces"]):
                    where = f"{name} {view['id']} surface {k}"
                    if (name, view["id"], k) in flat_moments:
                        continue
                    heights, widths = [], []
                    for piece in surface["pieces"]:
                        if piece.get("reference"):
                            continue
                        if "grid" in piece:
                            nodes = [P for row in grid_nodes(piece) for P in row]
                            heights += [P[2] for P in nodes]
                            widths.append(max(max(P[0] for P in nodes) - min(P[0] for P in nodes),
                                              max(P[1] for P in nodes) - min(P[1] for P in nodes)))
                        else:
                            heights += [P[2] for P in piece["points"]]
                            widths.append(2 * max(P[1] for P in piece["points"]))
                    self.assertTrue(all(math.isfinite(h) for h in heights), f"{where}: a height is not a number")
                    self.assertGreater(max(widths), 0, f"{where}: no width")
                    self.assertGreaterEqual((max(heights) - min(heights)) / max(widths), RELIEF,
                                            f"{where}: flatter than {RELIEF} of its width")

    def test_a_surface_drawn_taller_than_its_embedding_says_so_in_its_caption(self):
        """A view whose heights are drawn `vertical` times as tall as the embedding's says so in
        the first sentence of its caption, "(vertical scale $\\times N$)", and only Tolman-Bondi's
        cloud and Szekeres's are: at the scale the metric gives it Tolman-Bondi's cloud rises by a
        third of its width, and drawn twice as tall each moment stands in relief of at least 0.6 of
        its width; the dish of Szekeres's cloud is 4/15 deep at ct = 0 and is drawn three times as deep."""
        scaled = {name: view for name, data in self.embedding.items() for view in data["views"] if "vertical" in view}
        self.assertEqual(set(scaled), set(filter(chosen, {"tolman_bondi", "szekeres"})))
        if "szekeres" in scaled:
            self.assertEqual(scaled["szekeres"]["vertical"], 3)
            release = scaled["szekeres"]["surfaces"][1]
            self.assertEqual(release["time"], 0)
            for r, rho, z in release["pieces"][0]["points"]:
                self.assertAlmostEqual(rho, r, 6, f"Szekeres at ct = 0, rho at {r}")
                self.assertAlmostEqual(z / 3, 2 * r ** 3 / 3 - 2 * r ** 5 / 5 - 4 / 15, 6, f"Szekeres at ct = 0, z at {r}")
            for surface in scaled["szekeres"]["surfaces"]:
                self.assertTrue(all(P[2] == 0 for P in surface["pieces"][1]["points"]), "Szekeres's exterior is a plane")
        for name, view in scaled.items():
            self.assertIsInstance(view["vertical"], int, name)
            self.assertGreater(view["vertical"], 1, name)
            first = view["caption"][0].split(". ")[0]
            self.assertIn(f"(vertical scale $\\times {view['vertical']}$)", first, name)
        for data in self.embedding.values():
            for view in data["views"]:
                if "vertical" not in view:
                    self.assertNotIn("vertical scale", " ".join(view["caption"]), view["id"])
        if "tolman_bondi" in scaled:
            view = scaled["tolman_bondi"]
            self.assertEqual(view["vertical"], 2)
            for surface in view["surfaces"]:
                points = [P for piece in surface["pieces"] for P in piece["points"]]
                width = 2 * max(P[1] for P in points)
                relief = (max(P[2] for P in points) - min(P[2] for P in points)) / width
                self.assertGreaterEqual(relief, 0.6, f"Tolman-Bondi at ct = {surface['time']}")
                cloud = next(piece for piece in surface["pieces"] if piece["id"] == "cloud")["points"]
                self.assertGreaterEqual((cloud[-1][2] - cloud[0][2]) / (2 * cloud[-1][1]), 0.55,
                                        f"Tolman-Bondi's cloud at ct = {surface['time']}")

    @checks("krasnikov")
    def test_krasnikov_draws_how_far_the_tube_tips_the_light_cone_as_a_height(self):
        """1 - k, twice the published g_tx, of the tube the spacetime diagram declares, at
        ct = 5 rho_0, drawn over the plane of the tube's axis: 0 outside the tube, 1 where k = 0
        and 1.7977 on the axis in the middle of it."""
        view, piece = self.grid("krasnikov")
        grid = piece["grid"]
        self.assertEqual(grid["frame"], "cartesian")

        def step(q):
            return (1 + math.tanh(q / 0.15)) / 2

        def one_minus_k(X, Y):
            x = X + 2
            return 1.8 * step((1 - Y * Y) / 2) * step(5 - x) * step(x) * step(4 - x)
        nodes = [P for row in grid_nodes(piece) for P in row]
        for X, Y, Z in nodes:
            self.assertLess(abs(Z - one_minus_k(X, Y)), 1e-6, f"Krasnikov's height at {X}, {Y}")
            self.assertGreaterEqual(Z, 0, f"Krasnikov's height at {X}, {Y}")
            self.assertLess(Z, 1.8, f"Krasnikov's height at {X}, {Y}")
        u, v, z = grid["u"], grid["v"], grid["z"]
        for X, Y, want in ((0.0, 0.0, 1.7977122), (0.0, 1.0, 0.9), (0.0, 0.5, 1.7879529), (-2.0, 0.0, 0.8988561),
                           (2.0, 0.0, 0.8988546), (-2.5, 0.0, 0.0022849), (2.5, 0.0, 0.0022820), (0.0, 1.5, 0.0004326),
                           (0.0, 2.0, 0.0)):
            i, j = u.index(X), v.index(Y)
            self.assertLess(abs(z[i][j] - want), 1e-6, f"1 - k at x = {X + 2}, r = {abs(Y)}")
        # The ridge stands all along the path from x = 0 to D at the moment drawn: the tube is
        # complete, and flat beyond both ends and outside its wall.
        for X in (-1.5, -1.0, 0.0, 1.0, 1.5):
            self.assertGreater(z[u.index(X)][v.index(0.0)], 1.79, f"the ridge at x = {X + 2}")
        for X in (-3.0, 3.0):
            self.assertLess(z[u.index(X)][v.index(0.0)], 3e-6, f"beyond the tube at x = {X + 2}")
        curves = {c["class"]: c for c in view["surfaces"][0]["curves"]}
        for X, Y, Z in curves["wall"]["points"]:
            self.assertAlmostEqual(Z, 1.0, places=6, msg=f"k = 0 at {X}, {Y}")
            self.assertLess(abs(one_minus_k(X, Y) - 1), 0.02, f"k = 0 at {X}, {Y}")
        path = curves["path"]["points"]
        self.assertEqual((path[0][0], path[-1][0]), (-2.0, 2.0))
        self.assertTrue(all(Y == 0 for _, Y, _ in path))
        self.assertIn("1 - k", view["height"])

    def test_a_view_says_what_its_height_is_exactly_when_it_draws_a_grid(self):
        for name, data in self.embedding.items():
            for view in data["views"]:
                grids = any("grid" in p for s in view["surfaces"] for p in s["pieces"])
                self.assertEqual("height" in view, grids, f"{name} {view['id']}")
        self.assertEqual({name for name, data in self.embedding.items()
                          if any("height" in view for view in data["views"])},
                         {name for name, data in EXCEPTIONS.items() if "height" in data and chosen(name)})

    def test_a_grid_that_is_not_one_is_refused(self):
        def spoil(change, words):
            data = copy.deepcopy(self.named["krasnikov"])
            change(data["views"][0]["surfaces"][0]["pieces"][0])
            with self.assertRaises(build.DataError) as raised:
                self.load_folder({"krasnikov.json": data}, self.metrics)
            self.assertIn(words, str(raised.exception))
        spoil(lambda p: p["grid"]["u"].reverse(), "one way along its u")
        spoil(lambda p: p["grid"]["z"].pop(), "a height at every node")
        spoil(lambda p: p["grid"].__setitem__("frame", "spherical"), "neither a polar nor a Cartesian")
        spoil(lambda p: p.__setitem__("points", []), "a grid and a profile at once")
        spoil(lambda p: p.pop("edge"), "beyond its edge")

    def load_folder(self, files, metrics):
        with tempfile.TemporaryDirectory() as folder:
            for filename, data in files.items():
                (Path(folder) / filename).write_text(json.dumps(data), encoding="utf-8")
            with mock.patch.object(build, "EMBEDDING_DIR", Path(folder)):
                return build.load_embedding(metrics)


_turn_check = None


def turn_check(test):
    """One run of the page's own turning geometry in Node over every published figure that
    turns of the spacetimes this run checks, _tools/turn_check.cjs, shared by every test that
    reads it."""
    global _turn_check
    if shutil.which("node") is None:
        test.skipTest("Node is not installed, so the page's geometry cannot be run")
    if _turn_check is None:
        run = subprocess.run(["node", str(build.ROOT / "_tools" / "turn_check.cjs"), *sorted(CHOSEN)],
                             capture_output=True, text=True, timeout=600)
        test.assertEqual(run.returncode, 0, run.stderr[-2000:])
        _turn_check = json.loads(run.stdout)
    return _turn_check


class TurningEmbeddingDiagrams(unittest.TestCase):
    """The page turns an embedding diagram with MFS/assets/turn.js, which draws it again from the
    view's surfaces by what the figure's `turn`, its labels' `ring` and `clear` and its layers'
    `flat` say, as _tools/README.md, "Turning the figure", defines them. At the figure's own
    camera it must give back the published figure, and at any other it must keep inside the box,
    draw a surface of revolution the same from every side and hide what lies behind."""

    def setUp(self):
        self.embedding = embedding_files()
        self.views = [(name, view) for name, data in self.embedding.items() for view in data["views"]]
        found(self, self.views, "embedding view")

    def check(self):
        return turn_check(self)

    def test_every_figure_carries_what_turning_it_needs(self):
        page = (build.ROOT / "_layouts" / "mfs.html").read_text(encoding="utf-8")
        for name, view in self.views:
            figure, where = view["figure"], f"{name} {view['id']}"
            turn = figure["turn"]
            # A movie draws one frame at a time, all at one place.
            self.assertEqual(len(turn["origins"]), 1 if "movie" in view else len(view["surfaces"]), where)
            for origin in turn["origins"]:
                self.assertEqual(len(origin), 2, where)
            self.assertIsInstance(turn["meridians"], int, where)
            self.assertGreater(turn["meridians"], 0, where)
            classes = {p["class"] for s in view["surfaces"] for p in s["pieces"]}
            for cls, fill in turn["tint"].items():
                self.assertIn(cls, classes, where)
                self.assertIn(fill, {"cover", "star"}, where)
                self.assertRegex(page, r"\.em-" + fill + r"-fill\b", where)
            for mark in turn["marks"]:
                surface = view["surfaces"][mark["surface"]]
                self.assertIn(mark["piece"], {p["id"] for p in surface["pieces"]}, where)
                self.assertRegex(page, r"\.em-" + re.escape(mark["class"]) + r"\b", where)
            # Every grid piece names the lines of its grid the figure draws, by index, a stack of
            # ellipses in its own grid, since each frame of a movie holds its own rows.
            grids = {(k, p["id"]): p["grid"] for k, s in enumerate(view["surfaces"]) for p in s["pieces"]
                     if "grid" in p and p["grid"]["frame"] != "ellipses"}
            for surface in view["surfaces"] + view.get("movie", {}).get("frames", []):
                for p in surface["pieces"]:
                    if p.get("grid", {}).get("frame") != "ellipses":
                        continue
                    self.assertTrue(p["grid"]["lines"], where)
                    for line in p["grid"]["lines"]:
                        self.assertTrue(all(0 <= i < len(p["grid"]["u"]) for i in line["u"]), where)
                        self.assertTrue(all(0 <= j < len(p["grid"]["v"]) for j in line["v"]), where)
                        self.assertRegex(page, r"\.em-" + re.escape(line["class"]) + r"\b", where)
            self.assertEqual({(g["surface"], g["piece"]) for g in turn.get("grid", [])}, set(grids), where)
            for g in turn.get("grid", []):
                grid = grids[(g["surface"], g["piece"])]
                self.assertTrue(all(0 <= i < len(grid["u"]) for i in g["u"]), where)
                self.assertTrue(all(0 <= j < len(grid["v"]) for j in g["v"]), where)
                self.assertRegex(page, r"\.em-" + re.escape(g["class"]) + r"\b", where)
            for label in figure["labels"]:
                if "ring" in label:
                    ring = label["ring"]
                    self.assertLess(ring["ring"], len(view["surfaces"][ring["surface"]]["rings"]), where)
                    self.assertIn(ring["side"], (1, -1), where)
                if "clear" in label:
                    self.assertIn("ring", label, where)
                    self.assertGreater(label["clear"], 0, where)
            for layer in figure["layers"]:
                self.assertIn(layer.get("flat", False), (False, True), where)
        # Only a cone laid flat lies in the plane of the page: the cosmic string's and the point particle's.
        flat = {name for name, view in self.views if any(layer.get("flat") for layer in view["figure"]["layers"])}
        self.assertEqual(flat, set(filter(chosen, {"cosmic_string", "point_particle_2plus1"})))

    def test_a_label_that_names_a_circle_stands_beside_its_end(self):
        # The end of a circle on the page is (+-rho, z cos e) from where the axis meets z = 0, the
        # point of it farthest right or left, which is where a client turning the figure puts it.
        for name, view in self.views:
            figure, where = view["figure"], f"{name} {view['id']}"
            e = math.radians(figure["camera"]["elevation"])
            for label in figure["labels"]:
                if "ring" not in label:
                    continue
                ring = label["ring"]
                circle = view["surfaces"][ring["surface"]]["rings"][ring["ring"]]
                origin = figure["turn"]["origins"][ring["surface"]]
                x = origin[0] + ring["side"] * circle["rho"]
                y = origin[1] + circle["z"] * math.cos(e)
                self.assertAlmostEqual(label["at"][1], y, delta=2e-4, msg=f"{where} {label['text']}")
                if "clear" in label:
                    self.assertGreaterEqual(ring["side"] * (label["at"][0] - x), -2e-4, f"{where} {label['text']}")
                else:
                    self.assertAlmostEqual(label["at"][0], x, delta=2e-4, msg=f"{where} {label['text']}")

    def test_at_its_own_camera_the_page_draws_the_published_figure(self):
        views = {(v["metric"], v["view"]): v for v in self.check()["views"]}
        self.assertEqual(set(views), {(name, view["id"]) for name, view in self.views})
        for (name, vid), v in views.items():
            where = f"{name} {vid}"
            self.assertLess(v["labels"], 2e-4, f"{where}: a label is not where it was published")
            self.assertTrue(v["shown"], f"{where}: a label is hidden at the start")
            for cls, (published, drawn) in v["lines"].items():
                self.assertLessEqual(abs(published - drawn), 0.02 * max(published, drawn) + 0.01 * v["size"],
                                     f"{where}: {cls} runs {drawn:.3f} where it was published {published:.3f} long")
            for cls, (published, drawn) in v["fills"].items():
                self.assertLessEqual(abs(published - drawn), 3e-3 * v["boxArea"],
                                     f"{where}: {cls} covers {drawn:.4f} where it was published covering {published:.4f}")

    def test_a_turned_figure_keeps_to_its_box_at_one_scale(self):
        for v in self.check()["views"]:
            where = f"{v['metric']} {v['view']}"
            for step in v["sweep"]:
                at = f"{where} at azimuth {step['azimuth']}, elevation {step['elevation']}"
                self.assertEqual(step["bad"], 0, f"{at}: a point is not a number")
                self.assertLessEqual(step["outside"], 1e-9, f"{at}: drawn outside the box")
                self.assertGreater(step["scale"], 0, at)
                self.assertLessEqual(step["scale"], 1, at)
            start = v["sweep"][3]
            self.assertEqual(start["scale"], 1, f"{where} is not at its own scale at its own camera")

    def test_turning_a_surface_round_its_axis_leaves_its_outline_and_circles_in_place(self):
        # A height over a plane is no surface of revolution, and turns as the next test holds it.
        for v in self.check()["views"]:
            if v["moved"] is None:
                self.assertTrue("height" in v or v.get("stack"), f"{v['metric']} {v['view']}")
                continue
            self.assertLess(v["moved"], 1e-3, f"{v['metric']} {v['view']}: {v['moved']:.2e} of the lines moved")

    def test_a_height_over_a_plane_turns_under_the_hand(self):
        # Looked at along the vertical a height hides nothing of itself and its tint covers its
        # rim at the drawn scale; from just above the plane its relief hides part of its grid;
        # and turned all the way round it keeps to its box, the Krasnikov tube's rectangle drawn
        # smaller where it would stand wider or taller than it was published.
        heights = {v["metric"]: v["height"] for v in self.check()["views"] if "height" in v}
        self.assertEqual(set(heights), set(filter(chosen, {"alcubierre", "krasnikov", "natario", "tippett_tsang",
                                                          "moving_mirror"})))
        for metric_id, h in heights.items():
            for side in ("above", "below"):
                seen = h[side]
                self.assertEqual(seen["far"], 0, f"{metric_id} from {side}: part of the height is hidden")
                self.assertGreater(seen["seen"], 0, f"{metric_id} from {side}")
                self.assertLess(abs(seen["tint"] - seen["rim"]), 0.01 * seen["rim"],
                                f"{metric_id} from {side}: the tint covers {seen['tint']:.3f} of {seen['rim']:.3f}")
            self.assertGreater(h["low"]["far"], 0.02 * (h["low"]["far"] + h["low"]["seen"]),
                               f"{metric_id} from just above the plane: the relief hides none of the grid")
            self.assertLessEqual(h["turned"]["outside"], 1e-9, f"{metric_id} turned: drawn outside the box")
        if "krasnikov" in heights:
            self.assertLess(heights["krasnikov"]["turned"]["scale"], 1)

    def test_from_straight_above_the_upper_sheet_hides_the_lower_and_from_below_the_reverse(self):
        # Schwarzschild's circles run from 1.5 to 6 r_s on each sheet. Looking down the axis every
        # circle of the lower sheet lies under the upper sheet, all but the rim at 6 r_s, which
        # lies under the upper rim's edge, and nothing hides a circle of the upper sheet.
        seen = self.check()["schwarzschild"]
        for view, near, far in (("above", "r", "r2"), ("below", "r2", "r")):
            reach = seen[view]
            self.assertNotIn(near + "-far", reach, view)
            self.assertAlmostEqual(reach[near][0], 1.5, places=6, msg=view)
            self.assertAlmostEqual(reach[far + "-far"][0], 1.5, places=6, msg=view)
            self.assertGreater(reach.get(far, [6, 6])[0], 6 - 1e-3, f"{view}: a circle under the sheet is seen")

    def test_the_drawing_turns_under_a_drag_and_leaves_a_vertical_swipe_to_the_page(self):
        page = (build.ROOT / "_layouts" / "mfs.html").read_text(encoding="utf-8")
        self.assertTrue("{{ '/MFS/assets/turn.js' | relative_url }}" in page, "the page does not load turn.js")
        rule = re.search(r"\.mfs-turn \.mfs-cd-drawing \{([^}]*)\}", page)
        self.assertIsNotNone(rule)
        self.assertIn("touch-action: pan-y", rule.group(1))
        self.assertTrue(".mfs-print-body .mfs-turn-reset, .mfs-print-body .mfs-movie-play { display: none" in page,
                        "print shows the reset or the play button")
        # Both kinds of figure that turn are framed to turn, and the one controller wires them.
        self.assertTrue(re.search(r"function pjFigure\([\s\S]*?'mfs-pj-figure mfs-turn'", page),
                        "pjFigure() does not frame its figure to turn")
        self.assertTrue(re.search(r"function emFigure\([\s\S]*?'mfs-em-figure' \+ \(turns \? ' mfs-turn'", page),
                        "emFigure() does not frame its figure to turn")
        self.assertTrue("function wireTurning(root)" in page, "no wireTurning()")
        self.assertFalse("wireEmbedding" in page, "wireEmbedding() is left")


class StacksAndMovies(unittest.TestCase):
    """The ring views of Kasner's universe, Bianchi type I and the pp-wave stacked up their axis of
    time into the ring's world tube, each ellipse at the height of its time, and the embedding
    diagrams that change through a run of moments played as movies, as the captain asked on 30
    September 2026, from the numbers written and nothing else."""

    STACKS = chosen_entries({"kasner": 1.5, "kasner_scalar": 1.5, "kasner_magnetic": 1.5, "bianchi": 2.5, "pp_wave": 0.5, "aichelburg_sexl": 0.5, "khan_penrose": 3.0,
              "bell_szekeres": 3.0, "light_beam": 0.4, "chandrasekhar_xanthopoulos": 3.0, "belinski_zakharov": 2.0,
              "big_rip": 1.0,
              "murphy_viscous_universe": 1.0,
              "maximally_supersymmetric_plane_wave": 0.8,
              "cahen_wallach": 2.0,
              "holographic_superconductor": 1.0})   # the height of a unit of time
    # Every movie, by its spacetime and view, with its variable. The last nine stood as separate
    # pictures of their moments until the captain asked on 1 October 2026 for every one of them
    # to play, and TimeSlicedViewsAreMovies keeps any other from standing so again.
    MOVIES = chosen_entries({("chodos_detweiler", "cylinder"): "$t/t_0$",
              ("frw", "closed"): "$ct$", ("nordstrom_scalar", "dust"): "$ct$",
              ("kopczynski_trautman", "universe"): "$ct$", ("small_universes", "torus"): "$ct$",
              ("lqc_bounce", "universe"): "$t$",
              ("bekenstein_universes", "closed"): "$\\eta$", ("bekenstein_universes", "open"): "$\\eta$",
              ("malament_hogarth", "plane"): "$ct$", ("mixmaster", "sphere"): "$c\\tau$",
              ("oppenheimer_snyder", "collapse"): "$c\\tau$", ("white_hole", "explosion"): "$c\\tau$",
              ("vaidya", "shell"): "$v - r$",
              ("semiclosed_world", "bag"): "$c\\tau$", ("datt_ruban_t_models", "tsphere"): "$c\\tau$",
              ("lindquist_wheeler_lattice", "lattice"): "$c\\tau$",
              ("bonnor_vaidya", "shell"): "$v - r$", ("israel_shell", "shell"): "$v - r$", ("black_to_white_hole", "collapse"): "$ct$",
              ("charged_shell", "bounce"): "$v - r$",
              ("penrose_impulsive_wave", "snap"): "$ct$",
              ("ab_metrics", "neck"): "$\\tau$",
              ("cosmic_string", "unroll"): "$\\Delta\\phi$", ("point_particle_2plus1", "unroll"): "$\\Delta\\phi$",
              ("milne", "hyperboloids"): "$ct$",
              ("coleman_de_luccia", "hyperboloids"): "$c\\tau$",
              ("einstein_rosen_waves", "pulse"): "$ct$", ("nariai", "universe"): "$ct$",
              ("twisting_cylindrical_waves", "pulse"): "$ct$",
              ("domain_wall", "moments"): "$kct$", ("kantowski_sachs", "dust"): "$\\eta$",
              ("robinson_trautman", "fronts"): "$cu$", ("mcvittie", "flamm"): "$ct$",
              ("supertranslation_hair", "shells"): "$\\rho$",
              ("sultana_dyer", "paraboloid"): "$\\eta$",
              ("kastor_traschen", "two_holes"): "$c\\tau$",
              ("misner_brill_lindquist", "through"): "$a$", ("misner_brill_lindquist", "between"): "$a$",
              ("brill_waves", "strong"): "$a$",
              ("tolman_bondi", "cloud"): "$ct$", ("szekeres", "equators"): "$ct$", ("misner", "cylinders"): "$ct$",
              ("self_creating_universe", "cylinders"): "$c\\tau$",
              ("photon_rocket", "burn"): "$cu + r$", ("hayward", "history"): "$v - r$",
              ("mass_inflation", "tail"): "$v - r$",
              ("hiscock", "history"): "$v - r$",
              ("gott_time_machine", "cylinders"): "$c\\tau$", ("kantowski_sachs", "vacuum"): "$c\\tau$",
              ("ori_time_machine", "throat"): "$t$", ("senovilla", "universe"): "$act$",
              ("kasner", "ring"): "$t$", ("kasner_scalar", "ring"): "$t$", ("kasner_magnetic", "ring"): "$t$", ("bianchi", "ring"): "$c\\bar Ht$", ("pp_wave", "ring"): "$cu$",
              ("aichelburg_sexl", "ring"): "$u$", ("khan_penrose", "ring"): "$\\tau$",
              ("bell_szekeres", "ring"): "$\\xi$", ("chandrasekhar_xanthopoulos", "ring"): "$\\psi$",
              ("belinski_zakharov", "ring"): "$\\tau$", ("big_rip", "ring"): "$t$", ("murphy_viscous_universe", "ring"): "$t$",
              ("gowdy", "torus"): "$t$", ("light_beam", "ring"): "$u$",
              ("maximally_supersymmetric_plane_wave", "ring"): "$\\mu cu$",
              ("cahen_wallach", "ring"): "$cu/L$",
              ("string_wave", "ring"): "$u$", ("simpson_visser", "inside"): "$c\\tau$",
              ("hotta_tanaka", "ring"): "$\\tau$", ("holographic_superconductor", "ring"): "$-\\rho$",
              **{("roberts", case): "$ct$" for case in ("disperses", "threshold", "collapses")},
              ("black_string", "ripple"): "$v$"}, lambda key: key[0])

    def setUp(self):
        self.embedding = embedding_files()
        self.named = NamedFiles(build.EMBEDDING_DIR)

    def views(self, metric_id):
        return {v["id"]: v for v in self.named[metric_id]["views"]}

    def movies(self):
        """Every view with a movie, by its spacetime and view."""
        return {(name, v["id"]): v for name, data in self.embedding.items() for v in data["views"] if "movie" in v}

    def test_the_stacks_come_first_with_their_flat_rings_beside_them(self):
        for metric_id in self.STACKS:
            self.assertEqual([v["id"] for v in self.named[metric_id]["views"]], ["tube", "ring"], metric_id)
        stacked = {name for name, data in self.embedding.items() for v in data["views"]
                   for s in v["surfaces"] for p in s["pieces"] if p.get("grid", {}).get("frame") == "ellipses"}
        self.assertEqual(stacked, set(self.STACKS))

    def test_every_ellipse_of_a_stack_stands_at_the_height_of_its_time(self):
        for metric_id, lift in self.STACKS.items():
            views = self.views(metric_id)
            tube, = views["tube"]["surfaces"]
            grid = tube["pieces"][0]["grid"]
            rows = grid_nodes(tube["pieces"][0])
            # Every row is one moment, at the height of its time.
            for u, z in zip(grid["u"], grid["z"]):
                self.assertAlmostEqual(z, lift * u, delta=1e-6, msg=f"{metric_id} row at {u}")
            rings = [c for c in tube["curves"] if "time" in c]
            flat = views["ring"]["surfaces"]
            self.assertEqual([(c["label"], c["time"]) for c in rings], [(f["label"], f["time"]) for f in flat], metric_id)
            for ring, moment in zip(rings, flat):
                where = f"{metric_id} {ring['label']}"
                # The time is written to six decimals and the row at the double it was drawn at.
                i = min(range(len(grid["u"])), key=lambda k: abs(grid["u"][k] - ring["time"]))
                self.assertLess(abs(grid["u"][i] - ring["time"]), 1e-6, where)
                self.assertTrue(all(Z == grid["z"][i] for _, _, Z in ring["points"]), where)
                self.assertEqual(len(ring["points"]), len(grid["v"]), where)
                for P, Q in zip(ring["points"], rows[i]):
                    self.assertLess(math.dist(P, Q), 1e-6, where)
                # The same ellipse as the flat view's ring at that moment, only raised to its time.
                for P, Q in zip(ring["points"], moment["curves"][0]["points"]):
                    self.assertLess(math.hypot(P[0] - Q[0], P[1] - Q[1]), 1e-6, where)
            dots = tube["dots"]
            self.assertEqual(len(dots), 12 * len(rings), metric_id)
            axis = tube["axis"]
            self.assertEqual(axis["from"], grid["z"][0], metric_id)
            self.assertAlmostEqual(axis["to"], grid["z"][-1] + (grid["z"][-1] - grid["z"][0]) / 4, delta=1e-6, msg=metric_id)
            # The world lines of the twelve marked particles, every 30 degrees round the ring.
            self.assertEqual(grid["lines"], [{"class": "worldline", "u": [], "v": list(range(0, 360, 30))}], metric_id)
        # Kasner's rows reach t^p_1 along x and t^p_3 along z; the wave's close on a segment at the focus.
        if chosen("kasner"):
            kasner = self.views("kasner")["tube"]["surfaces"][0]["pieces"][0]["grid"]
            for u, a, b in zip(kasner["u"], kasner["a"], kasner["b"]):
                self.assertAlmostEqual(a, u ** (-2 / 7), delta=1e-6)
                self.assertAlmostEqual(b, u ** (6 / 7), delta=1e-6)
        if chosen("pp_wave"):
            wave = self.views("pp_wave")["tube"]["surfaces"][0]["pieces"][0]["grid"]
            self.assertEqual(wave["b"][-1], 0)

    def test_every_moment_label_of_a_stack_follows_its_ring_and_the_axis_is_named(self):
        for metric_id in self.STACKS:
            view = self.views(metric_id)["tube"]
            labels = view["figure"]["labels"]
            named = [L for L in labels if "curve" in L]
            self.assertEqual([L["curve"]["curve"] for L in named], list(range(len(named))), metric_id)
            self.assertEqual([L["text"] for L in named], [c["label"] for c in view["surfaces"][0]["curves"]], metric_id)
            # Right and left in turn from the first up, for any number of moments.
            self.assertEqual([L["curve"]["side"] for L in named], [1 if i % 2 == 0 else -1 for i in range(len(named))], metric_id)
            axis, = [L for L in labels if "axis" in L]
            self.assertEqual(axis["axis"], {"surface": 0}, metric_id)
            self.assertIn(("line", "axis"), [(k, c) for k, c, _ in view["figure"]["legend"]], metric_id)

    def test_the_movies_are_these(self):
        self.assertEqual({key: v["movie"]["variable"] for key, v in self.movies().items()}, self.MOVIES)
        self.assertEqual({name for (name, _), v in self.movies().items() if v["movie"].get("turns") is False},
                         set(filter(chosen, {"frw", "milne", "coleman_de_luccia"})))

    def test_every_movie_runs_through_its_frames_in_order_and_holds_its_moments(self):
        for (metric_id, view_id), view in self.movies().items():
            movie, where = view["movie"], f"{metric_id} {view_id}"
            values = [f["value"] for f in movie["frames"]]
            self.assertGreater(len(values), 30, where)
            self.assertTrue(all(b > a for a, b in zip(values, values[1:])), where)
            self.assertGreater(movie["seconds"], 0, where)
            self.assertEqual(len({f["label"] for f in movie["frames"]}), len(values), f"{where}: two frames share a name")
            frame_label, = [L for L in view["figure"]["labels"] if L.get("frame")]
            self.assertEqual(frame_label["text"], movie["frames"][0]["label"], where)
            if metric_id in ("cosmic_string", "point_particle_2plus1"):
                # A cone unrolling is one moment, moved without stretching.
                continue
            # A movie in time passes through every moment of its flat views, which are its frames.
            for surface in view["surfaces"]:
                frame = next(f for f in movie["frames"] if f["value"] == surface["time"])
                self.assertEqual(frame["label"], surface["label"], where)
                self.assertEqual(frame["pieces"], surface["pieces"], where)
                self.assertEqual(frame["rings"], surface["rings"], where)
                self.assertEqual(frame.get("curves"), surface.get("curves"), where)
                self.assertEqual(frame.get("dots"), surface.get("dots"), where)
            self.assertEqual((values[0], values[-1]), (view["surfaces"][0]["time"], view["surfaces"][-1]["time"]), where)

    @checks("frw")
    def test_the_closed_universe_plays_every_frame_as_its_sphere(self):
        # At ct = eta - sin(eta) the equator is the sphere of radius a = 1 - cos(eta).
        view = self.views("frw")["closed"]
        for frame in view["movie"]["frames"]:
            eta = bisect(lambda e: e - math.sin(e) - frame["value"], 0, 2 * math.pi)
            a = 1 - math.cos(eta)
            for piece, sign in (("near", -1), ("far", 1)):
                points = next(p for p in frame["pieces"] if p["id"] == piece)["points"]
                for r, rho, z in points:
                    self.assertLess(abs(rho - a * r), 1e-5, f"ct = {frame['value']}")
                    self.assertLess(abs(z - sign * a * math.sqrt(max(1 - r * r, 0))), 1e-5, f"ct = {frame['value']}")

    @checks("malament_hogarth", "oppenheimer_snyder", "vaidya", "israel_shell", "charged_shell", "tolman_bondi")
    def test_each_movie_with_a_rim_holds_it_still_while_the_rest_moves(self):
        # The flat plane beyond the Malament-Hogarth ball, the clocks released at 4 r_s outside
        # the collapsing dust, Vaidya's r = 4 r_s and the clocks released at 2.5 r_b outside
        # Tolman-Bondi's cloud stand at z = 0 in every frame, and the
        # Malament-Hogarth well only deepens as the removed event nears.
        for metric_id, piece in (("malament_hogarth", "flat"), ("oppenheimer_snyder", "exterior"), ("vaidya", None),
                                 ("israel_shell", None), ("charged_shell", None), ("tolman_bondi", "exterior")):
            view = self.named[metric_id]["views"][0]
            for frame in view["movie"]["frames"]:
                outer = frame["pieces"][-1] if piece is None else next(p for p in frame["pieces"] if p["id"] == piece)
                self.assertEqual(outer["points"][-1][2], 0, f"{metric_id} at {frame['value']}")
        depths = [-f["pieces"][0]["points"][0][2] for f in self.named["malament_hogarth"]["views"][0]["movie"]["frames"]]
        self.assertTrue(all(b >= a for a, b in zip(depths, depths[1:])), "the Malament-Hogarth well rises")

    @checks("cosmic_string")
    def test_the_cone_unrolls_without_stretching_to_the_deficit_angle(self):
        view = self.views("cosmic_string")["unroll"]
        frames = view["movie"]["frames"]
        delta = 8 * math.pi * (1 / 40)          # 8 pi G mu/c^2 at 4G mu/c^2 = 0.1
        for frame in frames:
            piece, = frame["pieces"]
            grid, where = piece["grid"], f"the cone at {frame['label']}"
            self.assertTrue(grid["open"], where)
            nodes = grid_nodes(piece)
            gap = 2 * math.pi - (grid["v"][-1] - grid["v"][0])
            self.assertAlmostEqual(gap, math.radians(frame["value"]), delta=1e-12, msg=where)
            # Along each line from the apex every distance is r, and round each circle every arc
            # is as long as on the cone, 0.9 r for each radian of phi.
            for i in range(len(grid["u"]) - 1):
                for j in range(len(grid["v"])):
                    self.assertAlmostEqual(math.dist(nodes[i][j], nodes[i + 1][j]), grid["u"][i + 1] - grid["u"][i],
                                           delta=1e-6, msg=where)
            step = 2 * math.pi / (len(grid["v"]) - 1)
            for u, a, b in zip(grid["u"], grid["a"], grid["b"]):
                self.assertEqual(a, b, where)
                for dv in (y - x for x, y in zip(grid["v"], grid["v"][1:])):
                    self.assertAlmostEqual(a * dv, 0.9 * u * step, delta=1e-6, msg=where)
        first, last = frames[0]["pieces"][0]["grid"], frames[-1]["pieces"][0]["grid"]
        self.assertEqual((first["v"][0], first["v"][-1]), (0, 2 * math.pi))
        self.assertAlmostEqual(first["a"][-1], 0.9 * first["u"][-1], delta=1e-6)
        # Laid flat, the cone is a plane missing the wedge of the deficit angle between its cut's edges.
        self.assertTrue(all(z == 0 for z in last["z"]))
        self.assertAlmostEqual(2 * math.pi - (last["v"][-1] - last["v"][0]), delta, delta=1e-12)
        self.assertAlmostEqual(math.degrees(delta), 36, delta=1e-12)
        self.assertEqual(frames[-1]["label"], "$\\Delta\\phi = 36° = \\delta$")

    def test_no_movie_says_how_it_loops(self):
        # Forward and back is the one way a movie plays, so no movie's data, the build and the
        # player carry a word for another.
        for key, view in self.movies().items():
            self.assertLessEqual(set(view["movie"]), {"variable", "seconds", "frames", "turns"}, str(key))
        view = copy.deepcopy(next(v for v in self.named["cosmic_string"]["views"] if "movie" in v))
        build.check_movie("cosmic_string", view)
        for loop in ("once", "pingpong", "restart"):
            view["movie"]["loop"] = loop
            with self.assertRaises(build.DataError) as raised:
                build.check_movie("cosmic_string", view)
            self.assertIn("every movie plays forward and back", str(raised.exception))
        player = (build.ROOT / "MFS" / "assets" / "turn.js").read_text(encoding="utf-8")
        body = re.search(r"\n  function movieFrame\(movie, ms\) \{[\s\S]*?\n  \}\n", player)
        self.assertIsNotNone(body, "no movieFrame()")
        self.assertEqual(set(re.findall(r"\bmovie\.(\w+)", body.group(0))), {"frames", "seconds"})
        page = (build.ROOT / "_layouts" / "mfs.html").read_text(encoding="utf-8")
        self.assertEqual(page.count("MfsTurn.movieFrame(movie, elapsed)"), 1, "the page plays a movie some other way")
        self.assertNotIn("movie.loop", page)

    def frames_shown(self, movies, times):
        """The frame MfsTurn.movieFrame() shows of each movie at each of `times` in milliseconds,
        run in Node as the page runs it."""
        if shutil.which("node") is None:
            self.skipTest("Node is not installed, so the page's movie player cannot be run")
        script = ("const t = require(process.argv[1]); const [movies, times] = JSON.parse(require('fs').readFileSync(0));"
                  "process.stdout.write(JSON.stringify(movies.map(m => times.map(ms => t.movieFrame(m, ms)))));")
        movies = [{k: m[k] for k in ("seconds", "loop") if k in m} | {"frames": [{"value": f["value"]} for f in m["frames"]]}
                  for m in movies]
        run = subprocess.run(["node", "-e", script, str(build.ROOT / "MFS" / "assets" / "turn.js")],
                             input=json.dumps([movies, times]), capture_output=True, text=True, timeout=120)
        self.assertEqual(run.returncode, 0, run.stderr[-2000:])
        return json.loads(run.stdout)

    def test_the_cosmic_string_turns_round_at_each_end_without_a_stutter(self):
        movie = next(v for v in self.named["cosmic_string"]["views"] if "movie" in v)["movie"]
        n, step = len(movie["frames"]), 1000 * movie["seconds"] / (len(movie["frames"]) - 1)
        # At the middle of each step: 0, 1, ..., n - 2, n - 1, n - 2, ..., 1, 0, 1, ..., each end once.
        [shown] = self.frames_shown([movie], [(j + 0.5) * step for j in range(4 * (n - 1) + 1)])
        up, down = list(range(n)), list(range(n - 2, 0, -1))
        self.assertEqual(shown, up + down + up + down + [0])
        # Millisecond by millisecond over two cycles, every frame, the ends and the first included,
        # stays on the screen for one step, so the turnarounds keep the pace.
        total = int(4 * (n - 1) * step)
        [shown] = self.frames_shown([movie], list(range(total)))
        runs = [len(list(g)) for _, g in itertools.groupby(shown)]
        self.assertEqual(len(runs), 4 * (n - 1))
        for run in runs[:-1]:
            self.assertLessEqual(abs(run - step), 1, runs)
        self.assertTrue(all(abs(a - b) == 1 for a, b in zip(shown, shown[1:]) if a != b))

    def test_every_movie_plays_forward_then_back_and_never_jumps_to_its_start(self):
        movies = {key: v["movie"] for key, v in self.movies().items()}
        self.assertEqual(set(movies), set(self.MOVIES))
        found(self, movies, "movie")
        # Three passes and a little more, forward, back and forward again, every 7 milliseconds.
        times = list(range(0, 3 * 1000 * max(m["seconds"] for m in movies.values()) + 2000, 7))
        for (name, movie), shown in zip(movies.items(), self.frames_shown(list(movies.values()), times)):
            values, period = [f["value"] for f in movie["frames"]], 1000 * movie["seconds"]
            n, span = len(values), values[-1] - values[0]
            half = period * (values[1] - values[0]) / span / 2
            # The frame shown is the one nearest in value, to within rounding where two are as near.
            for ms, k in zip(times, shown):
                p = (ms - half) % (2 * period)
                v = values[0] + span * (p if p < period else 2 * period - p) / period
                self.assertLessEqual(abs(values[k] - v), min(abs(x - v) for x in values) + 1e-6 * span, (name, ms))
            # It starts on its first frame and moves one frame at a time, so it never jumps, least
            # of all from its last frame to its first.
            self.assertEqual(shown[0], 0, str(name))
            self.assertTrue(all(abs(a - b) <= 1 for a, b in zip(shown, shown[1:])), str(name))
            # Its frames in the order shown: up to the last, down to the first, and up again.
            order = [k for k, _ in itertools.groupby(shown)]
            up, down = list(range(n)), list(range(n - 2, 0, -1))
            self.assertEqual(order, (up + down + up + down)[:len(order)], str(name))
            self.assertGreaterEqual(len(order), 3 * (n - 1), str(name))

    def frames_played(self, movies, ms, pause=None):
        """The frames the page's own wireMovie() hands its figure of each movie, in order, over
        `ms` milliseconds at sixty pictures a second, run in Node as the page runs it, with a
        clock and a play button stood in for the browser's. With `pause` as [at, for], the reader
        presses the button `at` milliseconds in and again `for` milliseconds later."""
        if shutil.which("node") is None:
            self.skipTest("Node is not installed, so the page's movie player cannot be run")
        page = (build.ROOT / "_layouts" / "mfs.html").read_text(encoding="utf-8")
        player = re.search(r"\n      var _movieObserver = null, _movieSeen = new Map\(\);\n"
                           r"      function wireMovie\([\s\S]*?\n      \}\n", page)
        self.assertIsNotNone(player, "no wireMovie()")
        script = ("const MfsTurn = require(process.argv[1]);"
                  "const [source, movies, ms, pause] = JSON.parse(require('fs').readFileSync(0));"
                  "let queue = [];"
                  "const window = {}, requestAnimationFrame = f => { queue.push(f); return queue.length; };"
                  "const wireMovie = new Function('MfsTurn', 'window', 'requestAnimationFrame',"
                  " source + '; return wireMovie;')(MfsTurn, window, requestAnimationFrame);"
                  "process.stdout.write(JSON.stringify(movies.map(movie => {"
                  "  queue = [];"
                  "  const played = [], play = { addEventListener(_, f) { this.click = f; } };"
                  "  const drawing = { isConnected: true, closest: () => ({ querySelector: () => play }) };"
                  "  wireMovie(drawing, movie, k => played.push(k), () => {});"
                  "  let pressed = 0;"
                  "  for (let now = 1000; now <= 1000 + ms; now += 1000 / 60) {"
                  "    if (pause && pressed < 2 && now - 1000 >= pause[0] + pressed * pause[1]) { play.click(); pressed++; }"
                  "    const due = queue; queue = []; due.forEach(f => f(now));"
                  "  }"
                  "  return played;"
                  "})));")
        movies = [{"seconds": m["seconds"], "frames": [{"value": f["value"]} for f in m["frames"]]} for m in movies]
        run = subprocess.run(["node", "-e", script, str(build.ROOT / "MFS" / "assets" / "turn.js")],
                             input=json.dumps([player.group(0), movies, ms, pause]), capture_output=True, text=True,
                             timeout=120)
        self.assertEqual(run.returncode, 0, run.stderr[-2000:])
        return json.loads(run.stdout)

    def test_the_page_itself_plays_every_movie_forward_then_back(self):
        # The path a reader's movie takes, wireMovie() in the page asking MfsTurn.movieFrame() at
        # every picture, for every movie the files hold and not only the frames the rule names: a
        # movie that played once and started again would hand its figure its first frame straight
        # after its last.
        movies = {key: v["movie"] for key, v in self.movies().items()}
        self.assertEqual(set(movies), set(self.MOVIES))
        found(self, movies, "movie")
        ms = 3 * 1000 * max(m["seconds"] for m in movies.values()) + 2000
        stood = 4321
        straight = self.frames_played(list(movies.values()), ms)
        shorter = self.frames_played(list(movies.values()), ms - stood)
        paused = self.frames_played(list(movies.values()), ms, [3210, stood])
        for (name, movie), played, short, held in zip(movies.items(), straight, shorter, paused):
            n = len(movie["frames"])
            up, down = list(range(n)), list(range(n - 2, 0, -1))
            # It starts on its first frame, so the first it is handed is the second, and in three
            # passes it is handed every frame on the way up, down and up again.
            self.assertEqual(played, (up + down + up + down)[1:len(played) + 1], str(name))
            self.assertGreaterEqual(len(played), 3 * (n - 1), str(name))
            # Paused, it takes up where it stopped, so it loses the time it stood and no more.
            self.assertEqual(held, played[:len(held)], str(name))
            self.assertLessEqual(abs(len(held) - len(short)), 1, str(name))
        # Every figure with a movie is wired to that one player, an embedding diagram's by
        # wireTurning() and a curvature figure's by wireCurvature(), and nothing else in the page
        # or the player keeps time for a movie.
        page = (build.ROOT / "_layouts" / "mfs.html").read_text(encoding="utf-8")
        self.assertEqual(re.findall(r"if \(movie\) (\w+)\(", page), ["wireMovie"])
        self.assertEqual(len(re.findall(r"\bwireMovie\(", page)), 3, "the page plays a movie some other way")
        curvature = re.search(r"function wireCurvature\([\s\S]*?\n      \}\n", page)
        self.assertIsNotNone(curvature, "no wireCurvature()")
        self.assertEqual(len(re.findall(r"\bwireMovie\(", curvature.group(0))), 1)
        player = (build.ROOT / "MFS" / "assets" / "turn.js").read_text(encoding="utf-8")
        for clock in ("requestAnimationFrame", "setInterval", "setTimeout", "performance.now", "Date.now"):
            self.assertNotIn(clock, player, "MfsTurn keeps time of its own")

    def test_the_page_plays_every_curvature_movie_forward_then_back(self):
        # A curvature figure's movie is played by the same player, by the same rule.
        movies = {(name, v["id"]): v["movie"] for name, data in curvature_files().items()
                  for v in data["views"] if "movie" in v}
        found(self, movies, "curvature figure holding a movie")
        ms = 3 * 1000 * max(m["seconds"] for m in movies.values()) + 2000
        for (name, movie), played in zip(movies.items(), self.frames_played(list(movies.values()), ms)):
            n = len(movie["frames"])
            up, down = list(range(n)), list(range(n - 2, 0, -1))
            # The run lasts three passes of the longest movie, so a shorter one goes up and down as
            # many times over as fit, every one of them the same cycle.
            cycles = len(played) // len(up + down) + 2
            self.assertEqual(played, ((up + down) * cycles)[1:len(played) + 1], str(name))
            self.assertGreaterEqual(len(played), 3 * (n - 1), str(name))

    def test_the_page_plays_a_movie_and_holds_it_still_for_a_reader_who_asks(self):
        page = (build.ROOT / "_layouts" / "mfs.html").read_text(encoding="utf-8")
        player = re.search(r"function wireMovie\([\s\S]*?\n      \}\n", page)
        self.assertIsNotNone(player, "no wireMovie()")
        self.assertIn("prefers-reduced-motion: reduce", player.group(0))
        self.assertIn("IntersectionObserver", player.group(0))
        self.assertRegex(page, r"\.mfs-movie-play:active \{ color: var\(--pink-light\); border-color: var\(--pink-light\); \}")


class SteadyMovieLabels(unittest.TestCase):
    """Nothing in the label that names a movie's frame moves as the movie plays, as the captain
    asked on 1 October 2026: until then the label was centred on its place, so the whole of it
    shifted whenever its number gained or lost a digit. Every movie the files hold is held to
    it, found by its `movie` and by the figure's label marked `frame`, and none by name: the
    label's place is the same at every frame, its names start at one left edge in a box as wide
    as the widest, and MfsTurn.steady() gives every frame's number one width with its point at
    one place, counting on nothing of a font but that its digits are of one width."""

    PADDED = re.compile(r"^(?P<lead>[^=]*=\s*)(?:\\hphantom\{(?P<front>(?:[0./]|\{-\})+)\})?(?P<sign>\{-\})?(?P<whole>\d+)"
                        r"(?P<rest>(?:\.\d+)?(?:/\d+)?)(?:\\hphantom\{(?P<back>(?:[0./]|\{-\})+)\})?(?P<tail>.*)$", re.S)

    def setUp(self):
        self.movies = {(name, v["id"]): v for name, data in embedding_files().items()
                       for v in data["views"] if "movie" in v}

    def steady(self, lists):
        """MfsTurn.steady() of each list of labels, run in Node as the page runs it."""
        if shutil.which("node") is None:
            self.skipTest("Node is not installed, so the page's labels cannot be made")
        script = ("const t = require(process.argv[1]); const lists = JSON.parse(require('fs').readFileSync(0));"
                  "process.stdout.write(JSON.stringify(lists.map(t.steady)));")
        run = subprocess.run(["node", "-e", script, str(build.ROOT / "MFS" / "assets" / "turn.js")],
                             input=json.dumps(lists), capture_output=True, text=True, timeout=120)
        self.assertEqual(run.returncode, 0, run.stderr[-2000:])
        return json.loads(run.stdout)

    @staticmethod
    def glyphs(*texts):
        """What a stretch of a number is made of, phantom or shown: how many of each kind of
        glyph, every digit one kind, which is its width in any font of digits of one width."""
        return collections.Counter(re.sub(r"\d", "0", "".join(t or "" for t in texts).replace("{-}", "-")))

    def test_every_frames_number_fills_one_box_with_its_point_at_one_place(self):
        labels = {key: [f["label"] for f in v["movie"]["frames"]] for key, v in self.movies.items()}
        # A curvature figure's movie names its frames by the same rule.
        labels.update({("curvature", name, v["id"]): [f["label"] for f in v["movie"]["frames"]]
                       for name, data in curvature_files().items() for v in data["views"] if "movie" in v})
        found(self, self.movies, "movie")
        for (key, plain), made in zip(labels.items(), self.steady(list(labels.values()))):
            self.assertEqual(len(made), len(plain), key)
            boxes, shown = set(), [collections.Counter(), collections.Counter()]
            for was, now in zip(plain, made):
                where = f"{key}: {was!r} as {now!r}"
                # The whole label is one formula, so the phantoms are set in its mathematics.
                self.assertRegex(was, r"^\$[^$]+\$$", where)
                # Nothing that shows is changed.
                self.assertEqual(re.sub(r"\\hphantom\{(?:[0./]|\{-\})+\}", "", now).replace("{-}", "-"), was, where)
                m = self.PADDED.match(now)
                self.assertIsNotNone(m, where)
                front = self.glyphs(m["front"], m["sign"], m["whole"])
                back = self.glyphs(m["rest"], m["back"])
                self.assertLessEqual(set(front), {"-", "0"}, where)
                self.assertLessEqual(set(back), {".", "/", "0"}, where)
                boxes.add((m["lead"], tuple(sorted(front.items())), tuple(sorted(back.items()))))
                shown[0] |= self.glyphs(m["sign"], m["whole"])
                shown[1] |= self.glyphs(m["rest"])
            # What stands before the number, the room before its point and the room from its
            # point on are the same in every frame.
            self.assertEqual(len(boxes), 1, f"{key}: the number's box changes from frame to frame, {sorted(boxes)}")
            # No room is kept that no frame fills: of each kind of glyph, the most any frame shows.
            [(_, front, back)] = boxes
            self.assertEqual([dict(front), dict(back)], [dict(c) for c in shown], key)

    def test_a_number_is_padded_with_what_it_lacks_and_nothing_else_is_touched(self):
        cases = [
            # A sign and a second digit: each frame gets the one it lacks.
            (["$cu + r = -2\\,m_0$", "$cu + r = 0.5\\,m_0$", "$cu + r = 24\\,m_0$"],
             ["$cu + r = \\hphantom{0}{-}2\\hphantom{.0}\\,m_0$", "$cu + r = \\hphantom{0{-}}0.5\\,m_0$",
              "$cu + r = \\hphantom{{-}}24\\hphantom{.0}\\,m_0$"]),
            # Decimals of different lengths, and what follows only the last frame's number.
            (["$\\Delta\\phi = 0°$", "$\\Delta\\phi = 36° = \\delta$"],
             ["$\\Delta\\phi = \\hphantom{0}0°$", "$\\Delta\\phi = 36° = \\delta$"]),
            (["$t = 1/4$", "$t = 0.35$", "$t = 2$"],
             ["$t = 1/4\\hphantom{.0}$", "$t = 0.35\\hphantom{/}$", "$t = 2\\hphantom{.00/}$"]),
            # Labels already of one shape are left as they are, a minus braced.
            (["$kct = -1.00$", "$kct = 0.00$"], ["$kct = {-}1.00$", "$kct = \\hphantom{{-}}0.00$"]),
            (["$\\eta = 0.00$", "$\\eta = 1.20$"], ["$\\eta = 0.00$", "$\\eta = 1.20$"]),
            # A digit in what stands before the "=" or after the number is no part of the value.
            (["$r_2 = 9\\,m_0$", "$r_2 = 10\\,m_0$"], ["$r_2 = \\hphantom{0}9\\,m_0$", "$r_2 = 10\\,m_0$"]),
            # A label with no number after an "=", or with its number outside the mathematics,
            # where a phantom would be printed as it is written, leaves every label as it came.
            (["$t = 1$", "the end"], ["$t = 1$", "the end"]),
            (["t = 1", "t = 10"], ["t = 1", "t = 10"]),
            (["$t$ = 1", "$t$ = 10"], ["$t$ = 1", "$t$ = 10"]),
            ([], []),
        ]
        self.assertEqual(self.steady([plain for plain, _ in cases]), [made for _, made in cases])

    def test_the_label_stands_at_one_place_whatever_the_frame(self):
        # One label names the frame, and the figure places it, never a frame.
        found(self, self.movies, "movie")
        for key, view in self.movies.items():
            named = [L for L in view["figure"]["labels"] if L.get("frame")]
            self.assertEqual(len(named), 1, key)
            self.assertEqual(named[0]["text"], view["movie"]["frames"][0]["label"], key)
            for frame in view["movie"]["frames"]:
                self.assertFalse({"at", "anchor", "dx", "dy", "labels"} & set(frame), key)
        # Turned, the page draws each frame's labels again: at the figure's own camera and four
        # others, every frame puts the label where the first does.
        checked = {(v["metric"], v["view"]): v["movie"] for v in turn_check(self)["views"] if "movie" in v}
        self.assertEqual(set(checked), set(self.movies))
        for key, movie in checked.items():
            self.assertEqual(movie["named"], 1, key)
            self.assertEqual(movie["label"], 0, f"{key}: the label moves {movie['label']:.2e} between frames")

    def test_the_page_sets_every_frames_name_from_one_left_edge_in_one_box(self):
        page = (build.ROOT / "_layouts" / "mfs.html").read_text(encoding="utf-8")
        figure = re.search(r"function emFigure\([\s\S]*?\n      \}\n", page)
        self.assertIsNotNone(figure, "no emFigure()")
        self.assertIn("MfsTurn.steady(names)", figure.group(0))
        self.assertIn("'<span class=\"mfs-frame-labels\">'", figure.group(0))
        curvature = re.search(r"function cvFigure\([\s\S]*?\n      \}\n", page)
        self.assertIsNotNone(curvature, "no cvFigure()")
        self.assertIn("MfsTurn.steady(names)", curvature.group(0))
        self.assertIn('<span class="cv-word mfs-frame-labels">', curvature.group(0))
        # Every name in one cell from its left edge, and a name not shown keeps its room, so the
        # label the place pins is as wide at every frame.
        rules = dict(re.findall(r"#mfs-content-panel (\.mfs-frame-label[^ {]*) \{([^}]*)\}", page))
        self.assertEqual(set(rules), {".mfs-frame-labels", ".mfs-frame-label", ".mfs-frame-label[hidden]"})
        for want in ("display: grid;", "justify-items: start;", "font-variant-numeric: tabular-nums;"):
            self.assertIn(want, rules[".mfs-frame-labels"])
        self.assertIn("grid-area: 1 / 1;", rules[".mfs-frame-label"])
        self.assertIn("visibility: hidden;", rules[".mfs-frame-label[hidden]"])
        self.assertNotIn("none", rules[".mfs-frame-label[hidden]"])
        # Playing changes which name shows and nothing else of the label.
        self.assertEqual(re.findall(r"names\.forEach\(function\(el, k\) \{([^}]*)\}", page), [" el.hidden = k !== shown; "])

    def test_the_value_under_the_pointer_keeps_its_room(self):
        page = (build.ROOT / "_layouts" / "mfs.html").read_text(encoding="utf-8")
        self.assertIn(".padStart(Math.max(lo.toFixed(3).length, hi.toFixed(3).length))", page)
        self.assertEqual(len(re.findall(r"r[xy]\.textContent = read\(", page)), 2)
        # A curvature figure's values are set by the two readers that keep their room, its
        # coordinates as a spacetime diagram's are and its values in ten places.
        setters = re.findall(r"\br[xyv]\.textContent = (\w+)\(", page)
        self.assertEqual(len(re.findall(r"\br[xyv]\.textContent = ", page)), len(setters))
        self.assertEqual(set(setters), {"read", "cvRead", "cvReadValue"})
        self.assertIn("return at.toFixed(3).replace('-', '−').padStart(Math.max(lo.toFixed(3).length, hi.toFixed(3).length));",
                      page[page.index("function cvRead("):])
        self.assertIn(".padStart(10);", page[page.index("function cvReadValue("):])
        self.assertRegex(page, r"\.mfs-nr-readout \.nr-ry \{ white-space: pre; font-variant-numeric: tabular-nums; \}")


class TimeSlicedViewsAreMovies(unittest.TestCase):
    """No embedding diagram sets the moments of a time out as separate pictures, as the captain
    asked on 1 October 2026, when Tolman-Bondi's cloud still stood as four of them beside the
    movies: a view that changes with a time, t, tau, eta, u or any other, plays as a movie with
    the page's play and pause button and a label naming the frame shown. The script that draws
    them, the build that reads them and the published files are each held to it."""

    def setUp(self):
        self.embedding = embedding_files()
        self.named = NamedFiles(build.EMBEDDING_DIR)
        self.views = [(name, view) for name, data in self.embedding.items() for view in data["views"]]

    def test_every_view_of_more_than_one_moment_is_a_movie(self):
        sliced = {(name, view["id"]) for name, view in self.views
                  if len(view["surfaces"]) > 1 or any("time" in s for s in view["surfaces"])}
        found(self, sliced, "view of more than one moment")
        self.assertEqual({key for key in sliced if key not in StacksAndMovies.MOVIES}, set())
        for name, view in self.views:
            if (name, view["id"]) in sliced:
                self.assertIn("movie", view, f"{name} {view['id']}: moments set out as separate pictures")

    def test_no_figure_stands_two_surfaces_side_by_side(self):
        # A movie draws one frame at a time on one axis, so a figure of more than one place holds
        # moments as separate pictures, and so does one with a label for each of them.
        for name, view in self.views:
            figure, where = view["figure"], f"{name} {view['id']}"
            self.assertEqual(len(figure["turn"]["origins"]), 1, where)
            moments = {s["label"] for s in view["surfaces"] if "label" in s}
            named = [L["text"] for L in figure["labels"] if L["text"] in moments and not L.get("frame")]
            self.assertEqual(named, [], where)

    def test_every_movie_of_moments_is_played_as_the_others_are(self):
        # The same shape as FRW's, Bianchi's and the cosmic string's: frames with a label and a
        # value, one label on the figure that names the frame shown, which is the first, and
        # nothing a client would have to treat on its own.
        keys = {"variable", "seconds", "frames", "turns"}
        for name, view in self.views:
            if "movie" not in view:
                continue
            movie, where = view["movie"], f"{name} {view['id']}"
            self.assertLessEqual(set(movie), keys, where)
            self.assertEqual(movie["seconds"], 5, where)
            for frame in movie["frames"]:
                self.assertLessEqual(set(frame), {"label", "value", "pieces", "rings", "curves", "dots"}, where)
                self.assertNotIn("time", frame, where)
            shown = [L for L in view["figure"]["labels"] if L.get("frame")]
            self.assertEqual([L["text"] for L in shown], [movie["frames"][0]["label"]], where)
        page = (build.ROOT / "_layouts" / "mfs.html").read_text(encoding="utf-8")
        self.assertIn("(movie ? MOVIE_PLAY : '')", page, "the page gives a movie no button")

    def test_the_build_refuses_moments_that_are_not_a_movie(self):
        view = copy.deepcopy(next(v for v in self.named["tolman_bondi"]["views"] if v["id"] == "cloud"))
        build.check_moments("tolman_bondi", view)
        still = {k: v for k, v in view.items() if k != "movie"}
        with self.assertRaises(build.DataError) as raised:
            build.check_moments("tolman_bondi", still)
        self.assertIn("separate pictures", str(raised.exception))
        # A movie that leaves a moment out, or stops short of the last, is refused too.
        short = copy.deepcopy(view)
        last = short["surfaces"][-1]["time"]
        short["movie"]["frames"] = [f for f in short["movie"]["frames"] if f["value"] != last]
        with self.assertRaises(build.DataError) as raised:
            build.check_moments("tolman_bondi", short)
        self.assertIn("does not hold its moment", str(raised.exception))
        changed = copy.deepcopy(view)
        changed["movie"]["frames"][0]["rings"] = []
        with self.assertRaises(build.DataError):
            build.check_moments("tolman_bondi", changed)
        beyond = copy.deepcopy(view)
        beyond["movie"]["frames"].append(dict(beyond["movie"]["frames"][-1], value=last + 1, label="later"))
        with self.assertRaises(build.DataError) as raised:
            build.check_moments("tolman_bondi", beyond)
        self.assertIn("from its first moment to its last", str(raised.exception))
        unordered = copy.deepcopy(view)
        unordered["surfaces"].reverse()
        with self.assertRaises(build.DataError):
            build.check_moments("tolman_bondi", unordered)
        # One moment is no sequence, with a movie, as the cosmic string's cone, or without.
        for metric_id in ("cosmic_string", "schwarzschild"):
            for one in self.named[metric_id]["views"]:
                build.check_moments(metric_id, one)

    def test_the_script_cannot_draw_moments_as_separate_pictures(self):
        # view() of the script is run as it stands, without the script's own imports: it refuses
        # more than one surface without a movie, and the script holds no other way to lay
        # surfaces out, every view being made by view().
        source = (build.ROOT / "_tools" / "derivations" / "embedding.py").read_text(encoding="utf-8")
        tree = ast.parse(source)
        functions = {node.name: node for node in tree.body if isinstance(node, ast.FunctionDef)}
        self.assertNotIn("sequence_figure", functions)
        scope = {}
        exec(compile(ast.Module(body=[functions["view"]], type_ignores=[]), "embedding.py", "exec"), scope)

        class Moment:
            def data(self):
                return {"pieces": []}
        made = scope["view"]("one", "One", "$r_s$", [Moment()], {})
        self.assertEqual(len(made["surfaces"]), 1)
        with self.assertRaises(AssertionError) as raised:
            scope["view"]("many", "Many", "$r_s$", [Moment(), Moment()], {})
        self.assertIn("no movie", str(raised.exception))
        with self.assertRaises(AssertionError):
            scope["view"]("many", "Many", "$r_s$", [Moment(), Moment()], {}, movie=None)
        made = scope["view"]("many", "Many", "$r_s$", [Moment(), Moment()], {}, movie={"frames": []})
        self.assertIn("movie", made)
        # Every view the script returns is made by view(), so none escapes that refusal: no
        # function builds the dictionary of a view by hand.
        builders = [node.name for node in ast.walk(tree) if isinstance(node, ast.FunctionDef) and node.name != "view"
                    and any(isinstance(d, ast.Dict) and {"surfaces", "figure"} <= {getattr(k, "value", None) for k in d.keys}
                            for d in ast.walk(node))]
        self.assertEqual(builders, [])


class TurningLightConeFigures(unittest.TestCase):
    """The page turns a figure of light cones in three dimensions with MFS/assets/turn.js, which
    draws it again from the figure's `turn`, its lines, cones, labels and slices in the drawing's
    (X, Y, T), as _tools/README.md, "Turning a figure of light cones", defines them. At the
    figure's own camera it must give back the published figure, layer for layer in the published
    order, and at any other it must keep inside the box at one scale."""

    def setUp(self):
        self.figures = {f"{name}/{figure['id']}": figure for name, data in diagram_files().items()
                        for figures in data.get("projections", {}).values() for figure in figures if "turn" in figure}
        found(self, self.figures, "figure of light cones that turns")

    def test_every_figure_of_light_cones_turns(self):
        checked = {f"{v['metric']}/{v['view']}" for v in turn_check(self)["figures"]}
        self.assertEqual(checked, set(self.figures))
        self.assertEqual(checked, {key for key in {"alcubierre/bubble", "godel/tipping", "gott_time_machine/loop", "kerr/dragging", "kerr_de_sitter/dragging",
                                   "kerr_newman/dragging", "kerr_taub_nut/dragging", "kerr_melvin/dragging", "kundt_waves/fronts",
                                   "near_horizon_extreme_kerr/dragging", "point_particle_2plus1/wedge",
                                   "som_raychaudhuri/tipping", "spinning_string/tipping", "stockum_dust/tipping",
                                   "bonnor_rotating_dust/tipping", "maitra_dust/tipping", "tippett_tsang/ring",
                                   "wormhole_time_machine/trip", "petrov_homogeneous/turning", "draining_bathtub/swirl",
                                   "vuorio_warped_ads/tipping", "cahen_wallach/turning"} if chosen(key.split("/")[0])})

    def test_at_its_own_camera_the_page_draws_the_published_figure(self):
        # Every point the generator does not thin is the published point to the published
        # rounding, a ten thousandth; a line it thins, 0.0005 across, may keep other points of
        # the same curve, so it lies within that of the published line, and a label is where
        # it was published to the rounding.
        for v in turn_check(self)["figures"]:
            where = f"{v['metric']}/{v['view']}"
            self.assertTrue(v["order"], f"{where}: the layers are not the published layers in the published order")
            self.assertLess(v["cones"], 1.5e-4, f"{where}: a cone strays {v['cones']:.1e} from its published place")
            self.assertLess(v["lines"], 7e-4, f"{where}: a line strays {v['lines']:.1e} from its published place")
            self.assertLess(v["slices"], 7e-4, f"{where}: a slice strays {v['slices']:.1e} from its published place")
            self.assertLess(v["labels"], 1e-4, f"{where}: a label strays {v['labels']:.1e} from its published place")
            self.assertTrue(v["shown"], f"{where}: a label is hidden at the start")

    def test_a_turned_figure_keeps_to_its_box_at_one_scale(self):
        for v in turn_check(self)["figures"]:
            where = f"{v['metric']}/{v['view']}"
            for step in v["sweep"]:
                at = f"{where} at azimuth {step['azimuth']}, elevation {step['elevation']}"
                self.assertEqual(step["bad"], 0, f"{at}: a point is not a number")
                self.assertLessEqual(step["outside"], 1e-9, f"{at}: drawn outside the box")
                self.assertGreater(step["scale"], 0, at)
                self.assertLessEqual(step["scale"], 1, at)
            self.assertEqual(v["sweep"][3]["scale"], 1, f"{where} is not at its own scale at its own camera")

    def test_the_turn_projects_onto_the_published_figure(self):
        # Independently of the page: the published apexes, in the order the cones are painted,
        # and the published labels are the figure's points seen from its own camera.
        for where, figure in self.figures.items():
            turn, camera = figure["turn"], figure["camera"]
            a, e = math.radians(camera["azimuth"]), math.radians(camera["elevation"])
            right = (-math.sin(a), math.cos(a), 0.0)
            up = (-math.sin(e) * math.cos(a), -math.sin(e) * math.sin(a), math.cos(e))

            def seen(P):
                return (sum(p * r for p, r in zip(P, right)), sum(p * u for p, u in zip(P, up)))
            apexes = [layer["at"] for layer in figure["layers"] if layer["class"] == "cone-apex"]
            self.assertEqual(len(apexes), len(turn["cones"]), where)
            for cone, at in zip(turn["cones"], apexes):
                self.assertLess(math.dist(seen(cone["apex"]), at), 1e-4, where)
                self.assertEqual(len(cone["rim"]) % turn["ribs"], 0, where)
            self.assertEqual(len(turn["labels"]), len(figure["labels"]), where)
            for place, label in zip(turn["labels"], figure["labels"]):
                self.assertEqual(len(set(place) & {"at", "circle"}), 1, f"{where} {label['text']}")
                if "circle" in place:
                    rho, t = place["circle"]
                    phi = math.radians(camera["azimuth"] + place["angle"])
                    P = (rho * math.cos(phi), rho * math.sin(phi), t)
                else:
                    P = place["at"]
                self.assertLess(math.dist(seen(P), label["at"]), 1e-4, f"{where} {label['text']}")
            # The figure turns about the axis, halfway up everything it draws.
            T = [p[2] for line in turn["lines"] for p in line["points"]]
            T += [p[2] for cone in turn["cones"] for p in [cone["apex"]] + cone["rim"]]
            self.assertEqual(turn["centre"][:2], [0, 0], where)
            self.assertAlmostEqual(turn["centre"][2], (min(T) + max(T)) / 2, delta=1e-6, msg=where)
            self.assertEqual([len(mark["fills"]) for mark in turn["slices"]],
                             [len(mark["fills"]) for mark in figure.get("slices", [])], where)
            self.assertEqual([len(mark["lines"]) for mark in turn["slices"]],
                             [len(mark["lines"]) for mark in figure.get("slices", [])], where)

    def test_only_the_flat_figures_of_light_rays_do_not_turn(self):
        # The cosmic string's beam lies in the plane t = 0 seen from straight above, drawn in the plane's
        # own flat coordinates, and the rays along Bonnor's beam of light in the plane y = 0 seen from the
        # side, t left out, so neither has another side; Lifshitz spacetime's rays lie in its plane y = 0
        # seen the same way, and the light swinging in the trap of Schrodinger spacetime's global chart
        # in its surface X = 0 seen from the side, V left out.
        still = {f"{name}/{figure['id']}" for name, data in diagram_files().items()
                 for figures in data.get("projections", {}).values() for figure in figures if "turn" not in figure}
        # Ernst and Wild's ergoregion is its meridional plane seen from the side, t and phi left out, and so
        # are the ergoregions of the swirling universe and of the black hole in it. The plane waves' light
        # lies in a surface of their Brinkmann and Rosen charts seen from the side, v left out.
        self.assertEqual(still, {key for key in {"cosmic_string/beam", "light_beam/lens", "lifshitz_spacetime/rays",
                                 "point_particle_2plus1/beam", "schrodinger_spacetime/trap", "kerr_melvin/tube",
                                 "maximally_supersymmetric_plane_wave/focus", "swirling_universe_black_hole/ergo",
                                 "cahen_wallach/pull", "cahen_wallach/push"} if chosen(key.split("/")[0])})


@page_source
class TurningUnderTheHand(unittest.TestCase):
    """The page turns every figure by one hand, MfsTurn.dragged(), keyed() and turned(), which
    _tools/turn_check.cjs drags across a surface of revolution, a height over a plane and a
    figure of light cones as the page does: a fifth of the drawing's width across is 36 degrees
    round the axis the way the hand moves, a whole width up or down stops looking straight up
    or straight down the axis, a whole turn and the Home key bring back the published figure,
    and the left arrow turns it 15 degrees."""

    def test_a_drag_turns_every_kind_of_figure(self):
        hand = turn_check(self)["hand"]
        self.assertEqual(set(hand), {"surface", "height", "cones"})
        for kind, steps in hand.items():
            home = steps["home"]
            start = (home["azimuth"], home["elevation"])
            for name, step in steps.items():
                self.assertEqual(step["bad"], 0, f"{kind} {name}")
                self.assertLessEqual(step["outside"], 1e-9, f"{kind} {name}: drawn outside the box")
            across, up, down = steps["across"], steps["up"], steps["down"]
            self.assertAlmostEqual(across["azimuth"], start[0] - 36, places=9, msg=kind)
            self.assertEqual(across["elevation"], start[1], kind)
            self.assertTrue(across["turned"], kind)
            self.assertGreater(across["moved"], 0.1, f"{kind}: a drag across left the drawing where it was")
            self.assertEqual((up["azimuth"], up["elevation"]), (start[0], -90), kind)
            self.assertEqual((down["azimuth"], down["elevation"]), (start[0], 90), kind)
            for tilt in (up, down):
                self.assertTrue(tilt["turned"], kind)
                self.assertGreater(tilt["moved"], 0.1, f"{kind}: a tilt left the drawing where it was")
            self.assertFalse(steps["round"]["turned"], f"{kind}: a whole turn is not the published figure")
            self.assertFalse(home["turned"], kind)
            self.assertEqual(home["moved"], 0, kind)
            self.assertEqual((steps["left"]["azimuth"], steps["left"]["elevation"]), (start[0] + 15, start[1]), kind)
        # The cone facing the reader goes the way the hand goes, by r_c sin 36 degrees, and the
        # left arrow takes it the other way.
        cones = hand["cones"]
        self.assertAlmostEqual(cones["across"]["facing"], math.asinh(1) * math.sin(math.radians(36)), delta=2e-3)
        self.assertLess(cones["left"]["facing"], 0)


# Simpson and Visser's three geometries at r_s = 1.
SV_A = {"bounce": 0.5, "null": 1.0, "wormhole": 2.0}


def sv_rstar(r, a):
    """The tortoise coordinate of Simpson and Visser's black bounce at r_s = 1, dr_*/dr = rho/(rho - 1) with
    rho = sqrt(r^2 + a^2), odd in r: a logarithm of the two horizons where a < 1, a pole at r = 0 where
    a = 1, and two arctangents where a > 1."""
    rho = math.sqrt(r * r + a * a)
    out = r + math.asinh(r / a)
    if a < 1:
        h = math.sqrt(1 - a * a)
        return out + math.log(abs((r - h) * (r - h * rho) / ((r + h) * (r + h * rho)))) / (2 * h)
    if a == 1:
        return out - (1 + rho) / r
    k = math.sqrt(a * a - 1)
    return out + (math.atan(r / k) + math.atan(r / (k * rho))) / k


def bisect(f, lo, hi, steps=200):
    """A root of f on [lo, hi], where f changes sign."""
    flo = f(lo)
    for _ in range(steps):
        mid = 0.5 * (lo + hi)
        if (f(mid) > 0) == (flo > 0):
            lo, flo = mid, f(mid)
        else:
            hi = mid
    return 0.5 * (lo + hi)


# Israel's shell of dust as its diagrams draw it, the shell that falls from rest at infinity, r_s = 1
# and mu = 1/2. In s = sqrt(1 + 8R) its radius, the flat time T inside it and the advanced time v
# outside it are R = (s^2 - 1)/8, T = -(s^3 + 3s - 4)/24 and
# v = -(s^3/3 - s^2 + 5s - 15 - 16 ln((s + 3)/6))/8, which _tools/test_israel_shell.py holds to the
# published metrics and the junction conditions. It crosses r_s at s = 3 and reaches R = 0 at s = 1.
def israel_radius(s):
    return (s * s - 1) / 8


def israel_inner_time(s):
    return -(s ** 3 + 3 * s - 4) / 24


def israel_advanced(s):
    return -(s ** 3 / 3 - s * s + 5 * s - 15 - 16 * math.log((s + 3) / 6)) / 8


ISRAEL_END = israel_advanced(1.0)       # the advanced time at which the shell reaches the centre


def israel_on_slice(w):
    """s where the slice v - r = w meets the shell, for w < ISRAEL_END."""
    return bisect(lambda s: israel_advanced(s) - israel_radius(s) - w, 1.0, 200.0)


def israel_kruskal(s):
    """Kruskal's U = (1 - R) e^(R - v/2) on the shell."""
    return (1 - israel_radius(s)) * math.exp(israel_radius(s) - israel_advanced(s) / 2)


def israel_event(p, q):
    """The v and r of the event outside the shell that the conformal diagram draws at (p, q). The
    outgoing ray p left the shell where the flat retarded time T - R = -((s + 1)^3 - 8)/24 is
    (7/3) tan p, and keeps u = v - 2r_* outside r_s and Kruskal's U inside it; the ingoing ray q
    entered the shell where T + R = -(s - 1)^3/24 is (7/3) tan q, or for q > 0 ends on r = 0 with the
    outgoing ray -q, at U = e^(-v/2)."""
    def left_at(angle):
        return (8 - 56 * math.tan(angle)) ** (1 / 3) - 1
    s_out = left_at(p)
    if q <= 0:
        v = israel_advanced(1 + (-56 * math.tan(q)) ** (1 / 3))
    else:
        v = -2 * math.log(israel_kruskal(left_at(-q)))
    if s_out > 3:
        R = israel_radius(s_out)
        u = israel_advanced(s_out) - 2 * (R + math.log(R - 1))
        r = bisect(lambda x: v - 2 * (x + math.log(x - 1)) - u, 1 + 1e-13, 1e3)
    else:
        U = israel_kruskal(s_out)
        # (1 - r) e^r falls from 1 at r = 0, where the moment ends on the singularity.
        f = lambda x: (1 - x) * math.exp(x - v / 2) - U
        r = 0.0 if f(0) < 1e-6 else bisect(f, 0, 1)
    return v, r


# The charged shell of dust as its diagrams draw it, r_s = 1, r_q = 12/25 and mu = 1/5, written here
# again from the junction condition alone. With E = 1/(2 mu), k = (r_q^2 - mu^2)/(2 mu),
# l = (r_q^2 + mu^2)/(2 mu) and a = E^2 - 1, the shell is R = k(E + cosh eta)/a at the flat time
# T = k(eta + E sinh eta)/a^(3/2) inside it, and its advanced time outside obeys
# dv/d(eta) = (R/sqrt(a))/(beta - Rdot), beta = E - l/R and Rdot = sqrt(a) sinh eta/(E + cosh eta),
# with v = r_* at the turn, eta = 0; on the way out v = -v(-eta) + 2 r_*.
CHARGED_RQ, CHARGED_MU, CHARGED_RP, CHARGED_RM = 12 / 25, 1 / 5, 16 / 25, 9 / 25
CHARGED_E = 1 / (2 * CHARGED_MU)
CHARGED_K = (CHARGED_RQ ** 2 - CHARGED_MU ** 2) / (2 * CHARGED_MU)
CHARGED_L = (CHARGED_RQ ** 2 + CHARGED_MU ** 2) / (2 * CHARGED_MU)
CHARGED_A = CHARGED_E ** 2 - 1
CHARGED_ETA_P = math.acosh(CHARGED_A * CHARGED_RP / CHARGED_K - CHARGED_E)
CHARGED_ETA_M = math.acosh(CHARGED_A * CHARGED_RM / CHARGED_K - CHARGED_E)


def charged_radius(eta):
    return CHARGED_K * (CHARGED_E + math.cosh(eta)) / CHARGED_A


def charged_inner_time(eta):
    return CHARGED_K * (eta + CHARGED_E * math.sinh(eta)) / CHARGED_A ** 1.5


def charged_rstar(r):
    """Reissner-Nordstrom's tortoise coordinate, zero at r = 0."""
    gap = CHARGED_RP - CHARGED_RM
    return (r + CHARGED_RP ** 2 / gap * math.log(max(abs(r / CHARGED_RP - 1), 1e-300))
            - CHARGED_RM ** 2 / gap * math.log(max(abs(r / CHARGED_RM - 1), 1e-300)))


def _charged_rate(eta):
    R = charged_radius(eta)
    speed = math.sqrt(CHARGED_A) * math.sinh(eta) / (CHARGED_E + math.cosh(eta))
    return R / math.sqrt(CHARGED_A) / (CHARGED_E - CHARGED_L / R - speed)


def charged_advanced(eta):
    """v on the shell, by Simpson's rule from the turn."""
    back = -abs(eta)
    n = 2 * max(8, int(abs(back) * 400))
    h = -back / n
    total = _charged_rate(back) + _charged_rate(0.0) + sum((4 if i % 2 else 2) * _charged_rate(back + i * h) for i in range(1, n))
    way_in = charged_rstar(charged_radius(0.0)) - total * h / 3
    return way_in if eta <= 0 else -way_in + 2 * charged_rstar(charged_radius(eta))


def charged_on_slice(w):
    """eta where the slice v - r = w meets the shell."""
    return bisect(lambda eta: charged_advanced(eta) - charged_radius(eta) - w, -12.0, CHARGED_ETA_M - 1e-12, 60)


def charged_event(p, q, length):
    """The v and r of the event outside the shell that the conformal diagram draws at (p, q), with
    length the L of p, q = arctan((cT -+ r)/L) inside the shell. The outgoing ray p left the shell
    where T - R = L tan p, and keeps v - 2r_* in whichever of the three stretches of r it left in;
    the ingoing ray q meets the shell where T + R = L tan q, and has the shell's v there."""
    left = bisect(lambda eta: charged_inner_time(eta) - charged_radius(eta) - length * math.tan(p), -12.0, 12.0, 80)
    met = bisect(lambda eta: charged_inner_time(eta) + charged_radius(eta) - length * math.tan(q), -12.0, 12.0, 80)
    v = charged_advanced(met)
    keeps = charged_advanced(left) - 2 * charged_rstar(charged_radius(left))
    if left < -CHARGED_ETA_P:
        lo, hi = CHARGED_RP + 1e-13, 1e3
    elif left < -CHARGED_ETA_M:
        lo, hi = CHARGED_RP - 1e-13, CHARGED_RM + 1e-13
    else:
        lo, hi = 1e-9, CHARGED_RM - 1e-13
    # v - 2r_* falls with r outside r_+ and inside r_-, and rises with it between the horizons.
    return v, bisect(lambda r: -(v - 2 * charged_rstar(r) - keeps), lo, hi, 200)


class Lukewarm:
    """The lukewarm charged black hole in de Sitter space as every diagram draws it, r_s = 1, r_q = 1/2 and
    Lambda = 27/64, H = 3/8: f = (1 - 1/(2r))^2 - 9r^2/64 has the roots 2, 2/3, (2 sqrt 7 - 4)/3 and
    -(2 sqrt 7 + 4)/3, and everything below is written from them and from nothing the generators hold."""
    ROOTS = (2.0, 2 / 3, (2 * math.sqrt(7) - 4) / 3, -(2 * math.sqrt(7) + 4) / 3)
    H = 3 / 8
    KAPPA = 3 / 16                  # the surface gravity of r_+ and of r_c

    @classmethod
    def rstar(cls, r):
        """sum_i ln|1 - r/r_i|/f'(r_i), with 1/f'(r_i) = -64 r_i^2/(9 prod_j (r_i - r_j))."""
        return sum(math.log(abs(1 - r / a)) * (-64 * a * a / (9 * math.prod(a - b for b in cls.ROOTS if b != a)))
                   for a in cls.ROOTS)

    @classmethod
    def static_t(cls, tau, r):
        """cT of the cosmological event at tau with the areal radius r: ln|H tau|/H + F(r) - F(1.298),
        F the sum of c_a ln|r - a| over a = 1/2 and the roots, c_a = -(8/3) a^4/prod_b (a - b)."""
        poles = (0.5,) + cls.ROOTS

        def F(x):
            return sum(math.log(abs(x - a)) * (-(8 / 3) * a ** 4 / math.prod(a - b for b in poles if b != a))
                       for a in poles)
        return math.log(abs(cls.H * tau)) / cls.H + F(r) - F(1.2979848366419)

    @classmethod
    def inside_t(cls):
        """The static time of the moment embedded inside r_-, on the cosmological chart's time."""
        return cls.static_t(-0.5 / cls.H, 0.35)

    @classmethod
    def tau_of(cls, r, at, sign):
        """tau of the event of static time `at` and radius r on the cosmological plane."""
        return sign * math.exp(cls.H * (at - cls.static_t(1 / cls.H, r))) / cls.H

    @classmethod
    def radius_at(cls, rho, at, sign, lo, hi):
        """The areal radius, between lo and hi, at which the moment `at` meets the comoving rho."""
        def miss(r):
            return (r - 0.5) / (cls.H * cls.tau_of(r, at, sign)) - rho
        lo, hi = lo + 1e-12 * (hi - lo), hi - 1e-12 * (hi - lo)
        if miss(lo) * miss(hi) > 0:
            return lo if abs(miss(lo)) < abs(miss(hi)) else hi
        return bisect(miss, lo, hi) if miss(lo) < 0 else bisect(lambda r: -miss(r), lo, hi)


# The angular radius of a cell of Lindquist and Wheeler's lattice of eight, 8 (2 psi - sin 2 psi) = 2 pi.
LW_PSI = bisect(lambda x: 8 * (2 * x - math.sin(2 * x)) - 2 * math.pi, 0.0, math.pi / 2)


def novikov_uv(s, tau):
    """Kruskal's U and V, r_s = 1, of Novikov's shell at rest at the areal radius s^2 + 1 on the far
    sheet, at its proper time tau since the release: slices.novikov_sheets, by hand."""
    R = s * s + 1
    eta = bisect(lambda e: 0.5 * R * math.sqrt(R) * (e + math.sin(e)) - tau, 0.0, math.pi) if tau else 0.0
    r = R * math.cos(eta / 2) ** 2
    phase = s * (eta + 0.5 * R * (eta + math.sin(eta)))
    return ((math.sin(eta / 2) - s * math.cos(eta / 2)) * math.exp((r - phase) / 2),
            (s * math.cos(eta / 2) + math.sin(eta / 2)) * math.exp((r + phase) / 2))


def ruban_sigma(eta, mu=1 / math.pi):
    """The conformal time of Ruban's dust at epsilon = 1 and r_s = 1 from its greatest expansion to eta,
    the integral of sin^2(e/2)/a with a = (1 - mu e) cot(e/2) + 2 mu, by Simpson's rule."""
    n = 2000
    f = [math.sin(e / 2) ** 2 / ((1 - mu * e) / math.tan(e / 2) + 2 * mu)
         for e in (math.pi + (eta - math.pi) * k / n for k in range(n + 1))]
    return (eta - math.pi) / n / 3 * (f[0] + f[-1] + 4 * sum(f[1:-1:2]) + 2 * sum(f[2:-1:2]))


def novikov_t(R, tau):
    """A shell of dust released from rest at areal radius R, r_s = 1, at its proper time tau:
    its r and Schwarzschild t, from the cycloid and Misner, Thorne and Wheeler's (31.10)."""
    eta = bisect(lambda e: 0.5 * R * math.sqrt(R) * (e + math.sin(e)) - tau, 0.0, math.pi)
    k, tan = math.sqrt(R - 1), math.tan(eta / 2)
    r = 0.5 * R * (1 + math.cos(eta))
    return r, math.log(abs((k + tan) / (k - tan))) + k * (eta + 0.5 * R * (eta + math.sin(eta)))


# Bonnor and Vaidya's charged shell as its drawings take it, in units of its mass: the charge, the
# horizons of the charged side, the surface gravity of the outer one, and the radius where the
# shell's energy density changes sign.
BV_Q = 0.96
BV_RP, BV_RM = 1.28, 0.72
BV_KAPPA = (BV_RP - BV_RM) / (2 * BV_RP ** 2)
BV_TURN = BV_Q ** 2 / 2


def bv_rstar(r):
    """Reissner-Nordstrom's tortoise coordinate at M = 1 and q = 24/25, zero at r = 0."""
    return r + (BV_RP ** 2 * math.log(abs(r / BV_RP - 1)) - BV_RM ** 2 * math.log(abs(r / BV_RM - 1))) / (BV_RP - BV_RM)


def bv_crossing(u):
    """p of the outgoing ray that crosses the falling shell v = 0 at r = -u/2, as the tower's ingoing
    map has it: -arctan e^(2 k r_*) outside r_+, arctan e^(2 k r_*) between the horizons, and
    pi - arctan e^(2 k r_*) inside r_-."""
    r = -u / 2
    x = 2 * BV_KAPPA * bv_rstar(r)
    g = math.pi / 2 if x > 700 else math.atan(math.exp(x))
    return -g if r > BV_RP else g if r > BV_RM else math.pi - g


def bv_radius(v, p):
    """r of the event of the charged side with the advanced time v > 0 and the tower's p, or None within
    1e-3 of a cell's edge, where the map is too steep to carry back."""
    for lo, hi, a, b in ((BV_RP, 60.0, -math.pi / 2, 0.0), (BV_RM, BV_RP, 0.0, math.pi / 2), (0.0, BV_RM, math.pi / 2, math.pi)):
        if a + 1e-3 < p < b - 1e-3:
            if lo == 0.0:
                star = (v + math.log(math.tan(math.pi - p)) / BV_KAPPA) / 2
            else:
                star = (v + math.log(math.tan(abs(p))) / BV_KAPPA) / 2
            return bisect(lambda r: bv_rstar(r) - star, lo + 1e-12, hi - 1e-12)
    return None


# p of each outcome of Roberts's collapse, as its drawings name them.
# The sphere r = 8 ell of Bartnik and McKinnon's soliton with one zero in its three other radial coordinates.
BARTNIK_MCKINNON_REACH = {"isotropic": 7.147807396, "tortoise": 24.219163201, "flow": 1.966805652}
ROBERTS_P = {"disperses": 0.9, "threshold": 1.0, "collapses": 2.0}


def murphy_scale(t):
    """The scale factor of Murphy's universe for radiation at the proper time t, in units of 1/H_s: the
    root of ln a + a^2/2 = t, by bisection on the logarithm of a."""
    lo, hi = -60.0, 10.0
    for _ in range(200):
        mid = (lo + hi) / 2
        if mid + math.exp(2 * mid) / 2 < t:
            lo = mid
        else:
            hi = mid
    return math.exp((lo + hi) / 2)


def lqc_bounce_eta(t, n=2000):
    """The conformal time of the bounce of loop quantum cosmology at the proper time t, in units of t_b:
    the integral of dt/(1 + t^2)^(1/6) from the bounce, by Simpson's rule."""
    h = t / n
    f = [(1 + (k * h) ** 2) ** (-1 / 6) for k in range(n + 1)]
    return h / 3 * (f[0] + f[-1] + 4 * sum(f[1:-1:2]) + 2 * sum(f[2:-1:2]))


def kopczynski_trautman_eta(t, n=2000):
    """The conformal time of Kopczynski and Trautman's universe at the proper time t, in units of l:
    the integral of dt/(1 + t^2)^(1/3) from the turn, by Simpson's rule."""
    h = t / n
    f = [(1 + (k * h) ** 2) ** (-1 / 3) for k in range(n + 1)]
    return h / 3 * (f[0] + f[-1] + 4 * sum(f[1:-1:2]) + 2 * sum(f[2:-1:2]))


def hotta_tanaka_moment(tau, system=None, fixed=None):
    """The sphere of Hotta and Tanaka's Kruskal chart through the embedded ring at the proper time
    tau, at 8GE/c^4 = a = 1 with the ring at 30 degrees, as its (u, v), or as the event where it
    meets the plane of another chart's view, in that view's drawn (X, Y): the ring keeps
    Z_0 + Z_1 = sinh(tau), and behind the shock its distance from the axis and its height along it
    are Podolsky and Ortaggio's, which fix Z_0 - Z_1 = (R^2 - 1)/(Z_0 + Z_1)."""
    theta0, g = math.pi / 6, 0.5
    U = math.sinh(tau)
    V = U
    if U > 0:
        Z = math.sin(theta0) * math.sqrt(1 + U * U) - g / math.sin(theta0) * U
        Z4 = math.cos(theta0) * math.sqrt(1 + U * U) + g * math.log(1 / math.tan(theta0 / 2)) * U
        V = (Z * Z + Z4 * Z4 - 1) / U
    x = (math.sqrt(1 + U * V) - 1) / (math.sqrt(1 + U * V) + 1)
    u, v = U * (1 - x) / 2, V * (1 - x) / 2
    if system is None:
        return u, v
    if system == "kruskal":
        return (v - u) / 2, (u + v) / 2
    if system == "global":
        return math.pi / 2 + math.atan(v) - math.atan(u), math.atan(u) + math.atan(v) + math.pi / 2
    if system == "null_cylindrical":
        k = 1 + fixed * fixed / 4
        omega = k if U * V == 0 else (-1 + math.sqrt(1 + U * V * k)) / (U * V / 2)
        return (V - U) * omega / 2, (U + V) * omega / 2
    # The conformally flat chart near a particle: eta^2 - rho^2 = 2 eta U and rho cos(theta) - 1 = eta (V - U)/2.
    c, k = math.cos(fixed), (V - U) / 2
    A, B = c * c - k * k, k + c * c * U
    eta = next(e for e in ((B + sign * math.sqrt(B * B + A)) / A for sign in (1, -1)) if e < 0 and 1 + k * e > 0)
    return (1 + k * eta) / c, eta


class Exceptions(unittest.TestCase):
    """Each spacetime's exceptions to the tests every spacetime is held to stand in a file of its own,
    _tools/derivations/exceptions/<id>.json, every one with its reason and every one naming a drawing
    and a moment that are published, so that none outlives what it excuses."""

    # The fields that hold groups of exceptions, each a reason and what it excuses, by drawing.
    GROUPS = {"hidden": list, "hidden_views": dict, "hidden_surfaces": dict, "hidden_moments": dict}

    def test_every_file_is_a_spacetimes_and_every_exception_carries_its_reason(self):
        spacetimes = {p.stem for p in build.METRICS_DIR.glob("*.json")}
        for name, data in EXCEPTIONS.items():
            if not chosen(name):
                continue
            where = f"exceptions/{name}.json"
            self.assertIn(name, spacetimes, f"{where} is of no spacetime")
            self.assertLessEqual(set(data), {*self.GROUPS, "flat", "flat_moments", "height"}, where)
            for field in ("flat", "height"):
                if field in data:
                    self.assertTrue(isinstance(data[field], str) and data[field].strip(), f"{where} {field}")
            for field, kind in self.GROUPS.items():
                places = []
                for group in data.get(field, []):
                    self.assertEqual(set(group), {"reason", "drawings"}, f"{where} {field}")
                    self.assertTrue(isinstance(group["reason"], str) and group["reason"].strip(), f"{where} {field}")
                    self.assertIsInstance(group["drawings"], kind, f"{where} {field}")
                    self.assertTrue(group["drawings"], f"{where} {field}")
                    places += list(group["drawings"])
                self.assertEqual(len(places), len(set(places)), f"{where} {field} names a drawing twice")
                for place in places:
                    self.assertEqual(place.removeprefix("conformal ").split("/")[0], name, f"{where} {field} {place}")
            for group in data.get("flat_moments", []):
                self.assertEqual(set(group), {"reason", "moments"}, f"{where} flat_moments")
                self.assertTrue(isinstance(group["reason"], str) and group["reason"].strip(), f"{where} flat_moments")
                self.assertTrue(group["moments"], f"{where} flat_moments")

    def test_every_exception_names_a_published_drawing_and_moment(self):
        for name, data in EXCEPTIONS.items():
            if not chosen(name):
                continue
            where = f"exceptions/{name}.json"
            drawings = set()
            if (build.DIAGRAMS_DIR / f"{name}.json").exists():
                drawings |= {f"{name}/{system}/{view['id']}" for part in ("systems", "projections")
                             for system, views in diagram_file(name).get(part, {}).items() for view in views}
            if (build.CONFORMAL_DIR / f"{name}.json").exists():
                drawings |= {f"conformal {name}/{view['id']}" for view in conformal_file(name)["views"]}
            views = ({view["id"]: view for view in embedding_file(name)["views"]}
                     if (build.EMBEDDING_DIR / f"{name}.json").exists() else {})
            for field in self.GROUPS:
                for group in data.get(field, []):
                    for place in group["drawings"]:
                        self.assertIn(place, drawings, f"{where} {field} names no published drawing")
                    for place, held in (group["drawings"].items() if field != "hidden" else ()):
                        for moment in held:
                            view, *at = [moment] if field == "hidden_views" else moment
                            self.assertIn(view, views, f"{where} {field} {place}: no embedding view {view}")
                            if at:
                                self.assertLess(at[0], len(views[view]["surfaces"]), f"{where} {field} {place} {moment}")
            for group in data.get("flat_moments", []):
                for view, k in group["moments"]:
                    self.assertIn(view, views, f"{where} flat_moments: no embedding view {view}")
                    self.assertLess(k, len(views[view]["surfaces"]), f"{where} flat_moments {view} {k}")
            if "flat" in data or "height" in data:
                self.assertTrue(views, f"{where} excuses an embedding diagram {name} does not have")


class Slices(unittest.TestCase):
    """Every moment an embedding diagram is cut from is drawn on its spacetime's other diagrams
    where it lies, from the numbers in the files alone, and on nothing else."""

    # The drawings of a spacetime on which no moment of its embedding lies, as other universes, another
    # cloud, the time reversed shell and cylinders where no surface of constant t is a moment of space,
    # from the `hidden` of each spacetime's exceptions.
    HIDDEN = {place for _, group in exceptions("hidden") for place in group["drawings"]}

    def setUp(self):
        self.diagrams, self.conformal, self.embedding = diagram_files(), conformal_files(), embedding_files()

    def drawings(self):
        """Every flat view, figure and conformal view, with its place and its spacetime."""
        for metric_id, data in self.diagrams.items():
            for part in ("systems", "projections"):
                for system_id, views in data.get(part, {}).items():
                    for view in views:
                        yield metric_id, f"{metric_id}/{system_id}/{view['id']}", view
        for metric_id, data in self.conformal.items():
            for view in data["views"]:
                yield metric_id, f"conformal {metric_id}/{view['id']}", view

    def moment(self, metric_id, mark):
        """The embedding view a slice marks and its surface, a ring of a stack standing for the
        surface at that ring's time."""
        view = next(v for v in self.embedding[metric_id]["views"] if v["id"] == mark["view"])
        surface = view["surfaces"][mark["surface"]]
        if "curve" in mark:
            ring = surface["curves"][mark["curve"]]
            surface = dict(surface, time=ring["time"], label=ring["label"])
        return view, surface

    # The moments of a spacetime's embedding views, a stack's each ring it marks at a time; a spacetime
    # WITHOUT_EMBEDDING has none.
    def every(self, metric_id):
        out = []
        for v in self.embedding.get(metric_id, {"views": []})["views"]:
            for i, surface in enumerate(v["surfaces"]):
                timed = [k for k, c in enumerate(surface.get("curves", [])) if "time" in c]
                out += [(v["id"], i, k) for k in timed] or [(v["id"], i, None)]
        return out

    # The embedding views a drawing does not mark though it marks others, as another spacetime of the
    # same line element is drawn there, from `hidden_views`.
    HIDDEN_VIEWS = {place: set(views) for _, group in exceptions("hidden_views")
                    for place, views in group["drawings"].items()}
    # The single moments a drawing does not hold, as (view, surface, ring) of every() with ring None for a
    # surface of its own, from `hidden_moments`.
    HIDDEN_MOMENTS = {place: {tuple(moment) for moment in moments} for _, group in exceptions("hidden_moments")
                      for place, moments in group["drawings"].items()}
    # The surfaces of a view that a drawing does not mark though it marks the view's others, as (view,
    # surface), from `hidden_surfaces`.
    HIDDEN_SURFACES = {place: {tuple(surface) for surface in surfaces} for _, group in exceptions("hidden_surfaces")
                       for place, surfaces in group["drawings"].items()}

    def reach(self, surface, system=None, reference=False):
        xs = [x for piece in surface["pieces"] if "points" in piece and (reference or not piece.get("reference"))
              and (system is None or piece["system"] == system) for x in (piece["points"][0][0], piece["points"][-1][0])]
        return min(xs), max(xs)

    def test_every_moment_appears_on_every_drawing_it_lies_on_and_on_no_other(self):
        drawn = 0
        for metric_id, where, view in self.drawings():
            marks = [(m["view"], m["surface"], m.get("curve")) for m in view.get("slices", [])]
            if where in self.HIDDEN:
                self.assertEqual(marks, [], where)
                continue
            every = [m for m in self.every(metric_id) if m[0] not in self.HIDDEN_VIEWS.get(where, set())
                     and m[:2] not in self.HIDDEN_SURFACES.get(where, set())
                     and m not in self.HIDDEN_MOMENTS.get(where, set())]
            self.assertEqual(marks, every, where)
            drawn += len(marks)
        counted(self, drawn, 100, "slices marked on a drawing")

    def test_every_slice_is_stamped_with_the_surface_it_marks(self):
        for metric_id, where, view in self.drawings():
            for mark in view.get("slices", []):
                target, _ = self.moment(metric_id, mark)
                self.assertEqual(build.embedding_moment_version(target, mark["surface"], mark.get("curve")),
                                 mark["version"], where)
                self.assertTrue(mark["lines"] or mark["points"] or mark["fills"], where)
                self.assertTrue(mark["label"].count("$") % 2 == 0 and mark["label"], where)

    def flat_expected(self, key, view, surface, mark):
        """The moment on a flat view as its drawn axes put it: Y as a function of X, and the
        ends a line of it may have short of the box."""
        t = surface.get("time")
        if key == "early_braneworlds/magnetic_general/radial":
            # The same magnetic surface: beta=(1+rho^2/a^2)^(1/3), with a=1.
            return (lambda X: 0.0), [(1+rho*rho)**(1/3) for rho in self.reach(surface)]
        if key.startswith("shell_void/"):
            # One fixed slice: exterior t=2t0, interior T=1.0043869706139898t0,
            # the matched clock reading in the independently recorded trajectory audit.
            inside = key.split("/")[1] == "interior"
            level = 1.0043869706139898 if inside else 2.0
            return (lambda X: level), list(self.reach(surface, "interior" if inside else "exterior"))
        if key.startswith("extended_inflation/"):
            # The fixed reference surface is t=t_0=1, eta=-1 at p=2.
            level = 1.0 if "/cosmic/" in key else -1.0
            return (lambda X: level), list(self.reach(surface))

        if key.startswith("electron_bubble/"):
            # At the shell r0=m=1, the two clocks agree at zero.
            # Its axial intersection is one pole in r or both poles in Cartesian z.
            return (lambda X: 0.0), ([-1.0, 1.0] if "/cartesian/" in key else [1.0])

        if key.startswith("exponential_scalar_universes/"):
            chart = key.split("/")[1]
            p = {"inflation": 2, "tracking": 2/3, "coasting": 1}[mark["view"]]
            clock = 1 if chart.startswith("cosmic") else 0 if p == 1 else 1/(1-p)
            lo, hi = self.reach(surface)
            return (lambda X: clock), [lo if chart == "cosmic_spherical" else -hi, hi]

        if key.startswith("hosotani_universe/"):
            # One sphere at maximum expansion: t=ell=1, eta=pi/2, over the full polar distance.
            height = 1.0 if "/cosmic_hyperspherical/" in key else math.pi / 2
            return (lambda X: height), list(self.reach(surface))
        if key.startswith("charach_malin/"):
            # The single embedded slice has t=1; logarithmic views plot -tau=ln(t)=0.
            return (lambda X: 0.0 if "_logarithmic/" in key else 1.0), list(self.reach(surface))
        if key.startswith("chodos_detweiler/"):
            # The cylinder cuts x at half a period and keeps the whole w circle.
            # X=2x and W=w/2 in the displayed rescaled chart; time is unchanged.
            _, system, plane = key.split("/")
            growing = plane in ("tx", "tX")
            end = math.pi if growing else 2 * math.pi
            if system == "rescaled":
                end *= 2 if growing else 0.5
            return (lambda X: t), [0.0, end]
        if key.startswith("thakurta/"):
            # The one embedded surface is eta=eta_0=1, with t=eta^3/3 and a=1.
            lo, hi = self.reach(surface, "conformal")
            system = key.split("/")[1]
            if system in ("conformal", "cosmic", "areal"):
                return (lambda X: 1.0 if system == "conformal" else 1 / 3), [lo, hi]
            if system == "isotropic":
                return (lambda X: 1.0), [(r - 1 + math.sqrt(r * (r - 2))) / 2 for r in (lo, hi)]
            ends = [r + 2 * math.log(r / 2 - 1) for r in (lo, hi)]
            if system == "tortoise":
                return (lambda X: 1.0), ends
            # U=-exp(-(1-x)/4), V=exp((1+x)/4), drawn in X=(V-U)/2, Y=(V+U)/2.
            return (lambda X: math.tanh(0.25) * X), [math.exp(x / 4) * math.cosh(0.25) for x in ends]

        if key.startswith("moncrief_taub_nut/"):
            # t=1 is T=1 and chi=psi for both logarithmic shifts.
            return (lambda X: 1.0), [0.0, 2*math.pi]
        if key.startswith("stiff_fluid_universes/"):
            # The one embedded symmetry plane is t=1/2,x=0, a point on each quotient.
            # chi=x+t^2/2=1/8 in the fluid labels; the characteristic view displays x,t.
            location = 0.125 if "/comoving/" in key else 0.0
            return (lambda X: 0.5), [location]
        if key.startswith("draining_bathtub/"):
            # The laboratory's t = 0 and the Kerr-like chart's T = 0, with T = t - ln(r^2 - 1)/2 for the
            # drain at |A| = c = 1: each a level line in its own chart and a curve in the other's.
            own = (mark["view"] == "plane") == (key.split("/")[1] == "laboratory")
            sign = 1 if mark["view"] == "funnel" else -1
            return (lambda X: 0.0 if own else sign * 0.5 * math.log(X * X - 1)), list(self.reach(surface))
        if key.startswith("unruh_acoustic_hole/"):
            # The laboratory's t = 0 and Unruh's tau = 0, with t - tau = ln|(r - 1)/(r + 1)|/4 + arctan(r)/2
            # at r_0 = c = 1: each a level line in its own chart and a curve in the other's.
            own = (mark["view"] == "space") == (key.split("/")[1] == "laboratory")
            sign = 1 if mark["view"] == "funnel" else -1
            # Unruh's tau = 0 lies outside the horizon alone, running off to minus infinity in t toward it, so
            # inside it the curve is taken at minus infinity and on it at a finite stand-in, which leaves a
            # slope taken across the curve's end as steep as the curve is there.
            def Y_of(X):
                if own:
                    return 0.0
                if X < 1:
                    return sign * -math.inf
                return sign * (math.log(max((X - 1) / (X + 1), 1e-300)) / 4 + math.atan(X) / 2)
            return Y_of, list(self.reach(surface))
        if key.startswith("bardeen_perturbations/"):
            # The perturbed Friedmann universe's moment eta = L, L = 1, where a = 1: a line in the longitudinal
            # chart, and in the synchronous chart of the same flat universe of dust the curve
            # eta (1 + e^(-x^2)/30) = 1, from eta -> eta + d_eta E with E = eta^2 e^(-x^2)/60.
            hi = self.reach(surface)[1]
            if key.endswith("/synchronous/well"):
                return (lambda X: 1 / (1 + math.exp(-X * X) / 30)), [-hi, hi]
            return (lambda X: 1.0), [-hi, hi]
        if key == "newman_tamburino/cylindrical/line":
            # The cylindrical limit's front, the surface of x and y at u = 0 and r = ell = 1, meets the
            # plane of u and r at x = 1 in the one event drawn at (r, cu + r) = (1, 1).
            return (lambda X: 1.0), [1.0]
        if key == "morgan_morgan/oblate_spheroidal/plane":
            # The plane z = 0 outside the rim, embedded out to Weyl's rho: xi = sqrt(rho^2/a^2 - 1).
            return (lambda X: 0.0), [math.sqrt(self.reach(surface)[1] ** 2 - 1)]
        if key.startswith("bach_weyl_ring/"):
            # The plane z = 0 on one side of the ring, embedded between two of Weyl's rho, in each chart's
            # own coordinate; the axis meets the disc inside the ring at its centre, sigma = pi in the
            # toroidal chart.
            lo, hi = self.reach(surface)
            of_rho = {"toroidal/axis": lambda r: math.pi, "toroidal/outer": lambda r: 2 * math.atanh(1 / r),
                      "toroidal/inner": lambda r: 2 * math.atanh(r),
                      "oblate_spheroidal/plane": lambda r: math.sqrt(r * r - 1),
                      "oblate_spheroidal/disc": lambda r: math.sqrt(1 - r * r)}.get(key.split("/", 1)[1], lambda r: r)
            return (lambda X: 0.0), sorted(of_rho(r) for r in (lo, hi))
        if key.startswith("bondi_sachs/"):
            # The sphere cu = 10 m_0, r = 10 m_0 is one event of each plane: r = 10 and cu + r = 20 on
            # Bondi's axes, 100 m_0 l = 10 and cu = 10 on the chart of l = 1/r.
            self.assertEqual((mark["lines"], len(mark["points"])), ([], 1), key)
            X0, X1 = view["box"][:2]
            self.assertAlmostEqual(X0 + mark["points"][0][0] * (X1 - X0), 10.0, delta=2e-3, msg=key)
            return (lambda X: 10.0 if key.endswith("compactified/equator") else 20.0), [10.0]
        if key.startswith("supertranslation_hair/"):
            # Compere and Long's shell t = 0, rho = rho_k at ell = 1 is one event of each equator:
            # (rho_k, 0) on the static plane, and r = sqrt((rho_k - 1)^2 - 9/4) with cu = -rho_k or
            # cv = rho_k on the Bondi planes, since t = u + W + 1 = v - W - 1 and rho = W + 1 there.
            self.assertEqual((mark["lines"], len(mark["points"])), ([], 1), key)
            rho = surface["time"]
            X_at, Y_at = {"static": (rho, 0.0), "bondi_retarded": (math.sqrt((rho - 1) ** 2 - 2.25), -rho),
                          "bondi_advanced": (math.sqrt((rho - 1) ** 2 - 2.25), rho)}[key.split("/")[1]]
            X0, X1 = view["box"][:2]
            self.assertAlmostEqual(X0 + mark["points"][0][0] * (X1 - X0), X_at, delta=2e-3, msg=key)
            return (lambda X: Y_at), [X_at]
        if key.startswith("hotta_tanaka/"):
            # The event where the sphere through the ring meets the view's plane, 15 degrees from a
            # particle or on the equator, rho = a/4 or 2a in the null cylindrical chart.
            _, system, plane = key.split("/")
            fixed = ({"equator": 2.0, "near": 0.25} if system == "null_cylindrical"
                     else {"equator": math.pi / 2, "near": math.pi / 12})[plane]
            X_at, Y_at = hotta_tanaka_moment(t, system, fixed)
            return (lambda X: Y_at + (X - X_at)), [X_at]
        if key.startswith("siklos/"):
            # Siklos's wave front u = v = 0, embedded on the disc out to the coordinate radius the profile
            # ends at: the event (0, 0) of a plane of u and v, and on Kaigorodov's planes the line of zero
            # v, t, U or V between the ends of the diameter eta = 0, where x = L(2L - xi)/(2L + xi), written
            # as x = L e^(-rho/L) = L e^(2Z) in the horospheric and homogeneous charts.
            if view["plane"][0] == "u":
                return (lambda X: 0.0), [0.0]
            top = self.reach(surface)[1]
            of_x = {"kaigorodov_horospheric": lambda x: -math.log(x),
                    "kaigorodov_homogeneous": lambda x: math.log(x) / 2}.get(key.split("/")[1], lambda x: x)
            return (lambda X: 0.0), sorted(of_x(x) for x in ((2 - top) / (2 + top), (2 + top) / (2 - top)))
        if key == "nordstrom_scalar/dust/radial":
            # A moment of the inertial time, from the centre to the galaxy at r = L.
            return (lambda X: t), list(self.reach(surface))
        if key.startswith("small_universes/torus"):
            # A moment of the torus meets the plane y = z = 0 along its whole circle of x: the line of its
            # time across every cell, at the conformal time 2 (3t/2)^(1/3) on the conformal chart.
            height = 2 * (1.5 * t) ** (1 / 3) if key.split("/")[1] == "torus_conformal" else t
            return (lambda X: height), None
        if key.startswith("tilted_universes/"):
            # Farnsworth's dust is embedded on its surface eta = 3, which is u = (sinh 3 - 3)/2 W, from its rim up.
            height = 3.0 if key.split("/")[1] == "farnsworth" else (math.sinh(3.0) - 3.0) / 2
            return (lambda X: height), list(self.reach(surface))
        if key == "poincare_dodecahedral/conformal/near":
            # Today in Luminet's universe, eta = 0.38882 of the conformal chart, along the line from us.
            return (lambda X: 0.38882), list(self.reach(surface))
        if key == "small_universes/horn/along":
            # The horn is embedded at the moment a = a_0, ct = (6/5 - (11/50) ln 11) a_0, from its rim up.
            return (lambda X: 1.2 - 0.22 * math.log(11)), list(self.reach(surface))
        if key.startswith("kopczynski_trautman/"):
            # A moment of the dust's time, out to the dust at r = l: the line of that time on the comoving
            # charts, through the centre on the Cartesian one, and of eta(t) on the conformal chart.
            lo, hi = self.reach(surface)
            system = key.split("/")[1]
            height = kopczynski_trautman_eta(t) if system == "conformal" else t
            return (lambda X: height), ([-hi, hi] if system == "comoving_cartesian" else [lo, hi])
        if key.startswith("lqc_bounce/"):
            # A moment of the comoving observers' time, out to the observer at r = c t_b: the line of that time
            # on the comoving charts, through the centre on the Cartesian ones, and of tau = asinh(t) on the
            # harmonic chart.
            lo, hi = self.reach(surface)
            system = key.split("/")[1]
            height = math.asinh(t) if system == "harmonic" else t
            return (lambda X: height), ([lo, hi] if system == "comoving_spherical" else [-hi, hi])
        if key.startswith("bekenstein_universes/"):
            # A moment of the conformal time eta, which each chart draws upward, from the centre to the
            # embedding's reach: the pole to the antipode in the closed chart, out to chi = 1 in the open one.
            return (lambda X: t), list(self.reach(surface))
        if key.startswith("murphy_viscous_universe/"):
            # A moment of the proper time t of the ring's world tube, a whole line: of that t on the cosmic
            # chart and of its scale factor, ln a + a^2/2 = t, on the scale factor charts, for q = 2.
            height = t if key.split("/")[1] == "cosmic" else murphy_scale(t)
            return (lambda X: height), None
        if key == "ppn_metric/cartesian/axis":
            # The line through the body's centre: the equator at t = 0 on both sides of the body.
            lo, hi = self.reach(surface)
            return (lambda X: 0.0), [-hi, -lo, lo, hi]
        if key == "exponential_metric/cartesian/axis":
            # The line through r = 0: the equator at t = 0 on both halves, from r = m/3 to 6m.
            lo, hi = self.reach(surface)
            return (lambda X: 0.0), [-hi, -lo, lo, hi]
        if key == "exponential_metric/areal/radial":
            # The areal radius R = r e^(m/r) covers the near side, from the throat R = e m out.
            return (lambda X: 0.0), [math.e, self.reach(surface)[1] * math.exp(1 / self.reach(surface)[1])]
        if key == "exponential_metric/harmonic/radial":
            # u = 1/r of the isotropic radius the embedding reads.
            return (lambda X: 0.0), sorted(1 / x for x in self.reach(surface))
        if key == "ppn_metric/areal/radial":
            # The areal radius is the isotropic radius the embedding reads plus gamma m, at gamma = 1.
            return (lambda X: 0.0), [x + 1 for x in self.reach(surface)]
        if key.startswith("schrodinger_spacetime/"):
            # Schrodinger spacetime's plane of x and r, t = 0 and xi = 0, or T = 0 and V = 0: the event
            # (0, 0) of a plane of the time and the null coordinate at one depth.
            return (lambda X: 0.0), [0.0]
        if key.startswith("robinson_trautman/"):
            # The fronts of one retarded time differ only in size, so a moment is the whole
            # outgoing ray u = u_k, drawn against r and cu + r, from r = 0 to the box.
            return (lambda X: t + X), [0.0]
        if key.startswith("distorted_schwarzschild/"):
            # The black hole in a tidal field at m = 1: the horizon's bifurcation surface is the point t = 0
            # of every plane at x = 1, r = 2m, z = m on Weyl's axis and rho = 0 in his equatorial plane, and
            # the equatorial plane's t = 0 runs over the r the embedding reaches, with x = r/m - 1 and
            # rho = m sqrt(x^2 - 1).
            _, system, plane = key.split("/")
            if not mark["lines"]:
                return (lambda X: 0.0), [{"prolate_spheroidal": 1.0, "spherical": 2.0,
                                          "weyl": 1.0 if plane.startswith("axis") else 0.0}[system]]
            to = {"prolate_spheroidal": lambda r: r - 1, "spherical": lambda r: r,
                  "weyl": lambda r: math.sqrt(max((r - 1) ** 2 - 1, 0.0))}[system]
            return (lambda X: 0.0), [to(r) for r in self.reach(surface)]
        if key.startswith("string_black_hole/") and not mark["lines"]:
            # The horizon's bifurcation sphere, the point t = 0, r = r_s of the static planes.
            return (lambda X: 0.0), [1.0]
        if key.startswith("string_black_hole/") and "eddington" not in key:
            return (lambda X: 0.0), list(self.reach(surface, "static"))
        if key == "black_string/kerr_schild/radial":
            # The plane of the time and r is Schwarzschild's at r_s = 1: static t = 0 is cT = v - r = ln(r - 1).
            return (lambda X: math.log(X - 1)), None
        if key.startswith(("schwarzschild/eddington_finkelstein", "string_black_hole/eddington_finkelstein",
                           "black_string/eddington_finkelstein")):
            sign = 1 if "ingoing" in key else -1
            finkelstein = key.endswith("finkelstein")
            return (lambda X: sign * math.log(X - 1) + (0 if finkelstein else sign * X)), None
        if key == "kaluza_klein_black_hole/eddington_finkelstein_ingoing/finkelstein":
            # The static t = 0 at q = 2 r_s = 2: v = sqrt 2 (r + ln(r - 1)), drawn against v - sqrt 2 r.
            return (lambda X: math.sqrt(2) * math.log(abs(X - 1))), list(self.reach(surface))
        if key == "kaluza_klein_black_hole/einstein_eddington_finkelstein/finkelstein":
            # The static t = 0 of the Einstein metric at q = 2 r_s = 2: v = r_*, drawn against v - r, with
            # r_* = s + (3/2) ln(2s + 2r + 1) - sqrt 2 ln((3r + 1 + 2 sqrt 2 s)/|r - 1|), s = sqrt(r(r + 1)).
            def above(X):
                root = math.sqrt(X * (X + 1))
                return (root + 1.5 * math.log(2 * root + 2 * X + 1)
                        - math.sqrt(2) * math.log((3 * X + 1 + 2 * math.sqrt(2) * root) / abs(X - 1)) - X)
            return above, list(self.reach(surface))
        if key.startswith("dilaton_black_hole/eddington_finkelstein"):
            # The plane of t and r is Schwarzschild's at r_s = 1: static t = 0 is v = r + ln(r - 1) and
            # u = -r - ln(r - 1), drawn against v - r and u + r or against v and u.
            sign = 1 if "ingoing" in key else -1
            finkelstein = key.endswith("finkelstein")
            return (lambda X: sign * math.log(X - 1) + (0 if finkelstein else sign * X)), list(self.reach(surface))
        if key.startswith("btz/eddington_finkelstein"):
            # v = ct + r_* and u = ct - r_*, r_* = (1/2) ln|(r - 1)/(r + 1)|, drawn against v - r and u + r.
            sign = 1 if "ingoing" in key else -1
            return (lambda X: sign * (0.5 * math.log(abs((X - 1) / (X + 1))) - X)), list(self.reach(surface))
        if key.startswith(("kerr_de_sitter/kerr_", "kerr_de_sitter/kerr_schild")):
            # Static t = 0 is v = r_*, u = -r_* and c tau = s, drawn against v - r, u + r and c tau:
            # dr_*/dr = (r^2 + a^2)/Delta_r and ds/dr = r_s r/((1 - Lambda r^2/3) Delta_r), both zero at
            # r = 0, each the real part of a sum over the simple poles of its integrand, of residue times
            # ln(1 - r/pole). The poles are the two horizons r_- and r_+, the other two roots of Delta_r,
            # which are the roots of the quadratic left when those are divided out, and +-sqrt(3/Lambda).
            rs, a, L, rm, rp = ((2.0, 0.5, -3.0, 0.1368868045598, 0.8594395618858) if key.endswith("_ads") else
                                (1.0, 0.45, 0.2, 0.2787501957398, 0.7847845522526))
            c = -L / 3
            pair = cmath.sqrt((rm + rp) ** 2 - 4 * a * a / (c * rm * rp))
            roots = [rm, rp, (-(rm + rp) + pair) / 2, (-(rm + rp) - pair) / 2]
            cut = cmath.sqrt(3 / L)

            def delta(z):
                return (z * z + a * a) * (1 - L * z * z / 3) - rs * z

            def slope(z):
                return 2 * z * (1 - L * z * z / 3) - 2 * L * z * (z * z + a * a) / 3 - rs

            def rstar(r):
                return sum((z * z + a * a) / slope(z) * cmath.log(1 - r / z) for z in roots).real

            def shift(r):
                horizons = sum(rs * z / ((1 - L * z * z / 3) * slope(z)) * cmath.log(1 - r / z) for z in roots)
                return (horizons + sum(rs * z / (-2 * L * z / 3 * delta(z)) * cmath.log(1 - r / z) for z in (cut, -cut))).real
            if "kerr_schild" in key:
                return shift, list(self.reach(surface))
            sign = 1 if "ingoing" in key else -1
            return (lambda X: sign * (rstar(X) - X)), list(self.reach(surface))
        if key.startswith("kiselev/eddington_finkelstein"):
            # At w = -2/3, r_s = 1 and r_q = 8 the roots of f are 4 -+ 2 sqrt 2 and
            # r_* = (4 sqrt 2 - 4) ln|1 - r/r_-| - (4 sqrt 2 + 4) ln|1 - r/r_+|, which vanishes at r = 0.
            # Static t = 0 is v = r_* and u = -r_*, drawn against v - r and u + r or against v and u.
            sign = 1 if "ingoing" in key else -1
            finkelstein = key.endswith("finkelstein")
            inner, outer, root = 4 - 2 * math.sqrt(2), 4 + 2 * math.sqrt(2), 4 * math.sqrt(2)

            def rstar(r):
                return (root - 4) * math.log(abs(1 - r / inner)) - (root + 4) * math.log(abs(1 - r / outer))
            return (lambda X: sign * (rstar(X) - (X if finkelstein else 0))), list(self.reach(surface))
        if key == "kiselev/hyperbolic/radial":
            # The static t = 0 of the matter alone is eta = 0, from the centre chi = 0 out of the drawing.
            return (lambda X: 0.0), [0.0]
        if key == "kiselev/conformally_flat/radial":
            # There tau = cosh(chi) and rho = sinh(chi), the hyperbola tau = sqrt(1 + rho^2) from the centre.
            return (lambda X: math.sqrt(1 + X * X)), [0.0]
        if key.startswith("schwarzschild_de_sitter/eddington_finkelstein"):
            # At r_s = 1 and Lambda = 1/5 the roots of f are those of r^3 - 15r + 15, by Viete's
            # trigonometric solution, and r_* = sum_i ln|1 - r/r_i|/f'(r_i), which vanishes at r = 0.
            # Static t = 0 is v = r_* and u = -r_*, drawn against v - r and u + r or against v and u.
            sign = 1 if "ingoing" in key else -1
            finkelstein = key.endswith("finkelstein")
            angle = math.acos(-1.5 * math.sqrt(0.2)) / 3
            roots = [2 * math.sqrt(5) * math.cos(angle - 2 * math.pi * k / 3) for k in range(3)]

            def rstar(r):
                return sum(math.log(abs(1 - r / ri)) / (1 / ri ** 2 - 0.4 * ri / 3) for ri in roots)
            return (lambda X: sign * (rstar(X) - (X if finkelstein else 0))), list(self.reach(surface))
        if key.startswith("reissner_nordstrom_de_sitter/"):
            # The lukewarm hole: static t = 0 between r_+ and r_c and inside r_-, and the moment
            # H tau = 1 of the cosmological chart, on which the areal radius is rho + 1/2. In an
            # Eddington-Finkelstein chart the static moments are v = r_* and u = -r_*, and the
            # cosmological one is v = cT + r_* or u = cT - r_*, drawn against v - r and u + r or against
            # v and u; on the cosmological plane a static moment is tau of Lukewarm.tau_of.
            L, chart, which = Lukewarm, key.split("/")[1], mark["view"]
            lo, hi = self.reach(surface)
            if chart == "cosmological":
                if which == "cosmological":
                    return (lambda X: 1 / L.H), [lo, hi]
                at, sign = (0.0, 1) if which == "between" else (L.inside_t(), -1)
                ends = [(r - 0.5) / (L.H * L.tau_of(r, at, sign)) for r in ((lo,) if which == "inside" else ())]
                return (lambda X: L.tau_of(L.radius_at(X, at, sign, lo, hi), at, sign)), ends
            sign = {"static": 0, "eddington_finkelstein_ingoing": 1, "eddington_finkelstein_outgoing": -1}[chart]
            lean = sign if key.endswith("finkelstein") else 0
            if which == "cosmological":
                return (lambda X: L.static_t(1 / L.H, X) + sign * L.rstar(X) - lean * X), [lo + 0.5, hi + 0.5]
            return (lambda X: sign * L.rstar(X) - lean * X), [lo, hi]
        if key.startswith("black_saturn/"):
            # The plane of Black Saturn's ring is embedded in two pieces, outside the ring and about the
            # hole, and t = 0 is marked on each plane along its own piece; on the polar chart's
            # theta = pi/2 the outside piece is z = -r^2/2, from r = 0 out.
            ends = {p["id"]: [p["points"][0][0], p["points"][-1][0]] for p in surface["pieces"]}
            if key.endswith("/far"):
                return (lambda X: 0.0), [0.0, math.sqrt(-2 * ends["outside"][0])]
            return (lambda X: 0.0), ends[key.rsplit("/", 1)[1]]
        if key.startswith("myers_perry/ingoing_kerr"):
            # One spin in five dimensions at mu = 1 and a = 3/5, r_+ = 4/5: Boyer-Lindquist t = 0 is v = r_*,
            # r_* = r + (5/8) ln((r - 4/5)/(r + 4/5)), drawn against v - r, on either plane.
            return (lambda X: 0.625 * math.log(abs((X - 0.8) / (X + 0.8)))), list(self.reach(surface))
        if key.startswith("tangherlini/eddington_finkelstein"):
            # Five dimensions at r_h = 1: r_* = r + ln((r - 1)/(r + 1))/2, and the static t = 0 is v = r_*
            # and u = -r_*, drawn against v - r and u + r or against v and u.
            sign = 1 if "ingoing" in key else -1
            finkelstein = key.endswith("finkelstein")
            return (lambda X: sign * (0.5 * math.log(abs((X - 1) / (X + 1))) + (0 if finkelstein else X))), \
                list(self.reach(surface))
        if key.startswith("boulware_deser/eddington_finkelstein"):
            # At r_0 = 13/12 and l = 5/12 the horizon is r_h = 1, W = sqrt(r^4 + (65/72)^2) is 97/72 on it, and
            # 1/f = 1 + (97/72)/(r^2 - 1) - (W - r^2 + 25/72)/(2 (W + 97/72)), so
            # r_* = r + (97/144) ln((r - 1)/(r + 1)) - J(r)/2 with J the integral of the last fraction's
            # numerator over W + 97/72 from 0, by Simpson's rule. Static t = 0 is v = r_* and u = -r_*, drawn
            # against v - r and u + r or against v and u.
            sign = 1 if "ingoing" in key else -1
            finkelstein = key.endswith("finkelstein")

            def rstar(r, n=2000):
                def of(x):
                    w = math.sqrt(x ** 4 + (65 / 72) ** 2)
                    return (w - x * x + 25 / 72) / (w + 97 / 72)
                h = r / n
                J = h / 3 * (of(0) + of(r) + sum((4 if k % 2 else 2) * of(k * h) for k in range(1, n)))
                return r + 97 / 144 * math.log(abs((r - 1) / (r + 1))) - J / 2
            return (lambda X: sign * (rstar(X) - (X if finkelstein else 0))), list(self.reach(surface))
        if key.startswith(("quantum_oppenheimer_snyder/eddington_finkelstein", "quantum_oppenheimer_snyder/painleve")):
            # At m = 1 and alpha = 5/4 the static t = 0 is v = r_* and u = -r_*, with dr_*/dr = 1/f and
            # r_* = 0 at r_b, drawn against v - r and u + r or against v and u, and the Painleve-Gullstrand
            # T = r_* - L(r), with L the integral of 1/(1 + N) from r_b, N = sqrt(2/r - alpha/r^4).
            sign = -1 if "outgoing" in key else 1
            finkelstein = key.endswith("finkelstein")
            if "painleve" in key:
                return (lambda X: qos_rstar(X) - qos_lead(X)), list(self.reach(surface))
            return (lambda X: sign * (qos_rstar(X) - (X if finkelstein else 0))), list(self.reach(surface))
        if key == "eye_of_the_storm/kerr_ingoing/axis":
            # At m = 1, a = 3m/5 and ell = m/5 the moment t = 0 outside r_+ is v = r_* on the axis,
            # drawn against v - r, running off toward the horizon.
            return (lambda X: eos_rstar(X) - X), list(self.reach(surface))
        if key.startswith("bardeen/eddington_finkelstein"):
            # At r_s = 1 and g = 1/3 the static t = 0 is v = r_* and u = -r_*, with dr_*/dr = 1/f and
            # r_* = 0 at the centre, drawn against v - r and u + r or against v and u. It runs off
            # toward the horizon its view ends on, so a line ends at the box or at the embedding's reach.
            sign = 1 if "ingoing" in key else -1
            finkelstein = key.endswith("finkelstein")
            return (lambda X: sign * (bardeen_rstar(X) - (X if finkelstein else 0))), list(self.reach(surface))
        if key.startswith("born_infeld_charge/eddington_finkelstein"):
            # In units of r_0 at r_q = r_0/2 the static t = 0 is v = r_* and u = -r_*, with dr_*/dr = 1/f and
            # r_* = 0 at the centre, drawn against v - r and u + r or against v and u: the whole of it for
            # Hoffmann's particle, and for the black hole a curve that runs off toward the horizon.
            sign = 1 if "ingoing" in key else -1
            case = "particle" if key.endswith("/particle") else "hole"
            lean = 0 if key.endswith("/chart") else 1
            return (lambda X: sign * (born_infeld_rstar(X, case) - lean * X)), list(self.reach(surface))
        if key in ("ads_soliton/horowitz_myers/radial", "ads_soliton/poincare/tz"):
            # The soliton's t = 0, embedded in the proper distance rho from the tip, at r_0 = L = 1:
            # r = cosh^(2/3)(3 rho/2) on Horowitz and Myers's plane and z = 1/r on the Poincare plane.
            lo, hi = (math.cosh(1.5 * rho) ** (2 / 3) for rho in self.reach(surface))
            return (lambda X: 0.0), [lo, hi] if key.endswith("radial") else [1 / hi, 1 / lo]
        if key.startswith("topological_black_hole/") and not mark["lines"]:
            # The hyperbolic horizon's bifurcation surface, the point t = 0, r = r_h = L of the static planes.
            return (lambda X: 0.0), [1.0]
        if key == "topological_black_hole/brane/tz":
            # The black string's t = 0 on the brane's plane, where z = L^2/r.
            lo, hi = self.reach(surface)
            return (lambda X: 0.0), [1 / hi, 1 / lo]
        if key.startswith("topological_black_hole/eddington_finkelstein"):
            # At mu = 1 and L = 1, 1/f = r/((r - 1)(r^2 + r + 1)), and r_* = (1/3) ln|1 - r| - (1/6) ln(r^2 + r + 1)
            # + (arctan((2r + 1)/sqrt 3) - pi/6)/sqrt 3, which vanishes at r = 0. Static t = 0 is v = r_* and
            # u = -r_*, drawn against v - r and u + r.
            sign = 1 if "ingoing" in key else -1
            w = math.sqrt(3)

            def rstar(r):
                return (math.log(abs(1 - r)) / 3 - math.log(r * r + r + 1) / 6
                        + (math.atan((2 * r + 1) / w) - math.pi / 6) / w)
            return (lambda X: sign * (rstar(X) - X)), list(self.reach(surface))
        if key.startswith("schwarzschild_ads/eddington_finkelstein"):
            # At r_s = 2 and L = 1, 1/f = r/((r - 1)(r^2 + r + 2)), and r_* = (1/4) ln|1 - r| - (1/8) ln((r^2 + r +
            # 2)/2) + (5/(4 sqrt 7))(arctan((2r + 1)/sqrt 7) - arctan(1/sqrt 7)), which vanishes at r = 0.
            # Static t = 0 is v = r_* and u = -r_*, drawn against v - r and u + r or against v and u.
            sign = 1 if "ingoing" in key else -1
            finkelstein = key.endswith("finkelstein")
            w = math.sqrt(7)

            def rstar(r):
                return (0.25 * math.log(abs(1 - r)) - 0.125 * math.log((r * r + r + 2) / 2)
                        + 5 / (4 * w) * (math.atan((2 * r + 1) / w) - math.atan(1 / w)))
            return (lambda X: sign * (rstar(X) - (X if finkelstein else 0))), list(self.reach(surface))
        if key.startswith("quantum_btz/eddington_finkelstein"):
            # At l_3 = 1, M = F = 1/4 and l = 15/16, 1/H = r/(r^3 - r/4 - 15/64) with the root 3/4 and a
            # complex pair -3/8 +- i sqrt(11)/8, and r_* = Re sum_i A_i ln(1 - r/r_i), A_i = r_i/(3r_i^2 - 1/4), vanishes at r = 0.
            # Static t = 0 is v = r_* and u = -r_*, drawn against v - r and u + r or against v and u.
            sign = 1 if "ingoing" in key else -1
            finkelstein = key.endswith("finkelstein")
            roots = [0.75, complex(-0.375, math.sqrt(11) / 8), complex(-0.375, -math.sqrt(11) / 8)]

            def rstar(r):
                return sum(ri / (3 * ri * ri - 0.25) * cmath.log(1 - r / ri) for ri in roots).real
            return (lambda X: sign * (rstar(X) - (X if finkelstein else 0))), list(self.reach(surface))
        if key.startswith("reissner_nordstrom_ads/eddington_finkelstein"):
            # At L = 1, r_s = 27/8 and r_q^2 = 11/8, r^2 f = (r - 1)(r - 1/2)(r^2 + 3r/2 + 11/4), and
            # r_* = Re sum_i ln(1 - r/r_i)/f'(r_i) over the four roots, with f' = r_s/r^2 - 2r_q^2/r^3 + 2r, which
            # vanishes at r = 0. Static t = 0 is v = r_* and u = -r_*, outside r_+ and inside r_-, drawn against
            # v - r and u + r or against v and u.
            sign = 1 if "ingoing" in key else -1
            finkelstein = key.endswith("finkelstein")
            roots = (1, 0.5, complex(-0.75, math.sqrt(35) / 4), complex(-0.75, -math.sqrt(35) / 4))

            def rstar(r):
                return sum(cmath.log(1 - r / z) / (27 / 8 / z ** 2 - 11 / 4 / z ** 3 + 2 * z) for z in roots).real
            return (lambda X: sign * (rstar(X) - (X if finkelstein else 0))), list(self.reach(surface))
        if key.startswith("hayward/eddington_finkelstein"):
            # At m = 1 and ell = 12/(7 sqrt 7), 1/F = 1 + 2r^2/((r - 12/7)(r - 6/7)(r + 4/7)), and
            # r_* = r + 3 ln|1 - 7r/12| - (6/5) ln|1 - 7r/6| + (1/5) ln(1 + 7r/4), which vanishes at r = 0.
            # Static t = 0 is v = r_* and u = -r_*, outside r_+ and inside r_-, drawn against v - r and u + r
            # or against v and u.
            sign = 1 if "ingoing" in key else -1
            finkelstein = key.endswith("finkelstein")

            def rstar(r):
                return (r + 3 * math.log(abs(1 - 7 * r / 12)) - 1.2 * math.log(abs(1 - 7 * r / 6))
                        + 0.2 * math.log(1 + 7 * r / 4))
            return (lambda X: sign * (rstar(X) - (X if finkelstein else 0))), list(self.reach(surface))
        if key.startswith("regular_black_hole_pure_gravity/eddington_finkelstein"):
            # At r_0 = 1 and l = 12/25, 1/f = 1 + r^2/((r^2 - 16/25)(r^2 - 9/25)), and
            # r_* = r + (10/7) ln|(r - 4/5)/(r + 4/5)| - (15/14) ln|(r - 3/5)/(r + 3/5)|, which vanishes at r = 0.
            # Static t = 0 is v = r_* and u = -r_*, outside r_+ and inside r_-, drawn against v - r and u + r
            # or against v and u.
            sign = 1 if "ingoing" in key else -1
            finkelstein = key.endswith("finkelstein")

            def rstar(r):
                return (r + 10 / 7 * math.log(abs((r - 0.8) / (r + 0.8)))
                        - 15 / 14 * math.log(abs((r - 0.6) / (r + 0.6))))
            return (lambda X: sign * (rstar(X) - (X if finkelstein else 0))), list(self.reach(surface))
        if key in ("hayward/evaporating/history", "mass_inflation/ingoing/tail"):
            # A slice of constant v - r, level against v - r.
            return (lambda X: t), list(self.reach(surface))
        if key == "hiscock/ingoing/history":
            # A slice of constant v - r, level against v - r, over the pieces read in the ingoing chart.
            return (lambda X: t), list(self.reach(surface, "ingoing"))
        if key in ("hiscock/outgoing/history", "hiscock/flat/after"):
            # Outside the surface of pair creation the moment is the slice u + r = T', level against
            # u + r, that meets v - r = t on the surface: there r = 3N and u = v - 18N with
            # N = cos^2(pi (v - 2)/12), so T' = v - 15N at the root of v - 3N = t, and T' = t once the
            # hole is gone, t >= 8. The flat space after the hole holds its part inside the last ray
            # u_0 = 8, r < (t - 8)/2, where (u + v)/2 = t.
            if key == "hiscock/flat/after":
                return (lambda X: t), [0.0, (t - 8) / 2]
            mass = lambda v: math.cos(math.pi * (v - 2) / 12) ** 2
            meets = bisect(lambda v: v - 3 * mass(v) - t, 2, 8) if t < 8 else None
            outer = t if meets is None else meets - 15 * mass(meets)
            return (lambda X: outer), list(self.reach(surface, "outgoing"))
        if key.startswith("global_monopole/eddington_finkelstein"):
            # Letelier's black hole at Delta = 0.19 and r_s = 1: r_* = r/0.81 + ln|0.81 r - 1|/0.81^2,
            # and the static t = 0 is v = r_* and u = -r_*, drawn against v - r and u + r or against v and u.
            sign = 1 if "ingoing" in key else -1
            finkelstein = key.endswith("finkelstein")
            return (lambda X: sign * (X / 0.81 + math.log(abs(0.81 * X - 1)) / 0.81 ** 2 - (X if finkelstein else 0))), \
                list(self.reach(surface))
        if key.startswith("kantowski_sachs/"):
            # A dust moment eta is level in the dust chart and at ct = pi/2 + eta + sin(eta) cos(eta)
            # in the comoving one, b_0 = 1; a vacuum moment is level at its T, its cylinder's radius.
            if "/schwarzschild_interior/" in key:
                return (lambda X: surface["pieces"][0]["points"][0][1]), None
            return (lambda X: t if "/dust/" in key else math.pi / 2 + t + math.sin(t) * math.cos(t)), None
        if key == "domain_wall/planar/tz":
            # A moment kct of the global chart meets the plane x = y = 0 along t = const, every z.
            return (lambda X: t), None
        if key == "domain_wall/inertial/through":
            # In the inertial chart of the side z < 0 it is cT = |R| tanh(kct), out to the wall at
            # R = cosh(kct), k = 1.
            return (lambda X: abs(X) * math.tanh(t)), [-math.cosh(t), math.cosh(t)]
        if key == "de_sitter/flat_slicing/tx":
            return (lambda X: -0.5 * math.log(1 + X * X)), None
        if key == "pp_wave/exact_plane_wave/tz" or key.startswith("aichelburg_sexl/null_cartesian"):
            return (lambda X: X + t), None
        if key.startswith("light_beam/"):
            # Bonnor's wave front u = u_k is ct - z = sqrt(2) u_k, drawn against z and ct.
            return (lambda X: X + math.sqrt(2) * t), None
        if key.startswith("string_wave/"):
            # The event at u = t on the surface v = 0, drawn against (v - u)/2 and (u + v)/2. In the
            # moving string chart its V is 2(X - A)A' + int_0^t A'^2 on the line X, for the pulse
            # A = exp(-4u^2)/2, whose A'^2 = 16u^2 exp(-8u^2).
            V = 0.0
            if "/moving_string/" in key:
                line = {"behind": -0.25, "ahead": 0.75}[key.rsplit("/", 1)[1]]
                A, dA = math.exp(-4 * t * t) / 2, -4 * t * math.exp(-4 * t * t)
                V = (2 * (line - A) * dA + math.sqrt(math.pi / 8) * math.erf(2 * math.sqrt(2) * t) / 2
                     - t * math.exp(-8 * t * t))
            return (lambda X: X + t), [(V - t) / 2]
        if key == "oppenheimer_snyder/exterior_schwarzschild/radial":
            lo, hi = self.reach(surface, "comoving_synchronous")

            def ct(X):
                R = bisect(lambda R: novikov_t(R, t)[0] - X, lo, hi) if t else X
                return novikov_t(R, t)[1] if t else 0.0
            # The curve runs from the star's surface, the shell released at lo, outward.
            return ct, [novikov_t(lo, t)[0] if t else lo, novikov_t(hi, t)[0] if t else hi]
        if key == "lindquist_wheeler_lattice/schwarzschild_cell/radial":
            lo, hi = self.reach(surface, "lindquist_wheeler")
            lo += 1e-9

            def ct(X):
                # Novikov's slice of the cell: the shell whose radius is X at the proper time t, held to
                # the shells the embedding reaches.
                if not t:
                    return 0.0
                if X <= novikov_t(lo, t)[0]:
                    return novikov_t(lo, t)[1]
                if X >= novikov_t(hi, t)[0]:
                    return novikov_t(hi, t)[1]
                return novikov_t(bisect(lambda R: novikov_t(R, t)[0] - X, lo, hi), t)[1]
            return ct, [novikov_t(hi, t)[0] if t else hi, 1.0]
        if key == "lindquist_wheeler_lattice/lindquist_wheeler/shells":
            return (lambda X: t), list(self.reach(surface, "lindquist_wheeler"))
        if key == "lindquist_wheeler_lattice/comparison_hypersphere/radial":
            # Drawn in units of a_m = r_s/sin^3(psi) at the eight cells' psi.
            return (lambda X: t * math.sin(LW_PSI) ** 3), list(self.reach(surface, "comparison_hypersphere"))
        if key == "oppenheimer_snyder/interior_comoving/through":
            return (lambda X: t / (2 * math.sqrt(2))), [0, math.pi / 4]
        if key in ("semiclosed_world/comoving/dust", "semiclosed_world/conformal/dust"):
            # A moment of the dust's proper time, in units of a_m = 2 sqrt 2 r_s, or its conformal
            # time, sqrt 2 (eta + sin eta) = c tau, across the dust to chi_0 = 3 pi/4.
            eta = bisect(lambda e: math.sqrt(2) * (e + math.sin(e)) - t, 0, math.pi) if t else 0.0
            return (lambda X: eta if "conformal" in key else t / (2 * math.sqrt(2))), [0, 3 * math.pi / 4]
        if key == "datt_ruban_t_models/exterior_kruskal/kruskal":
            # Novikov's shells outside a T-sphere, from its surface s = 0, on which V = U, to the shell
            # the embedding's rim rests at, drawn with (V - U)/2 across and (U + V)/2 up.
            s_top = self.reach(surface, "comoving_synchronous")[1]
            across = lambda s: (novikov_uv(s, t)[1] - novikov_uv(s, t)[0]) / 2

            def up(X):
                s = bisect(lambda x: across(x) - X, 0.0, s_top) if 0 < X < across(s_top) else (0.0 if X <= 0 else s_top)
                return sum(novikov_uv(s, t)) / 2
            return up, [0.0, across(s_top)]
        if key in ("semiclosed_world/schwarzschild/radial", "semiclosed_world/isotropic/radial"):
            # Novikov's shells of the far sheet, from the throat, which rests at r_s, to the shell the
            # embedding's rim rests at, r_s (s^2 + 1); the isotropic radius of the areal radius r on
            # that sheet is (r - 1/2 + sqrt(r (r - 1)))/2, and at tau = 0 the isotropic plane holds the
            # whole moment, from the surface of the dust behind the throat.
            s_lo, s_hi = self.reach(surface, "comoving_synchronous")
            top = s_hi * s_hi + 1
            iso = "isotropic" in key
            to_x = (lambda r: (r - 0.5 + math.sqrt(r * (r - 1))) / 2) if iso else (lambda r: r)
            to_r = (lambda X: X * (1 + 1 / (4 * X)) ** 2) if iso else (lambda X: X)
            if not t:
                behind = s_lo * s_lo + 1
                return (lambda X: 0.0), [(behind - 0.5 - math.sqrt(behind * (behind - 1))) / 2 if iso else 1.0, to_x(top)]

            def ct(X):
                r = min(max(to_r(X), novikov_t(1 + 1e-9, t)[0]), novikov_t(top, t)[0])
                return novikov_t(bisect(lambda R: novikov_t(R, t)[0] - r, 1 + 1e-9, top), t)[1]
            return ct, [to_x(novikov_t(top, t)[0])]
        if key.startswith("white_hole/"):
            # The collapse with the time reversed: the dust's moment tau is the collapse's moment
            # pi a_m/2 - tau before the rest, a_m = 2 sqrt 2 r_s.
            rest = math.pi * math.sqrt(2)
            if key == "white_hole/interior_comoving/through":
                return (lambda X: t / (2 * math.sqrt(2))), [0, math.pi / 4]
            if key == "white_hole/interior_conformal/through":
                eta = bisect(lambda e: math.sqrt(2) * (e - math.sin(e)) - t, 0, 2 * math.pi)
                return (lambda X: eta), [0, math.pi / 4]
            lo, hi = self.reach(surface, "comoving_synchronous")
            back = max(rest - t, 0.0)

            def shell(X):
                return bisect(lambda R: novikov_t(R, back)[0] - X, lo, hi) if back else X
            ends = [novikov_t(lo, back)[0] if back else lo, novikov_t(hi, back)[0] if back else hi]
            if key == "white_hole/exterior_schwarzschild/radial":
                return (lambda X: -novikov_t(shell(X), back)[1] if back else 0.0), ends

            def u_plus_r(X):
                # Kruskal's V of the collapse's shell is e^(-u/2) sqrt(2) e^(1 + pi/2) of the white hole's.
                # A point rounded past the end of the line, where the curve is steep beside the
                # surface, is read along the tangent there.
                if not ends[0] <= X <= ends[1]:
                    edge = min(max(X, ends[0]), ends[1])
                    inward = 1e-5 if edge == ends[0] else -1e-5
                    return u_plus_r(edge) + (u_plus_r(edge + inward) - u_plus_r(edge)) / inward * (X - edge)
                R = shell(X)
                eta = bisect(lambda e: 0.5 * R * math.sqrt(R) * (e + math.sin(e)) - back, 0.0, math.pi)
                k = math.sqrt(R - 1)
                log_V = (math.log(k * math.cos(eta / 2) + math.sin(eta / 2))
                         + (X + k * (eta + 0.5 * R * (eta + math.sin(eta)))) / 2)
                return math.pi + 2 + math.log(2) - 2 * log_V + X
            return u_plus_r, ends
        if key.startswith(("einstein_rosen_waves/", "twisting_cylindrical_waves/")):
            # A moment ct = T, drawn against rho and ct in both charts, out to the embedding's reach.
            return (lambda X: t), list(self.reach(surface))
        if key == "senovilla/cylindrical/radial":
            # A moment act = T, drawn against a rho and act, out to the embedding's reach.
            return (lambda X: t), list(self.reach(surface))
        if key == "gowdy/areal/plane":
            # A moment of the areal time, drawn against theta and t, once round the torus.
            return (lambda X: t), list(self.reach(surface))
        if key == "gowdy/logarithmic/plane":
            # The same moment against theta and -tau = ln t.
            return (lambda X: math.log(t)), list(self.reach(surface))
        if key in ("vaidya/eddington_finkelstein_ingoing/shell", "bonnor_vaidya/eddington_finkelstein_ingoing/shell"):
            return (lambda X: t), list(self.reach(surface))
        if key.startswith("charged_shell/"):
            # The balanced shell at rest: t = 0 in the isotropic chart, from the shell out. The falling
            # shell: the slice v - r = t outside the shell, level against v - r, the curve
            # ct = t + r - r_* in the static chart, outside r_+ and inside r_-, and u + r = t + 2r - 2r_*
            # in the outgoing chart inside r_-; inside the shell, the moment of the flat time T at which
            # that slice meets the shell.
            if key == "charged_shell/exterior_isotropic/point":
                return (lambda X: 0.0), list(self.reach(surface, "exterior_isotropic"))
            if key == "charged_shell/exterior_ingoing/shell":
                return (lambda X: t), list(self.reach(surface, "exterior_ingoing"))
            if key.startswith("charged_shell/exterior/"):
                return (lambda X: t + X - charged_rstar(X)), list(self.reach(surface, "exterior_ingoing"))
            if key == "charged_shell/exterior_outgoing/shell":
                return (lambda X: t + 2 * X - 2 * charged_rstar(X)), list(self.reach(surface, "exterior_ingoing"))
            level = charged_inner_time(charged_on_slice(t))
            return (lambda X: level), list(self.reach(surface, "interior"))
        if key.startswith("penrose_impulsive_wave/"):
            # Penrose's wave at beta = 1/2: the moment t behind the front and T = t/beta ahead of it. On the
            # retarded equator, drawn against r and u + r, that is ct out to the front and
            # r + 2(ct - r)/(1 + beta^2) beyond it. A null plane of constant rho meets it in the one event
            # V = ct/sqrt(2), U = (rho^2 - 1) ct/sqrt(2), drawn at (U + V)/2 and (V - U)/2.
            b = 0.5
            if key.endswith("behind/radial"):
                return (lambda X: t), list(self.reach(surface, "behind"))
            if key.endswith("ahead/radial"):
                return (lambda X: t / b), list(self.reach(surface, "ahead"))
            if key.endswith("retarded/equator"):
                return (lambda X: t if X <= t else X + 2 * (t - X) / (1 + b * b)), [0.0, t]
            rho = {"unit": 1.0, "near": 0.5}[key.rsplit("/", 1)[1]]
            at = rho * rho * t / (2 * math.sqrt(2))
            return (lambda X: (2 - rho * rho) * t / (2 * math.sqrt(2)) + (X - at)), [at]
        if key.startswith("black_to_white_hole/"):
            # Haggard and Rovelli's shell: Schwarzschild's slice ct = t outside the shell, from the shell
            # out, and inside it the flat moment ct = -R at which that slice meets the shell. Drawn
            # against r and ct, Kruskal's (V - U)/2 and (U + V)/2, the Painleve-Gullstrand time
            # t + 2 sqrt(r) + ln((sqrt(r) - 1)/(sqrt(r) + 1)), and Lemaitre's rho - c tau = 2 r^(3/2)/3.
            if key == "black_to_white_hole/interior/bounce":
                lo, hi = self.reach(surface, "interior")
                return (lambda X: lo), [0.0, hi - lo]
            lo, hi = self.reach(surface, "schwarzschild")

            def pg(r):
                return t + 2 * math.sqrt(r) + math.log((math.sqrt(r) - 1) / (math.sqrt(r) + 1))
            if key == "black_to_white_hole/schwarzschild/radial":
                return (lambda X: t), [lo, hi]
            if key == "black_to_white_hole/painleve_gullstrand_ingoing/falling":
                return pg, [lo, hi]
            if key == "black_to_white_hole/lemaitre/falling":
                return (lambda X: pg((1.5 * X) ** (2 / 3))), [2 * lo ** 1.5 / 3, 2 * hi ** 1.5 / 3]

            def uv(r):
                a = math.sqrt(r - 1) * math.exp(r / 2)
                return -a * math.exp(-t / 2), a * math.exp(t / 2)
            across = lambda r: (uv(r)[1] - uv(r)[0]) / 2  # noqa: E731

            def up(X):
                r = bisect(lambda x: across(x) - X, lo, hi) if across(lo) < X < across(hi) else (lo if X <= across(lo) else hi)
                return sum(uv(r)) / 2
            return up, [across(lo), across(hi)]
        if key.startswith("israel_shell/"):
            # Israel's shell: the slice v - r = t outside the shell, level against v - r, and in
            # Schwarzschild's chart the curve ct = t - ln(r - 1), which leaves the drawing toward r_s;
            # inside the shell, the moment of the flat time T at which that slice meets the shell.
            if key == "israel_shell/exterior_ingoing/shell":
                return (lambda X: t), list(self.reach(surface, "exterior_ingoing"))
            if key == "israel_shell/exterior/radial":
                return (lambda X: t - math.log(X - 1)), list(self.reach(surface, "exterior_ingoing"))
            level = israel_inner_time(israel_on_slice(t))
            return (lambda X: level), list(self.reach(surface, "interior"))
        if key == "tippett_tsang/cartesian/tx":
            # The moment ct = A/2 over the whole width of the height's plane.
            u = next(p for p in surface["pieces"] if "grid" in p)["grid"]["u"]
            return (lambda X: 0.5), [u[0], u[-1]]
        if key == "krasnikov/cylindrical/tx":
            path = next(c for c in surface["curves"] if c["class"] == "path")
            grid = next(p for p in surface["pieces"] if "grid" in p)["grid"]
            return (lambda X: 5.0), [grid["u"][0] - path["points"][0][0], grid["u"][-1] - path["points"][0][0]]
        if key.startswith("holographic_superconductor/inversion/"):
            # A moment of the inversion is the line rho = -u, drawn at the height -rho, its time u.
            return (lambda X: t), None
        if key == "kasner_scalar/logarithmic/taux":
            # A moment of t against x and -tau = ln t.
            return (lambda X: math.log(t)), None
        if key.startswith("kasner_scalar/kaluza_klein/"):
            # The same moment in the time of five dimensions, t = T^k/k with k = 1 + s_5/2.
            k = 1 + 10 / (13 * math.sqrt(6) - 10)
            return (lambda X: (k * t) ** (1 / k)), None
        if key == "big_rip/conformal/radial":
            # A moment of t at w = -3/2 is the line eta = -eta_0 (-t/t_0)^(7/3), in units of eta_0.
            return (lambda X: -(-t) ** (7 / 3)), None
        if key.startswith(("kasner", "bianchi", "malament_hogarth", "big_rip")):
            if key.startswith("malament_hogarth"):
                lo, hi = self.reach(surface)
                return (lambda X: t), [-hi, -lo, lo, hi]
            return (lambda X: t), None
        if key == "misner/misner/plane":
            return (lambda X: -t * t / 4), None
        if key.startswith("khan_penrose/"):
            # The event on sigma = 0 at tau = t, u = v = sin(t/2): drawn against v - u and u + v, or
            # against sigma and tau.
            return (lambda X: 2 * math.sin(t / 2) if "/double_null/" in key else t), [0.0]
        if key.startswith("belinski_zakharov/"):
            # The event on xi = 0 at tau = t, which is z = 0 at ct = w sinh(t): drawn against xi and tau, or
            # against z and ct, at w = 1.
            return (lambda X: math.sinh(t) if "/canonical/" in key else t), [0.0]
        if key.startswith("bell_szekeres/"):
            # The event on eta = 0 at xi = t, with a = b = 1: u = v = t/2 drawn against v - u and u + v;
            # xi against eta; chi = t - pi/2 against rho; U = V = -cos t/(sqrt 2 (1 + sin t)) drawn
            # against V - U and U + V; and t = tan xi at r = sec xi.
            chart = key.split("/")[1]
            if chart == "bertotti_robinson":
                return (lambda X: math.tan(t)), [1 / math.cos(t)]
            height = {"double_null": t, "time_space": t, "global": t - math.pi / 2,
                      "kruskal_szekeres": -math.sqrt(2) * math.cos(t) / (1 + math.sin(t))}[chart]
            return (lambda X: height), [0.0]
        if key.startswith("chandrasekhar_xanthopoulos/"):
            # The event on lambda = 0 at psi = t, at m = 1 and a = 4/5: eta = sin t against mu; psi against
            # lambda; -r = 0.6 sin t - 1 against theta/pi = 1/2; and v - r = r_*(r) - r against r on the
            # equator of the ingoing chart, with r_* = 0 at r = m.
            chart = key.split("/")[1]
            r = 1 - 0.6 * math.sin(t)
            if chart == "kerr_ingoing":
                star = (r - 1) - (2 / 3) * math.log((r - 0.4) / 0.6) + (8 / 3) * math.log((1.6 - r) / 0.6)
                return (lambda X: star - r), [r]
            if chart == "boyer_lindquist":
                return (lambda X: -r), [0.5]
            return (lambda X: math.sin(t) if chart == "prolate" else t), [0.0]
        if key in ("misner/milne/plane", "gott_time_machine/grant_milne/plane"):
            return (lambda X: t), None
        if key.startswith("self_creating_universe/"):
            # Gott and Li's moment c tau = t of the Kantowski-Sachs chart at r_0 = 1: level on its own plane,
            # x = cosh(t) exp(-tau) on the steady state plane, and rho = -eta cosh(t) on the conformal one.
            chart = key.split("/")[1]
            if chart == "steady_state":
                return (lambda X: math.log(math.cosh(t) / X)), None
            if chart == "conformal":
                return (lambda X: -X / math.cosh(t)), None
            return (lambda X: t), None
        if key.startswith("ori_time_machine/"):
            # Ori's moment t is level on every cylinder of the time and z: at T = t on the central circle
            # in both charts, where the charts agree, at the foliation's own t at x = 4, and there at
            # T = t + (a/2 - e) x^2 = t - 3/2 in the vacuum core's chart. On the Brinkmann plane it is
            # the hyperbola uv = -2t with both negative, drawn against X = (v - u)/2 and Y = (u + v)/2,
            # where uv = Y^2 - X^2.
            if key.endswith("/brinkmann/plane"):
                return (lambda X: -math.sqrt(X * X - 2 * t)), None
            return (lambda X: t - 1.5 if key.endswith("/vacuum_core/off_centre") else t), None
        if key.startswith(("godel/cylindrical", "stockum_dust/cylindrical", "som_raychaudhuri/cylindrical",
                           "vuorio_warped_ads/cylindrical",
                           "maitra_dust/cylindrical",
                           "bonnor_rotating_dust/cylindrical/outside",
                           "minkowski/rindler")):
            return (lambda X: 0.0), None
        if key == "spinning_string/helical/outside":
            # One turn of the helix c tau = a phi~/b = r_c phi~, drawn against r phi~ at r = 3 r_c/2.
            return (lambda X: X / 1.5), None
        if key.startswith("spinning_string/"):
            return (lambda X: 0.0), None
        if key.startswith("c_metric"):
            # The equator's t = 0 along the axis, and the horizon's bifurcation sphere as a point,
            # at r = 2m, or y = 1/(alpha r) = 3 on the Hong-Teo plane.
            hong_teo = "/hong_teo/" in key
            if not mark["lines"]:
                return (lambda X: 0.0), [3.0 if hong_teo else 2.0]
            lo, hi = self.reach(surface, "spherical")
            return (lambda X: 0.0), ([6 / hi, 6 / lo] if hong_teo else [lo, hi])
        if key.startswith("ab_metrics/"):
            # A moment tau of BI's de Sitter slicing, from the neck r = b out to the embedding's reach: a
            # line of constant time along r, through the neck at rho = +-sqrt(1 - b/r) and round the
            # circle of phi; on the static surface cos(theta) = cosh(tau) sin(phi) with
            # tanh(tau_s) = tanh(tau)/cos(phi); and on the inertial plane T = |X| tanh(tau).
            lo, hi = self.reach(surface)
            view_key = key.partition("/")[2]
            if view_key in ("b1_cone/radial", "b1_static/radial"):
                return (lambda X: t), [lo, hi]
            if view_key == "b1_neck/through":
                return (lambda X: t), [-math.sqrt(1 - 1 / hi), math.sqrt(1 - 1 / hi)]
            if view_key == "b1_cone/de_sitter":
                return (lambda X: t), None
            if view_key == "b1_static/surface":
                def static_time(X):
                    cos_phi = math.sqrt(max(1 - (math.cos(X) / math.cosh(t)) ** 2, 1e-300))
                    return math.atanh(max(-1 + 1e-15, min(1 - 1e-15, math.tanh(t) / cos_phi)))
                return static_time, None
            if view_key == "b1_cartesian/TX":
                ends = [sign * r * math.cosh(t) for sign in (-1, 1) for r in (lo, hi)]
                return (lambda X: abs(X) * math.tanh(t)), ends
            raise KeyError(key)
        if key.startswith("nariai"):
            # A global moment runs round the whole circle of chi, and on the static patch it is
            # sinh(ct) = sinh(ct_k)/sqrt(1 - r^2), Lambda = 1, from horizon to horizon; the sphere is
            # the event t = 0, r = 0, which is chi = pi/2 of the global chart.
            static = "/static/" in key
            if not mark["lines"]:
                return (lambda X: 0.0), [0.0 if static else math.pi / 2]
            if static:
                return (lambda X: math.asinh(math.sinh(t) / math.sqrt(max(1 - X * X, 1e-300)))), None
            return (lambda X: t), None
        if key == "cremmer_scherk/cartesian/tx":
            # The cylinder's moment t = 0 along x from -a to a, and the sphere at the event x = a.
            lo, hi = self.reach(surface) if mark["lines"] else (1, 1)
            return (lambda X: 0.0), [lo, hi]
        if key == "cremmer_scherk/cartesian/circle":
            # The cylinder's ring x = 0 at t = 0, the whole circle of phi.
            return (lambda X: 0.0), [0.0, 2 * math.pi]
        if key.startswith("plebanski_hacyan/sphere"):
            # The equator's moment t = 0 along z from -b to b, of which the Rindler chart covers
            # 0 < chi < b on tau = 0, and the sphere at the event z = b, chi = b.
            lo, hi = self.reach(surface) if mark["lines"] else (1, 1)
            return (lambda X: 0.0), [lo, hi]
        if key.startswith("plebanski_hacyan/anti_nariai"):
            # The hyperbolic plane at the event tau = 0, chi = 1, which is r = a cosh(1) of the static chart.
            return (lambda X: 0.0), [math.cosh(1.0) if "_static/" in key else 1.0]
        if key.startswith("bertotti_robinson"):
            lo, hi = self.reach(surface) if mark["lines"] else (1, 1)
            return (lambda X: 0.0), [lo, hi]
        if key.startswith("near_horizon_extreme_kerr/"):
            # The moment tau = 0 of the global chart along y, and the horizon's sphere at y = 0. On the
            # Poincare planes it is t = 0 with r = r_0 (sqrt(1 + y^2) + y) and x = r_0^2/r.
            lo, hi = self.reach(surface, "global") if mark["lines"] else (0.0, 0.0)
            if "/global/" not in key:
                lo, hi = (math.sqrt(1 + y * y) + y for y in (lo, hi))
            if "/inverse_radius/" in key:
                lo, hi = 1 / hi, 1 / lo
            return (lambda X: 0.0), [lo, hi]
        if key == "einstein_static/static_areal/radial":
            # r = R sin chi carries the near hemisphere of the moment, out to the equator r = R.
            return (lambda X: 0.0), [0.0, 1.0]
        if key.startswith("coleman_de_luccia/"):
            # A moment c tau of the open universe inside the anti-de Sitter bubble, l = 1: level in its own
            # chart, at eta = ln tan(c tau/2) in the conformal time, and in the static chart inside the wall
            # the curve sqrt(1 + r^2) cos(ct) = cos(c tau), out to r = sin(c tau) sinh(chi) of the chi reached.
            hi = self.reach(surface)[1]
            if key == "coleman_de_luccia/static_inside/out_of_flat":
                return (lambda X: math.atan2(math.sqrt(math.sin(t) ** 2 + X * X), math.cos(t))), [0.0, math.sin(t) * math.sinh(hi)]
            if key == "coleman_de_luccia/open_conformal/negative":
                return (lambda X: math.log(math.tan(t / 2))), [0.0, hi]
            return (lambda X: t), [0.0, hi]
        if key.startswith("milne/"):
            # The moment ct: level in the comoving charts, at c tau = ln(ct) in the logarithmic one
            # and the hyperbola cT = sqrt(c^2t^2 + R^2) in the inertial one, out to the chi reached.
            hi = self.reach(surface)[1]
            if key == "milne/inertial/through":
                return (lambda X: math.sqrt(t * t + X * X)), [0.0, t * math.sinh(hi)]
            if key == "milne/logarithmic_time/radial":
                return (lambda X: math.log(t)), [0.0, hi]
            return (lambda X: t), [0.0, math.sinh(hi) if "spherical" in key else hi]
        if key.startswith("kastor_traschen/"):
            # One hole's moment is c tau = -8m/3. A moment of the midplane between two holes is level
            # in tau, at ct = ln(H tau)/H in the comoving chart, H = -3/32, across the midplane out to
            # the rho the embedding reaches, and meets the axis through the holes at z = 0.
            if "/isotropic/" in key:
                return (lambda X: -8 / 3), list(self.reach(surface))
            level = math.log(-3 / 32 * t) / (-3 / 32) if "/comoving/" in key else t
            lo, hi = self.reach(surface)
            ends = [0.0] if key.endswith("/tz") else [-hi, hi] if key.endswith("/tx") else [lo, hi]
            return (lambda X: level), ends
        if key.startswith("mcvittie/"):
            # A moment of cosmic time is level in both charts, from the throat out to the comoving r the
            # embedding reaches, which the areal chart draws at R = ar(1 + 1/(4ar))^2, r_s = 1, with
            # a = sinh^(2/3)(3 H_0 t/2) and H_0 = 1/sqrt(15).
            ends = list(self.reach(surface))
            if "/areal/" in key:
                a = math.sinh(1.5 * t / math.sqrt(15)) ** (2 / 3)
                ends = [a * r * (1 + 1 / (4 * a * r)) ** 2 for r in ends]
            return (lambda X: t), ends
        if key.startswith("sultana_dyer/"):
            # A moment of the conformal time eta is level in the Kerr-Schild plane, the curve
            # ct = eta - ln(r - 1) of Schwarzschild's time, which leaves the drawing toward the horizon,
            # and the line v = eta + r of the advanced time, r_s = 1, out to the r the embedding reaches.
            lo, hi = self.reach(surface)
            if "/schwarzschild_time/" in key:
                return (lambda X: t - math.log(X - 1)), [hi]
            if "/eddington_finkelstein_ingoing/" in key:
                return (lambda X: t + X), [lo, hi]
            return (lambda X: t), [lo, hi]
        if key == "bonnor_magnetic_dipole/spheroidal/strut":
            # The equatorial plane meets the axis between the holes at the one event theta = pi/2.
            return (lambda X: 0.0), [math.pi / 2]
        if key.startswith(("zipoy_voorhees/prolate_spheroidal/", "erez_rosen/prolate_spheroidal/")):
            # The prolate spheroidal x is r/m - 1 of the circles the embedding reaches in r, at m = 1.
            return (lambda X: 0.0), [r - 1 for r in self.reach(surface)]
        if key.startswith("point_particle_2plus1/"):
            # The cone is read in the conical chart from its apex, at alpha = 3/4 and ell = 1: the wedge
            # chart's r is the same, the circumference radius is alpha r and the isotropic radius
            # (alpha r)^(1/alpha); the planet's own plane carries the planet, in chi.
            system = key.split("/")[1]
            if system == "planet":
                return (lambda X: 0.0), list(self.reach(surface, "planet"))
            of = {"conical": lambda r: r, "wedge": lambda r: r, "circumference": lambda r: 0.75 * r,
                  "isotropic": lambda r: (0.75 * r) ** (4 / 3)}[system]
            return (lambda X: 0.0), [of(r) for r in self.reach(surface, "conical", reference=True)]
        if key == "levi_civita/kasner/radial":
            # The Kasner form's r is the proper distance from the axis, rho^Sigma/Sigma with
            # Sigma = 3/4 at sigma = 1/4, of the circles the embedding reaches in Weyl's rho.
            return (lambda X: 0.0), [rho ** 0.75 / 0.75 for rho in self.reach(surface)]
        if key == "mclenaghan_tariq_tupper/upper_half_plane/surface":
            # The strip is embedded in r, while this view draws x = a exp(r/a), at a = 1.
            return (lambda X: 0.0), [math.exp(r) for r in self.reach(surface)]
        if key in ("ghost_neutrinos/areal/radial", "ghost_neutrinos/null/null"):
            # Taub's z = (2/3) r^(3/2) at L = 1, the strip reached in z; the null view draws
            # (mu + nu)/2 = r^2/4 across, and t = 0 is mu = nu, its zero up.
            of = (lambda r: r) if key.endswith("radial") else (lambda r: r * r / 4)
            return (lambda X: 0.0), [of((1.5 * z) ** (2 / 3)) for z in self.reach(surface)]
        if key == "kaluza_klein_monopole/taub_nut/radial":
            # The Taub-NUT radius is rho = r + 2m, at m = 1, of the circles the cigar reaches in r.
            return (lambda X: 0.0), [r + 2 for r in self.reach(surface)]
        if key.startswith("brans_dicke_sphere/"):
            # The embedding is read in the harmonic chart at k = r_0 = 1, where e^(-2u) = 1 - 2 r_0/r: the
            # harmonic plane draws that u, the spherical plane the radii r reached, and the isotropic
            # plane, in units of B = r_0/2, the radius r - 1 + sqrt(r(r - 2)).
            lo, hi = self.reach(surface)
            if key == "brans_dicke_sphere/harmonic/radial":
                return (lambda X: 0.0), [lo, hi]
            of_r = {"spherical": lambda r: r,
                    "isotropic": lambda r: r - 1 + math.sqrt(r * (r - 2))}[key.split("/")[1]]
            return (lambda X: 0.0), [of_r(2 / (1 - math.exp(-2 * u))) for u in (hi, lo)]
        if key.startswith("static_five_soliton/"):
            # Lake's h labels the embedded circles; r_s=1 in every drawn chart.
            system = key.split("/")[1]
            of_h = {"spherical": lambda h: 1/(1-h*h),
                    "isotropic": lambda h: (1+h)/(4*(1-h)),
                    "lake": lambda h: h}[system]
            return (lambda X: 0.0), [of_h(h) for h in self.reach(surface, "lake")]
        if key.startswith("fisher_jnw/"):
            # The embedding is read in the harmonic chart at k = 1/2 and b = 1, where e^(-u) = 1 - b/r: the
            # harmonic plane draws ku at k = 1, half of that u, and the other three the radii r reached, as
            # r itself, as R = r - 3b/4, and as the isotropic radius (r - 1/2 + sqrt(r(r - 1)))/2.
            lo, hi = self.reach(surface)
            if key == "fisher_jnw/harmonic/radial":
                return (lambda X: 0.0), [lo / 2, hi / 2]
            of_r = {"spherical": lambda r: r, "jnw": lambda r: r - 0.75,
                    "isotropic": lambda r: (r - 0.5 + math.sqrt(r * (r - 1))) / 2}[key.split("/")[1]]
            return (lambda X: 0.0), [of_r(1 / (1 - math.exp(-u))) for u in (hi, lo)]
        if key.startswith("penney_scalar/"):
            # The embedding is read in the spherical chart at a = 1, b = 1/4 and gamma = 1/2: the spherical plane
            # draws the radii r reached, and the harmonic plane, at k = 1 with a = 8k/3, the u with
            # e^(-2u) = (r - a)/(r - b), a ratio the scale leaves alone.
            lo, hi = self.reach(surface)
            if key == "penney_scalar/harmonic/radial":
                return (lambda X: 0.0), sorted(-math.log((r - 1) / (r - 0.25)) / 2 for r in (lo, hi))
            return (lambda X: 0.0), [lo, hi]
        if key.startswith("jackiw_teitelboim_black_hole/"):
            # The Euclidean disc is read in the proper distance rho at L = r_h = 1, and its meridians are the
            # moment t = 0: level in every chart, from the horizon r = cosh(rho) = 1 out in the static chart,
            # across the bifurcation point to V - U = +-2 tanh(rho/2) on the Kruskal plane, drawn against
            # (V - U)/2, to sigma = pi/2 +- arctan(sinh(rho)) in the global chart, and to z = 2 e^(-+rho) on
            # the Poincare plane.
            lo, hi = self.reach(surface)
            return (lambda X: 0.0), {
                "static": [math.cosh(lo), math.cosh(hi)], "proper_distance": [lo, hi],
                "kruskal": [-math.tanh(hi / 2), math.tanh(hi / 2)],
                "global": [math.pi / 2 - math.atan(math.sinh(hi)), math.pi / 2 + math.atan(math.sinh(hi))],
                "poincare": [2 * math.exp(-hi), 2 * math.exp(hi)]}[key.split("/")[1]]
        if key.startswith("eternal_traversable_wormhole/"):
            # The Euclidean cylinder is read in the global sigma at L = 1, and its meridians are the moment
            # t = 0: level in every chart, across the strip in the global chart, at z = tan(sigma/2) on the
            # Poincare plane, and in the thermal chart on the half sigma <= pi/2 alone, from its edge rho = 0
            # at the middle out to cosh(rho) = 1/sin(sigma).
            lo, hi = self.reach(surface)
            return (lambda X: 0.0), {
                "global": [lo, hi], "poincare": [math.tan(lo / 2), math.tan(hi / 2)],
                "thermal": [0.0, math.acosh(1 / math.sin(lo))]}[key.split("/")[1]]
        if key.startswith("witten_black_hole/"):
            # The cigar is read in Witten's proper distance r at lambda = m = 1, and its meridian theta = 0 is
            # the moment t = 0 from the horizon out, where e^(2x) = w = cosh^2 r and e^sigma = sinh r: level in
            # the charts of t; v - x = sigma - x and u + x = x - sigma, sigma = ln(e^(2x) - 1)/2, in the
            # Eddington-Finkelstein charts; and the line U + V = 0 of the Kruskal plane, drawn against
            # (V - U)/2, on both sides of the bifurcation point, the other side being the meridian theta = pi.
            lo, hi = self.reach(surface)
            chart = key.split("/")[1]
            if chart.startswith("eddington"):
                sign = 1 if chart.endswith("ingoing") else -1
                return (lambda X: sign * (0.5 * math.log(max(math.expm1(2 * X), 1e-300)) - X)), [math.log(math.cosh(hi))]
            return (lambda X: 0.0), {
                "witten": [lo, hi], "schwarzschild_gauge": [math.log(math.cosh(lo)), math.log(math.cosh(hi))],
                "dilaton": [math.cosh(lo) ** 2, math.cosh(hi) ** 2], "conformal": [math.log(math.sinh(hi))],
                "kruskal": [-math.sinh(hi), math.sinh(hi)]}[chart]
        if key.startswith("roberts/"):
            # A moment t of Roberts's diagonal chart, over the rho the embedding reaches, in each chart's
            # drawn axes, with u = s (t - rho), v = (t + rho)/s and s = sqrt(1 + p): level at t against
            # rho, level at s t against r = s rho, the straight line of u and v against (v - u)/2 and
            # (u + v)/2, the curve (1 + p/2) v - R against R = sqrt(rho (rho - p t)), and tau + x against
            # x in the scaling chart, where e^(2x) - 1 = 2 (t + rho)/((1 + p)(rho - t)).
            chart, case = key.split("/")[1:]
            p = ROBERTS_P[case]
            s = math.sqrt(1 + p)
            lo, hi = self.reach(surface)
            if chart == "diagonal":
                return (lambda X: t), [lo, hi]
            if chart == "advanced":
                return (lambda X: s * t), [s * lo, s * hi]
            if chart == "double_null":
                X_of = lambda rho: (t * (1 / s - s) + rho * (1 / s + s)) / 2
                rho_of = lambda X: (2 * X - t * (1 / s - s)) / (1 / s + s)
                return (lambda X: (t * (1 / s + s) + rho_of(X) * (1 / s - s)) / 2), [X_of(lo), X_of(hi)]
            if chart == "areal":
                rho_of = lambda X: (p * t + math.sqrt(p * p * t * t + 4 * X * X)) / 2
                R_of = lambda rho: math.sqrt(max(rho * (rho - p * t), 0.0))
                return (lambda X: (1 + p / 2) * (t + rho_of(X)) / s - X), [R_of(lo), R_of(hi)]

            def rho_of(X):
                E = math.expm1(2 * X) * (1 + p)
                return t * (E + 2) / (E - 2)
            x_of = lambda rho: math.log1p(2 * (t + rho) / ((1 + p) * (rho - t))) / 2 if rho > t else math.inf
            return (lambda X: X - math.log(s * (rho_of(X) - t) / 2)), [x_of(lo), x_of(hi)]
        if key in ("damour_solodukhin/isotropic/radial", "einstein_rosen_bridge/isotropic/radial"):
            # The isotropic radius of the circles the embedding reaches in the areal radius R, at r_s = 1:
            # r = (R - 1/2 + sqrt(R(R - 1)))/2 on one side and 1/(16 r) on the other.
            R = self.reach(surface)[1]
            r = (R - 0.5 + math.sqrt(R * (R - 1))) / 2
            return (lambda X: 0.0), [1 / (16 * r), r]
        if key.startswith("simpson_visser/"):
            # The moment t = 0 of one of Simpson and Visser's geometries at r_s = 1: level in the two
            # charts of t, v - r = r_* - r in the ingoing chart and u + r = r - r_* in the outgoing one. The
            # black bounce's exterior is embedded in the areal radius rho, the others in r, and
            # rho^2 = r^2 + a^2; a moment through the throat reaches rho = a.
            _, chart, case = key.split("/")
            a = SV_A[case]
            lo, hi = self.reach(surface)
            if case == "bounce":
                lo, hi = (math.sqrt(rho * rho - a * a) for rho in (lo, hi))
            if chart == "areal":
                near = 0.0 if lo < 0 < hi else min(abs(lo), abs(hi))
                return (lambda X: 0.0), [math.sqrt(near * near + a * a), math.sqrt(max(lo * lo, hi * hi) + a * a)]
            if chart == "spherical":
                return (lambda X: 0.0), [lo, hi]
            sign = 1 if chart.endswith("ingoing") else -1
            return (lambda X: sign * (sv_rstar(X, a) - X)), [lo, hi]
        if key.startswith("bartnik_mckinnon/") and "/areal/" not in key:
            # The sphere r = 8 ell the embedding reaches, on the soliton with one zero: its isotropic
            # radius, its tortoise coordinate and its tau = ln(rho/ell), which bartnik_mckinnon.py
            # integrates and test_bartnik_mckinnon.py holds to these numbers.
            return (lambda X: 0.0), [0.0, BARTNIK_MCKINNON_REACH[key.split("/")[1]]]
        if key == "gravastar/interior_tortoise/radial":
            # The tortoise coordinate of the circles the embedding reaches inside the shell, at L = 2 and
            # C = 64/195: x = (L/sqrt(C)) artanh(r/L).
            return (lambda X: 0.0), [2 / math.sqrt(64 / 195) * math.atanh(r / 2) for r in self.reach(surface, "interior")]
        if key == "randall_sundrum/conformal/tw":
            # The conformal distance of the circles the embedding reaches in the proper distance y, at k = 1.
            return (lambda X: 0.0), [math.copysign(math.expm1(abs(y)), y) for y in self.reach(surface)]
        if key == "randall_sundrum/poincare/tz":
            # One side, from the wall at z = 1/k to z = e^{ky}/k of the farthest circle.
            return (lambda X: 0.0), [1.0, math.exp(self.reach(surface)[1])]
        if key.startswith("lifshitz_spacetime/"):
            # Lifshitz spacetime's moment t = 0 at z = 2 and L = 1, over the embedding's reach in the proper
            # distance rho: r = e^rho, u = e^-rho, w = e^(-2 rho)/2 and s = e^(2 rho); in the two null charts
            # it is the curve v = -w, drawn against v - r and v - s/2; and on a plane of t and x the line
            # t = 0 over the strip 0 <= x < 2 pi L.
            lo, hi = self.reach(surface)
            view_id = key.rsplit("/", 1)[1]
            if view_id.startswith("tx_"):
                return (lambda X: 0.0), [0.0, 2 * math.pi]
            return {
                "tr": ((lambda X: 0.0), [math.exp(lo), math.exp(hi)]),
                "tu": ((lambda X: 0.0), [math.exp(-hi), math.exp(-lo)]),
                "trho": ((lambda X: 0.0), [lo, hi]),
                "tw": ((lambda X: 0.0), [0.5 * math.exp(-2 * hi), 0.5 * math.exp(-2 * lo)]),
                "vr": ((lambda X: -0.5 / max(X, 1e-9) ** 2 - X), [math.exp(lo), math.exp(hi)]),
                "vs": ((lambda X: -0.5 / max(X, 1e-9) - X / 2), [math.exp(2 * lo), math.exp(2 * hi)]),
            }[view_id]
        if key.startswith("elliptic_de_sitter/"):
            # Elliptic de Sitter space's moment t = 0 of the global chart at l = 1, over the embedding's
            # reach in chi: the line eta = 0; U = -V with V = tan(pi/4 - chi/2), drawn against (V - U)/2 = V
            # at (U + V)/2 = 0; the static t = 0 with r = sin chi; and t = -ln(1 + x^2)/2 with x = tan chi
            # in the planar chart, which leaves the box before the rim.
            lo, hi = self.reach(surface)
            return {
                "global": ((lambda X: 0.0), [lo, hi]),
                "conformal": ((lambda X: 0.0), [lo, hi]),
                "kruskal": ((lambda X: 0.0), [math.tan(math.pi / 4 - hi / 2), math.tan(math.pi / 4 - lo / 2)]),
                "static": ((lambda X: 0.0), [math.sin(lo), math.sin(hi)]),
                "planar": ((lambda X: -0.5 * math.log1p(X * X)), None),
            }[key.split("/")[1]]
        if key.startswith("rp3_geon/"):
            # The RP3 geon's moment t = 0 at r_s = 1, over the embedding's reach in the areal radius r: the
            # line T = 0 of Kruskal's chart, with X = sqrt(r - 1) e^(r/2), and the line t = 0 of the other
            # two, the isotropic radius being (r - 1/2 + sqrt(r(r - 1)))/2.
            lo, hi = self.reach(surface)
            return (lambda X: 0.0), {
                "kruskal": [math.sqrt(lo - 1) * math.exp(lo / 2), math.sqrt(hi - 1) * math.exp(hi / 2)],
                "schwarzschild": [lo, hi],
                "isotropic": [(r - 0.5 + math.sqrt(r * (r - 1))) / 2 for r in (lo, hi)],
            }[key.split("/")[1]]
        if key == "bonnor_charged_dust/spheroid_interior/axis":
            # The spheroid's equatorial plane meets its axis at the centre of the disc, the event t = 0, u = 0.
            return (lambda X: 0.0), [0.0]
        if key == "bonnor_charged_dust/exterior_areal/radial":
            # The areal radius is the isotropic radius the embedding reads plus m, at m = 1.
            return (lambda X: 0.0), [x + 1 for x in self.reach(surface, "exterior")]
        if key == "btz_multi_holes_wormholes/stereographic/fold":
            # The moment of time symmetry is U = 0, the hyperbola c tau = -sqrt(4 + x^2) at l = 1, out to
            # x = 4 rho/(1 - rho^2) at the embedding's reach in rho.
            hi = self.reach(surface)[1]
            return (lambda X: -math.sqrt(4 + X * X)), [4 * hi / (1 - hi * hi)]
        if key == "btz_multi_holes_wormholes/exterior/radial":
            # t = 0 from the wormhole's horizon r_+ = 2 arccosh(2)/pi to the edge of the sheet on phi = 0.
            hi = self.reach(surface)[1]
            horizon = 2 * math.acosh(2) / math.pi
            return (lambda X: 0.0), [horizon, horizon * (1 + hi * hi) / (1 - hi * hi)]
        if key in ("three_brane_throat/areal/radial", "three_brane_throat/horizon/radial"):
            # The embedding reads the isotropic radius rho at L = 1: the areal radius is (rho^4 + 1)^(1/4)
            # and Gibbons, Horowitz and Townsend's w is rho over it.
            areal = [(x ** 4 + 1) ** 0.25 for x in self.reach(surface, "isotropic")]
            if "/areal/" in key:
                return (lambda X: 0.0), areal
            return (lambda X: 0.0), [x / r for x, r in zip(self.reach(surface, "isotropic"), areal)]
        if key in ("freund_rubin/conformal/radial", "freund_rubin/static/radial"):
            # The embedding reads the global rho at L = 1: the static chart's r is sinh(rho) and the
            # conformal chart's chi is arctan(sinh(rho)).
            areal = [math.sinh(x) for x in self.reach(surface, "global")]
            if "/static/" in key:
                return (lambda X: 0.0), areal
            return (lambda X: 0.0), [math.atan(r) for r in areal]
        if key == "freund_rubin/poincare/radial":
            # The cylinder is read in the proper distance sigma, and r = L e^(sigma/L).
            return (lambda X: 0.0), [math.exp(x) for x in self.reach(surface, "proper")]
        if key == "three_brane_throat/throat/radial":
            # The throat's cylinder is read in the proper distance sigma, and r = L e^(sigma/L).
            return (lambda X: 0.0), [math.exp(x) for x in self.reach(surface, "throat_proper")]
        if key == "kodama_kink/areal/radial":
            # The embedding is read along the proper distance R, and each of its points carries the areal
            # radius rho = r(R): the side R > 0 runs from the throat, the least rho, out to the rim.
            rho = [p[1] for piece in surface["pieces"] if "points" in piece and not piece.get("reference")
                   for p in piece["points"]]
            return (lambda X: 0.0), [min(rho), max(rho)]
        if key == "einstein_dirac_maxwell_wormhole/compact/radial":
            # The compact x of the circles the embedding reaches in the areal radius r, at r_0 = 1:
            # x = +-sqrt(1 - 1/r), one on each side of the throat.
            x = math.sqrt(1 - 1 / self.reach(surface)[1])
            return (lambda X: 0.0), [-x, x]
        if key == "btz_shock_wave/kruskal/plane":
            # The moment u + v = -1/2 at alpha = 1, drawn against (v - u)/2 and (u + v)/2: level at -1/4,
            # over the embedding's reach in u, where (v - u)/2 = -1/4 - u.
            lo, hi = self.reach(surface)
            return (lambda X: -0.25), [-0.25 - hi, -0.25 - lo]
        if key == "btz_shock_wave/exterior/radial":
            # Outside the right horizon the moment has u < -1/2 and v = -1/2 - u, so uv = (1 - r)/(1 + r) at
            # l = R = 1 gives u, and ct = ln(v/(-u))/2; it reaches r = 3 at u = -1 and runs off down the horizon.
            def t_of(r):
                k = (r - 1) / (r + 1)
                u = (-0.5 - math.sqrt(0.25 + 4 * k)) / 2
                return 0.5 * math.log((-0.5 - u) / -u)
            lo = self.reach(surface)[0]
            return t_of, [(1 + lo * (lo + 0.5)) / (1 - lo * (lo + 0.5))]
        if key == "anti_de_sitter/poincare/tx":
            hi = self.reach(surface)[1]
            x = math.sqrt(2 * (math.sqrt(1 + hi * hi) - 1))
            return (lambda X: 0.0), [-x, x]
        if key == "godel/cartesian/tx":
            hi = self.reach(surface)[1]
            return (lambda X: 0.0), [-2 * hi, 2 * hi]
        if key == "vuorio_warped_ads/disc/plane":
            # Poincare's disc reads the cylindrical chart's proper distance r as R = tanh(r/2).
            hi = math.tanh(self.reach(surface)[1] / 2)
            return (lambda X: 0.0), [-hi, hi]
        if key == "vuorio_warped_ads/fibred/plane":
            hi = self.reach(surface)[1]
            return (lambda X: 0.0), [-hi, hi]
        if surface["pieces"][0].get("grid"):
            u = surface["pieces"][0]["grid"]["u"][-1]
            return (lambda X: 0.0), [-u, u]
        lo, hi = self.reach(surface)
        through = key.endswith("/tx") and not view["mirror"]
        return (lambda X: 0.0), ([-hi, -lo, lo, hi] if through else [lo, hi])

    def check_shock_moment(self, key, view, surface, mark):
        """The moment u + v = -alpha/2 of the BTZ shock in the discontinuous chart, V = v + alpha Theta(u), at
        alpha = 1: U + V = -1/2 ahead of the shock and 1/2 behind it, the two halves overlapping in (V - U)/2,
        so each line is held to its own half, from the shock U = 0 to the embedding's reach in u."""
        X0, X1, Y0, Y1 = view["box"]
        lo, hi = self.reach(surface)
        checked = 0
        self.assertEqual(len(mark["lines"]), 2, key)
        for line in mark["lines"]:
            points = [(X0 + u[0] * (X1 - X0), Y0 + u[1] * (Y1 - Y0)) for u in line]
            U = [Y - X for X, Y in points]
            behind = sum(U) > 0
            for (X, Y), u in zip(points, U):
                self.assertAlmostEqual(u + (X + Y) - behind, -0.5, delta=1e-3, msg=f"{key} at {X, Y}")
                checked += 1
            self.assertEqual(sorted(round(e, 3) for e in (U[0], U[-1])), sorted([0.0, round(hi if behind else lo, 3)]))
        return checked

    def test_every_slice_on_a_spacetime_diagram_lies_on_its_moment_and_ends_where_the_embedding_does(self):
        """Each point of a line or a point of a slice, carried back through the view's box and
        axes, lies on its moment, and a line stops at the embedding's reach or at the box."""
        checked = 0
        for metric_id, data in self.diagrams.items():
            for system_id, views in data["systems"].items():
                for view in views:
                    key = f"{metric_id}/{system_id}/{view['id']}"
                    X0, X1, Y0, Y1 = view["box"]
                    for mark in view.get("slices", []):
                        _, surface = self.moment(metric_id, mark)
                        if mark["fills"]:
                            self.check_region_from_above(key, view, surface, mark)
                            continue
                        if key.startswith("simpson_visser/") and mark["view"] == "inside":
                            checked += self.check_bounce_moment(key, view, surface, mark)
                            continue
                        if key == "btz_shock_wave/discontinuous/plane":
                            checked += self.check_shock_moment(key, view, surface, mark)
                            continue
                        if key == "tippett_tsang/polar/strip":
                            checked += self.check_polar_moment(key, view, surface, mark)
                            continue
                        if key.startswith("brane_wall_universe/"):
                            checked += self.check_wall_moment(key, view, surface, mark)
                            continue
                        Y_of, ends = self.flat_expected(key, view, surface, mark)
                        for line in mark["lines"] + [[p] for p in mark["points"]]:
                            for u in line:
                                X, Y = X0 + u[0] * (X1 - X0), Y0 + u[1] * (Y1 - Y0)
                                h = 2e-4 * (X1 - X0)
                                slope = (abs(Y_of(min(X + h, X1)) - Y_of(max(X - h, X0 + 1e-9 if key.startswith(
                                    ("schwarzschild/edd", "string_black_hole/edd", "oppenheimer_snyder/ext",
                                     "semiclosed_world/schwarzschild", "semiclosed_world/isotropic",
                                     "lindquist_wheeler_lattice/schwarzschild_cell",
                                     "black_string/edd", "black_string/kerr", "kaluza_klein_black_hole/edd",
                                     "kaluza_klein_black_hole/einstein_edd")) else X0))) / (2 * h)
                                         if X0 < X < X1 else 0)
                                tol = 1e-4 * (Y1 - Y0) + slope * 1e-4 * (X1 - X0) + 1e-9
                                self.assertLess(abs(Y - Y_of(X)), 3 * tol, f"{key} {mark['label']} at {u}")
                                checked += 1
                            for u in (line[0], line[-1]):
                                edge = min(u[0], 1 - u[0], u[1], 1 - u[1]) < 1.5e-4
                                X = X0 + u[0] * (X1 - X0)
                                at_reach = ends is not None and any(abs(X - e) < 1.5e-4 * (X1 - X0) for e in ends)
                                self.assertTrue(edge or at_reach, f"{key} {mark['label']} stops at {u}")
        counted(self, checked, 200, "slices on a spacetime diagram")

    def sv_inside_r(self, surface):
        """The r of a moment of the black bounce between its horizons, a = r_s/2: its cylinder's radius is
        sqrt(r^2 + a^2), and r is positive before the proper time of r = 0, the third moment's."""
        radius = surface["pieces"][0]["points"][0][1]
        inside = next(v for v in self.embedding["simpson_visser"]["views"] if v["id"] == "inside")
        turn = inside["surfaces"][2]["time"]
        return math.copysign(math.sqrt(max(radius * radius - 0.25, 0.0)), turn - surface["time"])

    def check_polar_moment(self, key, view, surface, mark):
        """Tippett and Tsang's moment ct = A/2 in the polar chart, xi sin(lambda) = A/2: one line, which
        crosses each xi twice and stops at both ends where the height's plane does, at |x| = 1.6 A,
        inside the strip, which runs to xi = 3.5 A so that its plotted region is no taller than 1:2."""
        X0, X1, Y0, Y1 = view["box"]
        (line,) = mark["lines"]
        for u in line:
            xi, lam = X0 + u[0] * (X1 - X0), Y0 + u[1] * (Y1 - Y0)
            self.assertAlmostEqual(xi * math.sin(lam), 0.5, delta=1e-3, msg=f"{key} {mark['label']} at {u}")
        grid = next(p for p in surface["pieces"] if "grid" in p)["grid"]
        reach = math.hypot(grid["u"][-1], 0.5)
        self.assertLess(reach, X1)
        for u in (line[0], line[-1]):
            self.assertAlmostEqual(X0 + u[0] * (X1 - X0), reach, delta=1.5e-4 * (X1 - X0),
                                   msg=f"{key} {mark['label']} stops at {u}")
        return len(line)

    @staticmethod
    def wall_moment():
        """The wall's own moment t = l/c of the braneworld cosmology, at l = mu = 1 and k = 0, as
        brane_wall_universe.Moment follows it along the normal geodesics from the wall, by the
        distance w along it."""
        sys.path.insert(0, str(build.ROOT / "_tools" / "derivations"))
        import brane_wall_universe
        return brane_wall_universe

    def check_wall_moment(self, key, view, surface, mark):
        """The wall's own moment on its two flat views. On the Gaussian normal plane it is the line
        ct = l from the wall out to w_* = ln(5)/2, where the lapse vanishes and the view is hatched.
        On the static plane every point lies on the curve (t(w), r(w)) the moment draws through the
        exterior and the black hole, and each line ends at the wall or on the box."""
        import numpy as np
        bwu = self.wall_moment()
        X0, X1, Y0, Y1 = view["box"]
        lo, hi = self.reach(surface)
        if key.split("/")[1] == "gaussian_normal":
            (line,) = mark["lines"]
            for u in line:
                self.assertAlmostEqual(Y0 + u[1] * (Y1 - Y0), 1.0, delta=1.5e-4 * (Y1 - Y0), msg=key)
            ends = sorted(X0 + u[0] * (X1 - X0) for u in (line[0], line[-1]))
            for got, want in zip(ends, (lo, min(hi, float(bwu.focal_distance(1.0))))):
                self.assertAlmostEqual(got, want, delta=1.5e-4 * (X1 - X0), msg=key)
            return len(line)
        moment = bwu.Moment(1.0)
        across = 0.5 * math.acosh((2 + moment.E ** 2) / moment.root)
        w_in, w_out = moment.w0 - across, moment.w0 + across
        crowd = 1 - np.geomspace(1e-12, 1.0, 4000)[::-1]
        w = np.concatenate([lo + (w_in - lo) * crowd, w_in + (moment.w0 - w_in) * np.geomspace(1e-12, 1.0, 4000),
                            w_out - (w_out - moment.w0) * np.geomspace(1e-12, 1.0, 4000)[::-1]])
        t, r = moment.static(w)
        curve = np.column_stack([(r - X0) / (X1 - X0), (t - Y0) / (Y1 - Y0)])
        checked = 0
        for line in mark["lines"]:
            for u in line:
                self.assertLess(float(np.min(np.hypot(*(curve - u).T))), 6e-4, f"{key} {mark['label']} at {u}")
                checked += 1
            for u in (line[0], line[-1]):
                edge = min(u[0], 1 - u[0], u[1], 1 - u[1]) < 1.5e-4
                wall = math.hypot(u[0] - curve[0][0], u[1] - curve[0][1]) < 1.5e-4
                self.assertTrue(edge or wall, f"{key} {mark['label']} stops at {u}")
        return checked

    def check_bounce_moment(self, key, view, surface, mark):
        """A moment of constant r between the horizons of the black bounce, a = r_s/2: one upright line
        at r, or at rho = sqrt(r^2 + a^2) in the areal chart, over the stretch |ct| <= r_s the cylinder
        reaches, about t = 0 or, in the ingoing chart, about v = r_*."""
        X0, X1, Y0, Y1 = view["box"]
        chart = key.split("/")[1]
        r = self.sv_inside_r(surface)
        at = math.sqrt(r * r + 0.25) if chart == "areal" else r
        middle = sv_rstar(r, 0.5) - r if chart.endswith("ingoing") else 0.0
        t0, t1 = self.reach(surface)
        self.assertEqual(len(mark["lines"]), 1, key)
        (line,) = mark["lines"]
        for u in line:
            self.assertAlmostEqual(X0 + u[0] * (X1 - X0), at, delta=1.5e-4 * (X1 - X0), msg=f"{key} {mark['label']}")
        ends = sorted(Y0 + u[1] * (Y1 - Y0) for u in (line[0], line[-1]))
        for got, want in zip(ends, (middle + t0, middle + t1)):
            self.assertAlmostEqual(got, max(Y0, min(Y1, want)), delta=1.5e-4 * (Y1 - Y0), msg=f"{key} {mark['label']}")
        return len(line)

    def check_region_from_above(self, key, view, surface, mark):
        """Kerr's equator from above: the plane outside the horizon, the box less the disc of r_+,
        or, where the embedding ends inside the box, as Kerr-de Sitter's does at its cosmological
        horizon, the ring out to the circle where it ends."""
        X0, X1, Y0, Y1 = view["box"]
        lo, hi = self.reach(surface)
        outer, hole = mark["fills"][0]
        if hi < min(X1, Y1) * (1 + 1e-9):
            for u in outer:
                self.assertAlmostEqual(math.hypot(X0 + u[0] * (X1 - X0), Y0 + u[1] * (Y1 - Y0)), hi, delta=2e-3, msg=key)
        else:
            self.assertEqual(sorted(map(tuple, outer)), sorted([(0, 0), (1, 0), (1, 1), (0, 1)]), key)
        for u in hole:
            self.assertAlmostEqual(math.hypot(X0 + u[0] * (X1 - X0), Y0 + u[1] * (Y1 - Y0)), lo, delta=2e-3, msg=key)

    def test_every_slice_on_a_figure_lies_on_its_floor(self):
        """A figure's slice is a region of its floor, t = 0, carried back from the page through
        its camera: a disc or a ring about the axis out to where the floor or the embedding
        ends, and a rim where the embedding stops short of the floor's edge."""
        for metric_id, data in self.diagrams.items():
            for figure in [f for views in data.get("projections", {}).values() for f in views]:
                a, e = (math.radians(figure["camera"][k]) for k in ("azimuth", "elevation"))

                def floor(p):
                    # (X, Y, 0) is seen at (-X sin a + Y cos a, -sin e (X cos a + Y sin a)).
                    along = -p[1] / math.sin(e)
                    return math.hypot(along * math.cos(a) - p[0] * math.sin(a), along * math.sin(a) + p[0] * math.cos(a))
                for mark in figure.get("slices", []):
                    _, surface = self.moment(metric_id, mark)
                    rings = [[floor(p) for p in ring] for rings in mark["fills"] for ring in rings]
                    rims = [[floor(p) for p in line] for line in mark["lines"]]
                    where = f"{metric_id}/{figure['id']}"
                    if metric_id in ("cosmic_string", "point_particle_2plus1"):
                        hi = self.reach(surface, "conical", reference=True)[1]
                        self.assertLessEqual(max(max(r) for r in rings), hi + 1e-3, where)
                        self.assertTrue(all(abs(r - hi) < 1e-3 for rim in rims for r in rim), where)
                        continue
                    if metric_id == "tippett_tsang":
                        # The moment ct = A/2 stands above the floor, the rectangle of the height's plane.
                        grid = next(p for p in surface["pieces"] if "grid" in p)["grid"]
                        (ring,) = [ring for rings in mark["fills"] for ring in rings]
                        corners = set()
                        for p in ring:
                            along = -(p[1] - 0.5 * math.cos(e)) / math.sin(e)
                            corners.add((round(along * math.cos(a) - p[0] * math.sin(a), 3),
                                         round(along * math.sin(a) + p[0] * math.cos(a), 3)))
                        self.assertEqual(corners, {(x, y) for x in (grid["u"][0], grid["u"][-1])
                                                   for y in (grid["v"][0], grid["v"][-1])}, where)
                        self.assertEqual(rims, [], where)
                        continue
                    for ring in rings + rims:
                        self.assertLess(max(ring) - min(ring), 2e-3 * max(ring), where)
                    if metric_id in ("godel", "som_raychaudhuri", "vuorio_warped_ads"):
                        self.assertAlmostEqual(rims[0][0], self.reach(surface)[1], delta=1e-3, msg=where)
                    if metric_id in ("kerr", "kerr_newman"):
                        self.assertAlmostEqual(min(r[0] for r in rings), self.reach(surface)[0], delta=1e-3, msg=where)

    def test_every_slice_on_a_conformal_diagram_lies_on_its_moment(self):
        """A moment of constant t through a bifurcation point or a centre is the line T = 0,
        Reissner-Nordstrom's inside r_- the line T = pi, the closed universe's T = eta, and the
        Malament-Hogarth and Einstein-Rosen moments ct = tan p + tan q over two, as
        p, q = arctan(ct -+ r) make them; Vaidya's are carried back through each side of the shell's own map."""
        for metric_id, data in self.conformal.items():
            for view in data["views"]:
                X0, X1, T0, T1 = view["box"]
                for mark in view.get("slices", []):
                    _, surface = self.moment(metric_id, mark)
                    where = f"conformal {metric_id}/{view['id']} {mark['label']}"
                    points = [p for line in mark["lines"] for p in line] + mark["points"]
                    self.assertTrue(all(X0 <= X <= X1 and T0 <= T <= T1 for X, T in points), where)
                    t = surface.get("time")
                    if metric_id == "extended_inflation":
                        # tan((T-X)/2)=eta-r and tan((T+X)/2)=eta+r, with eta=-1.
                        lo, hi = self.reach(surface)
                        for X, T in points:
                            tp, tq = math.tan((T-X)/2), math.tan((T+X)/2)
                            scale = 1 + tp*tp + tq*tq
                            self.assertLess(abs((tp+tq)/2+1), 2e-4*scale, where)
                            self.assertLessEqual(lo-2e-4*scale, (tq-tp)/2, where)
                            self.assertLessEqual((tq-tp)/2, hi+2e-4*scale, where)
                    elif metric_id == "exponential_scalar_universes":
                        p = {"inflation": 2, "tracking": 2/3, "coasting": 1}[mark["view"]]
                        eta = 0 if p == 1 else 1/(1-p)
                        lo, hi = self.reach(surface)
                        for X, T in points:
                            u, v = math.tan((T-X)/2), math.tan((T+X)/2)
                            tolerance = 2e-4*(1+u*u+v*v)
                            self.assertLess(abs((u+v)/2-eta), tolerance, where)
                            self.assertLessEqual(lo-tolerance, (v-u)/2, where)
                            self.assertLessEqual((v-u)/2, hi+tolerance, where)
                    elif metric_id == "charach_malin":
                        # p,q=atan(1 -+ theta) on each static matter slice.
                        lo, hi = self.reach(surface)
                        for X, T in points:
                            tp, tq = math.tan((T-X)/2), math.tan((T+X)/2)
                            scale = 1 + tp*tp + tq*tq
                            self.assertLess(abs((tp+tq)/2-1), 2e-4*scale, where)
                            self.assertLessEqual(lo-2e-4*scale, (tq-tp)/2, where)
                            self.assertLessEqual((tq-tp)/2, hi+2e-4*scale, where)
                    elif metric_id == "stiff_fluid_universes":
                        # At t=1/2,x=0, p=q=atan(2), hence X=0,T=2 atan(2).
                        self.assertEqual(len(points), 1, where)
                        X, T = points[0]
                        self.assertLess(abs(X), 2e-4, where)
                        self.assertLess(abs(T - 2 * math.atan(2)), 2e-4, where)
                    elif metric_id == "frw":
                        eta = bisect(lambda e: e - math.sin(e) - t, 0, 2 * math.pi)
                        self.assertTrue(all(abs(T - eta) < 2e-4 for _, T in points), where)
                        self.assertEqual(sorted(X for X, _ in points), [0, round(math.pi, 4)], where)
                    elif metric_id == "draining_bathtub":
                        # Kruskal's U = tan p and V = tan q with UV = (1 - r) e^(2r)/(1 + r) and V = e^v,
                        # v = t + r - ln(1 + r): the Kerr-like chart's T = 0 is the level line, and the
                        # laboratory's t = 0 is carried back through r.
                        for X, T in points:
                            if mark["view"] == "funnel":
                                self.assertLess(abs(T), 2e-4, where)
                                continue
                            tp, tq = math.tan((T - X) / 2), math.tan((T + X) / 2)
                            # Beside the drain UV = 1 - 2r^3/3, so four decimals of a point no longer fix r.
                            if not (1e-3 < tq < 1e3 and abs(tp) < 1e3 and tp * tq < 0.9):
                                continue
                            f = lambda r: (1 - r) * math.exp(2 * r) / (1 + r) - tp * tq
                            r = bisect(f, 0, 20)
                            self.assertLess(abs(math.log(tq) - r + math.log(1 + r)), 5e-3, f"{where} at {(X, T)}")
                    elif metric_id == "unruh_acoustic_hole":
                        # Kruskal's U = tan p and V = tan q with UV = (1 - r) e^(4r - 2 arctan r)/(1 + r) and
                        # V = e^(2v), v = t + r - arctan r: Unruh's tau = 0 is the level line, and the
                        # laboratory's t = 0 is carried back through r.
                        for X, T in points:
                            if mark["view"] == "funnel":
                                self.assertLess(abs(T), 2e-4, where)
                                continue
                            tp, tq = math.tan((T - X) / 2), math.tan((T + X) / 2)
                            # Beside the sink UV = 1 - 4r^5/5, so four decimals of a point no longer fix r.
                            if not (1e-3 < tq < 1e3 and abs(tp) < 1e3 and tp * tq < 0.9):
                                continue
                            f = lambda r: (1 - r) * math.exp(4 * r - 2 * math.atan(r)) / (1 + r) - tp * tq
                            r = bisect(f, 0, 20)
                            self.assertLess(abs(math.log(tq) / 2 - r + math.atan(r)), 5e-3, f"{where} at {(X, T)}")
                    elif metric_id == "moncrief_taub_nut":
                        # At t=T=1 the mapped compact circle has UV=exp(2w(1)).
                        import sys
                        sys.path.insert(0, str(build.ROOT / "_tools" / "derivations"))
                        import moncrief_taub_nut
                        expected = math.exp(2*float(moncrief_taub_nut.conformal_w(1.0)))
                        for X, T in points:
                            p, q = (T-X)/2, (T+X)/2
                            residual = math.sin(p)*math.sin(q)-expected*math.cos(p)*math.cos(q)
                            self.assertLess(abs(residual), 2e-4, where)
                    elif metric_id == "brane_wall_universe":
                        # The wall's own moment, p = arctan U and q = arctan V along it, as
                        # brane_wall_universe.Moment follows Kruskal's U and V out from the wall.
                        import numpy as np
                        lo, hi = self.reach(surface)
                        U, V = self.wall_moment().Moment(1.0).kruskal(np.linspace(lo, hi, 3000))
                        p, q = np.arctan(U), np.arctan(V)
                        a, b = np.column_stack([q - p, p + q])[:-1], np.column_stack([q - p, p + q])[1:]
                        for X, T in points:
                            # The distance from the point to the polyline through the moment's samples.
                            d = b - a
                            s = np.clip(((X - a[:, 0]) * d[:, 0] + (T - a[:, 1]) * d[:, 1]) / np.maximum((d * d).sum(1), 1e-300), 0, 1)
                            gap = np.hypot(a[:, 0] + s * d[:, 0] - X, a[:, 1] + s * d[:, 1] - T)
                            self.assertLess(float(np.min(gap)), 2e-4, f"{where} at {(X, T)}")
                    elif metric_id in ("malament_hogarth", "einstein_rosen_waves", "twisting_cylindrical_waves", "gowdy",
                                       "senovilla"):
                        lo, hi = self.reach(surface)
                        for X, T in points:
                            p, q = (T - X) / 2, (T + X) / 2
                            tp, tq = math.tan(p), math.tan(q)
                            scale = 1 + tp * tp + tq * tq
                            self.assertLess(abs((tp + tq) / 2 - t), 2e-4 * scale, where)
                            self.assertLessEqual(lo - 2e-4 * scale, (tq - tp) / 2, where)
                            self.assertLessEqual((tq - tp) / 2, hi + 2e-4 * scale, where)
                    elif metric_id == "big_rip":
                        # The moment t at w = -3/2 is eta = -(-t)^(7/3) in units of eta_0, carried by
                        # p, q = arctan(eta -+ r).
                        eta = -(-t) ** (7 / 3)
                        for X, T in points:
                            tp, tq = math.tan((T - X) / 2), math.tan((T + X) / 2)
                            scale = 1 + tp * tp + tq * tq
                            self.assertLess(abs((tp + tq) / 2 - eta), 2e-4 * scale, f"{where} at {(X, T)}")
                    elif metric_id == "bondi_sachs":
                        # The sphere cu = 10 m_0, r = 10 m_0 is one event of the axis: on the cone
                        # p = arctan(10/40) and on the world tube the view draws.
                        ((X, T),) = points
                        self.assertLess(abs((T - X) / 2 - math.atan(0.25)), 2e-4, where)
                        tube = next(layer["points"] for layer in view["layers"] if layer["class"] == "boundary")

                        def off(a, b):
                            (ax, ay), (bx, by) = a, b
                            k = ((X - ax) * (bx - ax) + (T - ay) * (by - ay)) / ((bx - ax) ** 2 + (by - ay) ** 2)
                            k = min(max(k, 0.0), 1.0)
                            return math.hypot(X - ax - k * (bx - ax), T - ay - k * (by - ay))
                        self.assertLess(min(off(a, b) for a, b in zip(tube, tube[1:])), 5e-3, where)
                    elif metric_id == "domain_wall":
                        # Each side is Minkowski's triangle cut at the wall X = pi/2, the side z > 0
                        # mirrored in it, and the moment kct runs along cT = R tanh(kct) on both.
                        for X, T in points:
                            p, q = (T - X) / 2, (T + X) / 2
                            a, b = (math.tan(p), math.tan(q)) if X <= math.pi / 2 else (math.tan(q - math.pi / 2),
                                                                                        math.tan(p + math.pi / 2))
                            cT, R = (a + b) / 2, (b - a) / 2
                            self.assertLess(abs(cT - R * math.tanh(t)), 2e-4 * (1 + a * a + b * b), f"{where} at {(X, T)}")
                    elif metric_id in ("misner", "gott_time_machine"):
                        # The hyperbola (ct - x)(ct + x) = c^2t^2 of the covering plane, every copy.
                        for X, T in points:
                            tp, tq = math.tan((T - X) / 2), math.tan((T + X) / 2)
                            self.assertLess(abs(tp * tq - t * t), 2e-4 * (1 + tp * tp + tq * tq), where)
                    elif metric_id == "self_creating_universe":
                        # On de Sitter's square V = tan T and W = cos(X)/cos T at r_0 = 1, and the moment
                        # c tau = t of the Kantowski-Sachs chart is V^2 - W^2 = sinh^2(t), every copy.
                        for X, T in points:
                            if abs(math.cos(T)) < 1e-3:
                                continue
                            V, W = math.tan(T), math.cos(X) / math.cos(T)
                            self.assertLess(abs(V * V - W * W - math.sinh(t) ** 2), 2e-3 * (1 + V * V + W * W), where)
                    elif metric_id == "ori_time_machine":
                        # The hyperbola of the covering plane on which T = t, every copy: its null
                        # coordinates -2 e^(-z/2) and 2T e^(z/2) multiply to -4t.
                        for X, T in points:
                            tp, tq = math.tan((T - X) / 2), math.tan((T + X) / 2)
                            self.assertLess(abs(tp * tq + 4 * t), 2e-4 * (1 + tp * tp + tq * tq), where)
                    elif metric_id == "roberts":
                        # arctan(u) and arctan(v), with ct = (u/s + s v)/2 and rho = (s v - u/s)/2 for
                        # s = sqrt(1 + p), over the rho the embedding reaches.
                        s = math.sqrt(1 + ROBERTS_P[view["id"].rsplit("_", 1)[1]])
                        lo, hi = self.reach(surface)
                        for X, T in points:
                            u, v = math.tan((T - X) / 2), math.tan((T + X) / 2)
                            scale = 1 + u * u + v * v
                            self.assertLess(abs((u / s + s * v) / 2 - t), 2e-4 * scale, f"{where} at {(X, T)}")
                            self.assertLessEqual(lo - 2e-4 * scale, (s * v - u / s) / 2, where)
                            self.assertLessEqual((s * v - u / s) / 2, hi + 2e-4 * scale, where)
                    elif metric_id == "hotta_tanaka":
                        # p = arctan(u) - pi/4 and q = arctan(v + Theta(u)/2) + pi/4 of the Kruskal chart.
                        (X, T), = points
                        u = math.tan((T - X) / 2 + math.pi / 4)
                        v = math.tan((T + X) / 2 - math.pi / 4) - (0.5 if u > 1e-9 else 0.0)
                        want = hotta_tanaka_moment(t)
                        self.assertLess(abs(u - want[0]) + abs(v - want[1]), 5e-4, where)
                    elif metric_id in ("maximally_supersymmetric_plane_wave", "cahen_wallach"):
                        # The wave front at mu c u = t meets the axis at the event v = 0, which
                        # p = arctan(2 mu v) - pi/2 and q = mu c u + pi/2 put at T = t, X = t + pi.
                        (X, T), = points
                        self.assertLess(abs(T - t) + abs(X - t - math.pi), 2e-4, where)
                    elif metric_id == "vaidya":
                        for X, T in points:
                            p, q = (T - X) / 2, (T + X) / 2
                            if q < math.pi / 4 - 1e-3:        # inside the shell, v < 0, flat
                                F = lambda u: math.atan((1 + u / 2) * math.exp(-u / 2))
                                v = bisect(lambda u: F(u) - q, -60, 0)
                                r = (v - bisect(lambda u: F(u) - p, -80, 0)) / 2
                            elif q > math.pi / 4 + 1e-3:      # outside, Kruskal's p = arctan U, q = arctan V
                                v = 2 * math.log(math.tan(q))
                                # (1 - r) e^r falls from its greatest value at r = 0, which is where
                                # the moment ends on the singularity.
                                f = lambda r: (1 - r) * math.exp(r - v / 2) - math.tan(p)
                                r = 0.0 if f(0) < 1e-4 else bisect(f, 0, 20)
                            else:
                                continue
                            self.assertLess(abs(v - r - t), 2e-3 * (1 + abs(v)), f"{where} at {(X, T)}")
                    elif metric_id == "black_to_white_hole":
                        # Inside the shell p, q = arctan(u, v), and the moment is the flat time
                        # (u + v)/2 = -R at which Schwarzschild's slice ct = t meets the shell. Outside it
                        # each point is carried back to Kruskal's U through the interior ray tan p and to
                        # V through the drawing's q, and lies on t = ln(-V/U).
                        sys.path.insert(0, str(build.ROOT / "_tools" / "derivations"))
                        import black_to_white_hole as bwh
                        level = -float(bwh.shell_radius_schwarzschild(t))
                        inside, outside = mark["lines"]
                        for X, T in inside:
                            a, b = math.tan((T - X) / 2), math.tan((T + X) / 2)
                            self.assertLess(abs((a + b) / 2 - level), 2e-3 * (1 + a * a + b * b), f"{where} at {(X, T)}")
                        self.assertLess(math.dist(inside[-1], outside[0]), 2e-4, where)
                        for X, T in outside:
                            p, q = (T - X) / 2, (T + X) / 2
                            if not (-1.4 < p < -0.2 and 0.05 < q < 1.4):
                                continue        # beside the infinities four decimals no longer fix the event
                            U = float(bwh.shell_U(math.tan(p)))
                            V = bisect(lambda x: float(bwh.q_of_V(x)) - q, bwh.V0, 60)
                            self.assertLess(abs(math.log(-V / U) - t), 2e-2, f"{where} at {(X, T)}")
                    elif metric_id == "charged_shell" and view["id"] == "exterior_isotropic":
                        # The balanced shell at rest: t = 0 is the level line T = 0 on both sides.
                        for line in mark["lines"]:
                            for X, T in line:
                                self.assertLess(abs(T), 1e-3, f"{where} at {(X, T)}")
                    elif metric_id == "charged_shell":
                        # Inside the shell p, q = arctan((cT -+ r)/L) with L the shell's T + R where it comes
                        # back to r_-: the moment is the flat time at which the slice v - r = t meets the
                        # shell, out to the shell. Outside it each point is carried back along its two rays
                        # to the shell, charged_event, and lies on v - r = t.
                        length = charged_inner_time(CHARGED_ETA_M) + CHARGED_RM
                        on_shell = charged_on_slice(t)
                        level, edge = charged_inner_time(on_shell), charged_radius(on_shell)
                        inside, outside = mark["lines"]
                        for X, T in inside:
                            a, b = length * math.tan((T - X) / 2), length * math.tan((T + X) / 2)
                            self.assertLess(abs((a + b) / 2 - level), 2e-3 * (1 + a * a + b * b), f"{where} at {(X, T)}")
                            self.assertLessEqual((b - a) / 2, edge + 2e-3 * (1 + a * a + b * b), where)
                        self.assertLess(math.dist(inside[-1], outside[0]), 2e-4, where)
                        for X, T in outside[1:]:
                            v, r = charged_event((T - X) / 2, (T + X) / 2, length)
                            if min(abs(r - CHARGED_RP), abs(r - CHARGED_RM)) < 0.02:
                                continue        # on a horizon the rounded point's ray is not well enough known
                            self.assertLess(abs(v - r - t), 5e-3 * (1 + abs(v)), f"{where} at {(X, T)}")
                    elif metric_id == "penrose_impulsive_wave":
                        # At beta = 1/2: behind the front p, q = arctan(ct -+ r), so ct = (tan p + tan q)/2 is t;
                        # ahead of it p = arctan((cT - R)/beta), q = arctan(beta (cT + R)), so cT is t/beta; the two
                        # lines meet on the front. On the null plane rho = 1 the event is U = 0, V = t/sqrt(2).
                        if view["id"] == "null":
                            (X, T), = mark["points"]
                            self.assertLess(abs(math.tan((T - X) / 2)) + abs(math.tan((T + X) / 2) - t / math.sqrt(2)),
                                            2e-3, where)
                        else:
                            inside, outside = mark["lines"]
                            for X, T in inside:
                                a, b = math.tan((T - X) / 2), math.tan((T + X) / 2)
                                self.assertLess(abs((a + b) / 2 - t), 2e-3 * (1 + a * a + b * b), f"{where} at {(X, T)}")
                            for X, T in outside:
                                a, b = math.tan((T - X) / 2), math.tan((T + X) / 2)
                                self.assertLess(abs((a / 2 + 2 * b) / 2 - 2 * t), 2e-3 * (1 + a * a + b * b),
                                                f"{where} at {(X, T)}")
                            self.assertLess(math.dist(inside[-1], outside[0]), 2e-4, where)
                    elif metric_id == "israel_shell":
                        # Inside the shell p, q = arctan(3(cT -+ r)/7): the moment is the flat time at which
                        # the slice v - r = t meets the shell, out to the shell. Outside it each point is
                        # carried back along its two rays to the shell, israel_event, and lies on v - r = t.
                        if t < ISRAEL_END:
                            on_shell = israel_on_slice(t)
                            level, edge = israel_inner_time(on_shell), israel_radius(on_shell)
                            inside, outside = mark["lines"]
                            for X, T in inside:
                                a, b = 7 / 3 * math.tan((T - X) / 2), 7 / 3 * math.tan((T + X) / 2)
                                self.assertLess(abs((a + b) / 2 - level), 2e-3 * (1 + a * a + b * b), f"{where} at {(X, T)}")
                                self.assertLessEqual((b - a) / 2, edge + 2e-3 * (1 + a * a + b * b), where)
                            self.assertLess(math.dist(inside[-1], outside[0]), 2e-4, where)
                        else:
                            (outside,) = mark["lines"]
                        for X, T in outside[1:]:
                            v, r = israel_event((T - X) / 2, (T + X) / 2)
                            self.assertLess(abs(v - r - t), 5e-3 * (1 + abs(v)), f"{where} at {(X, T)}")
                    elif metric_id == "bonnor_vaidya":
                        # Both views draw the moment before the shell turns: outside the shell through
                        # the tower's ingoing map, and inside it through the flat map, every p moved by pi/2.
                        for X, T in points:
                            p, q = (T - X) / 2 + math.pi / 2, (T + X) / 2
                            if q < math.pi / 4 - 1e-3:
                                v = bisect(lambda u: bv_crossing(u) - q - math.pi / 2, -80, -1e-12)
                                r = (v - bisect(lambda u: bv_crossing(u) - p, -80, -1e-12)) / 2
                            elif q > math.pi / 4 + 1e-3:
                                v = math.log(math.tan(q)) / BV_KAPPA
                                r = bv_radius(v, p)
                            else:
                                continue
                            if r is not None:
                                self.assertLess(abs(v - r - t), 2e-3 * (1 + abs(v)), f"{where} at {(X, T)}")
                    elif metric_id == "lindquist_wheeler_lattice" and view["id"] == "hypersphere":
                        # The rectangle of the cycloid parameter: T = eta with eta + sin eta = 2 c tau/a_m.
                        eta = bisect(lambda e: e + math.sin(e) - 2 * t * math.sin(LW_PSI) ** 3, 0, math.pi)
                        self.assertTrue(all(abs(T - eta) < 2e-4 for _, T in points), where)
                        self.assertEqual(sorted(round(X, 3) for X, _ in points),
                                         [round(LW_PSI, 3), round(math.pi - LW_PSI, 3)], where)
                    elif metric_id == "lindquist_wheeler_lattice":
                        # Kruskal's U V = (1 - r) e^r gives the radius of each point, the cycloids the shell
                        # that has that radius at the moment's proper time, and that shell its own V.
                        lo, hi = self.reach(surface, "lindquist_wheeler")
                        for X, T in points:
                            U, V = math.tan((T - X) / 2), math.tan((T + X) / 2)
                            if abs(T) > 1.5:
                                # Beside r = 0, the line T = pi/2, the tangents of the rounded point no longer
                                # fix the radius; a moment ends there once the throat has gone, c tau > pi/2.
                                self.assertGreater(t, math.pi / 2, where)
                                continue
                            r = bisect(lambda x: (1 - x) * math.exp(x) - U * V, 0.0, 3.0)

                            def radius(R):
                                eta = bisect(lambda e: e + math.sin(e) - 2 * t / R ** 1.5, 0.0, math.pi)
                                return R * math.cos(eta / 2) ** 2, eta
                            R = bisect(lambda R: radius(R)[0] - r, lo + 1e-9, hi) if t else r
                            eta = radius(R)[1] if t else 0.0
                            k = math.sqrt(max(R - 1, 0.0))
                            want = (k * math.cos(eta / 2) + math.sin(eta / 2)) * math.exp(
                                (radius(R)[0] + k * (eta + R / 2 * (eta + math.sin(eta)))) / 2)
                            self.assertLess(abs(math.atan(want) - math.atan(V)), 2e-3, f"{where} at {(X, T)}")
                    elif metric_id == "oppenheimer_snyder":
                        chi0 = math.pi / 4
                        eta = bisect(lambda e: math.sqrt(2) * (e + math.sin(e)) - t, 0, math.pi)
                        dust = [(X, T) for X, T in points if X <= chi0 + 1e-4]
                        self.assertTrue(dust and all(abs(T - eta) < 2e-4 for _, T in dust), where)
                        # The line crosses the star's surface at (chi_0, eta), where Novikov's
                        # slice outside meets the dust's moment.
                        inside, outside = mark["lines"]
                        self.assertEqual(inside, [[0, round(eta, 4)], [round(chi0, 4), round(eta, 4)]], where)
                        self.assertLess(math.dist(outside[0], [chi0, eta]), 2e-4, where)
                    elif metric_id == "semiclosed_world":
                        # The line T = eta across the dust to chi_0 = 3 pi/4, sqrt 2 (eta + sin eta) = c tau,
                        # and from the surface Novikov's curve outward through the throat.
                        chi0 = 3 * math.pi / 4
                        eta = bisect(lambda e: math.sqrt(2) * (e + math.sin(e)) - t, 0, math.pi) if t else 0.0
                        inside, outside = mark["lines"]
                        self.assertEqual(inside, [[0, round(eta, 4)], [round(chi0, 4), round(eta, 4)]], where)
                        self.assertLess(math.dist(outside[0], [chi0, eta]), 2e-4, where)
                        self.assertTrue(all(b[0] > a[0] for a, b in zip(outside, outside[1:])), where)
                        self.assertTrue(all(-2e-4 <= T < math.pi for _, T in outside), where)
                    elif metric_id == "datt_ruban_t_models":
                        # Inside, p, q = arctan(sigma -+ r) at the moment's eta = pi + e, (e + sin e)/2 = c tau,
                        # from r = -2 r_s to the surface, X = 0; outside, Novikov's curve from the surface out.
                        eta = math.pi + (bisect(lambda e: (e + math.sin(e)) / 2 - t, 0, math.pi) if t else 0.0)
                        sigma = ruban_sigma(eta)
                        inside, outside = mark["lines"]
                        for X, T in inside:
                            tp, tq = math.tan((T - X) / 2), math.tan((T + X) / 2)
                            self.assertLess(abs((tp + tq) / 2 - sigma), 2e-3 * (1 + tp * tp + tq * tq), f"{where} at {(X, T)}")
                            self.assertLessEqual((tq - tp) / 2, 1e-9)
                        tp, tq = math.tan((inside[0][1] - inside[0][0]) / 2), math.tan((inside[0][1] + inside[0][0]) / 2)
                        self.assertAlmostEqual((tq - tp) / 2, -2.0, delta=5e-3)
                        self.assertLess(math.dist(inside[-1], [0, 2 * math.atan(sigma)]), 2e-4, where)
                        self.assertLess(math.dist(outside[0], inside[-1]), 2e-4, where)
                        self.assertTrue(all(b[0] > a[0] for a, b in zip(outside, outside[1:])), where)
                        self.assertTrue(all(-2e-4 <= T < 2 * math.atan(ruban_sigma(2 * math.pi)) for _, T in outside), where)
                    elif metric_id == "white_hole":
                        chi0 = math.pi / 4
                        eta = bisect(lambda e: math.sqrt(2) * (e - math.sin(e)) - t, 0, 2 * math.pi) - math.pi
                        dust = [(X, T) for X, T in points if X <= chi0 + 1e-4]
                        self.assertTrue(dust and all(abs(T - eta) < 2e-4 for _, T in dust), where)
                        # The line crosses the core's surface at (chi_0, eta - pi), where Novikov's
                        # slice outside meets the dust's moment.
                        inside, outside = mark["lines"]
                        self.assertEqual(inside, [[0, round(eta, 4)], [round(chi0, 4), round(eta, 4)]], where)
                        self.assertLess(math.dist(outside[0], [chi0, eta]), 2e-4, where)
                    elif metric_id == "kantowski_sachs" and mark["view"] == "dust":
                        # p, q = arctan(tau -+ r), tau the integral of 2 cos^3(s)/(cos(s) + s sin(s))
                        # from 0 to eta, by Simpson's rule.
                        n = 2000
                        f = [2 * math.cos(t * k / n) ** 3 / (math.cos(t * k / n) + t * k / n * math.sin(t * k / n))
                             for k in range(n + 1)]
                        tau = t / n / 3 * (f[0] + f[-1] + 4 * sum(f[1:-1:2]) + 2 * sum(f[2:-1:2]))
                        for X, T in points:
                            tp, tq = math.tan((T - X) / 2), math.tan((T + X) / 2)
                            self.assertLess(abs((tp + tq) / 2 - tau), 2e-4 * (1 + tp * tp + tq * tq), f"{where} at {(X, T)}")
                    elif metric_id == "kantowski_sachs":
                        # Kruskal's square: tan p tan q = (1 - T/r_s) e^(T/r_s) inside the horizon.
                        radius = surface["pieces"][0]["points"][0][1]
                        for X, T in points:
                            tp, tq = math.tan((T - X) / 2), math.tan((T + X) / 2)
                            self.assertLess(abs(tp * tq - (1 - radius) * math.exp(radius)), 2e-3 * (1 + tp * tp) * (1 + tq * tq),
                                            f"{where} at {(X, T)}")
                    elif metric_id == "thakurta":
                        # tan(p)=(eta-x)/4, tan(q)=(eta+x)/4 on the fixed eta=1 surface.
                        for X, T in points:
                            tp, tq = math.tan((T - X) / 2), math.tan((T + X) / 2)
                            self.assertLess(abs(2 * (tp + tq) - 1), 5e-4 * (1 + tp * tp + tq * tq),
                                            f"{where} at {(X, T)}")
                    elif metric_id == "sultana_dyer":
                        # Kruskal's map of the ingoing chart at r_s = 1: tan q = e^((eta + r)/2) and
                        # tan p = (1 - r) e^((r - eta)/2), so each point gives its r, which runs from the
                        # centre out to the embedding's reach.
                        hi = self.reach(surface)[1]
                        for X, T in points:
                            tp, tq = math.tan((T - X) / 2), math.tan((T + X) / 2)
                            r = 2 * math.log(tq) - t
                            self.assertTrue(-2e-3 < r < hi + 2e-2, f"{where} at {(X, T)}")
                            self.assertLess(abs(tp - (1 - r) * math.exp((r - t) / 2)), 2e-3 * (1 + tp * tp) * (1 + abs(r)),
                                            f"{where} at {(X, T)}")
                    elif metric_id == "mcvittie":
                        # p = F(s_out) and q = -F(s_in), the times the event's two rays left R = r_s:
                        # each ray is run forward again from the throat x = 0 to the moment's time, by
                        # Runge and Kutta's rule in ln t, and both must arrive at one x, R = cosh^2 x.
                        h = 1 / math.sqrt(15)
                        kappa = (1 / 1.0851996154371 ** 2 - 2 * 1.0851996154371 / 15) / 2

                        def left_at(value):
                            return math.exp(bisect(lambda y: math.atan(y / 2 + math.log(2) / 2 + math.expm1(
                                kappa * math.exp(y)) / 20) - value, -60, math.log(60)))

                        def arrives(s, sign, n=1500):
                            def rate(y, x):
                                tt = math.exp(y)
                                return (h * tt / math.tanh(1.5 * h * tt) + sign * tt * math.tanh(x) / math.cosh(x) ** 2) / 2
                            y, x, d = math.log(s), 0.0, (math.log(t) - math.log(s)) / n
                            for _ in range(n):
                                k1 = rate(y, x)
                                k2 = rate(y + d / 2, x + d * k1 / 2)
                                k3 = rate(y + d / 2, x + d * k2 / 2)
                                k4 = rate(y + d, x + d * k3)
                                y, x = y + d, x + d * (k1 + 2 * k2 + 2 * k3 + k4) / 6
                            return x
                        met = 0
                        for X, T in points:
                            p_, q_ = (T - X) / 2, (T + X) / 2
                            if max(abs(p_), abs(q_)) > 1.2:
                                continue    # where F has flattened, four decimals no longer fix the time
                            out, back = arrives(left_at(p_), 1), arrives(left_at(-q_), -1)
                            self.assertLess(abs(out - back), 5e-3 * (1 + out), f"{where} at {(X, T)}")
                            met += 1
                        self.assertGreater(met, 5, where)
                    elif metric_id == "nordstrom_scalar" and mark["view"] == "dust":
                        # p, q = arctan((ct -+ r)/L), so tan p + tan q = 2ct/L on a moment of the inertial time.
                        for X, T in points:
                            self.assertLess(abs(math.tan((T - X) / 2) + math.tan((T + X) / 2) - 2 * t), 2e-3,
                                            f"{where} at {(X, T)}")
                    elif metric_id == "small_universes":
                        # p, q = arctan(eta -+ x), so tan p + tan q = 2 eta: 2 (3t/2)^(1/3) on a moment of the
                        # torus, and ln 11 on the moment the horn is embedded at.
                        eta = math.log(11) if mark["view"] == "horn" else 2 * (1.5 * t) ** (1 / 3)
                        for X, T in points:
                            self.assertLess(abs(math.tan((T - X) / 2) + math.tan((T + X) / 2) - 2 * eta), 2e-3,
                                            f"{where} at {(X, T)}")
                    elif metric_id == "poincare_dodecahedral":
                        # X = chi and T = eta: today in Luminet's universe, eta = 0.38882, across one cell along
                        # the line to the centre of a face, from chi = -pi/10 to pi/10.
                        self.assertTrue(all(abs(T - 0.3888) < 2e-4 for _, T in points), where)
                        self.assertEqual(sorted(X for X, _ in points), [round(-math.pi / 10, 4), round(math.pi / 10, 4)],
                                         where)
                    elif metric_id == "lqc_bounce":
                        # p, q = arctan((eta -+ r)/c t_b), so tan p + tan q = 2 eta(t) on a moment of the observers' time.
                        for X, T in points:
                            self.assertLess(abs(math.tan((T - X) / 2) + math.tan((T + X) / 2)
                                                - 2 * lqc_bounce_eta(t)), 2e-3, f"{where} at {(X, T)}")
                    elif metric_id == "bekenstein_universes":
                        # The closed strip draws T = eta; the open universe's p, q = gd(eta -+ chi), so
                        # ln tan(p/2 + pi/4) + ln tan(q/2 + pi/4) = 2 eta on a moment.
                        for X, T in points:
                            if mark["view"] == "closed":
                                self.assertLess(abs(T - t), 2e-4, f"{where} at {(X, T)}")
                            else:
                                p_, q_ = (T - X) / 2, (T + X) / 2
                                total = math.log(math.tan(p_ / 2 + math.pi / 4)) + math.log(math.tan(q_ / 2 + math.pi / 4))
                                self.assertLess(abs(total - 2 * t), 2e-3 / math.cos(q_) ** 2, f"{where} at {(X, T)}")
                    elif metric_id == "murphy_viscous_universe":
                        # p, q = arctan(eta -+ x), so tan p + tan q = 2 eta = 2(a - 1/a) on the moment t, for q = 2.
                        a = murphy_scale(t)
                        for X, T in points:
                            tp, tq = math.tan((T - X) / 2), math.tan((T + X) / 2)
                            scale = 1 + tp * tp + tq * tq
                            self.assertLess(abs((tp + tq) / 2 - (a - 1 / a)), 2e-4 * scale, f"{where} at {(X, T)}")
                    elif metric_id == "kopczynski_trautman":
                        # p, q = arctan((eta -+ r)/l), so tan p + tan q = 2 eta(t)/l on a moment of the dust's time.
                        for X, T in points:
                            self.assertLess(abs(math.tan((T - X) / 2) + math.tan((T + X) / 2)
                                                - 2 * kopczynski_trautman_eta(t)), 2e-3, f"{where} at {(X, T)}")
                    elif metric_id == "coleman_de_luccia":
                        # p, q = arctan(e^(eta -+ chi)) with eta = ln tan(c tau/2), so tan p tan q = tan^2(c tau/2).
                        for X, T in points:
                            tp, tq = math.tan((T - X) / 2), math.tan((T + X) / 2)
                            self.assertLess(abs(tp * tq - math.tan(t / 2) ** 2), 2e-3 * (1 + tp * tp) * (1 + tq * tq),
                                            f"{where} at {(X, T)}")
                    elif metric_id == "milne":
                        # p, q = arctan(ct e^-chi), arctan(ct e^chi), so tan p tan q = c^2t^2.
                        for X, T in points:
                            tp, tq = math.tan((T - X) / 2), math.tan((T + X) / 2)
                            self.assertLess(abs(tp * tq - t * t), 2e-3 * (1 + tp * tp) * (1 + tq * tq), f"{where} at {(X, T)}")
                    elif metric_id == "ab_metrics":
                        # p, q = arctan(tau -+ x), so tan p + tan q = 2 tau on a moment of the de Sitter slicing.
                        for X, T in points:
                            self.assertLess(abs(math.tan((T - X) / 2) + math.tan((T + X) / 2) - 2 * t), 2e-3,
                                            f"{where} at {(X, T)}")
                    elif metric_id == "nariai":
                        # A global moment is the line tan(eta) = sinh(ct), Lambda = 1; the sphere
                        # is the event at eta = 0.
                        eta = 0.0 if t is None else math.atan(math.sinh(t))
                        self.assertTrue(all(abs(T - eta) < 2e-4 for _, T in points), where)
                    elif metric_id == "robinson_trautman":
                        # The wave front u = u_k is the null line p = -arctan(e^{-cu_k/4m}).
                        self.assertTrue(all(abs((T - X) / 2 + math.atan(math.exp(-t / 4))) < 2e-4 for X, T in points), where)
                    elif metric_id == "aichelburg_sexl":
                        # The wave front u = u_k is the null line p = arctan u_k.
                        self.assertTrue(all(abs((T - X) / 2 - math.atan(t)) < 2e-4 for X, T in points), where)
                    elif metric_id == "light_beam":
                        # Bonnor's wave front u = u_k is the null line p = arctan(sqrt(2) u_k/l), l = 4R.
                        self.assertTrue(all(abs((T - X) / 2 - math.atan(math.sqrt(2) * t / 4)) < 2e-4 for X, T in points),
                                        where)
                    elif metric_id == "khan_penrose":
                        # The event u = v = sin(t/2) where both waves have passed, drawn with p, q = u, v.
                        for X, T in points:
                            self.assertLess(abs(X) + abs(T - 2 * math.sin(t / 2)), 2e-4, where)
                    elif metric_id == "belinski_zakharov":
                        # The event on xi = 0 at tau = t, z = 0 at ct = w sinh(t): X = 0 and T = 2 arctan(sinh t).
                        for X, T in points:
                            self.assertLess(abs(X) + abs(T - 2 * math.atan(math.sinh(t))), 2e-4, where)
                    elif metric_id == "bell_szekeres":
                        # The event on eta = 0 at xi = t: on the plane x = y = 0, p = q = t/2 at a = b = 1,
                        # and on the strip of the anti-de Sitter factor, rho = 0 at chi = t - pi/2.
                        height = t if view["id"] in ("double_null", "time_space") else t - math.pi / 2
                        for X, T in points:
                            self.assertLess(abs(X) + abs(T - height), 2e-4, where)
                    elif metric_id == "chandrasekhar_xanthopoulos":
                        # The event on lambda = 0 at psi = t: p = q = t/2 in the null coordinates u and v.
                        for X, T in points:
                            self.assertLess(abs(X) + abs(T - t), 2e-4, where)
                    elif metric_id == "reissner_nordstrom_de_sitter":
                        # Every region is placed by g(x) = arctan e^(-kappa x) of u and v, kappa = 3/16. The
                        # static moment between r_+ and r_c is the line T = 0; the one inside r_- is T = pi
                        # above the black hole and T = -pi below the white hole, and on the cosmological
                        # chart's own time it is cT = Lukewarm.inside_t() below the white hole; the
                        # cosmological moment is carried back to u, v and r through the white hole, the
                        # static region and the expanding one.
                        L = Lukewarm
                        H = math.pi / 2

                        def back(y):
                            return -math.log(math.tan(y)) / L.KAPPA
                        if mark["view"] == "between":
                            self.assertTrue(all(abs(T) < 2e-4 for _, T in points), where)
                            continue
                        if mark["view"] == "inside" and view["id"] != "cosmological":
                            level = math.pi if view["id"] in ("static", "ingoing") else -math.pi
                            self.assertTrue(all(abs(T - level) < 2e-4 for _, T in points), where)
                            continue
                        met = 0
                        for X, T in points:
                            p_, q_ = (T - X) / 2, (T + X) / 2
                            if min(abs(p_ / H - round(p_ / H)), abs(q_ / H - round(q_ / H))) < 0.03:
                                continue    # beside a horizon, where four decimals no longer fix u or v
                            if mark["view"] == "inside":
                                u, v = (back(p_ + math.pi), -back(-q_)) if p_ < -H else (back(-p_), -back(q_ + math.pi))
                                self.assertLess(abs((u + v) / 2 - L.inside_t()), 5e-2, f"{where} at {(X, T)}")
                            else:
                                u = back(-p_)
                                v, (lo, hi) = ((-back(-q_), (0.5, L.ROOTS[1])) if q_ < 0 else
                                               (-back(q_), (L.ROOTS[1], L.ROOTS[0])) if q_ < H else
                                               (-back(math.pi - q_), (L.ROOTS[0], 60.0)))
                                grows = L.rstar(lo + 0.3 * (hi - lo)) < L.rstar(lo + 0.6 * (hi - lo))
                                r = bisect(lambda r: (L.rstar(r) - (v - u) / 2) * (1 if grows else -1),
                                           lo + 1e-12, hi - 1e-12)
                                self.assertLess(abs((u + v) / 2 - L.static_t(1 / L.H, r)), 5e-2, f"{where} at {(X, T)}")
                            met += 1
                        self.assertGreater(met, 3, where)
                    elif metric_id == "kastor_traschen":
                        # One hole at m = 1 and H = -3/16 is the lukewarm hole at r_s = 2m run backward: the
                        # event (tau, r) is drawn where the lukewarm hole's cosmological chart draws
                        # (-tau/2, r/2), turned upside down. So each point turned back, (X, -T), is carried
                        # to u, v and the areal radius as the lukewarm hole's are, and lies on the moment
                        # tau' = 4/3 of that chart, which is c tau = -8m/3.
                        L = Lukewarm
                        H = math.pi / 2
                        met = 0
                        for X, T in points:
                            p_, q_ = (-T - X) / 2, (-T + X) / 2
                            if min(abs(p_ / H - round(p_ / H)), abs(q_ / H - round(q_ / H))) < 0.03:
                                continue    # beside a horizon, where four decimals no longer fix u or v
                            u = -math.log(math.tan(-p_)) / L.KAPPA
                            y, (lo, hi) = ((-q_, (0.5, L.ROOTS[1])) if q_ < 0 else (q_, (L.ROOTS[1], L.ROOTS[0]))
                                           if q_ < H else (math.pi - q_, (L.ROOTS[0], 60.0)))
                            v = math.log(math.tan(y)) / L.KAPPA
                            grows = L.rstar(lo + 0.3 * (hi - lo)) < L.rstar(lo + 0.6 * (hi - lo))
                            r = bisect(lambda r: (L.rstar(r) - (v - u) / 2) * (1 if grows else -1), lo + 1e-12, hi - 1e-12)
                            self.assertLess(abs((u + v) / 2 - L.static_t(4 / 3, r)), 5e-2, f"{where} at {(X, T)}")
                            met += 1
                        self.assertGreater(met, 3, where)
                    elif metric_id == "rn_metric" and mark["view"] == "inside":
                        self.assertTrue(all(abs(T - math.pi) < 2e-4 for _, T in points), where)
                    elif metric_id == "simpson_visser" and mark["view"] == "inside":
                        # Kruskal's square of the black bounce, a = r_s/2: p, q = arctan e^(-+k(ct -+ r_*)) with
                        # k = sqrt(3)/4, so tan p tan q = e^(2 k r_*) on the moment of constant r.
                        want = math.exp(2 * math.sqrt(3) / 4 * sv_rstar(self.sv_inside_r(surface), 0.5))
                        for X, T in points:
                            tp, tq = math.tan((T - X) / 2), math.tan((T + X) / 2)
                            self.assertLess(abs(tp * tq - want), 2e-3 * (1 + tp * tp) * (1 + tq * tq), f"{where} at {(X, T)}")
                    elif metric_id in ("bardeen", "reissner_nordstrom_ads", "quantum_oppenheimer_snyder"):
                        # Inside r- the moment runs through the inner bifurcation sphere, at T = pi, and
                        # outside r+ through the outer one, at T = 0, or a period up, at T = 2 pi, on the
                        # outgoing chart's view, whose exterior is the one above the white hole.
                        height = math.pi if mark["view"] == "inside" else 2 * math.pi if view["id"] == "outgoing" else 0.0
                        self.assertTrue(all(abs(T - height) < 2e-4 for _, T in points), where)
                    elif metric_id in ("hayward", "regular_black_hole_pure_gravity") and mark["view"] == "inside":
                        # Through the inner bifurcation point above the exteriors on the static and ingoing
                        # views, and through the one below them on the outgoing view, whose chart covers it.
                        height = -math.pi if view["id"] == "outgoing" else math.pi
                        self.assertTrue(all(abs(T - height) < 2e-4 for _, T in points), where)
                    elif metric_id == "hiscock":
                        # Both rays of every event are traced, so the slice is held to what needs no tracing:
                        # it runs outward with its ingoing rays arriving later and later, never left of the
                        # first centre, and a moment before the shell reaches the centre, t < 0, starts on
                        # that centre, X = 0, at p = q = arctan((t + 1)/3).
                        qs = [(T + X) / 2 for X, T in points]
                        self.assertTrue(all(b > a for a, b in zip(qs, qs[1:])), where)
                        self.assertTrue(all(X >= -1e-9 for X, _ in points), where)
                        if t < 0:
                            X, T = points[0]
                            self.assertLess(abs(X) + abs(T - 2 * math.atan((t + 1) / 3)), 5e-3, where)
                    elif metric_id == "hayward" and mark["view"] == "history":
                        # q = arctan((v - 4)/4) and p the same function of the advanced time v_0 at which the
                        # outgoing ray left the centre: the slice v - r = t starts on the centre, X = 0, at
                        # v = t, runs outward with v growing, and ends at r = 6, where v = t + 6.
                        X, T = points[0]
                        self.assertLess(abs(X) + abs(T - 2 * math.atan((t - 4) / 4)), 2e-4, where)
                        qs = [(T + X) / 2 for X, T in points]
                        self.assertTrue(all(b > a for a, b in zip(qs, qs[1:])), where)
                        self.assertTrue(all(X >= -1e-9 for X, _ in points), where)
                        v_end = 4 + 4 * math.tan(qs[-1])
                        self.assertLess(abs((v_end - t) - 6), 2e-2 * (1 + v_end * v_end / 16), where)
                    elif metric_id == "hosotani_universe":
                        # The one maximum-radius sphere is eta=pi/2 across 0<=chi<=pi.
                        self.assertTrue(all(abs(T - math.pi / 2) < 2e-4 for _, T in points), where)
                        self.assertAlmostEqual(min(X for X, _ in points), 0, delta=2e-4, msg=where)
                        self.assertAlmostEqual(max(X for X, _ in points), math.pi, delta=2e-4, msg=where)
                    elif metric_id == "tilted_universes":
                        # tan p = V and tan q = U are k e^(r/a) and k' e^(-r/a) on a surface of constant eta, so
                        # their product is one number along the slice, positive above the horizon, and the
                        # slice runs from the rim toward the whimper, q falling.
                        products = [math.tan((T - X) / 2) * math.tan((T + X) / 2) for X, T in points]
                        self.assertGreater(min(products), 0, where)
                        self.assertLess(max(products) - min(products), 2e-2 * max(products), where)
                        qs = [(T + X) / 2 for X, T in points]
                        self.assertTrue(all(b < a for a, b in zip(qs, qs[1:])), where)
                    else:
                        self.assertTrue(all(abs(T) < 2e-4 for _, T in points), where)

    def test_a_slice_of_a_view_or_surface_the_embedding_lacks_or_of_a_moved_moment_is_refused(self):
        metrics = build.load_metrics()
        diagrams, conformal, embedding = (build.load_diagrams(metrics), build.load_conformal(metrics),
                                          build.load_embedding(metrics))
        build.check_slices(diagrams, conformal, embedding)
        stray = copy.deepcopy(diagrams)
        stray["schwarzschild"]["systems"]["spherical"][0]["slices"][0]["view"] = "elsewhere"
        with self.assertRaises(build.DataError) as raised:
            build.check_slices(stray, conformal, embedding)
        self.assertIn("'elsewhere'", str(raised.exception))
        moved = copy.deepcopy(embedding)
        moved["kasner"]["views"][1]["surfaces"][1]["time"] = 0.6
        with self.assertRaises(build.DataError) as raised:
            build.check_slices(diagrams, conformal, moved)
        self.assertIn("kasner", str(raised.exception))
        self.assertIn("null_rays.py --slices", str(raised.exception))

    def test_every_slice_class_is_styled_on_the_page_and_in_print(self):
        _, rules = page_rules()
        for cls in (".nr-slice", ".cd-slice", ".pj-slice", ".nr-slice-fill", ".pj-slice-fill", ".mfs-slice-label"):
            for printed in (False, True):
                self.assertTrue(any(any(cls in selector for selector in selectors) and bool(p) == printed
                                    for selectors, _, p in rules), f"{cls}, print {printed}")


class ShadedRegions(unittest.TestCase):
    """While the Einstein static universe's conformal diagram shows Minkowski space, de Sitter
    space or anti-de Sitter space in its strip, the embedding diagram's sphere is shaded where
    that region covers it, as the captain asked on 30 September 2026. The region is Hawking and
    Ellis's at every conformal time the file carries, it is the region the conformal diagram
    draws, and the shade drawn is the part of the sphere between its circles, from the numbers in
    the files alone. _tools/README.md, "The shade of a conformal region", defines the fields."""

    # Hawking and Ellis's regions, the values of chi each spacetime covers at the conformal
    # time eta, or None where it covers none of the moment.
    HAWKING_ELLIS = {
        "minkowski": lambda eta: (0.0, math.pi - abs(eta)) if abs(eta) < math.pi else None,
        "de_sitter": lambda eta: (0.0, math.pi) if abs(eta) < math.pi / 2 else None,
        "anti_de_sitter": lambda eta: (0.0, math.pi / 2),
    }

    def setUp(self):
        self.conformal, self.embedding = conformal_files(), embedding_files()
        self.static = conformal_file("einstein_static")
        self.view = embedding_file("einstein_static")["views"][0]
        self.shades = {sh["view"]: sh for sh in self.view.get("shades", [])}

    @staticmethod
    def across(polygon, T):
        """The least and greatest X at which the line of constant T crosses a convex polygon."""
        xs = []
        for (x0, t0), (x1, t1) in zip(polygon, polygon[1:] + polygon[:1]):
            if min(t0, t1) <= T <= max(t0, t1):
                xs += [x0, x1] if t0 == t1 else [x0 + (T - t0) * (x1 - x0) / (t1 - t0)]
        return (min(xs), max(xs)) if xs else None

    def test_each_region_drawn_inside_the_strip_has_its_shade(self):
        if chosen("einstein_static"):
            views = [v["id"] for v in self.static["views"] if not v.get("system")]
            self.assertEqual(views, ["minkowski", "de_sitter", "anti_de_sitter"])
            self.assertEqual(list(self.shades), views)
        # Only the Einstein static universe's sphere is shaded, and every shade names a view its
        # spacetime's conformal diagram draws.
        for name, data in self.embedding.items():
            conformal = {v["id"] for v in self.conformal.get(name, {}).get("views", [])}
            for view in data["views"]:
                if name != "einstein_static":
                    self.assertNotIn("shades", view, name)
                for shade in view.get("shades", []):
                    self.assertIn(shade["view"], conformal, name)

    @checks("einstein_static")
    def test_the_shaded_region_matches_each_map_at_several_eta(self):
        for vid, shade in self.shades.items():
            times = [T for T, *_ in shade["reach"]]
            self.assertGreaterEqual(len(set(times)), 13, vid)
            self.assertIn(shade["T"], times, vid)
            cover = [L for L in next(v for v in self.static["views"] if v["id"] == vid)["layers"]
                     if L["kind"] == "fill" and L["class"] == "cover"]
            self.assertEqual(len(cover), 1, vid)
            for T, lo, hi in shade["reach"]:
                want = self.HAWKING_ELLIS[vid](T)
                drawn = self.across(cover[0]["points"], T)
                where = f"{vid} at eta = {T}"
                if want is None:
                    self.assertEqual((lo, hi), (None, None), where)
                    self.assertIsNone(drawn, where)
                    continue
                # Hawking and Ellis's region, and the region the conformal diagram draws.
                self.assertAlmostEqual(lo, want[0], delta=1e-7, msg=where)
                self.assertAlmostEqual(hi, want[1], delta=1e-7, msg=where)
                self.assertAlmostEqual(lo, drawn[0], delta=1e-4, msg=where)
                self.assertAlmostEqual(hi, drawn[1], delta=1e-4, msg=where)

    @checks("einstein_static")
    def test_the_shade_drawn_is_the_region_at_the_moment_drawn(self):
        surface = self.view["surfaces"][0]
        for vid, shade in self.shades.items():
            self.assertEqual(shade["T"], 0, vid)
            _, lo, hi = next(r for r in shade["reach"] if r[0] == shade["T"])
            self.assertAlmostEqual(shade["from"], lo, delta=1e-6, msg=vid)
            self.assertAlmostEqual(shade["to"], hi, delta=1e-6, msg=vid)
            piece = next(p for p in surface["pieces"] if p["id"] == shade["piece"])
            xs = [point[0] for point in piece["points"]]
            self.assertIn(shade["from"], xs, vid)
            self.assertIn(shade["to"], xs, vid)
            self.assertEqual(shade["legend"][:2], ["fill", "shade"], vid)
            self.assertTrue(shade["caption"].strip(), vid)
            # The shade takes the place of the sphere's tint: its fills are the shade alone.
            self.assertEqual({L["class"] for L in shade["layers"]}, {"shade"}, vid)
            self.assertTrue(all(L["kind"] == "fill" for L in shade["layers"]), vid)
        # The whole sphere, the sphere but one point and the hemisphere about the pole.
        self.assertEqual((self.shades["minkowski"]["from"], self.shades["minkowski"]["to"]), (0, math.pi))
        self.assertEqual((self.shades["de_sitter"]["from"], self.shades["de_sitter"]["to"]), (0, math.pi))
        self.assertEqual((self.shades["anti_de_sitter"]["from"], self.shades["anti_de_sitter"]["to"]), (0, math.pi / 2))

    @checks("einstein_static")
    def test_the_caption_says_the_moment_shaded(self):
        for vid, shade in self.shades.items():
            self.assertIn("pink, $\\eta = 0$", shade["caption"], vid)

    @checks("einstein_static")
    def test_the_page_turns_the_shade_as_it_was_published(self):
        views = {(v["metric"], v["view"]): v for v in turn_check(self)["views"]}
        v = views[("einstein_static", self.view["id"])]
        self.assertEqual([sh["view"] for sh in v["shades"]], list(self.shades))
        for sh in v["shades"]:
            for cls, (published, drawn) in sh["fills"].items():
                self.assertLessEqual(abs(published - drawn), 3e-3 * v["boxArea"], f"{sh['view']} {cls}")
            self.assertGreater(sh["away"], 0, sh["view"])
            self.assertEqual(sh["bad"], 0, sh["view"])
            self.assertLessEqual(sh["outside"], 1e-9, sh["view"])
        # Taking the shade away gives back the published tint.
        self.assertAlmostEqual(v["unshaded"]["cover"], v["fills"]["cover"][1], places=9)

    def test_the_page_shades_in_the_sites_pink_without_glow(self):
        page = (build.ROOT / "_layouts" / "mfs.html").read_text(encoding="utf-8")
        rule = re.search(r"\.em-shade-fill \{([^}]*)\}", page)
        self.assertIsNotNone(rule)
        self.assertIn("fill: var(--pink-light)", rule.group(1))
        self.assertNotRegex(rule.group(1), r"shadow|filter|transition")
        # The conformal diagram's buttons name their view, and choosing one shows its shade.
        self.assertIn("'data-cd-view': v.id", page)
        self.assertIn("showShade(root, button.getAttribute('data-cd-view'))", page)


class Relations(unittest.TestCase):
    """Every spacetime lists the spacetimes it is related to, each with a sentence on why, as the
    captain asked on 1 October 2026: every link resolves, and a relation written on one
    spacetime is written on the other too, in that one's own words, with the kind that answers
    it. _tools/README.md carries the format."""

    @staticmethod
    def metric(metric_id, related, references=()):
        return {"id": metric_id, "name": metric_id, "short_name": metric_id, "tags": ["t"],
                "related": related, "references": list(references)}

    def pair(self, kind="generalisation", back="special_case", text="Kerr with $a = 0$.", **extra):
        return [self.metric("a", [dict({"id": "b", "kind": kind, "text": text}, **extra)]),
                self.metric("b", [{"id": "a", "kind": back, "text": "The same, set spinning."}])]

    def problems(self, metrics):
        return build.relation_problems(metrics)

    def test_every_spacetime_on_disk_lists_at_least_one_related_spacetime(self):
        for metric in build.load_metrics():
            self.assertTrue(metric.get("related"), metric["id"])

    def test_every_link_on_disk_resolves_and_is_answered(self):
        """Stated here apart from the build's own check, so the two cannot go wrong together."""
        metrics = {m["id"]: m for m in build.load_metrics()}
        links = {(m["id"], entry["id"]) for m in metrics.values() for entry in m["related"]}
        self.assertGreater(len(links), len(metrics))
        for one, other in sorted(links):
            self.assertIn(other, metrics, f"{one}.json lists {other}, which has no metric file")
            self.assertNotEqual(one, other)
            self.assertIn((other, one), links, f"{one}.json lists {other}, and {other}.json does not list {one}")

    def test_each_side_of_a_relation_has_its_own_wording(self):
        metrics = {m["id"]: m for m in build.load_metrics()}
        for metric in metrics.values():
            for entry in metric["related"]:
                back = next(e for e in metrics[entry["id"]]["related"] if e["id"] == metric["id"])
                self.assertNotEqual(entry["text"], back["text"], f"{metric['id']} and {entry['id']}")

    def test_no_relation_states_how_many_spacetimes_there_are(self):
        count = len(build.load_metrics())
        for metric in build.load_metrics():
            for entry in metric["related"]:
                self.assertNotRegex(entry["text"], rf"(?i)\b{count}\b|\b(?:spacetimes|solutions) in (?:this|the) (?:collection|catalogue)\b",
                                    f"{metric['id']}.json: related[{entry['id']}]")

    def test_a_pair_that_answers_itself_passes(self):
        self.assertEqual(self.problems(self.pair()), [])
        self.assertEqual(self.problems(self.pair("family", "family")), [])

    def test_every_kind_has_an_answer_that_answers_back(self):
        for kind, back in build.RELATION_KINDS.items():
            self.assertEqual(build.RELATION_KINDS[back], kind)

    def limit_pair(self, back="limit_source", reached="The limit $m \\to \\infty$.",
                   source="The black hole this one is a limit of."):
        metrics = self.pair("limit", back, text=reached)
        metrics[1]["related"][0]["text"] = source
        return metrics

    def test_a_limit_is_answered_by_its_source_and_by_nothing_else(self):
        """A limit is its own pair of kinds, as the captain asked on 2 October 2026: the spacetime
        reached is a `limit` on the page of the one it is a limit of, which is its `limit_source`."""
        self.assertEqual(build.RELATION_KINDS["limit"], "limit_source")
        self.assertEqual(self.problems(self.limit_pair()), [])
        for back in ("generalisation", "special_case", "family", "limit"):
            found = self.problems(self.limit_pair(back))
            self.assertTrue(any("b.json has to list a as limit_source" in problem for problem in found), back)

    def test_a_limit_that_does_not_say_so_is_refused(self):
        self.assertEqual(self.problems(self.limit_pair(reached="Kerr with no spin.")),
                         ["a.json: related[b] is filed as a limit and does not say so, by the word limit "
                          "or by an arrow such as $m \\to \\infty$"])
        found = self.problems(self.limit_pair(source="The charged black hole."))
        self.assertEqual(len(found), 1)
        self.assertIn("b.json: related[a] is filed as a limit_source and does not say so", found[0])
        for text in ("The limit taken toward its horizon.", "The rod shrunk to a point, $q \\to \\infty$.",
                     "One of the limits Geroch took."):
            self.assertEqual(self.problems(self.limit_pair(reached=text, source=text + " Again.")), [], text)

    def test_every_limit_on_disk_says_so_on_both_pages(self):
        """Stated here apart from the build's own check, so the two cannot go wrong together."""
        limits = [(metric["id"], entry) for metric in build.load_metrics() for entry in metric["related"]
                  if entry["kind"] in ("limit", "limit_source")]
        self.assertGreater(len(limits), 20)
        self.assertEqual(sum(entry["kind"] == "limit" for _, entry in limits) * 2, len(limits))
        for metric_id, entry in limits:
            self.assertRegex(entry["text"], r"\blimits?\b|\\to\b", f"{metric_id}.json: related[{entry['id']}]")

    def test_no_special_case_on_disk_is_the_limit_of_a_boost_or_of_a_horizon(self):
        """The relations filed as special cases until 2 October 2026 that are limits, each by the
        pair of spacetimes: none may go back to `special_case`."""
        kinds = {(metric["id"], entry["id"]): entry["kind"] for metric in build.load_metrics()
                 for entry in metric["related"]}
        for source, limit in (("schwarzschild", "aichelburg_sexl"), ("rn_metric", "bertotti_robinson"),
                              ("majumdar_papapetrou", "bertotti_robinson"),
                              ("reissner_nordstrom_de_sitter", "bertotti_robinson"),
                              ("zipoy_voorhees", "curzon_chazy"), ("zipoy_voorhees", "levi_civita"),
                              ("schwarzschild_de_sitter", "nariai"), ("reissner_nordstrom_de_sitter", "nariai"),
                              ("hayward", "de_sitter"), ("bardeen", "de_sitter"), ("schwarzschild", "kasner"),
                              ("tolman_bondi", "vaidya")):
            self.assertEqual(kinds[source, limit], "limit", f"{source}.json: related[{limit}]")
            self.assertEqual(kinds[limit, source], "limit_source", f"{limit}.json: related[{source}]")

    def test_a_relation_written_one_way_only_is_refused(self):
        metrics = self.pair()
        metrics[1]["related"] = [{"id": "c", "kind": "family", "text": "A cousin."}]
        metrics.append(self.metric("c", [{"id": "b", "kind": "family", "text": "A cousin too."}]))
        self.assertEqual(self.problems(metrics), ["a.json lists b, and b.json does not list a"])

    def test_a_relation_answered_with_the_wrong_kind_is_refused(self):
        found = self.problems(self.pair("generalisation", "generalisation"))
        self.assertEqual(len(found), 2)
        self.assertIn("b.json has to list a as special_case", found[0])

    def test_a_link_to_no_spacetime_is_refused(self):
        metrics = self.pair()
        metrics[0]["related"].append({"id": "ghost", "kind": "family", "text": "Nobody."})
        self.assertEqual(self.problems(metrics), ["a.json: related[ghost] names a spacetime that has no metric file"])

    def test_a_spacetime_listing_itself_or_another_twice_is_refused(self):
        metrics = self.pair()
        metrics[0]["related"].append({"id": "a", "kind": "family", "text": "Itself."})
        metrics[0]["related"].append({"id": "b", "kind": "generalisation", "text": "Again."})
        self.assertEqual(self.problems(metrics), ["a.json lists itself as related", "a.json lists b twice"])

    def test_a_spacetime_with_no_relations_is_refused(self):
        metrics = self.pair()
        del metrics[1]["related"]
        found = self.problems(metrics)
        self.assertIn("b.json lists no related spacetimes", found)

    def test_an_unknown_kind_a_stray_field_and_an_empty_text_are_refused(self):
        self.assertIn("which is not one of", self.problems(self.pair("cousin"))[0])
        self.assertIn("and nothing else", self.problems(self.pair(note="x"))[0])
        self.assertIn("does not say why", self.problems(self.pair(text=" "))[0])

    def test_a_text_past_three_sentences_or_ending_mid_sentence_is_refused(self):
        self.assertIn("runs past 3 sentences", self.problems(self.pair(text="One. Two. Three. Four."))[0])
        self.assertIn("does not end at the end of a sentence", self.problems(self.pair(text="Kerr with no spin"))[0])
        self.assertEqual(self.problems(self.pair(text="He called it \"a counterexample\" [misner1967].")),
                         ["a.json: related[b] cites 'misner1967', which a.json does not list in its references"])

    def test_a_citation_outside_the_references_is_refused(self):
        metrics = self.pair(text="Kerr with $a = 0$ [kerr1963].")
        self.assertIn("cites 'kerr1963'", self.problems(metrics)[0])
        metrics[0]["references"] = ["kerr1963"]
        self.assertEqual(self.problems(metrics), [])

    def test_a_broken_relation_leaves_both_published_files_alone(self):
        before = {path: path.read_text(encoding="utf-8") for path in (build.INDEX_FILE, build.REFERENCES_FILE)}
        broken = build.load_metrics()
        broken[0] = dict(broken[0], related=broken[0]["related"][1:])
        with mock.patch.object(build, "load_metrics", return_value=broken):
            self.assertEqual(build.main([]), 2)
        for path, text in before.items():
            self.assertEqual(path.read_text(encoding="utf-8"), text)

    def test_the_page_lists_them_last_before_the_references(self):
        page = (build.ROOT / "_layouts" / "mfs.html").read_text(encoding="utf-8")
        body = page[page.index("var body = ["):]
        geodesics, related, references = (body.index(mark) for mark in (
            "mfs-section-label\">Geodesics<", "relatedSection(data),", "mfs-section-label\">references<"))
        self.assertLess(geodesics, related)
        self.assertLess(related, references)
        self.assertIn("mfs-section-label\">Related Spacetimes<", page)


@page_source
class SectionOrder(unittest.TestCase):
    """The sections of a spacetime's page stand in the order _tools/README.md gives under "Where the
    diagrams stand", each diagram beside the mathematics it draws: the spacetime diagram and the
    conformal diagram after the coordinates, their conventions and the line element, the embedding diagram after the
    metric, and the curvature figure's Ricci scalar after the Ricci scalar and its Kretschmann scalar
    after the Kretschmann scalar. The contents beside the page and its print copy read the page as
    drawn, so they follow it."""

    ORDER = (
        ("history", 'mfs-section-label">history<'),
        ("coordinates", 'mfs-section-label">coordinates<'),
        ("conventions", "conventionSection(data, coord),"),
        ("line element", 'mfs-section-label">line element<'),
        ("spacetime diagram", "diagramSection(data.id, coord.id)"),
        ("conformal diagram", "conformalSection(data.id, hasCoords ? coord.id : null),"),
        ("metric", 'mfs-section-label">metric<'),
        ("embedding diagram", "embeddingSection(data.id, hasCoords ? coord.id : null),"),
        ("inverse metric", 'mfs-section-label">inverse metric<'),
        ("Christoffel symbols", 'mfs-section-label">Christoffel symbols<'),
        ("Riemann tensor", 'mfs-section-label">Riemann tensor<'),
        ("Ricci tensor", 'mfs-section-label">Ricci tensor<'),
        ("Ricci scalar", 'mfs-section-label">Ricci scalar<'),
        ("Ricci scalar diagram", "curvatureSection(data.id, hasCoords ? coord.id : null, true),"),
        ("Kretschmann scalar", 'mfs-section-label">Kretschmann scalar<'),
        ("Kretschmann scalar diagram", "curvatureSection(data.id, hasCoords ? coord.id : null, false),"),
        ("Einstein tensor", 'mfs-section-label">Einstein tensor<'),
        ("Weyl tensor", 'mfs-section-label">Weyl tensor<'),
        ("geodesics", 'mfs-section-label">Geodesics<'),
        ("related spacetimes", "relatedSection(data),"),
        ("references", 'mfs-section-label">references<'),
    )

    def test_the_page_writes_its_sections_in_order(self):
        page = (build.ROOT / "_layouts" / "mfs.html").read_text(encoding="utf-8")
        body = page[page.index("var body = ["):page.index("window._mfsPrintBody = ")]
        for name, mark in self.ORDER:
            self.assertEqual(body.count(mark), 1, name)
        places = [body.index(mark) for _, mark in self.ORDER]
        self.assertEqual(places, sorted(places))
        # The slice and the shade first drawn are set before the body, so before every diagram.
        self.assertLess(page.index("firstViews(data.id, hasCoords ? coord.id : null);"), page.index("var body = ["))

    def test_the_readme_gives_the_same_order(self):
        readme = (build.ROOT / "_tools" / "README.md").read_text(encoding="utf-8")
        start = readme.index("The sections of a spacetime's page stand in this order")
        block = re.search(r"\n\n((?:    .*\n)+)", readme[start:]).group(1)
        names = [n.strip() for n in " ".join(block.split()).split(",")]
        self.assertEqual(names, [name for name, _ in self.ORDER])


class CurvatureCaptionSplit(unittest.TestCase):
    """A curvature view that draws both scalars stands in two places, and `ricci_sentences` says
    which sentences of its caption go under the Ricci scalar, every other under the Kretschmann
    scalar. The page counts a caption's sentences as build_mfs_data.caption_sentences() does, and
    every such view on disk names its Ricci scalar's sentences."""

    def test_each_place_opens_on_its_own_scalar_with_the_rest_as_written(self):
        opened = 0
        for name, data in curvature_files().items():
            for view in data["views"]:
                if "per_scalar" not in view:
                    continue
                first = build.caption_sentences(view["caption"][0])[0]
                rest = first[first.index(" of ") + 4:]
                kretschmann = view["per_scalar"]["kretschmann"]["opening"]
                ricci = view["per_scalar"]["ricci_scalar"]["opening"]
                self.assertEqual(kretschmann, "The Kretschmann scalar $K$ of " + rest, f"{name} {view['id']}")
                # The Ricci scalar's place stands it beside only what its panel draws, by that curve's
                # name in its legend, as Bardeen's beside Schwarzschild's $R = 0$.
                head = rest.rsplit(", beside ", 1)[0].rstrip(".")
                self.assertTrue(ricci.startswith("The Ricci scalar $R") and " of " + head in ricci, f"{name} {view['id']}")
                legend = {e[1]: e[2] for e in view["per_scalar"]["ricci_scalar"]["legend"]}
                self.assertEqual(", beside " in ricci, "reference" in legend, f"{name} {view['id']}")
                if "reference" in legend and legend["reference"].endswith("$R = 0$"):
                    self.assertIn(legend["reference"], ricci, f"{name} {view['id']}")
                opened += 1
        found(self, opened, "curvature view of both scalars")

    def test_no_sentence_of_a_place_names_the_other_scalar_but_as_a_relation(self):
        # Every sentence of a two-scalar view's caption stands with one scalar: one under the
        # Ricci scalar names no K, and one under the Kretschmann scalar names R only where it states
        # how K stands to R, as K = 5(R^mu_mu)^2/3 does.
        relations = ("K = 5(R^\\mu{}_\\mu)^2/3", "$K$ grows as $(R^\\mu{}_\\mu)^2/3$")
        for name, data in curvature_files().items():
            for view in data["views"]:
                if "ricci_sentences" not in view:
                    continue
                ricci_name = [p["name"] for p, r in zip(view["panels"], build.curvature_ricci_panels(view)) if r][0]
                symbol = "R^\\mu{}_\\mu" if "\\mu" in ricci_name else "R"
                named = {tuple(a) for a in view["ricci_sentences"]}
                for i, paragraph in enumerate(view["caption"]):
                    for j, sentence in enumerate(build.caption_sentences(paragraph)):
                        if (i, j) == (0, 0):
                            continue
                        maths = re.findall(r"\$[^$]*\$", sentence)
                        k = any(re.search(r"(^|[^A-Za-z\\])K([^A-Za-z]|$)", m) for m in maths)
                        r = (any(symbol in m for m in maths) if symbol != "R" else
                             any(re.search(r"(^|[^A-Za-z\\])R([^A-Za-z_]|$)", m) for m in maths))
                        where = f"{name} {view['id']} [{i},{j}]"
                        if (i, j) in named:
                            self.assertFalse(k, where + " names K under the Ricci scalar")
                        elif r:
                            self.assertTrue(any(rel in sentence for rel in relations), where + " names R under K")

    def test_the_page_counts_every_captions_sentences_as_the_build_does(self):
        if shutil.which("node") is None:
            self.skipTest("Node is not installed, so the page's counting cannot be run")
        page = (build.ROOT / "_layouts" / "mfs.html").read_text(encoding="utf-8")
        counting = re.search(r"      function cvSentences\([\s\S]*?\n      \}\n", page)
        self.assertIsNotNone(counting, "no cvSentences()")
        paragraphs = [p for data in curvature_files().values() for view in data["views"] for p in view["caption"]]
        found(self, paragraphs, "curvature caption")
        script = (counting.group(0) + "const texts = JSON.parse(require('fs').readFileSync(0, 'utf8'));"
                  "console.log(JSON.stringify(texts.map(cvSentences)));")
        run = subprocess.run(["node", "-e", script], input=json.dumps(paragraphs), capture_output=True, text=True,
                             check=True)
        for paragraph, sentences in zip(paragraphs, json.loads(run.stdout)):
            self.assertEqual(sentences, build.caption_sentences(paragraph), paragraph[:60])
            self.assertEqual(" ".join(sentences), paragraph)

    def test_every_view_of_both_scalars_names_its_ricci_scalars_sentences(self):
        both = 0
        for name, data in curvature_files().items():
            for view in data["views"]:
                ricci = build.curvature_ricci_panels(view)
                if any(ricci) and not all(ricci):
                    both += 1
                    self.assertTrue("ricci_sentences" in view, f"{name} {view['id']} names none")
                    build.check_ricci_sentences(f"{name} {view['id']}", view)
                    self.assertTrue("per_scalar" in view, f"{name} {view['id']} states no scale of each place")
                    build.check_per_scalar(f"{name} {view['id']}", view)
                    # A comparison that vanishes on the Ricci scalar's panel is named as what it is.
                    zero = any(c["class"] == "reference" and all(p[1] == 0 for p in c["points"])
                               for h in [view] + view.get("movie", {}).get("frames", [])
                               for c in h.get("curves", []) if c["panel"] == "R")
                    names = {e[1]: e[2] for e in view["per_scalar"]["ricci_scalar"]["legend"]}
                    if zero:
                        self.assertTrue(names["reference"].endswith("$R = 0$"), f"{name} {view['id']}")
                else:
                    self.assertFalse("ricci_sentences" in view, f"{name} {view['id']} draws one scalar")
                    self.assertFalse("per_scalar" in view, f"{name} {view['id']} draws one scalar")
        found(self, both, "curvature view of both scalars")

    def test_a_split_that_loses_or_repeats_a_sentence_is_refused(self):
        view = {"kind": "line", "panels": [{"id": "K"}, {"id": "R"}],
                "curves": [{"panel": "K", "scalar": "kretschmann"}, {"panel": "R", "scalar": "ricci_scalar"}],
                "caption": ["The Kretschmann scalar $K$ and the Ricci scalar $R$ of one $K = 1.5$. Two.", "Three."]}
        build.check_ricci_sentences("v", dict(view, ricci_sentences=[[0, 1], [1, 0]]))
        for bad in ([[0, 2]], [[1, 0], [0, 1]], [[0, 1], [0, 1]], [[0]], "all"):
            with self.assertRaises(build.DataError, msg=str(bad)):
                build.check_ricci_sentences("v", dict(view, ricci_sentences=bad))
        with self.assertRaises(build.DataError):
            build.check_ricci_sentences("v", dict(view, curves=view["curves"][:1], ricci_sentences=[]))
        with self.assertRaises(build.DataError):
            build.check_ricci_sentences("v", dict(view, caption=["The Kretschmann scalar $K$ of one. Two."],
                                                  ricci_sentences=[]))

    def test_a_place_that_names_what_it_does_not_draw_is_refused(self):
        view = {"kind": "line", "panels": [{"id": "K"}, {"id": "R"}],
                "curves": [{"panel": "K", "class": "K", "scalar": "kretschmann"},
                           {"panel": "R", "class": "R", "scalar": "ricci_scalar"},
                           {"panel": "K", "class": "reference", "scalar": "kretschmann"}],
                "marks": [{"class": "horizon", "x": 1}]}
        view["caption"] = ["The Kretschmann scalar $K$ and the Ricci scalar $R$ of one. Two."]
        good = {"kretschmann": {"opening": "The Kretschmann scalar $K$ of one.", "clip": "K.",
                                "legend": [["line", "K", "k"], ["line", "reference", "s"], ["line", "horizon", "h"]]},
                "ricci_scalar": {"opening": "The Ricci scalar $R$ of one.", "clip": "R.",
                                 "legend": [["line", "R", "r"], ["line", "horizon", "h"]]}}
        build.check_per_scalar("v", dict(view, per_scalar=good))
        bad = copy.deepcopy(good)
        bad["ricci_scalar"]["legend"].append(["line", "reference", "s"])
        missing = copy.deepcopy(good)
        del missing["kretschmann"]["clip"]
        reworded = copy.deepcopy(good)
        reworded["ricci_scalar"]["opening"] = "The Ricci scalar $R$ of another."
        for places in (bad, missing, reworded, {"kretschmann": good["kretschmann"]}):
            with self.assertRaises(build.DataError):
                build.check_per_scalar("v", dict(view, per_scalar=places))


class Bibliography(unittest.TestCase):
    def setUp(self):
        self.references = read(build.REFERENCES_FILE)
        self.entries = self.references["entries"]

    def test_every_entry_in_the_bib_file_is_served(self):
        source = build.BIB_FILE.read_text(encoding="utf-8")
        self.assertEqual(len(self.entries), source.count("\n@") + source.startswith("@"))

    def test_every_reference_the_metrics_cite_is_reachable_in_this_one_file(self):
        cited = {r for m in build.load_metrics() for r in m.get("references", [])}
        self.assertTrue(cited)
        self.assertEqual(cited - set(self.entries), set())

    def test_entries_carry_a_type_and_the_fields_a_citation_needs(self):
        for key, entry in self.entries.items():
            self.assertTrue(entry["type"], key)
            for field in ("title", "author", "year"):
                self.assertTrue(entry["fields"].get(field), f"{key} has no {field}")

    def test_braces_are_balanced_in_every_value(self):
        for key, entry in self.entries.items():
            for name, value in entry["fields"].items():
                self.assertEqual(value.count("{"), value.count("}"), f"{key}.{name}")

    def test_nested_and_quoted_values_are_read_whole(self):
        parsed = build.parse_bibtex(
            '@book{a, title = {The {ADM} Formalism}, author = "Arnowitt", year = {1962}}'
        )
        self.assertEqual(parsed["a"]["fields"]["title"], "The {ADM} Formalism")
        self.assertEqual(parsed["a"]["fields"]["author"], "Arnowitt")

    def test_a_repeated_key_is_refused(self):
        with self.assertRaises(build.DataError):
            build.parse_bibtex("@article{a, year = {1} }\n@article{a, year = {2} }")

    def test_a_field_left_without_a_value_is_refused_by_name(self):
        with self.assertRaises(build.DataError) as raised:
            build.parse_bibtex("@misc{lonely, url =\n}")
        self.assertIn("lonely", str(raised.exception))


class Citations(unittest.TestCase):
    def test_every_metric_lists_its_references_in_the_order_its_prose_first_cites_them(self):
        """The page numbers a citation by its place in `references`, so [1] is the first one read:
        the history's citations in order, then those of the related spacetimes, which stand
        right after it."""
        for metric in chosen_metrics():
            read = [metric.get("history") or ""] + [entry["text"] for entry in metric.get("related") or []]
            self.assertEqual(metric.get("references", []), build.cited_keys(" ".join(read)), metric["id"])

    def test_a_bracket_inside_mathematics_is_no_citation(self):
        self.assertEqual(build.cited_keys("where $r < \\sqrt[3]{3m/\\Lambda}$ [podolsky2009]"), ["podolsky2009"])

    def test_a_citation_the_bibliography_has_no_entry_for_is_refused(self):
        metric = {"id": "x", "name": "X", "short_name": "X", "tags": ["t"], "references": ["ghost"]}
        with self.assertRaises(build.DataError) as raised:
            build.check_citations([metric], {"kerr1963": {}})
        self.assertIn("x.json", str(raised.exception))
        self.assertIn("ghost", str(raised.exception))

    def test_a_dangling_citation_leaves_both_published_files_alone(self):
        before = {path: path.read_text(encoding="utf-8") for path in (build.INDEX_FILE, build.REFERENCES_FILE)}
        broken = build.load_metrics()
        broken[0] = dict(broken[0], references=["no_such_key"])
        with mock.patch.object(build, "load_metrics", return_value=broken):
            self.assertEqual(build.main([]), 2)
        for path, text in before.items():
            self.assertEqual(path.read_text(encoding="utf-8"), text)


class Dollars(unittest.TestCase):
    def test_a_stray_dollar_is_caught_and_balanced_mathematics_is_not(self):
        for text in ("holes at the points $z = z_i$ of the axis$, harmonic away from them",
                     "a mass $m and a charge",
                     "$$ds^2 = -dt^2$$ with $t$ and $x",
                     "one paragraph with $a$ and $b¶the next with c$"):
            with self.subTest(text):
                self.assertTrue(build.unbalanced_dollar(text))
        for text in ("no mathematics", "$a$ and $b$", "$$ds^2 = -dt^2$$ with $t$", "a price of \\$5",
                     "$a$¶$$b$$¶$c$"):
            with self.subTest(text):
                self.assertFalse(build.unbalanced_dollar(text))

    def test_a_stray_dollar_anywhere_in_a_file_is_named_by_its_place(self):
        data = {"id": "x", "history": "$a$", "coordinates": [{"parameters": [{"description": "the axis$, harmonic"}]}]}
        problems = build.dollar_problems("x.json", data)
        self.assertEqual(len(problems), 1)
        self.assertIn("x.json.coordinates[0].parameters[0].description", problems[0])

    def test_a_stray_dollar_leaves_both_published_files_alone(self):
        before = {path: path.read_text(encoding="utf-8") for path in (build.INDEX_FILE, build.REFERENCES_FILE)}
        with mock.patch.object(build, "unbalanced_dollar", return_value=True), \
                contextlib.redirect_stderr(io.StringIO()) as said:
            self.assertEqual(build.main([]), 2)
        self.assertIn("has an unbalanced $", said.getvalue())
        for path, text in before.items():
            self.assertEqual(path.read_text(encoding="utf-8"), text)


class TextRefusals(unittest.TestCase):
    """The build refuses a text a reader sees that breaks a rule of voice.py, naming its place."""

    @classmethod
    def setUpClass(cls):
        cls.metrics = build.load_metrics()

    def test_a_text_that_breaks_a_rule_is_named_by_its_file_and_its_place(self):
        kerr = copy.deepcopy(next(m for m in self.metrics if m["id"] == "kerr"))
        kerr["description"] += " It is important to note this entry."
        kerr["related"][0]["text"] = "A well-known limit."
        with self.assertRaises(build.DataError) as caught:
            build.check_texts([kerr], {}, {}, {}, {})
        said = str(caught.exception)
        self.assertIn("kerr.json: description has a word for the collection's machinery: 'entry'", said)
        self.assertIn(f"kerr.json: related {kerr['related'][0]['id']} has a compound spelled with a hyphen: "
                      "'well-known'", said)

    def test_a_caption_is_held_to_the_caption_rules_and_a_label_is_not(self):
        data = {"views": [{"id": "v", "label": "This is a label", "caption": ["This is the plane of $t$ and $r$."],
                           "legend": [], "labels": []}]}
        problems = build.text_problems(build.drawing_texts(build.conformal_prose, "x", data))
        self.assertEqual(problems, ["conformal/x.json v.caption[0] has a sentence opening \"This is\" in place of a "
                                    "noun phrase: 'This is'"])

    def test_a_history_may_say_published_and_no_other_field_may(self):
        metric = {"id": "x", "history": "Einstein published it in 1915.", "description": "Einstein published it."}
        problems = build.text_problems(build.metric_texts(metric))
        self.assertEqual(problems, ["x.json: description has a word for the collection's machinery: 'published'"])

    def test_one_spacetime_is_held_by_its_own_file_and_every_text_relating_another_to_it(self):
        metrics = copy.deepcopy(self.metrics)
        kerr = next(m for m in metrics if m["id"] == "kerr")
        other = next(m for m in metrics if m["id"] == kerr["related"][0]["id"])
        entry = next(e for e in other["related"] if e["id"] == "kerr")
        entry["text"] = "A well-known limit."
        far = next(m for m in metrics if m["id"] not in ("kerr", other["id"]) and
                   all(e["id"] != "kerr" for e in m["related"]))
        far["description"] += " Shown above."
        self.assertEqual(build.spacetime_text_problems("kerr", self.metrics), [])
        self.assertEqual(build.spacetime_text_problems("kerr", metrics),
                         [f"{other['id']}.json: related kerr has a compound spelled with a hyphen: 'well-known'"])


class EmbeddingRefusals(unittest.TestCase):
    """The build refuses a spacetime with no embedding diagram unless it is named as drawn without one."""

    def test_the_spacetimes_without_one_are_those_embedding_py_does_not_draw(self):
        tree = ast.parse((build.ROOT / "_tools" / "derivations" / "embedding.py").read_text(encoding="utf-8"))
        not_drawn = next(ast.literal_eval(node.value) for node in tree.body if isinstance(node, ast.Assign)
                         and [getattr(target, "id", None) for target in node.targets] == ["NOT_DRAWN"])
        self.assertEqual(set(build.WITHOUT_EMBEDDING), not_drawn)

    def test_a_missing_diagram_and_one_drawn_without_need_are_named(self):
        metrics = [{"id": "kerr"}, {"id": "lentz"}]
        self.assertEqual(build.embedding_problems(metrics, {"kerr": {}}), [])
        self.assertEqual(build.embedding_problems(metrics, {}), [
            "kerr.json has no embedding diagram, embedding/kerr.json; draw it with "
            "_tools/derivations/embedding.py --metric kerr"])
        self.assertEqual(build.embedding_problems(metrics, {"kerr": {}, "lentz": {}}), [
            "embedding/lentz.json draws a spacetime WITHOUT_EMBEDDING names as drawn without one"])


class TagRefusals(unittest.TestCase):
    """The build refuses tags that disagree with _tools/derivations/metric_tags.json."""

    @classmethod
    def setUpClass(cls):
        cls.metrics = build.load_metrics()
        cls.facts = build.metric_tags.load_facts()

    def problems(self, change):
        metrics = copy.deepcopy(self.metrics)
        change(next(m for m in metrics if m["id"] == "kerr"))
        return build.tag_problems(metrics, self.facts)

    def test_a_tag_the_charts_decide_against_is_named(self):
        self.assertEqual(self.problems(lambda kerr: kerr["tags"].append("static")),
                         ["kerr.json: tags carries 'static', which its charts decide it has not"])

    def test_a_retired_tag_is_named(self):
        self.assertIn("kerr.json: tags carries 'vacuum solution', which was merged into 'vacuum'",
                      self.problems(lambda kerr: kerr["tags"].append("vacuum solution")))

    def test_a_tag_one_spacetime_alone_carries_is_named(self):
        self.assertEqual(self.problems(lambda kerr: kerr["tags"].append("made up keyword")),
                         ["metrics: tags 'made up keyword' is carried by kerr alone"])

    def test_a_chart_changed_since_its_facts_were_computed_is_named(self):
        def change(kerr):
            kerr["coordinates"][0]["line_element"] += " "
        self.assertEqual(self.problems(change), [
            f"kerr.json: the chart {self.metrics_chart()!r} has changed since its facts were computed; "
            f"{build.TAG_RUN}"])

    def metrics_chart(self):
        return next(m for m in self.metrics if m["id"] == "kerr")["coordinates"][0]["id"]


class Aspect(unittest.TestCase):
    """The captain asked on 2 October 2026: "many of the spacetime diagrams have terrible aspect
    ratios. We should try to stick to 1:1 aspect ratios, allowing for up to 1:2 and 2:1 but no
    more than that." The region a view or a figure plots is its box, which the page and the
    application draw at one scale on both axes."""

    @staticmethod
    def file(box, mirror=False, part="systems"):
        return {"metric": "x", part: {"chart": [{"id": "view", "box": box, "mirror": mirror}]}}

    def test_the_aspect_is_the_width_of_the_box_over_its_height(self):
        self.assertEqual(build.view_aspect({"box": [0, 6, -3, 3], "mirror": False}), 1)
        self.assertEqual(build.view_aspect({"box": [1, 4, -3, 3], "mirror": False}), 0.5)
        self.assertEqual(build.view_aspect({"box": [-2, 2, 0, 1]}), 4)

    def test_a_view_through_a_centre_is_twice_as_wide_as_the_half_its_file_holds(self):
        self.assertEqual(build.view_aspect({"box": [0, 3, -3, 3], "mirror": True}), 1)
        self.assertEqual(build.aspect_problems("x", self.file([0, 1, 0, 4], mirror=True)), [])
        self.assertEqual(len(build.aspect_problems("x", self.file([0, 1, 0, 4]))), 1)

    def test_the_two_limits_are_allowed_and_a_box_just_past_either_is_not(self):
        for box in ([0, 1, 0, 2], [0, 2, 0, 1], [0, 1, 0, 1], [0, 1.3776, -1.3776, 1.3776],
                    [0, 2 * math.pi, 0, math.pi], [-0.3, 0.7, -0.5, 0.5]):
            with self.subTest(box=box):
                self.assertEqual(build.aspect_problems("x", self.file(box)), [])
        for box in ([0, 1, 0, 2.001], [0, 2.001, 0, 1], [0, 2, -6, 6], [-2, 2, 0, 1], [0, 1.3776, -25, 25]):
            with self.subTest(box=box):
                self.assertEqual(len(build.aspect_problems("x", self.file(box))), 1)

    def test_a_figure_in_three_dimensions_is_held_to_it_too(self):
        self.assertEqual(build.aspect_problems("x", self.file([-4, 4, -3, 1], part="projections")), [])
        problems = build.aspect_problems("x", self.file([-4.24, 4.24, -3.24, 0.24], part="projections"))
        self.assertEqual(len(problems), 1)

    def test_a_problem_names_the_spacetime_the_chart_the_view_and_the_shape(self):
        tall, = build.aspect_problems("kerr", self.file([0, 2, -6, 6]))
        self.assertIn("diagrams/kerr.json chart/view", tall)
        self.assertIn("1:6.00", tall)
        wide, = build.aspect_problems("kerr", self.file([-2, 2, 0, 1]))
        self.assertIn("4.00:1", wide)

    def test_every_kind_and_curvature_region_uses_the_same_limits(self):
        cases = [
            ("conformal", {"views": [{"id": "edge", "box": [0, 1, 0, 3]}]}),
            ("embedding", {"views": [{"id": "edge", "figure": {"box": [0, 1, 0, 3]}}]}),
            ("curvature", {"views": [{"id": "edge", "kind": "map",
                                     "x": {"from": 0, "to": 1}, "y": {"from": 0, "to": 3}}]}),
            ("curvature", {"views": [{"id": "edge", "kind": "line",
                                     "panels": [{"id": "p", "aspect": 3}]}]}),
            ("curvature", {"views": [{"id": "edge", "kind": "line",
                                     "panels": [{"id": str(i), "aspect": 1} for i in range(3)]}]}),
        ]
        for kind, data in cases:
            with self.subTest(kind=kind, data=data):
                problems = build.aspect_problems("x", data, kind)
                self.assertTrue(problems)
                self.assertIn(f"{kind}/x.json edge", problems[0])
                with self.assertRaises(build.DataError):
                    build.check_aspects({}, **{kind: {"x": data}})
        with self.assertRaisesRegex(build.DataError, "has no aspect measure"):
            build.check_aspects({}, future={"x": {"views": []}})

    def test_padding_only_enlarges_the_window_and_leaves_valid_windows_alone(self):
        for box in ([0, 1, 0, 4], [0, 4, 0, 1], [-1, 1, -1, 1]):
            padded = build.padded_window(box)
            self.assertLessEqual(padded[0], box[0])
            self.assertGreaterEqual(padded[1], box[1])
            self.assertLessEqual(padded[2], box[2])
            self.assertGreaterEqual(padded[3], box[3])
            self.assertTrue(0.5 <= build.view_aspect({"box": padded}) <= 2)
        self.assertEqual(build.padded_window([-1, 1, -1, 1]), [-1, 1, -1, 1])

    def test_an_enlarged_embedding_window_fits_labels_at_its_final_scale(self):
        sys.path.insert(0, str(build.ROOT / "_tools" / "derivations"))
        import embedding
        label = {"at": [0.25, 0], "text": "$" + "r" * 30 + "$", "class": "small",
                 "anchor": "l", "dx": 8, "dy": 0}
        box = embedding.fit_window([-0.3, 0.3, -2, 2], [label])
        x0, y0, x1, y1 = embedding.label_box(label, (box[1] - box[0]) / 560)
        self.assertLessEqual(box[0], x0)
        self.assertLessEqual(x1, box[1])
        self.assertLessEqual(box[2], y0)
        self.assertLessEqual(y1, box[3])
        self.assertTrue(0.5 <= build.view_aspect({"box": box}) <= 2)

    def test_a_diagram_outside_the_limits_stops_the_build_and_leaves_the_files_alone(self):
        diagrams = copy.deepcopy(diagram_files(every=True))
        name = sorted(diagrams)[0]
        view = next(iter(diagrams[name]["systems"].values()))[0]
        view["box"], view["mirror"] = [0, 1, 0, 3], False
        before = {path: path.read_text(encoding="utf-8") for path in (build.INDEX_FILE, build.REFERENCES_FILE)}
        for argv in (["--check"], []):
            with mock.patch.object(build, "load_diagrams", return_value=diagrams), \
                    contextlib.redirect_stderr(io.StringIO()) as said:
                self.assertEqual(build.main(argv), 2)
            self.assertIn(f"diagrams/{name}.json", said.getvalue())
            self.assertIn("outside 1:2 to 2:1", said.getvalue())
        for path, text in before.items():
            self.assertEqual(path.read_text(encoding="utf-8"), text)


if __name__ == "__main__":
    unittest.main()
