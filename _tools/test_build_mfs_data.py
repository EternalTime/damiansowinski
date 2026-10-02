#!/usr/bin/env python3
"""Tests for the MFS data generator.

    python3 -m unittest discover -s _tools
"""

import ast
import collections
import contextlib
import copy
import decimal
import io
import itertools
import json
import math
import re
import shutil
import subprocess
import tempfile
import unicodedata
import unittest
from decimal import ROUND_HALF_EVEN, Decimal
from fractions import Fraction
from pathlib import Path
from unittest import mock

import build_mfs_data as build

# The dashes the prose is not allowed to carry: U+2010 to U+2015, which are the
# hyphen, the non breaking hyphen, the figure dash, the en dash, the em dash and
# the horizontal bar, together with U+2212, the minus sign. The ASCII hyphen at
# U+002D is deliberately absent, because it spells Reissner-Nordstrom, Kerr-Newman,
# Lanczos-van Stockum, anti-de Sitter and pp-wave, where it is part of the name.
# _tools/README.md carries the rule this stands for.
DASHES = "\u2010\u2011\u2012\u2013\u2014\u2015\u2212"

# The fields written in sentences. The LaTeX and structural fields are mathematics
# and are left out on purpose, since a minus sign belongs in them.
PROSE_FIELDS = ("name", "short_name", "description", "history", "convention")

# Words that in this prose can only mean the collection's own machinery or the page
# describing itself, which a reader who came for a spacetime is never told about.
# _tools/README.md carries the rule this stands for.
MACHINERY = (
    r"(?i)\bentr(?:y|ies)\b",
    r"(?i)\b(?:this|the) collection\b",
    r"(?i)\b(?:displayed|published|printed|shown|listed) (?:here|above|below)\b",
    r"(?i)\bthe history describes\b",
    r"`",
)
MACHINERY_OUTSIDE_A_HISTORY = (
    r"(?i)\bpublish\w*",
    r"(?i)\bprint(?:ed|s)?\b",
    r"(?i)(?<!building )\bblocks?\b",
    r"(?i)\bvariants?\b",
)

# "is the distance the metric gives": a noun handed back by what gives it, which reads
# as artificial where "is the metric distance" says the same thing. The sentences and
# the templates say it the compact way. _tools/README.md carries the rule this stands for.
_WORD = r"(?:(?!(?:of|and|or|against|to|in|on|at|by)\b)[\w'-]+ )"
GIVES = rf"(?i)\bthe {_WORD}{{1,2}}the {_WORD}{{0,2}}(?:gives?|yields?|provides?|returns?)\b"
# Nor is the metric said to give anything: an equation is stated as it stands, "where $f = 1$,
# $ds^2 = -c^2dt^2$", never "where $f = 1$ and the metric gives $ds^2 = -c^2dt^2$".
METRIC_GIVES = r"(?i)\bthe metric gives?\b"

# The captain's voice, ~/VOICE.md, as far as a pattern can hold it, over every paragraph of every
# history and convention and every caption, note, restriction band and sentence stated in place
# of a drawing. Mathematics and quotations are left out first: a quotation is its speaker's own
# words. _tools/README.md carries the rules these stand for.
_NOUN = r"(?:(?!(?:of|and|or|against|to|in|on|at|by)\b)[\w'-]+ )"
VOICE = (
    ("an em dash", r"—"),
    ("a dash standing alone as punctuation", r"(?:^|\s)-(?:\s|$)"),
    ("a noun handed back to what gives it", rf"(?i)\bthe {_NOUN}{{1,2}}the {_NOUN}{{0,2}}"
                                            r"(?:gives?|yields?|provides?|returns?|affords?|supplies|supply"
                                            r"|produces?|delivers?|offers?|draws?|declares?)\b"),
    ("a paper, a result or a field cast as the actor",
     r"(?i)\b(?:papers?|articles?|literature|results?|stud(?:y|ies)|textbooks?|numerical relativity)\b "
     r"(?:reports?|shows?|landed|lands|argues?|finds?|reads?|tests?|demonstrates?|tells?|writes?|is written)\b"
     r"|\breads? as\b|\bwere read\b"),
    ("a drawing cast as the actor", r"(?i)\b(?:diagrams?|drawing|figure|caption|section|page)\b "
                                    r"(?:declares?|draws?|shows?|says|tells)\b"),
    ("a staged reveal", r"(?i)(?:^|[.;:,] |\b(?:and|but|so) )what\b[^.;:$]*?\b(?:is|are|was)\b"
                        r"|\bwhat comes out\b"),
    ("a sentence explaining what is not drawn",
     r"(?i)\b(?:is|are) (?:not|never) (?:drawn|shown|plotted)\b|\bnothing (?:is |are )?(?:drawn|shown|plotted)\b"
     r"|\b(?:drawn|shown) instead\b|\bno \w+(?: \w+)? (?:is|are) (?:drawn|shown|given|plotted)\b"),
    ("a thing said to be what it is not", r"(?i)\b(?:is|are) not\b[^.;:,]{1,80}\bbut\b"),
    ("an absolute for emphasis", r"(?i)\bevery single\b|\bat all\b|\bnever fails\b|\beach and every\b|\bwhatsoever\b"
                                 r"|\bactual(?:ly)?\b"),
    ("a machine's tell", r"(?i)\bit is important to note\b|\bit(?:'s| is) worth not(?:ing|e)\b|\bin essence\b"
                         r"|\binteresting question\b|\bdive in\b|\bfascinating\b"
                         r"|(?:^|\. )(?:Crucially|Importantly|Notably|Interestingly),|\bdelve\b|\butili[sz]e\b"
                         r"|\bleverage\b|\bbookkeeping\b|\bfunctionals?\b|\bin today's world\b|\bin modern times\b"
                         r"|\badditionally\b"),
    ("stacked transitions", r"(?i)\b(?:moreover|furthermore),? (?:moreover|furthermore)\b"),
    ("a sentence about the page", r"(?i)\bas (?:noted|shown|mentioned|described|stated|discussed) "
                                  r"(?:above|below|earlier)\b|\bthis (?:section|caption|note|page|entry)\b"),
    ("a slide into the future tense", r"(?i)\bwe will\b|\bwe'll\b"),
    ("a list without its Oxford comma", r"\$M\$, \$M\$ (?:and|or) \$M\$"),
    ("a figure of speech standing in for the claim",
     r"(?i)\bstands? on its own\b|\bhas its say\b|\bkeeps? its secrets\b|\ba question with answers\b"),
)
# A caption, note or restriction band is a compact noun phrase with its values in parentheses, as in
# the captain's own "A spherically symmetric distribution of dust collapsing from rest ($R_0 = 2\,r_s$),
# each point in the diagram a 2-sphere." It never opens "This is", never says "This is the whole of",
# and nothing in it "stands for" anything. _tools/README.md carries the rule these stand for.
CAPTION_VOICE = (
    ("\"This is the whole of\"", r"(?i)\bthis is the whole of\b"),
    ("\"stands for a sphere\"", r"(?i)\bstands? for (?:a|one|one such) (?:\w+ )?sphere\b"),
    ("\"stands for\"", r"(?i)\bstands? for\b"),
    ("a sentence opening \"This is\" in place of a noun phrase", r"^This is\b"),
)
# A hyphen joins two names, a name and a word, or a designation; beyond those it is part of a
# spelling only in these terms, and an ordinary compound is rewritten without it.
HYPHENATED_TERMS = {"anti-de", "anti-trapped", "plane-fronted", "pp-wave", "scalar-tensor"}

# The templates and pages whose words reach a reader, beside the generated files, and the
# data the site hands to agents.
TEMPLATES = ("_layouts/*.html", "_includes/*.html", "_includes/*.txt", "MFS/*.markdown", "*.markdown", "llms*.txt",
             "data/*.json")


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def prose(metric):
    """Yield every sentence carrying field of a metric, each with the name of its place."""
    for field in PROSE_FIELDS:
        if metric.get(field):
            yield field, metric[field]
    for entry in metric.get("related") or []:
        if isinstance(entry, dict) and entry.get("text"):
            yield f"related[{entry.get('id')}].text", entry["text"]
    for position, system in enumerate(metric.get("coordinates") or []):
        where = f"coordinates[{system.get('id') or position}]"
        if system.get("name"):
            yield f"{where}.name", system["name"]
        if system.get("convention"):
            yield f"{where}.convention", system["convention"]
        for parameter in system.get("parameters") or []:
            if parameter.get("description"):
                yield f"{where}.parameters[{parameter.get('symbol')}].description", parameter["description"]


def load_folder(files):
    """Load metrics from a scratch folder holding `files`, a map of file name to metric."""
    with tempfile.TemporaryDirectory() as folder:
        for filename, metric in files.items():
            (Path(folder) / filename).write_text(json.dumps(metric), encoding="utf-8")
        with mock.patch.object(build, "METRICS_DIR", Path(folder)):
            return build.load_metrics()


def diagram_files():
    return {p.stem: read(p) for p in sorted(build.DIAGRAMS_DIR.glob("*.json"))}


def conformal_files():
    return {p.stem: read(p) for p in sorted(build.CONFORMAL_DIR.glob("*.json"))}


def conformal_prose(name, conformal):
    """Yield every field of a conformal diagram file a reader sees, each with its place."""
    for view in conformal["views"]:
        where = f"conformal/{name}.json {view['id']}"
        yield f"{where}.label", view["label"]
        for position, paragraph in enumerate(view["caption"]):
            yield f"{where}.caption[{position}]", paragraph
        for field in ("restriction", "settings", "input"):
            if view.get(field):
                yield f"{where}.{field}", view[field]
        for position, (_, _, text) in enumerate(view["legend"]):
            yield f"{where}.legend[{position}]", text
        for position, label in enumerate(view["labels"]):
            yield f"{where}.labels[{position}]", label["text"]
        for position, mark in enumerate(view.get("slices", [])):
            yield f"{where}.slices[{position}].label", mark["label"]


def embedding_files():
    return {p.stem: read(p) for p in sorted(build.EMBEDDING_DIR.glob("*.json"))}


def grid_nodes(piece):
    """The nodes (X, Y, Z) of a grid piece, row i of u and column j of v, a polar grid's X and Y
    being u cos v and u sin v, as _tools/README.md defines them."""
    grid = piece["grid"]
    if grid["frame"] == "ellipses":
        # Row i the ellipse (a_i cos v, b_i sin v) at the height z_i.
        return [[(a * math.cos(v), b * math.sin(v), z) for v in grid["v"]]
                for a, b, z in zip(grid["a"], grid["b"], grid["z"])]
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


def embedding_prose(name, data):
    """Yield every field of an embedding diagram file a reader sees, each with its place."""
    for position, text in enumerate(data.get("stops", [])):
        yield f"embedding/{name}.json stops[{position}]", text
    for view in data["views"]:
        where = f"embedding/{name}.json {view['id']}"
        yield f"{where}.label", view["label"]
        yield f"{where}.unit", view["unit"]
        for position, paragraph in enumerate(view["caption"]):
            yield f"{where}.caption[{position}]", paragraph
        for field in ("settings", "input", "height"):
            if view.get(field):
                yield f"{where}.{field}", view[field]
        for position, text in enumerate(view.get("stops", [])):
            yield f"{where}.stops[{position}]", text
        for position, (_, _, text) in enumerate(view["figure"]["legend"]):
            yield f"{where}.figure.legend[{position}]", text
        for position, label in enumerate(view["figure"]["labels"]):
            yield f"{where}.figure.labels[{position}]", label["text"]
        for position, shade in enumerate(view.get("shades", [])):
            yield f"{where}.shades[{position}].legend", shade["legend"][2]
            yield f"{where}.shades[{position}].caption", shade["caption"]
        frames = view.get("movie", {}).get("frames", [])
        if frames:
            yield f"{where}.movie.variable", view["movie"]["variable"]
        for number, surface in [(f"surfaces[{k}]", x) for k, x in enumerate(view["surfaces"])] + \
                [(f"movie.frames[{k}]", x) for k, x in enumerate(frames)]:
            at = f"{where}.{number}"
            if surface.get("label"):
                yield f"{at}.label", surface["label"]
            for k, curve in enumerate(surface.get("curves", [])):
                if curve.get("label"):
                    yield f"{at}.curves[{k}].label", curve["label"]
            for piece in surface["pieces"]:
                for end in ("start", "end", "edge"):
                    if piece.get(end, {}).get("text"):
                        yield f"{at}.{piece['id']}.{end}", piece[end]["text"]
            for ring in surface["rings"]:
                if ring.get("label"):
                    yield f"{at}.rings[{ring['piece']} {ring['x']}]", ring["label"]


def diagram_prose(name, diagram):
    """Yield every sentence carrying field of a diagram file, each with the name of its place:
    its flat views and its figures in three dimensions."""
    for system, views in diagram["systems"].items():
        for view in views:
            where = f"diagrams/{name}.json {system}/{view['id']}"
            yield f"{where}.label", view["label"]
            for position, family in enumerate(view["families"]):
                yield f"{where}.families[{position}]", family
            for position, paragraph in enumerate(view["caption"]):
                yield f"{where}.caption[{position}]", paragraph
            if view.get("input"):
                yield f"{where}.input", view["input"]
            if view.get("cone"):
                yield f"{where}.cone", view["cone"]
            for position, marker in enumerate(view["markers"]):
                if marker.get("legend"):
                    yield f"{where}.markers[{position}].legend", marker["legend"]
            for position, mark in enumerate(view.get("slices", [])):
                yield f"{where}.slices[{position}].label", mark["label"]
    for system, figures in diagram.get("projections", {}).items():
        for figure in figures:
            where = f"diagrams/{name}.json {system}/{figure['id']}"
            yield f"{where}.label", figure["label"]
            for position, paragraph in enumerate(figure["caption"]):
                yield f"{where}.caption[{position}]", paragraph
            for field in ("settings", "input"):
                if figure.get(field):
                    yield f"{where}.{field}", figure[field]
            for position, (_, _, text) in enumerate(figure["legend"]):
                yield f"{where}.legend[{position}]", text
            for position, label in enumerate(figure["labels"]):
                yield f"{where}.labels[{position}]", label["text"]
            for position, mark in enumerate(figure.get("slices", [])):
                yield f"{where}.slices[{position}].label", mark["label"]



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


class PublishedFilesAreCurrent(unittest.TestCase):
    def test_check_mode_passes(self):
        self.assertEqual(build.main(["--check"]), 0)


class Index(unittest.TestCase):
    def setUp(self):
        self.metrics = build.load_metrics()
        self.index = read(build.INDEX_FILE)

    def test_one_entry_per_file_and_no_more(self):
        self.assertEqual(
            [e["id"] for e in self.index],
            [m["id"] for m in self.metrics],
        )
        self.assertEqual(
            sorted(e["id"] for e in self.index),
            sorted(p.stem for p in build.METRICS_DIR.glob("*.json")),
        )

    def test_entry_carries_the_search_fields_and_a_stamp(self):
        diagrams, conformal, embedding = diagram_files(), conformal_files(), embedding_files()
        for entry, metric in zip(self.index, self.metrics):
            expected = ({"id", "name", "tags", "version"} | ({"diagrams"} if metric["id"] in diagrams else set())
                        | ({"conformal"} if metric["id"] in conformal else set())
                        | ({"embedding"} if metric["id"] in embedding else set()))
            self.assertEqual(set(entry), expected)
            self.assertEqual(entry["name"], metric["short_name"])
            self.assertEqual(entry["tags"], metric["tags"])
            self.assertEqual(entry["version"], build.content_version(metric))

    def test_a_spacetime_with_diagrams_carries_their_stamp(self):
        diagrams = diagram_files()
        self.assertTrue(diagrams, "no diagram file was found, so nothing was checked")
        stamped = {e["id"]: e["diagrams"] for e in self.index if "diagrams" in e}
        self.assertEqual(set(stamped), set(diagrams))
        for metric_id, diagram in diagrams.items():
            self.assertEqual(stamped[metric_id], build.content_version(diagram))

    def test_a_spacetime_with_a_conformal_diagram_carries_its_stamp(self):
        conformal = conformal_files()
        self.assertTrue(conformal, "no conformal diagram file was found, so nothing was checked")
        stamped = {e["id"]: e["conformal"] for e in self.index if "conformal" in e}
        self.assertEqual(set(stamped), set(conformal))
        for metric_id, data in conformal.items():
            self.assertEqual(stamped[metric_id], build.content_version(data))

    def test_a_spacetime_with_an_embedding_diagram_carries_its_stamp(self):
        embedding = embedding_files()
        self.assertTrue(embedding, "no embedding diagram file was found, so nothing was checked")
        stamped = {e["id"]: e["embedding"] for e in self.index if "embedding" in e}
        self.assertEqual(set(stamped), set(embedding))
        for metric_id, data in embedding.items():
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
                  for m in build.load_metrics() for field, value in prose(m)]
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

    def test_no_metric_on_disk_spells_a_character_as_an_escape(self):
        # A \u escaped twice in the JSON reaches the page as the six characters of the
        # escape, as the î of Lemaître did in the Tolman-Bondi conventions, since neither
        # TeX nor the page reads it. The character itself is what the files carry.
        for metric in build.load_metrics():
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
        self.assertTrue(views)

    def test_no_prose_names_the_collection_or_the_page(self):
        """A reader is told about the spacetime, never about the collection that holds it."""
        fields = [(f"{m['id']}.json: {field}", field, value)
                  for m in build.load_metrics() for field, value in prose(m)]
        fields += [(field, field, value) for name, diagram in diagram_files().items()
                   for field, value in diagram_prose(name, diagram)]
        fields += [(field, field, value) for name, data in conformal_files().items()
                   for field, value in conformal_prose(name, data)]
        fields += [(field, field, value) for name, data in embedding_files().items()
                   for field, value in embedding_prose(name, data)]
        self.assertTrue(fields, "no prose was read, so nothing was checked")
        for where, field, value in fields:
            for pattern in MACHINERY:
                self.assertNotRegex(value, pattern, where)
            if field != "history":
                # A history tells of papers that were published and printed.
                for pattern in MACHINERY_OUTSIDE_A_HISTORY:
                    self.assertNotRegex(value, pattern, where)

    def test_no_prose_or_template_hands_a_noun_back_to_what_gives_it(self):
        """A reader is told the metric distance, never the distance the metric gives."""
        fields = [(f"{m['id']}.json: {field}", value)
                  for m in build.load_metrics() for field, value in prose(m)]
        fields += [(field, value) for name, diagram in diagram_files().items()
                   for field, value in diagram_prose(name, diagram)]
        fields += [(field, value) for name, data in conformal_files().items()
                   for field, value in conformal_prose(name, data)]
        fields += [(field, value) for name, data in embedding_files().items()
                   for field, value in embedding_prose(name, data)]
        templates = sorted({p for pattern in TEMPLATES for p in build.ROOT.glob(pattern)})
        self.assertTrue(templates, "no template was read, so none was checked")
        fields += [(str(p.relative_to(build.ROOT)), " ".join(p.read_text(encoding="utf-8").split()))
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
        self.assertTrue(fields, "no prose was read, so nothing was checked")
        for where, value in fields:
            for character in value:
                if character in DASHES:
                    self.fail(
                        f"{where} carries "
                        f"U+{ord(character):04X} {unicodedata.name(character)}. "
                        "Write a full stop, a semicolon or a comma instead, "
                        "or rewrite the sentence so it does not want the break."
                    )


def voiced_prose():
    """Yield every text held to the captain's voice, each with its place: the paragraphs of every
    history and convention, and every caption, note, restriction band and sentence stated in place
    of a drawing, and every sentence on why a spacetime is related to another. Labels and legends
    name things and are not sentences."""
    for metric in build.load_metrics():
        for number, paragraph in build.prose_paragraphs(metric.get("history") or ""):
            yield f"{metric['id']}.json: history paragraph {number}", paragraph
        for chart_id, convention in build.chart_conventions(metric):
            yield f"{metric['id']}.json: convention of {chart_id}", convention
        for entry in metric.get("related") or []:
            yield f"{metric['id']}.json: related {entry['id']}", entry["text"]
    sentences = re.compile(r"\.(?:caption\[\d+\]|input|settings|height|restriction|stops\[\d+\]|start|end|edge)$")
    fields = [field for name, diagram in diagram_files().items() for field in diagram_prose(name, diagram)]
    fields += [field for name, data in conformal_files().items() for field in conformal_prose(name, data)]
    fields += [field for name, data in embedding_files().items() for field in embedding_prose(name, data)]
    for where, value in fields:
        if sentences.search(where) or " stops[" in where:
            yield where, value


def without_mathematics_or_quotations(text):
    return re.sub(r'"[^"]*"', '"Q"', re.sub(r"\$\$.+?\$\$|\$[^$]+\$", "$M$", text, flags=re.DOTALL))


def opens_with_a_noun_phrase(caption):
    """Whether a caption opens with a noun phrase naming what is drawn, as "The plane of $t$ and
    $r$ ($\\theta = \\pi/2$, $\\phi = 0$)" does: it starts with a capital, never with "This", and
    its head, up to the first comma, parenthesis, stop or clause of its own, holds no "is" or "are"."""
    text = without_mathematics_or_quotations(caption)
    head = re.split(r"[,(.:;]|\b(?:which|where|when|that|as|whose|who|so|since)\b", text, maxsplit=1)[0]
    return (bool(re.match(r"[A-Z]", text)) and not text.startswith("This ")
            and not re.search(r"\b(?:is|are|was|were)\b", head))


class Voice(unittest.TestCase):
    """Every history, convention and caption keeps to the rules of the captain's voice that a
    pattern can hold."""

    def test_every_voiced_text_keeps_to_the_rules(self):
        fields = list(voiced_prose())
        self.assertGreater(len(fields), 600, "the voiced prose was not all read")
        for where, value in fields:
            text = without_mathematics_or_quotations(value)
            for rule, pattern in VOICE:
                found = re.search(pattern, text)
                self.assertIsNone(found, f"{where} has {rule}: {found and text[max(0, found.start() - 40):found.end() + 40]!r}")

    def test_every_caption_is_a_compact_noun_phrase(self):
        fields = [(where, value) for where, value in voiced_prose() if ".json: " not in where]
        self.assertGreater(len(fields), 300, "the captions were not all read")
        for where, value in fields:
            text = without_mathematics_or_quotations(value)
            for rule, pattern in CAPTION_VOICE:
                found = re.search(pattern, text)
                self.assertIsNone(found, f"{where} has {rule}: {found and text[max(0, found.start() - 40):found.end() + 40]!r}")

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

    def test_every_hyphen_joins_names_or_an_established_term(self):
        for where, value in voiced_prose():
            text = without_mathematics_or_quotations(value)
            for word in re.findall(r"[\w'’]+(?:-[\w'’]+)+", re.sub(r"\[[^\]]*\]", " ", text)):
                parts = word.split("-")
                named = any(part[0].isupper() or part[0].isdigit() for part in parts)
                self.assertTrue(named or word.lower() in HYPHENATED_TERMS,
                                f"{where} spells {word!r} with a hyphen; rewrite the compound without it")

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

    CONTRAST = (r"(?i)\brather than\b|\binstead\b"
                r"|\bnot (?:a|an|the|just|only|merely|simply)\b[^.;:]{0,80}?\bbut\b"
                r"|\b(?:is|are|was|were)\b[^.;:,]{1,60},\s*not (?!long\b|yet\b)\w"
                r"|\bless\b[^.;:,]{1,40},\s*more\b")
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
        fields = [(f"{m['id']}.json: {field}", value) for m in build.load_metrics() for field, value in prose(m)]
        fields += [field for name, diagram in diagram_files().items() for field in diagram_prose(name, diagram)]
        fields += [field for name, data in conformal_files().items() for field in conformal_prose(name, data)]
        fields += [field for name, data in embedding_files().items() for field in embedding_prose(name, data)]
        fields += [("_layouts/mfs.html", self.layout_text()), ("MFS/index.markdown",
                    (build.ROOT / "MFS" / "index.markdown").read_text(encoding="utf-8"))]
        return fields

    def caught(self, text):
        return [m.group(0) for m in re.finditer(self.CONTRAST, without_mathematics_or_quotations(text))]

    def test_no_text_a_reader_sees_defines_a_thing_by_what_it_is_not(self):
        fields = self.texts()
        self.assertGreater(len(fields), 1000, "the texts were not all read")
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
        for text in ("Cornelius Lanczos wrote down the rotating dust cylinder in 1924, not long after Einstein's field "
                     "equations had settled into their final form.",
                     "It was a clean result, not yet joined to any other solution.",
                     "the singularity theorems no longer apply", "each point in the diagram a single event.",
                     "The quote \"rather than\" is his own.", "where $a \\neq b$, not $a = b$"):
            with self.subTest(text):
                self.assertEqual(self.caught(text), [])


def strings(value, where):
    """Yield every string held anywhere in a parsed JSON value, each with its place."""
    if isinstance(value, str):
        yield where, value
    elif isinstance(value, dict):
        for key, item in value.items():
            yield from strings(item, f"{where}.{key}")
    elif isinstance(value, list):
        for position, item in enumerate(value):
            yield from strings(item, f"{where}[{position}]")


class NoCost(unittest.TestCase):
    """No text a reader of the spacetimes page sees says what physics requires in the words of
    money, as the captain asked on 29 September 2026: a shortcut that stays open requires exotic
    matter, and no field of any spacetime or diagram, nor the page, says cost, costs, costly or
    costing."""

    COST = r"(?i)\bcost(?:s|ly|ing)?\b"

    def texts(self):
        fields = Contrast().texts()
        for folder in (build.METRICS_DIR, build.DIAGRAMS_DIR, build.CONFORMAL_DIR, build.EMBEDDING_DIR):
            for path in sorted(folder.glob("*.json")):
                fields += strings(read(path), path.name)
        return fields

    def test_no_text_a_reader_sees_says_cost(self):
        fields = self.texts()
        self.assertGreater(len(fields), 10000, "the texts were not all read")
        for where, value in fields:
            found = re.search(self.COST, value)
            self.assertIsNone(found, f"{where} has {found and found.group(0)!r}; state what is required, "
                                     "what grows or what was spent")

    def test_the_rule_catches_every_form_and_leaves_other_words_alone(self):
        for text in ("a shortcut that stays open costs exotic matter", "at the cost of seven years of work",
                     "what the trip would cost", "a costly bubble", "Costing the trip in fuel"):
            with self.subTest(text):
                self.assertRegex(text, self.COST)
        for text in ("$\\cos t$", "the coast of the Baltic", "a cosmic string", "Costa Rica", "accosted",
                     "a shortcut that stays open requires exotic matter"):
            with self.subTest(text):
                self.assertNotRegex(text, self.COST)


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

    def test_every_history_on_disk_keeps_its_shape(self):
        metrics = build.load_metrics()
        self.assertTrue(all(m.get("history") for m in metrics), "a metric has no history")
        self.assertIsNone(build.check_prose_shape(metrics))

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
    WORDS = (".mfs-nr-figure", ".nr-tick", ".nr-name", ".nr-ref", ".mfs-slice-label", ".mfs-nr-readout",
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
        turn = (build.ROOT / "MFS" / "assets" / "turn.js").read_text(encoding="utf-8")
        found = {
            "conformal.py": re.search(r"^CD_LABEL_SIZE = \{([^}]*)\}", conformal, re.M).group(1),
            "embedding.py": re.search(r"^LAB = \{([^}]*)\}", embedding, re.M).group(1),
            "turn.js": re.search(r"var LAB = \{([^}]*)\}", turn).group(1),
        }
        for name, sizes in found.items():
            with self.subTest(name):
                self.assertEqual({int(v) for v in re.findall(r":\s*(\d+)", sizes)}, {page})

    def test_a_wide_letter_takes_its_width_in_a_label_s_box(self):
        # The margin of 1.5 em is kept from the box cdLabelSize() gives a label, and slices.py
        # places labels by the same box, so both give an m, an M and a W, which MathJax sets
        # about an em wide, 0.5 em over the 0.7 em of a letter: "$4m$" stood 1.4 em from the
        # edge of Majumdar-Papapetrou's embedding diagram while they did not.
        slices = (build.ROOT / "_tools" / "derivations" / "slices.py").read_text(encoding="utf-8")
        self.assertEqual(re.search(r'^WIDE = "(\w+)"$', slices, re.M).group(1),
                         re.search(r"var CD_WIDE = /\[(\w+)\]/g;", self.page).group(1))
        self.assertEqual(re.search(r'^WIDE = "(\w+)"$', slices, re.M).group(1), "mMW")
        self.assertIn("0.7 * len(math.replace(\" \", \"\")) + 0.5 * sum(math.count(c) for c in WIDE)", slices)
        self.assertIn("0.7 * math.replace(/ /g, '').length + 0.5 * (math.match(CD_WIDE) || []).length", self.page)


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

    def test_every_convention_on_disk_keeps_the_rule(self):
        for metric in build.load_metrics():
            with self.subTest(metric["id"]):
                self.assertEqual(build.convention_problems(metric), [])

    def test_every_chart_on_disk_carries_its_own_convention(self):
        charts = [(m["id"], c) for m in build.load_metrics() for c in m["coordinates"]]
        self.assertGreater(len(charts), 40, "the charts were not all read")
        for metric_id, chart in charts:
            with self.subTest(f"{metric_id}/{chart['id']}"):
                self.assertTrue(chart.get("convention", "").strip())

    def test_no_convention_on_disk_runs_past_five_sentences_or_says_what_the_captain_named(self):
        for metric in build.load_metrics():
            for chart_id, text in build.chart_conventions(metric):
                with self.subTest(f"{metric['id']}/{chart_id}"):
                    self.assertLessEqual(len(build.sentences(text)), 5)
                    for phrase in ("cost", "nowhere does the geometry break", "standing objection",
                                   "the Riemann tensor with the traces removed", "slots"):
                        self.assertNotIn(phrase, text.lower())

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
        self.assertTrue(self.diagrams, "no diagram file was found, so nothing was checked")

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
        self.assertEqual(set(build.load_diagrams(self.metrics)), set(self.diagrams))
        self.assertTrue(views)

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
        metric_id = sorted(self.diagrams)[0]
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
        self.assertTrue(views)

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
        metric_id = sorted(self.diagrams)[0]
        system_id = next(iter(self.diagrams[metric_id]["systems"]))
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
        self.assertTrue(self.figures, "no figure was found, so nothing was checked")

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


class StatedHorizons(unittest.TestCase):
    """Every horizon radius a Kerr or Kerr-Newman text prints is the root of Delta = r^2 - 2GMr/c^2
    + a^2 + r_Q^2 at the spin and charge its diagrams are drawn for, rounded once from that root,
    so that no value rounded to more places first, as 1.6245 for r_+ = 1.62449..., can be rounded
    again the wrong way."""

    PARAMETER = r"\$({}) = (?:\\frac\{{(\d+)\}}\{{(\d+)\}}|(\d+))\$"

    def parameters(self, name):
        """The values every diagram view of `name` is drawn for, which must agree."""
        stated = set()
        for views in diagram_files()[name]["systems"].values():
            for view in views:
                values = {m[1]: Fraction(int(m[2]), int(m[3])) if m[2] else Fraction(int(m[4]))
                          for m in re.finditer(self.PARAMETER.format("G|M|a|r_Q"), view["settings"])}
                stated.add(tuple(sorted(values.items())))
        self.assertEqual(len(stated), 1, f"{name}: its diagrams are drawn for different values: {stated}")
        values = dict(stated.pop())
        self.assertEqual((values["G"], values["M"]), (1, 1), name)
        return values["a"], values.get("r_Q", Fraction(0))

    def test_every_printed_horizon_radius_is_the_root_rounded_once(self):
        texts = {"kerr": [], "kerr_newman": []}
        for name in texts:
            texts[name] += diagram_prose(name, diagram_files()[name])
            texts[name] += conformal_prose(name, conformal_files()[name])
            texts[name] += embedding_prose(name, embedding_files()[name])
        found = []
        for name, fields in texts.items():
            a, r_Q = self.parameters(name)
            with decimal.localcontext(prec=50):
                root = (1 - a * a - r_Q * r_Q)
                root = (Decimal(root.numerator) / Decimal(root.denominator)).sqrt()
                radius = {"+": 1 + root, "-": 1 - root}
            for field, text in fields:
                for sign, printed in re.findall(r"r_([+-]) = (\d+\.\d+)", text):
                    places = Decimal(1).scaleb(-len(printed.split(".")[1]))
                    self.assertEqual(printed, str(radius[sign].quantize(places, ROUND_HALF_EVEN)),
                                     f"{field} prints r_{sign} = {printed}, where r_{sign} = {radius[sign]:.12f}")
                    found.append((name, field, sign))
        self.assertTrue({name for name, _, _ in found} == set(texts), f"a spacetime prints no horizon: {found}")
        # Kerr-Newman's r_+ = 1.62449... GM/c^2 is printed in the caption of the principal rays from
        # above and of the light cones on the equator, the two that once printed 1.625.
        for view in ("boyer_lindquist/above", "boyer_lindquist/dragging"):
            self.assertTrue(any(name == "kerr_newman" and f" {view}.caption" in field and sign == "+"
                                for name, field, sign in found), view)


class SzekeresAxis(unittest.TestCase):
    """The marked ray on each half of the axis of Szekeres's cloud, from the numbers in the file:
    inside the cloud it obeys c dt = (dR/dr +- R S'/S) dr, the upper sign toward theta = 0, with
    the declared R^(3/2) = r^(3/2) - (3/2) sqrt(2M) ct, 2M = r^3(5 - 3r^2)/4 and S'/S = 2r(1 - r^2),
    and outside it stays on R = r_s = 1/2, the horizon of the exterior."""

    @staticmethod
    def slope(r, t, sign):
        root = math.sqrt(5 - 3 * r * r)
        u = 1 - 3 * root * t / 4
        R, dR = r * u ** (2 / 3), u ** (2 / 3) + 3 * r * r * t / (2 * u ** (1 / 3) * root)
        return dR + sign * R * 2 * r * (1 - r * r)

    def ray(self, sign, r_to, steps=2000):
        """ct at r_to along the ray that crosses the surface r = 1 where R = 1/2, by Runge and Kutta."""
        r, t = 1.0, (1 - 0.5 ** 1.5) / (1.5 * math.sqrt(0.5))
        h = (r_to - r) / steps
        for _ in range(steps):
            k1 = self.slope(r, t, sign)
            k2 = self.slope(r + h / 2, t + h * k1 / 2, sign)
            k3 = self.slope(r + h / 2, t + h * k2 / 2, sign)
            k4 = self.slope(r + h, t + h * k3, sign)
            r, t = r + h, t + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6
        return t

    def test_the_marked_ray_is_null_in_the_cloud_and_the_horizon_outside(self):
        views = {view["id"]: view for view in diagram_files()["szekeres"]["systems"]["axisymmetric"]}
        self.assertEqual(set(views), {"north", "south"})
        for half, sign in (("north", 1), ("south", -1)):
            view = views[half]
            x0, x1, y0, y1 = view["box"]
            line, = next(m for m in view["markers"] if m["kind"] == "event")["lines"]
            points = [(x0 + (x1 - x0) * x, y0 + (y1 - y0) * y) for x, y in line]
            inside = [(r, t) for r, t in points if r < 1]
            outside = [(r, t) for r, t in points if r >= 1]
            self.assertGreater(len(inside), 3, half)
            self.assertGreater(len(outside), 1, half)
            for r, t in inside:
                self.assertAlmostEqual(t, self.ray(sign, r), delta=3e-4, msg=f"{half}, r = {r}")
            for r, t in outside:
                self.assertAlmostEqual((r ** 1.5 - 1.5 * math.sqrt(0.5) * t) ** (2 / 3), 0.5, delta=3e-4, msg=f"{half}, r = {r}")
            # The times the caption states: the ray leaves the centre at -0.81 toward theta = 0 and at -0.12 toward pi.
            self.assertEqual(f"{self.ray(sign, 1e-9):.2f}", {"north": "-0.81", "south": "-0.12"}[half])
            self.assertIn({"north": "$ct = -0.81\\,r_b$", "south": "$ct = -0.12\\,r_b$"}[half], " ".join(view["caption"]))


class StringWaveRing(unittest.TestCase):
    """The ring of free particles on the cone of the travelling wave on a string, from the numbers
    in the file: at b = 1/2 and ell = 1 the isotropic chart's geodesics across the string are

        x'' = -rho A'' + (x (x'^2 - y'^2) + 2 y x' y')/(2 rho^2),
        y'' =            (y (y'^2 - x'^2) + 2 x x' y')/(2 rho^2),

    with a prime for d/du, rho^2 = x^2 + y^2 and the declared pulse A = exp(-4u^2)/2, and the
    particle at (x, y) stands on the cone at r = 2 sqrt(rho), which the surface of half angle 30
    degrees draws at the distance r/2 from its axis and the height sqrt(3) r/2."""

    START, STEP = -4.0, 0.001

    @staticmethod
    def rate(u, w):
        x, y, vx, vy = w
        rho2 = x * x + y * y
        pulse = (32 * u * u - 4) * math.exp(-4 * u * u)
        return (vx, vy,
                -math.sqrt(rho2) * pulse + (x * (vx * vx - vy * vy) + 2 * y * vx * vy) / (2 * rho2),
                (y * (vy * vy - vx * vx) + 2 * x * vx * vy) / (2 * rho2))

    def run_to(self, angle, times):
        """The particle that starts at rest at the angle given on rho = 1, at each of the times, by
        Runge and Kutta."""
        w, u, out = (math.cos(angle), math.sin(angle), 0.0, 0.0), self.START, []
        for target in times:
            steps = round((target - u) / self.STEP)
            h = (target - u) / steps
            for _ in range(steps):
                k1 = self.rate(u, w)
                k2 = self.rate(u + h / 2, tuple(a + h * k / 2 for a, k in zip(w, k1)))
                k3 = self.rate(u + h / 2, tuple(a + h * k / 2 for a, k in zip(w, k2)))
                k4 = self.rate(u + h, tuple(a + h * k for a, k in zip(w, k3)))
                w = tuple(a + h * (p + 2 * q + 2 * r + s) / 6 for a, p, q, r, s in zip(w, k1, k2, k3, k4))
                u += h
            out.append(w)
        return out

    def setUp(self):
        self.view, = embedding_files()["string_wave"]["views"]
        self.frames = {frame["value"]: frame for frame in self.view["movie"]["frames"]}

    def test_every_marked_particle_is_where_its_geodesic_puts_it_on_the_cone(self):
        times = [-1.5, 0.0, 1.0, 2.5]
        self.assertEqual([surface["time"] for surface in self.view["surfaces"]], times)
        for k in range(12):
            for u, (x, y, _, _) in zip(times, self.run_to(math.radians(30 * k), times)):
                rho = math.hypot(x, y)
                r = 2 * math.sqrt(rho)
                want = (r / 2 * x / rho, r / 2 * y / rho, math.sqrt(3) * r / 2)
                for got, expected in zip(self.frames[u]["dots"][k]["at"], want):
                    self.assertAlmostEqual(got, expected, delta=2e-5, msg=f"particle {k} at u = {u}")

    def test_the_ring_stays_on_the_cone_in_every_frame(self):
        for u, frame in self.frames.items():
            curve, = frame["curves"]
            for X, Y, Z in curve["points"]:
                self.assertAlmostEqual(Z, math.sqrt(3) * math.hypot(X, Y), delta=2e-6, msg=f"u = {u}")

    def test_the_caption_states_what_the_ring_does(self):
        def radii(u):
            return [2 * math.hypot(X, Y) for X, Y, _ in self.frames[u]["curves"][0]["points"]]
        caption = " ".join(self.view["caption"])
        self.assertIn("reaching from $r = 1.5\\,\\ell$ to $2.5\\,\\ell$ on the crest", caption)
        self.assertAlmostEqual(min(radii(0.0)), 1.5, delta=1e-3)
        self.assertAlmostEqual(max(radii(0.0)), 2.5, delta=1e-3)
        self.assertIn("left falling toward the string", caption)
        self.assertAlmostEqual(min(radii(-1.5)), 2.0, delta=1e-3)
        self.assertLess(max(radii(1.0)), 1.75)
        self.assertLess(max(radii(2.5)), max(radii(1.0)) - 0.3)


class ConformalDiagrams(unittest.TestCase):
    """A conformal diagram is of a whole spacetime, drawn from what its metrics publish, and
    stops being published when any of that changes."""

    SLICES = {"kerr", "kerr_newman", "cosmic_string"}

    def setUp(self):
        self.metrics = build.load_metrics()
        self.conformal = conformal_files()
        self.assertTrue(self.conformal, "no conformal diagram file was found, so nothing was checked")

    def test_every_file_was_drawn_from_what_its_metrics_publish(self):
        by_id = {m["id"]: m for m in self.metrics}
        for metric_id, data in self.conformal.items():
            self.assertTrue(data["source"], metric_id)
            for source in data["source"]:
                system = next(s for s in by_id[source["metric"]]["coordinates"] if s["id"] == source["system"])
                self.assertEqual(build.diagram_source_version(system, source["fields"]), source["version"],
                                 f"{metric_id} from {source['metric']}/{source['system']}")
        self.assertEqual(set(build.load_conformal(self.metrics)), set(self.conformal))

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
            self.load_folder({"interior_schwarzschild.json": self.conformal["interior_schwarzschild"]}, changed)
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
        self.assertTrue(seen, "no tower with r = 0 outside a region r < r_- was found, so nothing was checked")

    def test_a_change_that_draws_nothing_leaves_the_diagrams_standing(self):
        changed = copy.deepcopy(self.metrics)
        for metric in changed:
            metric["history"] = metric.get("history", "") + " More."
            for system in metric.get("coordinates") or []:
                system["domains"] = (system.get("domains") or []) + ["x \\in \\mathbb{R}"]
                for parameter in system.get("parameters", []):
                    parameter["description"] = parameter.get("description", "") + " Reworded."
        self.assertEqual(set(build.load_conformal(changed)), set(self.conformal))

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
        kerr = self.conformal["kerr"]
        for data in (dict(kerr, views=[]), {key: value for key, value in kerr.items() if key != "views"}):
            with self.assertRaises(build.DataError) as raised:
                self.load_folder({"kerr.json": data}, self.metrics)
            self.assertIn("draws nothing", str(raised.exception))

    def test_every_caption_names_what_is_drawn(self):
        for metric_id, data in self.conformal.items():
            for view in data["views"]:
                self.assertTrue(view["caption"], f"{metric_id} {view['id']}")
                self.assertTrue(opens_with_a_noun_phrase(view["caption"][0]), f"{metric_id} {view['id']}")

    def test_the_spacetimes_drawn_on_a_slice_say_so_on_the_diagram(self):
        """A view of a surface that is not the whole spacetime carries its restriction, and the
        three spacetimes with no faithful picture of the whole carry it on every view."""
        self.assertLessEqual(self.SLICES, set(self.conformal))
        for metric_id in self.SLICES:
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
        self.assertTrue(self.embedding, "no embedding diagram file was found, so nothing was checked")

    def test_btz_is_one_surface_through_both_exteriors(self):
        """The BTZ moment t = 0 as the file holds it: each exterior a piece in flat space from the
        throat r_+ = l to sqrt(2) l and a piece in Minkowski space beyond, rho = r throughout, the
        two parts of each exterior meeting level in one circle marked `space`, and the two
        exteriors meeting at the throat."""
        views = self.embedding["btz"]["views"]
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

    def test_the_spinning_string_leaves_the_axis_along_the_light_cone_and_ends_on_its_cone(self):
        """The spinning string's moment t = 0 as the file holds it, at b = 0.9 and a = 0.9: a piece
        in Minkowski space from the null circle, where it runs along the light cone, to the level
        circle R = a b/sqrt(1 - b^2), and a piece in flat space beyond, whose last chord climbs at
        sqrt(g_RR - 1), a little below the slope sqrt(1/b^2 - 1) of the cosmic string's cone; every chord of both is the
        proper distance sqrt(g_RR) dR of the published metric, g_RR = R^2/(b^2 (R^2 + a^2))."""
        a = b = 0.9
        level = a * b / math.sqrt(1 - b * b)
        surface, = self.embedding["spinning_string"]["views"][0]["surfaces"]
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
        self.assertEqual(set(build.load_embedding(self.metrics)), set(self.embedding))

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
            self.load_folder({"interior_schwarzschild.json": self.embedding["interior_schwarzschild"]}, changed)
        self.assertIn("schwarzschild.json", str(raised.exception))

    def test_a_change_that_draws_nothing_leaves_the_diagrams_standing(self):
        changed = copy.deepcopy(self.metrics)
        for metric in changed:
            metric["history"] = metric.get("history", "") + " More."
            for system in metric.get("coordinates") or []:
                system["domains"] = (system.get("domains") or []) + ["x \\in \\mathbb{R}"]
                for parameter in system.get("parameters", []):
                    parameter["description"] = parameter.get("description", "") + " Reworded."
        self.assertEqual(set(build.load_embedding(changed)), set(self.embedding))

    def test_a_file_without_a_metric_or_a_source_or_a_surface_is_refused(self):
        metric = {"id": "x", "name": "X", "short_name": "X", "tags": ["t"], "coordinates": [{"id": "a"}]}
        with self.assertRaises(build.DataError) as raised:
            self.load_folder({"ghost.json": {"metric": "ghost", "source": [], "views": []}}, [metric])
        self.assertIn("ghost", str(raised.exception))
        with self.assertRaises(build.DataError) as raised:
            self.load_folder({"x.json": {"metric": "x", "source": [], "views": []}}, [metric])
        self.assertIn("drawn from", str(raised.exception))
        flamm = self.embedding["schwarzschild"]
        for data in (dict(flamm, views=[]), dict(flamm, views=[dict(flamm["views"][0], surfaces=[])])):
            with self.assertRaises(build.DataError) as raised:
                self.load_folder({"schwarzschild.json": data}, self.metrics)
            self.assertIn("draws nothing", str(raised.exception))

    def test_every_spacetime_draws_and_a_file_that_draws_nothing_says_why(self):
        # Every spacetime but Lentz's has a surface drawn, and none says beside its views that it
        # draws nothing.
        self.assertEqual(set(self.embedding), {m["id"] for m in self.metrics} - {"lentz"})
        for name, data in self.embedding.items():
            self.assertTrue(data["views"], name)
            self.assertNotIn("stops", data, f"{name} keeps what it does not draw in its views")
        # A file with no surface to draw says why in sentences, and a file with views may not.
        flamm = self.embedding["schwarzschild"]
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
        self.assertTrue(marked, "no curve or point was found, so nothing was checked")

    def test_a_piece_that_doubles_back_or_crosses_the_axis_is_refused(self):
        for spoil, words in ((lambda points: points.reverse() or points.insert(1, points[0]), "one way"),
                             (lambda points: points[3].__setitem__(1, -0.5), "below the axis")):
            data = copy.deepcopy(self.embedding["schwarzschild"])
            spoil(data["views"][0]["surfaces"][0]["pieces"][0]["points"])
            with self.assertRaises(build.DataError) as raised:
                self.load_folder({"schwarzschild.json": data}, self.metrics)
            self.assertIn(words, str(raised.exception))

    def test_no_file_writes_a_negative_zero(self):
        # Rounding a small negative number leaves -0.0, which a client would print as "-0".
        # The match alone is reported, since the file is too long to print.
        for path in sorted(build.EMBEDDING_DIR.glob("*.json")):
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

    def test_teos_equator_is_flamms_paraboloid_and_his_throat_a_dumbbell(self):
        """Teo's eq. (27), z = +-2 sqrt(b_0 (r - b_0)) on the equator for every spin, and the
        throat's circles of radius b_0 (1 + cos^2 theta) sin theta at a = 1/4, the same on either
        side of its equator, from the numbers written and nothing else."""
        views = {view["id"]: view["surfaces"][0] for view in self.embedding["teo_wormhole"]["views"]}
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

    def test_the_surfaces_are_the_ones_known_in_closed_form(self):
        """Flamm's paraboloid, the interior Schwarzschild cap, the catenoid, the cone, Gott's cap
        and the sphere, from the numbers written and nothing else, to their rounding."""
        def piece(metric_id, piece_id, number=0, view=0):
            surface = self.embedding[metric_id]["views"][view]["surfaces"][number]
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
                    middle = 0.5 * (r0 + r1)
                    if abs(middle - horizon) > 0.02 and r1 - r0 > 1e-9:
                        slope = math.sqrt(1 / bardeen_f(middle) - 1)
                        self.assertLess(abs(abs(z1 - z0) / (r1 - r0) - slope), 2e-3 * (1 + slope), f"Bardeen {pid} at {middle}")
        cap = piece("bardeen", "inside", view=1)
        a = (1 / 3) ** 1.5
        for r, _, z in cap[1:]:
            if r < 0.02:
                self.assertLess(abs(2 * a * (z - cap[0][2]) / r ** 2 - 1), 2e-2, f"Bardeen's cap at {r}")
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
        for number, surface in enumerate(self.embedding["vaidya"]["views"][0]["surfaces"]):
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
        # The photon rocket's slices of constant cu + r at theta = pi/2: before the light of the burn
        # reaches them, z^2 = 8 m_0 r, and once it has passed the rim, z^2 = 8 m r with m = m_0 e^(-6/5),
        # the mass the burn leaves; the rim r = 12 m_0 stands at z = 0.
        rocket = self.embedding["photon_rocket"]["views"][0]["surfaces"]
        for number, mass in ((0, 1.0), (len(rocket) - 1, math.exp(-1.2))):
            for r, rho, z in piece("photon_rocket", "whole", number):
                near(rho, r, f"photon rocket rho at {r}")
                near(z, 2 * math.sqrt(2 * mass) * (math.sqrt(r) - math.sqrt(12)), f"photon rocket z at {r}")
        # Oppenheimer-Snyder's dust is a cap of a sphere of radius a out to chi0 = pi/4 at every
        # moment, and at the release, the first moment, the outside is Flamm's paraboloid.
        for number in range(len(self.embedding["oppenheimer_snyder"]["views"][0]["surfaces"])):
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
        vertical = self.embedding["tolman_bondi"]["views"][0]["vertical"]
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
        for number in range(len(self.embedding["frw"]["views"][0]["surfaces"])):
            hemisphere = piece("frw", "near", number)
            a = -hemisphere[0][2]
            for sign, pid in ((-1, "near"), (1, "far")):
                for r, rho, z in piece("frw", pid, number):
                    near(rho, a * r, f"sphere rho at {r}")
                    near(z, sign * a * math.sqrt(max(1 - r * r, 0)), f"sphere z at {r}")

    def test_the_eleven_drawn_last_are_the_surfaces_known_in_closed_form(self):
        """The plane, the hyperboloid, the tube of the Malament-Hogarth toy, Taub's round
        moment, and the flat planes with what is marked on them, from the numbers written."""
        def view(metric_id):
            return self.embedding[metric_id]["views"][0]

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
            return next(v for v in self.embedding[metric_id]["views"] if v["id"] == "ring")
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
        view = self.embedding[metric_id]["views"][0]
        piece = view["surfaces"][0]["pieces"][0]
        self.assertIn("grid", piece, metric_id)
        return view, piece

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
        flat = {"minkowski", "kasner", "bianchi", "pp_wave", "aichelburg_sexl", "khan_penrose", "bell_szekeres", "light_beam"}
        # The domain wall's moment ct = 0, when the wall stops, is the flat disc of radius 1/k taken
        # twice and joined at its rim; the moments either side of it are the cones it opens into.
        flat_moments = {("domain_wall", "moments", 2)}
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
        self.assertEqual(set(scaled), {"tolman_bondi", "szekeres"})
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
                         {"alcubierre", "krasnikov", "natario", "kasner", "bianchi", "pp_wave", "aichelburg_sexl", "khan_penrose",
                          "bell_szekeres", "light_beam"})

    def test_a_grid_that_is_not_one_is_refused(self):
        def spoil(change, words):
            data = copy.deepcopy(self.embedding["krasnikov"])
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
    turns, _tools/turn_check.cjs, shared by every test that reads it."""
    global _turn_check
    if shutil.which("node") is None:
        test.skipTest("Node is not installed, so the page's geometry cannot be run")
    if _turn_check is None:
        run = subprocess.run(["node", str(build.ROOT / "_tools" / "turn_check.cjs")],
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
        self.assertTrue(self.views)

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
        # Only the cone laid flat lies in the plane of the page.
        flat = {name for name, view in self.views if any(layer.get("flat") for layer in view["figure"]["layers"])}
        self.assertEqual(flat, {"cosmic_string"})

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
        self.assertEqual(set(heights), {"alcubierre", "krasnikov", "natario"})
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

    STACKS = {"kasner": 1.5, "bianchi": 2.5, "pp_wave": 0.5, "aichelburg_sexl": 0.5, "khan_penrose": 3.0,
              "bell_szekeres": 3.0, "light_beam": 0.4}   # the height of a unit of time
    # Every movie, by its spacetime and view, with its variable. The last nine stood as separate
    # pictures of their moments until the captain asked on 1 October 2026 for every one of them
    # to play, and TimeSlicedViewsAreMovies keeps any other from standing so again.
    MOVIES = {("frw", "closed"): "$ct$", ("malament_hogarth", "plane"): "$ct$", ("mixmaster", "sphere"): "$c\\tau$",
              ("oppenheimer_snyder", "collapse"): "$c\\tau$", ("vaidya", "shell"): "$v - r$",
              ("cosmic_string", "unroll"): "$\\Delta\\phi$", ("milne", "hyperboloids"): "$ct$",
              ("einstein_rosen_waves", "pulse"): "$ct$", ("nariai", "universe"): "$ct$",
              ("domain_wall", "moments"): "$kct$", ("kantowski_sachs", "dust"): "$\\eta$",
              ("robinson_trautman", "fronts"): "$cu$", ("mcvittie", "flamm"): "$ct$",
              ("tolman_bondi", "cloud"): "$ct$", ("szekeres", "equators"): "$ct$", ("misner", "cylinders"): "$ct$",
              ("photon_rocket", "burn"): "$cu + r$",
              ("gott_time_machine", "cylinders"): "$c\\tau$", ("kantowski_sachs", "vacuum"): "$c\\tau$",
              ("ori_time_machine", "throat"): "$t$",
              ("kasner", "ring"): "$t$", ("bianchi", "ring"): "$c\\bar Ht$", ("pp_wave", "ring"): "$cu$",
              ("aichelburg_sexl", "ring"): "$u$", ("khan_penrose", "ring"): "$\\tau$",
              ("bell_szekeres", "ring"): "$\\xi$",
              ("gowdy", "torus"): "$t$", ("light_beam", "ring"): "$u$",
              ("string_wave", "ring"): "$u$", ("simpson_visser", "inside"): "$c\\tau$"}

    def setUp(self):
        self.embedding = embedding_files()

    def views(self, metric_id):
        return {v["id"]: v for v in self.embedding[metric_id]["views"]}

    def movies(self):
        """Every view with a movie, by its spacetime and view."""
        return {(name, v["id"]): v for name, data in self.embedding.items() for v in data["views"] if "movie" in v}

    def test_the_stacks_come_first_with_their_flat_rings_beside_them(self):
        for metric_id in self.STACKS:
            self.assertEqual([v["id"] for v in self.embedding[metric_id]["views"]], ["tube", "ring"], metric_id)
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
        kasner = self.views("kasner")["tube"]["surfaces"][0]["pieces"][0]["grid"]
        for u, a, b in zip(kasner["u"], kasner["a"], kasner["b"]):
            self.assertAlmostEqual(a, u ** (-2 / 7), delta=1e-6)
            self.assertAlmostEqual(b, u ** (6 / 7), delta=1e-6)
        wave = self.views("pp_wave")["tube"]["surfaces"][0]["pieces"][0]["grid"]
        self.assertEqual(wave["b"][-1], 0)

    def test_every_moment_label_of_a_stack_follows_its_ring_and_the_axis_is_named(self):
        for metric_id in self.STACKS:
            view = self.views(metric_id)["tube"]
            labels = view["figure"]["labels"]
            named = [L for L in labels if "curve" in L]
            self.assertEqual([L["curve"]["curve"] for L in named], list(range(len(named))), metric_id)
            self.assertEqual([L["text"] for L in named], [c["label"] for c in view["surfaces"][0]["curves"]], metric_id)
            self.assertEqual([L["curve"]["side"] for L in named], [1, -1] * (len(named) // 2), metric_id)
            axis, = [L for L in labels if "axis" in L]
            self.assertEqual(axis["axis"], {"surface": 0}, metric_id)
            self.assertIn(("line", "axis"), [(k, c) for k, c, _ in view["figure"]["legend"]], metric_id)

    def test_the_movies_are_these(self):
        self.assertEqual({key: v["movie"]["variable"] for key, v in self.movies().items()}, self.MOVIES)
        self.assertEqual({name for (name, _), v in self.movies().items() if v["movie"].get("turns") is False},
                         {"frw", "milne"})

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
            if metric_id == "cosmic_string":
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

    def test_each_movie_with_a_rim_holds_it_still_while_the_rest_moves(self):
        # The flat plane beyond the Malament-Hogarth ball, the clocks released at 4 r_s outside
        # the collapsing dust, Vaidya's r = 4 r_s and the clocks released at 2.5 r_b outside
        # Tolman-Bondi's cloud stand at z = 0 in every frame, and the
        # Malament-Hogarth well only deepens as the removed event nears.
        for metric_id, piece in (("malament_hogarth", "flat"), ("oppenheimer_snyder", "exterior"), ("vaidya", None),
                                 ("tolman_bondi", "exterior")):
            view = self.embedding[metric_id]["views"][0]
            for frame in view["movie"]["frames"]:
                outer = frame["pieces"][-1] if piece is None else next(p for p in frame["pieces"] if p["id"] == piece)
                self.assertEqual(outer["points"][-1][2], 0, f"{metric_id} at {frame['value']}")
        depths = [-f["pieces"][0]["points"][0][2] for f in self.embedding["malament_hogarth"]["views"][0]["movie"]["frames"]]
        self.assertTrue(all(b >= a for a, b in zip(depths, depths[1:])), "the Malament-Hogarth well rises")

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
        view = copy.deepcopy(next(v for v in self.embedding["cosmic_string"]["views"] if "movie" in v))
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
        movie = next(v for v in self.embedding["cosmic_string"]["views"] if "movie" in v)["movie"]
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
        self.assertTrue(self.movies)

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
        self.assertEqual(len(re.findall(r"\br[xy]\.textContent = ", page)), 2)
        self.assertRegex(page, r"\.mfs-nr-readout \.nr-ry \{ white-space: pre; font-variant-numeric: tabular-nums; \}")


class TimeSlicedViewsAreMovies(unittest.TestCase):
    """No embedding diagram sets the moments of a time out as separate pictures, as the captain
    asked on 1 October 2026, when Tolman-Bondi's cloud still stood as four of them beside the
    movies: a view that changes with a time, t, tau, eta, u or any other, plays as a movie with
    the page's play and pause button and a label naming the frame shown. The script that draws
    them, the build that reads them and the published files are each held to it."""

    def setUp(self):
        self.embedding = embedding_files()
        self.views = [(name, view) for name, data in self.embedding.items() for view in data["views"]]

    def test_every_view_of_more_than_one_moment_is_a_movie(self):
        sliced = {(name, view["id"]) for name, view in self.views
                  if len(view["surfaces"]) > 1 or any("time" in s for s in view["surfaces"])}
        self.assertTrue(sliced)
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
        view = copy.deepcopy(next(v for v in self.embedding["tolman_bondi"]["views"] if v["id"] == "cloud"))
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
            for one in self.embedding[metric_id]["views"]:
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
        self.assertTrue(self.figures)

    def test_every_figure_of_light_cones_turns(self):
        checked = {f"{v['metric']}/{v['view']}" for v in turn_check(self)["figures"]}
        self.assertEqual(checked, set(self.figures))
        self.assertEqual(checked, {"alcubierre/bubble", "godel/tipping", "gott_time_machine/loop", "kerr/dragging",
                                   "kerr_newman/dragging", "spinning_string/tipping", "stockum_dust/tipping", "wormhole_time_machine/trip"})

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
        # side, t left out, so neither has another side.
        still = {f"{name}/{figure['id']}" for name, data in diagram_files().items()
                 for figures in data.get("projections", {}).values() for figure in figures if "turn" not in figure}
        self.assertEqual(still, {"cosmic_string/beam", "light_beam/lens"})


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


# Simpson and Visser's three geometries at r_s = 1, and the embedding views each one's drawings mark.
SV_A = {"bounce": 0.5, "null": 1.0, "wormhole": 2.0}
SV_MOMENTS = ("outside", "inside", "null", "wormhole")
SV_MOMENTS_OF = {"bounce": ("outside", "inside"), "null": ("null",), "wormhole": ("wormhole",)}


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


def novikov_t(R, tau):
    """A shell of dust released from rest at areal radius R, r_s = 1, at its proper time tau:
    its r and Schwarzschild t, from the cycloid and Misner, Thorne and Wheeler's (31.10)."""
    eta = bisect(lambda e: 0.5 * R * math.sqrt(R) * (e + math.sin(e)) - tau, 0.0, math.pi)
    k, tan = math.sqrt(R - 1), math.tan(eta / 2)
    r = 0.5 * R * (1 + math.cos(eta))
    return r, math.log(abs((k + tan) / (k - tan))) + k * (eta + 0.5 * R * (eta + math.sin(eta)))


class Slices(unittest.TestCase):
    """Every moment an embedding diagram is cut from is drawn on its spacetime's other diagrams
    where it lies, from the numbers in the files alone, and on nothing else."""

    # The drawings on which no moment of the spacetime's embedding lies: other universes,
    # another cloud, the time reversed shell, and cylinders where no surface of constant t is
    # a moment of space.
    HIDDEN = {"btz/stationary/rotating", "btz/eddington_finkelstein_ingoing/rotating",
              "btz/eddington_finkelstein_outgoing/rotating", "conformal btz/rotating",
              "frw/comoving_spherical/radial", "frw/comoving_spherical/through", "frw/conformal_spherical/radial",
              "tolman_bondi/comoving_synchronous/collapse", "vaidya/eddington_finkelstein_outgoing/shell",
              "godel/cylindrical/beyond", "stockum_dust/cylindrical/beyond", "conformal frw/flat", "conformal frw/open",
              "misner/rindler/plane", "conformal misner/rindler",
              # Gott's region of closed timelike curves, which no moment of Grant's Milne time meets,
              # and the centre of momentum chart about the strings.
              "gott_time_machine/grant_rindler/plane", "conformal gott_time_machine/grant_rindler",
              "gott_time_machine/centre_of_momentum/loop",
              # Bell and Szekeres's regular chart, whose planes of T and Z lie off eta = 0, where the ring is.
              "bell_szekeres/regular/plane", "conformal bell_szekeres/regular",
              # The spinning string's cylinders inside r_c, whose circles are closed timelike curves.
              *[f"spinning_string/{s}/inside" for s in ("proper_radius", "rescaled_radius", "helical")],
              # Gowdy's sphere chart draws the inside of Schwarzschild's horizon; the moments
              # embedded are the torus universe's.
              "gowdy/sphere/plane", "conformal gowdy/sphere",
              # The wormhole time machine's axis of acceleration, which the embedded plane theta = pi/2 does
              # not meet, the mouth of its short throat, and the flat space outside the mouths.
              "wormhole_time_machine/wormhole/speeding", "wormhole_time_machine/wormhole/slowing",
              "wormhole_time_machine/short_throat/radial", "wormhole_time_machine/lorentz/trip",
              "conformal wormhole_time_machine/lorentz",
              # The axis of the Curzon-Chazy particle, which the embedded plane z = 0 meets only at R = 0.
              "curzon_chazy/weyl/axis", "curzon_chazy/spherical/axis",
              "conformal curzon_chazy/weyl_axis", "conformal curzon_chazy/spherical_axis",
              # The axis of Zipoy and Voorhees's metric, which the embedded equatorial plane does not meet.
              *[f"zipoy_voorhees/{s}/axis_{k}" for s in ("spherical", "prolate_spheroidal") for k in ("oblate", "prolate")],
              *[f"conformal zipoy_voorhees/{s}_axis_{k}" for s in ("spherical", "prolate_spheroidal")
                for k in ("oblate", "prolate")],
              # The axis of Szekeres's cloud, which the embedded surface through the shells' equators meets
              # only at the centre.
              "szekeres/axisymmetric/north", "szekeres/axisymmetric/south",
              # The axis of the photon rocket's flight, which the embedded surface of the rays that leave
              # the rocket sideways meets only at r = 0.
              *[f"photon_rocket/{s}/{half}" for s in ("rectilinear", "robinson_trautman") for half in ("behind", "ahead")],
              "conformal photon_rocket/behind", "conformal photon_rocket/ahead",
              # Two beams of light side by side, another spacetime than the single beam whose wave fronts
              # are embedded, and the plane y = 0 with t left out, which every wave front covers whole.
              "light_beam/null_cartesian/midway", "light_beam/null_cartesian/one", "conformal light_beam/midway",
              "light_beam/cartesian/lens"}

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

    # The moments of a spacetime's embedding views, a stack's each ring it marks at a time.
    def every(self, metric_id):
        out = []
        for v in self.embedding[metric_id]["views"]:
            for i, surface in enumerate(v["surfaces"]):
                timed = [k for k, c in enumerate(surface.get("curves", [])) if "time" in c]
                out += [(v["id"], i, k) for k in timed] or [(v["id"], i, None)]
        return out

    # The embedding views a drawing does not mark though it marks others: Gott's core replaces
    # the apex of the ideal string's cone, which is another spacetime than Gott's, and one of
    # Majumdar and Papapetrou's holes alone is another spacetime than two of them.
    # Letelier's black hole, r_s > 0, and the monopole with no mass at its centre are two spacetimes of
    # one line element: the static and Eddington-Finkelstein drawings are the black hole's, and the
    # Barriola-Vilenkin drawings the monopole's.
    HIDDEN_VIEWS = {"conformal cosmic_string/gott": {"unroll"},
                    **{f"global_monopole/{s}": {"monopole"} for s in (
                        "static/radial", "eddington_finkelstein_ingoing/finkelstein", "eddington_finkelstein_ingoing/chart",
                        "eddington_finkelstein_outgoing/finkelstein", "eddington_finkelstein_outgoing/chart")},
                    **{f"conformal global_monopole/{v}": {"monopole"} for v in ("static", "ingoing", "outgoing")},
                    # Tangherlini's black hole in five dimensions and in six are two spacetimes, each
                    # marked on the drawings of its own charts.
                    **{f"tangherlini/{s}": {"six"} for s in (
                        "spherical/radial", "eddington_finkelstein_ingoing/finkelstein", "eddington_finkelstein_ingoing/chart",
                        "eddington_finkelstein_outgoing/finkelstein", "eddington_finkelstein_outgoing/chart")},
                    **{f"conformal tangherlini/{v}": {"six"} for v in ("spherical", "ingoing", "outgoing")},
                    "tangherlini/spherical_six/radial": {"five"},
                    "conformal tangherlini/six": {"five"},
                    "global_monopole/conical/radial": {"black_hole"},
                    "conformal global_monopole/conical": {"black_hole"},
                    # The threaded black hole's bifurcation sphere, its horizon view, lies at v -> -infinity and
                    # u -> +infinity, off both Eddington-Finkelstein charts.
                    **{f"string_black_hole/eddington_finkelstein_{way}/{v}": {"horizon"}
                       for way in ("ingoing", "outgoing") for v in ("finkelstein", "chart")},
                    # The dilaton black hole's Einstein metric and its two string metrics measure one
                    # manifold three ways: the Einstein charts mark the Einstein metric's moment and each
                    # string chart its own.
                    **{f"dilaton_black_hole/{s}": {"string_magnetic", "string_electric"} for s in (
                        "static/radial", "eddington_finkelstein_ingoing/finkelstein", "eddington_finkelstein_ingoing/chart",
                        "eddington_finkelstein_outgoing/finkelstein", "eddington_finkelstein_outgoing/chart")},
                    **{f"conformal dilaton_black_hole/{v}": {"string_magnetic", "string_electric"}
                       for v in ("static", "ingoing", "outgoing")},
                    "dilaton_black_hole/string_magnetic/radial": {"einstein", "string_electric"},
                    "conformal dilaton_black_hole/string_magnetic": {"einstein", "string_electric"},
                    "dilaton_black_hole/string_electric/radial": {"einstein", "string_magnetic"},
                    "conformal dilaton_black_hole/string_electric": {"einstein", "string_magnetic"},
                    # Zipoy and Voorhees's oblate and prolate masses are two spacetimes of one line element,
                    # each equatorial drawing marking its own moment.
                    **{f"zipoy_voorhees/{s}/equator_{k}": {o} for s in ("spherical", "prolate_spheroidal")
                       for k, o in (("oblate", "prolate"), ("prolate", "oblate"))},
                    **{f"conformal zipoy_voorhees/{s}_equator_{k}": {o} for s in ("spherical", "prolate_spheroidal")
                       for k, o in (("oblate", "prolate"), ("prolate", "oblate"))},
                    "majumdar_papapetrou/cartesian/tz": {"one_hole"},
                    "majumdar_papapetrou/cartesian/tx": {"one_hole"},
                    "majumdar_papapetrou/cylindrical/radial": {"one_hole"},
                    "majumdar_papapetrou/isotropic/radial": {"two_holes"},
                    "conformal majumdar_papapetrou/one_hole": {"two_holes"},
                    # Melvin's universe with no hole and Ernst's hole inside it are two spacetimes of one entry,
                    # each chart's drawings marking its own moment.
                    "melvin/cylindrical/radial": {"ernst"}, "conformal melvin/cylindrical": {"ernst"},
                    "melvin/ernst/radial": {"universe"}, "conformal melvin/ernst": {"universe"},
                    # The dust universe and the inside of Schwarzschild's horizon are two members of
                    # the Kantowski-Sachs family, each marked on its own drawings.
                    "kantowski_sachs/comoving/tr": {"vacuum"}, "kantowski_sachs/dust/etar": {"vacuum"},
                    "conformal kantowski_sachs/dust": {"vacuum"},
                    "kantowski_sachs/schwarzschild_interior/Tr": {"dust"},
                    "conformal kantowski_sachs/vacuum": {"dust"},
                    # Teo's equatorial plane, embedded at a = 1, does not meet the axis, and his throat,
                    # embedded at a = 1/4, is marked where its poles meet the axis, drawn at that spin.
                    **{f"teo_wormhole/{s}/{surface}": {other} for s in ("spherical", "proper_radial")
                       for surface, other in (("axis", "equator"), ("equator", "throat"))},
                    **{f"conformal teo_wormhole/{s}_{surface}": {other} for s in ("spherical", "proper_radial")
                       for surface, other in (("axis", "equator"), ("equator", "throat"))},
                    # Simpson and Visser's black bounce, one way wormhole and traversable wormhole are three
                    # spacetimes of one line element, a = r_s/2, r_s and 2 r_s, each drawing marking its own
                    # moments. The outgoing chart's region between the horizons is the white hole, where the
                    # black hole's moments of constant r do not lie.
                    **{f"simpson_visser/{s}/{case}": set(SV_MOMENTS) - set(own) | (
                        {"inside"} if s.endswith("outgoing") else set())
                       for s in ("spherical", "areal", "eddington_finkelstein_ingoing", "eddington_finkelstein_outgoing")
                       for case, own in SV_MOMENTS_OF.items()},
                    **{f"conformal simpson_visser/{s}_{case}": set(SV_MOMENTS) - set(own)
                       for s in ("spherical", "areal", "ingoing", "outgoing") for case, own in SV_MOMENTS_OF.items()}}
    # The surfaces of a view that a drawing does not mark though it marks the view's others: the areal
    # chart of the black bounce is drawn on the side r > 0, where the moments of r < 0 do not lie.
    HIDDEN_SURFACES = {"simpson_visser/areal/bounce": {("inside", 3), ("inside", 4)}}

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
                     and m[:2] not in self.HIDDEN_SURFACES.get(where, set())]
            self.assertEqual(marks, every, where)
            drawn += len(marks)
        self.assertGreater(drawn, 100)

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
        if key.startswith("robinson_trautman/"):
            # The fronts of one retarded time differ only in size, so a moment is the whole
            # outgoing ray u = u_k, drawn against r and cu + r, from r = 0 to the box.
            return (lambda X: t + X), [0.0]
        if key.startswith("string_black_hole/") and not mark["lines"]:
            # The horizon's bifurcation sphere, the point t = 0, r = r_s of the static planes.
            return (lambda X: 0.0), [1.0]
        if key.startswith("string_black_hole/") and "eddington" not in key:
            return (lambda X: 0.0), list(self.reach(surface, "static"))
        if key.startswith(("schwarzschild/eddington_finkelstein", "string_black_hole/eddington_finkelstein")):
            sign = 1 if "ingoing" in key else -1
            finkelstein = key.endswith("finkelstein")
            return (lambda X: sign * math.log(X - 1) + (0 if finkelstein else sign * X)), None
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
        if key.startswith("tangherlini/eddington_finkelstein"):
            # Five dimensions at r_h = 1: r_* = r + ln((r - 1)/(r + 1))/2, and the static t = 0 is v = r_*
            # and u = -r_*, drawn against v - r and u + r or against v and u.
            sign = 1 if "ingoing" in key else -1
            finkelstein = key.endswith("finkelstein")
            return (lambda X: sign * (0.5 * math.log(abs((X - 1) / (X + 1))) + (0 if finkelstein else X))), \
                list(self.reach(surface))
        if key.startswith("bardeen/eddington_finkelstein"):
            # At r_s = 1 and g = 1/3 the static t = 0 is v = r_* and u = -r_*, with dr_*/dr = 1/f and
            # r_* = 0 at the centre, drawn against v - r and u + r or against v and u. It runs off
            # toward the horizon its view ends on, so a line ends at the box or at the embedding's reach.
            sign = 1 if "ingoing" in key else -1
            finkelstein = key.endswith("finkelstein")
            return (lambda X: sign * (bardeen_rstar(X) - (X if finkelstein else 0))), list(self.reach(surface))
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
        if key == "oppenheimer_snyder/interior_comoving/through":
            return (lambda X: t / (2 * math.sqrt(2))), [0, math.pi / 4]
        if key.startswith("einstein_rosen_waves/"):
            # A moment ct = T, drawn against rho and ct in both charts, out to the embedding's reach.
            return (lambda X: t), list(self.reach(surface))
        if key == "gowdy/areal/plane":
            # A moment of the areal time, drawn against theta and t, once round the torus.
            return (lambda X: t), list(self.reach(surface))
        if key == "gowdy/logarithmic/plane":
            # The same moment against theta and -tau = ln t.
            return (lambda X: math.log(t)), list(self.reach(surface))
        if key == "vaidya/eddington_finkelstein_ingoing/shell":
            return (lambda X: t), list(self.reach(surface))
        if key == "krasnikov/cylindrical/tx":
            path = next(c for c in surface["curves"] if c["class"] == "path")
            grid = next(p for p in surface["pieces"] if "grid" in p)["grid"]
            return (lambda X: 5.0), [grid["u"][0] - path["points"][0][0], grid["u"][-1] - path["points"][0][0]]
        if key.startswith(("kasner", "bianchi", "malament_hogarth")):
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
        if key in ("misner/milne/plane", "gott_time_machine/grant_milne/plane"):
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
        if key.startswith(("godel/cylindrical", "stockum_dust/cylindrical", "minkowski/rindler")):
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
        if key.startswith("bertotti_robinson"):
            lo, hi = self.reach(surface) if mark["lines"] else (1, 1)
            return (lambda X: 0.0), [lo, hi]
        if key == "einstein_static/static_areal/radial":
            # r = R sin chi carries the near hemisphere of the moment, out to the equator r = R.
            return (lambda X: 0.0), [0.0, 1.0]
        if key.startswith("milne/"):
            # The moment ct: level in the comoving charts, at c tau = ln(ct) in the logarithmic one
            # and the hyperbola cT = sqrt(c^2t^2 + R^2) in the inertial one, out to the chi reached.
            hi = self.reach(surface)[1]
            if key == "milne/inertial/through":
                return (lambda X: math.sqrt(t * t + X * X)), [0.0, t * math.sinh(hi)]
            if key == "milne/logarithmic_time/radial":
                return (lambda X: math.log(t)), [0.0, hi]
            return (lambda X: t), [0.0, math.sinh(hi) if "spherical" in key else hi]
        if key.startswith("mcvittie/"):
            # A moment of cosmic time is level in both charts, from the throat out to the comoving r the
            # embedding reaches, which the areal chart draws at R = ar(1 + 1/(4ar))^2, r_s = 1, with
            # a = sinh^(2/3)(3 H_0 t/2) and H_0 = 1/sqrt(15).
            ends = list(self.reach(surface))
            if "/areal/" in key:
                a = math.sinh(1.5 * t / math.sqrt(15)) ** (2 / 3)
                ends = [a * r * (1 + 1 / (4 * a * r)) ** 2 for r in ends]
            return (lambda X: t), ends
        if key.startswith("zipoy_voorhees/prolate_spheroidal/"):
            # The prolate spheroidal x is r/m - 1 of the circles the embedding reaches in r, at m = 1.
            return (lambda X: 0.0), [r - 1 for r in self.reach(surface)]
        if key == "levi_civita/kasner/radial":
            # The Kasner form's r is the proper distance from the axis, rho^Sigma/Sigma with
            # Sigma = 3/4 at sigma = 1/4, of the circles the embedding reaches in Weyl's rho.
            return (lambda X: 0.0), [rho ** 0.75 / 0.75 for rho in self.reach(surface)]
        if key == "kaluza_klein_monopole/taub_nut/radial":
            # The Taub-NUT radius is rho = r + 2m, at m = 1, of the circles the cigar reaches in r.
            return (lambda X: 0.0), [r + 2 for r in self.reach(surface)]
        if key == "damour_solodukhin/isotropic/radial":
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
        if key == "anti_de_sitter/poincare/tx":
            hi = self.reach(surface)[1]
            x = math.sqrt(2 * (math.sqrt(1 + hi * hi) - 1))
            return (lambda X: 0.0), [-x, x]
        if key == "godel/cartesian/tx":
            hi = self.reach(surface)[1]
            return (lambda X: 0.0), [-2 * hi, 2 * hi]
        if surface["pieces"][0].get("grid"):
            u = surface["pieces"][0]["grid"]["u"][-1]
            return (lambda X: 0.0), [-u, u]
        lo, hi = self.reach(surface)
        through = key.endswith("/tx") and not view["mirror"]
        return (lambda X: 0.0), ([-hi, -lo, lo, hi] if through else [lo, hi])

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
                        Y_of, ends = self.flat_expected(key, view, surface, mark)
                        for line in mark["lines"] + [[p] for p in mark["points"]]:
                            for u in line:
                                X, Y = X0 + u[0] * (X1 - X0), Y0 + u[1] * (Y1 - Y0)
                                h = 2e-4 * (X1 - X0)
                                slope = (abs(Y_of(min(X + h, X1)) - Y_of(max(X - h, X0 + 1e-9 if key.startswith(
                                    ("schwarzschild/edd", "string_black_hole/edd", "oppenheimer_snyder/ext")) else X0))) / (2 * h)
                                         if X0 < X < X1 else 0)
                                tol = 1e-4 * (Y1 - Y0) + slope * 1e-4 * (X1 - X0) + 1e-9
                                self.assertLess(abs(Y - Y_of(X)), 3 * tol, f"{key} {mark['label']} at {u}")
                                checked += 1
                            for u in (line[0], line[-1]):
                                edge = min(u[0], 1 - u[0], u[1], 1 - u[1]) < 1.5e-4
                                X = X0 + u[0] * (X1 - X0)
                                at_reach = ends is not None and any(abs(X - e) < 1.5e-4 * (X1 - X0) for e in ends)
                                self.assertTrue(edge or at_reach, f"{key} {mark['label']} stops at {u}")
        self.assertGreater(checked, 200)

    def sv_inside_r(self, surface):
        """The r of a moment of the black bounce between its horizons, a = r_s/2: its cylinder's radius is
        sqrt(r^2 + a^2), and r is positive before the proper time of r = 0, the third moment's."""
        radius = surface["pieces"][0]["points"][0][1]
        inside = next(v for v in self.embedding["simpson_visser"]["views"] if v["id"] == "inside")
        turn = inside["surfaces"][2]["time"]
        return math.copysign(math.sqrt(max(radius * radius - 0.25, 0.0)), turn - surface["time"])

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
        """Kerr's equator from above: the plane outside the horizon, the box less the disc of r_+."""
        X0, X1, Y0, Y1 = view["box"]
        lo = self.reach(surface)[0]
        outer, hole = mark["fills"][0]
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
                    if metric_id == "cosmic_string":
                        hi = self.reach(surface, "conical", reference=True)[1]
                        self.assertLessEqual(max(max(r) for r in rings), hi + 1e-3, where)
                        self.assertTrue(all(abs(r - hi) < 1e-3 for rim in rims for r in rim), where)
                        continue
                    for ring in rings + rims:
                        self.assertLess(max(ring) - min(ring), 2e-3 * max(ring), where)
                    if metric_id == "godel":
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
                    if metric_id == "frw":
                        eta = bisect(lambda e: e - math.sin(e) - t, 0, 2 * math.pi)
                        self.assertTrue(all(abs(T - eta) < 2e-4 for _, T in points), where)
                        self.assertEqual(sorted(X for X, _ in points), [0, round(math.pi, 4)], where)
                    elif metric_id in ("malament_hogarth", "einstein_rosen_waves", "gowdy"):
                        lo, hi = self.reach(surface)
                        for X, T in points:
                            p, q = (T - X) / 2, (T + X) / 2
                            tp, tq = math.tan(p), math.tan(q)
                            scale = 1 + tp * tp + tq * tq
                            self.assertLess(abs((tp + tq) / 2 - t), 2e-4 * scale, where)
                            self.assertLessEqual(lo - 2e-4 * scale, (tq - tp) / 2, where)
                            self.assertLessEqual((tq - tp) / 2, hi + 2e-4 * scale, where)
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
                    elif metric_id == "ori_time_machine":
                        # The hyperbola of the covering plane on which T = t, every copy: its null
                        # coordinates -2 e^(-z/2) and 2T e^(z/2) multiply to -4t.
                        for X, T in points:
                            tp, tq = math.tan((T - X) / 2), math.tan((T + X) / 2)
                            self.assertLess(abs(tp * tq + 4 * t), 2e-4 * (1 + tp * tp + tq * tq), where)
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
                    elif metric_id == "milne":
                        # p, q = arctan(ct e^-chi), arctan(ct e^chi), so tan p tan q = c^2t^2.
                        for X, T in points:
                            tp, tq = math.tan((T - X) / 2), math.tan((T + X) / 2)
                            self.assertLess(abs(tp * tq - t * t), 2e-3 * (1 + tp * tp) * (1 + tq * tq), f"{where} at {(X, T)}")
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
                    elif metric_id == "bell_szekeres":
                        # The event on eta = 0 at xi = t: on the plane x = y = 0, p = q = t/2 at a = b = 1,
                        # and on the strip of the anti-de Sitter factor, rho = 0 at chi = t - pi/2.
                        height = t if view["id"] in ("double_null", "time_space") else t - math.pi / 2
                        for X, T in points:
                            self.assertLess(abs(X) + abs(T - height), 2e-4, where)
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
                    elif metric_id == "rn_metric" and mark["view"] == "inside":
                        self.assertTrue(all(abs(T - math.pi) < 2e-4 for _, T in points), where)
                    elif metric_id == "simpson_visser" and mark["view"] == "inside":
                        # Kruskal's square of the black bounce, a = r_s/2: p, q = arctan e^(-+k(ct -+ r_*)) with
                        # k = sqrt(3)/4, so tan p tan q = e^(2 k r_*) on the moment of constant r.
                        want = math.exp(2 * math.sqrt(3) / 4 * sv_rstar(self.sv_inside_r(surface), 0.5))
                        for X, T in points:
                            tp, tq = math.tan((T - X) / 2), math.tan((T + X) / 2)
                            self.assertLess(abs(tp * tq - want), 2e-3 * (1 + tp * tp) * (1 + tq * tq), f"{where} at {(X, T)}")
                    elif metric_id == "bardeen":
                        # Inside r- the moment runs through the inner bifurcation sphere, at T = pi, and
                        # outside r+ through the outer one, at T = 0, or a period up, at T = 2 pi, on the
                        # outgoing chart's view, whose exterior is the one above the white hole.
                        height = math.pi if mark["view"] == "inside" else 2 * math.pi if view["id"] == "outgoing" else 0.0
                        self.assertTrue(all(abs(T - height) < 2e-4 for _, T in points), where)
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
        self.view = self.embedding["einstein_static"]["views"][0]
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
        views = [v["id"] for v in self.conformal["einstein_static"]["views"] if not v.get("system")]
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

    def test_the_shaded_region_matches_each_map_at_several_eta(self):
        for vid, shade in self.shades.items():
            times = [T for T, *_ in shade["reach"]]
            self.assertGreaterEqual(len(set(times)), 13, vid)
            self.assertIn(shade["T"], times, vid)
            cover = [L for L in next(v for v in self.conformal["einstein_static"]["views"] if v["id"] == vid)["layers"]
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

    def test_the_caption_says_the_moment_shaded(self):
        for vid, shade in self.shades.items():
            self.assertIn("pink, $\\eta = 0$", shade["caption"], vid)

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

    def test_the_collection_as_it_stands_holds_together(self):
        self.assertEqual(self.problems(build.load_metrics()), [])

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

    def test_the_page_lists_them_right_after_the_history(self):
        page = (build.ROOT / "_layouts" / "mfs.html").read_text(encoding="utf-8")
        body = page[page.index("var body = ["):]
        history, related, coordinates = (body.index(mark) for mark in (
            "mfs-section-label\">history<", "relatedSection(data),", "mfs-section-label\">coordinates<"))
        self.assertLess(history, related)
        self.assertLess(related, coordinates)
        self.assertIn("mfs-section-label\">Related Spacetimes<", page)


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
    def test_the_collection_as_it_stands_resolves(self):
        entries = build.build_references()["entries"]
        self.assertIsNone(build.check_citations(build.load_metrics(), entries))

    def test_every_metric_lists_its_references_in_the_order_its_prose_first_cites_them(self):
        """The page numbers a citation by its place in `references`, so [1] is the first one read:
        the history's citations in order, then those of the related spacetimes, which stand
        right after it."""
        for metric in build.load_metrics():
            read = [metric.get("history") or ""] + [entry["text"] for entry in metric.get("related") or []]
            self.assertEqual(metric.get("references", []), build.cited_keys(" ".join(read)), metric["id"])

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



class WormholeTrip(unittest.TestCase):
    """The round trip of the wormhole time machine's right mouth, as its figure in three
    dimensions publishes it in (Z, X, T) with c = 1: the physics the figure states, held to the
    numbers on disk."""

    def setUp(self):
        figure = diagram_files()["wormhole_time_machine"]["projections"]["lorentz"][0]
        self.lines = {}
        for line in figure["turn"]["lines"]:
            self.lines.setdefault(line["class"], []).append(line["points"])

    @staticmethod
    def interval(a, b):
        return (b[0] - a[0]) ** 2 + (b[1] - a[1]) ** 2 - (b[2] - a[2]) ** 2

    def test_both_mouths_move_slower_than_light_and_the_right_one_comes_home(self):
        left, right = self.lines["world"]
        self.assertEqual({(p[0], p[1]) for p in left}, {(0.0, 0.0)})
        for a, b in zip(right, right[1:]):
            self.assertGreater(b[2], a[2])
            self.assertLessEqual(self.interval(a, b), 1e-9)
        self.assertAlmostEqual(right[0][0], 10.0, places=5)
        self.assertAlmostEqual(right[-1][0], 10.0, places=5)
        self.assertAlmostEqual(max(p[0] for p in right), 31.48, places=2)

    def test_the_closed_geodesic_is_null_and_the_first_pair_of_ends_light_can_join(self):
        (start, end), = self.lines["ergo"]
        self.assertEqual(start[:2], [0.0, 0.0])
        self.assertAlmostEqual(self.interval(start, end), 0.0, delta=1e-6)
        self.assertAlmostEqual(start[2], 41.807, places=3)
        # Each dotted line joins the mouths at one proper time, the left mouth's being its T:
        # spacelike before the closed null geodesic and timelike after.
        pairs = self.lines["axis"]
        self.assertEqual([left[2] for left, _ in pairs], [9.0 * k for k in range(len(pairs))])
        for left, right in pairs:
            self.assertEqual(self.interval(left, right) < 0, left[2] > start[2], left)

    def test_the_closed_curve_is_timelike_and_inside_the_horizon(self):
        (start, end), = self.lines["ctc"]
        self.assertLess(self.interval(start, end), 0)
        self.assertGreater(end[2], start[2])
        (apex, _), = self.lines["ergo"]
        self.assertLess(self.interval(apex, start), 0)
        # Its far end is the right mouth at the proper time it left the left mouth at, after the
        # trip, when the right mouth's T runs ahead of its proper time by the time shift.
        self.assertAlmostEqual(end[2] - start[2], 21.718, places=3)
        self.assertAlmostEqual(end[0], 10.0, places=5)
        # Every generator drawn of the horizon is null, from the event the geodesic leaves.
        generators = [line for line in self.lines["horizon"] if len(line) == 2]
        self.assertEqual(len(generators), 12)
        for a, b in generators:
            self.assertEqual(a, apex)
            self.assertAlmostEqual(self.interval(a, b), 0.0, delta=1e-4)


if __name__ == "__main__":
    unittest.main()
