#!/usr/bin/env python3
"""Tests for the MFS data generator.

    python3 -m unittest discover -s _tools
"""

import contextlib
import copy
import decimal
import io
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
HYPHENATED_TERMS = {"anti-de", "anti-trapped", "plane-fronted", "scalar-tensor"}

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
    for position, system in enumerate(metric.get("coordinates") or []):
        where = f"coordinates[{system.get('id') or position}]"
        if system.get("name"):
            yield f"{where}.name", system["name"]
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
        for number, surface in enumerate(view["surfaces"]):
            at = f"{where}.surfaces[{number}]"
            if surface.get("label"):
                yield f"{at}.label", surface["label"]
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
    of a drawing. Labels and legends name things and are not sentences."""
    for metric in build.load_metrics():
        for field in ("history", "convention"):
            for number, paragraph in build.prose_paragraphs(metric.get(field) or ""):
                yield f"{metric['id']}.json: {field} paragraph {number}", paragraph
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
                                   ("#mfs-content-panel .mfs-turn-reset", True), (".mfs-result", False)):
            with self.subTest(button):
                self.assertEqual(self.last(button + ":active", "color"), "var(--pink-light)")
                if has_border:
                    self.assertEqual(self.last(button + ":active", "border-color"), "var(--pink-light)")
                places = [n for n, (selectors, _, printed) in enumerate(self.rules) if not printed
                          for s in selectors if s in (button + ":hover", button + ":active")]
                pressed = max(n for n, (selectors, _, _) in enumerate(self.rules) if button + ":active" in selectors)
                self.assertEqual(max(places), pressed)


class ConventionShape(unittest.TestCase):
    """Every convention is one paragraph or more, each of three to six sentences, broken
    only where a sentence ends."""

    HISTORY = "¶".join(HistoryShape.paragraph(3) for _ in range(5))

    @classmethod
    def metric(cls, convention, metric_id="x"):
        return dict(HistoryShape.metric(cls.HISTORY, metric_id), convention=convention)

    @staticmethod
    def convention(*counts):
        return "¶".join(HistoryShape.paragraph(count, "Signature") for count in counts)

    def test_every_convention_on_disk_keeps_its_shape(self):
        for metric in build.load_metrics():
            if metric.get("convention"):
                with self.subTest(metric["id"]):
                    self.assertEqual(
                        build.shape_problems(metric["id"], "convention", metric["convention"], 1), [])

    def test_a_convention_may_be_a_single_paragraph(self):
        for counts in ((3,), (6,), (3, 6), (4, 4, 5, 3)):
            self.assertIsNone(build.check_prose_shape([self.metric(self.convention(*counts))]))

    def test_a_spacetime_may_carry_no_convention(self):
        self.assertIsNone(build.check_prose_shape([HistoryShape.metric(self.HISTORY)]))

    def test_a_paragraph_too_short_or_too_long_is_refused_by_its_place(self):
        for counts, place in (((2,), 1), ((7,), 1), ((4, 4, 2), 3), ((5, 7, 4), 2)):
            with self.assertRaises(build.DataError) as raised:
                build.check_prose_shape([self.metric(self.convention(*counts))])
            self.assertIn(f"the convention, {list(counts)}, has", str(raised.exception))
            self.assertIn(f"paragraph {place}", str(raised.exception))

    def test_the_longest_paragraph_may_be_at_most_twice_the_shortest(self):
        with mock.patch.object(build, "PARAGRAPH_SENTENCES", (1, 10)):
            self.assertIsNone(build.check_prose_shape([self.metric(self.convention(2, 4))]))
            with self.assertRaises(build.DataError) as raised:
                build.check_prose_shape([self.metric(self.convention(2, 5))])
        self.assertIn("the convention, [2, 5], has a longest paragraph more than 2 times",
                      str(raised.exception))

    def test_every_convention_out_of_shape_is_named_at_once(self):
        with self.assertRaises(build.DataError) as raised:
            build.check_prose_shape([self.metric(self.convention(9), "one"),
                                     self.metric(self.convention(4, 4), "fine"),
                                     self.metric(self.convention(3, 1), "two")])
        self.assertIn("one.json: the convention", str(raised.exception))
        self.assertIn("two.json: the convention", str(raised.exception))
        self.assertNotIn("fine.json", str(raised.exception))

    def test_a_convention_out_of_shape_leaves_both_published_files_alone(self):
        before = {path: path.read_text(encoding="utf-8") for path in (build.INDEX_FILE, build.REFERENCES_FILE)}
        broken = build.load_metrics()
        broken[0] = dict(broken[0], convention=self.convention(8))
        for argv in (["--check"], []):
            with mock.patch.object(build, "load_metrics", return_value=broken), \
                    contextlib.redirect_stderr(io.StringIO()):
                self.assertEqual(build.main(argv), 2)
        for path, text in before.items():
            self.assertEqual(path.read_text(encoding="utf-8"), text)

    def test_every_break_stands_where_a_sentence_ends(self):
        # A break takes the place of the space between two sentences, so turning every break
        # back into a space gives the prose as it reads unbroken, and splitting at the breaks
        # gives each paragraph with nothing to trim. A paragraph that leads into a table may
        # end on a colon, since the table is not prose.
        for metric in build.load_metrics():
            for field in ("history", "convention"):
                text = metric.get(field)
                if not text:
                    continue
                with self.subTest(f"{metric['id']}.json {field}"):
                    self.assertNotRegex(text, r"\s¶|¶\s|^¶|¶$|¶¶")
                    paragraphs = text.split("¶")
                    for before, after in zip(paragraphs, paragraphs[1:]):
                        if "TABLE::" in (before[:7], after[:7]):
                            continue
                        self.assertEqual(len(build.sentences(f"{before} {after}")),
                                         len(build.sentences(before)) + len(build.sentences(after)),
                                         f"a break inside a sentence, before {after[:40]!r}")

    def test_the_breaks_change_no_word_of_main(self):
        """With every break turned back into a space, each convention reads exactly as main's.

        Cutting main's running conventions into paragraphs moves no word of them, and this
        holds the cut to that. Once main's own conventions are in paragraphs, the shape check
        holds them there and a later change may reword a convention on purpose, so this has
        nothing left to guard and stands aside.
        """
        def main_text(path):
            shown = subprocess.run(["git", "show", f"main:{path.relative_to(build.ROOT).as_posix()}"],
                                   cwd=build.ROOT, capture_output=True, text=True)
            return shown.stdout if shown.returncode == 0 else None

        if shutil.which("git") is None or subprocess.run(
                ["git", "rev-parse", "--verify", "--quiet", "main^{commit}"],
                cwd=build.ROOT, capture_output=True).returncode != 0:
            self.skipTest("no main to compare with")
        on_main = {}
        for metric in build.load_metrics():
            text = main_text(build.METRICS_DIR / f"{metric['id']}.json")
            if text is not None and json.loads(text).get("convention"):
                on_main[metric["id"]] = json.loads(text)["convention"]
        if not any(build.shape_problems(i, "convention", c, 1) for i, c in on_main.items()):
            self.skipTest("main's conventions are already in paragraphs")
        for metric in build.load_metrics():
            if metric["id"] in on_main:
                with self.subTest(metric["id"]):
                    self.assertEqual(metric.get("convention", "").replace("¶", " "),
                                     on_main[metric["id"]].replace("¶", " "))


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
        # Every spacetime has a surface drawn, and none says beside its views that it draws nothing.
        self.assertEqual(set(self.embedding), {m["id"] for m in self.metrics})
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
                            self.assertIn(grid["frame"], ("polar", "cartesian"), where)
                            self.assertEqual(len(grid["z"]), len(grid["u"]), where)
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
        # Reissner-Nordstrom's circles have their areal radius on both views, and inside r- each
        # side starts level at r_q^2/r_s, 0.2304 r_s, and ends at r- = 0.36 r_s.
        for view, pids in ((0, ("exterior", "other_exterior")), (1, ("inside", "other_inside"))):
            for pid in pids:
                points = piece("rn_metric", pid, view=view)
                for r, rho, z in points:
                    near(rho, r, f"Reissner-Nordstrom rho at {r}")
                if view:
                    self.assertEqual((points[0][0], points[-1][0]), (0.2304, 0.36), pid)
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
        # Vaidya's slices of constant v - r: a flat disc inside the shell, and outside it Flamm's
        # paraboloid moved in by r_s, z = 2 sqrt(r), from the shell at z = 0 or from the axis.
        for number, surface in enumerate(self.embedding["vaidya"]["views"][0]["surfaces"]):
            ids = [p["id"] for p in surface["pieces"]]
            if "inside" in ids:
                self.assertTrue(all(z == 0 for _, _, z in piece("vaidya", "inside", number)))
                outside = piece("vaidya", "outside", number)
                shell = outside[0][0]
                for r, rho, z in outside:
                    near(z, 2 * (math.sqrt(r) - math.sqrt(shell)), f"Vaidya z at {r}")
            else:
                for r, rho, z in piece("vaidya", "whole", number):
                    near(z, 2 * math.sqrt(r), f"Vaidya z at {r}")
        # Oppenheimer-Snyder's dust is a cap of a sphere of radius a out to chi0 = pi/4 at every
        # moment, and at the release, the first moment, the outside is Flamm's paraboloid.
        for number in range(len(self.embedding["oppenheimer_snyder"]["views"][0]["surfaces"])):
            dust = piece("oppenheimer_snyder", "dust", number)
            a = dust[-1][1] / math.sin(math.pi / 4)
            for chi, rho, z in dust:
                near(rho, a * math.sin(chi), f"Oppenheimer-Snyder cap rho at {chi}")
                near(z - dust[0][2], a * (1 - math.cos(chi)), f"Oppenheimer-Snyder cap z at {chi}")
        for r, rho, z in piece("oppenheimer_snyder", "exterior"):
            near(z, 2 * math.sqrt(r - 1) - 2, f"Oppenheimer-Snyder release z at {r}")
        # Tolman-Bondi's cloud at its release: every shell at its label, R = r, and outside it
        # Flamm's paraboloid of 2GM/c^2 = r_b/2 from the surface.
        cloud, outside = piece("tolman_bondi", "cloud"), piece("tolman_bondi", "exterior")
        for r, rho, z in cloud + outside:
            near(rho, r, f"Tolman-Bondi release rho at {r}")
        for r, rho, z in outside:
            near(z - outside[0][2], 2 * math.sqrt(0.5 * (r - 0.5)) - 1, f"Tolman-Bondi release z at {r}")
        # Bertotti-Robinson's equator is a cylinder of radius b, z = b ln r, and its sphere of radius b.
        for r, rho, z in piece("bertotti_robinson", "cylinder"):
            near(rho, 1.0, f"Bertotti-Robinson rho at {r}")
            near(z, math.log(r), f"Bertotti-Robinson z at {r}")
        for theta, rho, z in piece("bertotti_robinson", "sphere", view=1):
            near(rho, math.sin(theta), f"Bertotti-Robinson sphere rho at {theta}")
            near(z, 1 - math.cos(theta), f"Bertotti-Robinson sphere z at {theta}")
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
        for r, rho, z in piece("cosmic_string", "exterior"):
            near(rho, 0.9 * r, f"cone rho at {r}")
            near(z, math.sqrt(0.19) * r, f"cone z at {r}")
        core = piece("cosmic_string", "core")
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
        # The flat planes: every point level, and rho the distance along the profile.
        for metric_id, piece_id in (("minkowski", "plane"), ("natario", "plane"), ("lentz", "plane")):
            for x, rho, z in piece(view(metric_id)["surfaces"][0], piece_id):
                near(rho, x, f"{metric_id} rho at {x}")
                self.assertEqual(z, 0, f"{metric_id} z at {x}")
        # Anti-de Sitter's static equator on the hyperboloid Z = sqrt(L^2 + r^2) - L in Minkowski
        # space, and the light cone it nears, Z = rho - L.
        ads = view("anti_de_sitter")
        self.assertEqual(ads["space"], "minkowski")
        for r, rho, z in piece(ads["surfaces"][0], "sheet"):
            near(rho, r, f"anti-de Sitter rho at {r}")
            near(z, math.sqrt(1 + r * r) - 1, f"anti-de Sitter z at {r}")
        for r, rho, z in piece(ads["surfaces"][0], "cone"):
            near(z, rho - 1, f"the light cone at {r}")
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
        for surface in view("kasner")["surfaces"]:
            t = surface["time"]
            for X, Y, Z in surface["curves"][0]["points"]:
                near((X / t ** (-2 / 7)) ** 2 + (Y / t ** (6 / 7)) ** 2, 1, f"Kasner's ring at t = {t}", 1e-5)
        for metric_id in ("bianchi", "pp_wave"):
            for surface in view(metric_id)["surfaces"]:
                P = surface["curves"][0]["points"]
                A, B = max(abs(p[0]) for p in P), max(abs(p[1]) for p in P)
                for X, Y, Z in P:
                    near((X / A) ** 2 + ((Y / B) ** 2 if B else 1 - (X / A) ** 2), 1, f"{metric_id}'s ring", 1e-5)
        self.assertEqual(max(abs(p[1]) for p in view("pp_wave")["surfaces"][-1]["curves"][0]["points"]), 0)

        # Natario's lines of flow, closed and each on one level of the stream function n(r_s) y^2,
        # y being the drawing's X.
        n = lambda r: (math.tanh(4 * (r + 1)) - math.tanh(4 * (r - 1))) / (4 * math.tanh(4))  # noqa: E731
        flows = [c for c in view("natario")["surfaces"][0]["curves"] if c["class"] == "flow"]
        self.assertEqual(len(flows), 6)
        for curve in flows:
            self.assertTrue(curve["closed"])
            levels = [n(math.hypot(X, Y)) * X * X for X, Y, Z in curve["points"]]
            self.assertLess(max(levels) - min(levels), 1e-6 * max(levels), "a line of Natario's flow")

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
                          if any("height" in view for view in data["views"])}, {"alcubierre", "krasnikov"})

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
            self.assertEqual(len(turn["origins"]), len(view["surfaces"]), where)
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
            # Every grid piece names the lines of its grid the figure draws, by index.
            grids = {(k, p["id"]): p["grid"] for k, s in enumerate(view["surfaces"]) for p in s["pieces"] if "grid" in p}
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
                self.assertIn("height", v, f"{v['metric']} {v['view']}")
                continue
            self.assertLess(v["moved"], 1e-3, f"{v['metric']} {v['view']}: {v['moved']:.2e} of the lines moved")

    def test_a_height_over_a_plane_turns_under_the_hand(self):
        # Looked at along the vertical a height hides nothing of itself and its tint covers its
        # rim at the drawn scale; from just above the plane its relief hides part of its grid;
        # and turned all the way round it keeps to its box, the Krasnikov tube's rectangle drawn
        # smaller where it would stand wider or taller than it was published.
        heights = {v["metric"]: v["height"] for v in self.check()["views"] if "height" in v}
        self.assertEqual(set(heights), {"alcubierre", "krasnikov"})
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
        self.assertTrue(".mfs-print-body .mfs-turn-reset { display: none" in page, "print shows the reset button")
        # Both kinds of figure that turn are framed to turn, and the one controller wires them.
        for frame in ("emFigure", "pjFigure"):
            self.assertTrue(re.search(r"function " + frame + r"\([\s\S]*?'mfs-(em|pj)-figure mfs-turn'", page),
                            f"{frame}() does not frame its figure to turn")
        self.assertTrue("function wireTurning(root)" in page, "no wireTurning()")
        self.assertFalse("wireEmbedding" in page, "wireEmbedding() is left")


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
        self.assertEqual(checked, {"alcubierre/bubble", "godel/tipping", "kerr/dragging", "kerr_newman/dragging",
                                   "stockum_dust/tipping"})

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
                             [len(mark["fills"]) for mark in figure["slices"]], where)
            self.assertEqual([len(mark["lines"]) for mark in turn["slices"]],
                             [len(mark["lines"]) for mark in figure["slices"]], where)

    def test_only_the_cosmic_strings_flat_beam_does_not_turn(self):
        # The beam lies in the plane t = 0 seen from straight above, drawn in the plane's own
        # flat coordinates, so it has no other side.
        still = {f"{name}/{figure['id']}" for name, data in diagram_files().items()
                 for figures in data.get("projections", {}).values() for figure in figures if "turn" not in figure}
        self.assertEqual(still, {"cosmic_string/beam"})


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
    HIDDEN = {"frw/comoving_spherical/radial", "frw/comoving_spherical/through", "frw/conformal_spherical/radial",
              "tolman_bondi/comoving_synchronous/collapse", "vaidya/eddington_finkelstein_outgoing/shell",
              "godel/cylindrical/beyond", "stockum_dust/cylindrical/beyond", "conformal frw/flat", "conformal frw/open"}

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
        view = next(v for v in self.embedding[metric_id]["views"] if v["id"] == mark["view"])
        return view, view["surfaces"][mark["surface"]]

    def reach(self, surface, system=None, reference=False):
        xs = [x for piece in surface["pieces"] if "points" in piece and (reference or not piece.get("reference"))
              and (system is None or piece["system"] == system) for x in (piece["points"][0][0], piece["points"][-1][0])]
        return min(xs), max(xs)

    def test_every_moment_appears_on_every_drawing_it_lies_on_and_on_no_other(self):
        drawn = 0
        for metric_id, where, view in self.drawings():
            marks = [(m["view"], m["surface"]) for m in view.get("slices", [])]
            if where in self.HIDDEN:
                self.assertEqual(marks, [], where)
                continue
            every = [(v["id"], i) for v in self.embedding[metric_id]["views"] for i in range(len(v["surfaces"]))]
            self.assertEqual(marks, every, where)
            drawn += len(marks)
        self.assertGreater(drawn, 100)

    def test_every_slice_is_stamped_with_the_surface_it_marks(self):
        for metric_id, where, view in self.drawings():
            for mark in view.get("slices", []):
                target, _ = self.moment(metric_id, mark)
                self.assertEqual(build.embedding_moment_version(target, mark["surface"]), mark["version"], where)
                self.assertTrue(mark["lines"] or mark["points"] or mark["fills"], where)
                self.assertTrue(mark["label"].count("$") % 2 == 0 and mark["label"], where)

    def flat_expected(self, key, view, surface, mark):
        """The moment on a flat view as its drawn axes put it: Y as a function of X, and the
        ends a line of it may have short of the box."""
        t = surface.get("time")
        if key.startswith("schwarzschild/eddington_finkelstein"):
            sign = 1 if "ingoing" in key else -1
            finkelstein = key.endswith("finkelstein")
            return (lambda X: sign * math.log(X - 1) + (0 if finkelstein else sign * X)), None
        if key == "de_sitter/flat_slicing/tx":
            return (lambda X: -0.5 * math.log(1 + X * X)), None
        if key == "pp_wave/exact_plane_wave/tz":
            return (lambda X: X + t), None
        if key == "oppenheimer_snyder/exterior_schwarzschild/radial":
            lo, hi = self.reach(surface, "comoving_synchronous")

            def ct(X):
                R = bisect(lambda R: novikov_t(R, t)[0] - X, lo, hi) if t else X
                return novikov_t(R, t)[1] if t else 0.0
            # The curve runs from the star's surface, the shell released at lo, outward.
            return ct, [novikov_t(lo, t)[0] if t else lo, novikov_t(hi, t)[0] if t else hi]
        if key == "oppenheimer_snyder/interior_comoving/through":
            return (lambda X: t / (2 * math.sqrt(2))), [0, math.pi / 4]
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
        if key.startswith(("godel/cylindrical", "stockum_dust/cylindrical", "minkowski/rindler")):
            return (lambda X: 0.0), None
        if key.startswith("bertotti_robinson"):
            lo, hi = self.reach(surface) if mark["lines"] else (1, 1)
            return (lambda X: 0.0), [lo, hi]
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
                        Y_of, ends = self.flat_expected(key, view, surface, mark)
                        for line in mark["lines"] + [[p] for p in mark["points"]]:
                            for u in line:
                                X, Y = X0 + u[0] * (X1 - X0), Y0 + u[1] * (Y1 - Y0)
                                h = 2e-4 * (X1 - X0)
                                slope = (abs(Y_of(min(X + h, X1)) - Y_of(max(X - h, X0 + 1e-9 if key.startswith(
                                    ("schwarzschild/edd", "oppenheimer_snyder/ext")) else X0))) / (2 * h)
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
        Malament-Hogarth moments ct = tan p + tan q over two, as p, q = arctan(ct -+ r) make
        them; Vaidya's are carried back through each side of the shell's own map."""
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
                    elif metric_id == "malament_hogarth":
                        lo, hi = self.reach(surface)
                        for X, T in points:
                            p, q = (T - X) / 2, (T + X) / 2
                            tp, tq = math.tan(p), math.tan(q)
                            scale = 1 + tp * tp + tq * tq
                            self.assertLess(abs((tp + tq) / 2 - t), 2e-4 * scale, where)
                            self.assertLessEqual(lo - 2e-4 * scale, (tq - tp) / 2, where)
                            self.assertLessEqual((tq - tp) / 2, hi + 2e-4 * scale, where)
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
                    elif metric_id == "rn_metric" and mark["view"] == "inside":
                        self.assertTrue(all(abs(T - math.pi) < 2e-4 for _, T in points), where)
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
        moved["kasner"]["views"][0]["surfaces"][1]["time"] = 0.6
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

    def test_every_metric_lists_its_references_in_the_order_its_history_first_cites_them(self):
        """The page numbers a citation by its place in `references`, so [1] is the first one read."""
        for metric in build.load_metrics():
            cited = []
            for bracket in re.findall(r"\[([A-Za-z0-9_, ]+)\]", metric.get("history") or ""):
                for key in (k.strip() for k in bracket.split(",")):
                    if key not in cited:
                        cited.append(key)
            self.assertEqual(metric.get("references", []), cited, metric["id"])

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


if __name__ == "__main__":
    unittest.main()
