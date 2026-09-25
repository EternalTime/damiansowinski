#!/usr/bin/env python3
"""Tests for the MFS data generator.

    python3 -m unittest discover -s _tools
"""

import contextlib
import copy
import io
import json
import re
import tempfile
import unicodedata
import unittest
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


def diagram_prose(name, diagram):
    """Yield every sentence carrying field of a diagram file, each with the name of its place."""
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
        diagrams, conformal = diagram_files(), conformal_files()
        for entry, metric in zip(self.index, self.metrics):
            expected = ({"id", "name", "tags", "version"} | ({"diagrams"} if metric["id"] in diagrams else set())
                        | ({"conformal"} if metric["id"] in conformal else set()))
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

    def test_no_metric_on_disk_spells_a_character_as_an_escape(self):
        # A \u escaped twice in the JSON reaches the page as the six characters of the
        # escape, as the î of Lemaître did in the Tolman-Bondi conventions, since neither
        # TeX nor the page reads it. The character itself is what the files carry.
        for metric in build.load_metrics():
            for field, value in prose(metric):
                self.assertNotRegex(value, r"\\u[0-9a-fA-F]{4}",
                                    f"{metric['id']}.json: {field} spells a character as an escape")

    def test_every_caption_names_its_plane_and_never_a_block(self):
        """A reader is told which slice is drawn, in words, and never handed the shorthand."""
        views = 0
        for name, diagram in diagram_files().items():
            for system, drawn in diagram["systems"].items():
                for view in drawn:
                    views += 1
                    where = f"diagrams/{name}.json {system}/{view['id']}"
                    self.assertTrue(view["caption"][0].startswith("This is the "), where)
                    for field, value in diagram_prose(name, {"systems": {system: [view]}}):
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
        self.assertTrue(fields, "no prose was read, so nothing was checked")
        for where, field, value in fields:
            for pattern in MACHINERY:
                self.assertNotRegex(value, pattern, where)
            if field != "history":
                # A history tells of papers that were published and printed.
                for pattern in MACHINERY_OUTSIDE_A_HISTORY:
                    self.assertNotRegex(value, pattern, where)

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
        self.assertIsNone(build.check_history_shape(metrics))

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
        self.assertIsNone(build.check_history_shape([self.metric(self.history(3, 6, 3, 6, 3))]))

    def test_a_history_of_four_paragraphs_is_refused(self):
        with self.assertRaises(build.DataError) as raised:
            build.check_history_shape([self.metric(self.history(4, 4, 4, 4))])
        self.assertIn("4 paragraphs", str(raised.exception))

    def test_a_paragraph_too_short_or_too_long_is_refused_by_its_place(self):
        for counts, place in (((4, 4, 2, 4, 4), 3), ((4, 4, 4, 4, 7), 5)):
            with self.assertRaises(build.DataError) as raised:
                build.check_history_shape([self.metric(self.history(*counts))])
            self.assertIn(f"paragraph {place}", str(raised.exception))

    def test_the_longest_paragraph_may_be_at_most_twice_the_shortest(self):
        with mock.patch.object(build, "HISTORY_SENTENCES", (1, 10)):
            self.assertIsNone(build.check_history_shape([self.metric(self.history(2, 4, 3, 4, 2))]))
            with self.assertRaises(build.DataError) as raised:
                build.check_history_shape([self.metric(self.history(2, 5, 3, 4, 2))])
        self.assertIn("more than 2 times", str(raised.exception))

    def test_a_table_is_not_a_paragraph(self):
        history = self.history(3, 3, 3, 3) + "¶TABLE:: a | b ;; 1 | 2¶" + self.paragraph(3)
        self.assertEqual(build.history_shape(history), [3, 3, 3, 3, 3])
        self.assertIsNone(build.check_history_shape([self.metric(history)]))

    def test_every_history_out_of_shape_is_named_at_once(self):
        with self.assertRaises(build.DataError) as raised:
            build.check_history_shape([self.metric(self.history(3, 3), "one"),
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

    def test_every_caption_names_what_is_drawn(self):
        for metric_id, data in self.conformal.items():
            for view in data["views"]:
                self.assertTrue(view["caption"], f"{metric_id} {view['id']}")
                self.assertTrue(view["caption"][0].startswith("This is the "), f"{metric_id} {view['id']}")

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
