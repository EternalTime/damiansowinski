#!/usr/bin/env python3
"""Tests for the MFS data generator.

    python3 -m unittest discover -s _tools
"""

import copy
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import build_mfs_data as build


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


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
        for entry, metric in zip(self.index, self.metrics):
            self.assertEqual(set(entry), {"id", "name", "tags", "version"})
            self.assertEqual(entry["name"], metric["short_name"])
            self.assertEqual(entry["tags"], metric["tags"])
            self.assertEqual(entry["version"], build.content_version(metric))

    def test_a_metric_missing_its_short_name_is_refused(self):
        with self.assertRaises(build.DataError):
            self.load_one("x.json", {"id": "x", "name": "X", "tags": ["t"]})

    def test_a_metric_whose_id_is_not_its_filename_is_refused(self):
        with self.assertRaises(build.DataError):
            self.load_one("x.json", {"id": "y", "name": "X", "short_name": "X", "tags": ["t"]})

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
        with tempfile.TemporaryDirectory() as folder:
            for filename, metric in files.items():
                (Path(folder) / filename).write_text(json.dumps(metric), encoding="utf-8")
            with mock.patch.object(build, "METRICS_DIR", Path(folder)):
                return build.load_metrics()


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
