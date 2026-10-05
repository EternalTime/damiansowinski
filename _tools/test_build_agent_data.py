#!/usr/bin/env python3
"""Tests for the agent layer generator.

    .venv.noindex/bin/python -m unittest discover -s _tools
"""

import contextlib
import io
import json
import re
import unittest

import build_agent_data as agent
import build_mfs_data as mfs


def generated(path):
    return json.loads(path.read_text(encoding="utf-8"))


class PublishedFilesAreCurrent(unittest.TestCase):
    def test_check_mode_passes(self):
        err = io.StringIO()
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(err):
            status = agent.main(["--check"])
        self.assertEqual(status, 0, err.getvalue())


class Spacetimes(unittest.TestCase):
    def setUp(self):
        self.data = generated(agent.SPACETIMES_FILE)

    def test_every_metric_file_is_one_spacetime(self):
        ids = sorted(path.stem for path in mfs.METRICS_DIR.glob("*.json")
                     if not mfs.CONFLICT_COPY.search(path.stem))
        self.assertEqual(sorted(s["id"] for s in self.data["spacetimes"]), ids)
        self.assertEqual(self.data["count"], len(ids))

    def test_in_the_order_of_the_search_list(self):
        index = json.loads(mfs.INDEX_FILE.read_text(encoding="utf-8"))
        self.assertEqual([s["id"] for s in self.data["spacetimes"]], [e["id"] for e in index])

    def test_each_carries_its_metric_and_where_to_find_it(self):
        for s in self.data["spacetimes"]:
            with self.subTest(s["id"]):
                self.assertTrue(s["name"] and s["description"] and s["tags"])
                self.assertEqual(s["page_url"], "https://damiansowinski.com/MFS/")
                self.assertTrue((mfs.ROOT / s["data_url"].removeprefix(agent.SITE_URL + "/")).is_file())
                for chart in s["charts"]:
                    self.assertTrue(chart["coordinates"] and chart["line_element"] and chart["metric_components"])

    def test_each_chart_carries_its_own_convention_beside_the_shared_one(self):
        metrics = {m["id"]: m for m in mfs.load_metrics()}
        for s in self.data["spacetimes"]:
            with self.subTest(s["id"]):
                metric = metrics[s["id"]]
                self.assertEqual(s["convention"], metric["convention"])
                self.assertEqual([c["convention"] for c in s["charts"]],
                                 [c["convention"] for c in metric["coordinates"]])


class Publications(unittest.TestCase):
    def setUp(self):
        self.data = generated(agent.PUBLICATIONS_FILE)
        self.keys = generated(agent.DATA_DIR / "publication_keys.json")

    def test_in_the_order_the_page_lists_them(self):
        self.assertEqual([p["key"] for p in self.data["publications"]], self.keys)

    def test_each_has_what_an_agent_cites_by(self):
        for p in self.data["publications"]:
            with self.subTest(p["key"]):
                self.assertTrue(p["title"] and p["authors"] and p["venue"] and p["url"])
                self.assertIsInstance(p["year"], int)
                for author in p["authors"]:
                    self.assertNotRegex(author, r"[\\{}]|,")

    def test_every_pdf_on_the_site_exists(self):
        for p in self.data["publications"]:
            if p.get("pdf_url", "").startswith(agent.SITE_URL):
                with self.subTest(p["key"]):
                    self.assertTrue((mfs.ROOT / p["pdf_url"].removeprefix(agent.SITE_URL + "/")).is_file())

    def test_a_key_missing_from_the_bibliography_is_refused(self):
        with self.assertRaisesRegex(mfs.DataError, "nosuchkey"):
            agent.build_publications({}, ["nosuchkey"])

    def test_accents_authors_and_identifiers(self):
        entry = {"type": "article", "fields": {
            "title": "{Einstein}'s field", "author": "Pi\\~nero, Jordi and D\\\"obler, Niklas and others",
            "journal": "arXiv preprint arXiv:2503.02980", "year": "2025"}}
        p = agent.publication_entry("k", entry)
        self.assertEqual(p["title"], "Einstein's field")
        self.assertEqual(p["authors"], ["Jordi Piñero", "Niklas Döbler", "et al."])
        self.assertEqual(p["arxiv_url"], "https://arxiv.org/abs/2503.02980")
        self.assertEqual(p["url"], p["arxiv_url"])


class FullText(unittest.TestCase):
    def test_markup_becomes_markdown_with_absolute_links(self):
        source = (
            "---\ntitle: T\n---\n"
            "<a href=\"{{ '/' | relative_url }}\">← Back to Home</a>\n"
            "<h1 class=\"x\">Heading</h1>\n"
            "<script>var hidden = 1;</script>\n"
            "Text with <a href=\"{{ '/notes/' | relative_url }}\">notes </a>and an icon\n"
            "<svg><path d=\"M0\"/></svg>\nhere.\n"
            "<ul>\n  <li><a href=\"https://example.org\">One</a></li>\n</ul>\n"
            "<button class=\"b\" onclick=\"x()\">Ising Model</button>\n"
        )
        text = agent.page_markdown(source)
        self.assertNotIn("Back to Home", text)
        self.assertNotIn("hidden", text)
        self.assertIn("### Heading", text)
        self.assertIn("[notes](https://damiansowinski.com/notes/) and an icon here.", text)
        self.assertIn("- [One](https://example.org)", text)
        self.assertIn("- Ising Model (interactive applet)", text)
        self.assertNotRegex(text, r"<[a-z/]")

    def test_the_published_file_has_every_page_and_no_markup(self):
        text = agent.FULL_TEXT_FILE.read_text(encoding="utf-8")
        for _, title, url, _ in agent.PAGES:
            self.assertIn(f"## {title}\n\nURL: {agent.absolute(url)}", text)
        self.assertNotRegex(text, r"\{\{|\{%|<(div|span|script|svg|a )")
        self.assertNotIn("←", text)
        self.assertTrue(re.match(r"# .+\n\n> .+", text))


if __name__ == "__main__":
    unittest.main()
