"""The text a view of a drawing shows beside it, as `derivations/view_text.py` counts it for
every drawing program: .venv.noindex/bin/python -m unittest discover -s _tools

The count is held to views made up here, one of each shape the four drawing programs write, and
the warning to a view past the limit and to one within it. Every view of every published file is
counted, so that no program writes a view the count cannot read.
"""
import contextlib
import io
import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))

import view_text  # noqa: E402


def warnings(where, views):
    err = io.StringIO()
    with contextlib.redirect_stderr(err):
        view_text.warn(where, views)
    return err.getvalue()


class Count(unittest.TestCase):
    def test_a_curvature_view_counts_its_caption_clip_settings_and_the_text_of_its_legend(self):
        view = {"id": "K", "caption": ["abc", "de"], "clip": "fgh", "settings": "ij",
                "legend": [["line", "K", "klmn"], ["fill", "horizon", "o"]]}
        self.assertEqual(view_text.length(view), 3 + 2 + 3 + 2 + 4 + 1)

    def test_an_embedding_view_counts_the_legend_of_its_figure(self):
        view = {"id": "equator", "caption": ["abc"], "settings": "de",
                "figure": {"legend": [["line", "throat", "fghi"]]}}
        self.assertEqual(view_text.length(view), 3 + 2 + 4)

    def test_a_view_with_no_settings_counts_none(self):
        view = {"id": "plane", "caption": ["abc"], "settings": None, "legend": []}
        self.assertEqual(view_text.length(view), 3)

    def test_a_null_ray_view_counts_its_families_its_cones_and_each_kind_of_marker_once(self):
        view = {"id": "equator", "caption": ["abc"], "settings": "de", "families": ["ingoing", "outgoing"],
                "cones": [{"kind": "future"}], "cone": "cone",
                "markers": [{"kind": "shell", "legend": "the shell"}, {"kind": "shell", "legend": "the shell"},
                            {"kind": "grr"}]}
        self.assertEqual(view_text.length(view), 3 + 2 + 7 + 8 + 4 + 9)
        view["cones"] = []
        self.assertEqual(view_text.length(view), 3 + 2 + 7 + 8 + 9)


class Warn(unittest.TestCase):
    def view(self, n):
        return {"id": "K", "caption": ["x" * n], "clip": "", "settings": "", "legend": []}

    def test_a_view_past_the_limit_is_named_with_its_count(self):
        out = warnings("some_spacetime curvature", [self.view(view_text.LIMIT + 1)])
        self.assertIn("some_spacetime curvature K", out)
        self.assertIn(f"{view_text.LIMIT + 1} characters", out)

    def test_a_view_within_the_limit_says_nothing(self):
        self.assertEqual(warnings("some_spacetime curvature", [self.view(view_text.LIMIT)]), "")


class Published(unittest.TestCase):
    def test_every_published_view_is_counted(self):
        for folder in ("curvature", "conformal", "embedding", "diagrams"):
            for path in sorted((DATA / folder).glob("*.json")):
                data = json.loads(path.read_text(encoding="utf-8"))
                groups = [data["views"]] if "views" in data else [
                    *data["systems"].values(), *data.get("projections", {}).values()]
                for views in groups:
                    for view in views:
                        self.assertGreater(view_text.length(view), 0, f"{path.name} {view['id']}")
