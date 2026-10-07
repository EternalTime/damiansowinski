"""The page "More from Owl's Nest Creations": .venv.noindex/bin/python -m unittest discover -s _tools

The captain's order of 6 October 2026 is one page of the studio's other applications with his
two lines under them, opened by "Studio" beside the exit sign on the spacetimes page and
beside About at the top of the application's list. The applications are the studio site's own
list of its games, https://owlsnestcreations.com/data/games.json, read as the page opens, so a
release reaches the page with no change here; the application reads the same file. These run
the page's own reading of that list in Node over a list of the studio's shape, and hold the two
lines to his words.
"""
import json
import re
import shutil
import subprocess
import unittest

import build_mfs_data as build

PAGE = build.ROOT / "_layouts" / "mfs.html"
# My Favorite Spacetimes on the App Store, which the page never lists.
OWN_STORE_ID = "6810529195"

# The studio's list as it stood on 6 October 2026, with My Favorite Spacetimes as it will stand
# once it is released, and games a reader cannot get.
CATALOGUE = {
    "name": "Owl's Nest Creations games",
    "games": [
        {"name": "VoidFlux", "tagline": "A Puzzle in the depths of Electrostatics.", "description": "A puzzle game.",
         "released": True, "releaseDate": "2026-09-04", "platforms": ["iPhone", "iPad"],
         "image": "https://owlsnestcreations.com/assets/img/voidflux-icon.png",
         "pageUrl": "https://owlsnestcreations.com/games/voidflux/", "pressKitUrl": None,
         "appStoreUrl": "https://apps.apple.com/us/app/voidflux/id6778287432"},
        {"name": "The Verdant Engine", "tagline": "A serene Idle-RPG set in a world powered by Pattern Formation.",
         "released": False, "releaseNote": "Expected November release", "platforms": [], "image": None, "appStoreUrl": None},
        {"name": "My Favorite Spacetimes", "tagline": "A reference tool for General Relativity.", "released": True,
         "image": "https://owlsnestcreations.com/assets/img/mfs-icon.png",
         "appStoreUrl": f"https://apps.apple.com/us/app/my-favorite-spacetimes/id{OWN_STORE_ID}"},
        {"name": "Unreleased with an address", "tagline": "A line.", "released": False,
         "appStoreUrl": "https://apps.apple.com/app/id11"},
        {"name": "Released with no address", "tagline": "A line.", "released": True, "appStoreUrl": None},
        {"name": "Released elsewhere", "tagline": "A line.", "released": True, "appStoreUrl": "https://example.com/app/id12"},
        {"name": "", "tagline": "A line.", "released": True, "appStoreUrl": "https://apps.apple.com/app/id13"},
        {"name": "No picture", "tagline": "Its line.", "released": True, "image": "http://example.com/a.png",
         "appStoreUrl": "https://apps.apple.com/app/id14"},
    ],
}


def page_function(source, name):
    found = re.search(r"      function " + name + r"\(.*?\n      \}\n", source, re.S)
    if not found:
        raise AssertionError(f"the page no longer defines {name}")
    return found.group(0)


class ThePage(unittest.TestCase):
    """The page reads the studio's list as the application does, and carries the captain's lines."""

    @classmethod
    def setUpClass(cls):
        cls.source = PAGE.read_text(encoding="utf-8")

    def run_page(self, call, payload):
        if shutil.which("node") is None:
            self.skipTest("node is not installed")
        script = "".join([
            re.search(r"      var OWN_STORE_ID = .*?\n", self.source).group(0),
            page_function(self.source, "storeIdIn"),
            page_function(self.source, "studioApps"),
            page_function(self.source, "storeAddress"),
            "const input = JSON.parse(require('fs').readFileSync(0, 'utf8'));\n",
            f"process.stdout.write(JSON.stringify({call}));\n",
        ])
        run = subprocess.run(["node", "-e", script], input=json.dumps(payload), capture_output=True, text=True, timeout=60)
        self.assertEqual(run.returncode, 0, run.stderr)
        return json.loads(run.stdout)

    def test_studio_stands_before_the_exit_sign(self):
        self.assertRegex(self.source, r'<button type="button" id="mfs-studio-sign"[^>]*>Studio</button>\s*<a id="mfs-exit"')

    def test_the_list_is_the_studio_sites_own(self):
        self.assertIn("var STUDIO_LIST = 'https://owlsnestcreations.com/data/games.json';", self.source)
        self.assertIn(f"var OWN_STORE_ID = '{OWN_STORE_ID}';", self.source)

    def test_only_released_games_on_the_app_store_show_and_never_this_application(self):
        self.assertEqual(self.run_page("studioApps(input)", CATALOGUE), [
            {"store_id": "6778287432", "name": "VoidFlux", "line": "A Puzzle in the depths of Electrostatics.",
             "icon": "https://owlsnestcreations.com/assets/img/voidflux-icon.png"},
            {"store_id": "14", "name": "No picture", "line": "Its line.", "icon": None},
        ])

    def test_a_list_that_is_not_the_studios_shows_no_application(self):
        self.assertEqual(self.run_page("[studioApps(input[0]), studioApps(input[1]), studioApps(input[2])]",
                                       [None, [], {"name": "no games"}]), [[], [], []])

    def test_a_store_identifier_is_the_digits_after_id_and_its_address_carries_no_storefront(self):
        addresses = ["https://apps.apple.com/us/app/voidflux/id6778287432", "https://apps.apple.com/app/id42",
                     "https://itunes.apple.com/app/id7", "https://apps.apple.com/us/app/voidflux/id",
                     "https://apps.apple.com/us/app/voidflux/id12a", "http://apps.apple.com/app/id1",
                     "https://example.com/app/id1", None]
        self.assertEqual(self.run_page("input.map(storeIdIn)", addresses),
                         ["6778287432", "42", "7", None, None, None, None, None])
        self.assertEqual(self.run_page("storeAddress({store_id: '6778287432'})", None),
                         "https://apps.apple.com/app/id6778287432")

    def test_the_captains_two_lines_stand_under_the_applications_in_his_words(self):
        links = re.search(r"var STUDIO_LINKS = \[\n(.*?)\n      \];", self.source, re.S).group(1)
        self.assertEqual(re.findall(r"\{ words: \"(.*?)\", address: '(.*?)' \}", links), [
            ("See what's hatching in the nest...", "https://owlsnestcreations.com"),
            ("Check out Damian's research page...", "https://damiansowinski.com"),
        ])


if __name__ == "__main__":
    unittest.main()
