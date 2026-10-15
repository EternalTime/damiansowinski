"""Apple's App Store badge on the spacetimes page, the captain's order of 15 October 2026.

The badge links to the application's page on the App Store and stands where a visitor sees it
without hunting: before the exit sign in the title's band on a desktop, in a row of its own
under the title and the sign on a phone. It is Apple's artwork as the studio's site uses it,
unaltered. The site is built with Jekyll into a folder of its own, served on a free port, and
read in headless Chrome by `_tools/app_store_badge.mjs` on a wide desktop, a narrow one where
the title wraps to leave it room, and two phones.
"""

import functools
import http.server
import shutil
import subprocess
import tempfile
import threading
import unittest
from pathlib import Path

import build_mfs_data as build

BADGE = build.ROOT / "assets" / "images" / "app_store_button_white.svg"
SIZES = (("1440x900", False), ("700x800", False), ("390x844", True), ("320x568", True))


class Badge(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        for tool in ("bundle", "node"):
            if shutil.which(tool) is None:
                raise AssertionError(f"{tool} is not on the PATH")
        cls.home = tempfile.TemporaryDirectory(prefix="mfs-app-store-")
        site = Path(cls.home.name) / "site.noindex"
        built = subprocess.run(["bundle", "exec", "jekyll", "build", "-q", "-d", str(site)],
                               cwd=build.ROOT, capture_output=True, text=True, timeout=600)
        if built.returncode != 0:
            cls.home.cleanup()
            raise AssertionError(built.stderr[-3000:])
        handler = functools.partial(Quiet, directory=str(site))
        cls.server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)
        threading.Thread(target=cls.server.serve_forever, daemon=True).start()
        cls.base = f"http://127.0.0.1:{cls.server.server_address[1]}"

    @classmethod
    def tearDownClass(cls):
        cls.server.shutdown()
        cls.server.server_close()
        cls.home.cleanup()

    def test_it_is_apple_s_artwork(self):
        svg = BADGE.read_text(encoding="utf-8")
        self.assertIn("<title>Download_on_the_App_Store_Badge_US-UK_RGB_wht_092917</title>", svg)
        self.assertIn('viewBox="0 0 119.66407 40"', svg)

    def test_it_shows_and_links_to_the_application_on_every_size(self):
        for size, phone in SIZES:
            with self.subTest(size=size):
                args = ["node", str(build.ROOT / "_tools" / "app_store_badge.mjs"), self.base, f"--size={size}"]
                run = subprocess.run(args + (["--phone"] if phone else []), capture_output=True, text=True, timeout=120)
                self.assertEqual(run.returncode, 0, run.stdout + run.stderr)


class Quiet(http.server.SimpleHTTPRequestHandler):
    def log_message(self, *args):
        pass


if __name__ == "__main__":
    unittest.main()
