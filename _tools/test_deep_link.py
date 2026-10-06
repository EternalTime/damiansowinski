"""A spacetime opened at its own address, /MFS/?spacetime=<id>, comes up without a page error.

On 6 October 2026 such a link threw before the search's keywords were set up and the page
never came up. The site is built with Jekyll into a folder of its own, served on a free port,
and opened in headless Chrome by `_tools/deep_link.mjs`, which fails on any page error.
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

SPACETIMES = ("kerr", "godel")


class DeepLink(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        for tool in ("bundle", "node"):
            if shutil.which(tool) is None:
                raise AssertionError(f"{tool} is not on the PATH")
        cls.home = tempfile.TemporaryDirectory(prefix="mfs-deep-link-")
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

    def test_a_spacetime_comes_up_at_its_own_address(self):
        for spacetime in SPACETIMES:
            with self.subTest(spacetime=spacetime):
                run = subprocess.run(["node", str(build.ROOT / "_tools" / "deep_link.mjs"), self.base, spacetime],
                                     capture_output=True, text=True, timeout=120)
                self.assertEqual(run.returncode, 0, run.stdout + run.stderr)


class Quiet(http.server.SimpleHTTPRequestHandler):
    def log_message(self, *args):
        pass


if __name__ == "__main__":
    unittest.main()
