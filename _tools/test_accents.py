"""Accents on the spacetimes page and the publications page, the captain's order of 15 October 2026.

A search finds a spacetime whether its accents are typed or not, as the application's
alphabeticalKey folds them: "natario" finds Natário, "godel" Gödel, "lemaitre" Lemaître and
"nordstrom" Nordström. Every word written with a TeX accent in assets/data/references.bib shows
as that letter with that accent in every reference a spacetime lists and every publication,
which once read \\'e as an umlaut, so Poincaré came out Poincarë, and until this order printed
Stojkovi\\'c, Bizo\\'n, Chru\\'sciel and voja\\^go as typed. Both pages, and the applets' notes,
read the bibliography with _includes/unlatex.js, which build_agent_data.unlatex follows for the
agent layer.

The site is built with Jekyll into a folder of its own, served on a free port, and read in
headless Chrome by `_tools/accents.mjs`, which this file tells what to look for.
"""

import functools
import http.server
import json
import re
import shutil
import subprocess
import tempfile
import threading
import unicodedata
import unittest
from pathlib import Path

import build_agent_data as agent
import build_mfs_data as build

READER = build.ROOT / "_includes" / "unlatex.js"
PAGES = (build.ROOT / "_layouts" / "mfs.html", build.ROOT / "publications.markdown",
         build.ROOT / "assets" / "js" / "applet-shell.js")

# The searches the order named, each spelt without its accents first, and the spacetime each
# has to list.
SEARCHES = (
    (("natario", "Natário", "NATÁRIO"), "natario"),
    (("godel", "Gödel", "GÖDEL"), "godel"),
    (("lemaitre", "Lemaître"), "frw"),
    (("nordstrom", "Nordström", "nordström"), "nordstrom_scalar"),
    (("malament-hogarth", "Malament-Hogarth"), "malament_hogarth"),
)
# Keywords with accents, typed without them.
KEYWORDS = (("kahler", "Kähler"), ("painleve", "Painlevé-Gullstrand coordinates"))

# Each TeX accent by the name Unicode gives a letter carrying it, read from the character names
# rather than composed, so the expectation is not the readers' own arithmetic.
ACCENT_NAMES = {
    "'": "ACUTE", '"': "DIAERESIS", "`": "GRAVE", "^": "CIRCUMFLEX", "~": "TILDE", "=": "MACRON",
    ".": "DOT ABOVE", "u": "BREVE", "v": "CARON", "H": "DOUBLE ACUTE", "c": "CEDILLA",
    "k": "OGONEK", "r": "RING ABOVE",
}
ACCENT = re.compile(r"""\\(?:(['"`^~=.])|([uvHckr])(?=\s*\{|\s+[A-Za-z\\]))\s*"""
                    r"""(?:\{\s*(\\[ij](?![A-Za-z])|[A-Za-z])\s*\}|(\\[ij](?![A-Za-z])|[A-Za-z]))""")
TEX_ACCENT = re.compile(r"""\\(['"`^~=.]|[uvHckr]\s*\{)""")
# The fields a reference and a publication show with words in them.
SHOWN = ("author", "title")
VENUES = ("journal", "booktitle", "publisher")

# Names the order named, as a reader reads them, which the oracle below has to agree with.
NAMED = {"Nat\\'ario": "Natário", "Jos\\'e": "José", "Poincar\\'e": "Poincaré", "G\\\"odel": "Gödel",
         "Stojkovi\\'c": "Stojković", "Bizo\\'n": "Bizoń", "Chru\\'sciel": "Chruściel",
         "voja\\^go": "vojaĝo", "Fran\\c{c}ois": "François", "\\v{C}ech": "Čech", "\\'{\\i}": "í",
         "Sto{\\ss}vorg\\\"ange": "Stoßvorgänge"}


def accented(text):
    """The text with each TeX accent written as the letter carrying it, \\ss as ß, TeX's dashes and
    quotes written out as _tools/test_reference_tex.py holds them, and its braces gone."""
    def letter(match):
        sign, word, braced, bare = match.groups()
        base = (braced or bare).lstrip("\\")
        case = "CAPITAL" if base.isupper() else "SMALL"
        return unicodedata.lookup(f"LATIN {case} LETTER {base.upper()} WITH {ACCENT_NAMES[sign or word]}")
    text = re.sub(r"\\ss\b", "ß", ACCENT.sub(letter, text))
    for tex, written in (("---", "\u2013"), ("--", "\u2013"), ("``", '"'), ("''", '"'), ("`", "'")):
        text = text.replace(tex, written)
    return text.replace("{", "").replace("}", "")


def accented_words(fields):
    """The words a reference shows that the bibliography writes with a TeX accent."""
    venue = next((fields[v] for v in VENUES if fields.get(v)), "")
    words = []
    for value in [fields.get(f, "") for f in SHOWN] + [venue]:
        for token in value.split():
            if "$" not in token and ACCENT.search(token):
                words.append(accented(token).strip("()[],.;:!?\""))
    return words


def run_reader(values):
    """_includes/unlatex.js run in Node over each value."""
    script = ("const src = require('fs').readFileSync(process.argv[1], 'utf8');"
              "const unlatex = new Function(src + '; return unlatex;')();"
              "process.stdout.write(JSON.stringify(JSON.parse(require('fs').readFileSync(0, 'utf8')).map(unlatex)));")
    run = subprocess.run(["node", "-e", script, str(READER)], input=json.dumps(values),
                         capture_output=True, text=True, check=True)
    return json.loads(run.stdout)


class TheReader(unittest.TestCase):
    """The reader every page showing the bibliography includes, over every field of it."""

    @classmethod
    def setUpClass(cls):
        cls.entries = build.parse_bibtex(build.BIB_FILE.read_text(encoding="utf-8"))
        cls.values = [v for e in cls.entries.values() for v in e["fields"].values()]

    def test_the_oracle_reads_the_names_as_the_captain_does(self):
        for tex, reads in NAMED.items():
            self.assertEqual(accented(tex), reads)

    def test_every_page_includes_it_and_reads_no_accent_of_its_own(self):
        for page in PAGES:
            source = page.read_text(encoding="utf-8")
            with self.subTest(page.name):
                self.assertEqual(source.count("{% include unlatex.js %}"), 1)
                self.assertNotIn("function unlatex", source)

    def test_it_reads_the_named_accents(self):
        self.assertEqual(run_reader(list(NAMED)), list(NAMED.values()))

    def test_it_leaves_no_accent_as_tex_and_puts_each_on_its_letter(self):
        read = run_reader(self.values)
        for value, shown in zip(self.values, read):
            if not ACCENT.search(value):
                continue
            with self.subTest(value):
                self.assertIsNone(TEX_ACCENT.search(shown))
                for token in value.split():
                    if "$" not in token and ACCENT.search(token):
                        self.assertIn(accented(token).strip("()[],.;:!?\""), shown)

    def test_the_agent_layer_reads_every_field_as_the_pages_do(self):
        for value, shown in zip(self.values, run_reader(self.values)):
            with self.subTest(value):
                self.assertEqual(agent.unlatex(value), " ".join(shown.split()))


class BuiltSite:
    """The site built with Jekyll into a folder of its own and served on a free port at `base`."""

    @classmethod
    def setUpClass(cls):
        for tool in ("bundle", "node"):
            if shutil.which(tool) is None:
                raise AssertionError(f"{tool} is not on the PATH")
        cls.home = tempfile.TemporaryDirectory(prefix="mfs-site-")
        site = Path(cls.home.name) / "site.noindex"
        built = subprocess.run(["bundle", "exec", "jekyll", "build", "-q", "-d", str(site)],
                               cwd=build.ROOT, capture_output=True, text=True, timeout=600)
        if built.returncode != 0:
            cls.home.cleanup()
            raise AssertionError(built.stderr[-3000:])
        handler = functools.partial(Quiet, directory=str(site))
        cls.server = Server(("127.0.0.1", 0), handler)
        threading.Thread(target=cls.server.serve_forever, daemon=True).start()
        cls.base = f"http://127.0.0.1:{cls.server.server_address[1]}"

    @classmethod
    def tearDownClass(cls):
        cls.server.shutdown()
        cls.server.server_close()
        cls.home.cleanup()


class InTheBrowser(BuiltSite, unittest.TestCase):

    def spec(self):
        entries = build.parse_bibtex(build.BIB_FILE.read_text(encoding="utf-8"))
        index = json.loads(build.INDEX_FILE.read_text(encoding="utf-8"))
        listed = {e["id"] for e in index}
        # Every listed spacetime citing an entry with an accent, fewest pages first: each spacetime
        # in turn that shows the most of the entries not yet seen.
        citing = {m["id"]: {k for k in m.get("references", []) if accented_words(entries[k]["fields"])}
                  for m in build.load_metrics() if m["id"] in listed}
        left, spacetimes = set().union(*citing.values()), {}
        while left:
            best = max(sorted(citing), key=lambda i: len(citing[i] & left))
            spacetimes[best] = {k: accented_words(entries[k]["fields"]) for k in sorted(citing[best] & left)}
            left -= citing[best]
        keys = json.loads((build.ROOT / "_data" / "publication_keys.json").read_text(encoding="utf-8"))
        return {
            "search": [{"spellings": list(s), "finds": finds} for s, finds in SEARCHES],
            "keywords": [{"typed": typed, "offers": keyword,
                          "carrying": [e["id"] for e in index if keyword in e["tags"]]}
                         for typed, keyword in KEYWORDS],
            "spacetimes": spacetimes,
            "publications": [accented_words(entries[k]["fields"]) for k in keys],
        }

    def test_searches_fold_accents_and_every_accent_reads_as_its_letter(self):
        spec = self.spec()
        shown = {k: w for references in spec["spacetimes"].values() for k, w in references.items()}
        self.assertEqual(shown["frolov2003"], ["Stojković"])
        path = Path(self.home.name) / "accents.json"
        path.write_text(json.dumps(spec, ensure_ascii=False), encoding="utf-8")
        run = subprocess.run(["node", str(build.ROOT / "_tools" / "accents.mjs"), self.base, str(path)],
                             capture_output=True, text=True, timeout=600)
        self.assertEqual(run.returncode, 0, run.stdout + run.stderr)


class Server(http.server.ThreadingHTTPServer):
    # A page asks for dozens of files at once, and past the default backlog of five connections
    # macOS resets the rest, which Chrome reports as a console error: 40 clients at once had 380
    # of 400 requests reset at five and none at 128.
    request_queue_size = 128


class Quiet(http.server.SimpleHTTPRequestHandler):
    def log_message(self, *args):
        pass


if __name__ == "__main__":
    unittest.main()
