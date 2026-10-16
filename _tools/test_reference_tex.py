"""The TeX in assets/data/references.bib as a reader sees it, fixed on 16 October 2026.

The Gott-Li spacetime's 13th reference had a TeX error in it: its title is Elliptic de Sitter
space: $dS/\\mathbb{Z}_2$, and _includes/unlatex.js took the braces off inside the mathematics as
well as outside it, so MathJax was handed \\mathbbZ_2 and printed an undefined command in red.
Every title whose mathematics carries a brace was wrong the same way, $10^{-13}$ came out 10^-13
and $\\mathcal{N} = 4$ an error, and outside the mathematics TeX's dashes and quotes were printed
as typed. The reader now leaves mathematics as written, for MathJax, and writes everything else
out, on the spacetimes page, the publications page and the applets' notes alike.

The reader is run over every field of the bibliography in Node, and the site is built with Jekyll
into a folder of its own, served on a free port and read in headless Chrome by
`_tools/reference_tex.mjs`, which this file tells what to look for.
"""

import json
import re
import subprocess
import unittest
from pathlib import Path

import build_agent_data as agent
import build_mfs_data as build
import test_accents as accents

MATH = re.compile(r"(\$[^$]*\$)")
# TeX a reader must never see outside mathematics: a command, a dollar sign, a brace, TeX's quotes
# and its dashes. reference_tex.mjs looks for the same in the page.
TEX = re.compile(r"\\|\$|[{}]|``|''|`|--")
# The fields a reference shows, and those that are addresses rather than words.
SHOWN = ("author", "title", "journal", "booktitle", "publisher", "volume", "pages", "number")
VENUES = ("journal", "booktitle", "publisher")
ADDRESSES = ("url", "doi", "eprint")

# TeX the reader writes out, as a reader reads it.
WRITTEN = {
    "Elliptic de Sitter space: $dS/\\mathbb{Z}_2$": "Elliptic de Sitter space: $dS/\\mathbb{Z}_2$",
    "A Potential New Vector Boson of $8.7 \\times 10^{-13}$ {eV}": "A Potential New Vector Boson of $8.7 \\times 10^{-13}$ eV",
    "{Lense}--{Thirring}": "Lense\u2013Thirring",
    "The Thermodynamics of Computation---a Review": "The Thermodynamics of Computation\u2013a Review",
    "the game ``life''": 'the game "life"',
    "The unphysical nature of `warp drive'": "The unphysical nature of 'warp drive'",
    "G\\\"odel-type universes in $f(R)$ gravity": "Gödel-type universes in $f(R)$ gravity",
}


def mathematics(text):
    return MATH.findall(text)


def words(text):
    return [w for w in (" ".join(p.split()) for p in MATH.split(text)[::2]) if w]


def expected(fields, shown=SHOWN):
    """What a reference shows: its formulas, counted, and the words of its title as read."""
    venue = next((v for v in VENUES if fields.get(v)), None)
    counted = [f for f in shown if f not in VENUES or f == venue]
    return {"formulas": sum(len(mathematics(fields.get(f, ""))) for f in counted),
            "words": words(agent.unlatex(fields.get("title", "")))}


def carries_tex(fields):
    """Whether a reference shows mathematics, or TeX besides an accent, \\ss, \\& and a brace."""
    for name in SHOWN:
        value = fields.get(name, "")
        plain = re.sub(r"\\ss\b", "", accents.ACCENT.sub("", MATH.sub("", value))).replace("\\&", "")
        if MATH.search(value) or re.search(r"\\|``|''|`|--", plain):
            return True
    return False


def entry_texts(text):
    """Each entry of the bibliography as written, by its key, split as the pages split it."""
    texts = {}
    for chunk in re.split(r"(?=@\w+\s*\{)", text):
        key = re.match(r"@\w+\s*\{\s*([\w:]+)\s*,", chunk.strip())
        if key:
            texts[key.group(1)] = chunk.strip()
    return texts


class TheReader(unittest.TestCase):
    """_includes/unlatex.js over every field of the bibliography."""

    @classmethod
    def setUpClass(cls):
        cls.entries = build.parse_bibtex(build.BIB_FILE.read_text(encoding="utf-8"))
        cls.fields = [(k, f, v) for k, e in cls.entries.items() for f, v in e["fields"].items()]
        cls.read = accents.run_reader([v for _, _, v in cls.fields])

    def test_it_writes_out_the_tex_named(self):
        self.assertEqual(accents.run_reader(list(WRITTEN)), list(WRITTEN.values()))
        for tex, reads in WRITTEN.items():
            self.assertEqual(agent.unlatex(tex), reads)

    def test_every_dollar_sign_opens_or_closes_mathematics(self):
        for key, field, value in self.fields:
            with self.subTest(f"{key} {field}"):
                self.assertEqual(value.count("$") % 2, 0, value)

    def test_it_leaves_mathematics_as_written(self):
        for (key, field, value), shown in zip(self.fields, self.read):
            with self.subTest(f"{key} {field}"):
                self.assertEqual(mathematics(shown), mathematics(value))

    def test_it_leaves_no_tex_outside_mathematics(self):
        for (key, field, value), shown in zip(self.fields, self.read):
            if field in ADDRESSES:
                continue
            with self.subTest(f"{key} {field}"):
                self.assertIsNone(TEX.search(MATH.sub("", shown)), shown)

    def test_it_reads_every_address_as_written(self):
        for (key, field, value), shown in zip(self.fields, self.read):
            if field in ADDRESSES:
                with self.subTest(f"{key} {field}"):
                    self.assertEqual(shown, value)

    def test_the_gott_li_spacetimes_thirteenth_reference_is_elliptic_de_sitter_space(self):
        metric = json.loads((build.METRICS_DIR / "self_creating_universe.json").read_text(encoding="utf-8"))
        self.assertEqual(metric["references"][12], "parikh2003")
        self.assertEqual(self.entries["parikh2003"]["fields"]["title"], "Elliptic de Sitter space: $dS/\\mathbb{Z}_2$")


class InTheBrowser(accents.BuiltSite, unittest.TestCase):

    def spec(self):
        text = build.BIB_FILE.read_text(encoding="utf-8")
        entries = build.parse_bibtex(text)
        listed = {e["id"] for e in json.loads(build.INDEX_FILE.read_text(encoding="utf-8"))}
        metrics = {m["id"]: m["references"] for m in build.load_metrics()
                   if m["id"] in listed and m.get("references")}

        def references(spacetime, keys):
            order = metrics[spacetime]
            return {k: {"number": order.index(k) + 1, **expected(entries[k]["fields"])} for k in sorted(keys)}

        # The Gott-Li spacetime whole, then each spacetime in turn that shows the most of the
        # references carrying TeX not yet seen, until every one cited is.
        spacetimes = {"self_creating_universe": references("self_creating_universe", metrics["self_creating_universe"])}
        citing = {i: {k for k in keys if carries_tex(entries[k]["fields"])} for i, keys in metrics.items()}
        left = set().union(*citing.values()) - set(metrics["self_creating_universe"])
        while left:
            best = max(sorted(citing), key=lambda i: len(citing[i] & left))
            spacetimes[best] = references(best, citing[best] & left)
            left -= citing[best]

        keys = json.loads((build.ROOT / "_data" / "publication_keys.json").read_text(encoding="utf-8"))
        publications = [{"key": k, **expected(entries[k]["fields"])} for k in keys]
        # The publications' own bibliography, the first carrying the Gott-Li spacetime's 13th title.
        texts = entry_texts(text)
        title = entries["parikh2003"]["fields"]["title"]
        first = re.sub(r"(\btitle\s*=\s*)\{(?:[^{}]|\{[^{}]*\})*\}", lambda m: m.group(1) + "{" + title + "}",
                       texts[keys[0]], count=1)
        with_mathematics = [{"key": keys[0], **expected({**entries[keys[0]]["fields"], "title": title})}] + publications[1:]
        return {
            "spacetimes": spacetimes,
            "applets": [{"key": k, **expected(e["fields"], ("author", "title", *VENUES, "volume", "pages"))}
                        for k, e in entries.items()],
            "publications": publications,
            "withMathematics": {"bib": "\n\n".join([first] + [texts[k] for k in keys[1:]]),
                                "publications": with_mathematics},
        }

    def test_every_reference_shows_its_mathematics_set_and_no_tex(self):
        spec = self.spec()
        self.assertEqual(spec["spacetimes"]["self_creating_universe"]["parikh2003"],
                         {"number": 13, "formulas": 1, "words": ["Elliptic de Sitter space:"]})
        path = Path(self.home.name) / "references.json"
        path.write_text(json.dumps(spec, ensure_ascii=False), encoding="utf-8")
        run = subprocess.run(["node", str(build.ROOT / "_tools" / "reference_tex.mjs"), self.base, str(path)],
                             capture_output=True, text=True, timeout=1200)
        self.assertEqual(run.returncode, 0, run.stdout + run.stderr)


if __name__ == "__main__":
    unittest.main()
