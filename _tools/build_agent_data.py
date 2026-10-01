#!/usr/bin/env python3
"""Generate the machine-readable layer that AI agents read, from the site's own sources.

Run from anywhere:

    python3 _tools/build_agent_data.py

Reads   MFS/assets/data/metrics/*.json, assets/data/references.bib,
        _data/publication_keys.json, _data/libraries.json, _data/person.json
        and the page sources listed in PAGES
Writes  _data/generated/spacetimes.json, _data/generated/publications.json
        and llms-full.txt

Pass --check to verify the written files are up to date without changing them.

Jekyll publishes the two data files at /data/spacetimes.json and /data/publications.json
and builds /llms.txt from them, so everything an agent reads comes from the same sources
the pages read, and none of it is copied by hand.
"""

import argparse
import html
import json
import re
import sys
from pathlib import Path
from urllib.parse import urlsplit

import build_mfs_data as mfs

ROOT = mfs.ROOT
SITE_URL = "https://damiansowinski.com"
DATA_DIR = ROOT / "_data"
GENERATED_DIR = DATA_DIR / "generated"
SPACETIMES_FILE = GENERATED_DIR / "spacetimes.json"
PUBLICATIONS_FILE = GENERATED_DIR / "publications.json"
FULL_TEXT_FILE = ROOT / "llms-full.txt"

# The pages llms-full.txt carries, in the order it carries them. A page whose content is
# drawn by a script in the browser also names the section built here from that script's data.
PAGES = [
    ("index.markdown", "Home", "/", None),
    ("about.markdown", "About", "/about/", None),
    ("research.markdown", "Research", "/research/", None),
    ("publications.markdown", "Publications", "/publications/", "publications"),
    ("MFS/index.markdown", "My Favorite Spacetimes", "/MFS/", "spacetimes"),
    ("code.markdown", "Code", "/code.html", "libraries"),
    ("applets.markdown", "Applets", "/applets.html", None),
    ("notes.markdown", "Notes", "/notes/", None),
    ("teaching.markdown", "Teaching", "/teaching/", None),
    ("teaching/ice-dartmouth/Lectures.markdown", "From Epistemology to Information (lecture series)",
     "/teaching/ice-dartmouth/Lectures/", None),
    ("teaching/university-of-rochester/PHYSLABS.markdown", "Physics Lab Documentation",
     "/teaching/university-of-rochester/PHYSLABS/", None),
    ("teaching/university-of-rochester/PHYS100.markdown", "PHYS 100: Physics and the Natural World",
     "/teaching/university-of-rochester/PHYS100/", None),
    ("teaching/university-of-rochester/PHYS141.markdown", "PHYS 141: Honors Classical Mechanics",
     "/teaching/university-of-rochester/PHYS141/", None),
    ("press.markdown", "Press", "/press/", None),
]


def absolute(url):
    """The address `url` names on the live site, whatever form the page wrote it in."""
    url = url.strip()
    if not url or urlsplit(url).scheme:
        return url
    return SITE_URL + "/" + url.lstrip("/")


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


# ---------------------------------------------------------------------------------------
# Spacetimes


def spacetime_entry(metric):
    charts = []
    for chart in metric["coordinates"]:
        charts.append({
            "id": chart["id"],
            "name": chart["name"],
            "coordinates": chart["coords"],
            "domains": chart["domains"],
            "parameters": chart["parameters"],
            "convention": chart.get("convention"),
            "line_element": chart["line_element"],
            "metric_components": chart["metric_components"],
            "ricci_scalar": chart["ricci_scalar"],
            "kretschmann": chart["kretschmann"],
        })
    return {
        "id": metric["id"],
        "name": metric["name"],
        "short_name": metric["short_name"],
        "description": metric["description"],
        "tags": metric["tags"],
        "signature": metric.get("signature"),
        "convention": metric.get("convention"),
        "related": metric.get("related", []),
        "references": metric["references"],
        "page_url": absolute("/MFS/"),
        "data_url": absolute(f"/MFS/assets/data/metrics/{metric['id']}.json"),
        "charts": charts,
    }


def build_spacetimes(metrics):
    return {
        "name": "My Favorite Spacetimes",
        "description": (
            "A catalogue of exact solutions of Einstein's field equations. For each spacetime: "
            "its name, a description, tags, signature, the conventions shared by its charts and "
            "the spacetimes it is related to, each with its id, the kind of relation and why, "
            "and for each coordinate chart the coordinates, their domains, the parameters, the "
            "chart's own conventions, the line element, the nonzero metric components and the "
            "curvature invariants. Mathematics is LaTeX."
        ),
        "page_url": absolute("/MFS/"),
        "references_url": absolute("/MFS/assets/data/references.json"),
        "note": (
            "The catalogue page chooses a spacetime in the browser and has no address per "
            "spacetime, so page_url is the catalogue itself. data_url is the full record the "
            "page draws from: Christoffel symbols, Riemann, Ricci, Einstein and Weyl tensors, "
            "geodesic equations and the history with its references."
        ),
        "count": len(metrics),
        "spacetimes": [spacetime_entry(m) for m in sorted(metrics, key=mfs.sort_key)],
    }


# ---------------------------------------------------------------------------------------
# Publications

# The same accents the publications page's own reader turns into characters.
ACCENTS = {
    '"': dict(zip("aeiouAEIOUy", "äëïöüÄËÏÖÜÿ")),
    "'": dict(zip("aeiouAEIOUy", "áéíóúÁÉÍÓÚý")),
    "`": dict(zip("aeiouAEIOU", "àèìòùÀÈÌÒÙ")),
    "^": dict(zip("aeiouAEIOU", "âêîôûÂÊÎÔÛ")),
    "~": dict(zip("nNaAoO", "ñÑãÃõÕ")),
    "c": dict(zip("cCsS", "çÇşŞ")),
}
ACCENT = re.compile(r"""\\(["'`^~])\{?([A-Za-z])\}?|\\c\{?([cCsS])\}?""")
ARXIV = re.compile(r"arXiv[:\s.]*(\d{4}\.\d{4,5})", re.IGNORECASE)


def unlatex(text):
    def accent(match):
        mark, letter = (match.group(1), match.group(2)) if match.group(1) else ("c", match.group(3))
        return ACCENTS[mark].get(letter, letter)

    text = ACCENT.sub(accent, text)
    text = re.sub(r"\\ss\b", "ß", text).replace("\\&", "&")
    text = re.sub(r"\{([^{}]*)\}", r"\1", text)
    return " ".join(text.split())


def author_name(author):
    if author == "others":
        return "et al."
    if "," in author:
        last, first = (part.strip() for part in author.split(",", 1))
        return f"{first} {last}"
    return author


def publication_entry(key, entry):
    fields = {name: unlatex(value) for name, value in entry["fields"].items()}
    authors = [author_name(a.strip()) for a in re.split(r"\s+and\s+", fields.get("author", "")) if a.strip()]
    venue = fields.get("journal") or fields.get("booktitle") or fields.get("publisher")
    doi = fields.get("doi")
    arxiv_match = ARXIV.search(" ".join(filter(None, [venue, doi, fields.get("eprint")])))
    arxiv = arxiv_match.group(1) if arxiv_match else None
    pdf = absolute(fields["url"]) if fields.get("url") else None
    links = {
        "doi_url": f"https://doi.org/{doi}" if doi else None,
        "arxiv_url": f"https://arxiv.org/abs/{arxiv}" if arxiv else None,
        "pdf_url": pdf,
    }
    return {
        "key": key,
        "type": entry["type"],
        "title": fields.get("title"),
        "authors": authors,
        "year": int(fields["year"]) if fields.get("year", "").isdigit() else fields.get("year"),
        "venue": venue,
        "volume": fields.get("volume"),
        "number": fields.get("number"),
        "pages": fields["pages"].replace("--", "\u2013") if fields.get("pages") else None,
        "publisher": fields.get("publisher") if fields.get("publisher") != venue else None,
        "doi": doi,
        "arxiv": arxiv,
        **links,
        "url": links["doi_url"] or links["arxiv_url"] or pdf,
    }


def build_publications(entries, keys):
    missing = [key for key in keys if key not in entries]
    if missing:
        raise mfs.DataError(
            "_data/publication_keys.json names keys that are not in the bibliography: "
            + ", ".join(missing))
    publications = [publication_entry(key, entries[key]) for key in keys]
    return {
        "name": "Publications of Damian R. Sowinski",
        "page_url": absolute("/publications/"),
        "note": "In the order the publications page lists them. pdf_url is a copy with no pay-wall.",
        "count": len(publications),
        "publications": [{k: v for k, v in p.items() if v is not None} for p in publications],
    }


# ---------------------------------------------------------------------------------------
# llms-full.txt

LIQUID_URL = re.compile(r"""\{\{\s*["']([^"']*)["']\s*(\|\s*relative_url\s*)?\}\}""")
LIQUID_OTHER = re.compile(r"\{\{.*?\}\}|\{%.*?%\}", re.DOTALL)
DROPPED = re.compile(
    r"<(script|style|video|canvas|noscript)\b.*?</\1\s*>|<!--.*?-->", re.DOTALL | re.IGNORECASE)
# An icon sits inside a sentence, so the sentence closes over it.
ICON = re.compile(r"\s*<svg\b.*?</svg\s*>\s*", re.DOTALL | re.IGNORECASE)
# Blocks a script fills in the browser: their markup is an empty frame, so it is dropped
# and the page's section below is built from the script's data instead.
WIDGETS = ("repo-terminal",)
DIV = re.compile(r"<(/?)div\b[^>]*>", re.IGNORECASE)
IFRAME = re.compile(r"<iframe\b[^>]*\bsrc=\"([^\"]*)\"[^>]*>.*?</iframe>", re.DOTALL | re.IGNORECASE)
YOUTUBE = re.compile(r"youtube\.com/embed/([\w-]+)")
LINK = re.compile(r"<a\b([^>]*)>(.*?)</a\s*>", re.DOTALL | re.IGNORECASE)
HREF = re.compile(r"""\bhref\s*=\s*["']([^"']*)["']""", re.IGNORECASE)
HEADING = re.compile(r"<h([1-6])\b[^>]*>(.*?)</h\1\s*>", re.DOTALL | re.IGNORECASE)
WIP = re.compile(r"<span class=\"wip-label\">(.*?)</span>", re.DOTALL)
BUTTON = re.compile(r"<button\b[^>]*>(.*?)</button\s*>", re.DOTALL | re.IGNORECASE)
TAG = re.compile(r"<[^>]+>")


def inline_text(fragment):
    return " ".join(html.unescape(TAG.sub("", fragment)).split())


def drop_widgets(text):
    for name in WIDGETS:
        start = re.search(rf"<div\b[^>]*\bclass=\"[^\"]*\b{name}\b[^\"]*\"[^>]*>", text)
        if not start:
            raise mfs.DataError(f"no <div class=\"{name}\"> is left to drop; update WIDGETS")
        depth = 0
        for tag in DIV.finditer(text, start.start()):
            depth += -1 if tag.group(1) else 1
            if depth == 0:
                text = text[:start.start()] + text[tag.end():]
                break
    return text


def page_markdown(source):
    """The prose of a page source as Markdown, with every address absolute."""
    text = re.sub(r"\A---\n.*?\n---\n", "", source, count=1, flags=re.DOTALL)
    text = LIQUID_URL.sub(lambda m: absolute(m.group(1)) if m.group(2) else m.group(1), text)
    text = LIQUID_OTHER.sub("", text)
    text = DROPPED.sub("", text)
    text = ICON.sub(" ", text)
    if "repo-terminal" in text:
        text = drop_widgets(text)

    def iframe(match):
        video = YOUTUBE.search(match.group(1))
        if video:
            return f"\n[Video on YouTube](https://www.youtube.com/watch?v={video.group(1)})\n"
        return f"\n[Embedded player: {urlsplit(match.group(1)).hostname}]({match.group(1)})\n"

    def link(match):
        label = inline_text(match.group(2))
        href = HREF.search(match.group(1))
        target = absolute(href.group(1)) if href else ""
        if label.startswith("\u2190"):  # the "← Back to ..." links a reader clicks, not reads
            return ""
        if not label:
            return ""
        text = f"[{label}]({target})" if target else label
        # A space the page kept inside the link still divides it from the words beside it.
        raw = html.unescape(TAG.sub("", match.group(2)))
        return (" " if raw[:1].isspace() else "") + text + (" " if raw[-1:].isspace() else "")

    # The page's own heading levels, however it spelled them, nest under its "## " section.
    levels = sorted({int(level) for level, _ in HEADING.findall(text)})

    def heading(match):
        level = min(6, 3 + levels.index(int(match.group(1))))
        return f"\n\n{'#' * level} {inline_text(match.group(2))}\n\n"

    text = IFRAME.sub(iframe, text)
    text = LINK.sub(link, text)
    text = HEADING.sub(heading, text)
    text = WIP.sub(lambda m: f" ({inline_text(m.group(1)).lower()})", text)
    text = BUTTON.sub(lambda m: f"\n- {inline_text(m.group(1))} (interactive applet)", text)
    text = re.sub(r"<li\b[^>]*>", "\n- ", text, flags=re.IGNORECASE)
    text = re.sub(r"<br\s*/?>[ \t]*\n?|</(p|div|ul|ol|li)>", "\n", text, flags=re.IGNORECASE)
    text = re.sub(r"</?(strong|b)>", "**", text, flags=re.IGNORECASE)
    text = re.sub(r"</?(em|i)>", "*", text, flags=re.IGNORECASE)
    text = html.unescape(TAG.sub("", text))

    lines = []
    for line in text.splitlines():
        line = " ".join(line.split())
        if re.fullmatch(r"-{3,}|\*\*|\*", line):
            line = ""
        lines.append(line)
    text = "\n".join(lines)
    text = re.sub(r"\n{3,}", "\n\n", text).strip()
    # A list the page wrote as separate lines reads as one list.
    return re.sub(r"(?m)^(- .*)\n\n(?=- )", r"\1\n", text)


def publications_markdown(publications):
    lines = []
    for p in publications["publications"]:
        where = ", ".join(filter(None, [p.get("venue"), p.get("volume")]))
        line = f"- {', '.join(p['authors'])} ({p['year']}). [{p['title']}]({p['url']}). {where}."
        extra = [f"PDF: {p['pdf_url']}"] if p.get("pdf_url") and p["pdf_url"] != p["url"] else []
        lines.append(line + (" " + " ".join(extra) if extra else ""))
    return "\n".join(lines)


def spacetimes_markdown(spacetimes):
    parts = [
        f"Machine-readable catalogue: {absolute('/data/spacetimes.json')}",
    ]
    for s in spacetimes["spacetimes"]:
        block = [f"### {s['name']}", "", s["description"], "",
                 f"- id: `{s['id']}`", f"- tags: {', '.join(s['tags'])}"]
        if s.get("signature"):
            block.append(f"- signature: {s['signature']}")
        block.append(f"- full data: {s['data_url']}")
        for chart in s["charts"]:
            block.append(f"- chart {chart['name']}: coordinates $({', '.join(chart['coordinates'])})$, "
                         f"$${chart['line_element']}$$")
        parts.append("\n".join(block))
    return "\n\n".join(parts)


def libraries_markdown(libraries):
    lines = []
    for lib in libraries:
        repo = f"https://github.com/EternalTime/{lib.get('gh', lib['repo'])}"
        lines.append(f"- [{lib['name']}: {lib['label']}]({absolute(lib['docs'])}), source {repo}")
    return "\n".join(lines)


def build_full_text(person, spacetimes, publications, libraries):
    generated = {
        "publications": publications_markdown(publications),
        "spacetimes": spacetimes_markdown(spacetimes),
        "libraries": libraries_markdown(libraries),
    }
    out = [
        f"# {person['name']}",
        "",
        f"> {person['summary']}",
        "",
        "The full text of the main pages of damiansowinski.com as Markdown, for language models, "
        "generated from the page sources; the pages themselves are the reference.",
        f"Index: {absolute('/llms.txt')}",
    ]
    for source, title, url, extra in PAGES:
        body = page_markdown((ROOT / source).read_text(encoding="utf-8"))
        section = [f"## {title}", "", f"URL: {absolute(url)}"]
        if body:
            section += ["", body]
        if extra:
            section += ["", generated[extra]]
        out += ["", "\n".join(section)]
    return "\n".join(out).rstrip() + "\n"


# ---------------------------------------------------------------------------------------


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--check", action="store_true", help="report drift instead of writing")
    args = parser.parse_args(argv)

    try:
        entries = mfs.parse_bibtex(mfs.BIB_FILE.read_text(encoding="utf-8"))
        metrics = mfs.load_metrics()
        spacetimes = build_spacetimes(metrics)
        publications = build_publications(entries, read_json(DATA_DIR / "publication_keys.json"))
        full_text = build_full_text(read_json(DATA_DIR / "person.json"), spacetimes, publications,
                                    read_json(DATA_DIR / "libraries.json"))
        outputs = {
            SPACETIMES_FILE: mfs.serialise(spacetimes),
            PUBLICATIONS_FILE: mfs.serialise(publications),
            FULL_TEXT_FILE: full_text,
        }
    except mfs.DataError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    stale = []
    for path, text in outputs.items():
        current = path.read_text(encoding="utf-8") if path.exists() else None
        if args.check:
            if current != text:
                stale.append(path)
            continue
        if current != text:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(text, encoding="utf-8")
            print(f"wrote {path.relative_to(ROOT)}")
        else:
            print(f"unchanged {path.relative_to(ROOT)}")

    if stale:
        for path in stale:
            print(f"out of date: {path.relative_to(ROOT)}", file=sys.stderr)
        print("run python3 _tools/build_agent_data.py", file=sys.stderr)
        return 1
    if args.check:
        print("the agent layer is up to date")
    return 0


if __name__ == "__main__":
    sys.exit(main())
