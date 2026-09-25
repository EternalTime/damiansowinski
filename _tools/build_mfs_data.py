#!/usr/bin/env python3
"""Generate the My Favorite Spacetimes index and bibliography from the data folder.

Run from anywhere:

    python3 _tools/build_mfs_data.py

Reads   MFS/assets/data/metrics/*.json, MFS/assets/data/diagrams/*.json,
        MFS/assets/data/conformal/*.json and assets/data/references.bib
Writes  MFS/assets/data/metrics_index.json and MFS/assets/data/references.json

Pass --check to verify the written files are up to date without changing them.

The diagram files are drawn by _tools/derivations/null_rays.py and the conformal diagram
files by _tools/derivations/conformal.py, which need sympy and numpy; this script only
reads them. It refuses either kind drawn from components its metric no longer publishes,
and stamps each file's version into the index beside the metric's own.
"""

import argparse
import hashlib
import json
import re
import sys
import unicodedata
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
METRICS_DIR = ROOT / "MFS" / "assets" / "data" / "metrics"
INDEX_FILE = ROOT / "MFS" / "assets" / "data" / "metrics_index.json"
BIB_FILE = ROOT / "assets" / "data" / "references.bib"
REFERENCES_FILE = ROOT / "MFS" / "assets" / "data" / "references.json"
DIAGRAMS_DIR = ROOT / "MFS" / "assets" / "data" / "diagrams"
CONFORMAL_DIR = ROOT / "MFS" / "assets" / "data" / "conformal"

VERSION_LENGTH = 16


class DataError(Exception):
    pass


def content_version(value):
    """A stamp over the meaning of `value`, blind to key order and formatting."""
    canonical = json.dumps(value, sort_keys=True, ensure_ascii=False, separators=(",", ":"))
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()[:VERSION_LENGTH]


def name_key(name):
    """The key a displayed name sorts by: its letters as written, case and accents ignored.

    Unicode's compatibility decomposition (NFKD) splits a letter such as ö into o and a
    combining diaeresis, and dropping the combining marks leaves the base letter, so Gödel
    sorts as Godel and Natário as Natario without a table of letters kept here. casefold
    then puts de Sitter under D beside Anti-de Sitter under A, and pp-wave under P.
    """
    folded = "".join(c for c in unicodedata.normalize("NFKD", name) if not unicodedata.combining(c))
    return folded.casefold()


def sort_key(metric):
    """The order the search list is shown in: by the name the reader sees there.

    That name is `short_name`, the one the index publishes as `name`. There is no per entry
    override: a hand kept position is wrong the moment another spacetime is added. The id
    only breaks a tie between two equal names, so the order never depends on the folder.
    """
    return name_key(metric["short_name"]), metric["id"]


CONFLICT_COPY = re.compile(r" \d+$")


def load_metrics():
    metrics = []
    for path in sorted(METRICS_DIR.glob("*.json")):
        if CONFLICT_COPY.search(path.stem):
            continue
        try:
            metric = json.loads(path.read_text(encoding="utf-8"))
        except json.JSONDecodeError as exc:
            raise DataError(f"{path.name} is not valid JSON: {exc}") from exc
        for field in ("id", "name", "short_name", "tags"):
            if not metric.get(field):
                raise DataError(f"{path.name} has no {field}")
        if metric["id"] != path.stem:
            raise DataError(f"{path.name} carries the id {metric['id']!r}")
        if "sort_name" in metric:
            raise DataError(f"{path.name} has a sort_name; the list is ordered by short_name alone")
        metrics.append(metric)

    ids = [m["id"] for m in metrics]
    if len(set(ids)) != len(ids):
        raise DataError("two metric files share an id")
    return sorted(metrics, key=sort_key)


def diagram_source(system, fields):
    """The published fields of a coordinate system that a diagram was drawn from.

    Parameters count by their symbols alone, so rewording a parameter's description
    leaves every diagram standing.
    """
    source = {}
    for name in fields:
        if name == "parameters":
            source[name] = [p["symbol"] for p in system.get("parameters", [])]
        else:
            source[name] = system.get(name)
    return source


def diagram_source_version(system, fields):
    """The stamp a diagram view records of what it was drawn from, and is checked against."""
    return content_version(diagram_source(system, fields))


def check_diagram_stamp(path, system, what, source):
    if diagram_source_version(system, source["fields"]) != source["version"]:
        raise DataError(
            f"diagrams/{path.name}: {what} components {path.stem}.json no longer publishes; redraw it "
            f"with _tools/derivations/null_rays.py --metric {path.stem}")


def load_diagrams(metrics):
    """Every diagram file, refused if it was drawn from what its metric no longer publishes."""
    by_id = {m["id"]: m for m in metrics}
    diagrams = {}
    for path in sorted(DIAGRAMS_DIR.glob("*.json")):
        if CONFLICT_COPY.search(path.stem):
            continue
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except json.JSONDecodeError as exc:
            raise DataError(f"diagrams/{path.name} is not valid JSON: {exc}") from exc
        if data.get("metric") != path.stem:
            raise DataError(f"diagrams/{path.name} carries the metric {data.get('metric')!r}")
        if path.stem not in by_id:
            raise DataError(f"diagrams/{path.name} has no metric file to belong to")
        systems = {s["id"]: s for s in by_id[path.stem].get("coordinates") or []}
        # The flat views, the figures projected from three coordinates, and the reason a
        # system has neither, each stamped with what it was drawn from or speaks of.
        for part, what in (("systems", "view"), ("projections", "figure")):
            for system_id, views in data.get(part, {}).items():
                if system_id not in systems:
                    raise DataError(f"diagrams/{path.name} draws {system_id!r}, "
                                    f"which {path.stem}.json has no coordinate system for")
                for view in views:
                    check_diagram_stamp(path, systems[system_id], f"the {system_id} {what} {view['id']!r} was "
                                        "drawn from", view["source"])
        for system_id, reason in data.get("none", {}).items():
            if system_id not in systems:
                raise DataError(f"diagrams/{path.name} gives a reason for {system_id!r}, "
                                f"which {path.stem}.json has no coordinate system for")
            if system_id in data.get("systems", {}) or system_id in data.get("projections", {}):
                raise DataError(f"diagrams/{path.name} gives a reason for drawing nothing of {system_id!r}, "
                                "and draws it")
            check_diagram_stamp(path, systems[system_id], f"the reason {system_id} is not drawn speaks of",
                                reason["source"])
        diagrams[path.stem] = data
    return diagrams


def load_conformal(metrics):
    """Every conformal diagram file, refused if any field it was drawn from has changed.

    A conformal diagram is of a whole spacetime and may read more than one coordinate
    system, or another spacetime's, as the interior Schwarzschild star reads the exterior
    of schwarzschild.json, so the file lists every system it read under `source`.
    """
    by_id = {m["id"]: m for m in metrics}
    conformal = {}
    for path in sorted(CONFORMAL_DIR.glob("*.json")):
        if CONFLICT_COPY.search(path.stem):
            continue
        where = f"conformal/{path.name}"
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except json.JSONDecodeError as exc:
            raise DataError(f"{where} is not valid JSON: {exc}") from exc
        if data.get("metric") != path.stem:
            raise DataError(f"{where} carries the metric {data.get('metric')!r}")
        if path.stem not in by_id:
            raise DataError(f"{where} has no metric file to belong to")
        if not data.get("source"):
            raise DataError(f"{where} does not say what it was drawn from")
        for source in data["source"]:
            metric = by_id.get(source["metric"])
            if metric is None:
                raise DataError(f"{where} was drawn from {source['metric']}.json, which does not exist")
            systems = {s["id"]: s for s in metric.get("coordinates") or []}
            if source["system"] not in systems:
                raise DataError(f"{where} was drawn from {source['system']!r}, which "
                                f"{source['metric']}.json has no coordinate system for")
            if diagram_source_version(systems[source["system"]], source["fields"]) != source["version"]:
                raise DataError(
                    f"{where} was drawn from components of the {source['system']} system that "
                    f"{source['metric']}.json no longer publishes; redraw it with "
                    f"_tools/derivations/conformal.py --metric {path.stem}")
        for view in data.get("views", []):
            if view.get("system") and view["system"] not in {s["id"] for s in by_id[path.stem]["coordinates"]}:
                raise DataError(f"{where}: the view {view['id']!r} tints {view['system']!r}, which "
                                f"{path.stem}.json has no coordinate system for")
        conformal[path.stem] = data
    return conformal


def build_index(metrics, diagrams=None, conformal=None):
    index = []
    for m in metrics:
        entry = {
            "id": m["id"],
            "name": m["short_name"],
            "tags": m["tags"],
            "version": content_version(m),
        }
        if diagrams and m["id"] in diagrams:
            entry["diagrams"] = content_version(diagrams[m["id"]])
        if conformal and m["id"] in conformal:
            entry["conformal"] = content_version(conformal[m["id"]])
        index.append(entry)
    return index


ENTRY_START = re.compile(r"@(\w+)\s*\{\s*([^,\s]+)\s*,", re.MULTILINE)
FIELD = re.compile(r"(\w+)\s*=\s*")


def read_braced(text, start):
    """Return the text inside the braces that open at `start`, and the index after them."""
    depth = 0
    for i in range(start, len(text)):
        if text[i] == "{":
            depth += 1
        elif text[i] == "}":
            depth -= 1
            if depth == 0:
                return text[start + 1:i], i + 1
    raise DataError("unclosed brace in the bibliography")


def parse_bibtex(text):
    entries = {}
    for match in ENTRY_START.finditer(text):
        entry_type, key = match.group(1).lower(), match.group(2)
        if key in entries:
            raise DataError(f"the bibliography has two entries keyed {key!r}")
        body, _ = read_braced(text, text.index("{", match.start()))
        fields = {}
        position = body.index(",") + 1
        while True:
            field = FIELD.search(body, position)
            if not field:
                break
            after = field.end()
            if after >= len(body):
                raise DataError(f"the entry {key!r} ends with a field that has no value")
            if body[after] == "{":
                value, position = read_braced(body, after)
            elif body[after] == '"':
                end = body.index('"', after + 1)
                value, position = body[after + 1:end], end + 1
            else:
                end = body.find(",", after)
                end = len(body) if end == -1 else end
                value, position = body[after:end], end
            fields[field.group(1).lower()] = " ".join(value.split())
        entries[key] = {"type": entry_type, "fields": fields}
    return entries


def build_references():
    entries = parse_bibtex(BIB_FILE.read_text(encoding="utf-8"))
    if not entries:
        raise DataError("no entries were read out of the bibliography")
    return {"version": content_version(entries), "entries": entries}


def check_citations(metrics, entries):
    """Refuse to publish a citation the reader cannot resolve in the bibliography."""
    for metric in metrics:
        for key in metric.get("references", []):
            if key not in entries:
                raise DataError(
                    f"{metric['id']}.json cites {key!r}, which the bibliography has no entry for"
                )


# A history reads evenly: at least five paragraphs, each of three to six sentences, and the
# longest paragraph no more than twice the length of the shortest. A table is not a
# paragraph of prose and is left out of the count.
HISTORY_PARAGRAPHS = 5
HISTORY_SENTENCES = (3, 6)
HISTORY_SPREAD = 2

MATH = re.compile(r"\$\$.+?\$\$|\$[^$]+\$", re.DOTALL)
STOP = re.compile(r"[.?!][\"')\]]*\s+")


def sentences(paragraph):
    """The sentences of one paragraph of prose.

    A sentence ends at a full stop, a question mark or an exclamation mark, after any
    closing quote or bracket, where the next word begins with a capital letter or a digit.
    A capital standing alone before a full stop is an initial, as in J. Robert Oppenheimer,
    and ends nothing. Mathematics counts as a single word, which ends a sentence only when
    it closes with the stop itself, as a displayed equation at the end of a sentence does.
    """
    text = MATH.sub(lambda m: "MATH." if re.search(r"[.?!]\s*\$+$", m.group(0)) else "MATH",
                    paragraph).strip()
    found, start = [], 0
    for stop in STOP.finditer(text):
        after = text[stop.end():].lstrip("\"'([")
        initial = (text[stop.start()] == "." and stop.start() > 0 and text[stop.start() - 1].isupper()
                   and (stop.start() == 1 or not text[stop.start() - 2].isalnum()))
        if after and (after[0].isupper() or after[0].isdigit()) and not initial:
            found.append(text[start:stop.end()].strip())
            start = stop.end()
    if text[start:].strip():
        found.append(text[start:].strip())
    return found


def prose_paragraphs(history):
    """Each paragraph of a history that is prose, with its place among all of them."""
    return [(n, p) for n, p in enumerate(history.split("¶"), 1) if not p.startswith("TABLE::")]


def history_shape(history):
    """The number of sentences in each paragraph of a history, tables left out."""
    return [len(sentences(p)) for _, p in prose_paragraphs(history)]


def check_history_shape(metrics):
    """Refuse a history too short to read as a story, or with paragraphs of uneven length.

    Every history out of shape is named at once, with each paragraph that breaks the rule.
    """
    low, high = HISTORY_SENTENCES
    problems = []
    for metric in metrics:
        paragraphs = prose_paragraphs(metric.get("history") or "")
        shape = [len(sentences(p)) for _, p in paragraphs]
        where = f"{metric['id']}.json: the history, {shape},"
        if len(shape) < HISTORY_PARAGRAPHS:
            problems.append(f"{where} has {len(shape)} paragraphs and needs at least "
                            f"{HISTORY_PARAGRAPHS}")
        for (position, _), count in zip(paragraphs, shape):
            if not low <= count <= high:
                problems.append(f"{where} has {count} sentences in paragraph {position}, "
                                f"where a paragraph takes {low} to {high}")
        if shape and max(shape) > HISTORY_SPREAD * min(shape):
            problems.append(f"{where} has a longest paragraph more than {HISTORY_SPREAD} "
                            f"times its shortest")
    if problems:
        raise DataError("\n".join(problems))


def serialise(value):
    return json.dumps(value, indent=2, ensure_ascii=False) + "\n"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--check", action="store_true", help="report drift instead of writing")
    args = parser.parse_args(argv)

    try:
        references = build_references()
        metrics = load_metrics()
        check_citations(metrics, references["entries"])
        check_history_shape(metrics)
        diagrams = load_diagrams(metrics)
        conformal = load_conformal(metrics)
        outputs = {
            INDEX_FILE: serialise(build_index(metrics, diagrams, conformal)),
            REFERENCES_FILE: serialise(references),
        }
    except DataError as exc:
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
            path.write_text(text, encoding="utf-8")
            print(f"wrote {path.relative_to(ROOT)}")
        else:
            print(f"unchanged {path.relative_to(ROOT)}")

    if stale:
        for path in stale:
            print(f"out of date: {path.relative_to(ROOT)}", file=sys.stderr)
        print("run python3 _tools/build_mfs_data.py", file=sys.stderr)
        return 1
    if args.check:
        print("index, bibliography and diagram stamps are up to date")
    return 0


if __name__ == "__main__":
    sys.exit(main())
