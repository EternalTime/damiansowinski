#!/usr/bin/env python3
"""Generate the My Favorite Spacetimes index and bibliography from the data folder.

Run from anywhere:

    python3 _tools/build_mfs_data.py

Reads   MFS/assets/data/metrics/*.json, MFS/assets/data/diagrams/*.json and assets/data/references.bib
Writes  MFS/assets/data/metrics_index.json and MFS/assets/data/references.json

Pass --check to verify the written files are up to date without changing them.

The diagram files are drawn by _tools/derivations/null_rays.py, which needs sympy and
numpy; this script only reads them. It refuses a diagram drawn from components its metric
no longer publishes, and stamps each diagram file's version into the index beside the
metric's own.
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
        for system_id, views in data.get("systems", {}).items():
            if system_id not in systems:
                raise DataError(f"diagrams/{path.name} draws {system_id!r}, "
                                f"which {path.stem}.json has no coordinate system for")
            for view in views:
                source = view["source"]
                if diagram_source_version(systems[system_id], source["fields"]) != source["version"]:
                    raise DataError(
                        f"diagrams/{path.name}: the {system_id} view {view['id']!r} was drawn from "
                        f"components {path.stem}.json no longer publishes; redraw it with "
                        f"_tools/derivations/null_rays.py --metric {path.stem}")
        diagrams[path.stem] = data
    return diagrams


def build_index(metrics, diagrams=None):
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
        diagrams = load_diagrams(metrics)
        outputs = {
            INDEX_FILE: serialise(build_index(metrics, diagrams)),
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
