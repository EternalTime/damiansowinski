#!/usr/bin/env python3
"""Generate the My Favorite Spacetimes index and bibliography from the data folder.

Run from anywhere:

    python3 _tools/build_mfs_data.py

Reads   MFS/assets/data/metrics/*.json and assets/data/references.bib
Writes  MFS/assets/data/metrics_index.json and MFS/assets/data/references.json

Pass --check to verify the written files are up to date without changing them.
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

VERSION_LENGTH = 16


class DataError(Exception):
    pass


def content_version(value):
    """A stamp over the meaning of `value`, blind to key order and formatting."""
    canonical = json.dumps(value, sort_keys=True, ensure_ascii=False, separators=(",", ":"))
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()[:VERSION_LENGTH]


def sort_key(metric):
    """The order the search list is shown in: by name, accents folded away.

    `sort_name` exists for the spacetimes whose place in the list is not where
    their displayed name would put them.
    """
    name = metric.get("sort_name") or metric["short_name"]
    folded = "".join(c for c in unicodedata.normalize("NFKD", name) if not unicodedata.combining(c))
    return folded.casefold()


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
        metrics.append(metric)

    ids = [m["id"] for m in metrics]
    if len(set(ids)) != len(ids):
        raise DataError("two metric files share an id")
    return sorted(metrics, key=sort_key)


def build_index(metrics):
    return [
        {
            "id": m["id"],
            "name": m["short_name"],
            "tags": m["tags"],
            "version": content_version(m),
        }
        for m in metrics
    ]


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
        outputs = {
            INDEX_FILE: serialise(build_index(metrics)),
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
        print("index and bibliography are up to date")
    return 0


if __name__ == "__main__":
    sys.exit(main())
