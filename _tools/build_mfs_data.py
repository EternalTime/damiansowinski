#!/usr/bin/env python3
"""Generate the My Favorite Spacetimes index and bibliography from the data folder.

Run from anywhere:

    python3 _tools/build_mfs_data.py

Reads   MFS/assets/data/metrics/*.json, MFS/assets/data/diagrams/*.json,
        MFS/assets/data/conformal/*.json, MFS/assets/data/embedding/*.json
        and assets/data/references.bib
Writes  MFS/assets/data/metrics_index.json and MFS/assets/data/references.json

Pass --check to verify the written files are up to date without changing them.

The diagram files are drawn by _tools/derivations/null_rays.py, the conformal diagram
files by _tools/derivations/conformal.py and the embedding diagram files by
_tools/derivations/embedding.py, which need sympy and numpy; this script only reads them.
It refuses any of them drawn from components its metric no longer publishes, and stamps
each file's version into the index beside the metric's own.
"""

import argparse
import hashlib
import json
import math
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
EMBEDDING_DIR = ROOT / "MFS" / "assets" / "data" / "embedding"

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


def embedding_moment(view, index):
    """What an embedding surface says of the moment it is cut from: the view, its settings,
    the surface's label and time, and where each piece runs, which is the reach of the moment
    every other diagram draws. A grid is its frame and the ends of its rows and columns."""
    surface = view["surfaces"][index]
    pieces = []
    for piece in surface["pieces"]:
        if "grid" in piece:
            grid = piece["grid"]
            pieces.append([piece["id"], grid["frame"], grid["u"][0], grid["u"][-1], grid["v"][0], grid["v"][-1]])
        else:
            pieces.append([piece["id"], piece["class"], piece["metric"], piece["system"], piece["coordinate"],
                           piece["points"][0][0], piece["points"][-1][0]])
    return {"view": view["id"], "settings": view.get("settings"), "label": surface.get("label"),
            "time": surface.get("time"), "pieces": pieces,
            "curves": [[c["class"], c["points"][0], c["points"][-1]] for c in surface.get("curves", [])
                       if c["class"] == "path"]}


def embedding_moment_version(view, index):
    """The stamp a slice drawn on another diagram records of the embedding surface it marks."""
    return content_version(embedding_moment(view, index))


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
        # The flat views and the figures projected from three coordinates, each stamped with
        # what it was drawn from.
        for part, what in (("systems", "view"), ("projections", "figure")):
            for system_id, views in data.get(part, {}).items():
                if system_id not in systems:
                    raise DataError(f"diagrams/{path.name} draws {system_id!r}, "
                                    f"which {path.stem}.json has no coordinate system for")
                for view in views:
                    check_diagram_stamp(path, systems[system_id], f"the {system_id} {what} {view['id']!r} was "
                                        "drawn from", view["source"])
        if not any(data.get("systems", {}).values()) and not any(data.get("projections", {}).values()):
            raise DataError(f"diagrams/{path.name} draws nothing")
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
        if not data.get("views"):
            raise DataError(f"{where} draws nothing")
        for view in data["views"]:
            if view.get("system") and view["system"] not in {s["id"] for s in by_id[path.stem]["coordinates"]}:
                raise DataError(f"{where}: the view {view['id']!r} tints {view['system']!r}, which "
                                f"{path.stem}.json has no coordinate system for")
        conformal[path.stem] = data
    return conformal


def check_grid(at, piece):
    """A grid piece, a height over a plane, as _tools/README.md defines it: its values of u and v
    each strictly increasing, a polar grid's u never below the axis and its angles within one
    turn, a height at every node, one height where a polar grid meets its axis, and an edge that
    says what lies beyond it, with no profile beside it."""
    grid = piece["grid"]
    if grid.get("frame") not in ("polar", "cartesian"):
        raise DataError(f"{at} is a grid in neither a polar nor a Cartesian frame")
    u, v, z = grid.get("u") or [], grid.get("v") or [], grid.get("z") or []
    for name, values in (("u", u), ("v", v)):
        if len(values) < 2 or not all(b > a for a, b in zip(values, values[1:])):
            raise DataError(f"{at} does not run one way along its {name}")
    if grid["frame"] == "polar" and (u[0] < 0 or v[0] < 0 or v[-1] >= 2 * math.pi or len(v) < 3):
        raise DataError(f"{at} is a polar grid below the axis or beyond one turn")
    if len(z) != len(u) or any(len(row) != len(v) for row in z):
        raise DataError(f"{at} does not give a height at every node of its grid")
    if not all(isinstance(h, (int, float)) and math.isfinite(h) for row in z for h in row):
        raise DataError(f"{at} gives a height that is not a number")
    if grid["frame"] == "polar" and u[0] == 0 and len(set(z[0])) != 1:
        raise DataError(f"{at} meets its axis at more than one height")
    if any(key in piece for key in ("points", "start", "end", "coordinate")):
        raise DataError(f"{at} is a grid and a profile at once")
    if piece.get("edge", {}).get("kind") != "edge":
        raise DataError(f"{at} does not say what lies beyond its edge")


def load_embedding(metrics):
    """Every embedding diagram file, refused if any field it was drawn from has changed.

    Like a conformal diagram it is of a whole spacetime and may read more than one
    coordinate system, or another spacetime's, as the interior Schwarzschild star reads the
    exterior of schwarzschild.json, so the file lists every system it read under `source`.
    Each piece of each surface must run from its start to its end in one direction of its
    coordinate, never below the axis, and each grid piece hold a height at every node of its
    grid, as _tools/README.md defines them.
    """
    by_id = {m["id"]: m for m in metrics}
    embedding = {}
    for path in sorted(EMBEDDING_DIR.glob("*.json")):
        if CONFLICT_COPY.search(path.stem):
            continue
        where = f"embedding/{path.name}"
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
                    f"_tools/derivations/embedding.py --metric {path.stem}")
        # A spacetime with no surface to draw says why under `stops`, and draws no view; one
        # that draws keeps what it does not draw in its views.
        stated = data.get("views") == [] and bool(data.get("stops"))
        if stated and not all(isinstance(text, str) and text.strip() for text in data["stops"]):
            raise DataError(f"{where} states why it draws nothing in something other than sentences")
        if not stated and (not data.get("views") or not all(view.get("surfaces") for view in data["views"])):
            raise DataError(f"{where} draws nothing")
        if data.get("views") and "stops" in data:
            raise DataError(f"{where} draws views and says beside them that it draws nothing")
        read = {(s["metric"], s["system"]) for s in data["source"]}
        for view in data["views"]:
            if view.get("system") and view["system"] not in {s["id"] for s in by_id[path.stem]["coordinates"]}:
                raise DataError(f"{where}: the view {view['id']!r} names {view['system']!r}, which "
                                f"{path.stem}.json has no coordinate system for")
            for surface in view["surfaces"]:
                for piece in surface["pieces"]:
                    at = f"{where}: the view {view['id']!r}, piece {piece['id']!r}"
                    if (piece["metric"], piece["system"]) not in read:
                        raise DataError(f"{at} was drawn from {piece['metric']}/{piece['system']}, "
                                        "which the file's source does not stamp")
                    if "grid" in piece:
                        check_grid(at, piece)
                        continue
                    x = [point[0] for point in piece["points"]]
                    steps = [b - a for a, b in zip(x, x[1:])]
                    if len(x) < 2 or not (all(d > 0 for d in steps) or all(d < 0 for d in steps)):
                        raise DataError(f"{at} does not run one way along its coordinate")
                    if any(point[1] < 0 for point in piece["points"]):
                        raise DataError(f"{at} has a point below the axis, rho < 0")
                ids = {piece["id"] for piece in surface["pieces"]}
                grids = {piece["id"] for piece in surface["pieces"] if "grid" in piece}
                for ring in surface["rings"]:
                    if ring["piece"] not in ids:
                        raise DataError(f"{where}: the view {view['id']!r} marks a circle on "
                                        f"{ring['piece']!r}, which it does not draw")
                    if ring["piece"] in grids:
                        raise DataError(f"{where}: the view {view['id']!r} marks a circle of one height on "
                                        f"the grid {ring['piece']!r}, which has no circles of one height")
                for mark, kind in [(c, "a curve") for c in surface.get("curves", [])] + \
                                  [(d, "a point") for d in surface.get("dots", [])]:
                    if mark["piece"] not in ids:
                        raise DataError(f"{where}: the view {view['id']!r} marks {kind} on "
                                        f"{mark['piece']!r}, which it does not draw")
                for curve in surface.get("curves", []):
                    if len(curve["points"]) < 2 or any(len(point) != 3 for point in curve["points"]):
                        raise DataError(f"{where}: the view {view['id']!r} marks a curve that is not a "
                                        "polyline of points [X, Y, Z]")
        embedding[path.stem] = data
    return embedding


def check_slices(diagrams, conformal, embedding):
    """Every slice a spacetime diagram, a figure or a conformal diagram draws marks a surface
    of its spacetime's embedding diagram, by the view's id and the surface's place, and was
    drawn from that surface as the embedding file stands: a redrawn embedding whose moments
    or reach moved leaves the drawings that mark them unpublishable until they are redrawn."""
    def drawings():
        for metric_id, data in (diagrams or {}).items():
            for part in ("systems", "projections"):
                for system_id, views in data.get(part, {}).items():
                    for view in views:
                        yield (metric_id, f"diagrams/{metric_id}.json, the {system_id} view {view['id']!r}", view,
                               "null_rays.py --slices")
        for metric_id, data in (conformal or {}).items():
            for view in data["views"]:
                yield metric_id, f"conformal/{metric_id}.json, the view {view['id']!r}", view, "conformal.py"
    for metric_id, where, view, script in drawings():
        marks = view.get("slices", [])
        if not marks:
            continue
        views = {v["id"]: v for v in (embedding or {}).get(metric_id, {}).get("views", [])}
        for mark in marks:
            target = views.get(mark.get("view"))
            if target is None or not 0 <= mark.get("surface", -1) < len(target["surfaces"]):
                raise DataError(f"{where} marks the surface {mark.get('surface')!r} of the embedding view "
                                f"{mark.get('view')!r}, which embedding/{metric_id}.json does not draw")
            if embedding_moment_version(target, mark["surface"]) != mark.get("version"):
                raise DataError(f"{where} marks a moment of embedding/{metric_id}.json as it no longer "
                                f"stands; redraw it with _tools/derivations/{script} --metric {metric_id}")
            if not (mark.get("lines") or mark.get("points") or mark.get("fills")):
                raise DataError(f"{where} marks the moment {mark.get('label')!r} and draws nothing of it")


def build_index(metrics, diagrams=None, conformal=None, embedding=None):
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
        if embedding and m["id"] in embedding:
            entry["embedding"] = content_version(embedding[m["id"]])
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


# Prose set in paragraphs reads evenly: every paragraph of a history or a convention has
# three to six sentences, and the longest is no more than twice the length of the shortest.
# A history tells a story, so it runs to at least five paragraphs; a convention may be a
# single paragraph. A table is not a paragraph of prose and is left out of the count.
HISTORY_PARAGRAPHS = 5
PARAGRAPH_SENTENCES = (3, 6)
PARAGRAPH_SPREAD = 2

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


def prose_paragraphs(text):
    """Each paragraph of a history or a convention that is prose, with its place among all of them."""
    return [(n, p) for n, p in enumerate(text.split("¶"), 1) if not p.startswith("TABLE::")]


def paragraph_shape(text):
    """The number of sentences in each paragraph of a history or a convention, tables left out."""
    return [len(sentences(p)) for _, p in prose_paragraphs(text)]


def shape_problems(metric_id, field, text, fewest):
    """Each way one history or convention fails to read evenly, naming the paragraph at fault."""
    low, high = PARAGRAPH_SENTENCES
    paragraphs = prose_paragraphs(text)
    shape = [len(sentences(p)) for _, p in paragraphs]
    where = f"{metric_id}.json: the {field}, {shape},"
    problems = []
    if len(shape) < fewest:
        problems.append(f"{where} has {len(shape)} paragraphs and needs at least {fewest}")
    for (position, _), count in zip(paragraphs, shape):
        if not low <= count <= high:
            problems.append(f"{where} has {count} sentences in paragraph {position}, "
                            f"where a paragraph takes {low} to {high}")
    if shape and max(shape) > PARAGRAPH_SPREAD * min(shape):
        problems.append(f"{where} has a longest paragraph more than {PARAGRAPH_SPREAD} "
                        f"times its shortest")
    return problems


def check_prose_shape(metrics):
    """Refuse a history or a convention with paragraphs of uneven length, or a history too
    short to read as a story.

    Every one out of shape is named at once, with each paragraph that breaks the rule. A
    spacetime may carry no convention, and then there is nothing to count.
    """
    problems = []
    for metric in metrics:
        problems += shape_problems(metric["id"], "history", metric.get("history") or "",
                                   HISTORY_PARAGRAPHS)
        if metric.get("convention"):
            problems += shape_problems(metric["id"], "convention", metric["convention"], 1)
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
        check_prose_shape(metrics)
        diagrams = load_diagrams(metrics)
        conformal = load_conformal(metrics)
        embedding = load_embedding(metrics)
        check_slices(diagrams, conformal, embedding)
        outputs = {
            INDEX_FILE: serialise(build_index(metrics, diagrams, conformal, embedding)),
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
