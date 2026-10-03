#!/usr/bin/env python3
"""Generate the My Favorite Spacetimes index and bibliography from the data folder.

Run from anywhere:

    python3 _tools/build_mfs_data.py

Reads   MFS/assets/data/metrics/*.json, MFS/assets/data/diagrams/*.json,
        MFS/assets/data/conformal/*.json, MFS/assets/data/embedding/*.json
        and assets/data/references.bib
Writes  MFS/assets/data/metrics_index.json, MFS/assets/data/references.json and
        MFS/assets/data/relations.json, the relations the graph behind the list draws

Pass --check to verify the written files are up to date without changing them.

The diagram files are drawn by _tools/derivations/null_rays.py, the conformal diagram
files by _tools/derivations/conformal.py and the embedding diagram files by
_tools/derivations/embedding.py, which need sympy and numpy; this script only reads them.
It refuses any of them drawn from components its metric no longer publishes, and stamps
each file's version into the index beside the metric's own.

It also refuses a symbol that a chart's mathematics or a drawing uses and nothing defines;
symbols.py, beside this file, holds that rule.
"""

import argparse
import hashlib
import json
import math
import re
import sys
import unicodedata
from pathlib import Path

from symbols import symbol_problems

ROOT = Path(__file__).resolve().parent.parent
METRICS_DIR = ROOT / "MFS" / "assets" / "data" / "metrics"
INDEX_FILE = ROOT / "MFS" / "assets" / "data" / "metrics_index.json"
BIB_FILE = ROOT / "assets" / "data" / "references.bib"
REFERENCES_FILE = ROOT / "MFS" / "assets" / "data" / "references.json"
RELATIONS_FILE = ROOT / "MFS" / "assets" / "data" / "relations.json"
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


def embedding_moment(view, index, curve=None):
    """What an embedding surface says of the moment it is cut from: the view, its settings,
    the surface's label and time, and where each piece runs, which is the reach of the moment
    every other diagram draws. A grid is its frame and the ends of its rows and columns. A moment
    of a stack is the ring `curve` it marks at its time, with that ring's label and time."""
    surface = view["surfaces"][index]
    pieces = []
    for piece in surface["pieces"]:
        if "grid" in piece:
            grid = piece["grid"]
            pieces.append([piece["id"], grid["frame"], grid["u"][0], grid["u"][-1], grid["v"][0], grid["v"][-1]])
        else:
            pieces.append([piece["id"], piece["class"], piece["metric"], piece["system"], piece["coordinate"],
                           piece["points"][0][0], piece["points"][-1][0]])
    out = {"view": view["id"], "settings": view.get("settings"), "label": surface.get("label"),
           "time": surface.get("time"), "pieces": pieces,
           "curves": [[c["class"], c["points"][0], c["points"][-1]] for c in surface.get("curves", [])
                      if c["class"] == "path"]}
    if curve is not None:
        ring = surface["curves"][curve]
        out["ring"] = [curve, ring["label"], ring["time"]]
    return out


def embedding_moment_version(view, index, curve=None):
    """The stamp a slice drawn on another diagram records of the embedding surface it marks."""
    return content_version(embedding_moment(view, index, curve))


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
    if grid.get("frame") == "ellipses":
        check_ellipses(at, piece)
        return
    if grid.get("frame") not in ("polar", "cartesian"):
        raise DataError(f"{at} is a grid in neither a polar nor a Cartesian frame, nor a stack of ellipses")
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


def check_ellipses(at, piece):
    """A stack of ellipses, as _tools/README.md defines it: its rows u and its angles v each
    strictly increasing, the angles within one turn, or up to a whole turn where it is `open`,
    each row's semi-axes a and b, never negative, and height z, and the lines it names by the
    index of a row or an angle, with an edge that says what lies beyond it. A stack whose ellipses
    are sheared gives every row its sx and sy as well, the row being the ellipse
    (a cos v + sx sin v, sy cos v + b sin v), and then a and b are entries of a matrix of positive
    determinant and may take either sign."""
    grid = piece["grid"]
    u, v = grid.get("u") or [], grid.get("v") or []
    for name, values in (("u", u), ("v", v)):
        if len(values) < 2 or not all(b > a for a, b in zip(values, values[1:])):
            raise DataError(f"{at} does not run one way along its {name}")
    if v[0] < 0 or (v[-1] > 2 * math.pi if grid.get("open") else v[-1] >= 2 * math.pi) or len(v) < 3:
        raise DataError(f"{at} runs round beyond one turn")
    if ("sx" in grid) != ("sy" in grid):
        raise DataError(f"{at} gives one half of its shear without the other")
    sheared = "sx" in grid
    rows = [grid.get(k) or [] for k in ("a", "b", "z") + (("sx", "sy") if sheared else ())]
    if any(len(row) != len(u) for row in rows):
        raise DataError(f"{at} does not give each row its ellipse and height")
    if not all(isinstance(h, (int, float)) and math.isfinite(h) for row in rows for h in row):
        raise DataError(f"{at} gives a semi-axis or a height that is not a number")
    if not sheared and any(h < 0 for row in rows[:2] for h in row):
        raise DataError(f"{at} has a semi-axis below zero")
    if sheared and any(a * b - sx * sy < 0 for a, b, _, sx, sy in zip(*rows)):
        raise DataError(f"{at} has a row that turns its ring inside out")
    for line in grid.get("lines", []):
        if not (all(0 <= i < len(u) for i in line["u"]) and all(0 <= j < len(v) for j in line["v"])):
            raise DataError(f"{at} names a line its grid does not hold")
    if any(key in piece for key in ("points", "start", "end", "coordinate")):
        raise DataError(f"{at} is a grid and a profile at once")
    if piece.get("edge", {}).get("kind") != "edge":
        raise DataError(f"{at} does not say what lies beyond its edge")


def check_moments(where, view):
    """A view that changes from moment to moment plays as a movie and is never set out as
    separate pictures, as the captain asked on 1 October 2026: more than one surface is a run of
    moments, each with its label and its time in order, and the view then carries a movie that
    holds every one of them as a frame, from the first to the last."""
    surfaces = view["surfaces"]
    if len(surfaces) < 2:
        return
    times = [s.get("time") for s in surfaces]
    if not all(isinstance(t, (int, float)) for t in times) or not all(b > a for a, b in zip(times, times[1:])):
        raise DataError(f"{where}: the moments of the view {view['id']!r} are not in the order of their times")
    if "movie" not in view:
        raise DataError(f"{where}: the view {view['id']!r} holds {len(surfaces)} moments as separate pictures and "
                        "no movie to play them")
    frames = {f.get("value"): f for f in view["movie"].get("frames") or []}
    for surface in surfaces:
        frame = frames.get(surface["time"])
        if frame is None or any(frame.get(k) != surface.get(k) for k in ("label", "pieces", "rings", "curves", "dots")):
            raise DataError(f"{where}: the movie of the view {view['id']!r} does not hold its moment "
                            f"{surface.get('label')!r} as a frame")
    values = sorted(frames)
    if (values[0], values[-1]) != (times[0], times[-1]):
        raise DataError(f"{where}: the movie of the view {view['id']!r} does not run from its first moment to its last")


def check_movie(where, view):
    """A movie, as _tools/README.md defines it: at least two frames, each a surface with its label
    and a value of the movie's variable, the values strictly increasing, one pass taking a
    positive number of seconds, and no `loop`, since every movie plays forward and back."""
    movie = view["movie"]
    frames = movie.get("frames") or []
    values = [f.get("value") for f in frames]
    if len(frames) < 2 or not all(isinstance(x, (int, float)) for x in values) or \
            not all(b > a for a, b in zip(values, values[1:])):
        raise DataError(f"{where}: the movie of the view {view['id']!r} does not run through its frames in order")
    if not all(isinstance(f.get("label"), str) and f["label"] for f in frames) or not movie.get("variable"):
        raise DataError(f"{where}: the movie of the view {view['id']!r} leaves a frame or its variable unnamed")
    if not (isinstance(movie.get("seconds"), (int, float)) and movie["seconds"] > 0):
        raise DataError(f"{where}: the movie of the view {view['id']!r} takes no time")
    if "loop" in movie:
        raise DataError(f"{where}: the movie of the view {view['id']!r} says how it loops, which no movie does: "
                        "every movie plays forward and back")
    if movie.get("turns", False) is not False:
        raise DataError(f"{where}: the movie of the view {view['id']!r} says it turns, which every figure does unless "
                        "it says otherwise")


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
            check_moments(where, view)
            if "movie" in view:
                check_movie(where, view)
            for surface in view["surfaces"] + view.get("movie", {}).get("frames", []):
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
            rings = target["surfaces"][mark["surface"]].get("curves", [])
            if "curve" in mark and not (0 <= mark["curve"] < len(rings) and "time" in rings[mark["curve"]]):
                raise DataError(f"{where} marks the ring {mark['curve']!r} of the embedding view {mark['view']!r}, "
                                f"which embedding/{metric_id}.json does not mark at a time")
            if embedding_moment_version(target, mark["surface"], mark.get("curve")) != mark.get("version"):
                raise DataError(f"{where} marks a moment of embedding/{metric_id}.json as it no longer "
                                f"stands; redraw it with _tools/derivations/{script} --metric {metric_id}")
            if not (mark.get("lines") or mark.get("points") or mark.get("fills")):
                raise DataError(f"{where} marks the moment {mark.get('label')!r} and draws nothing of it")


def check_shades(conformal, embedding):
    """Every shade of an embedding view answers a view of its spacetime's conformal diagram, by
    that view's id, and tints a piece of the view's own surfaces between two points of its
    profile, so the page never shades for a conformal view it cannot show or a circle the surface
    does not have."""
    for metric_id, data in (embedding or {}).items():
        where = f"embedding/{metric_id}.json"
        views = {v["id"] for v in (conformal or {}).get(metric_id, {}).get("views", [])}
        for view in data.get("views", []):
            for shade in view.get("shades", []):
                at = f"{where}, the shade of the view {view['id']!r} for {shade.get('view')!r}"
                if shade.get("view") not in views:
                    raise DataError(f"{at} answers a view that conformal/{metric_id}.json does not draw")
                surfaces = view["surfaces"]
                if not 0 <= shade.get("surface", -1) < len(surfaces):
                    raise DataError(f"{at} names the surface {shade.get('surface')!r}, which the view does not draw")
                piece = {p["id"]: p for p in surfaces[shade["surface"]]["pieces"]}.get(shade.get("piece"))
                if piece is None or "points" not in piece:
                    raise DataError(f"{at} names the piece {shade.get('piece')!r}, which is no profile of its surface")
                xs = [point[0] for point in piece["points"]]
                if shade.get("from") not in xs or shade.get("to") not in xs:
                    raise DataError(f"{at} ends at {shade.get('from')!r} and {shade.get('to')!r}, "
                                    f"which are not both points of the profile {piece['id']!r}")


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


# What one spacetime is to another, as the `kind` of an entry of `related` names it, and the
# kind the other spacetime's entry has to carry back. A kind says what the spacetime named is
# to the one listing it: Kerr is a `generalisation` on Schwarzschild's page, and Schwarzschild
# a `special_case` on Kerr's. A `special_case` is reached at a value of a parameter, and a
# `limit` only as a parameter runs to the end of its range with the coordinates rescaled on
# the way, as the Aichelburg-Sexl shock is of Schwarzschild's field; the spacetime it is a
# limit of is its `limit_source`. _tools/README.md says what each one covers.
RELATION_KINDS = {
    "special_case": "generalisation",
    "generalisation": "special_case",
    "limit": "limit_source",
    "limit_source": "limit",
    "piece": "composite",
    "composite": "piece",
    "family": "family",
    "dual": "dual",
    "conformal": "conformal",
    "locally_same": "locally_same",
    "programme": "programme",
}
RELATION_SENTENCES = 3
# A relation that is a limit says so, by the word or by an arrow such as $m \to \infty$.
LIMIT_KINDS = ("limit", "limit_source")
LIMIT_SAID = re.compile(r"\blimits?\b|\\to\b")
CITATION = re.compile(r"\[([A-Za-z0-9_, ]+)\]")


def cited_keys(text):
    """The reference keys a prose field cites, in the order it first cites them."""
    keys = []
    for bracket in CITATION.findall(text or ""):
        for key in (k.strip() for k in bracket.split(",")):
            if key and key not in keys:
                keys.append(key)
    return keys


def relation_problems(metrics):
    """Each way the relations between the spacetimes fail to hold together."""
    by_id = {metric["id"]: metric for metric in metrics}
    kinds = {}
    problems = []
    for metric in metrics:
        name = f"{metric['id']}.json"
        related = metric.get("related")
        if not isinstance(related, list) or not related:
            problems.append(f"{name} lists no related spacetimes")
            continue
        for position, entry in enumerate(related):
            where = f"{name}: related[{position}]"
            if not isinstance(entry, dict) or set(entry) != {"id", "kind", "text"}:
                problems.append(f"{where} is not an id, a kind and a text and nothing else")
                continue
            other, kind, text = entry["id"], entry["kind"], entry["text"]
            where = f"{name}: related[{other}]"
            if other == metric["id"]:
                problems.append(f"{name} lists itself as related")
                continue
            if other not in by_id:
                problems.append(f"{where} names a spacetime that has no metric file")
                continue
            if (metric["id"], other) in kinds:
                problems.append(f"{name} lists {other} twice")
                continue
            if kind not in RELATION_KINDS:
                problems.append(f"{where} has the kind {kind!r}, which is not one of "
                                f"{', '.join(sorted(RELATION_KINDS))}")
                continue
            kinds[metric["id"], other] = kind
            if not isinstance(text, str) or not text.strip():
                problems.append(f"{where} does not say why")
                continue
            if not 1 <= len(sentences(CITATION.sub("", text))) <= RELATION_SENTENCES:
                problems.append(f"{where} runs past {RELATION_SENTENCES} sentences")
            if text != text.strip() or not re.search(r"[.?!][\"')]*(?: \[[A-Za-z0-9_, ]+\]\.?)?$", text):
                problems.append(f"{where} does not end at the end of a sentence")
            if kind in LIMIT_KINDS and not LIMIT_SAID.search(text):
                problems.append(f"{where} is filed as a {kind} and does not say so, by the word limit "
                                "or by an arrow such as $m \\to \\infty$")
            for key in cited_keys(text):
                if key not in metric.get("references", []):
                    problems.append(f"{where} cites {key!r}, which {name} does not list in its references")
    for (one, other), kind in sorted(kinds.items()):
        back = kinds.get((other, one))
        if back is None:
            problems.append(f"{one}.json lists {other}, and {other}.json does not list {one}")
        elif back != RELATION_KINDS[kind]:
            problems.append(f"{one}.json lists {other} as {kind}, so {other}.json has to list {one} as "
                            f"{RELATION_KINDS[kind]} and lists it as {back}")
    return problems


def build_relations(metrics):
    """The relations the graph behind the list draws, MFS/assets/graph.js: every two spacetimes
    that list each other under `related`, once, as their two ids in order, the pairs in order.
    check_relations holds every relation to running both ways, so a pair either file lists is
    a pair of the graph."""
    pairs = {tuple(sorted((metric["id"], entry["id"]))) for metric in metrics for entry in metric["related"]}
    return {"edges": [list(pair) for pair in sorted(pairs)]}


def serialise_relations(relations):
    """The relations with one pair to a line, which is short to read and to compare."""
    lines = ",\n".join("    " + json.dumps(pair, ensure_ascii=False) for pair in relations["edges"])
    return '{\n  "edges": [\n' + lines + "\n  ]\n}\n"


def check_relations(metrics):
    """Refuse a relation that leads nowhere, runs one way only, or disagrees with its reverse."""
    problems = relation_problems(metrics)
    if problems:
        raise DataError("\n".join(problems))


# Prose set in paragraphs reads evenly: every paragraph of a history has three to six
# sentences, and the longest is no more than twice the length of the shortest. A history
# tells a story, so it runs to at least five paragraphs. A table is not a paragraph of prose
# and is left out of the count.
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
    """Each paragraph of a history that is prose, with its place among all of them."""
    return [(n, p) for n, p in enumerate(text.split("¶"), 1) if not p.startswith("TABLE::")]


def paragraph_shape(text):
    """The number of sentences in each paragraph of a history, tables left out."""
    return [len(sentences(p)) for _, p in prose_paragraphs(text)]


def shape_problems(metric_id, field, text, fewest):
    """Each way one history fails to read evenly, naming the paragraph at fault."""
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
    """Refuse a history with paragraphs of uneven length, or one too short to read as a story.

    Every one out of shape is named at once, with each paragraph that breaks the rule.
    """
    problems = []
    for metric in metrics:
        problems += shape_problems(metric["id"], "history", metric.get("history") or "",
                                   HISTORY_PARAGRAPHS)
    if problems:
        raise DataError("\n".join(problems))


# A convention says only what a reader needs to read the mathematics: the coordinates and
# their units, what each parameter and function is, the factors of c and G, and any choice
# of chart, index, sign or normalisation the components depend on. Each chart carries its
# own, and the spacetime's `convention` holds what is true in every chart and names no
# coordinate. The reader sees the chart's text followed by the spacetime's, as one paragraph
# of at most five sentences. The captain cut every convention down to that on 29 September
# 2026, and named the words below, which never come back.
CONVENTION_SENTENCES = 5
CONVENTION_BANNED = (
    r"(?i)\bcost",
    r"(?i)nowhere does the geometry break",
    r"(?i)standing objection",
    r"(?i)the Riemann tensor with the traces removed",
    r"(?i)\bslots?\b",
)


def chart_conventions(metric):
    """Each chart's convention as the reader sees it, with the chart's id: the chart's own
    text, then the text the spacetime shares across its charts. A spacetime with no charts
    has the shared text alone, under no id."""
    shared = metric.get("convention") or ""
    charts = metric.get("coordinates") or [{"id": None}]
    return [(chart["id"], " ".join(text for text in (chart.get("convention") or "", shared) if text))
            for chart in charts]


def convention_problems(metric):
    """Each way one spacetime's conventions break the rule, naming the chart at fault."""
    where = f"{metric['id']}.json"
    problems = []
    if not isinstance(metric.get("convention"), str):
        problems.append(f"{where} has no convention shared by its charts")
    for chart in metric.get("coordinates") or []:
        if not isinstance(chart.get("convention"), str) or not chart["convention"].strip():
            problems.append(f"{where}: the chart {chart['id']!r} has no convention of its own")
    for chart_id, text in chart_conventions(metric):
        at = f"{where}: the convention" + (f" of the chart {chart_id!r}" if chart_id else "")
        if "¶" in text:
            problems.append(f"{at} runs to more than one paragraph")
        count = len(sentences(text))
        if count > CONVENTION_SENTENCES:
            problems.append(f"{at} has {count} sentences, where a convention takes at most "
                            f"{CONVENTION_SENTENCES}")
        for pattern in CONVENTION_BANNED:
            found = re.search(pattern, text)
            if found:
                problems.append(f"{at} says {found.group(0)!r}")
    return problems


def check_conventions(metrics):
    """Refuse a convention that says more than a reader needs to read the mathematics, as far
    as a count and a list of words can tell, naming every one at once."""
    problems = [problem for metric in metrics for problem in convention_problems(metric)]
    if problems:
        raise DataError("\n".join(problems))


def check_symbols(metrics, diagrams, conformal, embedding):
    """Refuse a symbol that a chart's mathematics or a drawing uses and nothing defines, naming
    every one at once. `symbols.py` holds the rule and the list of symbols that need no definition."""
    problems = [problem for metric in metrics
                for problem in symbol_problems(metric, diagrams.get(metric["id"]), conformal.get(metric["id"]),
                                               embedding.get(metric["id"]))]
    if problems:
        raise DataError("\n".join(problems))


# A `$` opens or closes mathematics wherever the page sets a text, so one left over, as the
# `$` after "axis" in a Kastor-Traschen parameter was on 2 October 2026, turns the words after
# it into mathematics and prints the mathematics after them as its source. No text holds a
# dollar that is not a delimiter, so every string of every file is held to it, paragraph by
# paragraph, since the page sets each paragraph of a history on its own.
DISPLAY_MATH = re.compile(r"\$\$.+?\$\$", re.DOTALL)
DOLLAR = re.compile(r"(?<!\\)\$")


def unbalanced_dollar(text):
    """Whether a text opens mathematics it does not close, or closes what it did not open."""
    return any(len(DOLLAR.findall(DISPLAY_MATH.sub(" ", paragraph))) % 2 for paragraph in text.split("¶"))


def strings(value, where):
    """Yield every string held anywhere in a parsed JSON value, each with its place."""
    if isinstance(value, str):
        yield where, value
    elif isinstance(value, dict):
        for key, item in value.items():
            yield from strings(item, f"{where}.{key}")
    elif isinstance(value, list):
        for position, item in enumerate(value):
            yield from strings(item, f"{where}[{position}]")


def dollar_problems(name, data):
    """Each string of one parsed file with an unbalanced `$`, named by its place."""
    return [f"{where} has an unbalanced $: {text}" for where, text in strings(data, name) if unbalanced_dollar(text)]


def check_dollars():
    """Refuse a metric or diagram file with an unbalanced `$` in any of its strings, naming each."""
    problems = []
    for folder in (METRICS_DIR, DIAGRAMS_DIR, CONFORMAL_DIR, EMBEDDING_DIR):
        for path in sorted(folder.glob("*.json")):
            if not CONFLICT_COPY.search(path.stem):
                problems += dollar_problems(path.name, json.loads(path.read_text(encoding="utf-8")))
    if problems:
        raise DataError("\n".join(problems))


# The widest and the tallest a spacetime diagram's plotted region may be, as width over height.
ASPECT_WIDEST = 2.0
ASPECT_TALLEST = 0.5


def view_aspect(view):
    """The width over the height of the region a view of a spacetime diagram plots, as the page
    and the application draw it: its box at one scale on both axes, a view drawn through a centre
    twice as wide as the half its file holds."""
    x0, x1, y0, y1 = view["box"]
    width = 2 * x1 if view.get("mirror") else x1 - x0
    return width / (y1 - y0)


def aspect_problems(name, data):
    """Each view and figure of one parsed diagram file whose plotted region is wider than 2:1 or
    taller than 1:2, named by its place. The captain asked on 2 October 2026 that the spacetime
    diagrams "stick to 1:1 aspect ratios, allowing for up to 1:2 and 2:1 but no more than that"."""
    problems = []
    for part in ("systems", "projections"):
        for system_id, views in data.get(part, {}).items():
            for view in views:
                aspect = view_aspect(view)
                if not ASPECT_TALLEST - 1e-9 <= aspect <= ASPECT_WIDEST + 1e-9:
                    shape = f"{aspect:.2f}:1" if aspect > 1 else f"1:{1 / aspect:.2f}"
                    problems.append(
                        f"diagrams/{name}.json {system_id}/{view['id']} plots a region of {shape}, outside "
                        "1:2 to 2:1; choose a window nearer 1:1 in its row and redraw it")
    return problems


def check_aspects(diagrams):
    """Refuse a spacetime diagram whose plotted region is outside 1:2 to 2:1, naming each."""
    problems = [problem for name, data in sorted(diagrams.items()) for problem in aspect_problems(name, data)]
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
        check_dollars()
        check_prose_shape(metrics)
        check_conventions(metrics)
        check_relations(metrics)
        diagrams = load_diagrams(metrics)
        check_aspects(diagrams)
        conformal = load_conformal(metrics)
        embedding = load_embedding(metrics)
        check_slices(diagrams, conformal, embedding)
        check_shades(conformal, embedding)
        check_symbols(metrics, diagrams, conformal, embedding)
        outputs = {
            INDEX_FILE: serialise(build_index(metrics, diagrams, conformal, embedding)),
            REFERENCES_FILE: serialise(references),
            RELATIONS_FILE: serialise_relations(build_relations(metrics)),
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
        print("index, bibliography, relations and diagram stamps are up to date")
    return 0


if __name__ == "__main__":
    sys.exit(main())
