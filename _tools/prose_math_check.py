#!/usr/bin/env python3
"""Show that a rewording of the prose changed no mathematics, number or citation.

    python3 _tools/prose_math_check.py [revision]

Every text a reader sees on the spacetimes page is read twice, from the working tree and from
`revision` (main by default): each paragraph of every history and convention, and every text of
every spacetime diagram, conformal diagram and embedding diagram file, its captions, notes,
restriction bands, sentences stated in place of a drawing, labels and legends. With its words
taken out, what is left of each text is its mathematics, every `$...$` and `$$...$$` character for
character, its numbers written in digits and its citations, in order, and that must be the same in
both. It prints each text where they differ and exits non-zero if any does, or if the two do not
hold the same texts in the same places.
"""

import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = "MFS/assets/data"
MATH = re.compile(r"\$\$.+?\$\$|\$[^$]+\$", re.DOTALL)


def skeleton(text):
    """The mathematics, the numbers in digits and the citations of a text, in order."""
    maths = MATH.findall(text)
    words = MATH.sub(" ", text)
    return maths, re.findall(r"\d+(?:\.\d+)?", words), re.findall(r"\[[^\]]+\]", words)


def texts(value, where):
    """Every string inside a JSON value, each with its place."""
    if isinstance(value, str):
        yield where, value
    elif isinstance(value, list):
        for position, item in enumerate(value):
            yield from texts(item, f"{where}[{position}]")
    elif isinstance(value, dict):
        for key, item in value.items():
            if key in ("points", "source", "layers", "surfaces", "rays", "cones", "markers", "ticks", "hatch",
                       "curves", "figure", "version", "fields", "box"):
                continue
            yield from texts(item, f"{where}.{key}")


def collect(read, names):
    """Map each place to its text, reading files through `read`."""
    found = {}
    for name in names:
        data = json.loads(read(name))
        if name.startswith(f"{DATA}/metrics/"):
            for field in ("history", "convention"):
                for number, paragraph in enumerate((data.get(field) or "").split("¶")):
                    if paragraph:
                        found[f"{name} {field} paragraph {number + 1}"] = paragraph
        else:
            for where, text in texts(data, name):
                found[where] = text
            for part in ("systems", "projections"):
                for system, views in (data.get(part) or {}).items():
                    for view in views:
                        for key in ("caption", "input", "settings", "legend", "labels"):
                            if key in view:
                                for where, text in texts(view[key], f"{name} {part}/{system}/{view['id']}.{key}"):
                                    found[where] = text
            for view in data.get("views") or []:
                for key in ("figure",):
                    if key in view:
                        for where, text in texts({"legend": view[key].get("legend"), "labels": view[key].get("labels")},
                                                 f"{name} {view['id']}.{key}"):
                            found[where] = text
                for surface in view.get("surfaces") or []:
                    for piece in surface.get("pieces") or []:
                        for end in ("start", "end", "edge"):
                            if (piece.get(end) or {}).get("text"):
                                found[f"{name} {view['id']}.{piece['id']}.{end}"] = piece[end]["text"]
    return found


def main(argv):
    revision = argv[0] if argv else "main"
    kinds = ("metrics", "diagrams", "conformal", "embedding")
    here = sorted(str(p.relative_to(ROOT)) for kind in kinds for p in (ROOT / DATA / kind).glob("*.json"))
    listed = subprocess.run(["git", "ls-tree", "-r", "--name-only", revision, "--"] + [f"{DATA}/{k}" for k in kinds],
                            cwd=ROOT, capture_output=True, text=True, check=True).stdout.split()
    there = sorted(n for n in listed if n.endswith(".json"))
    if here != there:
        print(f"the files differ: only here {sorted(set(here) - set(there))}, only in {revision} "
              f"{sorted(set(there) - set(here))}")
        return 1
    now = collect(lambda n: (ROOT / n).read_text(encoding="utf-8"), here)
    then = collect(lambda n: subprocess.run(["git", "show", f"{revision}:{n}"], cwd=ROOT, capture_output=True,
                                            text=True, check=True).stdout, there)
    if set(now) != set(then):
        print(f"the texts differ in place: only here {sorted(set(now) - set(then))[:5]}, only in {revision} "
              f"{sorted(set(then) - set(now))[:5]}")
        return 1
    reworded = sum(now[k] != then[k] for k in now)
    wrong = [k for k in sorted(now) if skeleton(now[k]) != skeleton(then[k])]
    for where in wrong:
        print(f"{where}\n  {revision}: {skeleton(then[where])}\n  here: {skeleton(now[where])}")
    print(f"{len(now)} texts, {reworded} reworded against {revision}, "
          f"{len(wrong)} with their mathematics, numbers or citations changed")
    return 1 if wrong else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
