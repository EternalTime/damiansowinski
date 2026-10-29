"""The text a view of a drawing shows beside it, counted one way for every drawing program.

A view's text is its caption, its clip note, the words of its legend and its settings, the
characters of each added together as the file writes them. A legend entry is [kind, class, text]
and only its text is shown; a view of the null rays has no legend of its own, and the page builds
one from the names of its two families of rays, the name of its cones where it draws any and the
legend each kind of marker carries. A view whose text runs past LIMIT characters is named in a warning, and the
program that drew it carries on and writes it.
"""
import sys

LIMIT = 400


def legend_words(view):
    """The words of the legend the page shows beside a view, in the order it shows them."""
    if "families" in view:
        words = list(view["families"])
        if view.get("cones") and view.get("cone"):
            words.append(view["cone"])
        # The page shows one entry for each kind of marker, so a kind drawn twice is counted once.
        legends = {}
        for marker in view.get("markers", []):
            if marker.get("legend"):
                legends[marker["kind"]] = marker["legend"]
        return words + list(legends.values())
    legend = view.get("legend", view.get("figure", {}).get("legend", []))
    return [entry[2] for entry in legend]


def length(view):
    """The characters of a view's caption, clip note, legend words and settings, added together."""
    texts = list(view["caption"]) + [view.get("clip") or "", view.get("settings") or ""] + legend_words(view)
    return sum(len(text) for text in texts)


def warn(where, views):
    """Warn of every view among `views` whose text runs past LIMIT, naming it as `where` and its id."""
    for view in views:
        count = length(view)
        if count > LIMIT:
            print(f"warning: {where} {view['id']} shows {count} characters of text, past {LIMIT}", file=sys.stderr)
