#!/usr/bin/env python3
"""The curvature figures hold the published scalars, value for value.

    .venv.noindex/bin/python -m unittest discover -s _tools

Every value a curvature figure carries is evaluated again here from the metric file it names,
the published `kretschmann` or `ricci_scalar` read through the checker's Reader at the point,
the other coordinates a curve holds fixed and the parameters the file gives, with each function
the chart leaves free put in as the closed form the file declares for it, and held to the
number written: a curve's value to the last few digits of a double, a map's cell to the six
significant figures it is written to. A map that declares the domain of its chart leaves null
every cell whose centre lies outside it, and no other. Every horizon marked on a line view is a
change of sign of the published g^xx along the drawn coordinate, and every change of sign on the
drawn range is marked, the range of a chart joined to another at a surface, as a star's to Schwarzschild's
vacuum, being the stretch its own curves cover. The file's shape is build_mfs_data.py's to
refuse, and the physics its captions state is curvature.py's to check before it writes.
"""

import json
import math
import sys
import unittest
from pathlib import Path

import build_mfs_data as build

try:
    import numpy as np
    import sympy as sp
except ImportError:  # pragma: no cover - the suite's environment check stops the run first
    np = sp = None

DERIVATIONS = Path(__file__).resolve().parent / "derivations"


def curvature_files():
    return {p.stem: json.loads(p.read_text(encoding="utf-8")) for p in sorted(build.CURVATURE_DIR.glob("*.json"))}


class Evaluator:
    """Published values of one chart, read through the checker's Reader as numpy functions of a
    chart's coordinates and parameters, kept by what they read."""

    def __init__(self):
        sys.path.insert(0, str(DERIVATIONS))
        import null_rays
        self.nr = null_rays
        self.made = {}

    def __call__(self, metric_id, system_id, field, coords, params, indices=None, functions=None, expression=None):
        functions = functions or {}
        key = (metric_id, system_id, field, tuple(coords), tuple(params), tuple(indices or ()),
               tuple(sorted(functions.items())), expression)
        if key not in self.made:
            _, entry, reader = self.nr.load(metric_id, system_id)
            if expression is not None:
                names = {reader._plain(n): s for n, s in reader.symbol.items()}
                names.update(reader.parameters)
                expr = sp.sympify(expression, locals={**self.nr.DECLARED_FUNCTIONS, **names})
            elif field == "inverse_metric_components":
                value = next(c["value"] for c in entry[field] if c["indices"] == list(indices))
                expr = reader(value)
            else:
                expr = reader(self.nr.strip_lhs(entry[field]))
            if reader.held:
                expr = expr.subs(reader.held).doit()
            names = {reader._plain(n): s for n, s in reader.symbol.items()}
            names.update(reader.parameters)
            for name, text in functions.items():
                rep = sp.sympify(text, locals={**self.nr.DECLARED_FUNCTIONS, **names})
                expr = expr.subs(reader.parameters[name], rep).doit()
            expr = expr.subs(reader.c, 1)
            by_plain = {reader._plain(name): symbol for name, symbol in reader.symbol.items()}
            symbols = [by_plain[reader._plain(c)] for c in coords] + [reader.parameters[reader._plain(p)] for p in params]
            stray = expr.free_symbols - set(symbols)
            free = {f for f in reader.parameters.values() if isinstance(f, sp.Function)} & expr.atoms(sp.Function)
            if free:
                raise AssertionError(f"{metric_id}/{system_id} {field} holds {free}, which the figure does not declare")
            if stray:
                raise AssertionError(f"{metric_id}/{system_id} {field} depends on {stray}, which the figure does not set")
            fn = sp.lambdify(symbols, sp.powsimp(expr, combine="exp"), self.nr.numeric_modules(expr))

            def evaluate(*values, fn=fn):
                with np.errstate(all="ignore"):
                    out = fn(*values)
                shape = np.broadcast(*[np.asarray(v, dtype=float) for v in values]).shape
                return np.array(np.broadcast_to(np.asarray(out, dtype=float), shape))
            self.made[key] = evaluate
        return self.made[key]

    def domain(self, metric_id, system_id, text, coords, params):
        """The domain a map declares for its chart, a closed form in sympy's spelling of the chart's
        coordinates and parameters, as a numpy function that is true where the chart holds."""
        _, _, reader = self.nr.load(metric_id, system_id)
        names = {reader._plain(n): s for n, s in reader.symbol.items()}
        names.update(reader.parameters)
        by_plain = {reader._plain(name): symbol for name, symbol in reader.symbol.items()}
        symbols = [by_plain[reader._plain(c)] for c in coords] + [reader.parameters[reader._plain(p)] for p in params]
        expr = sp.sympify(text, locals={**self.nr.DECLARED_FUNCTIONS, **names})
        stray = expr.free_symbols - set(symbols)
        if stray:
            raise AssertionError(f"{metric_id}/{system_id} the domain {text!r} depends on {stray}, which the map does not set")
        fn = sp.lambdify(symbols, expr, "numpy")
        return lambda *values: np.broadcast_to(np.asarray(fn(*values), dtype=bool),
                                               np.broadcast(*[np.asarray(v, dtype=float) for v in values]).shape)


def holders(view):
    """The view and each frame of its movie, with a name for each."""
    yield "", view
    for k, frame in enumerate(view.get("movie", {}).get("frames", [])):
        yield f" frame {k}", frame


def pictures(view):
    """What a line view draws at once, with a name for each: the view alone, or each frame of its
    movie with the curves and marks the view draws in every frame."""
    frames = view.get("movie", {}).get("frames", [])
    if not frames:
        yield "", view
    for k, frame in enumerate(frames):
        yield f" frame {k}", {"curves": view.get("curves", []) + frame.get("curves", []),
                              "marks": view.get("marks", []) + frame.get("marks", [])}


@unittest.skipIf(np is None, "numpy and sympy are not installed")
class CurvatureValues(unittest.TestCase):
    """Every value of every curvature figure is its published formula at its point."""

    @classmethod
    def setUpClass(cls):
        cls.files = curvature_files()
        cls.evaluate = Evaluator()

    def test_there_are_figures_to_check(self):
        self.assertTrue(self.files, "no curvature figure was found, so nothing was checked")

    def test_every_curve_is_its_published_scalar(self):
        checked = 0
        for name, data in self.files.items():
            for view in data["views"]:
                if view["kind"] != "line":
                    continue
                for place, holder in holders(view):
                    for k, curve in enumerate(holder.get("curves", [])):
                        where = f"curvature/{name}.json {view['id']}{place} curve {k}"
                        params, fixed = list(curve["params"]), curve.get("fixed", {})
                        fn = self.evaluate(curve["metric"], curve["system"], curve["scalar"],
                                           [curve["coordinate"]] + list(fixed), params,
                                           functions=curve.get("functions"))
                        points = np.array(curve["points"], dtype=float)
                        native = curve.get("horizontal", {}).get("values", points[:, 0])
                        arguments = [native, *fixed.values(), *[curve["params"][p] for p in params]]
                        expected = fn(*arguments)
                        if "horizontal" in curve:
                            horizontal = self.evaluate(curve["metric"], curve["system"], "horizontal",
                                                       [curve["coordinate"]] + list(fixed), params,
                                                       functions=curve.get("functions"),
                                                       expression=curve["horizontal"]["expression"])
                            with self.subTest(where + " horizontal"):
                                np.testing.assert_allclose(points[:, 0], horizontal(*arguments), rtol=1e-12, atol=1e-14)
                        with self.subTest(where):
                            np.testing.assert_allclose(points[:, 1], expected, rtol=1e-12, atol=0, err_msg=where)
                        checked += len(points)
        self.assertGreater(checked, 1000)

    def test_every_cell_is_its_published_scalar(self):
        checked = 0
        for name, data in self.files.items():
            for view in data["views"]:
                if view["kind"] != "map":
                    continue
                x, y = view["x"], view["y"]
                for place, holder in holders(view):
                    grid = holder.get("map")
                    if grid is None:
                        continue
                    where = f"curvature/{name}.json {view['id']}{place}"
                    params = list(grid["params"])
                    fn = self.evaluate(grid["metric"], grid["system"], grid["scalar"], [grid["x"], grid["y"]], params,
                                       functions=grid.get("functions"))
                    nx, ny = grid["nx"], grid["ny"]
                    xs = x["from"] + (np.arange(nx) + 0.5) * (x["to"] - x["from"]) / nx
                    ys = y["from"] + (np.arange(ny) + 0.5) * (y["to"] - y["from"]) / ny
                    X, Y = np.meshgrid(xs, ys)
                    if grid.get("frame") == "meridian":
                        # The box is the meridian plane, r sin theta across and r cos theta up.
                        X, Y = np.hypot(X, Y), np.arctan2(X, Y)
                    elif grid.get("frame") == "oblate":
                        # The box is the meridian plane of oblate spheroidal coordinates about a ring of
                        # radius f, sqrt(r^2 + f^2) sin theta across and r cos theta up, on the sheet r >= 0,
                        # r^2 taken inside the spheroid through the ring in the form that keeps its digits.
                        f = grid["params"][grid["focus"]]
                        s = X * X + Y * Y - f * f
                        q = np.sqrt(s * s + 4 * f * f * Y * Y)
                        r2 = np.where(s >= 0, (s + q) / 2, 2 * f * f * Y * Y / np.where(s >= 0, 1, q - s))
                        X, Y = np.sqrt(r2), np.arctan2(X * np.sqrt(r2), Y * np.sqrt(r2 + f * f))
                    values = [grid["params"][p] for p in params]
                    expected = fn(X, Y, *values)
                    written = np.array([[np.nan if v is None else v for v in row] for row in grid["cells"]], dtype=float)
                    drawn = ~np.isnan(written)
                    if "domain" in grid:
                        inside = self.evaluate.domain(grid["metric"], grid["system"], grid["domain"],
                                                      [grid["x"], grid["y"]], params)(X, Y, *values)
                        with self.subTest(where):
                            np.testing.assert_array_equal(drawn, inside, err_msg=f"{where}: a cell is null inside "
                                                          f"the domain {grid['domain']!r}, or a number outside it")
                    # Six significant figures: the value written is within half a unit of the sixth.
                    with self.subTest(where):
                        self.assertTrue(np.any(drawn), f"{where} leaves every cell null")
                        np.testing.assert_allclose(written[drawn], expected[drawn], rtol=5.1e-6, atol=0, err_msg=where)
                    checked += int(np.count_nonzero(drawn))
        self.assertGreater(checked, 1000)

    def test_every_horizon_marked_is_a_zero_of_the_published_inverse_metric_and_every_zero_is_marked(self):
        checked = 0
        for name, data in self.files.items():
            for view in data["views"]:
                if view["kind"] != "line":
                    continue
                for place, picture in pictures(view):
                    marks = [m["x"] for m in picture["marks"] if m["class"] == "horizon"]
                    # The charts the spacetime's own curves are drawn from, each at its parameters.
                    charts = {}
                    mapped = any("horizontal" in c for c in picture["curves"])
                    for curve in picture["curves"]:
                        # A mapped radial profile marks the vacuum comparison's horizon, including
                        # its continuation inside the matter. It is not a zero of g^chi_chi.
                        if ((mapped and curve["class"] == "reference") or
                                (not mapped and curve["class"] != "reference")):
                            key = (curve["metric"], curve["system"], curve["coordinate"],
                                   tuple(curve.get("fixed", {}).items()), tuple(sorted(curve["params"].items())),
                                   tuple(sorted((curve.get("functions") or {}).items())))
                            charts.setdefault(key, []).append(curve)
                    if not charts:
                        continue
                    where = f"curvature/{name}.json {view['id']}{place}"
                    changes = []
                    for (metric, system, coordinate, fixed, params, functions), curves in charts.items():
                        fn = self.evaluate(metric, system, "inverse_metric_components",
                                           [coordinate] + [c for c, _ in fixed], [p for p, _ in params],
                                           (coordinate, coordinate), dict(functions))
                        values = [v for _, v in fixed] + [v for _, v in params]
                        # A chart alone holds across the whole axis; a chart joined to another at a
                        # surface holds where its own curves are drawn.
                        lo, hi = view["x"]["from"], view["x"]["to"]
                        if len(charts) > 1:
                            lo = min(c["points"][0][0] for c in curves)
                            hi = max(c["points"][-1][0] for c in curves)
                        # The centres of 20000 even steps, so that no sample stands on a zero written
                        # exactly, as Schwarzschild's horizon r = r_s is at g = 0.
                        xs = lo + (np.arange(20000) + 0.5) * (hi - lo) / 20000
                        g = fn(xs, *values)
                        ok = np.isfinite(g[:-1]) & np.isfinite(g[1:]) & (np.sign(g[:-1]) * np.sign(g[1:]) < 0)
                        changes += list(xs[:-1][ok])
                        for at in marks:
                            if lo <= at <= hi:
                                below, above = fn(np.array([at - 1e-9, at + 1e-9]), *values)
                                with self.subTest(where):
                                    self.assertLess(below * above, 0, f"{where}: g^xx does not change sign at {at}")
                    with self.subTest(where):
                        self.assertEqual(len(changes), len(marks), f"{where} marks {marks}, g^xx changes sign near {changes}")
                        self.assertTrue(all(any(abs(at - c) < (view["x"]["to"] - view["x"]["from"]) / 10000 for c in changes)
                                            for at in marks), f"{where} marks {marks} away from the changes {changes}")
                    checked += 1
        self.assertGreater(checked, 0)


class CurvatureShape(unittest.TestCase):
    """What the build holds every file to, on the files on disk and on files made to break it."""

    def test_the_files_on_disk_load(self):
        self.assertEqual(set(build.load_curvature(build.load_metrics())), set(curvature_files()))

    def test_horizontal_coordinates_require_a_closed_expression_and_native_samples(self):
        import copy
        curve = {"coordinate": "\\chi", "fixed": {"\\tau": 0}, "params": {"a_m": 2},
                 "functions": {"a": "a_m - tau**2/4"}, "points": [[0, 1], [1, 1]],
                 "horizontal": {"expression": "a*sin(chi)", "values": [0, 0.5]}}
        build.check_horizontal("curve", curve, ["\\tau", "\\chi"])
        for change in (lambda h: h.update(expression="a*sin(theta)"),
                       lambda h: h.update(expression="a.__class__"),
                       lambda h: h.update(expression="sin(chi, 2)"),
                       lambda h: h.update(values=[0]),
                       lambda h: h.update(values=[0, 0]),
                       lambda h: h.update(values=[0, float("inf")])):
            broken = copy.deepcopy(curve)
            change(broken["horizontal"])
            with self.subTest(change), self.assertRaises(build.DataError):
                build.check_horizontal("curve", broken, ["\\tau", "\\chi"])

    def broken(self, change, name=None):
        import copy
        import tempfile
        from unittest import mock
        files = curvature_files()
        name = name or sorted(files)[0]
        data = copy.deepcopy(files[name])
        change(data)
        with tempfile.TemporaryDirectory() as folder:
            (Path(folder) / f"{name}.json").write_text(json.dumps(data), encoding="utf-8")
            with mock.patch.object(build, "CURVATURE_DIR", Path(folder)):
                with self.assertRaises(build.DataError) as raised:
                    build.load_curvature(build.load_metrics())
        return str(raised.exception)

    def test_a_changed_scalar_is_refused(self):
        metrics = build.load_metrics()
        for metric in metrics:
            if metric["id"] == "bardeen":
                for system in metric["coordinates"]:
                    system["kretschmann"] = system["kretschmann"].replace("47g^4", "48g^4")
        with self.assertRaises(build.DataError) as raised:
            build.load_curvature(metrics)
        self.assertIn("curvature/bardeen.json", str(raised.exception))
        self.assertIn("curvature.py --metric bardeen", str(raised.exception))

    def test_a_stretched_view_is_refused(self):
        def stretch(data):
            view = data["views"][0]
            if view["kind"] == "map":
                view["x"]["to"] *= 4
                view["x"]["ticks"] = [t for t in view["x"]["ticks"]]
            else:
                view["panels"][0]["aspect"] = 3
        self.assertIn("1:2 to 2:1", self.broken(stretch))

    def test_a_value_below_zero_on_a_logarithmic_scale_is_refused(self):
        def negate(data):
            view = data["views"][0]
            if view["kind"] == "map":
                view["map"]["cells"][0][0] = -1.0
            else:
                curve = next(c for c in view["curves"] if c["panel"] == "K")
                curve["points"][0][1] = -1.0
        self.assertIn("at or below zero", self.broken(negate))

    def test_a_map_read_in_a_frame_of_no_known_kind_is_refused(self):
        def reframed(data):
            view = next(v for v in data["views"] if v["kind"] == "map")
            view["map"]["frame"] = "polar"
        self.assertIn("frame", self.broken(reframed))

    def test_a_map_read_about_a_ring_of_no_named_radius_is_refused(self):
        files = curvature_files()
        self.assertIn("kerr_newman", files, "no map is read about a ring, so none was checked")
        for change in (lambda grid: grid.pop("focus"), lambda grid: grid.update(focus="b"),
                       lambda grid: grid.pop("frame")):
            def unfocused(data, change=change):
                change(next(v for v in data["views"] if v["kind"] == "map")["map"])
            with self.subTest(change):
                self.assertIn("ring", self.broken(unfocused, "kerr_newman"))

    def test_a_coordinate_held_at_no_number_is_refused(self):
        def unheld(data):
            view = next(v for v in data["views"] if v["kind"] == "line")
            curve = next(c for c in view["curves"] if "fixed" in c)
            curve["fixed"] = {name: "pi/2" for name in curve["fixed"]}
        self.assertIn("holds fixed", self.broken(unheld, "curzon_chazy"))

    def test_a_class_drawn_and_not_named_is_refused(self):
        def unnamed(data):
            data["views"][0]["legend"] = data["views"][0]["legend"][1:]
        self.assertIn("legend", self.broken(unnamed))

    def test_a_free_function_left_undeclared_is_refused(self):
        files = curvature_files()
        self.assertIn("bonnor_vaidya", files, "no figure declares a function, so none was checked")
        for change in (lambda grid: grid.pop("functions"), lambda grid: grid["functions"].pop("q"),
                       lambda grid: grid["functions"].update(q="")):
            def undeclared(data, change=change):
                change(data["views"][0]["map"])
            with self.subTest(change):
                self.assertIn("function", self.broken(undeclared, "bonnor_vaidya"))

    def test_a_null_cell_without_a_domain_is_refused(self):
        def emptied(data):
            view = next(v for v in data["views"] if v["kind"] == "map")
            view["map"].pop("domain", None)
            view["map"]["cells"][0][0] = None
        self.assertIn("domain", self.broken(emptied, "bach_weyl_ring"))

    def test_a_movie_that_says_how_it_loops_is_refused(self):
        def looped(data):
            view = next((v for v in data["views"] if "movie" in v), None)
            if view is None:
                view = data["views"][0]
                view["movie"] = {"variable": "$t$", "seconds": 5, "frames": [{"label": "$t = 0$", "value": 0},
                                                                             {"label": "$t = 1$", "value": 1}]}
            view["movie"]["loop"] = True
        self.assertIn("loops", self.broken(looped))


if __name__ == "__main__":
    unittest.main()
