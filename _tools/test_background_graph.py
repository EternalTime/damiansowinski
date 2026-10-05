"""The graph of relations drawn behind the spacetimes page, MFS/assets/graph.js and the script
beside its canvas in _layouts/mfs.html, as the captain asked on 3 October 2026.

The graph is drawn from the site's own data: the list's index and the relations
build_mfs_data.py writes from each spacetime's `related`. What the search finds is what the
graph gathers, since the list and the graph are handed one answer, mfsSearch's, and clearing
the search sends every spacetime back to its place in the whole collection. A layout is the
same for the same data every time, and a search's layout depends on the spacetimes it found
and the relations between them alone. The geometry is run here in Node as the page runs it;
`node _tools/background_graph.mjs` holds the page as drawn to the same rules.
"""

import json
import re
import shutil
import subprocess
import unittest

import build_mfs_data as build

GRAPH = build.ROOT / "MFS" / "assets" / "graph.js"
PAGE = build.ROOT / "_layouts" / "mfs.html"
# The room the graph has beside the list on a desktop of 1440 by 900 at the usual text size,
# and on narrower windows, where fewer names have the room.
ROOMS = ((1140, 720), (700, 600), (420, 520))


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def page_function(source, name):
    """A function of the page's search script, as the page defines it."""
    found = re.search(r"      function " + name + r"\(.*?\n      \}\n", source, re.S)
    if not found:
        raise AssertionError(f"the page no longer defines {name}")
    return found.group(0)


def node(script, payload):
    """`script` run in Node with the graph's geometry as G and `payload` as input, its output read as JSON."""
    if shutil.which("node") is None:
        raise unittest.SkipTest("node is not installed")
    prelude = f"const G = require({json.dumps(str(GRAPH))}); const input = JSON.parse(require('fs').readFileSync(0, 'utf8'));\n"
    run = subprocess.run(["node", "-e", prelude + script], input=json.dumps(payload),
                         capture_output=True, text=True, timeout=600)
    if run.returncode != 0:
        raise AssertionError(run.stderr[-3000:])
    return json.loads(run.stdout)


def word_matches(tag, q):
    """A tag matched from the start of one of its words, written apart from the page's mfsTagMatches."""
    t = tag.lower()
    return any(t.startswith(q, i) and (i == 0 or t[i - 1] in " -/") for i in range(len(t)))


class Relations(unittest.TestCase):
    """The graph's relations are every two spacetimes that list each other, once, and nothing else."""

    def test_every_relation_is_one_pair_of_the_graph(self):
        metrics = build.load_metrics()
        listed = {tuple(sorted((m["id"], r["id"]))) for m in metrics for r in m["related"]}
        edges = read(build.RELATIONS_FILE)["edges"]
        self.assertEqual(len(edges), len(listed))
        self.assertEqual({tuple(e) for e in edges}, listed)
        self.assertEqual(edges, sorted(edges))
        for one, other in edges:
            self.assertLess(one, other)

    def test_the_relations_file_is_what_the_build_writes(self):
        self.assertEqual(build.RELATIONS_FILE.read_text(encoding="utf-8"),
                         build.serialise_relations(build.build_relations(build.load_metrics())))

    def test_the_graph_counts_each_spacetime_s_relations(self):
        index, relations = read(build.INDEX_FILE), read(build.RELATIONS_FILE)
        graph = node("process.stdout.write(JSON.stringify(G.model(input[0], input[1])));", [index, relations])
        self.assertEqual([s["id"] for s in graph["spacetimes"]], [m["id"] for m in index])
        degree = {m["id"]: 0 for m in index}
        for one, other in relations["edges"]:
            degree[one] += 1
            degree[other] += 1
        self.assertEqual({s["id"]: s["degree"] for s in graph["spacetimes"]}, degree)
        self.assertEqual(len(graph["edges"]), len(relations["edges"]))


class SearchToGraph(unittest.TestCase):
    """Typing a keyword gathers exactly the spacetimes the list then shows, the spacetimes carrying
    it, and clearing the search brings back the whole graph."""

    @classmethod
    def setUpClass(cls):
        cls.source = PAGE.read_text(encoding="utf-8")
        cls.index = read(build.INDEX_FILE)
        cls.relations = read(build.RELATIONS_FILE)
        cls.keywords = sorted({t for m in cls.index for t in m["tags"]})
        search = page_function(cls.source, "mfsTagMatches") + page_function(cls.source, "mfsSearch")
        # Each keyword typed as the reader types it, and what the page's search and the graph make of it.
        script = search + """
        const [index, relations, queries, W, H] = input;
        const graph = G.model(index, relations), all = graph.spacetimes.map((s, i) => i);
        const home = G.gathered(graph, all, W, H).places.map(p => ({ x: p[0], y: p[1], z: p[2] }));
        const out = {};
        for (const q of queries) {
          const found = mfsSearch(index, q), ids = found && found.map(m => m.id);
          const members = ids && G.placesOf(graph, ids);
          const goal = G.goal(home, members, members && G.gathered(graph, members, W, H).places);
          out[q] = { listed: (found || index).map(m => m.id), front: graph.spacetimes.filter((s, i) => !goal.behind[i]).map(s => s.id),
                     home: JSON.stringify(goal.to) === JSON.stringify(home) };
        }
        process.stdout.write(JSON.stringify(out));
        """
        typed = [k for k in cls.keywords] + [k.upper() for k in cls.keywords[:20]] + ["  vacuum  ", "", "   "]
        cls.results = node(script, [cls.index, cls.relations, typed, *ROOMS[0]])

    def carrying(self, q):
        return [m["id"] for m in self.index if q in m["name"].lower() or any(word_matches(t, q) for t in m["tags"])]

    def test_typing_a_keyword_gathers_the_spacetimes_carrying_it(self):
        for keyword in self.keywords:
            with self.subTest(keyword):
                result = self.results[keyword]
                self.assertEqual(result["front"], self.carrying(keyword.lower()))
                self.assertFalse(result["home"])
                exact = [m["id"] for m in self.index if keyword in m["tags"]]
                self.assertTrue(set(exact) <= set(result["front"]))

    def test_the_graph_gathers_what_the_list_shows(self):
        for typed, result in self.results.items():
            with self.subTest(typed):
                self.assertEqual(result["front"], result["listed"])

    def test_a_keyword_no_name_holds_gathers_its_carriers_and_nothing_else(self):
        # "conformally flat" is the keyword the captain searched for on 2 October 2026.
        for keyword in ("conformally flat", "wormhole", "Petrov type D", "cosmological constant"):
            with self.subTest(keyword):
                carriers = [m["id"] for m in self.index
                            if any(word_matches(t, keyword.lower()) for t in m["tags"])]
                self.assertTrue(carriers)
                self.assertEqual(self.results[keyword]["front"], carriers)

    def test_capitals_and_spaces_round_a_keyword_change_nothing(self):
        for keyword in self.keywords[:20]:
            with self.subTest(keyword):
                self.assertEqual(self.results[keyword.upper()]["front"], self.results[keyword]["front"])
        self.assertEqual(self.results["  vacuum  "]["front"], self.results["vacuum"]["front"])

    def test_clearing_the_search_brings_back_the_whole_graph(self):
        everyone = [m["id"] for m in self.index]
        for cleared in ("", "   "):
            with self.subTest(repr(cleared)):
                self.assertTrue(self.results[cleared]["home"])
                self.assertEqual(self.results[cleared]["front"], everyone)
                self.assertEqual(self.results[cleared]["listed"], everyone)

    def test_the_list_hands_the_graph_the_search_s_answer(self):
        render = page_function(self.source, "renderResults")
        self.assertIn("var found = keyword ? mfsCarrying(METRIC_INDEX, keyword) : mfsSearch(METRIC_INDEX, query);", render)
        self.assertIn("window._mfsGraph.show(METRIC_INDEX, found && found.map(function(m) { return m.id; }));", render)
        self.assertIn("var matches = found || METRIC_INDEX;", render)
        self.assertEqual(self.source.count("window._mfsGraph.show("), 1)
        self.assertEqual(self.source.count("mfsSearch("), 2)


class Keywords(unittest.TestCase):
    """The keywords offered under the search field as the reader types, from the first letter on,
    as the captain asked on 5 October 2026: the index's tags that the search finds for what is
    typed, those it starts first, and a keyword chosen selects exactly the spacetimes carrying it,
    in the list and in the graph. `node _tools/background_graph.mjs` holds the page as drawn to
    the same, by the pointer and by the keyboard."""

    @classmethod
    def setUpClass(cls):
        cls.source = PAGE.read_text(encoding="utf-8")
        cls.index = read(build.INDEX_FILE)
        cls.relations = read(build.RELATIONS_FILE)
        cls.keywords = sorted({t for m in cls.index for t in m["tags"]})
        starts = {t.lower()[i:i + 2] for t in cls.keywords for i in range(len(t))
                  if (i == 0 or t[i - 1] in " -/") and " " not in t[i:i + 2]}
        letters = sorted({q[0] for q in starts} | set("abcdefghijklmnopqrstuvwxyz"))
        cls.typed = letters + sorted(q for q in starts if len(q) == 2) + ["zq", "W", "De", " v ", "", "  "]
        functions = "".join(page_function(cls.source, name) for name in ("mfsTagMatches", "mfsSearch", "mfsKeywords", "mfsCarrying"))
        script = functions + """
        const [index, relations, typed, keywords, W, H] = input;
        const graph = G.model(index, relations), all = graph.spacetimes.map((s, i) => i);
        const home = G.gathered(graph, all, W, H).places.map(p => ({ x: p[0], y: p[1], z: p[2] }));
        const offered = {}, chosen = {};
        for (const q of typed) offered[q] = mfsKeywords(index, q);
        for (const k of keywords) {
          const ids = mfsCarrying(index, k).map(m => m.id), members = G.placesOf(graph, ids);
          const goal = G.goal(home, members, G.gathered(graph, members, W, H).places);
          chosen[k] = { listed: ids, front: graph.spacetimes.filter((s, i) => !goal.behind[i]).map(s => s.id) };
        }
        process.stdout.write(JSON.stringify({ offered, chosen }));
        """
        cls.results = node(script, [cls.index, cls.relations, cls.typed, cls.keywords, *ROOMS[0]])

    def expected(self, typed):
        """The keywords offered for `typed`, worked out apart from the page."""
        q = typed.strip().lower()
        if not q:
            return []
        found = [t for t in self.keywords if word_matches(t, q)]
        order = lambda t: (t.lower(), t)
        return (sorted((t for t in found if t.lower().startswith(q)), key=order) +
                sorted((t for t in found if not t.lower().startswith(q)), key=order))

    def test_one_letter_offers_every_keyword_with_a_word_it_starts(self):
        for typed in self.typed:
            if len(typed) != 1:
                continue
            with self.subTest(typed):
                self.assertEqual(self.results["offered"][typed], self.expected(typed))
        self.assertEqual(self.results["offered"]["v"],
                         ["vacuum", "vacuum energy", "Vaidya", "vanishing scalar invariants", "void",
                          "energy condition violation", "null Killing vector"])
        self.assertEqual(self.results["offered"]["W"], self.results["offered"]["w"])

    def test_two_letters_offer_every_keyword_with_a_word_they_start(self):
        two = [t for t in self.typed if len(t.strip()) == 2]
        self.assertGreater(len(two), 100)
        for typed in two:
            with self.subTest(typed):
                self.assertEqual(self.results["offered"][typed], self.expected(typed))
                if typed != "zq":
                    self.assertTrue(self.results["offered"][typed])
        self.assertEqual(self.results["offered"]["zq"], [])
        self.assertEqual(self.results["offered"]["De"][:2], ["de Sitter", "de Sitter core"])
        self.assertIn("anti-de Sitter", self.results["offered"]["De"])

    def test_nothing_is_offered_while_nothing_is_typed(self):
        self.assertEqual(self.results["offered"][""], [])
        self.assertEqual(self.results["offered"]["  "], [])
        self.assertEqual(self.results["offered"][" v "], self.results["offered"]["v"])

    def test_every_keyword_is_offered_by_its_first_two_letters(self):
        for keyword in self.keywords:
            with self.subTest(keyword):
                self.assertIn(keyword, self.results["offered"][keyword[:2].lower()])

    def test_a_keyword_chosen_selects_exactly_the_spacetimes_carrying_it(self):
        for keyword in self.keywords:
            with self.subTest(keyword):
                carrying = [m["id"] for m in self.index if keyword in m["tags"]]
                self.assertEqual(self.results["chosen"][keyword]["listed"], carrying)
                self.assertEqual(self.results["chosen"][keyword]["front"], carrying)
        # Typed out, "de Sitter" finds every spacetime with a tag or a name it starts a word of.
        self.assertLess(len(self.results["chosen"]["de Sitter"]["listed"]),
                        len([m for m in self.index if "de sitter" in m["name"].lower() or any(word_matches(t, "de sitter") for t in m["tags"])]))

    def test_choosing_a_keyword_runs_the_search_with_it(self):
        choose = self.source[self.source.index("        function choose(i) {"):][:300]
        self.assertIn("input.value = keyword;", choose)
        self.assertIn("renderResults(index, keyword, keyword);", choose)
        self.assertEqual(self.source.count("renderResults("), 4)
        typing = self.source[self.source.index("        input.addEventListener('input', function() {"):][:200]
        self.assertIn("renderResults(index, this.value);", typing)
        self.assertIn("offer(mfsKeywords(index, this.value));", typing)


class SearchBar(unittest.TestCase):
    """The search field keeps the site's own look whatever is typed in it. On 5 October 2026 the
    captain saw it change colour and turn white while typing, which nothing of the page's own
    does, the field's colours read the same in every state: what does is the browser offering
    what was typed in the field before, under it and in its own colours, and filling the field
    in them. The browser now offers nothing of its own there, and no rule of the page sets the
    field apart in any state. `node _tools/background_graph.mjs` reads its colours in every state."""

    @classmethod
    def setUpClass(cls):
        cls.source = PAGE.read_text(encoding="utf-8")
        cls.field = re.search(r'<input id="mfs-search-input".*?>', cls.source, re.S).group(0)

    def test_the_browser_offers_nothing_of_its_own_in_the_field(self):
        for attribute in ('autocomplete="off"', 'autocorrect="off"', 'autocapitalize="off"', 'spellcheck="false"'):
            with self.subTest(attribute):
                self.assertIn(attribute, self.field)
        self.assertIn('type="text"', self.field)

    def test_the_field_is_a_combobox_of_the_keywords(self):
        for attribute in ('role="combobox"', 'aria-autocomplete="list"', 'aria-controls="mfs-keywords"', 'aria-expanded="false"'):
            with self.subTest(attribute):
                self.assertIn(attribute, self.field)
        self.assertIn('<div id="mfs-keywords" role="listbox" aria-label="Keywords" hidden></div>', self.source)

    def test_no_rule_sets_the_field_apart_in_any_state(self):
        sheets = [self.source] + [(build.ROOT / "assets" / "css" / name).read_text(encoding="utf-8") for name in ("mfs.css", "palette.css")]
        for sheet in sheets:
            self.assertEqual(re.findall(r"#mfs-search-input:(?!:placeholder)[^{]*\{", sheet), [])
            self.assertEqual(re.findall(r"(?<![\w-])input(?:\[[^\]]*\])?:[\w-]+[^{]*\{", sheet), [])

    def test_the_keywords_are_set_as_the_list_is_with_nothing_new(self):
        self.assertIn("\n    .mfs-result, .mfs-keyword {\n", self.source)
        self.assertIn("\n    .mfs-result:hover, .mfs-keyword.mfs-keyword-on {\n", self.source)
        rules = re.findall(r"\n    ([^{}\n]*mfs-keyword[^{}\n]*)\{(.*?)\}", self.source, re.S)
        self.assertTrue(rules)
        for selectors, body in rules:
            with self.subTest(selectors):
                for glow in ("shadow", "filter", "glow"):
                    self.assertNotIn(glow, body)


class Layout(unittest.TestCase):
    """A layout is the same for the same data, a search's layout is its spacetimes' alone, and
    every layout keeps the points apart and the names to the forty the room holds."""

    @classmethod
    def setUpClass(cls):
        cls.index = read(build.INDEX_FILE)
        cls.relations = read(build.RELATIONS_FILE)

    def layouts(self, members_by_name, rooms):
        script = """
        const [index, relations, groups, rooms] = input;
        const graph = G.model(index, relations), out = {};
        for (const [name, ids] of Object.entries(groups)) {
          const members = ids ? G.placesOf(graph, ids) : graph.spacetimes.map((s, i) => i);
          out[name] = { arranged: G.arranged(graph, members), rooms: rooms.map(([W, H]) => G.gathered(graph, members, W, H)) };
        }
        process.stdout.write(JSON.stringify(out));
        """
        return node(script, [self.index, self.relations, members_by_name, rooms])

    def test_the_same_data_gives_the_same_layout_every_time(self):
        groups = {"all": None, "vacuum": [m["id"] for m in self.index if "vacuum" in m["tags"]]}
        first, second = self.layouts(groups, ROOMS), self.layouts(groups, ROOMS)
        self.assertEqual(first, second)

    def test_a_search_s_layout_is_its_spacetimes_and_their_relations_alone(self):
        members = [m["id"] for m in self.index if "wormhole" in m["tags"]]
        script = """
        const [index, relations, members] = input;
        const whole = G.model(index, relations);
        const own = new Set(members);
        // The same spacetimes with every relation to a spacetime outside them taken away, and
        // with the rest of the collection taken away.
        const cut = G.model(index, { edges: relations.edges.filter(e => own.has(e[0]) && own.has(e[1])) });
        const alone = G.model(index.filter(m => own.has(m.id)), relations);
        const ids = g => G.placesOf(g, members);
        process.stdout.write(JSON.stringify([G.arranged(whole, ids(whole)), G.arranged(cut, ids(cut)), G.arranged(alone, ids(alone))]));
        """
        whole, cut, alone = node(script, [self.index, self.relations, members])
        self.assertEqual(len(whole), len(members))
        self.assertEqual(whole, cut)
        self.assertEqual(whole, alone)

    def test_the_whole_graph_fills_its_room_with_its_points_apart_and_forty_names_at_most(self):
        script = """
        const [index, relations, rooms] = input;
        const graph = G.model(index, relations), all = graph.spacetimes.map((s, i) => i);
        let most = 0;
        graph.spacetimes.forEach(s => { most = Math.max(most, s.degree); });
        const out = rooms.map(([W, H]) => {
          const laid = G.gathered(graph, all, W, H), camera = G.start();
          const items = laid.places.map((p, i) => {
            const at = G.seen(camera, { x: p[0], y: p[1], z: p[2] }, W, H), f = G.drawnAt(at.scale), s = graph.spacetimes[i];
            const lines = G.twoLines(s.name);
            return { x: at.x, y: at.y, near: at.near, weight: s.degree, f: f, r: G.radius(f, s.degree, most),
                     w: Math.max(...lines.map(l => l.length)) * G.LETTER * f, h: lines.length * G.LINE * f,
                     rank: 2, must: false, side: -1 };
          });
          const sides = G.labels(items, W, H);
          return { named: laid.named, written: sides.filter(s => s >= 0).length, apart: G.apartAll(items),
                   inside: items.every(it => it.x - it.r >= 0 && it.x + it.r <= W && it.y - it.r >= 0 && it.y + it.r <= H) };
        });
        process.stdout.write(JSON.stringify(out));
        """
        results = node(script, [self.index, self.relations, ROOMS])
        for (W, H), result in zip(ROOMS, results):
            with self.subTest(room=(W, H)):
                self.assertTrue(result["apart"])
                self.assertTrue(result["inside"])
                self.assertLessEqual(result["written"], 40)
                self.assertGreaterEqual(result["written"], result["named"])
        self.assertEqual(results[0]["named"], 40)

    def test_a_point_grows_linearly_with_its_relations(self):
        sizes = node("process.stdout.write(JSON.stringify([0, 5, 10, 15, 20].map(d => G.radius(1, d, 20))));", [])
        steps = [b - a for a, b in zip(sizes, sizes[1:])]
        for step in steps:
            self.assertAlmostEqual(step, steps[0])
        self.assertAlmostEqual(sizes[-1] / sizes[0], 3.0)

    def test_a_gathering_springs_to_its_place_and_stands_still(self):
        moved = node("process.stdout.write(JSON.stringify([0, 0.1, 0.3, 0.6, 1.0, 1.2, 2].map(G.sprung)));", [])
        self.assertEqual(moved[0], 0)
        self.assertEqual(moved[-2:], [1, 1])
        self.assertLess(max(moved), 1.02)
        self.assertGreater(moved[3], 0.9)


class OneThing(unittest.TestCase):
    """The list and the graph read as one thing, written and drawn, as the captain asked on
    3 October 2026: "The two interfaces need to feel like they're doing the same thing, one
    textually the other visually." A name in the list under the pointer or the keyboard lights
    its spacetime in the graph, a spacetime in the graph under the pointer turns its name in the
    list pink and scrolls the list to it, a press on either opens it, and only a spacetime the
    list shows is ever lit. `node _tools/background_graph.mjs` holds the
    page as drawn to the same."""

    @classmethod
    def setUpClass(cls):
        cls.source = PAGE.read_text(encoding="utf-8")
        cls.index = read(build.INDEX_FILE)

    def lit_by(self, cases):
        """The spacetime the graph draws lit, by the page's own inList and lit, for each case of
        the list's answer, the graph's pointer and the name the list lights."""
        script = page_function(self.source, "inList") + page_function(self.source, "lit") + """
        const [index, cases] = input;
        let graph = G.model(index, { edges: [] }), listed, hovered, named;
        const out = cases.map(([found, pointed, fromList]) => {
          listed = null;
          if (found) { listed = {}; found.forEach(id => { listed[id] = true; }); }
          hovered = pointed === null ? -1 : G.placesOf(graph, [pointed])[0];
          named = fromList;
          const i = lit();
          return i < 0 ? null : graph.spacetimes[i].id;
        });
        process.stdout.write(JSON.stringify(out));
        """
        return node(script, [self.index, cases])

    def test_a_name_in_the_list_lights_its_spacetime_and_only_it(self):
        ids = [m["id"] for m in self.index]
        cases = [[None, None, i] for i in ids] + [[None, None, None]]
        self.assertEqual(self.lit_by(cases), ids + [None])

    def test_the_graph_s_own_pointer_goes_before_the_list(self):
        a, b = self.index[0]["id"], self.index[1]["id"]
        self.assertEqual(self.lit_by([[None, a, b], [None, None, b], [None, a, None]]), [a, b, a])

    def test_only_a_spacetime_the_search_shows_is_lit(self):
        found = [m["id"] for m in self.index if any(word_matches(t, "wormhole") for t in m["tags"])]
        out = next(m["id"] for m in self.index if m["id"] not in found)
        lit = self.lit_by([[found, None, found[0]], [found, None, out], [found, out, None], [found, out, found[1]]])
        self.assertEqual(lit, [found[0], None, None, found[1]])

    def test_a_press_or_the_pointer_never_lands_on_a_spacetime_sent_behind(self):
        script = """
        const items = [{ x: 50, y: 50, r: 5, w: 40, h: 16, behind: true }, { x: 200, y: 50, r: 5, w: 40, h: 16, behind: false }];
        const sides = [0, 0];
        process.stdout.write(JSON.stringify([G.hit(items, sides, 50, 50), G.hit(items, sides, 75, 50), G.hit(items, sides, 54, 50),
                                             G.hit(items, sides, 200, 50), G.hit(items, sides, 225, 50)]));
        """
        self.assertEqual(node(script, []), [-1, -1, -1, 1, 1])

    def test_the_names_in_the_list_are_buttons_the_keyboard_reaches(self):
        render = page_function(self.source, "renderResults")
        self.assertIn("'<button type=\"button\" class=\"mfs-result'", render)
        self.assertIn("'</button>'", render)
        self.assertNotIn("<div class=\"mfs-result", self.source)
        self.assertIn(".mfs-result:focus-visible { outline: 2px solid var(--cyan); outline-offset: -2px; }", self.source)

    def test_the_pointer_and_the_keyboard_on_the_list_light_the_graph(self):
        for event in ("pointerover", "pointerleave", "focusin", "focusout"):
            with self.subTest(event):
                at = self.source.index(f"results.addEventListener('{event}'")
                self.assertIn("lightGraph();", self.source[at:at + 300])
        self.assertIn("window._mfsGraph.light(_pointed || _focused);", self.source)
        light = self.source[self.source.index("        light: function(id) {"):][:200]
        self.assertIn("named = id;", light)
        self.assertIn("draw();", light)

    def test_the_graph_turns_the_name_in_the_list_pink_and_nothing_else_and_scrolls_to_it(self):
        # The captain on 3 October 2026: its name turns the page's pink and its button gets no
        # background.
        self.assertIn("window._mfsListLight(i >= 0 ? graph.spacetimes[i].id : null);", page_function(self.source, "hover"))
        lit = [body for selectors, body in re.findall(r"\n    ([^{}\n]*mfs-result-lit[^{}\n]*)\{(.*?)\}", self.source, re.S)]
        self.assertEqual(lit, [" color: var(--pink-light); "])
        self.assertIn("\n    .mfs-result.mfs-result-lit { color: var(--pink-light); }\n", self.source)
        # The list's own hover keeps its background to itself.
        self.assertIn("\n    .mfs-result:hover, .mfs-keyword.mfs-keyword-on {\n", self.source)
        listing = self.source[self.source.index("window._mfsListLight = function(id) {"):][:900]
        self.assertIn("r.classList.toggle('mfs-result-lit', on);", listing)
        self.assertIn("results.scrollBy({ top: by, behavior: reducedMotion() ? 'auto' : 'smooth' });", listing)

    def test_a_press_on_a_spacetime_in_the_graph_opens_it_as_its_name_in_the_list_does(self):
        self.assertIn("if (i >= 0 && window._mfsOpen) window._mfsOpen(graph.spacetimes[i].id);", page_function(self.source, "letGo"))
        self.assertIn("var i = pointedAt(p);", page_function(self.source, "letGo"))
        self.assertIn("return inList(i) ? i : -1;", page_function(self.source, "pointedAt"))
        self.assertIn("el.addEventListener('click', function() { window._mfsOpen(this.dataset.id); });",
                      page_function(self.source, "renderResults"))


class OnThePage(unittest.TestCase):
    """The graph lies behind every panel, across the whole window, stays where it is under an
    open spacetime, and is not drawn on a phone."""

    @classmethod
    def setUpClass(cls):
        cls.source = PAGE.read_text(encoding="utf-8")

    def test_the_graph_lies_under_every_panel_across_the_whole_window(self):
        canvas = self.source.index('<canvas id="mfs-graph"')
        self.assertLess(self.source.index('<canvas id="wavy-grid"'), canvas)
        self.assertLess(canvas, self.source.index('<div id="mfs-search-panel"'))
        rule = re.search(r"#mfs-graph \{(.*?)\}", self.source, re.S).group(1)
        self.assertIn("z-index: 1;", rule)
        for panel in ("mfs-search-panel", "mfs-coffee-panel", "mfs-content-panel"):
            self.assertRegex(self.source, r'<div id="' + panel + r'" style="position:fixed;z-index:15;')
        # It is drawn to every edge of the window, as the captain asked on 3 October 2026, so
        # nothing of it is cut off short of the window's edge.
        self.assertIn("position: fixed; z-index: 1; top: 0; right: 0; bottom: 0; left: 0; width: 100%; height: 100%;", rule)
        # It is gathered in the place the spacetime's panel rests in, to the right of the list.
        room = re.search(r"#mfs-graph-room \{(.*?)\}", self.source, re.S).group(1)
        self.assertIn("top: var(--mfs-top); right: 0;", room)
        self.assertIn("width: calc(100vw - var(--mfs-left-w) - 20px); height: calc(var(--mfs-bottom) - var(--mfs-top));", room)
        self.assertIn("visibility: hidden; pointer-events: none;", room)

    def test_names_run_on_to_the_window_s_edge(self):
        self.assertIn("sides = G.labels(items, W, H, [-at[0], -at[1], CW - at[0], CH - at[1]]);", page_function(self.source, "render"))
        script = """
        const item = { x: 395, y: 100, r: 4, w: 60, h: 16, weight: 1, near: 1, rank: 2, must: false, side: -1 };
        process.stdout.write(JSON.stringify([G.labels([item], 400, 300), G.labels([item], 400, 300, [-500, -150, 600, 300])]));
        """
        inside, beyond = node(script, [])
        self.assertEqual(inside, [1])
        self.assertEqual(beyond, [0])

    def test_a_phone_draws_no_graph(self):
        phone = self.source[self.source.index("@media screen and (max-width: 37.5em)"):]
        self.assertIn("#mfs-graph, #mfs-graph-room { display: none !important; }", phone[:phone.index("</style>")])

    def test_the_graph_stays_under_an_open_spacetime(self):
        # The captain on 3 October 2026: "The graph should not disappear when a spacetime is
        # pressed, it should just get covered by the panel."
        self.assertNotIn("opened", page_function(self.source, "shown"))
        self.assertIn("open: function(on) { opened = on; if (on) hover(-1); },", self.source)
        for hook, on in (("window._mfsShowMetric = function", "true"), ("window._mfsPageBack = function", "true"),
                         ("window._mfsPageAway = function", "false")):
            with self.subTest(hook):
                body = self.source[self.source.index(hook):][:600]
                self.assertIn(f"window._mfsGraph.open({on});", body)
        self.assertIn("window._mfsGraph.leave();", self.source)

    def test_the_page_reads_its_own_build_of_the_graph(self):
        self.assertIn("""<script defer src="{{ '/MFS/assets/graph.js' | relative_url }}?v={{ site.time | date: "%s" }}"></script>""",
                      self.source)
        self.assertIn("""var RELATIONS = '{{ "/MFS/assets/data/relations.json" | relative_url }}?v=' + V;""", self.source)


if __name__ == "__main__":
    unittest.main()
