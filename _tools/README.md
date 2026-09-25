# Site tools

Jekyll ignores any directory whose name starts with an underscore, so nothing in here is published.

## Adding a spacetime to My Favorite Spacetimes

Write the new metric file at `MFS/assets/data/metrics/<id>.json`, with its `id` equal to the file name.
It needs `name` (the long title on the page), `short_name` (the label in the search list) and `tags`.
The search list is ordered by `short_name` alone, case and accents ignored, so Gödel sorts as Godel and de Sitter falls under D.
There is no per entry sort field, and the build refuses a file that carries `sort_name`; a new spacetime takes its place from its name.
Cite references by their key in `assets/data/references.bib`.

Every entry cites each of its references in its `history`, as `[key]` or `[key1, key2]` at the point the prose leans on it, and lists them in `references` in the order they are first cited.
The page turns those brackets into numbered links by looking the key up in `references`, so a square bracket in a history is always read as a citation and never printed.
A reason a reference is there, such as a novel beside the papers, goes in that prose; the `.bib` entries carry no annotations.
`godel` and `morris_thorne` are worked examples.

`_layouts/mfs.html` and `publications.markdown` do not read `references.json`; each parses `references.bib` in the browser with its own small reader.
That reader takes a value nested one brace deep, as in `{Einstein}'s` or `Rebou\c{c}as`, and turns the accent commands `\"`, `\'`, `` \` ``, `\^`, `\~`, `\c` and `\ss` and the escape `\&` into characters, and nothing else.
A value outside that set prints wrongly on the page without any error, so check a new entry's rendered line and not only the build.

Then run one command from the top of the repository:

    python3 _tools/build_mfs_data.py

That rewrites `MFS/assets/data/metrics_index.json` and `MFS/assets/data/references.json` from what is on disk.
The index is never edited by hand, so it cannot disagree with the folder.
Commit the new metric file together with whatever the command rewrote.

Every key a metric cites has to be an entry in `assets/data/references.bib`.
The command refuses to write anything if one is not, naming the metric file and the key it could not find, so a mistyped citation is caught here rather than published as a reference the reader cannot resolve.

Write a component value so that putting a minus sign in front of the whole string negates it.
The page groups a tensor's components by value and merges a value with its negation, printing `R^t{}_{\theta t\theta} = -R^t{}_{\theta\theta t} = \dots` on one line, and it recognises the negation by that leading minus alone.
So a value that is a bare sum wants collecting over a common denominator or wrapping in `\left(\right)` first, and its opposite wants writing as that string with a `-` in front rather than with the signs distributed through it.
Morris-Thorne is the worked example: its twenty four Riemann components print on six lines fully lowered and its twenty four Weyl components on six, and written the other way each would take twenty four.

A line too wide for the page scrolls sideways on its own, and nothing else on the page moves.
A tensor component keeps its first index label in place while the rest of the line scrolls, as the application does; any other formula scrolls whole.
The print copy cannot scroll, so there a wide line breaks before a top level `+`, `-` or any `=` but its first, and a line holding one piece too wide to break, such as a large fraction or a bracketed sum, is set small enough to fit the paper.
So a value whose whole length sits inside one `\left(\right)` prints small, which is one more reason to collect it rather than expand it.
`mathLine` and `printPrepare` in `_layouts/mfs.html` carry the details.

Nothing in the print copy may follow the last line of its body, not even blank space.
The repeated footer sits below the body's whole box, and closing space that pushes it past the foot of a page makes Chrome print one more page holding only the header and footer.
So printed sections are spaced from the section before them and never by a margin underneath, and the print stylesheet's comment beside that rule says why.

A formula inside prose, in a `history`, a `convention` or a caption, is set inline and cannot break, so one wider than its paragraph is given a line of its own that scrolls in the same way, and the prose around it wraps as before.
`fitProseMath` in `_layouts/mfs.html` measures that again whenever the paragraph changes width.

## The page on a phone

A screen narrower than 600px, or a touch screen under 500px tall, gets the same panels in one column that the page scrolls through: the title, the list with the coffee panel, then the spacetime.
That covers every iPhone upright and on its side, and the desktop layout is untouched by it.
The rules are one `@media screen` block in `_layouts/mfs.html`, after the main stylesheet, and its comment says why each choice was made.
They are for the screen alone, so the print copy is the same whatever screen it was printed from; a phone printing Kerr differs from a desktop only by MathJax's rounding at three pixels to the point, and did before this layout.

The desktop's panel scripts write their geometry into each panel's own style, so the phone block overrides it with `!important`, and it sets `--mfs-phone` on the root, which is how the scripts tell which layout is in force.
Anything new on the page has to hold at 390pt upright and 844 by 390 on its side: a line of mathematics, a coordinate domain or a formula in prose scrolls on its own there, and nothing may make the page wider than the screen, even for a frame, since a phone answers that by zooming the whole page out.
A spacetime diagram is drawn whole, as wide as the panel and never taller than the screen, with its labels scaling with it.

The page prints a spacetime's `signature` and `convention` under the heading "conventions", between the history and the coordinates, which is where the application reads them.
They belong to the spacetime and not to a chart, so the coordinate selector leaves them standing.
A `convention` is prose with inline TeX between dollar signs, split into paragraphs at `¶`, and never carries a citation or HTML.
An entry with neither field shows no conventions section at all.

A cloud-sync conflict copy dropped into the metrics folder, named like `kerr 2.json`, is passed over rather than read.
Those are the same names `.gitignore` already keeps out of the repository.
Any other file name is read, so a metric whose `id` does not match its file name still stops the command.

## Dashes in the prose

The prose fields carry no dashes as punctuation.
That covers `name`, `short_name`, `description`, `history`, `convention` and every other field written in sentences, including prose set between dollar signs.
Where a sentence wants a break, write a full stop, a semicolon, a comma or a pair of commas; an em dash, an en dash and a plain dash are all out.
If a sentence only reads with the break in it, rewrite the sentence rather than leave a stump.

A hyphen stays only where it is part of a spelling.
That means a pair of names, such as Reissner-Nordstrom, Kerr-Newman or Lanczos-van Stockum, and an established term, such as anti-de Sitter, pp-wave or Taub-NUT.
An ordinary English compound is rewritten so that the hyphen is not needed.
Minus signs inside the LaTeX fields are mathematics and are left alone.

The tests hold this rule over every metric file on disk, so a new spacetime carrying one of these dashes in its prose fails before it is published, named along with the field and the character that tripped it.

A character outside ASCII is written as itself, as the `î` of `Lemaître` is, and never as a JSON `\u` escape with its backslash doubled, since neither the page nor TeX reads such an escape and the reader sees all six characters of it.
The tests hold that rule over the same fields.

## Whom the prose is for

The prose on the page is for someone who came for a spacetime, and never for whoever builds the collection.
So it never names the collection's own machinery: no entry, no published block or published chart, no field in the JSON sense, no file, no variant, nothing said to be printed or published, and no metric identifier in backticks.
And no sentence has the page for its subject: nothing is above or below, no section shows anything, and no history describes anything.
A sentence whose subject is the collection, an entry, a chart as the collection lists it, a block, a field, a section or the page is rewritten so that its subject is the coordinate, the surface, the parameter, the horizon or the claim, as "Components are taken in the chart $x^0 = ct$" and "No component assumes a field equation" are.
A coordinate chart in the geometric sense stays, since it is mathematics, and so does a field that is physics.
The derivation notes follow the same rule, and where one has to say what the page states, it names the thing itself: the metric components, the Christoffel symbols or the Weyl tensor.
This file is written for whoever builds the collection, and is the one place its words belong.

The tests hold the words that can only ever mean the machinery out of every prose field of every metric and out of every caption, label and declared input of every diagram.
A history keeps "published" and "printed" for the papers it tells of, and "building block" stays wherever it is written, since neither names the collection.

## The shape of a history

A history has at least five paragraphs, and every paragraph has three to six sentences, so no paragraph runs to more than twice the length of another.
A table, the paragraph written as a `TABLE::` line, is not prose and is left out of the count.
Where an entry is too thin for five paragraphs, it wants more history, sourced and cited like the rest, never filler.

`sentences` in `build_mfs_data.py` does the counting.
A sentence ends at a full stop, a question mark or an exclamation mark, after any closing quote or bracket, where the next word begins with a capital or a digit.
A capital standing alone before a full stop is an initial, as in J. Robert Oppenheimer, and ends nothing.
Mathematics between dollar signs counts as one word, and ends a sentence only when the stop is inside it, as when a displayed equation closes one.

The command refuses to write anything, and `--check` fails, while any history is out of shape, naming every such history with the sentences in each of its paragraphs and the paragraph that breaks the rule.

## What the command writes

`MFS/assets/data/metrics_index.json` is one entry per metric file, in the order the search list shows them, carrying `id`, `name`, `tags` and `version`.
Its `name` is the metric's `short_name`, which is what the search list shows, not the long `name` the page titles the spacetime with.
An entry also carries `diagrams` when the spacetime has a file in `MFS/assets/data/diagrams/`, and `conformal` when it has one in `MFS/assets/data/conformal/`; the page fetches either file only for an entry that says so.

`MFS/assets/data/references.json` is the whole of `assets/data/references.bib` parsed into JSON, so a reader can pull every reference the collection cites in one fetch instead of walking every metric file.
The `.bib` file stays where it is and stays the one place a reference is written; the publications page still reads it.

## The version stamps

A metric's `version` is a hash of that metric's content.
It changes when the content changes and not otherwise, so reformatting a file or republishing the site leaves it alone.
This is what the iOS application uses to fetch only the spacetimes that actually changed, since GitHub Pages rewrites every file's own tag on every publish.
`references.json` carries one stamp of its own at the top for the same reason.
An index entry's `diagrams` and `conformal` are the same kind of hash over that spacetime's diagram file and conformal diagram file, so the application can fetch a redrawn diagram without fetching its metric again.

## Checking

    python3 _tools/build_mfs_data.py --check

reports whether the published files are still what the folder says they should be, and changes nothing.

    python3 -m unittest discover -s _tools

runs the tests, which include that check.

## Spacetime diagrams

`MFS/assets/data/diagrams/<metric_id>.json` holds the spacetime diagrams of one spacetime: null rays and future light cones on a plane of the time and one spatial coordinate, for each coordinate system that has any.
The page draws them in a section per coordinate system, after the geodesics, and the application displays the same files.
Nothing in them is drawn by eye or computed by the page.
`_tools/derivations/null_rays.py` computes every ray, cone and marker from the system's published `metric_components`, `inverse_metric_components`, `kretschmann` and `domains`, read through the checker's own `Reader`, and its docstring records the method.

A domain that holds only for some values of a parameter ends in `\;\text{for}\;` and a condition, as FRW's closed case writes `r \in [0, 1/\sqrt{k}) \;\text{for}\; k > 0` beside `r \in [0, \infty) \;\text{for}\; k \le 0`.
The script hatches a view by the domains whose conditions hold at that view's parameter values, and stops, naming the view, on a condition it cannot evaluate or on two domains for one coordinate that both hold.
A domain with no `\in`, such as a note on where a horizon sits, is printed and never read.

It needs more than sympy, and takes about nine minutes for the whole collection, the longest single view being FRW drawn through its observer at about a minute and a half:

    python3 -m venv /tmp/mfs-venv && /tmp/mfs-venv/bin/pip install sympy numpy scipy contourpy
    /tmp/mfs-venv/bin/python _tools/derivations/null_rays.py
    python3 _tools/build_mfs_data.py

`--metric <metric_id>` redraws one spacetime and is repeatable, which is what to use after editing one entry.
The second command stamps each diagram file's version into the index, as it does for the metrics.

`DIAGRAMS` in `null_rays.py` is the table of every view, one row each: the plane, the parameter values, the coordinates held fixed, the plot range, the orientation rule and any declared input.
`CAPTIONS` beside it carries each view's caption, which is prose, and the tests hold the captions, the labels and the declared inputs to the rules for prose as they hold the metrics.
A new view is a row in each, and the script refuses to run while one lacks the other.

A caption opens by naming its plane: the two coordinates drawn and the value of every coordinate held fixed, as "the plane of $t$ and $r$ at $\theta = \pi/2$ and $\phi = 0$".
It says what the drawing shows, and where the feature a reader comes looking for lies off the plane, it says where that feature is, as the Ellis-Bronnikov caption places the throat in $g_{\theta\theta}$ and the Gödel caption places the closed timelike curves.
A caption never stops at saying what a diagram leaves out.
Every curve drawn is a null curve; a caption calls it a null geodesic, the path light takes, only where no Christoffel symbol turns it out of the plane, and says so where one does, as for Gödel and for Kerr off the axis.

### Labels are TeX

Every label in a view is text with its mathematics in `$...$`, the form a caption takes: the button's name, the two axes, each tick and the line where a dust solution starts.
An axis names its coordinate and its unit as the entry's line element spells them, `$ct/r_s$` and not `ct / r_s`, so the page sets it with the same MathJax as the metric above it and the application with SwiftMath.
The script picks each axis's ticks, a step of 1, 2, 2.5 or 5 times a power of ten giving about six across, and writes them into the view as `ticks` with the value `at` in the box's units and the label as TeX, so neither reader chooses its own.
The page lays the labels over the SVG as HTML rather than drawing them as SVG text, and `_layouts/mfs.html` says why beside the renderer.

### A diagram is tied to what it was drawn from

Every view records the fields it was drawn from and a stamp over them, computed by `diagram_source_version` in `build_mfs_data.py`, the one function both scripts use.
The fields are the coordinates, the parameter symbols, the metric and its inverse, the Kretschmann scalar and the domains, and the Einstein tensor as well for a view whose input is solved as dust.
`build_mfs_data.py` recomputes each stamp from the metric file as it stands and refuses, naming the file, the system and the view, when one no longer matches, in `--check` and when writing alike.
So an edit to any of those fields leaves the collection unpublishable until its diagrams are redrawn, while an edit to a history, a reference or a parameter's description leaves them standing.

### Which way the cones point

Every cone is a future cone, and each row names how its chart is oriented.
`tau` takes a time function, the chart's $t$ unless the row says otherwise.
A static chart's $t$ stops being a time beyond a horizon, and there the chart alone cannot say which way is future, so those rows take a family instead.
`ingoing` makes the ingoing rays future directed toward smaller $r$ everywhere, which agrees with $t$ outside and reads the region inside $r_s$ as the black hole, as the ingoing Eddington-Finkelstein chart does; Schwarzschild's spherical chart, Reissner-Nordstrom, Taub-NUT, Kerr and Kerr-Newman use it.
`outgoing` does the same with the outgoing rays toward larger $r$, for de Sitter's static chart and the Krasnikov tube.
`vector` takes the chart's time direction, for Godel, whose published $g^{tt}$ is positive and which has no time function at all.

### Declared inputs

An entry that leaves a function free cannot be drawn without a choice, and every such choice is written in its row and printed beside its diagram.
FRW's scale factor and Bianchi I's three are dust, solved from each entry's own published $G^i{}_i = 0$.
Morris-Thorne is drawn as its Ellis-Bronnikov member, $\Phi = 0$ and $b = b_0^2/r$.
Alcubierre is drawn with his own profile, which the entry names, at $v_s = 2$; Natario with the same profile and his zero expansion field, whose divergence the script checks before using it; the Krasnikov tube with a tube built by a ship at the speed of light.
Kasner is drawn at the exponents $(-2/7, 3/7, 6/7)$.

### Checking the rays

    /tmp/mfs-venv/bin/python _tools/derivations/null_rays.py --verify

traces rays as the page's files are traced and measures, for every view with a closed form, how far the quantity each family should conserve drifts along a ray, together with the dust solutions and the equality of Natario's and Alcubierre's metrics on the plane of their axis.
It writes nothing and exits non-zero if anything fails; run it after changing the method.

### What is not drawn

Kerr and Kerr-Newman are drawn on the axis only.
Off it the null curves of fixed $\theta$ and $\phi$ are not null geodesics, since $\Gamma^\theta{}_{tt}$, $\Gamma^\theta{}_{rr}$ and $\Gamma^\phi{}_{tr}$ turn a light ray launched along one out of the plane, and inside the ergoregion the plane has no null direction at all.
The radial null geodesics they do have are the principal null congruence, which leaves the plane, and every null geodesic of van Stockum leaves its plane too.
Tolman-Bondi, Vaidya and the proper distance chart of Morris-Thorne leave functions free that no choice has been made for yet, and Oppenheimer-Snyder's collapse needs its two charts drawn together.
The Tolman-Oppenheimer-Volkoff star leaves its redshift and mass functions free, and its one closed form member, the star of uniform density, is drawn already as the interior Schwarzschild solution.
Lentz's soliton exists only as a numerical integral over his rhomboid sources, so no potential of his can be declared in closed form.
The Mixmaster's scale factors have no closed form to declare, being the chaotic solutions of its own field equations, and the one plane of time and an Euler angle that keeps its light rays, time against $\psi$, would show only $a_3$.
The Malament-Hogarth toy is Minkowski space times a conformal factor, which changes no null direction, so its diagram would be Minkowski's with one point missing; what it changes is proper time, which a diagram of light rays cannot show.
The cosmic string's two charts, Oppenheimer-Snyder's exterior, Natario's plane flow chart, the Brinkmann chart of the pp-wave and Minkowski's double null chart would each only repeat a plane drawn elsewhere, flat or the same as another chart's.

## Conformal diagrams

`MFS/assets/data/conformal/<metric_id>.json` holds the conformal diagram of one spacetime: the whole spacetime brought to a finite drawing with light at 45°, or, where no picture of the whole is faithful, a totally geodesic surface in it that says so.
The page draws it once per spacetime, under the heading "conformal diagram", between the conventions and the coordinates, since it belongs to the spacetime and not to a chart; its views are buttons, and the view shown first is the one that tints the region the selected coordinate system covers.
The application reads the same files.

`_tools/derivations/conformal.py` draws them, one function per spacetime, each carrying its derivation in its docstring, reading the published metric through the checker's `Reader` with the `load` and `published_matrix` the null rays use.
It needs the same environment as `null_rays.py` and takes a few seconds for the whole collection:

    /tmp/mfs-venv/bin/python _tools/derivations/conformal.py
    python3 _tools/build_mfs_data.py

`--metric <metric_id>` redraws one spacetime, and `--verify` prints every check and writes nothing.

Every map drawn is checked by the function that draws with it, at random points of the region it covers, against the published metric and inverse metric: lines of constant drawn null coordinate are light rays, the two families are distinct, the future is up, and the published inverse on the surface is the inverse of the published metric there.
Each spacetime adds the limits that place its horizons, infinities and singularities, and the published Kretschmann scalar must diverge wherever a line is drawn as a singularity and stay finite on every centre or throat drawn as regular.
The script refuses to write anything while one check fails; the docstring lists them.

A file records every coordinate system it read and a stamp over the fields it read from each, in `source`, computed by the same `diagram_source_version` as the null rays' stamps.
A diagram can read another spacetime's metric, as the interior Schwarzschild star reads the exterior of `schwarzschild.json`, and `build_mfs_data.py` refuses the file when any of those fields changes, naming the system and the command that redraws it.
The domains are not among the fields, since no conformal diagram reads one.

`DRAWN` in the script names the seventeen spacetimes that have a diagram, and `NOT_DRAWN` gives each of the others the reason it has none: every event's future is the whole spacetime for Gödel and van Stockum, the diagram is whatever a free function makes it for the warp drives, the Krasnikov tube, Tolman-Bondi and Bianchi, the Mixmaster has no surface that carries its causal structure, and so on.
A free function does not by itself rule a diagram out: where every choice of it gives the same shape, as for the Tolman-Oppenheimer-Volkoff star and the Morris-Thorne wormhole, the diagram is drawn with a declared choice that moves only the lines inside, and the Malament-Hogarth toy is drawn for every conformal factor at once, since a conformal factor changes no null direction.
Nothing is published for those, and the script stops if a metric file is in neither table, so a new spacetime needs a decision.

A view of a surface that is not the whole spacetime carries `restriction`, which the page prints in a band across the top of the figure, never in a footnote.
Kerr and Kerr-Newman are drawn on the symmetry axis and the cosmic string on the half plane of fixed $\phi$ and $z$, and the tests hold those three to carrying a restriction on every view.

Every text in a view is TeX in `$...$`: the labels on the drawing, the buttons, the legend, the caption, the restriction and the parameter values.
A caption opens by naming what is drawn, "This is the whole of ..." or "This is the plane ...", and its prose is for the reader of "Whom the prose is for" above; the tests hold every text in these files to that rule's words and to the dash rule.

`_tools/derivations/tex_check.cjs` typesets every TeX string the page sets, in the metrics, the diagram files and the conformal files, a published value together with its negation, through the TeX input MathJax loads on the page, and exits non-zero naming each one it cannot set.
sympy never reads the typesetting, so this is the check that catches a value that is right and prints as an error box:

    npm install --prefix /tmp/mfs-node mathjax-full
    node _tools/derivations/tex_check.cjs /tmp/mfs-node

## Checking the physics

`build_mfs_data.py` checks the shape of the data. It does not read the mathematics.
`_tools/derivations/verify_metrics.py` does, for every coordinate system of every metric file.

It builds $g_{\mu\nu}$ from the line element each entry prints, computes the inverse metric, both Christoffel variants, Riemann, Ricci, the Ricci scalar, Kretschmann, Einstein and Weyl in sympy, and compares all of it against every value the entry publishes.
A component the entry omits has to vanish, so a missing symbol is caught as well as a wrong one.
It names each disagreement by file, system and symbol, and exits non-zero if any remain.

It also checks that every term of every published expression carries the dimensions of its left hand side, which needs no algebra at all.
That pass names the file, the system, the field and the term it could not balance, and it is what catches a component printed in the bare chart rather than in $x^0 = cT$.

sympy is not installed system wide, and the virtual environment does not belong in the repository:

    python3 -m venv /tmp/mfs-venv && /tmp/mfs-venv/bin/pip install sympy
    /tmp/mfs-venv/bin/python _tools/derivations/verify_metrics.py

The whole collection took 259 seconds on 25 September 2026, and `--system <metric_id>/<system_id>` checks one system in seconds, which is what to use while editing a single entry.
The slowest systems are the general flow chart of `natario` and the potential flow of `lentz`, each under a minute, then `mixmaster` at about forty seconds; Kerr takes about seven seconds and the interior Schwarzschild solution about two.
`--budget <seconds>` changes how long sympy may spend on one tensor, and the default of 120 is several times what any tensor in the collection needs.
`--dimensions-only` runs the dimensional pass alone, which takes about a second over the whole collection, so there is no reason not to run it on every edit.

The speed is `norm`, the routine every tensor passes through, which puts an expression into a canonical form so that one that vanishes is exactly zero.
It reads the expression as a rational function of generators, reduces $\cos^2$ by $1-\sin^2$, and keeps every denominator as a product of irreducible factors, so it never takes a polynomial gcd; its docstring records what else was measured and why it lost.
Until 22 September 2026 it called sympy's general `simplify` instead, which took the whole collection about forty minutes and never finished Kerr's Riemann tensor at all.
It splits a radical into the square roots of the irreducible factors of its radicand, which is an identity only where each factor is positive.
Every radicand in the collection is positive factor by factor on the region its entry describes, and that reading is what lets the interior Schwarzschild radicals cancel.
An entry whose radicals change sign inside the region it claims would need that looked at again.

Both passes exit clean on the collection as it stands, and neither leaves anything `UNCHECKED`: since 22 September 2026 every published expression of every system is compared.
The last to be reached were the seventy Tolman-Bondi expressions that name a third derivative, which the reader now declares, and the inverse metric of Ellis-Bronnikov, which the entry now publishes.
`UNCHECKED` does not set the exit code, so a new one shows only in the output, and a clean run is one that prints none.
Every system publishes its inverse metric, because the page prints it beside the metric, so a system without one is reported `UNCHECKED` as missing it rather than excused.
Both passes print what failed on stderr and print the single line saying nothing failed on stdout, so a run piped through `2>&1 | tail` can look clean when it is not.
Read the exit code, and run the script before and after a change so that a failure you did not cause is not mistaken for one you did.

`derivations/audit-2026-09-18.md` groups and counts the disagreements the collection carried when the dimensional pass was added, and says which of them are the checker's fault rather than the physics'.
It is a dated snapshot of 433: the 320 that were Weyl components copied from the entry's own Riemann components have since been corrected, the 50 nested radicals cancel since the change to `norm`, and the rest were corrected entry by entry, leaving none.
`derivations/weyl.md` is the working behind those corrections, and is the thing to read before touching any `weyl_tensor`: a Weyl tensor equals Riemann only in a vacuum, it can be nonzero in a slot where Riemann vanishes, and two of the entries that publish one are conformally flat and so publish nothing.

Nothing is ever passed in silence. A value that cannot be parsed, a system with no time coordinate declaration, or a tensor sympy cannot finish in the budget is reported as `UNCHECKED` with the reason, separately from the disagreements.

## Printing a chart by machine

`tov`, `malament_hogarth`, `mixmaster` and `lentz` have their mathematics written by `_tools/derivations/print_charts.py`, which defines each of their charts, computes every tensor with the checker's own `Geometry`, and prints each value through `chart_printer.py` beside it:

    /tmp/mfs-venv/bin/python _tools/derivations/print_charts.py --metric tov

Every printed value is read back through the checker's `Reader` and compared with the value it came from before anything is written, and `verify_metrics.py` then checks the file like any other.
The script writes only the coordinate system's mathematics, the line element and the domains; the prose, the convention and the description of each parameter stay in the metric file and are carried over, and a parameter with no description stops the script.
A chart can pass a function that regroups the numerator of each value, which is how TOV's curvature is printed around $(\partial_r\Phi)^2 + \partial_r^2\Phi$ and Bianchi IX's around $a_1^2\cos^2\psi + a_2^2\sin^2\psi$, and a Ricci or Kretschmann scalar can be given in a structured form, which is checked against sympy before it is used.
`tov.md`, `malament_hogarth.md`, `mixmaster.md` and `lentz.md` in the same folder record why each chart is the one published.
The other entries were written by hand and are not touched by the script.

## The three conventions the checker encodes

All three are things the files do consistently, and only the newest of them say so in their own `convention` field, so they are recorded here and in the script's header as well.

The chart is $x^0 = cT$.
A coordinate carrying dimensions of time is not itself the chart coordinate; $c$ times it is, and every published component is a component in that chart even though the index is printed with the bare name.
Schwarzschild shows it plainly, publishing $g_{tt} = -(1-r_s/r)$ against a line element whose time term is $-(1-r_s/r)c^2dt^2$.
Because the rescaling is linear, a component in the chart is the one taken with the bare coordinate multiplied by $c$ once per upper time index and divided by $c$ once per lower one.

The Ricci tensor is the standard contraction $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$, on the first lower index, and it carries through to the Einstein tensor and the Ricci scalar.
That is the contraction the signature $(-,+,+,+)$ asks for, settled by the captain on 18 September 2026: the Ricci tensor should be whatever the metric sign convention has.
It is what makes ordinary matter come out with a positive energy density, so FRW publishes $G_{tt} = 3(\dot a^2+k)/a^2$ and the interior Schwarzschild solution $G^t{}_t = -3r_s/R^3 = -8\pi G\rho$.
The collection contracted on the last index before that, which is minus this in every slot, and every entry with a nonzero Ricci tensor was moved over in one pass.
The Weyl tensor is built from the same contraction and always was, because it is defined by removing the traces of Riemann and those traces are Riemann's own.

The dots in a geodesic equation are velocities of that same chart, so a dot on a time coordinate means $d(cT)/d\lambda$ and not $dT/d\lambda$, even though it is printed on the bare letter.
It has to be that one, because the equation is $\ddot x^\mu + \Gamma^\mu{}_{\nu\rho}\dot x^\nu\dot x^\rho = 0$ with the same printed $\Gamma$ the entry lists, and those are chart symbols.
The reading is checkable without sympy: on it every term of every equation carries the dimensions of its left hand side, and on the other reading a term mixing a time velocity with a space velocity comes out wrong by one factor of $c$.
That is what the dimensional pass does, over the geodesics and over every published component alike, and `_tools/derivations/kasner.md` Step 17 works the check through term by term for one entry.

## Adding a spacetime, as far as the checker is concerned

Telling a time coordinate from a length needs the dimensions of the parameters, which live in prose, so the script cannot read it off the file.
`DIMENSIONS` in `verify_metrics.py` declares the dimension of every coordinate and every parameter, per system.
A coordinate declared as a time is exactly a coordinate the chart multiplies by $c$, so that one table answers both which chart a component is in and what each term of it has to carry; there is no second table to keep in step with it.
A new coordinate system that is not listed there is reported `UNCHECKED` rather than guessed at, so adding a spacetime means adding its line, and forgetting to is visible rather than silent.
So is declaring a dimension for a parameter the system does not have, or leaving one out.

A dimension is written in $L$, $T$ and $M$.
Most entries never name a mass, because they fold it into a length such as $r_s = 2GM/c^2$ and quote that; Vaidya keeps $G$ and $m(u)$ explicit, so it declares $[G] = L^3M^{-1}T^{-2}$ and $[m] = M$ and its components balance through the length $Gm/c^2$.
Kerr does the same with a constant $M$, and adds the one declaration a reader cannot guess: its spin parameter $a = J/(Mc)$ is a length, not a dimensionless number, which is what makes $r^2 + a^2\cos^2\theta$ an area.

Some of what the table has to decide is a choice the entry leaves open.
FRW can be read with a dimensionless comoving $r$ and a scale factor carrying the length, or with $r$ a length, $a$ dimensionless and $k$ a curvature; the line element balances either way, and only the published Riemann tensor picks the second.
Where that happens, the table carries a comment saying which reading the entry's own values obey.

A parameter may be a function rather than a constant.
`a = a(t)` declares a function of one coordinate, which answers to a dot and to a prime, and `H = H(u,x,y)` declares a function of several, which answers to `\partial`: a published `\partial_x H`, `\partial_x^2 H` or `\partial_x\partial_y H` is read as that partial derivative, and so is one of any higher order.
Every such derivative is taken with respect to the chart coordinate, exactly as the dot and the prime are, so one along a time carries a factor of $1/c$ and one along a length does not.
`_tools/derivations/pp_wave.md` is the worked example, and its Step 15 is where that factor earns its keep.

The reader fixes no highest order.
A curvature is two derivatives of the metric, so it reaches the second derivative of a function only while the line element carries that function undifferentiated, and Tolman-Bondi is the entry where it does not: its $g_{rr}$ is built from $\partial_r R$, so its curvature reaches $\partial_r\partial_t^2R$.
`Reader._declare_partials` therefore declares a partial derivative when a published value names one, at whatever order it is written and in whatever order its factors are spelled, since mixed partials commute.
It still refuses a derivative of anything the system does not declare as a function, or along a coordinate the function is not declared to depend on, so a typo is an error naming the function and the coordinate rather than a new symbol.
`_tools/derivations/tolman_bondi.md` Step 18 is the worked example, and records that the entry's seventy third derivative expressions came back `UNCHECKED` until 22 September 2026.

The prime works only on a function whose name is not a LaTeX command.
The reader turns a prime into a name suffix before it turns `\Phi` into a name, so `b'` reads as the derivative of `b` while `\Phi'` leaves a stray backslash and is rejected as unhandled LaTeX.
`\partial` works on either kind of name, so an entry with a Greek named function writes every derivative that way; `_tools/derivations/morris_thorne.md` Step 17 says so of its own `\Phi` and `b`.

An entry whose parameters are not free needs a second line, in `PARAMETER_RELATIONS`.
Kasner prints three exponents bound by $\sum_i p_i = \sum_i p_i^2 = 1$ and claims its values only on the surface those two equations cut out; its Ricci tensor is zero there and nowhere else.
Checked against free symbols it would report seventy four disagreements that are not disagreements.
The table carries a rational parametrisation of the surface instead, and every constrained parameter is replaced by its parametrised value on both sides of each comparison, so the identity checked is the one the entry asserts.
That is exact rather than a sample, because a rational parametrisation covers a dense subset of an irreducible variety.
`_tools/derivations/kasner.md` derives the parametrisation the Kasner entry uses and why the constraints are the field equations.
