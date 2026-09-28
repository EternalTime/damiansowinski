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
Each of these observers measures everything it was handed before it changes any class, because a class changed between two measurements makes the second lay the whole page out again, and once per line that held Kerr-Newman still for seconds and Natario for nearly a minute.

A chosen spacetime's panel waits off the screen until its mathematics is set and laid out, and only then slides in.
While it waits, "Spacetime data loading..." shows where the panel rests, centred, under the panel so the panel slides in over it.
Most of that wait holds the main thread, so the line cannot be shown once it starts; it is asked for as the wait begins and fades in by a CSS transition with a delay, which the browser runs off the main thread, so a spacetime set within that delay never shows it.
`showNote` in `_layouts/mfs.html` carries the timings and why.

## The page on a phone

A screen narrower than 600px, or a touch screen under 500px tall, gets the same panels in one column that the page scrolls through: the title, the list with the coffee panel, then the spacetime.
That covers every iPhone upright and on its side, and the desktop layout is untouched by it.
The rules are one `@media screen` block in `_layouts/mfs.html`, after the main stylesheet, and its comment says why each choice was made.
They are for the screen alone, so the print copy is the same whatever screen it was printed from; a phone printing Kerr differs from a desktop only by MathJax's rounding at three pixels to the point, and did before this layout.

The desktop's panel scripts write their geometry into each panel's own style, so the phone block overrides it with `!important`, and it sets `--mfs-phone` on the root, which is how the scripts tell which layout is in force.
Anything new on the page has to hold at 390pt upright and 844 by 390 on its side: a line of mathematics, a coordinate domain or a formula in prose scrolls on its own there, and nothing may make the page wider than the screen, even for a frame, since a phone answers that by zooming the whole page out.
A spacetime diagram is drawn whole, as wide as the panel and never taller than the screen.

## The reader's text size

Every size a reader reads on the spacetimes page is in rem, so a reader who enlarges text in the browser, by Chrome's font size, a larger default font, or Safari's and Firefox's zoom of text only, enlarges every word on it: the name, the headings, the buttons, the domains, the prose, the mathematics, the list and every word on both kinds of diagram.
At the browser's usual 16px each is the pixel size it was when the page was set in pixels, and only the hairlines of borders stay in px.
The list's column grows with the text up to a quarter of the window, the panels start below the page's title however it wraps, the one-column layout's widths are in em so a large text size takes it on a narrower desktop window, and a sticky header that would cover more than a quarter of the view scrolls away instead.
A spacetime diagram's numbers and axis names are laid out around its plot, so they take the room they need and the numbers thin out where they would touch.
A conformal diagram's labels and a figure's stand at points inside the drawing, so the drawing keeps its proportion to them and is drawn larger with the text, scrolling sideways in its own frame where it is wider than the panel.
`nrFigure` and `cdFrameAround` in `_layouts/mfs.html` carry the details.
No word is broken inside itself at any size: the name, the headings, the choices, the page's title, the names in the list and the search hint are held to the size at which their widest word fits their line, which `fitWords` and `mfsWidestWord` measure.
A word of the prose wider than its whole line is hyphenated by `hyphenateWords` with TeX's English patterns, since the browser's own hyphenation never divides a capitalised name such as Schwarzschild.
On a phone the page's title and its exit sign share a row until they no longer fit side by side, when the sign takes a row of its own below the title.
The print copy's sizes are points and its prose 12pt, so it prints the same whatever the reader's text size on the screen.
Reproduce a text size fault the way a reader meets it, with a larger default font size and not with page zoom, which scales everything and hides it.

Every choice on the page, the chart, a view of a diagram and a tensor's index placement, is one control, a row of buttons built by `choiceButton`, and a row appears only where there is more than one thing to choose.
Every drawing, a spacetime diagram, a conformal diagram and a figure in three dimensions, stands on one dark ground, `--mfs-ground`, across the whole figure with its words, and prints on white.

The coordinates come first in the mathematics, right after the history, since the chart is chosen before anything that depends on it: their buttons, then each coordinate with its domain.
The page prints a spacetime's `signature` and `convention` under the heading "conventions" just below them, which is where the application reads them.
They belong to the spacetime and not to a chart, so choosing a chart leaves them standing.
A `convention` is prose with inline TeX between dollar signs, split into paragraphs at `¶`, and never carries a citation or HTML.
An entry with neither field shows no conventions section at all.

A cloud-sync conflict copy dropped into the metrics folder, named like `kerr 2.json`, is passed over rather than read.
Those are the same names `.gitignore` already keeps out of the repository.
Any other file name is read, so a metric whose `id` does not match its file name still stops the command.

## How long a spacetime takes to open

`_tools/page_timing.mjs` opens every spacetime and every chart in headless Chrome as a reader does, and times each from the click until its mathematics is typeset and the page answers again:

    bundle exec jekyll serve
    node _tools/page_timing.mjs http://127.0.0.1:4000

It prints each chart's time and the longest task that held the page, slowest last, and exits non-zero when a chart takes longer than `--budget`, 10 seconds by default.
`--metric <id>` times one spacetime, `--phone` lays the page out as an upright iPhone, and `--cpu 4` slows the processor four times over, as Chrome's own tools do to stand in for a slower phone.
It needs Chrome and Node 22 or later, and nothing installed.

Once the observers measured every line before marking any, what was left of a large chart's time went to laying out a whole tensor again and again.
A tensor with room for one column only was a multicol container of one column, which the browser lays out whole whenever anything inside it changes, so each index toggle, each fade marked at a line's edge and each change in the width of the prose above laid out all of Natário's Weyl tensor once more.
`reflowColumns` in `_layouts/mfs.html` now leaves a tensor of one column a plain block, as the print stylesheet already did, and counts every tensor's columns before it sets any.
`printPrepare` measures every line of the print copy before it sets the width of any, for the same reason the observers do.

On 25 September 2026, in Chrome 154 on an M4 Pro at 1440 by 900, Natário's general flow chart was ready in 2.9 seconds, 2.0 of them in MathJax, Lentz in 1.9, the Mixmaster in 1.6 and Kerr-Newman in 1.0, where with the multicol containers they had taken 4.7, 2.9, 2.4 and 1.4, and every other chart was ready in under a second.
Laid out as an iPhone with the processor slowed four times, Natário's took 10.0 seconds, 7.9 of them in MathJax, Lentz 6.5 and the Mixmaster 5.1, where they had taken 18.5, 13.7 and 10.0, and every other chart took under 5.
Once that chart was open, flipping the first index of its Weyl tensor took 0.5 or 1.1 seconds, where it had taken 1.5 or 2.9, and a change in the window's width 0.1, where it had taken 2.8.
Pressing print on it brought up the print dialog in 2.1 seconds, where it had taken about 3.
Those three are the largest by far, with a quarter of a million, 150,000 and 125,000 elements on the page, and the next, Kerr-Newman and Alcubierre, have about 60,000.
With every layout done once, the time left grows with the length of the published expressions and is mostly MathJax setting them, so a chart several times the size of Natário's would want its tensors typeset as the reader reaches them rather than all at once.

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
An entry also carries `diagrams` when the spacetime has a file in `MFS/assets/data/diagrams/`, `conformal` when it has one in `MFS/assets/data/conformal/`, and `embedding` when it has one in `MFS/assets/data/embedding/`; the page fetches each file only for an entry that says so.

`MFS/assets/data/references.json` is the whole of `assets/data/references.bib` parsed into JSON, so a reader can pull every reference the collection cites in one fetch instead of walking every metric file.
The `.bib` file stays where it is and stays the one place a reference is written; the publications page still reads it.

## The version stamps

A metric's `version` is a hash of that metric's content.
It changes when the content changes and not otherwise, so reformatting a file or republishing the site leaves it alone.
This is what the iOS application uses to fetch only the spacetimes that actually changed, since GitHub Pages rewrites every file's own tag on every publish.
`references.json` carries one stamp of its own at the top for the same reason.
An index entry's `diagrams`, `conformal` and `embedding` are the same kind of hash over that spacetime's diagram file, conformal diagram file and embedding diagram file, so the application can fetch a redrawn diagram without fetching its metric again.

## Checking

    python3 _tools/build_mfs_data.py --check

reports whether the published files are still what the folder says they should be, and changes nothing.

    python3 -m unittest discover -s _tools

runs the tests, which include that check.

## Spacetime diagrams

`MFS/assets/data/diagrams/<metric_id>.json` holds the spacetime diagrams of one spacetime: null rays and future light cones on a plane of the time and one spatial coordinate, for each coordinate system that has any, for Kerr and Kerr-Newman the principal null rays, which leave every such plane, and figures in three dimensions where no plane carries the causal structure, as "Figures in three dimensions" below describes.
The page draws them in a section per coordinate system, after the geodesics, and the application displays the same files.
Nothing in them is drawn by eye or computed by the page.
`_tools/derivations/null_rays.py` computes every ray, cone and marker from the system's published `metric_components`, `inverse_metric_components`, `kretschmann` and `domains`, and for the principal null rays its `weyl_tensor` and `christoffel` as well, read through the checker's own `Reader`, and its docstring records the method.

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
A view of rays that leave the plane names the surface they lie on and the coordinate left out of the drawing, as "the equatorial plane $\theta = \pi/2$ drawn in $t$ and $r$, with $\phi$ left out".
Every curve drawn is a null curve; a caption calls it a null geodesic, the path light takes, only where no Christoffel symbol turns it out of the plane, and says so where one does, as for Gödel and for Kerr off the axis.

A caption is prose in the register of the conventions and histories beside it, and it is read back whole, in the generator and on the page, since a sentence built from clauses that each keep the rules above can still read as machine prose.
It says what a thing is rather than what it is not: "its points are points, not spheres" says nothing that "each point of the diagram is a single event rather than a sphere of them" does not say plainly.
It gives the physics of a surface rather than how the figure lays it out, which the reader can see, and it ends on a fact rather than a flourish, so the lines inside Vaidya's shell crowd toward it because of the slope of $q$, not as "the price of a straight centre".

### Labels are TeX

Every label in a view is text with its mathematics in `$...$`, the form a caption takes: the button's name, the two axes, each tick and the line where a dust solution starts.
An axis names its coordinate and its unit as the entry's line element spells them, `$ct/r_s$` and not `ct / r_s`, so the page sets it with the same MathJax as the metric above it and the application with SwiftMath.
The script picks each axis's ticks, a step of 1, 2, 2.5 or 5 times a power of ten giving about six across, and writes them into the view as `ticks` with the value `at` in the box's units and the label as TeX, so neither reader chooses its own.
The page lays the labels over the SVG as HTML rather than drawing them as SVG text, and `_layouts/mfs.html` says why beside the renderer.

### A diagram is tied to what it was drawn from

Every view records the fields it was drawn from and a stamp over them, computed by `diagram_source_version` in `build_mfs_data.py`, the one function both scripts use.
The fields are the coordinates, the parameter symbols, the metric and its inverse, the Kretschmann scalar and the domains, the Einstein tensor as well for a view whose input is solved as dust, and the Weyl tensor and the Christoffel symbols for a view of the principal null rays.
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
Vaidya is drawn with a shell of null dust of mass $M$, $m$ jumping from $0$ to $M$ at $v = 0$ in the ingoing chart and its time reverse in the outgoing one, the shell the conformal diagram declares.
The Tolman-Oppenheimer-Volkoff star is the polytrope the conformal diagram declares, solved by the `StarSolver` the two scripts share and checked against the published $G^\theta{}_\theta$, which its construction never uses.
Tolman-Bondi is drawn for marginally bound dust whose density at $t = 0$ falls as $1 - r^2/r_b^2$, checked to solve the published $G^r{}_r = 0$, and Oppenheimer-Snyder's collapse is released from rest at twice its Schwarzschild radius, as its conformal diagram declares, in both of its charts.
The Malament-Hogarth toy is drawn for every conformal factor at once, since $\Omega^2$ drops out of the null condition, and checked with a declared factor that grows as $1/|ct|$ along the axis.
The cosmic string is drawn at the deficit its conformal diagram uses, $4G\mu/c^2 = 0.1$, Gödel at $\omega = 1$, and Kerr and Kerr-Newman's figures at the spins and charge of their flat views.

### Checking the rays

    /tmp/mfs-venv/bin/python _tools/derivations/null_rays.py --verify

traces rays as the page's files are traced and measures, for every view with a closed form, how far the quantity each family should conserve drifts along a ray, together with the dust solutions and the equality of Natario's and Alcubierre's metrics on the plane of their axis.
For the principal null rays it adds the checks the next section describes, and the closed forms $ct \mp r_*$ and $\phi \mp r_\sharp$ of Kerr and Kerr-Newman, which the drawing never uses.
On each cylinder of $t$ and $\phi$, van Stockum's and Gödel's, it checks the straight lines against the slopes the captions state, and checks what light launched along each family does against what the caption says, from the published Christoffel symbols: stays on the cylinder as a null geodesic, or is turned toward the axis or away from it.
It writes nothing and exits non-zero if anything fails; run it after changing the method.
It took 958 seconds on 25 September 2026.

### Rays that leave the plane

Off its axis, a light ray of Kerr that runs straight in or straight out turns in $\phi$ as it goes, so no plane of two coordinates holds it.
Such rays are the principal null congruence, and a row with `principal=True` draws it from the published Weyl tensor rather than from a formula for it.
At each point the script takes the Weyl tensor as an operator on bivectors, in a frame the published metric makes orthonormal, requires it to be of Petrov type D, and finds its principal plane: the timelike plane whose two null directions are the Weyl tensor's repeated principal null directions.
That plane's metric, in the basis whose projection onto the two drawn coordinates is their coordinate basis, takes the place of the coordinate plane's metric in the same null condition, so the tracing, the two families, the orientation rules and the cones are the plane method's own, and a cone is the future cone of the principal plane, as the legend says.
The row names in `leaves` the coordinates the rays move in off the drawn pair, and the script refuses the row if the published metric or the Weyl tensor depends on one of them, since only then are the drawn curves the projections of single rays.
It also refuses, naming the point, where the Weyl tensor is not of type D, where the principal plane has a part along a coordinate held fixed, or where it does not project one to one onto the drawn pair.

Kerr and Kerr-Newman each get two such views on the equator, where the ergoregion is widest and the ring singularity lies.
One is the projection onto $t$ and $r$, which carries the cones and both horizons; the other is the equatorial plane seen from above, $r$ and $\phi$ drawn as polar coordinates by `to_display=POLAR`, which carries the turning.
$r$ changes monotonically along every principal ray, so the two projections, which share it, fix each ray.
The view from above seeds twelve rays of each family evenly around the circle inscribed in its box, and `inside=True` stops each a thousandth of the drawing short of the edge of the domain, $r_+$, where $\phi$ winds without end.
Both views mark the ergosurface, $g_{tt} = 0$, with `mark_gtt`, whose marker carries its own legend.
Beside a curvature singularity the published components lose to rounding the digits the principal plane is found from, a direction error of rounding over $r^4$ beside Kerr's ring, so the plane is not taken where the Kretschmann scalar passes $10^{12}$ in the view's units, $r = 0.006\,GM/c^2$ beside the ring and under a pixel from it, and a ray stops there as it does wherever the metric is not finite.

Before a principal view is written, the script checks at points along the drawn $r$ that its rays satisfy the geodesic equation with the published Christoffel symbols and that their directions satisfy $C_{abc[d}k_{e]}k^bk^c = 0$ with the published Weyl tensor, and refuses to write the view if either misses.
`--verify` adds the closed forms, shows that the projection onto $t$ and $r$ is the same at $\theta = \pi/5$ as on the equator and equals the plane of $t$ and $r$ on the axis, shows that in Schwarzschild, Reissner-Nordstrom and Taub-NUT the principal plane is the plane of $t$ and $r$ their radial views already draw, and shows the script refusing van Stockum, whose Weyl tensor is of type I.

Kerr and Kerr-Newman are the only charts the extension draws something new for.
Every spherically symmetric chart that is not conformally flat has the plane of $t$ and $r$ as its principal plane, and so does Taub-NUT, whose twist sits in $g_{t\phi}$; Gödel's is the plane of $t$ and $z$ along its axis of rotation, which is flat.
van Stockum's Weyl tensor has four distinct principal directions and no repeated one, so its light rays, which frame dragging turns out of the plane of $t$ and $r$, form no principal congruence, and the pp-wave's one repeated direction, $\partial_v$, lies in the plane of its axis, already drawn.

### Figures in three dimensions

Where the causal structure a reader comes for turns in a direction no plane of two coordinates holds, a coordinate system also gets a figure in three dimensions: a slice of one time and two spatial coordinates, every other coordinate held fixed, drawn from the published metric by `_tools/derivations/projections.py` and projected once from a fixed camera.
The file holds the projection, polylines, polygons and points of the page's own plane in the order they are painted, farthest first, with TeX labels, which is the form a conformal diagram takes.
So the page draws a figure with the conformal diagram's frame and labels, it prints, and the application can draw the same data with its own renderer; nothing is left for a reader to turn, and there is no scene in the browser.
`FIGURES` and `CAPTIONS` in that script are its table, and `null_rays.py` writes each figure into its spacetime's diagram file under `projections`, keyed by coordinate system beside the flat views under `systems`, which an application that reads only the flat views passes over.
A figure is stamped with the fields it was drawn from, as a flat view is.

van Stockum's cylinder and the Gödel universe, in Gödel's cylindrical chart, each get the slice $z = 0$ of $t$, $r$ and $\phi$, with future light cones on the axis and at four places around each of the circles $r_c/2$, $r_c$ and $3r_c/2$, where $r_c$ is the zero of the published $g_{\phi\phi}$.
The cones tip over toward $+\phi$, one edge of each lies along the circle at $r_c$, and beyond it the circle is a closed timelike curve.
van Stockum's radius is drawn as the proper distance from the axis and Gödel's as $r$ itself, since its $g_{rr} = -g_{tt}$, so that in both the null directions straight out from the axis run at 45°.
Each also gets the cylinders of $t$ and $\phi$ at $r_c/2$ and $3r_c/2$ as flat views, drawn by the plane method with $\phi$ scaled by $r$ and `periodic` keeping the cylinder's two edges, which are one line, from being hatched.
Kerr and Kerr-Newman each get the equatorial slice of $t$, $r$ and $\phi$, drawn with $r$ as the radius down to the horizon, where the chart ends, with cones at four places around each of three circles: halfway through the ergoregion, on the ergosurface and at $3r_E/2$.
The horizon is the outer zero of the published $g^{rr}$ and the ergosurface the zero of the published $g_{tt}$ outside it, and before the figure is drawn every future generator on the inner circle is checked to move toward $+\phi$, none on the ergosurface to move back, and some on the outer circle to do so.
That is what no plane can show, since inside the ergoregion the plane of $t$ and $r$ at fixed $\phi$ has no null direction at all.
A figure of the principal null congruence alone would add nothing to the two flat views: $r$ changes monotonically along every principal ray, and its $t(r)$ and $\phi(r)$ are the same at every $\theta$, so the two projections fix every ray at every latitude.
In Boyer-Lindquist $t$ the cones near the hole are narrow, so they are drawn larger than van Stockum's, and with a cone every 45° round the ergosurface no label fits beside its circles, so the legend names them.
Alcubierre's bubble gets the slice $z = 0$ of $t$, $x$ and $y$ at the moment its declared profile centres it on $x = 0$, drawn cartesian with $t$ up: cones at the centre and around the circles $r_s = R/2$, $v_sf = 1$ and $r_s = 2R$, and the centre's world line $x = 2ct$.
The circle $v_sf = 1$ is the zero of the published $g_{tt}$, the same circle about $x = 2ct$ a little later is checked to be its zero there too, and the world line is checked against the published metric to have $g(u, u) = -1$, so that its clock keeps $t$.
The flat view on the axis shows the tilt along the direction of motion; the figure shows the same shift tipping every cone along $x$ off the axis too, so that the cones tip over in a ball about the ship.
Natário's drive has no such figure yet, Lentz's soliton cannot be written down, and the Krasnikov tube's causal structure lies along its axis, where the plane of $t$ and $x$ already draws it.
The cosmic string's figure is its plane $z = 0$ seen from straight overhead with $t$ left out: a beam of light rays run with the published Christoffel symbols past the string, laid flat with the angle $(1 - 4G\mu/c^2)\phi$ so that each ray is checked straight, and carried across the wedge the deficit removes.

A cone is the convex hull of its apex and its rim, every generator one Euclidean length in the drawing, so its size says nothing and its shape and tilt are the metric's; every generator is checked null against the published metric, and the published inverse is checked to be the inverse of the published metric on the slice.
Beyond $r_c$ the cones are so wide that from most directions the camera looks into their opening, and a cone projects to an oval with its apex inside; eight of its generators are therefore drawn faintly from the apex to the rim, which is what shows where the apex is and which way the cone opens.
A cone whose axis points along the camera's own direction is the worst case, so the cones on $r_c$ are turned by 45° from those on the other circles and none is placed where its neighbour's rim would cover it.

Every class a figure paints needs a style in `_layouts/mfs.html`, as `.pj-<class>` on the screen and under `.mfs-print-body` in print, and a fill as `.pj-<class>-fill`.
A path the stylesheet does not know is painted as nothing at all, so a test holds every class of every figure to having both.

A coordinate system with neither a flat view nor a figure has no section on the page, and a spacetime with neither has no diagram file, which a full redraw removes if one is left behind; `build_mfs_data.py` refuses a diagram file that draws nothing.

### What is not drawn

Kerr's and Kerr-Newman's planes of $t$ and $r$ are drawn on the axis only, and off it their principal null rays are drawn instead.
Off the axis the null curves of fixed $\theta$ and $\phi$ are not null geodesics, since $\Gamma^\theta{}_{tt}$, $\Gamma^\theta{}_{rr}$ and $\Gamma^\phi{}_{tr}$ turn a light ray launched along one out of the plane, and inside the ergoregion the plane has no null direction at all, which is why the equator is also drawn in three dimensions.
van Stockum's plane of $t$ and $r$ is not drawn: every null geodesic leaves it, and its Weyl tensor, of type I, gives no congruence to draw them by, so its cylinders of $t$ and $\phi$ and its figure are drawn instead.
The proper distance chart of Morris-Thorne leaves its functions free with no choice made for them, while its spherical chart draws the Ellis-Bronnikov member.
Natário's plane flow chart, the Brinkmann chart of the pp-wave and Minkowski's double null chart would each only repeat a plane drawn elsewhere, flat or the same as another chart's.

Nothing at all is drawn of Lentz's soliton, the Mixmaster universe or Gott's core.
Lentz's soliton exists only as a numerical integral over his rhomboid sources, so no soliton of the class can be written down, and a potential written in its place would draw another soliton of the class.
The Mixmaster's scale factors reach the singularity through an endless sequence of Kasner epochs, which a drawing of one solution follows only for a handful, and the one plane of time and an Euler angle that keeps its light rays, $t$ against $\psi$, would follow $a_3$ alone.
Gott's core, the cosmic string's interior, is flat on its plane of $t$ and $\chi$, and acts on light across it, in the cap of $\chi$ and $\phi$.

## Conformal diagrams

`MFS/assets/data/conformal/<metric_id>.json` holds the conformal diagram of one spacetime: the whole spacetime brought to a finite drawing with light at 45°, or, where no picture of the whole is faithful, a totally geodesic surface in it that says so.
The page draws it under the heading "conformal diagram", below the coordinates and the conventions, and it follows the chart chosen there: it shows the views that name that chart, then the views that name no chart, with buttons only when that makes more than one.
A chart none of whose views is its own shows no conformal diagram at all, as Vaidya's outgoing chart does, since the one drawn is the shell imploding in the ingoing chart.
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

`DRAWN` in the script names the seventeen spacetimes that have a diagram and `NOT_DRAWN` the twelve that have none: every event's future is the whole spacetime for Gödel and van Stockum, the diagram is whatever a free function makes it for the warp drives, the Krasnikov tube, Tolman-Bondi and Bianchi, the Mixmaster has no surface that carries its causal structure, and so on.
A free function does not by itself rule a diagram out: where every choice of it gives the same shape, as for the Tolman-Oppenheimer-Volkoff star and the Morris-Thorne wormhole, the diagram is drawn with a declared choice that moves only the lines inside, and the Malament-Hogarth toy is drawn for every conformal factor at once, since a conformal factor changes no null direction.
Those twelve have no file, which a full redraw removes if one is left behind, so the page has no conformal diagram section for them; `build_mfs_data.py` refuses a file that draws nothing, and the script stops if a metric file is in neither table, so a new spacetime needs a decision.

A view of a surface that is not the whole spacetime carries `restriction`, which the page prints in a band across the top of the figure, never in a footnote.
Kerr and Kerr-Newman are drawn on the symmetry axis and the cosmic string on the half plane of fixed $\phi$ and $z$, and the tests hold those three to carrying a restriction on every view.

Every text in a view is TeX in `$...$`: the labels on the drawing, the buttons, the legend, the caption, the restriction and the parameter values.
A caption opens by naming what is drawn, "This is the whole of ..." or "This is the plane ...", and its prose is for the reader of "Whom the prose is for" above; the tests hold every text in these files to that rule's words and to the dash rule.
It keeps the register the spacetime diagrams' captions keep, and it says what each point of the diagram stands for: a sphere of radius $r$ for the whole of a spherical spacetime, a single event on a slice such as Minkowski's plane $y = z = 0$, and a circle about the string times a line along it on the cosmic string's half plane.

`_tools/derivations/tex_check.cjs` typesets every TeX string the page sets, in the metrics, the diagram files and the conformal files, a published value together with its negation, through the TeX input MathJax loads on the page, and exits non-zero naming each one it cannot set.
sympy never reads the typesetting, so this is the check that catches a value that is right and prints as an error box:

    npm install --prefix /tmp/mfs-node mathjax-full
    node _tools/derivations/tex_check.cjs /tmp/mfs-node

## Embedding diagrams

`MFS/assets/data/embedding/<metric_id>.json` holds the embedding diagram of one spacetime: a slice of it, the equatorial plane at one moment, drawn as a surface in ordinary flat three dimensional space so that every distance along the surface is the distance the metric gives.
The page draws it under the heading "embedding diagram", just below the conformal diagram, and the application reads the same files.
The file is the definition in "The file, which the application reads" below, and the application is built against that section, so a change to the shape of the file is a change to it first.

`_tools/derivations/embedding.py` draws them, one function per spacetime with its derivation in its docstring, reading the published metric through the checker's `Reader` with the `load` and `published_matrix` the null rays use.
It needs the same environment as `null_rays.py` and took 41 seconds for the whole collection on 27 September 2026, most of it in the two collapsing dust clouds, whose every shell is solved on its cycloid by bisection:

    /tmp/mfs-venv/bin/python _tools/derivations/embedding.py
    python3 _tools/build_mfs_data.py

`--metric <metric_id>` redraws one spacetime, and `--verify` prints every check and writes nothing.

### The construction

Every slice drawn is the surface of one spatial coordinate $x$ and the angle $\phi$ of one coordinate system, every other coordinate held fixed, with the metric $g_{xx}(x)\,dx^2 + g_{\phi\phi}(x)\,d\phi^2$ read from the published `metric_components`.
The script refuses a slice with a cross term $g_{x\phi}$ or with a component that still depends on anything but $x$, since only then is it a surface of revolution.
A circle of constant $x$ has circumference $2\pi\sqrt{g_{\phi\phi}}$, so it is drawn at the radius $\rho = \sqrt{g_{\phi\phi}}$ from the axis, and the distance $\sqrt{g_{xx}}\,dx$ out to the next circle is the hypotenuse of $d\rho$ and $dz$, so $dz/dx = \sqrt{g_{xx} - (d\rho/dx)^2}$.
With an areal radius, $\rho = r$, that is the familiar $dz/dr = \sqrt{g_{rr} - 1}$.
The surface exists exactly where $g_{xx} \ge (d\rho/dx)^2$, and where a circle grows faster than the distance out to it no surface of revolution in flat space carries the slice.

$d\rho/dx$ is taken in sympy as $g_{\phi\phi}'/(2\sqrt{g_{\phi\phi}})$, so no absolute value is ever differentiated, and $z$ is the adaptive quadrature of $\sqrt{g_{xx} - (d\rho/dx)^2}$ between neighbouring points.
An end where the integrand diverges, as $g_{rr}$ does at a throat, is moved to the end of a new variable $s$ with $x = x_0 \pm (b - a)s^2$, which makes the integrand finite there.
A horizon found by `Slice.horizons()`, a root of $1/g_{xx}$, is kept exact, and next to it the metric is evaluated in $u = x - x_0$, expanded in sympy with the root exact, since the factor that vanishes there, as $r^2 - r_sr + r_q^2$, loses every digit to cancellation in floating point within $10^{-13}$ of an irrational root such as Kerr's $r_+ = 1 + \sqrt{19}/10$.
At such a horizon $1/g_{xx} = 0$ while $d\rho/dx$ is finite, which is checked, so the surface's tangent is vertical there exactly.
An interval is halved until the profile strays from its chord by less than $2 \times 10^{-5}$ of the drawing's size and by less than $0.004$ of the chord, and the chord is shorter than $1/90$ of the size, so a chord falls short of its arc by less than about $4 \times 10^{-5}$ of it.
No closed form is used to draw anything; the closed forms are only checked against.

### Checking the surface

The surface is measured as the application will draw it, from the rounded numbers the file holds, against the published metric, and the script refuses to write while any check fails:

- along: each chord of each profile, and each profile end to end, against the proper distance between the same two values of $x$;
- across: the straight line in space from each point to the next one $0.02$ further round the axis, against the length the metric gives the line that runs out at a steady proper distance while it turns steadily through the same angle;
- around: $\rho$ at every point against $\sqrt{g_{\phi\phi}}$;
- joins: where two pieces meet, as a star's surface meets the exterior, they meet at one point with one tangent, which says $g_{xx}$ agrees on both sides;
- forms: each surface against the closed form it is known by: Flamm's paraboloid, also as the exterior of the neutron star, of each collapse at its release and, moved in by $r_s$, of Vaidya's slices; the interior Schwarzschild cap, the catenoid in both its charts, the cone, Gott's cap, the spheres of FRW, de Sitter and Oppenheimer-Snyder's dust, and Bertotti-Robinson's cylinder and sphere; and the circumference radius of the Kerr, Kerr-Newman and Taub-NUT horizons;
- fields: a declared star, scale factor or dust cloud against the published Einstein tensor, the neutron star against $G^\theta{}_\theta = 8\pi p$, FRW's and Oppenheimer-Snyder's scale factors against their $G^r{}_r$ and $G^\chi{}_\chi$, and each cloud released from rest against Tolman-Bondi's $G^r{}_r = 0$ and its $G^t{}_t$, the density;
- stops: where a view or a file says a slice cannot be drawn, $g_{xx} - (d\rho/dx)^2$ is negative at every sample, or $g_{\phi\phi}$ is, where the circles are timelike; where it says a slice is a plane, it is zero, or the spatial metric has no cross term and no component that depends on a spatial coordinate.

On 27 September 2026 the worst chord missed its proper distance by $7.0 \times 10^{-5}$ of it, inside Reissner-Nordstrom's inner horizon, and the worst line across by $1.1 \times 10^{-4}$, on the sphere of radius $a = 1/2$ of FRW's first moment, where the rounding on chords a sixtieth of a unit long is most of it.
Every profile end to end was within $3.3 \times 10^{-5}$ of its proper length, every $\rho$ within $2.5 \times 10^{-8}$ of the drawing's size of $\sqrt{g_{\phi\phi}}$, every closed form within $2.5 \times 10^{-8}$ of the size, every join met to $3 \times 10^{-13}$ with tangents equal to $4 \times 10^{-16}$, and every declared solution satisfied the published field equations to $1.4 \times 10^{-14}$; 382 checks in all, none failing, and running the script again writes the same bytes.
The tests hold the files on disk to the same closed forms without sympy, from the numbers written and nothing else, so a redraw that changed a surface fails there too.

### Which spacetimes

`DRAWN` in the script names the spacetimes that have a diagram, `STATED` those with no surface to draw whose file says why, and `NOT_DRAWN` the others, which have no file; the script stops if a metric file is in none of the three or in more than one.
Each function in `STATED` checks what it says from the published metric before it says it.
Alcubierre's, Natario's and Lentz's warp drives, Kasner's universe, Bianchi type I and Minkowski space are flat on every slice of constant $t$, their spatial metric having no cross term and no component that depends on a spatial coordinate, so the equator is a plane and the physics is in how the slices are stacked.
Krasnikov's tube is flat outside, where $k = 1$, and deep inside it $k < 0$ makes the direction along the tube timelike at constant $t$, so a surface of constant $t$ is not a moment of space there.
Anti-de Sitter's static slice is the hyperbolic plane, $g_{rr} < (d\rho/dr)^2$ at every $r > 0$, which is checked.
Mixmaster's slices are squashed three spheres and a pp-wave's spacelike slices carry its profile, so neither is flat and neither has one surface that says anything, and a Malament-Hogarth spacetime's slices take whatever shape its arbitrary conformal factor gives them; those three have no file.
Schwarzschild is Flamm's paraboloid on both sheets of the Einstein-Rosen bridge, the slice of constant $t$ running through the bifurcation sphere into the other exterior.
The interior Schwarzschild star, at $R = 1.5\,r_s$ as its conformal diagram declares, is the cap of a sphere of radius $\sqrt{R^3/r_s}$ joined to the exterior of `schwarzschild.json`, which the file reads, with the vacuum paraboloid drawn on under the cap down to the throat the star replaces.
The Tolman-Oppenheimer-Volkoff star is the polytrope its other diagrams declare, $M = 1.40\,M_\odot$ and $R = 14.2$ km, with the mass function `null_rays.StarSolver` solves from the published Einstein tensor entering the slice as numbers; its exterior is checked to be Flamm's paraboloid of that mass, and the vacuum paraboloid is drawn under it as under Schwarzschild's star.
Morris-Thorne is drawn as its Ellis-Bronnikov member, $b = b_0^2/r$, the catenoid, as its other diagrams declare; the surface reads $b$ alone, and its proper radial chart, with $r(l) = \sqrt{l^2 + b_0^2}$, is checked to give the same surface.
Reissner-Nordstrom, at $r_q = 0.48\,r_s$ as its conformal diagram declares, has two views: outside $r_+$ the two sheets through the outer horizon's bifurcation sphere, and inside $r_-$ the two sides through the inner one, its widest circle there, each running in until $g_{rr} = 1$ at $r = r_q^2/r_s$, where the surface lies level; nearer the singularity and between the horizons no surface carries the slice, which the views state under `stops` and the script checks.
Kerr, at $a = 0.9\,GM/c^2$, and Kerr-Newman, at $a = 0.6\,GM/c^2$ and $r_Q = 0.5\,GM/c^2$, as their conformal diagrams declare, are the equator of a slice of constant Boyer-Lindquist $t$, where $g_{t\phi}$ drops out, drawn at the circumference radius $\rho = \sqrt{g_{\phi\phi}}$ through the bifurcation sphere at $r_+$ into a second exterior, with the edge of the ergosphere marked as the class `ergo`; Kerr's throat is checked to have $\rho = 2GM/c^2$ at any spin and Kerr-Newman's $2GM/c^2 - r_Q^2/r_+$.
De Sitter's static slice at $t = 0$ is the sphere of radius $\ell = \sqrt{3/\Lambda}$, drawn at $\Lambda = 3$: the static chart covers the hemisphere out to the horizon and the antipodal observer's patch the other, the waist of the hyperboloid; the flat slicing's slices are flat, which the view states and the script checks.
Vaidya is its imploding shell of radiation at four moments, slices of constant $v - r$, since a slice of constant $v$ is null; `Slice` pulls the published metric back along $v = T + r$ to $(1 + 2Gm/c^2r)\,dr^2 + r^2d\phi^2$, a flat disc inside the shell, which is checked, and $z^2 = 4r_sr$ outside, Flamm's paraboloid moved in by $r_s$, meeting at a `crease`; the figure sets the four in two rows.
Oppenheimer-Snyder is its dust released from rest at $R_0 = 2\,r_s$ at four moments of the dust's proper time: inside, the published interior chart at $a = (a_m/2)(1 + \cos\eta)$, checked to make the published $G^\chi{}_\chi$ vanish, a cap of a sphere of radius $a$; outside, the moment carried on as Novikov's slice, clocks released from rest at every radius with the dust, read from Tolman-Bondi's comoving chart with no dust in it, since a slice of constant Schwarzschild $t$ meets the dust at an angle after the release and cannot reach it inside $r_s$. `RestCloud` solves each shell's cycloid and is checked to make Tolman-Bondi's published $G^r{}_r$ vanish and its $G^t{}_t$ give the declared density; the two sides meet with one tangent, and at the release the outside is Flamm's paraboloid.
Tolman-Bondi is the cloud its spacetime diagram declares, density falling as $1 - r^2/r_b^2$ with $2GM/c^2 = r_b/2$, but released from rest, $E = -GM(r)/c^2r$, at four moments before its centre is crushed, from its own comoving chart with `RestCloud`, checked against its published field equations; the spacetime diagram's cloud is marginally bound, $E = 0$, and its slices are planes, which the view states and the script checks.
Bertotti-Robinson, a product of two dimensional anti-de Sitter space and a sphere of the one radius $b$, has two views: its equator at one moment, a line of the first factor times a great circle of the second, the cylinder $\rho = b$, $z = b\ln r$, and its sphere of $\theta$ and $\phi$ at one event, the second factor; the first factor is Lorentzian and has no surface, which the view states.
Van Stockum's rotating dust, at $R = 1$ as its spacetime diagrams declare, is the plane $z = 0$ at one moment, whose circles grow to $r = R/\sqrt{2}$ and shrink after, so the surface curls back toward the axis until $g_{rr} = (d\rho/dr)^2$ at $r = 0.83\,R$, found by sympy; the view states that nothing is drawn beyond, and that beyond $r = R$ the circles are closed timelike curves, both checked.
Taub-NUT, at $m = 1$ and $l = m/2$ as its spacetime diagram declares, is the equator of a slice of constant $t$, where $g_{t\phi} \propto \cos\theta$ vanishes, out of its horizon $r_+ = m + \sqrt{m^2 + l^2}$, whose circumference radius $\sqrt{2(mr_+ + l^2)}$ is checked; across the horizon lies Taub's cosmology, where $r$ is a time, which the view states and the script checks.
Godel's universe, at $\omega = 1$ as its spacetime diagrams declare, is the plane $z = 0$ about one world line of the dust at one moment of the cylindrical chart, which stops exactly where $\sinh^2 r = 1/\sqrt{2}$, since $g_{rr} = (d\rho/dr)^2$ there reduces to $4s^3 = 2s$; the view states that beyond it no surface carries the slice and that beyond $\sinh r = 1$ the circles are closed timelike curves, both checked.
Ellis-Bronnikov is the same catenoid in its own chart, whose $r$ is the proper distance from the throat, so it is one piece through the throat, checked against $z = \ell\,\mathrm{arcsinh}(r/\ell)$ and $\rho = \sqrt{r^2 + \ell^2}$.
The cosmic string is its cone at $4G\mu/c^2 = 0.1$ with Gott's core rounding the apex, and the figure lays the cone flat beside it, cut along one meridian, so that the missing wedge $\delta = 36°$ shows.
FRW is its closed universe of dust at five moments, $a = 1 - \cos\eta$ as its conformal diagram declares, checked to make the published $G^r{}_r$ vanish; its flat slices are planes and its open slices have no surface of revolution in flat space, and no surface at all as a whole, which the view states under `stops` and the script checks at the scale factors of dust.

### The file, which the application reads

The application downloads the file for a spacetime whose index entry carries `embedding`, the stamp of the whole file, which changes when the file does and not otherwise.
Every number is a plain JSON number, every text is TeX in `$...$` inside prose, as the conformal diagram's texts are, and nothing in the file needs any relativity to draw.

    {
      "metric": "schwarzschild",
      "source": [{"metric": "schwarzschild", "system": "spherical", "fields": [...], "version": "..."}],
      "views": [view, ...]
    }

`source` is what `build_mfs_data.py` checks the file against, one entry for every coordinate system it read, which may be another spacetime's, as the star reads `schwarzschild.json`; the application can ignore it.
`views` holds at least one view, unless the spacetime has no surface to draw.
Then `views` is empty and `stops`, a list of TeX sentences, says why, printed under the heading after "not drawn", as for the warp drives, whose slices are flat, and anti-de Sitter, whose slices are hyperbolic planes.
A file with a view never carries `stops` of its own; what a view does not draw is in the view's `stops`.
A view is:

- `id`: unique in the file.
- `label`: the name of its button, TeX, wanted only when a chart has more than one view to choose from.
- `system`, optional: the coordinate system the view belongs to. A view without it belongs to every chart. Show the views that name the chart being read and then the views that name none, as the conformal diagram does, and nothing if that leaves none.
- `unit`: TeX naming the length every number of the surface is measured in, as `$r_s$`, `$b_0$`, `$\ell$` or `$1/\sqrt{k}$`.
- `surfaces`: at least one surface. One surface is one moment; more than one is a sequence of moments in the order of their `time`.
- `figure`: the page's drawing of the view, described below; an application that draws the surfaces itself can ignore it.
- `caption`: a list of paragraphs.
- `settings`, optional: TeX prose giving the parameter values drawn at, printed after "drawn at". It names the length every number is measured in, as "$r_s = 1$, the unit of every length".
- `input`, optional: TeX prose naming a declared function or matter, printed after "drawn with".
- `stops`, optional: a list of TeX sentences, each saying where the construction stops or what in this spacetime has no surface, printed after "not drawn". FRW's view carries its flat and open universes here, which have no surface in the file.

A surface is:

- `label` and `time`, in a sequence only: the moment as TeX, and its time as a number in `unit`, the time the view's `settings` name: $ct$ for FRW and $v - r$ for Vaidya.
- `pieces`: its profile curves, at least one.
- `rings`: circles marked on it, possibly none.

Every surface is a surface of revolution about the vertical axis $z$, and every piece is a profile curve in a half plane through the axis.
Its `points` are `[x, rho, z]`, running from its `start` to its `end` with `x` strictly increasing or strictly decreasing, at least two of them:

- `rho` is the distance from the axis, never negative, and `z` the height, both in `unit`, rounded to below $10^{-7}$ of the piece's own extent and to at least six decimals, so a small piece, as a collapsing star's shrinking cap, keeps the precision of a large one;
- `x` is the value of the piece's `coordinate` at the point, as that coordinate system writes it with the view's `settings`, a length in `unit`, an angle in radians such as $\chi$, or a number such as FRW's comoving $r$, written as the double it was computed at, the shortest decimal that reads back as it, since next to an irrational horizon a rounded $x$ would fall inside it; a client needs it only to label a point or to find one.

No number is written as `-0.0`.
A piece holds the surface for the values of `x` from its first point to its last and says nothing beyond them; its `start` and `end` say what the surface does at each.
The pieces drawn on 27 September 2026 have 9 to 257 points each.

To draw a piece, turn every point about the axis, $(\rho\cos\phi, \rho\sin\phi, z)$ for $\phi$ from $0$ to $2\pi$, and join neighbouring points and neighbouring angles.
The points are close enough that straight segments between them are the surface to within $2 \times 10^{-5}$ of the drawing's size, the diameter of the widest circle drawn, so no smoothing is wanted, and a client may add angles as finely as it likes.
The pieces of a surface may be drawn in any order, since every piece but a reference piece is part of one surface and no two of them overlap.
Draw $\rho$ and $z$ at one scale, since a surface stretched along either axis no longer has the metric's distances.
The zero of $z$ is wherever the construction put it and means nothing; only differences in $z$ do, and every piece of one surface shares one $z$.
In a sequence every surface stands on its own axis at the origin, and placing them side by side or showing them one after another is the client's choice.

A piece also carries:

- `id`: unique in its surface.
- `class`: `sheet` for the slice a chart of the spacetime covers, `sheet2` for the same surface run on past where that chart ends, as the other exterior of Schwarzschild, the other side of a wormhole or the far hemisphere of a closed universe, `star` for the slice through matter, as a star or Gott's core, and `reference` for a surface that is not part of the slice at all, drawn for comparison.
- `reference`: `true` on a reference piece, and absent otherwise. A reference piece is drawn faintly and hides nothing, as the vacuum paraboloid under the star and the cone of an ideal string under Gott's core are.
- `metric` and `system`: the spacetime and the coordinate system the piece was read from, one of the file's `source` entries.
- `coordinate`: the coordinate `x` is, as that coordinate system spells it, as `r` or `\chi`.
- `start` and `end`: what the surface does at each end of the piece, `{"kind": ..., "text": ...}`, the text TeX prose for a reader and sometimes absent.

The kinds of end are:

- `axis`: the piece reaches the axis, $\rho = 0$, and the surface closes there smoothly, as at a star's centre or a sphere's pole.
- `apex`: the piece reaches the axis at a point where the surface is not smooth, as the cone's apex, where the string lies.
- `join`: the piece meets another piece of the same surface, which has the same $\rho$ and $z$ there and the same tangent.
- `crease`: the piece meets another piece of the same surface at the same $\rho$ and $z$ but at an angle, where a thin shell of matter folds the surface, as Vaidya's shell of radiation does.
- `throat`: the piece reaches a smallest circle, where its tangent is vertical and the surface turns back out as its own mirror image in the plane of that circle. The mirror image is the view's `sheet2` piece when the view draws the other side, as for Schwarzschild and the wormhole, and is not drawn otherwise, as under the star.
- `edge`: the drawing ends here and the surface runs on; draw the circle there as the rim of the drawing, or fade the surface out.
- `stops`: the construction ends here, because past this point no surface in flat space carries the slice, and the text says why.

A ring is a circle to mark on a surface, `{"piece", "class", "x", "rho", "z", "label"}`, of radius `rho` at height `z`, which is one of the points of its piece, and `label`, when present, is TeX to set beside it.
Its `class` is `r` for a circle of constant coordinate on a `sheet` or `star` piece, `r2` for one on a `sheet2` piece, `horizon` for the horizon, `throat` for a wormhole's throat, `surface` for the edge of matter, a star's surface or the edge of Gott's core, `chartedge` for the circle where a chart ends, as FRW's equator, and `reference` for a circle on a reference piece.

The `figure` is the drawing the page makes of the view, projected by the script once from a fixed camera, in the form a figure in three dimensions takes in the diagram files:

- `box`: `[Xmin, Xmax, Ymin, Ymax]` of the plane of the page, with $Y$ up, drawn at one scale.
- `camera`: `{"azimuth", "elevation"}`, the angles in degrees the surfaces were projected from, the azimuth measured round the axis from $\phi = 0$ toward $\phi = \pi/2$ and the elevation above the plane $z = 0$, looking at the origin.
- `layers`: painted in the order given, each with a `kind` and a `class`. A `fill` has `points`, a closed polygon, and optional `holes`, polygons painted with the even odd rule; a `line` has `points`, a polyline; a `point` has `at`. Every point is `[X, Y]` in the plane of the page, rounded to four decimals.
- `labels`: TeX at a point, `{"at", "text", "anchor", "class", "dx", "dy"}`, the anchor and the offset in units of a figure 628 wide, as a conformal diagram's are, `class` being `lab` or `small`.
  A `small` label names a circle and stands beside the end of it, on the surface and over its lines, so it is set on a ground of its own, dark on the screen and white in print, as the page sets it; a `lab` label stands where no line runs.
  No two labels overlap, which the script checks before it writes.
- `legend`: `[kind, class, text]` for each class it names, the kind being `fill`, `line` or `point`.

A `sheet` piece is tinted `cover` and a `star` piece `star`, while a `sheet2` or reference piece is left clear, so the part a chart covers stands out; each is tinted only where it is the surface nearest the camera, found by casting rays through the same truncated cones a client draws, and every line on the surfaces is split where another part of a surface hides it, the hidden part carrying its class with `-far` appended, which the page draws faint.
A reference piece hides nothing.
The fills are `cover`, `star` and `wedge`.
The lines are `outline`, where the surface turns edge on to the camera and the rim at an `edge` or `stops` end; `meridian`, the profile turned to evenly spaced angles, as the legend says; `reference`, the meridians and outline of a reference piece in dashes; the ring classes above, each ring drawn as its circle, or as a `point` where its $\rho$ is zero; and `cut`.
Each line class may also come with `-far`.
The page styles each class as `.em-<class>` for a line and `.em-<class>-fill` for a fill, on the screen and in print, and a test holds every class a figure paints to having both.
The cosmic string's figure also lays the cone flat beside it, a flat drawing in the plane of the page with the wedge the cone lacks painted as `wedge` and the two edges of the cut as `cut`.

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

`tov`, `malament_hogarth`, `mixmaster` and `lentz`, and the cylindrical chart of `godel`, have their mathematics written by `_tools/derivations/print_charts.py`, which defines each of those charts, computes every tensor with the checker's own `Geometry`, and prints each value through `chart_printer.py` beside it:

    /tmp/mfs-venv/bin/python _tools/derivations/print_charts.py --metric tov

Every printed value is read back through the checker's `Reader` and compared with the value it came from before anything is written, and `verify_metrics.py` then checks the file like any other.
The script writes only the coordinate system's mathematics, the line element and the domains; the prose, the convention and the description of each parameter stay in the metric file and are carried over, and a parameter with no description stops the script.
A chart can pass a function that regroups the numerator of each value, which is how TOV's curvature is printed around $(\partial_r\Phi)^2 + \partial_r^2\Phi$ and Bianchi IX's around $a_1^2\cos^2\psi + a_2^2\sin^2\psi$, and a Ricci or Kretschmann scalar can be given in a structured form, which is checked against sympy before it is used.
A chart can instead pass a whole `pretty`, as Gödel's cylindrical chart passes `chart_printer.hyperbolic`, which rewrites the exponentials the checker's `Geometry` hands back in $\sinh r$ and $\cosh r$.
A chart the script writes replaces the chart of its id, or joins the spacetime's other charts after them, so Gödel's Cartesian chart stays as it was written.
A parameter keeps the description it has in that chart, or in another chart of the same spacetime.
`tov.md`, `malament_hogarth.md`, `mixmaster.md`, `lentz.md` and `godel.md` in the same folder record why each chart is the one published; `godel.md` also records the transformation from the Cartesian chart, whose pullback sympy checks symbolically.
The other charts were written by hand and are not touched by the script.

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
