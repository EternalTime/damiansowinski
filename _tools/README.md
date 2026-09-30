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
Then run `python3 _tools/build_agent_data.py`, which adds the spacetime to `/data/spacetimes.json`, `/llms.txt` and `/llms-full.txt` for agents, and commit what it rewrote too.

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
`fitProseMath` in `_layouts/mfs.html` measures that again at the next frame after the paragraph changes width, never from inside the observer that sees the change, and passes over a displayed formula, which has its line already.
Each of these observers measures everything it was handed before it changes any class, because a class changed between two measurements makes the second lay the whole page out again, and once per line that held Kerr-Newman still for seconds and Natario for nearly a minute.

A chosen spacetime's panel waits off the screen until its mathematics is set and laid out, and only then slides in.
While it waits, "Spacetime data loading..." shows where the panel rests, centred, under the panel so the panel slides in over it.
Most of that wait holds the main thread, so the line cannot be shown once it starts; it is asked for as the wait begins and fades in by a CSS transition with a delay, which the browser runs off the main thread, so a spacetime set within that delay never shows it.
`showNote` in `_layouts/mfs.html` carries the timings and why.

## The panels on a desktop

Wherever the list and the spacetime stand side by side, the list's panel starts on the spacetime panel's top and the coffee panel ends on its bottom, with 6px between the two on the left, as the captain asked on 29 September 2026.
The spacetime panel's top is halfway between the foot of the Exit sign and 150px, or 30px below the title's foot where a large text size wraps the title lower.
`placeTop` and `syncPanel` in `_layouts/mfs.html` carry that geometry in `--mfs-top`, `--mfs-bottom` and `--mfs-coffee-h`, and `page_timing.mjs` holds every state of the page to it within a pixel.

## The page on a phone

A screen narrower than 600px, or a touch screen under 500px tall, gets the same panels in one column that the page scrolls through: the title, the list with the coffee panel, then the spacetime.
That covers every iPhone upright and on its side, and the desktop layout is untouched by it.
The rules are one `@media screen` block in `_layouts/mfs.html`, after the main stylesheet, and its comment says why each choice was made.
They are for the screen alone, so the print copy is the same whatever screen it was printed from; a phone printing Kerr differs from a desktop only by MathJax's rounding at three pixels to the point, and did before this layout.

The desktop's panel scripts write their geometry into each panel's own style, so the phone block overrides it with `!important`, and it sets `--mfs-phone` on the root, which is how the scripts tell which layout is in force.
Anything new on the page has to hold at 390pt upright and 844 by 390 on its side: a line of mathematics, a coordinate domain or a formula in prose scrolls on its own there, and nothing may make the page wider than the screen, even for a frame, since a phone answers that by zooming the whole page out.
A spacetime diagram is drawn whole, as wide as the panel and never taller than the screen.
An embedding diagram turns under a finger that sets off across it, while a finger that sets off up or down scrolls the page, which `touch-action: pan-y` on the drawing leaves to the browser.
Nothing that starts on that drawing scrolls it sideways, so two fingers moving across it scroll it where it is wider than its frame, as at a large text size, and two taps bring back the published figure.

## The reader's text size

The prose of the spacetimes page and its headings are sized in rem, so a reader who enlarges text in the browser, by Chrome's font size, a larger default font, or Safari's and Firefox's zoom of text only, enlarges every word of it: the headings, the domains, the prose, the mathematics and every word on both kinds of diagram.
At the browser's usual 16px each is the pixel size it was when the page was set in pixels.
The spacetime's name, the word "metric" above it and every button keep one size in px at any text size, as the captain asked on 29 September 2026: "When I said I wanted the title and subtitle fonts to increase along with the text, I meant the font in the actual prose, not in the buttons."
Those are the choices of a chart, a view and an index placement, the print button, the names in the list and the reset of a figure that turns, each at the size it had before the prose grew with the text: the name 40px, the word above it 21px, a choice and the print button 15px and a name in the list 18px, and on a phone 26px, 14px and 12px.
The reset is drawn at the drawing's scale at the usual text size, as wide as its frame.
`FixedSizes` in `_tools/test_build_mfs_data.py` holds those sizes and holds the headings of the prose in rem.
The list's column grows with the text up to a quarter of the window, the panels start below the page's title however it wraps, the one-column layout's widths are in em so a large text size takes it on a narrower desktop window, and a sticky header that would cover more than a quarter of the view scrolls away instead.
A spacetime diagram's numbers and axis names are laid out around its plot, so they take the room they need and the numbers thin out where they would touch.
Every number, axis name and label of a spacetime diagram stands at least 1.5 em of the diagram's words from the edge of its ground, 21px at the usual text size, one em of the prose beside it, as the captain asked on 29 September 2026: "the spacetime diagrams need larger margins".
It is measured from the farthest a word can reach: half a number's height above the plot and half the last number's width past its right hand edge, where those numbers stand centred on the plot's edge.
The name of the vertical axis is a formula turned on its side, so `fitAxisNames` makes its column as wide as the name is tall, on the screen and at the size the print copy sets it.
A conformal diagram's labels and a figure's stand at points inside the drawing, so the drawing keeps its proportion to them and is drawn larger with the text, scrolling sideways in its own frame where it is wider than the panel.
`nrFigure` and `cdFrameAround` in `_layouts/mfs.html` carry the details.
No word is broken inside itself at any size: the name, the headings, the choices, the page's title, the names in the list and the search hint are held to the size at which their widest word fits their line, which `fitWords` and `mfsWidestWord` measure.
A word of the prose wider than its whole line is hyphenated by `hyphenateWords` with TeX's English patterns, since the browser's own hyphenation never divides a capitalised name such as Schwarzschild.
On a phone the page's title and its exit sign share a row until they no longer fit side by side, when the sign takes a row of its own below the title.
The print copy's sizes are points and its prose 12pt, so it prints the same whatever the reader's text size on the screen.
Reproduce a text size fault the way a reader meets it, with a larger default font size and not with page zoom, which scales everything and hides it.

Every choice on the page, the chart, a view of a diagram and a tensor's index placement, is one control, a row of buttons built by `choiceButton`, and a row appears only where there is more than one thing to choose.
Nothing on the spacetimes page glows, as the captain asked on 29 September 2026: "Get rid of the glow on the pressed button in MFS. Only turn it pink. I already said to get rid of the glow earlier. No glow anywhere."
The chosen button of a row, the spacetime shown in the list and any button while it is pressed turn pink, `--pink-light`, at once and with no transition, so a press as short as a tap shows it, and nothing else about them changes.
The panels are drawn by their borders alone, the title and the exit sign flicker in by colour and opacity alone, and the grid behind the page is drawn without a shadow.
No text shadow, box shadow, blurring or shadowing filter, canvas shadow or SVG filter appears anywhere on the page, open or not, on the screen or on paper.
MathJax's right-click menu brings a stylesheet whose menus cast a grey 20px shadow, and `mfsMenuWithoutShadow` in `_layouts/mfs.html` sets every such rule to none once MathJax is ready.
`NoGlow` in `_tools/test_build_mfs_data.py` holds the page's sources to that, and `page_timing.mjs` holds the page as drawn to it, including each kind of button pressed.
No text a reader of the spacetimes page sees says cost, costs, costly or costing, as the captain asked on 29 September 2026: "Don't use the word cost."
Each such sentence states the physical fact itself, as in "a shortcut that stays open requires exotic matter", and `NoCost` in `_tools/test_build_mfs_data.py` holds every string of every metric and diagram file and the page to that.
Every drawing, a spacetime diagram, a conformal diagram and a figure in three dimensions, stands on one dark ground, `--mfs-ground`, across the whole figure with its words, and prints on white.

Every written area of the spacetimes page is set in the history's font, Source Code Pro at its regular weight and upright, at the history's size, `--mfs-prose`, which is 21px at the usual text size and 15px in the phone layout.
That is the history and its tables, the conventions, every diagram's captions, legends and notes, the sentences under "not drawn", the restriction bands, the references, "vanishing" where a tensor or a scalar is zero, the placeholders and the loading line.
The note beside the coffee button is in that font at its own smaller size, since it stands in the list's column.
The headings, the choices, the list of names and every word inside a drawing keep their own fonts and sizes.
Two rules in `_layouts/mfs.html` name them all, and `WrittenAreas` in `_tools/test_build_mfs_data.py` holds every rule that reaches into one of them to that font and size, so a new kind of prose goes into both.
The page declares Source Code Pro itself, in `assets/css/mfs.css`, since it does not load the site's `custom.css`.
The print copy is set whole in EB Garamond at 12pt, the history with it.

The coordinates come first in the mathematics, right after the history, since the chart is chosen before anything that depends on it: their buttons, then each coordinate with its domain.
The page prints a spacetime's `signature` and its conventions under the heading "conventions" just below them, which is where the application reads them.
The conventions follow the chart, as the captain asked on 29 September 2026: each chart in `coordinates` carries its own `convention`, the spacetime's top level `convention` holds what is true in every chart and names no coordinate, and the page shows the chosen chart's text followed by the spacetime's as one paragraph, drawn again with the tensors whenever the chart changes.
A spacetime with one chart reads the same way, and the top level field is present in every file, empty where nothing is shared, since the application reads it.
A spacetime has one chosen chart, its first until the reader chooses another and kept while the page stays open, and the buttons, domains, conventions, tensors and diagrams all read it.
`node _tools/chart_conventions.mjs http://127.0.0.1:4000`, and again with `--phone`, walks every chart of every spacetime with more than one through each control and holds each of those to the chart chosen.

GitHub Pages lets a browser keep a file for ten minutes, so the page fetches each spacetime's files at the stamp its index entry gives them, and the index and `references.bib` at the time the site was built.
A page therefore reads only the files of its own build, and a file that did not change stays in the browser's cache across publishes.
A `convention` is one string of prose with inline TeX between dollar signs, and never carries a citation, a table, a `¶` or HTML.
An entry with neither field shows no conventions section at all.

A cloud-sync conflict copy dropped into the metrics folder, named like `kerr 2.json`, is passed over rather than read.
Those are the same names `.gitignore` already keeps out of the repository.
Any other file name is read, so a metric whose `id` does not match its file name still stops the command.

## How long a spacetime takes to open

`_tools/page_timing.mjs` opens every spacetime and every chart in headless Chrome as a reader does, and times each from the click until its mathematics is typeset and the page answers again:

    bundle exec jekyll serve
    node _tools/page_timing.mjs http://127.0.0.1:4000

It prints each chart's time and the longest task that held the page, slowest last, then every page error and console error the page raised, and exits non-zero when a chart takes longer than `--budget`, 10 seconds by default, or when the page raised any error.
`--metric <id>` times one spacetime, `--phone` lays the page out as an upright iPhone, `--text 48` sets the browser's default font size to 48px, as a reader who enlarges text does, and `--cpu 4` slows the processor four times over, as Chrome's own tools do to stand in for a slower phone.
It needs Chrome and Node 22 or later, and nothing installed.
Every change to the page is run at least as a desktop at the usual text and as a phone at three times it, with neither run raising an error:

    node _tools/page_timing.mjs http://127.0.0.1:4000
    node _tools/page_timing.mjs http://127.0.0.1:4000 --phone --text 48

Chrome hands "ResizeObserver loop completed with undelivered notifications" to the page's error handlers and not to its console, so the page's errors are gathered by a handler installed before its scripts run.
That error meant an observer's callback had changed the size of something observed, which the browser cannot report in the same frame.
Opening any spacetime at three times the usual text on a phone raised it until 29 September 2026, since `fitWords` held a heading to its line and `fitProseMath` set a formula on a line of its own from inside the observers that watch them, and so did choosing a view of a spacetime diagram there in WebKit, since `fitPlots` held the plot to its width from inside the observer of the plot.
Those observers now only note what changed and fit it at the next frame through `nextFrame`, and the page fits at once where it changes what is set itself, once MathJax has set a spacetime and when a diagram's view is chosen, so no frame is drawn before its fit.

The script also watches every number and label of every drawing shown, and each spacetime diagram's plot, from the first frame a reader sees them until two frames after the last font has loaded, as each spacetime opens and as each chart and each view of a diagram is chosen, and counts any that moves by more than a twentieth of a pixel, or is shown or hidden, as an error.
In each of those states it also holds every number, axis name and label of a spacetime diagram to that margin of 1.5 em, and counts any nearer the edge of the diagram's ground as an error; `--print` does the same for the print copy of each chart, the print dialog held closed.
A browser fetches a font only once it lays out text set in it, and MathJax draws each alphabet in a font of its own, so until 29 September 2026 the first spacetime to use one, Bianchi's fraktur or Malament-Hogarth's small capitals, fitted its numbers and slice labels in the face drawn in its place and fitted them again once the font arrived, and in Chrome they moved by up to 4 pixels as the panel slid in.
Headless WebKit drew no frame between the two fits, but it takes a running transition's time afresh at each measurement, so within one frame a label measured after its figure seemed to have moved by as far as the panel slid in between, about 3 pixels a millisecond; measure the page there with the panel's transition turned off.
`fontsLoading` in `_layouts/mfs.html` now loads every font a spacetime's mathematics and labels are set in, a hidden view's among them, before anything is fitted, and the panel slides in only then.

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

A noun is never handed back to what gives it: "every distance along it is the metric distance", never "is the distance the metric gives", and the same for "the X the Y yields", "provides" and "returns".
Nor is the metric said to give anything; an equation is stated as it stands, "where $f = 1$, $ds^2 = -c^2dt^2$".
The tests hold both rules over every prose field, every caption, label and declared input of every diagram, the templates and pages a reader sees, and the data in `/data` handed to agents.

## The captain's voice

Every paragraph of every history and convention, and every caption, note, restriction band and sentence stated in place of a drawing, is written in the captain's voice, which `~/VOICE.md` sets out in full.
A convention or a caption takes his register for formal physics: an active "we" where an agentless passive would stand, as "we take components in the chart $x^0 = ct$", an equation introduced as it stands, and the Oxford comma in every list of three or more.
A history takes his narrative register and names the people who did the work, with a full name at first mention.
None of them casts a paper, a result, a field or a drawing as the actor, stages a reveal, as "what carries the curvature is the Kretschmann scalar" did, closes on an epigram, explains what is not drawn, says what a thing is not in order to say what it is, or leans on an absolute for emphasis.
A figure of speech never stands in for the claim it points at: "so the spacetime does not stand on its own" became "so anti-de Sitter space is not globally hyperbolic", as the captain asked on 28 September 2026.

`VOICE` and the class `Voice` in `_tools/test_build_mfs_data.py` hold the rules a pattern can catch over every one of those texts, with mathematics and quotations taken out first, since a quotation is its speaker's own words.
They are tested to catch the phrases the captain flagged and to pass plain physics, so a new rule goes into the table together with an example of each kind.
The same class holds every hyphen to joining two names, a name and a word, or a designation, or to one of the few established terms it lists, so an English compound such as "future directed" is written without one.
Labels and legends name things rather than state them, and stand outside these rules.
`_tools/derivations/voice_rewrites.md` records the rewrite of 28 September 2026 of the histories and captions text by text, before and after.
A rewording is checked against the text it replaces with

    python3 _tools/prose_math_check.py [revision]

which reads every history, convention and diagram text from the working tree and from `revision`, `main` by default, and exits non-zero where one's mathematics, numbers written in digits or citations differ.

## The shape of a history and a convention

A history has at least five paragraphs, and every paragraph has three to six sentences, so no paragraph runs to more than twice the length of another.
A table, the paragraph written as a `TABLE::` line, is not prose and is left out of the count.
Where an entry is too thin for five paragraphs, it wants more history, sourced and cited like the rest, never filler.

A convention says only what a reader needs to read the mathematics, with Schwarzschild's as the model: "We use coordinates $(t, r, \theta, \phi)$, with $t$ carrying dimensions of time. The Schwarzschild radius is $r_s = 2GM/c^2$. We keep factors of $c$ explicit."
That is the coordinates and their units, what each parameter and function is in a clause, the factors of $c$ and $G$, and any choice of chart, index, sign or normalisation the components depend on, such as $x^0 = ct$ or the Riemann and Ricci conventions.
The domains stand just above it, so it does not repeat the ranges.
It never describes what the tensors or diagrams show, explains a standard object such as the Weyl tensor, or carries physics, history, energy conditions or counts of components, and it never uses the word "cost"; the captain cut every convention down to that on 29 September 2026.
Each chart's text and the shared text together make one paragraph of at most five sentences.
The tests hold every break of a history to the end of a sentence, with no space on either side of the `¶`.

`sentences` in `build_mfs_data.py` does the counting.
A sentence ends at a full stop, a question mark or an exclamation mark, after any closing quote or bracket, where the next word begins with a capital or a digit.
A capital standing alone before a full stop is an initial, as in J. Robert Oppenheimer, and ends nothing.
Mathematics between dollar signs counts as one word, and ends a sentence only when the stop is inside it, as when a displayed equation closes one.

The command refuses to write anything, and `--check` fails, while any history is out of shape, naming each with the sentences in each of its paragraphs and the paragraph that breaks the rule.
It refuses in the same way a chart with no convention of its own, a spacetime with no shared `convention` field, and a chart whose convention, shared text included, runs past five sentences, runs to a second paragraph, or says "cost", "slots", "nowhere does the geometry break", "standing objection" or "the Riemann tensor with the traces removed", naming the chart at fault.

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
The page draws them in a section per coordinate system, the first of the three diagrams that close the mathematics, after the geodesics, and the application displays the same files.
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

A caption opens with a noun phrase naming its plane: the two coordinates drawn, with the value of every coordinate held fixed in parentheses, as "The plane of $t$ and $r$ ($\theta = \pi/2$, $\phi = 0$)".
It never opens "This is", and nothing in it "stands for" anything; `CAPTION_VOICE` in `_tools/test_build_mfs_data.py` holds every caption, note and restriction band on the site to that.
It says what the drawing shows, and where the feature a reader comes looking for lies off the plane, it says where that feature is, as the Ellis-Bronnikov caption places the throat in $g_{\theta\theta}$ and the Gödel caption places the closed timelike curves.
A caption never stops at saying what a diagram leaves out.
A view of rays that leave the plane names the surface they lie on and the coordinate left out of the drawing, as "The equatorial plane ($\theta = \pi/2$) drawn in $t$ and $r$, with $\phi$ left out".
Every curve drawn is a null curve; a caption calls it a null geodesic, the path light takes, only where no Christoffel symbol turns it out of the plane, and says so where one does, as for Gödel and for Kerr off the axis.

A caption is prose in the register of the conventions and histories beside it, and it is read back whole, in the generator and on the page, since a sentence built from clauses that each keep the rules above can still read as machine prose.
It says what a thing is and does, plainly, as "each point in the diagram a single event", and never defines it by a contrast with what it is not, as the captain asked on 29 September 2026: no "X rather than Y", "instead of", "not a X but Y", "is X, not Y" or "less X, more Y" anywhere a reader of the spacetimes page sees, in a history, a convention, a description, a caption, a note, a restriction band or the page's own words.
`Contrast` in `_tools/test_build_mfs_data.py` holds all of them to that, and a sentence where such a phrase states physics plainly may stand only as a named entry of its `ALLOWED`, with the reason.
It gives the physics of a surface, and the reader sees how the figure lays it out; it ends on a fact, so the lines inside Vaidya's shell crowd toward it because of the slope of $q$.

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
`outgoing` does the same with the outgoing rays toward larger $r$, for de Sitter's static chart and the Krasnikov tube, and in Misner's plane of $T$ and $\psi$, where $\partial_T$ is null and no chart time is timelike, with the family that tips over, toward larger $\psi$.
`split` is `ingoing` below the radius the row names in `split` and `outgoing` above it, for a static region between a black hole horizon and a cosmological one, where both agree with $t$: Schwarzschild-de Sitter's static chart uses it at its static radius, which reads the region inside $r_h$ as the black hole and the region beyond $r_c$ as the expanding universe.
`vector` takes the chart's time direction, for Godel, whose published $g^{tt}$ is positive and which has no time function at all.

### Declared inputs

An entry that leaves a function free cannot be drawn without a choice, and every such choice is written in its row and printed beside its diagram.
FRW's scale factor and Bianchi I's three are dust, solved from each entry's own published $G^i{}_i = 0$.
Morris-Thorne is drawn as its Ellis-Bronnikov member, $\Phi = 0$ and $b = b_0^2/r$.
Alcubierre is drawn with his own profile, which the entry names, at $v_s = 2$; Natario with the same profile and his zero expansion field, whose divergence the script checks before using it; the Krasnikov tube with a tube built by a ship at the speed of light.
Kasner is drawn at the exponents $(-2/7, 3/7, 6/7)$.
Misner space is drawn at Li and Gott's $\psi_0 = 4\pi$, a boost of rapidity $2\pi$, in each of its three charts, every plane drawn unrolled with its periodic edges one line; Misner's own plane marks $g^{TT} = 0$, the chronology horizon.
Vaidya is drawn with a shell of null dust of mass $M$, $m$ jumping from $0$ to $M$ at $v = 0$ in the ingoing chart and its time reverse in the outgoing one, the shell the conformal diagram declares.
The Tolman-Oppenheimer-Volkoff star is the polytrope the conformal diagram declares, solved by the `StarSolver` the two scripts share and checked against the published $G^\theta{}_\theta$, which its construction never uses.
Tolman-Bondi is drawn for marginally bound dust whose density at $t = 0$ falls as $1 - r^2/r_b^2$, checked to solve the published $G^r{}_r = 0$, and Oppenheimer-Snyder's collapse is released from rest at twice its Schwarzschild radius, as its conformal diagram declares, in both of its charts.
The Malament-Hogarth toy is drawn for every conformal factor at once, since $\Omega^2$ drops out of the null condition, and checked with a declared factor that grows as $1/|ct|$ along the axis.
The cosmic string is drawn at the deficit its conformal diagram uses, $4G\mu/c^2 = 0.1$, Gödel at $\omega = 1$, and Kerr and Kerr-Newman's figures at the spins and charge of their flat views.
The C-metric is drawn at $\alpha m = 1/6$ and $C = 3/4$, on the two halves of its axis, where $\sin\theta$ kills $\Gamma^\theta{}_{tt}$ and $\Gamma^\theta{}_{rr}$; beyond each horizon its time function is $r$, piecewise, taking the future as Griffiths, Krtouš and Podolský's extensions do, the black hole inside $2m$ and the region to the future of the acceleration horizon beyond $1/\alpha$.
The Einstein-Rosen waves are drawn with the pulse of Weber, Wheeler, and Bonnor at $C = a$, checked to solve the published $G^t{}_t$, $G^t{}_\rho$, $G^\phi{}_\phi$ and $G^z{}_z = 0$; on the plane of $t$ and $\rho$ the metric is conformally flat, so the rays are the same for every wave.

### Checking the rays

    /tmp/mfs-venv/bin/python _tools/derivations/null_rays.py --verify

traces rays as the page's files are traced and measures, for every view with a closed form, how far the quantity each family should conserve drifts along a ray, together with the dust solutions and the equality of Natario's and Alcubierre's metrics on the plane of their axis.
For the principal null rays it adds the checks the next section describes, and the closed forms $ct \mp r_*$ and $\phi \mp r_\sharp$ of Kerr and Kerr-Newman, which the drawing never uses.
On each cylinder of $t$ and $\phi$, van Stockum's and Gödel's, it checks the straight lines against the slopes the captions state, and checks what light launched along each family does against what the caption says, from the published Christoffel symbols: stays on the cylinder as a null geodesic, or is turned toward the axis or away from it.
It writes nothing and exits non-zero if anything fails; run it after changing the method.
`--metric <metric_id>` checks only that spacetime's closed forms and dust, and leaves out the checks of turning and of the principal null rays.
It took 570 seconds on 29 September 2026.

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

### Rays of no angular momentum

The rotating BTZ hole turns its light rays in $\phi$ as Kerr does, and between its ergosurface and its outer horizon the plane of $t$ and $r$ at fixed $\phi$ has no null direction; its Weyl tensor vanishes, as every Weyl tensor in three dimensions does, so it has no principal directions either.
A row with `quotient="phi"` divides the circles of $\phi$ out: the plane's metric is $h_{ab} = g_{ab} - g_{a\phi}g_{b\phi}/g_{\phi\phi}$, the part of the published metric orthogonal to the circles, and it takes the place of the coordinate plane's metric in the null condition, as the principal plane's does for Kerr.
Its rays are the shadows on $t$ and $r$ of the null geodesics with no angular momentum, each lifted into the spacetime by $d\phi = -(g_{\phi t}dt + g_{\phi r}dr)/g_{\phi\phi}$, and the upper left block of the published inverse metric is the inverse of $h$, so the horizons marked from $g^{rr}$ are the zeros of $h$'s own.
The script refuses such a row unless every coordinate is drawn, held fixed or divided out and nothing published depends on the one divided out, and before it writes the view it checks every lifted direction null against the published metric and geodesic by the published Christoffel symbols, which it stamps with the view; the header's "Rays of no angular momentum" gives the details.
The BTZ hole's three charts each draw the hole without rotation on the plane $\phi = 0$ and the rotating one with $\phi$ divided out, at $M = 1$ and $J = 4\ell/5$, where $r_\pm^2 = 4\ell^2/5$ and $\ell^2/5$.

### Figures in three dimensions

Where the causal structure a reader comes for turns in a direction no plane of two coordinates holds, a coordinate system also gets a figure in three dimensions: a slice of one time and two spatial coordinates, every other coordinate held fixed, drawn from the published metric by `_tools/derivations/projections.py` and projected once from a fixed camera.
The file holds the projection, polylines, polygons and points of the page's own plane in the order they are painted, farthest first, with TeX labels, which is the form a conformal diagram takes.
So the page draws a figure with the conformal diagram's frame and labels, it prints, and an application can draw the same data with its own renderer.
A figure of light cones also carries its pieces in three dimensions under `turn`, from which the page draws it again from whatever side the reader turns it to.
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

### Turning a figure of light cones

The figures of Gödel, van Stockum, Kerr, Kerr-Newman and Alcubierre carry `turn`, every piece of the figure in the drawing's $(X, Y, T)$, $T$ up, at six decimals, and the page turns them by the same hand as an embedding diagram.
The cosmic string's beam lies in the plane $T = 0$ seen from straight above, in the plane's own flat coordinates, so it has no other side and no `turn`.
`turn` holds:

- `lines`: every line the figure draws before its cones, in the order painted, each `{"class", "points"}`, a closed curve ending on its first point.
- `cones`: every future light cone, each `{"class", "apex", "rim"}`, its apex and the ends of its generators, all of one Euclidean length, in the order the figure paints them at its own camera, farthest first.
- `ribs`: how many of each cone's generators are drawn from the apex to the rim, evenly round it; every rim's length is a multiple of it.
- `labels`: where each of `labels` stands, in order: `{"at": [X, Y, T]}`, a point it keeps from every side, or `{"circle": [rho, T], "angle": delta}`, the circle of radius $\rho$ about the axis at height $T$ that it names, at the point of that circle $\delta$ degrees round from the camera's azimuth, counterclockwise seen from above, so it keeps to the side of the circle nearest the reader.
- `slices`: each of `slices` in three dimensions, in order, its `lines` and `fills` as the published slice has them.
- `centre`: the point the figure turns about, on the axis halfway between the lowest and the highest point it draws.

At another camera, azimuth $a$ and elevation $e$, with the vectors `right`, `up` and `toward` of `projections.Camera`:

- Every line is projected and thinned by Ramer-Douglas-Peucker to 0.0005 of the box's units.
- Every cone is painted after the lines, farthest first by `toward` at its apex, two cones whose apexes lie within $10^{-5}$ of the box's larger side of each other in depth keeping their order in `cones`.
  A cone is the convex hull of its projected apex and rim, each point taken to twelve decimals, filled as `cone`; its rim, closed, as `cone-rim`; every generator from the first round to the last a `ribs`th of the rim apart, from the apex, as `cone-rib`; where the apex is a corner of the hull, the two edges of the hull that meet there, as `cone`; and its apex, as the point `cone-apex`.
- Each slice's lines and rings are projected and thinned as a line is, and drawn under everything else, as the published slice is.
- A label stands at its point, or on its circle at $\phi = a + \delta$, with its anchor and offset as published; a label on a circle gives way while it would overlap a label before it, as a label with `ring` does on an embedding diagram.
- The figure turns about `centre`, which stays where the figure's own camera puts it on the page.
  It is drawn at one scale $s$, never above 1, the largest that keeps it within the height it had at its own camera and within the width it had on each side of its axis, and moved up or down only as far as that height needs; the `box` stays as it is, and every point drawn stays inside it.

At the figure's own camera these rules give back the published figure, the same layers in the same order, every point of a cone and every label to the published rounding and every line within its thinning.
`_tools/turn_check.cjs` holds the page's drawing of every figure to that, and to its box from every side.

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
The page draws it under the heading "conformal diagram", just below the spacetime diagram, and it follows the chart chosen at the top: it shows the views that name that chart, then the views that name no chart, with buttons only when that makes more than one.
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

`DRAWN` in the script names the twenty-four spacetimes that have a diagram and `NOT_DRAWN` the twelve that have none: every event's future is the whole spacetime for Gödel and van Stockum, the diagram is whatever a free function makes it for the warp drives, the Krasnikov tube, Tolman-Bondi and Bianchi, the Mixmaster has no surface that carries its causal structure, and so on.
Misner space is drawn in the plane $y = z = 0$ of the Minkowski space that covers it, one view per chart, each tinting one copy of the spacetime between two lines the boost carries onto one another; its views are drawn at $\psi_0 = 2$, since at the $\psi_0 = 4\pi$ of its other diagrams one copy stretches by $e^{2\pi}$ along each light ray and fills the drawing to within a pixel.
A free function does not by itself rule a diagram out: where every choice of it gives the same shape, as for the Tolman-Oppenheimer-Volkoff star and the Morris-Thorne wormhole, the diagram is drawn with a declared choice that moves only the lines inside, and the Malament-Hogarth toy is drawn for every conformal factor at once, since a conformal factor changes no null direction.
Those twelve have no file, which a full redraw removes if one is left behind, so the page has no conformal diagram section for them; `build_mfs_data.py` refuses a file that draws nothing, and the script stops if a metric file is in neither table, so a new spacetime needs a decision.

The Einstein-Rosen waves are drawn on the half plane of fixed $\phi$ and $z$, conformal to Minkowski's half diamond for every wave, since only the factor $e^{2(\gamma - \psi)}$ stands in front of $-c^2dt^2 + d\rho^2$; the maps are checked with random values of $\psi$ and $\gamma$ at every sample and with the pulse of Weber, Wheeler, and Bonnor, whose published Kretschmann scalar settles on the axis, at $\rho = 10^{-3}\,a$ and $10^{-4}\,a$ the same to a part in $10^4$ of its largest value, near $5.7\times10^5/a^4$ at $t = 0$.
The BTZ black hole's $1/N^2$ is a sum of simple poles with no constant term, so `BTZTower` sums its tortoise coordinate over every root of $N^2$, the negative ones included, which makes $r_*$ vanish both at $r = 0$ and as $r \to \infty$ and puts the conformal boundary and $r = 0$ on the vertical lines $X = \pm\pi/2$.
Without rotation that is the square Bañados, Henneaux, Teitelboim and Zanelli drew, and the rotating hole is drawn on its plane of $t$ and $r$ with $\phi$ divided out, `Plane(..., quotient="phi")`, as its spacetime diagram is, a tower between those lines cut to one period by `clip_in_t`.

A view of a surface that is not the whole spacetime carries `restriction`, which the page prints in a band across the top of the figure, never in a footnote.
Kerr and Kerr-Newman are drawn on the symmetry axis and the cosmic string on the half plane of fixed $\phi$ and $z$, and the tests hold those three to carrying a restriction on every view.

Every text in a view is TeX in `$...$`: the labels on the drawing, the buttons, the legend, the caption, the restriction and the parameter values.
A caption opens with a compact noun phrase naming what is drawn, its values in parentheses, in the captain's model: "A spherically symmetric distribution of dust collapsing from rest ($R_0 = 2\,r_s$), each point in the diagram a 2-sphere."
Its prose is for the reader of "Whom the prose is for" above; the tests hold every text in these files to that rule's words and to the dash rule.
It keeps the register the spacetime diagrams' captions keep, and it names what each point in the diagram is in that same construction: "each point in the diagram a 2-sphere of radius $r$" for a spherical spacetime, "each point in the diagram a single event rather than a 2-sphere of them" on a slice such as Minkowski's plane $y = z = 0$, and "each point in the diagram a circle around the string times a line along it" on the cosmic string's half plane.
No caption says "This is the whole of" or "stands for".

`_tools/derivations/tex_check.cjs` typesets every TeX string the page sets, in the metrics, the diagram files and the conformal files, a published value together with its negation, through the TeX input MathJax loads on the page, and exits non-zero naming each one it cannot set.
sympy never reads the typesetting, so this is the check that catches a value that is right and prints as an error box:

    npm install --prefix /tmp/mfs-node mathjax-full
    node _tools/derivations/tex_check.cjs /tmp/mfs-node

On 29 September 2026 it set all 24210 strings.
It does not read the embedding files yet, so a caption, setting or label of an embedding diagram that MathJax cannot set shows only on the page.

## Embedding diagrams

`MFS/assets/data/embedding/<metric_id>.json` holds the embedding diagram of one spacetime: a slice of it, the equatorial plane at one moment or another surface of one moment, drawn as a surface in ordinary flat three dimensional space so that every distance along the surface is the metric distance, or, where no surface in flat space carries it, in three dimensional Minkowski space.
Every spacetime has one but Lentz's soliton.
Every slice of constant $t$ in Lentz's class is flat for every potential $\phi$, and his soliton exists only as a numerical integral over rhomboid sources whose sizes, charges and places his paper leaves to a figure, so nothing of the soliton can be computed to draw on the plane, and the page shows no embedding diagram for it.
Three draw a quantity as a height over a plane in place of a slice, since their slices are flat and what makes each spacetime what it is shows in a number on the plane: Alcubierre's warp drive draws the expansion $\theta$ of the observers who ride its slices, Natário's the energy density $\varepsilon$ those observers measure, and the Krasnikov tube $1 - k$, how far it tips the light cone.
No height is a distance, and each of the three views says under `height` what its height stands for.
The page draws it under the heading "embedding diagram", just below the conformal diagram, and the application reads the same files.
The file is the definition in "The file, which the application reads" below, and the application is built against that section, so a change to the shape of the file is a change to it first.

`_tools/derivations/embedding.py` draws them, one function per spacetime with its derivation in its docstring, reading the published metric through the checker's `Reader` with the `load` and `published_matrix` the null rays use.
It needs the same environment as `null_rays.py` and took 409 seconds for the whole collection on 30 September 2026, most of it in the frames of the movies, above all Mixmaster's, whose every frame rewrites its metric in half angles in sympy, and the Malament-Hogarth plane's, whose conformal factor is a piecewise function evaluated whole at every point of every quadrature:

    /tmp/mfs-venv/bin/python _tools/derivations/embedding.py
    python3 _tools/build_mfs_data.py

`--metric <metric_id>` redraws one spacetime, and `--verify` prints every check and writes nothing.

### The construction

Every slice drawn, which is every surface but a height over a plane, is the surface of one spatial coordinate $x$ and the angle $\phi$ of one coordinate system, every other coordinate held fixed, with the metric $g_{xx}(x)\,dx^2 + g_{\phi\phi}(x)\,d\phi^2$ read from the published `metric_components`.
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

A slice may move one coordinate with the angle, as the Mixmaster great sphere's $\psi = -\phi$, where `Slice` pulls the published metric back along the move as it pulls Vaidya's back along $v = T + r$, and the result must depend on $x$ alone.
A plane whose angle is no coordinate of its chart, as the Malament-Hogarth plane of $x$ and $y$, is turned about the chart's origin: the metric is read on the profile, where the angle is zero, and the turn is checked at a thousand points of the plane to carry the pulled back metric onto itself to a part in $10^{12}$, since sympy seldom clears $\cos^2 + \sin^2$ inside a declared function.
A plane of two chart coordinates on which the published metric has constant coefficients and no cross term, as every slice of a warp drive, of Kasner's universe or of Bianchi type I, and a pp-wave's wave front, is `FlatPlane`: in the proper coordinates $\sqrt{g_{xx}}\,x$ and $\sqrt{g_{yy}}\,y$ it is Euclid's plane, drawn as a flat disc about the chart's origin with its profile along $x$.
The Mixmaster great sphere's $g_{\phi\phi}$ and $g_{\theta\theta} - (d\rho/d\theta)^2$ both vanish as $\theta^2$ at its poles, and sympy's forms of them subtract numbers that agree there, so `half_angles` writes each as a factored ratio of polynomials in $\sin^2(\theta/2)$, which keeps every digit.

Where every circle grows faster than the distance out to it, as on anti-de Sitter space's static slice, the slice is drawn in three dimensional Minkowski space, $dX^2 + dY^2 - dZ^2$.
There a profile turned about the $Z$ axis has the metric $\left((d\rho/dx)^2 - (dZ/dx)^2\right)dx^2 + \rho^2d\phi^2$, so it climbs at $dZ/dx = \sqrt{(d\rho/dx)^2 - g_{xx}}$, and the whole hyperbolic plane lies on one sheet of a hyperboloid.
Near its light cone a chord of such a surface is short in the metric against its length on the page, so it is written to below $10^{-9}$ of its extent, nine decimals where the others take seven.

A quantity drawn as a height over a plane is a `GridPiece`: the height $z(X, Y)$ at the nodes of a grid of the plane, polar or Cartesian, each cell cut into two flat triangles along its diagonal from $(i, j)$ to $(i + 1, j + 1)$.
The grid holds the values of $u$ and $v$ the drawing needs, as the lines it draws and the curves it marks along them, and every interval of $u$, and of a Cartesian grid's $v$ as well, whose cells hold a triangle further than $5 \times 10^{-4}$ of the drawing's size from the height, measured as the distance in space, is halved until none does.
That distance is the vertical miss divided by $\sqrt{1 + |\nabla z|^2}$: on a steep wall a triangle close to the height misses it vertically by several times as much.
Halving $u$ cannot close the gap a polar grid's straight chords leave where its angles are too few, and the script stops with that reason once a grid passes 200000 nodes rather than halve for ever.
The heights are rounded as they are written, and every triangle, line and check is taken from the rounded heights, so the surface drawn is exactly the surface in the file.
A curve marked on a grid is a run of its nodes or a level line of its triangles, found by marching over them, so every point of it lies on a triangle.

### Checking the surface

The surface is measured as the application will draw it, from the rounded numbers the file holds, against the published metric, and the script refuses to write while any check fails:

- along: each chord of each profile, and each profile end to end, against the proper distance between the same two values of $x$, a chord of a surface in Minkowski space measured as $\sqrt{d\rho^2 - dZ^2}$;
- across: the straight line in space from each point to the next one $0.02$ further round the axis, against the metric length of the line that runs out at a steady proper distance while it turns steadily through the same angle;
- around: $\rho$ at every point against $\sqrt{g_{\phi\phi}}$;
- joins: where two pieces meet, as a star's surface meets the exterior, they meet at one point with one tangent, which says $g_{xx}$ agrees on both sides;
- forms: each surface against the closed form it is known by: Flamm's paraboloid, also as the exterior of the neutron star, of each collapse at its release and, moved in by $r_s$, of Vaidya's slices; the interior Schwarzschild cap, the catenoid in both its charts, the cone, Gott's cap, the spheres of FRW, de Sitter, Oppenheimer-Snyder's dust and Taub's round moment, Bertotti-Robinson's cylinder and sphere, the planes of Minkowski space and of Natário's and Lentz's drives, and anti-de Sitter's hyperboloid; the circumference radius of the Kerr, Kerr-Newman and Taub-NUT horizons and of Taub's equator; and Kasner's, Bianchi's and the pp-wave's rings against the ellipses they are;
- heights: every node of a grid against the quantity it stands for, from the published metric and the declared functions, and every triangle within $5 \times 10^{-4}$ of the drawing's size of it in space: Alcubierre's $\theta = \nabla_\mu n^\mu$ against $v_s\,(x - x_s)/r_s\,df/dr_s$, negative ahead of the ship and positive behind it, with its level lines at half the greatest expansion and contraction and its slice checked flat; the Krasnikov tube's $1 - k$ against twice the published $g_{tx}$, with both null directions of its plane of $t$ and $x$ checked null for every $k$ and its level line $k = 0$ where the published $g_{xx}$ vanishes;
- fields: a declared star, scale factor or dust cloud against the published Einstein tensor, the neutron star against $G^\theta{}_\theta = 8\pi p$, FRW's and Oppenheimer-Snyder's scale factors against their $G^r{}_r$ and $G^\chi{}_\chi$, each cloud released from rest against Tolman-Bondi's $G^r{}_r = 0$ and its $G^t{}_t$, the density, and Taub's scale factors against every published Einstein component of the Mixmaster universe; the pp-wave's particles are run with its published Christoffel symbols, Kasner's and Bianchi's are held at rest by theirs, and the length of the Malament-Hogarth tube is checked against the proper time of the computer on its axis;
- marks: every curve and point marked on a surface lies on its piece, at the height the profile has at its distance from the axis or on the triangles of a grid, and each of Natário's lines of flow closes and holds Stokes's stream function to one value;
- stops: where a view or a file says a slice cannot be drawn, $g_{xx} - (d\rho/dx)^2$ is negative at every sample, or $g_{\phi\phi}$ is, where the circles are timelike; where it says a slice is a plane, it is zero, or the spatial metric has no cross term and no component that depends on a spatial coordinate.

On 28 September 2026 the worst chord missed its proper distance by $1.0 \times 10^{-4}$ of it, a chord $9 \times 10^{-4}$ long on the lip of the Malament-Hogarth well at $ct = -0.3$, where the rounding is most of it, and the worst line across by $1.1 \times 10^{-4}$, on the sphere of radius $a = 1/2$ of FRW's first moment, for the same reason.
Every profile end to end was within $3.6 \times 10^{-5}$ of its proper length, every $\rho$ within $2.5 \times 10^{-8}$ of the drawing's size of $\sqrt{g_{\phi\phi}}$, every closed form within $2.5 \times 10^{-8}$ of the size, every join met to $3 \times 10^{-13}$ with tangents equal to $3 \times 10^{-14}$, and every declared solution satisfied the published field equations to $1.4 \times 10^{-14}$, Taub's exactly; the triangles of Alcubierre's and the Krasnikov tube's grids lay within $4.9 \times 10^{-4}$ and $4.2 \times 10^{-4}$ of the size of their heights, in space; 645 checks in all, none failing, and running the script again writes the same bytes.
On 30 September 2026, with every frame of the movies checked as a moment is, it made 3854 checks, none failing, and the worst chord missed its proper distance by $1.4 \times 10^{-4}$ of it.
The tests hold the files on disk to the same closed forms without sympy, from the numbers written and nothing else, so a redraw that changed a surface fails there too.

### Which spacetimes

`DRAWN` in the script names every spacetime, and `STATED` and `NOT_DRAWN` are empty: a spacetime put in `STATED` would have a file with no surface that says why, checked from the published metric before it is said, and one in `NOT_DRAWN` no file, and the script stops if a metric file is in none of the three or in more than one.
Minkowski space is its equator from the spherical chart, the plane every other diagram is measured against, out to $r = 4\ell$ in any length $\ell$, and every slice of constant $t$ of its Cartesian chart is checked flat.
Misner space is its contracting region at $ct = -2$, $-1.5$, $-1$ and $-0.5$ in the Milne chart at $\psi_0 = 4\pi$, where $\chi$ runs once round $2\pi$: each moment's slice $z = 0$ is a flat cylinder of radius $c|t|$, checked against $\rho = c|t|$ and $z = y$, and its moments lie on the Misner and Milne planes and on none of the Rindler region's drawings, which `HIDDEN` records.
Anti-de Sitter's static equator has $g_{rr} < (d\rho/dr)^2$ at every $r > 0$, which is checked, so it is drawn in Minkowski space: the sheet $Z = \sqrt{L^2 + r^2} - L$ of a hyperboloid out to $r = 4L$, with the light cone it nears drawn as a reference piece, $Z = \rho - L$, given by its closed form since it is no slice of the spacetime.
The Malament-Hogarth plane $z = 0$, with the $\Omega$ its spacetime diagram declares, is turned about the removed event at $ct = -0.7$, $-0.3$, $-0.1$ and $0$, and played as a movie with a frame every $1/80$ of $ct$ up to the last whose well is no deeper than the tube at $ct = 0$, since the wells just before it reach deeper than the tube is drawn, the flat plane at $z = 0$ in every frame: a well inside the unit ball, flat beyond, and at $ct = 0$ a funnel into a tube of radius $1$ drawn down to $s = 0.03$, whose length from the rim down to $s$ is checked to be the proper time of the computer on the axis from $ct = -1$ to $-s$; the edge of the region where $\Omega > 1$ is $\sqrt{1 - c^2t^2}$ to the last digit, moved on by a unit in the last place where the declared function's two branches disagree within one of it, as they do at $ct = -0.7$.
The Mixmaster universe is the great two sphere of its three sphere, $\psi = -\phi$ and $\psi = 2\pi - \phi$ in the Euler angles, at five moments of Abraham Taub's universe, played as a movie with a frame about every $0.05\,m$ of $c\tau$, each at the $T$ whose proper time it is, $a_1 = a_2$, at $m = 1$ and $l = m/2$ as Taub-NUT declares: only there does the great sphere have an axis, through the identity, and where $a_3 > 2a_1/\sqrt{3}$ its curvature about the poles is negative and the drawing is the band about the equator, which the view states and the script checks.
Kasner's universe, Bianchi type I and the pp-wave are each a flat plane at four moments, Kasner's and Bianchi's the plane $y = 0$ at the exponents and the dust their spacetime diagrams declare, the pp-wave's its wave front at $cu = -3$, $-0.5$, $0$ and $0.661\,L$ with the pulse $A = e^{-u^2}$ its spacetime diagram declares, each with a ring of free particles marked on it as a curve and twelve of them as points: at rest in the chart for Kasner and Bianchi, and run from rest before the pulse for the pp-wave, which focuses every one onto the $x$ axis at once.
Each is drawn first as the ring's world tube, as the captain asked on 30 September 2026, and the flat moments are the view beside it: the ring at every time from the first moment to the last, each ellipse at the height of its time, $1.5\,\ell$ for a unit of Kasner's $t$, $2.5\,\ell$ for $1/\bar H$ of Bianchi's and $L/2$ for each $L$ of the wave's $cu$, as the view's `height` says.
`stack()` in the script builds it as a stack of ellipses whose columns are the ring's particles, one a degree, and checks every row against the ring at its time, from the published metric for Kasner and Bianchi and from the particles run with the published Christoffel symbols for the wave, and the tube's rows at the four moments against the flat views' rings.
The Krasnikov tube is $1 - k$, twice the published $g_{tx}$, as a height over the plane of its axis at $ct = 5$, a unit after the ship reached the far end of the declared tube, when the whole tube stands as a ridge along the path: over $-1 \le x \le 5$, as its spacetime diagram runs, and out to $r = 2$ on either side, on a Cartesian grid with its lines every $\rho_0/2$, with the path from $x = 0$ to $D$ and the level line $k = 0$, where the published $g_{xx}$ vanishes, marked as `wall`; `_tools/derivations/krasnikov_height.md` is the derivation.
Alcubierre's drive is the expansion $\theta = v_s\,\partial_xf$ of the observers who ride its flat slices, read from the published metric, as the height $\theta R^2/4c$ over the plane $z = 0$ of the path at $t = 0$, the ship heading along the drawing's $X$, on a polar grid of 72 angles with its circles every $R/4$: marked on it, the path, the circle $v_sf = 1$, the zero of the published $g_{tt}$, as `wall`, and the level lines where $\theta$ is half its greatest value, `contract` ahead of the ship and `expand` behind it; `_tools/derivations/alcubierre_expansion.md` is the derivation.
Natário's drive is the energy density $\varepsilon = -c^4K_{ij}K^{ij}/16\pi G$ of the observers who ride its flat slices, read from the published $G^{tt}$ with the declared field and checked against $-K_{ij}K^{ij}/2$ from the published shift, whose trace, the expansion, is checked to vanish, as the height $G^{tt}R^2/16$ over the plane $z = 0$ of the path at $t = 0$, the ship heading along the drawing's $X$ as Alcubierre's does, on a polar grid of 72 angles with its circles every $R/4$: a moat in the wall, deepest beside the ship at $r_s = 0.88R$, and level at the ship and far outside.
Marked on it, the path, the circle $r_s = R$ as `wall`, and six lines of the declared flow, each followed from the ship's plane $x = 0$ round to it again and set on the triangles of the grid.
Schwarzschild is Flamm's paraboloid on both sheets of the Einstein-Rosen bridge, the slice of constant $t$ running through the bifurcation sphere into the other exterior.
The interior Schwarzschild star, at $R = 1.5\,r_s$ as its conformal diagram declares, is the cap of a sphere of radius $\sqrt{R^3/r_s}$ joined to the exterior of `schwarzschild.json`, which the file reads, with the vacuum paraboloid drawn on under the cap down to the throat the star replaces.
The Tolman-Oppenheimer-Volkoff star is the polytrope its other diagrams declare, $M = 1.40\,M_\odot$ and $R = 14.2$ km, with the mass function `null_rays.StarSolver` solves from the published Einstein tensor entering the slice as numbers; its exterior is checked to be Flamm's paraboloid of that mass, and the vacuum paraboloid is drawn under it as under Schwarzschild's star.
Morris-Thorne is drawn as its Ellis-Bronnikov member, $b = b_0^2/r$, the catenoid, as its other diagrams declare; the surface reads $b$ alone, and its proper radial chart, with $r(l) = \sqrt{l^2 + b_0^2}$, is checked to give the same surface.
Reissner-Nordstrom, at $r_q = 0.48\,r_s$ as its conformal diagram declares, has two views: outside $r_+$ the two sheets through the outer horizon's bifurcation sphere, and inside $r_-$ the two sides through the inner one, its widest circle there, each running in until $g_{rr} = 1$ at $r = r_q^2/r_s$, where the surface lies level; nearer the singularity and between the horizons no surface carries the slice, which the views state under `stops` and the script checks.
The BTZ black hole is its moment $t = 0$ without rotation, $M = 1$ and $J = 0$, as its conformal diagram's square draws it, in one view through the throat $r_+ = \ell$ into both exteriors out to $r = 3\ell$: in flat space out to $r = \sqrt{2}\,\ell$, where $N^2 = 1$ and the surface lies level, and beyond that in Minkowski space, where it climbs from level toward a light cone as anti-de Sitter space's hyperboloid does.
Each exterior is two pieces, which meet in the ring of class `space`; the script checks that flat space carries no part beyond $\sqrt{2}\,\ell$ and Minkowski space none inside it, that the two parts of each exterior meet in one circle with one tangent, in the numbers and in the file, and that the two exteriors meet at the throat.
The rotating hole's diagrams list the moment in `slices.HIDDEN`, since their hole is another spacetime.
Kerr, at $a = 0.9\,GM/c^2$, and Kerr-Newman, at $a = 0.6\,GM/c^2$ and $r_Q = 0.5\,GM/c^2$, as their conformal diagrams declare, are the equator of a slice of constant Boyer-Lindquist $t$, where $g_{t\phi}$ drops out, drawn at the circumference radius $\rho = \sqrt{g_{\phi\phi}}$ through the bifurcation sphere at $r_+$ into a second exterior, with the edge of the ergosphere marked as the class `ergo`; Kerr's throat is checked to have $\rho = 2GM/c^2$ at any spin and Kerr-Newman's $2GM/c^2 - r_Q^2/r_+$.
Schwarzschild-de Sitter's static slice at $t = 0$, drawn at $\Lambda = 0.2/r_s^2$, runs from the throat $r_h = 1.085\,r_s$ to the widest circle $r_c = 3.215\,r_s$, standing vertical at both, and on through the cosmological bifurcation sphere into the next static region, the same surface turned over, up to the next throat: one period of a chain without end, whose horizons are the exact roots of the cubic in the published $g_{rr}$, each coefficient of the expansion about a root reduced modulo the root's minimal polynomial before it is evaluated.
De Sitter's static slice at $t = 0$ is the sphere of radius $\ell = \sqrt{3/\Lambda}$, drawn at $\Lambda = 3$: the static chart covers the hemisphere out to the horizon and the antipodal observer's patch the other, the waist of the hyperboloid; the flat slicing's slices are flat, which the view states and the script checks.
Vaidya is its imploding shell of radiation at four moments, played as a movie with a frame every $r_s/16$ of $v - r$ and the rim $r = 4\,r_s$ at $z = 0$ in every frame, slices of constant $v - r$, since a slice of constant $v$ is null; `Slice` pulls the published metric back along $v = T + r$ to $(1 + 2Gm/c^2r)\,dr^2 + r^2d\phi^2$, a flat disc inside the shell, which is checked, and $z^2 = 4r_sr$ outside, Flamm's paraboloid moved in by $r_s$, meeting at a `crease`.
Oppenheimer-Snyder is its dust released from rest at $R_0 = 2\,r_s$ at four moments of the dust's proper time, played as a movie with a frame about every $0.07\,r_s$ of $c\tau$ and the clocks released at $4\,r_s$, the rim of the drawing, at $z = 0$ in every frame: inside, the published interior chart at $a = (a_m/2)(1 + \cos\eta)$, checked to make the published $G^\chi{}_\chi$ vanish, a cap of a sphere of radius $a$; outside, the moment carried on as Novikov's slice, clocks released from rest at every radius with the dust, read from Tolman-Bondi's comoving chart with no dust in it, since a slice of constant Schwarzschild $t$ meets the dust at an angle after the release and cannot reach it inside $r_s$. `RestCloud` solves each shell's cycloid and is checked to make Tolman-Bondi's published $G^r{}_r$ vanish and its $G^t{}_t$ give the declared density; the two sides meet with one tangent, and at the release the outside is Flamm's paraboloid.
Tolman-Bondi is the cloud its spacetime diagram declares, density falling as $1 - r^2/r_b^2$ with $2GM/c^2 = r_b/2$, but released from rest, $E = -GM(r)/c^2r$, at four moments before its centre is crushed, from its own comoving chart with `RestCloud`, checked against its published field equations; the spacetime diagram's cloud is marginally bound, $E = 0$, and its slices are planes, which the view states and the script checks.
Every surface of the cloud is checked at the height the metric gives it and then drawn twice as tall, `vertical` 2, since at its own height the cloud rises by only a third of its width and reads on the page as a flat disc, as the captain found on 30 September 2026.
Bertotti-Robinson, a product of two dimensional anti-de Sitter space and a sphere of the one radius $b$, has two views: its equator at one moment, a line of the first factor times a great circle of the second, the cylinder $\rho = b$, $z = b\ln r$, and its sphere of $\theta$ and $\phi$ at one event, the second factor; the first factor is Lorentzian and has no surface, which the view states.
Van Stockum's rotating dust, at $R = 1$ as its spacetime diagrams declare, is the plane $z = 0$ at one moment, whose circles grow to $r = R/\sqrt{2}$ and shrink after, so the surface curls back toward the axis until $g_{rr} = (d\rho/dr)^2$ at $r = 0.83\,R$, found by sympy; the view states that nothing is drawn beyond, and that beyond $r = R$ the circles are closed timelike curves, both checked.
Taub-NUT, at $m = 1$ and $l = m/2$ as its spacetime diagram declares, is the equator of a slice of constant $t$, where $g_{t\phi} \propto \cos\theta$ vanishes, out of its horizon $r_+ = m + \sqrt{m^2 + l^2}$, whose circumference radius $\sqrt{2(mr_+ + l^2)}$ is checked; across the horizon lies Taub's cosmology, where $r$ is a time, which the view states and the script checks.
Godel's universe, at $\omega = 1$ as its spacetime diagrams declare, is the plane $z = 0$ about one world line of the dust at one moment of the cylindrical chart, which stops exactly where $\sinh^2 r = 1/\sqrt{2}$, since $g_{rr} = (d\rho/dr)^2$ there reduces to $4s^3 = 2s$; the view states that beyond it no surface carries the slice and that beyond $\sinh r = 1$ the circles are closed timelike curves, both checked.
Ellis-Bronnikov is the same catenoid in its own chart, whose $r$ is the proper distance from the throat, so it is one piece through the throat, checked against $z = \ell\,\mathrm{arcsinh}(r/\ell)$ and $\rho = \sqrt{r^2 + \ell^2}$.
The cosmic string is its cone at $4G\mu/c^2 = 0.1$ with Gott's core rounding the apex, and the figure lays the cone flat beside it, cut along one meridian, so that the missing wedge $\delta = 36°$ shows.
Its first view is a movie of the ideal string's cone unrolling onto the plane, as the captain asked on 30 September 2026, cut along $\phi = 0$: in each frame the cone of half angle $\alpha$ with $\sin\alpha = 2\pi(1 - 4G\mu/c^2)/(2\pi - \Delta\phi)$ carries the whole of it round $2\pi - \Delta\phi$ of the axis, so the cut's edges stand $\Delta\phi$ apart, from $0$, the cone, to the deficit angle $8\pi G\mu/c^2 = 36°$, where it lies flat.
Every distance from the apex and every circle keeps its length in every frame, which the script checks, as it checks the wedge of the last.
Gott's core is a cap of a sphere and does not unroll without stretching, so it stays on the second view.
The Einstein static universe is the equator of one moment of its hyperspherical chart, the sphere $\rho = R\sin\chi$, $z = -R\cos\chi$ from the pole to the antipode, the same at every moment, with the equator where the areal chart and Einstein's projection end marked as `chartedge`; the areal chart's slice is checked to be its lower hemisphere.
Its conformal diagram is the strip $0 \le \chi \le \pi$ of $\eta = ct/R$, and three views draw Minkowski, de Sitter and anti-de Sitter space inside it by the maps their own diagrams use, each checked to make the pulled back metric of the Einstein static universe one conformal factor times the other's, on its plane of $t$ and $r$ and on its spheres; `_tools/derivations/einstein_static.md` is the derivation.
The C-metric, at $\alpha m = 1/6$ and $C = 1/(1 + 2\alpha m) = 3/4$ as its other diagrams declare, has two views: the equator of $t = 0$, $dr^2/Q + C^2r^2\,d\phi^2$, from one black hole's horizon at $r = 2m$ out to the acceleration horizon at $r = 1/\alpha$, its widest circle, and on through that horizon's bifurcation into the second black hole's exterior, $\rho = Cr$ checked; and the black hole horizon itself, smooth at $\theta = 0$ and the apex of a cone at $\theta = \pi$, where the string meets it, with $d\rho/ds$ checked to be $C(1 \pm 2\alpha m)$ at its poles and its area against Griffiths, Krtouš and Podolský's $16\pi Cm^2/(1 - 4\alpha^2m^2)$.
The Milne universe is the equator of a moment of its comoving hyperbolic chart at $ct = 0.5$, $1$, $2$ and $3$ in any length $\ell$, played as a movie with a frame every $0.05$ of $ct$: its circles grow faster than the distance out to them at every $\chi > 0$, which is checked, so each moment is drawn in Minkowski space as the sheet $Z = ct\cosh\chi$, $\rho = ct\sinh\chi$ out to $\rho = 4\,\ell$, inside the light cone $Z = \rho$ drawn as a reference.
That Minkowski space is the plane $\theta = \pi/2$ of its own inertial chart, $Z = cT$, so every sheet stands where its moment lies; `_tools/derivations/milne.md` is the derivation.
The Einstein-Rosen waves are the pulse of Weber, Wheeler, and Bonnor at $C = a$ going out from the axis, the plane $z = 0$ of the cylindrical chart out to $\rho = 10\,a$ at the moments $ct = a$, $2a$, $4a$ and $8a$, played as a movie with a frame every $0.1\,a$ of $ct$ and the axis at $z = 0$ in every frame, with $\psi$ and $\gamma$ entering the slice as numbers and the pulse checked to make every published Einstein component vanish; each circle's radius is checked against $\rho e^{-\psi}$, and at $t = 0$ no surface of revolution in flat space carries the moment out to $\rho = 1.72\,a$, which the view states and the script checks.
FRW is its closed universe of dust at five moments, played as a movie with a frame about every $0.1$ of $ct$, $a = 1 - \cos\eta$ as its conformal diagram declares, checked to make the published $G^r{}_r$ vanish; its flat slices are planes and its open slices have no surface of revolution in flat space, and no surface at all as a whole, which the view states under `stops` and the script checks at the scale factors of dust.

### The file, which the application reads

The application downloads the file for a spacetime whose index entry carries `embedding`, the stamp of the whole file, which changes when the file does and not otherwise.
Every number is a plain JSON number, every text is TeX in `$...$` inside prose, as the conformal diagram's texts are, and nothing in the file needs any relativity to draw.

    {
      "metric": "schwarzschild",
      "source": [{"metric": "schwarzschild", "system": "spherical", "fields": [...], "version": "..."}],
      "views": [view, ...]
    }

`source` is what `build_mfs_data.py` checks the file against, one entry for every coordinate system it read, which may be another spacetime's, as the star reads `schwarzschild.json`; the application can ignore it.
`views` holds at least one view, unless the spacetime has no surface to draw, which no spacetime lacks today.
Then `views` is empty and `stops`, a list of TeX sentences, says why, printed under the heading after "not drawn".
A file with a view never carries `stops` of its own; what a view does not draw is in the view's `stops`.
A view is:

- `id`: unique in the file.
- `label`: the name of its button, TeX, wanted only when a chart has more than one view to choose from.
- `system`, optional: the coordinate system the view belongs to. A view without it belongs to every chart. Show the views that name the chart being read and then the views that name none, as the conformal diagram does, and nothing if that leaves none.
- `unit`: TeX naming the length every number of the surface is measured in, as `$r_s$`, `$b_0$`, `$\ell$` or `$1/\sqrt{k}$`.
- `surfaces`: at least one surface. One surface is one moment; more than one is a sequence of moments in the order of their `time`.
- `figure`: the page's drawing of the view, described below, with what it takes to draw it again from another camera; an application that draws the surfaces itself can ignore it.
- `caption`: a list of paragraphs.
- `settings`, optional: TeX prose giving the parameter values drawn at, printed after "drawn at". It names the length every number is measured in, as "$r_s = 1$, the unit of every length".
- `input`, optional: TeX prose naming a declared function or matter, printed after "drawn with".
- `stops`, optional: a list of TeX sentences, each saying where the construction stops or what in this spacetime has no surface, printed after "not drawn". FRW's view carries its flat and open universes here, which have no surface in the file.
- `height`, present exactly when a surface of the view has a grid piece: TeX prose saying what the height stands for and at what scale, printed after "height", as "$\theta$, a height of $R$ for $\theta = 4c/R$".
- `vertical`, optional: a whole number greater than 1, present when every `z` of the view's surfaces is that many times the embedding's height, so the surface stands out where the metric makes it too shallow to read on the page; the first sentence of the caption says "(vertical scale $\times N$)", and a client draws the numbers as they stand. Only Tolman-Bondi's cloud carries it, at 2.
- `movie`, optional: the surface played as a movie, described below.
- `shades`, optional: the parts of a piece tinted while a view of the conformal diagram is shown, described under "The shade of a conformal region" below. Only the Einstein static universe's sphere carries them.

A surface is:

- `label` and `time`, in a sequence only: the moment as TeX, and its time as a number in `unit`, the time the view's `settings` name: $ct$ for FRW and $v - r$ for Vaidya.
- `axis`, optional: `{"from", "to"}`, the axis of time a stack of moments stands on, drawn up its axis from the height `from` to the height `to` in the class `axis`.
- `pieces`: its profile curves, at least one.
- `rings`: circles marked on it, possibly none.
- `curves`, optional: curves marked on it that are no circles about its axis, as a ring of free particles stretched into an ellipse on a flat plane, each `{"piece", "class", "closed", "points"}`: the points `[X, Y, Z]` in the surface's own frame, $Z$ along its axis, $X$ along $\phi = 0$ and $Y$ along $\phi = \pi/2$, rounded as the piece's $\rho$ and $z$ are, and `closed` present and `true` when the last point joins the first. Every point lies on the piece, $\sqrt{X^2 + Y^2}$ a distance its profile reaches and $Z$ the height there, or on a grid piece its triangles, a node of the row at its height on a stack of ellipses. A curve turns with the surface, as its meridians do. A curve that is one moment of a stack, as a ring of particles at its time, also carries that moment's `label` and `time`, as a surface of a sequence does, and each is a moment the other diagrams mark.
- `dots`, optional: points marked on it, each `{"piece", "class", "at"}`, `at` an `[X, Y, Z]` on the piece, as twelve of a ring's particles.

Every piece but a grid piece, defined below, is a profile curve in a half plane through the vertical axis $z$, and a surface of such pieces is a surface of revolution about it.
Its `points` are `[x, rho, z]`, running from its `start` to its `end` with `x` strictly increasing or strictly decreasing, at least two of them:

- `rho` is the distance from the axis, never negative, and `z` the height, both in `unit`, rounded to below $10^{-7}$ of the piece's own extent and to at least six decimals, so a small piece, as a collapsing star's shrinking cap, keeps the precision of a large one; the Einstein-Rosen pulse, which at $ct = a$ bends its surface by the axis on a scale of $a/4$ against a drawing $20\,a$ across, is rounded to below $10^{-8}$;
- `x` is the value of the piece's `coordinate` at the point, as that coordinate system writes it with the view's `settings`, a length in `unit`, an angle in radians such as $\chi$, or a number such as FRW's comoving $r$, written as the double it was computed at, the shortest decimal that reads back as it, since next to an irrational horizon a rounded $x$ would fall inside it; a client needs it only to label a point or to find one.

No number is written as `-0.0`.
A piece holds the surface for the values of `x` from its first point to its last and says nothing beyond them; its `start` and `end` say what the surface does at each.
The pieces drawn on 28 September 2026 have 9 to 316 points each, and the curves 241 to 1102.

A grid piece carries `grid` and `edge` in place of `points`, `coordinate`, `start` and `end`, and draws a quantity as a height over the plane $z = 0$ of the surface's frame, which the view's `height` names:

- `frame`: `"polar"`, where `u` is the distance from the axis and `v` the angle from $X$ toward $Y$ in radians, the grid running round from its last angle back to its first, so that a node stands at $(u\cos v, u\sin v)$; `"cartesian"`, where `u` is $X$ and `v` is $Y$; or `"ellipses"`, a stack of ellipses, defined below.
- `u` and `v`: the values of the grid's rows and columns, each strictly increasing, at least two of each, a polar grid's `u` from $0$ or more and its `v`, at least three of them, within $[0, 2\pi)$, each written as the double it was computed at.
- `z`: a row for each value of `u`, each with a height for each value of `v`, in `unit`, rounded to below $10^{-7}$ of the grid's extent and to at least six decimals; where a polar grid's `u` is $0$ the row holds one height.

The surface is flat between its nodes.
Each cell between neighbouring rows $i$ and $i + 1$ and neighbouring columns $j$ and $j + 1$ is the two triangles $(i, j)$, $(i + 1, j)$, $(i + 1, j + 1)$ and $(i, j)$, $(i + 1, j + 1)$, $(i, j + 1)$, a polar grid's last column joining its first, so a line of the grid runs straight from node to node, a polar grid's circles included, and a client needs nothing but the nodes.
The triangles lie within $5 \times 10^{-4}$ of the drawing's size of the quantity they stand for, measured in space; on 28 September 2026 that took Alcubierre's polar grid 43 values of $u$ and 72 of $v$, and the Krasnikov tube's Cartesian grid 155 and 91.
`edge` says, as an end does, what lies beyond the grid's rim, and its kind is always `edge`.
A grid piece carries no rings, since no circle on it keeps one height, and every curve on it lies on its triangles.

A stack of ellipses, `"frame": "ellipses"`, draws a surface whose every row is an ellipse about the axis, and draws no height over a plane: the world tube of a ring of particles, each row the ring at the time `u`, and each frame of the cosmic string's cone unrolling, each row the circle at the distance `u` from its apex.
In place of `z` as rows of heights it carries a value for each row:

- `a` and `b`: the row's semi-axes along $X$ and $Y$, never negative, rounded as heights are.
- `z`: the row's height.
- `v`: the angles of its columns, within one turn, the grid running round from its last angle to its first unless it carries `"open": true`, as a cone cut along a meridian and opened does, whose `v` may run up to a whole turn.
- `lines`: the lines drawn on it, each `{"class", "u", "v"}`, the rows at the indices `u` and the columns at the indices `v`, as a ring particle's world line in the class `worldline` or the cut's edges in `cut`, since each frame of a movie holds its own rows.

The node of row $i$ and column $j$ stands at $(a_i\cos v_j, b_i\sin v_j, z_i)$, and each cell is two flat triangles as a polar grid's are, rows up and columns round, which turns every triangle one way, so the normal at a node is the sum of the triangles' normals, each $-(B - A) \times (C - A)$, pointing away from the axis, and a node that no triangle with any area meets, as a cone's apex at the end of its cut, takes the normal of the node before it in its row, or after it at the row's start.
It draws only the lines it names and its outline, and no rim: a tube's ends are its first and last moments, which it marks as curves.
A stack carries no `height` of its own meaning a quantity: the view's `height` says what its height is, as a time.

To draw a piece, turn every point about the axis, $(\rho\cos\phi, \rho\sin\phi, z)$ for $\phi$ from $0$ to $2\pi$, and join neighbouring points and neighbouring angles.
The points are close enough that straight segments between them are the surface to within $2 \times 10^{-5}$ of the drawing's size, the diameter of the widest circle drawn, so no smoothing is wanted, and a client may add angles as finely as it likes.
The pieces of a surface may be drawn in any order, since every piece but a reference piece is part of one surface and no two of them overlap.
Draw $\rho$ and $z$ at one scale, since a surface stretched along either axis no longer has the metric's distances.
The zero of $z$ is wherever the construction put it and means nothing; only differences in $z$ do, and every piece of one surface shares one $z$.
In a sequence every surface stands on its own axis at the origin, and placing them side by side or showing them one after another is the client's choice.

A view's `movie` plays a run of frames in turn, as the captain asked on 30 September 2026 for the embedding diagrams that change in time and for the cosmic string's cone unrolling:

- `frames`: at least two surfaces, each as `surfaces` defines one, with its `label` and in place of `time` its `value`, strictly increasing from frame to frame; a movie in time holds each surface of `surfaces` as one of its frames, with the same pieces and rings, and runs from the first to the last.
- `variable`: TeX naming what `value` is, as `$ct$` or `$\Delta\phi$`.
- `seconds`: how long one pass takes: `value` runs at a steady rate from the first frame to the last in that time, the last frame is held for one step of the mean, and the movie starts again.
- `turns`: `false` where turning the figure adds nothing, as for FRW's sphere, which looks the same from every side; absent otherwise.

Every frame stands on its own axis at the one origin, and `figure` draws the first frame alone in a box that holds every frame, with a label that carries `frame` and names the frame shown.
A client that does not play movies draws `surfaces` as a sequence.

A piece also carries:

- `id`: unique in its surface.
- `class`: `sheet` for the slice a chart of the spacetime covers, `sheet2` for the same surface run on past where that chart ends, as the other exterior of Schwarzschild, the other side of a wormhole or the far hemisphere of a closed universe, `star` for the slice through matter, as a star or Gott's core, and `reference` for a surface that is not part of the slice at all, drawn for comparison.
- `reference`: `true` on a reference piece, and absent otherwise. A reference piece is drawn faintly and hides nothing, as the vacuum paraboloid under the star and the cone of an ideal string under Gott's core are.
- `space`, optional: `"minkowski"` for a piece drawn in three dimensional Minkowski space, as anti-de Sitter's and the BTZ hole's beyond $r = \sqrt{2}\,\ell$, where every length along the piece is measured with $dX^2 + dY^2 - dZ^2$. A client draws such a piece as any other and tells the reader how its lengths are measured, as the view's `settings` and caption do. A piece without it is drawn in flat space.
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
Its `class` is `r` for a circle of constant coordinate on a `sheet` or `star` piece, `r2` for one on a `sheet2` piece, `horizon` for the horizon, `throat` for a wormhole's throat, `surface` for the edge of matter, a star's surface or the edge of Gott's core, `chartedge` for the circle where a chart ends or a drawing stops, as FRW's equator or where Reissner-Nordstrom's inner surface lies level, `space` for the circle where a surface passes from flat space into Minkowski space, as the BTZ hole's moment does where it lies level, `ergo` for the edge of an ergosphere on Kerr's and Kerr-Newman's equator, `wall` for the wall of a warp bubble, and `reference` for a circle on a reference piece.
A curve's `class` is `particles` for a ring of free particles, whose `dots` carry it too, `path` for the path of a warp bubble or the ship that builds the Krasnikov tube, `wall` for the circle $v_sf = 1$ on Alcubierre's height, the circle $r_s = R$ on Natário's and the level line $k = 0$ on the Krasnikov tube's, `contract` and `expand` for the level lines of Alcubierre's height where the expansion of the observers riding its slices is half its greatest value, ahead of the ship and behind it, and `flow` for a line of Natário's flow.

The `figure` is the drawing the page makes of the view, projected by the script once from a fixed camera, in the form a figure in three dimensions takes in the diagram files, and drawn again from the surfaces at whatever camera a reader turns it to, as "Turning the figure" below says:

- `box`: `[Xmin, Xmax, Ymin, Ymax]` of the plane of the page, with $Y$ up, drawn at one scale.
- `camera`: `{"azimuth", "elevation"}`, the angles in degrees the surfaces were projected from, the azimuth measured round the axis from $\phi = 0$ toward $\phi = \pi/2$ and the elevation above the plane $z = 0$, looking at the origin.
- `layers`: painted in the order given, each with a `kind` and a `class`. A `fill` has `points`, a closed polygon, and optional `holes`, polygons painted with the even odd rule; a `line` has `points`, a polyline; a `point` has `at`. Every point is `[X, Y]` in the plane of the page, rounded to four decimals. A layer drawn in the plane of the page rather than on a surface, as the cone laid flat is, carries `"flat": true`.
- `labels`: TeX at a point, `{"at", "text", "anchor", "class", "dx", "dy"}`, the anchor and the offset in units of a figure 628 wide, as a conformal diagram's are, `class` being `lab` or `small`, and for a label that names a circle `ring` and sometimes `clear`, for one that names a curve `curve`, for one that names an axis of time `axis`, and for a movie's name of the frame shown `frame`, which "Turning the figure" defines.
  A `small` label names a circle and stands beside the end of it, or names a moment of a sequence and stands below its surface; either way it can stand over the surface and its lines, so it is set on a ground of its own, dark on the screen and white in print, as the page sets it; a `lab` label stands where no line runs.
  No two labels overlap, which the script checks before it writes.
- `legend`: `[kind, class, text]` for each class it names, the kind being `fill`, `line` or `point`.
- `turn`: what a client needs to draw the figure again from another camera, which "Turning the figure" defines.

A `sheet` piece is tinted `cover` and a `star` piece `star`, while a `sheet2` or reference piece is left clear, so the part a chart covers stands out; each is tinted only where it is the surface nearest the camera, found by casting rays through the same truncated cones a client draws, and every line on the surfaces is split where another part of a surface hides it, the hidden part carrying its class with `-far` appended, which the page draws faint.
A reference piece hides nothing.
A point of a line is hidden when a ray from it toward the camera meets one of those cones further than $10^{-7}$ of the drawing's size along it.
A point of the outline is where the line of sight grazes the surface, and there the cone between two circles of a concave profile dips across that line by less than $10^{-5}$ of the drawing, so each point of the outline is judged by the two points $10^{-5}$ of the drawing off the surface on either side along its normal and is hidden only if both are.
The fills are `cover`, `star` and `wedge`.
The lines are `outline`, where the surface turns edge on to the camera and the rim at an `edge` or `stops` end or of a grid; `grid`, the rows and columns of a grid piece that `turn` names; `meridian`, the profile turned to evenly spaced angles, as the legend says; `reference`, the meridians and outline of a reference piece in dashes; the ring classes above, each ring drawn as its circle, or as a `point` where its $\rho$ is zero; the curve classes above, each curve drawn through its points, and `particles` also as a `point` at each of its `dots`; `worldline`, a particle's world line up a stack; `axis`, the axis of time a stack stands on; and `cut`.
Each line class may also come with `-far`.
The page styles each class as `.em-<class>` for a line and `.em-<class>-fill` for a fill, on the screen and in print, and a test holds every class a figure paints to having both.
The cosmic string's figure also lays the cone flat beside it, a flat drawing in the plane of the page with the wedge the cone lacks painted as `wedge` and the two edges of the cut as `cut`.

### Turning the figure

The page shows the published figure until a reader drags it, and then draws it again from the surfaces at the camera the drag has reached, with `MFS/assets/turn.js`, which is geometry alone and which the tests run against every published figure through `_tools/turn_check.cjs`.
The application turns its figures by the same rules, so the two agree, and at the figure's own camera the rules give back the published figure: the same lines split at the same points, the same labels at the same places, and the tint to within a point of its grid.

`turn` holds:

- `origins`: for each surface, in the order of `surfaces`, the point `[X, Y]` of the page where its axis meets $z = 0$ at the figure's camera.
- `meridians`: how many meridians each piece carries, at $\phi = 2\pi k/n$, a reference piece every other one.
- `tint`: the fill each piece class is tinted with, as `{"sheet": "cover", "star": "star"}`; a piece of any other class is left clear and still hides what lies behind it.
- `marks`: the lines a surface carries besides its circles, each `{"class", "surface", "piece", "phi"}`, the meridian at `phi` of that piece of that surface drawn in that class, as the cosmic string's `cut`.
- `grid`, present when a surface has a grid piece: the lines of each grid piece the figure draws, each `{"class", "surface", "piece", "u", "v"}`, the rows at the indices `u` and the columns at the indices `v` of that piece's grid, drawn in that class, `grid`.

A label that names a circle carries `ring`, `{"surface", "ring", "side"}`: the surface, the circle's place in that surface's `rings`, and `1` where the label stands to the right of the circle's end on the page or `-1` to the left.
It may also carry `clear`, a length in `unit`: the label stands past the outline wherever the outline crosses a height within `clear` of the circle's end outside that end, as the cone's label does below the cone's rim.
A label that names a curve carries `curve`, `{"surface", "curve", "side"}`, and stands beside the point of that curve farthest right on the page, `side` $1$, or farthest left, $-1$, the first of them where two tie, as a stack's moments are named beside their rings.
A label that names an axis of time carries `axis`, `{"surface"}`, and stands at the top of that surface's `axis`.
A movie's label that carries `frame` shows the `label` of the frame drawn, and stays where it is.
A label without `ring`, `curve` or `axis` stays where it is, as the moments of a sequence and the words on the cone laid flat do, and so does a layer with `flat`.

At another camera, azimuth $a$ and elevation $e$, with the vectors `right`, `up` and `toward` of `projections.Camera`:

- Each surface turns about its centre, the point of its axis at $z_c$ halfway between the lowest and the highest point of its pieces, which stays where the figure's camera put it, $(0, z_c\cos e_0)$ from its origin.
- A surface of revolution is as wide on the page as its widest circle from every camera, but its height, from $\min\,((z - z_c)\cos e - \rho|\sin e|)$ to $\max\,((z - z_c)\cos e + \rho|\sin e|)$ over the points of its pieces, changes with $e$.
  So all the surfaces are drawn at one scale $s$, never above 1, the largest that keeps each within the height it had at the figure's camera, and each is moved up or down only as far as that height needs; the `box` stays as it is, and every drawn point stays inside it.
- A grid piece is flat between its nodes, so its nodes bound it: its height on the page runs over $X u_0 + Y u_1 + (z - z_c)\cos e$ with `up` $= (u_0, u_1, \cos e)$, and since a grid need not be as wide from every side, a surface with one is also kept within the width it had at the figure's camera on each side of its axis, the largest of $X r_0 + Y r_1$ and of its negative over its nodes with `right` $= (r_0, r_1, 0)$, which draws the Krasnikov tube's rectangle smaller as it turns toward its diagonal.
- Each piece is drawn as the script draws it: its meridians, its outline, the circle at an `edge` or `stops` end, its rings, and the marks, every line split into the parts seen and hidden as above, cut halfway between points of the two kinds.
  Each curve of a surface is drawn through its points, and back to the first where it is `closed`, split the same way, and each of its `dots` as a `point` where it stands, hidden or not, as the script draws them.
- The outline of a piece runs through the angles on the circle of each of its points at which $\cos(\phi - a) = (d\rho/dz)\tan e$, with the tangent taken as numpy's `gradient` takes it, joined from point to point; a run of such points inside the piece turns at both of its ends, where its two sides meet, and a run that reaches an end of the piece runs on to that end on each side.
- A grid piece is drawn as the script draws it: the rows and columns `turn.grid` names, each from node to node and round to its first node where a polar grid runs round, its rim, and its outline, where the normal at its nodes is square to the line of sight.
  The normal at a node is the sum of the normals of the triangles that meet there, each pointing up and as long as twice the triangle's area, and the outline is found by marching squares over the grid as `isolines` finds it, a polar grid's first column repeated after its last, so each of its points lies on an edge of the grid.
  What hides a point is found from the triangles that cover its point of the page seen along `toward`, which a ray toward the camera meets after the difference of their depths, the triangles sorted into bins over the page by their boxes; every point of a line on a grid piece is judged by the two points $10^{-5}$ of the drawing above and below it and is hidden only if both are, since above and below are the two sides of a height.
- Each tinted piece is filled where it is the surface nearest the camera, found on a grid of 420 points across the longer side of the `box`, a grid piece by its triangles.
  While a shade is shown, the cones of its piece between the circles at its `from` and `to` are tinted `shade` and the rest of that piece is left clear, as `shade()` in `MFS/assets/turn.js` says.
- A label with `ring` stands beside the end of its circle on its side, $(\pm s\rho, s(z - z_c)\cos e)$ from the surface's centre, with its offset as published; one with `curve` beside its curve's point farthest to its side on the page, and one with `axis` at the top of the axis, each with its offset as published.
- A stack of ellipses is drawn as the script draws it: the lines it names, each from node to node and round to its first node where it runs round, its outline, and its axis of time, split at every point where a ray toward the camera meets its triangles; a point of a line on it is judged a hundred thousandth of its width straight out from the axis and back, the two sides of a tube or a cone, or above and below on the axis itself.
- A movie draws one frame at a time, the first until the client chooses another, and every frame turns about one centre, halfway between the lowest and the highest point any frame reaches, at one scale and one move up or down that keep all the frames together within the height and width they had at the figure's camera, so the frames keep one size and place as they change and as the reader turns them.
  It names its circle whether the circle runs in front of the surface or behind it, as the published labels do, and gives way only while it would overlap a label that stays put or a label before it in `labels`.
- The layers keep the published order: the flat fills, the tint, every line hidden and then every line seen, each in the order of the legend, and the points.

The page turns an embedding diagram and a figure of light cones by the same hand.
It turns a figure round its axis as far as the reader drags it, the drawing's width being half a turn, and tilts it from looking straight down the axis, $e = 90°$, to looking straight up it, $e = -90°$, never past, so the axis always stands up the page.
It stays where it is let go, and the reset button in the corner of the frame, a double click or a double tap, Home and Escape bring back the published figure, which is also what prints.
The arrow keys turn it by 15 degrees once the drawing has the focus, which a click on it gives.
`dragged()`, `keyed()` and `turned()` in `MFS/assets/turn.js` are those rules, and the tests drag a surface of revolution, a height over a plane and a figure of light cones with them.

The page plays a movie as `wireMovie()` in `_layouts/mfs.html` says: it plays once drawn and only while it is on the screen, keeping each frame once drawn at the figure's own camera, and a reader who asks for reduced motion finds it paused on its first frame.
The button in the top left corner of its frame, in the reset's size and colours and pink while pressed, plays and pauses it, and the label below the drawing names the frame shown.
FRW's sphere does not turn under the hand, and every other movie turns while it plays.
The print copy prints the first frame.

`_tools/turn_drag.mjs` turns one figure of each kind in headless Chrome as a reader does, Schwarzschild's surface, Krasnikov's height and Gödel's light cones, with the mouse, the keys, the reset button and, on a phone, a finger.
It holds every point drawn to the box, every label shown to the drawing, the reset button to the top right corner at its one size, the figure to having no glow, a drag up to stopping at straight below the axis, a finger drawn up to scrolling the page, and the print copy, turned or not, to the published figure:

    bundle exec jekyll serve
    node _tools/turn_drag.mjs http://127.0.0.1:4000 --phone --text 48

It starts Chrome with `_tools/chrome.mjs`, as `page_timing.mjs` does, and exports its checks, so a driver for another browser runs the same ones.

## The three diagrams together

A spacetime's spacetime diagram, conformal diagram and embedding diagram stand together, in that order, after the geodesics and before the references, as the captain asked on 28 September 2026.
That is where the application's Graphs section stands, between its Maths and its References, and the page has no row of section buttons, so the three need no heading of their own.
The spacetime diagram and the conformal diagram each draw the slice the embedding diagram is cut from, so the reader meets the line and then the surface cut along it.

### The slice of the embedding diagram

Every embedding diagram draws part of one moment of its spacetime, or of several moments in turn, and `_tools/derivations/embedding_slices.md` works out how each moment meets every other diagram of its spacetime.
`_tools/derivations/slices.py` is that table in code.
It reads each moment and how far the embedding reaches along it from the embedding file itself, and says what the moment is in the chart of each spacetime diagram's view: a line of constant time, the Eddington-Finkelstein curves of Schwarzschild's $t = 0$, de Sitter's static $t = 0$ in the flat slicing, the pp-wave's fronts $u = u_k$, Novikov's slice outside Oppenheimer and Snyder's dust, and the whole equator seen from above on Kerr's and Kerr-Newman's views of the principal rays.
`/tmp/mfs-venv/bin/python _tools/derivations/slices.py` checks every chart transformation a moment is carried through by pulling one published metric back onto the other, and Novikov's shells against the published exterior, sixteen checks in about a second.
`HIDDEN` in it names the drawings on which a moment of the spacetime lies and is not drawn, each with the reason: FRW's flat and open universes and Tolman-Bondi's marginally bound cloud are other spacetimes than the ones embedded, Vaidya's outgoing chart draws the time reverse of the shell embedded, and beyond $r_c$ Gödel's and van Stockum's circles are closed timelike curves.

Each drawing draws its slices with the map it draws everything else with: `null_rays.py` writes them into each flat view, `projections.py` into each figure and `conformal.py` into each conformal view.
`null_rays.py --slices` rewrites only the slices of the diagram files, from the embedding files as they stand, tracing no ray, and redraws the figures, which take seconds; run it, and `conformal.py`, after an embedding diagram moves a moment.
Each slice records a stamp over the embedding surface it marks, computed by `embedding_moment_version` in `build_mfs_data.py` over the view's settings and the surface's label, time and pieces, and the command refuses a slice whose stamp no longer matches, or that names an embedding view or surface the embedding file does not draw, in `--check` and when writing alike.
The `Slices` tests hold every point of every slice to its moment and every line's ends to the embedding's reach or the drawing's box, from the numbers in the files alone, and every moment to appearing on every drawing it lies on and on no other.

### The file, which the application reads

A flat view of a spacetime diagram, a figure in three dimensions and a view of a conformal diagram may each carry `slices`, a list with one entry for each surface of each embedding view whose moment lies on the drawing, in the order of the embedding file's views and surfaces:

- `view` and `surface`: the embedding view's `id` and the surface's place in its `surfaces`.
- `label`: the moment as TeX, the surface's own `label` in a sequence and otherwise the moment in the drawing's own words, as `$t = 0$`, `static $t = 0$` or `$t = 0$, $r = b$`.
- `version`: the stamp above; the application can ignore it.
- `lines`: polylines, `points`: points, and `fills`: regions, each a list of rings filled by the even odd rule, any of them possibly empty, all in the drawing's own coordinates: a flat view's unit square, as its rays, the right half of a view through the centre, which the reader mirrors as it mirrors the rays; a figure's plane of the page, as its layers; and a conformal view's $(X, T)$, as its layers.
- `place`, present where the drawing names the moment beside it: `at`, a point in the same coordinates, `anchor`, the corner of the label pinned there as a conformal label's anchor, and `dx` and `dy`, how far that corner stands from `at`, in the units the page draws the drawing in, a spacetime diagram's plot 520 wide and a conformal diagram 628 wide with its margin, $y$ down: to the edge of the line's stroke, 1.3 from the line, or of the point's dot, 5.1 from its centre.
  On a view through the centre `at` may have a negative first coordinate, the mirror of the point at its absolute value.

Draw a slice's regions under everything else on the drawing, and its lines and points over every other line and under the points that mark ends and events, as a solid line 2.6 wide of `--green-light`, #4db84a, on the dark ground and `--green-dark`, #2a6e2a, in print, a region tinted with the same green at 12%.
Green marks nothing else on these drawings but in dots, Kerr's ergosurface, so a solid green line reads at once as the new mark.
Set a label with `place` at its point, on a ground of its own as a conformal diagram's small labels are, in the same green, at the size of a small label on a conformal diagram and of a reference line's label on a spacetime diagram.
A slice without `place` is a region alone, the floor of a figure or the equator seen from above, and is named in the legend instead.

`slices.place` chooses each label's side for the drawing as it is drawn at the usual text size: at the right hand end of the slice's lines, on the line's upper edge over the rest of it, then on its lower edge, then past the end; then at the left hand end; then at a point 85% of the way along; for a point, round from above and to its right against the edge of its dot.
A label stands exactly on the edge of what it names, and labels may touch; nothing is added between them.
The first side that stays in the drawing and overlaps no other label wins, and among those, the first that covers none of its own line and none of another slice's.
A conformal diagram and a figure are drawn larger with the reader's text, their labels with them, so the place holds at every size.
A spacetime diagram's words grow around a plot that does not, so at a larger text size the page chooses, for all the labels of a plot together, the sides of their points and of the far ends of their lines that leave none past the plot and fewest over another label, in the same order, the file's own first, and a label wider than the room beside its point slides along to stay in the plot.
Where no choice keeps them clear, as at three times the usual text on a phone, the plot's labels are held smaller, a tenth at a time, until one does, as `fitWords` holds a heading to the size at which its widest word fits, and never below their size at the usual text; on 29 September 2026 that held Bianchi's at 0.9 and Oppenheimer and Snyder's dust at 0.73 of the reader's size, and every other plot kept it.
At the usual text size the choice is the file's place, and it returns there once the text is smaller again; `fitSliceLabels` in `_layouts/mfs.html` carries the details.

Each drawing's legend adds one entry for the slices of each embedding view it marks, shown with them: "the slice of the embedding diagram", "the slices of the embedding diagram" where the view is a sequence, and after it the moment where the drawing does not name it, as "the slice of the embedding diagram, $t = 0$".
Where two embedding views' entries would read the same, as Reissner-Nordström's outside $r_+$ and inside $r_-$, they are one entry.
Only the slices of the embedding view shown are drawn, so choosing another view of the embedding diagram draws its slices in place of the others, as a slider's moment would; the print copy prints every view of the embedding diagram, and so every slice.
A sequence's moments stand side by side in the embedding diagram, so all of them are drawn at once, and so are a movie's, whose frames pass through every one of them, and a stack's, which marks each moment as a ring at its time.
A moment of a stack is a slice that carries `curve` beside `surface`, the ring the stack marks at that time, and its stamp is taken over that ring's label and time as well.
The ideal string's cone unrolled is the same moment as the cone about Gott's core, marked on each drawing of the conical chart; Gott's own conformal diagram, a spacetime whose core replaces the apex of the ideal string's cone, marks only its own.

### The shade of a conformal region

The captain asked on 30 September 2026 that choosing Minkowski space, de Sitter space or anti-de Sitter space on the Einstein static universe's conformal diagram shade the part of the embedding diagram's sphere that region covers, in pink and without glow.
Each region is Hawking and Ellis's, found from the map each conformal view draws with and never written by hand: `conformal.covered()` takes the map from the other spacetime's own chart into the strip, finds by bisection in time the point of each radius at the conformal time $\eta$, and keeps the values of $\chi$ reached.
Minkowski space enters by its spherical chart through `mink_pq`, anti-de Sitter space by its static global chart through `ads_global_pq`, and de Sitter space by its closed slicing through `ds_closed_pq`, its hyperboloid embedding read backwards, since the static patch covers only the diamond $|\eta| + \chi < \pi/2$ about the pole.
`embedding.py` checks each region against Hawking and Ellis's, $\chi < \pi - |\eta|$, the whole sphere for $|\eta| < \pi/2$ and $\chi < \pi/2$, at $\eta$ from $-3$ to $3$ in steps of $1/2$, and `conformal.py` checks each region the strip draws against the same map, so the two diagrams agree.
The conformal diagram has no slider, so the shade is the region at the moment both diagrams draw, $t = 0$, which is $\eta = 0$, and its caption paragraph says so.

An embedding view's `shades` is a list with one entry for each view of the conformal diagram whose region the view shades, in the order of the conformal file's views:

- `view`: the `id` of the view in the spacetime's conformal file that the shade answers; the conformal file carries nothing new, and `build_mfs_data.py` refuses a shade whose view, surface, piece or circles do not exist.
- `surface` and `piece`: the surface's place in `surfaces` and the `id` of the piece shaded, a profile and never a grid piece.
- `from` and `to`: the values of the piece's `x` between which the moment is shaded, each written as the double of a point of its profile, so the shade ends on circles of the surface; at $\eta = 0$ Minkowski space shades $0$ to $\pi$, the sphere but its antipode $i^0$, which is one point, de Sitter space $0$ to $\pi$ and anti-de Sitter space $0$ to $\pi/2$.
- `T`: the conformal time of the moment shaded, the conformal view's own $T$, $0$ for the Einstein static universe.
- `reach`: the region at other conformal times, each `[T, from, to]` rounded to seven decimals, or `[T, null, null]` where the region covers none of the moment, as de Sitter space beyond $|\eta| = \pi/2$; a client with a slider of time can shade from it, and the tests hold it to Hawking and Ellis's.
- `layers`: the figure's fills with the shade in place of its piece's tint, at the figure's camera, in the form of the figure's own fill layers; the page paints them in place of the figure's fills and under all of its lines, and a client that draws the surfaces itself can ignore them.
- `legend`: the shade's legend entry, `["fill", "shade", text]`.
- `caption`: one TeX paragraph, read after the view's caption while the shade is shown.

While the conformal diagram shows a view that a shade answers, draw the shade and no other: the cones of its piece between `from` and `to` in `--pink-light`, #ff6da2, at 16% on the dark ground, the rest of that piece clear, every other piece as `tint` says, its legend entry first and every fill entry whose class the shade's `layers` do not paint hidden, and its caption paragraph after the view's caption.
Nothing about the shade glows, and it is shown at once with the view, with no transition.
A view the conformal diagram shows with no shade answering it leaves the figure as published.
A turned figure keeps its camera and takes the shade, as `shade()` in `MFS/assets/turn.js` draws it, and `turn_check.cjs` holds every shade drawn at the figure's camera to its published `layers`.
The print copy prints the page as it was first drawn, with the chart's own conformal view shown and so no shade.

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

The whole collection took 160 seconds on 29 September 2026, and `--system <metric_id>/<system_id>` checks one system in seconds, which is what to use while editing a single entry.
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

`tov`, `malament_hogarth`, `mixmaster`, `lentz`, `einstein_static`, `schwarzschild_de_sitter` and `milne`, the cylindrical chart of `godel`, and both charts of `c_metric` have their mathematics written by `_tools/derivations/print_charts.py`, which defines each of those charts, computes every tensor with the checker's own `Geometry`, and prints each value through `chart_printer.py` beside it:

    /tmp/mfs-venv/bin/python _tools/derivations/print_charts.py --metric tov

Every printed value is read back through the checker's `Reader` and compared with the value it came from before anything is written, and `verify_metrics.py` then checks the file like any other.
The script writes only the coordinate system's mathematics, the line element and the domains; the prose, the convention and the description of each parameter stay in the metric file and are carried over, and a parameter with no description stops the script.
A chart can pass a function that regroups the numerator of each value, which is how TOV's curvature is printed around $(\partial_r\Phi)^2 + \partial_r^2\Phi$ and Bianchi IX's around $a_1^2\cos^2\psi + a_2^2\sin^2\psi$, and a Ricci or Kretschmann scalar can be given in a structured form, which is checked against sympy before it is used.
A chart can instead pass a whole `pretty`, as Gödel's cylindrical chart passes `chart_printer.hyperbolic`, which rewrites the exponentials the checker's `Geometry` hands back in $\sinh r$ and $\cosh r$.
A chart the script writes replaces the chart of its id, or joins the spacetime's other charts after them, so Gödel's Cartesian chart stays as it was written.
A parameter keeps the description it has in that chart, or in another chart of the same spacetime.
`tov.md`, `malament_hogarth.md`, `mixmaster.md`, `lentz.md`, `godel.md` and `schwarzschild_de_sitter.md` in the same folder record why each chart is the one published; `godel.md` also records the transformation from the Cartesian chart, whose pullback sympy checks symbolically.
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
