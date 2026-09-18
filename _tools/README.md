# Site tools

Jekyll ignores any directory whose name starts with an underscore, so nothing in here is published.

## Adding a spacetime to My Favorite Spacetimes

Write the new metric file at `MFS/assets/data/metrics/<id>.json`, with its `id` equal to the file name.
It needs `name` (the long title on the page), `short_name` (the label in the search list) and `tags`.
Add `sort_name` only if the spacetime belongs somewhere else in the alphabetical list than its `short_name` puts it.
Cite references by their key in `assets/data/references.bib`.

Then run one command from the top of the repository:

    python3 _tools/build_mfs_data.py

That rewrites `MFS/assets/data/metrics_index.json` and `MFS/assets/data/references.json` from what is on disk.
The index is never edited by hand, so it cannot disagree with the folder.
Commit the new metric file together with whatever the command rewrote.

Every key a metric cites has to be an entry in `assets/data/references.bib`.
The command refuses to write anything if one is not, naming the metric file and the key it could not find, so a mistyped citation is caught here rather than published as a reference the reader cannot resolve.

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

## What the command writes

`MFS/assets/data/metrics_index.json` is one entry per metric file, in the order the search list shows them, carrying `id`, `name`, `tags` and `version`.
Its `name` is the metric's `short_name`, which is what the search list shows, not the long `name` the page titles the spacetime with.

`MFS/assets/data/references.json` is the whole of `assets/data/references.bib` parsed into JSON, so a reader can pull every reference the collection cites in one fetch instead of walking every metric file.
The `.bib` file stays where it is and stays the one place a reference is written; the publications page still reads it.

## The version stamps

A metric's `version` is a hash of that metric's content.
It changes when the content changes and not otherwise, so reformatting a file or republishing the site leaves it alone.
This is what the iOS application uses to fetch only the spacetimes that actually changed, since GitHub Pages rewrites every file's own tag on every publish.
`references.json` carries one stamp of its own at the top for the same reason.

## Checking

    python3 _tools/build_mfs_data.py --check

reports whether the published files are still what the folder says they should be, and changes nothing.

    python3 -m unittest discover -s _tools

runs the tests, which include that check.

## Checking the physics

`build_mfs_data.py` checks the shape of the data. It does not read the mathematics.
`_tools/derivations/verify_metrics.py` does, for every coordinate system of every metric file.

It builds $g_{\mu\nu}$ from the line element each entry prints, computes the inverse metric, both Christoffel variants, Riemann, Ricci, the Ricci scalar, Kretschmann, Einstein and Weyl in sympy, and compares all of it against every value the entry publishes.
A component the entry omits has to vanish, so a missing symbol is caught as well as a wrong one.
It names each disagreement by file, system and symbol, and exits non-zero if any remain.

sympy is not installed system wide, and the virtual environment does not belong in the repository:

    python3 -m venv /tmp/mfs-venv && /tmp/mfs-venv/bin/pip install sympy
    /tmp/mfs-venv/bin/python _tools/derivations/verify_metrics.py

The whole collection takes about ten minutes. `--system <metric_id>/<system_id>` checks one system and takes seconds, which is what to use while editing a single entry.
`--budget <seconds>` changes how long sympy may spend on one tensor.

Nothing is ever passed in silence. A value that cannot be parsed, a system with no time coordinate declaration, or a tensor sympy cannot finish in the budget is reported as `UNCHECKED` with the reason, separately from the disagreements.

## The three conventions the checker encodes

All three are things the files do consistently rather than things they write down, so they are recorded here and in the script's header.

The chart is $x^0 = cT$.
A coordinate carrying dimensions of time is not itself the chart coordinate; $c$ times it is, and every published component is a component in that chart even though the index is printed with the bare name.
Schwarzschild shows it plainly, publishing $g_{tt} = -(1-r_s/r)$ against a line element whose time term is $-(1-r_s/r)c^2dt^2$.
Because the rescaling is linear, a component in the chart is the one taken with the bare coordinate multiplied by $c$ once per upper time index and divided by $c$ once per lower one.

The Ricci tensor is contracted as $R_{\mu\nu} = R^\alpha{}_{\mu\nu\alpha}$, on the last lower index rather than the first.
That is the opposite sign from the commoner $R^\alpha{}_{\mu\alpha\nu}$, and it carries through to the Einstein tensor and the Ricci scalar.
Contracting each published Riemann tensor both ways reproduces the published Ricci in every slot this way and not the other, in Ellis-Bronnikov, Reissner-Nordstrom and Godel alike.
The Weyl tensor is the exception and is built from $R^\alpha{}_{\mu\alpha\nu}$, because it is defined by removing the traces of Riemann and those traces do not care what the file calls Ricci.

The dots in a geodesic equation are velocities of that same chart, so a dot on a time coordinate means $d(cT)/d\lambda$ and not $dT/d\lambda$, even though it is printed on the bare letter.
It has to be that one, because the equation is $\ddot x^\mu + \Gamma^\mu{}_{\nu\rho}\dot x^\nu\dot x^\rho = 0$ with the same printed $\Gamma$ the entry lists, and those are chart symbols.
The reading is checkable without sympy: on it every term of every equation carries the dimensions of its left hand side, and on the other reading a term mixing a time velocity with a space velocity comes out wrong by one factor of $c$.
A geodesic equation is the only thing an entry publishes that adds a time derivative to a space derivative, so it is the only place the mistake can hide, and `_tools/derivations/kasner.md` Step 17 works the check through term by term.

## Adding a spacetime, as far as the checker is concerned

Telling a time coordinate from a length needs the dimensions of the parameters, which live in prose, so the script cannot read it off the file.
`TIME_COORDINATES` in `verify_metrics.py` declares it per system instead.
A new coordinate system that is not listed there is reported `UNCHECKED` rather than guessed at, so adding a spacetime means adding its line, and forgetting to is visible rather than silent.

An entry whose parameters are not free needs a second line, in `PARAMETER_RELATIONS`.
Kasner prints three exponents bound by $\sum_i p_i = \sum_i p_i^2 = 1$ and claims its values only on the surface those two equations cut out; its Ricci tensor is zero there and nowhere else.
Checked against free symbols it would report seventy four disagreements that are not disagreements.
The table carries a rational parametrisation of the surface instead, and every constrained parameter is replaced by its parametrised value on both sides of each comparison, so the identity checked is the one the entry asserts.
That is exact rather than a sample, because a rational parametrisation covers a dense subset of an irreducible variety.
`_tools/derivations/kasner.md` derives the parametrisation the Kasner entry uses and why the constraints are the field equations.
