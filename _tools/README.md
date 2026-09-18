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

Write a component value so that putting a minus sign in front of the whole string negates it.
The page groups a tensor's components by value and merges a value with its negation, printing `R^t{}_{\theta t\theta} = -R^t{}_{\theta\theta t} = \dots` on one line, and it recognises the negation by that leading minus alone.
So a value that is a bare sum wants collecting over a common denominator or wrapping in `\left(\right)` first, and its opposite wants writing as that string with a `-` in front rather than with the signs distributed through it.
Morris-Thorne is the worked example: its twenty four Riemann components print on six lines fully lowered and its twenty four Weyl components on six, and written the other way each would take twenty four.

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

It also checks that every term of every published expression carries the dimensions of its left hand side, which needs no algebra at all.
That pass names the file, the system, the field and the term it could not balance, and it is what catches a component printed in the bare chart rather than in $x^0 = cT$.

sympy is not installed system wide, and the virtual environment does not belong in the repository:

    python3 -m venv /tmp/mfs-venv && /tmp/mfs-venv/bin/pip install sympy
    /tmp/mfs-venv/bin/python _tools/derivations/verify_metrics.py

The whole collection takes about a quarter of an hour. `--system <metric_id>/<system_id>` checks one system and takes seconds, which is what to use while editing a single entry.
`--budget <seconds>` changes how long sympy may spend on one tensor.
Kerr is the exception to the seconds, and so far it is the exception to the check.
Its metric is not diagonal and its components are rational functions with $(r^2+a^2\cos^2\theta)$ to a high power underneath, and `norm` calls sympy's general `simplify` on every one of them.
Its Christoffel symbols come out in about forty seconds; its Riemann tensor has run for over ninety minutes of processor time without finishing, and its Kretschmann scalar is further out again, so raising `--budget` does not rescue it.
`_tools/derivations/kerr.md` Step 17 records what was checked in place of a full pass, and names the change to `norm` that would bring Kerr inside the budget.
`--dimensions-only` runs the dimensional pass alone, which takes about a second over the whole collection, so there is no reason not to run it on every edit.

Neither pass exits clean on the collection as it stands.
The sympy pass reports the disagreements `derivations/audit-2026-09-18.md` counts, and the dimensional pass reports seven terms, in `frw` and in `stockum_dust`.
Both print what failed on stderr and print the single line saying nothing failed on stdout, so a run piped through `2>&1 | tail` can look clean when it is not.
Read the exit code, and run the script before and after a change so that a failure you did not cause is not mistaken for one you did.

`derivations/audit-2026-09-18.md` groups and counts the disagreements the collection carried when the dimensional pass was added, and says which of them are the checker's fault rather than the physics'.
It is a dated snapshot of 433, of which the 320 that were a Weyl block copied from the entry's own Riemann block have since been corrected, leaving 113.
`derivations/weyl.md` is the working behind those corrections, and is the thing to read before touching any `weyl_tensor`: a Weyl tensor equals Riemann only in a vacuum, it can be nonzero in a slot where Riemann vanishes, and two of the entries that publish one are conformally flat and so publish nothing.

Nothing is ever passed in silence. A value that cannot be parsed, a system with no time coordinate declaration, or a tensor sympy cannot finish in the budget is reported as `UNCHECKED` with the reason, separately from the disagreements.

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
`a = a(t)` declares a function of one coordinate, which answers to a dot and to a prime, and `H = H(u,x,y)` declares a function of several, which answers to `\partial`: a published `\partial_x H`, `\partial_x^2 H` or `\partial_x\partial_y H` is read as that first or second partial derivative, which is as far as any curvature tensor reaches.
Every such derivative is taken with respect to the chart coordinate, exactly as the dot and the prime are, so one along a time carries a factor of $1/c$ and one along a length does not.
`_tools/derivations/pp_wave.md` is the worked example, and its Step 15 is where that factor earns its keep.

The prime works only on a function whose name is not a LaTeX command.
The reader turns a prime into a name suffix before it turns `\Phi` into a name, so `b'` reads as the derivative of `b` while `\Phi'` leaves a stray backslash and is rejected as unhandled LaTeX.
`\partial` works on either kind of name, so an entry with a Greek named function writes every derivative that way; `_tools/derivations/morris_thorne.md` Step 17 says so of its own `\Phi` and `b`.

An entry whose parameters are not free needs a second line, in `PARAMETER_RELATIONS`.
Kasner prints three exponents bound by $\sum_i p_i = \sum_i p_i^2 = 1$ and claims its values only on the surface those two equations cut out; its Ricci tensor is zero there and nowhere else.
Checked against free symbols it would report seventy four disagreements that are not disagreements.
The table carries a rational parametrisation of the surface instead, and every constrained parameter is replaced by its parametrised value on both sides of each comparison, so the identity checked is the one the entry asserts.
That is exact rather than a sample, because a rational parametrisation covers a dense subset of an irreducible variety.
`_tools/derivations/kasner.md` derives the parametrisation the Kasner entry uses and why the constraints are the field equations.
