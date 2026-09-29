# Figures drawn in three dimensions on the spacetimes page

Every drawing on the spacetimes page, `_layouts/mfs.html` with the data under `MFS/assets/data/`, sorted by whether it has depth, with how each is drawn and what its data holds.
The captain asked on 29 September 2026 for every figure in three dimensions to turn as the embedding diagrams do, the light cones about Gödel's axis among them.

## Turned by the reader

### The embedding diagrams

Thirty-one views in twenty-nine files, `MFS/assets/data/embedding/<metric_id>.json`, every one of them carrying `figure.turn`.
`_tools/derivations/embedding.py` writes each surface in three dimensions under `surfaces` (a profile of $\rho$ against $z$ for a surface of revolution, a grid of heights over a plane for Alcubierre's and Krasnikov's, marked circles, curves and points) and projects it once from a fixed camera into `figure.layers`.
`emFigure()` draws the published layers, and `wireTurning()` draws the view again from its surfaces with `MfsTurn.draw()` from `MFS/assets/turn.js` once the reader drags it.
The data holds the three dimensional geometry.

### The light cones in three dimensions

Five figures under `projections` in `MFS/assets/data/diagrams/<metric_id>.json`, each a slice of one time and two spatial coordinates with every other coordinate held fixed, placed in Euclidean $(X, Y, T)$ with $T$ up and published as seen from azimuth $-90°$ and elevation $30°$.

| spacetime | coordinates | `id` | button | cones | lines besides the cones | labels |
| --- | --- | --- | --- | --- | --- | --- |
| Gödel | `cylindrical` | `tipping` | light cones about the axis | 13 | floor (12 spokes, circles $r_c/2$ and $2r_c$), `critical` $r_c$, `ctc` $3r_c/2$ and its 4 arrows, `axis` | $t$, $r = r_c$, $r = 3r_c/2$ |
| van Stockum | `cylindrical` | `tipping` | light cones about the axis | 13 | the same, at $R$ and $3R/2$ | $t$, $r = R$, $r = 3R/2$ |
| Kerr | `boyer_lindquist` | `dragging` | light cones on the equator | 12 | floor (12 spokes from $r_+$, circles $3r_E/2$ and $7r_E/4$), `horizon` $r_+$, `ergo` $r_E$, `axis` | $t$ |
| Kerr-Newman | `boyer_lindquist` | `dragging` | light cones on the equator | 12 | the same | $t$ |
| Alcubierre | `cartesian` | `bubble` | the bubble in three dimensions | 13 | floor (12 spokes, circles $R/2$, $2R$, $2.5R$), `ergo` $v_sf = 1$, `axis`, `world` $x = 2ct$ | $t$, $x = 2ct$ |

Each also carries one slice of its embedding diagram on its floor, the moment $t = 0$: a region, and for Gödel and van Stockum the rim where the embedding stops short of $r_c$.

`_tools/derivations/projections.py` builds every figure in $(X, Y, T)$ from the published metric: each future light cone as its apex and a rim of 96 generators of one Euclidean length, each checked null.
Its `Figure` class projects every piece through `Camera(-90, 30)` into the published layers:

- a line is projected, thinned by Ramer-Douglas-Peucker at 0.0005 of the page's units, and rounded to four decimals;
- a cone becomes the convex hull of its projected apex and rim, filled, its closed rim, 8 of its generators from the apex to the rim, the two generators that bound the hull where the apex lies on it, and the apex as a point;
- the cones are painted farthest first by their apexes' depth, after every other line;
- a label stands at the projection of a point of the slice;
- a slice's region and rim are projected and thinned as a line is.

The same class records every piece in $(X, Y, T)$ at six decimals under `turn`, beside the layers, which stay what the projection alone gives, byte for byte.
`pjFigure()` and `pjPaths()` in `_layouts/mfs.html` paint the published layers into an SVG in the conformal diagram's frame, with the TeX labels laid over it by `cdLabels()`, and `wireTurning()` draws the figure again from `turn` with `MfsTurn.draw()` once the reader drags it.
The data holds the projection and the three dimensional geometry.

### What the reader can do to either

- A drag across turns it round its axis, the drawing's whole width being half a turn, as far as the reader likes.
- A drag up or down tilts it, from looking straight down the axis, elevation $90°$, to looking straight up it, $-90°$, and never past, so the axis always stands up the page.
- A finger that sets off across the drawing turns it; one that sets off up or down scrolls the page (`touch-action: pan-y`), and two fingers scroll a drawing wider than its frame.
- The arrow keys turn it by $15°$ once the drawing has the focus.
- A double click, a double tap, Home, Escape or the reset button in the corner of the frame, shown only while the figure is turned, bring back the published figure.
- While a drag goes on an embedding diagram is made coarsely, and in full when it ends.
- Print shows the published figure, since the print copy is made from the page as first drawn.

`MfsTurn.dragged()`, `keyed()` and `turned()` are those rules, one hand for both kinds of figure.

## Not turned, and why

### Light passing the cosmic string

`cosmic_string.json`, `projections.conical`, `beam`: the plane $z = 0$ seen from straight overhead with $t$ left out, unrolled onto the page by the angle $(1 - 4G\mu/c^2)\phi$, in which the plane is flat.
Every point of it lies at $T = 0$ and the page's coordinates are the plane's own flat coordinates, so the figure has no depth, and from any other side it would only be the same map foreshortened.
The cone that plane is, seen in three dimensions, is the cosmic string's embedding diagram, which turns.

### The flat views of the spacetime diagrams

`MFS/assets/data/diagrams/<metric_id>.json`, `systems`: null rays and future light cones on a plane of two coordinates, drawn by `nrFigure()` in the plane's own coordinates.
Kerr's and Kerr-Newman's principal null rays leave every such plane, and are drawn as their two projections onto $t$ and $r$ and onto the equator seen from above, which together fix every ray, since $r$ changes monotonically along each.
Van Stockum's and Gödel's cylinders of $t$ and $\phi$ are drawn unrolled, $\phi$ scaled by $r$.
None of them is a picture of depth.

### The conformal diagrams

`MFS/assets/data/conformal/<metric_id>.json`: the plane $(X, T)$ of two null coordinates, light at $45°$, drawn by `cdSvg()`.

### Everything else drawn on the page

The legend swatches, a cone among them, are icons 34 by 16 pixels.
The grid behind the panels, `#wavy-grid`, is a flat grid displaced under each panel, drawn on a canvas.

## How the light cones turn

- `turn` holds the figure's lines in the order painted, its cones in the order painted at the published camera, how many generators of each rim are drawn as ribs, where each label stands, each slice in three dimensions, and the centre the figure turns about, on the axis halfway up the drawing; `_tools/README.md`, "Turning a figure of light cones", defines each field for the application.
- `drawCones()` in `MFS/assets/turn.js` draws the figure at any camera by the generator's own rules: projection, thinning, hull, rim, ribs, the generators that bound the hull and the farthest cone first, two cones at the same depth keeping their published order.
- A label that names a circle about the axis, $r = r_c$ and $r = 3r_c/2$ on Gödel's figure and $r = R$ and $r = 3R/2$ on van Stockum's, stands as far round its circle from the camera as it stood at the published camera, so it keeps to the side of the circle nearest the reader, and the outer one gives way while the two overlap, within $18.5°$ of seeing the floor edge on for Gödel and $25°$ for van Stockum, whose two circles stand closer; the legend names both circles.
- A turned figure keeps the box it was published in, drawn at one scale, never above 1, as far as its height or its width on either side of the axis needs, and moved up or down only as far as its height needs, down to between 0.50 and 0.52 of its size looking straight down or straight up the axis.
- `_tools/turn_check.cjs` holds every figure to giving back its published layers at its own camera, in the published order, each cone to the published rounding and each line within its thinning, and to staying inside its box from every side; it drags a surface of revolution, a height over a plane and a figure of light cones with the page's own hand.
- `_tools/turn_drag.mjs` drags one figure of each kind in a browser with a mouse, a finger and the keys, and holds the drawing, its reset button, its labels and its print copy to all of the above.
