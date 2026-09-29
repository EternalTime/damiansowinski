# Figures drawn in three dimensions on the spacetimes page

Every drawing on the spacetimes page, `_layouts/mfs.html` with the data under `MFS/assets/data/`, sorted by whether it has depth.
Surveyed on 29 September 2026 at commit 888485d, before the change that removes every glow.

## Turned by the reader already

### The embedding diagrams

Thirty-one views in twenty-nine files, `MFS/assets/data/embedding/<metric_id>.json`, every one of them carrying `figure.turn`.
`_tools/derivations/embedding.py` writes each surface in three dimensions under `surfaces` (a profile of $\rho$ against $z$ for a surface of revolution, a grid of heights over a plane for Alcubierre's and Krasnikov's, marked circles, curves and points) and projects it once from a fixed camera into `figure.layers`.
`emFigure()` draws the published layers; `wireEmbedding()` redraws them with `MfsTurn.draw()` from `MFS/assets/embedding-turn.js` once the reader drags.
The data holds the three dimensional geometry.

What the reader can do to one:

- A drag across turns it round its axis, the drawing's whole width being half a turn, as far as the reader likes.
- A drag up or down tilts it, from looking straight down the axis, elevation $90°$, to looking straight up it, $-90°$, and never past, so the axis always stands up the page.
- A finger that sets off across the drawing turns it; one that sets off up or down scrolls the page (`touch-action: pan-y`), and two fingers scroll a drawing wider than its frame.
- The arrow keys turn it by $15°$ once the drawing has the focus.
- A double click, a double tap, Home, Escape or the reset button in the corner of the frame, shown only while the figure is turned, bring back the published figure.
- While a drag goes on the drawing is made coarsely, and in full when it ends.
- Print shows the published figure, since the print copy is made from the page as first drawn.

## To be turned: the light cones in three dimensions

Five figures under `projections` in `MFS/assets/data/diagrams/<metric_id>.json`, each a slice of one time and two spatial coordinates with every other coordinate held fixed, placed in Euclidean $(X, Y, T)$ with $T$ up and seen from azimuth $-90°$ and elevation $30°$.

| spacetime | coordinates | `id` | button | cones | lines besides the cones | labels |
| --- | --- | --- | --- | --- | --- | --- |
| Gödel | `cylindrical` | `tipping` | light cones about the axis | 13 | floor (12 spokes, circles $r_c/2$ and $2r_c$), `critical` $r_c$, `ctc` $3r_c/2$ and its 4 arrows, `axis` | $t$, $r = r_c$, $r = 3r_c/2$ |
| van Stockum | `cylindrical` | `tipping` | light cones about the axis | 13 | the same, at $R$ and $3R/2$ | $t$, $r = R$, $r = 3R/2$ |
| Kerr | `boyer_lindquist` | `dragging` | light cones on the equator | 12 | floor (12 spokes from $r_+$, circles $3r_E/2$ and $7r_E/4$), `horizon` $r_+$, `ergo` $r_E$, `axis` | $t$ |
| Kerr-Newman | `boyer_lindquist` | `dragging` | light cones on the equator | 12 | the same | $t$ |
| Alcubierre | `cartesian` | `bubble` | the bubble in three dimensions | 13 | floor (12 spokes, circles $R/2$, $2R$, $2.5R$), `ergo` $v_sf = 1$, `axis`, `world` $x = 2ct$ | $t$, $x = 2ct$ |

Each also carries one slice of its embedding diagram on its floor, the moment $t = 0$: a region, and for Gödel and van Stockum the rim where the embedding stops short of $r_c$.

### How each is drawn

`_tools/derivations/projections.py` builds every figure in $(X, Y, T)$ from the published metric: each future light cone as its apex and a rim of 96 generators of one Euclidean length, each checked null.
Its `Figure` class projects every piece through `Camera(-90, 30)` the moment it is added and keeps only the result:

- a line is projected, thinned by Ramer-Douglas-Peucker at 0.0005 of the page's units, and rounded to four decimals;
- a cone becomes the convex hull of its projected apex and rim, filled, its closed rim, 8 of its generators from the apex to the rim, the two generators that bound the hull where the apex lies on it, and the apex as a point;
- the cones are painted farthest first by their apexes' depth at that one camera, after every other line;
- a label stands at the projection of a point of the slice;
- a slice's region and rim are projected and thinned as a line is.

`null_rays.py` writes the figure into `projections`, and `pjFigure()` and `pjSvg()` in `_layouts/mfs.html` paint the layers in order into an SVG in the conformal diagram's frame, with the TeX labels laid over it by `cdLabels()`.
Nothing on the page can turn them.

### What the data holds

The projection alone.
`layers`, `labels`, `slices` and `box` are all in the plane of the page at the one camera.
The three dimensional geometry, the apex and rim of every cone, the floor, the circles, the arrows, the axis, the world line, the points the labels hang from and the slice's floor, exists only inside `projections.py` while it runs.
The generator has to publish it for any reader, the page or the application, to draw the figure from another side.

## Not turned, and why

### Light passing the cosmic string

`cosmic_string.json`, `projections.conical`, `beam`: the plane $z = 0$ seen from straight overhead with $t$ left out, unrolled onto the page by the angle $(1 - 4G\mu/c^2)\phi$, in which the plane is flat.
Every point of it lies at $T = 0$ and the page's coordinates are the plane's own flat coordinates, so the figure has no depth, and from any other side it would only be the same map foreshortened.
The cone that plane is, seen in three dimensions, is the cosmic string's embedding diagram, which turns already.

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

## What turning the five light-cone figures takes

- `projections.py` records every piece of a figure in $(X, Y, T)$ as well as projecting it, and writes it into a `turn` block beside the layers: each line's points and class in the order painted, each cone's apex and rim, the point each label stands at and each slice's region and rim, all at the precision of the published figure.
  The published `layers`, `labels`, `box` and `slices` stay as they are, byte for byte.
- `MFS/assets/embedding-turn.js` becomes one module that turns both kinds of figure, sharing the camera, the Ramer-Douglas-Peucker thinning and the fitting of a turned figure into its box, with the light cones drawn by the generator's own rules: projection, thinning, hull, ribs and the farthest cone first.
- A turned light-cone figure keeps the box it was published in, drawn smaller at one scale only as far as a side needs more room than it had at the start and moved only as far as that room needs, as a turned embedding diagram is.
- `wireEmbedding()` becomes the one controller for every figure that turns, so the drag, the tilt limits, the keys, the double click, the double tap, the reset button and print are the same on both.
- `_tools/README.md` defines the new fields and the rules for drawing them at another camera, for the application.
- The tests hold every figure to giving back its published drawing at its own camera, to staying inside its box at every other, and to turning under a drag, for an embedding surface of revolution, an embedding height over a plane and a light-cone figure.
