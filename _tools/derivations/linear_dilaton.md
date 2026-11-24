# Linear dilaton vacuum

## Coordinates and sources

The two Einstein-frame charts use Cartesian spatial coordinates and natural units.
H. J. de Vega and I. L. Egusquiza, arXiv:hep-th/9309016, equation (3.2), give the linearly expanding metric in cosmic and conformal time and credit Robert C. Myers.
Roberta Angius, Matilda Delgado and Angel M. Uranga, arXiv:2207.13108, equations (3.1)-(3.8), give the flat string frame, linear dilaton, conformal transformation and cosmic time in the supercritical theory.
With the four-dimensional factor $g_E=e^{-2\Phi}g_S$ and $\Phi=-q\eta$, the relation $\tau=e^{q\eta}/q$ pulls the conformal metric back to the cosmic metric.
The coordinate selector holds Einstein-frame charts of one metric.
The string-frame metric is flat, while the Einstein-frame metric has nonzero curvature, so the conformal transformation changes the metric itself and cannot be included as another coordinate chart of that Einstein-frame geometry.
The flat string-frame line element and the conformal factor appear in the History, where the reader can compare both frames.
The sourced early boundary is not assigned a causal character here: de Vega and Egusquiza and Angius and coauthors describe it differently, and the distinction has not been resolved.

## Drawings

Both spacetime diagrams trace null geodesics of the published Einstein-frame metrics at $q=1$, on $y=z=0$.
The cosmic-time null coordinates are $\ln(q\tau)\mp qx$; the conformal-time null coordinates are $q(\eta\mp x)$.
Their marked moments are the four key moments of the embedding diagram, $q\tau=1/2$, $1$, $3/2$ and $2$.
The conformal views compactify the same plane with the arctangent of its two conformal null coordinates; the past boundary has divergent Einstein-frame curvature and the future null boundary has vanishing curvature along the sampled limits.
The curvature figure reads the stored scalars $R=6/\tau^2$ and $K=12/\tau^4$ from the cosmic chart.
All four drawing kinds are present.

`view_text.length` counts 124 characters on each spacetime view, 179 on each conformal view, 391 on the embedding's world tube, 311 on its ring, and 142 on the curvature view.

## Embedding diagram

### What the metric has

A moment of cosmic time is $q^2\tau^2(dx^2+dy^2+dz^2)$ in the cosmic chart, so every moment is flat, and its plane $z=0$ is Euclid's plane with every length scaled by $q\tau$.
A moment drawn alone is therefore the picture of empty flat space.
The conformal chart has the same moments: its plane at $\eta=\ln(q\tau)/q$ is $e^{2q\eta}(dx^2+dy^2)$, the same scale factor $q\tau$.
The dilaton $\Phi=-q\eta=-\ln(q\tau)$ takes one value over a whole moment, so a height of the dilaton over the plane would be level, and so would a height of either curvature scalar, which the curvature figure already plots through cosmic time.
What changes is the scale factor, linear in cosmic time.
The published Christoffel symbols of both charts have no $\Gamma^i{}_{\tau\tau}$ and no $\Gamma^i{}_{\eta\eta}$, so a free particle at rest in either chart stays at rest, and a ring of them on $x^2+y^2=1/q^2$ is at every moment the circle of proper radius $q\tau\cdot(1/q)=\tau$.
It grows at a steady rate, and its radius vanishes as $\tau$ falls to zero, where the curvature diverges.
The comoving distance light covers, $\int d\tau/(q\tau)$, diverges toward the past and toward the future, so the observer at the ring's centre has neither a particle horizon nor an event horizon.

### The decision

The embedding diagram has two views in the form of the Big Rip's ring of galaxies: the ring's world tube, and beside it the flat moments as a movie, both of the plane $z=0$ of the cosmic chart at $q=1$, with $1/q$ the unit of $\tau$ and of every length.
The ring of free particles stands at rest on $x^2+y^2=1/q^2$, the circle of radius $1/q$ at $q\tau=1$, where the scale factor is one, so its proper radius is $\tau$ itself.
The movie runs from $q\tau=1/2$ to $2$, a frame every $1/40$ of $q\tau$, 61 frames, with the key moments $q\tau=1/2$, $1$, $3/2$ and $2$, over which the ring grows fourfold by the same step from frame to frame.
Each moment is the disc of proper radius $2.5/q$ about the ring's centre, which holds the ring at its widest, $2/q$, with twelve particles marked every $30°$ and lines from the centre every $30°$.
The world tube stacks the ring at every $\tau$ from the first moment to the last, each circle at the height $2\tau$, with the world lines of the twelve particles and the axis of time.
Since the radius and the height are both proportional to $\tau$, every row lies on one cone, whose apex on the axis is the moment $\tau=0$, the curvature singularity; the tube is the part of that cone between $q\tau=1/2$ and $2$, and its world lines are straight lines through the apex.
The tube's caption says so, and the ring's caption gives the plane, the ring, the span, and the proper radius.
No quantity is drawn as a height over the plane, since every moment is homogeneous.

### What the generator checks

`embedding_views` in `linear_dilaton_figures.py` reads each moment as a `FlatPlane` of the published cosmic metric and checks its scale against $q\tau$, and reads the conformal chart's plane at $\eta=\ln(q\tau)/q$ and checks that its scale is the same.
`at_rest` checks the published Christoffel symbols of each chart for the missing $\Gamma^i{}_{\tau\tau}$ and $\Gamma^i{}_{\eta\eta}$.
Each frame is held to the published metric along, across, and around, as every slice is, its ring to the circle of proper radius $\tau$, the ring's growth from frame to frame to the step of $\tau$, and the disc to holding the ring at its widest.
The tube's rows are checked against the ring at every frame and at the four key moments, its triangles against its ellipses in space, and every row against the cone $\rho=z/2$ through $\tau=0$.
`test_embedding_is_a_ring_of_free_particles_of_proper_radius_tau` in `_tools/test_linear_dilaton.py` holds the file to the same, from its numbers and the published metric.

### The moments on the other drawings

Each key moment is the whole plane $z=0$, so on both spacetime diagrams it is the line of constant time across the drawing, $\tau$ in the cosmic chart and $\eta=\ln(q\tau)/q$ in the conformal chart.
On both conformal views it is the curve $p,v=\arctan(\eta\mp x)$ of that $\eta$ over every $x$.
Each is marked twice, once as a ring of the tube and once as a moment of the movie, as the Big Rip's moments are.
