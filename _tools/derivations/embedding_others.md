# Embedding diagrams for the last eleven spacetimes

Eighteen spacetimes have an embedding diagram drawn by `_tools/derivations/embedding.py`: the equatorial plane of one moment, a surface of one coordinate $x$ and one angle $\phi$ with metric $g_{xx}\,dx^2 + g_{\phi\phi}\,d\phi^2$, drawn as a surface of revolution in flat space at the radius $\rho = \sqrt{g_{\phi\phi}}$ and climbing at $dz/dx = \sqrt{g_{xx} - (d\rho/dx)^2}$.
Eleven had none.
Eight were stated in sentences, because their slices are flat or hyperbolic or not spacelike everywhere: Alcubierre, anti-de Sitter, Bianchi I, Kasner, Krasnikov, Lentz, Minkowski and Natário.
Three had no file at all, because their slices were thought to carry no one surface that says anything: the Malament-Hogarth toy, the Mixmaster universe and the pp-wave.

Each of the eleven now gets a drawing, and the drawings fall into three kinds.
Where a slice is curved and has an axis, it is drawn as the other eighteen are, a surface of revolution in flat space, and where the construction runs out the drawing stops at the circle where it does: the Malament-Hogarth toy and the Mixmaster universe.
Where a slice is curved the other way at every point, as the hyperbolic plane is, no surface in flat space carries it, and it is drawn in three dimensional Minkowski space instead, where the whole of it lies on one sheet of a hyperboloid: anti-de Sitter.
Where a slice is flat, the surface is a plane, and what makes the spacetime what it is lies in how the slices are stacked, not in the shape of any one; the plane is drawn with that feature marked on it, a ring of free particles as it is stretched, the wall of a warp bubble and the flow of space through it, or the circle inside which the direction along a tube is a time.

Every number was computed from the published metric components through the checker's `Reader`, with `null_rays.load` and `published_matrix` as the other diagrams read them, and every declared function is the one a diagram of the same spacetime already declares, unless it says otherwise.
`embedding.py` recomputes every one of them and checks it before it writes a file.

---

## Group 1. Curved slices

### The Malament-Hogarth toy

The line element is $\Omega^2(-c^2dt^2 + dx^2 + dy^2 + dz^2)$ with the event at the origin removed, and its spacetime diagram declares

$$\Omega = 1 + \frac{e^{1 - 1/(1 - \varrho^2)}}{\varrho},\qquad \varrho^2 = c^2t^2 + x^2 + y^2 + z^2 < 1,$$

and $\Omega = 1$ for $\varrho \ge 1$, which grows as $1/|ct|$ along the axis toward the removed event.

The slice drawn is the plane $z = 0$ at one moment $ct = T$, about the origin.
The published metric is diagonal, so on the plane, in polar coordinates $x = s\cos\phi$, $y = s\sin\phi$, it pulls back to

$$ds^2 = \Omega^2\left(ds^2 + s^2d\phi^2\right),\qquad \Omega = \Omega\left(\sqrt{T^2 + s^2}\right),$$

which sympy returns from the published components with no cross term, and the declared $\Omega$ depends on $x$ and $y$ only through $s$.
So the plane is a surface of revolution about the origin with $g_{ss} = \Omega^2$ and $\rho = s\,\Omega$, and

$$g_{ss} - \left(\frac{d\rho}{ds}\right)^2 = \Omega^2 - (\Omega + s\Omega')^2 = -s\Omega'\,(2\Omega + s\Omega').$$

The first factor is positive because $\Omega$ falls outward, and the second stays positive as well: its least value over $s$ is $0.1150$ at $T = 0$, reached at $s = 0.789$, and $0.1153$, $0.1530$, $0.4478$, $0.9814$ and $1.6129$ at $|cT| = 0.01$, $0.1$, $0.3$, $0.5$ and $0.7$.
The surface therefore exists at every $s$ at every moment, a well sunk into a plane, flat beyond the circle $s = \sqrt{1 - T^2}$ where $\Omega$ reaches $1$.

At a moment $T \ne 0$ the well has a smooth bottom, since $\Omega$ is finite at $s = 0$, $1 + e^{1 - 1/(1 - T^2)}/|T|$.
At $T = 0$ it has none: $\Omega \to 1/s$ as $s \to 0$, so the circles close in on the radius $\rho = s\Omega \to 1$ while the proper distance down to them, $\int\Omega\,ds$, grows as $\ln(1/s)$.
The slice through the removed event is a funnel that runs into a tube of radius $1$ and on down it for ever, as Bertotti and Robinson's cylinder runs on and as the throat of an extreme Reissner-Nordström hole does.

The depths below the flat rim, and the proper distance from the centre out to the rim, are:

| $cT$ | rim $s$ | $\Omega$ at the centre | depth | proper distance |
|---|---|---|---|---|
| $-0.9$ | $0.4359$ | $1.0156$ | $0.0363$ | $0.4382$ |
| $-0.7$ | $0.7141$ | $1.5466$ | $0.4468$ | $0.8975$ |
| $-0.5$ | $0.8660$ | $2.4331$ | $0.9244$ | $1.4651$ |
| $-0.3$ | $0.9539$ | $4.0194$ | $1.4991$ | $2.1683$ |
| $-0.1$ | $0.9950$ | $10.8995$ | $2.6596$ | $3.3879$ |

and at $T = 0$, from the rim at $s = 1$ down to the circle at $s$:

| $s$ | $\rho$ | depth | proper distance |
|---|---|---|---|
| $0.3$ | $1.205832$ | $1.1768$ | $1.3632$ |
| $0.1$ | $1.089950$ | $2.4290$ | $2.6208$ |
| $0.03$ | $1.029100$ | $3.6968$ | $3.8902$ |
| $0.01$ | $1.009900$ | $4.8148$ | $5.0084$ |
| $0.001$ | $1.000999$ | $7.1263$ | $7.3200$ |

The computer of Earman and Norton's toy rides the axis toward the removed event, and its clock reads $\int\Omega\,c\,dt$ along the axis.
From $ct = -1$ it reads $0.3470$ at $ct = -0.7$, $1.3632$ at $-0.3$, $2.6208$ at $-0.1$, $5.0084$ at $-0.01$ and $7.3200$ at $-0.001$, and without limit as $ct \to 0$.
At $T = 0$ the tube from its rim down to $s = 0.3$, $0.1$, $0.01$ and $0.001$ is $1.3632$, $2.6208$, $5.0084$ and $7.3200$ long: because $\Omega$ depends only on $c^2t^2 + x^2 + y^2 + z^2$, the tube of the slice $T = 0$ from its rim down to the circle $s$ is exactly as long as the computer's clock runs from $ct = -1$ to $ct = -s$.
The infinite proper time the computer spends reaching the removed event, the property the spacetime is named for, is the infinite length of the tube.

The drawing is a sequence of four moments, $cT = -0.7$, $-0.3$, $-0.1$ and $0$, each out to $s = 1.5$, where the plane is flat, with the last drawn down the tube to $s = 0.03$ and ended there with an edge, the tube running on.
The unit is the radius of the region where $\Omega \ne 1$, the unit the declared $\Omega$ is written in.
The edge of that region is $\sqrt{1 - T^2}$ correctly rounded, moved on by a unit in the last place where the declared function's two branches disagree there: within one unit of the edge the float $T^2 + s^2$ can fall below $1$ while $s^2 - (1 - T^2)$ rounds above $0$, and the inside branch then divides by a number of the wrong sign, which happens at $cT = -0.7$ and at no other moment drawn.
Nothing is checked against a field equation, since the spacetime has none, its stress-energy being defined as its Einstein tensor; the checks are the isometry of every piece, the flat rim, the closed form $\rho \to 1$ of the tube and the equality of the tube's length with the axis clock at the four moments.

### Anti-de Sitter

The static slice has, on its equator, $g_{rr} = 1/(1 + r^2/L^2)$ and $g_{\phi\phi} = r^2$, so $d\rho/dr = 1$ and $g_{rr} - (d\rho/dr)^2 = -r^2/(L^2 + r^2) < 0$ at every $r > 0$.
Every circle about the centre grows faster than the distance out to it, and no surface of revolution about the centre exists in flat space, from the centre outward.
The slice is the hyperbolic plane of curvature $-1/L^2$: with $r = L\sinh(\sigma/L)$ it is $d\sigma^2 + L^2\sinh^2(\sigma/L)\,d\phi^2$.

In three dimensional Minkowski space, with $dX^2 + dY^2 - dZ^2$, a profile $(\rho(r), Z(r))$ turned about the $Z$ axis has the metric $\left((d\rho/dr)^2 - (dZ/dr)^2\right)dr^2 + \rho^2d\phi^2$, so the construction runs with the sign of the defect turned over,

$$\frac{dZ}{dr} = \sqrt{\left(\frac{d\rho}{dr}\right)^2 - g_{rr}} = \frac{r}{\sqrt{L^2 + r^2}},\qquad Z = \sqrt{L^2 + r^2} - L,$$

which sympy returns from the published components.
The surface is one sheet of the hyperboloid $(Z + L)^2 - X^2 - Y^2 = L^2$, and the whole hyperbolic plane lies on it, every distance along it, measured with $dX^2 + dY^2 - dZ^2$, the metric distance.
Its tangent planes are spacelike everywhere, and it approaches the light cone $Z + L = \rho$ of the flat space it is drawn in without reaching it, so the conformal boundary at $r \to \infty$ lies along that cone at infinity.
Wilhelm Killing in 1880 and Henri Poincaré in 1881 each wrote the hyperbolic plane on this sheet.

The drawing is the static slice's equator out to $r = 4L$, the range of its spacetime diagrams, with circles of constant $r$ at $L$, $2L$, $3L$ and $4L$, and the light cone drawn faintly as a reference, down to its apex at $Z = -L$.
A distance on the drawing is not the distance a reader's eye measures: where the sheet is steep, $dZ$ nearly cancels $d\rho$, and the step from $r = 3L$ to $4L$, which looks longer than the step from $0$ to $L$, is $\mathrm{arcsinh}\,4 - \mathrm{arcsinh}\,3 = 0.2763\,L$ against $0.8814\,L$.
The caption and the view say so.

The rounding of the written surface needs care there.
At $r = 4L$ the chord of the profile has $\sqrt{d\rho^2 - dZ^2} = 0.2425\,d\rho$, so a chord a ninetieth of the drawing's width of $8L$ is $0.0218\,L$ long in the metric, and one half as long, as the halving next to a bend leaves some, is $0.0108\,L$.
Rounding $\rho$ and $Z$ to $2 \times 10^{-7}$, half the last place the file keeps for a piece $4L$ across, moves those chords by up to $7.4 \times 10^{-5}$ and $1.5 \times 10^{-4}$ of themselves, most of the $2 \times 10^{-4}$ the chords are held to.
A surface in Minkowski space is therefore written to nine decimals, which brings both below $2 \times 10^{-6}$.

The Poincaré patch's slice of constant $t$ is the same hyperbolic space, $L^2(dx^2 + dy^2 + dz^2)/z^2$, the upper half space, and its moment $t = 0$ is the global moment $t = 0$, since both are the surface fixed by $t \to -t$; so the one drawing serves both charts.
The checks are the isometry in Minkowski space of every chord along the profile and every chord across it, $\rho = \sqrt{g_{\phi\phi}}$, and the closed form $Z = \sqrt{L^2 + r^2} - L$.

### The Mixmaster universe

The slice of constant $t$ is a three sphere with metric $a_1^2(\omega^1)^2 + a_2^2(\omega^2)^2 + a_3^2(\omega^3)^2$ in the Euler angles $(\psi, \theta, \phi)$, and when all three scale factors are equal to $a$ it is the round sphere of radius $2a$, whose equator is a round two sphere of radius $2a$.

The two sphere drawn is the great sphere, the unit sphere of a three dimensional subspace when the three sphere is the unit quaternions.
Left translations are isometries of every Mixmaster slice, and they carry any great sphere onto any other, so at any moment every great sphere is congruent to every other and the great sphere is one surface.
Through the identity, the great sphere of the quaternions with no $k$ part is $\psi + \phi \equiv 0 \pmod{2\pi}$ in the Euler angles: two hemispheres, $\psi = -\phi$ and $\psi = 2\pi - \phi$, each with $\theta \in [0, \pi]$, meeting at $\theta = \pi$.
Pulled back from the published components along $\psi = -\phi$, sympy gives

$$g_{\theta\theta} = a_1^2\cos^2\phi + a_2^2\sin^2\phi,\qquad g_{\theta\phi} = \tfrac{1}{2}\left(a_2^2 - a_1^2\right)\sin 2\phi\sin\theta,$$

$$g_{\phi\phi} = \left(a_1^2\sin^2\phi + a_2^2\cos^2\phi\right)\sin^2\theta + a_3^2(1 - \cos\theta)^2,$$

and the same on the other hemisphere.
Three great circles cross the sphere at right angles: the meridian $\phi = 0$, which runs on through the other hemisphere at $\phi = \pi$, with $g_{\theta\theta} = a_1^2$ all the way round; the meridian $\phi = \pi/2$, with $g_{\theta\theta} = a_2^2$; and the equator $\theta = \pi$, with $g_{\phi\phi} = 4a_3^2$.
Their circumferences are $4\pi a_1$, $4\pi a_2$ and $4\pi a_3$: the three scale factors are the three principal circumferences of the great sphere divided by $4\pi$.

When all three differ, the metric depends on $\phi$ as well as $\theta$, so no turn about the axis through the identity carries the sphere onto itself, and no surface of revolution about that axis carries it.
When two are equal, $a_1 = a_2 = a$ and $a_3 = b$, the cross term vanishes and the rest loses $\phi$:

$$ds^2 = a^2d\theta^2 + \left(a^2\sin^2\theta + b^2(1 - \cos\theta)^2\right)d\phi^2,$$

a surface of revolution about the axis through the identity and its antipode, with $\phi$ closing after $2\pi$ on each hemisphere.
Its poles are the identity and its antipode, the points $\theta = 0$ of the two hemispheres, and its equator $\theta = \pi$ is a fibre of $\psi$, the circle of quaternions in the plane of $i$ and $j$, of circumference $4\pi b$, where the tangent is vertical and the two hemispheres join with one tangent; each meridian from pole to pole is $2\pi a$ long.

With $\rho^2 = a^2\sin^2\theta + b^2(1 - \cos\theta)^2$,

$$g_{\theta\theta} - \left(\frac{d\rho}{d\theta}\right)^2 = \left(1 - \frac{3b^2}{4a^2}\right)a^2\theta^2 + O(\theta^4)$$

near the pole, so the curvature of the great sphere there is $(1 - 3b^2/4a^2)/a^2$, $1/(2a)^2$ when $b = a$ as the round sphere's is.
For $b < 2a/\sqrt{3}$ the defect is positive everywhere and the whole sphere is drawn; for $b > 2a/\sqrt{3}$ it is negative near both poles, the curvature there is negative, no surface of revolution carries the sphere about its poles, and the drawing is the band about the equator between the two circles where the defect vanishes.

The moments are those of Abraham Taub's universe, the vacuum member of the family with $a_1 = a_2$, which Taub solved in closed form in 1951: in the Taub-NUT form, with Taub's time $T$,

$$a_1 = a_2 = \sqrt{T^2 + l^2},\qquad a_3 = 2l\sqrt{U},\qquad U = \frac{-T^2 + 2mT + l^2}{T^2 + l^2},\qquad c\,d\tau = \frac{dT}{\sqrt{U}},$$

with $\tau$ the proper time the published chart's $t$ is, at $m = 1$ and $l = m/2$ as Taub-NUT's spacetime diagram declares.
Substituted into the published Einstein tensor, with $a_i' = \sqrt{U}\,da_i/dT$ and $a_i'' = \sqrt{U}\,d(\sqrt{U}\,da_i/dT)/dT$, every component, $G^t{}_t$, the diagonal and the off diagonal $G^\psi{}_\phi$, $G^\psi{}_\theta$ and $G^\theta{}_\phi$, simplifies to zero in sympy.
This is the cosmological region that Taub-NUT's own embedding diagram stops at, across its horizon, where $r$ is a time: the same universe, now drawn.

It lives between $T_- = m - \sqrt{m^2 + l^2} = -0.118034$ and $T_+ = 2.118034$, where $a_3 = 0$ and the fibres close up, a proper time $c\tau = 3.821549\,m$ in all.
$b/a$ rises from $0$ through $1$ and $2/\sqrt{3}$ almost at once, by $T = -0.078211$, to its largest, $2.6938$ at $T = 0.1934$, falls back through $2/\sqrt{3}$ at $T = 0.841721$ and through $1$ at $T = 0.930523$, and returns to $0$.
The drawing is a sequence of five moments:

| $T$ | $c\tau$ | $a$ | $b$ | $b/a$ | drawn |
|---|---|---|---|---|---|
| $-0.1$ | $0.0922$ | $0.509902$ | $0.392232$ | $0.7692$ | the whole sphere |
| $0.2$ | $0.3937$ | $0.538516$ | $1.450327$ | $2.6932$ | the band $\theta \ge 2.731485$ on each side |
| $0.930523$ | $0.9479$ | $1.056349$ | $1.056349$ | $1$ | the round sphere of radius $2a$ |
| $1.5$ | $1.6422$ | $1.581139$ | $0.632456$ | $0.4$ | the whole sphere |
| $2.0$ | $2.8304$ | $2.061553$ | $0.242536$ | $0.1176$ | the whole sphere |

in units of $m$.
At $T = 0.2$ the band reaches $0.2208\,m$ of proper distance from the equator on each side, where its circles have radius $2.7887\,m$ against the equator's $2.9007\,m$.
At $T = 2$ the fibres are nearly closed, and the great sphere is two lobes, each nearly a sphere of radius $a$, joined by a waist of radius $2b = 0.4851\,m$: as $b \to 0$ the three sphere becomes the two sphere of radius $a$ that its fibres are drawn over, and the great sphere covers it twice.
Each $c\tau$ is the quadrature of $dT/\sqrt{U}$ from $T_-$, which converges there, since $U$ vanishes as $T - T_-$.

The checks are the isometry of each hemisphere, the join of the two with one tangent at the equator, $\rho = 2b$ there and the meridian's length $2\pi a$, the closed form of the round moment, a sphere of radius $2a$, the negative defect between each pole and the band at $T = 0.2$, and Taub's scale factors against the published field equations.
The view states that in every other Mixmaster universe the three scale factors differ, so its great sphere has no axis and no surface of revolution carries it, and that at a moment when two are equal and the third is more than $2/\sqrt{3}$ times them, as at $T = 0.2$, the great sphere curves negatively about its poles.

---

## Group 2. Flat slices, with what the spacetime does marked on them

On each of these the metric of the plane drawn has constant coefficients and no cross term, which is checked from the published components, so the plane is flat and its embedding is a flat disc.
A disc is a surface of revolution about any point of it, and is drawn about the point the spacetime turns on, the centre of the bubble or of the tube, with its profile running along one coordinate axis of the chart from that point, $\rho = \sqrt{g_{xx}}\,x$ and $z = 0$, and a circle or curve marked on it where the physics is.
The marks that are not circles about the centre are curves on the plane, given as points of the plane in the drawing's own coordinates.

### Minkowski space

The equator of a slice of constant $t$ in the spherical chart has $g_{rr} = 1$ and $g_{\phi\phi} = r^2$, so $\rho = r$ and $g_{rr} - (d\rho/dr)^2 = 0$: the plane, on which a circle of radius $r$ has circumference $2\pi r$.
It is the surface every other embedding diagram is measured against, and it is drawn as the others are, from the published spherical chart, out to $r = 4\ell$ as its spacetime diagram runs, with circles of constant $r$ at $\ell$, $2\ell$, $3\ell$ and $4\ell$.
Minkowski space has no length of its own, so $\ell$ is any length.
Every slice of constant $t$ of the Cartesian chart is the same flat space, and the plane is the equator of each; the checks are the isometry of the profile and $g_{rr} = (d\rho/dr)^2$ at every point, which the stated file already made.

### The Krasnikov tube

At constant $t$ the published metric is $k\,dx^2 + dr^2 + r^2d\phi^2$, with $k$ entering only $g_{xx}$ and $g_{tx}$.
The spacetime diagram declares the tube along $x$ from $0$ to $D = 4$, built by a ship that left $x = 0$ at $t = 0$ at the speed of light:

$$k = 1 - (2 - \delta)\,S\!\left(\frac{\rho_0^2 - r^2}{2\rho_0}\right)S(ct - x)\,S(x)\,S(D - x),\qquad S(w) = \frac{1 + \tanh(w/0.15)}{2},$$

with $\delta = 0.2$ and $\rho_0 = 1$.
Where $k < 0$ the direction along the tube is timelike at constant $t$, so the three dimensional slice of constant $t$ is not a moment of space inside the tube.
The surface of constant $t$ and $x$ across the tube is: its metric is $dr^2 + r^2d\phi^2$, whatever $k$ is, spacelike and flat inside the tube as outside.

The drawing is that cross-section halfway along the tube, $x = D/2 = 2$, at $ct = 3$, after the ship has passed, a flat disc out to $r = 2$.
There $k = -0.79771$ on the axis and $1$ far out, and it passes through $0$ at $r = 0.983122$, with $k = -1/2$ at $r = 0.8710$ and $k = 1/2$ at $r = 1.0693$.
The circle $k = 0$ is marked: inside it the direction along the tube, square to the drawing, is a time at constant $t$, which is how the tube carries a traveller back to within an instant of departure.
Before the ship passes, at $ct = 1$, $k = 1$ everywhere on the same disc, and the disc is the same; the tube is in the direction the drawing leaves out and in the time at which it is drawn, which the caption says.
The checks are the isometry of the disc, the published $g_{rr} = 1$ and $g_{\phi\phi} = r^2$ independent of $k$, and $k$ vanishing on the marked circle.
Since 28 September 2026 the figure draws $1 - k$ as a height over the plane of the tube's axis instead, which `krasnikov_height.md` derives.

### Alcubierre's warp drive

The published metric has unit lapse and flat slices, $dx^2 + dy^2 + dz^2$ at constant $t$, and the shift $v_sf$ along $x$.
The spacetime diagram declares $v_s = 2$ and Alcubierre's profile

$$f = \frac{\tanh\sigma(r_s + R) - \tanh\sigma(r_s - R)}{2\tanh\sigma R},\qquad r_s = \sqrt{(x - 2ct)^2 + y^2 + z^2},$$

with $R = 1$ and $\sigma = 4$, which centres the bubble on $x = 0$ at $t = 0$.
The drawing is the plane $z = 0$ at $t = 0$, the plane of the ship's path, a flat disc of radius $3R$ about the ship, as far as its spacetime diagram runs.

The path runs along the drawing's $Y$, away from the reader at the start, and the profile along $y$, so that the ends of the circles on the page, where their labels stand and stay as the drawing turns, lie clear of the path and of what is marked along it.
Marked on it, first, the circle $v_sf = 1$, which for $v_s = 2$ is $f = 1/2$, at $r_s = 1.000168\,R$: the zero of the published $g_{tt}$, the circle the figure in three dimensions draws its cones round, inside which space is carried past the chart faster than light.
The wall where $f$ falls from $0.9$ to $0.1$ lies between $r_s = 0.7262\,R$ and $1.2747\,R$.

Second, the expansion of the observers who ride the slices, $n^\mu = (1, v_sf, 0, 0)$ from the published inverse metric, whose volume changes at the rate

$$\theta = \nabla_\mu n^\mu = v_s\,\partial_xf = v_s\frac{x - x_s}{r_s}\frac{df}{dr_s},$$

which sympy returns from the published metric with $\sqrt{-g} = 1$.
On the plane it is $v_s\cos\varphi\,f'(r_s)$, with $\varphi$ the angle from the direction of travel, so space contracts ahead of the ship and expands behind it, most, $4.0027/R$, on the path at $r_s = R$.
The curves $\theta = \pm 2.0013/R$, half that, are marked: two crescents in the wall, one ahead and one behind, each crossing the path between $r_s = 0.7797\,R$ and $1.2203\,R$ and reaching $60°$ either side of it, where $\cos\varphi = 1/2$ at the steepest point of $f$.
This is the picture Miguel Alcubierre drew in 1994 as a surface of $\theta$ over the plane; the embedding diagram draws the plane as the flat surface it is, with $\theta$ marked on it.

The energy density the same observers measure is, from the published $G^{tt}$,

$$\frac{c^4}{8\pi G}G^{tt} = -\frac{c^4}{8\pi G}\frac{v_s^2}{4}\frac{y^2 + z^2}{r_s^2}\left(\frac{df}{dr_s}\right)^2,$$

checked on the plane against the declared profile to $3 \times 10^{-15}$: negative, zero on the path and largest at the sides of the wall, a ring about the direction of travel.
On a flat slice that density is all extrinsic, $16\pi G\rho/c^4 = K^2 - K_{ij}K^{ij}$, since the scalar curvature of the slice vanishes; the drawing's flatness and the negative energy are two faces of one fact.
Since 28 September 2026 the figure draws $\theta$ itself as a height over this plane, as Alcubierre did, with the path along the drawing's $X$, which `alcubierre_expansion.md` derives.

### Natário's warp drive

The published metric is Natário's flow chart, flat slices with the shift $(u, v, w)$, and the spacetime diagram declares $v_s = 2$, $n = f/2$ with Alcubierre's profile, and the zero expansion field

$$X = v_s\left[(2n + r_sn')\,e_x - n'\,x_r\,(x_r, y, z)/r_s\right],\qquad x_r = x - v_st,$$

whose divergence is checked to vanish.
The drawing is the plane $z = 0$ at $t = 0$, a flat disc of radius $3R$ about the ship, with the circle $r_s = R$ marked, the middle of the wall.

The field is axisymmetric about the path and divergence free, so in the plane $z = 0$ it flows along the level curves of Stokes's stream function,

$$\Psi = v_s\,n(r_s)\,y^2,$$

with $X_x = y^{-1}\partial_y\Psi$ and $X_y = -y^{-1}\partial_x\Psi$; the declared field is checked tangent to them, $X\cdot\nabla\Psi$ vanishing on the plane to $2 \times 10^{-16}$.
Inside the bubble $n = 1/2$ and the level curves are lines along the path: space there moves rigidly forward at $v_s$ with the ship.
Outside $n = 0$ and nothing moves.
In the wall $n\,r_s^2$ rises to $0.2802\,R^2$ at $r_s = 0.8837\,R$ and falls to zero, so every level curve below that closes: space runs forward through the bubble and back round it through the wall, compressed nowhere, as a fluid that cannot be squeezed.
The lines drawn are found by following the declared field itself, from where each crosses the ship's plane $x = 0$ round until it crosses it again, and each is checked to close there to $10^{-8}$ and to hold $\Psi$, taken from the declared field as $\int y\,X_x\,dy$ out from the path, to one value to $10^{-8}$ of it.
The path runs along the drawing's $Y$, as Alcubierre's did until its height was drawn.
The lines marked are the three pairs that cross the bubble at $y = \pm 0.2R$, $\pm 0.4R$ and $\pm 0.6R$, $\Psi = 0.04$, $0.16$ and $0.36$ in units of $R^2$ with $v_s = 2$, which come back through the wall at $y = \pm 1.5019R$, $\pm 1.2775R$ and $\pm 1.1111R$ on the plane $x = 0$.
José Natário built the drive in 2002 to show that the expansion Alcubierre's drive turns on is not needed: with $\theta = 0$ everywhere the ship is carried by the sliding of space, and the energy density the riding observers measure is again $K^2 - K_{ij}K^{ij}$ over $16\pi G/c^4$, now $-K_{ij}K^{ij}$ alone, negative.

### Lentz's soliton

The published metric is Lentz's, unit lapse, flat slices and the shift a gradient, $N_i = \partial_i\phi$, with $\phi$ left free, and the flat slices are one of the three things Lentz fixed to define the class.
His soliton exists only as a numerical integral of the wave equation $\partial_x^2\phi + \partial_y^2\phi - (2/v_h^2)\partial_z^2\phi = \rho_h$ over his rhomboid cells of source, so no member of the class can be written down, and a potential written in its place would draw another soliton, as the spacetime diagrams already say.
What the drawing holds for every soliton of the class is the slice itself: the plane $y = 0$ of the path along $z$, at one moment, flat, a disc about the soliton's centre with the path marked as a line through it.
Without a potential the published metric has no length in it, so the disc's radius is in any length $\ell$, and nothing else is marked.
The caption carries the physics: on the flat slice the energy density the riding observers measure is $(K^2 - K_{ij}K^{ij})c^4/16\pi G = \sigma_2(\partial_i\partial_j\phi)\,c^4/8\pi G$, the sum of the principal minors of the Hessian of $\phi$, which Lentz arranged to be positive wherever his soliton has any, and which Jessica Santiago, Sebastian Schuster and Matt Visser showed an observer moving fast enough through the slices measures as negative.
The check is the published spatial metric, $\delta_{ij}$ with no dependence on $\phi$.

---

## Group 3. Flat slices, with a ring of free particles stretched on them

In these three the plane drawn is flat at every moment, and what the spacetime does is carry free particles on it apart along one direction and together along another, which no one surface can show and a ring of free particles does.
Each drawing is a sequence of moments of one plane, a flat disc at each, with a ring of particles that starts as a circle marked on it as a curve, twelve of the particles as points, and the disc's own coordinates the proper distances along the two axes of the chart.
At each moment the plane's metric is $g_{11}\,dx_1^2 + g_{22}\,dx_2^2$ with constant coefficients, so the chart point $(x_1, x_2)$ stands at $(\sqrt{g_{11}}\,x_1, \sqrt{g_{22}}\,x_2)$ in the drawing, and a ring of particles is an ellipse whenever the particles' chart positions are a linear image of a circle, which in all three they are.

### Kasner's universe

The published metric is $-c^2dt^2 + t^{2p_1}dx^2 + t^{2p_2}dy^2 + t^{2p_3}dz^2$, with $t$ in the fixed unit the powers are read in, and the spacetime diagrams declare $(p_1, p_2, p_3) = (-2/7, 3/7, 6/7)$, which sum to $1$ as their squares do.
Every $\Gamma^i{}_{tt}$ of the published metric vanishes, so a particle at rest in the chart stays there: the ring of particles at $x^2 + z^2 = \ell^2$ in the plane $y = 0$ keeps its chart positions while the plane stretches under it.
At time $t$ it is the ellipse of semi-axes $t^{p_1}\ell$ along $x$ and $t^{p_3}\ell$ along $z$, the most and least expanding directions of the three, and it encloses $\pi\ell^2t^{p_1 + p_3} = \pi\ell^2t^{4/7}$.

| $t$ | along $x$ | along $z$ |
|---|---|---|
| $1/4$ | $1.4860\,\ell$ | $0.3048\,\ell$ |
| $1/2$ | $1.2190\,\ell$ | $0.5520\,\ell$ |
| $1$ | $\ell$ | $\ell$ |
| $2$ | $0.8203\,\ell$ | $1.8114\,\ell$ |

The discs are drawn out to $2\ell$.
Toward the singularity at $t = 0$ the ring grows without bound along $x$ while it shrinks to nothing along $z$ and along $y$, so every sphere of particles is drawn out into a needle, the cigar singularity of the solution Edward Kasner found in 1921; the exponent $p_1 < 0$ that does it is what every Kasner universe but $(1, 0, 0)$ has, one direction contracting while the other two expand.
The checks are the published spatial metric at each moment, the particles' chart positions held fixed by $\Gamma^i{}_{tt} = 0$, and each ellipse against its semi-axes.

### Bianchi type I

The published metric is $-c^2dt^2 + a_1^2dx^2 + a_2^2dy^2 + a_3^2dz^2$, and the spacetime diagram declares dust: the scale factors solved from the published $G^x{}_x = G^y{}_y = G^z{}_z = 0$ by `null_rays.DustSolver`, starting from $a_i = 1$ with rates $(-0.5, 1.5, 2.0)\,\bar H$ at the reference instant, $\bar H$ their mean, with the time shifted so that the singularity is $t = 0$ and lengths in $c/\bar H$.
The reference instant is $c\bar Ht = 0.377980$.
Every $\Gamma^i{}_{tt}$ of the published metric vanishes, so the dust itself is at rest in the chart, and the ring drawn is a ring of the dust, $x^2 + z^2 = \ell^2$ in the plane $y = 0$, the ellipse of semi-axes $a_1\ell$ and $a_3\ell$:

| $c\bar Ht$ | $a_1$ | $a_2$ | $a_3$ |
|---|---|---|---|
| $0.1$ | $1.383749$ | $0.474583$ | $0.363183$ |
| $0.377980$ | $1$ | $1$ | $1$ |
| $1$ | $0.890637$ | $1.749871$ | $2.071726$ |
| $2$ | $0.916967$ | $2.641095$ | $3.440658$ |

The discs are drawn out to $3.5\,\ell$.
Near the singularity the dust universe behaves as Kasner's, drawn out along $x$ and flattened along the other two; later the dust's own gravity turns the contraction along $x$ round, $a_1$ reaching its least, $0.8881$, at $c\bar Ht = 1.1780$ and growing after, as every direction of a dust universe expands at late times.
The checks are those of Kasner, with the scale factors from the solver the spacetime diagram uses, which solves the published field equations for them.

### The pp-wave

The published exact plane wave is $\left(A(x^2 - y^2) + 2Bxy\right)c^2du^2 - 2c\,du\,dv + dx^2 + dy^2$, and its spacetime diagram declares $A = e^{-u^2}$ and $B = 0$, a pulse of the plus polarisation with $A$ in units of $1/L^2$ and $cu$ in units of $L$.
A surface of constant $u$ is a wave front, and its metric is $dx^2 + dy^2$ whatever $v$ is on it, since $g_{vv}$ and $g_{xv}$ vanish: the plane drawn is the flat wave front, on which the distance between two particles is $\sqrt{\Delta x^2 + \Delta y^2}$.

The published Christoffel symbols have no $\Gamma^u{}_{\mu\nu}$, so $u$ is an affine parameter of every geodesic, and $\Gamma^x{}_{uu} = -Ax - By$, $\Gamma^y{}_{uu} = -Bx + Ay$ give

$$\frac{d^2x}{d(cu)^2} = A\,x,\qquad \frac{d^2y}{d(cu)^2} = -A\,y.$$

The ring is twelve particles at rest on the circle $x^2 + y^2 = L^2$ before the wave arrives, and since the equations are linear, each particle at $(L\cos\alpha, L\sin\alpha)$ is at $(X(u)L\cos\alpha, Y(u)L\sin\alpha)$ later, with $X$ and $Y$ the solutions from $X = Y = 1$ at rest, integrated from $cu = -8$:

| $cu$ | $X$ | $Y$ | enclosed area over $\pi L^2$ |
|---|---|---|---|
| $-3$ | $1.0000030$ | $0.9999970$ | $1.0000$ |
| $-0.5$ | $1.184116$ | $0.830017$ | $0.9828$ |
| $0$ | $1.556176$ | $0.551250$ | $0.8578$ |
| $0.660753$ | $2.645564$ | $0$ | $0$ |

The wave stretches the ring along $x$ and squeezes it along $y$, and the area it encloses falls although the profile is harmonic and the spacetime a vacuum, $G_{\mu\nu} = 0$: the Weyl curvature shears the ring, and the shear alone focuses it.
At $cu = 0.660753\,L$ every particle of the ring reaches the $x$ axis at once and the ring is a segment of half length $2.6456\,L$; the particles are still moving across it, and past it the ring turns inside out.
That focusing of every ray of a plane wave is what Roger Penrose used in 1965 to show that no plane wave spacetime is globally hyperbolic.

The discs are drawn out to $3L$, at $cu = -3$, $-0.5$, $0$ and $0.660753$.
The checks are the published metric on the wave front, the particles' equations against the published Christoffel symbols, and each ring against its closed form $X\cos\alpha$, $Y\sin\alpha$.

---

## What the drawings need that the other eighteen did not

Four things are new, and each is small.

A curve on a surface: the ellipses, crescents and loops on the flat slices are not circles about the axis, so each is written under the surface's `curves` as a polyline of points `[X, Y, Z]` in the surface's own frame, with a class, and is checked to lie on its piece.
A point on a surface: the twelve particles of each ring, under `dots`.
Both turn with the surface when a reader turns the drawing, and `MFS/assets/embedding-turn.js` draws them again at every camera by the rules the generator draws them by.
A surface in Minkowski space: anti-de Sitter's view carries `"space": "minkowski"`, and every length along it is measured with $dX^2 + dY^2 - dZ^2$, with the checks and the rounding changed to match.
A slice whose angle is not a chart coordinate: the Malament-Hogarth plane turns $x$ into $y$ about the removed event, and the Mixmaster great sphere moves $\psi$ with $\phi$, so the slice reads the published metric along the turn, as Vaidya's slice already reads it along $v = T + r$.
The Mixmaster great sphere's $g_{\phi\phi}$ and defect are also written in half angles before they are evaluated, since the forms sympy gives subtract numbers that agree at the poles and lose every digit there.

Every drawing keeps the rest of the construction: the profile from the published metric, the isometry checks along and across it, the closed forms where one exists, and the declared functions checked against the published field equations where there are any.
