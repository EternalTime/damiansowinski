# How far the Krasnikov tube tips the light cone, as a height

The captain's order of 28 September 2026 draws the Krasnikov tube's figure with the height standing for a physical quantity over the plane of the tube's axis, since the tube's cross section alone is a flat disc.
The quantity is $1 - k$, how far the tube tips the back edge of the light cone from where it stands in flat space.

## What $1 - k$ measures

The published line element, in Allen Everett and Thomas Roman's four dimensional form of Serguei Krasnikov's tube, is

$$ds^2 = -\left(c\,dt - dx\right)\left(c\,dt + k\,dx\right) + dr^2 + r^2d\phi^2,$$

and in the chart $x^0 = ct$ its components are $g_{tt} = -1$, $g_{tx} = (1 - k)/2$, $g_{xx} = k$, $g_{rr} = 1$ and $g_{\phi\phi} = r^2$.

The two null directions of the plane of $t$ and $x$ are

$$\ell_{\text{out}}^\mu = (1, 1, 0, 0),\qquad \ell_{\text{back}}^\mu = (k, -1, 0, 0),$$

and sympy finds both null with the published components for an arbitrary $k$.
The outbound edge of the light cone never moves.
Along the back edge $d(ct) = -k\,dx$, so a light signal sent back toward $x = 0$ over a length $L$ arrives after $\Delta(ct) = kL$, where in flat space it would take $L$.
So $1 - k$ is the coordinate time, as $ct$ per unit of length, that a signal sent home gains on flat space:

- $1 - k = 0$ outside the tube, where the metric is Minkowski's;
- $1 - k = 1$ where $k = 0$: the back edge lies along the axis, a signal home arrives at the $ct$ it left, and the published $g_{xx}$ vanishes;
- $1 - k > 1$ where $k < 0$: the back edge dips below the axis, a signal home arrives at an earlier $ct$ than it left, and the direction along the tube is a time;
- $1 - k = 2$ would be $k = -1$, where the determinant $-\frac{1}{4}r^2(1 + k)^2$ vanishes and the metric degenerates, and Krasnikov's $k = -1 + \delta$ keeps the tube short of it by $\delta$.

In the plane of $x$ and $ct$ the back edge of the forward light cone climbs $k$ in $ct$ for each unit it runs back along $x$, against $1$ in flat space, so $1 - k$ is how far it has been tipped down, as a change of slope.
A round trip over a length $D$, out at the speed of light and home along the back edge, takes $(1 + k)D = \left(2 - (1 - k)\right)D$ of $ct$, so $1 - k$ is also the part of the flat round trip, $2D$, that the tube takes away, per unit of length.
It is twice the published $g_{tx}$, exactly, which sympy confirms.

## The declared tube and the moment drawn

The spacetime diagram declares a tube along $x$ from $0$ to $D = 4\rho_0$, built by a ship that left $x = 0$ at $t = 0$ at the speed of light,

$$1 - k = (2 - \delta)\,S\!\left(\frac{\rho_0^2 - r^2}{2\rho_0}\right)S(ct - x)\,S(x)\,S(D - x),$$

with $\delta = 0.2$, $\rho_0 = 1$ and $S(q) = \left(1 + \tanh(q/0.15)\right)/2$, which agrees with the declared $k$ to $1.3 \times 10^{-15}$ at 5000 points of the plane drawn.
The factor $S(ct - x)$ holds the tube inside the future light cone of the departure, so nothing is built ahead of the ship.

The moment drawn is $ct = 5\rho_0$, the last moment of the spacetime diagram, which runs from $ct = -\rho_0$ to $5\rho_0$.
The ship reached $x = D$ at $ct = 4\rho_0$, and a unit of $ct$ later $S(ct - x) \ge S(\rho_0) = 1 - 1.6 \times 10^{-6}$ all along the tube, so the whole tube stands, from $x = 0$ to $x = D$, and the path home is open.
At an earlier moment the ridge would end at the ship: at $ct = 3\rho_0$, $1 - k$ is $1.797709$ on the axis at $x = 2\rho_0$ but $0.002282$ at $x = 3.5\rho_0$, ahead of the ship.

## The plane and the height

The plane is the plane of the axis, $\phi = 0$ and $\phi = \pi$, with $y = \pm r$ the distance from the axis on either side, at the moment $ct = 5\rho_0$.
It runs over $-\rho_0 \le x \le 5\rho_0$, as the spacetime diagram does, and out to $r = 2\rho_0$ on each side, as the cross section of the tube was drawn, where $1 - k$ has fallen to $3.7 \times 10^{-9}$.
Inside the tube, where $k < 0$, a surface of constant $t$ is no moment of space, so the height is a plot of $1 - k$ over the coordinates $x$ and $r$ and no embedding of any slice.

The height is $z = (1 - k)\rho_0$.
Along the axis in the middle of the tube it is $1.7977122\rho_0$, $k = -0.7977122$, the greatest anywhere, and the tube stands as a flat topped ridge $1.8\rho_0$ high, $2\rho_0$ wide and $D$ long.
At chosen points:

| $x$ | $r$ | $1 - k$ |
|---|---|---|
| $2\rho_0$ | $0$ | $1.7977122$ |
| $2\rho_0$ | $\rho_0/2$ | $1.7879529$ |
| $2\rho_0$ | $\rho_0$ | $0.9000000$ |
| $2\rho_0$ | $1.3\rho_0$ | $0.0179132$ |
| $2\rho_0$ | $1.5\rho_0$ | $0.0004326$ |
| $0$ | $0$ | $0.8988561$ |
| $0.2\rho_0$ | $0$ | $1.6809163$ |
| $3.8\rho_0$ | $0$ | $1.6809161$ |
| $4\rho_0$ | $0$ | $0.8988546$ |
| $-\rho_0/2$ | $0$ | $0.0022849$ |
| $4.5\rho_0$ | $0$ | $0.0022820$ |

At $r = \rho_0$ the radial factor is exactly $1/2$, so $1 - k = 0.9$ there to the eleventh decimal.
The ends differ in the sixth decimal because $S(ct - x)$ at $ct = 5\rho_0$ is $1 - 1.6 \times 10^{-6}$ at $x = 4\rho_0$ and closer to $1$ at $x = 0$.

The level line $1 - k = 1$, at height $\rho_0$, is the curve $k = 0$ that the cross section marked as a circle: in the middle of the tube it lies at $r = 0.9831218\rho_0$, and on the axis it crosses at $x = 0.0169506\rho_0$ and $3.9830492\rho_0$.
Inside it the direction along the tube is a time.
The walls climb at most $6.03$ in height per unit of distance across the tube at $x = 2\rho_0$ and $5.99$ along the axis at its ends, since each step is only $0.15\rho_0$ wide.

## The grid

The ridge is straight along $x$ and its walls run along the axes of the chart, so the height is sampled on a Cartesian grid of the plane, $z(x, y)$, and drawn as flat triangles.
The values of $x$ and of $y$ each hold every multiple of $\rho_0/2$, and $y = 0$, the axis, so each of those is a line of the grid.
Each cell between neighbouring values $x_i$, $x_{i+1}$ and $y_j$, $y_{j+1}$ is split into two triangles along the diagonal from $(i, j)$ to $(i + 1, j + 1)$, and every interval of $x$ or of $y$ whose cells hold a triangle further than $5 \times 10^{-4}$ of the drawing's size, $3 \times 10^{-3}\rho_0$ on the long side of $6\rho_0$, from the surface $z = (1 - k)\rho_0$, measured as the distance in space, is halved until none does.
The grid has 155 values of $x$ and 91 of $y$, $14105$ heights, about 122 kilobytes at six decimals, and the worst triangle misses by $2.52 \times 10^{-3}\rho_0$.
Sampled at 300000 points, the worst vertical miss is $1.28 \times 10^{-3}$ of the size, on the steepest wall.
A polar grid about the middle of the tube would need its finest spacing of angles everywhere the straight walls cross its circles.
