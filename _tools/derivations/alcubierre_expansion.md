# The expansion of Alcubierre's Eulerian observers as a height

The captain's order of 28 September 2026 draws Alcubierre's warp drive as a height plot of the expansion $\theta$ of the observers who ride the slices of constant $t$, over the plane of the ship's path, the way Miguel Alcubierre drew it in 1994.
The height is $\theta$, a rate of change of volume, and says nothing about the shape of the slice, which is flat.

## The observers and their expansion

The published line element is

$$ds^2 = -c^2dt^2 + \left(dx - v_sf\,c\,dt\right)^2 + dy^2 + dz^2,$$

and in the chart $x^0 = ct$ its components are $g_{tt} = -1 + v_s^2f^2$, $g_{tx} = -v_sf$ and $g_{xx} = g_{yy} = g_{zz} = 1$.
Read as a $3 + 1$ split, the lapse is $1$, the shift is $\beta^x = -v_sf$ and the slices of constant $t$ are flat.

The Eulerian observers move along the unit normal to those slices.
From the published inverse metric, $g^{tt} = -1$ and $g^{tx} = -v_sf$, so the lapse is $1/\sqrt{-g^{tt}} = 1$ and

$$n_\mu = (-1, 0, 0, 0),\qquad n^\mu = -g^{\mu t} = (1,\ v_sf,\ 0,\ 0).$$

The determinant of the published metric is $-1$ for every $f$ and every $v_s$, so the divergence needs no measure:

$$\theta = \nabla_\mu n^\mu = \frac{1}{\sqrt{-g}}\,\partial_\mu\!\left(\sqrt{-g}\,n^\mu\right) = \partial_t(1) + \partial_x(v_sf) = v_s\,\partial_xf.$$

sympy computes exactly $v_s\,\partial_xf$ from the published components with $f$ an arbitrary function of $t$, $x$, $y$ and $z$ and $v_s$ an arbitrary function of $t$.
That is a rate per unit of $ct$, an inverse length.
The lapse is $1$, so the observers' proper time is $t$, and per unit of it the rate is

$$\theta = c\,v_s\,\partial_xf,$$

the expression in the metric file's convention.

The same number follows from the motion alone.
An Eulerian observer moves along $x$ at $dx/dt = c\,v_sf$, since $n^x/n^t = v_sf$ in the chart, and the slices are flat and Cartesian, so the fractional rate at which a small volume of them changes is the divergence of that velocity, $c\,v_s\,\partial_xf$.
On the path, a centred difference of $v_sf$ agrees with the published $\theta$ to $1.8 \times 10^{-10}$.

For Alcubierre's radial profile $f(r_s)$, with $r_s = \sqrt{(x - x_s)^2 + y^2 + z^2}$, the derivative of $r_s$ along $x$ is $(x - x_s)/r_s$, so

$$\theta = c\,v_s\,\frac{x - x_s}{r_s}\,\frac{df}{dr_s},$$

the captain's expression with $c$ kept.

## The sign, against Alcubierre's paper

Miguel Alcubierre wrote the drive in Classical and Quantum Gravity 11 (1994) L73 with $G = c = 1$ as $\alpha = 1$, $\beta^x = -v_s(t)\,f(r_s(t))$, $\beta^y = \beta^z = 0$ and $\gamma_{ij} = \delta_{ij}$, his equations (2) to (5), which is the published metric.
He defines the extrinsic curvature as $K_{ij} = \frac{1}{2\alpha}\left(D_i\beta_j + D_j\beta_i - \partial_tg_{ij}\right)$, his equation (9), and the expansion as $\theta = -\alpha\,\mathrm{Tr}\,K$, his equation (11).
With flat slices and $\alpha = 1$ the trace is $\partial_x\beta_x = -v_s\,\partial_xf$, so his $\theta$ is $v_s\,\partial_xf$, and sympy finds it equal to the divergence above for an arbitrary $f$.

His equation (12) prints $\theta = v_s\,(x_s/r_s)\,df/dr_s$, with $x_s$ in the numerator where the derivative of $r_s$ along $x$ has $x - x_s$.
Only with $x - x_s$ does $\theta$ match his Figure 1, which he describes in these words: "We clearly see how the volume elements are expanding behind the spaceship, and contracting in front of it."
With $x_s$ in the numerator, $\theta$ would carry the one sign of $-x_s$ all the way round the wall.
He drew the figure over $x$ and $\rho = \sqrt{y^2 + z^2}$ with $\sigma = 8$ and $R = v_s = 1$.

In the collection the bubble's centre follows $dx_s/dt = c\,v_s$, and the spacetime diagram declares $v_s = 2$ and Alcubierre's profile

$$f = \frac{\tanh\sigma(r_s + R) - \tanh\sigma(r_s - R)}{2\tanh\sigma R},\qquad r_s = \sqrt{(x - 2ct)^2 + y^2 + z^2},$$

with $R = 1$ and $\sigma = 4$, so the ship moves toward $+x$: the greatest value of the declared $f$ on the axis sits at $x = 0$, $0.2R$ and $1.0R$ at $ct = 0$, $0.1R$ and $0.5R$.
Through the wall $f$ falls from $1$ to $0$, so $df/dr_s < 0$ there, and $\theta < 0$ where $x > x_s$, ahead of the ship, and $\theta > 0$ behind it: space contracts ahead and expands behind, as in Alcubierre's figure.

At $t = 0$, where $x_s = 0$, on the plane $z = 0$, $\theta$ from the published metric with the declared $f$ is, in units of $c/R$:

| $x$ | $y$ | $\theta$ |
|---|---|---|
| $R$ | $0$ | $-4.002683$ |
| $-R$ | $0$ | $4.002683$ |
| $R/2$ | $0$ | $-0.282695$ |
| $-R/2$ | $0$ | $0.282695$ |
| $0$ | $R$ | $0$ |
| $R/\sqrt{2}$ | $R/\sqrt{2}$ | $-2.830324$ |
| $-0.8R$ | $0.6R$ | $3.202146$ |
| $1.5R$ | $0$ | $-0.282793$ |
| $2R$ | $0$ | $-0.005367$ |
| $3R$ | $0$ | $-1.8 \times 10^{-6}$ |

and the closed form above agrees with the published metric to $3.1 \times 10^{-15}$ at 4000 points of the square $|x|, |y| \le 3R$.
On the plane, $(x - x_s)/r_s = \cos\varphi$, with $\varphi$ the angle from the direction of travel, so $\theta = c\,v_s\cos\varphi\,f'(r_s)$: zero across the ship's own plane $x = x_s$, and largest on the path.

The greatest expansion is $4.0026828\,c/R$, behind the ship on the path at $r_s = 1.0000001R$, and the fastest contraction the same amount ahead of it.
For $\sigma R \gg 1$ the extreme is close to $v_s\sigma/2$ in units of $c/R$, since $f' \approx -(\sigma/2)\,\mathrm{sech}^2\sigma(r_s - R)$ in the wall.
Alcubierre's $\sigma = 8$ with $v_s = 1$ and the collection's $\sigma = 4$ with $v_s = 2$ have the same product, so both reliefs rise to about $4c/R$, and his wall is half as thick.
The circle $v_sf = 1$, where the published $g_{tt}$ vanishes, lies at $r_s = 1.0001676R$, and on the path there $\theta = \pm 4.0026810\,c/R$.
The wall where $f$ falls from $0.9$ to $0.1$ runs from $r_s = 0.7262R$ to $1.2747R$.

## The plane and the height

The plane is $z = 0$ at $t = 0$, through the path and the ship at its centre, taken on both sides of the path, so that $\rho = |y|$ there.
$\theta$ depends on $y$ and $z$ only through $\rho$, so the plane holds all of it.
It is drawn over the disc of radius $3R$ about the ship, as far as the spacetime diagram runs, where $|\theta|$ has fallen below $2 \times 10^{-6}\,c/R$.

The height is $z = \theta R^2/4c$, so a height of $R$ stands for an expansion of $4c/R$.
The expansion behind the ship then rises $1.0006707R$ above the plane and the contraction ahead sinks as far below it, over a disc $6R$ across.
At a quarter the relief is already steep, climbing at most $3.08$ in height per unit of distance, at $r_s = 0.835R$ on the path, because the wall is only about $R/2$ thick; a larger factor would stand the two lobes up as thin blades, and a smaller one would flatten the relief of the wall's edges into the plane.
The curves $\theta = \pm 2.0013\,c/R$, half the greatest expansion and contraction, are the crescents that crossed the path between $r_s = 0.7797R$ and $1.2203R$ and reached $60°$ either side of it; on the relief they are its level lines at heights $\pm 0.50034R$.

## The grid

The relief is a surface over the disc, $z(\rho, \phi)$, with $\phi$ measured from the direction of travel, sampled on a polar grid and drawn as flat triangles.
The bubble is round, and $z = A(\rho)\cos\phi$ with $A = v_s f'(\rho)R^2/4$, so the grid's circles follow the wall and a coarse even spacing of angles carries the $\cos\phi$; a Cartesian grid would need the fine spacing of the wall across the whole disc in both directions.

There are 72 angles, $\phi_j = 5j°$, so that every third one is a meridian at the $15°$ of the other figures.
The radii run from $0$ to $3R$ and hold every multiple of $R/4$ and the circle $v_sf = 1$ at $1.0001676R$, so each of those is a line of the grid.
Each cell between neighbouring radii $i$, $i + 1$ and neighbouring angles $j$, $j + 1$ is split into two triangles along the diagonal from $(i, j)$ to $(i + 1, j + 1)$, and a band between two radii is halved until every triangle in it lies within $5 \times 10^{-4}$ of the drawing's size of the surface $z = \theta R^2/4c$, that is within $3 \times 10^{-3}R$, measured as the distance in space from the triangle to the surface.
That distance is the vertical miss divided by $\sqrt{1 + |\nabla z|^2}$, and it is the one to hold, since on the wall a triangle that lies on the surface to a thousandth of $R$ misses it vertically by three times as much where the relief climbs at $3$.
The worst triangle misses by $2.94 \times 10^{-3}R$, and the grid has 43 radii, $72 \times 43 = 3096$ heights, which at six decimals take about 27 kilobytes.
The rim is the polygon of the 72 angles at $3R$, within $4.8 \times 10^{-4}$ of the size of the circle.

A tolerance of $5 \times 10^{-4}$ of the size is about a third of a pixel on a figure 628 wide.
Halving it to $2.5 \times 10^{-4}$ asks for 96 angles and 63 radii, 6048 heights, and doubling it to $10^{-3}$ for 29 radii at 72 angles.
