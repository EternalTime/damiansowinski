# The Randall-Sundrum braneworld

The four charts of `randall_sundrum.json` are written by `print_charts.py --metric randall_sundrum`, which computes every tensor from each line element with the checker's own `Geometry` and reads every printed value back before it writes.
This note records what each chart is and where it comes from, why the tensors come out as printed, and what the diagrams draw.

## Step 1. The proper distance chart

Randall and Sundrum's solution, equation (4) of "An Alternative to Compactification" and equation (12) of "A Large Mass Hierarchy from a Small Extra Dimension" with $y = r_c\phi$, is

$$ds^2 = e^{-2k|y|}\left(-c^2dt^2 + dx_1^2 + dx_2^2 + dx_3^2\right) + dy^2 ,$$

with $y$ the proper distance from the wall at $y = 0$ and $-\infty < y < \infty$.
Chamblin, Hawking and Reall call these horospherical coordinates of anti-de Sitter space, their equation (2.1) with $l = 1/k$.
Every component is taken in the chart $x^0 = ct$.

## Step 2. The conformally flat chart

Randall and Sundrum's second paper changes variable to $z \equiv \mathrm{sgn}(y)\left(e^{k|y|} - 1\right)/k$ to turn the equation of a graviton into a Schrödinger problem.
The same coordinate is Chamblin, Hawking and Reall's $w = z - z_0$ with the wall at $z_0 = l$, and the entry writes it $w$, the letter the domain wall's conformal chart already uses:

$$ds^2 = \frac{-c^2dt^2 + dx_1^2 + dx_2^2 + dx_3^2 + dw^2}{\left(1 + k|w|\right)^2}, \qquad 1 + k|w| = e^{k|y|} .$$

## Step 3. The Poincaré chart of one side

With $kz = e^{ky}$ on the side $y > 0$ the metric is $\left(-c^2dt^2 + d\mathbf{x}^2 + dz^2\right)/k^2z^2$, the Poincaré chart of anti-de Sitter space of radius $1/k$, Chamblin, Hawking and Reall's equation (3.1).
The wall stands at $z = 1/k$, the side kept is $z \ge 1/k$, and the conformal boundary $z = 0$ lies in the part cut away.
It has no kink and every tensor is anti-de Sitter's: $R_{\mu\nu} = -4k^2g_{\mu\nu}$, $R = -20k^2$ and $K = 40k^4$.

## Step 4. The chart between two walls

The first paper's equation (12), $ds^2 = e^{-2kr_c|\phi|}\eta_{\mu\nu}dx^\mu dx^\nu + r_c^2d\phi^2$, has $-\pi < \phi \le \pi$ with $\phi$ and $-\phi$ one point, the wall of positive tension at $\phi = 0$ and the wall of negative tension at $\phi = \pi$.
The published components are those of the line element as it stands, which holds for $|\phi| < \pi$; taken as a periodic function of $\phi$ the warp factor has a second kink at $\phi = \pi$, whose delta has the opposite sign, as the first paper's equation (10) writes it, and the chart's convention says so.

## Step 5. The kink and the wall's tension

The checker's `norm` writes $|y|$ as $y\,\mathrm{sgn}(y)$, so that $\partial_y|y| = \mathrm{sgn}(y)$ and $\partial_y\mathrm{sgn}(y) = 2\delta(y)$, as `domain_wall.md` Step 4 sets out.
A curvature component is then anti-de Sitter's value plus $\delta(y)$ times a function continuous at the wall, and the chart's `pretty` prints the two apart, the delta with its coefficient taken on the wall, where $e^{-2k|y|} = 1$:

$$R^t{}_{yty} = -\left(k^2 - 2k\,\delta(y)\right), \qquad R = -\left(20k^2 - 16k\,\delta(y)\right), \qquad G^t{}_t = 6k^2 - 6k\,\delta(y), \qquad G^y{}_y = 6k^2 .$$

With $G_{\mu\nu} = -\Lambda g_{\mu\nu} + \left(8\pi G_5/c^4\right)T_{\mu\nu}$ the first term is $\Lambda = -6k^2$ and the delta is a wall with $T^\mu{}_\nu = -\sigma\,\delta(y)\,\delta^\mu{}_\nu$ along its own directions, of tension

$$\sigma = \frac{6kc^4}{8\pi G_5} = \frac{3kc^4}{4\pi G_5},$$

which is Chamblin, Hawking and Reall's equation (2.2), $\sigma = 6/\kappa^2l$ with $\kappa^2 = 8\pi G_5$.
The Kretschmann scalar holds the square of the delta, which is no distribution, so it is stated off the wall, $K = 40k^4$.
`randall_sundrum_check` holds each chart to $R_{\mu\nu} = -4k^2g_{\mu\nu}$ on each side of the wall, and each chart after the first to being the first pulled back through its map on each side.

## Step 6. Light rays

On the plane of $t$ and $y$ at fixed $x_i$ the metric is $-e^{-2k|y|}c^2dt^2 + dy^2$, so a ray keeps $kct \mp \mathrm{sgn}(y)\left(e^{k|y|} - 1\right)$, which is $k(ct \mp w)$.
The only Christoffel symbols with an upper $x_i$ are $\Gamma^{x_i}{}_{x_iy}$, so no ray is turned out of the plane.
$\partial_t$ is a Killing vector, so along a ray $-g_{tt}\,d(ct)/d\lambda$ is constant and $d\lambda \propto -g_{tt}\,dw = dw/(1 + k|w|)^2$: the affine parameter from the wall to the horizon $w = \infty$ is $1/k$, finite, while $t$ is infinite.
The horizon is therefore an edge of the chart and no infinity of the spacetime, and Chamblin, Hawking and Reall's section 2 gives the continuations beyond it.

## Step 7. The diagrams

The spacetime diagrams draw the plane of $t$ and the fifth coordinate in each chart at $k = 1$, the two walls at $kr_c = 1/2$, and `--verify` checks the rays against the closed forms of Step 6.

The conformal diagram is the whole diamond by $p, q = \arctan(k(ct \mp w))$, each point a flat space of three dimensions, with the wall on the axis, the horizons on the four edges, the Poincaré chart on the right half, and the two walls the strip $|w| < \left(e^{kr_c\pi} - 1\right)/k$ between the two images of the second wall.
It is checked against each chart's published metric on each side of the wall, and the affine parameter of Step 6 is integrated from the published $g_{tt}$.

The embedding diagram is the surface of $y$ and $x_1$ at one moment, $dy^2 + e^{-2k|y|}dx_1^2$, a piece of the hyperbolic plane of curvature $-k^2$ on each side.
A strip of width $2\pi/k$ along $x_1$, rolled up, is a surface of revolution with $\rho = e^{-k|y|}/k$ and $dZ/dy = \sqrt{1 - e^{-2k|y|}}$, the tractrix

$$kZ = \mathrm{sgn}(y)\left(\mathrm{arcosh}\,e^{k|y|} - \sqrt{1 - e^{-2k|y|}}\right),$$

so the two sides are the two halves of Beltrami's pseudosphere and the wall is its cuspidal rim.
The area between the wall and either horizon is $\int_0^\infty e^{-ky}dy \cdot 2\pi/k = 2\pi/k^2$.
Between two walls the same surface is cut at $\phi = \pm\pi$, and both copies are drawn, as the other two diagrams draw them.
