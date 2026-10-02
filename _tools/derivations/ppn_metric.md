# The parametrised post-Newtonian metric

Why each chart is the one published, what a post-Newtonian order is to the checker, and what `print_charts.py`, `null_rays.py`, `conformal.py` and `embedding.py` check before they write.

## Step 1: the charts and their sources

The isotropic line element is equation (11) of Thorne and Will (1971), "the original Eddington (1922)-Robertson (1962)-Schiff (1967) version of the PPN formalism": $ds^2 = -(1 - 2m/r + 2\beta m^2/r^2)c^2dt^2 + (1 + 2\gamma m/r)(dx^2 + dy^2 + dz^2)$ with $r^2 = x^2 + y^2 + z^2$ and $m = GM/c^2$, written here in the signature $(-,+,+,+)$.
It is the metric of Will (2014), Box 2, for a spherical body at rest, where the only potential is $U = m/r$.
Thorne and Will print it in Cartesian coordinates, which is the Cartesian chart; the isotropic chart is the same in polar coordinates, and `ppn_metric_check` pulls one back onto the other.
The areal chart is the form Eddington (1923) expands on page 105, $-(1 + a_1/r + \dots)^{-1}dr^2 - r^2d\Omega^2 + (1 + b_1/r + b_2/r^2 + \dots)dt^2$, in today's letters: with the areal radius $r + \gamma m$ his $a_1$ is $-2\gamma m$ and his $b_2$ is $2(\beta - \gamma)m^2$, and his three points, $b_1 = -2m$, $a_1 = -2m$ and $b_2 = 0$, are $\beta = \gamma = 1$.
`ppn_metric_check` carries the isotropic chart along $r \to r - \gamma m$ and holds the result to the areal chart, $g_{tt}$ through $m^2$ and the metric of space through $m$.
The rotating chart adds Will's (2014) equation (70), $g_{0i} = -\tfrac{1}{2}(4\gamma + 4 + \alpha_1)V_i$ with $\mathbf{V} = \tfrac{1}{2}\mathbf{J} \times \mathbf{x}/r^3$, for a body at rest in the preferred frame.
With $J = Mca$ that is $g_{t\phi} = -2\Delta am\sin^2\theta/r$, where $\Delta = \tfrac{1}{2}(1 + \gamma + \tfrac{1}{4}\alpha_1)$ is the factor of Will's (71) for the precession of a gyroscope, $1$ in general relativity.
For a stationary body the potentials $V_i$ and $W_i$ of the standard gauge are equal, so the term is the same in that gauge.
`ppn_metric_check` holds the chart to the isotropic one at $a = 0$ and to the published Lense-Thirring chart of `hartle_thorne` at $\beta = \gamma = 1$, $\alpha_1 = 0$ and first order in $m$.

## Step 2: what a post-Newtonian order is to the checker

The post-Newtonian metric is known to different orders in different components: with $v/c$ counted as first order and $m/r$ as second, $g_{tt}$ through the fourth, $g_{ti}$ through the third and $g_{ij}$ through the second (Will 2014, section 3.2).
So no single order cuts its tensors rightly: cut at $m$, the connection loses the term in $\beta$ that moves Mercury's perihelion, and cut at $m^2$, it prints terms of $\Gamma^i{}_{jk}$ that the metric's missing $m^2/r^2$ in $g_{ij}$ would change.
`ORDERS` names the time coordinate as a third item, `({"m": 2, "a": 1}, 2, "t")`, and `Reader.kept_order` then gives the order each component is kept to:

- a component of the connection, the Riemann tensor or the Ricci tensor with $n$ time indices, wherever they stand: $2 + n$, so $m^2$ with two, $am$ with one, and $m$ with none;
- a component of the Einstein or the Weyl tensor: $2$, or $3$ where $n$ is odd, since each takes a trace with $g_{tt}$, which multiplies the unknown fourth order of $R_{ij}$;
- the Ricci scalar: $2$;
- the Kretschmann scalar: $4$, the square of the lowest order of the Riemann tensor.

`Geometry` builds every tensor to the fourth order whole, and `PostNewtonian` wraps it and cuts each component as it hands the tensor over, so no product is cut before it is formed.
`geometry_of` returns the wrapped geometry to the checker, to `metric_tags.py` and, through `chart_printer.Chart`, to the printer.
The metric and its inverse are published whole, as the line element writes them, so the drawings read a metric and its exact inverse.

The rule is sound, and `_tools/test_ppn_metric.py` holds it so: with an arbitrary function of $r$ and $\theta$ added to each of the ten metric components at the first order the post-Newtonian metric leaves open, sixth in $g_{tt}$, fifth in $g_{ti}$ and fourth in $g_{ij}$, every tensor handed over is unchanged, while $G_{tt}$ cut at the fourth order, as $R_{tt}$ is, changes.

## Step 3: what the components say

Each value is printed order by order by `ppn_orders`, as $\Gamma^r{}_{tt} = m/r^2 - 2(\beta + \gamma)m^2/r^3$.
The Ricci tensor is $R_{tt} = -(1 - 2\beta + \gamma)m^2/r^4$ and $R_{rr} = 2(1 - \gamma)m/r^3$, with $R_{\theta\theta} = -(1 - \gamma)m/r$, and vanishes only at $\beta = \gamma = 1$: `ppn_metric_check` holds each chart to vanishing there and to not vanishing elsewhere.
The Weyl tensor carries $1 + \gamma$, the factor of the bending of light, and the Kretschmann scalar is $24(1 + \gamma^2)m^2/r^6$, Schwarzschild's $48m^2/r^6$ at $\gamma = 1$.
At $\beta = \gamma = 1$ the isotropic chart is Schwarzschild's metric in isotropic coordinates, $-((1 - m/2r)/(1 + m/2r))^2c^2dt^2 + (1 + m/2r)^4(dr^2 + r^2d\Omega^2)$, to these orders, which the check holds too.

## Step 4: the drawings

Every drawing takes $m = 1$, $\beta = \gamma = 1$ and a body whose surface is at $R = 10\,m$ in the isotropic radius, which is $11\,m$ in the areal radius; the rotating chart takes $a = 2m$ and $\alpha_1 = 0$.
The spacetime diagrams draw the plane of $t$ and $r$ of the isotropic and areal charts, the line through the body's centre of the Cartesian chart with the body hatched, and the rotating chart on its axis and on its equator with $\phi$ divided out.
`quotient_checks` in `null_rays.py` holds the lifted rays of that last view to the geodesic equation by the published Christoffel symbols, which are cut at post-Newtonian order, so it makes the check with $v/c$ a hundredth of what is drawn, where the terms left out are below its tolerance and the dragging is a hundred times above it.
`ppn_star` is the tortoise coordinate of each plane by quadrature, zero at the surface, which `--verify` holds the rays to and the conformal diagram builds its triangle from; to first order it is $r - R + (1 + \gamma)m\ln(r/R)$, whose logarithm is Shapiro's delay.
The embedding diagram draws the equator at $t = 0$, $(1 + 2\gamma m/r)(dr^2 + r^2d\phi^2)$, whose circles have the radius $\sqrt{r(r + 2\gamma m)}$ and which rises as $dz/dr = \sqrt{\gamma m(2r + 3\gamma m)/(r(r + 2\gamma m))}$; it is a plane at $\gamma = 0$, the same surface for every spin, and $0.30\,m$ below Flamm's paraboloid over the same circles by $3R$.
