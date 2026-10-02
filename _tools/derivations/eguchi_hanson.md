# The Eguchi-Hanson space

Eguchi and Hanson's gravitational instanton is a Riemannian space of four dimensions: every coordinate is spatial, the signature is $(+,+,+,+)$, and nothing in it is a time.
It is the collection's one space with no time, so this note records what that changes, the three charts with their sources, and what the embedding diagram draws.

## Step 1: no time

No coordinate is declared a time in `DIMENSIONS`, so the checker scales nothing by $c$, and each chart's line element is the chart line element as it stands.
The curvature conventions need no change: the Riemann tensor, the contraction $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$ and the Weyl tensor are the same expressions in any signature.
The owned tags follow from the metric: `vacuum` and `four-dimensional` hold, and `stationary` and `static` do not, since no coordinate is timelike anywhere.
A dot in a geodesic equation is a derivative with respect to arc length, which is an affine parameter on every geodesic of a Riemannian space.
There is no light cone, so `null_rays.py` has nothing to draw, and no null infinity, so the space stands in `NOT_DRAWN` in `conformal.py`.

## Step 2: Eguchi and Hanson's chart

$ds^2 = (1 - a^4/r^4)^{-1}dr^2 + \tfrac{r^2}{4}(d\theta^2 + \sin^2\theta\,d\phi^2) + \tfrac{r^2}{4}(1 - a^4/r^4)(d\psi + \cos\theta\,d\phi)^2$.
It is Eguchi and Hanson's metric II, (9) with (17) of their letter of 1978 and (2.26) of their paper of 1979, $[1 - (a/r)^4]^{-1}dr^2 + r^2(\sigma_x^2 + \sigma_y^2) + r^2[1 - (a/r)^4]\sigma_z^2$, written out with $\sigma_z = \tfrac{1}{2}(d\psi + \cos\theta\,d\phi)$ and $\sigma_x^2 + \sigma_y^2 = \tfrac{1}{4}(d\theta^2 + \sin^2\theta\,d\phi^2)$, as Gibbons and Hawking's (3.20) of 1979 has it.
The one parameter $a$ is a length, the constant of integration of the self-duality equation.
The range of $\psi$ is $[0, 2\pi)$, half an Euler angle's: near $r = a$ at fixed $\theta$ and $\phi$, with $u^2 = r^2(1 - a^4/r^4)$, the metric is $\tfrac{1}{4}(du^2 + u^2d\psi^2)$, Eguchi and Hanson's (2.38), flat in polar coordinates exactly when $\psi$ has the period $2\pi$.
With that period every surface of constant $r > a$ is a 3-sphere with opposite points identified, and $r = a$ is a 2-sphere of radius $a/2$, the bolt.
`eguchi_hanson_pretty` keeps $r^4 - a^4$ whole and writes each sum in whichever of $\sin^2\theta$ and $\cos^2\theta$ leaves it fewer terms.

## Step 3: the Kähler chart

$\rho^4 = r^4 - a^4$ is Eguchi and Hanson's (2.31) of 1979, and their (2.32) is $ds^2 = [1 + (a/\rho)^4]^{-1/2}(d\rho^2 + \rho^2\sigma_z^2) + [1 + (a/\rho)^4]^{1/2}\rho^2(\sigma_x^2 + \sigma_y^2)$.
With $W = \sqrt{\rho^4 + a^4} = r^2$ that is $\tfrac{\rho^2}{W}\left(d\rho^2 + \tfrac{\rho^2}{4}(d\psi + \cos\theta\,d\phi)^2\right) + \tfrac{W}{4}(d\theta^2 + \sin^2\theta\,d\phi^2)$.
$\rho$ is the radius of the flat space of two complex dimensions the metric is Kähler on, $\rho^2 = |z_1|^2 + |z_2|^2$, and the bolt is $\rho = 0$.
The complex chart itself is not published: its components are those of the Kähler form (2.33) in four real coordinates, and the polar form in $\rho$ carries the same metric in a chart the printer writes in $\rho$, $a$ and $W$.
The chart is checked to be Step 2's carried along $r = (\rho^4 + a^4)^{1/4}$, slot by slot.

## Step 4: the two-centre chart

Gibbons and Hawking's multi-centre metrics are $V\,d\mathbf{x}\cdot d\mathbf{x} + V^{-1}(d\tau + \boldsymbol{\omega}\cdot d\mathbf{x})^2$ with $V = \epsilon + \sum_i 2m/|\mathbf{x} - \mathbf{x}_i|$ and $\nabla\times\boldsymbol{\omega} = \nabla V$; with $\epsilon = 0$ and two centres the metric is Eguchi and Hanson's, as Eguchi and Hanson thought certain in 1979, Gibbons and Hawking state, and Prasad showed by giving the map.
Prasad's letter was not to hand, so the constants are Ishihara, Kimura, Matsuno and Tomizawa's of 2006, who cite him: centres at $z = \pm a$ and $V = \tfrac{a}{8}(1/R_1 + 1/R_2)$.
They write the flat space in spherical coordinates; the chart here is the same form in cylindrical coordinates $(\rho, z, \psi)$, in which the two distances are $R_{1,2} = \sqrt{\rho^2 + (z \mp a)^2}$ and $\omega = \tfrac{a}{8}\left((z - a)/R_1 + (z + a)/R_2\right)$ multiplies $d\psi$.
The map from Step 2 is $z = r^2\cos\theta/a$, $\rho = \sqrt{r^4 - a^4}\,\sin\theta/a$, $\tau = a\phi/4$, with the azimuth $\psi$ Eguchi and Hanson's $\psi$; on it $R_1 = (r^2 - a^2\cos\theta)/a$ and $R_2 = (r^2 + a^2\cos\theta)/a$, so $r^2 = a(R_1 + R_2)/2$.
Working it by hand: $V(d\rho^2 + dz^2) = dr^2/(1 - a^4/r^4) + \tfrac{r^2}{4}d\theta^2$, since the flat metric in these prolate spheroidal coordinates is $(R_1R_2)(4r^2dr^2/(r^4 - a^4) + d\theta^2)$, and the terms in $d\psi$ and $d\tau$ give the rest with $\tau = a\phi/4$.
So $\tau$ has the period $\pi a/2$, the bolt is the segment $\rho = 0$, $|z| < a$ with the circle of $\tau$ over it, and the centres are its poles, where the circle of $\tau$ closes and the space is regular.
`eguchi_hanson_pullback` checks the map at six random points in forty digits, since sympy does not clear the two radicals symbolically.
`eguchi_hanson_distances` prints every value in $\rho$, $a$, $R_1$ and $R_2$: $z = (R_2^2 - R_1^2)/4a$ and $\rho^2 = \left((R_1 + R_2)^2 - 4a^2\right)\left(4a^2 - (R_1 - R_2)^2\right)/16a^2$, after which no relation is left among the generators and a value that vanishes is exactly zero.
The Kretschmann scalar is $24576a^2/(R_1 + R_2)^6$, which is $384a^8/r^{12}$.
The chart took 56 seconds to print on 2 October 2026, the other two a second each.

## Step 5: the curvature

All three charts are Ricci flat, so the Weyl tensor is the Riemann tensor.
With the orientation $dr\wedge d\theta\wedge d\phi\wedge d\psi$ the Riemann tensor is self dual, $R_{\mu\nu\rho\sigma} = \tfrac{1}{2}\sqrt{g}\,\epsilon_{\rho\sigma\alpha\beta}R_{\mu\nu}{}^{\alpha\beta}$ with $\sqrt{g} = r^3\sin\theta/8$, which `_tools/test_eguchi_hanson.py` checks in every slot; with the other orientation it is anti-self dual, which is why Eguchi and Hanson's signature $\tau = -1$ is Gibbons and Hawking's $+1$.
$K = 384a^8/r^{12}$, which Ghezelbash quotes (arXiv:2108.07210) and which is the sum of the squares of Eguchi and Hanson's curvature components (2.28): $384/a^4$ on the bolt, finite.
For a Ricci flat metric the Euler density is $K/32\pi^2$, and its integral over $r \ge a$ with $\psi$ of period $2\pi$ is $3/2$; the boundary's term is $1/2$, so $\chi = 2$, Eguchi and Hanson's (2.47), the Euler number of one bolt.

## Step 6: the embedding diagram

A slice of constant time has no meaning here, so the diagram draws three surfaces of the space itself, each the fixed set of an isometry and so totally geodesic.

The fibre is the surface of $r$ and $\psi$ over one point of the bolt, the fixed set of the rotation $\phi \to \phi + \alpha$, $\psi \to \psi - \alpha$ over the pole $\theta = 0$.
On it $ds^2 = dr^2/(1 - a^4/r^4) + \tfrac{r^2}{4}(1 - a^4/r^4)d\psi^2$, so the circle at $r$ has radius $\varrho = \tfrac{r}{2}\sqrt{1 - a^4/r^4}$, and with $s$ the proper distance $d\varrho/ds = (1 + a^4/r^4)/2 \le 1$: the surface embeds everywhere, at $dz/dr = \sqrt{3 + a^4/r^4}/2$.
On the bolt $d\varrho/ds = 1$, the pole of a smooth surface, and far away $d\varrho/ds \to 1/2$, a cone of half angle $30°$ whose circles are $\pi r$ long, half a plane's.
That half is the identification of opposite points at infinity made visible: with $\psi$ of period $4\pi$ the far circles would be a plane's and the surface would have a cone's tip on the bolt, $d\varrho/ds = 2$ there.

The bolt is the surface of $\theta$ and $\phi$ at $r = a$, the fixed set of the rotation in $\psi$: $\tfrac{a^2}{4}(d\theta^2 + \sin^2\theta\,d\phi^2)$, the round sphere of radius $a/2$.

The third surface is $\theta = \pi/2$ at $\psi = 0$ and $\psi = \pi$, the fixed set of $\theta \to \pi - \theta$, $\psi \to -\psi$.
In each fibre the two rays $\psi = 0$ and $\psi = \pi$ are one straight line through the bolt, so the two halves join smoothly along the bolt's equator.
On each half $ds^2 = dr^2/(1 - a^4/r^4) + \tfrac{r^2}{4}d\phi^2$, so $\varrho = r/2$ and $dz/dr = \sqrt{(3r^4 + a^4)/(4(r^4 - a^4))}$, vertical on the equator, the smallest circle, and tending to $\sqrt{3}/2$, the same cone.
Both halves run out to the one boundary; the surface has the shape of a wormhole and one asymptotic region.

Hanson and Sha's isometric embedding of the whole space is in flat space of eleven dimensions, so no single surface in three carries more than a two dimensional piece of it.
