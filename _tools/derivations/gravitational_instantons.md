# The gravitational instantons of 1977 and 1978

Hawking, Gibbons, Page and Pope's instantons are Riemannian spaces of four dimensions: every coordinate is spatial, the signature is $(+,+,+,+)$, and nothing in any of them is a time.
The second mate ruled on 2 October 2026 that they are one page, centred on the Euclidean Schwarzschild solution, with the others as charts.
So a chart here is a space of its own, as the products of `plebanski_hacyan` are, and `REGIONS` in `metric_tags.py` says which charts are the same space.
This note records the source of each chart, what was checked, and what the embedding diagram draws.
`eguchi_hanson.md` Step 1 says what having no time changes for the checker, and all of it holds here.

## Step 1: the sources

Six papers carry the family, and only two of them could be read in full on 2 October 2026.
Hawking's letter of 1977, Gibbons and Hawking's paper of 1977 and letter of 1978, and Page's two letters of 1978 are behind their publishers' walls, with no open copy that INSPIRE, OSTI, CERN or the authors' pages hold; their abstracts are the publishers' own.
Gibbons and Pope's paper on $\mathbb{CP}^2$ (Commun. Math. Phys. 61, 239) and Gibbons and Hawking's classification of 1979 (Commun. Math. Phys. 66, 291) are open at Project Euclid, and the second writes out every metric of the family with its periods, its fixed points and a citation of where it came from.
Page's own account of the four instantons with $\Lambda > 0$, the talk he gave in Moscow in December 1978, is arXiv:0912.4922.
Eguchi, Gilkey and Hanson's review of 1980 lists all of them again in its Appendix D, from a copy on Hanson's own page.
Every line element below is therefore taken from the classification of 1979, by equation number, and checked against the review's list, and no wording is attributed to a paper that was not read.

## Step 2: the Euclidean Schwarzschild solution

$ds^2 = (1 - r_s/r)\,d\tau^2 + (1 - r_s/r)^{-1}dr^2 + r^2(d\theta^2 + \sin^2\theta\,d\phi^2)$ is Gibbons and Hawking's (3.1) with $r_s = 2M$, and item 6 of Eguchi, Gilkey and Hanson's list.
The page keeps $r_s$, as the Schwarzschild page does, so $\tau$ is a length, $t = -i\tau/c$, and its period is $4\pi r_s$, their $8\pi M$.
With that period $r = r_s$ is a 2-sphere of area $4\pi r_s^2$, their $16\pi M^2$, the bolt of $\partial_\tau$, and the manifold is $\mathbb{R}^2 \times S^2$.
The chart is checked to be the published metric of `schwarzschild` with $g_{tt}$ reversed, which is what $x^0 = i\tau$ does to a static metric, and to be Ricci flat.
$K = 12r_s^2/r^6$, as Schwarzschild's.

## Step 3: the chart regular on the bolt

$x = 2r_s\sqrt{1 - r_s/r}$ turns the first two terms into $(x/2r_s)^2d\tau^2 + (r/r_s)^4dx^2$: with $\tau/2r_s$ an angle of period $2\pi$ that is a plane in polar coordinates about $x = 0$, times a factor that is $1$ there.
So the bolt is the origin of polar coordinates and nothing more, and $x$ runs to $2r_s$ at infinity.
The areal radius is the name $r = 4r_s^3/(4r_s^2 - x^2)$, and `instanton_pretty` writes every value in $x$, $r_s$ and $r$.
The construction is the one Eguchi, Gilkey and Hanson use for the bolts of the Eguchi-Hanson space and of $\mathbb{CP}^2$ on their page 365, a radial coordinate in which the fibre is $du^2 + u^2d\psi^2$ to leading order; here it is exact in the circle's term.
A coordinate of this form is remembered from Gibbons and Hawking's paper of 1977, as $x = 4M(1 - 2M/r)^{1/2}$, but that paper could not be read, so nothing on the page attributes the chart to it.
The chart is checked to be Step 2's carried along the map, slot by slot.

## Step 4: the two Taub-NUT instantons

Taub-NUT's published metric with $t = -i\tau/c$ and $l = in$ is $f\,(d\tau + 2n\cos\theta\,d\phi)^2 + dr^2/f + (r^2 - n^2)(d\theta^2 + \sin^2\theta\,d\phi^2)$ with $f = (r^2 - 2mr + n^2)/(r^2 - n^2)$, real and positive definite for $r > n$.
Hawking's self-dual solution is $m = n$, $f = (r - n)/(r + n)$, Gibbons and Hawking's (3.9) and item 2 of the review's list in the radius and the one-forms $\sigma$ it uses there.
Page's is $m = 5n/4$, $f = (r - 2n)(r - n/2)/(r^2 - n^2)$, their (3.16), where the coefficient printed as $(r - 2r)$ is $(r - 2n)$, as item 11 of the review's list confirms: $r^2 - 2.5Nr + N^2$.
In both $\tau$ has the period $8\pi n$ and $(\tau/2n, \theta, \phi)$ are Euler angles on the 3-spheres of constant $r$.
The self-dual solution closes on the point $r = n$, a nut, and the manifold is $\mathbb{R}^4$; Page's closes on the 2-sphere $r = 2n$ of area $12\pi n^2$, a bolt, with the surface gravity $1/4n$.
Each chart is checked to be the published metric of `taub_nut` continued, at its own mass, and to be Ricci flat.
`_tools/test_gravitational_instantons.py` holds the self-dual chart's Riemann tensor to being self dual up to orientation and Page's to being neither, the nut to being a regular point and the bolt to its area and surface gravity.

## Step 5: the multi-centre metrics

Gibbons and Hawking's letter of 1978 is quoted by the review, item 5 of its list, as $V^{-1}(d\tau + \boldsymbol\omega\cdot d\mathbf{x})^2 + V\,d\mathbf{x}\cdot d\mathbf{x}$ with $\nabla V = \pm\nabla\times\boldsymbol\omega$ and $V = \epsilon + 2m\sum_i 1/|\mathbf{x} - \mathbf{x}_i|$.
The classification of 1979 writes the same metric with $V$ for the reciprocal, (3.17) to (3.19), and the nut parameter $n$ for $m$; the page takes the review's $V$ and the classification's $n$, so that the self-dual Taub-NUT chart is the case of one centre in the same letter.
With centres anywhere $\boldsymbol\omega$ has three components, and the Cartesian chart with $V$ and all three free printed 1.9 megabytes on 2 October 2026, a megabyte of it the Riemann and Weyl tensors, twice the largest chart of the collection.
The chart is therefore the family with its centres on one axis: cylindrical coordinates $(\rho, z, \phi)$, $V(\rho, z)$ and the one component $\omega = \omega_\phi(\rho, z)$, free, with no component assuming the field equation, as the cylindrical chart of `israel_wilson_perjes` is.
The field equation is then $\partial_\rho\omega = \rho\,\partial_zV$ and $\partial_z\omega = -\rho\,\partial_\rho V$, whose integrability is Laplace's equation.
`instanton_on_shell` writes every derivative of $\omega$ and every second derivative of $V$ along $z$ by them, and the chart is checked to be Ricci flat on every solution, to be the self-dual Taub-NUT chart at $V = 1 + 2n/R$ along $r = R + n$, and to be the published two-centre chart of `eguchi_hanson` at $\epsilon = 0$ with two centres.
The `vacuum` tag cannot be read from a chart with free functions, so `REGIONS` leaves it out.

## Step 6: the complex projective plane

$ds^2 = (1 + \Lambda r^2/6)^{-2}dr^2 + \tfrac{r^2}{4}(1 + \Lambda r^2/6)^{-2}(d\psi + \cos\theta\,d\phi)^2 + \tfrac{r^2}{4}(1 + \Lambda r^2/6)^{-1}(d\theta^2 + \sin^2\theta\,d\phi^2)$ is Gibbons and Pope's (25), Gibbons and Hawking's (3.24) and item 3 of the review's list.
$\psi$ has the period $4\pi$, Gibbons and Pope's (24); $r = 0$ is a nut and $r = \infty$ a bolt of area $6\pi/\Lambda$.
It solves $R_{\mu\nu} = \Lambda g_{\mu\nu}$, which is checked, and $K = 16\Lambda^2/3$, a constant, as it must be on a homogeneous space.
The chart in $r$ reaches the bolt only at infinity, so a second chart takes the distance from the nut: Gibbons and Pope's Table 1 writes the metric as $d\chi^2 + A(d\psi + \cos\theta\,d\phi)^2 + B(d\theta^2 + \sin^2\theta\,d\phi^2)$ with $A = \tfrac{1}{4}(6/\Lambda)\sin^2(\chi\sqrt{\Lambda/6})\cos^2(\chi\sqrt{\Lambda/6})$, $B = \tfrac{1}{4}(6/\Lambda)\sin^2(\chi\sqrt{\Lambda/6})$ and $\chi$ from $0$ to $\tfrac{\pi}{2}\sqrt{6/\Lambda}$.
That chart takes the length $L = \sqrt{6/\Lambda}$ for its one parameter, so that the argument of every sine is $\chi/L$, and is checked to be the chart in $r$ carried along $r = L\tan(\chi/L)$.
Page's (10) is the same chart in the angle $2\chi/L$.
The test file holds the Weyl tensor to being half flat and the volume to $18\pi^2/\Lambda^2$.

## Step 7: Page's space

Gibbons and Hawking's (3.25) is $ds^2 = \tfrac{3}{\Lambda}(1 + \nu^2)\{Q\,d\rho^2/P + Q\,(d\theta^2 + \sin^2\theta\,d\phi^2)/(3 + 6\nu^2 - \nu^4) + P\sin^2\rho\,(d\psi + \cos\theta\,d\phi)^2/(4(3 + \nu^2)^2Q)\}$ with $Q = 1 - \nu^2\cos^2\rho$, $P = 3 - \nu^2 - \nu^2(1 + \nu^2)\cos^2\rho$, and $\nu$ the positive root of $\nu^4 + 4\nu^3 - 6\nu^2 + 12\nu - 3 = 0$, $0.2817$.
Item 10 of the review's list is the same metric in $x = \cos\rho$, and Page's own (12) and (13) another form of it in the angle he calls $\chi$, the letter the chart takes.
Written that way the metric is Einstein at that one value of $\nu$ and at no other: sympy gives $R_{\mu\nu} - \Lambda g_{\mu\nu}$ proportional to the quartic.
A checker that compares rational functions would have to reduce every comparison by the quartic, and `PARAMETER_RELATIONS` has no rational parametrisation of a single algebraic number.
The Einstein equation itself says how to avoid that.
With the coefficient of $d\Omega^2$ written $Q/N$ and that of $(d\psi + \cos\theta\,d\phi)^2$ written $B\,P\sin^2\rho/4Q$, $R_{\mu\nu} = \Lambda g_{\mu\nu}$ holds for every $\nu$ exactly when $N = 3 + 6\nu^2 - \nu^4$ and $B = 16\nu^2/N^2$.
That is the local metric, a limit of Kerr-de Sitter's, with $\nu$ free; Gibbons and Hawking's $B = 1/(3 + \nu^2)^2$ equals it exactly where $4\nu(3 + \nu^2) = N$, which is the quartic.
What the quartic expresses is regularity: with $\psi$ of period $4\pi$ the circle of $\psi$ closes on each bolt at the rate $d\varrho/ds = 4\nu(3 + \nu^2)/N$, which is $1$, the rate of a plane, at the root and nowhere else.
So the chart publishes the local form, every component an identity in $\nu$, and the test file holds the rate to the quartic and the published metric to Gibbons and Hawking's at the root in forty digits.
The chart names $P$, $Q$ and $N$ and `instanton_pretty` writes every value in them.
The bolts $\chi = 0$ and $\chi = \pi$ are 2-spheres of area $12\pi(1 - \nu^4)/(N\Lambda)$, their (3.27), which the test holds too.

## Step 8: the tags

`REGIONS` counts five spaces: Schwarzschild's two charts, the self-dual Taub-NUT chart, Taub-bolt, the two charts of $\mathbb{CP}^2$, and Page's.
Three are Ricci flat and two are Einstein spaces with $\Lambda > 0$, so the page carries neither `vacuum` nor `Einstein space`, since a tag cannot say which chart it is true of; the prose and each chart's convention say it.
No chart is stationary, since no coordinate is timelike anywhere.

## Step 9: the embedding diagram

No chart has a moment of time, so each view is a surface of the space itself, the fixed set of an isometry and so totally geodesic, as the Eguchi-Hanson space's three are, and each belongs to its chart.
Where the angle of the surface is the fibre coordinate, the slice sweeps it at the rate that takes it once round its period while $\phi$ runs to $2\pi$.
At $\theta = 0$ the azimuth moves nothing by itself, and the form $d\tau + 2n\cos\theta\,d\phi$ there is $d(\tau + 2n\phi)$, so $\tau = 2n\phi$ takes $\tau + 2n\phi$ once round $8\pi n$.

The cigar is the surface of $r$ and $\tau$ of the Euclidean Schwarzschild solution: $\varrho = 2r_s\sqrt{1 - r_s/r}$, and with $s$ the proper distance $d\varrho/ds = r_s^2/r^2 \le 1$, so it embeds everywhere, at $dz/dr = \sqrt{(r + r_s)(r^2 + r_s^2)/r^3}$.
On the bolt the rate is $1$ and far away the surface is a cylinder of radius $2r_s$, whose circumference $4\pi r_s$ is the period of $\tau$.
In the regular chart the same surface has $\varrho = x$.
The surface $\theta = \pi/2$ at $\tau = 0$ and $\tau = 2\pi r_s$ is the fixed set of $\tau \to -\tau$: its two halves are one straight line through the bolt in each fibre, and it is Flamm's paraboloid on both sides of the bolt's equator, the moment of time symmetry of Kruskal's manifold.

The nut's surface has $\varrho = 4n\sqrt{(r - n)/(r + n)}$ and $d\varrho/ds = 4n^2/(r + n)^2$, and Taub-bolt's $\varrho = 4n\sqrt{f}$ and $d\varrho/ds = 2nf' = n(5r^2 - 8rn + 5n^2)/(r^2 - n^2)^2$, each $1$ where the circle closes and each a cylinder of radius $4n$ far away.

The multi-centre view is the axis $\rho = 0$ of two centres at $z = \pm 2n$ with $\epsilon = 1$: the circle of $\tau$ has the radius $4n/\sqrt{V}$, which closes at each centre, so the axis is a sphere between the centres, widest at $z = 0$ with radius $4n/\sqrt{3}$, and a cigar beyond each.
Between the centres $\omega = 0$ and $\tau = 4n\phi$; beyond them $\omega = \pm 4n$, the Dirac string, and $\phi$ alone runs once round the circle.
The sphere is the 2-cycle that two centres bound, which Gibbons and Hawking's signature counts.

$\mathbb{CP}^2$'s view is the surface of the radius and $\psi$ at $\theta = 0$ with $\psi = \phi$: $\varrho = (L/2)\sin(2\chi/L)$, a round sphere of radius $L/2$, a complex line, drawn from the distance chart for both charts, since the chart in $r$ reaches its far pole only at infinity.
Page's is the same surface of his space at $\Lambda = 1$ and $\nu$ the root to thirty digits: the fibre of the bundle, a closed surface from one bolt to the other.
