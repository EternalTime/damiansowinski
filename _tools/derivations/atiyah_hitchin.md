# The Atiyah-Hitchin space

Atiyah and Hitchin's space of two slowly moving magnetic monopoles is a Riemannian space of four dimensions with signature $(+,+,+,+)$ and no time in it, as Eguchi and Hanson's is.
Equations are named by their labels in the sources' TeX, which the arXiv copies carry.
This note records the sources of its two charts, the source errors met on the way, the chart that is left out and why, and what the drawings draw.
Every source named here is read in full text in the worker's literature review; Atiyah and Hitchin's letter of 1985, their Philosophical Transactions paper of 1985, their book of 1988 and Gibbons and Manton's paper of 1986 could not be read in full, and no equation below rests on them alone.

## Step 1: no time

No coordinate is declared a time in `DIMENSIONS`, so the checker scales nothing by $c$, and each chart's line element is the chart line element as it stands.
The Riemann tensor, the contraction $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$ and the Weyl tensor are the same expressions in any signature.
A dot in a geodesic equation is a derivative with respect to arc length.
There is no light cone and no null infinity, so `null_rays.py` records both charts in `NOT_DRAWN` and `conformal.py` names the space in its `NOT_DRAWN`.

## Step 2: the form of the metric

Every metric of four dimensions on which the rotations act with orbits of three dimensions is $ds^2 = f^2dr^2 + a^2\sigma_1^2 + b^2\sigma_2^2 + c^2\sigma_3^2$, as Tong writes it in the subsection "The Moduli Space" of the monopole lecture of his TASI lectures (the equation labelled `ahm` in his TeX), Schroers and Singer in the subsection "The Atiyah-Hitchin manifold and related spaces" (`AHmetric`) and Ghezelbash and Mann in their section on embedding the space in an M2-brane metric (`AHmetric`).
The forms are Hanany and Pioline's of section 3.1, $\sigma_1 = \sin\psi\,d\theta - \cos\psi\sin\theta\,d\phi$, $\sigma_2 = -\cos\psi\,d\theta - \sin\psi\sin\theta\,d\phi$, $\sigma_3 = d\psi + \cos\theta\,d\phi$, with $d\sigma_i = \tfrac12\epsilon_{ijk}\sigma_j\wedge\sigma_k$; Tong's $\sigma_1$ and $\sigma_2$ have the opposite signs, and the line element, which holds their squares, is the same.
Self duality is Ghezelbash and Mann's `conditions`, $a' = f((b - c)^2 - a^2)/(2bc)$ and its two cyclic permutations, which is Tong's $(2bc/f)\,da/dr = (b - c)^2 - a^2$ and Hanany and Pioline's `sdu` of section 3.1 with $\lambda = 1$ and $f = abc$.
The space is hyperkähler, so its curvature is self dual in one orientation and it is Ricci flat.

## Step 3: the modulus chart

Gibbons and Manton fix the radius by $f = -b/r$ (Weinberg and Yi, the subsection "Atiyah-Hitchin geometry for two identical monopoles", "Following Gibbons and Manton"; Schroers and Singer, the subsection "The Atiyah-Hitchin manifold and related spaces", "As explained in [GM]"), and write the solution in an angle $\beta$ with $r = 2K(\sin(\beta/2))$ (Weinberg and Yi, `rAndBeta`):
$ab = -\sin\beta\,r\,dr/d\beta + \tfrac12(1 - \cos\beta)r^2$, $bc = -\sin\beta\,r\,dr/d\beta - \tfrac12(1 + \cos\beta)r^2$, $ca = -\sin\beta\,r\,dr/d\beta$.
Ghezelbash and Mann's `abc`, `GAMMA` and `gama` are the same solution with the length $n$ that sets the scale, $K(\sin(\gamma/2)) = r/2n$, and $a$, $b$ and $c$ written out; at $n = 1$ their $a$, $b$ and $c$ equal Weinberg and Yi's products at 30 digits, and they solve the three equations of self duality at $10^{-28}$.
In the modulus $k = \sin(\beta/2)$, with $dr/d\beta = n(E - (1 - k^2)K)/(k\sqrt{1 - k^2})$ from $dK/dk = (E - (1 - k^2)K)/(k(1 - k^2))$, the three products become
$bc = -4n^2KE$, $ca = -4n^2K(E - (1 - k^2)K)$, $ab = 4n^2K(K - E)$,
with $K = \mathrm{K}(k^2)$ and $E = \mathrm{E}(k^2)$ the complete elliptic integrals at the parameter $k^2$, an algebraic rewriting checked at 30 digits; the elliptic integrals stay as they are, and no closed form is put for them.
Hanany and Pioline's appendix A gives Atiyah and Hitchin's own route to the same solution, $u = \sqrt{2\sin\theta}\,K(\sin(\theta/2))$, whose $w_i$ are these products halved, a constant scaling the equations allow.
So $a^2 = 4n^2K(K - E)(E - (1 - k^2)K)/E$, $b^2 = 4n^2KE(K - E)/(E - (1 - k^2)K)$, $c^2 = 4n^2KE(E - (1 - k^2)K)/(K - E)$, and $f^2dr^2 = (b^2/r^2)(dr/dk)^2dk^2$ with $g_{kk} = 4n^2E(K - E)(E - (1 - k^2)K)/(k^2(1 - k^2)^2K)$.
On $0 < k < 1$, $E > (1 - k^2)K$ and $K > E$, so each factor is positive and $a > 0$, $b > 0$, $c < 0$, Hanany and Pioline's signs.
The modulus as a coordinate on the space is Bakas's: his curve has modulus $k$, the Taub-NUT limit at $k \to 1$ and the bolt at $k \to 0$.
$K$ and $E$ are held as functions of $k$ with the rates $dK/dk = (E - (1 - k^2)K)/(k(1 - k^2))$ and $dE/dk = (E - K)/k$, which the reader checks against their definitions, so every value is a rational function of $k$, $K$ and $E$, among which no relation is left, and the Ricci tensor is exactly zero.
The check holds $g_{kk}$ and $g_{\psi\psi}$ to Weinberg and Yi's products evaluated independently with mpmath at three points, and the Ricci tensor to zero there.

## Step 4: the radius chart

Gibbons and Manton's radius $r$, with $f = -b/r$, is the chart every source writes the metric in, and in it $a$, $b$ and $c$ are fixed by $\beta(r)$, the inverse of $r = 2nK(\sin(\beta/2))$, for which neither sympy nor the reader has a closed form.
So the chart leaves $a$, $b$ and $c$ free functions of $r$, as Mixmaster's chart leaves its three scale factors, and no component assumes the equations of self duality, which its convention states:
with $f = -b/r$ they are $2rc\,a' = a^2 - (b - c)^2$, $2rac\,b' = b(b^2 - (c - a)^2)$ and $2rab\,c' = b(c^2 - (a - b)^2)$.
The check holds its Ricci tensor to zero at three points of the space, with $a$, $b$ and $c$ from Weinberg and Yi's products and their first and second derivatives carried along by those equations, and holds $g_{rr}(dr/dk)^2$ to the modulus chart's $g_{kk}$ there.
`metric_tags.py` sets `vacuum` in `OVERRULED`, as for Mixmaster, since the free chart is a vacuum only on the solution.
`metric_tags.py` samples a held name as a free positive number, so in the modulus chart $K$ and $E$ take values unrelated to $k$ and $\partial_\phi$ reads as timelike somewhere; with $K$ and $E$ the elliptic integrals of $k$ every coefficient is positive, and since $\phi$ is no translation no tag follows from it.

## Step 5: the core and the identifications

Near $r = \pi n$, $a = 2(r - \pi n)$, $b = \pi n + (r - \pi n)/2$ and $c = -\pi n + (r - \pi n)/2$ to first order, Ghezelbash and Mann's `abcofr`, checked numerically; Weinberg and Yi's `metricnearcenter` drops the factor $\pi^2$ of $\sigma_2^2 + \sigma_3^2$.
The circle of $\sigma_1$ closes smoothly only if its angle has period $\pi$, which is the identification $(\theta, \phi, \psi) \sim (\pi - \theta, \phi + \pi, -\psi)$, Weinberg and Yi's $I$, the exchange of the two monopoles; the core is then a 2-sphere of radius $\pi n$ and the space is the double cover.
The further identification $\psi \sim \psi + \pi$, Weinberg and Yi's $I'$, gives the Atiyah-Hitchin manifold $M_2^0$, with a real projective plane at the core and fundamental group $\mathbb{Z}_2$ (Seiberg and Witten, the passage "Comparison To Exact Metrics"); the domains state both, with $\psi \in [0, \pi)$.

## Step 6: far from the core

$a, b \to r\sqrt{1 - 2n/r}$ and $c \to -2n/\sqrt{1 - 2n/r}$ up to $O(e^{-r/n})$, Ghezelbash and Mann's `abcatinfinity`, Schroers and Singer's `AHasy`, Tong's `abc` in the subsection "The Moduli Space" of the monopole lecture at $n = 1$ with $a^2 - r^2(1 - 2/r) \approx -8r^3e^{-r}$, checked numerically.
That is the self dual Taub-NUT metric with a negative mass, Gibbons and Manton (1995) in their introduction and Hanany and Pioline after `rt`, singular at $r = 2n$, which lies inside the core radius $\pi n$.

## Step 7: the chart not printed

Hanany and Pioline's chart in $t$, $ds^2 = (abc)^2dt^2 + a^2\sigma_1^2 + b^2\sigma_2^2 + c^2\sigma_3^2$, their `ds`, carries the same space, with $w_1 = bc$, $w_2 = ca$ and $w_3 = ab$ "Following [AH]" and the Halphen system `dw`, $d(w_1 + w_2)/dt = -2w_1w_2$ and its cyclic permutations.
Their theta-function solution `ww123` is Gibbons and Manton's products divided by $2\pi^2$, at the modulus $k = \theta_2^2/\theta_3^2$, at 30 digits, and it solves $d(w_i + w_j)/dt = +2\pi w_iw_j$, not their `dw` as printed: the sign and the factor $\pi$ disagree, checked at three values of $t$ with mpmath's theta functions.
Their appendix A writes the infinite separation as $k \to \infty$, which is $k \to 1$.
A chart whose defining equations would have to be corrected from the source is left out; the modulus chart carries Atiyah and Hitchin's solution exactly, and the radius chart the form every source writes.
The Kähler coordinates of Ivanov and Rocek, Bakas, Ionas and Arai, Baba and Ionas give the metric through contour integrals and ratios of Hankel determinants, with no line element in closed form, and are not a chart here either.

## Step 8: the drawings

`atiyah_hitchin_drawings.py` draws them, and records in `NOT_DRAWN` for `null_rays.py` and `conformal.py` that a positive definite metric has no nonzero null tangent and no causal boundary.
Each of the three circles on the orbits closes after an angle $\pi$ once both identifications are made: $\psi$ by $I'$, the circle of $\sigma_1$ at $\theta = \pi/2$, $\psi = 0$ by $I$, which sends $\phi$ to $\phi + \pi$ there, and the circle of $\sigma_2$, $\theta$ at $\phi = \psi = 0$, by $I$ followed by $I'$, which sends the rotation by $\theta$ about the second axis to the rotation by $\theta + \pi$.
So each circle is drawn swept by twice its angle, and its radius in the flat space of the drawing is half its coefficient: $a/2$, $b/2$ or $|c|/2$.
The embedding diagrams are drawn in the radius chart with $a$, $b$ and $c$ given at each $r$ by Gibbons and Manton's solution, the root $k$ of $2K(k^2) = r/n$ found with mpmath and the products of Step 3, and their slopes by the equations of self duality, since the modulus chart's $g_{kk}$ vanishes at the core as $\pi^2n^2k^2/4$.
The first is the surface of $r$ and $\phi$ at $\theta = \pi/2$, $\psi = 0$, Atiyah and Hitchin's surface of two monopoles moving in one plane as Bachmaier, Dvali, Seitz and Valbuena-Bermudez describe it: $\rho = a/2$, which grows as the distance from the core does, $a = 2(r - \pi n)$ there, so the surface is smooth, and far away grows at half that rate, a cone of half angle $30°$, whose geodesic through the tip turns by $\pi/2$.
The second is the trumpet through the core at $\phi = 0$: below the core the surface of $r$ and $\psi$ at $\theta = \pi/2$, $\rho = |c|/2$ falling to $n$, a tube; above it the surface of $r$ and $\theta$ at $\psi = 0$, $\rho = b/2$ growing to the cone of $30°$; the two meet on the core's circle, $\rho = \pi n/2$, where $b = |c| = \pi n$ and both slopes are $1/4$ along the proper distance, so the surface is smooth there.
The curvature figure is the stored Kretschmann scalar of the modulus chart against $k$, from $k = 0.02$, where $r - \pi n = 3 \times 10^{-4}n$, to $k = 0.9995$, where $r = 9.7n$.
That scalar is a rational function of $k$, $K$ and $E$ whose numerator and denominator vanish to high order at the core, since $K - E$ and $E - (1 - k^2)K$ are of order $k^2$ there: evaluated in doubles it loses every digit below $k = 0.6$, so `Precise` evaluates the stored scalar itself in sixty digits, which holds it to a part in $10^{19}$ at $k = 0.02$.
It is $16(\lambda_1^2 + \lambda_2^2 + \lambda_3^2)$ with $\lambda_i$ the three frame components $R_{0i0i}$, whose sum vanishes, as self duality asks; that form is no shorter, so the scalar is printed as the printer writes it.
