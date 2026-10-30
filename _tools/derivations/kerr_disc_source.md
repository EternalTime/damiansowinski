# Israel's disc source of the Kerr metric

Kerr's metric in the coordinates of Boyer and Lindquist,

$$ds^2 = -\left(1 - \frac{2mr}{\Sigma}\right)c^2dt^2 - \frac{4mar\sin^2\theta}{\Sigma}c\,dt\,d\phi + \frac{\Sigma}{\Delta}dr^2 + \Sigma\,d\theta^2 + \left(r^2 + a^2 + \frac{2ma^2r\sin^2\theta}{\Sigma}\right)\sin^2\theta\,d\phi^2,$$

$$\Sigma = r^2 + a^2\cos^2\theta,\qquad \Delta = r^2 - 2mr + a^2,$$

kept on $r \ge 0$ alone, with the disc $r = 0$ a surface whose two faces, $\theta$ and $\pi - \theta$, are one.
That is Werner Israel's reading of the metric (`israel1970`): the jump in the normal derivatives across the glued disc is the field of a layer of matter on it.

Both charts are written by `_tools/derivations/print_charts.py --metric kerr_disc_source`, about fifteen seconds for the two on 29 October 2026, and `verify_metrics.py --metric kerr_disc_source` checks them in about ten.

## Step 1. The charts

The first chart is Boyer and Lindquist's, the `kerr` page's chart with $GM/c^2$ written $m$, as `kerr_taub_nut` writes it.
The second is Kerr and Schild's form in spheroidal coordinates,

$$ds^2 = -c^2d\tau^2 + \frac{\Sigma}{r^2 + a^2}dr^2 + \Sigma\,d\theta^2 + \left(r^2 + a^2\right)\sin^2\theta\,d\psi^2 + \frac{2mr}{\Sigma}\left(c\,d\tau + \frac{\Sigma}{r^2 + a^2}dr - a\sin^2\theta\,d\psi\right)^2,$$

which is the Kerr-Schild chart of `kerr_de_sitter` at $\Lambda = 0$, with $c\,d\tau = c\,dt + 2mr\,dr/\Delta$ and $d\psi = d\phi + 2mar\,dr/((r^2 + a^2)\Delta)$.
Its first four terms are flat space in oblate spheroidal coordinates, so at $m = 0$ the chart is Minkowski space and the disc $r = 0$ is the flat disc $z = 0$, $x^2 + y^2 < a^2$, with its rim the circle $x^2 + y^2 = a^2$.
The chart runs through both horizons and reaches the disc, which Boyer and Lindquist's does too for $a > m$.

`kerr_disc_source_check` holds each chart, at random rational points, to $R_{\mu\nu} = 0$ off the disc and to $\theta \to \pi - \theta$ carrying $g_{\mu\nu}$ and $\partial_r g_{\mu\nu}$ at $r = 0$ into themselves, so the gluing is an isometry and the two faces see the same extrinsic curvature; it holds the Kerr-Schild chart to Boyer and Lindquist's pulled back.
`KerrDiscForms` writes every value in $\Sigma$ and $\Delta$: Kerr-de Sitter's forms at $\Lambda = 0$, grouping a leftover sum by the powers of $m$ or of $\Delta$.

## Step 2. The layer on the disc

Take the unit normal on each face pointing into $r > 0$.
Then $n = \sqrt{g^{rr}}\,\partial_r$ with $g_{rr} = \Sigma/\Delta = \cos^2\theta$ at $r = 0$, and the extrinsic curvature is $K_{ij} = \partial_r g_{ij}/(2\sqrt{g_{rr}})$ on the face $\cos\theta > 0$.
Crossing the disc from the lower face to the upper, the normal of the lower face is reversed, so the jump is $[K_{ij}] = 2K_{ij}$, and Lanczos's equation $S_{ij} = -\frac{1}{8\pi}([K_{ij}] - h_{ij}[K])$ with $G = c = 1$ gives $S_{ij} = -\frac{1}{4\pi}(K_{ij} - h_{ij}K)$.

The induced metric is $h = \mathrm{diag}(-1, a^2\cos^2\theta, a^2\sin^2\theta)$ in $(t, \theta, \phi)$, a flat disc in polar coordinates with radius $\rho = a\sin\theta$, and

$$K_{tt} = \frac{m}{a^2\cos^3\theta},\qquad K_{t\phi} = -\frac{m\sin^2\theta}{a\cos^3\theta},\qquad K_{\phi\phi} = \frac{m\sin^4\theta}{\cos^3\theta},$$

the others zero, so $K = -m/(a^2\cos\theta)$.
The layer is

$$S^i{}_j = \lambda\left(\hat w^i\hat w_j + \hat e^i\hat e_j\right),\qquad \lambda = -\frac{m}{4\pi a^2\cos\theta} = -\frac{m}{4\pi a\sqrt{a^2 - \rho^2}},$$

with $\hat e$ the unit radial vector and $\hat w$ the unit spacelike vector along $\partial_t + \Omega\,\partial_\phi$, $\Omega = 1/(a\sin^2\theta) = a/\rho^2$.
The stream moves at $\rho\Omega = a/\rho$ times the speed of light: light speed on the rim, faster inside it, and without bound toward the centre.
The other eigenvector of the $(t, \phi)$ block is the timelike vector turning rigidly at $\Omega = 1/a$, with eigenvalue zero, so an observer turning with it measures zero energy density.

The Tolman mass, $\int(-S^t{}_t + S^\theta{}_\theta + S^\phi{}_\phi)\sqrt{-h}\,d\theta\,d\phi$, has density

$$\sigma = -\frac{ma}{2\pi\left(a^2 - \rho^2\right)^{3/2}}$$

per unit area of the flat disc, and the angular momentum $\int S^t{}_\phi\sqrt{-h}\,d\theta\,d\phi$ has $dJ/d\rho = -m\rho^3/(2(a^2 - \rho^2)^{3/2})$.
Both are negative and diverge at the rim, so the ring carries an infinite positive mass and angular momentum, and the totals are $m$ and $ma$.
These agree with equations (74) and (82) of Bianchi, Gambino, Pani and Riccioni (`bianchi2025`) at every point compared, to fifteen digits, and `_tools/test_kerr_disc_source.py` holds the layer to the charts.
At $m = 0$ the layer vanishes and the spacetime is flat, as Gibbons and Volkov note (`gibbons2017`).

## Step 3. What is left out

The sheet of negative $r$ is thrown away.
For $r > 0$, $g_{\phi\phi} = (r^2 + a^2 + 2ma^2r\sin^2\theta/\Sigma)\sin^2\theta > 0$, so no circle of constant $t$, $r$ and $\theta$ is timelike; Carter's closed timelike curves near the ring (`carter1968`) lie at $r < 0$.

## Step 4. The drawings

The page draws the case $a > m$ with $m = 1$, $a = 3/2$, where no horizon hides the disc, so every view shows the disc and its rim.
