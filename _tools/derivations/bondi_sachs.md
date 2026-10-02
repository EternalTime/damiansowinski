# Bondi and Sachs, the metric of an isolated body that radiates

Bondi, van der Burg and Metzner's metric of 1962 is, in the signature of this collection,

$$ds^2 = -\frac{V}{r}e^{2\beta}c^2du^2 - 2e^{2\beta}c\,du\,dr + r^2e^{2\gamma}\left(d\theta - U\,c\,du\right)^2 + r^2e^{-2\gamma}\sin^2\theta\,d\phi^2,$$

with $u$ constant on each outgoing light cone, $\theta$ and $\phi$ constant along each ray, and $r$ the luminosity distance, for which the sphere of constant $u$ and $r$ has the area $4\pi r^2$.
Its two charts are written by `_tools/derivations/print_charts.py --metric bondi_sachs` in about four minutes, and `verify_metrics.py --system bondi_sachs/<chart>` checks both in under a minute.

## Step 1. The charts

The chart of $u$, $r$, $\theta$ and $\phi$ is Bondi, van der Burg and Metzner's own, with their $V$, $\beta$, $U$ and $\gamma$ functions of $u$, $r$ and $\theta$ (doi:10.1098/rspa.1962.0161).
They write the signature $(+,-,-,-)$ and the term in $du^2$ multiplied out, $\left(Vr^{-1}e^{2\beta} - U^2r^2e^{2\gamma}\right)du^2$; the line element above is minus theirs with the square completed, which is the form Sachs's general metric takes.
Their source has an axis of symmetry and is unchanged by $\phi \to -\phi$, so it does not rotate and its waves have one polarisation.

The chart of $u$, $\ell = 1/r$, $\theta$ and $\phi$ is the one in which Penrose's rescaled metric is written, as in Mädler and Winicour's review, arXiv:1609.01731, section 4.
The chart printed is of the physical metric, and $\ell^2$ times it is the rescaled one, regular at $\ell = 0$.
The review's equation for the rescaled metric prints the term in $du^2$ without its minus sign; the printer checks the chart to be Bondi's pulled back along $r = 1/\ell$, which fixes the sign.

Sachs's general metric (doi:10.1098/rspa.1962.0206) has six functions of four coordinates,

$$ds^2 = -\frac{V}{r}e^{2\beta}du^2 - 2e^{2\beta}du\,dr + r^2h_{AB}\left(dx^A - U^Adu\right)\left(dx^B - U^Bdu\right),$$

with $2h_{AB}dx^Adx^B = \left(e^{2\gamma} + e^{2\delta}\right)d\theta^2 + 4\sin\theta\sinh(\gamma - \delta)\,d\theta\,d\phi + \sin^2\theta\left(e^{-2\gamma} + e^{-2\delta}\right)d\phi^2$.
It is not printed: its tensors did not finish in a quarter of an hour on 2 October 2026, and Bondi's chart, with four functions of three coordinates, already takes 870 kilobytes.
The History states Sachs's metric.
Newman and Unti's chart, with an affine distance for $r$, belongs to the spin coefficient formalism and is not printed either.

## Step 2. What the printer checks

Both charts leave $V$, $\beta$, $U$ and $\gamma$ free, so no component assumes a field equation.
Before anything is written `bondi_sachs_check` holds Bondi's chart to the members the relations name, each at its own functions:

- $V = r - r_s$, the rest zero: the published outgoing Eddington-Finkelstein chart of `schwarzschild`;
- $V = r - 2Gm(u)/c^2$: the published chart of `vaidya`;
- $U = -\alpha\sin\theta$, $V = r - 2m - 2\alpha r^2\cos\theta$, $\beta = \gamma = 0$: the published rectilinear chart of `photon_rocket`;
- $\gamma = 0$, $e^{2\beta} = f$, $U = \partial_\theta f/r$, $V = r\left(2H + 2r\,\partial_uf + (\partial_\theta f)^2\right)/f$: the published axisymmetric chart of `robinson_trautman` carried along $r \to f\,r$, so Robinson and Trautman's affine distance is $f$ times the luminosity distance, and with their vacuum $H$ the function $V = r\left(K + (\partial_\theta f)^2\right)/f - 2m/f^2$ is linear in $r$;
- $V = r$: flat.

It also holds $R_{rr} = 4\,\partial_r\beta/r - 2\left(\partial_r\gamma\right)^2$, the first of Bondi's main equations, which is the vacuum equation the description of $\beta$ states.
The Kretschmann scalar is printed by `bondi_sachs_scalar` as one expanded sum over its denominator, since sympy's `factor` did not finish on it in a quarter of an hour.

## Step 3. Bondi's expansion, and the burst drawn

With $\sigma(u, \theta)$ for Bondi's $c$, which this collection keeps for the speed of light, the expansion in $1/r$ of the vacuum solution is, in units $G = c = 1$,

$$\gamma = \frac{\sigma}{r} + \frac{C - \sigma^3/6}{r^3},\qquad \beta = -\frac{\sigma^2}{4r^2},\qquad U = -\frac{\partial_\theta\sigma + 2\sigma\cot\theta}{r^2} + \frac{2N + 3\sigma\,\partial_\theta\sigma + 4\sigma^2\cot\theta}{r^3},$$

$$V = r - 2M - \frac{\partial_\theta N + N\cot\theta - (\partial_\theta\sigma)^2 - 4\sigma\,\partial_\theta\sigma\cot\theta - \tfrac{1}{2}\sigma^2\left(1 + 8\cot^2\theta\right)}{r},$$

with the three supplementary equations

$$\partial_uM = -(\partial_u\sigma)^2 + \tfrac{1}{2}\partial_u\left(\partial_\theta^2\sigma + 3\,\partial_\theta\sigma\cot\theta - 2\sigma\right),\qquad -3\,\partial_uN = \partial_\theta M + 3\sigma\,\partial_u\partial_\theta\sigma + 4\sigma\,\partial_u\sigma\cot\theta + \partial_u\sigma\,\partial_\theta\sigma,$$

$$4\,\partial_uC = 2\sigma^2\partial_u\sigma + 2\sigma M + N\cot\theta - \partial_\theta N.$$

The paper itself was not to hand, so every coefficient was derived again on 2 October 2026: with the coefficients of the nonlinear terms left free, the series of $R_{rr}$, $R_{r\theta}$, $R_{ur}$, $R_{\theta\theta}$ and $R_{\phi\phi}$ in $1/r$ vanish through the orders kept only at the values above, and the equation for $\partial_uC$ is what $R_{\theta\theta}$ leaves at order $r^{-2}$.
The mean of the first supplementary equation over the sphere is the mass loss formula, $dm_B/du = -\tfrac{1}{2}\int_0^\pi(\partial_u\sigma)^2\sin\theta\,d\theta$, since the second term is a divergence.

The drawings take a weak burst, `bondi_burst` in `null_rays.py`: $\sigma = \partial_u^2Q\,\sin^2\theta$ with $Q = \tfrac{1}{4}m_0^3\sin^6(\pi cu/20m_0)$ from $u = 0$ to $cu = 20\,m_0$, a quadrupole moment that makes one swing and comes back.
For that shear the three equations integrate in closed form, with $E = \int(\partial_u^3Q)^2du$:

$$M = m_0 + 2\,\partial_u^2Q\left(3\cos^2\theta - 1\right) - E\sin^4\theta,\qquad N = 4\,\partial_uQ\,\sin\theta\cos\theta + \left(\tfrac{4}{3}\textstyle\int E\,du - 2(\partial_u^2Q)^2\right)\sin^3\theta\cos\theta,$$

and the coefficient of $r^{-3}$ in $\gamma$ follows by one more integration.
Before the burst the metric is Schwarzschild's of the mass $m_0$.
The Bondi mass falls by $\tfrac{8}{15}E = 0.00102\,m_0$.
After the burst the mass aspect is $m_0 - E\sin^4\theta$, uneven, and $N$ grows linearly in $u$ away from the axis and the equator's plane keeps a term in $V$ that grows with it: with no news to carry the unevenness off, the source the expansion describes is in one of the non-radiative motions Bondi, van der Burg and Metzner discuss.
On the axis $\sigma$ and that term vanish, and the plane is Schwarzschild's of the mass $m_0$ again, which is why the conformal diagram draws the axis.

A truncated expansion is no exact solution.
The Ricci tensor those orders leave, measured component by component in an orthonormal frame against the square root of the Kretschmann scalar, is 1.4% at $r = 10\,m_0$, 0.7% at $20\,m_0$ and 0.35% at $40\,m_0$, falling as $1/r$; `_tools/test_bondi_sachs.py` holds the published Ricci tensor to those numbers.
A burst strong enough to see by eye is out of reach of the expansion: at four times the amplitude the residue is 9.5% at $r = 10\,m_0$, and a strain of a tenth at a radius where the series holds would take more mass than the source has.
So the drawings show the causal picture, the cones, the world tube and null infinity, and state the mass loss and the strain in numbers.

## Step 4. The drawings

The spacetime diagrams are the plane of $u$ and $r$ on the equator and on the axis, where $U$ vanishes and the null curves drawn are null geodesics, from the world tube $r = 10\,m_0$ out, and the plane of $u$ and $\ell$ on the equator from $r = 8\,m_0$ to the edge $\ell = 0$, future null infinity.
In the chart of $\ell$ the ingoing rays obey $d\ell/d(cu) = \ell^3V/2$, so the future cone closes onto the edge.

The conformal diagram is the axis outside the world tube, `BurstPlane` in `conformal.py`: $p = \arctan(cu/40m_0)$, and $q = \arctan(cv/40m_0)$ with $v$ the advanced time an ingoing ray has before the burst, where the plane is Schwarzschild's, found for a ray met later by carrying it back.

The embedding diagram is the sphere of constant $u$ and $r$ at the world tube, $r = 10\,m_0$, in the middle of the burst, $cu = 10\,m_0$, where the shear is greatest.
Its metric is $r^2\left(e^{2\gamma}d\theta^2 + e^{-2\gamma}\sin^2\theta\,d\phi^2\right)$, a surface of revolution with circles of radius $r\,e^{-\gamma}\sin\theta$ and the area $4\pi r^2$ whatever $\gamma$ is.
There $\gamma = -0.0035$ on the equator, so the equator's radius is $1.0035\,r$ and the meridian from pole to pole is $0.9983\,\pi r$, an oblate sphere flattened by half a percent.
It is one moment and no movie, since the other moments differ from it by less than the eye can tell.
`embedding.py` hands $\gamma$ and its slope to the slice as numbers: sympy did not finish simplifying $e^{2\gamma}$ of the closed form in ten minutes.
