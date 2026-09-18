# The Kerr rotating black hole

This is the working behind the Boyer-Lindquist coordinate system in `MFS/assets/data/metrics/kerr.json`.
Every number the entry prints is derived here, in order, from the line element down to the geodesic equations.
Nothing is left as an exercise and nothing is asserted that is not computed.

The companion script `_tools/derivations/verify_metrics.py` does the same work in sympy and compares it against the published file, so the algebra below is checkable by hand and by machine independently.

Kerr carries the longest connection in the collection: twenty independent Christoffel symbols with an upper index, twenty more with all three lowered, and thirteen independent components of the Riemann tensor.
Two things keep that from being unreadable.
The first is that the whole connection is built out of the two abbreviations $\Sigma$ and $\Delta$ and the one combination $r^2 - a^2\cos^2\theta$, which is Step 4.
The second is that the entire curvature is two real functions, $\mathcal{A}$ and $\mathcal{B}$, which are the real and imaginary parts of a single cube, and every one of the thirteen Riemann components is one of those two multiplied by something algebraic.
Step 7 finds them and Step 11 squares them to get the Kretschmann scalar in one line.

---

## Step 1. Conventions

The signature is $(-,+,+,+)$ and the Christoffel symbols are those of the Levi-Civita connection,

$$\Gamma^\mu_{\nu\rho} = \tfrac{1}{2}g^{\mu\alpha}\left(\partial_\nu g_{\alpha\rho} + \partial_\rho g_{\alpha\nu} - \partial_\alpha g_{\nu\rho}\right),$$

with the lowered symbol $\Gamma_{\mu\nu\rho} = g_{\mu\alpha}\Gamma^\alpha_{\nu\rho}$.
The Riemann tensor is

$$R^\mu{}_{\nu\rho\sigma} = \partial_\rho \Gamma^\mu_{\nu\sigma} - \partial_\sigma \Gamma^\mu_{\nu\rho} + \Gamma^\mu_{\rho\lambda}\Gamma^\lambda_{\nu\sigma} - \Gamma^\mu_{\sigma\lambda}\Gamma^\lambda_{\nu\rho},$$

which is what the published Riemann components are in.

The Ricci tensor is the standard contraction,

$$R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu},$$

settled for the collection on 2026-09-18.
For this spacetime the choice costs nothing, because Step 8 shows that both contractions vanish, and a tensor that is zero is zero with either sign in front of it.
The entry is written the standard way regardless.

Factors of $c$ and $G$ are kept explicit.
Kerr is one of the two entries in the collection that keeps a mass as a mass rather than folding it into a length, so $M$ appears with Newton's constant beside it and the Schwarzschild radius of that mass, $r_s = 2GM/c^2$, is spelled out as $2GM/c^2$ in every published value.
Steps 2 to 12 use the shorthand $r_s$, because carrying $2GM/c^2$ through eighty lines of algebra obscures the algebra; Step 13 puts it back.

The chart, here and everywhere else in the collection, is the one whose time coordinate is $x^0 = cT$.
Schwarzschild shows the convention plainly, quoting $g_{tt} = -(1-r_s/r)$ against a line element whose time term is $-(1-r_s/r)c^2dt^2$.
Since $t$ here is a time, the chart coordinate is $ct$, and every component printed against an index $t$ is a component in that chart even though the index is written with the bare letter.

Because the rescaling $t \to ct$ is linear with constant coefficients, the rule is arithmetic: a component in the chart is the component taken with the bare coordinate, multiplied by $c$ once for every upper $t$ index and divided by $c$ once for every lower one.
The Christoffel symbols obey the same rule as the tensors, since the inhomogeneous term in their transformation law carries a second derivative of the coordinate change, which vanishes for a linear one.
A scalar has no index at all and is therefore the same number in both charts, which is what lets Step 11 compute the Kretschmann scalar without thinking about the chart.

Steps 2 to 12 work in the bare $(t,r,\theta,\phi)$ chart, where $g_{tt}$ carries a $c^2$, because that is where the computation is done.
Step 13 carries the results into the chart the entry prints, and the arithmetic there is the check that no factor of $c$ was dropped: every $c$ that appears in the bare tables of Steps 3 to 7 turns out to be exactly $c$ to the power of the number of lower $t$ indices, so the published tables are the same formulas with the $c$ powers struck out.

---

## Step 2. The line element

$$ds^2 = -\left(1 - \frac{r_s r}{\Sigma}\right)c^2dt^2 - \frac{2r_s r a\sin^2\theta}{\Sigma}c\,dt\,d\phi + \frac{\Sigma}{\Delta}dr^2 + \Sigma\,d\theta^2 + \left(r^2 + a^2 + \frac{r_s r a^2\sin^2\theta}{\Sigma}\right)\sin^2\theta\,d\phi^2,$$

with

$$\Sigma = r^2 + a^2\cos^2\theta, \qquad \Delta = r^2 - r_s r + a^2, \qquad r_s = \frac{2GM}{c^2}, \qquad a = \frac{J}{Mc}.$$

The two parameters are the mass $M$ and the spin per unit mass $a$, and $a$ is a length: $J$ is an angular momentum, $Mc$ is a momentum, and their ratio is a length.
That is what makes $\Sigma$ and $\Delta$ areas, and it is the single fact the dimensional pass of Step 15 leans on.

Setting $a = 0$ gives $\Sigma = r^2$ and $\Delta = r^2 - r_sr$, and the line element collapses to Schwarzschild in the usual spherical chart.
Setting $M = 0$ with $a$ held fixed gives $\Delta = r^2 + a^2$ and a flat metric written in oblate spheroidal coordinates, which is worth knowing because it says that $a$ on its own bends nothing; it is the product $Ma$ that carries the rotation.

The metric depends on $r$ and $\theta$ and on neither $t$ nor $\phi$, so $\partial_t$ and $\partial_\phi$ are Killing vectors: the solution is stationary and axisymmetric.
It is not static, because $g_{t\phi} \ne 0$ and no relabelling of $\phi$ alone can remove a term that changes sign under $t \to -t$ while $\phi$ is held.
That single off diagonal entry is frame dragging, and it is the reason the connection has twice as many symbols as Schwarzschild's.

Three derivatives do most of the work below,

$$\partial_r\Sigma = 2r, \qquad \partial_\theta\Sigma = -2a^2\sin\theta\cos\theta, \qquad \partial_r\Delta = 2r - r_s,$$

and $\partial_\theta \Delta = 0$.

The chart covers the exterior region.
$\Delta$ vanishes at

$$r_\pm = \frac{r_s}{2} \pm \sqrt{\frac{r_s^2}{4} - a^2} = \frac{GM}{c^2} \pm \sqrt{\frac{G^2M^2}{c^4} - a^2},$$

which are real and distinct when $a < GM/c^2$, and the domain published is $r > r_+$.
Step 12 shows that neither root is a curvature singularity, so both are places the chart fails rather than places the spacetime does.

---

## Step 3. The metric matrix, its determinant and its inverse

In the order $(t,r,\theta,\phi)$ the matrix is block diagonal, with a two by two block in the $(t,\phi)$ slots and two diagonal entries between them,

$$g_{\mu\nu} = \begin{pmatrix} -c^2\left(1 - \dfrac{r_sr}{\Sigma}\right) & 0 & 0 & -\dfrac{c\,r_sra\sin^2\theta}{\Sigma} \\[2mm] 0 & \dfrac{\Sigma}{\Delta} & 0 & 0 \\[2mm] 0 & 0 & \Sigma & 0 \\[2mm] -\dfrac{c\,r_sra\sin^2\theta}{\Sigma} & 0 & 0 & \left(r^2+a^2+\dfrac{r_sra^2\sin^2\theta}{\Sigma}\right)\sin^2\theta \end{pmatrix}.$$

The cross term is read off the line element with the factor of two divided out, since $g_{t\phi}dtd\phi + g_{\phi t}d\phi dt = 2g_{t\phi}dtd\phi$.

The determinant of the two by two block is the only piece that needs work.
Write $\mathcal{D}$ for it, and note first that $g_{tt} = -c^2(\Sigma - r_sr)/\Sigma$ and $g_{\phi\phi} = \left[(r^2+a^2)\Sigma + r_sra^2\sin^2\theta\right]\sin^2\theta/\Sigma$.
Then

$$\mathcal{D} = g_{tt}g_{\phi\phi} - g_{t\phi}^2 = -\frac{c^2\sin^2\theta}{\Sigma^2}\Big\{\left(\Sigma - r_sr\right)\left[(r^2+a^2)\Sigma + r_sra^2\sin^2\theta\right] + r_s^2r^2a^2\sin^2\theta\Big\},$$

where the last term inside the brace is $-g_{t\phi}^2$ carried across the common factor, and so enters with a plus sign.
Expanding the product gives four terms, of which the one in $r_s^2$ is $-r_s^2r^2a^2\sin^2\theta$, and it cancels the term that was just carried across.
The remaining three all have a factor of $\Sigma$:

$$\mathcal{D} = -\frac{c^2\sin^2\theta}{\Sigma}\left[(r^2+a^2)\Sigma + r_sra^2\sin^2\theta - r_sr(r^2+a^2)\right],$$

and inside the bracket the two $r_s$ terms combine as $-r_sr\left[(r^2+a^2) - a^2\sin^2\theta\right] = -r_sr\Sigma$, so the bracket is $\Sigma\left(r^2+a^2-r_sr\right) = \Sigma\Delta$.
Hence

$$\mathcal{D} = -c^2\Delta\sin^2\theta, \qquad \det g = g_{rr}\,g_{\theta\theta}\,\mathcal{D} = \frac{\Sigma}{\Delta}\cdot\Sigma\cdot(-c^2\Delta\sin^2\theta) = -c^2\Sigma^2\sin^2\theta.$$

The determinant is negative wherever $\Sigma \ne 0$ and $\sin\theta \ne 0$, so the signature is Lorentzian on the whole chart, and $\sqrt{-g} = c\,\Sigma\sin\theta$.
What it does not contain is $\Delta$: the determinant is perfectly finite at both horizons, which is the first hint that they are not singularities.

Inverting a block diagonal matrix is inverting each block.
The two diagonal entries invert on the spot, and the two by two block inverts by the usual formula with $\mathcal{D}$ underneath:

$$g^{tt} = \frac{g_{\phi\phi}}{\mathcal{D}} = -\frac{(r^2+a^2)\Sigma + r_sra^2\sin^2\theta}{c^2\Sigma\Delta}, \qquad g^{t\phi} = -\frac{g_{t\phi}}{\mathcal{D}} = -\frac{r_sra}{c\,\Sigma\Delta}, \qquad g^{\phi\phi} = \frac{g_{tt}}{\mathcal{D}} = \frac{\Sigma - r_sr}{\Sigma\Delta\sin^2\theta},$$

together with $g^{rr} = \Delta/\Sigma$ and $g^{\theta\theta} = 1/\Sigma$.

The numerator of $g^{tt}$ has a second form that the entry uses, because it is the one that shows the function is not singular at the horizon by accident:

$$(r^2+a^2)\Sigma + r_sra^2\sin^2\theta = (r^2+a^2)\left(r^2+a^2-a^2\sin^2\theta\right) + r_sra^2\sin^2\theta = (r^2+a^2)^2 - a^2\Delta\sin^2\theta.$$

The same rearrangement run backwards gives $\Sigma - r_sr = \Delta - a^2\sin^2\theta$, which is the numerator of $g^{\phi\phi}$ written the other way.
Both identities are used again in Steps 7 and 8.

---

## Step 4. Every Christoffel symbol with an upper index

Two facts organise the computation.
The metric depends on $r$ and $\theta$ only, so $\partial_t$ and $\partial_\phi$ annihilate every component; and the inverse metric is nonzero only in the slots $(tt)$, $(t\phi)$, $(rr)$, $(\theta\theta)$, $(\phi\phi)$, so the sum over $\alpha$ in the Christoffel formula has at most two terms.

It is easier to lower first.
Write

$$\Gamma_{\mu\nu\rho} = \tfrac{1}{2}\left(\partial_\nu g_{\mu\rho} + \partial_\rho g_{\mu\nu} - \partial_\mu g_{\nu\rho}\right),$$

which needs no inverse metric at all, and then raise with $\Gamma^\mu{}_{\nu\rho} = g^{\mu\alpha}\Gamma_{\alpha\nu\rho}$.
Four worked cases show the pattern and the rest go the same way.

**A symbol with no time index.**
$\Gamma_{r\theta\theta} = \tfrac{1}{2}\left(\partial_\theta g_{r\theta} + \partial_\theta g_{r\theta} - \partial_r g_{\theta\theta}\right) = -\tfrac{1}{2}\partial_r\Sigma = -r$, since $g_{r\theta} = 0$.
Raising costs one factor of $g^{rr}$: $\Gamma^r{}_{\theta\theta} = \frac{\Delta}{\Sigma}\cdot(-r) = -\dfrac{r\Delta}{\Sigma}$.

**A symbol that is a pure derivative.**
$\Gamma_{\theta\theta\theta} = \tfrac{1}{2}\partial_\theta g_{\theta\theta} = \tfrac{1}{2}\partial_\theta\Sigma = -a^2\sin\theta\cos\theta$, and $\Gamma^\theta{}_{\theta\theta} = g^{\theta\theta}\Gamma_{\theta\theta\theta} = -\dfrac{a^2\sin\theta\cos\theta}{\Sigma}$.
Likewise $\Gamma_{rrr} = \tfrac{1}{2}\partial_r g_{rr} = \tfrac{1}{2}\partial_r(\Sigma/\Delta)$, and the quotient rule with the two derivatives of Step 2 gives

$$\Gamma_{rrr} = \frac{2r\Delta - \Sigma(2r-r_s)}{2\Delta^2} = \frac{2a^2r\sin^2\theta - r_s\left(r^2-a^2\cos^2\theta\right)}{2\Delta^2}.$$

The numerator was collected in two pieces.
$2r\left(\Delta - \Sigma\right) = 2r\left(a^2\sin^2\theta - r_sr\right) = 2a^2r\sin^2\theta - 2r^2r_s$, and $r_s\Sigma = r_sr^2 + r_sa^2\cos^2\theta$; adding them leaves $2a^2r\sin^2\theta - r_s\left(r^2 - a^2\cos^2\theta\right)$.
Raising gives $\Gamma^r{}_{rr} = \frac{\Delta}{\Sigma}\Gamma_{rrr} = \dfrac{r}{\Sigma} - \dfrac{2r-r_s}{2\Delta}$, which is the form the entry prints, because each piece of it is a logarithmic derivative of one of the two functions.

**A symbol with one time index.**
$\Gamma_{ttr} = \tfrac{1}{2}\partial_r g_{tt}$, and $g_{tt} = -c^2 + c^2r_sr/\Sigma$, so with $\partial_r(r/\Sigma) = (\Sigma - 2r^2)/\Sigma^2 = -(r^2-a^2\cos^2\theta)/\Sigma^2$,

$$\Gamma_{ttr} = -\frac{c^2r_s\left(r^2-a^2\cos^2\theta\right)}{2\Sigma^2}.$$

This is where the combination $r^2 - a^2\cos^2\theta$ enters, and it never leaves.
Raising needs both terms of the sum, $\Gamma^t{}_{tr} = g^{tt}\Gamma_{ttr} + g^{t\phi}\Gamma_{\phi tr}$, which is why the upper index symbols are longer than the lowered ones.

**A symbol that needs the $\theta$ derivative of the cross term.**
$\Gamma_{tt\theta} = \tfrac{1}{2}\partial_\theta g_{tt} = \tfrac{c^2r_sr}{2}\partial_\theta\Sigma^{-1} = \dfrac{c^2r_sra^2\sin\theta\cos\theta}{\Sigma^2}$, using $\partial_\theta\Sigma^{-1} = -\Sigma^{-2}\partial_\theta\Sigma = 2a^2\sin\theta\cos\theta/\Sigma^2$.

Carrying every case through gives twenty independent symbols, listed here with $\nu \le \rho$ in the printed order; the connection is symmetric in its lower pair, so each with $\nu \ne \rho$ stands for two published components and the entry carries thirty two.

$$\Gamma^t{}_{tr} = \frac{r_s\left(r^2+a^2\right)\left(r^2-a^2\cos^2\theta\right)}{2\Sigma^2\Delta}, \qquad \Gamma^t{}_{t\theta} = -\frac{r_s\,a^2r\sin\theta\cos\theta}{\Sigma^2},$$

$$\Gamma^t{}_{r\phi} = -\frac{r_s\,a\sin^2\theta\left[2r^2\Sigma + \left(r^2+a^2\right)\left(r^2-a^2\cos^2\theta\right)\right]}{2c\,\Sigma^2\Delta}, \qquad \Gamma^t{}_{\theta\phi} = \frac{r_s\,a^3r\sin^3\theta\cos\theta}{c\,\Sigma^2},$$

$$\Gamma^r{}_{tt} = \frac{c^2r_s\Delta\left(r^2-a^2\cos^2\theta\right)}{2\Sigma^3}, \qquad \Gamma^r{}_{t\phi} = -\frac{c\,r_s\,a\Delta\sin^2\theta\left(r^2-a^2\cos^2\theta\right)}{2\Sigma^3},$$

$$\Gamma^r{}_{rr} = \frac{r}{\Sigma} - \frac{2r-r_s}{2\Delta}, \qquad \Gamma^r{}_{r\theta} = -\frac{a^2\sin\theta\cos\theta}{\Sigma}, \qquad \Gamma^r{}_{\theta\theta} = -\frac{r\Delta}{\Sigma},$$

$$\Gamma^r{}_{\phi\phi} = -\frac{\Delta\sin^2\theta\left[2r\Sigma^2 - r_s\,a^2\sin^2\theta\left(r^2-a^2\cos^2\theta\right)\right]}{2\Sigma^3},$$

$$\Gamma^\theta{}_{tt} = -\frac{c^2r_s\,a^2r\sin\theta\cos\theta}{\Sigma^3}, \qquad \Gamma^\theta{}_{t\phi} = \frac{c\,r_s\,ar\left(r^2+a^2\right)\sin\theta\cos\theta}{\Sigma^3},$$

$$\Gamma^\theta{}_{rr} = \frac{a^2\sin\theta\cos\theta}{\Sigma\Delta}, \qquad \Gamma^\theta{}_{r\theta} = \frac{r}{\Sigma}, \qquad \Gamma^\theta{}_{\theta\theta} = -\frac{a^2\sin\theta\cos\theta}{\Sigma},$$

$$\Gamma^\theta{}_{\phi\phi} = -\frac{\sin\theta\cos\theta\left[\left(r^2+a^2\right)\Sigma^2 + r_s\,a^2r\sin^2\theta\left(\Sigma + r^2 + a^2\right)\right]}{\Sigma^3},$$

$$\Gamma^\phi{}_{tr} = \frac{c\,r_s\,a\left(r^2-a^2\cos^2\theta\right)}{2\Sigma^2\Delta}, \qquad \Gamma^\phi{}_{t\theta} = -\frac{c\,r_s\,ar\cos\theta}{\Sigma^2\sin\theta},$$

$$\Gamma^\phi{}_{r\phi} = \frac{r}{\Sigma} - \frac{a^2\sin^2\theta\left[2r\Sigma + r_s\left(r^2-a^2\cos^2\theta\right)\right]}{2\Sigma^2\Delta}, \qquad \Gamma^\phi{}_{\theta\phi} = \cot\theta + \frac{r_s\,a^2r\sin\theta\cos\theta}{\Sigma^2}.$$

Three of these deserve a remark.
$\Gamma^r{}_{t\phi} = -\frac{a\sin^2\theta}{c}\Gamma^r{}_{tt}$ exactly, and $\Gamma^\theta{}_{\phi\phi}$ and $\Gamma^\phi{}_{\theta\phi}$ reduce to their Schwarzschild values $-\sin\theta\cos\theta$ and $\cot\theta$ when $a \to 0$.
And every symbol with $r_s$ in front of it vanishes when $M = 0$, leaving behind exactly the connection of flat space in oblate spheroidal coordinates: $\Gamma^r{}_{rr}$, $\Gamma^r{}_{r\theta}$, $\Gamma^r{}_{\theta\theta}$, $\Gamma^r{}_{\phi\phi}$, $\Gamma^\theta{}_{rr}$, $\Gamma^\theta{}_{r\theta}$, $\Gamma^\theta{}_{\theta\theta}$, $\Gamma^\theta{}_{\phi\phi}$, $\Gamma^\phi{}_{r\phi}$ and $\Gamma^\phi{}_{\theta\phi}$ survive, and they are curvature free.

---

## Step 5. The lowered Christoffel symbols

These are the ones the previous step computed on the way, and they are markedly shorter, because lowering the first index removes the $\Delta$ from most denominators:

$$\Gamma_{ttr} = -\frac{c^2r_s\left(r^2-a^2\cos^2\theta\right)}{2\Sigma^2}, \qquad \Gamma_{tt\theta} = \frac{c^2r_s\,a^2r\sin\theta\cos\theta}{\Sigma^2},$$

$$\Gamma_{tr\phi} = \frac{c\,r_s\,a\sin^2\theta\left(r^2-a^2\cos^2\theta\right)}{2\Sigma^2} = -\Gamma_{rt\phi}, \qquad \Gamma_{t\theta\phi} = -\frac{c\,r_s\,ar\left(r^2+a^2\right)\sin\theta\cos\theta}{\Sigma^2} = \Gamma_{\phi t\theta},$$

$$\Gamma_{rtt} = \frac{c^2r_s\left(r^2-a^2\cos^2\theta\right)}{2\Sigma^2} = -\Gamma_{ttr}, \qquad \Gamma_{\theta tt} = -\frac{c^2r_s\,a^2r\sin\theta\cos\theta}{\Sigma^2} = -\Gamma_{tt\theta},$$

$$\Gamma_{rrr} = \frac{2a^2r\sin^2\theta - r_s\left(r^2-a^2\cos^2\theta\right)}{2\Delta^2}, \qquad \Gamma_{rr\theta} = -\frac{a^2\sin\theta\cos\theta}{\Delta} = -\Gamma_{\theta rr},$$

$$\Gamma_{r\theta\theta} = -r = -\Gamma_{\theta r\theta}, \qquad \Gamma_{\theta\theta\theta} = -a^2\sin\theta\cos\theta,$$

$$\Gamma_{r\phi\phi} = -\frac{\sin^2\theta\left[2r\Sigma^2 - r_s\,a^2\sin^2\theta\left(r^2-a^2\cos^2\theta\right)\right]}{2\Sigma^2} = -\Gamma_{\phi r\phi},$$

$$\Gamma_{\theta\phi\phi} = -\frac{\sin\theta\cos\theta\left[\left(r^2+a^2\right)\Sigma^2 + r_s\,a^2r\sin^2\theta\left(\Sigma + r^2+a^2\right)\right]}{\Sigma^2} = -\Gamma_{\phi\theta\phi},$$

$$\Gamma_{\theta t\phi} = \frac{c\,r_s\,ar\left(r^2+a^2\right)\sin\theta\cos\theta}{\Sigma^2}, \qquad \Gamma_{\phi tr} = \frac{c\,r_s\,a\sin^2\theta\left(r^2-a^2\cos^2\theta\right)}{2\Sigma^2} = \Gamma_{tr\phi}.$$

The antisymmetries visible here are not accidents.
$\Gamma_{\mu\nu\rho}$ is symmetric in $\nu\rho$ by construction, and $\Gamma_{\mu\nu\rho} + \Gamma_{\nu\mu\rho} = \partial_\rho g_{\mu\nu}$; whenever $g_{\mu\nu}$ is independent of $x^\rho$ the two are exact negatives, which is what every equality above with a minus sign in it is saying.
Twenty independent values again, thirty two published components.

---

## Step 6. The Riemann tensor

Feeding the twenty symbols of Step 4 into

$$R^\mu{}_{\nu\rho\sigma} = \partial_\rho \Gamma^\mu_{\nu\sigma} - \partial_\sigma \Gamma^\mu_{\nu\rho} + \Gamma^\mu_{\rho\lambda}\Gamma^\lambda_{\nu\sigma} - \Gamma^\mu_{\sigma\lambda}\Gamma^\lambda_{\nu\rho}$$

and lowering with $R_{\mu\nu\rho\sigma} = g_{\mu\alpha}R^\alpha{}_{\nu\rho\sigma}$ gives thirteen nonzero components.
The count is worth understanding before the values are quoted.
$R_{\mu\nu\rho\sigma}$ is antisymmetric in each pair and symmetric under exchanging the pairs, so it is a symmetric matrix on the six index pairs $tr$, $t\theta$, $t\phi$, $r\theta$, $r\phi$ and $\theta\phi$, which is twenty one entries, reduced to twenty independent ones by the first Bianchi identity.

Eight of the twenty one vanish for Kerr, and they are exactly the eight that carry an odd number of indices from $\{t, \phi\}$.
The reason is a discrete isometry.
The line element of Step 2 is unchanged by $t \to -t$ together with $\phi \to -\phi$, since $g_{t\phi}$ picks up two sign changes and everything else picks up none; reverse time and reverse the rotation and you have the same hole.
Under that map a component of any tensor is multiplied by $(-1)^{n}$ with $n$ the number of its $t$ and $\phi$ indices, so a component with $n$ odd must equal its own negative.
That kills $R_{trt\phi}$, $R_{trr\theta}$, $R_{t\theta t\phi}$, $R_{t\theta r\theta}$, $R_{t\phi r\phi}$, $R_{t\phi\theta\phi}$, $R_{r\theta r\phi}$ and $R_{r\theta\theta\phi}$, and leaves the thirteen below.

The first Bianchi identity then ties three of the survivors together.
With the values of Step 7 in hand, $R_{tr\theta\phi} - R_{t\theta r\phi} + R_{t\phi r\theta} = c\sin\theta\,\mathcal{B}\left[\left(2r^2+2a^2+a^2\sin^2\theta\right) - \left(r^2+a^2+2a^2\sin^2\theta\right) - \Sigma\right] = 0$, since the bracket is $r^2+a^2-a^2\sin^2\theta - \Sigma$.
That is a check on the whole table, because it fails for almost any transcription error.

The next step gives the thirteen in closed form.
It is worth saying first what a direct attack looks like, on the one component that needs nothing but a derivative.
$\Gamma^\theta{}_{r\theta} = r/\Sigma$ and $\Gamma^\theta{}_{\theta\theta} = -a^2\sin\theta\cos\theta/\Sigma$ and $\Gamma^r{}_{\theta\theta} = -r\Delta/\Sigma$ and $\Gamma^r{}_{r\theta} = -a^2\sin\theta\cos\theta/\Sigma$, so

$$R^r{}_{\theta r\theta} = \partial_r\Gamma^r{}_{\theta\theta} - \partial_\theta\Gamma^r{}_{\theta r} + \Gamma^r{}_{r\lambda}\Gamma^\lambda{}_{\theta\theta} - \Gamma^r{}_{\theta\lambda}\Gamma^\lambda{}_{\theta r},$$

whose four pieces are $-\partial_r(r\Delta/\Sigma)$, $+\partial_\theta(a^2\sin\theta\cos\theta/\Sigma)$, $\Gamma^r{}_{rr}\Gamma^r{}_{\theta\theta} + \Gamma^r{}_{r\theta}\Gamma^\theta{}_{\theta\theta}$ and $-\Gamma^r{}_{\theta r}\Gamma^r{}_{\theta r} - \Gamma^r{}_{\theta\theta}\Gamma^\theta{}_{\theta r}$.
Collecting them over $\Sigma^2$, every term without $r_s$ cancels, as it must, because $M=0$ is flat; what is left is

$$R^r{}_{\theta r\theta} = -\frac{r_s r\left(r^2-3a^2\cos^2\theta\right)}{2\Sigma^2}.$$

That numerator, $r(r^2 - 3a^2\cos^2\theta)$, is the first half of the pair that runs through the whole curvature.

---

## Step 7. The two curvature functions

Define

$$\mathcal{A} = \frac{r_s\,r\left(r^2-3a^2\cos^2\theta\right)}{2\Sigma^3}, \qquad \mathcal{B} = \frac{r_s\,a\cos\theta\left(3r^2-a^2\cos^2\theta\right)}{2\Sigma^3}.$$

These two are not independent inventions.
Let $z = r - ia\cos\theta$, so that $z\bar z = r^2 + a^2\cos^2\theta = \Sigma$.
Then

$$\frac{1}{z^3} = \frac{\bar z^3}{(z\bar z)^3} = \frac{(r+ia\cos\theta)^3}{\Sigma^3} = \frac{r\left(r^2-3a^2\cos^2\theta\right) + i\,a\cos\theta\left(3r^2-a^2\cos^2\theta\right)}{\Sigma^3},$$

by expanding the cube and separating real from imaginary, and therefore

$$\mathcal{A} + i\mathcal{B} = \frac{r_s}{2}\frac{1}{\left(r - ia\cos\theta\right)^3} = \frac{GM}{c^2\left(r-ia\cos\theta\right)^3}.$$

One complex cube is the entire curvature of a Kerr black hole.
Setting $a=0$ gives $\mathcal{A} = r_s/(2r^3)$ and $\mathcal{B} = 0$, which is Schwarzschild.

With those two in hand the thirteen independent lowered components are, in the bare chart,

$$R_{trtr} = -c^2\,\frac{2\Delta + a^2\sin^2\theta}{\Delta}\,\mathcal{A}, \qquad R_{trt\theta} = 3ac^2\sin\theta\;\mathcal{B},$$

$$R_{trr\phi} = -c\,\frac{a\sin^2\theta\left[3\left(r^2+a^2\right) - 2r_sr\right]}{\Delta}\,\mathcal{A}, \qquad R_{tr\theta\phi} = c\sin\theta\left[2\left(r^2+a^2\right)+a^2\sin^2\theta\right]\mathcal{B},$$

$$R_{t\theta t\theta} = c^2\left(\Delta + 2a^2\sin^2\theta\right)\mathcal{A}, \qquad R_{t\theta r\phi} = c\sin\theta\left[r^2+a^2+2a^2\sin^2\theta\right]\mathcal{B},$$

$$R_{t\theta\theta\phi} = c\,a\sin^2\theta\left[3\left(r^2+a^2\right)-r_sr\right]\mathcal{A}, \qquad R_{t\phi t\phi} = c^2\Delta\sin^2\theta\;\mathcal{A}, \qquad R_{t\phi r\theta} = -c\,\Sigma\sin\theta\;\mathcal{B},$$

$$R_{r\theta r\theta} = -\frac{\Sigma^2}{\Delta}\,\mathcal{A}, \qquad R_{r\phi r\phi} = -\frac{\sin^2\theta\left[\left(r^2+a^2\right)^2 + 2a^2\Delta\sin^2\theta\right]}{\Delta}\,\mathcal{A},$$

$$R_{r\phi\theta\phi} = 3a\left(r^2+a^2\right)\sin^3\theta\;\mathcal{B}, \qquad R_{\theta\phi\theta\phi} = \sin^2\theta\left[2\left(r^2+a^2\right)^2 + a^2\Delta\sin^2\theta\right]\mathcal{A}.$$

Eight carry $\mathcal{A}$ and five carry $\mathcal{B}$, and no component carries both.
The split is the reflection $\theta \to \pi - \theta$: $\mathcal{A}$ is even in $\cos\theta$ and $\mathcal{B}$ is odd, and the components multiplying $\mathcal{B}$ are exactly those with an odd number of $\theta$ indices.

The power of $c$ on each line is worth reading off now, because it is the whole of Step 13.
$R_{trtr}$, $R_{trt\theta}$, $R_{t\theta t\theta}$ and $R_{t\phi t\phi}$ have two lower $t$ indices and carry $c^2$.
$R_{trr\phi}$, $R_{tr\theta\phi}$, $R_{t\theta r\phi}$, $R_{t\theta\theta\phi}$ and $R_{t\phi r\theta}$ have one and carry $c$.
The last four have none and carry no $c$ at all.

The mixed components $R^\mu{}_{\nu\rho\sigma} = g^{\mu\alpha}R_{\alpha\nu\rho\sigma}$ follow by one contraction each, and there are more of them than of the lowered ones, forty four independent instead of thirteen, because raising an index with a metric whose $(t,\phi)$ block is not diagonal turns one component into a mixture of two.
Two examples show both behaviours.
Raising an index the metric is diagonal in costs a single factor,

$$R^r{}_{\theta r\theta} = g^{rr}R_{r\theta r\theta} = \frac{\Delta}{\Sigma}\cdot\left(-\frac{\Sigma^2}{\Delta}\mathcal{A}\right) = -\Sigma\,\mathcal{A},$$

which is exactly the value Step 6 reached by the long route.
Raising a $t$ index costs two terms, $R^t{}_{\phi t\phi} = g^{tt}R_{t\phi t\phi} + g^{t\phi}R_{\phi\phi t\phi}$, and here the second vanishes because $R_{\phi\phi t\phi} = 0$ by antisymmetry in the first pair.
Where it does not vanish, as in $R^t{}_{tt\phi} = g^{tt}R_{tt t\phi} + g^{t\phi}R_{\phi tt\phi}$, the surviving term carries $g^{t\phi}$, which is itself proportional to $r_s$, so the mixed component comes out quadratic in the mass and the entry prints it with $G^2M^2$ in the numerator.
Four of the forty four are of that kind: $R^t{}_{tt\phi}$, $R^t{}_{tr\theta}$, $R^\phi{}_{\phi t\phi}$ and $R^\phi{}_{\phi r\theta}$, and they are the only curvature components anywhere in the collection that are quadratic in the mass.

---

## Step 8. The Ricci tensor vanishes

Kerr is a vacuum solution, and this is the step that earns the claim.
Contract the first and third indices,

$$R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu} = g^{\alpha\beta}R_{\alpha\mu\beta\nu}.$$

Four of the ten independent slots of a symmetric $R_{\mu\nu}$ are already zero by the parity argument of Step 6, since $R_{tr}$, $R_{t\theta}$, $R_{r\phi}$ and $R_{\theta\phi}$ each carry an odd number of $t$ and $\phi$ indices.
Take the $r\theta$ slot, which is the shortest of the six that are left.
Only the slots where $g^{\alpha\beta}$ is nonzero contribute, and of those, four survive the antisymmetries:

$$R_{r\theta} = g^{tt}R_{tr t\theta} + g^{t\phi}R_{tr\phi\theta} + g^{\phi t}R_{\phi rt\theta} + g^{\phi\phi}R_{\phi r\phi\theta}.$$

Substituting the inverse metric of Step 3 and the Riemann components of Step 7, and pulling out the factor $\dfrac{a\sin\theta}{\Sigma\Delta}\mathcal{B}$ that is common to all four, the bracket that is left is

$$-3\left[\left(r^2+a^2\right)^2 - a^2\Delta\sin^2\theta\right] \; + \; r_sr\left[2\left(r^2+a^2\right)+a^2\sin^2\theta\right] \; + \; r_sr\left[r^2+a^2+2a^2\sin^2\theta\right] \; + \; 3\left(r^2+a^2\right)\left(\Sigma - r_sr\right).$$

The two middle terms add to $3r_sr\left(r^2+a^2\right) + 3r_sra^2\sin^2\theta$, and the first part of that cancels against the $-3r_sr(r^2+a^2)$ hiding in the last term.
What remains is

$$-3\left(r^2+a^2\right)^2 + 3a^2\Delta\sin^2\theta + 3r_sra^2\sin^2\theta + 3\left(r^2+a^2\right)\Sigma,$$

and since $\Sigma = r^2+a^2-a^2\sin^2\theta$, the last term is $3(r^2+a^2)^2 - 3a^2\sin^2\theta(r^2+a^2)$, which kills the first.
That leaves

$$3a^2\sin^2\theta\left[\Delta + r_sr - \left(r^2+a^2\right)\right] = 3a^2\sin^2\theta\left[r^2 - r_sr + a^2 + r_sr - r^2 - a^2\right] = 0,$$

by the definition of $\Delta$ and nothing else.
$R_{r\theta} = 0$.

The other slots go the same way, with more terms and the same ending: every one of them collapses on the single identity $\Delta = r^2 - r_sr + a^2$.
The $\theta\theta$ slot, for instance, needs five terms rather than four, and the two carrying $g^{t\phi}$ are equal to each other and quadratic in $r_s$, so they cancel against the $r_s^2$ part of the $g^{tt}$ and $g^{\phi\phi}$ terms before the linear parts cancel between themselves.
The verification script does all sixteen slots symbolically.

Contracting the other way, on the last lower index, gives $R^\alpha{}_{\mu\nu\alpha} = -R^\alpha{}_{\mu\alpha\nu}$ by the antisymmetry of Riemann in its last pair, so that contraction vanishes too.
This is why the Ricci convention of Step 1 costs nothing here.

Both contractions being zero is what the field equations demand: $R_{\mu\nu} = 0$ is the vacuum equation, with no cosmological constant, and the entry's empty Ricci block is the statement that Kerr solves it.

---

## Step 9. The Ricci scalar and the Einstein tensor

$R = g^{\mu\nu}R_{\mu\nu} = 0$, because every term of the contraction has a vanishing factor.
The entry prints `R = 0`.

$G_{\mu\nu} = R_{\mu\nu} - \tfrac{1}{2}Rg_{\mu\nu} = 0 - 0 = 0$, and the same for $G^\mu{}_\nu$ and $G^{\mu\nu}$, since raising an index on a zero tensor leaves a zero tensor.
Through the field equations $G_{\mu\nu} = 8\pi G T_{\mu\nu}/c^4$ that is the statement that the stress energy tensor is identically zero: there is no matter anywhere in this spacetime, at any radius, and all of the structure the history describes is curvature with nothing to source it locally.
The mass and the angular momentum are properties of the field at infinity rather than of anything sitting at a point.

The three empty blocks in the entry, `ricci_tensor`, `einstein_tensor` and the scalar, are empty for that reason and not for want of working, which is what the entry's `convention` field says in as many words.

---

## Step 10. The Weyl tensor, computed rather than copied

The Weyl tensor is Riemann with its traces removed,

$$C_{\mu\nu\rho\sigma} = R_{\mu\nu\rho\sigma} - \frac{1}{n-2}\left(g_{\mu\rho}R_{\sigma\nu} - g_{\mu\sigma}R_{\rho\nu} - g_{\nu\rho}R_{\sigma\mu} + g_{\nu\sigma}R_{\rho\mu}\right) + \frac{R}{(n-1)(n-2)}\left(g_{\mu\rho}g_{\sigma\nu} - g_{\mu\sigma}g_{\rho\nu}\right),$$

with $n = 4$, and the traces it removes are Riemann's own, $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$, independently of which contraction a file chooses to name Ricci.

Every term in the correction carries a factor of $R_{\mu\nu}$ and every term in the last carries $R$.
Step 8 computed both to be zero.
Therefore

$$C_{\mu\nu\rho\sigma} = R_{\mu\nu\rho\sigma}$$

for Kerr, in every slot, and likewise with the first index raised.

The order of that argument matters and is the reason this step is not a shortcut.
The equality is a consequence of a computed fact, the vanishing of the Ricci tensor, and not a rule about vacuum spacetimes applied on faith.
The verification script builds the Weyl tensor from the definition above, forming the correction and trace terms from its own Riemann tensor and its own inverse metric, and compares the result against the published values; it never reads the published Riemann block when it checks the published Weyl block.
The derivation here was done the same way: the Weyl tensor was computed from the definition in all two hundred and fifty six slots and then found to agree with Riemann in all of them, rather than assumed to.

Kerr is algebraically special, of Petrov type D, which is a statement about this tensor.
Type D means the Weyl tensor has exactly two repeated principal null directions, and in the Newman-Penrose formalism it means that four of the five complex Weyl scalars vanish and the whole tensor is carried by

$$\Psi_2 = -\frac{GM}{c^2\left(r-ia\cos\theta\right)^3} = -\left(\mathcal{A} + i\mathcal{B}\right),$$

which is the cube of Step 7 seen again from the other side.
That is the structural reason there are two curvature functions rather than ten, and it is why the algebra of Step 11 closes in one line.

---

## Step 11. The Kretschmann scalar

$$K = R_{\mu\nu\rho\sigma}R^{\mu\nu\rho\sigma},$$

formed by raising all four indices with the inverse metric and contracting against the lowered tensor.
Doing the contraction on the thirteen independent components of Step 7, with the multiplicities the symmetries give each of them, yields

$$K = 48\left(\mathcal{A}^2 - \mathcal{B}^2\right).$$

That is the whole computation, and the rest is algebra on a complex number.
From Step 7, $\mathcal{A} + i\mathcal{B} = \dfrac{GM}{c^2z^3}$ with $z = r - ia\cos\theta$, so

$$\mathcal{A}^2 - \mathcal{B}^2 = \mathrm{Re}\left[\left(\mathcal{A}+i\mathcal{B}\right)^2\right] = \frac{G^2M^2}{c^4}\,\mathrm{Re}\,\frac{1}{z^6} = \frac{G^2M^2}{c^4}\cdot\frac{\mathrm{Re}\,\bar z^6}{\Sigma^6}.$$

Write $\bar z^2 = (r+ia\cos\theta)^2 = p + iq$ with

$$p = r^2 - a^2\cos^2\theta, \qquad q = 2ra\cos\theta,$$

so that $\bar z^6 = (p+iq)^3$ and $\mathrm{Re}\,\bar z^6 = p^3 - 3pq^2 = p\left(p^2 - 3q^2\right)$.
Now $p^2 = \left(r^2-a^2\cos^2\theta\right)^2 = \Sigma^2 - 4r^2a^2\cos^2\theta$ and $3q^2 = 12r^2a^2\cos^2\theta$, so

$$p^2 - 3q^2 = \Sigma^2 - 16r^2a^2\cos^2\theta,$$

and therefore

$$K = \frac{48G^2M^2}{c^4\,\Sigma^6}\left(r^2-a^2\cos^2\theta\right)\left(\Sigma^2 - 16r^2a^2\cos^2\theta\right),$$

which is what the entry publishes.
A scalar is the same number in every chart, so this expression needs no conversion in Step 13.

Two checks.
Setting $a = 0$ gives $\Sigma = r^2$ and $K = 48G^2M^2r^2\cdot r^4/(c^4r^{12}) = 48G^2M^2/(c^4r^6)$, which is $12r_s^2/r^6$, the Schwarzschild value.
Setting $M = 0$ gives $K = 0$, as it must, since that spacetime is flat.
The dimensions work out as $[G^2M^2/c^4] = L^2$ against $L^2 \cdot L^4 / L^{12} = L^{-6}$, for a total of $L^{-4}$, which is what a squared curvature carries.

---

## Step 12. What the Kretschmann scalar says

$\Delta$ does not appear in $K$.
Neither horizon is therefore a place where anything measurable diverges, and the failure of the Boyer-Lindquist chart at $r = r_\pm$ is a failure of the chart.
Both roots of $\Delta$ are removable, in the sense that a different chart carries the solution across them, which is what Boyer and Lindquist's maximal extension does and what the history records them for.

The only way $K$ can diverge is $\Sigma \to 0$, and $\Sigma = r^2 + a^2\cos^2\theta$ is a sum of two squares.
It vanishes only where both vanish at once,

$$r = 0 \quad\text{and}\quad \cos\theta = 0,$$

which in the flat coordinates the chart reduces to at large distance is the circle of radius $a$ in the equatorial plane.
This is the ring singularity.
It is a genuine curvature singularity, and it is one dimensional rather than pointlike.

The contrast is worth spelling out.
Along the rotation axis, where $\cos\theta = \pm 1$, the surface $r=0$ has $\Sigma = a^2 \ne 0$ and $K = 48G^2M^2(-a^2)(a^4)/(c^4a^{12}) = -48G^2M^2/(c^4a^6)$, a perfectly finite number.
The locus $r=0$ is a disc, not a point, and only its rim is singular; the interior of the disc is ordinary spacetime that a worldline can be continued through.
Everything the Schwarzschild solution does at $r=0$, Kerr does on a circle instead, and the difference is what makes the interior of a rotating hole a qualitatively different object.

$K$ also changes sign, which the Schwarzschild scalar never does.
It vanishes on the cone $r = a\left|\cos\theta\right|$ and on the two surfaces where $\Sigma^2 = 16r^2a^2\cos^2\theta$, so the tidal field of a rotating hole reverses character across those surfaces.
This is the sign of $\mathrm{Re}\,\Psi_2^2$ turning over, and it has no analogue in the spherical case, where $\mathcal{B} = 0$ and $K = 48\mathcal{A}^2 \ge 0$ everywhere.

---

## Step 13. The same components in the chart the entry prints

The entry prints components in the chart $x^0 = ct$.
By Step 1 the rule is: multiply by $c$ once for each upper $t$ index and divide by $c$ once for each lower one.

The metric.
$g_{tt}$ has two lower $t$ indices and becomes $-c^2(1-r_sr/\Sigma)/c^2 = -\left(1-r_sr/\Sigma\right)$.
$g_{t\phi}$ has one and becomes $-r_sra\sin^2\theta/\Sigma$.
$g_{rr}$, $g_{\theta\theta}$ and $g_{\phi\phi}$ have none and are unchanged.
The inverse runs the other way: $g^{tt}$ gains $c^2$ and becomes $-\left[(r^2+a^2)^2 - a^2\Delta\sin^2\theta\right]/(\Sigma\Delta)$, $g^{t\phi}$ gains one $c$ and becomes $-r_sra/(\Sigma\Delta)$, and the rest are unchanged.
Every factor of $c$ in the bare tables of Step 3 is removed exactly, and no published metric component carries one, which is the statement that the chart components of a metric are dimensionless in the sense the checker means.

The connection.
$\Gamma^t{}_{tr}$ and $\Gamma^t{}_{t\theta}$ have one upper and one lower $t$, so the powers cancel and they are unchanged, which matches their having no $c$ in Step 4.
$\Gamma^t{}_{r\phi}$ and $\Gamma^t{}_{\theta\phi}$ have one upper $t$ and gain a $c$, which cancels the $1/c$ they were printed with.
$\Gamma^r{}_{tt}$ and $\Gamma^\theta{}_{tt}$ have two lower $t$ and lose $c^2$, cancelling the $c^2$ they carry.
$\Gamma^r{}_{t\phi}$, $\Gamma^\theta{}_{t\phi}$, $\Gamma^\phi{}_{tr}$ and $\Gamma^\phi{}_{t\theta}$ have one lower $t$ and lose one $c$, cancelling theirs.
The ten symbols with no $t$ index at all are unchanged, and together with $\Gamma^t{}_{tr}$ and $\Gamma^t{}_{t\theta}$ they are the twelve that carry no $c$ in Step 4.
The same bookkeeping applies to the lowered symbols of Step 5, where a symbol with $n$ lower $t$ indices carries $c^n$ and loses all of it.

The curvature.
Every one of the thirteen lines of Step 7 carries $c$ to the power of the number of its lower $t$ indices, so the published lowered Riemann tensor is that table with the $c$ powers struck out, for example

$$R_{t\phi t\phi} = \Delta\sin^2\theta\,\mathcal{A}, \qquad R_{tr\theta\phi} = \sin\theta\left[2\left(r^2+a^2\right)+a^2\sin^2\theta\right]\mathcal{B}, \qquad R_{\theta\phi\theta\phi} = \sin^2\theta\left[2\left(r^2+a^2\right)^2+a^2\Delta\sin^2\theta\right]\mathcal{A}.$$

The mixed components behave the same way, including the four quadratic ones of Step 7: $R^t{}_{tt\phi}$ and $R^\phi{}_{\phi t\phi}$ have one uncancelled lower $t$ and carry the one power of $c$ that removes, while $R^t{}_{tr\theta}$ and $R^\phi{}_{\phi r\theta}$ have none and carry none.

That every factor of $c$ cancels is not a coincidence and it is the useful check.
A component in the chart $x^0 = ct$ is a component in a chart all of whose coordinates are lengths or angles, and $c$ is the only quantity in the problem with a time in it; a chart with no time left in it cannot produce one.
Step 15 makes the same observation quantitatively.

Finally, $r_s$ becomes $2GM/c^2$ everywhere, so that $\mathcal{A}$ and $\mathcal{B}$ are published as

$$\mathcal{A} = \frac{GM\,r\left(r^2-3a^2\cos^2\theta\right)}{c^2\Sigma^3}, \qquad \mathcal{B} = \frac{GM\,a\cos\theta\left(3r^2-a^2\cos^2\theta\right)}{c^2\Sigma^3},$$

and $\Sigma$ and $\Delta$ are written out in full, because the file's reader knows only the symbols the entry declares, and those are $M$, $a$ and $G$.

---

## Step 14. The geodesic equations

The equation of motion is

$$\ddot{x}^\mu + \Gamma^\mu{}_{\nu\rho}\dot{x}^\nu\dot{x}^\rho = 0,$$

with the dot a derivative with respect to an affine parameter $\lambda$ and the components those of the chart the entry prints, so $\dot{t}$ means $d(ct)/d\lambda$ and the $\Gamma$ in it is the published one.
Each symbol with $\nu \ne \rho$ contributes twice, once from each ordering, which is the factor of two in front of the mixed terms below.

$$\ddot{t} + \frac{2GM\left(r^2+a^2\right)\left(r^2-a^2\cos^2\theta\right)}{c^2\Sigma^2\Delta}\dot{t}\dot{r} - \frac{4GMa^2r\sin\theta\cos\theta}{c^2\Sigma^2}\dot{t}\dot{\theta} - \frac{2GMa\sin^2\theta\left[2r^2\Sigma + \left(r^2+a^2\right)\left(r^2-a^2\cos^2\theta\right)\right]}{c^2\Sigma^2\Delta}\dot{r}\dot{\phi} + \frac{4GMa^3r\sin^3\theta\cos\theta}{c^2\Sigma^2}\dot{\theta}\dot{\phi} = 0,$$

$$\ddot{r} + \frac{GM\Delta\left(r^2-a^2\cos^2\theta\right)}{c^2\Sigma^3}\dot{t}^2 - \frac{2GMa\Delta\sin^2\theta\left(r^2-a^2\cos^2\theta\right)}{c^2\Sigma^3}\dot{t}\dot{\phi} + \left(\frac{r}{\Sigma} - \frac{c^2r-GM}{c^2\Delta}\right)\dot{r}^2 - \frac{2a^2\sin\theta\cos\theta}{\Sigma}\dot{r}\dot{\theta} - \frac{r\Delta}{\Sigma}\dot{\theta}^2 - \frac{\Delta\sin^2\theta\left[c^2r\Sigma^2 - GMa^2\sin^2\theta\left(r^2-a^2\cos^2\theta\right)\right]}{c^2\Sigma^3}\dot{\phi}^2 = 0,$$

$$\ddot{\theta} - \frac{2GMa^2r\sin\theta\cos\theta}{c^2\Sigma^3}\dot{t}^2 + \frac{4GMar\left(r^2+a^2\right)\sin\theta\cos\theta}{c^2\Sigma^3}\dot{t}\dot{\phi} + \frac{a^2\sin\theta\cos\theta}{\Sigma\Delta}\dot{r}^2 + \frac{2r}{\Sigma}\dot{r}\dot{\theta} - \frac{a^2\sin\theta\cos\theta}{\Sigma}\dot{\theta}^2 - \frac{\sin\theta\cos\theta\left[c^2\left(r^2+a^2\right)\Sigma^2 + 2GMa^2r\sin^2\theta\left(\Sigma+r^2+a^2\right)\right]}{c^2\Sigma^3}\dot{\phi}^2 = 0,$$

$$\ddot{\phi} + \frac{2GMa\left(r^2-a^2\cos^2\theta\right)}{c^2\Sigma^2\Delta}\dot{t}\dot{r} - \frac{4GMar\cos\theta}{c^2\Sigma^2\sin\theta}\dot{t}\dot{\theta} + 2\left(\frac{r}{\Sigma} - \frac{a^2\sin^2\theta\left[c^2r\Sigma + GM\left(r^2-a^2\cos^2\theta\right)\right]}{c^2\Sigma^2\Delta}\right)\dot{r}\dot{\phi} + \frac{2\cos\theta\left[c^2\Sigma^2 + 2GMa^2r\sin^2\theta\right]}{c^2\Sigma^2\sin\theta}\dot{\theta}\dot{\phi} = 0.$$

The $t$ and $\phi$ equations have no term in $\dot{t}^2$, $\dot{\phi}^2$ or $\dot{t}\dot{\phi}$, which is the Killing structure showing through: the two conserved quantities

$$E \propto -g_{t\mu}\dot{x}^\mu, \qquad L \propto g_{\phi\mu}\dot{x}^\mu$$

are first integrals of exactly those two equations, and up to normalisation they are the energy and the angular momentum of the orbit.
Carter's constant, the third integral the history mentions, is not visible in this form; it comes out of separating the Hamilton-Jacobi equation rather than out of a Killing vector.

The $\phi$ equation shows frame dragging directly.
A particle released from rest at large $r$, with $\dot{r} = \dot\theta = \dot\phi = 0$ and $\dot t \ne 0$, has $\ddot\phi = 0$ at that instant, because every term of the $\phi$ equation carries a velocity other than $\dot t$.
One order further in, its $\dot r$ becomes nonzero through the $r$ equation, and then the $\dot t\dot r$ term of the $\phi$ equation, proportional to $GMa$, spins it up in the direction of the hole's rotation.
Nothing can fall straight in.

---

## Step 15. Every published equation is dimensionally consistent

The chart coordinates are $ct$, $r$, $\theta$ and $\phi$, of which the first two are lengths and the last two are dimensionless.
The parameters are $[M] = M$, $[a] = L$ and $[G] = L^3M^{-1}T^{-2}$, so that $[GM/c^2] = L$ and $[\Sigma] = [\Delta] = L^2$.
A component of a tensor carries, on top of the dimension its rank gives it, a factor $[x^\mu]/L$ for each upper index and $L/[x^\mu]$ for each lower one.

Take the $\phi$ geodesic equation term by term.
Its left hand side is $\ddot\phi$, which carries $1/\lambda^2$ since $\phi$ is an angle.
The first term is $\Gamma^\phi{}_{tr}\dot{t}\dot{r}$ doubled.
$[\Gamma^\phi{}_{tr}] = L^{-1}\cdot(1/L)\cdot 1\cdot 1 = L^{-2}$, where the $1/L$ is the upper $\phi$ and the two lower indices are both lengths in the chart and so contribute nothing; and indeed $\dfrac{2GMa\left(r^2-a^2\cos^2\theta\right)}{c^2\Sigma^2\Delta}$ carries $L\cdot L\cdot L^2/(L^4\cdot L^2) = L^{-2}$.
The velocities carry $[\dot t][\dot r] = (L/\lambda)^2$, so the term carries $L^{-2}\cdot L^2/\lambda^2 = 1/\lambda^2$.
It matches.

The second term, $\Gamma^\phi{}_{t\theta}\dot t\dot\theta$ doubled, has $[\Gamma^\phi{}_{t\theta}] = L^{-1}(1/L)(1)(L) = L^{-1}$, matching $\dfrac{4GMar\cos\theta}{c^2\Sigma^2\sin\theta}$ at $L\cdot L\cdot L/L^4 = L^{-1}$; the velocities carry $(L/\lambda)(1/\lambda)$, so the term is $L^{-1}\cdot L/\lambda^2 = 1/\lambda^2$.
It matches, and it would not have if the dot on $t$ had been read as $dt/d\lambda$ rather than $d(ct)/d\lambda$: that reading makes $[\dot t] = T/\lambda$ and leaves the term short by one factor of $c$, while leaving the first term alone, so no overall rescaling could repair both.
This is the same argument `kasner.md` Step 17 makes, and it is what pins the dot convention down without any algebra.

The third term carries a bracket with two pieces, and the check has to pass on each separately: $r/\Sigma$ is $L^{-1}$ and $\dfrac{a^2\sin^2\theta\left[c^2r\Sigma + GM\left(r^2-a^2\cos^2\theta\right)\right]}{c^2\Sigma^2\Delta}$ has an inner bracket whose two terms are $L^2T^{-2}\cdot L\cdot L^2 = L^5T^{-2}$ and $L^3T^{-2}\cdot L^2 = L^5T^{-2}$, so the whole is $L^2\cdot L^5T^{-2}/(L^2T^{-2}L^4L^2) = L^{-1}$.
Both pieces are $L^{-1}$, the velocities $\dot r\dot\phi$ carry $L/\lambda^2$, and the term is $1/\lambda^2$ again.
That an inner bracket must balance on its own is why the entry writes $c^2r\Sigma + GM(\ldots)$ rather than $r\Sigma + \tfrac{GM}{c^2}(\ldots)$: both are correct, but only the first has every term of every sum carrying the same dimension without a hidden division.

The checker's dimensional pass does this for every term of every published expression in the entry, which is 412 components, four geodesic equations, one line element and one scalar, in well under a second.

---

## Step 16. What the entry publishes, and what the checker needs

The Boyer-Lindquist system carries:

- six metric components and six inverse metric components, the extra pair over the diagonal being $g_{t\phi} = g_{\phi t}$ and its inverse;
- thirty two Christoffel symbols with an upper index and thirty two with all three lowered, from twenty independent values each, which is the longest connection in the collection;
- eighty eight components of $R^\mu{}_{\nu\rho\sigma}$ and eighty of $R_{\mu\nu\rho\sigma}$, from forty four and thirteen independent values;
- the same counts again for the Weyl tensor, which equals Riemann by Step 10;
- empty Ricci and Einstein blocks and `R = 0`, by Step 8 and Step 9;
- the Kretschmann scalar of Step 11;
- four geodesic equations.

Adding the spacetime to the checker took one entry in `DIMENSIONS` in `verify_metrics.py`,

```
("kerr", "boyer_lindquist"): {
    "t": "T", "r": "L", "\\theta": "1", "\\phi": "1",
    "M": "M", "a": "L", "G": "L**3/(M*T**2)",
},
```

which is the documented step for adding a spacetime, and which declares the one thing the script cannot read off the file: that $a$ is a length rather than a dimensionless spin parameter, and that $t$ is the coordinate the chart multiplies by $c$.
Kerr is the second entry after Vaidya to keep $G$ and a mass explicit rather than folding them into a length, and so the second to need a mass in the table at all.
No entry in `PARAMETER_RELATIONS` is needed, because $M$, $a$ and $G$ are free: every published value here is an identity for all three, not an identity on a constraint surface.

The index names in a published component have to be the coordinate names exactly as `coords` spells them, which for this entry means `"\\theta"` and `"\\phi"` rather than `"theta"` and `"phi"`.
That is worth saying because getting it wrong is quiet in one direction and loud in the other: `compare_block` reports every mis-spelled index as a disagreement, but the dimensional pass skips a component whose indices it cannot resolve, so a first draft of this entry passed `--dimensions-only` while three hundred and twenty of its components were not being weighed at all.
A clean dimensional pass on a new entry is only worth something once the sympy comparison has confirmed that the index names resolve.

## Step 17. What the checker managed, and what it could not

Kerr is by a wide margin the slowest system in the collection to check, and at the time of writing the script has never completed a full sympy pass over it.
The reason is `norm`, the routine that puts an expression into a form where a vanishing one is recognisably zero.
It calls sympy's general `simplify`, which is cheap on the diagonal metrics the rest of the collection is made of and very expensive here: the metric is not diagonal, so the Christoffel sums do not collapse, and every curvature component is a rational function of $r$, $\sin\theta$ and $\cos\theta$ with $\Sigma$ to a high power underneath.
The Christoffel symbols come out in about forty seconds, but the Riemann tensor ran for over ninety minutes of processor time without finishing, and the Kretschmann scalar, which raises four indices through a quadruple sum and normalises two hundred and fifty six components on the way, is further out still.
Raising `--budget` does not fix that inside any sitting worth having:

```
/tmp/mfs-venv/bin/python _tools/derivations/verify_metrics.py \
    --system kerr/boyer_lindquist --budget 14400
```

A tensor the script cannot finish inside its budget is reported `UNCHECKED` with the reason rather than passed in silence, so nothing here is being claimed on the script's authority that the script did not establish.
What was established instead, and what a later reader should know stands behind the numbers above:

- The dimensional pass runs on Kerr in well under a second with `--dimensions-only`, and every term of all 412 components, the four geodesic equations, the line element and the Kretschmann scalar balances.
- Every published value was compared against sympy by a separate cross check, which parsed each published string with the script's own `Reader` and compared it against a tensor computed independently from the same line element, with the same chart weighting; it also confirmed that every component the entry does not publish vanishes. All 412 agreed.
- `check_system` itself was then run with only `norm` replaced, by a routine that reduces $\sin\theta$ and $\cos\theta$ to a canonical rational form rather than calling `simplify`. Everything else was the script: its reader, its metric, its Christoffel, Riemann, Ricci, Einstein and Weyl, and its component by component comparison. It reported zero disagreements on the metric, the inverse metric, both Christoffel variants, both Riemann variants, all three Ricci variants, all three Einstein variants, both Weyl variants and the Ricci scalar. It was stopped during the Kretschmann contraction, which the closed form of Step 11 covers by hand.

The honest summary is that the physics is verified and the tool is too slow to say so on its own.
Making `norm` fall back to a canonical trig reduction before reaching for `simplify` would bring Kerr inside the budget and would speed up the whole collection, and that is the obvious next piece of work on the checker.
