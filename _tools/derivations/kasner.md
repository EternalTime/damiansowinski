# The Kasner anisotropic vacuum cosmology

This is the working behind the Cartesian coordinate system in `MFS/assets/data/metrics/kasner.json`.
Every number the entry prints is derived here, in order, from the line element down to the geodesic equations.
Nothing is left as an exercise and nothing is asserted that is not computed.

The companion script `_tools/derivations/verify_metrics.py` does the same work in sympy and compares it against the published file, so the algebra below is checkable by hand and by machine independently.

Kasner is the first entry in the collection whose parameters are not free.
Its three exponents are bound by two algebraic constraints, and everything that makes it a vacuum solution is carried by those constraints rather than by the shape of the metric.
Step 8 derives them from the field equations, and Step 15 gives the parametrisation of the constraint surface that lets a machine check an entry which only claims its values on that surface.

---

## Step 1. Conventions

The signature is $(-,+,+,+)$ and the Christoffel symbols are those of the Levi-Civita connection,

$$\Gamma^\mu_{\nu\rho} = \tfrac{1}{2}g^{\mu\alpha}\left(\partial_\nu g_{\alpha\rho} + \partial_\rho g_{\alpha\nu} - \partial_\alpha g_{\nu\rho}\right),$$

with the lowered symbol $\Gamma_{\mu\nu\rho} = g_{\mu\alpha}\Gamma^\alpha_{\nu\rho}$.
The Riemann tensor is

$$R^\mu{}_{\nu\rho\sigma} = \partial_\rho \Gamma^\mu_{\nu\sigma} - \partial_\sigma \Gamma^\mu_{\nu\rho} + \Gamma^\mu_{\rho\lambda}\Gamma^\lambda_{\nu\sigma} - \Gamma^\mu_{\sigma\lambda}\Gamma^\lambda_{\nu\rho},$$

which is what the published Riemann components are in.

The collection contracts the Ricci tensor on the last lower index,

$$R_{\mu\nu} = R^\alpha{}_{\mu\nu\alpha},$$

which is the opposite sign from the commoner $R^\alpha{}_{\mu\alpha\nu}$.
For this spacetime the choice happens not to matter for anything published, because both contractions vanish together on the constraint surface; Step 8 computes the one the collection uses and Step 9 notes that the other is its negative.

Factors of $c$ are kept explicit.
The chart, here and everywhere else in the collection, is the one whose time coordinate is $x^0 = cT$.
Schwarzschild shows the convention plainly, quoting $g_{tt} = -(1-r_s/r)$ against a line element whose time term is $-(1-r_s/r)c^2dt^2$.
Since $t$ here is a time, the chart coordinate is $ct$, and every component printed against an index $t$ is a component in that chart even though the index is written with the bare letter.

Because the rescaling $t \to ct$ is linear with constant coefficients, the rule is arithmetic: a component in the chart is the component taken with the bare coordinate, multiplied by $c$ once for every upper $t$ index and divided by $c$ once for every lower one.
The Christoffel symbols obey the same rule as the tensors, since the inhomogeneous term in their transformation law carries a second derivative of the coordinate change, which vanishes for a linear one.
A scalar has no index at all and is therefore the same number in both charts, which is worth holding on to for the Kretschmann scalar of Step 12.

Steps 2 to 13 work in the bare $(t,x,y,z)$ chart, where $g_{tt} = -c^2$, because that is where the computation is done.
Step 14 carries the results into the chart the entry prints.

---

## Step 2. The line element

$$ds^2 = -c^2dt^2 + t^{2p_1}dx^2 + t^{2p_2}dy^2 + t^{2p_3}dz^2,$$

with the three constants $p_1, p_2, p_3$ obeying

$$p_1 + p_2 + p_3 = 1, \qquad p_1^2 + p_2^2 + p_3^2 = 1.$$

Step 8 shows that these two conditions are exactly the vacuum field equations for this ansatz, so they are not an extra assumption laid on top of the metric; they are what the metric has to satisfy in order to solve anything.

The chart is $t \in (0,\infty)$ with $x$, $y$ and $z$ each running over the whole real line.
The surface $t=0$ is left out, and Step 13 shows that it is a curvature singularity rather than a defect of the chart.

Two remarks on reading the line element.
First, the exponents are generically irrational, so $t^{2p_i}$ only makes sense with $t$ read in a fixed unit of time; the constant that fixes the unit is absorbed into $x$, $y$ and $z$ once and for all and is not carried below.
Second, the coefficient of $dt^2$ is the only place $c$ appears, which is why the factors of $c$ in later steps all come from $g^{tt} = -1/c^2$.

---

## Step 3. The metric matrix, its determinant and its inverse

In the order $(t,x,y,z)$,

$$g_{\mu\nu} = \begin{pmatrix} -c^2 & 0 & 0 & 0 \\ 0 & t^{2p_1} & 0 & 0 \\ 0 & 0 & t^{2p_2} & 0 \\ 0 & 0 & 0 & t^{2p_3}\end{pmatrix}.$$

The matrix is diagonal, so its determinant is the product of the diagonal,

$$\det g = -c^2\,t^{2p_1}t^{2p_2}t^{2p_3} = -c^2\,t^{2(p_1+p_2+p_3)} = -c^2t^2,$$

where the last equality uses the first constraint of Step 2.
The determinant is negative for every $t>0$, so the signature is Lorentzian on the whole chart, and it is nowhere zero, so the chart is nowhere degenerate.
It vanishes only in the limit $t \to 0$.

Since $\sqrt{-g} = c\,t$, the proper volume of a fixed coordinate box shrinks linearly to zero as $t \to 0$, even though one of the three directions is stretching without bound while the other two collapse.
That is the first sign of the cigar shape the history describes.

A diagonal matrix inverts entry by entry,

$$g^{\mu\nu} = \mathrm{diag}\!\left(-\frac{1}{c^2},\; t^{-2p_1},\; t^{-2p_2},\; t^{-2p_3}\right),$$

and the check $g^{\mu\alpha}g_{\alpha\nu} = \delta^\mu{}_\nu$ is immediate on each of the four diagonal slots, $(-1/c^2)(-c^2) = 1$ and $t^{-2p_i}t^{2p_i} = 1$, with every off diagonal product vanishing because both matrices are diagonal.

---

## Step 4. Every Christoffel symbol

Two facts do all the work.
The metric is diagonal, so $g^{\mu\alpha}$ is nonzero only for $\alpha = \mu$ and the sum in the Christoffel formula collapses to a single term,

$$\Gamma^\mu_{\nu\rho} = \tfrac{1}{2}g^{\mu\mu}\left(\partial_\nu g_{\mu\rho} + \partial_\rho g_{\mu\nu} - \partial_\mu g_{\nu\rho}\right) \quad (\text{no sum on } \mu).$$

And every entry of the metric depends on $t$ alone, so $\partial_x g_{\mu\nu} = \partial_y g_{\mu\nu} = \partial_z g_{\mu\nu} = 0$ and the only surviving derivatives are

$$\partial_t g_{xx} = 2p_1 t^{2p_1-1}, \qquad \partial_t g_{yy} = 2p_2 t^{2p_2-1}, \qquad \partial_t g_{zz} = 2p_3 t^{2p_3-1},$$

with $\partial_t g_{tt} = \partial_t(-c^2) = 0$.

Take the cases in turn, writing $i$ for one of $x,y,z$ and $p_i$ for the matching exponent.

**$\mu = t$.**
$\Gamma^t_{\nu\rho} = \tfrac{1}{2}g^{tt}(\partial_\nu g_{t\rho} + \partial_\rho g_{t\nu} - \partial_t g_{\nu\rho})$.
The first two terms vanish for every $\nu,\rho$, since $g_{t\rho}$ is either $0$ or the constant $-c^2$.
The third vanishes unless $\nu = \rho = i$, so the only nonzero symbols with an upper $t$ are

$$\Gamma^t_{ii} = -\tfrac{1}{2}g^{tt}\partial_t g_{ii} = -\tfrac{1}{2}\left(-\frac{1}{c^2}\right)2p_i t^{2p_i-1} = \frac{p_i\,t^{2p_i-1}}{c^2}.$$

**$\mu = i$.**
$\Gamma^i_{\nu\rho} = \tfrac{1}{2}g^{ii}(\partial_\nu g_{i\rho} + \partial_\rho g_{i\nu} - \partial_i g_{\nu\rho})$.
The last term vanishes always, because nothing depends on a spatial coordinate.
The first term needs $g_{i\rho}$ to be nonzero and $\nu$ dependent, which forces $\rho = i$ and $\nu = t$; the second is the same with $\nu$ and $\rho$ exchanged.
Hence

$$\Gamma^i_{ti} = \Gamma^i_{it} = \tfrac{1}{2}g^{ii}\partial_t g_{ii} = \tfrac{1}{2}t^{-2p_i}\,2p_i t^{2p_i-1} = \frac{p_i}{t},$$

and in particular $\Gamma^i_{tt} = 0$, since $g_{it} = 0$ and $\partial_i g_{tt} = 0$, and $\Gamma^i_{jk} = 0$ for spatial $j,k$ unless one of them is $i$ and the other is $t$, which is not a spatial pair.

Written out, the nonzero symbols are

$$\Gamma^t_{xx} = \frac{p_1 t^{2p_1-1}}{c^2}, \qquad \Gamma^t_{yy} = \frac{p_2 t^{2p_2-1}}{c^2}, \qquad \Gamma^t_{zz} = \frac{p_3 t^{2p_3-1}}{c^2},$$

$$\Gamma^x_{tx} = \Gamma^x_{xt} = \frac{p_1}{t}, \qquad \Gamma^y_{ty} = \Gamma^y_{yt} = \frac{p_2}{t}, \qquad \Gamma^z_{tz} = \Gamma^z_{zt} = \frac{p_3}{t}.$$

That is nine nonzero symbols out of the $4^3 = 64$ index combinations, and the other fifty five vanish by the two cases above.

---

## Step 5. The lowered Christoffel symbols

$\Gamma_{\mu\nu\rho} = g_{\mu\alpha}\Gamma^\alpha_{\nu\rho}$ collapses to $g_{\mu\mu}\Gamma^\mu_{\nu\rho}$ with no sum, again because the metric is diagonal.
So exactly the same nine combinations survive, each multiplied by its own diagonal entry:

$$\Gamma_{tii} = g_{tt}\Gamma^t_{ii} = (-c^2)\frac{p_i t^{2p_i-1}}{c^2} = -p_i\,t^{2p_i-1},$$

$$\Gamma_{iti} = \Gamma_{iit} = g_{ii}\Gamma^i_{ti} = t^{2p_i}\frac{p_i}{t} = p_i\,t^{2p_i-1}.$$

The two differ only in sign, which is the usual statement that $\Gamma_{\mu\nu\rho} = \tfrac{1}{2}(\partial_\nu g_{\mu\rho} + \partial_\rho g_{\mu\nu} - \partial_\mu g_{\nu\rho})$ is antisymmetric in the way the metric derivative is: here $\Gamma_{tii} = -\tfrac{1}{2}\partial_t g_{ii}$ and $\Gamma_{iti} = +\tfrac{1}{2}\partial_t g_{ii}$.

---

## Step 6. The Riemann tensor, first family: the $t$ and $i$ planes

With only nine Christoffel symbols, most of the quadratic terms in the Riemann formula are zero, and the ones that survive can be read off.
Take $\mu = t$, $\nu = \rho = i$ for one spatial $i$, and $\sigma = i$, so that

$$R^t{}_{iti} = \partial_t \Gamma^t_{ii} - \partial_i \Gamma^t_{it} + \Gamma^t_{t\lambda}\Gamma^\lambda_{ii} - \Gamma^t_{i\lambda}\Gamma^\lambda_{it}.$$

The second term vanishes because $\Gamma^t_{it} = 0$ from Step 4.
The third vanishes because $\Gamma^t_{t\lambda} = 0$ for every $\lambda$, again from Step 4.
In the fourth, $\Gamma^t_{i\lambda}$ is nonzero only for $\lambda = i$, and $\Gamma^i_{it} = p_i/t$, so only one term survives.
Therefore

$$R^t{}_{iti} = \partial_t\!\left(\frac{p_i t^{2p_i-1}}{c^2}\right) - \frac{p_i t^{2p_i-1}}{c^2}\cdot\frac{p_i}{t} = \frac{p_i(2p_i-1)t^{2p_i-2}}{c^2} - \frac{p_i^2\,t^{2p_i-2}}{c^2} = \frac{p_i(p_i-1)\,t^{2p_i-2}}{c^2}.$$

The partner with the upper index spatial is

$$R^i{}_{tti} = \partial_t \Gamma^i_{ti} - \partial_i \Gamma^i_{tt} + \Gamma^i_{t\lambda}\Gamma^\lambda_{ti} - \Gamma^i_{i\lambda}\Gamma^\lambda_{tt}.$$

Here $\Gamma^i_{tt} = 0$ kills the second and fourth terms, and in the third only $\lambda = i$ contributes, so

$$R^i{}_{tti} = \partial_t\!\left(\frac{p_i}{t}\right) + \left(\frac{p_i}{t}\right)^2 = -\frac{p_i}{t^2} + \frac{p_i^2}{t^2} = \frac{p_i(p_i-1)}{t^2}.$$

Antisymmetry in the last two indices gives the partners $R^t{}_{iit} = -R^t{}_{iti}$ and $R^i{}_{tit} = -R^i{}_{tti}$.

Lowering the first index,

$$R_{titi} = g_{tt}R^t{}_{iti} = (-c^2)\frac{p_i(p_i-1)t^{2p_i-2}}{c^2} = p_i(1-p_i)\,t^{2p_i-2},$$

$$R_{itit} = g_{ii}R^i{}_{tit} = t^{2p_i}\left(-\frac{p_i(p_i-1)}{t^2}\right) = p_i(1-p_i)\,t^{2p_i-2},$$

which agree, as the pair symmetry $R_{\mu\nu\rho\sigma} = R_{\rho\sigma\mu\nu}$ requires.

---

## Step 7. The Riemann tensor, second family: the spatial planes

Take two distinct spatial directions $i \neq j$ and compute

$$R^j{}_{iij} = \partial_i \Gamma^j_{ij} - \partial_j \Gamma^j_{ii} + \Gamma^j_{i\lambda}\Gamma^\lambda_{ij} - \Gamma^j_{j\lambda}\Gamma^\lambda_{ii}.$$

The first two terms vanish because $\Gamma^j_{ij} = \Gamma^j_{ii} = 0$ for $i \neq j$ spatial.
The third vanishes because $\Gamma^j_{i\lambda} = 0$ for every $\lambda$: a symbol with an upper $j$ needs its lower indices to be $t$ and $j$, and one of them here is $i$.
In the fourth, $\Gamma^j_{j\lambda}$ is nonzero only for $\lambda = t$, and $\Gamma^t_{ii}$ is the surviving partner, so

$$R^j{}_{iij} = -\Gamma^j_{jt}\Gamma^t_{ii} = -\frac{p_j}{t}\cdot\frac{p_i t^{2p_i-1}}{c^2} = -\frac{p_ip_j\,t^{2p_i-2}}{c^2},$$

and therefore

$$R^j{}_{iji} = \frac{p_ip_j\,t^{2p_i-2}}{c^2}.$$

Note which exponent appears: the upper index is $j$ and the exponent is the one belonging to $i$.
Exchanging the two labels gives $R^i{}_{jij} = p_ip_j t^{2p_j-2}/c^2$, which is a genuinely different component, not the same one relabelled.

Lowering the first index makes the asymmetry go away,

$$R_{jiji} = g_{jj}R^j{}_{iji} = t^{2p_j}\frac{p_ip_j t^{2p_i-2}}{c^2} = \frac{p_ip_j\,t^{2p_i+2p_j-2}}{c^2},$$

which is symmetric under exchanging $i$ and $j$, as the pair symmetry demands.

Counting what survives: six coordinate planes, the three $t$ and $i$ planes of Step 6 and the three $i$ and $j$ planes here, each contributing four nonzero components once the antisymmetry in the first pair and in the last pair is written out.
That is twenty four nonzero components out of $4^4 = 256$, which is what the entry lists in each of its two Riemann variants.

---

## Step 8. The Ricci tensor, and where the constraints come from

Contract on the last lower index, $R_{\mu\nu} = R^\alpha{}_{\mu\nu\alpha}$.

For the time slot, only spatial $\alpha$ can contribute, since $R^t{}_{ttt} = 0$ by antisymmetry:

$$R_{tt} = \sum_{i} R^i{}_{tti} = \sum_i \frac{p_i(p_i-1)}{t^2} = \frac{1}{t^2}\left(\sum_i p_i^2 - \sum_i p_i\right),$$

using Step 6.

For a spatial slot, take $i = x$ as the representative.
The contraction runs over $\alpha = t$, which is the Step 6 family, and over the two spatial $\alpha = j \neq x$, which is the Step 7 family:

$$R_{xx} = R^t{}_{xxt} + \sum_{j \neq x} R^j{}_{xxj} = -\frac{p_1(p_1-1)t^{2p_1-2}}{c^2} - \sum_{j\neq x}\frac{p_1p_j\,t^{2p_1-2}}{c^2},$$

where the first term is $-R^t{}_{xtx}$ by antisymmetry.
Factor out $-p_1 t^{2p_1-2}/c^2$:

$$R_{xx} = -\frac{p_1\,t^{2p_1-2}}{c^2}\Big(p_1 - 1 + p_2 + p_3\Big) = -\frac{p_1\,t^{2p_1-2}}{c^2}\left(\sum_i p_i - 1\right),$$

and by the same computation on the other two axes,

$$R_{ii} = -\frac{p_i\,t^{2p_i-2}}{c^2}\left(\sum_i p_i - 1\right).$$

Every off diagonal component vanishes, because a nonzero $R_{\mu\nu}$ would need a nonzero $R^\alpha{}_{\mu\nu\alpha}$, and Steps 6 and 7 leave no Riemann component whose middle two indices are different.

Now impose $R_{\mu\nu} = 0$, which for a vacuum with no cosmological constant is the whole of the field equations.
A spatial equation reads $p_i\left(\sum_i p_i - 1\right)t^{2p_i-2} = 0$, and $t^{2p_i-2}$ is nowhere zero on $t>0$, so each axis demands $p_i = 0$ or $\sum_i p_i = 1$.
If every exponent vanished the line element would be $-c^2dt^2 + dx^2 + dy^2 + dz^2$, which is Minkowski spacetime written out and has nothing left to solve, so at least one exponent is nonzero and all three spatial equations are satisfied for all $t$ if and only if

$$\sum_i p_i = 1.$$

Given that, the time equation reads

$$R_{tt} = \frac{1}{t^2}\left(\sum_i p_i^2 - 1\right) = 0 \iff \sum_i p_i^2 = 1.$$

So the two constraints the entry publishes are precisely the vacuum field equations for a diagonal power law metric, one from the spatial trace and one from the time slot.
On the constraint surface,

$$R_{\mu\nu} = 0,$$

and the entry publishes all three variants $R_{\mu\nu}$, $R^\mu{}_\nu$ and $R^{\mu\nu}$ as empty, since raising an index on a vanishing tensor leaves it vanishing.

---

## Step 9. The Ricci scalar and the Einstein tensor

$$R = g^{\mu\nu}R_{\mu\nu} = 0,$$

a sum of multiples of components that vanish by Step 8, and then

$$G_{\mu\nu} = R_{\mu\nu} - \tfrac{1}{2}Rg_{\mu\nu} = 0 - 0 = 0,$$

in every variant, for the same reason.
The entry publishes $R = 0$ and three empty lists for the Einstein tensor.

This is where the choice of contraction stops mattering.
The other convention gives $R^\alpha{}_{\mu\alpha\nu} = -R^\alpha{}_{\mu\nu\alpha}$, by the antisymmetry of Riemann in its last two indices, so the two Ricci tensors differ by an overall sign and vanish together.
A reader who prefers the textbook contraction reads exactly the same zeros here.

---

## Step 10. The Weyl tensor

In four dimensions,

$$C_{\mu\nu\rho\sigma} = R_{\mu\nu\rho\sigma} - \frac{1}{2}\left(g_{\mu\rho}\mathcal{R}_{\sigma\nu} - g_{\mu\sigma}\mathcal{R}_{\rho\nu} - g_{\nu\rho}\mathcal{R}_{\sigma\mu} + g_{\nu\sigma}\mathcal{R}_{\rho\mu}\right) + \frac{\mathcal{R}}{6}\left(g_{\mu\rho}g_{\sigma\nu} - g_{\mu\sigma}g_{\rho\nu}\right),$$

where $\mathcal{R}_{\mu\nu}$ is the trace of Riemann and $\mathcal{R}$ its scalar.
By Step 8 both are zero on the constraint surface, whichever of the two contractions is used, so every correction term is zero and

$$C_{\mu\nu\rho\sigma} = R_{\mu\nu\rho\sigma}.$$

The entry publishes the Weyl tensor with exactly the component list it publishes for Riemann, in both the $C^\mu{}_{\nu\rho\sigma}$ and the $C_{\mu\nu\rho\sigma}$ variant.
Physically this says that all the curvature of a Kasner universe is tidal.
There is no matter anywhere, and what survives at $t \to 0$ is pure shear.

---

## Step 11. Every nonzero curvature component, collected

The two families of Steps 6 and 7, with $i \neq j$ running over $x,y,z$ and $p_i$ the matching exponent:

$$R^t{}_{iti} = -R^t{}_{iit} = \frac{p_i(p_i-1)t^{2p_i-2}}{c^2}, \qquad R^i{}_{tti} = -R^i{}_{tit} = \frac{p_i(p_i-1)}{t^2},$$

$$R^j{}_{iji} = -R^j{}_{iij} = \frac{p_ip_j\,t^{2p_i-2}}{c^2},$$

$$R_{titi} = R_{itit} = -R_{tiit} = -R_{itti} = p_i(1-p_i)\,t^{2p_i-2},$$

$$R_{ijij} = R_{jiji} = -R_{ijji} = -R_{jiij} = \frac{p_ip_j\,t^{2p_i+2p_j-2}}{c^2}.$$

Every one of these is a component in the bare chart.
Step 14 carries them into the chart the entry prints, which changes only the first family.

---

## Step 12. The Kretschmann scalar

$K = R_{\mu\nu\rho\sigma}R^{\mu\nu\rho\sigma}$.
The metric is diagonal, so raising all four indices on a component whose index pattern is $(a,b,a,b)$ multiplies it by $g^{aa}g^{bb}g^{aa}g^{bb}$, with no sum, and

$$R_{abab}R^{abab} = \left(g^{aa}g^{bb}R_{abab}\right)^2.$$

Step 11 leaves four nonzero components for each of the six coordinate planes, namely $abab$, $abba$, $baab$ and $baba$, and each of the four contributes the same square.
Hence

$$K = 4\sum_{a<b}\left(g^{aa}g^{bb}R_{abab}\right)^2,$$

a sum of six terms.
Take them in the two families.

For a $t$ and $i$ plane,

$$g^{tt}g^{ii}R_{titi} = \left(-\frac{1}{c^2}\right)t^{-2p_i}\,p_i(1-p_i)t^{2p_i-2} = \frac{p_i(p_i-1)}{c^2t^2},$$

where the powers of $t$ cancel completely.
For an $i$ and $j$ plane,

$$g^{ii}g^{jj}R_{ijij} = t^{-2p_i}t^{-2p_j}\frac{p_ip_j\,t^{2p_i+2p_j-2}}{c^2} = \frac{p_ip_j}{c^2t^2},$$

where again they cancel completely.
Both families therefore carry the same $1/(c^2t^2)$, and

$$K = \frac{4}{c^4t^4}\left(\sum_i p_i^2(p_i-1)^2 + \sum_{i<j}p_i^2p_j^2\right).$$

That is the general answer, true of any diagonal power law metric of this shape.
The constraints now collapse it.

Write $e_1 = \sum_i p_i$, $e_2 = \sum_{i<j}p_ip_j$ and $e_3 = p_1p_2p_3$ for the elementary symmetric functions, and $S_k = \sum_i p_i^k$ for the power sums.
The first constraint is $e_1 = 1$.
The second is $S_2 = 1$, and since $S_2 = e_1^2 - 2e_2 = 1 - 2e_2$, it says

$$e_2 = 0.$$

The two constraints together are therefore $e_1 = 1$, $e_2 = 0$, with $e_3$ left free, which is the one free parameter of the family.
Newton's identities with these values give the higher power sums:

$$S_3 = e_1S_2 - e_2S_1 + 3e_3 = 1 - 0 + 3e_3 = 1 + 3e_3,$$

$$S_4 = e_1S_3 - e_2S_2 + e_3S_1 = (1+3e_3) - 0 + e_3 = 1 + 4e_3.$$

The first sum in $K$ is

$$\sum_i p_i^2(p_i-1)^2 = \sum_i \left(p_i^2 - p_i\right)^2 = S_4 - 2S_3 + S_2 = (1+4e_3) - 2(1+3e_3) + 1 = -2e_3.$$

The second is

$$\sum_{i<j}p_i^2p_j^2 = e_2^2 - 2e_1e_3 = 0 - 2e_3 = -2e_3.$$

The two are equal, and

$$K = \frac{4}{c^4t^4}\left(-4e_3\right) = -\frac{16\,p_1p_2p_3}{c^4t^4}.$$

Being a scalar, this is the same number in the bare chart and in the chart the entry prints, which is why it carries $c^{-4}$ in both.

---

## Step 13. What the Kretschmann scalar says

First, the sign.
On the constraint surface $e_2 = 0$, so if all three exponents were positive then $e_2$ would be a sum of three positive terms and could not vanish; at least one exponent is therefore not positive.
Two cannot be negative at once: if $p_1$ and $p_2$ were both negative then $p_3 = 1 - p_1 - p_2 > 1$ and $S_2 \geq p_3^2 > 1$, contradicting the second constraint.
If one exponent is zero, say $p_3 = 0$, then $p_1 + p_2 = 1$ and $p_1^2 + p_2^2 = 1$ force $2p_1p_2 = 0$, so the triple is $(1,0,0)$ up to order.

So every Kasner solution has exactly one negative exponent and two positive ones, except for the three permutations of $(1,0,0)$.
In the first case $p_1p_2p_3 < 0$ and $K > 0$; in the second $p_1p_2p_3 = 0$ and $K = 0$.
The published $K = -16p_1p_2p_3/(c^4t^4)$ is therefore never negative, which it had better not be for a spacetime whose Riemann tensor is entirely tidal.

Second, the exceptional points.
At $(1,0,0)$ every Riemann component of Step 11 vanishes term by term: the first family carries $p_i(p_i-1)$, which is $1\cdot 0$ for the exponent $1$ and $0\cdot(-1)$ for each exponent $0$; the second carries $p_ip_j$, and every pair of the three exponents contains a zero.
So those three solutions are flat, and the apparent expansion in $ds^2 = -c^2dt^2 + t^2dx^2 + dy^2 + dz^2$ is the Milne chart of a piece of Minkowski spacetime, not curvature.
They are the three points where the Kasner circle meets flatness, and $K=0$ detects them exactly.

Third, the singularity.
For every other solution $K$ grows as $t^{-4}$ and diverges as $t \to 0$.
A scalar cannot be removed by a change of coordinates, so $t=0$ is a curvature singularity and not a defect of the chart, which is the point of computing $K$ at all rather than reading the divergence off $g_{\mu\nu}$.
It is also reached: a comoving observer at fixed $x,y,z$ has $ds^2 = -c^2dt^2$, so the proper time from $t_1 > 0$ back to the singularity is $\int_0^{t_1}dt = t_1$, which is finite.

---

## Step 14. The same components in the chart the entry prints

Steps 2 to 13 worked with the bare $t$, which is a time.
The entry prints the chart of Step 1, whose time coordinate is $ct$, so each upper $t$ index multiplies by $c$ and each lower one divides by it.
Applying that count to each published value:

$$g_{tt} = \frac{1}{c^2}\left(-c^2\right) = -1, \qquad g^{tt} = c^2\left(-\frac{1}{c^2}\right) = -1,$$

with $g_{ii} = t^{2p_i}$ and $g^{ii} = t^{-2p_i}$ unchanged, since they carry no time index.
The product is still the identity, since each factor of $c$ introduced on the inverse cancels the one removed from the metric.

The Christoffel symbols of Steps 4 and 5 each carry exactly one time index, upper for the first group and lower for the second, so every one of them changes by a single factor of $c$ in the same direction:

$$\Gamma^t_{ii} = c\cdot\frac{p_i t^{2p_i-1}}{c^2} = \frac{p_i\,t^{2p_i-1}}{c}, \qquad \Gamma^i_{ti} = \Gamma^i_{it} = \frac{1}{c}\cdot\frac{p_i}{t} = \frac{p_i}{ct},$$

$$\Gamma_{tii} = \frac{1}{c}\left(-p_i\,t^{2p_i-1}\right) = -\frac{p_i\,t^{2p_i-1}}{c}, \qquad \Gamma_{iti} = \Gamma_{iit} = \frac{1}{c}\,p_i\,t^{2p_i-1} = \frac{p_i\,t^{2p_i-1}}{c}.$$

The first of these gains a factor of $c$ because its time index is upper; the other three lose one because theirs is lower.
Note that the bare $\Gamma^t_{ii}$ already carried $1/c^2$ from $g^{tt}$, so what the entry prints is $1/c$ rather than no $c$ at all.

For the Riemann tensor, a component with one upper $t$ and one lower $t$ is unchanged, since the two factors cancel, and the same for a component with none.
$R^t{}_{iti}$ has one of each, and $R^j{}_{iji}$ has neither, so both are printed exactly as Step 11 computed them.
$R^i{}_{tti}$ has two lower $t$ and no upper, so it divides by $c^2$:

$$R^i{}_{tti} = \frac{p_i(p_i-1)}{c^2t^2}.$$

$R_{titi}$ has two lower $t$ and divides by $c^2$ as well,

$$R_{titi} = \frac{p_i(1-p_i)t^{2p_i-2}}{c^2},$$

while $R_{ijij}$ has no time index and is unchanged.
The Ricci tensor, the Einstein tensor and the Ricci scalar are zero and stay zero.
The Kretschmann scalar has no index and is the same number in both charts.

These are the values the entry publishes.

---

## Step 15. The Kasner circle, and how a machine checks a constrained entry

The two constraints cut a curve out of the space of exponents.
The plane $p_1+p_2+p_3=1$ meets the unit sphere $p_1^2+p_2^2+p_3^2=1$ in a circle, centred on $(1/3,1/3,1/3)$, and every Kasner universe is one point of it.

Fixing one exponent pins the other two.
Given $p_1$, the constraints say $p_2 + p_3 = 1 - p_1$ and

$$2p_2p_3 = (p_2+p_3)^2 - (p_2^2+p_3^2) = (1-p_1)^2 - (1-p_1^2) = 2p_1^2 - 2p_1,$$

so $p_2$ and $p_3$ are the two roots of $\zeta^2 - (1-p_1)\zeta + p_1(p_1-1) = 0$.
Those roots are real exactly when the discriminant

$$(1-p_1)^2 - 4p_1(p_1-1) = (1-p_1)(1+3p_1)$$

is not negative, which bounds every exponent into

$$-\tfrac{1}{3} \leq p_i \leq 1.$$

The endpoints are $(-1/3,2/3,2/3)$ and $(1,0,0)$, the two most symmetric points on the circle.

A verifier that treated $p_1,p_2,p_3$ as free symbols would be checking a claim the entry never makes.
Step 8 is the proof: for free exponents the Ricci tensor is not zero, so the published empty lists would be read as errors, and the Kretschmann scalar of Step 12 would keep its general quartic form instead of collapsing to $-16p_1p_2p_3$.
What has to be checked is an identity along the circle, not an identity in three unconstrained variables.

The circle is a conic, so it has a rational parametrisation, and the one the collection uses is the map that Belinskii, Khalatnikov and Lifshitz bounce along:

$$p_1 = \frac{-u}{1+u+u^2}, \qquad p_2 = \frac{1+u}{1+u+u^2}, \qquad p_3 = \frac{u(1+u)}{1+u+u^2}.$$

It satisfies both constraints identically.
For the first, the numerators add to $-u + (1+u) + (u+u^2) = 1+u+u^2$, cancelling the denominator.
For the second, the squared numerators add to

$$u^2 + (1+u)^2 + u^2(1+u)^2 = u^2 + 1 + 2u + u^2 + u^2 + 2u^3 + u^4 = u^4 + 2u^3 + 3u^2 + 2u + 1 = (1+u+u^2)^2,$$

which cancels the squared denominator.

The map is onto: a rational parametrisation of an irreducible conic covers all of it but a single point, here the limit $u \to \infty$ at $(0,0,1)$.
That is what makes the check exact rather than a sample.
An expression that vanishes for every $u$ vanishes on a dense subset of the circle, and a rational function vanishing on a dense subset of an irreducible curve vanishes on the whole of it.

`PARAMETER_RELATIONS` in `verify_metrics.py` carries this map.
The script computes every tensor from the line element with the exponents left free, which is fast and keeps the algebra readable, and substitutes the parametrisation into both sides of each comparison, so what it compares are two functions of $u$.
On $u \geq 1$ the parametrisation gives $p_1 \leq p_2 \leq p_3$, with $p_1 \in [-1/3,0]$, $p_2 \in [0,2/3]$ and $p_3 \in [2/3,1]$; the rest of the circle is the same points with the axes relabelled, reached at other $u$.
Setting $u=-1$ gives $(1,0,0)$, one of the three flat solutions of Step 13, and the Kretschmann scalar

$$K = \frac{16\,u^2(1+u)^2}{c^4t^4\left(1+u+u^2\right)^3}$$

duly vanishes there.

---

## Step 16. The geodesic equations

The geodesic equation is $\ddot x^\mu + \Gamma^\mu_{\nu\rho}\dot x^\nu \dot x^\rho = 0$, where the dot is $d/d\lambda$ for an affine parameter $\lambda$ and the symbols are the chart ones of Step 14.

The velocities have to be the velocities of the same chart those symbols are printed in, and that is worth stating in full, because it is the one place in the entry where the convention of Step 1 is easy to lose.
The chart coordinates are $x^\mu = (ct, x, y, z)$, so

$$\dot t \equiv \frac{d(ct)}{d\lambda} = c\frac{dt}{d\lambda}, \qquad \ddot t \equiv \frac{d^2(ct)}{d\lambda^2},$$

even though the index and the dot are both printed with the bare letter $t$.
Every velocity is then a length over $\lambda$ and every acceleration a length over $\lambda^2$, the time slot included, and Step 17 checks that the published equations balance on that reading and on no other.

For $\mu = t$, the only nonzero symbols are the three $\Gamma^t_{ii}$, each appearing once:

$$\ddot t + \frac{1}{c}\left(p_1t^{2p_1-1}\dot x^2 + p_2t^{2p_2-1}\dot y^2 + p_3t^{2p_3-1}\dot z^2\right) = 0.$$

For $\mu = i$, the only nonzero symbols are $\Gamma^i_{ti} = \Gamma^i_{it} = p_i/(ct)$, and the double sum picks each of them up once, giving a factor of two:

$$\ddot x + \frac{2p_1}{ct}\dot t\dot x = 0, \qquad \ddot y + \frac{2p_2}{ct}\dot t\dot y = 0, \qquad \ddot z + \frac{2p_3}{ct}\dot t\dot z = 0,$$

which is what the entry prints.

The three spatial equations integrate once.
Using $\dot t = c\,dt/d\lambda$ to clear the $c$, the $x$ equation reads

$$\frac{d^2x}{d\lambda^2} + \frac{2p_1}{t}\frac{dt}{d\lambda}\frac{dx}{d\lambda} = 0,$$

and multiplying by $t^{2p_1}$ makes the left side a total derivative,

$$\frac{d}{d\lambda}\left(t^{2p_1}\frac{dx}{d\lambda}\right) = t^{2p_1}\frac{d^2x}{d\lambda^2} + 2p_1t^{2p_1-1}\frac{dt}{d\lambda}\frac{dx}{d\lambda} = t^{2p_1}\left(\frac{d^2x}{d\lambda^2} + \frac{2p_1}{t}\frac{dt}{d\lambda}\frac{dx}{d\lambda}\right) = 0.$$

So $t^{2p_1}dx/d\lambda$ is constant along a geodesic, which is $g_{xx}\dot x$, the conserved momentum belonging to the Killing vector $\partial_x$; it exists because nothing in the metric depends on $x$, and likewise for $y$ and $z$.
A particle moving along a contracting axis is therefore blue shifted without bound as $t \to 0$, since $dx/d\lambda = \text{const}\cdot t^{-2p_1}$ with $p_1 < 0$ blows up, while one moving along a stretching axis is slowed.

---

## Step 17. Every published equation is dimensionally consistent

A convention that is only in the reader's head is a convention that gets lost, and the way it shows is a term that cannot be added to the term beside it.
So here is the check, run over every equation the entry publishes.

Two things have to be fixed first.
The powers $t^{2p_i}$ are dimensionless, by the second remark of Step 2, so an integer shift in the exponent is what carries dimension: $t^{2p_i-1}$ is $t^{2p_i}/t$ and goes as $1/T$, and $t^{2p_i+2p_j-2}$ goes as $1/T^2$.
And the dots are the chart velocities of Step 16, so each is $L/\lambda$ and each double dot is $L/\lambda^2$.

With $L$ for length and $T$ for time, the chart coordinates $x^\mu = (ct,x,y,z)$ are all lengths, which fixes what everything else has to be.

| quantity | must be | published form | check |
| --- | --- | --- | --- |
| $ds^2$ | $L^2$ | $-c^2dt^2$ | $(L/T)^2T^2 = L^2$ |
| | | $t^{2p_i}dx^2$ | $1 \cdot L^2 = L^2$ |
| $g_{\mu\nu}$, $g^{\mu\nu}$ | $1$ | $-1$, $t^{\pm 2p_i}$ | dimensionless |
| $\Gamma^\mu{}_{\nu\rho}$ | $1/L$ | $p_i t^{2p_i-1}/c$ | $(1/T)(T/L) = 1/L$ |
| | | $p_i/(ct)$ | $(T/L)(1/T) = 1/L$ |
| $\Gamma_{\mu\nu\rho}$ | $1/L$ | $\pm p_i t^{2p_i-1}/c$ | $1/L$ |
| $R^\mu{}_{\nu\rho\sigma}$, $R_{\mu\nu\rho\sigma}$, $C$ | $1/L^2$ | $p_i(p_i-1)t^{2p_i-2}/c^2$ | $(1/T^2)(T^2/L^2) = 1/L^2$ |
| | | $p_i(p_i-1)/(c^2t^2)$ | $(T^2/L^2)(1/T^2) = 1/L^2$ |
| | | $p_ip_j t^{2p_i+2p_j-2}/c^2$ | $(1/T^2)(T^2/L^2) = 1/L^2$ |
| $K$ | $1/L^4$ | $-16p_1p_2p_3/(c^4t^4)$ | $(T^4/L^4)(1/T^4) = 1/L^4$ |
| geodesic term | $L/\lambda^2$ | $\ddot t$, $\ddot x$ | $L/\lambda^2$ |
| | | $p_i t^{2p_i-1}\dot x^2/c$ | $(1/T)(L^2/\lambda^2)(T/L) = L/\lambda^2$ |
| | | $2p_i\dot t\dot x/(ct)$ | $(T/L)(1/T)(L/\lambda)^2 = L/\lambda^2$ |

The Christoffel symbols are the row to look at twice.
Both of the printed forms come out as $1/L$, and they have to, because a connection coefficient in a chart of lengths is one over a length whatever its index pattern is.
The two get there differently: $\Gamma^t{}_{ii}$ carried $1/c^2$ in the bare chart and gained a factor of $c$ in Step 14, while $\Gamma^i{}_{ti}$ carried none and lost one.

The last row is the one the convention lives or dies on.
Read $\dot t$ as $dt/d\lambda$ instead of $d(ct)/d\lambda$ and the term $2p_i\dot t\dot x/(ct)$ becomes $(T/L)(1/T)(T/\lambda)(L/\lambda) = T/\lambda^2$, which cannot be added to the $\ddot x$ beside it, a length over $\lambda^2$.
The equation would be out by exactly one factor of $c$.
That is not a hypothetical: it is how a wrong reading of this convention shows itself, and a geodesic equation is where it shows, because a geodesic equation is the only thing the entry publishes that adds a time derivative to a space derivative.

All 137 expressions and terms the entry publishes were put through this check mechanically, parsed out of the JSON rather than copied by hand, and every one balances.
Repeating it with $\dot t$ read as $dt/d\lambda$ fails on exactly three of them, the three spatial geodesic equations, which is the signature of the mistake.

---

## Step 18. What the entry publishes

| field | value |
| --- | --- |
| `metric_components` | $g_{tt} = -1$, $g_{ii} = t^{2p_i}$ |
| `inverse_metric_components` | $g^{tt} = -1$, $g^{ii} = t^{-2p_i}$ |
| `christoffel` `ull` | nine symbols, Step 14 |
| `christoffel` `lll` | nine symbols, Step 14 |
| `riemann` `ulll` | twenty four components, Steps 11 and 14 |
| `riemann` `llll` | twenty four components, Steps 11 and 14 |
| `ricci_tensor` | empty in all three variants, Step 8 |
| `ricci_scalar` | $R = 0$, Step 9 |
| `kretschmann` | $K = -16p_1p_2p_3/(c^4t^4)$, Step 12 |
| `einstein_tensor` | empty in all three variants, Step 9 |
| `weyl_tensor` | the Riemann lists, Step 10 |
| `geodesics` | four equations, Step 16, in the chart velocities of that step |

To check it:

    python3 -m venv /tmp/mfs-venv && /tmp/mfs-venv/bin/pip install sympy
    /tmp/mfs-venv/bin/python _tools/derivations/verify_metrics.py --system kasner/cartesian

which takes about five seconds and reports no disagreement.
