# The Kerr-Newman charged rotating black hole

This is the working behind the Boyer-Lindquist coordinate system in `MFS/assets/data/metrics/kerr_newman.json`.
Every number the entry prints is derived here, in order, from the line element down to the geodesic equations.
Nothing is left as an exercise and nothing is asserted that is not computed.

The companion script `_tools/derivations/verify_metrics.py` does the same work in sympy and compares it against the published file, so the algebra below is checkable by hand and by machine independently.

This entry is written to be read beside `_tools/derivations/kerr.md`, in the same chart and the same shapes, because the whole of the difference between the two spacetimes is one length.
Three facts organise the comparison, and each of them is established below rather than assumed.

The first is that the metric is Kerr's with the mass term $r_s r$ replaced by the mass function $\mathcal{M} = r_s r - r_Q^2$, which is Step 2.
The second is that this replacement commutes with both derivatives the connection needs, because $\partial_r\mathcal{M} = r_s$ is what $\partial_r(r_sr)$ was and $\partial_\theta\mathcal{M} = 0$ is what it was, so the whole connection is Kerr's with the same replacement and nothing else, which is Step 4.
The third is that the curvature is where the two solutions genuinely part company: Kerr's Ricci tensor vanishes and Kerr-Newman's does not, because a charged hole carries an electromagnetic field everywhere outside itself and that field has energy.
Steps 6 to 8 compute the field and its stress, Step 11 separates the Weyl tensor from it, and Step 12 shows the Kretschmann scalar picking up exactly two new terms.

Setting $r_Q = 0$ in any expression below returns the corresponding Kerr expression, and the places where that is worth saying out loud are said.

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
Unlike Kerr, this spacetime is one where the choice is visible in the published file, because the Ricci tensor here is not zero.
Step 7 computes it on this convention and Step 8 carries it into the Einstein tensor, where it comes out with a positive energy density for the electromagnetic field, which is the check that the convention and the signature agree with each other.

Factors of $c$ and $G$ are kept explicit.
Kerr-Newman keeps a mass as a mass, as Kerr does, so $M$ appears with Newton's constant beside it and the Schwarzschild radius of that mass, $r_s = 2GM/c^2$, is spelled out as $2GM/c^2$ in every published value.
The charge is the one quantity that is not kept in its own units, because the collection's dimensional table is written in $L$, $T$ and $M$ and there is no base dimension in it for charge.
It is folded instead into the charge radius $r_Q$, exactly as the Reissner-Nordstrom entry folds it, by

$$r_Q^2 = \frac{GQ^2}{4\pi\epsilon_0c^4},$$

which is a length squared: $Q^2/(4\pi\epsilon_0)$ is an energy times a length, $G/c^4$ is a length over an energy, and the product is an area.
Steps 2 to 13 use the shorthand $r_s$ and keep $r_Q$ as it stands; Step 14 puts $2GM/c^2$ back in place of $r_s$.

The chart, here and everywhere else in the collection, is the one whose time coordinate is $x^0 = cT$.
Since $t$ here is a time, the chart coordinate is $ct$, and every component printed against an index $t$ is a component in that chart even though the index is written with the bare letter.
Because the rescaling $t \to ct$ is linear with constant coefficients, the rule is arithmetic: a component in the chart is the component taken with the bare coordinate, multiplied by $c$ once for every upper $t$ index and divided by $c$ once for every lower one.
The Christoffel symbols obey the same rule as the tensors, since the inhomogeneous term in their transformation law carries a second derivative of the coordinate change, which vanishes for a linear one.
A scalar has no index at all and is therefore the same number in both charts.

Steps 2 to 13 work in the bare $(t,r,\theta,\phi)$ chart, where $g_{tt}$ carries a $c^2$, because that is where the computation is done.
Step 14 carries the results into the chart the entry prints, and the arithmetic there is the check that no factor of $c$ was dropped.

---

## Step 2. The line element, and the four blocks the charge enters through

$$ds^2 = -\left(1 - \frac{\mathcal{M}}{\Sigma}\right)c^2dt^2 - \frac{2\mathcal{M}a\sin^2\theta}{\Sigma}c\,dt\,d\phi + \frac{\Sigma}{\Delta}dr^2 + \Sigma\,d\theta^2 + \left(r^2 + a^2 + \frac{\mathcal{M}a^2\sin^2\theta}{\Sigma}\right)\sin^2\theta\,d\phi^2,$$

with

$$\Sigma = r^2 + a^2\cos^2\theta, \qquad \mathcal{M} = r_s r - r_Q^2, \qquad \Delta = r^2 - r_s r + a^2 + r_Q^2 = r^2 + a^2 - \mathcal{M},$$

and $r_s = 2GM/c^2$, $a = J/(Mc)$.

Every appearance of the charge in the metric is an appearance of $\mathcal{M}$.
Written against the Kerr line element of `kerr.md` Step 2 the statement is exact: the two differ by $r_sr \to \mathcal{M}$ in three places and by the same replacement inside $\Delta$, and by nothing else.
$\Sigma$ is untouched, so the ring singularity sits where it sat.

The three parameters are the mass $M$, the spin per unit mass $a$ and the charge radius $r_Q$, and $a$ and $r_Q$ are both lengths.
That is what makes $\Sigma$, $\Delta$ and $\mathcal{M}$ areas, and it is the fact the dimensional pass of Step 16 leans on.

Setting $r_Q = 0$ gives Kerr.
Setting $a = 0$ gives $\Sigma = r^2$ and $\Delta = r^2 - r_sr + r_Q^2$, and the line element collapses to Reissner-Nordstrom in the usual spherical chart, which is the entry `rn_metric` publishes.
Setting both to zero gives Schwarzschild.
Setting $M = 0$ and $r_Q = 0$ together gives flat space in oblate spheroidal coordinates; $M = 0$ with $r_Q \ne 0$ does not, because a charge without a mass still gravitates through the energy of its own field.

The metric depends on $r$ and $\theta$ and on neither $t$ nor $\phi$, so $\partial_t$ and $\partial_\phi$ are Killing vectors: the solution is stationary and axisymmetric, and it is not static, because $g_{t\phi} \ne 0$.

Four derivatives do the work below:

$$\partial_r\Sigma = 2r, \qquad \partial_\theta\Sigma = -2a^2\sin\theta\cos\theta, \qquad \partial_r\mathcal{M} = r_s, \qquad \partial_r\Delta = 2r - r_s,$$

with $\partial_\theta\mathcal{M} = \partial_\theta\Delta = 0$.

The third of these is the one that carries the whole of Step 4.
$\mathcal{M}$ differs from Kerr's $r_sr$ by a constant, so their radial derivatives are equal, and the replacement $r_sr \to \mathcal{M}$ therefore commutes with $\partial_r$ as well as with $\partial_\theta$.
Every quantity built from the metric by differentiating and multiplying is consequently the same function of $\left(r, \cos\theta, a, \mathcal{M}, r_s, \Sigma, \Delta\right)$ for both solutions, and the reader who has `kerr.md` open can read the two connections off one table.
The places where a bare $r_s$ survives are exactly the places where $\partial_r\Delta$ or $\partial_r\mathcal{M}$ was taken, and there the charge does not appear at all.

Two combinations recur often enough to name.
The first is what $\partial_r$ of the mass term produces,

$$\mathcal{P} = \tfrac{1}{2}\left(2r\mathcal{M} - r_s\Sigma\right) = \frac{r_s}{2}\left(r^2 - a^2\cos^2\theta\right) - r_Q^2 r, \qquad \partial_r\!\left(\frac{\mathcal{M}}{\Sigma}\right) = -\frac{2\mathcal{P}}{\Sigma^2},$$

which is Kerr's $\frac{r_s}{2}\left(r^2-a^2\cos^2\theta\right)$ with $-r_Q^2r$ added, and the second is what $\partial_\theta$ produces,

$$\partial_\theta\!\left(\frac{\mathcal{M}}{\Sigma}\right) = \frac{2a^2\mathcal{M}\sin\theta\cos\theta}{\Sigma^2},$$

which needs no new name.
So the connection is built out of $\mathcal{M}$, $\mathcal{P}$, $r_s$, $\Sigma$ and $\Delta$, and the charge enters through the first two.

The chart covers the exterior region.
$\Delta$ vanishes at

$$r_\pm = \frac{r_s}{2} \pm \sqrt{\frac{r_s^2}{4} - a^2 - r_Q^2} = \frac{GM}{c^2} \pm \sqrt{\frac{G^2M^2}{c^4} - a^2 - r_Q^2},$$

which are real and distinct when $a^2 + r_Q^2 < G^2M^2/c^4$, and the domain published is $r > r_+$.
Charge and spin enter that bound in the same way and on the same footing, which is the extremal condition of the solution; Step 13 shows that neither root is a curvature singularity.
The static limit, where $g_{tt}$ vanishes, is at $\Sigma = \mathcal{M}$, that is at $r_E = \frac{r_s}{2} + \sqrt{\frac{r_s^2}{4} - a^2\cos^2\theta - r_Q^2}$, and it lies outside $r_+$ except on the axis, so the ergosphere survives the charge and is thinned by it.

---

## Step 3. The metric matrix, its determinant and its inverse

In the order $(t,r,\theta,\phi)$ the matrix is block diagonal, with a two by two block in the $(t,\phi)$ slots and two diagonal entries between them,

$$g_{\mu\nu} = \begin{pmatrix} -c^2\left(1 - \dfrac{\mathcal{M}}{\Sigma}\right) & 0 & 0 & -\dfrac{c\,\mathcal{M}a\sin^2\theta}{\Sigma} \\[2mm] 0 & \dfrac{\Sigma}{\Delta} & 0 & 0 \\[2mm] 0 & 0 & \Sigma & 0 \\[2mm] -\dfrac{c\,\mathcal{M}a\sin^2\theta}{\Sigma} & 0 & 0 & \left(r^2+a^2+\dfrac{\mathcal{M}a^2\sin^2\theta}{\Sigma}\right)\sin^2\theta \end{pmatrix}.$$

The cross term is read off the line element with the factor of two divided out.

The determinant of the two by two block is the only piece that needs work, and it is worth doing in $\mathcal{M}$ rather than in $r_sr$, because the cancellation is then visibly independent of what the mass function is.
Write $\mathcal{D}$ for it, and note that $g_{tt} = -c^2(\Sigma - \mathcal{M})/\Sigma$ and $g_{\phi\phi} = \left[(r^2+a^2)\Sigma + \mathcal{M}a^2\sin^2\theta\right]\sin^2\theta/\Sigma$.
Then

$$\mathcal{D} = g_{tt}g_{\phi\phi} - g_{t\phi}^2 = -\frac{c^2\sin^2\theta}{\Sigma^2}\Big\{\left(\Sigma - \mathcal{M}\right)\left[(r^2+a^2)\Sigma + \mathcal{M}a^2\sin^2\theta\right] + \mathcal{M}^2a^2\sin^2\theta\Big\},$$

where the last term inside the brace is $-g_{t\phi}^2$ carried across the common factor.
Expanding, the two terms quadratic in $\mathcal{M}$ cancel, and the rest carries a factor of $\Sigma$:

$$\mathcal{D} = -\frac{c^2\sin^2\theta}{\Sigma}\left[(r^2+a^2)\Sigma + \mathcal{M}a^2\sin^2\theta - \mathcal{M}(r^2+a^2)\right],$$

and inside the bracket the two $\mathcal{M}$ terms combine as $-\mathcal{M}\left[(r^2+a^2) - a^2\sin^2\theta\right] = -\mathcal{M}\Sigma$, so the bracket is $\Sigma\left(r^2+a^2-\mathcal{M}\right) = \Sigma\Delta$.
Hence

$$\mathcal{D} = -c^2\Delta\sin^2\theta, \qquad \det g = \frac{\Sigma}{\Delta}\cdot\Sigma\cdot(-c^2\Delta\sin^2\theta) = -c^2\Sigma^2\sin^2\theta,$$

and $\sqrt{-g} = c\,\Sigma\sin\theta$.
The determinant is exactly Kerr's, with no trace of the charge in it at all.
That is not a coincidence: both solutions are Kerr-Schild, built on the same null congruence, and the volume element is the flat one in these coordinates whatever the mass function is.

Inverting a block diagonal matrix is inverting each block:

$$g^{tt} = -\frac{(r^2+a^2)\Sigma + \mathcal{M}a^2\sin^2\theta}{c^2\Sigma\Delta}, \qquad g^{t\phi} = -\frac{\mathcal{M}a}{c\,\Sigma\Delta}, \qquad g^{\phi\phi} = \frac{\Sigma - \mathcal{M}}{\Sigma\Delta\sin^2\theta},$$

together with $g^{rr} = \Delta/\Sigma$ and $g^{\theta\theta} = 1/\Sigma$.

Two rearrangements are used repeatedly below and are what the entry prints:

$$(r^2+a^2)\Sigma + \mathcal{M}a^2\sin^2\theta = \left(r^2+a^2\right)^2 - a^2\Delta\sin^2\theta, \qquad \Sigma - \mathcal{M} = \Delta - a^2\sin^2\theta.$$

Both follow from $\Sigma = r^2 + a^2 - a^2\sin^2\theta$ and $\Delta = r^2+a^2-\mathcal{M}$ and nothing else, which is why they hold for Kerr and for Kerr-Newman in the same words.

---

## Step 4. Every Christoffel symbol with an upper index

The metric depends on $r$ and $\theta$ only, so $\partial_t$ and $\partial_\phi$ annihilate every component, and the inverse metric is nonzero only in the slots $(tt)$, $(t\phi)$, $(rr)$, $(\theta\theta)$, $(\phi\phi)$, so the sum over $\alpha$ in the Christoffel formula has at most two terms.
It is easier to lower first, with

$$\Gamma_{\mu\nu\rho} = \tfrac{1}{2}\left(\partial_\nu g_{\mu\rho} + \partial_\rho g_{\mu\nu} - \partial_\mu g_{\nu\rho}\right),$$

and then raise with $\Gamma^\mu{}_{\nu\rho} = g^{\mu\alpha}\Gamma_{\alpha\nu\rho}$.
Three worked cases show where the charge goes and where it does not.

**A symbol the charge enters through $\mathcal{P}$.**
$\Gamma_{ttr} = \tfrac{1}{2}\partial_r g_{tt}$, and $g_{tt} = -c^2 + c^2\mathcal{M}/\Sigma$, so by the derivative recorded in Step 2,

$$\Gamma_{ttr} = \frac{c^2}{2}\,\partial_r\!\left(\frac{\mathcal{M}}{\Sigma}\right) = -\frac{c^2\mathcal{P}}{\Sigma^2}.$$

Kerr has $-\dfrac{c^2r_s\left(r^2-a^2\cos^2\theta\right)}{2\Sigma^2}$ here, and $\mathcal{P}$ is that numerator with $-r_Q^2r$ added.

**A symbol the charge enters through $\mathcal{M}$.**
$\Gamma_{tt\theta} = \tfrac{1}{2}\partial_\theta g_{tt} = \dfrac{c^2a^2\mathcal{M}\sin\theta\cos\theta}{\Sigma^2}$, by the other derivative of Step 2.
Kerr has $r_sr$ where this has $\mathcal{M}$.

**A symbol the charge does not enter at all.**
$\Gamma_{rrr} = \tfrac{1}{2}\partial_r(\Sigma/\Delta)$, and the quotient rule gives $\dfrac{2r\Delta - \Sigma(2r-r_s)}{2\Delta^2}$, whose numerator is $2r\left(\Delta-\Sigma\right) + r_s\Sigma = 2ra^2\sin^2\theta - 2\mathcal{P}$.
Raising gives

$$\Gamma^r{}_{rr} = \frac{\Delta}{\Sigma}\Gamma_{rrr} = \frac{r}{\Sigma} - \frac{2r-r_s}{2\Delta},$$

which is Kerr's expression unchanged, letter for letter, because the only charge in it is inside $\Delta$.
The $r_s$ standing alone in it is $\partial_r\Delta$ and is not a mass term; replacing it by anything charge dependent would be wrong.

Carrying every case through gives twenty independent symbols, listed with $\nu \le \rho$; the connection is symmetric in its lower pair, so each with $\nu \ne \rho$ stands for two published components and the entry carries thirty two.

$$\Gamma^t{}_{tr} = \frac{\left(r^2+a^2\right)\mathcal{P}}{\Sigma^2\Delta}, \qquad \Gamma^t{}_{t\theta} = -\frac{a^2\mathcal{M}\sin\theta\cos\theta}{\Sigma^2},$$

$$\Gamma^t{}_{r\phi} = -\frac{a\sin^2\theta\left[r\Sigma\mathcal{M} + \left(r^2+a^2\right)\mathcal{P}\right]}{c\,\Sigma^2\Delta}, \qquad \Gamma^t{}_{\theta\phi} = \frac{a^3\mathcal{M}\sin^3\theta\cos\theta}{c\,\Sigma^2},$$

$$\Gamma^r{}_{tt} = \frac{c^2\Delta\,\mathcal{P}}{\Sigma^3}, \qquad \Gamma^r{}_{t\phi} = -\frac{c\,a\Delta\sin^2\theta\,\mathcal{P}}{\Sigma^3},$$

$$\Gamma^r{}_{rr} = \frac{r}{\Sigma} - \frac{2r-r_s}{2\Delta}, \qquad \Gamma^r{}_{r\theta} = -\frac{a^2\sin\theta\cos\theta}{\Sigma}, \qquad \Gamma^r{}_{\theta\theta} = -\frac{r\Delta}{\Sigma},$$

$$\Gamma^r{}_{\phi\phi} = -\frac{\Delta\sin^2\theta\left[r\Sigma^2 - a^2\sin^2\theta\,\mathcal{P}\right]}{\Sigma^3},$$

$$\Gamma^\theta{}_{tt} = -\frac{c^2a^2\mathcal{M}\sin\theta\cos\theta}{\Sigma^3}, \qquad \Gamma^\theta{}_{t\phi} = \frac{c\,a\left(r^2+a^2\right)\mathcal{M}\sin\theta\cos\theta}{\Sigma^3},$$

$$\Gamma^\theta{}_{rr} = \frac{a^2\sin\theta\cos\theta}{\Sigma\Delta}, \qquad \Gamma^\theta{}_{r\theta} = \frac{r}{\Sigma}, \qquad \Gamma^\theta{}_{\theta\theta} = -\frac{a^2\sin\theta\cos\theta}{\Sigma},$$

$$\Gamma^\theta{}_{\phi\phi} = -\frac{\sin\theta\cos\theta\left[\left(r^2+a^2\right)\Sigma^2 + a^2\mathcal{M}\sin^2\theta\left(\Sigma + r^2 + a^2\right)\right]}{\Sigma^3},$$

$$\Gamma^\phi{}_{tr} = \frac{c\,a\,\mathcal{P}}{\Sigma^2\Delta}, \qquad \Gamma^\phi{}_{t\theta} = -\frac{c\,a\mathcal{M}\cos\theta}{\Sigma^2\sin\theta},$$

$$\Gamma^\phi{}_{r\phi} = \frac{r}{\Sigma} - \frac{a^2\sin^2\theta\left[r\Sigma + \mathcal{P}\right]}{\Sigma^2\Delta}, \qquad \Gamma^\phi{}_{\theta\phi} = \cot\theta + \frac{a^2\mathcal{M}\sin\theta\cos\theta}{\Sigma^2}.$$

Seven of the twenty carry no charge except through $\Delta$, and they are the seven with no $t$ index and no $\mathcal{M}$ or $\mathcal{P}$ in them.
Of the rest, six carry $\mathcal{M}$ and seven carry $\mathcal{P}$, and none carries both except $\Gamma^t{}_{r\phi}$, which is the one symbol where a radial derivative of the cross term and a radial derivative of $g_{\phi\phi}$ meet.
Putting $r_Q = 0$ turns $\mathcal{M}$ into $r_sr$ and $\mathcal{P}$ into $\frac{r_s}{2}(r^2-a^2\cos^2\theta)$ and reproduces the twenty symbols of `kerr.md` Step 4 exactly.

---

## Step 5. The lowered Christoffel symbols

Lowering the first index removes the $\Delta$ from most denominators:

$$\Gamma_{ttr} = -\frac{c^2\mathcal{P}}{\Sigma^2}, \qquad \Gamma_{tt\theta} = \frac{c^2a^2\mathcal{M}\sin\theta\cos\theta}{\Sigma^2},$$

$$\Gamma_{tr\phi} = \frac{c\,a\sin^2\theta\,\mathcal{P}}{\Sigma^2} = -\Gamma_{rt\phi}, \qquad \Gamma_{t\theta\phi} = -\frac{c\,a\left(r^2+a^2\right)\mathcal{M}\sin\theta\cos\theta}{\Sigma^2} = \Gamma_{\phi t\theta},$$

$$\Gamma_{rtt} = \frac{c^2\mathcal{P}}{\Sigma^2} = -\Gamma_{ttr}, \qquad \Gamma_{\theta tt} = -\frac{c^2a^2\mathcal{M}\sin\theta\cos\theta}{\Sigma^2} = -\Gamma_{tt\theta},$$

$$\Gamma_{rrr} = \frac{a^2r\sin^2\theta - \mathcal{P}}{\Delta^2}, \qquad \Gamma_{rr\theta} = -\frac{a^2\sin\theta\cos\theta}{\Delta} = -\Gamma_{\theta rr},$$

$$\Gamma_{r\theta\theta} = -r = -\Gamma_{\theta r\theta}, \qquad \Gamma_{\theta\theta\theta} = -a^2\sin\theta\cos\theta,$$

$$\Gamma_{r\phi\phi} = -\frac{\sin^2\theta\left[r\Sigma^2 - a^2\sin^2\theta\,\mathcal{P}\right]}{\Sigma^2} = -\Gamma_{\phi r\phi},$$

$$\Gamma_{\theta\phi\phi} = -\frac{\sin\theta\cos\theta\left[\left(r^2+a^2\right)\Sigma^2 + a^2\mathcal{M}\sin^2\theta\left(\Sigma + r^2+a^2\right)\right]}{\Sigma^2} = -\Gamma_{\phi\theta\phi}.$$

The antisymmetries are not accidents.
$\Gamma_{\mu\nu\rho}$ is symmetric in $\nu\rho$ by construction, and $\Gamma_{\mu\nu\rho} + \Gamma_{\nu\mu\rho} = \partial_\rho g_{\mu\nu}$; whenever $g_{\mu\nu}$ is independent of $x^\rho$ the two are exact negatives, which is what every equality with a minus sign above is saying.
Ten of these thirty two published components are equal up to sign to one of the others, which is why the entry's lowered block is built from ten written values rather than twenty.

---

## Step 6. The electromagnetic field, which is where the charge enters

Kerr-Newman is not a vacuum solution.
It solves the Einstein-Maxwell equations,

$$R_{\mu\nu} - \tfrac{1}{2}Rg_{\mu\nu} = \frac{8\pi G}{c^4}T_{\mu\nu}, \qquad T_{\mu\nu} = \frac{1}{\mu_0}\left(F_{\mu\alpha}F_\nu{}^\alpha - \tfrac{1}{4}g_{\mu\nu}F_{\alpha\beta}F^{\alpha\beta}\right),$$

with a field that has to be written down before the Ricci tensor of Step 7 means anything.
The four potential is the Coulomb potential of the charge dragged around by the rotation,

$$A_\mu dx^\mu = -\frac{Qr}{4\pi\epsilon_0c\,\Sigma}\left(c\,dt - a\sin^2\theta\,d\phi\right),$$

whose $a \to 0$ limit is $A_0 = -Q/(4\pi\epsilon_0cr)$, the Reissner-Nordstrom potential, and whose large $r$ limit is a point charge with a magnetic dipole of moment $Qa$.

Taking $F = dA$ and writing $\kappa = Q/(4\pi\epsilon_0c)$ for the constant in front gives, after the two derivatives,

$$F = \frac{\kappa}{\Sigma^2}\left[-\left(r^2-a^2\cos^2\theta\right)U + 2ar\cos\theta\,V\right],$$

where $U$ and $V$ are the two simple two forms

$$U = \left(c\,dt - a\sin^2\theta\,d\phi\right)\wedge dr, \qquad V = \sin\theta\,d\theta\wedge\left(\left(r^2+a^2\right)d\phi - a\,c\,dt\right).$$

Their nonzero components, which are all that Step 9 needs, are

$$U_{tr} = c, \qquad U_{r\phi} = a\sin^2\theta, \qquad V_{t\theta} = c\,a\sin\theta, \qquad V_{\theta\phi} = \left(r^2+a^2\right)\sin\theta,$$

in the bare chart these steps work in, together with the antisymmetric partners; Step 14 strikes the two factors of $c$ out along with every other one.
$U$ is the radial part of the field and $V$ the angular part: with $a = 0$ the second drops out and the first is the radial electric field of a point charge.

The pair $(U, V)$ is not an ad hoc bookkeeping device.
Write $z = r - ia\cos\theta$, as `kerr.md` Step 7 does, so that $z\bar z = \Sigma$ and $\bar z^2 = \left(r^2-a^2\cos^2\theta\right) + 2iar\cos\theta$.
Then the bracket above is $-\mathrm{Re}\left[\bar z^2\left(U + iV\right)\right]$, so the electromagnetic field and the curvature of Step 10 are built on the same complex structure and the same two principal null directions.
This is what makes Kerr-Newman algebraically special in both its curvature and its field, and it is the reason the whole solution closes in terms of one complex function.

Two facts are checked rather than assumed.
The field is closed, $dF = 0$, which holds because $F = dA$; and it is source free where the chart is defined,

$$\nabla_\mu F^{\mu\nu} = \frac{1}{\sqrt{-g}}\partial_\mu\left(\sqrt{-g}\,F^{\mu\nu}\right) = 0,$$

which is a genuine computation using $\sqrt{-g} = c\Sigma\sin\theta$ from Step 3, and which comes out zero in all four slots.
There is no charge anywhere in the region the chart covers; the charge is a property of the field at infinity, as the mass is.
The invariant is

$$F_{\alpha\beta}F^{\alpha\beta} = -\frac{2\kappa^2}{\Sigma^4}\left[\left(r^2-a^2\cos^2\theta\right)^2 - 4a^2r^2\cos^2\theta\right] = -\frac{2\kappa^2\,\mathrm{Re}\,\bar z^4}{\Sigma^4},$$

which is the usual $E^2 - B^2$ up to constants, and which changes sign, so there are surfaces outside a charged rotating hole where the magnetic part of the field dominates the electric part.

The one combination of constants that survives into the geometry is

$$\frac{8\pi G}{c^4}\cdot\frac{\kappa^2}{\mu_0} = \frac{8\pi G}{c^4}\epsilon_0c^2\cdot\frac{Q^2}{16\pi^2\epsilon_0^2c^2} = \frac{GQ^2}{2\pi\epsilon_0c^4} = 2r_Q^2,$$

using $1/\mu_0 = \epsilon_0c^2$ and the definition of $r_Q$ from Step 1.
That is the entire role of the charge in the gravitational field: it enters as the length $r_Q$ and in no other way, which is why the published entry never mentions $Q$, $\epsilon_0$ or $\mu_0$ at all.

---

## Step 7. The Ricci tensor is the Maxwell stress

The Maxwell stress tensor is trace free in four dimensions, since $g^{\mu\nu}\left(F_{\mu\alpha}F_\nu{}^\alpha - \tfrac14g_{\mu\nu}F^2\right) = F^2 - F^2 = 0$.
Taking the trace of the field equations therefore gives $-R = 0$, and with $R = 0$ the field equations read

$$R_{\mu\nu} = \frac{8\pi G}{c^4}T_{\mu\nu} = 2r_Q^2\left(\hat F_{\mu\alpha}\hat F_\nu{}^\alpha - \tfrac{1}{4}g_{\mu\nu}\hat F_{\alpha\beta}\hat F^{\alpha\beta}\right), \qquad \hat F = F/\kappa,$$

by the constant collected at the end of Step 6.
Carrying out that contraction with the field of Step 6 gives a remarkably compact answer.
Let

$$e^0 = \sqrt{\frac{\Delta}{\Sigma}}\left(c\,dt - a\sin^2\theta\,d\phi\right), \qquad e^1 = \sqrt{\frac{\Sigma}{\Delta}}\,dr,$$

the first two legs of the orthonormal coframe in which $ds^2 = -\left(e^0\right)^2 + \left(e^1\right)^2 + \left(e^2\right)^2 + \left(e^3\right)^2$, with $e^2 = \sqrt{\Sigma}\,d\theta$ and $e^3 = \sin\theta\left[\left(r^2+a^2\right)d\phi - a\,c\,dt\right]/\sqrt{\Sigma}$.
Then

$$R_{\mu\nu} = \mathcal{E}\left(g_{\mu\nu} + 2e^0_\mu e^0_\nu - 2e^1_\mu e^1_\nu\right), \qquad \mathcal{E} = \frac{r_Q^2}{\Sigma^2},$$

which says in the frame that

$$R_{\hat\mu\hat\nu} = \mathcal{E}\,\mathrm{diag}\left(1, -1, 1, 1\right), \qquad R^{\hat\mu}{}_{\hat\nu} = \mathcal{E}\,\mathrm{diag}\left(-1, -1, 1, 1\right).$$

That is the signature of a radial electric field: a positive energy density, a tension of the same size along the direction the field points, a pressure of the same size across it, and a trace of zero.
Reissner-Nordstrom has $R^{\hat\mu}{}_{\hat\nu} = \left(r_q^2/r^4\right)\mathrm{diag}(-1,-1,1,1)$, which is this with $\Sigma \to r^2$, and the `rn_metric` entry publishes exactly those four numbers.

Written out in the chart the entry prints, the six nonzero components are

$$R_{tt} = \frac{r_Q^2\left(\Delta + a^2\sin^2\theta\right)}{\Sigma^3}, \qquad R_{t\phi} = R_{\phi t} = -\frac{a\,r_Q^2\sin^2\theta\left(\Delta + r^2 + a^2\right)}{\Sigma^3},$$

$$R_{rr} = -\frac{r_Q^2}{\Sigma\Delta}, \qquad R_{\theta\theta} = \frac{r_Q^2}{\Sigma}, \qquad R_{\phi\phi} = \frac{r_Q^2\sin^2\theta\left[\left(r^2+a^2\right)^2 + a^2\Delta\sin^2\theta\right]}{\Sigma^3}.$$

The mixed form is shorter and shows the structure better:

$$R^t{}_t = -\frac{r_Q^2\left(r^2+a^2+a^2\sin^2\theta\right)}{\Sigma^3} = -R^\phi{}_\phi, \qquad R^r{}_r = -\frac{r_Q^2}{\Sigma^2} = -R^\theta{}_\theta,$$

$$R^t{}_\phi = \frac{2a\,r_Q^2\left(r^2+a^2\right)\sin^2\theta}{\Sigma^3}, \qquad R^\phi{}_t = -\frac{2a\,r_Q^2}{\Sigma^3},$$

and the upper form follows by one more contraction,

$$R^{tt} = \frac{r_Q^2\left[\left(r^2+a^2\right)^2 + a^2\Delta\sin^2\theta\right]}{\Sigma^3\Delta}, \qquad R^{t\phi} = \frac{a\,r_Q^2\left(\Delta+r^2+a^2\right)}{\Sigma^3\Delta},$$

$$R^{rr} = -\frac{r_Q^2\Delta}{\Sigma^3}, \qquad R^{\theta\theta} = \frac{r_Q^2}{\Sigma^3}, \qquad R^{\phi\phi} = \frac{r_Q^2\left(\Delta + a^2\sin^2\theta\right)}{\Sigma^3\Delta\sin^2\theta}.$$

Every one of the eighteen published Ricci components carries the factor $r_Q^2$, and every one of them vanishes when the charge does.
That is the whole of the difference between this entry's Ricci block and Kerr's empty one, and it is why the two entries are worth having beside each other.
The four slots that vanish, $R_{tr}$, $R_{t\theta}$, $R_{r\phi}$ and $R_{\theta\phi}$, vanish for the same parity reason Step 9 gives for the Riemann tensor, and they would vanish for any stationary axisymmetric solution of this kind.

The scale is worth reading physically.
$R_{\mu\nu}R^{\mu\nu} = 4r_Q^4/\Sigma^4$, so the one invariant the Ricci tensor has is $2\mathcal{E}$ up to a square root, and $\mathcal{E}$ itself is a curvature, $r_Q^2/\Sigma^2$, which at large $r$ falls off as $r^{-4}$, the energy density of a Coulomb field.

---

## Step 8. The Ricci scalar and the Einstein tensor

$R = g^{\mu\nu}R_{\mu\nu} = 0$, which Step 7 established from the trace of the field equations rather than by contracting the table, and which the contraction confirms: the frame form $\mathcal{E}\,\mathrm{diag}(-1,-1,1,1)$ has trace zero by inspection.
The entry prints `R = 0`, and it prints it for a different reason than Kerr does.
Kerr has $R = 0$ because $R_{\mu\nu} = 0$; Kerr-Newman has $R = 0$ because the electromagnetic field is trace free, which is the statement that the photon is massless and that Maxwell theory is conformally invariant in four dimensions.

Because the scalar vanishes,

$$G_{\mu\nu} = R_{\mu\nu} - \tfrac{1}{2}Rg_{\mu\nu} = R_{\mu\nu},$$

slot for slot, and the same holds for $G^\mu{}_\nu$ and $G^{\mu\nu}$.
The entry's three Einstein blocks are therefore the same eighteen values as its three Ricci blocks, and that equality is a fact about this solution rather than a copy: it is what a trace free source does.

Through $G_{\mu\nu} = 8\pi GT_{\mu\nu}/c^4$ the energy density the locally nonrotating observer measures, the one whose four velocity is the dual of $e^0$, is

$$T_{\hat 0\hat 0} = \frac{c^4}{8\pi G}G_{\hat 0\hat 0} = \frac{c^4}{8\pi G}\cdot\frac{r_Q^2}{\Sigma^2},$$

which is positive everywhere, as the energy density of an electromagnetic field must be.
That positivity is the check that the Ricci convention of Step 1 and the signature agree: with the other contraction every sign in Step 7 reverses and this observer measures a negative energy density.
The principal pressures are $-T_{\hat0\hat0}$ along $e^1$ and $+T_{\hat0\hat0}$ along $e^2$ and $e^3$, so the field is a tension along the field lines and a pressure across them, which is Faraday's picture of a stressed field and Maxwell's stress tensor in the form he wrote it.

---

## Step 9. The Riemann tensor, and the one term that separates it from Weyl

Feeding the twenty symbols of Step 4 into the definition of Step 1 and lowering with $R_{\mu\nu\rho\sigma} = g_{\mu\alpha}R^\alpha{}_{\nu\rho\sigma}$ gives thirteen nonzero independent components, the same count as Kerr.
The count follows the same way.
$R_{\mu\nu\rho\sigma}$ is a symmetric matrix on the six index pairs $tr$, $t\theta$, $t\phi$, $r\theta$, $r\phi$ and $\theta\phi$, which is twenty one entries, reduced to twenty by the first Bianchi identity, and eight of the twenty one vanish because the line element is unchanged by $t \to -t$ together with $\phi \to -\phi$, so a component with an odd number of indices from $\{t,\phi\}$ must equal its own negative.
The charge changes none of that, because $\mathcal{M}$ is even under the map exactly as $r_sr$ was.

The first Bianchi identity ties three of the survivors together in the same words as Kerr,

$$R_{tr\theta\phi} - R_{t\theta r\phi} + R_{t\phi r\theta} = c\sin\theta\,\mathcal{B}\left[\left(2r^2+2a^2+a^2\sin^2\theta\right) - \left(r^2+a^2+2a^2\sin^2\theta\right) - \Sigma\right] = 0,$$

with the $\mathcal{B}$ of Step 10, since the bracket is $r^2+a^2-a^2\sin^2\theta-\Sigma$.
It is a check on the whole table and it fails for almost any transcription error.

What is new is that the Riemann tensor is no longer the Weyl tensor.
The difference is fixed by the Ricci tensor of Step 7 through the definition Step 11 quotes, and because that Ricci tensor is $\mathcal{E}\,\mathrm{diag}(1,-1,1,1)$ in the frame, the difference is as simple as a difference of two curvature tensors ever gets:

$$R_{\mu\nu\rho\sigma} = C_{\mu\nu\rho\sigma} + \mathcal{E}\left(U_{\mu\nu}U_{\rho\sigma} + V_{\mu\nu}V_{\rho\sigma}\right),$$

with $U$ and $V$ the two Maxwell two forms of Step 6 and $\mathcal{E} = r_Q^2/\Sigma^2$.
The reason it collapses to two terms is worth seeing.
The correction the Weyl tensor removes is $\tfrac{1}{2}\left(g_{\mu\rho}R_{\sigma\nu} - g_{\mu\sigma}R_{\rho\nu} - g_{\nu\rho}R_{\sigma\mu} + g_{\nu\sigma}R_{\rho\mu}\right)$, and in the orthonormal frame, where both $g$ and $R$ are diagonal, its only nonzero entries are on the diagonal of the six by six matrix of index pairs, with the value $\tfrac{1}{2}\left(\eta_{aa}R_{bb} + \eta_{bb}R_{aa}\right)$ on the pair $ab$.
Substituting $R_{\hat\mu\hat\nu} = \mathcal{E}\,\mathrm{diag}(1,-1,1,1)$ gives $\mathcal{E}$ on the pair $01$, $\mathcal{E}$ on the pair $23$, and zero on the other four.
The two planes that survive are exactly the planes $U$ and $V$ span, and the square roots in the coframe cancel when the two are written in coordinates, which is why $U$ and $V$ carry no $\Delta$ or $\Sigma$ at all.

That leaves six of the thirteen independent components carrying a charge term and seven carrying none:

$$R_{trtr} = C_{trtr} + \frac{c^2r_Q^2}{\Sigma^2}, \qquad R_{trr\phi} = C_{trr\phi} + \frac{c\,a\,r_Q^2\sin^2\theta}{\Sigma^2}, \qquad R_{t\theta t\theta} = C_{t\theta t\theta} + \frac{c^2a^2r_Q^2\sin^2\theta}{\Sigma^2},$$

$$R_{t\theta\theta\phi} = C_{t\theta\theta\phi} + \frac{c\,a\left(r^2+a^2\right)r_Q^2\sin^2\theta}{\Sigma^2}, \qquad R_{r\phi r\phi} = C_{r\phi r\phi} + \frac{a^2r_Q^2\sin^4\theta}{\Sigma^2},$$

$$R_{\theta\phi\theta\phi} = C_{\theta\phi\theta\phi} + \frac{\left(r^2+a^2\right)^2r_Q^2\sin^2\theta}{\Sigma^2},$$

and $R_{trt\theta}$, $R_{tr\theta\phi}$, $R_{t\theta r\phi}$, $R_{t\phi t\phi}$, $R_{t\phi r\theta}$, $R_{r\theta r\theta}$ and $R_{r\phi\theta\phi}$ are equal to their Weyl counterparts.
The seven that are untouched are the five that carry $\mathcal{B}$, which lie in none of the two planes, together with $R_{t\phi t\phi}$ and $R_{r\theta r\theta}$, whose planes are $03$ and $12$.
The mixed components obey the same formula with the first index of each bivector raised,

$$R^\mu{}_{\nu\rho\sigma} = C^\mu{}_{\nu\rho\sigma} + \mathcal{E}\left(U^\mu{}_\nu U_{\rho\sigma} + V^\mu{}_\nu V_{\rho\sigma}\right),$$

where

$$U^t{}_r = -\frac{r^2+a^2}{c\,\Delta}, \quad U^r{}_t = -\frac{c\,\Delta}{\Sigma}, \quad U^r{}_\phi = \frac{a\Delta\sin^2\theta}{\Sigma}, \quad U^\phi{}_r = -\frac{a}{\Delta},$$

$$V^t{}_\theta = -\frac{a\sin\theta}{c}, \quad V^\theta{}_t = -\frac{c\,a\sin\theta}{\Sigma}, \quad V^\theta{}_\phi = \frac{\left(r^2+a^2\right)\sin\theta}{\Sigma}, \quad V^\phi{}_\theta = -\frac{1}{\sin\theta}.$$

The entry publishes each Riemann component as its Weyl component plus its charge term, in that order, so that the reader can see which part of the curvature is the free gravitational field and which part is the field of the charge sitting on top of it.

---

## Step 10. The two curvature functions

Define

$$\mathcal{A} = \frac{\dfrac{r_s}{2}r\left(r^2-3a^2\cos^2\theta\right) - r_Q^2\left(r^2-a^2\cos^2\theta\right)}{\Sigma^3}, \qquad \mathcal{B} = \frac{a\cos\theta\left[\dfrac{r_s}{2}\left(3r^2-a^2\cos^2\theta\right) - 2r_Q^2r\right]}{\Sigma^3}.$$

Setting $r_Q = 0$ gives the $\mathcal{A}$ and $\mathcal{B}$ of `kerr.md` Step 7 exactly.
As there, the two are one complex function.
With $z = r - ia\cos\theta$ and $z\bar z = \Sigma$,

$$\mathcal{A} + i\mathcal{B} = \frac{1}{z^3}\left(\frac{r_s}{2} - \frac{r_Q^2}{\bar z}\right) = \frac{\bar z^2\left(\dfrac{r_s}{2}\bar z - r_Q^2\right)}{\Sigma^3} = \frac{1}{z^3}\left(\frac{GM}{c^2} - \frac{r_Q^2}{\bar z}\right).$$

The middle form is the one to expand.
$\bar z^2 = p + iq$ with $p = r^2 - a^2\cos^2\theta$ and $q = 2ar\cos\theta$, and $\frac{r_s}{2}\bar z - r_Q^2 = \left(\frac{r_sr}{2} - r_Q^2\right) + i\frac{r_sa\cos\theta}{2}$, so the real part of the product is

$$p\left(\frac{r_sr}{2}-r_Q^2\right) - \frac{q\,r_sa\cos\theta}{2} = \frac{r_sr}{2}\left(r^2-3a^2\cos^2\theta\right) - r_Q^2\left(r^2-a^2\cos^2\theta\right),$$

since $-2a^2r\cos^2\theta \cdot \frac{r_s}{2} \cdot 2$ is what turns $r^2 - a^2\cos^2\theta$ into $r^2 - 3a^2\cos^2\theta$, and the imaginary part is

$$q\left(\frac{r_sr}{2}-r_Q^2\right) + \frac{p\,r_sa\cos\theta}{2} = a\cos\theta\left[\frac{r_s}{2}\left(3r^2-a^2\cos^2\theta\right) - 2r_Q^2r\right],$$

which are the numerators above.
Kerr's single cube $r_s/(2z^3)$ has become a cube with one extra pole in $\bar z$, and that single extra term is the entire effect of the charge on the free gravitational field.

With those two in hand the thirteen independent lowered Weyl components are, in the bare chart,

$$C_{trtr} = -c^2\,\frac{2\Delta + a^2\sin^2\theta}{\Delta}\,\mathcal{A}, \qquad C_{trt\theta} = 3ac^2\sin\theta\;\mathcal{B},$$

$$C_{trr\phi} = -c\,\frac{a\sin^2\theta\left[3\left(r^2+a^2\right) - 2\mathcal{M}\right]}{\Delta}\,\mathcal{A}, \qquad C_{tr\theta\phi} = c\sin\theta\left[2\left(r^2+a^2\right)+a^2\sin^2\theta\right]\mathcal{B},$$

$$C_{t\theta t\theta} = c^2\left(\Delta + 2a^2\sin^2\theta\right)\mathcal{A}, \qquad C_{t\theta r\phi} = c\sin\theta\left[r^2+a^2+2a^2\sin^2\theta\right]\mathcal{B},$$

$$C_{t\theta\theta\phi} = c\,a\sin^2\theta\left[3\left(r^2+a^2\right)-\mathcal{M}\right]\mathcal{A}, \qquad C_{t\phi t\phi} = c^2\Delta\sin^2\theta\;\mathcal{A}, \qquad C_{t\phi r\theta} = -c\,\Sigma\sin\theta\;\mathcal{B},$$

$$C_{r\theta r\theta} = -\frac{\Sigma^2}{\Delta}\,\mathcal{A}, \qquad C_{r\phi r\phi} = -\frac{\sin^2\theta\left[\left(r^2+a^2\right)^2 + 2a^2\Delta\sin^2\theta\right]}{\Delta}\,\mathcal{A},$$

$$C_{r\phi\theta\phi} = 3a\left(r^2+a^2\right)\sin^3\theta\;\mathcal{B}, \qquad C_{\theta\phi\theta\phi} = \sin^2\theta\left[2\left(r^2+a^2\right)^2 + a^2\Delta\sin^2\theta\right]\mathcal{A}.$$

These are Kerr's thirteen expressions, word for word, with the charged $\mathcal{A}$, $\mathcal{B}$, $\mathcal{M}$ and $\Delta$ in place of the uncharged ones.
Eight carry $\mathcal{A}$ and five carry $\mathcal{B}$, and no component carries both; the split is the reflection $\theta \to \pi - \theta$, under which $\mathcal{A}$ is even and $\mathcal{B}$ is odd.

The mixed components $C^\mu{}_{\nu\rho\sigma} = g^{\mu\alpha}C_{\alpha\nu\rho\sigma}$ follow by one contraction each, and there are more of them, eighty eight published from thirty eight distinct values, because raising an index with a metric whose $(t,\phi)$ block is not diagonal turns one component into a mixture of two.
Two of the thirty eight are quadratic in the mass function for that reason, $C^t{}_{tt\phi}$ and $C^t{}_{tr\theta}$ together with their partners, since the surviving term carries $g^{t\phi}$, which is itself proportional to $\mathcal{M}$.
Kerr prints those with $G^2M^2$ in the numerator; here they are printed as the product of the bracket $\mathcal{M}$ with the bracket $\Sigma^3\mathcal{A}$ or $\Sigma^3\mathcal{B}/(a\cos\theta)$, because the two factors no longer combine into a single power of the mass.

---

## Step 11. The Weyl tensor, computed rather than copied

The Weyl tensor is Riemann with its traces removed,

$$C_{\mu\nu\rho\sigma} = R_{\mu\nu\rho\sigma} - \frac{1}{n-2}\left(g_{\mu\rho}R_{\sigma\nu} - g_{\mu\sigma}R_{\rho\nu} - g_{\nu\rho}R_{\sigma\mu} + g_{\nu\sigma}R_{\rho\mu}\right) + \frac{R}{(n-1)(n-2)}\left(g_{\mu\rho}g_{\sigma\nu} - g_{\mu\sigma}g_{\rho\nu}\right),$$

with $n = 4$, and the traces it removes are Riemann's own.

For Kerr every term of the correction vanishes and the two tensors are equal.
Here the last term vanishes, because $R = 0$ by Step 8, and the middle term does not, because $R_{\mu\nu}$ does not.
Step 9 evaluated it: the correction is $\mathcal{E}\left(U\otimes U + V\otimes V\right)$ and it is nonzero in six of the thirteen independent slots.
The Weyl tensor of Kerr-Newman is therefore not its Riemann tensor, and the entry publishes two genuinely different blocks where Kerr publishes two identical ones.

This is the place where an entry can go wrong by copying, and the reason `_tools/derivations/weyl.md` tells a reader to compute rather than copy.
Had the Weyl block here been filled from the Riemann block, six of the thirteen independent values would have been wrong, and wrong by an amount that vanishes only in the limit that takes the solution back to Kerr.
The verification script builds the Weyl tensor from the definition above, forming the correction from its own Riemann tensor and its own Ricci tensor, and never reads the published Riemann block when it checks the published Weyl block.
The working here was done the same way, and the closed form of Step 9 was then recognised in the answer rather than assumed at the start.

Kerr-Newman is algebraically special, of Petrov type D, and the statement is about this tensor.
Four of the five complex Weyl scalars vanish and the whole tensor is carried by

$$\Psi_2 = -\left(\mathcal{A} + i\mathcal{B}\right) = -\frac{1}{\left(r-ia\cos\theta\right)^3}\left(\frac{GM}{c^2} - \frac{r_Q^2}{r+ia\cos\theta}\right),$$

which is the cube of Step 10 seen from the other side.
The electromagnetic field of Step 6 is aligned with the same two principal null directions, which is what a charged type D solution means and what makes the Newman-Janis trick of the entry's history work at all.

---

## Step 12. The Kretschmann scalar

$$K = R_{\mu\nu\rho\sigma}R^{\mu\nu\rho\sigma}.$$

The decomposition of Step 9 makes this two pieces rather than one.
Squaring $R = C + \mathcal{E}\left(U\otimes U + V\otimes V\right)$ and using the general identity

$$R_{\mu\nu\rho\sigma}R^{\mu\nu\rho\sigma} = C_{\mu\nu\rho\sigma}C^{\mu\nu\rho\sigma} + 2R_{\mu\nu}R^{\mu\nu} - \tfrac{1}{3}R^2,$$

which holds for any metric, together with $R = 0$ and $R_{\mu\nu}R^{\mu\nu} = 4r_Q^4/\Sigma^4$ from Step 7, gives

$$K = C_{\mu\nu\rho\sigma}C^{\mu\nu\rho\sigma} + \frac{8r_Q^4}{\Sigma^4}.$$

The Weyl square is the type D result, which for the thirteen components of Step 10 comes out as

$$C_{\mu\nu\rho\sigma}C^{\mu\nu\rho\sigma} = 48\left(\mathcal{A}^2 - \mathcal{B}^2\right) = 48\,\mathrm{Re}\left[\left(\mathcal{A}+i\mathcal{B}\right)^2\right] = \frac{48}{\Sigma^6}\,\mathrm{Re}\left[\bar z^4\left(\frac{r_s}{2}\bar z - r_Q^2\right)^2\right],$$

by Step 10, exactly as it does for Kerr, since the contraction never asks what $\mathcal{A}$ and $\mathcal{B}$ are made of.
Expanding the square leaves three powers of $\bar z$ to take the real part of,

$$\bar z^4\left(\frac{r_s}{2}\bar z - r_Q^2\right)^2 = \frac{r_s^2}{4}\bar z^6 - r_sr_Q^2\bar z^5 + r_Q^4\bar z^4,$$

and with $v = a\cos\theta$,

$$\mathrm{Re}\,\bar z^6 = r^6 - 15r^4v^2 + 15r^2v^4 - v^6, \qquad \mathrm{Re}\,\bar z^5 = r\left(r^4 - 10r^2v^2 + 5v^4\right), \qquad \mathrm{Re}\,\bar z^4 = r^4 - 6r^2v^2 + v^4.$$

Adding the Ricci piece, whose $\Sigma^2 = r^4 + 2r^2v^2 + v^4$ goes over the same denominator, collects the two terms in $r_Q^4$ into one:

$$48\left(r^4-6r^2v^2+v^4\right) + 8\left(r^4+2r^2v^2+v^4\right) = 8\left(7r^4 - 34r^2v^2 + 7v^4\right).$$

Therefore

$$K = \frac{12r_s^2\left(r^6-15a^2r^4\cos^2\theta+15a^4r^2\cos^4\theta-a^6\cos^6\theta\right) - 48r_sr_Q^2r\left(r^4-10a^2r^2\cos^2\theta+5a^4\cos^4\theta\right) + 8r_Q^4\left(7r^4-34a^2r^2\cos^2\theta+7a^4\cos^4\theta\right)}{\Sigma^6},$$

which with $r_s = 2GM/c^2$ is what the entry publishes.
A scalar is the same number in every chart, so this expression needs no conversion in Step 14.

Three checks.
Setting $r_Q = 0$ leaves the first term, and $r^6-15r^4v^2+15r^2v^4-v^6 = \left(r^2-v^2\right)\left(\Sigma^2-16r^2v^2\right)$, which is the form `kerr.md` Step 11 publishes.
Setting $a = 0$ gives $\Sigma = r^2$ and

$$K = \frac{12r_s^2r^6 - 48r_sr_Q^2r^5 + 56r_Q^4r^4}{r^{12}} = \frac{12r_s^2r^2 - 48r_sr_Q^2r + 56r_Q^4}{r^8},$$

which is exactly the Kretschmann scalar the `rn_metric` entry publishes for Reissner-Nordstrom.
The dimensions work out as $L^2\cdot L^6$, $L^2\cdot L^2\cdot L\cdot L^4$ and $L^4\cdot L^4$ over $L^{12}$, which is $L^{-4}$ in all three terms, what a squared curvature carries.

---

## Step 13. What the Kretschmann scalar says

$\Delta$ does not appear in $K$.
Neither horizon is therefore a place where anything measurable diverges, and the failure of the Boyer-Lindquist chart at $r = r_\pm$ is a failure of the chart, exactly as it is for Kerr and for Reissner-Nordstrom.

The only way $K$ can diverge is $\Sigma \to 0$, and $\Sigma = r^2 + a^2\cos^2\theta$ is a sum of two squares, so it vanishes only where $r = 0$ and $\cos\theta = 0$ together.
That is the ring singularity, unmoved by the charge, since the charge never enters $\Sigma$.
What the charge does change is how fast the divergence comes on.
Approaching the ring sends both $r$ and $v = a\cos\theta$ to zero together, so the term of lowest degree in the numerator wins, and that is the $r_Q^4$ term, which is quartic where the cross term is quintic and the mass term is sextic.
Kerr diverges as $\Sigma^{-3}$ and Kerr-Newman as $\Sigma^{-4}$: a charged rotating hole is singular in the same place as an uncharged one, but it is the energy of the electromagnetic field rather than the mass that sets the leading divergence there.

The sign structure is richer than Kerr's for the same reason.
$K$ changes sign across surfaces where the three terms balance, and in the uncharged limit those surfaces are the cone $r = a\left|\cos\theta\right|$ and the two surfaces $\Sigma^2 = 16r^2a^2\cos^2\theta$ of `kerr.md` Step 12.
On the axis, where $\cos\theta = \pm 1$ and $\Sigma = a^2$ at $r = 0$, $K$ is finite and equal to $\left(-12r_s^2a^6 + 56r_Q^4a^4\right)/a^{12}$, so the interior of the disc $r = 0$ is ordinary spacetime that a worldline can be continued through, and the charge shifts the value there without making it singular.

---

## Step 14. The same components in the chart the entry prints

The entry prints components in the chart $x^0 = ct$.
By Step 1 the rule is: multiply by $c$ once for each upper $t$ index and divide by $c$ once for each lower one.

The metric.
$g_{tt}$ has two lower $t$ indices and becomes $-\left(1-\mathcal{M}/\Sigma\right)$, $g_{t\phi}$ has one and becomes $-\mathcal{M}a\sin^2\theta/\Sigma$, and $g_{rr}$, $g_{\theta\theta}$, $g_{\phi\phi}$ have none and are unchanged.
The inverse runs the other way: $g^{tt}$ gains $c^2$ and becomes $-\left[\left(r^2+a^2\right)^2 - a^2\Delta\sin^2\theta\right]/(\Sigma\Delta)$, $g^{t\phi}$ gains one $c$ and becomes $-\mathcal{M}a/(\Sigma\Delta)$, and the rest are unchanged.
Every factor of $c$ in the bare tables of Step 3 is removed exactly, and no published metric component carries one.

The connection and the curvature.
Every symbol and every component in Steps 4, 5, 9 and 10 carries $c$ to the power of the number of its lower $t$ indices minus the number of its upper ones, so the published tables are those formulas with the $c$ powers struck out.
That is true of the charge terms as well: $U_{tr} = c$ has one lower $t$ and $U^t{}_r = -\left(r^2+a^2\right)/(c\Delta)$ has one upper, which is why the first is printed as $1$ and the second as $-\left(r^2+a^2\right)/\Delta$, with no $c$ in either.
A chart all of whose coordinates are lengths or angles cannot produce a factor of $c$, and $c$ is the only quantity in the problem with a time in it.

The parameters.
$r_s$ becomes $2GM/c^2$ everywhere, so the four blocks of Step 2 and Step 10 are published as

$$\mathcal{M} \to \left(\frac{2GMr}{c^2} - r_Q^2\right), \qquad \mathcal{P} \to \left(\frac{GM}{c^2}\left(r^2-a^2\cos^2\theta\right) - r_Q^2r\right),$$

$$\Sigma^3\mathcal{A} \to \left(\frac{GMr}{c^2}\left(r^2-3a^2\cos^2\theta\right) - r_Q^2\left(r^2-a^2\cos^2\theta\right)\right), \qquad \frac{\Sigma^3\mathcal{B}}{a\cos\theta} \to \left(\frac{GM}{c^2}\left(3r^2-a^2\cos^2\theta\right) - 2r_Q^2r\right),$$

and $\Sigma$ and $\Delta$ are written out in full, because the file's reader knows only the symbols the entry declares, which are $M$, $a$, $G$ and $r_Q$.
The numerators of $\mathcal{A}$ and $\mathcal{B}$ are what appear, rather than $\mathcal{A}$ and $\mathcal{B}$ themselves, because the entry cannot name a function; the powers of $\Sigma$ they carry are gathered into the denominator of each component.
Two identities between the four are worth recording, because they are what a reader checking a component by hand will want:

$$\Sigma^3\mathcal{A} = r\mathcal{P} - a^2\cos^2\theta\,\mathcal{M}, \qquad \frac{\Sigma^3\mathcal{B}}{a\cos\theta} = \mathcal{P} + r\mathcal{M}.$$

Both are immediate from the definitions and both fail if a single $r_Q^2$ is misplaced.

---

## Step 15. The geodesic equations

The equation of motion is

$$\ddot{x}^\mu + \Gamma^\mu{}_{\nu\rho}\dot{x}^\nu\dot{x}^\rho = 0,$$

with the dot a derivative with respect to an affine parameter $\lambda$ and the components those of the chart the entry prints, so $\dot{t}$ means $d(ct)/d\lambda$ and the $\Gamma$ in it is the published one.
Each symbol with $\nu \ne \rho$ contributes twice, once from each ordering, which is the factor of two in front of the mixed terms.
Written with the blocks of Step 2 and the symbols of Step 4:

$$\ddot{t} + \frac{2\left(r^2+a^2\right)\mathcal{P}}{\Sigma^2\Delta}\dot{t}\dot{r} - \frac{2a^2\mathcal{M}\sin\theta\cos\theta}{\Sigma^2}\dot{t}\dot{\theta} - \frac{2a\sin^2\theta\left[r\Sigma\mathcal{M} + \left(r^2+a^2\right)\mathcal{P}\right]}{\Sigma^2\Delta}\dot{r}\dot{\phi} + \frac{2a^3\mathcal{M}\sin^3\theta\cos\theta}{\Sigma^2}\dot{\theta}\dot{\phi} = 0,$$

$$\ddot{r} + \frac{\Delta\mathcal{P}}{\Sigma^3}\dot{t}^2 - \frac{2a\Delta\sin^2\theta\,\mathcal{P}}{\Sigma^3}\dot{t}\dot{\phi} + \left(\frac{r}{\Sigma} - \frac{2r - r_s}{2\Delta}\right)\dot{r}^2 - \frac{2a^2\sin\theta\cos\theta}{\Sigma}\dot{r}\dot{\theta} - \frac{r\Delta}{\Sigma}\dot{\theta}^2 - \frac{\Delta\sin^2\theta\left[r\Sigma^2 - a^2\sin^2\theta\,\mathcal{P}\right]}{\Sigma^3}\dot{\phi}^2 = 0,$$

$$\ddot{\theta} - \frac{a^2\mathcal{M}\sin\theta\cos\theta}{\Sigma^3}\dot{t}^2 + \frac{2a\left(r^2+a^2\right)\mathcal{M}\sin\theta\cos\theta}{\Sigma^3}\dot{t}\dot{\phi} + \frac{a^2\sin\theta\cos\theta}{\Sigma\Delta}\dot{r}^2 + \frac{2r}{\Sigma}\dot{r}\dot{\theta} - \frac{a^2\sin\theta\cos\theta}{\Sigma}\dot{\theta}^2 - \frac{\sin\theta\cos\theta\left[\left(r^2+a^2\right)\Sigma^2 + a^2\mathcal{M}\sin^2\theta\left(\Sigma+r^2+a^2\right)\right]}{\Sigma^3}\dot{\phi}^2 = 0,$$

$$\ddot{\phi} + \frac{2a\mathcal{P}}{\Sigma^2\Delta}\dot{t}\dot{r} - \frac{2a\mathcal{M}\cos\theta}{\Sigma^2\sin\theta}\dot{t}\dot{\theta} + 2\left(\frac{r}{\Sigma} - \frac{a^2\sin^2\theta\left[r\Sigma + \mathcal{P}\right]}{\Sigma^2\Delta}\right)\dot{r}\dot{\phi} + \frac{2\cos\theta\left[\Sigma^2 + a^2\mathcal{M}\sin^2\theta\right]}{\Sigma^2\sin\theta}\dot{\theta}\dot{\phi} = 0.$$

The entry prints the same four equations with $r_s$ replaced by $2GM/c^2$ and the blocks written out, and the $\dot{r}^2$ term of the $r$ equation is printed as $r/\Sigma - \left(c^2r - GM\right)/\left(c^2\Delta\right)$, which is the same number written so that each bracket balances dimensionally on its own.

These are the geodesics of an uncharged test particle.
A charged particle does not follow them, because it feels the field of Step 6 through the Lorentz force, and its equation of motion carries a term $\left(q/mc\right)F^\mu{}_\nu\dot{x}^\nu$ on the right hand side.
The entry publishes the geodesics, as every other entry in the collection does, and the distinction is worth stating because in this spacetime, unlike in Kerr, it is a real one.

The $t$ and $\phi$ equations have no term in $\dot{t}^2$, $\dot{\phi}^2$ or $\dot{t}\dot{\phi}$, which is the Killing structure showing through: the conserved energy and angular momentum are first integrals of exactly those two equations.
Carter's constant survives the charge, which is what makes Kerr-Newman geodesics separable in the same way Kerr's are, and it is again invisible in this form because it comes out of the Hamilton-Jacobi equation rather than out of a Killing vector.

---

## Step 16. Every published equation is dimensionally consistent

The chart coordinates are $ct$, $r$, $\theta$ and $\phi$, of which the first two are lengths and the last two are dimensionless.
The parameters are $[M] = M$, $[a] = L$, $[G] = L^3M^{-1}T^{-2}$ and $[r_Q] = L$, so that $[GM/c^2] = L$ and $[\Sigma] = [\Delta] = [\mathcal{M}] = L^2$, while $[\mathcal{P}] = L^3$ and $[\Sigma^3\mathcal{A}] = L^4$.
A component of a tensor carries, on top of the dimension its rank gives it, a factor $[x^\mu]/L$ for each upper index and $L/[x^\mu]$ for each lower one.

The four blocks are each a sum, so each has to balance internally, and that is where the charge declaration earns its keep.
In $\mathcal{M}$ the two terms are $2GMr/c^2$, which is $L \cdot L$, and $r_Q^2$, which is $L^2$.
In $\mathcal{P}$ they are $\left(GM/c^2\right)\left(r^2 - a^2\cos^2\theta\right)$, which is $L\cdot L^2$, and $r_Q^2r$, which is $L^3$.
In $\Sigma^3\mathcal{A}$ they are $L\cdot L\cdot L^2$ and $L^2\cdot L^2$, both $L^4$, and in $\Sigma^3\mathcal{B}/(a\cos\theta)$ they are $L \cdot L^2$ and $L^2 \cdot L$, both $L^3$.
In $\Delta$ all four terms are $L^2$.
Every one of them would fail if $r_Q$ were taken to be dimensionless, which is the one thing about this entry the script cannot read off the file and which its line in `DIMENSIONS` declares.

Take the Ricci tensor, which is the block Kerr does not have.
$\left[R_{\theta\theta}\right] = L^{-2}\cdot L \cdot L = 1$, since each lower $\theta$ contributes $L/[\theta] = L$; and $r_Q^2/\Sigma$ carries $L^2/L^2 = 1$.
$\left[R^r{}_r\right] = L^{-2}$, with the upper and lower $r$ cancelling, and $r_Q^2/\Sigma^2$ carries $L^2/L^4 = L^{-2}$.
$\left[R^{t\phi}\right] = L^{-2}\cdot\left(L/L\right)\cdot\left(1/L\right) = L^{-3}$, against $a r_Q^2\left(\Delta+r^2+a^2\right)/\left(\Sigma^3\Delta\right)$ at $L\cdot L^2\cdot L^2/\left(L^6\cdot L^2\right) = L^{-3}$.

Take the $\phi$ geodesic equation term by term, whose left hand side $\ddot\phi$ carries $1/\lambda^2$.
The first term is $2\Gamma^\phi{}_{tr}\dot t\dot r$, and $\left[\Gamma^\phi{}_{tr}\right] = L^{-1}\cdot(1/L) = L^{-2}$, matching $2a\mathcal{P}/\left(\Sigma^2\Delta\right)$ at $L\cdot L^3/\left(L^4\cdot L^2\right) = L^{-2}$; the velocities carry $\left(L/\lambda\right)^2$ and the term is $1/\lambda^2$.
The second is $2\Gamma^\phi{}_{t\theta}\dot t\dot\theta$, with $\left[\Gamma^\phi{}_{t\theta}\right] = L^{-1}(1/L)(1)(L) = L^{-1}$, matching $2a\mathcal{M}\cos\theta/\left(\Sigma^2\sin\theta\right)$ at $L\cdot L^2/L^4 = L^{-1}$, and the velocities carry $\left(L/\lambda\right)\left(1/\lambda\right)$, so the term is $1/\lambda^2$ again.
It would not have been, had the dot on $t$ been read as $dt/d\lambda$ rather than $d(ct)/d\lambda$: that reading leaves this term short by one factor of $c$ while leaving the first alone, so no overall rescaling could repair both.
This is the same argument `kasner.md` Step 17 and `kerr.md` Step 15 make, and it is what pins the dot convention down without any algebra.

The checker's dimensional pass does this for every term of every published expression in the entry, which is 448 components, four geodesic equations, one line element and one scalar, in well under a second.

---

## Step 17. What the entry publishes, and what the checker needs

The Boyer-Lindquist system carries:

- six metric components and six inverse metric components, the extra pair over the diagonal being $g_{t\phi} = g_{\phi t}$ and its inverse;
- thirty two Christoffel symbols with an upper index and thirty two with all three lowered, from twenty independent values each;
- eighty eight components of $R^\mu{}_{\nu\rho\sigma}$ and eighty of $R_{\mu\nu\rho\sigma}$, from thirty eight and thirteen distinct values, each printed as its Weyl component plus its charge term;
- the same counts again for the Weyl tensor, which here is not the Riemann tensor, by Step 11;
- six components each of $R_{\mu\nu}$, $R^\mu{}_\nu$ and $R^{\mu\nu}$, and the same eighteen again for the Einstein tensor, which equals the Ricci tensor because $R = 0$;
- `R = 0`, by Step 8;
- the Kretschmann scalar of Step 12;
- four geodesic equations.

That is 448 published components, against Kerr's 412; the difference is the thirty six Ricci and Einstein components that a vacuum solution does not have.

Adding the spacetime to the checker took one entry in `DIMENSIONS` in `verify_metrics.py`,

```
("kerr_newman", "boyer_lindquist"): {
    "t": "T", "r": "L", "\\theta": "1", "\\phi": "1",
    "M": "M", "a": "L", "G": "L**3/(M*T**2)", "r_Q": "L",
},
```

which is the documented step for adding a spacetime, and which declares the two things the script cannot read off the file: that $a$ and $r_Q$ are lengths rather than dimensionless numbers, and that $t$ is the coordinate the chart multiplies by $c$.
No entry in `PARAMETER_RELATIONS` is needed, because $M$, $a$, $G$ and $r_Q$ are free: every published value here is an identity for all four, not an identity on a constraint surface.
The extremality condition $a^2 + r_Q^2 \le G^2M^2/c^4$ is a condition for the chart's domain to be what the entry says it is, not a relation the algebra leans on.

Two practical points, inherited from the Kerr entry and confirmed again here.

The index names in a published component have to be the coordinate names exactly as `coords` spells them, which for this entry means `"\\theta"` and `"\\phi"` rather than `"theta"` and `"phi"`.
Getting it wrong is loud in one direction and quiet in the other: `compare_block` reports every misspelled index as a disagreement, but the dimensional pass skips a component whose indices it cannot resolve, so a clean `--dimensions-only` run means nothing on a new entry until the sympy comparison has confirmed that the index names resolve.

The reader splits a name it does not know into single letters, so `r_Q` has to be kept clear of the letter before it.
`ar_Q^2` parses as $a \cdot r \cdot {}_{-} \cdot Q^2$ and fails; `a\,r_Q^2` parses as intended, and every published value here that multiplies $r_Q$ by a preceding symbol carries that thin space.
A digit before it is safe, so `2r_Q^2` needs nothing, and `96GM\,r_Q^2r` needs one.

---

## Step 18. What the checker managed, and what it could not

Kerr-Newman is slow to check for the same reason Kerr is, and for one reason more.
`norm`, the routine that puts an expression into a form where a vanishing one is recognisably zero, calls sympy's general `simplify`, which is cheap on a diagonal metric and very expensive on this one: the metric is not diagonal, so the Christoffel sums do not collapse, and every curvature component is a rational function with $\Sigma$ to a high power underneath.
On top of that, the checker works in the bare chart, where the mass appears as $G$, $M$ and $c$ separately rather than as the one length $GM/c^2$, so its rational functions carry two more variables than the algebra needs.

Run as it stands, with `--system kerr_newman/boyer_lindquist`, the dimensional pass finishes in about a second and the sympy pass gets through the metric, the inverse metric and both Christoffel variants and then abandons the Riemann tensor on the budget; everything downstream of Riemann is abandoned with it, since the Ricci, Einstein and Weyl blocks are all built from it.
That is reported as `UNCHECKED` with the reason rather than passed in silence, so nothing below is being claimed on the script's authority that the script did not establish.

What was established, and what a later reader should know stands behind the numbers above:

- The dimensional pass runs on this entry in well under a second with `--dimensions-only`, and every term of all 448 components, the four geodesic equations, the line element and the Kretschmann scalar balances. Since the same pass skips silently over any component whose index names it cannot resolve, that is worth only as much as the sympy pass below, which confirms that all 448 resolve.
- `check_system` itself was then run with only `norm` replaced. Everything else was the script: its `Reader`, its metric read off the published line element, its Christoffel, Riemann, Ricci, Einstein and Weyl built from that metric, its chart weighting, and its component by component comparison, which also requires every component the entry does not publish to vanish. The replacement makes two changes of representation and no change of mathematics: it reduces $\sin\theta$ and $\cos\theta$ to a canonical rational form modulo $\sin^2+\cos^2=1$ instead of calling `simplify`, and it substitutes $M \to \mu c^2/G$ so that the one length the mass enters through is one symbol, which is a bijective change of a positive parameter and so preserves exactly which expressions are zero.
- On that run the script reported **zero disagreements and zero dimensional failures**, over the metric, the inverse metric, both Christoffel variants, both Riemann variants, all three Ricci variants, all three Einstein variants, both Weyl variants, the Ricci scalar and the four geodesic equations. The whole of it took 364 seconds of work: 81 seconds on $R^\mu{}_{\nu\rho\sigma}$, 60 on $R_{\mu\nu\rho\sigma}$, 88 on $C^\mu{}_{\nu\rho\sigma}$, 58 on $C_{\mu\nu\rho\sigma}$, and seconds or less on everything else.
- The Kretschmann scalar is the one thing that did not finish, exactly as it does not for Kerr. Forming it means raising four indices through a metric that is not diagonal and normalising 256 components on the way, and it ran past twenty minutes without finishing. It was checked instead by exact rational arithmetic: the published scalar and the contraction $R_{\mu\nu\rho\sigma}R^{\mu\nu\rho\sigma}$ of the computed Riemann tensor were evaluated at four random rational points, with $\theta$ carried by the Weierstrass parametrisation $\sin\theta = 2u/(1+u^2)$, $\cos\theta = (1-u^2)/(1+u^2)$ so that $\sin^2+\cos^2=1$ holds exactly, and agreed at every one. The same was done for $C_{\mu\nu\rho\sigma}C^{\mu\nu\rho\sigma}$ against $K - 8r_Q^4/\Sigma^4$, which is the decomposition Step 12 derives, and it agreed at every one as well. Both quantities are rational functions of $(r, a, r_s, r_Q, u)$, so agreement at random points is agreement of the functions unless the difference vanishes on a variety through all of them.
- The electromagnetic field of Step 6 was checked in the same way and symbolically: $dF = 0$ and $\nabla_\mu F^{\mu\nu} = 0$ in all four slots, $F$ equals the combination of $U$ and $V$ that Step 6 quotes, and the Ricci tensor of Step 7 equals both $2r_Q^2\left(\hat F_{\mu\alpha}\hat F_\nu{}^\alpha - \tfrac14 g_{\mu\nu}\hat F^2\right)$ and $\mathcal{E}\left(g + 2e^0e^0 - 2e^1e^1\right)$, in every slot.
- The split of Step 9, $R = C + \mathcal{E}\left(U\otimes U + V\otimes V\right)$, was checked in all 256 slots of both variants, and it is what the published Riemann block is built from.

The honest summary is that the physics is verified and the tool is slower than the physics.
The fix is the same one `kerr.md` Step 17 names: giving `norm` a canonical trig reduction to fall back on before it reaches for `simplify` would bring both of these entries inside the budget, and collapsing $GM/c^2$ to one symbol before the algebra would bring this one inside it comfortably.
That is a change to the checker, not to a spacetime, and it belongs to whoever takes the checker on next.
