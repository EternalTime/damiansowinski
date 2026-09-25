# The Weyl tensors of five spacetimes

We compute the Weyl tensors of FRW, the interior Schwarzschild solution, Godel's universe, Reissner-Nordstrom and the Lanczos-van Stockum dust.

In all five, the Weyl components had been copied from the Riemann components.
That is wrong in every case, because the Weyl tensor is Riemann with its traces taken out, and none of these five spacetimes has a traceless Riemann tensor.
Two of them are conformally flat, so their Weyl tensors are zero and every copied component was wrong.
The other three have nonzero Weyl tensors, and in two of them the copy was wrong in a second way as well: it had no component in a slot where Riemann happens to vanish but Weyl does not.

The script `_tools/derivations/verify_metrics.py` computes all of this in sympy and compares each component with its own result, so every value can be checked by hand and by machine independently.

---

## Step 1. Conventions

The signature is $(-,+,+,+)$ and the Riemann tensor is

$$R^\mu{}_{\nu\rho\sigma} = \partial_\rho \Gamma^\mu_{\nu\sigma} - \partial_\sigma \Gamma^\mu_{\nu\rho} + \Gamma^\mu_{\rho\lambda}\Gamma^\lambda_{\nu\sigma} - \Gamma^\mu_{\sigma\lambda}\Gamma^\lambda_{\nu\rho},$$

and every Riemann component is taken in this convention.

The Ricci tensor is the standard contraction on the first lower index,

$$R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu},$$

which is the one the signature asks for and the one every Ricci tensor, Einstein tensor and Ricci scalar is taken in.

The Weyl tensor wants exactly that trace, because Weyl is defined by subtracting from Riemann the parts of Riemann that its own traces carry, and those traces are a property of the Riemann tensor.
So the $R_{\mu\nu}$ and $R$ in the Weyl formula are the Ricci tensor and the Ricci scalar themselves, with no sign to keep track of.
Before 18 September 2026 the Ricci tensor was contracted on the last index, which is minus this.
On that older convention the trace in the Weyl formula was the negative of the Ricci tensor, and every Weyl value was already written for the trace rather than for the Ricci tensor.
The Weyl values are therefore unchanged by the move; only the sign relating them to the Ricci tensor has gone away.

The chart is $x^0 = cT$, as for every other spacetime.
Of the five spacetimes, only FRW in its comoving spherical chart has a coordinate that carries dimensions of time, and its Weyl tensor is zero, so no factor of $c$ survives in any Weyl component.
In the other four the chart does not rescale the time coordinate at all, Godel's being dimensionless and the rest already lengths, so their components are components of the coordinates as written in the line element.

---

## Step 2. The trace removal

In $n$ dimensions the Weyl tensor is

$$C_{\mu\nu\rho\sigma} = R_{\mu\nu\rho\sigma} - \frac{1}{n-2}\left(g_{\mu\rho}R_{\sigma\nu} - g_{\mu\sigma}R_{\rho\nu} - g_{\nu\rho}R_{\sigma\mu} + g_{\nu\sigma}R_{\rho\mu}\right) + \frac{R}{(n-1)(n-2)}\left(g_{\mu\rho}g_{\sigma\nu} - g_{\mu\sigma}g_{\rho\nu}\right),$$

and in the four dimensions all five spacetimes have, the two denominators are $2$ and $6$:

$$C_{\mu\nu\rho\sigma} = R_{\mu\nu\rho\sigma} - \tfrac{1}{2}\left(g_{\mu\rho}R_{\sigma\nu} - g_{\mu\sigma}R_{\rho\nu} - g_{\nu\rho}R_{\sigma\mu} + g_{\nu\sigma}R_{\rho\mu}\right) + \tfrac{1}{6}R\left(g_{\mu\rho}g_{\sigma\nu} - g_{\mu\sigma}g_{\rho\nu}\right).$$

This is the formula `verify_metrics.py` uses, in `Geometry.weyl_llll`, with $R_{\mu\nu}$ supplied by `ricci_ll`.

Three consequences recur throughout.

**The correction is not small.**
$C = R$ holds only when $R_{\mu\nu} = 0$, that is, only in a Ricci flat spacetime.
Of the five spacetimes, not one is Ricci flat: FRW and the interior Schwarzschild solution are filled with a fluid, Godel and Lanczos-van Stockum with rotating dust, and Reissner-Nordstrom with an electromagnetic field.
In Schwarzschild $C = R$ genuinely holds, because it is a vacuum, and carrying that relation over to a spacetime which is not a vacuum is exactly the mistake corrected in all five.

**Weyl can be nonzero where Riemann vanishes.**
The correction terms are built from the metric and the traces, not from the Riemann component in that slot, so a slot in which $R_{\mu\nu\rho\sigma} = 0$ generally has $C_{\mu\nu\rho\sigma} \neq 0$.
Godel and Lanczos-van Stockum both do this, which is why copying Riemann left components missing as well as wrong.

**Weyl has the symmetries of Riemann.**
It is antisymmetric in $\mu\nu$, antisymmetric in $\rho\sigma$, symmetric under exchange of the two pairs, and satisfies the first Bianchi identity, so one component in each symmetry class determines all of its images.

---

## Step 3. Friedmann Robertson Walker, which is conformally flat

FRW is written in two charts,

$$ds^2 = -c^2dt^2 + a^2\left(\frac{dr^2}{1-kr^2} + r^2d\Omega^2\right), \qquad ds^2 = a^2\left(-d\eta^2 + \frac{dr^2}{1-kr^2} + r^2d\Omega^2\right),$$

and the copy gave 24 nonzero Weyl components in each index position of each chart, 96 values, all of them copies of that chart's Riemann tensor.
Every one of them is wrong, and the Weyl tensor vanishes in both charts.

The argument is the standard isotropy one, and it needs no components at all.

Relative to the unit timelike vector $u^\mu$ of a comoving observer, the Weyl tensor of any four dimensional spacetime splits into two pieces,

$$E_{\mu\nu} = C_{\mu\alpha\nu\beta}u^\alpha u^\beta, \qquad H_{\mu\nu} = \tfrac{1}{2}\epsilon_{\mu\alpha\beta\gamma}C^{\alpha\beta}{}_{\nu\delta}u^\gamma u^\delta,$$

its electric and magnetic parts.
Both are symmetric, both are trace free, and both are orthogonal to $u^\mu$, so each is a symmetric trace free tensor on the three dimensional space the observer calls space.
Between them they carry all ten independent components of $C$, so $C = 0$ if and only if $E = 0$ and $H = 0$.

Now use the hypothesis that defines an FRW spacetime.
It is isotropic about the comoving worldline through each point, which means the full rotation group $SO(3)$ acting on that observer's three space is a symmetry of the geometry at that point, and so leaves every tensor built from the geometry alone invariant.
$E_{\mu\nu}$ and $H_{\mu\nu}$ are built from the geometry alone.
An $SO(3)$ invariant symmetric bilinear form on three dimensions is a multiple of $\delta_{ij}$, by Schur's lemma, and a multiple of $\delta_{ij}$ that is trace free is zero.
Hence $E = H = 0$ and the Weyl tensor vanishes at every point, for every $k$ and every scale factor $a$.

This is what "conformally flat" means here: $C = 0$ is the integrability condition, in four dimensions, for the metric to be $\Omega^2$ times a flat one.
For $k = 0$ the second chart displays the factor, since there $a^2(\eta)$ multiplies Minkowski in spherical coordinates; for $k = \pm 1$ that chart reduces the problem to the static metric $-d\eta^2 + d\Sigma_k^2$, which takes one further conformal map to reach flat space.

The vanishing is a fact about the spacetime, not about the chart, so it holds in both charts, and the conventions of FRW now say so.
The FRW Riemann tensor has its own separate errors, 41 of them, which are catalogued in `audit-2026-09-18.md` and which the Weyl correction leaves untouched.

---

## Step 4. The interior Schwarzschild solution, which is also conformally flat

$$ds^2 = -\frac{1}{4}\left(3\sqrt{1-\frac{r_s}{R}} - \sqrt{1-\frac{r^2r_s}{R^3}}\right)^{\!2}dt^2 + \frac{dr^2}{1-\dfrac{r^2r_s}{R^3}} + r^2d\Omega^2$$

is the uniform density static star, and its copied Weyl tensor had 24 components in each index position, 48 values, all of them wrong.
Here the isotropy argument of Step 3 is not available, since the solution is not homogeneous, so this one is done by computing.

Write any static spherically symmetric metric as

$$ds^2 = -e^{2\Phi(r)}dt^2 + \frac{dr^2}{1 - \dfrac{2m(r)}{r}} + r^2d\Omega^2,$$

which defines the mass function $m(r)$ from $g^{rr}$, in units with $G = c = 1$, the units of the line element.
Such a metric has one independent Weyl component, and carrying the formula of Step 2 through gives

$$C_{trtr} = \frac{e^{2\Phi}}{3r^2(r-2m)}\Big[(r^3 - 2r^2m)\left(\Phi'' + \Phi'^2\right) - r^2\Phi'm' - r^2\Phi' + 3rm\Phi' + rm' - 3m\Big].$$

The last two terms of the bracket are the interesting ones.
The $tt$ field equation for a static sphere is $m'(r) = 4\pi r^2\rho(r)$, so

$$rm' - 3m = -3r^3\left(\frac{m(r)}{r^3} - \frac{4\pi}{3}\rho(r)\right),$$

which is three $r^3$ times the difference between the mean density inside radius $r$ and the density at $r$.
That difference is what a tidal field measures, and for a star of uniform density it is zero at every radius.

Put the interior solution in.
Its $g^{rr}$ gives $2m/r = r_s r^2/R^3$, so with the abbreviation $k \equiv r_s/R^3$,

$$m(r) = \tfrac{1}{2}kr^3, \qquad m'(r) = \tfrac{3}{2}kr^2, \qquad rm' - 3m = \tfrac{3}{2}kr^3 - \tfrac{3}{2}kr^3 = 0,$$

which is the uniform density statement: $\rho = 3k/(8\pi)$ is a constant.
The two remaining terms carrying $m$ cancel against each other too,

$$-r^2\Phi'm' + 3rm\Phi' = -\tfrac{3}{2}kr^4\Phi' + \tfrac{3}{2}kr^4\Phi' = 0,$$

so the whole bracket collapses to

$$r^3\left(1-kr^2\right)\left(\Phi'' + \Phi'^2\right) - r^2\Phi',$$

and the Weyl tensor vanishes if and only if

$$r\left(1-kr^2\right)\left(\Phi'' + \Phi'^2\right) = \Phi'.$$

This is now a statement about $g_{tt}$ alone, and it is easiest checked on $e^\Phi$ rather than on $\Phi$, because $\left(e^\Phi\right)'' = \left(\Phi'' + \Phi'^2\right)e^\Phi$ turns it into the linear condition

$$r\left(1-kr^2\right)\left(e^\Phi\right)'' = \left(e^\Phi\right)'.$$

The interior solution's $g_{tt}$ is exactly $-e^{2\Phi}$ with

$$e^{\Phi} = \tfrac{1}{2}\left(3A - B\right), \qquad A \equiv \sqrt{1-kR^2} = \sqrt{1-\frac{r_s}{R}}, \qquad B \equiv \sqrt{1-kr^2} = \sqrt{1-\frac{r^2r_s}{R^3}},$$

in which $A$ is a constant and only $B$ depends on $r$.
Differentiating $B^2 = 1-kr^2$ gives $B' = -kr/B$, and differentiating again,

$$B'' = -\frac{k}{B} + \frac{kr B'}{B^2} = -\frac{k}{B} - \frac{k^2r^2}{B^3} = -\frac{k\left(1-kr^2\right) + k^2r^2}{B^3} = -\frac{k}{B^3}.$$

So $\left(e^\Phi\right)' = -\tfrac{1}{2}B' = \dfrac{kr}{2B}$ and $\left(e^\Phi\right)'' = -\tfrac{1}{2}B'' = \dfrac{k}{2B^3}$, and therefore

$$r\left(1-kr^2\right)\left(e^\Phi\right)'' = rB^2\cdot\frac{k}{2B^3} = \frac{kr}{2B} = \left(e^\Phi\right)',$$

which is the condition, identically in $r$.
The Weyl tensor of the interior Schwarzschild solution is zero in both index positions, and its conventions now say so.

Two remarks are worth keeping.
The constant $3A$ is fixed by requiring the pressure to vanish at $r = R$, so conformal flatness here is not an accident of a free constant: it holds for the physical uniform density star.
And the exterior it matches onto has $m(r) = r_s/2$ constant, so there $rm' - 3m = -3r_s/2 \neq 0$ and the Weyl tensor is not zero, which is Step 5 with $r_q = 0$.

---

## Step 5. Reissner-Nordstrom

$$ds^2 = -\frac{\Delta}{r^2}dt^2 + \frac{r^2}{\Delta}dr^2 + r^2d\Omega^2, \qquad \Delta \equiv r^2 - r_sr + r_q^2.$$

This is the mildest of the five.
The copy had all 24 nonzero components in each index position, and 16 of the 24 were already right, because in those slots the trace correction happens to vanish.
The 8 that disagreed in each index position are corrected, and no other component moved.

The same general formula serves.
Here $e^{2\Phi} = f \equiv \Delta/r^2 = 1 - 2m/r$, so the metric is one of the family of Step 4 with

$$m(r) = \frac{r_s}{2} - \frac{r_q^2}{2r},$$

and with $g_{tt}g_{rr} = -1$ the bracket simplifies a long way.
Substituting $\Phi = \tfrac{1}{2}\ln f$, so that

$$\Phi' = \frac{f'}{2f}, \qquad \Phi'' + \Phi'^2 = \frac{f''}{2f} - \frac{f'^2}{4f^2},$$

into the bracket of Step 4, and using $m = r(1-f)/2$ and $m' = (1-f)/2 - rf'/2$, every term carrying $f'^2$ cancels and what is left is

$$\tfrac{1}{2}r^3f'' - r^2f' + rf - r.$$

With $f = 1 - \dfrac{r_s}{r} + \dfrac{r_q^2}{r^2}$ the three pieces are

$$\tfrac{1}{2}r^3f'' = -r_s + \frac{3r_q^2}{r}, \qquad -r^2f' = -r_s + \frac{2r_q^2}{r}, \qquad r(f-1) = -r_s + \frac{r_q^2}{r},$$

which add to $-3\left(r_s - \dfrac{2r_q^2}{r}\right)$, and so

$$C_{trtr} = \frac{f}{3r^2\cdot rf}\cdot\left(-3\right)\left(r_s - \frac{2r_q^2}{r}\right) = \frac{2r_q^2 - r r_s}{r^4},$$

which is the corrected value.
Setting $r_q = 0$ gives $-r_s/r^3$, which is the same component of Schwarzschild, as it must be, since Schwarzschild is Ricci flat and there $C = R$.

The same answer falls out of the trace removal directly, which is the shorter check.
Reissner-Nordstrom has $R = 0$, so the whole $R$ term of Step 2 drops, and its Ricci tensor gives

$$R_{tt} = \frac{r_q^2\Delta}{r^6}, \qquad R_{rr} = -\frac{r_q^2}{r^2\Delta}.$$

Since $g_{tr} = 0$, only two of the four correction terms survive, and both come to the same thing:

$$g_{tt}R_{rr} = \left(-\frac{\Delta}{r^2}\right)\left(-\frac{r_q^2}{r^2\Delta}\right) = \frac{r_q^2}{r^4}, \qquad g_{rr}R_{tt} = \frac{r^2}{\Delta}\cdot\frac{r_q^2\Delta}{r^6} = \frac{r_q^2}{r^4}.$$

So the correction is $\tfrac{1}{2}\left(2r_q^2/r^4\right) = r_q^2/r^4$, and

$$C_{trtr} = R_{trtr} - \frac{r_q^2}{r^4} = \frac{3r_q^2 - rr_s}{r^4} - \frac{r_q^2}{r^4} = \frac{2r_q^2 - rr_s}{r^4}.$$

The Riemann component is Reissner-Nordstrom's own, so this is a check of the curvature against itself, and it is precisely the term the copy dropped.

This also says which 16 components the copy got right by luck.
The correction in the $t\theta t\theta$, $t\phi t\phi$, $r\theta r\theta$ and $r\phi r\phi$ classes is a difference of two equal terms and vanishes, for instance

$$\tfrac{1}{2}\left(g_{tt}R_{\theta\theta} + g_{\theta\theta}R_{tt}\right) = \tfrac{1}{2}\left(-\frac{\Delta}{r^2}\cdot\frac{r_q^2}{r^2} + r^2\cdot\frac{r_q^2\Delta}{r^6}\right) = 0,$$

so in those 16 slots Weyl and Riemann agree and the copy was harmless.
In the other two classes it does not vanish: the $trtr$ class is off by $r_q^2/r^4$ and the $\theta\phi\theta\phi$ class by

$$\tfrac{1}{2}\left(g_{\theta\theta}R_{\phi\phi} + g_{\phi\phi}R_{\theta\theta}\right) = r_q^2\sin^2\theta,$$

which turns $R_{\theta\phi\theta\phi} = \left(rr_s - r_q^2\right)\sin^2\theta$ into $C_{\theta\phi\theta\phi} = \left(rr_s - 2r_q^2\right)\sin^2\theta$.
Those two classes are the 8 corrected components of each index position.

The six independent components in each index position are

| class | $C_{\mu\nu\rho\sigma}$ | $C^\mu{}_{\nu\rho\sigma}$ |
| --- | --- | --- |
| $trtr$ | $\dfrac{2r_q^2 - rr_s}{r^4}$ | $\dfrac{rr_s - 2r_q^2}{r^2\Delta}$ |
| $t\theta t\theta$ | $\dfrac{\Delta\left(rr_s - 2r_q^2\right)}{2r^4}$ | $\dfrac{2r_q^2 - rr_s}{2r^2}$ |
| $t\phi t\phi$ | $\dfrac{\Delta\left(rr_s - 2r_q^2\right)\sin^2\theta}{2r^4}$ | $\dfrac{\left(2r_q^2 - rr_s\right)\sin^2\theta}{2r^2}$ |
| $r\theta r\theta$ | $\dfrac{2r_q^2 - rr_s}{2\Delta}$ | $\dfrac{2r_q^2 - rr_s}{2r^2}$ |
| $r\phi r\phi$ | $\dfrac{\left(2r_q^2 - rr_s\right)\sin^2\theta}{2\Delta}$ | $\dfrac{\left(2r_q^2 - rr_s\right)\sin^2\theta}{2r^2}$ |
| $\theta\phi\theta\phi$ | $\left(rr_s - 2r_q^2\right)\sin^2\theta$ | $\dfrac{\left(rr_s - 2r_q^2\right)\sin^2\theta}{r^2}$ |

Each is a multiple of the single quantity $rr_s - 2r_q^2$, which is the whole physical content: it is $2r^2$ times the effective tidal mass $m_W(r) = \dfrac{r_s}{2} - \dfrac{r_q^2}{r}$, which is the gravitational mass inside $r$ less the field energy that is not there.
The ratios $-2 : 1 : 1$ between the $tr$, $t\theta$ and $t\phi$ planes in the orthonormal frame are the statement that the electric Weyl tensor is trace free, which is Petrov type D.

---

## Step 6. Godel

$$ds^2 = \frac{1}{2\omega^2}\left[-dt^2 + dx^2 - 2e^x\,dt\,dy - \tfrac{1}{2}e^{2x}dy^2 + dz^2\right],$$

with the metric, inverse metric and Ricci tensor

$$g_{tt} = -\frac{1}{2\omega^2}, \quad g_{ty} = -\frac{e^x}{2\omega^2}, \quad g_{xx} = g_{zz} = \frac{1}{2\omega^2}, \quad g_{yy} = -\frac{e^{2x}}{4\omega^2},$$

$$R_{tt} = 1, \quad R_{ty} = e^x, \quad R_{yy} = e^{2x}, \quad R = -2\omega^2,$$

and every other component of $R_{\mu\nu}$ is zero.
It is rank one, as the Ricci tensor of dust has to be: with the dust four velocity $u^\mu = \sqrt{2}\,\omega\,\delta^\mu_t$, which is a unit vector because $g_{tt} = -\dfrac{1}{2\omega^2}$, its covariant form is $u_\mu = -\dfrac{1}{\sqrt{2}\,\omega}\left(1,0,e^x,0\right)$ and

$$R_{\mu\nu} = 2\omega^2u_\mu u_\nu,$$

which reproduces all three components and, contracted with $u^\mu u^\nu = -1$, gives $R = -2\omega^2$ as well.
In particular $R_{xx} = R_{zz} = 0$, which keeps the trace corrections short.

**A copied component.**
Take $\mu\nu\rho\sigma = txtx$.
The correction term is

$$\tfrac{1}{2}\left(g_{tt}R_{xx} - g_{tx}R_{tx} - g_{xt}R_{xt} + g_{xx}R_{tt}\right) = \tfrac{1}{2}\left(0 - 0 - 0 + \frac{1}{2\omega^2}\cdot 1\right) = \frac{1}{4\omega^2},$$

since $R_{xx} = 0$ and $g_{tx} = 0$.
The trace term is

$$\tfrac{1}{6}R\left(g_{tt}g_{xx} - g_{tx}g_{xt}\right) = \tfrac{1}{6}\left(-2\omega^2\right)\left(-\frac{1}{2\omega^2}\cdot\frac{1}{2\omega^2}\right) = \frac{1}{12\omega^2}.$$

Godel's Riemann tensor has $R_{txtx} = \dfrac{1}{4\omega^2}$, so

$$C_{txtx} = \frac{1}{4\omega^2} - \frac{1}{4\omega^2} + \frac{1}{12\omega^2} = \frac{1}{12\omega^2},$$

where the copy said $\dfrac{1}{4\omega^2}$.
The correction happens to cancel the Riemann term exactly here, so the whole of this component is the Ricci scalar term.

**A component the copy left out.**
Take $tztz$.
Godel's Riemann tensor has no $R_{tztz}$, and correctly so: it is zero.
But $R_{zz} = 0$ and $g_{tz} = 0$ give the same correction as before,

$$\tfrac{1}{2}g_{zz}R_{tt} = \frac{1}{4\omega^2}, \qquad \tfrac{1}{6}R\,g_{tt}g_{zz} = \frac{1}{12\omega^2},$$

so

$$C_{tztz} = 0 - \frac{1}{4\omega^2} + \frac{1}{12\omega^2} = -\frac{1}{6\omega^2},$$

which is nonzero.
Godel's $z$ axis is the rotation axis, and its Riemann tensor has no component at all carrying a $z$, so the copy had none either: that is where 20 of the 40 lowered components and 18 of the mixed ones went missing.
This is the clearest demonstration that Weyl is not Riemann: the trace removal can put curvature into a slot that Riemann leaves empty.

The eight independent lowered components, and the ten mixed ones, are

| class | $C_{\mu\nu\rho\sigma}$ | class | $C^\mu{}_{\nu\rho\sigma}$ |
| --- | --- | --- | --- |
| $txtx$ | $\dfrac{1}{12\omega^2}$ | $ttty$ | $\dfrac{e^x}{6}$ |
| $txxy$ | $-\dfrac{e^x}{12\omega^2}$ | $txtx$ | $-\dfrac{1}{6}$ |
| $tyty$ | $\dfrac{e^{2x}}{24\omega^2}$ | $txxy$ | $\dfrac{e^x}{2}$ |
| $tztz$ | $-\dfrac{1}{6\omega^2}$ | $tyty$ | $\dfrac{e^{2x}}{12}$ |
| $tzyz$ | $-\dfrac{e^x}{6\omega^2}$ | $tztz$ | $\dfrac{1}{3}$ |
| $xyxy$ | $\dfrac{e^{2x}}{6\omega^2}$ | $tzyz$ | $\dfrac{e^x}{2}$ |
| $xzxz$ | $-\dfrac{1}{12\omega^2}$ | $xyxy$ | $\dfrac{e^{2x}}{3}$ |
| $yzyz$ | $-\dfrac{5e^{2x}}{24\omega^2}$ | $xzxz$ | $-\dfrac{1}{6}$ |
| | | $yyty$ | $-\dfrac{e^x}{6}$ |
| | | $yzyz$ | $-\dfrac{1}{6}$ |

Each index position has 40 nonzero components, the images of these under the symmetries of Step 2.
The mixed components have two extra classes because raising the first index with $g^{ty} \neq 0$ turns a $y$ into a $t$, so $C^t{}_{tty}$ and $C^y{}_{yty}$ are nonzero although no lowered component has a repeated index.

The Weyl tensor is dimensionally uniform: with $x$, $y$, $z$ and $t$ all dimensionless and $[\omega] = 1/L$, every lowered component carries $L^2$ and every mixed one carries nothing, which is what `--dimensions-only` checks.

---

## Step 7. Lanczos-van Stockum

$$ds^2 = -dt^2 - \frac{2r^2}{R}\,dt\,d\phi + e^{-r^2/R^2}dr^2 + r^2\left(1-\frac{r^2}{R^2}\right)d\phi^2 + e^{-r^2/R^2}dz^2,$$

the rigidly rotating dust cylinder, with the Ricci tensor

$$R_{tt} = \frac{2E}{R^2}, \quad R_{t\phi} = \frac{2Er^2}{R^3}, \quad R_{rr} = R_{zz} = \frac{2}{R^2}, \quad R_{\phi\phi} = \frac{2Er^2\left(r^2+R^2\right)}{R^4}, \quad R = \frac{4E}{R^2},$$

writing $E \equiv e^{r^2/R^2}$ throughout.

**A worked component.**
Take $trtr$, where $g_{tr} = 0$, $g_{tt} = -1$ and $g_{rr} = 1/E$.
The correction is

$$\tfrac{1}{2}\left(g_{tt}R_{rr} + g_{rr}R_{tt}\right) = \tfrac{1}{2}\left(\left(-1\right)\frac{2}{R^2} + \frac{1}{E}\cdot\frac{2E}{R^2}\right) = \tfrac{1}{2}\left(-\frac{2}{R^2} + \frac{2}{R^2}\right) = 0,$$

and the trace term is

$$\tfrac{1}{6}R\,g_{tt}g_{rr} = \tfrac{1}{6}\cdot\frac{4E}{R^2}\cdot\left(-1\right)\cdot\frac{1}{E} = -\frac{2}{3R^2}.$$

Its Riemann tensor has $R_{trtr} = \dfrac{1}{R^2}$, so

$$C_{trtr} = \frac{1}{R^2} - 0 - \frac{2}{3R^2} = \frac{1}{3R^2},$$

which is the corrected value, where the copy said $\dfrac{1}{R^2}$.
Here the correction term vanishes on its own and the Ricci scalar term does all the work, taking exactly two thirds of the Riemann component away.

As in Godel, the axial direction is where components went missing: $R_{tztz} = 0$ while

$$C_{tztz} = 0 - \tfrac{1}{2}\left(g_{tt}R_{zz} + g_{zz}R_{tt}\right) + \tfrac{1}{6}Rg_{tt}g_{zz} = -\tfrac{1}{2}\left(-\frac{2}{R^2} + \frac{2}{R^2}\right) - \frac{2}{3R^2} = -\frac{2}{3R^2},$$

which the copy left out.
That one class is the whole of what was missing, four lowered components and two mixed ones; the other 36 and 38 were present and wrong.

The independent components are

| class | $C_{\mu\nu\rho\sigma}$ | class | $C^\mu{}_{\nu\rho\sigma}$ |
| --- | --- | --- | --- |
| $trtr$ | $\dfrac{1}{3R^2}$ | $ttt\phi$ | $\dfrac{Er^2}{3R^3}$ |
| $trr\phi$ | $-\dfrac{4r^2}{3R^3}$ | $trtr$ | $-\dfrac{R^2+3r^2}{3R^4}$ |
| $t\phi t\phi$ | $\dfrac{Er^2}{3R^2}$ | $trr\phi$ | $\dfrac{r^2\left(2R^2+r^2\right)}{R^5}$ |
| $tztz$ | $-\dfrac{2}{3R^2}$ | $t\phi t\phi$ | $-\dfrac{Er^2\left(R^2-r^2\right)}{3R^4}$ |
| $tz\phi z$ | $-\dfrac{5r^2}{3R^3}$ | $tztz$ | $\dfrac{2R^2+3r^2}{3R^4}$ |
| $r\phi r\phi$ | $\dfrac{r^2\left(2R^2+7r^2\right)}{3R^4}$ | $tz\phi z$ | $\dfrac{r^2\left(2R^2+r^2\right)}{R^5}$ |
| $rzrz$ | $-\dfrac{1}{3R^2 E}$ | $r\phi r\phi$ | $\dfrac{Er^2\left(2R^2+7r^2\right)}{3R^4}$ |
| $\phi z\phi z$ | $-\dfrac{r^2\left(R^2+8r^2\right)}{3R^4}$ | $rzrz$ | $-\dfrac{1}{3R^2}$ |
| | | $\phi\phi t\phi$ | $-\dfrac{Er^2}{3R^3}$ |
| | | $\phi z\phi z$ | $-\dfrac{R^2+3r^2}{3R^4}$ |

There are 40 nonzero lowered components and 44 mixed ones, the extra four being the images of $C^t{}_{tt\phi}$ and $C^\phi{}_{\phi t\phi}$, which exist for the same reason as Godel's: $g^{t\phi} \neq 0$, so raising an index mixes $t$ with $\phi$.

Two checks on these components can be done by eye.
Every lowered component carries $1/L^2$ with $[r] = [R] = [t] = [z] = L$ and $\phi$ dimensionless, and $E$ is dimensionless because $r^2/R^2$ is.
And every component vanishes as $R \to \infty$ at fixed $r$, which is the flat limit: the dust density that the Ricci scalar carries is $\rho = \dfrac{E}{2\pi R^2}$, on the sign convention in which the field equation reads $G_{\mu\nu} = -8\pi T_{\mu\nu}$, and that goes to zero there too.

---

## Step 8. What the checker says

`verify_metrics.py` reported 433 disagreements over all spacetimes before the Weyl correction and 113 after it.
The two runs were diffed line by line rather than compared as totals: 320 lines were removed, no line was added, and every one of the 320 names a Weyl component of one of these five spacetimes.

| spacetime and chart | Weyl disagreements removed |
| --- | ---: |
| FRW, comoving spherical | 48 |
| FRW, conformal spherical | 48 |
| Godel, Cartesian | 80 |
| interior Schwarzschild, spherical | 48 |
| Reissner-Nordstrom, spherical | 16 |
| Lanczos-van Stockum, cylindrical | 80 |
| | **320** |

Each 80 is 40 components in each of the two index positions, each 48 is 24 in each, and the 16 is 8 in each.

The 113 that remain are the rest of `audit-2026-09-18.md`, none of it touched by the Weyl correction.

| spacetime | left | what it is |
| --- | ---: | --- |
| interior Schwarzschild | 51 | the 50 nested radicals, which are the checker's limitation and not the physics, and the inverse metric component that does not invert the line element |
| FRW | 41 | the Riemann tensor, its index positions and signs and omissions, and the doubly raised Ricci tensor |
| Ellis-Bronnikov | 13 | curvature components it leaves out, 10 of them in its own Weyl tensor, and one geodesic whose sign contradicts its own Christoffel symbol |
| Godel | 4 | the doubly raised Einstein tensor, which is written unsymmetric |
| Lanczos-van Stockum | 4 | the same |

Ellis-Bronnikov is the sixth spacetime with Weyl disagreements, and the Weyl correction deliberately leaves it alone.
Its Weyl tensor is the one the previous pass got right: every component present is correct, and its 10 disagreements are components it omits, such as $C_{\theta t t\theta} = -\dfrac{\ell^2}{3\left(\ell^2+r^2\right)}$.
Supplying missing components is a different job from correcting copied ones.

The dimensional pass reports the same 7 terms before and after, in the FRW Riemann tensor and the Lanczos-van Stockum Einstein tensor, and none of them in a Weyl tensor.
