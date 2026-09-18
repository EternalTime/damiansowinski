# The Morris-Thorne traversable wormhole

This is the working behind the two coordinate systems in `MFS/assets/data/metrics/morris_thorne.json`.
Every number the entry prints is derived here, in order, from the line element down to the geodesic equations.
Nothing is left as an exercise and nothing is asserted that is not computed.

The companion script `_tools/derivations/verify_metrics.py` does the same work in sympy and compares it against the published file, so the algebra below is checkable by hand and by machine independently.

This entry is the first in the collection whose whole point is an inequality.
The geometry is written down first and the matter is read off it afterwards, which is the direction Morris and Thorne ran the problem in, and what comes back is a demand no ordinary matter can meet: at the throat the radial tension has to exceed the energy density.
Step 11 is where that inequality appears, Step 12 turns it into the failure of the null energy condition, and Step 16 shows it again in a chart where the throat is an ordinary point rather than a place two of the published components have to be read as limits.

The entry publishes the general case.
The redshift function $\Phi$ and the shape function $b$ are left as declared functions of the radial coordinate, so every component below is an identity in them, and the inequality of Step 11 is a statement about every wormhole of this form rather than about a chosen example.

---

## Step 1. Conventions

The signature is $(-,+,+,+)$ and the Christoffel symbols are those of the Levi-Civita connection,

$$\Gamma^\mu_{\nu\rho} = \tfrac{1}{2}g^{\mu\alpha}\left(\partial_\nu g_{\alpha\rho} + \partial_\rho g_{\alpha\nu} - \partial_\alpha g_{\nu\rho}\right),$$

with the lowered symbol $\Gamma_{\mu\nu\rho} = g_{\mu\alpha}\Gamma^\alpha_{\nu\rho}$.
The Riemann tensor is

$$R^\mu{}_{\nu\rho\sigma} = \partial_\rho \Gamma^\mu_{\nu\sigma} - \partial_\sigma \Gamma^\mu_{\nu\rho} + \Gamma^\mu_{\rho\lambda}\Gamma^\lambda_{\nu\sigma} - \Gamma^\mu_{\sigma\lambda}\Gamma^\lambda_{\nu\rho},$$

which is what the published Riemann components are in.

The Ricci tensor is the standard contraction on the first lower index,

$$R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu},$$

settled for the collection on 2026-09-18.
Here the choice matters in every slot, because this spacetime is not a vacuum and its Ricci tensor is nowhere zero.
Contracting on the last index instead would give the negative of every Ricci and Einstein component published, and it would flip the sign of the inequality Step 11 exists to state, so the convention is carried explicitly through Step 10 where the field equations are written down.
On the convention used here the field equations read

$$G_{\mu\nu} = \frac{8\pi G}{c^4}T_{\mu\nu},$$

and ordinary matter comes out with a positive energy density.

Factors of $G$ and $c$ are kept explicit.
The chart is the collection's, the one whose zeroth coordinate is

$$x^0 = ct.$$

The index is printed with the bare letter $t$, but the component printed against it is a component of that chart.
Everything from Step 3 onward is computed directly in it, so nothing has to be converted at the end.
Both radial coordinates used below, the areal radius $r$ and the proper distance $l$, are lengths already, so a radial derivative carries no factor of $c$ and no published component of either chart carries one either.

The dots in the geodesic equations are velocities of that same chart, so $\dot{t}$ means $d(ct)/d\lambda$ and not $dt/d\lambda$, and the equations are

$$\ddot{x}^\mu + \Gamma^\mu{}_{\nu\rho}\dot{x}^\nu\dot{x}^\rho = 0$$

with the same printed connection the entry lists.
Step 15 checks that this is the reading on which every term of every equation carries the dimensions of its left hand side.

Two abbreviations run through the areal chart below.

$$f(r) = 1 - \frac{b(r)}{r}, \qquad \Phi' = \partial_r\Phi, \qquad b' = \partial_r b.$$

They are shorthand for this document only.
The published file writes $1 - b/r$ out wherever it appears and spells every derivative $\partial_r\Phi$, $\partial_r^2\Phi$ or $\partial_r b$, because the checker reads the file symbol by symbol and has no way to be told what an abbreviation means.
One identity in that shorthand is used constantly and is worth recording once:

$$f' = \partial_r\left(1 - \frac{b}{r}\right) = \frac{b - rb'}{r^2}.$$

---

## Step 2. The line element, and what each function is for

The solution is

$$ds^2 = -e^{2\Phi(r)}c^2dt^2 + \frac{dr^2}{1 - \dfrac{b(r)}{r}} + r^2\left(d\theta^2 + \sin^2\theta\,d\phi^2\right).$$

The coordinate $r$ is the areal radius, so that the sphere labelled $r$ has area $4\pi r^2$ and circumference $2\pi r$.
It is not a distance from anything.
There is no center in this spacetime, and $r$ does not run down to zero: it runs from the throat radius $b_0$ outward, and the far side of the wormhole is a second copy of the same chart joined to this one at $r = b_0$.
That is the one thing about the areal chart a reader has to hold on to, and it is why Step 16 publishes a second chart in which the throat is an ordinary point and one coordinate covers both sides.

The two functions do two different jobs.

$\Phi(r)$ is the redshift function.
It fixes how much a signal climbing out of the wormhole is reddened, since a static clock at $r$ runs at $e^{\Phi(r)}$ times the rate of a clock at infinity.
Morris and Thorne's requirement on it is that it be finite everywhere: $e^{2\Phi}$ must never reach zero, because a zero of $g_{tt}$ is a horizon, and a horizon is a one way door, which is exactly what a traveller cannot be given [morris1988].
$\Phi$ is dimensionless, which it has to be to sit inside an exponential.

$b(r)$ is the shape function.
It fixes the spatial geometry alone, and it is a length beside $r$, which is what leaves $1 - b/r$ dimensionless.
Its name is earned in the embedding of Step 9: the shape of the funnel a slice of this spacetime makes in flat space depends on $b$ and not at all on $\Phi$.
It obeys $b(b_0) = b_0$ at the throat and $b(r) < r$ outside it, so that $f > 0$ away from the throat and the radial part of the metric stays positive.

At $r = b_0$ the metric component $g_{rr} = 1/f$ diverges.
That is a coordinate degeneracy and not a curvature singularity, in the same way that Schwarzschild's $r = r_s$ is, and Step 17 shows the curvature is perfectly finite there.
It does mean that several components published in this chart carry $r - b$ in a denominator, and that their value at the throat is the limit $r \to b_0$ rather than a value the chart holds.

---

## Step 3. The metric and its inverse

In the chart $x^0 = ct$, with the ordering $(t, r, \theta, \phi)$, the metric is diagonal:

$$g_{\mu\nu} = \operatorname{diag}\left(-e^{2\Phi},\ \frac{1}{f},\ r^2,\ r^2\sin^2\theta\right).$$

The $c^2$ of the line element is absorbed by the rescaling $t \mapsto ct$, which is why $g_{tt}$ is printed as $-e^{2\Phi}$ and carries no dimensions.
A diagonal metric inverts termwise:

$$g^{\mu\nu} = \operatorname{diag}\left(-e^{-2\Phi},\ f,\ \frac{1}{r^2},\ \frac{\csc^2\theta}{r^2}\right).$$

The determinant is $g = -e^{2\Phi}r^4\sin^2\theta/f$, which is negative wherever $f > 0$, so the signature is $(-,+,+,+)$ on the whole domain and the chart is Lorentzian everywhere except at the throat itself.

---

## Step 4. The Christoffel symbols

Everything depends on $r$ alone except the two $\theta$ dependences of the sphere, so only $\partial_r$ and $\partial_\theta$ can contribute.

$$\Gamma^t{}_{tr} = \tfrac{1}{2}g^{tt}\partial_r g_{tt} = \tfrac{1}{2}\left(-e^{-2\Phi}\right)\left(-2\Phi'e^{2\Phi}\right) = \Phi'.$$

$$\Gamma^r{}_{tt} = -\tfrac{1}{2}g^{rr}\partial_r g_{tt} = -\tfrac{1}{2}f\left(-2\Phi'e^{2\Phi}\right) = f\,e^{2\Phi}\Phi' = \left(1 - \frac{b}{r}\right)e^{2\Phi}\partial_r\Phi.$$

$$\Gamma^r{}_{rr} = \tfrac{1}{2}g^{rr}\partial_r g_{rr} = \tfrac{1}{2}f\,\partial_r\!\left(\frac{1}{f}\right) = -\frac{f'}{2f} = \frac{r\,\partial_r b - b}{2r\left(r - b\right)},$$

where the last step puts $f = (r-b)/r$ and $f' = (b - rb')/r^2$ in.

$$\Gamma^r{}_{\theta\theta} = -\tfrac{1}{2}g^{rr}\partial_r g_{\theta\theta} = -rf = -\left(r - b\right), \qquad \Gamma^r{}_{\phi\phi} = -\left(r - b\right)\sin^2\theta.$$

$$\Gamma^\theta{}_{r\theta} = \Gamma^\phi{}_{r\phi} = \tfrac{1}{2}g^{\theta\theta}\partial_r g_{\theta\theta} = \frac{1}{r}, \qquad \Gamma^\theta{}_{\phi\phi} = -\sin\theta\cos\theta, \qquad \Gamma^\phi{}_{\theta\phi} = \cot\theta.$$

Those are the thirteen nonzero symbols the entry publishes, counting the two orderings of each mixed pair separately.
The two of them that carry a time index, $\Gamma^t{}_{tr}$ and $\Gamma^r{}_{tt}$, illustrate the chart rule of Step 1.
$\Gamma^t{}_{tr}$ has one upper time index and one lower one, so the two factors of $c$ cancel and the chart symbol equals the bare one.
$\Gamma^r{}_{tt}$ has two lower time indices, so the chart symbol is $c^{-2}$ times the symbol computed with the bare $t$, and the $c^2$ that would otherwise sit in front of $f e^{2\Phi}\Phi'$ is gone.
Neither prints a $c$, which is the whole point of the convention.

Lowering the first index with the diagonal metric is a single multiplication each time:

$$\Gamma_{ttr} = g_{tt}\Gamma^t{}_{tr} = -e^{2\Phi}\partial_r\Phi, \qquad \Gamma_{rtt} = g_{rr}\Gamma^r{}_{tt} = \frac{1}{f}\cdot f e^{2\Phi}\Phi' = e^{2\Phi}\partial_r\Phi,$$

$$\Gamma_{rrr} = \frac{1}{f}\cdot\frac{r b' - b}{2r\left(r-b\right)} = \frac{r\,\partial_r b - b}{2\left(r - b\right)^2}, \qquad \Gamma_{r\theta\theta} = \frac{1}{f}\left(-\left(r-b\right)\right) = -r,$$

$$\Gamma_{\theta r\theta} = r, \qquad \Gamma_{\theta\phi\phi} = -r^2\sin\theta\cos\theta, \qquad \Gamma_{\phi r\phi} = r\sin^2\theta, \qquad \Gamma_{\phi\theta\phi} = r^2\sin\theta\cos\theta.$$

Note that $\Gamma_{rrr}$ is the one place where the factor $1/f$ turns a single power of $r - b$ in the denominator into a square.

---

## Step 5. The Riemann tensor, in the frame that makes it readable

The coordinate components of Riemann for this metric are six independent numbers dressed in metric factors.
The clean way to get them, and the way that pays for itself twice over in Steps 10 and 13, is to compute them once in the orthonormal frame of a static observer and then dress them.

Take the orthonormal coframe

$$\omega^{\hat{t}} = e^{\Phi}d(ct), \qquad \omega^{\hat{r}} = \frac{dr}{\sqrt{f}}, \qquad \omega^{\hat{\theta}} = r\,d\theta, \qquad \omega^{\hat{\phi}} = r\sin\theta\,d\phi,$$

whose dual frame is

$$e_{\hat{t}} = e^{-\Phi}\partial_{ct}, \qquad e_{\hat{r}} = \sqrt{f}\,\partial_r, \qquad e_{\hat{\theta}} = \frac{1}{r}\partial_\theta, \qquad e_{\hat{\phi}} = \frac{1}{r\sin\theta}\partial_\phi.$$

In that frame the metric is $\eta_{\hat{\mu}\hat{\nu}} = \operatorname{diag}(-1,1,1,1)$ and the curvature has exactly four independent components, which the spherical symmetry and the staticity force to be

$$E_1 \equiv R_{\hat{t}\hat{r}\hat{t}\hat{r}} = f\left(\Phi'' + \Phi'^2\right) + \tfrac{1}{2}f'\Phi',$$

$$E_2 \equiv R_{\hat{t}\hat{\theta}\hat{t}\hat{\theta}} = R_{\hat{t}\hat{\phi}\hat{t}\hat{\phi}} = \frac{f\Phi'}{r},$$

$$E_3 \equiv R_{\hat{r}\hat{\theta}\hat{r}\hat{\theta}} = R_{\hat{r}\hat{\phi}\hat{r}\hat{\phi}} = -\frac{f'}{2r} = \frac{r\,\partial_r b - b}{2r^3},$$

$$E_4 \equiv R_{\hat{\theta}\hat{\phi}\hat{\theta}\hat{\phi}} = \frac{1 - f}{r^2} = \frac{b}{r^3}.$$

$E_1$ is the tidal stretch along the radius, $E_2$ the tidal squeeze transverse to it, $E_3$ the curvature of the radial plane and $E_4$ the intrinsic curvature of the sphere, which is $1/r^2$ in flat space and is $b/r^3$ here because the sphere is embedded in a space that is not flat.

Two things are already visible.
$E_1$ and $E_2$ are the only pieces $\Phi$ enters at all; $E_3$ and $E_4$ know nothing about the redshift function.
And $E_2$ carries a bare factor of $f$, which vanishes at the throat, so the transverse tidal squeeze goes to zero there however $\Phi$ is chosen.
That second statement needs one caveat, made good in Step 11: $\partial_r\Phi$ itself diverges at the throat for any wormhole worth the name, and it is the absence of a horizon that keeps the divergence slow enough for the product $f\,\partial_r\Phi$ to vanish anyway.
$E_1$ does not vanish there.
Its three terms each blow up or vanish separately and the sum is finite, which Step 16 computes cleanly in a chart with no such cancellation in it.

Every component the entry publishes is one of these four with metric factors on it.
Converting a frame index back into a coordinate one multiplies by the corresponding coframe factor, so

$$R_{trtr} = \frac{e^{2\Phi}}{f}E_1 = \frac{e^{2\Phi}\left(2r\left(r-b\right)\left(\left(\partial_r\Phi\right)^2 + \partial_r^2\Phi\right) - \left(r\,\partial_r b - b\right)\partial_r\Phi\right)}{2r\left(r-b\right)},$$

$$R_{t\theta t\theta} = r^2 e^{2\Phi}E_2 = \left(r - b\right)e^{2\Phi}\partial_r\Phi, \qquad R_{t\phi t\phi} = R_{t\theta t\theta}\sin^2\theta,$$

$$R_{r\theta r\theta} = \frac{r^2}{f}E_3 = \frac{r\,\partial_r b - b}{2\left(r - b\right)}, \qquad R_{r\phi r\phi} = R_{r\theta r\theta}\sin^2\theta,$$

$$R_{\theta\phi\theta\phi} = r^4\sin^2\theta\,E_4 = b\,r\sin^2\theta.$$

The remaining published components are these with one index raised, or with a pair swapped, and the antisymmetry $R_{\mu\nu\rho\sigma} = -R_{\mu\nu\sigma\rho}$ supplies the sign.
Twenty four components in each variant, all of them one of four independent curvatures in a dress, which is why the entry's page shows the fully lowered variant on six lines and the mixed one on ten once the equal and opposite components are grouped together.

Everything that is not one of those six vanishes.
The reason is symmetry rather than algebra.
The spacetime is invariant under $t \mapsto -t$ and under reflection of either angle, and a Riemann component whose two index pairs are not the same pair, such as $R_{trt\theta}$, is odd under one of those reflections and so has to be zero.

---

## Step 6. The Ricci tensor

Contract on the first lower index.
In the orthonormal frame the contraction $R_{\hat{\mu}\hat{\nu}} = \eta^{\hat{\alpha}\hat{\beta}}R_{\hat{\beta}\hat{\mu}\hat{\alpha}\hat{\nu}}$ is a sum of at most three terms, and the $\eta^{\hat{t}\hat{t}} = -1$ is where the signature earns its keep:

$$R_{\hat{t}\hat{t}} = E_1 + 2E_2, \qquad R_{\hat{r}\hat{r}} = -E_1 + 2E_3, \qquad R_{\hat{\theta}\hat{\theta}} = R_{\hat{\phi}\hat{\phi}} = -E_2 + E_3 + E_4.$$

Dressing those back into the chart gives what the entry publishes:

$$R_{tt} = e^{2\Phi}\left(E_1 + 2E_2\right) = \frac{e^{2\Phi}\left(2r\left(r - b\right)\left(\left(\partial_r\Phi\right)^2 + \partial_r^2\Phi\right) + \left(4r - 3b - r\,\partial_r b\right)\partial_r\Phi\right)}{2r^2},$$

$$R_{rr} = \frac{-E_1 + 2E_3}{f} = \frac{\left(r\,\partial_r b - b\right)\left(2 + r\,\partial_r\Phi\right)}{2r^2\left(r - b\right)} - \left(\partial_r\Phi\right)^2 - \partial_r^2\Phi,$$

$$R_{\theta\theta} = r^2\left(-E_2 + E_3 + E_4\right) = \frac{b}{2r} + \frac{\partial_r b}{2} - \left(r - b\right)\partial_r\Phi, \qquad R_{\phi\phi} = R_{\theta\theta}\sin^2\theta.$$

The off diagonal components vanish, by the same reflection argument as in Step 5.
None of the four diagonal ones vanishes for a general pair of functions, which is the plainest statement that this spacetime is not a vacuum: something has to be there.

---

## Step 7. The Ricci scalar

$$R = \eta^{\hat{\mu}\hat{\nu}}R_{\hat{\mu}\hat{\nu}} = -R_{\hat{t}\hat{t}} + R_{\hat{r}\hat{r}} + 2R_{\hat{\theta}\hat{\theta}} = -2E_1 - 4E_2 + 4E_3 + 2E_4,$$

which written out is

$$R = \frac{2\partial_r b + \left(3b + r\,\partial_r b - 4r\right)\partial_r\Phi - 2r\left(r - b\right)\left(\left(\partial_r\Phi\right)^2 + \partial_r^2\Phi\right)}{r^2}.$$

---

## Step 8. The Einstein tensor

$G_{\hat{\mu}\hat{\nu}} = R_{\hat{\mu}\hat{\nu}} - \tfrac{1}{2}\eta_{\hat{\mu}\hat{\nu}}R$ with $R/2 = -E_1 - 2E_2 + 2E_3 + E_4$ collapses to something much simpler than the Ricci tensor it came from:

$$G_{\hat{t}\hat{t}} = 2E_3 + E_4 = \frac{r\,\partial_r b - b}{r^3} + \frac{b}{r^3} = \frac{\partial_r b}{r^2},$$

$$G_{\hat{r}\hat{r}} = 2E_2 - E_4 = \frac{2f\,\partial_r\Phi}{r} - \frac{b}{r^3},$$

$$G_{\hat{\theta}\hat{\theta}} = G_{\hat{\phi}\hat{\phi}} = E_1 + E_2 - E_3.$$

The first of those is the single cleanest formula in the entry.
Every trace of $\Phi$ has cancelled out of it, and what is left says that the energy density is fixed by the slope of the shape function and by nothing else.
The second has lost its dependence on $\Phi''$.
Only the transverse pressure still carries the second derivative of the redshift function.

Dressed into the chart:

$$G_{tt} = e^{2\Phi}G_{\hat{t}\hat{t}} = \frac{e^{2\Phi}\,\partial_r b}{r^2}, \qquad G_{rr} = \frac{G_{\hat{r}\hat{r}}}{f} = \frac{2r\left(r - b\right)\partial_r\Phi - b}{r^2\left(r - b\right)},$$

$$G_{\theta\theta} = r^2\left(E_1 + E_2 - E_3\right) = r\left(r - b\right)\left(\left(\partial_r\Phi\right)^2 + \partial_r^2\Phi\right) + \frac{\left(2r - b - r\,\partial_r b\right)\partial_r\Phi}{2} + \frac{b}{2r} - \frac{\partial_r b}{2},$$

with $G_{\phi\phi} = G_{\theta\theta}\sin^2\theta$.
The entry also publishes the mixed and fully raised variants, of which the one worth quoting is

$$G^t{}_t = -\frac{\partial_r b}{r^2},$$

negative wherever the shape function is rising, which on the sign convention of Step 1 is $-8\pi G\rho/c^2$ and so says the same thing as $G_{tt}$.

---

## Step 9. The embedding, and where the throat is

Take a moment of time and the equatorial plane, $t$ constant and $\theta = \pi/2$.
What is left is a two dimensional surface with

$$ds^2 = \frac{dr^2}{1 - \dfrac{b}{r}} + r^2d\phi^2,$$

and the question of what it looks like is the question of embedding it in flat three dimensional space, $ds^2 = dz^2 + dr^2 + r^2d\phi^2$, as a surface of revolution $z = z(r)$.
Matching the two gives

$$\left(\frac{dz}{dr}\right)^2 + 1 = \frac{1}{1 - \dfrac{b}{r}}, \qquad \text{so} \qquad \frac{dz}{dr} = \pm\left(\frac{r}{b(r)} - 1\right)^{-1/2}.$$

This is where the shape function gets its name: the profile of the funnel depends on $b$ and not at all on $\Phi$.
At $r = b_0$, where $b = r$, the slope $dz/dr$ is infinite: the funnel is vertical there, which is the throat.
Away from it $r$ increases and the surface opens out.

For the surface to be a wormhole rather than a spike it has to flare outward at the throat, which is the statement that $r$ is a minimum as a function of the embedding height $z$:

$$\frac{d^2r}{dz^2} = \frac{b - r\,\partial_r b}{2b^2} > 0 \qquad \text{at } r = b_0.$$

With $b(b_0) = b_0$ this is exactly

$$\partial_r b(b_0) < 1,$$

the flare out condition [morris1988].
It is published in the entry as part of the description of $b$, since it is a constraint on the function rather than an equation between published components, and the checker has no inequalities to test it with.
Step 16 shows it again in a form that needs no embedding at all.

---

## Step 10. What the geometry demands of its matter

Read the field equations of Step 1 in the static orthonormal frame, where the stress energy of a static, spherically symmetric source is diagonal:

$$T_{\hat{t}\hat{t}} = \rho c^2, \qquad T_{\hat{r}\hat{r}} = p_r = -\tau, \qquad T_{\hat{\theta}\hat{\theta}} = T_{\hat{\phi}\hat{\phi}} = p.$$

$\rho$ is the mass density a static observer measures, $\tau$ is the radial tension, which is minus the radial pressure, and $p$ is the transverse pressure.
Tension rather than pressure is the natural variable here, and Step 11 is why.

With Step 8 in hand the three equations are immediate:

$$\rho = \frac{c^2}{8\pi G}\frac{\partial_r b}{r^2}, \qquad \tau = \frac{c^4}{8\pi G}\left(\frac{b}{r^3} - \frac{2\left(1 - \dfrac{b}{r}\right)\partial_r\Phi}{r}\right), \qquad p = \frac{c^4}{8\pi G}\left(E_1 + E_2 - E_3\right).$$

These are Morris and Thorne's equations (12), with $G$ and $c$ kept rather than set to one [morris1988].

Two readings are worth making before the inequality.
The density is fixed by $b$ alone, so a shape function that is flat, $\partial_r b = 0$, describes a geometry threaded by matter of zero density, which is what the Schwarzschild limit of Step 14 turns out to be.
The redshift function never appears in the density at all.
It appears in the tension only through the combination $f\,\partial_r\Phi$, and that combination is the subject of the first half of Step 11.

---

## Step 11. The inequality

Evaluate the two of those at $r = b_0$, where $b(b_0) = b_0$ and therefore $f(b_0) = 0$.

The tension loses its $\Phi$ term there, and it is worth being careful about why, because the factor $f$ is multiplying something that is not finite.
In the proper distance of Step 16 the throat sits at $l = 0$ with $r(l) = b_0 + \tfrac{1}{2}\partial_l^2r(0)\,l^2 + \dots$, so $\partial_l r$ vanishes linearly in $l$ and

$$f = \left(\partial_l r\right)^2 \sim \left(\partial_l^2r(0)\right)^2l^2, \qquad \partial_r\Phi = \frac{\partial_l\Phi}{\partial_l r} \sim \frac{\partial_l\Phi(0)}{\partial_l^2r(0)\,l}.$$

$\partial_r\Phi$ does diverge at the throat, but only as $1/l$, so $f\,\partial_r\Phi$ vanishes as $l$.
What bounds $\partial_l\Phi(0)$ is the no horizon requirement of Step 2: $\Phi$ is finite at the throat, so it is a smooth function of the proper distance there and its first derivative is a number.
A $\Phi$ that was allowed a horizon could diverge fast enough for the product to survive, and Step 14 shows the one case where it does.

So, for a wormhole with no horizon,

$$\tau_0 = \frac{c^4}{8\pi G}\frac{b_0}{b_0^3} = \frac{c^4}{8\pi Gb_0^2}.$$

This is a number.
It does not depend on the redshift function, and it does not depend on the shape function beyond the throat radius itself.
Every wormhole of this form with a throat of radius $b_0$ has exactly this radial tension at the throat.

The density there is

$$\rho_0c^2 = \frac{c^4}{8\pi Gb_0^2}\partial_r b(b_0),$$

so the difference is

$$\tau_0 - \rho_0c^2 = \frac{c^4}{8\pi Gb_0^2}\left(1 - \partial_r b(b_0)\right).$$

The flare out condition of Step 9 is $\partial_r b(b_0) < 1$, so the bracket is positive and

$$\boxed{\ \tau_0 > \rho_0c^2\ }$$

at the throat of every traversable wormhole of this form.

The radial tension beats the energy density.
Morris and Thorne put a number on it for a throat of three kilometers: about $10^{37}$ dynes per square centimeter, which is the order of the pressure at the center of the heaviest neutron stars, and which no material anybody has a theory of can supply as a tension [morris1988].
They called the stuff exotic.

The inequality is not an artefact of the ansatz.
It is the flare out condition, rearranged.
Anything that makes the funnel of Step 9 open out instead of pinching to a point makes the bracket positive, and the bracket positive is the inequality.

---

## Step 12. The null energy condition, and why it fails

The null energy condition asks that $T_{\mu\nu}k^\mu k^\nu \ge 0$ for every null $k^\mu$.
It is the weakest of the standard energy conditions, implied by all the others, and it is the one the singularity theorems and topological censorship are built on.

Take the radial null vector, which in the orthonormal frame is $k^{\hat{\mu}} = (1,1,0,0)$: a light ray crossing the throat straight through.
Then

$$T_{\hat{\mu}\hat{\nu}}k^{\hat{\mu}}k^{\hat{\nu}} = T_{\hat{t}\hat{t}} + T_{\hat{r}\hat{r}} = \rho c^2 + p_r = \rho c^2 - \tau,$$

which at the throat is $\rho_0c^2 - \tau_0 < 0$ by Step 11.
The same thing read geometrically, without ever naming a stress tensor, is

$$G_{\hat{\mu}\hat{\nu}}k^{\hat{\mu}}k^{\hat{\nu}} = G_{\hat{t}\hat{t}} + G_{\hat{r}\hat{r}} = \left(2E_3 + E_4\right) + \left(2E_2 - E_4\right) = 2\left(E_2 + E_3\right) = \frac{2f\,\partial_r\Phi}{r} - \frac{f'}{r},$$

where the $E_4$ terms cancelled, which is the algebraic reason the sphere's own curvature never enters.
At the throat the first term vanishes, by the argument of Step 11, and $f' = \left(1 - \partial_r b(b_0)\right)/b_0$, so

$$G_{\hat{\mu}\hat{\nu}}k^{\hat{\mu}}k^{\hat{\nu}}\Big|_{r=b_0} = -\frac{1 - \partial_r b(b_0)}{b_0^2} < 0.$$

So an observer riding that light ray measures a negative energy density.
The failure is not marginal and it is not tunable: it is forced by the same inequality that makes the thing a wormhole at all, and it holds however $\Phi$ is chosen, since the redshift function leaves no trace at the throat once the no horizon requirement is in force.

This is the local version of a theorem.
Friedman, Schleich and Witt's topological censorship says that under the null energy condition any topological shortcut collapses too fast for light to get through it [friedman1993].
The calculation above is that theorem's contrapositive in one spacetime: a shortcut light can cross is a shortcut where the condition fails.

---

## Step 13. The Weyl tensor

The Weyl tensor is Riemann with its traces removed,

$$C_{\mu\nu\rho\sigma} = R_{\mu\nu\rho\sigma} - \left(g_{\mu[\rho}R_{\sigma]\nu} - g_{\nu[\rho}R_{\sigma]\mu}\right) + \tfrac{1}{3}R\,g_{\mu[\rho}g_{\sigma]\nu},$$

built from the same contraction as everywhere else.
This spacetime is not a vacuum, so the Weyl tensor is not a copy of the Riemann tensor, and it has to be computed.

In the orthonormal frame the subtraction is arithmetic on the four numbers of Step 5.
The result is a single scalar,

$$\Psi = \tfrac{1}{3}\left(E_2 + E_4 - E_1 - E_3\right),$$

with the whole tensor being

$$C_{\hat{t}\hat{r}\hat{t}\hat{r}} = -\Psi, \qquad C_{\hat{t}\hat{\theta}\hat{t}\hat{\theta}} = C_{\hat{t}\hat{\phi}\hat{t}\hat{\phi}} = \tfrac{1}{2}\Psi, \qquad C_{\hat{r}\hat{\theta}\hat{r}\hat{\theta}} = C_{\hat{r}\hat{\phi}\hat{r}\hat{\phi}} = -\tfrac{1}{2}\Psi, \qquad C_{\hat{\theta}\hat{\phi}\hat{\theta}\hat{\phi}} = \Psi.$$

The pattern $-1, \tfrac{1}{2}, \tfrac{1}{2}, -\tfrac{1}{2}, -\tfrac{1}{2}, 1$ sums to zero as it must, and it is the signature of Petrov type D, which is what every static spherically symmetric spacetime is away from the points where $\Psi$ vanishes.

Written out, $6r^3\Psi$ is the bracket every published Weyl component is a multiple of:

$$6r^3\Psi = 3b - r\,\partial_r b + r\left(2r - 3b + r\,\partial_r b\right)\partial_r\Phi - 2r^2\left(r - b\right)\left(\left(\partial_r\Phi\right)^2 + \partial_r^2\Phi\right).$$

So, for instance,

$$C^t{}_{rtr} = \frac{6r^3\Psi}{6r^2\left(r-b\right)}, \qquad C_{\theta\phi\theta\phi} = \frac{r\sin^2\theta\cdot 6r^3\Psi}{6} = r^4\sin^2\theta\,\Psi,$$

and the other twenty two components of each variant are this one with metric factors and a sign.
The Weyl tensor vanishes only where $\Psi$ does, which is the condition for the geometry to be conformally flat, and that is a differential equation on $\Phi$ and $b$ rather than something a general pair of them satisfies.

---

## Step 14. The Kretschmann scalar

For a Riemann tensor whose only nonzero orthonormal components are the six pairs of Step 5, every term of $R_{\mu\nu\rho\sigma}R^{\mu\nu\rho\sigma}$ is one of those six squared.
Each unordered index pair contributes four terms, one for each ordering of the two pairs, and raising indices with $\eta$ contributes an even number of minus signs each time, so

$$K = 4\left(E_1^2 + 2E_2^2 + 2E_3^2 + E_4^2\right).$$

It is a sum of squares, so it is positive wherever the geometry is curved at all, and it can only be made to vanish by making all four of the $E$ vanish together, which is flat space.
Written out in the published functions,

$$K = 4\left(\frac{\left(r - b\right)\left(\left(\partial_r\Phi\right)^2 + \partial_r^2\Phi\right)}{r} + \frac{\left(b - r\,\partial_r b\right)\partial_r\Phi}{2r^2}\right)^2 + \frac{8\left(r-b\right)^2\left(\partial_r\Phi\right)^2}{r^4} + \frac{2\left(r\,\partial_r b - b\right)^2}{r^6} + \frac{4b^2}{r^6}.$$

Nothing in it diverges at the throat, though the areal chart does not show that as plainly as one would like.
The $E_2$ term vanishes there, and the $E_3$ and $E_4$ terms are manifestly finite.
$E_1$ is the one to watch.
It carries a factor $f$ that vanishes at the throat against a factor $\partial_r\Phi$ that in general diverges there, because a redshift function that is smooth in the proper distance of Step 16 has $\partial_r\Phi = \partial_l\Phi/\partial_l r$ with a denominator going to zero.
The product is finite, and Step 16 computes it with no indeterminate form anywhere, giving $E_1(b_0) = \left(\left(\partial_l\Phi\right)^2 + \partial_l^2\Phi\right)\big|_{l=0}$.
With that,

$$K\big|_{r=b_0} = 4E_1(b_0)^2 + \frac{2\left(1 - \partial_r b(b_0)\right)^2 + 4}{b_0^4},$$

which is positive and bounded, so the degeneracy of $g_{rr}$ at the throat is a fault of the chart and not of the spacetime.
This is the sense in which a Morris-Thorne wormhole has no singularity to hide: both $\Phi$ and $b$ are assumed finite and smooth as functions of the proper distance on the whole domain, and there is nowhere for a curvature blow up to come from.

The check on all of this is the Schwarzschild limit.
Put $b = r_s$ constant, so $\partial_r b = 0$, and $e^{2\Phi} = 1 - r_s/r = f$.
Then $\partial_r\Phi = r_s/(2r^2f)$, and the four curvatures come out as $E_1 = -r_s/r^3$, $E_2 = r_s/(2r^3)$, $E_3 = -r_s/(2r^3)$ and $E_4 = r_s/r^3$.
So $G_{\hat{t}\hat{t}} = 2E_3 + E_4 = 0$, the solution is a vacuum as it must be, while

$$K = 4\left(\frac{r_s^2}{r^6} + \frac{2r_s^2}{4r^6} + \frac{2r_s^2}{4r^6} + \frac{r_s^2}{r^6}\right) = \frac{12r_s^2}{r^6},$$

which is the Schwarzschild value the collection publishes elsewhere.
$\Psi = \tfrac{1}{3}(E_2 + E_4 - E_1 - E_3) = r_s/r^3$, which is the Schwarzschild Weyl scalar, and there the Weyl tensor is equal to the Riemann tensor, as it must be in a vacuum.

This limit is also the exception that proves Step 11.
A constant shape function passes the flare out condition, $\partial_r b(b_0) = 0 < 1$, and it puts $b_0 = r_s$, so it looks like a wormhole and it is the Einstein-Rosen bridge.
But $e^{2\Phi} = f$ vanishes at $r = r_s$, which is a horizon, so $\partial_r\Phi = r_s/(2r^2f)$ diverges as $1/f$ rather than as $1/\sqrt{f}$ and the product $f\,\partial_r\Phi$ survives at the throat instead of vanishing.
It survives at exactly the value that cancels the other term: $\tau_0 = \dfrac{c^4}{8\pi G}\left(\dfrac{r_s}{r_s^3} - \dfrac{r_s}{r_s^3}\right) = 0$, and the vacuum needs no exotic matter.
What it has instead is the horizon, and a horizon is a one way door.
That is the trade the whole subject turns on, and it is why Morris and Thorne put the no horizon requirement first.

---

## Step 15. The geodesic equations, and the dimensional check

With the connection of Step 4 the equations $\ddot{x}^\mu + \Gamma^\mu{}_{\nu\rho}\dot{x}^\nu\dot{x}^\rho = 0$ read

$$\ddot{t} + 2\partial_r\Phi\,\dot{r}\dot{t} = 0,$$

$$\ddot{r} + \left(1 - \frac{b}{r}\right)e^{2\Phi}\partial_r\Phi\,\dot{t}^2 + \frac{r\,\partial_r b - b}{2r\left(r - b\right)}\dot{r}^2 - \left(r - b\right)\dot{\theta}^2 - \left(r - b\right)\sin^2\theta\,\dot{\phi}^2 = 0,$$

$$\ddot{\theta} + \frac{2}{r}\dot{r}\dot{\theta} - \sin\theta\cos\theta\,\dot{\phi}^2 = 0, \qquad \ddot{\phi} + \frac{2}{r}\dot{r}\dot{\phi} + 2\cot\theta\,\dot{\theta}\dot{\phi} = 0,$$

with the dots as declared in Step 1.
The first integrates once to $e^{2\Phi}\dot{t} = \text{constant}$, which is the conserved energy of the static Killing vector, and the last to $r^2\sin^2\theta\,\dot{\phi} = \text{constant}$, the angular momentum.

The dimensional check is the one that pins the reading of the dots down.
Every term of an equation is measured against its own second derivative, so the terms of the $r$ equation all have to carry $L/\lambda^2$ with $\lambda$ the affine parameter.

$\ddot{r}$ carries $L/\lambda^2$ by definition.

$\left(1 - b/r\right)e^{2\Phi}\partial_r\Phi\,\dot{t}^2$: the bracket and the exponential are dimensionless, $\partial_r\Phi$ carries $1/L$ because $\Phi$ is dimensionless and $r$ is a length, and $\dot{t} = d(ct)/d\lambda$ carries $L/\lambda$.
The product is $\tfrac{1}{L}\cdot\tfrac{L^2}{\lambda^2} = \tfrac{L}{\lambda^2}$.

Read the other way, with $\dot{t} = dt/d\lambda$ carrying $T/\lambda$, that term would come out short of $\ddot{r}$ by a factor of $c^2$.
That is the failure the checker's dimensional pass exists to catch, and it is why the chart convention has to be stated rather than assumed.

$\dfrac{r\,\partial_r b - b}{2r\left(r-b\right)}\dot{r}^2$: the fraction is a length over a length squared, so $1/L$, times $L^2/\lambda^2$.

$\left(r - b\right)\dot{\theta}^2$: a length times $1/\lambda^2$, since an angle is dimensionless and so $\dot\theta$ carries $1/\lambda$.

All four carry $L/\lambda^2$.

---

## Step 16. The proper distance chart

The areal chart has two faults, both of them the same fault.
It is degenerate at $r = b_0$, and it covers one side of the wormhole only.
Both are cured by using the proper radial distance as the coordinate instead.

Define $l$ by

$$dl = \pm\frac{dr}{\sqrt{1 - \dfrac{b(r)}{r}}},$$

with $l = 0$ at the throat, the upper sign on one side and the lower on the other.
$l$ is the distance a traveller with a ruler actually measures, it runs over the whole real line, and it passes through the throat without noticing anything.
The line element becomes

$$ds^2 = -e^{2\Phi(l)}c^2dt^2 + dl^2 + r(l)^2\left(d\theta^2 + \sin^2\theta\,d\phi^2\right),$$

with the areal radius $r(l)$ now a declared function of the coordinate rather than the coordinate itself.
It has a minimum $r(0) = b_0$ at the throat and grows without bound on both sides.

This is the form the Ellis-Bronnikov drainhole is published in, and that entry is the special case

$$\Phi = 0, \qquad r(l) = \sqrt{l^2 + \ell^2}, \qquad b_0 = \ell,$$

which makes it the best available check on this chart, since it was computed independently and years earlier.
Substituting those two functions into the ninety nine values `ellis_bronnikov.json` publishes and that this chart publishes as well reproduces every one of them exactly.
Three by hand:

$$\partial_l r = \frac{l}{\sqrt{l^2+\ell^2}}, \qquad \partial_l^2 r = \frac{\ell^2}{\left(l^2+\ell^2\right)^{3/2}}, \qquad 1 - \left(\partial_l r\right)^2 = \frac{\ell^2}{l^2+\ell^2},$$

so $R_{ll} = -2\partial_l^2r/r = -2\ell^2/\left(l^2+\ell^2\right)^2$, which is that entry's one nonzero Ricci component; $G_{tt} = \left(1 - \left(\partial_l r\right)^2 - 2r\partial_l^2 r\right)/r^2 = -\ell^2/\left(l^2+\ell^2\right)^2$, which is the negative energy density it names in its own convention field; and $K = 8\left(\partial_l^2r\right)^2/r^2 + 4\left(1 - \left(\partial_l r\right)^2\right)^2/r^4 = 12\ell^4/\left(l^2+\ell^2\right)^4$, which is its Kretschmann scalar.
The differences that remain are that entry's and not this chart's: it publishes no inverse metric, and it omits twelve nonzero Riemann and Weyl components that this chart carries.
Those twelve are part of the disagreements the checker already reports against that file, and they are left alone here.

The dictionary between the charts is one identity and its derivative.
Squaring the definition of $l$,

$$\left(\partial_l r\right)^2 = 1 - \frac{b}{r}, \qquad \text{so} \qquad b = r\left(1 - \left(\partial_l r\right)^2\right),$$

and differentiating that along $l$ gives, wherever $\partial_l r \ne 0$,

$$\partial_l^2 r = \frac{b - r\,\partial_r b}{2r^2}.$$

At the throat $\partial_l r = 0$ and the limit of the right hand side is

$$\partial_l^2 r(0) = \frac{1 - \partial_r b(b_0)}{2b_0}.$$

So the flare out condition $\partial_r b(b_0) < 1$ is exactly

$$\partial_l^2 r(0) > 0,$$

which is the statement that the areal radius really has a minimum at the throat rather than an inflection.
That is the whole of Step 9 without an embedding diagram, and it is the form the entry publishes in this chart's parameter description.

Everything else is easier here than in the areal chart, because the metric is polynomial in $r$ and its derivatives with no $1/f$ anywhere.
The four curvatures of Step 5 become

$$E_1 = \left(\partial_l\Phi\right)^2 + \partial_l^2\Phi, \qquad E_2 = \frac{\partial_l\Phi\,\partial_l r}{r}, \qquad E_3 = -\frac{\partial_l^2 r}{r}, \qquad E_4 = \frac{1 - \left(\partial_l r\right)^2}{r^2},$$

which agree with the areal expressions under the dictionary, by the chain rule $\partial_l = \partial_l r\,\partial_r$.
The Einstein tensor is

$$G_{\hat{t}\hat{t}} = \frac{1 - \left(\partial_l r\right)^2 - 2r\,\partial_l^2 r}{r^2}, \qquad G_{\hat{r}\hat{r}} = \frac{\left(\partial_l r\right)^2 - 1 + 2r\,\partial_l\Phi\,\partial_l r}{r^2},$$

$$G_{\hat{\theta}\hat{\theta}} = G_{\hat{\phi}\hat{\phi}} = \left(\partial_l\Phi\right)^2 + \partial_l^2\Phi + \frac{\partial_l\Phi\,\partial_l r + \partial_l^2 r}{r},$$

and the null energy condition test of Step 12 comes out in one line, everywhere and not only at the throat:

$$G_{\hat{\mu}\hat{\nu}}k^{\hat{\mu}}k^{\hat{\nu}} = G_{\hat{t}\hat{t}} + G_{\hat{r}\hat{r}} = \frac{2\left(\partial_l\Phi\,\partial_l r - \partial_l^2 r\right)}{r}.$$

At the throat $\partial_l r = 0$, so the first term drops and

$$G_{\hat{\mu}\hat{\nu}}k^{\hat{\mu}}k^{\hat{\nu}}\Big|_{l=0} = -\frac{2\partial_l^2 r(0)}{b_0} = -\frac{1 - \partial_r b(b_0)}{b_0^2},$$

which is the areal chart's answer, reached without taking a limit and without the redshift function ever being touched.
That agreement between two independently computed charts is the strongest check in this document.

The Weyl scalar in this chart is

$$3r^2\Psi = -\left(r^2\left(\left(\partial_l\Phi\right)^2 + \partial_l^2\Phi\right) - r\,\partial_l\Phi\,\partial_l r + \left(\partial_l r\right)^2 - r\,\partial_l^2 r - 1\right),$$

the bracket every published Weyl component of this chart is a multiple of, and the Kretschmann scalar is the same sum of four squares,

$$K = 4\left(\left(\partial_l\Phi\right)^2 + \partial_l^2\Phi\right)^2 + \frac{8\left(\partial_l\Phi\right)^2\left(\partial_l r\right)^2}{r^2} + \frac{8\left(\partial_l^2 r\right)^2}{r^2} + \frac{4\left(1 - \left(\partial_l r\right)^2\right)^2}{r^4}.$$

Both are manifestly finite at $l = 0$, with nothing to cancel and no limit to take.

---

## Step 17. What the entry publishes, and what the checker needs

For each of the two charts the entry publishes the line element, the metric and its inverse, both Christoffel variants, both Riemann variants, the three Ricci variants, the Ricci scalar, the Kretschmann scalar, the three Einstein variants, both Weyl variants and the four geodesic equations.

No block is empty and nothing vanishes identically.
Within each block the components that are zero are zero for one reason, given in Step 5: the metric is diagonal and depends on the radial coordinate alone, and the spacetime is invariant under reversal of $t$ and under reflection of either angle, so every off diagonal Ricci and Einstein component vanishes and the only surviving Riemann and Weyl components are those whose two index pairs are the same pair.

`verify_metrics.py` needed one declaration per chart.
`DIMENSIONS` gained an entry for each, declaring $t$ as a time, so the checker knows the chart multiplies it by $c$, and declaring $\Phi$ dimensionless, $b$ a length and $b_0$ a length in the areal chart, and $\Phi$ dimensionless and $r$ a length in the proper distance chart, where $r$ is a declared function rather than a coordinate.

The two functions are declared to the reader as `\Phi = \Phi(r)` and `b = b(r)`, and every published derivative of them is spelled `\partial_r\Phi`, `\partial_r^2\Phi` or `\partial_r b` rather than with a prime.
A function of one coordinate does answer to a prime, and `b'` reads correctly, but a function whose name is a LaTeX command does not: the reader turns primes into suffixes before it turns `\Phi` into a name, so `\Phi'` leaves a stray backslash behind and is rejected as unhandled LaTeX.
The `\partial` spelling works for either kind of name, and it is what this entry uses throughout so that the two functions are written the same way.
Both radial coordinates are lengths, so none of those derivatives carries a factor of $c$, unlike the retarded time derivatives of the Vaidya and pp-wave entries.

No entry in `PARAMETER_RELATIONS` is needed.
$\Phi$ and $b$ are free: every published component is an identity in them, and the entry claims nothing that holds only on some surface in the space of functions.
The two constraints that do exist, $b(b_0) = b_0$ and $\partial_r b(b_0) < 1$, are a boundary condition and an inequality.
Neither is a relation the checker's machinery can carry, which takes rational parametrisations of equalities between free parameters, and neither is needed for any published value to be correct.
They are published in the parameter descriptions instead, where a reader meets them next to the function they constrain.
