# The Tolman-Bondi inhomogeneous dust universe

This is the working behind the comoving synchronous coordinate system in `MFS/assets/data/metrics/tolman_bondi.json`.
Every number the entry prints is derived here, in order, from the line element down to the geodesic equations.
Nothing is left as an exercise and nothing is asserted that is not computed.

The companion script `_tools/derivations/verify_metrics.py` does the same work in sympy and compares it against the published file, so the algebra below is checkable by hand and by machine independently.
Step 18 says what that comparison covers, and how the third derivative of Step 6 came to be read by it.

Like FRW and like Bianchi type I, this entry's metric is not a formula.
Its areal radius $R(r,t)$ and its energy function $E(r)$ are undetermined, and no field equation is imposed on them anywhere in this file.
That is the whole point here, more than it is for either of those two.
The solution is named for a dust universe, and a reader is entitled to ask where the dust is; the answer is that it is one equation away, and the entry publishes the geometry on both sides of that equation rather than only on the far side of it.
Step 10 reads a density and two pressures off the Einstein tensor, Step 11 sets the pressures to zero and gets the evolution equation and the mass function $M(r)$ out of them, Step 12 turns the time component into $\rho = M'/(4\pi R^2\partial_r R)$, and Step 13 sets $R = a(t)r$ and recovers the FRW entry component for component.

---

## Step 1. Conventions

The signature is $(-,+,+,+)$ and the Christoffel symbols are those of the Levi-Civita connection,

$$\Gamma^\mu{}_{\nu\rho} = \tfrac{1}{2}g^{\mu\alpha}\left(\partial_\nu g_{\alpha\rho} + \partial_\rho g_{\alpha\nu} - \partial_\alpha g_{\nu\rho}\right),$$

with the lowered symbol $\Gamma_{\mu\nu\rho} = g_{\mu\alpha}\Gamma^\alpha{}_{\nu\rho}$.
The Riemann tensor is

$$R^\mu{}_{\nu\rho\sigma} = \partial_\rho \Gamma^\mu{}_{\nu\sigma} - \partial_\sigma \Gamma^\mu{}_{\nu\rho} + \Gamma^\mu{}_{\rho\lambda}\Gamma^\lambda{}_{\nu\sigma} - \Gamma^\mu{}_{\sigma\lambda}\Gamma^\lambda{}_{\nu\rho},$$

which is what the published Riemann components are in.

The Ricci tensor is the standard contraction,

$$R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu},$$

settled for the collection on 18 September 2026.
It is the contraction the signature asks for, in that it is the one that gives ordinary matter a positive energy density, and Step 10 is where that matters here: the published $G_{tt}$ is $8\pi G\rho/c^2$ with $\rho$ positive, and on the other contraction every sign in this file would turn over and the dust would weigh a negative amount.

The Weyl tensor is

$$C_{\mu\nu\rho\sigma} = R_{\mu\nu\rho\sigma} - \tfrac{1}{2}\left(g_{\mu\rho}R_{\sigma\nu} - g_{\mu\sigma}R_{\rho\nu} - g_{\nu\rho}R_{\sigma\mu} + g_{\nu\sigma}R_{\rho\mu}\right) + \tfrac{1}{6}R\left(g_{\mu\rho}g_{\sigma\nu} - g_{\mu\sigma}g_{\rho\nu}\right),$$

built from the same contraction, and Step 14 computes it from that definition rather than copying the Riemann tensor across.
This spacetime is not a vacuum, so the two are not equal, and the difference is the physical content of Step 14.

Factors of $c$ and $G$ are kept explicit.
The chart, here and everywhere else in the collection, is the one whose time coordinate is $x^0 = ct$.
Schwarzschild shows the convention plainly, quoting $g_{tt} = -(1-r_s/r)$ against a line element whose time term is $-(1-r_s/r)c^2dt^2$.
Since $t$ here is a time, the chart coordinate is $ct$, and every component printed against an index $t$ is a component in that chart even though the index is written with the bare letter.

Because the rescaling $t \to ct$ is linear with constant coefficients, the rule is arithmetic: a component in the chart is the component taken with the bare coordinate, multiplied by $c$ once for every upper $t$ index and divided by $c$ once for every lower one.
The Christoffel symbols obey the same rule as the tensors, since the inhomogeneous term in their transformation law carries a second derivative of the coordinate change, which vanishes for a linear one.
A scalar has no index at all and is therefore the same number in both charts.

This entry does the whole computation in the chart from the start, as Bianchi type I and Alcubierre do, because its metric functions are undetermined and nothing is lost by differentiating with respect to $x^0$ directly.
A printed $\partial_t$ is that derivative and no other:

$$\partial_t \equiv \frac{\partial}{\partial(ct)} = \frac{1}{c}\frac{\partial}{\partial t}, \qquad \partial_r \equiv \frac{\partial}{\partial r},$$

so that $\partial_t R$ and $\partial_r R$ are both dimensionless, $\partial_t^2R$, $\partial_r^2R$ and $\partial_r\partial_t R$ each carry an inverse length, and $\partial_r\partial_t^2R$ carries an inverse length squared.
Every factor of $c$ the chart demands is in that one definition, and apart from the $-c^2dt^2$ of the line element no factor of $c$ appears anywhere in the entry.
The energy function is a function of $r$ alone, so its derivative is written with a prime, $E' = dE/dr$.
The relation to the notation a cosmologist uses is $\dot{R} = c\,\partial_t R$, so a printed $\partial_t R$ of $10^{-3}$ is an expansion velocity of $3\times10^5$ metres per second.

The dots in the geodesic equations are velocities of that same chart with respect to an affine parameter, so $\dot{t}$ means $d(ct)/d\lambda$ and not $dt/d\lambda$, and the equations are

$$\ddot{x}^\mu + \Gamma^\mu{}_{\nu\rho}\dot{x}^\nu\dot{x}^\rho = 0$$

with the same printed $\Gamma$ the entry lists above them.
Step 17 checks that reading term by term.

Three abbreviations are used below and never in the published file, which has no way to declare them:

$$f^2 \equiv 1 + 2E, \qquad \mathcal{K} \equiv (\partial_t R)^2 - 2E, \qquad \mathcal{J} \equiv \partial_t R\,\partial_r\partial_t R - E'.$$

$\mathcal{J}$ is one half of $\partial_r\mathcal{K}$, which is the one relation between them, and everywhere a formula below carries any of the three the published file carries it written out.

---

## Step 2. The line element, and what the Lemaître-Tolman-Bondi form is

$$ds^2 = -c^2dt^2 + \frac{(\partial_r R)^2}{1+2E}dr^2 + R(r,t)^2\left(d\theta^2 + \sin^2\theta\,d\phi^2\right),$$

with $R$ positive and twice differentiable in both arguments, $E$ a function of $r$ alone with $1 + 2E > 0$, and $\theta$ and $\phi$ the usual angles on the sphere.

Three separate decisions are in that line element and it is worth taking them apart, because only the third is a choice of gauge.

The first is that the coordinates are comoving.
There is no $dt\,dr$ term and no shift of any kind, so the curves $r,\theta,\phi = \text{const}$ are worldlines, and Step 16 shows that they are geodesics.
The matter, when Step 11 puts matter in, is at rest in these coordinates for all time, which is what a dust can do and a fluid with pressure cannot.

The second is that the chart is synchronous.
The coefficient of $-c^2dt^2$ is exactly one, so $t$ is proper time along every one of those worldlines at once, and the surfaces $t = \text{const}$ are the slices orthogonal to them.
That is what makes $t$ a cosmic time in the sense the history of the entry means: every observer riding with the dust can set their clock by it.

The third is the radial coefficient, and it is the one that carries the name.
The general spherically symmetric comoving synchronous metric is

$$ds^2 = -c^2dt^2 + X(r,t)^2dr^2 + R(r,t)^2d\Omega^2$$

for two free functions $X$ and $R$, and Step 10 computes the mixed Einstein component of that general form,

$$G_{tr} = -\frac{2X}{R}\,\partial_t\!\left(\frac{\partial_r R}{X}\right).$$

So $G_{tr}$ vanishes exactly when $\partial_r R/X$ does not depend on $t$, that is when $X = \partial_r R/f(r)$ for some function $f$ of $r$ alone, and writing $f^2 = 1 + 2E$ gives the line element above.
The Lemaître-Tolman-Bondi form is therefore not an ansatz picked for its convenience.
It is the general solution of one of the field equations, the one that says there is no radial energy flux, and a spherically symmetric comoving dust has to satisfy it.
What is left free after that is $E(r)$, and what is left to impose is the rest of the field equations, which Step 11 does.

$R$ is the areal radius: the sphere labelled by $r$ at time $t$ has area $4\pi R^2$, by the angular part of the line element.
It is a length, and it is a function of both coordinates, since each shell of dust moves.
$r$ is a label for a shell and nothing more.
It is a length here only because the chart has to give $dr$ some dimension and $\partial_r R$ is then dimensionless, which is the reading Step 17 and the checker's `DIMENSIONS` table use; relabelling the shells by any increasing $\tilde{r}(r)$ leaves the line element in the same form with $\tilde{E}(\tilde{r}) = E(r)$, so nothing physical depends on which labelling is chosen.

$E$ is dimensionless and has two readings.
Locally it is the energy per unit mass of the shell $r$ in units of $c^2$, as Step 11 makes precise: the evolution equation it satisfies is the Newtonian energy equation of a particle falling radially, with $E$ the constant total energy.
Geometrically $1 + 2E$ carries the curvature of the spatial slice at that shell, exactly as $1 - kr^2$ does in FRW, which is why Step 13 finds $E = -kr^2/2$ in the homogeneous limit.
A shell with $E > 0$ has enough energy to expand forever, a shell with $E < 0$ must turn around and recollapse, and $E = 0$ is the marginally bound case that this entry deliberately does not specialise to.

The chart is regular where $R > 0$, $\partial_r R > 0$ and $1 + 2E > 0$.
Each of the three failures is a different thing.
$R \to 0$ with matter present is a curvature singularity, by Step 15, and it is the big bang of this model or its final crunch.
$1 + 2E \to 0$ is a place where the spatial slice closes off, the analogue of $r = 1/\sqrt{k}$ in a closed FRW universe, and the geometry continues through it in another chart.
$\partial_r R \to 0$ is a shell crossing, where two labels have reached the same areal radius, and Step 12 shows the density blowing up there while Step 15 shows the Kretschmann scalar blowing up with it.

---

## Step 3. The metric matrix, its determinant and its inverse

In the chart $x^\mu = (ct, r, \theta, \phi)$ the line element gives

$$g_{\mu\nu} = \begin{pmatrix} -1 & 0 & 0 & 0 \\ 0 & \dfrac{(\partial_r R)^2}{1+2E} & 0 & 0 \\ 0 & 0 & R^2 & 0 \\ 0 & 0 & 0 & R^2\sin^2\theta \end{pmatrix}.$$

Every entry is dimensionless, which is what a chart component of a metric has to be once all four coordinates are lengths.
The matrix is diagonal, so the inverse is the reciprocal entry by entry,

$$g^{tt} = -1, \qquad g^{rr} = \frac{1+2E}{(\partial_r R)^2}, \qquad g^{\theta\theta} = \frac{1}{R^2}, \qquad g^{\phi\phi} = \frac{1}{R^2\sin^2\theta},$$

and the determinant is

$$\det g = -\frac{(\partial_r R)^2R^4\sin^2\theta}{1+2E}, \qquad \sqrt{-g} = \frac{\partial_r R\,R^2\sin\theta}{\sqrt{1+2E}}.$$

The determinant is nonzero wherever the three regularity conditions of Step 2 hold, and the volume element $\sqrt{-g}\,dr\,d\theta\,d\phi$ is the one Step 12 integrates to get the mass inside a shell.

---

## Step 4. Every Christoffel symbol

The metric is diagonal, so the sum over $\alpha$ collapses to the one term $\alpha = \mu$, and

$$\Gamma^\mu{}_{\nu\rho} = \frac{1}{2g_{\mu\mu}}\left(\partial_\nu g_{\mu\rho} + \partial_\rho g_{\mu\nu} - \partial_\mu g_{\nu\rho}\right)$$

with no sum on $\mu$.
Four of them are worth writing out, and the rest go the same way.

The symbol with an upper $t$ and two lower $r$ needs only the time derivative of $g_{rr}$, since $g_{tt}$ is a constant and $g_{tr}$ is zero:

$$\Gamma^t{}_{rr} = -\frac{1}{2g_{tt}}\partial_t g_{rr} = \frac{1}{2}\partial_t\frac{(\partial_r R)^2}{1+2E} = \frac{\partial_r R\,\partial_r\partial_t R}{1+2E},$$

where $E$ does not depend on $t$ and so passes through the derivative.
Its angular counterparts are the same computation on $g_{\theta\theta}$ and $g_{\phi\phi}$,

$$\Gamma^t{}_{\theta\theta} = \frac{1}{2}\partial_t R^2 = R\,\partial_t R, \qquad \Gamma^t{}_{\phi\phi} = R\,\partial_t R\sin^2\theta.$$

The symbol with an upper $r$ and a mixed pair is

$$\Gamma^r{}_{tr} = \frac{1}{2g_{rr}}\partial_t g_{rr} = \frac{1+2E}{2(\partial_r R)^2}\cdot\frac{2\,\partial_r R\,\partial_r\partial_t R}{1+2E} = \frac{\partial_r\partial_t R}{\partial_r R},$$

which is the fractional rate at which the proper radial distance between two neighbouring shells is stretching.
The purely radial one carries the only appearance of $E'$ in the connection:

$$\Gamma^r{}_{rr} = \frac{1}{2g_{rr}}\partial_r g_{rr} = \frac{1+2E}{2(\partial_r R)^2}\left(\frac{2\,\partial_r R\,\partial_r^2R}{1+2E} - \frac{2E'(\partial_r R)^2}{(1+2E)^2}\right) = \frac{(1+2E)\partial_r^2R - E'\,\partial_r R}{(1+2E)\partial_r R}.$$

The remaining ones are

$$\Gamma^r{}_{\theta\theta} = -\frac{1}{2g_{rr}}\partial_r g_{\theta\theta} = -\frac{(1+2E)R}{\partial_r R}, \qquad \Gamma^r{}_{\phi\phi} = -\frac{(1+2E)R\sin^2\theta}{\partial_r R},$$

$$\Gamma^\theta{}_{t\theta} = \Gamma^\phi{}_{t\phi} = \frac{\partial_t R}{R}, \qquad \Gamma^\theta{}_{r\theta} = \Gamma^\phi{}_{r\phi} = \frac{\partial_r R}{R},$$

$$\Gamma^\theta{}_{\phi\phi} = -\sin\theta\cos\theta, \qquad \Gamma^\phi{}_{\theta\phi} = \cot\theta.$$

The last two are the connection of the unit two sphere and know nothing about the rest of the geometry.
The pair $\partial_t R/R$ and $\partial_r R/R$ is the transverse expansion rate and the transverse gradient, and their equality across $\theta$ and $\phi$ is the spherical symmetry.

Counting the symmetric partners $\Gamma^\mu{}_{\nu\rho} = \Gamma^\mu{}_{\rho\nu}$, the entry publishes nineteen symbols with an upper index, from eleven distinct values.

---

## Step 5. The lowered Christoffel symbols

$\Gamma_{\mu\nu\rho} = g_{\mu\alpha}\Gamma^\alpha{}_{\nu\rho}$ is one multiplication per symbol, since the metric is diagonal.
The time row picks up the minus sign of $g_{tt} = -1$,

$$\Gamma_{trr} = -\frac{\partial_r R\,\partial_r\partial_t R}{1+2E}, \qquad \Gamma_{t\theta\theta} = -R\,\partial_t R, \qquad \Gamma_{t\phi\phi} = -R\,\partial_t R\sin^2\theta,$$

while the radial row multiplies by $g_{rr}$,

$$\Gamma_{rtr} = \frac{\partial_r R\,\partial_r\partial_t R}{1+2E}, \qquad \Gamma_{rrr} = \frac{\partial_r R\left((1+2E)\partial_r^2R - E'\,\partial_r R\right)}{(1+2E)^2},$$

$$\Gamma_{r\theta\theta} = -R\,\partial_r R, \qquad \Gamma_{r\phi\phi} = -R\,\partial_r R\sin^2\theta,$$

and the angular rows multiply by $R^2$ and $R^2\sin^2\theta$,

$$\Gamma_{\theta t\theta} = R\,\partial_t R, \qquad \Gamma_{\theta r\theta} = R\,\partial_r R, \qquad \Gamma_{\theta\phi\phi} = -R^2\sin\theta\cos\theta,$$

$$\Gamma_{\phi t\phi} = R\,\partial_t R\sin^2\theta, \qquad \Gamma_{\phi r\phi} = R\,\partial_r R\sin^2\theta, \qquad \Gamma_{\phi\theta\phi} = R^2\sin\theta\cos\theta.$$

That is nineteen more published symbols, from seven distinct values, since lowering the index has merged several that differed only by a factor of the metric.

---

## Step 6. The Riemann tensor in the plane of time and radius

This is the one component family where the third derivative of $R$ appears, and it is worth doing in full because everything later either inherits it or cancels it.

$$R^t{}_{rtr} = \partial_t\Gamma^t{}_{rr} - \partial_r\Gamma^t{}_{rt} + \Gamma^t{}_{t\lambda}\Gamma^\lambda{}_{rr} - \Gamma^t{}_{r\lambda}\Gamma^\lambda{}_{rt}.$$

The second and third terms vanish, because $\Gamma^t{}_{rt}$ and $\Gamma^t{}_{t\lambda}$ are zero for every $\lambda$, and the fourth has only $\lambda = t$, so

$$R^t{}_{rtr} = \partial_t\left(\frac{\partial_r R\,\partial_r\partial_t R}{1+2E}\right) - \frac{\partial_r R\,\partial_r\partial_t R}{1+2E}\cdot\frac{\partial_r\partial_t R}{\partial_r R} = \frac{(\partial_r\partial_t R)^2 + \partial_r R\,\partial_r\partial_t^2R}{1+2E} - \frac{(\partial_r\partial_t R)^2}{1+2E},$$

$$R^t{}_{rtr} = \frac{\partial_r R\,\partial_r\partial_t^2R}{1+2E}.$$

The square of the first derivative cancels and a third derivative of $R$ is left standing.
There is nothing to be done about it.
The radial metric component already carries one derivative of $R$, and a curvature carries two derivatives of the metric, so the curvature of this line element reaches the third derivative of the function the line element is written with, which is one order past every other entry in the collection.
Step 18 says what that cost at the checker, until the checker learned to read it.

The partner with the indices the other way up is the same expression divided by $g_{rr}$ and multiplied by $g_{tt}$,

$$R^r{}_{ttr} = \frac{\partial_r\partial_t^2R}{\partial_r R},$$

and the fully lowered one is $R_{trtr} = g_{tt}R^t{}_{rtr} = -\partial_r R\,\partial_r\partial_t^2R/(1+2E)$.

---

## Step 7. The Riemann tensor in the remaining planes

The plane of time and an angle is the tidal acceleration of the transverse directions:

$$R^t{}_{\theta t\theta} = \partial_t\Gamma^t{}_{\theta\theta} - \Gamma^t{}_{\theta\theta}\Gamma^\theta{}_{\theta t} = \left((\partial_t R)^2 + R\,\partial_t^2R\right) - R\,\partial_t R\cdot\frac{\partial_t R}{R} = R\,\partial_t^2R,$$

and the same with $\phi$ in place of $\theta$ carries the extra $\sin^2\theta$ of $g_{\phi\phi}$,

$$R^t{}_{\phi t\phi} = R\,\partial_t^2R\sin^2\theta.$$

The mixed versions follow by raising and lowering: $R^\theta{}_{tt\theta} = \partial_t^2R/R$ and $R_{t\theta t\theta} = -R\,\partial_t^2R$.
Only the second time derivative of $R$ is here, which is what a tidal acceleration is.

The plane of radius and an angle brings in $E'$ for the first time:

$$R^r{}_{\theta r\theta} = \partial_r\Gamma^r{}_{\theta\theta} + \Gamma^r{}_{rt}\Gamma^t{}_{\theta\theta} + \Gamma^r{}_{rr}\Gamma^r{}_{\theta\theta} - \Gamma^r{}_{\theta\theta}\Gamma^\theta{}_{\theta r}.$$

Writing the four terms out,

$$\partial_r\Gamma^r{}_{\theta\theta} = -\frac{2E'R + (1+2E)\partial_r R}{\partial_r R} + \frac{(1+2E)R\,\partial_r^2R}{(\partial_r R)^2},$$

$$\Gamma^r{}_{rt}\Gamma^t{}_{\theta\theta} = \frac{R\,\partial_t R\,\partial_r\partial_t R}{\partial_r R}, \qquad \Gamma^r{}_{rr}\Gamma^r{}_{\theta\theta} = -\frac{(1+2E)R\,\partial_r^2R}{(\partial_r R)^2} + \frac{E'R}{\partial_r R},$$

$$-\Gamma^r{}_{\theta\theta}\Gamma^\theta{}_{\theta r} = (1+2E).$$

The two terms in $\partial_r^2R$ cancel against each other and the two in $(1+2E)$ cancel against each other, leaving

$$R^r{}_{\theta r\theta} = \frac{R\left(\partial_t R\,\partial_r\partial_t R - E'\right)}{\partial_r R} = \frac{R\,\mathcal{J}}{\partial_r R},$$

with $R^r{}_{\phi r\phi}$ the same times $\sin^2\theta$ and $R^\theta{}_{rr\theta} = -\partial_r R\,\mathcal{J}/\left((1+2E)R\right)$.

The plane of the two angles is the intrinsic curvature of the sphere, corrected by the motion:

$$R^\theta{}_{\phi\theta\phi} = \partial_\theta\Gamma^\theta{}_{\phi\phi} + \Gamma^\theta{}_{\theta t}\Gamma^t{}_{\phi\phi} + \Gamma^\theta{}_{\theta r}\Gamma^r{}_{\phi\phi} - \Gamma^\theta{}_{\phi\phi}\Gamma^\phi{}_{\phi\theta},$$

$$R^\theta{}_{\phi\theta\phi} = -\cos 2\theta + (\partial_t R)^2\sin^2\theta - (1+2E)\sin^2\theta + \cos^2\theta = \left((\partial_t R)^2 - 2E\right)\sin^2\theta = \mathcal{K}\sin^2\theta.$$

That combination is the one the whole solution turns on.
It is the sectional curvature of the sphere, which would be $1/R^2$ for a sphere sitting still in flat space, reduced by the expansion and by the spatial curvature of the slice, and Step 11 will find that the field equations make it a mass divided by a radius.

Every component not in one of these four families vanishes.
The entry publishes twenty four components of $R^\mu{}_{\nu\rho\sigma}$ from ten distinct values and twenty four of $R_{\mu\nu\rho\sigma}$ from six, the fully lowered block being shorter because lowering the first index merges pairs that the mixed block keeps apart.

---

## Step 8. The Ricci tensor

Contracting on the first lower index, $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$, each diagonal entry is a sum of three of the components above.

$$R_{tt} = R^r{}_{trt} + R^\theta{}_{t\theta t} + R^\phi{}_{t\phi t} = -\frac{\partial_r\partial_t^2R}{\partial_r R} - \frac{2\,\partial_t^2R}{R} = -\frac{R\,\partial_r\partial_t^2R + 2\,\partial_r R\,\partial_t^2R}{R\,\partial_r R}.$$

This is the combination Raychaudhuri's equation is about.
The comoving congruence of Step 16 is geodesic and has no rotation, so its expansion is $\Theta = \partial_t\ln\sqrt{-g} = \partial_r\partial_t R/\partial_r R + 2\,\partial_t R/R$, the sum of one radial and two transverse rates, and the identity above rearranges to

$$R_{tt} = -\partial_t\Theta - \left(\frac{\partial_r\partial_t R}{\partial_r R}\right)^2 - 2\left(\frac{\partial_t R}{R}\right)^2,$$

which is Raychaudhuri's equation with the quadratic terms written out as the squares of those three rates.
Step 12 makes the left hand side $4\pi G\rho/c^2$, so a positive density drives $\partial_t\Theta$ down and the expansion decelerates.

$$R_{rr} = R^t{}_{rtr} + R^\theta{}_{r\theta r} + R^\phi{}_{r\phi r} = \frac{\partial_r R\left(R\,\partial_r\partial_t^2R + 2\mathcal{J}\right)}{(1+2E)R},$$

$$R_{\theta\theta} = R^t{}_{\theta t\theta} + R^r{}_{\theta r\theta} + R^\phi{}_{\theta\phi\theta} = \frac{\partial_r R\left(R\,\partial_t^2R + \mathcal{K}\right) + R\,\mathcal{J}}{\partial_r R},$$

with $R_{\phi\phi} = R_{\theta\theta}\sin^2\theta$.
The mixed components all vanish, including $R_{tr}$, which Step 10 makes the most of.

The mixed and doubly raised variants are one or two factors of the inverse metric away, and the entry publishes all three, four components each.

---

## Step 9. The Ricci scalar and the Einstein tensor

$$R = g^{\mu\nu}R_{\mu\nu} = -R_{tt} + \frac{(1+2E)R_{rr}}{(\partial_r R)^2} + \frac{2R_{\theta\theta}}{R^2},$$

$$R = \frac{2\left(R^2\,\partial_r\partial_t^2R + 2R\,\partial_r R\,\partial_t^2R + 2R\,\mathcal{J} + \partial_r R\,\mathcal{K}\right)}{R^2\,\partial_r R}.$$

The Einstein tensor is $G_{\mu\nu} = R_{\mu\nu} - \tfrac{1}{2}Rg_{\mu\nu}$, and its time component is where the third derivative goes away.
With $g_{tt} = -1$,

$$G_{tt} = R_{tt} + \frac{R}{2} = \left(-\frac{\partial_r\partial_t^2R}{\partial_r R} - \frac{2\,\partial_t^2R}{R}\right) + \left(\frac{\partial_r\partial_t^2R}{\partial_r R} + \frac{2\,\partial_t^2R}{R} + \frac{2\mathcal{J}}{R\,\partial_r R} + \frac{\mathcal{K}}{R^2}\right),$$

$$G_{tt} = \frac{\partial_r R\,\mathcal{K} + 2R\,\mathcal{J}}{R^2\,\partial_r R}.$$

Everything with two time derivatives in it has cancelled, and so has the third derivative with it.
That is not luck.
$G_{tt}$ is the Hamiltonian constraint of the slicing, and a constraint by definition involves only the slice and its first time derivative; here the slice carries $R$ and $\partial_r R$ and its first time derivative carries $\partial_t R$ and $\partial_r\partial_t R$, which is exactly the list above.

The same cancellation happens in the radial component, which is the other constraint direction in disguise, and

$$G^r{}_r = -\frac{2R\,\partial_t^2R + \mathcal{K}}{R^2}.$$

It does not happen in the angular components, which are genuine evolution equations:

$$G^\theta{}_\theta = G^\phi{}_\phi = -\frac{R\,\partial_r\partial_t^2R + \partial_r R\,\partial_t^2R + \mathcal{J}}{R\,\partial_r R}.$$

The entry publishes all three variants of $G_{\mu\nu}$, four components each.

---

## Step 10. What the Einstein tensor says: no flux, a density and two pressures

Read through the field equations $G_{\mu\nu} = 8\pi G\,T_{\mu\nu}/c^4$, the four components above are a statement about the matter, and it is worth being careful about which four numbers they are.

The observer to read them with is the comoving one, with four velocity $u^\mu = (1,0,0,0)$ in this chart and $u_\mu = (-1,0,0,0)$, which is a unit timelike vector because $g_{tt} = -1$.
For that observer the energy density is $T_{\mu\nu}u^\mu u^\nu = T_{tt}$, the energy flux is $T_{tr}$, and the principal pressures are $c^4G^r{}_r/8\pi G$ radially and $c^4G^\theta{}_\theta/8\pi G$ in the two transverse directions, which are equal by spherical symmetry.

The first thing the Einstein tensor says is that the flux is zero.
$G_{tr}$ vanishes identically, for every $R$ and every $E$, so

$$T_{tr} = 0$$

with nothing imposed.
That is the content of the third decision of Step 2, and it is worth seeing where it comes from rather than reading it off a table of components.
For the general comoving synchronous metric $ds^2 = -c^2dt^2 + X^2dr^2 + R^2d\Omega^2$ the same computation as Step 6 and Step 7, with $X$ in place of $\partial_r R/\sqrt{1+2E}$, gives

$$G_{tr} = R_{tr} = -\frac{2X}{R}\,\partial_t\!\left(\frac{\partial_r R}{X}\right),$$

which vanishes if and only if $\partial_r R/X$ is a function of $r$ alone.
Calling that function $\sqrt{1+2E}$ is the line element of this entry.
So the entry's radial coefficient is exactly the statement that no energy crosses a comoving shell, and the matter it can hold is therefore matter that stays where it is put.

The remaining three components give a density and two pressures:

$$\rho = \frac{c^2}{8\pi G}G_{tt} = \frac{c^2}{8\pi G}\cdot\frac{\partial_r R\,\mathcal{K} + 2R\,\mathcal{J}}{R^2\,\partial_r R},$$

$$p_r = \frac{c^4}{8\pi G}G^r{}_r = -\frac{c^4}{8\pi G}\cdot\frac{2R\,\partial_t^2R + \mathcal{K}}{R^2}, \qquad p_\perp = \frac{c^4}{8\pi G}G^\theta{}_\theta = -\frac{c^4}{8\pi G}\cdot\frac{R\,\partial_r\partial_t^2R + \partial_r R\,\partial_t^2R + \mathcal{J}}{R\,\partial_r R}.$$

The density here is a mass density, so that $T_{tt} = \rho c^2$ and $G_{tt} = 8\pi G\rho/c^2$, which is the normalisation the FRW entry uses as well.
Nothing so far is dust.
An arbitrary $R$ and an arbitrary $E$ give an anisotropic fluid at rest, with a radial pressure and a transverse pressure that need not agree with each other or vanish.

One relation among those three holds whatever $R$ and $E$ are, and it is the one that does the work in the next step.
Differentiating $R^2G^r{}_r$ with respect to $r$ and comparing against $G^\theta{}_\theta$ gives the identity

$$G^\theta{}_\theta = \frac{\partial_r\left(R^2G^r{}_r\right)}{2R\,\partial_r R}, \qquad \text{equivalently} \qquad \partial_r p_r = \frac{2\left(p_\perp - p_r\right)\partial_r R}{R}.$$

The second form is the radial equation of hydrostatic balance for an anisotropic fluid, with no gravitational term in it because the flow is geodesic, and it is nothing more than the Bianchi identity $\nabla_\mu G^{\mu}{}_{r} = 0$ written out.
It says that the transverse pressure is not independent: it is fixed by the radial pressure and its gradient.

---

## Step 11. Dust, which is one equation and not three

A dust is matter with no pressure at all, so the condition is $p_r = p_\perp = 0$.
By the identity at the end of Step 10, the second of those follows from the first.
If $p_r$ vanishes at every $r$, then $\partial_r p_r$ vanishes too, and the identity forces $p_\perp = p_r = 0$ wherever $\partial_r R \neq 0$.
That is the sense in which the vanishing of the transverse pressure is not an assumption imposed on this solution: it is a consequence of the radial one vanishing, and no second condition is available to impose.

So the whole of dust is $G^r{}_r = 0$, which by Step 9 is

$$2R\,\partial_t^2R + (\partial_t R)^2 - 2E = 0.$$

This integrates once in $t$ at fixed $r$.
Multiply by $\partial_t R$ and notice that the left hand side is then an exact derivative,

$$\partial_t\left[R\left((\partial_t R)^2 - 2E\right)\right] = \partial_t R\left(2R\,\partial_t^2R + (\partial_t R)^2 - 2E\right),$$

so the bracket is constant in time.
Its value is a function of $r$ alone, and writing that function as $2GM(r)/c^2$, which is a length because $M$ is a mass,

$$R\,\mathcal{K} = R\left((\partial_t R)^2 - 2E\right) = \frac{2GM(r)}{c^2}, \qquad \text{that is} \qquad (\partial_t R)^2 = 2E + \frac{2GM(r)}{c^2R}.$$

That is the Lemaître-Tolman-Bondi evolution equation, and the arbitrary function it introduces is the mass function.
Multiplying by $c^2/2$ turns it into $\tfrac{1}{2}\dot{R}^2 = c^2E + GM/R$, the energy equation of a Newtonian particle falling radially in the field of a mass $M$, with $c^2E$ its total energy per unit mass.
That is the reading of $E$ promised in Step 2, and it is why $E < 0$ recollapses and $E \geq 0$ does not.
Differentiating the evolution equation once more gives the acceleration in the same Newtonian form,

$$\partial_t^2R = -\frac{GM(r)}{c^2R^2},$$

and each shell falls in the field of the mass interior to it and is deaf to everything outside, which is Birkhoff's theorem doing its work shell by shell.
That is the independence the entry's history describes: the model is a stack of noninteracting Newtonian shells, and its exact solution in general relativity is the same one Newton would have written.

The three cases of the evolution equation are the three signs of $E$, and the general solution is parametric rather than explicit.
For $E > 0$ it is $R = \left(GM/2c^2E\right)\left(\cosh\eta - 1\right)$ with $c(t - t_B(r)) = \left(GM/c^2\right)\left(2E\right)^{-3/2}\left(\sinh\eta - \eta\right)$, for $E < 0$ the same with hyperbolic functions replaced by circular ones and $2E$ by $-2E$, and only the marginally bound case $E = 0$ closes in elementary form as $R = \left(9GM(r)(t - t_B)^2/2\right)^{1/3}$.
The third arbitrary function $t_B(r)$ is the bang time, the moment at which each shell had zero radius.
None of the three is assumed anywhere in this file, and the entry publishes the general case for exactly that reason: writing the parabolic case would have bought a closed form for $R$ at the price of a solution that no longer has a bound shell in it.

---

## Step 12. The density through the mass function

With the evolution equation in hand the time component of Step 9 collapses.
Substituting $\mathcal{K} = 2GM/(c^2R)$ and its radial derivative $\partial_r\mathcal{K} = 2\mathcal{J} = 2GM'/(c^2R) - 2GM\,\partial_r R/(c^2R^2)$ into

$$G_{tt} = \frac{\partial_r R\,\mathcal{K} + 2R\,\mathcal{J}}{R^2\,\partial_r R}$$

gives, after the two terms in $GM\,\partial_r R/R$ cancel,

$$G_{tt} = \frac{2GM'(r)}{c^2R^2\,\partial_r R} = \frac{8\pi G}{c^2}\rho, \qquad \rho(r,t) = \frac{M'(r)}{4\pi R^2\,\partial_r R}.$$

That is the density of the dust, and it is the published $G_{tt}$ read through the field equations and nothing else.
It is worth seeing why it has to take that shape.
The proper volume of the shell between $r$ and $r + dr$ is $4\pi R^2\,\partial_r R\,dr/\sqrt{1+2E}$, by the volume element of Step 3, so the rest mass it contains is $\rho$ times that, namely $M'dr/\sqrt{1+2E}$.
$M(r)$ is therefore not quite the rest mass inside the shell: it is the active gravitating mass, the rest mass corrected by the binding energy through the factor $\sqrt{1+2E}$, which is the same correction that separates the two in a static star.
What matters for the geometry is that $M$ does not depend on $t$.
The dust is comoving, no shell crosses another while $\partial_r R > 0$, and so the mass inside a given label is carried along unchanged; that conservation is precisely the equation $G^r{}_r = 0$ of Step 11, since

$$G^r{}_r = -\frac{2}{R^2\,\partial_t R}\,\partial_t\!\left[\frac{R\,\mathcal{K}}{2}\right]$$

and the bracket is $GM/c^2$.
The combination $\tfrac{1}{2}R\mathcal{K}$ in it is the Misner-Sharp mass of the sphere, defined off the shell as well as on it, and the two Einstein components of Step 9 that carry no third derivative are exactly its radial and time derivatives:

$$G_{tt} = \frac{2\,\partial_r\mu}{R^2\,\partial_r R}, \qquad G^r{}_r = -\frac{2\,\partial_t\mu}{R^2\,\partial_t R}, \qquad \mu \equiv \frac{R\left((\partial_t R)^2 - 2E\right)}{2}.$$

Vanishing radial pressure says the Misner-Sharp mass is conserved along the flow, and the density is its radial gradient per unit areal volume.

Two more readings come free.
The Ricci scalar of Step 9, with the same substitution, becomes

$$R = \frac{2GM'}{c^2R^2\,\partial_r R} = \frac{8\pi G\rho}{c^2},$$

which is the trace of the field equations for a pressureless source, and the Ricci tensor becomes $R_{\mu\nu} = \left(4\pi G\rho/c^2\right)\left(g_{\mu\nu} + 2u_\mu u_\nu\right)$, the dust form.
The density diverges in two different ways, and they are different singularities.
Where $R \to 0$ with $M' \neq 0$ the whole shell has collapsed to a point and the geometry is singular with it, by Step 15.
Where $\partial_r R \to 0$ with $M' \neq 0$ two neighbouring shells have arrived at the same areal radius, the proper distance between them has gone to zero, and a finite mass sits in no volume; that is a shell crossing, and it is a real curvature singularity in this model although a mild one, in the sense that the solution can often be continued through it once the dust is allowed to be something a little more realistic.

---

## Step 13. The homogeneous limit, which is FRW

The collection publishes the FRW universe as an entry of its own, and a reader who meets an inhomogeneous dust will want to know exactly how the two are related.
Set

$$R(r,t) = a(t)\,r, \qquad E(r) = -\frac{kr^2}{2}.$$

The line element becomes

$$ds^2 = -c^2dt^2 + \frac{a^2dr^2}{1-kr^2} + a^2r^2d\Omega^2,$$

which is the FRW comoving spherical line element exactly as that entry prints it.
Every published component of this entry reduces to the one the FRW entry prints beside it, which is forced rather than lucky: the substitution turns this line element into that one, and each entry's components are the curvature of its own line element and nothing else.
Writing $a' = \partial_t a$ for the chart derivative, as the Bianchi entry does,

$$G_{tt} = \frac{3\left(a'^2 + k\right)}{a^2}, \qquad G^r{}_r = G^\theta{}_\theta = -\frac{2a\,a'' + a'^2 + k}{a^2}, \qquad R = \frac{6\left(a\,a'' + a'^2 + k\right)}{a^2},$$

and the first of those is the Friedmann equation with the same positive sign that the FRW entry publishes, which is a check on the contraction convention of Step 1 as much as on the algebra.
The Kretschmann scalar of Step 15 becomes

$$K = \frac{12a''^2}{a^2} + \frac{12\left(a'^2+k\right)^2}{a^4},$$

which is the FRW entry's Kretschmann scalar, and the Weyl tensor of Step 14 becomes zero in every one of its twenty four slots, which is the FRW entry's empty Weyl block and the statement that a homogeneous isotropic universe is conformally flat.

The dust condition in the limit is the Friedmann equation, and the mass function is

$$\frac{2GM(r)}{c^2} = R\left((\partial_t R)^2 - 2E\right) = a r^3\left(a'^2 + k\right), \qquad \rho = \frac{M'}{4\pi R^2\partial_r R} = \frac{3c^2\left(a'^2+k\right)}{8\pi G a^2},$$

a density that does not depend on $r$, which is homogeneity, and a mass that grows as the cube of the label, which is a uniform ball.
The inhomogeneity of the general case is exactly the failure of $M$ to be proportional to $r^3$ at fixed $a$, or of $E$ to be proportional to $r^2$, and Step 14 shows that the Weyl tensor measures the first of those two directly.

---

## Step 14. The Weyl tensor, computed rather than copied

This spacetime is not a vacuum, so the Weyl tensor is not the Riemann tensor and it is computed here from the definition of Step 1 by removing the traces.
Doing that leaves every nonzero component proportional to a single scalar,

$$\mathcal{W} \equiv R\left(\partial_r R\,\partial_t^2R + \partial_t R\,\partial_r\partial_t R - E' - R\,\partial_r\partial_t^2R\right) - \partial_r R\left((\partial_t R)^2 - 2E\right),$$

with

$$C_{trtr} = \frac{\partial_r R\,\mathcal{W}}{3R^2(1+2E)}, \qquad C_{t\theta t\theta} = -\frac{\mathcal{W}}{6\,\partial_r R}, \qquad C_{r\theta r\theta} = \frac{\partial_r R\,\mathcal{W}}{6(1+2E)}, \qquad C_{\theta\phi\theta\phi} = -\frac{R^2\mathcal{W}\sin^2\theta}{3\,\partial_r R},$$

and the rest from these by the symmetries.
In the orthonormal frame carried by the comoving observer, whose radial leg is $\partial_r$ divided by $\sqrt{g_{rr}}$ and whose transverse legs are $\partial_\theta$ and $\partial_\phi$ divided by $R$ and $R\sin\theta$, those four numbers become

$$C_{\hat t\hat r\hat t\hat r} = \frac{\mathcal{W}}{3R^2\,\partial_r R} = -2\,C_{\hat t\hat\theta\hat t\hat\theta},$$

so the radial tide is minus twice the transverse one.
That is the tidal pattern of a point mass, stretching along the radius and squeezing across it, and a single scalar carrying all of it is what makes the spacetime Petrov type D, which every spherically symmetric spacetime is.
The twenty four components of each variant print on six and eight lines respectively, because so many of them are the same value with a sign.

That the Weyl tensor is not the Riemann tensor is visible in any one slot.
Step 7 has $R_{t\theta t\theta} = -R\,\partial_t^2R$ while this step has $C_{t\theta t\theta} = -\mathcal{W}/(6\,\partial_r R)$, and the two agree only where $\mathcal{W} = 6R\,\partial_r R\,\partial_t^2R$, which on the dust shell of Step 11 is the condition $M' = 0$, that is a region with no matter in it.
That is the general rule recorded in `weyl.md` arriving in a particular case: the two tensors coincide in vacuum because the traces removed are the ones that vanish there, and nowhere else.
Here they differ in all twenty four slots wherever the density is nonzero, so the entry's Weyl block is what the definition gives rather than a copy of anything.

On the dust shell of Step 11 the scalar becomes something worth reading.
Substituting the evolution equation and writing $m \equiv GM/c^2$ for the mass in length units,

$$\mathcal{W} = 2m' - \frac{6m\,\partial_r R}{R}, \qquad C^t{}_{\theta t\theta} = \frac{m'}{3\,\partial_r R} - \frac{m}{R} = \frac{4\pi G\rho R^2}{3c^2} - \frac{GM}{c^2R}.$$

The second term is the mass actually inside the shell, and the first is the mass a uniform ball at the local density would have inside that same shell.
The Weyl tensor is the difference between them.
It vanishes when and only when $M = 4\pi\rho R^3/3$, which is the homogeneous case of Step 13, and it is what a tidal force measures: an observer inside a uniform ball feels no tide, and an observer inside a lump feels the difference between the lump and the uniform ball that would fit around them.
That is also why the Oppenheimer-Snyder interior, which is a uniform ball of this dust, is conformally flat and joins onto a Schwarzschild exterior whose Weyl tensor is all there is.

---

## Step 15. The Kretschmann scalar

$K = R_{\mu\nu\rho\sigma}R^{\mu\nu\rho\sigma}$ is a sum over the four planes of Step 6 and Step 7, each contributing the square of its sectional curvature with multiplicity, and for this metric it collects into four squares:

$$K = \frac{4\left[(\partial_r R)^2\mathcal{K}^2 + 2R^2\mathcal{J}^2 + R^4\left(\partial_r\partial_t^2R\right)^2 + 2R^2(\partial_r R)^2\left(\partial_t^2R\right)^2\right]}{R^4(\partial_r R)^2}.$$

Every term is a square with a positive coefficient, so $K \geq 0$ for this whole family, and it vanishes only where all four vanish together, which is flat space.

On the dust shell the third and fourth squares turn into the mass and its gradient, and the whole thing reduces to

$$K = \frac{48m^2}{R^6} - \frac{32m\,m'}{R^5\,\partial_r R} + \frac{12m'^2}{R^4(\partial_r R)^2}, \qquad m = \frac{GM}{c^2}.$$

Three readings follow from those three terms.
Set $m' = 0$, which is a shell with vacuum around it, and $K = 48m^2/R^6 = 12r_s^2/R^6$ is the Schwarzschild Kretschmann scalar, as Birkhoff's theorem requires.
Set $m = 4\pi\rho R^3/3$ with $\rho$ homogeneous, which is the FRW case, and the three terms collect into the FRW scalar quoted in Step 13.
Let $R \to 0$ at fixed $r$ with $m \neq 0$ and the first term diverges as $R^{-6}$: that is the big bang of this model, or the central singularity of its collapse, and it is a curvature singularity and not a coordinate one.
Let $\partial_r R \to 0$ with $m' \neq 0$ and the last term diverges: that is the shell crossing of Step 12, and the density of Step 12 diverges with it.

---

## Step 16. The geodesic equations

With the Christoffel symbols of Step 4 and the dots of Step 1, $\ddot{x}^\mu + \Gamma^\mu{}_{\nu\rho}\dot{x}^\nu\dot{x}^\rho = 0$ reads

$$\ddot{t} + \frac{\partial_r R\,\partial_r\partial_t R}{1+2E}\dot{r}^2 + R\,\partial_t R\left(\dot{\theta}^2 + \sin^2\theta\,\dot{\phi}^2\right) = 0,$$

$$\ddot{r} + \frac{2\,\partial_r\partial_t R}{\partial_r R}\dot{t}\dot{r} + \frac{(1+2E)\partial_r^2R - E'\,\partial_r R}{(1+2E)\partial_r R}\dot{r}^2 - \frac{(1+2E)R}{\partial_r R}\left(\dot{\theta}^2 + \sin^2\theta\,\dot{\phi}^2\right) = 0,$$

$$\ddot{\theta} + \frac{2\,\partial_t R}{R}\dot{t}\dot{\theta} + \frac{2\,\partial_r R}{R}\dot{r}\dot{\theta} - \sin\theta\cos\theta\,\dot{\phi}^2 = 0,$$

$$\ddot{\phi} + \frac{2\,\partial_t R}{R}\dot{t}\dot{\phi} + \frac{2\,\partial_r R}{R}\dot{r}\dot{\phi} + 2\cot\theta\,\dot{\theta}\dot{\phi} = 0.$$

The factors of two on the mixed terms are the two orderings of the same symmetric Christoffel symbol, which the sum counts twice and the entry lists twice.

The first thing to check is the comoving worldline, $r,\theta,\phi$ constant.
Every term of the last three equations vanishes on it, and the first reduces to $\ddot{t} = 0$, so it is a geodesic and its affine parameter is $ct$ up to a linear change.
The dust is therefore in free fall, which is the only way a pressureless fluid can move, and the coordinate $t$ is proper time along it, which is what synchronous means.
No pressure gradient exists to push a dust element off that path, and Step 11 is that statement in the language of the field equations.

One first integral survives and the one a reader will look for first does not.
The metric does not depend on $\phi$, so $L = R^2\sin^2\theta\,\dot{\phi}$ is conserved and the orbit can be put in the plane $\theta = \pi/2$, where $\dot{\theta} = 0$ solves the third equation identically.
What does not survive is a conserved energy.
The metric depends on $t$ through $R$, so $\partial_t$ is not a Killing vector and there is no analogue of the $E$ of a static spacetime, which is why light climbing out of this geometry is redshifted by the expansion rather than by a potential.
That absence is worth naming: it is the same absence that makes the FRW redshift a stretching of wavelengths, and it is the reason an inhomogeneous void model can imitate an accelerating universe in the way the entry's history describes.

---

## Step 17. Every published expression is dimensionally consistent

The declarations are $[t] = T$, $[r] = L$, $[\theta] = [\phi] = 1$, $[R] = L$ and $[E] = 1$.
Because the chart coordinate is $ct$, all four chart coordinates are lengths, which makes every chart component of the metric dimensionless and fixes the rest.
A chart derivative of $R$ divides a length by a length, so

$$[\partial_r R] = [\partial_t R] = 1, \qquad [\partial_r^2R] = [\partial_t^2R] = [\partial_r\partial_t R] = [E'] = \frac{1}{L}, \qquad [\partial_r\partial_t^2R] = \frac{1}{L^2},$$

each order of differentiation costing one inverse length.
The abbreviations of Step 1 then carry $[f^2] = 1$, $[\mathcal{K}] = 1$ and $[\mathcal{J}] = 1/L$.

The rule a published component is measured against is that a field with every index down carries its own dimension times $L$ for each lower index divided by that index's own dimension, and the inverse for an upper one; the field itself carries $1$ for the metric, $1/L$ for a Christoffel symbol, $1/L^2$ for Riemann, Ricci, Einstein, Weyl and the Ricci scalar, and $1/L^4$ for the Kretschmann scalar.

Four worked cases cover every pattern in the file.

$\Gamma^t{}_{rr} = \partial_r R\,\partial_r\partial_t R/(1+2E)$ carries $1 \cdot L^{-1} = L^{-1}$, and its index weight is $1$ because $t$ and $r$ are both lengths in the chart, so it matches the $1/L$ a Christoffel symbol wants.

$R^\theta{}_{\phi\theta\phi} = \mathcal{K}\sin^2\theta$ is dimensionless, and its index weight is $L^{-1}$ for the upper angle and $L$ for each of the three lower ones, so the target is $L^{-2}\cdot L^2 = 1$, which is what it carries.
The angles are the one place where the weights are not all one, and this is the component that shows it.

$G_{tt} = \left(\partial_r R\,\mathcal{K} + 2R\,\mathcal{J}\right)/(R^2\partial_r R)$ has a first term carrying $1/L^2$ and a second carrying $L\cdot L^{-1}/L^2 = 1/L^2$, which is the target for a doubly lowered curvature with two length indices.
Both terms carry it, which is the requirement, and a component that mixed a $\partial_t^2R$ with an $E'$ without a compensating $R$ would fail it here.

The $t$ geodesic equation is measured against $\ddot{t}$, which carries $L/\lambda^2$.
The term $\Gamma^t{}_{rr}\dot{r}^2$ carries $L^{-1}\left(L/\lambda\right)^2 = L/\lambda^2$, and the term $R\,\partial_t R\,\dot{\theta}^2$ carries $L\cdot 1\cdot\left(1/\lambda\right)^2 = L/\lambda^2$, since the chart velocity of a dimensionless angle is $1/\lambda$.
Both match.
On the other reading of the dots, in which $\dot{t}$ meant $dt/d\lambda$, the second of those two terms would come out short by a factor of $c$ against the first, which is the check the dimensional pass makes over the whole collection.

The Kretschmann scalar wants $1/L^4$, and each of its four squares delivers it: $(\partial_r R)^2\mathcal{K}^2/\left(R^4(\partial_r R)^2\right)$ is $L^{-4}$, and so are $2R^2\mathcal{J}^2$, $R^4(\partial_r\partial_t^2R)^2$ and $2R^2(\partial_r R)^2(\partial_t^2R)^2$ over the same denominator, since $\mathcal{J}$, $\partial_r\partial_t^2R$ and $\partial_t^2R$ carry $L^{-1}$, $L^{-2}$ and $L^{-1}$ respectively.

---

## Step 18. What the entry publishes, the declaration the checker needs, and the third derivative

The one system publishes one hundred and seventy three expressions:

- the line element, four metric components and four inverse metric components;
- nineteen Christoffel symbols with an upper index, from eleven distinct values, and nineteen with all three lowered, from seven;
- twenty four components of $R^\mu{}_{\nu\rho\sigma}$ from ten distinct values and twenty four of $R_{\mu\nu\rho\sigma}$ from six;
- four Ricci components in each of three variants, and the Ricci scalar;
- four Einstein components in each of three variants;
- twenty four Weyl components with an upper index from eight distinct values and twenty four with all four down from six;
- the Kretschmann scalar and four geodesic equations.

For the checker to read it, `DIMENSIONS` in `verify_metrics.py` needs one line:

```
("tolman_bondi", "comoving_synchronous"): {
    "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "R": "L", "E": "1",
},
```

The declaration `"t": "T"` is what tells the checker that $t$ is a time and so that the chart coordinate is $ct$, and every factor of $c$ in the comparison follows from it, as does the reading of $\partial_t$ as a chart derivative.
The two that a reader cannot guess are `"R": "L"` and `"E": "1"`.
The areal radius is a length and the shell label beside it is a length too, which is what makes $\partial_r R$ dimensionless; had $r$ been declared dimensionless instead, as a comoving label reasonably might be, $\partial_r R$ would carry a length and not one published component would balance.
The energy function has to be dimensionless because it is added to $1$.

The entry needs no line in `PARAMETER_RELATIONS`, and that is the most important thing to say about it.
Every published value here is an identity in a free $R(r,t)$ and a free $E(r)$.
The evolution equation of Step 11 is not assumed anywhere in the file, and if it were the Einstein block would print a single component instead of four and the entry would be claiming something on a constraint surface that the checker compares against free functions.
That is the same discipline the Alcubierre entry keeps when it declines to use $\partial_t f = -v_s\partial_x f$: the published components are what the geometry says, and the field equations are derived from them in Step 11 rather than fed into them.
The cost is that a reader meeting the Einstein block sees an anisotropic fluid rather than a dust, and the entry's own prose and this file are where that is resolved.

The index names in a published component have to be the coordinate names exactly as `coords` spells them, so `"\\theta"` and `"\\phi"` and not `"theta"` and `"phi"`.
Getting it wrong is loud in the sympy pass and silent in the dimensional one, which skips a component whose indices it cannot resolve, so a clean `--dimensions-only` run means nothing until the sympy comparison has confirmed that the index names resolve.

The third derivative is where this entry first met the checker's limits.
When the entry landed, `Reader._declare_parameter` in `verify_metrics.py` declared, for a parameter that is a function of several coordinates, every first partial derivative and every second one and nothing beyond, and `_tools/README.md` said that second derivatives are as far as any curvature tensor reaches.
That was true of every other entry in the collection and is false here, for the reason Step 6 gives: the radial metric component already carries one derivative of $R$, so a curvature, which is two derivatives of the metric, reaches $\partial_r\partial_t^2R$.
Seventy of the one hundred and seventy three published expressions name that third derivative:

- the four components of each Riemann variant in the plane of $t$ and $r$, eight in all;
- the $tt$ and $rr$ components of all three Ricci variants, six in all, and the Ricci scalar;
- the $\theta\theta$ and $\phi\phi$ components of all three Einstein variants, six in all;
- all twenty four components of each Weyl variant, forty eight in all;
- the Kretschmann scalar.

The script reported each of them `UNCHECKED` rather than passing it in silence, once from the dimensional pass and once from the sympy pass, one hundred and forty lines for seventy expressions, and until 22 September 2026 they were the last expressions in the collection it could not read.
They were checked in the meantime by a script outside the repository that patched a third loop into `_declare_parameter` and called the script's own `check_system` on this entry, and it disagreed with none of them.

Since 22 September 2026 the reader has no fixed order at all.
`Reader._declare_partials` declares a partial derivative when a published value names one, at whatever order it is written, provided the function is declared and every coordinate it is differentiated along is one the function is declared to depend on.
A third loop would have closed this entry and left the next entry one order further on to meet the same wall, while declaring on demand costs nothing for the orders nobody writes and still turns a typo into an error rather than a new symbol: $\partial_\theta R$ is refused with the reason that $R$ is declared a function of $r$ and $t$ only.
The spelling does not matter either, because mixed partials commute and every spelling of $\partial_r\partial_t^2R$ is read as the one object.

The whole system, all one hundred and seventy three expressions, is now compared against sympy in the ordinary way and agrees:

```
python3 _tools/derivations/verify_metrics.py \
    --system tolman_bondi/comoving_synchronous
```

reports no disagreement, no dimensional failure and nothing `UNCHECKED`, in about two seconds, and exits zero.
That the seventy are really compared and not merely parsed was confirmed by breaking them: one run doubled the third derivative term of the Kretschmann scalar and the third derivative in one Weyl component, and came back with exactly those two disagreements, and a second run spelled the third derivative $\partial_t\partial_r\partial_t R$ throughout the Weyl block and came back clean.
