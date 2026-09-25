# The curvature components missing from Ellis-Bronnikov, and the sign of its radial geodesic

The Ellis-Bronnikov wormhole needed three corrections, and each is worked out by hand.
Only the corrected components are derived in full, and not the whole spacetime as elsewhere.
No component of the Ellis-Bronnikov wormhole had a wrong value.
Its metric, its eighteen Christoffel symbols, its Riemann tensor with the first index up, its Ricci tensor with both indices down, with one up and with both up, its Ricci scalar, its Kretschmann scalar, its Einstein tensor in the same three index positions and its Weyl tensor with the first index up were all correct before the corrections, and none of them changes.

What was wrong was an omission and a sign.

1. Of the twelve nonzero components of the Riemann tensor with every index lowered, ten were given and two were left out.
2. Of the twenty four nonzero components of the Weyl tensor with every index lowered, fourteen were given and ten were left out.
3. The radial geodesic equation had the sign opposite to the one its own $\Gamma^r{}_{\theta\theta}$ forces.

Every component left out of a list of nonzero components is read as a claim that it vanishes, so the two omissions were not silence.
They asserted that twelve components of curvature are zero where they are not, and the checker reports an omission for exactly that reason.

We compute every one of the twelve added components from the connection and from the definition of the Weyl tensor, by hand, so that the captain can check any one of them without running anything.
The companion script `_tools/derivations/verify_metrics.py` does the same work in sympy and compares it against the stored components; on this one system it runs in about three seconds, which makes Ellis-Bronnikov one of the cheapest spacetimes to check again.

Two independent confirmations are carried alongside the algebra.
Minkowski in spherical coordinates is this spacetime at $\ell = 0$ and settles the geodesic sign on its own, in Step 9.
Morris-Thorne's wormhole, worked out the same day, contains this spacetime as the special case $\Phi = 0$, $r(l) = \sqrt{l^2+\ell^2}$, and specialised in full in Step 11 it reproduces every component of this spacetime, including all twelve of the added ones.

The inverse metric, the one set of components still missing, was derived on 22 September 2026, after the rest, in Step 13.
Every other coordinate system of every spacetime had one, and its absence was the last line the checker reported `UNCHECKED` for want of anything to check.

---

## Step 1. Conventions

The signature is $(-,+,+,+)$.
The connection is the Levi-Civita one,

$$\Gamma^\mu{}_{\nu\rho} = \tfrac{1}{2}g^{\mu\alpha}\left(\partial_\nu g_{\alpha\rho} + \partial_\rho g_{\alpha\nu} - \partial_\alpha g_{\nu\rho}\right),$$

with the lowered symbol $\Gamma_{\mu\nu\rho} = g_{\mu\alpha}\Gamma^\alpha{}_{\nu\rho}$, and the Riemann tensor is

$$R^\mu{}_{\nu\rho\sigma} = \partial_\rho \Gamma^\mu{}_{\nu\sigma} - \partial_\sigma \Gamma^\mu{}_{\nu\rho} + \Gamma^\mu{}_{\rho\lambda}\Gamma^\lambda{}_{\nu\sigma} - \Gamma^\mu{}_{\sigma\lambda}\Gamma^\lambda{}_{\nu\rho},$$

which fixes the sign of every Riemann component.

The Ricci tensor is the standard contraction on the first lower index,

$$R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu},$$

the contraction every spacetime has used since 18 September 2026.
Here the choice is visible in exactly one slot, because this spacetime has a single nonzero Ricci component, and that component already has the right sign: $R_{rr} = -2\ell^2/(r^2+\ell^2)^2$ is negative, and with it $G_{tt} = -\ell^2/(r^2+\ell^2)^2$ is negative, which is the negative energy density the throat has to be held open with.
The conventions of the spacetime say so, and none of the three corrections changes it.

The Weyl tensor is built from that same Ricci tensor and that same Ricci scalar through

$$C_{\mu\nu\rho\sigma} = R_{\mu\nu\rho\sigma} - \frac{1}{n-2}\left(g_{\mu\rho}R_{\sigma\nu} - g_{\mu\sigma}R_{\rho\nu} - g_{\nu\rho}R_{\sigma\mu} + g_{\nu\sigma}R_{\rho\mu}\right) + \frac{R}{(n-1)(n-2)}\left(g_{\mu\rho}g_{\nu\sigma} - g_{\mu\sigma}g_{\nu\rho}\right),$$

which in four dimensions is

$$C_{\mu\nu\rho\sigma} = R_{\mu\nu\rho\sigma} - \tfrac{1}{2}\left(g_{\mu\rho}R_{\sigma\nu} - g_{\mu\sigma}R_{\rho\nu} - g_{\nu\rho}R_{\sigma\mu} + g_{\nu\sigma}R_{\rho\mu}\right) + \tfrac{1}{6}R\left(g_{\mu\rho}g_{\nu\sigma} - g_{\mu\sigma}g_{\nu\rho}\right).$$

This is not a vacuum spacetime.
Its Ricci tensor does not vanish, so its Weyl tensor is not its Riemann tensor, and every Weyl component in Step 7 is computed from that formula rather than copied from anywhere.
The difference is not a decoration: of the ten added components, the four carrying a pair of $r$ indices come out a third of the corresponding Riemann component, the two purely angular ones come out two thirds of theirs, and the four carrying a pair of $t$ indices are nonzero in Weyl where Riemann is identically zero.

Factors of $c$ and $G$ are explicit.
The chart is the one whose zeroth coordinate is

$$x^0 = ct,$$

which is why $g_{tt} = -1$ stands against a line element whose time term is $-c^2dt^2$.
Because the rescaling $t \to ct$ is linear with constant coefficients, a component in this chart is the component taken with the bare coordinate, multiplied by $c$ once for every upper time index and divided by $c$ once for every lower one, and the Christoffel symbols follow the same rule, since the inhomogeneous term in their transformation law carries a second derivative of the coordinate change and that vanishes for a linear one.

We compute everything directly in the $x^0$ chart, so no factor of $c$ appears anywhere.
That is not an accident but a property of this metric: $g_{tt} = -1$ is a constant, so every derivative of it vanishes, and $t$ never meets a metric function.
The dimensions of the twelve added components in that chart are checked in Step 10.

The dots in the geodesic equations are velocities of that same chart, so $\dot{t}$ means $d(ct)/d\lambda$ and not $dt/d\lambda$, and the equations are

$$\ddot{x}^\mu + \Gamma^\mu{}_{\nu\rho}\dot{x}^\nu\dot{x}^\rho = 0$$

with the same Christoffel symbols as Step 3.
The third correction is to the radial one of these equations, and it is made in Step 9.

We use one abbreviation throughout, purely to keep the lines short:

$$\rho^2 \equiv r^2 + \ell^2.$$

It is shorthand for the derivation only, and never enters the components themselves.
Every stored component writes $r^2+\ell^2$ out wherever it appears, because the checker reads each one symbol by symbol and has no way to be told what an abbreviation means.
Note that $\rho$ is the areal radius of this geometry and not a radial coordinate: the sphere labelled $r$ has area $4\pi\rho^2$, and $\rho$ has its minimum $\ell$ at the throat $r = 0$ and grows without bound on both sides.
The one derivative needed constantly is

$$\partial_r \rho^2 = 2r.$$

---

## Step 2. The line element and the metric

The solution is

$$ds^2 = -c^2dt^2 + dr^2 + \left(r^2+\ell^2\right)\left(d\theta^2 + \sin^2\theta\,d\phi^2\right),$$

which in the chart coordinate $x^0 = ct$ is

$$ds^2 = -(dx^0)^2 + dr^2 + \rho^2\left(d\theta^2 + \sin^2\theta\,d\phi^2\right).$$

The radial coordinate is the proper radial distance and runs over the whole line, $r \in (-\infty,\infty)$.
There is no center and no singularity.
The two ends $r \to +\infty$ and $r \to -\infty$ are separate asymptotically flat regions joined through the throat at $r = 0$, where the areal radius takes its minimum value $\ell$.

The metric is diagonal, and in the order $(t,r,\theta,\phi)$ its components are

$$g_{tt} = -1, \qquad g_{rr} = 1, \qquad g_{\theta\theta} = \rho^2, \qquad g_{\phi\phi} = \rho^2\sin^2\theta,$$

with the inverse

$$g^{tt} = -1, \qquad g^{rr} = 1, \qquad g^{\theta\theta} = \frac{1}{\rho^2}, \qquad g^{\phi\phi} = \frac{1}{\rho^2\sin^2\theta}.$$

When the three corrections were made the inverse metric was missing, and the checker reported one `UNCHECKED` line against this system for it.
It is in place now, derived in Step 13.

Two facts about this metric do most of the work and are worth stating once.

**Nothing depends on $t$, and $g_{tt}$ is constant.**
So no Christoffel symbol carries a $t$ index at all: $\Gamma^t{}_{\mu\nu} = 0$ and $\Gamma^\mu{}_{t\nu} = 0$ for every $\mu$ and $\nu$.
Every Riemann component with a $t$ index is therefore built entirely out of terms that are zero, and the whole $t$ row and $t$ column of the Riemann tensor vanish identically.
That is the reason the Riemann tensor has only twelve nonzero components where a generic spherically symmetric metric has twenty four, and it is also the reason the Weyl tensor has components the Riemann tensor does not: the trace terms of Step 7 do not care that $\Gamma$ has no $t$ index, only that $g_{tt} = -1$ is not zero.

**The angular part is $\rho^2$ times the unit sphere, and $\rho$ depends on $r$ alone.**
So the angular sector behaves exactly as it does in flat spherical coordinates, with $r$ replaced by $\rho$ in the metric but not in the derivative: $\partial_r\rho = r/\rho$, which is smaller than one in magnitude everywhere and vanishes at the throat, and that single mismatch is the whole curvature of this spacetime.

---

## Step 3. The connection

Every curvature component is built from these Christoffel symbols, which were right and are restated unchanged.
The checker agrees with all nine symbols with the first index up and all nine with every index lowered, and none was touched by this correction.
The symmetric partner $\Gamma^\mu{}_{\rho\nu} = \Gamma^\mu{}_{\nu\rho}$ is left implicit, though the checker is given both orderings.

$$\Gamma^r{}_{\theta\theta} = -r, \qquad \Gamma^r{}_{\phi\phi} = -r\sin^2\theta,$$

$$\Gamma^\theta{}_{r\theta} = \Gamma^\phi{}_{r\phi} = \frac{r}{\rho^2}, \qquad \Gamma^\theta{}_{\phi\phi} = -\sin\theta\cos\theta, \qquad \Gamma^\phi{}_{\theta\phi} = \cot\theta.$$

Two of them are worth deriving, since the whole of Step 5 and the whole of Step 9 turn on them.

$$\Gamma^r{}_{\theta\theta} = \tfrac{1}{2}g^{rr}\left(-\partial_r g_{\theta\theta}\right) = \tfrac{1}{2}\left(-\partial_r \rho^2\right) = -r,$$

$$\Gamma^\theta{}_{r\theta} = \tfrac{1}{2}g^{\theta\theta}\,\partial_r g_{\theta\theta} = \frac{1}{2\rho^2}\,\partial_r\rho^2 = \frac{r}{\rho^2}.$$

The first is the flat space value exactly, with no $\ell$ anywhere in it, because $g_{rr} = 1$ and $\partial_r\rho^2 = \partial_r r^2$.
The second is not: $r/\rho^2$ is the flat space $1/r$ only when $\ell$ vanishes, and the gap between $1/r$ and $r/\rho^2$ is where the curvature of this spacetime comes from.

Three features of these symbols do all the work in the curvature and the geodesics.
There is no $\Gamma$ with a $t$ index of any kind, for the reason given in Step 2.
There is no $\Gamma^\mu{}_{rr}$ of any kind, because $g_{rr} = 1$ is constant, so the curve of constant $\theta$ and $\phi$ is a geodesic and the radial direction is as simple as it can be.
And $\Gamma^\theta{}_{r\theta} = \Gamma^\phi{}_{r\phi}$, which is why the $\theta$ and $\phi$ rows of the Riemann tensor come out the same up to a factor of $\sin^2\theta$.

---

## Step 4. The mixed Riemann tensor, restated

The tensor $R^\mu{}_{\nu\rho\sigma}$ has twelve nonzero components, and the checker agrees with all twelve.
They are not being changed, but two of them are the input to Step 5, so all twelve are restated in full.

$$R^r{}_{\theta r\theta} = -\frac{\ell^2}{\rho^2}, \qquad R^r{}_{\theta\theta r} = \frac{\ell^2}{\rho^2}, \qquad R^r{}_{\phi r\phi} = -\frac{\ell^2}{\rho^2}\sin^2\theta, \qquad R^r{}_{\phi\phi r} = \frac{\ell^2}{\rho^2}\sin^2\theta,$$

$$R^\theta{}_{rr\theta} = \frac{\ell^2}{\rho^4}, \qquad R^\theta{}_{r\theta r} = -\frac{\ell^2}{\rho^4}, \qquad R^\theta{}_{\phi\theta\phi} = \frac{\ell^2}{\rho^2}\sin^2\theta, \qquad R^\theta{}_{\phi\phi\theta} = -\frac{\ell^2}{\rho^2}\sin^2\theta,$$

$$R^\phi{}_{rr\phi} = \frac{\ell^2}{\rho^4}, \qquad R^\phi{}_{r\phi r} = -\frac{\ell^2}{\rho^4}, \qquad R^\phi{}_{\theta\theta\phi} = -\frac{\ell^2}{\rho^2}, \qquad R^\phi{}_{\theta\phi\theta} = \frac{\ell^2}{\rho^2}.$$

Every nonzero component of this tensor carries $\ell^2$ as a factor, which is the statement that $\ell = 0$ is flat space.

The last two are the input to Step 5, so we compute one of them rather than quote it.
Take $\mu = \phi$, $\nu = \theta$, $\rho = \theta$, $\sigma = \phi$ in the definition:

$$R^\phi{}_{\theta\theta\phi} = \partial_\theta\Gamma^\phi{}_{\theta\phi} - \partial_\phi\Gamma^\phi{}_{\theta\theta} + \Gamma^\phi{}_{\theta\lambda}\Gamma^\lambda{}_{\theta\phi} - \Gamma^\phi{}_{\phi\lambda}\Gamma^\lambda{}_{\theta\theta}.$$

The four terms, one at a time.

The first is $\partial_\theta\cot\theta = -\csc^2\theta = -1/\sin^2\theta$.

The second vanishes, since $\Gamma^\phi{}_{\theta\theta} = 0$.

In the third, the only $\lambda$ for which $\Gamma^\phi{}_{\theta\lambda}$ is nonzero is $\lambda = \phi$, giving

$$\Gamma^\phi{}_{\theta\phi}\Gamma^\phi{}_{\theta\phi} = \cot^2\theta.$$

In the fourth, $\Gamma^\phi{}_{\phi\lambda}$ is nonzero for $\lambda = r$ and $\lambda = \theta$, but $\Gamma^\lambda{}_{\theta\theta}$ is nonzero only for $\lambda = r$, so only $\lambda = r$ survives:

$$-\Gamma^\phi{}_{\phi r}\Gamma^r{}_{\theta\theta} = -\frac{r}{\rho^2}\cdot(-r) = \frac{r^2}{\rho^2}.$$

Adding them,

$$R^\phi{}_{\theta\theta\phi} = -\frac{1}{\sin^2\theta} + \frac{\cos^2\theta}{\sin^2\theta} + \frac{r^2}{\rho^2} = -\frac{1-\cos^2\theta}{\sin^2\theta} + \frac{r^2}{\rho^2} = -1 + \frac{r^2}{r^2+\ell^2},$$

so that

$$\boxed{R^\phi{}_{\theta\theta\phi} = -\frac{\ell^2}{r^2+\ell^2}, \qquad R^\phi{}_{\theta\phi\theta} = +\frac{\ell^2}{r^2+\ell^2},}$$

the second by the antisymmetry $R^\mu{}_{\nu\rho\sigma} = -R^\mu{}_{\nu\sigma\rho}$ in the last two indices.
The whole angular dependence cancelled, as it must: $-\csc^2\theta + \cot^2\theta = -1$ exactly.
What is left is the deficit $1 - r^2/\rho^2$, which is $\ell^2/\rho^2$, and that single expression is the curvature of the throat.
At $\ell = 0$ it is zero, which is flat space, and at the throat $r = 0$ it is $1$, its largest value anywhere.

---

## Step 5. The two Riemann components that were missing

### The lowering

The fully lowered Riemann tensor is

$$R_{\mu\nu\rho\sigma} = g_{\mu\alpha}R^\alpha{}_{\nu\rho\sigma},$$

and the metric is diagonal, so the sum has one term and the lowering multiplies by a single metric component chosen by the first index.
For the three indices that appear in a nonzero component that is

$$g_{rr} = 1, \qquad g_{\theta\theta} = \rho^2, \qquad g_{\phi\phi} = \rho^2\sin^2\theta.$$

So an $r$ row is unchanged by the lowering, a $\theta$ row is multiplied by $\rho^2$ and a $\phi$ row by $\rho^2\sin^2\theta$.
Applying that to the twelve of Step 4 gives the twelve fully lowered components, and the two missing ones are the last pair.

### The two that were missing

$$R_{\phi\theta\theta\phi} = g_{\phi\phi}R^\phi{}_{\theta\theta\phi} = \rho^2\sin^2\theta\cdot\left(-\frac{\ell^2}{\rho^2}\right) = -\ell^2\sin^2\theta,$$

$$R_{\phi\theta\phi\theta} = g_{\phi\phi}R^\phi{}_{\theta\phi\theta} = \rho^2\sin^2\theta\cdot\frac{\ell^2}{\rho^2} = +\ell^2\sin^2\theta.$$

$$\boxed{R_{\phi\theta\theta\phi} = -\ell^2\sin^2\theta, \qquad R_{\phi\theta\phi\theta} = +\ell^2\sin^2\theta.}$$

The factor $\rho^2$ cancels completely.
That is why these two, together with the partners $R_{\theta\phi\theta\phi}$ and $R_{\theta\phi\phi\theta}$ that were never missing, are the four components of the fully lowered Riemann tensor that are pure $\ell^2$ with no $\rho$ left in them.

### A one line check by antisymmetry

The two partners that were never missing are

$$R_{\theta\phi\theta\phi} = \ell^2\sin^2\theta, \qquad R_{\theta\phi\phi\theta} = -\ell^2\sin^2\theta,$$

and the Riemann tensor is antisymmetric in its first pair,

$$R_{\mu\nu\rho\sigma} = -R_{\nu\mu\rho\sigma},$$

so $R_{\phi\theta\theta\phi} = -R_{\theta\phi\theta\phi}$ follows without touching the connection.
Reading off the value already given yields $-\ell^2\sin^2\theta$, which is what the lowering gave.
The two added components were therefore forced by two already given, and the components contradicted themselves: one half of an antisymmetric pair was given as nonzero and the other half, by omission, as zero.

### The census

The fully lowered Riemann tensor of this spacetime has twelve nonzero components, three pairs of indices carrying four each, and all twelve stand in the table with the two added ones marked.

| first pair | components |
| --- | --- |
| $\{r,\theta\}$ | $R_{r\theta r\theta} = -\dfrac{\ell^2}{\rho^2}$, $\;R_{r\theta\theta r} = \dfrac{\ell^2}{\rho^2}$, $\;R_{\theta r r\theta} = \dfrac{\ell^2}{\rho^2}$, $\;R_{\theta r\theta r} = -\dfrac{\ell^2}{\rho^2}$ |
| $\{r,\phi\}$ | $R_{r\phi r\phi} = -\dfrac{\ell^2}{\rho^2}\sin^2\theta$, $\;R_{r\phi\phi r} = \dfrac{\ell^2}{\rho^2}\sin^2\theta$, $\;R_{\phi r r\phi} = \dfrac{\ell^2}{\rho^2}\sin^2\theta$, $\;R_{\phi r\phi r} = -\dfrac{\ell^2}{\rho^2}\sin^2\theta$ |
| $\{\theta,\phi\}$ | $R_{\theta\phi\theta\phi} = \ell^2\sin^2\theta$, $\;R_{\theta\phi\phi\theta} = -\ell^2\sin^2\theta$, **$R_{\phi\theta\theta\phi} = -\ell^2\sin^2\theta$**, **$R_{\phi\theta\phi\theta} = \ell^2\sin^2\theta$** |

The first two rows were complete before the correction.
The third carried only its $\theta\phi$ half.
Everything with a $t$ index is zero, for the reason given in Step 2, so there are no further rows.

---

## Step 6. The Ricci tensor and the Ricci scalar

The Weyl tensor of Step 7 needs both, so we restate and confirm them.
Both were already given and both are right; neither is changed.

Contracting on the first lower index,

$$R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu},$$

and the only nonzero component is the radial one.
Reading the four $\alpha$ values off the components of Step 4 for $\mu = \nu = r$,

$$R_{rr} = R^t{}_{rtr} + R^r{}_{rrr} + R^\theta{}_{r\theta r} + R^\phi{}_{r\phi r} = 0 + 0 - \frac{\ell^2}{\rho^4} - \frac{\ell^2}{\rho^4},$$

so

$$R_{rr} = -\frac{2\ell^2}{\rho^4} = -\frac{2\ell^2}{(r^2+\ell^2)^2},$$

which is the value already given.
That the $\theta\theta$ and $\phi\phi$ components vanish is the cancellation that makes this solution what it is: for $\mu = \nu = \theta$ the surviving terms are

$$R_{\theta\theta} = R^r{}_{\theta r\theta} + R^\phi{}_{\theta\phi\theta} = -\frac{\ell^2}{\rho^2} + \frac{\ell^2}{\rho^2} = 0,$$

using $R^\phi{}_{\theta\phi\theta}$ from Step 4, the component whose lowered partner Step 5 had to add.

The Ricci scalar is then

$$R = g^{\mu\nu}R_{\mu\nu} = g^{rr}R_{rr} = -\frac{2\ell^2}{(r^2+\ell^2)^2},$$

numerically equal to $R_{rr}$ because $g^{rr} = 1$ and there is nothing else to sum.
It is worth flagging that the two are the same number here, since $R$ and $R_{rr}$ both appear in Step 7 and a reader checking a cancellation there needs to know which of them a given $-2\ell^2/\rho^4$ came from.

Two quantities recur in the Weyl tensor:

$$R_{rr} = R = -\frac{2\ell^2}{\rho^4}, \qquad \frac{R}{6} = -\frac{\ell^2}{3\rho^4}.$$

---

## Step 7. The ten Weyl components that were missing

### What the trace terms can do

With the Weyl tensor as written in Step 1,

$$C_{\mu\nu\rho\sigma} = R_{\mu\nu\rho\sigma} - \tfrac{1}{2}\left(g_{\mu\rho}R_{\sigma\nu} - g_{\mu\sigma}R_{\rho\nu} - g_{\nu\rho}R_{\sigma\mu} + g_{\nu\sigma}R_{\rho\mu}\right) + \tfrac{1}{6}R\left(g_{\mu\rho}g_{\nu\sigma} - g_{\mu\sigma}g_{\nu\rho}\right),$$

the middle bracket has only one way to be nonzero.
The Ricci tensor has the single component $R_{rr}$, so a term of that bracket survives only when its Ricci factor is $R_{rr}$, which means two of the four indices are $r$; and the metric is diagonal, so the metric factor beside it needs its own two indices equal.
Together that pins the bracket down: it contributes only to components whose index set is $\{r,r,X,X\}$ for a single other coordinate $X$.

The last bracket is easier.
The metric is diagonal, so $g_{\mu\rho}g_{\nu\sigma}$ survives only when $\mu = \rho$ and $\nu = \sigma$, and $g_{\mu\sigma}g_{\nu\rho}$ only when $\mu = \sigma$ and $\nu = \rho$.
It therefore contributes to every component of the form $C_{XYXY}$ and $C_{XYYX}$, whatever $X$ and $Y$ are, and in particular to the ones carrying a $t$ index, where both the Riemann tensor and the middle bracket are zero.

That last bracket is the whole reason the Weyl tensor of this spacetime has twenty four nonzero components where the Riemann tensor has twelve.
The Weyl tensor of this spacetime is nonzero in the time direction even though its Riemann tensor is not, and the only thing carrying it there is $\tfrac{1}{6}R\,g_{tt}g_{XX}$.
It is the trace of the curvature, not the curvature itself, that reaches the $t$ row.

### The four with a $t$ index

Take $\mu = \theta$, $\nu = t$, $\rho = t$, $\sigma = \theta$.
The Riemann term is zero by Step 2.
In the middle bracket, $g_{\theta t}R_{\theta t}$, $g_{\theta\theta}R_{tt}$, $g_{tt}R_{\theta\theta}$ and $g_{t\theta}R_{t\theta}$ are all zero, the first and last because the metric is diagonal and the other two because the Ricci tensor has only $R_{rr}$.
So the whole component is the last bracket:

$$C_{\theta t t\theta} = \tfrac{1}{6}R\left(g_{\theta t}g_{t\theta} - g_{\theta\theta}g_{tt}\right) = \tfrac{1}{6}R\left(0 - \rho^2\cdot(-1)\right) = \tfrac{1}{6}R\,\rho^2 = -\frac{2\ell^2}{6\rho^4}\cdot\rho^2,$$

$$\boxed{C_{\theta t t\theta} = -\frac{\ell^2}{3\rho^2} = -\frac{\ell^2}{3(r^2+\ell^2)}, \qquad C_{\theta t\theta t} = +\frac{\ell^2}{3(r^2+\ell^2)}.}$$

This is the component the audit of 18 September 2026 named, and it is the value sympy gives.

The $\phi$ version is the same computation with $g_{\phi\phi} = \rho^2\sin^2\theta$ in place of $g_{\theta\theta} = \rho^2$:

$$C_{\phi t t\phi} = \tfrac{1}{6}R\left(-g_{\phi\phi}g_{tt}\right) = \tfrac{1}{6}R\,\rho^2\sin^2\theta,$$

$$\boxed{C_{\phi t t\phi} = -\frac{\ell^2}{3\rho^2}\sin^2\theta, \qquad C_{\phi t\phi t} = +\frac{\ell^2}{3\rho^2}\sin^2\theta.}$$

### The four with a pair of $r$ indices

Take $\mu = \theta$, $\nu = r$, $\rho = r$, $\sigma = \theta$.
Here all three pieces contribute, which happens only for the eight components whose index pair is $\{r,\theta\}$ or $\{r,\phi\}$, four of them already given and four added.
On the $\{t,\theta\}$ and $\{t,\phi\}$ pairs only the Ricci scalar term is left, on $\{t,r\}$ and $\{\theta,\phi\}$ exactly two of the three survive.

The Riemann term is $R_{\theta r r\theta} = \ell^2/\rho^2$, already given.

The middle bracket is

$$-\tfrac{1}{2}\left(g_{\theta r}R_{\theta r} - g_{\theta\theta}R_{rr} - g_{rr}R_{\theta\theta} + g_{r\theta}R_{r\theta}\right) = -\tfrac{1}{2}\left(0 - \rho^2\cdot\left(-\frac{2\ell^2}{\rho^4}\right) - 0 + 0\right) = -\frac{\ell^2}{\rho^2}.$$

The last bracket is

$$\tfrac{1}{6}R\left(g_{\theta r}g_{r\theta} - g_{\theta\theta}g_{rr}\right) = \tfrac{1}{6}R\left(0 - \rho^2\right) = -\frac{1}{6}\cdot\left(-\frac{2\ell^2}{\rho^4}\right)\rho^2 = +\frac{\ell^2}{3\rho^2}.$$

Adding the three,

$$C_{\theta r r\theta} = \frac{\ell^2}{\rho^2} - \frac{\ell^2}{\rho^2} + \frac{\ell^2}{3\rho^2},$$

$$\boxed{C_{\theta r r\theta} = +\frac{\ell^2}{3(r^2+\ell^2)}, \qquad C_{\theta r\theta r} = -\frac{\ell^2}{3(r^2+\ell^2)}.}$$

The Riemann term and the middle bracket cancelled each other exactly, and what survives is the Ricci scalar term alone.
That is the sense in which these Weyl components could not have been copies of the Riemann components even in the slots where both are nonzero: here $C_{\theta rr\theta}$ is $\tfrac{1}{3}$ of $R_{\theta rr\theta} = \ell^2/\rho^2$, while in the previous family $C_{\theta t t\theta}$ sits in a slot where Riemann is identically zero.

The $\phi$ version again carries the extra $\sin^2\theta$ throughout, with $R_{\phi rr\phi} = (\ell^2/\rho^2)\sin^2\theta$, $g_{\phi\phi} = \rho^2\sin^2\theta$:

$$C_{\phi r r\phi} = \frac{\ell^2}{\rho^2}\sin^2\theta - \frac{\ell^2}{\rho^2}\sin^2\theta + \frac{\ell^2}{3\rho^2}\sin^2\theta,$$

$$\boxed{C_{\phi r r\phi} = +\frac{\ell^2}{3\rho^2}\sin^2\theta, \qquad C_{\phi r\phi r} = -\frac{\ell^2}{3\rho^2}\sin^2\theta.}$$

### The two purely angular ones

Take $\mu = \phi$, $\nu = \theta$, $\rho = \theta$, $\sigma = \phi$.
This is the one family whose Riemann input is itself a component added in Step 5, which is why the two corrections cannot be checked apart.

The Riemann term is $R_{\phi\theta\theta\phi} = -\ell^2\sin^2\theta$, from Step 5.

The middle bracket vanishes: it would need $R_{\theta\theta}$ or $R_{\phi\phi}$, both of which are zero, or an off diagonal metric component, and there are none.

The last bracket is

$$\tfrac{1}{6}R\left(g_{\phi\theta}g_{\theta\phi} - g_{\phi\phi}g_{\theta\theta}\right) = \tfrac{1}{6}R\left(0 - \rho^2\sin^2\theta\cdot\rho^2\right) = -\frac{1}{6}\cdot\left(-\frac{2\ell^2}{\rho^4}\right)\rho^4\sin^2\theta = +\frac{\ell^2}{3}\sin^2\theta.$$

Adding the two,

$$C_{\phi\theta\theta\phi} = -\ell^2\sin^2\theta + \frac{\ell^2}{3}\sin^2\theta,$$

$$\boxed{C_{\phi\theta\theta\phi} = -\frac{2}{3}\ell^2\sin^2\theta, \qquad C_{\phi\theta\phi\theta} = +\frac{2}{3}\ell^2\sin^2\theta.}$$

This is the check that ties the two corrections together.
The component $C_{\theta\phi\theta\phi} = \tfrac{2}{3}\ell^2\sin^2\theta$ was already given, and it is the same computation with $R_{\theta\phi\theta\phi} = \ell^2\sin^2\theta$, never missing, in place of the added $R_{\phi\theta\theta\phi}$, and the value that comes out is minus it, as the first pair antisymmetry requires.
Had the Riemann component of Step 5 been omitted from the calculation rather than merely from the list of components, this Weyl component would have come out $+\tfrac{1}{3}\ell^2\sin^2\theta$ and not matched its own antisymmetric partner.

### The census

The fully lowered Weyl tensor has twenty four nonzero components, six pairs of indices carrying four each.
All twenty four stand in the table, with the ten added ones in bold.

| pair | already given | added |
| --- | --- | --- |
| $\{t,r\}$ | $C_{trtr} = -\dfrac{2\ell^2}{3\rho^4}$, $\;C_{trrt} = \dfrac{2\ell^2}{3\rho^4}$, $\;C_{rttr} = \dfrac{2\ell^2}{3\rho^4}$, $\;C_{rtrt} = -\dfrac{2\ell^2}{3\rho^4}$ | none, this pair was complete |
| $\{t,\theta\}$ | $C_{t\theta t\theta} = \dfrac{\ell^2}{3\rho^2}$, $\;C_{t\theta\theta t} = -\dfrac{\ell^2}{3\rho^2}$ | **$C_{\theta t t\theta} = -\dfrac{\ell^2}{3\rho^2}$**, **$C_{\theta t\theta t} = \dfrac{\ell^2}{3\rho^2}$** |
| $\{t,\phi\}$ | $C_{t\phi t\phi} = \dfrac{\ell^2}{3\rho^2}\sin^2\theta$, $\;C_{t\phi\phi t} = -\dfrac{\ell^2}{3\rho^2}\sin^2\theta$ | **$C_{\phi t t\phi} = -\dfrac{\ell^2}{3\rho^2}\sin^2\theta$**, **$C_{\phi t\phi t} = \dfrac{\ell^2}{3\rho^2}\sin^2\theta$** |
| $\{r,\theta\}$ | $C_{r\theta r\theta} = -\dfrac{\ell^2}{3\rho^2}$, $\;C_{r\theta\theta r} = \dfrac{\ell^2}{3\rho^2}$ | **$C_{\theta r r\theta} = \dfrac{\ell^2}{3\rho^2}$**, **$C_{\theta r\theta r} = -\dfrac{\ell^2}{3\rho^2}$** |
| $\{r,\phi\}$ | $C_{r\phi r\phi} = -\dfrac{\ell^2}{3\rho^2}\sin^2\theta$, $\;C_{r\phi\phi r} = \dfrac{\ell^2}{3\rho^2}\sin^2\theta$ | **$C_{\phi r r\phi} = \dfrac{\ell^2}{3\rho^2}\sin^2\theta$**, **$C_{\phi r\phi r} = -\dfrac{\ell^2}{3\rho^2}\sin^2\theta$** |
| $\{\theta,\phi\}$ | $C_{\theta\phi\theta\phi} = \dfrac{2}{3}\ell^2\sin^2\theta$, $\;C_{\theta\phi\phi\theta} = -\dfrac{2}{3}\ell^2\sin^2\theta$ | **$C_{\phi\theta\theta\phi} = -\dfrac{2}{3}\ell^2\sin^2\theta$**, **$C_{\phi\theta\phi\theta} = \dfrac{2}{3}\ell^2\sin^2\theta$** |

---

## Step 8. The shape of both omissions

Both omissions have the same shape, and it is worth naming.

Every one of the twelve added components is the first pair swap $\mu \leftrightarrow \nu$ of a component already given.
The components with every index lowered were complete on the index pair $\{t,r\}$ in Weyl and on $\{r,\theta\}$ and $\{r,\phi\}$ in Riemann, and carried exactly one of the two orderings for every other pair.
Nothing was ever wrong; half of each list was simply never written down.

That is a consistent failure mode for components listed by hand and worth watching for in other spacetimes.
A tensor antisymmetric in its first pair has twice as many nonzero components as a reader counts when working through the independent ones, and a list of nonzero components has to carry both, since the checker reads anything absent as zero.
The Ellis-Bronnikov Riemann and Weyl tensors with the first index up did not have the problem, because there the first index is upstairs and the second down, and the swap is not a symmetry, so nothing invited the reader to think one of the pair stood for both.

---

## Step 9. The radial geodesic

### The equation as it stood

The radial equation stood as

$$\ddot{r} = -r\dot{\theta}^2 - r\sin^2\theta\,\dot{\phi}^2$$

beside its own $\Gamma^r{}_{\theta\theta} = -r$ and $\Gamma^r{}_{\phi\phi} = -r\sin^2\theta$.
Those two connection symbols were right.
The equation was not.

### What the connection gives

The geodesic equation is

$$\ddot{x}^\mu + \Gamma^\mu{}_{\nu\rho}\dot{x}^\nu\dot{x}^\rho = 0,$$

and for $\mu = r$ the only nonzero symbols with an upper $r$ are $\Gamma^r{}_{\theta\theta}$ and $\Gamma^r{}_{\phi\phi}$, so

$$\ddot{r} + \Gamma^r{}_{\theta\theta}\dot{\theta}^2 + \Gamma^r{}_{\phi\phi}\dot{\phi}^2 = 0,$$

$$\ddot{r} - r\dot{\theta}^2 - r\sin^2\theta\,\dot{\phi}^2 = 0,$$

$$\boxed{\ddot{r} = r\dot{\theta}^2 + r\sin^2\theta\,\dot{\phi}^2.}$$

The equation as it stood was the negative of this.
The other three geodesic equations are right and are not touched: $\ddot{t} = 0$, because no connection symbol carries a $t$ index, and the $\theta$ and $\phi$ equations already solve $\ddot{x}^\mu = -\Gamma^\mu{}_{\nu\rho}\dot{x}^\nu\dot{x}^\rho$ correctly against the same Christoffel symbols.

### The first integral confirms the sign, with no connection at all

Put the motion in the equatorial plane $\theta = \pi/2$, so that $\dot\theta = 0$ and $\sin^2\theta = 1$.
The $\theta$ geodesic equation keeps it there, since both of its terms carry either $\dot\theta$ or $\cos\theta$ and both vanish on that plane.
The metric has the two Killing vectors $\partial_t$ and $\partial_\phi$, so

$$E = \dot{t}, \qquad L = \rho^2\dot{\phi} = (r^2+\ell^2)\dot{\phi}$$

are constant along a geodesic, and the normalisation $g_{\mu\nu}\dot{x}^\mu\dot{x}^\nu = -\kappa$, with $\kappa = 1$ for a timelike geodesic and $0$ for a null one, reads

$$-\dot{t}^2 + \dot{r}^2 + \rho^2\dot{\phi}^2 = -\kappa,$$

so that

$$\dot{r}^2 = E^2 - \kappa - \frac{L^2}{r^2+\ell^2}.$$

Differentiating with respect to the affine parameter,

$$2\dot{r}\ddot{r} = L^2\cdot\frac{2r}{(r^2+\ell^2)^2}\dot{r}, \qquad \ddot{r} = \frac{L^2 r}{(r^2+\ell^2)^2} = r\left(\frac{L}{r^2+\ell^2}\right)^2 = r\dot{\phi}^2,$$

which is the boxed equation restricted to the equatorial plane.
The sign is positive, and no Christoffel symbol was used to get it.

This also says what the sign means physically.
The effective potential $L^2/(r^2+\ell^2)$ has its **maximum** at the throat $r = 0$, because $r^2+\ell^2$ has its minimum there, so a particle with angular momentum is pushed away from the throat rather than towards it.
The wrong sign asserted the opposite, that the throat attracts orbiting matter, which is exactly backwards for a geometry whose defining property is that its areal radius bottoms out at the throat.
Nothing is falling inwards in this spacetime: there is no mass, no horizon, and $\ddot{r} = 0$ for purely radial motion at every $r$, including at the throat.

### Two other spacetimes already had the right sign

**Minkowski.** At $\ell = 0$ this metric is flat space in spherical coordinates, whose radial geodesic is

$$\ddot{r} - r\dot{\theta}^2 - r\sin^2\theta\,\dot{\phi}^2 = 0,$$

the same equation in the other arrangement.
It is the familiar statement that a free particle passing the origin in a straight line has $r$ at a minimum and increasing, with the centrifugal term positive.
The Ellis-Bronnikov radial geodesic had the opposite sign to an equation its own $\ell = 0$ limit already carries.

**Morris-Thorne.** Set $\Phi = 0$ and $r(l) = \sqrt{l^2+\ell^2}$ in the proper radial chart of Morris-Thorne's wormhole, whose radial geodesic is

$$\ddot{l} + e^{2\Phi}\partial_l\Phi\,\dot{t}^2 - r\,\partial_l r\,\dot{\theta}^2 - r\,\partial_l r\sin^2\theta\,\dot{\phi}^2 = 0.$$

With $\Phi = 0$ the $\dot{t}^2$ term goes, and $\partial_l r = l/\sqrt{l^2+\ell^2}$ gives $r\,\partial_l r = l$, so the equation becomes

$$\ddot{l} - l\dot{\theta}^2 - l\sin^2\theta\,\dot{\phi}^2 = 0,$$

which is the boxed equation with $l$ for $r$.

---

## Step 10. Dimensions

The dimensional pass in `verify_metrics.py` reads a component of a rank four tensor with all indices down as carrying

$$\frac{1}{L^2}\cdot\frac{L^4}{[x^\mu][x^\nu][x^\rho][x^\sigma]},$$

where $[x^\mu]$ is the dimension of the chart coordinate, a length for $t$ because the chart coordinate is $ct$, a length for $r$ because it is a proper distance, and nothing at all for $\theta$ and $\phi$, which are angles.
The parameter $\ell$ is a length.

The twelve added components fall into two families by that rule, and both balance.

**Eight carry two angular indices and two length indices**, the eight Weyl components whose index pair is $\{t,\theta\}$, $\{t,\phi\}$, $\{r,\theta\}$ or $\{r,\phi\}$.
Such a component carries $L^4/(1\cdot L\cdot L\cdot 1)\cdot L^{-2} = 1$, so it has to be dimensionless, and $\ell^2/\left(3(r^2+\ell^2)\right)$ is a length squared over a length squared.
The four of them with a $\phi$ index carry the extra $\sin^2\theta$, which is dimensionless too.

**Four carry four angular indices**, the two Riemann components and the two Weyl components on the pair $\{\theta,\phi\}$.
Such a component carries $L^4/1\cdot L^{-2} = L^2$, and $\ell^2\sin^2\theta$ and $\tfrac{2}{3}\ell^2\sin^2\theta$ both carry $L^2$.

The pass reports no failure for this system, before or after.
The three failures it still reports across all the spacetimes are in the Einstein tensor with both indices up of the Lanczos-van Stockum dust in its cylindrical chart, and the corrections do not touch them.

---

## Step 11. Morris-Thorne, specialised, as an independent check

Morris-Thorne's spacetime, worked out on 18 September 2026, is the general static spherically symmetric wormhole with a free redshift function and a free shape function.
Its proper radial chart is

$$ds^2 = -e^{2\Phi(l)}c^2dt^2 + dl^2 + r(l)^2\left(d\theta^2 + \sin^2\theta\,d\phi^2\right),$$

and setting

$$\Phi = 0, \qquad r(l) = \sqrt{l^2+\ell^2}, \qquad l \to r$$

turns it into this spacetime exactly.
Setting $\Phi = 0$ is what is meant by calling Ellis-Bronnikov the wormhole with zero redshift function, and $r(l) = \sqrt{l^2+\ell^2}$ is the areal radius $\rho$ of Step 1 written along the proper distance.

That specialisation is a real check rather than a restatement, because the two spacetimes were derived independently and neither was written from the other.
It was carried out on every component of both spacetimes at once: each Morris-Thorne component was parsed, specialised and matched against the Ellis-Bronnikov component with the same indices, with $l$ renamed to $r$.

Against the Ellis-Bronnikov components as they stood, Morris-Thorne reproduced **99** values and produced **12** more that were missing.
Those 12 are exactly the 12 added, with exactly the values derived in Step 5 and Step 7.
Against the corrected components it reproduces all **111** and produces nothing that is not there.

Three of the specialisations are short enough to do by hand, one from each family.
Two derivatives are needed and both are elementary:

$$\partial_l r = \frac{l}{r}, \qquad \partial_l^2 r = \frac{1}{r} - \frac{l^2}{r^3} = \frac{r^2-l^2}{r^3} = \frac{\ell^2}{r^3},$$

so that $(\partial_l r)^2 = l^2/(l^2+\ell^2)$ and $r\,\partial_l^2 r = \ell^2/(l^2+\ell^2)$.

**$R_{\phi\theta\theta\phi}$.** In Morris-Thorne it is $-r^2\left(1 - (\partial_l r)^2\right)\sin^2\theta$.
Here $1 - (\partial_l r)^2 = \ell^2/(l^2+\ell^2)$ and $r^2 = l^2+\ell^2$, so the product is $-\ell^2\sin^2\theta$, which is Step 5.

**$C_{\theta t t\theta}$.** In Morris-Thorne it is

$$\frac{e^{2\Phi}\left(r^2\left((\partial_l\Phi)^2 + \partial_l^2\Phi\right) - r\,\partial_l\Phi\,\partial_l r + (\partial_l r)^2 - r\,\partial_l^2 r - 1\right)}{6}.$$

Every term carrying $\Phi$ drops, and $e^{2\Phi} = 1$, so what is left is

$$\frac{(\partial_l r)^2 - r\,\partial_l^2 r - 1}{6} = \frac{1}{6}\left(\frac{l^2}{l^2+\ell^2} - \frac{\ell^2}{l^2+\ell^2} - 1\right) = \frac{1}{6}\cdot\frac{-2\ell^2}{l^2+\ell^2} = -\frac{\ell^2}{3(l^2+\ell^2)},$$

which with $l \to r$ is Step 7.

**$C_{\phi\theta\theta\phi}$.** In Morris-Thorne it is $r^2\sin^2\theta$ times the same bracket over $3$, so the specialisation is

$$\frac{(l^2+\ell^2)\sin^2\theta}{3}\cdot\frac{-2\ell^2}{l^2+\ell^2} = -\frac{2}{3}\ell^2\sin^2\theta,$$

which is Step 7.

The check cuts both ways, and it is worth recording which way it cut here.
It could have found a disagreement in a value of either spacetime, and it found none in 99 of them, so the two spacetimes agree wherever both spoke.
Where they did not agree was where one of them was silent.

---

## Step 12. What the checker says

`verify_metrics.py` compares each stored component against sympy in the $x^0 = cT$ chart and separately checks that every term carries the dimension its left hand side fixes.

Before this change the checker reported 13 disagreements for Ellis-Bronnikov: the 2 Riemann omissions of Step 5, the 10 Weyl omissions of Step 7 and the 1 geodesic of Step 9.
After it, none.

```
$ python3 _tools/derivations/verify_metrics.py --system ellis_bronnikov/spherical

  ellis_bronnikov/spherical

29 metric files, 1 coordinate systems, 1 checked.

1 UNCHECKED:
  UNCHECKED ellis_bronnikov/spherical.inverse_metric_components: the entry does not publish one

Every published value balances dimensionally and agrees with sympy.
```

The one `UNCHECKED` line is unchanged and is not a disagreement.
It stands for the missing inverse metric, which the audit of 18 September 2026 already listed among its six `UNCHECKED` lines across all the spacetimes.
Supplying one is a separate job from these three corrections and was left alone.

The full sweep was run twice to show that nothing else moved, once against the tree as it stood and once against the corrected tree, and the two reports were diffed line by line.
Both runs were made against the same tree, at commit b7bb85a with only the three Ellis-Bronnikov corrections between them, so the only thing that differs between the two reports is this change.
That matters because several spacetimes were landing in parallel that day, and a comparison taken across two different trees would not have isolated anything.
The disagreements across all the spacetimes went from **72** to **59**, over 29 metrics and 39 coordinate systems both times.
The 13 lines that disappeared are exactly the 13 of Ellis-Bronnikov, and no line appeared that was not there before.

Before, the 72 stood at 51 in the interior Schwarzschild star, 13 in Ellis-Bronnikov, 4 in Gödel's universe and 4 in the Lanczos-van Stockum dust.
After, the 59 stand at 51, 4 and 4, with nothing left in Ellis-Bronnikov.
The 51 are the 50 nested radicals sympy will not combine, which the audit of 18 September 2026 records are not disagreements at all, plus the one inverse metric component of the interior Schwarzschild star that does not invert its own line element; the 8 are all in the Einstein tensors with both indices up of the two rotating dust spacetimes, Gödel's and Lanczos-van Stockum's.
None of them is in Ellis-Bronnikov and none of them was touched.

The dimensional pass reports 3 terms before and 3 after, the same three, all of them in the Einstein tensor with both indices up of the Lanczos-van Stockum dust in its cylindrical chart.

One difference between the two reports is not a disagreement and is worth naming so that the next reader is not puzzled by it.
The `UNCHECKED` count went from 171 to 160, and all 11 lines that went are Taub-NUT in its spherical chart, which timed out on the 120 second budget in the first run and finished inside it in the second.
The budget is wall clock and the machine was carrying a load average above 11 from other work during the first run.
That moves the comparison in the safe direction rather than the dangerous one: the second run checked eleven sets of components the first never compared at all, and still found nothing new to disagree with.

The totals of 72 and 59 across all the spacetimes belong to that one tree, and they are not what a sweep run later will report.
Krasnikov and Kerr-Newman landed during the same work, and both had until then been without any coordinate system at all, so the checkable coordinate systems went from 39 to 41 underneath these corrections.
What carries over from the pair of runs is the difference, which is 13 and is all of it in Ellis-Bronnikov.

---

## Step 13. The inverse metric

Every other coordinate system of every spacetime has its inverse metric, and a reader finds it beside the metric.
Ellis-Bronnikov had none, and that left the one slot that the checker reported `UNCHECKED` because there was nothing in it to compare.
We fill the slot rather than excuse it in the checker, because every spacetime gives its reader an inverse metric, and because the one thing a checker could say about an empty slot is that it is empty.

The metric of Step 2 is diagonal, so its inverse is diagonal too, with each component the reciprocal of the metric component in the same slot:

$$g^{\mu\nu} = \operatorname{diag}\left(\frac{1}{g_{tt}},\ \frac{1}{g_{rr}},\ \frac{1}{g_{\theta\theta}},\ \frac{1}{g_{\phi\phi}}\right).$$

In the chart $x^0 = ct$ that is

$$g^{tt} = -1, \qquad g^{rr} = 1, \qquad g^{\theta\theta} = \frac{1}{r^2+\ell^2}, \qquad g^{\phi\phi} = \frac{1}{\left(r^2+\ell^2\right)\sin^2\theta},$$

and every off diagonal component vanishes.
The product $g^{\mu\alpha}g_{\alpha\nu}$ is diagonal with entries $(-1)(-1)$, $1\cdot1$, $\rho^2/\rho^2$ and $\rho^2\sin^2\theta/(\rho^2\sin^2\theta)$, each of them $1$, so it is the identity wherever the metric is defined, which is everywhere off the poles $\sin\theta = 0$, the same coordinate singularity flat spherical coordinates have.

The time component follows the chart rule of Step 1 in the other direction from the metric's.
With the bare coordinate the time component of the line element is $-c^2$, whose reciprocal is $-1/c^2$, and a component with two upper time indices is multiplied by $c$ twice on the way into the chart, which gives $-1$.

The dimensions balance by the rule of Step 10 with the indices raised.
An upper index contributes $[x^\mu]/L$, so $g^{tt}$ and $g^{rr}$ are dimensionless, as $-1$ and $1$ are, and $g^{\theta\theta}$ and $g^{\phi\phi}$ carry $1/L^2$, as $1/(r^2+\ell^2)$ does.

The two angular values are written as in every other spacetime, the $\phi\phi$ one as Taub-NUT writes its own $1/\left((r^2+l^2)\sin^2\theta\right)$, with the sum in parentheses.
The inverse metric is set as a four by four matrix under the metric's own, one value per slot and zero in every slot left out, so no value is grouped with its negation.

With the inverse metric in place the system reports nothing `UNCHECKED`:

```
$ python3 _tools/derivations/verify_metrics.py --system ellis_bronnikov/spherical

  ellis_bronnikov/spherical

29 metric files, 1 coordinate systems, 1 checked.

Every published value balances dimensionally and agrees with sympy.
```

A wrong value in the inverse metric is caught, which was confirmed by giving $g^{\phi\phi}$ without its $\sin^2\theta$ for one run: the checker reported it as a disagreement against sympy's $1/\left((\ell^2+r^2)\sin^2\theta\right)$ and exited non-zero.
