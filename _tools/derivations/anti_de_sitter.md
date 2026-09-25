# The anti-de Sitter spacetime

We work anti-de Sitter space in two charts, the global static chart and the Poincaré patch.
We derive every component in order, from the line element down to the geodesic equations.
We leave nothing as an exercise and assert nothing we do not compute.

The companion script `_tools/derivations/verify_metrics.py` does the same work in sympy and compares the result component by component, so the algebra is checkable by hand and by machine independently.

Anti-de Sitter space is the third maximally symmetric solution, alongside Minkowski space and de Sitter space, and the only one of the three with a timelike boundary.
Almost all of its curvature is settled in one line, the Riemann tensor of Step 6 written as an antisymmetrised product of metrics; the Ricci tensor and scalar, the Einstein tensor, the Weyl tensor and the Kretschmann scalar of Steps 7 to 10 are that line contracted, traced and stripped of its traces.
Two things are not mere consequences of maximal symmetry: in Step 12 the global static chart has no horizon and a boundary reachable in finite coordinate time, and in Step 13 the Poincaré patch is built as a second chart on the same quadric.

---

## Step 1. Conventions

The signature is $(-,+,+,+)$ and the Christoffel symbols are those of the Levi-Civita connection,

$$\Gamma^\mu{}_{\nu\rho} = \tfrac{1}{2}g^{\mu\alpha}\left(\partial_\nu g_{\alpha\rho} + \partial_\rho g_{\alpha\nu} - \partial_\alpha g_{\nu\rho}\right),$$

with the lowered symbol $\Gamma_{\mu\nu\rho} = g_{\mu\alpha}\Gamma^\alpha{}_{\nu\rho}$.
The Riemann tensor is

$$R^\mu{}_{\nu\rho\sigma} = \partial_\rho \Gamma^\mu{}_{\nu\sigma} - \partial_\sigma \Gamma^\mu{}_{\nu\rho} + \Gamma^\mu{}_{\rho\lambda}\Gamma^\lambda{}_{\nu\sigma} - \Gamma^\mu{}_{\sigma\lambda}\Gamma^\lambda{}_{\nu\rho},$$

and the Riemann components follow it.

The Ricci tensor is the standard contraction,

$$R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu},$$

the contraction in use for every spacetime since 2026-09-18.
For anti-de Sitter space the choice of contraction is visible rather than invisible.
This spacetime has a nonzero Ricci tensor, and on this contraction it comes out as $R_{\mu\nu} = \Lambda g_{\mu\nu}$ with the sign of $\Lambda$ itself, so a negative cosmological constant gives a Ricci tensor whose $tt$ component is positive and whose spatial components are negative.
Contracting on the last index instead would give every component of Steps 7, 8 and 9 the opposite sign and would put a negative $\Lambda$ where a positive one belongs.

Factors of $G$ and $c$ are kept explicit.
The solution has one parameter, the anti-de Sitter radius

$$L, \qquad \Lambda = -\frac{3}{L^2},$$

which is a length, and no mass enters anywhere, so $G$ never appears in a component of the metric or the curvature: it appears only in the field equations of Step 8, where it multiplies a stress energy tensor that is zero.

The chart is the one used for every spacetime, whose zeroth coordinate is

$$x^0 = ct.$$

The index is written with the bare letter $t$, but the component written against it is a component of that chart.
Because the rescaling $t \to ct$ is linear with a constant coefficient, a chart component is the component taken with the bare coordinate multiplied by $c$ once per upper time index and divided by $c$ once per lower one, and the Christoffel symbols follow the same rule as the tensors.
We compute directly in the chart from Step 3 onward, so no conversion is needed at the end, and neither of the two metrics depends on $t$, so no factor of $c$ survives in any component of either.

The dots in the geodesic equations are velocities of that same chart, so $\dot{t}$ means $d(ct)/d\lambda$ and the equations are

$$\ddot{x}^\mu + \Gamma^\mu{}_{\nu\rho}\dot{x}^\nu\dot{x}^\rho = 0$$

with the same $\Gamma$.
This is the reading on which every term balances, as we check in Step 14.

One abbreviation runs through the global chart,

$$f(r) = 1 + \frac{r^2}{L^2} = \frac{L^2 + r^2}{L^2},$$

and it is only a shorthand.
Every component is written out in $r$ and $L$, because the checker reads each one symbol by symbol and has no way to be told what an abbreviation means.

---

## Step 2. The line element, and where the radius comes from

The static global form is

$$ds^2 = -\left(1 + \frac{r^2}{L^2}\right)c^2dt^2 + \frac{dr^2}{1 + r^2/L^2} + r^2\left(d\theta^2 + \sin^2\theta\,d\phi^2\right),$$

with $t \in (-\infty,\infty)$, $r \in [0,\infty)$ and the usual ranges for the two angles.

It is the induced metric on a quadric.
Take a flat five dimensional space of signature $(-,-,+,+,+)$,

$$ds_5^2 = -dU^2 - dV^2 + dX^2 + dY^2 + dZ^2,$$

and the hypersurface

$$-U^2 - V^2 + X^2 + Y^2 + Z^2 = -L^2,$$

It is connected, since the constraint reads $U^2 + V^2 = L^2 + X^2 + Y^2 + Z^2$ and so keeps the pair $(U,V)$ outside a circle of radius $L$, and its topology is that circle times the three remaining directions.
What makes it anti-de Sitter rather than de Sitter is that the two directions entering with a minus sign are both timelike.
Parametrise it by

$$U = \sqrt{L^2+r^2}\,\cos\frac{ct}{L}, \qquad V = \sqrt{L^2+r^2}\,\sin\frac{ct}{L}, \qquad (X,Y,Z) = r\,\hat{n}(\theta,\phi),$$

with $\hat{n}$ the unit radial vector.
The constraint holds identically, since $-(L^2+r^2) + r^2 = -L^2$, and

$$dU^2 + dV^2 = \frac{r^2}{L^2+r^2}dr^2 + \frac{L^2+r^2}{L^2}c^2dt^2,$$

so that

$$ds^2 = -\frac{L^2+r^2}{L^2}c^2dt^2 + \left(1 - \frac{r^2}{L^2+r^2}\right)dr^2 + r^2d\Omega^2 = -f\,c^2dt^2 + \frac{dr^2}{f} + r^2d\Omega^2,$$

which is the line element of the global static chart.

Two things follow immediately from the embedding and are used later.
The isometry group of the quadric is the group preserving the flat five dimensional form, $SO(3,2)$, which has ten generators, and ten is the maximum a four dimensional spacetime can carry, so this spacetime is maximally symmetric.
And $(U,V)$ traces a circle as $ct$ runs over $[0, 2\pi L)$, so on the quadric itself $t$ is periodic with period $2\pi L/c$, and since $\partial_t$ is timelike everywhere, every one of those circles is a closed timelike curve.
Taking $t$ over the whole real line, as the domain of the coordinates does, is passing to the universal cover, which unwraps that circle and removes every closed timelike curve at once.

De Sitter space is the same construction with one timelike direction fewer, on $-U^2 + X^2 + Y^2 + Z^2 + W^2 = +L^2$, and its static form carries $1 - r^2/L^2$ where this one carries $1 + r^2/L^2$.
The whole of the difference between the two spacetimes is that sign, and its consequences are paid out in Step 12.

---

## Step 3. The metric matrix, its determinant and its inverse

In the chart $x^0 = ct$ the line element gives a diagonal metric,

$$g_{\mu\nu} = \mathrm{diag}\left(-f,\ f^{-1},\ r^2,\ r^2\sin^2\theta\right), \qquad f = \frac{L^2+r^2}{L^2}.$$

Every entry is dimensionless in the first two slots and carries $L^2$ in the angular ones, which is what the chart asks of a metric whose zeroth coordinate is a length and whose angles are pure numbers.

The determinant is

$$\det g = (-f)\left(f^{-1}\right)\left(r^2\right)\left(r^2\sin^2\theta\right) = -r^4\sin^2\theta,$$

so $\sqrt{-g} = r^2\sin\theta$, exactly as in flat space in spherical coordinates: the two factors of $f$ cancel, which is the same accident that makes Schwarzschild's volume element flat.
The metric is nowhere degenerate for $r > 0$, and $r = 0$ is the ordinary coordinate axis of spherical coordinates rather than anything geometric.

The inverse is the elementwise reciprocal,

$$g^{\mu\nu} = \mathrm{diag}\left(-f^{-1},\ f,\ r^{-2},\ r^{-2}\sin^{-2}\theta\right),$$

which in components is $g^{tt} = -\left(1 + r^2/L^2\right)^{-1}$, $g^{rr} = 1 + r^2/L^2$, $g^{\theta\theta} = 1/r^2$ and $g^{\phi\phi} = 1/r^2\sin^2\theta$.

---

## Step 4. The Christoffel symbols with an upper index

The metric is static and diagonal and depends on $r$ and $\theta$ alone, so the surviving symbols are those of any metric of the form $-f\,c^2dt^2 + f^{-1}dr^2 + r^2d\Omega^2$, with

$$f' = \frac{2r}{L^2}.$$

The time symbols come from $\partial_r g_{tt}$ alone,

$$\Gamma^t{}_{tr} = \Gamma^t{}_{rt} = \frac{f'}{2f} = \frac{r}{L^2+r^2}, \qquad \Gamma^r{}_{tt} = \frac{ff'}{2} = \frac{r\left(L^2+r^2\right)}{L^4}.$$

The first of these is the one component where a factor of $c$ could have survived and does not: it carries one upper time index and one lower time index, so the chart weight is $c^0$.
The second carries two lower time indices and no upper one, so its chart value is $c^{-2}$ times the value computed with the bare $t$, and the $c^2$ sitting in the bare $g_{tt}$ is exactly what that cancels.

The radial symbols are

$$\Gamma^r{}_{rr} = -\frac{f'}{2f} = -\frac{r}{L^2+r^2}, \qquad \Gamma^r{}_{\theta\theta} = -rf = -\frac{r\left(L^2+r^2\right)}{L^2}, \qquad \Gamma^r{}_{\phi\phi} = -rf\sin^2\theta,$$

and the angular ones are those of the round sphere of radius $r$,

$$\Gamma^\theta{}_{r\theta} = \Gamma^\theta{}_{\theta r} = \Gamma^\phi{}_{r\phi} = \Gamma^\phi{}_{\phi r} = \frac{1}{r}, \qquad \Gamma^\theta{}_{\phi\phi} = -\cos\theta\sin\theta, \qquad \Gamma^\phi{}_{\theta\phi} = \Gamma^\phi{}_{\phi\theta} = \cot\theta.$$

That is thirteen nonzero symbols counting both orderings of each symmetric pair, and they are all the Christoffel symbols with the first index up.
Setting $L \to \infty$ sends $f \to 1$ and leaves only the four angular symbols, which are flat space in spherical coordinates; that limit is the first check to make on any later expression.

---

## Step 5. The lowered Christoffel symbols

Lowering the first index with the diagonal metric is multiplication by one entry,

$$\Gamma_{\mu\nu\rho} = g_{\mu\mu}\Gamma^\mu{}_{\nu\rho} \quad (\text{no sum}),$$

and the factors of $f$ cancel in most slots.

$$\Gamma_{ttr} = \Gamma_{trt} = (-f)\frac{r}{L^2+r^2} = -\frac{r}{L^2}, \qquad \Gamma_{rtt} = f^{-1}\frac{r\left(L^2+r^2\right)}{L^4} = \frac{r}{L^2},$$

$$\Gamma_{rrr} = f^{-1}\left(-\frac{r}{L^2+r^2}\right) = -\frac{L^2r}{\left(L^2+r^2\right)^2}, \qquad \Gamma_{r\theta\theta} = -r, \qquad \Gamma_{r\phi\phi} = -r\sin^2\theta,$$

$$\Gamma_{\theta r\theta} = \Gamma_{\theta\theta r} = r, \qquad \Gamma_{\theta\phi\phi} = -r^2\cos\theta\sin\theta, \qquad \Gamma_{\phi r\phi} = \Gamma_{\phi\phi r} = r\sin^2\theta, \qquad \Gamma_{\phi\theta\phi} = \Gamma_{\phi\phi\theta} = r^2\cos\theta\sin\theta.$$

Thirteen again, since lowering an index cannot create or destroy a zero in a diagonal metric, and these are all the Christoffel symbols with every index lowered.
In the pair $\Gamma_{ttr} = -r/L^2$ and $\Gamma_{rtt} = +r/L^2$ the magnitudes are equal and the signs are opposite, which is the statement that the only thing bending a static worldline here is the growth of $g_{tt}$ with radius.

---

## Step 6. The Riemann tensor, in one line

A maximally symmetric spacetime has a Riemann tensor built from the metric alone, because there is no other tensor available that has its symmetries and is invariant under the whole isometry group:

$$R_{\mu\nu\rho\sigma} = K_0\left(g_{\mu\rho}g_{\nu\sigma} - g_{\mu\sigma}g_{\nu\rho}\right),$$

with $K_0$ a constant, the sectional curvature of every two plane at every point.
For anti-de Sitter space

$$K_0 = -\frac{1}{L^2} = \frac{\Lambda}{3},$$

negative, which is the whole content of the word anti in its name.

Rather than assume that, here is the component that fixes it, computed from Step 4 by the definition in Step 1:

$$R^t{}_{rtr} = \partial_t\Gamma^t{}_{rr} - \partial_r\Gamma^t{}_{rt} + \Gamma^t{}_{t\lambda}\Gamma^\lambda{}_{rr} - \Gamma^t{}_{r\lambda}\Gamma^\lambda{}_{rt}.$$

The first term is zero because nothing depends on $t$, and $\Gamma^t{}_{rr}$ is zero besides.
The rest are

$$-\partial_r\frac{r}{L^2+r^2} = -\frac{L^2-r^2}{\left(L^2+r^2\right)^2}, \qquad \Gamma^t{}_{tr}\Gamma^r{}_{rr} = -\frac{r^2}{\left(L^2+r^2\right)^2}, \qquad -\Gamma^t{}_{rt}\Gamma^t{}_{rt} = -\frac{r^2}{\left(L^2+r^2\right)^2},$$

and their sum is

$$R^t{}_{rtr} = \frac{-L^2 + r^2 - r^2 - r^2}{\left(L^2+r^2\right)^2} = -\frac{L^2+r^2}{\left(L^2+r^2\right)^2} = -\frac{1}{L^2+r^2}.$$

Against the maximally symmetric form this is $K_0 g_{rr} = K_0 L^2/(L^2+r^2)$, so $K_0 = -1/L^2$, with no $r$ left in it.
A second component, in a plane with no time in it at all, confirms the same constant:

$$R^\theta{}_{\phi\theta\phi} = \partial_\theta\Gamma^\theta{}_{\phi\phi} + \Gamma^\theta{}_{\theta r}\Gamma^r{}_{\phi\phi} - \Gamma^\theta{}_{\phi\phi}\Gamma^\phi{}_{\phi\theta} = \left(\sin^2\theta - \cos^2\theta\right) - \frac{\left(L^2+r^2\right)\sin^2\theta}{L^2} + \cos^2\theta = -\frac{r^2\sin^2\theta}{L^2},$$

which is $K_0 g_{\phi\phi}$ with the same $K_0$.
The two planes are as unlike as this spacetime has to offer, one containing the time direction and one not, and they agree, which is maximal symmetry showing itself rather than being assumed.

With $K_0$ fixed, every component follows by index gymnastics.
Raising the first index of the maximally symmetric form gives

$$R^\mu{}_{\nu\rho\sigma} = -\frac{1}{L^2}\left(\delta^\mu{}_\rho\, g_{\nu\sigma} - \delta^\mu{}_\sigma\, g_{\nu\rho}\right),$$

so in a diagonal metric the nonzero components are exactly

$$R^\mu{}_{\nu\mu\nu} = -\frac{g_{\nu\nu}}{L^2}, \qquad R^\mu{}_{\nu\nu\mu} = +\frac{g_{\nu\nu}}{L^2} \qquad (\mu \neq \nu, \text{ no sum}),$$

which is two components for each of the twelve ordered pairs of distinct indices, twenty four in all.
The fully lowered tensor is nonzero on the same planes,

$$R_{\mu\nu\mu\nu} = R_{\nu\mu\nu\mu} = -\frac{g_{\mu\mu}g_{\nu\nu}}{L^2}, \qquad R_{\mu\nu\nu\mu} = R_{\nu\mu\mu\nu} = +\frac{g_{\mu\mu}g_{\nu\nu}}{L^2},$$

four components for each of the six unordered pairs, twenty four again.
Both forms are given in full, and the checker requires every other component to vanish, which holds the count to exactly these.

The six planes of the global chart carry

$$R_{trtr} = \frac{1}{L^2}, \qquad R_{t\theta t\theta} = \frac{r^2\left(L^2+r^2\right)}{L^4}, \qquad R_{t\phi t\phi} = \frac{r^2\left(L^2+r^2\right)\sin^2\theta}{L^4},$$

$$R_{r\theta r\theta} = -\frac{r^2}{L^2+r^2}, \qquad R_{r\phi r\phi} = -\frac{r^2\sin^2\theta}{L^2+r^2}, \qquad R_{\theta\phi\theta\phi} = -\frac{r^4\sin^2\theta}{L^2}.$$

The three timelike planes come out positive and the three spacelike ones negative, which is the signature at work rather than a change of curvature: every sectional curvature, formed with the right normalisation for the plane's own signature, is the same $-1/L^2$.

---

## Step 7. The Ricci tensor and the Ricci scalar

Contract the first index of the lowered form with the third, which is the standard contraction of Step 1:

$$R_{\mu\nu} = g^{\alpha\beta}R_{\alpha\mu\beta\nu} = -\frac{1}{L^2}g^{\alpha\beta}\left(g_{\alpha\beta}g_{\mu\nu} - g_{\alpha\nu}g_{\mu\beta}\right) = -\frac{1}{L^2}\left(4g_{\mu\nu} - g_{\mu\nu}\right) = -\frac{3}{L^2}g_{\mu\nu}.$$

So

$$R_{\mu\nu} = \Lambda g_{\mu\nu}, \qquad \Lambda = -\frac{3}{L^2},$$

which is the defining property of an Einstein space, here with a negative constant.
In the global chart that is

$$R_{tt} = \frac{3\left(L^2+r^2\right)}{L^4}, \qquad R_{rr} = -\frac{3}{L^2+r^2}, \qquad R_{\theta\theta} = -\frac{3r^2}{L^2}, \qquad R_{\phi\phi} = -\frac{3r^2\sin^2\theta}{L^2}.$$

The mixed and fully raised forms are the same statement with indices moved,

$$R^\mu{}_\nu = \Lambda\,\delta^\mu{}_\nu = -\frac{3}{L^2}\delta^\mu{}_\nu, \qquad R^{\mu\nu} = \Lambda\,g^{\mu\nu},$$

and the mixed form is the plainest, because all four of its components are the same negative constant and no coordinate appears in any of them.
That is the clearest form of the statement that this spacetime is the same at every point and in every direction.

The scalar is the trace,

$$R = g^{\mu\nu}R_{\mu\nu} = 4\Lambda = -\frac{12}{L^2},$$

constant and negative.
De Sitter space, run through the same four lines with $K_0 = +1/L^2$, gives $R = +12/L^2$.

---

## Step 8. The Einstein tensor, and why it is not zero

$$G_{\mu\nu} = R_{\mu\nu} - \tfrac{1}{2}Rg_{\mu\nu} = \Lambda g_{\mu\nu} - 2\Lambda g_{\mu\nu} = -\Lambda g_{\mu\nu} = \frac{3}{L^2}g_{\mu\nu}.$$

A vacuum with a nonzero Einstein tensor looks like a contradiction, and is not.
The field equations with a cosmological constant are

$$G_{\mu\nu} + \Lambda g_{\mu\nu} = \frac{8\pi G}{c^4}T_{\mu\nu},$$

and by the Einstein tensor just computed the left hand side vanishes identically, so $T_{\mu\nu} = 0$: there is no matter anywhere in this spacetime.
The Einstein tensor is nonzero only because the cosmological term has been kept on the left, where it is geometry, rather than moved to the right and called a vacuum energy.
Moving it is a relabelling and changes no physics, but it is worth naming what it would say here: a vacuum energy density of $\Lambda c^4/8\pi G$, which for $\Lambda < 0$ is negative, and negative vacuum energy is what makes this spacetime a box instead of an expansion.

Component by component the Einstein tensor is minus the Ricci tensor,

$$G_{tt} = -\frac{3\left(L^2+r^2\right)}{L^4}, \qquad G_{rr} = \frac{3}{L^2+r^2}, \qquad G_{\theta\theta} = \frac{3r^2}{L^2}, \qquad G_{\phi\phi} = \frac{3r^2\sin^2\theta}{L^2}, \qquad G^\mu{}_\nu = \frac{3}{L^2}\delta^\mu{}_\nu,$$

and that relation, $G_{\mu\nu} = -R_{\mu\nu}$, holds in any four dimensional Einstein space and is the quickest check on the two tensors.

---

## Step 9. The Weyl tensor, computed rather than copied

The Weyl tensor is the Riemann tensor with all of its traces removed,

$$C_{\mu\nu\rho\sigma} = R_{\mu\nu\rho\sigma} - \frac{1}{n-2}\left(g_{\mu\rho}R_{\sigma\nu} - g_{\mu\sigma}R_{\rho\nu} - g_{\nu\rho}R_{\sigma\mu} + g_{\nu\sigma}R_{\rho\mu}\right) + \frac{R}{(n-1)(n-2)}\left(g_{\mu\rho}g_{\sigma\nu} - g_{\mu\sigma}g_{\rho\nu}\right),$$

with $n = 4$ here, so the two coefficients are $1/2$ and $1/6$.

Substituting Steps 6 and 7 into it, with $P_{\mu\nu\rho\sigma} = g_{\mu\rho}g_{\nu\sigma} - g_{\mu\sigma}g_{\nu\rho}$ for the antisymmetrised product:

$$R_{\mu\nu\rho\sigma} = \Lambda P_{\mu\nu\rho\sigma}/3, \qquad R_{\mu\nu} = \Lambda g_{\mu\nu}, \qquad R = 4\Lambda.$$

The middle bracket becomes $\Lambda\left(g_{\mu\rho}g_{\sigma\nu} - g_{\mu\sigma}g_{\rho\nu} - g_{\nu\rho}g_{\sigma\mu} + g_{\nu\sigma}g_{\rho\mu}\right) = 2\Lambda P_{\mu\nu\rho\sigma}$, and the last is $4\Lambda P_{\mu\nu\rho\sigma}/6$, so

$$C_{\mu\nu\rho\sigma} = \frac{\Lambda}{3}P - \frac{1}{2}\cdot 2\Lambda P + \frac{2\Lambda}{3}P = \left(\frac{\Lambda}{3} - \Lambda + \frac{2\Lambda}{3}\right)P = 0,$$

identically, in every slot, for either sign of $\Lambda$.

The reason is worth stating plainly, because it is not the reason a vacuum's Weyl tensor equals its Riemann tensor.
There the Ricci tensor is zero and so the subtraction removes nothing.
Here the Ricci tensor is emphatically not zero; what happens instead is that the Riemann tensor of a maximally symmetric space is built out of its own traces and nothing else, so removing the traces removes all of it.
The general distinction is worked in `_tools/derivations/weyl.md`, and two other spacetimes have a vanishing Weyl tensor, FRW and the interior Schwarzschild solution, which are conformally flat without being maximally symmetric.
Conformal flatness is what a vanishing Weyl tensor means in four dimensions: it is the integrability condition for the metric to be a position dependent multiple of a flat one in some chart.
That chart is the Poincaré patch of Step 13.

The Weyl tensor therefore vanishes both with the first index up and with every index lowered, and the checker turns each of the two into twenty four assertions that a component is zero.

---

## Step 10. The Kretschmann scalar

Raise all four indices of the maximally symmetric form and contract with itself:

$$K = R_{\mu\nu\rho\sigma}R^{\mu\nu\rho\sigma} = K_0^2\left(g_{\mu\rho}g_{\nu\sigma} - g_{\mu\sigma}g_{\nu\rho}\right)\left(g^{\mu\rho}g^{\nu\sigma} - g^{\mu\sigma}g^{\nu\rho}\right).$$

The four products are $n^2$, $-n$, $-n$ and $n^2$, so

$$K = 2n(n-1)K_0^2 = 24K_0^2 = \frac{24}{L^4},$$

and in terms of the cosmological constant

$$K = \frac{8\Lambda^2}{3}.$$

It is a constant, positive, and blind to the sign of $\Lambda$: de Sitter space of the same radius has the same Kretschmann scalar.
That constancy is the statement that there is no curvature singularity anywhere in this spacetime, and no distinguished place in it at all.
Any place a chart appears to single out, such as the Poincaré horizon of Step 13, is a fact about the chart, and $K$ is the quickest way to see it.

---

## Step 11. The geodesic equations

With the symbols of Step 4 in

$$\ddot{x}^\mu + \Gamma^\mu{}_{\nu\rho}\dot{x}^\nu\dot{x}^\rho = 0,$$

and remembering that each off diagonal pair appears twice in the sum,

$$\ddot{t} + \frac{2r}{L^2+r^2}\dot{t}\dot{r} = 0,$$

$$\ddot{r} + \frac{r\left(L^2+r^2\right)}{L^4}\dot{t}^2 - \frac{r}{L^2+r^2}\dot{r}^2 - \frac{r\left(L^2+r^2\right)}{L^2}\left(\dot{\theta}^2 + \sin^2\theta\,\dot{\phi}^2\right) = 0,$$

$$\ddot{\theta} + \frac{2}{r}\dot{r}\dot{\theta} - \cos\theta\sin\theta\,\dot{\phi}^2 = 0,$$

$$\ddot{\phi} + \frac{2}{r}\dot{r}\dot{\phi} + 2\cot\theta\,\dot{\theta}\dot{\phi} = 0.$$

The first equation is $\frac{d}{d\lambda}\left(f\dot{t}\right) = 0$ in disguise, which is the conserved energy of the static Killing vector, and the last two are those of any spherically symmetric metric, so motion in a plane stays in that plane.

The sign of the $\dot{t}^2$ term in the radial equation is the one that matters.
It is positive, so the coordinate acceleration it produces is inward at every radius, for every worldline with $\dot{t} \neq 0$, and it grows without bound as $r$ grows.
Schwarzschild's corresponding term falls off as $1/r^2$ and de Sitter's is negative.
Anti-de Sitter space pulls everything back toward the origin harder the further out it goes, which is Step 12 in one line.

---

## Step 12. No horizon, and a boundary that light can reach

Two facts about this chart are the mirror of de Sitter's.

The first is that

$$f = 1 + \frac{r^2}{L^2} \geq 1$$

never vanishes and never changes sign.
So $g_{tt} < 0$ everywhere, $\partial_t$ is timelike at every point, and there is no Killing horizon anywhere in the spacetime.
The static chart is global: unlike Schwarzschild's, which stops at $r_s$, and unlike de Sitter's, whose $1 - r^2/L^2$ vanishes at $r = L$ and hides everything beyond a cosmological horizon, this one covers the whole spacetime, and observers at rest at every radius share one and the same notion of time.

The second is that the edge of the spacetime is close, measured in that time.
A radial null ray has $ds^2 = 0$ with $d\theta = d\phi = 0$, so

$$c\,dt = \frac{dr}{f} = \frac{dr}{1 + r^2/L^2},$$

and the coordinate time to go from the origin out to infinity is

$$c\,\Delta t = \int_0^\infty\frac{dr}{1 + r^2/L^2} = L\arctan\frac{r}{L}\bigg|_0^\infty = \frac{\pi L}{2},$$

which is finite.
Proper radial distance out to infinity is not finite,

$$\int_0^\infty\frac{dr}{\sqrt{f}} = \infty,$$

so the boundary is infinitely far away and reached in finite time, which is only possible because the coordinate speed of light $dr/dt = cf$ grows without bound.
Light leaves the origin, arrives at spatial infinity at $ct = \pi L/2$, and can be reflected back to arrive home at $ct = \pi L$.

The consequence is the one the history of anti-de Sitter space is built on.
The conformal boundary at $r \to \infty$ is timelike, not null as it is in an asymptotically flat spacetime, so information can enter through it in finite time.
No spacelike surface is a Cauchy surface: the future of any such surface depends on what comes in from the edge, and anti-de Sitter space is not globally hyperbolic.
Evolving a field in it is not an initial value problem alone; it is an initial value problem plus a boundary condition at infinity.
That is a defect if one wants a self contained cosmology, and it is the entire point if one wants a boundary for a field theory to live on.

Massive bodies, by contrast, never get out.
For a radial timelike geodesic, write $E = f\,d(ct)/d\tau$ for the conserved quantity of the first geodesic equation and normalise $g_{\mu\nu}\dot{x}^\mu\dot{x}^\nu = -c^2$ with $\tau$ the proper time.
Then

$$-\frac{E^2}{f} + \frac{1}{f}\left(\frac{dr}{d\tau}\right)^2 = -c^2 \quad\Longrightarrow\quad \left(\frac{dr}{d\tau}\right)^2 = E^2 - c^2f = \frac{c^2}{L^2}\left(A^2 - r^2\right), \qquad A^2 = \frac{L^2\left(E^2-c^2\right)}{c^2}.$$

The motion is a harmonic oscillation between $r = 0$ and the turning point $r = A$, which is finite for every finite $E$, so no timelike geodesic reaches the boundary.
Dividing the two rates gives the coordinate time as an integral over the radius,

$$\frac{d(ct)}{dr} = \frac{E}{f}\cdot\frac{L}{c\sqrt{A^2-r^2}} = \frac{EL^3}{c\left(L^2+r^2\right)\sqrt{A^2-r^2}},$$

so the round trip out to the turning point and back takes

$$c\,\Delta t = \frac{2EL^3}{c}\int_0^A\frac{dr}{\left(L^2+r^2\right)\sqrt{A^2-r^2}} = \frac{2EL^3}{c}\cdot\frac{\pi}{2L\sqrt{L^2+A^2}} = \frac{\pi EL^2}{c\sqrt{L^2+A^2}},$$

and with $\sqrt{L^2+A^2} = LE/c$ it collapses to

$$c\,\Delta t = \pi L,$$

independent of $E$ altogether.
Every radial timelike geodesic, however energetically it is thrown, returns to the origin after the same coordinate time $ct = \pi L$, exactly twice the time light takes to reach the boundary.
Anti-de Sitter space focuses its own geodesics, and that is the sense in which the negative curvature acts as a confining box.

---

## Step 13. The Poincaré patch

The second chart is the one every holography paper is written in,

$$ds^2 = \frac{L^2}{z^2}\left(-c^2dt^2 + dx^2 + dy^2 + dz^2\right), \qquad z \in (0,\infty),$$

a position dependent multiple of the Minkowski metric, which is possible at all only because the Weyl tensor of Step 9 vanishes.

It is a chart on the same quadric.
Take

$$U = \frac{z^2 + L^2 + x^2 + y^2 - c^2t^2}{2z}, \qquad W = \frac{z^2 - L^2 + x^2 + y^2 - c^2t^2}{2z}, \qquad V = \frac{Lct}{z}, \qquad X = \frac{Lx}{z}, \qquad Y = \frac{Ly}{z},$$

and the constraint $-U^2 - V^2 + X^2 + Y^2 + W^2 = -L^2$ holds identically while the induced metric is the line element of the Poincaré patch.
Since $U - W = L^2/z$, the chart covers only the half of the quadric where $U > W$, and the null surface $U = W$, which is $z \to \infty$, bounds it.
That surface is the Poincaré horizon.
Nothing happens there: the Kretschmann scalar of Step 10 is the same constant $24/L^4$ on it as everywhere else, and the global chart of Step 2 runs straight through it without noticing.
It is a horizon of the chart, in the same way that the Rindler horizon is a horizon of the Rindler chart on flat space.

The other end is the interesting one.
As $z \to 0$ the conformal factor $L^2/z^2$ diverges, and what it multiplies is Minkowski space in $(t,x,y)$.
Stripping off the factor leaves a finite metric on a three dimensional flat spacetime, and that is the conformal boundary: the same timelike boundary the global chart reaches at $r \to \infty$, now presented as a flat spacetime one dimension down.
That presentation is what makes the chart useful, since the field theory a holographic duality puts there is a conformal field theory on ordinary Minkowski space.

The curvature is quickest from the conformal form.
Write $g_{\mu\nu} = e^{2\omega}\eta_{\mu\nu}$ with $\omega = \ln(L/z)$, so that $\partial_z\omega = -1/z$ and every other derivative of $\omega$ vanishes.
The connection of a conformally rescaled flat metric is

$$\Gamma^\mu{}_{\nu\rho} = \delta^\mu{}_\nu\partial_\rho\omega + \delta^\mu{}_\rho\partial_\nu\omega - \eta_{\nu\rho}\eta^{\mu\sigma}\partial_\sigma\omega,$$

which with only $\partial_z\omega$ surviving gives

$$\Gamma^t{}_{tz} = \Gamma^x{}_{xz} = \Gamma^y{}_{yz} = \Gamma^z{}_{zz} = -\frac{1}{z}, \qquad \Gamma^z{}_{tt} = -\frac{1}{z}, \qquad \Gamma^z{}_{xx} = \Gamma^z{}_{yy} = +\frac{1}{z},$$

together with the symmetric partners of the first three, ten nonzero symbols in all.
The sign flip between $\Gamma^z{}_{tt}$ and $\Gamma^z{}_{xx}$ is the $\eta_{\nu\rho}$ in the last term, and it is the only place the signature enters.
Lowering the first index multiplies each by $\pm L^2/z^2$ and gives the ten Christoffel symbols with every index lowered, all of them $\pm L^2/z^3$.

The curvature is Step 6 again, because the spacetime is the same one:

$$R_{\mu\nu\rho\sigma} = -\frac{1}{L^2}\left(g_{\mu\rho}g_{\nu\sigma} - g_{\mu\sigma}g_{\nu\rho}\right), \qquad R_{\mu\nu} = -\frac{3}{L^2}g_{\mu\nu}, \qquad R = -\frac{12}{L^2}, \qquad K = \frac{24}{L^4}, \qquad C_{\mu\nu\rho\sigma} = 0,$$

so with $g_{tt} = -L^2/z^2$ and $g_{xx} = g_{yy} = g_{zz} = L^2/z^2$ the components are

$$R^\mu{}_{t\mu t} = \frac{1}{z^2}, \qquad R^\mu{}_{i\mu i} = -\frac{1}{z^2}, \qquad R_{titi} = \frac{L^2}{z^4}, \qquad R_{ijij} = -\frac{L^2}{z^4},$$

$$R_{tt} = \frac{3}{z^2}, \qquad R_{ii} = -\frac{3}{z^2}, \qquad R^\mu{}_\nu = -\frac{3}{L^2}\delta^\mu{}_\nu, \qquad G^\mu{}_\nu = \frac{3}{L^2}\delta^\mu{}_\nu,$$

with $i$ and $j$ distinct spatial indices and no sum anywhere.
Swapping the last two indices of a Riemann component flips its sign, which doubles each of these into twenty four nonzero components, with the first index up and with every index lowered alike, the same count as the global chart since it depends on the dimension and not on the chart.

The geodesic equations are

$$\ddot{t} - \frac{2}{z}\dot{t}\dot{z} = 0, \qquad \ddot{x} - \frac{2}{z}\dot{x}\dot{z} = 0, \qquad \ddot{y} - \frac{2}{z}\dot{y}\dot{z} = 0, \qquad \ddot{z} - \frac{1}{z}\left(\dot{t}^2 - \dot{x}^2 - \dot{y}^2 + \dot{z}^2\right) = 0,$$

and the first three integrate at once to $\dot{t}, \dot{x}, \dot{y} \propto z^2$, which are the three conserved momenta of the boundary directions.
A ray of light sent straight out at the boundary has $\dot{x} = \dot{y} = 0$ and $\dot{z} = \pm\dot{t}$, so $dz/d(ct) = \pm 1$, and it covers the whole depth of the patch at unit coordinate speed: it reaches $z = 0$ from any depth in finite coordinate time, which is the statement of Step 12 in the global chart.

---

## Step 14. Every expression is dimensionally consistent

The chart coordinates of the global system are $(ct, r, \theta, \phi)$, of dimensions $L$, $L$, $1$, $1$, and those of the Poincaré patch are $(ct, x, y, z)$, all four lengths.
The one parameter $L$ is a length in both.
A metric component carries

$$[g_{\mu\nu}] = \frac{L^2}{[x^\mu][x^\nu]},$$

a Christoffel symbol carries $L^{-1}$ corrected by $[x^\mu]/L$ per upper index and $L/[x^\mu]$ per lower one, a Riemann, Ricci, Einstein or Weyl component carries $L^{-2}$ corrected the same way, and $K$ carries $L^{-4}$.

Four samples, one per rank and one per chart.

$g_{\phi\phi} = r^2\sin^2\theta$ must carry $L^2/(1\cdot 1) = L^2$, and it does, while $g_{tt} = -(1 + r^2/L^2)$ must carry $L^2/(L\cdot L) = 1$, and the ratio $r^2/L^2$ is what makes it dimensionless.
This is where the anti-de Sitter radius earns its dimension: if $L$ were not a length there would be no way to add $1$ to $r^2/L^2$ at all.

$\Gamma^r{}_{tt} = r\left(L^2+r^2\right)/L^4$ carries $L^{-1}$ from the field, times $L/L = 1$ for the upper $r$, times $L/L = 1$ for each lower $t$, so $L^{-1}$ in all; and $L\cdot L^2/L^4$ is $L^{-1}$.

$R_{\theta\phi\theta\phi} = -r^4\sin^2\theta/L^2$ carries $L^{-2}$ from the field times $L/1 = L$ for each of its four lower angular indices, so $L^2$; and $L^4/L^2$ is $L^2$.

$R^{tt} = 3z^2/L^4$ in the Poincaré patch carries $L^{-2}$ times $L/L = 1$ for each upper index, so $L^{-2}$; and $L^2/L^4$ is $L^{-2}$.

The geodesic equations are measured against their own second derivatives.
In the radial equation of Step 11 $\ddot{r}$ carries $L/\lambda^2$, and the term $r\left(L^2+r^2\right)\dot{t}^2/L^4$ carries

$$\frac{L\cdot L^2}{L^4}\cdot\left(\frac{L}{\lambda}\right)^2 = \frac{L}{\lambda^2},$$

which matches, and matches only because $\dot{t}$ is $d(ct)/d\lambda$ and not $dt/d\lambda$.
On the other reading that term would come out one factor of $c$ away from $\ddot{r}$, and the dimensional pass in `verify_metrics.py` would name it.
It does not, on any expression of either chart.

---

## Step 15. The components, and what the checker needs

There are two coordinate systems, the global static chart and the Poincaré patch.

The global static chart has the line element; four $g_{\mu\nu}$ and four $g^{\mu\nu}$; thirteen $\Gamma^\mu{}_{\nu\rho}$ and thirteen $\Gamma_{\mu\nu\rho}$; twenty four $R^\mu{}_{\nu\rho\sigma}$ and twenty four $R_{\mu\nu\rho\sigma}$; four each of $R_{\mu\nu}$, $R^\mu{}_\nu$, $R^{\mu\nu}$, $G_{\mu\nu}$, $G^\mu{}_\nu$ and $G^{\mu\nu}$; $R = -12/L^2$; $K = 24/L^4$; a vanishing $C^\mu{}_{\nu\rho\sigma}$ and $C_{\mu\nu\rho\sigma}$; and four geodesic equations.

The Poincaré patch has the same shape with ten Christoffel symbols in each index position instead of thirteen, its transverse directions being flat coordinates rather than angles.

For `verify_metrics.py` to read them, the two systems need two lines in `DIMENSIONS`,

    ("anti_de_sitter", "static_global"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "L": "L",
    },
    ("anti_de_sitter", "poincare"): {
        "t": "T", "x": "L", "y": "L", "z": "L", "L": "L",
    },

declaring $t$ a time, which is what tells the script that the chart multiplies it by $c$, and declaring the radius a length.
The radius being named $L$ is the one thing worth a second look, because the dimensional pass writes its own lengths with the letter $L$ as well.
The two never meet: the parameter is a symbol of the system's reader and the dimension is a symbol of the dimension table, they live in different expressions, and the pass reports the dimension of the parameter as $L$ because that is what the table was told to give it.
No relation in `PARAMETER_RELATIONS` is needed, since $L$ is free and every component is an identity in it.

Both systems check fast, the metrics being diagonal and their components short:

    /tmp/mfs-venv/bin/python _tools/derivations/verify_metrics.py \
        --system anti_de_sitter/static_global --system anti_de_sitter/poincare

runs in about three seconds and reports no disagreements, no dimensional failures and nothing unchecked.
