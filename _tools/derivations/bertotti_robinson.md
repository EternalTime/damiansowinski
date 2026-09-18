# The Bertotti-Robinson electrovacuum

This is the working behind the two coordinate systems in `MFS/assets/data/metrics/bertotti_robinson.json`.
Every number the entry prints is derived here, in order, from the line element down to the geodesic equations.
Nothing is left as an exercise and nothing is asserted that is not computed.

The companion script `_tools/derivations/verify_metrics.py` does the same work in sympy and compares it against the published file, so the algebra below is checkable by hand and by machine independently.

This is an electrovacuum: the Ricci scalar vanishes, the Ricci tensor does not, and the Einstein tensor is the Maxwell stress of a field that is the same at every event.
Reissner-Nordstrom is the other entry in the collection carrying a stress tensor of that kind, and the difference between the two is worth holding on to while reading.
There the field falls off with the radius and the curvature diverges at the centre; here the field does not fall off at all and the curvature is the same number everywhere.
Step 8 is where that correspondence is made, and it is the whole physical content of the entry.
Step 9 computes the Weyl tensor and gets zero, which is not a shortcut but a result, and finds that the same single condition on the two radii is responsible for the vanishing of the Ricci scalar and for the vanishing of Weyl.

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

settled for the collection on 2026-09-18.
Here the settlement costs something, because this Ricci tensor is not zero and the other contraction is minus it in every slot.
On the convention published, the energy density of the electromagnetic field comes out positive, which is the reason the convention was chosen: Step 8 gets $G^t{}_t = -1/b^2$ and reads it as $-8\pi G u/c^4$ with $u > 0$.
The Weyl tensor of Step 9 is built from the same contraction.

Factors of $G$ and $c$ are kept explicit, and so is $\epsilon_0$ where the electromagnetic field itself is written down.
The solution has one parameter, the length $b$, and no published component of either chart carries a factor of $c$, of $G$ or of $\epsilon_0$ at all.
Those constants enter only in Step 8, where $b$ is traded for the strength of the field that holds the geometry up.

The chart is the collection's, the one whose zeroth coordinate is

$$x^0 = ct.$$

The index is printed with the bare letter $t$, but the component printed against it is a component of that chart.
Because neither published metric depends on $t$, the rescaling leaves nothing behind, and every component below is a function of the radial coordinate and of $\theta$ alone.
The dots in the geodesic equations are velocities of that same chart, so $\dot{t}$ means $d(ct)/d\lambda$ with $\lambda$ an affine parameter, and Step 11 checks that this is the reading on which every published term balances.

---

## Step 2. Which chart, and why

The solution is the product of a two dimensional anti-de Sitter spacetime with a round two sphere, both of the same radius $b$,

$$\mathrm{AdS}_2(b) \times S^2(b),$$

and the choice of chart is a choice of chart for the $\mathrm{AdS}_2$ factor alone, since the sphere has only one.
The entry publishes two, and the first is the one the page opens on.

The first is the static chart,

$$ds^2 = -\frac{r^2}{b^2}c^2dt^2 + \frac{b^2}{r^2}dr^2 + b^2\left(d\theta^2 + \sin^2\theta\,d\phi^2\right),$$

which is the form the solution is usually met in.
Three things recommend it.
It is static and diagonal, with $\partial_t$ a Killing vector that is timelike on the whole of $r > 0$, so the components below are the ones a static observer measures.
It is the chart in which the spacetime appears as the near horizon limit of an extremal Reissner-Nordstrom black hole, worked out in Step 13, so $r$ is the coordinate that came from the parent hole's radial coordinate and $t$ the one that came from its time.
And it puts the whole of the electromagnetic field into a single component that does not depend on position, which is Step 8.

One warning has to go with it, and the entry's `convention` field carries it.
Neither $r$ here nor $x$ in the second chart is an areal radius.
Every sphere of this spacetime has the same area $4\pi b^2$, because $g_{\theta\theta} = b^2$ is a constant, so the areal radius is $b$ everywhere and is not a coordinate at all.
What $r$ measures is position along the $\mathrm{AdS}_2$ factor, and the surface $r = 0$ is a degenerate Killing horizon rather than a centre.

The second chart is the Poincaré chart,

$$ds^2 = \frac{b^2}{x^2}\left(-c^2dt^2 + dx^2\right) + b^2\left(d\theta^2 + \sin^2\theta\,d\phi^2\right),$$

reached from the first by $rx = b^2$, which Step 12 carries out.
It is published for three reasons.
The $\mathrm{AdS}_2$ factor is manifestly a conformal multiple of a flat two dimensional metric in it, which is as close as any chart comes to displaying by eye the conformal flatness that Step 9 proves of the whole four dimensional spacetime.
It is the chart the solution is used in, since the Poincaré patch of $\mathrm{AdS}_2$ is where near horizon black hole thermodynamics and two dimensional holography do their work.
And it gives a second and independent reading of the same invariants: $R = 0$, $C_{\mu\nu\rho\sigma} = 0$ and $K = 8/b^4$ come out of two different sets of Christoffel symbols, which is a check on both.

Neither chart is complete.
Both cover a Poincaré patch of the $\mathrm{AdS}_2$ factor, bounded by a horizon at $r = 0$, equivalently $x = \infty$, and the maximal extension is the global $\mathrm{AdS}_2$ cylinder.
Nothing is lost by stopping there, because Step 10 finds the curvature constant: there is no singularity anywhere in this spacetime to extend towards.

---

## Step 3. The metric, its determinant and its inverse

Order the coordinates $(t, r, \theta, \phi)$ and write the static chart with $x^0 = ct$:

$$g_{\mu\nu} = \mathrm{diag}\left(-\frac{r^2}{b^2},\ \frac{b^2}{r^2},\ b^2,\ b^2\sin^2\theta\right).$$

The first two entries are dimensionless, as they must be in a chart whose zeroth coordinate is $ct$ and whose first is a length, and the last two carry $L^2$ against the two angles.
Note the shape of it: the first two entries are reciprocals of each other up to sign, and the last two do not mention $r$.
Both facts are the product structure showing itself, and both are used below.

The determinant is a product of the four:

$$\det g = -\frac{r^2}{b^2}\cdot\frac{b^2}{r^2}\cdot b^2\cdot b^2\sin^2\theta = -b^4\sin^2\theta, \qquad \sqrt{-g} = b^2\sin\theta.$$

It does not depend on $r$ at all.
That is the volume element of a spacetime whose spheres never change size, and it is the first sign that $r$ is not a radius.

The inverse is the entrywise reciprocal, since the metric is diagonal:

$$g^{\mu\nu} = \mathrm{diag}\left(-\frac{b^2}{r^2},\ \frac{r^2}{b^2},\ \frac{1}{b^2},\ \frac{1}{b^2\sin^2\theta}\right).$$

The entry publishes both, and the checker confirms that the published pair multiply to the identity.

---

## Step 4. The Christoffel symbols

Nothing depends on $t$ or on $\phi$, so only $\partial_r$ and $\partial_\theta$ ever act, and the metric is diagonal, so each symbol is one term.
For a diagonal metric,

$$\Gamma^\mu{}_{\mu\nu} = \tfrac{1}{2}\partial_\nu\ln|g_{\mu\mu}|, \qquad \Gamma^\mu{}_{\nu\nu} = -\tfrac{1}{2}g^{\mu\mu}\partial_\mu g_{\nu\nu} \quad (\mu \neq \nu),$$

with no summation in either.

Take them in turn.
From $g_{tt} = -r^2/b^2$,

$$\Gamma^t{}_{tr} = \Gamma^t{}_{rt} = \tfrac{1}{2}\partial_r\ln\frac{r^2}{b^2} = \frac{1}{r},$$

$$\Gamma^r{}_{tt} = -\tfrac{1}{2}g^{rr}\partial_r g_{tt} = -\tfrac{1}{2}\cdot\frac{r^2}{b^2}\cdot\left(-\frac{2r}{b^2}\right) = \frac{r^3}{b^4}.$$

From $g_{rr} = b^2/r^2$,

$$\Gamma^r{}_{rr} = \tfrac{1}{2}\partial_r\ln\frac{b^2}{r^2} = -\frac{1}{r}.$$

The two symbols that would carry the sphere into the radial direction are

$$\Gamma^r{}_{\theta\theta} = -\tfrac{1}{2}g^{rr}\partial_r g_{\theta\theta} = 0, \qquad \Gamma^r{}_{\phi\phi} = -\tfrac{1}{2}g^{rr}\partial_r g_{\phi\phi} = 0,$$

because $g_{\theta\theta}$ and $g_{\phi\phi}$ do not depend on $r$, and for the same reason

$$\Gamma^\theta{}_{r\theta} = \tfrac{1}{2}\partial_r\ln b^2 = 0, \qquad \Gamma^\phi{}_{r\phi} = \tfrac{1}{2}\partial_r\ln\left(b^2\sin^2\theta\right) = 0.$$

Those four zeros are the ones a reader who knows Schwarzschild will look for, where they are $-rf$, $-rf\sin^2\theta$ and $1/r$ twice.
Their absence here is the product structure, and Step 5 says what follows from it.

The sphere keeps its own two symbols,

$$\Gamma^\theta{}_{\phi\phi} = -\sin\theta\cos\theta, \qquad \Gamma^\phi{}_{\theta\phi} = \Gamma^\phi{}_{\phi\theta} = \cot\theta,$$

which are those of a round sphere of any radius, since a constant factor in $g_{\theta\theta}$ and $g_{\phi\phi}$ cancels between the metric and its inverse.

So seven symbols are nonzero, counting the reflections in the lower pair, and the entry publishes all seven.
Lowering the first index with $\Gamma_{\mu\nu\rho} = g_{\mu\mu}\Gamma^\mu{}_{\nu\rho}$ multiplies each by one diagonal entry:

$$\Gamma_{ttr} = \Gamma_{trt} = -\frac{r^2}{b^2}\cdot\frac{1}{r} = -\frac{r}{b^2}, \qquad \Gamma_{rtt} = \frac{b^2}{r^2}\cdot\frac{r^3}{b^4} = \frac{r}{b^2}, \qquad \Gamma_{rrr} = \frac{b^2}{r^2}\cdot\left(-\frac{1}{r}\right) = -\frac{b^2}{r^3},$$

$$\Gamma_{\theta\phi\phi} = -b^2\sin\theta\cos\theta, \qquad \Gamma_{\phi\theta\phi} = \Gamma_{\phi\phi\theta} = b^2\sin^2\theta\cot\theta = b^2\sin\theta\cos\theta.$$

Both variants are in the entry.

The factors of $c$ are worth one line, because they are where an entry in this collection is most easily got wrong.
Computed with the bare $t$, the symbol $\Gamma^r{}_{tt}$ carries $c^2$ and $\Gamma^t{}_{tr}$ carries none.
In the chart $x^0 = ct$ a symbol is multiplied by $c$ once per upper time index and divided by $c$ once per lower one, so $\Gamma^r{}_{tt}$ loses its $c^2$ and $\Gamma^t{}_{tr}$, which has one upper time index and one lower, is unchanged.
Every $c$ cancels, which is why none appears in the list.

---

## Step 5. The product structure, read off the connection

Split the indices into the two factors: let $a, b$ run over $(t, r)$ and $i, j$ over $(\theta, \phi)$.
Step 4 found that every mixed symbol vanishes,

$$\Gamma^a{}_{ij} = 0, \qquad \Gamma^i{}_{aj} = 0, \qquad \Gamma^a{}_{ib} = 0, \qquad \Gamma^i{}_{ab} = 0,$$

and that the surviving symbols are exactly those of the two factors taken separately.
The connection is therefore the direct sum of the connection of $\mathrm{AdS}_2$ and the connection of $S^2$, and so is everything built from it.
Concretely, Riemann has no component with indices drawn from both blocks, the two blocks below are the Riemann tensors of the factors, and a curve that starts in one factor stays in it.

That is what makes this spacetime easy and is also what makes it strange.
The sphere is rigid: no motion, no wave and no observer changes its area, because there is no connection component that could tilt a radial direction into an angular one.

Each factor is a two dimensional space of constant curvature, and in two dimensions the whole Riemann tensor is one number, the Gaussian curvature $K$:

$$R_{abcd} = K\left(g_{ac}g_{bd} - g_{ad}g_{bc}\right).$$

For the $\mathrm{AdS}_2$ factor $K_1 = -1/b^2$ and for the sphere $K_2 = +1/b^2$.
Step 6 confirms both by direct computation rather than by citation.

---

## Step 6. The Riemann tensor

Only two independent components survive, one per factor.

For the $\mathrm{AdS}_2$ block take $\mu = t$, $\nu = r$, $\rho = t$, $\sigma = r$ in the definition.
The derivative terms are

$$\partial_t\Gamma^t{}_{rr} - \partial_r\Gamma^t{}_{rt} = 0 - \partial_r\frac{1}{r} = \frac{1}{r^2},$$

and the quadratic terms are

$$\Gamma^t{}_{t\lambda}\Gamma^\lambda{}_{rr} - \Gamma^t{}_{r\lambda}\Gamma^\lambda{}_{rt} = \Gamma^t{}_{tr}\Gamma^r{}_{rr} - \Gamma^t{}_{rt}\Gamma^t{}_{tr} = \frac{1}{r}\cdot\left(-\frac{1}{r}\right) - \frac{1}{r}\cdot\frac{1}{r} = -\frac{2}{r^2}.$$

Adding them,

$$R^t{}_{rtr} = \frac{1}{r^2} - \frac{2}{r^2} = -\frac{1}{r^2}.$$

Lowering with $g_{tt}$,

$$R_{trtr} = -\frac{r^2}{b^2}\cdot\left(-\frac{1}{r^2}\right) = \frac{1}{b^2},$$

and since $g_{tt}g_{rr} = -1$ in this chart, the two dimensional formula gives

$$R_{trtr} = K_1\left(g_{tt}g_{rr}\right) = -K_1,$$

so $K_1 = -1/b^2$ as claimed.
Raising the other way,

$$R^r{}_{trt} = g^{rr}R_{rtrt} = \frac{r^2}{b^2}\cdot\frac{1}{b^2} = \frac{r^2}{b^4}.$$

For the sphere block the two dimensional formula and the standard symbols give

$$R_{\theta\phi\theta\phi} = K_2\,g_{\theta\theta}g_{\phi\phi} = \frac{1}{b^2}\cdot b^2\cdot b^2\sin^2\theta = b^2\sin^2\theta,$$

$$R^\theta{}_{\phi\theta\phi} = g^{\theta\theta}R_{\theta\phi\theta\phi} = \sin^2\theta, \qquad R^\phi{}_{\theta\phi\theta} = g^{\phi\phi}R_{\phi\theta\phi\theta} = 1.$$

Every component not obtainable from these two by the symmetries of Riemann is zero, and in particular every component with one index in each block is zero, which is Step 5.

The published blocks list the eight nonzero components of each variant, which are these two together with the sign flips forced by antisymmetry in the last pair and by the pair exchange:

$$R^t{}_{rtr} = -\frac{1}{r^2} = -R^t{}_{rrt}, \qquad R^r{}_{trt} = \frac{r^2}{b^4} = -R^r{}_{ttr},$$

$$R_{trtr} = R_{rtrt} = \frac{1}{b^2} = -R_{trrt} = -R_{rttr},$$

$$R^\theta{}_{\phi\theta\phi} = \sin^2\theta = -R^\theta{}_{\phi\phi\theta}, \qquad R^\phi{}_{\theta\phi\theta} = 1 = -R^\phi{}_{\theta\theta\phi},$$

$$R_{\theta\phi\theta\phi} = R_{\phi\theta\phi\theta} = b^2\sin^2\theta = -R_{\theta\phi\phi\theta} = -R_{\phi\theta\theta\phi}.$$

Two remarks before going on.
The lowered components carry no $r$ at all: the curvature of this spacetime is the same at every event, which is the homogeneity the history section claims, made arithmetic.
And the two blocks carry opposite signs in the same normalisation, $-1/b^2$ against $+1/b^2$, which is the balance the next three steps live on.

---

## Step 7. The Ricci tensor, and why the scalar vanishes

Contract on the first lower index.
Because Riemann is block diagonal, the sum for an index in one factor runs over that factor alone.

In the $\mathrm{AdS}_2$ block,

$$R_{tt} = R^\alpha{}_{t\alpha t} = R^r{}_{trt} = \frac{r^2}{b^4}, \qquad R_{rr} = R^t{}_{rtr} = -\frac{1}{r^2},$$

and both are $K_1$ times the metric,

$$R_{tt} = -\frac{1}{b^2}\left(-\frac{r^2}{b^2}\right), \qquad R_{rr} = -\frac{1}{b^2}\cdot\frac{b^2}{r^2}.$$

In the sphere block,

$$R_{\theta\theta} = R^\phi{}_{\theta\phi\theta} = 1 = \frac{1}{b^2}\cdot b^2, \qquad R_{\phi\phi} = R^\theta{}_{\phi\theta\phi} = \sin^2\theta = \frac{1}{b^2}\cdot b^2\sin^2\theta.$$

So in one line,

$$R_{ab} = K_1 g_{ab} = -\frac{1}{b^2}g_{ab}, \qquad R_{ij} = K_2 g_{ij} = +\frac{1}{b^2}g_{ij}, \qquad R_{ai} = 0,$$

which is the general statement that a product of two dimensional factors is an Einstein space in each factor separately but not as a whole.

The mixed form is the cleanest way to see it, and the entry publishes it:

$$R^t{}_t = R^r{}_r = -\frac{1}{b^2}, \qquad R^\theta{}_\theta = R^\phi{}_\phi = +\frac{1}{b^2}.$$

The Ricci scalar is the trace of that,

$$R = 2K_1 + 2K_2 = -\frac{2}{b^2} + \frac{2}{b^2} = 0.$$

It vanishes because the two radii are equal, and for no other reason.
Had the sphere had radius $b_2 \neq b_1$ the scalar would have been $2(1/b_2^2 - 1/b_1^2)$, which is the Nariai family of products and is not an electrovacuum.
The entry publishes $R = 0$ and the front end prints it as vanishing.

The Ricci tensor itself is emphatically not zero, and this is the place to say why the collection's contraction convention matters here when it did not for the vacuum entries.
The other contraction, on the last index, gives $-R_{\mu\nu}$, and with it $R^t{}_t$ would read $+1/b^2$ and the energy density of Step 8 would come out negative.
Nothing in the geometry changes, but the reading of it does, and the collection reads it the way that makes ordinary matter positive.

---

## Step 8. The Einstein tensor is the Maxwell stress of a uniform field

Since $R = 0$, the Einstein tensor is the Ricci tensor:

$$G_{\mu\nu} = R_{\mu\nu} - \tfrac{1}{2}Rg_{\mu\nu} = R_{\mu\nu},$$

so the entry's two blocks are componentwise identical, and they are published that way rather than one of them omitted.
In mixed form,

$$G^t{}_t = G^r{}_r = -\frac{1}{b^2}, \qquad G^\theta{}_\theta = G^\phi{}_\phi = +\frac{1}{b^2}, \qquad G^\mu{}_\mu = 0.$$

That pattern, $\mathrm{diag}(-1,-1,+1,+1)$, is the signature of a radial electromagnetic field and of nothing else.
Read it as a stress tensor through $G_{\mu\nu} = 8\pi G T_{\mu\nu}/c^4$, with $T^t{}_t = -u$ the energy density, $T^r{}_r = p_r$ the radial pressure and $T^\theta{}_\theta = T^\phi{}_\phi = p_\perp$ the transverse pressure:

$$u = \frac{c^4}{8\pi G b^2} > 0, \qquad p_r = -u, \qquad p_\perp = +u.$$

A positive energy density, a radial tension of exactly the same size, a transverse pressure of the same size again, and a vanishing trace.
That is the stress of a static radial electric field, and the trace is zero because the Maxwell stress tensor is traceless in four dimensions, which is the same fact as $R = 0$ in Step 7.

Now produce the field itself rather than argue from the pattern.
Take

$$F = -\frac{E}{c}\,d(ct)\wedge dr, \qquad \text{that is} \qquad F_{0r} = -F_{r0} = -\frac{E}{c}$$

in the chart $x^0 = ct$, with $E$ a constant and every other component zero.
It comes from the potential $A_0 = Er/c$, which is a potential growing linearly with position: the uniform field of a parallel plate capacitor, not the field of a point charge.

Three properties, in order.

First, the field a static observer measures is $E$ at every event.
The orthonormal coframe of the static chart is

$$e^0 = \frac{r}{b}d(ct), \qquad e^1 = \frac{b}{r}dr, \qquad e^2 = b\,d\theta, \qquad e^3 = b\sin\theta\,d\phi,$$

so that $g = -e^0\otimes e^0 + e^1\otimes e^1 + e^2\otimes e^2 + e^3\otimes e^3$, and since $d(ct)\wedge dr = e^0\wedge e^1$ exactly, the factors of $r$ cancelling between the two legs,

$$F = -\frac{E}{c}\,e^0\wedge e^1.$$

The frame component is the constant $-E/c$ wherever the observer stands.
This is the precise sense in which the field is uniform, and it is why the geometry can be homogeneous: a field that fell off with position would need a position for it to fall off from.

Second, the field is covariantly constant,

$$\nabla_\mu F_{\nu\rho} = \partial_\mu F_{\nu\rho} - \Gamma^\alpha{}_{\mu\nu}F_{\alpha\rho} - \Gamma^\alpha{}_{\mu\rho}F_{\nu\alpha} = 0$$

in all sixty four slots.
Only $\mu = r$ and $\mu = t$ can give anything: the components of $F$ are constants, so no derivative term survives, the only occupied rows of $F$ are $t$ and $r$, and Step 5 has already removed every symbol that could connect those rows to the sphere.
For $\mu = r$,

$$\nabla_r F_{tr} = -\Gamma^t{}_{rt}F_{tr} - \Gamma^r{}_{rr}F_{tr} = -\frac{1}{r}F_{tr} + \frac{1}{r}F_{tr} = 0,$$

the two symbols of Step 4 cancelling because $\Gamma^t{}_{tr} = -\Gamma^r{}_{rr}$, which is in turn because $g_{tt}g_{rr}$ is a constant.
For $\mu = t$,

$$\nabla_t F_{tr} = -\Gamma^\alpha{}_{tt}F_{\alpha r} - \Gamma^\alpha{}_{tr}F_{t\alpha} = -\Gamma^r{}_{tt}F_{rr} - \Gamma^t{}_{tr}F_{tt} = 0,$$

both terms dying on the antisymmetry of $F$.
A covariantly constant field is the strongest sense in which a field can be the same everywhere, and it is the defining property of this solution.
Its two invariants follow at once and are constants:

$$F_{\mu\nu}F^{\mu\nu} = -\frac{2E^2}{c^2}, \qquad F_{\mu\nu}\tilde{F}^{\mu\nu} = 0,$$

the first negative because the field is purely electric, the second zero for the same reason.

Third, its stress tensor is the one the Einstein tensor asked for.
With

$$T_{\mu\nu} = \frac{1}{\mu_0}\left(F_{\mu\alpha}F_\nu{}^\alpha - \tfrac{1}{4}g_{\mu\nu}F_{\alpha\beta}F^{\alpha\beta}\right), \qquad \mu_0 = \frac{1}{\epsilon_0c^2},$$

the frame components are

$$T_{\hat{0}\hat{0}} = \frac{1}{\mu_0}\left(\frac{E^2}{c^2} - \frac{1}{4}\cdot\frac{2E^2}{c^2}\right) = \frac{\epsilon_0E^2}{2}, \qquad T_{\hat{1}\hat{1}} = -\frac{\epsilon_0E^2}{2}, \qquad T_{\hat{2}\hat{2}} = T_{\hat{3}\hat{3}} = \frac{\epsilon_0E^2}{2},$$

that is

$$T^\mu{}_\nu = \frac{\epsilon_0E^2}{2}\,\mathrm{diag}(-1,-1,+1,+1),$$

with the familiar $u = \epsilon_0E^2/2$ in front.
Comparing slot by slot with the Einstein tensor above, all four equations of $G^\mu{}_\nu = 8\pi G T^\mu{}_\nu/c^4$ reduce to the single condition

$$\frac{1}{b^2} = \frac{4\pi G\epsilon_0E^2}{c^4}, \qquad\text{that is}\qquad b = \frac{c^2}{E\sqrt{4\pi G\epsilon_0}}, \qquad b^2 = \frac{c^4}{8\pi G u}.$$

That is the whole content of the entry's parameter.
$b$ is not free once the field is given: the radius of the sphere and the radius of the anti-de Sitter factor are both fixed by the strength of the field, and a stronger field makes a tighter throat.
Run the numbers and the field has to be violent, because $c^4/G$ is enormous: a throat a kilometre across takes about $10^{24}$ volts per metre, and an energy density of some $5\times10^{36}$ joules per cubic metre to go with it.

Two footnotes to this step.
The magnetic case is the same geometry.
Replacing $F$ by its dual, $F = B b^2\sin\theta\,d\theta\wedge d\phi$, which is a monopole field threading the sphere with constant flux, gives the same stress tensor with $\epsilon_0E^2/2$ replaced by $B^2/2\mu_0$, so it changes nothing above and the entry's geometry covers both.
And the electrovacuum condition is what forced the two radii to agree in the first place.
For a product $\mathrm{AdS}_2(b_1)\times S^2(b_2)$ with the radii left free, the same computation gives

$$G^t{}_t = G^r{}_r = -\frac{1}{b_2^2}, \qquad G^\theta{}_\theta = G^\phi{}_\phi = +\frac{1}{b_1^2},$$

whose trace vanishes only when $b_1 = b_2$.
A Maxwell field has a traceless stress tensor, so it can only source a product whose radii agree, and the equal radii the entry publishes are a consequence of the source being electromagnetic and not an extra assumption.

---

## Step 9. The Weyl tensor, computed and zero

In four dimensions,

$$C_{\mu\nu\rho\sigma} = R_{\mu\nu\rho\sigma} - \tfrac{1}{2}\left(g_{\mu\rho}R_{\sigma\nu} - g_{\mu\sigma}R_{\rho\nu} - g_{\nu\rho}R_{\sigma\mu} + g_{\nu\sigma}R_{\rho\mu}\right) + \tfrac{1}{6}R\left(g_{\mu\rho}g_{\sigma\nu} - g_{\mu\sigma}g_{\rho\nu}\right).$$

This is not a vacuum, so Weyl is not Riemann and cannot be copied from it.
It is computed here in the general case, with the two radii left free, so that the answer says something rather than merely being zero.

Write $K_1$ and $K_2$ for the two Gaussian curvatures, $h$ for the metric of the first factor and $k$ for that of the second, so that Step 5 and Step 7 give

$$R_{abcd} = K_1\left(h_{ac}h_{bd} - h_{ad}h_{bc}\right), \qquad R_{ijkl} = K_2\left(k_{ik}k_{jl} - k_{il}k_{jk}\right), \qquad R_{ab} = K_1h_{ab}, \qquad R_{ij} = K_2k_{ij},$$

with every mixed component of Riemann zero and $R = 2(K_1+K_2)$.

Take the three kinds of slot in turn.

Both index pairs in the first factor.
The trace correction is

$$\tfrac{1}{2}\left(h_{ac}R_{db} - h_{ad}R_{cb} - h_{bc}R_{da} + h_{bd}R_{ca}\right) = K_1\left(h_{ac}h_{bd} - h_{ad}h_{bc}\right),$$

which is the whole of $R_{abcd}$, and the scalar term leaves

$$C_{abcd} = \tfrac{1}{3}\left(K_1+K_2\right)\left(h_{ac}h_{bd} - h_{ad}h_{bc}\right).$$

Both index pairs in the second factor.
The same computation with $k$ and $K_2$ gives

$$C_{ijkl} = \tfrac{1}{3}\left(K_1+K_2\right)\left(k_{ik}k_{jl} - k_{il}k_{jk}\right).$$

One pair in each.
Here Riemann contributes nothing, and the two surviving terms of the correction carry opposite signs:

$$C_{aibj} = 0 - \tfrac{1}{2}\left(h_{ab}R_{ji} + k_{ij}R_{ba}\right) + \tfrac{1}{6}\cdot 2\left(K_1+K_2\right)h_{ab}k_{ij} = -\tfrac{1}{6}\left(K_1+K_2\right)h_{ab}k_{ij}.$$

So every component of the Weyl tensor of a product of two dimensional factors is proportional to $K_1 + K_2$, and the spacetime is conformally flat exactly when the two curvatures are equal and opposite.
For Bertotti-Robinson $K_1 = -1/b^2$ and $K_2 = +1/b^2$, so

$$C_{\mu\nu\rho\sigma} = 0$$

in every slot, and both published Weyl blocks are empty, which the front end prints as all vanish.

The condition $K_1 + K_2 = 0$ is the same condition that made the Ricci scalar vanish in Step 7 and the same condition the traceless Maxwell stress forced in Step 8.
So for this family the three statements are one statement: the source is electromagnetic, the radii agree, and the spacetime is conformally flat.
Take the radii apart and all three fail together.
That is worth saying because conformal flatness is rare, and because everywhere else in the collection it is bought with symmetry or with a fluid.
Minkowski publishes an empty Weyl tensor by being flat outright and anti-de Sitter by being maximally symmetric, so neither has anything left to remove; Friedmann Robertson Walker and the interior Schwarzschild solution publish one on the strength of a perfect fluid and a preferred rest frame.
This entry is the only one that gets there through its source being electromagnetic, and the only one where the condition for it is an equality between two curvatures that were free to differ.

An empty Weyl block is a claim, not an omission, and it is the claim the checker is strictest about: it requires every component the entry does not print to vanish, so a Weyl tensor published as empty is checked in all two hundred and fifty six slots.

---

## Step 10. The Kretschmann scalar

For a two dimensional space of Gaussian curvature $K$, the fully raised Riemann tensor is $R^{abcd} = K(h^{ac}h^{bd} - h^{ad}h^{bc})$, and contracting it with the lowered one gives

$$R_{abcd}R^{abcd} = K^2\left(4 - 2 - 2 + 4\right) = 4K^2,$$

where each $4$ is $\delta^a_a\delta^b_b$ over a two dimensional block and each $2$ is $\delta^a_b\delta^b_a$.
Riemann is block diagonal, so the two blocks add with no cross term:

$$K_{\text{Kretschmann}} = R_{\mu\nu\rho\sigma}R^{\mu\nu\rho\sigma} = 4K_1^2 + 4K_2^2 = \frac{4}{b^4} + \frac{4}{b^4} = \frac{8}{b^4}.$$

The entry publishes

$$K = \frac{8}{b^4}.$$

It is a constant, positive, finite everywhere and free of every coordinate.
There is no curvature singularity anywhere in this spacetime, and the surface $r = 0$ where $g_{tt}$ vanishes is a horizon and nothing worse, as one expects of the throat of an extremal hole.
Note also that the two factors contribute equally, $4/b^4$ each, and that the scalar cannot tell the negative curvature from the positive one, since it squares them: $K$ alone does not distinguish this spacetime from a product of two spheres of radius $b$.
It takes the Ricci tensor of Step 7, which is not squared, to see the difference.

---

## Step 11. The geodesic equations

The equation is

$$\ddot{x}^\mu + \Gamma^\mu{}_{\nu\rho}\dot{x}^\nu\dot{x}^\rho = 0,$$

with the dot the derivative along an affine parameter $\lambda$ and every coordinate the chart's, so $\dot{t} = d(ct)/d\lambda$.
Substituting the seven symbols of Step 4, with the factor of two coming from the sum over the two orderings of a symmetric pair,

$$\ddot{t} + \frac{2}{r}\dot{t}\dot{r} = 0,$$

$$\ddot{r} + \frac{r^3}{b^4}\dot{t}^2 - \frac{1}{r}\dot{r}^2 = 0,$$

$$\ddot{\theta} - \sin\theta\cos\theta\,\dot{\phi}^2 = 0,$$

$$\ddot{\phi} + 2\cot\theta\,\dot{\theta}\dot{\phi} = 0,$$

which is what the entry publishes.

The two halves do not talk to each other, which is Step 5 again and is the most striking thing about motion here.
The last two equations are the geodesic equations of a sphere of radius $b$ with no reference to $r$ or $t$ whatever, so the angular motion of any free particle is uniform motion along a great circle of a fixed sphere, forever.
The first two are the geodesic equations of $\mathrm{AdS}_2$ with no reference to the angles, so the radial motion is unaffected by how fast the particle is going around.
There are no orbits in this spacetime in the sense Schwarzschild has them: nothing circles a centre, because there is no centre, and the angular momentum of a particle exerts no centrifugal effect on its radial motion at all.

The dimensional check is the one the collection insists on, and the second equation is where the reading of the dot is forced.
Every term of the first equation carries $L/\lambda^2$: $\ddot{t}$ does because the chart coordinate $ct$ is a length, and $\dot{t}\dot{r}/r$ does because it is $(L/\lambda)(L/\lambda)/L$.
That equation would balance on the bare reading too, both of its terms simply carrying $T/\lambda^2$ instead, so it decides nothing.
The second does.
There $\ddot{r}$ and $\dot{r}^2/r$ carry $L/\lambda^2$ whatever the dot on $t$ means, while $r^3\dot{t}^2/b^4$ carries $(L^3/L^4)(L^2/\lambda^2) = L/\lambda^2$ only because $\dot{t}$ is $d(ct)/d\lambda$ and so carries a length over the affine parameter.
Read as a bare $dt/d\lambda$ that one term would come out as $T^2/(L\lambda^2)$, short of its neighbours by exactly the $c^2$ the chart supplies, and the equation would not balance.
In the last two equations the left hand side carries $1/\lambda^2$, since an angle is dimensionless, and so does every term beside it.

---

## Step 12. The Poincaré chart

Substitute

$$r = \frac{b^2}{x}, \qquad dr = -\frac{b^2}{x^2}dx,$$

which is an orientation reversing map of $r > 0$ onto $x > 0$, so the deep throat $r \to 0$ is $x \to \infty$.
The two radial terms become

$$-\frac{r^2}{b^2}c^2dt^2 = -\frac{b^2}{x^2}c^2dt^2, \qquad \frac{b^2}{r^2}dr^2 = \frac{b^2}{b^4/x^2}\cdot\frac{b^4}{x^4}dx^2 = \frac{b^2}{x^2}dx^2,$$

and the sphere is untouched, so

$$ds^2 = \frac{b^2}{x^2}\left(-c^2dt^2 + dx^2\right) + b^2\left(d\theta^2 + \sin^2\theta\,d\phi^2\right).$$

The $\mathrm{AdS}_2$ factor is now visibly $b^2/x^2$ times the flat two dimensional metric, which is the Poincaré form.

Everything is recomputed from scratch in this chart rather than transported, so that the two systems are independent checks on each other.
The metric is $\mathrm{diag}(-b^2/x^2,\ b^2/x^2,\ b^2,\ b^2\sin^2\theta)$ and the inverse the entrywise reciprocal.
The connection comes out unusually uniform:

$$\Gamma^t{}_{tx} = \Gamma^x{}_{tt} = \Gamma^x{}_{xx} = -\frac{1}{x},$$

all three equal, with the sphere's two symbols unchanged, and lowering gives

$$\Gamma_{ttx} = \frac{b^2}{x^3}, \qquad \Gamma_{xtt} = \Gamma_{xxx} = -\frac{b^2}{x^3}.$$

The curvature follows as before:

$$R^t{}_{xtx} = R^x{}_{ttx} = -\frac{1}{x^2}, \qquad R_{txtx} = \frac{b^2}{x^4}, \qquad R_{tt} = -R_{xx} = \frac{1}{x^2},$$

with the sphere block identical to Step 6, and the mixed Ricci tensor is the same four numbers as before,

$$R^t{}_t = R^x{}_x = -\frac{1}{b^2}, \qquad R^\theta{}_\theta = R^\phi{}_\phi = +\frac{1}{b^2}.$$

The geodesic equations are

$$\ddot{t} - \frac{2}{x}\dot{t}\dot{x} = 0, \qquad \ddot{x} - \frac{1}{x}\dot{t}^2 - \frac{1}{x}\dot{x}^2 = 0,$$

with the two angular equations unchanged.

The invariants agree with Step 7, Step 9 and Step 10 exactly: $R = 0$, every component of Weyl vanishes, and $K = 8/b^4$.
They had to, being invariants, and the point of publishing the second chart is that they do so through a different set of Christoffel symbols and so check the first.

---

## Step 13. Where the geometry comes from

The entry's history calls this spacetime the throat of an extremal Reissner-Nordstrom hole.
That is a statement about a limit and is worth doing, because it is where the parameter $b$ gets its physical meaning.

The extremal Reissner-Nordstrom solution has $r_s = 2r_q$ and so a double horizon at $\rho = r_q$:

$$ds^2 = -\left(1 - \frac{r_q}{\rho}\right)^2c^2dt^2 + \frac{d\rho^2}{\left(1 - \frac{r_q}{\rho}\right)^2} + \rho^2\left(d\theta^2 + \sin^2\theta\,d\phi^2\right).$$

Zoom in on the horizon by

$$\rho = r_q + \epsilon r, \qquad t = \frac{\tau}{\epsilon},$$

and let $\epsilon \to 0$.
Since

$$1 - \frac{r_q}{\rho} = \frac{\epsilon r}{r_q + \epsilon r} \longrightarrow \frac{\epsilon r}{r_q},$$

the time term goes to $-(r^2/r_q^2)c^2d\tau^2$, the two factors of $\epsilon$ cancelling against the two in $dt^2$; the radial term goes to $(r_q^2/r^2)dr^2$, the two factors of $\epsilon$ cancelling the other way; and the sphere goes to $r_q^2d\Omega^2$, frozen at the horizon radius.
The limit is exactly the static chart of Step 2 with

$$b = r_q, \qquad r_q^2 = \frac{GQ^2}{4\pi\epsilon_0c^4},$$

the charge radius the Reissner-Nordstrom entry defines.
The two readings of $b$ agree: the electric field at the horizon of an extremal hole is $E = Q/4\pi\epsilon_0r_q^2$, and putting that into the $b$ of Step 8 returns $b = r_q$ identically.

So the throat is not an analogy.
An extremal charged hole has an infinitely long neck of constant circumference, and the geometry of that neck, on its own, is this entry.
The distance down it is the coordinate $r$, the constant circumference is the $4\pi b^2$ of Step 3, and the uniform field of Step 8 is the hole's own field, which stops falling off once the sphere stops shrinking.

---

## Step 14. What the checker checks

Run

    python3 -m venv /tmp/mfs-venv && /tmp/mfs-venv/bin/pip install sympy
    /tmp/mfs-venv/bin/python _tools/derivations/verify_metrics.py --system bertotti_robinson/static --system bertotti_robinson/poincare

and it reports both systems checked with no disagreement and no dimensional failure, in about three seconds.

It reads $g_{\mu\nu}$ off each line element, computes the inverse, both Christoffel variants, Riemann in both variants, Ricci in three, the Ricci scalar, Kretschmann, Einstein in three and Weyl in two, weights every component into the chart $x^0 = ct$, and compares it against what the entry prints.
Anything the entry omits has to vanish, which is what makes the two empty Weyl blocks and the many zero slots of the other blocks into checked claims rather than silence.
It also confirms that the published metric and the published inverse multiply to the identity, and that each published geodesic equation is $\ddot{x}^\mu + \Gamma^\mu{}_{\nu\rho}\dot{x}^\nu\dot{x}^\rho = 0$ for the published symbols.

Both systems are declared in `DIMENSIONS` in that script, as $t$ a time, $r$ and $x$ lengths, $\theta$ and $\phi$ angles and $b$ a length.
That declaration is what tells the checker which coordinate the chart multiplies by $c$, and it is also what the dimensional pass weighs every published term against.
There is no `PARAMETER_RELATIONS` line, because $b$ is a free parameter of the family: fixing it needs the field of Step 8, which is not part of the metric, so every published value here is an identity in $b$.
