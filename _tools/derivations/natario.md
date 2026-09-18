# The Natário warp drive

This is the working behind the two coordinate systems in `MFS/assets/data/metrics/natario.json`.
Every number the entry prints is derived here, in order, from the line element down to the geodesic equations.
Nothing is left as an exercise and nothing is asserted that is not computed.

The companion script `_tools/derivations/verify_metrics.py` does the same work in sympy and compares it against the published file, so the algebra below is checkable by hand and by machine independently.

This entry publishes a construction rather than a solution.
There is no field equation to solve here and no parameter to fit: a flow field is chosen on flat space, the line element is built from it, and whatever stress energy the Einstein tensor then reports is what that flow would cost to maintain.
That is the honest way to read every warp drive, and it is what makes the arithmetic below worth doing, because the cost comes out negative and stays negative no matter which flow is chosen.

Two things are settled here that the Alcubierre entry cannot settle.
Step 4 shows that the expansion of the drive is one number, the divergence of the flow, and that Natário's condition is the vanishing of exactly that number.
Step 9 then computes the energy density with the divergence left free and finds it to be a difference of two squares, of which Natário's condition removes the positive one.
Taking the expansion away is therefore not a repair.
It removes the only term in the energy density that was ever positive and leaves the negative one untouched, which is why the drive is interesting rather than fixed.

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
The choice is not free here.
Neither chart is a vacuum, so the two contractions differ in every slot of the Ricci tensor, of the Ricci scalar and of the Einstein tensor, and the sign of the energy density in Step 9 is the sign this convention gives.
It is the convention on which ordinary matter has a positive energy density, which is what makes the negative one found below a statement about the drive rather than about the bookkeeping.

The chart is the collection's, the one whose zeroth coordinate is

$$x^0 = ct.$$

The index is printed with the bare letter $t$, but the component printed against it is a component of that chart.
A printed partial derivative is taken with respect to the same chart coordinate,

$$\partial_t = \frac{1}{c}\frac{\partial}{\partial t},$$

while $\partial_x$, $\partial_y$ and $\partial_z$ are what they look like.
That single factor is what makes $\partial_t u$ and $\partial_x u$ carry the same dimension, $1/T$, so that they can be added.
The dots in the geodesic equations are velocities of that same chart, so $\dot{t}$ means $d(ct)/d\lambda$, and Step 21 checks that this is the reading on which every published term balances.

Two coordinate systems are published, and they are not two charts of one spacetime.
The Cartesian flow chart carries an arbitrary flow and so carries a family of spacetimes.
The plane symmetric chart is the member of that family whose flow points one way and is carried by the two coordinates across it.
The collection has done this once before, in `pp_wave`, whose exact plane wave chart is the subcase of its Brinkmann chart with a quadratic profile, and the reason is the same: the general member is where the definition lives and the special member is where the curvature closes.

---

## Step 2. The line element, and what a warp drive is made of

Write flat space in Cartesian coordinates and let it flow.
At each moment the slice of constant $t$ is ordinary Euclidean three space, and on it lives a velocity field

$$\vec{V} = u\,\partial_x + v\,\partial_y + w\,\partial_z,$$

whose three components are arbitrary functions of all four coordinates.
The line element of the drive is the flat one written in coordinates that are dragged along by that flow,

$$ds^2 = -c^2dt^2 + \left(dx - u\,dt\right)^2 + \left(dy - v\,dt\right)^2 + \left(dz - w\,dt\right)^2,$$

which is the first published line element.
This is the 3+1 form with lapse one, flat spatial slices and shift $-\vec{V}$, and it is the only structure a warp drive has.
Nothing has been assumed about $\vec{V}$ yet, and nothing below assumes anything about it until Step 4.

The reading is the one the history field gives.
A curve with $dx^i/dt = V^i$ has $ds^2 = -c^2dt^2$, so an observer carried by the flow ages at the rate the coordinate $t$ counts, whatever the flow is doing.
Those are the Eulerian observers, the ones at rest with respect to the slices, and their four velocity is

$$n^\mu = \left(1, \frac{u}{c}, \frac{v}{c}, \frac{w}{c}\right), \qquad n_\mu = \left(-1, 0, 0, 0\right) = -\partial_\mu(ct),$$

as Step 3 confirms once the metric is written down.
A ship is a point held at rest in the flow, and the whole trick of a warp drive is that $\vec{V}$ has no speed limit while the ship never moves through the space around it.

Alcubierre's choice is a flow with one component,

$$\vec{V} = c\,v_s f\,\partial_x,$$

with $v_s$ the speed of the bubble in units of $c$ and $f$ a shape function that is one at the ship and zero outside the bubble.
His own profile takes $f$ to be a function of the distance $r_s = \sqrt{\left(x - x_s(t)\right)^2 + y^2 + z^2}$ from the ship, and the entry that publishes that drive leaves $f$ an arbitrary function of all four coordinates, which is the form used in Step 10.
Natário's choice is any $\vec{V}$ at all whose divergence vanishes.
Both are this line element; Step 4 is where they part.

---

## Step 3. The metric matrix, its determinant and its inverse

Expanding the squares in the line element and dividing by the factor of $c$ that each $dt$ carries into the chart $x^0 = ct$,

$$g_{tt} = -\left(1 - \frac{u^2 + v^2 + w^2}{c^2}\right), \qquad g_{ti} = -\frac{V_i}{c}, \qquad g_{ij} = \delta_{ij},$$

which are the published components, with $V_i$ standing for $u$, $v$ and $w$ in turn.
Every one of them is dimensionless, as a chart component must be.

The determinant is worth a line of its own.
The matrix has the block form with $\delta_{ij}$ in its spatial corner, so its determinant is the Schur complement $\det(\delta_{ij})\left(g_{tt} - g_{ti}\delta^{ij}g_{jt}\right)$, which is

$$\det g = -\left(1 - \frac{V^2}{c^2}\right) - \frac{V^2}{c^2} = -1,$$

exactly, for every flow.
So $\sqrt{-g} = 1$ and the coordinate four volume element is the flat one.
That is a property of the 3+1 form with unit lapse and flat slices and not of Natário's condition, and the two are easy to confuse.
Natário's condition is about the three volume carried along the flow, which Step 4 computes, not about this determinant.

The inverse is

$$g^{tt} = -1, \qquad g^{ti} = -\frac{V^i}{c}, \qquad g^{ij} = \delta^{ij} - \frac{V^iV^j}{c^2},$$

as published, and it is checked in one line:

$$g^{t\alpha}g_{\alpha t} = (-1)\left(-1 + \frac{V^2}{c^2}\right) + \left(-\frac{V^i}{c}\right)\left(-\frac{V_i}{c}\right) = 1 - \frac{V^2}{c^2} + \frac{V^2}{c^2} = 1.$$

Lowering $n^\mu$ with this metric gives $n_t = g_{tt} + g_{ti}V^i/c = -1 + V^2/c^2 - V^2/c^2 = -1$ and $n_i = g_{it} + g_{ij}V^j/c = -V_i/c + V_i/c = 0$, which is the claim of Step 2.
The Eulerian observers are the ones whose one form is $-d(ct)$, so they are hypersurface orthogonal and the slices are their instantaneous spaces.

---

## Step 4. The expansion, and the one number that separates Natário from Alcubierre

The expansion of a congruence is the divergence of its four velocity, and with $\sqrt{-g} = 1$ from Step 3 that divergence is the flat one:

$$\nabla_\mu n^\mu = \frac{1}{\sqrt{-g}}\partial_\mu\left(\sqrt{-g}\,n^\mu\right) = \partial_t(1) + \partial_i\!\left(\frac{V^i}{c}\right) = \frac{1}{c}\left(\partial_x u + \partial_y v + \partial_z w\right).$$

This is the rate at which a small ball of Eulerian observers changes volume per unit of proper length along their worldlines, and it carries $1/L$ as an expansion must.
Natário's condition is that it vanish,

$$\partial_x u + \partial_y v + \partial_z w = 0,$$

and that is the whole of his construction.
A divergence free field is the velocity field of an incompressible flow, which is why the history can say that space slides around the ship like water around a stone with no parcel compressed on the way past.

Alcubierre's field fails it, and by how much is worth writing down.
With $\vec{V} = c\,v_sf\,\partial_x$,

$$\partial_xu + \partial_yv + \partial_zw = c\,v_s\,\partial_xf,$$

which is what the Alcubierre entry publishes for $\nabla_\mu n^\mu$.
For his own radial profile that is $c\,v_s\,f'(r_s)\left(x - x_s\right)/r_s$, negative in front of the ship and positive behind it for a bump that falls off outward.
That is the contraction ahead and expansion behind of the story, and it is a property of the one particular field, not of the line element the two drives share.

The rate of strain of the flow,

$$\theta_{ij} = \tfrac{1}{2}\left(\partial_iV_j + \partial_jV_i\right), \qquad \theta = \delta^{ij}\theta_{ij} = \partial_xu + \partial_yv + \partial_zw,$$

is the object every tensor below is built from, and Natário's condition is exactly that its trace vanish while its trace free part is left alone.
Nothing forces the trace free part to be small, and Step 9 is the price of that.
The antisymmetric part,

$$\omega_{ij} = \tfrac{1}{2}\left(\partial_iV_j - \partial_jV_i\right),$$

is the vorticity, and it appears in the connection at Step 5 and in the geodesics at Step 11 without ever appearing in the energy density.

---

## Step 5. The Christoffel symbols of a general flow

With $g_{ij} = \delta_{ij}$ constant and $g_{tt}$ and $g_{ti}$ built from the flow, the connection is a short computation whose answer is entirely in $\theta_{ij}$ and $\omega_{ij}$.
Writing $\beta^i = V^i/c$ for the dimensionless flow, the published symbols are

$$\Gamma^t{}_{ij} = \frac{\theta_{ij}}{c}, \qquad \Gamma^t{}_{ti} = -\frac{\beta^j\theta_{ji}}{c}, \qquad \Gamma^t{}_{tt} = \frac{\beta^i\beta^j\theta_{ij}}{c},$$

$$\Gamma^i{}_{jk} = \beta^i\,\Gamma^t{}_{jk}, \qquad \Gamma^i{}_{tj} = \beta^i\,\Gamma^t{}_{tj} + \frac{\omega_{ij}}{c}, \qquad \Gamma^i{}_{tt} = \beta^i\,\Gamma^t{}_{tt} - \frac{1}{c}\left(\partial_tV_i + \frac{V^j\partial_iV_j}{c}\right).$$

Three of these are worth reading rather than checking.
The first says that the whole of the connection among the spatial directions is the rate of strain of the flow, so a flow with no strain at all, a rigid motion of flat space, has no connection and no curvature.
The fourth says that every spatial symbol is the corresponding time symbol carried along by the flow, which is the algebraic shadow of the fact that the slices are being dragged rigidly.
The last carries $\partial_tV_i + \partial_i(V^2/2)/c$, the combination that appears in the Bernoulli equation of an ordinary fluid, and it is what accelerates a particle that is not moving with the flow.

Two published values check the pattern.
For $\Gamma^t{}_{xy}$ the formula gives $\theta_{xy}/c = \left(\partial_yu + \partial_xv\right)/2c$, which is what the file prints.
For $\Gamma^x{}_{yz}$ it gives $\beta^x\theta_{yz}/c = u\left(\partial_zv + \partial_yw\right)/2c^2$, which is also what the file prints.
The file writes all forty of them out, because the checker reads it symbol by symbol and has no way to be told what $\theta_{ij}$ means.

---

## Step 6. The lowered Christoffel symbols

Lowering with $g_{\mu\alpha}$ mixes the time and space rows, since $g_{ti}$ does not vanish, and the result is published in full.
The compact pattern of Step 5 does not survive the lowering: $\Gamma_{\mu\nu\rho} = g_{\mu\alpha}\Gamma^\alpha{}_{\nu\rho}$ puts a factor of $\left(1 - V^2/c^2\right)$ on one term and $-\beta_i$ on another, and the sum does not collect into $\theta_{ij}$ and $\omega_{ij}$ again.
The file prints the components rather than a formula, and what is worth stating is the check they satisfy, which is that they are symmetric in their last two indices and reproduce the metric under

$$\partial_\rho g_{\mu\nu} = \Gamma_{\mu\nu\rho} + \Gamma_{\nu\mu\rho},$$

which sympy confirms component by component in the published file.
For the plane symmetric chart of Step 14 the lowered symbols are short enough to read, and there the identity can be checked by eye: $\Gamma_{xty} + \Gamma_{txy} = -\partial_yu/2c - \partial_yu/2c = -\partial_yu/c$, which is $\partial_y g_{tx} = \partial_y(-u/c)$.

---

## Step 7. The Ricci scalar of a general flow

The Ricci scalar of this family is

$$R = \frac{2}{c}\left(\partial_t\theta + \frac{V^i\partial_i\theta}{c}\right) + \frac{\theta^2 + \theta_{ij}\theta^{ij}}{c^2},$$

which the file prints with $\theta$ written out in its three terms.
The first bracket is the material derivative of the expansion, the rate at which the expansion changes along the flow, and Natário's condition kills it twice over: $\theta$ is zero everywhere, so both its time derivative and its gradient are zero.
What is left on that condition is

$$R = \frac{\theta_{ij}\theta^{ij}}{c^2} \ge 0,$$

the squared rate of strain, which vanishes only where the flow is a rigid motion.
So a Natário drive has a positive Ricci scalar wherever it is doing anything at all.
The sign is the first hint of Step 9: the Ricci scalar of ordinary matter with positive energy density and small pressure is negative on this convention, since $R = -8\pi G T/c^4$ and $T \approx -\rho c^2$.

---

## Step 8. The Einstein tensor, and where its three pieces come from

The entry publishes the Einstein tensor of the general flow with both indices up, which is the variant in which its pieces are the constraint equations of the 3+1 split and can be read off rather than dug out.
Because $n_\mu = -\partial_\mu(ct)$ has only a time component, projecting on $n$ is the same as taking a $t$ index:

$$G^{\mu\nu}n_\mu n_\nu = G^{tt}, \qquad G^{\mu\nu}n_\mu = -G^{t\nu}.$$

The Gauss equation for a slice of intrinsic curvature $^{(3)}R$ and extrinsic curvature $K_{ij}$ gives the Hamiltonian constraint

$$2\,G_{\mu\nu}n^\mu n^\nu = {}^{(3)}R + K^2 - K_{ij}K^{ij}.$$

Here the slices are flat, so $^{(3)}R = 0$, and with unit lapse and flat slices the extrinsic curvature is the rate of strain of the flow, $K_{ij} = -\theta_{ij}/c$, with $K = -\theta/c$.
The constraint is quadratic in $K$, so the sign convention for it does not enter the answer.
What comes out is

$$G^{tt} = \frac{\theta^2 - \theta_{ij}\theta^{ij}}{2c^2},$$

which is the published component, written there as

$$G^{tt} = \frac{\partial_xu\,\partial_yv + \partial_yv\,\partial_zw + \partial_zw\,\partial_xu}{c^2} - \frac{\left(\partial_yu + \partial_xv\right)^2 + \left(\partial_zu + \partial_xw\right)^2 + \left(\partial_zv + \partial_yw\right)^2}{4c^2},$$

the two forms being the same expression multiplied out.

The Codazzi equation gives the momentum constraint, and it comes out as the piece of $G^{ti}$ that is not simply $G^{tt}$ carried along by the flow:

$$G^{ti} = \frac{V^i}{c}\,G^{tt} + \frac{1}{2c}\left(\partial_i\theta - \nabla^2V^i\right),$$

where $\nabla^2$ is the flat Laplacian of the slice.
On Natário's condition the first term inside the bracket is gone and the momentum constraint is the Laplacian of the flow alone.

The spatial components are the evolution equations, and they carry the second time derivatives of the flow:

$$G^{ij} = S^{ij} + \frac{V^i}{c}\,\frac{\partial_j\theta - \nabla^2V^j}{2c} + \frac{V^j}{c}\,\frac{\partial_i\theta - \nabla^2V^i}{2c} + \frac{V^iV^j}{c^2}\,G^{tt},$$

with $S^{ij}$ the part that is genuinely new, printed in the file term by term.
Every one of the ten published components is one of these three objects plus the flow carrying the others along, which is why each of them repeats the expression for $G^{tt}$ inside it.

---

## Step 9. The energy density, and why removing the expansion does not help

The field equations are $G_{\mu\nu} = 8\pi G\,T_{\mu\nu}/c^4$.
The energy density an Eulerian observer measures is

$$\rho = T_{\mu\nu}n^\mu n^\nu = \frac{c^4}{8\pi G}\,G_{\mu\nu}n^\mu n^\nu = \frac{c^4}{8\pi G}\,G^{tt},$$

so Step 8 has already computed it:

$$\rho = \frac{c^2}{16\pi G}\left(\theta^2 - \theta_{ij}\theta^{ij}\right).$$

Read that expression before imposing anything.
It is a difference of two squares, and the Cauchy-Schwarz inequality on a symmetric three by three matrix gives $\theta^2 \le 3\,\theta_{ij}\theta^{ij}$, with equality only when $\theta_{ij}$ is a multiple of $\delta_{ij}$, that is only when the flow is a pure uniform expansion at that point.
So the positive term can never be larger than three times the negative one, and it can reach that only for a flow that is expanding isotropically everywhere, which no bubble is.

Now impose Natário's condition, $\theta = 0$.
The positive term is gone and what is left is

$$\rho = -\frac{c^2}{16\pi G}\,\theta_{ij}\theta^{ij} \le 0,$$

with equality only where $\theta_{ij}$ vanishes entirely, that is only where the flow is a rigid motion of flat space and the drive is doing nothing.
A bubble that carries a ship has a flow that is at rest far away and moving at the ship, so $\theta_{ij}$ cannot vanish everywhere, and the energy density is strictly negative somewhere.
The weak energy condition fails, and with it the dominant.

This is the result the entry exists to show, and the way it arrives matters.
Alcubierre's drive violates the energy conditions too, and a natural first thought is that the expansion is to blame, since the expansion is the visible strangeness of his construction and the thing his story is told about.
The computation says otherwise.
The expansion enters the energy density only through $+\theta^2$, the one term that is positive, so it was never the source of the violation; it was the only thing that could ever have worked against it.
Natário's drive removes it and the negative term is untouched.
Step 10 sharpens this from an inequality into an identity on Alcubierre's own field, where the expansion cancels out of the energy density exactly.
That is what the history means by saying the expansion was scenery.

The null energy condition, which is weaker than the weak one and so could in principle survive it, fails as well.
Step 17 checks that on an explicit null vector in the plane symmetric chart, where the arithmetic closes, and finds the same sign.
Lobo and Visser reached the same conclusion for both drives without any approximation, and Santiago, Schuster and Visser later showed that no warp drive a physicist would call reasonable escapes it.

---

## Step 10. The same energy density, computed on Alcubierre's own drive

The comparison can be made exact, because `MFS/assets/data/metrics/alcubierre.json` publishes the other drive in the same chart and the same conventions.
Its line element is $ds^2 = -c^2dt^2 + \left(dx - v_sf\,c\,dt\right)^2 + dy^2 + dz^2$ with $v_s$ the bubble speed in units of $c$ and $f$ a shape function left arbitrary in all four coordinates.
In the notation used here that is the flow

$$\vec{V} = c\,v_s f\,\partial_x,$$

a flow with one component which does depend on the coordinate it points along.
Its expansion is $\theta = \partial_xu = c\,v_s\,\partial_xf$, which is what that entry publishes for $\nabla_\mu n^\mu$, so the two entries agree on the number Step 4 says they differ by.

Now put that flow into the energy density of Step 9 rather than into the story.
With only one component the rate of strain has

$$\theta_{xx} = \partial_xu, \qquad \theta_{xy} = \frac{\partial_yu}{2}, \qquad \theta_{xz} = \frac{\partial_zu}{2},$$

and nothing else, so $\theta^2 = \left(\partial_xu\right)^2$ while $\theta_{ij}\theta^{ij} = \left(\partial_xu\right)^2 + \left(\left(\partial_yu\right)^2 + \left(\partial_zu\right)^2\right)/2$, and

$$\theta^2 - \theta_{ij}\theta^{ij} = -\frac{\left(\partial_yu\right)^2 + \left(\partial_zu\right)^2}{2}.$$

The expansion has cancelled identically, against the longitudinal part of the rate of strain that always accompanies it.
So

$$\rho = -\frac{c^2}{32\pi G}\left(\left(\partial_yu\right)^2 + \left(\partial_zu\right)^2\right) = -\frac{c^4v_s^2}{32\pi G}\left(\left(\partial_yf\right)^2 + \left(\partial_zf\right)^2\right),$$

which is the energy density the Alcubierre entry publishes, and $G^{tt}$ and the Ricci scalar of Step 7 reproduce that entry's published components term for term when the same substitution is made.
Two entries built independently agree, which is worth more than either of them alone.

Read what that cancellation says.
The expansion never entered the energy density of Alcubierre's drive at all.
It is not that Natário's condition removes a cost; the cost was never charged to the expansion, and the negative energy of the Alcubierre drive is already the transverse gradient of its flow and nothing else.
The two drives have the same expression for the energy density, and the only difference between them is whether the profile is allowed to depend on the coordinate the flow points along.

What Natário's construction does buy is real but smaller than the story suggests.
Alcubierre's drive stretches and squeezes volume elements and Natário's does not, which is a genuine geometric difference and the one the two fields actually differ by.
Everything the drive needs in order to be a drive survives the change: the ship still rides in a flat region, still moves at $\vec{V}$ with respect to the outside, and still has no speed limit.
Everything the drive needs in order to exist fails the same way in both.
The negative energy is not reduced, not rearranged into a more comfortable place and not made to depend on the bubble being slow: it is the trace free rate of strain squared, and a flow that is at rest at infinity and moving at the ship must be strained in between.

---

## Step 11. The geodesic equations of a general flow

The geodesic equation is

$$\ddot{x}^\mu + \Gamma^\mu{}_{\nu\rho}\dot{x}^\nu\dot{x}^\rho = 0,$$

with the dots the chart velocities of Step 1.
Substituting the connection of Step 5 and writing

$$A^i = \dot{x}^i - \frac{V^i}{c}\dot{t}$$

for the velocity of the particle relative to the flow it is sitting in, the time equation collapses to

$$\ddot{t} + \frac{1}{c}\,\theta_{ij}A^iA^j = 0,$$

which is the first published equation with $\theta_{ij}$ written out in its six terms.
The spatial equations are

$$\ddot{x}^i + \frac{V^i}{c^2}\,\theta_{jk}A^jA^k - \frac{\dot{t}}{c}\,\left(\partial_jV_i - \partial_iV_j\right)A^j - \frac{\dot{t}^2}{c}\left(\partial_tV_i + \frac{V^j\partial_jV_i}{c}\right) = 0,$$

which are the other three, again with everything written out.

Three readings, in order of importance.
The last bracket is the material derivative of the flow, $\partial_tV_i + \left(\vec{V}\cdot\nabla\right)V_i/c$, the acceleration of the flow itself, and it is the only term that survives when the particle moves with the flow.
The middle term is twice the vorticity of Step 4 acting on the relative velocity, a Coriolis force, and being antisymmetric it does no work on that relative velocity.
Every other term carries $A^i$ at least once, so a particle with $A^i = 0$ everywhere along its worldline, one that is exactly carried by the flow, has

$$\ddot{t} = 0, \qquad \ddot{x}^i = \frac{\dot{t}^2}{c}\left(\partial_tV_i + \frac{V^j\partial_jV_i}{c}\right),$$

and the second of these is just $d(V^i\dot{t}/c)/d\lambda$, the statement that the particle keeps moving with the flow.
So the Eulerian worldlines are geodesics.
That is the sense in which a ship at rest in the drive feels nothing, and it is checked here on the published equations rather than asserted.

---

## Step 12. The plane symmetric member, and why it is divergence free of itself

The curvature of a general flow does not fit on a page.
With three free functions of four coordinates a single component of the Riemann tensor runs to dozens of terms and the Kretschmann scalar to thousands, so the entry publishes the curvature in the member of the family where it closes.

Take the flow to point along $x$ and to be carried by the two coordinates across it,

$$\vec{V} = u(t,y,z)\,\partial_x.$$

Then

$$\partial_xu + \partial_yv + \partial_zw = \partial_xu = 0$$

identically, because $u$ does not depend on $x$, and the other two components are not there.
So this is a Natário drive for every profile $u$ whatever, with no condition left to impose, and every number published in this chart is a number of Natário's drive rather than of a general flow.
The line element is

$$ds^2 = -c^2dt^2 + \left(dx - u\,dt\right)^2 + dy^2 + dz^2.$$

What it describes is a warp corridor.
Where $u$ is constant every Christoffel symbol of Step 14 vanishes, since each carries a derivative of $u$, so the geometry there is flat; and it is flat in the obvious way, because $x' = x - ut$ with constant $u$ turns the line element into $-c^2dt^2 + dx'^2 + dy^2 + dz^2$ exactly.
Take $u$ to be a constant $v_s$ inside a tube of any cross section, zero outside it, and smooth in between.
Inside, the ship sits at rest in a flat region and moves at $v_s$ with respect to the outside; outside, the universe is at rest; and the wall between them carries the whole of the curvature.
The map between the two flat regions is $x' = x - v_st$, which is a Galilean shear and not a Lorentz boost, and that is exactly why $v_s$ may exceed $c$.
They are two patches joined through a wall, not two inertial frames of one flat spacetime.

What this flow cannot do is close.
A flow along $x$ that is divergence free cannot depend on $x$, so the corridor runs to infinity in both directions and there is no front and no back.
Closing it off means letting $u$ depend on $x$, and then the divergence forces transverse components, which is precisely the structure of Natário's own field: the flow parts ahead of the ship, slides around it and closes behind.
The Cartesian chart of Steps 2 to 11 is where that general case lives.

This chart is also exactly where the two drives meet.
Step 10 wrote the Alcubierre drive as the one component flow $\vec{V} = c\,v_sf\,\partial_x$, and Natário's condition on a flow with one component is the single requirement that it not depend on $x$, which for that entry reads $\partial_xf = 0$.
So the plane symmetric chart is the divergence free member of the Alcubierre family, published there with an arbitrary shape function and here with the one condition imposed.
The two entries carry the same line element with the same conventions, and everything below can be read against the corresponding step of `_tools/derivations/alcubierre.md`.

---

## Step 13. The plane symmetric metric and its inverse

Expanding the line element in the chart $x^0 = ct$,

$$g_{tt} = -\left(1 - \frac{u^2}{c^2}\right), \qquad g_{tx} = g_{xt} = -\frac{u}{c}, \qquad g_{xx} = g_{yy} = g_{zz} = 1,$$

with every other component zero, as published.
The determinant is $-1$ by the computation of Step 3, and the inverse is

$$g^{tt} = -1, \qquad g^{tx} = g^{xt} = -\frac{u}{c}, \qquad g^{xx} = 1 - \frac{u^2}{c^2}, \qquad g^{yy} = g^{zz} = 1.$$

The single check is $g^{tt}g_{tt} + g^{tx}g_{xt} = \left(1 - u^2/c^2\right) + u^2/c^2 = 1$.

---

## Step 14. The plane symmetric Christoffel symbols

Only $\partial_tu$, $\partial_yu$ and $\partial_zu$ can appear, since $u$ does not depend on $x$.
The published symbols are

$$\Gamma^t{}_{ty} = -\frac{u\,\partial_yu}{2c^2}, \qquad \Gamma^t{}_{xy} = \frac{\partial_yu}{2c}, \qquad \Gamma^x{}_{tt} = -\frac{\partial_tu}{c}, \qquad \Gamma^x{}_{ty} = -\frac{\partial_yu}{2c}\left(1 + \frac{u^2}{c^2}\right),$$

$$\Gamma^x{}_{xy} = \frac{u\,\partial_yu}{2c^2}, \qquad \Gamma^y{}_{tt} = -\frac{u\,\partial_yu}{c^2}, \qquad \Gamma^y{}_{tx} = \frac{\partial_yu}{2c},$$

together with everything that follows from swapping $y$ for $z$ and from the symmetry in the last two indices.
They agree with the general pattern of Step 5 on the one nonzero rate of strain and vorticity this flow has,

$$\theta_{xy} = \frac{\partial_yu}{2}, \qquad \theta_{xz} = \frac{\partial_zu}{2}, \qquad \omega_{yx} = \frac{\partial_yu}{2}, \qquad \omega_{zx} = \frac{\partial_zu}{2},$$

with $\theta = 0$ as Step 12 requires.
For instance Step 5 gives $\Gamma^y{}_{tx} = \beta^y\Gamma^t{}_{tx} + \omega_{yx}/c$, and with $\beta^y = 0$ that is $\partial_yu/2c$, which is the printed value.

Every symbol carries at least one derivative of $u$, which is the statement that a corridor of constant speed is flat.

---

## Step 15. The plane symmetric Riemann tensor

The independent components of the lowered Riemann tensor, taking the pair symmetries into account, are

$$R_{txtx} = \frac{\left(\partial_yu\right)^2 + \left(\partial_zu\right)^2}{4c^2}, \qquad R_{txty} = -\frac{\partial_t\partial_yu}{2c}, \qquad R_{txtz} = -\frac{\partial_t\partial_zu}{2c},$$

$$R_{tyty} = -\frac{u\,\partial_y^2u}{c^2} - \frac{3\left(\partial_yu\right)^2}{4c^2} - \frac{u^2\left(\partial_yu\right)^2}{4c^4}, \qquad R_{tytz} = -\frac{u\,\partial_y\partial_zu}{c^2} - \frac{3\,\partial_yu\,\partial_zu}{4c^2} - \frac{u^2\,\partial_yu\,\partial_zu}{4c^4},$$

$$R_{tyxy} = \frac{\partial_y^2u}{2c} + \frac{u\left(\partial_yu\right)^2}{4c^3}, \qquad R_{tyxz} = R_{tzxy} = \frac{\partial_y\partial_zu}{2c} + \frac{u\,\partial_yu\,\partial_zu}{4c^3},$$

$$R_{xyxy} = -\frac{\left(\partial_yu\right)^2}{4c^2}, \qquad R_{xyxz} = -\frac{\partial_yu\,\partial_zu}{4c^2}, \qquad R_{xzxz} = -\frac{\left(\partial_zu\right)^2}{4c^2},$$

together with the ones obtained by swapping $y$ for $z$.
The file publishes the mixed form $R^\mu{}_{\nu\rho\sigma}$ in a hundred slots and the lowered form in as many, since each independent value fills several of them through the antisymmetries.

Two of these are worth a remark.
The purely spatial block is the Gauss equation with nothing else in it.
The slice is flat, so $R_{ijkl} = K_{ik}K_{jl} - K_{il}K_{jk}$, and with the only nonzero extrinsic curvatures $K_{xy} = -\partial_yu/2c$ and $K_{xz} = -\partial_zu/2c$ from Step 8 that gives

$$R_{xyxy} = K_{xx}K_{yy} - K_{xy}K_{xy} = -\frac{\left(\partial_yu\right)^2}{4c^2},$$

which is the value above.
The curvature of the spatial block is therefore the shear of the flow squared, with a minus sign, and it is a statement about how the slice sits in spacetime rather than about the slice, which is flat.
And $R_{yzyz}$ is absent from the list because it vanishes, which Step 19 needs.

---

## Step 16. The plane symmetric Ricci tensor and Ricci scalar

Contracting on the first index, $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$, gives the published values, of which the three that carry the argument are

$$R_{xx} = -\frac{\left(\partial_yu\right)^2 + \left(\partial_zu\right)^2}{2c^2}, \qquad R_{yy} = \frac{\left(\partial_yu\right)^2}{2c^2}, \qquad R_{yz} = \frac{\partial_yu\,\partial_zu}{2c^2},$$

with $R_{tx}$, $R_{tt}$ and the mixed time components following from these by the flow carrying them along.
The Ricci scalar is

$$R = g^{\mu\nu}R_{\mu\nu} = \frac{\left(\partial_yu\right)^2 + \left(\partial_zu\right)^2}{2c^2},$$

which is the published value, and it is exactly the general result of Step 7 on Natário's condition:

$$\theta_{ij}\theta^{ij} = 2\left(\theta_{xy}^2 + \theta_{xz}^2\right) = \frac{\left(\partial_yu\right)^2 + \left(\partial_zu\right)^2}{2},$$

and $R = \theta_{ij}\theta^{ij}/c^2$ as Step 7 says it must be.
The general formula and the special chart agree, which is the one place the two published systems can be checked against each other.

---

## Step 17. The plane symmetric Einstein tensor, the matter and the energy conditions

The Einstein tensor is published in all three variants here, since the components are short.
With both indices up the components are

$$G^{tt} = -\frac{\left(\partial_yu\right)^2 + \left(\partial_zu\right)^2}{4c^2}, \qquad G^{xy} = \frac{\partial_t\partial_yu}{2c}, \qquad G^{yy} = -G^{zz} = \frac{\left(\partial_yu\right)^2 - \left(\partial_zu\right)^2}{4c^2}, \qquad G^{yz} = \frac{\partial_yu\,\partial_zu}{2c^2},$$

with $G^{tx}$ and $G^{xx}$ carrying these along with the flow as Step 8 describes.

The energy density is immediate:

$$\rho = \frac{c^4}{8\pi G}G^{tt} = -\frac{c^2}{32\pi G}\left(\left(\partial_yu\right)^2 + \left(\partial_zu\right)^2\right),$$

negative wherever the profile has any transverse gradient at all, which is everywhere in the wall of the corridor and nowhere else.
It agrees with the general formula of Step 9, $\rho = -c^2\theta_{ij}\theta^{ij}/16\pi G$, on the value of $\theta_{ij}\theta^{ij}$ computed in Step 16.
The expansion is exactly zero here, by construction and not by approximation, and the energy density is still negative.
That is the entry's claim, stated on an exact solution with no small parameter anywhere in it.

The null energy condition fails too, and on this chart the failure can be exhibited rather than cited.
Take the null vector $k^\mu = n^\mu + e^\mu$ with $n^\mu$ the Eulerian four velocity of Step 3 and $e^\mu = \partial_y$, which is a unit spatial vector orthogonal to it.
Then $G_{\mu\nu}n^\mu e^\nu = 0$, because $G_{ty}n^t + G_{xy}n^x = -u\,\partial_t\partial_yu/2c^2 + \left(\partial_t\partial_yu/2c\right)\left(u/c\right)$ cancels exactly, and so

$$G_{\mu\nu}k^\mu k^\nu = G^{tt} + G_{yy} = -\frac{\left(\partial_yu\right)^2 + \left(\partial_zu\right)^2}{4c^2} + \frac{\left(\partial_yu\right)^2 - \left(\partial_zu\right)^2}{4c^2} = -\frac{\left(\partial_zu\right)^2}{2c^2},$$

which is negative unless $\partial_zu$ vanishes.
Taking $e^\mu = \partial_z$ instead gives $-\left(\partial_yu\right)^2/2c^2$ by the same cancellation, so unless the profile is constant one of the two is strictly negative and the null energy condition is violated.
Since the null condition is the weakest of the four, its failure carries the weak, the strong and the dominant with it.

---

## Step 18. The plane symmetric Kretschmann scalar

The Kretschmann scalar is

$$K = R_{\mu\nu\rho\sigma}R^{\mu\nu\rho\sigma} = \frac{11\left(\left(\partial_yu\right)^2 + \left(\partial_zu\right)^2\right)^2}{4c^4} - \frac{2\left(\left(\partial_y^2u\right)^2 + 2\left(\partial_y\partial_zu\right)^2 + \left(\partial_z^2u\right)^2\right)}{c^2} + \frac{2\left(\left(\partial_t\partial_yu\right)^2 + \left(\partial_t\partial_zu\right)^2\right)}{c^2},$$

which is the published value, and every profile dependence in it is transverse.
It is free of $u$ itself, as a curvature invariant of this geometry has to be, since adding a constant to $u$ is the coordinate change $x' = x - \text{constant}\times t$ of Step 12 and cannot move an invariant.

Three things follow.
Where the profile is flat, inside the corridor and outside it, $K$ vanishes and the geometry is flat, which is the statement Step 12 made with the Christoffel symbols.
The scalar is not positive definite, since the middle term can dominate: a profile with a large transverse Hessian and a small gradient has $K < 0$.
And nothing in it blows up for any smooth profile, so the corridor has no singularity anywhere.
The obstructions to building one are the negative energy of Step 17 and the horizons the history describes, not a curvature that runs away.

---

## Step 19. The plane symmetric Weyl tensor, computed and not copied

The Weyl tensor is Riemann with its traces removed,

$$C_{\mu\nu\rho\sigma} = R_{\mu\nu\rho\sigma} - \frac{1}{2}\left(g_{\mu\rho}R_{\sigma\nu} - g_{\mu\sigma}R_{\rho\nu} - g_{\nu\rho}R_{\sigma\mu} + g_{\nu\sigma}R_{\rho\mu}\right) + \frac{R}{6}\left(g_{\mu\rho}g_{\sigma\nu} - g_{\mu\sigma}g_{\rho\nu}\right)$$

in four dimensions, with the traces built from the same contraction the entry publishes.
This spacetime is not a vacuum, so the correction terms do not vanish and the Weyl tensor is not a copy of the Riemann tensor.
`_tools/derivations/weyl.md` is the collection's working on exactly this point and it is what the computation here follows.

The difference is visible in slots where the two disagree about whether anything is happening at all.
Step 15 recorded that $R_{yzyz} = 0$.
The Weyl tensor in the same slot is

$$C_{yzyz} = -\frac{\left(\partial_yu\right)^2 + \left(\partial_zu\right)^2}{6c^2} = -\frac{R}{3},$$

which is nonzero for every profile that does anything.
So is $C_{xyyz} = \partial_t\partial_zu/4c$, whose Riemann counterpart $R_{xyyz}$ is also zero.
Copying Riemann into the Weyl block would have published zero in both places and would have been wrong in both.

The purely spatial block of Weyl,

$$C_{xyxy} = \frac{\left(\partial_zu\right)^2}{3c^2} - \frac{\left(\partial_yu\right)^2}{6c^2}, \qquad C_{xzxz} = \frac{\left(\partial_yu\right)^2}{3c^2} - \frac{\left(\partial_zu\right)^2}{6c^2}, \qquad C_{yzyz} = -\frac{\left(\partial_yu\right)^2 + \left(\partial_zu\right)^2}{6c^2},$$

sums to zero, which is the trace free condition doing its work and a useful check on the arithmetic.

---

## Step 20. The plane symmetric geodesic equations

From the connection of Step 14 the four equations are

$$\ddot{t} + \frac{\partial_yu}{c}\left(\dot{x} - \frac{u}{c}\dot{t}\right)\dot{y} + \frac{\partial_zu}{c}\left(\dot{x} - \frac{u}{c}\dot{t}\right)\dot{z} = 0,$$

$$\ddot{x} - \frac{\partial_tu}{c}\dot{t}^2 - \frac{\partial_yu}{c}\left(\dot{t} - \frac{u}{c}\dot{x} + \frac{u^2}{c^2}\dot{t}\right)\dot{y} - \frac{\partial_zu}{c}\left(\dot{t} - \frac{u}{c}\dot{x} + \frac{u^2}{c^2}\dot{t}\right)\dot{z} = 0,$$

$$\ddot{y} + \frac{\partial_yu}{c}\left(\dot{x} - \frac{u}{c}\dot{t}\right)\dot{t} = 0, \qquad \ddot{z} + \frac{\partial_zu}{c}\left(\dot{x} - \frac{u}{c}\dot{t}\right)\dot{t} = 0,$$

which are the published equations.
Every term carries the combination $\dot{x} - u\dot{t}/c$, the velocity relative to the flow, except the one term in $\partial_tu$, which is the flow changing under the particle.

So a particle that starts out moving with the flow, $\dot{x} = u\dot{t}/c$ with $\dot{y} = \dot{z} = 0$, has $\ddot{t} = \ddot{y} = \ddot{z} = 0$ and $\ddot{x} = \partial_tu\,\dot{t}^2/c$, which is exactly the rate at which $u\dot{t}/c$ itself changes along the worldline.
It keeps moving with the flow, as Step 11 found in general.
A ship at rest in the corridor is in free fall, is carried at $u$ with respect to the outside universe, and reads its own proper time as the coordinate $t$.

The transverse equations say the rest.
A particle crossing the wall with any velocity relative to the flow is pushed sideways at a rate set by the transverse gradient of the profile, which is the same gradient that carries the negative energy density of Step 17 and the whole of the curvature of Step 18.
Everything this spacetime does, it does in the wall.

---

## Step 21. Every published equation is dimensionally consistent

The declarations are $[t] = T$, $[x] = [y] = [z] = L$ and $[u] = [v] = [w] = L/T$, which is `DIMENSIONS` in `verify_metrics.py` for both systems.
Since the chart is $x^0 = ct$, a chart component with $p$ upper time indices and $q$ lower ones carries the bare component times $c^{p-q}$, and a printed $\partial_t$ is $c^{-1}\partial/\partial t$.

The line element: $c^2dt^2$ carries $L^2$ and so does $\left(dx - u\,dt\right)^2$, since $u\,dt$ carries $(L/T)(T) = L$.
The metric: $u/c$ and $u^2/c^2$ are dimensionless, as every chart component of a metric must be.
The connection: $\left[\partial_yu/c\right] = (1/T)(T/L) = 1/L$, which is what a Christoffel symbol carries, and $\left[u\,\partial_yu/c^2\right] = (L/T)(1/T)(T^2/L^2) = 1/L$ as well.
The curvature: $\left[\left(\partial_yu\right)^2/c^2\right] = (1/T^2)(T^2/L^2) = 1/L^2$ for Riemann, Ricci, Einstein and the Ricci scalar alike, and $\left[\partial_y^2u/c\right] = \left((L/T)/L^2\right)(T/L) = 1/L^2$ for the terms carrying second derivatives.
The Kretschmann scalar: each of its three groups carries $1/L^4$, the first as $\left(1/L^2\right)^2$ and the other two as $\left(1/L^2\right)\left(1/L^2\right)$ through $\left[\partial_y^2u\right] = 1/(LT)$ and $\left[\partial_t\partial_yu\right] = 1/(LT)$.

The geodesics are the case where the chart convention earns its keep.
Every term must carry $[x^\mu]$ over an affine parameter squared, and with $\dot{t} = d(ct)/d\lambda$ carrying $L/\lambda$ exactly as $\dot{x}$ does, the term $\left(\partial_yu/c\right)\left(\dot{x} - u\dot{t}/c\right)\dot{y}$ carries $(1/L)(L/\lambda)(L/\lambda) = L/\lambda^2$, which matches $\ddot{t}$.
On the other reading, with $\dot{t} = dt/d\lambda$, the same term would be short by a factor of $c$ against $\ddot{t}$, and the mixed terms in the $x$ equation would each be wrong by one factor in the opposite direction.
The dimensional pass of `verify_metrics.py` makes this check over every term of every published expression and reports nothing for either system.

---

## Step 22. What the entry publishes, and what the checker did

The Cartesian flow chart publishes the line element, the metric and its inverse, both Christoffel variants, the Ricci scalar, the Einstein tensor with both indices up and the four geodesic equations, for an arbitrary flow.
It stops there on purpose.
Its Riemann tensor has a hundred and ninety two nonvanishing components whose longest runs to fifty eight terms, and its Kretschmann scalar is longer than the rest of the entry put together; publishing them would be a wall of algebra and not a piece of knowledge, and Step 12 says where to find the curvature instead.

The plane symmetric chart publishes the full set: line element, metric, inverse, both Christoffel variants, both Riemann variants, all three Ricci variants, the Ricci scalar, the Kretschmann scalar, all three Einstein variants, both Weyl variants and the geodesics.

`verify_metrics.py` rebuilds all of it from the two line elements and compares it symbol by symbol, with a component the file omits required to vanish, and reports no disagreement and nothing unchecked for either system.
The Cartesian chart takes about ninety seconds, most of it in the Riemann tensor that the Einstein tensor is built from and that the file does not print; the plane symmetric chart takes about twenty.
The dimensional pass reports nothing for either, which is the check of Step 21 run over every term.
