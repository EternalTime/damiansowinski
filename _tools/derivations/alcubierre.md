# The Alcubierre warp drive

This is the working behind the coordinate system in `MFS/assets/data/metrics/alcubierre.json`.
Every value the entry prints is derived here, in order, from the line element down to the geodesic equations.
Nothing is left as an exercise and nothing is asserted that is not computed.

The companion script `_tools/derivations/verify_metrics.py` does the same work in sympy and compares it against the published file, so the algebra below is checkable by hand and by machine independently.

The entry exists for one number.
The history says that the warp drive needs exotic matter and that the amount of it is the objection, and Steps 10 and 11 are where the mathematics says the same thing: the energy density an Eulerian observer measures is

$$\rho = -\frac{c^4}{32\pi G}v_s^2\left[\left(\partial_y f\right)^2 + \left(\partial_z f\right)^2\right],$$

minus a sum of squares, negative wherever the wall of the bubble has any transverse structure at all, growing as the square of the speed and as the square of the steepness of the wall.
Step 11 derives the same expression a second way, from the Hamiltonian constraint, where the negative sign turns out to be structural rather than an accident of the profile.

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
The settlement is not free here.
This spacetime is not a vacuum anywhere the bubble wall is, so the Ricci tensor, the Ricci scalar and the Einstein tensor all change sign with the choice, and the sign the collection now carries is the one that puts ordinary matter at positive energy density.
That matters more for this entry than for any other, because the whole question it answers is the sign of an energy density, and the answer would read backwards under the other contraction.

The chart is the collection's, the one whose zeroth coordinate is

$$x^0 = ct.$$

The index is printed with the bare letter $t$, but the component printed against it is a component of that chart, and everything from Step 4 onward is computed directly in it.
A printed $\partial_t$ is the derivative along that chart coordinate,

$$\partial_t = \frac{1}{c}\frac{\partial}{\partial t},$$

and the prime on $v_s$ is the same derivative, $v_s' = dv_s/d(ct)$.
Each printed derivative therefore carries one inverse length per order, which is what Step 17 checks.
The dots in the geodesic equations are velocities of that same chart, so $\dot{t}$ means $d(ct)/d\lambda$.

---

## Step 2. The line element, and the two things it is built from

The warp drive metric is

$$ds^2 = -c^2dt^2 + \left(dx - v_s f\,c\,dt\right)^2 + dy^2 + dz^2.$$

It carries two functions and no constants.
The first is $v_s(t)$, the velocity of the bubble in units of $c$, so that the centre of the bubble follows the world line $x = x_s(t)$ with $dx_s/dt = c\,v_s$.
Nothing in the solution bounds it; $v_s > 1$ is allowed, and that is the entire point of the construction.
The second is $f$, the shape function, which is $1$ at the centre of the bubble and $0$ far outside it, so that the wall of the bubble is wherever $f$ falls between the two.

The entry leaves $f$ an arbitrary function of all four coordinates.
Alcubierre's own choice is the radial profile

$$f(r_s) = \frac{\tanh\left(\sigma(r_s + R)\right) - \tanh\left(\sigma(r_s - R)\right)}{2\tanh(\sigma R)},\qquad r_s = \sqrt{\left(x - x_s(t)\right)^2 + y^2 + z^2},$$

a bubble of radius $R$ whose wall has thickness of order $1/\sigma$, but no component published here needs it.
A profile of that kind depends on $t$ and $x$ only through $x - x_s(t)$ and so obeys

$$\partial_t f = -v_s\,\partial_x f$$

in this chart, and that is exactly the one relation the published components may not use, because they are identities for an arbitrary $f$ and the checker tests them as such.
Keeping $f$ general costs nothing and buys the whole family: every warp drive of this form, Alcubierre's profile among them, reads its curvature off the same table.

---

## Step 3. The $3+1$ reading, and the shift function

Written as a $3+1$ split in the chart of Step 1,

$$ds^2 = -N^2\left(dx^0\right)^2 + \gamma_{ij}\left(dx^i - \beta^i dx^0\right)\left(dx^j - \beta^j dx^0\right),$$

the line element of Step 2 has

$$N = 1,\qquad \gamma_{ij} = \delta_{ij},\qquad \beta^i = \left(v_s f, 0, 0\right).$$

The shift is dimensionless in this chart, being measured per unit of $x^0$; as a velocity it is $c\,v_s f$.

The lapse is one, the slices of constant $t$ are ordinary flat Euclidean space, and the whole of the geometry sits in a single component of the shift.
That is the sentence the rest of the document unpacks.
Space is not curved on any slice; what the warp drive does is slide the slices past each other, and the curvature is entirely in how they are stacked.

Because the metric depends on $v_s$ and $f$ only through their product, it is worth naming that product.
Write

$$\beta = v_s f,$$

dimensionless, the shift per unit of $x^0$.
Every curvature quantity below is a function of $\beta$ and its derivatives alone, so the computation is done once, in $\beta$, and the published components are what it becomes under

$$\partial_i\beta = v_s\,\partial_i f \quad (i = x, y, z),\qquad \partial_t\beta = v_s'f + v_s\,\partial_t f.$$

The second of those is where $v_s'$ enters the geometry, and the only place it does: the bubble's acceleration reaches the curvature only through $\partial_t\beta$.
This $\beta$ is shorthand for this document alone.
The published file writes every component out in $v_s$ and $f$, because the checker reads the file symbol by symbol and has no way to be told what an abbreviation means.

---

## Step 4. The metric matrix, its determinant and its inverse

Expanding the square in Step 2 and dividing the time terms by $c^2$ to reach the chart $x^0 = ct$,

$$ds^2 = \left(-1 + \beta^2\right)\left(dx^0\right)^2 - 2\beta\,dx^0dx + dx^2 + dy^2 + dz^2,$$

so in the order $(t, x, y, z)$

$$g_{\mu\nu} = \begin{pmatrix} -1 + \beta^2 & -\beta & 0 & 0 \\ -\beta & 1 & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1\end{pmatrix}.$$

The lower right block is the identity and the upper left block has determinant $\left(-1+\beta^2\right) - \beta^2 = -1$, so

$$\det g = -1$$

everywhere, whatever the bubble does.
The chart is therefore never degenerate and carries no coordinate singularity at all: these four coordinates cover the whole spacetime, and there is no second chart to be forced into by a place the first one fails.

Inverting the upper left block by hand,

$$g^{\mu\nu} = \begin{pmatrix} -1 & -\beta & 0 & 0 \\ -\beta & 1 - \beta^2 & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1\end{pmatrix},$$

which the entry publishes as $g^{tt} = -1$, $g^{tx} = g^{xt} = -v_s f$ and $g^{xx} = 1 - v_s^2f^2$.
Note that $g^{tt} = -1$ exactly: the lapse is one, so proper time and coordinate time agree for the observer at rest in the slicing, and that observer is the one Step 10 is about.

---

## Step 5. The Christoffel symbols

With $\Gamma_{\mu\nu\rho} = \tfrac{1}{2}\left(\partial_\nu g_{\mu\rho} + \partial_\rho g_{\mu\nu} - \partial_\mu g_{\nu\rho}\right)$ and the matrix of Step 4, each symbol is a single derivative of $\beta$, sometimes multiplied by $\beta$ itself.
Four of them worked out:

$$\Gamma_{ttt} = \tfrac{1}{2}\partial_t g_{tt} = \tfrac{1}{2}\partial_t\left(\beta^2\right) = \beta\,\partial_t\beta,$$

$$\Gamma_{ttx} = \tfrac{1}{2}\partial_x g_{tt} = \beta\,\partial_x\beta,$$

$$\Gamma_{txx} = \partial_x g_{tx} - \tfrac{1}{2}\partial_t g_{xx} = -\partial_x\beta,$$

$$\Gamma_{xtt} = \partial_t g_{xt} - \tfrac{1}{2}\partial_x g_{tt} = -\partial_t\beta - \beta\,\partial_x\beta.$$

The complete list of the independent lowered symbols is

$$\Gamma_{ttt} = \beta\,\partial_t\beta,\qquad \Gamma_{tti} = \beta\,\partial_i\beta,$$
$$\Gamma_{txx} = -\partial_x\beta,\qquad \Gamma_{txy} = -\tfrac{1}{2}\partial_y\beta,\qquad \Gamma_{txz} = -\tfrac{1}{2}\partial_z\beta,$$
$$\Gamma_{xtt} = -\partial_t\beta - \beta\,\partial_x\beta,\qquad \Gamma_{xty} = -\tfrac{1}{2}\partial_y\beta,\qquad \Gamma_{xtz} = -\tfrac{1}{2}\partial_z\beta,$$
$$\Gamma_{ytt} = -\beta\,\partial_y\beta,\qquad \Gamma_{ytx} = \tfrac{1}{2}\partial_y\beta,\qquad \Gamma_{ztt} = -\beta\,\partial_z\beta,\qquad \Gamma_{ztx} = \tfrac{1}{2}\partial_z\beta,$$

with $i$ running over $x$, $y$, $z$ in the second of them, and with the symmetry in the last two indices supplying the rest.
Raising the first index with $g^{\mu\alpha}$ of Step 4 mixes the $t$ and $x$ rows, since $g^{tx} = -\beta$ is not zero.
Two examples:

$$\Gamma^t{}_{tt} = g^{tt}\Gamma_{ttt} + g^{tx}\Gamma_{xtt} = -\beta\,\partial_t\beta + \left(-\beta\right)\left(-\partial_t\beta - \beta\,\partial_x\beta\right) = \beta^2\,\partial_x\beta,$$

$$\Gamma^x{}_{tt} = g^{xt}\Gamma_{ttt} + g^{xx}\Gamma_{xtt} = -\beta^2\partial_t\beta + \left(1 - \beta^2\right)\left(-\partial_t\beta - \beta\,\partial_x\beta\right) = -\partial_t\beta - \beta\,\partial_x\beta + \beta^3\,\partial_x\beta.$$

The first of those is worth a second look.
$\Gamma^t{}_{tt} = \beta^2\partial_x\beta$ says that coordinate time is not an affine parameter along the world lines of constant $x$, $y$, $z$ once the bubble is moving, and it is the reason the published $\Gamma^t{}_{tt} = v_s^3f^2\partial_x f$ carries three powers of the speed.

The whole upper index list, in $\beta$:

$$\Gamma^t{}_{tt} = \beta^2\partial_x\beta,\qquad \Gamma^t{}_{tx} = -\beta\,\partial_x\beta,\qquad \Gamma^t{}_{ty} = -\tfrac{1}{2}\beta\,\partial_y\beta,\qquad \Gamma^t{}_{tz} = -\tfrac{1}{2}\beta\,\partial_z\beta,$$
$$\Gamma^t{}_{xx} = \partial_x\beta,\qquad \Gamma^t{}_{xy} = \tfrac{1}{2}\partial_y\beta,\qquad \Gamma^t{}_{xz} = \tfrac{1}{2}\partial_z\beta,$$
$$\Gamma^x{}_{tt} = -\partial_t\beta - \beta\,\partial_x\beta + \beta^3\partial_x\beta,\qquad \Gamma^x{}_{tx} = -\beta^2\partial_x\beta,$$
$$\Gamma^x{}_{ty} = -\tfrac{1}{2}\left(1 + \beta^2\right)\partial_y\beta,\qquad \Gamma^x{}_{tz} = -\tfrac{1}{2}\left(1 + \beta^2\right)\partial_z\beta,$$
$$\Gamma^x{}_{xx} = \beta\,\partial_x\beta,\qquad \Gamma^x{}_{xy} = \tfrac{1}{2}\beta\,\partial_y\beta,\qquad \Gamma^x{}_{xz} = \tfrac{1}{2}\beta\,\partial_z\beta,$$
$$\Gamma^y{}_{tt} = -\beta\,\partial_y\beta,\qquad \Gamma^y{}_{tx} = \tfrac{1}{2}\partial_y\beta,\qquad \Gamma^z{}_{tt} = -\beta\,\partial_z\beta,\qquad \Gamma^z{}_{tx} = \tfrac{1}{2}\partial_z\beta.$$

Substituting $\beta = v_s f$ turns these eighteen into the thirty the entry publishes with an upper index, since a term in $\partial_t\beta$ splits into one in $v_s'$ and one in $\partial_t f$.
Every symbol vanishes where $f$ is constant, which is everywhere except the wall: the interior of the bubble and the space outside it are flat, and the ship floats in a patch of Minkowski space.

---

## Step 6. The Riemann tensor

With the connection in hand the curvature is the definition of Step 1 applied component by component.
Two worked examples, both in $\beta$.

The first is an off diagonal component with one upper index:

$$R^t{}_{xxy} = \partial_x\Gamma^t{}_{xy} - \partial_y\Gamma^t{}_{xx} + \Gamma^t{}_{x\lambda}\Gamma^\lambda{}_{xy} - \Gamma^t{}_{y\lambda}\Gamma^\lambda{}_{xx}.$$

The derivative terms give $\tfrac{1}{2}\partial_x\partial_y\beta - \partial_x\partial_y\beta = -\tfrac{1}{2}\partial_x\partial_y\beta$.
In the first quadratic term only $\lambda = t$ and $\lambda = x$ survive, and they cancel:

$$\left(-\beta\,\partial_x\beta\right)\left(\tfrac{1}{2}\partial_y\beta\right) + \left(\partial_x\beta\right)\left(\tfrac{1}{2}\beta\,\partial_y\beta\right) = 0.$$

The second cancels the same way, so

$$R^t{}_{xxy} = -\tfrac{1}{2}\partial_x\partial_y\beta,$$

published as $-\tfrac{1}{2}v_s\,\partial_x\partial_y f$.

The second example is the transverse component, which is the one the energy density comes from:

$$R^x{}_{yxy} = \partial_x\Gamma^x{}_{yy} - \partial_y\Gamma^x{}_{yx} + \Gamma^x{}_{x\lambda}\Gamma^\lambda{}_{yy} - \Gamma^x{}_{y\lambda}\Gamma^\lambda{}_{yx}.$$

Here $\Gamma^x{}_{yy} = 0$ and $\Gamma^\lambda{}_{yy} = 0$ for every $\lambda$, so the first and third terms are absent.
The second is $-\partial_y\left(\tfrac{1}{2}\beta\,\partial_y\beta\right) = -\tfrac{1}{2}\left(\partial_y\beta\right)^2 - \tfrac{1}{2}\beta\,\partial_y^2\beta$, and the fourth is

$$-\left[\left(-\tfrac{1}{2}\left(1+\beta^2\right)\partial_y\beta\right)\left(\tfrac{1}{2}\partial_y\beta\right) + \left(\tfrac{1}{2}\beta\,\partial_y\beta\right)^2\right] = \tfrac{1}{4}\left(\partial_y\beta\right)^2,$$

so

$$R^x{}_{yxy} = -\tfrac{1}{4}\left(\partial_y\beta\right)^2 - \tfrac{1}{2}\beta\,\partial_y^2\beta.$$

Lowering the first index has to go through both $g_{xt}$ and $g_{xx}$, and the second derivative cancels against $R^t{}_{yxy} = -\tfrac{1}{2}\partial_y^2\beta$:

$$R_{xyxy} = g_{xt}R^t{}_{yxy} + g_{xx}R^x{}_{yxy} = \tfrac{1}{2}\beta\,\partial_y^2\beta - \tfrac{1}{4}\left(\partial_y\beta\right)^2 - \tfrac{1}{2}\beta\,\partial_y^2\beta = -\tfrac{1}{4}\left(\partial_y\beta\right)^2.$$

A perfect square with a minus sign in front of it, left standing after the second derivatives have gone.
That is the first appearance of the object the whole entry turns on.

The published tables carry one hundred and twenty nonzero components with an upper index and one hundred with all four down, each of them this computation with different indices.

---

## Step 7. The Ricci tensor

Contracting on the first lower index as Step 1 says, the transverse diagonal component is the shortest:

$$R_{yy} = R^t{}_{yty} + R^x{}_{yxy} + R^y{}_{yyy} + R^z{}_{yzy} = \left(\tfrac{1}{2}\beta\,\partial_y^2\beta + \tfrac{3}{4}\left(\partial_y\beta\right)^2\right) + \left(-\tfrac{1}{4}\left(\partial_y\beta\right)^2 - \tfrac{1}{2}\beta\,\partial_y^2\beta\right) + 0 + 0,$$

so

$$R_{yy} = \tfrac{1}{2}\left(\partial_y\beta\right)^2,$$

and likewise $R_{zz} = \tfrac{1}{2}\left(\partial_z\beta\right)^2$ and $R_{yz} = \tfrac{1}{2}\partial_y\beta\,\partial_z\beta$.
The $\beta\,\partial_y^2\beta$ pieces cancel between the two surviving terms and leave a square.
The rest of the tensor is longer but the same work:

$$R_{xx} = \partial_t\partial_x\beta + \left(\partial_x\beta\right)^2 + \beta\,\partial_x^2\beta - \tfrac{1}{2}\left(\partial_y\beta\right)^2 - \tfrac{1}{2}\left(\partial_z\beta\right)^2,$$

$$R_{xy} = \tfrac{1}{2}\partial_t\partial_y\beta + \tfrac{1}{2}\beta\,\partial_x\partial_y\beta + \partial_x\beta\,\partial_y\beta,$$

$$R_{tx} = \tfrac{1}{2}\partial_y^2\beta + \tfrac{1}{2}\partial_z^2\beta - \beta\,\partial_t\partial_x\beta - \beta\left(\partial_x\beta\right)^2 + \tfrac{1}{2}\beta\left(\partial_y\beta\right)^2 + \tfrac{1}{2}\beta\left(\partial_z\beta\right)^2 - \beta^2\,\partial_x^2\beta,$$

with $R_{tt}$, $R_{ty}$ and $R_{tz}$ longer still and printed in the entry.
Contracting once more with $g^{\mu\nu}$ of Step 4,

$$R = 2\,\partial_t\partial_x\beta + 2\left(\partial_x\beta\right)^2 + 2\beta\,\partial_x^2\beta + \tfrac{1}{2}\left(\partial_y\beta\right)^2 + \tfrac{1}{2}\left(\partial_z\beta\right)^2,$$

which under $\beta = v_s f$ is the published

$$R = 2v_s\,\partial_t\partial_x f + 2v_s'\,\partial_x f + 2v_s^2\left(\partial_x f\right)^2 + \tfrac{1}{2}v_s^2\left(\partial_y f\right)^2 + \tfrac{1}{2}v_s^2\left(\partial_z f\right)^2 + 2v_s^2f\,\partial_x^2 f.$$

The Ricci scalar is not sign definite and it is not the thing to look at.
The longitudinal terms in it can be made either sign by choosing how the bubble accelerates, which is exactly why the objection to the warp drive has to be stated as a contraction of the Einstein tensor with an observer rather than as a statement about $R$.

---

## Step 8. The Einstein tensor

Subtracting the trace,

$$G_{\mu\nu} = R_{\mu\nu} - \tfrac{1}{2}Rg_{\mu\nu},$$

gives, in $\beta$,

$$G_{tt} = -\beta\left(\partial_y^2\beta + \partial_z^2\beta\right) - \tfrac{1}{4}\left[\left(\partial_y\beta\right)^2 + \left(\partial_z\beta\right)^2\right] - \tfrac{3}{4}\beta^2\left[\left(\partial_y\beta\right)^2 + \left(\partial_z\beta\right)^2\right],$$

$$G_{tx} = \tfrac{1}{2}\left(\partial_y^2\beta + \partial_z^2\beta\right) + \tfrac{3}{4}\beta\left[\left(\partial_y\beta\right)^2 + \left(\partial_z\beta\right)^2\right],$$

$$G_{xx} = -\tfrac{3}{4}\left[\left(\partial_y\beta\right)^2 + \left(\partial_z\beta\right)^2\right],$$

together with $G_{ty}$, $G_{tz}$, $G_{xy}$, $G_{xz}$, $G_{yy}$, $G_{yz}$ and $G_{zz}$, all printed in the entry in $v_s$ and $f$.
Three of these are enough for what comes next, and it is worth noticing already that $G_{xx}$ is minus a sum of squares on its own.

---

## Step 9. The Eulerian observer

The observer this entry is about is the one at rest in the slicing: the one whose world line is orthogonal to every slice of constant $t$, the Eulerian observer of the $3+1$ split of Step 3.
Its four velocity is the unit normal to the slices,

$$n_\mu = \left(-1, 0, 0, 0\right),\qquad n^\mu = g^{\mu\nu}n_\nu = \left(1, \beta, 0, 0\right),$$

using $g^{tt} = -1$ and $g^{tx} = -\beta$ from Step 4, and it is a unit timelike vector,

$$n^\mu n_\mu = -1.$$

This observer is not at rest in the coordinates: it drifts in $x$ at $dx/d(ct) = \beta$, carried along by the shift.
Far from the bubble $\beta = 0$ and it is an ordinary inertial observer of the flat space the bubble travels through.
Inside the bubble $\beta = v_s$ and it moves with the ship.
Its proper time is $t$ itself, because the lapse is one.

Two facts about this congruence come out of the connection of Step 5 with no further work.
Its expansion is

$$\theta = \nabla_\mu n^\mu = \frac{1}{\sqrt{-g}}\partial_\mu\left(\sqrt{-g}\,n^\mu\right) = \partial_x\beta = v_s\,\partial_x f,$$

using $\det g = -1$ from Step 4.
For a bubble profile that falls off with distance from the centre, $\partial_x f$ is negative ahead of the ship and positive behind it, so space contracts in front and expands behind: the picture the history describes is this one line of algebra.
The congruence is also geodesic, since $n^\mu$ is the unit normal of a unit lapse foliation, which is why the ship rides for free.

---

## Step 10. The energy density, and the sign it cannot escape

The energy density this observer measures is

$$\rho = T_{\mu\nu}n^\mu n^\nu = \frac{c^4}{8\pi G}G_{\mu\nu}n^\mu n^\nu = \frac{c^4}{8\pi G}G^{\mu\nu}n_\mu n_\nu = \frac{c^4}{8\pi G}G^{tt},$$

the last step because $n_\mu$ has only a time component and that component is $-1$.
So the published $G^{tt}$ is the energy density, up to the constant $c^4/8\pi G$, and nothing else in the entry needs to be read to answer the question the warp drive raises.

Raising both indices with Step 4,

$$G^{tt} = \left(g^{tt}\right)^2G_{tt} + 2g^{tt}g^{tx}G_{tx} + \left(g^{tx}\right)^2G_{xx} = G_{tt} + 2\beta\,G_{tx} + \beta^2G_{xx},$$

and putting in the three components of Step 8, the terms in $\beta\left(\partial_y^2\beta + \partial_z^2\beta\right)$ cancel between the first two, and the coefficients of $\beta^2\left[\left(\partial_y\beta\right)^2 + \left(\partial_z\beta\right)^2\right]$ come to $-\tfrac{3}{4} + \tfrac{3}{2} - \tfrac{3}{4} = 0$.
What is left is

$$G^{tt} = -\tfrac{1}{4}\left[\left(\partial_y\beta\right)^2 + \left(\partial_z\beta\right)^2\right] = -\tfrac{1}{4}v_s^2\left[\left(\partial_y f\right)^2 + \left(\partial_z f\right)^2\right],$$

which is what the entry publishes, and therefore

$$\rho = -\frac{c^4}{32\pi G}v_s^2\left[\left(\partial_y f\right)^2 + \left(\partial_z f\right)^2\right].$$

Read it slowly, because every clause of the objection to the warp drive is in it.

It is minus a sum of squares, so it is never positive, for any shape function and any speed.
It is proportional to $v_s^2$, so it does not care which way the bubble goes, it grows as the square of the speed, and it shrinks without ever vanishing as the bubble is slowed: an arbitrarily slow warp bubble still needs exotic matter, which is the Lobo and Visser result the history cites.
It involves only the transverse gradient of $f$, the variation across the direction of travel, so it vanishes only for a bubble with no transverse structure at all, which is to say for no bubble.
It says nothing about $\partial_t f$, nor about $v_s'$, so it is untouched by how the bubble accelerates: a warp drive that starts gently is exotic in exactly the way one that does not is.
And the weak energy condition, which asks that $T_{\mu\nu}u^\mu u^\nu \geq 0$ for every timelike $u$, is violated at every point of the wall by the observer who is most at home there.

For Alcubierre's radial profile the transverse gradient is

$$\left(\partial_y f\right)^2 + \left(\partial_z f\right)^2 = \frac{y^2 + z^2}{r_s^2}\left(\frac{df}{dr_s}\right)^2,$$

and $\rho$ becomes his expression of 1994: the negative energy sits in a belt around the bubble, thickest on the equator where $y^2 + z^2$ is largest and vanishing on the axis of travel, where the transverse gradient is zero.

---

## Step 11. The same density from the Hamiltonian constraint

The result of Step 10 is worth deriving a second way, because the second way explains the sign rather than merely producing it.

For any spacelike slicing, the Gauss equation contracted twice gives the Hamiltonian constraint,

$$R^{(3)} + K^2 - K_{ij}K^{ij} = 2G_{\mu\nu}n^\mu n^\nu = \frac{16\pi G}{c^4}\rho,$$

with $R^{(3)}$ the Ricci scalar of the slice, $K_{ij}$ its extrinsic curvature and $K$ its trace.

Both pieces on the left are known here.
The slices are flat, $\gamma_{ij} = \delta_{ij}$, so

$$R^{(3)} = 0,$$

and the whole of the curvature is in the way the slices are embedded.
The extrinsic curvature needs no new machinery either, since with $n_j = 0$ for spatial $j$ and $n_t = -1$,

$$K_{ij} = -\nabla_i n_j = -\left(\partial_i n_j - \Gamma^\lambda{}_{ij}n_\lambda\right) = -\Gamma^t{}_{ij},$$

so the table of Step 5 has it already.
Reading the three symbols $\Gamma^t{}_{xx} = \partial_x\beta$, $\Gamma^t{}_{xy} = \tfrac{1}{2}\partial_y\beta$ and $\Gamma^t{}_{xz} = \tfrac{1}{2}\partial_z\beta$ off it, and noting that $\Gamma^t{}_{yy}$, $\Gamma^t{}_{yz}$ and $\Gamma^t{}_{zz}$ are all absent from it,

$$K_{xx} = -\partial_x\beta,\qquad K_{xy} = -\tfrac{1}{2}\partial_y\beta,\qquad K_{xz} = -\tfrac{1}{2}\partial_z\beta,$$

with every purely transverse component zero.
This is the symmetrised gradient of the shift, as a unit lapse and a time independent spatial metric require.
Only one of these components is diagonal, and that is the whole story.
The trace is $K = -\partial_x\beta$, so

$$K^2 = \left(\partial_x\beta\right)^2,\qquad K_{ij}K^{ij} = \left(\partial_x\beta\right)^2 + \tfrac{1}{2}\left(\partial_y\beta\right)^2 + \tfrac{1}{2}\left(\partial_z\beta\right)^2,$$

and the longitudinal parts cancel exactly:

$$K^2 - K_{ij}K^{ij} = -2\left(K_{xy}^2 + K_{xz}^2\right) = -\tfrac{1}{2}\left[\left(\partial_y\beta\right)^2 + \left(\partial_z\beta\right)^2\right].$$

Hence

$$\rho = \frac{c^4}{16\pi G}\left(K^2 - K_{ij}K^{ij}\right) = -\frac{c^4}{32\pi G}v_s^2\left[\left(\partial_y f\right)^2 + \left(\partial_z f\right)^2\right],$$

the same expression, reached without computing a single Riemann component.
The sign convention for $K_{ij}$ does not enter, since both terms are quadratic in it.

Now the structure is visible.
On a flat slice the energy density is $K^2 - K_{ij}K^{ij}$, and that combination is the difference between the square of a trace and the sum of the squares of all the entries.
A shift with one component and a gradient across it puts everything it has into off diagonal entries $K_{xy}$ and $K_{xz}$, which contribute to $K_{ij}K^{ij}$ and not to $K$.
The only part that could have given a positive contribution, $K_{xx}^2$, appears in both terms and cancels.
The negative energy of the warp drive is not a feature of Alcubierre's profile or of any other: it is what flat slices and a one component shift make inevitable, and the only way out of it is to give up one of those two, which is what every later variant of the warp drive does.

---

## Step 12. How much exotic matter, and why that is the objection

Integrating the density of Step 10 over a slice, for a radial profile $f(r_s)$, and averaging $\left(y^2 + z^2\right)/r_s^2$ over the sphere at fixed $r_s$, which gives $2/3$,

$$E = \int\rho\,d^3x = -\frac{c^4}{32\pi G}v_s^2\cdot\frac{8\pi}{3}\int_0^\infty r_s^2\left(\frac{df}{dr_s}\right)^2dr_s = -\frac{c^4}{12G}v_s^2\int_0^\infty r_s^2\left(\frac{df}{dr_s}\right)^2dr_s.$$

For a bubble of radius $R$ whose wall has thickness $\Delta$, the derivative is of order $1/\Delta$ across a range of order $\Delta$, so the integral is of order $R^2/\Delta$ and

$$E \sim -\frac{c^4}{12G}\frac{v_s^2R^2}{\Delta}.$$

The coefficient is $c^4/12G \approx 1.0\times10^{43}$ joules per metre, so a bubble of radius a hundred metres with a wall a metre thick, moving at the speed of light, needs about $10^{47}$ joules of negative energy, which is around half the rest energy of the Sun.
That is already absurd, and it is the optimistic case.
The quantum inequalities Pfenning and Ford applied to this geometry force the wall to be thinner than about a hundred Planck lengths, and at $\Delta \sim 10^{-33}$ metres the same formula gives something of order $10^{80}$ joules, ten orders of magnitude beyond the mass energy of the visible universe.
The $1/\Delta$ is what does the damage: the requirement is not that some exotic matter exist but that an amount of it scaling inversely with the thinness of a wall that quantum field theory insists must be thin.

This is the arithmetic behind the sentence in the entry's history, and it is worth seeing that it follows from one published component and a single integral.

---

## Step 13. The light cones, and the horizon at $f = 1 - 1/v_s$

Everything so far has been about the source.
The kinematics is worth one step, because it is where the claim that nothing moves faster than light locally is settled.

For a photon travelling along the axis, with $dy = dz = 0$, the line element of Step 2 gives $\left(dx - \beta\,dx^0\right)^2 = \left(dx^0\right)^2$, so

$$\frac{dx}{d(ct)} = \beta \pm 1.$$

The local light cone is the ordinary one, tipped over by the shift.
The ship sits at the centre of the bubble, where $f = 1$ and so $\beta = v_s$, and its own coordinate velocity is $dx/d(ct) = v_s$, which lies strictly between $v_s - 1$ and $v_s + 1$.
It is inside its own light cone at every point of its world line, whatever $v_s$ is.
Far outside the bubble $\beta = 0$ and the cones are upright.
Nothing anywhere moves faster than light with respect to the local geometry; what moves faster than light is the bubble, with respect to the flat space it is passing through, and the coordinate speed of a distant observer is not a velocity that special relativity ever bounded.

Now measure the same rays from the bubble, by using the displacement $x - x_s(t)$ instead of $x$.
Subtracting the bubble's own $v_s$,

$$\frac{d\left(x - x_s\right)}{d(ct)} = v_s\left(f - 1\right) \pm 1.$$

The forward ray, the one with $+1$, moves forward relative to the bubble as long as $v_s(f-1) + 1 > 0$, which is to say as long as

$$f > 1 - \frac{1}{v_s}.$$

For $v_s \leq 1$ the right hand side is not positive and every ray gets through: a subluminal bubble is transparent to its own occupants.
For $v_s > 1$ the surface $f = 1 - 1/v_s$ lies inside the front wall, and a light ray sent forward from the ship cannot cross it.
For a bubble that is unchanging in its own frame this surface is a horizon, and the ship, sitting at the centre in flat space, cannot signal to the front of its own bubble, let alone steer it.
Whatever creates and sustains the wall has to have been put in place beforehand, along the whole path, which is a second objection to the warp drive and one that costs nothing to see: it is two lines of the line element.

---

## Step 14. The Weyl tensor, computed rather than copied

A Weyl tensor equals a Riemann tensor only where the Ricci tensor vanishes, and here it vanishes nowhere on the wall, so the trace free part has to be built:

$$C_{\mu\nu\rho\sigma} = R_{\mu\nu\rho\sigma} - \left(g_{\mu[\rho}R_{\sigma]\nu} - g_{\nu[\rho}R_{\sigma]\mu}\right) + \tfrac{1}{3}R\,g_{\mu[\rho}g_{\sigma]\nu},$$

in the four dimensional form the checker uses.

The difference is not cosmetic.
Take the purely transverse slot, the one with no $x$ and no $t$ in it.
Step 6 found $R_{xyxy} = -\tfrac{1}{4}\left(\partial_y\beta\right)^2$, but drop the $x$ as well and the Riemann component is empty,

$$R_{yzyz} = 0,$$

while the same slot of the Weyl tensor is

$$C_{yzyz} = \tfrac{1}{3}\left[\partial_t\partial_x\beta + \left(\partial_x\beta\right)^2 + \beta\,\partial_x^2\beta\right] - \tfrac{1}{6}\left[\left(\partial_y\beta\right)^2 + \left(\partial_z\beta\right)^2\right],$$

which is not zero, and which depends on the longitudinal structure of the bubble that $R_{yzyz}$ knows nothing about.
Thirty six slots of each published Weyl variant are nonzero where the corresponding Riemann slot vanishes, and none of them would exist in a file that had copied one block into the other.
There is no slot the other way round: every nonzero Riemann component of this spacetime leaves a nonzero Weyl component behind.

The entry publishes one hundred and fifty six nonzero components of $C^\mu{}_{\nu\rho\sigma}$ and one hundred and thirty six of $C_{\mu\nu\rho\sigma}$, all of them computed from the formula above.

---

## Step 15. The Kretschmann scalar

The invariant

$$K = R_{\mu\nu\rho\sigma}R^{\mu\nu\rho\sigma}$$

is a polynomial of forty two terms in $\beta$ and its first and second derivatives when it is expanded, and it collapses to six when the right groupings are taken.
Write

$$A = \partial_t\partial_x\beta + \left(\partial_x\beta\right)^2 + \beta\,\partial_x^2\beta = \partial_x\left(\partial_t\beta + \beta\,\partial_x\beta\right),$$

$$B_i = \partial_t\partial_i\beta + \beta\,\partial_x\partial_i\beta + 2\,\partial_x\beta\,\partial_i\beta \quad (i = y, z),$$

$$Q = \left(\partial_y\beta\right)^2 + \left(\partial_z\beta\right)^2,$$

$$H = \left(\partial_x\partial_y\beta\right)^2 + \left(\partial_x\partial_z\beta\right)^2 + \left(\partial_y^2\beta\right)^2 + \left(\partial_z^2\beta\right)^2 + 2\left(\partial_y\partial_z\beta\right)^2.$$

The quantity $\partial_t\beta + \beta\,\partial_x\beta$ inside $A$ is $n^\mu\partial_\mu\beta$, the rate of change of the shift along the Eulerian congruence of Step 9.
In these,

$$K = 4A^2 + 2B_y^2 + 2B_z^2 - 2AQ + \tfrac{11}{4}Q^2 - 2H,$$

which is what the entry publishes, written out in $v_s$ and $f$, where $Q = v_s^2\left[\left(\partial_y f\right)^2 + \left(\partial_z f\right)^2\right]$ is the same transverse gradient the energy density of Step 10 is built from.
The grouped form and the expanded one are the same expression, which is what the checker tests; the grouped one is published because it renders in under half the width and because it shows where the pieces come from.

Two things are worth saying about it.

It is a polynomial, with no denominators anywhere, so wherever the shape function is smooth the curvature invariant is finite.
The warp drive has no curvature singularity at all: unlike every black hole entry in this collection, there is nowhere the geometry breaks down.
What makes it impossible is entirely the source it demands, never the geometry it produces, and this scalar is where that is visible.

It also vanishes identically wherever $f$ is constant, which is to say everywhere except the wall, confirming from the invariant side what Step 5 said from the connection: the inside of the bubble and the outside are flat pieces of Minkowski space, and all of the curvature lives in the shell between them.

---

## Step 16. The geodesic equations

The equations are

$$\ddot{x}^\mu + \Gamma^\mu{}_{\nu\rho}\dot{x}^\nu\dot{x}^\rho = 0$$

with the symbols of Step 5 and the dots of Step 1, that is, derivatives of the chart coordinates $\left(ct, x, y, z\right)$ with respect to an affine parameter.
Writing them out with $\beta = v_s f$, the two transverse ones are short:

$$\ddot{y} + v_s\,\partial_y f\,\dot{t}\dot{x} - v_s^2f\,\partial_y f\,\dot{t}^2 = 0,\qquad \ddot{z} + v_s\,\partial_z f\,\dot{t}\dot{x} - v_s^2f\,\partial_z f\,\dot{t}^2 = 0,$$

and the $t$ and $x$ equations are the longer ones the entry prints.
Each cross term appears once for each ordering of its two indices, which is where the factor of two in $-2v_s^2f\,\partial_xf\,\dot{t}\dot{x}$ comes from.

The equations say the thing Step 9 said: where $f$ is constant every $\Gamma$ vanishes and the equations become $\ddot{x}^\mu = 0$, so an observer inside the bubble, or far outside it, moves in a straight line at constant speed through a flat patch.
Only a world line that crosses the wall feels anything at all.

---

## Step 17. Every published expression is dimensionally consistent

The chart coordinates are $\left(ct, x, y, z\right)$, all of dimension $L$, and both parameters are dimensionless.
A metric component then carries

$$[g_{\mu\nu}] = \frac{L^2}{[x^\mu][x^\nu]} = 1,$$

a Christoffel symbol carries $L^{-1}$, a Riemann, Ricci, Einstein or Weyl component carries $L^{-2}$, and $K$ carries $L^{-4}$, each corrected by $[x^\mu]/L = 1$ per index, which is to say not corrected at all: in this chart every index is a length and the bookkeeping is as simple as it gets.

What does the work here is the derivative convention.
A printed $\partial_i f$ is $\partial f/\partial x^i$ with $f$ dimensionless, so it carries $L^{-1}$; a printed $\partial_t f$ is $c^{-1}\partial f/\partial t$, which carries $T^{-1}\cdot T/L = L^{-1}$ as well; and $v_s'$ is $dv_s/d(ct)$, again $L^{-1}$.
Every derivative of every order therefore carries one inverse length per order, and a published component is balanced exactly when each of its terms carries as many derivatives as its rank demands.

Four samples.

$g_{tt} = -1 + v_s^2f^2$ is dimensionless in both terms, since both parameters are.

$\Gamma^t{}_{xy} = \tfrac{1}{2}v_s\,\partial_y f$ carries one derivative and so $L^{-1}$, as a Christoffel symbol with any three indices of this chart must.

$G^{tt} = -\tfrac{1}{4}v_s^2\left[\left(\partial_y f\right)^2 + \left(\partial_z f\right)^2\right]$ carries two derivatives and so $L^{-2}$, and the energy density built from it carries $\left(c^4/G\right)L^{-2} = MLT^{-2}\cdot L^{-2} = ML^{-1}T^{-2}$, which is an energy per unit volume.
That is the check that the factor of $c^4/8\pi G$ in Step 10 is the right one.

The Kretschmann scalar carries four derivatives in every term, grouped as Step 15 publishes it or expanded into its forty two, which is $L^{-4}$.

The geodesic equations are measured against their own second derivatives.
In the $y$ equation of Step 16, $\ddot{y}$ carries $L/\lambda^2$, and the term $v_s\,\partial_yf\,\dot{t}\dot{x}$ carries

$$\frac{1}{L}\cdot\frac{L}{\lambda}\cdot\frac{L}{\lambda} = \frac{L}{\lambda^2},$$

which matches, and matches only because $\dot{t}$ is $d(ct)/d\lambda$ rather than $dt/d\lambda$.
On the other reading that term would be one factor of $c$ away from its left hand side, and so would every one of the twelve terms of the $x$ equation.

---

## Step 18. The bubble frame is this chart with $f$ replaced by $f - 1$

The natural second chart for a warp drive is the one in which the ship is at rest, reached by measuring $x$ from the centre of the bubble,

$$x' = x - x_s(t),\qquad dx = dx' + c\,v_s\,dt.$$

Substituting into the line element of Step 2,

$$dx - v_s f\,c\,dt = dx' + c\,v_s\left(1 - f\right)dt = dx' - v_s\left(f - 1\right)c\,dt,$$

so

$$ds^2 = -c^2dt^2 + \left(dx' - v_s\left(f-1\right)c\,dt\right)^2 + dy^2 + dz^2.$$

That is the same line element with the shape function shifted by a constant.
Since the entry's $f$ is arbitrary, and since $f$ and $f - 1$ have the same derivatives of every order, every component published in this entry is already a component of the bubble frame under $f \to f - 1$.
Nothing new would be learned by printing the second chart, and the entry does not: it would be six hundred and eighty expressions restating the six hundred and eighty already there.

The substitution is a good check of the physics, too.
The energy density of Step 10 depends only on the transverse derivatives of $f$, which are unchanged by it, so the Eulerian observers of the two charts, who are the same observers, measure the same negative density.
The exotic matter is not an artefact of describing the bubble from the frame it is flying through.

---

## Step 19. What the entry publishes, and the declarations the checker needs

The one system publishes six hundred and eighty expressions: the line element, six metric components and six inverse ones, thirty Christoffel symbols with an upper index and twenty three with all three down, one hundred and twenty Riemann components with an upper index and one hundred with all four down, sixteen Ricci components in each of three variants, the Ricci scalar, the Kretschmann scalar, sixteen Einstein components in each of three variants, one hundred and fifty six Weyl components with an upper index and one hundred and thirty six with all four down, and four geodesic equations.

For the checker to read it, `DIMENSIONS` in `verify_metrics.py` needs one line:

    ("alcubierre", "cartesian"): {
        "t": "T", "x": "L", "y": "L", "z": "L", "v_s": "1", "f": "1",
    },

The declaration `"t": "T"` is what tells the checker that $t$ is a time and so that the chart coordinate is $ct$; every factor of $c$ in the comparison follows from it, as does the reading of $\partial_t$ and of the prime on $v_s$ as chart derivatives.
The prime rather than a dot is deliberate, and the next entry that names a function with a subscript will want to know why: the reader in `verify_metrics.py` binds `\dot{...}` and `\partial_x` to a name of letters only, so `\dot{v_s}` is not a rate it can read, while `v_s'` is, and a parameter whose name carries an underscore therefore takes its time derivative as a prime.
Declaring both parameters dimensionless is what puts the bubble velocity in units of $c$: the shift is $c\,v_s f$, a velocity, and the $c$ is written out in the line element rather than hidden in $v_s$.

The entry needs no line in `PARAMETER_RELATIONS`, and it is worth being clear about why, since it is the one place this spacetime could have been published with an identity that is not one.
Alcubierre's own shape function satisfies $\partial_t f = -v_s\,\partial_x f$, and using it would shorten a good many of the published components.
It would also be false for a bubble whose profile changes shape as it travels, which the entry allows, so no published value uses it.
The checker compares against a free $f$ of four coordinates and a free $v_s$ of one, which is the claim the entry actually makes.

Running

    /tmp/mfs-venv/bin/python _tools/derivations/verify_metrics.py --system alcubierre/cartesian

reports no disagreements and no dimensional failures.
