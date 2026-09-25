# The Bianchi type I anisotropic cosmology

The Bianchi type I anisotropic cosmology is worked here in one coordinate system, the type I Cartesian chart.
Every component of the metric and of its curvature is derived in order, from the line element down to the geodesic equations.
We leave nothing as an exercise and assert nothing we do not compute.

The companion script `_tools/derivations/verify_metrics.py` does the same work in sympy and checks every component against its own result, so the algebra is checkable by hand and by machine independently.

Like FRW, this metric is not a formula: its scale factors are undetermined functions, and no field equation is imposed on them in computing any tensor.
Unlike FRW there are three of them, and their differences are the whole of the anisotropy.
That is deliberate: type I is at its most useful with matter in it, so the Ricci and Einstein tensors are computed for arbitrary $a_i(t)$ and the matter is read off them afterwards.
The density and the three principal pressures follow in Step 10, and emptying the spacetime out in Step 12 recovers the Kasner solution, which is a spacetime in its own right.

---

## Step 1. Conventions

The signature is $(-,+,+,+)$ and the Christoffel symbols are those of the Levi-Civita connection,

$$\Gamma^\mu_{\nu\rho} = \tfrac{1}{2}g^{\mu\alpha}\left(\partial_\nu g_{\alpha\rho} + \partial_\rho g_{\alpha\nu} - \partial_\alpha g_{\nu\rho}\right),$$

with the lowered symbol $\Gamma_{\mu\nu\rho} = g_{\mu\alpha}\Gamma^\alpha_{\nu\rho}$.
The Riemann tensor is

$$R^\mu{}_{\nu\rho\sigma} = \partial_\rho \Gamma^\mu_{\nu\sigma} - \partial_\sigma \Gamma^\mu_{\nu\rho} + \Gamma^\mu_{\rho\lambda}\Gamma^\lambda_{\nu\sigma} - \Gamma^\mu_{\sigma\lambda}\Gamma^\lambda_{\nu\rho},$$

which fixes the sign of every Riemann component.

The Ricci tensor is the standard contraction,

$$R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu},$$

settled on 18 September 2026 for every spacetime alike.
Before that the Ricci tensor was contracted on the last lower index, which is the same tensor with the opposite sign, and every Ricci and Einstein tensor has since been brought over, in `verify_metrics.py` as well; for a spacetime with matter in it the choice is not a matter of taste, because only the standard contraction gives an expanding universe a positive energy density, which is Step 11.

Factors of $c$ and $G$ are kept explicit.
The chart, here as for every other spacetime, is the one whose time coordinate is $x^0 = ct$.
Schwarzschild shows the convention plainly, quoting $g_{tt} = -(1-r_s/r)$ against a line element whose time term is $-(1-r_s/r)c^2dt^2$.
Since $t$ here is a time, the chart coordinate is $ct$, and every component written against an index $t$ is a component in that chart even though the index is written with the bare letter.

Because the rescaling $t \to ct$ is linear with constant coefficients, the rule is arithmetic: a component in the chart is the component taken with the bare coordinate, multiplied by $c$ once for every upper $t$ index and divided by $c$ once for every lower one.
The Christoffel symbols obey the same rule as the tensors, since the inhomogeneous term in their transformation law carries a second derivative of the coordinate change, which vanishes for a linear one.
A scalar has no index at all and is therefore the same number in both charts.

Unlike Kasner and Vaidya, type I is computed in the chart from the start: their metric functions are explicit in $t$ and in $u$, and rewriting them in the chart coordinate would bury a power of $c$ inside each one.
Here the metric functions are undetermined, so nothing is lost by differentiating with respect to $x^0$ directly, and a prime is defined to mean exactly that:

$$a_i' \equiv \frac{da_i}{d(ct)} = \frac{1}{c}\frac{da_i}{dt}, \qquad a_i'' \equiv \frac{d^2a_i}{d(ct)^2} = \frac{1}{c^2}\frac{d^2a_i}{dt^2}.$$

Every factor of $c$ the chart demands is in those two definitions, and no other factor of $c$ appears in any component.
That is not a way of hiding them: the scale factors are dimensionless, so $a_i'$ carries an inverse length and $a_i''$ an inverse length squared, which is what a chart component with those index patterns has to carry, and every term is held to that requirement in Step 16.
The relation to the notation a cosmologist uses is $\dot{a}_i/a_i = c\,a_i'/a_i$, so the Hubble rate along the $x$ axis is $c$ times the rate per unit chart length.

Two abbreviations shorten the working:

$$H_i \equiv \frac{a_i'}{a_i}, \qquad u_i \equiv \frac{a_i''}{a_i}.$$

$H_i$ is the expansion rate along the $i$ axis measured per unit chart length, and $u_i$ its acceleration.
The checker has no way to declare them, so every component it checks carries $H_i$ and $u_i$ written out.

---

## Step 2. The line element, and what type I means

$$ds^2 = -c^2dt^2 + a_1(t)^2dx^2 + a_2(t)^2dy^2 + a_3(t)^2dz^2,$$

with $a_1, a_2, a_3$ positive and twice differentiable, and $x$, $y$, $z$ each running over the whole real line.

The Bianchi classification is a statement about the symmetry group of the spatial slices, and type I is the case where that group is the abelian $\mathbb{R}^3$ of translations.
Concretely: the slices $t = \text{const}$ are ordinary flat three dimensional space, every point of a slice is equivalent to every other, and the only thing that distinguishes one direction from another is how fast the slice is being stretched along it.
Three functions of time carry all of that, and there is nothing else in the geometry to carry anything.
The three invariants $(n_1, n_2, n_3)$ of the Bianchi classification are all zero here, which is the algebraic way of saying the same thing.

A general type I metric carries a full symmetric spatial matrix $g_{ij}(t)$, and the diagonal form used here is the standard restriction of it: it is the case where the principal axes of the expansion do not rotate, so one constant linear change of the comoving coordinates lines them up with $x$, $y$ and $z$ once and for all.
Only that case is treated, and it is the case the literature means when it says Bianchi I.

Writing $a_1 = a_2 = a_3$ gives the spatially flat FLRW model, so this line element contains the isotropic universe as the one point where the three functions coincide.
Everything anisotropic is carried by their differences.

The chart is written $t \in (0,\infty)$ because the models of interest reach a singularity at a finite time, which is placed at $t = 0$.
No result depends on that choice; the results hold on any interval where the three scale factors are smooth and positive.

In the chart $x^0 = ct$ the line element reads

$$ds^2 = -(dx^0)^2 + a_1^2dx^2 + a_2^2dy^2 + a_3^2dz^2,$$

and every component from Step 3 to Step 15 is computed there.

---

## Step 3. The metric matrix, its determinant and its inverse

In the order $(t,x,y,z)$, where the time slot is the chart's $x^0 = ct$,

$$g_{\mu\nu} = \begin{pmatrix} -1 & 0 & 0 & 0 \\ 0 & a_1^2 & 0 & 0 \\ 0 & 0 & a_2^2 & 0 \\ 0 & 0 & 0 & a_3^2 \end{pmatrix}, \qquad g^{\mu\nu} = \begin{pmatrix} -1 & 0 & 0 & 0 \\ 0 & a_1^{-2} & 0 & 0 \\ 0 & 0 & a_2^{-2} & 0 \\ 0 & 0 & 0 & a_3^{-2} \end{pmatrix}.$$

The matrix is diagonal, so the inverse is the reciprocal entry by entry, and $g_{\mu\alpha}g^{\alpha\nu} = \delta_\mu^\nu$ by inspection.

The determinant is

$$\det g = -a_1^2a_2^2a_3^2, \qquad \sqrt{-g} = a_1a_2a_3 \equiv V.$$

$V$ is the comoving volume element, the volume of a coordinate box of unit side, and it is the one combination of the three scale factors that behaves like a single isotropic scale factor.
It runs through the whole of what follows, and a singularity of this spacetime is a time where $V \to 0$.

---

## Step 4. Every Christoffel symbol

Only the three spatial diagonal components of $g$ depend on anything, and each depends only on $x^0$.
Throughout, $\partial_0 = \partial/\partial x^0$ is the derivative the prime denotes, while the index on a component is written $t$, since $t$ is the letter that names the chart slot $x^0$:

$$\partial_0 g_{ii} = 2a_ia_i', \qquad \text{no sum on } i,$$

and every other derivative of the metric vanishes.
Feed that into the connection formula.

A symbol with three spatial indices needs a spatial derivative of the metric and gets none, so it vanishes.
A symbol with three time indices needs $\partial_0 g_{00}$, which is zero because $g_{00} = -1$ is constant, so $\Gamma^t{}_{tt} = 0$.
Two families survive.

For $\Gamma^t{}_{ii}$ the only surviving term is the one with the minus sign, since $g_{ti} = 0$ kills the first two:

$$\Gamma^t{}_{ii} = \tfrac{1}{2}g^{tt}\left(-\partial_0 g_{ii}\right) = \tfrac{1}{2}(-1)(-2a_ia_i') = a_ia_i'.$$

For $\Gamma^i{}_{ti}$ the surviving term is the first, and $g^{ii} = a_i^{-2}$:

$$\Gamma^i{}_{ti} = \Gamma^i{}_{it} = \tfrac{1}{2}g^{ii}\,\partial_0 g_{ii} = \tfrac{1}{2}a_i^{-2}\,2a_ia_i' = \frac{a_i'}{a_i} = H_i.$$

There is no $\Gamma^i{}_{jk}$ with $j \neq k$ and no $\Gamma^t{}_{ij}$ with $i \neq j$, because the metric is diagonal and each spatial component depends on the time alone.
The complete list is nine nonzero symbols, counting the two index orders separately:

$$\Gamma^t{}_{xx} = a_1a_1', \quad \Gamma^t{}_{yy} = a_2a_2', \quad \Gamma^t{}_{zz} = a_3a_3',$$

$$\Gamma^x{}_{tx} = \Gamma^x{}_{xt} = \frac{a_1'}{a_1}, \quad \Gamma^y{}_{ty} = \Gamma^y{}_{yt} = \frac{a_2'}{a_2}, \quad \Gamma^z{}_{tz} = \Gamma^z{}_{zt} = \frac{a_3'}{a_3}.$$

That is every Christoffel symbol with its first index up.

---

## Step 5. The lowered Christoffel symbols

$\Gamma_{\mu\nu\rho} = g_{\mu\alpha}\Gamma^\alpha{}_{\nu\rho}$, and the metric is diagonal, so each symbol is multiplied by the single diagonal entry matching its first index:

$$\Gamma_{tii} = g_{tt}\Gamma^t{}_{ii} = -a_ia_i', \qquad \Gamma_{iti} = \Gamma_{iit} = g_{ii}\Gamma^i{}_{ti} = a_i^2\frac{a_i'}{a_i} = a_ia_i'.$$

Nine nonzero symbols again, differing from Step 4 only in the time row, where the sign flips.
That is every Christoffel symbol with every index lowered.

---

## Step 6. The Riemann tensor, first family: the time and space planes

Take $i$ to be one of $x,y,z$ with $a_i$ the matching scale factor, and no sum on $i$ anywhere in this step.

$$R^t{}_{iti} = \partial_0\Gamma^t{}_{ii} - \partial_i\Gamma^t{}_{it} + \Gamma^t{}_{t\lambda}\Gamma^\lambda{}_{ii} - \Gamma^t{}_{i\lambda}\Gamma^\lambda{}_{it}.$$

The second term vanishes because $\Gamma^t{}_{it} = 0$.
The third vanishes because $\Gamma^t{}_{t\lambda} = 0$ for every $\lambda$.
The fourth has one surviving value of $\lambda$, namely $\lambda = i$:

$$\Gamma^t{}_{ii}\Gamma^i{}_{it} = a_ia_i'\cdot\frac{a_i'}{a_i} = (a_i')^2.$$

The first term is

$$\partial_0\left(a_ia_i'\right) = (a_i')^2 + a_ia_i'',$$

so the two squares cancel and

$$R^t{}_{iti} = a_ia_i'', \qquad R^t{}_{iit} = -a_ia_i''.$$

The same calculation with the spatial index raised instead gives the other member of the family:

$$R^i{}_{tti} = \partial_0\Gamma^i{}_{ti} - \partial_i\Gamma^i{}_{tt} + \Gamma^i{}_{t\lambda}\Gamma^\lambda{}_{ti} - \Gamma^i{}_{i\lambda}\Gamma^\lambda{}_{tt}.$$

The second and fourth terms vanish because $\Gamma^i{}_{tt} = 0$ and $\Gamma^\lambda{}_{tt} = 0$.
The third has $\lambda = i$ only and contributes $H_i^2$.
The first is

$$\partial_0\left(\frac{a_i'}{a_i}\right) = \frac{a_i''}{a_i} - H_i^2,$$

so again the squares cancel:

$$R^i{}_{tti} = \frac{a_i''}{a_i} = u_i, \qquad R^i{}_{tit} = -u_i.$$

Twelve components in all across the three axes, which are the components of the Riemann tensor with its first index up that carry a time index.
Both members are built from $a_i''$ alone: in a type I universe the curvature of a plane containing the time direction is the acceleration of the scale factor in that direction, and an axis expanding at a constant rate contributes nothing.

---

## Step 7. The Riemann tensor, second family: the spatial planes

Now take $i \neq j$, both spatial, with no sum on either.

$$R^i{}_{jij} = \partial_i\Gamma^i{}_{jj} - \partial_j\Gamma^i{}_{ji} + \Gamma^i{}_{i\lambda}\Gamma^\lambda{}_{jj} - \Gamma^i{}_{j\lambda}\Gamma^\lambda{}_{ji}.$$

The first two terms are spatial derivatives of quantities that depend on time alone, so both vanish.
The fourth vanishes because $\Gamma^i{}_{j\lambda}$ is zero for every $\lambda$ when $i \neq j$.
The third has one surviving value, $\lambda = t$:

$$\Gamma^i{}_{it}\Gamma^t{}_{jj} = \frac{a_i'}{a_i}\,a_ja_j' = \frac{a_j\,a_i'a_j'}{a_i},$$

so

$$R^i{}_{jij} = \frac{a_j\,a_i'a_j'}{a_i}, \qquad R^i{}_{jji} = -\frac{a_j\,a_i'a_j'}{a_i}.$$

Twelve more components, four for each of the three coordinate planes, which completes the Riemann tensor with its first index up at twenty four.
This family is built from the product of two expansion rates and carries no second derivative at all: the curvature of a spatial plane is the shear between the two axes spanning it.

Lowering the first index with the diagonal metric turns the two families into

$$R_{titi} = -a_ia_i'' = -R_{tiit} = -R_{itti} = R_{itit}, \qquad R_{ijij} = a_ia_j\,a_i'a_j' = -R_{ijji},$$

and $R_{ijij} = R_{jiji}$, symmetric under exchanging the two axes as the pair symmetry of Riemann requires.
That is the Riemann tensor with every index lowered, twenty four components again.

Six independent numbers describe the whole curvature, $u_1, u_2, u_3$ and $H_1H_2, H_1H_3, H_2H_3$, which is what a diagonal metric in four dimensions with no off diagonal curvature can carry.

---

## Step 8. The Ricci tensor

Contract on the first and third indices, $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$.

For the time slot, $\alpha = t$ contributes nothing by the antisymmetry of Riemann in its last two indices, so only the three spatial values survive and Step 6 supplies them:

$$R_{tt} = \sum_i R^i{}_{tit} = -\sum_i u_i = -\left(\frac{a_1''}{a_1} + \frac{a_2''}{a_2} + \frac{a_3''}{a_3}\right).$$

For a spatial slot take $i = x$ as the representative.
The contraction runs over $\alpha = t$, which is the Step 6 family, and over the two spatial $\alpha = j \neq x$, which is the Step 7 family:

$$R_{xx} = R^t{}_{xtx} + \sum_{j\neq x} R^j{}_{xjx} = a_1a_1'' + \sum_{j \neq x}\frac{a_1\,a_j'a_1'}{a_j} = a_1a_1'' + a_1a_1'\left(\frac{a_2'}{a_2} + \frac{a_3'}{a_3}\right),$$

using $R^j{}_{xjx}$, which is Step 7 with the roles of $i$ and $j$ exchanged.
The same on the other two axes gives

$$R_{ii} = a_ia_i'' + a_ia_i'\sum_{j\neq i}\frac{a_j'}{a_j}, \qquad \text{equivalently} \qquad R^i{}_i = u_i + H_i\sum_{j\neq i}H_j.$$

Every off diagonal component vanishes, because a nonzero $R_{\mu\nu}$ with $\mu \neq \nu$ would need a Riemann component whose second and fourth indices differ, and Steps 6 and 7 leave none.

Raising indices on a diagonal metric is division:

$$R^t{}_t = -R_{tt} = \sum_i u_i, \qquad R^i{}_i = \frac{R_{ii}}{a_i^2}, \qquad R^{tt} = R_{tt}, \qquad R^{ii} = \frac{R_{ii}}{a_i^4},$$

which give the Ricci tensor in all three index positions.

The mixed spatial component has a compact form, on which the shear of Step 10 and the vacuum solution of Step 12 both rest.
With $V = a_1a_2a_3$,

$$\frac{1}{V}\frac{d}{d(ct)}\left(V H_i\right) = \frac{1}{V}\left(V'H_i + VH_i'\right) = H_i\sum_k H_k + \left(u_i - H_i^2\right) = u_i + H_i\sum_{j\neq i}H_j = R^i{}_i,$$

using $V'/V = \sum_k H_k$ and $H_i' = u_i - H_i^2$.
So the mixed spatial Ricci component is the rate of change of $V H_i$, per unit volume.

---

## Step 9. The Ricci scalar and the Einstein tensor

$$R = g^{\mu\nu}R_{\mu\nu} = -R_{tt} + \sum_i \frac{R_{ii}}{a_i^2} = \sum_i u_i + \sum_i\left(u_i + H_i\sum_{j\neq i}H_j\right).$$

The double sum counts each unordered pair twice, so

$$R = 2\left(\sum_i u_i + \sum_{i<j}H_iH_j\right) = 2\left(\frac{a_1''}{a_1} + \frac{a_2''}{a_2} + \frac{a_3''}{a_3} + \frac{a_1'a_2'}{a_1a_2} + \frac{a_1'a_3'}{a_1a_3} + \frac{a_2'a_3'}{a_2a_3}\right),$$

which holds for any three scale factors.

For the Einstein tensor, $G_{\mu\nu} = R_{\mu\nu} - \tfrac{1}{2}Rg_{\mu\nu}$.
In the time slot the second derivatives cancel exactly:

$$G_{tt} = -\sum_i u_i + \tfrac{1}{2}\cdot 2\left(\sum_i u_i + \sum_{i<j}H_iH_j\right) = \sum_{i<j}H_iH_j = \frac{a_1'a_2'}{a_1a_2} + \frac{a_1'a_3'}{a_1a_3} + \frac{a_2'a_3'}{a_2a_3}.$$

That cancellation is the general fact that the time time Einstein equation is a constraint on initial data, carrying first derivatives only.

In a spatial slot, take $i = x$ and let $j,k$ be the other two axes:

$$G^x{}_x = R^x{}_x - \tfrac{1}{2}R = u_1 + H_1(H_2+H_3) - \left(u_1+u_2+u_3\right) - \left(H_1H_2+H_1H_3+H_2H_3\right),$$

and the $H_1$ terms cancel against each other, leaving

$$G^x{}_x = -\left(u_2 + u_3 + H_2H_3\right) = -\left(\frac{a_2''}{a_2} + \frac{a_3''}{a_3} + \frac{a_2'a_3'}{a_2a_3}\right).$$

The equation along one axis involves only the other two, which is the anisotropic generalisation of the FLRW acceleration equation and the reason the three axes are coupled at all.
The fully lowered and doubly raised versions follow by the diagonal rule of Step 8:

$$G_{ii} = a_i^2\,G^i{}_i, \qquad G^{ii} = \frac{G^i{}_i}{a_i^2}, \qquad G^{tt} = G_{tt},$$

which give the Einstein tensor in all three index positions.

---

## Step 10. What the Einstein tensor says: a density and three pressures

Nothing so far has used a field equation.
Imposing one is what turns the three scale factors from arbitrary functions into a cosmology, and the tensors are kept general rather than solved precisely so that any matter model can be inserted here.

Take a comoving perfect fluid, at rest in these coordinates, whose four velocity in this chart is $u^\mu = (c,0,0,0)$ and whose stress energy tensor is

$$T_{\mu\nu} = \left(\rho + \frac{p}{c^2}\right)u_\mu u_\nu + p\,g_{\mu\nu},$$

whose components in this chart are $T_{tt} = \rho c^2$ and $T^i{}_i = p$, since $u_t = -c$ and $u^\mu u_\mu = -c^2$.
An anisotropic universe has no reason to carry a single pressure, so allow three, one per axis, and write the mixed stress energy tensor as

$$T^\mu{}_\nu = \mathrm{diag}\left(-\rho c^2,\, p_x,\, p_y,\, p_z\right),$$

which is the comoving perfect fluid exactly when the three agree.
Then $G_{\mu\nu} = \dfrac{8\pi G}{c^4}T_{\mu\nu}$ gives four equations.

The time equation is the constraint of Step 9:

$$\frac{a_1'a_2'}{a_1a_2} + \frac{a_1'a_3'}{a_1a_3} + \frac{a_2'a_3'}{a_2a_3} = \frac{8\pi G}{c^2}\rho,$$

which, written with dots for $d/dt$ using $H_i = \dot{a}_i/(ca_i)$, is

$$\frac{\dot{a}_1\dot{a}_2}{a_1a_2} + \frac{\dot{a}_1\dot{a}_3}{a_1a_3} + \frac{\dot{a}_2\dot{a}_3}{a_2a_3} = 8\pi G\rho.$$

Setting $a_1 = a_2 = a_3 = a$ turns the left side into $3(\dot{a}/a)^2$ and recovers the flat Friedmann equation exactly.

The three spatial equations are

$$\frac{\ddot{a}_2}{a_2} + \frac{\ddot{a}_3}{a_3} + \frac{\dot{a}_2\dot{a}_3}{a_2a_3} = -\frac{8\pi G}{c^2}p_x,$$

and its two cyclic partners, which for $a_1=a_2=a_3$ collapse to $2\ddot{a}/a + \dot{a}^2/a^2 = -8\pi Gp/c^2$, the flat FLRW acceleration equation.

Two consequences are worth stating because they are what the anisotropy buys.

First, subtract two spatial equations.
By Step 8 the difference of two mixed Ricci components is a total derivative,

$$R^i{}_i - R^j{}_j = \frac{1}{V}\frac{d}{d(ct)}\Big(V\left(H_i - H_j\right)\Big),$$

and $G^i{}_i - G^j{}_j = R^i{}_i - R^j{}_j$ because the trace term is the same in both.
So if the pressure is isotropic, $p_x = p_y = p_z$, then

$$V\left(H_i - H_j\right) = \text{constant},$$

and the shear between any two axes dies off as $1/V$.
An expanding type I universe filled with an ordinary fluid isotropises on its own, which is the quantitative form of the question the whole Bianchi family exists to ask.

Second, the three spatial equations are independent, so three unequal scale factors are matched slot for slot by three unequal pressures, and $T^x{}_x - T^y{}_y$ measures the anisotropy of the matter the way $H_1 - H_2$ measures the anisotropy of the geometry.

---

## Step 11. Why the contraction has to be the standard one

The two conventions for the Ricci tensor differ by a sign, since $R^\alpha{}_{\mu\nu\alpha} = -R^\alpha{}_{\mu\alpha\nu}$ by the antisymmetry of Riemann in its last two indices.
The Einstein tensor, built from Ricci and its trace, flips with it.

For a vacuum solution the choice is invisible, which is why the zeros of Kasner hold under either reading.
Here it is not invisible.
With the standard contraction, Step 9 gives

$$G_{tt} = \sum_{i<j}H_iH_j,$$

which for a universe expanding along every axis is a sum of three positive terms, so $\rho > 0$.
With the other contraction $G_{tt}$ is the negative of that, and the same expanding universe would be filled with matter of negative energy density.
The check is not a matter of taste: a diagonal type I model with the right hand side of the field equations fixed to $8\pi G T_{\mu\nu}/c^4$ has a positive density if and only if the Ricci tensor is contracted the standard way.

Everything other than the Ricci tensor, the Ricci scalar and the Einstein tensor is untouched by the convention.
The Christoffel symbols, the Riemann tensor, the Kretschmann scalar and the geodesic equations do not involve a contraction at all, and the Weyl tensor is defined by removing the traces of Riemann, which are Riemann's own regardless of which sign is called Ricci.

---

## Step 12. The vacuum specialisation, which is Kasner

Set $\rho = 0$ and $p_x = p_y = p_z = 0$.
The field equations become $R_{\mu\nu} = 0$, and Step 8 turns the three spatial ones into

$$\frac{d}{d(ct)}\left(VH_i\right) = V R^i{}_i = 0, \qquad i = 1,2,3,$$

so each $VH_i$ is a constant $q_i$.
Summing, and using $V'/V = \sum_i H_i$,

$$V' = V\sum_i H_i = \sum_i q_i \equiv q,$$

a constant, so the comoving volume is linear in the time:

$$V = q\left(ct - ct_0\right).$$

Two cases.
If $q = 0$ the volume is constant, every $H_i$ is the constant $q_i/V$, and the time equation $\sum_{i<j}H_iH_j = 0$ together with $\sum_i H_i = 0$ forces $\sum_i H_i^2 = (\sum_i H_i)^2 - 2\sum_{i<j}H_iH_j = 0$, hence every $H_i = 0$.
That is Minkowski spacetime in disguise.

Otherwise put the singularity at $t_0 = 0$, so $V = qct$ and

$$H_i = \frac{q_i}{V} = \frac{p_i}{ct}, \qquad p_i \equiv \frac{q_i}{q}, \qquad \sum_i p_i = 1$$

by construction.
Integrating $a_i'/a_i = p_i/(ct)$, which in ordinary time is $\dot{a}_i/a_i = p_i/t$, gives

$$a_i \propto t^{p_i}.$$

The time equation has not been used yet, and it is the second constraint:

$$\sum_{i<j}H_iH_j = \frac{1}{c^2t^2}\sum_{i<j}p_ip_j = 0 \implies \left(\sum_i p_i\right)^2 - \sum_i p_i^2 = 0 \implies \sum_i p_i^2 = 1.$$

So the vacuum type I model is either Minkowski or

$$ds^2 = -c^2dt^2 + t^{2p_1}dx^2 + t^{2p_2}dy^2 + t^{2p_3}dz^2, \qquad \sum_i p_i = \sum_i p_i^2 = 1,$$

which is the Kasner line element, line for line.
The two constraints that accompany the Kasner metric are derived here as the content of the field equations, and the derivation used all four of them: the three spatial equations gave the power law and the sum rule, and the time equation gave the sum of squares.

The relation between the two spacetimes is therefore exact.
Type I with a fluid and Kasner are one family read at two settings of the matter, Kasner is the setting where the matter is gone, and every component of Kasner can be obtained from those of type I by substituting $a_i = t^{p_i}$ and imposing the two constraints.
The substitution is checked in Step 14 on the Kretschmann scalar, the sharpest of the invariants.

---

## Step 13. The Weyl tensor

In four dimensions,

$$C_{\mu\nu\rho\sigma} = R_{\mu\nu\rho\sigma} - \frac{1}{2}\left(g_{\mu\rho}R_{\sigma\nu} - g_{\mu\sigma}R_{\rho\nu} - g_{\nu\rho}R_{\sigma\mu} + g_{\nu\sigma}R_{\rho\mu}\right) + \frac{R}{6}\left(g_{\mu\rho}g_{\sigma\nu} - g_{\mu\sigma}g_{\rho\nu}\right),$$

where $R_{\mu\nu}$ here is the trace of Riemann in the standard sense and $R$ its scalar, which by Step 11 is the same tensor as the Ricci tensor of Step 8.

Unlike in vacuum, there is no shortcut available: the Ricci tensor does not vanish, so Weyl is not Riemann and every component has to be computed.

Take $C_{txtx}$.
The surviving terms are the ones whose metric factors are nonzero, namely $g_{tt}$ and $g_{xx}$:

$$C_{txtx} = R_{txtx} - \frac{1}{2}\left(g_{tt}R_{xx} + g_{xx}R_{tt}\right) + \frac{R}{6}g_{tt}g_{xx}.$$

Substituting Steps 7, 8 and 9 and factoring $a_1^2$ out of the result,

$$C_{txtx} = a_1^2\left[-u_1 + \frac{1}{2}\left(u_1 + H_1(H_2+H_3)\right) + \frac{1}{2}\sum_i u_i - \frac{1}{6}\cdot 2\left(\sum_i u_i + \sum_{i<j}H_iH_j\right)\right],$$

and the bracket simplifies to one sixth of

$$E_1 \equiv \frac{a_2''}{a_2} + \frac{a_3''}{a_3} - 2\frac{a_1''}{a_1} + \frac{a_1'a_2'}{a_1a_2} + \frac{a_1'a_3'}{a_1a_3} - 2\frac{a_2'a_3'}{a_2a_3},$$

with $E_2$ and $E_3$ defined by cycling the labels.
The whole Weyl tensor is built from those three:

$$C_{titi} = \frac{a_i^2}{6}E_i = -C_{tiit} = -C_{itti} = C_{itit}, \qquad C_{ijij} = -\frac{a_i^2a_j^2}{6}E_k = -C_{ijji},$$

where in the second family $k$ is the axis that is neither $i$ nor $j$.
Raising the first index divides by the matching diagonal entry, which gives the Weyl tensor with its first index up:

$$C^t{}_{iti} = -\frac{a_i^2}{6}E_i, \qquad C^i{}_{tti} = -\frac{1}{6}E_i, \qquad C^i{}_{jij} = -\frac{a_j^2}{6}E_k.$$

Three structural facts about $E_i$, each checked in sympy.

First, $E_1 + E_2 + E_3 = 0$, which is the tracefree property of Weyl showing up in the only place this geometry leaves for it.
Two independent functions describe the entire conformal curvature of a type I universe.

Second, rearranged,

$$E_i = \left(\sum_k u_k + \sum_{k<l}H_kH_l\right) - 3\left(u_i + H_jH_k\right),$$

so the Weyl tensor vanishes if and only if the three combinations $u_i + H_jH_k$ agree.
The isotropic case $a_1 = a_2 = a_3$ satisfies that identically, which is the statement that the flat FLRW model is conformally flat, and that is the only place the Weyl tensor of this spacetime is allowed to vanish.

Third, in vacuum the Ricci tensor vanishes and $C_{\mu\nu\rho\sigma} = R_{\mu\nu\rho\sigma}$, which is the fact stated for Kasner.
It is a fact, not a shortcut, and substituting $a_i = t^{p_i}$ on the Kasner circle shows it component by component.
There

$$E_i = \frac{6\,p_i(1-p_i)}{c^2t^2}, \qquad \frac{a_i^2}{6}E_i = \frac{p_i(1-p_i)\,t^{2p_i-2}}{c^2},$$

and the right hand side is exactly the $R_{titi}$ of Kasner, so $C_{titi} = R_{titi}$ there as claimed.
All three $E_i$ vanish together only when every exponent is $0$ or $1$, which on the circle means $(1,0,0)$ up to permutation: the three isolated points where the Kasner metric is flat space in disguise.
All the curvature of an anisotropic vacuum is tidal, and the Weyl tensor is where it lives.

---

## Step 14. The Kretschmann scalar

$K = R_{\mu\nu\rho\sigma}R^{\mu\nu\rho\sigma}$.
The metric is diagonal, so raising all four indices on a component whose index pattern is $(a,b,a,b)$ multiplies it by $g^{aa}g^{bb}g^{aa}g^{bb}$ with no sum, and each coordinate plane contributes four components, $abab$, $abba$, $baab$ and $baba$, all of the same square.
Hence

$$K = 4\sum_{a<b}\left(g^{aa}g^{bb}R_{abab}\right)^2.$$

For a time and space plane, $g^{tt}g^{ii}R_{titi} = (-1)a_i^{-2}(-a_ia_i'') = u_i$.
For a spatial plane, $g^{ii}g^{jj}R_{ijij} = a_i^{-2}a_j^{-2}a_ia_j a_i'a_j' = H_iH_j$.
So

$$K = 4\left[\left(\frac{a_1''}{a_1}\right)^2 + \left(\frac{a_2''}{a_2}\right)^2 + \left(\frac{a_3''}{a_3}\right)^2 + \left(\frac{a_1'a_2'}{a_1a_2}\right)^2 + \left(\frac{a_1'a_3'}{a_1a_3}\right)^2 + \left(\frac{a_2'a_3'}{a_2a_3}\right)^2\right],$$

which holds for every type I model.

Two things follow immediately.
It is a sum of squares, so $K \geq 0$ for every type I model whatever the matter, and it vanishes only when every $u_i$ and every product $H_iH_j$ vanishes, which is flat space.
And it diverges wherever an acceleration or a shear does, which is what happens as the comoving volume $V \to 0$.

The Kasner check.
Substituting $a_i = t^{p_i}$ gives $u_i = p_i(p_i-1)/(c^2t^2)$ and $H_iH_j = p_ip_j/(c^2t^2)$, so

$$K = \frac{4}{c^4t^4}\left[\sum_i p_i^2(p_i-1)^2 + \sum_{i<j}p_i^2p_j^2\right].$$

On the Kasner circle that bracket collapses to $-4p_1p_2p_3$, so

$$K = -\frac{16\,p_1p_2p_3}{c^4t^4},$$

which is the Kretschmann scalar of Kasner, exactly.
The identity was verified in sympy along the rational parametrisation of the circle that `verify_metrics.py` uses, which covers a dense subset of it and therefore proves it there.
It also says something: since $K \geq 0$ always, the product $p_1p_2p_3$ cannot be positive on the circle, which is the algebraic shadow of the fact that one Kasner exponent is always negative.

---

## Step 15. The geodesic equations and their first integrals

The geodesic equation is

$$\ddot{x}^\mu + \Gamma^\mu{}_{\nu\rho}\dot{x}^\nu\dot{x}^\rho = 0,$$

with the dot the derivative with respect to an affine parameter $\lambda$ and each $x^\mu$ a chart coordinate, so $\dot{t}$ means $d(ct)/d\lambda$ and matches the Christoffel symbols of Step 4 slot for slot.

The time equation collects the three $\Gamma^t{}_{ii}$:

$$\ddot{t} + a_1a_1'\dot{x}^2 + a_2a_2'\dot{y}^2 + a_3a_3'\dot{z}^2 = 0.$$

Each spatial equation collects the two orderings of $\Gamma^i{}_{ti}$, which is where the factor of two comes from:

$$\ddot{x} + 2\frac{a_1'}{a_1}\dot{t}\dot{x} = 0, \qquad \ddot{y} + 2\frac{a_2'}{a_2}\dot{t}\dot{y} = 0, \qquad \ddot{z} + 2\frac{a_3'}{a_3}\dot{t}\dot{z} = 0.$$

These are the four geodesic equations.

Each spatial equation integrates once.
Since $d/d\lambda = \dot{t}\,d/d(ct)$ acting on a function of the time alone,

$$\frac{d}{d\lambda}\left(a_i^2\dot{x}^i\right) = 2a_ia_i'\dot{t}\dot{x}^i + a_i^2\ddot{x}^i = a_i^2\left(\ddot{x}^i + 2\frac{a_i'}{a_i}\dot{t}\dot{x}^i\right) = 0,$$

so

$$a_i^2\dot{x}^i = \text{constant}, \qquad i = 1,2,3,$$

three conserved momenta, one per axis, which is exactly what the three translation symmetries of the type I slices are supposed to give by Noether's theorem.
A particle's comoving velocity decays as $a_i^{-2}$ along each axis separately, and the anisotropy of the expansion is therefore imprinted on the momentum distribution of anything moving through it.

The comoving worldlines $x = y = z = \text{const}$ are geodesics: the three spatial equations are satisfied trivially, and the time equation reduces to $\ddot{t} = 0$, so $t$ is an affine parameter and, up to a constant factor, the proper time of a comoving observer.
That is why the coordinate $t$ can be called cosmic time at all.

---

## Step 16. Every equation is dimensionally consistent

In the chart $x^0 = ct$ every coordinate is a length, since $t$ is a time and the chart carries $ct$.
The scale factors are dimensionless by declaration, which is what makes the line element balance: $a_i^2dx^2$ carries $L^2$ and so does $c^2dt^2$.

With every chart coordinate a length, a component of the metric carries nothing, a Christoffel symbol carries $1/L$, a Riemann, Ricci, Einstein or Weyl component carries $1/L^2$, and the Kretschmann scalar carries $1/L^4$.
The primes supply exactly that: $a_i' = c^{-1}da_i/dt$ carries $1/L$ and $a_i''$ carries $1/L^2$.

Run down the list.
$\Gamma^t{}_{xx} = a_1a_1'$ carries $1/L$.
$\Gamma^x{}_{tx} = a_1'/a_1$ carries $1/L$.
$R^t{}_{xtx} = a_1a_1''$ carries $1/L^2$, and so does $R^x{}_{ttx} = a_1''/a_1$ and $R^x{}_{yxy} = a_2a_1'a_2'/a_1$, the last as a product of two rates.
$R_{xx}$, $G_{xx}$ and $C_{txtx}$ all carry $1/L^2$, since the scale factor prefactors are dimensionless.
$K$ is a sum of squares of quantities carrying $1/L^2$ and carries $1/L^4$.

A geodesic equation is measured against its own second derivative, and the dots in it are chart velocities, so $\dot{x}^\mu$ carries $L/\lambda$ and $\ddot{x}^\mu$ carries $L/\lambda^2$.
In the time equation, $\ddot{t}$ carries $L/\lambda^2$ and $a_1a_1'\dot{x}^2$ carries $(1/L)(L/\lambda)^2 = L/\lambda^2$ as well.
In a spatial equation, $\ddot{x}$ and $2(a_1'/a_1)\dot{t}\dot{x}$ balance the same way.

The dimensional pass in `verify_metrics.py` does this term by term over every expression, and every term of type I passes it with no exceptions.

---

## Step 17. The components, and the declaration the checker needs

In the type I Cartesian chart the spacetime has the line element, the metric and its inverse, nine Christoffel symbols in each of two index positions, twenty four Riemann components in each of two index positions, the Ricci tensor in three index positions, the Ricci scalar, the Kretschmann scalar, the Einstein tensor in three index positions, the Weyl tensor in twenty four components in each of two index positions, and the four geodesic equations.

No tensor vanishes.
The Ricci tensor, the Ricci scalar and the Einstein tensor are nonzero because no field equation has been imposed, and the Weyl tensor is nonzero because the model is not conformally flat unless the three scale factors conspire as Step 13 describes.
The two specialisations each remove a different part of the curvature: the Ricci and Einstein tensors vanish in vacuum, which is Kasner, and the Weyl tensor vanishes under isotropy, which is flat FLRW.
No setting of the scale factors empties both at once without emptying the spacetime into Minkowski, since a vanishing Ricci tensor and a vanishing Weyl tensor together leave no Riemann tensor at all.

`verify_metrics.py` checks a system only if it is declared in `DIMENSIONS`, and the declaration for type I is

```python
("bianchi", "type_i_cartesian"): {
    "t": "T", "x": "L", "y": "L", "z": "L", "a_1": "1", "a_2": "1", "a_3": "1",
},
```

which says that $t$ is the coordinate the chart multiplies by $c$, that the three comoving coordinates are lengths, and that the three scale factors are dimensionless.
No `PARAMETER_RELATIONS` line is needed, because unlike Kasner type I constrains nothing: every component holds for arbitrary scale factors.

Two facts matter to whoever lands that declaration.

The reader in `verify_metrics.py` gives a declared function of a time coordinate a prime and a dot that both mean $d/d(ct)$, dividing the bare derivative by $c$, which is the convention of Step 1 and is why no component of type I carries an explicit $c$.
The parser reaches a name like `a_1` through its prime spelling only, since its dot spelling accepts letters alone, so the components are written with primes throughout.

When type I was first checked, the checker contracted the Ricci tensor on the last lower index, the older convention, and against that reading it reported twenty five disagreements: the four components of the Ricci tensor in each of its three index positions, the same for the Einstein tensor, and the Ricci scalar, every one of them a pure sign.
None was an error in the components, which are written with the standard contraction the captain settled on 18 September 2026.
The checker was brought onto that contraction in the pass that flipped the sign of the older Ricci and Einstein tensors, and all twenty five went to zero with no component of type I changed.
