# The pp-wave plane gravitational waves

We work the pp-wave in two charts, the Brinkmann chart with an arbitrary profile and the exact plane wave chart.
We derive every component in order, from the line element down to the geodesic equations.
We leave nothing as an exercise and assert nothing we do not compute.

The companion script `_tools/derivations/verify_metrics.py` does the same work in sympy and compares the result component by component, so the algebra is checkable by hand and by machine independently.

Two things set the pp-wave apart from the spacetimes worked before it.
Its metric carries an arbitrary function rather than a parameter, so the line element describes a whole family at once and the field equations become a condition on that function rather than a statement that some component vanishes.
That condition appears in Step 9, and it is the one thing about this family worth remembering: a pp-wave is a vacuum exactly where its profile is a harmonic function of the two transverse coordinates.
And its curvature is invisible to every scalar built from it.
The Ricci scalar and the Kretschmann scalar are zero for every profile, vacuum or not, while the Weyl tensor is not; the reason is in Step 12, and that fact is the whole of the argument in the history of the pp-wave, in which the waves were nearly reasoned out of existence on paper and were put back by asking what an instrument would measure.

---

## Step 1. Conventions

The signature is $(-,+,+,+)$ and the Christoffel symbols are those of the Levi-Civita connection,

$$\Gamma^\mu_{\nu\rho} = \tfrac{1}{2}g^{\mu\alpha}\left(\partial_\nu g_{\alpha\rho} + \partial_\rho g_{\alpha\nu} - \partial_\alpha g_{\nu\rho}\right),$$

with the lowered symbol $\Gamma_{\mu\nu\rho} = g_{\mu\alpha}\Gamma^\alpha_{\nu\rho}$.
The Riemann tensor is

$$R^\mu{}_{\nu\rho\sigma} = \partial_\rho \Gamma^\mu_{\nu\sigma} - \partial_\sigma \Gamma^\mu_{\nu\rho} + \Gamma^\mu_{\rho\lambda}\Gamma^\lambda_{\nu\sigma} - \Gamma^\mu_{\sigma\lambda}\Gamma^\lambda_{\nu\rho},$$

and the Riemann components follow it.

The Ricci tensor is contracted on the middle index,

$$R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu},$$

which is the standard contraction.
The choice matters here, because the Brinkmann chart carries an arbitrary profile and its Ricci tensor does not vanish.
Contracting on the last lower index instead, $R^\alpha{}_{\mu\nu\alpha}$, gives the negative of everything in Step 9 and Step 10, and the sign carries through to the field equations in Step 10.
The vacuum condition itself is the statement that one component is zero, so it is the same equation on either convention, as at the end of Step 9.

The Weyl tensor is built by removing the traces of Riemann, and those traces are this same contraction, $S_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu} = R_{\mu\nu}$.
So the trace removed in Step 11 is the Ricci tensor itself, with no sign to carry between them.

The chart, here and for every other spacetime, is the one whose time coordinate is $x^0 = cT$, so here it is $x^0 = cu$.
The index is written with the bare letter $u$, but the component written against it is a component of the chart whose zeroth coordinate is $cu$ and therefore carries a length.
Throughout, write

$$U \equiv cu$$

for that chart coordinate, so that the component written with an index $u$ is the component computed with an index $U$.

This has one consequence worth stating before any algebra is done.
A partial derivative in the components is taken with respect to the chart coordinate, not with respect to the bare one:

$$\partial_u H \equiv \frac{\partial H}{\partial (cu)} = \frac{1}{c}\frac{\partial H}{\partial u},$$

which is the same reading already given to the dot in Vaidya's $\dot{m}$ and the prime in the conformal FRW chart.
The transverse derivatives $\partial_x$ and $\partial_y$ need no such factor, because $x$ and $y$ are their own chart coordinates.
The same reading applies to the dots in the geodesic equations, where $\dot{u}$ means $d(cu)/d\lambda$.
It is the reading that makes every term of every equation carry the dimensions of its left hand side, and we check it term by term in Step 15.

We work directly in the chart $x^0 = cu$ from Step 2 to Step 15, so the components we compute need no conversion at the end.

---

## Step 2. The line element

In Brinkmann form,

$$ds^2 = H(u,x,y)\,c^2du^2 - 2c\,du\,dv + dx^2 + dy^2.$$

The coordinate $u$ is the retarded time labelling the wave fronts, $v$ is the coordinate along the rays, and $x$ and $y$ are the two transverse directions.
The profile $H$ is an arbitrary function of $u$, $x$ and $y$.
The one thing asked of it is that it does not depend on $v$, and by Step 6 this is exactly what makes the rays parallel.

In the chart coordinate $U = cu$ the line element is

$$ds^2 = H\,dU^2 - 2\,dU\,dv + dx^2 + dy^2,$$

which carries no $c$ at all.
Every coordinate of the chart is a length, $H$ is dimensionless, and $v$ carries a length because it is multiplied by $c\,du$ rather than by $c^2du^2$.
That is the same arrangement the outgoing Vaidya chart has, where the mixed term is $-2c\,du\,dr$ and $r$ is a length, and it is not the arrangement of the Minkowski double null chart, whose mixed term is $-c^2du\,dv$ and whose $u$ and $v$ are both times.

Surfaces of constant $u$ are null: the vector $\partial_\mu U$ has $g^{\mu\nu}\partial_\mu U\,\partial_\nu U = g^{UU} = 0$ by Step 3.
They are the wave fronts, and the whole solution is the statement that the profile is carried along them unchanged.

---

## Step 3. The metric matrix, its determinant and its inverse

In the order $(U, v, x, y)$ the line element gives

$$g_{\mu\nu} = \begin{pmatrix} H & -1 & 0 & 0 \\ -1 & 0 & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1\end{pmatrix},$$

which gives the metric components $g_{uu} = H$, $g_{uv} = g_{vu} = -1$ and $g_{xx} = g_{yy} = 1$.
The mixed term of the line element is $-2\,dU\,dv$ and contributes half of its coefficient to each of the two off diagonal slots, which is where the $-1$ comes from.

The determinant is the product of the two blocks,

$$\det g = \left(H\cdot 0 - (-1)(-1)\right)\cdot 1 = -1,$$

so $\sqrt{-g} = 1$ for every profile.
The chart is unimodular: the volume element is the flat one no matter how strong the wave is, which already says that no amount of profile can make the coordinate volume of a region change.

The determinant is negative and the transverse block is positive definite, so the $(U,v)$ block has one positive and one negative eigenvalue and the signature is $(-,+,+,+)$ as claimed.

Inverting the $(U,v)$ block,

$$\begin{pmatrix} H & -1 \\ -1 & 0\end{pmatrix}^{-1} = \frac{1}{-1}\begin{pmatrix} 0 & 1 \\ 1 & H\end{pmatrix} = \begin{pmatrix} 0 & -1 \\ -1 & -H\end{pmatrix},$$

so

$$g^{\mu\nu} = \begin{pmatrix} 0 & -1 & 0 & 0 \\ -1 & -H & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1\end{pmatrix},$$

which gives the inverse metric components $g^{uv} = g^{vu} = -1$, $g^{vv} = -H$ and $g^{xx} = g^{yy} = 1$, together with $g^{uu} = 0$.
The check is immediate: $g^{U\alpha}g_{\alpha U} = (-1)(-1) = 1$, $g^{U\alpha}g_{\alpha v} = (-1)\cdot 0 = 0$, and $g^{v\alpha}g_{\alpha v} = (-1)(-1) + (-H)\cdot 0 = 1$.

Two facts from this matrix recur throughout and are worth naming now.
The first is that $g^{UU} = 0$, so an upper $U$ index can never be produced from a lower one: $g^{\mu U}$ is nonzero only for $\mu = v$.
The second is that $H$ appears in the inverse in one slot only, $g^{vv}$, and every tensor computed from it turns out to have no $v$ index at all, so that slot never gets used.

---

## Step 4. The lowered Christoffel symbols

With $\Gamma_{\mu\nu\rho} = \tfrac{1}{2}\left(\partial_\nu g_{\mu\rho} + \partial_\rho g_{\mu\nu} - \partial_\mu g_{\nu\rho}\right)$, the only entry of the metric that is not constant is $g_{UU} = H$, and $H$ does not depend on $v$.
So the only nonvanishing first derivatives of the metric are

$$\partial_U g_{UU} = \partial_u H, \qquad \partial_x g_{UU} = \partial_x H, \qquad \partial_y g_{UU} = \partial_y H,$$

where the first is a derivative with respect to $U = cu$ and is written $\partial_u H$ by the convention of Step 1.

A symbol is therefore nonzero only if two of its three indices are the $U$ and $U$ of $g_{UU}$, the remaining one naming the direction differentiated in.
There are two ways to place them.

**Both metric indices in the first and second slots.**
Taking $\mu = \nu = U$ and $\rho$ free,

$$\Gamma_{UU\rho} = \tfrac{1}{2}\left(\partial_U g_{U\rho} + \partial_\rho g_{UU} - \partial_U g_{U\rho}\right) = \tfrac{1}{2}\partial_\rho H,$$

the outer two terms cancelling whatever $\rho$ is.
This gives $\Gamma_{UUU} = \tfrac{1}{2}\partial_u H$, $\Gamma_{UUx} = \tfrac{1}{2}\partial_x H$ and $\Gamma_{UUy} = \tfrac{1}{2}\partial_y H$, and $\Gamma_{UUv} = 0$ because $\partial_v H = 0$.
The symbol is symmetric in its last two indices, so $\Gamma_{UxU} = \Gamma_{UUx}$ and $\Gamma_{UyU} = \Gamma_{UUy}$ as well.

**Both metric indices in the last two slots.**
Taking $\nu = \rho = U$ and $\mu$ free,

$$\Gamma_{\mu UU} = \tfrac{1}{2}\left(2\,\partial_U g_{\mu U} - \partial_\mu g_{UU}\right) = -\tfrac{1}{2}\partial_\mu H,$$

since $g_{\mu U}$ is constant for every $\mu$.
This gives $\Gamma_{xUU} = -\tfrac{1}{2}\partial_x H$ and $\Gamma_{yUU} = -\tfrac{1}{2}\partial_y H$, while $\Gamma_{vUU} = -\tfrac{1}{2}\partial_v H = 0$ and $\Gamma_{UUU}$ is the one already found.

Those seven are the whole of the Christoffel symbols with every index lowered:

$$\Gamma_{uuu} = \tfrac{1}{2}\partial_u H, \qquad \Gamma_{uux} = \Gamma_{uxu} = \tfrac{1}{2}\partial_x H, \qquad \Gamma_{uuy} = \Gamma_{uyu} = \tfrac{1}{2}\partial_y H,$$

$$\Gamma_{xuu} = -\tfrac{1}{2}\partial_x H, \qquad \Gamma_{yuu} = -\tfrac{1}{2}\partial_y H.$$

Every one of them carries at least two $u$ indices, and none of them carries a $v$.

---

## Step 5. The Christoffel symbols with an upper index

Raise the first index with $\Gamma^\mu{}_{\nu\rho} = g^{\mu\alpha}\Gamma_{\alpha\nu\rho}$, taking the rows of the inverse metric one at a time.

**The $U$ row.**
$g^{U\alpha}$ is nonzero only for $\alpha = v$, and every lowered symbol with a first index $v$ vanishes by Step 4.
So

$$\Gamma^u{}_{\nu\rho} = 0$$

identically, for every profile.
No Christoffel symbol in this chart carries an upper $u$ index.

**The $v$ row.**
$g^{v\alpha}$ is nonzero for $\alpha = U$ and $\alpha = v$, so $\Gamma^v{}_{\nu\rho} = -\Gamma_{U\nu\rho} - H\Gamma_{v\nu\rho} = -\Gamma_{U\nu\rho}$, the second term dropping because every lowered symbol with a first index $v$ vanishes.
This gives

$$\Gamma^v{}_{uu} = -\tfrac{1}{2}\partial_u H, \qquad \Gamma^v{}_{ux} = \Gamma^v{}_{xu} = -\tfrac{1}{2}\partial_x H, \qquad \Gamma^v{}_{uy} = \Gamma^v{}_{yu} = -\tfrac{1}{2}\partial_y H.$$

**The transverse rows.**
$g^{xx} = g^{yy} = 1$, so $\Gamma^x{}_{\nu\rho} = \Gamma_{x\nu\rho}$ and $\Gamma^y{}_{\nu\rho} = \Gamma_{y\nu\rho}$, giving

$$\Gamma^x{}_{uu} = -\tfrac{1}{2}\partial_x H, \qquad \Gamma^y{}_{uu} = -\tfrac{1}{2}\partial_y H.$$

Those seven are the whole of the Christoffel symbols with the first index up.
A symbol with an upper $v$ index and a symbol with an upper transverse index carry the same $-\tfrac{1}{2}\partial_a H$: the first moves a ray along itself, the second bends a transverse direction, and they are the same number because the wave shears the transverse plane and pushes along the ray by the same amount.

---

## Step 6. The parallel ray, and why the fronts are planes

Let $k$ be the vector field $\partial_v$, so that $k^\mu = \delta^\mu_v$.
Lowering it with the metric of Step 3,

$$k_\mu = g_{\mu v} = (-1, 0, 0, 0),$$

so $k_\mu = -\partial_\mu U = -\partial_\mu(cu)$, the one form normal to the wave fronts.
It is null, $k^\mu k_\mu = g_{vv} = 0$.

It is also covariantly constant.
Since $k_\mu$ has constant components,

$$\nabla_\mu k_\nu = \partial_\mu k_\nu - \Gamma^\alpha{}_{\mu\nu}k_\alpha = -\Gamma^U{}_{\mu\nu}k_U = \Gamma^u{}_{\mu\nu} = 0$$

by the $U$ row of Step 5.
This is the defining property of the family and the reason for its name.
A covariantly constant null vector field means the rays, the curves with $u$, $x$ and $y$ held constant, are null geodesics with $v$ an affine parameter, and that they are parallel to one another everywhere rather than only at a point; the wave fronts they are orthogonal to are then flat planes, since the transverse metric $dx^2 + dy^2$ is the flat one for every $u$.
Plane fronted waves with parallel rays is what pp abbreviates.

The profile's independence of the ray coordinate enters at Step 4: $H$ has no $v$ in it, so $g_{\mu\nu}$ has no $v$ in it, so no derivative of the metric can produce the $\Gamma_{v\nu\rho}$ that would have spoiled the $U$ row.

---

## Step 7. The Riemann tensor, and why it is linear in the profile

The quadratic terms of the Riemann tensor cancel identically here, which is what makes every curvature component a single second derivative of $H$ rather than a sum of squares.
The argument is a count over the seven symbols of Step 5.

Consider $\Gamma^\mu{}_{\rho\lambda}\Gamma^\lambda{}_{\nu\sigma}$.
The second factor is nonzero only if its upper index $\lambda$ is $v$, $x$ or $y$.
In the first factor $\lambda$ sits in a lower slot, and the lower slots of the nonzero symbols hold only $U$, $x$ and $y$, never $v$.
So $\lambda$ is $x$ or $y$, which forces the second factor to be $\Gamma^x{}_{UU}$ or $\Gamma^y{}_{UU}$ and therefore forces $\nu = \sigma = U$.
With $\lambda$ transverse, the first factor must be $\Gamma^v{}_{Ux}$ or $\Gamma^v{}_{Uy}$, which forces $\mu = v$ and $\rho = U$.
The only surviving product is therefore the one with $\rho = \sigma = U$, and the other quadratic term $\Gamma^\mu{}_{\sigma\lambda}\Gamma^\lambda{}_{\nu\rho}$ is the same expression with $\rho$ and $\sigma$ exchanged.
With $\rho = \sigma$ the two are equal and the difference is zero.

So for every profile

$$R^\mu{}_{\nu\rho\sigma} = \partial_\rho \Gamma^\mu{}_{\nu\sigma} - \partial_\sigma \Gamma^\mu{}_{\nu\rho},$$

and the Riemann tensor is linear in $H$.
This is why pp-waves superpose, which almost nothing else in general relativity does.

Write $a$ and $b$ for transverse indices, each running over $x$ and $y$.
There are three rows to take.

**The $v$ row with a transverse second index.**
From $\Gamma^v{}_{aU} = -\tfrac{1}{2}\partial_a H$ and $\Gamma^v{}_{ab} = 0$,

$$R^v{}_{aUb} = \partial_U\Gamma^v{}_{ab} - \partial_b\Gamma^v{}_{aU} = \tfrac{1}{2}\partial_a\partial_b H, \qquad R^v{}_{abU} = -\tfrac{1}{2}\partial_a\partial_b H.$$

**The $v$ row with $U$ as the second index.**
Here $\Gamma^v{}_{U\sigma} = -\tfrac{1}{2}\partial_\sigma H$ for every $\sigma$ at once, including $\sigma = U$ and $\sigma = v$, so

$$R^v{}_{U\rho\sigma} = -\tfrac{1}{2}\left(\partial_\rho\partial_\sigma H - \partial_\sigma\partial_\rho H\right) = 0.$$

Mixed partial derivatives commute, so this whole row vanishes, and with it every component that would have carried a derivative of $H$ with respect to $u$.
No component of the curvature knows how fast the profile changes in retarded time.

**The transverse rows.**
The only symbol with an upper transverse index is $\Gamma^a{}_{UU} = -\tfrac{1}{2}\partial_a H$, so the second index must be $U$, and

$$R^a{}_{UUb} = \partial_U\Gamma^a{}_{Ub} - \partial_b\Gamma^a{}_{UU} = \tfrac{1}{2}\partial_a\partial_b H, \qquad R^a{}_{UbU} = -\tfrac{1}{2}\partial_a\partial_b H.$$

Those sixteen components are the whole of the Riemann tensor with the first index up: four choices of the transverse pair $(a,b)$ in each of the four slots

$$R^v{}_{aub} = \tfrac{1}{2}\partial_a\partial_b H, \quad R^v{}_{abu} = -\tfrac{1}{2}\partial_a\partial_b H, \quad R^a{}_{uub} = \tfrac{1}{2}\partial_a\partial_b H, \quad R^a{}_{ubu} = -\tfrac{1}{2}\partial_a\partial_b H,$$

and nothing else is nonzero.
The whole curvature is the transverse Hessian of the profile,

$$\mathcal{H}_{ab} \equiv \partial_a\partial_b H,$$

a symmetric two by two matrix of functions of $u$, $x$ and $y$.

---

## Step 8. The lowered Riemann tensor

Lower the first index with $R_{\mu\nu\rho\sigma} = g_{\mu\alpha}R^\alpha{}_{\nu\rho\sigma}$, again by rows of the metric.

The $v$ row of $g_{\mu\nu}$ has only $g_{vU} = -1$, so $R_{v\nu\rho\sigma} = -R^u{}_{\nu\rho\sigma} = 0$: the lowered tensor has no component with a first index $v$, and by the pair symmetries it has no component with a $v$ anywhere.
The $U$ row has $g_{UU} = H$ and $g_{Uv} = -1$, so $R_{U\nu\rho\sigma} = H R^u{}_{\nu\rho\sigma} - R^v{}_{\nu\rho\sigma} = -R^v{}_{\nu\rho\sigma}$, the profile dropping out with the vanishing upper $u$ row.
The transverse rows are unchanged.

The result is one formula,

$$R_{uaub} = -\tfrac{1}{2}\partial_a\partial_b H,$$

together with everything the symmetries $R_{\mu\nu\rho\sigma} = -R_{\nu\mu\rho\sigma} = -R_{\mu\nu\sigma\rho} = R_{\rho\sigma\mu\nu}$ generate from it, which is the sixteen components of the Riemann tensor with every index lowered.
The pair symmetry is visible in it: $R_{uaub}$ and $R_{aubu}$ are the same number, as they must be.

Two properties of the lowered Riemann tensor carry everything that follows.
Every nonzero component has exactly two $u$ indices, one in each antisymmetric pair, and two transverse ones.
And no nonzero component has a $v$ index in any slot.

---

## Step 9. The Ricci tensor, and where the harmonic condition comes from

Contract on the middle index, $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$, using the Riemann tensor of Step 7 with the first index up.

For the contraction to be nonzero the component must have its first and third indices equal.
The $v$ row of Step 7 has $R^v{}_{a\rho\sigma}$ with $\rho$ drawn from $U$ and the transverse directions, so its third index is never $v$ and it contributes nothing.
The transverse rows have $R^a{}_{U\rho\nu}$, whose third index equals the first when $\rho = a$, and by Step 7, $R^a{}_{Ua\nu}$ is nonzero only for $\nu = U$, where it is $R^a{}_{UaU} = -\tfrac{1}{2}\partial_a\partial_a H$ with no sum.
The $u$ row is zero.

So the Ricci tensor has one component,

$$R_{uu} = R^x{}_{uxu} + R^y{}_{uyu} = -\tfrac{1}{2}\left(\partial_x^2 H + \partial_y^2 H\right),$$

and every other component vanishes.
On the other contraction it is the negative of this, $+\tfrac{1}{2}(\partial_x^2 H + \partial_y^2 H)$, and the two agree on when it is zero.

Write $\Delta_\perp \equiv \partial_x^2 + \partial_y^2$ for the Laplacian of the transverse plane, so that the whole Ricci tensor of a pp-wave is

$$R_{\mu\nu} = -\tfrac{1}{2}\left(\Delta_\perp H\right)k_\mu k_\nu,$$

with $k_\mu$ the covariantly constant null one form of Step 6, whose only component is $k_u = -1$.
The Ricci tensor is null, in the same sense Vaidya's is: it is a multiple of $k_\mu k_\nu$, and the multiple is minus half the trace of the transverse Hessian.

The vacuum field equations are $R_{\mu\nu} = 0$, and with one component that is one equation:

$$\boxed{\ \partial_x^2 H + \partial_y^2 H = 0.\ }$$

A pp-wave is a vacuum exactly where its profile is a harmonic function of the two transverse coordinates, at each moment of retarded time, with the dependence on $u$ left completely free.
The equation is linear, which is the superposition of Step 7 seen again, and it is an equation in two dimensions, so its solutions are the real parts of holomorphic functions of $\zeta = x + iy$ with $u$ carried along as a parameter.

One consequence is worth drawing out, because it explains why the exact plane wave chart looks the way it does.
A function harmonic on the whole transverse plane and bounded there is constant, and a constant profile has vanishing Hessian and therefore vanishing curvature by Step 7.
So a nontrivial vacuum pp-wave has to grow without bound in the transverse directions, or else be singular somewhere in the plane, or be defined on less than the whole plane.
These spacetimes are never asymptotically flat sideways, and the exact plane waves of Step 14 are quadratic for exactly this reason.

---

## Step 10. The Ricci scalar, the Einstein tensor and the matter

The Ricci scalar is $R = g^{\mu\nu}R_{\mu\nu}$, and the Ricci tensor has only a $uu$ component, so only $g^{uu}$ can reach it.
By Step 3 that entry of the inverse metric is zero, so

$$R = g^{uu}R_{uu} = 0$$

for every profile, whether or not the spacetime is a vacuum.
This is not the vacuum condition in disguise; it is the null structure, and its general version is in Step 12.

With $R = 0$ the Einstein tensor is the Ricci tensor,

$$G_{\mu\nu} = R_{\mu\nu} - \tfrac{1}{2}Rg_{\mu\nu} = R_{\mu\nu}, \qquad G_{uu} = -\tfrac{1}{2}\left(\partial_x^2 H + \partial_y^2 H\right),$$

which is the whole of the Einstein tensor with both indices down.
Raising indices with Step 3 sends a lower $u$ to an upper $v$ and multiplies by $-1$ each time, since $g^{vu} = -1$, so

$$G^u{}_v = R^u{}_v = 0, \qquad G^v{}_u = R^v{}_u = \tfrac{1}{2}\left(\partial_x^2 H + \partial_y^2 H\right), \qquad G^{vv} = R^{vv} = -\tfrac{1}{2}\left(\partial_x^2 H + \partial_y^2 H\right),$$

which are the Einstein tensor with mixed indices and with both indices up.

Where the profile is not harmonic the spacetime is not empty, and the stress energy it needs is pure radiation.
With the standard contraction the field equations carry their usual sign, $G_{\mu\nu} = 8\pi G\,T_{\mu\nu}/c^4$, so

$$T_{\mu\nu} = -\frac{c^4}{16\pi G}\left(\Delta_\perp H\right)k_\mu k_\nu,$$

a null fluid streaming along the same rays the wave does.
An observer with four velocity $\xi$ measures the energy density $T_{\mu\nu}\xi^\mu\xi^\nu = -\frac{c^4}{16\pi G}\left(\Delta_\perp H\right)(\xi^u)^2$, so the energy conditions hold exactly where $\Delta_\perp H \le 0$, and the profile of a physical non vacuum pp-wave is superharmonic rather than harmonic.
Electromagnetic plane waves and beams of null dust are the usual sources, and the vacuum case is the boundary between them and the profiles no matter can support.

---

## Step 11. The Weyl tensor, and why it is not a copy of Riemann

The Weyl tensor in four dimensions is Riemann with its traces removed,

$$C_{\mu\nu\rho\sigma} = R_{\mu\nu\rho\sigma} - \tfrac{1}{2}\left(g_{\mu\rho}S_{\sigma\nu} - g_{\mu\sigma}S_{\rho\nu} - g_{\nu\rho}S_{\sigma\mu} + g_{\nu\sigma}S_{\rho\mu}\right) + \tfrac{S}{6}\left(g_{\mu\rho}g_{\sigma\nu} - g_{\mu\sigma}g_{\rho\nu}\right),$$

built from $S_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$, which by Step 1 is the Ricci tensor itself here, and from $S = g^{\mu\nu}S_{\mu\nu} = R = 0$.
So the trace term drops and the only surviving piece of $S_{\mu\nu}$ is

$$S_{uu} = R_{uu} = -\tfrac{1}{2}\Delta_\perp H.$$

Take the family Riemann populates, $C_{uaub}$ with $a$ and $b$ transverse:

$$C_{uaub} = R_{uaub} - \tfrac{1}{2}\left(g_{uu}S_{ba} - g_{ub}S_{ua} - g_{au}S_{bu} + g_{ab}S_{uu}\right).$$

Of the four correction terms, $S_{ba} = 0$ kills the first, and $g_{ub} = g_{au} = 0$ kill the middle two, leaving $g_{ab} = \delta_{ab}$ times $S_{uu}$.
So

$$C_{uaub} = -\tfrac{1}{2}\partial_a\partial_b H + \tfrac{1}{4}\delta_{ab}\Delta_\perp H = -\tfrac{1}{2}\left(\partial_a\partial_b - \tfrac{1}{2}\delta_{ab}\Delta_\perp\right)H.$$

The Weyl tensor is the trace free part of the transverse Hessian, exactly as the Ricci tensor was its trace.
Written out, that is

$$C_{uxux} = -C_{uyuy} = \tfrac{1}{4}\left(\partial_y^2 H - \partial_x^2 H\right), \qquad C_{uxuy} = -\tfrac{1}{2}\partial_x\partial_y H,$$

and the sixteen components of the Weyl tensor with every index lowered are these and what the symmetries of Step 8 generate.
Raising the first index sends $u$ to $v$ and flips the sign, giving the Weyl tensor with the first index up.

Nothing new appears in a slot where Riemann vanished, and that has to be checked rather than assumed, because it is exactly where a Weyl tensor differs from a Riemann tensor in general.
A correction term is nonzero only when the two indices landing on $S$ are both $u$ and the two landing on the metric pick out a nonzero entry, which by Step 3 means $g_{uu}$, $g_{uv}$, $g_{vu}$, $g_{xx}$ or $g_{yy}$.
The transverse choices are the family just computed.
The choices $g_{uu}$, $g_{uv}$ and $g_{vu}$ put a $u$ in the same antisymmetric pair as one of the two $u$ indices on $S$, and the four correction terms then cancel in pairs: $C_{uuvu}$, for instance, has $-\tfrac{1}{2}(g_{uv}S_{uu} - g_{uu}S_{vu} - g_{uv}S_{uu} + g_{uu}S_{vu}) = 0$.
So the Weyl tensor has the same sixteen nonzero slots as the Riemann tensor and no others.

The curvature is of Petrov type N, and that can be read off the components without any of the Newman-Penrose apparatus.
No nonzero component of $C$ carries a $v$ index, so contracting any slot with the ray $k^\mu = \delta^\mu_v$ of Step 6 gives zero:

$$C_{\mu\nu\rho\sigma}k^\sigma = 0,$$

and by the symmetries the same holds in every slot.
A null direction that annihilates the Weyl tensor outright is a principal null direction of multiplicity four, which is the definition of type N, the algebraically special type of pure radiation.
The two independent numbers left in it, $\tfrac{1}{2}(\partial_x^2 - \partial_y^2)H$ and $\partial_x\partial_y H$, are the two polarisations.

In a vacuum $\Delta_\perp H = 0$ and the correction term vanishes, so there the Weyl tensor is the Riemann tensor exactly.
That is why the Weyl and Riemann tensors of the exact plane wave chart of Step 14 are identical, and it is a genuine equality rather than the copying mistake undone in `_tools/derivations/weyl.md`: it holds because that chart is a vacuum, and sympy confirms it from the line element without being told.

---

## Step 12. Why every curvature scalar vanishes while the curvature does not

The Kretschmann scalar is $K = R_{\mu\nu\rho\sigma}R^{\mu\nu\rho\sigma}$, and it is zero for every profile.
The reason is the two properties of Step 8, and it is worth doing slowly because the history of the pp-wave turns on it.

Raise all four indices of the second factor.
For a term to contribute, the lowered tensor inside it must be nonzero, so by Step 8 each of its pairs holds one $u$ index and one transverse index.
Raising a transverse index leaves it transverse, since $g^{xx} = g^{yy} = 1$ and the transverse rows of the inverse metric have nothing else in them.
Raising a $u$ index goes through $g^{\mu u}$, which by Step 3 is nonzero only for $\mu = v$.
So every nonzero component of $R^{\mu\nu\rho\sigma}$ carries two $v$ indices, one in each pair, and two transverse ones.

But the first factor $R_{\mu\nu\rho\sigma}$ has no nonzero component with a $v$ index anywhere.
Every term of the sum therefore has a zero in it, and

$$K = 0$$

identically, for every profile, vacuum or not.
The same argument applied to the Weyl tensor, whose components sit in the same slots by Step 11, gives $C_{\mu\nu\rho\sigma}C^{\mu\nu\rho\sigma} = 0$, and applied to the Ricci tensor gives $R_{\mu\nu}R^{\mu\nu} = R_{vv}R^{vv}\cdot 0 = 0$, since $R^{\mu\nu}$ is nonzero only in the $vv$ slot and $R_{vv} = 0$.
The Ricci scalar was already zero in Step 10 for the same reason in its simplest form, $g^{uu} = 0$.

The general statement is that every scalar polynomial in the curvature vanishes here.
Each contraction in such a polynomial must pair an index of one curvature factor with an index of another; the pairing goes through the inverse metric; and the inverse metric can only pair a $u$ with a $v$, while the curvature has no $v$ slot to offer.
These are the vanishing scalar invariant spacetimes, and the flat metric is not the only one among them.

So the situation is the one the history of the pp-wave describes.
The curvature is not zero: the Weyl tensor of Step 11 is nonzero wherever the transverse Hessian of the profile has a trace free part, and in Step 13 it does measurable work on nearby particles.
Yet every scalar one can build from that curvature is zero, so no scalar can be used to tell a wave from flat space, and no argument conducted in scalars can settle whether a wave is there at all.
The settlement had to come from the relative acceleration of two freely falling particles, which is the Riemann tensor with its indices still on it.

---

## Step 13. The geodesic equations

The geodesic equation is $\ddot{x}^\mu + \Gamma^\mu{}_{\nu\rho}\dot{x}^\nu\dot{x}^\rho = 0$, with the dots denoting derivatives of the chart coordinates with respect to an affine parameter, so $\dot{u}$ means $d(cu)/d\lambda$ and the symbols are the chart symbols of Step 5.
Summing over the pairs, and doubling the symbols that are symmetric in two different indices,

$$\ddot{u} = 0,$$

$$\ddot{v} - \tfrac{1}{2}\partial_u H\,\dot{u}^2 - \partial_x H\,\dot{u}\dot{x} - \partial_y H\,\dot{u}\dot{y} = 0,$$

$$\ddot{x} - \tfrac{1}{2}\partial_x H\,\dot{u}^2 = 0, \qquad \ddot{y} - \tfrac{1}{2}\partial_y H\,\dot{u}^2 = 0,$$

which are the four geodesic equations.
The factors of two in the $v$ equation come from $\Gamma^v{}_{ux}$ and $\Gamma^v{}_{xu}$ being the same symbol counted twice, and the $x$ and $y$ equations have no such factor because $\Gamma^x{}_{uu}$ is a single term.

The first equation is the whole of the $u$ row of Step 5 being empty: $cu$ is an affine parameter along every geodesic, so

$$\dot{u} = \beta, \qquad cu = \beta\lambda + \text{constant}$$

for a constant $\beta$ fixed by the initial data, and a particle's retarded time advances uniformly no matter what the wave does.
For $\beta = 0$ the transverse equations say $\ddot{x} = \ddot{y} = 0$: a particle riding the wave front at the speed of light is not deflected at all.

For $\beta \ne 0$ the transverse equations are

$$\ddot{x}^a = \tfrac{1}{2}\beta^2\,\partial_a H,$$

the profile acting as a potential in the transverse plane with the retarded time as the clock.
Two nearby particles separated by a transverse $\xi^a$ obey the difference of two such equations, which to first order is

$$\ddot{\xi}^a = \tfrac{1}{2}\beta^2\,\partial_a\partial_b H\,\xi^b = -\beta^2\,C^a{}_{ubu}\,\xi^b + \tfrac{1}{4}\beta^2\left(\Delta_\perp H\right)\xi^a,$$

the geodesic deviation driven by the transverse Hessian of Step 7, split into the trace free part that is the Weyl tensor of Step 11 and the trace that is the Ricci tensor of Step 9.
In a vacuum only the first survives, so the driving matrix is trace free: what one transverse direction gains the other loses, and a ring of particles is squeezed into an ellipse of the same area.
The two polarisations of Step 11 are the two ways of doing that, at forty five degrees to each other, and this is the equation an interferometer integrates.

---

## Step 14. The exact plane wave

A profile with no dependence on $x$ and $y$ at all has vanishing Hessian and so, by Step 7, vanishing curvature: it is flat space in a strange chart.
The first profile that is not flat is a quadratic one, and it is the one Bondi, Pirani and Robinson treated.

Write the general quadratic,

$$H = \alpha(u)x^2 + 2\beta(u)xy + \gamma(u)y^2 + \text{terms linear and constant in }x\text{ and }y.$$

The linear and constant terms have vanishing second transverse derivatives, so by Step 7 they contribute nothing to any curvature component and can be dropped without losing a solution.
The vacuum condition of Step 9 is

$$\Delta_\perp H = 2\alpha + 2\gamma = 0,$$

so $\gamma = -\alpha$, and writing $A = \alpha$ and $B = \beta$ the general quadratic vacuum profile is

$$H = A(u)\left(x^2 - y^2\right) + 2B(u)xy = \operatorname{Re}\left[\left(A - iB\right)\left(x + iy\right)^2\right],$$

with $A$ and $B$ arbitrary functions of the retarded time.
Both carry $1/L^2$, since $H$ is dimensionless and $x^2$ is a length squared, and they are the amplitudes of the two polarisations of Step 11.
This is the second chart of the pp-wave, and every component in it comes from the Brinkmann formulas by substituting the four derivatives

$$\partial_x H = 2\left(Ax + By\right), \qquad \partial_y H = 2\left(Bx - Ay\right), \qquad \partial_u H = A'\left(x^2 - y^2\right) + 2B'xy,$$

$$\partial_x^2 H = 2A, \qquad \partial_y^2 H = -2A, \qquad \partial_x\partial_y H = 2B,$$

where a prime is a derivative with respect to the chart coordinate $cu$, by the convention of Step 1.

The Christoffel symbols of Step 5 become

$$\Gamma^v{}_{uu} = -\tfrac{1}{2}\left(A'\left(x^2 - y^2\right) + 2B'xy\right), \qquad \Gamma^v{}_{ux} = \Gamma^x{}_{uu} = -\left(Ax + By\right), \qquad \Gamma^v{}_{uy} = \Gamma^y{}_{uu} = Ay - Bx,$$

with the lowered symbols of Step 4 equal to these in the $x$ and $y$ rows and to the negative of the $v$ row in its $u$ row, which is the raising of Step 5 in reverse.
The Riemann tensor of Step 7 becomes

$$R^v{}_{xux} = R^x{}_{uux} = A, \qquad R^v{}_{yuy} = R^y{}_{uuy} = -A, \qquad R^v{}_{xuy} = R^v{}_{yux} = B,$$

and the rest by the symmetries, with the sign flips in the slots of Step 7.
The Ricci tensor, the Ricci scalar, the Einstein tensor and the Kretschmann scalar are all zero: the first two because the profile is harmonic by construction, the last by the argument of Step 12, which needed no vacuum.
The Weyl tensor equals the Riemann tensor, for the reason given in Step 11.

The curvature is a pair of functions of $u$ alone, so the tidal field is the same everywhere on a wave front, which is what makes this the exact analogue of a plane wave in electromagnetism and what the word plane in its name means beyond the flatness of the fronts.
Constant $A$ and $B$ give a homogeneous wave that never turns off.
Amplitudes vanishing outside an interval of $u$ give a sandwich wave, flat before and flat after, with a burst in between; the geodesics of Step 13 enter parallel and leave converging, which is the focusing Penrose used to show that a plane wave admits no Cauchy surface.
That is the last piece of the pp-wave's history, and it is visible in the oscillator equation $\ddot{x} = \beta^2(Ax + By)$, whose solutions are focused by any pulse of $A$ of one sign no matter how brief.

---

## Step 15. Every equation is dimensionally consistent

The dimensions are $[u] = T$, $[v] = [x] = [y] = L$ and $[H] = 1$, with $[A] = [B] = 1/L^2$ in the second chart.
Since the chart multiplies $u$ by $c$, every chart coordinate carries a length, so a chart component of the metric is dimensionless, a Christoffel symbol carries $1/L$, and a Riemann, Ricci, Einstein or Weyl component carries $1/L^2$.

The line element: $H\,c^2du^2$ carries $1\cdot(L/T)^2T^2 = L^2$, the mixed term $2c\,du\,dv$ carries $(L/T)\cdot T\cdot L = L^2$, and $dx^2$ carries $L^2$.
All three agree, which is what fixes $v$ as a length.

The profile derivatives: $\partial_u H = c^{-1}\partial H/\partial u$ carries $(T/L)(1/T) = 1/L$, and $\partial_x H$ carries $1/L$ directly.
So every Christoffel symbol of Step 5 carries $1/L$, as a Christoffel symbol must.
Had $\partial_u H$ been read as the bare $\partial H/\partial u$ it would carry $1/T$, and $\Gamma^v{}_{uu}$ would have been wrong by one factor of $c$ against its own transverse neighbours $\Gamma^v{}_{ux}$, which is the error the dimensional pass of the checker exists to catch.

The second derivatives $\partial_a\partial_b H$ carry $1/L^2$, so every Riemann, Ricci, Einstein and Weyl component of Steps 7 to 11 carries $1/L^2$.

The geodesics: each term is measured against $\ddot{x}^\mu$, which carries $[x^\mu_{\text{chart}}]/\lambda^2 = L/\lambda^2$, and a chart velocity $\dot{x}^\mu$ carries $L/\lambda$.
In the $x$ equation, $\partial_x H\,\dot{u}^2$ carries $(1/L)(L/\lambda)^2 = L/\lambda^2$.
In the $v$ equation, $\partial_u H\,\dot{u}^2$ and $\partial_x H\,\dot{u}\dot{x}$ both carry $L/\lambda^2$, and they do so only because the written $\partial_u$ carries the factor of $c$ that the dot on $\dot u$ also carries.

The second chart: $Ax$ carries $(1/L^2)L = 1/L$, matching a Christoffel symbol; $A$ carries $1/L^2$, matching a Riemann component; and $A' = c^{-1}dA/du$ carries $(T/L)(1/L^2)(1/T) = 1/L^3$, so $A'(x^2 - y^2)$ carries $1/L$ and sits correctly in $\Gamma^v{}_{uu}$.

---

## Step 16. The components, and what the checker needed

In each of the two charts we compute the line element, the metric and its inverse, the Christoffel symbols with the first index up and with every index lowered, the Riemann tensor with the first index up and with every index lowered, the Ricci tensor with both indices down, with mixed indices and with both up, the Ricci scalar, the Kretschmann scalar, the Einstein tensor in the same three index positions, the Weyl tensor with the first index up and with every index lowered, and the four geodesic equations.
In the Brinkmann chart they hold for an arbitrary profile, so its Ricci and Einstein tensors do not vanish, and the vacuum condition is stated in the conventions of the pp-wave and in the description of $H$ rather than imposed on the components.
The exact plane wave chart has the condition already solved, so its Ricci and Einstein tensors vanish and its components hold for any amplitudes at all.

`verify_metrics.py` needed one declaration per chart in `DIMENSIONS`, each declaring $u$ a time and $v$, $x$ and $y$ lengths, with $[H] = 1$ in the first and $[A] = [B] = 1/L^2$ in the second.
No relation in `PARAMETER_RELATIONS` is needed for either.
Nothing in the Brinkmann chart requires the profile to be harmonic, so its parameters are free; the components of the plane wave chart hold for every $A$ and $B$, because the harmonic condition is already built into the shape of its profile rather than left as a constraint between parameters.
That is the opposite of Kasner, whose exponents are bound by two equations and whose components are true only on the surface those equations cut out.

Two things in the checker's reader are new with the pp-wave, and both are general rather than particular to it.
A parameter may now be declared as a function of several coordinates, as `H = H(u,x,y)`, where before only a function of one coordinate was understood.
And a partial derivative may be written with $\partial$, as $\partial_x H$ or $\partial_x^2 H$ or $\partial_x\partial_y H$, which the reader expands into the name it has declared for that derivative.
Every first and second partial derivative of a declared function was declared, which was as far as any curvature tensor of any spacetime then reached, and a function of one coordinate still answers to a dot and to a prime as it did before.
Tolman-Bondi later reached a third, and since 22 September 2026 the reader declares a partial derivative of any order when a component names it; the reason is in Step 18 of `tolman_bondi.md`.
Like the dot and the prime, the $\partial$ is taken with respect to the chart coordinate, so $\partial_u$ carries the $1/c$ of Step 1 and $\partial_x$ does not, and that convention is what makes every term of Step 15 balance.
