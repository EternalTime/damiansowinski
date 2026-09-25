# Null coordinates on Minkowski spacetime

Minkowski spacetime carries two null coordinate systems, Double Null and Spherical Null.
We derive every component of each, in order, from the transformation down to the geodesic equations.
We leave nothing as an exercise and assert nothing we do not compute.

The companion script `_tools/derivations/verify_metrics.py` does the same work in sympy and compares every component it computes with the stated one, so the algebra is checkable by hand and by machine independently.
That script checks every coordinate system of every spacetime the same way, not only these two.

---

## Step 1. Conventions

The signature is $(-,+,+,+)$ and the Riemann tensor is

$$R^\mu{}_{\nu\rho\sigma} = \partial_\rho \Gamma^\mu_{\nu\sigma} - \partial_\sigma \Gamma^\mu_{\nu\rho} + \Gamma^\mu_{\rho\lambda}\Gamma^\lambda_{\nu\sigma} - \Gamma^\mu_{\sigma\lambda}\Gamma^\lambda_{\nu\rho},$$

with the Christoffel symbols of the Levi-Civita connection

$$\Gamma^\mu_{\nu\rho} = \tfrac{1}{2}g^{\mu\alpha}\left(\partial_\nu g_{\alpha\rho} + \partial_\rho g_{\alpha\nu} - \partial_\alpha g_{\nu\rho}\right),$$

and the lowered symbol defined by $\Gamma_{\mu\nu\rho} = g_{\mu\alpha}\Gamma^\alpha_{\nu\rho}$.
The Ricci tensor is the contraction $R_{\nu\sigma} = R^\mu{}_{\nu\mu\sigma}$.

Factors of $c$ are kept explicit.
The coordinates $u$ and $v$ carry dimensions of time, exactly as $t$ does, since they are built as $t \pm (\text{length})/c$.

The chart, here and for every other spacetime, is the one whose time coordinate is $x^0 = cT$.
The Cartesian and spherical charts of Minkowski spacetime have $g_{tt} = -1$ against a line element $-c^2dt^2$, which is that chart; so does Schwarzschild, with $g_{tt} = -(1-r_s/r)$ against a time term $-(1-r_s/r)c^2dt^2$.
The two null systems follow it too, and that is where the factors of $c$ in their components come from.
Since $u$ and $v$ are times, the chart coordinates are the rescaled $cu$ and $cv$, and every component with an index $u$ or $v$ is a component in that rescaled chart even though the index is written with the bare letter.

Because the rescaling is linear with constant coefficients, the rule is arithmetic: a component in the chart is the component taken with the bare coordinate, multiplied by $c$ once for every upper $u$ or $v$ index and divided by $c$ once for every lower one.
The Christoffel symbols obey the same rule as the tensors, since the inhomogeneous term in their transformation law carries a second derivative of the coordinate change, which vanishes for a linear one.
This is why $\Gamma^u{}_{\theta\theta}$ carries a $c$ and $\Gamma^\theta{}_{u\theta}$ carries a $1/c$, while $\Gamma^\theta{}_{\phi\phi}$ and $\Gamma^\phi{}_{\theta\phi}$, which have no null index at all, carry none.
The line element itself is written in $du$ and $dv$ throughout, so it is unchanged by any of this; only the components sit in the $cu, cv$ chart.

In Steps 4 to 7 and Steps 14 to 18 we work in the bare $(u,v)$ chart, because that is where the transformation from Cartesian naturally lands.
In Steps 7a, 17a and 22a we then carry the results into the rescaled chart.

The starting point in both parts is the Cartesian line element of Minkowski spacetime,

$$ds^2 = -c^2dt^2 + dx^2 + dy^2 + dz^2,$$

so the Cartesian metric in the chart $(t,x,y,z)$ is $\eta_{ab} = \mathrm{diag}(-c^2, 1, 1, 1)$.

---

# Part I. Double null coordinates

## Step 2. The transformation and its inverse

Define

$$u = t - \frac{x}{c}, \qquad v = t + \frac{x}{c},$$

keeping $y$ and $z$ untouched.
Adding and subtracting the two definitions inverts them,

$$t = \frac{u+v}{2}, \qquad x = \frac{c(v-u)}{2}.$$

Both $u$ and $v$ run over the whole real line, and the map is a linear bijection of $\mathbb{R}^2$ onto $\mathbb{R}^2$, so the chart covers all of Minkowski spacetime with no boundary and no excluded point.
A surface of constant $u$ is the worldsheet of a light front moving in the $+x$ direction, and a surface of constant $v$ is one moving in the $-x$ direction.

## Step 3. The Jacobian

Write $x^a = (t,x,y,z)$ for the Cartesian chart and $\tilde x^\mu = (u,v,y,z)$ for the new one.
Differentiating the inverse transformation of Step 2 entry by entry gives

$$J^a{}_\mu = \frac{\partial x^a}{\partial \tilde x^\mu} = \begin{pmatrix} \tfrac{1}{2} & \tfrac{1}{2} & 0 & 0 \\[2pt] -\tfrac{c}{2} & \tfrac{c}{2} & 0 & 0 \\[2pt] 0 & 0 & 1 & 0 \\[2pt] 0 & 0 & 0 & 1 \end{pmatrix},$$

rows labelled by $a = t,x,y,z$ and columns by $\mu = u,v,y,z$.
Its determinant is $\left(\tfrac{1}{2}\cdot\tfrac{c}{2} - \tfrac{1}{2}\cdot\left(-\tfrac{c}{2}\right)\right)\cdot 1 \cdot 1 = \tfrac{c}{2}$, which is nonzero, so the chart is regular everywhere.

## Step 4. The metric by substituting the differentials

From Step 2, $dt = \tfrac{1}{2}(du + dv)$ and $dx = \tfrac{c}{2}(dv - du)$.
Substituting both into the Cartesian line element,

$$-c^2dt^2 = -\frac{c^2}{4}\left(du^2 + 2\,du\,dv + dv^2\right),$$

$$dx^2 = \frac{c^2}{4}\left(du^2 - 2\,du\,dv + dv^2\right).$$

The $du^2$ terms cancel, the $dv^2$ terms cancel, and the two cross terms add,

$$-c^2dt^2 + dx^2 = -\frac{c^2}{4}\left(2\,du\,dv\right) - \frac{c^2}{4}\left(2\,du\,dv\right) = -c^2\,du\,dv.$$

Therefore

$$ds^2 = -c^2\,du\,dv + dy^2 + dz^2,$$

which is the double null line element.

## Step 5. The metric by contracting the Jacobian

The same result follows from $g_{\mu\nu} = J^a{}_\mu J^b{}_\nu \eta_{ab}$, which for $\eta_{ab} = \mathrm{diag}(-c^2,1,1,1)$ reads

$$g_{\mu\nu} = -c^2 \frac{\partial t}{\partial \tilde x^\mu}\frac{\partial t}{\partial \tilde x^\nu} + \frac{\partial x}{\partial \tilde x^\mu}\frac{\partial x}{\partial \tilde x^\nu} + \frac{\partial y}{\partial \tilde x^\mu}\frac{\partial y}{\partial \tilde x^\nu} + \frac{\partial z}{\partial \tilde x^\mu}\frac{\partial z}{\partial \tilde x^\nu}.$$

Taking the ten independent components in turn, using the columns of Step 3:

$$g_{uu} = -c^2\left(\tfrac{1}{2}\right)^2 + \left(-\tfrac{c}{2}\right)^2 = -\frac{c^2}{4} + \frac{c^2}{4} = 0,$$

$$g_{vv} = -c^2\left(\tfrac{1}{2}\right)^2 + \left(\tfrac{c}{2}\right)^2 = -\frac{c^2}{4} + \frac{c^2}{4} = 0,$$

$$g_{uv} = -c^2\left(\tfrac{1}{2}\right)\left(\tfrac{1}{2}\right) + \left(-\tfrac{c}{2}\right)\left(\tfrac{c}{2}\right) = -\frac{c^2}{4} - \frac{c^2}{4} = -\frac{c^2}{2},$$

$$g_{uy} = g_{uz} = g_{vy} = g_{vz} = 0,$$

since the $u$ and $v$ columns have no $y$ or $z$ entry and the $y$ and $z$ columns have no $t$ or $x$ entry, and finally

$$g_{yy} = g_{zz} = 1, \qquad g_{yz} = 0.$$

## Step 6. The metric matrix, its determinant and the null directions

In the order $(u,v,y,z)$,

$$g_{\mu\nu} = \begin{pmatrix} 0 & -\tfrac{c^2}{2} & 0 & 0 \\[2pt] -\tfrac{c^2}{2} & 0 & 0 & 0 \\[2pt] 0 & 0 & 1 & 0 \\[2pt] 0 & 0 & 0 & 1\end{pmatrix}, \qquad \det g = \left(0 \cdot 0 - \tfrac{c^4}{4}\right)\cdot 1 \cdot 1 = -\frac{c^4}{4}.$$

The determinant is negative, so the signature is still Lorentzian, and it is nowhere zero, so the chart is nowhere degenerate.

The two diagonal zeros are the whole point of the system.
The coordinate vector fields are

$$\partial_u = \frac{\partial t}{\partial u}\partial_t + \frac{\partial x}{\partial u}\partial_x = \frac{1}{2}\partial_t - \frac{c}{2}\partial_x, \qquad \partial_v = \frac{1}{2}\partial_t + \frac{c}{2}\partial_x,$$

and their norms in the Cartesian chart are

$$g(\partial_u,\partial_u) = -c^2\left(\tfrac{1}{2}\right)^2 + \left(-\tfrac{c}{2}\right)^2 = 0, \qquad g(\partial_v,\partial_v) = -c^2\left(\tfrac{1}{2}\right)^2 + \left(\tfrac{c}{2}\right)^2 = 0.$$

Both coordinate directions are null.
They are not parallel, because $g(\partial_u,\partial_v) = -c^2/2 \neq 0$, so together they span a timelike plane.
The inverse metric of Step 7 has $g^{uu} = g^{vv} = 0$, which says that $\nabla u$ and $\nabla v$ are null covectors, so the level sets $u = \text{const}$ and $v = \text{const}$ are null hypersurfaces.

## Step 7. The inverse metric

The $2 \times 2$ block $\begin{pmatrix} 0 & k \\ k & 0\end{pmatrix}$ has inverse $\begin{pmatrix} 0 & 1/k \\ 1/k & 0\end{pmatrix}$, since the product of the two is $\begin{pmatrix} k/k & 0 \\ 0 & k/k \end{pmatrix}$.
With $k = -c^2/2$ this gives $1/k = -2/c^2$, and the $y,z$ block is its own inverse, so

$$g^{\mu\nu} = \begin{pmatrix} 0 & -\tfrac{2}{c^2} & 0 & 0 \\[2pt] -\tfrac{2}{c^2} & 0 & 0 & 0 \\[2pt] 0 & 0 & 1 & 0 \\[2pt] 0 & 0 & 0 & 1\end{pmatrix}.$$

Checking $g^{\mu\alpha}g_{\alpha\nu} = \delta^\mu{}_\nu$ component by component, with $\alpha$ summed over all four values:

$$g^{u\alpha}g_{\alpha u} = g^{uv}g_{vu} = \left(-\frac{2}{c^2}\right)\left(-\frac{c^2}{2}\right) = 1,$$

$$g^{u\alpha}g_{\alpha v} = g^{uv}g_{vv} = \left(-\frac{2}{c^2}\right)(0) = 0,$$

$$g^{v\alpha}g_{\alpha v} = g^{vu}g_{uv} = \left(-\frac{2}{c^2}\right)\left(-\frac{c^2}{2}\right) = 1,$$

$$g^{v\alpha}g_{\alpha u} = g^{vu}g_{uu} = \left(-\frac{2}{c^2}\right)(0) = 0,$$

$$g^{y\alpha}g_{\alpha y} = g^{yy}g_{yy} = 1, \qquad g^{z\alpha}g_{\alpha z} = g^{zz}g_{zz} = 1,$$

and every remaining mixed component, such as $g^{u\alpha}g_{\alpha y}$, vanishes because the only nonzero $g^{u\alpha}$ is $g^{uv}$ and $g_{vy} = 0$.
The product is the identity.

## Step 7a. The same metric in the rescaled chart

In Steps 4 to 7 we worked with the bare $u$ and $v$, which are times.
The chart of Step 1 has the coordinates $cu$ and $cv$, so each lower null index divides by $c$ and each upper one multiplies by it.
Only $g_{uv}$ and $g^{uv}$ carry null indices at all, and each carries two:

$$g_{uv} = g_{vu} = \frac{1}{c^2}\left(-\frac{c^2}{2}\right) = -\frac{1}{2}, \qquad g^{uv} = g^{vu} = c^2\left(-\frac{2}{c^2}\right) = -2.$$

The components $g_{yy} = g_{zz} = 1$ have no null index and are unchanged, and the product is still the identity, since each factor of $c$ introduced on the inverse is cancelled by the one removed from the metric.
These four numbers are the components of the metric and its inverse in the rescaled chart.
In the rescaled chart the flat metric in double null coordinates carries no $c$ at all, which is the point of using that chart.

## Step 8. Every Christoffel symbol

Every entry of $g_{\mu\nu}$ in Step 6 is a constant: $0$, $-c^2/2$ and $1$ are all independent of $u$, $v$, $y$ and $z$.
Therefore

$$\partial_\nu g_{\alpha\rho} = 0 \quad \text{for every choice of } \nu, \alpha, \rho.$$

The Christoffel formula is built from exactly three such derivatives,

$$\Gamma^\mu_{\nu\rho} = \tfrac{1}{2}g^{\mu\alpha}\left(\partial_\nu g_{\alpha\rho} + \partial_\rho g_{\alpha\nu} - \partial_\alpha g_{\nu\rho}\right) = \tfrac{1}{2}g^{\mu\alpha}\left(0 + 0 - 0\right) = 0,$$

for all $4^3 = 64$ index combinations.
Every Christoffel symbol vanishes, both $\Gamma^\mu_{\nu\rho}$ and $\Gamma_{\mu\nu\rho}$, the second because $\Gamma_{\mu\nu\rho} = g_{\mu\alpha}\Gamma^\alpha_{\nu\rho}$ is a sum of multiples of symbols that are already zero.
The checker therefore expects no nonzero Christoffel symbol of either kind.

Note that this is a stronger statement than it looks.
The metric is not diagonal, and the coordinates are not orthogonal, yet the connection still vanishes, because what the connection sees is not whether the metric is diagonal but whether it is constant.

## Step 9. Riemann, Ricci, Einstein and Weyl

With every $\Gamma^\mu_{\nu\rho} = 0$ from Step 8, each of the four terms of the Riemann tensor is zero separately:

$$R^\mu{}_{\nu\rho\sigma} = \partial_\rho(0) - \partial_\sigma(0) + (0)(0) - (0)(0) = 0.$$

The contraction $R_{\nu\sigma} = R^\mu{}_{\nu\mu\sigma}$ is then a sum of zeros, so the Ricci tensor vanishes, and $R = g^{\nu\sigma}R_{\nu\sigma} = 0$.
The Kretschmann scalar $K = R_{\mu\nu\rho\sigma}R^{\mu\nu\rho\sigma}$ is a sum of products of vanishing components, so $K = 0$.
The Einstein tensor $G_{\mu\nu} = R_{\mu\nu} - \tfrac{1}{2}Rg_{\mu\nu} = 0 - 0 = 0$.
The Weyl tensor is Riemann minus terms built entirely from Ricci and $R$, all of which are zero, so it vanishes too.
The Riemann, Ricci, Einstein and Weyl tensors all vanish, and $R = 0$, $K = 0$ for the two scalars.

## Step 10. The geodesic equations

The geodesic equation is $\ddot x^\mu + \Gamma^\mu_{\nu\rho}\dot x^\nu \dot x^\rho = 0$, where the dot is $d/d\lambda$ for an affine parameter $\lambda$.
The vanishing symbols of Step 8 kill the second term for every $\mu$, leaving

$$\ddot u = 0, \qquad \ddot v = 0, \qquad \ddot y = 0, \qquad \ddot z = 0,$$

which are the geodesic equations in double null coordinates.
The solutions are affine in $\lambda$ in each of the four coordinates, as they must be, since they are linear combinations of the Cartesian straight lines of Step 2.

---

# Part II. Spherical null coordinates

## Step 11. The transformation and its inverse

Let $r = \sqrt{x^2+y^2+z^2}$ be the Cartesian radius and define the retarded and advanced times

$$u = t - \frac{r}{c}, \qquad v = t + \frac{r}{c},$$

keeping the angles $\theta$ and $\phi$.
Adding and subtracting inverts them,

$$t = \frac{u+v}{2}, \qquad r = \frac{c(v-u)}{2}.$$

Since $r \geq 0$, the chart requires $v \geq u$, and the surface $v = u$ is the worldline $r=0$ where the angles degenerate, exactly as they do at the origin of the ordinary spherical chart.
An outgoing light ray travelling radially has $dr = c\,dt$ and therefore $du = 0$, so $u$ labels outgoing rays; an ingoing ray has $dr = -c\,dt$ and therefore $dv = 0$, so $v$ labels ingoing rays.
This is the pair of labels a Penrose diagram is drawn in.

## Step 12. The Jacobian

Write the unit radial vector

$$\hat n(\theta,\phi) = \left(\sin\theta\cos\phi,\; \sin\theta\sin\phi,\; \cos\theta\right),$$

so that the Cartesian position is $(x,y,z) = r\,\hat n$ with $r = c(v-u)/2$.
Differentiating,

$$\frac{\partial t}{\partial u} = \frac{1}{2}, \quad \frac{\partial t}{\partial v} = \frac{1}{2}, \quad \frac{\partial t}{\partial \theta} = \frac{\partial t}{\partial \phi} = 0,$$

$$\frac{\partial (x,y,z)}{\partial u} = -\frac{c}{2}\hat n, \quad \frac{\partial (x,y,z)}{\partial v} = \frac{c}{2}\hat n, \quad \frac{\partial (x,y,z)}{\partial \theta} = r\,\partial_\theta\hat n, \quad \frac{\partial (x,y,z)}{\partial \phi} = r\,\partial_\phi\hat n,$$

with

$$\partial_\theta \hat n = \left(\cos\theta\cos\phi,\; \cos\theta\sin\phi,\; -\sin\theta\right), \qquad \partial_\phi \hat n = \left(-\sin\theta\sin\phi,\; \sin\theta\cos\phi,\; 0\right).$$

Written as a matrix with rows $(t,x,y,z)$ and columns $(u,v,\theta,\phi)$,

$$J^a{}_\mu = \begin{pmatrix} \tfrac{1}{2} & \tfrac{1}{2} & 0 & 0 \\[2pt] -\tfrac{c}{2}\hat n^1 & \tfrac{c}{2}\hat n^1 & r\,\partial_\theta \hat n^1 & r\,\partial_\phi \hat n^1 \\[2pt] -\tfrac{c}{2}\hat n^2 & \tfrac{c}{2}\hat n^2 & r\,\partial_\theta \hat n^2 & r\,\partial_\phi \hat n^2 \\[2pt] -\tfrac{c}{2}\hat n^3 & \tfrac{c}{2}\hat n^3 & r\,\partial_\theta \hat n^3 & r\,\partial_\phi \hat n^3 \end{pmatrix}.$$

## Step 13. The six dot products the transformation needs

Every component of the transformed metric is one of six dot products of the three vectors of Step 12, so we compute all six explicitly.

$$\hat n \cdot \hat n = \sin^2\theta\cos^2\phi + \sin^2\theta\sin^2\phi + \cos^2\theta = \sin^2\theta + \cos^2\theta = 1,$$

$$\hat n \cdot \partial_\theta \hat n = \sin\theta\cos\theta\cos^2\phi + \sin\theta\cos\theta\sin^2\phi - \cos\theta\sin\theta = \sin\theta\cos\theta - \sin\theta\cos\theta = 0,$$

$$\hat n \cdot \partial_\phi \hat n = -\sin^2\theta\cos\phi\sin\phi + \sin^2\theta\sin\phi\cos\phi + 0 = 0,$$

$$\partial_\theta \hat n \cdot \partial_\theta \hat n = \cos^2\theta\cos^2\phi + \cos^2\theta\sin^2\phi + \sin^2\theta = \cos^2\theta + \sin^2\theta = 1,$$

$$\partial_\theta \hat n \cdot \partial_\phi \hat n = -\sin\theta\cos\theta\cos\phi\sin\phi + \sin\theta\cos\theta\sin\phi\cos\phi + 0 = 0,$$

$$\partial_\phi \hat n \cdot \partial_\phi \hat n = \sin^2\theta\sin^2\phi + \sin^2\theta\cos^2\phi = \sin^2\theta.$$

## Step 14. The metric

Contract the Jacobian of Step 12 with $\eta_{ab} = \mathrm{diag}(-c^2,1,1,1)$, so that

$$g_{\mu\nu} = -c^2\frac{\partial t}{\partial \tilde x^\mu}\frac{\partial t}{\partial \tilde x^\nu} + \frac{\partial (x,y,z)}{\partial \tilde x^\mu}\cdot\frac{\partial (x,y,z)}{\partial \tilde x^\nu},$$

and take the ten independent components using Step 13.

$$g_{uu} = -c^2\left(\tfrac{1}{2}\right)^2 + \left(-\tfrac{c}{2}\right)^2 \hat n\cdot\hat n = -\frac{c^2}{4} + \frac{c^2}{4} = 0,$$

$$g_{vv} = -c^2\left(\tfrac{1}{2}\right)^2 + \left(\tfrac{c}{2}\right)^2 \hat n\cdot\hat n = -\frac{c^2}{4} + \frac{c^2}{4} = 0,$$

$$g_{uv} = -c^2\left(\tfrac{1}{2}\right)\left(\tfrac{1}{2}\right) + \left(-\tfrac{c}{2}\right)\left(\tfrac{c}{2}\right)\hat n\cdot\hat n = -\frac{c^2}{4} - \frac{c^2}{4} = -\frac{c^2}{2},$$

$$g_{u\theta} = -c^2\left(\tfrac{1}{2}\right)(0) + \left(-\tfrac{c}{2}\right)r\;\hat n\cdot\partial_\theta\hat n = 0,$$

$$g_{u\phi} = \left(-\tfrac{c}{2}\right)r\;\hat n\cdot\partial_\phi\hat n = 0, \qquad g_{v\theta} = \left(\tfrac{c}{2}\right)r\;\hat n\cdot\partial_\theta\hat n = 0, \qquad g_{v\phi} = \left(\tfrac{c}{2}\right)r\;\hat n\cdot\partial_\phi\hat n = 0,$$

$$g_{\theta\theta} = r^2\;\partial_\theta\hat n\cdot\partial_\theta\hat n = r^2, \qquad g_{\theta\phi} = r^2\;\partial_\theta\hat n\cdot\partial_\phi\hat n = 0, \qquad g_{\phi\phi} = r^2\;\partial_\phi\hat n\cdot\partial_\phi\hat n = r^2\sin^2\theta.$$

Substituting $r = c(v-u)/2$, so that $r^2 = c^2(v-u)^2/4$, the line element is

$$ds^2 = -c^2\,du\,dv + \frac{c^2(v-u)^2}{4}\left(d\theta^2 + \sin^2\theta\,d\phi^2\right),$$

which is the spherical null line element.

## Step 15. The same result starting from the spherical chart

As a second route, start instead from the spherical line element of Minkowski spacetime,

$$ds^2 = -c^2dt^2 + dr^2 + r^2d\theta^2 + r^2\sin^2\theta\,d\phi^2,$$

and use only $t = (u+v)/2$ and $r = c(v-u)/2$ with the angles fixed.
Then $dt = \tfrac{1}{2}(du+dv)$ and $dr = \tfrac{c}{2}(dv-du)$, which are the same two differentials as in Step 4, so

$$-c^2dt^2 + dr^2 = -\frac{c^2}{4}\left(du+dv\right)^2 + \frac{c^2}{4}\left(dv-du\right)^2 = -c^2\,du\,dv,$$

and the angular part is carried over unchanged with $r^2 = c^2(v-u)^2/4$.
The two routes agree.

## Step 16. The metric matrix, its determinant and the null directions

In the order $(u,v,\theta,\phi)$, writing $r = c(v-u)/2$ throughout,

$$g_{\mu\nu} = \begin{pmatrix} 0 & -\tfrac{c^2}{2} & 0 & 0 \\[2pt] -\tfrac{c^2}{2} & 0 & 0 & 0 \\[2pt] 0 & 0 & \tfrac{c^2(v-u)^2}{4} & 0 \\[2pt] 0 & 0 & 0 & \tfrac{c^2(v-u)^2}{4}\sin^2\theta\end{pmatrix}.$$

The matrix is block diagonal, so its determinant is the product of the two blocks,

$$\det g = \left(-\frac{c^4}{4}\right)\cdot \frac{c^2(v-u)^2}{4}\cdot\frac{c^2(v-u)^2\sin^2\theta}{4} = -\frac{c^8(v-u)^4\sin^2\theta}{64}.$$

This is negative wherever $v > u$ and $\sin\theta \neq 0$, so the signature is Lorentzian on the whole chart except at the axis and at $r=0$, which are the same two degeneracies the ordinary spherical chart has.

As in Step 6, $g_{uu} = g_{vv} = 0$ makes $\partial_u$ and $\partial_v$ null, and $g^{uu} = g^{vv} = 0$ in Step 17 makes the light cones $u = \text{const}$ and $v = \text{const}$ null hypersurfaces.
Here they are the outgoing and ingoing light cones of the events on the worldline $r=0$.

## Step 17. The inverse metric

The block structure of Step 16 gives the inverse block by block.
The $(u,v)$ block is inverted exactly as in Step 7, giving $g^{uv} = g^{vu} = -2/c^2$ and $g^{uu} = g^{vv} = 0$.
The two diagonal entries invert as reciprocals,

$$g^{\theta\theta} = \frac{4}{c^2(v-u)^2}, \qquad g^{\phi\phi} = \frac{4}{c^2(v-u)^2\sin^2\theta}.$$

Checking $g^{\mu\alpha}g_{\alpha\nu} = \delta^\mu{}_\nu$ component by component:

$$g^{u\alpha}g_{\alpha u} = g^{uv}g_{vu} = \left(-\frac{2}{c^2}\right)\left(-\frac{c^2}{2}\right) = 1, \qquad g^{u\alpha}g_{\alpha v} = g^{uv}g_{vv} = 0,$$

$$g^{v\alpha}g_{\alpha v} = g^{vu}g_{uv} = \left(-\frac{2}{c^2}\right)\left(-\frac{c^2}{2}\right) = 1, \qquad g^{v\alpha}g_{\alpha u} = g^{vu}g_{uu} = 0,$$

$$g^{\theta\alpha}g_{\alpha\theta} = g^{\theta\theta}g_{\theta\theta} = \frac{4}{c^2(v-u)^2}\cdot\frac{c^2(v-u)^2}{4} = 1,$$

$$g^{\phi\alpha}g_{\alpha\phi} = g^{\phi\phi}g_{\phi\phi} = \frac{4}{c^2(v-u)^2\sin^2\theta}\cdot\frac{c^2(v-u)^2\sin^2\theta}{4} = 1.$$

Every off block component vanishes, because the only nonzero $g^{u\alpha}$ is $g^{uv}$ and $g_{v\theta} = g_{v\phi} = 0$, the only nonzero $g^{\theta\alpha}$ is $g^{\theta\theta}$ and $g_{\theta u} = g_{\theta v} = g_{\theta\phi} = 0$, and likewise for $v$ and $\phi$.
The product is the identity.

## Step 17a. The same metric in the rescaled chart

As in Step 7a, we pass to the chart whose null coordinates are $cu$ and $cv$.
The two components carrying null indices carry two apiece,

$$g_{uv} = g_{vu} = -\frac{1}{2}, \qquad g^{uv} = g^{vu} = -2,$$

while $g_{\theta\theta}$, $g_{\phi\phi}$ and their inverses carry none, so they keep exactly the values of Steps 16 and 17,

$$g_{\theta\theta} = \frac{c^2(v-u)^2}{4}, \qquad g_{\phi\phi} = \frac{c^2(v-u)^2\sin^2\theta}{4}.$$

The factor $c^2(v-u)^2/4$ is just $r^2$ written out, and $r$ is a length however the null coordinates are scaled, which is why no $c$ moves here.

## Step 18. The derivatives of the metric

There are only five nonvanishing first derivatives of the metric in this chart, and the whole of the connection is built from them.

The component $g_{uv} = g_{vu} = -c^2/2$ is a constant, so $\partial_\lambda g_{uv} = 0$ for every $\lambda$.
The components $g_{uu}$, $g_{vv}$, $g_{u\theta}$, $g_{u\phi}$, $g_{v\theta}$, $g_{v\phi}$ and $g_{\theta\phi}$ are identically zero, so their derivatives vanish too.
That leaves $g_{\theta\theta}$ and $g_{\phi\phi}$:

$$\partial_u g_{\theta\theta} = \partial_u\left(\frac{c^2(v-u)^2}{4}\right) = -\frac{c^2(v-u)}{2}, \qquad \partial_v g_{\theta\theta} = \frac{c^2(v-u)}{2}, \qquad \partial_\theta g_{\theta\theta} = \partial_\phi g_{\theta\theta} = 0,$$

$$\partial_u g_{\phi\phi} = -\frac{c^2(v-u)\sin^2\theta}{2}, \qquad \partial_v g_{\phi\phi} = \frac{c^2(v-u)\sin^2\theta}{2},$$

$$\partial_\theta g_{\phi\phi} = \frac{c^2(v-u)^2}{4}\cdot 2\sin\theta\cos\theta = \frac{c^2(v-u)^2\sin\theta\cos\theta}{2}, \qquad \partial_\phi g_{\phi\phi} = 0.$$

Nothing in this chart depends on $\phi$, so $\partial_\phi$ of every metric component is zero.

## Step 19. The symbols with an upper $u$ or $v$

Because $g^{u\alpha}$ is nonzero only for $\alpha = v$, the formula collapses to

$$\Gamma^u_{\nu\rho} = \tfrac{1}{2}g^{uv}\left(\partial_\nu g_{v\rho} + \partial_\rho g_{v\nu} - \partial_v g_{\nu\rho}\right).$$

The first two terms carry a derivative of $g_{v\rho}$ or $g_{v\nu}$.
By Step 18 every component of the metric with a $v$ index is either identically zero or the constant $-c^2/2$, so both of those terms vanish for every $\nu$ and $\rho$, and only the third survives:

$$\Gamma^u_{\nu\rho} = -\tfrac{1}{2}g^{uv}\,\partial_v g_{\nu\rho} = \frac{1}{c^2}\,\partial_v g_{\nu\rho}.$$

By Step 18 the derivative $\partial_v g_{\nu\rho}$ is nonzero only for $(\nu\rho) = (\theta\theta)$ and $(\nu\rho) = (\phi\phi)$, which leaves exactly two symbols:

$$\Gamma^u_{\theta\theta} = \frac{1}{c^2}\cdot\frac{c^2(v-u)}{2} = \frac{v-u}{2}, \qquad \Gamma^u_{\phi\phi} = \frac{1}{c^2}\cdot\frac{c^2(v-u)\sin^2\theta}{2} = \frac{(v-u)\sin^2\theta}{2}.$$

The same argument with $u$ and $v$ exchanged gives

$$\Gamma^v_{\nu\rho} = -\tfrac{1}{2}g^{vu}\,\partial_u g_{\nu\rho} = \frac{1}{c^2}\,\partial_u g_{\nu\rho},$$

$$\Gamma^v_{\theta\theta} = \frac{1}{c^2}\left(-\frac{c^2(v-u)}{2}\right) = -\frac{v-u}{2}, \qquad \Gamma^v_{\phi\phi} = -\frac{(v-u)\sin^2\theta}{2}.$$

All other $\Gamma^u_{\nu\rho}$ and $\Gamma^v_{\nu\rho}$ are zero.
Note that the factors of $c$ cancel exactly here: the $c^2$ in the inverse metric is cancelled by the $c^2$ in $\partial g_{\theta\theta}$, which is why these four symbols carry no $c$ at all.

## Step 20. The symbols with an upper $\theta$

Because $g^{\theta\alpha}$ is nonzero only for $\alpha = \theta$,

$$\Gamma^\theta_{\nu\rho} = \tfrac{1}{2}g^{\theta\theta}\left(\partial_\nu g_{\theta\rho} + \partial_\rho g_{\theta\nu} - \partial_\theta g_{\nu\rho}\right).$$

The only nonzero component of the form $g_{\theta\rho}$ is $g_{\theta\theta}$, so the first term needs $\rho = \theta$ and the second needs $\nu = \theta$.
The third term needs $\partial_\theta g_{\nu\rho} \neq 0$, which by Step 18 happens only for $(\nu\rho) = (\phi\phi)$.
Three families survive.

First, $\nu \in \{u,v\}$ and $\rho = \theta$, where only the first term lives:

$$\Gamma^\theta_{u\theta} = \tfrac{1}{2}g^{\theta\theta}\,\partial_u g_{\theta\theta} = \frac{1}{2}\cdot\frac{4}{c^2(v-u)^2}\cdot\left(-\frac{c^2(v-u)}{2}\right) = -\frac{1}{v-u},$$

$$\Gamma^\theta_{v\theta} = \frac{1}{2}\cdot\frac{4}{c^2(v-u)^2}\cdot\frac{c^2(v-u)}{2} = \frac{1}{v-u}.$$

Second, $\nu = \theta$ and $\rho \in \{u,v\}$, where only the second term lives and gives the same numbers, as the symmetry $\Gamma^\mu_{\nu\rho} = \Gamma^\mu_{\rho\nu}$ requires:

$$\Gamma^\theta_{\theta u} = -\frac{1}{v-u}, \qquad \Gamma^\theta_{\theta v} = \frac{1}{v-u}.$$

Third, $(\nu\rho) = (\phi\phi)$, where only the third term lives:

$$\Gamma^\theta_{\phi\phi} = -\tfrac{1}{2}g^{\theta\theta}\,\partial_\theta g_{\phi\phi} = -\frac{1}{2}\cdot\frac{4}{c^2(v-u)^2}\cdot\frac{c^2(v-u)^2\sin\theta\cos\theta}{2} = -\sin\theta\cos\theta.$$

The case $\nu = \rho = \theta$ would take the first two terms and the third together, giving $\tfrac{1}{2}g^{\theta\theta}\left(\partial_\theta g_{\theta\theta} + \partial_\theta g_{\theta\theta} - \partial_\theta g_{\theta\theta}\right) = \tfrac{1}{2}g^{\theta\theta}\partial_\theta g_{\theta\theta} = 0$, since $g_{\theta\theta}$ does not depend on $\theta$.
Every other $\Gamma^\theta_{\nu\rho}$ is zero.

## Step 21. The symbols with an upper $\phi$

Because $g^{\phi\alpha}$ is nonzero only for $\alpha = \phi$,

$$\Gamma^\phi_{\nu\rho} = \tfrac{1}{2}g^{\phi\phi}\left(\partial_\nu g_{\phi\rho} + \partial_\rho g_{\phi\nu} - \partial_\phi g_{\nu\rho}\right).$$

The third term vanishes for every $\nu$ and $\rho$, because nothing depends on $\phi$.
The only nonzero component of the form $g_{\phi\rho}$ is $g_{\phi\phi}$, so the first term needs $\rho = \phi$ and the second needs $\nu = \phi$.
With $\rho = \phi$ and $\nu \in \{u,v,\theta\}$,

$$\Gamma^\phi_{u\phi} = \tfrac{1}{2}g^{\phi\phi}\,\partial_u g_{\phi\phi} = \frac{1}{2}\cdot\frac{4}{c^2(v-u)^2\sin^2\theta}\cdot\left(-\frac{c^2(v-u)\sin^2\theta}{2}\right) = -\frac{1}{v-u},$$

$$\Gamma^\phi_{v\phi} = \frac{1}{v-u},$$

$$\Gamma^\phi_{\theta\phi} = \frac{1}{2}\cdot\frac{4}{c^2(v-u)^2\sin^2\theta}\cdot\frac{c^2(v-u)^2\sin\theta\cos\theta}{2} = \frac{\cos\theta}{\sin\theta} = \cot\theta,$$

and the mirror cases $\nu = \phi$ with $\rho \in \{u,v,\theta\}$ give the same three numbers,

$$\Gamma^\phi_{\phi u} = -\frac{1}{v-u}, \qquad \Gamma^\phi_{\phi v} = \frac{1}{v-u}, \qquad \Gamma^\phi_{\phi\theta} = \cot\theta.$$

The case $\nu = \rho = \phi$ gives $\tfrac{1}{2}g^{\phi\phi}\left(\partial_\phi g_{\phi\phi} + \partial_\phi g_{\phi\phi} - \partial_\phi g_{\phi\phi}\right) = 0$.
Every other $\Gamma^\phi_{\nu\rho}$ is zero.

The connection has fifteen nonzero symbols in total in the bare $(u,v)$ chart, those of Steps 19 to 21.
In Step 22a we carry them into the rescaled chart.

## Step 22. The same symbols with the first index lowered

Lowering uses $\Gamma_{\mu\nu\rho} = g_{\mu\alpha}\Gamma^\alpha_{\nu\rho}$, and the metric of Step 16 makes each of these a single product.
Since the only nonzero $g_{u\alpha}$ is $g_{uv} = -c^2/2$, the first index $u$ picks out the $v$ symbol, and the other way round:

$$\Gamma_{u\nu\rho} = -\frac{c^2}{2}\Gamma^v_{\nu\rho}, \qquad \Gamma_{v\nu\rho} = -\frac{c^2}{2}\Gamma^u_{\nu\rho},$$

$$\Gamma_{\theta\nu\rho} = \frac{c^2(v-u)^2}{4}\Gamma^\theta_{\nu\rho}, \qquad \Gamma_{\phi\nu\rho} = \frac{c^2(v-u)^2\sin^2\theta}{4}\Gamma^\phi_{\nu\rho}.$$

Carrying the fifteen symbols of Steps 19 to 21 through:

$$\Gamma_{u\theta\theta} = -\frac{c^2}{2}\left(-\frac{v-u}{2}\right) = \frac{c^2(v-u)}{4}, \qquad \Gamma_{u\phi\phi} = \frac{c^2(v-u)\sin^2\theta}{4},$$

$$\Gamma_{v\theta\theta} = -\frac{c^2}{2}\left(\frac{v-u}{2}\right) = -\frac{c^2(v-u)}{4}, \qquad \Gamma_{v\phi\phi} = -\frac{c^2(v-u)\sin^2\theta}{4},$$

$$\Gamma_{\theta u\theta} = \Gamma_{\theta\theta u} = \frac{c^2(v-u)^2}{4}\left(-\frac{1}{v-u}\right) = -\frac{c^2(v-u)}{4}, \qquad \Gamma_{\theta v\theta} = \Gamma_{\theta\theta v} = \frac{c^2(v-u)}{4},$$

$$\Gamma_{\theta\phi\phi} = \frac{c^2(v-u)^2}{4}\left(-\sin\theta\cos\theta\right) = -\frac{c^2(v-u)^2\sin\theta\cos\theta}{4},$$

$$\Gamma_{\phi u\phi} = \Gamma_{\phi\phi u} = \frac{c^2(v-u)^2\sin^2\theta}{4}\left(-\frac{1}{v-u}\right) = -\frac{c^2(v-u)\sin^2\theta}{4}, \qquad \Gamma_{\phi v\phi} = \Gamma_{\phi\phi v} = \frac{c^2(v-u)\sin^2\theta}{4},$$

$$\Gamma_{\phi\theta\phi} = \Gamma_{\phi\phi\theta} = \frac{c^2(v-u)^2\sin^2\theta}{4}\cot\theta = \frac{c^2(v-u)^2\sin\theta\cos\theta}{4}.$$

These fifteen are the lowered symbols in the bare chart, and in Step 22a we carry them over too.

## Step 22a. Both kinds of symbol in the rescaled chart

In the chart of Step 1, whose null coordinates are $cu$ and $cv$, each symbol is multiplied by $c$ once for every upper $u$ or $v$ index and divided by $c$ once for every lower one.
For $\Gamma^\mu{}_{\nu\rho}$ that means the four with an upper null index gain a $c$, the eight with one lower null index lose one, and the three purely angular symbols are untouched:

$$\Gamma^u_{\theta\theta} = \frac{c(v-u)}{2}, \qquad \Gamma^u_{\phi\phi} = \frac{c(v-u)\sin^2\theta}{2}, \qquad \Gamma^v_{\theta\theta} = -\frac{c(v-u)}{2}, \qquad \Gamma^v_{\phi\phi} = -\frac{c(v-u)\sin^2\theta}{2},$$

$$\Gamma^\theta_{u\theta} = \Gamma^\theta_{\theta u} = -\frac{1}{c(v-u)}, \qquad \Gamma^\theta_{v\theta} = \Gamma^\theta_{\theta v} = \frac{1}{c(v-u)},$$

$$\Gamma^\phi_{u\phi} = \Gamma^\phi_{\phi u} = -\frac{1}{c(v-u)}, \qquad \Gamma^\phi_{v\phi} = \Gamma^\phi_{\phi v} = \frac{1}{c(v-u)},$$

$$\Gamma^\theta_{\phi\phi} = -\sin\theta\cos\theta, \qquad \Gamma^\phi_{\theta\phi} = \Gamma^\phi_{\phi\theta} = \cot\theta.$$

For $\Gamma_{\mu\nu\rho}$ every index is lower, so each symbol simply loses one factor of $c$ for each null index it carries.
The twelve with exactly one lose a single $c$, and the three with none keep what Step 22 gave them:

$$\Gamma_{u\theta\theta} = \frac{c(v-u)}{4}, \qquad \Gamma_{u\phi\phi} = \frac{c(v-u)\sin^2\theta}{4}, \qquad \Gamma_{v\theta\theta} = -\frac{c(v-u)}{4}, \qquad \Gamma_{v\phi\phi} = -\frac{c(v-u)\sin^2\theta}{4},$$

$$\Gamma_{\theta u\theta} = \Gamma_{\theta\theta u} = -\frac{c(v-u)}{4}, \qquad \Gamma_{\theta v\theta} = \Gamma_{\theta\theta v} = \frac{c(v-u)}{4},$$

$$\Gamma_{\phi u\phi} = \Gamma_{\phi\phi u} = -\frac{c(v-u)\sin^2\theta}{4}, \qquad \Gamma_{\phi v\phi} = \Gamma_{\phi\phi v} = \frac{c(v-u)\sin^2\theta}{4},$$

$$\Gamma_{\theta\phi\phi} = -\frac{c^2(v-u)^2\sin\theta\cos\theta}{4}, \qquad \Gamma_{\phi\theta\phi} = \Gamma_{\phi\phi\theta} = \frac{c^2(v-u)^2\sin\theta\cos\theta}{4}.$$

These are the fifteen nonzero Christoffel symbols of each kind in the rescaled chart, with the first index up and with every index lowered.
The lowering identity survives the change of chart, as it must: $\Gamma_{u\theta\theta} = g_{uv}\Gamma^v_{\theta\theta} = \left(-\tfrac{1}{2}\right)\left(-\tfrac{c(v-u)}{2}\right) = \tfrac{c(v-u)}{4}$, using the $g_{uv}$ of Step 17a rather than the one of Step 16.

## Step 23. The shape of the connection

Before computing the curvature it pays to name the pattern, because it is what makes the cancellations visible.
For the curvature, in Steps 23 to 27, we use the bare chart symbols of Steps 19 to 21 rather than the rescaled ones of Step 22a, because the pattern is cleaner without the factors of $c$ and the conclusion is that every curvature component vanishes.
That conclusion is the same in either chart: rescaling a coordinate multiplies a component by a nonzero power of $c$, which cannot turn a zero into anything else.

Let $a,b,e,f$ run over the two null indices $\{u,v\}$ and let $i,j,k,l,m$ run over the two angular indices $\{\theta,\phi\}$.
Write $w = v-u$, and attach to each null index the sign

$$\sigma_u = -1, \qquad \sigma_v = +1,$$

which is a label on the letter $u$ or $v$, not a tensor.
With that label, $\partial_a w = \sigma_a$, and let $\hat\gamma_{ij} = \mathrm{diag}\left(1, \sin^2\theta\right)$ be the metric of the unit two-sphere, so that $g_{ij} = \tfrac{c^2w^2}{4}\hat\gamma_{ij}$.

In this notation the symbols of Steps 19 to 21 are exactly these, and nothing else is nonzero:

$$\Gamma^a_{ij} = -\frac{w}{2}\sigma_a\,\hat\gamma_{ij}, \qquad \Gamma^i_{aj} = \Gamma^i_{ja} = \frac{\sigma_a}{w}\delta^i{}_j, \qquad \Gamma^i_{jk} = \hat\Gamma^i_{jk},$$

$$\Gamma^a_{bc} = 0, \qquad \Gamma^a_{bi} = 0, \qquad \Gamma^i_{ab} = 0,$$

where $\hat\Gamma^i_{jk}$ are the Christoffel symbols of the unit two-sphere, namely $\hat\Gamma^\theta_{\phi\phi} = -\sin\theta\cos\theta$ and $\hat\Gamma^\phi_{\theta\phi} = \hat\Gamma^\phi_{\phi\theta} = \cot\theta$, with the other four zero.
Check the signs against Step 19: $\Gamma^u_{\theta\theta} = -\tfrac{w}{2}(-1)(1) = \tfrac{w}{2}$ and $\Gamma^v_{\theta\theta} = -\tfrac{w}{2}(+1)(1) = -\tfrac{w}{2}$, as computed there, and against Step 20: $\Gamma^\theta_{u\theta} = \tfrac{-1}{w}\delta^\theta{}_\theta = -\tfrac{1}{w}$.

Two facts will be used repeatedly.
An upper null index forces both lower indices to be angular, and an upper angular index forbids both lower indices from being null.
And $\hat\Gamma^i_{jk}$ and $\hat\gamma_{ij}$ depend on $\theta$ alone, so $\partial_a$ annihilates them.

## Step 24. The Riemann tensor, case by case

Each index of $R^\mu{}_{\nu\rho\sigma}$ is either null or angular, giving four combinations for $(\mu,\nu)$ and, since the tensor is antisymmetric in the last pair, three for $\{\rho,\sigma\}$.
We work through all twelve cases, and every one of them is zero.
Throughout, $\lambda$ in the quadratic terms is summed over all four indices, and the rule from Step 23 decides which values contribute.

**(1) $R^a{}_{bef}$.**
The derivative terms need $\Gamma^a_{bf}$ and $\Gamma^a_{be}$, both zero by Step 23.
The quadratic terms are $\Gamma^a_{e\lambda}\Gamma^\lambda_{bf} - \Gamma^a_{f\lambda}\Gamma^\lambda_{be}$, and every $\Gamma^\lambda_{bf}$ and $\Gamma^\lambda_{be}$ with two lower null indices is zero.
The whole component is zero.

**(2) $R^a{}_{bej}$.**
The derivative terms need $\Gamma^a_{bj}$ and $\Gamma^a_{be}$, and an upper null index with a lower null index is zero, so both vanish.
The quadratic terms are $\Gamma^a_{e\lambda}\Gamma^\lambda_{bj} - \Gamma^a_{j\lambda}\Gamma^\lambda_{be}$; in the first, $\Gamma^a_{e\lambda}$ carries a lower null index against an upper null index and is zero for every $\lambda$; in the second, $\Gamma^\lambda_{be}$ has two lower null indices and is zero.
The whole component is zero.

**(3) $R^a{}_{bij}$.**
The derivative terms need $\Gamma^a_{bj}$ and $\Gamma^a_{bi}$, both zero.
The quadratic terms are $\Gamma^a_{i\lambda}\Gamma^\lambda_{bj} - \Gamma^a_{j\lambda}\Gamma^\lambda_{bi}$.
An upper null index on $\Gamma^a_{i\lambda}$ forces $\lambda = k$ angular, and $\Gamma^k_{bj} = \tfrac{\sigma_b}{w}\delta^k{}_j$, so the first term is $\tfrac{\sigma_b}{w}\Gamma^a_{ij}$ and the second is $\tfrac{\sigma_b}{w}\Gamma^a_{ji}$.
Since $\Gamma^a_{ij}$ is symmetric in $ij$, the two cancel.

**(4) $R^a{}_{jef}$.**
The derivative terms need $\Gamma^a_{jf}$ and $\Gamma^a_{je}$, which mix one angular and one null lower index against an upper null index and are zero.
The quadratic terms both contain a factor $\Gamma^a_{e\lambda}$ or $\Gamma^a_{f\lambda}$, which is zero for every $\lambda$ by the same rule.
The whole component is zero.

**(5) $R^a{}_{jek}$.**
The derivative terms are $\partial_e\Gamma^a_{jk} - \partial_k\Gamma^a_{je}$.
The second is zero because $\Gamma^a_{je}$ is zero; the first is

$$\partial_e\left(-\frac{w}{2}\sigma_a\hat\gamma_{jk}\right) = -\frac{\sigma_e\sigma_a}{2}\hat\gamma_{jk},$$

using $\partial_e w = \sigma_e$ and $\partial_e \hat\gamma_{jk} = 0$.
The quadratic terms are $\Gamma^a_{e\lambda}\Gamma^\lambda_{jk} - \Gamma^a_{k\lambda}\Gamma^\lambda_{je}$.
The first is zero for every $\lambda$.
In the second, $\Gamma^a_{k\lambda}$ forces $\lambda = l$ angular, and $\Gamma^l_{je} = \tfrac{\sigma_e}{w}\delta^l{}_j$, so

$$-\Gamma^a_{kl}\,\frac{\sigma_e}{w}\delta^l{}_j = -\left(-\frac{w}{2}\sigma_a\hat\gamma_{kj}\right)\frac{\sigma_e}{w} = +\frac{\sigma_a\sigma_e}{2}\hat\gamma_{kj}.$$

The two surviving pieces are $-\tfrac{\sigma_a\sigma_e}{2}\hat\gamma_{jk}$ and $+\tfrac{\sigma_a\sigma_e}{2}\hat\gamma_{kj}$, and $\hat\gamma$ is symmetric, so they cancel.

**(6) $R^a{}_{jkl}$.**
Antisymmetry in the last pair leaves $(k,l) = (\theta,\phi)$ as the only case, and $j$ is $\theta$ or $\phi$.

For $j = \theta$, the component is

$$R^a{}_{\theta\theta\phi} = \partial_\theta\Gamma^a_{\theta\phi} - \partial_\phi\Gamma^a_{\theta\theta} + \Gamma^a_{\theta m}\hat\Gamma^m_{\theta\phi} - \Gamma^a_{\phi m}\hat\Gamma^m_{\theta\theta},$$

where $\lambda$ has already been restricted to angular values $m$ because the upper index is null.
Now $\Gamma^a_{\theta\phi} \propto \hat\gamma_{\theta\phi} = 0$, and $\partial_\phi$ annihilates everything, and $\hat\Gamma^m_{\theta\phi}$ is nonzero only for $m=\phi$ where it multiplies $\Gamma^a_{\theta\phi} = 0$, and $\hat\Gamma^m_{\theta\theta} = 0$ for both $m$.
Every term is zero.

For $j = \phi$, the component is

$$R^a{}_{\phi\theta\phi} = \partial_\theta\Gamma^a_{\phi\phi} - \partial_\phi\Gamma^a_{\phi\theta} + \Gamma^a_{\theta m}\hat\Gamma^m_{\phi\phi} - \Gamma^a_{\phi m}\hat\Gamma^m_{\phi\theta}.$$

Term by term, with $\Gamma^a_{\phi\phi} = -\tfrac{w}{2}\sigma_a\sin^2\theta$ and $\Gamma^a_{\theta\theta} = -\tfrac{w}{2}\sigma_a$:

$$\partial_\theta\Gamma^a_{\phi\phi} = -\frac{w\sigma_a}{2}\cdot 2\sin\theta\cos\theta = -w\sigma_a\sin\theta\cos\theta,$$

$$\partial_\phi\Gamma^a_{\phi\theta} = 0,$$

$$\Gamma^a_{\theta m}\hat\Gamma^m_{\phi\phi} = \Gamma^a_{\theta\theta}\hat\Gamma^\theta_{\phi\phi} = \left(-\frac{w\sigma_a}{2}\right)\left(-\sin\theta\cos\theta\right) = +\frac{w\sigma_a}{2}\sin\theta\cos\theta,$$

$$-\Gamma^a_{\phi m}\hat\Gamma^m_{\phi\theta} = -\Gamma^a_{\phi\phi}\hat\Gamma^\phi_{\phi\theta} = -\left(-\frac{w\sigma_a}{2}\sin^2\theta\right)\cot\theta = +\frac{w\sigma_a}{2}\sin\theta\cos\theta.$$

The sum is $-w\sigma_a\sin\theta\cos\theta + \tfrac{w\sigma_a}{2}\sin\theta\cos\theta + \tfrac{w\sigma_a}{2}\sin\theta\cos\theta = 0$.

**(7) $R^i{}_{bef}$.**
The derivative terms need $\Gamma^i_{bf}$ and $\Gamma^i_{be}$, which have an upper angular index against two lower null indices and are zero.
The quadratic terms contain $\Gamma^\lambda_{bf}$ and $\Gamma^\lambda_{be}$, zero for every $\lambda$.
The whole component is zero.

**(8) $R^i{}_{bej}$.**
The derivative terms are $\partial_e\Gamma^i_{bj} - \partial_j\Gamma^i_{be}$.
The second vanishes since $\Gamma^i_{be} = 0$, and the first is

$$\partial_e\left(\frac{\sigma_b}{w}\delta^i{}_j\right) = -\frac{\sigma_b\sigma_e}{w^2}\delta^i{}_j.$$

The quadratic terms are $\Gamma^i_{e\lambda}\Gamma^\lambda_{bj} - \Gamma^i_{j\lambda}\Gamma^\lambda_{be}$.
The second vanishes since $\Gamma^\lambda_{be} = 0$.
In the first, $\Gamma^i_{e\lambda}$ forces $\lambda = k$ angular, so it is $\tfrac{\sigma_e}{w}\delta^i{}_k \cdot \tfrac{\sigma_b}{w}\delta^k{}_j = \tfrac{\sigma_e\sigma_b}{w^2}\delta^i{}_j$.
The two surviving pieces cancel.

**(9) $R^i{}_{bjk}$.**
The derivative terms are $\partial_j\left(\tfrac{\sigma_b}{w}\delta^i{}_k\right) - \partial_k\left(\tfrac{\sigma_b}{w}\delta^i{}_j\right) = 0$, because an angular derivative does not see $w$ and the Kronecker deltas are constants.
The quadratic terms are $\Gamma^i_{j\lambda}\Gamma^\lambda_{bk} - \Gamma^i_{k\lambda}\Gamma^\lambda_{bj}$.
A lower null index on $\Gamma^\lambda_{bk}$ forces $\lambda = l$ angular, giving $\hat\Gamma^i_{jl}\tfrac{\sigma_b}{w}\delta^l{}_k = \tfrac{\sigma_b}{w}\hat\Gamma^i_{jk}$ for the first and $\tfrac{\sigma_b}{w}\hat\Gamma^i_{kj}$ for the second.
The Christoffel symbols of the sphere are symmetric in their lower pair, so the two cancel.

**(10) $R^i{}_{jef}$.**
The derivative terms are

$$\partial_e\left(\frac{\sigma_f}{w}\delta^i{}_j\right) - \partial_f\left(\frac{\sigma_e}{w}\delta^i{}_j\right) = \left(-\frac{\sigma_f\sigma_e}{w^2} + \frac{\sigma_e\sigma_f}{w^2}\right)\delta^i{}_j = 0.$$

The quadratic terms are $\Gamma^i_{e\lambda}\Gamma^\lambda_{jf} - \Gamma^i_{f\lambda}\Gamma^\lambda_{je}$; in both, the upper angular index against a lower null index forces $\lambda$ angular, giving $\tfrac{\sigma_e\sigma_f}{w^2}\delta^i{}_j$ and $\tfrac{\sigma_f\sigma_e}{w^2}\delta^i{}_j$, which cancel.

**(11) $R^i{}_{jek}$.**
The derivative terms are $\partial_e\Gamma^i_{jk} - \partial_k\Gamma^i_{je} = \partial_e\hat\Gamma^i_{jk} - \partial_k\left(\tfrac{\sigma_e}{w}\delta^i{}_j\right) = 0 - 0 = 0$, the first because the sphere symbols do not depend on $u$ or $v$, the second because $w$ does not depend on the angles.
The quadratic terms are $\Gamma^i_{e\lambda}\Gamma^\lambda_{jk} - \Gamma^i_{k\lambda}\Gamma^\lambda_{je}$.
In the first, $\Gamma^i_{e\lambda}$ forces $\lambda = l$ angular, giving $\tfrac{\sigma_e}{w}\delta^i{}_l\hat\Gamma^l_{jk} = \tfrac{\sigma_e}{w}\hat\Gamma^i_{jk}$; the null value $\lambda = b$ would need $\Gamma^i_{eb}$, which is zero.
In the second, $\Gamma^\lambda_{je}$ forces $\lambda = l$ angular, giving $\hat\Gamma^i_{kl}\tfrac{\sigma_e}{w}\delta^l{}_j = \tfrac{\sigma_e}{w}\hat\Gamma^i_{kj}$; the null value $\lambda = a$ would need $\Gamma^a_{je}$, which is zero.
The two cancel by the symmetry of $\hat\Gamma^i_{jk}$.

**(12) $R^i{}_{jkl}$.**
This is the only case where something has to cancel against something else rather than term by term.
Split the sum over $\lambda$ into its angular and null parts.
The angular part is exactly the Riemann tensor of the unit two-sphere,

$$\partial_k\hat\Gamma^i_{jl} - \partial_l\hat\Gamma^i_{jk} + \hat\Gamma^i_{km}\hat\Gamma^m_{jl} - \hat\Gamma^i_{lm}\hat\Gamma^m_{jk} = \hat R^i{}_{jkl}.$$

The null part is

$$\Gamma^i_{ka}\Gamma^a_{jl} - \Gamma^i_{la}\Gamma^a_{jk} = \sum_{a \in \{u,v\}}\left[\frac{\sigma_a}{w}\delta^i{}_k\left(-\frac{w}{2}\sigma_a\hat\gamma_{jl}\right) - \frac{\sigma_a}{w}\delta^i{}_l\left(-\frac{w}{2}\sigma_a\hat\gamma_{jk}\right)\right].$$

Each bracket carries $\sigma_a^2 = 1$ and the $w$ cancels, so the sum over the two null values of $a$ gives a factor of two,

$$\Gamma^i_{ka}\Gamma^a_{jl} - \Gamma^i_{la}\Gamma^a_{jk} = 2 \cdot \left(-\frac{1}{2}\right)\left(\delta^i{}_k\hat\gamma_{jl} - \delta^i{}_l\hat\gamma_{jk}\right) = -\left(\delta^i{}_k\hat\gamma_{jl} - \delta^i{}_l\hat\gamma_{jk}\right).$$

So the whole component is

$$R^i{}_{jkl} = \hat R^i{}_{jkl} - \left(\delta^i{}_k\hat\gamma_{jl} - \delta^i{}_l\hat\gamma_{jk}\right),$$

and the two pieces are equal by the curvature of the unit two-sphere in Step 25.

## Step 25. The curvature of the unit two-sphere

The unit two-sphere has constant curvature $+1$, meaning

$$\hat R^i{}_{jkl} = \delta^i{}_k\hat\gamma_{jl} - \delta^i{}_l\hat\gamma_{jk}.$$

We verify that identity rather than quote it.
Antisymmetry in $kl$ leaves $(k,l) = (\theta,\phi)$, and two free indices give four components, of which two are computed and two are immediate.

$$\hat R^\theta{}_{\phi\theta\phi} = \partial_\theta\hat\Gamma^\theta_{\phi\phi} - \partial_\phi\hat\Gamma^\theta_{\phi\theta} + \hat\Gamma^\theta_{\theta m}\hat\Gamma^m_{\phi\phi} - \hat\Gamma^\theta_{\phi m}\hat\Gamma^m_{\phi\theta}.$$

With $\hat\Gamma^\theta_{\phi\phi} = -\sin\theta\cos\theta$, $\hat\Gamma^\theta_{\phi\theta} = 0$, $\hat\Gamma^\theta_{\theta m} = 0$ for both $m$, and $\hat\Gamma^\theta_{\phi m}$ nonzero only for $m = \phi$,

$$\hat R^\theta{}_{\phi\theta\phi} = -\left(\cos^2\theta - \sin^2\theta\right) - 0 + 0 - \left(-\sin\theta\cos\theta\right)\cot\theta = \sin^2\theta - \cos^2\theta + \cos^2\theta = \sin^2\theta,$$

which matches $\delta^\theta{}_\theta\hat\gamma_{\phi\phi} - \delta^\theta{}_\phi\hat\gamma_{\phi\theta} = \sin^2\theta - 0$.

$$\hat R^\phi{}_{\theta\theta\phi} = \partial_\theta\hat\Gamma^\phi_{\theta\phi} - \partial_\phi\hat\Gamma^\phi_{\theta\theta} + \hat\Gamma^\phi_{\theta m}\hat\Gamma^m_{\theta\phi} - \hat\Gamma^\phi_{\phi m}\hat\Gamma^m_{\theta\theta}.$$

With $\hat\Gamma^\phi_{\theta\phi} = \cot\theta$, $\hat\Gamma^\phi_{\theta\theta} = 0$, $\hat\Gamma^m_{\theta\theta} = 0$ for both $m$, and $\hat\Gamma^\phi_{\theta m}$ nonzero only for $m = \phi$,

$$\hat R^\phi{}_{\theta\theta\phi} = -\frac{1}{\sin^2\theta} - 0 + \cot\theta\cdot\cot\theta - 0 = \frac{\cos^2\theta - 1}{\sin^2\theta} = -1,$$

which matches $\delta^\phi{}_\theta\hat\gamma_{\theta\phi} - \delta^\phi{}_\phi\hat\gamma_{\theta\theta} = 0 - 1$.

The remaining two vanish on both sides.
On the right, $\delta^\theta{}_\theta\hat\gamma_{\theta\phi} - \delta^\theta{}_\phi\hat\gamma_{\theta\theta} = 0 - 0 = 0$ and $\delta^\phi{}_\theta\hat\gamma_{\phi\phi} - \delta^\phi{}_\phi\hat\gamma_{\phi\theta} = 0 - 0 = 0$.
On the left,

$$\hat R^\theta{}_{\theta\theta\phi} = \partial_\theta\hat\Gamma^\theta_{\theta\phi} - \partial_\phi\hat\Gamma^\theta_{\theta\theta} + \hat\Gamma^\theta_{\theta m}\hat\Gamma^m_{\theta\phi} - \hat\Gamma^\theta_{\phi m}\hat\Gamma^m_{\theta\theta} = 0,$$

because $\hat\Gamma^\theta_{\theta\phi}$, $\hat\Gamma^\theta_{\theta\theta}$, $\hat\Gamma^\theta_{\theta m}$ for both $m$, and $\hat\Gamma^m_{\theta\theta}$ for both $m$ are each zero, and

$$\hat R^\phi{}_{\phi\theta\phi} = \partial_\theta\hat\Gamma^\phi_{\phi\phi} - \partial_\phi\hat\Gamma^\phi_{\phi\theta} + \hat\Gamma^\phi_{\theta m}\hat\Gamma^m_{\phi\phi} - \hat\Gamma^\phi_{\phi m}\hat\Gamma^m_{\phi\theta} = 0,$$

because $\hat\Gamma^\phi_{\phi\phi} = 0$, because $\partial_\phi\cot\theta = 0$, because the only nonzero $\hat\Gamma^\phi_{\theta m}$ is $\hat\Gamma^\phi_{\theta\phi}$ and it meets $\hat\Gamma^\phi_{\phi\phi} = 0$, and because the only nonzero $\hat\Gamma^\phi_{\phi m}$ is $\hat\Gamma^\phi_{\phi\theta}$ and it meets $\hat\Gamma^\theta_{\phi\theta} = 0$.

Substituting the identity into Step 24 case (12),

$$R^i{}_{jkl} = \left(\delta^i{}_k\hat\gamma_{jl} - \delta^i{}_l\hat\gamma_{jk}\right) - \left(\delta^i{}_k\hat\gamma_{jl} - \delta^i{}_l\hat\gamma_{jk}\right) = 0.$$

All twelve cases of Step 24 are zero, so $R^\mu{}_{\nu\rho\sigma} = 0$ identically on the chart.

The cancellation is not an accident of the algebra.
The coefficient of the sphere term is $1$ because $g^{\mu\nu}\partial_\mu r\,\partial_\nu r = 2g^{uv}\partial_u r\,\partial_v r = 2\left(-\tfrac{2}{c^2}\right)\left(-\tfrac{c}{2}\right)\left(\tfrac{c}{2}\right) = 1$, that is, because $r$ measured this way really is the proper radius of the sphere it labels.

## Step 26. The cross-check by tensor transformation

Independently of Step 24, the Riemann tensor is a tensor, so it transforms as

$$R^\mu{}_{\nu\rho\sigma} = \frac{\partial \tilde x^\mu}{\partial x^a}\frac{\partial x^b}{\partial \tilde x^\nu}\frac{\partial x^c}{\partial \tilde x^\rho}\frac{\partial x^d}{\partial \tilde x^\sigma}\,R^a{}_{bcd}.$$

In the Cartesian chart every metric component is constant, so by the argument of Step 8 every Christoffel symbol vanishes there and $R^a{}_{bcd} = 0$ at every point.
The right hand side is then a sum of finite Jacobian factors times zero, so the left hand side is zero at every point of the null chart.
The transformation law reproduces the result of Step 24 without touching a single index, and the two together establish that the Riemann tensor vanishes.

## Step 27. Ricci, Ricci scalar, Kretschmann, Einstein and Weyl

With $R^\mu{}_{\nu\rho\sigma} = 0$ from Steps 24 to 26:

$$R_{\nu\sigma} = R^\mu{}_{\nu\mu\sigma} = 0, \qquad R = g^{\nu\sigma}R_{\nu\sigma} = 0,$$

$$K = R_{\mu\nu\rho\sigma}R^{\mu\nu\rho\sigma} = 0, \qquad G_{\mu\nu} = R_{\mu\nu} - \tfrac{1}{2}Rg_{\mu\nu} = 0.$$

The Weyl tensor in four dimensions is

$$C_{\mu\nu\rho\sigma} = R_{\mu\nu\rho\sigma} - \left(g_{\mu[\rho}R_{\sigma]\nu} - g_{\nu[\rho}R_{\sigma]\mu}\right) + \tfrac{1}{3}R\,g_{\mu[\rho}g_{\sigma]\nu},$$

every term of which is built from Riemann, Ricci or $R$, all of which vanish, so $C_{\mu\nu\rho\sigma} = 0$.
All four tensors vanish, and $R = 0$, $K = 0$ for the two scalars, in both null systems.
This is the flat spacetime it started as, written in a chart where the metric is no longer constant, and no change of chart can give a tensor that vanishes in one chart a nonzero component in another.

## Step 28. The geodesic equations

Take $\ddot x^\mu = -\Gamma^\mu_{\nu\rho}\dot x^\nu \dot x^\rho$ with the fifteen symbols of Step 22a, and remember that the sum over $\nu$ and $\rho$ runs over both orderings, which is what turns each pair of equal mixed symbols into a factor of two.
The equations are taken in the same chart as those symbols, so $\dot u$ and $\dot v$ are the rates of the chart coordinates $cu$ and $cv$ along the curve, matching the indices the symbols carry.

For $\mu = u$, the only contributions are $\Gamma^u_{\theta\theta}$ and $\Gamma^u_{\phi\phi}$:

$$\ddot u = -\frac{c(v-u)}{2}\dot\theta^2 - \frac{c(v-u)\sin^2\theta}{2}\dot\phi^2.$$

For $\mu = v$, the same two symbols with the opposite sign:

$$\ddot v = \frac{c(v-u)}{2}\dot\theta^2 + \frac{c(v-u)\sin^2\theta}{2}\dot\phi^2.$$

For $\mu = \theta$, the contributions are $\Gamma^\theta_{u\theta}$ and $\Gamma^\theta_{\theta u}$ together, $\Gamma^\theta_{v\theta}$ and $\Gamma^\theta_{\theta v}$ together, and $\Gamma^\theta_{\phi\phi}$:

$$\ddot\theta = \frac{2}{c(v-u)}\dot u\dot\theta - \frac{2}{c(v-u)}\dot v\dot\theta + \sin\theta\cos\theta\,\dot\phi^2 = -\frac{2\left(\dot v - \dot u\right)}{c(v-u)}\dot\theta + \sin\theta\cos\theta\,\dot\phi^2.$$

For $\mu = \phi$, the contributions are the four mixed symbols $\Gamma^\phi_{u\phi}$, $\Gamma^\phi_{\phi u}$, $\Gamma^\phi_{v\phi}$, $\Gamma^\phi_{\phi v}$ and the pair $\Gamma^\phi_{\theta\phi}$, $\Gamma^\phi_{\phi\theta}$:

$$\ddot\phi = -\frac{2\left(\dot v - \dot u\right)}{c(v-u)}\dot\phi - 2\cot\theta\,\dot\theta\dot\phi.$$

These four are the geodesic equations in spherical null coordinates.

## Step 29. Consistency with the spherical chart

In the spherical chart of Minkowski spacetime the geodesic equations are $\ddot t = 0$, $\ddot r = r\dot\theta^2 + r\sin^2\theta\,\dot\phi^2$, and the two angular equations with $-\tfrac{2}{r}\dot r$ in front.
Those hold in the same chart, so the $t$ in them is the chart coordinate $ct$, and the null chart coordinates are $cu = ct - r$ and $cv = ct + r$.
Differentiating $cu = ct - r$ twice along the curve gives $\ddot u = \ddot t - \ddot r$ with every dot now a chart rate, so

$$\ddot u = 0 - \left(r\dot\theta^2 + r\sin^2\theta\,\dot\phi^2\right) = -r\dot\theta^2 - r\sin^2\theta\,\dot\phi^2,$$

and $r = c(v-u)/2$ turns that into the $\ddot u$ equation of Step 28 exactly.
The same step with $cv = ct + r$ gives the $\ddot v$ equation.
For the angles, $\dot r = \tfrac{1}{2}\left(\dot v - \dot u\right)$ and $r = \tfrac{c}{2}(v-u)$ give

$$\frac{2\dot r}{r} = \frac{2 \cdot \tfrac{1}{2}\left(\dot v - \dot u\right)}{\tfrac{c}{2}(v-u)} = \frac{2\left(\dot v - \dot u\right)}{c(v-u)},$$

which is the coefficient standing in front of $\dot\theta$ and $\dot\phi$ in Step 28.
The two coordinate systems describe the same geodesics, as they must.

---

## What the machine checks

`_tools/derivations/verify_metrics.py` builds the metric of every coordinate system from its stated line element, and then computes the inverse metric, the Christoffel symbols with the first index up and with every index lowered, Riemann, Ricci, the Ricci scalar, Kretschmann, Einstein and Weyl in sympy.
It reads every stated component, translates the LaTeX back into an expression and asserts equality, and it also requires every other component to vanish, which catches an omission as well as a wrong number.
It weights each component into the $x^0 = cT$ chart of Step 1 before comparing, which is the arithmetic of Steps 7a, 17a and 22a done by machine, and it checks each geodesic equation against $\ddot x^\mu + \Gamma^\mu_{\nu\rho}\dot x^\nu\dot x^\rho = 0$ built from the symbols it computed.
Anything it cannot parse, or cannot finish simplifying inside its time budget, it reports as UNCHECKED and names, rather than passing it in silence.
It exits non-zero and names what disagreed.

It needs sympy, in a virtual environment outside the repository:

    python3 -m venv /tmp/mfs-venv && /tmp/mfs-venv/bin/pip install sympy
    /tmp/mfs-venv/bin/python _tools/derivations/verify_metrics.py

Both null systems pass it.
Pass `--system minkowski/spherical_null` to check just one.
