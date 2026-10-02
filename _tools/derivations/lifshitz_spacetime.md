# Lifshitz spacetime

The six charts of `lifshitz_spacetime.json` are written by `print_charts.py --metric lifshitz_spacetime`, which computes every tensor from each line element with the checker's own `Geometry` and reads every printed value back before it writes.
This note records what each chart is and where it comes from, what the checks hold it to, and what the diagrams draw.

## Step 1. Kachru, Liu and Mulligan's chart

Kachru, Liu and Mulligan's metric, their (2.1), is

$$ds^2 = L^2\left(-r^{2z}dt^2 + r^2d\mathbf{x}^2 + \frac{dr^2}{r^2}\right), \qquad 0 < r < \infty ,$$

with every coordinate a pure number, since they set the Planck length to 1.
It is the most general metric with their scaling $t \to \lambda^zt$, $\mathbf{x} \to \lambda\mathbf{x}$, $r \to r/\lambda$, their (2.2), together with translations, rotations and the reflections of space and time, once the scale is taken for the radial coordinate.
They work in four dimensions, $\mathbf{x} = (x, y)$, and so does the entry.
The entry gives the coordinates their dimensions, as Copsey and Mann's (1.2) does and as the Poincaré chart of `anti_de_sitter` does: $t \to ct/L$, $\mathbf{x} \to \mathbf{x}/L$ and $r \to r/L$, so that

$$ds^2 = -\left(\frac{r}{L}\right)^{2z}c^2dt^2 + \frac{r^2}{L^2}\left(dx^2 + dy^2\right) + \frac{L^2}{r^2}dr^2 .$$

Every component is taken in the chart $x^0 = ct$.

## Step 2. What supports it

The mixed Einstein tensor is diagonal and constant, $G^\mu{}_\nu = \mathrm{diag}(3,\ z^2 + z + 1,\ z^2 + z + 1,\ 2z + 1)/L^2$, the Ricci scalar is $-2(z^2 + 2z + 3)/L^2$, Kachru, Liu and Mulligan's value, and the Kretschmann scalar is $4(z^4 + 2z^2 + 3)/L^4$.
Kachru, Liu and Mulligan support the metric with a 2-form and a 3-form field strength joined by a topological coupling, their (2.3) to (2.11), with $\Lambda = -(z^2 + z + 4)/2L^2$.
Taylor's massive vector field does the same, and Copsey and Mann show the two actions dual, their (1.7) to (1.13).
`lifshitz_check` uses the vector, in the normalisation of Copsey and Mann's (1.3) to (1.6) at $d = 4$:

$$G_{ab} + \Lambda g_{ab} = \tfrac{1}{2}\left(F_{ac}F_b{}^c - \tfrac{1}{4}g_{ab}F^2\right) + \tfrac{1}{2}m^2\left(A_aA_b - \tfrac{1}{2}g_{ab}A^2\right), \qquad \nabla_aF^{ab} = m^2A^b ,$$

with $A = q\,(r/L)^z\,c\,dt$, $q^2 = 2(z - 1)/z$ and $m^2 = 2z/L^2$.
Both equations hold slot by slot, and $q$ is real exactly when $z \ge 1$.
On the null vectors along $r$ and along $x$ the Einstein tensor gives $2(z - 1)/L^2$ and $(z - 1)(z + 2)/L^2$, so the null energy condition holds exactly for $z \ge 1$, which `_tools/test_lifshitz.py` holds on the published file.
The Weyl tensor is $z(z - 1)$ times a tensor that does not vanish, so the spacetime is conformally flat only at $z = 1$, where it is anti-de Sitter space.

## Step 3. The inverse radius

With $u = L^2/r$,

$$ds^2 = -\left(\frac{L}{u}\right)^{2z}c^2dt^2 + \frac{L^2}{u^2}\left(dx^2 + dy^2 + du^2\right) ,$$

Kachru, Liu and Mulligan's (3.1) in Lorentzian signature, Hořava and Melby-Thompson's (2.6) and Keeler, Knodel and Liu's (1.1).
At $z = 1$ it is the published Poincaré chart of `anti_de_sitter` with $u$ for its $z$, which the check holds component by component.
The letter $z$ is the exponent on this page, so the coordinate is $u$, as in Kachru, Liu and Mulligan.

## Step 4. The proper distance

With $r = Le^{\rho/L}$,

$$ds^2 = -e^{2z\rho/L}c^2dt^2 + e^{2\rho/L}\left(dx^2 + dy^2\right) + d\rho^2 ,$$

the domain wall coordinates of Taylor's (2.1), with her exponents $\alpha_t = z$ and $\alpha_x = 1$.
Koroteev and Libanov's braneworld of December 2007 has the same form in five dimensions, $e^{-2\xi k|z|}dt^2 - e^{-2\zeta k|z|}d\mathbf{x}^2 - dz^2$ with an anisotropic fluid in the bulk, their solution $a = \xi k|z|$, $b = \zeta k|z|$; the name Lifshitz scaling entered that paper in its version of December 2009.

## Step 5. The tortoise coordinate

On the plane of $t$ and $r$ a light ray has $c\,dt = \pm L^{z+1}dr/r^{z+1}$, so with

$$w = \frac{L^{z+1}}{z\,r^z}, \qquad ds^2 = \frac{L^2}{z^2w^2}\left(-c^2dt^2 + dw^2\right) + \left(\frac{L}{zw}\right)^{2/z}\left(dx^2 + dy^2\right) .$$

It is the gauge $A = C$ of Keeler, Knodel and Liu's (3.4), which they call a tortoise coordinate.
The plane of $t$ and $w$ is that of anti-de Sitter space's Poincaré chart with the radius $L/z$.

## Step 6. The two null charts

The advanced time $v = ct - w$ is constant along each ingoing ray of fixed $x$ and $y$, and

$$ds^2 = -\left(\frac{r}{L}\right)^{2z}dv^2 + 2\left(\frac{r}{L}\right)^{z-1}dv\,dr + \frac{r^2}{L^2}\left(dx^2 + dy^2\right) ,$$

the Eddington-Finkelstein chart of Keränen, Keski-Vakkuri and Thorlacius's (9) with $b = 1$, their Lifshitz vacuum.
Along such a ray $g_{vr}\,dr$ is the step of an affine parameter, so $s = (2L/z)(r/L)^z$ is one, and

$$ds^2 = -\frac{z^2s^2}{4L^2}dv^2 + dv\,ds + \left(\frac{zs}{2L}\right)^{2/z}\left(dx^2 + dy^2\right) ,$$

Copsey and Mann's (2.16) with their $\tau = v$ and $u = s$.
Its block of $v$ and $s$ is finite and invertible at $s = 0$, a null surface reached at a finite affine parameter, their (2.7), and $g_{xx}$ vanishes there.
Each chart after the first is checked to be the first pulled back, $J^{\mathsf T}gJ$.

## Step 7. The singularity

Every scalar built from the Riemann tensor is constant.
Along the timelike geodesic of energy $E$ with no momentum along $x$ or $y$, $u = (E(L/r)^{2z},\ 0,\ 0,\ -(r/L)\sqrt{E^2(L/r)^{2z} - 1})$ in the chart $x^0 = ct$, the tidal tensor $E_{ab} = R_{acbd}u^cu^d$ has the component $(1 + (z - 1)E^2(L/r)^{2z})/L^2$ along $x$ in the falling frame, Copsey and Mann's (2.12), which grows without bound toward $r = 0$ for $z \neq 1$.
The proper time to $r = 0$ is finite: at $z = 2$, $E = 1$ and $L = 1$ it is $\tfrac{1}{2}\arcsin r_0^2$ from $r_0$.
`conformal.py` checks both on the published Riemann tensor and Christoffel symbols before it draws the edges as a singularity, and `_tools/test_lifshitz.py` holds the component for a general $z$.

## Step 8. What the diagrams draw

Every diagram is drawn at $z = 2$, Kachru, Liu and Mulligan's own case, and $L = 1$.

The spacetime diagrams draw the plane of the time and the radial coordinate in each chart, where $w = L^3/2r^2 = u^2/2L = (L/2)e^{-2\rho/L} = L^2/2s$ and a ray keeps $ct \mp w$.
The null charts are drawn against $v - r$ and $v - s/2$, which are time functions where $r < 2^{1/3}L$ and $s < \sqrt{2}\,L$: $\tau = v - ar$ has $|\nabla\tau|^2 = -2a(L/r)^{z-1} + a^2r^2/L^2$.
The inverse radius chart adds three planes of $t$ and $x$ at fixed depth, whose null lines have the slope $(L/u)^{z-1}$.
They are no geodesics for $z > 1$: on them $\Gamma^u{}_{tt}\dot t^2 + \Gamma^u{}_{xx}\dot x^2 = -(z - 1)(L/u)^{2z}u\,\dot t^2/L^2$, so light launched along one is turned toward larger $u$.

The figure follows that turning.
A null geodesic with energy $E$ and momentum $p$ along $x$ has $\dot u^2 = (u/L)^2\left(E^2(u/L)^{2z} - p^2(u/L)^2\right)$ in units with $c = 1$, Keeler, Knodel and Liu's effective potential (1.3), and turns back where $(u/L)^{z-1} = p/E$.
At $z = 2$ its path is $dx/du = p/\sqrt{E^2u^2/L^2 - p^2}$, the catenary $u = u_0\cosh((x - x_0)/u_0)$ with $u_0 = Lp/E$, and a ray that leaves the depth $u_s$ at the angle $\alpha$ from the straight way to the boundary, in the frame of an observer at rest, has $u_0 = u_s\sin\alpha$.
At $z = 1$ the spacetime is conformally flat and the same rays are straight lines of the half plane, which all reach the boundary.
The boundary lies at an infinite affine distance, so the rays are integrated in $\sigma$ with $d\lambda = d\sigma/u^2$.

The conformal diagram is Minkowski's triangle by $p, q = \arctan((ct \mp w)/L)$, the boundary on $X = 0$ and the two null edges $r = 0$, Copsey and Mann's Figure 1: the Poincaré patch with singularities where its horizons would be.

The embedding diagram is the moment $t = 0$, $y = 0$ of the proper distance chart, $d\rho^2 + e^{2\rho/L}dx^2$ for every $z$, with a strip of width $2\pi L$ along $x$ rolled up: the circle at $\rho$ has the radius $Le^{\rho/L}$, which grows more slowly than the distance out for $\rho < 0$, the pseudosphere with $dZ/d\rho = -\sqrt{1 - e^{2\rho/L}}$ measured down from the rim, and faster for $\rho > 0$, where the strip is drawn in Minkowski space with $dZ/d\rho = \sqrt{e^{2\rho/L} - 1}$.
The geometry is static, so there is no movie.
