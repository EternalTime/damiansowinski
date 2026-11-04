# Rotation halts cylindrical collapse

## Claims to verify

1. The primary record identifies the authors, title, date, and journal fields.
2. The nonzero particle angular momentum imposes a lower circumference bound of 8 pi alpha Lambda, with the source's definitions and hypotheses.
3. The shell's two counter-rotating dust streams determine its surface stress and both junction conditions.
4. The interior's flatness needs its precise scope checked against the general dynamical solution.
5. The non-rotating collapse and the four related geometries have distinct, sourced physical connections.

## Sources and conventions

We use $G = c = 1$.
The primary records of Apostolatos and Thorne (1992), DOI 10.1103/PhysRevD.46.2435, and Echeverria (1993), DOI 10.1103/PhysRevD.47.2271, agree with Crossref and INSPIRE.
The two founding full texts remain unread: APS PDFs returned HTTP 403, APS Harvest fulltext returned HTTP 401, and the CaltechAUTHORS fulltext returned HTTP 403.
The bound is verified against the founding abstract, without a quotation.
The bulk and shell equations are read in Gleiser and Ramirez, arXiv:1106.3122v1, equations (1)-(13) and (15), which attribute the model to Apostolatos and Thorne.
Kurita and Nakao, arXiv:1112.4252, equations (2.13)-(2.15), independently state the stress and the junction convention.
The regular interior of a dynamical cylindrical shell can contain waves.
Flatness belongs to static equilibrium and to the restricted dynamical family in section III B of Gleiser and Ramirez.

## Charts and tensor checks

The canonical coordinates on the two vacuum regions are $(t,\rho,\phi,z)$, with independent times on the two sides and $\rho=R$ at the shell.
They are equation (1) of Gleiser and Ramirez, with their $r$ renamed $\rho$ and with the angular coordinate before the axial coordinate.
The junction coordinates $(\tau,\phi,z)$ are their equation (4).
Flat interior and static exterior are special choices of the functions in these coordinates.
Every tensor is retained with the functions free; the vacuum differential constraints are stated with those functions.
The bulk checks evaluate two independent vacuum profiles at rational points, and the interior check also evaluates the flat specialization.
The intrinsic shell metric is three dimensional, but the bulk spacetime is four dimensional.

## Bound and shell matching

The canonical bulk metric is
$$ds^2 = e^{2(\gamma-\psi)}(-dt^2+dr^2)+r^2e^{-2\psi}d\phi^2+e^{2\psi}dz^2.$$
The vacuum equations are $\psi_{,tt}=\psi_{,rr}+\psi_{,r}/r$, $\gamma_{,r}=r(\psi_{,r}^2+\psi_{,t}^2)$, and $\gamma_{,t}=2r\psi_{,r}\psi_{,t}$.
Write $V=dR/d\tau$, $X_\pm=dt_\pm/d\tau=\sqrt{e^{2\psi_\Sigma-2\gamma_\pm}+V^2}$, and $\psi_{,n}^\pm=X_\pm\psi_{,r}^\pm+V\psi_{,t}^\pm$.
Continuity of the first fundamental form requires $\psi_+=\psi_-=\psi_\Sigma$ and
$$h=-d\tau^2+R^2e^{-2\psi_\Sigma}d\phi^2+e^{2\psi_\Sigma}dz^2.$$
For $u=\alpha e^{\psi_\Sigma}/R$, the two dust streams have orthonormal stress $S_{ab}=\operatorname{diag}(\sigma,p,0)$, where $\sigma=\lambda\sqrt{1+u^2}/(2\pi R)$ and $p=\lambda u^2/(2\pi R\sqrt{1+u^2})$.
With an outward normal on both faces, the full Israel condition is $[K_{ab}]=-8\pi(S_{ab}-h_{ab}S/2)$.
Its two independent jumps are
$$[\psi_{,n}]=-\frac{2\lambda}{\sqrt{R^2+\alpha^2e^{2\psi_\Sigma}}},\qquad [X]=-\frac{4\lambda\sqrt{R^2+\alpha^2e^{2\psi_\Sigma}}}{R}.$$
These are Gleiser and Ramirez equations (5)-(8), with their $J$ renamed $\alpha$.

A regular axis has $\gamma_-(t,0)=0$.
The radial vacuum constraint implies $\gamma_-\geq0$, and therefore $X_-\leq\sqrt{e^{2\psi_\Sigma}+V^2}$.
Future oriented timelike evolution in this canonical region has $X_+>|V|$.
Consequently
$$\frac{4\lambda\sqrt{R^2+\alpha^2e^{2\psi_\Sigma}}}{R}=X_--X_+<\sqrt{e^{2\psi_\Sigma}+V^2}-|V|\leq e^{\psi_\Sigma}.$$
Since $\sqrt{R^2+\alpha^2e^{2\psi_\Sigma}}\geq\alpha e^{\psi_\Sigma}$ for $\alpha>0$, we obtain $R>4\alpha\lambda$.
The circumference $C=2\pi R e^{-\psi_\Sigma}$ and proper axial rest mass density $\Lambda=\lambda e^{-\psi_\Sigma}$ then obey $C>8\pi\alpha\Lambda$ at regular points.
The weak inequality states the lower bound without claiming that its endpoint is attained.
The turning circumference also depends on initial data.

For the restricted flat interior, $\psi_-=\gamma_-=0$, $X_-=\sqrt{1+V^2}$, and the jumps and equation (9) imply
$$\ddot R=\frac{(1+V^2)\alpha^2}{R(R^2+\alpha^2)}-\frac{\lambda\sqrt{1+V^2}(R^2+2\alpha^2)^2}{R^2(R^2+\alpha^2)^{3/2}}.$$
This reproduces equation (15) independently from the general shell equations.
Equation (14) in the preprint has a different expression for $X_+$ from equation (6); the general junction expression above is retained and the reduced equation is rederived.
An admissible local bounce has $\alpha=1$, $\lambda=0.08$, $R=2$, and $V=0$.
Here $X_+=1-0.16\sqrt5>0$, $\ddot R=0.1-0.72/5^{3/2}>0$, and $R>0.32$.
Its exterior Cauchy data must satisfy the moving-boundary compatibility conditions; a static logarithmic exterior is valid only at equilibrium.

## Drawings and moments

The radial interior drawing follows a local bounce of the restricted flat-interior family, with $\alpha=1$, $\lambda=0.08$, $R(0)=2$, and $V(0)=0$.
The shell trajectory is integrated in proper time, including $dt_-/d\tau=\sqrt{1+V^2}$.
Its conversion to interior time is interpolated and independently checked against the integrated trajectory at 701 interleaved points within $2\times10^{-12}$.
The axial shell drawing uses the same trajectory and the induced proper-time metric.
Both null families on these planes preserve $t\pm\rho$ or $\tau\pm z$.
The exterior drawing is the exact equilibrium at $R_0=2$, $\alpha=1$, $\lambda=\sqrt5/18$, with $\kappa=1/2$ and $\gamma_0=2\ln(3/2)$.
Its radial null families preserve $t\pm\rho$.
It is a different member from the bouncing example.

The turning interior disc and the material cylinder are embedded separately at their common turning event.
On the radial interior drawing, the disc is the interval $0\leq\rho\leq2$ and the cylinder intersects it at $\rho=2$.
On the axial shell drawing, the disc intersects at $z=0$ and the cylinder is the interval $-2\leq z\leq2$.
Every intersection names its embedding view, surface, and moment through the shared slice registry.
The equilibrium exterior surface reaches $2\leq\rho\leq5$ and appears only on the equilibrium exterior drawing.
Its circles have circumferential radius $\rho\sqrt{\rho/2}$.
The three embeddings satisfy the stored spatial metrics along, across, and around their surfaces.
Only the existing cover, meridian, r, shell, mark, K, and front styles are painted.
The displayed source text counts are 204, 205, and 198 for the spacetime views, 158, 186, and 175 for the embeddings, and 230 for the curvature view, using view_text.length.
There is one moment per embedding view and no movie.

NOT_DRAWN: a global conformal diagram of the dynamical family requires complete wave data and a globally specified shell history.
Those data are free in the sourced equations, so a unique global causal completion is undetermined.
NOT_DRAWN: a dynamical exterior ray view for the turning example requires compatible exterior Cauchy data beyond the shell jets.
The equilibrium exterior is the exact exterior example drawn here.
No exterior waveform is inferred from the shell trajectory alone.

## Curvature

The curvature curves read the stored Kretschmann scalars of the two bulk charts, with the equilibrium exterior functions and the flat interior specialization.
The exterior scalar is positive and falls as $\rho^{-7}$ for $\kappa=1/2$.
The regular bulk Ricci scalar vanishes on both sides.
The exterior and interior curves stop on opposite sides of $\rho=2$.
NOT_DRAWN: a bulk Kretschmann scalar on the ideal shell would square its delta distribution, which has no canonical distributional value.
The shell stress is verified through the extrinsic-curvature jump.

## Verification

The dedicated numerical check reads the stored metric and Christoffel symbols on both faces.
It contracts the outward normal with the shell acceleration and angular connection components, and compares all three orthonormal components with the two-stream stress.
Five restricted dynamical states and two equilibria pass at 80 decimal digits, with a maximum restricted Israel-jump residual of $7.9066\times10^{-82}$.
A separate state with nonzero interior wave derivatives, $\psi_\Sigma=0.2$ and $\gamma_-=0.1$, satisfies the full jump with a residual of $1.2519\times10^{-81}$.
The flat-interior specialization has zero Riemann components; an independent regular cylindrical vacuum wave has nonzero Riemann components.
Admissible random shell states obey the strict lower bound, and candidate states at or below it fail the positive-lapse condition.
