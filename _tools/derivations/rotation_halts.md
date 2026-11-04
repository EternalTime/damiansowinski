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

Pending the source equations.

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

Pending chart selection.

## Curvature

Pending stored scalars.

## Verification

Pending scoped checks.
