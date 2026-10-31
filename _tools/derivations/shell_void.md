# A cosmological void surrounded by a shell

## Sources and scope

Kei-ichi Maeda and Humitaka Sato treated the expansion of a compressed shell in a dust universe in 1983, doi:10.1143/PTP.70.772.
Nobuyuki Sakai and Paul Haines recover the homogeneous exterior when their McVittie mass vanishes, arXiv:astro-ph/9909183, section 2 and equation (3-1).
Norimasa Sugiura, Ken-ichi Nakao, Daisuke Ida, Nobuyuki Sakai, and Hideki Ishihara take a Minkowski interior, an exterior filled with dust, and a dust shell, arXiv:astro-ph/9912414, sections 2 and 3.1.
The construction considered here has an empty interior and zero cosmological constant.

## The two smooth regions

The interior is Minkowski spacetime in spherical coordinates.
The exterior is a Friedmann universe filled with dust.
The shell has intrinsic metric $ds_\Sigma^2=-d\tau^2+R(\tau)^2d\Omega^2$ and surface stress $S^i{}_j=\operatorname{diag}(-\sigma,0,0)$ in units $c=1$.
The radius $R$ evolves under the junction conditions.

## The shell equations

The normal points from the empty interior toward the cosmological exterior.
With $K_{ij}=h_i{}^a h_j{}^b\nabla_a n_b$ and $[K]=K_+-K_-$, the junction equation is $[K^i{}_j]=-8\pi G(S^i{}_j-\delta^i{}_j S^k{}_k/2)$.
These are equations (2-8) through (2-10) of Sakai and Haines, with the surface stress of (2-18) at zero surface pressure.
The smooth tensors on either side and the distribution supported on the shell must be treated separately.

## The flat dust exterior

In units $G=c=1$, take exterior cosmic time $t$, $a=(t/t_0)^{2/3}$, $H=2/(3t)$, and dust density $\rho=1/(6\pi t^2)$.
The shell has exterior comoving radius $r=R/a$ and peculiar speed $v=a\,dr/dt$.
Define $x=HR$, $w=x+v=dR/dt$, $\gamma=(1-v^2)^{-1/2}$, and $B=\sqrt{1+2xv+x^2}$.
Induced metric continuity fixes $dT/dt=B$ for the interior Minkowski time $T$, and $d\tau/dt=1/\gamma$.
The identity $B^2-w^2=1-v^2$ ensures that $-dT^2+dR^2=-d\tau^2$ on the shell.

The angular curvatures on its two sides are

$$K^\theta{}_{\theta,-}=\frac{\gamma B}{R},\qquad K^\theta{}_{\theta,+}=\frac{\gamma(1+xv)}{R}.$$

Their jump fixes the surface density:

$$\sigma=\frac{\gamma[B-(1+xv)]}{4\pi R}=\frac{H^2R}{4\pi\gamma[B+1+xv]}.$$

The rationalized expression is positive for the expanding ordinary shell and avoids subtracting nearly equal quantities when $HR$ is small.
The proper shell mass is $m=4\pi R^2\sigma$.
Its Newtonian limit is the removed background mass $4\pi\rho R^3/3$.

## Evolution and accretion

For surface dust the temporal jump is $[K^\tau{}_\tau]=4\pi\sigma$.
The normal stress equation is $-\sigma\overline{K^\tau{}_\tau}=\rho\gamma^2v^2$, with the bar denoting the mean of the two limiting values.
Consequently $K^\tau{}_{\tau,+}=2\pi\sigma-\rho\gamma^2v^2/\sigma$.
The homogeneous exterior has $K^\tau{}_{\tau,+}=\gamma^3dv/dt+\gamma vH$, so the motion obeys

$$\frac{dR}{dt}=HR+v,\qquad\frac{dv}{dt}=\frac{-\gamma vH+2\pi\sigma-\rho\gamma^2v^2/\sigma}{\gamma^3},\qquad\frac{dT}{dt}=B.$$

These follow directly from the junction equations (2-8)-(2-10), the stress (2-18), and the normal flux (2-20) of Sakai and Haines at zero McVittie mass and zero surface pressure.
The angular constraint fixes $\sigma$ at each moment.
Tangential conservation independently requires

$$\frac{d\sigma}{dt}+\frac{2\sigma}{R}\frac{dR}{dt}=\rho\gamma v,\qquad\frac{dm}{dt}=4\pi R^2\rho\gamma v.$$

The influx raises the proper mass when $v>0$.

## Initial data for the drawings

Take $t_0=1$, $a(t_0)=1$, $R(t_0)=0.15$, $v(t_0)=0$, and $T(t_0)=0$ in the chosen time unit.
The initial size is $H_0R_0=0.1$, the surface density is approximately $0.002645984$, and the interior speed is approximately $0.09950372$.
The shell starts with positive peculiar acceleration $dv/dt\simeq0.01662521$.
The small void power law $R\propto t^{(15+\sqrt{17})/24}$ is an asymptotic Newtonian result; the trajectory obeys the relativistic differential equations above.

## Bulk and surface curvature

Minkowski's bulk Riemann tensor vanishes.
The dust exterior has its ordinary smooth Friedmann curvature.
The singular Einstein tensor on the shell is $8\pi S_{ab}\delta(n)$, where $n$ is signed Gaussian normal distance and $S_{ab}=\sigma u_au_b$ for shell velocity $u^a$.
The bulk scalar invariants are defined away from the shell.
Products of the shell's delta distribution do not define an ordinary curvature squared scalar.
Sugiura and colleagues, section 2.3 and equations (2.18)-(2.27), derive finite optical jumps by integrating the singular curvature through the shell.

## Numerical motion and the junction check

The drawing interval is $1\le t/t_0\le4$.
DOP853 integrates the relativistic equations at relative tolerance $2\times10^{-12}$ and absolute tolerance $2\times10^{-14}$, with maximum step $0.025t_0$.
At $t=2t_0$ the shell has $R=0.24387876888560137t_0$, $v=0.007872607582882827$, and $T=1.0043869706139898t_0$.
At $t=4t_0$ it has $R=0.4106499137590503t_0$, $v=0.010403560061993835$, and $T=3.011244482443109t_0$.

The normal covector in each radial plane is $n_a=\sqrt{-\det h}(-u^r,u^0)$, with $h$ its two dimensional metric and $u^a$ the shell velocity.
Contracting its acceleration with this normal determines $K^\tau{}_\tau=n_aa^a$.
The angular curvature follows independently from $K^\theta{}_\theta=-n_a\Gamma^a_{\theta\theta}/g_{\theta\theta}$.
The metrics and their Christoffel symbols determine these quantities on each side.
At $t/t_0=1,2,4$, both jumps have residuals smaller than $10^{-15}/t_0$.
The dedicated junction checks use additional time samples, differentiate the algebraic proper shell mass to check its influx, and perturb the acceleration to confirm that the temporal jump rejects an altered trajectory.

## Spatial slice and global causal boundary

The equatorial spatial slice meets the shell at $t=2t_0$ outside and $T=1.0043869706139898t_0$ inside.
Each side has intrinsic metric $d\rho^2+\rho^2d\phi^2$, with $\rho=r$ inside and $\rho=ar$ outside.
Their circumferences and radial distances agree at $\rho=R$.
The time normal of this piecewise spacelike slice has a corner at the shell, while its intrinsic geometry is a plane.
A circle marks the shell on that plane.

The freely chosen dust Friedmann exterior and shell initial data describe a family with different global causal boundaries.
The numerical initial value solution specifies a finite future interval and leaves its past extension free.
A global conformal diagram is therefore recorded in NOT_DRAWN, with no guessed cosmological or shell boundary.

## Surface Weyl curvature

The smooth Minkowski and Friedmann bulks are conformally flat.
The singular Riemann tensor has $R_{nanb}^{\mathrm{sing}}=-[K_{ab}]\delta(n)$ in Gaussian normal coordinates.
For dust, $[K^\tau{}_\tau]=4\pi\sigma$ and $[K^\theta{}_\theta]=-4\pi\sigma$.
Taking its traces leaves a nonzero singular Weyl tensor, so conformal flatness of the bulks is not a conformal flatness claim about the glued spacetime.
Its curvature figure samples only the smooth exterior Kretschmann scalar and marks the shell trajectory separately.

## The proper radial coordinate of the source

Sakai and Haines equations (2-1) and (2-2), at zero McVittie mass, use $ds^2=-dt^2+a^2[d\chi^2+f(\chi)^2d\Omega^2]$.
Here $\chi$ and the curvature sign $k$ are dimensionless, $a$ has dimensions of length, and $f$ is $\sin\chi$, $\chi$, or $\sinh\chi$ for $k=1,0,-1$.
The shell position is $\chi_s(t)$, with $R=a f(\chi_s)$.
This coordinate covers the closed exterior through its equator to its antipodal centre.
The radial coordinate $r=f(\chi)$ covers a hemisphere in the closed case and obeys $dr^2/(1-kr^2)=d\chi^2$ on that hemisphere.
A constant rescaling makes $r$ a length and $a$ dimensionless in the alternative comoving coordinate.
The flat specialization of the proper radial coordinate has $a=t_0(t/t_0)^{2/3}$ and horizontal plotting distance $t_0\chi$.
The two flat exterior diagrams represent the same shell in their respective radial coordinates.

The interior null coordinates used to check the ray families are $(T+r)/R_0$ and $(T-r)/R_0$, with $R_0=0.15t_0$ the initial cavity radius.
This normalization makes each coordinate dimensionless on the scale of the cavity.
Their constant-coordinate curves are the Minkowski null rays in the interior.
