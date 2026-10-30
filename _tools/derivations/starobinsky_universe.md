# Starobinsky universe

## Semiclassical geometry

The 1980 construction is a family of homogeneous isotropic solutions of semiclassical field equations with local one loop vacuum polarization.
Its evolving scale factor is defined by a differential equation.
The constant-curvature de Sitter solutions are special members, given explicitly in Eq. (7) of Starobinsky (1980).
The curvature of the evolving family varies with time.

## Primary sources

Starobinsky, A.A., A new type of isotropic cosmological models without singularity, Physics Letters B 91 (1980), 99-102, doi:10.1016/0370-2693(80)90670-X.
The author-uploaded full-text HTML is at https://www.researchgate.net/publication/222625284_A_New_Type_of_Isotropic_Cosmological_Models_Without_Singularity.
Mukhanov and Chibisov (1981), JETP Letters 33, 532-535, Eq. (1), state the same covariant equation and Eq. (2) states the conformal de Sitter metric.
Hirano (2026), arXiv:2610.00401, Appendix A, Eq. (13), independently states the exact flat homogeneous constraint.
Yoshida and Quintin (2018), doi:10.1088/1361-6382/aacf4b, Secs. II-III, derive the past affine length and the parallel transported curvature diagnostic.

## Homogeneous constraint

Let $\tau=ct$, $h=a'/a$, $Q=h^2+k/a^2$, with primes denoting $d/d\tau$.
The scale factor $a$ is dimensionless, the comoving radius is a length, and $k$ has dimensions $L^{-2}$.
The constants $H_0$ and $M$ are inverse lengths.
The semiclassical constraint is

$$Q-\frac{Q^2}{H_0^2}+\frac{1}{M^2}\left(2hh''+6h^2h'-(h')^2-\frac{2kh^2}{a^2}+\frac{k^2}{a^4}\right)=0.$$

For the flat expanding branch, $x=(h/H_0)^{3/2}$ and $y=h'/(2H_0^2x^{1/3})$ obey Starobinsky's Eq. (12),

$$y\frac{dy}{dx}=\frac{M^2}{12H_0^2}(x-x^{-1/3})-y.$$

The decreasing phase curve A departs from the de Sitter saddle $(x,y)=(1,0)$.

## Exact equations and approximations

The field equations and the phase-plane equation are exact within the stated semiclassical model.
The perturbations of Eqs. (8)-(10) retain the first order in their small amplitude.
The oscillatory scale factor in Eq. (14) is a late-time asymptotic formula.
The local scalaron reduction neglects $Q^2/H_0^2$ when $|Q|/H_0^2\ll1$.
In a flat expanding phase with $h^2\gg M^2$, leading slow roll gives $h'\simeq-M^2/6$, as in Hirano (2026), Eq. (5).
Integrating this approximate relation gives $a\simeq a_0\exp[-M^2(\tau_s-\tau)^2/12]$ for $M(\tau_s-\tau)\gg1$.
The Gaussian scale factor is approximate.

## Past boundary

The flat branch approaching $a\propto e^{H_0\tau}$ has finite past null affine length $\int_{-\infty}^{\tau}a\,d\tau$.
Its leading small-amplitude departure has exponent $q=\sqrt{9/4+M^2/H_0^2}-3/2$, with $h-H_0\propto a^q$.
The tidal diagnostic $h'/a^2\propto a^{q-2}$ diverges for a nonzero leading mode with $0<M^2/H_0^2<10$.
Scalar curvature invariants approach their de Sitter values in the same limit.
This is an application of the geometric diagnostic of Yoshida and Quintin; their scalar-field example uses the Einstein frame.

## Chart normalization

The cosmic and conformal charts use a dimensionless scale factor, a radial length and $k$ with dimensions $L^{-2}$.
The closed and open exact charts use a dimensionless radial coordinate $\chi$, so their spatial metrics include $H_0^{-2}$.
Their dimensionless scale factors are $\cosh(H_0\tau)$ and $\sinh(H_0\tau)$, with $k=H_0^2$ and $k=-H_0^2$ respectively.
The flat exact chart uses $a=\exp(H_0\tau)$ and $k=0$.
The three choices obey $Q=D=H_0^2$, where $D=a''/a$.
The two evolving charts retain every scale-factor derivative in their tensors, without a slow-roll reduction.

The $R^2$ contribution to the covariant 00 equation is $18(Q^2-D^2)+6hR'$.
Dividing by 18 gives the higher-derivative bracket in the homogeneous constraint.
The dedicated test checks this normalization independently of the chart printer.
It also recovers the founding phase-plane equation by substituting $h=H_0x^{2/3}$ and $h'=2H_0^2x^{1/3}y$.
The coefficient of $h^2h'$ is 6.
The coefficient 2 printed in Ellis, Kaloper, Olive and Yokoyama (1999), Eq. (5), disagrees with the covariant expression, the founding phase-plane equation and Hirano's Eq. (13), and is not used.

## Regular trace evolution

Write $S=R'$, $R=6(h'+2h^2+k/a^2)$ and $Q=h^2+k/a^2$.
The regular equations are

$$a'=ah,\qquad h'=R/6-2h^2-k/a^2,\qquad R'=S,$$
$$S'=-3hS-M^2R+\frac{12M^2}{H_0^2}Q(R/6-Q).$$

The constraint multiplied by $M^2$ is

$$C=M^2(Q-Q^2/H_0^2)+QR/3-R^2/36+hS/3.$$

The trace system obeys $C'=-4hC$ for all three spatial curvatures, with $(k/a^2)'=-2hk/a^2$.
Initial data with $C=0$ therefore satisfy the homogeneous constraint throughout the evolution.
The dedicated test checks this polynomial identity and rejects a changed anomaly coefficient as a control.
No equation in the numerical evolution is divided by $h$.

## Numerical realization

The illustrative trajectory uses $M/H_0=1/2$, without a claim of observational parameter fitting.
Set $H_0=1$ for the integration and define $u=h-1$.
The growing unstable manifold near the saddle has the second-order expansion

$$h'=qu+bu^2+O(u^3),\qquad b=\frac{q(4q+3)}{6(q+1)}.$$

The finite initial displacement is $u=-10^{-4}$.
The initial $R$ follows from this approximate slope, and the initial $S$ is fixed by the exact constraint.
The second-order slope is an initialization approximation; the subsequent differential equation retains the anomaly term and every higher derivative.
Time and the scale factor are normalized by $\tau=0$ and $a=1$ at the crossing $h=1/\sqrt2$.
DOP853 integrates through $H_0\tau=160$ with relative tolerance $2\times10^{-12}$, absolute tolerance $2\times10^{-14}$ and maximum late-time step $0.2$.
Conformal time is integrated from that same crossing, with $\eta=0$ there, avoiding subtraction of the enormous early-time conformal integral.
Its numerical inverse is a cubic Hermite interpolation with exact slope $d\tau/d\eta=a$ and cosmic-time sample spacing $0.004$.
The scoped check compares a displacement halved and tolerances quartered after the common normalization.

## Figures and causal limits

The spacetime diagrams trace the published metric on a radial plane at $\theta=\pi/2$ and $\phi=0$.
The numerical scale factor and its first two derivatives are shared by the evolving rays, conformal-map checks and curvature samples.
The conserved radial null coordinates are $\eta\pm r$ for the flat branches and $H_0\eta\pm\chi$ for the curved exact branches.
The exact primitives are $H_0\eta=\arctan\sinh(H_0\tau)$, $-\exp(-H_0\tau)$ and $\ln\tanh(H_0\tau/2)$ for the closed, flat and open de Sitter charts.

The evolving flat branch has $\eta\to-\infty$ from its past exponential asymptotic and $\eta\to+\infty$ from the sourced late-time $a\propto\tau^{2/3}$ asymptotic.
Its compactification is the whole Minkowski-shaped radial triangle under $p,q=\arctan[H_0(\eta\mp r)]$.
The past edge is labeled as a null tidal boundary; it is not labeled as ordinary complete null infinity.
The de Sitter special cases use their own full conformal domains and regular past continuations.
No geometry is truncated at a reheating time.

The flat spatial equator at the reference epoch is the plane $z=0$, with circumferential radius $\rho=r$.
The whole closed de Sitter equator at $\tau=0$ has $\rho=H_0^{-1}\sin\chi$ and $z=-H_0^{-1}\cos\chi$.
Each embedding view contains one moment, so no movie or stack registration is needed.
The curvature figure evaluates the stored Ricci and Kretschmann scalars on the same numerical trajectory from $H_0\tau=-6$ to 80.
Its panels use linear scales and no new frame kind.
