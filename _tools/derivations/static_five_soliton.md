# Static soliton of five dimensional gravity

## Sources

Dietrich Kramer's original article, *Axialsymmetrische stationäre Lösungen der projektiven Feldtheorie*, appears in Acta Physica Polonica B 2(6), 807-811 (1971).
The scan preserved by INSPIRE, record 71611, identifies the year on page 807 and the static spherical solution in equation (22), page 810.
Kayll Lake's *Static Ricci-flat 5-manifolds admitting the 2-sphere*, Classical and Quantum Gravity 23, 5871-5881 (2006), supplies a geometrical treatment of the full neutral class.
His reference [10] says 1970; the original scan establishes 1971.
The original articles of Chodos and Detweiler (1982), Gross and Perry (1983), and Davidson and Owen (1985) have verified metadata, with full text still owed.

## Charts

Write $f=1-r_s/r$, with $r_s>0$ a length, and take a spacelike coordinate $y$.
With signature $(-,+,+,+,+)$ the family is

$$ds_5^2=-f^a c^2dt^2+f^{-a-b}dr^2+r^2f^{1-a-b}d\Omega^2+f^bdy^2,\qquad a^2+ab+b^2=1.$$

In Kramer's equation (22), $a=A-B$, $b=2B$, and $r_s=2m$.
His scalar $S=e^{2\Phi}=f^B$ supplies the extra term $S^2dy^2$.
The restriction $A^2+3B^2=1$ is therefore the ellipse above.
Lake's equations (14)-(17) use $h=\sqrt f$, $\delta=2a$, $\alpha=2b$, and $C=1/r_s$.
The isotropic radius $\rho$ obeys $r=\rho(1+r_s/4\rho)^2$, with $h=(4\rho-r_s)/(4\rho+r_s)$.
Lake's equation (19) uses its reciprocal; the simultaneous change $(\alpha,\delta)\mapsto(-\alpha,-\delta)$ in his equation (21) relates the conventions.

## Parameter range and limits

The whole ellipse is retained, including negative exponents and the exceptional points.
The exterior is $r>r_s$; $r$ is not the areal radius, which is $R=rf^{(1-a-b)/2}$.
As $r\to\infty$ the geometry is locally Minkowski space, or asymptotically a product with a circle if $y$ is periodic.
The point $(a,b)=(1,0)$ is the Schwarzschild black string.
The point $(a,b)=(0,1)$ is time times Euclidean Schwarzschild; smooth closure of the circle requires period $4\pi r_s$ for $y$.
At the other points the inner endpoint is a curvature singularity, as Lake's equations (26)-(32) establish.
The negative exceptional exponents have a singular inner endpoint on this exterior branch.
Reduction along $y$ with $g_E=f^{b/2}g_4$ yields a massless scalar solution with exponent $\gamma=a+b/2$, so $\gamma^2=1-3b^2/4$.

## Drawings

The drawings take $r_s=1$, $a=8/7$, and $b=-3/7$, a generic singular member of the ellipse.
The symmetry spheres shrink to zero radius at the inner endpoint, while the fifth circle expands.
Each radial null plane fixes $y=0$, $\theta=\pi/2$, and $\phi=0$.
Its optical radius satisfies $dr_*/dr=f^{-13/14}$ and $r_*(r_s)=0$.
Writing $w=f^{1/14}$ makes this integral $14r_s\int_0^w(1-v^{14})^{-2}dv$.
The hypergeometric antiderivative is finite at $r=r_s$, so the radial conformal quotient is a triangle with a timelike singular left boundary.
The compactification uses $p,q=\arctan((ct\mp r_*)/16r_s)$.

The spatial equator has metric $f^{-5/7}dr^2+r^2f^{2/7}d\phi^2$.
Its circles have radius $R=rf^{1/7}$.
The Euclidean height satisfies $(dz/dr)^2=f^{-12/7}[5r_s/(7r)-36r_s^2/(49r^2)]$.
It lies level at $r=36r_s/35$, where $h=1/6$.
Inside that circle the ambient signature is Minkowski, with $(dZ/dr)^2=-(dz/dr)^2$.
The two parts meet with one tangent plane.
The inner surface continues toward the singular endpoint; the numerical profile ends at $h=10^{-4}$, whose circle radius is $0.072r_s$.
The outer surface ends at $r=4r_s$ and continues to spatial infinity.
Both parts are retained, and their heights are checked against independent quadratures of the slice.
The fixed figure window contains both parts through rotation.

The curvature figure evaluates the stored Kretschmann scalar.
It compares the generic member to the coincident curves of the regular black string and bubble, $K=12r_s^2/r^6$.
For the generic member $K$ grows as $(r-r_s)^{-18/7}$ near the inner endpoint.

## Algebra on the parameter ellipse

The rational parametrisation $a=(u^2-1)/(u^2+u+1)$ and $b=-u(u+2)/(u^2+u+1)$ covers a dense subset of the full ellipse.
The excluded value $(1,-1)$ is its limit as $u\to\infty$.
Polynomial remainders modulo $a^2+ab+b^2-1$ reduce the tensor coefficients while leaving radial powers intact.
The numerical source checks compare the unreduced Ricci tensor with zero and each radial transformation with Kramer's spherical form, including both regular members and negative exponents.
Positive exterior factors are restored before numerical substitution where the canonical algebra splits an isotropic ratio into signed irreducible factors.

## Checks

The dedicated checks compare the Einstein reduction with the massless scalar geometry, the bubble circumference with proper radius, and the symmetry sphere's turning point with the full family.
The scoped page check follows completion of the tensors, diagrams, and reciprocal relations.
