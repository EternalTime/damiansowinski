# McLenaghan-Tariq-Tupper electrovacuum

## Sources and access

Raymond G. McLenaghan and Nessim Tariq published their solution in Journal of Mathematical Physics 16, 2306-2312 (November 1975), DOI 10.1063/1.522461.
N. Tariq and B. O. J. Tupper published a class of algebraically general non-null electrovacua in General Relativity and Gravitation 6, 345-360 (August 1975), DOI 10.1007/BF00761969.
Both original full texts remain unread; their primary bibliographic records were fetched.
Tupper's sequel, General Relativity and Gravitation 7, 479-486 (June 1976), DOI 10.1007/BF00766405, acknowledges the independent McLenaghan-Tariq solution on its publisher record.
No quotation or equation is attributed to an unread original.

Alexey Anatolievich Magazev's 2021 article, DOI 10.1007/s11040-021-09385-3, is available as arXiv:1908.10035v3.
Section 5, equation (29), supplies an explicit chart, equations (31)-(33) its invariant frame and coframe, and equation (42) the integrated Hamilton-Jacobi equation.
The entire preprint was read.
Anatoly Leonidovich Evtyagin and Alexey Anatolievich Magazev's 2021 article, DOI 10.24147/1812-3996.2021.26(1).21-26, was read in the university-hosted journal issue, pages 21-26.
Its Section 2 states the difference between the metric and electromagnetic symmetries, and Sections 3-4 treat the Klein-Gordon equation.

Ian Anderson and Charles Torre's 2022 article, DOI 10.1007/s10714-022-02913-8, was read as arXiv:2109.02558v3.
Its introduction identifies the McLenaghan-Tariq electrovacuum as homogeneous and non-inheriting, and equations (2.4)-(2.14) explain duality rotations and symmetry inheritance.
Its new spacetime has a positive cosmological constant and is a distinct solution.
Its global hyperbolicity and geodesic completeness results are not transferred to this electrovacuum.

## Chart definitions

Magazev's equation (29), with the overall sign changed to the site's convention, is

\[
ds^2=-(dX^1+dX^2-kX^3dX^4)^2+(dX^2)^2+e^{-kX^2}(dX^3)^2+e^{kX^2}(dX^4)^2,
\qquad k>0.
\]

The substitution

\[
ct=X^1+X^2,\quad r=X^2,\quad Y=X^3,\quad Z=X^4,\quad a=2/k
\]

produces the exponential form

\[
ds^2=-(c\,dt-2Y\,dZ/a)^2+dr^2+e^{-2r/a}dY^2+e^{2r/a}dZ^2.
\]

Then

\[
x=a e^{r/a},\quad y=Y,\quad \varphi=Z/a
\]

produces the upper-half-plane form

\[
ds^2=-(c\,dt-2y\,d\varphi)^2+\frac{a^2}{x^2}(dx^2+dy^2)+x^2d\varphi^2.
\]

This last form is also explicitly displayed by Frausto, Kolář, Málek and Torre (2025), equation (4.16), with attribution to the textbook's equation (12.21).
The symbol \(\varphi\) is a real coordinate; no periodic identification is imposed by these transformations.
The changes of coordinates agree at the numerical points used for the Einstein-Maxwell equations.
The History and the reciprocal relation were approved on 2026-10-30.

## Maxwell field and inherited motions

Write \(T=ct\) and take the oriented orthonormal coframe
\[
 e^0=dT-2y\,dz/a,\qquad e^1=dr,\qquad e^2=e^{-r/a}dy,\qquad e^3=e^{r/a}dz.
\]
Its exterior derivatives are \(de^0=-2e^2\wedge e^3/a\), \(de^1=0\), \(de^2=-e^1\wedge e^2/a\), and \(de^3=e^1\wedge e^3/a\).
The Einstein tensor in this frame is \(2\operatorname{diag}(1,-1,1,1)/a^2\).
In the normalization \(G_{ab}=2(F_{ac}F_b{}^c-g_{ab}F_{cd}F^{cd}/4)\), a compatible field is
\[
 F=\frac{\sqrt2}{a}\left[\cos(2r/a)e^0\wedge e^1+\sin(2r/a)e^2\wedge e^3\right].
\]
The orientation obeys \({}^{\star}(e^0\wedge e^1)=-e^2\wedge e^3\).
Both \(dF\) and \(d{}^{\star}F\) vanish by the four exterior derivatives above.
A constant duality phase may be added.
The two electromagnetic invariants cannot vanish together, since their squared sum is \(16/a^4\).

The metric Killing fields in the exponential coordinates are
\[
 K_0=\partial_T,\quad K_1=\partial_z,\quad K_2=\partial_y+2z\partial_T/a,\quad
 K_3=\partial_r+(y/a)\partial_y-(z/a)\partial_z.
\]
All four preserve the coframe.
The first three preserve \(F\); the fourth obeys \(\mathcal L_{K_3}F=-(2/a){}^{\star}F\).
Thus the connected inherited subgroup of this transitive group is three dimensional.
Its generators obey \([K_1,K_2]=2K_0/a\), while \([K_3,K_1]=K_1/a\) and \([K_3,K_2]=-K_2/a\).
The full motion group and the field's inherited subgroup are distinct.
The stored scalar values are \(R=0\) and \(K=64/a^4\).

## Tensor computation

All three charts obey the same Einstein-Maxwell equations.
The source and both source-free Maxwell equations are checked at three numerical points per chart.
Each chart is also compared at those points with the exponential metric pulled back by its stated map.
The original-coordinate names are \((u,r,y,z)=(X^1,X^2,X^3,X^4)\).
The extracted Evtyagin-Magazev matrix is not used.

## Diagram slices

Every drawing uses \(a=1\).
The spacetime views follow the totally geodesic surface \(y=z=0\), with both coordinates real, and its images in the other charts.
Its induced metric is \(-dT^2+dr^2\).
The conformal diamond compactifies this surface alone, using \(p,q=\arctan((T\mp r)/a)\).
No global causal-boundary claim follows from this surface.
The spatial surface \(t=z=0\) has metric \(dr^2+e^{-2r/a}dy^2\), a hyperbolic plane.
An isometric immersion of its strip of width \(2\pi a\) along \(y\), for \(0\le r\le2a\), is Beltrami's pseudosphere with radius \(a e^{-r/a}\).
The strip is rolled once about the drawing's axis; the coordinates of spacetime retain their real domains.
Its radial ends \(r=0\) and \(r=2a\) lie at \(x=a\) and \(x=ae^2\) in the upper half plane.
The drawing is a finite immersed strip, with no global embedding or periodic spacetime identification asserted.
The curvature figure evaluates the stored scalars on \(-2a\le r\le2a\); both are constant, so one still view suffices.
