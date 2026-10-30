# Early braneworlds

## Sources and scope

Keiichi Akama presented Pregeometry at Nara in 1982; Springer published the proceedings in 1983, in Lecture Notes in Physics 176, pages 267 to 271.
The chapter DOI is 10.1007/3-540-11994-9_41, and hep-th/0001113 reproduces its text with a cover added in 2000.
The six dimensional bulk is flat; the Higgs and gauge fields form a vortex, and the effective gravitational action concerns the gently curved four dimensional core.
Akama's equation (4) introduces tubular coordinates around an unspecified embedding, so it specifies no unique additional curved bulk geometry.

Valery Rubakov and Mikhail Shaposhnikov's two papers in Physics Letters B 125 (1983) concern different constructions.
The paper on pages 136 to 138, DOI 10.1016/0370-2693(83)91253-4, concerns field theoretic trapping on a scalar wall in a flat five dimensional background.
The paper on pages 139 to 143, DOI 10.1016/0370-2693(83)91254-6, concerns six dimensional gravity with a cosmological constant.
The wall's fermion modes follow Rubakov's own derivation in hep-ph/0104152, section 3.1, equations (2) to (5).
The six dimensional gravitational solution follows George Lavrelashvili's explicit reconstruction in arXiv:1002.4566, equations (1) to (8), also published in Proceedings of A. Razmadze Mathematical Institute 159 (2012), pages 65 to 74.

Matt Visser's 1985 model, DOI 10.1016/0370-2693(85)90112-1, is reproduced in hep-th/9910093.
The equations on PDF pages 5 to 7 specify an electric Einstein-Maxwell solution with a cosmological constant, with a varying lapse and flat spatial slices.
The 1999 note is distinct from the 1985 body.

Euan J. Squires's 1986 paper, DOI 10.1016/0370-2693(86)90347-3, concerns vacuum Einstein geometries with a cosmological constant and scalar trapping.
Peter Anthony Bannister's 1987 Durham thesis, The cosmological constant and dimensional reduction, derives the three dimensional example in chapter 3, equations (3) to (6), and the higher dimensional family in chapter 4, equations (27) to (28).
Bannister's regular axis is a specified member of the family, selected in equations (35) to (36).

Gary Gibbons and David Wiltshire's 1987 Einstein-Maxwell membranes, DOI 10.1016/0550-3213(87)90125-8, are reproduced in hep-th/0109093.
The general membrane is their equations (3.1) to (3.6); the flat membrane with vanishing cosmological constants is (3.7).
The six dimensional case has a four dimensional membrane and two transverse directions.
The four dimensional specialization of their arbitrary dimensional construction is Melvin's magnetic universe.

## Source geometries

We take the mostly positive signature and write the time coordinate as $x^0 = ct$ when forming components.

Akama's Cartesian six dimensional background and the polar coordinates of the vortex are

$$ds^2=-c^2dt^2+dx^2+dy^2+dz^2+d\rho^2+\rho^2d\phi^2.$$

The scalar wall background has one transverse direction,

$$ds^2=-c^2dt^2+dx^2+dy^2+dz^2+dw^2.$$

The wall profile affects the fermion equation; gravitational backreaction is neglected in this construction.

For the six dimensional gravity solution, let $L=\sqrt{8/(5\Lambda)}$ in Lavrelashvili's notation, and set $u=\rho/L$.
His equations (7) and (8), after reversing the overall signature, read

$$ds^2=\cos^{4/5}u\,(-c^2dt^2+dx^2+dy^2+dz^2)+d\rho^2+L^2\sin^2u\cos^{-6/5}u\,d\phi^2.$$

The radial interval is $0\leq\rho<\pi L/2$, with $\phi$ of period $2\pi$.
The endpoint at $\pi L/2$ has a vanishing four dimensional warp factor and diverging angular circumference.
Lavrelashvili's perturbations with angular harmonic $m\geq2$ include an unstable scalar mode.

Visser's explicit electric solution is

$$ds^2=-\cosh^2(w/L)c^2dt^2+dx^2+dy^2+dz^2+dw^2.$$

In his gravitationally normalized units the electric field is $E=1/L$, the potential is $A_t=\sinh(w/L)$, and his field equations are $G_{AB}=\lambda g_{AB}+T_{AB}$ with $\lambda=E^2/2$.
Thus the cosmological constant in $G_{AB}+\Lambda g_{AB}=T_{AB}$ is $\Lambda=-E^2/2$.
His more general lapse-only ansatz leaves its supporting stress unspecified until these field equations are imposed.

Squires's family in Bannister's notation has $n$ tangential coordinates and one proper radial coordinate $w$.
With $\alpha^2=\lambda n/4$, its scale factors are

$$a_A(w)=\sinh^{1/n+q_A}(\alpha w)\cosh^{1/n-q_A}(\alpha w),\qquad \sum_Aq_A=0,\qquad\sum_Aq_A^2=1-1/n.$$

The sign of the time factor is negative and the other signs are positive.
The three dimensional member has $a_0=\cosh(\alpha w)$ and $a_1=\sinh(\alpha w)$, chapter 3 equations (6.1) and (6.2).
With the periodic spatial coordinate scaled for a regular axis, this is global anti-de Sitter space in three dimensions.
Bannister's regular axis in five dimensions has three factors $\cosh^{1/2}(\alpha w)$ and one factor $\sinh(\alpha w)\cosh^{-1/2}(\alpha w)$.
It is the five dimensional anti-de Sitter soliton under an explicit radial change of coordinate.
The general exponent family carries the two constraints above; it is not an arbitrary common warp factor.

For the flat six dimensional magnetic membrane, set $H=1+\rho^2/a^2$.
Gibbons and Wiltshire's (3.7) reads

$$ds^2=H^{2/3}(-c^2dt^2+dx^2+dy^2+dz^2+d\rho^2)+\rho^2H^{-2}d\phi^2,$$

$$F=B_0\rho H^{-2}d\rho\wedge d\phi,\qquad a^2=\frac{2}{3\pi G B_0^2}.$$

Their normalization is $\kappa^2=4\pi G$ in equation (2.1), and must be retained or explicitly converted when checking the Maxwell stress.
The transverse area is infinite and radial proper distance is unbounded.
Their general membrane uses the radial warp coordinate $\beta$ and the function

$$\Delta(\beta)=\bar\lambda+\frac{2Gk}{\beta^{d-3}}-\frac{8\pi GB_0^2}{(d-2)(d-3)\beta^{2(d-3)}}-\frac{2\Lambda\beta^2}{(d-1)(d-2)}.$$

The line element is $ds^2=\beta^2d\bar s^2+d\beta^2/\Delta+\Delta d\tilde\phi^2$ after reversing the overall signature, with the membrane Einstein metric carrying the corresponding reversed signature.
Its domains depend on the roots of $\Delta$ and the angular period required at a regular axis.
The flat member, curved membrane, and bulk cosmological constant cases retain their distinct assumptions.

## Localization and stability

Akama's low energy fields are suppressed outside the vortex; his Einstein action is induced by fluctuations, equations (27) to (29).
Rubakov's fermion zero mode is proportional to $\exp[-h\int_0^w\varphi_c(s)\,ds]$ times a four dimensional chiral spinor, hep-ph/0104152 equation (5).
Visser's massive geodesics have transverse turning points, and the scalar wave equation has the explicit bound spectrum on PDF page 7.
Gibbons and Wiltshire find flux controlled chiral fermion zero modes in equations (5.2) to (5.11).
For their flat magnetic membrane with both cosmological constants zero, the axisymmetric neutral scalar spectrum extends continuously to arbitrarily small positive mass, equations (6.1) and (6.2).
The mass gap in the curved cases depends on the boundary conditions discussed after (6.9); the constant scalar zero mode is not square integrable in that discussion.
Their stability argument concerns radial perturbations of the flat magnetic membrane, section 4; section 7 leaves stability of the curved extensions open.

## Coordinates and normalizations

The two vortex coordinate systems are the Cartesian background and its transverse polar coordinates, with $u=\rho\cos\phi$ and $v=\rho\sin\phi$.
The scalar wall retains its five dimensional flat background.
The six dimensional vacuum uses the two scale factors $f=\cos^{2/5}(\rho/L)$ and $h=L\sin(\rho/L)\cos^{-3/5}(\rho/L)$.
Its source cosmological constant has the positive value $8/(5L^2)$ after converting to $G_{AB}+\Lambda g_{AB}=0$.

Visser's free lapse and his electric solution are separate coordinate systems.
The electric solution has $R_{tt}=\cosh^2(w/L)/L^2$, $R_{ww}=-1/L^2$, and flat ordinary spatial directions.
His field normalization is absorbed into the electric field as stated above.

The Squires family has five dimensions, four tangential scale factors, and the explicit constraints $\sum q_A=0$ and $\sum q_A^2=3/4$.
The scale factors and their first and second derivatives appear explicitly in the tensors, with the source definitions above.
The three dimensional member uses its regular angular coordinate, with the spatial coordinate in Bannister's chapter 3 rescaled by $L$.
For the five dimensional regular axis, set the anti-de Sitter soliton's radius to $2L$, its tip radius to $2L$, its areal coordinate to $r=2L\sqrt{\cosh(w/L)}$, and its circle coordinate to $\tau=L\phi$.
This pulls the soliton's five dimensional metric back to the regular axis form.

In the general magnetic coordinates, $\beta$ is dimensionless, the angle $\phi$ has period $2\pi$, and $\tilde\phi=C\phi$.
We use $\mu=2Gk$ and $Q^2=2\pi G B_0^2/3$, so $D=K+\mu/\beta^3-Q^2/\beta^6-\Lambda\beta^2/10$ has inverse squared length units and $C$ has squared length units.
The flat membrane sets $K=0$; the curved membrane uses a maximally symmetric four dimensional Einstein factor of sectional curvature $K$ in its static spherical coordinates.
Gibbons and Wiltshire also permit other Einstein factors with the same Ricci tensor.
Intervals require $D>0$, and a simple regular axis at $\beta_0$ requires $C=2/|D'(\beta_0)|$.
The Maxwell two form is $F=CB_0\beta^{-4}d\beta\wedge d\phi$ with coupling $8\pi G$.
In six dimensions the trace reversed Einstein equation is $R_{AB}=\Lambda g_{AB}/2+8\pi G(F_{AC}F_B{}^C-g_{AB}F_{CD}F^{CD}/8)$.

For the flat magnetic solution, the change $\beta=(1+\rho^2/a^2)^{1/3}$, $\mu=Q^2=4/(9a^2)$, $C=3a^2/2$, and $\Lambda=K=0$ relates the two radial coordinates.
The scalar radial kinetic and mass measures both reduce to $\rho\,d\rho$.
Thus the angular harmonic zero has the Bessel equation and a continuum down to arbitrarily small positive mass; the constant function has divergent norm.

## Drawings

Every coordinate system has one radial light section with its fixed coordinates stated in the caption.
All views are stills.
The free lapse uses Visser's electric solution, and the general Squires view uses its regular axis exponents.
The general flat magnetic view uses $\mu=Q^2=4/9$, $C=3/2$, and $\Lambda=0$, in the interval $\beta>1$.
The curved magnetic section has $K=1/4$, $\mu=2$, $Q=1/2$, $C=1$, and $\Lambda=-1/10$, at $r=1/2$.
Its null curves belong to that fixed radius section; they are not asserted to be six dimensional null geodesics.

The magnetic transverse surface at $a=1$ is embedded from its regular axis to $\rho=5a$, with circle radius $\rho/(1+\rho^2/a^2)$ and its full circumference retained.
The conformal triangle describes the totally geodesic magnetic radial section, with $p,q=\arctan((ct\mp\rho)/a)$, and makes no assertion of a conformal completion of all six dimensions.
The embedding's $t=0$ slice is carried onto both magnetic radial charts by the coordinate map above and onto the conformal section by these null coordinates.
The other light sections belong to different localization backgrounds and carry no slice of this magnetic surface.
The curvature line uses the stored magnetic Kretschmann scalar alone.
