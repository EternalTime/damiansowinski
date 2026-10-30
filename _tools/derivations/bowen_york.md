# Bowen-York initial data

## A maximal spatial slice

The initial data are a Riemannian metric $\gamma_{ij}$ and an extrinsic curvature $K_{ij}$ on a spatial hypersurface, with $K = \gamma^{ij}K_{ij} = 0$.
In geometrized units $G = c = 1$, the conformal decomposition is $\gamma_{ij} = \psi^4\delta_{ij}$ and $K_{ij} = \psi^{-2}\widehat A_{ij}$, with $\psi > 0$.
Brandt and Brügmann, Phys. Rev. Lett. 78, 3606 (1997), equations (1)-(3), distinguish these physical and conformal tensors.
All curvature of $\gamma_{ij}$ is intrinsic curvature in three dimensions.

## The momentum constraint

For $r = |\boldsymbol{x}|$ and $n^i = x^i/r$, the conformal tensor is

$$\widehat A^{ij} = \frac{3}{2r^2}\left(P^i n^j + P^j n^i - (\delta^{ij} - n^i n^j)P^k n_k\right) + \frac{3}{r^3}\left(\epsilon^{ikl}S_k n_l n^j + \epsilon^{jkl}S_k n_l n^i\right).$$

Here $P^i$ and $S^i$ are constant linear and angular momentum parameters, and all indices on $\widehat A^{ij}$ are raised and lowered with the Euclidean metric $\delta_{ij}$.
Brandt and Brügmann's equation (5) states this tensor, and equations (2) and (7) state its divergence constraint and its superposition over punctures.
The Hamiltonian equation is

$$\nabla^2\psi = -\frac{1}{8}\psi^{-7}\widehat A_{ij}\widehat A^{ij},$$

where $\nabla^2$ is the Euclidean Laplacian.

## Punctures

Brandt and Brügmann's equations (8)-(10) set $\psi = \sum_A m_A/(2r_A) + u$, with $u\to1$ at spatial infinity and $u$ of class $C^2$ at each puncture.
Here $r_A = |\boldsymbol{x} - \boldsymbol{x}_A|$, $\boldsymbol{x}_A$ is the puncture position, and $m_A>0$ is its bare mass parameter.
Each puncture is another asymptotically flat end, as their equations (11)-(14) show.
Their individual momentum and spin parameters agree with individual ADM quantities in the limit of large separation, as the text following equation (7) specifies.

## Inversion-symmetric throats

Cook, Living Rev. Relativ. 3, 5 (2000), section 3.2.1, equations (67)-(75), states the original throat construction with its assumptions and boundary conditions.
For a single throat of coordinate radius $a>0$, the inversion is $r\mapsto a^2/r$ and the conformal factor obeys $\psi(r,\theta,\phi) = (a/r)\psi(a^2/r,\theta,\phi)$.
At the throat this implies $\partial_r\psi = -\psi/(2a)$, and asymptotic flatness requires $\psi\to1$ at infinity.
Cook's equation (72) adds the image tensor

$$\mp\frac{3a^2}{2r^4}\left(P_i n_j + P_j n_i + (\delta_{ij} - 5n_i n_j)P^k n_k\right)$$

to the momentum and spin tensor of equation (70).
The sign labels the two inversion parities of the conformal tensor.
The extension to several throats uses tensor images and the corresponding boundary condition at each throat.
The 1980 publisher abstract verifies the original inversion and two-ended geometry; the detailed equations here follow Cook's readable review.

## Spatial coordinate systems

Cartesian coordinates $(x,y,z)$ have $d\ell^2 = \psi^4(dx^2+dy^2+dz^2)$ and $\psi = \psi(x,y,z)$, with each puncture or excised throat removed from the domain as its construction requires.
Brandt and Brügmann state Cartesian and spherical coordinates explicitly after equation (5).
In spherical coordinates $(r,\theta,\phi)$, $d\ell^2 = \psi^4[dr^2+r^2(d\theta^2+\sin^2\theta\,d\phi^2)]$, with $\psi = \psi(r,\theta,\phi)$ in the general single-hole case.
For the axisymmetric spinning hole, Gleiser, Nicasio, Price, and Pullin, Phys. Rev. D 57, 3401 (1998), section II, equations (2)-(8), use the spherical Euclidean background and $\psi = \psi(r,\theta)$.
In their spin orientation the conformal components are $\widehat A_{r\phi}=\widehat A_{\phi r}=3S\sin^2\theta/r^2$ and $\widehat A_{ij}\widehat A^{ij}=18S^2\sin^2\theta/r^6$.
The equivalent cylindrical coordinates $\rho=r\sin\theta$, $z=r\cos\theta$ have $d\ell^2=\psi^4(d\rho^2+dz^2+\rho^2d\phi^2)$, with $\psi=\psi(\rho,z)$.

## Intrinsic curvature

For a Euclidean conformal background the scalar curvature is ${}^{(3)}R=-8\psi^{-5}\nabla^2\psi$.
The vacuum maximal constraint therefore reads ${}^{(3)}R=K_{ij}K^{ij}=\psi^{-12}\widehat A_{ij}\widehat A^{ij}$.
For the single spinning hole this becomes ${}^{(3)}R=18S^2\sin^2\theta/(r^6\psi^{12})$.
The metric remains conformally flat while its intrinsic Ricci tensor is generally nonzero.

## The angle-reduced algebraic model

Altas and Tekin, arXiv:2007.14279, section II, equation (20), retain a source proportional to $(\boldsymbol{P}\cdot\boldsymbol{n})^2$ in the Hamiltonian equation for nonzero linear momentum.
Their equation (21) discards angular dependence to obtain a radial ordinary differential equation.
A radial conformal factor has a radial Euclidean Laplacian, while that full source depends on the polar angle for nonzero momentum.
The algebraic radial solution of the reduced equation consequently does not specify the general vacuum initial data in three dimensions.

## Numerical spinning throat

The drawings use one inversion throat with $P^i=0$, $S^i=(0,0,a^2)$ and $a$ as the unit of length.
The solver retains the entire angular source $18S^2\sin^2\theta/r^6$.
With $t=a/r$ and $\mu=\cos\theta$, its equation is

$$t^2\partial_t^2\psi+\partial_\mu[(1-\mu^2)\partial_\mu\psi]=-\frac94\left(\frac{S}{a^2}\right)^2t^4(1-\mu^2)\psi^{-7}.$$

The compact interval $0\le t\le1$ includes spatial infinity, where $\psi=1$.
The throat condition is $\partial_t\psi=\psi/2$ at $t=1$.
The correction to $1+t$ is collocated at 33 Chebyshev points in $t$, with eight even Legendre modes in $\mu$.
Gauss-Legendre projection resolves the nonlinear angular source at every fixed-point iteration.
The dedicated test samples the Hamiltonian equation between collocation points, checks the throat and infinity conditions, and compares against 41 radial points and ten angular modes.
The zero-spin solution is $\psi=1+a/r$.

The equatorial embedding reads the published spherical metric with this numerical conformal factor.
Its circumferential radius is $\mathcal R=r\psi(r,\pi/2)^2$, and its height obeys $Z'^2=\psi^4-\mathcal R'^2$.
The generator checks that this radicand is nonnegative throughout the drawn interval $a\le r\le8a$ and measures both drawn sheets against their intrinsic metric.
The lower sheet is the inversion image of the exterior, parametrized by its own exterior radial coordinate.
Both sheets meet at $r=a$ and $Z=0$.
The figure window encloses the whole surface with a square frame.

The curvature generator evaluates the stored intrinsic Ricci scalar in the supported meridian frame.
Its function declaration is a finite polynomial approximation in $a/r$ and $\cos\theta$, written through the existing closed-expression function representation.
This declaration records the numerical solution; it is not an exact algebraic solution of the constraint.
The map covers the exterior $r\ge a$ and marks the inversion throat.
There is no time coordinate, movie or stack.

## Absent causal diagrams

A three-dimensional Riemannian initial-data slice has no null directions or intrinsic causal compactification.
The null and conformal generators therefore draw no causal diagram for any of its three spatial charts.
The extrinsic data and constraints do not specify an explicit four-dimensional spacetime development.
