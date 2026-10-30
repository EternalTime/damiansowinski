# Gowdy universes with scalar and electromagnetic waves

## Sources and scope

The entry covers the minimally coupled scalar-wave and polarised electromagnetic families associated with the 1979 founding works.
Both founding bibliographic records and publisher abstracts were fetched and verified.
Their original full texts remain unread after access refusals.
The equations are sourced to readable later primary treatments, under the instruction to preserve that debt explicitly.
The history, proposed references, reciprocal Gowdy wording and literature ledger are preserved in the task data directory.

The scalar source is Szybka and Wyrębowski, Physical Review D 94, 024059 (2016), arXiv:1509.00346v3, equations (1)-(7) and (21)-(24), and sections 6.1, 6.3 and 6.4.
The electromagnetic source is Nungesser and Rendall, Classical and Quantum Gravity 26, 105019 (2009), arXiv:0901.4000v2, equations (2)-(16), with the global constraint on page 7.
The 1980 combined-field work supports a historical attribution from its verified abstract and supplies no formula to this addition.
Conformally coupled scalar extensions and the combined-field family are outside the proposed two-matter-system scope.

## Reference inventory

`charach_malin1979` is Cosmological model with gravitational and scalar waves, Physical Review D 19(4), 1058-1064, published 15 February 1979, DOI 10.1103/PhysRevD.19.1058.
`charach1979` is Electromagnetic Gowdy universe, Physical Review D 19(12), 3516-3523, published 15 June 1979, DOI 10.1103/PhysRevD.19.3516.
`charach_malin1980` is Cosmological model with gravitational, electromagnetic, and scalar waves, Physical Review D 21(12), 3284-3294, published 15 June 1980, DOI 10.1103/PhysRevD.21.3284.
`szybka_wyrebowski2016` is Backreaction for Einstein-Rosen waves coupled to a massless scalar field, Physical Review D 94(2), 024059, published 29 July 2016, DOI 10.1103/PhysRevD.94.024059.
`nungesser_rendall2009` is Strong cosmic censorship for solutions of the Einstein-Maxwell field equations with polarized Gowdy symmetry, Classical and Quantum Gravity 26(10), 105019, published 29 April 2009, DOI 10.1088/0264-9381/26/10/105019.
Each record was compared with its fetched publisher-deposited Crossref metadata, with the primary text or publisher abstract used for claim attribution.

## Coordinates and toroidal identifications

The areal coordinates are $(t,\theta,x,y)$ with $t>0$ and each spatial coordinate of period $2\pi$.
For both matter systems the areal line element is
$$ds^2=L^2\left\{t^{-1/2}e^{\lambda/2}(-dt^2+d\theta^2)+t(e^Pdx^2+e^{-P}dy^2)\right\}.$$
All coordinates are dimensionless, $L$ is a constant length, and $P$ and $\lambda$ depend on $t$ and $\theta$.
The scalar source writes the longitudinal coefficient as $e^f$, so $f=(\lambda-\ln t)/2$.

The logarithmic coordinates use $t=e^{-\tau}$ and the future runs towards decreasing $\tau$.
Their coefficients will be obtained by pulling back the areal metric, with $dt^2=e^{-2\tau}d\tau^2$.
The planned chart identifiers are `scalar_areal`, `scalar_logarithmic`, `electromagnetic_areal` and `electromagnetic_logarithmic`.
The electromagnetic logarithmic representation occurs in Nungesser and Lim (2013), but equations (1) and (3) there have typographical problems confirmed in the fetched TeX source.
The Lorentzian signs and the placement of the spatial differentials will therefore come from the 2009 areal equation (2) and the explicit pullback.

Periodicity of the longitudinal coefficient requires the integral of its spatial constraint around the circle to vanish.
This condition must be checked alongside the integer spatial wave numbers used in the drawings.

## Matter fields and equations

For the scalar system define the dimensionless field $\Phi=\sqrt{8\pi}\,\varphi$, where $\varphi$ is the scalar in the 2016 cosmological equations and $G=c=1$.
The matter equation is $R_{ab}=\partial_a\Phi\,\partial_b\Phi$.
Both $P$ and $\Phi$ obey $u_{tt}+t^{-1}u_t-u_{\theta\theta}=0$.
The longitudinal constraints in the common notation are
$$\lambda_t=t(P_t^2+P_\theta^2)+2t(\Phi_t^2+\Phi_\theta^2),\qquad
\lambda_\theta=2tP_tP_\theta+4t\Phi_t\Phi_\theta.$$
The chart printer checks these rescaled equations against every Ricci component at three independent rational derivative data in each chart.
The scalar field and its first and spatial second derivatives vary independently of the gravitational polarisation.

For the electromagnetic system use the potential $A=L\chi(t,\theta)\,dy$, with $F=dA$ and
$$R_{ab}=2\left(F_{ac}F_b{}^c-\tfrac14g_{ab}F_{cd}F^{cd}\right).$$
Nungesser and Rendall's equations (3), (5)-(7) with $\omega=0$ become
$$P_{tt}+t^{-1}P_t-P_{\theta\theta}=\frac{2e^P}{t}(\chi_t^2-\chi_\theta^2),\qquad
\chi_{tt}-\chi_{\theta\theta}=-P_t\chi_t+P_\theta\chi_\theta,$$
$$\lambda_t=t(P_t^2+P_\theta^2)+4e^P(\chi_t^2+\chi_\theta^2),\qquad
\lambda_\theta=2tP_tP_\theta+8e^P\chi_t\chi_\theta.$$
The factor $L$ in the potential accompanies the constant rescaling of the metric.
The chart printer checks every Einstein-Maxwell component, all source-free Maxwell equations and the Ricci trace at three independent rational derivative data in both electromagnetic charts.
The constant length varies between these data, which tests the factor of $L$ in the potential.

## Charts and curvature

All four charts publish the off-shell tensors generated from the line element by `print_charts.py`.
Keeping $P$ and $\lambda$ free makes every component independently readable and checkable through the existing metric reader.
The scalar and electromagnetic equations appear beside their own chart's parameters, with the logarithmic equations obtained by $\partial_\tau=-t\partial_t$.
No vacuum Ricci or Weyl components are substituted into either matter chart.
The dedicated test checks both logarithmic metrics against the areal metric pulled back by $\mathrm{diag}(-e^{-\tau},1,1,1)$ at three times.

## Drawings

The drawings keep the torus identifications, with integer spatial wave number one and zero integrated longitudinal constraint.
The plane of $t$ and $\theta$ is conformal to a flat cylinder, with null curves $t\pm\theta$ constant before the periodic identifications.
The Killing momenta along $x$ and $y$ vanish on these light planes and remain zero along their null geodesics.
The logarithmic ray diagrams plot $-\tau$ upwards, in the future direction.
Each conformal view uses $p,q=\arctan(t\mp\theta)$ and identifies the two spatial edges at the same $t$.

Define the Bessel primitive
$$B(t,\theta)=\tfrac12t^2(J_0(t)^2+J_1(t)^2)-tJ_0(t)J_1(t)\cos^2\theta.$$
The scalar drawings use $P=\tfrac12J_0(t)\cos\theta$, $\Phi=\tfrac7{20}J_0(t)\cos\theta$ and $\lambda=\tfrac{99}{200}B$.
The coefficient $99/200=(1/2)^2+2(7/20)^2$ keeps both scalar and gravitational contributions to the sourced constraints.
The electromagnetic drawings use the geodesic wave map
$$U=\tfrac12\ln t+\tfrac25J_0(t)\cos\theta,\qquad \bar P=\ln\cosh U,\qquad \chi=\tanh U.$$
Here $d\bar P^2+e^{2\bar P}d\chi^2=dU^2$, so the nonlinear wave-map equations reduce to the linear equation for $U$.
The longitudinal primitive is $\bar\lambda=\tfrac14\ln t+\tfrac25J_0(t)\cos\theta+\tfrac4{25}B$.
Equations (9)-(10) of the 2009 source then give $P=2\bar P+\ln t$ and $\lambda=4(\bar\lambda+\bar P)+\ln t$.
These are derived specialisations of the readable sources, with no claim that these particular parameters occur in the unavailable founding papers.
The dedicated tests check both local wave equations and both longitudinal constraints at five times, periodicity of every field and vanishing integrated momentum around the spatial circle.

Each matter system has one static spatial embedding at $L=1$, $t=1$ and $y=0$.
The radius is $\rho=L\sqrt{t}e^{P/2}$ and the height obeys $z_\theta^2=L^2t^{-1/2}e^{\lambda/2}-\rho_\theta^2$.
Both selected slices have positive embedding radicand around the complete period $0\leq\theta\leq2\pi$.
They are cut-open tubes with their end circles identified; their nonzero total height change precludes interpreting the displayed tube as a closed surface in Euclidean space.
The surface checks use the stored spatial components and adaptive quadrature of the height.
The square display windows enlarge the narrower dimension and retain the complete geometry through rotation.
There is one moment per view, so neither a movie nor a stack registration applies.

The scalar profiles have finite $P$, $\Phi$ and $\lambda$ at $t=0$, and the leading metric has Kasner exponents $(-1/3,2/3,2/3)$.
For the electromagnetic profile $U\sim\tfrac12\ln t+\tfrac25\cos\theta$, $P$ and $\lambda$ have finite limits there and the same exponents follow.
Both selected profiles therefore have a past curvature singularity.
The conformal boundary at future time infinity is one point after compactification, with winding rays continued across the identified edges.

The curvature file contains scalar-family $R$ and $K$ maps and an electromagnetic-family $K$ map over $0\leq\theta\leq2\pi$, $0.2\leq t\leq5$.
`curvature.py` reads the stored scalar expressions and substitutes exactly the declared profiles; it computes no curvature tensor.
Each map uses the existing rectangular frame, with signed logarithmic colour values.
The electromagnetic Ricci scalar vanishes by its field equations and has no separate map.

## Validation

The page-scoped check runs after the page is complete.
The History and reciprocal relation retain the exact approved wording.
The verified founding-paper metadata and abstracts remain the sources for historical attribution, with their refused full texts still recorded as owed in the literature ledger.
