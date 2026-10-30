# Gowdy universes with scalar and electromagnetic waves

## Sources and scope

The entry covers the minimally coupled scalar-wave and polarised electromagnetic families associated with the 1979 founding works.
Both founding bibliographic records and publisher abstracts were fetched and verified.
Their original full texts remain unread after access refusals.
The equations are sourced to readable later primary treatments, under the instruction to preserve that debt explicitly.
The history, proposed references, related wording and literature ledger await review in the task data directory.

The scalar source is Szybka and Wyrębowski, Physical Review D 94, 024059 (2016), arXiv:1509.00346v3, equations (1)-(7) and (21)-(24), and sections 6.1, 6.3 and 6.4.
The electromagnetic source is Nungesser and Rendall, Classical and Quantum Gravity 26, 105019 (2009), arXiv:0901.4000v2, equations (2)-(16), with the global constraint on page 7.
The 1980 combined-field work supports a historical attribution from its verified abstract and supplies no formula to this addition.
Conformally coupled scalar extensions and the combined-field family are outside the proposed two-matter-system scope.

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
These rescaled equations will receive an independent tensor check after history approval.

For the electromagnetic system use the potential $A=L\chi(t,\theta)\,dy$, with $F=dA$ and
$$R_{ab}=2\left(F_{ac}F_b{}^c-\tfrac14g_{ab}F_{cd}F^{cd}\right).$$
Nungesser and Rendall's equations (3), (5)-(7) with $\omega=0$ become
$$P_{tt}+t^{-1}P_t-P_{\theta\theta}=\frac{2e^P}{t}(\chi_t^2-\chi_\theta^2),\qquad
\chi_{tt}-\chi_{\theta\theta}=-P_t\chi_t+P_\theta\chi_\theta,$$
$$\lambda_t=t(P_t^2+P_\theta^2)+4e^P(\chi_t^2+\chi_\theta^2),\qquad
\lambda_\theta=2tP_tP_\theta+8e^P\chi_t\chi_\theta.$$
The factor $L$ in the potential accompanies the constant rescaling of the metric.
The Maxwell equation, trace-free Ricci tensor and these constraints will be checked for both electromagnetic charts.

## Charts and curvature

Tensor generation and independent field-equation checks remain to be done after review.
The scalar normalization, electromagnetic normalization and pullbacks must be checked explicitly before the page is published.

## Drawings

The drawings will keep the torus identifications, with integer spatial wave number and zero integrated longitudinal constraint.
The plane of $t$ and $\theta$ is conformal to a flat cylinder, with null curves $t\pm\theta$ constant before the periodic identifications.
Parameters, embedding domains and curvature scalars have not yet been checked or drawn.

## Validation

The page-scoped check will run after the page is complete.
