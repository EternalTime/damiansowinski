# Cylindrical waves around a string

## Coordinates and source

Gürtuğ and Sakallı, equations (1)-(3), give the vacuum metric in cylindrical coordinates and the auxiliary coordinates $\tau$ and $\sigma$.
Their equation (55) introduces $\tau=\sinh\psi$ and $\sigma=\cosh\theta$.
The two stored charts cover $\sigma>1$ and $\theta>0$ respectively, excluding the axis where polar coordinates fail.
The source uses $q^2-p^2=1$; the stored charts choose $p>0$ and $q=\sqrt{1+p^2}$, with the other signs obtained by reversing coordinate orientations.
The length $\ell$ restores an overall scale to the source's dimensionless coordinates.

The cylindrical time and radius satisfy $t=\tau\sigma$ and $\rho^2=(1+\tau^2)(\sigma^2-1)$.
For $t\in\mathbb R$ and $\rho>0$, put $D=\sqrt{(1+\rho^2-t^2)^2+4t^2}$, $\sigma=\sqrt{(D+1+\rho^2-t^2)/2}$ and $\tau=t/\sigma$.
The relation has a unique nonnegative solution for $\tau^2$, since the product of its two roots is $-t^2$.
The cylindrical line element is given in the History, and the explicit inverse above fixes it at every event in this domain.
The cylindrical tensors are not stored: the explicit radical chart did not finish the printer's 300 second chart step.

## Checks and drawings

The two stored charts were printed from the source expressions and checked against the source metric and the vacuum equations at two numerical points each.
The page-local tests take the circumference and radial distance from the stored auxiliary line element, at fixed $t=0$, and check both angle deficits and the regular-axis case $\alpha p=1$.

Both radial diagrams show the surface of fixed $\phi$ and $z$ at $p=1$, $\alpha=2$ and $\ell=1$.
The null directions have $\operatorname{arsinh}\tau\pm\operatorname{arcosh}\sigma$ constant, or $\psi\pm\theta$ constant.
The conformal views compactify $t\pm\rho$ and mark the string axis.
The embedding shows the surface $\psi=0$, $z=0$ for $0\leq\theta\leq2$, and its isometry is checked by the embedding generator.
The curvature map reads the stored Kretschmann scalar on $1\leq\sigma\leq4$ and $|\tau|\leq3/2$, at the same parameter values.
