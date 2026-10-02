# The static sphere of Brans and Dicke's theory

The static, spherically symmetric, asymptotically flat vacuum of Brans and Dicke's theory of gravity, Brans's class I.
This note records the three charts, where each comes from, the field equations each is held to, and what each diagram draws.
Units here have $G = c = 1$ where a source does; the metric file keeps $c$ and writes every mass as a length.

## Step 1: Brans and Dicke's isotropic chart

Brans and Dicke's (31) and (32), with their constants $\alpha_0 = \beta_0 = 0$: $ds^2 = -h^{2/\lambda}dt^2 + (1 + B/\rho)^4h^{2(\lambda - C - 1)/\lambda}(d\rho^2 + \rho^2d\Omega^2)$ with $h = (\rho - B)/(\rho + B)$, and their (33), $\lambda^2 = (C + 1)^2 - C(1 - \omega C/2)$.
They write $r$ for the isotropic radius; the metric file writes $\rho$, as Fisher's page does, and keeps $r$ for the radius of Step 3.
Their footnote 22 says the form of the solution was suggested to Brans by C. Misner.
The paper prints the scalar field as $\phi = \phi_0h^{-C/\lambda}$, and the field equations hold with the other sign, $\phi = \phi_0h^{C/\lambda}$, which is the sign of Brans's Scholarpedia article, his (35), and the one that agrees with their own weak field (28), where $\phi$ grows toward the mass and $C = -1/(2 + \omega)$ is negative.
`print_charts.brans_dicke_check` holds the chart to the equations with that sign, so the sign in (32) is a misprint.
The metric file takes $B$, $C$, and $\lambda$ as free parameters and reads $\omega = 2(\lambda^2 - C^2 - C - 1)/C^2$ off them, which is (33) solved for $\omega$.
Far from the body $g_{tt} = -(1 - 4B/\lambda\rho)$ and $g_{\rho\rho} = 1 + 4B(C + 1)/\lambda\rho$, so $GM/c^2 = 2B/\lambda$ and the post-Newtonian parameter is $\gamma = C + 1$.
Their (34) for a body whose gravity is weak throughout is $C = -1/(2 + \omega)$, which gives $\lambda = \sqrt{(2\omega + 3)/(2\omega + 4)}$ and $\gamma = (1 + \omega)/(2 + \omega)$.

## Step 2: the field equations

With no matter, Brans and Dicke's (11) and (13) are $R_{\mu\nu} = \omega\,\partial_\mu\phi\,\partial_\nu\phi/\phi^2 + \nabla_\mu\nabla_\nu\phi/\phi$ and $\nabla^\mu\nabla_\mu\phi = 0$, which are Campanelli and Lousto's (5) and (6).
`brans_dicke_check` holds every chart to both, the first in every slot with the Hessian written through the chart's own Christoffel symbols, and the second through logarithmic derivatives of the determinant so that no root is taken.
It also holds each chart to a vanishing Ricci tensor without the field.
The reader holds $h^{2/\lambda}$ as a power of $B - \rho$ with $(-1)^{2/\lambda}$ beside it, a phase of its generator, so the chart's powers are counted on that base and the pullback of Step 3 compares moduli at rational points to thirty digits.
The printer writes an exponent over $\lambda$ on one line, $h^{(2C + 2)/\lambda}$.

## Step 3: Campanelli and Lousto's chart

Campanelli and Lousto's (7), "a power generalization of the Schwarzschild metric": $ds^2 = -A^{m+1}dt^2 + A^{n-1}dr^2 + r^2A^{n}d\Omega^2$ with $A = 1 - 2r_0/r$ and $\phi = \phi_0A^{-(m+n)/2}$.
Agnese and La Camera, and Vanzo, Zerbini, and Faraoni, use the same chart.
The map from Step 1 is $r = \rho(1 + B/\rho)^2$, for which $A = h^2$ and $(1 + B/\rho)^4(d\rho^2 + \rho^2d\Omega^2) = dr^2/A + r^2d\Omega^2$, so $r_0 = 2B$, $m = 1/\lambda - 1$, and $n = 1 - (C + 1)/\lambda$.
In these letters $\omega = -2(m^2 + n^2 + mn + m - n)/(m + n)^2$, $GM/c^2 = (m + 1)r_0$, and $\gamma = (1 - n)/(1 + m)$.
The Kretschmann scalar diverges on $r = 2r_0$ as $(r - 2r_0)^{-2 - 2n}$ times a polynomial, and the spheres there have the area $4\pi r^2A^{n}$, which vanishes for $n > 0$.

## Step 4: Dicke's units and Fisher's metric

Dicke's transformation of units is the conformal map $\tilde g_{\mu\nu} = (\phi/\phi_0)g_{\mu\nu}$, under which the theory is Einstein's with a massless scalar field.
In the chart of Step 3, $(\phi/\phi_0)g_{\mu\nu}$ has $-A^{1 + (m - n)/2}$, $A^{-1 - (m - n)/2}$, and $r^2A^{(n - m)/2}$, which is the spherical chart of `fisher_jnw` with $b = 2r_0$ and $\gamma = 1 + (m - n)/2 = (C + 2)/2\lambda$.
`brans_dicke_check` compares the two component by component at rational points, against the chart `print_charts.fisher_jnw` writes.
Faraoni, Hammad, Cardini, and Gobeil state the same of the general solution: it is conformal to the Fisher-Wyman geometry.

## Step 5: Bronnikov's harmonic coordinate

Bronnikov, Constantinidis, Evangelista, and Fabris's (9), from Bronnikov's paper of 1973: $ds^2 = \phi^{-1}\left(e^{-2bu}dt^2 - e^{2bu}s^{-2}(k, u)\left(du^2s^{-2}(k, u) + d\Omega^2\right)\right)$ with $s(k, u) = k^{-1}\sinh ku$ for $k > 0$, and their (13), $e^{-2ku} = 1 - 2k/r$.
The bracket is Fisher's metric in its harmonic chart with the mass $b$, and their (10) makes $\ln\phi$ linear in $u$.
The metric file writes $\phi = \phi_0e^{-su}$ with $s$ a length, so that $ds^2 = e^{su}\left(-e^{-2bu}c^2dt^2 + k^2e^{2bu}\sinh^{-2}(ku)\left(k^2du^2\sinh^{-2}(ku) + d\Omega^2\right)\right)$.
Then $k = r_0 = 2B$, $b = (C + 2)B/\lambda$, $s = 2BC/\lambda$, and their (12), $2k^2 = 2b^2 + S$ with $S = (\omega + 3/2)s^2$, is $\omega = 2(k^2 - b^2)/s^2 - 3/2$.
The coordinate $u$ is harmonic for Fisher's metric, the one in Dicke's units, and is an inverse length, which `DIMENSIONS` declares as `L**(-1)`.
`brans_hyperbolic` prints the chart in $\sinh(ku)$, $\cosh(ku)$, and one exponential.

## Step 6: the diagrams

Every diagram is drawn at $\omega = 6$, the least value Brans and Dicke's bound from Mercury allows, and at $C = -1/4$, where $\lambda = 1$ and every power is rational: $m = 0$, $n = 1/4$, $b = 7k/8$, $s = -k/4$.
That is twice the scalar charge of a body of weak gravity, whose $C = -1/8$ gives $n = 0.096$ and confines the difference from Schwarzschild's metric to within a few thousandths of $2r_0$: its level circle of the embedding is $r = 2.005\,r_0$, and light from $r = 4r_0$ needs $65\,r_0/c$ to reach the singularity.
At $C = -1/4$ the tortoise coordinate has $dr_*/dr = A^{-7/8}$ and $r_* = 16r_0A^{1/8}\,{}_2F_1(2, 1/8; 9/8; A)$, which vanishes at $r = 2r_0$; `null_rays.py --verify` checks the rays of all three charts against $ct \mp r_*$.

Since $r_*$ is finite at $r = 2r_0$, $p, q = \arctan((ct \mp r_*)/\ell)$ bring the spacetime into Minkowski's triangle with the singularity a timelike line on $X = 0$, the diagram of Fisher's metric, to which this one is conformal on every plane of $t$ and the radius.

The equator of a moment of $t$ has circles of radius $rA^{n/2}$ and is Fisher's slice at $\gamma = 1 - n$, since the two metrics differ by the factor $\phi/\phi_0$ alone, which the slice of Fisher's at $\gamma = (C + 2)/2\lambda$ does not share: the two pages draw different surfaces.
The surface lies level on $r_e = (2 - n)^2r_0/2(1 - n) = 49r_0/24$, in flat space outside it and in three dimensional Minkowski space inside.
With $w = A^{1/8} = e^{-u/4}$ the radius is $2w/(1 - w^8)$, the level circle is $w^8 = 1/49$, and the height is the integral of $2\sqrt{|49w^8 - 1|}/(1 - w^8)^{3/2}\,dw$ from $7^{-1/4}$, which `embedding.py` checks by a quadrature of its own.
It enters the singular point along the light cone, with $1 - (dZ/d\rho)^2 = 64A$, while the radius falls only as $A^{1/8}$, so the difference of $d\rho$ and $dZ$ is lost to rounding before the circles are small.
The profile is written to sixteen decimals and stops at $u = 12$, a circle of radius $2e^{-3}r_0 = 0.1\,r_0$.
The geometry is static, so it is one surface and no movie.
