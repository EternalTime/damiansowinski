# The Bell-Szekeres colliding electromagnetic waves

Where both waves have passed, $u \ge 0$, $v \ge 0$ and $au + bv < \pi/2$, Bell and Szekeres's metric is

$$ds^2 = -2\,du\,dv + \cos^2(au - bv)\,dx^2 + \cos^2(au + bv)\,dy^2 .$$

All six charts are written by `_tools/derivations/print_charts.py --metric bell_szekeres`, and `verify_metrics.py --system bell_szekeres/<chart>` checks each in seconds.
Every source writes the signature $(+,-,-,-)$, and every chart here is its source's with the signature flipped.

## Step 1. Sources

The double null chart is eq. (15.7) of J. B. Griffiths, *Colliding Plane Waves in General Relativity* (Clarendon Press, 1991), after P. Bell and P. Szekeres, Gen. Rel. Grav. 5, 275 (1974), and eq. (4) of A. Feinstein and M. A. Pérez Sebastián, Class. Quantum Grav. 12, 2723 (1995).
Ahead of one wave or both the metric is the same with $u$ or $v$ set to zero, Griffiths's eq. (15.3).
The chart of $\xi = au + bv$ and $\eta = bv - au$ is Feinstein and Pérez Sebastián's eq. (9).
The regular chart is the one Bell and Szekeres removed the coordinate singularity with, Griffiths's eqs. (15.11) and (15.12).
The global chart is C. J. S. Clarke and S. A. Hayward's, Class. Quantum Grav. 6, 615 (1989), as Griffiths gives it in eqs. (15.13) to (15.15).
The Kruskal-Szekeres chart is Feinstein and Pérez Sebastián's eqs. (10) to (13), after M. Dorca and E. Verdaguer, Nucl. Phys. B 403, 770 (1993).
The Bertotti-Robinson chart is the conformally flat form of Griffiths's eq. (15.9) and Feinstein and Pérez Sebastián's eq. (8).

## Step 2. Units

$u$, $v$, $x$ and $y$ are lengths and $a$ and $b$ inverse lengths, so $au \pm bv$ is a pure number.
$q = 1/\sqrt{2ab}$ is the radius of both factors of AdS$_2 \times$ S$^2$, and $k = \sqrt{2ab} = 1/q$.
$\xi$, $\eta$ and the coordinates of the regular, global and Bertotti-Robinson charts are pure numbers; the Kruskal-Szekeres $U$ and $V$ are lengths.

## Step 3. The maps

With $\xi = au + bv$ and $\eta = bv - au$, $-2\,du\,dv = (-d\xi^2 + d\eta^2)/(2ab)$, so

$$ds^2 = \frac{-d\xi^2 + d\eta^2}{2ab} + \cos^2\eta\,dx^2 + \cos^2\xi\,dy^2 :$$

$(\xi, y)$ is an anti-de Sitter space of two dimensions and $(\eta, x)$ a sphere, with $\eta$ the latitude and $kx$ the longitude.
The regular chart is $T = -\cos\xi\cosh ky$, $Z = \cos\xi\sinh ky$, $X = \cos\eta\cos kx$, $Y = \cos\eta\sin kx$, so $T^2 - Z^2 = \cos^2\xi$ and $X^2 + Y^2 = \cos^2\eta$; it covers one sign of $\eta$, with $\eta = 0$ on its rim $X^2 + Y^2 = 1$.
The global chart is $\sinh\rho = Z$, $\sin\chi = T/\cosh\rho$, $\theta = \pi/2 - \eta$, $\phi = kx$, so $\cos\chi\cosh\rho = \sin\xi$: the collision $\xi = 0$ is $\chi = -\pi/2$ and the horizon is $\cos\chi\cosh\rho = 1$.
Griffiths's eq. (15.13) writes $\tilde\theta$ and $\tilde\phi$ through $\tilde X$ and $\tilde Y$; the ranges of his eq. (15.15), $|\cos\tilde\theta| < \cos\chi\cosh\rho$ and $-\infty < \tilde\phi < \infty$, are those of $\theta = \pi/2 + au - bv$ and $\phi = kx$, the angles of his eq. (15.8), which are the ones used here.
The Kruskal-Szekeres chart is $U = -q\,e^{ky}\cos\xi/(1 + \sin\xi)$, $V = -q\,e^{-ky}\cos\xi/(1 + \sin\xi)$, so $2abUV = (1 - \sin\xi)/(1 + \sin\xi)$ and $(1 + \sin\xi)^2 = 4/(1 + 2abUV)^2$.
The Bertotti-Robinson chart is $t = e^{ky}\tan\xi$, $r = e^{ky}\sec\xi$, $\theta = \pi/2 - \eta$, $\phi = kx$, which puts the region on $0 \le t < r$.
Griffiths's eq. (15.8) and Feinstein and Pérez Sebastián's eq. (7) reach the same line element through $t + r = \coth(\tfrac12\mathrm{sech}^{-1}\cos\xi - y/2q)$ and $t - r = -\tanh(\tfrac12\mathrm{sech}^{-1}\cos\xi + y/2q)$, which differs from the map here by an isometry of the anti-de Sitter factor and sends half of the region to $r < 0$ with time running down; the map here keeps $r > 0$ and $t$ increasing to the future.
`bell_szekeres_check` pulls every chart back through its map onto the double null chart at random points, to forty digits.

## Step 4. The field equations

In units with $G = c = 1$ and Gaussian fields, $G_{\mu\nu} = 2\left(F_{\mu\alpha}F_\nu{}^\alpha - \tfrac14 g_{\mu\nu}F_{\alpha\beta}F^{\alpha\beta}\right)$ with $F = dA$, $A = \sin(au + bv)\,dy$, which is $k$ times the area form of the anti-de Sitter factor, a uniform electric field.
`bell_szekeres_check` checks this in every slot of every chart, with $dF = 0$ and $\partial_\mu(\sqrt{-g}F^{\mu\nu}) = 0$, and that the Weyl tensor vanishes.
$G_{uu} = 2a^2$ and $G_{vv} = 2b^2$ are the two waves, $R = 0$, and $K = 32a^2b^2 = 8/q^4$ everywhere, the Bertotti-Robinson universe's $8/b^4$ at $b = q$.

## Step 5. The horizon

$g_{yy} = \cos^2\xi$ vanishes at $\xi = \pi/2$, where the Killing vector $\partial_y$ is null: a Killing-Cauchy horizon.
In the regular and Kruskal-Szekeres charts $\partial_y$ is a boost, $-k(Z\partial_T + T\partial_Z)$ and $k(U\partial_U - V\partial_V)$, the horizon is $T = -|Z|$ or $UV = 0$, and every event of it at finite $y$ is the one event $T = Z = 0$, $U = V = 0$ where its two null branches cross.
An observer at rest in $\eta$, $x$ and $y$ has no $\Gamma^\mu{}_{\xi\xi}$ to turn it, so it is a geodesic, and reaches the horizon after the proper time $q\pi/2$.

## Step 6. The diagrams

The spacetime diagrams draw the plane of each chart's first two coordinates at $a = b = 1$; `null_rays.py --verify` compares their rays with the closed forms: $u$ and $v$; $\xi \pm \eta$; $(-T \mp Z)/(1 + \sqrt{1 - T^2 + Z^2})$; $\chi \pm \arctan\sinh\rho$; $U$ and $V$; and $t \pm r$.
The conformal diagram draws the plane $x = y = 0$ with all four regions, as Khan and Penrose's, and the plane $\eta = 0$, $x = 0$ in the strip of the anti-de Sitter factor, $p, q = (\chi \mp \arctan\sinh\rho)/2$, where $\tan p = kU = t - r$ and $\tan q = kV = -1/(t + r)$.
The embedding diagram is the plane of $x$ and $y$ on $\eta = 0$, flat at each $\xi$ with the metric $dx^2 + \cos^2\xi\,dy^2$, with a ring of free particles at rest, the ellipse of semi-axes $\ell$ and $\ell\cos\xi$, stacked into its world tube.
