# The Tolman-Oppenheimer-Volkoff star

The Tolman-Oppenheimer-Volkoff star takes the line element of its two founding papers of 1939, in the areal chart.

The components were computed and written by `_tools/derivations/print_charts.py`, which reads every value back through the checker's own reader before writing it.
`_tools/derivations/verify_metrics.py --system tov/spherical` then recomputes all of it from the line element, in about three seconds.

---

## Step 1. Which line element, and why

Tolman's paper of 1939 and the Oppenheimer and Volkoff paper that follows it in the same issue of the Physical Review both write the static sphere of fluid as

$$ds^2 = -e^{\lambda}dr^2 - r^2d\theta^2 - r^2\sin^2\theta\,d\phi^2 + e^{\nu}dt^2,$$

in the signature $(+,-,-,-)$ and in units with $G = c = 1$.
Oppenheimer and Volkoff then write the radial function through a mass, $e^{-\lambda} = 1 - 2u/r$, and it is in $u$ and $\nu$ that the equation now carrying their names is stated.

In the signature $(-,+,+,+)$ it reads, with $e^{\nu}$ written $e^{2\Phi}$ and $u$ written $m$:

$$ds^2 = -e^{2\Phi}c^2dt^2 + \frac{dr^2}{1 - 2m/r} + r^2\left(d\theta^2 + \sin^2\theta\,d\phi^2\right).$$

This is also the form of Misner, Thorne and Wheeler's chapter on stellar structure, so a reader coming from either source recognises it.
$m(r)$ is a length, $m = GM(r)/c^2$, as Schwarzschild's $r_s = 2GM/c^2$ folds a mass into a length, and $\Phi$ is dimensionless.
Both are left as free functions of $r$, because the equation of state that would fix them is exactly what the history says the equation cannot supply by itself.

It is the Morris-Thorne line element with $b = 2m$, and the two are meant to be read against each other: the same geometry, with a regular centre $m(0) = 0$ in one and a throat $b(b_0) = b_0$ in the other.

## Step 2. What the curvature says

With $T^\mu{}_\nu = \operatorname{diag}(-\rho c^2, p, p, p)$ and $G^\mu{}_\nu = 8\pi Gc^{-4}T^\mu{}_\nu$, the Einstein tensor gives three equations.

$$G^t{}_t = -\frac{2\partial_r m}{r^2} \quad\Longrightarrow\quad \partial_r m = \frac{4\pi G}{c^2}r^2\rho,$$

$$G^r{}_r = \frac{2r(r - 2m)\partial_r\Phi - 2m}{r^3} \quad\Longrightarrow\quad \partial_r\Phi = \frac{m + 4\pi Gr^3p/c^4}{r(r - 2m)},$$

and the angular equation $G^\theta{}_\theta = G^r{}_r$, which says the pressure is the same in every direction.

## Step 3. Why the angular equation is the hydrostatic one

For this metric the contracted Bianchi identity $\nabla_\mu G^\mu{}_r = 0$ reads

$$\partial_r G^r{}_r + \left(G^r{}_r - G^t{}_t\right)\partial_r\Phi + \frac{2}{r}\left(G^r{}_r - G^\theta{}_\theta\right) = 0$$

identically, for every $\Phi$ and every $m$; sympy reduces the left hand side to zero.
Read through the field equations it is $\partial_r p_r + (\rho c^2 + p_r)\partial_r\Phi + \frac{2}{r}(p_r - p_\perp) = 0$, so once $\rho$ and $p$ are defined by the first two equations, isotropy $p_\perp = p_r$ is exactly the statement $\partial_r p = -(\rho c^2 + p)\partial_r\Phi$.
Substituting Step 2's $\partial_r\Phi$ gives the Tolman-Oppenheimer-Volkoff equation

$$\partial_r p = -\frac{\left(\rho c^2 + p\right)\left(m + 4\pi Gr^3p/c^4\right)}{r(r - 2m)},$$

as the conventions state it.
With $p \ll \rho c^2$, $4\pi r^3p \ll Mc^2$ and $m \ll r$ it becomes $\partial_r p = -GM(r)\rho/r^2$, Newton's hydrostatic equation.

## Step 4. A check against the interior Schwarzschild solution

Put the uniform density star into the components: $m = r_s r^3/(2R^3)$ and $e^{\Phi} = \frac{1}{2}\left(3\sqrt{1 - r_s/R} - \sqrt{1 - r_s r^2/R^3}\right)$, which is the interior Schwarzschild line element.
Then $G^r{}_r - G^\theta{}_\theta$ vanishes identically and $G^t{}_t = -3r_s/R^3$, the value of the interior Schwarzschild solution.
So the general chart contains the interior Schwarzschild solution as its member with constant density, and agrees with it.

## Step 5. How the values are written

Every Riemann component of this chart carries the combination $(\partial_r\Phi)^2 + \partial_r^2\Phi$, which is $e^{-\Phi}\partial_r^2 e^{\Phi}$, the tidal stretching along the radius.
`chart_printer.py` collects each numerator around that combination and around $\partial_r\Phi$, with the coefficients factored, so the components read as Morris-Thorne's hand written ones do.

The Kretschmann scalar is written as the sum of squares of the four independent components in the orthonormal frame of a static observer,

$$K = 4A^2 + 8B^2 + 8C^2 + 4D^2,$$

with $A = R_{\hat t\hat r\hat t\hat r}$, $B = R_{\hat t\hat\theta\hat t\hat\theta}$, $C = R_{\hat r\hat\theta\hat r\hat\theta}$ and $D = R_{\hat\theta\hat\phi\hat\theta\hat\phi}$, the factors counting how often each appears among the 256 components.
`print_charts.py` checks that form against sympy's contraction before writing it.
Outside the star, $m = GM/c^2$ and $e^{2\Phi} = 1 - 2m/r$, and it reduces to Schwarzschild's $48m^2/r^6$.

## Step 6. Dimensions

`DIMENSIONS` declares $t$ a time, $r$ a length, $\Phi$ a pure number and $m$ a length.
Every term of every component balances on that declaration.
