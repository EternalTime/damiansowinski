# Topological black holes in anti-de Sitter space

The static vacuum solutions with $\Lambda = -3/L^2$ over a surface of constant curvature $k$,

$$ds^2 = -f\,c^2dt^2 + \frac{dr^2}{f} + r^2\,d\Sigma_k^2,\qquad f = k - \frac{\mu}{r} + \frac{r^2}{L^2}.$$

Its six charts are written by `_tools/derivations/print_charts.py --metric topological_black_hole`, and `verify_metrics.py --system topological_black_hole/<chart>` checks each in seconds.

## Step 1. The parameters

The parameters are anti-de Sitter's radius $L$, as in the anti-de Sitter and Schwarzschild-anti-de Sitter entries, and a length $\mu$.
At $k = 1$ the line element is Schwarzschild-anti-de Sitter's with $\mu = r_s$.
Brill, Louko and Peldán write $F = k - 2M/R + Q^2/R^2 - \Lambda R^2/3$, their (2.1), so $\mu = 2M$ at $Q = 0$, and a closed horizon whose $d\Sigma_k^2$ has area $V$ has the mass $(V/4\pi)M$.
Birmingham writes $f = k - \omega_dM/r^{d-3} + r^2/l^2$ in $d$ dimensions, his (3) and (4), and Vanzo $V = -1 - 2\eta/r + r^2/l^2$, his (1.4).

## Step 2. The charts and their sources

- `static`: the family as Brill, Louko and Peldán and Birmingham write it, with $d\Sigma_k^2 = d\rho^2/(1 - k\rho^2) + \rho^2d\phi^2$, one form for all three curvatures, so that $k$ is a parameter the checker keeps free. $\rho = \sin\theta$, $\theta$ and $\sinh\theta$ carry it to the three forms of their (2.2).
- `eddington_finkelstein_ingoing` and `eddington_finkelstein_outgoing`: $v, u = ct \pm r_*$ with $dr_*/dr = 1/f$ and $r_*(0) = 0$, as the Schwarzschild-anti-de Sitter entry builds them.
- `black_string`: Lemos's (2), Phys. Lett. B 353, 46, $ds^2 = -(\alpha^2r^2 - b/(\alpha r))dt^2 + dr^2/(\alpha^2r^2 - b/(\alpha r)) + r^2d\phi^2 + \alpha^2r^2dz^2$, with $\alpha = 1/L$ and $b = \mu/L$. His horizon $\alpha r = b^{1/3}$ is $r_h = (\mu L^2)^{1/3}$, and his Kretschmann scalar $24\alpha^4(1 + b^2/(2\alpha^6r^6))$ is $24/L^4 + 12\mu^2/r^6$.
- `brane`: Hartnoll's (42) and (43), Class. Quantum Grav. 26, 224002, $ds^2 = (L^2/r^2)(-f\,dt^2 + dr^2/f + dx^idx^i)$ with $f = 1 - (r/r_+)^d$, at $d = 3$ with his $r$ written $z$. It is the static chart at $k = 0$ under $z = L^2/r$, $x = L\rho\cos\phi$, $y = L\rho\sin\phi$, with $z_h = L^2/r_h$ and $\mu = L^4/z_h^3$.
- `hyperbolic`: Mann's (6), (7) and (10), Class. Quantum Grav. 14, L109, at $q = 0$, and his (2) of Class. Quantum Grav. 14, 2927, with $d\Omega^2 = d\theta^2 + \sinh^2\theta\,d\phi^2$. It is the static chart at $k = -1$ under $\rho = \sinh\theta$.

`topological_check` in `print_charts.py` holds each chart to $R_{\mu\nu} = -(3/L^2)g_{\mu\nu}$ and each after the first to being the static chart pulled back through these maps.

## Step 3. The horizons

$L^2rf = r^3 + kL^2r - L^2\mu$.
For $k = 0$ the one real root is $r_h = (\mu L^2)^{1/3}$.
For $k = -1$ the cubic $r^3 - L^2r - L^2\mu$ has one positive root for $\mu \ge 0$, which is $L$ at $\mu = 0$, and two for $-2L/(3\sqrt{3}) < \mu < 0$, which meet at $r = L/\sqrt{3}$ at the lower bound, Vanzo's and Brill, Louko and Peldán's $M_{\rm crit} = -l/(3\sqrt{3})$.
The temperature is $T = \hbar c\,f'(r_h)/(4\pi k_B) = \hbar c(3r_h^2 + kL^2)/(4\pi k_BL^2r_h)$, Brill, Louko and Peldán's (3.4).

## Step 4. What the diagrams draw

Three spacetimes, in units of $L$.

- The flat hole at $\mu = L$: $r_h = L$, $\kappa = 3/(2L)$, and $r_* = \frac{1}{3}\ln|1 - r| - \frac{1}{6}\ln(r^2 + r + 1) + (\arctan((2r + 1)/\sqrt{3}) - \pi/6)/\sqrt{3}$, which tends to $R = \pi/(3\sqrt{3})$. Kruskal's $UV = (1 - r)(r^2 + r + 1)^{-1/2}\exp(\sqrt{3}(\arctan((2r + 1)/\sqrt{3}) - \pi/6))$ is $-e^{\pi/\sqrt{3}} = -6.134$ on the boundary, and the ray that leaves a boundary at $t = 0$ meets $r = 0$ at $X = 2\arctan e^{\pi/(2\sqrt{3})} - \pi/2 = 0.803$.
- The hyperbolic hole without mass: $f = r^2/L^2 - 1$, the BTZ hole's $N^2$ at $M = 1$ and $J = 0$, so its plane of $t$ and $r$ and its conformal diagram are that hole's.
- The hyperbolic hole at $\mu = -120L/343$: with roots $a$, $b$ and $-(a + b)$ the cubic needs $a^2 + ab + b^2 = L^2$ and $\mu = -ab(a + b)/L^2$, and $a = 5L/7$, $b = 3L/7$ solve the first in rationals. Then $f'(r_i) = 26/35$, $-22/21$ and $-143/56$, $r_* = \frac{35}{26}\ln|1 - 7r/5| - \frac{21}{22}\ln|1 - 7r/3| - \frac{56}{143}\ln(1 + 7r/8)$, which vanishes at $r = 0$ and tends to $R = -0.3035$, and $\kappa_+ = 13/35$.

## Step 5. The embedding

On the plane $z = 0$ of the black string at $t = 0$ the metric is $dr^2/f + r^2d\phi^2$ with $f = r^2/L^2 - \mu/r$.
In flat space $dz/dr = \sqrt{1/f - 1}$, real from the throat $r_h$ to the root of $r^3 - L^2r - \mu L^2$, which at $\mu = L$ is $1.3247\,L$, and beyond it the slice is drawn in Minkowski space at $dZ/dr = \sqrt{1 - 1/f}$, as Schwarzschild-anti-de Sitter's is.
The horizon of the hyperbolic hole at $\mu = 0$ has the metric $L^2(d\theta^2 + \sinh^2\theta\,d\phi^2)$, whose circles grow as $\cosh\theta$, so it is the sheet $Z = L\cosh\theta - L$ of a hyperboloid in Minkowski space.
