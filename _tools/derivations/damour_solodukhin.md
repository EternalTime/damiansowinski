# The Damour-Solodukhin wormhole

Damour and Solodukhin's wormhole is Schwarzschild's metric with a constant $\lambda^2$ added to $-g_{tt}$.
This note records the five charts, where each comes from, what the curvature is at the throat, and what each diagram draws.
Every chart writes $r_s = 2GM/c^2$ for the areal radius of the throat, with $M$ Damour and Solodukhin's mass parameter, and $a = 1 + \lambda^2$ below.

## Step 1: Damour and Solodukhin's chart

Their equation (2.1) is $ds^2 = -(1 - r_s/r + \lambda^2)c^2dt^2 + dr^2/(1 - r_s/r) + r^2d\Omega^2$, with $r \ge r_s$ on each of two sides.
At $r = r_s$ the component $g_{rr}$ diverges and $g_{tt} = -\lambda^2$ does not vanish, so $\partial_t$ is timelike everywhere and $r = r_s$ is a throat, the smallest sphere, where the chart ends.
It is the Morris-Thorne metric with the constant shape function $b = r_s$ and $e^{2\Phi} = 1 - r_s/r + \lambda^2$.

## Step 2: the rescaled time

Far away $g_{tt} \to -a$, so $t$ is not the proper time of a clock at rest there.
Bueno, Cano, Goelen, Hertog and Vercnocke redefine $t \to t/\sqrt{a}$ and $M \to M(1 + \lambda^2)$ and write, as their (35) and (36), $ds^2 = -(1 - 2M/r)dt^2 + dr^2/(1 - 2M(1 + \lambda^2)/r) + r^2d\Omega^2$.
Their $M$ is $1/a$ of Damour and Solodukhin's, so their throat $r = 2M(1 + \lambda^2)$ is the same sphere $r = r_s$, and with $r_s$ kept as the throat's radius the chart reads $ds^2 = -(1 - r_s/ar)c^2dt^2 + dr^2/(1 - r_s/r) + r^2d\Omega^2$.
Bronnikov, Konoplya and Pappas write it in exactly that form, their (4.14), $f = 1 - 2m/r(1 + \lambda^2)$ with the throat at $r_0 = 2m$.
That keeps one meaning of $r_s$ on every chart and every drawing.
Tsukamoto's (A6) is the same statement about the two masses read from the other side: in the rescaled form the ADM mass is $M(1 + \lambda^2)$.
The mass a distant orbit measures is $r_sc^2/2Ga$, from $g_{tt}$, while $g_{rr}$ falls off as for the mass $r_sc^2/2G$; the two differ because the matter's pressures reach to infinity, falling as $\lambda^2r_s/r^3$.

## Step 3: through the throat

Bueno and his collaborators write the relation between $r$ and the tortoise coordinate of their time with an auxiliary variable $\rho$, their (40) and (41): $r/M = 2 + \lambda^2(1 + \cosh\rho)$ and $r_*/M = (2 + \lambda^2)\rho + \lambda^2\sinh\rho$, with $\rho$ over the whole line and the throat at $\rho = 0$.
With $r_s$ the throat's radius that is $r = r_s(2 + \lambda^2(1 + \cosh\rho))/2a$ and $x = r_s((2 + \lambda^2)\rho + \lambda^2\sinh\rho)/2a$.
Then $dr/d\rho = r_s\lambda^2\sinh\rho/2a$, $r - r_s/a = r_s\lambda^2(1 + \cosh\rho)/2a$ and $r - r_s = r_s\lambda^2(\cosh\rho - 1)/2a$, so $dr^2/(1 - r_s/r) = (1 - r_s/ar)\,r^2d\rho^2$ and $dx = r\,d\rho$.
Taking $\rho$ as the coordinate, the metric is $(1 - r_s/ar)(-c^2dt^2 + r^2d\rho^2) + r^2d\Omega^2$, regular at $\rho = 0$, where $1 - r_s/ar = \lambda^2/a$.

Bueno and his collaborators make their rotating wormhole regular at its throat with $r = r_+ + \rho^2/M$, whose static case is Einstein and Rosen's coordinate.
Einstein and Rosen removed the same divergence of $g_{rr}$ from Schwarzschild's metric in 1935 with $u^2 = r - 2m$, their (5a): $ds^2 = -4(u^2 + 2m)du^2 - (u^2 + 2m)^2d\Omega^2 + u^2/(u^2 + 2m)\,dt^2$ in their signature.
The spatial part of Damour and Solodukhin's metric is Schwarzschild's, so the same $u$ gives $ds^2 = -(u^2/(u^2 + r_s) + \lambda^2)c^2dt^2 + 4(u^2 + r_s)du^2 + (u^2 + r_s)^2d\Omega^2$, a polynomial chart through the throat $u = 0$.
Damour and Solodukhin continue through the throat with the proper radial distance $y$, for which $r = r_s + y^2/4r_s$ to leading order; $y = 2\sqrt{r_s}\,u$ to that order, and exactly $2\sqrt{r_s}\,u$ is the height of Flamm's paraboloid.
The exact $y$ and Damour and Solodukhin's tortoise coordinate $z$ have no closed inverse, so neither is a chart with components to print.
$u$ carries the dimension of the square root of a length, which `DIMENSIONS` declares as `L**(1/2)`.

Nandi, Karimov, Izmailov and Potapov put the metric in isotropic form, their (21) and (22): with the areal radius $r(1 + r_s/4r)^2$, $ds^2 = -\left(\left(\frac{4r - r_s}{4r + r_s}\right)^2 + \lambda^2\right)c^2dt^2 + (1 + r_s/4r)^4(dr^2 + r^2d\Omega^2)$.
The throat is $r = r_s/4$, every component is finite there, and $r \to r_s^2/16r$ exchanges the two sides, so $0 < r < r_s/4$ is the whole of the other side and $r = 0$ its far end.

`print_charts.damour_solodukhin_check` pulls Damour and Solodukhin's metric back through each of the four maps and compares it with the chart's own, slot by slot.

## Step 4: the curvature

$G^t{}_t = 0$ in every chart, since $g_{rr}$ is Schwarzschild's: the matter has no energy density.
$G^r{}_r = -\lambda^2r_s/r^2(ar - r_s)$, which is $-1/r_s^2$ at the throat for every $\lambda$, a radial tension $c^4/8\pi Gr_s^2$, the value every Morris-Thorne throat of that radius has.
$G^\theta{}_\theta = \lambda^2r_s(2ar - r_s)/4r^2(ar - r_s)^2$, which is $(1 + 2\lambda^2)/4\lambda^2r_s^2$ at the throat, a tangential pressure that grows as $1/\lambda^2$.
Both fall to the order of $\lambda^2r_s/r^3$ once $r - r_s$ is much more than $\lambda^2r_s$, so the matter sits in a shell of that thickness in $r$.
The Kretschmann scalar is the sum $4A^2 + 8B^2 + 8C^2 + 4D^2$ over the frame components of a static observer, $A = r_s(4(r - r_s)(ar - r_s) - \lambda^2rr_s)/4r^3(ar - r_s)^2$, $B = r_s(r - r_s)/2r^3(ar - r_s)$, $C = -r_s/2r^3$ and $D = r_s/r^3$.
At the throat it is $(1 + 24\lambda^4)/4\lambda^4r_s^4$, and at $\lambda = 0$ it is Schwarzschild's $12r_s^2/r^6$; the script checks the throat's three values and that the Ricci tensor vanishes at $\lambda = 0$.

## Step 5: the diagrams

Every diagram is drawn at $r_s = 1$ and $\lambda = 1/5$.
On the plane of $t$ and $r$ the metric is $(1 - r_s/r + \lambda^2)(-c^2dt^2 + dr_*^2)$ with $dr_*/dr = r/\sqrt{(r - r_s)(ar - r_s)}$, and
$r_* = \sqrt{(r - r_s)(ar - r_s)}/a + \frac{(1 + a)r_s}{2a^{3/2}}\ln\frac{2ar - (1 + a)r_s + 2\sqrt{a(r - r_s)(ar - r_s)}}{\lambda^2r_s}$,
which vanishes at the throat and whose derivative is checked against the published metric.
The tortoise coordinate of the rescaled time is $\sqrt{a}\,r_*$, which is the $x$ of Step 3.
`null_rays.py --verify` checks the rays of each chart against $ct \mp r_*$, $ct \mp \sqrt{a}\,r_*$, $ct \mp x$, $ct \mp r_*\,\mathrm{sgn}(4r - r_s)$ and $ct \mp r_*\,\mathrm{sgn}(u)$.
For small $\lambda$, $\sqrt{a}\,r_* = r + r_s\ln(r/r_s - 1) + r_s(\ln(4/\lambda^2) - 1)$, Bueno and his collaborators' (42), and the constant is Damour and Solodukhin's delay $r_s\ln(1/\lambda^2)$ to leading order.

Through the throat $r_*\,\mathrm{sgn}$ runs over the whole line, so $p, q = \arctan((ct \mp r_*\,\mathrm{sgn})/\ell)$ give the full diamond with the throat on its axis and no horizon; the drawing takes $\ell = 4r_s$, and $\sqrt{a}\,\ell$ in the two charts of the rescaled time, so that one event is one point in all five views, which `conformal.py` checks.

The equator of a moment of $t$ is $dr^2/(1 - r_s/r) + r^2d\phi^2$ for every $\lambda$, Flamm's paraboloid $z^2 = 4r_s(r - r_s)$ on both sides of the throat, which `embedding.py` checks in Damour and Solodukhin's chart, in the isotropic chart and, as $z = 2\sqrt{r_s}\,u$, in Einstein and Rosen's.
The geometry is static, so it is one surface and no movie.
