# Penrose's spherical impulsive wave

The charts of `penrose_impulsive_wave`, why each is the one published, and the maps between them.
`print_charts.py` writes all four, `piw_check` there holds each to what this note derives, and `_tools/test_penrose_impulsive_wave.py` holds the published files to the same maps written out afresh.

## Step 1: the source

Penrose's wave is flat spacetime cut along the future light cone of one event and glued back with a warp, a holomorphic map $Z \to h(Z)$ of the stereographic plane of the wave front.
The continuous form of the metric is equation (15) of J. Podolský and J. B. Griffiths, Class. Quantum Grav. 17, 1401 (2000), arXiv:gr-qc/0001049:
$ds^2 = 2\left|\dfrac{V}{p}\,dZ + U\,\Theta(U)\,p\,\bar H\,d\bar Z\right|^2 + 2\,dU\,dV - 2\epsilon\,dU^2$, with $p = 1 + \epsilon Z\bar Z$, $\epsilon = -1, 0, 1$, and $H$ half the Schwarzian derivative of $h$, their equation (8).
The region $U < 0$ is inside the cone and $U > 0$ outside it.

## Step 2: why the string

The checker reads real coordinates, and a free holomorphic $H$ cannot be written as free real functions without the Cauchy-Riemann equations between them.
Left free, the two real functions would publish a curvature ahead of the wave that the solution does not have.
So the charts take the one warp with a worked physical reading, $h = Z^\beta$, Podolský and Griffiths's equation (16) with $\beta = 1 - \delta$, for which $H = k/Z^2$ with $k = (1 - \beta^2)/4$, their equation (17).
Ahead of the wave the flat space is the cone of a cosmic string with $1 - 4G\mu/c^2 = \beta$, and behind it there is no string.
The History states the general metric.

## Step 3: the null chart, $\epsilon = 0$

With $Z = \rho e^{i\phi}$, $dZ = e^{i\phi}(d\rho + i\rho\,d\phi)$ and $\bar H\,d\bar Z = (k/\rho^2)e^{i\phi}(d\rho - i\rho\,d\phi)$.
So $V\,dZ + U\Theta\bar H\,d\bar Z = e^{i\phi}\left[(V + kU\Theta/\rho^2)\,d\rho + i\rho(V - kU\Theta/\rho^2)\,d\phi\right]$, and the squared modulus is diagonal:
$ds^2 = 2\,dU\,dV + 2(V + kU\Theta/\rho^2)^2d\rho^2 + 2\rho^2(V - kU\Theta/\rho^2)^2d\phi^2$.
The chart defines $\Theta = \tfrac{1}{2}(1 + \mathrm{sgn}(U))$ as a name, so the checker's kink machinery reads it: the Christoffel symbols carry $\Theta$, the Riemann tensor is $\delta(U)$ times a function on the front, $R_{U\rho U\rho} = -2kV\delta(U)/\rho^2 = -R_{U\phi U\phi}/\rho^2$, and the Ricci tensor and the Kretschmann scalar vanish identically.
`chart_printer.step` prints every value once for both sides: as it stands where the sides agree, with $U$ written $U\Theta$ where the value behind the wave is the value ahead of it at $U = 0$, and otherwise as the value behind plus $\Theta$ times the difference.
The chart covers $V > 0$, and ahead of the wave it ends on the string, $V\rho^2 = kU$, where $g_{\phi\phi}$ vanishes.

## Step 4: the retarded chart, $\epsilon = 1$

With $\epsilon = 1$ the cones $U = $ const have their vertices on the timelike line through the centre, as Hogan observed (Phys. Rev. D 49, 6521).
Behind the wave Podolský and Griffiths's (3) gives $ct = (V - 2U)/\sqrt2$, $r = V/\sqrt2$, and $Z = \cot(\theta/2)e^{i\phi}$.
So put $V = \sqrt2\,r$, $U = -u/\sqrt2$, and $\rho = \cot(\theta/2)$ on both sides.
Then $p = 1/\sin^2(\theta/2)$, $d\rho = -d\theta/(2\sin^2(\theta/2))$, $\rho/p = \sin\theta/2$, and $p/\rho = 2/\sin\theta$, and the line element is
$ds^2 = -du^2 - 2\,du\,dr + (r - 2ku\Theta/\sin^2\theta)^2d\theta^2 + \sin^2\theta\,(r + 2ku\Theta/\sin^2\theta)^2d\phi^2$,
with $\Theta = \tfrac{1}{2}(1 - \mathrm{sgn}(u))$, which is $1$ ahead of the wave.
Behind the wave $u = ct - r$ and the metric is Minkowski's in retarded time.
Ahead of it the chart ends on the string, $r\sin^2\theta = -2ku$.

## Step 5: the flat space on either side

Podolský and Griffiths's (4) and (5) with $h = Z^\beta$ have $h''/h' = (\beta - 1)/Z$ and $h''/h' - 2h'/h = -(1 + \beta)/Z$, so
$A = \rho^{1-\beta}/\beta p$, $B = \rho^{1+\beta}/\beta p$, $C = h/(\beta p\rho^{\beta-1})$,
$D = \rho^{1-\beta}\left(p(1-\beta)^2/4\rho^2 + \epsilon\beta\right)/\beta$, $E = \rho^{1+\beta}\left(p(1+\beta)^2/4\rho^2 - \epsilon\beta\right)/\beta$, and $F = hpk/(\beta\rho^{1+\beta})$,
and $v = AV - DU$, $u = BV - EU$, $\eta = CV - FU$ are null coordinates of the cone ahead of the wave, with the angle of $\eta$ equal to $\beta\phi$.
`piw_flat` writes them, and `piw_check` pulls $-2\,du\,dv + 2|d\eta|^2$ back onto each continuous chart on each side at random points in forty digits.
The chart ahead of the wave is the published conical chart of `cosmic_string` at $1 - 4G\mu/c^2 = \beta$, and the chart behind it the published spherical chart of `minkowski`; both are checked.

## Step 6: the plane across the string

On the equator, $\rho = 1$ and $p = 2$, the map of Step 5 with $\epsilon = 1$ is
$cT = (2r + (1 + \beta^2)u)/2\beta$ and $R = (2r + (1 - \beta^2)u)/2\beta$, with $z = 0$.
On the front, $u = 0$, $cT = R = r/\beta$: the moment $t$ behind the wave meets the moment $T = t/\beta$ of the string's rest frame there, and the circle has the circumference $2\pi ct$ from both sides.
That is the surface the embedding diagram draws, a flat disc of radius $ct$ joined to the cone from $R = ct/\beta$ out.
A particle at rest ahead of the wave, $R = R_0$, is the line $dr/du = -(1 - \beta^2)/2$ of the retarded chart, whose plane of $u$ and $r$ has constant metric components and no Christoffel symbol turning a ray out of it, so the line runs on straight behind the front, where $ct = u + r$: $dr/d(ct) = -(1 - \beta^2)/(1 + \beta^2)$.
A ring of such particles falls toward the centre at that speed, $3c/5$ at $\beta = 1/2$.
J. Podolský and R. Steinbauer, Phys. Rev. D 67, 064013 (2003), give $\delta(1 - \delta/2)/(1 - \delta)$ per unit proper time, their equation (ring), which is $\gamma v$ for this $v$ at $\beta = 1 - \delta$.

## Step 7: the drawings

Every drawing takes $\beta = 1/2$, $k = 3/16$.
The conformal diagram of the retarded, behind and ahead charts is one picture: $p = \arctan u$ and $q = \arctan(u + 2r)$ on the equator, which behind the wave is Minkowski's own map and ahead of it is $p = \arctan((cT - R)/\beta)$, $q = \arctan(\beta(cT + R))$, with the string on $\tan q = \beta^2\tan p$.
The null chart's plane $\rho = 1$, $\phi = 0$ meets each moment of the embedding in one event, on the front, at $V = ct/\sqrt2$; the plane $\rho = 1/2$ meets it behind the front at $U = -3ct/4\sqrt2$.
