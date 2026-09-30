# The Milne universe

The four charts of `milne.json` are written by `print_charts.py --metric milne`, which computes every tensor from each line element with the checker's own `Geometry` and reads every printed value back before it writes.
This note derives what the charts are and what the diagrams draw.

## Step 1. The comoving hyperbolic chart

Inside the future light cone of the event $T = R = 0$ of Minkowski spacetime, $cT > R$, write

$$T = t\cosh\chi, \qquad R = ct\sinh\chi,$$

with $t > 0$ and $\chi \ge 0$, and keep the angles $\theta$ and $\phi$.
Then $c\,dT = c\cosh\chi\,dt + ct\sinh\chi\,d\chi$ and $dR = c\sinh\chi\,dt + ct\cosh\chi\,d\chi$, so

$$-c^2dT^2 + dR^2 = -c^2dt^2 + c^2t^2d\chi^2,$$

the cross terms cancelling, and $R^2d\Omega^2 = c^2t^2\sinh^2\chi\,d\Omega^2$, which together give

$$ds^2 = -c^2dt^2 + c^2t^2\left(d\chi^2 + \sinh^2\chi\,d\theta^2 + \sinh^2\chi\sin^2\theta\,d\phi^2\right).$$

This is the FLRW metric with $k = -1$ and scale factor $ct$, the hyperbolic space of radius $ct$ at each moment.
Inverting, $c^2t^2 = c^2T^2 - R^2$ and $\tanh\chi = R/cT$: $t$ is the proper time from the event along the straight world line of constant $\chi$, and $\chi$ is the rapidity of that world line relative to the one at $\chi = 0$, which moves at $R/T = c\tanh\chi$.
Every value of $\chi$ at $t = 0$ is the one event $T = R = 0$.

Every component is taken in the chart $x^0 = ct$.
Since the components depend on $t$, `print_charts.py` names $t$ as the chart's `time`: `Geometry` computes with the symbol standing for $x^0 = ct$, and every value is printed and read back with that symbol written as $c$ times the $t$ the line element prints, which is why $g_{\chi\chi} = c^2t^2$ and $\Gamma^t{}_{\chi\chi} = ct$.

## Step 2. The curvature

The chart covers an open set of Minkowski spacetime, so the Riemann tensor vanishes in every chart, and with it the Ricci tensor, the Ricci and Kretschmann scalars, the Einstein tensor and the Weyl tensor.
`print_charts.py` confirms each component by component in each of the four charts.
As an FLRW model, $\dot a = 1$ in $a = ct$ and $x^0 = ct$, so the Friedmann equation $\dot a^2 + k = (8\pi G/3c^2)\rho a^2$ with $k = -1$ gives $\rho = 0$: it is the empty open universe.

## Step 3. The comoving spherical chart

With $r = \sinh\chi$, $dr = \cosh\chi\,d\chi$ and $d\chi^2 = dr^2/(1 + r^2)$, so

$$ds^2 = -c^2dt^2 + c^2t^2\left(\frac{dr^2}{1 + r^2} + r^2d\Omega^2\right),$$

FRW's comoving chart with $k = -1$ and $a = ct$, and $ctr = R$ is the areal radius.
Its radial rays obey $dr/d(ct) = \pm\sqrt{1 + r^2}/ct$, so $\ln t \pm \operatorname{arcsinh} r$ is constant along them, as $\ln t \pm \chi$ is in Step 1.

## Step 4. The logarithmic time

With $t = t_0e^{\tau/t_0}$ for any positive time $t_0$, $dt = e^{\tau/t_0}d\tau$ and

$$ds^2 = e^{2\tau/t_0}\left(-c^2d\tau^2 + c^2t_0^2\left(d\chi^2 + \sinh^2\chi\,d\Omega^2\right)\right),$$

the static universe with hyperbolic space of radius $ct_0$ times the conformal factor $e^{2\tau/t_0}$, so that $\tau/t_0$ is the conformal time.
$\tau$ runs over the whole line, the event $t = 0$ lying at $\tau \to -\infty$, and $t_0$ fixes only its origin, $\tau = 0$ at $t = t_0$.
In the chart $x^0 = c\tau$ the factor is $e^{2x^0/ct_0}$, which the chart line element in `print_charts.py` writes, and `chart_printer.Printer` sets the exponent inline as $e^{2\tau/t_0}$, the way the collection writes $e^{r^2/R^2}$.

## Step 5. The inertial chart

The chart of Step 1 before the change of coordinates, restricted to the Milne universe: $ds^2 = -c^2dT^2 + dR^2 + R^2d\Omega^2$ with $T > 0$ and $R < cT$.
The published domain `R \in [0, cT)` has an end that depends on $T$; `null_rays.parse_domains` reads such an end as a function of the plane's other coordinate, so the spacetime diagram hatches $R > cT$.

## Step 6. The spacetime diagrams

- `comoving_hyperbolic/through`: rays $\chi = \pm\ln(t/t_1)$, checked by `--verify` to keep $\ln t \pm \chi$.
- `comoving_spherical/radial`: rays keeping $\ln t \pm \operatorname{arcsinh} r$; $|\nabla R|^2 = 1$ for the areal radius $R = ctr$, since $g^{tt}(\partial_{ct}R)^2 + g^{rr}(\partial_rR)^2 = -r^2 + (1 + r^2)$, so no sphere is trapped and no apparent horizon is marked.
- `logarithmic_time/radial`, at $t_0 = 1$: straight rays keeping $\tau \pm t_0\chi$.
- `inertial/through`: straight rays keeping $cT \pm R$, with $R > cT$ hatched.

## Step 7. The conformal diagram

Minkowski's own diagram uses $p = \arctan(u/\ell)$ and $q = \arctan(v/\ell)$ with $u = cT - R$ and $v = cT + R$.
By Step 1, $u = ct\,e^{-\chi}$ and $v = ct\,e^{\chi}$, so the Milne universe, $u > 0$ and $v > 0$, is $p > 0$ and $q > 0$ inside Minkowski's triangle: the wedge with corners $(X, T) = (0, 0)$, $(\pi/2, \pi/2)$ and $(0, \pi)$, bounded by the centre, by the light cone $p = 0$ and by the half of $\mathscr{I}^+$ with $p > 0$.
Each chart enters through its own map: $\chi = \operatorname{arcsinh} r$, $t = e^{\tau}$ at $t_0 = \ell/c = 1$, and Minkowski's map as it stands for the inertial chart.
`conformal.py` checks every map against the published metric, and checks the limits: the line $t = 0$ goes to the event $p = q = 0$, $\chi \to \infty$ at fixed $t$ goes to $\mathscr{I}^+$ at $q = \pi/2$, and the light cone $R = cT$ goes to $p = 0$.
A moment of constant $t$ satisfies $\tan p\tan q = c^2t^2/\ell^2$.

## Step 8. The embedding diagram

On the equator $\theta = \pi/2$ of a moment, $g_{\chi\chi} = c^2t^2$ and $\rho = \sqrt{g_{\phi\phi}} = ct\sinh\chi$, so

$$g_{\chi\chi} - \left(\frac{d\rho}{d\chi}\right)^2 = c^2t^2\left(1 - \cosh^2\chi\right) = -c^2t^2\sinh^2\chi < 0$$

at every $\chi > 0$, and no surface of revolution in flat space carries the slice.
In three dimensional Minkowski space, $dX^2 + dY^2 - dZ^2$, the profile climbs at $dZ/d\chi = \sqrt{(d\rho/d\chi)^2 - g_{\chi\chi}} = ct\sinh\chi$, and started at $Z = ct$ on the axis it is $Z = ct\cosh\chi$: the sheet $Z^2 - \rho^2 = c^2t^2$ of a hyperboloid.
That is $cT$ and $R$ of Step 1 on the plane $\theta = \pi/2$, so the space the sheet is drawn in is the plane $\theta = \pi/2$ of the inertial chart and each moment stands where it lies in spacetime.
The sheets for all $t$ nest inside the light cone $Z = \rho$ of the event $T = R = 0$, drawn as a reference, and the circle of constant $\chi$ moves out along the line $Z = \rho\coth\chi$ through its apex.
`embedding.py` draws the moments $ct = 0.5$, $1$, $2$ and $3$, each out to $\rho = 4$, and a movie with a frame every $0.05$ of $ct$, and checks every frame against $\rho = ct\sinh\chi$ and $Z = ct\cosh\chi$.
