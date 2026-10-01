# McVittie

Why the two charts of `mcvittie.json` are the ones published, how their values are grouped, and how its diagrams are built.
`print_charts.py --metric mcvittie` writes both charts and checks everything stated here that is algebra; `conformal.py` and `embedding.py` check the rest where they draw.

## The mass is $r_s$

McVittie, and every paper since, writes the mass as $m = GM/c^2$, so that the lapse is $(1 - m/2ar)/(1 + m/2ar)$.
The collection writes Schwarzschild, Schwarzschild-de Sitter and the global monopole in $r_s = 2GM/c^2$, and McVittie's metric is Schwarzschild's at constant $a$ and Kottler's at constant $H$, so it is written in $r_s$ as well, with $\mu = r_s/4ar$.
There is a second reason, in the checker.
`norm` splits a radical into the roots of the irreducible factors of its radicand, and sympy's `factor_list` fixes the sign of each factor by the order of its symbols.
With the names $R$ and $m$ it returns $\sqrt{R - 2m}$ as $i\sqrt{2m - R}$; with $R$ and $r_s$ it returns $\sqrt{R - r_s}$, which is the reading the areal chart needs outside the throat.

## The isotropic chart

$$ds^2 = -\left(\frac{1 - \mu}{1 + \mu}\right)^2c^2dt^2 + a^2(1 + \mu)^4\left(dr^2 + r^2d\Omega^2\right), \qquad \mu = \frac{r_s}{4ar}.$$

This is McVittie's equation (29), with his $e^{\beta/2}$ written $a$ and his $\mu(t)$, bound by $\dot\beta/2 = -\dot\mu/\mu$, written $m/a$.
The dot on $a$ is $d/d(ct)$, as in FRW's comoving chart, and the chart printer's option `dotted` prints it as FRW writes it.
Every curvature value is collected by $\ddot a$ and $\dot a^2$, and the coefficients then factor: the numerator of $R^t{}_{rtr}$ is $a(4ar + r_s)^7\ddot a - 2r_s(4ar + r_s)^6\dot a^2 + 4096a^5r^3r_s(4ar - r_s)$, the last term Schwarzschild's.
The domain is $r > r_s/(4a)$.
The other half, $0 < r < r_s/(4a)$, is a second copy of the same areal radii that ends on spacelike singularities in the past and in the future, and Kaloper, Kleban and Martin set it aside.

## The areal chart

With $R = ar(1 + \mu)^2$ and the same $t$,

$$ds^2 = -\left(1 - \frac{r_s}{R} - \frac{H^2R^2}{c^2}\right)c^2dt^2 - \frac{2HR}{\sqrt{1 - r_s/R}}\,dt\,dR + \frac{dR^2}{1 - r_s/R} + R^2d\Omega^2,$$

since $\partial_tR|_r = HR(1 - \mu)/(1 + \mu)$, $\partial_rR|_t = a(1 + \mu)(1 - \mu)$ and $\sqrt{1 - r_s/R} = (1 - \mu)/(1 + \mu)$ for $\mu < 1$.
`mcvittie_areal` in `print_charts.py` checks the pullback slot by slot, and `slices.py` checks it again in numbers at the scale factor the diagrams use.
$H = (da/dt)/a$ is a frequency here, with $c$ explicit, as de Sitter's flat slicing has it; the History, the Maths and the diagrams use that one $H$.
The published dot on it is the chart's, $\dot H = dH/d(ct)$, which is why $\dot H$ stands over one power of $c$.
The determinant of the $t$, $R$ block is $-1$, so $g^{RR} = 1 - r_s/R - H^2R^2/c^2$, whose zeros are the apparent horizons.

## The fluid, and the two scalars

$G^t{}_t = -3H^2/c^2$ in both charts: the density is that of the unperturbed universe, uniform on each moment, $8\pi G\rho = 3H^2$.
$G^\theta{}_\theta = -(3H^2/c^2 + 2\dot H/(c\sqrt{1 - r_s/R}))$, the pressure, $8\pi Gp/c^2 = -3H^2 - 2(dH/dt)/\sqrt{1 - r_s/R}$.
The Weyl tensor is Schwarzschild's at the areal radius, $C_{abcd}C^{abcd} = 12r_s^2/R^6$.
With $q = 1/\sqrt{1 - r_s/R}$ and $c = 1$, $R^t{}_t = 3H^2 + 3\dot Hq$ and $R^r{}_r = R^\theta{}_\theta = R^\phi{}_\phi = 3H^2 + \dot Hq$ in the fluid's frame, so

$$R = 12H^2 + 6\dot Hq, \qquad K = C^2 + 2R_{ab}R^{ab} - \tfrac13R^2 = \frac{12r_s^2}{R^6} + 12H^4 + 12\left(H^2 + \dot Hq\right)^2.$$

Both are written that way in both charts and checked against sympy.
At $\dot H = 0$ they are Kottler's $4\Lambda$ and $12r_s^2/R^6 + 8\Lambda^2/3$ with $\Lambda = 3H^2/c^2$.
Where $\dot H \neq 0$ both diverge on $R = r_s$, the Kretschmann scalar as $\dot H^2/(1 - r_s/R)$, the first power only.

## The declared universe

Every diagram takes dust and a cosmological constant, $a = \sinh^{2/3}(3H_0t/2)$ and $H = H_0\coth(3H_0t/2)$, the expansion Lake and Abdelqader chose, at $r_s = 1$ and $H_0 = c/(\sqrt{15}\,r_s)$.
That is $\Lambda r_s^2 = 3H_0^2r_s^2/c^2 = 1/5$, the value Schwarzschild-de Sitter is drawn at, so the late horizons are its $r_- = 1.0852\,r_s$ and $r_+ = 3.2146\,r_s$.
$g^{RR} = 0$ has no root while $Hr_s/c > 2/(3\sqrt3)$, and its two roots appear together on $R = 3r_s/2$ at $ct = 2.0972\,r_s$.

## The singular curve on the spacetime diagrams

`null_rays.py` marks a singular edge where the Kretschmann scalar passes $10^8$ and grows fiftyfold between $10^{-4}$ and $10^{-5}$ of the drawing.
McVittie's grows tenfold per decade and its coefficient $\dot H^2$ falls off as $e^{-6H_0t}$, so in double precision that test sees nothing after the first instants.
A row may therefore declare the curve, `singular_zero="R - r_s"` or `"4*a*r - r_s"`, and `Plot.weak_singularity` draws the zero set and checks the same two conditions at $10^{-20}$ and $10^{-30}$ of the chart's unit from it, in 60 digits, on the side the published domains claim.
Nothing is marked beyond a declared singular curve, since the isotropic chart's formulas run on into the other half.

## The conformal diagram

Put $R = r_s\cosh^2x$, so that $\sqrt{1 - r_s/R} = \tanh x$.
The radial null condition of the areal chart, $dR/d(ct) = \sqrt{1 - r_s/R}\,(HR/c \pm \sqrt{1 - r_s/R})$, becomes, at $r_s = c = 1$,

$$\frac{dx}{dt} = \frac{H(t) \pm \tanh x\,\mathrm{sech}^2x}{2},$$

which is regular from the throat $x = 0$ to $x \to \infty$.
At $x = 0$ both families have $dx/dt = H/2 > 0$: every radial ray leaves the singular sphere, which therefore lies in the past of every event, as Kaloper, Kleban and Martin showed.
So every event has two times, $s_{\rm out}$ and $s_{\rm in}$, at which its outgoing and its ingoing ray left $R = r_s$, with $s_{\rm in} \le s_{\rm out} \le t$, and they are null coordinates.
An outgoing ray gains $x$ all the way, so $s_{\rm out}$ is integrated in $x$, every point in one call; an ingoing ray rises and may fall, so $s_{\rm in}$ is traced back in time to the event $x = 0$.
Both are integrated for $\ln t$, where $tH \to 2/3$ as $t \to 0$, because a ray from far out left at a time exponentially small in $x$.

The drawing takes $p = F(s_{\rm out})$ and $q = -F(s_{\rm in})$, so that $T = p + q \ge 0$ and the singular sphere is the line $T = 0$, with late times on the left.
$F(s) = \arctan(\tfrac12\ln 2s + (e^{\kappa s} - 1)/20)$ spreads the early times by their logarithm and closes on $\pi/2$ as a Kruskal coordinate does, with $\kappa = (1/r_-^2 - 2r_-/15)/2 = 0.352$ the surface gravity of Kottler's $r_-$.
With the logarithm alone the label of $ct = 60\,r_s$ sits at $p = 1.18$, and every late curve stops visibly short of the edge $p = \pi/2$.

What the rays then say, each checked by `conformal.py`:

- An ingoing ray that left before $T_c = 0.0298\,r_s/c$ escapes to $R \to \infty$, and the outgoing rays that cross it stop at a last one, whose label has settled to $10^{-10}$ by $ct = 120\,r_s$. The curve of those last labels is future infinity, spacelike.
- The ingoing ray that left at $T_c$, found by bisection, holds $R = r_+$: the cosmological event horizon.
- Every later ingoing ray ends at $R = r_-$ as $t \to \infty$, where $s_{\rm out} \to \infty$: the null edge $p = \pi/2$, on which the published Kretschmann scalar tends to Kottler's $12/r_-^6 + 24H_0^4$.
- Those that left before about $1.28\,r_s/c$ rise above $r_-$, cross the inner branch of $g^{RR} = 0$ and fall back to it, and reach the edge from the region between the branches, the black hole horizon. Later ones reach it from below, on the part Lake and Abdelqader read as the white hole horizon of the Schwarzschild-de Sitter spacetime. The inner branch ends on the edge between the two, at $q = -F(1.28)$, which is where its own label settles; the overshoot of a ray above $r_-$ falls below $10^{-12}$ for rays that left after $1.3\,r_s/c$, so that time is not found by bisection and no point is marked there.

The tests check each slice on this diagram independently: they invert $F$ at each point, run both rays forward again from the throat by Runge and Kutta's rule, and require them to arrive at one $x$ at the moment's time.

## The embedding

At one $t$ the equator of the isotropic chart is $a^2(1 + \mu)^4(dr^2 + r^2d\phi^2)$, which in $x = ar$ is the equator of Schwarzschild's space in isotropic coordinates.
So every moment is Flamm's paraboloid $z^2 = 4r_s(R - r_s)$ over the areal radius, the same surface for every $a$, from the throat $R = r_s$, where the spacetime is singular and the moment ends.
The movie runs $ct$ from $1$ to $7\,r_s$: the surface stands still, the circles of constant comoving $r$ slide outward along it, and from $ct = 2.10\,r_s$ the two circles of $g^{RR} = 0$ stand on it.
