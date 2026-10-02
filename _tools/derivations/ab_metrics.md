# The A- and B-metrics of Ehlers and Kundt

The static vacuum fields of type D are seven: Schwarzschild's, which Ehlers and Kundt call A I, the other two A-metrics, the three B-metrics and the C-metric.
The collection has A I as `schwarzschild` and the C-metric as `c_metric`, and this spacetime holds the other five, each in the charts its literature uses.
Every chart takes Ehlers and Kundt's constant $b$, a length, for its one parameter; in Hruška and Podolský's letters $b = 2M$, and $b = 2|n|$.

## Step 1: the letters

Ehlers and Kundt, Martins, and Hruška and Podolský each letter the coordinates their own way, and Ehlers and Kundt use $r$ for a dimensionless coordinate of the surface and $\varphi$ for a length.
The charts here keep one lettering across the family: $r$ or $\sigma$ for the radius, a length; $\chi$, $\theta$, $\phi$ and $x$ on the surfaces, dimensionless; $t$ for a time in seconds; $\tau$ for a dimensionless time on a surface; and $z$ for the length that stands where A I has its time.
A time that is a length is not multiplied by $c$, so only the two charts of $t$ take components in $x^0 = ct$.

## Step 2: the charts and their sources

- `a2_static`: Ehlers and Kundt's A II, their table 2-3.1 as Hruška and Podolský quote it in (2) of arXiv:1808.03508, whose equation numbers the rest of this list uses, $r < b$.
- `a2_cone`: the same line element beyond $r = b$, where the radius is a time, Hruška and Podolský's (32) with $\sigma$ for the radius and $2M = b$.
  `ab_metrics_check` holds it to being the static chart with $r \to \sigma$ and $ct \to z$.
- `a2_kruskal`: $-(4b^3/r)e^{-r/b}\,dU\,dV$ with $(r/b - 1)e^{r/b} = UV$.
  A II's block of $t$ and $r$ is minus Schwarzschild's, so Kruskal's construction goes through with the roles of the regions exchanged: $UV < 0$ is static and $UV > 0$ is not, and the singularity $UV = -1$ is timelike.
  Gott gives extensions of this kind; the chart is checked to carry the static chart's block under $U = -\sqrt{1 - r/b}\,e^{(r - ct)/2b}$, $V = \sqrt{1 - r/b}\,e^{(r + ct)/2b}$.
  The radius is held as a function, as the white hole's Kruskal chart holds it, with `white_hole_kruskal(..., inside=True)` for the changed sign.
- `a2_cartesian`: Hruška and Podolský's (18), the inertial coordinates of the flat spacetime A II is at $b = 0$, inside the cone $T^2 > X^2 + Y^2$, checked to be `a2_cone` carried along $T = \sigma\cosh\chi$, $X = \sigma\sinh\chi\cos\phi$, $Y = \sigma\sinh\chi\sin\phi$.
- `b1_static`: Ehlers and Kundt's B I, Hruška and Podolský's (5).
- `b1_cone`: the same field on the de Sitter slicing of the surfaces of constant $r$ and $z$, Hruška and Podolský's (33), the form Gott and Plebański use, checked to be the static chart carried along Ehlers and Kundt's map $\cos\theta = \cosh\tau\sin\phi$, $\tanh\tau_s = \tanh\tau/\cos\phi$.
- `b1_neck`: Ehlers and Kundt's (2-3.47), Hruška and Podolský's (27) and (39), $r = b/(1 - \rho^2)$, which runs through $r = b$.
- `b1_cartesian`: Hruška and Podolský's (20), the inertial coordinates outside the cone, checked to be `b1_cone` carried along $T = r\sinh\tau$, $X = r\cosh\tau\cos\phi$, $Y = r\cosh\tau\sin\phi$.
- `a3`: Ehlers and Kundt's A III, with the flat plane in polar coordinates $r\chi$ and $\phi$ so that the three A-metrics differ in one function of $\chi$.
- `b2_static`: Ehlers and Kundt's B II, Hruška and Podolský's (6).
- `b2_neck`: Ehlers and Kundt's (2-3.48), Hruška and Podolský's (30), $r = b/(1 + \rho^2)$, with the anti-de Sitter surface in its static coordinates, checked to be the static chart carried along $\sinh\chi' = \sinh\chi\cosh\tau$, $\tan\tau' = \tanh\chi\sinh\tau$.
- `b3`: Hruška and Podolský's (31) with $2n = b$.

Every chart is checked to be a vacuum with $K = 12b^2/r^6$ in its own radius.

## Step 3: the time of the inertial charts

The inertial charts take their time $T$ as a length, the time multiplied by $c$.
Their components hold $T$ outside any function, and a chart whose time is in seconds has to print each such $T$ as $cT$ beside a name whose definition holds $c^2T^2$, which the printer's read-back does not carry through.
With $T$ a length no value needs $c$.

## Step 4: the root outside the cone

`norm` writes $\sqrt{X^2 + Y^2 - T^2}$, whose radicand leads with a negative term in the chart's order of coordinates, as $i\sqrt{T^2 - X^2 - Y^2}$.
Every comparison is sound on that branch and a number taken from it is the wrong one, as `wahlquist.md` Step 5 records.
So `ab_inertial` reads $\sqrt{T^2 - X^2 - Y^2}$ back as $-ir$ before it prints, and `ab_metrics_check` puts the root back on its branch before the pullback is evaluated.
The drawings read the file through the `Reader`, which keeps the root whole.

## Step 5: the drawings

All at $b = 1$.

- Spacetime diagrams, sixteen views, one or two for each chart; `null_rays.py --verify` holds every ray to the null coordinate its caption states.
  On the inertial planes the time function is $T^2 - X^2$ inside the cone and $T/|X|$ outside it, since $T$ itself is no time function near the horizon and the neck.
- Conformal diagrams, eight views. A II is the whole spacetime, each point a hyperbolic plane: Kruskal's hexagon on its side, by $p = \arctan U$, $q = \arctan V$, with the singularity on $X = \pm\pi/2$. B I is its plane $\phi = 0$, $z = 0$, the full diamond by $p, q = \arctan(\tau \mp x)$ with $x = \pm 2\,\mathrm{arcosh}\sqrt{r/b}$, the neck on its axis.
- Embedding diagram: the surface of $r$ and $\phi$ of B I on the de Sitter slicing, Flamm's paraboloid at $\tau = 0$, a movie from $\tau = -0.8$ to $0.8$. Away from $\tau = 0$ the surface lies level at $r = b\cosh^2\tau/\sinh^2\tau$ and stops, Hruška and Podolský's (41) and Fig. 3.

A III, B II and B III have no conformal diagram and no embedding diagram: the literature draws none, and the page's drawings of them are the spacetime diagrams.
