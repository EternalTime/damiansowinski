# The black hole threaded by a cosmic string

Aryal, Ford and Vilenkin's metric of 1986 is Schwarzschild's with the angle about the axis shortened by $b = 1 - 4G\mu/c^2$:

$$
ds^2 = -\left(1 - \frac{r_s}{r}\right)c^2dt^2 + \frac{dr^2}{1 - \dfrac{r_s}{r}} + r^2\left(d\theta^2 + b^2\sin^2\theta\,d\phi^2\right).
$$

## Step 1. The four charts

The static chart keeps $\phi$ periodic in $2\pi$ and carries the deficit in $g_{\phi\phi}$.
The chart with the wedge removed takes $\tilde\phi = b\phi$, periodic in $2\pi b$, and its line element is Schwarzschild's own; `print_charts.py` checks it slot by slot as the static chart pulled back through $\phi = \tilde\phi/b$.
The two Eddington-Finkelstein charts take Schwarzschild's $v, u = ct \pm (r + r_s\ln|r/r_s - 1|)$, since the string enters $g_{\phi\phi}$ alone and the plane of $t$ and $r$ is Schwarzschild's.
Every curvature component is Schwarzschild's with a factor $b^2$ beside each $\sin^2\theta$ of a lowered $\phi$, the Ricci tensor vanishes off the axis, and $K = 12r_s^2/r^6$.
The parameters are $r_s$, a length, and $b$, a number, as `DIMENSIONS` in `verify_metrics.py` declares.

## Step 2. The axis

Near $\theta = 0$ a sphere of constant $t$ and $r$ has the metric $r^2(d\vartheta^2 + b^2\vartheta^2d\phi^2)$, a cone with the deficit angle $2\pi(1 - b) = 8\pi G\mu/c^2$, and the same holds at $\theta = \pi$.
The two halves of the axis carry the same tension, so no chart of the family has an acceleration, in contrast with the C-metric, where they differ.

## Step 3. The embedding diagram

On the equator at $t = 0$ the metric is $dr^2/(1 - r_s/r) + b^2r^2d\phi^2$, so $\rho = br$ and

$$
\frac{dz}{dr} = \sqrt{\frac{r}{r - r_s} - b^2}, \qquad
z = r_s\left(w\sqrt{1 + k^2w^2} + \frac{\mathrm{arsinh}(kw)}{k}\right), \quad w = \sqrt{r/r_s - 1}, \quad k = \sqrt{1 - b^2}.
$$

At $b = 1$ that is Flamm's $2r_sw$, and far out $dz/d\rho \to k/b$, the cone of the string, of half angle $\arcsin b$.
The horizon at one moment has the metric $r_s^2(d\theta^2 + b^2\sin^2\theta\,d\phi^2)$: $\rho = br_s\sin\theta$ and $dz/d\theta = r_s\sqrt{1 - b^2\cos^2\theta}$, a spindle of length $\pi r_s$ along a meridian and of area $4\pi br_s^2$, which is the $A$ of $S = A/4$.

## Step 4. The other two diagrams

The spacetime diagrams are Schwarzschild's planes with the string's parameter held at $b = 0.9$, and the conformal diagram is Kruskal and Szekeres's hexagon, each point a sphere of radius $r$ with the wedge missing.
The bifurcation sphere is the point $t = 0$, $r = r_s$ of the static planes and the point $(0, 0)$ of the hexagon, and lies off both Eddington-Finkelstein charts.
