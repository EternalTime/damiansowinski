# The Gowdy cosmologies

Gowdy's vacuum spacetimes have two commuting spacelike Killing vectors, $\partial_\sigma$ and $\partial_\delta$, whose orbits are closed surfaces, and compact spatial sections.
In Gowdy's own notation, which the charts keep,
$$ds^2 = L^2\left\{e^{2a}\left(-dt^2 + d\theta^2\right) + R\left[e^{P}\left(d\sigma + Q\,d\delta\right)^2 + e^{-P}d\delta^2\right]\right\},$$
with $a$, $R$, $P$ and $Q$ functions of $t$ and $\theta$, and $L$ a constant length.
Every coordinate is a pure number, so no coordinate is a time the checker multiplies by $c$, and $L$ is the one dimensionful symbol.

## Step 1: the orbit area

The orbit through $(t, \theta)$ has the area $4\pi^2L^2R$, and in a vacuum $R$ obeys $\partial_\theta^2R - \partial_t^2R = 0$.
The solution $R = t$ with $\theta$ periodic is the torus $T^3$, and $R = \sin t\sin\theta$ is both $S^3$ and $S^2 \times S^1$, which differ in the boundary conditions on the rotation axes $\theta = 0$ and $\theta = \pi$.
Gowdy's Scholarpedia article (2014) is the source for both, and for $L$.

## Step 2: the three charts

`areal` is the torus in the areal time, $e^{2a} = t^{-1/2}e^{\lambda/2}$, the form Ringström's Living Review (2010) works in.
`logarithmic` is the same universe in $\tau = -\ln t$, which reaches $+\infty$ at the singularity: $ds^2 = L^2\{e^{(\lambda + \tau)/2}(-e^{-2\tau}d\tau^2 + d\theta^2) + e^{-\tau}[\dots]\}$, Rendall and Weaver's equation (1) (2001).
`sphere` is $R = \sin t\sin\theta$ with $e^{2a}$ left as it stands, Garfinkle's equation (1) (1999) with $M = 2a$ and $L = P$.
`print_charts.py` checks the logarithmic chart to be the areal one pulled back through $t = e^{-\tau}$.

## Step 3: the sign of $\lambda$

Berger and Moncrief's equation (3.1) (1993) prints $e^{\lambda/2}$ beside constraint equations that hold for $e^{-\lambda/2}$, which Berger and Garfinkle (1998) correct by writing $e^{-\lambda/2}$ in the metric.
The charts keep $e^{+\lambda/2}$, as Ringström, Gowdy, and Rendall and Weaver do, and the constraints are then
$$\partial_t\lambda = t\left[(\partial_tP)^2 + (\partial_\theta P)^2 + e^{2P}\left((\partial_tQ)^2 + (\partial_\theta Q)^2\right)\right], \qquad \partial_\theta\lambda = 2t\left(\partial_tP\,\partial_\theta P + e^{2P}\partial_tQ\,\partial_\theta Q\right),$$
with the opposite overall sign in $\tau$, since $\partial_\tau = -t\,\partial_t$.
`gowdy_check` substitutes the stated equations into every Ricci component of each chart and requires zero, so the signs are sympy's and no paper's.

## Step 4: the field equations of the sphere chart

With $R = \sin t\sin\theta$ the two wave equations pick up $\cot t\,\partial_t - \cot\theta\,\partial_\theta$ in place of $\partial_t/t$, and $R_{t\theta} = 0$ and $R_{tt} + R_{\theta\theta} = 0$ are two linear equations for $\partial_ta$ and $\partial_\theta a$:
$$\cot\theta\,\partial_ta + \cot t\,\partial_\theta a = \tfrac12\left(\cot t\cot\theta + \partial_tP\,\partial_\theta P + e^{2P}\partial_tQ\,\partial_\theta Q\right),$$
$$\cot t\,\partial_ta + \cot\theta\,\partial_\theta a = \tfrac14\left((\partial_tP)^2 + (\partial_\theta P)^2 + e^{2P}\left((\partial_tQ)^2 + (\partial_\theta Q)^2\right) - \cot^2t - \cot^2\theta\right) - 1.$$
Their determinant $\cot^2\theta - \cot^2t$ vanishes on the diagonals $t = \theta$ and $t + \theta = \pi$, where the gradient of $R$ is null and Gowdy's matching constraints take their place.
The inside of Schwarzschild's horizon is a member: with $r = L(1 - \cos t)$, $r_s = 2L$ and $T = L\delta$, $e^{2a} = (1 - \cos t)^2$, $e^{P} = (1 - \cos t)^2\sin\theta/\sin t$ and $Q = 0$, which `gowdy_check` holds to being a vacuum and to being Schwarzschild's interior metric.

## Step 5: the declared wave

The polarised equation $\partial_t^2P + \partial_tP/t - \partial_\theta^2P = 0$ separates into Bessel functions of order zero.
The diagrams declare $Q = 0$, $P = A\,Y_0(t)\cos\theta$ with $A = -\pi/4$, for which
$$\lambda = A^2\left(\tfrac12t^2\left(Y_0^2 + Y_1^2\right) - t\,Y_0Y_1\cos^2\theta\right),$$
by $\tfrac{d}{dt}(tY_0Y_1) = t(Y_0^2 - Y_1^2)$ and $\tfrac{d}{dt}\left[t^2(Y_0^2 + Y_1^2) - tY_0Y_1\right] = t(Y_0^2 + Y_1^2)$.
Each generator checks the wave against the published Einstein tensor.
As $t \to 0$, $Y_0 \to (2/\pi)\ln t$, so $P \to v(\theta)\tau$ with the asymptotic velocity $v = \tfrac12\cos\theta$, between $-\tfrac12$ and $\tfrac12$, and the Kretschmann scalar diverges at every $\theta$.

## Step 6: the Kasner universe at each $\theta$

A homogeneous polarised torus has $P = v\tau = -v\ln t$ and $\lambda = v^2\ln t$, so
$$ds^2 = L^2\left\{t^{(v^2 - 1)/2}\left(-dt^2 + d\theta^2\right) + t^{1 - v}d\sigma^2 + t^{1 + v}d\delta^2\right\}.$$
The proper time is $T \propto t^{(v^2 + 3)/4}$, which makes it Kasner's universe with exponents
$$p_\theta = \frac{v^2 - 1}{v^2 + 3}, \qquad p_\sigma = \frac{2(1 - v)}{v^2 + 3}, \qquad p_\delta = \frac{2(1 + v)}{v^2 + 3},$$
whose sum and whose squares' sum are both 1.
At $v = \pm 1$ they are $(0, 0, 1)$ and $(0, 1, 0)$, the flat Kasner universe.

## Step 7: the diagrams

The plane of the time and $\theta$ is conformally flat in the areal time, so its rays are $t \pm \theta = $ const for every wave, and $e^{-\tau} \pm \theta = $ const in the logarithmic chart; a geodesic that starts in the plane with no momentum along $\sigma$ or $\delta$ keeps none, since both are Killing directions and the metric has no cross term with them, so the rays are null geodesics.
The conformal diagram is the part of Minkowski's triangle between $\theta = 0$ and $\theta = 2\pi$ under $p, q = \arctan(t \mp \theta)$: the singularity is the segment $T = 0$, and since $\theta$ is bounded every ray ends at the one point $i^+$.
The embedding diagram is the torus of $\theta$ and $\sigma$ at one $t$ and $\delta$, a surface of revolution about an axis along $\theta$ with circles of radius $L\sqrt{t}\,e^{P/2}$, drawn as a tube whose two ends are one circle and played as a movie in $t$.
The sphere chart is drawn for the inside of Schwarzschild's horizon of Step 4: its plane of $t$ and $\theta$ is $L^2e^{2a}(-dt^2 + d\theta^2)$, with rays $t \pm \theta = $ const, and since both coordinates are bounded the conformal diagram is the square itself, $p, q = (t \mp \theta)/2$.
The singularity is the edge $t = 0$, where $r = 0$, the edge $t = \pi$ is the horizon, where $K = 12r_s^2/r^6 = 3/(4L^4)$, and the sides are the poles of the sphere.
The gradient of $R = \sin t\sin\theta$ has the norm $\sin(t + \theta)\sin(t - \theta)$ up to a positive factor, so it is null on the two diagonals, timelike in the lower and upper triangles and spacelike in the left and right ones.
The logarithmic chart is drawn against $-\tau$ so that the future is up, which mirrors the plane, and `null_rays.py` swaps the two families there so that the one named first still moves left.

## Step 8: how the printer tells values apart

A chart of three free functions of two coordinates has mixed curvature tensors with some ninety nonzero components, and `chart_printer.py` once compared every pair of them exactly to find the equal and the opposite ones, which took over twelve minutes for one tensor.
`Chart.fingerprint` now evaluates each value at one random point, every symbol, function and derivative a number of its own, and the exact comparison is asked only of two values whose numbers agree or are opposite.
The three charts print in under five minutes.
