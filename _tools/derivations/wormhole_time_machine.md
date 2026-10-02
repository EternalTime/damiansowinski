# The wormhole time machine: the charts and what the diagrams draw

The three charts of `wormhole_time_machine.json` are written by `print_charts.py --metric wormhole_time_machine`.
This note records where each chart comes from, why its tensors read as they do, and the numbers every diagram carries.
Units here have $c = 1$ and lengths in the throat radius $r_0$; the file keeps $c$.

## Step 1. The three charts and their sources

`lorentz` is the flat space outside the mouths in Lorentz coordinates $(T, X, Y, Z)$.
It is the model of Friedman, Morris, Novikov, Echeverria, Klinkhammer, Thorne and Yurtsever (1990), Section II A: choose two timelike world lines in Minkowski space, mark each with its proper time $\tau$, cut the ball of radius $b$ out of the slice of simultaneity through each point, and identify the two spheres at equal $\tau$.
Every tensor of the chart vanishes, and the physics is in its domains.

`wormhole` is Morris, Thorne and Yurtsever's own chart (1988), inside the wormhole and just outside its mouths,

$$ds^2 = -(1 + glF\cos\theta)^2e^{2\Phi}dt^2 + dl^2 + r^2(d\theta^2 + \sin^2\theta\,d\phi^2).$$

$\Phi(l)$ and $r(l)$ are the static wormhole's, $g(t)$ is the right mouth's acceleration in its own asymptotic rest frame, and $F(l)$ vanishes for $l \le 0$ and rises to 1 at the right mouth.
With $g = 0$ it is the proper radial chart of `morris_thorne.json`.

`short_throat` is one mouth of the 1990 model in its rest frame, $r = b + |l|$:

$$ds^2 = -dt^2 + dl^2 + (b + |l|)^2(d\theta^2 + \sin^2\theta\,d\phi^2).$$

It is the chart through the throat of `thin_shell_wormhole.json` at $r_s = 0$, and the checker reads its kink as `thin_shell_wormhole.md` Step 3 describes.
Its curvature is $R_{\theta l\theta l} = -2b\,\delta(l)$, which is the 1990 paper's equation (2), $R_{alcl} = -(2/b)\gamma_{ac}\delta(l)$ with $\gamma_{\theta\theta} = b^2$, and $G_{tt} = -4\delta(l)/b$, a negative surface energy density.

## Step 2. The lapse as a name

The chart `wormhole` declares $N = 1 + glF\cos\theta/c^2$ as a name, and the checker holds it as a function of $t$, $l$ and $\theta$ while the tensors are built, as it holds Szekeres's $E$.
Every value is then written in $N$, $\Phi$, $r$ and their derivatives, and `Reader.surface` writes the definition out on both sides of each comparison.
Written out in $g$ and $F$ a single Riemann component runs to thirteen terms; in $N$ the longest is four, $R_{tltl} = N(N(\partial_l\Phi)^2 + N\partial_l^2\Phi + 2\partial_l N\,\partial_l\Phi + \partial_l^2N)e^{2\Phi}$.
No component uses $\partial_\theta^2N = \cot\theta\,\partial_\theta N$, which holds for this $N$, so every value stands for any lapse of $t$, $l$ and $\theta$.

## Step 3. Outside the right mouth the chart is flat

Just outside the right mouth $F = 1$, $\Phi = 0$ and $r = l$, and the chart is $-(1 + gl\cos\theta)^2dt^2 + dl^2 + l^2d\Omega^2$, the frame of an accelerated observer.
Morris, Thorne and Yurtsever give the transformation to the Lorentz chart,

$$T = T_R + v\gamma\,l\cos\theta, \quad Z = Z_R + \gamma\,l\cos\theta, \quad X = l\sin\theta\cos\phi, \quad Y = l\sin\theta\sin\phi,$$

with $dT_R/dt = \gamma$ and $dZ_R/dt = v\gamma$ along the mouth's world line.
With the rapidity $\eta$, $v = \tanh\eta$, the pair $(T, Z)$ is $(T_R, Z_R)$ plus $l\cos\theta$ times $(\sinh\eta, \cosh\eta)$, so $-dT^2 + dZ^2 = -(1 + \dot\eta\,l\cos\theta)^2dt^2 + d(l\cos\theta)^2$ and $g = \dot\eta = \gamma^2dv/dt$.
`wormhole_accelerated_frame` in `print_charts.py` checks the pullback in every slot before the chart is written.
Just outside the left mouth $F = 0$ and the same map with $v = 0$ is the identity on $T = t$.
So $t$ is the proper time of both mouths, and that is the whole of the time machine: the wormhole joins the mouths at equal $t$, while the Lorentz time of the right mouth runs ahead of its $t$ after the trip.

## Step 4. The trip the diagrams declare

Every diagram draws one round trip, held by `WormholeTrip` in `null_rays.py`.
The left mouth rests at $Z = 0$ with $\tau = T$.
The right mouth starts at $Z = D = 10$ with $\tau = T = 0$ and has rapidity

$$\eta(\tau) = \tfrac{3}{2}\sin^3(2\pi\tau/P), \qquad 0 \le \tau \le P = 54,$$

and none before or after.
$\eta$ is odd about $P/2$, so $Z_R(P) = D + \int_0^P\sinh\eta\,d\tau = D$: the mouth comes home.
Its greatest speed is $\tanh(3/2) = 0.905$ and its greatest distance $Z = 31.48$.
Its acceleration $g = d\eta/d\tau = (\pi/6)\sin^2(\pi\tau/27)\cos(\pi\tau/27)$ starts and ends at zero and is greatest, $\pi/(9\sqrt{3}) = 0.2015$, at $\tau = (P/2\pi)\arctan\sqrt{2} = 8.21$, and least, $-0.2015$, at $P/2 - 8.21 = 18.79$.
The time shift is $T_R(P) - P = \int_0^P(\cosh\eta - 1)\,d\tau = 21.72$, which exceeds $D$.

Morris, Thorne and Yurtsever ask for $g_{\max}S \ll 1$ over the length $S$ of the wormhole, which makes the lapse differ from 1 by a part too small to draw.
The flat views draw $|l| \le 1.5$, where $gl$ reaches $0.3$, so that the opening of the cones shows.

## Step 5. The closed null geodesic and the horizon

A light ray leaving the left mouth at $\tau$ along $+Z$ meets the right mouth at $T = \tau + Z$.
It arrives at the right mouth's own proper time $\tau$ when $T_R(\tau) - \tau = Z_R(\tau)$, and the first root is $\tau_c = 41.807$, with the right mouth on its way home at $(T, Z) = (59.816, 18.009)$.
Entering there it leaves the left mouth at $\tau_c$, the event it started from: the closed null geodesic $\mathcal{C}$.
Before $\tau_c$ the mouths at one $\tau$ are spacelike separated, which `WormholeTrip` checks on a grid of 5401 times, so no causal curve closes through the throat earlier.
After it they are timelike separated at every time of the same grid, so $\tau_c$ is the only root.

With the mouths small beside $D$, let $f(\tau)$ be the proper time at which light from the left mouth at $\tau$ reaches the right mouth.
For $\tau > \tau_c$ it arrives early, $f(\tau) < \tau$, and steps out of the left mouth at $f(\tau)$, so from the left mouth at any $\tau > \tau_c$ a causal curve reaches the left mouth at every time after $\tau_c$.
An event $p$ therefore lies on a closed timelike curve exactly when it is in the future of the left mouth at $\tau_c$: from there to $p$, on to the right mouth late enough, and back through the throat.
So the Cauchy horizon is the future light cone of the event $\mathcal{C}$ leaves from, as Kim and Thorne describe it ("roughly speaking, a future light cone"), its generators peeling off $\mathcal{C}$.
The right mouth at $\tau_c$ lies on that cone, so its own future adds nothing.
The mouths' finite size bends the generators that pass through the throat, which this picture leaves out.

The figure in three dimensions draws the cone's rim and twelve generators, and a closed timelike curve from the left mouth at $\tau = 60$ straight to the right mouth at the same $\tau$, after the trip, where the two are $D = 10$ apart in space and the time shift $21.72$ apart in $T$: a speed of $0.46$.
The conformal diagram tints the future of the same event on the plane $X = Y = 0$.

## Step 6. The embedding

On $\theta = \pi/2$ the lapse is 1 and the slice $t$ constant has the metric $dl^2 + r^2d\phi^2$, in which neither $g$ nor $F$ appears.
With $r = \sqrt{1 + l^2}$, $dz/dl = \sqrt{1 - (dr/dl)^2} = 1/r$ and $z = \mathrm{arcsinh}\,l$, the catenoid $r = \cosh z$ of `morris_thorne.json`.
`embedding.wormhole_time_machine` checks the surface against that closed form and checks it the same at $t = 0$, $8.21$ and $18.79$.
Both sides open into the one flat space of the chart `lorentz`, a distance apart that no surface of revolution can show, so the drawing is the wormhole alone.
The mouth of the short throat is two flat sheets joined along a circle, with no height to draw.
