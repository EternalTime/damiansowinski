# The Malament-Hogarth toy spacetime

Malament-Hogarth is a property of a spacetime's causal structure, not a line element.
The line element is the one its references build to show that the property can be had at all.

The components were computed and written by `_tools/derivations/print_charts.py`.
`_tools/derivations/verify_metrics.py --system malament_hogarth/cartesian` recomputes them from the line element in about seven seconds.

---

## Step 1. What the references build

The three papers behind the property build three things.

Etesi and Németi (2002) work in the Kerr black hole, in Boyer-Lindquist coordinates: an observer falls through the outer horizon and meets the Malament-Hogarth event on the inner one, while the computer circles outside on a stable orbit.
Kerr has its own chart in those coordinates.

Hogarth (1992) and Earman and Norton (1993) both use anti-de Sitter space as the example that solves the field equations, and Earman and Norton name Reissner-Nordström as one that meets the energy conditions.
Both have their own charts.

What neither of those carries is the construction Earman and Norton introduce, on page 28 of their paper, as "a useful concrete example of a M-H spacetime", and draw as their Figure 1, "A Toy Malament-Hogarth Spacetime":

> Start with Minkowski spacetime $\mathbb{R}^4$, $\eta_{ab}$ and choose a scalar field $\Omega$ which is 1 outside of a compact set $C$ and which goes rapidly to $+\infty$ as the point $r$ is approached.
> The M-H spacetime is then $M$, $g_{ab}$ where $M = \mathbb{R}^4 - r$ and $g_{ab} = \Omega^2\eta_{ab}$.

Welch (2008) reproduces the same toy with Earman and Norton's figure.
It is the line element here: the removed point $r$ at the origin of Minkowski's inertial chart,

$$ds^2 = \Omega^2\left(-c^2dt^2 + dx^2 + dy^2 + dz^2\right),$$

with $\Omega(t, x, y, z)$ left free, since Earman and Norton fix nothing about it beyond those two conditions.
No coordinates were invented: the chart is Minkowski's own and the only function is theirs.

## Step 2. Why it has the property

A conformal factor leaves every null direction where it was, so the light cones and the causal structure are Minkowski's with one point missing.
An observer at rest in the chart ages $c\,d\tau = \Omega\,c\,dt$.
Let the computer ride the $t$ axis from $t = -T$ toward the origin; that worldline has no future endpoint in $M$, and its proper time $\int_{-T}^0\Omega\,dt$ diverges when $\Omega$ grows at least as fast as $1/|t|$ along it.
An event $p$ on the axis above the origin has the whole of that worldline in its chronological past, because a timelike curve from any point on it can pass around the missing point.
That is the definition the history states.

Earman and Norton add that $\Omega$ can be chosen so that the computer is in free fall, "if $\gamma_1$ is a geodesic of $\eta_{ab}$, choose an $\Omega$ with $\gamma_1$ as an axis of symmetry".
The geodesic equations show it: on the axis $\partial_x\Omega = \partial_y\Omega = \partial_z\Omega = 0$, the three spatial equations are solved by $\dot x = \dot y = \dot z = 0$, and the time equation becomes $\ddot t + (\partial_t\Omega/\Omega)\dot t^2 = 0$, which fixes the affine parameter.

Their equation (5.2) gives the blueshift of a signal from a sender at rest where the factor is $\Omega_1$ to a receiver at rest where it is $\Omega_2$ as $\omega_2/\omega_1 = \Omega_1/\Omega_2$.
It follows from the conformal map alone: a null geodesic of $g$ is one of $\eta$ with tangent $k_g = \Omega^{-2}k_\eta$, a static observer's velocity is $\Omega^{-1}\partial_t$, so the frequency it measures is $\Omega^{-1}$ times the Minkowski one, which is the same at both ends.
With the receiver outside $C$, where $\Omega = 1$, the ratio is $\Omega$ at the sender, and it diverges as the computer nears the missing point.

## Step 3. What the curvature is

The Weyl tensor of any metric conformal to a flat one vanishes, and every Weyl component is zero for that reason, which the checker confirms.
The Ricci scalar is

$$R = -\frac{6\,\Box\Omega}{\Omega^3},\qquad \Box = \partial_x^2 + \partial_y^2 + \partial_z^2 - \partial_t^2.$$

Earman and Norton, on page 34: the toy "can be regarded as a solution to Einstein's field equations with vanishing cosmological constant by computing the Einstein tensor $G_{ab}(g)$ and then defining $T_{ab} \equiv (1/8\pi)G_{ab}$. However, there is no guarantee that even the weak energy condition ... will be satisfied."
The convention says the same.

## Step 4. How the values are written

Each numerator is collected by powers of $\Omega$, so a component reads as $\Omega$ times its second derivatives plus the squares of its first.

The Kretschmann scalar is written through the Schouten tensor $P_{ab}$.
Since the Weyl tensor vanishes, $R_{abcd} = g_{ac}P_{bd} + g_{bd}P_{ac} - g_{ad}P_{bc} - g_{bc}P_{ad}$, and in four dimensions

$$K = 8P_{ab}P^{ab} + 4\left(P^a{}_a\right)^2.$$

For $g = \Omega^2\eta$, $\Omega^2P_{ab} = 2\partial_a\Omega\,\partial_b\Omega - \Omega\,\partial_a\partial_b\Omega - \frac{1}{2}\eta_{ab}\eta^{cd}\partial_c\Omega\,\partial_d\Omega$ and $P^a{}_a = -\Box\Omega/\Omega^3$.
Raising with $\Omega^{-2}\eta^{ab}$ turns $8P_{ab}P^{ab}$ into a signed sum of ten squares over $\Omega^8$, the three with one time index entering with a minus sign.
That form has 1420 characters against 2132 for the expanded numerator, and `print_charts.py` checks it against sympy's contraction.
