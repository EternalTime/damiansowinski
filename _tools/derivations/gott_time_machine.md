# Gott's time machine: the charts and what the diagrams draw

J. Richard Gott's two strings (1991) are flat everywhere off the strings, so every chart of `gott_time_machine.json` has vanishing curvature and the physics sits in the domains: which wedges are removed and which events are identified.
This note derives the numbers those domains and the diagrams carry.
Units here have $c = 1$; the file keeps $c$ and $G$.

## Step 1. One string

A string of mass per unit length $\mu$ at rest is Minkowski space with a wedge of angle $2\alpha = 8\pi G\mu/c^2$ removed along it and the faces identified by the rotation through $2\alpha$ about the string.
In the chart `string_rest` the deficit is carried by $g_{\phi\phi} = (1 - 4G\mu/c^2)^2r^2$ with $\phi$ periodic in $2\pi$, the conical chart of `cosmic_string.json`.
The flat angle about the string is $(1 - 4G\mu/c^2)(\phi - \pi)$ measured from the direction $\phi = \pi$, which the chart points at the other string, so the plane $y = 0$ on which Gott joins the two halves is $r\cos[(1 - 4G\mu/c^2)(\phi - \pi)] = d$.

## Step 2. Two strings

Gott takes the string on $x = 0$, $y = d$ with its wedge $|x| < (y - d)\tan\alpha$ opening away from the plane $y = 0$, boosts the half $y \ge 0$ to speed $v$ along $+x$, and joins it along $y = 0$ to its mirror image, boosted to $-v$.
The plane $y = 0$ is carried into itself by both boosts, so the join is smooth.
In the centre of momentum chart the strings are on $x = \pm vt$, $y = \pm d$, and the wedges are $\gamma|x \mp vt| < \pm(y \mp d)\tan\alpha$, with $\gamma = (1 - v^2)^{-1/2}$.
The faces of each wedge are identified at equal times of that string's rest frame, so two identified events differ in the chart's $t$ by $2\gamma v(y \mp d)\tan\alpha$.

## Step 3. The closed timelike curve

Take $A = (t, x, y) = (0, x_0, 0)$ and $B = (0, -x_0, 0)$.
In the rest frame of the upper string, $x_1 = \gamma(x - vt)$ and $t_1 = \gamma(t - vx)$, they are $(\mp\gamma vx_0, \pm\gamma x_0, 0)$.
Close the wedge by turning one side through $2\alpha$ about the string: $A$ and $B$ then lie at the same distance from the string, the face is the bisector of the angle between them, and the straight line joining them crosses the face at right angles.
So the path meets each face at the foot of the perpendicular from the event, at the distance $s = \gamma x_0\sin\alpha - d\cos\alpha$ along the face, and its length is $2(\gamma x_0\cos\alpha + d\sin\alpha)$.
The rocket has the time $2\gamma vx_0$ for it, so its speed in that frame is $(\cos\alpha + (d/\gamma x_0)\sin\alpha)/v$, less than one for large $x_0$ exactly when $v > \cos\alpha$, which is $\gamma\sin\alpha > 1$, Gott's condition.
The crossing events $C$ and $D$ have $t_1 = 0$, $x_1 = \pm s\sin\alpha$, $y = d + s\cos\alpha$, and in the chart $t = \pm\gamma vs\sin\alpha$: the rocket reaches $C$ after $A$ and leaves $D$, the same event, before $B$.
The return round the lower string is the same path turned by $\pi$ about the $t$ axis.
`projections.gott_loop` draws this at $\alpha = \pi/3$, $v = 4/5$, $d = 1/2$, and $x_0 = 3/2$, and checks every stretch timelike against the published metric, $C$ and $D$ on the faces at one rest time and carried onto one another by the rotation, and the curve closed.

## Step 4. The transformation round both strings

Round the upper string the developing map is $P_1 = B(v)\,T(d)\,R(2\alpha)\,T(-d)\,B(-v)$, the rotation about the string in its rest frame, and round the lower one $P_2 = B(-v)\,T(-d)\,R(2\alpha)\,T(d)\,B(v)$, with $B$ a boost along $x$, $T$ a translation along $y$, and $R$ a rotation in the plane of $x$ and $y$.
With $\tan(\alpha/2) = u$ and $e^{\xi} = w$, $v = \tanh\xi$, every entry of $P_2P_1$ is a rational function, and sympy gives, with $X = \gamma^2\sin^2\alpha$:

- the Lorentz part $H$ has $(\mathrm{tr}\,H - 1)/2 = 1 + 8X(X - 1)$, so for $X > 1$ it is a boost of rapidity $a$ with $\cosh a = 1 + 8X(X - 1)$, that is $\cosh(a/2) = 2X - 1$ and $\cosh(a/4) = \gamma\sin\alpha$;
- its axis, the vector $H$ leaves alone, is $(t, x, y) \propto (\cot\alpha/\sinh\xi, 0, 1)$, the $Y$ axis of a frame moving along $y$ at speed $\cot\alpha/\sinh\xi = c^2\cot\alpha/(\gamma v)$, which is below $c$ exactly when $X > 1$;
- the part of the translation along that axis has $b^2 = (4d\sinh\xi\sin\alpha)^2/(X - 1)$, so $b = 4\gamma vd\sin\alpha/(c\sinh(a/4))$.

These are James Grant's (3.12) and (3.13).
The other order, $P_1P_2$, is the same transformation seen from the other side of the strings and has the axis $(-\cot\alpha/\sinh\xi, 0, 1)$.
So away from the strings the spacetime is Minkowski space $(T, X, Y, z)$ with $(T, X, Y)$ identified with its image under the boost of rapidity $a$ in the plane of $T$ and $X$ and the shift $b$ along $Y$: Grant's generalisation of Misner space, which is Misner space at $b = 0$.
The charts `grant_rindler` and `grant_milne` are the Rindler chart of the wedge $X > c|T|$ and the Milne chart of the two light cones of the origin, with $(\eta, Y) \sim (\eta + a, Y + b)$ and $(\chi, Y) \sim (\chi + a, Y + b)$.

## Step 5. Where the closed timelike curves are

The squared interval between $(\eta, \xi, Y)$ and its $n$th image is $n^2b^2 - 4\xi^2\sinh^2(na/2)$.
It vanishes on $\xi_n = nb/(2\sinh(na/2))$, Grant's $n$th polarised hypersurface, his (3.16), and is negative beyond it.
$\xi_n \to 0$ as $n$ grows, so a closed timelike curve passes through every event with $\xi > 0$, and the chronology horizon is $\xi = 0$, the null surfaces $X = \pm cT$.
In the Milne chart the same interval is $n^2b^2 + 4c^2\tau^2\sinh^2(na/2) > 0$, so no closed timelike curve enters the light cones.
The null generators of the horizon keep $Y$, and the identification moves $Y$ by $b$, so no closed null geodesic lies in it.
Dustin Laurence (1994) showed that the region of closed timelike curves of Gott's spacetime is one of the two wedges.

## Step 6. The values the diagrams are drawn at

Every diagram uses $4G\mu/c^2 = 1/3$ ($\alpha = \pi/3$, a wedge of $120°$ on each string), $v = 4c/5$ ($\gamma = 5/3$), and $d = \ell/2$.
Then $\gamma\sin\alpha = 5/(2\sqrt{3}) = 1.443$, $a = 4\,\mathrm{arccosh}(5/(2\sqrt{3})) = 3.640$, $\sinh(a/4) = \sqrt{13/12}$, and $b = 8\ell/\sqrt{13} = 2.219\,\ell$.
The first three polarised hypersurfaces are at $\xi = 0.369$, $0.117$, and $0.028\,\ell$, and the cylinders of Step 7 at $c\tau = -2$, $-1$, $-1/2$, and $-1/4$ times $\ell$ have radii $1.211$, $0.678$, $0.457$, and $0.382\,\ell$, above $b/2\pi = 0.353\,\ell$.

## Step 7. The embedding

A moment $\tau$ of the Milne chart at $z = 0$ has the flat metric $c^2\tau^2d\chi^2 + dY^2$ with $(\chi, Y) \sim (\chi + a, Y + b)$, a flat cylinder of circumference $L = \sqrt{a^2c^2\tau^2 + b^2}$.
The map $\chi = a\phi/2\pi - bx/(c|\tau|L)$, $Y = b\phi/2\pi + ac|\tau|x/L$ carries $\phi$ once round the circle and $x$ along the cylinder, and pulls the metric back to $dx^2 + (L/2\pi)^2d\phi^2$ with no cross term.
`embedding.Slice` takes it through `swept`, which may name the surface's own two coordinates for this.
The circumference shrinks to $b$ at the horizon, where Misner's cylinders, with $b = 0$, close up.
A moment of the centre of momentum time that passes through both strings is a flat sheet with two conical points, which is no surface of revolution, so the construction of `embedding.py` does not draw it.
