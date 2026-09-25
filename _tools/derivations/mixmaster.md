# The Mixmaster universe

The Mixmaster universe is the Bianchi type IX line element with its three scale factors left free functions of time, as Bianchi type I leaves its own.
That the Mixmaster evolution is chaotic is a fact about the solutions of the field equations, not about the line element, which can be written down whatever those solutions do.

The components were computed and written by `_tools/derivations/print_charts.py` in about four minutes, and `_tools/derivations/verify_metrics.py --system mixmaster/euler_angles` recomputes them from the line element in about forty seconds.

---

## Step 1. The line element

Take Euler angles $(\psi, \theta, \phi)$ on the three sphere, $\psi \in [0, 4\pi)$, $\theta \in [0, \pi]$, $\phi \in [0, 2\pi)$, and the three one forms

$$\omega^1 = \cos\psi\,d\theta + \sin\psi\sin\theta\,d\phi,\qquad \omega^2 = \sin\psi\,d\theta - \cos\psi\sin\theta\,d\phi,\qquad \omega^3 = d\psi + \cos\theta\,d\phi.$$

sympy confirms $d\omega^1 = \omega^2\wedge\omega^3$ with its cyclic partners, component by component, so the structure constants are $C^i{}_{jk} = \epsilon_{ijk}$, which is type IX in Bianchi's classification.

The line element is

$$ds^2 = -c^2dt^2 + a_1^2(\omega^1)^2 + a_2^2(\omega^2)^2 + a_3^2(\omega^3)^2,$$

with $a_1$, $a_2$ and $a_3$ free positive functions of $t$.
They are lengths, because the angles are pure numbers; `DIMENSIONS` declares them so.
With all three equal to $a$ the spatial part is $a^2$ times four times the round unit three sphere, so the slice is a sphere of radius $2a$ and the metric is closed Friedmann.
This is the diagonal type IX metric, which is what Misner (1969) and Belinskii, Khalatnikov and Lifshitz (1970) study; the general type IX metric also carries off diagonal terms, and the diagonal case is the one both papers take as their starting point.

## Step 2. The field equations it carries

In the orthonormal frame $c\,dt$, $a_1\omega^1$, $a_2\omega^2$, $a_3\omega^3$ the Ricci tensor is diagonal, and sympy gives

$$R^{\hat 1}{}_{\hat 1} = \frac{\left(a_1'a_2a_3\right)'}{a_1a_2a_3} + \frac{a_1^4 - \left(a_2^2 - a_3^2\right)^2}{2a_1^2a_2^2a_3^2}$$

and its cyclic partners, and

$$G^t{}_t = -\left(\frac{a_1'a_2'}{a_1a_2} + \frac{a_1'a_3'}{a_1a_3} + \frac{a_2'a_3'}{a_2a_3}\right) + \frac{a_1^4 + a_2^4 + a_3^4 - 2a_1^2a_2^2 - 2a_1^2a_3^2 - 2a_2^2a_3^2}{4a_1^2a_2^2a_3^2},$$

with the prime the derivative along $x^0 = ct$.
Setting these to zero is the vacuum, and the conventions state them.
Without the terms free of derivatives they are the Bianchi type I equations whose vacuum is Kasner; those terms are the spatial curvature of the closed slice, Misner's potential walls.

## Step 3. How the values are written

The components are coordinate components, so $\psi$ and $\theta$ appear in them even though the geometry is homogeneous: $\omega^1$ and $\omega^2$ turn with $\psi$.
A canonical form writes $\cos^2\psi$ as $1 - \sin^2\psi$, which hides the structure: $g_{\theta\theta}$ comes out as $a_1^2 - (a_1^2 - a_2^2)\sin^2\psi$ rather than $a_1^2\cos^2\psi + a_2^2\sin^2\psi$.
`chart_printer.py` restores it.
It groups each numerator by powers of $\sin\theta$ and $\cos\theta$, raising every term to the highest degree of its parity with $\sin^2 + \cos^2 = 1$, then does the same inside each group for $\psi$, and keeps $(a_1 - a_2)(a_1 + a_2)$ as $a_1^2 - a_2^2$.

The Ricci scalar is written as Bianchi I's plus the curvature of the slice,

$$R = 2\left(\sum_i\frac{a_i''}{a_i} + \sum_{i<j}\frac{a_i'a_j'}{a_ia_j}\right) - \frac{a_1^4 + a_2^4 + a_3^4 - 2a_1^2a_2^2 - 2a_1^2a_3^2 - 2a_2^2a_3^2}{2a_1^2a_2^2a_3^2},$$

which for $a_1 = a_2 = a_3 = a$ at rest gives $3/(2a^2) = 6/(2a)^2$, the scalar curvature of a sphere of radius $2a$.

The Kretschmann scalar is written through the orthonormal frame.
Nine frame components are independent and nonzero: $R_{\hat 0\hat i\hat 0\hat i} = -a_i''/a_i$, the three $R_{\hat i\hat j\hat i\hat j}$, and three that pair $(\hat 0\hat i)$ with $(\hat j\hat k)$.
The last three vanish in Bianchi I and exist here because the frame is not a coordinate frame, and each carries one time index on one side, so

$$K = 4\left(\sum_i R_{\hat 0\hat i\hat 0\hat i}^2 + \sum_{i<j}R_{\hat i\hat j\hat i\hat j}^2\right) - 8\sum_i R_{\hat 0\hat i\hat j\hat k}^2.$$

That form has 1042 characters against 1666 for the expanded one, and `print_charts.py` checks it against sympy's contraction before writing it.

## Step 4. What is not drawn, and why

A diagram of null rays needs the scale factors, and the only honest ones are numerical solutions of the equations of Step 2, which are chaotic.
`_tools/README.md` records that with the other undrawn spacetimes.
