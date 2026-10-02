# Tippett and Tsang's time machine

Benjamin Tippett and David Tsang's bubble of 2017 is Minkowski's metric plus a top hat function $h$ times the difference between Rindler's metric and Minkowski's, with Rindler's time read as the polar angle of the plane of $t$ and $x$.
This note records the four charts, what each is checked against, the bubble the diagrams declare, how the drawings take the future in a spacetime that has no arrow of time, what the embedding diagram draws, and why no conformal diagram is drawn.

## Step 1: the sources

The published paper is Class. Quantum Grav. 34, 095006 (2017), doi:10.1088/1361-6382/aa6549, "Traversable acausal retrograde domains in spacetime".
The preprint, arXiv:1310.7985 of October 2013, is titled "Traversable Achronal Retrograde Domains In Spacetime" and differs in its analysis.
The preprint proposes a term $4t^3hW\,dx\,dt$ to give the walls an arrow of time; the published paper drops it and states that the spacetime is not time orientable and has naked singularities.
The charts follow the published paper, whose equations (1) to (6) are the ones cited below.
The lay companion, arXiv:1310.7983, "The Blue Box White Paper", carries the same line element in a footnote.

## Step 2: the Cartesian chart

Equation (1) is $ds^2 = [1 - h\,2t^2/(x^2 + t^2)](-dt^2 + dx^2) + h\,[4xt/(x^2 + t^2)]\,dx\,dt + dy^2 + dz^2$ with $c = 1$.
The chart publishes it with $c$ restored and $h$ left an arbitrary function of the four coordinates, as Alcubierre's chart leaves $f$, so every component is an identity that assumes nothing about the bubble.
With $ct = \xi\sin\lambda$ and $x = \xi\cos\lambda$ the bracket is $1 - 2h\sin^2\lambda$ and the cross term is $2h\sin 2\lambda\,c\,dt\,dx$, so $g = \eta + h\,(g_R - \eta)$, where $g_R = -\xi^2d\lambda^2 + d\xi^2 + dy^2 + dz^2$.
`print_charts.tippett_tsang_bubble` checks that $h = 0$ is Minkowski's metric and that $\det g = -(1 - 4h(1 - h)c^2t^2/(x^2 + c^2t^2))$.
The determinant vanishes exactly where $x = 0$ and $h = 1/2$, and there $g_{tt}$, $g_{xx}$ and $g_{tx}$ all vanish, as the paper's section 2.3 says.
With the paper's bubble the same function checks, at random points in forty digits, that on the plane $y = z = 0$ only $G_{yy}$ and $G_{zz}$ do not vanish, as the preprint's section on the stress tensor says, that $G_{tt} = 0$ on the moment $t = 0$, that $G_{ab}N^aN^b < 0$ for some null $N$ at a point of the wall off that plane, and that the Kretschmann scalar grows without bound toward $x = 0$, $h = 1/2$.

A chart whose components hold the time explicitly names it as `time`, and `chart_printer.Chart.named_time` writes that time as $ct$ while leaving $h$ and its derivatives as the file prints them.
Before this chart no chart had both an explicit time and a free function of it.

## Step 3: the polar chart

Equation (2), $t = \xi\sin\lambda$ and $x = \xi\cos\lambda$, carries the whole metric to $ds^2 = [h + (1 - h)\cos 2\lambda](-\xi^2d\lambda^2 + d\xi^2) - 2(1 - h)\,\xi\sin 2\lambda\,d\lambda\,d\xi + dy^2 + dz^2$.
The paper writes its equation (6) for the null directions in these coordinates: $d\xi = \xi\,d\lambda\,[\sin 2\lambda \pm \sqrt{1 + 2b\cos 2\lambda + b^2}]/(\cos 2\lambda + b)$ with $b = h/(1 - h)$.
`print_charts.tippett_tsang_polar` checks the chart to be the Cartesian one pulled back, to be Rindler's at $h = 1$, and to have equation (6) for its null directions.
The paper's bubble does not depend on $\lambda$ in this chart, and equation (6) writes $h(\xi, y, z)$, so the chart declares $h$ a function of those three: a bubble that runs round a circle of the plane of $t$ and $x$.
Both charts write their Kretschmann scalar as an expanded numerator over the fourth power of the determinant's numerator, since sympy does not finish factoring the numerator of the Cartesian chart's.
The surface $\cos 2\lambda + b = 0$, where the paper's light rays kink, is $g_{\xi\xi} = 0$ here and $g_{xx} = g_{tt} = 0$ in the Cartesian chart, and the diagrams mark it as the zero of $g^{\xi\xi}$ and of $g^{xx}$.

## Step 4: the interior and Rindler's chart

Inside the bubble $h = 1$, and the interior chart is the Cartesian one there, with $\det g = -1$; it is flat, and singular only at the origin of the plane, which the bubble does not reach.
Rindler's chart is equation (3), and `print_charts.tippett_tsang_rindler` checks it to be the interior chart pulled back and the circle of constant $\xi$, $y$ and $z$ to be timelike with the acceleration $1/\xi$ of the paper's section 2.1.
$\lambda$ is an angle of the plane of $t$ and $x$, so it has the period $2\pi$: the interior is the Rindler region of Misner space for a boost of rapidity $2\pi$, the published `rindler` chart of `misner` with $\eta$ identified at $2\pi$.

## Step 5: the declared bubble

The diagrams draw the paper's own bubble, equations (4) and (5), $h = \tfrac12 + \tfrac12\tanh\alpha[R^4 - y^4 - z^4 - (x^2 + t^2 - A^2)^2]$ with $A = 100$, $R = 70$ and $\alpha = 1/6000000$.
In units of $A$ that is $R = 7/10$ and $\alpha = 50/3$, so the argument of the tanh is $4.0$ at the centre of the box, where $h = 0.9993$.
On the plane $y = z = 0$ the walls are the circles $\xi^2 = A^2 \mp R^2$, $\xi = 0.714\,A$ and $1.221\,A$, and the four singular events are where those circles cross $x = 0$.

## Step 6: which way the future is taken

The spacetime with the singular surfaces removed is not time orientable.
Carry the future round the loop that runs inside the ring from $\lambda = 0$ to $\lambda = \pi$ along $\partial_\lambda$, out through the wall along $t = 0$, where the metric is Minkowski's whatever $h$ is, round the outside to the right, and back in along $t = 0$: it returns reversed, since $g(\partial_t, \partial_\lambda) = -x$ changes sign with $x$.
A drawing therefore states its choice, which is Tippett and Tsang's in their figure 3: counterclockwise, along $\partial_\lambda$, where $h > 1/2$, and upward, along $\partial_t$, where $h < 1/2$.
$\partial_\lambda$ is timelike wherever $h > 1/2$, since $g_{\lambda\lambda} = -\xi^2[h + (1 - h)\cos 2\lambda]$, and $\partial_t$ wherever $h < 1/2$, since $g_{tt} = -(1 - 2h\sin^2\lambda)$, so each rule gives a cone wherever it is used.
In `null_rays.py` that is a row's `tau`, a piecewise time function, $\lambda$ inside the ring and $t$ outside; a coordinate named `lambda` is written `lambda_` there, since sympy cannot parse the bare word.
The two null line fields are each continuous round a singular event, and each turns through half a turn there, which is how the arrow is lost while the families stay distinct.
The family named moving left is $x + ct$ constant outside the bubble on both sides.

## Step 7: the diagrams

The plane $y = z = 0$ holds its light rays, since $h$ is even in $y$ and in $z$ and every $\Gamma^y{}_{ab}$ and $\Gamma^z{}_{ab}$ with $a$, $b$ in the plane carries a first derivative of $h$ along $y$ or $z$.
Four flat views draw it: the plane of $t$ and $x$ of the Cartesian chart, the paper's figure 4; the same plane unrolled in the polar chart, where the bubble is a band; the interior chart's plane, whose rays are the spirals $\xi = \xi_0e^{\pm\lambda}$; and Rindler's strip.
The figure in three dimensions, `projections.time_machine_ring`, draws the slice $z = 0$ of $t$, $x$ and $y$: the wall as a tube round the circle $\xi = A$, that circle, a closed timelike curve of proper length $2\pi A$, and cones round it, on the axis and outside.

No conformal diagram is drawn.
The whole spacetime has no arrow of time, so it has no future and no past null infinity, as the paper says of its figure 5, and no conformal diagram of Minkowski's kind.
The interior alone is one copy of Rindler's wedge under a boost of rapidity $2\pi$, which stretches by $e^{2\pi}$ along each light ray and fills the wedge to within a pixel, as Misner space does at $\psi_0 = 4\pi$.

## Step 8: the turn of the cone as a height

On the plane of $t$ and $x$ the metric is the block $\begin{pmatrix}-a & b\\ b & a\end{pmatrix}$ with $a = 1 - 2h\sin^2\lambda$ and $b = h\sin 2\lambda$, which is $N\,[-(\cos\chi\,c\,dt + \sin\chi\,dx)^2 + (\cos\chi\,dx - \sin\chi\,c\,dt)^2]$ with $N = \sqrt{a^2 + b^2}$ and $\tan 2\chi = b/a$.
So the two null directions stand $45°$ on either side of an axis turned by $\chi = \tfrac12\,\mathrm{atan2}(g_{tx}, g_{xx})$ from the $t$ axis toward $-x$: $0$ outside the bubble, and inside it $\lambda$ counted from the nearer half of the $x$ axis.
The embedding diagram draws $\chi$ as a height over the plane $z = 0$ at $ct = A/2$, as the Krasnikov tube draws $1 - k$: no slice is embedded, and the height is the angle alone.
At that moment the boxes stand about $x = \pm\sqrt{3}A/2$, where $\lambda = 30°$ and $150°$, so the centres are turned by $\pm\pi/6$.
A moment of $t$ is spacelike through the declared bubble only while $c^2t^2 < (A^2 - R^2)/2 = 0.255\,A^2$, which $ct = A/2$ meets, and `embedding.tippett_tsang` checks $g_{xx} > 0$ on it.
The moment $t = 0$ itself carries Minkowski's metric exactly, whatever $h$ is, and its extrinsic curvature vanishes, since the metric is even in $t$: the cone is upright everywhere on it.
The moment is marked on the Cartesian and polar planes and on the figure in three dimensions; the interior chart's and Rindler's planes draw the flat interior continued, another spacetime, and do not carry it.
