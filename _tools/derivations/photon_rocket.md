# Kinnersley's photon rocket

Why each chart is the one published, and what `print_charts.py`, `null_rays.py`, `conformal.py` and `embedding.py` check before they write.

## Step 1: the rectilinear chart

The line element is von der Gönna and Kramer's (58) and Podolský's (12) of 2008 with $\Lambda = 0$, in the signature $(-,+,+,+)$:
$ds^2 = -(1 - 2m/r - 2\alpha r\cos\theta)c^2du^2 - 2c\,du\,dr + r^2(d\theta + \alpha\sin\theta\,c\,du)^2 + r^2\sin^2\theta\,d\phi^2$.
Expanding the square gives Podolský's $g_{uu} = -(1 - 2m/r - 2\alpha r\cos\theta - \alpha^2r^2\sin^2\theta)$ and $g_{u\theta} = \alpha r^2\sin\theta$.
$m = GM/c^2$ is a length and $\alpha$ an inverse length, and both are free functions of $u$ in every tensor.
Podolský's transformation (22) to Minkowski coordinates at $m = 0$ has $X = x(u) + r(-\cos\theta\cosh w + \sinh w)$ with $\dot x = \sinh w$ and $w = \int\alpha\,du$, so $\theta = 0$ points against the acceleration.

## Step 2: the Robinson-Trautman chart

The line element is Podolský's (52) and (53) of 2011 in four dimensions with $\Lambda = 0$:
$ds^2 = -(1 - 2m/r - 2r\,\partial_u p/p)c^2du^2 - 2c\,du\,dr + r^2p^{-2}(d\theta^2 + \sin^2\theta\,d\phi^2)$, with $p = U_0 - U_1\cos\theta - U_2\sin\theta\cos\phi - U_3\sin\theta\sin\phi$.
$U_i$ are the spatial components of the rocket's four-velocity over $c$ in the background inertial frame, Podolský's $\dot z^i$, and $U_0 = \sqrt{1 + U_1^2 + U_2^2 + U_3^2}$, so the normalisation holds identically and the three $U_i$ with $m$ are Kinnersley's four arbitrary functions of time.
The checker holds $p$ as a function of $u$, $\theta$ and $\phi$ while the tensors are built.
`RocketForms.reduce` uses that $p$ and each of its derivatives along $u$ is a constant on the sphere less a function linear in the direction of the ray, $q = A - \mathbf{B}\cdot\mathbf{n}$, for which $\partial_\theta^2 q = A - q$, $\partial_\theta\partial_\phi q = \cot\theta\,\partial_\phi q$ and $\partial_\phi^2 q = (A - q)\sin^2\theta - \sin\theta\cos\theta\,\partial_\theta q$.
The constant of $p$ is $U_0$, and the unit four-velocity gives $2U_0p = p^2 + 1 + (\partial_\theta p)^2 + (\partial_\phi p)^2/\sin^2\theta$, whose derivatives along $u$ are the constants of the derivatives of $p$.
That leaves $p$, $\partial_\theta p$, $\partial_\phi p$ and their derivatives along $u$ as generators with no relation among them, so a value that vanishes for the rocket's $p$ is exactly zero.

## Step 3: what is checked before a chart is written

`photon_rocket_check` holds each chart to $R_{uu} = 2(3m\,\partial_u\ln p - \partial_u m)/r^2$ as the one component of the Ricci tensor, with $\partial_u\ln p = \alpha\cos\theta$ in the rectilinear chart, Podolský's (13) and (7); to the Kretschmann scalar $48m^2/r^6$; and to a vanishing Riemann tensor at $m = 0$.
`photon_rocket_pullback` holds the rectilinear chart to being the Robinson-Trautman chart with $U = (\sinh w, 0, 0)$ and $\alpha = \partial_u w$, pulled back through the aberration $\tan(\theta'/2) = e^{-w}\cot(\theta/2)$, von der Gönna and Kramer's (59) with the angle measured from the other pole.

## Step 4: the declared burn

Every diagram is drawn for one burn, in units of the mass $m_0$ before it: $\alpha = 2\sin^2(\pi cu/10m_0)/(25\,m_0)$ for $0 < cu < 10\,m_0$, so $w = (cu - 5\sin(\pi cu/5)/\pi)/25$ ends at $2/5$, a speed of $\tanh(2/5) = 0.380\,c$.
The density of the radiation is proportional to $3\alpha m\cos\theta - \partial_u m$, positive in every direction where $\partial_u m \le -3|\alpha|m$, and the burn takes equality, $m = m_0e^{-3w}$, which ends at $m_0e^{-6/5} = 0.301\,m_0$.
The greatest $\alpha m$ is $0.0501$ at $cu = 3.95\,m_0$, below $1/16$, so on the axis behind the rocket $g^{rr} = 1 - 2m/r - 2\alpha r$ keeps two separate zeros: the outer comes in to $4.73\,m_0$ at $cu = 4.59\,m_0$ and the inner swells to $2.17\,m_0$ at $cu = 2.70\,m_0$ before it shrinks to $0.60\,m_0$.
Ahead of the rocket the one zero of $1 - 2m/r + 2\alpha r$ lies inside $2m$, at $0.854$ of it where $\alpha m$ is greatest.

## Step 5: the three diagrams

On the axis $\sin\theta = 0$ removes $g_{u\theta}$ and every $\Gamma^\theta$ of the plane of $u$ and $r$, so the null curves `null_rays.py` traces there are null geodesics, in both charts.
The conformal diagram's map is derived in `RocketPlane`'s docstring in `conformal.py` and checked there against the published metric of both charts.
The embedding diagram is the surface $\theta = \pi/2$ at constant $cu + r = T$: with $du = -dr/c$ the rectilinear line element gives $(1 + 2m/r + \alpha^2r^2)\,dr^2 + r^2d\phi^2$, a surface of revolution with $dz/dr = \sqrt{2m/r + \alpha^2r^2}$.
