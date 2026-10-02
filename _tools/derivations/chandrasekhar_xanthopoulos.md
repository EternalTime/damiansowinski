# The Chandrasekhar-Xanthopoulos colliding waves

Where both waves have passed, $0 \le |\mu| \le \eta < 1$, the metric is

$$ds^2 = m^2X\left(-\frac{d\eta^2}{1 - \eta^2} + \frac{d\mu^2}{1 - \mu^2}\right) + \frac{Y}{X}\left(dx - \omega\,dy\right)^2 + \frac{(1 - \eta^2)(1 - \mu^2)X}{Y}dy^2 ,$$

with $X = (1 - p\eta)^2 + q^2\mu^2$, $Y = 1 - p^2\eta^2 - q^2\mu^2$, $\omega = 2q\left(\eta(1 - \mu^2) - p(\eta^2 - \mu^2)\right)/Y$ and $p^2 + q^2 = 1$.

All four charts are written by `_tools/derivations/print_charts.py --metric chandrasekhar_xanthopoulos`, and `verify_metrics.py --system chandrasekhar_xanthopoulos/<chart>` checks each.
Every source writes the signature $(+,-,-,-)$, and every chart here is its source's with the signature flipped.

## Step 1. Sources

The solution is S. Chandrasekhar and B. C. Xanthopoulos's, Proc. R. Soc. Lond. A 408, 175 (1986).
The chart of $\eta$ and $\mu$ is theirs, as J. B. Griffiths gives it in eq. (13.22) of *Colliding Plane Waves in General Relativity* (Clarendon Press, 1991), where $\eta$ and $\mu$ are written $t$ and $z$, and as O. Gurtug and M. Halilsoy give it in eq. (77) of arXiv:1509.05174.
The chart of $\psi$ and $\lambda$, $\eta = \sin\psi$ and $\mu = \sin\lambda$, is Griffiths's eq. (A.3), which he ascribes to Chandrasekhar and Ferrari, Proc. R. Soc. Lond. A 396, 55 (1984).
The Boyer-Lindquist chart is Griffiths's eqs. (13.23) to (13.27), with the signs of his $t$, $p$ and $q$ changed as his text after eq. (13.29) says, which is the case of Chandrasekhar and Xanthopoulos: the collision at $r = m$ is followed by the inner horizon.
The ingoing chart is Kerr's own, which is regular on that horizon, in the form R. H. Boyer and R. W. Lindquist give, J. Math. Phys. 8, 265 (1967).

## Step 2. The coordinate x, and units

Chandrasekhar and Xanthopoulos fix the constant in $\omega$ so that it vanishes on $\eta = 1$.
Griffiths's eq. (13.31) shifts $x$ by $2qy/(1 + p)$, after which the cross term vanishes on the collision, and the metric there is $m^2(-d\eta^2 + d\mu^2) + dx^2 + dy^2$, continuous with flat space ahead of the waves.
The charts here use that shifted $x$.
Its $\omega$ is theirs plus $2q/(1 + p)$, which is $2q\left(\eta(1 - \mu^2) - p(\eta^2 - \mu^2)\right)/Y$, and behind one wave alone, $\eta = \mu$, it is $2q\eta$, Griffiths's eq. (13.32).

The sources write every coordinate as a pure number.
Here $\eta$, $\mu$, $\psi$ and $\lambda$ are pure numbers, $x$ and $y$ are lengths, and $m$ is the length that multiplies the part of the metric across the fronts; in Kerr's charts $m$ is $GM/c^2$.
Griffiths's line element carries $X/2$ where this one carries $m^2X$, a rescaling of his $x$ and $y$ by $\sqrt{2}$ that his eq. (13.25) undoes with its factors of $\sqrt{2}M$.

## Step 3. The constraint

The metric is a vacuum only on $p^2 + q^2 = 1$.
The charts take one free parameter, the angle $\alpha$, and define $p = \cos\alpha$ and $q = \sin\alpha$ as names, so that the checker's canonical form, which writes $\cos^2\alpha$ as $1 - \sin^2\alpha$, makes every identity that holds on the circle exact, and no table of relations is needed.
The names $\rho = 1 - p\eta$, $X$ and $Y$ are defined the same way, and every value is printed in them by `CxForms`, which writes a value in $\rho$, reduces the even powers of $p$, factors, and names the factors.
$\rho^2 - 2\rho + q^2 = -p^2(1 - \eta^2)$ is Kerr's $\Delta/m^2$, and $2\rho - X = Y$.

## Step 4. The map onto Kerr

With $a = mq$,

$$r = m(1 - p\eta), \qquad \cos\theta = \mu, \qquad t = x - \frac{2q}{p}y, \qquad \phi = -\frac{y}{mp}$$

carries the chart of $\eta$ and $\mu$ onto Kerr's metric in Boyer and Lindquist's coordinates, with $\Sigma = m^2X$ and $\Delta = -m^2p^2(1 - \eta^2)$.
So $r$ falls from $m$ at the collision to $r_- = m(1 - p)$ on $\eta = 1$, the Killing vectors $\partial_x$ and $\partial_y$ are combinations of $\partial_t$ and $\partial_\phi$, both spacelike there, and $\phi$ runs over the whole line.
The region where both waves have passed is $|\cos\theta| \le (m - r)/\sqrt{m^2 - a^2}$.
The ingoing chart is $dv = dt + (r^2 + a^2)\,dr/\Delta$ and $d\tilde\phi = d\phi + a\,dr/\Delta$, with $r_* = 0$ at $r = m$, so that the event $t = 0$ at the collision is $v = 0$.
`cx_check` holds each chart to a vanishing Ricci tensor, to Kerr's Kretschmann scalar $48m^2(r^2 - a^2\cos^2\theta)(\Sigma^2 - 16a^2r^2\cos^2\theta)/\Sigma^6$, and each chart after the first to being the first carried along its map, at random points to forty digits.

## Step 5. Why there is no double null chart

With null coordinates $u$ and $v$ that vanish on the fronts, $\psi = u + v$ and $\lambda = u - v$, so the double null chart is the chart of $\psi$ and $\lambda$ with $-d\psi^2 + d\lambda^2 = -4\,du\,dv$, and the sines of $u$ and $v$ are the null coordinates of Khan and Penrose that Griffiths's eq. (13.21) uses.
The checker's canonical form spreads $\sin(u + v)$ and $\sin(u - v)$ over the sines and cosines of $u$ and of $v$, and neither that chart nor the one with Khan and Penrose's square roots finished its curvature in a quarter of an hour, where the chart of $\psi$ and $\lambda$ takes four minutes.
The conventions of the charts of $\eta$ and $\mu$ and of $\psi$ and $\lambda$ define $u$ and $v$, and the regions ahead of the waves are the metric with $u$ or $v$ set to zero.

## Step 6. The diagrams

The spacetime diagrams draw the plane $x = y = 0$, which the reflection of $x$ and $y$ together keeps fixed, at $p = 3/5$ and $q = 4/5$; `null_rays.py --verify` compares their rays with the closed forms $\arcsin\eta \pm \arcsin\mu$, $\psi \pm \lambda$ and $\arcsin((m - r)/\sqrt{m^2 - a^2}) \pm \theta$.
The Boyer-Lindquist chart is drawn with $-r$ up, since $r$ falls toward the future, and the ingoing chart draws Kerr's principal null rays on the equator, $\mu = 0$, through the horizon to the ring singularity.
The conformal diagram draws the same plane with all four regions, as Khan and Penrose's and Bell and Szekeres's are drawn, and checks that the metric with $v$ set to zero is a vacuum and curved.
An observer at rest on $\mu = 0$ reaches the horizon after the proper time $m\int_0^1(1 - p\eta)\,d\eta/\sqrt{1 - \eta^2} = m(\pi/2 - p)$.

The embedding diagram is the plane of $x$ and $y$ on $\lambda = 0$, flat at each $\psi$, with a ring of free particles at rest.
Its metric there is

$$h = \frac{1}{1 - p\eta}\begin{pmatrix} 1 + p\eta & -2q\eta \\ -2q\eta & \dfrac{4q^2\eta^2 + (1 - \eta^2)(1 - p\eta)^2}{1 + p\eta} \end{pmatrix}, \qquad \det h = 1 - \eta^2 ,$$

with a cross term, so the ring's ellipse has no fixed axes, and which way it is drawn has to be said.
It is drawn in the frame that parallel transport carries along the world line of the ring's centre: the covectors $E_i$ of that frame obey $dE_{ia}/d\psi = E_{ib}\Gamma^b{}_{\psi a}$ with $\Gamma^b{}_{\psi a} = \tfrac12 h^{bc}\partial_\psi h_{ca}$, which keeps $E^TE = h$, and the ring is $E(\cos\varphi, \sin\varphi)$ for $\varphi$ round the circle.
The ring first stretches along the direction at $-\alpha/2$ from the $x$ axis, the eigenvector of $\partial_\eta h$ at the collision, its area is $\pi\ell^2\cos\psi$, and on the horizon $h$ has rank one and the ring closes onto a segment of half length $\ell\sqrt{(1 + p)(5 - 3p)}/q$, which is $2\sqrt2\,\ell$ at $p = 3/5$.
Between the two the long axis turns through $63°$ at $p = 3/5$, which `embedding.py --verify` checks.
The moments are stacked into the ring's world tube as the other rings of the collection are, and since these ellipses turn, each row of the stack carries the whole frame: the row at $\psi$ is $(a\cos v + s_x\sin v,\; s_y\cos v + b\sin v)$ with $E = \begin{pmatrix} a & s_x \\ s_y & b \end{pmatrix}$, the file's `a`, `sx`, `sy` and `b`.
`_tools/test_chandrasekhar_xanthopoulos.py` holds every row to $E^TE = h$ and $\det E = \cos\psi$ and runs the transport again by itself.
