# The symmetric plane waves of Cahen and Wallach

Cahen and Wallach's Theorem 5(a), Bull. Am. Math. Soc. 76 (1970) 585, writes the solvable indecomposable Lorentzian symmetric space of dimension $n + 2$ as $\mathbb{R} \times \mathbb{R} \times W$ with the metric $2\,ds\,dt + A(w,w)\,ds^2 + (dw,dw)$ at the point $(s,t,w)$, their $W$ negative definite.
In the signature $(-,+,+,+)$ and four dimensions, with $A$ diagonal, this is the Brinkmann chart
$$ds^2 = \left(\lambda_1x^2 + \lambda_2y^2\right)c^2du^2 - 2c\,du\,dv + dx^2 + dy^2,$$
Figueroa-O'Farrill's (12) of arXiv:1112.4967 and Blau and O'Loughlin's (2.5) of hep-th/0212135, with their $A_{ij}$ equal to our $\lambda_i\delta_{ij}$.
Their own $A$ carries the opposite sign, since their transverse metric is negative.
The scaling $(u, v) \mapsto (\mu u, v/\mu)$ multiplies both eigenvalues by $\mu^2$, so a member is the ratio $\lambda_1 : \lambda_2$ with its sign, and $\lambda_1\lambda_2 \neq 0$ keeps it from splitting off a flat line.

## Curvature

In the chart $x^0 = cu$ the only Riemann components are $R_{uxux} = -\lambda_1$ and $R_{uyuy} = -\lambda_2$ and their symmetries, constants, and every term of $\nabla R$ pairs a Christoffel symbol with an index $u$ or $v$ against a component that has none to spare, so $\nabla R = 0$.
The Ricci tensor is $R_{uu} = -(\lambda_1 + \lambda_2)$, the Ricci scalar and the Kretschmann scalar vanish, and the Weyl tensor is $C_{uxux} = -C_{uyuy} = -(\lambda_1 - \lambda_2)/2$.
So the wave is a vacuum at $\lambda_1 + \lambda_2 = 0$, conformally flat at $\lambda_1 = \lambda_2$, and otherwise null dust of energy density $-c^4(\lambda_1 + \lambda_2)/(8\pi G)$ moving along $\partial_v$, positive when the pull outweighs the push.
`cw_check` in `print_charts.py` holds the chart to each of these at rational points.

## Geodesics

The transverse equations are $\ddot{x} = \lambda_1x\,\dot{u}^2$ and $\ddot{y} = \lambda_2y\,\dot{u}^2$ with $u$ affine.
A negative eigenvalue $\lambda = -\mu^2$ gives $x = x_0\cos(\mu cu) + (\dot{x}_0/\mu c)\sin(\mu cu)$, so light leaving the axis at $u = 0$ comes back to it at $u = \pi/\mu c$, Marolf and Ross's first caustic, (3.4) of hep-th/0208197.
A positive one $\lambda = k^2$ gives $\cosh$ and $\sinh$ of $kcu$.

## The rotating chart

For $\lambda_1 = \lambda_2 = -\omega^2$ put $X + iY = e^{i\omega cu}(x + iy)$.
Then $dX^2 + dY^2 = dx^2 + dy^2 + 2\omega c\,(x\,dy - y\,dx)\,du + \omega^2c^2(x^2 + y^2)\,du^2$, and the last term cancels the profile:
$$ds^2 = -2c\,du\,dv + dx^2 + dy^2 + 2\omega c\left(x\,dy - y\,dx\right)du.$$
This is Nappi and Witten's (17) of hep-th/9310112 at $b = 0$, where $\epsilon_{jk}\,da_j\,a_k\,du$ is $\omega = 1/2$ with the orientation reversed, and Halliday and Szabo's (2.7) of hep-th/0502054.
Its Ricci tensor is $R_{uu} = 2\omega^2$, the conformally flat member's $-(\lambda_1 + \lambda_2)$.

## The Rosen chart

For the vacuum member $\lambda_1 = k^2$, $\lambda_2 = -k^2$ take $a_1 = \sinh(kcu)$ and $a_2 = \sin(kcu)$, solutions of $a'' = \lambda a$ in $cu$ that vanish at $u = 0$, Blau and O'Loughlin's (2.9).
With $x = p\,a_1$, $y = q\,a_2$ and $v_B = v + \tfrac{1}{2}\left(p^2a_1a_1' + q^2a_2a_2'\right)$, the cross terms cancel and $(a'^2 + a\,a'')$ cancels the profile, leaving
$$ds^2 = -2c\,du\,dv + \sinh^2(kcu)\,dp^2 + \sin^2(kcu)\,dq^2$$
on $0 < u < \pi/kc$, where $\sin(kcu)$ closes again: the rays leaving the axis at $u = 0$ meet in $y$ at $u = \pi/kc$ and spread in $x$ without end.
`cw_check` holds both charts to the Brinkmann chart pulled back, slot by slot, at rational points.

## Tags

Only the Brinkmann chart describes the whole family, and $\partial_u$ there has norm $(\lambda_1x^2 + \lambda_2y^2)c^2$ of either sign, so the charts cannot decide `stationary`.
$\partial_u + a\,\partial_v$ with $a > 0$ is Killing and timelike about the axis for every member, and `OVERRULED` in `metric_tags.py` carries that reason.

## Drawings

The values drawn are $\lambda_1 = 1/L^2$, $\lambda_2 = -1/L^2$ for the Brinkmann chart, the vacuum member that the photon sphere of a Schwarzschild black hole sees (Fransen, arXiv:2301.06999, section 3), $k = 1/L$ for the Rosen chart and $\omega = 1/L$ for the rotating chart.
