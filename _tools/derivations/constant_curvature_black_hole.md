# The black holes of constant curvature

Anti-de Sitter space of four or more dimensions with its points identified along a boost, Bañados, gr-qc/9703040, Phys. Rev. D 57, 1068 (1998), and in four dimensions Åminneborg, Bengtsson, Holst and Peldán, gr-qc/9604005, Holst and Peldán, gr-qc/9705067, and Bañados, Gomberoff and Martínez, hep-th/9805087.
The charts are written and checked by `constant_curvature_black_hole.py`, and the drawings by `constant_curvature_black_hole_drawings.py`, both beside this file.

## The charts and their sources

Every chart holds the anti-de Sitter radius $l$, a length, with $\Lambda = -3/l^2$ in four dimensions and $-6/l^2$ in five, and Bañados's $r_+$, a length, the radius of each circle on the horizon, where the squared norm of the boost $\xi = \partial_\varphi$ is $r_+^2$.
No coordinate has dimensions of time, so no chart holds $c$.

- `kruskal`, coordinates $(y_0, y_1, y_2, \varphi)$, is Bañados, Gomberoff and Martínez's eqs. (K), (r-y) and (ds/krus): $ds^2 = l^2(r + r_+)^2/r_+^2\,(-dy_0^2 + dy_1^2 + dy_2^2) + r^2d\varphi^2$ with $r = r_+(1 + y^2)/(1 - y^2)$, $y^2 = -y_0^2 + y_1^2 + y_2^2$, $-1 < y^2 < 1$ and $\varphi \sim \varphi + 2\pi$.
  The event horizon is the cone $y^2 = 0$, the singularity the hyperboloid $y^2 = -1$, where $r = 0$, and infinity $y^2 = 1$.
- `static`, coordinates $(t, r, \theta, \varphi)$, is their eqs. (schw1), (ds/sch) and (sph2), $y_0 = f\sin\theta\sinh(r_+t/l)$, $y_1 = f\sin\theta\cosh(r_+t/l)$, $y_2 = f\cos\theta$ with $f = \sqrt{(r - r_+)/(r + r_+)}$, $0 < \theta < \pi$ and $r > r_+$: $ds^2 = (r^2 - r_+^2)(-\sin^2\theta\,dt^2 + (l^2/r_+^2)d\theta^2) + l^2dr^2/(r^2 - r_+^2) + r^2d\varphi^2$.
  Their $t$ is dimensionless, as $\sinh(r_+t/l)$ requires; Creighton and Mann's eq. (1), gr-qc/9710042, is the same chart with their $t$ equal to $r_+$ times this one, a length.
- `exterior`, coordinates $(t, r, \theta, \varphi)$, is their eq. (ds2/sch), $y_0 = f\sinh(r_+t/l)$, $y_1 = f\cos\theta\cosh(r_+t/l)$, $y_2 = f\sin\theta\cosh(r_+t/l)$, which covers the whole exterior: $ds^2 = (r^2 - r_+^2)(-dt^2 + (l^2/r_+^2)\cosh^2(r_+t/l)\,d\theta^2) + l^2dr^2/(r^2 - r_+^2) + r^2d\varphi^2$.
  The arXiv text of (ds2/sch) prints $\cosh(r_+t/l)$ without the square; the pullback of the Kruskal chart along their own map has the square, as Cai's eq. (2.15), hep-th/0206223, and Anand, Campos Delgado and Samart's eq. (2.5), arXiv:2410.16210, print it, and the unsquared form fails the field equations.
- `interior`, coordinates $(r, t, \psi, \varphi)$, is Holst and Peldán's eqs. (SchwarzschildCoord) and (SchwarzschildMetric), with their unit radius restored as $l$ and their period $2\pi\alpha$ of $\psi$ written $2\pi r_+/l$: $ds^2 = l^2(-dr^2/(1 - r^2) + (1 - r^2)dt^2 + r^2d\psi^2 + (1 - r^2)\sinh^2t\,d\varphi^2)$, $0 < r < 1$, $t \ge 0$.
  Their $r$ is the radius of the Kruskal chart divided by $r_+$, read off the length $2\pi r_+r$ of the circle of $\psi$; their $\psi$ is $r_+\varphi/l$ in the other charts' angle; their $\varphi$, the rotation of their $X$ and $Z$, is the exterior chart's $\theta$.
  The singularity is $r = 0$, the event horizon $r = 1$, and the apparent horizon of the slices of constant $r$ their $r = 1/2$, the root of their eq. (0ExpansionEquation).
- `null`, coordinates $(u, v, \psi, \varphi)$, is their eqs. (uvCoordinates) and (uvMetric), $\tan u = e^{-t}\sqrt{(1 - r)/(1 + r)}$, $\tan v = e^t\sqrt{(1 - r)/(1 + r)}$: $ds^2 = l^2\sec^2(u - v)(-4\,du\,dv + \cos^2(u + v)\,d\psi^2 + \sin^2(u - v)\,d\varphi^2)$.
  They state $0 < u < \pi/4$, $0 < v < \pi/2$ and $v > u$; the interior they map is the part $u + v < \pi/2$, since $\tan u\tan v = (1 - r)/(1 + r) < 1$, so the domain is written $u \in (0, \pi/4)$, $v \in [u, \pi/2 - u)$, the singularity $u + v = \pi/2$ and the event horizon $u = 0$.
  The chart is printed with $\alpha = u + v$ and $\beta = v - u$ held as names with the rates $\partial_u\alpha = \partial_v\alpha = 1$, $-\partial_u\beta = \partial_v\beta = 1$, two coordinates with no relation between them, so that every value is a rational function of their sines and cosines.
- `kruskal_five`, coordinates $(y_0, y_1, y_2, y_3, \varphi)$, is Bañados's eqs. (y) and (ds/krus) at $n = 5$, the same as Bañados, Gomberoff and Martínez's eq. (ds/krus5) and Cai's eq. (2.9).
- `static_five`, coordinates $(t, r, \theta, \chi, \varphi)$, is Bañados's eqs. (cc), (ds/sch) and (sph), $y_0 = f\cos\theta\sinh(r_+t/l)$, $y_1 = f\cos\theta\cosh(r_+t/l)$, $y_2 = f\sin\theta\sin\chi$, $y_3 = f\sin\theta\cos\chi$: $ds^2 = (r^2 - r_+^2)(-\cos^2\theta\,dt^2 + (l^2/r_+^2)(d\theta^2 + \sin^2\theta\,d\chi^2)) + l^2dr^2/(r^2 - r_+^2) + r^2d\varphi^2$.
  Bañados prints $0 < \theta < \pi$; the chart's $\partial_t$ vanishes on $\theta = \pi/2$, and $\theta \in [0, \pi/2)$ is the static patch $y_1 > |y_0|$, the range of Bañados, Gomberoff and Martínez's eq. (sph3).
  Their eq. (sph3) prints $-\sin^2\theta\,dt^2$ with $\sin^2\theta\,d\chi^2$, which fails the field equations; Bañados's $\cos^2\theta$ is the pullback of the Kruskal chart.
- `exterior_five`, coordinates $(t, r, \theta, \chi, \varphi)$, is Cai's eqs. (2.13) to (2.15), $y_0 = f\sinh(r_+t/l)$, $y_1 = f\cos\theta\cosh(r_+t/l)$, $y_2 = f\sin\theta\cos\chi\cosh(r_+t/l)$, $y_3 = f\sin\theta\sin\chi\cosh(r_+t/l)$, with $0 \le \theta \le \pi$; Cai says its $t$ is dimensionless.

Left off: Bañados, Gomberoff and Martínez's rotating metric of five dimensions, eq. (dsam), and its extreme limit, eq. (dsext).
Figueroa-O'Farrill, Madden, Ross and Simón, hep-th/0402094, sec. 5.2, show the quotient by a boost combined with a rotation has no event horizon in five dimensions, and the page states nothing that needs the rotating metric.

The reader takes $r_+$ under the spelling `r_plus`, `Reader` in `verify_metrics.py`, since the parser takes no sign in a name.
The radius $r$ of the Kruskal charts is a defined name written out wherever it stands: held with its rates $\partial r/\partial y^2 = (r + r_+)^2/(2r_+)$, the curvature came out constant only through the relation of $r$ with $y^2$, which the rates leave out, and the printed Kretschmann scalar did not read back as $24/l^4$.

## The checks

`check_chart` holds each chart at three events inside its domain, at $l = 7/5$ and $r_+ = 3/5$, in thirty digits.

- Every component of the Riemann tensor is $R_{abcd} = -(g_{ac}g_{bd} - g_{ad}g_{bc})/l^2$, anti-de Sitter space at every point, as Bañados's construction requires.
- The Kruskal charts are the pullback of $-dx_0^2 + dx_1^2 + \dots + dx_{n-1}^2 - dx_n^2$ along Bañados's eq. (y), $x_\alpha = 2ly_\alpha/(1 - y^2)$, $x_{n-1} = (lr/r_+)\sinh(r_+\varphi/l)$, $x_n = (lr/r_+)\cosh(r_+\varphi/l)$, on his hyperboloid (ads/n), $-x_0^2 + x_1^2 + \dots + x_{n-1}^2 - x_n^2 = -l^2$.
- The static and exterior charts of both dimensions are the pullback of the Kruskal chart along their source's map of $y$.
- The interior chart is the pullback of $l^2(dX^2 + dY^2 + dZ^2 - dU^2 - dV^2)$ along Holst and Peldán's eq. (SchwarzschildCoord) at radius $l$, on $X^2 + Y^2 + Z^2 - U^2 - V^2 = -l^2$.
- The null chart is the pullback of the interior chart along $r = (1 - \tan u\tan v)/(1 + \tan u\tan v)$ and $t = \tfrac12\ln(\tan v/\tan u)$.

The printed scalars are those of anti-de Sitter space, $R = -12/l^2$ and $K = 24/l^4$ in four dimensions and $R = -20/l^2$ and $K = 40/l^4$ in five, and every Weyl component vanishes.
`verify_metrics.py` then checks every printed value of every chart.

`metric_tags.py` finds the static charts stationary and static, since $\partial_t$ is timelike where they cover; `OVERRULED` there sets both against them, since that part of the exterior is all they cover and Holst and Peldán show no Killing field is timelike throughout the exterior.

## The drawings

Every drawing takes $l = 1$, the unit of every length, and $r_+ = l$.

- The Kruskal charts are drawn on the plane $y_2 = 0$ (and $y_3 = 0$ in five dimensions) at $\varphi = 0$, the fixed set of two reflections and so totally geodesic, on which the conformal factor $4l^2/(1 - y^2)^2$ leaves light at 45°.
  The row's `where`, $(1 + y_0^2 - y_1^2)(1 - y_0^2 + y_1^2)$, hatches what lies beyond infinity and beyond the singularity.
  The singularity is no curvature singularity, $K = 24/l^4$ there, so it is the edge of the hatching and carries no singular marker, as the spinless hole of three dimensions has none at its $r = 0$.
- The static chart at $\theta = \pi/2$ and the exterior chart at $\theta = 0$ draw the one wedge $y_1 > |y_0|$ of that plane, with the rays $dt = \pm l\,dr/(r^2 - r_+^2)$; the five dimensional charts are drawn on their axis $\theta = 0$, where their metric is the same.
- Holst and Peldán's interior is drawn in $t$ across and $-r$ up, since the future of the black hole runs toward its singularity $r = 0$, with the apparent horizon $r = 1/2$ marked; their null chart is drawn in $(v - u)/2$ across and $(u + v)/2$ up, with the apparent horizon $\tan u\tan v = 1/3$.
- The embedding diagram is the moment $t = 0$ of the exterior chart through $\theta = 0$ and $\theta = \pi$, the line $y_0 = y_2 = 0$ of the Kruskal chart times the circle of $\varphi$, totally geodesic as the fixed set of $t \to -t$ and $y_2 \to -y_2$.
  Its metric $l^2dr^2/(r^2 - r_+^2) + r^2d\varphi^2$ is the moment $t = 0$ of the spinless hole of three dimensions at $M = r_+^2/l^2$, so it is drawn as that hole's is: in flat space from the throat $r = r_+$ to $r = \sqrt{r_+^2 + l^2}$, where it lies level, and in Minkowski space of three dimensions beyond, out to $r = 3l$.
  The throat is the circle in which the moment meets the event horizon, and the two sheets are the halves $\theta = 0$ and $\theta = \pi$ of the one exterior, joined round the axis $r = r_+$ by the angle $\theta$.
- The conformal diagram is the half plane of $y_0$ and $\rho = \sqrt{y_1^2 + y_2^2}$ of the Kruskal chart, $4l^2(-dy_0^2 + d\rho^2)/(1 - \rho^2 + y_0^2)^2$, each point a torus of $\theta$ and $\varphi$, under $p = \arctan(y_0 - \rho)$ and $q = \arctan(y_0 + \rho)$.
  Since $\tan p\tan q = y_0^2 - \rho^2$, infinity lies on $X = \pi/2$, the singularities on $T = \pm\pi/2$, the axis on $X = 0$ and the event horizon on $X = \pm T$: one exterior and one infinity, as Holst and Peldán's Fig. 10 and Åman, Åminneborg, Bengtsson and Pidokrajt's Fig. 1, arXiv:0801.4214, draw it.
  The curves $\rho^2 - y_0^2 = 1/3$ and $3/5$ are $r = 2r_+$ and $4r_+$, and $\rho^2 - y_0^2 = -1/3$ in the black hole is $r = r_+/2$, Holst and Peldán's apparent horizon.
  The view names no chart, since the same rectangle holds in five dimensions with each point a 2-sphere times a circle.
- The moment $t = 0$ is drawn on the Kruskal, static and exterior charts of four dimensions and on the conformal diagram; the exceptions file hides it on the charts of five dimensions, another spacetime, and on the interior charts, which it meets only in the circle $r = r_+$ at their edge.
- No curvature figure is drawn: every chart is anti-de Sitter space at every point, so the Kretschmann and Ricci scalars are the constants above.
