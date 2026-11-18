# Kaluza-Klein Melvin universe

## Coordinates and sources

The spacetime is five dimensional Minkowski space with each point identified with the one reached by a step of $2\pi R$ along the fifth dimension and a turn of $2\pi RB$ about an axis.
Fay Dowker, Jerome P. Gauntlett, Gary W. Gibbons and Gary T. Horowitz, arXiv:hep-th/9507143, give the flat metric in cylindrical coordinates, equation (2.1), the identifications, equation (2.2), and the metric in the coordinate $\tilde\varphi = \varphi - Bx^5$, constant along the identified direction, equation (2.3).
Fay Dowker, Jerome P. Gauntlett, Steven B. Giddings and Gary T. Horowitz, arXiv:hep-th/9312172, write the same metric as their equation (5.6), from the Kaluza-Klein form (5.5) of Gibbons and Maeda's solution at the coupling $a = \sqrt3$, and state that it is flat.

The adapted chart is (2.3) and (5.6), with the sources' $x^5$ and $x_5$ written $y$, since the checker reads a superscript as a power, and their $\tilde\varphi$ written $\varphi$.
In it $\varphi$ has period $2\pi$ and $y$ has period $2\pi R$, each identified on its own.
The Minkowski chart is (2.1) with the identification (2.2) as a line of its domains, and its angle is $\phi = \varphi + By$.
`kaluza_klein_melvin.py` beside this file holds both charts; its check holds the adapted chart to the Minkowski chart pulled back by $\phi = \varphi + By$ and every Riemann component of both to zero, at rational points.
Factors of $c$ are restored in the time term, and the components are taken in the chart $x^0 = ct$, as on every page.

The fields of four dimensions, Dowker and coauthors' (2.4)-(2.6), are a metric of another manifold, so they are no chart of this spacetime: the History states them, and `_tools/test_kaluza_klein_melvin.py` rebuilds the adapted chart from them and checks $F_{\mu\nu}F^{\mu\nu} = 2B^2/\Lambda^4$, the statement after (1.5) of hep-th/9312172 at $a = \sqrt3$.

## Drawings

The drawings take $R = 1$ and $B = 1/2R$, the largest field of a distinct spacetime, $-1/2R < B \le 1/2R$ by (2.2).

The embedding diagram is the surface of $\rho$ and $y$ through the axis at $t = 0$ and $z = 0$, on the half planes $\varphi = 0$ and $\varphi = \pi$ of the adapted chart.
In flat space it is a helicoid; its metric $d\rho^2 + (1 + B^2\rho^2)dy^2$ is that of a surface of revolution whose circle, the fifth dimension, has the radius $R\sqrt{1 + B^2\rho^2}$, so it is drawn as one.
Along it $dR/ds = B^2\rho R/\sqrt{1 + B^2\rho^2}$ stays below $BR \le 1/2$, so flat space carries all of it, and the two half planes meet on the axis, the smallest circle, where both sheets stand vertical.
It is drawn out to $\rho = 4R$, where its box is $0.56$ wide for every unit of height.

The spacetime diagram of each chart is the plane of $t$ and $\rho$ through the axis, $z = 0$ and $y = 0$, with $x = \rho$ on $\varphi = 0$ and $x = -\rho$ on $\varphi = \pi$.
The metric on it is $-c^2dt^2 + d\rho^2$ in both charts, and no Christoffel symbol turns a ray with $\dot z = \dot\varphi = \dot y = 0$ out of it, so the plane is totally geodesic, each point a single event, and every ray is at 45°.
The moment of the embedding, $t = 0$, crosses each view from $x = -4R$ to $4R$.

The conformal diagram of each chart brings the same plane into the whole diamond with $p, q = \arctan((ct \mp x)/R)$, as Minkowski's plane $y = z = 0$ is brought, and marks the same moment.

The space is flat, $K = R = 0$ in both charts, so it has no curvature figure: as "Curvature figures" in `_tools/README.md` says, a spacetime without one has no file, no section and no table naming it.

The shared count of `view_text.py` gives each spacetime view, each conformal view and the embedding view fewer than 400 characters.
