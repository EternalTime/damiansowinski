# Dust collapse in anti-de Sitter space

The charts are written by `ads_dust_collapse.py` through `print_charts.py`, and the drawings by `ads_dust_collapse_drawings.py`.

## The charts and their sources

Three dimensions, $\Lambda = -1/\ell^2$, from R. B. Mann and S. F. Ross, Phys. Rev. D 47, 3319 (1993).

- `disc_interior` is their (3), with the scale factor left a function of $t$: their (4) is $\ddot a = -a/\ell^2$, their (6) from rest is $a = \cos(ct/\ell)$, and their (9) gives $k = 8\pi G_3\rho_0/c^2 - 1/\ell^2$.
- `disc_exterior` is their (2), the hole of Bañados, Teitelboim and Zanelli without rotation, with their (21), $M = 8\pi G_3\rho_0r_0^2/c^2 - 1$.
- `disc_exterior_ingoing` is the same exterior on the advanced time, $dv = c\,dT + dR/(R^2/\ell^2 - M)$, the form of section 2 of S. Carlip, Class. Quantum Grav. 12, 2853 (1995).
  It is the chart in which the edge crosses the horizon, and $v = 0$ is set where the edge is at rest.

Four dimensions, $\Lambda = -3/L^2$, from W. L. Smith and R. B. Mann, Phys. Rev. D 56, 4942 (1997), and F. C. Mena, J. Natário and P. Tod, Adv. Theor. Math. Phys. 12, 1163 (2008).

- `hyperbolic_interior` is Smith and Mann's (4), with their $k$ written $1/b^2$, and their Friedmann equation (6).
- `hyperbolic_slab` is the same cloud in $r = b\cosh(\chi/b)$, Mena, Natário and Tod's Table 1 (c), from one face through the smallest surface to the other.
- `hyperbolic_exterior` is Smith and Mann's (2), with their (28) for the mass, $\mu = 8\pi G\rho_0r_0^3/(3c^2)$.
- `flat_interior` is Mena, Natário and Tod's section 3.2.1, case (iii), the open universe sliced by horospheres; their $R^2(d\rho^2 + e^{2\rho}(d\theta^2 + d\varphi^2))$ is this chart with $\chi = b\rho$, $x = b\theta$ and $y = b\varphi$.
- `flat_exterior` is the static flat hole of their section 3.1, written with $(R^2/b^2)(dx^2 + dy^2)$ so that $x$ and $y$ are the cloud's own lengths on both sides of the face.

The letter $k$ is the curvature of the disc's space alone, and $b$ the curvature radius of the clouds of four dimensions, so that neither has two meanings.
The scale factor is a pure number with $a = 1$ at the release, and the clouds of four dimensions obey $\dot a^2 = 1/b^2 + (1/L^2 - 1/b^2)/a - a^2/L^2$, with $8\pi G\rho_0/c^2 = 3(1/L^2 - 1/b^2)$.

## The interior that is not a chart

J. P. S. Lemos, Phys. Rev. D 57, 4600 (1998), prints the interior $-dt^2 + a^2(dl^2 + l^2(d\theta^2 + d\phi^2))$ with flat $\theta$ and $\phi$.
Its space is a cone over a plane, and $G^l{}_l - G^\theta{}_\theta = 1/(l^2a^2)$, so it is no dust solution.
Mena, Natário and Tod list the spatially flat universe with the flat hole as a pair that cannot be joined, and F. C. Mena and J. M. Oliveira, Ann. Phys. 387, 135 (2017), prove the same for a radiating exterior.
No chart is taken from it; the flat cloud is the horospherical one, and `test_ads_dust_collapse.py` holds the statement above.

## The junction

Across the surface the metric of the surface is continuous, $R_s = r_0a$, and the surface is the exterior's radial geodesic, $(dR_s/d(c\tau))^2 + f(R_s) = E^2$.
The energy is the interior's: $E^2 = 1 - kr_0^2$ for the disc, $r_0^2/b^2 - 1$ for the hyperbolic cloud and $R_0^2/b^2$ for the flat one, which is the continuity of the extrinsic curvature.
`check` in `ads_dust_collapse.py` holds each interior to dust with its cosmological constant, the slab to Smith and Mann's chart pulled back, and each exterior to its vacuum equations, to the published chart of `btz` or of `topological_black_hole` it is, and to this junction, at points.

## The drawings

The disc is drawn at $k = 1/(4\ell^2)$ and $r_0 = \ell$, so that $M = 1/4$, the horizon is $R = \ell/2$, the edge is $R_s = \ell\cos(ct/\ell)$, it crosses the horizon at $ct = \pi\ell/3$, and the crunch is at $ct = \pi\ell/2$.
In the static time $\tan(ct/\ell) = \sqrt3\tanh(cT/2\ell)$, and in the advanced time $\tan(ct/2\ell) = \sqrt3\,(e^{v/2\ell} - 1)/(3 - e^{v/2\ell})$, so the edge crosses the horizon at $v = \ell\ln(9/4)$ and reaches $R = 0$ at $v = \ell\ln 3$.
The hyperbolic cloud is drawn at $b = 8L/7$ and $r_0 = 2L$, so that $\mu = 15L/8$ and $R_h = 3L/2$, and the flat cloud at $b = \sqrt{8/7}\,L$ and $R_0 = 2L$, so that $\mu = L$ and $R_h = L$.
The static rows of the clouds declare the face as a function of the static time, since a geodesic integrated in that chart loses its digits where the surface gravity is $1.9/L$.

The embedding diagrams are the moments of rest.
The disc is a cap of a sphere of radius $1/\sqrt k$ joined with one tangent to the slice $T = 0$ of its exterior, which lies in flat space out to $R = \sqrt{1 + M}\,\ell$ and in Minkowski space beyond.
The flat cloud is drawn with $x$ closed into a circle, as its faces are closed into tori: its surface $y = 0$ is a pseudosphere of radius $b$, the horn that narrows to the cusp, and the circle's length is chosen so that the exterior lies level at $R = 4L$.
The hyperbolic cloud has no surface of revolution among its totally geodesic sections, which are hyperbolic planes with the slab a band about a geodesic, so it is not embedded, and `exceptions/ads_dust_collapse.json` lists its drawings as holding no embedded moment.

The conformal diagram is the disc's, by the construction of the `oppenheimer_snyder` diagram.
The dust is drawn in $p, q = \arctan(\eta \mp \chi)$, with $\eta = \operatorname{artanh}\sin(ct/\ell)$ its conformal time, which has no end, so the crunch is one point.
Outside, Kruskal's $U$ and $V$ are drawn as $P(U)$ and $Q(V)$, fixed by agreement on the edge, by the moment of rest being the line $T = 0$ and by $R = 0$ being the line $T = \pi$.
The clouds of four dimensions end on a spacelike crunch at a finite conformal time, and their diagram is the one Mena, Natário and Tod draw in Ann. Henri Poincaré 10, 1359 (2010); the spacetime diagram of `hyperbolic_slab` shows the slab between its two faces.

The curvature figure is the disc's stored Kretschmann and Ricci scalars at $a = \cos(ct/\ell)$, uniform over the disc, against its proper time.
