# A shell of light falling in anti-de Sitter space

The spacetime is Vaidya's ingoing metric with a negative cosmological constant: light falls inward along the rays of constant advanced time $v$, and the mass function $m(v)$ grows as it arrives.
`ads_vaidya.py` beside this file defines the six charts and checks each one, and `print_charts.py` writes their tensors.

## Sources of the charts

Each chart is read from a fetched source, and the units are restored from it: $c$ multiplies $v$, the anti-de Sitter radius is $\ell$ in every dimension, and $m$ keeps the dimension its source gives it.

1. Three dimensions, ingoing, $(v, r, \theta)$.
   Husain, Phys. Rev. D 50, R2361 (1994), equations (1), (5) and (6) without charge: $ds^2 = -f\,dv^2 + 2\,dv\,dr + r^2d\theta^2$ with $f = -(\Lambda r^2 + g(v))$, $\Lambda = -1/l^2$ and $dg/dv = \rho(v)$, the luminosity of the light.
   The mass function $m$ is Husain's $g$; his $m(v)$ is the radius $l\sqrt{g}$ of the apparent horizon, his equation (9).
   With $g$ a positive constant $M$ the metric is the hole of Bañados, Teitelboim and Zanelli without spin, and anti-de Sitter space is $g = -1$, which is not Husain's vacuum $g = 0$, the hole of no mass.
2. Three dimensions, Poincaré form, $(v, z, x)$.
   Balasubramanian and coauthors, Phys. Rev. Lett. 106, 191601 (2011), equation (1) with $d = 2$: $ds^2 = z^{-2}\left(-(1 - m(v)z^d)dv^2 - 2\,dz\,dv + d\mathbf{x}^2\right)$, with the radius set to 1.
   It is the first chart along $z = \ell^2/r$ and $x = \ell\theta$ with the circle unrolled, and its mass function is $1/\ell^2$ times the first chart's.
3. Four dimensions, round front, $(v, r, \theta, \phi)$.
   Lemos, Phys. Rev. D 59, 044020 (1999), equations (2) and (3), labelled 2.2 and 2.3 in his source: $f = 1 + \alpha^2r^2 - 2m(v)/r$ with $\alpha^2 = -\Lambda/3$, in units $G = c = 1$, and $T_{ab} = (4\pi r^2)^{-1}(dm/dv)k_ak_b$.
   The chart restores $G$ and $c$, so $m$ is a mass, as on the page of Vaidya's metric.
4. Four dimensions, flat front, $(v, r, x, y)$.
   Lemos's equations (6) and (7), labelled 3.1 and 3.2: $f = \alpha^2r^2 - qm(v)/r$ with $r^2(d\theta^2 + d\phi^2)$ on the front, and $q$ a constant fixed by whether the front is a torus, a cylinder or a plane.
   The chart writes his $qm$ as one length $m$ and his pure numbers $\theta$ and $\phi$ as the lengths $x = \ell\theta$ and $y = \ell\phi$, which is also Abajo-Arrastia, Aparício and López's equation (3.1) at $d = 3$ with the radius restored, JHEP 11 (2010) 149.
5. Four dimensions, Poincaré form, $(v, z, x, y)$.
   Balasubramanian and coauthors' equation (1) with $d = 3$.
   It is the fourth chart along $z = \ell^2/r$, and its mass function is $1/\ell^4$ times the fourth chart's.
6. Five dimensions, round front, $(v, r, \chi, \theta, \phi)$.
   Hubeny, Liu and Rangamani, JHEP 01 (2007) 009, equation (5.12) and appendix F: $ds^2 = -f\,dv^2 + 2\,dv\,dr + r^2d\Omega_3^2$ with $f = r^2 + 1 - m(v)/r^2$, $T_{vv} = q(v)/r^3$ and $m' = \tfrac{2}{3}q$ at $\Lambda = -6$.

The article of Danielsson, Keski-Vakkuri and Kruczenski (1999) has no chart here: its shell is of slowly moving matter and it writes no Vaidya metric.

## Geometry

In $D = d + 1$ dimensions the ingoing charts have $f = k + r^2/\ell^2 - \mu(v)/r^{d-2}$, with $k = 1$ for a round front and $0$ for a flat one, and $\Lambda = -d(d - 1)/2\ell^2$.
The Einstein tensor with the cosmological term has one component, $(G + \Lambda g)_{vv} = \tfrac{1}{2}(d - 1)\mu'/r^{d-1}$, with the prime a derivative along $x^0 = cv$, which is Husain's $\rho/2r$, Lemos's $2m'/r^2$ and Hubeny, Liu and Rangamani's $q/r^3$.
The Kretschmann scalar is anti-de Sitter's $2d(d + 1)/\ell^4$ plus the hole's term, $(d - 1)^2(d - 2)\,d\,\mu^2/r^{2d}$, which vanishes in three dimensions: there the metric has constant curvature away from the light, and the singularity at $r = 0$ of the hole is a singularity of the causal structure.
The apparent horizon is where $f = 0$ on each surface of constant $v$ (Husain's equation (9); Abajo-Arrastia and coauthors' (3.5)).
The event horizon is the outgoing ray $dr/d(cv) = f/2$ that ends on the final horizon (their equation (3.7)).

## Checks

`check_chart` in `ads_vaidya.py` holds every chart to the one stress component above.
With a constant mass it holds the ingoing chart of three dimensions to the published ingoing chart of `btz` without spin and, at $m = -1$, to anti-de Sitter space with a regular centre; the chart of four dimensions to the published ingoing chart of `schwarzschild_ads` and, without the cosmological term, to the published ingoing chart of `vaidya`; the flat front to the published ingoing chart of `topological_black_hole` at $k = 0$, whose plane is in polar coordinates; and the chart of five dimensions to the published static chart of `brane_wall_universe` at $k = 1$ carried to advanced time and, without the cosmological term, to the published ingoing chart of `tangherlini`.
Each Poincaré chart is held to the flat chart of its dimension pulled back along $z = \ell^2/r$.
`verify_metrics.py --metric ads_vaidya` then checks every printed value and its dimension.
