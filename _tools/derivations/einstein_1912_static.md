# Einstein's static field of 1912

A. Einstein's Prague theory of static gravity: space is flat and the speed of light varies from place to place, serving as the gravitational potential.
This note records the four charts, the source of each, what each is checked against, and what the drawings draw.

## Step 1: the sources

The three papers are Annalen der Physik 340, 898 (1911), doi:10.1002/andp.19113401005, Annalen der Physik 343, 355 (1912), doi:10.1002/andp.19123430704, and Annalen der Physik 343, 443 (1912), doi:10.1002/andp.19123430709; each record was checked on Crossref on 2 October 2026, and each names A. Einstein as its one author.
His reply to Max Abraham is Annalen der Physik 343, 1059 (1912), doi:10.1002/andp.19123431014, checked the same way.
Domenico Giulini's study of the second equation is doi:10.1007/978-3-319-06761-2_10, arXiv:1306.5966, whose source was read in full.
Norbert Straumann's account of the Zurich notebook is Annalen der Physik 523, 488 (2011), doi:10.1002/andp.201110467, arXiv:1106.0900, and Galina Weinstein's of the letters and of the polemic with Abraham is arXiv:1202.2791; both were read in full, and every quotation of the page's history is in the words of one of them.
The texts of Einstein's own papers are not openly readable, so the history cites these three for what the papers say.

## Step 2: the letter N

Einstein writes $c$ for the speed of light at a place and $c_0$ for a constant.
The collection keeps $c$ for the constant, so the charts write the speed of light at a place as $cN$, with $N$ a pure number, and the history keeps Einstein's letters.
His first equation, $\Delta c = kc\rho$ with $k = 4\pi G/c^2$, is $\nabla^2N = 4\pi G\rho N/c^2$, and his second, Giulini's (3), is $N\nabla^2N - \tfrac{1}{2}(\nabla N)^2 = 4\pi G\rho N^2/c^2$, which is linear in $\sqrt{N}$, Giulini's $\Psi$.

## Step 3: the free chart

The second paper's law of motion is that $\int\sqrt{c^2dt^2 - dx^2 - dy^2 - dz^2}$ is stationary, so a body moves on a geodesic of $ds^2 = -N^2c^2dt^2 + dx^2 + dy^2 + dz^2$.
The chart `static` leaves $N$ a free function of the three lengths.
`print_charts.einstein_1912_check` holds its Ricci scalar to $-2\nabla^2N/N$, its $R_{tt}$ to $N\nabla^2N$, and its Einstein tensor to no $tt$ component, which is what flat space means for it.

## Step 4: the uniform field

The first paper's accelerated frame has $c = c_0 + ax$, a solution of Laplace's equation.
The chart `uniform` counts the height from the plane where it vanishes, $N = az/c^2$, so that a body at rest at $z = c^2/a$ has the proper acceleration $a$.
It is checked to be flat and is the published Rindler chart of `minkowski` under $T \to t$, $X \to z$, which the test file holds.
It does not solve the second equation, $N\nabla^2N - \tfrac{1}{2}(\nabla N)^2 = -a^2/2c^4$, which is why Einstein wrote that the equivalence principle holds only for infinitely small fields.

## Step 5: the two bodies

Outside a spherical body the first equation is Laplace's for $N$, so $N = 1 - m/r$ with $m = GM/c^2$, the chart `february`; the paper of 1911 has the same $c = c_0(1 + \Phi/c^2)$.
The second equation is Laplace's for $\sqrt{N}$, so $\sqrt{N} = 1 - R_g/r$ with Giulini's gravitational radius $R_g = GM/2c^2 = m/2$, the chart `march`.
Each is checked to be the free chart at its own $N$ and to solve its equation, $R = 0$ in February and $R = -(\nabla N)^2/N^2$ in March.
Space is flat, so $\gamma = 0$ and light bends by half of general relativity's angle; $-g_{tt} = 1 - 2U + 2\beta U^2$ with $U = m/r$ has $\beta = 1/2$ in February and $3/4$ in March, so a perihelion advances by $(2 + 2\gamma - \beta)/3 = 1/2$ and $5/12$ of general relativity's advance.
`_tools/test_einstein_1912_static.py` runs a ray and a planet with the published Christoffel symbols to those numbers.
The lapse vanishes on $r = m$ and on $r = m/2$, where the Kretschmann scalar diverges as $(r - m)^{-2}$ and $(r - m/2)^{-4}$; light takes an infinite time $t$ to reach either sphere and a falling body a finite proper time.
Giulini's theorem is $R_g \le R$ for every star, so in March the surface of a star covers that sphere.

## Step 6: the drawings

The free chart is drawn through the centre of Giulini's star of uniform density, $\sqrt{N} = \sinh(\omega r)/(\omega r\cosh(\omega R))$ inside and $1 - R_g/r$ outside, with $R_g = R(1 - \tanh(\omega R)/\omega R)$, at $\omega R = 3/2$, where $m = 0.793\,R$, $N = 0.181$ at the centre and $0.364$ on the surface.
`null_rays.e12_star` declares it as a function of $r^2$, in which it is smooth through the centre, with its first two derivatives.
On a plane of $t$ and one length the metric is $N^2(-c^2dt^2 + dx_*^2)$ with $dx_* = dx/N$, so each conformal diagram is the part of Minkowski's that $x_*$ reaches: the whole diamond for the star, Rindler's wedge for the uniform field, and a full diamond for each body, whose two left edges are the sphere where $N$ vanishes, null and singular.
The embedding diagram draws the equator outside the body of March, a flat plane, under Flamm's paraboloid of the same mass.
The star, the uniform field and the body of February are other spacetimes than that body, so their drawings do not mark its moment; `HIDDEN` in `slices.py` says so.
`REGIONS` in `metric_tags.py` counts the free chart alone, since the three solutions have symmetries the general static field lacks.
