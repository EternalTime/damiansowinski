# Spheroidal stars

## Sources and stellar families

Vaidya and Tikekar, doi:10.1007/BF02714870, published in September 1982, interpret the spatial geometry as a spheroid.
Durgapal, doi:10.1088/0305-4470/15/8/039, published in August 1982, chooses a positive integer power for the squared lapse.
Durgapal and Bannerji, doi:10.1103/PhysRevD.27.328, published in January 1983, have an erratum, doi:10.1103/PhysRevD.28.2695, published in November 1983.
Finch and Skea, doi:10.1088/0264-9381/6/4/007, published in April 1989, solve the isotropy equation for a paraboloidal spatial geometry.
All four original records have been verified against Crossref.
Their original full texts remain unread: the available Academy repository address returns 404, the APS text endpoint refuses access, and the IOP PDF returns a bot-management response.
No quotation or equation is attributed to an unread page.

The equations below follow the completely read texts of Delgaty and Lake, arXiv:gr-qc/9809013v1, Thirukkanesh and Maharaj, arXiv:0904.0781v1, and Khugaev, Dadhich, and Molina, arXiv:1603.07118v1.
The first reproduces the named metrics and records their equivalences.
The second derives the neutral field equations and the complete Finch-Skea and Durgapal-Bannerji solutions.
The third derives the spheroidal spatial geometry and identifies the original Vaidya-Tikekar choice as $K=2$ in its convention $g_{rr}=(1+Kr^2/L^2)/(1-r^2/L^2)$.

## Metric potentials and matter

Take $x=Cr^2$, $C>0$, and write

$$ds^2=-y(x)^2c^2dt^2+Z(x)^{-1}dr^2+r^2(d\theta^2+\sin^2\theta\,d\phi^2).$$

With $\mathcal{E}=8\pi G\epsilon/c^4$ and $\mathcal{P}=8\pi Gp/c^4$, where $\epsilon$ is energy density, Thirukkanesh and Maharaj's (4a)-(4c) are

$$\frac{\mathcal{E}}C=\frac{1-Z}x-2Z_x,
\qquad \frac{\mathcal{P}}C=4Z\frac{y_x}y+\frac{Z-1}x,$$

$$4x^2Zy_{xx}+2x^2Z_xy_x+(xZ_x-Z+1)y=0.$$

The last equation enforces equality of radial and tangential pressures.
Specifying $Z$ alone leaves the lapse and pressure undetermined.

### Vaidya-Tikekar and Durgapal-Bannerji

Thirukkanesh and Maharaj's (31), together with their $a=-1/2$, $b=1$, and zero electric field, is

$$Z=\frac{2-x}{2(1+x)},\qquad
y=A(1+x)^{3/2}+B(5+2x)\sqrt{2-x},\qquad
\frac{\mathcal{E}}C=\frac{3(3+x)}{2(1+x)^2}.$$

Both independent lapse solutions are retained.
Khugaev, Dadhich, and Molina's $K=2$ has $x=2r^2/L^2$, so its radial potential is exactly this $Z$.
Substitution of each lapse basis into their pressure-isotropy equation, or into the equation above, gives zero.
Thus the complete isotropic solution spaces agree, including their two integration constants.

### Durgapal III, IV, and V

Write $y=a(1+x)^{n/2}$, with $a>0$ fixing the normalization of time.
Delgaty and Lake reproduce the following radial potentials, with a dimensionless integration constant $k$:

$$Z_{\rm III}=1-\frac{3x}{2(1+x)}\left(1+\frac{k}{\sqrt{1+4x}}\right),\qquad n=3,$$

$$Z_{\rm IV}=\frac{7-10x-x^2}{7(1+x)^2}
+\frac{kx}{(1+x)^2(1+5x)^{2/5}},\qquad n=4,$$

$$Z_{\rm V}=\frac{1-\frac{x}{112}(309+54x+8x^2)
+\frac{kx}{(1+6x)^{1/3}}}{(1+x)^3},\qquad n=5.$$

Their pressure and density are, for each stated $Z$,

$$\frac{\mathcal{P}}C=\frac{2nZ}{1+x}+\frac{Z-1}x,
\qquad \frac{\mathcal{E}}C=\frac{1-Z}x-2Z_x.$$

These formulas retain the radial integration constant.
Taking only the lapse power would omit the rest of each solution.

### Finch-Skea

Thirukkanesh and Maharaj's (30) has $s=\sqrt{1+x}$ and

$$Z=\frac1{1+x},\qquad
y=A(\cos s+s\sin s)+B(\sin s-s\cos s),\qquad
\frac{\mathcal{E}}C=\frac{3+x}{(1+x)^2}.$$

The pressure follows from the common equation above.
In particular, $y_x=(A\cos s+B\sin s)/2$.
An algebraic lapse chosen independently of this equation would generally describe anisotropic matter.

## Centre, surface and exterior

## Equivalences

Delgaty and Lake's Appendix C identifies Durgapal-Bannerji with Buchdahl's 1959 finite star, called Buch1 in their tables.
Their page 6 has exactly the $Z$ and two lapse basis functions stated above, with the constants interchanged.
The $K=2$ Vaidya-Tikekar star belongs to this same solution space.
It is represented once.

The existing `buchdahl_sphere` is Buchdahl's 1964 sphere, called Buch2 by Delgaty and Lake, with isotropic potential $f=GM/(2c^2\sqrt{r^2+R^2})$ and fluid of infinite extent.
It is a different geometry from the finite 1959 star.
The $K=0$ member of the general spheroidal ansatz is the constant-density star already represented by `interior_schwarzschild`.
It needs a relation and an explicit limiting identification, not another named star.

Delgaty and Lake also identify Durgapal III with Heintzmann IIa, and Durgapal II with Wyman IIa at its $n=1$.
The first lapse power is Tolman IV, whose $g_{tt}$ is proportional to $1+Cr^2$.
The separate Durgapal-Pande-Phuloria II solution of 1984 is Tolman VII.
That 1984 result is not Durgapal's 1982 IV or V.
Neither IV nor V nor Finch-Skea has Tolman VII's quadratic density profile.
The `tov` geometry has free mass and lapse functions; each isotropic star here specializes those functions and obeys the equilibrium equation.
## Drawing parameters

## Charts and diagrams

## Validation
