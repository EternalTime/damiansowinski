# Tilted universes and the whimper

Spatially homogeneous universes whose matter does not move at right angles to the surfaces of homogeneity, on the simplest of them: the locally rotationally symmetric universes of Bianchi's type V, Farnsworth's dust, and the flat model in which Ellis and King drew the whimper.
This note records the sources, the four charts, what each is checked against, and what the drawings draw.

## Step 1: the sources

The candidate list named three papers, and each was verified on Crossref before anything was built.

- A. R. King and G. F. R. Ellis, "Tilted homogeneous cosmological models", Communications in Mathematical Physics 31, 209-242 (1973), doi:10.1007/BF01646266, received 27 November 1972. Read in full from Project Euclid.
- G. F. R. Ellis and A. R. King, "Was the big bang a whimper?", Communications in Mathematical Physics 38, 119-156 (1974), doi:10.1007/BF01651508, received 15 December 1973. Read in full from Project Euclid.
- David L. Farnsworth, "Some New General Relativistic Dust Metrics Possessing Isometries", Journal of Mathematical Physics 8, 2315-2317 (1967), doi:10.1063/1.1705157. The record was verified; the paper is behind a paywall and was not read.

Farnsworth's line element is therefore taken from two papers that quote it and were read: Anton and Clifton, arXiv:2302.05715, their (51) to (56), and Mena, Natário and Tod, arXiv:0707.2519, their (15) to (18), which is the same solution with $\rho = -r$ and a cosmological constant.
Anton and Clifton's is $ds^2 = -dt^2 + X^2(t + Cr)\,dr^2 + e^{-2r}Y^2(t + Cr)\left(dy^2 + dz^2\right)$ with $X = k^{-1}(CY' - Y)$, $Y'^2 = D/3Y + k^2$, and in parametric form $Y = \tfrac{1}{2}Wk(\cosh\eta - 1)$, $t + Cr = \tfrac{1}{2}W(\sinh\eta - \eta)$, $W = D/3k^3$.
A rescaling of $y$ and $z$ sets $k = 1$, and the page writes $X = Y - C\,dY/du$, the opposite sign, which the metric does not see and which is positive where the density is.

The names are as the papers print them: A. R. King, C. B. Collins, L. C. Shepley and S. T. C. Siklos have initials alone, George Ellis is from the address block that closes the paper of 1974, and David L. Farnsworth is Crossref's author line.

What was read for each statement of the history:

- Eliot's line: the last line of "The Hollow Men", first printed whole in Poems 1909-1925 (Faber and Gwyer, London, 1925); the wording was checked against the text as Wikipedia's article on the poem quotes it, and the book against Open Library's record.
- King and Ellis 1973: the abstract, the introduction (the quotations on a fundamental observer), (1.6) for the hyperbolic angle of tilt, Theorem 3.1 (no tilted models of type I), Theorem 4.1 and the paragraph after it (the tilted locally rotationally symmetric models are of type V, and Farnsworth found the dust solution), and the conclusion (Farnsworth's model "particularly simple").
- Ellis and King 1974: the abstract, section 1 (Shepley's example), Theorem 4.1 (a Cauchy horizon in every group type of class B and in none of class A), Lemma 4.4 and Theorem 4.5 (the generators of the horizon are incomplete in the past and end at an intermediate singularity, with the curvature finite in the frame of the fluid), Theorem 5.4 (the stationary region), page 143 ("whimpers"), section 7 (the model on Minkowski's plane, "crosses an infinite number of matter world lines in a finite affine distance", "an open set of models", "not negligible", the infinite blueshift and the "little bang").
- Shepley 1969, Physics Letters A 28, 695: the record and the title on Crossref; what it shows is cited as Ellis and King state it. The letter was not read.
- Collins 1974, read from Project Euclid: the abstract and the introduction ("wishful thinking and dangerous extrapolation", and its first footnote on Taub-NUT space).
- Siklos 1978, read from Project Euclid: the abstract ("two fewer parameters", "not stable features"), the introduction and the close of section 4.
- Secrest and others, arXiv:2009.14826, Krishnan, Mondol and Sheikh-Jabbari, arXiv:2209.14918, and Allahyari and others, arXiv:2307.15791: the abstracts, and of the last its sections 1 and 3.2.

## Step 2: the charts

Four charts, each written by `print_charts.py` and checked by `verify_metrics.py`.

`homogeneous` is Farnsworth's comoving chart with $u = ct + Cr$ for its time: $-(du - C\,dr)^2 + X(u)^2dr^2 + e^{-2r}Y(u)^2(dy^2 + dz^2)$, with $X$ and $Y$ free, so no component assumes a field equation.
The surfaces of constant $u$ are the orbits of $\partial_y$, $\partial_z$ and $\partial_r + y\partial_y + z\partial_z$, whose brackets are those of Bianchi's type V, and the rotation of the plane of $y$ and $z$ is a fourth isometry.
The lines of constant $r$, $y$ and $z$ are geodesics with $u/c$ their proper time, and $g^{uu} = -(X^2 - C^2)/X^2$, so the hyperbolic angle between them and the normal of the surfaces has $\tanh\beta = C/X$, Anton and Clifton's $v = -C/X$.
With $C = 0$ and $X = Y$ the chart is the horn of `small_universes`.

`farnsworth` is the dust solution in the parameter $\eta$: $Y = W\sinh^2(\eta/2)$, $du = Y\,d\eta$, $X = W\sinh^2(\eta/2) - C\coth(\eta/2)$.
Every value is a rational function of $t = \coth(\eta/2)$, $W$ and $C$, and `TiltedParametric` writes each factor in $X$ and $C$ where it is a polynomial in them, as $X$, $X - C$ and $X + C$ are, and otherwise in $\sinh(\eta/2)$ and $\cosh(\eta/2)$.
The check holds it to the homogeneous chart along $du = Y\,d\eta$, to $G_{\mu\nu} = \rho\,u_\mu u_\nu$ with $u_\mu = (-Y, C, 0, 0)$ and $\rho = 3/(WX\sinh^4(\eta/2))$, and to a vanishing Weyl tensor at $C = 0$.
The density is positive where $X > 0$, which is $\eta > \eta_s$, and the Kretschmann scalar is $3(5W^2\sinh^6(\eta/2) + 4C^2\cosh^2(\eta/2))/(W^4X^2\sinh^{14}(\eta/2))$: it diverges on $X = 0$ and is finite on the horizon $X = C$.

`flat_model` is the member with no dust: $dY/du = 1$, so $Y = u + C$ and $X = Y - C = u$ with the origin of $u$ chosen on $X = 0$.
Its Riemann tensor vanishes, and it is Minkowski space along $cT + x = (u + C)e^{-r}$, $cT - x = (u - C)e^{r} + (u + C)e^{-r}(y^2 + z^2)$, $\xi = (u + C)e^{-r}y$, $\zeta = (u + C)e^{-r}z$, the point $u\,n - C\,\partial_r n$ for $n$ the unit hyperboloid in horospheres, mirrored in $x$.
The interval from the origin is $C^2 - u^2$, so the horizon $u = C$ is the light cone of the origin, the lines of matter are straight lines tangent to the hyperboloid at the distance $C$, and on the plane $y = z = 0$ they are $x = C\,\mathrm{sech}\,r - cT\tanh r$.
This is the model of Ellis and King's section 7 with a straight line painted on and the two other dimensions kept; that Farnsworth's solution contains it is this note's own computation, which `tilted_universes_check` holds.

`inertial` is the same space in those inertial coordinates.

## Step 3: the horizon and the whimper

A ray of the plane $y = z = 0$ has $dr/d\eta = Y/(C \pm X)$.
The Killing vector $\partial_r$ gives the first integral $XY\,d\eta/d\lambda$ for both families, so the affine parameter is the integral of $XY\,d\eta$: finite toward $X = 0$ and infinite toward $\eta \to \infty$.
On the horizon the ray of the second family has $\eta$ constant and $\ddot r = X'\dot r^2$ in the homogeneous chart, so $\dot r = 1/(k - X'\lambda)$: it reaches $r \to +\infty$, its past end, at the finite parameter $k/X'$, and runs on for ever toward $r \to -\infty$.
That end is the whimper, Ellis and King's Lemma 4.4.

## Step 4: the drawings

Farnsworth's dust is drawn at $C = W$, where $X = 0$ at $\eta_s = 1.9684$, $u = 0.7707\,W$, and $X = C$ at $\eta_H = 2.3730$, $u = 1.4725\,W$; $Y$ at the singularity is the root of $Y^3 = Y + 1$.
The spacetime diagrams draw the plane of the time and $r$ in each chart, the cones oriented by the velocity of the matter, since under the horizon no coordinate is a time function; the homogeneous view starts above the singularity, and the inertial view draws the lines of matter, the hyperbola they touch and the surfaces of homogeneity.
The rays are checked against the null coordinates $V = (\eta - \eta_H)/(\eta_H - \eta_s)\,e^{(r - G)/a}$ and $r - F_+$ of `tilted_null_tables`, which are regular on the horizon.

The conformal diagrams use $p = \arctan V$ and $q = \arctan e^{-(r - F_+)/a}$, and the flat model's $p, q = \arctan((cT \mp x)/C)$.
Both fill the quadrilateral with corners $(\pi/2, -\pi/2)$, $(-\pi/2, \pi/2)$, $(0, \pi)$ and $(\pi/2, \pi/2)$: the right edge is $X = 0$, since $-VU = 1$ there, the lower left edge is $r \to \infty$, and the horizon is $p = 0$.
It is the shape of Ellis and King's figure 2c.

The embedding diagram is the slice $z = 0$ of the surface of homogeneity $\eta = 3$, where $Y = 4.534\,W$ and $X = 3.429\,W$.
Its metric is $(X^2 - C^2)\,dr^2 + Y^2e^{-2r}dy^2$, a hyperbolic plane of radius $\sqrt{X^2 - C^2} = 3.280\,W$ ruled by horocycles, and no surface of revolution carries such a plane whole, so the strip $0 \le y < 2\pi$ is rolled up, as the horn of `small_universes` is: Beltrami's pseudosphere of that radius, from the rim $r = \ln(Y/\sqrt{X^2 - C^2}) = 0.324$.
In this spacetime $y$ is not periodic, and the caption says that the plane runs on.
The moment is marked on the two drawings of the dust in each kind, and the drawings of the flat model, another spacetime, carry none.
