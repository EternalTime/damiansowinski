# Pre-big bang dilaton vacuum

## Coordinates and sources

Maurizio Gasperini and Gabriele Veneziano, Phys. Rep. 373, 1 (2003), arXiv:hep-th/0207130, equation (3.1) with its exponents constant, give the synchronous chart, $ds^2 = -dt^2 + \sum_i(-t)^{2\beta_i}dx_i^2$ with $\phi = (\sum_i\beta_i - 1)\ln(-t)$ and $\sum_i\beta_i^2 = 1$ on $t < 0$.
It is the four-dimensional case of Veneziano's equation (4.9) in CERN-TH-6077/91, which writes the shifted dilaton as $\bar\phi = -\ln(t/t_0)$ with the sign of $-t$ dropped from the argument.
The field equations are Veneziano's (3.7) with no potential, $R_{\mu\nu} + \nabla_\mu\nabla_\nu\phi = 0$ and $R + 2\Box\phi - (\nabla\phi)^2 = 0$, whose homogeneous form is Gasperini and Veneziano's (2.10) to (2.12).
With the dilaton above, the space components of the first vanish for any exponents, its $tt$ component is $-(\sum_i\beta_i^2 - 1)/t^2$, and the second is $(\sum_i\beta_i^2 - 1)/t^2$, so both hold exactly on the sphere $\sum_i\beta_i^2 = 1$.
On the sphere the Ricci tensor is $R_{\mu\nu} = -\nabla_\mu\nabla_\nu\phi$, so $R_{tt} = (\sum_i\beta_i - 1)/t^2$, $R_{ii} = \beta_i(\sum_j\beta_j - 1)(-t)^{2\beta_i - 2}$, and $R = (\sum_i\beta_i - 1)^2/t^2$.
The Kretschmann scalar is $K = 4[\sum_i\beta_i^2(\beta_i - 1)^2 + \sum_{i<j}\beta_i^2\beta_j^2]/t^4$ for any exponents, the sum of the squares of $R^{ti}{}_{ti} = \beta_i(\beta_i - 1)/t^2$ and $R^{ij}{}_{ij} = \beta_i\beta_j/t^2$.
It vanishes only where each exponent is zero or one and at most one is nonzero, which on the sphere is $(1, 0, 0)$ and its permutations, the flat members; everywhere else it grows as $1/t^4$ toward $t = 0$.
The Weyl components are written in their shortest form modulo $\sum_i\beta_i^2 = 1$, eight terms in place of ten.
The values on the sphere are compared through its stereographic projection from $(0, 0, 1)$, $\beta = (2u, 2v, u^2 + v^2 - 1)/(1 + u^2 + v^2)$, which reaches every point but that flat member.

The conformal chart takes the two members with $\beta_1 = \beta_2 = \beta_3 = \beta$, where the sphere gives $\beta = \mp 1/\sqrt3$.
Veneziano's (2.3) defines the conformal time by $a\,d\eta = dt$, and Gasperini and Veneziano's (2.94) give $a_\pm \sim |\eta|^{\pm 1/(\sqrt3 \mp 1)}$ for three dimensions of space.
With $a = (-t)^\beta$ the time is $\eta = -(-t)^{1 - \beta}/(1 - \beta)$, so $-t = (-\eta/(1 + k))^{1 + k}$ and $a = (-\eta/(1 + k))^k$ with $k = \beta/(1 - \beta)$.
So $k = -1/(1 + \sqrt3)$ is $a_-$, the expanding member, and $k = 1/(\sqrt3 - 1)$ is $a_+$, the contracting one.
The two values are the roots of $2k^2 - 2k - 1$, which is $3\beta^2 = 1$ written in $k$, and the dilaton is $\phi = (3\beta - 1)\ln(-t) = (2k - 1)\ln(-\eta/(1 + k))$.
Both field equations in this chart are proportional to $2k^2 - 2k - 1$, and its curvature holds for every $k$, so none of its components assumes a value of $k$.
The constant $1 + k$ keeps $x$, $y$, and $z$ the coordinates of the synchronous chart and puts $t = -1$ at $\eta = -(1 + k)$, where $a = 1$.
The two charts' Ricci scalars agree: at $k = -1/(1 + \sqrt3)$, $6k(k - 1)/(\eta^2a^2) = (1 + \sqrt3)^2/t^2$, since $\eta^2a^2 = (1 + k)^2t^2$.

The Einstein frame metric $\tilde g_{\mu\nu} = e^{-\phi}g_{\mu\nu}$ is a second metric, conformal to this one, Gasperini and Veneziano's (1.8) to (1.14) and (3.51).
With $B = \sum_i\beta_i$, its time is $-\tilde t \propto (-t)^{(3 - B)/2}$ and its scale factors are $(-\tilde t)^{p_i}$ with $p_i = (1 - B + 2\beta_i)/(3 - B)$, so $\sum_i p_i = 1$ and $\sum_i p_i^2 = 1 - 2(1 - B)^2/(3 - B)^2$.
That is Kasner's universe with a massless scalar field, its $q^2 = 2(1 - B)^2/(3 - B)^2$, run toward $\tilde t = 0$ from below; the isotropic expanding member has $p_i = 1/3$, the massless scalar universe $\tilde a = (-\tilde t)^{1/3}$.
It appears in the History and in the relation to Kasner's universe with a scalar field, and no chart of this spacetime carries it, since its line element is another metric.

## Drawings

Every drawing is of the expanding isotropic member, $\beta_1 = \beta_2 = \beta_3 = -1/\sqrt3$, which is $k = -1/(1 + \sqrt3)$ in conformal time.
The plane $y = z = 0$ is totally geodesic, since the reflections $y \to -y$ and $z \to -z$ fix it, and its light rays see only $\beta_1$.
Along $x$ a ray keeps $x \mp (-t)^s/s$ with $s = 1 - \beta_1 = 1 + 1/\sqrt3$, which is $x \pm \eta$ in conformal time, so the cones close up as $t \to 0$, where the scale factor diverges.
Both spacetime diagrams span $x \in [-2.2, 2.2]$ and a time from $-4.4$ to $0$, the curvature singularity along the top edge.
The marked moment is $t = -1$, where every scale factor is one, which is $\eta = -1/s = -(1 + k)$ in conformal time.
The embedding turns the plane $z = 0$ at that moment about its origin: it is Euclid's plane, $dx^2 + dy^2$, out to the radius $2$ of the slices.
The conformal diagram compactifies the plane $y = z = 0$ in conformal time with $p, q = \arctan(\eta \mp x)$, the half of Minkowski's diamond below $T = 0$: the singularity $\eta = 0$ is its top edge, where the Kretschmann scalar grows as $(-\eta)^{-4(1 + k)}$, and its lower edges are past null infinity, where the scalar falls to zero.
The curvature figure reads the stored scalars of the synchronous chart, $K = (20 + 8\sqrt3)/(3t^4)$ and $R = (1 + \sqrt3)^2/t^2$ at this member, from $t = -4.4$ to $t = -0.2$.
All four drawing kinds are present.

The collection's `view_text.length` counts 223 characters on the synchronous spacetime view, 156 on the conformal time one, 246 on the conformal diagram, 212 on the embedding, and 163 on the curvature view.
