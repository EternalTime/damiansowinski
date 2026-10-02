# Senovilla's universe with no singularity

Senovilla's universe in his own cylindrical chart is

$$ds^2 = \cosh^4(act)\cosh^2(3a\rho)\left(-c^2dt^2 + d\rho^2\right) + \frac{\cosh^4(act)\sinh^2(3a\rho)}{9a^2\cosh^{2/3}(3a\rho)}\,d\phi^2 + \frac{dz^2}{\cosh^2(act)\cosh^{2/3}(3a\rho)}.$$

The chart is written by `_tools/derivations/print_charts.py --metric senovilla`, and `verify_metrics.py --system senovilla/cylindrical` checks it in about thirteen seconds.

## Step 1. The source of the chart

The line element is (1) of Chinea, Fernández-Jambrina and Senovilla, Phys. Rev. D 45, 481 (1992), which restates Senovilla, Phys. Rev. Lett. 64, 2219 (1990), with $ct$ written for their $t$ and $\rho$ for their $r$.
Senovilla's review of 1998, section 7.6, and his paper of 2007, its appendix, write the same line element in the same coordinates.
It is the one chart the literature uses for this metric, so it is the one chart published.
Chinea, Fernández-Jambrina and Senovilla remark that Cartesian-like coordinates exist which cover the axis, and write none down.
Ruiz and Senovilla's family, (44) of the review, carries a second constant $K$ and is this metric at $K = 1$; it is another spacetime and is told in the History.
`DIMENSIONS` declares $[a] = L^{-1}$.

## Step 2. The cube root

$g_{\phi\phi}$ and $g_{zz}$ carry $\cosh^{-2/3}(3a\rho)$.
The reader takes a function raised to a negative or fractional power, `\cosh^{-2/3}(3a\rho)`, since 2 October 2026.
`norm` writes the hyperbolic functions in exponentials, so the cube root becomes the cube roots of the irreducible factors of $e^{6a\rho} + 1$, $e^{2a\rho} + 1$ and $e^{4a\rho} - e^{2a\rho} + 1$, and of $2$.
Since 2 October 2026 `_canonical` takes a root above the square as a generator $w$ with $w^q = b$: powers of $w$ in a numerator are folded below $q$, and below the line $w$ may stand only as a factor of its own, $1/w = w^{q-1}/b$.
Every component of this chart is one power of the cube root times a rational function, so no sum below a line holds one, and the whole geometry takes about seven seconds.
`_tools/test_roots.py` holds the reader and `norm` to both.

## Step 3. The printing

`senovilla_pretty` finds the power $\cosh^{-k/3}(3a\rho)$, $k = 0$, $1$ or $2$, that a value carries by trying the three, writes what is left in $\sinh$ and $\cosh$ of $act$ and of $3a\rho$, reduces by $\sinh^2 = \cosh^2 - 1$, and factors.
The printer's `arguments` sets those two arguments tight, as the line element writes them, keeps a fractional power on the name, $\cosh^{2/3}(3a\rho)$, and writes $\sinh^k/\cosh^k$ of one argument as $\tanh^k$.
The metric and its inverse are written as the line element writes them.

## Step 4. The checks against the sources

`senovilla_fluid` checks, exactly and before anything is written:

- $G_{\mu\nu} = (\varepsilon + p)u_\mu u_\nu + p\,g_{\mu\nu}$ in every slot, with $u_\mu dx^\mu = -\cosh^2(act)\cosh(3a\rho)\,c\,dt$, $\varepsilon = 15a^2/\cosh^4(act)\cosh^4(3a\rho)$ and $p = \varepsilon/3$ in units of $c^4/8\pi G$, their (2), (5) and (6);
- the fluid's expansion, $3a\sinh(act)/\cosh^3(act)\cosh(3a\rho)$, their (3);
- its acceleration, $a_\rho = 3a\tanh(3a\rho)$, their (4);
- the circular light ray: a circle about the axis is a null geodesic where $a(\cosh^2(3a\rho) - 4)/\sinh(3a\rho)\cosh(3a\rho)$ vanishes, at $\cosh(3a\rho) = 2$, their section 2.

The geodesic equations published are their (7) and (8), with the two that their Killing vectors integrate, (9) and (10), written out, and each is read back against the one printed from the Christoffel symbols.
The Ricci scalar vanishes, since the stress of radiation is trace free, and the Kretschmann scalar is

$$K = \frac{24a^4\left(24\cosh^4(act)\cosh^2(3a\rho) + 33\cosh^4(act) - 32\cosh^2(act)\cosh^2(3a\rho) + 8\cosh^4(3a\rho)\right)}{\cosh^{12}(act)\cosh^8(3a\rho)},$$

which is $792\,a^4$ at the bounce on the axis and smaller everywhere else.

## Step 5. The diagrams

The plane of $t$ and $\rho$ at fixed $\phi$ and $z$ is totally geodesic and conformally flat, so the spacetime diagram's rays run at 45° and the conformal diagram is Minkowski's half diamond by $p, q = \arctan((ct \mp \rho)a)$.
Along a radial light ray $d\rho/d\lambda = h/\cosh^4(act)\cosh^2(3a\rho)$, case 3 of their section 2, so the affine parameter grows without bound toward every far edge, and no edge is singular.
The slice of constant $t$ and $z$ has $g_{\rho\rho} = \cosh^4(act)\cosh^2(3a\rho)$ and circles of radius $\cosh^2(act)\sinh(3a\rho)/3a\cosh^{1/3}(3a\rho)$, whose derivative along $\rho$ is $\cosh^2(act)(2\cosh^2 + 1)/3\cosh^{4/3}$.
It embeds in flat space at every $\rho$, since $3\cosh^{7/3} \ge 2\cosh^2 + 1$ with equality on the axis, and time enters only through the factor $\cosh^2(act)$, so the movie from $act = -1$ to $1$ is one surface shrinking and growing.
The spacetime diagram, the conformal diagram and the embedding diagram are all drawn, each at $a = 1$.
