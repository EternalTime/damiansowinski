# Schrödinger spacetime

Anti-de Sitter space in null coordinates with one term in $dt^2$ added,

$$ds^2 = \frac{L^2}{r^2}\left(-\frac{\beta^2}{r^2}c^2dt^2 - 2c\,dt\,d\xi + d\mathbf{x}^2 + dr^2\right),$$

whose isometries are the Schrödinger group of the $d$ dimensions of $\mathbf{x}$.
Its six charts are written by `_tools/derivations/print_charts.py --metric schrodinger_spacetime`, which took two seconds on 2 October 2026, and `verify_metrics.py --system schrodinger_spacetime/<chart>` checks each in under a second.

## Step 1. The line element, its sign and its units

Son's (19), Phys. Rev. D 78, 046003, is $ds^2 = -2(dx^+)^2/z^4 + (-2\,dx^+dx^- + dx^idx^i + dz^2)/z^2$, with the anti-de Sitter radius set to 1.
Balasubramanian and McGreevy's (2.1), Phys. Rev. Lett. 101, 061601, is $L^2(-dt^2/r^{2z} + (d\vec{x}^2 + 2\,d\xi\,dt)/r^2 + dr^2/r^2)$ with the dynamical exponent $z$, so their radial coordinate is Son's $z$ and their $\xi$ is minus his $x^-$.
The charts take Balasubramanian and McGreevy's names, $t$, $\xi$ and $r$, which leaves $z$ free for the exponent, and Son's sign of the cross term, which is also Blau, Hartong and Rollier's.
The coefficient of the first term is written $\beta^2$, as in Herzog, Rangamani and Ross's (2.3), JHEP 11 (2008) 080, Adams, Balasubramanian and McGreevy's (1.1) and Blau, Hartong and Rollier's (1.1) of their second paper: the boost $t \to \lambda t$, $\xi \to \xi/\lambda$ turns $\beta$ into $\lambda\beta$, so it has no invariant meaning, and $\beta = 0$ is anti-de Sitter space.
Son's metric is $\beta^2 = 2$.
$t$ carries dimensions of time, as the time of the nonrelativistic theory, and $\xi$, $x$, $r$, $L$ and $\beta$ are lengths, so the chart is $x^0 = ct$ as everywhere in the collection.

## Step 2. The matter

With $d = n - 3$ spatial directions in $n$ dimensions and $\Lambda = -(d + 1)(d + 2)/2L^2$, the tensor $G_{ab} + \Lambda g_{ab}$ vanishes in every slot but $tt$, where it is $(z - 1)(2z + d)\,h/r^2$ with $h = (\beta/r)^{2z - 2}$.
That is Balasubramanian and McGreevy's pressureless dust, whose density for $d = 3$ they give as $(2z^2 + z - 3)/L^2 = (z - 1)(2z + 3)/L^2$, and Herzog, Rangamani and Ross's null dust $T_{uu} \propto r^{2\nu + 2}$ in the inverse radius.
Son's toy model, his (21) to (23), is a vector field of mass $m^2 = 2(d + 2)$ with $C^- = 1$; for every $z$ it is
$$C = -L\sqrt{\frac{z - 1}{z}}\,\frac{\beta^{z - 1}}{r^z}\,c\,dt,\qquad m^2 = \frac{z(z + d)}{L^2},$$
which is Balasubramanian and McGreevy's $m_A^2 = z(z + d)/L^2$, and at $z = 2$ and $L = 1$ Son's $2(d + 2)$.
`schrodinger_check` holds every chart to Proca's equation $\nabla_aH^{ab} = m^2C^b$ and to $G_{ab} + \Lambda g_{ab} = H_{ac}H_b{}^c - \tfrac{1}{4}g_{ab}H^2 + m^2(C_aC_b - \tfrac{1}{2}g_{ab}C^2)$, slot by slot.
In five dimensions the Einstein tensor itself has no $tt$ component, since $R_{tt} = 10\beta^2/r^4$ and $\tfrac{1}{2}Rg_{tt} = 10\beta^2/r^4$ cancel, and the dust is all in $\Lambda g_{tt}$.

## Step 3. The charts and their sources

- `poincare`: Son's (19) and Balasubramanian and McGreevy's (2.1) at $z = 2$, with one spatial direction, four dimensions, so that it is Siklos's published line element at $H = -\beta^2/x^2$, with his $u$, $v$, $x$ and $y$ its $ct$, $\xi$, $r$ and $x$, and anti-de Sitter space's published Poincaré chart at $\beta = 0$ along $\sqrt{2}\,ct' = ct + \xi$, $\sqrt{2}\,x' = \xi - ct$. `schrodinger_check` holds both.
- `inverse_radius`: Herzog, Rangamani and Ross's (2.1) and (2.3) at $\nu = 1$, $r^2(-2\,du\,dv - \beta^2r^2du^2 + d\mathbf{x}^2) + dr^2/r^2$, with $\rho = L^2/r$ for their $r$, also the form of Maldacena, Martelli and Tachikawa.
- `global`: Blau, Hartong and Rollier's (3.20), JHEP 07 (2009) 027, under their (3.19): $\omega t = \tan\omega T$, $r = R/\cos\omega T$, $x = X/\cos\omega T$ and $\xi = V + (\omega/2c)(R^2 + X^2)\tan\omega T$, with $\omega$ a frequency. The metric gains the one term $-\omega^2(R^2 + X^2)\,dT^2$ in the bracket. It is geodesically complete for every $\omega > 0$, their section 3.3, and the Poincaré chart is its part $|\omega T| < \pi/2$.
- `dynamical_exponent`: Balasubramanian and McGreevy's (2.1) with $z$ free, written with the term $h = (\beta/r)^{2z - 2}$ named, so that every value is a polynomial in $h$ over $r$, $L$ and $z$: each derivative of $h$ is $h$ times a rational function, as with Kiselev's term. $z = 1$ is anti-de Sitter space and $z = 2$ the Poincaré chart.
- `poincare_5d`: the Poincaré chart with two spatial directions, the metric Herzog, Rangamani and Ross, Maldacena, Martelli and Tachikawa, and Adams, Balasubramanian and McGreevy found in type IIB supergravity.
- `poincare_6d`: with three, the case Son proposes for fermions at unitarity in three dimensions of space.

Each chart after the first is checked to be the first pulled back, at $z = 2$, or with every spatial direction but the first left out.
In every chart $\partial_\xi$, or $\partial_V$, is a null Killing vector and a repeated principal null direction of the Weyl tensor, $C_{abcd}k^d = 0$, which in four dimensions is Petrov type N.

## Step 4. What is not a chart

The black holes with Schrödinger asymptotics of the three papers of July 2008 carry a dilaton and are another spacetime.
Blau, Hartong and Rollier's embedding in a flat space of two more dimensions, their second paper's (2.11), is no chart.
The three dimensional case, $d = 0$, has no line of atoms, and its Weyl tensor vanishes identically.

## Step 5. The planes drawn

On a plane of the time and the null coordinate at one depth the metric is $-(L^2/r^2)(2c\,dt\,d\xi + h\,c^2dt^2)$, with $h = \beta^2/r^2$, $\beta^2\rho^2/L^4$, $\beta^2/R^2 + \omega^2R^2/c^2$ on the axis $X = 0$, or $(\beta/r)^{2z - 2}$.
Its null directions are $\partial_\xi$ and $d\xi = -\tfrac{1}{2}h\,c\,dt$, so one edge of every cone lies in the surface of constant time and the cone opens to the half plane as $h \to \infty$, toward the boundary.
$g^{tt} = 0$ in every chart, so no chart's time is a time function, and by Blau, Hartong and Rollier's second paper the spacetime has none: the rows of `null_rays.py` take their future from the timelike Killing vector of the time, `orient="vector"`, as Gödel's do.
The rays along $\xi$ are null geodesics, the lightlike lines of that paper.
Along the other family the radial equation leaves $\ddot{r} = (z - 1)(h/r)\,\dot{t}^2$, which is $(\beta^2/r^3)\,\dot{t}^2$ at $z = 2$, so light launched along it is turned toward larger $r$; in the global chart $\ddot{R} = (\beta^2/R^3 - \omega^2R/c^2)\,\dot{T}^2$, which vanishes at $R^2 = \beta c/\omega$, the bottom of the trap, where that family is geodesic too.
`_tools/test_schrodinger_spacetime.py` holds both accelerations to the published Christoffel symbols.
The views are drawn at $L = \beta = 1$ and $\omega = c/\beta$: $r = \beta/2$, $\beta$ and $2\beta$, the same events in $\rho$, the same depths in $R$, the exponents $z = 1$, $3/2$ and $3$ at $r = \beta/\sqrt{2}$, where $h = 2^{z - 1}$, and $r = \beta$ in five and six dimensions.
At $r = \beta/2$ the exponent $z = 3$ has $h = 16$, a cone too wide for `--verify` to tell its second family from the first across the box, which is why the exponents are drawn a little deeper.

## Step 6. Light in the trap

The surface $X = 0$ of the global chart is totally geodesic, by the reflection of $X$.
At $L = \beta = 1$ and $c = 1$ a null geodesic of the Poincaré chart with momenta $E$ and $P_\xi$ and none along $x$ has $r^2 = r_0^2 + (t - t_0)^2/r_0^2$ with $r_0^2 = P_\xi/2E$ the least depth, from Blau, Hartong and Rollier's (2.2): it never reaches $r = 0$ and reaches $r = \infty$ at a finite affine parameter.
Carried to the global chart at $t_0 = 0$, $\beta = 1$ and $\omega = c$ it is $R^2 = R_0^2\cos^2T + \sin^2T/R_0^2$, which swings between $R_0$ and $1/R_0$ with the period $\pi$ and is constant at $R_0 = 1$.
`trap_rays` in `projections.py` integrates five such rays from $T = 0$ with the published Christoffel symbols, both ways to $|T| = \pi$, and checks each null, of constant momentum along $V$, and on that closed form.

## Step 7. No conformal diagram

Every event on a surface of constant $T$ has the same future, the whole of the spacetime after it, and the same past, Blau, Hartong and Rollier's (3.7), so no map of a surface onto a plane with light at 45° is faithful to it; the spacetime is in `NOT_DRAWN` of `conformal.py`, beside the pp-wave it is conformal to.
No plane of two coordinates through $\partial_\xi$ is a spacetime of its own with a conformal boundary to draw: Balasubramanian and McGreevy note that the boundary is one dimensional.

## Step 8. The plane of $x$ and $r$

The surface $t = 0$, $\xi = 0$ has the metric $(L^2/r^2)(dx^2 + dr^2)$ for every $\beta$, every $z$ and every $\omega$, the hyperbolic plane of curvature $-1/L^2$ as Poincaré's half plane.
It is no surface of revolution in the chart's coordinates, so `Slice` sweeps it about its point $x = 0$, $r = L$: with $\zeta = \tau e^{i\phi}$ on the unit disc, $x + ir = iL(1 + \zeta)/(1 - \zeta)$, and along the half line $x = 0$, $r \ge L$, which is the profile, $\tau = (r - L)/(r + L)$ and the proper distance is $s = L\ln(r/L)$.
The pulled back metric is $L^2dr^2/r^2 + ((r^2 - L^2)/2r)^2d\phi^2$, whose circles have the radius $L\sinh(s/L)$, so the plane is the sheet $Z = L\cosh(s/L) - L$ of a hyperboloid in Minkowski space, as Siklos's wave front is.
The lines $r = L/2$, $L$ and $2L$ are horocycles, each checked to lie on the sheet and to have chords equal to the arcs $L\,dx/r$ of the published metric.
