# The Hotta-Tanaka shock wave

The five charts of `hotta_tanaka.json` are written by `print_charts.py --metric hotta_tanaka`, which computes every tensor from each line element with the checker's own `Geometry` and reads every printed value back before it writes.
This note derives what the charts are, how a delta on a curved background is read, and what the diagrams draw.

## Step 1. The five-dimensional form

De Sitter space of radius $a$ is the hyperboloid $-Z_0^2 + Z_1^2 + Z_2^2 + Z_3^2 + Z_4^2 = a^2$ in a flat space of five dimensions, with $\Lambda = 3/a^2$.
Hotta and Tanaka write the Schwarzschild-de Sitter metric on it to first order in the mass $m$, boost it along $Z_1$ with speed $v$, and let $v \to 1$ with $p = m/\sqrt{1 - v^2}$ fixed.
Their equation (18) is the limit,

$$ds^2 = ds^2_{dS} + \frac{4GE}{c^4}\left(z\ln\frac{1 + z}{1 - z} - 2\right)\delta(Z_0 + Z_1)\left(dZ_0 + dZ_1\right)^2, \qquad z = \frac{Z_4}{a},$$

with $GE/c^4$ for their $p$.
The shock is on the null cone $Z_0 + Z_1 = 0$, a cosmological horizon, where $Z_2^2 + Z_3^2 + Z_4^2 = a^2$: a sphere of radius $a$ at every $Z_0 - Z_1$.
The profile diverges at $z = \pm 1$, the two null particles, each of energy $E$ by their (43).
Every chart is this form pulled back, and `hotta_tanaka_check` holds each to it: the chart without the shock is the hyperboloid's, and the shock's term is the profile at $z = Z_4/a$ on the shock times $\delta(Z_0 + Z_1)(dZ_0 + dZ_1)^2$.

A boost along $Z_1$ multiplies $Z_0 + Z_1$ by a number $\kappa$ and the shock's term by the same $\kappa$, so $E$ is the energy in the frame of the five coordinates, the one in which the black hole was at rest before the boost.

## Step 2. The charts

`hotta_tanaka_embedding` holds each chart's map into the hyperboloid.

The conformally flat chart is Hotta and Tanaka's (22) to (27), with their conformal time written $-\eta$ so that $\eta$ increases to the future: $Z_0 + Z_1 = (\eta^2 - \rho^2)/2\eta$, which has unit slope across the shock $\rho = -\eta$, and $Z_4 = a\cos\theta$ on it.
Its origin $\rho = 0$ is the free observer whose event horizon the shock is, and the shock is the sphere of comoving radius $-\eta$ about that observer, of proper radius $a$.

The global chart is Podolský and Griffiths's (10) and (11): $Z_0 = -a\cot\eta$ and $Z_1 = a\cos\chi/\sin\eta$, so $Z_0 + Z_1 = a(\cos\chi - \cos\eta)/\sin\eta$, with slope $a$ across the shock $\chi = \eta$, which is why its term carries $4GEa/c^4$.

The Kruskal chart is the form of Dray and 't Hooft that Hotta and Tanaka use in their (29) to (36), with $v$ turned to increase to the future and both null coordinates scaled alike: $Z_0 \pm Z_1 = 2a^2(u, v)/(a^2 - uv)$ and $R = a(a^2 + uv)/(a^2 - uv)$.
$Z_0 + Z_1$ has slope 2 across $u = 0$, so the shock's term is $(8GE/c^4)(\dots)\delta(u)\,du^2$.
In the global chart's coordinates $u = a\tan((\eta - \chi)/2)$ and $v = -a\cot((\eta + \chi)/2)$.
It covers the whole hyperboloid, $-a^2 < uv < a^2$, with the poles on $uv = -a^2$ and infinity on $uv = a^2$.

The null cylindrical chart is Podolský and Griffiths's unified form of 1999, their (7), in Podolský and Ortaggio's parametrisation (8): $Z_0 \pm Z_1 = (u, v)/\Omega$, $Z_2 + iZ_3 = \rho e^{i\phi}/\Omega$ and $Z_4 = a(2/\Omega - 1)$ with $\Omega = 1 + (\rho^2 - uv)/4a^2$.
On the shock $z = (4a^2 - \rho^2)/(4a^2 + \rho^2)$, so $(1 + z)/(1 - z) = 4a^2/\rho^2$ and the profile inside the conformal bracket is $-(4GE/c^4)\left((1 - \rho^2/4a^2)\ln(\rho^2/4a^2) + 2 + \rho^2/2a^2\right)$.
One particle is at $\rho = 0$ and the other at $\rho \to \infty$.
With $a \to \infty$ at a fixed $\rho$ the chart is Aichelburg and Sexl's null cylindrical chart with $\rho_0 = 2a/e$, Hotta and Tanaka's (21).
The Cartesian form of the same chart, in $x = \rho\cos\phi$ and $y = \rho\sin\phi$, is not printed: its Riemann tensor is the same sum of the de Sitter part and the wave written over one denominator in $x$ and $y$, some sixteen terms a component.

The Kundt chart is Podolský and Griffiths's (18) and (19) with $w = a^2/t$ and $u = \rho - t$ for their $t$ and $\rho$.
Its line element holds $u$ only in the delta, so the delta may be any profile of $u$ and every tensor is linear in it: the chart is the impulsive member of the Kundt class $KN(\Lambda)$ of Ozsváth, Robinson, and Rózga, and it is the one chart whose delta needs no rule.
It covers the half $\cos\phi > 0$ of each wave front.

## Step 3. A delta on a curved background

In every chart but Kundt's the coefficient of the delta, or what multiplies it in the inverse metric, varies across the shock, and the curvature holds products of two deltas.
A product of deltas is no distribution, so the checker reads the delta of a system listed in `IMPULSES` as a smooth pulse $d(x)$ while the tensors are built, `Pulse` in `verify_metrics.py`, and takes each value to the limit of a narrow pulse once it stands whole, `on_the_shock`:

- a term with one factor is a distribution, and $f(x)\,d^{(k)}(x) = \sum_j (-1)^j\binom{k}{j}f^{(j)}(0)\,d^{(k-j)}(x)$ exactly;
- a product of $N$ factors of total weight $W$, each $d^{(k)}$ counting $k + 1$, times $x^n$ scales as the width of the pulse to the power $n + 1 - W$: it goes to zero for $n \ge W$, and for $n = W - 1$ to a multiple of $d(x)$ whose factor is the integral of $s^n$ times the pulses, which for two factors of a pulse even about the shock is odd and vanishes;
- any other product is left standing, and the comparison it enters fails.

In the four charts every product that arises is $d^2$ times a function that vanishes on the shock, $u\,K^2d^2$ in $\Gamma^v{}_{uu}$ and its like in the Riemann tensor, so every value has a limit, and the limit needs only that the pulse be even.
The convention on the page says so.
The limit is taken at the end and never between two steps, since the limit of a product is not the product of the limits: `Reduced(raw=True)` in `chart_printer.py` hands the unreduced tensor back in where an index is raised.

The argument of a delta may be a sum of two coordinates with a unit coefficient on one, as $\eta - \chi$; the first such coordinate is written as the argument less the rest, so a coefficient on the shock is written in $\chi$ in the global chart and in $\rho$ in the conformally flat one.

## Step 4. What is checked

`hotta_tanaka_check` holds each chart, before it is written, to:

- the hyperboloid, and the five-dimensional form of Step 1 for the shock's term;
- $R_{\mu\nu} = (3/a^2)g_{\mu\nu}$, the delta of $g_{uu}$ with it, so the shock adds no source off the particles;
- the Weyl tensor Gauss's equation makes of the five-dimensional wave.

The last is independent of Step 3.
The five-dimensional metric is a pp-wave, whose only curvature is $R_{UpUq} = -\tfrac12\partial_p\partial_qH\,\delta(U)$, and the hyperboloid's extrinsic curvature is $(g_{ab} + \tfrac12\mathcal{G}\,\delta(U)\,\ell_a\ell_b)/a$ with $\mathcal{G} = Z_4H' - H$ and $\ell = dU$, as Podolský and Ortaggio's $T\cdot\nabla_TN$ gives it.
Gauss's equation then makes the Weyl tensor

$$C_{abcd} = \left[e^A_ae^B_be^C_ce^D_dR^{(5)}_{ABCD} + \frac{\mathcal{G}}{2a^2}\left(g_{ac}\ell_b\ell_d + g_{bd}\ell_a\ell_c - g_{ad}\ell_b\ell_c - g_{bc}\ell_a\ell_d\right)\right]\delta(U),$$

and the Riemann tensor is that plus $(g_{ac}g_{bd} - g_{ad}g_{bc})/a^2$ with the whole metric, shock included.
Each chart's Weyl tensor, computed through the pulse, is compared with this at three random points in forty digits, since the profile's logarithm is spelled in more than one way.
The Kretschmann scalar is $24/a^4$ exactly, the wave being of Petrov type N.

## Step 5. The jump of a ray

On a plane of the two null directions at a fixed place on the wave front the rays of one family run beside the shock, and each ray of the other crosses it.
In the Kruskal chart $-4a^4du\,dv/(a^2 - uv)^2 + (8GE/c^4)P\,\delta(u)\,du^2 = 0$ with $P = \cos\theta\ln\frac{1 + \cos\theta}{1 - \cos\theta} - 2$ gives $dv/du = (2GE/c^4)P\,\delta(u)$ at $u = 0$, so

$$\Delta v = \frac{2GE}{c^4}P(\theta),$$

the same at every $v$, and equal to $-\Gamma^v{}_{uu}$'s coefficient of $\delta'(u)$ integrated twice.
$P = -2$ on the equator, where the ray comes out at a smaller $v$, advanced, and $P$ is positive within $33.53°$ of either particle, where it is delayed; Sfetsos prints the angle as $33.52°$.
In the other charts the jump is the same statement: $-a\cot((\eta + \chi)/2)$ and $-2a^2/(\eta - \rho)$ are the Kruskal $v$, the null cylindrical chart's $v$ jumps by its own profile, and the Kundt chart's $w$ by $(2GE/c^4)P/\sin\theta$ on $\phi = 0$.
A linear reading, $\Delta(\eta + \chi) = (4GE/c^4a)\sin^2\chi\,P$, holds only for a small jump, since $\sin^2\eta$ changes across it.

## Step 6. The spacetime diagrams

Two planes of each chart at $8GE/c^4 = a$: through the equator of the wave front, $P = -2$, and $15°$ from a particle, $P = 1.917$; the null cylindrical chart's are $\rho = 2a$, the equator, and $\rho = a/4$.
`null_rays.py` integrates the rays, so each view declares the pulse $e^{-s^2/w^2}/(w\sqrt\pi)$, $w = 0.05$, in place of the delta.
Inside the pulse $g^{rr}$ takes the pulse's sign, which marks nothing of the spacetime, so a view with a declared pulse marks no zero of $g^{rr}$.
Each box is placed so that no cone of the lattice stands within the pulse.
The equatorial plane is totally geodesic by the reflection $Z_4 \to -Z_4$; on the other, $\Gamma^\theta{}_{uu}$ turns a crossing ray toward the particle, and the captions say so.

## Step 7. The conformal diagram

On the plane $\theta = \pi/2$, with $8GE/c^4 = a = 1$, $p = \arctan u - \pi/4$ and $q = \arctan(v - \Delta v\,\Theta(u)) + \pi/4$ with $\Delta v = -1/2$.
Ahead of the shock $X = q - p = \chi$ and $T = p + q = \eta - \pi/2$, de Sitter's square under its diagonal.
Behind it each crossing ray keeps its $q$, so future infinity $uv = 1$ is $T = \arctan u + \arctan(1/u + \tfrac12)$, which rises from $\pi/2$ at the corner the shock ends on to $\pi/2 + \arctan\tfrac12$, and the pole $\chi = 0$ is $X = \pi/2 + \arctan(\tfrac12 - 1/u) - \arctan u$.
A ray that leaves the pole $\chi = \pi$ at $\eta_s$ has $v = \tan(\eta_s/2)$ ahead of the shock and $v - \tfrac12$ behind it, and reaches the pole $\chi = 0$, where $v < 0$, when $\eta_s < 2\arctan\tfrac12$.
In de Sitter space no such ray arrives, which is Gao and Wald's statement that the past of an observer holds no Cauchy surface, and the taller drawing is Leblond, Marolf, and Myers's.
`conformal.py` checks both maps null against the published metric of the Kruskal and the global charts, the jump against $\Gamma^v{}_{uu}$, the plane totally geodesic in both, and the two statements above.

## Step 8. The embedding diagram

A sphere of constant $u$ and $v$ of the Kruskal chart is round, of radius $R = a(a^2 + uv)/(a^2 - uv)$.
The dust at rest in the closed slicing on the great sphere $\chi = \pi/2$ has $u = v = a\tanh(\tau/2a)$ and meets the shock all at once, at $\tau = 0$, where its sphere is the wave front.
The ring $\theta_0 = 30°$ of that dust is run from $\tau = -0.6\,a$ with the published Christoffel symbols and a pulse of width $0.005$, and each frame of the movie is the sphere through it.
By Podolský and Ortaggio's (37) and (43), with $s = \sinh(\tau/a)$ and $g = 4GE/c^4a$, behind the shock

$$Z = a\sin\theta_0\sqrt{1 + s^2} - \frac{g\,a\,s}{\sin\theta_0}, \qquad Z_4 = a\cos\theta_0\sqrt{1 + s^2} + g\,a\,s\ln\cot\frac{\theta_0}{2},$$

so the ring closes on the axis at $s = \sin^2\theta_0/\sqrt{g^2 - \sin^4\theta_0}$, which is $1/\sqrt3$ here, and a ring with $\sin^2\theta_0 > g$ never does.
The kicks follow from the five-dimensional geodesic equation: $\Delta(dZ_p/dU) = \tfrac12(H_{,p} - Z_p\mathcal{G}/a^2)$ with $\mathcal{G} = (8GE/c^4)/\sin^2\theta_0$.
The script checks the run against these to twice the width of the pulse, the error falling in proportion to the width.

## Step 9. The moments on the other drawings

Each key moment of the movie is the sphere $Z_0 + Z_1 = U$, $Z_0 - Z_1 = V$ with $U = a\sinh(\tau/a)$ and $V = (R^2 - a^2)/U$, or $V = U$ up to the shock.
On a plane of a chart it is the event with those two values, `slices.hotta_tanaka_event`, at $a = 1$:

- Kruskal: $u, v = (U, V)(1 - x)/2$ with $x = (\sqrt{1 + UV} - 1)/(\sqrt{1 + UV} + 1)$; global: $\eta, \chi = \pi/2 + \arctan v \pm \arctan u$;
- null cylindrical: $u, v = (U, V)\Omega$ with $UV\Omega^2/4 + \Omega = 1 + \rho^2/4$;
- conformally flat: $\eta^2 - \rho^2 = 2\eta U$ and $\rho\cos\theta - 1 = \eta(V - U)/2$, a quadratic for $\eta$;
- Kundt: $w = (V - U + 2\cos\theta)/2\sin\theta$ and $wu^2 + 2u = 2U/\sin\theta$.

Up to the shock $V = U$, and then the conformally flat chart's plane $\theta = \pi/2$ has no solution with $\eta$ finite and the Kundt chart's has $w = 0$, its edge; those views are in `HIDDEN`.
The conformal diagram marks each moment as a point by its own map.
