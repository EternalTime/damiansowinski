# Belinski and Zakharov's gravitational solitons

The spacetime is a family, every vacuum solution the inverse scattering method of Belinski and Zakharov makes by adding solitons to a seed.
One member is written out and drawn: the wave of two solitons on the Kasner universe with $s_1 = s_2 = 1/2$, their (5.10) to (5.12) in Sov. Phys. JETP 48, 985 (1978).
The charts are written by `_tools/derivations/print_charts.py --metric belinski_zakharov`, and `verify_metrics.py --system belinski_zakharov/pole` and `belinski_zakharov/canonical` check them.

## Step 1. Which soliton solution, and why

The 1978 paper works three examples: one soliton on Kasner's universe, (4.17), two solitons on the isotropic Kasner universe, (5.10), and two solitons on flat space, (5.17).
The one soliton solution holds only outside the light cone $z^2 \ge t^2$, joins the seed across it with a jump in its first derivatives, and its hump moves faster than light, so the authors themselves decline to call it a soliton.
The two solitons on flat space are the Kerr-NUT solution continued to an imaginary polar angle, their (5.19), and Kerr's black hole is in the collection already.
The wave of (5.10) is the one the authors describe as two localized disturbances, it is regular everywhere after the singularity, it is in closed form, and its seed, Kasner's universe, and its setting, Gowdy's, are both in the collection.
So it is the member drawn.

## Step 2. The pole coordinates

With $\alpha = ct$ and $\beta = z$ the pole trajectory of a constant $w_1 - iw_2$ solves $\mu^2 + 2(z - w_1 + iw_2)\mu + c^2t^2 = 0$.
Belinski and Zakharov's (5.14) and (5.15) choose the coordinates in which the root is a perfect square: with $w = w_2$ and $w_1 = 0$,

$$ct = w\sinh\tau\cosh\xi, \qquad z = w\cosh\tau\sinh\xi,$$

which is their (5.15) with their new $z$ written $\tau$ and their new $t$ written $\xi$, since the first is the time: $-c^2dt^2 + dz^2 = w^2(\sinh^2\tau + \cosh^2\xi)(-d\tau^2 + d\xi^2)$.
The map is one to one from $\tau > 0$ onto $t > 0$: $Y = \cosh^2\tau$ is the one root above 1 of $w^2Y^2 + (z^2 - w^2 - c^2t^2)Y - z^2 = 0$.
In these coordinates their (5.16) reads $\sigma = \rho^2/\alpha^2 = \tanh^2(\tau/2)$ for the modulus of the pole and $\sin^2\varphi = 1/\cosh^2\xi$ for its phase.

## Step 3. The metric

Substituting (5.16) into (5.10), with $p_1 = \cosh\beta$ and $p_2 = \sinh\beta$ for their two constants bound by $p_1^2 - p_2^2 = 1$, every factor of $(\cosh\tau + 1)^2\cosh^2\xi$ cancels, and with

$$N = \cosh^2\beta\sinh^2\tau + \sinh^2\beta\cosh^2\xi$$

their $Q$ is $4N/((\cosh\tau + 1)^2\cosh^2\xi)$ and the block is $g_{xx}, g_{yy} = \sinh\tau\cosh\xi\,(N + 2 \pm 2\cosh\beta\cosh\tau)/N$ and $g_{xy} = -2\sinh\tau\cosh\xi\sinh\beta\sinh\xi/N$.
The conformal factor as (5.10) prints it, $C_1\alpha^{3/2}\sigma^{-1}Q$, is the coefficient of $-d\tau^2 + d\xi^2$ in these coordinates, and with $C_1 = 1/4$ it is $w^2N/\sqrt{\sinh\tau\cosh\xi}$.
A direct solution of the two equations for the conformal factor confirms the power $3/2$ and that it stands in the pole coordinates: in the chart $\alpha = ct$, $\beta = z$ the factor is that one divided by $w^2(\sinh^2\tau + \cosh^2\xi)$.
The constant $C_1$ is no parameter, since $t, z, w \to \lambda(t, z, w)$ with $x, y \to \lambda^{-1/2}(x, y)$ rescales it, so $w$ carries the one scale and $x$ and $y$ are lengths.
`belinski_zakharov_check` holds the chart to a vanishing Ricci tensor, exactly, to (5.10) as printed at three random points in forty digits, to $\det g_{ab} = \alpha^2$, and to the diagonal metric $\mathrm{diag}(\alpha/\sigma, \alpha\sigma)$ at $\beta = 0$.

## Step 4. What the wave does

On $\xi = \pm\tau$, which is $z = \pm ct$, the light cone of the event $t = 0$, $z = 0$, sit the two pulses.
Far inside the cone, $\tau \to \infty$ at fixed $z$, $N \to \cosh^2\beta\sinh^2\tau$ and the block tends to $(ct/w)\,\mathrm{diag}(1, 1)$, with $f \to \cosh^2\beta\,(ct/w)^{-1/2}$.
Far outside it, $\tau \to 0$ at fixed $t$, $N \to \sinh^2\beta\cosh^2\xi$, the block tends to the same, and $f \to \sinh^2\beta\,(ct/w)^{-1/2}$.
Both are the Kasner universe with $s_1 = s_2 = 1/2$, and the constant in $f$ differs by $\coth^2\beta$ across the pulses.
At $\beta = 0$ the metric is diagonal and $f$ falls as $1/z^2$ outside the cone.
On $\xi = 0$ the block is diagonal with $g_{xx}/g_{yy} = ((\cosh\beta\cosh\tau + 1)/(\cosh\beta\cosh\tau - 1))^2$, which is what the ring of particles in the embedding diagram shows.

## Step 5. The canonical chart

Belinski and Zakharov's block form (1.1), with $\alpha = ct$, is $f(-c^2dt^2 + dz^2) + g_{ab}dx^adx^b$ with $\det g = c^2t^2$.
The block is written as Gowdy's is, $g_{ab}dx^adx^b = (ct/w)(e^P(dx + Q\,dy)^2 + e^{-P}dy^2)$, which has that determinant for every $P$ and $Q$, with $w$ a length so that $x$ and $y$ are lengths.
The wave's own $f$, $P$ and $Q$ in $t$ and $z$ hold a root inside a root, $\cosh\tau = \sqrt{(S + c^2t^2 - z^2 + w^2)/2w^2}$ with $S^2 = (c^2t^2 - z^2 - w^2)^2 + 4w^2c^2t^2$, so the chart leaves the three functions free, as Gowdy's and black Saturn's charts do, and no component assumes a field equation.
The parameters state the vacuum equations, their (1.3) to (1.5) in $P$ and $Q$: two wave equations and the two first derivatives of $\ln f$.
`belinski_zakharov_check` replaces $\partial_t^2P$, $\partial_t^2Q$ and every derivative of $f$ by them and holds each Ricci component to vanishing, holds the wave's $f = N/((\sinh^2\tau + \cosh^2\xi)\sqrt{ct/w})$, $e^P = (N + 2 + 2\cosh\beta\cosh\tau)/N$ and $Q = -2\sinh\beta\sinh\xi/(N + 2 + 2\cosh\beta\cosh\tau)$ to solving them at three random points in forty digits, and holds the chart at those functions to being the pole chart carried along the map.
`OVERRULED` in `metric_tags.py` sets the tag `vacuum`, which the free functions cannot show.

## Step 6. The drawings

Every drawing takes $w = 1$ and $\cosh\beta = 5/4$, $\sinh\beta = 3/4$, $\beta = \ln 2$.
Both planes of the time and the direction of travel are conformally flat, so their rays are at 45° and `--verify` checks them against $\tau \pm \xi$ and $ct \pm z$; the canonical plane's row declares the wave's functions and is checked against the published Einstein tensor.
The conformal diagram is the upper half of Minkowski's diamond, $p, q = \arctan((ct \mp z)/w)$, with the singularity along the bottom.
The embedding diagram is a ring of free particles at rest on $\xi = 0$, where $\Gamma^\xi{}_{\tau\tau}$ vanishes by the symmetry $\xi \to -\xi$, $y \to -y$, drawn as its world tube and as a movie.
`_tools/test_belinski_zakharov.py` holds the files to all of this.
