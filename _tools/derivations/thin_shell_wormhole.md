# Visser's thin shell wormhole

The two charts of `thin_shell_wormhole.json` are written by `print_charts.py --metric thin_shell_wormhole`, which computes every tensor from each line element with the checker's own `Geometry` and reads every printed value back before it writes.
This note derives what the charts are, why their tensors come out as printed, and what the diagrams draw.

## Step 1. The spacetime

Matt Visser takes two copies of Schwarzschild's spacetime, removes from each the region $r \le a$ with $a > r_s$, and identifies the two spheres $r = a$.
Each side is a piece of Schwarzschild's exterior,

$$ds^2 = -\left(1 - \frac{r_s}{r}\right)c^2dt^2 + \frac{dr^2}{1 - r_s/r} + r^2\,d\Omega^2, \qquad r \ge a,$$

so off the throat every tensor is Schwarzschild's, and the Schwarzschild chart publishes them for either side.
The throat is a timelike hypersurface, and $a > r_s$ keeps every point of it outside the horizon, so $\partial_t$ is timelike everywhere and there is no horizon at all.

## Step 2. The chart through the throat

With $\ell = r - a$ on one side and $a - r$ on the other, $r = a + |\ell|$ and $dr^2 = d\ell^2$ on both, so

$$ds^2 = -\left(1 - \frac{r_s}{a + |\ell|}\right)c^2dt^2 + \frac{d\ell^2}{1 - r_s/(a + |\ell|)} + (a + |\ell|)^2\,d\Omega^2, \qquad \ell \in \mathbb{R}.$$

The metric is continuous at $\ell = 0$ and its first derivative jumps there: $\partial_\ell r = \operatorname{sgn}\ell$ and $\partial_\ell^2 r = 2\delta(\ell)$.
`thin_shell_pullback` in `print_charts.py` checks that the chart is Schwarzschild's with $r = a + \ell$ on the side $\operatorname{sgn}\ell = 1$ and with $r = a - \ell$ on the other.
$\ell$ is not the proper distance from the throat, which is $\int d\ell/\sqrt{1 - r_s/r}$ and has no inverse in closed form; $\ell$ is the one coordinate through the throat whose metric is algebraic, so the checker can compare every component exactly.

## Step 3. How the checker reads the kink

The reader writes $|\ell|$ as sympy's `Abs` and $\operatorname{sgn}\ell$ as `sign`, and sympy differentiates them to $\operatorname{sgn}\ell$ and $2\delta(\ell)$.
`_kinked` in `verify_metrics.py` writes $|\ell| = \ell\operatorname{sgn}\ell$, so every value is a rational function of $\ell$ and $\operatorname{sgn}\ell$, and `_canonical` folds $\operatorname{sgn}^2\ell = 1$ exactly as it folds a square root into its radicand.
That form is canonical: a numerator $A + B\operatorname{sgn}\ell$ over a denominator free of $\operatorname{sgn}\ell$ vanishes for every $\ell$ exactly when it vanishes on both branches $\operatorname{sgn}\ell = \pm 1$ as rational functions, which forces $A = B = 0$.
A delta multiplies a function continuous at $\ell = 0$ by that function's value there, so the coefficient of $\delta(\ell)$ is evaluated at $\ell = 0$, after checking that its odd part in $\operatorname{sgn}\ell$ vanishes there, since $\delta(\ell)\operatorname{sgn}\ell$ has no value.
Every Christoffel symbol below carries $\operatorname{sgn}\ell$ to the first power and no delta, every curvature component is linear in its second derivatives, so no product of a delta with $\operatorname{sgn}\ell$ arises, and the check never stops on one.
`KinkForms` in `print_charts.py` prints the canonical form back in $|\ell|$ and $\operatorname{sgn}\ell$: it puts $\ell = |\ell|\operatorname{sgn}\ell$, folds the square of $\operatorname{sgn}\ell$ again, multiplies out the conjugate $(a + |\ell|)(a - |\ell|)$ that canonical form leaves in a denominator, and factors the even and odd parts back into powers of $a + |\ell|$.

## Step 4. The Christoffel symbols

With $f = 1 - r_s/r$ and $\partial_\ell f = (r_s/r^2)\operatorname{sgn}\ell$ they are Schwarzschild's times $\operatorname{sgn}\ell$ wherever an index is $\ell$ an odd number of times:

$$\Gamma^\ell{}_{tt} = \frac{r_s(r - r_s)}{2r^3}\operatorname{sgn}\ell, \quad \Gamma^t{}_{t\ell} = -\Gamma^\ell{}_{\ell\ell} = \frac{r_s\operatorname{sgn}\ell}{2r(r - r_s)}, \quad \Gamma^\ell{}_{\theta\theta} = -(r - r_s)\operatorname{sgn}\ell, \quad \Gamma^\theta{}_{\ell\theta} = \frac{\operatorname{sgn}\ell}{r},$$

with $r = a + |\ell|$, and the angular symbols unchanged.
They jump at the throat and carry no delta, so every geodesic crosses it with its velocity continuous and its acceleration jumping: a radial light ray keeps $d\ell/d(ct) = \pm f$, whose derivative along $\ell$ changes sign at $\ell = 0$.

## Step 5. The curvature

Each Riemann component is Schwarzschild's off the throat plus a delta from $\partial_\ell\Gamma$, since $\partial_\ell\operatorname{sgn}\ell = 2\delta(\ell)$ and every coefficient of it is evaluated at $r = a$.
For $\Gamma^\theta{}_{\ell\theta} = \operatorname{sgn}\ell/r$ the delta is $2\delta(\ell)/a$, which is why

$$R^\theta{}_{\ell\ell\theta} = \frac{r_s}{2r^2(r - r_s)} + \frac{2\delta(\ell)}{a},$$

and the others follow the same way.
The Ricci and Einstein tensors vanish off the throat, as Schwarzschild's do, and are pure deltas on it:

$$G^t{}_t = \frac{4(a - r_s)}{a^2}\delta(\ell), \qquad G^\theta{}_\theta = G^\phi{}_\phi = \frac{2a - r_s}{a^2}\delta(\ell), \qquad G^\ell{}_\ell = 0, \qquad R = -\frac{2(4a - 3r_s)}{a^2}\delta(\ell).$$

## Step 6. The shell

The proper distance $\eta$ across the throat has $d\eta = d\ell/\sqrt{f(a)}$ there, so $\delta(\ell) = \delta(\eta)/\sqrt{f(a)}$.
A shell of surface energy density $\sigma$ and surface pressure $p$ has $T^\mu{}_\nu = \operatorname{diag}(-\sigma, 0, p, p)\,\delta(\eta)$, and $G^\mu{}_\nu = (8\pi G/c^4)T^\mu{}_\nu$ gives

$$\sigma = -\frac{c^4}{2\pi Ga}\sqrt{1 - \frac{r_s}{a}}, \qquad p = \frac{c^4}{4\pi Ga}\,\frac{1 - r_s/2a}{\sqrt{1 - r_s/a}},$$

Visser's equation (4.3) and Poisson and Visser's equations (18) and (19) with $M = r_s/2$ and $G = c = 1$.
$G^\ell{}_\ell = 0$ says the shell carries no stress across itself, as a thin shell must, and $\sigma < 0$ for every $a > r_s$.
The null energy condition on the shell, $\sigma + p \ge 0$, reads $(3r_s - 2a)/2a \ge 0$ times a positive factor, so it holds exactly when the throat is at or inside the photon sphere $a = 3r_s/2$, as Cardoso, Franzin and Pani state.

## Step 7. The Weyl tensor

The Weyl tensor carries a delta as well, and every delta in it is a multiple of $(2a - 3r_s)\delta(\ell)$, for instance $C^\theta{}_{\phi\theta\phi} = r_s\sin^2\theta/r + (2a - 3r_s)\sin^2\theta\,\delta(\ell)/3$.
At $a = 3r_s/2$ it vanishes: there $\sigma + p = 0$, so the surface stress is $-\sigma$ times the metric of the throat, the same for every observer on the shell, and a stress with no preferred direction along the shell has no traceless part to put into the Weyl tensor.

## Step 8. The Kretschmann scalar

Off the throat the Kretschmann scalar is Schwarzschild's, $12r_s^2/(a + |\ell|)^6$.
On it the square of the curvature holds $\delta(\ell)^2$, which has no value as a distribution, so the scalar is defined only for $\ell \ne 0$, as the chart's convention says.
`off_support` in `verify_metrics.py` compares it there, with every delta at the kink set to zero, and `print_charts.py` does the same before it writes.

## Step 9. The spacetime diagrams

Both are drawn at $r_s = 1$ and $a = 5/4$, inside the photon sphere, where $f(a) = 1/5$.
Through the throat the rays are $ct \mp \ell_* = $ const with $\ell_* = \ell + r_s\ln(1 + |\ell|/(a - r_s))\operatorname{sgn}\ell = \ell + \ln(1 + 4|\ell|)\operatorname{sgn}\ell$, since $d\ell_*/d\ell = 1/f$ and $\ell_*(0) = 0$; `--verify` holds each ray to that closed form.
On one side the rays are Schwarzschild's, $ct \mp r_*$ with $r_* = r + \ln(r - 1)$, from the throat out, and the box starts at $r = a$, where the chart ends.
The throat is drawn as the shell, the line $\ell = 0$ and $r = a$.

## Step 10. The conformal diagram

On the plane of $t$ and $\ell$ the metric is $f\left(-c^2dt^2 + d\ell_*^2\right)$, and $\ell_*$ runs over the whole line, so $p, q = \arctan((ct \mp \ell_*)/r_s)$ bring the plane into Minkowski's full diamond with the throat on its axis, and each side of the throat is half the diamond, as each side of the Ellis-Bronnikov wormhole is.
`conformal.py` checks both maps null against the published metric of each chart, and Schwarzschild's chart on one side is the right half, entered through $\ell = r - a$.

## Step 11. The embedding diagram

The equator of $t = 0$ on either side is Flamm's paraboloid, $dz/dr = \sqrt{r_s/(r - r_s)}$, and from the throat $z = \pm 2(\sqrt{r_s(r - r_s)} - \sqrt{r_s(a - r_s)})$.
At the throat both sides climb at $dz/dr = \sqrt{r_s/(a - r_s)} = 2$, one up and one down, so the surface has a crease there, the embedding of the jump in extrinsic curvature that the shell is.
The script checks each side against that closed form, the two sides for one circle at the throat and for the slope $2$ on either side of it, and the chart through the throat for the same surface.
