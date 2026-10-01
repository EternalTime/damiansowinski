# Teo's rotating wormhole

The two charts of `teo_wormhole.json` are written by `print_charts.py --metric teo_wormhole`, which computes every tensor from each line element with the checker's own `Geometry` and reads every printed value back before it writes.
This note says what the charts are, why the example is published in place of the canonical metric, and what the diagrams draw.
The source is E. Teo, "Rotating traversable wormholes", Phys. Rev. D 58, 024014 (1998), arXiv:gr-qc/9803098, and equation numbers below are his.

## Step 1. The canonical metric and the example

Teo's canonical metric (19) is

$$ds^2 = -N^2dt^2 + \left(1 - \frac{b}{r}\right)^{-1}dr^2 + r^2K^2\left[d\theta^2 + \sin^2\theta\,(d\phi - \omega\,dt)^2\right],$$

with $N$, $b$, $K$ and $\omega$ functions of $r$ and $\theta$.
His example (26), in units $G = c = b = 1$, is $N = K = 1 + (4a\cos\theta)^2/r$, $b = 1$ and $\omega = 2a/r^3$, with $a$ the angular momentum.
The canonical metric with four free functions of two coordinates was tried first: the checker's `Geometry` had not finished its Riemann tensor after ten minutes, and what it would print is far longer than a page holds.
So the file publishes the example, whose every tensor is a rational function of $r$, $\cos\theta$ and $\sin\theta$, and the history carries the canonical metric.

## Step 2. Units

With the throat radius $b_0$ restored, lengths are in units of $b_0$, and Teo's $a$ is the pure number $a = GJ/(c^3b_0^2)$ for an angular momentum $J$.
Then $N = K = 1 + 16a^2b_0\cos^2\theta/r$, $b = b_0$ and $\omega = 2ab_0^2c/r^3$, an angular velocity, so the line element is

$$ds^2 = -N^2c^2dt^2 + \frac{dr^2}{1 - b_0/r} + r^2N^2\left(d\theta^2 + \sin^2\theta\left(d\phi - \frac{2ab_0^2}{r^3}c\,dt\right)^2\right).$$

Tsukamoto and Bambi, Phys. Rev. D 91, 084013 (2015), restore a length $d$ in the same place, $N = K = 1 + 16a^2d\cos^2\theta/r$, and keep $a$ an area; with $b = d$ theirs is this metric with their $a$ equal to $ab_0^2$ here.
In the chart $x^0 = ct$ the cross term is $g_{t\phi} = -2ab_0^2N^2\sin^2\theta/r$ and $g_{tt} = -N^2\left(1 - 4a^2b_0^4\sin^2\theta/r^4\right)$.
The ergosurface $g_{tt} = 0$ is $r^2 = 2|a|b_0^2\sin\theta$, which lies outside the throat for some $\theta$ exactly when $|a| > 1/2$, Teo's condition.

## Step 3. The proper distance chart

Teo's (15) defines $dl/dr = \pm(1 - b/r)^{-1/2}$, and (16) is the metric in $(t, l, \theta, \phi)$, with $dl^2$ in place of the radial term.
For $b = b_0$ his (28) gives $l = \pm\left(\sqrt{r(r - b_0)} + b_0\ln\left(\sqrt{r/b_0} + \sqrt{r/b_0 - 1}\right)\right)$, which has no inverse in closed form.
The chart therefore declares $r = r(l)$ a function, as the proper distance chart of `morris_thorne` does, and its values hold for every $r(l)$; `teo_wormhole_pullback` checks that with $\partial_l r = \sqrt{1 - b_0/r}$ it is the spherical chart pulled back.
The first printing named $N$ in the proper distance chart as well, but a name defined through a declared function is held as a function while the tensors are built, which the diagram generators do not read, so both charts write $rN = r + 16a^2b_0\cos^2\theta$ out, the factor every value is printed around.

## Step 4. The inverse of (28) for the diagrams

With $r = b_0\cosh^2w$, (28) is $l/b_0 = w + \tfrac12\sinh 2w$, Kepler's kind of equation.
`teo_inverse` in `null_rays.py` solves it by Newton's method from $w = \tfrac12\,\mathrm{asinh}(2l/b_0)$, and returns $\rho = r/b_0 = \cosh^2w$ and $\sigma = d\rho/d(l/b_0) = \tanh w$, which is $\pm\sqrt{1 - b_0/r}$, odd and smooth through the throat.
`teo_rho` and `teo_sigma` are sympy functions with $\rho' = \sigma$ and $\sigma' = 1/(2\rho^2)$, since $d\tanh w/dw = \mathrm{sech}^2w$ and $dw/d(l/b_0) = 1/(2\cosh^2w)$, and `lambdify` evaluates them with `teo_inverse`; a row of `null_rays.py` names them through `DECLARED_FUNCTIONS`.
`conformal.py` and `embedding.py` hand the same numbers to their `numeric`.

## Step 5. What the diagrams draw

On the axis $\sin\theta = 0$ removes the dragging term and $\Gamma^\theta{}_{tt}$, $\Gamma^\theta{}_{rr}$ vanish, so the plane of $t$ and $r$ is totally geodesic, with the metric $-N^2c^2dt^2 + dl^2$ and $N = 1 + 16a^2b_0/r$.
At $a = 1/4$, $N = 1 + b_0/r$ and $dl_* = dl/N$ integrates, with $r = b_0\cosh^2w$, to

$$l_* = \pm b_0\left(\sqrt{\rho(\rho - 1)} - \ln\left(\sqrt{\rho} + \sqrt{\rho - 1}\right) + \sqrt{2}\,\mathrm{artanh}\sqrt{\frac{\rho - 1}{2\rho}}\right), \qquad \rho = r/b_0,$$

since $2\cosh^4w/(\cosh^2w + 1) = 2\sinh^2w + 2/(\cosh^2w + 1)$.
The rays keep $ct \mp l_*$, which `null_rays.py --verify` checks, and $\arctan((ct \mp l_*)/b_0)$ give the conformal diagram's diamond.
On the equator $N = 1$, and with the circles of $\phi$ divided out, `quotient="phi"`, the metric is $-c^2dt^2 + dl^2$ for every spin: the rays of no angular momentum keep $ct \mp l$.
The equatorial slice of constant $t$ has $g_{rr} = r/(r - b_0)$ and $g_{\phi\phi} = r^2$, so it is Flamm's paraboloid $z = \pm 2\sqrt{b_0(r - b_0)}$, Teo's (27).
The throat $r = b_0$ at constant $t$ has the metric $b_0^2N^2(d\theta^2 + \sin^2\theta\,d\phi^2)$ with $N = 1 + 16a^2\cos^2\theta$; at $a = 1/4$ its circles have radius $b_0(1 + \cos^2\theta)\sin\theta$, least at the equator and greatest where $\cos^2\theta = 1/3$, and $g_{\theta\theta} - (d\rho/d\theta)^2 = b_0^2\left((1 + c^2)^2 - c^2(3c^2 - 1)^2\right) \ge 0$ with $c = \cos\theta$, so it embeds whole, Teo's dumbbell.

## Step 6. The reader and a subscript

Printing the line element first returned a diagonal metric: the reader took the `ab_0` of $2ab_0^2$ as one unknown token, split it letter by letter, and the stray `_0` multiplied the dragging term by zero.
`Reader._part_subscripted` in `verify_metrics.py` now parts a declared name with a subscript from the letters written before it, and refuses any other unknown token that holds an underscore, so such a term can no longer vanish without a word.
The whole collection was checked again after the change, and every other system reads as before.
