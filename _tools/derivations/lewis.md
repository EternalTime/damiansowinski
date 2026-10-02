# The Lewis metrics: the charts and what the diagrams draw

Lewis's family is the general stationary vacuum solution with cylindrical symmetry.
This note records where each chart of `lewis.json` comes from, how the charts are held to one another, and what the diagrams take from them.
Units here have $c = 1$; the file keeps $c$.

## Step 1. Lewis's chart and the unit of length

In the notation of da Silva, Herrera, Paiva, and Santos (1995), which Bronnikov, Santos, and Wang's review (2020, their section on the Lewis vacuum) and Costa, Natário, and Santos (2021) keep, the line element is

$$ds^2 = -f\,dt^2 + 2k\,dt\,d\phi + e^\mu(dr^2 + dz^2) + l\,d\phi^2,$$

$$f = ar^{1 - n} - \frac{c^2r^{n + 1}}{n^2a}, \qquad k = -Af, \qquad l = \frac{r^2}{f} - A^2f, \qquad A = \frac{cr^{n + 1}}{naf} + b, \qquad e^\mu = r^{(n^2 - 1)/2}.$$

Their coordinates $t$, $r$, and $z$ are pure numbers, and their $c$ is a constant of the solution.
The file writes that constant $q$, since $c$ is the speed of light on every page, and measures $r$ in a length $\ell$, the scale Costa, Natário, and Santos call $\mathcal{R}$ where they join the vacuum to van Stockum's cylinder: $a$ and $q$ are then numbers and $b$ is a length.
Every power whose exponent holds $n$ is a power of $r/\ell$, which is what lets the dimensional pass balance each term, and the chart names $u = (r/\ell)^n$ and $h = (r/\ell)^{(n^2 - 1)/2}$.
Since $fl + k^2 = r^2$, the functions $k$ and $l$ have forms with no $f$ below the line,

$$k = -\frac{q\,r\,u}{n\,a} - bf, \qquad l = \frac{r\,u}{a}\left(\ell - \frac{2bq}{n}\right) - b^2f,$$

which the file states, and `lewis_check` holds them to Lewis's $-Af$ and $r^2/f - A^2f$.
With all four constants real this is the Weyl class.
At $b = q = 0$ it is checked to be Levi-Civita's published Weyl chart with $\sigma = (1 - n)/4$ and $C^2 = a$.

## Step 2. The canonical chart

Costa, Natário, and Santos carry Lewis's chart along $\phi \to \phi + \Omega t$ with $\Omega = q/(n\ell - bq)$, and reach their canonical form,

$$ds^2 = -\frac{r^{4\sigma}}{\alpha}\left(dt + \frac{4j}{1 - 4\sigma}d\phi\right)^2 + r^{4\sigma(2\sigma - 1)}(dr^2 + dz^2) + \alpha r^{2(1 - 2\sigma)}d\phi^2,$$

with the Komar mass per unit length $\sigma = (1 - n)/4$, their $\lambda_m$, the Komar angular momentum per unit length $j = b(n\ell - bq)/4\ell$, and $\alpha = (n\ell - bq)^2/(n^2\ell^2a)$.
The file writes the powers on $r/\ell$ and names $u = (r/\ell)^{4\sigma}$ and $h = (r/\ell)^{4\sigma(2\sigma - 1)}$.
`lewis_check` pulls Lewis's chart back along that map and compares every slot at random points in thirty digits, and holds the chart at $j = 0$ to Levi-Civita's published chart with $C^2 = 1/\alpha$.
This is their branch $a > 0$; the branch $a < 0$ gives the same line element with $n \to -n$.

## Step 3. The Lewis class

With $n = im$ the metric is real for the constants the review gives,

$$a = \tfrac{1}{2}(a_1 + ib_1)^2, \qquad b = \frac{a_2 + ib_2}{a_1 + ib_1}, \qquad c = \tfrac{1}{2}m(a_1^2 + b_1^2), \qquad a_1b_2 - a_2b_1 = 1,$$

and the functions become

$$f = r\left((a_1^2 - b_1^2)\cos\psi + 2a_1b_1\sin\psi\right), \quad k = -r\left((a_1a_2 - b_1b_2)\cos\psi + (a_1b_2 + a_2b_1)\sin\psi\right), \quad l = -r\left((a_2^2 - b_2^2)\cos\psi + 2a_2b_2\sin\psi\right),$$

with $\psi = m\ln r$ and $e^\mu = r^{-(m^2 + 1)/2}$.
The review's arXiv text prints $a = \tfrac{1}{2}(a_1 + b_1)^2$, without the $i$, which does not make $f$ real; the form above does, and `lewis_check` continues Lewis's chart to these complex constants and compares it with the Lewis class's chart at random points.
The file takes $b_2$ as the name $(1 + a_2b_1)/a_1$, so no published value holds a relation it does not state, and the chart asks $a_1 \neq 0$.

## Step 4. van Stockum's exteriors

Van Stockum's exterior for $wR < 1/2$ is, as Costa, Natário, and Santos write it with $N = n/2$:

$$F = \frac{(n - 1)x^{n + 1} + (n + 1)x^{1 - n}}{2n}, \quad M = wR^2\frac{(n + 1)x^{n + 1} + (n - 1)x^{1 - n}}{2n}, \quad L = R^2\frac{(n + 1)^3x^{n + 1} + (n - 1)^3x^{1 - n}}{8n},$$

$$H = e^{-w^2R^2}x^{-2w^2R^2}, \qquad x = r/R, \qquad n = \sqrt{1 - 4w^2R^2}.$$

The file takes $n$ and $R$ as the constants and $w = \sqrt{1 - n^2}/2R$ as a name, and writes the cross term as $-2M\,dt\,d\phi$, the sense of rotation of the published chart of the dust in `stockum_dust.json`, whose own $R$ is $1/w$.
$H$ is $(r/\ell)^{-(1 - n^2)/2}$ with $\ell = R/\sqrt{e}$, and the file keeps $\ell$ as a constant of the chart: the checker has no generator for the number $e^{-1/4}$, and the vacuum equations hold for every $\ell$, the junction alone fixing it.
The exterior for $wR > 1/2$ is the same with $n = im$, and the exterior at $wR = 1/2$ is its limit $n \to 0$:

$$F = x\left(\cos\psi - \frac{\sin\psi}{m}\right), \quad M = wR^2x\left(\cos\psi + \frac{\sin\psi}{m}\right), \quad L = \frac{R^2x}{4}\left((3 - m^2)\cos\psi + \frac{(1 - 3m^2)\sin\psi}{m}\right), \quad \psi = m\ln x,$$

$$F = x(1 - \ln x), \qquad M = \tfrac{1}{2}Rx(1 + \ln x), \qquad L = \tfrac{1}{4}R^2x(3 + \ln x),$$

which are the forms Tipler (1974) gives.
`lewis_check` holds the light cylinder to being Lewis's chart at $a = (n + 1)(\ell/R)^{1 - n}/2n$, $b = -(1 - n)wR^2/(1 + n)$, and $q = w\ell$, which are Costa, Natário, and Santos's substitutions with the signs of $b$ and $q$ reversed for the reversed rotation, the critical cylinder to the limit, the heavy one to the continuation and to being a member of the Lewis class, and each of the three, at $\ell = R/\sqrt{e}$, to carrying the published dust's metric and its first derivative along $r$ at $r = R$.

## Step 5. The radical of the light cylinder

The checker splits a radical over the irreducible factors of its radicand, each with a positive leading coefficient, so it writes $\sqrt{1 - n^2}$ as $i\sqrt{n - 1}\sqrt{n + 1}$.
That is the same algebraic number, and every comparison the checker makes is an identity in it, but taken with principal roots at $0 < n < 1$ it is minus the radical.
So a value of this chart that has been through `norm` is not evaluated numerically as it stands: `lewis_check` puts the radical back on the branch $0 < n < 1$ first, and the diagram scripts read the published values through the `Reader` alone, which keeps $\sqrt{1 - n^2}$ whole.
A rational parametrisation of the circle $n^2 + 4w^2R^2 = 1$ in `PARAMETER_RELATIONS` was tried and left: it puts a rational function of the parameter into every exponent, which the checker's generators do not take.

## Step 6. What the diagrams draw

The Weyl class is drawn at $n = 1/2$: in Lewis's chart with $a = 4/9$, $b = 3\ell/2$, and $q = 1/9$, and in the canonical chart with $\sigma = 1/8$, $j = \ell/8$, and $\alpha = 1$, which Step 2's map makes one spacetime.
Its circle $r = \ell$ is null, with closed timelike curves inside, and Lewis's $\partial_t$ is spacelike beyond $r = 4\ell$, where the rigidly turning coordinates outrun light.
On a cylinder of $t$ and $\phi$ the null curves are $dt = (k \pm r)\,d\phi/f$, and in van Stockum's letters $dt = (-M \pm r)\,d\phi/F$.
The canonical chart's plane of $t$ and $r$ has $d\phi = 0$, so its metric is Levi-Civita's, $u(-dt^2 + dr_*^2)/\alpha$ with $r_* = (16/9)r^{9/16}$ at $\sigma = 1/8$; its null rays are null geodesics, and its conformal diagram is Levi-Civita's triangle.
The future is taken from the static time $t + C\phi$ in the Weyl class's own charts, from $t$ outside the light and critical cylinders, whose $L$ is positive at every radius, and from $\phi$ in the bands where $f$ or $F$ is negative; there $g_{t\phi} < 0$ puts $+\partial_\phi$ in the cone that holds $+\partial_t$ nearer the cylinder.
At $wR = 1/2$ the curve $dt = (R/2)\,d\phi$ is a null geodesic at every radius, the null limit of the circular orbits of Levi-Civita's cylinder at $\sigma = 1/4$.
The embedding diagram is the moment $t = 0$ of the canonical chart outside the null circle, $h\,dr^2 + (\alpha r^2/u - C^2u/\alpha)\,d\phi^2$: in Minkowski space from $r = \ell$ to the circle where $g_{rr} = (\partial_r\sqrt{g_{\phi\phi}})^2$, the root $r = 1.472\,\ell$ of $16s^{17} - 9s^{16} - 16s^9 + 6s^8 - 1$ in $s = r^{1/8}$, and in flat space beyond, as the spinning string's is.
A figure of light cones in three dimensions is not drawn: the builder of the spinning string's measures the proper distance from an axis where the metric is regular, and the two cylinders of the canonical chart show the same tipping.
