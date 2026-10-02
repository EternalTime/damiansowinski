# Maitra's rotating dust: the chart, the source and what the diagrams draw

Maitra's solution has one parameter, a length $a$, and this note records where the chart of `maitra_dust.json` comes from, what its source is, and what the diagrams take from it.
Units here have $c = 1$ unless a line keeps it; the file keeps $c$.

## Step 1. What was read

Maitra's paper (J. Math. Phys. 7, 1025, 1966) was verified on Crossref and its abstract read; the paper itself was not opened.
Its functions were taken from two later papers that print them and agree.
Krasiński (1998), under his (6.21), writes Maitra's solution in the signature $(+,-,-,-)$ with $g_{t\phi} = m$, $g_{\phi\phi} = m^2 - r^2$, the timelike Killing field of unit length, and his (6.24) and (6.25),

$$m = -\frac{a}{2}\left(\sqrt{1 + u^2} - 1 - \ln\frac{\sqrt{1 + u^2} + 1}{2}\right), \qquad \ln\frac{dy}{dr} = -\frac{\sqrt{1 + u^2} - 1}{4u^2} + \frac{1}{8} - \frac{1}{4}\ln\frac{\sqrt{1 + u^2} + 1}{2}, \qquad u = \frac{2r}{a},$$

with $y$ the proper distance from the axis.
Chan and Santos (2024) derive the same solution by integrating the field equations: their (1) with $f = 1$, $l = r^2 - k^2$ (45), $k$ in (48) and $\gamma$ in (50), with $x = 2c_1r$.
The two agree with $k = -m$, $\gamma = 2\ln(dy/dr)$ and $c_1 = 1/a$.

## Step 2. The chart

In the signature $(-,+,+,+)$ the line element is

$$ds^2 = -\left(c\,dt - k\,d\phi\right)^2 + r^2d\phi^2 + e^{\gamma}\left(dr^2 + dz^2\right),$$

with $s = \sqrt{1 + 4r^2/a^2}$, $k = \tfrac{a}{2}\left(s - 1 - \ln\tfrac{s + 1}{2}\right)$ and $\gamma = \tfrac{1}{4} - \tfrac{1}{2(s + 1)} - \tfrac{1}{2}\ln\tfrac{s + 1}{2}$.
The slopes are $k' = 2r/a(s + 1)$, which is Chan and Santos's (47), and $\gamma' = -k'^2/2r$, their (49).
The checker holds $k$ and $\gamma$ as functions with those slopes, `HELD` and `RATES` in `verify_metrics.py`, having checked each against its definition, so every tensor is a rational function of $r$, $a$, $s$, $k$ and $e^{\gamma}$, and the chart prints in seconds.
`maitra_pretty` writes every even power of $r$ by $r^2 = a^2(s^2 - 1)/4$, after which no relation holds among what is left and the value factors.

## Step 3. The source

`maitra_dust_check` holds the chart to dust: $G_{ab} = R\,u_au_b$ with every other slot zero, where the Ricci scalar

$$R = \frac{4e^{-\gamma}}{a^2s(s + 1)^2} = \frac{8\pi G\rho}{c^2}$$

is $1/a^2$ on the axis and falls as $r^{-5/2}$ far away.
The dust's velocity is $u = v\,\partial_t + \Omega\,\partial_\phi$ with $\Omega = -\tfrac{1}{a}\sqrt{2/(s + 1)}$ and $v = \sqrt{(s + 1)/2} + k\Omega$, Chan and Santos's (4), (6) and (41) at $f = 1$: it circles toward $-\phi$, more slowly farther out, and its speed past the lines of constant $r$, $\phi$ and $z$ is $k'$, which is $0.618\,c$ at $r = a$ and tends to $c$.
The check is made on the three products $u_t^2$, $u_tu_\phi$ and $u_\phi^2$, which hold no nested root, and on $\Gamma^r{}_{ab}u^au^b = 0$ with the root cleared.
It also holds the determinant to $-r^2e^{2\gamma}$, $g^{tt}$ to $-(1 - k^2/r^2)$, $k'^2$ to $(s - 1)/(s + 1) < 1$, Whittaker's mass per unit length to $(s - 1)/8(s + 1)$, which tends to $1/8$, their (51), and the chart to the published dust of `stockum_dust` at $R = 2a$ with $\phi$ reversed through the second order in $r$.

Since $k' < 1$ and $k(0) = 0$, $k < r$ at every radius.
So $g_{\phi\phi} = r^2 - k^2$ is positive and the gradient of $t$ is timelike everywhere: $t$ rises along every timelike curve, and none closes, which is Maitra's title.

## Step 4. The diagrams

- The plane of $t$ and $r$ with $\phi$ divided out, whose curves are the null geodesics of no angular momentum, $c\,dt/dr = \pm e^{\gamma/2}\sqrt{1 - k^2/r^2}$, checked against the quadrature `_maitra_plane`.
- The cylinders of $t$ and $\phi$ at $r = a$ and $r = 5a$, where the null curves are $c\,dt = (k \pm r)\,d\phi$; light launched along either family is turned away from the axis at both radii, $\Gamma^r{}_{ab}k^ak^b < 0$.
- The figure of light cones about the axis, `projections.maitra_dust`, with the proper distance as its radius and no null circle.
- The embedding of the plane $z = 0$, `embedding.maitra_dust`: the circle at $r$ has radius $\sqrt{r^2 - k^2}$, and $e^{\gamma} - (\partial_r\sqrt{r^2 - k^2})^2$ is positive at every radius, checked out to $1000\,a$, so the surface never stops; it is drawn to $6a$.

No conformal diagram is drawn: the planes of $t$ and $r$ and of $t$ and $z$ off the axis do not hold their own light rays, as for van Stockum's dust and Gödel's universe.
