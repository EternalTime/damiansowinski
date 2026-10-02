# Tolman's solution VII: the charts, the source and what the diagrams draw

Tolman's seventh solution is a static sphere of perfect fluid whose density falls as the square of the areal radius.
This note records where each chart of `tolman_vii.json` comes from, how the checker reads its phase, and what the diagrams take from it.
Units here have $G = c = 1$ unless a line keeps them; the file keeps both.

## Step 1. Tolman's chart

Tolman (1939) writes a static sphere as $ds^2 = -e^\lambda dr^2 - r^2d\theta^2 - r^2\sin^2\theta\,d\phi^2 + e^\nu dt^2$, his (2.2), and assumes for solution VII

$$e^{-\lambda} = 1 - \frac{r^2}{R^2} + \frac{4r^4}{A^4},$$

with $R$ and $A$ constants.
His (4.7) is the result,

$$e^\nu = B^2\sin^2\ln\left(\frac{e^{-\lambda/2} + 2r^2/A^2 - A^2/4R^2}{C}\right)^{1/2},
\qquad 8\pi\rho = \frac{3}{R^2} - \frac{20r^2}{A^4},
\qquad 8\pi p = -\frac{1}{R^2} + \frac{4r^2}{A^4} + \frac{4e^{-\lambda/2}}{A^2}\left(B^2e^{-\nu} - 1\right)^{1/2},$$

with $B$ and $C$ two more constants and the cosmological constant left out.
The file writes $Z = e^{-\lambda} = g^{rr}$ and $\psi$ for the argument of the sine,

$$\psi = \frac{1}{2}\ln\left(\frac{1}{C}\left(\sqrt{Z} + \frac{2r^2}{A^2} - \frac{A^2}{4R^2}\right)\right),
\qquad g_{tt} = -B^2\sin^2\psi,$$

in the signature $(-,+,+,+)$.
Then $(B^2e^{-\nu} - 1)^{1/2} = \cot\psi$ and the pressure is $8\pi p = -1/R^2 + 4r^2/A^4 + (4\sqrt{Z}/A^2)\cot\psi$.
Tolman's $R$ is the length in the central density, $8\pi\rho_c = 3/R^2$, and is no radius of the star: the surface is at his $r_b$, where the pressure vanishes, his (5.3).
His (5.4) fixes $B$ there, $B^2\sin^2\psi = Z$, and $C$ fixes $r_b$.
`print_charts.tolman_vii("tolman")` builds the chart, and `tolman_vii_check` holds it to the density and the pressure above.

## Step 2. The derivative of the phase

The phase is half a logarithm of radicals, and it stands only inside a sine.
Its derivative is algebraic.
With $X = \sqrt{Z} + 2r^2/A^2 - A^2/4R^2$,

$$\frac{dX}{dr} = \frac{Z'}{2\sqrt{Z}} + \frac{4r}{A^2} = \frac{4r}{A^2}\,\frac{X}{\sqrt{Z}},
\qquad\text{so}\qquad \frac{d\psi}{dr} = \frac{2r}{A^2\sqrt{Z}},$$

since $Z' = (4r/A^2)(4r^2/A^2 - A^2/2R^2)$.
Every tensor is therefore a rational function of $r$, $\sqrt{Z}$ and $\cot\psi$, with $\sin^2\psi$ a factor wherever a time index is lowered.

## Step 3. Lattimer and Prakash's chart

Lattimer and Prakash (2001) take the member whose density vanishes at the surface, $\rho = \rho_c(1 - r^2/R^2)$ with $R$ now the radius of the star, their (16).
With $x = r^2/R^2$ and the compactness $\beta = GM/(Rc^2)$, their (17) is

$$e^{-\lambda} = 1 - \beta x(5 - 3x), \qquad e^\nu = \left(1 - \frac{5\beta}{3}\right)\cos^2\phi,
\qquad \phi = \phi_1 + \frac{w_1 - w}{2}, \qquad w = \ln\left(x - \frac{5}{6} + \sqrt{\frac{e^{-\lambda}}{3\beta}}\right),$$

with $\phi_1 = \arctan\sqrt{\beta/(3(1 - 2\beta))}$ and $w_1 = w(x = 1)$, and

$$P = \frac{c^4}{4\pi R^2G}\left(\sqrt{3\beta e^{-\lambda}}\tan\phi - \frac{\beta}{2}(5 - 3x)\right), \qquad \rho_c = \frac{15\beta c^2}{8\pi GR^2}.$$

The file writes their $\phi$ as $\psi$, since $\phi$ is the azimuth, and $Z$ for $e^{-\lambda}$, and puts the two logarithms under one, with numerator and denominator multiplied by six.
As in Step 2, $dw/dx = \sqrt{3\beta/Z}$, so $d\psi/dr = -(r/R^2)\sqrt{3\beta/Z}$.
At the surface $\psi = \phi_1$, the pressure vanishes, $Z = 1 - 2\beta$ and $g_{tt} = -(1 - 5\beta/3)\cos^2\phi_1 = -(1 - 2\beta)$: Schwarzschild's exterior of the mass $M$ on both counts.
It is Tolman's chart at $R^2 \to R^2/5\beta$, $A^4 = 4R^4/3\beta$, $B^2 = 1 - 5\beta/3$ and $C = \sqrt{3\beta}\,(1/6 + \sqrt{(1 - 2\beta)/3\beta})\,e^{2\phi_1 - \pi}$, which makes his phase $\pi/2$ less theirs.
`tolman_vii_check` holds the chart to the density, the pressure, the three statements at the surface and that map.

The central values are their (18): $\psi_c = \pi/2$, where the central pressure is infinite, at $\beta \approx 0.3862$, and the speed of sound at the centre, $c_s^2/c^2 = \tan\psi_c\,(\tan\psi_c/5 + \sqrt{\beta/3})$, is that of light at $\beta \approx 0.2698$.
`_tools/test_tolman_vii.py` holds both numbers to the published chart.

## Step 4. Why the parameters are $R$ and $r_s$

The checker splits a radical into the roots of the irreducible factors of its radicand, each read as positive and each normalised to lead with a positive coefficient.
$1 - 2\beta$ is normalised to $-(2\beta - 1)$, and the root of that comes out as $i\sqrt{2\beta - 1}$ above the line and below it alike, so $\sqrt{1 - 2\beta}\cdot\sqrt{1/(1 - 2\beta)}$ reads as $-1$.
With the star's radius and Schwarzschild radius for parameters, $1 - 2\beta = (R - r_s)/R$ and every factor under a root leads with a positive coefficient.
So the chart declares $R$ and $r_s$ and names $\beta = r_s/2R$, and every published value is written in $\beta$.

## Step 5. A name held as a function

Written out, $\psi$ puts a logarithm of radicals inside every cosine, and the checker's canonical form spells such a logarithm differently on a second pass than on the first, since it factors the argument before it splits the radicals.
`HELD` in `verify_metrics.py` therefore lists $\psi$ for both charts: the reader holds it as a function of $r$ while the tensors are built, as it holds a defined name that carries a declared function, and `Reader.surface` writes the definition out on both sides of each comparison, once, so that every logarithm in a comparison has one spelling.
The comparison is still made on the definition the file publishes.
`print_charts.TolmanForms.reduce` writes each derivative of the held $\psi$ by Step 2, which `tolman_vii_check` proves of the definition first, and `pretty` writes a value in $\tan\psi$ or $\cot\psi$, in $Z$ and $\sqrt{Z}$, and in $\sqrt{3\beta}$.

## Step 6. What the diagrams draw

Every diagram draws the star of Step 3 at $R = 2r_s$, $\beta = 1/4$, which is below both limits of Step 3.
The spacetime diagrams of Tolman's chart draw the same star in his constants, $R^2 = 16r_s^2/5$, $A^4 = 256r_s^4/3$, $B^2 = 7/12$ and $C \approx 0.0799$, and `conformal.tolman_vii` checks that the two charts have one metric on the plane of $t$ and $r$.
At the centre $g_{tt} = -0.199$.

On a slice of constant $t$ the equator has $g_{rr} = 1/Z$ with $1 - Z = 2m/r$ and $m = \beta r^3(5R^2 - 3r^2)/2R^4$, so the embedded surface climbs at $dz/dr = \sqrt{2m/(r - 2m)}$, which `embedding.tolman_vii` checks by quadrature.
Its Gaussian curvature is $m'/r^2 - m/r^3 = \beta(5R^2 - 6r^2)/R^4$, positive inside $r = R\sqrt{5/6}$ and negative outside, and at $r = R$ it is $-\beta/R^2 = -M/R^3$, Flamm's value, so the curvature is continuous across the surface of the star as well as the tangent.
