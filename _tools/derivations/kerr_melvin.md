# Kerr's black hole in Melvin's universe

Ernst and Wild's metric, Harrison's transformation of Kerr's, in Boyer and Lindquist's coordinates of the seed is

$$ds^2 = H\left(-\frac{\Delta\Sigma}{A}c^2dt^2 + \Sigma\left(\frac{dr^2}{\Delta} + d\theta^2\right)\right) + \frac{A\sin^2\theta}{H\Sigma}\left(k\,d\phi - \omega\,c\,dt\right)^2,$$

with Kerr's $\Sigma = r^2 + a^2\cos^2\theta$, $\Delta = r^2 - 2mr + a^2$ and $A = (r^2 + a^2)^2 - a^2\Delta\sin^2\theta$, and

$$H = \left(1 + \frac{B^2A\sin^2\theta}{4\Sigma}\right)^2 + \frac{B^4m^2a^2\cos^2\theta}{4}\left(3 - \cos^2\theta + \frac{a^2\sin^4\theta}{\Sigma}\right)^2,$$

$$\omega = \frac{a}{A}\left(2mr + \frac{mB^4}{8}\left(\left(3r^2 + 6mr - a^2\right)\left(r^3 + a^2r + 2a^2m\right) + \Delta\cos^2\theta\left(6r\left(r^2 + a^2\right) - \left(r^3 - 3a^2r + 2a^2m\right)\cos^2\theta\right)\right)\right),$$

$$k = 1 + a^2m^2B^4.$$

Both charts are written by `_tools/derivations/print_charts.py --metric kerr_melvin`, one chart to a run with `--system`, and `verify_metrics.py --system kerr_melvin/boyer_lindquist` and `--system kerr_melvin/rotating` check them.

## Step 1. The sources

Ernst and Wild (1976) state in their abstract that they give the electromagnetic field of a Kerr-Newman hole of any charge in a magnetic universe and the metric "in the particular case of charge $Q = 2B_0J$", $J = ma$.
The paper itself was not read: its publisher refuses the request, and the abstract is the record on INSPIRE and Crossref.
That case is the transformation of the uncharged Kerr seed, as Bičák and Hejda (2015, section II B) say, and they print its metric as $|\Lambda|^2\Sigma[-(\Delta/\mathcal{A})dt^2 + dr^2/\Delta + d\vartheta^2] + (\mathcal{A}/(\Sigma|\Lambda|^2))\sin^2\vartheta\,(d\varphi - \omega\,dt)^2$ with

$$\Lambda = 1 + \frac{B^2}{4}\frac{\mathcal{A}}{\Sigma}\sin^2\vartheta - \frac{i}{2}B^2Ma\cos\vartheta\left(3 - \cos^2\vartheta + \frac{a^2}{\Sigma}\sin^4\vartheta\right).$$

$H$ is $|\Lambda|^2$, the square of the real part plus the square of the imaginary part.
Gibbons, Mujtaba, and Pope (2013, appendix B) print the metric of the magnetized Kerr-Newman hole whole, with $H = 1 + (H_1B + H_2B^2 + H_3B^3 + H_4B^4)/R^2$ and $\omega = ((2mr - \tilde q^2)a + \omega_1B + \dots + \omega_4B^4)/\Sigma$; their $R^2$ is $\Sigma$ here and their $\Sigma$ is $A$.
At $q = p = 0$ only $H_2$, $H_4$ and $\omega_4$ are left, and

$$1 + \frac{H_2B^2 + H_4B^4}{\Sigma} = |\Lambda|^2,
\qquad
\frac{8\,\omega_4}{am} = \left(3r^2 + 6mr - a^2\right)\left(r^3 + a^2r + 2a^2m\right) + \Delta\cos^2\theta\left(6r\left(r^2 + a^2\right) - \left(r^3 - 3a^2r + 2a^2m\right)\cos^2\theta\right),$$

both identities of polynomials, checked in sympy when the chart was first written; the second is how $\omega$ is printed.
`_tools/test_kerr_melvin.py` holds the published $H$ and $\omega$ to Gibbons, Mujtaba, and Pope's polynomials as they print them.

## Step 2. The period of the azimuth

On the axis $\sin\theta = 0$ and $\cos^2\theta = 1$, so $H = 1 + B^4m^2a^2 = k$ at both poles.
With the seed's angle $\phi_s$ of period $2\pi$ the circle of latitude $\theta$ near a pole has radius $\sqrt{g_{\phi_s\phi_s}} \to (r^2 + a^2)\theta/(k\sqrt{\Sigma})$ at the distance $\sqrt{g_{\theta\theta}}\,\theta = \sqrt{k\Sigma}\,\theta$ from it, a ratio of $1/k$ on the axis, where $\Sigma = r^2 + a^2$: a cone.
Hiscock (1981) found it, and gave the angle the period $2\pi k$ instead.
Booth, Hunt, Palomo-Lozano, and Kunduri (2015, their equation for the line element) write the same thing with an angle of period $2\pi$, $\phi_s = k\phi$, and the charts here follow them, so that $g_{\phi\phi} = k^2A\sin^2\theta/(H\Sigma)$ is the square of the circle's radius and every drawing reads a circumference from the published metric as it does for Kerr.
`kerr_melvin_check` holds $H$ to $k$ at both poles.

## Step 3. Five names each chart holds

Written out, a component of the Riemann tensor is a rational function whose numerator has thousands of terms: with $H$ and $\omega$ alone held as functions of $r$ and $\theta$ the Riemann tensor did not print in five minutes.
Held together with Kerr's $\Sigma$, $\Delta$ and $A$, in Ernst and Wild's form of the line element, the chart printed in 75 seconds and ran to a megabyte, since every metric component is then a product of three or four functions and each derivative spreads over all of them.
Each chart therefore writes the line element in the four functions of a stationary field with an axis,

$$ds^2 = -N\,c^2dt^2 + F\left(\frac{dr^2}{\Delta} + d\theta^2\right) + P\left(k\,d\phi - \omega\,c\,dt\right)^2,
\qquad
N = \frac{H\Delta\Sigma}{A},
\quad
F = H\Sigma,
\quad
P = \frac{A\sin^2\theta}{H\Sigma},$$

declares $\Sigma$, $\Delta$, $A$, $H$, $\omega$, $N$, $F$ and $P$ among its parameters, and holds $\Delta$, $\omega$, $N$, $F$ and $P$ as functions while the tensors are built, `HELD` in `verify_metrics.py`, so that every value is written in the five and their derivatives.
$a$, $m$ and $B$ then enter a value only through $k$, which `kerr_melvin_pretty` writes back as $k$.
The reader writes each held definition out in full as it reads it, since $N$ holds $H$, which holds $A$, which holds $\Delta$, and `agree_held` passes two values that already agree as functions of the held names without writing them out: equal in the names, they are equal for every definition of them.
The Ricci scalar is published in the names as every other value is, a sum of their second derivatives that no relation among five free functions makes zero.
For the definitions it vanishes, the trace of Maxwell's stress.
Written out from the five names its canonical form took 115 seconds, against the 120 the checker allows one tensor, so the checker is not asked for it: `kerr_melvin_traceless` builds the scalar from Ernst and Wild's own form of the line element with $H$ and $\omega$ alone held, writes it out and finds the canonical form zero, in about a minute and a half each time the first chart is printed, and `_tools/test_kerr_melvin.py` holds the published expression to zero at a point of rational coordinates.

## Step 4. What the check holds

`kerr_melvin_check` writes the five names out and holds each chart to the Einstein-Maxwell equations, $G_{\mu\nu} = 2\left(F_{\mu\alpha}F_\nu{}^\alpha - \tfrac{1}{4}g_{\mu\nu}F^2\right)$ in geometric units, at two points of rational coordinates and parameters, exactly, with the potential of Gibbons, Mujtaba, and Pope's appendix B at $q = p = 0$,

$$\mathcal{A} = \Phi_0\,c\,dt + \Phi_3\left(k\,d\phi - \omega\,c\,dt\right),
\qquad
\Phi_3 = \frac{\partial_BH}{2H},
\qquad
\Phi_0 = \frac{2}{B}\left(\omega - \frac{2mra}{A}\right).$$

Their $\Phi_3 = (\chi_1B + \chi_3B^3)/(R^2H)$ has $\chi_1 = H_2$ and $\chi_3 = 2H_4$, which is the derivative of $H$ along $B$ over $2H$, and their $\Phi_0 = \Phi_0^{(3)}B^3/(4\Sigma)$ has $\Phi_0^{(3)} = 8\omega_4$, which is the second.
It holds $\omega$ on the horizon to $a/(2mr_+) + amB^4(r_+ + m)/2$ at every $\theta$, the first chart at $B = 0$ to the published Boyer-Lindquist chart of `kerr` with $GM/c^2 = m$ and at $a = 0$ to the published chart of Ernst's hole in `melvin` with $r_s = 2m$, and the second chart to the first pulled back along $\phi = \tilde\phi + \Omega\,ct$.

## Step 5. The charge, and the rate of the horizon

The charge inside a sphere round the hole is $Q = \tfrac{1}{4\pi}\oint {*F}$, and with Gibbons, Mujtaba, and Pope's $\psi = am\left(3r^2 + a^2 - (r^2 - a^2)\cos^2\theta\right)\cos\theta\,B/(\Sigma H)$ for the uncharged seed it is $\tfrac{k}{2}\left[\psi\right]_{\theta = \pi}^{\theta = 0}$ in size: $\psi = \pm 2amB/k$ at the poles, so $|Q| = 2amB$, Wald's charge with the seed's $J = am$.
`_tools/test_kerr_melvin.py` integrates the flux of the published field's dual over a sphere and holds it to $2amB$.
On the horizon $\Delta = 0$, $A = (r_+^2 + a^2)^2 = 4m^2r_+^2$, and $\omega$ is the same at every $\theta$, as the rigidity of a Killing horizon asks: $\omega_+ = a/(2mr_+) + amB^4(r_+ + m)/2$.
The horizon turns at $d\phi/d(ct) = \omega_+/k$, and the second chart at $\Omega = \omega_+/k$ has $\partial_t$ for the horizon's null generator.
