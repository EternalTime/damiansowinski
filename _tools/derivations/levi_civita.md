# Levi-Civita's static cylinder

Levi-Civita's static vacuum with cylindrical symmetry, in Weyl's coordinates and in its Kasner form, is

$$ds^2 = -\rho^{4\sigma}c^2dt^2 + \rho^{4\sigma(2\sigma - 1)}\left(d\rho^2 + dz^2\right) + \frac{\rho^{2(1 - 2\sigma)}}{C^2}\,d\phi^2,$$

$$ds^2 = -r^{2p_0}c^2dt^2 + dr^2 + \ell^2r^{2p_2}\,d\phi^2 + r^{2p_3}\,dz^2,\qquad p_0 + p_2 + p_3 = p_0^2 + p_2^2 + p_3^2 = 1.$$

Both charts are written by `_tools/derivations/print_charts.py --metric levi_civita`, and `verify_metrics.py --system levi_civita/weyl` and `--system levi_civita/kasner` check each in seconds.

## Step 1. Units and the date

Levi-Civita's ninth note is *Rendiconti della Reale Accademia dei Lincei* **28**, 101-109 (1919); volume 26, page 307 (1917), which da Silva, Herrera, Paiva and Santos cite, is the first note of the series, as MacCallum's editorial note says.
The mass parameter is $\sigma = G\lambda/c^2$ for a line of mass $\lambda$ per unit length: $g_{tt} = -\rho^{4\sigma} \approx -(1 + 4\sigma\ln\rho)$ is $-(1 + 2\Phi/c^2)$ with Newton's $\Phi = 2G\lambda\ln\rho$.
A power of the radius whose exponent holds $\sigma$ is read with the radius in a fixed unit of length, as Kasner's $t^{2p_i}$ is, so only the numeric part of an exponent carries a dimension: $\rho^{2 - 4\sigma}$ is an area, and $C$ is a pure number.
`DIMENSIONS` declares that, and $[\ell] = L$ in the Kasner form.

## Step 2. The Kasner form

With $\Sigma = 4\sigma^2 - 2\sigma + 1$, the proper distance from the axis is $r = \int\rho^{2\sigma(2\sigma - 1)}d\rho = \rho^\Sigma/\Sigma$.
Then $\rho^{4\sigma} = (\Sigma r)^{4\sigma/\Sigma}$, and absorbing the constants into $t$ and $z$ leaves the Kasner form with

$$p_0 = \frac{2\sigma}{\Sigma},\qquad p_2 = \frac{1 - 2\sigma}{\Sigma},\qquad p_3 = \frac{2\sigma(2\sigma - 1)}{\Sigma},\qquad \ell = \frac{\Sigma^{p_2}}{C}.$$

Their sum is $\Sigma/\Sigma = 1$, and the sum of their squares is $\left(4\sigma^2 + (2\sigma - 1)^2(4\sigma^2 + 1)\right)/\Sigma^2 = 1$.
`PARAMETER_RELATIONS` carries this parametrisation of Kasner's circle, so the checker compares the Kasner form's values on the circle, where alone they are a vacuum.
`levi_civita_on_kasner_circle` in `print_charts.py` checks the Ricci tensor to vanish there, writes the Weyl tensor as the Riemann tensor, and checks $K = -16p_0p_2p_3/r^4$, which is $64\sigma^2(2\sigma - 1)^2/(\Sigma^3r^4)$, Herrera, Santos, Teixeira and Wang's eq. (12).

## Step 3. Printing a power of the radius

`chart_printer.py` writes rational exponents only.
`radial_powers` in `print_charts.py` gathers the powers of the radius in a value, which is one monomial in it, and hands the printer a placeholder whose text is the radius with its exponent.
The whole part of a negative exponent is written below the line, as $2\sigma\rho^{8\sigma - 8\sigma^2}/\rho$, so the dimension of a value is the dimension of what stands outside the powers that hold $\sigma$.
Every printed value is read back and compared as for any other chart.

## Step 4. The geometry

The Kretschmann scalar is $64\sigma^2(2\sigma - 1)^2\Sigma\,\rho^{-4\Sigma}$, zero at $\sigma = 0$ and $\sigma = 1/2$ and divergent on the axis otherwise.
On the plane of $t$ and $\rho$ the null curves are $ct = \pm\rho_* + $ const with $\rho_* = \rho^{(2\sigma - 1)^2}/(2\sigma - 1)^2$, which is zero on the axis, so a ray leaves the axis and reaches any $\rho$ in a finite $t$: no horizon.
A circular geodesic has $\dot\phi^2/\dot t^2 = -\partial_\rho g_{tt}/\partial_\rho g_{\phi\phi}$, and its speed past the static observers is $W^2 = g_{\phi\phi}\dot\phi^2/(-g_{tt}\dot t^2) = 2\sigma/(1 - 2\sigma)$ at every radius, the speed of light at $\sigma = 1/4$.
On the slice $t = 0$, $z = 0$ the circle of radius $\rho^{1 - 2\sigma}/C$ grows at $(1 - 2\sigma)\rho^{-4\sigma^2}/C$ times the distance out to it, which is below 1 only beyond $\rho_0 = ((1 - 2\sigma)/C)^{1/(4\sigma^2)}$.

## Step 5. The diagrams

Every diagram is drawn at $\sigma = 1/4$ and $C = 1$, where the Kasner exponents are $(2/3, 2/3, -1/3)$ and $\ell = (3/4)^{2/3}$.
There the circles grow as the power $2/3$ of the distance out to them, a horn the eye tells from a cone; at a small $\sigma$ the surface is within a few percent of the cosmic string's cone.
The spacetime diagrams draw the plane of $t$ and $\rho$ and the plane of $t$ and $r$, and `null_rays.py --verify` compares their rays with $ct \pm 4\rho^{1/4}$ and $ct \pm 3r^{1/3}$.
The conformal diagram is Minkowski's triangle, $p, q = \arctan(ct \mp 4\rho^{1/4})$, with the axis a timelike singularity on $X = 0$, one view for each chart, the same triangle since $4\rho^{1/4} = (4/3)^{1/3}\,3r^{1/3}$.
The embedding diagram is the plane $z = 0$ from $\rho_0 = 1/16$ to $\rho = 4$, the horn $z = (4\sqrt\rho - 1)^{3/2}/6$ with circles of radius $\sqrt\rho$, checked against both closed forms; inside $\rho_0$ no surface of revolution in flat space carries the slice, which the view states and the script checks.
