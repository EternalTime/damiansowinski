# Kiselev

Kiselev's black hole is his Eq. (18) with one term in the sum,

$$ds^2 = -f\,c^2dt^2 + \frac{dr^2}{f} + r^2\,d\Omega^2,\qquad f = 1 - \frac{r_s}{r} - \left(\frac{r_q}{r}\right)^{3w + 1},$$

with $r_s$ his $r_g$, the Schwarzschild radius, $r_q$ his normalisation length, and $w$ his state parameter $w_q$.
Its six charts are written by `_tools/derivations/print_charts.py --metric kiselev`, and `verify_metrics.py --system kiselev/<chart>` checks each in seconds.
Equation numbers are those of arXiv:gr-qc/0210040v3, the text of Classical and Quantum Gravity 20, 1187, counted from its source, and Visser's are those of arXiv:1908.11058.

## Step 1. The parameters and the charts

The parameters are Kiselev's, and both lengths keep the metric function a pure number for every $w$, which a constant $K$ in $K/r^{3w+1}$, Visser's form, does not: $K = r_q^{3w+1}$ carries a power of a length that changes with $w$.
Visser's $K$ may have either sign, and Kiselev's form holds the sign that gives a positive energy density for $w < 0$, his $c\,w_q \ge 0$.
Each chart has a source:

- `static`, Kiselev's (18) with one term.
- `eddington_finkelstein_outgoing`, the retarded time of Toshmatov, Stuchlík and Ahmedov's (5) and (6), $du = dt - dr/f$ and $ds^2 = -f\,du^2 - 2\,du\,dr + r^2d\Omega^2$, from which they and Ghosh build the rotating metric.
- `eddington_finkelstein_ingoing`, its time reverse, the advanced time $v = ct + r_*$.
- `linear`, Kiselev's (22) with no charge and no de Sitter term, his example $w = -2/3$, $f = 1 - r_s/r - r/r_q$, the metric Cvetič, Gibbons and Pope write as their (6.5) and Fernando as her (7).
- `hyperbolic`, Kiselev's (25): the matter alone at $w = -2/3$, under the map he gives before it, $\chi = -\tfrac{1}{2}\ln(1 - r/r_q)$ and $t \to t/2r_q$.
- `conformally_flat`, Kiselev's (29), which his (27) and (28) lead to, "the Fock transformation" $\tau = e^{t}\cosh\chi$ and $\rho = e^{t}\sinh\chi$.

Kiselev keeps the letter $t$ for the rescaled time of his (25); it is $\eta$ here, since $t$ is the static time of the other charts, and every drawing is held to all the charts together.
The charts that keep $w$ name the term, $h = (r_q/r)^{3w + 1}$.
Every derivative of $h$ is $h$ times a rational function, $dh/dr = -(3w + 1)h/r$, so every value is a polynomial in $h$ over $r$, $r_s$ and $w$, and it is printed as one, with no power of $r$ that holds $w$ left standing.

## Step 2. The stress tensor

`kiselev_check` in `print_charts.py` holds every chart to Kiselev's stress tensor before it is written.
In the mixed components of this collection, with $G = 8\pi T$,

$$G^t{}_t = G^r{}_r = \frac{3w\,h}{r^2},\qquad G^\theta{}_\theta = G^\phi{}_\phi = -\frac{3w(3w + 1)\,h}{2r^2},$$

which is his (13) and (14), with $\rho_q = \tfrac{1}{2}c\,3w/r^{3(1 + w)}$ in his units $4\pi G = 1$ and with his $c = -r_q^{3w+1}$.
So the energy density is $\rho = -3w\,h/8\pi r^2$, positive for $w < 0$, the radial pressure is $p_r = -\rho$, and the tangential pressure is $p_t = \tfrac{1}{2}(3w + 1)\rho$.
Their ratio is $p_t/p_r = -(1 + 3w)/2$ and their average over directions is $(p_r + 2p_t)/3 = w\rho$, Visser's (3) to (5): the matter has two pressures for every $w$ but $-1$, and $w$ is the average.
The Ricci scalar is $R = 3w(3w - 1)h/r^2$, Kiselev's (17), and at $w = -2/3$ it is his $6/(r_q r)$.

The same function checks the static chart against three published metrics: Schwarzschild's of radius $r_s + r_q$ at $w = 0$, Kottler's with $\Lambda = 3/r_q^2$ at $w = -1$, and Reissner and Nordström's with $r_q^2 \to -r_q^2$ at $w = 1/3$.
Each Eddington-Finkelstein chart is the static one pulled back, the chart of the linear term is the static one at $w = -2/3$, the hyperbolic chart is that one at $r_s = 0$ pulled back along his map to it, and the conformally flat chart is the hyperbolic one pulled back along the inverse of his (27) and (28).

## Step 3. The example $w = -2/3$

With $f = 1 - r_s/r - r/r_q = -(r - r_-)(r - r_+)/(r_q r)$ the horizons are his (36) and (37),

$$r_\pm = \tfrac{1}{2}\left(r_q \pm \sqrt{r_q^2 - 4r_qr_s}\right),$$

real for $r_q \ge 4r_s$, with $f$ greatest at $r = \sqrt{r_sr_q}$.
The drawings take $r_q = 8\,r_s$, where $r_\pm = (4 \pm 2\sqrt{2})\,r_s$.
$1/f$ has no polynomial part, so

$$r_* = \sum_\pm \frac{\ln|1 - r/r_\pm|}{f'(r_\pm)},\qquad \sum_\pm \frac{1}{f'(r_\pm)} = -r_q,$$

which vanishes at $r = 0$ and falls as $-r_q\ln r$ as $r \to \infty$.
That is what sets this spacetime apart from Kottler's, whose $r_*$ tends to a constant: beyond $r_+$ both null edges of the region are $r = \infty$, and since $g_{tt}g_{rr} = -1$ the area radius is an affine parameter on every radial ray, Jacobson's result, so those edges are null infinity.
`conformal.KiselevMap` draws every region through Kottler's one function with the surface gravities $\kappa_\mp = \pm f'(r_\mp)/2$.

## Step 4. The matter alone

At $r_s = 0$ and $w = -2/3$, $f = 1 - r/r_q$, the one horizon is $r = r_q$ and $r_* = -r_q\ln|1 - r/r_q|$.
Kiselev's map is $\chi = r_*/2r_q$, so $u, v = 2r_q(\eta \mp \chi)$, and

$$ds^2 = 4r_q^2e^{-2\chi}\left(-d\eta^2 + d\chi^2 + \sinh^2\chi\,d\Omega^2\right) = \frac{4r_q^2}{(\tau + \rho)^2}\left(-d\tau^2 + d\rho^2 + \rho^2d\Omega^2\right),$$

with $\tau \pm \rho = e^{\eta \pm \chi}$ and $r = 2r_q\rho/(\tau + \rho)$.
The Weyl tensor vanishes, and the Ricci scalar is his $(3/r_q^2)(1 + \coth\chi) = (3/r_q^2)(\tau + \rho)/\rho$.
The conformally flat chart runs on past the horizon $\tau = \rho$ into $\tau < \rho$, where $r > r_q$, as far as $\tau + \rho = 0$, where $r \to \infty$: along a ray of constant $\tau - \rho$ the affine parameter is $r$ itself, so that edge is past null infinity, and the chart's own infinity $\tau + \rho \to \infty$ is the future horizon, reached at a finite affine parameter.

The static slice has $g_{rr} = r_q/(r_q - r)$, so its embedding climbs at $dz/dr = \sqrt{r/(r_q - r)}$, and

$$z = r_q\left(\arcsin\sqrt{r/r_q} - \sqrt{(r/r_q)(1 - r/r_q)}\right),$$

which with $r = \tfrac{1}{2}r_q(1 - \cos\psi)$ is $z = \tfrac{1}{2}r_q(\psi - \sin\psi)$: the cycloid of a circle of diameter $r_q$.
`embedding.kiselev` draws it by quadrature and checks it against this form.
