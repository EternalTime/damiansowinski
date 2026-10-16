# Kodama's kink

Kodama's kink is a static, spherically symmetric wormhole held open by a ghost scalar field with a quartic potential.
No closed form is known for its metric.
The entry therefore gives two charts with their metric functions left free, as `boson_star` does, and every drawing is made from one numerical solution, which `_tools/derivations/kodama_kink.py` computes and `_tools/test_kodama_kink.py` holds to the physics.

The charts are written by `_tools/derivations/print_charts.py --metric kodama_kink`, and `verify_metrics.py --system kodama_kink/<chart>` checks each in about a second.

## Step 1: the sources

The published paper, T. Kodama, Phys. Rev. D 18, 3529 (1978), doi:10.1103/PhysRevD.18.3529, was reached only through its abstract.
Its preprint, CBPF-NF-003/78 of January 1978, was read in full from the scan the Centro Brasileiro de Pesquisas Físicas keeps at cbpfindex.cbpf.br, and so was the preprint of T. Kodama, L. C. S. de Oliveira and F. C. Santos, CBPF-NF-036/78, published as Phys. Rev. D 19, 3576 (1979).
The equation numbers below are the preprints'.

- `areal`: Kodama's (2), $ds^2 = e^{2\eta}(dx^0)^2 - e^{2\alpha}dr^2 - r^2d\Omega^2$ in the signature $(+,-,-,-)$, here $ds^2 = -e^{2\eta}c^2dt^2 + e^{2\alpha}dr^2 + r^2d\Omega^2$, on either side of the throat $r = r_0$, where $e^{-2\alpha}$ vanishes.
- `proper_radial`: Kodama, de Oliveira and Santos's (2), $ds^2 = e^{2h}(dx^0)^2 - dR^2 - r^2d\Omega^2$, with $(dr/dR)^2 = e^{-2\alpha}$, their (3), here $ds^2 = -e^{2h}c^2dt^2 + dR^2 + r(R)^2d\Omega^2$, with $R$ over the whole line and the throat at $R = 0$.

## Step 2: the field equations

Kodama's Lagrangian, his (1), is $\kappa\mathcal{L} = \sqrt{-g}\left[\tfrac{1}{2}R + \epsilon\left(S_{,\alpha}S_{,\beta}g^{\alpha\beta} - V(S^2)\right)\right]$ with $\kappa = 8\pi G/c^4$, $\epsilon = -1$ for the ghost, and his (7), $V = \tfrac{1}{2}(\mu/f)^2(1 - f^2S^2)^2$.
In the signature $(-,+,+,+)$ the stress tensor is $\kappa T_{ab} = 2\epsilon\,\partial_aS\,\partial_bS - \epsilon g_{ab}\left(g^{cd}\partial_cS\,\partial_dS + V\right)$, so for a static $S$ and $\epsilon = -1$

$$G^t{}_t = e^{-2\alpha}(\partial_rS)^2 + V,\qquad G^r{}_r = V - e^{-2\alpha}(\partial_rS)^2,\qquad G^\theta{}_\theta = e^{-2\alpha}(\partial_rS)^2 + V,$$

and along the proper distance the same with $\alpha = 0$.
These are Kodama's (3) and (4), $\eta_1 = \epsilon rS_1^2 - \alpha_1$ and $2r\alpha_1 = \epsilon r^2S_1^2 + 1 - (1 - \epsilon r^2V)e^{2\alpha}$, as `kodama_kink_check` in `print_charts.py` confirms on the published $G^t{}_t$ and $G^r{}_r$.
The field obeys $\Box S = S\,dV/d(S^2)$, Kodama's (5), in which $\epsilon$ cancels:

$$\partial_R^2S + \left(\partial_Rh + \frac{2\,\partial_Rr}{r}\right)\partial_RS = -\mu^2S\left(1 - f^2S^2\right).$$

`kodama_kink_check` solves the first two for $\partial_r\eta$ and $\partial_r\alpha$ (or $\partial_Rh$ and $\partial_R^2r$), takes the wave equation for $\partial^2S$, and finds the published $G^\theta{}_\theta$ equal to the field's stress across the radius identically, which is the contracted Bianchi identity.
No Einstein component stands off the diagonal.
At $S = 1/f$, a vacuum, the areal chart with $e^{2\eta} = e^{-2\alpha} = 1 - r_s/r$ is the published metric of `schwarzschild`, and the proper chart is the areal one pulled back along $r = r(R)$ with $e^{-\alpha} = \partial_Rr$ and $\eta = h$.

Kodama, de Oliveira and Santos's (19) and (20), $\bar y'(0) = (2f^2 + x_0^2)^{1/2}/x_0$ and $\bar x''(0) = f^2(2f^2 + x_0^2)/2x_0$, disagree with Kodama's own (25) to (28) and with the equations above; they read as misprints of $f^{-2}$ for $f^2$ and a dropped $\sqrt2$.
At the throat, in $x = \mu r$ and $y = fS$, the $G^R{}_R$ equation gives $y'(0)^2 = \tfrac{1}{2} + f^2/x_0^2$ and the $G^t{}_t$ equation $x''(0) = (2f^2 + x_0^2)/(2x_0f^2)$, which are Kodama's (27) and (28) carried from his areal $\rho = x - x_0$ to the proper distance.

## Step 3: the numerical solution

`kodama_kink.Kink` works in units with $c = \mu = 1$.
It reads the published $G^t{}_t$ and $G^\theta{}_\theta$ of the proper chart from the metric file, which give $\partial_R^2r$ and $\partial_R^2h$ regularly at the throat, and leaves the published $G^R{}_R$, which there reads $0/0$ for $\partial_Rh$, for a check.

The solution is even in $R$ for $r$ and $h$ and odd for $S$, so it starts at the throat from $r(0) = r_0$, $r'(0) = h'(0) = S(0) = 0$ and the $S'(0)$ above.
The field reaches its vacuum $S = 1/f$ far away for one $r_0$ alone, the eigenvalue of Kodama's section 4, found by bisection: too small and the field turns back before it gets there, too large and it overshoots.
The field of the last radius that turns back is followed to its turning point, the solver's `edge`, at $\mu R = 13.2$ for the declared kink, where $1 - fS = 2\times10^{-9}$.
From there on the metric is Schwarzschild's of the mass $m = r(1 - r'^2)/2$, with $R - R_{\rm edge} = F(r) - F(r_{\rm edge})$, $F(r) = \sqrt{r(r - 2m)} + 2m\ln(\sqrt r + \sqrt{r - 2m})$, and $h$ is shifted to meet $\tfrac{1}{2}\ln(1 - 2m/r)$ there, so that $h \to 0$ far away.

The areal chart's $\eta$ and $\alpha$ are read from the proper chart at the same sphere, `kodama_kink.Areal`, by inverting $r(R)$ from a table in $\sqrt{r - r_0}$, since $r - r_0$ grows as $R^2$ from the throat.

## Step 4: the declared kink and the checks

Every drawing declares $f = 1.25$, the value of all three figures of Kodama, de Oliveira and Santos and of the third curve of Kodama's Fig. 1:

| | |
| --- | --- |
| throat radius $r_0$ | $1.475/\mu$ |
| mass $GM/c^2$ on either side | $0.390/\mu$ |
| $e^h$ at the throat | $0.703$ |
| field's edge | $\mu R = 13.2$ |

Kodama's Fig. 1 has the throats of $f = 0.77$, $1$ and $1.25$ at about $\mu r_0 = 0.25$, $0.78$ and $1.47$; the solver gives $0.257$, $0.776$ and $1.475$.
His and his colleagues' Fig. 1c has $h \approx -0.36$ at the throat of $f = 1.25$; the solver gives $-0.353$.
Below $f \approx 0.645$ he found no solution and could not confirm it; the solver follows the throat down to $\mu r_0 = 0.006$ at $f = 0.64$, where it closes linearly toward $f \approx 0.634$.

For large $f$ he quotes $x_0 \approx 1.30f^2$; the solver gives $r_0/f^2 = 1.2245$ at $f = 40$, and the limit is $\sqrt6/2 = 1.22474$.
In that limit the kink is a wall of width $1/\mu$ on a throat of radius about $f^2/\mu$, with surface energy $-\tfrac{4\sqrt2}{3}/f^2$ and tangential pressure the opposite of it, the flat kink $y = \tanh(u/\sqrt2)$ integrated across.
Israel's junction conditions for a wall joining two Schwarzschild exteriors of mass $m$ at radius $a$, with $8\pi\sigma = -(4/a)\sqrt{1 - 2m/a}$ and $8\pi p = (2/a)(1 - m/a)/\sqrt{1 - 2m/a}$, give $p = -\sigma$ only at $a = 3m$, the photon sphere, and then $a = \sqrt6f^2/2$.
The solver gives $r_0/m = 3.002$ at $f = 20$.
So Kodama's $1.30$ was his integration's accuracy, and the History does not repeat it.

`_tools/test_kodama_kink.py` holds the solver to the published $G^R{}_R$, to $10^{-10}$ of the energy at the throat; to Kodama's figures and his throat values (25) to (28); to the thin wall; to parity; to Schwarzschild's metric outside the edge; and the areal chart to its own published $G^t{}_t$.
It holds the drawings to the solver from the numbers in their files: every ray of both charts is a null curve, $c\,dt = \pm e^{-h}\,dR$, the embedded surface climbs at $\sqrt{1 - (\partial_Rr)^2}$ with the throat its narrowest circle, and the conformal diagram lies in the diamond.

## Step 5: stability

Kodama, de Oliveira and Santos's radial perturbations, their case (b), hold $X(0) = 0$, the radius of the throat.
V. Dzhunushaliev, V. Folomeev, B. Kleihaus and J. Kunz, Phys. Rev. D 87, 104036 (2013), arXiv:1302.5217, let it move, and for their case without fluid, which "corresponds to the problem considered in" Kodama, de Oliveira and Santos, find the lowest $\omega^2$ "negative for any size of the throat".
K. A. Bronnikov, J. C. Fabris and A. Zhidenko, Eur. Phys. J. C 71, 1791 (2011), prove instability for $V = 0$ only, and V. Dzhunushaliev, V. Folomeev, B. Kleihaus and J. Kunz, arXiv:1710.01884, treat topologically trivial ghost fields "in contrast to" Kodama's.

## Step 6: dimensions

`DIMENSIONS` declares $t$ a time, $r$, $R$ and $r_0$ lengths, and $\eta$, $\alpha$ and $h$ pure numbers, since each sits in an exponential.
