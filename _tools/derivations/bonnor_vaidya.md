# Bonnor-Vaidya

William Bonnor and Prahlad Vaidya's charged radiating star is Reissner-Nordström's metric with its mass and its charge functions of a null time,

$$ds^2 = -F\,dw^2 \pm 2\,dw\,dr + r^2\,d\Omega^2,\qquad F = 1 - \frac{2m(w)}{r} + \frac{q(w)^2}{r^2},$$

with the upper sign for the advanced time $w = v$ and the lower for the retarded time $w = u$.
Its three charts are written by `_tools/derivations/print_charts.py --metric bonnor_vaidya`, and `verify_metrics.py --system bonnor_vaidya/<chart>` checks each in seconds.

## Step 1. The parameters and the charts

The mass and the charge are lengths, $m = GM/c^2$ and $q^2 = GQ^2/(4\pi\epsilon_0c^4)$, and the null time is a length too, so no factor of $G$ or $c$ stands in a component.
They are Reissner-Nordström's lengths as the collection writes them, $r_s = 2m$ and $r_q = q$, and Vaidya's $Gm/c^2$.
The charge is written $q$ and never $e$, which the reader of the checker takes for the base of the exponential.
Each chart has a source:

- `eddington_finkelstein_outgoing`: Bonnor and Vaidya's own retarded coordinates, in which the star sends the dust out along the rays of constant $u$.
  Only the first two pages of their paper, General Relativity and Gravitation 1, 127 (1970), were read: its flat line element (2.1), $ds^2 = -r^2(d\theta^2 + \sin^2\theta\,d\phi^2) + 2\,du\,dr + du^2$ in the signature $(+,-,-,-)$, the potential (2.3) and the current (2.5).
  The curved line element is taken from Booth's (55), Physical Review D 93, 084005, and Kehle and Unger's (4.10), arXiv:2402.10190, which are the form published here.
- `eddington_finkelstein_ingoing`: Ori's (2) and (3), Classical and Quantum Gravity 8, 1559 (1991), Booth's (27), and Kehle and Unger's (4.5).
- `homothetic`: Koh, Park and Sherif's (29) to (35), JHEP 02 (2024) 028, for Lake and Zannias's homothetic case $m = \mu v$, $q = Qm/M$.
  With $R = Mr/m$ and $dv = (m/M)\,dV$, so that $v = (M/\mu)e^{\mu V/M}$,

  $$ds^2 = e^{2\mu V/M}\left(-\left(1 - \frac{2M}{R} + \frac{Q^2}{R^2} - \frac{2\mu R}{M}\right)dV^2 + 2\,dV\,dR + R^2\,d\Omega^2\right).$$

Not published as a chart: Chirenti and Saa's double null coordinates, arXiv:1012.5110, in which $r(u, v)$ has no closed form, and Lasky and Lun's coordinates of Painlevé and Gullstrand's kind, arXiv:0704.3634, whose coefficients are given by differential equations.
No diagonal chart with a varying mass and charge was found in the literature.

## Step 2. The field equations

`bonnor_vaidya_check` holds each Eddington-Finkelstein chart to the Einstein-Maxwell equations with a null current, in units $G = c = 4\pi\epsilon_0 = 1$.
The field is the Coulomb field of the charge the null time has reached, $A = \pm(q/r)\,dw$, whose divergence is a current along the rays, $4\pi J^r = \pm\partial_wq/r^2$, which is null.
The Einstein tensor is $8\pi$ times that field's Maxwell stress, $G^w{}_w = G^r{}_r = -q^2/r^4$ and $G^\theta{}_\theta = G^\phi{}_\phi = q^2/r^4$, plus the dust,

$$G_{ww}\Big|_{\text{dust}} = \pm\frac{2\left(r\,\partial_wm - q\,\partial_wq\right)}{r^3},$$

Ori's (6).
Its sign changes on $r_c = q\,\partial_wq/\partial_wm$, Ori's (11).
With $m$ and $q$ constant the chart is the published metric of `rn_metric` pulled back along $t = w \mp r_*$, and with $q = 0$ it is the published chart of the same name of `vaidya`.
The homothetic chart is checked to be the ingoing one pulled back, and its dust to be $G_{VV} = 2\mu(MR - Q^2)/(MR^3)$, which changes sign on $R = Q^2/M$.

## Step 3. The shell the drawings take

Every drawing but the homothetic one takes a thin shell of mass $M$ and charge $q = 24M/25$, the ratio of Reissner-Nordström's drawings, falling along $v = 0$ onto flat space, or its time reverse leaving along $u = 0$.
Outside the shell $F = (r - r_+)(r - r_-)/r^2$ with $r_+ = 32M/25$ and $r_- = 18M/25$, surface gravities $\kappa_\pm = (r_+ - r_-)/2r_\pm^2$, and

$$r_* = r + \frac{r_+^2\ln|r/r_+ - 1| - r_-^2\ln|r/r_- - 1|}{r_+ - r_-}.$$

Across the shell $G_{vv}$ is a delta function of $v$ times $\int 2(r\,\partial_vm - q\,\partial_vq)/r^3\,dv = 2(Mr - q^2/2)/r^3$, exactly, since $\int q\,\partial_vq\,dv = q^2/2$.
So the shell's energy density is positive outside $r_b = q^2/2M = 288M/625$ and negative inside it, and $r_b < r_-$ for every $q < M$.
On the sphere $r_b$ itself $F = 1$ on both sides of the shell.
Sympy's Heaviside function is $1/2$ at zero, which on the shell would make $q^2$ a quarter where $m$ is a half, so the rows of `null_rays.py` give both steps the flat side's value there.

## Step 4. The slices of the embedding

A slice of constant $v$ is null, so the moments are slices of constant $v - r = T$, on which the ingoing metric is $(1 + 2m/r - q^2/r^2)\,dr^2 + r^2d\phi^2$ on the equator.
Inside the shell that is a flat disc.
Outside it $dz/dr = \sqrt{2M/r - q^2/r^2}$, and with $s = \sqrt{2Mr - q^2}$

$$z = 2s - 2q\arctan(s/q),$$

which lies level on $r = q^2/2M$ and has no real height inside it.
The fold at the shell's radius $R = -T$ has the slope $\sqrt{2M/R - q^2/R^2}$, greatest at $R = q^2/M$ and zero at $R = r_b$, the last moment drawn.
Inside $r = \sqrt{M^2 + q^2} - M$ a slice of constant $v - r$ is timelike on the charged side, so the cones of the spacetime diagrams are oriented by the family that keeps the null time.

## Step 5. The conformal diagrams

The charged side is a `Tower` of the roots $r_\pm$, and the ingoing chart covers its cells I, II and III' with $q = \arctan e^{\kappa_+v}$ throughout.
On the flat side $u = v - 2r$, each outgoing ray keeps the $p$ it has where it crosses the shell at $r = -u/2$, $p = F(u)$, and $q = F(v) - \pi/2$.
The centre $u = v$ is then the line $q - p = -\pi/2$, on which the tower puts $r = 0$ of III' as well, and the two sides agree on the shell since $F(0) = 3\pi/4$.
Every $p$ is moved by $\pi/2$, which puts that line on $X = 0$.
The leaving shell is the same drawing under $(p, q) \to (-q, -p)$.

Ori's reading turns the shell round at $r_b$.
It leaves along the outgoing ray through that event, $p = p_b = F(-2r_b)$, through III', the reflected cell II and the reflected cell I, which the outgoing chart covers with the cell's own time $-u - r_*$; on it $u = u_s$ with $\arctan e^{-\kappa_+u_s} = \pi - p_b$.
The flat region is $v < 0$ together with $u > u_b = -2r_b$.
Inside the leaving shell an ingoing ray $v > 0$ takes the $q$ of the event where it meets the shell, at $r = (v - u_b)/2$, $q = Q(v)$, and $p = Q(u) + \pi/2$ for $u > 0$, which keeps the centre on $X = 0$ and agrees with $F$ at $u = 0$.
`conformal.py --verify` checks every map against the published metric of its chart, both shells to be one place from their two sides, and the centre to be straight.
Ori treats a continuous stream, for which the turning surface is $r_c(v)$; the thin shell is its limit, and his closing note cites Dray's and Barrabès and Israel's thin shells as reaching the same conclusion.
