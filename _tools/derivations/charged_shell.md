# The charged shell of dust

The five charts of `charged_shell.json` are written by `print_charts.py --metric charged_shell`, which computes every tensor from each line element with the checker's own `Geometry` and reads every printed value back before it writes.
This note records the source of each chart, the motion of the shell between them, and what the diagrams draw.

## Step 1. The spacetime

A thin spherical shell of dust of rest mass $m$ and charge $Q$, with flat spacetime inside it and Reissner-Nordström's field of mass $M$ and charge $Q$ outside (`delacruz1967`, `kuchar1968`, `boulware1973`, `gao2008`).
The shell encloses no charge, so there is no field inside it.
Three lengths stand for the two masses and the charge, $r_s = 2GM/c^2$, $\mu = Gm/c^2$ and $r_q$, with $r_q^2 = GQ^2/(4\pi\epsilon_0c^4)$, the lengths of `rn_metric.json`.
With $r_q < r_s/2$ the field outside has the horizons $r_\pm = r_s/2 \pm \sqrt{r_s^2/4 - r_q^2}$.

- `interior`: Minkowski's spherical chart, $-c^2dT^2 + dr^2 + r^2d\Omega^2$, with $T$ the time of an observer at rest at the centre, $r \le R$.
- `exterior`: Reissner-Nordström's static chart, $f = 1 - r_s/r + r_q^2/r^2$ (`reissner1916`, `nordstrom1918`), $r \ge R$. It ends on each horizon.
- `exterior_ingoing`: the ingoing Eddington-Finkelstein chart, $v = ct + r_*$ with $dr_*/dr = 1/f$, which follows the shell in through $r_+$ and $r_-$ to its turn and back up to $r_-$, where $v \to \infty$. It is Bonnor and Vaidya's ingoing chart with a constant mass and charge (`bonnor1970`, `ori1991charged`).
- `exterior_outgoing`: the outgoing chart, $u = ct - r_*$, which begins on the $r_-$ the shell went in through, at $u \to -\infty$, and follows it from the turn out through $r_-$ and $r_+$ into the next exterior.
- `exterior_isotropic`: the isotropic chart of Arnowitt, Deser and Misner (`arnowitt1960`, `arnowitt1962` section 7), in which the spatial metric is $\chi^4\delta_{ij}$ with $\chi^2 = \psi^2 - \phi^2$, $\psi = 1 + m/32\pi r$ and $\phi = e/16\pi r$ in their units. With $a = r_s/4 + r_q/2$ and $b = r_s/4 - r_q/2$ that is $\chi^2 = (1 + a/\rho)(1 + b/\rho)$, the areal radius is $r = (\rho + a)(\rho + b)/\rho$, and the lapse is $(\rho^2 - ab)/((\rho + a)(\rho + b))$. The chart takes $a$ and $b$ for its parameters, since every value is short in them and long in $r_s$ and $r_q$. At $b = 0$ it is the single hole of `majumdar_papapetrou.json` with $m = a$.

`charged_shell_check` holds the interior to the published spherical chart of `minkowski.json`, the static chart to the published metric of `rn_metric.json`, slot by slot, every exterior chart to $G^a{}_b = (r_q^2/r^4)\,\mathrm{diag}(-1, -1, 1, 1)$, to a vanishing Ricci scalar and to being the static chart pulled back, and the isotropic chart at $b = 0$ to Majumdar and Papapetrou's published single hole.

## Step 2. Why no chart runs through the shell

As for `israel_shell.md` Step 2: on the shell $dT/d\tau$ and $dt/d\tau$ differ and $g_{rr}$ jumps from $1$ to $1/f$.
Each chart declares $R$, or $\epsilon$ in the isotropic chart, as a function of its own time, and `REGIONS` and `OVERRULED` in `metric_tags.py` say for the tags that the two regions are one spacetime that is neither static nor stationary.

## Step 3. The motion

With a dot for $d/d(c\tau)$ along the shell, the junction conditions for dust, $S_{ab} = \sigma u_au_b$ with $4\pi R^2\sigma = m$, are those of `israel_shell.md` Step 3 with $f$ now Reissner-Nordström's:

$$\sqrt{1 + \dot{R}^2} - \sqrt{f + \dot{R}^2} = \frac{\mu}{R}, \qquad \frac{r_s}{2} = \mu\sqrt{1 + \dot{R}^2} + \frac{r_q^2 - \mu^2}{2R},$$

which is $M = m\sqrt{1 + \dot{R}^2} + (Q^2 - m^2)/2R$ in $G = c = 1$, Kuchař's law of conservation of energy (`kuchar1968`) and Gao and Lemos's (44) (`gao2008`).
So

$$\gamma = \frac{dT}{d\tau} = \frac{r_s}{2\mu} - \frac{r_q^2 - \mu^2}{2\mu R}, \qquad \beta = f\,\frac{dt}{d\tau} = \frac{r_s}{2\mu} - \frac{r_q^2 + \mu^2}{2\mu R}, \qquad \dot{R}^2 = \gamma^2 - 1 = \beta^2 - f,$$

and $dv/d(c\tau) = 1/(\beta - \dot{R})$, $du/d(c\tau) = 1/(\beta + \dot{R})$.
The shell is at rest where $\gamma = 1$, at $R = (r_q^2 - \mu^2)/(r_s - 2\mu)$, Gao and Lemos's (46).
For $r_q < r_s/2$ and $\mu < r_q$ that radius is at most $r_-$, with equality at $\mu = r_-$: its derivative along $\mu$ vanishes at $\mu = r_\pm$, and at $\mu = r_-$ it is $r_-(r_+ - r_-)/(r_+ + r_- - 2r_-) = r_-$.
$\beta$ at the turn is $1 - \mu/R$, positive for $\mu < r_-$, so the shell turns in the region inside $r_-$ that the ingoing chart reaches.

At rest, $\gamma = 1$ and $\beta = \sqrt{f}$, the first integral is $M = m + (Q^2 - m^2)/2R$.
In the isotropic radius $\epsilon$ of the shell, $R = (\epsilon + a)(\epsilon + b)/\epsilon$ and $R\sqrt{f} = \epsilon - ab/\epsilon$, it reads $\mu = a + b + 2ab/\epsilon$, which is $M = -\epsilon + \sqrt{\epsilon^2 + 2m\epsilon + Q^2}$, Arnowitt, Deser and Misner's (7.1): their bare mass is the rest mass of the dust and their $\epsilon$ the isotropic radius.
With $\mu = r_q = r_s/2$ the shell is at rest at every radius outside $r = a$, Gao and Lemos's case (iib).

`_tools/test_charged_shell.py` takes $u^a$, $n_a$ and the acceleration from the published metric components and Christoffel symbols of each chart and holds the unit speed and both jumps at five points, outside $r_+$, between the horizons, inside $r_-$, and for a field with no horizon, and holds the isotropic chart to the relation at rest.

## Step 4. The shell drawn

Every diagram but the isotropic chart's draws $r_s = 1$, $r_q = 12/25$ and $\mu = 1/5$: $r_+ = 16/25$, $r_- = 9/25$, $E = r_s/2\mu = 5/2$, $k = (r_q^2 - \mu^2)/2\mu = 119/250$ and $a = E^2 - 1 = 21/4$.
Inside, the motion is that of a charge in a repulsive Coulomb field in flat space, $\gamma = E - k/R$, and in the parameter $\eta$, zero at the turn,

$$R = \frac{k(E + \cosh\eta)}{a}, \qquad c\tau = \frac{k(E\eta + \sinh\eta)}{a^{3/2}}, \qquad cT = \frac{k(\eta + E\sinh\eta)}{a^{3/2}}.$$

The turn is at $R = k/(E - 1) = 119/375$, and the shell is on $r_+$ at $\eta = \mp 2.198$ and on $r_-$ at $\eta = \mp 0.936$.
The advanced time is one quadrature, $dv/d\eta = (R/\sqrt{a})/(\beta - \dot{R})$ over $\eta \le 0$, held in `charged_shell.py` as a table summed from the turn outward, with the static time zero at the turn, $v(0) = r_*(R)$ and $r_* = r + A_+\ln|r/r_+ - 1| + A_-\ln|r/r_- - 1|$, $A_\pm = \pm r_\pm^2/(r_+ - r_-)$.
The way out is the way in turned over in time, $u(\eta) = -v(-\eta)$, and $v = u + 2r_*$ there, which runs off to infinity where the shell comes back to $r_-$.

The numbers the captions state: the shell crosses $r_+$ at $cT = -0.527$, $v = -0.022$ and $v - r = -0.662$, and $r_-$ at $cT = -0.144$, $v = 0.106$ and $v - r = -0.254$; it turns at $v = 0.303$, $v - r = -0.015$; its proper time from $r_+$ to the turn is $0.393$ and from $r_-$ to the turn $0.135$; it comes in and leaves at $\sqrt{a}/E = 0.917$ of the speed of light.
The event horizon is $r = r_+$ outside the shell and the outgoing ray $cT - r = -1.167$ inside it.
The Cauchy horizon is the branch of $r_-$ at $v \to \infty$ outside the shell and the ingoing ray $cT + r = 0.504$ inside it.

The isotropic chart's diagrams draw the balanced shell, $b = 0$ and $\mu = a = 1$, at rest at $\epsilon = 1$, where its areal radius is $2$ and the lapse $\epsilon/(\epsilon + a) = 1/2$.
The embedding diagrams mark the circle of each shell a thousandth of the unit inside the fold, since on the fold itself whether the wall in front hides it hangs on the last digit kept, and `turn.js` and the script then disagree.

## Step 5. The diagrams

- Spacetime diagrams, seven views: the plane of $T$ and $r$ inside the shell and the line through the centre, with the shell, the event horizon and the Cauchy horizon; the static chart outside $r_+$ and inside $r_-$, where the shell turns at $t = 0$; the ingoing chart against $v - r$ and the outgoing chart against $u + r$; and the isotropic chart's plane of $t$ and $\rho$ with the balanced shell at rest. `--verify` checks the rays against $cT \pm r$, $ct \pm r_*$, $v$ and $v - 2r_*$, $u + 2r_*$ and $u$, and $ct \pm (\rho + 2\ln\rho - 1/\rho)$.
- Conformal diagram, five views. Inside the falling shell, $p, q = \arctan((cT \mp r)/L)$, Minkowski's compactification with $L = 0.504$, the shell's $cT + R$ where it comes back to $r_-$. Outside, an outgoing ray keeps the $p$ it has where it left the shell and an ingoing ray the $q$ it has where it meets it; every ray outside the shell crosses it once, so the whole of the outside lies in Minkowski's triangle, with $r_+$ on $p = -1.163$, the inner horizon's two branches on $p = -\pi/4$ and $q = \pi/4$, and the next $r_+$ on $q = 1.163$. A ray is traced to the shell by the quantity it keeps, $v$ or $v - 2r_*$, on the stretch of the shell's way it left in. The region inside $r_-$ beyond the Cauchy horizon, which neither null chart covers, is entered through its own static time by the symmetry of the whole about the turn. The balanced shell is Minkowski's triangle by the tortoise coordinate, as the gravastar's is.
- Embedding diagram, two views. The falling shell on Vaidya's slices: outside the shell $v - r = w$, on which the metric is $(1 + r_s/r - r_q^2/r^2)dr^2 + r^2d\phi^2$ and $z = 2s - 2r_q\arctan(s/r_q)$ with $s = \sqrt{r_sr - r_q^2}$; inside, the moment of $T$ at which that slice meets the shell, a flat disc. A movie of 81 frames from $w = -1.5$ to $1$. The balanced shell at rest: the slice $t = 0$ of the isotropic chart, $z = 2s + \ln((s - 1)/(s + 1))$ with $s = \sqrt{2r/a - 1}$, cut off at the shell by a flat disc, with the throat below the shell drawn on as a reference.
