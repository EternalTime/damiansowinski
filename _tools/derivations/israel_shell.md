# Israel's collapsing shell of dust

The three charts of `israel_shell.json` are written by `print_charts.py --metric israel_shell`, which computes every tensor from each line element with the checker's own `Geometry` and reads every printed value back before it writes.
This note records the source of each chart, the motion of the shell between them, and what the diagrams draw.

## Step 1. The spacetime

A thin spherical shell of dust of rest mass $m$, with flat spacetime inside it and Schwarzschild's vacuum of mass $M$ outside (`israel1966`, `poisson2004`).
The areal radius $r$ is one coordinate on both sides, and the shell is $r = R$.
Two lengths stand for the two masses, $r_s = 2GM/c^2$ and $\mu = Gm/c^2$.

- `interior`: Minkowski's spherical chart, $-c^2dT^2 + dr^2 + r^2d\Omega^2$, with $T$ the time of an observer at rest at the centre, $r \le R$.
- `exterior`: Schwarzschild's chart, $r \ge R$ (`schwarzschild1916`). It ends on $r_s$, which the shell reaches only as $t \to \infty$.
- `exterior_ingoing`: the ingoing Eddington-Finkelstein chart, $v = ct + r + r_s\ln|r/r_s - 1|$ (`eddington1924`, `finkelstein1958`), which follows the shell through $r_s$ to $r = 0$.

`israel_shell_check` holds the interior to a vanishing Riemann tensor and to the published spherical chart of `minkowski.json`, each exterior chart to a vacuum, to $K = 12r_s^2/r^6$ and to the published chart of `schwarzschild.json` it names, slot by slot, and the ingoing chart to being the exterior pulled back.

## Step 2. Why no chart runs through the shell

The time of one side is not the time of the other: on the shell $dT/d\tau$ and $dt/d\tau$ differ, and $g_{rr}$ jumps from $1$ to $1/(1 - r_s/R)$.
A chart with a continuous metric needs the proper distance from the shell and the shell's own proper time, in which neither side is in closed form.
So the shell is in no chart, as the gravastar's is not, each chart declares $R$ as a function of its own time, and the motion is derived in Step 3 and held to the published metrics by `_tools/test_israel_shell.py`.
`REGIONS` and `OVERRULED` in `metric_tags.py` say so for the tags: both regions are empty and static, and the spacetime is neither a vacuum nor static.

## Step 3. The motion

With a dot for $d/d(c\tau)$ along the shell and $f = 1 - r_s/R$, the unit normal has $n^r = \sqrt{1 + \dot{R}^2}$ inside and $\pm\sqrt{f + \dot{R}^2}$ outside, and

$$K^\theta{}_\theta = \frac{n^r}{R}, \qquad K^\tau{}_\tau = n_a a^a = \frac{d n^r}{dR}.$$

For dust, $S_{ab} = \sigma u_au_b$, the junction conditions are $[[K^\theta{}_\theta]] = -4\pi G\sigma/c^2$ and $[[K^\tau{}_\tau]] = +4\pi G\sigma/c^2$.
The two together give $4\pi R^2\sigma = m$, a constant, and

$$\sqrt{1 + \dot{R}^2} - \sqrt{1 - r_s/R + \dot{R}^2} = \frac{\mu}{R}, \qquad \frac{r_s}{2} = \mu\sqrt{1 + \dot{R}^2} - \frac{\mu^2}{2R},$$

which is $M = m\sqrt{1 + \dot{R}^2} - m^2/2R$ in $G = c = 1$.
With $E = r_s/2\mu = M/m$,

$$\gamma = \frac{dT}{d\tau} = E + \frac{\mu}{2R}, \qquad \beta = f\,\frac{dt}{d\tau} = E - \frac{\mu}{2R}, \qquad \dot{R}^2 = \gamma^2 - 1 = \beta^2 - f, \qquad \frac{dv}{d(c\tau)} = \frac{1}{\beta - \dot{R}}.$$

$\beta$ changes sign at $R = \mu/2E$, inside $r_s$, where the outward normal turns to point down Kruskal's $U$.
A shell with $E < 1$ is at rest at $R = \mu^2/(2\mu - r_s)$ and one with $E = 1$ falls from rest at infinity.

`_tools/test_israel_shell.py` takes $u^a$, $n_a$ and the acceleration from the published metric components and Christoffel symbols of each chart and holds the unit speed, both jumps and the first integral at five points, two of them inside $r_s$.

## Step 4. The shell drawn

Every diagram draws $E = 1$, $\mu = r_s/2$, in units of $r_s$, where $\dot{R}^2 = (8R + 1)/16R^2$.
In $s = \sqrt{1 + 8R}$ every time along the shell is elementary, and `israel_shell.py` holds them:

$$R = \frac{s^2 - 1}{8}, \quad c\tau = -\frac{(s - 1)^2(s + 2)}{24}, \quad cT = -\frac{s^3 + 3s - 4}{24}, \quad v = -\frac{1}{8}\left(\frac{s^3}{3} - s^2 + 5s - 15 - 16\ln\frac{s + 3}{6}\right),$$

since $\dot{R} = -2s/(s^2 - 1)$, $\gamma = (s^2 + 1)/(s^2 - 1)$, $\beta = (s^2 - 3)/(s^2 - 1)$, $f = (s^2 - 9)/(s^2 - 1)$ and $dv/d(c\tau) = (s + 1)/(s + 3)$.
The flat null times on the shell are $cT + R = -(s - 1)^3/24$ and $cT - R = -((s + 1)^3 - 8)/24$, and $R(T)$ is Cardano's root $s = 2\sinh(\tfrac{1}{3}\mathrm{arsinh}(2 - 12cT))$.

The shell crosses $r_s$ at $s = 3$: $c\tau = -5/6$, $cT = -4/3$, $v = 0$, $\dot{R} = -3/4$.
It reaches $R = 0$ at $s = 1$: $\tau = T = 0$ and $v = (32/3 - 16\ln\tfrac{3}{2})/8 = 0.5224$.
The event horizon is $r = r_s$ outside the shell and the outgoing ray $cT - r = -7/3$ inside it, so it leaves the centre at $cT = -7/3$, before the shell arrives.

## Step 5. The diagrams

- Spacetime diagrams: the plane of $T$ and $r$ inside the shell and the line through the centre, with the shell and the event horizon; Schwarzschild's plane outside it; and the ingoing chart against $v - r$, where the shell runs through $r_s$ to $r = 0$ and the singular edge begins at $v - r = 0.52$. A row's `curves` draws the shell, a world line $r = R(x^0)$ the row declares. `--verify` checks the rays against $cT \pm r$, $ct \pm r_*$, and $v$ and $v - 2r_*$.
- Conformal diagram: inside, $p, q = \arctan(3(cT \mp r)/7)$, Minkowski's compactification with the length $7/3$, so that the centre is $X = 0$ and the horizon $p = -\pi/4$. Outside, an outgoing ray keeps the $p$ it has where it left the shell and an ingoing ray with $v < 0.5224$ the $q$ it has where it enters, which makes both continuous on the shell; an ingoing ray with $v \ge 0.5224$ ends on the singularity and takes $q = -p$ of the outgoing ray that ends there with it, so the singularity is $T = 0$. Near the last ingoing ray to meet the shell $q$ goes as $-(0.5224 - v)^{3/2}$ below it and as $(v - 0.5224)^{1/2}$ above it, so the lines outside bend there.
- Embedding diagram: Vaidya's slices. Outside the shell $v - r = w$, on which the metric is $(1 + r_s/r)dr^2 + r^2d\phi^2$ and $z = 2\sqrt{r_sr}$; inside, the moment of $T$ at which that slice meets the shell, a flat disc. The fold at the shell has slope $\sqrt{r_s/R}$. A movie of 97 frames from $w = -5$ to $1$.
