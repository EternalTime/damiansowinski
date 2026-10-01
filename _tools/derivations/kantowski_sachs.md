# The Kantowski-Sachs cosmologies

The three charts of `kantowski_sachs.json` are written by `print_charts.py --metric kantowski_sachs`, which computes every tensor from each line element with the checker's own `Geometry` and reads every printed value back before it writes.
This note derives what the charts are, why their tensors come out as printed, and what the diagrams draw.

## Step 1. The comoving chart

The group $\mathbb{R} \times SO(3)$ acts on moments $\mathbb{R} \times S^2$, translations along the axis $r$ and rotations of the sphere, so a metric it leaves alone has one scale factor for the axis and one for the spheres,

$$ds^2 = -c^2dt^2 + a(t)^2dr^2 + b(t)^2\left(d\theta^2 + \sin^2\theta\,d\phi^2\right).$$

$r$ is a length, so $a$ is a pure number and $b$ carries the length.
Every component is taken in the chart $x^0 = ct$, and a prime is $d/d(ct)$.
In the orthonormal frame $c\,dt$, $a\,dr$, $b\,d\theta$, $b\sin\theta\,d\phi$ the curvature has four independent components, $a''/a$ on the plane of $t$ and $r$, $b''/b$ on each plane of $t$ and an angle, $a'b'/(ab)$ on each plane of $r$ and an angle, and $(b'^2 + 1)/b^2$ on the sphere, the last holding the sphere's own curvature $1/b^2$.
The Ricci scalar is twice their sum with those multiplicities, and the Kretschmann scalar four times the sum of their squares with the same multiplicities, which is how both are printed.
The Einstein tensor is

$$G^t{}_t = -\frac{b'^2 + 1 + 2ba'b'/a}{b^2}, \qquad G^r{}_r = -\frac{2bb'' + b'^2 + 1}{b^2}, \qquad G^\theta{}_\theta = G^\phi{}_\phi = -\frac{ab'' + ba'' + a'b'}{ab}.$$

## Step 2. Dust

Dust at rest in the chart has $G^r{}_r = G^\theta{}_\theta = 0$.
The first is the equation of a cycloid for $b$ alone: with $b = b_0\cos^2\eta$ and $c\,dt = 2b_0\cos^2\eta\,d\eta$, so that $ct = b_0(\eta + \sin\eta\cos\eta)$, $b' = -\tan\eta$ and $2bb'' + b'^2 + 1 = 0$.
The second is then linear in $a$, and its two solutions are $\tan\eta$ and $1 + \eta\tan\eta$.
The first alone is the vacuum of Step 3, and the second carries the dust, so with the scale of $r$ chosen to make its coefficient one,

$$a = 1 + (\eta + \kappa)\tan\eta, \qquad G^\eta{}_\eta = -\frac{1}{b_0^2\,a\cos^4\eta} = -\frac{8\pi G\rho}{c^2},$$

and $\rho\,a\,b^2 = c^2/8\pi G$ is constant, as it must be for dust in a volume $a\,b^2$.
For $|\kappa| < \pi/2$ the product $a\cos\eta = \cos\eta + (\eta + \kappa)\sin\eta$ has its least value, $\cos\kappa$, at $\eta = -\kappa$ and is positive at both ends, so $a > 0$ on the whole interval and both singularities are cigarlike, $b \to 0$ with $a \to \infty$.
Every tensor of the dust chart is a rational function of $\tan\eta$ and $\eta + \kappa$, since $\cos^4\eta = 1/(1 + \tan^2\eta)^2$ and $d\tan\eta/d\eta = 1 + \tan^2\eta$, and `kantowski_sachs_dust` in `print_charts.py` factors each value there, where factoring is unique, and prints $1 + \tan^2\eta$ as $1/\cos^2\eta$ and the factor $1 + (\eta + \kappa)\tan\eta$ as the line element writes it.
`kantowski_sachs_member` checks, slot by slot, that the chart's metric is the comoving one with those $a$, $b$ and $c\,dt$.

## Step 3. The vacuum

With no dust, $a = \tan\eta$ and $b = r_s\cos^2\eta$.
Writing $T = r_s\cos^2\eta$ gives $a^2 = r_s/T - 1$ and $c\,dt = dT/a$, so

$$ds^2 = -\frac{dT^2}{r_s/T - 1} + \left(\frac{r_s}{T} - 1\right)dr^2 + T^2d\Omega^2, \qquad 0 < T < r_s,$$

which is Schwarzschild's line element inside the horizon with the areal radius called $T$ and $ct$ called $r$.
Its Kretschmann scalar is $12r_s^2/T^6$, finite at the horizon $T = r_s$, where $a = 0$, and divergent at $T = 0$.
The future lies toward smaller $T$ in the black hole and toward larger $T$ in the white hole.

## Step 4. The spacetime diagrams

`null_rays.py` draws three planes of time and $r$.
The comoving plane has $a$ and $b$ solved by `DustSolver` from the published $G^r{}_r = G^\theta{}_\theta = 0$, from rest at $a = 1$, $b = b_0$, which is $\kappa = 0$, and `--verify` checks the solution against Step 2's closed forms and the first singularity against $ct = -\pi b_0/2$.
The dust chart's plane is the same universe in $\eta$, whose rays keep $r \mp \tau$ with $\tau = \int_0^\eta 2\cos^2s\,ds/(1 + s\tan s)$, and the vacuum plane's rays keep $r \pm (T + r_s\ln(1 - T/r_s))$; `--verify` checks both.
$\tau$ runs only to $\pm 1.2189\,b_0$, so light crosses $2.44\,b_0$ of $r$ in the whole life of the dust universe.

## Step 5. The conformal diagram

On the plane of $T$ and $r$ the vacuum metric is $-F\,d(r)^2 + dT^2/F$ with $F = 1 - r_s/T$ and the roles of the two coordinates exchanged, so `Tower` draws it as cell II of Kruskal's square, the triangle between the horizon and the singularity.
The dust universe's plane is conformal to the strip $|\tau| < 1.2189\,b_0$ of Minkowski's plane, and Minkowski's maps $p, q = \arctan(\tau \mp r)$ send it to the lens between the curves $\tau = \pm 1.2189\,b_0$, which meet at the two ends of the axis.

## Step 6. The embedding diagram

The equator of a moment has the metric $a^2dr^2 + b^2d\phi^2$ with $a$ and $b$ constant on it, a flat cylinder of radius $b$ on which the stretch $|r| \le 1$ is $2a$ long.
The dust universe is played as a movie in $\eta$ from $-1.1$ to $1.1$, and the vacuum is set as five cylinders from $T = 0.9\,r_s$ to $0.1\,r_s$, in the order of the proper time $c\tau = r_s(\eta + \sin\eta\cos\eta)$ of an observer at fixed $r$ since the horizon.
The dust's moments are marked on the comoving and dust planes and on the lens, and the vacuum's on the plane of $T$ and $r$ and in Kruskal's square; `slices.py` checks that Step 2's $a$, $b$ and $t$ pull the comoving chart back onto the dust chart.
