# The gravastar

The three charts of `gravastar.json` are written by `print_charts.py --metric gravastar`, which computes every tensor from each line element with the checker's own `Geometry` and reads every printed value back before it writes.
This note records the source of each chart, what the shell between them carries, and what the diagrams draw.

## Step 1. The spacetime

Pawel Mazur and Emil Mottola write a static, spherically symmetric collapsed star as their (1),

$$ds^2 = -f(r)\,c^2dt^2 + \frac{dr^2}{h(r)} + r^2\,d\Omega^2,$$

with de Sitter space inside, their (8), $f = C\,h = C(1 - H_0^2r^2)$, and Schwarzschild's vacuum outside, their (9), $f = h = 1 - r_s/r$ (`mazur2004`).
Between the two they put a layer of fluid with $p = \rho$ held in by two thin shells.
Matt Visser and David Wiltshire merged that layer and its two shells into one thin shell at a radius $a$ outside $r_s$, which leaves three layers: de Sitter space, a shell, Schwarzschild (`visser2004`, section 2).
That three layer star is the spacetime published, with the shell at $r = R$, $L = 1/H_0$ the radius of the de Sitter space, and

$$r_s < R < L .$$

The first inequality keeps the shell outside Schwarzschild's horizon and the second inside de Sitter's, so $f > 0$ and $h > 0$ at every radius and $\partial_t$ is timelike everywhere.

## Step 2. The constant C

The induced metric on the shell is $-f(R)\,c^2dt^2 + R^2d\Omega^2$ from either side, and it has to be one metric, so $f$ is continuous there:

$$C\left(1 - \frac{R^2}{L^2}\right) = 1 - \frac{r_s}{R}, \qquad C = \frac{1 - r_s/R}{1 - R^2/L^2}.$$

Mazur and Mottola call $C$ "an arbitrary constant, corresponding to the freedom to redefine the interior time coordinate", and Pani, Berti, Cardoso, Chen and Norte fix it by this condition, their $\alpha$ (`pani2009`, section II).
With it one coordinate $t$ serves both sides.
The interior charts declare $C$ as a name for that expression, so every published value that writes $C$ is read as the expression, and `gravastar` in `print_charts.py` prints with a `pretty` that writes $r_s$ back as $R - CR(L^2 - R^2)/L^2$, after which $R$ cancels: the interior knows of the shell only through $C$.
Visser and Wiltshire write both sides with $f = h = 1 - 2m(r)/r$, their (3), each side with a time of its own; theirs inside is $\sqrt{C}\,t$.

## Step 3. The three charts

- `interior`: Mazur and Mottola's (1) with their (8), $r \le R$. It is the static patch of de Sitter space with its time rescaled, so $R_{\mu\nu} = (3/L^2)g_{\mu\nu}$, the Weyl tensor vanishes, $R = 12/L^2$ and $K = 24/L^4$, and only $g_{tt}$, $g^{tt}$, $\Gamma^r{}_{tt}$ and the components with a lowered pair of $t$ carry $C$.
- `interior_tortoise`: the tortoise coordinate Pani and his collaborators use for the perturbations of the interior, $dr/dx = \sqrt{fh} = \sqrt{C}\,(1 - r^2/L^2)$, so $x = (L/\sqrt{C})\,\mathrm{artanh}(r/L)$ and $r = L\tanh(\sqrt{C}\,x/L)$ (`pani2009`, section III A, where $C = 1$ and the coordinate is $r_*$). Then $f = C/\cosh^2(\sqrt{C}\,x/L)$ and $dr^2/h = f\,dx^2$, which is the line element published, conformally flat on the plane of $t$ and $x$. The shell is at $x_R = (L/2\sqrt{C})\ln((L + R)/(L - R))$.
- `exterior`: Schwarzschild's chart from the shell out, Mazur and Mottola's (9).

`gravastar_check` holds the two interior charts to $R_{\mu\nu} = (3/L^2)g_{\mu\nu}$ and to a vanishing Weyl tensor, the tortoise chart to being the interior chart pulled back along $r = L\tanh(\sqrt{C}\,x/L)$, the exterior to a vanishing Ricci tensor and to being the published metric of `schwarzschild.json` slot by slot, and $g_{tt}$ on the shell to $-(1 - r_s/R)$ from both sides.
The exterior's tortoise coordinate $r + r_s\ln(r/r_s - 1)$ has no inverse in closed form, so the exterior has no tortoise chart.

## Step 4. Why no chart runs through the shell

In the areal radius $g_{rr} = 1/h$ jumps at $r = R$, from $1/(1 - R^2/L^2)$ to $1/(1 - r_s/R)$, so the Christoffel symbols of a chart through the shell hold a delta times a function that jumps, which no reading as a distribution fixes; `_on_a_kink` in `verify_metrics.py` stops on exactly that.
A chart with a continuous metric needs the proper distance from the shell, Visser and Wiltshire's Gaussian normal coordinate $\eta$ with $d\eta = dr/\sqrt{h}$, their (9), or the tortoise coordinate, and neither has $r$ in closed form outside.
So the shell is in no chart, the two sides are published each in its own, as the collapsing star of Oppenheimer and Snyder is, and what the shell carries is derived in Step 5 and held to the published metrics by `_tools/test_gravastar.py`.

## Step 5. The shell

With the unit normal pointing outward, the extrinsic curvature of the sphere $r = R$ in a chart of this form is

$$K^t{}_t = \frac{\sqrt{h}}{2f}\frac{df}{dr}, \qquad K^\theta{}_\theta = K^\phi{}_\phi = \frac{\sqrt{h}}{r}.$$

A shell of surface density $\sigma$ and surface tension $\vartheta$ has, by the junction conditions of Israel (`israel1966`) as Visser and Wiltshire write them for a static shell, their (12) and (13), with $[[X]]$ the value outside less the value inside and $G = c = 1$,

$$\left[\left[\frac{\sqrt{h}}{R}\right]\right] = -4\pi\sigma, \qquad \left[\left[\frac{1 - m/R - m'}{R\sqrt{1 - 2m/R}}\right]\right] = -8\pi\vartheta,$$

where $h = 1 - 2m(r)/r$, so $m = r_s/2$ outside and $m = r^3/2L^2$ inside, and the second bracket is $K^t{}_t + K^\theta{}_\theta$.
With the units restored,

$$\sigma = \frac{c^2}{4\pi GR}\left(\sqrt{1 - \frac{R^2}{L^2}} - \sqrt{1 - \frac{r_s}{R}}\right), \qquad \vartheta = \frac{c^4}{8\pi GR}\left(\frac{1 - 2R^2/L^2}{\sqrt{1 - R^2/L^2}} - \frac{1 - r_s/2R}{\sqrt{1 - r_s/R}}\right).$$

Three things follow.

- $\sigma = 0$ exactly when $L^2 = R^3/r_s$, which is $\tfrac{4}{3}\pi R^3\rho = Mc^2$: the vacuum inside accounts for the whole mass. Then $h$ is continuous as well as $f$, $C = 1$, and the shell carries a tension alone, which is negative, a pressure. That is the thin shell gravastar of Pani and his collaborators and of the papers on echoes (`pani2009`, `cardoso2016prd`).
- For $L^2 > R^3/r_s$ the vacuum accounts for less than the mass and $\sigma > 0$. The diagrams draw such a star, $R = 5r_s/4$ and $L = 2r_s$, where $C = 64/195$ and $4\pi G R\sigma/c^2 = 0.33$.
- As $R \to r_s$ the second term of $\vartheta$ diverges, which is why Visser and Wiltshire advise a shell at a finite distance outside $r_s$.

Mazur and Mottola's star of 2015 lies on the edge of the family, $R = L = r_s$, where $f$ vanishes on the shell from both sides and continuity no longer fixes $C$: Schwarzschild's star of uniform density at $R = r_s$ has $-g_{tt} = \tfrac{1}{4}(1 - r^2/r_s^2)$ and $1/g_{rr} = 1 - r^2/r_s^2$, so $C = 1/4$ (`mazur2015`).

`_tools/test_gravastar.py` takes $K^t{}_t$ and $K^\theta{}_\theta$ from the published $g_{tt}$ and $g_{rr}$ of the two charts and holds the two closed forms above, the three statements about $\sigma$, and the limit $C = 1/4$ against the published metric of `interior_schwarzschild.json`.

## Step 6. The diagrams

Every diagram draws the star of Step 5 with $r_s = 1$, $R = 5/4$ and $L = 2$.

- Spacetime diagrams: the plane of $t$ and $r$ inside the shell and the line through the centre, where $dr/d(ct) = \pm\sqrt{C}(1 - r^2/L^2)$, $0.57$ at the centre and $0.35$ on the shell; the plane of $t$ and $x$, where every ray runs at 45° and the shell is at $x_R = 2.56\,r_s$; and Schwarzschild's plane from the shell out, where $dr/d(ct) = \pm(1 - r_s/r)$ is $\pm 0.2$ on the shell. The slope of a ray jumps across the shell, from $0.35$ to $0.2$, since $\sqrt{fh}$ does. `--verify` checks the rays against $ct \mp x$ inside and $ct \mp r_*$ outside.
- Conformal diagram: $x$ runs on through the shell, $x = x_R + r_* - r_*(R)$ outside, and $p, q = \arctan((ct \mp x)/x_R)$ give Minkowski's triangle with the shell a timelike line through $(X, T) = (\pi/2, 0)$, one view for each chart.
- Embedding diagram: the equator at one moment of $t$ is a cap of the sphere of radius $L$ inside the shell and Flamm's paraboloid outside. The cap reaches the shell at $dz/dr = (R/L)/\sqrt{1 - R^2/L^2} = 0.80$ and the paraboloid leaves it at $\sqrt{r_s/(R - r_s)} = 2$, so they meet at a crease, and since $\sqrt{h} = 1/\sqrt{1 + (dz/dr)^2}$ the crease is the jump of $\sqrt{h}$ that $\sigma$ measures. The geometry is static, so it is one surface.
