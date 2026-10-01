# The dilaton black hole of Gibbons and Maeda and of Garfinkle, Horowitz and Strominger

The Einstein metric of the charged black hole of low energy string theory is

$$ds^2 = -\left(1 - \frac{r_s}{r}\right)c^2dt^2 + \frac{dr^2}{1 - \dfrac{r_s}{r}} + r\left(r - r_d\right)\left(d\theta^2 + \sin^2\theta\,d\phi^2\right),$$

with $r_s = 2GM/c^2$ and $r_d = 2r_q^2/r_s$, $r_q^2 = Q^2G/(4\pi\epsilon_0c^4)$, which is $Q^2/M$ in units with $G = c = 1$.
Its five charts are written by `_tools/derivations/print_charts.py --metric dilaton_black_hole` in about 8 seconds, and `verify_metrics.py --system dilaton_black_hole/<chart>` checks each.

## Step 1. The parameters

Garfinkle, Horowitz and Strominger write the metric with $2M$ and $Q^2/M$, with the dilaton's value at infinity set to zero.
The collection keeps Schwarzschild's $r_s$ for the first and names the second $r_d$, the radius at which the area $4\pi r(r - r_d)$ of the spheres vanishes, so every value reduces to Schwarzschild's at $r_d = 0$ term by term.
The Reissner-Nordström entry's charge radius $r_q$ fixes it, $r_d = 2r_q^2/r_s$, and the extremal charge $Q^2 = 2M^2$ is $r_d = r_s$.
Gibbons and Maeda's $r_+$ and $r_-$ are $r_s$ and $r_d$.

## Step 2. The charts

The static chart is the one above.
The two Eddington-Finkelstein charts are built on Schwarzschild's tortoise coordinate $r_* = r + r_s\ln|r/r_s - 1|$, since the plane of $t$ and $r$ is Schwarzschild's.
The two string charts print the string metric $e^{2\varphi}g_{\mu\nu}$, in the same $t$ and $r$: with $e^{-2\varphi} = 1 - r_d/r$ for the magnetically charged hole,

$$ds^2 = -\frac{1 - \dfrac{r_s}{r}}{1 - \dfrac{r_d}{r}}c^2dt^2 + \frac{dr^2}{\left(1 - \dfrac{r_s}{r}\right)\left(1 - \dfrac{r_d}{r}\right)} + r^2d\Omega^2,$$

and with $e^{2\varphi} = 1 - r_d/r$ for the electrically charged one,

$$ds^2 = \left(1 - \frac{r_d}{r}\right)\left(-\left(1 - \frac{r_s}{r}\right)c^2dt^2 + \frac{dr^2}{1 - \dfrac{r_s}{r}}\right) + \left(r - r_d\right)^2d\Omega^2.$$

They are other metrics on the same manifold, the ones strings move in, and the literature on the extremal holes is written in them, so they stand beside the Einstein charts.
At $r_d = r_s$ the first is $-c^2dt^2 + dr^2/(1 - r_s/r)^2 + r^2d\Omega^2$, the infinitely long throat, and the second $-(1 - r_s/r)^2c^2dt^2 + dr^2 + (r - r_s)^2d\Omega^2$, with flat moments.
The printer checks each string chart to be that power of $1 - r_d/r$ times the static chart.

## Step 3. The field equations

For the action $\int\sqrt{-g}\,[R - 2(\nabla\varphi)^2 - e^{-2\varphi}F^2]$ and the magnetically charged hole, $F = Q\sin\theta\,d\theta\wedge d\phi$ with $Q^2 = r_sr_d/2$, the printer checks in each of the three Einstein charts, slot by slot, that

$$G_{\mu\nu} = 2\partial_\mu\varphi\,\partial_\nu\varphi - g_{\mu\nu}(\nabla\varphi)^2 + 2e^{-2\varphi}\left(F_{\mu\alpha}F_\nu{}^\alpha - \tfrac14 g_{\mu\nu}F^2\right),$$

that $\partial_\mu(\sqrt{-g}\,e^{-2\varphi}F^{\mu\nu}) = 0$, and that $\Box\varphi = -\tfrac12 e^{-2\varphi}F^2$.
The Ricci scalar is $R = 2(\nabla\varphi)^2 = r_d^2(r - r_s)/(2r^3(r - r_d)^2)$, which vanishes on the horizon.

## Step 4. The diagrams

Every drawing is at $r_d = r_s/2$, the charge $Q = M$ at which the Reissner-Nordström black hole is extremal.
The spacetime diagrams start at the singularity $r = r_d$, their left edge; on the plane of $t$ and $r$ all three metrics have the null rays $ct \pm r_*$, which `--verify` checks.
The conformal diagram is Kruskal's square with $U$ and $V$ each divided by $\sqrt{k}$, $k = (1 - r_d/r_s)e^{r_d/r_s}$, which moves the tortoise coordinate by $r_s\ln k$ and puts the singularity $UV = k$ on the straight lines $T = \pm\pi/2$; `ShiftedTower` in `conformal.py` does it.
The embedding diagram has one view for each metric, the equator of the moment $t = 0$ through the bifurcation sphere into the second exterior: the Einstein metric's throat has the radius $\sqrt{r_s(r_s - r_d)}$, the magnetic string metric's the radius $r_s$, and the electric string metric's surface is Flamm's paraboloid $z^2 = 4(r_s - r_d)(r - r_s)$ in the areal radius $r - r_d$, which is checked.
