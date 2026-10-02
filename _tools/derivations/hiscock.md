# Hiscock's evaporating black hole

The three charts of `hiscock.json` are written by `print_charts.py --metric hiscock`, which computes every tensor from each line element with the checker's own `Geometry` and reads every printed value back before it writes.
This note records the source of each chart, what the surface between the two Vaidya charts carries, and what the diagrams draw.

## Step 1. The spacetime and its sources

William Hiscock's two papers of 1981 are `hiscock1981` and `hiscock1981outgoing`.
The first builds evaporating holes from the ingoing Vaidya metric alone, the second joins the ingoing metric to the outgoing one along a timelike surface "near (but outside) the apparent horizon".
Both are behind the publisher's wall, so every statement about them here and in the History is taken from their abstracts, as the publisher and INSPIRE print them, and the construction from the papers that write it out:

- Christian Lübbe and Paul Tod, `lubbe2009`, section 2 and figure 1: the regions, the matching surfaces, the ingoing metric as their (1) and the outgoing metric as their (2).
- Pierre Martin-Dussaud and Carlo Rovelli, `martindussaud2019`, their (17) to (21) and Appendix A: five patches, flat space in double null coordinates, Schwarzschild's metric, the ingoing Vaidya metric with a mass $N(v)$, the outgoing Vaidya metric with a mass $M(u)$, and flat space again.
- Chang-Zhong Guo, Wen-Cong Gan and Fu-Wen Shu, `guo2023`, their (2.1) to (2.3): the simplest model, two shells of mass $m_0$ and $-m_0$.
- Sean Hayward, `hayward2006`: the surface of pair creation at a fixed radius, and the surface layer on it.

The units are Hiscock's, $G = c = 1$, with the mass a length $m = GM/c^2$ and both null times lengths.

## Step 2. The charts

The ingoing chart is $ds^2 = -(1 - 2m/r)\,dv^2 + 2\,dv\,dr + r^2d\Omega^2$ with $m = m(v)$.
It is flat space where $m = 0$, Schwarzschild's metric in ingoing Eddington-Finkelstein coordinates where $m$ is constant, which `hiscock_check` holds to the published chart of `schwarzschild.json` at $r_s = 2m$ slot by slot, and the evaporating hole where $m$ falls.
Its one Einstein component is $G_{vv} = 2\,\partial_v m/r^2$, null dust of negative energy density moving inward while the mass falls.
The outgoing chart is $ds^2 = -(1 - 2m/r)\,du^2 - 2\,du\,dr + r^2d\Omega^2$ with $m = m(u)$, and $G_{uu} = -2\,\partial_u m/r^2$, null dust of positive energy density moving outward.
The flat chart is $ds^2 = -du\,dv + \tfrac{1}{4}(v - u)^2d\Omega^2$, the ingoing chart without mass along $r = (v - u)/2$, which `hiscock_check` checks by pulling it back.
Martin-Dussaud and Rovelli also write Schwarzschild's region in double null coordinates, $r = 2m(1 + W(e^{(v - u)/4m - 1}))$ with Lambert's $W$; the checker's `Reader` has no such function, so that chart is not published, and the region is the ingoing chart with a constant mass.

Each Vaidya chart declares $R$, the areal radius of the surface of pair creation, as a function of its null time; it stands in the domains alone, $r \le R$ in the ingoing chart and $r \ge R$ in the outgoing one.

## Step 3. The surface of pair creation

Let the surface be $r = R(v)$ in the ingoing chart, with the mass $N(v)$ inside and $M(u)$ outside.
Its induced metric is $(-F + 2R')\,dv^2 + R^2d\Omega^2$ from inside, $F = 1 - 2N/R$, and $(-F_+ - 2\,dR/du)\,du^2 + R^2d\Omega^2$ from outside.
Israel's first condition makes them one metric, and his second, with the normal toward larger $r$, gives the surface energy density and pressure

$$\sigma = -\frac{[[K^\theta{}_\theta]]}{4\pi}, \qquad p = \frac{[[K^\tau{}_\tau + K^\theta{}_\theta]]}{8\pi}.$$

At a fixed radius $r_0$, $K^\theta{}_\theta = \sqrt{F}/r_0$ on each side, so $\sigma = (\sqrt{F_-} - \sqrt{F_+})/4\pi r_0$, which vanishes exactly when the two masses agree on the surface.
With $M = N$ there, the first condition is $du/dv = 1 - 2R'/F$, and what is left is

$$p = -\frac{N'}{4\pi r_0F^{3/2}} \quad\text{at a fixed radius } r_0,$$

a pressure while the mass falls: no surface energy density and a negative tension, as Hayward found and as James Bardeen had pointed out to him.
The two streams carry no net energy away from the surface and a net momentum outward, and the layer is what balances it.
`_tools/test_hiscock.py` holds all of this to the published metrics and Christoffel symbols of the two charts, for a surface at a fixed radius and for the one drawn.

## Step 4. The model drawn

Every diagram draws one model, `HISCOCK_MASS` in `null_rays.py`, in units of the mass $m_0$ of a shell of light that falls in along $v = 0$ and makes the hole.
From $v_1 = 2$ to $v_0 = 8$ the mass inside is $N(v) = \cos^2(\pi(v - 2)/12)$, which starts and ends with no slope, so $dN/dv \to 0$ as $N \to 0$, the condition Hiscock found necessary for a finite flux on the Cauchy horizon.
The surface of pair creation is $R = 3N(v)$, outside the apparent horizon $r = 2N$ by half its radius and timelike while the mass falls, since $ds^2 = (-\tfrac13 + 6N')\,dv^2$ on it.
That choice is ours: Hiscock's abstract fixes only that the surface is timelike, near the apparent horizon and outside it.
On it $F = 1/3$, so $du/dv = 1 - 18N'$ and

$$u = v - 18N(v),$$

normalised to $u = v - 2r$ in the flat space left at the end, which puts the first of the radiation on $u_1 = -16$ and the last on $u_0 = v_0 = 8$.
The mass outside is $M(u) = N(v(u))$, with $dM/du = N'/(1 - 18N') \to 0$ as $M \to 0$, Hiscock's condition in retarded time.
Before $u_1$ the outgoing chart is Schwarzschild's with $u = v - 2r_* - 12 - 4\ln 2$, $r_* = r + 2\ln(r/2 - 1)$, and begins on the shell.

The event horizon is the outgoing ray that reaches $r = 0$ at $v_0$.
Outgoing rays converge on it toward the past, so it is found by integrating $dr/dv = (1 - 2N/r)/2$ back from $r = 0.0007$ at $v = 7.9$, inside the apparent horizon $r = 0.0014$: it crosses the shell at $r = 1.706$ and left the centre at $v = -3.412$.
It lies inside the apparent horizon from the shell to the end, so light that starts between the two escapes.

## Step 5. The diagrams

The spacetime diagrams draw the ingoing chart inside the surface, the outgoing chart beyond it, the flat space after the last ray, and the simplest model, two shells with $v_1 = 6$, whose event horizon and Cauchy horizon are the rays into and out of the last point of the singularity.

The conformal diagram draws the whole spacetime by tracing both rays of every event, `hiscock` in `conformal.py`.
Its $q$ is a function of the advanced time each ingoing ray had at past null infinity, which is the ingoing chart's $v$ only until $v_1$: a ray that crosses the radiation is delayed, and the last ray into the hole, $v_0 = 8$ inside, left past null infinity at $23.79$.
Its $p$ is a function of the advanced time $w$ at which each outgoing ray left the centre of the flat space inside the shell, the same function, which puts that centre on $X = 0$; the rays of the flat space after the hole take $q$'s function of their own starting time, moved down to be continuous across $u_0$.
The singularity runs from the event where the shell reaches the centre to the last point, where the event horizon, the apparent horizon, the surface of pair creation and the Cauchy horizon all meet.

The embedding diagram is a movie of the moments $v - r = T$ inside the surface and $u + r = T'$ outside it, with $T' = v - 15N(v)$ at the $v$ where the two meet, $v - 3N(v) = T$.
Both carry $(1 + 2m/r)\,dr^2 + r^2d\phi^2$, and the mass is the same on the two sides of the surface, so the profile $dz/dr = \sqrt{2m/r}$ has no fold there.
