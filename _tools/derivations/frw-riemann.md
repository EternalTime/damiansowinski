# The Riemann tensor of the Friedmann-Lemaître-Robertson-Walker universe

This is the working behind the Riemann block of `MFS/assets/data/metrics/frw.json`, and behind the two components it forces elsewhere in that file.
It is narrower than the other files in this directory, which derive a whole entry.
The FLRW Christoffel symbols, Ricci tensor, Einstein tensor, Ricci scalar and Kretschmann scalar were already right when this was written, and the Weyl block was already empty, which is what a conformally flat spacetime requires.
What was wrong was the curvature tensor sitting between the connection and the Ricci tensor, so that the entry published a Riemann tensor which does not contract to its own Ricci tensor.
Every component this file changed is computed here from the connection, by hand, so that the captain can check any one of them without running anything.

The companion script `_tools/derivations/verify_metrics.py` does the same work in sympy and compares it against the published file.
Running it on `frw/comoving_spherical` and `frw/conformal_spherical` takes about thirteen seconds, which makes this the cheapest entry in the collection to re-check.

---

## Step 1. Conventions

The signature is $(-,+,+,+)$.
The connection is the Levi-Civita one,

$$\Gamma^\mu{}_{\nu\rho} = \tfrac{1}{2}g^{\mu\alpha}\left(\partial_\nu g_{\alpha\rho} + \partial_\rho g_{\alpha\nu} - \partial_\alpha g_{\nu\rho}\right),$$

and the Riemann tensor is

$$R^\mu{}_{\nu\rho\sigma} = \partial_\rho \Gamma^\mu{}_{\nu\sigma} - \partial_\sigma \Gamma^\mu{}_{\nu\rho} + \Gamma^\mu{}_{\rho\lambda}\Gamma^\lambda{}_{\nu\sigma} - \Gamma^\mu{}_{\sigma\lambda}\Gamma^\lambda{}_{\nu\rho}.$$

The Ricci tensor is the standard contraction on the first lower index,

$$R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu},$$

which the collection settled on for every entry on 18 September 2026.
On that contraction FLRW reads $R_{tt} = -3\ddot{a}/a$ and $G_{tt} = 3(\dot{a}^2+k)/a^2$, the second of which is positive for ordinary matter, and both of which the entry already published.
Step 7 shows that the corrected Riemann tensor of this file contracts to exactly those, which is the whole point of the correction: the published Ricci tensor was never in doubt, and the Riemann tensor had to be brought into line with it rather than the other way round.

Factors of $c$ and $G$ are explicit.
The chart is the one whose time coordinate is $x^0 = ct$, which is why the entry prints $g_{tt} = -1$ against a line element whose time term is $-c^2dt^2$.
Because the rescaling $t \to ct$ is linear with constant coefficients, a component in this chart is the component taken with the bare coordinate multiplied by $c$ once for every upper time index and divided by $c$ once for every lower one, and the Christoffel symbols follow the same rule because the inhomogeneous term in their transformation law carries a second derivative of the coordinate change, which vanishes for a linear one.

Every dot in this file and in the entry is a derivative with respect to that chart coordinate,

$$\dot{a} \equiv \frac{da}{dx^0} = \frac{1}{c}\frac{da}{dt}, \qquad \ddot{a} \equiv \frac{d^2a}{(dx^0)^2} = \frac{1}{c^2}\frac{d^2a}{dt^2}.$$

This is what keeps the printed components free of $c$ while every one of them is a component in the $x^0 = ct$ chart.
It is also what makes them dimensionally honest.
The scale factor $a$ is dimensionless and $r$ is a comoving length, so $k$ is a curvature carrying $1/L^2$, and $\dot{a}$, being a derivative with respect to a length, carries $1/L$.
Then $\dot{a}^2 + k$ carries $1/L^2$ and $\ddot{a}/a$ carries $1/L^2$, which is what a Riemann component with one upper and three lower indices, all of them lengths, has to carry.
That last sentence is the whole of the dimensional argument that condemned four of the old components on sight, and Step 8 returns to it.

Working in the $x^0$ chart from the start, rather than in the bare $t$ chart and converting at the end, means no factor of $c$ appears anywhere below.
Every expression in Steps 2 to 9 is already in the form the entry prints.

---

## Step 2. The line element and the metric

$$ds^2 = -c^2dt^2 + a^2\left(\frac{dr^2}{1-kr^2} + r^2d\theta^2 + r^2\sin^2\theta\,d\phi^2\right),$$

which in the chart coordinate $x^0 = ct$ is

$$ds^2 = -(dx^0)^2 + a^2\left(\frac{dr^2}{1-kr^2} + r^2d\theta^2 + r^2\sin^2\theta\,d\phi^2\right).$$

The metric is diagonal, and in the order $(t,r,\theta,\phi)$ its components are

$$g_{tt} = -1, \qquad g_{rr} = \frac{a^2}{1-kr^2}, \qquad g_{\theta\theta} = a^2r^2, \qquad g_{\phi\phi} = a^2r^2\sin^2\theta,$$

with the inverse

$$g^{tt} = -1, \qquad g^{rr} = \frac{1-kr^2}{a^2}, \qquad g^{\theta\theta} = \frac{1}{a^2r^2}, \qquad g^{\phi\phi} = \frac{1}{a^2r^2\sin^2\theta}.$$

Two abbreviations are used throughout, purely to keep the lines short:

$$f \equiv 1 - kr^2, \qquad f' = \frac{df}{dr} = -2kr.$$

The scale factor depends on $x^0$ alone and on nothing else, and $f$ depends on $r$ alone and on nothing else.
That single fact kills most of the terms below before they are written.

---

## Step 3. The connection

These are the entry's own published Christoffel symbols, restated because every Riemann component is built from them.
The checker agrees with all nineteen, and none of them was touched by this correction.
They are listed once, with the symmetric partner $\Gamma^\mu{}_{\rho\nu} = \Gamma^\mu{}_{\nu\rho}$ left implicit.

$$\Gamma^t{}_{rr} = \frac{a\dot{a}}{f}, \qquad \Gamma^t{}_{\theta\theta} = r^2a\dot{a}, \qquad \Gamma^t{}_{\phi\phi} = r^2a\dot{a}\sin^2\theta,$$

$$\Gamma^r{}_{tr} = \Gamma^\theta{}_{t\theta} = \Gamma^\phi{}_{t\phi} = \frac{\dot{a}}{a},$$

$$\Gamma^r{}_{rr} = \frac{kr}{f}, \qquad \Gamma^r{}_{\theta\theta} = -rf, \qquad \Gamma^r{}_{\phi\phi} = -rf\sin^2\theta,$$

$$\Gamma^\theta{}_{r\theta} = \Gamma^\phi{}_{r\phi} = \frac{1}{r}, \qquad \Gamma^\theta{}_{\phi\phi} = -\sin\theta\cos\theta, \qquad \Gamma^\phi{}_{\theta\phi} = \cot\theta.$$

The entry writes $\Gamma^r{}_{\theta\theta} = r(kr^2-1)$ and $\Gamma^r{}_{\phi\phi} = r(kr^2-1)\sin^2\theta$, which are $-rf$ and $-rf\sin^2\theta$.

Three features of this list do all the work below.
There is no $\Gamma^\mu{}_{tt}$ of any kind, because the comoving observers are geodesic, so every term carrying $\Gamma^\lambda{}_{tt}$ drops.
There is no $\Gamma^t{}_{t\mu}$ of any kind, because $g_{tt}$ is the constant $-1$, so every term carrying $\Gamma^t{}_{t\lambda}$ drops.
And the three symbols $\Gamma^r{}_{tr}$, $\Gamma^\theta{}_{t\theta}$ and $\Gamma^\phi{}_{t\phi}$ are the same number $\dot{a}/a$, which is why the three components $R^r{}_{ttr}$, $R^\theta{}_{tt\theta}$ and $R^\phi{}_{tt\phi}$ come out equal in Step 4.

---

## Step 4. The six independent Riemann components, computed

A Riemann tensor in four dimensions has twenty independent components in general.
FLRW is isotropic and homogeneous, so almost all of them vanish and the survivors fall into two families: the ones with a time index pair, which carry $\ddot{a}$, and the purely spatial ones, which carry $\dot{a}^2 + k$.
Six computations fix the whole tensor, and the other eighteen nonzero components follow from them by the antisymmetry $R^\mu{}_{\nu\rho\sigma} = -R^\mu{}_{\nu\sigma\rho}$ and by relabelling.

### 4a. $R^r{}_{ttr}$, and with it $R^\theta{}_{tt\theta}$ and $R^\phi{}_{tt\phi}$

Take $\mu = r$, $\nu = t$, $\rho = t$, $\sigma = r$ in the definition:

$$R^r{}_{ttr} = \partial_t\Gamma^r{}_{tr} - \partial_r\Gamma^r{}_{tt} + \Gamma^r{}_{t\lambda}\Gamma^\lambda{}_{tr} - \Gamma^r{}_{r\lambda}\Gamma^\lambda{}_{tt}.$$

The second and fourth terms vanish because no $\Gamma^\mu{}_{tt}$ exists.
The first is

$$\partial_t\frac{\dot{a}}{a} = \frac{\ddot{a}}{a} - \frac{\dot{a}^2}{a^2},$$

and the third has the single surviving value $\lambda = r$,

$$\Gamma^r{}_{tr}\Gamma^r{}_{tr} = \frac{\dot{a}^2}{a^2}.$$

The two $\dot{a}^2/a^2$ cancel and

$$\boxed{R^r{}_{ttr} = \frac{\ddot{a}}{a}.}$$

The same three lines with $r$ replaced by $\theta$ or by $\phi$ give the same answer, because $\Gamma^\theta{}_{t\theta}$ and $\Gamma^\phi{}_{t\phi}$ are the same $\dot{a}/a$ and the sum again collapses to its diagonal term:

$$R^\theta{}_{tt\theta} = R^\phi{}_{tt\phi} = \frac{\ddot{a}}{a}.$$

These three are the components the old file got most conspicuously wrong, and the value $\ddot{a}/a$ is the same for all three no matter what $r$ or $\theta$ is.
That is isotropy, and it is worth noticing that it is impossible to reconcile with a published value carrying an $r^2$ or a $\sin^2\theta$.

### 4b. $R^t{}_{rtr}$

Take $\mu = t$, $\nu = r$, $\rho = t$, $\sigma = r$:

$$R^t{}_{rtr} = \partial_t\Gamma^t{}_{rr} - \partial_r\Gamma^t{}_{rt} + \Gamma^t{}_{t\lambda}\Gamma^\lambda{}_{rr} - \Gamma^t{}_{r\lambda}\Gamma^\lambda{}_{rt}.$$

The second and third terms vanish because no $\Gamma^t{}_{t\mu}$ exists.
The first is

$$\partial_t\frac{a\dot{a}}{f} = \frac{\dot{a}^2 + a\ddot{a}}{f},$$

since $f$ does not depend on time.
The fourth has the single surviving value $\lambda = r$,

$$\Gamma^t{}_{rr}\Gamma^r{}_{rt} = \frac{a\dot{a}}{f}\cdot\frac{\dot{a}}{a} = \frac{\dot{a}^2}{f}.$$

Subtracting,

$$\boxed{R^t{}_{rtr} = \frac{a\ddot{a}}{f} = \frac{a\ddot{a}}{1-kr^2}.}$$

### 4c. $R^t{}_{\theta t\theta}$ and $R^t{}_{\phi t\phi}$

The same pattern with $r$ replaced by $\theta$:

$$R^t{}_{\theta t\theta} = \partial_t\Gamma^t{}_{\theta\theta} - \Gamma^t{}_{\theta\theta}\Gamma^\theta{}_{\theta t} = \partial_t\left(r^2a\dot{a}\right) - r^2a\dot{a}\cdot\frac{\dot{a}}{a} = r^2\left(\dot{a}^2 + a\ddot{a}\right) - r^2\dot{a}^2,$$

so

$$\boxed{R^t{}_{\theta t\theta} = r^2a\ddot{a}, \qquad R^t{}_{\phi t\phi} = r^2a\ddot{a}\sin^2\theta,}$$

the second by the same computation with the extra constant factor $\sin^2\theta$ riding along in $\Gamma^t{}_{\phi\phi}$.

### 4d. $R^r{}_{\theta r\theta}$, and with it $R^r{}_{\phi r\phi}$

This is the first purely spatial one, and the first that needs more than two terms.
Take $\mu = r$, $\nu = \theta$, $\rho = r$, $\sigma = \theta$:

$$R^r{}_{\theta r\theta} = \partial_r\Gamma^r{}_{\theta\theta} - \partial_\theta\Gamma^r{}_{\theta r} + \Gamma^r{}_{r\lambda}\Gamma^\lambda{}_{\theta\theta} - \Gamma^r{}_{\theta\lambda}\Gamma^\lambda{}_{\theta r}.$$

The second term vanishes because $\Gamma^r{}_{\theta r} = 0$.
The first is

$$\partial_r(-rf) = -f - rf' = -(1-kr^2) + 2kr^2 = -1 + 3kr^2.$$

The third runs over $\lambda = t$ and $\lambda = r$:

$$\Gamma^r{}_{rt}\Gamma^t{}_{\theta\theta} = \frac{\dot{a}}{a}\cdot r^2a\dot{a} = r^2\dot{a}^2, \qquad \Gamma^r{}_{rr}\Gamma^r{}_{\theta\theta} = \frac{kr}{f}\cdot(-rf) = -kr^2.$$

The fourth has the single surviving value $\lambda = \theta$:

$$\Gamma^r{}_{\theta\theta}\Gamma^\theta{}_{\theta r} = (-rf)\cdot\frac{1}{r} = -(1-kr^2).$$

Adding, with the fourth term subtracted,

$$R^r{}_{\theta r\theta} = \left(-1 + 3kr^2\right) + r^2\dot{a}^2 - kr^2 + \left(1 - kr^2\right) = r^2\dot{a}^2 + kr^2,$$

so

$$\boxed{R^r{}_{\theta r\theta} = r^2\left(\dot{a}^2+k\right), \qquad R^r{}_{\phi r\phi} = r^2\left(\dot{a}^2+k\right)\sin^2\theta.}$$

The $\phi$ version is the same four lines with $\sin^2\theta$ multiplying $\Gamma^t{}_{\phi\phi}$, $\Gamma^r{}_{\phi\phi}$ and nothing else, so it factors straight out.

### 4e. $R^\theta{}_{r\theta r}$, and with it $R^\phi{}_{r\phi r}$

Take $\mu = \theta$, $\nu = r$, $\rho = \theta$, $\sigma = r$:

$$R^\theta{}_{r\theta r} = \partial_\theta\Gamma^\theta{}_{rr} - \partial_r\Gamma^\theta{}_{r\theta} + \Gamma^\theta{}_{\theta\lambda}\Gamma^\lambda{}_{rr} - \Gamma^\theta{}_{r\lambda}\Gamma^\lambda{}_{r\theta}.$$

The first term vanishes because $\Gamma^\theta{}_{rr} = 0$.
The second is

$$-\partial_r\frac{1}{r} = \frac{1}{r^2}.$$

The third runs over $\lambda = t$ and $\lambda = r$:

$$\Gamma^\theta{}_{\theta t}\Gamma^t{}_{rr} = \frac{\dot{a}}{a}\cdot\frac{a\dot{a}}{f} = \frac{\dot{a}^2}{f}, \qquad \Gamma^\theta{}_{\theta r}\Gamma^r{}_{rr} = \frac{1}{r}\cdot\frac{kr}{f} = \frac{k}{f}.$$

The fourth has the single surviving value $\lambda = \theta$:

$$\Gamma^\theta{}_{r\theta}\Gamma^\theta{}_{r\theta} = \frac{1}{r^2}.$$

The two $1/r^2$ cancel and

$$\boxed{R^\theta{}_{r\theta r} = \frac{\dot{a}^2+k}{f} = \frac{\dot{a}^2+k}{1-kr^2}, \qquad R^\phi{}_{r\phi r} = \frac{\dot{a}^2+k}{1-kr^2}.}$$

The $\phi$ version is the same computation with $\Gamma^\phi{}_{\phi t} = \dot{a}/a$ and $\Gamma^\phi{}_{\phi r} = 1/r$ in place of the $\theta$ ones, and every one of those is the same number as its $\theta$ counterpart.
No $\sin^2\theta$ enters, because none of the four symbols used carries one.
This is the component the old file published with a spurious $\sin^2\theta$, and Step 6 explains where the reader's intuition for that factor goes wrong.

### 4f. $R^\theta{}_{\phi\theta\phi}$ and $R^\phi{}_{\theta\phi\theta}$, which are not equal

These two are the pair worth doing slowly, because the correction turns on them in both charts.

Take $\mu = \theta$, $\nu = \phi$, $\rho = \theta$, $\sigma = \phi$:

$$R^\theta{}_{\phi\theta\phi} = \partial_\theta\Gamma^\theta{}_{\phi\phi} - \partial_\phi\Gamma^\theta{}_{\phi\theta} + \Gamma^\theta{}_{\theta\lambda}\Gamma^\lambda{}_{\phi\phi} - \Gamma^\theta{}_{\phi\lambda}\Gamma^\lambda{}_{\phi\theta}.$$

The second term vanishes because $\Gamma^\theta{}_{\phi\theta} = 0$.
The first is

$$\partial_\theta(-\sin\theta\cos\theta) = \sin^2\theta - \cos^2\theta.$$

The third runs over $\lambda = t$ and $\lambda = r$:

$$\Gamma^\theta{}_{\theta t}\Gamma^t{}_{\phi\phi} = \frac{\dot{a}}{a}\cdot r^2a\dot{a}\sin^2\theta = r^2\dot{a}^2\sin^2\theta, \qquad \Gamma^\theta{}_{\theta r}\Gamma^r{}_{\phi\phi} = \frac{1}{r}\cdot(-rf\sin^2\theta) = -(1-kr^2)\sin^2\theta.$$

The fourth has the single surviving value $\lambda = \phi$:

$$\Gamma^\theta{}_{\phi\phi}\Gamma^\phi{}_{\phi\theta} = (-\sin\theta\cos\theta)(\cot\theta) = -\cos^2\theta.$$

Adding, with the fourth subtracted, the two $\cos^2\theta$ cancel and the bare $\sin^2\theta$ cancels against the $-\sin^2\theta$ inside $-(1-kr^2)\sin^2\theta$:

$$R^\theta{}_{\phi\theta\phi} = \sin^2\theta + r^2\dot{a}^2\sin^2\theta - \sin^2\theta + kr^2\sin^2\theta,$$

so

$$\boxed{R^\theta{}_{\phi\theta\phi} = r^2\left(\dot{a}^2+k\right)\sin^2\theta.}$$

Now take $\mu = \phi$, $\nu = \theta$, $\rho = \phi$, $\sigma = \theta$, which is the same four indices with the raised one swapped:

$$R^\phi{}_{\theta\phi\theta} = \partial_\phi\Gamma^\phi{}_{\theta\theta} - \partial_\theta\Gamma^\phi{}_{\theta\phi} + \Gamma^\phi{}_{\phi\lambda}\Gamma^\lambda{}_{\theta\theta} - \Gamma^\phi{}_{\theta\lambda}\Gamma^\lambda{}_{\theta\phi}.$$

The first term vanishes because $\Gamma^\phi{}_{\theta\theta} = 0$.
The second is

$$-\partial_\theta\cot\theta = \csc^2\theta.$$

The third runs over $\lambda = t$ and $\lambda = r$, and the value $\lambda = \theta$ contributes nothing because $\Gamma^\theta{}_{\theta\theta} = 0$:

$$\Gamma^\phi{}_{\phi t}\Gamma^t{}_{\theta\theta} = \frac{\dot{a}}{a}\cdot r^2a\dot{a} = r^2\dot{a}^2, \qquad \Gamma^\phi{}_{\phi r}\Gamma^r{}_{\theta\theta} = \frac{1}{r}\cdot(-rf) = -(1-kr^2).$$

The fourth has the single surviving value $\lambda = \phi$:

$$\Gamma^\phi{}_{\theta\phi}\Gamma^\phi{}_{\theta\phi} = \cot^2\theta.$$

Adding, with the fourth subtracted, and using $\csc^2\theta - \cot^2\theta = 1$:

$$R^\phi{}_{\theta\phi\theta} = \csc^2\theta + r^2\dot{a}^2 - 1 + kr^2 - \cot^2\theta = 1 + r^2\dot{a}^2 - 1 + kr^2,$$

so

$$\boxed{R^\phi{}_{\theta\phi\theta} = r^2\left(\dot{a}^2+k\right).}$$

The two differ by exactly $\sin^2\theta$, and neither is a typo for the other.
Step 6 says why in one line.

---

## Step 5. The mixed block in full

Antisymmetry in the last two indices, $R^\mu{}_{\nu\rho\sigma} = -R^\mu{}_{\nu\sigma\rho}$, doubles the six boxed values of Step 4 into twelve, and the relabelling already noted in Step 4 carries them to all twenty-four nonzero components.
This is the corrected `frw/comoving_spherical.riemann.ulll` block, in the order the file publishes it.

| $R^\mu{}_{\nu\rho\sigma}$ | value | | $R^\mu{}_{\nu\rho\sigma}$ | value |
| --- | --- | --- | --- | --- |
| $R^t{}_{rtr}$ | $\dfrac{a\ddot{a}}{1-kr^2}$ | | $R^\theta{}_{tt\theta}$ | $\dfrac{\ddot{a}}{a}$ |
| $R^t{}_{rrt}$ | $-\dfrac{a\ddot{a}}{1-kr^2}$ | | $R^\theta{}_{t\theta t}$ | $-\dfrac{\ddot{a}}{a}$ |
| $R^t{}_{\theta t\theta}$ | $r^2a\ddot{a}$ | | $R^\theta{}_{rr\theta}$ | $-\dfrac{\dot{a}^2+k}{1-kr^2}$ |
| $R^t{}_{\theta\theta t}$ | $-r^2a\ddot{a}$ | | $R^\theta{}_{r\theta r}$ | $\dfrac{\dot{a}^2+k}{1-kr^2}$ |
| $R^t{}_{\phi t\phi}$ | $r^2a\ddot{a}\sin^2\theta$ | | $R^\theta{}_{\phi\theta\phi}$ | $r^2(\dot{a}^2+k)\sin^2\theta$ |
| $R^t{}_{\phi\phi t}$ | $-r^2a\ddot{a}\sin^2\theta$ | | $R^\theta{}_{\phi\phi\theta}$ | $-r^2(\dot{a}^2+k)\sin^2\theta$ |
| $R^r{}_{ttr}$ | $\dfrac{\ddot{a}}{a}$ | | $R^\phi{}_{tt\phi}$ | $\dfrac{\ddot{a}}{a}$ |
| $R^r{}_{trt}$ | $-\dfrac{\ddot{a}}{a}$ | | $R^\phi{}_{t\phi t}$ | $-\dfrac{\ddot{a}}{a}$ |
| $R^r{}_{\theta r\theta}$ | $r^2(\dot{a}^2+k)$ | | $R^\phi{}_{rr\phi}$ | $-\dfrac{\dot{a}^2+k}{1-kr^2}$ |
| $R^r{}_{\theta\theta r}$ | $-r^2(\dot{a}^2+k)$ | | $R^\phi{}_{r\phi r}$ | $\dfrac{\dot{a}^2+k}{1-kr^2}$ |
| $R^r{}_{\phi r\phi}$ | $r^2(\dot{a}^2+k)\sin^2\theta$ | | $R^\phi{}_{\theta\theta\phi}$ | $-r^2(\dot{a}^2+k)$ |
| $R^r{}_{\phi\phi r}$ | $-r^2(\dot{a}^2+k)\sin^2\theta$ | | $R^\phi{}_{\theta\phi\theta}$ | $r^2(\dot{a}^2+k)$ |

Eight of these twenty-four were already right in the file, and they are the eight in the $R^r{}_{\theta\cdot\cdot}$, $R^r{}_{\phi\cdot\cdot}$, $R^\theta{}_{r\cdot\cdot}$ and $R^\theta{}_{\phi\cdot\cdot}$ rows.
The other sixteen are Step 9.

---

## Step 6. Why $\sin^2\theta$ rides with the lower index and not the upper one

The rule the old file appears to have followed is that a $\phi$ index anywhere brings a $\sin^2\theta$.
It does not, and the reason is worth stating once because it decides four components in the comoving chart and two in the conformal one.

The metric carries the factor: $g_{\phi\phi} = g_{\theta\theta}\sin^2\theta$.
So lowering a $\phi$ index multiplies by $\sin^2\theta$ relative to lowering a $\theta$ index, and raising one divides by it.
Take the pair of Step 4f and lower the first index of each, using $g_{\theta\theta} = a^2r^2$ and $g_{\phi\phi} = a^2r^2\sin^2\theta$:

$$R_{\theta\phi\theta\phi} = g_{\theta\theta}R^\theta{}_{\phi\theta\phi} = a^2r^2\cdot r^2(\dot{a}^2+k)\sin^2\theta = r^4a^2(\dot{a}^2+k)\sin^2\theta,$$

$$R_{\phi\theta\phi\theta} = g_{\phi\phi}R^\phi{}_{\theta\phi\theta} = a^2r^2\sin^2\theta\cdot r^2(\dot{a}^2+k) = r^4a^2(\dot{a}^2+k)\sin^2\theta.$$

They agree, as the pair exchange symmetry $R_{\mu\nu\rho\sigma} = R_{\rho\sigma\mu\nu}$ demands.
And they can only agree if the mixed components differ by exactly the $\sin^2\theta$ that separates $g_{\theta\theta}$ from $g_{\phi\phi}$.
So the two mixed components are obliged to be unequal, and the file publishing them equal was the error, not the inequality.

The same argument settles $R^\phi{}_{rr\phi}$ against $R^r{}_{\phi\phi r}$, which the pair exchange also forces to share a lowered value:

$$R_{\phi rr\phi} = g_{\phi\phi}R^\phi{}_{rr\phi} = a^2r^2\sin^2\theta\cdot\left(-\frac{\dot{a}^2+k}{1-kr^2}\right) = -\frac{r^2a^2(\dot{a}^2+k)\sin^2\theta}{1-kr^2},$$

$$R_{r\phi\phi r} = g_{rr}R^r{}_{\phi\phi r} = \frac{a^2}{1-kr^2}\cdot\left(-r^2(\dot{a}^2+k)\sin^2\theta\right) = -\frac{r^2a^2(\dot{a}^2+k)\sin^2\theta}{1-kr^2}.$$

The one $\sin^2\theta$ they both carry arrives from $g_{\phi\phi}$ in the first and from the mixed component itself in the second.
It is counted once, not twice, which is exactly what the old $R^\phi{}_{rr\phi}$ got wrong by carrying its own $\sin^2\theta$ on top of the one $g_{\phi\phi}$ would supply.

The short version, good for checking any component by eye: each $\phi$ index counts $+\tfrac{1}{2}$ when it sits below and $-\tfrac{1}{2}$ when it sits above, and the total is the power of $\sin^2\theta$ the component carries.
For $R^\theta{}_{\phi\theta\phi}$ that is $+\tfrac{1}{2}+\tfrac{1}{2} = 1$, for $R^\phi{}_{\theta\phi\theta}$ it is $-\tfrac{1}{2}+\tfrac{1}{2} = 0$, and the same count reads $g_{\phi\phi}$ and $g^{\phi\phi}$ correctly too.

---

## Step 7. The fully lowered block, and the check that closes it

Lowering the first index of a diagonal metric is a single multiplication with no sum,

$$R_{\mu\nu\rho\sigma} = g_{\mu\mu}R^\mu{}_{\nu\rho\sigma} \quad (\text{no sum over } \mu),$$

which gives the corrected `frw/comoving_spherical.riemann.llll` block.
The file previously published twelve of these twenty-four and left twelve out, so this is the block that grew.

The six independent values are

$$R_{trtr} = g_{tt}R^t{}_{rtr} = -\frac{a\ddot{a}}{1-kr^2}, \qquad R_{rttr} = g_{rr}R^r{}_{ttr} = \frac{a^2}{1-kr^2}\cdot\frac{\ddot{a}}{a} = \frac{a\ddot{a}}{1-kr^2},$$

$$R_{t\theta t\theta} = -r^2a\ddot{a}, \qquad R_{\theta tt\theta} = a^2r^2\cdot\frac{\ddot{a}}{a} = r^2a\ddot{a},$$

$$R_{t\phi t\phi} = -r^2a\ddot{a}\sin^2\theta, \qquad R_{\phi tt\phi} = a^2r^2\sin^2\theta\cdot\frac{\ddot{a}}{a} = r^2a\ddot{a}\sin^2\theta,$$

$$R_{r\theta r\theta} = \frac{a^2}{1-kr^2}\cdot r^2(\dot{a}^2+k) = \frac{r^2a^2(\dot{a}^2+k)}{1-kr^2}, \qquad R_{\theta r\theta r} = a^2r^2\cdot\frac{\dot{a}^2+k}{1-kr^2} = \frac{r^2a^2(\dot{a}^2+k)}{1-kr^2},$$

$$R_{r\phi r\phi} = \frac{r^2a^2(\dot{a}^2+k)}{1-kr^2}\sin^2\theta, \qquad R_{\phi r\phi r} = a^2r^2\sin^2\theta\cdot\frac{\dot{a}^2+k}{1-kr^2} = \frac{r^2a^2(\dot{a}^2+k)}{1-kr^2}\sin^2\theta,$$

$$R_{\theta\phi\theta\phi} = R_{\phi\theta\phi\theta} = r^4a^2(\dot{a}^2+k)\sin^2\theta,$$

and the remaining twelve are these with the last two indices swapped and the sign reversed.

The three pairs that had to come out equal did: $R_{trtr}$ against $R_{rtrt} = -R_{rttr}$, $R_{r\theta r\theta}$ against $R_{\theta r\theta r}$, and $R_{\theta\phi\theta\phi}$ against $R_{\phi\theta\phi\theta}$.
None of those equalities was imposed.
Each is a consequence of six separately computed mixed components meeting four separately written metric factors, so the pair exchange symmetry is a real check on the whole of Step 4 and not a restatement of it.

### The contraction to the published Ricci tensor

This is the check the old block failed, and the reason the correction was necessary at all.
Contract the corrected mixed block on the first lower index, $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$:

$$R_{tt} = R^r{}_{trt} + R^\theta{}_{t\theta t} + R^\phi{}_{t\phi t} = -\frac{\ddot{a}}{a} - \frac{\ddot{a}}{a} - \frac{\ddot{a}}{a} = -3\frac{\ddot{a}}{a},$$

$$R_{rr} = R^t{}_{rtr} + R^\theta{}_{r\theta r} + R^\phi{}_{r\phi r} = \frac{a\ddot{a}}{1-kr^2} + \frac{\dot{a}^2+k}{1-kr^2} + \frac{\dot{a}^2+k}{1-kr^2} = \frac{a\ddot{a} + 2(\dot{a}^2+k)}{1-kr^2},$$

$$R_{\theta\theta} = R^t{}_{\theta t\theta} + R^r{}_{\theta r\theta} + R^\phi{}_{\theta\phi\theta} = r^2a\ddot{a} + r^2(\dot{a}^2+k) + r^2(\dot{a}^2+k) = r^2\left(a\ddot{a} + 2(\dot{a}^2+k)\right),$$

$$R_{\phi\phi} = R^t{}_{\phi t\phi} + R^r{}_{\phi r\phi} + R^\theta{}_{\phi\theta\phi} = r^2\left(a\ddot{a} + 2(\dot{a}^2+k)\right)\sin^2\theta.$$

All four are exactly what the entry already published, and they were left untouched.
Notice that $R_{\theta\theta}$ needs $R^\phi{}_{\theta\phi\theta}$ without a $\sin^2\theta$ and $R_{\phi\phi}$ needs $R^\theta{}_{\phi\theta\phi}$ with one.
The old file had both carrying $\sin^2\theta$, which cannot produce both of those Ricci components, so the old Riemann block was inconsistent with the same entry's own Ricci tensor whichever way the contraction was taken.

---

## Step 8. The dimensional argument, which condemns four components without any algebra

The dimensional pass in `verify_metrics.py` flagged four terms in the old block, all in the same family, and they are worth seeing because they need no derivative at all.

From Step 1, $a$ is dimensionless, $r$ is a length, $k$ carries $1/L^2$, and $\dot{a}$ and $\ddot{a}$, being derivatives with respect to the length $x^0$, carry $1/L$ and $1/L^2$.

The bookkeeping for an indexed component is the one the checker uses: a curvature carries a base $1/L^2$, and each index then contributes the dimension of its own chart coordinate divided by $L$ when it sits above, and $L$ divided by that dimension when it sits below.
A chart coordinate that is already a length contributes nothing either way, and that covers $t$, whose chart coordinate is $x^0 = ct$, and $r$.
An angle contributes $1/L$ above and $L$ below, so an angle appearing once above and once below contributes nothing either.
That is exactly the shape of $R^\theta{}_{tt\theta}$, whose two $t$ indices contribute nothing and whose two $\theta$ indices cancel each other, so it must carry the bare $1/L^2$.
Now compare:

$$\frac{\ddot{a}}{a} \ \text{carries}\ \frac{1}{L^2}, \qquad r^2a\ddot{a} \ \text{carries}\ L^2\cdot\frac{1}{L^2} = 1.$$

So the old $R^\theta{}_{tt\theta} = r^2a\ddot{a}$ was dimensionless where a curvature was wanted, and no choice of sign or convention could have rescued it.
The same holds for $R^\theta{}_{t\theta t}$, $R^\phi{}_{tt\phi}$ and $R^\phi{}_{t\phi t}$.
The value $r^2a\ddot{a}$ is the correct $R_{\theta tt\theta}$, which carries no net dimension because its one upper index has been traded for a lower one, and that is the fingerprint of the mistake: a fully lowered value printed in a mixed slot.

---

## Step 9. What was wrong, component by component

The sixteen mixed components that changed fall into three groups, and the audit of 18 September 2026 named all three.

**Six were the right value with the wrong sign.**
These are the six $R^t{}_{\cdot\cdot\cdot}$ components.
The file published $R^t{}_{rtr} = -a\ddot{a}/(1-kr^2)$ where Step 4b gives $+a\ddot{a}/(1-kr^2)$, and similarly for the $\theta$ and $\phi$ rows and their antisymmetric partners.
These six are exactly the lowered values in disguise as well, since $g_{tt} = -1$ makes a sign flip and a lowering the same operation on a $t$ index, which is why they read as a sign error rather than as the placement error the next group is.

**Six carried the fully lowered value where the mixed one belongs.**
These are the $R^r{}_{tt r}$, $R^\theta{}_{tt\theta}$ and $R^\phi{}_{tt\phi}$ components and their partners.
The file published $R^r{}_{ttr} = a\ddot{a}/(1-kr^2)$, which is $R_{rttr}$ from Step 7, where Step 4a gives $R^r{}_{ttr} = \ddot{a}/a$.
The two differ by the factor $g_{rr} = a^2/(1-kr^2)$, which is precisely the lowering.
The $\theta$ and $\phi$ ones were published as $r^2a\ddot{a}$ and $r^2a\ddot{a}\sin^2\theta$, which are $R_{\theta tt\theta}$ and $R_{\phi tt\phi}$, and those are the four the dimensional pass caught in Step 8.

**Four carried a $\sin^2\theta$ that belongs to a different slot.**
These are $R^\phi{}_{rr\phi}$, $R^\phi{}_{r\phi r}$, $R^\phi{}_{\theta\theta\phi}$ and $R^\phi{}_{\theta\phi\theta}$, published with $\sin^2\theta$ where Steps 4e and 4f give none.
Step 6 is the whole of the explanation.

The fully lowered block was in a different state: of its twenty-four nonzero components the file published twelve and omitted twelve, and of the twelve published, ten carried the wrong sign.
The two that were right are $R_{\theta\phi\theta\phi}$ and $R_{\theta\phi\phi\theta}$.
The twelve omitted are every component whose first index is $r$ paired with $t$, every component whose first index is $\theta$ paired with $t$ or $r$, and every component whose first index is $\phi$.
Since the block is a list of nonzero components and anything absent from it is read as zero, omitting them asserted that they vanish, which Step 7 shows they do not.

---

## Step 10. The two components this forces elsewhere in the comoving chart

### $R^{rr}$

The entry publishes the Ricci tensor in three variants, and the doubly raised one had $R^{rr}$ raised with $g_{rr}$ where $g^{rr}$ was wanted.
The correct raise is

$$R^{rr} = g^{rr}g^{rr}R_{rr} = \left(\frac{1-kr^2}{a^2}\right)^2\cdot\frac{a\ddot{a}+2(\dot{a}^2+k)}{1-kr^2} = \frac{(1-kr^2)\left(a\ddot{a}+2(\dot{a}^2+k)\right)}{a^4}.$$

The file had $\left(a\ddot{a}+2(\dot{a}^2+k)\right)/\left((1-kr^2)a^4\right)$, which is the above divided by $(1-kr^2)^2$, exactly the error of raising with the metric instead of its inverse twice over.
The other three components of that block were already right, which is what one expects: $g^{tt} = g_{tt} = -1$, so $R^{tt}$ cannot tell the difference, and the $\theta$ and $\phi$ ones were written out correctly.

### Nothing else moves

The Einstein tensor, the Ricci scalar, the Kretschmann scalar and the Weyl block are untouched, and the checker confirms all four.
The Weyl block is empty in both variants, which is correct and was fixed by an earlier job: FLRW is conformally flat, so $C_{\mu\nu\rho\sigma} = 0$ identically, and that is a fact about the spacetime rather than a copy of anything.
It is worth saying plainly that the vanishing was checked rather than assumed: `verify_metrics.py` builds the Weyl tensor from the corrected Riemann tensor, the Ricci tensor and the Ricci scalar through

$$C_{\mu\nu\rho\sigma} = R_{\mu\nu\rho\sigma} - \frac{1}{n-2}\left(g_{\mu\rho}R_{\sigma\nu} - g_{\mu\sigma}R_{\rho\nu} - g_{\nu\rho}R_{\sigma\mu} + g_{\nu\sigma}R_{\rho\mu}\right) + \frac{R}{(n-1)(n-2)}\left(g_{\mu\rho}g_{\sigma\nu} - g_{\mu\sigma}g_{\rho\nu}\right),$$

and finds every one of its 256 components zero, against an empty published block.

---

## Step 11. The conformal chart

The second coordinate system writes the same spacetime as

$$ds^2 = a^2\left(-d\eta^2 + \frac{dr^2}{1-kr^2} + r^2d\theta^2 + r^2\sin^2\theta\,d\phi^2\right),$$

with $\eta$ a length, so that no factor of $c$ appears in the chart at all, and with the prime denoting $d/d\eta$.
Twenty-two of its twenty-four mixed components and all twenty-four of its fully lowered components were already right.
The two that were wrong are the pair of Step 4f, published without their $\sin^2\theta$.

The connection symbols the computation needs are

$$\Gamma^\theta{}_{\theta\eta} = \frac{a'}{a}, \quad \Gamma^\theta{}_{\theta r} = \frac{1}{r}, \quad \Gamma^\theta{}_{\phi\phi} = -\sin\theta\cos\theta, \quad \Gamma^\phi{}_{\theta\phi} = \cot\theta,$$

$$\Gamma^\eta{}_{\phi\phi} = \frac{1}{2}g^{\eta\eta}\left(-\partial_\eta g_{\phi\phi}\right) = \frac{1}{2}\left(-\frac{1}{a^2}\right)\left(-2aa'r^2\sin^2\theta\right) = \frac{a'}{a}r^2\sin^2\theta, \qquad \Gamma^r{}_{\phi\phi} = -rf\sin^2\theta.$$

Then

$$R^\theta{}_{\phi\theta\phi} = \partial_\theta\Gamma^\theta{}_{\phi\phi} + \Gamma^\theta{}_{\theta\eta}\Gamma^\eta{}_{\phi\phi} + \Gamma^\theta{}_{\theta r}\Gamma^r{}_{\phi\phi} - \Gamma^\theta{}_{\phi\phi}\Gamma^\phi{}_{\phi\theta},$$

term by term,

$$\left(\sin^2\theta - \cos^2\theta\right) + \frac{a'^2}{a^2}r^2\sin^2\theta - (1-kr^2)\sin^2\theta + \cos^2\theta = r^2\left(\frac{a'^2}{a^2}+k\right)\sin^2\theta,$$

which is the same cancellation as Step 4f with $\dot{a}$ replaced by $a'/a$, and

$$\boxed{R^\theta{}_{\phi\theta\phi} = r^2\left(\frac{a'^2}{a^2}+k\right)\sin^2\theta, \qquad R^\theta{}_{\phi\phi\theta} = -r^2\left(\frac{a'^2}{a^2}+k\right)\sin^2\theta.}$$

There is a one line check on this inside the entry itself, which needs no derivative at all.
The conformal fully lowered block already published

$$R_{\theta\phi\theta\phi} = r^4a^2\left(\frac{a'^2}{a^2}+k\right)\sin^2\theta,$$

and $g_{\theta\theta} = a^2r^2$, so

$$R^\theta{}_{\phi\theta\phi} = \frac{R_{\theta\phi\theta\phi}}{g_{\theta\theta}} = r^2\left(\frac{a'^2}{a^2}+k\right)\sin^2\theta.$$

The entry was therefore already contradicting itself across its own two variants, and the fully lowered one was the one telling the truth.
Its sibling $R^\phi{}_{\theta\phi\theta} = r^2(a'^2/a^2+k)$ was published without a $\sin^2\theta$ and was already right, which is the same asymmetry Step 6 describes.

---

## Step 12. What the checker says

`verify_metrics.py` compares each published component against sympy in the $x^0 = cT$ chart and separately checks that every term carries the dimension its left hand side fixes.
Before this change the two FLRW systems reported 41 disagreements and 4 terms whose dimensions do not balance.
After it they report none of either:

```
$ python3 _tools/derivations/verify_metrics.py \
      --system frw/comoving_spherical --system frw/conformal_spherical

  frw/comoving_spherical
  frw/conformal_spherical

29 metric files, 2 coordinate systems, 2 checked.

Every published value balances dimensionally and agrees with sympy.
```

The 41 break down as the 16 mixed components of Step 9, the 22 fully lowered ones of the same step, which are 10 wrong values and 12 omissions, the 1 doubly raised Ricci component of Step 10 and the 2 conformal components of Step 11.
That is $16 + 22 + 1 + 2 = 41$.

The full sweep was run twice to show that nothing else moved, once against the tree as it stood and once against the corrected tree, and the two reports were diffed line by line.
The collection went from 113 disagreements to 72.
Every one of the 41 lines that disappeared names `frw`, and no line appeared that was not there before.
The dimensional pass went from 7 terms to 3, and the 4 that went are the four of Step 8; the 3 that remain are the Lanczos-van Stockum doubly raised Einstein tensor, which this change does not touch.
What is left of the 72 is 51 in `interior_schwarzschild`, 13 in `ellis_bronnikov`, 4 in `godel` and 4 in `stockum_dust`, which is the audit's own accounting of what survives once the Weyl group and the FLRW group are gone.

One difference between the two reports is not a disagreement and is worth naming so that the next reader is not puzzled by it.
The `UNCHECKED` count went from 20 to 22, because `kerr/boyer_lindquist.christoffel.ull` and `.lll` finished inside the 120 second budget in the first run and timed out in the second.
The machine was carrying a load average near 150 from other work at the time, and the budget is wall clock.
Neither block reported a disagreement in the run that did finish, so this moved the `UNCHECKED` count and left the disagreement count alone.
`_tools/README.md` already records that Kerr is the entry a full sweep cannot be relied on to complete.
