# The geodesic equations of Godel, Reissner-Nordstrom and the Schwarzschild interior

This is the working behind the `geodesics` field of three entries:

* `godel/cartesian` in `MFS/assets/data/metrics/godel.json`,
* `rn_metric/spherical` in `MFS/assets/data/metrics/rn_metric.json`,
* `interior_schwarzschild/spherical` in `MFS/assets/data/metrics/interior_schwarzschild.json`.

A sweep of every chart in the collection found these three, and only these three, carrying a full set of curvature blocks and no geodesic equations at all.
Every other chart had them.
This file derives the twelve missing equations, four per chart, and nothing else: the curvature blocks of all three entries are left exactly as they were.

The three are unrelated as spacetimes and they are together here because they were missing the same thing.
They are also a useful set to do at once, because each one exercises a different part of the convention.
Godel has an off diagonal metric and a time coordinate carrying no dimensions at all, Reissner-Nordstrom has two horizons and a term whose sign changes inside them, and the Schwarzschild interior has square roots that a computer algebra system will not simplify without being told the signs, which is what decides the form the equations are printed in.

---

## Step 1. Conventions

The signature is $(-,+,+,+)$ and the Christoffel symbols are those of the Levi-Civita connection,

$$\Gamma^\mu{}_{\nu\rho} = \tfrac{1}{2}g^{\mu\alpha}\left(\partial_\nu g_{\alpha\rho} + \partial_\rho g_{\alpha\nu} - \partial_\alpha g_{\nu\rho}\right),$$

which is what the `christoffel` block of each entry publishes in its `ull` variant.

The geodesic equation is the one the whole collection uses,

$$\ddot{x}^\mu + \Gamma^\mu{}_{\nu\rho}\dot{x}^\nu\dot{x}^\rho = 0,$$

with the same $\Gamma$ that is printed above it in the entry, and with the dot the derivative with respect to an affine parameter $\lambda$.

The dots are velocities of the chart, not of the bare coordinate letters, and each entry now says so in its own `convention` field, in the way `kasner.json` and `vaidya.json` say it.
The chart is $x^0 = ct$, so

$$\dot{t} \equiv \frac{d(ct)}{d\lambda}, \qquad \ddot{t} \equiv \frac{d^2(ct)}{d\lambda^2},$$

even though the index and the dot are both printed with the bare letter $t$.
This is the one place in an entry where the chart convention can be lost silently, because a geodesic equation is the only thing an entry publishes that adds a derivative of the time coordinate to a derivative of a space coordinate; Step 3 checks all twelve equations on that reading.

All three entries here use geometric units, so the chart coordinate and the printed letter happen to coincide.
Godel sets $G = c = 1$, and the other two set $c = 1$.
That does not make the convention idle: it fixes what the dimensional pass has to find, and the three entries answer it differently.
Reissner-Nordstrom and the Schwarzschild interior both carry $t$ as a length, beside $r$ and beside the length parameters, so no published component of either carries a factor of $c$.
Godel carries all four coordinates as pure numbers, because $e^x$ needs a dimensionless $x$ and the metric then needs the same of $t$, $y$ and $z$, which leaves the entire length of the solution in $1/\omega$.

The Ricci contraction is the standard one, $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$, settled across the collection on 2026-09-18.
It does not reach anything in this file.
A geodesic is a statement about the connection alone, and the connection is fixed by the metric before any curvature is contracted at all, so the twelve equations below would be the same under either convention.
They are quoted here only so that the file is readable beside the entries, whose curvature blocks are on the standard contraction.

---

## Step 2. What is being built on, and how it was checked

Each entry already publishes its own Christoffel symbols, and the geodesic equations below are assembled from those published symbols by the double sum in Step 1.
That is the fastest route and it is the one taken, so the equations inherit whatever the connection blocks say.
Inheriting is not the same as trusting, and the connection of each entry was therefore recomputed from that entry's own line element, independently of the file, before any equation was written.

The recomputation was done twice.
Once in sympy, from the line element alone, which reproduced all nine nonzero symbols of Godel, all thirteen of Reissner-Nordstrom and all thirteen of the Schwarzschild interior, with nothing missing and nothing extra.
Once by hand, for two symbols of each entry, worked out in full below at Steps A2, B2 and C2.
The hand checks were chosen to be the symbols the equations lean on hardest: the coefficient of $\dot{t}^2$ in the radial equation and the coefficient of $\dot{t}\dot{r}$ or $\dot{t}\dot{x}$ in the time equation, which are the two places a sign or a factor of two would go unnoticed.

All three connections came back clean, in the sense that matters: every published symbol is the right number.
Not one of them had to be corrected, and the equations below are the entries' own connections read back as equations of motion.
What is not clean is the written form of three of the Schwarzschild interior's symbols, which the checker cannot reconcile with its own rebuild even though they agree to forty digits; Step C4 diagnoses that and says what it decided about the twelve equations.

One consequence of building on the published connection is worth stating plainly, because it shows up in Part C.
An expression can be the right number and still be written in a form a computer algebra system will not recognise as equal to another form of the same number.
The Schwarzschild interior is exactly that case, and Step C4 says which form its geodesics are printed in and why.

---

## Part A. The Godel universe

### A1. The line element

The entry publishes, in the chart $(t,x,y,z)$ with all four coordinates running over $\mathbb{R}$,

$$ds^2 = -\frac{1}{2\omega^2}dt^2 - \frac{e^x}{\omega^2}dt\,dy + \frac{1}{2\omega^2}dx^2 - \frac{e^{2x}}{4\omega^2}dy^2 + \frac{1}{2\omega^2}dz^2,$$

which gathers into

$$ds^2 = \frac{1}{2\omega^2}\left[-\left(dt + e^x dy\right)^2 + dx^2 + \tfrac{1}{2}e^{2x}dy^2 + dz^2\right].$$

The second form is the one to read the geometry off.
The metric matrix, in the order $(t,x,y,z)$, is

$$g_{\mu\nu} = \frac{1}{2\omega^2}\begin{pmatrix} -1 & 0 & -e^x & 0 \\ 0 & 1 & 0 & 0 \\ -e^x & 0 & -\tfrac{1}{2}e^{2x} & 0 \\ 0 & 0 & 0 & 1\end{pmatrix}, \qquad \det g = -\frac{e^{2x}}{32\omega^8},$$

and the inverse the entry publishes is

$$g^{\mu\nu} = 2\omega^2\begin{pmatrix} 1 & 0 & -2e^{-x} & 0 \\ 0 & 1 & 0 & 0 \\ -2e^{-x} & 0 & 2e^{-2x} & 0 \\ 0 & 0 & 0 & 1\end{pmatrix}.$$

Two facts about the chart are worth having before the connection.

The determinant never vanishes, so the chart is nowhere degenerate, and it is negative everywhere, so the signature is Lorentzian on the whole of $\mathbb{R}^4$.

And $g_{yy} = -e^{2x}/(4\omega^2)$ is negative at every point, so $\partial_y$ is a timelike vector field everywhere, not merely somewhere.
That is the Godel pathology in its simplest chart form.
The curves it generates are not closed here, because $y$ runs over the whole line, but the sign is the same sign that closes them in the cylindrical chart.

### A2. The connection, with two symbols checked by hand

The entry publishes nine nonzero symbols:

$$\Gamma^t{}_{tx} = \Gamma^t{}_{xt} = 1, \qquad \Gamma^t{}_{xy} = \Gamma^t{}_{yx} = \frac{e^x}{2},$$

$$\Gamma^x{}_{ty} = \Gamma^x{}_{yt} = \frac{e^x}{2}, \qquad \Gamma^x{}_{yy} = \frac{e^{2x}}{2},$$

$$\Gamma^y{}_{tx} = \Gamma^y{}_{xt} = -\cosh x + \sinh x.$$

The last one is printed in hyperbolic functions in the file and is the plain exponential $-e^{-x}$, since $\cosh x - \sinh x = e^{-x}$.
The equations below use the exponential, which is also the form the inverse metric is printed in.

Here is $\Gamma^y{}_{tx}$ by hand.
Writing the definition out,

$$\Gamma^y{}_{tx} = \tfrac{1}{2}g^{y\alpha}\left(\partial_t g_{\alpha x} + \partial_x g_{\alpha t} - \partial_\alpha g_{tx}\right),$$

and the sum over $\alpha$ has only the two terms $\alpha = t$ and $\alpha = y$, since $g^{yx}$ and $g^{yz}$ vanish.
The $\alpha = t$ term is zero: $g_{tx} = 0$ kills the first and third pieces, and $g_{tt}$ is constant so $\partial_x g_{tt} = 0$.
The $\alpha = y$ term keeps one piece only,

$$\partial_t g_{yx} + \partial_x g_{yt} - \partial_y g_{tx} = 0 + \partial_x\left(-\frac{e^x}{2\omega^2}\right) - 0 = -\frac{e^x}{2\omega^2},$$

so

$$\Gamma^y{}_{tx} = \tfrac{1}{2}\left(4e^{-2x}\omega^2\right)\left(-\frac{e^x}{2\omega^2}\right) = -e^{-x},$$

which is the published value.

Here is $\Gamma^t{}_{tx}$ by hand.
The same two brackets appear, because the index pattern $tx$ is the same and only the raised index has changed,

$$\Gamma^t{}_{tx} = \tfrac{1}{2}g^{t\alpha}\left(\partial_t g_{\alpha x} + \partial_x g_{\alpha t} - \partial_\alpha g_{tx}\right) = \tfrac{1}{2}\left(-4e^{-x}\omega^2\right)\left(-\frac{e^x}{2\omega^2}\right) = 1,$$

which is the published value, and the exponentials cancelling against each other is what leaves the bare $1$.

One structural remark, which the equations of Step A3 make visible.
No symbol carries $\omega$.
The reason is that $\omega$ enters the metric only as the constant overall factor $1/(2\omega^2)$, and a constant rescaling of $g$ leaves $g^{-1}\partial g$ untouched, so it cannot reach the connection.
The parameter $\omega$ therefore sets the unit of length in which $ds^2$ is measured and does not appear in a single geodesic equation.

### A3. The four equations

The double sum runs over all sixteen pairs, and every symbol with two distinct lower indices is picked up twice.

For $\mu = t$, the symbols available are $\Gamma^t{}_{tx}$ and $\Gamma^t{}_{xy}$, each counted twice,

$$\Gamma^t{}_{\nu\rho}\dot{x}^\nu\dot{x}^\rho = 2\dot{t}\dot{x} + e^x\dot{x}\dot{y}.$$

For $\mu = x$, the symbols are $\Gamma^x{}_{ty}$, counted twice, and $\Gamma^x{}_{yy}$, counted once,

$$\Gamma^x{}_{\nu\rho}\dot{x}^\nu\dot{x}^\rho = e^x\dot{t}\dot{y} + \frac{e^{2x}}{2}\dot{y}^2.$$

For $\mu = y$, the only symbol is $\Gamma^y{}_{tx}$, counted twice,

$$\Gamma^y{}_{\nu\rho}\dot{x}^\nu\dot{x}^\rho = -2e^{-x}\dot{t}\dot{x}.$$

For $\mu = z$ there is nothing at all.

So the entry publishes

$$\ddot{t} + 2\dot{t}\dot{x} + e^x\dot{x}\dot{y} = 0,$$

$$\ddot{x} + e^x\dot{t}\dot{y} + \frac{e^{2x}}{2}\dot{y}^2 = 0,$$

$$\ddot{y} - 2e^{-x}\dot{t}\dot{x} = 0,$$

$$\ddot{z} = 0.$$

### A4. Three first integrals, and a check that uses two equations at once

Nothing in the metric depends on $t$, on $y$ or on $z$, so $\partial_t$, $\partial_y$ and $\partial_z$ are Killing vectors and each gives a constant of the motion.
Lowering the velocity on each in turn,

$$E \equiv -g_{t\mu}\dot{x}^\mu = \frac{1}{2\omega^2}\left(\dot{t} + e^x\dot{y}\right),$$

$$L \equiv g_{y\mu}\dot{x}^\mu = -\frac{e^x}{2\omega^2}\dot{t} - \frac{e^{2x}}{4\omega^2}\dot{y},$$

$$P \equiv g_{z\mu}\dot{x}^\mu = \frac{\dot{z}}{2\omega^2}.$$

That $P$ is constant is the fourth equation and nothing more.
That $E$ is constant is a genuine check, because it needs two of the equations at once and they have to cancel each other exactly.
Differentiating,

$$2\omega^2\frac{dE}{d\lambda} = \ddot{t} + e^x\dot{x}\dot{y} + e^x\ddot{y},$$

and substituting the first and third equations,

$$= \left(-2\dot{t}\dot{x} - e^x\dot{x}\dot{y}\right) + e^x\dot{x}\dot{y} + e^x\left(2e^{-x}\dot{t}\dot{x}\right) = -2\dot{t}\dot{x} + 2\dot{t}\dot{x} = 0.$$

The $e^x\dot{x}\dot{y}$ terms cancel between the first equation and the derivative of the cross term, and the $\dot{t}\dot{x}$ terms cancel between the first and the third, with the factor $e^{x}$ of the cross term meeting the factor $e^{-x}$ of $\Gamma^y{}_{tx}$.
Either the coefficient $2$ in the first equation or the coefficient $-2$ in the third would break it on its own.
The same substitution kills $dL/d\lambda$, with the $\dot{x}\dot{y}$ terms cancelling against each other and the $\dot{t}\dot{x}$ terms doing the same.

### A5. What the equations say

The dust is in free fall.
The Godel solution is sourced by pressureless dust comoving with $\partial_t$, and a pressureless fluid has to move on geodesics, since there is no pressure gradient available to push it off one.
The equations say so directly: setting $\dot{x} = \dot{y} = \dot{z} = 0$ leaves $\ddot{t} = 0$ in the first equation and nothing at all in the other three, because $\Gamma^\mu{}_{tt} = 0$ for every $\mu$.
The worldlines of constant $x$, $y$ and $z$ are therefore geodesics, affinely parametrised by $t$ itself, and since $g_{tt} = -1/(2\omega^2)$ the unit four velocity along them is $u^\mu = \sqrt{2}\,\omega\,\delta^\mu_t$.
That is the $u^\mu$ the entry's own `convention` field names when it writes $R_{\mu\nu} = 2\omega^2u_\mu u_\nu$.

The timelike direction that closes is not free.
Section A1 noted that $\partial_y$ is timelike everywhere.
The second equation says what it costs to follow it: a curve with $\dot{y}$ alone nonzero has

$$\ddot{x} = -\frac{e^{2x}}{2}\dot{y}^2 < 0,$$

so no orbit of $\partial_y$ is a geodesic, and anyone who insists on following one is being pushed toward smaller $x$ the whole way.
That is the family which closes on itself once the chart is wrapped into the cylindrical form, so the closed timelike curves Godel is famous for cost fuel: they have to be driven, and the second equation is where the bill appears.

The $z$ direction is flat and free.
The fourth equation is $\ddot{z} = 0$ with no coupling to anything, which is the geodesic statement of the fact that Godel's solution is a product of a three dimensional Lorentzian factor with a line.

---

## Part B. Reissner-Nordstrom

### B1. The line element

Write

$$\Delta \equiv r^2 - r_s r + r_q^2,$$

as the entry's `convention` field does, with $r_s = 2GM/c^2$ the Schwarzschild radius and $r_q^2 = Q^2G/(4\pi\epsilon_0c^4)$ the square of the charge radius.
Then

$$ds^2 = -\frac{\Delta}{r^2}dt^2 + \frac{r^2}{\Delta}dr^2 + r^2\left(d\theta^2 + \sin^2\theta\,d\phi^2\right),$$

which is the usual $-\left(1 - \frac{r_s}{r} + \frac{r_q^2}{r^2}\right)dt^2 + \ldots$ written over a common denominator, since $\Delta/r^2$ is exactly that bracket.
The horizons are the roots of $\Delta$,

$$r_\pm = \frac{r_s}{2} \pm \sqrt{\frac{r_s^2}{4} - r_q^2},$$

real when $r_s \geq 2r_q$.

The metric is diagonal, which makes the connection short, and the only function of $r$ in it is $\Delta/r^2$ and its reciprocal.
Its derivative is used three times below, so it is worth having once:

$$\partial_r\frac{\Delta}{r^2} = \partial_r\left(1 - \frac{r_s}{r} + \frac{r_q^2}{r^2}\right) = \frac{r_s}{r^2} - \frac{2r_q^2}{r^3} = \frac{r r_s - 2r_q^2}{r^3}.$$

The combination $r r_s - 2r_q^2$ that appears in every one of the entry's radial symbols is nothing but this derivative, cleared of its denominator.

### B2. The connection, with two symbols checked by hand

The entry publishes thirteen nonzero symbols, and the radial sector is

$$\Gamma^t{}_{tr} = \Gamma^t{}_{rt} = \frac{r r_s - 2r_q^2}{2r\Delta}, \qquad \Gamma^r{}_{tt} = \frac{\Delta\left(r r_s - 2r_q^2\right)}{2r^5},$$

$$\Gamma^r{}_{rr} = -\frac{r r_s - 2r_q^2}{2r\Delta}, \qquad \Gamma^r{}_{\theta\theta} = -\frac{\Delta}{r}, \qquad \Gamma^r{}_{\phi\phi} = -\frac{\Delta}{r}\sin^2\theta,$$

with the angular sector the usual $\Gamma^\theta{}_{r\theta} = \Gamma^\phi{}_{r\phi} = 1/r$, $\Gamma^\theta{}_{\phi\phi} = -\sin\theta\cos\theta$ and $\Gamma^\phi{}_{\theta\phi} = \cot\theta$ of any spherically symmetric chart with areal radius $r$.
The file prints $\Gamma^r{}_{\theta\theta}$ as $r_s - r - r_q^2/r$, which is $-\Delta/r$ expanded, and $\Gamma^r{}_{rr}$ over the cleared denominator $2r^3 + 2rr_q^2 - 2r^2r_s$, which is $2r\Delta$.

Here is $\Gamma^r{}_{tt}$ by hand.
The metric is diagonal and $g_{tt}$ depends on $r$ alone, so only one term of the definition survives,

$$\Gamma^r{}_{tt} = -\tfrac{1}{2}g^{rr}\,\partial_r g_{tt} = -\tfrac{1}{2}\cdot\frac{\Delta}{r^2}\cdot\partial_r\left(-\frac{\Delta}{r^2}\right) = \frac{\Delta}{2r^2}\cdot\frac{r r_s - 2r_q^2}{r^3} = \frac{\Delta\left(r r_s - 2r_q^2\right)}{2r^5},$$

using the derivative of Step B1, and that is the published value.

Here is $\Gamma^t{}_{tr}$ by hand.
The same derivative appears, now divided by $g_{tt}$ rather than multiplied by $g^{rr}$,

$$\Gamma^t{}_{tr} = \tfrac{1}{2}g^{tt}\,\partial_r g_{tt} = \tfrac{1}{2}\,\partial_r\ln\left(\frac{\Delta}{r^2}\right) = \frac{r^2}{2\Delta}\cdot\frac{r r_s - 2r_q^2}{r^3} = \frac{r r_s - 2r_q^2}{2r\Delta},$$

which is the published value.
The pairing is the one every static spherically symmetric chart has: $\Gamma^t{}_{tr}$ is half the logarithmic derivative of the redshift factor, and $\Gamma^r{}_{tt}$ is that same derivative weighted by the factor itself.

Setting $r_q = 0$ turns $\Delta$ into $r(r-r_s)$ and every symbol above into the Schwarzschild one the `schwarzschild/spherical` entry publishes, which is a third check and a free one.

### B3. The four equations

For $\mu = t$ the only symbol is $\Gamma^t{}_{tr}$, counted twice.
For $\mu = r$ the four symbols $\Gamma^r{}_{tt}$, $\Gamma^r{}_{rr}$, $\Gamma^r{}_{\theta\theta}$ and $\Gamma^r{}_{\phi\phi}$ are each counted once, since each has a repeated lower index.
For $\mu = \theta$ and $\mu = \phi$ the angular sector gives the standard pattern, with $\Gamma^\theta{}_{r\theta}$, $\Gamma^\phi{}_{r\phi}$ and $\Gamma^\phi{}_{\theta\phi}$ counted twice and $\Gamma^\theta{}_{\phi\phi}$ once.
So the entry publishes

$$\ddot{t} + \frac{r r_s - 2r_q^2}{r\left(r^2 - r r_s + r_q^2\right)}\dot{r}\dot{t} = 0,$$

$$\ddot{r} + \frac{\left(r^2 - r r_s + r_q^2\right)\left(r r_s - 2r_q^2\right)}{2r^5}\dot{t}^2 - \frac{r r_s - 2r_q^2}{2r\left(r^2 - r r_s + r_q^2\right)}\dot{r}^2 - \frac{r^2 - r r_s + r_q^2}{r}\dot{\theta}^2 - \frac{r^2 - r r_s + r_q^2}{r}\sin^2\theta\,\dot{\phi}^2 = 0,$$

$$\ddot{\theta} + \frac{2}{r}\dot{r}\dot{\theta} - \cos\theta\sin\theta\,\dot{\phi}^2 = 0,$$

$$\ddot{\phi} + \frac{2}{r}\dot{r}\dot{\phi} + 2\cot\theta\,\dot{\theta}\dot{\phi} = 0.$$

$\Delta$ is written out in full rather than abbreviated, so that each equation stands on its own and so that the machine check of Step 4 reads the same symbols the rest of the entry uses.
Setting $r_q = 0$ returns the `schwarzschild/spherical` equations term by term.

### B4. The first integrals and the equatorial plane

The two Killing vectors $\partial_t$ and $\partial_\phi$ give

$$E \equiv -g_{t\mu}\dot{x}^\mu = \frac{\Delta}{r^2}\dot{t}, \qquad L \equiv g_{\phi\mu}\dot{x}^\mu = r^2\sin^2\theta\,\dot{\phi}.$$

The first of these is the time equation, rearranged.
Differentiating,

$$\frac{dE}{d\lambda} = \frac{\Delta}{r^2}\ddot{t} + \dot{t}\dot{r}\,\partial_r\frac{\Delta}{r^2} = \frac{\Delta}{r^2}\left(\ddot{t} + \frac{r^2}{\Delta}\cdot\frac{r r_s - 2r_q^2}{r^3}\dot{r}\dot{t}\right) = \frac{\Delta}{r^2}\left(\ddot{t} + \frac{r r_s - 2r_q^2}{r\Delta}\dot{r}\dot{t}\right),$$

and the bracket is the published time equation exactly, so $E$ is constant if and only if that equation holds.
That is the cleanest possible confirmation of the coefficient: no factor of two is free in it.

The third equation is what makes the equatorial plane usable.
Its only term without a $\dot\theta$ is $-\cos\theta\sin\theta\,\dot{\phi}^2$, which vanishes at $\theta = \pi/2$, so a geodesic starting in that plane with $\dot\theta = 0$ has $\ddot\theta = 0$ and stays in it.
That is why the two integrals above, together with the normalisation of the velocity, are enough to reduce the problem to one radial equation, and why $L$ simplifies to $r^2\dot\phi$ in every textbook treatment.

### B5. Where the attraction reverses

The coefficient of $\dot{t}^2$ in the radial equation is the gravitational field a static observer feels, and it factorises,

$$\Gamma^r{}_{tt} = \frac{\Delta\left(r r_s - 2r_q^2\right)}{2r^5},$$

into the metric function $\Delta$ and the second factor $r r_s - 2r_q^2$, which is the one that carries the sign of the source.

Outside the outer horizon both factors are positive, so the term is positive, so $\ddot{r}$ is negative for a static particle, which is an attraction.
$\Delta > 0$ there by definition.
And $r r_s - 2r_q^2 > 0$ there because $r \geq r_+ \geq r_s/2$ gives $r r_s \geq r_s^2/2 \geq 2r_q^2$, the last step being the condition $r_s \geq 2r_q$ for horizons to exist at all, which the entry states.
Equality holds only in the extremal case $r_s = 2r_q$ evaluated on the horizon itself.

Inside, the second factor changes sign at

$$r_\ast = \frac{2r_q^2}{r_s},$$

and the inequality just used says $r_\ast \leq r_s/2 \leq r_+$, so the reversal is always hidden behind the outer horizon.

The place to read it is the inner static region $r < r_-$, where $\Delta$ has turned positive again and a static observer exists once more.
Since $r_+r_- = r_q^2$ and $r_+ + r_- = r_s$, the two radii compare as

$$\frac{r_\ast}{r_-} = \frac{2r_q^2/r_s}{r_q^2/r_+} = \frac{2r_+}{r_+ + r_-} > 1,$$

so $r_- < r_\ast$ strictly whenever the two horizons are distinct, and the whole of the inner static region lies below $r_\ast$.
There $\Delta > 0$ while $r r_s - 2r_q^2 < 0$, so $\Gamma^r{}_{tt} < 0$ and the radial equation gives $\ddot{r} > 0$ for a static particle: the charge term has overtaken the mass term and the field points outward.
That is the repulsion of the timelike Reissner-Nordstrom singularity, the reason a radially infalling particle is turned around before reaching $r = 0$ rather than crushed on it.
The whole of that statement is one sign in one coefficient of the second equation, and it is the physical reason the equation is worth publishing rather than leaving to the reader.

Setting $r_q = 0$ removes $r_\ast$ altogether, and the attraction never reverses, which is Schwarzschild.

---

## Part C. The Schwarzschild interior

### C1. The line element and two abbreviations

The entry publishes, on $0 \leq r \leq R$ with $r$ the areal radius and $R$ the stellar radius,

$$ds^2 = -\frac{1}{4}\left(\sqrt{1 - \frac{r^2r_s}{R^3}} - 3\sqrt{1 - \frac{r_s}{R}}\right)^2 dt^2 + \frac{dr^2}{1 - \dfrac{r^2r_s}{R^3}} + r^2\left(d\theta^2 + \sin^2\theta\,d\phi^2\right).$$

Two abbreviations shorten everything below and are used only in this file, never in the entry:

$$f(r) \equiv \sqrt{1 - \frac{r^2r_s}{R^3}}, \qquad h \equiv \sqrt{1 - \frac{r_s}{R}}.$$

The second is a constant and it is the first evaluated at the surface, $h = f(R)$.
In these,

$$g_{tt} = -\tfrac{1}{4}\left(f - 3h\right)^2, \qquad g_{rr} = \frac{1}{f^2}, \qquad g_{\theta\theta} = r^2, \qquad g_{\phi\phi} = r^2\sin^2\theta.$$

The derivative used throughout is

$$f' = \frac{d}{dr}\sqrt{1 - \frac{r^2r_s}{R^3}} = \frac{-r r_s/R^3}{f},$$

and $h' = 0$.

The bracket $3h - f$ is positive on the whole star exactly when the star is allowed to exist.
At the surface it is $3h - h = 2h > 0$ for any $r_s < R$.
At the centre it is $3h - 1$, which is positive when $9\left(1 - r_s/R\right) > 1$, that is when

$$r_s < \frac{8}{9}R,$$

or equivalently $R > \tfrac{9}{8}r_s$, which is the compactness limit the entry's own `history` field names.
At $R = \tfrac{9}{8}r_s$ the central value of $g_{tt}$ reaches zero, the central redshift and the central pressure both diverge, and no static uniform sphere exists past it; Buchdahl later showed the same number bounds every equation of state.
So $3h - f > 0$ everywhere on a star the entry describes, and every sign below is settled by that.

### C2. The connection, with two symbols checked by hand

Here is $\Gamma^t{}_{tr}$ by hand.
As in Part B it is half a logarithmic derivative,

$$\Gamma^t{}_{tr} = \tfrac{1}{2}\,\partial_r\ln\left(-g_{tt}\right) = \partial_r\ln\left|f - 3h\right| = \frac{f'}{f - 3h} = \frac{-r r_s/R^3}{f\left(f - 3h\right)} = \frac{r r_s}{R^3 f\left(3h - f\right)}.$$

Here is $\Gamma^r{}_{tt}$ by hand,

$$\Gamma^r{}_{tt} = -\tfrac{1}{2}g^{rr}\,\partial_r g_{tt} = -\tfrac{1}{2}f^2\cdot\left(-\tfrac{1}{2}\left(f - 3h\right)f'\right) = \frac{f^2\left(f - 3h\right)}{4}\cdot\frac{-r r_s/R^3}{f} = \frac{r r_s}{4R^3}f\left(3h - f\right).$$

The other two radial symbols come out the same way,

$$\Gamma^r{}_{rr} = \tfrac{1}{2}g^{rr}\,\partial_r g_{rr} = -\frac{f'}{f} = \frac{r r_s/R^3}{f^2} = \frac{r r_s}{R^3 - r^2r_s},$$

$$\Gamma^r{}_{\theta\theta} = -\tfrac{1}{2}g^{rr}\,\partial_r\left(r^2\right) = -rf^2 = r\left(\frac{r^2r_s}{R^3} - 1\right), \qquad \Gamma^r{}_{\phi\phi} = \Gamma^r{}_{\theta\theta}\sin^2\theta,$$

and the last two are exactly the forms the entry prints.
The angular sector is again the standard one of an areal radius chart.

The two hand values are tied to each other by the relation every static diagonal chart obeys,

$$\Gamma^r{}_{tt} = -g^{rr}g_{tt}\,\Gamma^t{}_{tr},$$

which follows at once from $\Gamma^t{}_{tr} = \partial_r g_{tt}/(2g_{tt})$ and $\Gamma^r{}_{tt} = -g^{rr}\partial_r g_{tt}/2$.
Here $-g^{rr}g_{tt} = f^2\left(3h-f\right)^2/4$, and multiplying $\Gamma^t{}_{tr}$ by it does give $\Gamma^r{}_{tt}$, so the two hand computations check each other.
The same relation holds in Part B, where $-g^{rr}g_{tt} = \Delta^2/r^4$ turns $\Gamma^t{}_{tr}$ into $\Gamma^r{}_{tt}$ term for term, and it is the reason the same bracket appears in both, once in a denominator and once in a numerator.

### C3. The four equations

Counting the double sum as before,

$$\ddot{t} + \frac{2 r r_s}{R^3\sqrt{1 - \dfrac{r^2 r_s}{R^3}}\left(3\sqrt{1 - \dfrac{r_s}{R}} - \sqrt{1 - \dfrac{r^2 r_s}{R^3}}\right)}\dot{r}\dot{t} = 0,$$

$$\ddot{r} + \frac{r r_s}{4R^3}\sqrt{1 - \frac{r^2 r_s}{R^3}}\left(3\sqrt{1 - \frac{r_s}{R}} - \sqrt{1 - \frac{r^2 r_s}{R^3}}\right)\dot{t}^2 + \frac{r r_s}{R^3 - r^2 r_s}\dot{r}^2 + r\left(\frac{r^2 r_s}{R^3} - 1\right)\dot{\theta}^2 + r\left(\frac{r^2 r_s}{R^3} - 1\right)\sin^2\theta\,\dot{\phi}^2 = 0,$$

$$\ddot{\theta} + \frac{2}{r}\dot{r}\dot{\theta} - \cos\theta\sin\theta\,\dot{\phi}^2 = 0,$$

$$\ddot{\phi} + \frac{2}{r}\dot{r}\dot{\phi} + 2\cot\theta\,\dot{\theta}\dot{\phi} = 0.$$

The first two are the hand values of Step C2 with $f$ and $h$ written out, the first doubled because its lower indices are distinct.
The angular pair is identical to Reissner-Nordstrom's and to Schwarzschild's, as it has to be, since those two equations see only $g_{\theta\theta} = r^2$ and $g_{\phi\phi} = r^2\sin^2\theta$.

### C4. Why these equations are printed in the radicals of the metric

The entry's own `christoffel` block prints the time symbol differently,

$$\Gamma^t{}_{tr} = \frac{r r_s}{-R^3 + r^2r_s + 3R\sqrt{\left(R - r_s\right)\left(R^3 - r^2r_s\right)}},$$

with the two radicals gathered under one root sign.
That is the same number as the form of Step C2, on the whole of the chart, and here is why.
The denominator rewrites as

$$-R^3 + r^2r_s = -R^3\left(1 - \frac{r^2r_s}{R^3}\right) = -R^3f^2,$$

$$\sqrt{\left(R - r_s\right)\left(R^3 - r^2r_s\right)} = \sqrt{Rh^2\cdot R^3f^2} = R^2fh,$$

the second line using $R > 0$ and $f, h \geq 0$, which the domain $0 \leq r \leq R$ and the bound $r_s < 8R/9$ of Step C1 guarantee.
So the denominator is $-R^3f^2 + 3R^3fh = R^3f\left(3h - f\right)$, and the two forms agree.

The published geodesics use the $\sqrt{1 - \cdot}$ form rather than the gathered one, and the reason is worth recording rather than leaving as a matter of taste.
The gathered form needs the step $\sqrt{AB} = \sqrt{A}\sqrt{B}$, which is only true when neither factor is negative.
A computer algebra system carrying $R$, $r$ and $r_s$ as real symbols will not take that step, since it is false for real arguments in general, so it cannot reduce the difference of the two forms to zero.
`_tools/derivations/verify_metrics.py` carries exactly such symbols, and it rebuilds the connection from the line element before comparing, so the geodesic equations have to be in a form whose difference from the rebuilt one it can close.
Written in $\sqrt{1 - r^2r_s/R^3}$ and $\sqrt{1 - r_s/R}$, which are the radicals the line element is itself written in, the difference collapses without any assumption about signs and the comparison passes.
Written in the gathered form, it does not, and the checker reports a disagreement for an equation that is perfectly correct.

The identity was confirmed three ways before choosing: symbolically under a simplifier told that the factors are positive, numerically at several interior points to forty digits, and by the hand derivation of Step C2, which reaches the $\sqrt{1 - \cdot}$ form straight from the line element without passing through the gathered one at all.

This is not a hypothetical limitation, and the entry already shows it.
Run the checker against `interior_schwarzschild/spherical` and it reports the three gathered symbols of the `christoffel` block, $\Gamma^t{}_{tr}$, $\Gamma^t{}_{rt}$ and $\Gamma^r{}_{tt}$, as disagreements in both variants, along with the Riemann, Ricci and Einstein components built on them.
Those reports are about the form and not about the numbers: each of those symbols is the right value, as Step C2 derives independently and as a forty digit evaluation confirms.
The `christoffel` block was not touched here, because the curvature of this entry is being corrected by separate work and moving it would collide.
What this file could control is the form of its own twelve equations, and choosing the radicals the line element already uses is what lets the geodesic comparison close where the blocks above it do not.

### C5. The surface, checked against the exterior

The entry's `convention` field says the solution matches the exterior Schwarzschild metric at $r = R$, and the connection is where that can be checked directly.
Evaluating the interior symbols at $r = R$, where $f = h$ and so $3h - f = 2h$ and $h^2 = 1 - r_s/R$:

| symbol | interior at $r = R$ | exterior Schwarzschild at $r = R$ | |
| --- | --- | --- | --- |
| $\Gamma^t{}_{tr}$ | $\dfrac{r_s}{2R\left(R - r_s\right)}$ | $\dfrac{r_s}{2R\left(R-r_s\right)}$ | matches |
| $\Gamma^r{}_{tt}$ | $\dfrac{r_s\left(R - r_s\right)}{2R^3}$ | $\dfrac{r_s\left(R-r_s\right)}{2R^3}$ | matches |
| $\Gamma^r{}_{\theta\theta}$ | $r_s - R$ | $r_s - R$ | matches |
| $\Gamma^r{}_{rr}$ | $\dfrac{r_s}{R\left(R - r_s\right)}$ | $-\dfrac{r_s}{2R\left(R - r_s\right)}$ | jumps |

The exterior column is the `schwarzschild/spherical` entry's own connection at $r = R$, which the reader can check against `schwarzschild.json` directly.

Three of the four match, and the fourth is expected to jump.
The junction conditions constrain the induced metric on the surface and its extrinsic curvature, which between them involve $g_{tt}$, $g_{\theta\theta}$ and their radial derivatives, and those are precisely the three that match.
They do not involve $\partial_r g_{rr}$, and $\partial_r g_{rr}$ is where the uniform density model has its discontinuity: the density drops from a constant to zero across the surface, so the mass function has a kink there, so $g_{rr}$ is continuous and its radial derivative is not.
$\Gamma^r{}_{rr}$ is built out of that derivative alone, and it is the only symbol of the four that is.

The geodesic equations carry this in one term.
A particle crossing the surface has its $\dot{t}\dot{r}$, $\dot{t}^2$ and $\dot\theta^2$ coefficients continuous and its $\dot{r}^2$ coefficient discontinuous, which is the exact geometric signature of a star with a sharp edge, and is an artefact of the uniform density idealisation rather than of the chart.

### C6. What the radial equation says

A free particle falls inward everywhere in the star, and at the centre it does not fall at all.
A particle momentarily at rest at radius $r$ has $\dot{r} = \dot\theta = \dot\phi = 0$, so the second equation leaves

$$\ddot{r} = -\frac{r r_s}{4R^3}f\left(3h - f\right)\dot{t}^2,$$

and every factor on the right is positive on $0 < r \leq R$ by Step C1, so $\ddot{r} < 0$ throughout the star.
At $r = 0$ the coefficient vanishes, so the centre is a free fall point, which spherical symmetry demands and which is a check on the factor of $r$ in the numerator.

The Newtonian limit is exact and it is the one to expect.
For $r_s \ll R$ both radicals go to $1$, the bracket goes to $2$, and

$$\Gamma^r{}_{tt} \longrightarrow \frac{r r_s}{2R^3} = \frac{GMr}{c^2R^3}$$

to first order in $r_s$.
For a slowly moving particle with $\lambda = \tau$ the chart velocity $\dot{t} = d(ct)/d\tau$ tends to $c$, so

$$\frac{d^2r}{d\tau^2} \longrightarrow -\frac{GMr}{R^3},$$

which is the field inside a uniform Newtonian sphere, linear in $r$ and vanishing at the centre.
That the factors of $c$ land correctly is a consequence of the chart convention of Step 1 and would fail by exactly one factor of $c$ on any other reading of the dots.

---

## Step 3. Every published equation is dimensionally consistent

A geodesic equation is measured against its own second derivative: every term of it has to carry $[x^\mu]/\lambda^2$, where $[x^\mu]$ is the dimension of the chart coordinate the equation belongs to.
The dots are the chart velocities of Step 1, so each carries $[x^\mu]/\lambda$ and each double dot $[x^\mu]/\lambda^2$.

Godel carries all four coordinates as pure numbers, so every dot is $1/\lambda$, every double dot is $1/\lambda^2$, and every coefficient has to be dimensionless.

| equation | term | carries |
| --- | --- | --- |
| $t$ | $\ddot{t}$, $2\dot{t}\dot{x}$, $e^x\dot{x}\dot{y}$ | $1/\lambda^2$ |
| $x$ | $\ddot{x}$, $e^x\dot{t}\dot{y}$, $\tfrac{1}{2}e^{2x}\dot{y}^2$ | $1/\lambda^2$ |
| $y$ | $\ddot{y}$, $2e^{-x}\dot{t}\dot{x}$ | $1/\lambda^2$ |
| $z$ | $\ddot{z}$ | $1/\lambda^2$ |

The exponentials are dimensionless because $x$ is, which is what forced $x$ to be dimensionless in the first place, and $\omega$ does not appear at all by Step A2, so there is nothing else to balance.

Reissner-Nordstrom and the Schwarzschild interior carry $t$ and $r$ as lengths, beside the length parameters $r_s$, $r_q$ and $R$, with $\theta$ and $\phi$ dimensionless.
So $\dot{t}$ and $\dot{r}$ are $L/\lambda$, $\dot\theta$ and $\dot\phi$ are $1/\lambda$, and the two equations that matter read

| equation | term | carries | check |
| --- | --- | --- | --- |
| $t$ | $\ddot{t}$ | $L/\lambda^2$ | |
| | $\dfrac{r r_s - 2r_q^2}{r\Delta}\dot{r}\dot{t}$ | $L/\lambda^2$ | $\dfrac{L^2}{L\cdot L^2}\cdot\dfrac{L^2}{\lambda^2}$ |
| $r$ | $\ddot{r}$ | $L/\lambda^2$ | |
| | $\dfrac{\Delta\left(r r_s - 2r_q^2\right)}{2r^5}\dot{t}^2$ | $L/\lambda^2$ | $\dfrac{L^2\cdot L^2}{L^5}\cdot\dfrac{L^2}{\lambda^2}$ |
| | $\dfrac{\Delta}{r}\dot\theta^2$ | $L/\lambda^2$ | $\dfrac{L^2}{L}\cdot\dfrac{1}{\lambda^2}$ |
| $\theta$ | $\ddot\theta$, $\dfrac{2}{r}\dot{r}\dot\theta$, $\cos\theta\sin\theta\,\dot\phi^2$ | $1/\lambda^2$ | $\dfrac{1}{L}\cdot\dfrac{L}{\lambda}\cdot\dfrac{1}{\lambda}$ |

The interior Schwarzschild coefficients balance the same way.
Its radicals are dimensionless, since $r^2r_s/R^3$ and $r_s/R$ are, so the time coefficient is $L^2/L^3 = 1/L$ and the $\dot{t}^2$ coefficient of the radial equation is $L^2/L^3 = 1/L$ as well, each meeting an $L^2/\lambda^2$ from the pair of dots.
The $\dot{r}^2$ coefficient is $L^2/L^3 = 1/L$ and the $\dot\theta^2$ coefficient is $L\cdot 1 = L$, which are the same two patterns the table above shows.

The angular equations of the last two entries are identical and both balance at $1/\lambda^2$, which is the check that catches a missing factor of $r$: a term $\dot{r}\dot\theta/r$ is $\left(1/L\right)\left(L/\lambda\right)\left(1/\lambda\right) = 1/\lambda^2$, and the same term without its $1/r$ would be a length over $\lambda^2$ and could not be added to $\ddot\theta$.

All twelve equations were put through this check mechanically, parsed out of the JSON rather than copied by hand, and every term balances.

---

## Step 4. What the three entries publish, and how to check it

| entry | field | value |
| --- | --- | --- |
| `godel/cartesian` | `geodesics` | four equations, Step A3 |
| `rn_metric/spherical` | `geodesics` | four equations, Step B3 |
| `interior_schwarzschild/spherical` | `geodesics` | four equations, Step C3 |

Each entry's `convention` field gained one sentence naming the chart, giving the geodesic equation in the form of Step 1, and fixing the dots as derivatives of that chart, in the way `kasner.json` and `vaidya.json` do.
Nothing else in the three files was touched.
In particular every curvature block is byte for byte what it was, since all three were being corrected by separate work at the same time as this.

To check:

    python3 -m venv /tmp/mfs-venv && /tmp/mfs-venv/bin/pip install sympy
    /tmp/mfs-venv/bin/python _tools/derivations/verify_metrics.py --system godel/cartesian
    /tmp/mfs-venv/bin/python _tools/derivations/verify_metrics.py --system rn_metric/spherical
    /tmp/mfs-venv/bin/python _tools/derivations/verify_metrics.py --system interior_schwarzschild/spherical

The first two take a few seconds each and the third takes about twelve minutes, almost all of it in the curvature comparison rather than in the geodesics.
The geodesic comparison itself is well under a second for all three, and it reports no disagreement for any of them.
The other blocks do: at the time of writing the checker reports four disagreements in Godel's `einstein_tensor` and a long list in the Schwarzschild interior's curvature, most of the second list being the form problem of Step C4.
Those are the subject of separate work and none of them is touched here.
Reissner-Nordstrom comes back clean throughout.
Each run rebuilds the connection from the entry's own line element and measures every published equation against

$$\ddot{x}^\mu + \Gamma^\mu{}_{\nu\rho}\dot{x}^\nu\dot{x}^\rho,$$

so it is an independent check of the equations rather than a restatement of the file.
The dimensional pass of Step 3 is the same script with `--dimensions-only`, which runs over the whole collection in about eight seconds.
