# The de Sitter spacetime

We work out de Sitter space in two coordinate systems, the static patch and the flat slicing.
Every component is derived in order, from the line element down to the geodesic equations.
We leave nothing as an exercise and assert nothing we do not compute.

The companion script `_tools/derivations/verify_metrics.py` does the same work in sympy and compares it against the stored components, so every line of the algebra is checkable by hand and by machine independently.

de Sitter is the first of the spacetimes that is a vacuum and yet has a curvature the Ricci tensor carries.
Every other vacuum among them empties the Ricci tensor and leaves the whole of the curvature in Weyl.
This one does the opposite, and does it completely: its Ricci tensor is $R_{\mu\nu} = \Lambda g_{\mu\nu}$ in Step 9, and its Weyl tensor is empty in every slot in Step 11.
That is what maximal symmetry costs and what it buys, and it is why the Weyl tensor of Step 11 is computed rather than copied.

---

## Step 1. Conventions

The signature is $(-,+,+,+)$ and the Christoffel symbols are those of the Levi-Civita connection,

$$\Gamma^\mu_{\nu\rho} = \tfrac{1}{2}g^{\mu\alpha}\left(\partial_\nu g_{\alpha\rho} + \partial_\rho g_{\alpha\nu} - \partial_\alpha g_{\nu\rho}\right),$$

with the lowered symbol $\Gamma_{\mu\nu\rho} = g_{\mu\alpha}\Gamma^\alpha_{\nu\rho}$.
The Riemann tensor is

$$R^\mu{}_{\nu\rho\sigma} = \partial_\rho \Gamma^\mu_{\nu\sigma} - \partial_\sigma \Gamma^\mu_{\nu\rho} + \Gamma^\mu_{\rho\lambda}\Gamma^\lambda_{\nu\sigma} - \Gamma^\mu_{\sigma\lambda}\Gamma^\lambda_{\nu\rho},$$

which fixes the sign of every Riemann component.

The Ricci tensor is the standard contraction,

$$R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu},$$

the contraction every spacetime has used since 2026-09-18.
Here the choice is not free, as it was for the vacuum spacetimes whose Ricci tensor vanishes either way.
This spacetime has a nonzero Ricci tensor, so the two contractions differ in every slot of it, of the Einstein tensor and of the Ricci scalar, and which sign the components carry, and why that is the one that reads correctly, comes out in Step 9.

The chart is the one whose zeroth coordinate is

$$x^0 = ct.$$

The index is written with the bare letter $t$, but the component written against it is a component of that chart.
We compute everything from Step 4 onward directly in it.
The dots in the geodesic equations are velocities of that same chart, so $\dot{t}$ means $d(ct)/d\lambda$, and it is on this reading that every term balances in Step 19.

The one parameter of the static patch is the cosmological constant $\Lambda$, taken positive, carrying an inverse length squared.
It sets one length,

$$\ell = \sqrt{\frac{3}{\Lambda}},$$

which is the horizon radius of Step 13 and the Hubble radius $c/H$ of Step 15 at once, since the Hubble rate of the flat slicing is $H = c\sqrt{\Lambda/3} = c/\ell$.
The abbreviation

$$f = 1 - \frac{\Lambda r^2}{3} = \frac{3 - \Lambda r^2}{3} = 1 - \frac{r^2}{\ell^2}$$

runs through Steps 3 to 14.
It is shorthand for the derivation only, and never enters the components themselves.
Every stored component is written out in $r$ and $\Lambda$, because the checker reads each one symbol by symbol and has no way to be told what an abbreviation means.

---

## Step 2. The field equation with a cosmological constant, and what it forces

The field equation this spacetime solves is the empty one with a $\Lambda$ term,

$$G_{\mu\nu} + \Lambda g_{\mu\nu} = 0, \qquad G_{\mu\nu} = R_{\mu\nu} - \tfrac{1}{2}Rg_{\mu\nu}.$$

Its content can be read off before any metric is written down.
Take the trace with $g^{\mu\nu}$ in four dimensions:

$$R - \tfrac{1}{2}R\cdot 4 + 4\Lambda = -R + 4\Lambda = 0,$$

so

$$R = 4\Lambda.$$

Put that back into the equation:

$$R_{\mu\nu} = \tfrac{1}{2}Rg_{\mu\nu} - \Lambda g_{\mu\nu} = 2\Lambda g_{\mu\nu} - \Lambda g_{\mu\nu} = \Lambda g_{\mu\nu},$$

and therefore

$$G_{\mu\nu} = \Lambda g_{\mu\nu} - 2\Lambda g_{\mu\nu} = -\Lambda g_{\mu\nu}.$$

Those three identities hold for any solution of the $\Lambda$ vacuum equation in four dimensions, whatever its symmetry.
Both charts satisfy them, and each of them is recovered from the metric rather than from this argument in Steps 9, 10 and 16.
The word vacuum is doing narrower work than usual here: there is no matter, so $T_{\mu\nu} = 0$, but the Ricci tensor is not zero, because $\Lambda$ sits on the geometry side of the equation and is a curvature in its own right.

What is special about de Sitter, and does not follow from the field equation alone, is maximal symmetry: it is the solution whose Riemann tensor is built from the metric and nothing else,

$$R_{\mu\nu\rho\sigma} = \frac{\Lambda}{3}\left(g_{\mu\rho}g_{\nu\sigma} - g_{\mu\sigma}g_{\nu\rho}\right).$$

Four components of the static patch yield that form in Step 8, and the flat slicing yields the same form with the same $\Lambda$ in Step 16.

---

## Step 3. The static patch line element

The primary chart is the static one,

$$ds^2 = -\left(1 - \frac{\Lambda r^2}{3}\right)c^2dt^2 + \frac{dr^2}{1 - \dfrac{\Lambda r^2}{3}} + r^2\left(d\theta^2 + \sin^2\theta\,d\phi^2\right),$$

that is $ds^2 = -fc^2dt^2 + f^{-1}dr^2 + r^2d\Omega^2$.

Two things about it matter before any algebra.

The coordinate $r$ is the areal radius, exactly as in Schwarzschild: the sphere at fixed $t$ and $r$ has area $4\pi r^2$, since the angular part of the line element is $r^2d\Omega^2$ untouched.
That is a definition of $r$, not a result, and it is what makes the horizon radius of Step 13 a statement about an area rather than about a distance.

The chart is static only where $f > 0$.
At $r = \ell$ the coefficient of $dt^2$ vanishes and the coefficient of $dr^2$ blows up, and past it the two swap sign, so $t$ stops being a time.
Nothing happens to the curvature there, by Step 13.

The form is the same one Schwarzschild wears, with $f = 1 - r_s/r$ replaced by $f = 1 - r^2/\ell^2$.
Every Christoffel symbol can be written once in terms of $f$ and $f'$ and then specialised, and the specialisation is

$$f' = \frac{df}{dr} = -\frac{2\Lambda r}{3}, \qquad f'' = -\frac{2\Lambda}{3}.$$

The second derivative is a constant, which is the whole of the difference between this spacetime and a black hole.

---

## Step 4. The metric matrix, its determinant and its inverse

In the chart $x^0 = ct$ the time term is $-f(c\,dt)^2$, so the matrix of the metric in the coordinate order $(t, r, \theta, \phi)$ is

$$g_{\mu\nu} = \begin{pmatrix} -f & 0 & 0 & 0 \\ 0 & f^{-1} & 0 & 0 \\ 0 & 0 & r^2 & 0 \\ 0 & 0 & 0 & r^2\sin^2\theta \end{pmatrix},$$

which written out is

$$g_{tt} = -\left(1 - \frac{\Lambda r^2}{3}\right), \qquad g_{rr} = \frac{3}{3 - \Lambda r^2}, \qquad g_{\theta\theta} = r^2, \qquad g_{\phi\phi} = r^2\sin^2\theta.$$

The chart component $g_{tt}$ is dimensionless because the $c^2$ of the line element has been absorbed into $(c\,dt)^2$.

The determinant is

$$\det g = (-f)\cdot f^{-1}\cdot r^2\cdot r^2\sin^2\theta = -r^4\sin^2\theta,$$

so $\sqrt{-g} = r^2\sin\theta$.
The two metric functions cancel each other exactly, which is a property of every metric of the form $-f\,,f^{-1}$ and the reason the volume element of the static patch is the flat one.

The matrix is diagonal, so the inverse is the reciprocal slot by slot,

$$g^{tt} = -\frac{3}{3 - \Lambda r^2}, \qquad g^{rr} = 1 - \frac{\Lambda r^2}{3}, \qquad g^{\theta\theta} = \frac{1}{r^2}, \qquad g^{\phi\phi} = \frac{1}{r^2\sin^2\theta},$$

which is the inverse metric in full.

---

## Step 5. The Christoffel symbols with an upper index

Every symbol comes from the connection formula of Step 1 with the diagonal metric of Step 4.

The two that mix $t$ and $r$:

$$\Gamma^t{}_{tr} = \Gamma^t{}_{rt} = \tfrac{1}{2}g^{tt}\partial_r g_{tt} = \tfrac{1}{2}\left(-\frac{1}{f}\right)(-f') = \frac{f'}{2f} = -\frac{\Lambda r}{3 - \Lambda r^2},$$

$$\Gamma^r{}_{tt} = -\tfrac{1}{2}g^{rr}\partial_r g_{tt} = -\tfrac{1}{2}f\cdot(-f') = \frac{ff'}{2} = -\frac{\Lambda r\left(3 - \Lambda r^2\right)}{9}.$$

In the first, $f'/(2f)$ becomes $(-2\Lambda r/3)/(2(3-\Lambda r^2)/3)$, and the threes cancel.
In the second, $ff'/2$ becomes $\left((3-\Lambda r^2)/3\right)\left(-2\Lambda r/3\right)/2$, which is where the nine in the denominator comes from.

The purely radial one:

$$\Gamma^r{}_{rr} = \tfrac{1}{2}g^{rr}\partial_r g_{rr} = \tfrac{1}{2}f\cdot\partial_r\!\left(\frac{1}{f}\right) = \tfrac{1}{2}f\cdot\left(-\frac{f'}{f^2}\right) = -\frac{f'}{2f} = \frac{\Lambda r}{3 - \Lambda r^2}.$$

It is minus $\Gamma^t{}_{tr}$, which is again the $-f\,,f^{-1}$ pairing.

The two that push the spheres outward:

$$\Gamma^r{}_{\theta\theta} = -\tfrac{1}{2}g^{rr}\partial_r g_{\theta\theta} = -\tfrac{1}{2}f\cdot 2r = -rf = -\frac{r\left(3 - \Lambda r^2\right)}{3},$$

$$\Gamma^r{}_{\phi\phi} = -rf\sin^2\theta = -\frac{r\left(3 - \Lambda r^2\right)\sin^2\theta}{3}.$$

And the four the round sphere carries whatever the radial functions are:

$$\Gamma^\theta{}_{r\theta} = \Gamma^\theta{}_{\theta r} = \Gamma^\phi{}_{r\phi} = \Gamma^\phi{}_{\phi r} = \frac{1}{r},$$

$$\Gamma^\theta{}_{\phi\phi} = -\cos\theta\sin\theta, \qquad \Gamma^\phi{}_{\theta\phi} = \Gamma^\phi{}_{\phi\theta} = \cot\theta.$$

That is thirteen nonzero symbols counting both orderings of each mixed pair, and it is all of the Christoffel symbols with an upper index.
Not one of them carries a factor of $c$: the metric has no $c$ in this chart and does not depend on $t$, so nothing can produce one.

---

## Step 6. The lowered Christoffel symbols

Lowering the first index multiplies by the diagonal metric component that matches it.

$$\Gamma_{ttr} = \Gamma_{trt} = g_{tt}\Gamma^t{}_{tr} = (-f)\frac{f'}{2f} = -\frac{f'}{2} = \frac{\Lambda r}{3},$$

$$\Gamma_{rtt} = g_{rr}\Gamma^r{}_{tt} = \frac{1}{f}\cdot\frac{ff'}{2} = \frac{f'}{2} = -\frac{\Lambda r}{3},$$

$$\Gamma_{rrr} = g_{rr}\Gamma^r{}_{rr} = \frac{1}{f}\cdot\left(-\frac{f'}{2f}\right) = -\frac{f'}{2f^2} = \frac{3\Lambda r}{\left(3 - \Lambda r^2\right)^2}.$$

The last one is where the $f^{-2}$ of the denominator comes from, and it is the only lowered symbol whose form differs from the corresponding one with an upper index by more than a sign.

The rest lose their metric functions entirely:

$$\Gamma_{r\theta\theta} = \frac{1}{f}(-rf) = -r, \qquad \Gamma_{r\phi\phi} = -r\sin^2\theta,$$

$$\Gamma_{\theta r\theta} = \Gamma_{\theta\theta r} = r^2\cdot\frac{1}{r} = r, \qquad \Gamma_{\phi r\phi} = \Gamma_{\phi\phi r} = r^2\sin^2\theta\cdot\frac{1}{r} = r\sin^2\theta,$$

$$\Gamma_{\theta\phi\phi} = r^2(-\cos\theta\sin\theta) = -r^2\cos\theta\sin\theta, \qquad \Gamma_{\phi\theta\phi} = \Gamma_{\phi\phi\theta} = r^2\sin^2\theta\cot\theta = r^2\cos\theta\sin\theta.$$

That is every Christoffel symbol with every index lowered, thirteen symbols again.

---

## Step 7. Four components of the Riemann tensor, computed

The Riemann tensor of a static spherically symmetric metric has four independent components up to the symmetries, one for each plane the geometry distinguishes.
We compute all four from the definition of Step 1, and by the argument of Step 8 they are enough.

**The $tr$ plane.**
With $\mu = t$, $\nu = r$, $\rho = t$, $\sigma = r$, and using $\Gamma^t{}_{rr} = 0$,

$$R^t{}_{rtr} = -\partial_r\Gamma^t{}_{rt} + \Gamma^t{}_{t\lambda}\Gamma^\lambda{}_{rr} - \Gamma^t{}_{r\lambda}\Gamma^\lambda{}_{rt}.$$

The derivative is $\partial_r\left(f'/2f\right) = f''/2f - (f')^2/2f^2$.
The first sum has one surviving term, $\Gamma^t{}_{tr}\Gamma^r{}_{rr} = (f'/2f)(-f'/2f) = -(f')^2/4f^2$.
The second has one, $\Gamma^t{}_{rt}\Gamma^t{}_{rt} = (f')^2/4f^2$.
Adding,

$$R^t{}_{rtr} = -\frac{f''}{2f} + \frac{(f')^2}{2f^2} - \frac{(f')^2}{4f^2} - \frac{(f')^2}{4f^2} = -\frac{f''}{2f} = \frac{\Lambda}{3 - \Lambda r^2},$$

using $f'' = -2\Lambda/3$ and $f = (3-\Lambda r^2)/3$.
The terms in $(f')^2$ cancel exactly, which they do for any $f$, so this component is $-f''/2f$ and nothing else.

**A $t$ and angle plane.**
With $\mu = t$, $\nu = \theta$, $\rho = t$, $\sigma = \theta$, only one product survives:

$$R^t{}_{\theta t\theta} = \Gamma^t{}_{tr}\Gamma^r{}_{\theta\theta} = \frac{f'}{2f}\cdot(-rf) = -\frac{rf'}{2} = \frac{\Lambda r^2}{3}.$$

**An $r$ and angle plane.**

$$R^r{}_{\theta r\theta} = \partial_r\Gamma^r{}_{\theta\theta} + \Gamma^r{}_{rr}\Gamma^r{}_{\theta\theta} - \Gamma^r{}_{\theta\theta}\Gamma^\theta{}_{\theta r} = \left(-f - rf'\right) + \left(-\frac{f'}{2f}\right)(-rf) - (-rf)\frac{1}{r},$$

$$R^r{}_{\theta r\theta} = -f - rf' + \frac{rf'}{2} + f = -\frac{rf'}{2} = \frac{\Lambda r^2}{3}.$$

It equals the previous one, which is the first hint of the symmetry.

**The sphere itself.**

$$R^\theta{}_{\phi\theta\phi} = \partial_\theta\Gamma^\theta{}_{\phi\phi} + \Gamma^\theta{}_{\theta r}\Gamma^r{}_{\phi\phi} - \Gamma^\theta{}_{\phi\phi}\Gamma^\phi{}_{\phi\theta}.$$

The derivative is $\partial_\theta(-\cos\theta\sin\theta) = \sin^2\theta - \cos^2\theta$, the second term is $(1/r)(-rf\sin^2\theta) = -f\sin^2\theta$, and the third is $-(-\cos\theta\sin\theta)(\cot\theta) = \cos^2\theta$.
Adding,

$$R^\theta{}_{\phi\theta\phi} = \sin^2\theta - \cos^2\theta - f\sin^2\theta + \cos^2\theta = (1-f)\sin^2\theta = \frac{\Lambda r^2\sin^2\theta}{3}.$$

The $\cos^2\theta$ terms cancel, and what is left is the curvature of a sphere of radius $r$ corrected by $f$, which for this $f$ is $\Lambda r^2/3$ exactly.

---

## Step 8. The maximally symmetric form, and the rest of the Riemann tensor

Each of the four results of Step 7 is reproduced by the single expression

$$R^\mu{}_{\nu\rho\sigma} = \frac{\Lambda}{3}\left(\delta^\mu{}_\rho g_{\nu\sigma} - \delta^\mu{}_\sigma g_{\nu\rho}\right).$$

Checking them in turn:

$$\frac{\Lambda}{3}\delta^t{}_t g_{rr} = \frac{\Lambda}{3}\cdot\frac{1}{f} = \frac{\Lambda}{3 - \Lambda r^2}, \qquad \frac{\Lambda}{3}\delta^t{}_t g_{\theta\theta} = \frac{\Lambda r^2}{3},$$

$$\frac{\Lambda}{3}\delta^r{}_r g_{\theta\theta} = \frac{\Lambda r^2}{3}, \qquad \frac{\Lambda}{3}\delta^\theta{}_\theta g_{\phi\phi} = \frac{\Lambda r^2\sin^2\theta}{3},$$

which are the four values computed.
Since a static spherically symmetric Riemann tensor has no independent components beyond those four, and that expression has the right symmetries, it is the whole tensor.

Lowering the first index gives the form quoted in Step 2,

$$R_{\mu\nu\rho\sigma} = \frac{\Lambda}{3}\left(g_{\mu\rho}g_{\nu\sigma} - g_{\mu\sigma}g_{\nu\rho}\right),$$

and every one of the twenty four nonzero components, with the first index up or with every index lowered, is one slot of it.
Three examples, one per pattern:

$$R_{trtr} = \frac{\Lambda}{3}g_{tt}g_{rr} = \frac{\Lambda}{3}(-f)\frac{1}{f} = -\frac{\Lambda}{3},$$

$$R_{t\theta t\theta} = \frac{\Lambda}{3}g_{tt}g_{\theta\theta} = -\frac{\Lambda r^2 f}{3} = -\frac{\Lambda r^2\left(3-\Lambda r^2\right)}{9}, \qquad R_{\theta\phi\theta\phi} = \frac{\Lambda}{3}g_{\theta\theta}g_{\phi\phi} = \frac{\Lambda r^4\sin^2\theta}{3}.$$

The first is a constant, with both metric functions cancelled, and it is the cleanest statement that the curvature does not vary from point to point.

The antisymmetry $R_{\mu\nu\rho\sigma} = -R_{\mu\nu\sigma\rho}$ is what pairs the components: every slot appears twice, once with the last two indices in each order and with opposite signs, which is how the count reaches twenty four in a tensor with six independent planes.

---

## Step 9. The Ricci tensor, which is not zero

Contract the form of Step 8 on the first index, as Step 1 requires:

$$R_{\nu\sigma} = R^\mu{}_{\nu\mu\sigma} = \frac{\Lambda}{3}\left(\delta^\mu{}_\mu g_{\nu\sigma} - \delta^\mu{}_\sigma g_{\nu\mu}\right) = \frac{\Lambda}{3}\left(4g_{\nu\sigma} - g_{\nu\sigma}\right) = \Lambda g_{\nu\sigma}.$$

This is the identity Step 2 got from the field equation, now recovered from the metric.
Slot by slot, in the chart,

$$R_{tt} = -\Lambda f = -\frac{\Lambda\left(3 - \Lambda r^2\right)}{3}, \qquad R_{rr} = \frac{\Lambda}{f} = \frac{3\Lambda}{3 - \Lambda r^2}, \qquad R_{\theta\theta} = \Lambda r^2, \qquad R_{\phi\phi} = \Lambda r^2\sin^2\theta,$$

and with one index raised the metric disappears entirely,

$$R^\mu{}_\nu = \Lambda\delta^\mu{}_\nu,$$

so all four diagonal components of the Ricci tensor with one index up are the bare $\Lambda$.
That form shows the isotropy directly: a curvature that is a multiple of the identity in an orthonormal frame is the definition of isotropy, and it holds at every point of the chart, including at the horizon where the other two index positions have a factor that misbehaves.

The contraction matters here in a way it did not for the vacuum spacetimes.
Contracting on the last index instead would give $R_{\mu\nu} = -\Lambda g_{\mu\nu}$, $R = -4\Lambda$ and $G_{\mu\nu} = +\Lambda g_{\mu\nu}$, and the field equation would have to be written $G_{\mu\nu} = \Lambda g_{\mu\nu}$ to compensate.
On the standard contraction the sign works out the way the rest of physics reads it: a positive $\Lambda$ acts as a fluid of positive energy density $\rho_\Lambda = \Lambda c^4/(8\pi G)$ and pressure $p_\Lambda = -\rho_\Lambda$, and $G_{tt}$ comes out positive where that density is, exactly as it does for ordinary matter in FRW.

---

## Step 10. The Ricci scalar and the Einstein tensor

The scalar is the trace of Step 9,

$$R = g^{\nu\sigma}R_{\nu\sigma} = \Lambda g^{\nu\sigma}g_{\nu\sigma} = 4\Lambda,$$

a constant, written $R = 4\Lambda$.

The Einstein tensor follows without any further computation:

$$G_{\mu\nu} = R_{\mu\nu} - \tfrac{1}{2}Rg_{\mu\nu} = \Lambda g_{\mu\nu} - 2\Lambda g_{\mu\nu} = -\Lambda g_{\mu\nu},$$

so every Einstein component is minus the Ricci component in the same slot and the same index position, and the mixed Einstein tensor is $G^\mu{}_\nu = -\Lambda\delta^\mu{}_\nu$.
Written out, $G_{\mu\nu} + \Lambda g_{\mu\nu} = 0$, which is the field equation of Step 2 and the sense in which this spacetime is empty.

---

## Step 11. The Weyl tensor, computed rather than assumed

The Weyl tensor in four dimensions is Riemann with its traces removed,

$$C_{\mu\nu\rho\sigma} = R_{\mu\nu\rho\sigma} - \frac{1}{2}\left(g_{\mu\rho}R_{\sigma\nu} - g_{\mu\sigma}R_{\rho\nu} - g_{\nu\rho}R_{\sigma\mu} + g_{\nu\sigma}R_{\rho\mu}\right) + \frac{R}{6}\left(g_{\mu\rho}g_{\sigma\nu} - g_{\mu\sigma}g_{\rho\nu}\right),$$

which is the expression the checker computes, the one derived in `_tools/derivations/weyl.md`.
It is not Riemann here, because the Ricci tensor is not zero, and it is not the Ricci tensor either.
It has to be computed.

Write $P_{\mu\nu\rho\sigma} = g_{\mu\rho}g_{\nu\sigma} - g_{\mu\sigma}g_{\nu\rho}$ for the tensor that appears in all three terms.
From Step 8, $R_{\mu\nu\rho\sigma} = \tfrac{\Lambda}{3}P_{\mu\nu\rho\sigma}$.
From Step 9, $R_{\mu\nu} = \Lambda g_{\mu\nu}$, so the middle bracket is

$$g_{\mu\rho}\Lambda g_{\sigma\nu} - g_{\mu\sigma}\Lambda g_{\rho\nu} - g_{\nu\rho}\Lambda g_{\sigma\mu} + g_{\nu\sigma}\Lambda g_{\rho\mu} = 2\Lambda\left(g_{\mu\rho}g_{\nu\sigma} - g_{\mu\sigma}g_{\nu\rho}\right) = 2\Lambda P_{\mu\nu\rho\sigma},$$

where the four terms have collapsed to two because $g$ is symmetric.
From Step 10, $R = 4\Lambda$, so the last term is $\tfrac{4\Lambda}{6}P_{\mu\nu\rho\sigma}$.
Adding the three,

$$C_{\mu\nu\rho\sigma} = \frac{\Lambda}{3}P_{\mu\nu\rho\sigma} - \frac{1}{2}\cdot 2\Lambda P_{\mu\nu\rho\sigma} + \frac{2\Lambda}{3}P_{\mu\nu\rho\sigma} = \Lambda\left(\frac{1}{3} - 1 + \frac{2}{3}\right)P_{\mu\nu\rho\sigma} = 0.$$

The three rational numbers cancel identically, with no appeal to the value of $\Lambda$, to the coordinates or to the chart.
So the Weyl tensor vanishes in every slot, with the first index up and with every index lowered, and the checker requires every component of both to vanish, as it does for FRW for the same reason.

This is the arithmetic behind the word conformally flat.
A maximally symmetric spacetime has no tidal curvature at all: the whole of its Riemann tensor is trace, so removing the traces leaves nothing.
A vanishing Weyl tensor here is not the vanishing Ricci tensor of a vacuum: Schwarzschild has $R_{\mu\nu} = 0$ and $C = R \neq 0$, de Sitter has $C = 0$ and $R_{\mu\nu} \neq 0$, and the two are opposite corners of the same decomposition.

---

## Step 12. The Kretschmann scalar

The Kretschmann scalar is $K = R_{\mu\nu\rho\sigma}R^{\mu\nu\rho\sigma}$, and with Riemann proportional to $P$ it needs only the square of $P$.
Raising all four indices and contracting,

$$P_{\mu\nu\rho\sigma}P^{\mu\nu\rho\sigma} = \left(g_{\mu\rho}g_{\nu\sigma} - g_{\mu\sigma}g_{\nu\rho}\right)\left(g^{\mu\rho}g^{\nu\sigma} - g^{\mu\sigma}g^{\nu\rho}\right).$$

The first term against the first gives $\delta^\mu{}_\mu\delta^\nu{}_\nu = 4\cdot 4 = 16$.
The first against the second gives $\delta^\mu{}_\sigma\delta^\sigma{}_\mu = 4$, with a minus sign, and the second against the first the same.
The second against the second gives $16$.
So

$$P_{\mu\nu\rho\sigma}P^{\mu\nu\rho\sigma} = 16 - 4 - 4 + 16 = 24,$$

which is $2n(n-1)$ at $n = 4$, and

$$K = \left(\frac{\Lambda}{3}\right)^2\cdot 24 = \frac{24\Lambda^2}{9} = \frac{8\Lambda^2}{3}.$$

That value is a constant: the same number at the origin, at the horizon and at every point of the flat slicing, and the other chart gives it again as $24H^4/c^4$ in Step 16.
Nothing in this spacetime is singular, and $K$ is the shortest proof of it.

---

## Step 13. The horizon at $r = \sqrt{3/\Lambda}$ is a coordinate fact

At $r = \ell = \sqrt{3/\Lambda}$ the function $f$ vanishes, so $g_{tt} \to 0$ and $g_{rr} \to \infty$.
Three quantities say at once that this is the chart failing and not the geometry.

$K = 8\Lambda^2/3$ is finite there, and so is every other curvature invariant, because all of them are built from $\Lambda$ and the metric alone.

$R^\mu{}_\nu = \Lambda\delta^\mu{}_\nu$ has no $r$ in it at all, so the Ricci tensor is as regular at $r = \ell$ as it is anywhere.

The Riemann tensor with its first index up has components such as $R^t{}_{rtr} = \Lambda/(3-\Lambda r^2)$ that do blow up, and components such as $R^t{}_{\theta t\theta} = \Lambda r^2/3$ that do not, in the same tensor at the same point.
A tensor cannot be singular in one slot and regular in another at the same point in a good chart, so the chart is bad.

What the surface is instead is a Killing horizon of $\partial_t$, the boundary of the region a single observer at $r = 0$ can ever receive a signal from.
Its surface gravity is

$$\kappa = \frac{c^2}{2}\left|f'(\ell)\right| = \frac{c^2}{2}\cdot\frac{2\Lambda\ell}{3} = \frac{c^2\Lambda}{3}\sqrt{\frac{3}{\Lambda}} = c^2\sqrt{\frac{\Lambda}{3}} = cH,$$

the one combination the geometry can make, and it is what fixes the Gibbons-Hawking temperature $T = \hbar H/2\pi k_B$ that every de Sitter observer sees.
Unlike a black hole horizon it is observer dependent: it is centred wherever the observer at $r = 0$ is placed, and the flat slicing of Step 15 crosses it without noticing.
The domain of the static patch states the radius, alongside the range of $r$, rather than leaving a reader to find the zero of $g_{tt}$.

---

## Step 14. The geodesic equations of the static patch

The geodesic equation is

$$\ddot{x}^\mu + \Gamma^\mu{}_{\nu\rho}\dot{x}^\nu\dot{x}^\rho = 0,$$

with the $\Gamma$ of Step 5 and the dots velocities of the chart, $\dot{t} = d(ct)/d\lambda$.

For $\mu = t$ the only symbol is $\Gamma^t{}_{tr}$, which appears twice in the sum:

$$\ddot{t} - \frac{2\Lambda r}{3 - \Lambda r^2}\dot{t}\dot{r} = 0.$$

For $\mu = r$ the four diagonal symbols contribute:

$$\ddot{r} - \frac{\Lambda r\left(3 - \Lambda r^2\right)}{9}\dot{t}^2 + \frac{\Lambda r}{3 - \Lambda r^2}\dot{r}^2 - \frac{r\left(3 - \Lambda r^2\right)}{3}\left(\dot\theta^2 + \sin^2\theta\,\dot\phi^2\right) = 0,$$

where the last bracket collects $\Gamma^r{}_{\theta\theta}$ and $\Gamma^r{}_{\phi\phi}$, which differ only by the $\sin^2\theta$.

The two angular equations are the ones any spherically symmetric metric carries:

$$\ddot\theta + \frac{2}{r}\dot{r}\dot\theta - \cos\theta\sin\theta\,\dot\phi^2 = 0, \qquad \ddot\phi + \frac{2}{r}\dot{r}\dot\phi + 2\cot\theta\,\dot\theta\dot\phi = 0.$$

The first equation integrates once, because $\partial_t$ is a Killing vector: multiplying it by $f$ gives $\tfrac{d}{d\lambda}\left(f\dot{t}\right) = 0$, so

$$E = \left(1 - \frac{\Lambda r^2}{3}\right)\dot{t}$$

is conserved, and it is the energy the timelike Killing vector $\partial_t$ defines, per unit rest mass and in units of $c^2$.
The fourth equation integrates the same way on $\partial_\phi$, giving the angular momentum $L = r^2\sin^2\theta\,\dot\phi$.
Neither constant appears among the geodesic equations, which are equations rather than their solutions, but both are one line from them and worth knowing.

---

## Step 15. The flat slicing, and the change of coordinates that produces it

The second chart is the one cosmology uses,

$$ds^2 = -c^2dt^2 + e^{2Ht}\left(dx^2 + dy^2 + dz^2\right), \qquad H = c\sqrt{\frac{\Lambda}{3}} = \frac{c}{\ell}.$$

Its parameter is $H$ rather than $\Lambda$, carrying an inverse time, and the two are the same thing through

$$\frac{3H^2}{c^2} = \Lambda.$$

That relation is what turns every component of this chart into a statement about $\Lambda$: $R = 12H^2/c^2$ is $4\Lambda$, and $K = 24H^4/c^4$ is $8\Lambda^2/3$.

The change of coordinates from the flat slicing to the static patch is worth doing in full, because it produces the static $f$ out of nothing but a completed square.
Write the spatial part in polar form, $\rho^2 = x^2+y^2+z^2$, and set

$$r = \rho\,e^{Ht}.$$

Then $dr = e^{Ht}d\rho + Hr\,dt$, so $e^{Ht}d\rho = dr - Hr\,dt$ and

$$e^{2Ht}d\rho^2 = dr^2 - 2Hr\,dt\,dr + H^2r^2dt^2.$$

Substituting into the line element,

$$ds^2 = -c^2\left(1 - \frac{H^2r^2}{c^2}\right)dt^2 - 2Hr\,dt\,dr + dr^2 + r^2d\Omega^2,$$

and the bracket is already $f$, since $H^2/c^2 = \Lambda/3$.
Complete the square on the first two terms:

$$-c^2f\,dt^2 - 2Hr\,dt\,dr = -c^2f\left(dt + \frac{Hr\,dr}{c^2f}\right)^2 + \frac{H^2r^2}{c^2f}dr^2.$$

The leftover joins $dr^2$ to give $\left(H^2r^2/c^2f + 1\right)dr^2 = f^{-1}dr^2$, because $H^2r^2/c^2 + f = 1$.
So with the static time

$$T = t + \int\frac{Hr\,dr}{c^2f} = t - \frac{1}{2H}\ln\left(1 - \frac{\Lambda r^2}{3}\right)$$

the line element becomes $-fc^2dT^2 + f^{-1}dr^2 + r^2d\Omega^2$, which is Step 3.

Two things follow from the shape of that transformation.
The logarithm is real only for $r < \ell$, so the static patch is the part of the flat slicing inside the horizon, and the flat slicing carries on past it with nothing happening to any of its components.
And $T \to \infty$ as $r \to \ell$ at fixed $t$, which is the usual way a Killing horizon looks from a chart that crosses it: infinitely far away in the static time, a finite step in the slicing time.

---

## Step 16. The flat slicing's connection and curvature

In the chart $x^0 = ct$ the metric of the flat slicing is

$$g_{\mu\nu} = \mathrm{diag}\left(-1,\, a^2,\, a^2,\, a^2\right), \qquad a = e^{Ht},$$

with $g_{tt} = -1$ exactly, since the whole of the $c^2$ went into $(c\,dt)^2$.
Derivatives are taken with respect to the chart time, so write

$$a' = \frac{da}{d(ct)} = \frac{H}{c}a, \qquad \frac{a'}{a} = \frac{H}{c},$$

an inverse length, which is where every factor of $c$ in this chart comes from and the only place one appears.

The connection of a spatially flat expanding metric has two families:

$$\Gamma^t{}_{ij} = a a'\delta_{ij} = \frac{H}{c}e^{2Ht}\delta_{ij}, \qquad \Gamma^i{}_{tj} = \Gamma^i{}_{jt} = \frac{a'}{a}\delta^i{}_j = \frac{H}{c}\delta^i{}_j,$$

with $i$ and $j$ running over $x$, $y$ and $z$.
That is nine nonzero symbols, three of the first kind and six of the second, and it is all of the Christoffel symbols with an upper index.
Lowering the first index gives $\Gamma_{tij} = -aa'\delta_{ij}$ and $\Gamma_{itj} = \Gamma_{ijt} = aa'\delta_{ij}$, the only change being the sign that $g_{tt} = -1$ puts on the first family.

The Riemann tensor takes two computations.
For the time and space planes,

$$R^t{}_{itj} = \partial_{ct}\Gamma^t{}_{ij} - \Gamma^t{}_{jk}\Gamma^k{}_{it} = \left((a')^2 + aa''\right)\delta_{ij} - aa'\cdot\frac{a'}{a}\delta_{ij} = aa''\delta_{ij},$$

and with $a = e^{Ht}$ the chart second derivative is $a'' = (H/c)^2a$, so

$$R^t{}_{itj} = \frac{H^2}{c^2}e^{2Ht}\delta_{ij}.$$

For the purely spatial planes, with $i \neq j$, both derivative terms vanish because $\Gamma^i{}_{jj}$ and $\Gamma^i{}_{ji}$ are zero when the indices differ, and only one product survives:

$$R^i{}_{jij} = \Gamma^i{}_{it}\Gamma^t{}_{jj} = \frac{a'}{a}\cdot aa' = \left(a'\right)^2 = \frac{H^2}{c^2}e^{2Ht},$$

the expansion of the slicing making the spatial sections curve in spacetime even though each of them is flat.

Both results are the maximally symmetric form of Step 8 with $\Lambda$ replaced by $3H^2/c^2$:

$$R^\mu{}_{\nu\rho\sigma} = \frac{H^2}{c^2}\left(\delta^\mu{}_\rho g_{\nu\sigma} - \delta^\mu{}_\sigma g_{\nu\rho}\right),$$

since $\tfrac{H^2}{c^2}g_{ij} = \tfrac{H^2}{c^2}a^2\delta_{ij}$ and $-\tfrac{H^2}{c^2}g_{tt} = \tfrac{H^2}{c^2}$, which is the component $R^x{}_{ttx}$.
So the same $\Lambda$ appears in a chart that never mentions it, and everything from Step 9 to Step 12 carries over verbatim:

$$R_{\mu\nu} = \frac{3H^2}{c^2}g_{\mu\nu}, \qquad R^\mu{}_\nu = \frac{3H^2}{c^2}\delta^\mu{}_\nu, \qquad R = \frac{12H^2}{c^2}, \qquad G_{\mu\nu} = -\frac{3H^2}{c^2}g_{\mu\nu}, \qquad C_{\mu\nu\rho\sigma} = 0,$$

$$K = \left(\frac{H^2}{c^2}\right)^2\cdot 24 = \frac{24H^4}{c^4} = \frac{8\Lambda^2}{3}.$$

The components carry the scale factor where the metric does: $R_{xx} = 3H^2e^{2Ht}/c^2$ has it, $R^x{}_x = 3H^2/c^2$ does not, and $R^{xx} = 3H^2e^{-2Ht}/c^2$ has its inverse.

---

## Step 17. Both charts sit on one hyperboloid, and the flat slicing covers half of it

de Sitter space is the hyperboloid

$$-\left(X^0\right)^2 + \left(X^1\right)^2 + \left(X^2\right)^2 + \left(X^3\right)^2 + \left(X^4\right)^2 = \ell^2$$

in five dimensional Minkowski space with metric $\mathrm{diag}(-1,1,1,1,1)$, and both charts are charts on it.
No component needs the hyperboloid, but it is what makes the relation between the two charts a fact about regions rather than about formulas.

The static patch is

$$X^0 = \sqrt{\ell^2 - r^2}\,\sinh\frac{cT}{\ell}, \qquad X^4 = \sqrt{\ell^2 - r^2}\,\cosh\frac{cT}{\ell}, \qquad X^i = r\,n^i(\theta,\phi),$$

with $n^i$ the unit radial vector.
The induced metric is $-fc^2dT^2 + f^{-1}dr^2 + r^2d\Omega^2$ exactly, and the constraint holds identically.
Since $\left(X^4\right)^2 - \left(X^0\right)^2 = \ell^2 - r^2 > 0$ with $X^4 > 0$, the chart covers the wedge $X^4 > \left|X^0\right|$, one of four such wedges, which is the sense in which a static observer sees a quarter of the spacetime.

The flat slicing is

$$X^0 = \ell\sinh Ht + \frac{\rho^2}{2\ell}e^{Ht}, \qquad X^4 = \ell\cosh Ht - \frac{\rho^2}{2\ell}e^{Ht}, \qquad X^i = e^{Ht}x^i,$$

whose induced metric is $-c^2dt^2 + e^{2Ht}\left(dx^2+dy^2+dz^2\right)$ and which satisfies the constraint identically as well.
Adding the first two,

$$X^0 + X^4 = \ell\left(\sinh Ht + \cosh Ht\right) = \ell e^{Ht} > 0,$$

so the flat slicing covers exactly the half of the hyperboloid on which $X^0 + X^4 > 0$, and no point outside it, for any $t$ and any $x^i$.
That half is the expanding one, bounded by the past horizon $X^0 + X^4 = 0$, and the static patch is the part of it with $r < \ell$.
The half covering is easy to miss, because every coordinate of the flat slicing runs over the whole real line, which by itself suggests the whole spacetime.

---

## Step 18. The geodesic equations of the flat slicing

With the connection of Step 16 the four equations are

$$\ddot{t} + \frac{H}{c}e^{2Ht}\left(\dot{x}^2 + \dot{y}^2 + \dot{z}^2\right) = 0,$$

$$\ddot{x} + \frac{2H}{c}\dot{t}\dot{x} = 0, \qquad \ddot{y} + \frac{2H}{c}\dot{t}\dot{y} = 0, \qquad \ddot{z} + \frac{2H}{c}\dot{t}\dot{z} = 0,$$

the factor of two in the last three coming from $\Gamma^i{}_{tj}$ appearing in both orderings.
A particle at rest in these coordinates has $\dot{x} = \dot{y} = \dot{z} = 0$, and then the first equation gives $\ddot{t} = 0$ and the other three are satisfied: the comoving worldlines are geodesics, which is what makes $x$, $y$ and $z$ comoving labels.
Any peculiar velocity decays, since the spatial equations integrate exactly to $\dot{x}^i \propto e^{-2Ht}$, so the proper peculiar momentum $a\dot{x}^i$ falls as $1/a$, which is the redshifting of momentum in an expanding universe read straight off the connection.

---

## Step 19. Every expression is dimensionally consistent

The chart coordinates of the static patch are $(ct, r, \theta, \phi)$, of dimensions $L$, $L$, $1$, $1$, and its parameter is $[\Lambda] = L^{-2}$.
Those of the flat slicing are $(ct, x, y, z)$, all $L$, with $[H] = T^{-1}$.
A metric component then carries

$$[g_{\mu\nu}] = \frac{L^2}{[x^\mu][x^\nu]},$$

a Christoffel symbol carries $L^{-1}$ corrected by $[x^\mu]/L$ per upper index and $L/[x^\mu]$ per lower one, a Riemann or Weyl component carries $L^{-2}$ corrected the same way, and $K$ carries $L^{-4}$.

Four samples.

$g_{tt} = -\left(1 - \Lambda r^2/3\right)$ must be dimensionless, and $\Lambda r^2$ is $L^{-2}\cdot L^2 = 1$, so both terms are.
This is the check that decides that $\Lambda$ is a curvature: no other dimension lets it stand beside the $1$.

$\Gamma^r{}_{tt} = -\Lambda r(3-\Lambda r^2)/9$ carries $L^{-1}$ from the field and $L/L = 1$ from each of its three indices, so it must carry $L^{-1}$, and $\Lambda r$ does.

$R = 4\Lambda$ carries $L^{-2}$ and $K = 8\Lambda^2/3$ carries $L^{-4}$, as a scalar of each rank must.

In the flat slicing every factor of $c$ earns its place in the same way.
$\Gamma^x{}_{tx} = H/c$ must carry $L^{-1}$, and $H/c$ is $T^{-1}\cdot T/L = L^{-1}$.
$R_{xx} = 3H^2e^{2Ht}/c^2$ must carry $L^{-2}$, which $H^2/c^2$ does, and the exponential is dimensionless because $Ht$ is.
Had the components been written in the bare chart rather than in $x^0 = ct$, those factors of $c$ would not be there and the dimensional pass would name every one of them.

The geodesic equations are measured against their own second derivatives.
In the first equation of Step 18, $\ddot{t}$ carries $L/\lambda^2$, and the term $\left(H/c\right)e^{2Ht}\dot{x}^2$ carries

$$\frac{1}{L}\cdot\left(\frac{L}{\lambda}\right)^2 = \frac{L}{\lambda^2},$$

which matches, and matches only because $\dot{t}$ is $d(ct)/d\lambda$ and $\dot{x}$ is $dx/d\lambda$, both velocities of the same chart.
On the other reading, with $\dot t = dt/d\lambda$, that term would be one factor of $c$ away from its left hand side.
It is not, on any of the 218 expressions of the two systems between them.

---

## Step 20. The expressions, and the declarations the checker needs

The static patch has 113 expressions: the line element, four metric components and four inverse ones, thirteen Christoffel symbols in each of two index positions, twenty four Riemann components in each of two, four Ricci components in each of three, four Einstein components in each of three, two scalars, a vanishing Weyl tensor in each of two index positions, and four geodesic equations.
The flat slicing has 105, differing only in having nine Christoffel symbols per index position rather than thirteen.

For the checker to read either of them, `DIMENSIONS` in `verify_metrics.py` needs one line per system:

    ("de_sitter", "static_spherical"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "\\Lambda": "1/L**2",
    },
    ("de_sitter", "flat_slicing"): {
        "t": "T", "x": "L", "y": "L", "z": "L", "H": "1/T",
    },

The declaration `"t": "T"` is what tells the checker that $t$ is a time and so that the chart coordinate is $ct$; every factor of $c$ in the comparison follows from it.
Declaring $\Lambda$ as $L^{-2}$ and $H$ as $T^{-1}$ is what makes the two charts balance without either of them naming the other's parameter.
Neither system needs a line in `PARAMETER_RELATIONS`, because neither parameter is constrained: every positive $\Lambda$ gives a de Sitter space, and the components hold for all of them.

Running

    /tmp/mfs-venv/bin/python _tools/derivations/verify_metrics.py --system de_sitter/static_spherical --system de_sitter/flat_slicing

reports no disagreements and no dimensional failures.
