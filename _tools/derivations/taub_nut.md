# The Taub-NUT vacuum spacetime

This is the working behind the coordinate system in `MFS/assets/data/metrics/taub_nut.json`.
Every number the entry prints is derived here, in order, from the line element down to the geodesic equations.
Nothing is left as an exercise and nothing is asserted that is not computed.

The companion script `_tools/derivations/verify_metrics.py` does the same work in sympy and compares it against the published file, so the algebra below is checkable by hand and by machine independently.

Taub-NUT is the first entry in the collection whose metric is not diagonal.
The whole of what makes it interesting sits in one off diagonal component, $g_{t\phi}$, and the entry publishes it rather than rotating it away, because it cannot be rotated away.
Step 2 says what that term is, Step 8 shows it does not disturb the vacuum condition, Step 11 shows it doubling the Kretschmann scalar into a real part and an imaginary part, and Step 14 says where the price is paid.

---

## Step 1. Conventions

The signature is $(-,+,+,+)$ and the Christoffel symbols are those of the Levi-Civita connection,

$$\Gamma^\mu_{\nu\rho} = \tfrac{1}{2}g^{\mu\alpha}\left(\partial_\nu g_{\alpha\rho} + \partial_\rho g_{\alpha\nu} - \partial_\alpha g_{\nu\rho}\right),$$

with the lowered symbol $\Gamma_{\mu\nu\rho} = g_{\mu\alpha}\Gamma^\alpha_{\nu\rho}$.
The Riemann tensor is

$$R^\mu{}_{\nu\rho\sigma} = \partial_\rho \Gamma^\mu_{\nu\sigma} - \partial_\sigma \Gamma^\mu_{\nu\rho} + \Gamma^\mu_{\rho\lambda}\Gamma^\lambda_{\nu\sigma} - \Gamma^\mu_{\sigma\lambda}\Gamma^\lambda_{\nu\rho},$$

which is what the published Riemann components are in.

The Ricci tensor is the standard contraction,

$$R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu},$$

settled for the collection on 2026-09-18.
For this spacetime the settlement costs nothing, because the contraction vanishes and so does its negative: Step 8 computes it and gets zero in every slot, so an entry written with either convention would print the same empty block.
The Weyl tensor is built from the same contraction in Step 10, which is the only place the choice could have shown.

Factors of $G$ and $c$ are kept explicit.
Both parameters of the solution are lengths.
The mass enters only through the geometric mass

$$m = \frac{GM}{c^2},$$

which is half the Schwarzschild radius, and the NUT parameter $l$ is a length beside it.
That is the same bargain Schwarzschild strikes with $r_s = 2GM/c^2$, and it is why no published component of this entry carries a factor of $c$ or of $G$ at all.

The chart is the collection's, the one whose zeroth coordinate is

$$x^0 = ct.$$

The index is printed with the bare letter $t$, but the component printed against it is a component of that chart.
Everything from Step 3 onward is computed directly in it, so no conversion is needed at the end.
Because the metric of this spacetime does not depend on $t$, the rescaling leaves nothing behind: in the chart $x^0 = ct$ every component of the metric, of the connection and of every curvature tensor is a function of $r$ and $\theta$ alone, with no $c$ anywhere.
The dots in the geodesic equations are velocities of that same chart, so $\dot{t}$ means $d(ct)/d\lambda$, and Step 15 checks that this is the reading on which every published term balances.

Three abbreviations run through everything below.

$$\rho^2 = r^2 + l^2, \qquad \Delta = r^2 - 2mr - l^2, \qquad f = \frac{\Delta}{\rho^2}.$$

They are shorthand for this document only.
The published file writes every component out in $r$, $\theta$, $m$ and $l$, because the checker reads the file symbol by symbol and has no way to be told what an abbreviation means.

---

## Step 2. The line element, and the one form that is the point of it

The solution is

$$ds^2 = -f(r)\left(c\,dt + 2l\cos\theta\,d\phi\right)^2 + \frac{dr^2}{f(r)} + \left(r^2+l^2\right)\left(d\theta^2 + \sin^2\theta\,d\phi^2\right),$$

with

$$f(r) = \frac{r^2 - 2mr - l^2}{r^2 + l^2}.$$

Setting $l = 0$ gives $f = 1 - 2m/r = 1 - r_s/r$, kills the cross term and returns Schwarzschild exactly, which is the first check to make on any expression below.

Write the one form inside the bracket as

$$\sigma = c\,dt + A, \qquad A = 2l\cos\theta\,d\phi.$$

$A$ is the vector potential of a magnetic monopole of strength $2l$, written in the gauge that is symmetric between the two poles, and

$$dA = -2l\sin\theta\,d\theta\wedge d\phi$$

is $-2l$ times the area form of the unit sphere, the monopole field itself.
That is the content of calling $l$ a gravitomagnetic charge: the NUT parameter enters the metric the way a magnetic charge enters electromagnetism, through a potential whose curl is a monopole field, and a monopole potential cannot be made regular everywhere on a sphere.
Step 14 collects what that costs.

The horizon is where $f$ vanishes, that is where $\Delta = 0$:

$$r_\pm = m \pm \sqrt{m^2 + l^2}.$$

Both roots are real for every $m$ and every $l$, and $r_- < 0 < r_+$, so unlike Reissner-Nordstrom there is no extremal case and no naked case.
The entry publishes the exterior, $r > r_+$.

---

## Step 3. The metric matrix, its determinant and its inverse

Expanding the square in Step 2 with $x^0 = ct$ and ordering the coordinates $(t, r, \theta, \phi)$,

$$g_{\mu\nu} = \begin{pmatrix} -f & 0 & 0 & -2lf\cos\theta \\ 0 & \dfrac{1}{f} & 0 & 0 \\ 0 & 0 & \rho^2 & 0 \\ -2lf\cos\theta & 0 & 0 & \rho^2\sin^2\theta - 4l^2f\cos^2\theta\end{pmatrix}.$$

Every entry is dimensionless, as it must be in a chart whose four coordinates are a length, an angle and an angle, once the zeroth is $ct$.
The two mixed entries are the whole of the difference from Schwarzschild.

The determinant comes out of the $2\times 2$ block in $(t,\phi)$ and the two diagonal entries between.
For the block,

$$g_{tt}g_{\phi\phi} - g_{t\phi}^2 = -f\left(\rho^2\sin^2\theta - 4l^2f\cos^2\theta\right) - 4l^2f^2\cos^2\theta = -f\rho^2\sin^2\theta,$$

where the two terms in $4l^2f^2\cos^2\theta$ cancel exactly.
Multiplying by $g_{rr}g_{\theta\theta} = \rho^2/f$,

$$\det g = -\rho^4\sin^2\theta, \qquad \sqrt{-g} = \left(r^2+l^2\right)\sin\theta.$$

The volume element knows nothing about $m$ and nothing about the cross term.
It is the flat measure of a sphere of areal radius $\rho$, which is the first hint that the NUT parameter is not a source of anything.

Inverting the same block,

$$g^{tt} = \frac{g_{\phi\phi}}{-f\rho^2\sin^2\theta} = -\frac{1}{f} + \frac{4l^2\cos^2\theta}{\rho^2\sin^2\theta}, \qquad g^{t\phi} = \frac{-g_{t\phi}}{-f\rho^2\sin^2\theta} = -\frac{2l\cos\theta}{\rho^2\sin^2\theta}, \qquad g^{\phi\phi} = \frac{1}{\rho^2\sin^2\theta},$$

and the rest are reciprocals:

$$g^{rr} = f, \qquad g^{\theta\theta} = \frac{1}{\rho^2}.$$

Two things are worth reading off before going on.
The inverse metric carries $1/\sin^2\theta$ in three of its entries, which is the axis making its first appearance, and $g^{\phi\phi}$ is the same as it would be on a sphere, so the cross term has been pushed entirely into $g^{tt}$ and $g^{t\phi}$.

---

## Step 4. The Christoffel symbols with an upper index

Nothing depends on $t$ or on $\phi$, so only $\partial_r$ and $\partial_\theta$ ever act.
Writing $f' = df/dr$ and

$$P = mr^2 + 2l^2r - ml^2, \qquad f' = \frac{2P}{\rho^4},$$

which is the derivative of $\Delta/\rho^2$ with the quotient rule and nothing else, the twenty five nonzero symbols are the following fifteen and their reflections in the lower pair.

$$\Gamma^t{}_{tr} = \frac{f'}{2f}, \qquad \Gamma^t{}_{t\theta} = -\frac{2l^2f}{\rho^2}\cot\theta, \qquad \Gamma^t{}_{r\phi} = -\frac{2W\cos\theta}{\rho^2\Delta}, $$

$$\Gamma^t{}_{\theta\phi} = \frac{l\left[\left(\rho^4+4l^2\Delta\right)\sin^2\theta - 2\left(\rho^4+2l^2\Delta\right)\right]}{\rho^4\sin\theta},$$

$$\Gamma^r{}_{tt} = \frac{ff'}{2}, \qquad \Gamma^r{}_{t\phi} = lff'\cos\theta, \qquad \Gamma^r{}_{rr} = -\frac{f'}{2f}, \qquad \Gamma^r{}_{\theta\theta} = -rf,$$

$$\Gamma^r{}_{\phi\phi} = -rf\sin^2\theta + 2l^2ff'\cos^2\theta,$$

$$\Gamma^\theta{}_{t\phi} = -\frac{lf\sin\theta}{\rho^2}, \qquad \Gamma^\theta{}_{r\theta} = \frac{r}{\rho^2}, \qquad \Gamma^\theta{}_{\phi\phi} = -\left(1 + \frac{4l^2f}{\rho^2}\right)\sin\theta\cos\theta,$$

$$\Gamma^\phi{}_{t\theta} = \frac{lf}{\rho^2\sin\theta}, \qquad \Gamma^\phi{}_{r\phi} = \frac{r}{\rho^2}, \qquad \Gamma^\phi{}_{\theta\phi} = \left(1 + \frac{2l^2f}{\rho^2}\right)\cot\theta.$$

Here $W$ is the polynomial that Step 7 will meet again,

$$W = l\left(r^3 - 3mr^2 - 3l^2r + l^2m\right).$$

Set $l = 0$ and the list collapses to Schwarzschild's nine: $\Gamma^t{}_{tr} = f'/2f$, $\Gamma^r{}_{tt} = ff'/2$, $\Gamma^r{}_{rr} = -f'/2f$, $\Gamma^r{}_{\theta\theta} = -rf$, $\Gamma^r{}_{\phi\phi} = -rf\sin^2\theta$, $\Gamma^\theta{}_{r\theta} = \Gamma^\phi{}_{r\phi} = 1/r$, $\Gamma^\theta{}_{\phi\phi} = -\sin\theta\cos\theta$ and $\Gamma^\phi{}_{\theta\phi} = \cot\theta$.
Six symbols survive only because $l \neq 0$, and every one of them carries an odd power of $l$: those are the ones the cross term paid for.

Two of them deserve a sentence.
$\Gamma^\theta{}_{t\phi}$ and $\Gamma^\phi{}_{t\theta}$ are the gravitomagnetic terms proper.
They couple a velocity in time to a velocity in angle, which is exactly how a magnetic field couples in the Lorentz force, and they are why a particle dropped radially in this spacetime does not fall radially.

The published file writes each of these twenty five with $\rho^2$, $\Delta$, $f$ and $P$ substituted out, so for instance the first of them appears as

$$\Gamma^t{}_{tr} = \frac{mr^2 + 2l^2r - l^2m}{\left(r^2+l^2\right)\left(r^2 - 2mr - l^2\right)}.$$

---

## Step 5. The lowered Christoffel symbols

$\Gamma_{\mu\nu\rho} = g_{\mu\alpha}\Gamma^\alpha{}_{\nu\rho}$, and because the metric mixes $t$ with $\phi$ the lowering is not a single multiplication in those two slots:

$$\Gamma_{t\nu\rho} = g_{tt}\Gamma^t{}_{\nu\rho} + g_{t\phi}\Gamma^\phi{}_{\nu\rho}, \qquad \Gamma_{\phi\nu\rho} = g_{\phi t}\Gamma^t{}_{\nu\rho} + g_{\phi\phi}\Gamma^\phi{}_{\nu\rho},$$

while the $r$ and $\theta$ rows lower with one factor each.
The count is unchanged at twenty five, since the mixing is invertible and no symbol is annihilated by it.
The entry publishes both variants in full.

---

## Step 6. The orthonormal coframe

The Riemann tensor of this metric has eighty nonzero components with one index up and seventy two with all four down, and writing them out in the coordinate basis explains nothing.
In an orthonormal coframe there are eight, and they are two numbers.

Take

$$e^0 = \sqrt{f}\left(dx^0 + A\right), \qquad e^1 = \frac{dr}{\sqrt{f}}, \qquad e^2 = \rho\,d\theta, \qquad e^3 = \rho\sin\theta\,d\phi,$$

which is the line element of Step 2 read straight off, so that

$$g = -e^0\otimes e^0 + e^1\otimes e^1 + e^2\otimes e^2 + e^3\otimes e^3$$

and the frame metric is $\eta_{ab} = \mathrm{diag}(-1,1,1,1)$ exactly.
This is a frame and not a chart: $e^0$ is not the differential of anything, and that is the same fact as the Misner string.

---

## Step 7. The Riemann tensor in the coframe

Transforming the seventy two lowered components into this frame leaves eight, and they are built from two polynomials:

$$U = mr^3 + 3l^2r^2 - 3l^2mr - l^4, \qquad W = lr^3 - 3lmr^2 - 3l^3r + l^3m.$$

The eight are

$$R_{0101} = -\frac{2U}{\rho^6}, \qquad R_{0202} = R_{0303} = \frac{U}{\rho^6}, \qquad R_{1212} = R_{1313} = -\frac{U}{\rho^6}, \qquad R_{2323} = \frac{2U}{\rho^6},$$

$$R_{0123} = -\frac{2W}{\rho^6}, \qquad R_{0213} = -\frac{W}{\rho^6}, \qquad R_{0312} = \frac{W}{\rho^6},$$

and every frame component not obtainable from these by the symmetries of Riemann is zero, $R_{1323}$ among them.

Now split the curvature the way an observer along $e^0$ sees it.
The electric part is the tidal tensor $E_{ij} = R_{0i0j}$ and the magnetic part is its dual, $B_{ij} = \tfrac{1}{2}\varepsilon_i{}^{kl}R_{0jkl}$.
From the eight components above,

$$E_{ij} = \frac{U}{\rho^6}\,\mathrm{diag}(-2,1,1), \qquad B_{ij} = \frac{W}{\rho^6}\,\mathrm{diag}(-2,1,1).$$

Both are the same matrix, $\mathrm{diag}(-2,1,1)$, scaled by one number each, which is the signature of a Petrov type D field: one radial direction stretched twice as hard as the two transverse ones are squeezed.
Schwarzschild has the electric half of this and nothing else.
The NUT parameter turns on the magnetic half, and it is genuinely a second, independent tidal field: no observer can transform $B$ away, because both halves are built from the same Weyl tensor.

The two numbers are one complex number.
Since

$$U + iW = (m + il)(r - il)^3 \qquad\text{and}\qquad \rho^6 = \left[(r+il)(r-il)\right]^3,$$

the ratio is

$$\frac{U + iW}{\rho^6} = \frac{m + il}{(r + il)^3},$$

so the entire curvature of Taub-NUT is the Schwarzschild expression $m/r^3$ with the mass complexified to $m + il$ and the radius complexified to $r + il$.
That identity is not decoration.
It is the reason $l$ is called a gravitomagnetic charge, it is why every polynomial in this document has coefficients $1, 3, 3, 1$ up to signs, and Step 11 gets the Kretschmann scalar out of it in one line.

---

## Step 8. The Ricci tensor vanishes

In an orthonormal frame the Ricci tensor is $R_{ab} = \eta^{cd}R_{cadb} = -R_{0a0b} + \sum_k R_{kakb}$.

The diagonal entries are four sums, and all four use only the tracelessness of $\mathrm{diag}(-2,1,1)$:

$$R_{00} = R_{0101} + R_{0202} + R_{0303} = \frac{U}{\rho^6}\left(-2 + 1 + 1\right) = 0,$$

$$R_{11} = -R_{0101} + R_{1212} + R_{1313} = \frac{U}{\rho^6}\left(2 - 1 - 1\right) = 0,$$

$$R_{22} = -R_{0202} + R_{1212} + R_{2323} = \frac{U}{\rho^6}\left(-1 - 1 + 2\right) = 0,$$

and $R_{33}$ the same as $R_{22}$ by the symmetry between $e^2$ and $e^3$.

The six off diagonal entries are shorter still.
Each of them is a sum of frame components whose two index pairs are different and share one index, such as $(20)(21)$ or $(21)(23)$, and not one of those appears among the eight: for example

$$R_{01} = R_{2021} + R_{3031} = 0 + 0, \qquad R_{13} = -R_{0103} + R_{2123} = 0 + 0.$$

The three components carrying $W$ pair $(01)$ with $(23)$, $(02)$ with $(13)$ and $(03)$ with $(12)$, which share no index at all, so they never reach an off diagonal slot of Ricci.
That is the statement that $B_{ij}$ is symmetric, in components.
Ricci is zero in all ten slots.

This is the field equation, and it is what fixes $f$.
Had we written $f = (r^2 - 2mr + k)/(r^2+l^2)$ with $k$ free and left it to the vacuum condition, $R_{00} = 0$ would have come out as a rational function of $r$ whose numerator has exactly two roots in $k$: the constant $k = -l^2$, and $k = r(2m-r)$, which is not a constant and only says $f \equiv 0$.
The function in Step 2 is not a guess that happens to work but the only one of its form that does.
So in the chart,

$$R_{\mu\nu} = 0,$$

and the entry publishes an empty block for it in all three variants, $R_{\mu\nu}$, $R^\mu{}_\nu$ and $R^{\mu\nu}$.
The sign convention of the contraction cannot be read off an empty block, which is why Step 1 records it in prose instead.

---

## Step 9. The Ricci scalar and the Einstein tensor

$$R = g^{\mu\nu}R_{\mu\nu} = 0, \qquad G_{\mu\nu} = R_{\mu\nu} - \tfrac{1}{2}Rg_{\mu\nu} = 0.$$

Both follow from Step 8 with no further work, and the Einstein tensor is published empty in all three variants as well.
There is no stress energy anywhere: $T_{\mu\nu} = 0$, the spacetime is empty, and everything the entry prints below is curvature of the vacuum sort, sourced by nothing at any point of the manifold the chart covers.

---

## Step 10. The Weyl tensor, computed rather than copied

The Weyl tensor in four dimensions is Riemann with its traces removed:

$$C_{\mu\nu\rho\sigma} = R_{\mu\nu\rho\sigma} - \frac{1}{2}\left(g_{\mu\rho}\mathcal{R}_{\sigma\nu} - g_{\mu\sigma}\mathcal{R}_{\rho\nu} - g_{\nu\rho}\mathcal{R}_{\sigma\mu} + g_{\nu\sigma}\mathcal{R}_{\rho\mu}\right) + \frac{\mathcal{R}}{6}\left(g_{\mu\rho}g_{\sigma\nu} - g_{\mu\sigma}g_{\rho\nu}\right),$$

where $\mathcal{R}_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$ and $\mathcal{R} = g^{\mu\nu}\mathcal{R}_{\mu\nu}$ are the traces of Riemann itself.
Those traces are the ones Weyl is defined by, whatever a file chooses to call its Ricci tensor, so this formula is convention free.

The computation was run as written, on all two hundred and fifty six components, with $\mathcal{R}_{\mu\nu}$ and $\mathcal{R}$ recomputed from the published Riemann tensor rather than assumed.
Both correction terms came out identically zero, because both are built from $\mathcal{R}_{\mu\nu}$, which Step 8 showed vanishes.
Hence

$$C_{\mu\nu\rho\sigma} = R_{\mu\nu\rho\sigma}$$

in every one of the seventy two nonzero slots, and $C^\mu{}_{\nu\rho\sigma} = R^\mu{}_{\nu\rho\sigma}$ in every one of the eighty.

The equality is a theorem about vacuum solutions and not a shortcut through the arithmetic.
It says something: all of the curvature of Taub-NUT is Weyl curvature, which is to say tidal and free, none of it traceable to matter at the point where it is felt.
The published Weyl block is therefore numerically identical to the published Riemann block, and the check that it should be was the computation above, not the copy.

In the frame of Step 7 the same statement reads $E_{ij} + iB_{ij} = \mathfrak{c}\,\mathrm{diag}(-2,1,1)$ with

$$\mathfrak{c} = \frac{m+il}{(r+il)^3},$$

a single complex Weyl scalar, which is the compact form of Petrov type D.

---

## Step 11. The Kretschmann scalar

$K = R_{\mu\nu\rho\sigma}R^{\mu\nu\rho\sigma}$ is easiest in the frame, where raising an index costs a minus sign for every $0$ it carries.

The purely electric part contributes twice.
Each component $R_{0i0i}$ appears four times in the unrestricted sum over the four indices, and $R^{0i0i} = (-1)^2R_{0i0i} = R_{0i0i}$, so

$$4\left(R_{0101}^2 + R_{0202}^2 + R_{0303}^2\right) = 4\left(4+1+1\right)\frac{U^2}{\rho^{12}} = \frac{24U^2}{\rho^{12}},$$

and the purely spatial components $R_{1212}$, $R_{1313}$, $R_{2323}$ appear four times each with no sign at all,

$$4\left(1 + 1 + 4\right)\frac{U^2}{\rho^{12}} = \frac{24U^2}{\rho^{12}}.$$

The magnetic part contributes with the opposite sign.
Each of $R_{0123}$, $R_{0213}$, $R_{0312}$ appears eight times, and raising one $0$ index gives $R^{0ijk} = -R_{0ijk}$, so

$$-8\left(4 + 1 + 1\right)\frac{W^2}{\rho^{12}} = -\frac{48W^2}{\rho^{12}}.$$

Adding the three,

$$K = \frac{48\left(U^2 - W^2\right)}{\rho^{12}} = \frac{48\left(U-W\right)\left(U+W\right)}{\left(r^2+l^2\right)^6},$$

which by Step 7 is

$$K = 48\,\mathrm{Re}\left[\left(\frac{m+il}{(r+il)^3}\right)^2\right] = 48\,\mathrm{Re}\left[\frac{(m+il)^2}{(r+il)^6}\right].$$

The two factors, written out, are what the entry publishes:

$$U - W = mr^3 - lr^3 + 3lmr^2 + 3l^2r^2 - 3l^2mr + 3l^3r - l^3m - l^4,$$

$$U + W = mr^3 + lr^3 - 3lmr^2 + 3l^2r^2 - 3l^2mr - 3l^3r + l^3m - l^4.$$

At $l = 0$ both factors are $mr^3$ and $\rho^{12} = r^{12}$, so

$$K \to \frac{48m^2}{r^6} = \frac{48G^2M^2}{c^4r^6},$$

which is Schwarzschild's, as it has to be.

---

## Step 12. What the Kretschmann scalar says

$K$ has denominator $\left(r^2+l^2\right)^6$, which for $l \neq 0$ is never zero for any real $r$.
There is no curvature singularity anywhere in Taub-NUT.
That is the first of the spacetime's several refusals to behave: a solution with a mass parameter, a horizon at $r_+ = m + \sqrt{m^2+l^2}$, and nothing at the centre for the mass to be.
Continue $r$ down through the horizon and through zero and out the far side and the curvature stays bounded the whole way.

$K$ is also finite on the axis, at $\theta = 0$ and $\theta = \pi$, being independent of $\theta$ altogether.
The Misner string of Step 14 is therefore not a curvature singularity and cannot be one; it is a defect of the chart, and Step 14 says which.

The second thing $K$ says is that it changes sign.
Expanding the quotient of Step 11 for large $r$,

$$K = \frac{48\left(m^2 - l^2\right)}{r^6} + O\!\left(\frac{1}{r^7}\right),$$

so once $l > m$ the Kretschmann scalar is negative far from the source, which $48m^2/r^6$ never is.
Both $U-W$ and $U+W$ are cubics in $r$ with real roots for any $l \neq 0$, so $K$ passes through zero at finite radius in every case; whether any of those radii lies in the exterior is what $l$ against $m$ decides.
For $l > m$ exactly one of them does, and for $l < m$ they all sit at or inside the horizon and $K$ stays positive over the whole published domain.
At a radius where $K$ vanishes the magnetic part has caught the electric part up, $|E| = |B|$, and the tidal field is null in the same algebraic sense a null electromagnetic field is.
Nothing of the kind happens in Schwarzschild, and the difference is entirely the NUT parameter's.

---

## Step 13. The geodesic equations

The equations are

$$\ddot{x}^\mu + \Gamma^\mu{}_{\nu\rho}\dot{x}^\nu\dot{x}^\rho = 0$$

with the $\Gamma$ of Step 4 and the dot the derivative with respect to an affine parameter $\lambda$ along the curve, taken of the chart coordinates, so $\dot{t} = d(ct)/d\lambda$ and $\dot{r} = dr/d\lambda$.
Summing over both orders of $\nu$ and $\rho$ doubles every symbol whose lower pair is unequal, which is where the factors of two below come from.

$$\ddot{t} + \frac{f'}{f}\dot{t}\dot{r} - \frac{4l^2f\cos\theta}{\rho^2\sin\theta}\dot{t}\dot{\theta} - \frac{4W\cos\theta}{\rho^2\Delta}\dot{r}\dot{\phi} + \frac{2l\left[\left(\rho^4+4l^2\Delta\right)\sin^2\theta - 2\left(\rho^4+2l^2\Delta\right)\right]}{\rho^4\sin\theta}\dot{\theta}\dot{\phi} = 0,$$

$$\ddot{r} + \frac{ff'}{2}\dot{t}^2 + 2lff'\cos\theta\,\dot{t}\dot{\phi} - \frac{f'}{2f}\dot{r}^2 - rf\dot{\theta}^2 + \left(2l^2ff'\cos^2\theta - rf\sin^2\theta\right)\dot{\phi}^2 = 0,$$

$$\ddot{\theta} - \frac{2lf\sin\theta}{\rho^2}\dot{t}\dot{\phi} + \frac{2r}{\rho^2}\dot{r}\dot{\theta} - \left(1 + \frac{4l^2f}{\rho^2}\right)\sin\theta\cos\theta\,\dot{\phi}^2 = 0,$$

$$\ddot{\phi} + \frac{2lf}{\rho^2\sin\theta}\dot{t}\dot{\theta} + \frac{2r}{\rho^2}\dot{r}\dot{\phi} + 2\left(1 + \frac{2l^2f}{\rho^2}\right)\cot\theta\,\dot{\theta}\dot{\phi} = 0.$$

The entry prints these with $f$, $\rho^2$, $\Delta$ and $W$ substituted out.

The third equation is the one to read.
Set $l = 0$ and its first term goes, leaving the Schwarzschild statement that a particle launched in the equatorial plane with $\dot{\theta} = 0$ and $\dot{\phi} \neq 0$ stays in it.
With $l \neq 0$ the term $-2lf\sin\theta\,\dot{t}\dot{\phi}/\rho^2$ survives, and any particle with $\dot{\phi} \neq 0$ is pushed off the plane at once.
There are no planar orbits in Taub-NUT.
Orbits lie instead on cones about the axis, which is the gravitational copy of the result that a charge moving in a magnetic monopole field is confined to a cone, and the same conserved quantity is behind both.

---

## Step 14. The Misner string

$\theta$ runs over the open interval $(0, \pi)$ in the published domains, and the exclusion of the two endpoints is not tidiness.

On a sphere $d\phi$ is undefined at both poles, so a one form $h(\theta)d\phi$ extends to the pole only if $h$ vanishes there.
The potential of Step 2 has $h = 2l\cos\theta$, which is $+2l$ at $\theta = 0$ and $-2l$ at $\theta = \pi$.
It vanishes at neither.
So in this gauge the metric is singular along the whole axis, and the singularity is a coordinate one: Step 12 showed $K$ is finite there.
The singular half axis is the Misner string, the gravitational counterpart of the Dirac string that trails a magnetic monopole, and it is there for the same reason.

The gauge can be moved but not removed.
Adding the exact form $\mp 2l\,d\phi$ gives

$$A_\pm = 2l\left(\cos\theta \mp 1\right)d\phi,$$

and $A_+$ vanishes at $\theta = 0$ while $A_-$ vanishes at $\theta = \pi$.
Either choice clears one half of the axis and doubles the trouble on the other.
This is exactly the two patch construction Wu and Yang use for the monopole, and it works there because the transition function between the patches is a phase.
Here the transition is a shift of the time coordinate, forced by asking the two patches to describe the same $\sigma = c\,dt + A$,

$$c\,t_+ = c\,t_- + 4l\phi,$$

and asking it to be single valued as $\phi$ goes around by $2\pi$ forces $t$ to be periodic,

$$t \sim t + \frac{8\pi l}{c}.$$

That is Misner's repair, and it is worse than the disease.
A periodic time coordinate whose orbits are timelike somewhere makes closed timelike curves through every point of the region where they are, so the spacetime is regular and causally hopeless, or causally sensible and singular along a half axis, and there is no third option.
The entry publishes the second, with $\theta \in (0,\pi)$ in the domains and the whole story in the convention field, because the components printed are the components of a chart, and a chart that covers the axis in this gauge does not exist.

---

## Step 15. Every published expression is dimensionally consistent

The chart coordinates are $(ct, r, \theta, \phi)$, of dimensions $L$, $L$, $1$, $1$.
The parameters $m$ and $l$ are both lengths.
A metric component then carries

$$[g_{\mu\nu}] = \frac{L^2}{[x^\mu][x^\nu]},$$

a Christoffel symbol carries $L^{-1}$ corrected by $[x^\mu]/L$ per upper index and $L/[x^\mu]$ per lower one, a Riemann or Weyl component carries $L^{-2}$ corrected the same way, and $K$ carries $L^{-4}$.

Three samples, one per rank.

$g_{t\phi} = -2lf\cos\theta$ must carry $L^2/(L\cdot 1) = L$, and $f$ is dimensionless while $l$ is a length, so it carries $L$.
This is the component that decides the whole scheme: if $m$ and $l$ were not lengths there would be no way for a $dt\,d\phi$ term to balance against a $d\theta^2$ term in the same line element.

$\Gamma^\theta{}_{t\phi} = -lf\sin\theta/\rho^2$ carries $L^{-1}$ from the field, times $1/L$ for the upper $\theta$, times $L/L = 1$ for the lower $t$, times $L/1 = L$ for the lower $\phi$, so $L^{-1}$ in all; and $l/\rho^2$ is $L/L^2 = L^{-1}$.

$K = 48(U^2-W^2)/\rho^{12}$ must carry $L^{-4}$, and $U$ and $W$ are quartics in lengths, so $U^2$ carries $L^8$ against $\rho^{12} = L^{12}$.

The geodesic equations are measured against their own second derivatives.
In the first of them $\ddot{t}$ carries $L/\lambda^2$, and the term $-4W\cos\theta\,\dot{r}\dot{\phi}/\rho^2\Delta$ carries

$$\frac{L^4}{L^2\cdot L^2}\cdot\frac{L}{\lambda}\cdot\frac{1}{\lambda} = \frac{L}{\lambda^2},$$

which matches, and matches only because $\dot{t}$ is $d(ct)/d\lambda$ and not $dt/d\lambda$.
On the other reading the first term of that equation would come out one factor of $c$ away from the rest, and the dimensional pass in `verify_metrics.py` would name it.
It does not, on any of the three hundred and seventy three expressions the entry publishes.

---

## Step 16. What the entry publishes, and what the checker needs

The file carries one coordinate system, `spherical`, and in it:

the line element; the six nonzero $g_{\mu\nu}$ and the six nonzero $g^{\mu\nu}$; the twenty five $\Gamma^\mu{}_{\nu\rho}$ and the twenty five $\Gamma_{\mu\nu\rho}$; the eighty $R^\mu{}_{\nu\rho\sigma}$ and the seventy two $R_{\mu\nu\rho\sigma}$; empty blocks for $R_{\mu\nu}$, $R^\mu{}_\nu$, $R^{\mu\nu}$, $G_{\mu\nu}$, $G^\mu{}_\nu$ and $G^{\mu\nu}$; $R = 0$; the Kretschmann scalar; the eighty $C^\mu{}_{\nu\rho\sigma}$ and the seventy two $C_{\mu\nu\rho\sigma}$; and the four geodesic equations.

For `verify_metrics.py` to read it, the system needs one line in `DIMENSIONS`,

    ("taub_nut", "spherical"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "m": "L", "l": "L",
    },

declaring $t$ a time, which is what tells the script that the chart multiplies it by $c$, and declaring both parameters lengths, which is what lets the dimensional pass weigh the cross term.
No entry in `PARAMETER_RELATIONS` is needed: $m$ and $l$ are free, and every value published is an identity in both.

The one thing to know before rerunning the check is that this system is slow.
Every published value agrees with sympy, but the metric has an off diagonal block, so the inverse is dense in the $t$ and $\phi$ rows and every contraction below Riemann carries the extra terms.
None of the tensors from Riemann down reliably finishes inside the default per tensor budget of one hundred and twenty seconds, and the Kretschmann scalar never does.
The budget is wall clock, so how much of the entry comes back `UNCHECKED` at the default depends on what else the machine is doing.
Verify it at a larger budget instead:

    /tmp/mfs-venv/bin/python _tools/derivations/verify_metrics.py \
        --system taub_nut/spherical --budget 2400

which reports no disagreements and no dimensional failures, and takes from nine minutes to about half an hour depending on the load.
The dimensional pass is unaffected and runs over the whole entry in a moment, as it does everywhere else.
