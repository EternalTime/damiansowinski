# The Vaidya radiating star

This is the working behind the two coordinate systems in `MFS/assets/data/metrics/vaidya.json`.
Every number the entry prints is derived here, in order, from the line element down to the geodesic equations.
Nothing is left as an exercise and nothing is asserted that is not computed.

The companion script `_tools/derivations/verify_metrics.py` does the same work in sympy and compares it against the published file, so the algebra below is checkable by hand and by machine independently.

Vaidya is the first entry in the collection that is not a vacuum.
Its Ricci tensor and its Einstein tensor do not vanish, and the single component each of them has is the whole physical content of the solution: a flux of pure radiation streaming along null rays.
Step 8 is where that component appears, Step 10 reads it as a stress energy tensor, and Step 12 shows why no scalar invariant can see it at all.

---

## Step 1. Conventions

The signature is $(-,+,+,+)$ and the Christoffel symbols are those of the Levi-Civita connection,

$$\Gamma^\mu_{\nu\rho} = \tfrac{1}{2}g^{\mu\alpha}\left(\partial_\nu g_{\alpha\rho} + \partial_\rho g_{\alpha\nu} - \partial_\alpha g_{\nu\rho}\right),$$

with the lowered symbol $\Gamma_{\mu\nu\rho} = g_{\mu\alpha}\Gamma^\alpha_{\nu\rho}$.
The Riemann tensor is

$$R^\mu{}_{\nu\rho\sigma} = \partial_\rho \Gamma^\mu_{\nu\sigma} - \partial_\sigma \Gamma^\mu_{\nu\rho} + \Gamma^\mu_{\rho\lambda}\Gamma^\lambda_{\nu\sigma} - \Gamma^\mu_{\sigma\lambda}\Gamma^\lambda_{\nu\rho},$$

which is what the published Riemann components are in.

The Ricci tensor is the standard contraction on the first lower index,

$$R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu},$$

which is what the collection publishes everywhere.
Here the choice matters, because the Ricci tensor does not vanish.
Everything below is computed with it, contracting on the last index instead gives the negative of every Ricci and Einstein component, and Step 10 carries the sign through to the field equations, which on this convention read $G_{\mu\nu} = 8\pi G\,T_{\mu\nu}/c^4$.

Factors of $G$ and $c$ are kept explicit.
The chart, here and everywhere else in the collection, is the one whose time coordinate is $x^0 = cT$, so here it is $x^0 = cu$ in the outgoing chart and $x^0 = cv$ in the ingoing one.
The index is printed with the bare letter $u$, but the component printed against it is a component of the chart in which the zeroth coordinate is $cu$ and therefore carries a length.

This has one consequence worth stating before any algebra is done.
A dot is a derivative with respect to the chart coordinate, not with respect to the bare one:

$$\dot{m} \equiv \frac{dm}{d(cu)} = \frac{1}{c}\frac{dm}{du},$$

so $\dot{m}$ is a mass per unit length rather than a mass per unit time.
The same reading applies to the dots in the geodesic equations, where $\dot{u}$ means $d(cu)/d\lambda$.
It is the reading that makes every term of every published equation carry the dimensions of its left hand side, which Step 15 checks term by term.

Everything from Step 2 to Step 15 is done directly in the chart $x^0 = cu$, so no conversion step is needed at the end: the components computed below are the components the entry prints.
Where an index is written $0$ it is the chart index the entry prints as $u$.

---

## Step 2. The line element

In outgoing, or retarded, Eddington-Finkelstein form,

$$ds^2 = -\left(1 - \frac{2Gm(u)}{c^2r}\right)c^2du^2 - 2c\,du\,dr + r^2d\theta^2 + r^2\sin^2\theta\,d\phi^2.$$

The coordinate $u$ is the retarded time, labelling the outgoing light fronts; $r$ is the areal radius, so that the sphere of coordinate radius $r$ has area $4\pi r^2$; and $m(u)$ is the mass the metric still reports after what has already radiated away.

Two abbreviations carry the rest of the derivation.
Write the instantaneous Schwarzschild radius and the redshift factor as

$$r_s(u) \equiv \frac{2Gm(u)}{c^2}, \qquad f(u,r) \equiv 1 - \frac{r_s(u)}{r} = 1 - \frac{2Gm(u)}{c^2r},$$

so that the line element is $ds^2 = -f\,(c\,du)^2 - 2(c\,du)\,dr + r^2d\Omega^2$.
Both are lengths and dimensionless respectively, and the derivative of the first with respect to the chart coordinate is

$$\dot{r}_s = \frac{dr_s}{d(cu)} = \frac{2G\dot{m}}{c^2},$$

which is dimensionless, since $\dot m$ is a mass per unit length.

Holding $m$ constant returns the outgoing Eddington-Finkelstein chart of the Schwarzschild entry, with $r_s = 2Gm/c^2$.
Every formula below reduces to that entry's when $\dot{m} = 0$, and every term carrying $\dot m$ is a term Schwarzschild does not have.
That is the practical form of the remark in the entry's history: Birkhoff's theorem forces an empty spherical exterior to be static, so a spherical exterior that is not static cannot be empty.

The chart is $u \in (-\infty,\infty)$, $r \in (0,\infty)$, $\theta \in [0,\pi]$, $\phi \in [0,2\pi)$.
Unlike the Schwarzschild chart in $(t,r)$, nothing goes wrong at $r = r_s$; Step 3 shows the determinant does not know that surface is there.

---

## Step 3. The metric matrix, its determinant and its inverse

Reading the coefficients off the line element in the order $(x^0, r, \theta, \phi)$ with $x^0 = cu$, and remembering that the cross term $-2(c\,du)dr$ contributes $g_{0r} = g_{r0} = -1$ because it appears twice in the quadratic form,

$$g_{\mu\nu} = \begin{pmatrix} -f & -1 & 0 & 0 \\ -1 & 0 & 0 & 0 \\ 0 & 0 & r^2 & 0 \\ 0 & 0 & 0 & r^2\sin^2\theta \end{pmatrix}.$$

These four are the entry's `metric_components`: $g_{uu} = -\left(1-\dfrac{2Gm}{c^2r}\right)$, $g_{ur} = g_{ru} = -1$, $g_{\theta\theta} = r^2$ and $g_{\phi\phi} = r^2\sin^2\theta$.

The matrix is block diagonal, so its determinant is the product of the determinants of the two blocks,

$$\det g = \left[(-f)(0) - (-1)(-1)\right]\cdot r^2\cdot r^2\sin^2\theta = -r^4\sin^2\theta.$$

The mass function has cancelled out of it entirely.
So $\sqrt{-g} = r^2\sin\theta$, the same volume element as flat space in spherical coordinates, and it is nonzero for every $r > 0$ away from the poles: the chart is degenerate nowhere, and in particular not at $r = r_s$.

The inverse follows from inverting the upper block, $\begin{pmatrix} -f & -1 \\ -1 & 0\end{pmatrix}^{-1} = \begin{pmatrix} 0 & -1 \\ -1 & f \end{pmatrix}$, and inverting the two diagonal entries:

$$g^{\mu\nu} = \begin{pmatrix} 0 & -1 & 0 & 0 \\ -1 & f & 0 & 0 \\ 0 & 0 & r^{-2} & 0 \\ 0 & 0 & 0 & r^{-2}\sin^{-2}\theta \end{pmatrix}.$$

Three of the four checks of $g^{\mu\alpha}g_{\alpha\nu} = \delta^\mu{}_\nu$ on the upper block are worth writing out, since the vanishing $g^{00}$ is what most of the later algebra turns on:

$$g^{0\alpha}g_{\alpha 0} = g^{0r}g_{r0} = (-1)(-1) = 1, \qquad g^{0\alpha}g_{\alpha r} = g^{00}g_{0r} + g^{0r}g_{rr} = 0,$$

$$g^{r\alpha}g_{\alpha 0} = g^{r0}g_{00} + g^{rr}g_{r0} = (-1)(-f) + f(-1) = 0, \qquad g^{r\alpha}g_{\alpha r} = g^{r0}g_{0r} = 1.$$

So the entry's `inverse_metric_components` are $g^{ur} = g^{ru} = -1$, $g^{rr} = 1 - \dfrac{2Gm}{c^2r}$, $g^{\theta\theta} = r^{-2}$ and $g^{\phi\phi} = r^{-2}\sin^{-2}\theta$, with $g^{uu} = 0$.

That last fact, $g^{00} = 0$, says that the surfaces of constant $u$ are null.
It is the reason the Ricci scalar vanishes in Step 9 even though the Ricci tensor does not.

---

## Step 4. The lowered Christoffel symbols

The lowered symbol is

$$\Gamma_{\mu\nu\rho} = \tfrac{1}{2}\left(\partial_\nu g_{\mu\rho} + \partial_\rho g_{\mu\nu} - \partial_\mu g_{\nu\rho}\right),$$

symmetric in its last two indices.
Only five derivatives of the metric are nonzero:

$$\partial_0 g_{00} = -\partial_0 f = \frac{\dot{r}_s}{r}, \qquad \partial_r g_{00} = -\partial_r f = -\frac{r_s}{r^2},$$

$$\partial_r g_{\theta\theta} = 2r, \qquad \partial_r g_{\phi\phi} = 2r\sin^2\theta, \qquad \partial_\theta g_{\phi\phi} = 2r^2\sin\theta\cos\theta.$$

Note that $\partial_0$ reaches only $g_{00}$, since $r_s$ is the only place the retarded time appears, and that $\partial_0 f = -\dot r_s/r$.

Taking the cases in turn:

$$\Gamma_{000} = \tfrac{1}{2}\partial_0 g_{00} = \frac{\dot{r}_s}{2r} = \frac{G\dot{m}}{c^2r},$$

$$\Gamma_{00r} = \Gamma_{0r0} = \tfrac{1}{2}\left(\partial_0 g_{0r} + \partial_r g_{00} - \partial_0 g_{0r}\right) = \tfrac{1}{2}\partial_r g_{00} = -\frac{r_s}{2r^2} = -\frac{Gm}{c^2r^2},$$

$$\Gamma_{r00} = \tfrac{1}{2}\left(2\partial_0 g_{r0} - \partial_r g_{00}\right) = \frac{r_s}{2r^2} = \frac{Gm}{c^2r^2},$$

$$\Gamma_{0rr} = \tfrac{1}{2}\left(2\partial_r g_{0r} - \partial_0 g_{rr}\right) = 0, \qquad \Gamma_{rrr} = \tfrac{1}{2}\partial_r g_{rr} = 0,$$

$$\Gamma_{r\theta\theta} = -\tfrac{1}{2}\partial_r g_{\theta\theta} = -r, \qquad \Gamma_{\theta r\theta} = \Gamma_{\theta\theta r} = \tfrac{1}{2}\partial_r g_{\theta\theta} = r,$$

$$\Gamma_{r\phi\phi} = -\tfrac{1}{2}\partial_r g_{\phi\phi} = -r\sin^2\theta, \qquad \Gamma_{\phi r\phi} = \Gamma_{\phi\phi r} = \tfrac{1}{2}\partial_r g_{\phi\phi} = r\sin^2\theta,$$

$$\Gamma_{\theta\phi\phi} = -\tfrac{1}{2}\partial_\theta g_{\phi\phi} = -r^2\cos\theta\sin\theta, \qquad \Gamma_{\phi\theta\phi} = \Gamma_{\phi\phi\theta} = \tfrac{1}{2}\partial_\theta g_{\phi\phi} = r^2\cos\theta\sin\theta.$$

Every other lowered symbol vanishes, because it would need a derivative that is not on the list.
These thirteen are exactly the entry's `christoffel.lll` block.

The only symbol carrying $\dot{m}$ is $\Gamma_{000}$.
That is the whole of the difference from Schwarzschild at this level, and it is where the radiation will come from.

---

## Step 5. The Christoffel symbols with an upper index

Raise with $\Gamma^\mu_{\nu\rho} = g^{\mu\alpha}\Gamma_{\alpha\nu\rho}$.
Because $g^{0\alpha}$ is nonzero only for $\alpha = r$, and $g^{r\alpha}$ only for $\alpha \in \{0,r\}$, the four rules are

$$\Gamma^0_{\nu\rho} = -\Gamma_{r\nu\rho}, \qquad \Gamma^r_{\nu\rho} = -\Gamma_{0\nu\rho} + f\,\Gamma_{r\nu\rho}, \qquad \Gamma^\theta_{\nu\rho} = \frac{\Gamma_{\theta\nu\rho}}{r^2}, \qquad \Gamma^\phi_{\nu\rho} = \frac{\Gamma_{\phi\nu\rho}}{r^2\sin^2\theta}.$$

Applying them to the thirteen symbols of Step 4:

$$\Gamma^0_{00} = -\Gamma_{r00} = -\frac{r_s}{2r^2} = -\frac{Gm}{c^2r^2}, \qquad \Gamma^0_{\theta\theta} = -\Gamma_{r\theta\theta} = r, \qquad \Gamma^0_{\phi\phi} = -\Gamma_{r\phi\phi} = r\sin^2\theta,$$

$$\Gamma^r_{00} = -\Gamma_{000} + f\,\Gamma_{r00} = \frac{r_s}{2r^2}f - \frac{\dot{r}_s}{2r} = \frac{Gm}{c^2r^2}\left(1-\frac{2Gm}{c^2r}\right) - \frac{G\dot{m}}{c^2r},$$

$$\Gamma^r_{0r} = \Gamma^r_{r0} = -\Gamma_{00r} + f\,\Gamma_{r0r} = \frac{r_s}{2r^2} = \frac{Gm}{c^2r^2},$$

$$\Gamma^r_{\theta\theta} = f\,\Gamma_{r\theta\theta} = -rf = -r\left(1-\frac{2Gm}{c^2r}\right), \qquad \Gamma^r_{\phi\phi} = -rf\sin^2\theta,$$

$$\Gamma^\theta_{r\theta} = \Gamma^\theta_{\theta r} = \frac{1}{r}, \qquad \Gamma^\theta_{\phi\phi} = -\cos\theta\sin\theta, \qquad \Gamma^\phi_{r\phi} = \Gamma^\phi_{\phi r} = \frac{1}{r}, \qquad \Gamma^\phi_{\theta\phi} = \Gamma^\phi_{\phi\theta} = \cot\theta.$$

Those fifteen are the entry's `christoffel.ull` block, and there are no others.

Two features of the list are used again and again below.
There is no symbol of the form $\Gamma^\mu_{rr}$, for any $\mu$, which Step 14 turns into the statement that the outgoing rays are affinely parametrised geodesics.
And $\Gamma^r_{00}$ is the only symbol carrying $\dot{m}$, the trace of $\Gamma_{000}$ after the raising.

---

## Step 6. The Riemann tensor, computed component by component

The independent components of a spherically symmetric metric of this form sit in four planes: the $(0,r)$ plane, the two mixed planes $(0,\theta)$ and $(r,\theta)$ with their $\phi$ copies, and the $(\theta,\phi)$ plane.
Six direct computations from

$$R^\mu{}_{\nu\rho\sigma} = \partial_\rho \Gamma^\mu_{\nu\sigma} - \partial_\sigma \Gamma^\mu_{\nu\rho} + \Gamma^\mu_{\rho\lambda}\Gamma^\lambda_{\nu\sigma} - \Gamma^\mu_{\sigma\lambda}\Gamma^\lambda_{\nu\rho}$$

give everything, and the rest of the published components follow by the symmetries.

**$R^0{}_{00r}$.**
The two connection terms are $\Gamma^0_{0\lambda}\Gamma^\lambda_{0r}$, which needs $\Gamma^0_{0r}$ or $\Gamma^0_{0\theta}$ and finds neither, and $\Gamma^0_{r\lambda}\Gamma^\lambda_{00}$, which needs a symbol $\Gamma^0_{r\lambda}$ and finds none.
With $\Gamma^0_{0r} = 0$ the first derivative term goes too, leaving

$$R^0{}_{00r} = -\partial_r\Gamma^0_{00} = -\partial_r\left(-\frac{r_s}{2r^2}\right) = -\frac{r_s}{r^3} = -\frac{2Gm}{c^2r^3}.$$

**$R^r{}_{r0r}$.**
Here $\Gamma^r_{rr} = 0$ kills the first derivative term and every term $\Gamma^\lambda_{rr}$ vanishes, so

$$R^r{}_{r0r} = -\partial_r\Gamma^r_{r0} - \Gamma^r_{r\lambda}\Gamma^\lambda_{r0} = -\partial_r\left(\frac{r_s}{2r^2}\right) - \Gamma^r_{r0}\Gamma^0_{r0} = \frac{r_s}{r^3} = \frac{2Gm}{c^2r^3},$$

using $\Gamma^0_{r0} = 0$ for the last term.

**$R^r{}_{00r}$.**
This is the one with cancellations.
The derivative terms are

$$\partial_0\Gamma^r_{0r} = \frac{\dot{r}_s}{2r^2}, \qquad -\partial_r\Gamma^r_{00} = -\partial_r\left(\frac{r_s}{2r^2} - \frac{r_s^2}{2r^3} - \frac{\dot{r}_s}{2r}\right) = \frac{r_s}{r^3} - \frac{3r_s^2}{2r^4} - \frac{\dot{r}_s}{2r^2},$$

and the connection terms are $\Gamma^r_{0\lambda}\Gamma^\lambda_{0r} = \left(\Gamma^r_{0r}\right)^2 = \dfrac{r_s^2}{4r^4}$ and $-\Gamma^r_{r\lambda}\Gamma^\lambda_{00} = -\Gamma^r_{r0}\Gamma^0_{00} = \dfrac{r_s^2}{4r^4}$.
The two terms in $\dot{r}_s$ cancel against each other and the three in $r_s^2$ combine, leaving

$$R^r{}_{00r} = \frac{r_s}{r^3} - \frac{r_s^2}{r^4} = \frac{r_s}{r^3}f = \frac{2Gm}{c^2r^3}\left(1-\frac{2Gm}{c^2r}\right).$$

**$R^0{}_{\theta 0\theta}$.**
Only one term survives, the connection term $\Gamma^0_{00}\Gamma^0_{\theta\theta}$:

$$R^0{}_{\theta 0\theta} = \Gamma^0_{0\lambda}\Gamma^\lambda_{\theta\theta} = \left(-\frac{r_s}{2r^2}\right)(r) = -\frac{r_s}{2r} = -\frac{Gm}{c^2r}.$$

**$R^\theta{}_{00\theta}$ and $R^\theta{}_{0r\theta}$.**
Every term with an upper $\theta$ and a lower $0$ vanishes, so both reduce to a single product with $\Gamma^\theta_{\theta r} = 1/r$:

$$R^\theta{}_{00\theta} = -\Gamma^\theta_{\theta\lambda}\Gamma^\lambda_{00} = -\frac{1}{r}\Gamma^r_{00} = -\frac{r_s}{2r^3}f + \frac{\dot{r}_s}{2r^2} = -\frac{Gm}{c^2r^3}\left(1-\frac{2Gm}{c^2r}\right) + \frac{G\dot{m}}{c^2r^2},$$

$$R^\theta{}_{0r\theta} = -\Gamma^\theta_{\theta\lambda}\Gamma^\lambda_{0r} = -\frac{1}{r}\Gamma^r_{0r} = -\frac{r_s}{2r^3} = -\frac{Gm}{c^2r^3}.$$

**$R^r{}_{\theta 0\theta}$ and $R^r{}_{\theta r\theta}$.**
For the first, $\partial_0\Gamma^r_{\theta\theta} = \partial_0(-r + r_s) = \dot{r}_s$, and the connection terms are

$$\Gamma^r_{0\lambda}\Gamma^\lambda_{\theta\theta} = \Gamma^r_{00}\Gamma^0_{\theta\theta} + \Gamma^r_{0r}\Gamma^r_{\theta\theta} = \left(\frac{r_s}{2r^2}f - \frac{\dot r_s}{2r}\right)r + \frac{r_s}{2r^2}(-rf) = -\frac{\dot{r}_s}{2},$$

with $-\Gamma^r_{\theta\lambda}\Gamma^\lambda_{\theta 0} = 0$ since $\Gamma^\theta_{\theta 0} = 0$.
Hence

$$R^r{}_{\theta 0\theta} = \dot{r}_s - \frac{\dot{r}_s}{2} = \frac{\dot{r}_s}{2} = \frac{G\dot{m}}{c^2}.$$

For the second, $\partial_r\Gamma^r_{\theta\theta} = \partial_r(-r+r_s) = -1$, the term $\Gamma^r_{r\lambda}\Gamma^\lambda_{\theta\theta} = \Gamma^r_{r0}\Gamma^0_{\theta\theta} = \dfrac{r_s}{2r}$, and $-\Gamma^r_{\theta\theta}\Gamma^\theta_{\theta r} = f$, so

$$R^r{}_{\theta r\theta} = -1 + \frac{r_s}{2r} + f = -\frac{r_s}{2r} = -\frac{Gm}{c^2r}.$$

**$R^\theta{}_{\phi\theta\phi}$.**
Here $\partial_\theta\Gamma^\theta_{\phi\phi} = \sin^2\theta - \cos^2\theta$, the term $\Gamma^\theta_{\theta r}\Gamma^r_{\phi\phi} = -f\sin^2\theta$, and $-\Gamma^\theta_{\phi\phi}\Gamma^\phi_{\phi\theta} = \cos^2\theta$, so

$$R^\theta{}_{\phi\theta\phi} = \sin^2\theta - \cos^2\theta - f\sin^2\theta + \cos^2\theta = (1-f)\sin^2\theta = \frac{r_s\sin^2\theta}{r} = \frac{2Gm\sin^2\theta}{c^2r}.$$

Each $\phi$ component is its $\theta$ partner times $\sin^2\theta$ where the index pattern calls for it, by spherical symmetry, and the entry's `riemann.ulll` block lists all thirty four nonzero components that the antisymmetry $R^\mu{}_{\nu\rho\sigma} = -R^\mu{}_{\nu\sigma\rho}$ generates from these.

---

## Step 7. The lowered Riemann tensor

Lowering with $R_{\mu\nu\rho\sigma} = g_{\mu\alpha}R^\alpha{}_{\nu\rho\sigma}$ uses only the first row of the metric, $g_{0\alpha}R^\alpha = -fR^0 - R^r$ and $g_{r\alpha}R^\alpha = -R^0$, together with $g_{\theta\theta} = r^2$ and $g_{\phi\phi} = r^2\sin^2\theta$.
The four independent components are

$$R_{0r0r} = -f R^0{}_{r0r} - R^r{}_{r0r} = -\frac{r_s}{r^3} = -\frac{2Gm}{c^2r^3},$$

$$R_{0\theta 0\theta} = -f R^0{}_{\theta 0\theta} - R^r{}_{\theta 0\theta} = \frac{r_s}{2r}f - \frac{\dot{r}_s}{2} = \frac{Gm}{c^2r}\left(1-\frac{2Gm}{c^2r}\right) - \frac{G\dot{m}}{c^2},$$

$$R_{0\theta r\theta} = -f R^0{}_{\theta r\theta} - R^r{}_{\theta r\theta} = \frac{r_s}{2r} = \frac{Gm}{c^2r},$$

$$R_{\theta\phi\theta\phi} = r^2 R^\theta{}_{\phi\theta\phi} = r_s r\sin^2\theta = \frac{2Gmr\sin^2\theta}{c^2},$$

using $R^0{}_{r0r} = 0$ and $R^0{}_{\theta r\theta} = 0$, neither of which has any surviving term in the Riemann formula.
Their $\phi$ partners are $R_{0\phi 0\phi} = \sin^2\theta\,R_{0\theta 0\theta}$ and $R_{0\phi r\phi} = \sin^2\theta\,R_{0\theta r\theta}$, and

$$R_{r\theta r\theta} = -R^0{}_{\theta r\theta} = 0,$$

which is the component a static chart would carry and this one does not.
The pair symmetry $R_{\mu\nu\rho\sigma} = R_{\rho\sigma\mu\nu}$ and the two antisymmetries generate the thirty two components the entry's `riemann.llll` block lists.

Only one of them, $R_{0\theta0\theta}$ and its $\phi$ copy, carries $\dot{m}$.
Every component with both pairs in the $(0,r)$ plane or both in the $(\theta,\phi)$ plane is exactly the Schwarzschild value with $r_s$ read at the retarded time.

---

## Step 8. The Ricci tensor, and where the radiation enters

With the contraction $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$, the $00$ component is a sum of four terms, of which the first vanishes by the antisymmetry in the last two indices:

$$R_{00} = R^0{}_{000} + R^r{}_{0r0} + R^\theta{}_{0\theta0} + R^\phi{}_{0\phi0}.$$

Step 6 computed all three of the survivors with their last two indices the other way about, so each enters with a minus, and the $\phi$ term equals the $\theta$ term:

$$R_{00} = -\frac{r_s}{r^3}f - 2\left(-\frac{r_s}{2r^3}f + \frac{\dot{r}_s}{2r^2}\right) = -\frac{\dot{r}_s}{r^2} = -\frac{2G\dot{m}}{c^2r^2}.$$

This is the heart of the entry.
The three terms carrying $f$ cancel exactly, which is Schwarzschild's vacuum condition surviving inside the radiating solution, and what is left is the single term the time dependence of the mass produced.
The Ricci tensor is proportional to $\dot{m}$ and to nothing else, so it vanishes at any retarded time when the star is not radiating and nowhere else.

Every other component vanishes, and the two cancellations are worth showing rather than asserting.
For $R_{0r}$,

$$R_{0r} = R^0{}_{00r} + R^r{}_{0rr} + R^\theta{}_{0\theta r} + R^\phi{}_{0\phi r} = -\frac{r_s}{r^3} + 0 + \frac{r_s}{2r^3} + \frac{r_s}{2r^3} = 0,$$

where $R^0{}_{00r} = -r_s/r^3$ and $R^r{}_{0rr} = 0$ by the antisymmetry in its last two indices.
For $R_{\theta\theta}$, the first two terms are the Step 6 values of $R^0{}_{\theta 0\theta}$ and $R^r{}_{\theta r\theta}$ as they stand, and the third is $R^\phi{}_{\theta\phi\theta} = R_{\theta\phi\theta\phi}/(r^2\sin^2\theta)$:

$$R_{\theta\theta} = R^0{}_{\theta 0\theta} + R^r{}_{\theta r\theta} + R^\phi{}_{\theta\phi\theta} = -\frac{r_s}{2r} - \frac{r_s}{2r} + \frac{r_s}{r} = 0.$$

The same cancellation kills $R_{\phi\phi}$, and $R_{rr} = 0$ because every one of its four terms would need a symbol $\Gamma^\mu_{rr}$, which Step 5 found does not exist.
The components with mixed angular indices vanish by spherical symmetry.

So the entry's `ricci_tensor.ll` block has exactly one entry,

$$R_{uu} = -\frac{2G\dot{m}}{c^2r^2},$$

and the two raised variants follow from $g^{0\alpha}$ having only the one nonzero entry $g^{0r} = -1$:

$$R^r{}_u = g^{r0}R_{00} = \frac{2G\dot{m}}{c^2r^2}, \qquad R^{rr} = \left(g^{r0}\right)^2R_{00} = -\frac{2G\dot{m}}{c^2r^2},$$

with $R^u{}_u = g^{u\alpha}R_{\alpha u} = g^{uu}R_{uu} = 0$ because $g^{uu} = 0$.
An index raised on a null tensor moves to a different slot rather than staying where it was, which is why the `ul` variant is printed against $[r,u]$ and the `uu` variant against $[r,r]$.

---

## Step 9. The Ricci scalar and the Einstein tensor

The Ricci scalar is

$$R = g^{\mu\nu}R_{\mu\nu} = g^{00}R_{00} = 0,$$

since $R_{00}$ is the only nonzero component and $g^{00} = 0$ from Step 3.
The scalar curvature of a radiating star is zero however hard it is radiating, which is the first sign that the source is null.

With a vanishing scalar the Einstein tensor is the Ricci tensor,

$$G_{\mu\nu} = R_{\mu\nu} - \tfrac{1}{2}Rg_{\mu\nu} = R_{\mu\nu},$$

so the entry publishes

$$G_{uu} = -\frac{2G\dot{m}}{c^2r^2}, \qquad G^r{}_u = \frac{2G\dot{m}}{c^2r^2}, \qquad G^{rr} = -\frac{2G\dot{m}}{c^2r^2},$$

and nothing else in any of the three variants.

---

## Step 10. The null dust, and reading $-\dot{m}$ as a luminosity

Let $k_\mu$ be minus the gradient of the chart's time coordinate,

$$k_\mu dx^\mu = -d(cu), \qquad k_\mu = (-1,0,0,0).$$

It is null, since $g^{\mu\nu}k_\mu k_\nu = g^{00} = 0$, and raising it gives

$$k^\mu = g^{\mu\nu}k_\nu = -g^{\mu 0} = (0,1,0,0), \qquad k^\mu\partial_\mu = \partial_r,$$

which points along increasing $r$ at fixed $u$: outgoing, as the retarded time demands.
It is a geodesic field with $r$ as an affine parameter, because

$$k^\nu\nabla_\nu k^\mu = \Gamma^\mu_{rr} = 0$$

for every $\mu$, which is the second use of the missing $\Gamma^\mu_{rr}$ noted in Step 5.

Since $k_\mu k_\nu$ has the single nonzero component $k_0k_0 = 1$, Step 9 can be rewritten with no loss as

$$G_{\mu\nu} = -\frac{2G\dot{m}}{c^2r^2}\,k_\mu k_\nu.$$

That is the pure radiation form: the Einstein tensor is a scalar times the outer square of one null covector.

The standard contraction makes the field equations $G_{\mu\nu} = 8\pi G\,T_{\mu\nu}/c^4$, as Step 1 recorded, so

$$T_{\mu\nu} = \frac{c^4}{8\pi G}G_{\mu\nu} = -\frac{c^2\dot{m}}{4\pi r^2}\,k_\mu k_\nu = -\frac{c}{4\pi r^2}\frac{dm}{du}k_\mu k_\nu.$$

Both signs moved together, so the stress energy is what it was: it is a measurable thing and cannot depend on which trace of Riemann is named Ricci.

This is null dust: energy streaming along the null geodesics of $k$ with no rest mass, no pressure transverse to the flow and no shear.

The sign is the physics.
An observer with four velocity $U^\mu$ measures the energy density

$$T_{\mu\nu}U^\mu U^\nu = -\frac{c^2\dot{m}}{4\pi r^2}\left(k_\mu U^\mu\right)^2,$$

which is nonnegative for every observer if and only if $\dot{m} \le 0$.
The weak energy condition therefore holds exactly when the star is losing mass, which is why the outgoing chart is the radiating one and the entry's parameter description says $\dot m \le 0$.

The luminosity follows from the same expression.
A static observer far from the star has $U^\mu = (1/\sqrt{f},0,0,0)$ in the chart, from $g_{\mu\nu}U^\mu U^\nu = -1$, so $k_\mu U^\mu = -1/\sqrt{f}$ and the energy density that observer measures is

$$\rho = -\frac{c^2\dot{m}}{4\pi r^2 f} \longrightarrow -\frac{c^2\dot{m}}{4\pi r^2} \quad (r \to \infty).$$

The radiation crosses the sphere of radius $r$ at the speed of light, so the power crossing it is $4\pi r^2c\rho$, and in the limit

$$L = -c^3\dot{m} = -c^2\frac{dm}{du}.$$

That is the reading Lindquist, Schwartz and Misner gave the mass function, quoted in the entry's history: the rate at which $m$ falls with retarded time is the luminosity, in units where an energy is $c^2$ times a mass.

---

## Step 11. The Weyl tensor

In four dimensions, with $R = 0$ here,

$$C_{\mu\nu\rho\sigma} = R_{\mu\nu\rho\sigma} - \tfrac{1}{2}\left(g_{\mu\rho}R_{\sigma\nu} - g_{\mu\sigma}R_{\rho\nu} - g_{\nu\rho}R_{\sigma\mu} + g_{\nu\sigma}R_{\rho\mu}\right).$$

The trace the formula wants is the published Ricci tensor itself, which is what the standard contraction buys.
Its only nonzero component is $R_{00} = -\dot{r}_s/r^2$, so a correction term survives only when one index of each pair is $0$ and the metric factor multiplying $R_{00}$ is nonzero.

For the $(0,r)$ plane every candidate correction vanishes:

$$C_{0r0r} = R_{0r0r} - \tfrac{1}{2}\left(g_{00}R_{rr} - g_{0r}R_{0r} - g_{r0}R_{r0} + g_{rr}R_{00}\right) = R_{0r0r} = -\frac{2Gm}{c^2r^3},$$

since $R_{rr} = R_{0r} = 0$ and $g_{rr} = 0$.
For the mixed plane the last term survives:

$$C_{0\theta 0\theta} = R_{0\theta 0\theta} - \tfrac{1}{2}g_{\theta\theta}R_{00} = R_{0\theta 0\theta} + \frac{\dot{r}_s}{2} = \frac{r_s}{2r}f = \frac{Gm}{c^2r}\left(1-\frac{2Gm}{c^2r}\right),$$

and the term in $\dot{m}$ has been removed exactly.
The remaining independent components, $C_{0\theta r\theta} = R_{0\theta r\theta}$ and $C_{\theta\phi\theta\phi} = R_{\theta\phi\theta\phi}$, are uncorrected, because $g_{\theta 0} = 0$ and $g_{\phi 0} = 0$ leave nothing for $R_{00}$ to multiply.

So the Weyl tensor is the Riemann tensor of Schwarzschild with $r_s$ read at the retarded time, and it carries no $\dot{m}$ anywhere.
The radiation lives entirely in the Ricci part of the curvature and the tidal part is Coulomb like, the same algebraically special form Schwarzschild has.
The entry's two Weyl variants are exactly its two Riemann variants with every $\dot{m}$ struck out, which is visible by eye in the published file.

---

## Step 12. The Kretschmann scalar

Write $P = R_{0r0r}$, $Q = R_{0\theta0\theta}$, $C = R_{0\theta r\theta}$ and $D = R_{\theta\phi\theta\phi}$ for the four independent components of Step 7, and recall $R_{r\theta r\theta} = 0$.

Group the sum $K = R_{\mu\nu\rho\sigma}R^{\mu\nu\rho\sigma}$ by index pairs.
Each unordered pair $\{\mu,\nu\}$ with $\mu \ne \nu$ contributes the same product in each of its two orderings, and the same holds for the second pair, so

$$K = 4\sum_{I,J}R_{IJ}R^{IJ},$$

the sum running over the six pairs $[0r], [0\theta], [0\phi], [r\theta], [r\phi], [\theta\phi]$ in each slot.
No published component mixes one plane with another, so only four blocks contribute.

**The $(0,r)$ plane.**
Raising an index $0$ replaces it by $-1$ times an index $r$, and raising an index $r$ gives $-1$ times an index $0$ plus $f$ times an index $r$.
Applied to the first pair, and using $R_{rr0r} = R_{000r} = 0$,

$$R^{0r}{}_{0r} = -g^{r\beta}R_{r\beta 0r} = R_{r00r} = -R_{0r0r},$$

and the same operation on the second pair gives $R^{0r0r} = R_{0r0r} = P$.
The block contributes $P^2 = r_s^2/r^6$.

**The $(0,\theta)$ and $(r,\theta)$ block.**
Collect the four components into the symmetric matrix $M = \begin{pmatrix} Q & C \\ C & 0\end{pmatrix}$ over the pair indices $[0\theta]$ and $[r\theta]$.
Raising such a pair index is the matrix $N = \dfrac{1}{r^2}\begin{pmatrix} 0 & -1 \\ -1 & f \end{pmatrix}$, the inverse metric of the $(0,r)$ block times the $g^{\theta\theta} = 1/r^2$ of the $\theta$ slot.
The block therefore contributes

$$\sum_{a,b}R_{a\theta b\theta}R^{a\theta b\theta} = \operatorname{tr}\left(MNMN\right) = \frac{2C^2}{r^4} = \frac{r_s^2}{2r^6},$$

where the trace is computed from $NM = \dfrac{1}{r^2}\begin{pmatrix} -C & 0 \\ fC - Q & -C\end{pmatrix}$, whose square has $C^2/r^4$ on both diagonal entries.

This is the step that explains Step 11's conclusion arithmetically.
$Q$ is the only component carrying $\dot{m}$, and it appears in $MNMN$ only off the diagonal, so it cannot reach the trace.
The radiation is invisible to this scalar not by cancellation but because the null structure never lets $Q$ multiply anything but itself, and $R_{r\theta r\theta} = 0$ is what it would have needed.

**The $(0,\phi)$ and $(r,\phi)$ block.**
Every entry carries one factor of $\sin^2\theta$ and $g^{\phi\phi}$ carries $1/\sin^2\theta$, so the factors cancel and the contribution is again $r_s^2/2r^6$.

**The $(\theta,\phi)$ plane.**

$$R_{\theta\phi\theta\phi}R^{\theta\phi\theta\phi} = D\cdot\left(g^{\theta\theta}\right)^2\left(g^{\phi\phi}\right)^2D = r_sr\sin^2\theta\cdot\frac{r_s}{r^7\sin^2\theta} = \frac{r_s^2}{r^6}.$$

Adding the four blocks and multiplying by the four orderings,

$$K = 4\left(\frac{r_s^2}{r^6} + \frac{r_s^2}{2r^6} + \frac{r_s^2}{2r^6} + \frac{r_s^2}{r^6}\right) = \frac{12r_s^2}{r^6} = \frac{48G^2m^2}{c^4r^6},$$

which is what the entry publishes.

Three things follow.
It is the Schwarzschild value with the instantaneous mass, so the curvature blows up as $r \to 0$ whenever the mass function has not already reached zero there, and that singularity is real rather than a defect of the chart.
It carries no $\dot m$, so no amount of radiation shows up in it; the same holds for every scalar built from the Ricci tensor, since $R = 0$ and $R_{\mu\nu}R^{\mu\nu} = R_{00}R^{00} = \left(g^{0r}\right)^2R_{00}R_{rr} = 0$, the Ricci tensor being null.
And it is finite at $r = r_s$, so the surface where $f$ vanishes is not a singularity here any more than it is in Schwarzschild.

---

## Step 13. The geodesic equations

The geodesic equation is $\ddot{x}^\mu + \Gamma^\mu{}_{\nu\rho}\dot{x}^\nu\dot{x}^\rho = 0$ with the dots denoting derivatives of the chart coordinates with respect to an affine parameter, so $\dot{u}$ is $d(cu)/d\lambda$ and the symbols are the chart symbols of Step 5.
Summing over the pairs, and doubling the symbols that are symmetric in two different indices,

$$\ddot{u} - \frac{Gm}{c^2r^2}\dot{u}^2 + r\dot{\theta}^2 + r\sin^2\theta\,\dot{\phi}^2 = 0,$$

$$\ddot{r} + \left(\frac{Gm}{c^2r^2}\left(1-\frac{2Gm}{c^2r}\right) - \frac{G\dot{m}}{c^2r}\right)\dot{u}^2 + \frac{2Gm}{c^2r^2}\dot{u}\dot{r} - r\left(1-\frac{2Gm}{c^2r}\right)\left(\dot{\theta}^2 + \sin^2\theta\,\dot{\phi}^2\right) = 0,$$

$$\ddot{\theta} + \frac{2}{r}\dot{r}\dot{\theta} - \cos\theta\sin\theta\,\dot{\phi}^2 = 0,$$

$$\ddot{\phi} + \frac{2}{r}\dot{r}\dot{\phi} + 2\cot\theta\,\dot{\theta}\dot{\phi} = 0.$$

The first equation comes from $\Gamma^0_{00}, \Gamma^0_{\theta\theta}, \Gamma^0_{\phi\phi}$; the second from $\Gamma^r_{00}$, the pair $\Gamma^r_{0r} = \Gamma^r_{r0}$, and $\Gamma^r_{\theta\theta}, \Gamma^r_{\phi\phi}$; the third from the pair $\Gamma^\theta_{r\theta}$ and from $\Gamma^\theta_{\phi\phi}$; the fourth from the pairs $\Gamma^\phi_{r\phi}$ and $\Gamma^\phi_{\theta\phi}$.

The radiation enters the motion of a test body through exactly one term, the $-\dfrac{G\dot{m}}{c^2r}\dot{u}^2$ of the radial equation.
It is proportional to $\dot{u}^2$, so it acts only on a body crossing the outgoing light fronts, and it flips sign with $\dot{m}$, so a radiating star and an accreting one bend the radial motion opposite ways.
Set $\dot{m} = 0$ and every equation here is the Schwarzschild one in the same chart.

The last two equations are the ordinary spherical ones, so the equatorial plane $\theta = \pi/2$ is preserved exactly as in Schwarzschild.

---

## Step 14. The radial null geodesics and the horizon

Setting $d\theta = d\phi = 0$ and $ds^2 = 0$ in the line element factorises it,

$$(c\,du)\left[f\,(c\,du) + 2\,dr\right] = 0,$$

so there are two families.
The first is $u = \text{const}$, whose tangent is exactly the $k^\mu = \partial_r$ of Step 10: these are the outgoing rays the radiation itself travels along, and Step 10 showed $r$ is an affine parameter on them.
The second obeys

$$\frac{d(cu)}{dr} = -\frac{2}{f},$$

the ingoing family, which is the one the chart bends in order to make the outgoing rays straight.

The expansion of a spherically symmetric null congruence is $\vartheta = \dfrac{2}{r}\dfrac{dr}{d\lambda}$, since the cross sections are spheres of area $4\pi r^2$.
The outgoing family has $dr/d\lambda = 1$ and therefore $\vartheta_{\text{out}} = 2/r$, which never vanishes.
The ingoing family, parametrised so that $d(cu)/d\lambda = 1$ and hence $dr/d\lambda = -f/2$, has

$$\vartheta_{\text{in}} = -\frac{f}{r} = -\frac{1}{r}\left(1 - \frac{2Gm(u)}{c^2r}\right),$$

which vanishes at $r = 2Gm(u)/c^2$ and is positive inside it.
Inside that radius both families are expanding, so the spheres there are anti trapped: the retarded solution describes a white hole exterior, or a star shining outward, rather than a black hole.
Step 16 shows the mirror statement for the advanced chart, where the same surface is the apparent horizon of a black hole.

The surface moves.
It sits at $r = 2Gm(u)/c^2$, which changes with retarded time whenever $\dot{m} \ne 0$, so it is a horizon that recedes as the star radiates rather than a null surface fixed once and for all.

---

## Step 15. Every published equation is dimensionally consistent

The declarations are $[r] = L$, $[u] = T$, $[\theta] = [\phi] = 1$, $[c] = LT^{-1}$, $[G] = L^3M^{-1}T^{-2}$ and $[m] = M$.
The chart coordinate is $x^0 = cu$, which carries $L$ like the other three, and

$$\left[\frac{Gm}{c^2}\right] = \frac{L^3M^{-1}T^{-2}\cdot M}{L^2T^{-2}} = L, \qquad [\dot{m}] = \frac{M}{L}, \qquad \left[\frac{G\dot{m}}{c^2}\right] = 1.$$

So $r_s$ is a length and $f$ and $\dot{r}_s$ are dimensionless, as Step 2 said.

The line element: $c^2du^2$ carries $L^2$, the cross term $c\,du\,dr$ carries $LT^{-1}\cdot T\cdot L = L^2$, and $r^2d\theta^2$ carries $L^2$.

The metric components: with every chart coordinate a length, $g_{\mu\nu}$ is dimensionless for the first two indices and carries $L^2$ against $\theta$ or $\phi$.
$g_{uu} = -f$ and $g_{ur} = -1$ are dimensionless, $g_{\theta\theta} = r^2$ carries $L^2$.

The Christoffel symbols carry $1/L$ with every index a length: $\Gamma^0_{00} = -Gm/(c^2r^2)$ carries $L/L^2 = 1/L$, and both terms of $\Gamma^r_{00}$ do too, the second being $G\dot{m}/(c^2r)$, which is $1/L$ because $G\dot m/c^2$ is dimensionless.
$\Gamma^0_{\theta\theta} = r$ carries $L$, which is what an index pattern with two angles and one length asks for, $1/L \cdot L^2 = L$.

The curvature carries $1/L^2$ on all length indices: $R_{0r0r} = -2Gm/(c^2r^3)$ is $L/L^3$, and both terms of $R_{0\theta0\theta}$ carry $L^0$, namely $(Gm/c^2)/r$ and $G\dot{m}/c^2$, which is $1/L^2$ times the $L^2$ two angular indices ask for.
The Ricci and Einstein components carry $1/L^2$: $2G\dot{m}/(c^2r^2)$ is a dimensionless numerator over $L^2$.
The Kretschmann scalar carries $1/L^4$: $48G^2m^2/(c^4r^6)$ is $L^2/L^6$.

The geodesic equations are measured against their own second derivative, so every term carries $L/\lambda^2$ with $\lambda$ the affine parameter.
In the radial equation $\ddot r$ carries $L/\lambda^2$; the term $\dfrac{Gm}{c^2r^2}\dot{u}^2$ carries $\dfrac{L}{L^2}\cdot\dfrac{L^2}{\lambda^2} = \dfrac{L}{\lambda^2}$, using that $\dot u = d(cu)/d\lambda$ carries $L/\lambda$ and not $T/\lambda$; the term $\dfrac{G\dot{m}}{c^2r}\dot{u}^2$ carries $\dfrac{1}{L}\cdot\dfrac{L^2}{\lambda^2}$, the same; and $r\dot\theta^2$ carries $L/\lambda^2$ because an angle is dimensionless.
Read the other way, with $\dot{u} = du/d\lambda$ carrying $T/\lambda$, the first of those terms would come out short by a factor of $c^2$ against $\ddot{r}$, which is the failure the checker's dimensional pass is built to catch.

---

## Step 16. The ingoing chart

The advanced, or ingoing, form is

$$ds^2 = -\left(1 - \frac{2Gm(v)}{c^2r}\right)c^2dv^2 + 2c\,dv\,dr + r^2d\theta^2 + r^2\sin^2\theta\,d\phi^2,$$

with $v$ the advanced time and $m$ now a function of it.
It differs from Step 2 in the sign of the cross term alone, which means it is the same computation under the reflection

$$cu = -cv, \qquad \text{so} \qquad \frac{\partial}{\partial(cu)} = -\frac{\partial}{\partial(cv)}, \qquad \dot{m} \longmapsto -\dot{m}.$$

Under a reflection of one chart coordinate a tensor component picks up a factor of $-1$ for each index in that slot, upper or lower, and the Christoffel symbols follow the same rule because the transformation is linear.
So every published value of the ingoing chart is the outgoing one with

$$(-1)^{\#\text{ of } 0 \text{ indices}} \quad\text{and}\quad \dot{m} \to -\dot{m}$$

applied together.
Three checks against the published file:

$\Gamma^v_{vv} = +\dfrac{Gm}{c^2r^2}$, three indices in the reflected slot and no $\dot m$, so the sign flips from the outgoing $-\dfrac{Gm}{c^2r^2}$.

$\Gamma_{vvv} = \dfrac{G\dot{m}}{c^2r}$, three indices and one $\dot{m}$, so the two sign flips cancel and the value is unchanged.

$R_{vv} = +\dfrac{2G\dot{m}}{c^2r^2}$, two indices and one $\dot{m}$, so the value flips from the outgoing $-\dfrac{2G\dot{m}}{c^2r^2}$.

The physics the reflection carries is the interesting part.
The stress energy becomes

$$T_{\mu\nu} = +\frac{c^2\dot{m}}{4\pi r^2}k_\mu k_\nu, \qquad k_\mu dx^\mu = -d(cv), \qquad k^\mu\partial_\mu = -\partial_r,$$

so the null rays now run inward and the weak energy condition holds when $\dot{m} \ge 0$: the advanced chart is the accreting one, and the entry's parameter description says so.
Redoing Step 14 in it, the outgoing family obeys $d(cv)/dr = 2/f$ and, parametrised so that $dr/d\lambda = f$, has expansion

$$\vartheta_{\text{out}} = \frac{2f}{r},$$

which vanishes at $r = 2Gm(v)/c^2$ and is negative inside it, where the ingoing expansion $-2/r$ is negative too.
The spheres inside that radius are trapped, so the advanced solution is the black hole one, growing as it swallows the dust.
That is the chart in which collapse is studied, and it is the setting of the naked singularity result the entry's history cites.

Everything else is unchanged, including the Kretschmann scalar $K = 48G^2m^2/(c^4r^6)$, which is even in the number of reflected indices and free of $\dot{m}$.

---

## Step 17. What the entry publishes, and what the checker needs

For each of the two charts the entry publishes the line element, the metric and its inverse, both Christoffel variants, both Riemann variants, the three Ricci variants, the Ricci scalar, the Kretschmann scalar, the three Einstein variants, both Weyl variants and the four geodesic equations.

`verify_metrics.py` needed one declaration per chart for it.
`DIMENSIONS` gained an entry for each, declaring $u$ and $v$ as times, so the checker knows the chart multiplies them by $c$, along with $[G] = L^3M^{-1}T^{-2}$ and $[m] = M$.
That is the first mass in the table, because every other entry in the collection folds its mass into a length such as $r_s$ and never has to name one, so `BASE_DIMENSIONS` gained $M$ alongside $L$ and $T$.

The mass function is declared to the reader as `m = m(u)`, which is the same declaration FRW makes for its scale factor, and it is what tells the checker that $\dot{m}$ means a derivative with respect to the chart coordinate $cu$ rather than with respect to $u$.
No entry in `PARAMETER_RELATIONS` is needed, because $G$ and $m$ are free: nothing in this solution constrains the mass function, and that is exactly why it models a star whose luminosity is whatever its physics makes it.
