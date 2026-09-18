# The Krasnikov tube

This is the working behind the coordinate system in `MFS/assets/data/metrics/krasnikov.json`.
Every value the entry prints is derived here, in order, from the line element down to the geodesic equations.
Nothing is left as an exercise and nothing is asserted that is not computed.

The companion script `_tools/derivations/verify_metrics.py` does the same work in sympy and compares it against the published file, so the algebra below is checkable by hand and by machine independently.

The entry exists for one contraction.
The tube is built to tip one null ray, the one that comes home along the axis, and Step 11 is where the mathematics hands back the bill for tipping it.
Writing $w = \ln\left(1 + k\right)$, the Einstein tensor contracted twice with that ray is

$$G_{\mu\nu}\ell_{\text{back}}^\mu\ell_{\text{back}}^\nu = -\frac{\left(1 + k\right)^2}{2r}\,\partial_r\left(r\,\partial_r w\right),$$

minus the flat Laplacian of $w$ in the transverse plane, against a weight that is positive wherever the metric is not degenerate.
A Laplacian integrates to zero over a plane on which its argument is regular at the origin and constant far out, so the contraction cannot have one sign unless it vanishes identically, and it vanishes identically only when $k$ does not depend on $r$, which is to say only when there is no tube.
The null energy condition therefore fails in the wall of every Krasnikov tube of finite radius, for every shape function, and it fails on exactly the ray the construction exists to open.
Step 10 reads the same fact as an energy density and gets a sum of squares with a minus sign in front of it, which is the same shape the Alcubierre entry arrives at.

---

## Step 1. Conventions

The signature is $(-,+,+,+)$ and the Christoffel symbols are those of the Levi-Civita connection,

$$\Gamma^\mu{}_{\nu\sigma} = \tfrac{1}{2}g^{\mu\alpha}\left(\partial_\nu g_{\alpha\sigma} + \partial_\sigma g_{\alpha\nu} - \partial_\alpha g_{\nu\sigma}\right),$$

with the lowered symbol $\Gamma_{\mu\nu\sigma} = g_{\mu\alpha}\Gamma^\alpha{}_{\nu\sigma}$.
The Riemann tensor is

$$R^\mu{}_{\nu\sigma\tau} = \partial_\sigma\Gamma^\mu{}_{\nu\tau} - \partial_\tau\Gamma^\mu{}_{\nu\sigma} + \Gamma^\mu{}_{\sigma\lambda}\Gamma^\lambda{}_{\nu\tau} - \Gamma^\mu{}_{\tau\lambda}\Gamma^\lambda{}_{\nu\sigma},$$

which is what the published Riemann components are in.
The dummy indices are written $\sigma$ and $\tau$ rather than $\rho$ and $\sigma$ throughout, because $r$ is a coordinate here and $\rho$ is reserved for the energy density.

The Ricci tensor is the standard contraction,

$$R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu},$$

settled for the collection on 2026-09-18.
The settlement is not free here.
The wall of the tube is not a vacuum, so the Ricci tensor, the Ricci scalar and the Einstein tensor all change sign with the choice, and the sign the collection now carries is the one that puts ordinary matter at positive energy density.
The entire question this entry answers is the sign of an energy density, so the answer would read backwards under the other contraction.

The chart is the collection's, the one whose zeroth coordinate is

$$x^0 = ct.$$

The index is printed with the bare letter $t$, but the component printed against it is a component of that chart, and everything from Step 4 onward is computed directly in it.
A printed $\partial_t$ is the derivative along that chart coordinate,

$$\partial_t = \frac{1}{c}\frac{\partial}{\partial t},$$

so each printed derivative carries one inverse length per order, which is what Step 16 checks.
The dots in the geodesic equations are velocities of that same chart, so $\dot{t}$ means $d(ct)/d\lambda$.

---

## Step 2. What Krasnikov built, and why the entry is four dimensional

Krasnikov's problem is not how to go fast.
It is that the distance to Deneb is fixed by the geometry, and the half of his 1998 paper that fails proves you cannot shorten the outbound leg from inside a globally hyperbolic spacetime.
What is left is the way back.
Nothing forbids modifying the metric in a region the ship has already passed through, and a modification that swings the light cones far enough open lets the return leg run backward through the clock left at home.

He wrote that modification in two dimensions:

$$ds^2 = -\left(c\,dt - dx\right)\left(c\,dt + k\,dx\right),$$

with $k$ a function of $t$ and $x$ built out of smoothed step functions, equal to $1$ outside the modified region and to a constant just above $-1$ inside it.
The two null directions of this line element are read straight off the two brackets.
Setting the first to zero gives $d(ct) = dx$, the outbound ray, which does not move no matter what $k$ does.
Setting the second to zero gives $d(ct) = -k\,dx$, and that is the whole of the construction: at $k = 1$ it is the ordinary homeward ray, at $k = 0$ it is the spatial direction $-x$, and at $k < 0$ it runs toward decreasing $x$ with $d(ct)$ negative as well.

A two dimensional entry would stop there, and stopping there is exactly what this collection cannot do.
In two dimensions the Einstein tensor vanishes identically for every metric there is, because $R_{\mu\nu} = \tfrac{1}{2}Rg_{\mu\nu}$ is an identity when the Riemann tensor has one independent component.
So a two dimensional Krasnikov entry could print a line element, a connection, a curvature scalar and a set of tipped cones, and it could print nothing at all about the source, because in two dimensions there is no source to print.
Since the source is the whole of the physics, and since the history of this spacetime is the history of an argument about how much exotic matter it needs, the entry carries the four dimensional form instead.

That form is Everett and Roman's, and it is the same product with a flat transverse plane appended,

$$ds^2 = -\left(c\,dt - dx\right)\left(c\,dt + k\,dx\right) + dr^2 + r^2d\phi^2,$$

with $k$ now a function of $t$, $x$ and $r$.
What the fourth dimension adds is not decoration.
In two dimensions there is no direction for $k$ to fall off in, so the modified region is a slab that fills all of space at the values of $x$ it covers, and it has no wall.
In four dimensions $k$ returns to $1$ beyond some radius, the modified region becomes a tube of finite cross section, and there is a wall where $k$ turns.
Every curvature component in this entry is supported on that wall, and every one of them vanishes when $k$ is independent of $r$, which is the two dimensional case dressed in four coordinates.

---

## Step 3. The line element, and the domain it lives on

The coordinates are $\left(t, x, r, \phi\right)$: $t$ carries dimensions of time, $x$ is length along the axis the tube is laid along, $r$ is distance from that axis and $\phi$ is the angle around it.
The one parameter is $k$, dimensionless, an arbitrary function of $t$, $x$ and $r$ and not of $\phi$.
Expanding the product,

$$ds^2 = -c^2dt^2 + \left(1 - k\right)c\,dt\,dx + k\,dx^2 + dr^2 + r^2d\phi^2,$$

which is the form Everett and Roman print.
At $k = 1$ the cross term goes and the $dx^2$ term becomes $+dx^2$, leaving Minkowski space in cylindrical coordinates, so the tube is asymptotically flat in every direction and is flat wherever $k$ is constant, whatever that constant is.

The published $k$ is left arbitrary.
Krasnikov assembles his own out of smoothed step functions: one closes the tube off at a radius, one holds it inside the future light cone of the departure so that nothing is ever built ahead of the ship, and a pair opens it between the start of the journey and its end.
Nothing below uses any of that, so every component published here is an identity for an arbitrary $k$.

One restriction is not optional.
The determinant computed in Step 4 is $\det g = -\tfrac{1}{4}r^2\left(1 + k\right)^2$, which vanishes at $k = -1$, and $1 + k$ divides every curvature component in the entry.
So the solution lives on $k > -1$, and Krasnikov's interior value $k = -1 + \delta$ approaches the boundary of the chart from above as $\delta$ goes to zero.
That is not a technicality, it is the mechanism.
The return leg runs along $d(ct) = -k\,dx$, so at $k = -1 + \delta$ a homeward leg of $\Delta x = -D$ takes $\Delta(ct) = -\left(1 - \delta\right)D$, and a round trip that went out at nearly the speed of light costs $\delta D$ in all.
Making the round trip shorter means making $\delta$ smaller, which means driving the metric toward degeneracy, and Step 12 is where the price of that shows up.
The trip time never reaches zero and never goes negative, which is why one tube is not a time machine and two of them, laid along different routes, are.

---

## Step 4. The metric matrix, its determinant and its inverse

Reading the coefficients off the line element in the chart $x^0 = ct$, where the metric carries no dimensions,

$$g_{\mu\nu} = \begin{pmatrix} -1 & \tfrac{1}{2}\left(1-k\right) & 0 & 0 \\ \tfrac{1}{2}\left(1-k\right) & k & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & r^2\end{pmatrix}.$$

Two things about this matrix are worth noticing before any derivative is taken.
The first is that $g_{tt} = -1$ exactly, for every $k$: the lapse of this slicing is one, and the observer sitting at fixed $x$, $r$ and $\phi$ has $t$ for a proper time.
The second is that $g_{xx} = k$ is negative inside the tube, so $x$ is a timelike coordinate there.
Neither of those makes the chart bad, and the second is the tipping of the cones stated as a component.

The block is two by two plus a flat transverse plane, so the determinant is

$$\det g = r^2\left[\left(-1\right)k - \tfrac{1}{4}\left(1-k\right)^2\right] = -\tfrac{1}{4}r^2\left(k^2 + 2k + 1\right) = -\tfrac{1}{4}r^2\left(1 + k\right)^2.$$

It is negative for every $k \ne -1$, as a Lorentzian metric's determinant must be in four dimensions, and it degenerates only at $k = -1$.
Inverting the two by two block with that determinant,

$$g^{\mu\nu} = \begin{pmatrix} -\dfrac{4k}{\left(1+k\right)^2} & \dfrac{2\left(1-k\right)}{\left(1+k\right)^2} & 0 & 0 \\[4pt] \dfrac{2\left(1-k\right)}{\left(1+k\right)^2} & \dfrac{4}{\left(1+k\right)^2} & 0 & 0 \\[4pt] 0 & 0 & 1 & 0 \\[4pt] 0 & 0 & 0 & \dfrac{1}{r^2}\end{pmatrix},$$

which at $k=1$ is $\operatorname{diag}\left(-1, 1, 1, 1/r^2\right)$ as it should be.
Every one of these six entries is published, and the factor $\left(1+k\right)^{-2}$ in three of them is where the powers of $1 + k$ in the whole entry come from.

---

## Step 5. The lowered Christoffel symbols

The lowered symbol is

$$\Gamma_{\mu\nu\sigma} = \tfrac{1}{2}\left(\partial_\nu g_{\mu\sigma} + \partial_\sigma g_{\mu\nu} - \partial_\mu g_{\nu\sigma}\right),$$

and the matrix of Step 4 makes most of these short, because only $g_{tx}$ and $g_{xx}$ depend on $k$ at all and $g_{tt}$ depends on nothing.

That last fact kills a whole family at once:

$$\Gamma_{tt\sigma} = \tfrac{1}{2}\left(\partial_t g_{t\sigma} + \partial_\sigma g_{tt} - \partial_t g_{t\sigma}\right) = \tfrac{1}{2}\partial_\sigma g_{tt} = 0,$$

so every symbol with two leading $t$ indices vanishes.
Four more, worked:

$$\Gamma_{txx} = \partial_x g_{tx} - \tfrac{1}{2}\partial_t g_{xx} = -\tfrac{1}{2}\partial_x k - \tfrac{1}{2}\partial_t k,$$

$$\Gamma_{xtt} = \partial_t g_{xt} - \tfrac{1}{2}\partial_x g_{tt} = -\tfrac{1}{2}\partial_t k,$$

$$\Gamma_{rtx} = \tfrac{1}{2}\left(\partial_t g_{rx} + \partial_x g_{rt} - \partial_r g_{tx}\right) = -\tfrac{1}{2}\partial_r\!\left(\tfrac{1}{2}\left(1-k\right)\right) = \tfrac{1}{4}\partial_r k,$$

$$\Gamma_{rxx} = \partial_x g_{rx} - \tfrac{1}{2}\partial_r g_{xx} = -\tfrac{1}{2}\partial_r k.$$

The complete list of the seventeen nonzero lowered symbols is

$$\Gamma_{txx} = -\tfrac{1}{2}\left(\partial_t k + \partial_x k\right),\qquad \Gamma_{txr} = \Gamma_{trx} = -\tfrac{1}{4}\partial_r k,$$
$$\Gamma_{xtt} = -\tfrac{1}{2}\partial_t k,\qquad \Gamma_{xtx} = \Gamma_{xxt} = \tfrac{1}{2}\partial_t k,\qquad \Gamma_{xtr} = \Gamma_{xrt} = -\tfrac{1}{4}\partial_r k,$$
$$\Gamma_{xxx} = \tfrac{1}{2}\partial_x k,\qquad \Gamma_{xxr} = \Gamma_{xrx} = \tfrac{1}{2}\partial_r k,$$
$$\Gamma_{rtx} = \Gamma_{rxt} = \tfrac{1}{4}\partial_r k,\qquad \Gamma_{rxx} = -\tfrac{1}{2}\partial_r k,\qquad \Gamma_{r\phi\phi} = -r,$$
$$\Gamma_{\phi r\phi} = \Gamma_{\phi\phi r} = r.$$

Every one of them is a single first derivative of $k$ against a rational number, with no $k$ in front of it anywhere.
The last three carry no $k$ at all: they are the symbols of the flat transverse plane in polar coordinates, and they would be there in Minkowski space written the same way.

---

## Step 6. The Christoffel symbols with an upper index

Raising the first index with $g^{\mu\alpha}$ of Step 4 mixes the $t$ and $x$ rows only, since $g^{tr}$, $g^{t\phi}$, $g^{xr}$ and $g^{x\phi}$ all vanish.
The $r$ and $\phi$ rows are raised by $g^{rr}=1$ and $g^{\phi\phi} = r^{-2}$, so they barely change:

$$\Gamma^r{}_{tx} = \Gamma_{rtx} = \tfrac{1}{4}\partial_r k,\qquad \Gamma^r{}_{xx} = -\tfrac{1}{2}\partial_r k,\qquad \Gamma^r{}_{\phi\phi} = -r,\qquad \Gamma^\phi{}_{r\phi} = \frac{1}{r^2}\cdot r = \frac{1}{r}.$$

The first two of these are the only place a first derivative of $k$ survives with no $1+k$ underneath it, and they are the reason the radial geodesic equation of Step 15 is as short as it is.

The $t$ and $x$ rows pick up the inverse metric in full.
Two worked:

$$\Gamma^t{}_{tt} = g^{tt}\Gamma_{ttt} + g^{tx}\Gamma_{xtt} = 0 + \frac{2\left(1-k\right)}{\left(1+k\right)^2}\left(-\tfrac{1}{2}\partial_t k\right) = \frac{\left(k-1\right)\partial_t k}{\left(1+k\right)^2},$$

$$\Gamma^x{}_{tt} = g^{xt}\Gamma_{ttt} + g^{xx}\Gamma_{xtt} = 0 + \frac{4}{\left(1+k\right)^2}\left(-\tfrac{1}{2}\partial_t k\right) = -\frac{2\,\partial_t k}{\left(1+k\right)^2}.$$

The first is worth a second look.
$\Gamma^t{}_{tt}$ is not zero, so coordinate time is not an affine parameter along the world lines of constant $x$, $r$ and $\phi$ once the tube is switched on, even though those world lines have $t$ for a proper time.
That is not a contradiction: they are not geodesics.
The observer who holds station at fixed $x$ inside the tube is accelerating, and the acceleration is what $\partial_t k$ measures.

The whole upper index list is the twenty two symbols

$$\Gamma^t{}_{tt} = \frac{\left(k-1\right)\partial_t k}{\left(1+k\right)^2},\qquad \Gamma^t{}_{tx} = -\frac{\left(k-1\right)\partial_t k}{\left(1+k\right)^2},\qquad \Gamma^t{}_{tr} = \frac{\left(k-1\right)\partial_r k}{2\left(1+k\right)^2},$$
$$\Gamma^t{}_{xx} = \frac{2k\,\partial_t k + \left(1+k\right)\partial_x k}{\left(1+k\right)^2},\qquad \Gamma^t{}_{xr} = \frac{\partial_r k}{\left(1+k\right)^2},$$
$$\Gamma^x{}_{tt} = -\frac{2\,\partial_t k}{\left(1+k\right)^2},\qquad \Gamma^x{}_{tx} = \frac{2\,\partial_t k}{\left(1+k\right)^2},\qquad \Gamma^x{}_{tr} = -\frac{\partial_r k}{\left(1+k\right)^2},$$
$$\Gamma^x{}_{xx} = \frac{\left(1+k\right)\left(\partial_t k + \partial_x k\right) - 2\,\partial_t k}{\left(1+k\right)^2},\qquad \Gamma^x{}_{xr} = \frac{\left(k+3\right)\partial_r k}{2\left(1+k\right)^2},$$
$$\Gamma^r{}_{tx} = \tfrac{1}{4}\partial_r k,\qquad \Gamma^r{}_{xx} = -\tfrac{1}{2}\partial_r k,\qquad \Gamma^r{}_{\phi\phi} = -r,\qquad \Gamma^\phi{}_{r\phi} = \frac{1}{r},$$

with the symmetry in the last two indices supplying the rest.
Every symbol that involves $k$ vanishes where $k$ is constant, which is everywhere except the wall: the inside of the tube and the space outside it are both flat, and the tube is a seam between two flat regions.

---

## Step 7. The Riemann tensor

With the connection in hand the curvature is the definition of Step 1 applied component by component.
Eight blocks come out nonzero, and the shortest of them is the one that does the most work.

Take $R_{trtr}$, in the plane spanned by the time and the radius.
Working with the upper index first,

$$R^t{}_{rtr} = \partial_t\Gamma^t{}_{rr} - \partial_r\Gamma^t{}_{rt} + \Gamma^t{}_{t\lambda}\Gamma^\lambda{}_{rr} - \Gamma^t{}_{r\lambda}\Gamma^\lambda{}_{rt}.$$

Here $\Gamma^t{}_{rr} = 0$ and $\Gamma^\lambda{}_{rr} = 0$ for every $\lambda$, so the first and third terms are absent, and in the fourth only $\lambda = t$ and $\lambda = x$ survive.
What is left carries a $\partial_r^2 k$ out of the surviving derivative term.
Lowering the first index has to go through both $g_{tt}$ and $g_{tx}$, since $g^{tx}$ is not zero, and the second derivative cancels between the two contributions, leaving

$$R_{trtr} = \frac{\left(\partial_r k\right)^2}{4\left(1+k\right)^2} = \tfrac{1}{4}\left(\partial_r\ln\left(1+k\right)\right)^2.$$

A perfect square, with the whole $k$ dependence packaged into the logarithm.
That is the first appearance of $w = \ln\left(1+k\right)$, and from here on it is the variable the geometry is actually written in.

The rest of the independent blocks, in the same variable where they shorten:

$$R_{txtx} = \frac{\left(1+k\right)\left(\partial_r k\right)^2 - 8\left(1+k\right)\left(\partial_t^2 k + \partial_t\partial_x k\right) + 8\,\partial_t k\left(\partial_t k + \partial_x k\right)}{16\left(1+k\right)},$$

$$R_{txtr} = -\frac{\left(1+k\right)\partial_t\partial_r k - \partial_t k\,\partial_r k}{4\left(1+k\right)},\qquad R_{t\phi x\phi} = \tfrac{1}{4}r\,\partial_r k,$$

$$R_{trxr} = \frac{2\left(1+k\right)^2\partial_r^2 k - \left(1+k\right)\left(\partial_r k\right)^2 - 2\left(\partial_r k\right)^2}{8\left(1+k\right)^2},$$

together with $R_{txxr}$, $R_{xrxr}$ and $R_{x\phi x\phi}$, all printed in the entry.
Four pairings that a metric of this shape might have been expected to curve vanish identically for every $k$: $R_{t\phi t\phi}$, $R_{t\phi r\phi}$, $R_{x\phi r\phi}$ and $R_{r\phi r\phi}$ are all zero.
They will not stay zero in Step 13, and that is the point of Step 13.

The published tables carry sixty two nonzero components with an upper index and forty eight with all four down, each of them this computation with different indices.
Every one of them vanishes when $k$ is constant, which says again that the tube is flat inside and out.

---

## Step 8. The Ricci tensor and the Ricci scalar

Contracting on the first lower index as Step 1 says, the two transverse components are the short ones:

$$R_{rr} = -\frac{2\left(1+k\right)\partial_r^2 k - \left(\partial_r k\right)^2}{2\left(1+k\right)^2} = -\partial_r^2 w - \tfrac{1}{2}\left(\partial_r w\right)^2,\qquad R_{\phi\phi} = -\frac{r\,\partial_r k}{1+k} = -r\,\partial_r w,$$

using $\partial_r w = \partial_r k/\left(1+k\right)$ and $\partial_r^2 w = \partial_r^2 k/\left(1+k\right) - \left(\partial_r k\right)^2/\left(1+k\right)^2$ to go between the two forms.
The $\phi\phi$ component is the cleanest statement in the entry that curvature here is the radial gradient of $w$ and nothing else.

The time component is

$$R_{tt} = \frac{\left(1+k\right)\left(\partial_r k\right)^2 - 4\left(1+k\right)\left(\partial_t^2 k + \partial_t\partial_x k\right) + 4\,\partial_t k\left(\partial_t k + \partial_x k\right)}{2\left(1+k\right)^3},$$

with $R_{tx}$, $R_{tr}$, $R_{xx}$ and $R_{xr}$ longer and printed in the entry, in three variants.
The combinations $\partial_t k + \partial_x k$ and $\partial_t^2 k + \partial_t\partial_x k$ appear together everywhere in this tensor, and they are not an accident: $\partial_t + \partial_x$ is the derivative along the outbound null ray $\ell_{\text{out}}$ of Step 2, and the tensor is organised by it.

Contracting once more with $g^{\mu\nu}$,

$$R = -\frac{4r\left(1+k\right)^2\partial_r^2 k + 4\left(1+k\right)^2\partial_r k - r\left(1+k\right)\left(\partial_r k\right)^2 - 8r\left(1+k\right)\left(\partial_t^2 k + \partial_t\partial_x k\right) + 8r\,\partial_t k\left(\partial_t k + \partial_x k\right)}{2r\left(1+k\right)^3},$$

which is the published scalar.
It is not sign definite and it is not the thing to look at.
Its longitudinal terms can be made either sign by choosing how sharply the tube switches on in $t$ and $x$, which is why the objection to the tube, like the objection to the warp drive, has to be stated as a contraction of the Einstein tensor with a vector rather than as a statement about $R$.

---

## Step 9. The Einstein tensor

Subtracting the trace,

$$G_{\mu\nu} = R_{\mu\nu} - \tfrac{1}{2}Rg_{\mu\nu},$$

gives ten nonzero components, published in all three variants.
One of them is the entry, and it is worth writing out before the others:

$$G_{tt} = -\frac{4r\left(1+k\right)\partial_r^2 k + 4\left(1+k\right)\partial_r k - 3r\left(\partial_r k\right)^2}{4r\left(1+k\right)^2}.$$

Two things have happened to it that did not happen to $R_{tt}$.
Every derivative along $t$ and along $x$ has cancelled out of it, so $G_{tt}$ depends on the radial profile of $k$ alone, at each point, however violently the tube is switching on there.
And what is left reorganises completely in $w$.
Using $\partial_r k = \left(1+k\right)\partial_r w$ and $\partial_r^2 k = \left(1+k\right)\left(\partial_r^2 w + \left(\partial_r w\right)^2\right)$,

$$G_{tt} = -\left(\partial_r^2 w + \left(\partial_r w\right)^2\right) - \frac{\partial_r w}{r} + \tfrac{3}{4}\left(\partial_r w\right)^2 = -\frac{1}{r}\partial_r\left(r\,\partial_r w\right) - \tfrac{1}{4}\left(\partial_r w\right)^2.$$

The first term is the flat two dimensional Laplacian of $w$ in the transverse plane and the second is minus a square.
The rest of the tensor is longer and printed in the entry: $G_{tx}$, $G_{tr}$, $G_{xx}$, $G_{xr}$, $G_{rr}$ and $G_{\phi\phi}$.
Of those, $G_{tx}$ and $G_{xx}$ are free of derivatives along $t$ and along $x$ as $G_{tt}$ is, and the other four are not.

---

## Step 10. The energy density, and the sign it cannot escape

The observer to read $G_{tt}$ with is the one already picked out by Step 4.
Since $g_{tt} = -1$ in this chart for every $k$, the vector

$$u^\mu = \left(1, 0, 0, 0\right),\qquad g_{\mu\nu}u^\mu u^\nu = -1,$$

is a unit timelike vector everywhere on $k > -1$, and its world lines are the observers who hold station at fixed $x$, $r$ and $\phi$ with $t$ for a proper time.
Far outside the tube they are ordinary inertial observers of the flat space the tube is laid through.
They are future directed inside it too, which takes one line to check: writing $u$ in the basis of the two null directions of Step 2,

$$\left(1,0\right) = \frac{1}{1+k}\left(1,1\right) + \frac{1}{1+k}\left(k,-1\right),$$

and both coefficients are positive exactly when $k > -1$.

Through the field equations the energy density this observer measures is

$$\rho = T_{\mu\nu}u^\mu u^\nu = \frac{c^4}{8\pi G}G_{\mu\nu}u^\mu u^\nu = \frac{c^4}{8\pi G}G_{tt} = -\frac{c^4}{8\pi G}\left[\frac{1}{r}\partial_r\left(r\,\partial_r w\right) + \tfrac{1}{4}\left(\partial_r w\right)^2\right].$$

The second term is negative everywhere.
The first has both signs, and it integrates to nothing.
The surfaces of constant $t$ and $x$ carry the metric $dr^2 + r^2d\phi^2$, which is the Euclidean plane, so the area element is $r\,dr\,d\phi$ and

$$\int_0^\infty \frac{1}{r}\partial_r\left(r\,\partial_r w\right)r\,dr = \left[r\,\partial_r w\right]_0^\infty = 0$$

for any $k$ that is regular on the axis and constant outside the tube.
The energy per unit length of the tube, integrated over that plane, is therefore

$$\int\rho\,dA = -\frac{c^4}{16 G}\int_0^\infty\left(\partial_r w\right)^2 r\,dr,$$

minus a sum of squares, which is the same shape the Alcubierre entry arrives at by a different route.
It is strictly negative for every $k$ that depends on $r$ at all, and no choice of shape function and no choice of tube depth can turn it positive.
A tube that dodged the negative energy would have to be one with $\partial_r w = 0$ everywhere, and that is a tube of infinite radius, which is to say the two dimensional case of Step 2 and not a tube.

---

## Step 11. The null energy condition, and the ray it fails on

The density of Step 10 is one observer's reading and it is a weak energy condition statement.
The sharper statement is the null one, and it is sharper here because it can be aimed.

The two null directions of the $\left(t,x\right)$ plane, as vectors in this chart, are

$$\ell_{\text{out}}^\mu = \left(1, 1, 0, 0\right),\qquad \ell_{\text{back}}^\mu = \left(k, -1, 0, 0\right),$$

and both are null for every $k$, which the matrix of Step 4 confirms in two lines:

$$g_{\mu\nu}\ell_{\text{out}}^\mu\ell_{\text{out}}^\nu = -1 + \left(1-k\right) + k = 0,$$

$$g_{\mu\nu}\ell_{\text{back}}^\mu\ell_{\text{back}}^\nu = -k^2 - k\left(1-k\right) + k = 0.$$

The first is the outbound ray, which no tube ever moves.
The second is the ray the tube was built to tip, the one that comes home.

Contracting the Einstein tensor of Step 9 with each of them gives the two results this entry exists for.
On the outbound ray,

$$G_{\mu\nu}\ell_{\text{out}}^\mu\ell_{\text{out}}^\nu = G_{tt} + 2G_{tx} + G_{xx} = 0$$

identically, for every $k$, with no assumption about the shape function whatever.
The tube asks nothing of the source along the direction it did not alter.

On the homeward ray,

$$G_{\mu\nu}\ell_{\text{back}}^\mu\ell_{\text{back}}^\nu = k^2G_{tt} - 2kG_{tx} + G_{xx} = -\frac{\left(1+k\right)^2}{2r}\,\partial_r\left(r\,\partial_r w\right).$$

Every derivative along $t$ and along $x$ has cancelled here too, and what is left is the transverse Laplacian of $w$ again, now against the weight $\left(1+k\right)^2/2$, which is positive wherever the chart is good.

That settles it.
Suppose the null energy condition held everywhere on some slice, so that $G_{\mu\nu}\ell_{\text{back}}^\mu\ell_{\text{back}}^\nu \ge 0$ at every point of it.
Then $\partial_r\left(r\,\partial_r w\right) \le 0$ for every $r$.
But its integral over $r$ is $\left[r\,\partial_r w\right]_0^\infty$, which is zero for the same reason as in Step 10, and a function of one sign with zero integral is zero.
So $r\,\partial_r w$ is constant, and regularity on the axis makes the constant zero, so $\partial_r w = 0$, so $k$ does not depend on $r$, so there is no tube.

The null energy condition fails in the wall of every Krasnikov tube of finite radius, for every shape function, at some radius on every slice.
This is Visser, Bassett and Liberati's superluminal censorship for this particular geometry, written as an identity rather than as a theorem with hypotheses: tipping a light cone in Einstein gravity is the same act as putting negative energy along the direction you tipped it.

---

## Step 12. How much, and why a fast tube costs more

Step 10 gives the amount as well as the sign.
Take a slice through the middle of a tube of radius $r_0$ whose wall has thickness $\Delta$, with $k = -1+\delta$ inside and $k = 1$ outside.
Then $w = \ln\left(1+k\right)$ runs from $\ln\delta$ to $\ln 2$ across the wall, so $\partial_r w$ is of order $\ln\left(2/\delta\right)/\Delta$ there, and

$$\int\rho\,dA \sim -\frac{c^4}{16G}\,r_0\Delta\left(\frac{\ln\left(2/\delta\right)}{\Delta}\right)^2 = -\frac{c^4}{16G}\,\frac{r_0}{\Delta}\left(\ln\frac{2}{\delta}\right)^2,$$

since the integrand is supported on an annulus of radius $r_0$ and width $\Delta$.

Three readings come out of that one expression.
A thinner wall costs more, as $1/\Delta$.
A wider tube costs more, as $r_0$.
And a faster round trip costs more, as the square of the logarithm of $1/\delta$, because $\delta$ is at once the fraction of the distance the round trip costs in time and the distance the metric sits from degenerate.

The factor out front is $c^4/16G$, about $7.6\times10^{42}$ newtons, and it is a constant, so what sets the size of the answer is the ratio $r_0/\Delta$ of the radius of the tube to the thickness of its wall.
Everett and Roman's numbers come from that ratio: run it against the quantum inequalities, which cap how thin a wall of negative energy may be and for how long, and the wall comes out unphysically thin, which drives $r_0/\Delta$ to a figure with tens of zeroes in it and the total to their $10^{16}$ galactic masses for a tube a metre across.
That is why the argument about the Krasnikov tube moved, in Krasnikov's own 2003 reply, off the size of the number and onto whether the quantum inequality that bounds it is derived under assumptions a shortcut builder has any reason to respect.

---

## Step 13. The Weyl tensor, computed rather than copied

The Weyl tensor is Riemann with its traces removed,

$$C_{\mu\nu\sigma\tau} = R_{\mu\nu\sigma\tau} - \frac{1}{2}\left(g_{\mu\sigma}R_{\tau\nu} - g_{\mu\tau}R_{\sigma\nu} - g_{\nu\sigma}R_{\tau\mu} + g_{\nu\tau}R_{\sigma\mu}\right) + \frac{R}{6}\left(g_{\mu\sigma}g_{\tau\nu} - g_{\mu\tau}g_{\sigma\nu}\right)$$

in four dimensions, with the same Ricci tensor as Step 8, which is the same contraction as everywhere else in the collection.
This spacetime is not a vacuum anywhere the wall is, so $C$ is not $R$ and copying one into the other would be wrong.
Step 7 already named four blocks of Riemann that vanish identically, and the quickest check that the subtraction is doing real work is to look at one of them.

$R_{t\phi t\phi} = 0$ for every $k$.
The correction terms at those indices are not zero, because $g_{tt} = -1$, $g_{\phi\phi} = r^2$ and $R_{\phi\phi} = -r\,\partial_r w$ are all nonzero, and what comes out is

$$C_{t\phi t\phi} = \frac{2r^2\left(1+k\right)^2\partial_r^2 k - r\left(1+k\right)^2\partial_r k - 2r^2\left(1+k\right)\left(\partial_r k\right)^2 + 2r^2\left(1+k\right)\left(\partial_t^2 k + \partial_t\partial_x k\right) - 2r^2\,\partial_t k\left(\partial_t k + \partial_x k\right)}{6\left(1+k\right)^3},$$

which is not zero.
Twenty four of the slots the entry publishes for $C_{\mu\nu\sigma\tau}$ are like that, nonzero where the corresponding Riemann component vanishes, and no slot goes the other way.
The entry publishes eighty four components with an upper index and seventy two with all four down, and the Weyl tensor is the one tensor in the file that could not have been obtained by reading another block off the page.

The physical reading is the usual one.
Riemann is what the matter does plus what the geometry does on its own, Ricci is what the matter does, and Weyl is the rest.
A tube whose wall carries stress leaves a tidal field in the vacuum on either side of it, and these are its components.

---

## Step 14. The Kretschmann scalar

Contracting Riemann with itself,

$$K = R_{\mu\nu\sigma\tau}R^{\mu\nu\sigma\tau},$$

gives a polynomial in $k$ and its first and second derivatives over $4r^2\left(1+k\right)^6$.
The entry prints it with the combination $\partial_t k + \partial_x k$ collected wherever it appears, which is the derivative along the outbound null ray of Step 11 and turns thirty three separate monomials into five terms grouped by their power of $1+k$.
That is a change of writing and not of content, and its shape says two things that matter.

The first is that it is finite wherever the shape function is smooth and $k$ stays above $-1$.
The tube has no curvature singularity anywhere.
What stands in the way of building one is the source it demands, not a place where the geometry comes apart, which is exactly the situation the warp drives are in as well.

The second is the power of $1+k$ underneath.
$K$ carries a curvature squared, so $\left(1+k\right)^{-6}$ against the $\left(1+k\right)^{-3}$ of the Ricci scalar is the expected bookkeeping, and both of them blow up in the same limit $\delta\to 0$ that makes the round trip fast.
The numerator carries $1+k$ to the fourth power at most, so at fixed derivatives of $k$ the scalar grows like $\left(1+k\right)^{-2}$ as the metric is driven toward degeneracy.
A tube made fast by taking $\delta$ small is a tube whose wall is strongly curved, and Step 12 charges for the same thing in energy.

---

## Step 15. The geodesic equations

With the symbols of Step 6 and the chart dots of Step 1, the equation

$$\ddot{x}^\mu + \Gamma^\mu{}_{\nu\sigma}\dot{x}^\nu\dot{x}^\sigma = 0$$

gives four equations, all printed in the entry.
The two short ones are the transverse pair:

$$\ddot{r} + \tfrac{1}{2}\partial_r k\,\dot{t}\,\dot{x} - \tfrac{1}{2}\partial_r k\,\dot{x}^2 - r\,\dot{\phi}^2 = 0,$$

$$\ddot{\phi} + \frac{2}{r}\dot{r}\,\dot{\phi} = 0.$$

The $\phi$ equation is the flat one, so angular momentum about the axis is conserved exactly, whatever the tube does.
The $r$ equation is more interesting than it looks.
Its two $k$ terms combine into $\tfrac{1}{2}\partial_r k\,\dot{x}\left(\dot{t} - \dot{x}\right)$, which vanishes on the outbound null ray of Step 11, where $\dot{t} = \dot{x}$.
A ray launched along the axis in the outbound direction feels no radial force from the wall at all, even passing through it, which is the same fact as $G_{\mu\nu}\ell_{\text{out}}^\mu\ell_{\text{out}}^\nu = 0$ seen from the other side.
A ray coming home has $\dot{t} = -k\dot{x}$, and the combination becomes $-\tfrac{1}{2}\left(1+k\right)\partial_r k\,\dot{x}^2$, in which the factor $1 + k$ that everything else in the entry divides by appears in the numerator instead.

The $t$ and $x$ equations are longer and carry every symbol of Step 6 with a $1+k$ under it.
The entry prints all four.

---

## Step 16. Every published expression is dimensionally consistent

The dimensions are declared as $\left[t\right] = T$ and $\left[x\right] = \left[r\right] = L$, with $\phi$ and $k$ dimensionless.
Since $t$ is a time, the chart coordinate is $ct$ and the chart dimensions are $L$, $L$, $L$ and $1$.
Every printed derivative is along a chart coordinate, so $\left[\partial_t k\right] = \left[\partial_x k\right] = \left[\partial_r k\right] = L^{-1}$ and every second derivative carries $L^{-2}$.

The line element balances term by term: $c^2dt^2$ is $L^2$, the cross term $\left(1-k\right)c\,dt\,dx$ is $\left(L/T\right)\cdot T\cdot L = L^2$, $k\,dx^2$ and $dr^2$ are $L^2$, and $r^2d\phi^2$ is $L^2$ because $\phi$ is a pure number.

A metric component with two lower indices carries $L^2/\left(\left[x^\mu\right]\left[x^\nu\right]\right)$ in the chart, which is $1$ for every pair here except $\phi\phi$, where it is $L^2$; $g_{\phi\phi} = r^2$ obliges.
A Christoffel symbol carries $L^{-1}$ on top of the index weights, so $\Gamma^r{}_{\phi\phi}$ has to carry $L^{-1}\cdot 1\cdot L^2 = L$, and it is $-r$.
The Einstein tensor carries $L^{-2}$ on top of the index weights, so $G_{tt}$ has to be $L^{-2}$, and each of its three terms is: $\left(1+k\right)\partial_r^2 k$ over $\left(1+k\right)^2$ is $L^{-2}$, $\left(1+k\right)\partial_r k$ over $r\left(1+k\right)^2$ is $L^{-1}/L$, and $r\left(\partial_r k\right)^2$ over $r\left(1+k\right)^2$ is $L^{-2}$.

The geodesic equations are where the dot convention earns its keep, and the $t$ equation shows it in one line.
Its five terms pair $\dot{t}^2$, $\dot{t}\dot{x}$, $\dot{x}^2$, $\dot{t}\dot{r}$ and $\dot{x}\dot{r}$ with coefficients that are all of dimension $L^{-1}$, being single derivatives of $k$ over a dimensionless denominator.
For those five terms to be addable at all, $\dot{t}$, $\dot{x}$ and $\dot{r}$ must carry the same dimension, and the left hand side $\ddot{t}$ then fixes it: with $\dot{t} = d(ct)/d\lambda$ every velocity is $L/\lambda$, every term is $L^{-1}\left(L/\lambda\right)^2 = L/\lambda^2$, and $\ddot{t} = d^2(ct)/d\lambda^2$ is $L/\lambda^2$ as well.
On the other reading, with $\dot{t} = dt/d\lambda$, the first term would carry $T^2/\left(L\lambda^2\right)$ and the third $L/\lambda^2$, and the equation would be adding a time to a length.

Running the dimensional pass alone,

    /tmp/mfs-venv/bin/python _tools/derivations/verify_metrics.py --system krasnikov/cylindrical --dimensions-only

reports nothing, over all three hundred and eighty four published expressions.

---

## Step 17. What the tube is, against the warp drives beside it

Alcubierre's bubble and Natario's are published in this collection alongside this entry, and all three are superluminal geometries bought with negative energy.
The difference is worth stating in components, because it is not a difference of degree.

A warp drive is a shift.
Its line element is $-c^2dt^2 + \left(dx - \beta\,c\,dt\right)^2 + dy^2 + dz^2$, the lapse is one, the spatial slices are flat, and the whole of the geometry sits in $\beta$, a velocity field that carries the flat patch the ship is sitting in.
The ship does not move through the coordinates; the coordinates move it.
That is the source of the causal objection Everett and Roman raise, which is not about energy at all: the pilot at the centre of a bubble is outside the past light cone of its leading wall, so he cannot raise it, steer it or turn it off once he is inside.

The tube is not a shift.
There is no $\beta$ in it, its $x$ is an ordinary coordinate, and nothing in the geometry transports anybody.
The ship crosses $x$ under its own power at whatever speed a rocket can manage, and what it leaves behind it is a region where $g_{xx} = k$ has been driven negative, which is to say a region where the light cones have been opened.
Every component in this entry is a component of a geometry along a path that has already been travelled.
The traveller builds the tube, which answers the control objection, and the price of the answer is that the outbound leg is no faster than light and the whole of the saving is on the way home.

The two also fail the energy conditions in different places.
The warp drive's density is minus the squared transverse gradient of the shape function, so it is negative wherever the wall has any structure and the bubble is symmetric front to back about it.
The tube's failure, by Step 11, is on one null direction and not the other: $\ell_{\text{out}}$ costs nothing and $\ell_{\text{back}}$ costs everything.
That asymmetry is the whole content of the construction, written as a contraction.

---

## Step 18. What the entry publishes, and the declarations the checker needs

The one system publishes three hundred and eighty four expressions: the line element, six metric components and six inverse ones, twenty two Christoffel symbols with an upper index and seventeen with all three down, sixty two Riemann components with an upper index and forty eight with all four down, ten Ricci components in each of three variants, the Ricci scalar, the Kretschmann scalar, ten Einstein components in each of three variants, eighty four Weyl components with an upper index and seventy two with all four down, and four geodesic equations.

Every one of them is a rational function of $k$ and its partial derivatives whose denominator is a power of $1+k$, and each is printed in whichever of the two readings, in $k$ or in $1+k$, is shorter.
The choice is cosmetic and the two are equal as expressions, but it is worth knowing which one is being read: $\Gamma^t{}_{tt}$ prints as $\left(k-1\right)\partial_t k$ over $\left(1+k\right)^2$ because that is one term, while $G_{tt}$ prints its numerator in powers of $1+k$ because grouping it that way turns five terms into three.

For the checker to read the entry, `DIMENSIONS` in `verify_metrics.py` needs one line:

    ("krasnikov", "cylindrical"): {
        "t": "T", "x": "L", "r": "L", "\\phi": "1", "k": "1",
    },

The declaration `"t": "T"` is what tells the checker that $t$ is a time and so that the chart coordinate is $ct$; every factor of $c$ in the comparison follows from it, as does the reading of $\partial_t k$ as a chart derivative.
The declaration `"r": "L"` is the one a reader might guess wrong: $r$ here is the distance from the axis of a cylindrical chart, not an areal radius, and the entry has no spherical coordinate at all.
The solution names no length of its own, so there is nothing else in the line.
The radius of the tube and the thickness of its wall are both hidden inside $k$, which is why $r_0$ and $\Delta$ appear in Step 12 as features of a particular profile and nowhere in the published components.

The entry needs no line in `PARAMETER_RELATIONS`, and it is worth saying why, since this is the one place the entry could have been published with an identity that is not one.
A tube whose switching profile depended on $t$ and $x$ only through $ct - x$ would obey $\partial_t k = -\partial_x k$, and that one relation would kill every term in $\partial_t k + \partial_x k$ in Steps 8 and 14 outright and shorten a good half of the file.
It is not a relation the construction has.
Krasnikov's own $k$ carries a separate factor that opens the tube between the start of the journey and its end, which depends on $x$ by itself, and in any case the relation would be false for the arbitrary $k$ the entry allows.
So no published value uses it, and the checker compares against a free $k$ of three coordinates, which is the claim the entry actually makes.

Running

    /tmp/mfs-venv/bin/python _tools/derivations/verify_metrics.py --system krasnikov/cylindrical

reports no disagreements, no dimensional failures and nothing unchecked, in a little under two minutes.
