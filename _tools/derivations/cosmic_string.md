# The Vilenkin-Gott cosmic string

The Vilenkin-Gott cosmic string has two coordinate systems, the conical exterior of a straight string and Gott's interior cap.
Every component of both, from the line element down to the geodesic equations, follows in order from the one before it.
None of them is asserted without its computation.

The script `_tools/derivations/verify_metrics.py` repeats every computation in sympy and compares each result with the expression derived by hand, so the algebra is checkable by hand and by machine independently.

Outside the string every curvature tensor vanishes, and that vanishing is the result, not a computation left undone.
The exterior of a straight string is flat: the Riemann tensor, the Ricci tensor, the Ricci scalar, the Einstein tensor, the Kretschmann scalar and the Weyl tensor are every one of them zero by the computation of Steps 8 and 9, at every point where the chart is defined.
The gravitation is still there, in the global structure of Steps 10 and 11.
It is the global statement that a circle around the string closes after $2\pi$ of angle has bought only $2\pi(1 - 4G\mu/c^2)$ of circumference, so a wedge of angle $\delta = 8\pi G\mu/c^2$ has been cut out of every plane the string crosses and the cut edges sewn together.
A cone is flat everywhere except at its tip, and this is that statement in four dimensions.

Gott's 1985 interior, the second chart and the subject of Steps 13 to 21, is the tip opened out into a source of finite size.
There the curvature is not zero, and the Einstein tensor of Step 17 is the stress of the string itself.
Matching the interior to the exterior in Step 20 forces the deficit angle to be the total Gaussian curvature of the string's cross section, which is the Gauss-Bonnet theorem, and is why a flat exterior and a curved interior are two readings of the same number.

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

the contraction in use since 2026-09-18.
For the conical exterior the settlement costs nothing, because a tensor that vanishes vanishes on either contraction.
For the interior of Step 17 it decides a sign, and the sign it picks there is the one that gives a sphere positive curvature.

The zeroth coordinate of each chart is

$$x^0 = ct.$$

The index is written with the bare letter $t$, but a component carrying it is a component of that chart.
The metric of Step 5 and everything built from it are computed directly in it.
The dots in the geodesic equations are velocities of that same chart, so $\dot{t}$ means $d(ct)/d\lambda$, and that reading is the one on which every term balances in Step 22.

The string is described by one number, its mass per unit length $\mu$, which for a string is also its tension.
The only dimensionless combination $\mu$ forms with $G$ and $c$ is

$$\frac{G\mu}{c^2},$$

since $[\mu] = M/L$ and $[G] = L^3M^{-1}T^{-2}$ give $[G\mu] = L^2T^{-2} = [c^2]$.
Every property of the string depends on its mass per unit length only through that number.
Two abbreviations run through the exterior,

$$\alpha = 1 - \frac{4G\mu}{c^2}, \qquad \delta = 2\pi(1 - \alpha) = \frac{8\pi G\mu}{c^2},$$

the first the factor the metric carries and the second the deficit angle.
They abbreviate the derivation only.
The checker compares each component written out in $G$, $\mu$ and $c$, because it reads symbol by symbol and has no way to be told what an abbreviation means.
The deficit angle $\delta$ is nevertheless declared a parameter of the spacetime, with its value given in the declaration, so that a reader meets it where the spacetime's other constants are listed rather than buried in a component.

The interior of Steps 13 to 21 carries two numbers instead, the curvature radius $\ell$ of the cap and the polar angle $\chi_0$ it reaches, and by the matching of Step 20 they carry between them no more information than $\mu$ and a choice of core size.

---

## Step 2. The source, and what makes a string a string

A straight static string along the $z$ axis has a stress tensor supported on that axis and invariant under boosts along it.
That last requirement is what distinguishes a string from a line of dust.
A boost along $z$ must leave the source unchanged, because the string is a one dimensional object with no internal structure along its length and nothing marks a preferred rest frame on it.
The only diagonal tensor with that symmetry has equal entries in the $t$ and $z$ slots,

$$T^\mu{}_\nu = \operatorname{diag}\left(-\rho c^2,\, 0,\, 0,\, -\rho c^2\right),$$

which is an energy density $\rho c^2$ together with a tension of the same size along the string, and no pressure across it.
The energy density and the tension being equal is the statement $\mu = $ (tension): the string's tension equals its mass per unit length.

Integrating the energy density across the string's cross section gives the mass per unit length,

$$\mu = \int \rho \, dA,$$

which is the one number the exterior can depend on.
For the idealized string the integrand is a delta function on the axis and the cross section has no size at all.
Gott's interior gives it a size in Step 13, and that is the whole of its content.

The trace of that stress tensor is $T = -2\rho c^2$, so the trace reversed source is

$$T^\mu{}_\nu - \tfrac{1}{2}T\delta^\mu{}_\nu = \operatorname{diag}\left(-\rho c^2,\, 0,\, 0,\, -\rho c^2\right) + \operatorname{diag}\left(\rho c^2, \rho c^2, \rho c^2, \rho c^2\right) = \operatorname{diag}\left(0,\, \rho c^2,\, \rho c^2,\, 0\right).$$

The $t$ slot of the trace reversed source is empty, and that one line is the reason everything that follows is flat.

---

## Step 3. Why the exterior is flat, before any exact solution

The linearized field equation in harmonic gauge is

$$\Box \bar{h}_{\mu\nu} = -\frac{16\pi G}{c^4}T_{\mu\nu},$$

and for a static source the Newtonian potential is read off $h_{tt}$, whose source is the trace reversed stress tensor computed at the end of Step 2.
That component is zero.
A straight string exerts no Newtonian attraction on a particle at rest beside it, because its tension contributes to the active gravitational mass exactly as much as its energy density does and with the opposite sign.
This is Vilenkin's 1981 result, and it is why a string passing through the solar system would not be noticed by anything that watches orbits [vilenkin1981].

What is left is the spatial part, sourced by $\rho c^2$ in the $x$ and $y$ slots.
Solving the two dimensional Poisson equation gives a logarithm, and the resulting metric is

$$ds^2 = -c^2dt^2 + dz^2 + \left(1 - \frac{8G\mu}{c^2}\ln\frac{r}{r_0}\right)\left(dr^2 + r^2d\phi^2\right)$$

to first order in $G\mu/c^2$.
A logarithm looks like curvature, but the transverse part of that line element is conformally flat in two dimensions, which every two dimensional metric is, and the conformal factor here is exactly the one a flat cone acquires when it is written in isothermal coordinates.
Substituting

$$\bar{r} = \left(1 - \frac{4G\mu}{c^2}\right)^{-1} r\left(\frac{r}{r_0}\right)^{-4G\mu/c^2}$$

and dropping terms of second order turns it into

$$ds^2 = -c^2dt^2 + dz^2 + d\bar{r}^2 + \left(1 - \frac{4G\mu}{c^2}\right)^2\bar{r}^2 d\phi^2,$$

which is Step 4's line element.
The form arrived at by linearizing is exact: it solves the full field equations, because the field equations outside the source are $R_{\mu\nu} = 0$ and this metric has $R^\mu{}_{\nu\rho\sigma} = 0$ by Step 8, which is a great deal more than the vacuum equations ask for.
There is nothing left for a higher order correction to correct.

---

## Step 4. The conical line element, and which of its two forms to use

The exterior is

$$ds^2 = -c^2dt^2 + dr^2 + \left(1 - \frac{4G\mu}{c^2}\right)^2 r^2 d\phi^2 + dz^2, \qquad \phi \in [0, 2\pi),$$

with $t$ a time, $r$ the proper distance from the string, $\phi$ an angle around it and $z$ the proper distance along it.
Every component of the exterior is computed in this form.

The other form in circulation writes the same geometry as

$$ds^2 = -c^2dt^2 + dr^2 + r^2 d\tilde\phi^2 + dz^2, \qquad \tilde\phi \in \left[0, 2\pi\alpha\right),$$

with $\tilde\phi = \alpha\phi$.
The two are the same spacetime, carried one into the other by the substitution of Step 10.
The first is the better choice, for three reasons.

Its angle is periodic in $2\pi$, so the chart has the same periodicity as every other cylindrical or spherical chart, and a reader does not have to carry a nonstandard identification around with the components.

Its deviation from Minkowski space is visible in a metric component.
In the second form every metric component, every Christoffel symbol and every curvature component is identical to flat space's, and the entire physics has been moved into the range of a coordinate, which is a fact about the chart's domain rather than about any component.
A spacetime whose whole content is a footnote to the domain of a coordinate looks empty, and the exterior already has eight curvature tensors, counting each index position separately, that vanish identically.

The tension appears where the field equation put it.
$g_{\phi\phi} = \alpha^2r^2$ is the only component in which $\mu$ occurs at all, in either chart, which says as plainly as it can be said that a string bends nothing but the closing up of angles around it.

---

## Step 5. The metric matrix, its determinant and its inverse

In the chart $(x^0, x^1, x^2, x^3) = (ct, r, \phi, z)$ the line element of Step 4 is diagonal,

$$g_{\mu\nu} = \operatorname{diag}\left(-1,\, 1,\, \alpha^2r^2,\, 1\right),$$

which written out is

$$g_{tt} = -1, \qquad g_{rr} = 1, \qquad g_{\phi\phi} = \left(1 - \frac{4G\mu}{c^2}\right)^2r^2, \qquad g_{zz} = 1.$$

The factor of $c^2$ that multiplies $dt^2$ in the line element is the one the chart absorbs: $-c^2dt^2 = -(d(ct))^2$, so the chart component is $-1$ and carries no dimensions, as every chart component of a metric must.

The determinant is

$$\det g = -\alpha^2r^2,$$

negative for every $r > 0$ and every $\alpha > 0$, so the chart is a good Lorentzian one away from the axis.
It degenerates at $r = 0$, where the string is, and at $\alpha = 0$, which is $G\mu/c^2 = 1/4$ and is where the deficit angle would eat the whole of the plane.
Observationally $G\mu/c^2$ is smaller than about $10^{-7}$, so $\alpha$ is within a part in ten million of one and the cone is a very shallow one.

The inverse of a diagonal matrix is the diagonal matrix of reciprocals,

$$g^{\mu\nu} = \operatorname{diag}\left(-1,\, 1,\, \frac{1}{\alpha^2r^2},\, 1\right),$$

written out as

$$g^{tt} = -1, \qquad g^{rr} = 1, \qquad g^{\phi\phi} = \frac{1}{\left(1 - \dfrac{4G\mu}{c^2}\right)^2r^2}, \qquad g^{zz} = 1.$$

---

## Step 6. The Christoffel symbols of the exterior, with an upper index

Only one metric component depends on any coordinate:

$$\partial_r g_{\phi\phi} = 2\alpha^2 r,$$

and every other first derivative of every component is zero.
A Christoffel symbol is a sum of such derivatives, so any symbol without at least one $r$ index and two $\phi$ indices among its three slots vanishes.
That leaves two.

$$\Gamma^r{}_{\phi\phi} = \tfrac{1}{2}g^{rr}\left(\partial_\phi g_{r\phi} + \partial_\phi g_{r\phi} - \partial_r g_{\phi\phi}\right) = -\tfrac{1}{2}\left(2\alpha^2r\right) = -\alpha^2 r.$$

$$\Gamma^\phi{}_{r\phi} = \tfrac{1}{2}g^{\phi\phi}\left(\partial_r g_{\phi\phi} + \partial_\phi g_{\phi r} - \partial_\phi g_{r\phi}\right) = \frac{1}{2\alpha^2r^2}\left(2\alpha^2 r\right) = \frac{1}{r}.$$

The symbol is symmetric in its lower pair, so $\Gamma^\phi{}_{\phi r} = \Gamma^\phi{}_{r\phi}$, and the checker compares it twice, once for each ordering.
Written out, the nonzero symbols are

$$\Gamma^r{}_{\phi\phi} = -\left(1 - \frac{4G\mu}{c^2}\right)^2r, \qquad \Gamma^\phi{}_{r\phi} = \Gamma^\phi{}_{\phi r} = \frac{1}{r}.$$

These are the symbols of flat space in cylindrical coordinates with the single factor $\alpha^2$ in the first of them.
Nothing carries a factor of $c$, because no symbol has a time index: the time part of the metric is a constant, so the chart's factors of $c$ never have anything to multiply.

---

## Step 7. The lowered Christoffel symbols

Lowering the upper index on a diagonal metric multiplies by one component:

$$\Gamma_{r\phi\phi} = g_{rr}\Gamma^r{}_{\phi\phi} = -\alpha^2r,$$

$$\Gamma_{\phi r\phi} = \Gamma_{\phi\phi r} = g_{\phi\phi}\Gamma^\phi{}_{r\phi} = \alpha^2r^2\cdot\frac{1}{r} = \alpha^2 r,$$

which written out are

$$\Gamma_{r\phi\phi} = -\left(1 - \frac{4G\mu}{c^2}\right)^2r, \qquad \Gamma_{\phi r\phi} = \Gamma_{\phi\phi r} = \left(1 - \frac{4G\mu}{c^2}\right)^2r.$$

The lowered symbol is symmetric in its last two indices only, and the sign flip between $\Gamma_{r\phi\phi}$ and $\Gamma_{\phi r\phi}$ is the usual antisymmetry of the first two indices of $\Gamma_{\mu\nu\rho} + \Gamma_{\nu\mu\rho} = \partial_\rho g_{\mu\nu} = 0$ here, since $g_{r\phi}$ does not depend on anything.

---

## Step 8. The Riemann tensor of the exterior, computed and zero

Two Christoffel symbols is few enough that the Riemann tensor can be done by hand in full.
The only pair of directions in which anything depends on anything is $(r, \phi)$, so the only components that can be nonzero are the ones built from those two.

The first:

$$R^r{}_{\phi r\phi} = \partial_r \Gamma^r{}_{\phi\phi} - \partial_\phi \Gamma^r{}_{\phi r} + \Gamma^r{}_{r\lambda}\Gamma^\lambda{}_{\phi\phi} - \Gamma^r{}_{\phi\lambda}\Gamma^\lambda{}_{\phi r}.$$

The second term is zero because $\Gamma^r{}_{\phi r}$ is zero, and the third is zero because $\Gamma^r{}_{r\lambda}$ is zero for every $\lambda$.
The fourth runs over $\lambda$ and survives only for $\lambda = \phi$.
So

$$R^r{}_{\phi r\phi} = \partial_r\left(-\alpha^2 r\right) - \Gamma^r{}_{\phi\phi}\Gamma^\phi{}_{\phi r} = -\alpha^2 - \left(-\alpha^2 r\right)\left(\frac{1}{r}\right) = -\alpha^2 + \alpha^2 = 0.$$

The derivative term and the quadratic term cancel exactly, for every value of $\alpha$.

The second:

$$R^\phi{}_{r\phi r} = \partial_\phi \Gamma^\phi{}_{rr} - \partial_r \Gamma^\phi{}_{r\phi} + \Gamma^\phi{}_{\phi\lambda}\Gamma^\lambda{}_{rr} - \Gamma^\phi{}_{r\lambda}\Gamma^\lambda{}_{r\phi}.$$

The first and third terms vanish because $\Gamma^\lambda{}_{rr} = 0$ for every $\lambda$, and the last survives only for $\lambda = \phi$.
So

$$R^\phi{}_{r\phi r} = -\partial_r\left(\frac{1}{r}\right) - \left(\frac{1}{r}\right)\left(\frac{1}{r}\right) = \frac{1}{r^2} - \frac{1}{r^2} = 0.$$

Every other component is zero for one of three reasons.
A component with a $t$ or a $z$ index needs a Christoffel symbol carrying that index, and there is none.
A component with three or four distinct spatial indices cannot be formed from the two indices $r$ and $\phi$.
A component with the index pattern of one of those two, in any of the orderings the symmetries of Riemann relate them by, is that one up to a sign and a factor of the metric, and is therefore zero as well.

So

$$R^\mu{}_{\nu\rho\sigma} = 0, \qquad R_{\mu\nu\rho\sigma} = 0$$

identically, everywhere on $r > 0$, in both index positions.
This is the exact statement that the linearized calculation of Step 3 could only suggest.

---

## Step 9. Everything built from Riemann, including Weyl

The Ricci tensor is a contraction of Riemann, so

$$R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu} = 0$$

in all three index positions, $R_{\mu\nu}$, $R^\mu{}_\nu$ and $R^{\mu\nu}$.
The Ricci scalar is a contraction of that,

$$R = g^{\mu\nu}R_{\mu\nu} = 0,$$

and the Einstein tensor is built from both,

$$G_{\mu\nu} = R_{\mu\nu} - \tfrac{1}{2}Rg_{\mu\nu} = 0,$$

again in all three index positions.
The Kretschmann scalar is quadratic in Riemann,

$$K = R_{\mu\nu\rho\sigma}R^{\mu\nu\rho\sigma} = 0.$$

That implication runs one way only.
A flat metric has $K = 0$, which is this line, but a spacetime with $K = 0$ need not be flat, because in Lorentzian signature $K$ is not a sum of squares: a pp wave has every curvature invariant zero and a Riemann tensor that is not.
The exterior's $K = 0$ is a consequence of Step 8 and never an argument for it.

The Weyl tensor is the one curvature tensor that is not a contraction of Riemann but a subtraction from it,

$$C_{\mu\nu\rho\sigma} = R_{\mu\nu\rho\sigma} - \tfrac{1}{2}\left(g_{\mu\rho}R_{\nu\sigma} - g_{\mu\sigma}R_{\nu\rho} - g_{\nu\rho}R_{\mu\sigma} + g_{\nu\sigma}R_{\mu\rho}\right) + \tfrac{1}{6}R\left(g_{\mu\rho}g_{\nu\sigma} - g_{\mu\sigma}g_{\nu\rho}\right),$$

and it is never filled in by copying Riemann into it, a warning set out in full in `_tools/derivations/weyl.md`.
Here every term on the right is separately zero, by Step 8 for the first and by the vanishing Ricci tensor and Ricci scalar for the other two, so

$$C_{\mu\nu\rho\sigma} = 0, \qquad C^\mu{}_{\nu\rho\sigma} = 0.$$

The exterior is conformally flat, which is not news for a spacetime that is flat outright, but it follows from the computation rather than from an assumption.

Eight curvature tensors, counting each index position separately, and two scalars vanish, and that vanishing is the physics of the exterior rather than an omission.

---

## Step 10. The cone is flat, and is not Minkowski space

Define

$$\tilde\phi = \alpha\phi = \left(1 - \frac{4G\mu}{c^2}\right)\phi.$$

Then $\alpha^2r^2d\phi^2 = r^2d\tilde\phi^2$ and the line element becomes

$$ds^2 = -c^2dt^2 + dr^2 + r^2d\tilde\phi^2 + dz^2,$$

which is Minkowski space in cylindrical coordinates, component for component.
That is the second proof of Step 8 and it takes one line.

What the substitution does not do is make the spacetime Minkowski space.
The coordinate $\phi$ is an angle and the points $\phi$ and $\phi + 2\pi$ are the same point, so $\tilde\phi$ and $\tilde\phi + 2\pi\alpha$ are the same point, and the new angle closes early:

$$\tilde\phi \in [0, 2\pi\alpha), \qquad 2\pi\alpha = 2\pi - \frac{8\pi G\mu}{c^2} = 2\pi - \delta.$$

A wedge of opening angle $\delta$ is missing, and the two faces of the cut are the same points.
Take a flat sheet of paper, cut out a wedge of angle $\delta$ with its apex at the origin, tape the cut edges together, and the result is a cone: intrinsically flat everywhere the paper is, with all of the geometry concentrated at a point that carries no paper at all.
Multiply that cone by a flat $t$ and a flat $z$ and the spacetime is built.

The invariant statement is about holonomy.
Parallel transport a vector once around any loop encircling the string, at any radius, and it comes back rotated through the angle $\delta$, in the plane transverse to the string and always by the same amount.
Around a loop that does not encircle the string it comes back unrotated, because the region the loop bounds is flat and the rotation is the integral of the curvature over it.
Curvature is a local quantity, holonomy is not, and the string is entirely holonomy.
Nothing measured inside any one simply connected neighbourhood of the exterior can tell it from empty space.
A tidal measurement returns zero, a gyroscope returns nothing, a gravimeter returns nothing.
Only a trip the whole way around, or a comparison between two travellers who went opposite ways, turns the missing wedge up.

The string itself is the apex.
At $r = 0$ the chart fails and the curvature, read as a distribution rather than as a function, is a delta function on the axis with total strength $\delta$, which is Step 20's statement in the exterior's language.
It is a conical singularity and not a curvature singularity: every scalar built from the curvature is zero on approach, $K$ included, and the geometry is extendable as a cone in a way a black hole's centre is not.

---

## Step 11. What a deficit angle does, since it is all the string does

Three consequences, each a global statement and each measurable in principle.

A circle around the string at proper radius $r$, at fixed $t$ and $z$, has circumference

$$\oint\sqrt{g_{\phi\phi}}\,d\phi = \int_0^{2\pi}\alpha r\,d\phi = 2\pi\alpha r = 2\pi r - \delta r,$$

against a proper radius of exactly $r$, since $g_{rr} = 1$ makes $r$ the distance to the axis.
The ratio of circumference to radius is short of $2\pi$ by $\delta$ at every radius, which is the sense in which the deficit is not a feature of any particular place.

Two null geodesics leaving a source behind the string and passing it on opposite sides both travel in straight lines, because Step 10 makes the exterior flat, and yet they meet the observer from directions that differ.
The angle between them is $\delta$ in the limit of a source infinitely far behind a string perpendicular to the line of sight, and in general it is reduced by the ratio of distances,

$$\Delta\theta = \frac{8\pi G\mu}{c^2}\cdot\frac{d_{\mathrm{LS}}}{d_{\mathrm{S}}}\sin\vartheta,$$

with $d_{\mathrm{LS}}$ the distance from string to source, $d_{\mathrm{S}}$ from observer to source and $\vartheta$ the angle between the string and the line of sight [gott1985].
The lensing is achromatic, the two images are of equal brightness and neither is magnified, because each ray has travelled through flat space the whole way and nothing has focused anything.
That signature is unlike every other lens in astronomy, and it is why a string would be unmistakable if one were ever found.
For $G\mu/c^2 = 10^{-6}$, the value a grand unified scale string would have, $\delta$ is $2.5\times 10^{-5}$ radians, or about five arcseconds.

The third leads to closed timelike curves.
A circuit around the string is shorter than a circuit through flat space would be, by the missing wedge, so a traveller going around can arrive before a light signal sent across.
Around one static string that costs nothing, since the saving is in distance and not in time and no closed causal curve results.
Gott's 1991 construction takes two such strings moving past one another at high speed and arranges for a circuit around both to be a closed timelike curve [gott1991].
Nothing about a single straight string misbehaves.
The misbehaviour is what two such strings can be made to do to each other.

---

## Step 12. The geodesic equations of the exterior

The geodesic equation is

$$\ddot{x}^\mu + \Gamma^\mu{}_{\nu\rho}\dot{x}^\nu\dot{x}^\rho = 0,$$

with $\Gamma$ the symbols of Step 6 and the dot a derivative with respect to an affine parameter $\lambda$ of the chart coordinates, so that $\dot{t}$ is $d(ct)/d\lambda$.
Two of the four equations have no Christoffel symbol to collect:

$$\ddot{t} = 0, \qquad \ddot{z} = 0.$$

Coordinate time and height along the string are both affine, which is the statement that $\partial_t$ and $\partial_z$ are covariantly constant Killing vectors here.

The radial equation collects $\Gamma^r{}_{\phi\phi}$:

$$\ddot{r} - \left(1 - \frac{4G\mu}{c^2}\right)^2r\dot{\phi}^2 = 0.$$

The angular equation collects $\Gamma^\phi{}_{r\phi}$ twice, once for each ordering of the lower indices, which is where the factor of two comes from:

$$\ddot{\phi} + \frac{2}{r}\dot{r}\dot{\phi} = 0.$$

The second of these is the conservation of

$$L = \alpha^2 r^2 \dot\phi,$$

since $\dot L = \alpha^2(2r\dot r\dot\phi + r^2\ddot\phi) = \alpha^2r^2(\ddot\phi + 2\dot r\dot\phi/r) = 0$, and $L$ is the angular momentum per unit mass conjugate to the Killing vector $\partial_\phi$.
Substituting it into the radial equation gives $\ddot r = L^2/(\alpha^2r^3)$, which integrates to a straight line in the transverse plane, as Step 10 says it must.

Every one of these four is the flat space equation, with $\alpha^2$ in the one place the metric carries it.
A particle near a cosmic string moves exactly as a free particle moves, and only a particle that goes around it, or two that go around it different ways, ever learns otherwise.

---

## Step 13. Gott's interior, which is where the mass is

The exterior of Step 4 has all of its source at $r = 0$, in a line of zero thickness, which is an idealization and not a spacetime.
Gott's 1985 solution replaces the apex of the cone with a smooth cap and joins the two [gott1985].

The interior is a piece of a two sphere, carried rigidly along the string and along time.
Write $\chi$ for the polar angle measured from the axis and $\ell$ for the radius of the sphere, so the transverse section of the string is a spherical cap of radius $\ell$ reaching out to $\chi = \chi_0$.
The line element is

$$ds^2 = -c^2dt^2 + \ell^2d\chi^2 + \ell^2\sin^2\chi\,d\phi^2 + dz^2, \qquad \chi \in [0, \chi_0], \quad \phi \in [0, 2\pi).$$

The proper distance from the axis is $\ell\chi$, so the core has proper radius $\ell\chi_0$.
The chart is regular at $\chi = 0$ in the way polar coordinates are: the circumference of a small circle is $2\pi\ell\sin\chi \to 2\pi\ell\chi$, which is $2\pi$ times the proper radius, so there is no deficit and no apex inside the core.

Writing the cap with a dimensionless polar angle rather than with a proper radius $\varrho = \ell\chi$, in which the same metric reads $-c^2dt^2 + d\varrho^2 + \ell^2\sin^2(\varrho/\ell)d\phi^2 + dz^2$, is a presentational choice and nothing more.
The angle is the natural coordinate of a sphere and it keeps every trigonometric argument a bare symbol, which is what spherical charts do with $\theta$.
It is also the form the checker can read: the LaTeX dialect `verify_metrics.py` parses takes $\sin^2\chi$, a squared trigonometric function applied to a bare name, and does not take $\sin^2(\varrho/\ell)$, whose argument is an expression, so the second form would have to be spelled out as a parenthesized power and would be reported `UNCHECKED` if it were not.
A spacetime whose natural radial variable is an angle should be written in that angle.

The spacetime is a metric product of the flat two dimensional $(t, z)$ plane along the string with the two sphere of radius $\ell$ across it.
Nothing in either factor depends on anything in the other, and that single fact decides the connection and the curvature of Steps 15 to 19.

---

## Step 14. The metric of the cap and its inverse

In the chart $(ct, \chi, \phi, z)$,

$$g_{\mu\nu} = \operatorname{diag}\left(-1,\, \ell^2,\, \ell^2\sin^2\chi,\, 1\right),$$

written out as

$$g_{tt} = -1, \qquad g_{\chi\chi} = \ell^2, \qquad g_{\phi\phi} = \ell^2\sin^2\chi, \qquad g_{zz} = 1.$$

Both $\chi$ and $\phi$ are dimensionless, so both of their metric components carry $\ell^2$ and both are still dimensionless as chart components must be, since $[g_{\mu\nu}] = L^2/([x^\mu][x^\nu])$ and each of those coordinates contributes nothing.

The determinant is $\det g = -\ell^4\sin^2\chi$, which is negative for $0 < \chi < \pi$, so the chart is good on the whole cap except at the axis $\chi = 0$, where it fails exactly as polar coordinates fail at their origin and for the same harmless reason.

The inverse is

$$g^{\mu\nu} = \operatorname{diag}\left(-1,\, \frac{1}{\ell^2},\, \frac{1}{\ell^2\sin^2\chi},\, 1\right).$$

---

## Step 15. The Christoffel symbols of the cap

Only $g_{\phi\phi}$ depends on a coordinate:

$$\partial_\chi g_{\phi\phi} = 2\ell^2\sin\chi\cos\chi = \ell^2\sin 2\chi.$$

As in Step 6 this leaves two symbols with an upper index.

$$\Gamma^\chi{}_{\phi\phi} = -\tfrac{1}{2}g^{\chi\chi}\partial_\chi g_{\phi\phi} = -\frac{1}{2\ell^2}\left(2\ell^2\sin\chi\cos\chi\right) = -\sin\chi\cos\chi.$$

$$\Gamma^\phi{}_{\chi\phi} = \Gamma^\phi{}_{\phi\chi} = \tfrac{1}{2}g^{\phi\phi}\partial_\chi g_{\phi\phi} = \frac{2\ell^2\sin\chi\cos\chi}{2\ell^2\sin^2\chi} = \cot\chi.$$

These are the symbols of the unit two sphere, with $\chi$ in the role $\theta$ usually plays and no dependence on $\ell$ at all, which is what it means for the radius to be a scale rather than a shape.
Lowering the index:

$$\Gamma_{\chi\phi\phi} = g_{\chi\chi}\Gamma^\chi{}_{\phi\phi} = -\ell^2\sin\chi\cos\chi,$$

$$\Gamma_{\phi\chi\phi} = \Gamma_{\phi\phi\chi} = g_{\phi\phi}\Gamma^\phi{}_{\chi\phi} = \ell^2\sin^2\chi\cot\chi = \ell^2\sin\chi\cos\chi.$$

With the index up or down there are three nonzero symbols, and none carries a $t$ or a $z$ index, because the flat factor of the product contributes no connection.

---

## Step 16. The Riemann tensor of the cap

The product structure of Step 13 says the answer before the computation: the Riemann tensor of a metric product is the sum of the Riemann tensors of the factors, one of which is flat, so everything sits in the sphere.
The computation is two components, as in Step 8, and this time they do not cancel.

$$R^\chi{}_{\phi\chi\phi} = \partial_\chi\Gamma^\chi{}_{\phi\phi} - \Gamma^\chi{}_{\phi\phi}\Gamma^\phi{}_{\phi\chi} = \partial_\chi\left(-\sin\chi\cos\chi\right) - \left(-\sin\chi\cos\chi\right)\cot\chi$$

$$= \left(\sin^2\chi - \cos^2\chi\right) + \cos^2\chi = \sin^2\chi.$$

$$R^\phi{}_{\chi\phi\chi} = -\partial_\chi\Gamma^\phi{}_{\chi\phi} - \Gamma^\phi{}_{\chi\phi}\Gamma^\phi{}_{\chi\phi} = \csc^2\chi - \cot^2\chi = 1.$$

The two surviving orderings follow from the antisymmetry in the last pair, $R^\mu{}_{\nu\rho\sigma} = -R^\mu{}_{\nu\sigma\rho}$, so there are four nonzero components with an upper index:

$$R^\chi{}_{\phi\chi\phi} = \sin^2\chi, \quad R^\chi{}_{\phi\phi\chi} = -\sin^2\chi, \quad R^\phi{}_{\chi\phi\chi} = 1, \quad R^\phi{}_{\chi\chi\phi} = -1.$$

Lowering the first index gives the fully covariant components,

$$R_{\chi\phi\chi\phi} = g_{\chi\chi}R^\chi{}_{\phi\chi\phi} = \ell^2\sin^2\chi, \qquad R_{\phi\chi\phi\chi} = g_{\phi\phi}R^\phi{}_{\chi\phi\chi} = \ell^2\sin^2\chi,$$

with the two sign flipped orderings beside them, and the equality of those two is the pair symmetry $R_{\mu\nu\rho\sigma} = R_{\rho\sigma\mu\nu}$ checking itself.

The Gaussian curvature of the transverse section can be read straight off:

$$K_{\mathrm{G}} = \frac{R_{\chi\phi\chi\phi}}{g_{\chi\chi}g_{\phi\phi} - g_{\chi\phi}^2} = \frac{\ell^2\sin^2\chi}{\ell^4\sin^2\chi} = \frac{1}{\ell^2},$$

which is the curvature of a sphere of radius $\ell$, as it had better be.
Its integral over the cap is the deficit angle of Step 20.

---

## Step 17. The Ricci tensor, the scalar and the Einstein tensor

Contract on the first index, the convention of Step 1:

$$R_{\chi\chi} = R^\alpha{}_{\chi\alpha\chi} = R^\phi{}_{\chi\phi\chi} = 1,$$

$$R_{\phi\phi} = R^\alpha{}_{\phi\alpha\phi} = R^\chi{}_{\phi\chi\phi} = \sin^2\chi,$$

with $R_{tt} = R_{zz} = 0$ because no Riemann component carries a $t$ or a $z$ index.
Raising gives the other two index positions,

$$R^\chi{}_\chi = R^\phi{}_\phi = \frac{1}{\ell^2}, \qquad R^{\chi\chi} = \frac{1}{\ell^4}, \qquad R^{\phi\phi} = \frac{1}{\ell^4\sin^2\chi},$$

and the mixed form is the clearest of the three: the Ricci tensor is $1/\ell^2$ on the sphere directions and zero on the flat ones.

The Ricci scalar is

$$R = g^{\mu\nu}R_{\mu\nu} = \frac{1}{\ell^2} + \frac{\sin^2\chi}{\ell^2\sin^2\chi} = \frac{2}{\ell^2},$$

which is twice the Gaussian curvature of Step 16, as the scalar curvature of a two dimensional factor always is.
The sign is where the contraction convention earns its keep.
On the other contraction every value in this step reverses, and a sphere would be reported as having negative scalar curvature, which is wrong in the plainest possible sense.

The Einstein tensor is $G_{\mu\nu} = R_{\mu\nu} - \tfrac{1}{2}Rg_{\mu\nu}$, and the two sphere directions cancel exactly:

$$G_{\chi\chi} = 1 - \frac{1}{\ell^2}\cdot\ell^2 = 0, \qquad G_{\phi\phi} = \sin^2\chi - \frac{1}{\ell^2}\cdot\ell^2\sin^2\chi = 0,$$

while the two flat directions, which had no Ricci at all, pick up the trace term:

$$G_{tt} = 0 - \frac{1}{\ell^2}(-1) = \frac{1}{\ell^2}, \qquad G_{zz} = 0 - \frac{1}{\ell^2}(1) = -\frac{1}{\ell^2}.$$

In mixed form,

$$G^t{}_t = G^z{}_z = -\frac{1}{\ell^2},$$

and the other two slots are empty.
Put that beside the field equation $G^\mu{}_\nu = (8\pi G/c^4)T^\mu{}_\nu$ and the stress tensor of Step 2:

$$-\frac{1}{\ell^2} = \frac{8\pi G}{c^4}\left(-\rho c^2\right) \implies \frac{1}{\ell^2} = \frac{8\pi G\rho}{c^2}.$$

The interior is a uniform energy density $\rho c^2$ with an equal tension along $z$ and no pressure across the string, which is precisely the source a string must have by Step 2, and the equality of the $t$ and $z$ slots of $G^\mu{}_\nu$ is not something that was arranged: it fell out of the geometry.
A cap of a sphere multiplied by a flat plane is a cosmic string interior, and it is one for no other reason than that its Einstein tensor has two equal negative entries in the flat directions and nothing anywhere else.
The radius is fixed by the density,

$$\ell = \frac{c}{\sqrt{8\pi G\rho}},$$

and that is the only place $\rho$ enters, since every component of the interior is written in $\ell$.

---

## Step 18. The Kretschmann scalar of the cap

The Riemann tensor has one independent component, which by the symmetries appears in four index orderings, each with the same square.
So

$$K = R_{\mu\nu\rho\sigma}R^{\mu\nu\rho\sigma} = 4\left(R_{\chi\phi\chi\phi}\right)^2\left(g^{\chi\chi}g^{\phi\phi}\right)^2 = 4\left(\ell^2\sin^2\chi\right)^2\left(\frac{1}{\ell^2}\cdot\frac{1}{\ell^2\sin^2\chi}\right)^2 = \frac{4}{\ell^4}.$$

It is constant over the cap, as it must be for a homogeneous source, and equals $R^2$, which is the relation a product with a single curved two dimensional factor always satisfies.
Nowhere on the cap is there a singularity, and the axis $\chi = 0$ that the chart loses is an ordinary interior point of the string.

Against the exterior's $K = 0$, this says that all of the local curvature this spacetime has is inside the string, and outside it there is none at any order.

---

## Step 19. The Weyl tensor of the cap, computed rather than copied

The interior is the one part of the spacetime with a nonzero Weyl tensor, and it is emphatically not a copy of Riemann.
Riemann has four nonzero components, all of them in the $(\chi, \phi)$ directions.
Weyl has twenty four, and they occupy every pair of directions the spacetime has, including the two in which Riemann is empty.

With $R = 2/\ell^2$ and the Ricci tensor of Step 17, the definition

$$C_{\mu\nu\rho\sigma} = R_{\mu\nu\rho\sigma} - \tfrac{1}{2}\left(g_{\mu\rho}R_{\nu\sigma} - g_{\mu\sigma}R_{\nu\rho} - g_{\nu\rho}R_{\mu\sigma} + g_{\nu\sigma}R_{\mu\rho}\right) + \tfrac{1}{6}R\left(g_{\mu\rho}g_{\nu\sigma} - g_{\mu\sigma}g_{\nu\rho}\right)$$

gives six independent components, one for each pair of coordinate directions.
Four of them start from a Riemann component that is zero.

$$C_{t\chi t\chi} = 0 - \tfrac{1}{2}\left(g_{tt}R_{\chi\chi} + g_{\chi\chi}R_{tt}\right) + \tfrac{1}{6}R\,g_{tt}g_{\chi\chi} = \tfrac{1}{2} - \frac{1}{3\ell^2}\ell^2 = \frac{1}{6}.$$

$$C_{t\phi t\phi} = -\tfrac{1}{2}\left(g_{tt}R_{\phi\phi}\right) + \tfrac{1}{6}R\,g_{tt}g_{\phi\phi} = \frac{\sin^2\chi}{2} - \frac{\sin^2\chi}{3} = \frac{\sin^2\chi}{6}.$$

$$C_{tztz} = 0 - 0 + \tfrac{1}{6}\cdot\frac{2}{\ell^2}\cdot(-1)(1) = -\frac{1}{3\ell^2}.$$

$$C_{\chi z\chi z} = -\tfrac{1}{2}\left(g_{zz}R_{\chi\chi}\right) + \tfrac{1}{6}R\,g_{\chi\chi}g_{zz} = -\tfrac{1}{2} + \frac{1}{3} = -\frac{1}{6},$$

and likewise $C_{\phi z\phi z} = -\sin^2\chi/6$.
The one that starts from a nonzero Riemann component is

$$C_{\chi\phi\chi\phi} = \ell^2\sin^2\chi - \tfrac{1}{2}\left(g_{\chi\chi}R_{\phi\phi} + g_{\phi\phi}R_{\chi\chi}\right) + \tfrac{1}{6}R\,g_{\chi\chi}g_{\phi\phi} = \ell^2\sin^2\chi - \ell^2\sin^2\chi + \frac{\ell^2\sin^2\chi}{3} = \frac{\ell^2\sin^2\chi}{3},$$

a third of what Riemann has there, not all of it.
The twenty four nonzero components are these six, each in the four orderings the antisymmetry in the first pair and in the last pair relate, with the sign that antisymmetry demands, and the mixed form $C^\mu{}_{\nu\rho\sigma}$ has the same components with the first index raised.

Two checks.
The Weyl tensor is traceless on every pair, $g^{\mu\rho}C_{\mu\nu\rho\sigma} = 0$, and the $\nu\sigma = tt$ trace reads

$$g^{\chi\chi}C_{\chi t\chi t} + g^{\phi\phi}C_{\phi t\phi t} + g^{zz}C_{ztzt} = \frac{1}{\ell^2}\cdot\frac{1}{6} + \frac{1}{\ell^2\sin^2\chi}\cdot\frac{\sin^2\chi}{6} - \frac{1}{3\ell^2} = 0,$$

which it is.
And in the null tetrad

$$l = \tfrac{1}{\sqrt 2}\left(\partial_{ct} + \partial_z\right), \quad n = \tfrac{1}{\sqrt 2}\left(\partial_{ct} - \partial_z\right), \quad m = \tfrac{1}{\sqrt 2}\left(\frac{1}{\ell}\partial_\chi + \frac{i}{\ell\sin\chi}\partial_\phi\right),$$

normalized by $l\cdot n = -1$ and $m\cdot\bar m = 1$, the five Weyl scalars come out

$$\Psi_0 = \Psi_1 = \Psi_3 = \Psi_4 = 0, \qquad \Psi_2 = -\frac{1}{6\ell^2},$$

so the interior is Petrov type D, with its two repeated principal null directions along the string.
That is the same algebraic type as the Schwarzschild and Kerr spacetimes, reached here by a source of finite density rather than by a vacuum.

The Bertotti-Robinson spacetime makes the contrast, since it too is a metric product of two two dimensional factors.
There the two factors are an anti-de Sitter plane and a sphere of the same radius, whose curvatures are equal and opposite, the two contributions to Weyl cancel and the spacetime is conformally flat.
Here the flat factor contributes nothing to cancel against the sphere's $1/\ell^2$, and what is left over is that $\Psi_2$.
A Weyl tensor is a difference, and it takes the whole of both factors to know it.

---

## Step 20. Matching the cap to the cone, and the Gauss-Bonnet reading

The cap of Step 13 and the cone of Step 4 are joined at the edge of the core.
Continuity of the induced metric and of its normal derivative across that surface, which are the Israel conditions for a boundary carrying no surface layer, reduce here to two statements about the circumference function.

Let $\varrho$ be the proper distance from the axis in both charts, $\varrho = \ell\chi$ inside and $\varrho = r + \text{const}$ outside.
Write $C(\varrho)/2\pi = \sqrt{g_{\phi\phi}}$ for the circumference radius.
Inside it is $\ell\sin\chi$, whose derivative with respect to $\varrho = \ell\chi$ is $\cos\chi$.
Outside it is $\alpha r$, whose derivative with respect to $\varrho$ is $\alpha$.
Matching the function fixes the exterior's radial offset and matching its derivative gives the condition that matters:

$$\alpha = \cos\chi_0, \qquad \text{that is} \qquad \cos\chi_0 = 1 - \frac{4G\mu}{c^2}.$$

The cone's opening is the slope of the cap at its rim, which is geometrically obvious once seen: the cone is the tangent developable of the sphere along the circle $\chi = \chi_0$, the shape a sheet of paper takes when wrapped tangent to a ball.

Now compute $\mu$ independently, by integrating the density of Step 17 over the cap.
The area element is $\ell^2\sin\chi\,d\chi\,d\phi$, so

$$\mu = \int\rho\,dA = \rho\int_0^{2\pi}\!\!\int_0^{\chi_0}\ell^2\sin\chi\,d\chi\,d\phi = 2\pi\rho\ell^2\left(1 - \cos\chi_0\right).$$

Substituting $\rho = c^2/(8\pi G\ell^2)$ from Step 17,

$$\mu = \frac{c^2}{8\pi G\ell^2}\cdot 2\pi\ell^2\left(1 - \cos\chi_0\right) = \frac{c^2}{4G}\left(1 - \cos\chi_0\right),$$

which rearranges to $\cos\chi_0 = 1 - 4G\mu/c^2$.
That is the matching condition again, arrived at from the source rather than from the junction, and the agreement is the solution's consistency.

The deficit angle is therefore

$$\delta = 2\pi\left(1 - \cos\chi_0\right) = \frac{8\pi G\mu}{c^2},$$

which is the solid angle the cap subtends at the centre of its sphere, and also, by Step 16, the total Gaussian curvature of the transverse section:

$$\int_{\mathrm{cap}}K_{\mathrm{G}}\,dA = \frac{1}{\ell^2}\cdot 2\pi\ell^2\left(1 - \cos\chi_0\right) = 2\pi\left(1 - \cos\chi_0\right) = \delta.$$

This is the Gauss-Bonnet theorem for the transverse two dimensional surface, and it is the sharpest statement about the string.
The deficit angle of a cone is the integral of the curvature over whatever was used to cap it, and the answer does not depend on how the capping was done.
Make the core smaller and denser at fixed $\mu$ and $\ell$ shrinks with $\chi_0$ growing to keep the product fixed; replace the cap with some other profile of the same total mass and the exterior does not notice.
In the limit of zero thickness the curvature becomes a delta function on the axis of total strength $\delta$, which is the distributional statement at the end of Step 10.
The exterior is flat because all of the curvature has been pushed into a region that shrank to nothing, and the one thing that survived the shrinking is its integral.

---

## Step 21. The geodesic equations of the cap

With the symbols of Step 15 and the same convention as Step 12,

$$\ddot{t} = 0, \qquad \ddot{z} = 0,$$

$$\ddot{\chi} - \sin\chi\cos\chi\,\dot{\phi}^2 = 0,$$

$$\ddot{\phi} + 2\cot\chi\,\dot{\chi}\dot{\phi} = 0.$$

These are the great circle equations of the two sphere, with $t$ and $z$ carried along freely, which is what the product structure requires.
A particle at rest inside the string stays at rest, and one moving across the core follows a great circle of the cap until it reaches the rim and hands over to the straight lines of Step 12.

The conserved quantity of the angular equation is $\ell^2\sin^2\chi\,\dot\phi$, which matches the exterior's $\alpha^2r^2\dot\phi$ at the rim once $\alpha = \cos\chi_0$ is used, since both are $\sqrt{g_{\phi\phi}}^2\dot\phi$ and $g_{\phi\phi}$ is continuous there by Step 20.

---

## Step 22. Every expression is dimensionally consistent

The chart coordinates of the exterior are $(ct, r, \phi, z)$, of dimensions $L$, $L$, $1$, $L$, and its parameters are $[\mu] = M/L$, $[G] = L^3M^{-1}T^{-2}$ and $[\delta] = 1$.
Those of the interior are $(ct, \chi, \phi, z)$, of dimensions $L$, $1$, $1$, $L$, with $[\ell] = L$, $[\chi_0] = 1$, $[\rho] = M/L^3$ and the same $G$.
A metric component then carries

$$[g_{\mu\nu}] = \frac{L^2}{[x^\mu][x^\nu]},$$

a Christoffel symbol carries $L^{-1}$ corrected by $[x^\mu]/L$ per upper index and $L/[x^\mu]$ per lower one, a Riemann, Ricci, Einstein or Weyl component carries $L^{-2}$ corrected the same way, and $K$ carries $L^{-4}$.

Five samples.

$g_{\phi\phi} = (1 - 4G\mu/c^2)^2r^2$ must carry $L^2/(1\cdot 1) = L^2$.
The bracket is dimensionless because $[G\mu/c^2] = L^2T^{-2}/(L^2T^{-2}) = 1$, as in Step 1, and $r^2$ carries $L^2$.
This is the check that decides that $\mu$ is a mass per unit length: no other dimension lets $G\mu/c^2$ stand beside the $1$.

$\Gamma^\phi{}_{r\phi} = 1/r$ carries $L^{-1}$ from the field, $1/L$ from the upper $\phi$, $1$ from the lower $r$ and $L$ from the lower $\phi$, so it must carry $L^{-1}$, and it does.

$R_{\chi\chi} = 1$ carries $L^{-2}$ from the field and $L/1$ twice from its two dimensionless lower indices, so it must be dimensionless, and it is.
This is why a dimensionless radial coordinate is not a problem for the checker: the index weights absorb it exactly.

$C_{t\chi t\chi} = 1/6$ carries $L^{-2}$ from the field, $L/L = 1$ from each $t$ index and $L/1 = L$ from each $\chi$, so it must be dimensionless, while $C_{tztz} = -1/(3\ell^2)$ has two length indices in place of the two angles and must carry $L^{-2}$, which it does.
Those two components of one tensor carrying different dimensions is the index weighting working, not an inconsistency.

$R = 2/\ell^2$ carries $L^{-2}$ and $K = 4/\ell^4$ carries $L^{-4}$, as a scalar of each rank must.

The geodesic equations are measured against their own second derivatives.
In the radial equation of Step 12, $\ddot{r}$ carries $L/\lambda^2$, and the term $\alpha^2r\dot\phi^2$ carries

$$L\cdot\left(\frac{1}{\lambda}\right)^2 = \frac{L}{\lambda^2},$$

which matches, because $\dot\phi$ is $d\phi/d\lambda$ with $\phi$ dimensionless.
In the angular equation of Step 21, $\ddot\phi$ carries $\lambda^{-2}$ and $\cot\chi\,\dot\chi\dot\phi$ carries $1\cdot\lambda^{-1}\cdot\lambda^{-1}$, which matches too.
The equations $\ddot t = 0$ are where the chart convention shows: $\ddot t$ carries $L/\lambda^2$ only because $\dot t$ means $d(ct)/d\lambda$, and on the other reading it would carry $T/\lambda^2$ and the four equations of a single geodesic would not be four equations of the same kind.
Neither system has a term mixing a time velocity with a space velocity, since no Christoffel symbol of either carries a time index, so the cosmic string is a weak test of that convention and a clean one.

---

## Step 23. The components, and the declarations the checker needs

The conical exterior has twenty one expressions for the checker: the line element, four metric components and four inverse ones, three Christoffel symbols in each of two index positions, two scalars and four geodesic equations.
Its Riemann, Ricci, Einstein and Weyl tensors vanish in two index positions each, eight zero tensors in all, and the checker tests them as hard as any nonzero component, because it requires every other component to vanish.

The interior has eighty nine: the line element, four metric components and four inverse ones, three Christoffel symbols in each of two index positions, four Riemann components in each of two, two Ricci components in each of three, two Einstein components in each of three, two scalars, twenty four Weyl components in each of two index positions and four geodesic equations.

For the checker to read either of them, `DIMENSIONS` in `verify_metrics.py` needs one line per system:

    ("cosmic_string", "conical"): {
        "t": "T", "r": "L", "\\phi": "1", "z": "L",
        "\\mu": "M/L", "G": "L**3/(M*T**2)", "\\delta": "1",
    },
    ("cosmic_string", "interior_cap"): {
        "t": "T", "\\chi": "1", "\\phi": "1", "z": "L",
        "\\ell": "L", "\\chi_0": "1", "\\rho": "M/L**3", "\\mu": "M/L",
        "G": "L**3/(M*T**2)",
    },

The declaration `"t": "T"` is what tells the checker that $t$ is a time and so that the chart coordinate is $ct$.
Declaring $\mu$ as a mass per unit length and $G$ with its own dimensions, rather than folding them into a length as most other spacetimes do with $r_s = 2GM/c^2$, is what makes the deficit come out as the dimensionless group the physics is usually quoted in; the Vaidya spacetime is the only other one that keeps $G$ and a mass apart in this way.
The declarations of $\delta$, $\chi_0$, $\rho$ and $\mu$ in the interior name parameters no component of the interior uses, because every one of its components is written in $\ell$.
They are declared because they are the spacetime's parameters all the same, and the checker requires its declaration table to agree exactly with the parameters of the spacetime, which is what keeps the two from drifting.

The domain of each system includes a fact rather than an interval, the exterior's flatness and the interior's matching condition, and those are written inside `\text{}`.
Every domain condition is typeset as mathematics, so a word left outside `\text{}` comes out as a row of italic symbols with the spaces stripped out; the horizon condition of de Sitter space is the pattern.

Neither system needs a line in `PARAMETER_RELATIONS`.
The exterior's $\mu$ is free, every value of it below $c^2/4G$ giving a cone, and the interior's $\ell$ and $\chi_0$ are free in the same way, with Step 20's $\cos\chi_0 = 1 - 4G\mu/c^2$ relating them to the exterior's parameter rather than constraining either of them within its own system.
Every component of each system holds for all values of its own parameters, which is the condition for the table not to be needed.

Running

    /tmp/mfs-venv/bin/python _tools/derivations/verify_metrics.py --system cosmic_string/conical --system cosmic_string/interior_cap

reports no disagreements and no dimensional failures.
