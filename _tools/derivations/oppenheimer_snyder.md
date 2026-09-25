# Oppenheimer-Snyder collapse

This is the working behind the two coordinate systems in `MFS/assets/data/metrics/oppenheimer_snyder.json`.
Every value the entry prints is derived here, in order, from the line elements down to the geodesic equations.
Nothing is left as an exercise and nothing is asserted that is not computed.

The companion script `_tools/derivations/verify_metrics.py` does the same work in sympy and compares it against the published file, so the algebra below is checkable by hand and by machine independently.

This entry is the only one in the collection that is two spacetimes rather than one.
Neither half is new.
The interior is the closed case of the Friedmann-Lemaître-Robertson-Walker dust cosmology, which `MFS/assets/data/metrics/frw.json` already publishes, written here in the chart that makes the surface of the star a coordinate surface.
The exterior is the Schwarzschild vacuum, which `MFS/assets/data/metrics/schwarzschild.json` already publishes, in the same spherical chart and with the same $r_s$.
The work of this document is therefore not in either curvature block but in the seam between them, and Steps 11 to 15 are the seam.

One equation comes out of it,

$$r_s = a_m\sin^3\chi_0,$$

binding the exterior's mass to the interior's scale factor, and it is the whole of the matching.
Step 15 reads it three ways: as $M = \tfrac{4}{3}\pi\rho R^3$ holding at every moment of the collapse, as the compactness $\sin^2\chi_0 = r_s/R_0$ of the star at release, and as the redshift $\cos\chi_0 = \sqrt{1 - r_s/R_0}$ a distant observer sees from its surface before it starts to fall.
Steps 16 to 18 are what the joined solution then says: the star falls to a singularity in a finite proper time, its surface crosses the Schwarzschild radius at a conformal time $\eta = \pi - 2\chi_0$ with time to spare before the singularity, and the same crossing takes forever in the exterior chart, the star fading with an e folding time of $4GM/c^3$.

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
The settlement costs nothing in the exterior, where every one of these tensors vanishes on either contraction, and it decides a sign in the interior, where the dust is.
On this contraction the interior's Einstein tensor comes out as

$$G_{\tau\tau} = 3\frac{\dot{a}^2+1}{a^2} = \frac{8\pi G\rho}{c^2},$$

which is the Friedmann equation with a positive density on the right.
On the other contraction it would be the same equation with a negative one, and Step 7 is where that is checked rather than assumed.

The chart is the collection's, the one whose zeroth coordinate is

$$x^0 = c\tau \quad\text{inside}, \qquad x^0 = ct \quad\text{outside}.$$

The index is printed with the bare letter, but the component printed against it is a component of that chart, and everything from Step 4 onward is computed directly in it.
A dot on the scale factor is a derivative with respect to the chart time,

$$\dot{a} = \frac{da}{d(c\tau)} = \frac{1}{c}\frac{da}{d\tau},$$

so that $\dot{a}$ is dimensionless and $\ddot{a}$ carries one inverse length.
The dots in the geodesic equations are velocities of the same chart, so $\dot{\tau}$ means $d(c\tau)/d\lambda$ and $\dot{t}$ means $d(ct)/d\lambda$.
Step 21 checks that this is the reading on which every published term balances.

Two conventions of this document only, never of the published file.
A subscript $\Sigma$ marks a quantity evaluated on the surface of the star.
A square bracket marks the jump of a quantity across that surface, $[X] = X_{\text{out}} - X_{\text{in}}$, both sides taken in the same coordinates on the surface.

---

## Step 2. The two pieces, and why there have to be two

A star of dust has a boundary, and a boundary is where a solution of the Einstein equations stops being the solution it is.
Inside, the source is a ball of pressureless matter of uniform density, and the field equations there are the field equations of a dust cosmology.
Outside, the source is nothing at all, and Birkhoff's theorem leaves no freedom: the only spherically symmetric vacuum is Schwarzschild, static whether or not the star inside it is [birkhoff1923].
That last point is the one worth pausing on, because it is what makes the problem solvable at all.
The star may be collapsing at any rate it likes, and the geometry outside it does not care, does not radiate, and does not change.

So the spacetime is two exact solutions with a common boundary, and the only question is what the field equations demand where they meet.
They demand two things, and Steps 11 to 14 impose them:

* the induced metric on the surface is the same computed from either side, so the surface has one geometry rather than two;
* the extrinsic curvature of the surface is the same computed from either side, so there is no thin shell of matter sitting at the join.

The first is Darmois' condition and the second is Darmois' too, though it is most often met in Israel's form, which reads the failure of the second as the stress tensor of a surface layer,

$$S_{ij} = -\frac{c^4}{8\pi G}\left(\left[K_{ij}\right] - h_{ij}\left[K\right]\right).$$

Asking for no surface layer is asking for $[K_{ij}] = 0$.
The Oppenheimer-Snyder star has all of its matter in its interior and none of it painted on its boundary, so that is the condition to impose [oppenheimersnyder1939].

The static version of the same construction is already in the collection.
`MFS/assets/data/metrics/interior_schwarzschild.json` is a ball of uniform density held up by pressure and matched to the same exterior, and the pressure it needs is what diverges when the ball is too compact.
This entry is what happens when there is no pressure to begin with.

---

## Step 3. The interior line element

The interior is homogeneous and isotropic about its centre, so its spatial slices are the three geometries that are, and the collapse of a ball of dust released from rest picks the closed one.
The reason is the field equation rather than a choice, and Step 8 shows it: a solution with $\dot{a} = 0$ at some moment needs $k > 0$ to satisfy the Friedmann equation with a positive density, because the equation at that moment reads $k = 8\pi G\rho a^2/(3c^2) > 0$.
Writing the closed slice as a three sphere of radius $a$ and a dimensionless polar angle $\chi$,

$$ds^2 = -c^2d\tau^2 + a(\tau)^2\left[d\chi^2 + \sin^2\chi\left(d\theta^2 + \sin^2\theta\,d\phi^2\right)\right].$$

This is the same geometry `frw.json` publishes at $k = +1$, under $r = \sin\chi$, which turns $dr^2/(1-r^2)$ into $d\chi^2$.
The two charts differ in where the length sits.
There $r$ is a length and $a$ is dimensionless; here $\chi$ is an angle and $a$ carries the length, which is the normalisation a collapse wants, since $a$ is then a radius that shrinks to zero and $a_m$ is a radius the star is released from.

The star is the ball $\chi \le \chi_0$, and $\chi_0$ is a constant because the dust is comoving: a particle of dust sits at fixed $\chi$ for all time, and so does the outermost one.
The areal radius of the sphere at $\chi$, meaning the radius read off the area $4\pi\tilde{r}^2$ of that sphere, is

$$\tilde{r}(\tau,\chi) = a(\tau)\sin\chi,$$

and the surface of the star has areal radius

$$R(\tau) = a(\tau)\sin\chi_0.$$

The whole of the matching is written in that one function.

Two bounds on $\chi_0$ are worth stating before they are used.
The chart needs $\chi_0 < \pi$ for the ball to be a ball.
The physics needs $\chi_0 < \pi/2$, because $\tilde{r} = a\sin\chi$ has its maximum at $\chi = \pi/2$, and a star whose surface sat beyond the equator of the three sphere would be one whose spheres get smaller as you move outward through the matter, which is a star already inside a trapped region at the moment it is released.
Step 15 gets the same bound from the exterior, where it reads $R_0 > r_s$.

---

## Step 4. The metric and its inverse

In the chart $x^0 = c\tau$ the line element of Step 3 is

$$ds^2 = -\left(dx^0\right)^2 + a^2d\chi^2 + a^2\sin^2\chi\,d\theta^2 + a^2\sin^2\chi\sin^2\theta\,d\phi^2,$$

so the metric is diagonal with

$$g_{\tau\tau} = -1, \qquad g_{\chi\chi} = a^2, \qquad g_{\theta\theta} = a^2\sin^2\chi, \qquad g_{\phi\phi} = a^2\sin^2\chi\sin^2\theta,$$

and the inverse is the reciprocal of each,

$$g^{\tau\tau} = -1, \qquad g^{\chi\chi} = \frac{1}{a^2}, \qquad g^{\theta\theta} = \frac{1}{a^2\sin^2\chi}, \qquad g^{\phi\phi} = \frac{1}{a^2\sin^2\chi\sin^2\theta}.$$

The factor of $c^2$ that stood in front of $d\tau^2$ is in the chart coordinate now, which is what makes $g_{\tau\tau}$ the dimensionless $-1$ rather than $-c^2$.
The determinant is $-a^6\sin^4\chi\sin^2\theta$, which is nonzero on $0 < \chi < \pi$, $0 < \theta < \pi$ and $a > 0$, so the chart is good everywhere except at the two axes, where a spherical chart always fails, and at $a = 0$, which is not a coordinate failure at all but Step 9's singularity.

---

## Step 5. The Christoffel symbols

Only $a$ depends on the chart time and only the angular blocks depend on $\chi$ and $\theta$, so the symbols come in three families.
From the definition in Step 1, with $\partial_0 = \partial/\partial(c\tau)$,

$$\Gamma^\tau{}_{\chi\chi} = -\tfrac{1}{2}g^{\tau\tau}\partial_0 g_{\chi\chi} = a\dot{a},$$

and the same computation on the two angular slots gives

$$\Gamma^\tau{}_{\theta\theta} = a\dot{a}\sin^2\chi, \qquad \Gamma^\tau{}_{\phi\phi} = a\dot{a}\sin^2\chi\sin^2\theta.$$

The mixed symbols are the logarithmic rate of the scale factor,

$$\Gamma^\chi{}_{\tau\chi} = \Gamma^\theta{}_{\tau\theta} = \Gamma^\phi{}_{\tau\phi} = \frac{\dot{a}}{a},$$

each symmetric in its lower pair, and the purely spatial ones are those of the unit three sphere, untouched by $a$ because it cancels between the metric and its inverse,

$$\Gamma^\chi{}_{\theta\theta} = -\sin\chi\cos\chi, \qquad \Gamma^\chi{}_{\phi\phi} = -\sin\chi\cos\chi\sin^2\theta, \qquad \Gamma^\theta{}_{\chi\theta} = \Gamma^\phi{}_{\chi\phi} = \cot\chi,$$

$$\Gamma^\theta{}_{\phi\phi} = -\sin\theta\cos\theta, \qquad \Gamma^\phi{}_{\theta\phi} = \cot\theta.$$

Counting the two orderings of each mixed lower pair, that is the eighteen symbols the entry publishes.

The lowered variant follows by one multiplication by the diagonal metric, which flips the sign of the three with a time index and multiplies the spatial ones by $a^2$,

$$\Gamma_{\tau\chi\chi} = -a\dot{a}, \qquad \Gamma_{\chi\tau\chi} = a\dot{a}, \qquad \Gamma_{\chi\theta\theta} = -a^2\sin\chi\cos\chi, \qquad \Gamma_{\theta\chi\theta} = a^2\sin\chi\cos\chi,$$

and so on through the list, eighteen again.

$\Gamma^\chi{}_{\tau\tau} = 0$ is worth naming, because Step 12 leans on it: a curve of constant $\chi$, $\theta$ and $\phi$ is a geodesic, so every particle of the dust, including every particle on the surface, is in free fall.
That is what pressureless means, and it is built into the chart.

---

## Step 6. The Riemann tensor

Two families again.
The time-space components come from differentiating $\Gamma^\chi{}_{\tau\chi}$,

$$R^\tau{}_{\chi\tau\chi} = \partial_0\Gamma^\tau{}_{\chi\chi} - \Gamma^\tau{}_{\chi\lambda}\Gamma^\lambda{}_{\chi\tau} = \left(\dot{a}^2 + a\ddot{a}\right) - \dot{a}^2 = a\ddot{a},$$

and the purely spatial ones combine the curvature of the three sphere with the square of the expansion rate,

$$R^\chi{}_{\theta\chi\theta} = \left(\dot{a}^2+1\right)\sin^2\chi.$$

The $1$ in that bracket is the curvature of the unit three sphere and would be $k$ in the chart `frw.json` uses.
Everything else follows from these two by the symmetries and by the factors of $\sin^2\chi$ and $\sin^2\theta$ the angular metric carries:

$$R^\tau{}_{\theta\tau\theta} = a\ddot{a}\sin^2\chi, \qquad R^\tau{}_{\phi\tau\phi} = a\ddot{a}\sin^2\chi\sin^2\theta, \qquad R^\chi{}_{\tau\chi\tau} = R^\theta{}_{\tau\theta\tau} = R^\phi{}_{\tau\phi\tau} = -\frac{\ddot{a}}{a},$$

$$R^\theta{}_{\chi\theta\chi} = R^\phi{}_{\chi\phi\chi} = \dot{a}^2+1, \qquad R^\theta{}_{\phi\theta\phi} = \left(\dot{a}^2+1\right)\sin^2\chi\sin^2\theta, \qquad R^\phi{}_{\theta\phi\theta} = \left(\dot{a}^2+1\right)\sin^2\chi.$$

With the antisymmetry in the last pair that is twenty four nonzero components, and the entry publishes all twenty four rather than half of them, because the checker requires every component it does not find to vanish.
The fully lowered variant is these multiplied by $g_{\mu\mu}$ on the first index, which flips the sign of the six with a leading time index and multiplies the spatial ones by $a^2$ or $a^2\sin^2\chi$:

$$R_{\tau\chi\tau\chi} = -a\ddot{a}, \qquad R_{\chi\theta\chi\theta} = a^2\left(\dot{a}^2+1\right)\sin^2\chi, \qquad R_{\theta\phi\theta\phi} = a^2\left(\dot{a}^2+1\right)\sin^4\chi\sin^2\theta,$$

twenty four again.

---

## Step 7. Ricci, the scalar, Einstein, and the field equations

Contracting on the first lower index as Step 1 says,

$$R_{\tau\tau} = R^\chi{}_{\tau\chi\tau} + R^\theta{}_{\tau\theta\tau} + R^\phi{}_{\tau\phi\tau} = -3\frac{\ddot{a}}{a},$$

$$R_{\chi\chi} = R^\tau{}_{\chi\tau\chi} + R^\theta{}_{\chi\theta\chi} + R^\phi{}_{\chi\phi\chi} = a\ddot{a} + 2\dot{a}^2 + 2,$$

with $R_{\theta\theta}$ and $R_{\phi\phi}$ the same bracket carrying $\sin^2\chi$ and $\sin^2\chi\sin^2\theta$.
The scalar is

$$R = g^{\mu\nu}R_{\mu\nu} = 3\frac{\ddot{a}}{a} + 3\frac{a\ddot{a}+2\dot{a}^2+2}{a^2} = 6\frac{a\ddot{a}+\dot{a}^2+1}{a^2},$$

and the Einstein tensor is

$$G_{\tau\tau} = R_{\tau\tau} + \tfrac{1}{2}R = 3\frac{\dot{a}^2+1}{a^2}, \qquad G_{\chi\chi} = -\left(2a\ddot{a}+\dot{a}^2+1\right),$$

the second derivative cancelling out of the time slot, as it must, because $G_{\tau\tau}$ is a constraint on a slice rather than an equation of motion.

The source is pressureless dust of density $\rho$, at rest in this chart, so its four velocity is $u^\mu = (c,0,0,0)$ and

$$T^{\mu\nu} = \rho\, u^\mu u^\nu, \qquad T_{\tau\tau} = \rho c^2, \qquad T_{\chi\chi} = T_{\theta\theta} = T_{\phi\phi} = 0.$$

With $G_{\mu\nu} = 8\pi G T_{\mu\nu}/c^4$ the two independent equations are

$$3\frac{\dot{a}^2+1}{a^2} = \frac{8\pi G\rho}{c^2}, \qquad 2a\ddot{a}+\dot{a}^2+1 = 0.$$

The first is the Friedmann equation of a closed dust universe and the second is the vanishing of the pressure.
Both signs are the contraction's doing.
Had the Ricci tensor been contracted on the last index instead, every component of this step would carry the opposite sign, the first equation would read $3(\dot{a}^2+1)/a^2 = -8\pi G\rho/c^2$, and a ball of ordinary matter would have to have a negative density to collapse.
That is the concrete content of the settlement named in Step 1.

The entry publishes these components with $a(\tau)$ left free, not with the dust condition imposed.
That is deliberate, and it is the same choice `frw.json` makes.
A published $G_{\chi\chi} = 0$ would be a true statement about the Oppenheimer-Snyder solution and a false statement about the line element it is printed under, and the checker compares against the line element.

---

## Step 8. The first integral, and the cycloid

Multiply the pressureless condition by $\dot{a}$ and read it backwards:

$$2a\dot{a}\ddot{a} + \dot{a}^3 + \dot{a} = \frac{d}{d(c\tau)}\left[a\left(\dot{a}^2+1\right)\right] = 0,$$

so the bracket is a constant of the collapse.
Calling that constant $a_m$,

$$\dot{a}^2 + 1 = \frac{a_m}{a}.$$

The name is chosen by what the equation says: $\dot{a}^2 \ge 0$ forces $a \le a_m$, with equality exactly where $\dot{a} = 0$, so $a_m$ is the largest the star ever is, and a star released from rest is released at $a = a_m$.
Comparing with the Friedmann equation of Step 7 gives the same constant as a statement about the matter,

$$a_m = \frac{8\pi G\rho a^3}{3c^2},$$

which says that $\rho a^3$ does not change, the dust being neither created nor destroyed as the ball shrinks.

The solution is the cycloid, most easily written with a parameter $\eta$,

$$a = \frac{a_m}{2}\left(1+\cos\eta\right), \qquad c\tau = \frac{a_m}{2}\left(\eta+\sin\eta\right),$$

running from $\eta = 0$, where $a = a_m$ and $\dot{a} = 0$, to $\eta = \pi$, where $a = 0$.
It is a solution because

$$\frac{da}{d(c\tau)} = \frac{da/d\eta}{d(c\tau)/d\eta} = \frac{-\tfrac{1}{2}a_m\sin\eta}{\tfrac{1}{2}a_m(1+\cos\eta)} = -\frac{\sin\eta}{1+\cos\eta},$$

and then

$$\dot{a}^2 + 1 = \frac{\sin^2\eta + (1+\cos\eta)^2}{(1+\cos\eta)^2} = \frac{2+2\cos\eta}{(1+\cos\eta)^2} = \frac{2}{1+\cos\eta} = \frac{a_m}{a}.$$

The parameter is conformal time: $c\,d\tau = a\,d\eta$ by the second equation, so the line element of Step 3 is $a(\eta)^2$ times the metric of a static three sphere,

$$ds^2 = a(\eta)^2\left[-d\eta^2 + d\chi^2 + \sin^2\chi\,d\Omega^2\right],$$

which is why $\eta$ is the right parameter for anything to do with light.
A radial null ray has $d\chi = \pm\,d\eta$ in it, exactly, with no reference to $a$ at all.
Step 17 is entirely that observation.

The entry does not publish this chart.
It publishes the comoving one and declares $a_m$ as a parameter, because the components of the conformal chart would be the same geometry written twice, and because the free function $a(\tau)$ is the honest statement of what the published curvature depends on.

---

## Step 9. Kretschmann, the singularity, and a Weyl tensor that is zero

The Kretschmann scalar of the interior, from the fully lowered Riemann tensor of Step 6 and its fully raised partner, is

$$K = R_{\mu\nu\rho\sigma}R^{\mu\nu\rho\sigma} = 12\frac{a^2\ddot{a}^2 + \left(\dot{a}^2+1\right)^2}{a^4},$$

the first term coming from the three independent blocks with a time index and the second from the three purely spatial ones, each block counted four times over by the two antisymmetries, which is where both twelves come from.
On the dust solution of Step 8 it collapses to a single power.
Differentiating the first integral gives $\ddot{a} = -a_m/(2a^2)$, and substituting both it and $\dot{a}^2+1 = a_m/a$,

$$K = \frac{12}{a^4}\left[\frac{a_m^2}{4a^2} + \frac{a_m^2}{a^2}\right] = \frac{15\,a_m^2}{a^6}.$$

That diverges as $a \to 0$, so $\eta = \pi$ is a curvature singularity and not a failure of the chart.
Every particle of the dust reaches it at the same proper time, since $a$ depends on $\tau$ alone, which is the one feature of this model that a real star would not share.

The Weyl tensor of the interior vanishes identically,

$$C^\mu{}_{\nu\rho\sigma} = 0,$$

in every component and at every moment of the collapse, and the entry publishes both variants as empty blocks.
This is computed rather than quoted: subtracting the Ricci traces from the Riemann tensor of Step 6,

$$C_{\mu\nu\rho\sigma} = R_{\mu\nu\rho\sigma} - \left(g_{\mu[\rho}R_{\sigma]\nu} - g_{\nu[\rho}R_{\sigma]\mu}\right) + \tfrac{1}{3}R\,g_{\mu[\rho}g_{\sigma]\nu},$$

every slot cancels.
The reason it must is that a homogeneous isotropic slicing is conformally flat, which Step 8 showed in the most direct possible way: the line element is a function of $\eta$ times a static metric, and the static metric it multiplies, the Einstein universe, is itself conformally flat.
So all of the interior's curvature is Ricci, which is to say all of it is the dust, and none of it is tidal.
A small ball of test particles falling with the dust is squeezed but not distorted.

One practical note for whoever writes the next entry with a three sphere in it.
Two of the Weyl components come out of sympy's `simplify` as

$$\frac{a\dot{a}\left(\sin 2\chi\tan\chi + \cos 2\chi - 1\right)}{2\tan\chi},$$

which is zero, since $\sin2\chi\tan\chi = 2\sin^2\chi$ and $\cos2\chi - 1 = -2\sin^2\chi$, but is not recognisably zero to `simplify` alone.
It is the exact form the docstring of `norm` in `verify_metrics.py` names as its reason for existing, and `norm` does reduce it, which is why the empty published block passes.
A derivation that trusted a bare `simplify` would have published two components that are not there.

---

## Step 10. The exterior

Outside the surface the solution is Schwarzschild,

$$ds^2 = -\left(1-\frac{r_s}{r}\right)c^2dt^2 + \frac{dr^2}{1-\dfrac{r_s}{r}} + r^2d\theta^2 + r^2\sin^2\theta\,d\phi^2, \qquad r_s = \frac{2GM}{c^2},$$

in the standard spherical chart, and the entry publishes it: thirteen Christoffel symbols in each variant, twenty four Riemann components in each, an empty Ricci tensor and an empty Einstein tensor in all three variants each, $R = 0$, $K = 12r_s^2/r^6$, twenty four Weyl components in each variant and four geodesic equations.

None of that is derived here, because `MFS/assets/data/metrics/schwarzschild.json` publishes the same chart and `verify_metrics.py` checks it, and repeating the derivation would be repeating a page the collection already carries.
Two things about it do belong here.

The first is the domain.
The Schwarzschild entry publishes the exterior region $r > r_s$ of an eternal black hole.
This entry publishes $r \ge R(\tau)$ of a collapsing star, which is a different region of the same local geometry: while the star is large the chart never reaches $r_s$ at all, and the part of the eternal solution that this spacetime does not contain, the white hole and the second asymptotic region, is replaced by the star.
The published components are identical because the geometry is; the domain is the whole of the difference, and it is the reason this entry does not simply point at the other one.

The exterior time runs over $t \ge 0$ and not the whole line; `oppenheimer_snyder.json` gave $t \in (-\infty, \infty)$ until 24 September 2026.
The interior starts at the moment of rest, $\tau \ge 0$ with $a(0) = a_m$ and $\dot{a}(0) = 0$, and the model says nothing about what came before: a star held at rest by something the model leaves out, or the time reverse of the collapse, are equally consistent with it, and it contains neither.
So before release there is no $R(t)$ and no surface to be outside of, and a domain claiming the whole of $t$ claimed more than the model supports.
The origin of Schwarzschild time is free, and it is put at the release, so that $t = 0$ and $\tau = 0$ are the same instant on the surface.
That instant is a slice of time symmetry on both sides at once: inside because $\dot{a} = 0$ there, outside because every slice of constant $t$ in the static chart has vanishing extrinsic curvature, and on the surface because it is at rest, $dR/dt = 0$, the geodesic of energy $\tilde{E} = c^2\cos\chi_0$ that Steps 14 and 15 put it on being at its turning point $R_0$.
The upper end needs no bound: the surface takes an infinite $t$ to reach $r_s$, as Step 18 shows, so $R(t) > r_s$ at every $t$ the coordinates reach, and the part of the black hole outside the star, $r < r_s$ after the surface has crossed, lies in neither set of coordinates.

The second is the Weyl tensor.
It is equal to the Riemann tensor here, component for component, and that is a fact about a vacuum rather than a licence to copy a block: with $R_{\mu\nu} = 0$ and $R = 0$ every correction term in the definition quoted in Step 9 is zero, so $C_{\mu\nu\rho\sigma} = R_{\mu\nu\rho\sigma}$ identically.
It was computed from Riemann minus its traces for this entry, in sympy, and it came out equal, all twenty four components of it, which is what a vacuum requires.
So the two halves of this spacetime carry complementary curvature: the interior is all Ricci and no Weyl, the exterior is all Weyl and no Ricci.

---

## Step 11. The surface, and the first junction condition

The surface $\Sigma$ is a timelike hypersurface, three dimensional, and it is described twice.

From inside it is $\chi = \chi_0$, and the coordinates that run along it are $(\tau,\theta,\phi)$.
Reading the interior line element of Step 3 at fixed $\chi$,

$$ds^2_\Sigma\big|_{\text{in}} = -c^2d\tau^2 + a(\tau)^2\sin^2\chi_0\left(d\theta^2+\sin^2\theta\,d\phi^2\right).$$

From outside it is a curve in the $(t,r)$ plane crossed with the spheres, parametrised as $r = R(\lambda)$, $t = T(\lambda)$.
Reading the exterior line element of Step 10 along it,

$$ds^2_\Sigma\big|_{\text{out}} = \left[-\left(1-\frac{r_s}{R}\right)c^2\left(\frac{dT}{d\lambda}\right)^2 + \frac{1}{1-\dfrac{r_s}{R}}\left(\frac{dR}{d\lambda}\right)^2\right]d\lambda^2 + R^2\left(d\theta^2+\sin^2\theta\,d\phi^2\right).$$

Setting the two equal, slot by slot, is the first junction condition.

The angular slots give

$$R = a\sin\chi_0,$$

which is the statement already made in Step 3, that the areal radius of the surface is the areal radius of the surface.
It has content anyway: it says the exterior's $r$ coordinate of the surface is not free, but is the interior's scale factor times a fixed sine, so the entire history $R(\tau)$ is known as soon as $a(\tau)$ is.

The time slot gives

$$\left(1-\frac{r_s}{R}\right)c^2\left(\frac{dT}{d\lambda}\right)^2 - \frac{1}{1-\dfrac{r_s}{R}}\left(\frac{dR}{d\lambda}\right)^2 = c^2\left(\frac{d\tau}{d\lambda}\right)^2,$$

which, taking $\lambda = \tau$, is exactly the statement that the exterior worldline of the surface is timelike and parametrised by its own proper time.
So $\tau$ is not merely the interior's time coordinate: it is the proper time of the surface read from either side, and the two descriptions share it.
That is what makes the next step's comparison of extrinsic curvatures a comparison of two tensors on one manifold rather than an equivocation.

Writing $f(r) = 1 - r_s/r$ and a dot for $d/d\tau$ for the length of this step and the next two, the time slot solves to

$$\frac{dT}{d\tau} = \frac{1}{f(R)}\sqrt{f(R) + \frac{1}{c^2}\left(\frac{dR}{d\tau}\right)^2},$$

the root taken positive so that the exterior time runs forward.
Nothing has been imposed yet beyond continuity: any $R(\tau)$ at all can be fitted with this $T(\tau)$.
The dynamics is in the second condition.

---

## Step 12. The extrinsic curvature from inside

The extrinsic curvature of $\Sigma$ is the rate at which its unit normal turns as it is carried along the surface,

$$K_{ij} = \nabla_\mu n_\nu\, e^\mu{}_i e^\nu{}_j,$$

with $e^\mu{}_i$ the three vectors tangent to $\Sigma$ and $n$ the outward unit normal, outward on both sides so that the two are comparable.

Inside, the surface is a level set of $\chi$, so the normal is proportional to $d\chi$, and the normalisation $g^{\chi\chi}n_\chi n_\chi = 1$ gives

$$n_\mu dx^\mu = a\,d\chi, \qquad n^\mu\partial_\mu = \frac{1}{a}\partial_\chi.$$

The tangent vectors are $\partial_\tau$, $\partial_\theta$ and $\partial_\phi$, and since the normal has only a $\chi$ component and does not depend on $\tau$, $\theta$ or $\phi$, each component of $K$ is a single Christoffel symbol from Step 5,

$$K_{ij} = -\Gamma^\chi{}_{ij}\,n_\chi = -a\,\Gamma^\chi{}_{ij}.$$

The three that matter are

$$K_{\tau\tau} = -a\,\Gamma^\chi{}_{\tau\tau} = 0, \qquad K_{\theta\theta} = -a\,\Gamma^\chi{}_{\theta\theta} = a\sin\chi_0\cos\chi_0, \qquad K_{\phi\phi} = K_{\theta\theta}\sin^2\theta,$$

the first because $\Gamma^\chi{}_{\tau\tau}$ vanishes, which Step 5 already flagged as the free fall of the dust.
The mixed components vanish with their Christoffel symbols.
In mixed form, which is the form that compares most cleanly,

$$K^\tau{}_\tau\big|_{\text{in}} = 0, \qquad K^\theta{}_\theta\big|_{\text{in}} = K^\phi{}_\phi\big|_{\text{in}} = \frac{\cos\chi_0}{a\sin\chi_0} = \frac{\cos\chi_0}{R}.$$

Note what has happened to $\chi_0$.
It entered as a label for which sphere of the three sphere is the edge of the star, and it has come out as the one number the extrinsic curvature of that edge depends on.

---

## Step 13. The extrinsic curvature from outside

Outside, the surface is not a coordinate surface, so the normal has to be built.
The tangent along the worldline is the four velocity

$$u^\mu = \left(\frac{dT}{d\tau},\ \frac{dR}{d\tau},\ 0,\ 0\right), \qquad u\cdot u = -c^2,$$

and the outward normal is the unit spacelike vector in the $(t,r)$ plane orthogonal to it,

$$n_\mu dx^\mu = -\frac{dR}{d\tau}\,dt + \frac{dT}{d\tau}\,dr.$$

That it is orthogonal is immediate, $n_\mu u^\mu = -\dot{R}\dot{T} + \dot{T}\dot{R} = 0$, and that it is a unit vector uses the normalisation of Step 11:

$$g^{\mu\nu}n_\mu n_\nu = -\frac{1}{c^2f}\left(\frac{dR}{d\tau}\right)^2 + f\left(\frac{dT}{d\tau}\right)^2 = 1.$$

The angular component of the extrinsic curvature is again a single Christoffel symbol, now the exterior's $\Gamma^r{}_{\theta\theta} = -rf$,

$$K_{\theta\theta}\big|_{\text{out}} = -\Gamma^r{}_{\theta\theta}n_r = R\,f(R)\frac{dT}{d\tau} = R\sqrt{f(R) + \frac{1}{c^2}\left(\frac{dR}{d\tau}\right)^2},$$

the last equality substituting the $dT/d\tau$ of Step 11.
In mixed form

$$K^\theta{}_\theta\big|_{\text{out}} = K^\phi{}_\phi\big|_{\text{out}} = \frac{1}{R}\sqrt{f(R)+\frac{1}{c^2}\left(\frac{dR}{d\tau}\right)^2}.$$

The remaining component is the normal projection of the four acceleration,

$$K^\tau{}_\tau\big|_{\text{out}} = -n_\mu A^\mu, \qquad A^\mu = u^\nu\nabla_\nu u^\mu,$$

which is zero exactly when the surface falls freely in the exterior geometry.

---

## Step 14. The second junction condition, and the one equation it gives

Now compare.

The $\tau\tau$ slot gives $n_\mu A^\mu = 0$, so the surface is a radial timelike geodesic of Schwarzschild.
Nothing has been assumed about the matter to reach that, beyond the pressurelessness that put $\Gamma^\chi{}_{\tau\tau} = 0$ into Step 5: the outermost shell of dust is in free fall in its own exterior field, which is Newton's statement that a shell falls as a particle would.

The $\theta\theta$ slot gives, cancelling one factor of $R = a\sin\chi_0$,

$$\cos\chi_0 = \sqrt{f(R) + \frac{1}{c^2}\left(\frac{dR}{d\tau}\right)^2}.$$

This is the equation with content.
Squaring it and writing $\dot{R}/c = \sin\chi_0\,\dot{a}$ in the chart derivative of Step 1,

$$\cos^2\chi_0 = 1 - \frac{r_s}{a\sin\chi_0} + \sin^2\chi_0\,\dot{a}^2.$$

Everything in it is known except the relation between $r_s$ and the interior.
Use the interior's own first integral from Step 8, $\dot{a}^2 = a_m/a - 1$:

$$\cos^2\chi_0 = 1 - \frac{r_s}{a\sin\chi_0} + \sin^2\chi_0\left(\frac{a_m}{a}-1\right) = \cos^2\chi_0 + \frac{1}{a}\left(a_m\sin^2\chi_0 - \frac{r_s}{\sin\chi_0}\right).$$

The $\cos^2\chi_0$ cancels from both sides, and what is left has to vanish at every value of $a$, which is to say at every moment of the collapse.
So

$$\boxed{\ r_s = a_m\sin^3\chi_0.\ }$$

Two things are worth noticing about how that came out.

It is a constraint on constants, not a differential equation, and it had to be: the two field equations have already fixed $a(\tau)$ up to the constant $a_m$, and the matching has only to say which exterior mass goes with which interior.
Had the $a$ dependence not cancelled, the junction conditions would have been inconsistent with the field equations and the construction would have failed.

And the matching used the interior's first integral, which is the pressureless condition integrated once.
That is the sense in which this is a dust result.
A ball with pressure has a different $a$ equation, its surface is not in free fall, and $K^\tau{}_\tau$ does not vanish on either side; matching it is the calculation `interior_schwarzschild.json` does at the other extreme, where the ball does not move at all.

---

## Step 15. What that equation says

**It is the mass of the ball.** Substituting $a_m = 8\pi G\rho a^3/(3c^2)$ from Step 8 and $R = a\sin\chi_0$,

$$r_s = \frac{8\pi G\rho a^3\sin^3\chi_0}{3c^2} = \frac{8\pi G\rho R^3}{3c^2}, \qquad\text{that is}\qquad M = \frac{4\pi}{3}\rho R^3.$$

The Schwarzschild mass of the exterior is the density times the Euclidean volume of a ball of the surface's areal radius, at every moment of the collapse, not only at the start.
Both sides of that equation change as the star falls, and they change together, because $\rho R^3$ is constant for the same reason $\rho a^3$ is.

**It is the compactness at release.** The star is released from rest at $a = a_m$, where its surface has areal radius $R_0 = a_m\sin\chi_0$, so

$$\sin^2\chi_0 = \frac{a_m\sin^3\chi_0}{a_m\sin\chi_0} = \frac{r_s}{R_0}.$$

The comoving angle of the surface, which looked like a piece of bookkeeping about the three sphere, is the compactness of the star at the moment it starts to fall.
A star released far outside its own Schwarzschild radius occupies a small polar cap of the three sphere; a star released close to it occupies a large one.
And $R_0 > r_s$ is $\sin^2\chi_0 < 1$, which is the bound $\chi_0 < \pi/2$ Step 3 got from the geometry of the three sphere.
The two readings of that bound are the same bound.

**It is a redshift, and an energy.** Since $\cos\chi_0 = \sqrt{1-\sin^2\chi_0}$,

$$\cos\chi_0 = \sqrt{1-\frac{r_s}{R_0}},$$

which is the redshift factor of a static observer sitting on the surface of the star before the collapse begins.
It is also, by Step 14, the conserved energy per unit rest mass of the surface's geodesic, $\tilde{E} = c^2\cos\chi_0$, which is what the energy of a particle released from rest at $R_0$ has to be.
A star with $\chi_0$ close to $\pi/2$ is one whose surface is already deeply redshifted before anything happens.

**And it is short of the mass the star is made of.** The proper volume of the interior at a given $\tau$ is not $\tfrac{4}{3}\pi R^3$ but the volume of a polar cap of a three sphere,

$$V = 4\pi a^3\int_0^{\chi_0}\sin^2\chi\,d\chi = 2\pi a^3\left(\chi_0 - \sin\chi_0\cos\chi_0\right),$$

so the proper mass, the mass got by counting the dust, is $M_p = \rho V$, and

$$\frac{M}{M_p} = \frac{\tfrac{4}{3}\pi\rho a^3\sin^3\chi_0}{2\pi\rho a^3\left(\chi_0-\sin\chi_0\cos\chi_0\right)} = \frac{2\sin^3\chi_0}{3\left(\chi_0-\sin\chi_0\cos\chi_0\right)} = 1 - \frac{3}{10}\chi_0^2 + O(\chi_0^4).$$

The ratio is less than one for every $\chi_0 > 0$, and the deficit is the gravitational binding energy.
Its leading term is worth writing out, using $\chi_0^2 \approx r_s/R_0$ from the second reading above:

$$\left(M_p - M\right)c^2 \approx \frac{3}{10}\frac{r_s}{R_0}Mc^2 = \frac{3}{5}\frac{GM^2}{R_0},$$

which is exactly the Newtonian binding energy of a uniform sphere.
The exterior mass is the rest mass of the dust less the energy it took to assemble it, and the relativistic calculation reproduces the textbook Newtonian number in the weak field limit without being told to.

---

## Step 16. The proper time to the singularity

From Step 8 the collapse runs from $\eta = 0$ to $\eta = \pi$, so it takes the proper time

$$\tau_s = \frac{a_m}{2c}\left(\pi + \sin\pi\right) = \frac{\pi a_m}{2c}.$$

In terms of the exterior's variables, using $a_m = R_0/\sin\chi_0$ and $\sin^2\chi_0 = r_s/R_0$ from Step 15,

$$\tau_s = \frac{\pi}{2c}\frac{R_0}{\sin\chi_0} = \frac{\pi}{2c}\frac{R_0^{3/2}}{r_s^{1/2}} = \frac{\pi}{2}\sqrt{\frac{R_0^3}{2GM}},$$

which is the Newtonian free fall time of a particle dropped from $R_0$ onto a point mass $M$, unchanged.
Written with the density instead, using $M = \tfrac{4}{3}\pi\rho_0R_0^3$,

$$\tau_s = \sqrt{\frac{3\pi}{32G\rho_0}},$$

depending on nothing but the density the ball started at.
A ball of the Sun's mean density, about $1.4\times10^3\ \mathrm{kg\,m^{-3}}$, falls in for about half an hour whatever its size; a ball of nuclear density falls in in a fraction of a millisecond.

This is the finite proper time of the history: the surface, and every other particle of the dust, reaches the singularity in a span of its own time that is not merely finite but unremarkable, with nothing locally to mark the horizon it passed on the way.
Step 18 is the other half of that sentence.

---

## Step 17. The horizon, from both sides

**Where the surface crosses.** The surface has areal radius $R = a\sin\chi_0$ and the horizon of the exterior sits at $r_s = a_m\sin^3\chi_0$, so the crossing $R = r_s$ happens at

$$a = a_m\sin^2\chi_0, \qquad\text{that is}\qquad 1+\cos\eta = 2\sin^2\chi_0, \qquad \cos\eta = -\cos2\chi_0,$$

$$\eta_H = \pi - 2\chi_0.$$

The proper time it happens at, and the proper time left afterwards, follow from the cycloid,

$$c\tau_H = \frac{a_m}{2}\left(\pi - 2\chi_0 + \sin2\chi_0\right), \qquad c\left(\tau_s - \tau_H\right) = \frac{a_m}{2}\left(2\chi_0 - \sin2\chi_0\right).$$

The second is the interesting one.
For a star that was large and diffuse when it was released, $\chi_0$ is small, and expanding with $a_m = r_s/\sin^3\chi_0$,

$$c\left(\tau_s-\tau_H\right) = \frac{r_s}{2\sin^3\chi_0}\left(2\chi_0 - \sin2\chi_0\right) = \frac{2r_s}{3} + \frac{r_s}{5}\chi_0^2 + O(\chi_0^4).$$

So whatever the star was, once its surface is at the Schwarzschild radius it has about $2r_s/3c = 4GM/3c^3$ of proper time left.
For ten solar masses that is about sixty six microseconds.
The limit $\chi_0 \to \pi/2$ is the other extreme, where $c(\tau_s-\tau_H) \to \pi a_m/2$ is the whole collapse, because a star released at $R_0 = r_s$ crosses the horizon at the moment it is released.

**The apparent horizon.** Inside the star the horizon is not at any fixed areal radius, and the thing that is locally defined is the trapped region: the set of spheres whose area decreases along their own outgoing light rays.
Step 8 said that a radial outgoing ray has $d\chi = d\eta$, so along one the areal radius changes at

$$\frac{d}{d\eta}\left[a(\eta)\sin\chi\right]_{\chi = \eta + \text{const}} = \frac{da}{d\eta}\sin\chi + a\cos\chi = \frac{a_m}{2}\left[\cos\chi + \cos\left(\chi+\eta\right)\right],$$

which vanishes on

$$\chi_{AH}(\eta) = \frac{\pi-\eta}{2}$$

and is negative for $\chi > \chi_{AH}$.
So the trapped region is the outer part of the star, not the inner, and its boundary sweeps inward as $\eta$ grows, reaching the centre exactly at $\eta = \pi$, the singularity.
Where does it first touch the surface? At $\chi_{AH} = \chi_0$, which is $\eta = \pi - 2\chi_0 = \eta_H$.
The apparent horizon reaches the surface of the star at precisely the moment the surface crosses the Schwarzschild radius, which is a consistency check on the whole construction: the interior and the exterior have to agree about when the star becomes trapped, and they do, having been told nothing about each other except the two junction conditions.

**The event horizon.** The event horizon is the boundary of what can still get out, so it is the outgoing null ray that reaches the surface exactly when the surface reaches $r_s$.
A ray reaching the surface a moment earlier leaves it while it is still outside $r_s$ and climbs away through the vacuum; one reaching it a moment later leaves from inside $r_s$ and does not.
With $\chi = \eta + \text{const}$ through $(\eta_H,\chi_0)$,

$$\chi_{EH}(\eta) = \chi_0 - \left(\eta_H - \eta\right) = \eta - \pi + 3\chi_0 .$$

It reaches the centre of the star at

$$\eta = \pi - 3\chi_0,$$

which is the moment the event horizon is born, and from then until $\eta_H$ it grows outward through the matter at the speed of light in comoving angle, meeting the apparent horizon at the surface.
For $\chi_0 > \pi/3$, which is $R_0 < \tfrac{4}{3}r_s$, that moment is before $\eta = 0$: such a star already has an event horizon inside it at the instant it is released.

Comparing the two through the collapse, at $\eta = \eta_H - \delta$,

$$\chi_{EH} = \chi_0 - \delta, \qquad \chi_{AH} = \chi_0 + \frac{\delta}{2},$$

so the trapped region $\chi > \chi_{AH}$ sits strictly inside the black hole region $\chi > \chi_{EH}$, which is the general theorem in this particular case.
The event horizon is a global object that knows the star's whole future, and it starts growing before anything local has gone wrong; the apparent horizon is local and lags behind it.

---

## Step 18. The view from outside, and why the star appears to freeze

Everything in Step 16 and Step 17 was written in $\tau$ or $\eta$, the time of somebody falling with the star.
The distant observer uses $t$, and Step 11 already has the conversion.
Substituting $f(R)\,dT/d\tau = \cos\chi_0$, which is what the matching of Step 14 made of it,

$$\frac{dT}{d\tau} = \frac{\cos\chi_0}{1-\dfrac{r_s}{R}}.$$

As $R \to r_s$ this diverges, and it diverges like $1/(R-r_s)$ while $R$ approaches $r_s$ linearly in $\tau$, so $T$ diverges logarithmically:

$$T \sim -\frac{r_s}{c}\ln\left(\tau_H - \tau\right).$$

The crossing that takes the surface a finite proper time to reach takes infinite Schwarzschild time.
The star never crosses, as far as the exterior chart is concerned, and this is the freezing the 1939 paper describes.
It is a statement about the chart and not about the star: $t$ is the proper time of a static observer far away, there is no static observer at $r_s$ for it to be the time of, and the surface itself passes $r_s$ at the perfectly ordinary $\eta_H$ of Step 17.

What the distant observer actually sees is light, so the sharper statement is about the redshift.
Let $u = ct - r_*$ be the retarded time, with $r_* = r + r_s\ln\left|r/r_s - 1\right|$ the tortoise coordinate, so that an outgoing radial ray has $u$ constant.
A ray leaving the surface at proper time $\tau$ carries

$$\frac{du}{d\tau} = c\frac{dT}{d\tau} - \frac{1}{f(R)}\frac{dR}{d\tau} = \frac{c\cos\chi_0 - dR/d\tau}{f(R)}.$$

Near the crossing $dR/d\tau \to -c\cos\chi_0$, from the geodesic first integral of Step 14 at $f = 0$, so the numerator tends to $2c\cos\chi_0$, while

$$f(R) = \frac{R-r_s}{R} \approx \frac{c\cos\chi_0\left(\tau_H-\tau\right)}{r_s}.$$

Therefore

$$\frac{du}{d\tau} \approx \frac{2r_s}{\tau_H-\tau}, \qquad \tau_H - \tau \propto e^{-u/2r_s}.$$

The emitter's clock and the receiver's are related by exactly that factor, so a spectral line of proper frequency $\omega_0$ arrives at

$$\omega_\infty = \omega_0\frac{d\tau}{du} \propto \omega_0\,e^{-u/2r_s},$$

falling by a factor of $e$ every

$$\Delta t = \frac{2r_s}{c} = \frac{4GM}{c^3}$$

of the distant observer's time.
For a solar mass that is about twenty microseconds, and for ten solar masses about two hundred.
The received flux falls faster still, since the photons arrive rarer as well as redder and the cone of directions that escapes closes down as the surface approaches the horizon; this entry computes only the frequency factor, which is enough to make the point.
The star does not so much freeze as vanish: after a few of those e foldings there is nothing left to see, and the eternal hovering image the chart suggests is far too faint to detect.
An infalling observer crossing with the surface, meanwhile, measures the ordinary local physics of Step 9 and has the sixty six microseconds of Step 17 left to think about it.

---

## Step 19. The curvature does not match, and does not have to

The junction conditions asked for the continuity of two things, the induced metric and the extrinsic curvature.
They did not ask for the continuity of the curvature tensors, and the curvature tensors are not continuous.

The Ricci tensor obviously jumps, since the density jumps from $\rho$ to zero at the surface, and the field equations tie the Ricci tensor to the density pointwise.
The Weyl tensor jumps too, and less obviously, since Step 9 computed it to be identically zero inside and Step 10 has it nonzero outside.
The cleanest statement of both is the Kretschmann scalar.
Just inside the surface, from Step 9 at $a = R/\sin\chi_0$ and $a_m = r_s/\sin^3\chi_0$,

$$K_{\text{in}} = \frac{15a_m^2}{a^6} = \frac{15\,r_s^2}{R^6},$$

and just outside, from Step 10 at $r = R$,

$$K_{\text{out}} = \frac{12\,r_s^2}{R^6}.$$

So crossing the surface outward multiplies the Kretschmann scalar by $4/5$, everywhere and at every moment, which is as clean a way as there is to say that the two curvatures are of entirely different kinds and happen to be of comparable size.

This is allowed because the Riemann tensor at a surface is not determined by the induced metric and the extrinsic curvature alone.
The Gauss-Codazzi relations fix the components with two or fewer normal indices from $h_{ij}$ and $K_{ij}$, and those components do match; the component with two normal indices involves the normal derivative of $K_{ij}$, which jumps whenever the source does.
A discontinuous density gives a discontinuous curvature and no surface layer, which is exactly the physical situation: a sharp edge to the star, and nothing painted on the edge.

---

## Step 20. The geodesic equations

With the dots the chart velocities of Step 1, the geodesic equation $\ddot{x}^\mu + \Gamma^\mu{}_{\nu\rho}\dot{x}^\nu\dot{x}^\rho = 0$ and the Christoffel symbols of Step 5 give, for the interior,

$$\ddot{\tau} + a\dot{a}\left(\dot{\chi}^2 + \sin^2\chi\,\dot{\theta}^2 + \sin^2\chi\sin^2\theta\,\dot{\phi}^2\right) = 0,$$

$$\ddot{\chi} + \frac{2\dot{a}}{a}\dot{\tau}\dot{\chi} - \sin\chi\cos\chi\left(\dot{\theta}^2+\sin^2\theta\,\dot{\phi}^2\right) = 0,$$

$$\ddot{\theta} + \frac{2\dot{a}}{a}\dot{\tau}\dot{\theta} + 2\cot\chi\,\dot{\chi}\dot{\theta} - \sin\theta\cos\theta\,\dot{\phi}^2 = 0,$$

$$\ddot{\phi} + \frac{2\dot{a}}{a}\dot{\tau}\dot{\phi} + 2\cot\chi\,\dot{\chi}\dot{\phi} + 2\cot\theta\,\dot{\theta}\dot{\phi} = 0,$$

the factors of two on the mixed terms being the two orderings of each mixed Christoffel index pair.
The exterior's four are Schwarzschild's and the entry publishes them unchanged.

Two solutions of the interior equations are used above.
A comoving particle, $\dot{\chi} = \dot{\theta} = \dot{\phi} = 0$, satisfies all four with $\ddot{\tau} = 0$, so the dust is in free fall and $\tau$ is its proper time, which is Step 5's remark and Step 12's $K_{\tau\tau} = 0$.
A radial null ray has $\dot{\theta} = \dot{\phi} = 0$ and $a\dot{\chi} = \pm\dot{\tau}$ from the line element, which integrates to the $\chi = \pm\eta + \text{const}$ of Step 8, and Step 17 is built on it.

---

## Step 21. Every published expression is dimensionally consistent

The interior's chart coordinates are $\left(c\tau, \chi, \theta, \phi\right)$, one length and three dimensionless angles, and the parameters are $a$ and $a_m$, both lengths, and $\chi_0$, an angle.
A component then carries

$$[g_{\mu\nu}] = \frac{L^2}{[x^\mu][x^\nu]},$$

with an upper index contributing $[x^\mu]/L$ and a lower one $L/[x^\mu]$ on top of the $1$, $L^{-1}$, $L^{-2}$ or $L^{-4}$ the field itself carries.
Since the angles are dimensionless, each angular index contributes a whole factor of $L$, which is what makes the bookkeeping here less trivial than in a Cartesian chart.

The derivative convention does the rest.
A dot is $d/d(c\tau)$, so $\dot{a}$ is dimensionless and $\ddot{a}$ carries $L^{-1}$.

Five samples.

$g_{\theta\theta} = a^2\sin^2\chi$ carries $L^2$, which is $L^2/([\theta][\theta]) = L^2$.

$\Gamma^\tau{}_{\chi\chi} = a\dot{a}$ carries $L$, and the requirement is $L^{-1}\cdot\left([\tau]/L\right)\cdot\left(L/[\chi]\right)^2 = L^{-1}\cdot 1\cdot L^2 = L$.
A published $\Gamma^\chi{}_{\tau\chi} = \dot{a}/a$ carries $L^{-1}$, and its requirement is $L^{-1}\cdot L^{-1}\cdot 1\cdot L = L^{-1}$.
Both work only because $\dot{a}$ is the chart derivative; on the reading $\dot{a} = da/d\tau$ each would be one factor of $c$ out.

$R^\chi{}_{\theta\chi\theta} = \left(\dot{a}^2+1\right)\sin^2\chi$ is dimensionless, and so is its requirement, $L^{-2}\cdot L^{-1}\cdot L\cdot L\cdot L = 1$.
The two terms of the bracket balance because $\dot{a}$ is dimensionless, which is the whole reason this chart puts the length in $a$ rather than in $\chi$.

$K = 12\left(a^2\ddot{a}^2 + \left(\dot{a}^2+1\right)^2\right)/a^4$ carries $\left(L^2L^{-2} + 1\right)/L^4 = L^{-4}$, in both terms.

The geodesic equations are measured against their own second derivatives.
In the $\tau$ equation $\ddot{\tau}$ carries $L/\lambda^2$, since the chart coordinate is a length, and the term $a\dot{a}\sin^2\chi\,\dot{\theta}^2$ carries

$$L\cdot 1\cdot\frac{1}{\lambda^2} = \frac{L}{\lambda^2},$$

because $\dot{\theta}$ is an angle per affine parameter.
In the $\chi$ equation $\ddot{\chi}$ carries $\lambda^{-2}$ and the mixed term $2\left(\dot{a}/a\right)\dot{\tau}\dot{\chi}$ carries $L^{-1}\cdot\left(L/\lambda\right)\cdot\left(1/\lambda\right) = \lambda^{-2}$, which matches only on the reading that $\dot{\tau}$ is $d(c\tau)/d\lambda$.

The exterior's balance is Schwarzschild's and holds for the same reasons; `verify_metrics.py` checks both systems term by term.

---

## Step 22. What the entry publishes, and the declarations the checker needs

The interior publishes the line element, four metric components and four inverse ones, eighteen Christoffel symbols in each of two variants, twenty four Riemann components in each of two variants, four Ricci components in each of three variants, the Ricci scalar, the Kretschmann scalar, four Einstein components in each of three variants, two empty Weyl blocks and four geodesic equations.
The exterior publishes Schwarzschild's spherical chart entire, as Step 10 lists it.

For the checker to read them, `DIMENSIONS` in `verify_metrics.py` needs one line each:

    ("oppenheimer_snyder", "interior_comoving"): {
        "\\tau": "T", "\\chi": "1", "\\theta": "1", "\\phi": "1",
        "a": "L", "\\chi_0": "1", "a_m": "L",
    },
    ("oppenheimer_snyder", "exterior_schwarzschild"): {
        "t": "T", "r": "L", "\\theta": "1", "\\phi": "1", "r_s": "L",
    },

The declaration `"\\tau": "T"` is what tells the checker that $\tau$ is a time and so that the chart coordinate is $c\tau$; every factor of $c$ in the comparison follows from it, as does the reading of the dot on $a$ as a chart derivative.
Declaring `"a": "L"` with `"\\chi": "1"` is the choice of Step 3, the opposite of the one `frw.json` declares, and it is the declaration that makes $\dot{a}$ dimensionless and $R = a\sin\chi_0$ a length.
$\chi_0$ and $a_m$ are declared although no published component contains them, because the entry states the matching in them and the checker requires the declarations and the parameter list to agree exactly.
Neither appears in a component value on purpose, quite apart from taste: the reader in `verify_metrics.py` turns `\chi` into a name before it looks at what follows, so a `\chi_0` inside a published value would be read as $\chi$ times an unknown, and rejected.

The entry needs no line in `PARAMETER_RELATIONS`.
It could have wanted one, and the reason it does not is Step 7's last paragraph: the published components are printed for a free $a(\tau)$, so every one of them is an identity that holds for any line element of that form, and the dust solution, the cycloid and the matching are stated in the entry's prose and derived here rather than substituted into a curvature block.
The relation $r_s = a_m\sin^3\chi_0$ ties a parameter of one system to a parameter of another, which `PARAMETER_RELATIONS` has no way to express and no need to, since neither symbol appears in a component.

Running

    /tmp/mfs-venv-oppenheimer/bin/python _tools/derivations/verify_metrics.py \
        --system oppenheimer_snyder/interior_comoving \
        --system oppenheimer_snyder/exterior_schwarzschild

reports no disagreements and no dimensional failures for either system.
