# Voice rewrites of the captions and conventions

Every paragraph of every History section, every paragraph of every conventions section, and every caption, note, restriction band and sentence stated in place of a drawing on the spacetimes page was read against ~/VOICE.md and rewritten here where it strayed from it.
The captain added one rule on 28 September 2026: a vague figurative phrase never stands in for a precise claim, as "so the spacetime does not stand on its own" stood for "so anti-de Sitter space is not globally hyperbolic".
Each entry names the source the text is written from, then gives the text before and after, one sentence to a line; joining the lines of a block with single spaces gives the text exactly.
Every `$...$` of every entry is unchanged character for character and in order, and so is every number written in digits outside it and every citation.
Of 675 distinct texts, 274 are rewritten and 401 stand as they were.

## Alcubierre Warp Drive

### history: `alcubierre/history[1]`

Source: `MFS/assets/data/metrics/alcubierre.json`, `history`, paragraph 2.

Before:

> Locally nothing moves faster than light, since the ship sits at rest in its own little patch, yet the bubble itself glides through spacetime at any speed one likes.
> Special relativity forbids outrunning a light beam locally; it says nothing about how fast a region of space may be made to move.
> The catch is in the source.
> To bend spacetime this way the Einstein equations demand a region of negative energy density, exotic matter that violates the energy conditions every known form of matter obeys.

After:

> Locally nothing moves faster than light, since the ship sits at rest in its own little patch, yet the bubble itself glides through spacetime at any speed one likes.
> Special relativity forbids outrunning a light beam locally; it says nothing about how fast a region of space may be made to move.
> The difficulty is the source: to bend spacetime this way the Einstein equations demand a region of negative energy density, exotic matter that violates the energy conditions every known form of matter obeys.

### history: `alcubierre/history[2]`

Source: `MFS/assets/data/metrics/alcubierre.json`, `history`, paragraph 3.

Before:

> Worse, the amount is not a technicality.
> In 1997 Michael Pfenning and Larry Ford applied quantum inequalities to the bubble and found that its walls must be only a few hundred Planck lengths thick, and that the total negative energy required exceeds, by many orders of magnitude, anything physically conceivable, comparable in scale to all the mass and energy of the visible universe [pfenning1997].
> So the warp drive lives as a thought experiment rather than an engineering proposal, but a remarkably productive one.

After:

> The negative energy needed is enormous.
> In 1997 Michael Pfenning and Larry Ford applied quantum inequalities to the bubble and found that its walls must be only a few hundred Planck lengths thick, and that the total negative energy required exceeds, by many orders of magnitude, anything physically conceivable, comparable in scale to all the mass and energy of the visible universe [pfenning1997].
> The warp drive therefore remains a thought experiment, though a remarkably productive one.

### history: `alcubierre/history[3]`

Source: `MFS/assets/data/metrics/alcubierre.json`, `history`, paragraph 4.

Before:

> It sharpened the question of how far the energy conditions can be pushed, and it spawned a family of variants.
> In 1999 Chris Van Den Broeck showed that a clever change of the bubble's shape, a tiny throat opening into a large interior, could cut the required negative energy from the scale of a universe down to a few solar masses [vandenbroeck1999].
> José Natário demonstrated in 2002 that the picture of contraction ahead and expansion behind is not even essential, constructing a warp drive with no volume expansion at all [natario2002].
> Yet the central obstruction proved stubborn: in 2004 Francisco Lobo and Matt Visser showed that the violations of the energy conditions are not merely an extravagance of high speed but persist even for arbitrarily slow bubbles [lobo2004].

After:

> The drive sharpened the question of how far the energy conditions can be pushed, and it spawned a family of variants.
> In 1999 Chris Van Den Broeck showed that a clever change of the bubble's shape, a tiny throat opening into a large interior, could cut the required negative energy from the scale of a universe down to a few solar masses [vandenbroeck1999].
> José Natário showed in 2002 that the ship needs no contraction ahead and expansion behind, constructing a warp drive with no volume expansion [natario2002].
> The central obstruction held: in 2004 Francisco Lobo and Matt Visser showed that the energy conditions fail for bubbles of every speed, however slow [lobo2004].

### history: `alcubierre/history[4]`

Source: `MFS/assets/data/metrics/alcubierre.json`, `history`, paragraph 5.

Before:

> Alexey Bobrick, Gianni Martire, and Erik Lentz took up the question of the source again in 2021.
> Bobrick and Martire gave a general framework that encompasses every warp drive metric yet proposed and clarified which of their features are negotiable, describing physical, subluminal warp drives free of the worst pathologies [bobrick2021]; Lentz exhibited soliton solutions that can in principle be sourced by purely positive energy densities drawn from a conducting plasma and classical electromagnetic fields, the first warp geometries built from familiar matter rather than the forbidden kind [lentz2021].
> Whether anyone could ever build such a drive is another matter.
> But as a clean demonstration of what general relativity permits, and of exactly where it digs in its heels, the warp drive has few equals.

After:

> Alexey Bobrick, Gianni Martire, and Erik Lentz took up the question of the source again in 2021.
> Bobrick and Martire gave a general framework that encompasses every warp drive metric yet proposed and clarified which of their features are negotiable, describing physical, subluminal warp drives free of the worst pathologies [bobrick2021]; Lentz exhibited soliton solutions that can in principle be sourced by purely positive energy densities drawn from a conducting plasma and classical electromagnetic fields, the first warp geometries built from familiar matter rather than the forbidden kind [lentz2021].
> Whether anyone could build such a drive remains open.
> General relativity permits the geometry, and the energy conditions mark where it resists.

### conventions: `alcubierre/convention[0]`

Source: `MFS/assets/data/metrics/alcubierre.json`, `convention`, paragraph 1.

Before:

> Coordinates are $(t, x, y, z)$ with $t$ carrying dimensions of time and $x$, $y$ and $z$ lengths, and factors of $c$ are kept explicit.
> The line element is the standard warp drive form, which in the language of a $3+1$ split has unit lapse, flat spatial slices and a single shift component $c\,v_s f$ along $x$, so the whole of the geometry sits in the shift and the spacetime is flat wherever the shape function is constant.
> Both parameters are dimensionless.
> $v_s$ is the velocity of the bubble in units of $c$, so its centre follows $x = x_s(t)$ with $dx_s/dt = c\,v_s$, and $f$ is the shape function, left an arbitrary function of all four coordinates rather than fixed to Alcubierre's profile of hyperbolic tangents, so every component holds for any bubble of this form.

After:

> We use coordinates $(t, x, y, z)$, with $t$ carrying dimensions of time and $x$, $y$, and $z$ lengths, and keep every factor of $c$ explicit.
> The line element is the standard warp drive form: in the language of a $3+1$ split it has unit lapse, flat spatial slices, and a single shift component $c\,v_s f$ along $x$, so the geometry sits wholly in the shift, and the spacetime is flat wherever the shape function is constant.
> Both parameters are dimensionless.
> The first, $v_s$, is the velocity of the bubble in units of $c$, so its centre follows $x = x_s(t)$ with $dx_s/dt = c\,v_s$; the second, $f$, is the shape function, which we leave an arbitrary function of all four coordinates rather than fix to Alcubierre's profile of hyperbolic tangents, so every component holds for any bubble of this form.

### conventions: `alcubierre/convention[1]`

Source: `MFS/assets/data/metrics/alcubierre.json`, `convention`, paragraph 2.

Before:

> Components are taken in the chart $x^0 = ct$, where a chart component of the metric carries no dimensions.
> Every $\partial_t$ is the derivative along that chart coordinate, $\partial_t = c^{-1}\partial/\partial t$, and the prime on $v_s$ is the same derivative, $v_s' = dv_s/d(ct)$, so each derivative carries one inverse length per order.
> The dots in the geodesic equations are velocities of that same chart, so $\dot{t}$ means $d(ct)/d\lambda$, and the equations are $\ddot{x}^\mu + \Gamma^\mu{}_{\nu\rho}\dot{x}^\nu\dot{x}^\rho = 0$ with $\Gamma$ the Christoffel symbols of that same chart.
> The signature is $(-,+,+,+)$, the Riemann tensor is $R^\mu{}_{\nu\rho\sigma} = \partial_\rho\Gamma^\mu{}_{\nu\sigma} - \partial_\sigma\Gamma^\mu{}_{\nu\rho} + \Gamma^\mu{}_{\rho\lambda}\Gamma^\lambda{}_{\nu\sigma} - \Gamma^\mu{}_{\sigma\lambda}\Gamma^\lambda{}_{\nu\rho}$, and the Ricci tensor is the standard contraction $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$.

After:

> We take components in the chart $x^0 = ct$, where a chart component of the metric carries no dimensions.
> Every $\partial_t$ is the derivative along that chart coordinate, $\partial_t = c^{-1}\partial/\partial t$, and a prime on $v_s$ is the same derivative, $v_s' = dv_s/d(ct)$, so each order of derivative carries one inverse length.
> The dots in the geodesic equations are velocities in the same chart, so $\dot{t}$ means $d(ct)/d\lambda$, and the equations read $\ddot{x}^\mu + \Gamma^\mu{}_{\nu\rho}\dot{x}^\nu\dot{x}^\rho = 0$, with $\Gamma$ the Christoffel symbols of that chart.
> We work in the signature $(-,+,+,+)$, with the Riemann tensor $R^\mu{}_{\nu\rho\sigma} = \partial_\rho\Gamma^\mu{}_{\nu\sigma} - \partial_\sigma\Gamma^\mu{}_{\nu\rho} + \Gamma^\mu{}_{\rho\lambda}\Gamma^\lambda{}_{\nu\sigma} - \Gamma^\mu{}_{\sigma\lambda}\Gamma^\lambda{}_{\nu\rho}$ and the Ricci tensor its standard contraction, $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$.

### conventions: `alcubierre/convention[2]`

Source: `MFS/assets/data/metrics/alcubierre.json`, `convention`, paragraph 3.

Before:

> One component carries the physics of the drive.
> The observer at rest in the slicing, the Eulerian observer with $n_\mu = (-1, 0, 0, 0)$ and $n^\mu = (1, v_s f, 0, 0)$, measures the energy density $\rho = T_{\mu\nu}n^\mu n^\nu = \dfrac{c^4}{8\pi G}G^{tt} = -\dfrac{c^4}{32\pi G}v_s^2\left[\left(\partial_y f\right)^2 + \left(\partial_z f\right)^2\right]$, which is $G^{tt}$ read through the field equations, since the lapse is $1$ and $n_\mu$ has only a time component.
> It is minus a sum of squares, so it is negative wherever the wall of the bubble has any transverse gradient at all, and no choice of shape function and no choice of speed can turn it positive; it grows as the square of the speed and as the square of the steepness of the wall, which is why a thin fast bubble costs so much.
> That is the standing objection to the drive, written as a component: the weak energy condition fails everywhere the wall is.
> For Alcubierre's radial profile the transverse gradient is $\left(\partial_y f\right)^2 + \left(\partial_z f\right)^2 = \dfrac{y^2 + z^2}{r_s^2}\left(\dfrac{df}{dr_s}\right)^2$, and the density is his expression of 1994.
> The expansion of those same observers is $\theta = \nabla_\mu n^\mu = c\,v_s\,\partial_x f$, negative ahead of the ship and positive behind it, which is the contraction of space in front and the expansion behind.

After:

> We turn now to the energy density, which carries the physics of the drive.
> The observer at rest in the slicing, the Eulerian observer with $n_\mu = (-1, 0, 0, 0)$ and $n^\mu = (1, v_s f, 0, 0)$, measures the energy density $\rho = T_{\mu\nu}n^\mu n^\nu = \dfrac{c^4}{8\pi G}G^{tt} = -\dfrac{c^4}{32\pi G}v_s^2\left[\left(\partial_y f\right)^2 + \left(\partial_z f\right)^2\right]$, by the field equations a multiple of $G^{tt}$ alone, since the lapse is $1$ and $n_\mu$ has only a time component.
> The density is minus a sum of squares, negative wherever the wall of the bubble has any transverse gradient, and no choice of shape function or of speed turns it positive; it grows as the square of the speed and as the square of the steepness of the wall, which is why a thin, fast bubble costs so much.
> The weak energy condition therefore fails everywhere in the wall, and this is the standing objection to the drive.
> For Alcubierre's radial profile the transverse gradient is $\left(\partial_y f\right)^2 + \left(\partial_z f\right)^2 = \dfrac{y^2 + z^2}{r_s^2}\left(\dfrac{df}{dr_s}\right)^2$, and the density is Alcubierre's own expression of 1994.
> The expansion of the same observers is $\theta = \nabla_\mu n^\mu = c\,v_s\,\partial_x f$, negative ahead of the ship and positive behind it: space contracts in front of the ship and expands behind.

### conventions: `alcubierre/convention[3]`

Source: `MFS/assets/data/metrics/alcubierre.json`, `convention`, paragraph 4.

Before:

> This is not a vacuum, so the Weyl tensor is not the Riemann tensor but the Riemann tensor with its traces removed: $C_{yzyz}$ is nonzero while $R_{yzyz}$ vanishes, and thirty six slots of $C_{\mu\nu\rho\sigma}$ are like that.
> The Kretschmann scalar is a polynomial in $v_s$ and in the first and second derivatives of $f$, so it stays finite wherever the shape function is smooth: the bubble carries no curvature singularity anywhere, and what stands in its way is the source it demands rather than a place where the geometry breaks.
> The bubble frame, in which the ship sits at rest at the origin, is reached by $x \to x - x_s(t)$, and the line element takes the same form there with $f$ replaced by $f - 1$, so the bubble frame is not a second chart: every component in the chart $x^0 = ct$ is already a component of that frame under the substitution.

After:

> The spacetime is not a vacuum, and its Weyl tensor is the Riemann tensor with the traces removed, traces that here do not vanish: $C_{yzyz}$ is nonzero while $R_{yzyz}$ vanishes, and thirty six slots of $C_{\mu\nu\rho\sigma}$ behave the same way.
> The Kretschmann scalar is a polynomial in $v_s$ and in the first and second derivatives of $f$, so it stays finite wherever the shape function is smooth; the bubble has no curvature singularity, and the one obstacle to it is the source it demands, since nowhere does the geometry break.
> The bubble frame, in which the ship sits at rest at the origin, follows from $x \to x - x_s(t)$, and there the line element takes the same form with $f$ replaced by $f - 1$; every component in the chart $x^0 = ct$ is therefore already a component of that frame under the substitution, and the bubble frame needs no chart of its own.

### caption of a figure in three dimensions: `diagrams/alcubierre/projections/cartesian/bubble.caption[0]`

Source: `_tools/derivations/projections.py` line 577, `CAPTIONS ("alcubierre", "cartesian", "bubble")`.

Before:

> This is the slice $z = 0$ of $t$, $x$ and $y$ through a bubble moving at twice the speed of light along $x$, with $t$ up, $ct$ and $x$ drawn at one scale, at the moment $t = 0$ when the bubble is centred on $x = 0$.
> Far from the bubble $f = 0$ and the cones stand upright, as Minkowski's do.
> Inside it $f$ is close to 1, and the shift $v_sf$ tips every cone forward along $x$ so far that the vertical lies outside it: nothing inside can stay at fixed $x$.
> The tilt is along $x$ everywhere, since the shift points along $x$, and it depends only on the distance from the centre, as $f$ does, so the cones tip over in a ball about the centre.

After:

> This is the slice $z = 0$ of $t$, $x$, and $y$ through a bubble moving at twice the speed of light along $x$, with $t$ up, $ct$ and $x$ drawn at one scale, at the moment $t = 0$ when the bubble is centred on $x = 0$.
> Far from the bubble $f = 0$ and the cones stand upright, as Minkowski's do.
> Inside it $f$ is close to 1, and the shift $v_sf$ tips every cone forward along $x$ so far that the vertical lies outside it: nothing inside can stay at fixed $x$.
> The tilt is along $x$ everywhere, since the shift points along $x$, and it depends only on the distance from the centre, as $f$ does, so the cones tip over in a ball about the centre.

### embedding diagram caption: `embedding/alcubierre/plane.caption[0]`

Source: `_tools/derivations/embedding.py` line 3674, `CAPTIONS ("alcubierre", "plane")`.

Before:

> This is the expansion $\theta$ of the observers who ride the slices of Alcubierre's warp drive, drawn as a height over the plane $z = 0$ of the ship's path at the moment $t = 0$: the height is $\theta$, a height of $R$ for an expansion of $4c/R$, and no surface of the spacetime, whose slices are flat.
> The observers are carried along $x$ at $v_sf$ times the speed of light, faster than light inside the circle $v_sf = 1$, and a small volume of them changes at the rate $\theta = c\,v_s\,\partial_xf = c\,v_s\,\frac{x - x_s}{r_s}\frac{df}{dr_s}$.

After:

> This is the expansion $\theta$ of the observers who ride the slices of Alcubierre's warp drive, drawn as a height over the plane $z = 0$ of the ship's path at the moment $t = 0$: the height stands for $\theta$ alone, a height of $R$ for an expansion of $4c/R$, and the slices themselves are flat.
> The observers are carried along $x$ at $v_sf$ times the speed of light, faster than light inside the circle $v_sf = 1$, and a small volume of them changes at the rate $\theta = c\,v_s\,\partial_xf = c\,v_s\,\frac{x - x_s}{r_s}\frac{df}{dr_s}$.

### embedding diagram note: `embedding/alcubierre/plane.input`

Source: `_tools/derivations/embedding.py` line 2850, `in alcubierre() input=`.

Before:

> $v_s = 2$, and Alcubierre's own profile, $f = [\tanh\sigma(r_s + R) - \tanh\sigma(r_s - R)]/(2\tanh\sigma R)$ with $R = 1$ and $\sigma = 4$, as the spacetime diagram declares.

After:

> $v_s = 2$, and Alcubierre's own profile, $f = [\tanh\sigma(r_s + R) - \tanh\sigma(r_s - R)]/(2\tanh\sigma R)$ with $R = 1$ and $\sigma = 4$, as in the spacetime diagram.

## Anti-de Sitter

### history: `anti_de_sitter/history[0]`

Source: `MFS/assets/data/metrics/anti_de_sitter.json`, `history`, paragraph 1.

Before:

> Anti-de Sitter space was born as the other branch of the same idea.
> When Willem de Sitter solved Einstein's amended equations for an empty universe in 1917, the cosmological constant $\Lambda$ was free to take either sign [desitter1917, einstein1917].
> Positive $\Lambda$ gave the exponentially expanding world that bears his name; negative $\Lambda$ gave its mirror image, a maximally symmetric vacuum with constant negative curvature.
> For decades this second solution had no champion and no name beyond the obvious one: the anti-de Sitter universe, de Sitter's geometry run in reverse.

After:

> Anti-de Sitter space was born as the other branch of the same idea.
> When Willem de Sitter solved Einstein's amended equations for an empty universe in 1917, the cosmological constant $\Lambda$ was free to take either sign [desitter1917, einstein1917].
> Positive $\Lambda$ gave the exponentially expanding world that bears his name; negative $\Lambda$ gave its mirror image, a maximally symmetric vacuum with constant negative curvature.
> For decades this second solution had no champion and no name beyond the obvious one: the anti-de Sitter universe, de Sitter's geometry with its curvature reversed in sign.

### history: `anti_de_sitter/history[1]`

Source: `MFS/assets/data/metrics/anti_de_sitter.json`, `history`, paragraph 2.

Before:

> As a cosmology it was a nonstarter.
> Anti-de Sitter space is static, not expanding, so it cannot describe a universe like ours.
> Its spatial sections are hyperbolic, and in its natural form it is riddled with closed timelike curves, and one must pass to the universal cover to be rid of them [hawking1973].
> Strangest of all, it has a timelike boundary at spatial infinity: a light ray can reach the edge of the universe and return in finite time, so the spacetime does not stand on its own.
> To evolve anything inside it, one must say what happens at that boundary [hawking1973].
> For sixty years this was a catalogue of pathologies more than a subject.

After:

> As a cosmology it was a nonstarter.
> Anti-de Sitter space is static, so it cannot describe an expanding universe like ours.
> Its spatial sections are hyperbolic, and in its natural form it is riddled with closed timelike curves, and one must pass to the universal cover to be rid of them [hawking1973].
> Strangest of all, it has a timelike boundary at spatial infinity: a light ray can reach the edge of the universe and return in finite time, so anti-de Sitter space is not globally hyperbolic.
> To evolve anything inside it, one must specify what happens at that boundary [hawking1973].
> For sixty years it was studied mostly for these pathologies.

### history: `anti_de_sitter/history[3]`

Source: `MFS/assets/data/metrics/anti_de_sitter.json`, `history`, paragraph 4.

Before:

> Then in 1997 Juan Maldacena turned the pathology into the point [maldacena1998].
> That troublesome timelike boundary, he argued, is where a quantum field theory lives, and the gravitational physics of the anti-de Sitter interior is exactly equivalent to a conformal field theory on its edge, which carries one dimension fewer; within months Gubser, Klebanov, and Polyakov, and independently Witten, sharpened the proposal into a precise dictionary between the two sides [gubser1998, witten1998].
> This is the AdS/CFT correspondence, the sharpest realization of the holographic principle anticipated by 't Hooft and Susskind [thooft1993, susskind1995], and it made anti-de Sitter space one of the most studied geometries in theoretical physics, not because the universe is anti-de Sitter, but because gravity there can be read off a boundary that has no gravity at all.
> Edward Witten used it in 1998 to relate the thermodynamics of the gauge theory on the boundary to that of Schwarzschild black holes in anti-de Sitter space, so that confinement and the mass gap of the gauge theory are coded in classical geometry [witten1998thermal].
> Shinsei Ryu and Tadashi Takayanagi proposed in 2006 that the entanglement entropy of a region of the boundary is the area of a minimal surface in the interior, by analogy with the Bekenstein-Hawking formula for the entropy of a black hole [ryu2006].

After:

> Then in 1997 Juan Maldacena found a use for the boundary [maldacena1998].
> That troublesome timelike boundary, he argued, is where a quantum field theory lives, and the gravitational physics of the anti-de Sitter interior is exactly equivalent to a conformal field theory on its edge, which carries one dimension fewer; within months Steven Gubser, Igor Klebanov, and Alexander Polyakov, and independently Edward Witten, sharpened the proposal into a precise dictionary between the two sides [gubser1998, witten1998].
> This is the AdS/CFT correspondence, the sharpest realization of the holographic principle anticipated by Gerard 't Hooft and Leonard Susskind [thooft1993, susskind1995], and it made anti-de Sitter space one of the most studied geometries in theoretical physics, though the universe is not anti-de Sitter, because there gravity is encoded in a theory on the boundary that has no gravity.
> Witten used it in 1998 to relate the thermodynamics of the gauge theory on the boundary to that of Schwarzschild black holes in anti-de Sitter space, so that confinement and the mass gap of the gauge theory are coded in classical geometry [witten1998thermal].
> Shinsei Ryu and Tadashi Takayanagi proposed in 2006 that the entanglement entropy of a region of the boundary is the area of a minimal surface in the interior, by analogy with the Bekenstein-Hawking formula for the entropy of a black hole [ryu2006].

### conventions: `anti_de_sitter/convention[0]`

Source: `MFS/assets/data/metrics/anti_de_sitter.json`, `convention`, paragraph 1.

Before:

> Coordinates are $(t, r, \theta, \phi)$ in the global chart, with $t$ carrying dimensions of time and $r$ the areal radius, and $(t, x, y, z)$ in the Poincaré patch, with $z$ the direction off the boundary.
> Factors of $c$ and $G$ are kept explicit, and the one parameter is the anti-de Sitter radius $L$, a length, which fixes the cosmological constant through $\Lambda = -3/L^2$; no mass enters, so $G$ appears in no component.
> Components are taken in the chart $x^0 = ct$, where the metric is dimensionless and no component carries a factor of $c$, and the dots in the geodesic equations are derivatives of that same chart with respect to an affine parameter, so $\dot{t}$ means $d(ct)/d\lambda$ and every term of every equation carries the dimensions of its left hand side.
> The geodesic equations are $\ddot{x}^\mu + \Gamma^\mu{}_{\nu\rho}\dot{x}^\nu\dot{x}^\rho = 0$ with $\Gamma$ the Christoffel symbols of that same chart.

After:

> We use coordinates $(t, r, \theta, \phi)$ in the global chart, with $t$ carrying dimensions of time and $r$ the areal radius, and $(t, x, y, z)$ in the Poincaré patch, with $z$ the direction off the boundary.
> We keep factors of $c$ and $G$ explicit, and the one parameter is the anti-de Sitter radius $L$, a length, which fixes the cosmological constant through $\Lambda = -3/L^2$; no mass enters, so $G$ appears in no component.
> We take components in the chart $x^0 = ct$, where the metric is dimensionless and no component carries a factor of $c$, and a dot in the geodesic equations is a derivative of that chart's coordinates with respect to an affine parameter, so $\dot{t}$ means $d(ct)/d\lambda$, and every term of every equation carries the dimensions of its left hand side.
> The geodesic equations read $\ddot{x}^\mu + \Gamma^\mu{}_{\nu\rho}\dot{x}^\nu\dot{x}^\rho = 0$, with $\Gamma$ the Christoffel symbols of the same chart.

### conventions: `anti_de_sitter/convention[1]`

Source: `MFS/assets/data/metrics/anti_de_sitter.json`, `convention`, paragraph 2.

Before:

> The signature is $(-,+,+,+)$, the Riemann tensor is $R^\mu{}_{\nu\rho\sigma} = \partial_\rho\Gamma^\mu{}_{\nu\sigma} - \partial_\sigma\Gamma^\mu{}_{\nu\rho} + \Gamma^\mu{}_{\rho\lambda}\Gamma^\lambda{}_{\nu\sigma} - \Gamma^\mu{}_{\sigma\lambda}\Gamma^\lambda{}_{\nu\rho}$, and the Ricci tensor is the standard contraction $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$.
> On that contraction the Ricci tensor carries the sign of $\Lambda$, and here that sign is negative.
> The spacetime is maximally symmetric, with ten Killing vectors and every sectional curvature equal to $-1/L^2$, so its Riemann tensor is $R_{\mu\nu\rho\sigma} = -\left(g_{\mu\rho}g_{\nu\sigma} - g_{\mu\sigma}g_{\nu\rho}\right)/L^2$ and every other curvature follows from it.
> Contracting it gives $R_{\mu\nu} = \Lambda g_{\mu\nu}$, $R = 4\Lambda = -12/L^2$ and the Kretschmann scalar $K = 8\Lambda^2/3 = 24/L^4$, all of them constant, so nothing curves more sharply anywhere than it does everywhere and there is no singularity to find.

After:

> We work in the signature $(-,+,+,+)$, with the Riemann tensor $R^\mu{}_{\nu\rho\sigma} = \partial_\rho\Gamma^\mu{}_{\nu\sigma} - \partial_\sigma\Gamma^\mu{}_{\nu\rho} + \Gamma^\mu{}_{\rho\lambda}\Gamma^\lambda{}_{\nu\sigma} - \Gamma^\mu{}_{\sigma\lambda}\Gamma^\lambda{}_{\nu\rho}$ and the Ricci tensor its standard contraction, $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$.
> With that contraction the Ricci tensor carries the sign of $\Lambda$, which here is negative.
> The spacetime is maximally symmetric, with ten Killing vectors and every sectional curvature equal to $-1/L^2$, so its Riemann tensor is $R_{\mu\nu\rho\sigma} = -\left(g_{\mu\rho}g_{\nu\sigma} - g_{\mu\sigma}g_{\nu\rho}\right)/L^2$, and every other curvature follows from it.
> Contracting it gives $R_{\mu\nu} = \Lambda g_{\mu\nu}$, $R = 4\Lambda = -12/L^2$, and the Kretschmann scalar $K = 8\Lambda^2/3 = 24/L^4$, all of them constant, so the curvature is the same at every point and there is no singularity.

### conventions: `anti_de_sitter/convention[2]`

Source: `MFS/assets/data/metrics/anti_de_sitter.json`, `convention`, paragraph 3.

Before:

> The stress energy vanishes and the field equations read $G_{\mu\nu} + \Lambda g_{\mu\nu} = 0$, so the Einstein tensor is not zero here: it is $G_{\mu\nu} = -\Lambda g_{\mu\nu} = 3g_{\mu\nu}/L^2$, which is minus the Ricci tensor component by component.
> The Weyl tensor vanishes identically, in both charts and in every slot, and not for the reason a vacuum's does: a maximally symmetric Riemann tensor is built from its own traces and nothing else, so removing them leaves nothing behind.
> Anti-de Sitter space is therefore conformally flat, which the Poincaré patch shows outright by writing the whole metric as a position dependent multiple of Minkowski space.

After:

> The stress energy vanishes, and the field equations read $G_{\mu\nu} + \Lambda g_{\mu\nu} = 0$, so the Einstein tensor does not: it is $G_{\mu\nu} = -\Lambda g_{\mu\nu} = 3g_{\mu\nu}/L^2$, which is minus the Ricci tensor component by component.
> The Weyl tensor vanishes identically, in both charts and in every slot, and the reason is the symmetry: a maximally symmetric Riemann tensor is built from its own traces and nothing else, so removing them leaves nothing behind.
> Anti-de Sitter space is therefore conformally flat, and in the Poincaré patch the whole metric is a position dependent multiple of Minkowski space.

### conventions: `anti_de_sitter/convention[3]`

Source: `MFS/assets/data/metrics/anti_de_sitter.json`, `convention`, paragraph 4.

Before:

> Two coordinate facts are the mirror of de Sitter's.
> The first is that $1 + r^2/L^2$ never vanishes, so $\partial_t$ is timelike everywhere, the static chart covers the whole spacetime and there is no horizon of any kind; de Sitter's $1 - r^2/L^2$ vanishes at $r = L$ and gives that spacetime its cosmological horizon, and this is exactly where the two solutions part.
> The second is that the edge is reachable.
> A radial light ray obeys $c\,dt = dr/(1 + r^2/L^2)$, whose integral out to $r \to \infty$ is finite, so light leaves the origin and arrives at spatial infinity at the coordinate time $ct = \pi L/2$, and can come back.
> The conformal boundary is timelike, no spacelike surface determines the future of the spacetime, and anti-de Sitter space is not globally hyperbolic: evolving anything inside it means stating a boundary condition at that edge, which is what makes the boundary a place a field theory can live.

After:

> Two coordinate facts mirror de Sitter's.
> The first is that $1 + r^2/L^2$ never vanishes, so $\partial_t$ is timelike everywhere, the static chart covers the whole spacetime, and there is no horizon of any kind; de Sitter's $1 - r^2/L^2$ vanishes at $r = L$ and gives that spacetime its cosmological horizon.
> The second is that the edge lies within reach.
> A radial light ray obeys $c\,dt = dr/(1 + r^2/L^2)$, whose integral out to $r \to \infty$ is finite, so light leaves the origin, arrives at spatial infinity at the coordinate time $ct = \pi L/2$, and can come back.
> The conformal boundary is timelike, no spacelike surface determines the future of the spacetime, and anti-de Sitter space is not globally hyperbolic: evolving anything inside it means stating a boundary condition at that edge, which is why a field theory can live on the boundary.

### conventions: `anti_de_sitter/convention[4]`

Source: `MFS/assets/data/metrics/anti_de_sitter.json`, `convention`, paragraph 5.

Before:

> Massive bodies never get there.
> A radial timelike geodesic turns around at a finite radius and refocuses at the origin after the coordinate time $ct = \pi L$, whatever energy it was thrown with.
> The chart takes $t$ over the whole real line, which is the universal cover: the quadric this metric is induced on closes $t$ into a circle of period $2\pi L/c$ and so carries closed timelike curves through every point, and unwrapping that circle is what removes them.

After:

> Massive bodies never reach the boundary.
> A radial timelike geodesic turns around at a finite radius and refocuses at the origin after the coordinate time $ct = \pi L$, whatever energy it was thrown with.
> The chart takes $t$ over the whole real line, which makes it the universal cover: the quadric on which this metric is induced closes $t$ into a circle of period $2\pi L/c$ and so carries closed timelike curves through every point, and unwrapping that circle removes them.

### conventions: `anti_de_sitter/convention[5]`

Source: `MFS/assets/data/metrics/anti_de_sitter.json`, `convention`, paragraph 6.

Before:

> The Poincaré patch is the chart the holography literature is written in, and it covers a wedge of the same spacetime rather than all of it.
> The conformal factor $L^2/z^2$ blows up at $z \to 0$, which is the conformal boundary, and the metric it induces there is Minkowski space in $(t, x, y)$ up to that factor, which is the flat spacetime the dual conformal field theory lives on.
> The other end, $z \to \infty$, is the Poincaré horizon, and it is a horizon of the chart and not of the spacetime: every curvature scalar takes the same constant value there as everywhere else, and the global chart runs straight through it.

After:

> The Poincaré patch, the chart of work on holography, covers only a wedge of the same spacetime.
> The conformal factor $L^2/z^2$ blows up at $z \to 0$, the conformal boundary, and up to that factor the metric induced there is Minkowski space in $(t, x, y)$, the flat spacetime on which the dual conformal field theory lives.
> The other end, $z \to \infty$, is the Poincaré horizon, a horizon of the chart alone: every curvature scalar takes the same constant value there as everywhere else, and the global chart runs straight through it.

### conformal diagram caption: `conformal/anti_de_sitter/poincare.caption[1]`

Source: `_tools/derivations/conformal.py` line 2310, `CAPTIONS ("anti_de_sitter", "poincare")`.

Before:

> They cover a wedge of it.
> Their $z \to 0$ is the conformal boundary, and $z \to \infty$ is the Poincaré horizon, the pair of null lines from $(\sigma, ct/L) = (-\pi/2, 0)$.
> The curvature there is the same as everywhere else, and the global coordinates run smoothly across it.

After:

> The Poincaré coordinates cover a wedge of the strip.
> Their $z \to 0$ is the conformal boundary, and $z \to \infty$ is the Poincaré horizon, the pair of null lines from $(\sigma, ct/L) = (-\pi/2, 0)$.
> The curvature there is the same as everywhere else, and the global coordinates run smoothly across it.

### embedding diagram caption: `embedding/anti_de_sitter/hyperboloid.caption[1]`

Source: `_tools/derivations/embedding.py` line 3582, `CAPTIONS ("anti_de_sitter", "hyperboloid")`.

Before:

> Where the sheet is steep a step along it is shorter than it looks: the circles at $L$, $2L$, $3L$ and $4L$ stand $0.88$, $0.56$, $0.37$ and $0.28\,L$ apart.
> The sheet nears the light cone of the space it is drawn in, dashed, without ever reaching it, and the conformal boundary of anti-de Sitter space lies along that cone at infinity.
> Wilhelm Killing in 1880 and Henri Poincaré in 1881 each described the hyperbolic plane as this sheet, and David Hilbert proved in 1901 that no surface in flat space carries the whole of it.

After:

> Where the sheet is steep a step along it is shorter than it looks: the circles at $L$, $2L$, $3L$, and $4L$ stand $0.88$, $0.56$, $0.37$, and $0.28\,L$ apart.
> The sheet nears the light cone of the space it is drawn in, dashed, without ever reaching it, and the conformal boundary of anti-de Sitter space lies along that cone at infinity.
> Wilhelm Killing in 1880 and Henri Poincaré in 1881 each described the hyperbolic plane as this sheet, and David Hilbert proved in 1901 that no surface in flat space carries the whole of it.

### embedding diagram, not drawn: `embedding/anti_de_sitter/hyperboloid.stops[0]`

Source: `_tools/derivations/embedding.py` line 3183, `in anti_de_sitter() stops=`.

Before:

> At every $r > 0$ the circles grow faster than the distance out to them, $g_{rr} < (\partial_r\sqrt{g_{\phi\phi}})^2$, so no surface of revolution in flat space carries the slice, and it is drawn in Minkowski space instead.

After:

> At every $r > 0$ the circles grow faster than the distance out to them, $g_{rr} < (\partial_r\sqrt{g_{\phi\phi}})^2$, and no surface of revolution in flat space carries the slice; Minkowski space carries it.

## Bertotti-Robinson Electrovacuum

### history: `bertotti_robinson/history[0]`

Source: `MFS/assets/data/metrics/bertotti_robinson.json`, `history`, paragraph 1.

Before:

> In 1959 two relativists, working independently, wrote down the same unusual solution.
> Bruno Bertotti in Princeton and Ivor Robinson in the same year each solved the coupled Einstein-Maxwell equations for the case of a uniform electromagnetic field: a field that is the same, in a precise covariant sense, at every event [bertotti1959, robinson1959].
> Most charged solutions of general relativity have a field that falls off with distance from a source; this one does not, and the price of that uniformity is paid in the geometry.

After:

> Two relativists, working independently, wrote down the same unusual solution in 1959.
> Bruno Bertotti in Princeton and Ivor Robinson in the same year each solved the coupled Einstein-Maxwell equations for the case of a uniform electromagnetic field: a field that is the same, in a precise covariant sense, at every event [bertotti1959, robinson1959].
> Most charged solutions of general relativity have a field that falls off with distance from a source; this one does not, and its uniform energy holds every sphere of the spacetime at one radius.

### history: `bertotti_robinson/history[3]`

Source: `MFS/assets/data/metrics/bertotti_robinson.json`, `history`, paragraph 4.

Before:

> For years the Bertotti-Robinson universe was a tidy oddity.
> Its importance became clear once black holes were understood: it is the limit near the horizon of an extremal Reissner-Nordström black hole, the geometry one finds deep in the throat of a charged hole pushed to the edge of having a horizon at all [reissner1916, nordstrom1918].
> Gary Gibbons and Paul Townsend described the extreme Reissner-Nordström black hole in 1993, as a solution of $N = 2$ supergravity, as a soliton that runs from flat spacetime at infinity down an infinite throat to AdS$_2 \times$ S$^2$, both of them vacua that keep all of the theory's supersymmetry [gibbons1993].
> Sergio Ferrara, Renata Kallosh, and Andrew Strominger found in 1995 that at the horizon of the extremal black holes of $N = 2$ supergravity the scalar fields flow to values fixed by the charges alone, and that the geometry there is the Bertotti-Robinson universe [ferrara1995].

After:

> For years the Bertotti-Robinson universe was a curiosity.
> Its importance became clear once black holes were understood: it is the limit near the horizon of an extremal Reissner-Nordström black hole, the geometry one finds deep in the throat of a charged hole whose charge has reached its mass, the largest charge at which it keeps a horizon [reissner1916, nordstrom1918].
> Gary Gibbons and Paul Townsend described the extreme Reissner-Nordström black hole in 1993, as a solution of $N = 2$ supergravity, as a soliton that runs from flat spacetime at infinity down an infinite throat to AdS$_2 \times$ S$^2$, both of them vacua that keep all of the theory's supersymmetry [gibbons1993].
> Sergio Ferrara, Renata Kallosh, and Andrew Strominger found in 1995 that at the horizon of the extremal black holes of $N = 2$ supergravity the scalar fields flow to values fixed by the charges alone, and that the geometry there is the Bertotti-Robinson universe [ferrara1995].

### conventions: `bertotti_robinson/convention[0]`

Source: `MFS/assets/data/metrics/bertotti_robinson.json`, `convention`, paragraph 1.

Before:

> Coordinates are $(t, r, \theta, \phi)$ in the static chart and $(t, x, \theta, \phi)$ in the Poincaré chart, with $t$ a time and $r$ and $x$ lengths, and factors of $c$ and $G$ are kept explicit.
> Neither $r$ nor $x$ is an areal radius: every sphere of the spacetime has the same areal radius $b$, and $r$ and $x$ run along the AdS$_2$ factor, related by $rx = b^2$.
> Neither set of coordinates is global.
> Under $x = b^2/r$ the static coordinates become the Poincaré coordinates, so the two cover one and the same Poincaré patch of AdS$_2$, and the static $r = 0$, which is the Poincaré $x \to \infty$, is a degenerate Killing horizon at the edge of that patch rather than the end of the spacetime, which continues through it into the rest of global AdS$_2$.

After:

> We use coordinates $(t, r, \theta, \phi)$ in the static chart and $(t, x, \theta, \phi)$ in the Poincaré chart, with $t$ a time and $r$ and $x$ lengths, and keep factors of $c$ and $G$ explicit.
> Neither $r$ nor $x$ is an areal radius: every sphere of the spacetime has the same areal radius $b$, and $r$ and $x$ run along the AdS$_2$ factor, related by $rx = b^2$.
> Neither set of coordinates is global.
> Under $x = b^2/r$ the static coordinates become the Poincaré coordinates, so the two cover one and the same Poincaré patch of AdS$_2$; the static $r = 0$, the Poincaré $x \to \infty$, is a degenerate Killing horizon at the edge of that patch, and the spacetime continues through it into the rest of global AdS$_2$.

### conventions: `bertotti_robinson/convention[1]`

Source: `MFS/assets/data/metrics/bertotti_robinson.json`, `convention`, paragraph 2.

Before:

> The single parameter $b$ is the common radius of the two factors and is fixed by the uniform electromagnetic field.
> Components are taken in the chart $x^0 = ct$, so the metric components carry no dimensions, and the dots in the geodesic equations are derivatives of that same chart with respect to an affine parameter, so $\dot{t}$ means $d(ct)/d\lambda$.
> Nothing in the metric depends on $t$, so no component carries a factor of $c$ at all.
> The Ricci tensor is the standard contraction $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$.
> The Ricci scalar vanishes and the Weyl tensor vanishes, the first because the two factors have equal and opposite curvature and the second because the spacetime is conformally flat, but the Ricci tensor does not: this is an electrovacuum, and the Einstein tensor is the Maxwell stress of a covariantly constant field, $G^t{}_t = G^r{}_r = -\dfrac{1}{b^2}$ and $G^\theta{}_\theta = G^\phi{}_\phi = +\dfrac{1}{b^2}$, which is a positive energy density $u = \dfrac{c^4}{8\pi G b^2}$ with a radial tension of the same size and an equal transverse pressure, the trace free stress of a radial electric field of strength $E = \dfrac{c^2}{b\sqrt{4\pi G\epsilon_0}}$.

After:

> The single parameter $b$ is the common radius of the two factors, and the uniform electromagnetic field fixes it.
> We take components in the chart $x^0 = ct$, so the metric components carry no dimensions, and a dot in the geodesic equations is a derivative of that chart's coordinates with respect to an affine parameter, so $\dot{t}$ means $d(ct)/d\lambda$.
> Nothing in the metric depends on $t$, so no component carries a factor of $c$.
> The Ricci tensor is the standard contraction, $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$.
> The Ricci scalar vanishes, since the two factors have equal and opposite curvature, and so does the Weyl tensor, since the spacetime is conformally flat, while the Ricci tensor does not: this is an electrovacuum, and the Einstein tensor is the Maxwell stress of a covariantly constant field, $G^t{}_t = G^r{}_r = -\dfrac{1}{b^2}$ and $G^\theta{}_\theta = G^\phi{}_\phi = +\dfrac{1}{b^2}$, a positive energy density $u = \dfrac{c^4}{8\pi G b^2}$ with a radial tension of the same size and an equal transverse pressure, the trace free stress of a radial electric field of strength $E = \dfrac{c^2}{b\sqrt{4\pi G\epsilon_0}}$.

## Bianchi Anisotropic Cosmologies

### history: `bianchi/history[1]`

Source: `MFS/assets/data/metrics/bianchi.json`, `history`, paragraph 2.

Before:

> Its relevance to physics rests on a single observation: a universe is spatially homogeneous if some symmetry group slides each spatial slice rigidly onto itself, and those groups are precisely Bianchi's nine.
> The Bianchi types therefore catalogue every way a universe can be the same at every point yet not the same in every direction: anisotropic generalizations of the homogeneous, isotropic Friedmann-Lemaître-Robertson-Walker model.
> Abraham Taub wrote down the vacuum members in 1951 [taub1951], and the subject was put on its modern footing by George Ellis and Malcolm MacCallum in 1969, whose division of the types into "class A" and "class B", according to whether a certain trace of the structure constants vanishes, remains the organizing scheme of the field [ellismaccallum1969].

After:

> Bianchi's classification matters to physics because of a single observation: a universe is spatially homogeneous if some symmetry group slides each spatial slice rigidly onto itself, and those groups are precisely Bianchi's nine.
> The Bianchi types therefore catalogue every way a universe can be the same at every point yet not the same in every direction: anisotropic generalizations of the homogeneous, isotropic Friedmann-Lemaître-Robertson-Walker model.
> Abraham Taub wrote down the vacuum members in 1951 [taub1951], and George Ellis and Malcolm MacCallum put the subject on its modern footing in 1969, dividing the types into "class A" and "class B" according to whether a certain trace of the structure constants vanishes, the organizing scheme of the field ever since [ellismaccallum1969].

### history: `bianchi/history[2]`

Source: `MFS/assets/data/metrics/bianchi.json`, `history`, paragraph 3.

Before:

> The nine types are told apart by their structure constants, which for each symmetry group reduce to a symmetric matrix with eigenvalues $(n_1, n_2, n_3)$, each free to be positive, negative, or zero, together with a vector whose single independent component is $a$.
> Two rules then fix everything: a type is class A when $a = 0$ and class B when $a \neq 0$, and the type itself is read off from the signs of $(n_1, n_2, n_3)$.
> The class and the type are therefore not labels but consequences of the two invariants, and the nine types, with the group and the isotropic FLRW limit each admits, are these:

After:

> The nine types are told apart by their structure constants, which for each symmetry group reduce to a symmetric matrix with eigenvalues $(n_1, n_2, n_3)$, each free to be positive, negative, or zero, together with a vector whose single independent component is $a$.
> Two rules then fix the class and the type: a type is class A when $a = 0$ and class B when $a \neq 0$, and the type itself follows from the signs of $(n_1, n_2, n_3)$.
> The class and the type therefore follow from the two invariants, and the nine types, with the group and the isotropic FLRW limit each admits, are these:

### history: `bianchi/history[4]`

Source: `MFS/assets/data/metrics/bianchi.json`, `history`, paragraph 5.

Before:

> Let us define these groups precisely.
> The Abelian group is flat space under addition, $\mathbb{R}^3 = \{\,(x,y,z) : x,y,z \in \mathbb{R}\,\}$ with group law $(\mathbf{u}, \mathbf{v}) \mapsto \mathbf{u} + \mathbf{v}$; every element commutes, which is what gives type I its flat slices.
> The Heisenberg group is the simplest case that does not commute, the unipotent upper triangular matrices $$H_3 = \left\{\, \begin{pmatrix} 1 & a & c \\ 0 & 1 & b \\ 0 & 0 & 1 \end{pmatrix} : a, b, c \in \mathbb{R} \,\right\},$$ whose single nontrivial commutator is the canonical relation of quantum mechanics.

After:

> Let us define these groups precisely.
> The Abelian group is flat space under addition, $\mathbb{R}^3 = \{\,(x,y,z) : x,y,z \in \mathbb{R}\,\}$ with group law $(\mathbf{u}, \mathbf{v}) \mapsto \mathbf{u} + \mathbf{v}$; every element commutes, which gives type I its flat slices.
> The Heisenberg group is the simplest case that does not commute, the unipotent upper triangular matrices $$H_3 = \left\{\, \begin{pmatrix} 1 & a & c \\ 0 & 1 & b \\ 0 & 0 & 1 \end{pmatrix} : a, b, c \in \mathbb{R} \,\right\},$$ whose single nontrivial commutator is the canonical relation of quantum mechanics.

### history: `bianchi/history[5]`

Source: `MFS/assets/data/metrics/bianchi.json`, `history`, paragraph 6.

Before:

> The two that matter most for cosmology are those of types IX and VIII.
> The first is the special unitary group $$SU(2) = \left\{\, U \in \mathrm{M}_2(\mathbb{C}) : U^\dagger U = \mathbb{1},\ \det U = 1 \,\right\},$$ the double cover of $SO(3)$ and the isometry group of the round $S^3$, which is why type IX has spherical slices and a closed FLRW limit.
> Its noncompact counterpart is the special linear group $$SL(2, \mathbb{R}) = \left\{\, M \in \mathrm{M}_2(\mathbb{R}) : \det M = 1 \,\right\},$$ whose universal cover gives the type VIII group.
> Because $SL(2,\mathbb{R})$ is not simply connected, since it retracts onto its $SO(2)$ rotation subgroup and so has fundamental group $\pi_1 = \mathbb{Z}$, that cover, written $\widetilde{SL(2,\mathbb{R})}$, is defined as the unique simply connected Lie group sharing the Lie algebra $\mathfrak{sl}(2,\mathbb{R})$: an infinite cover with one sheet for each integer, one that, unlike every other group of the nine, admits no faithful matrix representation in finitely many dimensions.
> It plays the role for hyperbolic geometry that $SU(2)$ plays for the sphere.

After:

> The groups of types IX and VIII matter most for cosmology.
> The first is the special unitary group $$SU(2) = \left\{\, U \in \mathrm{M}_2(\mathbb{C}) : U^\dagger U = \mathbb{1},\ \det U = 1 \,\right\},$$ the double cover of $SO(3)$ and the isometry group of the round $S^3$, which is why type IX has spherical slices and a closed FLRW limit.
> Its noncompact counterpart is the special linear group $$SL(2, \mathbb{R}) = \left\{\, M \in \mathrm{M}_2(\mathbb{R}) : \det M = 1 \,\right\},$$ whose universal cover gives the type VIII group.
> Because $SL(2,\mathbb{R})$ is not simply connected, since it retracts onto its $SO(2)$ rotation subgroup and so has fundamental group $\pi_1 = \mathbb{Z}$, that cover, written $\widetilde{SL(2,\mathbb{R})}$, is defined as the unique simply connected Lie group sharing the Lie algebra $\mathfrak{sl}(2,\mathbb{R})$: an infinite cover with one sheet for each integer, one that, unlike every other group of the nine, admits no faithful matrix representation in finitely many dimensions.
> It plays the role for hyperbolic geometry that $SU(2)$ plays for the sphere.

### history: `bianchi/history[6]`

Source: `MFS/assets/data/metrics/bianchi.json`, `history`, paragraph 7.

Before:

> The dynamics are where the Bianchi cosmologies turn strange.
> The class A models can be cast as a particle moving in a potential, and their behavior on the approach to the singularity ranges from the simple to the chaotic.
> Type I, the Kasner solution, expands along two axes while contracting along the third, a single fixed anisotropic motion.
> Type IX is another matter entirely: this is the Mixmaster universe, introduced by Charles Misner in 1969, whose approach to the singularity is an unending, chaotic sequence of Kasner epochs, the axes of expansion and contraction trading places in a sequence that never settles [misner1969].
> The same year and shortly after, Vladimir Belinskii, Isaak Khalatnikov, and Evgeny Lifshitz advanced the conjecture that gives the Bianchi models their deepest significance: that this Mixmaster oscillation is not peculiar to type IX but is the generic behavior of any spacetime near a singularity, with neighboring points decoupling so that each evolves like its own independent Bianchi IX universe [bkl1970].
> If the BKL conjecture holds, the structure of the cosmic singularity is, locally and everywhere, a Bianchi cosmology.

After:

> The Bianchi cosmologies turn strange in their dynamics.
> The class A models can be cast as a particle moving in a potential, and their behavior on the approach to the singularity ranges from the simple to the chaotic.
> Type I, the Kasner solution, expands along two axes while contracting along the third, a single fixed anisotropic motion.
> Type IX is another matter entirely: this is the Mixmaster universe, introduced by Charles Misner in 1969, whose approach to the singularity is an unending, chaotic sequence of Kasner epochs, the axes of expansion and contraction trading places in a sequence that never settles [misner1969].
> The same year and shortly after, Vladimir Belinskii, Isaak Khalatnikov, and Evgeny Lifshitz advanced the conjecture that gives the Bianchi models their deepest significance: that this Mixmaster oscillation is the generic behavior of any spacetime near a singularity, type IX being one instance of it, with neighboring points decoupling so that each evolves like its own independent Bianchi IX universe [bkl1970].
> If the BKL conjecture holds, the structure of the cosmic singularity is, locally and everywhere, a Bianchi cosmology.

### history: `bianchi/history[7]`

Source: `MFS/assets/data/metrics/bianchi.json`, `history`, paragraph 8.

Before:

> The diagonal type I model, $ds^2 = -c^2dt^2 + a_1(t)^2dx^2 + a_2(t)^2dy^2 + a_3(t)^2dz^2$, is the member the rest of the family is measured against.
> Its slices are ordinary flat space, so the three scale factors are the whole of the geometry and the whole of the anisotropy, and the shape of the curvature is easy to say: expansion along an axis is carried by $a_i''/a_i$ and the shear between two axes by the product of their expansion rates.
> The Ricci and Einstein tensors are given for arbitrary scale factors rather than for a vacuum, because the useful case is type I with matter in it, where the time component of the Einstein tensor returns the density and the three spatial components return three separate pressures, one per axis.
> Take the matter away and the same four equations force the scale factors to be powers of the time whose exponents sum to one along with their squares, which is the Kasner solution.
> Type I with a fluid and Kasner are one family at two settings of the matter, and type IX is where the family stops being solvable in closed form and becomes the Mixmaster universe.

After:

> The diagonal type I model, $ds^2 = -c^2dt^2 + a_1(t)^2dx^2 + a_2(t)^2dy^2 + a_3(t)^2dz^2$, is the member against which the rest of the family is measured.
> Its slices are ordinary flat space, so the three scale factors are the whole of the geometry and the whole of the anisotropy, and the shape of the curvature is easy to say: expansion along an axis is carried by $a_i''/a_i$ and the shear between two axes by the product of their expansion rates.
> We take the Ricci and Einstein tensors for arbitrary scale factors, because the useful case is type I with matter in it, where the time component of the Einstein tensor returns the density and the three spatial components return three separate pressures, one per axis.
> Without the matter, the same four equations force the scale factors to be powers of the time whose exponents sum to one along with their squares, the Kasner solution.
> Type I with a fluid and Kasner's vacuum belong to one family, and in type IX the family can no longer be solved in closed form and becomes the Mixmaster universe.

### history: `bianchi/history[8]`

Source: `MFS/assets/data/metrics/bianchi.json`, `history`, paragraph 9.

Before:

> A classification of symmetric spaces from the nineteenth century became, in turn, the framework for anisotropic cosmology, the testing ground for the cosmological principle, and the key to the singularity at the beginning of time.
> The very nearly isotropic universe we observe sits inside the Bianchi family as a special, finely tuned member, and that is exactly why the family is indispensable.
> It supplies the space of alternatives against which isotropy can be measured, the answer to "how anisotropic could the early universe have been, and how do we know it wasn't?" The Bianchi universes are less a single spacetime than a periodic table of them, the complete enumeration of homogeneous worlds, within which our own isotropic cosmos is one carefully balanced element.

After:

> A classification of symmetric spaces from the nineteenth century became the framework for anisotropic cosmology, for tests of the cosmological principle, and for the study of the singularity at the beginning of time.
> It is the complete list of homogeneous cosmologies, and the very nearly isotropic universe we observe sits inside it as a special, finely tuned member, which is why cosmologists need the whole family.
> The family supplies the alternatives against which isotropy can be measured, and with them the question "how anisotropic could the early universe have been, and how do we know it wasn't?"

### conventions: `bianchi/convention[0]`

Source: `MFS/assets/data/metrics/bianchi.json`, `convention`, paragraph 1.

Before:

> Coordinates are $(t, x, y, z)$ with $t$ a cosmic time and $x$, $y$ and $z$ comoving Cartesian lengths, and factors of $c$ and $G$ are kept explicit.
> The three scale factors $a_1(t)$, $a_2(t)$ and $a_3(t)$ are positive and dimensionless, so every length sits in the coordinates; the chart covers any interval on which all three are smooth and positive, and it is written as $t > 0$ because the models of interest run into a singularity at a finite time, which is placed at $t = 0$.
> Components are taken in the chart $x^0 = ct$, where a chart component of the metric carries no dimensions.
> A prime is the derivative with respect to that chart time, $a_i' = \dfrac{da_i}{d(ct)} = \dfrac{1}{c}\dfrac{da_i}{dt}$ and $a_i'' = \dfrac{1}{c^2}\dfrac{d^2a_i}{dt^2}$, so every rate carries an inverse length and the factors of $c$ the chart demands sit in the definition of the prime rather than scattered through the components.
> The dots in the geodesic equations are derivatives of those same chart coordinates with respect to an affine parameter, so $\dot{t}$ means $d(ct)/d\lambda$.

After:

> We use coordinates $(t, x, y, z)$, with $t$ a cosmic time and $x$, $y$, and $z$ comoving Cartesian lengths, and keep factors of $c$ and $G$ explicit.
> The three scale factors $a_1(t)$, $a_2(t)$, and $a_3(t)$ are positive and dimensionless, so every length sits in the coordinates; the chart covers any interval on which all three are smooth and positive, and we write it as $t > 0$ because the models of interest run into a singularity at a finite time, which we place at $t = 0$.
> We take components in the chart $x^0 = ct$, where a chart component of the metric carries no dimensions.
> A prime is the derivative with respect to that chart time, $a_i' = \dfrac{da_i}{d(ct)} = \dfrac{1}{c}\dfrac{da_i}{dt}$ and $a_i'' = \dfrac{1}{c^2}\dfrac{d^2a_i}{dt^2}$, so every rate carries an inverse length, and the chart's factors of $c$ sit in the definition of the prime instead of being scattered through the components.
> A dot in the geodesic equations is a derivative of the same chart coordinates with respect to an affine parameter, so $\dot{t}$ means $d(ct)/d\lambda$.

### conventions: `bianchi/convention[1]`

Source: `MFS/assets/data/metrics/bianchi.json`, `convention`, paragraph 2.

Before:

> The Ricci tensor is the standard contraction $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$.
> No component assumes a field equation, because type I is at its most useful with matter in it: the scale factors are left arbitrary, and a comoving perfect fluid is read straight off the Einstein tensor, whose time component gives $G_{tt} = \dfrac{8\pi G}{c^2}\rho$ and whose three spatial components give one principal pressure per axis, $G^x{}_x = \dfrac{8\pi G}{c^4}p_x$, which is what an anisotropic universe carries instead of a single $p$.
> Empty it out, $\rho = p_x = p_y = p_z = 0$, and those same four equations force $a_i = t^{p_i}$ with $p_1 + p_2 + p_3 = p_1^2 + p_2^2 + p_3^2 = 1$: Kasner's solution is the vacuum specialisation of type I.

After:

> The Ricci tensor is the standard contraction, $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$.
> No component assumes a field equation, because type I is at its most useful with matter in it: we leave the scale factors arbitrary and take a comoving perfect fluid straight from the Einstein tensor, whose time component gives $G_{tt} = \dfrac{8\pi G}{c^2}\rho$ and whose three spatial components give one principal pressure per axis, $G^x{}_x = \dfrac{8\pi G}{c^4}p_x$, which an anisotropic universe carries in place of a single $p$.
> Emptied of matter, $\rho = p_x = p_y = p_z = 0$, those same four equations force $a_i = t^{p_i}$ with $p_1 + p_2 + p_3 = p_1^2 + p_2^2 + p_3^2 = 1$, and Kasner's solution is the vacuum specialisation of type I.

### spacetime diagram caption: `diagrams/bianchi/systems/type_i_cartesian/tx.caption[1]`

Source: `_tools/derivations/null_rays.py` line 852, `CAPTIONS ("bianchi", "type_i_cartesian", "tx")`.

Before:

> So along $x$ the cones close toward $t = 0$, as in Kasner's contracting direction, and open again later as $a_1$ turns round.

After:

> Along $x$ the cones therefore close toward $t = 0$, as in Kasner's contracting direction, and open again later as $a_1$ turns round.

### embedding diagram caption: `embedding/bianchi/ring.caption[0]`

Source: `_tools/derivations/embedding.py` line 3632, `CAPTIONS ("bianchi", "ring")`.

Before:

> This is the plane $y = 0$ of a Bianchi type I universe of dust at four moments of cosmic time, each drawn as a surface in flat space so that every distance along it is the metric distance.
> At every moment the plane is flat, $a_1^2dx^2 + a_3^2dz^2$ being Euclid's plane with its axes scaled, so the drawing is a flat disc, and what the geometry does shows in a ring of the dust itself, whose grains stay at rest in the chart.

After:

> This is the plane $y = 0$ of a Bianchi type I universe of dust at four moments of cosmic time, each drawn as a surface in flat space so that every distance along it is the metric distance.
> At every moment the plane is flat, $a_1^2dx^2 + a_3^2dz^2$ being Euclid's plane with its axes scaled, so the drawing is a flat disc, and the uneven expansion shows in a ring of the dust itself, whose grains stay at rest in the chart.

### embedding diagram note: `embedding/bianchi/ring.input`

Source: `_tools/derivations/embedding.py` line 3025, `in bianchi() input=`.

Before:

> Dust: the three scale factors solved from this spacetime's own $G^x{}_x = G^y{}_y = G^z{}_z = 0$, starting from $a_i = 1$ with rates $(-0.5, 1.5, 2.0)\,\bar H$, $\bar H$ their mean, as the spacetime diagram declares.

After:

> Dust: the three scale factors solved from this spacetime's own $G^x{}_x = G^y{}_y = G^z{}_z = 0$, starting from $a_i = 1$ with rates $(-0.5, 1.5, 2.0)\,\bar H$, $\bar H$ their mean, as in the spacetime diagram.

## Vilenkin-Gott Cosmic String

### history: `cosmic_string/history[0]`

Source: `MFS/assets/data/metrics/cosmic_string.json`, `history`, paragraph 1.

Before:

> Tom Kibble spent the middle of the 1970s asking what a cooling universe does with a symmetry it can no longer afford, and answered in 1976 that the vacuum will not agree with itself.
> Along certain lines the old symmetric phase stays trapped, held there by the topology of the vacuum manifold and not by anything local [kibble1976].
> The defect left behind is a thread thinner than a proton, drawn across cosmological distances, its tension equal to its mass per unit length.

After:

> Tom Kibble spent the middle of the 1970s asking what a cooling universe does with a symmetry it breaks, and answered in 1976 that regions too far apart to communicate break it differently and cannot always be matched where they meet.
> Along certain lines the old symmetric phase stays trapped, held there by the topology of the vacuum manifold alone [kibble1976].
> The defect left behind is a thread thinner than a proton, drawn across cosmological distances, its tension equal to its mass per unit length.

### history: `cosmic_string/history[2]`

Source: `MFS/assets/data/metrics/cosmic_string.json`, `history`, paragraph 3.

Before:

> J. Richard Gott gave the exact solution in 1985, that cone outside and a spherical cap within, and read the optics off it [gott1985].
> The missing wedge lenses by itself, so one quasar behind the string shows up twice at equal brightness, and the microwave sky steps in temperature across its track.
> Robert Brandenberger and Neil Turok computed in 1986 what a network of strings would stamp on the microwave background [brandenberger1986], and in a 1994 review Brandenberger set out the whole defect theory of structure formation [brandenberger1994].
> The acoustic peaks later measured so cleanly settled it against: they are the sound of a plasma driven in step, and strings do not drive it that way.

After:

> J. Richard Gott gave the exact solution in 1985, that cone outside and a spherical cap within, and worked out its optics [gott1985].
> The missing wedge lenses by itself, so one quasar behind the string shows up twice at equal brightness, and the microwave sky steps in temperature across its track.
> Robert Brandenberger and Neil Turok computed in 1986 what a network of strings would stamp on the microwave background [brandenberger1986], and in a 1994 review Brandenberger set out the whole defect theory of structure formation [brandenberger1994].
> The acoustic peaks measured later decided against it: they are the sound of a plasma set oscillating in step, and strings do not set it oscillating that way.

### history: `cosmic_string/history[3]`

Source: `MFS/assets/data/metrics/cosmic_string.json`, `history`, paragraph 4.

Before:

> Then Gott came back to his own solution and used it to build a time machine.
> In 1991 he gave the exact geometry of two straight strings passing one another without touching, and showed that if each moves fast enough, a traveler who loops around both while they pass returns to his own past [gott1991].
> The condition he gives is $\gamma_s > (\sin 4\pi\mu)^{-1}$ for each string in the laboratory frame, with $\mu$ the mass per unit length in Planck units and $8\pi\mu$ the deficit angle.
> Neither string misbehaves.

After:

> Then Gott came back to his own solution and used it to build a time machine.
> In 1991 he gave the exact geometry of two straight strings passing one another without touching, and showed that if each moves fast enough, a traveler who loops around both while they pass returns to their own past [gott1991].
> The condition he gives is $\gamma_s > (\sin 4\pi\mu)^{-1}$ for each string in the laboratory frame, with $\mu$ the mass per unit length in Planck units and $8\pi\mu$ the deficit angle.
> Each string on its own is an ordinary source of positive energy.

### history: `cosmic_string/history[5]`

Source: `MFS/assets/data/metrics/cosmic_string.json`, `history`, paragraph 6.

Before:

> The objections came fast.
> Stanley Deser, Roman Jackiw and Gerard 't Hooft, who had made gravity in three dimensions a flat space problem in 1984 [deser1984], added up the energy and momentum of Gott's pair and found it tachyonic: together the strings move faster than light though neither does alone, so the pair cannot be realized by physical, timelike sources [deser1992].
> Then came the attempts to build one anyway.
> Four pages earlier in that same issue of Physical Review Letters, Sean Carroll, Edward Farhi and Alan Guth had built the spacetime of one gravitating particle decaying into two and found that an open universe never holds enough energy to reach Gott's threshold [carroll1992].
> A closed universe will let you assemble the pair, and 't Hooft shut that door too: such a universe has a finite lifetime and crunches before anybody finishes the circuit [thooft1992].
> Amos Ori had come at it earlier still: if such a spacetime has closed curves at all, every surface of constant time in it meets one, so the machine was never switched on [ori1991].

After:

> The objections came fast.
> Stanley Deser, Roman Jackiw, and Gerard 't Hooft, who had made gravity in three dimensions a flat space problem in 1984 [deser1984], added up the energy and momentum of Gott's pair and found it tachyonic: together the strings move faster than light though neither does alone, so the pair cannot be realized by physical, timelike sources [deser1992].
> Then came the attempts to build one anyway.
> Four pages earlier in that same issue of Physical Review Letters, Sean Carroll, Edward Farhi, and Alan Guth had built the spacetime of one gravitating particle decaying into two and found that an open universe never holds enough energy to reach Gott's threshold [carroll1992].
> A closed universe will let you assemble the pair, and 't Hooft showed that this fails too: such a universe has a finite lifetime and crunches before anybody finishes the circuit [thooft1992].
> Amos Ori had come at it earlier still: if such a spacetime has any closed curves, every surface of constant time in it meets one, so the machine was never switched on [ori1991].

### history: `cosmic_string/history[6]`

Source: `MFS/assets/data/metrics/cosmic_string.json`, `history`, paragraph 7.

Before:

> Hawking generalized the quarrel in 1992.
> Causality violation growing in a finite region out of an infinite initial surface must break the averaged weak energy condition on its Cauchy horizon, which disposes of closed curves made from finite lengths of string; and as a curve nears closure the stress tensor's expectation value runs away, wrecking the horizon before any traveler arrives [hawking1992].
> In his words, "It seems that there is a Chronology Protection Agency which prevents the appearance of closed timelike curves and so makes the universe safe for historians."

After:

> Stephen Hawking generalized the quarrel in 1992.
> Causality violation growing in a finite region out of an infinite initial surface must break the averaged weak energy condition on its Cauchy horizon, which disposes of closed curves made from finite lengths of string; and as a curve nears closure the stress tensor's expectation value runs away, wrecking the horizon before any traveler arrives [hawking1992].
> In his words, "It seems that there is a Chronology Protection Agency which prevents the appearance of closed timelike curves and so makes the universe safe for historians."

### history: `cosmic_string/history[7]`

Source: `MFS/assets/data/metrics/cosmic_string.json`, `history`, paragraph 8.

Before:

> These arguments rule out manufacturing such a spacetime, but chronology protection is still a conjecture, since the divergence it leans on arrives where the semiclassical approximation gives out [visser2003].
> Gott's spacetime is still there, an exact solution made of two threads of positive energy obeying every energy condition, with no wormhole, no exotic matter and no singularity in it.
> Ian Stewart told the whole story as fiction in his novel Flatterland of 2001, whose characters learn to make the string's cone by cutting a wedge from the plane and gluing its edges, meet Gott's time machine of two strings passing each other at nearly light speed, and are told that their universe holds too little energy to build one [stewart2001flatterland].
> That is why an argument about time travel so often ends up back at two cosmic strings.

After:

> These arguments rule out manufacturing such a spacetime, but chronology protection is still a conjecture, since the divergence it leans on arrives where the semiclassical approximation gives out [visser2003].
> Gott's spacetime is still there, an exact solution made of two threads of positive energy obeying every energy condition, with no wormhole, no exotic matter, and no singularity in it.
> Arguments about time travel often return to it, since it needs nothing exotic.
> Ian Stewart told the whole story as fiction in his novel Flatterland of 2001, whose characters learn to make the string's cone by cutting a wedge from the plane and gluing its edges, meet Gott's time machine of two strings passing each other at nearly light speed, and are told that their universe holds too little energy to build one [stewart2001flatterland].

### conventions: `cosmic_string/convention[0]`

Source: `MFS/assets/data/metrics/cosmic_string.json`, `convention`, paragraph 1.

Before:

> Coordinates of the conical exterior are $(t, r, \phi, z)$ with $t$ carrying dimensions of time, $r$ the proper distance from the string, $\phi$ running the full $2\pi$ and $z$ the proper length along the string, and factors of $c$ and $G$ are kept explicit.
> The string's mass per unit length is $\mu$, equal to its tension, the one dimensionless combination it forms is $4G\mu/c^2$, and the deficit angle is $\delta = 8\pi G\mu/c^2$.
> The exterior is written with $\phi$ periodic in $2\pi$ and the deficit carried by the metric, in $g_{\phi\phi} = (1 - 4G\mu/c^2)^2r^2$, rather than with $g_{\phi\phi} = r^2$ and $\phi$ closing short of $2\pi$; the two are the same geometry.

After:

> For the conical exterior we use coordinates $(t, r, \phi, z)$, with $t$ carrying dimensions of time, $r$ the proper distance from the string, $\phi$ running the full $2\pi$, and $z$ the proper length along the string, and we keep factors of $c$ and $G$ explicit.
> The string's mass per unit length is $\mu$, equal to its tension; the one dimensionless combination it forms is $4G\mu/c^2$, and the deficit angle is $\delta = 8\pi G\mu/c^2$.
> We write the exterior with $\phi$ periodic in $2\pi$ and the deficit carried by the metric, in $g_{\phi\phi} = (1 - 4G\mu/c^2)^2r^2$, where one could equally take $g_{\phi\phi} = r^2$ and let $\phi$ close short of $2\pi$; the two are the same geometry.

### conventions: `cosmic_string/convention[1]`

Source: `MFS/assets/data/metrics/cosmic_string.json`, `convention`, paragraph 2.

Before:

> Components are taken in the chart $x^0 = ct$, where a chart component of the metric carries no dimensions, and the dots in the geodesic equations are derivatives of that same chart with respect to an affine parameter, so $\dot{t}$ means $d(ct)/d\lambda$ and every term of every equation carries the dimensions of its left hand side.
> The geodesic equations are $\ddot{x}^\mu + \Gamma^\mu{}_{\nu\rho}\dot{x}^\nu\dot{x}^\rho = 0$ with $\Gamma$ the Christoffel symbols of that same chart.
> The signature is $(-,+,+,+)$, the Riemann tensor is $R^\mu{}_{\nu\rho\sigma} = \partial_\rho\Gamma^\mu{}_{\nu\sigma} - \partial_\sigma\Gamma^\mu{}_{\nu\rho} + \Gamma^\mu{}_{\rho\lambda}\Gamma^\lambda{}_{\nu\sigma} - \Gamma^\mu{}_{\sigma\lambda}\Gamma^\lambda{}_{\nu\rho}$, and the Ricci tensor is the standard contraction $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$.

After:

> We take components in the chart $x^0 = ct$, where a chart component of the metric carries no dimensions, and a dot in the geodesic equations is a derivative of that chart's coordinates with respect to an affine parameter, so $\dot{t}$ means $d(ct)/d\lambda$, and every term of every equation carries the dimensions of its left hand side.
> The geodesic equations read $\ddot{x}^\mu + \Gamma^\mu{}_{\nu\rho}\dot{x}^\nu\dot{x}^\rho = 0$, with $\Gamma$ the Christoffel symbols of the same chart.
> We work in the signature $(-,+,+,+)$, with the Riemann tensor $R^\mu{}_{\nu\rho\sigma} = \partial_\rho\Gamma^\mu{}_{\nu\sigma} - \partial_\sigma\Gamma^\mu{}_{\nu\rho} + \Gamma^\mu{}_{\rho\lambda}\Gamma^\lambda{}_{\nu\sigma} - \Gamma^\mu{}_{\sigma\lambda}\Gamma^\lambda{}_{\nu\rho}$ and the Ricci tensor its standard contraction, $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$.

### conventions: `cosmic_string/convention[2]`

Source: `MFS/assets/data/metrics/cosmic_string.json`, `convention`, paragraph 3.

Before:

> Every curvature tensor of the exterior vanishes, and that emptiness is the physics of this spacetime.
> The Riemann tensor vanishes identically, and with it the Ricci tensor, the Ricci scalar, the Einstein tensor, the Kretschmann scalar and the Weyl tensor, because the cone is flat at every point where it is defined: the substitution $\tilde\phi = (1 - 4G\mu/c^2)\phi$ turns the line element into that of Minkowski space in cylindrical coordinates, and only the range of the new angle remembers the string, since $\tilde\phi$ closes at $2\pi(1 - 4G\mu/c^2)$ instead of at $2\pi$.
> No measurement made inside one neighbourhood tells the exterior from empty space.
> The whole of the gravitation is the global fact that a circle of proper radius $r$ about the string has circumference $2\pi(1 - 4G\mu/c^2)r$, a wedge of angle $\delta$ short of a full one, and that missing wedge is what lenses: two rays passing the string on opposite sides arrive $\delta$ apart, so one quasar behind it is seen twice, at equal brightness and with neither image magnified.

After:

> Every curvature tensor of the exterior vanishes.
> The Riemann tensor vanishes identically, and with it the Ricci tensor, the Ricci scalar, the Einstein tensor, the Kretschmann scalar, and the Weyl tensor, because the cone is flat at every point where it is defined: the substitution $\tilde\phi = (1 - 4G\mu/c^2)\phi$ turns the line element into that of Minkowski space in cylindrical coordinates, and only the range of the new angle remembers the string, since $\tilde\phi$ closes at $2\pi(1 - 4G\mu/c^2)$ instead of at $2\pi$.
> No measurement made inside one neighbourhood tells the exterior from empty space.
> The string's gravitation lies wholly in the global fact that a circle of proper radius $r$ about it has circumference $2\pi(1 - 4G\mu/c^2)r$, a wedge of angle $\delta$ short of a full one, and the missing wedge lenses light: two rays passing the string on opposite sides arrive $\delta$ apart, so one quasar behind it is seen twice, at equal brightness and with neither image magnified.

### conventions: `cosmic_string/convention[3]`

Source: `MFS/assets/data/metrics/cosmic_string.json`, `convention`, paragraph 4.

Before:

> The second chart is Gott's 1985 interior, the source the cone is the outside of, and the one place in this spacetime where the curvature is not zero.
> It is a cap of a sphere of radius $\ell$ carried rigidly along the string, written with the polar angle $\chi$ from the axis, so the proper distance from the axis is $\ell\chi$ and the core reaches $\chi_0$.
> Its Einstein tensor is nonzero, $G^t{}_t = G^z{}_z = -1/\ell^2 = -8\pi G\rho/c^2$, which is a uniform energy density $\rho c^2$ with a tension of the same size along $z$ and no pressure across the string, exactly the stress a straight string carries.

After:

> The second chart is Richard Gott's interior of 1985, the source whose outside is the cone, and the one place in this spacetime where the curvature is not zero.
> It is a cap of a sphere of radius $\ell$ carried rigidly along the string, written with the polar angle $\chi$ from the axis, so the proper distance from the axis is $\ell\chi$, and the core reaches $\chi_0$.
> Its Einstein tensor is nonzero, $G^t{}_t = G^z{}_z = -1/\ell^2 = -8\pi G\rho/c^2$: a uniform energy density $\rho c^2$ with a tension of the same size along $z$ and no pressure across the string, which is the stress of a straight string.

### conventions: `cosmic_string/convention[4]`

Source: `MFS/assets/data/metrics/cosmic_string.json`, `convention`, paragraph 5.

Before:

> Its Kretschmann scalar is the constant $4/\ell^4$, so the cap is as smooth as the cone is flat, and its Weyl tensor is nonzero, because a flat plane multiplied by a sphere is not conformally flat.
> That Weyl tensor fills the $t$ and $z$ slots, where the Riemann tensor is empty, so it is neither a copy of Riemann nor confined to the two directions Riemann lives in.
> Matching the cap to the cone at $\chi_0$ gives $\cos\chi_0 = 1 - 4G\mu/c^2$, and integrating $\rho$ over the cap gives $\mu = 2\pi\rho\ell^2(1 - \cos\chi_0)$, which is the same equation read the other way.
> Gauss and Bonnet make the agreement inevitable: the deficit angle is the total Gaussian curvature of the transverse section, $\delta = \int K_{\mathrm{G}}\,dA = 2\pi(1 - \cos\chi_0) = 8\pi G\mu/c^2$, so the flat exterior and the curved interior are two readings of one number.

After:

> The core's Kretschmann scalar is the constant $4/\ell^4$, so the cap is smooth everywhere, and its Weyl tensor is nonzero, because a flat plane multiplied by a sphere is not conformally flat.
> That Weyl tensor fills the $t$ and $z$ slots, where the Riemann tensor is empty, so it differs from Riemann and reaches beyond the two directions in which Riemann lives.
> Matching the cap to the cone at $\chi_0$ gives $\cos\chi_0 = 1 - 4G\mu/c^2$, and integrating $\rho$ over the cap gives $\mu = 2\pi\rho\ell^2(1 - \cos\chi_0)$, the same equation solved the other way round.
> The theorem of Gauss and Bonnet guarantees the agreement: the deficit angle is the total Gaussian curvature of the transverse section, $\delta = \int K_{\mathrm{G}}\,dA = 2\pi(1 - \cos\chi_0) = 8\pi G\mu/c^2$, so the deficit of the flat exterior equals the integrated curvature of the interior.

## de Sitter

### history: `de_sitter/history[0]`

Source: `MFS/assets/data/metrics/de_sitter.json`, `history`, paragraph 1.

Before:

> In 1917, having just published the general theory two years earlier, Albert Einstein added a term to his field equations, the cosmological constant $\Lambda$, in an attempt to hold the universe static against its own gravity [einstein1917].
> He had a particular cosmos in mind: ordinary matter spread uniformly through a closed, spherical space, with $\Lambda$'s repulsion finely tuned to balance the pull of the dust on itself exactly, so that nothing moved for all time.
> Within months the Dutch astronomer Willem de Sitter found a second solution to the very same amended equations, and it was everything Einstein's was not [desitter1917].
> He kept the same $\Lambda$ and threw out the matter, and the empty world that remained refused to sit still.

After:

> Two years after publishing the general theory, Albert Einstein added a term to his field equations in 1917, the cosmological constant $\Lambda$, in an attempt to hold the universe static against its own gravity [einstein1917].
> He had a particular cosmos in mind: ordinary matter spread uniformly through a closed, spherical space, with $\Lambda$'s repulsion finely tuned to balance the pull of the dust on itself exactly, so that nothing moved for all time.
> Within months the Dutch astronomer Willem de Sitter found a second solution to the very same amended equations [desitter1917].
> He kept the same $\Lambda$ and threw out the matter, and the empty world that remained was not static.

### history: `de_sitter/history[1]`

Source: `MFS/assets/data/metrics/de_sitter.json`, `history`, paragraph 2.

Before:

> The distinction was never about which side of the equation $\Lambda$ was written on; that choice, geometry or vacuum energy, is a convention and changes no physics.
> What separated the two solutions was their matter content: Einstein's a balance of dust against $\Lambda$ on a knife edge, de Sitter's the empty limit where only $\Lambda$ remains.
> Eddington would later sum up the pair with a line that has outlived most of the physics around it: Einstein's was matter without motion, de Sitter's motion without matter [eddington1933].
> Peter Phillips built his story "Manna", in Astounding Science Fiction of February 1949, on that emptiness: the ghosts of two medieval monks leave a physicist a de Sitter version of the Riemann-Christoffel tensor, and he has to work out how anything can exist in a universe with neither matter nor radiation [phillips1949manna].

After:

> Which side of the equation $\Lambda$ is written on, geometry or vacuum energy, is a convention and changes no physics.
> The two solutions differ in their matter content: Einstein's a balance of dust against $\Lambda$ on a knife edge, de Sitter's the empty limit where only $\Lambda$ remains.
> Arthur Eddington later summed up the pair: Einstein's was matter without motion, de Sitter's motion without matter [eddington1933].
> Peter Phillips built his story "Manna", in Astounding Science Fiction of February 1949, on that emptiness: the ghosts of two medieval monks leave a physicist a de Sitter version of the Riemann-Christoffel tensor, and he has to work out how anything can exist in a universe with neither matter nor radiation [phillips1949manna].

### history: `de_sitter/history[3]`

Source: `MFS/assets/data/metrics/de_sitter.json`, `history`, paragraph 4.

Before:

> De Sitter wrote his solution in static coordinates, and in that form it looked like a strange, empty, eternal world.
> The strangeness turned out to be an artifact of the coordinates.
> In 1925 Georges Lemaître, and independently Howard Robertson three years later, showed that the de Sitter universe could be rewritten in coordinates that are spatially flat and homogeneous, at the cost of a scale factor growing exponentially in time [lemaitre1925, robertson1928].
> When Edwin Hubble published his relation between the distances and velocities of the nebulae in 1929, he called its outstanding feature the possibility that it represents the de Sitter effect [hubble1929].
> The empty, static world was in fact expanding, and expanding forever; "motion without matter" had been hiding an accelerating universe.

After:

> De Sitter wrote his solution in static coordinates, and in that form it looked like a strange, empty, eternal world.
> The strangeness turned out to be an artifact of the coordinates.
> In 1925 Georges Lemaître, and independently Howard Robertson three years later, showed that the de Sitter universe could be rewritten in coordinates that are spatially flat and homogeneous, at the cost of a scale factor growing exponentially in time [lemaitre1925, robertson1928].
> When Edwin Hubble published his relation between the distances and velocities of the nebulae in 1929, he called its outstanding feature the possibility that it represents the de Sitter effect [hubble1929].
> The empty world looked static only in de Sitter's coordinates; it expands forever, and at an accelerating rate.

### history: `de_sitter/history[4]`

Source: `MFS/assets/data/metrics/de_sitter.json`, `history`, paragraph 5.

Before:

> That accelerating expansion is why de Sitter space refuses to stay a historical curiosity.
> It is the unique maximally symmetric solution with positive $\Lambda$, and it sits at both ends of modern cosmology: as the nearly exponential expansion of the inflationary epoch, and as the attractor at late times toward which a universe like our own, dominated by dark energy, appears to be heading [misner1973, wald1984].
> In Jonathan Lethem's novel As She Climbed Across the Table of 1997 a physicist opens a Farhi-Guth universe, a bubble of de Sitter space, and explains that to adhere it to the Schwarzschild space outside he had to develop a pair of anti-trapped surfaces in an asymptotically Minkowskian background [lethem1997table].
> Gary Gibbons and Stephen Hawking showed in 1977 that an observer in de Sitter space has an event horizon whose area can be interpreted as entropy, and that a particle detector the observer carries registers thermal radiation coming apparently from that horizon [gibbons1977].
> Robert Wald proved in 1983 that expanding homogeneous universes with a positive cosmological constant evolve exponentially toward the de Sitter solution, for every Bianchi type but IX, and for type IX when $\Lambda$ is large enough [wald1983].
> Two teams measuring distant supernovae of type Ia, one with Adam Riess and one with Saul Perlmutter as first author, found in 1998 and 1999 that the expansion of our own universe is accelerating, with a positive cosmological constant [riess1998, perlmutter1999].

After:

> The accelerating expansion is why cosmologists still study de Sitter space.
> It is the unique maximally symmetric solution with positive $\Lambda$, and it sits at both ends of modern cosmology: as the nearly exponential expansion of the inflationary epoch, and as the attractor at late times toward which a universe like our own, dominated by dark energy, appears to be heading [misner1973, wald1984].
> In Jonathan Lethem's novel As She Climbed Across the Table of 1997 a physicist opens a Farhi-Guth universe, a bubble of de Sitter space, and explains that to adhere it to the Schwarzschild space outside he had to develop a pair of anti-trapped surfaces in an asymptotically Minkowskian background [lethem1997table].
> Gary Gibbons and Stephen Hawking showed in 1977 that an observer in de Sitter space has an event horizon whose area can be interpreted as entropy, and that a particle detector the observer carries registers thermal radiation coming apparently from that horizon [gibbons1977].
> Robert Wald proved in 1983 that expanding homogeneous universes with a positive cosmological constant evolve exponentially toward the de Sitter solution, for every Bianchi type but IX, and for type IX when $\Lambda$ is large enough [wald1983].
> Two teams measuring distant supernovae of type Ia, one with Adam Riess and one with Saul Perlmutter as first author, found in 1998 and 1999 that the expansion of our own universe is accelerating, with a positive cosmological constant [riess1998, perlmutter1999].

### conventions: `de_sitter/convention[0]`

Source: `MFS/assets/data/metrics/de_sitter.json`, `convention`, paragraph 1.

Before:

> Coordinates of the static patch are $(t, r, \theta, \phi)$ with $t$ having dimensions of time and $r$ the areal radius, and factors of $c$ are kept explicit.
> The cosmological constant $\Lambda$ is positive and carries an inverse length squared, and the two lengths it sets are the horizon radius $\sqrt{3/\Lambda}$ and the Hubble radius $c/H$ with $H = c\sqrt{\Lambda/3}$, which are the same length.
> Components are taken in the chart $x^0 = ct$, where a chart component of the metric carries no dimensions, and the dots in the geodesic equations are derivatives of that same chart with respect to an affine parameter, so $\dot{t}$ means $d(ct)/d\lambda$ and every term of every equation carries the dimensions of its left hand side.
> The geodesic equations are $\ddot{x}^\mu + \Gamma^\mu{}_{\nu\rho}\dot{x}^\nu\dot{x}^\rho = 0$ with $\Gamma$ the Christoffel symbols of that same chart.
> The signature is $(-,+,+,+)$, the Riemann tensor is $R^\mu{}_{\nu\rho\sigma} = \partial_\rho\Gamma^\mu{}_{\nu\sigma} - \partial_\sigma\Gamma^\mu{}_{\nu\rho} + \Gamma^\mu{}_{\rho\lambda}\Gamma^\lambda{}_{\nu\sigma} - \Gamma^\mu{}_{\sigma\lambda}\Gamma^\lambda{}_{\nu\rho}$, and the Ricci tensor is the standard contraction $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$.

After:

> For the static patch we use coordinates $(t, r, \theta, \phi)$, with $t$ carrying dimensions of time and $r$ the areal radius, and keep factors of $c$ explicit.
> The cosmological constant $\Lambda$ is positive and carries an inverse length squared; it sets the horizon radius $\sqrt{3/\Lambda}$ and the Hubble radius $c/H$ with $H = c\sqrt{\Lambda/3}$, and the two are the same length.
> We take components in the chart $x^0 = ct$, where a chart component of the metric carries no dimensions, and a dot in the geodesic equations is a derivative of that chart's coordinates with respect to an affine parameter, so $\dot{t}$ means $d(ct)/d\lambda$, and every term of every equation carries the dimensions of its left hand side.
> The geodesic equations read $\ddot{x}^\mu + \Gamma^\mu{}_{\nu\rho}\dot{x}^\nu\dot{x}^\rho = 0$, with $\Gamma$ the Christoffel symbols of the same chart.
> We work in the signature $(-,+,+,+)$, with the Riemann tensor $R^\mu{}_{\nu\rho\sigma} = \partial_\rho\Gamma^\mu{}_{\nu\sigma} - \partial_\sigma\Gamma^\mu{}_{\nu\rho} + \Gamma^\mu{}_{\rho\lambda}\Gamma^\lambda{}_{\nu\sigma} - \Gamma^\mu{}_{\sigma\lambda}\Gamma^\lambda{}_{\nu\rho}$ and the Ricci tensor its standard contraction, $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$.

### conventions: `de_sitter/convention[1]`

Source: `MFS/assets/data/metrics/de_sitter.json`, `convention`, paragraph 2.

Before:

> This is a vacuum with a cosmological constant rather than an empty one, so the Ricci tensor does not vanish: $R_{\mu\nu} = \Lambda g_{\mu\nu}$, $R = 4\Lambda$ and $G_{\mu\nu} = -\Lambda g_{\mu\nu}$, which is $G_{\mu\nu} + \Lambda g_{\mu\nu} = 0$ written the other way round.
> The spacetime is maximally symmetric, its Riemann tensor being $R_{\mu\nu\rho\sigma} = \dfrac{\Lambda}{3}\left(g_{\mu\rho}g_{\nu\sigma} - g_{\mu\sigma}g_{\nu\rho}\right)$, so the Weyl tensor vanishes identically in both charts: every trace free piece of Riemann is empty and the curvature is carried entirely by the Ricci tensor.
> The Kretschmann scalar is the constant $8\Lambda^2/3$ everywhere, so neither the horizon of the static patch nor the edge of the flat slicing is a singularity; both are failures of their chart.
> The second chart is the flat slicing $(t, x, y, z)$, whose parameter is the Hubble rate $H$ rather than $\Lambda$, related by $3H^2/c^2 = \Lambda$; its invariants are the same, $R = 12H^2/c^2 = 4\Lambda$ and $K = 24H^4/c^4 = 8\Lambda^2/3$, and covers only the expanding half of the spacetime.

After:

> De Sitter space is a vacuum with a cosmological constant, so its Ricci tensor does not vanish: $R_{\mu\nu} = \Lambda g_{\mu\nu}$, $R = 4\Lambda$, and $G_{\mu\nu} = -\Lambda g_{\mu\nu}$, which is the field equation $G_{\mu\nu} + \Lambda g_{\mu\nu} = 0$ rearranged.
> The spacetime is maximally symmetric, with the Riemann tensor $R_{\mu\nu\rho\sigma} = \dfrac{\Lambda}{3}\left(g_{\mu\rho}g_{\nu\sigma} - g_{\mu\sigma}g_{\nu\rho}\right)$, so the Weyl tensor vanishes identically in both charts: every trace free piece of Riemann is empty, and the Ricci tensor carries the whole curvature.
> The Kretschmann scalar is the constant $8\Lambda^2/3$ everywhere, so neither the horizon of the static patch nor the edge of the flat slicing is a singularity; each is a failure of its chart.
> The second chart is the flat slicing $(t, x, y, z)$, whose parameter is the Hubble rate $H$ in place of $\Lambda$, related by $3H^2/c^2 = \Lambda$; its invariants are the same, $R = 12H^2/c^2 = 4\Lambda$ and $K = 24H^4/c^4 = 8\Lambda^2/3$, and it covers only the expanding half of the spacetime.

### conformal diagram caption: `conformal/de_sitter/static.caption[0]`

Source: `_tools/derivations/conformal.py` line 2274, `CAPTIONS ("de_sitter", "static")`.

Before:

> This is the whole of de Sitter spacetime, the hyperboloid $-X_0^2 + X_1^2 + \dots + X_4^2 = L^2$ with $L = \sqrt{3/\Lambda}$, and each point of the diagram stands for a sphere.
> Its global coordinates give the metric $\frac{L^2}{\cos^2 T}(-dT^2 + d\chi^2 + \sin^2\chi\,d\Omega^2)$ on the square $|T| < \pi/2$, $0 \le \chi \le \pi$, with $X = \chi$ across.

After:

> This is the whole of de Sitter spacetime, the hyperboloid $-X_0^2 + X_1^2 + \dots + X_4^2 = L^2$ with $L = \sqrt{3/\Lambda}$, and each point of the diagram stands for a sphere.
> In its global coordinates the metric is $\frac{L^2}{\cos^2 T}(-dT^2 + d\chi^2 + \sin^2\chi\,d\Omega^2)$ on the square $|T| < \pi/2$, $0 \le \chi \le \pi$, with $X = \chi$ across.

## Ellis-Bronnikov Traversable Wormhole

### history: `ellis_bronnikov/history[0]`

Source: `MFS/assets/data/metrics/ellis_bronnikov.json`, `history`, paragraph 1.

Before:

> In 1973 two physicists on two sides of an iron curtain, neither aware of the other, found the same geometry.
> Homer Ellis, working at the University of Colorado, wrote down what he called a "drainhole", a solution to the Einstein equations coupled to a massless phantom scalar field, characterized by a throat connecting two asymptotically flat regions [ellis1973].
> He wanted a model of a particle free of the point singularity at the center of the Schwarzschild manifold, and opened that singularity up with a scalar field coupled to the geometry with the opposite of the usual sign [ellis1973].
> In his general drainhole an ether flows through the hole and the two sides carry masses of opposite sign; with the ether at rest, the case of zero redshift function, the particle has no mass on either side, yet the hole stays open [ellis1973].
> His paper, first submitted in December 1969, appeared in January 1973 [ellis1973].

After:

> Two physicists on two sides of an iron curtain, neither aware of the other, found the same geometry in 1973.
> Homer Ellis, working at the University of Colorado, wrote down what he called a "drainhole", a solution to the Einstein equations coupled to a massless phantom scalar field, characterized by a throat connecting two asymptotically flat regions [ellis1973].
> He wanted a model of a particle free of the point singularity at the center of the Schwarzschild manifold, and opened that singularity up with a scalar field coupled to the geometry with the opposite of the usual sign [ellis1973].
> In his general drainhole an ether flows through the hole and the two sides carry masses of opposite sign; with the ether at rest, the case of zero redshift function, the particle has no mass on either side, yet the hole stays open [ellis1973].
> His paper, first submitted in December 1969, appeared in January 1973 [ellis1973].

### history: `ellis_bronnikov/history[2]`

Source: `MFS/assets/data/metrics/ellis_bronnikov.json`, `history`, paragraph 3.

Before:

> Neither author was thinking primarily about traversability; that framing would come fifteen years later.
> The metric Morris and Thorne set out in the second box of their 1988 paper as the simplest traversable wormhole is exactly Ellis's, taken with small changes from the final exam of a course Thorne taught at Caltech in the autumn of 1985 [morris1988].
> Thorne returned to it in 2015 with Oliver James, Eugénie von Tunzelmann, and Paul Franklin of the effects studio Double Negative, who used it as the starting point for the wormhole of the film Interstellar and noted that Morris and Thorne, unaware of Ellis's paper, had not credited him, for which they apologized [james2015].

After:

> Neither author was thinking primarily about traversability; that framing would come fifteen years later.
> The metric Michael Morris and Kip Thorne set out in the second box of their 1988 paper as the simplest traversable wormhole is exactly Ellis's, taken with small changes from the final exam of a course Thorne taught at Caltech in the autumn of 1985 [morris1988].
> Thorne returned to it in 2015 with Oliver James, Eugénie von Tunzelmann, and Paul Franklin of the effects studio Double Negative, who used it as the starting point for the wormhole of the film Interstellar and noted that Morris and Thorne, unaware of Ellis's paper, had not credited him, for which they apologized [james2015].

### history: `ellis_bronnikov/history[3]`

Source: `MFS/assets/data/metrics/ellis_bronnikov.json`, `history`, paragraph 4.

Before:

> The solution Ellis and Bronnikov found was remarkable on its own terms.
> The metric is geodesically complete: there is no singularity anywhere, and the two asymptotically flat ends are joined by a smooth throat of radius $\ell$.
> The price is a phantom scalar field, a matter source whose kinetic term carries the wrong sign and violates the null energy condition.
> This was the first clear instance of what Morris and Thorne would later identify as the exotic matter requirement [morris1988]: the wormhole geometry is geometrically clean but physically costly.

After:

> The solution Ellis and Bronnikov found was remarkable on its own terms.
> The metric is geodesically complete: there is no singularity anywhere, and the two asymptotically flat ends are joined by a smooth throat of radius $\ell$.
> The price is a phantom scalar field, a matter source whose kinetic term carries the wrong sign and violates the null energy condition.
> This was the first clear instance of what Morris and Thorne would later identify as the exotic matter requirement [morris1988]: a wormhole with a smooth throat needs a source that violates the null energy condition.

### history: `ellis_bronnikov/history[4]`

Source: `MFS/assets/data/metrics/ellis_bronnikov.json`, `history`, paragraph 5.

Before:

> The solution sat quietly in the literature until the late 1980s, when the traversable wormhole program brought it renewed attention [visser1995].
> Cristian Armendáriz-Picón argued in 2002 that the massless wormholes of a phantom scalar field are stable [armendarizpicon2002].
> The same year Hisa-aki Shinkai and Sean Hayward evolved the Ellis wormhole numerically and found it unstable: a pulse carrying net negative energy makes the throat blow up into an inflating universe, positive energy makes it collapse to a black hole, and ordinary matter, a traveller included, always brings collapse [shinkai2002].
> José Antonio González, Francisco Siddhartha Guzmán, and Olivier Sarbach proved in 2009 that every static wormhole of a massless ghost scalar field has exactly one unstable mode, growing on the time light takes to cross the throat [gonzalez2009a, gonzalez2009b].

After:

> The solution drew little attention until the late 1980s, when the program of traversable wormholes brought it renewed attention [visser1995].
> Cristian Armendáriz-Picón argued in 2002 that the massless wormholes of a phantom scalar field are stable [armendarizpicon2002].
> The same year Hisa-aki Shinkai and Sean Hayward evolved the Ellis wormhole numerically and found it unstable: a pulse carrying net negative energy makes the throat blow up into an inflating universe, positive energy makes it collapse to a black hole, and ordinary matter, a traveller included, always brings collapse [shinkai2002].
> José Antonio González, Francisco Siddhartha Guzmán, and Olivier Sarbach proved in 2009 that every static wormhole of a massless ghost scalar field has exactly one unstable mode, growing on the time light takes to cross the throat [gonzalez2009a, gonzalez2009b].

### history: `ellis_bronnikov/history[5]`

Source: `MFS/assets/data/metrics/ellis_bronnikov.json`, `history`, paragraph 6.

Before:

> Fumio Abe worked out in 2010 the light curves an Ellis wormhole would produce as a gravitational microlens, with a dip of about four percent just outside the crossing of the Einstein ring [abe2010].
> Ryuichi Takahashi and Hideki Asada found no lensing by Ellis wormholes among about fifty thousand quasars of the Sloan Digital Sky Survey Quasar Lens Search and set an upper bound on their abundance [takahashi2013].
> It is now a standard example in graduate GR, simple enough to work with by hand, sharp enough to expose the essential tension between wormhole geometry and the energy conditions: the geometry asks for a bridge, and the field equations send the bill in exotic matter.

After:

> Fumio Abe worked out in 2010 the light curves an Ellis wormhole would produce as a gravitational microlens, with a dip of about four percent just outside the crossing of the Einstein ring [abe2010].
> Ryuichi Takahashi and Hideki Asada found no lensing by Ellis wormholes among about fifty thousand quasars of the Sloan Digital Sky Survey Quasar Lens Search and set an upper bound on their abundance [takahashi2013].
> The Ellis wormhole is now a standard example in graduate courses on general relativity: simple enough to work with by hand, it shows that a throat joining two flat ends needs a source that violates the energy conditions.

### conventions: `ellis_bronnikov/convention[0]`

Source: `MFS/assets/data/metrics/ellis_bronnikov.json`, `convention`, paragraph 1.

Before:

> The time coordinate is $ct$.
> The coordinate $r$ is the proper radial distance, $r \in (-\infty,\infty)$, with the wormhole throat at $r = 0$.
> The parameter $\ell > 0$ is the throat radius.
> The Ricci tensor is the standard contraction $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$, on which $G_{tt} = -\dfrac{\ell^2}{(r^2+\ell^2)^2}$ is negative, which is the negative energy density the throat needs to stay open.

After:

> We take $ct$ as the time coordinate.
> The coordinate $r$ is the proper radial distance, $r \in (-\infty,\infty)$, with the wormhole throat at $r = 0$.
> The parameter $\ell > 0$ is the radius of the throat.
> With the Ricci tensor the standard contraction, $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$, the component $G_{tt} = -\dfrac{\ell^2}{(r^2+\ell^2)^2}$ is negative, and this is the negative energy density the throat needs to stay open.

### spacetime diagram caption: `diagrams/ellis_bronnikov/systems/spherical/radial.caption[1]`

Source: `_tools/derivations/null_rays.py` line 606, `CAPTIONS ("ellis_bronnikov", "spherical", "radial")`.

Before:

> The throat is where the spheres are smallest.
> The angular part of the metric is $g_{\theta\theta} = r^2 + \ell^2$, so the areal radius $R = \sqrt{r^2 + \ell^2}$ takes its least value, $\ell$, at $r = 0$, where $\partial_r R = 0$.
> The faint vertical lines are the spheres $R = 1.5\ell$, $2\ell$ and $3\ell$, one of each size on either side of the throat.
> The Kretschmann scalar $12\ell^4/(r^2 + \ell^2)^4$ is finite everywhere.

After:

> The throat is where the spheres are smallest.
> The angular part of the metric is $g_{\theta\theta} = r^2 + \ell^2$, so the areal radius $R = \sqrt{r^2 + \ell^2}$ takes its least value, $\ell$, at $r = 0$, where $\partial_r R = 0$.
> The faint vertical lines are the spheres $R = 1.5\ell$, $2\ell$, and $3\ell$, one of each size on either side of the throat.
> The Kretschmann scalar $12\ell^4/(r^2 + \ell^2)^4$ is finite everywhere.

### embedding diagram caption: `embedding/ellis_bronnikov/wormhole.caption[0]`

Source: `_tools/derivations/embedding.py` line 3541, `CAPTIONS ("ellis_bronnikov", "wormhole")`.

Before:

> This is the equatorial plane $\theta = \pi/2$ of the Ellis-Bronnikov wormhole at one moment of $t$, drawn as a surface in flat space so that every distance along it is the metric distance.
> Its $r$ is the proper distance from the throat, running from $-\infty$ on one side to $\infty$ on the other, so the circles of constant $r$ stand at equal steps along the surface, and the circle at $r$ has circumference $2\pi\sqrt{r^2 + \ell^2}$.
> The surface that carries both is the catenoid $\sqrt{r^2 + \ell^2} = \ell\cosh(z/\ell)$, one piece through the throat at $r = 0$, where the circles are smallest and the surface stands vertical.

After:

> This is the equatorial plane $\theta = \pi/2$ of the Ellis-Bronnikov wormhole at one moment of $t$, drawn as a surface in flat space so that every distance along it is the metric distance.
> Its $r$ is the proper distance from the throat, running from $-\infty$ on one side to $\infty$ on the other, so the circles of constant $r$ stand at equal steps along the surface, and the circle at $r$ has circumference $2\pi\sqrt{r^2 + \ell^2}$.
> Both hold on the catenoid $\sqrt{r^2 + \ell^2} = \ell\cosh(z/\ell)$, one piece through the throat at $r = 0$, where the circles are smallest and the surface stands vertical.

### embedding diagram caption: `embedding/ellis_bronnikov/wormhole.caption[1]`

Source: `_tools/derivations/embedding.py` line 3548, `CAPTIONS ("ellis_bronnikov", "wormhole")`.

Before:

> It is the surface the Morris-Thorne wormhole draws, and the same metric: Michael Morris and Kip Thorne set it out in 1988 as the simplest traversable wormhole, with the shape function $b = \ell^2/R$ of the areal radius $R = \sqrt{r^2 + \ell^2}$, unaware that Homer Ellis and Kirill Bronnikov had each found it in 1973.
> The areal radius turns back at the throat, so it covers one side at a time; the proper $r$ runs straight through.

After:

> The same surface, with the same metric, is the Morris-Thorne wormhole of Michael Morris and Kip Thorne, who set it out in 1988 as the simplest traversable wormhole, with the shape function $b = \ell^2/R$ of the areal radius $R = \sqrt{r^2 + \ell^2}$, unaware that Homer Ellis and Kirill Bronnikov had each found it in 1973.
> The areal radius turns back at the throat, so it covers one side at a time, while the proper $r$ runs straight through.

## Friedmann-Lemaître-Robertson-Walker Expanding Universe

### history: `frw/history[1]`

Source: `MFS/assets/data/metrics/frw.json`, `history`, paragraph 2.

Before:

> Five years later, Georges Lemaître independently derived the expanding solution and connected it to the growing observational record [lemaitre1927]: Vesto Slipher had been measuring the redshifts of spiral nebulae since 1913 [slipher1913], Carl Wilhelm Wirtz had noted a systematic recession trend in 1918 and again in 1922 [wirtz1918, wirtz1922], and Edwin Hubble would quantify the relationship between distance and recession in 1929 [hubble1929].
> Lemaître's 1927 paper had already estimated the rate of expansion, but that passage is missing from the English translation the Monthly Notices printed in 1931 [lemaitre1931]; Mario Livio found in 2011 that Lemaître had made the translation himself and left out his provisional discussion of radial velocities [livio2011].
> When Lemaître discussed his paper with Einstein at the Solvay conference of 1927, Einstein, as Lemaître recalled, made favourable technical remarks and then called the physics abominable [nussbaumer2014].
> Arthur Eddington, who had been examining with G. C. McVittie whether Einstein's static world is stable, learned of Lemaître's paper and showed in 1930 that the static world is unstable [eddington1930].
> Einstein gave up his static universe in April 1931, in a note to the Prussian Academy that took Friedmann's equations as its starting point and dropped the cosmological constant as theoretically unsatisfactory [oraifeartaigh2014].

After:

> Five years later, Georges Lemaître independently derived the expanding solution and connected it to the growing observational record [lemaitre1927]: Vesto Slipher had been measuring the redshifts of spiral nebulae since 1913 [slipher1913], Carl Wilhelm Wirtz had noted a systematic recession trend in 1918 and again in 1922 [wirtz1918, wirtz1922], and Edwin Hubble would quantify the relationship between distance and recession in 1929 [hubble1929].
> Lemaître's 1927 paper had already estimated the rate of expansion, but that passage is missing from the English translation printed in the Monthly Notices in 1931 [lemaitre1931]; Mario Livio found in 2011 that Lemaître had made the translation himself and left out his provisional discussion of radial velocities [livio2011].
> When Lemaître discussed his paper with Einstein at the Solvay conference of 1927, Einstein, as Lemaître recalled, made favourable technical remarks and then called the physics abominable [nussbaumer2014].
> Arthur Eddington, who had been examining with G. C. McVittie whether Einstein's static world is stable, learned of Lemaître's paper and showed in 1930 that the static world is unstable [eddington1930].
> Einstein gave up his static universe in April 1931, in a note to the Prussian Academy that started from Friedmann's equations and dropped the cosmological constant as theoretically unsatisfactory [oraifeartaigh2014].

### history: `frw/history[2]`

Source: `MFS/assets/data/metrics/frw.json`, `history`, paragraph 3.

Before:

> Howard Robertson and Arthur Walker placed the whole framework on rigorous geometric footing in the 1930s, proving that this metric is the unique spacetime consistent with spatial homogeneity and isotropy [robertson1935, walker1937], what physicists now call the cosmological principle.
> Homogeneity and isotropy together make each spatial slice maximally symmetric, carrying the largest symmetry group a geometry of three dimensions can.
> The spacetime itself is not, because the scale factor $a(t)$ changes with time and breaks the temporal symmetries that the maximally symmetric vacua, Minkowski and de Sitter and anti-de Sitter, retain; FLRW is the cosmology one gets by demanding maximal symmetry of space alone, and letting time do the rest.

After:

> Howard Robertson and Arthur Walker placed the whole framework on rigorous geometric footing in the 1930s, proving that this metric is the unique spacetime consistent with spatial homogeneity and isotropy [robertson1935, walker1937], what physicists now call the cosmological principle.
> Homogeneity and isotropy together make each spatial slice maximally symmetric, carrying the largest symmetry group a geometry of three dimensions can.
> The spacetime itself is not, because the scale factor $a(t)$ changes with time and breaks the temporal symmetries that the maximally symmetric vacua, Minkowski, de Sitter, and anti-de Sitter, retain; FLRW is the cosmology of maximal symmetry in space alone, with the scale factor free to change in time.

### conventions: `frw/convention[0]`

Source: `MFS/assets/data/metrics/frw.json`, `convention`, paragraph 1.

Before:

> Coordinates are cosmic time $t$ and a comoving radial coordinate $r$, the areal radius of the sphere through a point divided by $a$, with $\theta \in [0,\pi]$ and $\phi \in [0,2\pi)$.
> The scale factor $a(t) > 0$ describes the expansion of the universe.
> The curvature parameter $k \in \{-1, 0, +1\}$ encodes the spatial geometry: $k=-1$ (open), $k=0$ (flat), $k=+1$ (closed), in units of the inverse square of the curvature radius, since $r$ is a length.
> Dots denote derivatives with respect to $t$.
> The Ricci tensor is the standard contraction $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$, so the field equations read $G_{\mu\nu} = 8\pi G\,T_{\mu\nu}/c^4$ and $G_{tt} = 3\dfrac{\dot{a}^2+k}{a^2}$ is the Friedmann equation with a positive energy density on the other side of it.
> The spacetime is conformally flat, so its Weyl tensor vanishes identically.

After:

> We use cosmic time $t$ and a comoving radial coordinate $r$, the areal radius of the sphere through a point divided by $a$, with $\theta \in [0,\pi]$ and $\phi \in [0,2\pi)$.
> The scale factor $a(t) > 0$ describes the expansion of the universe.
> The curvature parameter $k \in \{-1, 0, +1\}$ sets the spatial geometry, $k=-1$ (open), $k=0$ (flat), and $k=+1$ (closed), in units of the inverse square of the curvature radius, since $r$ is a length.
> Dots denote derivatives with respect to $t$.
> The Ricci tensor is the standard contraction, $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$, so the field equations read $G_{\mu\nu} = 8\pi G\,T_{\mu\nu}/c^4$, and $G_{tt} = 3\dfrac{\dot{a}^2+k}{a^2}$ is the Friedmann equation with a positive energy density on the other side of it.
> The spacetime is conformally flat, so its Weyl tensor vanishes identically.

### conventions: `frw/convention[1]`

Source: `MFS/assets/data/metrics/frw.json`, `convention`, paragraph 2.

Before:

> For $k \le 0$ the coordinates cover the whole of each spatial slice as $r$ runs over $[0, \infty)$.
> For $k = +1$ they do not.
> The slice is then a three sphere $S^3$ of radius $a/\sqrt{k}$, and $g_{rr} = a^2/(1-kr^2)$ diverges at $r = 1/\sqrt{k}$, which is its equator, so $r$ ends there having covered one hemisphere, and the other hemisphere would repeat the same values of $r$.
> Writing $r = \sin\chi/\sqrt{k}$ turns $dr^2/(1-kr^2)$ into $d\chi^2/k$ and covers the whole three sphere as the polar angle $\chi$ runs over $[0, \pi]$, which is how Friedmann wrote the closed case in 1922, how Misner, Thorne and Wheeler and Wald write it in their textbooks, and how the interior of the Oppenheimer-Snyder star is written.

After:

> For $k \le 0$ the coordinates cover the whole of each spatial slice as $r$ runs over $[0, \infty)$.
> For $k = +1$ they do not.
> The slice is then a three sphere $S^3$ of radius $a/\sqrt{k}$, and $g_{rr} = a^2/(1-kr^2)$ diverges at $r = 1/\sqrt{k}$, its equator, so $r$ ends there having covered one hemisphere, and the other hemisphere would repeat the same values of $r$.
> Writing $r = \sin\chi/\sqrt{k}$ turns $dr^2/(1-kr^2)$ into $d\chi^2/k$ and covers the whole three sphere as the polar angle $\chi$ runs over $[0, \pi]$.
> Friedmann wrote the closed case this way in 1922; Misner, Thorne, and Wheeler write it so in their textbook, as Wald does in his, and the interior of the Oppenheimer-Snyder star takes the same form.

### conformal diagram note: `conformal/frw/flat.input`

Source: `_tools/derivations/conformal.py` line 1734, `in frw() input=`.

Shown in: `conformal/frw/flat.input`, `conformal/frw/closed.input`, `conformal/frw/open.input`.

Before:

> Dust: $a \propto \eta^2$, $1 - \cos\eta$ and $\cosh\eta - 1$ for $k = 0$, $+1$ and $-1$, each solved from this spacetime's own $G^r{}_r = 0$, with lengths in units of $1/\sqrt{|k|}$ where $k$ is not zero.

After:

> Dust: $a \propto \eta^2$, $1 - \cos\eta$, and $\cosh\eta - 1$ for $k = 0$, $+1$, and $-1$, each solved from this spacetime's own $G^r{}_r = 0$, with lengths in units of $1/\sqrt{|k|}$ where $k$ is not zero.

## Gödel Rotating Universe

### history: `godel/history[0]`

Source: `MFS/assets/data/metrics/godel.json`, `history`, paragraph 1.

Before:

> Kurt Gödel presented his rotating universe in the July 1949 issue of Reviews of Modern Physics, the issue dedicated to Einstein's seventieth birthday: a gift, as it were, that unsettled the honoree [godel1949].
> The solution describes a universe filled with pressureless dust in uniform rotation, homogeneous in space and time, with no expansion.
> Its cosmological constant is negative, tuned to the density of the dust.
> About any world line of the dust the light cones tip as one moves outward, until at a critical radius they become tangent to a closed null curve, and beyond that radius closed timelike curves exist through every event: trajectories that loop back to their own past.
> Rotating dust had been solved before, by Cornelius Lanczos in 1924 and by Willem Jacob van Stockum in 1937 [lanczos1924, vanstockum1937], but Gödel's paper cites neither, and his universe is homogeneous where theirs is a cylinder turning about one axis.

After:

> Kurt Gödel presented his rotating universe in the July 1949 issue of Reviews of Modern Physics, the issue dedicated to Einstein's seventieth birthday: a gift, as it were, that unsettled the honoree [godel1949].
> The solution describes a universe filled with pressureless dust in uniform rotation, homogeneous in space and time, with no expansion.
> Its cosmological constant is negative, tuned to the density of the dust.
> About any world line of the dust the light cones tip as one moves outward, until at a critical radius they become tangent to a closed null curve, and beyond that radius closed timelike curves exist through every event: trajectories that loop back to their own past.
> Rotating dust had been solved before, by Cornelius Lanczos in 1924 and by Willem Jacob van Stockum in 1937 [lanczos1924, vanstockum1937], but Gödel cites neither, and his universe is homogeneous where theirs is a cylinder turning about one axis.

### history: `godel/history[1]`

Source: `MFS/assets/data/metrics/godel.json`, `history`, paragraph 2.

Before:

> Gödel stated plainly in the paper that it is theoretically possible in these worlds to travel into the past [godel1949].
> Gödel said more in the essay he wrote for the same birthday, in the volume Paul Arthur Schilpp gathered in Einstein's honor [godel1949idealistic].
> By a round trip on a rocket ship in a sufficiently wide curve, he wrote, one could reach any region of the past, present and future, though the velocities and the fuel such a trip would need were far beyond anything that could be expected ever to become a practical possibility.
> The trip was not his point.

After:

> Gödel stated plainly in the paper that it is theoretically possible in these worlds to travel into the past [godel1949].
> Gödel said more in the essay he wrote for the same birthday, in the volume Paul Arthur Schilpp gathered in Einstein's honor [godel1949idealistic].
> By a round trip on a rocket ship in a sufficiently wide curve, he wrote, one could reach any region of the past, present, and future, though the velocities and the fuel such a trip would need were far beyond anything that could be expected ever to become a practical possibility.
> The trip was not his point.

### history: `godel/history[2]`

Source: `MFS/assets/data/metrics/godel.json`, `history`, paragraph 3.

Before:

> Relativity already seemed to him to prove the view of Parmenides, Kant and the modern idealists that change is an illusion or an appearance due to our mode of perception.
> James Jeans had concluded from the cosmic time of the known cosmological solutions that there was no reason to abandon an absolute time lapsing objectively, and Gödel's answer was his own solutions, since the mere compatibility with the laws of nature of worlds with no distinguished absolute time throws light on the meaning of time in worlds where one can be defined.
> Einstein answered in the same volume that the essay was an important contribution to the analysis of the concept of time, that the problem had disturbed him already while he was building the general theory, and that it would be interesting to weigh whether such solutions are to be excluded on physical grounds [einstein1949remarks].
> At the International Congress of Mathematicians in 1950 Gödel went on to rotating universes that also expand, as ours does, and showed that for them the absence of closed timelike lines is equivalent to the existence of a world time [godel1952].

After:

> Relativity already seemed to him to prove the view of Parmenides, Kant, and the modern idealists that change is an illusion or an appearance due to our mode of perception.
> James Jeans had concluded from the cosmic time of the known cosmological solutions that there was no reason to abandon an absolute time lapsing objectively, and Gödel's answer was his own solutions, since the mere compatibility with the laws of nature of worlds with no distinguished absolute time throws light on the meaning of time in worlds where one can be defined.
> Einstein answered in the same volume that the essay was an important contribution to the analysis of the concept of time, that the problem had disturbed him already while he was building the general theory, and that it would be interesting to weigh whether such solutions are to be excluded on physical grounds [einstein1949remarks].
> At the International Congress of Mathematicians in 1950 Gödel went on to rotating universes that also expand, as ours does, and showed that for them the absence of closed timelike lines is equivalent to the existence of a world time [godel1952].

### history: `godel/history[5]`

Source: `MFS/assets/data/metrics/godel.json`, `history`, paragraph 6.

Before:

> Neal Stephenson names it outright in Anathem: when the avout of Arbre and the Laterran Jules Verne Durand work out between them that a spaceship in a universe that happens to be rotating, if it travels far and fast enough, will travel backwards in time, a result the avout know and had always thought little more than a curiosity, Durand tells them that on his world it was discovered by "a kind of Saunt named Gödel: a friend of the Saunt who had earlier discovered geometrodynamics" [stephenson2008].
> The Gödel universe is not a serious model of the cosmos we inhabit, since it does not expand and we observe expansion, but it remains one of the most important exact solutions in the theory.
> It is the solution that put the causal structure of spacetime on trial, and Einstein's field equations permit it, closed timelike curves and all.

After:

> Neal Stephenson names it outright in Anathem: when the avout of Arbre and the Laterran Jules Verne Durand work out between them that a spaceship in a universe that happens to be rotating, if it travels far and fast enough, will travel backwards in time, a result the avout know and had always thought little more than a curiosity, Durand tells them that on his world it was discovered by "a kind of Saunt named Gödel: a friend of the Saunt who had earlier discovered geometrodynamics" [stephenson2008].
> The Gödel universe is not a serious model of the cosmos we inhabit, since it does not expand and we observe expansion, but it remains one of the most important exact solutions in the theory.
> It showed that Einstein's field equations permit closed timelike curves through every event of a universe filled with rotating dust.

### conventions: `godel/convention[0]`

Source: `MFS/assets/data/metrics/godel.json`, `convention`, paragraph 1.

Before:

> Coordinates $(t, x, y, z)$ all range over $\mathbb{R}$.
> The parameter $\omega$ is the angular velocity of the dust.
> The Ricci tensor is the standard contraction $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$, on which $R_{\mu\nu} = 2\omega^2u_\mu u_\nu$ with $u^\mu$ the four velocity of the dust.
> Units with $G = c = 1$ are used throughout, and the chart coordinate is $x^0 = ct$, which those units make $t$ itself.
> The argument of $e^x$ forces $x$ to be dimensionless and the metric forces $t$, $y$ and $z$ to be dimensionless with it, so the one length of the solution is $1/\omega$.
> The geodesic equations are $\ddot{x}^\mu + \Gamma^\mu{}_{\nu\rho}\dot{x}^\nu\dot{x}^\rho = 0$ with $\Gamma$ the Christoffel symbols of that same chart, and the dots are derivatives of those chart coordinates with respect to an affine parameter $\lambda$, so $\dot{t}$ means $d(ct)/d\lambda$ and every term of every equation carries the dimensions of its left hand side.

After:

> The coordinates $(t, x, y, z)$ all range over $\mathbb{R}$.
> The parameter $\omega$ is the angular velocity of the dust.
> With the Ricci tensor the standard contraction, $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$, we have $R_{\mu\nu} = 2\omega^2u_\mu u_\nu$, with $u^\mu$ the four velocity of the dust.
> We work in units with $G = c = 1$ throughout, and the chart coordinate $x^0 = ct$ is then $t$ itself.
> The argument of $e^x$ forces $x$ to be dimensionless, and the metric forces $t$, $y$, and $z$ to be dimensionless with it, so the one length of the solution is $1/\omega$.
> The geodesic equations read $\ddot{x}^\mu + \Gamma^\mu{}_{\nu\rho}\dot{x}^\nu\dot{x}^\rho = 0$, with $\Gamma$ the Christoffel symbols of the same chart, and a dot is a derivative of those chart coordinates with respect to an affine parameter $\lambda$, so $\dot{t}$ means $d(ct)/d\lambda$, and every term of every equation carries the dimensions of its left hand side.

### conventions: `godel/convention[1]`

Source: `MFS/assets/data/metrics/godel.json`, `convention`, paragraph 2.

Before:

> Gödel's own coordinates of 1949, $(t, r, \phi, z)$, are cylindrical about one world line of the dust, with $t$, $r$ and $z$ dimensionless as before and $\phi$ an angle.
> The Cartesian coordinates are carried to them by $e^x = \cosh 2r + \cos\phi\sinh 2r$, $ye^x = \sqrt{2}\sin\phi\sinh 2r$, $\tan\tfrac{1}{2}\left(\phi + (t_x - 2t)/\sqrt{2}\right) = e^{-2r}\tan\tfrac{1}{2}\phi$ and $z_x = 2z$, where $t_x$ and $z_x$ are the Cartesian $t$ and $z$, and the line element becomes $\dfrac{2}{\omega^2}\left(-dt^2 + dr^2 + \left(\sinh^2 r - \sinh^4 r\right)d\phi^2 - 2\sqrt{2}\sinh^2 r\,dt\,d\phi + dz^2\right)$, with the same curvature, $K = 12\omega^4$.
> The circles of constant $t$, $r$ and $z$ have $g_{\phi\phi} = 2(1 - \sinh^2 r)\sinh^2 r/\omega^2$, so they are spacelike inside $\sinh r = 1$, null on it and timelike beyond it, where each is a closed timelike curve.
> Every world line of the dust is equivalent to every other, so such a curve passes through every event.

After:

> Kurt Gödel's own coordinates of 1949, $(t, r, \phi, z)$, are cylindrical about one world line of the dust, with $t$, $r$, and $z$ dimensionless as before and $\phi$ an angle.
> The Cartesian coordinates are carried to them by $e^x = \cosh 2r + \cos\phi\sinh 2r$, $ye^x = \sqrt{2}\sin\phi\sinh 2r$, $\tan\tfrac{1}{2}\left(\phi + (t_x - 2t)/\sqrt{2}\right) = e^{-2r}\tan\tfrac{1}{2}\phi$, and $z_x = 2z$, where $t_x$ and $z_x$ are the Cartesian $t$ and $z$, and the line element becomes $\dfrac{2}{\omega^2}\left(-dt^2 + dr^2 + \left(\sinh^2 r - \sinh^4 r\right)d\phi^2 - 2\sqrt{2}\sinh^2 r\,dt\,d\phi + dz^2\right)$, with the same curvature, $K = 12\omega^4$.
> The circles of constant $t$, $r$, and $z$ have $g_{\phi\phi} = 2(1 - \sinh^2 r)\sinh^2 r/\omega^2$, so they are spacelike inside $\sinh r = 1$, null on it, and timelike beyond it, where each is a closed timelike curve.
> Every world line of the dust is equivalent to every other, so such a curve passes through every event.

### spacetime diagram caption: `diagrams/godel/systems/cartesian/tx.caption[1]`

Source: `_tools/derivations/null_rays.py` line 860, `CAPTIONS ("godel", "cartesian", "tx")`.

Before:

> The metric on this plane is $(-dt^2 + dx^2)/2\omega^2$, so the curves drawn are null and run at 45°.
> They are not null geodesics, though.
> $\Gamma^y{}_{tx} = -e^{-x}$ is not zero, so a light ray launched along one of these curves is turned out of the plane into $y$.
> The closed timelike curves of the Gödel universe circle each world line of the dust beyond a critical radius, through $y$ as well as $x$, so they cross this plane rather than lie in it.

After:

> The metric on this plane is $(-dt^2 + dx^2)/2\omega^2$, so the curves drawn are null and run at 45°.
> They are not null geodesics, though.
> Since $\Gamma^y{}_{tx} = -e^{-x}$ is not zero, a light ray launched along one of these curves is turned out of the plane into $y$.
> The closed timelike curves of the Gödel universe circle each world line of the dust beyond a critical radius, through $y$ as well as $x$, so they cross this plane rather than lie in it.

### spacetime diagram caption: `diagrams/godel/systems/cylindrical/beyond.caption[0]`

Source: `_tools/derivations/null_rays.py` line 951, `CAPTIONS ("godel", "cylindrical", "beyond")`.

Before:

> This is the cylinder of $t$ and $\phi$ at $r = 3r_c/2$ and $z = 0$, opened along $\phi = \pm\pi$ in the same way.
> Beyond $r_c$ the coefficient $g_{\phi\phi} = 2\sinh^2 r\,(1 - \sinh^2 r)/\omega^2$ is negative, and the cones have tipped over past the horizontal: the null curve moving to $+\phi$, $dt = -\tfrac{1}{2}(3 - \sqrt{2})\,d\phi$, goes down in $t$, while the one moving to $-\phi$, $dt = -\tfrac{1}{2}(17 - \sqrt{2})\,d\phi$, climbs steeply.
> Every horizontal line, run toward $+\phi$, points into the future cones, so the circle of constant $t$, $r$ and $z$ is a closed timelike curve.

After:

> This is the cylinder of $t$ and $\phi$ at $r = 3r_c/2$ and $z = 0$, opened along $\phi = \pm\pi$ in the same way.
> Beyond $r_c$ the coefficient $g_{\phi\phi} = 2\sinh^2 r\,(1 - \sinh^2 r)/\omega^2$ is negative, and the cones have tipped over past the horizontal: the null curve moving to $+\phi$, $dt = -\tfrac{1}{2}(3 - \sqrt{2})\,d\phi$, goes down in $t$, while the one moving to $-\phi$, $dt = -\tfrac{1}{2}(17 - \sqrt{2})\,d\phi$, climbs steeply.
> Every horizontal line, run toward $+\phi$, points into the future cones, so the circle of constant $t$, $r$, and $z$ is a closed timelike curve.

### caption of a figure in three dimensions: `diagrams/godel/projections/cylindrical/tipping.caption[0]`

Source: `_tools/derivations/projections.py` line 564, `CAPTIONS ("godel", "cylindrical", "tipping")`.

Before:

> This is the slice $z = 0$ of $t$, $r$ and $\phi$ about the axis $r = 0$, the world line of one particle of the dust, with $t$ up and $r$ as the radius, which puts the null directions $dt = \pm dr$ at 45°.
> The cones stand at $t = 0$ on the axis and around the circles $r = r_c/2$, $r_c$ and $3r_c/2$, with $\sinh r_c = 1$.
> On the axis they are upright, and farther out the cross term $g_{t\phi} = -2\sqrt{2}\sinh^2 r/\omega^2$ tips them over toward $+\phi$, counterclockwise seen from above.

After:

> This is the slice $z = 0$ of $t$, $r$, and $\phi$ about the axis $r = 0$, the world line of one particle of the dust, with $t$ up and $r$ as the radius, which puts the null directions $dt = \pm dr$ at 45°.
> The cones stand at $t = 0$ on the axis and around the circles $r = r_c/2$, $r_c$, and $3r_c/2$, with $\sinh r_c = 1$.
> On the axis they are upright, and farther out the cross term $g_{t\phi} = -2\sqrt{2}\sinh^2 r/\omega^2$ tips them over toward $+\phi$, counterclockwise seen from above.

## Schwarzschild's stellar interior metric

### history: `interior_schwarzschild/history[1]`

Source: `MFS/assets/data/metrics/interior_schwarzschild.json`, `history`, paragraph 2.

Before:

> Schwarzschild wrote his own exterior solution in an unusual coordinate.
> His original radial coordinate was defined so that the surface $r = r_s$ mapped to the coordinate origin, so the singularity appeared to him to sit there.
> The coordinate singularity we now associate with the Schwarzschild radius was, in his frame, sitting at the center of the geometry.
> Droste rederived the exterior solution in 1916 using the full coordinate range, making the structure of the solution considerably clearer [droste1917], though his paper was not published until 1917, almost certainly after Schwarzschild's death.

After:

> Schwarzschild wrote his own exterior solution in an unusual coordinate.
> His original radial coordinate was defined so that the surface $r = r_s$ mapped to the coordinate origin, so the singularity appeared to him to sit there.
> The coordinate singularity we now associate with the Schwarzschild radius was, in his frame, sitting at the center of the geometry.
> Johannes Droste rederived the exterior solution in 1916 using the full coordinate range, making the structure of the solution considerably clearer [droste1917], though his paper was not published until 1917, almost certainly after Schwarzschild's death.

### history: `interior_schwarzschild/history[2]`

Source: `MFS/assets/data/metrics/interior_schwarzschild.json`, `history`, paragraph 3.

Before:

> Schwarzschild also found that at an instant the space inside the star has the geometry of a spherical space, a cap of a three sphere whose radius of curvature he put at about 500 times the radius of the Sun, and he called it an interesting result of Einstein's theory that this geometry, until then "a mere possibility," must be real inside gravitating spheres [schwarzschild1916b].
> He noticed too that the mass a distant planet feels is smaller than the substance of the sphere, the more so the more compact it is, and read the difference as energy radiated away as a sphere of fixed mass contracts [schwarzschild1916b].
> Ludwig Flamm drew the picture later in 1916: a plane through the center of the star has the geometry of a spherical cap that meets the paraboloid of the exterior tangentially at the surface, the two together forming a kind of funnel [flamm1916].

After:

> Schwarzschild also found that at an instant the space inside the star has the geometry of a spherical space, a cap of a three sphere whose radius of curvature he put at about 500 times the radius of the Sun, and he called it an interesting result of Einstein's theory that this geometry, until then "a mere possibility," must be real inside gravitating spheres [schwarzschild1916b].
> He noticed too that the mass a distant planet feels is smaller than the substance of the sphere, the more so the more compact it is, and took the difference to be energy radiated away as a sphere of fixed mass contracts [schwarzschild1916b].
> Ludwig Flamm drew the picture later in 1916: a plane through the center of the star has the geometry of a spherical cap that meets the paraboloid of the exterior tangentially at the surface, the two together forming a kind of funnel [flamm1916].

### history: `interior_schwarzschild/history[3]`

Source: `MFS/assets/data/metrics/interior_schwarzschild.json`, `history`, paragraph 4.

Before:

> The interior solution carries a sharp physical result: the central pressure of the star diverges as the stellar radius $R$ approaches $\frac{9}{8}r_s$ from above [schwarzschild1916b].
> Schwarzschild drew the conclusion himself: "there is a limit to the concentration, above which a sphere of incompressible fluid can not exist" [schwarzschild1916b].
> Since the solution assumes an incompressible fluid, with a fixed equation of state, this sets an absolute compactness limit for this particular matter model.
> Buchdahl later generalized the bound in 1959 to arbitrary equations of state, showing that no static, spherically symmetric star can have $R < \frac{9}{4}GM/c^2$ [buchdahl1959].
> Subrahmanyan Chandrasekhar showed in 1964 that a star in general relativity becomes dynamically unstable before it can contract to that limiting radius, for any finite stiffness of its matter, and for the homogeneous sphere of this solution the $\frac{9}{8}$ limit is reached only as the ratio of specific heats becomes infinite [chandrasekhar1964].

After:

> The interior solution carries a sharp physical result: the central pressure of the star diverges as the stellar radius $R$ approaches $\frac{9}{8}r_s$ from above [schwarzschild1916b].
> Schwarzschild drew the conclusion himself: "there is a limit to the concentration, above which a sphere of incompressible fluid can not exist" [schwarzschild1916b].
> Since the solution assumes an incompressible fluid, with a fixed equation of state, this sets an absolute compactness limit for this particular matter model.
> Hans Buchdahl generalized the bound in 1959 to arbitrary equations of state, showing that no static, spherically symmetric star can have $R < \frac{9}{4}GM/c^2$ [buchdahl1959].
> Subrahmanyan Chandrasekhar showed in 1964 that a star in general relativity becomes dynamically unstable before it can contract to that limiting radius, for any finite stiffness of its matter, and for the homogeneous sphere of this solution the $\frac{9}{8}$ limit is reached only as the ratio of specific heats becomes infinite [chandrasekhar1964].

### conventions: `interior_schwarzschild/convention[0]`

Source: `MFS/assets/data/metrics/interior_schwarzschild.json`, `convention`, paragraph 1.

Before:

> The coordinate $r$ is the areal radius, $r \in [0, R]$, where $R$ is the stellar radius.
> The parameter $r_s = 2GM/c^2$ is the Schwarzschild radius of the star.
> The solution assumes a uniform (incompressible) perfect fluid and matches the exterior Schwarzschild metric at $r = R$.
> The interior is conformally flat, so its Weyl tensor vanishes identically, and the exterior it matches onto is not.

After:

> The coordinate $r$ is the areal radius, $r \in [0, R]$, where $R$ is the stellar radius.
> The parameter $r_s = 2GM/c^2$ is the Schwarzschild radius of the star.
> The star is a uniform (incompressible) perfect fluid, joined to the exterior Schwarzschild metric at $r = R$.
> The interior is conformally flat, so its Weyl tensor vanishes identically, while the exterior it joins is not conformally flat.

### conventions: `interior_schwarzschild/convention[1]`

Source: `MFS/assets/data/metrics/interior_schwarzschild.json`, `convention`, paragraph 2.

Before:

> The Ricci tensor is the standard contraction $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$, on which the Einstein tensor is that of a fluid of positive density, $G^t{}_t = -\dfrac{3r_s}{R^3} = -8\pi G\rho$.
> Units with $c = 1$ are used throughout, and components are taken in the chart $x^0 = ct$, which those units make $t$ itself, so $t$ is a length beside $r$, $r_s$ and $R$ and no component carries a factor of $c$.
> The geodesic equations are $\ddot{x}^\mu + \Gamma^\mu{}_{\nu\rho}\dot{x}^\nu\dot{x}^\rho = 0$ with $\Gamma$ the Christoffel symbols of that same chart, and the dots are derivatives of that chart with respect to an affine parameter $\lambda$, so $\dot{t}$ means $d(ct)/d\lambda$ and every term of every equation carries the dimensions of its left hand side.

After:

> With the Ricci tensor the standard contraction, $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$, the Einstein tensor is that of a fluid of positive density, $G^t{}_t = -\dfrac{3r_s}{R^3} = -8\pi G\rho$.
> We work in units with $c = 1$ throughout and take components in the chart $x^0 = ct$, which is then $t$ itself, so $t$ is a length beside $r$, $r_s$, and $R$, and no component carries a factor of $c$.
> The geodesic equations read $\ddot{x}^\mu + \Gamma^\mu{}_{\nu\rho}\dot{x}^\nu\dot{x}^\rho = 0$, with $\Gamma$ the Christoffel symbols of the same chart, and a dot is a derivative of that chart's coordinates with respect to an affine parameter $\lambda$, so $\dot{t}$ means $d(ct)/d\lambda$, and every term of every equation carries the dimensions of its left hand side.

### conformal diagram caption: `conformal/interior_schwarzschild/spherical.caption[1]`

Source: `_tools/derivations/conformal.py` line 2378, `CAPTIONS ("interior_schwarzschild", "spherical")`.

Before:

> Both sides give $g_{tt} = -(1 - r_s/R)$ at the surface, so $t$ is one coordinate throughout.
> The tortoise coordinate $r_* = \int\sqrt{g_{rr}/(-g_{tt})}\,dr$ runs from the centre through the surface, and $p, q = \arctan((t \mp r_*)/R)$ bring the spacetime into Minkowski's triangle, the causal structure of empty space, with the star a timelike tube from $i^-$ to $i^+$.

After:

> At the surface $g_{tt} = -(1 - r_s/R)$ on both sides, so $t$ is one coordinate throughout.
> The tortoise coordinate $r_* = \int\sqrt{g_{rr}/(-g_{tt})}\,dr$ runs from the centre through the surface, and $p, q = \arctan((t \mp r_*)/R)$ bring the spacetime into Minkowski's triangle, the causal structure of empty space, with the star a timelike tube from $i^-$ to $i^+$.

### embedding diagram caption: `embedding/interior_schwarzschild/star.caption[1]`

Source: `_tools/derivations/embedding.py` line 3359, `CAPTIONS ("interior_schwarzschild", "star")`.

Before:

> At $r = R$ both sides give $g_{rr} = 1/(1 - r_s/R)$, so the cap meets the paraboloid in one circle with one tangent plane, and the surface is smooth across the surface of the star.
> The vacuum paraboloid would run on down to a throat at $r_s$.
> The star, at $R = 1.5\,r_s$, ends it above there, and its circles shrink to a point at the centre instead.

After:

> At $r = R$, $g_{rr} = 1/(1 - r_s/R)$ on both sides, so the cap meets the paraboloid in one circle with one tangent plane, and the surface is smooth across the surface of the star.
> The vacuum paraboloid would run on down to a throat at $r_s$.
> The star, at $R = 1.5\,r_s$, ends it above there, and its circles shrink to a point at the centre instead.

## Kasner Anisotropic Vacuum Cosmology

### history: `kasner/history[2]`

Source: `MFS/assets/data/metrics/kasner.json`, `history`, paragraph 3.

Before:

> Those two conditions are a plane and a sphere.
> The plane $p_1 + p_2 + p_3 = 1$ cuts the unit sphere $p_1^2 + p_2^2 + p_3^2 = 1$ in a circle, and every Kasner universe is one point on it.
> Save for three isolated points where the geometry is secretly flat, one exponent comes out negative: two axes stretch, the third squeezes, and the squeeze is compulsory.
> Run it back to $t = 0$: the stretching directions crush to nothing while the contracting one runs away, so the singularity is met not as a point but as a cigar, unboundedly long and thin.

After:

> Those two conditions are a plane and a sphere.
> The plane $p_1 + p_2 + p_3 = 1$ cuts the unit sphere $p_1^2 + p_2^2 + p_3^2 = 1$ in a circle, and every Kasner universe is one point on it.
> Save for three isolated points where the spacetime is flat, one exponent comes out negative: two axes stretch, and the third squeezes.
> Run it back to $t = 0$: the stretching directions crush to nothing while the contracting one runs away, so the singularity is met as a cigar, unboundedly long and thin.

### history: `kasner/history[3]`

Source: `MFS/assets/data/metrics/kasner.json`, `history`, paragraph 4.

Before:

> In 1963 Evgeny Lifshitz and Isaak Khalatnikov took Kasner's exponents as the leading behavior of a region collapsing toward a singularity, counted the arbitrary functions such a solution can carry, and came up short.
> A singularity, they concluded, is what special initial data buys you, an artifact of Friedmann's symmetry rather than a verdict of the theory [lifshitz1963].
> In 1965 Roger Penrose took that premise apart by assuming no symmetry whatever: once a trapped surface forms, geodesics end, and generality is no protection [penrose1965].
> The Soviet position held anyway, until its own authors, and not Penrose, overturned it five years later.

After:

> In 1963 Evgeny Lifshitz and Isaak Khalatnikov took Kasner's exponents as the leading behavior of a region collapsing toward a singularity, counted the arbitrary functions such a solution can carry, and came up short.
> A singularity, they concluded, arises only from special initial data, an artifact of Friedmann's symmetry that general initial data avoid [lifshitz1963].
> Roger Penrose refuted that premise in 1965 without assuming any symmetry: once a trapped surface forms, geodesics end, however general the initial data [penrose1965].
> The Soviet position held anyway, until its own authors, rather than Penrose, overturned it five years later.

### history: `kasner/history[4]`

Source: `MFS/assets/data/metrics/kasner.json`, `history`, paragraph 5.

Before:

> Hunting their own error, Khalatnikov and Lifshitz announced in January 1970 that the general solution does carry a singularity [khalatnikov1970], and that summer Vladimir Belinskii joined them to set out the regime they had all missed.
> Near a singularity neighboring points stop speaking to one another, and each one alone sits in a Kasner solution for a while, hands the negative exponent to a different axis, then sits in another [bkl1970].
> The axis collapsing this epoch is expanding the next.
> Charles Misner had met the same oscillation in 1969 in the closed Bianchi type IX model [misner1969].

After:

> Hunting their own error, Khalatnikov and Lifshitz announced in January 1970 that the general solution does carry a singularity [khalatnikov1970], and that summer Vladimir Belinskii joined them to set out the regime they had all missed.
> Near a singularity neighboring points fall out of causal contact, and each one alone sits in a Kasner solution for a while, hands the negative exponent to a different axis, then sits in another [bkl1970].
> The axis collapsing this epoch is expanding the next.
> Charles Misner had met the same oscillation in 1969 in the closed Bianchi type IX model [misner1969].

### history: `kasner/history[5]`

Source: `MFS/assets/data/metrics/kasner.json`, `history`, paragraph 6.

Before:

> One number fixes the order in which the axes take their turn.
> Parametrize the exponents by $u$: each bounce sends $u$ to $u - 1$ until it falls below one, whereupon the era closes and the next opens at $1/u$, the map that manufactures continued fractions.
> An irrational start never repeats, so the descent is chaotic.
> So the simple example of section three became the local building block from which every modern account of a generic singularity is assembled: one Kasner solution per point of space, swapped over and over on the way down.
> Kasner thought the best thing in his paper was the section four metric, "apparently the simplest solution of Einstein's equations which has thus far been obtained." It is the section three metric that relativists remember.

After:

> One number fixes the order in which the axes take their turn.
> Parametrize the exponents by $u$: each bounce sends $u$ to $u - 1$ until it falls below one, whereupon the era closes and the next opens at $1/u$, the map that generates continued fractions.
> An irrational start never repeats, so the descent is chaotic.
> The simple example of section three thereby became the local building block from which every modern account of a generic singularity is assembled: one Kasner solution per point of space, swapped over and over on the way down.
> Kasner thought the best thing in his paper was the section four metric, "apparently the simplest solution of Einstein's equations which has thus far been obtained." Relativists remember the metric of section three.

### conventions: `kasner/convention[0]`

Source: `MFS/assets/data/metrics/kasner.json`, `convention`, paragraph 1.

Before:

> Coordinates are $(t, x, y, z)$ with $t$ having dimensions of time, and factors of $c$ are kept explicit.
> The expansion exponents obey $p_1 + p_2 + p_3 = 1$ and $p_1^2 + p_2^2 + p_3^2 = 1$, which leaves one free exponent; every component assumes both.
> The powers $t^{2p_i}$ are read with $t$ in a fixed unit of time, so they are dimensionless and the lengths sit in $x$, $y$ and $z$.
> Components are taken in the chart $x^0 = ct$, and the dots in the geodesic equations are derivatives of that same chart with respect to an affine parameter, so $\dot{t}$ means $d(ct)/d\lambda$ and every term of every equation carries the dimensions of its left hand side.
> The singularity is at $t = 0$ and the chart covers $t > 0$.

After:

> We use coordinates $(t, x, y, z)$, with $t$ carrying dimensions of time, and keep factors of $c$ explicit.
> The expansion exponents obey $p_1 + p_2 + p_3 = 1$ and $p_1^2 + p_2^2 + p_3^2 = 1$, which leaves one exponent free; every component assumes both conditions.
> We evaluate the powers $t^{2p_i}$ with $t$ in a fixed unit of time, so they are dimensionless, and the lengths sit in $x$, $y$, and $z$.
> We take components in the chart $x^0 = ct$, and a dot in the geodesic equations is a derivative of that chart's coordinates with respect to an affine parameter, so $\dot{t}$ means $d(ct)/d\lambda$, and every term of every equation carries the dimensions of its left hand side.
> The singularity is at $t = 0$, and the chart covers $t > 0$.

### embedding diagram caption: `embedding/kasner/ring.caption[0]`

Source: `_tools/derivations/embedding.py` line 3620, `CAPTIONS ("kasner", "ring")`.

Before:

> This is the plane $y = 0$ of Kasner's universe at four moments of $t$, each drawn as a surface in flat space so that every distance along it is the metric distance.
> At every moment the plane is flat, $t^{2p_1}dx^2 + t^{2p_3}dz^2$ being Euclid's plane with its axes scaled, so the drawing is a flat disc, and what the geometry does shows in a ring of particles at rest in the chart, which stay at rest because the metric has no $\Gamma^i{}_{tt}$.

After:

> This is the plane $y = 0$ of Kasner's universe at four moments of $t$, each drawn as a surface in flat space so that every distance along it is the metric distance.
> At every moment the plane is flat, $t^{2p_1}dx^2 + t^{2p_3}dz^2$ being Euclid's plane with its axes scaled, so the drawing is a flat disc, and the uneven expansion shows in a ring of particles at rest in the chart, which stay at rest because the metric has no $\Gamma^i{}_{tt}$.

### embedding diagram caption: `embedding/kasner/ring.caption[1]`

Source: `_tools/derivations/embedding.py` line 3625, `CAPTIONS ("kasner", "ring")`.

Before:

> The ring is the circle $x^2 + z^2 = \ell^2$ at $t = 1$, and at time $t$ the ellipse reaching $t^{p_1}\ell$ along $x$ and $t^{p_3}\ell$ along $z$.
> With $(p_1, p_2, p_3) = (-2/7, 3/7, 6/7)$ the direction $x$ contracts while $y$ and $z$ expand, so toward the singularity at $t = 0$ every sphere of particles is drawn out into a needle along $x$.
> Edward Kasner found the solution in 1921: its exponents sum to $1$, so volumes grow as $t$, and the vacuum asks that their squares sum to $1$ as well.

After:

> The ring is the circle $x^2 + z^2 = \ell^2$ at $t = 1$, and at time $t$ the ellipse reaching $t^{p_1}\ell$ along $x$ and $t^{p_3}\ell$ along $z$.
> With $(p_1, p_2, p_3) = (-2/7, 3/7, 6/7)$ the direction $x$ contracts while $y$ and $z$ expand, so toward the singularity at $t = 0$ every sphere of particles is drawn out into a needle along $x$.
> Edward Kasner found the solution in 1921: its exponents sum to $1$, so volumes grow as $t$, and the vacuum field equations require their squares to sum to $1$ as well.

### embedding diagram note: `embedding/kasner/ring.input`

Source: `_tools/derivations/embedding.py` line 2994, `in kasner() input=`.

Before:

> Exponents $(p_1, p_2, p_3) = (-2/7, 3/7, 6/7)$, a point on the Kasner circle, as the spacetime diagrams declare.

After:

> Exponents $(p_1, p_2, p_3) = (-2/7, 3/7, 6/7)$, a point on the Kasner circle, as in the spacetime diagrams.

### embedding diagram note: `embedding/kasner/ring.settings`

Source: `_tools/derivations/embedding.py` line 2992, `in kasner() settings=`.

Before:

> $(p_1, p_2, p_3) = (-2/7, 3/7, 6/7)$ and $t$ in the unit the powers are read in, with $\ell$ the ring's radius at $t = 1$, the unit of every length.

After:

> $(p_1, p_2, p_3) = (-2/7, 3/7, 6/7)$ and $t$ in the unit of time in which the powers are evaluated, with $\ell$ the ring's radius at $t = 1$, the unit of every length.

## Kerr Rotating Black Hole

### history: `kerr/history[0]`

Source: `MFS/assets/data/metrics/kerr.json`, `history`, paragraph 1.

Before:

> For nearly half a century after Schwarzschild, the rotating case resisted everyone.
> Schwarzschild's 1916 solution described a static, spherical mass, but real stars spin, and the search for the exterior field of a rotating body became one of the notorious unsolved problems of general relativity, attempted and abandoned by a generation of relativists.
> The breakthrough came in 1963 from Roy Kerr, a New Zealander then at the University of Texas, who found it not by brute force but by restricting attention to metrics of a special algebraic type and letting that structure guide him to the answer [kerr1963].
> The paper is two pages long.
> It is, by a common reckoning, the most important exact solution ever found to Einstein's equations.

After:

> For nearly half a century after Schwarzschild, the rotating case resisted everyone.
> Schwarzschild's 1916 solution described a static, spherical mass, but real stars spin, and the search for the exterior field of a rotating body became one of the notorious unsolved problems of general relativity, attempted and abandoned by a generation of relativists.
> The breakthrough came in 1963 from Roy Kerr, a New Zealander then at the University of Texas, who found it by restricting attention to metrics of a special algebraic type and letting that structure guide him to the answer, where brute force had failed [kerr1963].
> The paper is two pages long.
> Its solution is, by a common reckoning, the most important exact solution ever found to Einstein's equations.

### history: `kerr/history[2]`

Source: `MFS/assets/data/metrics/kerr.json`, `history`, paragraph 3.

Before:

> The geometry Kerr uncovered is far richer than Schwarzschild's.
> A rotating black hole drags spacetime around with it, so inertial frames near the hole are swept along in the direction of spin, and this frame dragging becomes inescapable inside a region called the ergosphere, where nothing can remain at rest relative to the distant stars.
> The solution has two horizons rather than one, and at its core the point singularity of Schwarzschild is stretched into a ring.
> Others made much of this structure clear only later: in 1967 Robert Boyer and Richard Lindquist gave the coordinates now universally used to write the metric and built its maximal analytic extension [boyer1967], and in 1968 Brandon Carter discovered a hidden constant of motion that renders the geodesic equations integrable, which is the reason orbits around a Kerr hole can be solved at all [carter1968].
> Carter found the price of the extension as well: in every case but the spherical one the extended geometry contains closed timelike lines, and the only geodesics that reach the ring singularity lie in its equatorial plane [carter1968].
> Carl Sagan ruled that extension out as a route between the stars in Contact, where Eda, one of the travellers, explains that the exact Kerr solution has an interior tunnel but that it is unstable, the slightest perturbation sealing it off into a singularity through which nothing can pass [sagan1985].

After:

> The geometry Kerr uncovered is far richer than Schwarzschild's.
> A rotating black hole drags spacetime around with it, so inertial frames near the hole are swept along in the direction of spin, and this frame dragging becomes inescapable inside a region called the ergosphere, where nothing can remain at rest relative to the distant stars.
> The solution has two horizons rather than one, and at its core the point singularity of Schwarzschild is stretched into a ring.
> Others made much of this structure clear only later: in 1967 Robert Boyer and Richard Lindquist gave the coordinates now universally used to write the metric and built its maximal analytic extension [boyer1967], and in 1968 Brandon Carter discovered a hidden constant of motion that renders the geodesic equations integrable, so that orbits around a Kerr hole can be solved exactly [carter1968].
> Carter found the price of the extension as well: in every case but the spherical one the extended geometry contains closed timelike lines, and the only geodesics that reach the ring singularity lie in its equatorial plane [carter1968].
> Carl Sagan ruled that extension out as a route between the stars in Contact, where Eda, one of the travellers, explains that the exact Kerr solution has an interior tunnel but that it is unstable, the slightest perturbation sealing it off into a singularity through which nothing can pass [sagan1985].

### history: `kerr/history[3]`

Source: `MFS/assets/data/metrics/kerr.json`, `history`, paragraph 4.

Before:

> The ergosphere is not merely a curiosity.
> In 1969 Roger Penrose showed that a particle entering it can split so that one fragment falls in carrying negative energy while the other escapes with more energy than the original possessed, energy mined directly from the hole's rotation [penrose1969].
> Roger Blandford and Roman Znajek showed in 1977 that magnetic field lines threading a spinning black hole can carry off its energy and angular momentum electromagnetically, and proposed the mechanism as the engine of active galactic nuclei [blandford1977].

After:

> The ergosphere stores energy that can be drawn out.
> In 1969 Roger Penrose showed that a particle entering it can split so that one fragment falls in carrying negative energy while the other escapes with more energy than the original possessed, energy mined directly from the hole's rotation [penrose1969].
> Roger Blandford and Roman Znajek showed in 1977 that magnetic field lines threading a spinning black hole can carry off its energy and angular momentum electromagnetically, and proposed the mechanism as the engine of active galactic nuclei [blandford1977].

### history: `kerr/history[4]`

Source: `MFS/assets/data/metrics/kerr.json`, `history`, paragraph 5.

Before:

> Saul Teukolsky found in 1972 that the gravitational and electromagnetic perturbations of a Kerr black hole obey wave equations that separate completely [teukolsky1972, teukolsky1973], and with William Press he used them in 1973 to test the hole's stability, finding no unstable mode for any spin, even close to the extreme [press1973].
> And the Kerr metric carries a significance beyond any single phenomenon: uniqueness theorems establish that an isolated, stationary black hole in our universe is characterized entirely by its mass and spin, which is to say that it is a Kerr black hole.
> Werner Israel proved in 1967 that the only static black hole is Schwarzschild's [israel1967], Brandon Carter showed in 1971 that an axisymmetric black hole has only two degrees of freedom, its mass and its angular momentum [carter1971], Stephen Hawking showed in 1972 that a stationary black hole that rotates must be axisymmetric [hawking1972], and David Robinson completed the argument in 1975 by proving the Kerr family the only one [robinson1975].
> Subrahmanyan Chandrasekhar closed his Nobel lecture of 1983 on this point: stationary, isolated black holes "are all, every single one of them, described exactly by the Kerr solution," which he called "the only instance we have of an exact description of a macroscopic object" [chandrasekhar1984].

After:

> Saul Teukolsky found in 1972 that the gravitational and electromagnetic perturbations of a Kerr black hole obey wave equations that separate completely [teukolsky1972, teukolsky1973], and with William Press he used them in 1973 to test the hole's stability, finding no unstable mode for any spin, even close to the extreme [press1973].
> Uniqueness theorems establish that an isolated, stationary black hole in our universe is characterized entirely by its mass and spin, which is to say that it is a Kerr black hole.
> Werner Israel proved in 1967 that the only static black hole is Schwarzschild's [israel1967], Brandon Carter showed in 1971 that an axisymmetric black hole has only two degrees of freedom, its mass and its angular momentum [carter1971], Stephen Hawking showed in 1972 that a stationary black hole that rotates must be axisymmetric [hawking1972], and David Robinson completed the argument in 1975 by proving the Kerr family the only one [robinson1975].
> Subrahmanyan Chandrasekhar closed his Nobel lecture of 1983 on this point: stationary, isolated black holes "are all, every single one of them, described exactly by the Kerr solution," which he called "the only instance we have of an exact description of a macroscopic object" [chandrasekhar1984].

### history: `kerr/history[5]`

Source: `MFS/assets/data/metrics/kerr.json`, `history`, paragraph 6.

Before:

> Kip Thorne and the visual effects team at Double Negative drew Gargantua, the black hole of Christopher Nolan's film Interstellar of 2014, by tracing bundles of light through the Kerr metric in Boyer-Lindquist coordinates, with its spin set at 0.6 of the maximum for the images although the story's time dilations need a hole spinning almost as fast as it can [nolan2014interstellar, james2015lensing].
> The first detection of gravitational waves, on 14 September 2015, recorded two black holes of about 36 and 29 solar masses merging into one of 62 with a spin about two thirds of the maximum, the fading end of the signal matching a hole relaxing to a final stationary Kerr configuration [abbott2016].
> The Event Horizon Telescope collaboration published the first image of a black hole's shadow in 2019, a ring 42 microarcseconds across at the core of the galaxy M87 around a mass of about 6.5 billion Suns, and found it consistent with the shadow of a Kerr black hole as general relativity predicts [eht2019].
> Every quiescent black hole astronomers observe, from stellar remnants to the giants at the centers of galaxies, is described by these two pages of 1963.

After:

> Kip Thorne and the visual effects team at Double Negative drew Gargantua, the black hole of Christopher Nolan's film Interstellar of 2014, by tracing bundles of light through the Kerr metric in Boyer-Lindquist coordinates, with its spin set at 0.6 of the maximum for the images although the story's time dilations need a hole spinning almost as fast as it can [nolan2014interstellar, james2015lensing].
> The first detection of gravitational waves, on 14 September 2015, recorded two black holes of about 36 and 29 solar masses merging into one of 62 with a spin about two thirds of the maximum, the fading end of the signal matching a hole relaxing to a final stationary Kerr configuration [abbott2016].
> The Event Horizon Telescope collaboration published the first image of a black hole's shadow in 2019, a ring 42 microarcseconds across at the core of the galaxy M87 around a mass of about 6.5 billion Suns, and found it consistent with the shadow of a Kerr black hole as general relativity predicts [eht2019].
> Every quiescent black hole astronomers observe, from stellar remnants to the giants at the centers of galaxies, is described by the metric of those two pages of 1963.

### conventions: `kerr/convention[0]`

Source: `MFS/assets/data/metrics/kerr.json`, `convention`, paragraph 1.

Before:

> Boyer-Lindquist coordinates are $(t, r, \theta, \phi)$ with $t$ having dimensions of time, and factors of $c$ and $G$ are kept explicit.
> The two parameters are the mass $M$ and the spin per unit mass $a = J/(Mc)$, which carries dimensions of length; the Schwarzschild radius of the same mass is $r_s = 2GM/c^2$, and every component is written with $2GM/c^2$ in place of it.
> The abbreviations the literature reads the solution through are $\Sigma = r^2 + a^2\cos^2\theta$ and $\Delta = r^2 - r_s r + a^2$, and both are written out in full wherever they appear, so that no component leans on a symbol left undefined.
> The horizons are the two roots of $\Delta = 0$, at $r_\pm = GM/c^2 \pm \sqrt{G^2M^2/c^4 - a^2}$, and the chart covers the exterior $r > r_+$.
> Components are taken in the chart $x^0 = ct$, so a component carrying an index $t$ is a component in that chart even though the index is written with the bare letter, and the dots in the geodesic equations are derivatives of that same chart with respect to an affine parameter, so $\dot{t}$ means $d(ct)/d\lambda$ and every term of every equation carries the dimensions of its left hand side.

After:

> We use the coordinates of Robert Boyer and Richard Lindquist, $(t, r, \theta, \phi)$, with $t$ carrying dimensions of time, and keep factors of $c$ and $G$ explicit.
> The two parameters are the mass $M$ and the spin per unit mass $a = J/(Mc)$, which carries dimensions of length; the Schwarzschild radius of the same mass is $r_s = 2GM/c^2$, and we write every component with $2GM/c^2$ in its place.
> The abbreviations common in work on the solution are $\Sigma = r^2 + a^2\cos^2\theta$ and $\Delta = r^2 - r_s r + a^2$, and we write both out in full wherever they appear, so that no component depends on a symbol left undefined.
> The horizons are the two roots of $\Delta = 0$, at $r_\pm = GM/c^2 \pm \sqrt{G^2M^2/c^4 - a^2}$, and the chart covers the exterior $r > r_+$.
> We take components in the chart $x^0 = ct$, so a component carrying an index $t$ is a component in that chart even though we write the index with the bare letter, and a dot in the geodesic equations is a derivative of that chart's coordinates with respect to an affine parameter, so $\dot{t}$ means $d(ct)/d\lambda$, and every term of every equation carries the dimensions of its left hand side.

### conventions: `kerr/convention[1]`

Source: `MFS/assets/data/metrics/kerr.json`, `convention`, paragraph 2.

Before:

> The Ricci tensor is the standard contraction $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$, though for Kerr the choice makes no difference at all, because both contractions of the Riemann tensor vanish.
> Kerr is a vacuum solution, so the Ricci tensor, the Ricci scalar and the Einstein tensor are identically zero.
> Because the Ricci terms the Weyl tensor subtracts from Riemann are built out of those same vanishing traces, $C^\mu{}_{\nu\rho\sigma} = R^\mu{}_{\nu\rho\sigma}$ here, and the Weyl tensor evaluated from its definition agrees with Riemann component by component.
> What carries the curvature is the Kretschmann scalar, which is finite at both horizons, so neither is a singularity, and which diverges only where $\Sigma = 0$, that is on the ring $r = 0$, $\theta = \pi/2$.

After:

> The Ricci tensor is the standard contraction, $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$, though for Kerr the choice makes no difference, because both contractions of the Riemann tensor vanish.
> Kerr is a vacuum solution, so the Ricci tensor, the Ricci scalar, and the Einstein tensor are identically zero.
> Because the Ricci terms subtracted from Riemann to form the Weyl tensor are built from those same vanishing traces, $C^\mu{}_{\nu\rho\sigma} = R^\mu{}_{\nu\rho\sigma}$ here, and the Weyl tensor evaluated from its definition agrees with Riemann component by component.
> The curvature shows in the Kretschmann scalar, which is finite at both horizons, so neither is a singularity, and which diverges only where $\Sigma = 0$, that is, on the ring $r = 0$, $\theta = \pi/2$.

### spacetime diagram caption: `diagrams/kerr/systems/boyer_lindquist/radial.caption[0]`

Source: `_tools/derivations/null_rays.py` line 742, `CAPTIONS ("kerr", "boyer_lindquist", "radial")`.

Before:

> This is the plane of $t$ and $r$ on the rotation axis, $\theta = 0$, drawn for $a = 0.9\,GM/c^2$.
> The curves drawn are null, and on the axis they are also null geodesics, the paths light takes.
> Off the axis a light ray launched along a curve of fixed $\theta$ and $\phi$ is turned out of the plane by $\Gamma^\theta{}_{tt}$, $\Gamma^\theta{}_{rr}$ and $\Gamma^\phi{}_{tr}$.
> On the axis $g^{rr} = \Delta/\Sigma$ vanishes at both roots of $\Delta$, $r_\pm = GM/c^2 \pm \sqrt{(GM/c^2)^2 - a^2}$, and the cones close at both; between them they point to smaller $r$, following the ingoing family.

After:

> This is the plane of $t$ and $r$ on the rotation axis, $\theta = 0$, drawn for $a = 0.9\,GM/c^2$.
> The curves drawn are null, and on the axis they are also null geodesics, the paths light takes.
> Off the axis a light ray launched along a curve of fixed $\theta$ and $\phi$ is turned out of the plane by $\Gamma^\theta{}_{tt}$, $\Gamma^\theta{}_{rr}$, and $\Gamma^\phi{}_{tr}$.
> On the axis $g^{rr} = \Delta/\Sigma$ vanishes at both roots of $\Delta$, $r_\pm = GM/c^2 \pm \sqrt{(GM/c^2)^2 - a^2}$, and the cones close at both; between them they point to smaller $r$, following the ingoing family.

### conformal diagram caption: `conformal/kerr/axis.caption[0]`

Source: `_tools/derivations/conformal.py` line 2253, `CAPTIONS ("kerr", "axis")`.

Before:

> This is the symmetry axis $\theta = 0$ of the maximally extended Kerr spacetime, the surface the rotations leave fixed and so totally geodesic.
> On it the metric is $-\frac{\Delta}{r^2 + a^2}c^2dt^2 + \frac{r^2 + a^2}{\Delta}dr^2$ with $\Delta = r^2 - 2GMr/c^2 + a^2$, and its two simple roots give it the tower of Reissner-Nordström.
> Carter extended the axis this way in 1966.

After:

> This is the symmetry axis $\theta = 0$ of the maximally extended Kerr spacetime, the surface the rotations leave fixed and so totally geodesic.
> On it the metric is $-\frac{\Delta}{r^2 + a^2}c^2dt^2 + \frac{r^2 + a^2}{\Delta}dr^2$ with $\Delta = r^2 - 2GMr/c^2 + a^2$, and its two simple roots give it the tower of Reissner-Nordström.
> Brandon Carter extended the axis this way in 1966.

### embedding diagram caption: `embedding/kerr/equator.caption[1]`

Source: `_tools/derivations/embedding.py` line 3418, `CAPTIONS ("kerr", "equator")`.

Before:

> On the equator the throat's circumference is $4\pi GM/c^2$ whatever the spin, since $r_+^2 + a^2 = 2GMr_+/c^2$ there.
> The dotted circle is the edge of the ergosphere, $r = 2GM/c^2$ on the equator, inside which nothing can stand still against the rotation.
> Nothing in the shape of the slice marks it: the ergosphere lies in how the slices are stacked, the rotation dragging each one round past the next.

After:

> On the equator the throat's circumference is $4\pi GM/c^2$ whatever the spin, since $r_+^2 + a^2 = 2GMr_+/c^2$ there.
> The dotted circle is the edge of the ergosphere, $r = 2GM/c^2$ on the equator, inside which nothing can stand still against the rotation.
> The ergosphere leaves the shape of the slice unmarked and lies in how the slices are stacked, the rotation dragging each one round past the next.

## Kerr-Newman Charged Rotating Black Hole

### history: `kerr_newman/history[0]`

Source: `MFS/assets/data/metrics/kerr_newman.json`, `history`, paragraph 1.

Before:

> By the middle of the 1960s two black hole solutions carried structure beyond bare mass.
> Kerr had just written down the rotating hole in 1963 [kerr1963], and the charged case without rotation had been known since the work of Reissner and Nordström in 1916 and 1918 [reissner1916, nordstrom1918].
> The obvious remaining case, a hole that both spins and carries electric charge, fell almost immediately, and by an unexpected route.
> Ezra Newman and Allen Janis noticed that Kerr's metric could be conjured out of Schwarzschild's by a peculiar trick: shifting the radial coordinate into the complex plane and taking a particular slice [newman1965janis].
> The maneuver had no rigorous justification, yet it worked.

After:

> By the middle of the 1960s two black hole solutions carried structure beyond bare mass.
> Kerr had just written down the rotating hole in 1963 [kerr1963], and the charged case without rotation had been known since the work of Hans Reissner and Gunnar Nordström in 1916 and 1918 [reissner1916, nordstrom1918].
> The obvious remaining case, a hole that both spins and carries electric charge, fell almost immediately, and by an unexpected route.
> Ezra Newman and Allen Janis noticed that Kerr's metric could be conjured out of Schwarzschild's by a peculiar trick: shifting the radial coordinate into the complex plane and taking a particular slice [newman1965janis].
> The maneuver had no rigorous justification, yet it worked.

### history: `kerr_newman/history[1]`

Source: `MFS/assets/data/metrics/kerr_newman.json`, `history`, paragraph 2.

Before:

> In 1965 Newman, together with five collaborators, applied the same complex shift to the charged Reissner-Nordström solution and read off the result, a metric describing a mass that is at once rotating and charged [newman1965kn].
> His collaborators were E. Couch, K. Chinnapared, A. Exton, A. Prakash, and R. Torrence, all, like Newman, of the University of Pittsburgh [newman1965kn].
> This is the Kerr-Newman solution.
> It unifies every black hole metric that came before it: set the charge to zero and it returns to Kerr, switch off the spin and it becomes Reissner-Nordström, do both and Schwarzschild reappears.

After:

> In 1965 Newman, together with five collaborators, applied the same complex shift to the charged Reissner-Nordström solution and obtained a metric describing a mass that is at once rotating and charged [newman1965kn].
> His collaborators were E. Couch, K. Chinnapared, A. Exton, A. Prakash, and R. Torrence, all, like Newman, of the University of Pittsburgh [newman1965kn].
> This is the Kerr-Newman solution.
> It contains every black hole metric that came before it: with the charge set to zero it returns to Kerr, with the spin switched off it becomes Reissner-Nordström, and with both it is Schwarzschild.

### history: `kerr_newman/history[4]`

Source: `MFS/assets/data/metrics/kerr_newman.json`, `history`, paragraph 5.

Before:

> Its importance is that it is the end of the line.
> Kerr-Newman is the most general stationary, asymptotically flat black hole in Einstein-Maxwell theory, and the uniqueness theorems make this precise: an isolated black hole at equilibrium is fixed completely by three numbers, namely its mass, its angular momentum, and its electric charge [carter1968].
> Werner Israel proved the static cases in 1967 and 1968 [israel1967, israel1968], and Brandon Carter and Stephen Hawking carried the argument to rotating holes in 1971 and 1972 [carter1971, hawking1972].
> David Robinson proved in 1975 that the Kerr family is the unique family of vacuum black holes with a horizon that is not degenerate [robinson1975], and Paweł Mazur settled the charged case in 1982, proving that the only exterior of a stationary, rotating electrovacuum black hole with such a horizon is Kerr-Newman [mazur1982].

After:

> Kerr-Newman is the most general stationary, asymptotically flat black hole in Einstein-Maxwell theory, and the uniqueness theorems make this precise: an isolated black hole at equilibrium is fixed completely by three numbers, namely its mass, its angular momentum, and its electric charge [carter1968].
> Werner Israel proved the static cases in 1967 and 1968 [israel1967, israel1968], and Brandon Carter and Stephen Hawking carried the argument to rotating holes in 1971 and 1972 [carter1971, hawking1972].
> David Robinson proved in 1975 that the Kerr family is the unique family of vacuum black holes with a horizon that is not degenerate [robinson1975], and Paweł Mazur settled the charged case in 1982, proving that the only exterior of a stationary, rotating electrovacuum black hole with such a horizon is Kerr-Newman [mazur1982].

### history: `kerr_newman/history[5]`

Source: `MFS/assets/data/metrics/kerr_newman.json`, `history`, paragraph 6.

Before:

> Whatever rich and complicated star collapsed to form it, the remnant remembers only these three.
> This is the content of the dictum that black holes have no hair [ruffini1971].
> Gary Gibbons showed in 1975 that a charged black hole also loses its charge spontaneously by emitting particles [gibbons1975].
> In the universe we actually observe, charge is quickly neutralized, so real black holes are described by Kerr; Kerr-Newman stands as the theoretical capstone of the family rather than its astrophysical workhorse.

After:

> Whatever rich and complicated star collapsed to form it, the remnant remembers only these three.
> This is the content of the dictum that black holes have no hair [ruffini1971].
> Gary Gibbons showed in 1975 that a charged black hole also loses its charge spontaneously by emitting particles [gibbons1975].
> In the universe we observe, charge is quickly neutralized, so real black holes are described by Kerr, and Kerr-Newman matters for the theory of the family more than for astronomy.

### conventions: `kerr_newman/convention[0]`

Source: `MFS/assets/data/metrics/kerr_newman.json`, `convention`, paragraph 1.

Before:

> Boyer-Lindquist coordinates are $(t, r, \theta, \phi)$ with $t$ having dimensions of time, and factors of $c$ and $G$ are kept explicit.
> The parameters are the mass $M$, the spin per unit mass $a = J/(Mc)$, which carries dimensions of length, and the charge radius $r_Q$, fixed by $r_Q^2 = GQ^2/(4\pi\epsilon_0 c^4)$ from the total electric charge $Q$; the Schwarzschild radius of the same mass is $r_s = 2GM/c^2$, and every component is written with $2GM/c^2$ in place of it.
> The abbreviations the literature reads the solution through are $\Sigma = r^2 + a^2\cos^2\theta$ and $\Delta = r^2 - r_s r + a^2 + r_Q^2$, and both are written out in full wherever they appear, so that no component leans on a symbol left undefined.

After:

> We use the coordinates of Robert Boyer and Richard Lindquist, $(t, r, \theta, \phi)$, with $t$ carrying dimensions of time, and keep factors of $c$ and $G$ explicit.
> The parameters are the mass $M$, the spin per unit mass $a = J/(Mc)$, which carries dimensions of length, and the charge radius $r_Q$, fixed by $r_Q^2 = GQ^2/(4\pi\epsilon_0 c^4)$ from the total electric charge $Q$; the Schwarzschild radius of the same mass is $r_s = 2GM/c^2$, and we write every component with $2GM/c^2$ in its place.
> The abbreviations common in work on the solution are $\Sigma = r^2 + a^2\cos^2\theta$ and $\Delta = r^2 - r_s r + a^2 + r_Q^2$, and we write both out in full wherever they appear, so that no component depends on a symbol left undefined.

### conventions: `kerr_newman/convention[1]`

Source: `MFS/assets/data/metrics/kerr_newman.json`, `convention`, paragraph 2.

Before:

> Setting $r_Q = 0$ returns every component to Kerr's, in the same chart and the same form.
> The horizons are the two roots of $\Delta = 0$, at $r_\pm = GM/c^2 \pm \sqrt{G^2M^2/c^4 - a^2 - r_Q^2}$, and the chart covers the exterior $r > r_+$.
> Components are taken in the chart $x^0 = ct$, so a component carrying an index $t$ is a component in that chart even though the index is written with the bare letter, and the dots in the geodesic equations are derivatives of that same chart with respect to an affine parameter, so $\dot{t}$ means $d(ct)/d\lambda$ and every term of every equation carries the dimensions of its left hand side.

After:

> Setting $r_Q = 0$ returns every component to Kerr's, in the same chart and the same form.
> The horizons are the two roots of $\Delta = 0$, at $r_\pm = GM/c^2 \pm \sqrt{G^2M^2/c^4 - a^2 - r_Q^2}$, and the chart covers the exterior $r > r_+$.
> We take components in the chart $x^0 = ct$, so a component carrying an index $t$ is a component in that chart even though we write the index with the bare letter, and a dot in the geodesic equations is a derivative of that chart's coordinates with respect to an affine parameter, so $\dot{t}$ means $d(ct)/d\lambda$, and every term of every equation carries the dimensions of its left hand side.

### conventions: `kerr_newman/convention[2]`

Source: `MFS/assets/data/metrics/kerr_newman.json`, `convention`, paragraph 3.

Before:

> The Ricci tensor is the standard contraction $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$.
> Kerr-Newman is an electrovacuum rather than a vacuum, so unlike Kerr its Ricci tensor does not vanish: it is the Maxwell stress of the charged rotating field, every component of it carries the factor $r_Q^2$, and it is traceless, $R = 0$, which is why the Einstein tensor equals the Ricci tensor slot for slot in all three index positions.
> Its scale is $r_Q^2/\Sigma^2$, the energy density the locally nonrotating observer measures being $c^4r_Q^2/(8\pi G\Sigma^2)$.
> The Weyl tensor is therefore not the Riemann tensor here, as it is for Kerr: the difference between the two is exactly $r_Q^2/\Sigma^2$ times the Maxwell bivectors.
> The Kretschmann scalar is finite at both horizons, so neither is a singularity, and it diverges only where $\Sigma = 0$, that is on the ring $r = 0$, $\theta = \pi/2$.

After:

> The Ricci tensor is the standard contraction, $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$.
> Kerr-Newman is an electrovacuum, so, unlike Kerr's, its Ricci tensor does not vanish: it is the Maxwell stress of the charged rotating field, every component of it carries the factor $r_Q^2$, and it is traceless, $R = 0$, which is why the Einstein tensor equals the Ricci tensor slot for slot in all three index positions.
> Its scale is $r_Q^2/\Sigma^2$, the energy density measured by the locally nonrotating observer being $c^4r_Q^2/(8\pi G\Sigma^2)$.
> The Weyl tensor therefore differs from the Riemann tensor here, where for Kerr the two agree: the difference between them is exactly $r_Q^2/\Sigma^2$ times the Maxwell bivectors.
> The Kretschmann scalar is finite at both horizons, so neither is a singularity, and it diverges only where $\Sigma = 0$, that is, on the ring $r = 0$, $\theta = \pi/2$.

## Krasnikov Tube

### history: `krasnikov/history[1]`

Source: `MFS/assets/data/metrics/krasnikov.json`, `history`, paragraph 2.

Before:

> Speed will not do it.
> Geometry might, because in general relativity a traveller partly controls how far he has to go.
> Working at the Pulkovo observatory in St Petersburg, Krasnikov killed half the hope: in a globally hyperbolic spacetime, under assumptions he makes precise, nothing you do on the way gets you to Deneb sooner than light does.
> The other half of the hope survives, because nothing forbids you from rebuilding the road behind you.

After:

> Speed will not do it.
> Geometry might, because in general relativity a traveller partly controls how far they have to go.
> Working at the Pulkovo observatory in St Petersburg, Krasnikov ruled out one route: in a globally hyperbolic spacetime, under assumptions he makes precise, nothing you do on the way gets you to Deneb sooner than light does.
> The other route stays open, because nothing forbids you from rebuilding the road behind you.

### history: `krasnikov/history[2]`

Source: `MFS/assets/data/metrics/krasnikov.json`, `history`, paragraph 3.

Before:

> So he builds it.
> The metric is altered in a narrow region along the ship's worldline, and inside that region spacetime stays flat while the light cones swing open, far enough that the homeward edge of the cone tips toward decreasing coordinate time.
> The traveller returns running backward through the clock on Earth, his own worldline future-directed throughout.
> Krasnikov did this in two dimensions, and his preprint had been circulating since November 1995, so Allen Everett and Thomas Roman got there ahead of the journal: in 1997 they built the four-dimensional version, named it the Krasnikov tube, and titled the paper a superluminal subway [everett1997].
> Their tube is Minkowski outside and flat inside with all the curvature in thin walls, and the ship lays it down as it goes at ordinary sublight speed.
> That control is the point: a pilot at the center of an Alcubierre warp bubble is causally cut off from the bubble wall and can neither raise it on demand nor steer it [alcubierre1994, everett1997].

After:

> So he builds it.
> The metric is altered in a narrow region along the ship's worldline, and inside that region spacetime stays flat while the light cones swing open, far enough that the homeward edge of the cone tips toward decreasing coordinate time.
> The traveller returns running backward through the clock on Earth, their own worldline directed to the future throughout.
> Krasnikov did this in two dimensions, and his preprint had been circulating since November 1995, so Allen Everett and Thomas Roman got there ahead of the journal: in 1997 they built the version in four dimensions, named it the Krasnikov tube, and titled the paper a superluminal subway [everett1997].
> Their tube is Minkowski outside and flat inside with all the curvature in thin walls, and the ship lays it down as it goes at ordinary sublight speed.
> Here the tube differs from a warp drive: a pilot at the center of an Alcubierre warp bubble is causally cut off from the bubble wall and can neither raise it on demand nor steer it [alcubierre1994, everett1997].

### history: `krasnikov/history[3]`

Source: `MFS/assets/data/metrics/krasnikov.json`, `history`, paragraph 4.

Before:

> Krasnikov gave the payoff as a calendar.
> An astronaut leaves Earth in 2000, flies toward Deneb at close to light speed, and arrives in 3600.
> Coming home down the tube he has just laid, he lands on Earth in 2002 [krasnikov1998].
> Krasnikov called machines of this sort space machines and saw that they are square roots of time machines.
> Run the trick once and causality survives, because the far end only ever moves away from you; run it twice along two routes and you have closed a timelike curve.
> Everett and Roman made that precise: one tube carries no closed timelike curves, and two non-overlapping tubes make a time machine [everett1997].

After:

> Krasnikov gave the payoff as a calendar.
> An astronaut leaves Earth in 2000, flies toward Deneb at close to light speed, and arrives in 3600.
> Coming home down the tube just laid, the astronaut lands on Earth in 2002 [krasnikov1998].
> Krasnikov called machines of this sort space machines and saw that they are square roots of time machines.
> Run the trick once and causality survives, because the far end only ever moves away from you; run it twice along two routes and you have closed a timelike curve.
> Everett and Roman made that precise: one tube carries no closed timelike curves, and two tubes that do not overlap make a time machine [everett1997].

### history: `krasnikov/history[4]`

Source: `MFS/assets/data/metrics/krasnikov.json`, `history`, paragraph 5.

Before:

> Then Everett and Roman computed the stress-energy of their tube and ran it against the quantum inequalities, which cap how much negative energy density a region may carry and for how long.
> The walls come out unphysically thin and the totals absurd: a tube one meter long and one meter wide asks for negative energy of order $10^{16}$ galactic masses, and a tube reaching the nearest star asks for $10^{32}$ of them [everett1997].
> The obstruction then hardened into theorems.
> Ken Olum pinned down what superluminal travel means in a curved spacetime and proved, granting the generic condition, that any spacetime permitting it violates the weak energy condition [olum1998].
> Matt Visser, Bruce Bassett and Stefano Liberati tied the tipping of light cones in Einstein gravity to violation of the null energy condition and named the pattern superluminal censorship [visser2000].

After:

> Then Everett and Roman computed the stress energy of their tube and ran it against the quantum inequalities, which cap how much negative energy density a region may carry and for how long.
> The walls come out unphysically thin and the totals absurd: a tube one meter long and one meter wide asks for negative energy of order $10^{16}$ galactic masses, and a tube reaching the nearest star asks for $10^{32}$ of them [everett1997].
> The obstruction then hardened into theorems.
> Ken Olum pinned down what superluminal travel means in a curved spacetime and proved, granting the generic condition, that any spacetime permitting it violates the weak energy condition [olum1998].
> Matt Visser, Bruce Bassett, and Stefano Liberati tied the tipping of light cones in Einstein gravity to violation of the null energy condition and named the pattern superluminal censorship [visser2000].

### history: `krasnikov/history[5]`

Source: `MFS/assets/data/metrics/krasnikov.json`, `history`, paragraph 6.

Before:

> Krasnikov did not concede.
> In 2003 he attacked the quantum inequality itself, which is derived assuming spacetime looks approximately Minkowski below the radius of curvature, the very assumption a shortcut builder should want to break.
> He showed by example that the inequality need not force large densities, that large densities need not sum to a large total, and that the total is meaningless in some of the situations where it gets invoked [krasnikov2003].
> The case against the tube now rests on the energy conditions, and they are what anyone wanting to build one would have to break.
> Michael Flynn joined the worlds of his Spiral Arm novels by Krasnikov tubes, into which ships slide to cross between the stars, and one of his peoples, the Confederation, names its tubes after rivers [flynn2008january, flynn2010jim, flynn2012lion].

After:

> Krasnikov did not concede.
> In 2003 he attacked the quantum inequality itself, which is derived assuming spacetime looks approximately Minkowski below the radius of curvature, the very assumption a shortcut builder should want to break.
> He showed by example that the inequality need not force large densities, that large densities need not sum to a large total, and that the total is meaningless in some of the situations where it gets invoked [krasnikov2003].
> The case against the tube now rests on the energy conditions, which anyone building one would have to break.
> Michael Flynn joined the worlds of his Spiral Arm novels by Krasnikov tubes, into which ships slide to cross between the stars, and one of his peoples, the Confederation, names its tubes after rivers [flynn2008january, flynn2010jim, flynn2012lion].

### conventions: `krasnikov/convention[0]`

Source: `MFS/assets/data/metrics/krasnikov.json`, `convention`, paragraph 1.

Before:

> Coordinates are $(t, x, r, \phi)$, cylindrical about the axis the tube is laid along: $t$ carries dimensions of time, $x$ is length measured along the axis, $r$ is distance from the axis and $\phi$ is the angle around it, and factors of $c$ and $G$ are kept explicit.
> The tube is taken in four dimensions, in the form Everett and Roman gave Krasnikov's construction rather than Krasnikov's own line element, and the choice is forced.
> Krasnikov wrote the tube in two dimensions, as $ds^2 = -\left(c\,dt - dx\right)\left(c\,dt + k\,dx\right)$, and in two dimensions the Einstein tensor vanishes identically for every metric there is.

After:

> We use coordinates $(t, x, r, \phi)$, cylindrical about the axis along which the tube is laid: $t$ carries dimensions of time, $x$ is length measured along the axis, $r$ is distance from the axis, and $\phi$ is the angle around it, and we keep factors of $c$ and $G$ explicit.
> We take the tube in four dimensions, as Allen Everett and Thomas Roman wrote Serguei Krasnikov's construction, and the choice is forced.
> Krasnikov's own line element is two dimensional, $ds^2 = -\left(c\,dt - dx\right)\left(c\,dt + k\,dx\right)$, and in two dimensions the Einstein tensor vanishes identically for every metric there is.

### conventions: `krasnikov/convention[1]`

Source: `MFS/assets/data/metrics/krasnikov.json`, `convention`, paragraph 2.

Before:

> In two dimensions the tube has a line element, a connection, one curvature scalar and a set of tipped light cones, and not one component of its source, so nothing in two dimensions says what the tube costs; since what the tube costs is the whole of the physics, the two dimensional form cannot serve.
> That form is the same product with the flat transverse plane $dr^2 + r^2d\phi^2$ appended and with $k$ allowed to depend on $r$ as well, and the dependence on $r$ is exactly what makes the thing a tube.
> In two dimensions there is no direction for $k$ to fall off in, so the modified region is a slab filling all of space at the values of $x$ it covers and it has no wall; in four dimensions $k$ returns to $1$ beyond some radius, the modified region has a finite cross section, and the wall where $k$ turns is where every curvature component lives.

After:

> In two dimensions the tube has a line element, a connection, one curvature scalar, and a set of tipped light cones, but not one component of its source, so the two dimensional form cannot say what the tube costs, and without the cost it cannot serve.
> The four dimensional form is the same product with the flat transverse plane $dr^2 + r^2d\phi^2$ appended and with $k$ allowed to depend on $r$ as well, and the dependence on $r$ makes the thing a tube.
> In two dimensions $k$ has no direction to fall off in, so the modified region is a slab filling all of space at the values of $x$ it covers, with no wall; in four dimensions $k$ returns to $1$ beyond some radius, the modified region has a finite cross section, and every curvature component lives in the wall where $k$ turns.

### conventions: `krasnikov/convention[2]`

Source: `MFS/assets/data/metrics/krasnikov.json`, `convention`, paragraph 3.

Before:

> Components are taken in the chart $x^0 = ct$, where a chart component of the metric carries no dimensions.
> Every $\partial_t$ is the derivative along that chart coordinate, $\partial_t = c^{-1}\partial/\partial t$, so each derivative carries one inverse length per order.
> The dots in the geodesic equations are velocities of that same chart, so $\dot{t}$ means $d(ct)/d\lambda$, and the equations are $\ddot{x}^\mu + \Gamma^\mu{}_{\nu\sigma}\dot{x}^\nu\dot{x}^\sigma = 0$ with $\Gamma$ the Christoffel symbols of that same chart.
> The signature is $(-,+,+,+)$, the Riemann tensor is $R^\mu{}_{\nu\sigma\tau} = \partial_\sigma\Gamma^\mu{}_{\nu\tau} - \partial_\tau\Gamma^\mu{}_{\nu\sigma} + \Gamma^\mu{}_{\sigma\lambda}\Gamma^\lambda{}_{\nu\tau} - \Gamma^\mu{}_{\tau\lambda}\Gamma^\lambda{}_{\nu\sigma}$, and the Ricci tensor is the standard contraction $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$.
> The choice of contraction matters here: the wall of the tube is not a vacuum, so the Ricci tensor, the Ricci scalar and the Einstein tensor all change sign with it, and the whole question the tube raises is the sign of an energy density.

After:

> We take components in the chart $x^0 = ct$, where a chart component of the metric carries no dimensions.
> Every $\partial_t$ is the derivative along that chart coordinate, $\partial_t = c^{-1}\partial/\partial t$, so each order of derivative carries one inverse length.
> The dots in the geodesic equations are velocities in the same chart, so $\dot{t}$ means $d(ct)/d\lambda$, and the equations read $\ddot{x}^\mu + \Gamma^\mu{}_{\nu\sigma}\dot{x}^\nu\dot{x}^\sigma = 0$, with $\Gamma$ the Christoffel symbols of that chart.
> We work in the signature $(-,+,+,+)$, with the Riemann tensor $R^\mu{}_{\nu\sigma\tau} = \partial_\sigma\Gamma^\mu{}_{\nu\tau} - \partial_\tau\Gamma^\mu{}_{\nu\sigma} + \Gamma^\mu{}_{\sigma\lambda}\Gamma^\lambda{}_{\nu\tau} - \Gamma^\mu{}_{\tau\lambda}\Gamma^\lambda{}_{\nu\sigma}$ and the Ricci tensor its standard contraction, $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$.
> The choice of contraction matters here: the wall of the tube is not a vacuum, so the Ricci tensor, the Ricci scalar, and the Einstein tensor all change sign with it, and the tube turns on the sign of an energy density.

### conventions: `krasnikov/convention[4]`

Source: `MFS/assets/data/metrics/krasnikov.json`, `convention`, paragraph 5.

Before:

> The first is the outbound ray and it never moves, whatever $k$ does.
> The second is the whole of the construction: at $k = 1$ it is the ordinary ray running back down the axis at the speed of light, as $k$ falls it swings, at $k = 0$ it is purely spatial, and for $k < 0$ it has $d(ct) < 0$ along $dx < 0$, which is a future directed ray coming home while the coordinate time decreases.
> Along it $d(ct) = -k\,dx$, so Krasnikov's $k = -1 + \delta$ brings the ship home over $\Delta x = -D$ in $\Delta(ct) = -\left(1 - \delta\right)D$, and a round trip that went out at nearly the speed of light costs $\delta D$ in all.
> It can be made as short as one likes and it never comes out negative, which is why one tube is not a time machine and two of them are, and it is why $\delta$ cannot be taken to zero: the metric degenerates there.

After:

> The outbound ray never moves, whatever $k$ does.
> The tube tips the returning ray: at $k = 1$ it is the ordinary ray running back down the axis at the speed of light, as $k$ falls it swings, at $k = 0$ it is purely spatial, and for $k < 0$ it has $d(ct) < 0$ along $dx < 0$, a future directed ray coming home while the coordinate time decreases.
> Along it $d(ct) = -k\,dx$, so Krasnikov's $k = -1 + \delta$ brings the ship home over $\Delta x = -D$ in $\Delta(ct) = -\left(1 - \delta\right)D$, and a round trip that went out at nearly the speed of light costs $\delta D$ in all.
> The round trip can be made as short as one likes but never negative, so one tube is no time machine, though two of them make one, and $\delta$ cannot be taken to zero, since the metric degenerates there.

### conventions: `krasnikov/convention[5]`

Source: `MFS/assets/data/metrics/krasnikov.json`, `convention`, paragraph 6.

Before:

> One contraction carries the physics of the tube.
> On the outbound ray the source is asked for nothing at all, $G_{\mu\nu}\ell_{\text{out}}^\mu\ell_{\text{out}}^\nu = 0$ identically and for every $k$.
> On the ray that comes home, the one the tube was built to tip, $G_{\mu\nu}\ell_{\text{back}}^\mu\ell_{\text{back}}^\nu = -\dfrac{\left(1 + k\right)^2}{2r}\,\partial_r\left(r\,\partial_r\ln\left(1 + k\right)\right)$, which is minus the flat Laplacian of $\ln\left(1 + k\right)$ in the transverse plane against a positive weight.
> Its radial integral $\int_0^\infty \partial_r\left(r\,\partial_r\ln\left(1 + k\right)\right)dr = \left[r\,\partial_r\ln\left(1 + k\right)\right]_0^\infty$ vanishes for any $k$ regular on the axis and constant outside the tube.
> So if the null energy condition held everywhere on a slice, an integrand of one sign would integrate to zero and therefore vanish, $r\,\partial_r\ln\left(1 + k\right)$ would be constant and hence zero, $k$ would not depend on $r$, and there would be no tube.
> The null energy condition fails in the wall of every Krasnikov tube of finite radius, for every shape function, and it fails on exactly the ray the construction exists to open.

After:

> We now contract the Einstein tensor with the two null directions.
> Along the outbound ray the contraction vanishes, $G_{\mu\nu}\ell_{\text{out}}^\mu\ell_{\text{out}}^\nu = 0$ identically and for every $k$.
> Along the ray that comes home, which the tube was built to tip, $G_{\mu\nu}\ell_{\text{back}}^\mu\ell_{\text{back}}^\nu = -\dfrac{\left(1 + k\right)^2}{2r}\,\partial_r\left(r\,\partial_r\ln\left(1 + k\right)\right)$, minus the flat Laplacian of $\ln\left(1 + k\right)$ in the transverse plane against a positive weight.
> Its radial integral $\int_0^\infty \partial_r\left(r\,\partial_r\ln\left(1 + k\right)\right)dr = \left[r\,\partial_r\ln\left(1 + k\right)\right]_0^\infty$ vanishes for any $k$ regular on the axis and constant outside the tube.
> If the null energy condition held everywhere on a slice, an integrand of one sign would integrate to zero and so vanish, $r\,\partial_r\ln\left(1 + k\right)$ would be constant and hence zero, $k$ would not depend on $r$, and there would be no tube.
> The null energy condition therefore fails in the wall of every Krasnikov tube of finite radius, for every shape function, along the returning ray that the construction opens.

### conventions: `krasnikov/convention[6]`

Source: `MFS/assets/data/metrics/krasnikov.json`, `convention`, paragraph 7.

Before:

> The energy density is that same statement read by an observer.
> Since $g_{tt} = -1$ in this chart, $u^\mu = (1,0,0,0)$ is a unit timelike vector for every $k > -1$, so the observer sitting at fixed $x$, $r$ and $\phi$ measures $\rho = T_{\mu\nu}u^\mu u^\nu = \dfrac{c^4}{8\pi G}G_{tt}$, and $G_{tt}$ is $-\dfrac{1}{r}\partial_r\left(r\,\partial_r\ln\left(1 + k\right)\right) - \dfrac{1}{4}\left(\partial_r\ln\left(1 + k\right)\right)^2$, which carries no derivative along $t$ or along $x$ at all.
> The transverse plane is flat, its area element is $r\,dr\,d\phi$, and the first term is a total divergence in it, so the energy per unit length of the tube is $\displaystyle\int \rho\,dA = -\frac{c^4}{16G}\int_0^\infty \left(\partial_r\ln\left(1 + k\right)\right)^2 r\,dr$, minus a sum of squares exactly as Alcubierre's is.
> It is negative for every tube that has a wall and no shape function can turn it positive.
> A thin wall makes $\partial_r k$ large and a deep tube makes $1 + k$ small, and the density grows as the square of their ratio, which is the arithmetic behind the galactic masses Everett and Roman arrived at.

After:

> An observer at rest in the chart sees the same failure as a negative energy density.
> Since $g_{tt} = -1$ in this chart, $u^\mu = (1,0,0,0)$ is a unit timelike vector for every $k > -1$, so the observer sitting at fixed $x$, $r$, and $\phi$ measures $\rho = T_{\mu\nu}u^\mu u^\nu = \dfrac{c^4}{8\pi G}G_{tt}$, and $G_{tt}$ is $-\dfrac{1}{r}\partial_r\left(r\,\partial_r\ln\left(1 + k\right)\right) - \dfrac{1}{4}\left(\partial_r\ln\left(1 + k\right)\right)^2$, which carries no derivative along $t$ or along $x$.
> The transverse plane is flat, its area element is $r\,dr\,d\phi$, and the first term is a total divergence in it, so the energy per unit length of the tube is $\displaystyle\int \rho\,dA = -\frac{c^4}{16G}\int_0^\infty \left(\partial_r\ln\left(1 + k\right)\right)^2 r\,dr$, minus a sum of squares, as Alcubierre's is.
> It is negative for every tube that has a wall, and no shape function can turn it positive.
> A thin wall makes $\partial_r k$ large and a deep tube makes $1 + k$ small, and the density grows as the square of their ratio; this is the reason Everett and Roman found masses of galactic size.

### conventions: `krasnikov/convention[7]`

Source: `MFS/assets/data/metrics/krasnikov.json`, `convention`, paragraph 8.

Before:

> This is not a vacuum, so the Weyl tensor is not the Riemann tensor but the Riemann tensor with its traces removed: $C_{t\phi t\phi}$ is nonzero where $R_{t\phi t\phi}$ vanishes, and twenty four slots of $C_{\mu\nu\sigma\tau}$ are like that.
> The Kretschmann scalar is a polynomial in $k$ and its first and second derivatives over $\left(1 + k\right)^6$, so it stays finite wherever the shape function is smooth and $k$ stays above $-1$: the tube carries no curvature singularity anywhere, and what stands in its way is the source it demands rather than a place where the geometry breaks.
> Alcubierre's bubble and Natario's are the tube's nearest relatives, and the difference is what the geometry does with the traveller.
> A warp bubble carries him: he is at rest in a flat patch, the patch is what moves, and he is cut off from the wall he would have to raise.

After:

> The tube is not a vacuum, and its Weyl tensor is the Riemann tensor with the traces removed, traces that here do not vanish: $C_{t\phi t\phi}$ is nonzero where $R_{t\phi t\phi}$ vanishes, and twenty four slots of $C_{\mu\nu\sigma\tau}$ behave the same way.
> The Kretschmann scalar is a polynomial in $k$ and its first and second derivatives over $\left(1 + k\right)^6$, so it stays finite wherever the shape function is smooth and $k$ stays above $-1$; the tube has no curvature singularity, and the one obstacle to it is the source it demands, since nowhere does the geometry break.
> Alcubierre's bubble and Natário's are the tube's nearest relatives, and they differ from it in what the geometry does with the traveller.
> A warp bubble carries its traveller: they sit at rest in a flat patch, the patch moves, and they are cut off from the wall they would have to raise.

### conventions: `krasnikov/convention[8]`

Source: `MFS/assets/data/metrics/krasnikov.json`, `convention`, paragraph 9.

Before:

> The tube carries nobody.
> Its $x$ is an ordinary coordinate and the ship crosses it under its own power at ordinary speed, laying the modification down behind itself as it goes, so every component is a component of a geometry along a path already travelled.
> What the tube alters is the light cone rather than the position, and the whole of the saving is on the way home.

After:

> The tube, by contrast, carries no one.
> Its $x$ is an ordinary coordinate, and the ship crosses it under its own power at ordinary speed, laying the modification down behind itself as it goes, so every component belongs to a geometry along a path already travelled.
> The tube alters the light cone and leaves the position alone, and the time is saved on the way home.

### spacetime diagram note: `diagrams/krasnikov/systems/cylindrical/tx.input`

Source: `_tools/derivations/null_rays.py` line 443, `DIAGRAMS input=`.

Before:

> A tube along $x$ from $0$ to $D = 4$, built by a ship that left $x = 0$ at $t = 0$ at the speed of light: $k = 1 - (2 - \delta)\,S(\tfrac{\rho_0^2 - r^2}{2\rho_0})\,S(t - x)\,S(x)\,S(D - x)$ with $\delta = 0.2$, $\rho_0 = 1$ and $S$ a step of width $0.15$ built from $\tanh$.

After:

> A tube along $x$ from $0$ to $D = 4$, built by a ship that left $x = 0$ at $t = 0$ at the speed of light: $k = 1 - (2 - \delta)\,S(\tfrac{\rho_0^2 - r^2}{2\rho_0})\,S(t - x)\,S(x)\,S(D - x)$ with $\delta = 0.2$, $\rho_0 = 1$, and $S$ a step of width $0.15$ built from $\tanh$.

### embedding diagram caption: `embedding/krasnikov/plane.caption[0]`

Source: `_tools/derivations/embedding.py` line 3658, `CAPTIONS ("krasnikov", "plane")`.

Before:

> This is the tilt of the light cone in the Krasnikov tube, $1 - k$, drawn as a height over the plane of the tube's axis at the moment $ct = 5\rho_0$: a height of $\rho_0$ for $1 - k = 1$, and no surface of the spacetime.
> Along the back edge of the light cone $c\,dt = -k\,dx$, so a light signal sent home over a length $L$ arrives $(1 - k)L/c$ sooner than in flat space: $1 - k$ is $0$ outside the tube, $1$ on the curve $k = 0$, where the signal arrives at the moment it left, and close to $1.8$ deep inside, where it arrives $0.8L/c$ before it left.

After:

> This is the tilt of the light cone in the Krasnikov tube, $1 - k$, drawn as a height over the plane of the tube's axis at the moment $ct = 5\rho_0$, a height of $\rho_0$ for $1 - k = 1$; the height stands for the tilt alone.
> Along the back edge of the light cone $c\,dt = -k\,dx$, so a light signal sent home over a length $L$ arrives $(1 - k)L/c$ sooner than in flat space: $1 - k$ is $0$ outside the tube, $1$ on the curve $k = 0$, where the signal arrives at the moment it left, and close to $1.8$ deep inside, where it arrives $0.8L/c$ before it left.

### embedding diagram note: `embedding/krasnikov/plane.input`

Source: `_tools/derivations/embedding.py` line 2728, `in krasnikov() input=`.

Before:

> A tube along $x$ from $0$ to $D = 4$, built by a ship that left $x = 0$ at $t = 0$ at the speed of light: $k = 1 - (2 - \delta)\,S(\tfrac{\rho_0^2 - r^2}{2\rho_0})\,S(ct - x)\,S(x)\,S(D - x)$ with $\delta = 0.2$, $\rho_0 = 1$ and $S$ a step of width $0.15$ built from $\tanh$, as the spacetime diagram declares.

After:

> A tube along $x$ from $0$ to $D = 4$, built by a ship that left $x = 0$ at $t = 0$ at the speed of light: $k = 1 - (2 - \delta)\,S(\tfrac{\rho_0^2 - r^2}{2\rho_0})\,S(ct - x)\,S(x)\,S(D - x)$ with $\delta = 0.2$, $\rho_0 = 1$, and $S$ a step of width $0.15$ built from $\tanh$, as in the spacetime diagram.

## Lentz Positive-Energy Warp Drive

### history: `lentz/history[0]`

Source: `MFS/assets/data/metrics/lentz.json`, `history`, paragraph 1.

Before:

> Every warp drive from 1994 onward arrived with the same bill attached: negative energy density, exotic matter, the sort of source particle physics has no classical example of [alcubierre1994].
> Erik Lentz, then at Göttingen, went looking for the line in the derivation where that bill actually gets written, and found a choice nobody had thought to question.
> A warp geometry is fixed by its shift vector, the field that says how space flows past the ship, and every construction to that point had tied the components of that vector together in one of two ways.

After:

> Every warp drive from 1994 onward arrived with the same bill attached: negative energy density, exotic matter, the sort of source particle physics has no classical example of [alcubierre1994].
> Erik Lentz, then at Göttingen, went looking for the step in the derivation that forces the negative energy, and found a choice nobody had thought to question.
> A warp geometry is fixed by its shift vector, the field that says how space flows past the ship, and every construction to that point had tied the components of that vector together in one of two ways.

### history: `lentz/history[1]`

Source: `MFS/assets/data/metrics/lentz.json`, `history`, paragraph 2.

Before:

> Alcubierre used a linear relation, which produces the celebrated toroid of negative energy density around the bubble; Natário used an expansionless elliptic one, which forces the energy density into a negative square of the extrinsic curvature [natario2002].
> Parabolic and hyperbolic relations were sitting there, unexplored.
> Lentz took the hyperbolic one, a wave equation on each spatial slice, and the sign flipped [lentz2021].

After:

> Alcubierre used a linear relation, which produces the celebrated toroid of negative energy density around the bubble; Natário used an expansionless elliptic one, which forces the energy density into a negative square of the extrinsic curvature [natario2002].
> Parabolic and hyperbolic relations were sitting there, unexplored.
> Lentz took the hyperbolic one, a wave equation on each spatial slice, and the energy density came out positive [lentz2021].

### history: `lentz/history[2]`

Source: `MFS/assets/data/metrics/lentz.json`, `history`, paragraph 3.

Before:

> Lentz built his soliton from cells of hyperbolic source.
> Each cell carries positive energy density on its own and keeps it positive in the company of neighbours of the same size and orientation, so bigger solutions get assembled from them like bricks.
> There is a central region of mild tidal forces where proper time keeps step with coordinate time far away, so a passenger ages at the rate the rest of the universe does and stays put with respect to the soliton.
> And the stress-energy Lentz's geometry demands is that of a conducting plasma together with classical electromagnetic fields, which is familiar matter.
> It is not cheap.
> For a soliton a hundred metres in radius with a metre of source thickness, running at the speed of light, Lentz put the mass equivalent at a few tenths of a solar mass, which he called immense, and which sits at the same order as the estimate Michael Pfenning and Larry Ford had made for an Alcubierre bubble of those dimensions [pfenning1997].

After:

> Lentz built his soliton from cells of hyperbolic source.
> Each cell carries positive energy density on its own and keeps it positive in the company of neighbours of the same size and orientation, so bigger solutions get assembled from them like bricks.
> There is a central region of mild tidal forces where proper time keeps step with coordinate time far away, so a passenger ages at the rate the rest of the universe does and stays put with respect to the soliton.
> And the stress energy of Lentz's geometry is that of a conducting plasma together with classical electromagnetic fields, which is familiar matter.
> It is not cheap.
> For a soliton a hundred metres in radius with a metre of source thickness, running at the speed of light, Lentz put the mass equivalent at a few tenths of a solar mass, which he called immense, and which sits at the same order as the estimate Michael Pfenning and Larry Ford had made for an Alcubierre bubble of those dimensions [pfenning1997].

### history: `lentz/history[3]`

Source: `MFS/assets/data/metrics/lentz.json`, `history`, paragraph 4.

Before:

> Two other groups proposed warp drives with positive energy that same year.
> Alexey Bobrick and Gianni Martire gave a general framework for the whole family [bobrick2021]; Shaun Fell and Lavinia Heisenberg found hidden geometric structure in the energy that yielded superluminal solitons with positive semidefinite energy, four orders of magnitude below a solar mass [fell2021].
> Then Jessica Santiago, Sebastian Schuster, and Matt Visser pointed out an elementary and awkward fact about all three constructions: the positive energy density that Lentz, Bobrick and Martire, and Fell and Heisenberg each found is what one family of observers measures, the comoving Eulerian ones.
> The weak energy condition asks every timelike observer, and for a stress-energy of this shape an observer moving quickly enough, while still below light speed, measures the density as negative.
> They went further: generic warp drives violate the null energy condition, which drags the weak, strong and dominant conditions down with it [santiago2022].

After:

> Two other groups proposed warp drives with positive energy that same year.
> Alexey Bobrick and Gianni Martire gave a general framework for the whole family [bobrick2021]; Shaun Fell and Lavinia Heisenberg found hidden geometric structure in the energy that yielded superluminal solitons with positive semidefinite energy, four orders of magnitude below a solar mass [fell2021].
> Then Jessica Santiago, Sebastian Schuster, and Matt Visser pointed out an elementary and awkward fact about all three constructions: the positive energy density that Lentz, Bobrick and Martire, and Fell and Heisenberg each found is the density measured by one family of observers, the comoving Eulerian ones.
> The weak energy condition applies to every timelike observer, and for a stress energy of this shape an observer moving quickly enough, while still below light speed, measures the density as negative.
> They went further: generic warp drives violate the null energy condition, and so violate the weak, strong, and dominant conditions as well [santiago2022].

### history: `lentz/history[4]`

Source: `MFS/assets/data/metrics/lentz.json`, `history`, paragraph 5.

Before:

> Santiago, Schuster, and Visser had refuted the claim of positive energy, but Lentz had already changed what the argument was about.
> Before Lentz the argument was whether familiar matter could source a warp geometry at all.
> After him it was about which particular condition has to give, and where.
> In 2024 Jared Fuchs, Christopher Helmerich, Alexey Bobrick, Luke Sellers, Brandon Melcher, and Gianni Martire gave one answer: they wrapped a shift distribution like Alcubierre's around a stable matter shell of positive ADM mass and got a warp drive at constant velocity satisfying all of the energy conditions [fuchs2024].
> Their drive is subluminal, so nobody outruns light with it.
> It is still a warp bubble that ordinary matter pays for, which is what Lentz set out to find with his hyperbolic soliton.

After:

> Santiago, Schuster, and Visser had refuted the claim of positive energy, but Lentz had already changed what the argument was about.
> Before Lentz the argument was whether familiar matter could source a warp geometry.
> After him it was about which energy condition has to fail, and where.
> In 2024 Jared Fuchs, Christopher Helmerich, Alexey Bobrick, Luke Sellers, Brandon Melcher, and Gianni Martire gave one answer: they wrapped a shift distribution like Alcubierre's around a stable matter shell of positive ADM mass and got a warp drive at constant velocity satisfying all of the energy conditions [fuchs2024].
> Their drive is subluminal, so nobody outruns light with it.
> It is still a warp bubble sourced by ordinary matter, the thing Lentz set out to find with his hyperbolic soliton.

### conventions: `lentz/convention[0]`

Source: `MFS/assets/data/metrics/lentz.json`, `convention`, paragraph 1.

Before:

> Coordinates are $(t, x, y, z)$ with $t$ carrying dimensions of time and $x$, $y$ and $z$ lengths, and factors of $c$ are kept explicit.
> The line element is Lentz's own, the 3+1 form with unit lapse and flat slices, $ds^2 = -\left(N^2 - N^iN_i\right)c^2dt^2 - 2N_i\,dx^i\,c\,dt + \delta_{ij}dx^idx^j$ with $N = 1$, and with his hyperbolic class built in: the shift is the gradient of a potential, $N_i = \partial_i\phi$.
> The potential $\phi$ is a length and is left a free function of all four coordinates, so every component holds for any soliton of the class, in any state of motion.
> Lentz's own solitons move rigidly along $z$ at $c\,v_s$, with $\phi$ a function of $x$, $y$ and $z - z_s(t)$, so that $\partial_t\phi = -v_s\,\partial_z\phi$ in this chart; he ties the three second derivatives together with a wave equation on each slice, $\partial_x^2\phi + \partial_y^2\phi - \dfrac{2}{v_h^2}\partial_z^2\phi = \rho_h$, whose source $\rho_h$ he lays out in rhomboid cells and whose wave fronts cross a slice at $v_h/\sqrt{2}$, and he sets $v_s$ equal to the shift at the centre so that the passengers there ride with the soliton.
> None of that is assumed: it is a choice of $\phi$, and every component is an identity without it.

After:

> We use coordinates $(t, x, y, z)$, with $t$ carrying dimensions of time and $x$, $y$, and $z$ lengths, and keep every factor of $c$ explicit.
> The line element is Erik Lentz's own, the 3+1 form with unit lapse and flat slices, $ds^2 = -\left(N^2 - N^iN_i\right)c^2dt^2 - 2N_i\,dx^i\,c\,dt + \delta_{ij}dx^idx^j$ with $N = 1$, with his hyperbolic class built in: the shift is the gradient of a potential, $N_i = \partial_i\phi$.
> The potential $\phi$ is a length, and we leave it a free function of all four coordinates, so every component holds for any soliton of the class, in any state of motion.
> Lentz's own solitons move rigidly along $z$ at $c\,v_s$, with $\phi$ a function of $x$, $y$, and $z - z_s(t)$, so that $\partial_t\phi = -v_s\,\partial_z\phi$ in this chart; he ties the three second derivatives together with a wave equation on each slice, $\partial_x^2\phi + \partial_y^2\phi - \dfrac{2}{v_h^2}\partial_z^2\phi = \rho_h$, whose source $\rho_h$ he lays out in rhomboid cells and whose wave fronts cross a slice at $v_h/\sqrt{2}$, and he sets $v_s$ equal to the shift at the centre so that the passengers there ride with the soliton.
> We assume none of this: it is a choice of $\phi$, and every component is an identity without it.

### conventions: `lentz/convention[1]`

Source: `MFS/assets/data/metrics/lentz.json`, `convention`, paragraph 2.

Before:

> Components are taken in the chart $x^0 = ct$, where a chart component of the metric carries no dimensions.
> $\partial_t$ is the derivative along that chart coordinate, $\partial_t = c^{-1}\partial/\partial t$, so each derivative carries one inverse length per order, and the dots in the geodesic equations are velocities of the same chart, so $\dot{t}$ means $d(ct)/d\lambda$.
> The signature is $(-,+,+,+)$ and the Ricci tensor is the standard contraction $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$, which the Ricci scalar and the Einstein tensor follow, and that contraction fixes the sign of every energy density.

After:

> We take components in the chart $x^0 = ct$, where a chart component of the metric carries no dimensions.
> The derivative $\partial_t$ runs along that chart coordinate, $\partial_t = c^{-1}\partial/\partial t$, so each order of derivative carries one inverse length, and the dots in the geodesic equations are velocities in the same chart, so $\dot{t}$ means $d(ct)/d\lambda$.
> We work in the signature $(-,+,+,+)$, and the Ricci tensor is the standard contraction, $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$, from which the Ricci scalar and the Einstein tensor follow; that contraction fixes the sign of every energy density.

### conventions: `lentz/convention[2]`

Source: `MFS/assets/data/metrics/lentz.json`, `convention`, paragraph 3.

Before:

> The observers who ride the slices, with $n_\mu = (-1, 0, 0, 0)$, measure the energy density $\dfrac{c^4}{8\pi G}G^{tt}$, and for any potential $G^{tt}$ is the sum of the three principal minors of the Hessian of $\phi$, $G^{tt} = \partial_x^2\phi\,\partial_y^2\phi - \left(\partial_x\partial_y\phi\right)^2 + \partial_x^2\phi\,\partial_z^2\phi - \left(\partial_x\partial_z\phi\right)^2 + \partial_y^2\phi\,\partial_z^2\phi - \left(\partial_y\partial_z\phi\right)^2$, which is Lentz's $\tfrac{1}{2}\left(K^2 - K^i{}_jK^j{}_i\right)$ with the extrinsic curvature $K_{ij} = -\partial_i\partial_j\phi$.
> Unlike Alcubierre's and Natário's, this density has no fixed sign, and Lentz's hyperbolic relation and his rhomboids are how he makes it positive everywhere his soliton has any.
> The momentum the same observers measure, $\dfrac{c^4}{8\pi G}G^t{}_i$, vanishes identically for every potential, because the flow is a gradient, and that is Lentz's $J_i = 0$.
> What Santiago, Schuster and Visser add is visible in the rest of $G^{\mu\nu}$: the density another observer measures mixes in the stresses, and an observer moving fast enough relative to the slices finds it negative.
> The gradient flow also leaves the Kretschmann scalar a sum of squares, four times the squares of the cofactors of the Hessian and four times the squares of the tidal field of the riding observers, since the part of the curvature that would mix the two vanishes for any potential.

After:

> The observers who ride the slices, with $n_\mu = (-1, 0, 0, 0)$, measure the energy density $\dfrac{c^4}{8\pi G}G^{tt}$, and for any potential $G^{tt}$ is the sum of the three principal minors of the Hessian of $\phi$, $G^{tt} = \partial_x^2\phi\,\partial_y^2\phi - \left(\partial_x\partial_y\phi\right)^2 + \partial_x^2\phi\,\partial_z^2\phi - \left(\partial_x\partial_z\phi\right)^2 + \partial_y^2\phi\,\partial_z^2\phi - \left(\partial_y\partial_z\phi\right)^2$, which is Lentz's $\tfrac{1}{2}\left(K^2 - K^i{}_jK^j{}_i\right)$ with the extrinsic curvature $K_{ij} = -\partial_i\partial_j\phi$.
> Unlike Alcubierre's and Natário's, this density has no fixed sign, and Lentz makes it positive, wherever his soliton has any, with his hyperbolic relation and his rhomboids.
> The momentum measured by the same observers, $\dfrac{c^4}{8\pi G}G^t{}_i$, vanishes identically for every potential, because the flow is a gradient; this is Lentz's $J_i = 0$.
> Jessica Santiago, Sebastian Schuster, and Matt Visser turned to the rest of $G^{\mu\nu}$ and showed that the density another observer measures mixes in the stresses, and that an observer moving fast enough relative to the slices finds it negative.
> The gradient flow also leaves the Kretschmann scalar a sum of squares, four times the squares of the cofactors of the Hessian and four times the squares of the tidal field of the riding observers, since the part of the curvature that would mix the two vanishes for any potential.

### embedding diagram caption: `embedding/lentz/plane.caption[0]`

Source: `_tools/derivations/embedding.py` line 3700, `CAPTIONS ("lentz", "plane")`.

Before:

> This is the plane $y = 0$ of the path of Lentz's soliton at one moment, drawn as a surface in flat space so that every distance along it is the metric distance.
> Every slice of constant $t$ is flat, $dx^2 + dy^2 + dz^2$, for every potential $\phi$, because flat slices are one of the three things Erik Lentz fixed in 2021 to define his class, with a unit lapse and a shift that is the gradient of $\phi$.

After:

> This is the plane $y = 0$ of the path of Lentz's soliton at one moment, drawn as a surface in flat space so that every distance along it is the metric distance.
> Every slice of constant $t$ is flat, $dx^2 + dy^2 + dz^2$, for every potential $\phi$, because flat slices are one of the three conditions Erik Lentz imposed in 2021 to define his class, with a unit lapse and a shift that is the gradient of $\phi$.

### embedding diagram caption: `embedding/lentz/plane.caption[1]`

Source: `_tools/derivations/embedding.py` line 3704, `CAPTIONS ("lentz", "plane")`.

Before:

> His soliton exists only as a numerical integral over rhomboid cells of source, so no potential is drawn in its place and the plane carries only its path.
> On a flat slice the energy density the riding observers measure is $\sigma_2(\partial_i\partial_j\phi)\,c^4/8\pi G$, the sum of the principal minors of the Hessian of $\phi$, which Lentz arranged to be positive; Jessica Santiago, Sebastian Schuster and Matt Visser showed in 2022 that an observer moving fast enough through the slices measures it negative.

After:

> His soliton exists only as a numerical integral over rhomboid cells of source, and the plane carries its path alone.
> On a flat slice the energy density measured by the riding observers is $\sigma_2(\partial_i\partial_j\phi)\,c^4/8\pi G$, the sum of the principal minors of the Hessian of $\phi$, which Lentz arranged to be positive; Jessica Santiago, Sebastian Schuster, and Matt Visser showed in 2022 that an observer moving fast enough through the slices measures it negative.

## Malament-Hogarth Super-Turing Spacetimes

### history: `malament_hogarth/history[0]`

Source: `MFS/assets/data/metrics/malament_hogarth.json`, `history`, paragraph 1.

Before:

> Malament-Hogarth is not a metric.
> It is a property a spacetime may or may not possess, and the spacetimes that possess it are strange enough to deserve a name.
> A spacetime is Malament-Hogarth if it contains a worldline of infinite proper duration whose entire history lies within the causal past of some single event.
> Picture a computer that runs forever along that worldline, and an observer who travels to the special event in finite time of their own: every step the computer will ever take has, by the time the observer arrives, already had the chance to send its result ahead.
> The observer can therefore learn the outcome of an unending computation in a finite span of their own life.

After:

> Malament-Hogarth names a property of spacetimes rather than a metric, and the spacetimes that possess it are strange enough to deserve a name.
> A spacetime is Malament-Hogarth if it contains a worldline of infinite proper duration whose entire history lies within the causal past of some single event.
> Picture a computer that runs forever along that worldline, and an observer who travels to the special event in finite time of their own: every step the computer will ever take has, by the time the observer arrives, already had the chance to send its result ahead.
> The observer can therefore learn the outcome of an unending computation in a finite span of their own life.

### history: `malament_hogarth/history[1]`

Source: `MFS/assets/data/metrics/malament_hogarth.json`, `history`, paragraph 2.

Before:

> The idea grew out of work by Itamar Pitowsky in 1990 [pitowsky1990] and was put on firm footing by Mark Hogarth, who in 1992 proved that ordinary spacetimes cannot do this but that some can, anti-de Sitter space among them [hogarth1992].
> Pitowsky's version was set in Minkowski spacetime: a mathematician circles ever faster, so that his own clock runs only a finite time while his graduate students check Fermat's conjecture one case after another forever [pitowsky1990, earman1993].
> John Earman and John Norton showed that the mathematician's acceleration must grow without bound, so that any real traveler would be crushed [earman1993].
> Hogarth told his version as the story of the immortal computer HAL testing every even number against the Goldbach conjecture while the mortal Dave follows a worldline one hour long [hogarth1992].

After:

> The idea grew out of work by Itamar Pitowsky in 1990 [pitowsky1990], and Mark Hogarth put it on firm footing, proving in 1992 that ordinary spacetimes cannot do this but that some can, anti-de Sitter space among them [hogarth1992].
> Pitowsky's version was set in Minkowski spacetime: a mathematician circles ever faster, so that their own clock runs only a finite time while their graduate students check Fermat's conjecture one case after another forever [pitowsky1990, earman1993].
> John Earman and John Norton showed that the mathematician's acceleration must grow without bound, so that any real traveler would be crushed [earman1993].
> Hogarth told his version as the story of the immortal computer HAL testing every even number against the Goldbach conjecture while the mortal Dave follows a worldline one hour long [hogarth1992].

### history: `malament_hogarth/history[2]`

Source: `MFS/assets/data/metrics/malament_hogarth.json`, `history`, paragraph 3.

Before:

> The name was fixed a year later by Earman and Norton, who analyzed these supertask geometries and credited the construction to both Hogarth and David Malament [earman1993].
> Malament supplied the quick proof that such a spacetime cannot be globally hyperbolic, so that no Cauchy surface determines its whole history [earman1993].
> The property is a feature of global causal structure, not of any particular line element: Malament-Hogarth is a category of spacetimes, and its members include some of the best known exact solutions.
> Earman and Norton noted that in the Gödel universe every event is a Malament-Hogarth event, and that the Reissner-Nordström spacetime of the charged black hole has the property too [earman1993].

After:

> Earman and Norton fixed the name a year later, analyzing these supertask geometries and crediting the construction to both Hogarth and David Malament [earman1993].
> Malament supplied the quick proof that such a spacetime cannot be globally hyperbolic, so that no Cauchy surface determines its whole history [earman1993].
> The property belongs to the global causal structure, whatever the line element: Malament-Hogarth is a category of spacetimes, and its members include some of the best known exact solutions.
> Earman and Norton noted that in the Gödel universe every event is a Malament-Hogarth event, and that the Reissner-Nordström spacetime of the charged black hole has the property too [earman1993].

### history: `malament_hogarth/history[3]`

Source: `MFS/assets/data/metrics/malament_hogarth.json`, `history`, paragraph 4.

Before:

> The category is more than a curiosity because it bears on the limits of computation.
> Gábor Etesi and István Németi showed in 2002 that an observer in a Malament-Hogarth spacetime could, in principle, settle questions no Turing machine can decide, such as the halting problem or the consistency of arithmetic, by reading off the result of an infinite computation [etesi2002].
> They also noted that the requisite causal structure is not confined to exotic toy models: the interior of a rotating or charged black hole, past its inner horizon, is Malament-Hogarth.
> The obstacles to ever building such a device are severe, since the computer must survive an eternity and the observer must cross a horizon where the infalling history of the universe arrives infinitely blueshifted [earman1993], but the conclusion stands as a matter of principle: whether the Church-Turing thesis holds is not a fact of logic alone but depends on the shape of spacetime.

After:

> The category matters because it bears on the limits of computation.
> Gábor Etesi and István Németi showed in 2002 that an observer in a Malament-Hogarth spacetime could, in principle, settle questions no Turing machine can decide, such as the halting problem or the consistency of arithmetic, by receiving the result of an infinite computation [etesi2002].
> They also noted that the requisite causal structure is not confined to exotic toy models: the interior of a rotating or charged black hole, past its inner horizon, is Malament-Hogarth.
> The obstacles to ever building such a device are severe, since the computer must survive an eternity and the observer must cross a horizon where the infalling history of the universe arrives infinitely blueshifted [earman1993], but the conclusion stands as a matter of principle: whether the Church-Turing thesis holds depends on the shape of spacetime as well as on logic.

### conventions: `malament_hogarth/convention[0]`

Source: `MFS/assets/data/metrics/malament_hogarth.json`, `convention`, paragraph 1.

Before:

> Coordinates are $(t, x, y, z)$, the inertial coordinates of Minkowski space, with $t$ carrying dimensions of time and $x$, $y$ and $z$ lengths, and factors of $c$ are kept explicit.
> The chart is the toy spacetime Earman and Norton used to show that the property can be had at all, the construction they credit to Malament and to Hogarth: take Minkowski space, remove a single event, placed at the origin, and multiply the metric by $\Omega^2$.
> The conformal factor $\Omega$ is a dimensionless function, equal to $1$ outside a compact region $C$ around the removed event and growing without bound as that event is approached, so the chart covers every event but the origin and is Minkowski space wherever $\Omega = 1$.

After:

> We use coordinates $(t, x, y, z)$, the inertial coordinates of Minkowski space, with $t$ carrying dimensions of time and $x$, $y$, and $z$ lengths, and keep every factor of $c$ explicit.
> The chart is the toy spacetime John Earman and John Norton used to show that the property can be had, a construction they credit to David Malament and to Mark Hogarth: it starts from Minkowski space, removes a single event, placed at the origin, and multiplies the metric by $\Omega^2$.
> The conformal factor $\Omega$ is a dimensionless function, equal to $1$ outside a compact region $C$ around the removed event and growing without bound as that event is approached, so the chart covers every event but the origin and is Minkowski space wherever $\Omega = 1$.

### conventions: `malament_hogarth/convention[1]`

Source: `MFS/assets/data/metrics/malament_hogarth.json`, `convention`, paragraph 2.

Before:

> Components are taken in the chart $x^0 = ct$, where a chart component of the metric carries no dimensions.
> $\partial_t$ is the derivative along that chart coordinate, $\partial_t = c^{-1}\partial/\partial t$, so each derivative carries one inverse length per order, and the dots in the geodesic equations are velocities of the same chart, so $\dot{t}$ means $d(ct)/d\lambda$.
> The signature is $(-,+,+,+)$ and the Ricci tensor is the standard contraction $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$, which the Ricci scalar and the Einstein tensor follow.

After:

> We take components in the chart $x^0 = ct$, where a chart component of the metric carries no dimensions.
> The derivative $\partial_t$ runs along that chart coordinate, $\partial_t = c^{-1}\partial/\partial t$, so each order of derivative carries one inverse length, and the dots in the geodesic equations are velocities in the same chart, so $\dot{t}$ means $d(ct)/d\lambda$.
> We work in the signature $(-,+,+,+)$, and the Ricci tensor is the standard contraction, $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$, from which the Ricci scalar and the Einstein tensor follow.

### conventions: `malament_hogarth/convention[2]`

Source: `MFS/assets/data/metrics/malament_hogarth.json`, `convention`, paragraph 3.

Before:

> A conformal factor changes no null direction, so the light cones of this spacetime are Minkowski's, and so is its causal structure, one point short.
> What the factor changes is the proper time: an observer at rest in the chart ages $d\tau = \Omega\,dt$.
> Put the computer on the $t$ axis, $x = y = z = 0$, from $t = -T$ up to the missing origin.
> That worldline has no future endpoint in this spacetime, because its endpoint was removed, and its proper time $\int_{-T}^{0}\Omega\,dt$ is finite in Minkowski space and infinite with the factor as soon as $\Omega$ grows at least as fast as $1/|t|$ along it.

After:

> A conformal factor changes no null direction, so the light cones of this spacetime are Minkowski's, and so is its causal structure, one point short.
> The factor changes the proper time instead: an observer at rest in the chart ages $d\tau = \Omega\,dt$.
> We put the computer on the $t$ axis, $x = y = z = 0$, from $t = -T$ up to the missing origin.
> That worldline has no future endpoint in this spacetime, because its endpoint was removed, and its proper time $\int_{-T}^{0}\Omega\,dt$ is finite in Minkowski space and infinite with the factor as soon as $\Omega$ grows at least as fast as $1/|t|$ along it.

### conventions: `malament_hogarth/convention[3]`

Source: `MFS/assets/data/metrics/malament_hogarth.json`, `convention`, paragraph 4.

Before:

> Any event $p$ in the future light cone of the origin, such as $(t, 0, 0, 0)$ with $t > 0$, has the whole of that worldline in its past, since a signal sent from anywhere on it can pass around the missing point, and an observer can reach $p$ in finite proper time along any worldline that keeps away from that point, where $\Omega$ stays bounded.
> That is the Malament-Hogarth property, and $p$ is the event at which an eternity of computation has been delivered.
> It comes at a price: light the computer sends from where the factor is $\Omega$ reaches an observer at rest outside $C$ with its frequency raised by that same factor, $\omega_{\text{received}}/\omega_{\text{sent}} = \Omega$, which grows without bound as the computer nears the missing point.
> Earman and Norton add that $\Omega$ can be chosen so that the computer falls freely, by making the $t$ axis an axis of symmetry of $\Omega$, and the geodesic equations show why: on the axis every $\partial_x\Omega$, $\partial_y\Omega$ and $\partial_z\Omega$ vanishes, the spatial equations are solved by $\dot{x} = \dot{y} = \dot{z} = 0$, and the time equation, $\ddot{t} + \left(\partial_t\Omega/\Omega\right)\dot{t}^2 = 0$, only fixes the affine parameter.

After:

> Any event $p$ in the future light cone of the origin, such as $(t, 0, 0, 0)$ with $t > 0$, has the whole of that worldline in its past, since a signal sent from anywhere on it can pass around the missing point, and an observer can reach $p$ in finite proper time along any worldline that keeps away from that point, where $\Omega$ stays bounded.
> That is the Malament-Hogarth property, and $p$ is the event at which an eternity of computation has been delivered.
> It comes at a price: light the computer sends from where the factor is $\Omega$ reaches an observer at rest outside $C$ with its frequency raised by that same factor, $\omega_{\text{received}}/\omega_{\text{sent}} = \Omega$, which grows without bound as the computer nears the missing point.
> Earman and Norton add that $\Omega$ can be chosen so that the computer falls freely, by making the $t$ axis an axis of symmetry of $\Omega$, and the reason is in the geodesic equations: on the axis every $\partial_x\Omega$, $\partial_y\Omega$, and $\partial_z\Omega$ vanishes, the spatial equations are solved by $\dot{x} = \dot{y} = \dot{z} = 0$, and the time equation, $\ddot{t} + \left(\partial_t\Omega/\Omega\right)\dot{t}^2 = 0$, only fixes the affine parameter.

### conventions: `malament_hogarth/convention[4]`

Source: `MFS/assets/data/metrics/malament_hogarth.json`, `convention`, paragraph 5.

Before:

> The Weyl tensor vanishes identically, as it does for every metric conformal to a flat one, so all of the curvature is Ricci curvature, and the Ricci scalar is $R = -6\,\Box\Omega/\Omega^3$ with $\Box = \partial_x^2 + \partial_y^2 + \partial_z^2 - \partial_t^2$ the flat wave operator of the chart.
> Nothing in the construction asks that curvature to be the Einstein tensor of anything reasonable.
> Read through $G_{\mu\nu} = 8\pi Gc^{-4}T_{\mu\nu}$, the matter it calls for lives wherever $\Omega$ varies and is whatever the derivatives of $\Omega$ make it, and Earman and Norton say so plainly: the toy can be read as a solution by defining its source to be $G_{\mu\nu}$, with no guarantee that even the weak energy condition holds.
> It shows what the property is rather than where it might be found; anti-de Sitter space and the charged and rotating black holes have the property while solving the field equations with honest sources.

After:

> The Weyl tensor vanishes identically, as it does for every metric conformal to a flat one, so all of the curvature is Ricci curvature, and the Ricci scalar is $R = -6\,\Box\Omega/\Omega^3$ with $\Box = \partial_x^2 + \partial_y^2 + \partial_z^2 - \partial_t^2$ the flat wave operator of the chart.
> The construction does not require that curvature to be the Einstein tensor of anything reasonable.
> Through $G_{\mu\nu} = 8\pi Gc^{-4}T_{\mu\nu}$, the matter it calls for lives wherever $\Omega$ varies and is whatever the derivatives of $\Omega$ make it, as Earman and Norton say plainly: the toy becomes a solution once its source is defined to be $G_{\mu\nu}$, with no guarantee that even the weak energy condition holds.
> The toy exhibits the property, and the spacetimes where it is found are others: anti-de Sitter space and the charged and rotating black holes have the property while solving the field equations with honest sources.

### embedding diagram note: `embedding/malament_hogarth/plane.input`

Source: `_tools/derivations/embedding.py` line 3143, `in malament_hogarth() input=`.

Before:

> $\Omega = 1 + e^{1 - 1/(1 - \varrho^2)}/\varrho$ for $\varrho^2 = c^2t^2 + x^2 + y^2 + z^2 < 1$ and $\Omega = 1$ beyond, as the spacetime diagram declares.

After:

> $\Omega = 1 + e^{1 - 1/(1 - \varrho^2)}/\varrho$ for $\varrho^2 = c^2t^2 + x^2 + y^2 + z^2 < 1$ and $\Omega = 1$ beyond, as in the spacetime diagram.

### embedding diagram note: `embedding/malament_hogarth/plane.settings`

Source: `_tools/derivations/embedding.py` line 3141, `in malament_hogarth() settings=`.

Before:

> $c\,t$ and every length in the unit the declared $\Omega$ is written in, the radius of the region where $\Omega > 1$; each moment is the plane $z = 0$ of one $t$.

After:

> $c\,t$ and every length in the unit in which the declared $\Omega$ is written, the radius of the region where $\Omega > 1$; each moment is the plane $z = 0$ of one $t$.

## Minkowski Flat Spacetime

### history: `minkowski/history[0]`

Source: `MFS/assets/data/metrics/minkowski.json`, `history`, paragraph 1.

Before:

> Hermann Minkowski was Einstein's former mathematics professor at Zürich, and he understood what Einstein had done in 1905 better than Einstein did, at least geometrically.
> In 1908 he showed that the Lorentz transformations of special relativity are not ad hoc kinematic rules but rotations in a pseudo-Riemannian manifold of four dimensions with metric signature $(-,+,+,+)$ [minkowski1908]. "Henceforth," he declared in his famous Cologne lecture [minkowski1909], "space by itself and time by itself are doomed to fade away into mere shadows, and only a kind of union of the two will preserve an independent reality." Greg Egan built the universe of his Orthogonal trilogy by turning the one minus sign in Minkowski's signature into a plus [egan2010plusminus, egan2011clockwork, egan2012eternal, egan2013arrows].

After:

> Hermann Minkowski was Einstein's former mathematics professor at Zürich, and he understood what Einstein had done in 1905 better than Einstein did, at least geometrically.
> In 1908 he showed that the Lorentz transformations of special relativity, which had looked like ad hoc kinematic rules, are rotations in a pseudo-Riemannian manifold of four dimensions with metric signature $(-,+,+,+)$ [minkowski1908]. "Henceforth," he declared in his famous Cologne lecture [minkowski1909], "space by itself and time by itself are doomed to fade away into mere shadows, and only a kind of union of the two will preserve an independent reality." Greg Egan built the universe of his Orthogonal trilogy by turning the one minus sign in Minkowski's signature into a plus [egan2010plusminus, egan2011clockwork, egan2012eternal, egan2013arrows].

### history: `minkowski/history[4]`

Source: `MFS/assets/data/metrics/minkowski.json`, `history`, paragraph 5.

Before:

> The deep significance of this horizon became clear with the Unruh effect, derived independently by Davies in 1975 and Unruh in 1976 [davies1975, unruh1976]: an observer accelerating through the Minkowski vacuum detects a thermal bath of particles at temperature $T = \hbar a / 2\pi c k_B$.
> This result, that acceleration and temperature are two sides of the same coin, established the connection between horizons and thermodynamics that continues to shape our understanding of quantum gravity.
> Stephen Fulling had shown in 1973 that quantizing a field in the Rindler wedge, where the metric is static in the accelerated coordinates, gives a particle interpretation different from the standard one of flat spacetime [fulling1973].
> Geoffrey Sewell proved in 1982, from an analogue of a theorem of Joseph Bisognano and Eyvind Wichmann, that any quantum field, interacting or not, restricted to a region bounded by event horizons is in thermal equilibrium, a rigorous and model independent proof of a generalized Hawking-Unruh effect [bisognano1976, sewell1982].

After:

> This horizon turned out to be thermal, through the Unruh effect, derived independently by Paul Davies in 1975 and William Unruh in 1976 [davies1975, unruh1976]: an observer accelerating through the Minkowski vacuum detects a thermal bath of particles at temperature $T = \hbar a / 2\pi c k_B$.
> This result, which ties temperature to acceleration, established the connection between horizons and thermodynamics that continues to shape our understanding of quantum gravity.
> Stephen Fulling had shown in 1973 that quantizing a field in the Rindler wedge, where the metric is static in the accelerated coordinates, gives a particle interpretation different from the standard one of flat spacetime [fulling1973].
> Geoffrey Sewell proved in 1982, from an analogue of a theorem of Joseph Bisognano and Eyvind Wichmann, that any quantum field, interacting or not, restricted to a region bounded by event horizons is in thermal equilibrium, a rigorous and model independent proof of a generalized Hawking-Unruh effect [bisognano1976, sewell1982].

### conventions: `minkowski/convention[0]`

Source: `MFS/assets/data/metrics/minkowski.json`, `convention`, paragraph 1.

Before:

> Factors of $c$ are kept explicit throughout.
> The Cartesian and spherical charts use a time coordinate $t$, the null charts a retarded time $u$ and an advanced time $v$, and the Rindler chart a time $T$; all four carry dimensions of time, and $a$ is a proper acceleration.
> Components are taken in the chart $x^0 = ct$, so a component whose index is written on one of those letters is a component of that chart rather than of the bare coordinate.
> The dots in the geodesic equations are derivatives of that same chart with respect to an affine parameter, so $\dot{t}$ means $d(ct)/d\lambda$ and $\dot{u}$ means $d(cu)/d\lambda$, and every term of every equation carries the dimensions of its left hand side.

After:

> We keep factors of $c$ explicit throughout.
> The Cartesian and spherical charts use a time coordinate $t$, the null charts a retarded time $u$ and an advanced time $v$, and the Rindler chart a time $T$; all four carry dimensions of time, and $a$ is a proper acceleration.
> We take components in the chart $x^0 = ct$, so a component whose index we write with one of those letters is a component of that chart and not of the bare coordinate.
> A dot in the geodesic equations is a derivative of that chart's coordinates with respect to an affine parameter, so $\dot{t}$ means $d(ct)/d\lambda$ and $\dot{u}$ means $d(cu)/d\lambda$, and every term of every equation carries the dimensions of its left hand side.

### spacetime diagram caption: `diagrams/minkowski/systems/rindler/tx.caption[1]`

Source: `_tools/derivations/null_rays.py` line 646, `CAPTIONS ("minkowski", "rindler", "tx")`.

Before:

> That line is the Rindler horizon, which only the accelerated observers have.
> The spacetime is flat, its Kretschmann scalar is zero, and the rays carry on across $X = 0$ into the rest of Minkowski spacetime, which the Cartesian chart covers whole.

After:

> The edge of the wedge is the Rindler horizon, which only the accelerated observers have.
> The spacetime is flat, its Kretschmann scalar is zero, and the rays carry on across $X = 0$ into the rest of Minkowski spacetime, which the Cartesian chart covers whole.

### embedding diagram caption: `embedding/minkowski/plane.caption[1]`

Source: `_tools/derivations/embedding.py` line 3569, `CAPTIONS ("minkowski", "plane")`.

Before:

> Every other embedding diagram is measured against this one.
> Where a circle's circumference falls short of $2\pi$ times the distance out to it, as around a star, the plane curves into a bowl; where it exceeds it, as in anti-de Sitter space, no surface of revolution in flat space carries the plane at all.
> Hermann Minkowski set out in 1908 the geometry in which space at one moment of any inertial observer is this flat space of Euclid.

After:

> Every other slice is measured against this flat plane.
> Where a circle's circumference falls short of $2\pi$ times the distance out to it, as around a star, the plane curves into a bowl; where it exceeds it, as in anti-de Sitter space, no surface of revolution in flat space carries the plane.
> Hermann Minkowski set out in 1908 the geometry in which space at one moment of any inertial observer is this flat space of Euclid.

## Mixmaster Chaotic Cosmology

### history: `mixmaster/history[1]`

Source: `MFS/assets/data/metrics/mixmaster.json`, `history`, paragraph 2.

Before:

> Evgeny Lifshitz and Isaak Khalatnikov reached the same oscillation in Moscow from the other direction, at the cost of seven years of work.
> They had argued in 1963 that the general solution of the Einstein equations carries no singularity, the known ones being a penalty of symmetry [lifshitz1963].
> Penrose's theorem of 1965 made that position uninhabitable [penrose1965].
> With Vladimir Belinskii they went back and in 1970 published retraction and replacement together: a general solution that does reach a singularity, through Mixmaster oscillations [khalatnikov1970, bkl1970].
> The replacement was the bolder claim.
> Near a spacelike singularity the time derivatives drown out the spatial ones, so every point of a generic spacetime collapses on its own private Bianchi IX schedule, set by one number that sheds a unit at each bounce until it drops below one, whereupon the next era opens at its reciprocal.

After:

> Evgeny Lifshitz and Isaak Khalatnikov reached the same oscillation in Moscow from the other direction, at the cost of seven years of work.
> They had argued in 1963 that the general solution of the Einstein equations carries no singularity, the known singularities being consequences of symmetry [lifshitz1963].
> Roger Penrose's theorem of 1965 refuted that position [penrose1965].
> With Vladimir Belinskii they went back and in 1970 published retraction and replacement together: a general solution that does reach a singularity, through Mixmaster oscillations [khalatnikov1970, bkl1970].
> The replacement was the bolder claim.
> Near a spacelike singularity the time derivatives dominate the spatial ones, so every point of a generic spacetime collapses on its own private Bianchi IX schedule, set by one number that sheds a unit at each bounce until it drops below one, whereupon the next era opens at its reciprocal.

### history: `mixmaster/history[3]`

Source: `MFS/assets/data/metrics/mixmaster.json`, `history`, paragraph 4.

Before:

> By then the model had failed its own purpose: the oscillations do not mix the sky into causal contact [doroshkevich1971], and the horizon problem waited for inflation [guth1981].
> John Barrow gave the bounce map a positive Lyapunov exponent in 1982 [barrow1982], but an exponent measures separation per unit of time and relativity supplies no preferred unit: change the clock and other coordinates return zero [francisco1988, berger1991].
> Neil Cornish and Janna Levin ended a decade of argument in 1997 by declining to measure time at all: the phase space carries a fractal, and a fractal does not care what clock you hold [cornish1997].
> Adilson Motter settled the principle in 2003: positive Lyapunov exponents stay positive under a change of coordinates [motter2003].

After:

> By then the model had failed its own purpose: the oscillations do not mix the sky into causal contact [doroshkevich1971], and the horizon problem waited for inflation [guth1981].
> John Barrow gave the bounce map a positive Lyapunov exponent in 1982 [barrow1982], but an exponent measures separation per unit of time and relativity supplies no preferred unit: change the clock and other coordinates return zero [francisco1988, berger1991].
> Neil Cornish and Janna Levin ended a decade of argument in 1997 by characterizing the chaos without any measure of time: the phase space carries a fractal, and a fractal is the same whatever clock one uses [cornish1997].
> Adilson Motter settled the principle in 2003: positive Lyapunov exponents stay positive under a change of coordinates [motter2003].

### history: `mixmaster/history[4]`

Source: `MFS/assets/data/metrics/mixmaster.json`, `history`, paragraph 5.

Before:

> Misner's billiard reaches beyond four dimensions.
> Pure gravity stops bouncing in eleven spacetime dimensions and above [demaret1985], awkwardly, since eleven is where supergravity lives; restore its three form field and the walls come back, the table now the Weyl chamber of the hyperbolic Kac-Moody algebra $E_{10}$ [damour2003].
> The appliance Misner named in 1969 is the shape of the singularity in M-theory.

After:

> Misner's billiard reaches beyond four dimensions.
> Pure gravity stops bouncing in eleven spacetime dimensions and above [demaret1985], awkwardly, since eleven is where supergravity lives; restore its three form field and the walls come back, the table now the Weyl chamber of the hyperbolic Kac-Moody algebra $E_{10}$ [damour2003].
> The billiard of the Mixmaster, which Misner named in 1969, thus describes the approach to the singularity in M-theory.

### conventions: `mixmaster/convention[0]`

Source: `MFS/assets/data/metrics/mixmaster.json`, `convention`, paragraph 1.

Before:

> Coordinates are $(t, \psi, \theta, \phi)$ with $t$ a cosmic time and $\psi$, $\theta$ and $\phi$ Euler angles on a three sphere, $\psi \in [0, 4\pi)$, $\theta \in [0, \pi]$ and $\phi \in [0, 2\pi)$, and factors of $c$ are kept explicit.
> Each slice of constant $t$ is closed, and the line element is written in the three one forms that the symmetry of the slice leaves unchanged, $\omega^1 = \cos\psi\,d\theta + \sin\psi\sin\theta\,d\phi$, $\omega^2 = \sin\psi\,d\theta - \cos\psi\sin\theta\,d\phi$ and $\omega^3 = d\psi + \cos\theta\,d\phi$, as $ds^2 = -c^2dt^2 + a_1^2(\omega^1)^2 + a_2^2(\omega^2)^2 + a_3^2(\omega^3)^2$.
> They obey $d\omega^1 = \omega^2\wedge\omega^3$ and its two cyclic partners, which is what makes the slices type IX in Bianchi's list, the type whose group is that of rotations.
> This is the form Misner wrote in 1969 and Belinskii, Khalatnikov and Lifshitz in 1970, with their $a$, $b$ and $c$ written $a_1$, $a_2$ and $a_3$; Misner's variables are the logarithms, $a_1 = e^{\alpha + \beta_+ + \sqrt{3}\beta_-}$, $a_2 = e^{\alpha + \beta_+ - \sqrt{3}\beta_-}$ and $a_3 = e^{\alpha - 2\beta_+}$, with $e^{\alpha}$ carrying the length.

After:

> We use coordinates $(t, \psi, \theta, \phi)$, with $t$ a cosmic time and $\psi$, $\theta$, and $\phi$ Euler angles on a three sphere, $\psi \in [0, 4\pi)$, $\theta \in [0, \pi]$, and $\phi \in [0, 2\pi)$, and keep factors of $c$ explicit.
> Each slice of constant $t$ is closed, and we write the line element in the three one forms that the symmetry of the slice leaves unchanged, $\omega^1 = \cos\psi\,d\theta + \sin\psi\sin\theta\,d\phi$, $\omega^2 = \sin\psi\,d\theta - \cos\psi\sin\theta\,d\phi$, and $\omega^3 = d\psi + \cos\theta\,d\phi$, as $ds^2 = -c^2dt^2 + a_1^2(\omega^1)^2 + a_2^2(\omega^2)^2 + a_3^2(\omega^3)^2$.
> They obey $d\omega^1 = \omega^2\wedge\omega^3$ and its two cyclic partners, which makes the slices type IX in Bianchi's list, the type whose group is that of rotations.
> This is the form Charles Misner wrote in 1969, and Vladimir Belinskii, Isaak Khalatnikov, and Evgeny Lifshitz in 1970, with their $a$, $b$, and $c$ written $a_1$, $a_2$, and $a_3$; Misner's variables are the logarithms, $a_1 = e^{\alpha + \beta_+ + \sqrt{3}\beta_-}$, $a_2 = e^{\alpha + \beta_+ - \sqrt{3}\beta_-}$, and $a_3 = e^{\alpha - 2\beta_+}$, with $e^{\alpha}$ carrying the length.

### conventions: `mixmaster/convention[1]`

Source: `MFS/assets/data/metrics/mixmaster.json`, `convention`, paragraph 2.

Before:

> The three scale factors are left free functions of $t$ and are lengths, since the angles carry none.
> Where all three equal $a$ the slice is a round three sphere of radius $2a$ and the metric is closed Friedmann, and where two agree the slices keep a further rotation about the third direction, the symmetry of Taub's vacuum solution of 1951.
> In the coordinate basis the components carry $\psi$ only because $\omega^1$ and $\omega^2$ turn with it; the geometry is the same at every point of a slice.

After:

> We leave the three scale factors free functions of $t$; they are lengths, since the angles carry none.
> Where all three equal $a$ the slice is a round three sphere of radius $2a$ and the metric is closed Friedmann, and where two agree the slices keep a further rotation about the third direction, the symmetry of Abraham Taub's vacuum solution of 1951.
> In the coordinate basis the components carry $\psi$ only because $\omega^1$ and $\omega^2$ turn with it; the geometry is the same at every point of a slice.

### conventions: `mixmaster/convention[2]`

Source: `MFS/assets/data/metrics/mixmaster.json`, `convention`, paragraph 3.

Before:

> Components are taken in the chart $x^0 = ct$, where a chart component of the metric carries no dimensions.
> A prime is the derivative along that chart time, $a_i' = da_i/d(ct)$, so $a_i'$ is dimensionless and $a_i''$ carries an inverse length, and the dots in the geodesic equations are velocities of the same chart, so $\dot{t}$ means $d(ct)/d\lambda$.
> The signature is $(-,+,+,+)$ and the Ricci tensor is the standard contraction $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$, which the Ricci scalar and the Einstein tensor follow.

After:

> We take components in the chart $x^0 = ct$, where a chart component of the metric carries no dimensions.
> A prime is the derivative along that chart time, $a_i' = da_i/d(ct)$, so $a_i'$ is dimensionless and $a_i''$ carries an inverse length, and the dots in the geodesic equations are velocities in the same chart, so $\dot{t}$ means $d(ct)/d\lambda$.
> We work in the signature $(-,+,+,+)$, and the Ricci tensor is the standard contraction, $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$, from which the Ricci scalar and the Einstein tensor follow.

### conventions: `mixmaster/convention[3]`

Source: `MFS/assets/data/metrics/mixmaster.json`, `convention`, paragraph 4.

Before:

> No component assumes a field equation.
> Every component is an identity for any three scale factors, and the Mixmaster is the vacuum among them.
> In the frame $c\,dt$, $a_1\omega^1$, $a_2\omega^2$ and $a_3\omega^3$ the Ricci tensor is diagonal, and the vacuum equations read $\dfrac{\left(a_1'a_2a_3\right)'}{a_1a_2a_3} = \dfrac{\left(a_2^2 - a_3^2\right)^2 - a_1^4}{2a_1^2a_2^2a_3^2}$ with its two cyclic partners, together with the constraint $G^t{}_t = 0$, $\dfrac{a_1'a_2'}{a_1a_2} + \dfrac{a_1'a_3'}{a_1a_3} + \dfrac{a_2'a_3'}{a_2a_3} = \dfrac{a_1^4 + a_2^4 + a_3^4 - 2a_1^2a_2^2 - 2a_1^2a_3^2 - 2a_2^2a_3^2}{4a_1^2a_2^2a_3^2}$.
> Drop the right hand sides and these are the Kasner equations, whose solutions are the Kasner epochs; the right hand sides are the curvature of the closed slice, the walls of Misner's billiard, negligible while one scale factor is much larger than the other two and decisive when the growing one catches them, which is the bounce that deals out the next set of exponents.
> The line element can be written down whatever the solution does; the chaos is in the solutions of these equations, not in the metric that carries them.

After:

> No component assumes a field equation.
> Every component is an identity for any three scale factors, and the Mixmaster is the vacuum among them.
> In the frame $c\,dt$, $a_1\omega^1$, $a_2\omega^2$, and $a_3\omega^3$ the Ricci tensor is diagonal, and the vacuum equations read $\dfrac{\left(a_1'a_2a_3\right)'}{a_1a_2a_3} = \dfrac{\left(a_2^2 - a_3^2\right)^2 - a_1^4}{2a_1^2a_2^2a_3^2}$ with its two cyclic partners, together with the constraint $G^t{}_t = 0$, $\dfrac{a_1'a_2'}{a_1a_2} + \dfrac{a_1'a_3'}{a_1a_3} + \dfrac{a_2'a_3'}{a_2a_3} = \dfrac{a_1^4 + a_2^4 + a_3^4 - 2a_1^2a_2^2 - 2a_1^2a_3^2 - 2a_2^2a_3^2}{4a_1^2a_2^2a_3^2}$.
> Without the right hand sides these are the Kasner equations, whose solutions are the Kasner epochs; the right hand sides are the curvature of the closed slice, the walls of Misner's billiard, negligible while one scale factor is much larger than the other two and decisive when the growing one catches them, in the bounce that deals out the next set of exponents.
> The line element holds whatever the solution does, and the chaos lives in the solutions of these equations.

### embedding diagram caption: `embedding/mixmaster/sphere.caption[0]`

Source: `_tools/derivations/embedding.py` line 3603, `CAPTIONS ("mixmaster", "sphere")`.

Before:

> This is the great two sphere of the Mixmaster universe's three sphere at five moments of its proper time $\tau$, each drawn as a surface in flat space so that every distance along it is the metric distance.
> Every great sphere of a Mixmaster slice is congruent to every other, and the three great circles in which it meets its planes of symmetry have circumferences $4\pi a_1$, $4\pi a_2$ and $4\pi a_3$, so the scale factors can be read off it.
> When all three agree it is a round sphere of radius $2a$, the equator of a round three sphere.

After:

> This is the great two sphere of the Mixmaster universe's three sphere at five moments of its proper time $\tau$, each drawn as a surface in flat space so that every distance along it is the metric distance.
> Every great sphere of a Mixmaster slice is congruent to every other, and the three great circles in which it meets its planes of symmetry have circumferences $4\pi a_1$, $4\pi a_2$, and $4\pi a_3$, so the sphere carries all three scale factors.
> When all three agree it is a round sphere of radius $2a$, the equator of a round three sphere.

### embedding diagram, not drawn: `embedding/mixmaster/sphere.stops[1]`

Source: `_tools/derivations/embedding.py` line 3277, `in mixmaster() stops=`.

Before:

> At the second moment $a_3$ is more than $2/\sqrt{3}$ times $a_1$, the curvature about the poles is negative, and only the band about the equator is drawn.

After:

> At the second moment $a_3$ is more than $2/\sqrt{3}$ times $a_1$, the curvature about the poles is negative, and only the band about the equator has a surface of revolution in flat space.

## Morris-Thorne Traversable Wormhole

### history: `morris_thorne/history[0]`

Source: `MFS/assets/data/metrics/morris_thorne.json`, `history`, paragraph 1.

Before:

> Carl Sagan wanted his heroine's route across the galaxy to survive a physicist's reading.
> In the summer of 1985 he sent Kip Thorne a prepublication draft of Contact and asked him to make the gravitational physics in it as accurate as possible [sagan1985, morris1988].
> Science fiction moved people by dropping them into black holes, and a black hole is no route at all: tidal forces kill the traveler before the horizon unless the hole outweighs ten thousand suns, and the horizon is a membrane that admits and never returns.
> Wormholes were no better off.
> The bridge Einstein and Nathan Rosen ran between two sheets of spacetime in 1935 [einstein1935] pinches off faster than light can cross it, as Robert Fuller and John Wheeler showed in 1962 [fuller1962].
> There was nothing on the shelf to hand a novelist.

After:

> Carl Sagan wanted his heroine's route across the galaxy to survive a physicist's reading.
> In the summer of 1985 he sent Kip Thorne a prepublication draft of Contact and asked him to make the gravitational physics in it as accurate as possible [sagan1985, morris1988].
> Science fiction moved people by dropping them into black holes, and a black hole is no route: tidal forces kill the traveler before the horizon unless the hole outweighs ten thousand suns, and the horizon is a membrane that admits and never returns.
> Wormholes were no better off.
> The bridge Einstein and Nathan Rosen ran between two sheets of spacetime in 1935 [einstein1935] pinches off faster than light can cross it, as Robert Fuller and John Wheeler showed in 1962 [fuller1962].
> There was nothing on the shelf to hand a novelist.

### history: `morris_thorne/history[1]`

Source: `MFS/assets/data/metrics/morris_thorne.json`, `history`, paragraph 2.

Before:

> So Thorne and his graduate student Michael Morris ran the problem backwards.
> The usual direction is source first: write the stress-energy, solve for the geometry.
> They began instead from what a traveler would insist on.

After:

> So Thorne and his graduate student Michael Morris ran the problem backwards.
> The usual direction is source first: write the stress energy, solve for the geometry.
> They began instead from what a traveler would insist on.

### history: `morris_thorne/history[3]`

Source: `MFS/assets/data/metrics/morris_thorne.json`, `history`, paragraph 4.

Before:

> At the throat the radial tension has to beat the energy density, and Morris and Thorne worked the number for a throat of three kilometers radius: about $10^{37}$ dynes per square centimeter, the same magnitude as the pressure at the center of the most massive of neutron stars [morris1988].
> Tension beating energy density means that an observer crossing the throat fast enough, near the speed of light, measures a negative energy density, so the null energy condition fails and the singularity theorems lose their footing.
> They called the stuff exotic and declined to be embarrassed, since quantum field theory, they wrote, gives tantalizing hints that such material might, in fact, be possible.
> Nor was the requirement a quirk of their ansatz.
> Topological censorship, proved in 1993 by John Friedman, Kristin Schleich and Donald Witt, holds that under the null energy condition any topological shortcut collapses too fast for light to cross, so a shortcut that stays open costs exotic matter [friedman1993].

After:

> At the throat the radial tension has to beat the energy density, and Morris and Thorne worked the number for a throat of three kilometers radius: about $10^{37}$ dynes per square centimeter, the same magnitude as the pressure at the center of the most massive of neutron stars [morris1988].
> Tension beating energy density means that an observer crossing the throat fast enough, near the speed of light, measures a negative energy density, so the null energy condition fails and the singularity theorems no longer apply.
> They called the stuff exotic and declined to be embarrassed, since quantum field theory, they wrote, gives tantalizing hints that such material might, in fact, be possible.
> Nor was the requirement a quirk of their ansatz.
> Topological censorship, proved in 1993 by John Friedman, Kristin Schleich, and Donald Witt, holds that under the null energy condition any topological shortcut collapses too fast for light to cross, so a shortcut that stays open costs exotic matter [friedman1993].

### history: `morris_thorne/history[4]`

Source: `MFS/assets/data/metrics/morris_thorne.json`, `history`, paragraph 5.

Before:

> Nobody had yet asked what happens if you take one mouth of the shortcut and move it.
> Put one mouth on a ship, run it out at speed, bring it home.
> It has aged less than the mouth that stayed, the tunnel keeps its ends in step as seen from inside, and a walk through the throat is a walk into your own past.
> Morris, Thorne and Ulvi Yurtsever published that the same year, and the arguing began [morris1988prl].

After:

> Nobody had yet asked what happens if you take one mouth of the shortcut and move it.
> Put one mouth on a ship, run it out at speed, bring it home.
> It has aged less than the mouth that stayed, the tunnel keeps its ends in step as seen from inside, and a walk through the throat is a walk into your own past.
> Morris, Thorne, and Ulvi Yurtsever published that the same year, and the arguing began [morris1988prl].

### history: `morris_thorne/history[5]`

Source: `MFS/assets/data/metrics/morris_thorne.json`, `history`, paragraph 6.

Before:

> Fernando Echeverria, Gunnar Klinkhammer and Thorne threw billiard balls in so that each came back out to knock its younger self off course, and found in 1991 that the grandfather paradox fails in the direction nobody expected: not too few consistent histories but too many [echeverria1991].
> Sung-Won Kim and Thorne found the quantum stress-energy diverging as the machine switches on, but so feebly that quantum gravity takes over first, they argued, and the wormhole carries on regardless [kim1991].
> Stephen Hawking held that the divergence wins, and answered in 1992 with chronology protection, his experimental evidence being that we have not been overrun by tourists from the future [hawking1992].
> They never agreed on what counts as close to a horizon, which is really a question about quantum gravity.
> Sagan's request for a way across the galaxy had become a research program.

After:

> Fernando Echeverria, Gunnar Klinkhammer, and Thorne threw billiard balls in so that each came back out to knock its younger self off course, and found in 1991 that the grandfather paradox fails in the direction nobody expected, with too many consistent histories where too few had been feared [echeverria1991].
> Sung-Won Kim and Thorne found the quantum stress energy diverging as the machine switches on, but so feebly that quantum gravity takes over first, they argued, and the wormhole carries on regardless [kim1991].
> Stephen Hawking held that the divergence wins, and answered in 1992 with chronology protection, his experimental evidence being that we have not been overrun by tourists from the future [hawking1992].
> They never agreed on what counts as close to a horizon, which is a question about quantum gravity.
> Sagan's request for a way across the galaxy had become a research program.

### conventions: `morris_thorne/convention[0]`

Source: `MFS/assets/data/metrics/morris_thorne.json`, `convention`, paragraph 1.

Before:

> Coordinates are $(t, r, \theta, \phi)$ in the areal chart and $(t, l, \theta, \phi)$ in the proper distance chart, with $t$ carrying dimensions of time, $r$ the areal radius, whose sphere has circumference $2\pi r$, and $l$ the proper radial distance, both lengths.
> The redshift function $\Phi$ is dimensionless, since it sits inside an exponential, and the shape function $b$ is a length beside $r$, which is what leaves $1 - b/r$ dimensionless.
> Components are taken in the chart $x^0 = ct$, so every component is a component of that chart even though its index is written with the bare letter $t$; a radial derivative carries no factor of $c$, because $r$ and $l$ are lengths already.
> The dots in the geodesic equations are derivatives of that same chart with respect to an affine parameter, so $\dot{t}$ means $d(ct)/d\lambda$, and the equations are $\ddot{x}^\mu + \Gamma^\mu{}_{\nu\rho}\dot{x}^\nu\dot{x}^\rho = 0$ with the connection of that chart.
> The signature is $(-,+,+,+)$ and the Ricci tensor is the standard contraction $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$, which the Ricci scalar and the Einstein tensor follow.

After:

> We use coordinates $(t, r, \theta, \phi)$ in the areal chart and $(t, l, \theta, \phi)$ in the proper distance chart, with $t$ carrying dimensions of time, $r$ the areal radius, whose sphere has circumference $2\pi r$, and $l$ the proper radial distance, both lengths.
> The redshift function $\Phi$ is dimensionless, since it sits inside an exponential, and the shape function $b$ is a length beside $r$, so $1 - b/r$ is dimensionless.
> We take components in the chart $x^0 = ct$, so every component is a component of that chart even though we write its index with the bare letter $t$; a radial derivative carries no factor of $c$, because $r$ and $l$ are lengths already.
> A dot in the geodesic equations is a derivative of that chart's coordinates with respect to an affine parameter, so $\dot{t}$ means $d(ct)/d\lambda$, and the equations read $\ddot{x}^\mu + \Gamma^\mu{}_{\nu\rho}\dot{x}^\nu\dot{x}^\rho = 0$, with the connection of that chart.
> We work in the signature $(-,+,+,+)$, and the Ricci tensor is the standard contraction, $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$, from which the Ricci scalar and the Einstein tensor follow.

### conventions: `morris_thorne/convention[1]`

Source: `MFS/assets/data/metrics/morris_thorne.json`, `convention`, paragraph 2.

Before:

> The geometry is the general case: $\Phi$ and $b$ are left as arbitrary functions, so every component is an identity for any pair of them and no example is chosen.
> The price of the throat is read off the Einstein tensor.
> In the orthonormal frame of a static observer, $e_{\hat{t}} = e^{-\Phi}\partial_{ct}$, $e_{\hat{r}} = \sqrt{1 - b/r}\,\partial_r$, $e_{\hat{\theta}} = r^{-1}\partial_\theta$ and $e_{\hat{\phi}} = (r\sin\theta)^{-1}\partial_\phi$, the field equations $G_{\mu\nu} = 8\pi Gc^{-4}T_{\mu\nu}$ read the energy density and the radial tension straight out of $G_{\hat{t}\hat{t}} = \dfrac{\partial_r b}{r^2}$ and $G_{\hat{r}\hat{r}} = \dfrac{2(1 - b/r)\partial_r\Phi}{r} - \dfrac{b}{r^3}$, giving $\rho c^2 = \dfrac{c^4}{8\pi G}\dfrac{\partial_r b}{r^2}$ and $\tau = -p_r = \dfrac{c^4}{8\pi G}\left(\dfrac{b}{r^3} - \dfrac{2(1 - b/r)\partial_r\Phi}{r}\right)$.

After:

> We keep the general case: $\Phi$ and $b$ stay arbitrary functions, so every component is an identity for any pair of them, and we choose no example.
> The Einstein tensor sets the price of the throat.
> In the orthonormal frame of a static observer, $e_{\hat{t}} = e^{-\Phi}\partial_{ct}$, $e_{\hat{r}} = \sqrt{1 - b/r}\,\partial_r$, $e_{\hat{\theta}} = r^{-1}\partial_\theta$, and $e_{\hat{\phi}} = (r\sin\theta)^{-1}\partial_\phi$, the field equations $G_{\mu\nu} = 8\pi Gc^{-4}T_{\mu\nu}$ take the energy density and the radial tension straight from $G_{\hat{t}\hat{t}} = \dfrac{\partial_r b}{r^2}$ and $G_{\hat{r}\hat{r}} = \dfrac{2(1 - b/r)\partial_r\Phi}{r} - \dfrac{b}{r^3}$, giving $\rho c^2 = \dfrac{c^4}{8\pi G}\dfrac{\partial_r b}{r^2}$ and $\tau = -p_r = \dfrac{c^4}{8\pi G}\left(\dfrac{b}{r^3} - \dfrac{2(1 - b/r)\partial_r\Phi}{r}\right)$.

### conventions: `morris_thorne/convention[2]`

Source: `MFS/assets/data/metrics/morris_thorne.json`, `convention`, paragraph 3.

Before:

> At the throat $r = b_0$ the factor $1 - b/r$ vanishes and the $\partial_r\Phi$ term goes with it.
> That needs one word of care, because $\partial_r\Phi$ diverges at the throat: written against the proper distance $l$ of the second chart it goes as $\partial_l\Phi(0)/\left(\partial_l^2r(0)\,l\right)$ while $1 - b/r$ goes as $l^2$, so the product vanishes, and what bounds $\partial_l\Phi(0)$ is precisely the requirement that there be no horizon.
> So the redshift function leaves no trace at the throat, and $\tau_0 = \dfrac{c^4}{8\pi Gb_0^2}$ stands against $\rho_0c^2 = \dfrac{c^4}{8\pi Gb_0^2}\partial_r b(b_0)$, giving $\tau_0 - \rho_0c^2 = \dfrac{c^4}{8\pi Gb_0^2}\left(1 - \partial_r b(b_0)\right) > 0$ by the flare out condition.
> The radial tension at the throat beats the energy density there.
> That is the same inequality as the failure of the null energy condition: a radial null vector $k^\mu$ measures $T_{\mu\nu}k^\mu k^\nu \propto \rho c^2 + p_r = \rho c^2 - \tau < 0$, a negative energy density in the frame of a light beam crossing the throat.
> It holds for every wormhole of this form, not for some choice of $\Phi$ and $b$.

After:

> At the throat $r = b_0$ the factor $1 - b/r$ vanishes, and the $\partial_r\Phi$ term goes with it.
> That needs one word of care, because $\partial_r\Phi$ diverges at the throat: written against the proper distance $l$ of the second chart it goes as $\partial_l\Phi(0)/\left(\partial_l^2r(0)\,l\right)$ while $1 - b/r$ goes as $l^2$, so the product vanishes, and the requirement that there be no horizon bounds $\partial_l\Phi(0)$.
> The redshift function therefore leaves no trace at the throat, and $\tau_0 = \dfrac{c^4}{8\pi Gb_0^2}$ stands against $\rho_0c^2 = \dfrac{c^4}{8\pi Gb_0^2}\partial_r b(b_0)$, giving $\tau_0 - \rho_0c^2 = \dfrac{c^4}{8\pi Gb_0^2}\left(1 - \partial_r b(b_0)\right) > 0$ by the flare out condition.
> The radial tension at the throat beats the energy density there.
> The same inequality is the failure of the null energy condition: a radial null vector $k^\mu$ measures $T_{\mu\nu}k^\mu k^\nu \propto \rho c^2 + p_r = \rho c^2 - \tau < 0$, a negative energy density in the frame of a light beam crossing the throat.
> It holds for every wormhole of this form, whatever $\Phi$ and $b$ are.

### conventions: `morris_thorne/convention[3]`

Source: `MFS/assets/data/metrics/morris_thorne.json`, `convention`, paragraph 4.

Before:

> The areal chart is the one the field equations are cleanest in, but it is a chart of one side and it is degenerate at the throat, where $g_{rr}$ diverges and where the components carrying $r - b$ in a denominator have to be read as limits.
> The proper distance chart holds the throat as an ordinary point: $dl = \pm\dfrac{dr}{\sqrt{1 - b/r}}$ runs over the whole real line and covers both sides at once, $l = 0$ is the throat, and the areal radius $r(l)$ has its minimum $b_0$ there.
> In that chart flare out is the statement that the minimum is a real one, $\partial_l^2 r(0) > 0$, and the energy condition fails in a single line, because $G_{\hat{t}\hat{t}} + G_{\hat{r}\hat{r}} = \dfrac{2\left(\partial_l\Phi\,\partial_l r - \partial_l^2 r\right)}{r}$ everywhere, which at the throat, where $\partial_l r = 0$, is $-\dfrac{2\partial_l^2 r(0)}{b_0} < 0$.
> Nothing in either chart vanishes identically.
> Both charts are static and spherically symmetric, so the metric is diagonal and depends on the radial coordinate alone: every off diagonal component of the Ricci and Einstein tensors is zero, and the only Riemann and Weyl components that survive are those whose two index pairs are the same pair drawn from $tr$, $t\theta$, $t\phi$, $r\theta$, $r\phi$ and $\theta\phi$, since a component mixing two different pairs would have to change sign under the time reversal or the reflection the spacetime is invariant under.

After:

> The field equations are cleanest in the areal chart, but it is a chart of one side, and it is degenerate at the throat, where $g_{rr}$ diverges and where the components carrying $r - b$ in a denominator have to be taken as limits.
> The proper distance chart holds the throat as an ordinary point: $dl = \pm\dfrac{dr}{\sqrt{1 - b/r}}$ runs over the whole real line and covers both sides at once, $l = 0$ is the throat, and the areal radius $r(l)$ has its minimum $b_0$ there.
> In that chart flare out is the statement that the minimum is a real one, $\partial_l^2 r(0) > 0$, and the energy condition fails in a single line, because $G_{\hat{t}\hat{t}} + G_{\hat{r}\hat{r}} = \dfrac{2\left(\partial_l\Phi\,\partial_l r - \partial_l^2 r\right)}{r}$ everywhere, which at the throat, where $\partial_l r = 0$, is $-\dfrac{2\partial_l^2 r(0)}{b_0} < 0$.
> Nothing in either chart vanishes identically.
> Both charts are static and spherically symmetric, so the metric is diagonal and depends on the radial coordinate alone: every off diagonal component of the Ricci and Einstein tensors is zero, and the only Riemann and Weyl components that survive are those whose two index pairs are the same pair drawn from $tr$, $t\theta$, $t\phi$, $r\theta$, $r\phi$, and $\theta\phi$, since a component mixing two different pairs would have to change sign under the time reversal or the reflection that leaves the spacetime invariant.

### conventions: `morris_thorne/convention[4]`

Source: `MFS/assets/data/metrics/morris_thorne.json`, `convention`, paragraph 5.

Before:

> The spacetime is not a vacuum, so the Weyl tensor is not a copy of the Riemann tensor.
> It is of Petrov type D and carries a single scalar, the bracket every Weyl component is a multiple of, which vanishes only for the $\Phi$ and $b$ that make the geometry conformally flat.
> The Kretschmann scalar is a sum of four squares, one for each independent orthonormal curvature, so it is positive wherever the geometry is curved at all, and it stays finite at the throat, where the first of those squares is an indeterminate product in the areal chart and an ordinary number in the proper distance one.
> That is the sense in which this spacetime has no singularity to hide.

After:

> The spacetime is not a vacuum, so its Weyl tensor differs from the Riemann tensor.
> It is of Petrov type D and carries a single scalar, the bracket every Weyl component is a multiple of, which vanishes only for the $\Phi$ and $b$ that make the geometry conformally flat.
> The Kretschmann scalar is a sum of four squares, one for each independent orthonormal curvature, so it is positive wherever the geometry is curved, and it stays finite at the throat, where the first of those squares is an indeterminate product in the areal chart and an ordinary number in the proper distance one.
> This wormhole therefore has no curvature singularity.

### spacetime diagram caption: `diagrams/morris_thorne/systems/spherical/radial.caption[1]`

Source: `_tools/derivations/null_rays.py` line 619, `CAPTIONS ("morris_thorne", "spherical", "radial")`.

Before:

> In this areal chart the cones close toward the throat at $r = b_0$, as they would at a horizon, because $g_{rr} = (1 - b_0^2/r^2)^{-1}$ diverges there.
> But $g_{tt} = -1$ stays finite, so $\partial_t$ is timelike right up to the throat, and $r = b_0$ is only the edge of this chart.
> The rays reach the throat in finite $t$ and pass into the other mouth, which the Ellis-Bronnikov chart, running through the throat, covers in full.
> Below $b_0$ the formula gives a metric on this plane with no null directions at all.

After:

> In this areal chart the cones close toward the throat at $r = b_0$, as they would at a horizon, because $g_{rr} = (1 - b_0^2/r^2)^{-1}$ diverges there.
> But $g_{tt} = -1$ stays finite, so $\partial_t$ is timelike right up to the throat, and $r = b_0$ is only the edge of this chart.
> The rays reach the throat in finite $t$ and pass into the other mouth, which the Ellis-Bronnikov chart, running through the throat, covers in full.
> Below $b_0$ the formula gives a metric on this plane with no null directions.

## Natário Warp Drive

### history: `natario/history[0]`

Source: `MFS/assets/data/metrics/natario.json`, `history`, paragraph 1.

Before:

> The Alcubierre warp drive comes with a story attached, and the story is the part that gets remembered: space contracts in front of the ship, expands behind it, and the ship rides the wave from here to there [alcubierre1994].
> José Natário asked whether the expansion does any of the work.
> It does not.
> In 2002 he published a warp drive whose volume elements are neither stretched nor squeezed anywhere in the spacetime, and which carries a ship exactly as Alcubierre's does [natario2002].
> The contraction and expansion turned out to be a marginal consequence of one particular choice; make a different choice and they go away, while the drive keeps working.

After:

> The Alcubierre warp drive comes with a story attached, and the story is the part that gets remembered: space contracts in front of the ship, expands behind it, and the ship rides the wave from here to there [alcubierre1994].
> José Natário asked whether the expansion does any of the work, and showed that it does none.
> In 2002 he published a warp drive whose volume elements are neither stretched nor squeezed anywhere in the spacetime, and which carries a ship exactly as Alcubierre's does [natario2002].
> The contraction and expansion turned out to be a marginal consequence of one particular choice; make a different choice and they go away, while the drive keeps working.

### history: `natario/history[1]`

Source: `MFS/assets/data/metrics/natario.json`, `history`, paragraph 2.

Before:

> Natário put a prettier picture in their place.
> Hold a stone fixed in a river and watch the water: the flow parts ahead of it, slides around the sides, and closes up behind, and no parcel of water is compressed on the way past.
> That is the Natário drive.
> Space flows around a ship that sits still in it, and the flow is arranged so that volumes are carried along untouched.
> Alcubierre's treadmill of stretched and squeezed space can be dismantled without disturbing anything the drive actually needs.

After:

> Natário put a prettier picture in their place.
> Hold a stone fixed in a river and watch the water: the flow parts ahead of it, slides around the sides, and closes up behind, and no parcel of water is compressed on the way past.
> That is the Natário drive.
> Space flows around a ship that sits still in it, and the flow is arranged so that volumes are carried along untouched.
> The stretching and squeezing of Alcubierre's space can be removed without disturbing anything the drive needs.

### history: `natario/history[2]`

Source: `MFS/assets/data/metrics/natario.json`, `history`, paragraph 3.

Before:

> His drive still needs what Alcubierre's needed: negative energy.
> Francisco Lobo and Matt Visser examined both drives in 2004 and found the classical energy conditions violated with no perturbative approximation doing the damage, and violated still when the bubble is made to crawl; slowness buys nothing, and the negative energy stored in the warp field has to be a sizable fraction of the mass of the ship it carries [lobo2004].
> Jessica Santiago, Sebastian Schuster and Matt Visser sharpened this in 2022 into a claim about the whole family.
> Any warp drive a physicist would call physically reasonable violates the null energy condition, and so violates the weak, strong and dominant conditions along with it; move to modified gravity and the violation returns as a violation of the purely geometric null convergence condition [santiago2022].

After:

> His drive still needs what Alcubierre's needed: negative energy.
> Francisco Lobo and Matt Visser examined both drives in 2004 and found the classical energy conditions violated with no perturbative approximation doing the damage, and violated still when the bubble is made to crawl, and the negative energy stored in the warp field has to be a sizable fraction of the mass of the ship it carries [lobo2004].
> Jessica Santiago, Sebastian Schuster, and Matt Visser sharpened this in 2022 into a claim about the whole family.
> Any warp drive a physicist would call physically reasonable violates the null energy condition, and so violates the weak, strong, and dominant conditions along with it; move to modified gravity and the violation returns as a violation of the purely geometric null convergence condition [santiago2022].

### history: `natario/history[4]`

Source: `MFS/assets/data/metrics/natario.json`, `history`, paragraph 5.

Before:

> Brandon Mattingly and colleagues drew both drives in 2021, plotting curvature invariants for the Alcubierre and Natário spacetimes, which is a way of looking at a geometry with no coordinate system standing in front of it.
> In their plots a safe harbor sits at the center where the crew rides flat, with the shape function drawing out the bubble wall around them.
> At constant velocity the Natário drive leaves no wake; the wake appears once the bubble accelerates [mattingly2021].

After:

> Brandon Mattingly and colleagues drew both drives in 2021, plotting curvature invariants for the Alcubierre and Natário spacetimes, a way of looking at a geometry that does not depend on its coordinates.
> In their plots a safe harbor sits at the center where the crew rides flat, with the shape function drawing out the bubble wall around them.
> At constant velocity the Natário drive leaves no wake; the wake appears once the bubble accelerates [mattingly2021].

### history: `natario/history[5]`

Source: `MFS/assets/data/metrics/natario.json`, `history`, paragraph 6.

Before:

> Natário did not settle whether anyone can travel this way.
> He settled which part of the story we were entitled to believe.
> The expansion of space was scenery, and taking it down leaves the warp drive standing where it stood, with the same debt in negative energy still unpaid.

After:

> Natário did not settle whether anyone can travel this way.
> He showed that the expansion of space does none of the work of a warp drive.
> Without it the drive still carries the ship, and it still needs negative energy.

### conventions: `natario/convention[0]`

Source: `MFS/assets/data/metrics/natario.json`, `convention`, paragraph 1.

Before:

> Coordinates are $(t,x,y,z)$ with $t$ carrying dimensions of time and $x$, $y$ and $z$ the three lengths of a slice of constant $t$, which is ordinary flat space.
> The drive is written in the 3+1 form: the lapse is one, the slices are flat, and the whole of the geometry is a flow field $\vec{V} = u\,\partial_x + v\,\partial_y + w\,\partial_z$ whose components are velocities and which enters the line element as the shift.
> The observers who ride the slices, the Eulerian observers with $n_\mu = -\partial_\mu(ct)$, are carried through the coordinates at $\vec{V}$ and age at the rate $t$ counts, so a ship held at rest in the flow is in free fall and feels nothing.

After:

> We use coordinates $(t,x,y,z)$, with $t$ carrying dimensions of time and $x$, $y$, and $z$ the three lengths of a slice of constant $t$, which is ordinary flat space.
> We write the drive in the 3+1 form: the lapse is one, the slices are flat, and the geometry is wholly a flow field $\vec{V} = u\,\partial_x + v\,\partial_y + w\,\partial_z$, whose components are velocities and which enters the line element as the shift.
> The observers who ride the slices, the Eulerian observers with $n_\mu = -\partial_\mu(ct)$, are carried through the coordinates at $\vec{V}$ and age at the rate $t$ counts, so a ship held at rest in the flow is in free fall and feels nothing.

### conventions: `natario/convention[1]`

Source: `MFS/assets/data/metrics/natario.json`, `convention`, paragraph 2.

Before:

> Natário's condition on the field is that it be divergence free, $\partial_x u + \partial_y v + \partial_z w = 0$; that quantity is the expansion of those observers, so the condition says exactly that volume elements are carried along unchanged.
> Alcubierre's drive is the same line element with the one component flow $\vec{V} = c\,v_s f\,\partial_x$, whose divergence is $c\,v_s\,\partial_x f$ and is not zero, and that one number is the whole difference between the two spacetimes.
> It is a smaller difference than the story suggests.
> Putting a flow with one component into the energy density of the Eulerian observers gives $\theta^2 - \theta_{ij}\theta^{ij} = -\left(\left(\partial_y u\right)^2 + \left(\partial_z u\right)^2\right)/2$, in which the expansion has cancelled identically against the part of the rate of strain that always accompanies it, so the expansion never entered the energy density of Alcubierre's drive either and the two drives carry the same expression for it.

After:

> José Natário asks that the field be divergence free, $\partial_x u + \partial_y v + \partial_z w = 0$; that quantity is the expansion of those observers, so the condition says that volume elements are carried along unchanged.
> Alcubierre's drive is the same line element with the one component flow $\vec{V} = c\,v_s f\,\partial_x$, whose divergence is $c\,v_s\,\partial_x f$ and is not zero, and the two spacetimes differ in that one number alone.
> The difference does not reach the energy density.
> Putting a flow with one component into the energy density of the Eulerian observers gives $\theta^2 - \theta_{ij}\theta^{ij} = -\left(\left(\partial_y u\right)^2 + \left(\partial_z u\right)^2\right)/2$, in which the expansion has cancelled identically against the part of the rate of strain that always accompanies it, so the expansion never entered the energy density of Alcubierre's drive either, and the two drives carry the same expression for it.

### conventions: `natario/convention[2]`

Source: `MFS/assets/data/metrics/natario.json`, `convention`, paragraph 3.

Before:

> The Cartesian chart assumes nothing about the field, so every component in it holds for either drive and Natário's case is read off it by setting the divergence to zero.
> Every piece of the geometry is taken for the general flow: the metric and its inverse, the connection with its first index up and with every index down, the Riemann tensor in both of those positions, the Ricci tensor in all three index positions, the Ricci scalar, the Kretschmann scalar, the Einstein tensor in all three, the Weyl tensor in both positions and the geodesic equations.
> That generality costs length.
> Three free functions of four coordinates put 672 components into the Riemann and Weyl tensors, two index positions each, carrying 18738 terms between them, and the longest of those components runs to 5005 characters.
> The same geometry closes in a few lines in the plane symmetric chart, the simplest member of the family, and Alcubierre's drive is the member whose flow has one component.

After:

> The Cartesian chart assumes nothing about the field, so every component in it holds for either drive, and Natário's case follows from it on setting the divergence to zero.
> We take every piece of the geometry for the general flow: the metric and its inverse, the connection with its first index up and with every index down, the Riemann tensor in both of those positions, the Ricci tensor in all three index positions, the Ricci scalar, the Kretschmann scalar, the Einstein tensor in all three, the Weyl tensor in both positions, and the geodesic equations.
> The generality makes the components long.
> Three free functions of four coordinates put 672 components into the Riemann and Weyl tensors, two index positions each, carrying 18738 terms between them, and the longest of those components runs to 5005 characters.
> The same geometry closes in a few lines in the plane symmetric chart, the simplest member of the family, and Alcubierre's drive is the member whose flow has one component.

### conventions: `natario/convention[3]`

Source: `MFS/assets/data/metrics/natario.json`, `convention`, paragraph 4.

Before:

> The Weyl tensor of this chart is nowhere a copy of its Riemann tensor, since neither chart is a vacuum: the two are nonzero in the same 144 slots with every index lowered and disagree in all 144 of them.
> The Kretschmann scalar closes in one form only, the signed sum of squares of the three pieces of the 3+1 split, the Gauss, the Codazzi and the Ricci equation, taken in the orthonormal frame of the Eulerian observers where each piece is short; expanded flat it is 498 terms.
> The Einstein tensor with both indices up is where the content sits: $G^{tt}$ is the Hamiltonian constraint of a flat slice and $G^{tx}$, $G^{ty}$ and $G^{tz}$ are the momentum constraint, each of them carried along by the flow.
> The mixed Einstein tensor is the shortest of the three, because an upper time index is minus a projection on those observers, $g^{t\alpha} = -n^\alpha$, so $R^t{}_x$, $R^t{}_y$ and $R^t{}_z$ are the momentum constraint with the flow stripped off them, $\left(\partial_i\theta - \nabla^2V_i\right)/2c$; with both indices down it is the form the field equations are written in, and it is the longest, because a lower time index drags the whole spatial part along with the flow.

After:

> The Weyl tensor of this chart differs from its Riemann tensor everywhere, since neither chart is a vacuum: the two are nonzero in the same 144 slots with every index lowered and disagree in all 144 of them.
> The Kretschmann scalar closes in one form only, the signed sum of squares of the three pieces of the 3+1 split, the Gauss, the Codazzi, and the Ricci equation, taken in the orthonormal frame of the Eulerian observers, where each piece is short; expanded flat it is 498 terms.
> The Einstein tensor with both indices up holds the content: $G^{tt}$ is the Hamiltonian constraint of a flat slice, and $G^{tx}$, $G^{ty}$, and $G^{tz}$ are the momentum constraint, each of them carried along by the flow.
> The mixed Einstein tensor is the shortest of the three, because an upper time index is minus a projection on those observers, $g^{t\alpha} = -n^\alpha$, so $R^t{}_x$, $R^t{}_y$, and $R^t{}_z$ are the momentum constraint with the flow stripped off them, $\left(\partial_i\theta - \nabla^2V_i\right)/2c$; with both indices down, the form in which the field equations are written, it is the longest, because a lower time index drags the whole spatial part along with the flow.

### conventions: `natario/convention[4]`

Source: `MFS/assets/data/metrics/natario.json`, `convention`, paragraph 5.

Before:

> The plane symmetric chart is the member of the family whose flow points along $x$ and is carried by the two coordinates across it, $\vec{V} = u(t,y,z)\,\partial_x$.
> It is divergence free for every profile $u$, because a component that is free of its own coordinate contributes nothing to the divergence, so every component in it is a component of Natário's drive rather than of a general flow.
> It is the exact geometry of a warp corridor: a tube of space of any cross section sliding along $x$ at the speed $u$, flat inside where the ship rides, flat outside where the universe is at rest, and with the whole of the curvature in the wall between them.
> What such a flow cannot do is close, since being divergence free forbids it to depend on the coordinate it points along, so the corridor runs to infinity at both ends.
> Closing it off is exactly what calls for the transverse components that Natário's own field carries, and the Cartesian chart is where those live.

After:

> The plane symmetric chart is the member of the family whose flow points along $x$ and is carried by the two coordinates across it, $\vec{V} = u(t,y,z)\,\partial_x$.
> It is divergence free for every profile $u$, because a component that is free of its own coordinate contributes nothing to the divergence, so every component in it is a component of Natário's drive rather than of a general flow.
> It is the exact geometry of a warp corridor: a tube of space of any cross section sliding along $x$ at the speed $u$, flat inside where the ship rides, flat outside where the universe is at rest, and with the whole of the curvature in the wall between them.
> Such a flow cannot close, since being divergence free forbids it to depend on the coordinate it points along, so the corridor runs to infinity at both ends.
> Closing it off calls for the transverse components of Natário's own field, which live in the Cartesian chart.

### conventions: `natario/convention[5]`

Source: `MFS/assets/data/metrics/natario.json`, `convention`, paragraph 6.

Before:

> Components are taken in the chart $x^0 = ct$, so every partial derivative is taken with respect to that chart coordinate: $\partial_t u = c^{-1}\partial u/\partial t$, while $\partial_x$, $\partial_y$ and $\partial_z$ need no such factor.
> The dots in the geodesic equations are derivatives of the same chart with respect to an affine parameter, so $\dot{t}$ means $d(ct)/d\lambda$, which is what makes every term of every equation carry the dimensions of its left hand side.
> The combination $\dot{x} - c^{-1}u\,\dot{t}$ that runs through those equations is the velocity of a particle relative to the flow it is sitting in, and a geodesic is bent only through that combination and through the rate at which the flow itself changes.

After:

> We take components in the chart $x^0 = ct$, so every partial derivative is taken with respect to that chart coordinate: $\partial_t u = c^{-1}\partial u/\partial t$, while $\partial_x$, $\partial_y$, and $\partial_z$ need no such factor.
> A dot in the geodesic equations is a derivative of the same chart's coordinates with respect to an affine parameter, so $\dot{t}$ means $d(ct)/d\lambda$, and every term of every equation then carries the dimensions of its left hand side.
> The combination $\dot{x} - c^{-1}u\,\dot{t}$ that runs through those equations is the velocity of a particle relative to the flow it sits in, and a geodesic bends only through that combination and through the rate at which the flow itself changes.

### conventions: `natario/convention[6]`

Source: `MFS/assets/data/metrics/natario.json`, `convention`, paragraph 7.

Before:

> The Ricci tensor is the standard contraction $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$, so the field equations read $G_{\mu\nu} = 8\pi G\,T_{\mu\nu}/c^4$.
> Neither chart is a vacuum.
> Because $n_\mu = -\partial_\mu(ct)$, the energy density the Eulerian observers measure is the single component $\rho = T_{\mu\nu}n^\mu n^\nu = c^4G^{tt}/8\pi G$, and $G^{tt}$ is $\left(\theta^2 - \theta_{ij}\theta^{ij}\right)/2c^2$, with $\theta_{ij} = \left(\partial_iV_j + \partial_jV_i\right)/2$ the rate of strain of the flow and $\theta$ its trace.
> Natário's condition kills the trace and leaves minus a sum of squares, so $\rho = -c^2\theta_{ij}\theta^{ij}/16\pi G$, which is negative wherever the flow is strained at all and zero only where the flow is a rigid motion of flat space.
> Taking the expansion away therefore does not pay off the debt in negative energy; it only removes the one term that could ever have paid it.
> In the plane symmetric chart that statement is the closed expression $\rho = -c^2\left(\left(\partial_y u\right)^2 + \left(\partial_z u\right)^2\right)/32\pi G$.

After:

> The Ricci tensor is the standard contraction, $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$, so the field equations read $G_{\mu\nu} = 8\pi G\,T_{\mu\nu}/c^4$.
> Neither chart is a vacuum.
> Because $n_\mu = -\partial_\mu(ct)$, the energy density measured by the Eulerian observers is the single component $\rho = T_{\mu\nu}n^\mu n^\nu = c^4G^{tt}/8\pi G$, and $G^{tt}$ is $\left(\theta^2 - \theta_{ij}\theta^{ij}\right)/2c^2$, with $\theta_{ij} = \left(\partial_iV_j + \partial_jV_i\right)/2$ the rate of strain of the flow and $\theta$ its trace.
> Natário's condition kills the trace and leaves minus a sum of squares, so $\rho = -c^2\theta_{ij}\theta^{ij}/16\pi G$, which is negative wherever the flow is strained and zero only where the flow is a rigid motion of flat space.
> Taking the expansion away therefore leaves the negative energy in place and removes the one positive term, the square of the expansion.
> In the plane symmetric chart the density has the closed form $\rho = -c^2\left(\left(\partial_y u\right)^2 + \left(\partial_z u\right)^2\right)/32\pi G$.

### spacetime diagram caption: `diagrams/natario/systems/cartesian_flow/tx.caption[0]`

Source: `_tools/derivations/null_rays.py` line 879, `CAPTIONS ("natario", "cartesian_flow", "tx")`.

Before:

> This is the plane of $t$ and $x$ on the axis of motion, $y = z = 0$.
> There the field reduces to $u = 2nv_s = v_s f$, which is exactly Alcubierre's shift, so on this plane the two metrics are the same and so are their light rays.

After:

> This is the plane of $t$ and $x$ on the axis of motion, $y = z = 0$.
> There the field reduces to $u = 2nv_s = v_s f$, which is Alcubierre's shift, so on this plane the two metrics are the same, and so are their light rays.

### embedding diagram caption: `embedding/natario/plane.caption[1]`

Source: `_tools/derivations/embedding.py` line 3692, `CAPTIONS ("natario", "plane")`.

Before:

> The lines marked are lines of that flow.
> Inside the bubble space moves forward at $v_s$ with the ship, outside it is at rest, and every line closes back through the wall, as the flow of a fluid that cannot be compressed closes round an obstacle.
> José Natário built the drive in 2002 to show that Alcubierre's contraction ahead and expansion behind are not what carries the ship; the energy density the riding observers measure is still negative, $-c^4K_{ij}K^{ij}/16\pi G$, since the trace $K$ vanishes with the expansion.

After:

> The lines marked are lines of that flow.
> Inside the bubble space moves forward at $v_s$ with the ship, outside it is at rest, and every line closes back through the wall, as the flow of a fluid that cannot be compressed closes round an obstacle.
> José Natário built the drive in 2002 to show that the ship is carried without Alcubierre's contraction ahead and expansion behind; the energy density measured by the riding observers is still negative, $-c^4K_{ij}K^{ij}/16\pi G$, since the trace $K$ vanishes with the expansion.

### embedding diagram note: `embedding/natario/plane.input`

Source: `_tools/derivations/embedding.py` line 2917, `in natario() input=`.

Before:

> $v_s = 2$, $n = f/2$ with Alcubierre's profile, and the zero expansion field $X = v_s[(2n + \rho n')\,e_x - n'\,x_r\,(x_r, y, z)/\rho]$, $x_r = x - v_s t$, as the spacetime diagram declares.

After:

> $v_s = 2$, $n = f/2$ with Alcubierre's profile, and the zero expansion field $X = v_s[(2n + \rho n')\,e_x - n'\,x_r\,(x_r, y, z)/\rho]$, $x_r = x - v_s t$, as in the spacetime diagram.

## Oppenheimer-Snyder Gravitational Collapse

### history: `oppenheimer_snyder/history[0]`

Source: `MFS/assets/data/metrics/oppenheimer_snyder.json`, `history`, paragraph 1.

Before:

> Earlier in 1939, Oppenheimer and Volkoff had shown that a cold star above a certain mass has nothing left to hold it up [oppenheimer1939].
> What happens to it then?
> Oppenheimer answered the obvious next question that same year with his student Hartland Snyder, in a paper now recognized as the first description of a black hole forming in real time [oppenheimersnyder1939].
> They took the simplest collapse imaginable: a uniform sphere of pressureless dust, with no internal pressure to slow its fall, and solved it exactly by matching a contracting interior, itself a piece of a closed dust cosmology, onto the static Schwarzschild geometry outside [schwarzschild1916].

After:

> Earlier in 1939, J. Robert Oppenheimer and George Volkoff had shown that a cold star above a certain mass has nothing left to hold it up [oppenheimer1939].
> What happens to it then?
> Oppenheimer answered the obvious next question that same year with his student Hartland Snyder, in a paper now recognized as the first description of a black hole forming in real time [oppenheimersnyder1939].
> They took the simplest collapse imaginable: a uniform sphere of pressureless dust, with no internal pressure to slow its fall, and solved it exactly by matching a contracting interior, itself a piece of a closed dust cosmology, onto the static Schwarzschild geometry outside [schwarzschild1916].

### history: `oppenheimer_snyder/history[2]`

Source: `MFS/assets/data/metrics/oppenheimer_snyder.json`, `history`, paragraph 3.

Before:

> The paper was published on the day Germany invaded Poland, and it went largely unregarded for two decades, set aside by the war, and by a lingering disbelief that nature would permit such an object.
> Christopher Nolan put the coincidence on film in Oppenheimer, where Oppenheimer tells Snyder that their paper is in print and Snyder holds up a newspaper announcing that Hitler has invaded Poland [nolan2023oppenheimer].
> The paper's vindication came in the 1960s.
> The worry had always been that the singularity was an artifact of the model's perfect symmetry, an idealization realistic, lumpy collapse would avoid; Roger Penrose's 1965 theorem dispelled that hope, proving that once a trapped surface forms a singularity is unavoidable, symmetry or no [penrose1965].
> Oppenheimer and Snyder had not described a special case.
> They had described the generic fate of massive matter.

After:

> The paper was published on the day Germany invaded Poland, and it went largely unregarded for two decades, set aside by the war, and by a lingering disbelief that nature would permit such an object.
> Christopher Nolan put the coincidence on film in Oppenheimer, where Oppenheimer tells Snyder that their paper is in print and Snyder holds up a newspaper announcing that Hitler has invaded Poland [nolan2023oppenheimer].
> The paper's vindication came in the 1960s.
> The worry had always been that the singularity was an artifact of the model's perfect symmetry, an idealization realistic, lumpy collapse would avoid; Roger Penrose's 1965 theorem dispelled that hope, proving that once a trapped surface forms a singularity is unavoidable, symmetry or no [penrose1965].
> Oppenheimer and Snyder had described the generic fate of massive matter.

### history: `oppenheimer_snyder/history[4]`

Source: `MFS/assets/data/metrics/oppenheimer_snyder.json`, `history`, paragraph 5.

Before:

> The correspondence is suggestive rather than exact: the real star shed turbulent outer layers as it went, nothing like the pressureless, perfectly symmetric ball of the 1939 model, and the interpretation is contested.
> Emma Beasor and colleagues, pointing to a mid-infrared source that lingers a decade on and the absence of X-rays from accretion, argue the event was a dusty stellar eruption rather than a black hole's birth [beasor2026].
> Whichever way it resolves, the question Oppenheimer and Snyder opened, whether massive matter truly collapses without limit, is at last one that astronomers can test with telescopes.

After:

> The correspondence is suggestive rather than exact: the real star shed turbulent outer layers as it went, nothing like the pressureless, perfectly symmetric ball of the 1939 model, and the interpretation is contested.
> Emma Beasor and colleagues, pointing to a source in the middle infrared that lingers a decade on and the absence of X-rays from accretion, argue the event was a dusty stellar eruption rather than a black hole's birth [beasor2026].
> Whichever way it resolves, the question Oppenheimer and Snyder opened, whether massive matter collapses without limit, is at last one that astronomers can test with telescopes.

### conventions: `oppenheimer_snyder/convention[0]`

Source: `MFS/assets/data/metrics/oppenheimer_snyder.json`, `convention`, paragraph 1.

Before:

> The spacetime is two charts joined at the surface of the star.
> The interior is a closed dust cosmology in comoving coordinates $(\tau, \chi, \theta, \phi)$, where $\tau$ is the proper time of the dust, $\chi$ is a dimensionless comoving polar angle and the scale factor $a(\tau)$ carries the length, so the areal radius of the sphere at $\chi$ is $a\sin\chi$ and the star is the ball $\chi \le \chi_0$.
> It is the $k = +1$ case of the Friedmann-Lemaître-Robertson-Walker metric, written with the curvature in the coordinate rather than in a parameter.

After:

> We cover the spacetime with two charts joined at the surface of the star.
> The interior is a closed dust cosmology in comoving coordinates $(\tau, \chi, \theta, \phi)$, where $\tau$ is the proper time of the dust, $\chi$ is a dimensionless comoving polar angle, and the scale factor $a(\tau)$ carries the length, so the areal radius of the sphere at $\chi$ is $a\sin\chi$, and the star is the ball $\chi \le \chi_0$.
> It is the $k = +1$ case of the Friedmann-Lemaître-Robertson-Walker metric, written with the curvature in the coordinate rather than in a parameter.

### conventions: `oppenheimer_snyder/convention[1]`

Source: `MFS/assets/data/metrics/oppenheimer_snyder.json`, `convention`, paragraph 2.

Before:

> The exterior is the Schwarzschild vacuum in the standard spherical chart, with $r_s = 2GM/c^2$.
> It is valid from the moment of rest on and outside the surface, $t \ge 0$ and $r \ge R(t)$, where $R = a\sin\chi_0$ is the areal radius of the surface read along the exterior time and $t = 0$ is set at the moment of rest, when the star and the static exterior share one slice of time symmetry.
> The model begins there, so before it there is no surface to be outside of.
> $R(t)$ stays above $r_s$ at every finite $t$, since the surface takes an infinite Schwarzschild time to reach $r_s$, so the part of the black hole outside the star lies beyond both sets of coordinates.

After:

> The exterior is the Schwarzschild vacuum in the standard spherical chart, with $r_s = 2GM/c^2$.
> It is valid from the moment of rest on and outside the surface, $t \ge 0$ and $r \ge R(t)$, where $R = a\sin\chi_0$ is the areal radius of the surface taken along the exterior time, and we set $t = 0$ at the moment of rest, when the star and the static exterior share one slice of time symmetry.
> The model begins there, so before it there is no surface to be outside of.
> The areal radius $R(t)$ stays above $r_s$ at every finite $t$, since the surface takes an infinite Schwarzschild time to reach $r_s$, so the part of the black hole outside the star lies beyond both sets of coordinates.

### conventions: `oppenheimer_snyder/convention[2]`

Source: `MFS/assets/data/metrics/oppenheimer_snyder.json`, `convention`, paragraph 3.

Before:

> Factors of $c$ and $G$ are kept explicit.
> Components are taken in the chart $x^0 = c\tau$ inside and $x^0 = ct$ outside, so a chart component of either metric carries no dimensions, and a dot on the scale factor is a derivative with respect to that chart time, $\dot{a} = \dfrac{da}{d(c\tau)} = \dfrac{1}{c}\dfrac{da}{d\tau}$, which is dimensionless.
> The dots in the geodesic equations are velocities of the same chart with respect to an affine parameter, so $\dot{\tau}$ means $d(c\tau)/d\lambda$ and $\dot{t}$ means $d(ct)/d\lambda$, and the equations are $\ddot{x}^\mu + \Gamma^\mu{}_{\nu\rho}\dot{x}^\nu\dot{x}^\rho = 0$ with the Christoffel symbols of that chart.

After:

> We keep factors of $c$ and $G$ explicit.
> We take components in the chart $x^0 = c\tau$ inside and $x^0 = ct$ outside, so a chart component of either metric carries no dimensions, and a dot on the scale factor is a derivative with respect to that chart time, $\dot{a} = \dfrac{da}{d(c\tau)} = \dfrac{1}{c}\dfrac{da}{d\tau}$, which is dimensionless.
> The dots in the geodesic equations are velocities in the same chart with respect to an affine parameter, so $\dot{\tau}$ means $d(c\tau)/d\lambda$ and $\dot{t}$ means $d(ct)/d\lambda$, and the equations read $\ddot{x}^\mu + \Gamma^\mu{}_{\nu\rho}\dot{x}^\nu\dot{x}^\rho = 0$, with the Christoffel symbols of that chart.

### conventions: `oppenheimer_snyder/convention[3]`

Source: `MFS/assets/data/metrics/oppenheimer_snyder.json`, `convention`, paragraph 4.

Before:

> The Ricci tensor is the standard contraction $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$, so the interior has $G_{\tau\tau} = 3\dfrac{\dot{a}^2+1}{a^2} = \dfrac{8\pi G\rho}{c^2}$, the Friedmann equation with a positive energy density on the other side of it, and the pressureless condition $G_{\chi\chi} = 0$ is $2a\ddot{a}+\dot{a}^2+1 = 0$.
> That equation has the first integral $\dot{a}^2 + 1 = \dfrac{a_m}{a}$, whose solution is the cycloid $a = \dfrac{a_m}{2}(1+\cos\eta)$, $c\tau = \dfrac{a_m}{2}(\eta+\sin\eta)$, collapsing from rest at $a = a_m$ to $a = 0$ in the proper time $\tau_s = \dfrac{\pi a_m}{2c}$.
> The components leave $a(\tau)$ a free function, exactly as those of the Friedmann-Lemaître-Robertson-Walker metric do, so that every component is an identity rather than a claim about the dust solution alone.

After:

> The Ricci tensor is the standard contraction, $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$, so the interior has $G_{\tau\tau} = 3\dfrac{\dot{a}^2+1}{a^2} = \dfrac{8\pi G\rho}{c^2}$, the Friedmann equation with a positive energy density on the other side of it, and the pressureless condition $G_{\chi\chi} = 0$ is $2a\ddot{a}+\dot{a}^2+1 = 0$.
> That equation has the first integral $\dot{a}^2 + 1 = \dfrac{a_m}{a}$, whose solution is the cycloid $a = \dfrac{a_m}{2}(1+\cos\eta)$, $c\tau = \dfrac{a_m}{2}(\eta+\sin\eta)$, collapsing from rest at $a = a_m$ to $a = 0$ in the proper time $\tau_s = \dfrac{\pi a_m}{2c}$.
> We leave $a(\tau)$ a free function in the components, as the components of the Friedmann-Lemaître-Robertson-Walker metric do, so every component is an identity rather than a claim about the dust solution alone.

### conventions: `oppenheimer_snyder/convention[4]`

Source: `MFS/assets/data/metrics/oppenheimer_snyder.json`, `convention`, paragraph 5.

Before:

> The two charts are joined by the Darmois conditions, that the induced metric and the extrinsic curvature agree across the surface.
> The first gives $R = a\sin\chi_0$ and makes $\tau$ the proper time of the surface on both sides; the second gives $K_{\tau\tau} = 0$, which is the statement that the surface falls freely, and $\cos\chi_0 = \sqrt{1-\dfrac{r_s}{R}+\left(\dfrac{dR}{d(c\tau)}\right)^2}$, which together with the interior's first integral forces $r_s = a_m\sin^3\chi_0$.
> That single equation is the whole of the matching.
> It says $M = \dfrac{4\pi}{3}\rho R^3$ at every moment of the collapse, it reads $\chi_0$ as a compactness through $\sin^2\chi_0 = r_s/R_0$ and as a redshift through $\cos\chi_0 = \sqrt{1-r_s/R_0}$, and it puts the surface on the radial Schwarzschild geodesic of energy $\tilde{E} = c^2\cos\chi_0$.

After:

> We join the two charts by the conditions of Georges Darmois, that the induced metric and the extrinsic curvature agree across the surface.
> The first gives $R = a\sin\chi_0$ and makes $\tau$ the proper time of the surface on both sides; the second gives $K_{\tau\tau} = 0$, which is the statement that the surface falls freely, and $\cos\chi_0 = \sqrt{1-\dfrac{r_s}{R}+\left(\dfrac{dR}{d(c\tau)}\right)^2}$, which together with the interior's first integral forces $r_s = a_m\sin^3\chi_0$.
> The matching comes down to that single equation.
> It says $M = \dfrac{4\pi}{3}\rho R^3$ at every moment of the collapse, it makes $\chi_0$ a compactness, through $\sin^2\chi_0 = r_s/R_0$, and a redshift, through $\cos\chi_0 = \sqrt{1-r_s/R_0}$, and it puts the surface on the radial Schwarzschild geodesic of energy $\tilde{E} = c^2\cos\chi_0$.

### conventions: `oppenheimer_snyder/convention[5]`

Source: `MFS/assets/data/metrics/oppenheimer_snyder.json`, `convention`, paragraph 6.

Before:

> The surface reaches $r = r_s$ at $\eta = \pi - 2\chi_0$, a finite proper time before the singularity, and takes an infinite Schwarzschild time $t$ to get there, which is the same crossing seen from the two sides.
> The interior is conformally flat, so its Weyl tensor vanishes identically and all of its curvature is the dust; the exterior is Ricci flat, so all of its curvature is Weyl.
> Neither tensor is continuous at the surface, which the junction conditions do not ask for: the Kretschmann scalar is $15r_s^2/R^6$ just inside and $12r_s^2/R^6$ just outside.

After:

> The surface reaches $r = r_s$ at $\eta = \pi - 2\chi_0$, a finite proper time before the singularity, and takes an infinite Schwarzschild time $t$ to get there; the two times belong to one crossing, seen from the two sides.
> The interior is conformally flat, so its Weyl tensor vanishes identically and all of its curvature is the dust; the exterior is Ricci flat, so all of its curvature is Weyl.
> Neither tensor is continuous at the surface, and the junction conditions do not require it to be: the Kretschmann scalar is $15r_s^2/R^6$ just inside and $12r_s^2/R^6$ just outside.

### spacetime diagram caption: `diagrams/oppenheimer_snyder/systems/exterior_schwarzschild/radial.caption[0]`

Source: `_tools/derivations/null_rays.py` line 1009, `CAPTIONS ("oppenheimer_snyder", "exterior_schwarzschild", "radial")`.

Before:

> This is the plane of $t$ and $r$ at $\theta = \pi/2$ and $\phi = 0$ outside the collapsing star, where the metric is Schwarzschild's.
> The surface falls freely from rest at $r = 2r_s$ at $t = 0$, along its radial geodesic, and what lies inside it, $r < R(t)$, is the star, which these coordinates do not cover.

After:

> This is the plane of $t$ and $r$ at $\theta = \pi/2$ and $\phi = 0$ outside the collapsing star, where the metric is Schwarzschild's.
> The surface falls freely from rest at $r = 2r_s$ at $t = 0$, along its radial geodesic, and inside it, $r < R(t)$, lies the star, which these coordinates do not cover.

## pp-wave Plane Gravitational Waves

### history: `pp_wave/history[0]`

Source: `MFS/assets/data/metrics/pp_wave.json`, `history`, paragraph 1.

Before:

> H. W. Brinkmann was not looking for a gravitational wave.
> In 1925, asking which Einstein spaces map conformally onto one another, he was handed a family carrying a null vector field parallel to itself everywhere [brinkmann1925].
> He did not call them waves.
> Nobody did for thirty six years, until Kundt and Ehlers named them the plane-fronted waves with parallel rays [kundt1961, ehlerskundt1962].
> Relativists were busy with something more basic: Eddington sorted the plane waves into three kinds in 1922 and found two of them travelling at no fixed speed at all, mere sinuosities in the coordinate system whose only relevant speed, he wrote, was the speed of thought [eddington1922].
> The remedy, to follow the curvature and not the potentials, he left to somebody else.

After:

> Hans Brinkmann was not looking for a gravitational wave.
> In 1925, asking which Einstein spaces map conformally onto one another, he was handed a family carrying a null vector field parallel to itself everywhere [brinkmann1925].
> He did not call them waves.
> Nobody did for thirty six years, until Wolfgang Kundt and Jürgen Ehlers named them the plane-fronted waves with parallel rays [kundt1961, ehlerskundt1962].
> Relativists were busy with something more basic: Arthur Eddington sorted the plane waves into three kinds in 1922 and found two of them travelling at no fixed speed, mere sinuosities in the coordinate system whose only relevant speed, he wrote, was the speed of thought [eddington1922].
> The remedy, to follow the curvature and not the potentials, he left to somebody else.

### history: `pp_wave/history[4]`

Source: `MFS/assets/data/metrics/pp_wave.json`, `history`, paragraph 5.

Before:

> At Chapel Hill in January 1957 Bondi asked from the floor whether it could be made into an absorber of gravitational energy.
> Pirani had put in no absorption term, only a spring.
> Feynman supplied the friction: let one mass carry a stick running past touching a second, the wave slides them against each other, the rubbing makes heat, and heat is energy [feynman2011].
> He had hesitated to say it, having missed the session on gravitational waves entirely.

After:

> At Chapel Hill in January 1957 Hermann Bondi asked from the floor whether it could be made into an absorber of gravitational energy.
> Pirani had put in no absorption term, only a spring.
> Richard Feynman supplied the friction: let one mass carry a stick running past touching a second, the wave slides them against each other, the rubbing makes heat, and heat is energy [feynman2011].
> He had hesitated to say it, having missed the session on gravitational waves entirely.

### history: `pp_wave/history[5]`

Source: `MFS/assets/data/metrics/pp_wave.json`, `history`, paragraph 6.

Before:

> Bondi came a skeptic and left converted, with a version of the argument in Nature within months [bondi1957].
> In 1959, with Pirani and Ivor Robinson, he gave the exact plane waves, the most symmetric of Brinkmann's family, their full treatment [bondipirani1959].
> Then Penrose found the sting: plane waves are geodesically complete and about as well behaved as a spacetime gets, yet they focus light so relentlessly that no spacelike surface in one can carry Cauchy data for the whole of it [penrose1965plane].

After:

> Bondi came a skeptic and left converted, with a version of the argument in Nature within months [bondi1957].
> In 1959, with Pirani and Ivor Robinson, he gave the exact plane waves, the most symmetric of Brinkmann's family, their full treatment [bondipirani1959].
> Then Roger Penrose found the sting: plane waves are geodesically complete and about as well behaved as a spacetime gets, yet they focus light so relentlessly that no spacelike surface in one can carry Cauchy data for the whole of it [penrose1965plane].

### history: `pp_wave/history[6]`

Source: `MFS/assets/data/metrics/pp_wave.json`, `history`, paragraph 7.

Before:

> Joseph Weber, also in that audience, went home to Maryland and built the first antenna out of aluminum [weber1960].
> His own detections did not survive scrutiny; the enterprise did.
> Hulse and Taylor's binary pulsar was losing orbital period at exactly the radiated rate [hulse1975, taylor1979], and on 14 September 2015 two detectors caught the same fifth of a second of chirp [abbott2016].
> LIGO is Feynman's stick at scale, with hanging mirrors for the rubbing masses and a laser where the friction was, and what sweeps through it is, to any accuracy the instrument can tell, one of Brinkmann's metrics.

After:

> Joseph Weber, also in that audience, went home to Maryland and built the first antenna out of aluminum [weber1960].
> His own detections did not survive scrutiny; the enterprise did.
> Russell Hulse and Joseph Taylor's binary pulsar was losing orbital period at exactly the radiated rate [hulse1975, taylor1979], and on 14 September 2015 two detectors caught the same fifth of a second of chirp [abbott2016].
> LIGO is Feynman's stick at scale, with hanging mirrors for the rubbing masses and a laser where the friction was, and the wave sweeping through it is, to any accuracy the instrument can tell, one of Brinkmann's metrics.

### conventions: `pp_wave/convention[0]`

Source: `MFS/assets/data/metrics/pp_wave.json`, `convention`, paragraph 1.

Before:

> Coordinates are $(u, v, x, y)$ with $u$ the retarded time labelling the wave fronts and carrying dimensions of time, $v$ an affine parameter along the null rays and carrying a length, and $x$ and $y$ the two transverse lengths.
> The profile $H$ is an arbitrary dimensionless function of $u$, $x$ and $y$; that $v$ is absent from it is what leaves $\partial_v$ covariantly constant and null, which is the parallel ray the family is named for.
> Components are taken in the chart $x^0 = cu$, so every partial derivative is taken with respect to that chart coordinate, $\partial_u H = c^{-1}\partial H/\partial u$, and a prime on an amplitude is the same derivative; the transverse $\partial_x$ and $\partial_y$ need no such factor.
> The dots in the geodesic equations are derivatives of that same chart with respect to an affine parameter, so $\dot{u}$ means $d(cu)/d\lambda$, which is what makes every term of every equation carry the dimensions of its left hand side.
> The Ricci tensor is contracted on the middle index, $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$, which is the standard contraction and the one every Ricci and Einstein component is in.

After:

> We use coordinates $(u, v, x, y)$, with $u$ the retarded time labelling the wave fronts and carrying dimensions of time, $v$ an affine parameter along the null rays and carrying a length, and $x$ and $y$ the two transverse lengths.
> The profile $H$ is an arbitrary dimensionless function of $u$, $x$, and $y$; since $v$ is absent from it, $\partial_v$ is covariantly constant and null, the parallel ray for which the family is named.
> We take components in the chart $x^0 = cu$, so every partial derivative is taken with respect to that chart coordinate, $\partial_u H = c^{-1}\partial H/\partial u$, and a prime on an amplitude is the same derivative; the transverse $\partial_x$ and $\partial_y$ need no such factor.
> A dot in the geodesic equations is a derivative of that chart's coordinates with respect to an affine parameter, so $\dot{u}$ means $d(cu)/d\lambda$, and every term of every equation then carries the dimensions of its left hand side.
> We contract the Ricci tensor on the middle index, $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$, the standard contraction, and every Ricci and Einstein component uses it.

### conventions: `pp_wave/convention[1]`

Source: `MFS/assets/data/metrics/pp_wave.json`, `convention`, paragraph 2.

Before:

> In the Brinkmann chart the profile is arbitrary, so its Ricci and Einstein tensors are not assumed to vanish: the one surviving component of each is $R_{uu} = G_{uu} = -\dfrac{1}{2}\left(\partial_x^2 H + \partial_y^2 H\right)$, and the field equations reduce to the single condition that $H$ be harmonic in the transverse plane, $\partial_x^2 H + \partial_y^2 H = 0$, at every moment of retarded time.
> The exact plane wave chart is the subcase whose profile is quadratic in $x$ and $y$, $H = A(u)(x^2 - y^2) + 2B(u)xy$, which is the general quadratic solution of that condition and so is a vacuum for any amplitudes.
> The Weyl tensor carries the trace free part of the transverse Hessian of $H$ and is of Petrov type N, the algebraically special type of pure radiation, so it vanishes only where those second derivatives are isotropic, $\partial_x^2 H = \partial_y^2 H$ and $\partial_x\partial_y H = 0$, and is nonzero for every other profile.
> Both the Ricci scalar and the Kretschmann scalar vanish identically, for every profile and whether or not the spacetime is a vacuum, because raising a $u$ index turns it into a $v$ index and no curvature component of this family carries one, so every full contraction meets a zero.

After:

> In Hans Brinkmann's chart the profile is arbitrary, so we do not assume its Ricci and Einstein tensors vanish: the one surviving component of each is $R_{uu} = G_{uu} = -\dfrac{1}{2}\left(\partial_x^2 H + \partial_y^2 H\right)$, and the field equations reduce to the single condition that $H$ be harmonic in the transverse plane, $\partial_x^2 H + \partial_y^2 H = 0$, at every moment of retarded time.
> The exact plane wave chart is the subcase whose profile is quadratic in $x$ and $y$, $H = A(u)(x^2 - y^2) + 2B(u)xy$, which is the general quadratic solution of that condition and so is a vacuum for any amplitudes.
> The Weyl tensor carries the trace free part of the transverse Hessian of $H$ and is of Petrov type N, the algebraically special type of pure radiation, so it vanishes only where those second derivatives are isotropic, $\partial_x^2 H = \partial_y^2 H$ and $\partial_x\partial_y H = 0$, and is nonzero for every other profile.
> Both the Ricci scalar and the Kretschmann scalar vanish identically, for every profile and whether or not the spacetime is a vacuum, because raising a $u$ index turns it into a $v$ index and no curvature component of this family carries one, so every full contraction meets a zero.

### spacetime diagram caption: `diagrams/pp_wave/systems/exact_plane_wave/tz.caption[0]`

Source: `_tools/derivations/null_rays.py` line 897, `CAPTIONS ("pp_wave", "exact_plane_wave", "tz")`.

Before:

> This is the plane the wave travels in, on its axis $x = y = 0$, drawn with $u = t - z$ and $v = (t + z)/2$ so that the axes read as $t$ and $z$; the chart's own $u$ and $v$ are both null.
> On the axis the profile $A(x^2 - y^2) + 2Bxy$ vanishes whatever $A$ and $B$ are, so the metric on this plane is flat and the rays are at 45°.

After:

> This is the plane the wave travels in, on its axis $x = y = 0$, drawn with $u = t - z$ and $v = (t + z)/2$ so that the axes stand for $t$ and $z$; the chart's own $u$ and $v$ are both null.
> On the axis the profile $A(x^2 - y^2) + 2Bxy$ vanishes whatever $A$ and $B$ are, so the metric on this plane is flat and the rays are at 45°.

### spacetime diagram note: `diagrams/pp_wave/systems/exact_plane_wave/tz.input`

Source: `_tools/derivations/null_rays.py` line 450, `DIAGRAMS input=`.

Before:

> The amplitudes $A(u)$ and $B(u)$ do not enter on the axis, so any profile gives this diagram.

After:

> The amplitudes $A(u)$ and $B(u)$ do not enter on the axis, so the diagram is the same for every profile.

### embedding diagram caption: `embedding/pp_wave/ring.caption[0]`

Source: `_tools/derivations/embedding.py` line 3645, `CAPTIONS ("pp_wave", "ring")`.

Before:

> This is the wave front of a plane gravitational wave at four values of its retarded time $u$, each drawn as a surface in flat space so that every distance along it is the metric distance.
> A surface of constant $u$ has the metric $dx^2 + dy^2$ whatever $v$ is on it, so the drawing is a flat disc, and what the wave does shows in a ring of free particles at rest on the circle $x^2 + y^2 = L^2$ before the pulse arrives.

After:

> This is the wave front of a plane gravitational wave at four values of its retarded time $u$, each drawn as a surface in flat space so that every distance along it is the metric distance.
> A surface of constant $u$ has the metric $dx^2 + dy^2$ whatever $v$ is on it, so the drawing is a flat disc, and the wave shows in a ring of free particles at rest on the circle $x^2 + y^2 = L^2$ before the pulse arrives.

### embedding diagram note: `embedding/pp_wave/ring.input`

Source: `_tools/derivations/embedding.py` line 3082, `in pp_wave() input=`.

Before:

> A pulse of the plus polarisation, $A = e^{-u^2}/L^2$ and $B = 0$, as the spacetime diagram declares.

After:

> A pulse of the plus polarisation, $A = e^{-u^2}/L^2$ and $B = 0$, as in the spacetime diagram.

## Reissner-Nordström Charged Black Hole

### history: `rn_metric/history[0]`

Source: `MFS/assets/data/metrics/rn_metric.json`, `history`, paragraph 1.

Before:

> The charged generalization of the Schwarzschild metric arrived in pieces, and the story parallels the uncharged case almost beat for beat.
> In May 1916, Hans Reissner gave the gravitational field of a point charge, solving the Einstein equations in the presence of a Coulombic electric field [reissner1916].
> His derivation used Schwarzschild's auxiliary radius $R$, so the solution was not yet in its modern form, and the structure of its singularities was interpreted through the same lens that had obscured Schwarzschild's.
> Nordström was himself no stranger to the problem of gravity: his scalar theory of gravitation had been a serious contender before Einstein's field equations settled the matter in 1915.

After:

> The charged generalization of the Schwarzschild metric arrived in pieces, and the story parallels the uncharged case almost beat for beat.
> In May 1916, Hans Reissner gave the gravitational field of a point charge, solving the Einstein equations in the presence of a Coulombic electric field [reissner1916].
> His derivation used Schwarzschild's auxiliary radius $R$, so the solution was not yet in its modern form, and its singularities were interpreted in the same auxiliary radius that had hidden the structure of Schwarzschild's.
> Gunnar Nordström was himself no stranger to the problem of gravity: his scalar theory of gravitation had been a serious contender before Einstein's field equations settled the matter in 1915.

### history: `rn_metric/history[1]`

Source: `MFS/assets/data/metrics/rn_metric.json`, `history`, paragraph 2.

Before:

> He turned to the victorious theory anyway, and in 1918 did for Reissner's metric what Droste had done for Schwarzschild's: imposing isotropy directly on the auxiliary radius, he extended the areal coordinate to $r = 0$, giving the solution its modern form [nordstrom1918].
> Hermann Weyl had arrived at the same metric independently in 1917, as one member of a broad class of static fields in which the gravitational and electrostatic potentials depend on one another [weyl1917, bicak2000].
> George Barker Jeffery derived it again in 1921, solving the field equations directly for the field of an electron, and learned of Nordström's work only after his own paper was written [jeffery1921].

After:

> Nordström turned to the victorious theory anyway, and in 1918 did for Reissner's metric what Johannes Droste had done for Schwarzschild's: imposing isotropy directly on the auxiliary radius, he extended the areal coordinate to $r = 0$, giving the solution its modern form [nordstrom1918].
> Hermann Weyl had arrived at the same metric independently in 1917, as one member of a broad class of static fields in which the gravitational and electrostatic potentials depend on one another [weyl1917, bicak2000].
> George Barker Jeffery derived it again in 1921, solving the field equations directly for the field of an electron, and learned of Nordström's work only after his own paper was written [jeffery1921].

### history: `rn_metric/history[5]`

Source: `MFS/assets/data/metrics/rn_metric.json`, `history`, paragraph 6.

Before:

> Michael Simpson and Roger Penrose found in 1973, from computer calculations of test electromagnetic fields, that instabilities arise at the inner horizon though not at the event horizon, and inferred that in the full theory the inner horizon becomes a curvature singularity once perturbations without symmetry are present [simpson1973].
> Eric Poisson and Werner Israel made the mechanism explicit in 1989 and 1990: the slowly decaying radiative tail of a collapse is blueshifted without limit at the inner horizon and inflates the internal mass parameter of the hole without bound, while the mass seen from outside stays finite, which they called mass inflation [poisson1989, poisson1990].
> Mihalis Dafermos proved in 2005 that their picture is right for a spherical charged black hole coupled to a scalar field: the spacetime extends continuously across a boundary on which the Hawking mass blows up, so it cannot be extended as a differentiable metric, and in Christodoulou's formulation by continuous extension the strong cosmic censorship conjecture fails for that system [dafermos2005].

After:

> Michael Simpson and Roger Penrose found in 1973, from computer calculations of test electromagnetic fields, that instabilities arise at the inner horizon though not at the event horizon, and inferred that in the full theory the inner horizon becomes a curvature singularity once perturbations without symmetry are present [simpson1973].
> Eric Poisson and Werner Israel made the mechanism explicit in 1989 and 1990: the slowly decaying radiative tail of a collapse is blueshifted without limit at the inner horizon and inflates the internal mass parameter of the hole without bound, while the mass seen from outside stays finite, which they called mass inflation [poisson1989, poisson1990].
> Mihalis Dafermos proved in 2005 that their picture is right for a spherical charged black hole coupled to a scalar field: the spacetime extends continuously across a boundary on which the Hawking mass blows up, so it cannot be extended as a differentiable metric, and in Demetrios Christodoulou's formulation by continuous extension the strong cosmic censorship conjecture fails for that system [dafermos2005].

### conventions: `rn_metric/convention[0]`

Source: `MFS/assets/data/metrics/rn_metric.json`, `convention`, paragraph 1.

Before:

> The coordinate $r$ is the areal radius.
> The parameter $r_s = 2GM/c^2$ is the Schwarzschild radius and $r_q^2 = Q^2 G / (4\pi\epsilon_0 c^4)$ is the charge radius squared, where $Q$ is the total electric charge.
> The metric function is $\Delta = r^2 - r_s r + r_q^2$.
> Horizons occur at $r_\pm = \dfrac{r_s}{2} \pm \sqrt{\dfrac{r_s^2}{4} - r_q^2}$ when $r_s \geq 2r_q$.

After:

> The coordinate $r$ is the areal radius.
> The parameter $r_s = 2GM/c^2$ is the Schwarzschild radius, and $r_q^2 = Q^2 G / (4\pi\epsilon_0 c^4)$ is the square of the charge radius, where $Q$ is the total electric charge.
> We write the metric with $\Delta = r^2 - r_s r + r_q^2$.
> The horizons sit at $r_\pm = \dfrac{r_s}{2} \pm \sqrt{\dfrac{r_s^2}{4} - r_q^2}$ when $r_s \geq 2r_q$.

### conventions: `rn_metric/convention[1]`

Source: `MFS/assets/data/metrics/rn_metric.json`, `convention`, paragraph 2.

Before:

> When $r_s > 2r_q$ the coordinates cover three regions of the maximally extended spacetime, one of each kind: an exterior $r > r_+$, a region $r_- < r < r_+$ where $r$ is the time, and an inner region $0 < r < r_-$ where $t$ is a time again and the singularity at $r = 0$ is timelike.
> The metric is unchanged under $t \to -t$, so between the horizons the coordinates cannot say whether $r$ decreases toward the future, the black hole, or increases, the white hole.
> Here the region between the horizons is the black hole, the one an observer falling in from the exterior enters, and the inner region lies to its future.

After:

> When $r_s > 2r_q$ the coordinates cover three regions of the maximally extended spacetime, one of each kind: an exterior $r > r_+$, a region $r_- < r < r_+$ where $r$ is the time, and an inner region $0 < r < r_-$ where $t$ is a time again and the singularity at $r = 0$ is timelike.
> The metric is unchanged under $t \to -t$, so between the horizons the coordinates alone do not fix whether $r$ decreases toward the future, the black hole, or increases, the white hole.
> We take the region between the horizons to be the black hole, the one an observer falling in from the exterior enters, and the inner region lies to its future.

### conventions: `rn_metric/convention[2]`

Source: `MFS/assets/data/metrics/rn_metric.json`, `convention`, paragraph 3.

Before:

> The Ricci tensor is the standard contraction $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$, on which the Einstein tensor is that of the radial electric field: $G^t{}_t = G^r{}_r = -\dfrac{r_q^2}{r^4}$ carries a positive energy density and the two angular components carry the opposite sign.
> Units with $c = 1$ are used throughout, and components are taken in the chart $x^0 = ct$, which those units make $t$ itself, so $t$ is a length beside $r$, $r_s$ and $r_q$ and no component carries a factor of $c$.
> The geodesic equations are $\ddot{x}^\mu + \Gamma^\mu{}_{\nu\rho}\dot{x}^\nu\dot{x}^\rho = 0$ with $\Gamma$ the Christoffel symbols of that same chart, and the dots are derivatives of that chart with respect to an affine parameter $\lambda$, so $\dot{t}$ means $d(ct)/d\lambda$ and every term of every equation carries the dimensions of its left hand side.

After:

> With the Ricci tensor the standard contraction, $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$, the Einstein tensor is that of the radial electric field: $G^t{}_t = G^r{}_r = -\dfrac{r_q^2}{r^4}$ carries a positive energy density, and the two angular components carry the opposite sign.
> We work in units with $c = 1$ throughout and take components in the chart $x^0 = ct$, which is then $t$ itself, so $t$ is a length beside $r$, $r_s$, and $r_q$, and no component carries a factor of $c$.
> The geodesic equations read $\ddot{x}^\mu + \Gamma^\mu{}_{\nu\rho}\dot{x}^\nu\dot{x}^\rho = 0$, with $\Gamma$ the Christoffel symbols of the same chart, and a dot is a derivative of that chart's coordinates with respect to an affine parameter $\lambda$, so $\dot{t}$ means $d(ct)/d\lambda$, and every term of every equation carries the dimensions of its left hand side.

### spacetime diagram caption: `diagrams/rn_metric/systems/spherical/radial.caption[1]`

Source: `_tools/derivations/null_rays.py` line 698, `CAPTIONS ("rn_metric", "spherical", "radial")`.

Before:

> The chart alone cannot say which way is future in the two inner regions.
> An ingoing chart runs smoothly through both horizons, and taking the future from it makes the region between the horizons the black hole.
> The Kretschmann scalar diverges at $r = 0$.

After:

> The chart alone does not fix which way is future in the two inner regions.
> An ingoing chart runs smoothly through both horizons, and we take the future from it, which makes the region between the horizons the black hole.
> The Kretschmann scalar diverges at $r = 0$.

### conformal diagram caption: `conformal/rn_metric/tower.caption[1]`

Source: `_tools/derivations/conformal.py` line 2231, `CAPTIONS ("rn_metric", "tower")`.

Before:

> Every region is placed by the Kruskal coordinate of the outer horizon, $p = \pm\arctan e^{-\kappa_+ u}$ and $q = \pm\arctan e^{\kappa_+ v}$ with $u, v = t \mp r_*$, and the regions above the inner horizon are the reflection $(p, q) \to (\pi - q, \pi - p)$ of those below.
> This map is smooth across $r_+$ and puts the singularity $r = 0$ exactly on the vertical lines $X = \pm\pi/2$, where it is timelike.
> Across $r_-$ it is continuous and cannot also be smooth, because the late light rays that reach $\mathscr{I}^+$ are the rays that pile up at the Cauchy horizon $r_-$, and one function of the ray has to serve both.

After:

> We place every region by the Kruskal coordinate of the outer horizon, $p = \pm\arctan e^{-\kappa_+ u}$ and $q = \pm\arctan e^{\kappa_+ v}$ with $u, v = t \mp r_*$, and the regions above the inner horizon are the reflection $(p, q) \to (\pi - q, \pi - p)$ of those below.
> This map is smooth across $r_+$ and puts the singularity $r = 0$ exactly on the vertical lines $X = \pm\pi/2$, where it is timelike.
> Across $r_-$ it is continuous and cannot also be smooth, because the late light rays that reach $\mathscr{I}^+$ are the rays that pile up at the Cauchy horizon $r_-$, and one function of the ray has to serve both.

### conformal diagram caption: `conformal/rn_metric/tower.caption[2]`

Source: `_tools/derivations/conformal.py` line 2239, `CAPTIONS ("rn_metric", "tower")`.

Before:

> The coordinates $t$ and $r > 0$ cover one region of each kind: an exterior, a region between the horizons and a region inside $r_-$.
> Between the horizons their $t$ alone cannot tell the black hole from the white hole, and the region is taken to be the black hole an infalling observer enters.

After:

> The coordinates $t$ and $r > 0$ cover one region of each kind: an exterior, a region between the horizons, and a region inside $r_-$.
> Between the horizons their $t$ alone cannot tell the black hole from the white hole, and we take the region to be the black hole an infalling observer enters.

### conformal diagram note: `conformal/rn_metric/tower.settings`

Source: `_tools/derivations/conformal.py` line 1076, `in reissner_nordstrom()`.

Shown in: `conformal/rn_metric/tower.settings`, `conformal/rn_metric/malament_hogarth.settings`.

Before:

> $r_q = 0.48\,r_s$, so that $r_+ = 0.64\,r_s$, $r_- = 0.36\,r_s$ and $\kappa_-/\kappa_+ = 3.2$; at $r_q = 0.4\,r_s$ the ratio is 16 and every line inside $r_-$ would lie within $10^{-6}$ of the singularity.

After:

> $r_q = 0.48\,r_s$, so that $r_+ = 0.64\,r_s$, $r_- = 0.36\,r_s$, and $\kappa_-/\kappa_+ = 3.2$; at $r_q = 0.4\,r_s$ the ratio is 16, and every line inside $r_-$ would lie within $10^{-6}$ of the singularity.

### embedding diagram, not drawn: `embedding/rn_metric/outside.stops[0]`

Source: `_tools/derivations/embedding.py` line 1839, `in rn_metric()`.

Shown in: `embedding/rn_metric/outside.stops[0]`, `embedding/rn_metric/inside.stops[1]`.

Before:

> Between the horizons, $r_- < r < r_+$, $g_{rr} < 0$: $r$ is a time there and a slice of constant $t$ is not a moment of space, so nothing is drawn.

After:

> Between the horizons, $r_- < r < r_+$, $g_{rr} < 0$: $r$ is a time there, and a slice of constant $t$ is not a moment of space.

## Schwarzschild Black Hole

### history: `schwarzschild/history[0]`

Source: `MFS/assets/data/metrics/schwarzschild.json`, `history`, paragraph 1.

Before:

> In December 1915, weeks after Einstein presented his field equations to the Prussian Academy, Karl Schwarzschild wrote to Einstein from the Russian front with an exact solution.
> The war, he told Einstein, was kindly disposed toward him, allowing him, "despite fierce gunfire at a decidedly terrestrial distance, to take this walk into this your land of ideas" [schulmann1998].
> He had found the gravitational field outside a spherical mass that does not rotate, working under artillery fire in what must rank among the more remarkable acts of mathematical focus in the history of physics.
> Einstein presented the result to the Academy on January 13, 1916 [schwarzschild1916].

After:

> Weeks after Einstein presented his field equations to the Prussian Academy, in December 1915, Karl Schwarzschild wrote to Einstein from the Russian front with an exact solution.
> The war, he told Einstein, was kindly disposed toward him, allowing him, "despite fierce gunfire at a decidedly terrestrial distance, to take this walk into this your land of ideas" [schulmann1998].
> He had found the gravitational field outside a spherical mass that does not rotate, working under artillery fire in what must rank among the more remarkable acts of mathematical focus in the history of physics.
> Einstein presented the result to the Academy on January 13, 1916 [schwarzschild1916].

### history: `schwarzschild/history[2]`

Source: `MFS/assets/data/metrics/schwarzschild.json`, `history`, paragraph 3.

Before:

> Schwarzschild's original coordinates covered only $r > r_s$.
> Ludwig Flamm gave the exterior its first picture in 1916: a plane through the center, measured with the metric's own rulers, has the geometry of the surface made by rotating a parabola about its directrix, the funnel now called Flamm's paraboloid [flamm1916].
> Droste in 1916 and Hilbert in 1917 each rederived the solution with the coordinate domain extended to all $r > 0$, giving the metric the form in which it is used today [droste1917, hilbert1917].
> That the metric carries Schwarzschild's name alone is one of those small injustices history commits without malice.

After:

> Schwarzschild's original coordinates covered only $r > r_s$.
> Ludwig Flamm gave the exterior its first picture in 1916: a plane through the center, measured with the metric's own rulers, has the geometry of the surface made by rotating a parabola about its directrix, the funnel now called Flamm's paraboloid [flamm1916].
> Johannes Droste in 1916 and David Hilbert in 1917 each rederived the solution with the coordinate domain extended to all $r > 0$, giving the metric the form in which it is used today [droste1917, hilbert1917].
> The metric nonetheless carries Schwarzschild's name alone.

### history: `schwarzschild/history[3]`

Source: `MFS/assets/data/metrics/schwarzschild.json`, `history`, paragraph 4.

Before:

> The solution harbors a singularity at $r = r_s = 2GM/c^2$, which puzzled everyone for decades.
> Eddington introduced a coordinate change in 1924 that removed it from the metric components [eddington1924], and Lemaître recognized in 1933 that the singularity was a coordinate artifact, and that freely falling observers cross it without incident [lemaitre1933].
> The full picture came in a burst: Finkelstein in 1958 [finkelstein1958], then Kruskal and Szekeres independently in 1960 [kruskal1960, szekeres1960], who each found maximally extended coordinate systems revealing that $r = r_s$ is not a singularity at all, but a horizon, a membrane that lets matter and light fall in but never back out.
> In Poul Anderson's story "Kyrie" of 1968 a being of plasma falls with a collapsing star, and since by distant clocks the star takes infinite time to shrink to the Schwarzschild radius, its dying cry goes on forever in the mind of the telepath who hears it [anderson1968kyrie].
> In Geoffrey Landis's "Approaching Perimelasma" of 1998 the narrator falls into a black hole that technicians have checked is not rotating, past three Schwarzschild radii, inside which no orbit is stable, and one and a half, where orbiting takes the speed of light, to the horizon where the directions of space and time change identity [landis1998perimelasma].
> Birkhoff had already shown in 1923 that the exterior Schwarzschild metric is the unique spherically symmetric vacuum solution, so any pulsating or collapsing sphere, as long as it remains spherical, wears this geometry outside [birkhoff1923].

After:

> The solution harbors a singularity at $r = r_s = 2GM/c^2$, which puzzled relativists for decades.
> Arthur Eddington introduced a coordinate change in 1924 that removed it from the metric components [eddington1924], and Georges Lemaître recognized in 1933 that the singularity was a coordinate artifact, and that freely falling observers cross it without incident [lemaitre1933].
> The full picture came in a burst: David Finkelstein in 1958 [finkelstein1958], then Martin Kruskal and George Szekeres independently in 1960 [kruskal1960, szekeres1960], who each found maximally extended coordinate systems revealing that $r = r_s$ is a horizon, a membrane that lets matter and light fall in but never back out.
> In Poul Anderson's story "Kyrie" of 1968 a being of plasma falls with a collapsing star, and since by distant clocks the star takes infinite time to shrink to the Schwarzschild radius, its dying cry goes on forever in the mind of the telepath who hears it [anderson1968kyrie].
> In Geoffrey Landis's "Approaching Perimelasma" of 1998 the narrator falls into a black hole that technicians have checked is not rotating, past three Schwarzschild radii, inside which no orbit is stable, and one and a half, where orbiting takes the speed of light, to the horizon where the directions of space and time change identity [landis1998perimelasma].
> George Birkhoff had already shown in 1923 that the exterior Schwarzschild metric is the unique spherically symmetric vacuum solution, so any pulsating or collapsing sphere, as long as it remains spherical, wears this geometry outside [birkhoff1923].

### conventions: `schwarzschild/convention[0]`

Source: `MFS/assets/data/metrics/schwarzschild.json`, `convention`, paragraph 1.

Before:

> The Schwarzschild radius is $r_s = 2GM/c^2$.
> Coordinates are $(t, r, \theta, \phi)$ with $t$ having dimensions of time.
> Factors of $c$ are kept explicit.

After:

> The Schwarzschild radius is $r_s = 2GM/c^2$.
> We use coordinates $(t, r, \theta, \phi)$, with $t$ carrying dimensions of time.
> We keep factors of $c$ explicit.

### spacetime diagram caption: `diagrams/schwarzschild/systems/spherical/radial.caption[1]`

Source: `_tools/derivations/null_rays.py` line 538, `CAPTIONS ("schwarzschild", "spherical", "radial")`.

Before:

> The domain of this chart stops at $r_s$.
> Evaluated inside it, the same components give cones lying on their side, because there it is $r$ that is the time.
> They cannot say whether that region is the black hole or the white hole; the ingoing Eddington-Finkelstein chart, which runs smoothly across $r_s$, makes it the black hole, and there every cone points to $r = 0$.
> The Kretschmann scalar $12r_s^2/r^6$ is finite at $r_s$ and diverges only at $r = 0$.

After:

> The domain of this chart stops at $r_s$.
> Evaluated inside it, the same components give cones lying on their side, because there $r$ is the time.
> The components alone do not fix whether that region is the black hole or the white hole; we take the future from the ingoing Eddington-Finkelstein chart, which runs smoothly across $r_s$, and this makes it the black hole, where every cone points to $r = 0$.
> The Kretschmann scalar $12r_s^2/r^6$ is finite at $r_s$ and diverges only at $r = 0$.

### conformal diagram caption: `conformal/schwarzschild/spherical.caption[0]`

Source: `_tools/derivations/conformal.py` line 2198, `CAPTIONS ("schwarzschild", "spherical")`.

Before:

> This is the whole of the Schwarzschild spacetime, maximally extended, and each point of the diagram stands for a sphere of radius $r$.
> Kruskal and Szekeres's $U = -e^{-u/2r_s}$ and $V = e^{v/2r_s}$, with $u, v = ct \mp r_*$ and $r_* = r + r_s\ln|r/r_s - 1|$, make the metric regular through $r = r_s$, where $UV = (1 - r/r_s)e^{r/r_s}$ vanishes.
> With $p = \arctan U$ and $q = \arctan V$ the singularity $UV = 1$ lies exactly on the straight lines $T = \pm\pi/2$, since $\tan(p + q) = (U + V)/(1 - UV)$ diverges there.

After:

> This is the whole of the Schwarzschild spacetime, maximally extended, and each point of the diagram stands for a sphere of radius $r$.
> The coordinates of Martin Kruskal and George Szekeres, $U = -e^{-u/2r_s}$ and $V = e^{v/2r_s}$, with $u, v = ct \mp r_*$ and $r_* = r + r_s\ln|r/r_s - 1|$, make the metric regular through $r = r_s$, where $UV = (1 - r/r_s)e^{r/r_s}$ vanishes.
> With $p = \arctan U$ and $q = \arctan V$ the singularity $UV = 1$ lies exactly on the straight lines $T = \pm\pi/2$, since $\tan(p + q) = (U + V)/(1 - UV)$ diverges there.

### embedding diagram caption: `embedding/schwarzschild/flamm.caption[0]`

Source: `_tools/derivations/embedding.py` line 3341, `CAPTIONS ("schwarzschild", "flamm")`.

Before:

> This is the equatorial plane $\theta = \pi/2$ of the Schwarzschild spacetime at one moment of $t$, drawn as a surface in flat space so that every distance along it is the metric distance.
> On it the metric is $dr^2/(1 - r_s/r) + r^2d\phi^2$: the circle of radius $r$ has circumference $2\pi r$, while the distance out to the next circle, $dr/\sqrt{1 - r_s/r}$, is longer than $dr$.
> The surface of revolution that carries both is the paraboloid $z^2 = 4r_s(r - r_s)$, which Ludwig Flamm found in 1916.

After:

> This is the equatorial plane $\theta = \pi/2$ of the Schwarzschild spacetime at one moment of $t$, drawn as a surface in flat space so that every distance along it is the metric distance.
> On it the metric is $dr^2/(1 - r_s/r) + r^2d\phi^2$: the circle of radius $r$ has circumference $2\pi r$, while the distance out to the next circle, $dr/\sqrt{1 - r_s/r}$, is longer than $dr$.
> Both hold on the paraboloid $z^2 = 4r_s(r - r_s)$, a surface of revolution that Ludwig Flamm found in 1916.

### embedding diagram caption: `embedding/schwarzschild/flamm.caption[1]`

Source: `_tools/derivations/embedding.py` line 3347, `CAPTIONS ("schwarzschild", "flamm")`.

Before:

> Every slice of constant $t$ passes through the bifurcation sphere $r = r_s$, where the circles are smallest, and runs on through it into a second exterior, the same paraboloid turned over.
> Einstein and Rosen took this bridge between the two sheets as a model of a particle in 1935.
> Robert Fuller and John Wheeler showed in 1962 that its throat closes before light can cross it, so nothing passes from one exterior to the other.

After:

> Every slice of constant $t$ passes through the bifurcation sphere $r = r_s$, where the circles are smallest, and runs on through it into a second exterior, the same paraboloid turned over.
> Albert Einstein and Nathan Rosen took this bridge between the two sheets as a model of a particle in 1935.
> Robert Fuller and John Wheeler showed in 1962 that its throat closes before light can cross it, so nothing passes from one exterior to the other.

## Lanczos-van Stockum Rotating Dust Universe

### history: `stockum_dust/history[0]`

Source: `MFS/assets/data/metrics/stockum_dust.json`, `history`, paragraph 1.

Before:

> Cornelius Lanczos wrote down the rotating dust cylinder in 1924, not long after Einstein's field equations had settled into their final form [lanczos1924].
> He was not hunting for physical models.
> The solution describes an infinite cylinder of pressureless dust in rigid rotation, a geometry with no reasonable astrophysical realization, but exact solutions were scarce, and one took what one could get.
> It was largely set aside, its sting unnoticed.

After:

> Cornelius Lanczos wrote down the rotating dust cylinder in 1924, not long after Einstein's field equations had settled into their final form [lanczos1924].
> He was not hunting for physical models.
> The solution describes an infinite cylinder of pressureless dust in rigid rotation, a geometry with no reasonable astrophysical realization, but exact solutions were scarce, and one took what one could get.
> It was largely set aside, its closed timelike curves unnoticed.

### history: `stockum_dust/history[1]`

Source: `MFS/assets/data/metrics/stockum_dust.json`, `history`, paragraph 2.

Before:

> Willem Jacob van Stockum returned to it in 1937 and worked it through carefully, mapping the causal structure of the spacetime [vanstockum1937].
> He was working at the Mathematical Institute of the University of Edinburgh, and E. T. Whittaker communicated his paper to the Royal Society of Edinburgh [vanstockum1937].
> He had set out to extend Weyl's static fields with an axis of symmetry to rotating ones, and solved the case of dust rotating rigidly about an axis, both inside a cylinder and in the empty space around it [vanstockum1937, bonnor1980].
> He noticed that for a dense enough cylinder a light signal sent against the rotation runs along a circle on its surface, so that an observer there could look right round the cylinder, the light returning after about forty minutes [vanstockum1937].

After:

> Willem Jacob van Stockum returned to it in 1937 and worked it through carefully, mapping the causal structure of the spacetime [vanstockum1937].
> He was working at the Mathematical Institute of the University of Edinburgh, and Edmund Whittaker communicated his paper to the Royal Society of Edinburgh [vanstockum1937].
> He had set out to extend Weyl's static fields with an axis of symmetry to rotating ones, and solved the case of dust rotating rigidly about an axis, both inside a cylinder and in the empty space around it [vanstockum1937, bonnor1980].
> He noticed that for a dense enough cylinder a light signal sent against the rotation runs along a circle on its surface, so that an observer there could look right round the cylinder, the light returning after about forty minutes [vanstockum1937].

### history: `stockum_dust/history[2]`

Source: `MFS/assets/data/metrics/stockum_dust.json`, `history`, paragraph 3.

Before:

> He found the sting, buried in the geometry all along: at $r = R$, where $R$ is the cylinder radius, a closed null curve first appears.
> Beyond that radius the spacetime harbors closed timelike curves, trajectories that loop back to their own past.
> A particle following such a curve would return to the event it started from, having never left the forward flow of its own proper time.
> The Lanczos-van Stockum spacetime was, by more than a decade, the first exact solution of general relativity known to contain them, predating Gödel's rotating universe of 1949.

After:

> He found the closed curves the geometry had held all along: at $r = R$, where $R$ is the cylinder radius, a closed null curve first appears.
> Beyond that radius the spacetime harbors closed timelike curves, trajectories that loop back to their own past.
> A particle following such a curve would return to the event it started from, having never left the forward flow of its own proper time.
> The Lanczos-van Stockum spacetime was, by more than a decade, the first exact solution of general relativity known to contain them, predating Gödel's rotating universe of 1949.

### history: `stockum_dust/history[3]`

Source: `MFS/assets/data/metrics/stockum_dust.json`, `history`, paragraph 4.

Before:

> Frank Tipler took up the solution in 1974 and showed that the field of the rotating cylinder allows a closed timelike line to connect any two events, which led him to suggest that a finite rotating cylinder would also act as a time machine [tipler1974].
> Two years later he proved that a region holding closed timelike lines cannot grow out of regular initial data in an asymptotically flat spacetime free of singularities [tipler1976].
> Larry Niven took the title of Tipler's paper for a story in Analog in 1977, in which a mathematician offers the emperor of seventy worlds a time machine made by putting a rapid spin on a massive cylinder [niven1977cylinders].
> William Bonnor showed in 1980 that van Stockum's cylinder has three kinds of exterior according to its mass per unit length, that the exterior is static only for the lightest, and that the heaviest, which contains closed timelike lines, can occur at physically possible densities and radii [bonnor1980].
> The solution remains a theoretical laboratory rather than a physical proposal, but as a precise arena for questions about causality, cosmic censorship, and the permissiveness of Einstein's equations, it has lost none of its bite [misner1973].

After:

> Frank Tipler took up the solution in 1974 and showed that the field of the rotating cylinder allows a closed timelike line to connect any two events, which led him to suggest that a finite rotating cylinder would also act as a time machine [tipler1974].
> Two years later he proved that a region holding closed timelike lines cannot grow out of regular initial data in an asymptotically flat spacetime free of singularities [tipler1976].
> Larry Niven took the title of Tipler's paper for a story in Analog in 1977, in which a mathematician offers the emperor of seventy worlds a time machine made by putting a rapid spin on a massive cylinder [niven1977cylinders].
> William Bonnor showed in 1980 that van Stockum's cylinder has three kinds of exterior according to its mass per unit length, that the exterior is static only for the lightest, and that the heaviest, which contains closed timelike lines, can occur at physically possible densities and radii [bonnor1980].
> The solution remains a theoretical laboratory rather than a physical proposal, and it is still used to pose precise questions about causality, cosmic censorship, and what Einstein's equations permit [misner1973].

### conventions: `stockum_dust/convention[0]`

Source: `MFS/assets/data/metrics/stockum_dust.json`, `convention`, paragraph 1.

Before:

> Cylindrical coordinates $(t, r, \phi, z)$ are used.
> The parameter $R$ sets the scale of the dust cylinder; $r = R$ is the critical radius at which closed null curves first appear.
> The metric function $e^{-r^2/R^2}$ appears in the $g_{rr}$ and $g_{zz}$ components.
> The $g_{t\phi}$ term off the diagonal encodes the frame dragging due to rotation.
> The Ricci tensor is the standard contraction $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$, on which the Einstein tensor is $G_{\mu\nu} = 8\pi\rho\,u_\mu u_\nu$ for the comoving dust, of positive density $\rho = \dfrac{e^{r^2/R^2}}{2\pi R^2}$.
> Units with $G = c = 1$ are used throughout.

After:

> We use cylindrical coordinates $(t, r, \phi, z)$.
> The parameter $R$ sets the scale of the dust cylinder, and $r = R$ is the critical radius at which closed null curves first appear.
> The function $e^{-r^2/R^2}$ enters the metric in the $g_{rr}$ and $g_{zz}$ components.
> The term $g_{t\phi}$, off the diagonal, carries the dragging of frames by the rotation.
> With the Ricci tensor the standard contraction, $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$, the Einstein tensor is $G_{\mu\nu} = 8\pi\rho\,u_\mu u_\nu$ for the comoving dust, of positive density $\rho = \dfrac{e^{r^2/R^2}}{2\pi R^2}$.
> We work in units with $G = c = 1$ throughout.

### spacetime diagram caption: `diagrams/stockum_dust/systems/cylindrical/beyond.caption[0]`

Source: `_tools/derivations/null_rays.py` line 922, `CAPTIONS ("stockum_dust", "cylindrical", "beyond")`.

Before:

> This is the cylinder of $t$ and $\phi$ at $r = 3R/2$ and $z = 0$, opened along $\phi = \pm\pi$ in the same way.
> Beyond $r = R$ the coefficient $g_{\phi\phi} = r^2(1 - r^2/R^2)$ is negative, and the cones have tipped over past the horizontal: the null curve moving to $+\phi$, $dt = r(1 - r/R)\,d\phi$, goes down in $t$, while the one moving to $-\phi$, $dt = -r(1 + r/R)\,d\phi$, climbs steeply.
> Every horizontal line, run toward $+\phi$, points into the future cones, so the circle of constant $t$, $r$ and $z$ is a closed timelike curve.

After:

> This is the cylinder of $t$ and $\phi$ at $r = 3R/2$ and $z = 0$, opened along $\phi = \pm\pi$ in the same way.
> Beyond $r = R$ the coefficient $g_{\phi\phi} = r^2(1 - r^2/R^2)$ is negative, and the cones have tipped over past the horizontal: the null curve moving to $+\phi$, $dt = r(1 - r/R)\,d\phi$, goes down in $t$, while the one moving to $-\phi$, $dt = -r(1 + r/R)\,d\phi$, climbs steeply.
> Every horizontal line, run toward $+\phi$, points into the future cones, so the circle of constant $t$, $r$, and $z$ is a closed timelike curve.

### caption of a figure in three dimensions: `diagrams/stockum_dust/projections/cylindrical/tipping.caption[0]`

Source: `_tools/derivations/projections.py` line 553, `CAPTIONS ("stockum_dust", "cylindrical", "tipping")`.

Before:

> This is the slice $z = 0$ of $t$, $r$ and $\phi$, with $t$ up and the proper distance from the axis, $\int e^{-r^2/2R^2}dr$, as the radius, which puts the null directions straight out from the axis at 45°.
> The cones stand at $t = 0$ on the axis and around the circles $r = R/2$, $R$ and $3R/2$.
> On the axis they are upright, and farther out the cross term $g_{t\phi} = -r^2/R$ tips them over toward $+\phi$, counterclockwise seen from above.

After:

> This is the slice $z = 0$ of $t$, $r$, and $\phi$, with $t$ up and the proper distance from the axis, $\int e^{-r^2/2R^2}dr$, as the radius, which puts the null directions straight out from the axis at 45°.
> The cones stand at $t = 0$ on the axis and around the circles $r = R/2$, $R$, and $3R/2$.
> On the axis they are upright, and farther out the cross term $g_{t\phi} = -r^2/R$ tips them over toward $+\phi$, counterclockwise seen from above.

## Taub-Newman-Unti-Tamburino Vacuum Spacetime

### history: `taub_nut/history[0]`

Source: `MFS/assets/data/metrics/taub_nut.json`, `history`, paragraph 1.

Before:

> Taub-NUT is two discoveries that turned out to be one spacetime.
> In 1951 Abraham Taub, looking for vacuum solutions with a high degree of symmetry, wrote down an empty cosmology whose spatial slices carry the symmetry group of the sphere $S^3$, a homogeneous but anisotropic universe that expands and recollapses [taub1951].
> It was a clean result that stood on its own, and for a decade it sat among the catalogue of exact solutions without causing trouble.

After:

> Taub-NUT is two discoveries that turned out to be one spacetime.
> In 1951 Abraham Taub, looking for vacuum solutions with a high degree of symmetry, wrote down an empty cosmology whose spatial slices carry the symmetry group of the sphere $S^3$, a homogeneous but anisotropic universe that expands and recollapses [taub1951].
> It was a clean result, not yet joined to any other solution, and for a decade it sat among the catalogue of exact solutions without causing trouble.

### history: `taub_nut/history[1]`

Source: `MFS/assets/data/metrics/taub_nut.json`, `history`, paragraph 2.

Before:

> The trouble arrived in 1963.
> Ezra Newman, Louis Tamburino, and Theodore Unti set out to generalize the Schwarzschild metric and found a new vacuum solution carrying one extra parameter beyond the mass, the NUT parameter, after their initials [nut1963].
> That same year Charles Misner showed that their "generalized Schwarzschild" space and Taub's cosmology were not two solutions but one: Taub's universe is the interior, the NUT metric the exterior, joined across a horizon into a single analytic whole [misner1963].
> Hence Taub-NUT.

After:

> The trouble arrived in 1963.
> Ezra Newman, Louis Tamburino, and Theodore Unti set out to generalize the Schwarzschild metric and found a new vacuum solution carrying one extra parameter beyond the mass, the NUT parameter, after their initials [nut1963].
> That same year Charles Misner showed that their "generalized Schwarzschild" space and Taub's cosmology were one solution: Taub's universe is the interior, the NUT metric the exterior, joined across a horizon into a single analytic whole [misner1963].
> Hence Taub-NUT.

### history: `taub_nut/history[2]`

Source: `MFS/assets/data/metrics/taub_nut.json`, `history`, paragraph 3.

Before:

> Misner also showed that the joined spacetime is spectacularly badly behaved.
> The NUT parameter acts like a magnetic charge for gravity, a "gravitomagnetic monopole", and like the Dirac string trailing a magnetic monopole, it drags along an unavoidable singular axis, the Misner string, that no choice of coordinates can remove without introducing closed timelike curves through every point.
> The spacetime is not asymptotically flat in any ordinary sense [misner1963] and is not globally hyperbolic [hawking1973].
> Misner's own verdict was that it serves as "a counterexample to almost any conjecture one might make," and that is exactly why it endures: Taub-NUT is the geometry relativists reach for when they want to know whether a theorem that sounds plausible is actually true.

After:

> Misner also showed that the joined spacetime is spectacularly badly behaved.
> The NUT parameter acts like a magnetic charge for gravity, a "gravitomagnetic monopole", and like the Dirac string trailing a magnetic monopole, it drags along an unavoidable singular axis, the Misner string, that no choice of coordinates can remove without introducing closed timelike curves through every point.
> The spacetime is not asymptotically flat in any ordinary sense [misner1963] and is not globally hyperbolic [hawking1973].
> Misner's own verdict was that it serves as "a counterexample to almost any conjecture one might make," and relativists still turn to it to test whether a plausible theorem holds.

### history: `taub_nut/history[3]`

Source: `MFS/assets/data/metrics/taub_nut.json`, `history`, paragraph 4.

Before:

> William Bonnor read the singular half axis physically in 1969, as the field of a mass together with a massless source of angular momentum running out to infinity along half the axis [bonnor1969].
> Stephen Hawking and George Ellis gave Taub-NUT a section of its own in their 1973 monograph, as a spacetime of interest for its pathological global properties [hawking1973].
> They used it as their example of a spacetime whose incomplete geodesics are imprisoned in a compact region where the curvature stays regular, so that an observer on one of them would wind round and round forever yet never get beyond a certain time in their life [hawking1973].

After:

> William Bonnor interpreted the singular half axis physically in 1969, as the field of a mass together with a massless source of angular momentum running out to infinity along half the axis [bonnor1969].
> Stephen Hawking and George Ellis gave Taub-NUT a section of its own in their 1973 monograph, as a spacetime of interest for its pathological global properties [hawking1973].
> They used it as their example of a spacetime whose incomplete geodesics are imprisoned in a compact region where the curvature stays regular, so that an observer on one of them would wind round and round forever yet never get beyond a certain time in their life [hawking1973].

### conventions: `taub_nut/convention[0]`

Source: `MFS/assets/data/metrics/taub_nut.json`, `convention`, paragraph 1.

Before:

> Coordinates are $(t, r, \theta, \phi)$ with $t$ having dimensions of time, and factors of $c$ and $G$ are kept explicit.
> The mass parameter is the geometric mass $m = GM/c^2$ and the NUT parameter $l$ is a length as well, so $G$ and the mass $M$ enter only through $m$.
> Components are taken in the chart $x^0 = ct$, where the metric is dimensionless and no component carries a factor of $c$, and the dots in the geodesic equations are derivatives of that same chart with respect to an affine parameter, so $\dot{t}$ means $d(ct)/d\lambda$ and every term of every equation carries the dimensions of its left hand side.

After:

> We use coordinates $(t, r, \theta, \phi)$, with $t$ carrying dimensions of time, and keep factors of $c$ and $G$ explicit.
> The mass parameter is the geometric mass $m = GM/c^2$, and the NUT parameter $l$ is a length as well, so $G$ and the mass $M$ enter only through $m$.
> We take components in the chart $x^0 = ct$, where the metric is dimensionless and no component carries a factor of $c$, and a dot in the geodesic equations is a derivative of that chart's coordinates with respect to an affine parameter, so $\dot{t}$ means $d(ct)/d\lambda$, and every term of every equation carries the dimensions of its left hand side.

### conventions: `taub_nut/convention[1]`

Source: `MFS/assets/data/metrics/taub_nut.json`, `convention`, paragraph 2.

Before:

> The geodesic equations are $\ddot{x}^\mu + \Gamma^\mu{}_{\nu\rho}\dot{x}^\nu\dot{x}^\rho = 0$ with $\Gamma$ the Christoffel symbols of that same chart.
> The signature is $(-,+,+,+)$, the Riemann tensor is $R^\mu{}_{\nu\rho\sigma} = \partial_\rho\Gamma^\mu{}_{\nu\sigma} - \partial_\sigma\Gamma^\mu{}_{\nu\rho} + \Gamma^\mu{}_{\rho\lambda}\Gamma^\lambda{}_{\nu\sigma} - \Gamma^\mu{}_{\sigma\lambda}\Gamma^\lambda{}_{\nu\rho}$, and the Ricci tensor is the standard contraction $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$, a choice that is invisible here because the Ricci tensor vanishes either way.
> The solution is a vacuum one: $R_{\mu\nu} = 0$, $R = 0$ and $G_{\mu\nu} = 0$, and the Weyl tensor is therefore equal to the Riemann tensor component by component, which is a consequence of the field equations.

After:

> The geodesic equations read $\ddot{x}^\mu + \Gamma^\mu{}_{\nu\rho}\dot{x}^\nu\dot{x}^\rho = 0$, with $\Gamma$ the Christoffel symbols of the same chart.
> We work in the signature $(-,+,+,+)$, with the Riemann tensor $R^\mu{}_{\nu\rho\sigma} = \partial_\rho\Gamma^\mu{}_{\nu\sigma} - \partial_\sigma\Gamma^\mu{}_{\nu\rho} + \Gamma^\mu{}_{\rho\lambda}\Gamma^\lambda{}_{\nu\sigma} - \Gamma^\mu{}_{\sigma\lambda}\Gamma^\lambda{}_{\nu\rho}$ and the Ricci tensor its standard contraction, $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$, a choice that makes no difference here because the Ricci tensor vanishes either way.
> The solution is a vacuum: $R_{\mu\nu} = 0$, $R = 0$, and $G_{\mu\nu} = 0$, and the Weyl tensor is therefore equal to the Riemann tensor component by component, as a consequence of the field equations.

### conventions: `taub_nut/convention[2]`

Source: `MFS/assets/data/metrics/taub_nut.json`, `convention`, paragraph 3.

Before:

> The chart stops short of the axis, and that is the Misner string.
> The one form $c\,dt + 2l\cos\theta\,d\phi$ needs $\cos\theta$ to vanish where $d\phi$ is undefined, and $\cos\theta$ is $+1$ at $\theta = 0$ and $-1$ at $\theta = \pi$, so in this gauge the potential is singular along both halves of the axis and $\theta$ is restricted to $(0, \pi)$.
> The gauge $2l(\cos\theta \mp 1)d\phi$ clears one half and doubles the other, and Misner's own repair, making $t$ periodic with period $8\pi l/c$, clears both at the price of closed timelike curves through every point.
> Nothing about the string is curvature: the Kretschmann scalar is finite on the axis and blows up only at $r^2 + l^2 = 0$, which no real $r$ reaches while $l \neq 0$.

After:

> The chart stops short of the axis, where the Misner string lies.
> The one form $c\,dt + 2l\cos\theta\,d\phi$ needs $\cos\theta$ to vanish where $d\phi$ is undefined, and $\cos\theta$ is $+1$ at $\theta = 0$ and $-1$ at $\theta = \pi$, so in this gauge the potential is singular along both halves of the axis, and $\theta$ is restricted to $(0, \pi)$.
> The gauge $2l(\cos\theta \mp 1)d\phi$ clears one half and doubles the other, and Charles Misner's own repair, making $t$ periodic with period $8\pi l/c$, clears both at the price of closed timelike curves through every point.
> The string carries no curvature: the Kretschmann scalar is finite on the axis and blows up only at $r^2 + l^2 = 0$, which no real $r$ reaches while $l \neq 0$.

### spacetime diagram caption: `diagrams/taub_nut/systems/spherical/radial.caption[1]`

Source: `_tools/derivations/null_rays.py` line 708, `CAPTIONS ("taub_nut", "spherical", "radial")`.

Before:

> The twist the NUT parameter brings sits in the cross term $g_{t\phi}$, which drops out of the metric on this plane.
> No Christoffel symbol turns these null curves out of the plane either, so they are null geodesics, the paths light takes, even though the solution is not spherically symmetric.

After:

> The NUT parameter's twist sits in the cross term $g_{t\phi}$, which drops out of the metric on this plane.
> No Christoffel symbol turns these null curves out of the plane either, so they are null geodesics, the paths light takes, even though the solution is not spherically symmetric.

## Tolman-Bondi Inhomogeneous Dust Universe

### history: `tolman_bondi/history[1]`

Source: `MFS/assets/data/metrics/tolman_bondi.json`, `history`, paragraph 2.

Before:

> The solution's charm is that it is exactly solvable despite being inhomogeneous, and the reason is a kind of independence: each concentric shell of dust evolves on its own, like a miniature FLRW universe with its own spatial curvature and its own clock, neighboring shells neither knowing nor caring what the others do.
> Set the radial profiles to constants and ordinary FLRW reappears; take a single uniform ball and match it to an exterior vacuum and one recovers the Oppenheimer-Snyder model of collapse [oppenheimersnyder1939, christodoulou1984].
> W. B. Bonnor showed in 1974 that an open Robertson-Walker universe can grow from many initial states, since there are inhomogeneous models of this kind whose initial density is an arbitrary function of radius and which approach the Robertson-Walker form as time goes on [bonnor1974].
> LTB is the connective tissue between the homogeneous cosmos and the collapsing star.

After:

> The solution's charm is that it is exactly solvable despite being inhomogeneous, and the reason is a kind of independence: each concentric shell of dust evolves on its own, like a miniature FLRW universe with its own spatial curvature and its own clock, neighboring shells neither knowing nor caring what the others do.
> Set the radial profiles to constants and ordinary FLRW reappears; take a single uniform ball and match it to an exterior vacuum and one recovers the Oppenheimer-Snyder model of collapse [oppenheimersnyder1939, christodoulou1984].
> William Bonnor showed in 1974 that an open Robertson-Walker universe can grow from many initial states, since there are inhomogeneous models of this kind whose initial density is an arbitrary function of radius and which approach the Robertson-Walker form as time goes on [bonnor1974].
> LTB thus holds both the homogeneous cosmos and the collapsing star as special cases.

### history: `tolman_bondi/history[3]`

Source: `MFS/assets/data/metrics/tolman_bondi.json`, `history`, paragraph 4.

Before:

> That flexibility has made it the standard tool for inhomogeneity in general relativity.
> Lemaître had already used it to discuss how clusters of nebulae might form, an attempt within the exact theory that Krasiński judged at least twenty years ahead of its time [krasinski1997].
> It models the growth of cosmic voids and the nonlinear formation of structure, and it has been pressed into more provocative service: a sufficiently large local void, described by an LTB metric, can mimic the dimming of distant supernovae that is usually read as cosmic acceleration, raising the question of whether some of what we attribute to dark energy is instead a sign that we sit near the center of an underdensity [tomita2001].

After:

> The solution's flexibility has made it the standard tool for inhomogeneity in general relativity.
> Lemaître had already used it to discuss how clusters of nebulae might form, an attempt within the exact theory that Krasiński judged at least twenty years ahead of its time [krasinski1997].
> It models the growth of cosmic voids and the nonlinear formation of structure, and it has been pressed into more provocative service: a sufficiently large local void, described by an LTB metric, can mimic the dimming of distant supernovae that is usually taken as cosmic acceleration, raising the question of whether some of what we attribute to dark energy is instead a sign that we sit near the center of an underdensity [tomita2001].

### history: `tolman_bondi/history[4]`

Source: `MFS/assets/data/metrics/tolman_bondi.json`, `history`, paragraph 5.

Before:

> Marie-Noëlle Célérier gave in 2000 an inhomogeneous model with no cosmological constant that reproduced the supernova data [celerier2000], and Håvard Alnes, Morad Amarzguioui, and Øyvind Grøn found in 2006 that an underdense bubble centered near the observer fits the supernova distances without dark energy [alnes2006].
> Adam Moss, James Zibin, and Douglas Scott found void models in severe tension with the data in 2011: they predict too low a local Hubble rate and much less local structure than is observed [moss2011].
> Pengjie Zhang and Albert Stebbins showed the same year that a void large enough to stand in for dark energy would raise a kinetic Sunyaev-Zel'dovich signal far above the observed limits [zhang2011].
> The proposal is constrained and largely disfavored [moss2011, zhang2011], but that it can even be posed is a measure of how much room inhomogeneity leaves in the field equations.

After:

> Marie-Noëlle Célérier gave in 2000 an inhomogeneous model with no cosmological constant that reproduced the supernova data [celerier2000], and Håvard Alnes, Morad Amarzguioui, and Øyvind Grøn found in 2006 that an underdense bubble centered near the observer fits the supernova distances without dark energy [alnes2006].
> Adam Moss, James Zibin, and Douglas Scott found void models in severe tension with the data in 2011: they predict too low a local Hubble rate and much less local structure than is observed [moss2011].
> Pengjie Zhang and Albert Stebbins showed the same year that a void large enough to stand in for dark energy would raise a kinetic Sunyaev-Zel'dovich signal far above the observed limits [zhang2011].
> The proposal is constrained and largely disfavored [moss2011, zhang2011], by observation alone, since the field equations allow such a void.

### conventions: `tolman_bondi/convention[0]`

Source: `MFS/assets/data/metrics/tolman_bondi.json`, `convention`, paragraph 1.

Before:

> Coordinates are $(t, r, \theta, \phi)$ with $t$ a cosmic time and $r$ a comoving label for a shell of matter, and the areal radius $R(r,t)$ carries the size of that shell, whose sphere has area $4\pi R^2$.
> Factors of $c$ and $G$ are kept explicit.
> Components are taken in the chart $x^0 = ct$, where a chart component of the metric carries no dimensions, and every $\partial_t$ is the derivative along that same chart coordinate, $\partial_t = c^{-1}\partial/\partial t$, so $\partial_t R$ and $\partial_r R$ are dimensionless and each further derivative carries one more inverse length; the energy function depends on $r$ alone, so its derivative is written $E' = dE/dr$.
> Apart from the $-c^2dt^2$ of the line element no factor of $c$ appears in any component, because that one definition holds them all.
> The dots in the geodesic equations are velocities of that same chart, so $\dot{t}$ means $d(ct)/d\lambda$, and the equations are $\ddot{x}^\mu + \Gamma^\mu{}_{\nu\rho}\dot{x}^\nu\dot{x}^\rho = 0$ with $\Gamma$ the Christoffel symbols of that same chart.
> The signature is $(-,+,+,+)$, the Riemann tensor is $R^\mu{}_{\nu\rho\sigma} = \partial_\rho\Gamma^\mu{}_{\nu\sigma} - \partial_\sigma\Gamma^\mu{}_{\nu\rho} + \Gamma^\mu{}_{\rho\lambda}\Gamma^\lambda{}_{\nu\sigma} - \Gamma^\mu{}_{\sigma\lambda}\Gamma^\lambda{}_{\nu\rho}$, and the Ricci tensor is the standard contraction $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$.

After:

> We use coordinates $(t, r, \theta, \phi)$, with $t$ a cosmic time and $r$ a comoving label for a shell of matter, and the areal radius $R(r,t)$ carries the size of that shell, whose sphere has area $4\pi R^2$.
> We keep factors of $c$ and $G$ explicit.
> We take components in the chart $x^0 = ct$, where a chart component of the metric carries no dimensions, and every $\partial_t$ is the derivative along that same chart coordinate, $\partial_t = c^{-1}\partial/\partial t$, so $\partial_t R$ and $\partial_r R$ are dimensionless, and each further derivative carries one more inverse length; the energy function depends on $r$ alone, so we write its derivative $E' = dE/dr$.
> Apart from the $-c^2dt^2$ of the line element no factor of $c$ appears in any component, because that one definition holds them all.
> The dots in the geodesic equations are velocities in the same chart, so $\dot{t}$ means $d(ct)/d\lambda$, and the equations read $\ddot{x}^\mu + \Gamma^\mu{}_{\nu\rho}\dot{x}^\nu\dot{x}^\rho = 0$, with $\Gamma$ the Christoffel symbols of that chart.
> We work in the signature $(-,+,+,+)$, with the Riemann tensor $R^\mu{}_{\nu\rho\sigma} = \partial_\rho\Gamma^\mu{}_{\nu\sigma} - \partial_\sigma\Gamma^\mu{}_{\nu\rho} + \Gamma^\mu{}_{\rho\lambda}\Gamma^\lambda{}_{\nu\sigma} - \Gamma^\mu{}_{\sigma\lambda}\Gamma^\lambda{}_{\nu\rho}$ and the Ricci tensor its standard contraction, $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$.

### conventions: `tolman_bondi/convention[1]`

Source: `MFS/assets/data/metrics/tolman_bondi.json`, `convention`, paragraph 2.

Before:

> The radial coefficient is not a guess that happens to work.
> For the general comoving synchronous metric $ds^2 = -c^2dt^2 + X^2dr^2 + R^2d\Omega^2$ the mixed Einstein component is $G_{tr} = -\dfrac{2X}{R}\partial_t\!\left(\dfrac{\partial_r R}{X}\right)$, which vanishes exactly when $X = \partial_r R/\sqrt{1+2E}$ for some $E(r)$, so this form is the general solution of the one field equation that says no energy crosses a comoving shell, and $G_{tr}$ is identically zero.
> No other component assumes a field equation.
> Every component is an identity for an arbitrary $R(r,t)$ and an arbitrary $E(r)$, and the matter is read straight off the Einstein tensor as a fluid at rest with density $\rho = \dfrac{c^2}{8\pi G}G_{tt}$, radial pressure $p_r = \dfrac{c^4}{8\pi G}G^r{}_r$ and transverse pressure $p_\perp = \dfrac{c^4}{8\pi G}G^\theta{}_\theta$.

After:

> One field equation forces the radial coefficient.
> For the general comoving synchronous metric $ds^2 = -c^2dt^2 + X^2dr^2 + R^2d\Omega^2$ the mixed Einstein component is $G_{tr} = -\dfrac{2X}{R}\partial_t\!\left(\dfrac{\partial_r R}{X}\right)$, which vanishes exactly when $X = \partial_r R/\sqrt{1+2E}$ for some $E(r)$, so this form is the general solution of the one field equation that says no energy crosses a comoving shell, and $G_{tr}$ is identically zero.
> No other component assumes a field equation.
> Every component is an identity for an arbitrary $R(r,t)$ and an arbitrary $E(r)$, and the matter comes straight from the Einstein tensor as a fluid at rest with density $\rho = \dfrac{c^2}{8\pi G}G_{tt}$, radial pressure $p_r = \dfrac{c^4}{8\pi G}G^r{}_r$, and transverse pressure $p_\perp = \dfrac{c^4}{8\pi G}G^\theta{}_\theta$.

### conventions: `tolman_bondi/convention[2]`

Source: `MFS/assets/data/metrics/tolman_bondi.json`, `convention`, paragraph 3.

Before:

> Dust is then one equation and not three.
> Setting $p_r = 0$ gives $2R\,\partial_t^2R + (\partial_t R)^2 - 2E = 0$, which integrates in time to the Lemaître-Tolman-Bondi evolution equation $(\partial_t R)^2 = 2E + \dfrac{2GM(r)}{c^2R}$ with the mass function $M(r)$ as its constant of integration, and the transverse pressure vanishes along with it rather than by a second assumption, because $G^\theta{}_\theta = \dfrac{\partial_r\left(R^2G^r{}_r\right)}{2R\,\partial_r R}$ is an identity in this chart.
> On that equation the time component becomes $G_{tt} = \dfrac{2GM'}{c^2R^2\,\partial_r R} = \dfrac{8\pi G}{c^2}\rho$, so the density is $\rho(r,t) = \dfrac{M'(r)}{4\pi R^2\,\partial_r R}$: the gravitating mass inside a shell is carried along unchanged, and the density is its radial gradient per unit areal volume, which diverges both where a shell collapses to $R = 0$ and where two shells meet at $\partial_r R = 0$.
> FLRW is the homogeneous limit, reached by $R = a(t)r$ and $E = -kr^2/2$: the line element becomes the comoving spherical FRW one, $G_{tt}$ becomes the Friedmann equation $3(a'^2+k)/a^2$, the mass function becomes $2GM/c^2 = ar^3(a'^2+k)$ and the density becomes independent of $r$.

After:

> Dust then needs a single condition.
> Setting $p_r = 0$ gives $2R\,\partial_t^2R + (\partial_t R)^2 - 2E = 0$, which integrates in time to the Lemaître-Tolman-Bondi evolution equation $(\partial_t R)^2 = 2E + \dfrac{2GM(r)}{c^2R}$ with the mass function $M(r)$ as its constant of integration, and the transverse pressure vanishes along with it rather than by a second assumption, because $G^\theta{}_\theta = \dfrac{\partial_r\left(R^2G^r{}_r\right)}{2R\,\partial_r R}$ is an identity in this chart.
> With that equation the time component becomes $G_{tt} = \dfrac{2GM'}{c^2R^2\,\partial_r R} = \dfrac{8\pi G}{c^2}\rho$, so the density is $\rho(r,t) = \dfrac{M'(r)}{4\pi R^2\,\partial_r R}$: the gravitating mass inside a shell is carried along unchanged, and the density is its radial gradient per unit areal volume, which diverges both where a shell collapses to $R = 0$ and where two shells meet at $\partial_r R = 0$.
> FLRW is the homogeneous limit, reached by $R = a(t)r$ and $E = -kr^2/2$: the line element becomes the comoving spherical FRW one, $G_{tt}$ becomes the Friedmann equation $3(a'^2+k)/a^2$, the mass function becomes $2GM/c^2 = ar^3(a'^2+k)$, and the density becomes independent of $r$.

### conventions: `tolman_bondi/convention[3]`

Source: `MFS/assets/data/metrics/tolman_bondi.json`, `convention`, paragraph 4.

Before:

> This is not a vacuum, so the Weyl tensor is not the Riemann tensor but the Riemann tensor with its traces removed.
> Every nonzero component is proportional to one scalar, which on the dust equation is $C^t{}_{\theta t\theta} = \dfrac{4\pi G\rho R^2}{3c^2} - \dfrac{GM}{c^2R}$, the difference between the mass a uniform ball at the local density would hold inside the shell and the mass actually inside it.
> That is the tidal field of the inhomogeneity alone, and it is why all twenty four Weyl components vanish in the FLRW limit while the Riemann tensor does not.
> The chart is regular where $R > 0$, $\partial_r R > 0$ and $1 + 2E > 0$.

After:

> The spacetime is not a vacuum, and its Weyl tensor is the Riemann tensor with the traces removed.
> Every nonzero component is proportional to one scalar, which with the dust equation is $C^t{}_{\theta t\theta} = \dfrac{4\pi G\rho R^2}{3c^2} - \dfrac{GM}{c^2R}$, the difference between the mass a uniform ball at the local density would hold inside the shell and the mass inside it.
> This scalar is the tidal field of the inhomogeneity alone, which is why all twenty four Weyl components vanish in the FLRW limit while the Riemann tensor does not.
> The chart is regular where $R > 0$, $\partial_r R > 0$, and $1 + 2E > 0$.

### spacetime diagram caption: `diagrams/tolman_bondi/systems/comoving_synchronous/collapse.caption[0]`

Source: `_tools/derivations/null_rays.py` line 1019, `CAPTIONS ("tolman_bondi", "comoving_synchronous", "collapse")`.

Before:

> This is the plane of $t$ and the comoving $r$ at $\theta = \pi/2$ and $\phi = 0$, through a cloud of dust whose density falls from its centre to zero at its surface $r_b$, with vacuum outside.
> Every shell falls on its own clock, $R^{3/2} = r^{3/2} - \tfrac{3}{2}\sqrt{2GM(r)/c^2}\,ct$, and reaches $R = 0$ at its own time: the centre first, at $ct = 0.60\,r_b$, and the surface at $0.94\,r_b$.
> The singularity, where the Kretschmann scalar diverges, is that curve, and beyond it there is no spacetime.
> The rays obey $c\,dt = \pm\partial_r R\,dr$, and outside the cloud the same coordinates are Lemaître's for Schwarzschild's exterior, carried by observers who fall freely from rest at infinity.

After:

> This is the plane of $t$ and the comoving $r$ at $\theta = \pi/2$ and $\phi = 0$, through a cloud of dust whose density falls from its centre to zero at its surface $r_b$, with vacuum outside.
> Every shell falls on its own clock, $R^{3/2} = r^{3/2} - \tfrac{3}{2}\sqrt{2GM(r)/c^2}\,ct$, and reaches $R = 0$ at its own time: the centre first, at $ct = 0.60\,r_b$, and the surface at $0.94\,r_b$.
> The singularity, where the Kretschmann scalar diverges, is that curve, and beyond it there is no spacetime.
> The rays obey $c\,dt = \pm\partial_r R\,dr$, and outside the cloud the same coordinates are Georges Lemaître's for Schwarzschild's exterior, carried by observers who fall freely from rest at infinity.

### embedding diagram caption: `embedding/tolman_bondi/cloud.caption[1]`

Source: `_tools/derivations/embedding.py` line 3482, `CAPTIONS ("tolman_bondi", "cloud")`.

Before:

> Richard Tolman found these solutions in 1934 and Hermann Bondi took them up in 1947.
> The spacetime diagram draws this cloud marginally bound, $E = 0$, falling from rest at infinity, and then every slice of constant $t$ is flat; drawn here released from rest from the same density at $t = 0$, its slices curve.

After:

> Richard Tolman found these solutions in 1934, and Hermann Bondi took them up in 1947.
> In the spacetime diagram the cloud is marginally bound, $E = 0$, falling from rest at infinity, and then every slice of constant $t$ is flat; released from rest from the same density at $t = 0$, as it is here, its slices curve.

### embedding diagram note: `embedding/tolman_bondi/cloud.input`

Source: `_tools/derivations/embedding.py` line 2240, `in tolman_bondi() input=`.

Before:

> The cloud the spacetime diagram draws, its density falling as $1 - r^2/r_b^2$ to zero at $r_b$ with $R(r, 0) = r$, but released from rest, $E = -GM(r)/c^2r$, each shell falling on its own cycloid, checked to solve this spacetime's own $G^r{}_r = 0$ and to give its density.

After:

> The cloud of the spacetime diagram, its density falling as $1 - r^2/r_b^2$ to zero at $r_b$ with $R(r, 0) = r$, but released from rest, $E = -GM(r)/c^2r$, each shell falling on its own cycloid, checked to solve this spacetime's own $G^r{}_r = 0$ and to give its density.

## Tolman-Oppenheimer-Volkoff Relativistic Star

### history: `tov/history[0]`

Source: `MFS/assets/data/metrics/tov.json`, `history`, paragraph 1.

Before:

> A star is a truce between gravity, which pulls its matter inward, and pressure, which pushes back.
> Writing that truce down in Newtonian gravity gives the classical equation of stellar structure, and closing it requires one more ingredient: a relation between pressure and density, an equation of state [phillips1999].
> The oldest useful choice is a polytrope, $p = K\rho^{\gamma}$, which reduces the whole problem to the Lane-Emden equation studied since the nineteenth century [lane1870].
> The Tolman-Oppenheimer-Volkoff equation is the version of that truce in general relativity, and it is the polytrope, or some richer equation of state in its place, that turns the equation into an actual star.

After:

> A star is a truce between gravity, which pulls its matter inward, and pressure, which pushes back.
> Writing that truce down in Newtonian gravity gives the classical equation of stellar structure, and closing it requires one more ingredient: a relation between pressure and density, an equation of state [phillips1999].
> The oldest useful choice is a polytrope, $p = K\rho^{\gamma}$, which reduces the whole problem to the Lane-Emden equation studied since the nineteenth century [lane1870].
> The Tolman-Oppenheimer-Volkoff equation is the version of that truce in general relativity, and a polytrope, or some richer equation of state in its place, turns the equation into a model of a star.

### history: `tov/history[1]`

Source: `MFS/assets/data/metrics/tov.json`, `history`, paragraph 2.

Before:

> The relativistic story has a sharp lineage.
> In 1931 Subrahmanyan Chandrasekhar fed a relativistically degenerate electron gas, itself a polytrope, into the Newtonian equations and discovered that white dwarfs possess a maximum mass, above which electron degeneracy pressure simply cannot hold [chandrasekhar1931].
> Lev Landau reached a limit of about 1.5 solar masses independently, in a paper he wrote in Zurich in February 1931 and published the next year, and guessed that heavier stars hold regions where "atomic nuclei come in close contact, forming one gigantic nucleus" [landau1932, yakovlev2013].
> Walter Baade and Fritz Zwicky went further in 1934, proposing that a supernova marks "the transition of an ordinary star into a neutron star, consisting mainly of neutrons" [baade1934].

After:

> The relativistic story has a sharp lineage.
> In 1931 Subrahmanyan Chandrasekhar fed a relativistically degenerate electron gas, itself a polytrope, into the Newtonian equations and discovered that white dwarfs possess a maximum mass, above which electron degeneracy pressure cannot hold the star up [chandrasekhar1931].
> Lev Landau reached a limit of about 1.5 solar masses independently, in a paper he wrote in Zurich in February 1931 and published the next year, and guessed that heavier stars hold regions where "atomic nuclei come in close contact, forming one gigantic nucleus" [landau1932, yakovlev2013].
> Walter Baade and Fritz Zwicky went further in 1934, proposing that a supernova marks "the transition of an ordinary star into a neutron star, consisting mainly of neutrons" [baade1934].

### history: `tov/history[2]`

Source: `MFS/assets/data/metrics/tov.json`, `history`, paragraph 3.

Before:

> Eight years after Chandrasekhar, Richard Tolman gave a systematic method for static spheres of perfect fluid in full general relativity [tolman1939], and in the immediately following paper J. Robert Oppenheimer and George Volkoff applied the relativistic equilibrium equation to a cold gas of neutrons, deriving what is now called the TOV equation and finding that neutron cores, too, have a maximum mass [oppenheimer1939].
> Above three quarters of a solar mass they found no static solution at all, and concluded that a massive enough star, its nuclear fuel spent, will "contract indefinitely, although more and more slowly, never reaching true equilibrium" [oppenheimer1939].
> Their estimate was crude, but the conclusion was momentous: beyond it, nothing known could stop collapse.

After:

> Eight years after Chandrasekhar, Richard Tolman gave a systematic method for static spheres of perfect fluid in full general relativity [tolman1939], and in the immediately following paper J. Robert Oppenheimer and George Volkoff applied the relativistic equilibrium equation to a cold gas of neutrons, deriving what is now called the TOV equation and finding that neutron cores, too, have a maximum mass [oppenheimer1939].
> Above three quarters of a solar mass they found no static solution, and concluded that a massive enough star, its nuclear fuel spent, will "contract indefinitely, although more and more slowly, never reaching true equilibrium" [oppenheimer1939].
> Their estimate was crude, but the conclusion was momentous: beyond it, nothing known could stop collapse.

### history: `tov/history[3]`

Source: `MFS/assets/data/metrics/tov.json`, `history`, paragraph 4.

Before:

> This is why the TOV equation, fed a polytrope, is the workhorse of the physics of compact stars.
> By itself the equation is underdetermined: it relates pressure, density, and the enclosed mass, but cannot say how pressure responds to compression.
> Supply that closure, whether an idealized polytrope or a detailed nuclear equation of state, and the TOV equation returns a family of stars with one parameter, a curve of mass against radius, and a maximum mass past which static equilibrium is impossible and a black hole must form.
> Clifford Rhoades and Remo Ruffini showed in 1974 that causality and general relativity alone cap that maximum for a neutron star at 3.2 solar masses, whatever the equation of state does where it is unknown, which gave observers a way to tell a neutron star from a black hole by weighing it [rhoades1974].

After:

> The TOV equation, fed a polytrope, is therefore the workhorse of the physics of compact stars.
> By itself the equation is underdetermined: it relates pressure, density, and the enclosed mass, but cannot say how pressure responds to compression.
> Supply that closure, whether an idealized polytrope or a detailed nuclear equation of state, and the TOV equation returns a family of stars with one parameter, a curve of mass against radius, and a maximum mass past which static equilibrium is impossible and a black hole must form.
> Clifford Rhoades and Remo Ruffini showed in 1974 that causality and general relativity alone cap that maximum for a neutron star at 3.2 solar masses, whatever the equation of state does where it is unknown, which gave observers a way to tell a neutron star from a black hole by weighing it [rhoades1974].

### history: `tov/history[5]`

Source: `MFS/assets/data/metrics/tov.json`, `history`, paragraph 6.

Before:

> Timing the delay of a pulsar's pulses as they pass its companion, Demorest and colleagues weighed PSR J1614-2230 at $1.97 \pm 0.04$ solar masses in 2010, which ruled out almost every proposed equation of state with hyperons or boson condensates [demorest2010]; Antoniadis and colleagues found $2.01 \pm 0.04$ solar masses for a pulsar in a 2.46 hour orbit about a white dwarf in 2013 [antoniadis2013], and Cromartie and colleagues about 2.14 for PSR J0740+6620 [cromartie2020].
> Two teams modelling the X-rays that the NICER telescope recorded from the hot spots of PSR J0030+0451 each found in 2019 a star of about 1.4 solar masses and a radius of about 13 kilometers [riley2019, miller2019].
> The LIGO and Virgo detectors recorded two neutron stars spiralling together on 17 August 2017 [abbott2017], and the tidal distortion of the stars in their last orbits, together with the demand that the equation of state hold up 1.97 solar masses, gave both radii as $11.9 \pm 1.4$ kilometers [abbott2018].

After:

> Timing the delay of a pulsar's pulses as they pass its companion, Paul Demorest and colleagues weighed PSR J1614-2230 at $1.97 \pm 0.04$ solar masses in 2010, which ruled out almost every proposed equation of state with hyperons or boson condensates [demorest2010]; John Antoniadis and colleagues found $2.01 \pm 0.04$ solar masses for a pulsar in a 2.46 hour orbit about a white dwarf in 2013 [antoniadis2013], and Thankful Cromartie and colleagues about 2.14 for PSR J0740+6620 [cromartie2020].
> Two teams modelling the X-rays that the NICER telescope recorded from the hot spots of PSR J0030+0451 each found in 2019 a star of about 1.4 solar masses and a radius of about 13 kilometers [riley2019, miller2019].
> The LIGO and Virgo detectors recorded two neutron stars spiralling together on 17 August 2017 [abbott2017], and the tidal distortion of the stars in their last orbits, together with the demand that the equation of state hold up 1.97 solar masses, gave both radii as $11.9 \pm 1.4$ kilometers [abbott2018].

### conventions: `tov/convention[0]`

Source: `MFS/assets/data/metrics/tov.json`, `convention`, paragraph 1.

Before:

> Coordinates are $(t, r, \theta, \phi)$ with $t$ carrying dimensions of time and $r$ the areal radius, whose sphere has area $4\pi r^2$, and factors of $c$ and $G$ are kept explicit.
> The line element is the one Tolman, and then Oppenheimer and Volkoff, wrote down in 1939 for a static sphere of fluid: their $e^{\nu}$ is $e^{2\Phi}$, and their $e^{-\lambda}$ is written through the mass function, $e^{-\lambda} = 1 - 2m/r$, as Oppenheimer and Volkoff wrote it with their $u$ in units where $G = c = 1$.
> Both functions are left free.

After:

> We use coordinates $(t, r, \theta, \phi)$, with $t$ carrying dimensions of time and $r$ the areal radius, whose sphere has area $4\pi r^2$, and keep factors of $c$ and $G$ explicit.
> The line element is the one Richard Tolman, and then J. Robert Oppenheimer and George Volkoff, wrote down in 1939 for a static sphere of fluid: their $e^{\nu}$ is $e^{2\Phi}$, and we write their $e^{-\lambda}$ through the mass function, $e^{-\lambda} = 1 - 2m/r$, as Oppenheimer and Volkoff wrote it with their $u$ in units where $G = c = 1$.
> We leave both functions free.

### conventions: `tov/convention[2]`

Source: `MFS/assets/data/metrics/tov.json`, `convention`, paragraph 3.

Before:

> Components are taken in the chart $x^0 = ct$, where a chart component of the metric carries no dimensions; a radial derivative needs no factor of $c$, and each order of it carries one inverse length.
> The dots in the geodesic equations are velocities of that same chart, so $\dot{t}$ means $d(ct)/d\lambda$, and the equations are $\ddot{x}^\mu + \Gamma^\mu{}_{\nu\rho}\dot{x}^\nu\dot{x}^\rho = 0$ with $\Gamma$ the Christoffel symbols of that same chart.
> The signature is $(-,+,+,+)$ and the Ricci tensor is the standard contraction $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$, which the Ricci scalar and the Einstein tensor follow.

After:

> We take components in the chart $x^0 = ct$, where a chart component of the metric carries no dimensions; a radial derivative needs no factor of $c$, and each order of it carries one inverse length.
> The dots in the geodesic equations are velocities in the same chart, so $\dot{t}$ means $d(ct)/d\lambda$, and the equations read $\ddot{x}^\mu + \Gamma^\mu{}_{\nu\rho}\dot{x}^\nu\dot{x}^\rho = 0$, with $\Gamma$ the Christoffel symbols of that chart.
> We work in the signature $(-,+,+,+)$, and the Ricci tensor is the standard contraction, $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$, from which the Ricci scalar and the Einstein tensor follow.

### conventions: `tov/convention[3]`

Source: `MFS/assets/data/metrics/tov.json`, `convention`, paragraph 4.

Before:

> No component assumes a field equation.
> Every component is an identity for any pair of functions, and the star is read off the Einstein tensor.
> A perfect fluid at rest has $T^\mu{}_\nu = \mathrm{diag}\left(-\rho c^2, p, p, p\right)$, and the three independent components of $G^\mu{}_\nu = \dfrac{8\pi G}{c^4}T^\mu{}_\nu$ are the equations of stellar structure.

After:

> No component assumes a field equation.
> Every component is an identity for any pair of functions, and the star comes from the Einstein tensor.
> A perfect fluid at rest has $T^\mu{}_\nu = \mathrm{diag}\left(-\rho c^2, p, p, p\right)$, and the three independent components of $G^\mu{}_\nu = \dfrac{8\pi G}{c^4}T^\mu{}_\nu$ are the equations of stellar structure.

### conventions: `tov/convention[4]`

Source: `MFS/assets/data/metrics/tov.json`, `convention`, paragraph 5.

Before:

> The time component, $G^t{}_t = -\dfrac{2\partial_r m}{r^2}$, gives $\partial_r m = \dfrac{4\pi G}{c^2}r^2\rho$: the mass is summed shell by shell over $4\pi r^2dr$ rather than over the proper volume, and the difference between the two sums is the binding energy of the star.
> The radial component, $G^r{}_r = \dfrac{2r\left(r - 2m\right)\partial_r\Phi - 2m}{r^3}$, gives $\partial_r\Phi = \dfrac{m + 4\pi Gr^3p/c^4}{r\left(r - 2m\right)}$, in which the pressure gravitates alongside the mass.
> The angular component asks that the pressure be the same in every direction, $G^\theta{}_\theta = G^r{}_r$, and once the first two have defined $\rho$ and $p$ the Bianchi identity makes that the same statement as $\nabla_\mu T^\mu{}_r = 0$, which is $\partial_r p = -\left(\rho c^2 + p\right)\partial_r\Phi$.
> Put together they are the Tolman-Oppenheimer-Volkoff equation, $\partial_r p = -\dfrac{\left(\rho c^2 + p\right)\left(m + 4\pi Gr^3p/c^4\right)}{r\left(r - 2m\right)}$, which is Newton's $\partial_r p = -GM(r)\rho/r^2$ when $p \ll \rho c^2$ and $m \ll r$.

After:

> The time component, $G^t{}_t = -\dfrac{2\partial_r m}{r^2}$, gives $\partial_r m = \dfrac{4\pi G}{c^2}r^2\rho$: the mass is summed shell by shell over $4\pi r^2dr$ rather than over the proper volume, and the difference between the two sums is the binding energy of the star.
> The radial component, $G^r{}_r = \dfrac{2r\left(r - 2m\right)\partial_r\Phi - 2m}{r^3}$, gives $\partial_r\Phi = \dfrac{m + 4\pi Gr^3p/c^4}{r\left(r - 2m\right)}$, in which the pressure gravitates alongside the mass.
> The angular component requires the pressure to be the same in every direction, $G^\theta{}_\theta = G^r{}_r$, and once the first two have defined $\rho$ and $p$ the Bianchi identity makes that the same statement as $\nabla_\mu T^\mu{}_r = 0$, which is $\partial_r p = -\left(\rho c^2 + p\right)\partial_r\Phi$.
> Together they give the Tolman-Oppenheimer-Volkoff equation, $\partial_r p = -\dfrac{\left(\rho c^2 + p\right)\left(m + 4\pi Gr^3p/c^4\right)}{r\left(r - 2m\right)}$, which is Newton's $\partial_r p = -GM(r)\rho/r^2$ when $p \ll \rho c^2$ and $m \ll r$.

### conventions: `tov/convention[5]`

Source: `MFS/assets/data/metrics/tov.json`, `convention`, paragraph 6.

Before:

> Each of the three corrections it makes to Newton steepens the pressure gradient the star needs, and the pressure that holds the star up is itself a source of gravity, so past some mass squeezing harder only adds weight; that is the maximum mass Oppenheimer and Volkoff found.
> Three equations for the four unknowns $\rho$, $p$, $m$ and $\Phi$ leave one short, and an equation of state closes them.
> A uniform density gives the interior Schwarzschild solution, whose Weyl tensor vanishes, and a polytrope gives stars whose Weyl tensor does not.

After:

> Each of the three corrections the Tolman-Oppenheimer-Volkoff equation makes to Newton's steepens the pressure gradient the star needs, and the pressure that holds the star up is itself a source of gravity, so past some mass squeezing harder only adds weight; this is why Oppenheimer and Volkoff found a maximum mass.
> Three equations for the four unknowns $\rho$, $p$, $m$, and $\Phi$ leave one short, and an equation of state closes them.
> A uniform density gives the interior Schwarzschild solution, whose Weyl tensor vanishes, and a polytrope gives stars whose Weyl tensor does not.

### conformal diagram note: `conformal/tov/spherical.input`

Source: `_tools/derivations/conformal.py` line 2077, `in tov() input=`, an f-string.

Before:

> A polytrope, $p = K\rho_0^2$ with rest mass density $\rho_0$ and energy density $\rho c^2 = \rho_0c^2 + p$, at $K = 100$ and a central $\rho_0 = 1.28\times10^{-3}$ in units where $G = c = M_\odot = 1$, solved from this spacetime's own $G^t{}_t$ and $G^r{}_r$: a star of $M = 1.40\,M_\odot$ and $R = 14.2$ km, the one numerical relativity tests its codes on.

After:

> A polytrope, $p = K\rho_0^2$ with rest mass density $\rho_0$ and energy density $\rho c^2 = \rho_0c^2 + p$, at $K = 100$ and a central $\rho_0 = 1.28\times10^{-3}$ in units where $G = c = M_\odot = 1$, solved from this spacetime's own $G^t{}_t$ and $G^r{}_r$: a star of $M = 1.40\,M_\odot$ and $R = 14.2$ km, the star on which numerical relativists test their codes.

### embedding diagram note: `embedding/tov/star.input`

Source: `_tools/derivations/embedding.py` line 1741, `in tov() input=`, an f-string.

Before:

> A polytrope, $p = K\rho_0^2$ with rest mass density $\rho_0$ and energy density $\rho c^2 = \rho_0c^2 + p$, at $K = 100$ and a central $\rho_0 = 1.28\times10^{-3}$, solved from this spacetime's own $G^t{}_t$ and $G^r{}_r$: a star of $M = 1.40\,M_\odot$ and $R = 14.2$ km, the one numerical relativity tests its codes on.

After:

> A polytrope, $p = K\rho_0^2$ with rest mass density $\rho_0$ and energy density $\rho c^2 = \rho_0c^2 + p$, at $K = 100$ and a central $\rho_0 = 1.28\times10^{-3}$, solved from this spacetime's own $G^t{}_t$ and $G^r{}_r$: a star of $M = 1.40\,M_\odot$ and $R = 14.2$ km, the star on which numerical relativists test their codes.

## Vaidya Radiating Star

### history: `vaidya/history[0]`

Source: `MFS/assets/data/metrics/vaidya.json`, `history`, paragraph 1.

Before:

> Schwarzschild's solution describes the empty space outside a static star, but a star that shines is not surrounded by empty space; it is surrounded by the energy it is pouring out.
> In 1951 the Indian physicist Prahalad Chunilal Vaidya found the metric for exactly this situation: the gravitational field outside a radiating spherical mass [vaidya1951].
> He had come to the problem in 1942 as the private research student of V. V. Narlikar at Banaras, who offered him two outstanding problems, the field of a radiating star and the field of a rotating star; since a rotating star must also radiate, they chose the radiating one [vaidya1997].

After:

> Schwarzschild's solution describes the empty space outside a static star, but a star that shines is surrounded by the energy it pours out.
> In 1951 the Indian physicist Prahalad Chunilal Vaidya found the metric for this situation: the gravitational field outside a radiating spherical mass [vaidya1951].
> He had come to the problem in 1942 as the private research student of V. V. Narlikar at Banaras, who offered him two outstanding problems, the field of a radiating star and the field of a rotating star; since a rotating star must also radiate, they chose the radiating one [vaidya1997].

### history: `vaidya/history[1]`

Source: `MFS/assets/data/metrics/vaidya.json`, `history`, paragraph 2.

Before:

> The trick is disarmingly simple in hindsight.
> Take the Schwarzschild solution and let the mass $M$, instead of being constant, become a function $m(u)$ of retarded time, the time labelling outgoing light fronts, so that as radiation streams away to infinity the mass in the metric steadily decreases.
> By Birkhoff's theorem a truly empty spherical exterior must be static Schwarzschild; Vaidya's is not empty, and so it need not be static [birkhoff1923].
> Vaidya himself reckoned that for any known star the term his solution adds to Schwarzschild's is of order $10^{-20}$, far below anything observable [vaidya1997].

After:

> The trick is disarmingly simple in hindsight.
> Take the Schwarzschild solution and let the mass $M$, instead of being constant, become a function $m(u)$ of retarded time, the time labelling outgoing light fronts, so that as radiation streams away to infinity the mass in the metric steadily decreases.
> By Birkhoff's theorem an empty spherical exterior must be static Schwarzschild; Vaidya's is not empty, and so it need not be static [birkhoff1923].
> Vaidya himself reckoned that for any known star the term his solution adds to Schwarzschild's is of order $10^{-20}$, far below anything observable [vaidya1997].

### history: `vaidya/history[2]`

Source: `MFS/assets/data/metrics/vaidya.json`, `history`, paragraph 3.

Before:

> It is filled with null dust: a flux of pure radiation, energy streaming along light rays with no rest mass of its own.
> The solution comes in two faces, an outgoing version for a star losing energy and an ingoing version for one accreting it.
> In 1965 Richard Lindquist, R. A. Schwartz, and Charles Misner clarified what the dependence on time means physically, showing that $-dm/du$ is precisely the power the star radiates, tied by redshift and Doppler factors to the luminosity a moving observer would actually measure [lindquist1965].

After:

> Vaidya's exterior is filled with null dust: a flux of pure radiation, energy streaming along light rays with no rest mass of its own.
> The solution comes in two faces, an outgoing version for a star losing energy and an ingoing version for one accreting it.
> In 1965 Richard Lindquist, R. A. Schwartz, and Charles Misner clarified what the dependence on time means physically, showing that $-dm/du$ is precisely the power the star radiates, tied by redshift and Doppler factors to the luminosity a moving observer would measure [lindquist1965].

### history: `vaidya/history[3]`

Source: `MFS/assets/data/metrics/vaidya.json`, `history`, paragraph 4.

Before:

> In 1951 Vaidya also published in the Physical Review solutions for spheres of fluid that radiate, with matter and outflowing radiation inside, an expanding zone of pure radiation around them, and empty space beyond [vaidya1951pr].
> When he sent the Astrophysical Journal the paper it printed in 1966, which carries Oppenheimer and Snyder's collapse over to matter that radiates as it falls, its editor, Subrahmanyan Chandrasekhar, answered that it must be published and followed with three pages on how to make it readable to astronomers [vaidya1997, vaidya1966].
> William Kinnersley generalized the shining star in 1969 to a source that accelerates arbitrarily [kinnersley1969], William Bonnor and Vaidya gave its charged version in 1970 [bonnor1970], and Vaidya and L. K. Patel found a radiating Kerr metric in 1973 [vaidya1973].

After:

> Vaidya also published in 1951, in the Physical Review, solutions for spheres of fluid that radiate, with matter and outflowing radiation inside, an expanding zone of pure radiation around them, and empty space beyond [vaidya1951pr].
> When he sent the Astrophysical Journal the paper that appeared there in 1966, which carries Oppenheimer and Snyder's collapse over to matter that radiates as it falls, its editor, Subrahmanyan Chandrasekhar, answered that it must be published and followed with three pages on how to make it readable to astronomers [vaidya1997, vaidya1966].
> William Kinnersley generalized the shining star in 1969 to a source that accelerates arbitrarily [kinnersley1969], William Bonnor and Vaidya gave its charged version in 1970 [bonnor1970], and Vaidya and L. K. Patel found a radiating Kerr metric in 1973 [vaidya1973].

### history: `vaidya/history[5]`

Source: `MFS/assets/data/metrics/vaidya.json`, `history`, paragraph 6.

Before:

> Most consequentially, the metric became a testing ground for cosmic censorship, the conjecture that gravitational collapse always hides its singularities behind horizons.
> Yuhji Kuroda showed in 1984 that when radiation falls in slowly enough, the singularity at the centre can be globally naked, "a counterexample to the cosmic censorship hypothesis" [kuroda1984].
> Achilles Papapetrou showed that the collapse of ingoing null dust can instead produce a naked singularity, visible to the outside universe, one of the first concrete counterexamples to the conjecture [papapetrou1985].
> Dwivedi and Joshi showed in 1989 that the naked singularity of this model is a strong curvature singularity [dwivedi1989].
> A metric built to describe a shining star turned out to be one of the sharpest probes of when general relativity keeps its singularities decently out of sight.

After:

> Most consequentially, the metric became a testing ground for cosmic censorship, the conjecture that gravitational collapse always hides its singularities behind horizons.
> Yuhji Kuroda showed in 1984 that when radiation falls in slowly enough, the singularity at the centre can be globally naked, "a counterexample to the cosmic censorship hypothesis" [kuroda1984].
> Achilles Papapetrou showed that the collapse of ingoing null dust can instead produce a naked singularity, visible to the outside universe, one of the first concrete counterexamples to the conjecture [papapetrou1985].
> I. H. Dwivedi and Pankaj Joshi showed in 1989 that the naked singularity of this model is a strong curvature singularity [dwivedi1989].
> A metric built to describe a shining star became one of the sharpest tests of whether collapse hides its singularities behind horizons.

### conventions: `vaidya/convention[0]`

Source: `MFS/assets/data/metrics/vaidya.json`, `convention`, paragraph 1.

Before:

> Coordinates are $(u, r, \theta, \phi)$ in the outgoing chart and $(v, r, \theta, \phi)$ in the ingoing one, with $u$ the retarded time and $v$ the advanced time, both carrying dimensions of time, and $r$ the areal radius, so that a sphere of coordinate radius $r$ has area $4\pi r^2$.
> Factors of $G$ and $c$ are kept explicit; the mass function $m$ depends on the null time alone, and holding it constant returns Schwarzschild with $r_s = 2Gm/c^2$.
> Components are taken in the chart $x^0 = cu$, or $x^0 = cv$, so a dot is a derivative with respect to that chart coordinate: $\dot{m} = c^{-1}dm/du$ is a mass per unit length, and the dots in the geodesic equations are $d(cu)/d\lambda$ and the like, which is what makes every term of an equation carry the dimensions of its left hand side.

After:

> We use coordinates $(u, r, \theta, \phi)$ in the outgoing chart and $(v, r, \theta, \phi)$ in the ingoing one, with $u$ the retarded time and $v$ the advanced time, both carrying dimensions of time, and $r$ the areal radius, so that a sphere of coordinate radius $r$ has area $4\pi r^2$.
> We keep factors of $G$ and $c$ explicit; the mass function $m$ depends on the null time alone, and holding it constant returns Schwarzschild with $r_s = 2Gm/c^2$.
> We take components in the chart $x^0 = cu$, or $x^0 = cv$, so a dot is a derivative with respect to that chart coordinate: $\dot{m} = c^{-1}dm/du$ is a mass per unit length, and the dots in the geodesic equations are $d(cu)/d\lambda$ and the like, so every term of an equation carries the dimensions of its left hand side.

### conventions: `vaidya/convention[1]`

Source: `MFS/assets/data/metrics/vaidya.json`, `convention`, paragraph 2.

Before:

> The luminosity the star radiates to infinity is $L = -c^2 dm/du = -c^3\dot{m}$.
> The Ricci tensor is the standard contraction $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$, so the field equations read $G_{\mu\nu} = 8\pi G\,T_{\mu\nu}/c^4$.
> This spacetime is not a vacuum: its stress energy is pure radiation, $T_{\mu\nu} = -\dfrac{c^2\dot{m}}{4\pi r^2}k_\mu k_\nu$ in the outgoing chart and $T_{\mu\nu} = +\dfrac{c^2\dot{m}}{4\pi r^2}k_\mu k_\nu$ in the ingoing one, with $k_\mu = -\partial_\mu(cu)$ or $k_\mu = -\partial_\mu(cv)$ the null one form the radiation travels along, so the Ricci and Einstein tensors vanish only where $\dot{m}$ does.

After:

> The star radiates to infinity the luminosity $L = -c^2 dm/du = -c^3\dot{m}$.
> The Ricci tensor is the standard contraction, $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$, so the field equations read $G_{\mu\nu} = 8\pi G\,T_{\mu\nu}/c^4$.
> The spacetime is not a vacuum: its stress energy is pure radiation, $T_{\mu\nu} = -\dfrac{c^2\dot{m}}{4\pi r^2}k_\mu k_\nu$ in the outgoing chart and $T_{\mu\nu} = +\dfrac{c^2\dot{m}}{4\pi r^2}k_\mu k_\nu$ in the ingoing one, with $k_\mu = -\partial_\mu(cu)$ or $k_\mu = -\partial_\mu(cv)$ the null one form along which the radiation travels, so the Ricci and Einstein tensors vanish only where $\dot{m}$ does.

### embedding diagram note: `embedding/vaidya/shell.input`

Source: `_tools/derivations/embedding.py` line 2022, `in vaidya() input=`.

Before:

> An imploding shell of radiation, $m = 0$ for $v < 0$ and $m = M$ for $v > 0$, as the conformal diagram draws it.

After:

> An imploding shell of radiation, $m = 0$ for $v < 0$ and $m = M$ for $v > 0$, as in the conformal diagram.
