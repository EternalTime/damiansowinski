# Voice rewrites of the captions and conventions

Every caption, note, restriction band and sentence stated in place of a drawing on the spacetimes page, and every paragraph of every conventions section, was read against ~/VOICE.md and rewritten here where it strayed from it.
Each entry names the source the text is written from, then gives the text before and after, one sentence to a line; joining the lines of a block with single spaces gives the text exactly.
Every `$...$` of every entry is unchanged character for character and in order, and so is every number written in digits outside it.
Of 511 distinct texts, 158 are rewritten and 353 stand as they were.

## Alcubierre Warp Drive

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
