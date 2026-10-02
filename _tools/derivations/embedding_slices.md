# Where each embedding diagram is cut from

Every embedding diagram draws part of one moment of its spacetime, or of several moments in turn: a surface of one spatial coordinate $x$ and the angle $\phi$ at one value of a time, every other coordinate held fixed.
The moment is a hypersurface of three dimensions, spacelike everywhere it is drawn, save the wave fronts of the pp-wave and of the Aichelburg-Sexl shock, which are null.
The surface drawn is its equator, or for the flat slices a plane through it, and each point of a spherically symmetric moment's line on a conformal diagram stands for a sphere whose equator is one circle of the surface.

On a drawing of a plane of two coordinates, the moment appears where the plane meets it.
For a plane of $t$ and $r$ at $\theta = \pi/2$ and $\phi = 0$ and a moment of constant $t$, that is the line $t = t_0$, and it is the profile of the embedding at $\phi = 0$ itself.
For a plane that is not in the embedding's surface, as the axis of Kerr's black hole or the Poincaré plane of anti-de Sitter space, it is still the line where the plane crosses the moment, and the moment's equator is the surface drawn.
On a figure in three dimensions of one time and two spatial coordinates, the moment is the surface of constant time, the figure's floor wherever the figure has one.

The moment is drawn over the part of it the embedding's pieces reach, since every point of that part lies on a circle of the drawn surface.
Each piece holds its coordinate $x$ from its first point to its last, so the reach is read from the file: Schwarzschild's $r$ from $r_s$ to $6\,r_s$ on both sheets, the Malament-Hogarth plane's $s$ from $0.03$ to $1.5$ at $ct = 0$.
A reference piece, as the vacuum paraboloid under a star, is no part of the moment and adds nothing to the reach, save for the ideal cosmic string, whose own moment is the reference cone together with the sheet outside Gott's core.
Where the slice is homogeneous and the unit of the embedding is any length, as for Minkowski space, Kasner's and Bianchi's planes and the fronts of the pp-wave and of the Aichelburg-Sexl shock, the reach is the whole drawing.
The moment is then cut to the drawing's box, and where it lies on an edge of the box it is drawn on that edge.

Each moment appears on a drawing in one of six ways: a line, a curve, a point, a surface, the whole drawing, or not at all, and the last always has a reason in the physics.

---

## Spacetime by spacetime

A flat view is named by its keys in the diagram file, `system/view`, a figure in three dimensions by its id, and a conformal view by its id, and every number is in its drawing's own units.

### Aichelburg-Sexl

- The embedding is the wave front $x$, $y$ of four values of $u = ct - z$, $u = -1$, $0.5$, $1$ and $1.5$ in units of $8GE/c^4$.
  The moment is the null hypersurface $u = u_k$, every $v$ on it carrying the same flat front.
- `null_cartesian/half`, `/eighth` and `/thirtysecond`, the planes $x = \rho_0/2$, $\rho_0/8$ and $\rho_0/32$ at $y = 0$ drawn in $z$ and $ct$ with $u = ct - z$ and $v = ct + z$: four null lines $ct = z + u_k$ at 45°, parallel to the shock $u = 0$; the first, $u = -1$, runs before the shock and the other three behind it.
- `shock`, the conformal diagram of the plane $x = \rho_0/8$, $y = 0$: the four lines of constant $p = \arctan u_k$, from $\mathscr{I}^-$ to $\mathscr{I}^+$, parallel to the shock.

### Bonnor's beam of light

- The embedding is the wave front $x$, $y$ of six values of Bonnor's $u = (ct - z)/\sqrt{2}$, $u = 0$, $2$, $4$, $6$, $8$ and $10$ in units of the beam's radius $R$, for the single uniform beam.
  The moment is the null hypersurface $u = u_k$, every $v$ on it carrying the same flat front.
- `cartesian/axis` and `/beside`, the planes $x = 0$ and $x = 3R$ at $y = 0$ drawn in $z$ and $ct$: six null lines $ct = z + \sqrt{2}\,u_k$ at 45°, each a ray moving with the beam.
- `null_cylindrical_interior/edge`, `null_cylindrical_exterior/twice` and `/four`, the planes $\rho = R$, $2R$ and $4R$ drawn in $z$ and $ct$: the same six lines, $u = u_k$.
- `null_cartesian/midway` and `/one` draw two beams side by side, another spacetime than the single beam embedded, and carry no slice.
- The figure `lens`, the plane $y = 0$ with $t$ left out: every wave front covers the whole of it, so no slice is drawn.
- `cartesian_axis` and `interior_axis`, the conformal diagrams of the beam's axis: the six lines of constant $p = \arctan(\sqrt{2}\,u_k/\ell)$, $\ell = 4R$, from $\mathscr{I}^-$ to $\mathscr{I}^+$; `midway` is the two beams' and carries none.

### Alcubierre

- The embedding is the plane $z = 0$ of the Cartesian coordinates at $t = 0$, when the bubble is centred on $x = 0$, a disc of radius $3R$.
- `cartesian/tx`, the plane of $t$ and $x$ at $y = z = 0$: the line $ct = 0$ across the whole box, $x$ from $-3R$ to $3R$, which is the disc's diameter along the path.
- The figure `bubble`, the slice $z = 0$ of $t$, $x$ and $y$: its floor, the plane $t = 0$ out to $2.5R$, on which every cone stands.
- No conformal diagram.
- The height plot now being built draws the expansion $\theta$ over the plane of $x$ and $\rho$ at one moment of $t$; that plane is again the plane $z = 0$ at that $t$, so both appearances stay, at whatever $t$ that change declares.

### Anti-de Sitter

- The embedding is the equator $\theta = \pi/2$ of the static global chart at $t = 0$, $r$ from $0$ to $4L$, drawn in Minkowski space; its reference cone is no part of the moment.
- `static_global/radial`: the line $ct = 0$ from $r = 0$ to $4L$, the whole width.
- `static_global/through`, the line through the centre: the line $ct = 0$ from $x = -4L$ to $4L$, the whole width.
- `poincare/tx`, the plane of $t$ and $x$ at $y = 0$ and $z = L$: the line $ct = 0$ across the whole box.
  In the embedding space the static moment $t = 0$ is the Poincaré moment $t = 0$, and all of it lies in the Poincaré patch, since the patch's edge $U - X_3 = \sqrt{L^2 + r^2} - r\cos\theta$ is positive on it.
  On this plane $r^2 = x^2 + x^4/4L^2$, so the box's edge $|x| = 2L$ is $r = 2\sqrt{2}\,L$, inside the reach.
  The plane crosses the equator itself only at $x = 0$, the centre; the rest of the line lies in the moment off its equator.
- `global`, the strip, each point a sphere: the line $T = 0$ from the centre $\sigma = 0$ to $\sigma = \arctan 4$.
- `poincare`, the plane $x = y = 0$: the line $T = 0$ from $\sigma = -\arctan 4$ to $\arctan 4$, the diameter through the centre, which is $z$ from $(\sqrt{17} - 4)L$ to $(\sqrt{17} + 4)L$.

### Bertotti-Robinson

- The embedding has two views at one moment of $t$, drawn at $t = 0$: the equator, $r$ from $b/e$ to $eb$ with $\phi$, and the sphere of $\theta$ and $\phi$ at one event, drawn at $r = b$.
- `static/radial`: the equator is the line $ct = 0$ from $r = b/e$ to $eb$; the sphere is the point $ct = 0$, $r = b$ on it.
- `poincare/tx`: the Poincaré chart has the same $t$ and $x = b^2/r$, which carries one line element onto the other, so the equator is the line $ct = 0$ from $x = b/e$ to $eb$ and the sphere the point $x = b$.
- `static` and `poincare`, the strip of the first factor, each point a sphere of radius $b$: the line $T = 0$ over the same stretch of $r$ and the point at $r = b$, the same on both, since the two charts cover the same wedge.

### Bianchi type I

- The embedding is the plane $y = 0$ at four moments of cosmic time, $c\bar Ht = 0.10$, $0.38$, $1.00$ and $2.00$, counted from the singularity.
- `type_i_cartesian/tx`, whose $ct$ is counted from the singularity on its bottom edge and whose dashed line $a_i = 1$ stands at $c\bar Ht = 0.378$, the second moment: four lines of constant $ct$ across the whole box, the last on its top edge.
- No conformal diagram.

### The cosmic string

- The embedding is the plane $z = 0$ at one moment of $t$, drawn at $t = 0$: the cone outside Gott's core from $r = \ell\tan\chi_0$ to $3\ell$, the core, and under the core the ideal string's cone as a reference.
- No flat views.
- The figure `beam` is the plane $z = 0$ itself seen from overhead with $t$ left out, so the moment is the whole drawing.
- `conical`, the half plane of the ideal string: the line $T = 0$ from the string to $r = 3\ell$, which is the reference cone together with the sheet.
- `gott`, the half plane with the core, in the proper distance $\rho$: the line $T = 0$ from the axis to $\rho = \ell\chi_0 + 3\ell - \ell\tan\chi_0$, through the edge of the core.

### The spinning string

- The embedding is the plane $z = 0$ at the moment $t = 0$ outside the null circle, from $r = r_c$ to $5\,r_c$, read in the circumference radius chart, whose $R = \sqrt{b^2r^2 - a^2}$.
- `outside` in the proper radius, rescaled radius and circumference radius charts, cylinders of $t$ and $\phi$ at a radius the embedding reaches: the whole line $t = 0$.
- `helical/outside`: the circle of constant $t$ is $c\tau = a\tilde\phi/b$, one turn of the helix, from the left edge to the same event on the right edge.
- `inside` in each chart: inside $r_c$ the circles are closed timelike curves and no surface of constant $t$ is a moment of space, so nothing is drawn.
- The figure `tipping`: the floor from $r_c$ out to its edge at $2r_c$.
- The four conformal views: the line $T = 0$ from $r = r_c$ to $5\,r_c$.

### The travelling wave on a string

- The embedding is the surface of constant $u$ and $v$ across the string, the cone out to $r = 3\ell$, at the four moments $u = -1.5\,\ell$, $0$, $\ell$ and $2.5\,\ell$ on $v = 0$, with the ring of free particles on it.
- `null_conical/toward`, `null_conical/away` and `isotropic/beside`, planes of $u$ and $v$ at one place across the string, which each moment meets at one event: four points, $(u_k, 0)$.
- `moving_string/behind` and `moving_string/ahead`: the same events, whose $V = 2(X - A)A' + \int_0^{u_k}A'^2\,du$ on the line $X$, which `slices.checks()` holds to the map between the two charts.
- No conformal diagram.

### de Sitter

- The embedding is the equator of the static chart at $t = 0$: the hemisphere the static chart covers, $r$ from $0$ to $\ell$, and the antipodal observer's hemisphere beyond the horizon.
- `static_spherical/radial`: the line $ct = 0$ from $r = 0$ to the horizon $r = \ell$; beyond it $t$ is no time and the moment does not continue in this chart.
- `static_spherical/through`: the line $ct = 0$ from $x = -\ell$ to $\ell$.
- `flat_slicing/tx`: a curve.
  With $X_0 = \ell\sinh Ht + \tfrac{1}{2}H\rho^2e^{Ht}$ in the embedding space, the static moment is $X_0 = 0$, which in the flat chart is $Ht = -\tfrac{1}{2}\ln(1 + H^2\rho^2/c^2)$, and with $\rho = |x|$ on this plane it runs from $ct = 0$ at $x = 0$ down to $ct = -\tfrac{1}{2}\ln 5 = -0.80\,c/H$ at the box's edges.
  All of the observer's hemisphere lies on it, $r = \rho e^{Ht}$ reaching $\ell$ only as $\rho \to \infty$; the antipodal hemisphere lies outside the flat chart.
- `static` and `flat`, the square, each point a sphere: the line $T = 0$ across the whole square, the observer's hemisphere $\chi \le \pi/2$ and the antipodal one beyond it, the same on both, since both are the whole spacetime.

### The domain wall

- The embedding is the equator $\theta = \pi/2$ of the global chart at five moments, $kct = -1$, $-0.5$, $0$, $0.5$ and $1$, $z$ from $-1/k$ to $1/k$, the centre of one side through the wall to the centre of the other, drawn in Minkowski space.
- `planar/tz`, the plane $x = y = 0$: the moment $kct = t_k$ of the global chart meets it along the line $t = t_k$, every $z$, since both charts put the point of each side at $R = (1 - k|z|)\cosh(kct)$, $cT = (1 - k|z|)\sinh(kct)$ of its inertial chart; five lines across the whole box.
- `inertial/through`, the line through the centre of the side $z < 0$: the moment is the cone $cT = R\tanh(kct)$, so five lines $cT = |x|\tanh(kct_k)$ through the centre, out to the wall at $|x| = \cosh(kct_k)/k$.
- `planar`, `global`, `conformal` and `inertial`, the whole spacetime, each point a 2-sphere: five curves from the centre of one side at $T = 0$ through the wall to the centre of the other, each the image of $cT = R\tanh(kct_k)$ on both sides, the one at $kct = 0$ the line $T = 0$.

### The Kantowski-Sachs cosmologies

- The embedding has two views, two members of one family: the dust universe symmetric in time at five moments of the dust chart, $\eta = 0$, $0.3$, $0.6$, $0.9$ and $1.2$, and the inside of Schwarzschild's horizon at $T = 0.9$, $0.7$, $0.5$, $0.3$ and $0.1\,r_s$, each the equator over $|r| \le 1$, a cylinder that runs on along $r$.
- `dust/etar`: each dust moment is the line $\eta = \eta_k$, every $r$; five lines across the whole box, and no vacuum moment.
- `comoving/tr`, drawn as the same dust from rest: the moment $\eta_k$ is the line $ct = \pi/2 + \eta_k + \sin\eta_k\cos\eta_k$, counted from the first singularity at $b_0 = 1$; five lines across the whole box.
- `schwarzschild_interior/Tr`: each vacuum moment is the line $T = T_k$, every $r$; five lines across the whole box, and no dust moment.
- `dust`, the lens: five curves $\tau = \tau(\eta_k)$ from one end of the axis to the other.
- `vacuum`, Kruskal's square: five curves of constant $T$ across the black hole's triangle, from one end of the horizon to the other, $\tan p\tan q = (1 - T_k/r_s)e^{T_k/r_s}$.

### Einstein-Rosen waves

- The embedding is four moments of the cylindrical chart, $ct = a$, $2a$, $4a$ and $8a$, of the pulse of Weber, Wheeler, and Bonnor at $C = a$ going out from the axis, each the plane $z = 0$ from the axis to $\rho = 10\,a$.
- `cylindrical/radial`, the plane of $t$ and $\rho$: horizontal lines $ct = a$, $2a$, $4a$ and $8a$ from $\rho = 0$ to $10\,a$, inside the box, which runs from $ct = -a$ to $9a$.
- `null/radial`, drawn against $(v - u)/2$ and $(u + v)/2$: the same lines, $u = ct - \rho$ and $v = ct + \rho$ along each.
- `cylindrical` and `null`, the conformal diagrams: four curves $p, q = \arctan((ct \mp \rho)/a)$ from the axis to $\rho = 10\,a$.

### Ellis-Bronnikov

- The embedding is the equator at one moment of $t$, drawn at $t = 0$, $r$ from $-5\ell$ to $5\ell$ through the throat.
- `spherical/radial`: the line $ct = 0$ across the whole box, $-3\ell$ to $3\ell$.
- `spherical`, the diamond: the line $T = 0$ from $r = -5\ell$ to $5\ell$, $X = 2\arctan(r/\ell)$.

### FRW

- The embedding is the equator of the closed universe of dust, $k = +1$, at five moments $ct = 0.18$, $1.23$, $3.14$, $5.05$ and $6.10$, which are $\eta = \pi/3$, $2\pi/3$, $\pi$, $4\pi/3$ and $5\pi/3$, each the near hemisphere and the far one.
- Every flat view, `comoving_spherical/radial`, `comoving_spherical/through` and `conformal_spherical/radial`, draws the flat universe, $k = 0$: not visible, since the moments drawn are the closed universe's, and the flat universe's moments are planes, which the embedding states and does not draw.
- `closed`, the rectangle $0 \le \chi \le \pi$: five lines $T = \eta$ across its whole width, the near hemisphere $\chi \le \pi/2$ and the far one.
- `flat` and `open`: not visible, other universes.

### The global monopole

- The embedding has two views of two spacetimes of one line element: the monopole with no mass at its centre, the equator of a moment of the Barriola-Vilenkin chart, drawn at $t = 0$, $r$ from the centre to $3\ell$; and Letelier's black hole at $r_s = 1$, the equator of the static moment $t = 0$, $r$ from the throat $r_h = r_s/(1 - \Delta) = 1.235\,r_s$ to $6\,r_s$ on both sheets through the bifurcation sphere, both at $\Delta = 0.19$.
- `conical/radial`: the monopole's line $ct = 0$ from the centre to $3\ell$; the black hole is another spacetime and is not marked.
- `static/radial`: the black hole's line $ct = 0$ from $r_h$ to the box's edge at $6\,r_s$; the monopole is not marked.
- `eddington_finkelstein_ingoing/finkelstein` and `/chart`: with $v = ct + r_*$, $r_* = r/(1 - \Delta) + r_s\ln|(1 - \Delta)r/r_s - 1|/(1 - \Delta)^2$, the moment is the curve $v = r_*$ outside $r_h$, which runs off to $v \to -\infty$ at the horizon, as Schwarzschild's does.
- `eddington_finkelstein_outgoing/finkelstein` and `/chart`: the time reverse, $u = -r_*$.
- `static`, `ingoing` and `outgoing`, Kruskal's diagram of the black hole: the line $T = 0$ through the bifurcation point across both exteriors, $r$ from $r_h$ to $6\,r_s$ on each side, the same on all three.
- `conical`, Minkowski's triangle with its centre singular: the monopole's line $T = 0$ from the centre to $r = 3\ell$.

### The dilaton black hole

- The embedding has three views at $r_d = r_s/2$, the equator of the static moment $t = 0$ as each of three metrics on one manifold measures it: the Einstein metric, the string metric of the magnetically charged hole, and the string metric of the electrically charged one, each from the throat $r = r_s$ to $6\,r_s$ on both sheets through the bifurcation sphere.
- `static/radial`: the Einstein metric's line $ct = 0$ from $r_s$ to $6\,r_s$; the string metrics' moments are not marked.
- `eddington_finkelstein_ingoing/finkelstein` and `/chart`: with $v = ct + r_*$, $r_* = r + r_s\ln|r/r_s - 1|$, the Einstein metric's moment is the curve $v = r_*$ outside $r_s$, as Schwarzschild's is.
- `eddington_finkelstein_outgoing/finkelstein` and `/chart`: the time reverse, $u = -r_*$.
- `string_magnetic/radial` and `string_electric/radial`: each string metric's own line $ct = 0$ from $r_s$ to $6\,r_s$.
- `static`, `ingoing` and `outgoing`, Kruskal's diagram: the Einstein metric's line $T = 0$ through the bifurcation point across both exteriors; `string_magnetic` and `string_electric`: the same line, marked as that string metric's moment.

### Gödel

- The embedding is the plane $z = 0$ of the cylindrical chart at $t = 0$ about the world line $r = 0$, $r$ from $0$ to $\operatorname{arcsinh} 2^{-1/4} = 0.764$, where it stops.
- `cartesian/tx`, the plane of Gödel's own $t$ and $x$ at $y = z = 0$: the line $t = 0$ from $x = -1.529$ to $1.529$.
  Gödel's transformation, $e^x = \cosh 2r + \cos\phi\sinh 2r$, $ye^x = \sqrt{2}\sin\phi\sinh 2r$ and $\tan\left(\tfrac{\phi}{2} + \tfrac{t - 2t'}{2\sqrt{2}}\right) = e^{-2r}\tan\tfrac{\phi}{2}$ with $t'$ the cylindrical time, puts the plane $y = 0$ at $\phi = 0$ and $\pi$, where $x = \pm 2r$ and $t = 2t'$, so the moment $t' = 0$ is $t = 0$ out to $|x| = 2 \times 0.764$; the circles of constant $t'$ become closed timelike curves at $|x| = 2\operatorname{arcsinh} 1 = 1.763$, inside the box.
  The published charts are to be pulled back one onto the other before this is drawn.
- `cylindrical/inside`, the cylinder $r = r_c/2 = 0.441$: the line $t = 0$ across its whole width, every $\phi$.
- `cylindrical/beyond`, the cylinder $r = 3r_c/2$: not visible; the circles there are closed timelike curves, no surface of constant $t$ is a moment of space there, and the embedding stops at $0.764$.
- The figure `tipping`: its floor, the plane $t = 0$, out to the circle $r = 0.764$, inside the critical circle $r_c = 0.881$ of the null circles; the floor runs on to $2r_c$.
- No conformal diagram.

### Interior Schwarzschild

- The embedding is the equator at one moment of $t$, drawn at $t = 0$: the star, $r$ from $0$ to $R = 1.5\,r_s$, and the exterior from $R$ to $4\,r_s$; the vacuum paraboloid under the cap is a reference.
- `spherical/radial`, which draws the star alone: the line $t = 0$ across the whole box, $r$ from $0$ to $R$.
- `spherical/through`: the line $t = 0$ from $x = -R$ to $R$.
- `spherical`, the half diamond: the line $T = 0$ from the centre through the surface to $r = 4\,r_s$.

### Kasner

- The embedding is the plane $y = 0$ at $t = 1/4$, $1/2$, $1$ and $2$.
- `cartesian/tx`, the plane $y = z = 0$, and `cartesian/tz`, the plane $x = y = 0$: each meets the plane $y = 0$ of a moment along its whole line of constant $t$, so each shows four lines $ct = 1/4$, $1/2$, $1$ and $2$ across the whole box, the last on its top edge.
- No conformal diagram.

### Kerr

- The embedding is the equator of a moment of constant Boyer-Lindquist $t$, drawn at $t = 0$, $r$ from $r_+ = 1.436\,GM/c^2$ to $8\,GM/c^2$ on both sheets through the bifurcation sphere.
- `boyer_lindquist/radial`, the plane of $t$ and $r$ on the axis $\theta = 0$: the line $ct = 0$ from $r_+$ to $4\,GM/c^2$, the box's edge; the axis crosses the moment there, off its equator.
- `boyer_lindquist/principal`, the principal null rays projected on $t$ and $r$ at $\theta = \pi/2$: the line $ct = 0$ from $r_+$ to the box's edge.
- `boyer_lindquist/above`, the equatorial plane from above with $t$ left out: every moment projects onto the whole plane, so the moment is the whole drawing outside the horizon, $r_+ < r \le 3\,GM/c^2$.
- The figure `dragging`: its floor, the plane $t = 0$ from $r_+$ to $\tfrac{7}{4}r_E$, carrying the ergosurface.
- `axis`, the tower of the symmetry axis: the line $T = 0$ through the outer bifurcation point across the two exteriors, $r$ from $r_+$ to $8\,GM/c^2$ on each side.

### Kerr-Newman

- As Kerr, with $r_+ = 1.624\,GM/c^2$, on `boyer_lindquist/radial`, `principal`, `above`, the figure `dragging` and the tower `axis`.

### The Krasnikov tube

- The embedding is the cross section of constant $t$ and $x$ at $ct = 3$ and $x = 2$, $r$ from $0$ to $2\rho_0$.
- `cylindrical/tx`, the plane of $t$ and $x$ on the axis $r = 0$: the point $ct = 3$, $x = 2$, where the axis passes through the cross section's centre.
- No conformal diagram.
- The height plot now being built draws the tube over the plane of $x$ and $\rho$ at one moment of $t$; on `cylindrical/tx` that moment is the line $ct = t_0$ across the $x$ it reaches, at whatever $t_0$ that change declares, and the point goes.

### Lentz

- No embedding diagram, and nothing else is drawn of Lentz's soliton.

### Malament-Hogarth

- The embedding is the plane $z = 0$ at $ct = -0.7$, $-0.3$, $-0.1$ and $0$, $s$ from the removed event out to $1.5$, the last moment from $s = 0.03$.
- `cartesian/tx`: four lines of constant $ct$ from $x = -1.5$ to $1.5$, the one at $ct = 0$ broken where the embedding stops, $|x| < 0.03$, about the removed event.
- `cartesian`, the triangle, each point a sphere: four curves.
  With $p, q = \arctan(ct \mp r)$ a moment of $t \ne 0$ bows from the axis toward $i^0$, and the moment $ct = 0$ is the line $T = 0$ from $r = 0.03$ to $1.5$, ending short of the removed event on the axis.

### Melvin

- The embedding has two views, each one moment: `universe`, Melvin's plane $z = 0$ at $t = 0$ from the axis to $\rho = 4/B$, and `ernst`, the equator of Ernst's hole at $t = 0$ from the throat $r = r_s$ to $5\,r_s$ on both sheets.
- `cylindrical/radial`, the plane of $t$ and $\rho$: the line $ct = 0$ from $\rho = 0$ to $4/B$, the edge of the box.
- `ernst/radial`, the equator's plane of $t$ and $r$: the line $ct = 0$ from $r = r_s$ to $5\,r_s$; the other sheet lies in the other exterior, which the chart does not cover.
- `cylindrical`, the half diamond: the curve $p, q = \arctan(\mp B\rho)$ from the axis to $\rho = 4/B$.
- `ernst`, the hexagon: the line $T = 0$ through the bifurcation point, from $r = 5\,r_s$ in one exterior to $5\,r_s$ in the other.

### Levi-Civita

- The embedding is one moment: the plane $z = 0$ at $t = 0$ in Weyl's coordinates, from $\rho = 1/16$ to $\rho = 4$, at $\sigma = 1/4$ and $C = 1$.
- `weyl/radial`, the plane of $t$ and $\rho$: the line $ct = 0$ from $\rho = 1/16$ to $4$, the edge of the box.
- `kasner/radial`, the plane of $t$ and $r$: the line $ct = 0$ from $r = \rho^\Sigma/\Sigma$ at $\rho = 1/16$ to its value at $\rho = 4$, with $\Sigma = 3/4$, that is from $r = 1/6$ to $3.77$; the Kasner form's $t$ is Weyl's times a constant, so the moment is the same.
- `weyl` and `kasner`, the triangle: the curve $p, q = \arctan(\mp 4\rho^{1/4})$ between the same two circles.

### Damour-Solodukhin wormhole

- The embedding is Flamm's paraboloid, the equator at one moment of $t$, drawn at $t = 0$ at $r_s = 1$, $r$ from the throat $r_s$ to $6\,r_s$ on both sides.
- `spherical/radial` and `rescaled/radial`, one side of the throat: the line $ct = 0$ from $r = r_s$ to the box's edge at $6\,r_s$, which the embedding also reaches; a moment of one time is a moment of the other.
- `throat/radial`, both sides in $\rho$ with $\cosh\rho = 52r/r_s - 51$ at $\lambda = 1/5$: the line $ct = 0$ from $\rho = -6.26$ to $6.26$, cut to the box at $|\rho| = 6$.
- `isotropic/radial`, both sides in the isotropic radius: the line $ct = 0$ from $r = 0.0114\,r_s$, the sphere of areal radius $6\,r_s$ on the other side, to $5.49\,r_s$, cut to the box at $2\,r_s$.
- `einstein_rosen/radial`, both sides in $u = \pm\sqrt{r - r_s}$: the line $ct = 0$ from $u = -\sqrt{5}$ to $\sqrt{5}$, cut to the box at $|u| = 2$.
- All five conformal views, the diamond in the tortoise coordinate: the line $T = 0$ through the throat from $6\,r_s$ on one side to $6\,r_s$ on the other, the same on each.

### Simpson-Visser black bounce

Three spacetimes of one line element at $r_s = 1$, each marked on its own drawings alone: the black bounce $a = r_s/2$, the one way wormhole $a = r_s$ and the traversable wormhole $a = 2\,r_s$.

- `outside`, the black bounce's equator at $t = 0$, read in the areal radius from the horizon $\rho = r_s$ to $6\,r_s$ on both sheets.
  On `spherical/bounce` it is the line $ct = 0$ from $r = \sqrt{3}/2$ to the box's edge, on `areal/bounce` the line from $\rho = r_s$ to $6\,r_s$, and on the Eddington-Finkelstein planes the curve $v = r_*$ or $u = -r_*$, which runs off the box toward the horizon.
  On the four conformal views of $a = r_s/2$ it is the line $T = 0$ through the bifurcation point.
- `inside`, five moments of constant $r$ between the horizons, $r = 0.75$, $0.4$, $0$, $-0.4$ and $-0.75\,r_s$, each the stretch $|ct| \le r_s$ of its cylinder.
  On `spherical/bounce` each is the upright segment at $r$ from $ct = -r_s$ to $r_s$, and on the ingoing plane the same segment about $v = r_*(r)$.
  `areal/bounce` is drawn on the side $r > 0$, so it marks the first three, at $\rho = \sqrt{r^2 + a^2}$, the third on its edge $\rho = a$.
  The outgoing chart's region between the horizons is the white hole, so `eddington_finkelstein_outgoing/bounce` marks none of them.
  On the conformal views each is the curve $\tan p\tan q = e^{2\kappa r_*(r)}$ of the black hole, $r = 0$ on the line $T = \pi/2$.
- `null`, the one way wormhole's equator at $t = 0$ from $r = r_s/50$ to $6\,r_s$: the line $ct = 0$ on `spherical/null`, from $\rho = 1.0002\,r_s$ on `areal/null`, the curve $v = r_*$ or $u = -r_*$ on the Eddington-Finkelstein planes, and the line $T = 0$ of one region $r > 0$ on the conformal views.
- `wormhole`, the traversable wormhole's equator at $t = 0$ from $r = -6\,r_s$ to $6\,r_s$: the line $ct = 0$ across `spherical/wormhole`, from the throat $\rho = a$ out on `areal/wormhole`, the curve $v = r_*$ or $u = -r_*$ on the Eddington-Finkelstein planes, and the line $T = 0$ through the throat on the conformal views.

### Fisher-Janis-Newman-Winicour scalar field

- The embedding is the equator at one moment of $t$, drawn at $t = 0$ at $\gamma = 1/2$ and $b = 1$, read in the harmonic chart at $k = 1/2$ from $u = \ln(4/3)$, the sphere $r = 4b$, to $u = 17.5$, where $r - b = 2.5 \times 10^{-8}\,b$.
- `spherical/radial`: the line $ct = 0$ from $r = b$ to $4b$, with $r = b/(1 - e^{-bu})$.
- `jnw/radial`: the same line in $R = r - 3b/4$, from $b/4$ to $13b/4$.
- `isotropic/radial`: the same line in the isotropic radius $(r - b/2 + \sqrt{r(r - b)})/2$, from $b/4$ to $3.48\,b$.
- `harmonic/radial`, drawn in units of $1/k$ at $k = 1$, so at half the embedding's $u$: the line $ct = 0$ from $ku = 0.144$ to $8.75$, cut to the box at $ku = 4$.
- All four conformal views, the triangle: the line $T = 0$ from the singularity on $X = 0$ out to the sphere $r = 4b$, the same on each.

### Kaluza-Klein monopole

The cigar is the surface of $r$ and $x_5$ on the half axis $\theta = 0$ at $t = 0$ in Gross and Perry's chart, from the nut $r = 0$ out to $r = 16m$ at $m = 1$.
Every plane of $t$ and the radius is drawn on that half axis, so the moment is the line $t = 0$ from the nut to $16m$ on each: $r$ itself in Gross and Perry's chart and the Hopf chart, and $\rho = r + 2m$, from $2m$ to $18m$, with the Taub-NUT radius.
On the view through the nut it is the same line, which the page mirrors, and on each conformal view it is the curve $T = 0$ from the centre out to $r_*(16m)$.

### Myers-Perry

- The embedding is four moments, each at $t = 0$ from the throat out on both sheets: with one spin in five dimensions the plane of rotation, $\theta = \pi/2$, and the transverse plane, $\theta = 0$, both from $r_+ = 0.8\sqrt{\mu}$ to $4\sqrt{\mu}$; with equal spins the surface of $\rho$ and the Hopf fibre, from $\rho_+$ to $4\sqrt{\mu}$; and in six dimensions the transverse plane, from $r_+$ to $3\,\mu^{1/3}$.
- `boyer_lindquist/transverse` and `boyer_lindquist/rotation`: the line $ct = 0$ from $r_+$ to the edge of the box, each plane marking the moment embedded on it.
- `ingoing_kerr/transverse` and `ingoing_kerr/rotation`: $v = r_*$ with $r_* = r + \tfrac58\sqrt{\mu}\ln((r - r_+)/(r + r_+))$, the same on both planes, which runs off the drawing toward the horizon.
- `equal_spins/radial` and `boyer_lindquist_six/transverse`: the line $ct = 0$ from the throat out.
- `boyer_lindquist_six/rotation` marks nothing: the embedded transverse plane meets the plane of rotation nowhere outside the horizon, which `HIDDEN` records.
- The conformal views `transverse`, `ingoing`, `rotation`, `equal` and `six`: the line $T = 0$ through the bifurcation surface, from one exterior to the other.
- The four moments are four surfaces, three of them in other spacetimes than the fourth, so each drawing marks its own alone, which `HIDDEN_VIEWS` records.

### Schwarzschild-Tangherlini

- The embedding is two moments, one in each dimension: the plane of $r$ and $\phi$ at $t = 0$ in the static chart, every other angle at $\pi/2$, from the throat $r_h$ to $6\,r_h$ on both sheets.
- `spherical/radial` and `spherical_six/radial`, the plane of $t$ and $r$: the line $ct = 0$ from $r_h$ to the edge of the box at $6\,r_h$, each chart marking its own dimension's moment.
- The Eddington-Finkelstein planes of five dimensions: $v = r_*$ and $u = -r_*$ with $r_* = r + \tfrac{1}{2}r_h\ln((r - r_h)/(r + r_h))$, which runs off the drawing toward the horizon.
- `spherical`, `ingoing`, `outgoing` and `six`, the square: the line $T = 0$ through the bifurcation sphere, from $6\,r_h$ in one exterior to $6\,r_h$ in the other.
- Five dimensions and six are two spacetimes, so no drawing of one marks the other's moment, which `HIDDEN_VIEWS` records.

### Curzon-Chazy

- The embedding is one moment: the plane $z = 0$ at $t = 0$ in Weyl's chart, from $\rho = 0.7226\,m$, where the surface starts, to $5\,m$.
- `weyl/equator` and `spherical/equator`, the plane of $t$ and $\rho$ or $r$: the line $ct = 0$ from $0.7226\,m$ to the edge of the box at $4\,m$.
- `weyl_equator` and `spherical_equator`, the triangle: the curve $p, q = \arctan(\mp\rho_*/\ell)$ over the same stretch.
- `weyl/axis`, `spherical/axis`, `weyl_axis` and `spherical_axis`: nothing, since the plane meets the axis only at $\rho = 0$, inside where the embedding stops, which `HIDDEN` records.

### Zipoy-Voorhees

- The embedding is two moments of two spacetimes: the equatorial plane at $t = 0$ in the spherical chart for the oblate $q = 1$, from $r = 2.5161\,m$ to $6\,m$, and for the prolate $q = -1/2$, from $r = 2.0020\,m$ to $6\,m$.
- `spherical/equator_oblate` and `spherical/equator_prolate`, the plane of $t$ and $r$: the line $ct = 0$ from where that deformation's surface starts to the edge of the box at $6\,m$.
- `prolate_spheroidal/equator_oblate` and `prolate_spheroidal/equator_prolate`, the plane of $t$ and $x$: the same line at $x = r/m - 1$.
- The four equatorial triangles of the conformal diagram: the curve $u, v = \arctan(\mp r_*/\ell)$ over the same stretch.
- Each equatorial drawing marks its own deformation's moment and none of the other's, which the tests list in `HIDDEN_VIEWS`.
- The eight drawings of the axis: nothing, since the equatorial plane does not meet the axis, which `HIDDEN` records.

### Robinson-Trautman

- The embedding is a movie of one wave front, the surface of $\theta$ and $\phi$ at one $u$ and one $r$ of the axisymmetric chart, at $cu = 0$, $0.25$, $0.5$, $1$ and $2\,m$; the fronts of one $u$ differ only in size and each is drawn with its own $r$ as the unit, so a moment is every front of its $u$.
- `axisymmetric/axis` and `axisymmetric/equator`, the planes of $u$ and $r$ drawn against $r$ and $cu + r$: five outgoing rays $u = u_k$, each from $r = 0$ to the edge of the box.
- `axis`, the conformal diagram: the five null lines $p = -\arctan e^{-cu_k/4m}$ from the singularity to future null infinity, the first of them the first front.
### McVittie

- The embedding is a movie of the moments $ct = 1$ to $7\,r_s$ of the isotropic chart, with the moments $ct = 1$, $3$, $5$ and $7\,r_s$ as its surfaces, each from the throat $r = r_s/4a$ out to the comoving $r$ whose areal radius is $6\,r_s$.
- `isotropic/radial`, the plane of $t$ and $r$: the level line $ct$ from the throat to that $r$, or to the edge of the box at $2\,r_s$.
- `areal/radial`, the plane of $t$ and $R$: the same cosmic time, so the level line $ct$ from $R = r_s$ to the edge of the box at $5\,r_s$; `slices.py` checks that the areal chart pulls back onto the isotropic plane and that `mcvittie_areal` is the areal radius.
- `dust_lambda`, the conformal diagram: each moment through the diagram's own labels, from the singular sphere on $T = 0$, where both rays left at that moment's time, out to $R = 6\,r_s$.

### Milne

- The embedding is the equator of the comoving hyperbolic chart at $ct = 0.5$, $1$, $2$ and $3$, each out to $\chi = \operatorname{arcsinh}(4/ct)$, where it is $4$ from the axis, in any length $\ell$; the reference cone is no part of any moment.
- `comoving_hyperbolic/through`: four lines of constant $ct$ from the centre to that $\chi$, the widest two cut by the box's edge $\chi = 3$.
- `comoving_spherical/radial`: four lines of constant $ct$, $r = \sinh\chi$ from $0$ to $4/ct$, the first two cut by the box's edge $r = 4$.
- `logarithmic_time/radial`, at $t_0 = 1$: four lines $c\tau = \ln(ct)$, from $\chi = 0$ to the same $\chi$.
- `inertial/through`: four hyperbolae $c^2T^2 - R^2 = c^2t^2$ across the centre, $|R|$ out to $4$, since $R = ct\sinh\chi$.
- `comoving_hyperbolic`, `comoving_spherical`, `logarithmic_time` and `inertial`, the wedge of Minkowski's triangle: four curves, $p, q = \arctan(ct\,e^{\mp\chi})$, from the centre toward the corner where the light cone of the event $T = R = 0$ meets $\mathscr{I}^+$.

### Minkowski

- The embedding is the equator of the spherical chart at one moment of $t$, drawn at $t = 0$, out to $r = 4\ell$ in any length $\ell$.
- `spherical/radial`: the line $ct = 0$ across the whole box, $r$ from $0$ to $4$.
- `spherical_null/radial`, drawn with $(v - u)/2$ across and $(u + v)/2$ up: the moment is $u + v = 0$, the line $(u + v)/2 = 0$ from $0$ to $4$.
- `cartesian/tx`: the line $ct = 0$ across the whole box, the equator's diameter along $x$.
- `rindler/tx`: $ct = X\sinh(aT/c)$ vanishes only at $T = 0$, so the line $cT = 0$ across the whole box, the half of the moment the wedge holds.
- `spherical` and `spherical_null`, the half diamond: the line $T = 0$ from the centre to $r = 4\ell$.
- `cartesian`, `double_null` and `rindler`, the diamond of the plane $y = z = 0$: the line $T = 0$ from $x = -4\ell$ to $4\ell$, on each, since each is the same plane.

### Mixmaster

- The embedding is the great two sphere of Taub's universe at five moments of its proper time.
- No spacetime diagram and no conformal diagram, so the moments appear nowhere else.

### Morris-Thorne

- The embedding is the equator at one moment of $t$, drawn at $t = 0$, of the Ellis-Bronnikov member, $r$ from the throat $b_0$ to $5b_0$ on both sheets.
- `spherical/radial`, one side of the throat in the areal $r$: the line $ct = 0$ from $r = b_0$ to the box's edge, $4b_0$.
- `spherical` and `proper_radial`, the diamond in $l = \pm\sqrt{r^2 - b_0^2}$: the line $T = 0$ from $l = -\sqrt{24}\,b_0$ to $\sqrt{24}\,b_0$, both sides through the throat, the same on both.

### Natário

- The embedding is the energy density of the riding observers as a height over the plane $z = 0$ at $t = 0$, on a disc of radius $3R$.
- `cartesian_flow/tx`: the line $ct = 0$ across the whole box, as far as the grid's rim reaches, $3R$ either side of the ship.
- No conformal diagram.

### Oppenheimer-Snyder

- The embedding is four moments of the dust's proper time since the release from rest at $R_0 = 2\,r_s$, $c\tau = 0$, $2.48$, $4.01$ and $4.39\,r_s$, which are $\eta = 0$, $3\pi/10$, $3\pi/5$ and $4\pi/5$ of $a = \tfrac{a_m}{2}(1 + \cos\eta)$.
  Inside it is the interior chart's moment of constant $\tau$, $\chi$ from $0$ to $\pi/4$; outside it is Igor Novikov's slice, the moment of clocks released from rest at every areal radius $R$ from $2\,r_s$ to $4\,r_s$ with the dust.
- `interior_comoving/through`, drawn in units of $a_m$: four lines of constant $c\tau/a_m = 0$, $0.876$, $1.418$ and $1.551$ across the whole star, the first on the bottom edge and the last just under the crunch at $\pi/2$.
- `exterior_schwarzschild/radial`, the plane of Schwarzschild's $t$ and $r$ outside the star: curves.
  A shell released from rest at $R$ is at $r = \tfrac{R}{2}(1 + \cos\eta)$ at $\tau = \tfrac{R}{2}\sqrt{R/r_s}(\eta + \sin\eta)/c$, with Schwarzschild's $t$ in the closed form of Misner, Thorne and Wheeler's (31.10), so each moment is the curve of the shells from $R = 2\,r_s$ to $4\,r_s$ at that $\tau$, drawn where $r > r_s$, since the chart ends at the horizon.
  At $\tau = 0$ it is the line $ct = 0$ from $r = 2\,r_s$ to the box's edge at $3\,r_s$, on the bottom edge.
  At $c\tau = 2.48\,r_s$ the star's surface is still outside its horizon, and the curve runs from the surface's world line at $r = 1.59\,r_s$, $ct = 3.82\,r_s$, down to $ct = 3.03\,r_s$ at the box's edge.
  At $4.01$ and $4.39\,r_s$ the inner shells are inside the horizon, and each curve rises toward $t = \infty$ as its shells near $r_s$, leaving the box through its top edge, $ct = 6\,r_s$, near $r = 1.6$ and $1.9\,r_s$.
- `collapse`, the dust and the exterior joined at the star's surface: four curves, each a line $T = \eta$ across the dust and Novikov's curve outside it through the Kruskal coordinates of each shell, met at the surface.

### pp-wave

- The embedding is the wave front $x$, $y$ of four values of the retarded time, $cu = -3$, $-0.5$, $0$ and $0.661\,L$.
  The moment is the null hypersurface $u = u_k$, every $v$ on it carrying the same flat front.
- `exact_plane_wave/tz`, the plane $x = y = 0$ drawn in $z$ and $ct$ with $u = ct - z$ and $v = (ct + z)/2$: four null lines $ct = z + cu_k$ at 45°, the first crossing the box's corner from $z = L$ to $2L$; the drawing is scale free on the axis, so it is read in units of $L$.
- No conformal diagram.

### Reissner-Nordström

- The embedding has two views at one moment of $t$, drawn at $t = 0$, at $r_q = 0.48\,r_s$: outside, $r$ from $r_+ = 0.64\,r_s$ to $6\,r_s$ on both sheets through the outer bifurcation sphere; inside, $r$ from $r_q^2/r_s = 0.2304\,r_s$ to $r_- = 0.36\,r_s$ on both sides of the inner bifurcation sphere.
- `spherical/radial` is drawn today at $r_q = 0.4\,r_s$, with horizons at $0.8$ and $0.2\,r_s$: not visible as it stands, since that is another black hole.
  Redrawn at $r_q = 0.48\,r_s$, the charge its conformal diagram and embedding declare, it shows both views on the line $t = 0$: from $r_+$ to the box's edge at $2\,r_s$, and from $0.2304$ to $0.36\,r_s$, with the region between the horizons, where $t$ is no time, between them.
- `tower`, the maximal extension: the outside view is the line $T = 0$ through the outer bifurcation point across the exterior the chart covers, tinted as its cover, and its mirror, and the inside view is the line through the inner bifurcation point above them, across the two regions inside $r_-$.
- `malament_hogarth`, the same tower with the event beyond the Cauchy horizon marked: the same two lines.

### Bardeen

- The embedding has two views at one moment of $t$, drawn at $t = 0$, at $g = r_s/3$: outside, $r$ from $r_+ = 0.775\,r_s$ to $6\,r_s$ on both sheets through the outer bifurcation sphere; inside, $r$ from the centre to $r_- = 0.301\,r_s$ on both sides of the inner bifurcation sphere.
- `static/radial`: both views on the line $t = 0$, from $r_+$ to the box's edge at $2\,r_s$ and from $0$ to $r_-$, with the region between the horizons, where $t$ is no time, between them.
- The Eddington-Finkelstein views: the moment is $v = r_*$ in the ingoing chart and $u = -r_*$ in the outgoing one, with $r_* = 0$ at the centre, each curve running off toward the horizon its view ends on.
- `tower`, `ingoing` and `outgoing`: the outside view is the line $T = 0$ through the outer bifurcation point, on the outgoing view the line $T = 2\pi$ of the exterior that chart covers, and the inside view is the line $T = \pi$ through the inner bifurcation point, from one centre to the other.

### Schwarzschild

- The embedding is Flamm's paraboloid, the equator of a moment of constant $t$, drawn at $t = 0$, $r$ from $r_s$ to $6\,r_s$ on both sheets through the bifurcation sphere.
- `spherical/radial`: the line $ct = 0$ from $r = r_s$ to the box's edge at $6\,r_s$, which the embedding also reaches; inside $r_s$ the surfaces of constant $t$ are timelike.
- `eddington_finkelstein_ingoing/finkelstein`, drawn against $v - r$: with $v = ct + r + r_s\ln|r/r_s - 1|$ the moment is the curve $v - r = r_s\ln(r/r_s - 1)$, from the bottom edge at $r = (1 + e^{-3})r_s = 1.050\,r_s$ up to $1.609\,r_s$ at $r = 6\,r_s$; it never reaches the horizon, which the ingoing chart crosses at finite $v$.
- `eddington_finkelstein_ingoing/chart`, drawn against $v$: the curve $v = r + r_s\ln(r/r_s - 1)$, from $v = 0$ at $r = 1.278\,r_s$ to $v = 6\,r_s$ at $r = 4.693\,r_s$.
- `eddington_finkelstein_outgoing/finkelstein` and `/chart`: the time reverse, $u + r = -r_s\ln(r/r_s - 1)$ from $3\,r_s$ at $1.050\,r_s$ down to $-1.609\,r_s$ at $6\,r_s$, and $u = -r - r_s\ln(r/r_s - 1)$ from $0$ at $1.278\,r_s$ to $-6\,r_s$ at $4.693\,r_s$.
- `spherical`, `ingoing` and `outgoing`, Kruskal's diagram with $p = \arctan U$ and $q = \arctan V$: the moment $t = 0$ is $U = -V$, the line $T = 0$ through the bifurcation point across both exteriors, $r$ from $r_s$ to $6\,r_s$ on each side, the same on all three, since all three are the whole spacetime.

### Szekeres

- The embedding is the surface $\theta = \pi/2$ through the equators of the shells of the marginally bound cloud, $f = 0$ and $S'/S = 2r(1 - r^2)$ inside $r_b = 1$, at $ct = -0.4$, $0$, $0.3$ and $0.55\,r_b$, played as a movie.
- `axisymmetric/north` and `axisymmetric/south` draw the two halves of the axis of symmetry, $\theta = 0$ and $\theta = \pi$, which that surface meets only at the centre $r = 0$: not visible.
- No conformal diagram.

### Reissner-Nordström-de Sitter

- The embedding has three views of the lukewarm hole, $r_q = r_s/2$ and $\Lambda = 27/(64\,r_s^2)$: the equator of the static moment $t = 0$ from the throat $r_+ = 2r_s/3$ to the widest circle $r_c = 2\,r_s$ and back to the next throat; the static moment inside $r_-$, from $r = 0.2495\,r_s$ to $r_- = 0.4305\,r_s$ and back; and the moment $\tau = 1/H$ of the cosmological chart, $\rho$ from $r_s/100$ to $2\,r_s$, on which the areal radius is $\rho + r_s/2$.
- `static/radial`: the lines $ct = 0$ from $r_+$ to $r_c$ and from $0.2495\,r_s$ to $r_-$, and the cosmological moment as the curve $cT = F(r) - F(1.298\,r_s)$ of `slices.rnds_static_t` from $r_+$ out; between $r_s/2$ and $r_+$ that moment lies in the white hole, and the static plane reads that range of $r$ as the black hole.
- `eddington_finkelstein_ingoing/finkelstein` and `/chart`: the static moments are $v = r_*$, and the cosmological one is $v = cT + r_*$ between $r_+$ and $r_c$, the part of it the ingoing chart covers.
- `eddington_finkelstein_outgoing/finkelstein` and `/chart`: the static moments are $u = -r_*$, and the cosmological one is $u = cT - r_*$ from $r_s/2 + r_s/100$ out, one curve through $r_+$ and $r_c$, where the logarithms of $T$ and of $r_*$ cancel.
- `cosmological/plane`: the cosmological moment is the line $c\tau = 8r_s/3$; the static $t = 0$ between the horizons is $H\tau = e^{-HT_1(r)/c}$, $\rho = (r - r_s/2)/H\tau$, with $T_1$ the static time at $H\tau = 1$, and the moment inside $r_-$, where $\tau < 0$, is the one through $H\tau = -1/2$ on $r = 0.35\,r_s$, since the static time there is fixed only up to a constant.
- `static`, `ingoing`, `outgoing` and `cosmological`, the maximal extension: the static moment between the horizons is the line $T = 0$ from $X = 0$ to $X = 2\pi$ on all four; the moment inside $r_-$ is the line $T = \pi$ above the black hole on the static and ingoing views and $T = -\pi$ below the white hole on the outgoing one, and on the cosmological view the curve of constant static time $cT = $ `slices.RNDS_INSIDE_T` below the white hole, through the inner bifurcation sphere; the cosmological moment is one curve through the white hole, the static region and the expanding region on all four.

### The black string

- The embedding has three views: across the string, the equator of the static moment $t = 0$ at $z = 0$, Flamm's paraboloid from the throat $r_s$ to $6\,r_s$ on both sheets; along the string, the rippling horizon at the advanced times $v = 0$, $14$, $28$ and $42\,r_s$; and across the string of six dimensions, the catenoid from $r_h$ to $6\,r_h$.
- `static/radial`: the line $ct = 0$ from $r_s$ to $6\,r_s$.
- `eddington_finkelstein_ingoing/finkelstein` and `/chart`: the curve $v = r + r_s\ln(r/r_s - 1)$, Schwarzschild's.
- `kerr_schild/radial`: the curve $cT = v - r = r_s\ln(r/r_s - 1)$.
- `static_six/radial`: the line $ct = 0$ from $r_h$ to $6\,r_h$, the moment of six dimensions alone.
- `static`, `ingoing` and `kerr_schild`, the maximal extension: the line $T = 0$ through the bifurcation point from one exterior to the other; `six` the same for the moment of six dimensions.
- The rippled horizon is the perturbed string, another spacetime than the one whose planes are drawn: not visible on any of them, and five dimensions and six are marked each on its own charts.

### Schwarzschild-de Sitter

- The embedding is the equator of Kottler's static moment $t = 0$ at $\Lambda = 0.2/r_s^2$, $r$ from the throat $r_h = 1.085\,r_s$ to the widest circle $r_c = 3.215\,r_s$ and back to the next throat, through the cosmological bifurcation sphere into the next static region.
- `static/radial`: the line $ct = 0$ from $r_h$ to $r_c$; inside $r_h$ and beyond $r_c$ the surfaces of constant $t$ are timelike.
- `eddington_finkelstein_ingoing/finkelstein` and `/chart`: with $v = ct + r_*$ and $r_*(0) = 0$ the moment is the curve $v = r_*$ between the horizons, which runs off to $v \to -\infty$ at $r_h$ and to $v \to +\infty$ at $r_c$, since the ingoing chart crosses only the future half of the black hole horizon and the past half of the cosmological one.
- `eddington_finkelstein_outgoing/finkelstein` and `/chart`: the time reverse, $u = -r_*$.
- `static`, `ingoing` and `outgoing`, the chain of the maximal extension: the moment is the line $T = 0$ from the black hole's bifurcation point at $X = 0$ through the static region, the cosmological bifurcation point at $X = \pi$ and the next static region to the next black hole's bifurcation point at $X = 2\pi$, the same on all three.

### van Stockum

- The embedding is the plane $z = 0$ at one moment of $t$, drawn at $t = 0$, $r$ from $0$ to $0.834\,R$, where it stops.
- `cylindrical/inside`, the cylinder $r = R/2$: the line $t = 0$ across its whole width.
- `cylindrical/beyond`, the cylinder $r = 3R/2$: not visible; beyond $r = R$ the circles are closed timelike curves, and the embedding stops at $0.834\,R$.
- The figure `tipping`, drawn with the proper distance from the axis as its radius: its floor, the plane $t = 0$, out to the proper distance of $r = 0.834\,R$.
- No conformal diagram.

### Taub-NUT

- The embedding is the equator of a moment of constant $t$, drawn at $t = 0$, $r$ from the horizon $r_+ = 2.118\,m$ to $8\,m$, one sheet only, since across the horizon lies Taub's cosmology.
- `spherical/radial`: the line $ct = 0$ from $r_+$ to the box's edge at $6\,m$.
- No conformal diagram.

### Thin shell wormhole

- The embedding is the equator at one moment of $t$, drawn at $t = 0$ at $a = 1.25\,r_s$, $r$ from the throat $a$ to $5\,r_s$ on both sides.
- `spherical/radial`, one side of the throat: the line $ct = 0$ from $r = a$ to the box's edge, $4.25\,r_s$.
- `throat/radial`, both sides through the throat with $\ell = \pm(r - a)$: the line $ct = 0$ from $\ell = -3.75\,r_s$ to $3.75\,r_s$, cut to the box at $|\ell| = 3\,r_s$.
- `throat` and `spherical`, the diamond in $\ell_*$: the line $T = 0$ from $\ell = -3.75\,r_s$ to $3.75\,r_s$, both sides through the throat, the same on both.

### Schwarzschild-anti-de Sitter

- The embedding is the equator of the static moment $t = 0$ at $r_s = 2L$, $r$ from the throat $r_h = L$ to $3L$ on both exteriors through the bifurcation sphere, in flat space out to $2^{1/3}L$ and in Minkowski space beyond.
- `static/radial`: the line $ct = 0$ from $r_h$ to the box's edge at $3L$; inside $r_h$ the surfaces of constant $t$ are timelike.
- `eddington_finkelstein_ingoing/finkelstein` and `/chart`: with $v = ct + r_*$ and $r_*(0) = 0$ the moment is the curve $v = r_*$ outside $r_h$, which runs off to $v \to -\infty$ at the horizon and rises to $r_* = 0.331\,L$ at $3L$, short of its limit $R = 0.658\,L$ at the conformal boundary.
- `eddington_finkelstein_outgoing/finkelstein` and `/chart`: the time reverse, $u = -r_*$.
- `static`, `ingoing` and `outgoing`, the conformal diagram: the moment $t = 0$ is $U = -V$, the line $T = 0$ through the bifurcation sphere across both exteriors, $r$ from $L$ to $3L$ on each side, the same on all three.

### Hayward

- The embedding has three views at $\ell = 12m/7\sqrt{7}$: the equator of the static moment $t = 0$ outside $r_+ = 12m/7$, $r$ from $r_+$ to $6\,m$ on both sheets through the outer bifurcation sphere; the same moment inside $r_- = 6m/7$, $r$ from $0$ to $r_-$ on both sides of the inner bifurcation sphere; and, for the forming and evaporating chart, a movie of the slices $v - r = T$ from $T = -6$ to $8\,m_0$, $r$ from $0$ to $6\,m_0$.
- `static/radial`: the line $ct = 0$ from $r_+$ to the box's edge at $4\,m$, and from $0$ to $r_-$, with the region between the horizons, where $t$ is no time, between them.
- `eddington_finkelstein_ingoing/finkelstein` and `/chart`: with $v = ct + r_*$ and $r_*(0) = 0$ the outside moment is the curve $v = r_*$, which runs off to $v \to -\infty$ at $r_+$, and the inside moment the curve $v = r_*$ from the centre, which runs off to $v \to +\infty$ at $r_-$.
- `eddington_finkelstein_outgoing/finkelstein` and `/chart`: the time reverse, $u = -r_*$.
- `evaporating/history`: each of the six moments the view names is the line $v - r = T$, level on the drawing, from $r = 0$ to the box's edge at $6\,m_0$.
- `static`, `ingoing` and `outgoing`, the conformal diagram: the outside moment is the line $T = 0$ through the outer bifurcation point across both exteriors; the inside moment is the line through the inner bifurcation point above them, $T = \pi$, on the static and ingoing views, whose charts cover an inner region there, and the line $T = -\pi$ through the inner bifurcation point below on the outgoing view, whose chart covers an inner region there.
- `history`, the conformal diagram: each moment is the curve $v - r = T$ from the centre on $X = 0$ out to $r = 6\,m_0$.

### The threaded black hole

- The embedding has two views at $b = 0.9$ and $r_s = 1$: the equator of the static chart's $t = 0$, $r$ from $r_s$ to $6\,r_s$ on both sheets through the bifurcation sphere, and the horizon itself, the bifurcation sphere, $\theta$ from $0$ to $\pi$.
- `static/radial` and `wedge/radial`: the line $ct = 0$ from $r = r_s$ to the box's edge at $6\,r_s$, and the bifurcation sphere as the point $ct = 0$, $r = r_s$.
- The Eddington-Finkelstein views: the string enters $g_{\phi\phi}$ alone, so the equator's moment is Schwarzschild's curve on each, $v - r = r_s\ln(r/r_s - 1)$ and its time reverse; the bifurcation sphere lies at $v \to -\infty$ and $u \to +\infty$, off both charts, so those views mark the equator alone, which `HIDDEN_VIEWS` in the tests records.
- `static`, `wedge`, `ingoing` and `outgoing`, Kruskal's diagram: the equator's moment is the line $T = 0$ across both exteriors, and the bifurcation sphere the point $(X, T) = (0, 0)$, on all four.

### Tolman-Bondi

- The embedding is a cloud released from rest, $E = -GM(r)/c^2r$, at $ct = 0$, $0.6$, $1$ and $1.3\,r_b$.
- `comoving_synchronous/collapse` draws the marginally bound cloud, $E = 0$: not visible, since that is another cloud, whose moments are planes, as the embedding states.
- No conformal diagram.

### Tolman-Oppenheimer-Volkoff

- The embedding is the equator at one moment of $t$, drawn at $t = 0$: the star, $r$ from $0$ to $R = 9.59\,GM_\odot/c^2$, and the exterior to $28.76\,GM_\odot/c^2$; the vacuum paraboloid is a reference.
- `spherical/radial`: the line $ct = 0$ across the whole box, $r$ from $0$ to $16$.
- `spherical/through`: the line $ct = 0$ from $x = -16$ to $16$.
- `spherical`, the half diamond: the line $T = 0$ from the centre through the surface to $r = 28.76$.

### Vaidya

- The embedding is four moments of constant $v - r$ in the ingoing chart, $-3$, $-1.5$, $-0.5$ and $1\,r_s$, since a moment of constant $v$ is a light cone; each runs from $r = 0$ to $4\,r_s$, flat inside the shell, which falls along $v = 0$ and so crosses the moment at $r = -(v - r)$, and Flamm's paraboloid moved in by $r_s$ outside it.
- `eddington_finkelstein_ingoing/shell`, drawn against $v - r$: horizontal lines $cv - r = -1.5$, $-0.5$ and $1\,r_s$ from $r = 0$ to $4\,r_s$.
  The first moment, $-3\,r_s$, lies below the box, whose lowest $cv - r$ is $-2.5\,r_s$; the box is to be taken down to $-3\,r_s$, so that the first moment lies on its bottom edge, and the outgoing view's box up to $3\,r_s$, which keeps the two views each other's time reverse.
- `eddington_finkelstein_outgoing/shell` draws the exploding shell, the time reverse: not visible, since its moments of constant $u + r$ belong to that other spacetime.
- `shell`, the conformal diagram of the imploding shell: four curves, each flat Minkowski time $ct = v - r$ inside the shell, through the Minkowski map there, and Schwarzschild's $v - r$ outside it, through Kruskal's map, met at the shell; the last, $v - r = r_s$, lies wholly outside the shell and crosses the horizon to end on the singularity at $r = 0$.

---

## In a line each

| spacetime | flat views | figure | conformal diagram |
| --- | --- | --- | --- |
| Aichelburg-Sexl | four null lines, three times | none | four null lines |
| Bonnor's beam of light | six null lines, five times | none, every front covers it | six null lines, twice |
| Alcubierre | line | floor | none drawn |
| anti-de Sitter | line, line, line | none | line, line |
| Bell-Szekeres | four events, five times; not visible on the regular view, which lies off $\eta = 0$ | none | four events, five times |
| Bertotti-Robinson | line and point, twice | none | line and point, twice |
| Bianchi I | four lines | none | none drawn |
| cosmic string | none | whole drawing | line, line |
| de Sitter | line, line, curve | none | line, line |
| dilaton black hole | line, four curves, line, line | none | line, five times |
| domain wall | five lines, twice | none | five curves, four times |
| Ellis-Bronnikov | line | none | line |
| FRW | not visible, three times | none | five lines on the closed universe |
| global monopole | line, line, four curves | none | line, four times |
| Kantowski-Sachs | five lines, three times | none | five curves, twice |
| Gödel | line, line, not visible | floor | none drawn |
| interior Schwarzschild | line, line | none | line |
| Kasner | four lines, twice | none | none drawn |
| Kerr | line, line, whole drawing | floor | line |
| Kerr-Newman | line, line, whole drawing | floor | line |
| Krasnikov | point, a line after the height plot | none | none drawn |
| Lentz | none drawn | none | none drawn |
| Malament-Hogarth | four lines | none | four curves |
| Milne | four lines, three times, and four hyperbolae | none | four curves, four times |
| Minkowski | line, four times | none | line, five times |
| Mixmaster | none drawn | none | none drawn |
| Morris-Thorne | line | none | line, line |
| Natário | line | none | none drawn |
| Oppenheimer-Snyder | four lines, four curves | none | four curves |
| pp-wave | four null lines | none | none drawn |
| Reissner-Nordström | two lines once redrawn at $r_q = 0.48\,r_s$ | none | two lines, twice |
| Bardeen | two lines, two curves four times | none | two lines, three times |
| Schwarzschild | line, four curves | none | line, three times |
| Schwarzschild-de Sitter | line, four curves | none | line, three times |
| van Stockum | line, not visible | floor | none drawn |
| Szekeres | not visible, twice | none | none drawn |
| Taub-NUT | line | none | none drawn |
| thin shell wormhole | line, line | none | line, line |
| Tolman-Bondi | not visible | none | none drawn |
| TOV | line, line | none | line |
| Vaidya | three lines, four with the box taken to $-3\,r_s$; not visible on the outgoing view | none | four curves |

As drawn in 4de26cc, the Krasnikov tube's slice is a line, Reissner-Nordström's are two lines, and Vaidya's ingoing view carries all four lines.

---

## The data

Each drawing draws its own slices with the map it draws everything else with, in the same run, so a slice cannot drift from the drawing it stands on.

- `_tools/derivations/slices.py`, new, declares every moment: for each embedding view and surface, the coordinate system, the coordinate or expression held constant and its value, the coordinates held fixed on the surface, and the reach of each piece.
  Where a moment is drawn in a second chart, it declares the transformation, as Schwarzschild's $v = ct + r + r_s\ln|r/r_s - 1|$, Gödel's and the flat slicing's, and the script checks each one by pulling the published metric of one chart back onto the other.
- `null_rays.py` writes `slices` into each flat view it draws a moment on, and `projections.py` into each figure: each `{"view", "surface", "label", "lines", "points", "fills"}`, the embedding view's `id`, the index of the surface in its `surfaces`, the moment as TeX, and the geometry in the drawing's own coordinates, the ones its rays and layers use.
- `conformal.py` does the same for each conformal view, in the coordinates of its layers.
- Each file records in `source` a stamp over the embedding file's moments, and `build_mfs_data.py` refuses a drawing whose stamp no longer matches, or whose slice names a view or surface the embedding file lacks, in `--check` and when writing, as it does for every other source.
- A drawing on which a moment is not visible carries no entry for it, and `slices.py` says why in a sentence, which the tests read.
- `_tools/README.md` defines `slices` in all three file contracts, since the application reads them.

As built in 4de26cc, the stamp is on each slice rather than in the file's `source`: a slice's `version` is taken over the embedding view's settings and the surface's label, time and pieces by `embedding_moment_version` in `build_mfs_data.py`.

The tests, which run without sympy, hold every slice drawn to where its embedding says it is, from the numbers in the files alone:

- on every flat view, each point of a slice mapped back through `to_display` and the box lies on its moment: $ct = t_0$ for a moment of constant $t$, $v - r - r_s\ln|r/r_s - 1| = 0$ on Schwarzschild's Eddington-Finkelstein views, $Ht = -\tfrac{1}{2}\ln(1 + H^2x^2/c^2)$ on the flat slicing, $ct - z = cu_k$ on the pp-wave's axis, and Novikov's closed form on Oppenheimer-Snyder's exterior;
- each slice's ends are the reach of its pieces, or the box's edge where the box cuts it first;
- on every conformal view, each point lies on its moment through the view's closed form, as $T = 0$ on Kruskal's diagram and on every diamond and strip at $t = 0$, $T = \eta_k$ on the closed universe, and $p, q = \arctan(ct \mp r)$ at the Malament-Hogarth moments;
- on every figure, the floor's fill lies in its plane of constant time;
- every moment of every embedding view appears on every drawing `slices.py` declares it on, and on no other.

---

## The layout

The three diagrams of a spacetime stand together as one run, in the place the spacetime diagram has now: after the geodesics and before the references.
That is where the mathematics ends and the reading ends, which is the order of the application's buttons, the history, the maths, the new Graphs and the references, and it keeps the chart chosen at the top in force over all three, as it is over every quantity between.
The website's spacetime page has no row of section buttons, so the run has no button of its own and no heading above the three headings it holds.

The order within the run is the spacetime diagram, the conformal diagram, then the embedding diagram.
The two that carry the slice come first, from the chart's own plane to the whole spacetime, and the surface they point to comes last, so the eye reads the line and then the surface cut along it.

At desktop width each drawing stands beside its caption, as every drawing on the page does now, and the three follow one another down the page.
Three abreast would set each drawing near a third of the panel, about 310 pixels at 1440 pixels wide, and a conformal diagram's labels are sized in units of its width, 14 pixels at its composed width of 628, so they would fall to about 7.
At phone width, and wherever a large text size narrows the room, each drawing stands above its caption, as now, and the run is one column.
In print the three follow one another as they do on the screen.

    desktop                                     phone
    +------------------------------------+      +--------------------+
    | spacetime diagram   [chart views]  |      | spacetime diagram  |
    | +--------------+  caption          |      | [views]            |
    | | t=0 ======== |  legend           |      | +----------------+ |
    | |   (green)    |  drawn at ...     |      | | t=0 ========== | |
    | +--------------+                   |      | +----------------+ |
    | conformal diagram   [views]        |      | caption, legend    |
    | +--------------+  caption          |      | conformal diagram  |
    | |  <== T=0 ==> |  legend           |      | +----------------+ |
    | +--------------+                   |      | |  <=== T=0 ===> | |
    | embedding diagram   [views]        |      | +----------------+ |
    | +--------------+  caption          |      | caption, legend    |
    | |   surface    |  legend           |      | embedding diagram  |
    | +--------------+  drawn at ...     |      | +----------------+ |
    +------------------------------------+      | |    surface     | |
                                                | +----------------+ |
                                                | caption, legend    |
                                                +--------------------+

The slice is drawn in `--green-light`, #4db84a, a solid line of 2.6 pixels on the screen, and `--green-dark`, #2a6e2a, in print.
Every other colour of the palette already stands for something on these drawings: the bright cyan `--mfs-cyan-light`, #37dfff, solid, for the rays of one family and for lines of constant $t$, pink for the other family and for lines of constant $r$, gold for cones and light rays, the site's `--cyan`, #2abed9, dashed, for horizons, amber for the surfaces of matter, pale blue for the apparent horizon and pink red for singularities.
Green appears on them only dotted, for Kerr's ergosurface, so a solid green line is the one mark that reads as new at a glance.
A floor or a whole drawing that is the moment is tinted the same green at 12%, as the covered region is tinted the bright cyan, with its rim in the solid green.

Each slice carries its moment as a label beside its right hand end or its point, in the same green: "$t = 0$", or for a sequence the label its surface carries in the embedding diagram, "$ct = 0.18$", "$v - r = -1.5\,r_s$", so a line and its surface are matched by the same words.
Each drawing's legend adds one entry, the green line, "the moment the embedding diagram draws", or "the moments the embedding diagram draws" for a sequence.
Where the embedding has two views, as Reissner-Nordström's and Bertotti-Robinson's, the slices shown are those of the view its buttons have chosen, and the slices follow the buttons as a slider's moment would; in print every view is printed, so every slice is.
The site's sequences play as movies that pass through every moment rather than standing behind a slider, so every moment of a sequence is drawn at once.

As built in 9ce16f6, the legend reads "the slice of the embedding diagram", or "the slices of the embedding diagram" for a sequence, with the moment after it where the drawing does not name it.
Since 4cd25e4 and 3fd8142 a label stands on whichever side `slices.place` chooses, from the right hand end of its line round to a point 85% of the way along, exactly on the edge of what it names, and `_tools/README.md` gives the order.

---

## What waits on the two changes before this one

The font change to `_layouts/mfs.html` touches the stylesheet the run is laid out in, so the layout is built on it once it lands.
The Alcubierre and Krasnikov height plots change `embedding.py`, `turn.js`, `_tools/README.md` and both spacetimes' embedding files: Alcubierre's moment stays the plane $z = 0$ at one $t$, drawn as the line of constant $t$ on its flat view and the floor of its figure, and Krasnikov's becomes the plane of $x$ and $\rho$ at one $t$, so its point on the flat view becomes a line.
Both are drawn from the files as those changes leave them.
Both landed before any slice was drawn, the font in 84167d9 and 3b0ff21 and the height plots in b912642.
