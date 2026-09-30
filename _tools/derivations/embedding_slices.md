# Where each embedding diagram is cut from

Every embedding diagram draws part of one moment of its spacetime, or of several moments in turn: a surface of one spatial coordinate $x$ and the angle $\phi$ at one value of a time, every other coordinate held fixed.
The moment is a hypersurface of three dimensions, spacelike everywhere it is drawn, save the pp-wave's wave fronts, which are null.
The surface drawn is its equator, or for the flat slices a plane through it, and each point of a spherically symmetric moment's line on a conformal diagram stands for a sphere whose equator is one circle of the surface.

On a drawing of a plane of two coordinates, the moment appears where the plane meets it.
For a plane of $t$ and $r$ at $\theta = \pi/2$ and $\phi = 0$ and a moment of constant $t$, that is the line $t = t_0$, and it is the profile of the embedding at $\phi = 0$ itself.
For a plane that is not in the embedding's surface, as the axis of Kerr's black hole or the Poincaré plane of anti-de Sitter space, it is still the line where the plane crosses the moment, and the moment's equator is the surface drawn.
On a figure in three dimensions of one time and two spatial coordinates, the moment is the surface of constant time, the figure's floor wherever the figure has one.

The moment is drawn over the part of it the embedding's pieces reach, since every point of that part lies on a circle of the drawn surface.
Each piece holds its coordinate $x$ from its first point to its last, so the reach is read from the file: Schwarzschild's $r$ from $r_s$ to $6\,r_s$ on both sheets, the Malament-Hogarth plane's $s$ from $0.03$ to $1.5$ at $ct = 0$.
A reference piece, as the vacuum paraboloid under a star, is no part of the moment and adds nothing to the reach, save for the ideal cosmic string, whose own moment is the reference cone together with the sheet outside Gott's core.
Where the slice is homogeneous and the unit of the embedding is any length, as for Minkowski space, Kasner's and Bianchi's planes and the pp-wave's fronts, the reach is the whole drawing.
The moment is then cut to the drawing's box, and where it lies on an edge of the box it is drawn on that edge.

Each moment appears on a drawing in one of six ways: a line, a curve, a point, a surface, the whole drawing, or not at all, and the last always has a reason in the physics.

---

## Spacetime by spacetime

A flat view is named by its keys in the diagram file, `system/view`, a figure in three dimensions by its id, and a conformal view by its id, and every number is in its drawing's own units.

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

### de Sitter

- The embedding is the equator of the static chart at $t = 0$: the hemisphere the static chart covers, $r$ from $0$ to $\ell$, and the antipodal observer's hemisphere beyond the horizon.
- `static_spherical/radial`: the line $ct = 0$ from $r = 0$ to the horizon $r = \ell$; beyond it $t$ is no time and the moment does not continue in this chart.
- `static_spherical/through`: the line $ct = 0$ from $x = -\ell$ to $\ell$.
- `flat_slicing/tx`: a curve.
  With $X_0 = \ell\sinh Ht + \tfrac{1}{2}H\rho^2e^{Ht}$ in the embedding space, the static moment is $X_0 = 0$, which in the flat chart is $Ht = -\tfrac{1}{2}\ln(1 + H^2\rho^2/c^2)$, and with $\rho = |x|$ on this plane it runs from $ct = 0$ at $x = 0$ down to $ct = -\tfrac{1}{2}\ln 5 = -0.80\,c/H$ at the box's edges.
  All of the observer's hemisphere lies on it, $r = \rho e^{Ht}$ reaching $\ell$ only as $\rho \to \infty$; the antipodal hemisphere lies outside the flat chart.
- `static` and `flat`, the square, each point a sphere: the line $T = 0$ across the whole square, the observer's hemisphere $\chi \le \pi/2$ and the antipodal one beyond it, the same on both, since both are the whole spacetime.

### Einstein-Rosen waves

- The embedding is four moments of the cylindrical chart, $ct = a$, $2a$, $4a$ and $8a$, of the pulse of Weber, Wheeler and Bonnor at $C = a$ going out from the axis, each the plane $z = 0$ from the axis to $\rho = 10\,a$.
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

### Schwarzschild

- The embedding is Flamm's paraboloid, the equator of a moment of constant $t$, drawn at $t = 0$, $r$ from $r_s$ to $6\,r_s$ on both sheets through the bifurcation sphere.
- `spherical/radial`: the line $ct = 0$ from $r = r_s$ to the box's edge at $6\,r_s$, which the embedding also reaches; inside $r_s$ the surfaces of constant $t$ are timelike.
- `eddington_finkelstein_ingoing/finkelstein`, drawn against $v - r$: with $v = ct + r + r_s\ln|r/r_s - 1|$ the moment is the curve $v - r = r_s\ln(r/r_s - 1)$, from the bottom edge at $r = (1 + e^{-3})r_s = 1.050\,r_s$ up to $1.609\,r_s$ at $r = 6\,r_s$; it never reaches the horizon, which the ingoing chart crosses at finite $v$.
- `eddington_finkelstein_ingoing/chart`, drawn against $v$: the curve $v = r + r_s\ln(r/r_s - 1)$, from $v = 0$ at $r = 1.278\,r_s$ to $v = 6\,r_s$ at $r = 4.693\,r_s$.
- `eddington_finkelstein_outgoing/finkelstein` and `/chart`: the time reverse, $u + r = -r_s\ln(r/r_s - 1)$ from $3\,r_s$ at $1.050\,r_s$ down to $-1.609\,r_s$ at $6\,r_s$, and $u = -r - r_s\ln(r/r_s - 1)$ from $0$ at $1.278\,r_s$ to $-6\,r_s$ at $4.693\,r_s$.
- `spherical`, `ingoing` and `outgoing`, Kruskal's diagram with $p = \arctan U$ and $q = \arctan V$: the moment $t = 0$ is $U = -V$, the line $T = 0$ through the bifurcation point across both exteriors, $r$ from $r_s$ to $6\,r_s$ on each side, the same on all three, since all three are the whole spacetime.

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
| Alcubierre | line | floor | none drawn |
| anti-de Sitter | line, line, line | none | line, line |
| Bertotti-Robinson | line and point, twice | none | line and point, twice |
| Bianchi I | four lines | none | none drawn |
| cosmic string | none | whole drawing | line, line |
| de Sitter | line, line, curve | none | line, line |
| Ellis-Bronnikov | line | none | line |
| FRW | not visible, three times | none | five lines on the closed universe |
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
| Schwarzschild | line, four curves | none | line, three times |
| Schwarzschild-de Sitter | line, four curves | none | line, three times |
| van Stockum | line, not visible | floor | none drawn |
| Taub-NUT | line | none | none drawn |
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
Every other colour of the palette already stands for something on these drawings: teal for the rays of one family and for lines of constant $t$, pink for the other family and for lines of constant $r$, gold for cones and light rays, cyan for horizons, amber for the surfaces of matter, pale blue for the apparent horizon and pink red for singularities.
Green appears on them only dotted, for Kerr's ergosurface, so a solid green line is the one mark that reads as new at a glance.
A floor or a whole drawing that is the moment is tinted the same green at 12%, as the covered region is tinted teal, with its rim in the solid green.

Each slice carries its moment as a label beside its right hand end or its point, in the same green: "$t = 0$", or for a sequence the label its surface carries in the embedding diagram, "$ct = 0.18$", "$v - r = -1.5\,r_s$", so a line and its surface are matched by the same words.
Each drawing's legend adds one entry, the green line, "the moment the embedding diagram draws", or "the moments the embedding diagram draws" for a sequence.
Where the embedding has two views, as Reissner-Nordström's and Bertotti-Robinson's, the slices shown are those of the view its buttons have chosen, and the slices follow the buttons as a slider's moment would; in print every view is printed, so every slice is.
The site's sequences stand side by side in one figure rather than behind a slider, so every moment of a sequence is drawn at once.

As built in 9ce16f6, the legend reads "the slice of the embedding diagram", or "the slices of the embedding diagram" for a sequence, with the moment after it where the drawing does not name it.
Since 4cd25e4 and 3fd8142 a label stands on whichever side `slices.place` chooses, from the right hand end of its line round to a point 85% of the way along, exactly on the edge of what it names, and `_tools/README.md` gives the order.

---

## What waits on the two changes before this one

The font change to `_layouts/mfs.html` touches the stylesheet the run is laid out in, so the layout is built on it once it lands.
The Alcubierre and Krasnikov height plots change `embedding.py`, `turn.js`, `_tools/README.md` and both spacetimes' embedding files: Alcubierre's moment stays the plane $z = 0$ at one $t$, drawn as the line of constant $t$ on its flat view and the floor of its figure, and Krasnikov's becomes the plane of $x$ and $\rho$ at one $t$, so its point on the flat view becomes a line.
Both are drawn from the files as those changes leave them.
Both landed before any slice was drawn, the font in 84167d9 and 3b0ff21 and the height plots in b912642.
