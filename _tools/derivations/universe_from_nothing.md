# The universe from nothing

The five charts of `universe_from_nothing.json` are written by `print_charts.py --metric universe_from_nothing`, which computes every tensor from each line element with the checker's own `Geometry` and reads every printed value back before it writes.
This note records the source of each chart and what the diagrams draw.

## Step 1. The geometry

`vilenkin1982`, Alexander Vilenkin, "Creation of universes from nothing", Phys. Lett. B 117, 25 (1982), writes the closed Robertson-Walker metric (2) with the scale factor of a vacuum of energy density $\rho_v$, $\dot a^2 + 1 = H^2a^2$ (3), solved by $a = H^{-1}\cosh(Ht)$ (4), with $H = (8\pi G\rho_v/3)^{1/2}$ and $c = 1$.
Changing $t$ to $-it$ gives $-\dot a^2 + 1 = H^2a^2$ and $a = H^{-1}\cos(Ht)$ (5), the four sphere, defined for $|t| < \pi/2H$.
The page writes $\ell = c/H = \sqrt{3/\Lambda}$, so the Lorentzian half is $a = \ell\cosh(ct/\ell)$ for $t \ge 0$ and the Riemannian half $a = \ell\cos(\tau/\ell)$ for $-\pi\ell/2 \le \tau \le 0$, with $ct = -i\tau$.
Both are the surface of radius $\ell$ in flat space of five dimensions, the hyperboloid $-X_0^2 + X_1^2 + \dots + X_4^2 = \ell^2$ and the sphere $X_0^2 + \dots + X_4^2 = \ell^2$, so every local quantity is that of constant curvature: $R_{\mu\nu} = (3/\ell^2)g_{\mu\nu}$, no Weyl tensor, $R = 12/\ell^2$ and $K = 24/\ell^4$.
At $t = 0$ and $\tau = 0$ each half ends on a 3-sphere of radius $\ell$ with $da/dt = da/d\tau = 0$, so the extrinsic curvature of the join vanishes from both sides and nothing is needed there; this is Vilenkin's "finite size ($a = H^{-1}$) and zero 'velocity'".
The saddle point of Hartle and Hawking's no boundary integral for this model is the same geometry, "a four sphere, continued at its equator to de Sitter spacetime", in `feldbrugge2017`'s words.

## Step 2. The charts

With $n$ the unit vector of $\chi$, $\theta$ and $\phi$ on $S^3$, $(X_1, \dots, X_4) = a\,n$ in every chart:

- `closed`: Vilenkin's (2) and (4) in the angle $\chi$, $X_0 = \ell\sinh(ct/\ell)$, $a = \ell\cosh(ct/\ell)$; it is also the global chart (7) to (9) of `spradlin2003`, Spradlin, Strominger and Volovich, "Les Houches lectures on de Sitter space", hep-th/0110007, on its half $t \ge 0$.
- `four_sphere`: Vilenkin's (5), $X_0 = \ell\sin(\tau/\ell)$, $a = \ell\cos(\tau/\ell)$, on the southern half $X_0 \le 0$.
- `scale_factor`: Vilenkin's (3) and its Euclidean form solved for the time: $d\tau = da/\sqrt{1 - a^2/\ell^2}$ below $a = \ell$ and $c\,dt = da/\sqrt{a^2/\ell^2 - 1}$ above it, so one line element, $da^2/(1 - a^2/\ell^2) + a^2\,d\Omega_3^2$, is the four sphere below the join and de Sitter space above it, with $a$ the "particle coordinate" of his barrier; $X_0 = -\sqrt{\ell^2 - a^2}$ on the sphere and $X_0 = \sqrt{a^2 - \ell^2}$ on the hyperboloid.
- `conformal`: `spradlin2003`'s (10) and (11), $\cosh(ct/\ell) = 1/\cos\eta$, on $0 \le \eta < \pi/2$; $X_0 = \ell\tan\eta$, $a = \ell/\cos\eta$. It is the Einstein static universe's metric times $\ell^2/\cos^2\eta$.
- `lapse`: `feldbrugge2017`'s (14) with the lapse $N \to N/a$, $ds^2 = -(N^2/q)\,dt^2 + q\,d\Omega_3^2$ with $q = a^2$, whose equation of motion and constraint (19) make $q$ a quadratic in $t$, their (20). Solved with $\dot q = 0$ at the waist, $q = \ell^2$ there, and the lapse set to $\ell$ with the time a length, it is $q = \ell^2 + c^2t^2$, so $X_0 = ct$ and $ct = \ell\sinh(ct_c/\ell)$ for the proper time $t_c$ of the closed slicing.

`universe_from_nothing_check` holds each chart to the sphere or the hyperboloid pulled back along these, both sides of $a = \ell$ for the scale factor chart, to $R_{\mu\nu} = (3/\ell^2)g_{\mu\nu}$, to a vanishing Weyl tensor, to the signature its convention states, and the four charts that reach the join to reaching it at rest on a 3-sphere of radius $\ell$.

## Step 3. The numbers the page quotes

The Euclidean action of the whole four sphere is $S_E = -3/(8G^2\rho_v)$ with $\hbar = c = 1$, Vilenkin's $-3m_p^4/8\rho_v$ and `vilenkin2002`'s (7); the tunnelling probability of `vilenkin1984` is $\exp(-3/8G^2\rho_v)$, the inverse of $e^{-S_E}$.

## Step 4. The drawings

Every drawing is at $\ell = 1$.
The spacetime diagrams are the plane of the time and $\chi$ at $\theta = \pi/2$, $\phi = 0$ of each Lorentzian chart, from the waist upward, with light rays $\chi \pm \eta$ in conformal time: $\eta = \arctan\sinh(ct/\ell)$ in the closed slicing, $\arctan(ct/\ell)$ in the quadratic gauge, and $\mathrm{arcsec}(a/\ell)$ above the join in the scale factor chart.
The conformal diagram is the upper half of de Sitter's square, $0 \le \eta < \pi/2$, with the waist on its lower edge, where the four sphere is joined on.
The embedding diagram is the slice $\chi = \theta = \pi/2$ of the scale factor chart, $da^2/(1 - a^2/\ell^2) + a^2\,d\phi^2$: below $a = \ell$ it is a hemisphere of radius $\ell$ in flat space, $z = -\sqrt{\ell^2 - a^2}$, and above it a timelike surface in Minkowski space $dX^2 + dY^2 - dZ^2$, the hyperboloid $Z = \sqrt{a^2 - \ell^2}$, the bowl and its skirt meeting on the circle $a = \ell$ with one vertical tangent.
`_tools/test_universe_from_nothing.py` holds the published files to all of this.
