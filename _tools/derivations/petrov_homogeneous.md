# Petrov's homogeneous vacuum: the charts and what the diagrams draw

Petrov's solution of 1962 is the one vacuum field without a cosmological constant whose group of motions has four parameters and carries any event to any other in exactly one way.
This note records where each chart of `petrov_homogeneous.json` comes from, how the two charts are held to one another and to their neighbours, and what the diagrams take from them.
Units here have $c = 1$; the file keeps $c$.

## Step 1. Petrov's chart

Petrov's chapter lists seven line elements, and the seventh, as the scan of the book prints it, is

$$ds^2 = -e^{x^4}\left[\cos(\sqrt{3}x^4)\,(dx^1)^2 - 2\sin(\sqrt{3}x^4)\,dx^1dx^2 - \cos(\sqrt{3}x^4)\,(dx^2)^2\right] - (dx^4)^2 - e^{-2x^4}(dx^3)^2$$

in the signature $(+,-,-,-)$, with an overall constant $k$.
Stephani, Kramer, MacCallum, Hoenselaers, and Herlt give it as their (12.14), and Gibbons and Gielen (2008) relabel it as their (2.1),

$$k^2ds^2 = dr^2 + e^{-2r}dz^2 + e^r\left(\cos\sqrt{3}r\,(d\phi^2 - dt^2) - 2\sin\sqrt{3}r\,d\phi\,dt\right),$$

which is the chart the file states, with the length $\ell = 1/k$ restored: $r$, $z$, $\phi$, and $ct$ are lengths, each over the whole real line, and the phase is $\psi = \sqrt{3}r/\ell$.
`petrov_homogeneous_check` holds it to a vanishing Ricci tensor, to the determinant $-1$, to Gibbons and Gielen's four Killing vectors (2.2), which are independent at every event and close into their algebra (2.3), and to their left-invariant frame (3.11).

## Step 2. The scalars of the curvature

The Kretschmann scalar and the Chern-Pontryagin scalar both vanish, and the cubic scalar $R_{ab}{}^{cd}R_{cd}{}^{ef}R_{ef}{}^{ab}$ is $-48/\ell^6$; `petrov_invariants` computes the three and the check holds them.
So the two quadratic invariants of the Weyl tensor vanish and the cubic one does not, which is Petrov type I with eigenvalues proportional to the three cube roots of unity.
It agrees with Ferrando and Sáez's characterization, a vacuum of type I whose Weyl eigenvalues are constant, and with the Kretschmann scalar of the Lewis class, proportional to $3 - m^2$, at $m = \sqrt{3}$.

## Step 3. The light cones turn

On a plane of constant $r$ and $z$ the metric is $e^{r/\ell}\,\mathrm{Re}\left[e^{i\psi}(d\phi + i\,dt)^2\right]$, a flat Lorentzian plane whose null directions make the angles $\pi/4 - \psi/2$ and $3\pi/4 - \psi/2$ with the $\phi$ axis.
The cones turn through $\psi/2$, and the vector $\cos(\psi/2)\,\partial_t + \sin(\psi/2)\,\partial_\phi$ is timelike of norm $-e^{r/\ell}$ at every $r$: it is Gibbons and Gielen's arrow of time, and the convention takes it for the future.
The flat views are drawn at $\psi = 0$, $\pi/2$, $\pi$, and $2\pi$, with the time functions $t$, $t + \phi$, $\phi$, and $-t$.
Along a null direction at the angle $\theta$, $\Gamma^r{}_{ab}k^ak^b = -e^{r/\ell}\cos(\psi + \pi/3 + 2\theta)/\ell$, which is $\pm\sqrt{3}e^{r/\ell}/2\ell$ on the two null directions, so neither null line of a plane is a null geodesic.

## Step 4. The closed timelike curve

No coordinate is periodic, and the spacetime still holds closed timelike curves, as Gibbons and Gielen state.
`projections.petrov_loop_events` builds one at $\ell = 1$.
With $w = e^{2\pi i/3}$, a curve that moves in $r$ at unit rate and along the cones' axis at the rate $\kappa e^{-r/2}$ has $d(t + i\phi) = \kappa e^{wr}|dr|$ and a tangent of norm $1 - \kappa^2$, so $t + i\phi$ moves by $Z(r) = \kappa e^{wr}/w$, along a logarithmic spiral.
On the plane $r_1 = \pi/\sqrt{3}$ the future is $+\phi$ and on $-r_1$ it is $-\phi$, so the curve runs along $+\phi$ from $-A_1$ to $A_1$ on the first with $dt = -\beta\,d\phi$, crosses, runs along $-\phi$ from $A_3$ to $-A_3$ on the second, and crosses back.
It closes when $A_3 - A_1 = \mathrm{Im}\,I$ and $\beta(A_1 + A_3) = \mathrm{Re}\,I$ with $I = Z(r_1) - Z(-r_1)$.
The figure draws it at $\kappa = 5/4$ and $\beta = 4/5$, where its speed past the observers who move along the cones' axis is $0.8\,c$ on every stretch, and checks every step timelike and future directed against the published metric.

## Step 5. The chart outside the dust

Bonnor (1979) read the solution as the vacuum outside van Stockum's cylinder of dust at $aR = 1$, and Gibbons and Gielen give the functions of the Weyl-Papapetrou form, from Tipler (1974), and the map.
The file writes it as the heavy cylinder of `lewis` at $m = \sqrt{3}$, where $w = 1/R$:

$$ds^2 = -F\,dt^2 - 2M\,dt\,d\phi + \frac{\ell^2}{r^2}(dr^2 + dz^2) + L\,d\phi^2,$$

$$F = \frac{r}{R}\left(\cos\psi - \frac{\sin\psi}{\sqrt{3}}\right), \qquad M = r\left(\cos\psi + \frac{\sin\psi}{\sqrt{3}}\right), \qquad L = -\frac{2Rr\sin\psi}{\sqrt{3}}, \qquad \psi = \sqrt{3}\ln\frac{r}{R}.$$

Any $\ell$ gives a vacuum, and the vacuum joins the dust smoothly at $\ell = R/\sqrt{e}$.
The map is $r_P = \ell\ln(r/R)$, $z_P = \ell z/R$,

$$t_P - \phi_P = 3^{1/4}\,t, \qquad \sqrt{2 + \sqrt{3}}\,\phi_P + \sqrt{2 - \sqrt{3}}\,t_P = 3^{1/4}\sqrt{2}\,R\,\phi,$$

which is Gibbons and Gielen's with the time reversed: theirs has the cross term $+2M$, and with the reversal the two charts share a future on the surface and the dust turns as the published dust of `stockum_dust` does.
`petrov_homogeneous_check` holds the chart to a vacuum, to $FL + M^2 = r^2$, to the published heavy cylinder at $m = \sqrt{3}$, to Petrov's chart carried along the map, and to the published dust and its first derivative on $r = R$.
On the surface $L = 0$, so the circle there is a closed null curve, the image of one null direction of Petrov's plane $r = 0$; $L$ is negative out to $r = e^{\pi/\sqrt{3}}R$, and $F$ is negative from $\psi = \pi/3$ to $4\pi/3$.
At $r = 12R$ both are positive again and the cones have turned through $123°$, so the future is the side of falling $t$; on the surface $g^{tt}$ vanishes with $L$, and the flat view takes $t + R\phi/c$ for its time function.

## Step 6. The embedding and what is not drawn

The plane of $r$ and $z$ at fixed $t$ and $\phi$ has the metric $dr^2 + e^{-2r/\ell}dz^2$, a hyperbolic plane with horocycles for its lines of constant $r$, and it is totally geodesic.
A strip of width $2\pi\ell$ along $z$, rolled up, is half of Beltrami's pseudosphere above $r = 0$ and a sheet in Minkowski space below it, the mirror image of the moment of `lifshitz_spacetime`.
It meets each plane of $t$ and $\phi$ in one event and the figure's slice along its axis, so no moment is marked on them, and `slices.HIDDEN` says so.
No conformal diagram is drawn: the planes of $t$ and $\phi$ are flat and their null lines are not geodesics, and no plane of a time and $r$ holds its own light rays, since $\Gamma^\phi{}_{tr}$ turns them out of it.
