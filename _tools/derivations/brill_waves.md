# Brill's time-symmetric gravitational waves

Brill's waves are one moment of a spacetime: a Riemannian space of three dimensions, signature $(+,+,+)$, with no time in it.
No spacetime around the slice is known in closed form, so the page states the slice and nothing of its evolution.
This note records the three charts with their sources, the solver the drawings use, and what the embedding diagram draws.

## Step 1: one moment

At a moment of time symmetry the extrinsic curvature vanishes, the momentum constraint holds identically, and the Hamiltonian constraint in a vacuum is $R = 0$ (Cook, Living Rev. Relativ. 3, 5, section 2.3).
No coordinate is declared a time in `DIMENSIONS`, as for Misner's and Brill and Lindquist's data, so the checker scales nothing by $c$ and a dot in a geodesic equation is a derivative with respect to arc length.
In three dimensions the Weyl tensor vanishes identically, and the charts print none.
The slice is not conformally flat, since its Cotton tensor does not vanish for a free $q$, so the one owned tag is `three-dimensional`.
There is no light cone, so `null_rays.py` draws nothing, and the slice stands in `NOT_DRAWN` in `conformal.py`.

## Step 2: Brill's chart

$ds^2 = \psi^4\left[e^{2q}(d\rho^2 + dz^2) + \rho^2 d\phi^2\right]$ is Brill's metric (Ann. Phys. 7, 466, 1959), as Ó Murchadha writes it (arXiv:gr-qc/9302023, his (1) with $Aq$ written $q$) and Alcubierre and others (Class. Quantum Grav. 17, 2159, section V).
Brill's paper was not to hand; the formulas are taken from those two and each is checked here.
The base metric $e^{2q}(d\rho^2 + dz^2) + \rho^2 d\phi^2$ has the scalar curvature $-2e^{-2q}(\partial_\rho^2 q + \partial_z^2 q)$, Ó Murchadha's (3), and $\psi^4$ times it has
$$R = -\frac{8\nabla^2\psi + 2\psi\left(\partial_\rho^2 q + \partial_z^2 q\right)}{\psi^5 e^{2q}},$$
with $\nabla^2$ the flat Laplacian, so the constraint is Brill's equation, $\nabla^2\psi + \tfrac{1}{4}\psi(\partial_\rho^2 q + \partial_z^2 q) = 0$, Ó Murchadha's (24).
$\psi$ and $q$ are left free functions of $\rho$ and $z$, as Brill and Lindquist's $\psi$ is, so no component assumes the constraint.
`brill_waves_check` holds the Ricci scalar to that form, the base metric to its curvature, and the chart at $\psi = 1$, $q = 0$ to flat space.

$q$ must vanish on the axis with its first derivative along $\rho$ and fall off faster than $1/r$.
Then $\int(\partial_\rho^2 q + \partial_z^2 q)\rho\,d\rho\,dz = 0$, Ó Murchadha's (4), and with $\psi \to 1 + M/2r$ Green's identity gives
$$2\pi M = \int\frac{(\nabla\psi)^2}{\psi^2}\,dV,$$
his (8), Brill's positive mass, where $M$ is the length $GM/c^2$.
The two shapes the page names are Holz, Miller, Wakano and Wheeler's, $q = a\rho^2e^{-r^2}$, and Eppley's, $q = a\rho^2/(1 + r^n)$ with $n \geq 4$, in units of the width $\lambda$ (Alcubierre and others, (38) and (39)); the check holds both to vanishing on the axis.

## Step 3: the spherical chart and the chart without axial symmetry

Eppley solved for $\psi$ in spherical coordinates, $\rho = r\sin\theta$ and $z = r\cos\theta$, where the metric is $\psi^4\left[e^{2q}(dr^2 + r^2d\theta^2) + r^2\sin^2\theta\,d\phi^2\right]$ (Phys. Rev. D 16, 1609).
The check holds it to being the cylindrical chart pulled back and its Ricci scalar to the same form with $\partial_r^2 q + \partial_r q/r + \partial_\theta^2 q/r^2$ for the source.

Alcubierre, Allen, Brügmann, Lanfermann, Seidel, Suen and Tobias kept Brill's form and let $q$ depend on $\phi$, $q = a\rho^2e^{-r^2}\left[1 + c\,\rho^2\cos^2(n\phi)/(1 + \rho^2)\right]$ (Phys. Rev. D 61, 041501, their (1) and (2)).
The page writes their $c$ as $b$, since $c$ is the speed of light on every page.
The third chart is Brill's with both functions free in all three coordinates; its Ricci scalar holds the terms in $\partial_\phi q$ and $\partial_\phi\psi$, and the check holds it to the cylindrical chart's where neither function depends on $\phi$.

## Step 4: the solver

`brill_wave.py` solves Brill's equation for Holz and others' wave.
In units of $\lambda$ the potential is $V = \tfrac{1}{4}(\partial_\rho^2 q + \partial_z^2 q) = \tfrac{a}{4}e^{-r^2}(2 - 12\rho^2 + 4\rho^2r^2)$, which `brill_waves_check` holds to the stated $q$.
With $\mu = \cos\theta$ it is $V_0(r) + V_2(r)P_2(\mu)$, so $\psi = 1 + \sum_l u_l(r)P_l(\mu)$ over even $l$ obeys ordinary equations coupled through $(2l + 1)/2\int P_lP_2P_{l'}\,d\mu$.
Each $u_l$ is collocated on Chebyshev's points of $[0, 8]$ and continued beyond as the multipole $u_l(8)(8/r)^{l+1}$, since $e^{-64}$ is below $10^{-27}$.
Ninety-six points and sixteen polynomials give the mass to a part in $10^9$ at $a = 12$.

The masses are $0.0339$, $0.1263$, $0.699$, $2.913$ and $4.669$ at $a = 1$, $2$, $5$, $10$ and $12$, where Alcubierre and others' Table V has $0.0338 \pm 0.0004$, $0.1262 \pm 0.0009$, $0.696 \pm 0.003$, $2.912 \pm 0.008$ and $4.67 \pm 0.01$.
Brill's integral of $(\nabla\psi/\psi)^2$ agrees with the coefficient of $1/r$ to a part in $10^{11}$.
The mass is positive for negative $a$ as well, $3.21$ at $a = -5$.

Past $a = 12$ the mass grows fast, $7.63$ at $14$, $13.6$ at $16$ and $31.5$ at $18$, and it diverges just above $a = 20$, beyond which $\psi$ has a node and the coefficient of $1/r$ is negative.
That is the critical amplitude of Cantor and Brill and of Beig and Ó Murchadha, where the potential first holds a bound state at zero energy (Ó Murchadha, sections 3 to 5).
The page states no number for it, since no source gives one for this wave.

## Step 5: the apparent horizon

A surface of revolution is a curve $(\rho(s), z(s))$ of the half plane, and its area is $2\pi\int\psi^4e^q\rho\,ds$, so it is minimal when the curve is a geodesic of $(\psi^4e^q\rho)^2(d\rho^2 + dz^2)$: $d\theta/ds = n\cdot\nabla\ln(\psi^4e^q\rho)$, with $\theta$ the angle of the tangent.
On the axis $q$ vanishes and $d\theta/ds = 2\partial_z\psi/\psi$.
A curve shot from the axis at right angles closes into a surface when it meets the plane $z = 0$ at right angles, and at a moment of time symmetry a minimal surface is marginally trapped, the outermost being the apparent horizon (Gibbons, Commun. Math. Phys. 27, 87).

A pair of such surfaces, an inner and an outer, first appears at $a = 11.816$.
Alcubierre and others bracket it between $11.81$ and $11.82$, where Holz and others had $7.5$.
At $a = 12$ the outer one has the area $1090.6$ and $16\pi M^2/A = 1.0049$, so the Penrose inequality holds; Alcubierre and others found $A \sim 1.1 \times 10^3$ and $0.997$, within their error.
`_tools/test_brill_waves.py` holds all of it.

## Step 6: the embedding diagram

The slice is a moment already, so both views are of the plane $z = 0$, the fixed set of $z \to -z$.
On it $q = a\rho^2e^{-\rho^2}$ and the metric is $\psi^4(e^{2q}d\rho^2 + \rho^2d\phi^2)$: the circle at $\rho$ has the radius $\rho\psi^2$ and $dz/d\rho = \psi\sqrt{\psi^2e^{2q} - (\psi + 2\rho\,\partial_\rho\psi)^2}$.

For $a$ below $10.8$ that is imaginary in a band about $\rho = 2$, where $q$ has died away and $\psi$ is still rising: from $1.59$ to $2.43$ at $a = 4$, narrowing to nothing between $a = 10.5$ and $11$.
So the weaker wave, $a = 4$, which Alcubierre and others saw disperse, is drawn in two pieces, as Black Saturn's plane is, with the two circles where the surface stops set at one height.

For $a$ from $11$ to $13$ the plane is one piece and is played as a movie in $a$.
At $\rho = 1$, where $q = a/e$, the distance between neighbouring circles is stretched by $e^q$, $57$ at $a = 11$ and $119$ at $a = 13$, while the circles there have radii below $\lambda$: the plane is a long narrow stalk, and the surface rises $38$ to $89$ from its centre to the rim at $\rho = 6$.
The first chords at the centre are a ten thousandth of that height, so each frame is rounded to a part in $10^9$ of its extent.
The circle where the apparent horizon cuts the plane is marked from $a = 11.82$.
From $a = 12.19$ the plane's own circles have a widest and then a narrowest, a neck above the stalk; the frames at $12.2$ and above show it.

The charts took 19 seconds to print and the embedding diagram 195 seconds to draw on 2 October 2026.
