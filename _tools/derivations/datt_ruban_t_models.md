# The T-models of Datt and Ruban

Dust whose spheres all have one radius at each moment, $ds^2 = -c^2dt^2 + a(t, r)^2dr^2 + b(t)^2d\Omega^2$: the spherical dust solution that Lemaître, Tolman and Bondi's family leaves out, and the Kantowski-Sachs dust universe with its scale factor along the axis free to differ from shell to shell.
This note records the five charts, where each comes from, the one place a source is corrected, and what each diagram draws.

## Step 1: the sources

Every paper was verified before anything was built.
B. Datt, Z. Phys. 108, 314 (1938), doi:10.1007/BF01374951, is on Crossref with the author line "B. Datt"; its English reprint is Gen. Rel. Grav. 31, 1619 (1999).
V. A. Ruban's letter, JETP Lett. 8, 414-417 (1968), was read from the scan in the journal's own archive at `jetpletters.ru`, issue 12 of volume 8, and his paper, Sov. Phys. JETP 29, 1027-1034 (1969), from `jetp.ras.ru`; the reprints are Gen. Rel. Grav. 33, 369 and 375 (2001).
Andrzej Krasiński's editor's notes on Datt (1999) and on Ruban (2001) were read in full from his own page at `users.camk.edu.pl/akr`, the second with Irina Dymnikova's biography of Ruban, which is where his name Vladimir and his life come from.
Krasiński and Giono's paper on Ruban's charged solution, Gen. Rel. Grav. 44, 239 (2012), was read from arXiv:1107.4897.
The letter of 1968 does not cite Datt; the paper of 1969 does, as its reference 4.

## Step 2: the solution, and a factor of 2

With $b$ independent of $r$ the equation $G^r{}_r = 0$ is $2b\,b'' + b'^2 + 1 = 0$, the closed Friedmann universe's, so $b = r_s\sin^2(\eta/2)$ with $c\,dt = b\,d\eta$, and $G^\theta{}_\theta = 0$ is linear in $a$.
Its two solutions are $\cot(\eta/2)$, the inside of Schwarzschild's horizon, and $1 - \tfrac{1}{2}\eta\cot(\eta/2)$.
The pages write $a = \epsilon\cot(\eta/2) + 2\mu\left(1 - \tfrac{1}{2}\eta\cot(\eta/2)\right)$ with $\mu(r) = (G/c^2)\,d\mathcal{M}/dr$, for which $G^\eta{}_\eta = -2\mu/(a\,b^2)$: the density is $8\pi G\rho/c^2 = 2\mu/(a\,b^2)$ and the rest mass between two shells is $\rho\,4\pi b^2a\,dr = (c^2/G)\,\mu\,dr$.
That is Krasiński's form, his (3) and (4) of 1999, with his $X = \mu$ and $Y = \epsilon$.
Ruban's (2) of 1968 and (18) of 1969 print $\mathcal{M}'$ where the field equations need $2\mathcal{M}'$ beside his density $\mathcal{M}'/(4\pi r^2e^{\omega/2})$; his (19), the de Sitter case, is consistent as printed.
`datt_ruban_check` holds the chart to the density, and `_tools/test_datt_ruban_t_models.py` holds the factor: with one $\mu$ in place of two the density does not come out.
Rescaling $r$ multiplies $\epsilon$ and $\mu$ together, which is how $\epsilon$ is brought to $1$, $0$ or $-1$.

With $\epsilon = 1$ the scale factor is $a = (1 - \mu\eta)\cot(\eta/2) + 2\mu$, positive near the bang and positive near the crunch only if $\mu \ge 1/2\pi$; below that $a$ reaches zero before the crunch and neighbouring shells meet.
Ruban's sentence after his (18) gives the threshold as $\mathcal{M}' > \tfrac{1}{2}\pi$ in the scan's typography, which is $1/2\pi$ in the normalisation the field equations fix.
Krasiński and Giono's (6.6), $0 < Y < \pi X/2$, is not used: their (6.2) changes the sign of $Y$ between the expansion and the collapse, which leaves $a$ with a corner at the greatest expansion, while the solution in $\eta$ is one analytic function across it.
What is taken from them is that shell crossings can be avoided without charge and cannot with it, and that the surface of a T-sphere stays inside the horizon and touches it once.

## Step 3: the five charts

`comoving` leaves $a(t, r)$ and $b(t)$ free, so no component assumes a field equation; it is checked to be the published comoving chart of `kantowski_sachs` slot for slot, with $a$ free in $r$.

`ruban` is Ruban's chart in the cycloid's parameter.
$a$ is a name the chart defines, which the reader holds as a function of $\eta$ and $r$ because its definition holds $\mu(r)$.
`DattRubanCycloid` writes every value as a rational function of $t = \cot(\eta/2)$, $q = \epsilon - \mu\eta$ and $\mu$, among which nothing is related, so a value that vanishes is exactly zero; it then writes each curvature component in $a$ itself by $q = (a - 2\mu)/t$, which leaves no $\cos(\eta/2)$ below the line, and the connection in $q$.
The checker reads the functions of $\eta/2$ alone, so $\mu\sin\eta$ is printed $2\mu\sin(\eta/2)\cos(\eta/2)$.
The chart is checked to be dust of density $2\mu/(a\,b^2)$ and, at $\mu = 1/2$, to be the published dust chart of `kantowski_sachs` along $\eta = 2\eta_K + \pi$ with $\epsilon = \pi/2 - \kappa$.

`areal` takes the radius of the spheres for its time, $T = r_s\sin^2(\eta/2)$ on $0 < \eta < \pi$: Datt's own form as Krasiński writes it, Krasiński and Giono's (6.2), and Plebański and Krasiński's (19.101) by their footnote 5.
There $\partial_Ta = (2\mu T - r_sa)/(2T(r_s - T))$ is rational in $a$, so every value is printed in $a$, $\mu$, $T$ and $r_s$ by `datt_ruban_rates`, with no arcsine anywhere.
It is checked against the derivative of its own $a$, to be dust, to be Ruban's chart carried along its map, and at $\mu = 0$, $\epsilon = 1$ to be the published inside of Schwarzschild's horizon.

`de_sitter` is Ruban's (19): $b = \ell\cosh(ct/\ell)$ and $a = \epsilon\sinh(ct/\ell) + \mu\left(\left(\tfrac{\pi}{2} - \arctan\sinh(ct/\ell)\right)\sinh(ct/\ell) - 1\right)$, with $\partial_{ct}a = (a\cosh^2 + \mu)/(\ell\sinh\cosh)$.
It is held to $G^\mu{}_\nu + 3/\ell^2 = \mathrm{diag}(-2\mu/(a\,b^2), 0, 0, 0)$.
$a(0) = -\mu$ is negative, so the dust lives where $a > 0$, after the curve $a = 0$ for $\epsilon = 1$, and the rate's $\sinh$ below the line never vanishes there.

`exterior_kruskal` is the vacuum outside a T-sphere, the white hole's chart of the same name on the side $V \ge U$ of the surface, with the areal radius held as Lambert's function.

## Step 4: the junction

In the areal chart the block of $T$, $\theta$ and $\phi$ holds no $a$, $\mu$ or $\epsilon$ and does not change along $r$.
So a surface of constant $r$ has the induced metric $-dT^2/(r_s/T - 1) + T^2d\Omega^2$ of the inside of Schwarzschild's horizon and no extrinsic curvature, $K_{ab} = \partial_rg_{ab}/2a = 0$, whatever the dust is, and the same holds from the vacuum's side for a line of constant Schwarzschild $t$ inside the horizon.
A T-sphere of any rest mass therefore joins Schwarzschild's vacuum of the same $r_s$ with no shell on the surface, which is Ruban's constant active mass.
In Kruskal's coordinates the surface is $V = U = -\cos(\eta/2)\,e^{\sin^2(\eta/2)/2}$, from $-1$ on the past singularity through the crossing of the horizons at $\eta = \pi$ to $1$ on the future one.

## Step 5: the drawings

Spacetime diagrams, one for each chart.
The comoving, Ruban and areal planes draw one tube that runs on in both directions, $\epsilon = 1$ and $\mu = (2 + \tanh(r/r_s))/2\pi$, above the threshold on every shell; the comoving chart's free functions are declared from Ruban's solution through `dr_eta`, the cycloid's parameter at the dust's proper time, and checked to solve the published $G^r{}_r = 0$ and $G^\theta{}_\theta = 0$.
The de Sitter plane draws the same $\mu$ at $\ell = 1$, hatched below the curve $a = 0$, which is marked as the singularity it is.
Kruskal's plane draws the vacuum outside a T-sphere with the surface on its left edge; a row's `curves` may now say that the drawn coordinate is signed, since Kruskal's $V$ is negative on half of the surface.

The conformal diagram and the embedding draw a T-sphere, $\epsilon = 1$ and $\mu = 1/\pi$ on $r \le 0$, the member symmetric in time.
`conformal.TSphere` sends the dust to the left half of the Kantowski-Sachs lens by $p, q = \arctan(\sigma \mp r)$ with $d\sigma = b\,d\eta/a$, $\sigma_m = (\pi/2)\,1.2189 = 1.9147\,r_s$, and fits Kruskal's manifold to the surface by $P(U)$ and $Q(V)$, each singularity outside a level line; the event horizon is then the ray $p = 0$, which runs back through the dust to the bang at $r = -\sigma_m$.
The embedding is a movie of the dust's proper time from the greatest expansion to $c\tau = 1.5\,r_s$: the tube of radius $b$, of which the stretch $-2r_s \le r \le 0$ is drawn, $a|r|$ long, and outside it Novikov's slice from its throat, which is the surface of the dust, Flamm's paraboloid at $\tau = 0$.
The four planes of the tube carry no moment of the embedding, which draws another member, and `slices.HIDDEN` says so.
