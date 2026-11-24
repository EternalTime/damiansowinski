# The wormholes of a conformally coupled scalar field

Carlos Barceló and Matt Visser, "Traversable wormholes from massless conformally coupled scalar fields", Physics Letters B 466, 127-134 (1999), arXiv gr-qc/9908029, read in full from its text.
Sourya Ray, "Revisiting the static spherically symmetric solutions of gravity with a conformally coupled scalar field", Classical and Quantum Gravity 42, 017002 (2025), arXiv 2406.14345, section 4 at $D = 4$, for the parameters $\mu$ and $b$.

## Step 1. The field equations

Barceló and Visser write Einstein's equations as $\kappa G_{\mu\nu} = T_{\mu\nu}$ with $\kappa = c^4/8\pi G$ and the stress tensor of their (2.1), which has no trace, so $R = 0$.
Their (2.2) and (2.3) are $R_{\mu\nu} = (\kappa - \phi_c^2/6)^{-1}\left(\tfrac23\partial_\mu\phi_c\partial_\nu\phi_c - \tfrac16 g_{\mu\nu}(\nabla\phi_c)^2 - \tfrac13\phi_c\nabla_\mu\nabla_\nu\phi_c\right)$ and $\Box\phi_c = 0$.
Their solutions are the Janis-Newman-Winicour-Wyman metrics (2.6) times $\Omega^2$, with $\Omega$ their (2.17), $\Omega = \alpha_+f^{\sigma} + \alpha_-f^{-\sigma}$, $f = 1 - 2\eta/r$ and $\sigma = \sin\chi/2\sqrt3$, and the field their (2.14), $\phi_c = \pm\sqrt{6\kappa}\,(\alpha_+f^{\sigma} - \alpha_-f^{-\sigma})/(\alpha_+f^{\sigma} + \alpha_-f^{-\sigma})$.
Their angle $\Delta$ has $\tan(\Delta/2) = (\alpha_+ - \alpha_-)/(\alpha_+ + \alpha_-)$; with $\alpha_+ + \alpha_- = 1$ the page writes $b = \tan(\Delta/2)$, so $\alpha_\pm = (1 \pm b)/2$, and $\mu = \eta/2$.
The page keeps $c$: each time term is $c^2dt^2$, $t$ a time, and every component is in the chart $x^0 = ct$.

## Step 2. The charts

`isotropic` is their (3.22) in the isotropic radius $\bar r$ of their (3.21), $r = \bar r(1 + \eta/2\bar r)^2$, written $\rho$ here, since the checker's reader takes no barred coordinate; Fisher-Janis-Newman-Winicour's isotropic chart writes the same radius $\rho$.
With $\alpha_\pm = (1 \pm b)/2$ the bracket of (3.22) is $(\rho - b\mu)/(\rho + \mu)$, so $ds^2 = (1 - b\mu/\rho)^2\left(-(1 + \mu/\rho)^{-2}c^2dt^2 + (1 + \mu/\rho)^2(d\rho^2 + \rho^2d\Omega^2)\right)$, Ray's form, and the field (3.23) is $\phi_c = \sqrt{6\kappa}\,(b\rho - \mu)/(\rho - b\mu)$.
For $b < 0$ the throat is $\rho = \mu\sqrt{-b}$, their $\bar r = (\eta/2)|\tan(\Delta/2)|^{1/2}$, of areal radius $\mu(1 + \sqrt{-b})^2$, and $\rho \to 0$ is a second flat region.
`three_parameter` is the whole class in their $r$, (2.6) times $\Omega^2$, with $f = 1 - 4\mu/r$, the line element $\Omega^2\left(-F\,c^2dt^2 + dr^2/F + (f/F)r^2d\Omega^2\right)$ with $F = f^{\gamma}$ and $\gamma = \cos\chi$, and $\Psi = ((1 + b)f^{\sigma} - (1 - b)f^{-\sigma})/2$, so that $\phi_c = \sqrt{6\kappa}\,\Psi/\Omega$.
$F$, $\Omega$ and $\Psi$ are held as names with declared slopes, $dF/dr = 4\mu\gamma F/r(r - 4\mu)$, $d\Omega/dr = 4\mu\sigma\Psi/r(r - 4\mu)$ and $d\Psi/dr = 4\mu\sigma\Omega/r(r - 4\mu)$, so every value is a rational function of them, $r$ and the parameters.
$\gamma$ and $\sigma$ lie on the conic $\gamma^2 + 12\sigma^2 = 1$, which $R = 0$ needs, and the chart defines them through its rational parametrisation by $w = \tan(\chi/2)/2\sqrt3$: $\gamma = (1 - 12w^2)/(1 + 12w^2)$ and $\sigma = 2w/(1 + 12w^2)$, covering $\chi \in [0, \pi)$.
A definition by a root fails: the checker's canonical form writes the root of a radicand with a negative leading term as $i$ times the root of its negative, which evaluates to minus the root, and $\sqrt{1 - 12\sigma^2}$ came back as $-\gamma$, the member at $\pi - \chi$; a defined $\sqrt{(1 - \gamma^2)/12}$ split into $\sqrt{\gamma - 1}\sqrt{\gamma + 1}$ with the wrong sign the same way.
`three_parameter_pretty` writes each value back in $\gamma$ and $\sigma$ by $w = \sigma/(1 + \gamma)$, reduced by $12\sigma^2 = 1 - \gamma^2$ to the first power of $\sigma$.
Barceló and Visser note that their areal radius (2.18) cannot be inverted in closed form, so no chart is given in it.

## Step 3. The checks

`check_chart` holds the isotropic chart, in forty digits at three points of $(\mu, b, \rho)$, to (2.2) in every slot, to (2.3) and to $R = 0$; the areal radius to its least value $\mu(1 + \sqrt{-b})^2$ at $\rho = \mu\sqrt{-b}$; far out to $g_{tt} = -1 + 2\mu(1 + b)/\rho$, so $G_{\rm eff}M = \mu(1 + b)$, their (3.26); at $b = 0$ to the published isotropic chart of `majumdar_papapetrou` with $m = \mu$; and on $\rho > \mu$ to $\Omega^2$ times the published isotropic chart of `fisher_jnw` at $\gamma = 1/2$ and Fisher's $b = 4\mu$.
`check_three` holds the chart in $r$, at four points of $(\mu, b, w, r)$, one a wormhole and one with $\gamma < 0$, to (2.2), (2.3) and $R = 0$; to $\Omega^2$ times the published spherical chart of `fisher_jnw` with Fisher's $b = 4\mu$; far out to $r_s = 4\mu(\gamma + 2b\sigma)$, twice their $G_{\rm eff}M = \eta(\cos\chi + \tan(\Delta/2)\sin\chi/\sqrt3)$ of footnote 5; and at $w = 1/6$, where $\chi = \pi/3$, to the isotropic chart along $r = \rho(1 + \mu/\rho)^2$ on $\rho > \mu$.
The first version of `check_three`, with $\sigma$ defined by its root, failed (2.2) at its first point, which is how the wrong branch was found.

## Step 4. The drawings

Every drawing of the wormholes is of one, $\mu = 1$ and $b = -1/4$, with its throat at $\rho = \mu/2$, inside the sphere $\rho = \mu$ where the chart in $r$ ends, so that the drawings show the continuation through that sphere.
On the plane of $t$ and $\rho$ the metric is $(1 - b\mu/\rho)^2(1 + \mu/\rho)^{-2}\left(-c^2dt^2 + (1 + \mu/\rho)^4d\rho^2\right)$, so the rays keep $ct \mp \rho_*$ with $\rho_* = \rho + 2\mu\ln(\rho/\mu) - \mu^2/\rho$ for every $b$, the closed form `null_rays.py --verify` checks.
The spacetime diagram of the isotropic chart runs from $\rho = 0$ to $3\mu$ and marks the throat and the sphere $\rho = \mu$; $g^{\rho\rho}$ vanishes at $\rho \to 0$ with the diverging conformal factor, the second flat region, so the row declares `where_is_infinity` and no horizon is marked there.
The chart in $r$ is drawn at another member, $w = 1/4$ and $b = -1/4$, $\chi = 81.8°$, $\gamma = 1/7$ and $\sigma = 2/7$, whose spheres shrink to a point at the naked singularity $r = 4\mu$, as the areal radius goes as $f^{1/7}$ there.
At $w = 1/10$, a member whose spheres grow toward $r = 4\mu$, the throat stands within $0.1\mu$ of the singularity, too close to read, and at $w = 1/6$ the chart in $r$ draws only the side $\rho > \mu$ of the wormhole the isotropic chart draws whole.
The moment $t = 0$ of the wormhole does not lie on that view, which `exceptions/conformal_scalar_wormholes.json` and `HIDDEN` say.
The conformal diagram is the full diamond, $p, q = \arctan((ct \mp (\rho_* - \rho_{*t}))/\mu)$, with the throat on its axis, two flat regions and no horizon, the same for every $b < 0$; the chart in $r$ has no view of its own, since its member drawn would be Minkowski's triangle with the singularity on its left edge, and the page's conformal diagram is the wormhole's.
The embedding diagram is the slice $t = 0$ on the equator, one piece through the throat from $\rho = \mu/20$ to $5\mu$, two spheres each the other's image under $\rho \to \mu^2(-b)/\rho$, which keeps the areal radius $\rho + \mu(1 - b) - b\mu^2/\rho$ and the proper distance, so the surface is the same on both sides of the throat for every $b < 0$; the generator checks equal circumferences and heights.
The curvature figure plots the stored Kretschmann scalar of the isotropic chart against the proper distance from the throat, $l = \rho - \rho_t + \mu(1 - b)\ln(\rho/\rho_t) + b\mu^2(1/\rho - 1/\rho_t)$, a closed horizontal expression, from $-7.5\mu$ to $7.5\mu$, greatest within $\mu/10$ of the throat, and that of the chart in $r$ at the naked member against $r$ from $4.1\mu$ to $8\mu$, growing without bound toward $r = 4\mu$; both stored Ricci scalars vanish, so each view has $K$ alone.
