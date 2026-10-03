# The Big Rip

The flat Friedmann universe filled with phantom energy alone, a perfect fluid with $p = w\rho c^2$ and $w < -1$.
`big_rip` in `print_charts.py` prints its three charts and `big_rip_check` holds each to its source before anything is written.

## Step 1: the scale factor

Friedmann's equations for flat space with one fluid of constant $w$ give $\rho \propto a^{-3(1+w)}$, which grows with $a$ when $w < -1$.
The solution that reaches infinity at $t = 0$ is
$$a(t) = \left(-\frac{t}{t_0}\right)^{2/3(1+w)}, \qquad t < 0,$$
which is case E of Chiba, Takahashi and Sugiyama, arXiv:astro-ph/0501661 section 2, $a = (-t)^{2/3(1+w)}$, with time measured in the unit $t_0$.
At $t = -t_0$ the scale factor is $1$, and $t_0 = 2/(3|1+w|H_0)$ with $H_0$ the Hubble rate there, which is Caldwell, Kamionkowski and Weinberg's time left, $t_{\rm rip} - t_0 \simeq \tfrac{2}{3}|1+w|^{-1}H_0^{-1}(1-\Omega_m)^{-1/2}$ (arXiv:astro-ph/0302506), at $\Omega_m = 0$.
The Hubble rate is $H = 2/(3(1+w)t)$, positive for $t < 0$.

`big_rip_check` holds every chart to a perfect fluid at rest, $G^i{}_i = -w\,G^0{}_0$ on each spatial axis and no flux, and the comoving charts to Friedmann's equation $G^0{}_0 = -3H^2/c^2$.

## Step 2: the charts

- `cartesian`: $ds^2 = -c^2dt^2 + a^2(dx^2 + dy^2 + dz^2)$, the spatially flat Friedmann universe of Caldwell's phantom cosmology (arXiv:astro-ph/9908168) and of Nojiri, Odintsov and Tsujikawa (arXiv:hep-th/0501025, section II). Checked against the comoving chart pulled back along $x = r\sin\theta\cos\phi$, $y = r\sin\theta\sin\phi$, $z = r\cos\theta$.
- `comoving`: $ds^2 = -c^2dt^2 + a^2(dr^2 + r^2d\Omega^2)$, the left side of Chiba, Takahashi and Sugiyama's (1).
- `conformal`: $ds^2 = a^2(-d\eta^2 + dr^2 + r^2d\Omega^2)$, the right side of their (1), with $d\eta = c\,dt/a$. Integrating, $\eta = -\eta_0(-t/t_0)^{(1+3w)/3(1+w)}$ with $\eta_0 = 3(1+w)ct_0/(1+3w)$, and then $a = (-\eta/\eta_0)^{2/(1+3w)}$, their $\eta \propto -(-t)^{(1+3w)/3(1+w)}$. Checked against the comoving chart pulled back along $\eta(t)$.

Each chart names the scale factor $a$ as a parameter defined by its expression, as Weyl's chart of the Curzon-Chazy particle names $R$.
`big_rip_scale` prints every value as a rational function of $w$ times a whole power of the time and a whole power of $a$: the checker's canonical form writes $(-t/t_0)^q$ as $(-1)^q t^q t_0^{-q}$, and the printer gathers those factors back.

## Step 3: curvature

In the comoving charts
$$R = \frac{4(1 - 3w)}{3c^2t^2(1+w)^2}, \qquad K = \frac{16(5 + 6w + 9w^2)}{27c^4t^4(1+w)^4},$$
both diverging as $t \to 0$, and the Weyl tensor vanishes, since every Friedmann universe is conformally flat.
The Ricci scalar is $8\pi G(\rho - 3p/c^2)/c^2 = 3H^2(1 - 3w)/c^2$, which is the first expression with $H = 2/(3(1+w)t)$.

## Step 4: the event horizon

Light from the centre covers the comoving distance $\int_t^0 c\,dt'/a = -\eta(t)$ before the rip, so the event horizon of the observer at $r = 0$ is $r = -\eta$, with proper radius
$$a(-\eta) = \frac{3(1+w)}{1+3w}\,c(-t),$$
Chiba, Takahashi and Sugiyama's $R_c$ for $w < -1$, and Caldwell, Kamionkowski and Weinberg's $3\delta t(1+w)/(1+3w)$ with $c = 1$.
At $w = -3/2$, the value every drawing takes, $a = (-t/t_0)^{-4/3}$, $\eta_0 = 3ct_0/7$, the horizon's comoving radius is $\eta_0(-t/t_0)^{7/3}$ and its proper radius $3c(-t)/7$.
The Hubble sphere $R = c/H = 3c(-t)/4$ lies outside it.

## Step 5: the drawings

The spacetime diagrams (`null_rays.py`) draw the planes of $t$ and $x$, of $t$ and $r$ and of $\eta$ and $r$ at $w = -3/2$, each marking the event horizon as the ray through $r = \eta_0$ at $t = -t_0$ (or through $r = \eta_0$ at $\eta = -\eta_0$).
The conformal diagram (`conformal.py`) is the lower half of Minkowski space's, $p, q = \arctan((\eta \mp r)/\eta_0)$.
The embedding diagram (`embedding.py`) draws the flat plane $y = 0$ with a ring of galaxies at rest on $x^2 + z^2 = (ct_0/2)^2$, of proper radius $act_0/2$, and the event horizon, of proper radius $3c(-t)/7$; the ring leaves the horizon at $t = -(7/6)^{3/7}t_0$.
