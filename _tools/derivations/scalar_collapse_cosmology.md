# The Husain-Martinez-Núñez scalar field universe

Husain, Martinez and Núñez's solution is Fisher's static solution for a mass with a massless scalar field, at the exponent $\gamma = \pm\sqrt{3}/2$, multiplied by the conformal time,

$$ds^2 = \frac{\eta}{\eta_0}\left[-f^{\gamma}d\eta^2 + \frac{dr^2}{f^{\gamma}} + f^{1-\gamma}r^2d\Omega^2\right], \qquad f = 1 - \frac{b}{r}.$$

Its four charts are written by `_tools/derivations/print_charts.py --metric scalar_collapse_cosmology`, from `_tools/derivations/scalar_collapse_cosmology.py`, and `verify_metrics.py --metric scalar_collapse_cosmology` checks all four in seconds.

## Step 1: the parameters and the units

Husain, Martinez and Núñez, Phys. Rev. D 50, 3783, equations (3) to (5), write the metric with $G = c = 1$ as $(at + b)$ times Fisher's bracket, with $f^2 = (1 - 2c/r)^\alpha$, $\alpha = \pm\sqrt{3}/2$, and the field $\phi = \pm(4\sqrt{\pi})^{-1}\ln[d\,(1 - 2c/r)^{\alpha/\sqrt{3}}(at + b)^{\sqrt{3}}]$ for $R_{\mu\nu} = 8\pi\phi_\mu\phi_\nu$.
The charts shift the time so that the factor is $\eta/\eta_0$, with $\eta = t + b/a$ and $\eta_0 = 1/a$, write their $2c$ as $b$, as the chart of `fisher_jnw` does, and their $\alpha$ as $\gamma$.
The conformal time and $\eta_0$ are lengths, as on `sultana_dyer`, so the conformal chart has no $c$.
The expanding branch is $\eta > 0$ with $\eta_0 > 0$ and its time reverse $\eta < 0$ with $\eta_0 < 0$; in both the factor $\eta/\eta_0$ is positive.
The field is printed in the scaling of `roberts` and `fisher_jnw`, $R_{\mu\nu} = (8\pi G/c^4)\partial_\mu\varphi\,\partial_\nu\varphi$, so $\varphi = \sqrt{c^4/4\pi G}\,\phi$ for their $\phi$ with $d = 1$.

## Step 2: the sources of the charts

- `conformal_time`: equations (3) and (4) of Husain, Martinez and Núñez, shifted as above; Husain, Phys. Rev. D 69, 084002, equation (10), and Kastor and Traschen, Class. Quantum Grav. 34, 035012, their (4), write it in the same form.
- `comoving_time`: Faraoni, Galaxies 1, 114, the solution in comoving time, equation (3.162) in the arXiv text (arXiv:1309.4915, labelled `HMNcomoving`), with $c\,dt = a\,d\eta$, $a = (\eta/\eta_0)^{1/2} = (t/t_0)^{1/3}$ and $t_0 = 2\eta_0/3c$.
- `fonarev_conformal`: Fonarev's family as Maeda, arXiv:0704.2731, writes it in equations (6) to (11), $a(\eta)^2$ times the bracket with any exponent, and as Kastor and Traschen write it in their (12); Maeda's $2/(\lambda^2 - 2)$ and Kastor and Traschen's exponent are the $s$ of the chart, $s = (1 - \gamma^2)/(2\gamma^2 - 1)$ on the solutions, and $s = 1/2$ at $\gamma^2 = 3/4$ is the first chart.
- `fonarev`: Fonarev, Class. Quantum Grav. 12, 1739, equation (2.18) in the arXiv text (arXiv:gr-qc/9409020), signature $(+,-,-,-)$, $e^{8\alpha^2aT}(1 - 2m/r)^{\delta}dT^2 - e^{2aT}[\dots]$, with Fonarev's $a$ written $k/c$, $2m$ written $b$, $\delta$ written $\gamma$ and $4\alpha^2$ written $n$; equation (2.7) gives $\delta = 2\alpha/\sqrt{4\alpha^2 + 1}$, so $n = \gamma^2/(1 - \gamma^2)$, and $V_o = a^2(3 - 4\alpha^2)/8\pi G$, which vanishes at $n = 3$, $\gamma^2 = 3/4$.

Three more forms are not charts of the page.
The advanced time $v$ of Wang, Wu and Gao, Eur. Phys. J. C 79, 823, needs $\int f^{-\gamma}dr$, which has no closed form in elementary functions, so no component of it can be printed.
Faraoni's description in the comoving time and the areal radius $R = a\,rf^{(1-\gamma)/2}$ needs $r$ as a function of $t$ and $R$, which is implicit.
Fonarev's solution for a scalar field coupled conformally to the curvature, equations (3.29) and (3.36), solves other field equations, and equation (3.36) writes $dr^2 + r^2d\theta^2$ where the isotropic radius is meant while equation (3.41) writes $2m/\bar r$ where (3.36) has $m/2\bar r$; the History tells it and no chart is printed.

## Step 3: what is checked

`chart_check` in `scalar_collapse_cosmology.py` holds, at three points in forty digits:

- each chart to $R_{\mu\nu} = 2\,\partial_\mu\varphi\,\partial_\nu\varphi + Wg_{\mu\nu}$, in the scaling of step 1, and $\varphi$ to its wave equation $\Box\varphi = W'(\varphi)/2$, for both signs of the massless exponent and, in the two charts of Fonarev's family, for $\gamma = 4/5$ and $\gamma = -3/5$ as well, with $W = 0$ where $\gamma^2 = 3/4$;
- each chart off the solutions to failing the first, so the check can fail;
- the conformal chart to Faraoni's Ricci scalar, $R = \gamma^2b^2f^{\gamma-2}\eta_0/(6r^4\eta) - 3\eta_0/(2\eta^3f^{\gamma})$, and to a null gradient of the areal radius on the apparent horizon of Husain, Martinez and Núñez's (11), $\eta = r^2f^{1-\gamma}/(2r - b(1 + \gamma))$;
- each chart after the first to being the chart before it pulled back: the comoving chart onto the conformal one along $\eta = \eta_0(t/t_0)^{2/3}$, the first chart of the family onto the conformal one at $s = 1/2$, and Fonarev's chart onto the family's along $\eta = \eta_0e^{kT/sc}$ with $s = 1/(n - 1)$ and $\eta_0 = sc/k$.

`_tools/test_scalar_collapse_cosmology.py` holds the conventions' printed fields and potential to the published Ricci tensor, and the drawings to the numbers their captions state.

## Step 4: what the diagrams draw

Every drawing takes $b = 1$, the unit of every length, and $c = 1$; `scalar_collapse_cosmology_drawings.py` holds them all.
The massless members are drawn at $\eta_0 = 2b$, so the comoving chart has $t_0 = 4b/3c$, and Fonarev's family at $\gamma = 4/5$, where $s = 9/7$ and $n = 16/9$, with $k = c/b$, so that $T = 0$ is $\eta = \eta_0 = 9b/7$.
On the plane of the time and $r$ every chart is a factor times $f^{\gamma}(-d\eta^2 + dr_*^2)$, with the tortoise radius $r_* = b f^{1-\gamma}\,{}_2F_1(2, 1-\gamma; 2-\gamma; f)/(1-\gamma)$, the integral of $f^{-\gamma}dr$ from $r = b$, finite there for $|\gamma| < 1$.

- Spacetime diagrams: one plane for each chart, both massless members in the conformal chart. The rays keep $\eta \pm r_*$, which `null_rays.py --verify` checks, held beyond $r = 1.01\,b$, where the rays stand nearly upright, and in the comoving chart beyond $ct = 0.1\,b$, where $\eta$ is steep in $t$. The apparent horizon is the marker `areal=True` finds, where $|\nabla R|^2$ vanishes for $R^2 = g_{\theta\theta}$; at $\gamma = \sqrt{3}/2$ it rises from the big bang to $\eta = 3.58\,b$ within $0.0125\,b$ of $r = b$, a quarter of a cell of the grid it is found on, so the marker starts beside $r = b$ at that time, falls to $1.67\,b$ at $r = 1.72\,b$ and grows as $r \approx 2\eta$. Each plane that reaches $\eta = 0$ declares the time positive in `where`: beneath the big bang the chart runs on into the time reverse, where its signature turns over and the two families trade orientation, so a ray seeded beside $\eta = 0$ ends on it where the declaration makes the chart's null directions undefined. At $\gamma = 4/5$ the Kretschmann scalar passes $10^8$ beside $r = b$ too slowly for the test at $10^{-5}$ of the drawing, so Fonarev's two planes declare $r = b$ in `singular_zero`.
- Conformal diagram: $p, q = \arctan((\eta \mp r_*)/\ell)$ with $\ell = 4b$ carries each chart onto the triangle $T > 0$ of Minkowski's compactification: the big bang on $T = 0$ from the corner to $i^0$, the timelike singularity $r = b$ on $X = 0$ up to $i^+$, and future null infinity between them, the picture of Husain's Fig. 1 of 2004 for the expanding branch. One view for each chart, both massless members in the conformal chart.
- Embedding diagram: the equatorial plane of the moment $\eta = \eta_0$, where $a = 1$, for each massless member, from $\rho = a\,rf^{(1-\gamma)/2}$ and $g_{rr} - (d\rho/dr)^2 = a^2f^{-1-\gamma}(\gamma/r - (1+\gamma)^2/4r^2)$. At $\gamma = \sqrt{3}/2$ this vanishes on $r_e = (1+\gamma)^2b/4\gamma = 1.00518\,b$; beyond it the plane is drawn in flat space, inside it in three dimensional Minkowski space, as Fisher's plane is drawn on `fisher_jnw`. At $\gamma = -\sqrt{3}/2$ it is negative everywhere, and the whole plane is drawn in Minkowski space. Toward $r = b$ both surfaces close on the light cone, a chord's interval being about $15(r - b)^{1/2}$ of its extent in $\rho$ at $\gamma = \sqrt{3}/2$ and $(r - b)^{1/2}$ at $\gamma = -\sqrt{3}/2$, so fourteen decimals hold the isometry to the checks' part in 5000 only down to $r - b = 10^{-7}b$, a circle of radius $0.34\,b$, and $6 \times 10^{-4}b$, a circle of radius $10^{-3}b$; there the profiles stop and their edges say the surface runs on into the axis. The moment is the same shape at every $\eta$, magnified by $a$, so one moment is drawn.
- Curvature figure: $K$ and $R$ of the moment $\eta = \eta_0$ against $r$, one line view for each massless member, from the conformal chart with $\eta$ held. Far out both tend to the flat Friedmann universe's $K = 15\eta_0^2/4\eta^6$ and $R = -3\eta_0/2\eta^3$, which `curvature.py` checks at $b = 0$ and at $r = 10^8b$; $R$ changes sign at $r = 1.014\,b$ and $1.72\,b$, and at $\gamma = -\sqrt{3}/2$ $K$ is least, $0.153/b^4$, at $r = 4.01\,b$.

The moments lie on the drawings of their own member alone: Fonarev's four drawings are of another member of the family, and each massless member's drawings hold only its own moment, which `exceptions/scalar_collapse_cosmology.json` and `HIDDEN` in `slices.py` record.

## Step 5: the names in the History

Each name is as its article's author line prints it, in the Crossref or INSPIRE record each entry of `assets/data/references.bib` is stamped with, and the people the live pages already name are written as those pages write them.
