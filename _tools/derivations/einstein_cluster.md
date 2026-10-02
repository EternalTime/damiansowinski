# Einstein's cluster

Einstein's cluster of 1939 is a static ball of particles on circular geodesics about its centre, the planes of the orbits spread evenly over every direction, with no pressure along the radius.
Its five charts are written by `_tools/derivations/print_charts.py --metric einstein_cluster`, and `verify_metrics.py --system einstein_cluster/<chart>` checks each in seconds.
Equation numbers of Einstein's are those of Annals of Mathematics 40, 922.

## Step 1. The charts and their sources

- `areal`: the areal radius with the mass function free.
  C. Gilbert rewrote the cluster in these coordinates in 1954, with $\nu' = (e^\lambda - 1)/r$ and the particles' speed $V^2 = (e^\lambda - 1)/2$, his (8) and (10).
  With $e^{-\lambda} = 1 - 2m/r$ those are $\partial_r\Phi = m/(r(r - 2m))$ and $V^2 = m/(r - 2m)$, the form of Florides (1974), Böhmer and Harko (2007), Böhmer and Lobo's (9), and Cardoso and his colleagues' (3) to (5).
  Lake's (8) is the same relation solved for the mass, $m = r^2\partial_r\Phi/(1 + 2r\partial_r\Phi)$.
- `constant_speed`: Einstein's own example, his section 7, in the areal radius.
  His $\sigma = \sigma_0$ constant makes $V$ constant by his (23), $V^2 = 2\sigma/(1 - \sigma)^2$, and the density fall as $1/r^2$.
  Then $m = V^2r/(1 + 2V^2)$, $g_{rr} = 1 + 2V^2$ and $e^{2\Phi} \propto r^{2V^2}$, normalised so that $t$ is Schwarzschild's time outside $r = R$.
- `isotropic`: the same cluster in Einstein's chart, his (2a), $ds^2 = -a(dr^2 + r^2d\Omega^2) + b\,dt^2$, with his (17) and (18b): $a = (1 + \sigma)^4(\rho/\rho_0)^{-4\sigma/(1 + \sigma)}$ and $b = ((1 - \sigma)/(1 + \sigma))^2(\rho/\rho_0)^{4\sigma/(1 - \sigma^2)}$.
  The isotropic radius is written $\rho$ here, since $r$ is the areal radius in the other charts, and his $r_0$ is $\rho_0$.
- `uniform`: Florides's cluster of uniform density, $m = r_sr^3/(2R^3)$, whose redshift function integrates to $e^{2\Phi} = (1 - r_s/R)^{3/2}(1 - r_sr^2/R^3)^{-1/2}$, Böhmer and Lobo's (10).
- `hyperspherical`: the same cluster with $r = a\sin\chi$, $a = \sqrt{R^3/r_s}$, Böhmer and Lobo's (13) with the time kept as Schwarzschild's: $ds^2 = -(\cos^3\chi_0/\cos\chi)\,c^2dt^2 + a^2(d\chi^2 + \sin^2\chi\,d\Omega^2)$.

## Step 2. The mass function is a name

A mass function left free cannot be published beside its redshift function, since $\Phi$ is then an integral of $m$.
So the areal chart declares $\Phi(r)$ as its function and $m$ as a name for $r^2\partial_r\Phi/(1 + 2r\partial_r\Phi)$, which is $G^r{}_r = 0$ solved for $m$.
The reader holds a defined name whose definition holds a declared function as a function of the coordinates, as it holds Szekeres's $E$, so the tensors are built with $m(r)$ and $\Phi(r)$ side by side.
The chart's `reduce` then writes $\partial_r\Phi$ and $\partial_r^2\Phi$ in $m$ and $\partial_r m$, and every value is printed in $m$, $\partial_r m$ and $e^{2\Phi}$.
`verify_metrics.py` writes the name out on both sides of each comparison, so what it checks is the line element with one free function.

## Step 3. What the curvature says

With $T^\mu{}_\nu = \operatorname{diag}(-\varepsilon, 0, p_\perp, p_\perp)$,

$$G^t{}_t = -\frac{2\partial_r m}{r^2},\qquad G^r{}_r = 0,\qquad G^\theta{}_\theta = G^\phi{}_\phi = \frac{m\,\partial_r m}{r^2(r - 2m)}.$$

So $p_\perp = \varepsilon\,m/(2(r - 2m)) = \tfrac{1}{2}\varepsilon V^2$, Gilbert's (5) and (6): the stress across the radius is the momentum the particles carry sideways.
$V^2 = m/(r - 2m)$ reaches 1 at $r = 3m$, Einstein's bound, his (6a) and (19b), which reads $\sigma < 2 - \sqrt{3}$ in his isotropic chart.
A circular orbit is stable where $r > 6m(r)$, Zapolsky's result, which is $V^2 < 1/4$ and $\sigma < 5 - 2\sqrt{6}$.

`einstein_cluster_check` holds every chart to these: the areal radius is $\sqrt{g_{\theta\theta}}$, the mass function is the Misner-Sharp mass, $G^r{}_r$ and every off diagonal component vanish, $G^t{}_t = -2m'/r^2$, and $2(r - 2m)G^\theta{}_\theta = -m\,G^t{}_t$.
Each closed form has Schwarzschild's lapse at its surface, and the isotropic and hyperspherical charts are the constant speed and uniform charts pulled back.

## Step 4. How the values are printed

The constant speed chart's values are each a rational function of $V$ times a power of $r/R$, printed by `named_powers`.
Einstein's chart names his $a$ and $b$ as parameters, and `einstein_powers` writes every value as a rational function of $\sigma$ times whole powers of $\rho$, $a$ and $b$: the canonical form does not take a power whose exponent is a rational function of $\sigma$ through a product, so the pullback of that chart is compared at six points to forty digits.
The uniform chart's values each hold at most one power of $\sqrt{R - r_s}$ and of $\sqrt{R^3 - r_sr^2}$, and `einstein_radicals` writes them so.
The checker's canonical form holds the second radical as $i\sqrt{r_sr^2 - R^3}$, a generator with the right square, which is why the pullback of the hyperspherical chart reads the lapse as the line element writes it.
The hyperspherical chart's denominators are written in $\cos\chi$.

## Step 5. What is drawn

Every drawing is in units of the cluster's Schwarzschild radius.
Einstein's cluster is drawn at $V = 1/2$, where $m = r/6$ at every radius, the edge of stability, and $R = 3\,r_s$; in his chart that is $\sigma = 5 - 2\sqrt{6}$ and $\rho_0 = R/(1 + \sigma)^2 = 2.47\,r_s$.
Florides's is drawn with the same mass and radius, $a = \sqrt{27}\,r_s$ and $\sin\chi_0 = 1/\sqrt{3}$.
The areal chart is drawn for a declared cluster with no surface, $e^{2\Phi} = 1 - r_s/\sqrt{r^2 + b^2}$ at $b = 2\,r_s$: its density is positive everywhere, its mass approaches $r_s/2$, and its fastest particles move at $0.37\,c$.
The slice of Einstein's cluster is a cone, since $g_{rr}$ is constant, and the slice of Florides's is the cap of Schwarzschild's star of the same density.
