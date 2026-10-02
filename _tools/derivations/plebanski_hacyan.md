# Plebański-Hacyan and anti-Nariai

The products of two surfaces of constant curvature, $K_1$ for the Lorentzian one and $K_2$ for the other, solve the Einstein-Maxwell equations with $\Lambda = (K_1 + K_2)/2$ and a uniform field of energy density $\varepsilon$, $8\pi G\varepsilon/c^4 = (K_2 - K_1)/2$.
Ortaggio and Podolský's table 1 (Class. Quantum Grav. 19, 5221, gr-qc/0209068) lists the six with $K_2 \ge K_1$: Minkowski, Nariai, anti-Nariai, Bertotti-Robinson, and Plebański and Hacyan's two.
The charts are written by `_tools/derivations/print_charts.py --metric plebanski_hacyan`, which took 3 seconds on 2 October 2026, and `verify_metrics.py --system plebanski_hacyan/<chart>` checks each in under a second.

## Step 1. The source of each chart

Plebański and Hacyan's paper of 1979, J. Math. Phys. 20, 1004, is behind a paywall, so each line element is taken from a paper that quotes it or derives it.

- `sphere`, $K_1 = 0$, $K_2 = 1/b^2$: Cardoso, Dias, and Lemos's (28) at $D = 4$, $-dt^2 + dx^2 + \rho_u^2d\Omega^2$, with their $x$ written $z$ and $\rho_u = b$.
- `sphere_rindler`: the chart their limit lands in, $-\chi^2dT^2 + d\chi^2 + \rho_u^2d\Omega^2$, with their $T$ written $\tau$.
- `plane`, $K_1 = -1/a^2$, $K_2 = 0$: Plebański and Hacyan's own chart as Podolský and Ortaggio write it, their (37) at $A_1 = 0$, $2\,d\zeta\,d\bar\zeta - 2\,du\,dv + 2\Lambda v^2du^2$, with $\zeta = (x + iy)/\sqrt{2}$, $\Lambda = -1/2a^2$, and their $v$ written $w$, the letter Ortaggio and Podolský's (14) has.
- `plane_null`: Ortaggio and Podolský's (11), $2\,du\,dv/(1 + \Lambda uv)^2$, with $v$ reversed so that the signature is the collection's; $w = v/(1 + uv/2a^2)$ carries it onto `plane`.
- `plane_static`: Cardoso, Dias, and Lemos's (37) at $k = 0$ and $D = 4$, with $1/A = a^2$.
- `anti_nariai`, $K_1 = K_2 = -1/a^2$: Dias and Lemos's (17) at $K_0 = 1$, the neutral solution, with $R_0 = a$.
- `anti_nariai_static`: their (21) and (22), with $R$ written $r$.
- `exceptional`: Podolský and Ortaggio's (37), with $A_1\zeta + \bar A_1\bar\zeta$ written $f\,x + g\,y$ for two real functions of $u$.

## Step 2. What the script checks

`plebanski_hacyan_check` holds $G^\mu{}_\nu$ to $-K_2$ on the block of the first two coordinates and $-K_1$ on the other, which with $\Lambda = (K_1 + K_2)/2$ is the stress of a uniform field along the first surface, $G^\mu{}_\nu + \Lambda\delta^\mu{}_\nu = \tfrac{1}{2}(K_2 - K_1)\,\mathrm{diag}(-1, -1, 1, 1)$.
The exceptional chart passes the same check whatever $f$ and $g$ are: they enter the Christoffel symbols and no curvature component.
The Rindler chart is held to being `sphere` pulled back, the null chart to being `plane` pulled back, `anti_nariai` to being its static chart pulled back along $r = a\cosh\chi$ and $ct = a\tau$, and the exceptional chart to being `plane` at $f = g = 0$.

## Step 3. The drawings

Each product is drawn in units of its one radius.
The spacetime diagrams hold the other surface at one point, the equator of the sphere, the origin of the flat plane, or $\theta = 1$ of the hyperbolic plane, off the pole where $g^{\phi\phi}$ diverges.
The exceptional chart has no drawing: $\Gamma^x{}_{uu} = f/2$ carries a ray out of the plane of $u$ and $w$, and a drawing would need a choice of $f$ and $g$.
The conformal diagrams place Plebański and Hacyan's chart in the strip of the anti-de Sitter factor through the Poincaré patch, $ct = u - a^2/w$ and $z = a^2/w$ on $w > 0$, with $q$ lowered by $\pi$ on $w < 0$, and the static charts through $\tan q = e^{\tau + \chi_*}$ and $\tan p = -e^{-(\tau - \chi_*)}$ with $\chi_* = \ln\tanh(\chi/2)$.
The moment of anti-de Sitter space times a flat plane is a flat plane, so it has no embedding diagram, and anti-Nariai's hyperbolic plane is drawn in Minkowski space as the hyperbolic horizon of `topological_black_hole` is.

## Step 4. The limits

`_tools/test_plebanski_hacyan.py` holds the published static chart of `reissner_nordstrom_de_sitter` at $9r_s^2/4 = 8r_q^2 = 2/\Lambda$ to a triple horizon at $r = 1/\sqrt{2\Lambda} = b$ with $r_q^2/b^4 = \Lambda$, and the published hyperbolic chart of `topological_black_hole` at $\mu = -2L/(3\sqrt{3})$ to a double horizon at $r_h = L/\sqrt{3}$ with $f''(r_h)/2 = 1/r_h^2$, which is anti-Nariai at $a = r_h$.
