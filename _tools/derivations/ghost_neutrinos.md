# Davis-Ray ghost neutrinos

Talmadge M. Davis and John R. Ray, Phys. Rev. D 9, 331 (1974): a static plane symmetric solution of the Einstein-Dirac equations whose neutrino field has a current along the symmetry axis and a vanishing energy-momentum tensor, so that the geometry is a vacuum.
Their full text is behind a subscription and was not read; their abstract, and the reference list the publisher deposited with Crossref, which cites Taub's vacuum of 1951, were.
The geometry is taken from M. L. Bedran, M. O. Calvão, F. M. Paiva and I. Damião Soares, Phys. Rev. D 55, 3431 (1997), arXiv:gr-qc/9608058, read in full.

## Charts

All three come from Bedran et al., whose signature is $(+,-,-,-)$ with $c = 1$ and no length scale; each is turned to $(-,+,+,+)$, with $c$ and a length $L$ put back.
$L$ is no physical parameter: rescaling $t$, $x$ and $y$ by constants removes it, as Bedran et al. say Taub's solution has no essential parameter.

- `taub`, their (2.1): $ds^2 = -(z/L)^{-2/3}c^2dt^2 + dz^2 + (z/L)^{4/3}(dx^2 + dy^2)$, with $z > 0$ the proper distance from the singularity.
- `areal`, their (2.2) with Table 1's $r = (3z/2)^{2/3}$: $z = \tfrac{2}{3}L(r/L)^{3/2}$ gives $dz^2 = (r/L)dr^2$, $(z/L)^{-2/3} = (2/3)^{-2/3}L/r$ and $(z/L)^{4/3} = (2/3)^{4/3}(r/L)^2$, and the constants go into $t$, $x$ and $y$: $t_{\rm Taub} = (2/3)^{1/3}t$, $x_{\rm Taub} = (2/3)^{-2/3}x$. Then $ds^2 = -(L/r)c^2dt^2 + (r/L)dr^2 + (r/L)^2(dx^2 + dy^2)$.
- `null`, their (3.4) with Table 1. In the areal chart the tortoise radius is $r_* = r^2/2L$, since $(r/L)dr^2 = (L/r)dr_*^2$, so the $(t, r)$ part is $(L/r)(-c^2dt^2 + dr_*^2)$. With $\kappa = 1/2L$, Table 1's $u = e^{-\kappa r_*}\sinh(\kappa ct)$ and $v = e^{-\kappa r_*}\cosh(\kappa ct)$ give $dv^2 - du^2 = \kappa^2(v^2 - u^2)(dr_*^2 - c^2dt^2)$ and $v^2 - u^2 = e^{-r^2/2L^2}$. Then $\mu = -\ln(v - u)$, $\nu = -\ln(v + u)$ give $dv^2 - du^2 = (v^2 - u^2)\,d\mu\,d\nu$ and $r = L\sqrt{2(\mu + \nu)}$, so $ds^2 = 2\sqrt{2}L^2\,d\mu\,d\nu/\sqrt{\mu + \nu} + 2(\mu + \nu)(dx^2 + dy^2)$, their (3.4) with the sign of the signature turned.

Bedran et al.'s (3.2), the same geometry in $u$ and $v$, carries $\sqrt{\ln(1/(v^2 - u^2))}$; the chart printer reads its two logarithms back in two equal forms that the symbolic comparison cannot match, so $u$ and $v$ are stated in the null chart's convention and used by the conformal diagram, and are not printed as a chart of their own.

`ghost_neutrinos_check` in `print_charts.py` holds every chart to $R_{\mu\nu} = 0$, to $K = 64/27z^4$ (Bedran et al.'s (2.7)) and, but for Taub's own, to Taub's metric pulled back through the maps above, all at points.

## The neutrino field

Derived here and held by `_tools/test_ghost_neutrinos.py`.
With the tetrad $\mathrm{diag}(cA, C, C, 1)$, $A = (z/L)^{-1/3}$, $C = (z/L)^{2/3}$, the spinor $\psi = (AC^2)^{-1/2}\chi = (L/z)^{1/2}\chi$ solves the massless Dirac equation for every constant $\chi$.
Its symmetric energy-momentum tensor has only $tx$ and $ty$ components, each proportional to the $y$ or $x$ component of the current $\bar\psi\gamma^\mu\psi$, so it vanishes exactly when the current points along $z$, that is for $\chi$ of definite chirality and spin along $z$.
The current is null, with $j^z \propto L/z$, and $\sqrt{-g} = cz/L$, so the flux $\sqrt{-g}\,j^z$ is the same on every plane: $z^{-1}\cdot z^{4/3}\cdot z^{-1/3} = z^0$.
