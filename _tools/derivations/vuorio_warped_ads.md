# Vuorio's universe and warped anti-de Sitter space: the charts, the field equations and what the diagrams draw

The family has two parameters, the twist $\Omega$ of the Killing vector $\partial_t$ and the curvature scale $m$ of the plane it turns about, both inverse lengths.
This note records where each chart of `vuorio_warped_ads.json` comes from, what it is checked against, and what the diagrams take from it.
Units here have $c = 1$ unless a line keeps it; the file keeps $c$.

## Step 1. The papers, and how Vuorio's own was read

Every paper the History cites was verified on Crossref and INSPIRE before anything was built, with its author line and DOI.
Vuorio's own paper (Phys. Lett. B 163, 91, 1985, MIT-CTP-1267) is on Elsevier's site alone, which refuses this machine, and KEK holds no scan of it.
Its metric was read from two accounts that quote it: Chow, Pope and Sezgin's (0906.3559) quotation of his equation (2.21),
$$ds^2 = \frac{9}{\mu^2}\left[-\left(dt + 2\,d\theta - 2\cosh\sigma\,d\theta\right)^2 + d\sigma^2 + \sinh^2\sigma\,d\theta^2\right],$$
and Ait Moussa, Clement and Leygnac's (gr-qc/0303042) equation (2), the same solution at $\mu = 3$ with the integration constant Vuorio set to $-2$ for regularity on the axis.
Percacci, Sodano and Vuorio's paper of 1987 was read whole from KEK's scan of its preprint, KEKSCAN 2000-36-627, `https://lib-extopc.kek.jp/preprints/PDF/2000/0036/0036627.pdf`; it calls Vuorio's solution "the well-known Gödel-like solution" and shows that the untwisted stationary vacua are flat.
Deser, Jackiw and Templeton's and Nutku's papers are cited for their abstracts alone, read on INSPIRE.

## Step 2. The family and its field equations

Every chart has the form $-(c\,dt + A)^2 + h$, with $h$ the metric of a hyperbolic plane of Gaussian curvature $-m^2$ and $dA = 2\Omega$ times its area.
These are Rebouças and Tiomno's homogeneous metrics of Gödel's type with $m^2 > 0$ (Rebouças and Santos, 0906.5354, their $H(r)$ and $D(r)$, signature reversed).
Gürses (0812.2576, Proposition 13) shows that such a metric solves topologically massive gravity, $G_{\mu\nu} + \Lambda g_{\mu\nu} + C_{\mu\nu}/\mu = 0$, with $\mu = 3\Omega$ and a background curvature fixed by the cosmological constant; his $\lambda$ is $-\Lambda$, and his $r_2 = -2(w^2 + 3\lambda)$ with $r_2 = -2m^2$ gives $\Lambda = (\Omega^2 - m^2)/3$.
`vuorio_check` in `print_charts.py` holds every chart to that equation exactly, with $C_{\mu\nu} = \epsilon_\mu{}^{\alpha\beta}\nabla_\alpha(R_{\beta\nu} - \tfrac14 R g_{\beta\nu})$ and $\epsilon^{012} = +1/\sqrt{-g}$ in the chart's own order of coordinates, the convention of Anninos, Li, Padi, Song and Strominger; the sign of $\mu$ follows the orientation, and every chart is written so that it comes out $+3\Omega$.
The trace of the equation gives $R = 6\Lambda = 2(\Omega^2 - m^2)$, which every chart publishes, and the Kretschmann scalar is $44\Omega^4 - 24\Omega^2m^2 + 4m^4$, anti-de Sitter space's $12/\ell^4$ at $m = 2\Omega$ with $\ell = 1/\Omega$.
Vuorio's universe is $m = \Omega$, $\Lambda = 0$; Gödel's slice is $m^2 = 2\Omega^2$ (Rebouças and Santos); anti-de Sitter space is $m = 2\Omega$, which is Rooman and Spindel's $\mu = 1$ and Boyda, Ganguli, Horava and Varadarajan's $m^2 = 8\Omega^2$ in their normalisation with $4\sqrt{2}\Omega/m^2$.
In Rooman and Spindel's notation the stretching is $2\Omega/m$, and in Anninos et al.'s timelike warping $4\nu^2/(\nu^2 + 3) = 4\Omega^2/m^2$ with $\ell^2/(\nu^2 + 3) = 1/m^2$.

## Step 3. The four charts

cylindrical: Rebouças and Tiomno's $H = (4\Omega/m^2)\sinh^2(mr/2)$, $D = \sinh(mr)/m$, with $r$ the proper distance from the axis. `vuorio_against_vuorio` checks that at $m = \Omega$ it is Vuorio's (2.21) pulled back through $t_V = \Omega ct$, $\sigma = \Omega r$, $\theta = -\phi$, exactly. The circles of constant $t$ and $r$ are null at $\tanh(mr/2) = m/2\Omega$, Rebouças and Santos's critical radius, which is $r_c = \ln 3/\Omega$ at Vuorio's member.

disc: Bengtsson and Sandin's chart on Poincaré's disc (gr-qc/0509076, section 4), $R = \tanh(mr/2)$, in which the null circle is $R = m/2\Omega$, their $1/\lambda$.

fibred: Anninos et al.'s timelike warped anti-de Sitter space, their (3.4), written in $\Omega$ and $m$: the plane in coordinates about its geodesic $\sigma = 0$, $dA = (2\Omega/m^2)\cosh\sigma\,d\sigma\wedge du$.

horospherical: Rooman and Spindel's (4), Gödel's own coordinates carried to every member, $A = (2\Omega/m)e^{mx}dy$.

Each chart after the first is checked against the cylindrical one in `vuorio_check`: the plane pulled back through the map onto the hyperboloid $X_0^2 - X_1^2 - X_2^2 = 1$, and the difference of the twisting forms closed, so that a shift of $t$ carries one chart onto the other; both at six random points in forty digits.
`slices.py` checks Rooman and Spindel's transformation with its shift, $t_h = t + (2\Omega/m^2)(2\arctan(e^{-mr}\tan(\phi/2)) - \phi)$, $e^{mx} = \cosh mr + \cos\phi\sinh mr$, $y\,e^{mx} = \sinh(mr)\sin(\phi)/m$, and the disc's $R = \tanh(mr/2)$, by pulling the published metrics back.

## Step 4. The diagrams

Every drawing is of Vuorio's member, $m = \Omega = 1$.
Spacetime diagrams, five views: the cylinders of $t$ and $\phi$ at $r_c/2$ and $3r_c/2$, where the null lines are $c\,dt = (-H \pm D)\,d\phi$; the plane of $t$ and $R$ through the disc's centre, with $ct \pm 2\,\mathrm{artanh}\,R$ constant; and the planes $u = 0$ and $y = 0$, each $-c^2dt^2 + d\sigma^2$ or $-c^2dt^2 + dx^2$. At $r_c/2$ the curve moving to $+\phi$ is a null geodesic, since $2\lambda + 6\sinh^2(r/2) - 1 = 0$ there.
The turning figure of light cones about one axis, which in two dimensions of space is the whole spacetime out to $2r_c$.
The embedding diagram of the moment $t = 0$ about one axis: $\rho = 2s\sqrt{1 - 3s^2}$ with $s = \sinh(r/2)$, widest at $s^2 = 1/6$, $\rho = 1/\sqrt{3}$, and stopping where $9s^4 + 6s^2 - 2 = 0$, $s^2 = (\sqrt{3} - 1)/3$, inside the null circle $s^2 = 1/3$.
No conformal diagram: closed timelike curves pass through every event, as in Gödel's universe.
