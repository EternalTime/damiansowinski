# The draining bathtub

Visser's acoustic metric of a fluid that swirls down a drain, $ds^2 = -c^2dt^2 + (dr - A\,dt/r)^2 + (r\,d\theta - B\,dt/r)^2$, with $c$ the speed of sound.
The fluid moves at $(A\,\hat r + B\,\hat\theta)/r$, inward for $A < 0$, and sound moves at $c$ past it, which is all the line element says.
`print_charts.draining_bathtub` writes the three charts, `draining_bathtub_check` holds each to its source, and `_tools/test_draining_bathtub.py` holds the published files and drawings to the same facts.

## Step 1. The sources, each verified before anything was built

- Matt Visser, "Acoustic black holes: horizons, ergospheres and Hawking radiation", Class. Quantum Grav. 15, 1767 (1998), doi:10.1088/0264-9381/15/6/024, arXiv:gr-qc/9712010. Crossref; read in full from the arXiv source.
- W. G. Unruh, "Experimental Black-Hole Evaporation?", Phys. Rev. Lett. 46, 1351 (1981), doi:10.1103/PhysRevLett.46.1351. Crossref, and the abstract from INSPIRE's record, whose source is APS.
- W. G. Unruh, "Sonic analogue of black holes and the effects of high frequencies on black hole evaporation", Phys. Rev. D 51, 2827 (1995), doi:10.1103/PhysRevD.51.2827, arXiv:gr-qc/9409008. Crossref; abstract and text from the arXiv source.
- Carlos Barceló, Stefano Liberati and Matt Visser, "Analogue Gravity", Living Rev. Relativ. 8, 12 (2005), doi:10.12942/lrr-2005-12, arXiv:gr-qc/0505065. Crossref; the history chapter read in version 1 of the arXiv source, which is the review of 2005.
- Ralf Schützhold and William G. Unruh, "Gravity wave analogues of black holes", Phys. Rev. D 66, 044019 (2002), doi:10.1103/PhysRevD.66.044019, arXiv:gr-qc/0205099. Crossref; read from the arXiv source.
- Soumen Basak and Parthasarathi Majumdar, "'Superresonance' from a rotating acoustic black hole", Class. Quantum Grav. 20, 3907 (2003), doi:10.1088/0264-9381/20/18/304, arXiv:gr-qc/0203059. Crossref; read from the arXiv source.
- Emanuele Berti, Vitor Cardoso and José P. S. Lemos, "Quasinormal modes and classical wave propagation in analogue black holes", Phys. Rev. D 70, 124006 (2004), doi:10.1103/PhysRevD.70.124006, arXiv:gr-qc/0408099. Crossref; read from the arXiv source.
- Carlos Barceló, Stefano Liberati, Sebastiano Sonego and Matt Visser, "Causal structure of analogue spacetimes", New J. Phys. 6, 186 (2004), doi:10.1088/1367-2630/6/1/186, arXiv:gr-qc/0408022. Crossref; read from the arXiv source.
- Silke Weinfurtner, Edmund W. Tedford, Matthew C. J. Penrice, William G. Unruh and Gregory A. Lawrence, "Measurement of Stimulated Hawking Emission in an Analogue System", Phys. Rev. Lett. 106, 021302 (2011), doi:10.1103/PhysRevLett.106.021302, arXiv:1008.1911. Crossref; abstract from the arXiv source.
- Theo Torres, Sam Patrick, Antonin Coutant, Maurício Richartz, Edmund W. Tedford and Silke Weinfurtner, "Rotational superradiant scattering in a vortex flow", Nature Phys. 13, 833 (2017), doi:10.1038/nphys4151, arXiv:1612.06180. Crossref; the preprint read from the arXiv source and the published abstract from the journal's page.

The candidate list gave the gain Torres and the others measured as 20 percent.
That is the preprint's abstract; the published abstract says "14% ± 8%", and the History quotes the published figure.
The tank's size, the hole and the depth are from the preprint's text.

## Step 2. The charts and where each comes from

| chart | source | line element |
| --- | --- | --- |
| `laboratory` | Visser's equation for the draining bathtub, section "Vortex geometries" | $-c^2dt^2 + (dr - A\,dt/r)^2 + (r\,d\theta - B\,dt/r)^2$ |
| `kerr_like` | Basak and Majumdar's (6), as Berti, Cardoso and Lemos correct it in their (16) and (17) | $-(1 - (A^2 + B^2)/(c^2r^2))c^2dT^2 + dr^2/(1 - A^2/(c^2r^2)) - 2B\,dT\,d\phi + r^2d\phi^2$ |
| `vortex_filament` | Visser's last equation of the same section | the laboratory chart with $dz^2$ added |

Berti, Cardoso and Lemos take $A > 0$ for a drain, so their $A$ is minus Visser's; the page keeps Visser's sign.
Their $\tilde t$ and $\tilde\phi$ are written $T$ and $\phi$, since the reader takes a tilde on a command and not on a Latin letter.
The map is $dT = dt + A\,r\,dr/(c^2r^2 - A^2)$ and $d\phi = d\theta + A\,B\,dr/(r(c^2r^2 - A^2))$; Basak and Majumdar's second line lacks the $1/r$, which is one of the typos Berti, Cardoso and Lemos correct.
The symbol $c$ is the speed of sound on this page alone, and the shared convention says so; the checker's chart $x^0 = ct$ takes it as it stands.

## Step 3. What the check holds

- The laboratory charts against Visser's expanded form, slot by slot.
- A ray with $d\vec x/dt = \vec v + c\,\hat n$ is null for every unit $\hat n$.
- $g^{rr} = 1 - A^2/(c^2r^2)$, zero on the horizon $r = |A|/c$; $g_{tt} = 0$ on $r = \sqrt{A^2 + B^2}/c$; $g^{tt} = -1$, so $t$ is a time everywhere.
- The Kerr-like chart is the laboratory chart pulled back along the map above.
- $R = 2(A^2 + B^2)/(c^2r^4)$ and $K = 44(A^2 + B^2)^2/(c^4r^8)$ in all three charts: the curvature depends on the speed of the fluid alone.

## Step 4. The drawings

All are of the drain $A = -1$, $B = \sqrt 3$ at $c = 1$: the horizon at 1 and the ergosurface at 2.

- Spacetime diagrams: the plane of time and $r$ with the angle divided out, the rays of zero angular momentum, in each chart, and the spring $A = 1$ in the laboratory chart. The closed forms are $t + r - \ln(1 + r)$ and $t - r - \ln|r - 1|$ for the drain, their time reverse for the spring, and $T \pm r_*$ with $r_* = r + \tfrac12\ln((r - 1)/(r + 1))$.
- A figure in three dimensions, `projections.bathtub`: cones of sound on the horizon, on the ergosurface and at $r = 3.2$, on a floor carrying six streamlines. In the laboratory's angle the cones tip past the vertical inside the ergosurface; they do not all turn one way there, as Kerr's do in Boyer and Lindquist's angle, since $\theta$ and $\phi$ differ by a function of $r$.
- Conformal diagram: $f = 1 - 1/r^2$ has surface gravity 1 and $UV = (1 - r)e^{2r}/(1 + r)$, so $r = 0$ is the line $T = \pi/2$, as in Kruskal's diagram. The laboratory's $v = t + r - \ln(1 + r)$ is regular, so the laboratory chart is the ingoing one and covers the outside and the black hole; the rest of Kruskal's diagram is left out, since the water is all there is, and its edge is $t \to -\infty$.
- Embedding diagram: the laboratory's $t = 0$ is a flat plane, and the Kerr-like chart's $T = 0$, which is $t = \tfrac12\ln(r^2 - 1)$, is the catenoid $z = \operatorname{arcosh} r$. Both moments are marked on every drawing of the drain but the figure, which carries the plane as its floor.
