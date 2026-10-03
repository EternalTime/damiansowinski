# Unruh's acoustic black hole

Unruh's acoustic metric of a spherically symmetric, steady, convergent flow, for the flow Visser calls the canonical acoustic black hole: a fluid of constant density, so that continuity makes the radial velocity $v = -c\,r_0^2/r^2$.
With $c$ the speed of sound, the laboratory chart is $ds^2 = -c^2dt^2 + (dr + c\,r_0^2\,dt/r^2)^2 + r^2d\Omega^2$, up to the constant factor $\rho/c$, which is dropped.
`print_charts.unruh_acoustic_hole` writes the two charts, `unruh_acoustic_hole_check` holds each to its source, and `_tools/test_unruh_acoustic_hole.py` holds the published files and drawings to the same facts.

## Step 1. The sources, each verified before anything was built

- W. G. Unruh, "Experimental Black-Hole Evaporation?", Phys. Rev. Lett. 46, 1351 (1981), doi:10.1103/PhysRevLett.46.1351. Crossref and the publisher's abstract page; the published pages read in full from a teaching copy, with the exponents of the last paragraph read off the scanned page.
- V. Moncrief, "Stability of stationary, spherical accretion onto a Schwarzschild black hole", Astrophys. J. 235, 1038 (1980), doi:10.1086/157707. Crossref; read in full from the NASA ADS scan, whose first page prints "Vincent Moncrief".
- Matt Visser, "Acoustic black holes: horizons, ergospheres and Hawking radiation", Class. Quantum Grav. 15, 1767 (1998), doi:10.1088/0264-9381/15/6/024, arXiv:gr-qc/9712010. Crossref; read in full from the arXiv source.
- Carlos Barceló, Stefano Liberati and Matt Visser, "Analogue Gravity", Living Rev. Relativ. 14, 3 (2011), doi:10.12942/lrr-2011-3, arXiv:gr-qc/0505065. Crossref; the history chapter read in version 3 of the arXiv source, which is the review of 2011.
- Theodore Jacobson, "Black-hole evaporation and ultrashort distances", Phys. Rev. D 44, 1731 (1991), doi:10.1103/PhysRevD.44.1731. Crossref and the publisher's abstract page.
- Jeff Steinhauer, "Observation of quantum Hawking radiation and its entanglement in an analogue black hole", Nature Phys. 12, 959 (2016), doi:10.1038/nphys3863, arXiv:1510.00621. Crossref; the published abstract from the journal's page and the speeds from the arXiv text.
- Juan Ramón Muñoz de Nova, Katrine Golubkov, Victor I. Kolobov and Jeff Steinhauer, "Observation of thermal Hawking radiation and its temperature in an analogue black hole", Nature 569, 688 (2019), doi:10.1038/s41586-019-1241-0, arXiv:1809.00913. Crossref; the published abstract from the journal's page.

The candidate list also names Weinfurtner, Tedford, Penrice, Unruh and Lawrence (2011), verified at Crossref; the draining bathtub tells that experiment, and this History leaves it to that page.
Unruh's metric in his time prints $c\,dr^2/(c^2 - v^2)$, which balances only as $c^2dr^2/(c^2 - v^2)$, Visser's (57); the History writes the balanced form and cites both.

## Step 2. The charts and where each comes from

| chart | source | line element |
| --- | --- | --- |
| `laboratory` | Visser's (55), with the lower sign of a flow falling inward | $-c^2dt^2 + (dr + c\,r_0^2\,dt/r^2)^2 + r^2d\Omega^2$ |
| `unruh` | Unruh's time $\tau = t + \int v\,dr/(c^2 - v^2)$, Visser's (56) and (57) | $-(1 - r_0^4/r^4)\,c^2d\tau^2 + dr^2/(1 - r_0^4/r^4) + r^2d\Omega^2$ |

With $v = -c\,r_0^2/r^2$ the map is $c\,d\tau = c\,dt - r_0^2r^2\,dr/(r^4 - r_0^4)$, and integrated, $c(\tau - t) = -\tfrac14 r_0\ln|(r - r_0)/(r + r_0)| - \tfrac12 r_0\arctan(r/r_0)$, with the constant chosen to vanish at $r = 0$.

## Step 3. What the check holds

- The laboratory chart against Visser's (55) expanded, slot by slot.
- A ray with $d\vec x/dt = \vec v + c\,\hat n$ is null for every unit $\hat n$.
- $g^{rr} = 1 - r_0^4/r^4$, zero on the horizon $r = r_0$, where $g_{tt}$ vanishes too, so the horizon and the ergosurface coincide, as Moncrief notes of every spherical flow; each moment of constant $t$ is flat space.
- Unruh's chart is the laboratory chart pulled back along the map above.
- $R = 6r_0^4/r^6$ and $K = 468\,r_0^8/r^{12}$ in both charts, singular at the sink alone.
- Unruh's $T = (\hbar/2\pi k_B)\,\partial v/\partial r$ on the horizon is $\hbar c/(\pi k_Br_0)$, and Visser's surface gravity $\tfrac12\partial(c^2 - v^2)/\partial r$ is $2c^2/r_0$.

## Step 4. The drawings

All are at $r_0 = c = 1$, with the horizon at 1.

- Spacetime diagrams: the radial plane on the equator in each chart. The closed forms are $t + r - \arctan r$ and $t - r - \tfrac12\ln|(r - 1)/(r + 1)|$ in the laboratory, and $\tau \pm r_*$ with $r_* = r + \tfrac14\ln|(r - 1)/(r + 1)| - \tfrac12\arctan r$ in Unruh's time.
- Conformal diagram: $f = 1 - 1/r^4$ has surface gravity 2 and $UV = (1 - r)e^{4r - 2\arctan r}/(1 + r)$, so $r = 0$ is the line $T = \pi/2$, as in Kruskal's diagram. Two roots of $f$ are imaginary, so `conformal.Tower` does not apply, and the map is written out. The laboratory's $v = t + r - \arctan r$ is regular, so the laboratory chart is the ingoing one and covers the outside and the black hole; the fluid is all there is, and the edge is $t \to -\infty$.
- Embedding diagram: the laboratory's $t = 0$ is flat, so its equator is a plane, and Unruh's $\tau = 0$ climbs at $dz/dr = 1/\sqrt{r^4 - 1}$, so $z = \mathrm{F}(\arccos(1/r) \mid 1/2)/\sqrt 2$, which levels off at $\mathrm{K}(1/2)/\sqrt 2 \approx 1.311$. Both moments are marked on every spacetime and conformal diagram.
