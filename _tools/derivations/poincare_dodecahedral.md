# The Poincaré dodecahedral universe

A closed Friedmann universe whose space is the 3-sphere divided by the binary icosahedral group $I^*$, so that space is one regular dodecahedron with each face glued to the face opposite it after a tenth of a turn.
This note records the sources, the three charts, what each is checked against, and what the drawings draw.

## Step 1: the sources

Every paper the candidate list names was verified on Crossref and on arXiv, authors, title, volume, pages and DOI, and read in its arXiv version:

- Luminet, Weeks, Riazuelo, Lehoucq and Uzan, *Nature* 425, 593-595 (2003), doi:10.1038/nature01944, arXiv:astro-ph/0310253, whose first page says the arXiv text is "A slightly edited version" of the published one. The quiet quadrupole ("about 1/7"), the vanishing of correlations beyond 60°, the bell, the gluing and the "illusion", the edge angles of 120° and about 117°, the tiling by 120 and the volume 120 times smaller, the ranges of $\Omega_0$ for the octopole and the quadrupole, "no free parameters in its construction", global homogeneity, the inradius 0.31, outradius 0.39 and horizon radius 0.38 at $\Omega_0 \simeq 1.013$, and the six pairs of circles of about 35°.
- Cornish, Spergel, Starkman and Komatsu, *Physical Review Letters* 92, 201302 (2004), arXiv:astro-ph/0310233: no matched circles of radius above 25°, and "the search excludes the Poincaré Dodecohedron suggested in [7] as this model predicts back-to-back circles of radius 35°".
- Roukema, Lew, Cechowska, Marecki and Bajtlik, *Astronomy & Astrophysics* 423, 821-831 (2004), arXiv:astro-ph/0402608: $\Omega_{\rm tot} = 1.009$ brings the circles' radius nearly to zero, the strongest match at $11 \pm 1°$ for the left handed screw motion of 36° and not the right handed one or none, and the screw motion of $\pi/5$, "equal to the in-diameter of the fundamental domain". The names are the arXiv author line; Crossref prints initials.
- Aurich, Lustig and Steiner, *Classical and Quantum Gravity* 22, 2061-2083 (2005), arXiv:astro-ph/0412569: the line element (6) with (8), the coordinates (7) and (29), the quaternions and the Clifford translations (18), the group $I^*$ of order 120 and its two generators, the wave numbers (26), the 59 eigenfunctions of Luminet's group and their own 10521 up to $\beta = 155$, $\Omega_{\rm tot}$ from 1.016 to 1.020, circles of 40° to 50°, and "the signal of the six pairs of matched circles could be missed".
- Key, Cornish, Spergel and Starkman, *Physical Review D* 75, 084034 (2007), arXiv:astro-ph/0604616: the "soccer ball universe", the peak at 11° reproduced and gone once the mean of each circle is removed, the reach to circles of about 5°, and "effectively rules out the Poincaré Dodecahedral model as an interesting shape for the Universe".

Two more were read for what came before 2003:

- Gott, *Monthly Notices of the Royal Astronomical Society* 193, 153-169 (1980), read on the ADS scan, section 4: "The group $I^*$ divides $S^3$ into 120 regular dodecahedral cells, the nearest images are at $\chi = \pi/5$", the cell inscribed in a sphere of radius $0.388\,a_0$, no space form with a smaller largest radius of a cell, and, from the summary, "the proper radius of the universe cannot be much smaller than its radius of curvature". The wish for "cells that are small in all dimensions relative to the radius of curvature of the sphere" is section 4 too. Lachièze-Rey and Luminet quote the radius as $0.338$, a slip: Gott's scan and the computation below give $0.388$.
- Lachièze-Rey and Luminet, *Physics Reports* 254, 135-214 (1995), arXiv:gr-qc/9605010, section 7.4.4: the Poincaré manifold as $S^3/I$, its faces "identified after rotating by 1/10th turn".

Poincaré's paper of 1904 and Weber and Seifert's "Die beiden Dodekaederräume" of 1933 were verified on Crossref but are behind paywalls and were not read, so neither is cited.

## Step 2: the charts

The local geometry is Friedmann's closed universe, so every chart is his line element with the scale factor left free, and the identification stands in the domains and in the shared convention.
The scale factor is the radius of curvature of space, a length, as Aurich, Lustig and Steiner's $R(t)$ is.

- `comoving` is their (6) with (8), $ds^2 = -c^2dt^2 + a^2(d\chi^2 + \sin^2\chi\,d\Omega^2)$, in the angles (7) about one observer.
- `conformal` is the same in their conformal time, $a\,d\eta = c\,dt$, a pure number, $ds^2 = a^2(-d\eta^2 + d\chi^2 + \sin^2\chi\,d\Omega^2)$.
- `toroidal` is their (29), $dv^2/(2v(1 - 2v)) + (1 - 2v)d\alpha^2 + 2v\,d\gamma^2$ on the unit 3-sphere at $(\sqrt{1 - 2v}\cos\alpha, \sqrt{2v}\cos\gamma, \sqrt{2v}\sin\gamma, \sqrt{1 - 2v}\sin\alpha)$.

`poincare_dodecahedral_check` in `print_charts.py` holds each chart to a vanishing Weyl tensor, to a time orthogonal to space with Friedmann's lapse, and to its space being $a^2$ times the unit 3-sphere of $\mathbb{R}^4$ pulled back through (7) or (29).
In the toroidal chart it holds the shift $(\alpha, \gamma) \to (\alpha + s, \gamma + s)$ to the left multiplication by $\cos s + \sin s\,ij$, their (18) with $a_k = \cos s$ and $d_k = \sin s$.
Every element of $I^*$ of real part $\cos(\pi/5)$ is conjugate in $S^3$ to $\cos(\pi/5) + \sin(\pi/5)\,ij$, and conjugation is a rotation, so with the group turned to hold that element the toroidal chart's identification includes the shift by $\pi/5$; `test_poincare_dodecahedral.py` checks both.
The charts in cosmic time print $a'^2 + 1$ whole, `poincare_dodecahedral_pretty`.
The tag `spherically symmetric` is overruled in `metric_tags.py`, since the group leaves only a finite group of rotations about any point.

## Step 3: the group, computed

`dodecahedral.py` builds $I^*$ as the closure of Aurich, Lustig and Steiner's generators $j$ and $\sigma/2 + i/2\sigma + j/2$, and `test_poincare_dodecahedral.py` holds it to what the texts state:

- 120 elements, none but 1 with a fixed point, each moving every point the same distance;
- twelve nearest images at exactly $\pi/5$, so the inradius is $\pi/10$;
- twenty corners of the Dirichlet cell about 1, each the point equidistant from 1 and three nearest images that no image is nearer to, all at $0.38814$, Gott's $0.388$;
- the multiplicity of each wave number $\beta$, $\beta$ times the invariants of $I^*$ in the representation of $SU(2)$ of dimension $\beta$, nonzero exactly at Aurich, Lustig and Steiner's list (26), never at an even $\beta$, and 59 in all at 13, 21 and 25.

## Step 4: the drawings

The comoving and toroidal charts are drawn with closed dust at rest at its largest radius $a_m$, $a = a_m(1 - \cos\psi)/2$ and $ct = a_m(\psi - \sin\psi)/2$, solved by `null_rays.DustSolver` from the chart's own $G^\chi{}_\chi = 0$ or $G^\alpha{}_\alpha = 0$.
Both draw the great circle through us and the centre of a face, with the group turned so that the face lies at $\theta = \pi/2$, $\phi = 0$: a Clifford translation through $\pi/5$ carries it onto itself, so it passes through ten copies of the cell.
In the toroidal chart that circle is $v = 0$, and its identification is the chart's own shift of $\alpha$.

The conformal chart is drawn with Luminet's universe, $\Omega_m = 0.28$ and $\Omega_\Lambda = \Omega_0 - \Omega_m$ with $\Omega_0 = 1.013$, their Figure 4, in units of the radius of curvature $a_0$ today, with $H_0a_0/c = 1/\sqrt{\Omega_0 - 1}$ and $\Lambda = 3\Omega_\Lambda H_0^2/c^2 = 169.15/a_0^2$; radiation is left out.
`DustSolver` solves the published $G^\chi{}_\chi = -\Lambda$, and `conformal.py` checks the same scale factor against $a' = H_0\sqrt{\Omega_m a + \Omega_\Lambda a^4 - (\Omega_0 - 1)a^2}$.
Today is $\eta_0 = 0.38882$, the last scattering at $a = a_0/1101$ is $\eta = 0.012988$, and future infinity is $\eta = 0.5168$, so the sphere of last scattering has the radius $0.3758$, Luminet's "about 0.38", and light sent today reaches $\chi = 0.128$.
Cornish, Spergel and Starkman's $d = 2R_c\arctan(\tan(R_{\rm lss}/R_c)\cos\alpha)$ with $d = \pi a_0/5$ gives matched circles of $\alpha = 34.6°$, Luminet's "about 35°".
The sphere of last scattering lies inside our cell and its twelve neighbours, and reaches into all twelve, which a sample of 400000 of its points confirms.

The embedding diagram is the equator $\theta = \pi/2$ of space today in Luminet's universe, a sphere of radius $a_0$, drawn as two hemispheres so that each is a graph over its circles.
On it the faces between adjacent cells $g$ and $h$ are the great circles $n\cdot(g - h) = 0$, kept where $g$ and $h$ are both nearer than every other image; the sphere crosses 46 of the 120 cells, and the section of our own cell reaches from $\chi = \pi/10$ to $0.3879$.
Each curve's height is read off the drawn profile, which near the equator stands up to $10^{-4}\,a_0$ from the round sphere.
The moment today is marked on the conformal chart's spacetime diagram and on the conformal diagram, and hidden on the closed dust, another member of the family.
