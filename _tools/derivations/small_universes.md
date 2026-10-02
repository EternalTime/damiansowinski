# Small universes

A Friedmann universe whose space is a quotient of flat or of hyperbolic space by a discrete group of isometries that moves every point: a torus, another closed flat form, a closed hyperbolic manifold, or a horn.
This note records the sources, the four charts, what each is checked against, and what the drawings draw.

## Step 1: the sources

G. F. R. Ellis, "Topology and cosmology", General Relativity and Gravitation 2, 7-21 (1971), doi:10.1007/BF02450512, checked on Crossref and INSPIRE.
Springer's free preview shows pages 7 and 8, which were read: the paper was "Read on 18 May 1970 at the Gwatt Seminar on the Bearings of Topology upon General Relativity", and the quotation in the history is from its second paragraph.
The rest is behind the paywall and was not read, so what the paper covers is cited to Luminet and Roukema, arXiv:astro-ph/9901364, section 1: "Ellis detailed the classification of 3-dimensional Riemannian manifolds useful for cosmology and started to explore the observational consequences of a toroidal universe."
The name George Ellis is from the arXiv author list of gr-qc/0302094; the 1971 paper prints G. F. R. Ellis.

The candidate list's entry gave one paper and it verifies.
Three identifiers looked up along the way were wrong in the brief handed to the search and are corrected in `references.bib`: the translation of Schwarzschild's lecture is doi:10.1088/0264-9381/15/9/003, Ellis and Schreiber is doi:10.1016/0375-9601(86)90032-0, and Sokolov and Starobinskii's English translation carries the ADS year 1976 for a Russian original of 1975.

What was read for each statement of the history:

- Levin, arXiv:gr-qc/0108043, section 2: the video game torus, the lattice of images, "globally anisotropic", "can be made infinitesimally small", Mostow's rigidity.
- Schwarzschild's lecture of 9 August 1900: the date and place from Kragh, arXiv:1205.4909; the quotation from the 1998 translation as Luminet and Roukema's abstract quotes it. Neither the original nor the translation was read.
- Friedmann 1924: the closing remark as Luminet quotes it, arXiv:0704.3579. Neither the paper nor its 1999 translation was read.
- Lachièze-Rey and Luminet, arXiv:gr-qc/9605010: Einstein's equations as local, the 18 flat forms with ten closed and six orientable (section 6), ghosts (section 11.1), rigidity (section 8.3), and "small universes" credited to Ellis and Schreiber (section 10.1, reference 40).
- Sokolov and Shvartsman 1974: the abstract and section 2, jetp.ras.ru. Gott 1980: the summary, ADS scan. Ellis and Schreiber 1986: the abstract.
- Zel'dovich and Starobinskii 1984, ADS scan, page 135: "The simplest example is a 3-torus, that is, three-dimensional Euclidean space with the identification $x \equiv x + L_1$, $y \equiv y + L_2$, $z \equiv z + L_3$", of volume $L_1L_2L_3$, which leaves space homogeneous and violates its isotropy.
- Sokolov and Starobinskii, ADS scan, page 630: the horn's metric, their (2), $dl^2 = R^2[dx^2 + e^{-2x}(dy^2 + dz^2)]$ with $R$ the radius of curvature, the toroidal horn's identifications, their (11), the map (5) to the spherical form, and "every spliced hyperbolic universe is globally inhomogeneous".
- Gabai, Meyerhoff and Milley, arXiv:0705.4325, abstract. The Weeks manifold's volume 0.9427, in-radius 0.5192 and out-radius 0.7525 are Luminet and Roukema's table 1.
- Stevens, Scott and Silk 1993; Cornish, Spergel and Starkman 1998; Cornish, Spergel, Starkman and Komatsu 2004; Planck 2015 XVIII; Petersen and others 2023: abstracts.

Initials stand where no record gives a first name: G. Schreiber, D. D. Sokolov, V. F. Shvartsman, A. A. Starobinskii as his papers of 1975 and 1984 print it, and G. D. Mostow.

## Step 2: the charts

The local geometry is Friedmann's, so every chart is his line element with the scale factor left free, and the identifications stand in the domains.

`torus` is the flat universe in comoving Cartesian coordinates with $x$, $y$ and $z$ periodic in $L_1$, $L_2$ and $L_3$, Zel'dovich and Starobinskii's identification; Gott writes the same for the Einstein-de Sitter universe.
`torus_conformal` is the same in the conformal time, $ds^2 = a^2(-d\eta^2 + dx^2 + dy^2 + dz^2)$, Lachièze-Rey and Luminet's (46) and (49) at zero curvature.
`hyperbolic` is the open universe about one observer, $ds^2 = -c^2dt^2 + a^2(d\chi^2 + \sinh^2\chi\,d\Omega^2)$, their (37), whose scale factor is the radius of curvature of space, a length, as their $R(t)$ is; two points related by an element of the group are one point.
`horn` is Sokolov and Starobinskii's (2) with their (11), written with the periods $b_2$ and $b_3$ since $a$ is the scale factor here.

`small_universes_check` in `print_charts.py` holds every chart to a vanishing Weyl tensor.
The two torus charts are held to the published charts of `frw` at $k = 0$ pulled back along $x = r\sin\theta\cos\phi$, $y = r\sin\theta\sin\phi$, $z = r\cos\theta$, and to having no component that depends on $x$, $y$ or $z$, so that every translation is an isometry.
The hyperbolic chart is held to the published comoving chart of `frw` at $k = -1/\ell^2$ carried along $r = \ell\sinh\chi$, with its scale factor $\ell$ times Friedmann's, for every length $\ell$.
The horn is held to the Einstein tensor of the open universe, to having no component that depends on $y$ or $z$, and to being the spherical form of hyperbolic space carried along Sokolov and Starobinskii's (5).

The hyperbolic charts print the factor $a'^2 - 1$ whole, `small_universes_pretty`: factoring would split it into $a' + 1$ and $a' - 1$.
`chart_printer.hyperbolic` passes over the `Tuple` a derivative of a free function holds, which it met here for the first time.

The tag `spherically symmetric` is overruled in `metric_tags.py`: the hyperbolic chart is spherical about one observer, and the identifications leave only a discrete group of rotations.

## Step 3: the scale factors drawn

The torus is drawn cubic, $L_1 = L_2 = L_3 = L$, and filled with dust, $a = (3ct/2L)^{2/3}$, which is one when the Hubble radius $3ct/2$ is $L$.
Its conformal time is $\eta = 2L\sqrt{a}$, the comoving distance light has crossed since the bang, so light has been once round at $ct = L/12$, twice at $2L/3$ and three times at $9L/4$.
An image at the comoving distance $d$, seen from $\eta = 3L$, shows the galaxy at $\eta = 3L - d$, at the time $ct = (3L - d)^3/12L^2$ and the redshift $1 + z = 9L^2/(3L - d)^2$.

The hyperbolic charts are drawn with open dust, $a = A(\cosh\eta - 1)$ and $ct = A(\sinh\eta - \eta)$, which `null_rays.DustSolver` solves from the chart's own $G^\chi{}_\chi = 0$ or $G^x{}_x = 0$ starting from $a = a_0$ and $a' = 6/5$.
There $a'^2 = 1 + 2A/a$ gives $A = 11a_0/50$, $\cosh\eta = 61/11$ and $\sinh\eta = 60/11$, so $\eta = \ln 11$ and $ct = (6/5 - (11/50)\ln 11)a_0 = 0.6725\,a_0$, and the density is $1 - 1/a'^2 = 11/36$ of the critical density.
The rays are checked against $\eta(t) \pm \chi$ with $\eta(t)$ found from the parametric solution by Newton's method, which the tracing does not use.

On the horn two images a period of $y$ apart at $x$ are the distance $d$ apart with $\sinh(d/2a) = b_2e^{-x}/2$, the chord of a horocycle.
With $b_2 = 2\pi$ that is $a\ln 11$, the distance light has crossed by the dashed line, at $x = \ln(\pi\sqrt{11}/5) = 0.734$.

## Step 4: the drawings

Spacetime diagrams: one cell of the torus in each of its charts, the plane of the time and $x$ with its edges one line, and the same plane unrolled over six cells with one galaxy at each of its places and our past light cone from $\eta = 3L$; the hyperbolic chart about one observer with the in-radius and the out-radius of the Weeks manifold's cell; and the horn along its length.
The plane of the time and $y$ of the horn is not drawn: $\Gamma^x{}_{yy} = e^{-2x}$ does not vanish, so its null curves are not geodesics.

Conformal diagram: the cell of the torus in Minkowski's triangle by $p, q = \arctan((\eta \mp x)/L)$, as Gowdy's torus is drawn, in the comoving and in the conformal time, with a ray from the bang on its first three laps; and the horn on the upper half of Minkowski's diamond.

Embedding diagram: the slice $z = 0$ of the torus is a flat torus, which no smooth surface of flat space carries whole, so the circle of $x$ is closed and the circle of $y$ is cut open, a cylinder of circumference $aL$ and length $aL$ whose ends are one circle, played as a movie through $ct = L/96$, $L/12$ and $9L/32$, where light has been half way, once and one and a half times round.
The slice $z = 0$ of the horn at one moment is Beltrami's pseudosphere from its rim $x = 0$ up, at $b_2 = 2\pi$.
The closed hyperbolic forms have no embedding drawn: the equator of the open universe is the hyperbolic plane, whose embedding `frw` records as impossible about a point.

`_tools/test_small_universes.py` holds the published files and the drawings to all of the above.
