# Nordström's scalar theory of gravity

Gunnar Nordström's second theory of 1913, as A. Einstein and A. D. Fokker wrote it in 1914, is a conformal factor $\Phi^2$ on Minkowski's metric in a preferred chart, with the field equation $R = 24\pi G\,T/c^4$.
This note records the four charts, the source of each, what each is checked against, and what the drawings draw.

## Step 1: the sources

Nordström's second theory is Annalen der Physik 347, 533 (1913), doi:10.1002/andp.19133471303, and Einstein and Fokker's paper is Annalen der Physik 349, 321 (1914), doi:10.1002/andp.19143491009; both records were checked on Crossref.
John Norton's history is Archive for History of Exact Sciences 45, 17 (1992), doi:10.1007/BF00375886, reprinted in The Genesis of General Relativity (2007), doi:10.1007/978-1-4020-4000-9_27.
The quotations of the page's history were checked against the reprint, which is the text Norton hosts.
Nordström's first paper of 1912, Einstein's Vienna lecture of 1913 and Laue's review of 1917 have no record on Crossref or INSPIRE, so the bibliography does not hold them and the history cites Norton for what they say.

## Step 2: the free chart

Einstein and Fokker's first assumption is $ds^2 = \Phi^2(dx^2 + dy^2 + dz^2 - c^2dt^2)$ in preferred coordinates, Norton's equation (56).
The chart `conformal` leaves $\Phi$ a free function of the four coordinates, so every component is an identity.
`print_charts.nordstrom_scalar_check` holds its Weyl tensor to vanishing and its Ricci scalar to $-6\,\Box\Phi/\Phi^3$, Ravndal's equation (15) and Deruelle's (8.7) with $1 + \Phi$ written $\Phi$.
With $g = \Phi^2\eta$ and $\Phi = 1 + U/c^2$ the Newtonian limit is $R = -6\nabla^2U/c^2 = -24\pi G\rho/c^2$, which is $24\pi G\,T/c^4$ for $T = -\rho c^2$, Deruelle's (7.5).
A vacuum is therefore a solution of the wave equation, and superposition holds outside matter.

## Step 3: the point mass

Outside a static spherical body $\Phi = 1 - GM/c^2r$, Deruelle's (9.15) and Deruelle and Sasaki's (3.8).
The chart `spherical` writes $m = GM/c^2$.
It is checked to be the free chart at that $\Phi$, to have $R = 0$ and a Ricci tensor that does not vanish, $R_{\mu\nu}R^{\mu\nu} = 4m^2(m^2 - 4mr + 6r^2)/(r - m)^8$, Deruelle and Sasaki's (3.9), and $K = 2R_{\mu\nu}R^{\mu\nu}$.
In the isotropic radius $-g_{tt} = 1 - 2U + U^2$ and $g_{rr} = 1 - 2U + U^2$ with $U = m/r$, so $\beta = 1/2$ and $\gamma = -1$, as Deruelle states.
The orbit equation is $u'' = m/\ell^2 - (1 + m^2/\ell^2)u$ for $u = 1/r$ and angular momentum $\ell$ per unit mass, so the perihelion turns back by $\pi m^2/\ell^2$ a revolution, a sixth of Einstein's $6\pi m^2/\ell^2$; `_tools/test_nordstrom_scalar.py` runs a planet with the published Christoffel symbols to that number.
The areal radius is $r - m$, which vanishes where $\Phi$ does, and the Kretschmann scalar diverges there as $(r - m)^{-8}$.

## Step 4: the uniform field

Giulini's section 4.1 takes $\Phi$ linear in the height, $1 + gz/c^2$, a vacuum of the theory.
The chart `uniform` counts $z$ from the plane where $\Phi$ vanishes, $\Phi = az/c^2$, so that a body at rest at $z = c^2/a$, where $\Phi = 1$, has the proper acceleration $a$.
Giulini's equations of motion are written in the proper time $s$ of the flat background, and his (26a) is $z = (c^2/a)\sqrt{1 - a^2s^2/c^2}$ in this chart, which reaches the singular plane at $s = c/a$, his (29a).
A clock carried by the body reads $d\tau = \Phi\,ds$, so it reads $\pi c/4a$ on arrival, which is the number the caption states.
The test file runs the fall with the published Christoffel symbols and holds it to both.

## Step 5: the dust universe

Sundrum's section 9 solves the theory for homogeneous dust: $\partial_t^2\Phi$ is a negative constant, his (9.2) to (9.4), so $\Phi$ is a parabola in the inertial time and the universe runs from a bang to a crunch.
The chart `dust` writes it $\Phi = 1 - c^2t^2/L^2$, with $R\Phi^3 = -12/L^2$, so the density is $\rho_0/\Phi^3$ with $2\pi G\rho_0 = c^2/L^2$.
Deruelle and Sasaki's (2.16) carries the opposite sign and no factor $\Phi^3$, and is not followed.
The dust is at rest and in free fall, and a galaxy's clock runs $\int\Phi\,dt = 4L/3c$ between the two singularities.
Every value is printed around the factor $L^2 - c^2t^2$.

## Step 6: the drawings

Every cone in every chart is Minkowski's, since $\Phi^2$ drops out of the null condition; the free chart is drawn for any $\Phi$, with a plane wave $1 + \tfrac{1}{2}\cos(k(x - ct))$ declared.
The conformal diagrams are the parts of Minkowski's diagram where $\Phi > 0$: the whole diamond for the wave, the triangle with a timelike singularity on its edge for the point mass and the uniform field, and the slab between two spacelike singularities for the dust.
The embedding diagram draws the point mass's equator in three dimensional Minkowski space, where $dZ/dr = \sqrt{2mr - m^2}/r$ and $Z = 2s - 2m\arctan(s/m)$ with $s = \sqrt{2mr - m^2}$, and the dust universe's flat moments as a movie.
The uniform field's plane of $x$ and $z$ has no rotation to turn it about, and the wave's has none either, so neither is embedded.
`REGIONS` in `metric_tags.py` counts the free chart alone, since the three solutions have symmetries the theory's general spacetime lacks.
