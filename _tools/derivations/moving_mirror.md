# The moving mirror of Fulling and Davies

Flat spacetime of two dimensions to the right of a perfectly reflecting mirror, with a massless quantum field on it.
This note records the six charts, where each comes from, what each is held to, and what each diagram draws.
The sources have $\hbar = c = 1$; the metric file keeps $c$, writes the null coordinates as the lengths $u = ct - x$ and $v = ct + x$, and takes $\kappa$ as an inverse length.

## Step 1: the papers

The two papers the candidate list names were verified on Crossref before anything was built.
S. A. Fulling and P. C. W. Davies, "Radiation from a moving mirror in two dimensional space-time: conformal anomaly", Proceedings of the Royal Society of London A 348, 393-414 (1976), doi:10.1098/rspa.1976.0045.
P. C. W. Davies and S. A. Fulling, "Radiation from moving mirrors and from black holes", Proceedings of the Royal Society of London A 356, 237-257 (1977), doi:10.1098/rspa.1977.0130.
Neither is open, so what each did is taken from its abstract, which Crossref carries in full, and every formula from open papers that attribute it to them.
Those are Good, Anderson and Evans, Physical Review D 88, 025023 (2013), arXiv:1303.6756, and Physical Review D 94, 065010 (2016), arXiv:1605.06635, and Akal, Kusuki, Shiba, Takayanagi and Wei, Physical Review Letters 126, 061604 (2021), arXiv:2011.12005.
The first names of Fulling and Davies are from Crossref's record of Fulling's paper of 1973 and INSPIRE's author record of Davies; C. M. Wilson's is given by no record found and stays as initials.

## Step 2: the inertial and null charts

Good, Anderson and Evans (2016), (2.1) and (2.2): $ds^2 = -dt^2 + dr^2 = -du\,dv$ with $u = t - r$, $v = t + r$, the mirror on $r = z(t)$, and only the part of the spacetime to its right kept.
Their (2.9a) gives the modes $e^{-i\omega v} - e^{-i\omega p(u)}$, with the ray tracing function $p$ defined so that $v = p(u)$ on the mirror, and $f$ its inverse, $u = f(v)$ there, their 2013 paper's (2.13) and (2.15).
The inertial chart names $z(t)$ and the null chart $p(u)$ as declared functions; neither enters the metric, and each chart's domain states the side of the mirror kept.

## Step 3: the chart that brings the mirror to rest

Akal and his coauthors' (1), which they attribute to Fulling and Davies: $\tilde u = p(u)$, $\tilde v = v$ "maps this into a simple setup with a static mirror $\tilde u - \tilde v = 0$".
The chart writes $U$ for their $\tilde u$, since the checker reads a coordinate's name as one symbol.
With $u = f(U)$, $ds^2 = -du\,dv = -f'(U)\,dU\,dv$, a factor times the flat metric of $U$ and $v$, which is the "conformally static nature of the problem" of Fulling and Davies's abstract.
Its one Christoffel symbol is $\Gamma^U{}_{UU} = f''/f'$.

## Step 4: the flux

Fulling and Davies's result, as Good, Anderson and Evans (2016) write it in (4.1): $F(u) = \langle T_{uu}\rangle = \frac{1}{24\pi}\left[\frac{3}{2}\left(\frac{p''}{p'}\right)^2 - \frac{p'''}{p'}\right]$, which with $u$ a length is $\frac{\hbar c^2}{24\pi}[\ldots]$, an energy per unit time.
Davies, Fulling and Unruh (1976) give the same quantity from the conformal factor: for $ds^2 = -C\,dU\,dv$ with no curvature and the vacuum of the modes of $U$ and $v$, $T_{UU} = -\frac{1}{12\pi}C^{1/2}\partial_U^2C^{-1/2}$, and $T_{uu} = T_{UU}/C^2$ since $dU/du = 1/C$.
With $C = f' = 1/p'$ the two agree identically, which `print_charts.moving_mirror_check` holds for a free $f$, writing $p''/p' = -f''/f'^2$ and $p'''/p' = 3f''^2/f'^4 - f'''/f'^3$.
So the flux of every mirror is read from the published metric of its own chart, and no closed form is drawn.

## Step 5: Carlitz and Willey's mirror

Good, Anderson and Evans (2013), Table I: $z = -t - \frac{1}{\kappa}\mathrm{W}(e^{-2\kappa t})$, $p = -\frac{1}{\kappa}e^{-\kappa u}$, $f = -\frac{1}{\kappa}\ln(-\kappa v)$.
The chart `thermal` is Step 3's in $cT = (v + U)/2$, $X = (v - U)/2$: $f' = -1/\kappa U = 1/\kappa(X - cT)$, so $ds^2 = (-c^2dT^2 + dX^2)/\kappa(X - cT)$, with the mirror on $X = 0$ and the chart ending on $X = cT$, which is $u \to \infty$.
The check pulls the inertial chart back along that map, holds the world line to $U = v$ at four times, and holds the flux to $\kappa^2/48\pi$ at every $u$, their (4.3) and the sentence after it.
The mirror turns where $dz/dt = 0$, which is $\mathrm{W} = 1$: $ct = -1/2\kappa$, $z = -1/2\kappa$.

## Step 6: the mirror that imitates a collapse

Good, Anderson and Evans (2016), (2.14) and (2.18) at $v_H = 0$: $z = -t - \mathrm{W}(2e^{-2\kappa t})/2\kappa$ and $u = h(v) = v - \frac{1}{\kappa}\ln(-\kappa v)$.
Their (2.26) is the relation between the retarded times inside and outside a shell of light collapsing to a Schwarzschild black hole, $u_s = u - 4M\ln((v_H - u)/4M)$, the same function at $\kappa = 1/4M$, which with $G$ and $c$ restored is $\kappa = 1/2r_s$.
The chart `collapse` has $f' = 1 - 1/\kappa U = 1 + 1/\kappa(X - cT)$.
Its flux is their (4.2), $\frac{\kappa^2}{48\pi}(4\mathrm{W} + 1)/(\mathrm{W} + 1)^4$ with $\mathrm{W} = \mathrm{W}(e^{-\kappa u})$, which is $-\kappa U$ on the ray, and their (4.5) puts its fastest rise on $\kappa u = \ln 2 - 1/2$, where $\mathrm{W} = 1/2$ and the flux is $16/27$ of its last value.
Far in the past $\mathrm{W}(2e^{-2\kappa t}) \approx -2\kappa t - \ln(-\kappa t)$, so $z \approx \ln(-\kappa ct)/2\kappa$ and the mirror's speed falls as $1/2\kappa t$.

## Step 7: the uniformly accelerating mirror

The hyperbola $x^2 - c^2t^2 = 1/\kappa^2$ has $uv = -1/\kappa^2$, so $p = -1/\kappa^2u$, a ratio of linear functions, for which the Schwarzian form of Step 4 vanishes: "A uniformly accelerating mirror does not radiate", in the words of Fulling and Davies's abstract.
Rindler's coordinates in their conformal form, $ct = e^{\kappa\xi}\sinh(\kappa c\eta)/\kappa$ and $x = e^{\kappa\xi}\cosh(\kappa c\eta)/\kappa$, bring it to rest on $\xi = 0$ with $ds^2 = e^{2\kappa\xi}(-c^2d\eta^2 + d\xi^2)$, a static chart: the "static universe" of the same abstract and the fixed surface of Davies (1975).
The spacetime to the right of this mirror lies inside the wedge $x > c|t|$, so the chart covers all of it.
`metric_tags.py` counts only the inertial, null and mirror at rest charts toward the tags, since this static chart is one mirror's.

## Step 8: the diagrams

Spacetime diagrams, ten views at $\kappa = 1$, every window square.
The inertial chart draws the three mirrors, the null chart and the chart of $U$ and $v$ the two that recede, and each of the last three charts its own mirror.
Each view marks the mirror's world line, hatches what lies behind it, and marks rays seeded a millionth to the right of the mirror and traced both ways, so that each is drawn arriving and leaving: at $\kappa v = -2, -1, -1/2, -1/4, -1/8$ for the receding mirrors and $1/2, 1, 2, 4$ for the hyperbola.
The receding mirrors' views mark the last ray, $v = 0$.
`Diagram.signed_curves` lets a declared world line take either sign, as a mirror's $x$ does.

Conformal diagrams, eight views, with $p = \arctan(\kappa u)$ and $q = \arctan(\kappa v)$: the part of Minkowski's diamond to the right of each mirror.
The receding mirrors run from $i^-$ to the point $(X, T) = (-\pi/2, \pi/2)$ of the left future null infinity, where the ray $v = 0$ ends, and the hyperbola is the vertical line $X = \pi/2$, since $\arctan(-1/u) = \pi/2 + \arctan u$ for $u < 0$.

The embedding diagram's place is taken by a height, as it is for the Krasnikov tube: a moment of this spacetime is a half line, and what makes it what it is is the radiation.
The height is the flux of Step 6 over the plane $0.5 \le \kappa x \le 6.5$, $-3 \le \kappa ct \le 3$, twice its ratio to the thermal flux, read from the published $g_{XX}$ of the `collapse` chart by the formula of Step 4.

`_tools/test_moving_mirror.py` holds the published files to all of this.
