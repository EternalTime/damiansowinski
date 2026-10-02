# The white hole of Novikov and Ne'eman

The six charts of `white_hole.json` are written by `print_charts.py --metric white_hole`, which computes every tensor from each line element with the checker's own `Geometry` and reads every printed value back before it writes.
This note records the source of each chart, the junction at the surface of the core, and what the diagrams draw.

## Step 1. The spacetime and its sources

The spacetime is Oppenheimer and Snyder's ball of dust with the time reversed: a ball cut from a closed Friedmann universe of dust that leaves a singularity, crosses the Schwarzschild sphere of the vacuum around it outward, comes to rest, and falls back.
I. D. Novikov proposed such cores in 1964 as pieces of the big bang that start late, `novikov1964`, read here in its translation `novikov1965`, all seven pages.
His own exact solution is the marginally bound one: Tolman's dust with no energy function, set in Einstein and Straus's vacuole, with the moment each shell leaves the singularity a function $t_0(r)$ of the shell.
He lists the bound core as his case 2, one that will "expand to some finite dimensions greater than the gravitational radius, after which it will begin to contract again".
Yuval Ne'eman's note of 1965, `neeman1965`, has no scan at ADS and was not opened; what is said of it is taken from Zel'dovich and Novikov, `zeldovichnovikov1967`, and from Ne'eman and Tauber, `neemantauber1967`, who draw the lagging core on Kruskal's diagram.

## Step 2. The charts

- `interior_comoving`: $-c^2d\tau^2 + a^2(d\chi^2 + \sin^2\chi\,d\Omega^2)$ with $a(\tau)$ free, the chart of `oppenheimer_snyder` with $\tau$ counted from the first singularity.
- `interior_conformal`: the same along the cycloid, $a = a_m\sin^2(\eta/2)$, $c\tau = \tfrac{a_m}{2}(\eta - \sin\eta)$, which is dust exactly, $G^\eta{}_\eta = -3/(a_m^2\sin^6(\eta/2))$; the half angle is printed tight.
- `exterior_schwarzschild` and `exterior_eddington_finkelstein`: the published charts of `schwarzschild`, slot by slot. The retarded time $u$ counts from the ray the surface sends as it crosses $r_s$ on its way out.
- `exterior_kruskal`: Kruskal's null coordinates, $ds^2 = -(4r_s^3/r)e^{-r/r_s}\,dU\,dV + r^2d\Omega^2$ with $UV = (1 - r/r_s)e^{r/r_s}$, `kruskal1960`. The areal radius is $r = r_s(1 + \mathrm{W}(-UV/e))$ with Lambert's function, which the checker's `Reader` reads as `\mathrm{W}` and `norm` reduces by $e^{\mathrm{W}(x)} = x/\mathrm{W}(x)$, so that every value is a rational function of $\mathrm{W}$, $U$ and $V$. The printer writes $\partial_U r = -r_s^2Ve^{-r/r_s}/r$ and its companion, and prefers the form of a value that does not divide by $U$ or $V$, which would be $0/0$ on a horizon.
- `novikov_comoving`: Novikov's solution as he writes it, $g_{rr} = (\partial_r R)^2$ with $R = (3c/2)^{2/3}F^{1/3}(t - t_0)^{2/3}$ and $8\pi G\rho/c^2 = F'/(R^2\partial_r R)$. Here $b = ct_0$ and $q = ct - b$, so no value holds the time outside $q$, and $\partial_r R = R(qF' - 2Fb')/3Fq$ makes every value a power of $R$ times a rational function.

The origin of $u$ is at the crossing because Kruskal's coordinates crowd the white hole out of sight otherwise: with the origin at the moment of rest the surface crosses the past horizon at $V = 0$, $U = -\sqrt{2}\,e^{1 + \pi/2} = -18.5$ and leaves the singularity at $U = -e^{\pi} = -23.1$, $V = -0.043$.

## Step 3. The junction

On the surface $\chi = \chi_0$ the areal radius is $R = a\sin\chi_0$, and with $r_s = a_m\sin^3\chi_0$ it is a radial geodesic of the exterior of energy $\cos\chi_0$.
$K^\theta{}_\theta = \cos\chi_0/R$ from both sides, so the surface carries no layer.
It crosses $r_s$ at $\eta = 2\chi_0$ and $2\pi - 2\chi_0$.
Along it on the way out $du/d\tau = 1/(E + \sqrt{E^2 - f})$, and for $\chi_0 = \pi/4$ the advanced time of the collapse's crossing is $v_h = (\pi + 2 + \ln 2)r_s$ in closed form, so the surface leaves $r = 0$ at $u = (2 + \ln 2 - \pi)r_s$ and rests at $u = (\pi + \ln 2)r_s$.
`_tools/test_white_hole.py` holds all of this to the published files.

## Step 4. The drawings

Every drawing takes the core that rests at $R_0 = 2r_s$, $\chi_0 = \pi/4$, as `oppenheimer_snyder` does.
The spacetime diagrams draw the surface before its moment of rest as well as after, `surface_whole` in `null_rays.py`, and Novikov's vacuole with $r_1 = r_s$, $r_2 = 4r_s$, whose singularity lies in the past of the drawing, which `crunch_curves` now finds from either side.
The conformal diagram is the collapse's after the moment of rest and its mirror image before it.
The embedding movie is the collapse's read backwards, from $c\tau = 0.06\,r_s$ to the rest.
Its moments are not marked on Kruskal's view, where the later ones lie a dozen widths of the drawing away, nor on Novikov's vacuole, another core than the one embedded.
