# Point particles in three dimensions: the charts and what the diagrams draw

Gravity in three spacetime dimensions has no field outside its sources, so every chart of `point_particle_2plus1.json` but one is flat, and the physics sits in one number, $\alpha = 1 - 4Gm/c^2$, which removes a wedge of angle $2\pi(1 - \alpha) = 8\pi Gm/c^2$.
This note records where each chart comes from and what the diagrams take from it.
Units here have $c = 1$; the file keeps $c$.
Equation numbers are those of Staruszkiewicz (1963), of Deser, Jackiw and 't Hooft (1984), cited DJT, and of Deser and Jackiw (1988), each read from the paper itself.

## Step 1. The conical chart

DJT's (2.8b) and (2.8c) bring the space of one particle to $d\rho^2 + \rho^2d\theta'^2$ with $0 \le \theta' \le 2\pi\alpha$.
With the angle left periodic in $2\pi$, $\theta' = \alpha\phi$, that is

$$ds^2 = -dt^2 + dr^2 + \alpha^2r^2d\phi^2,$$

the form Gott gave the cosmic string in 1985, with $r$ the proper distance from the particle.
`print_charts.point_particle_2plus1("conical")` builds it and the checker's `Geometry` finds every curvature tensor zero.

## Step 2. The wedge and circumference charts

The wedge chart keeps $\theta = \alpha\phi$, so the line element is Minkowski's and the particle shows only in the range of the angle: Staruszkiewicz's (8) and (9), DJT's (2.8c), and Deser and Jackiw's (2.3).
The circumference chart takes the radius of the circle itself, $R = \alpha r$, so that $ds^2 = -dt^2 + dR^2/\alpha^2 + R^2d\phi^2$: Staruszkiewicz's (6), with his $e^N = 1/\alpha$, and Deser and Jackiw's (2.4), which they call the imbedded form.
`point_particle_check` pulls the conical chart back through $\phi = \theta/\alpha$ and through $r = R/\alpha$ and compares every slot exactly.

## Step 3. The isotropic chart and two bodies at rest

DJT's static solution for any number of particles, their (2.7), is $dl^2 = C\prod_n|\mathbf{r} - \mathbf{r}_n|^{-8Gm_n}(dx^2 + dy^2)$.
For one particle it is their (2.8a), $r^{-8Gm}(dr^2 + r^2d\theta^2)$, and $-8Gm = 2\alpha - 2$.
The file writes the coordinate radius $\rho$ and an arbitrary length $\ell$, so that every term balances: $(\rho/\ell)^{2\alpha - 2}(d\rho^2 + \rho^2d\phi^2)$, which is the conical chart at $r = (\ell/\alpha)(\rho/\ell)^\alpha$, DJT's (2.8b).
That map carries the power $\alpha$, so the pullback is checked at six random points in forty digits.

For two particles at $(x, y) = (\pm d, 0)$ the factor is $\Omega = (\rho_1/\ell)^{2\alpha_1 - 2}(\rho_2/\ell)^{2\alpha_2 - 2}$, which Staruszkiewicz had as his (17) by a Schwarz-Christoffel map.
The chart is checked to be flat, to be the isotropic chart about either particle when the other has no mass, and to fall off far away as one particle with $\alpha = \alpha_1 + \alpha_2 - 1$, which is DJT's (3.2), the masses add.
`conformal_cones` prints every value as a rational function times whole powers of the conformal factors, and the two bodies' connection as the sum of each particle's own part, $(\alpha_i - 1)$ times a rational function over $\rho_i^2$, since the two potentials add.

## Step 4. The moving particle

Staruszkiewicz's recipe for bodies in motion, and DJT's (5.6), is to cut the wedge in the particle's rest frame and boost.
With $x_1 = \gamma(x - \beta t)$, $t_1 = \gamma(t - \beta x)$ and the wedge about the negative $x_1$ axis, $|y| < -x_1\tan(\delta/2)$ with $\delta = 2\pi(1 - \alpha)$, two identified events share $t_1$ and $x_1$ and differ by the sign of $y$.
So they share $t$ and $x$ as well: in the frame where the particle moves the faces are identified at equal $t$, and the wedge is $|y| < \gamma(\beta t - x)\tan(\delta/2)$, of half angle $\arctan(\gamma\tan(\delta/2))$.
`point_particle_check` pulls the wedge chart back through the boost exactly, and `projections.moving_wedge` checks each of those statements on the events it draws.

## Step 5. Gott and Alpert's planet

Gott and Alpert (1984) give a planet of uniform density a spherical cap inside and the cone outside.
The cap is $a^2(d\chi^2 + \sin^2\chi\,d\phi^2)$, and the checker's Einstein tensor for it is $G^t{}_t = -1/a^2$ and nothing else: dust at rest with no pressure, of mass per unit area $\sigma = c^2/8\pi Ga^2$ under DJT's $G_{\mu\nu} = 8\pi GT_{\mu\nu}$.
The cap out to $\chi_0$ has area $2\pi a^2(1 - \cos\chi_0)$, so its mass is $(c^2/4G)(1 - \cos\chi_0)$ and $\cos\chi_0 = \alpha$: at the edge the circle $a\sin\chi_0$ and its rate of growth $\cos\chi_0$ are the cone's $\alpha r$ and $\alpha$ at $r = a\tan\chi_0$.
The formulas are this note's, checked by `point_particle_check`; of Gott and Alpert's paper only the abstract's statement was read.

## Step 6. The values drawn

Every diagram is drawn at $\alpha = 3/4$, a particle of mass $c^2/16G$, whose wedge is a right angle; eight of them at rest close space into the surface of a cube.

The plane of $t$ and the radius at fixed angle is flat in every chart, so its rays are $t \pm r$, $t \pm R/\alpha$ and $t \pm (\ell/\alpha)(\rho/\ell)^\alpha$, which `--verify` checks.
Two particles of $\alpha = 3/4$ at $x = \pm d$ with $\ell = d$ are drawn on the plane $y = 0$ through both, where $\Omega = |x^2/d^2 - 1|^{-1/2}$ and the proper distance from the midpoint is `null_rays._two_particles_distance`: $x\,F(\tfrac14, \tfrac12; \tfrac32; x^2)$ between them and a smooth quadrature beyond.
The particles are $B(\tfrac12, \tfrac34)\,d = 2.396\,d$ apart, which a test holds.
A ray that reaches a particle ends there, so the plane is drawn in two views, each ending on a particle.

The figure in three dimensions of the conical chart is the cosmic string's beam of light, `projections.string_rays`, at this deficit, and that of the moving chart is the wedge at $v = 3c/5$, $103°$ wide.
The conformal diagram is Minkowski's half diamond in the proper distance from the particle, one view for each chart of one particle, the planet's with its edge marked, and Minkowski's whole diamond for two bodies in the proper distance along the line through them.
The embedding diagram is the cone of half angle $\arcsin\alpha$ with the planet's cap at its apex, and the ideal particle's cone unrolling onto the plane, the cosmic string's `unrolling` with the particle's names.
