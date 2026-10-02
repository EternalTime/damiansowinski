# Bonnor's beam of light

The four charts of `light_beam.json` are written by `print_charts.py --metric light_beam`, which computes every tensor from each line element with the checker's own `Geometry` and reads every printed value back before it writes.
This note derives what the charts are, why their tensors come out as printed, and what the diagrams draw.
The source of every chart is W. B. Bonnor, "The gravitational field of light", *Communications in Mathematical Physics* **13**, 163-174 (1969), cited below by its equation numbers.

## Step 1. The Cartesian chart

Bonnor's equation (2.13), in his signature $(+,-,-,-)$ and with $G = c = 1$, is $ds^2 = (-dx^2 - dy^2 - dz^2 + dt^2) + A\,(dt - dz)^2$ with $A \ge 0$, his (2.2).
In the signature $(-,+,+,+)$ with $c$ restored it is

$$ds^2 = -c^2dt^2 + dx^2 + dy^2 + dz^2 - A\left(c\,dt - dz\right)^2 ,$$

so his $A$ enters with a minus sign.
His section 4 takes $A = A(x, y)$, which makes the field stationary, and that is the free function the chart declares: a pure number, since it multiplies a squared length.
With $A$ depending on $u$ as well the metric is the general pp-wave, which is the chart of `pp_wave.json`.

Every component is taken in the chart $x^0 = ct$.
The metric and its inverse are written by hand as Minkowski's plus the beam, $g_{tt} = -1 - A$, $g_{tz} = A$, $g_{zz} = 1 - A$ and $g^{tt} = -1 + A$, $g^{tz} = A$, $g^{zz} = 1 + A$, and both are checked against sympy.
The inverse is linear in $A$ because $k = c\,dt - dz$ is null: $(\eta - A\,k \otimes k)^{-1} = \eta^{-1} + A\,k^\sharp \otimes k^\sharp$.

## Step 2. The null chart

Bonnor's (2.12) is $\sqrt{2}\,u = t - z$ and $\sqrt{2}\,v = t + z$, both lengths once $t$ is $ct$, and his (2.1) is $ds^2 = -dx^2 - dy^2 + 2\,du\,dv + 2A\,du^2$, which in this signature is

$$ds^2 = -2\,du\,dv + dx^2 + dy^2 - 2A\,du^2 .$$

$2\,du\,dv = c^2dt^2 - dz^2$ and $2\,du^2 = (c\,dt - dz)^2$, so it is the Cartesian chart's metric.
No coordinate is a time, so no component carries a factor of $c$ from the chart.
The inverse has $g^{uv} = -1$, $g^{vv} = 2A$ and $g^{uu} = 0$.

## Step 3. The Christoffel symbols and the geodesics

The only nonconstant component is $g_{uu} = -2A(x, y)$, so

$$\Gamma^x{}_{uu} = \partial_xA, \qquad \Gamma^y{}_{uu} = \partial_yA, \qquad \Gamma^v{}_{ux} = \partial_xA, \qquad \Gamma^v{}_{uy} = \partial_yA ,$$

which are Bonnor's (8.1) for a profile that does not depend on $u$.
No $\Gamma^u$ is nonzero, so $\ddot u = 0$, his (8.4), and $u$ is an affine parameter on every geodesic that is not a ray of the beam.
The geodesic equations are his (8.2) and (8.3): $\ddot x + \partial_xA\,\dot u^2 = 0$ and $\ddot v + 2(\partial_xA\,\dot x + \partial_yA\,\dot y)\dot u = 0$.
In the Cartesian chart every symbol multiplies $\dot t - \dot z$, the rate of $ct - z$, so the equations are printed grouped around it, and each is read back against the one computed from the Christoffel symbols.
A ray with $\dot t = \dot z$ and $\dot x = \dot y = 0$ solves them for every profile, his (8.5): light moving with the beam is undeflected.

## Step 4. The curvature and the source

The Riemann tensor of the null chart is $R_{uiuj} = \partial_i\partial_jA$ for $i, j = x, y$, and the Ricci tensor is $R_{uu} = \partial_x^2A + \partial_y^2A$, with no nonlinear term, which is why Bonnor finds the exact solution and the linear approximation identical (his section 4).
The Ricci scalar vanishes, since $g^{uu} = 0$, and so does the Kretschmann scalar, since every Riemann component carries two lower $u$ indices.
So $G_{\mu\nu} = R_{\mu\nu} = \tfrac{1}{2}\left(\partial_x^2A + \partial_y^2A\right)k_\mu k_\nu$ with $k = c\,dt - dz = \sqrt{2}\,du$, and Einstein's equations with $T_{\mu\nu} = \epsilon\,k_\mu k_\nu$, light of energy density $\epsilon$ moving along $+z$, give

$$\partial_x^2A + \partial_y^2A = \frac{16\pi G}{c^4}\,\epsilon ,$$

Bonnor's (2.10).
`light_beam_null_dust` in `print_charts.py` checks that in every slot of every chart, and that $k^\mu$ is null and covariantly constant, which is what makes the metric a pp-wave.
The Weyl tensor carries the trace free part of $\partial_i\partial_jA$, so it is of Petrov type N wherever it does not vanish.

## Step 5. The uniform beam

Bonnor's (5.6) and (5.5), with $r$ the distance from the axis, $a$ the beam's radius and $m = \pi\rho_0a^2$ the energy per unit length, are $A_i = 4mr^2/a^2$ for $r \le a$ and $A_e = 4m + 8m\ln(r/a)$ for $r \ge a$.
With $G$ and $c$ restored $m$ is the pure number $\pi G\epsilon R^2/c^4$, where $\epsilon$ is the energy density and $R$ the radius, so

$$A = \frac{4\pi G\epsilon}{c^4}\,\rho^2 \quad (\rho \le R), \qquad A = \frac{4\pi G\epsilon R^2}{c^4}\left(1 + 2\ln\frac{\rho}{R}\right) \quad (\rho \ge R),$$

and $g_{uu} = -2A$ is what the two null cylindrical charts print, in the polar coordinates of his (8.10).
$A$ and $dA/d\rho$ are continuous at $\rho = R$, both $4\pi G\epsilon R^2/c^4$ and $8\pi G\epsilon R/c^4$, his (5.4).
Inside, $\partial_i\partial_jA$ is a multiple of $\delta_{ij}$, so the Weyl tensor vanishes and the beam's interior is conformally flat, with $G_{uu} = 16\pi G\epsilon/c^4 = (8\pi G/c^4)\,\epsilon\,k_uk_u$.
Outside, $\ln\rho$ is harmonic, the Ricci tensor vanishes and the Weyl tensor equals the Riemann tensor.

The radial equation is $\ddot\rho + \Gamma^\rho{}_{uu}\dot u^2 - \rho\dot\phi^2 = 0$ with $\Gamma^\rho{}_{uu} = dA/d\rho$, his (8.11): $8\pi G\epsilon\rho/c^4$ inside, a pendulum, and $8\pi G\epsilon R^2/(c^4\rho)$ outside, the pull of a line of matter in Newton's theory.
With no angular momentum, $\rho'^2 + 2A$ is constant along a geodesic, a prime a derivative along $u$, which is his (8.16).

## Step 6. The spacetime diagrams

Every diagram is drawn at $m = 1/32$ in units of $R$, so $A = \rho^2/8$ inside the beam and $(1 + 2\ln\rho)/8$ outside it, $1/8$ at the edge.
On a plane of the time and $z$ at a fixed place across the beam $A$ is a constant, and the null condition $-2\,du\,dv - 2A\,du^2 = 0$ has the roots $du = 0$, a ray with the beam, and $dv = -A\,du$, a null curve against it, which keeps $v + Au$, or $ct + z + A(ct - z)$.
`--verify` holds the rays of each view to those two quantities.
The curve against the beam is a null geodesic only where the gradient of $A$ vanishes, his (8.6): on the axis of the uniform beam, and midway between two beams.
A constant added to $A$ is the change of coordinates $v \to v + \text{const}\cdot u$, so the tilt of a cone on one plane belongs to the chart and only differences of $A$ between two places do not.

The Cartesian chart's free profile is the uniform beam, drawn on its axis and at $x = 3R$.
The null Cartesian chart's is two uniform beams with their axes at $x = \pm 2R$, the sum of the two profiles, his (6.1) with each beam's interior: midway between them $A = (1 + 2\ln 2)/4$ and its gradient vanishes by symmetry, and on the axis of one $A = (1 + 2\ln 4)/8$ is the other's alone.
The null cylindrical charts are drawn at $\rho = R$, $2R$ and $4R$, where $g_{uu} = -1/4$, $-(1 + 2\ln 2)/4$ and $-(1 + 2\ln 4)/4$.
The box runs from $ct = -2R$ to $14R$ so that it holds the six wave fronts of the embedding diagram.

## Step 7. Light sent against the beam

The figure `lens` is the plane $y = 0$ of the Cartesian chart seen from the side with $t$ left out, which the reflection $y \to -y$ keeps every geodesic launched in it inside.
A ray sent with the beam starts at $z = -7R$ with $dz = c\,dt$ and is checked to keep its $x$.
A ray sent against it starts at $z = 7R$ with $dx = 0$ along the other null direction of the plane of $t$ and $z$, $c\,dt = -(1 - A)\,dz/(1 + A)$, and is integrated in $t$, $x$ and $z$ with the published Christoffel symbols.
It is checked null against the published metric, to keep $\dot t - \dot z$, and to keep $\dot x^2 + A(\dot t - \dot z)^2$, which is (8.16) with $2k^2 = (\dot t - \dot z)^2$.
Inside the beam $x'' = -x/4R^2$ along Bonnor's $u$, so a ray that starts at $x = b$ and stays inside is $b\cos(u/2R)$, checked to $10^{-8}$, and every such ray reaches the axis when $u$ has grown by $\pi R$: the beam focuses light sent the other way.

## Step 8. The conformal diagrams

Where no Christoffel symbol turns a ray of a plane of $u$ and $v$ out of it, the plane is totally geodesic and its metric $-2\,du\,d(v + Au)$ is flat, so $p = \arctan(\sqrt{2}\,u/\ell)$ and $q = \arctan(\sqrt{2}(v + Au)/\ell)$ bring it into Minkowski's diamond.
That holds on the axis of the uniform beam, where $A = 0$, drawn in the Cartesian chart and in the null cylindrical chart inside the beam, and midway between the two beams, where $A = (1 + 2\ln 2)/4$; `conformal.py` checks each map null against the published metric and each plane against the published Christoffel symbols.
Outside a single beam $dA/d\rho$ vanishes nowhere, so that chart has no view.

## Step 9. The embedding diagram

A surface of constant $u$ and $v$ has the metric $dx^2 + dy^2$, flat, so each moment is a flat disc with the beam's edge marked on it, and the beam shows in what it does to free particles.
360 particles at rest at $u = 0$ on the circle $\rho = 2R$ are run with the published $\Gamma^x{}_{uu}$ and $\Gamma^y{}_{uu}$ of the null Cartesian chart and the uniform beam's profile.
By Step 5 the ring stays a circle with $\rho'^2 = 2(A(2R) - A(\rho))$, which outside the beam is $\ln(2R/\rho)/2$ in units of $R$.
With $\rho = 2e^{-s^2}$ the time to reach $\rho$ is $u = 4\sqrt{2}\int_0^s e^{-s^2}ds = 2\sqrt{2\pi}\,\mathrm{erf}\sqrt{\ln(2/\rho)}$, so the ring reaches the edge at $u_1 = 2\sqrt{2\pi}\,\mathrm{erf}\sqrt{\ln 2} = 3.815$ with the speed $v_1 = \sqrt{\ln 2/2} = 0.589$.
Inside, $\rho = \cos(s/2) - 2v_1\sin(s/2)$ with $s = u - u_1$, which vanishes at $s = 2\arctan(1/2v_1)$: the ring closes on the axis at $u = 5.223$.
Every particle passes through the axis to the far side, the ring opens again as it closed, and it is back at rest on its first circle at $u = 10.446$, each particle opposite where it started.
The script checks the edge, the focus, the interior form, the constancy of $\rho'^2 + 2A$ and the return to rest, and the tests hold the rings in the file to the same closed forms.
The moments are $u = 0$, $2$, $4$, $6$, $8$ and $10$, each the null hypersurface $u = u_k$, and they are marked as null lines on the spacetime and conformal diagrams of the single beam.
Past the focus a ring is read round from the particle that then stands on the positive $x$ axis, so the tube's semi-axes stay positive and its twelve marked lines are still the world lines of the twelve marked particles, each crossing the axis.
