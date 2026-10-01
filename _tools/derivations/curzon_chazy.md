# The Curzon-Chazy particle

The Curzon-Chazy particle in Weyl's canonical chart and in the spherical chart $\rho = r\sin\theta$, $z = r\cos\theta$ is

$$ds^2 = -e^{-2m/R}c^2dt^2 + e^{2m/R}\left(e^{-m^2\rho^2/R^4}\left(d\rho^2 + dz^2\right) + \rho^2d\phi^2\right),\qquad R = \sqrt{\rho^2 + z^2},$$

$$ds^2 = -e^{-2m/r}c^2dt^2 + e^{2m/r}\left(e^{-m^2\sin^2\theta/r^2}\left(dr^2 + r^2d\theta^2\right) + r^2\sin^2\theta\,d\phi^2\right).$$

Both charts are written by `_tools/derivations/print_charts.py --metric curzon_chazy`, and `verify_metrics.py --system curzon_chazy/weyl` and `--system curzon_chazy/spherical` check each in seconds.

## Step 1. Weyl's two functions

A static, axially symmetric vacuum metric is $-e^{2\psi}c^2dt^2 + e^{-2\psi}\left(e^{2\gamma}(d\rho^2 + dz^2) + \rho^2d\phi^2\right)$, with $\psi$ harmonic in the flat space of $\rho$, $\phi$ and $z$ and $\gamma$ its quadrature, $\partial_\rho\gamma = \rho\left((\partial_\rho\psi)^2 - (\partial_z\psi)^2\right)$ and $\partial_z\gamma = 2\rho\,\partial_\rho\psi\,\partial_z\psi$.
Newton's potential of a point mass, $\psi = -m/R$, gives $\gamma = -m^2\rho^2/2R^4$, which vanishes on the axis, so the axis carries no conical singularity.
`curzon_chazy_check` in `print_charts.py` checks Laplace's equation, both quadratures and the vanishing of the Ricci tensor before anything is written, and `curzon_chazy_pullback` checks that the spherical chart is the pullback of Weyl's.
$m = GM/c^2$ is a length, and Abdelqader and Lake show that $M$ is the mass at spatial infinity.

## Step 2. A name the chart defines

Every value of Weyl's chart is a polynomial in $\rho$, $z$ and $R$, so the chart declares `R = \sqrt{\rho^2 + z^2}` among its parameters and prints around it.
`Reader._declare_parameter` in `verify_metrics.py` reads a declaration whose right side is no function of coordinates as a definition, and every published `R` as that expression, so the checker compares the same values it would with the root written out.
`DIMENSIONS` declares $[R] = L$, and `Dimensions` checks that the definition carries the declared dimension.
`weyl_distance` in `print_charts.py` writes each power of $\rho^2 + z^2$ as a power of $R$ and each even power of $z$ left over as $R^2 - \rho^2$, and gathers the exponentials as the line element writes them.
The Ricci scalar is also written $R$, as it is for Tolman-Bondi, whose areal radius has the same letter; the page prints it as $R = 0$.

## Step 3. The curvature

The Kretschmann scalar is $16m^2\left(m^4\rho^2 - 3m^3\rho^2R + 3m^2\rho^2R^2 + 3m^2R^4 - 6mR^5 + 3R^6\right)e^{-4m/R + 2m^2\rho^2/R^4}/R^{12}$, eight times Abdelqader and Lake's $w1R$.
On the axis it is $48m^2(|z| - m)^2e^{-4m/|z|}/z^8$, which vanishes at $z = \pm m$ and goes to zero as $z \to 0$; in the plane $z = 0$ it grows as $e^{2m^2/\rho^2}$.
So the curvature at $R = 0$ depends on the direction of approach, as Gautreau and Anderson found.
The circle of Weyl's radius $\rho$ in the plane $z = 0$ has circumference $2\pi\rho\,e^{m/\rho}$, least at $\rho = m$ and unbounded as $\rho \to 0$, and Stachel's surfaces of constant $t$ and $r$ have least area at $r = 0.5389\,m$, which a quadrature of $2\pi r^2e^{2m/r}\int_0^\pi e^{-m^2\sin^2\theta/2r^2}\sin\theta\,d\theta$ confirms.

## Step 4. The two planes

The rotations about the axis fix the axis, and the reflections $z \to -z$ and $\phi \to -\phi$ fix the half plane $z = 0$, $\phi = 0$, so both are totally geodesic and their null curves are null geodesics.
On the axis the metric is $-e^{-2m/z}c^2dt^2 + e^{2m/z}dz^2$ for $z > 0$, with $g_{tt}g_{zz} = -1$, so $z$ is an affine parameter along every ray, and $z_* = z\,e^{2m/z} - 2m\,\mathrm{Ei}(2m/z)$ falls to $-\infty$ as $z \to 0$: a ray reaches $z = 0$ after a finite affine distance and at $t = \pm\infty$.
In the plane the metric is $-e^{-2m/\rho}c^2dt^2 + e^{2m/\rho - m^2/\rho^2}d\rho^2$, and $\rho_* = \int_0^\rho e^{2m/s - m^2/2s^2}ds$ is finite at $\rho = 0$, so a ray reaches the ring in a finite time; $c\,dt/d\rho$ is greatest, $e^2$, at $\rho = m/2$.

## Step 5. The diagrams

The spacetime diagrams draw the two planes in each chart at $m = 1$, and `null_rays.py --verify` compares their rays with $ct \pm z_*$ and $ct \pm \rho_*$.
The conformal diagrams take $p, q = \arctan((ct \mp x_*)/\ell)$ with $\ell = 4m$: the half axis is the whole diamond, null infinity on the right and the edge $z = 0$ on the left, drawn as the edge of the chart with the Kretschmann scalar checked to go to zero on it, and the half plane is Minkowski's triangle with the ring on $X = 0$, drawn as a singularity with the Kretschmann scalar checked to diverge.
The embedding diagram is the plane $z = 0$ at $t = 0$: $g_{\rho\rho} \ge \left(\partial_\rho(\rho\,e^{m/\rho})\right)^2$ reduces to $e^{-m^2/\rho^2} \ge (1 - m/\rho)^2$, which holds from $\rho = 0.7226\,m$ out, so the surface is a funnel with its neck, of radius $m\,e$, at $\rho = m$, widening again below the neck until it lies level where the drawing stops.
The plane's moment is marked on the drawings of the plane and on none of the axis, which the plane meets only at $\rho = 0$; `slices.HIDDEN` records it.
