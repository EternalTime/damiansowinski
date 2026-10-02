# Bonnor's rotating dust cloud: the charts, the source and what the diagrams draw

Bonnor's cloud has one parameter, a length $a$ with $a^2 = 2h$ for Bonnor's rotation parameter $h$, and its angular momentum is $J = c^3a^2/2G$.
This note records where each chart of `bonnor_rotating_dust.json` comes from, what its source is, and what the diagrams take from it.
Units here have $G = c = 1$ unless a line keeps them; the file keeps $c$.

## Step 1. The cylindrical chart

Van Stockum's class of rigidly rotating dust is, in comoving coordinates and the signature $(-,+,+,+)$,

$$ds^2 = -\left(dt - K\,d\phi\right)^2 + \rho^2d\phi^2 + e^{\mu}\left(d\rho^2 + dz^2\right),$$

with $K(\rho, z)$ a solution of $K_{\rho\rho} - K_\rho/\rho + K_{zz} = 0$ and $\mu$ its quadrature, $\mu_\rho = (K_z^2 - K_\rho^2)/2\rho$ and $\mu_z = -K_\rho K_z/\rho$.
Bratek, Jałocha and Kutschera's (2) to (4) are these equations with $e^{2\Psi}$ for $e^\mu$, and Astesiano, Bini, Geralico and Ruggiero's (3.3) to (3.5) with $A = -K$.
Writing $K = \rho\,\partial_\rho\xi$ turns the first into Laplace's equation for $\xi$ in the flat space of $\rho$, $\phi$ and $z$, and Bonnor's cloud is $\xi = -2h/r$, the potential of a point, with $r = \sqrt{\rho^2 + z^2}$: Ilyas, Yang, Malafarina and Bambi's section 2.1.
Then

$$K = \frac{2h\rho^2}{r^3} = \frac{a^2\rho^2}{r^3}, \qquad \mu = \frac{h^2\rho^2(\rho^2 - 8z^2)}{2r^8} = \frac{a^4\rho^2(\rho^2 - 8z^2)}{8r^8},$$

which is Collas and Klein's (1) and (2) and Astesiano and his coauthors' (3.9) and (3.10) with $m = a$.
Bonnor, Collas and Klein, and Astesiano and his coauthors write $r$ for the distance from the axis and $R$ for the distance from the centre.
The file takes Bratek, Jałocha and Kutschera's and Ilyas and his coauthors' letters, $\rho$ and $r$, so that $R$ is free for the Ricci scalar, and takes $a$ for its parameter, the radius of the null circle in the equatorial plane and Bratek, Jałocha and Kutschera's dipole amplitude, $K_E^{(1)} = a^2\sin^2\theta/r$.
`print_charts.bonnor_rotating_dust("cylindrical")` builds it, with $r$ a name the chart defines.

## Step 2. The source

`bonnor_rotating_dust_source` holds the chart to van Stockum's three equations and to its source.
The Einstein tensor is $G^{tt} = D$ and nothing else, dust at rest in the chart, with

$$D = e^{-\mu}\,\frac{K_\rho^2 + K_z^2}{\rho^2} = \frac{a^4(\rho^2 + 4z^2)}{r^8}\,e^{-\mu},$$

which is $8\pi G/c^2$ times the density and equals the Ricci scalar: Collas and Klein's (4).
The four-velocity $\partial_t$ is a Killing vector of unit length, so the dust's world lines are geodesics that neither expand nor shear, and the square of their vorticity vector is $D/4$, $2\pi G$ times the density, as on the axis of van Stockum's cylinder.
The determinant is $-\rho^2e^{2\mu}$, so the signature holds wherever the metric is finite.
The Kretschmann scalar at $z = 0$ tends to $-36a^4/\rho^8$ far out and on the axis it is $4a^4(8a^4 - 27z^4)/z^{12}$, Astesiano and his coauthors' (3.15) and (3.16).

## Step 3. The mass

In the dust's rest space the volume element is $e^{\mu}\rho\,d\rho\,d\phi\,dz$, which cancels the $e^{-\mu}$ of the density, so the dust outside the sphere $r = r_0$ has the mass

$$\frac{c^2}{8\pi G}\int D\,e^{\mu}\rho\,d\rho\,d\phi\,dz = \frac{c^2a^4}{4G}\int_{r_0}^\infty\frac{dr}{r^4}\int_0^\pi(1 + 3\cos^2\theta)\sin\theta\,d\theta = \frac{c^2a^4}{3Gr_0^3}.$$

It grows without bound as $r_0 \to 0$.
Bratek, Jałocha and Kutschera's surface integral for the Komar mass inside the sphere, $(1/8\pi)\oint K\,\partial_rK\,d\theta\,d\phi/\sin\theta$, is $-a^4/3r_0^3$, minus the mass of the dust outside, at every $r_0$.
So the total mass is zero, as $g_{tt} = -1$ says, and the centre holds a negative mass without bound: Bonnor's reading, as de Araujo and Wang report it.
`_tools/test_bonnor_rotating_dust.py` holds both integrals to the published metric and Ricci scalar.

## Step 4. The closed timelike curves

$g_{\phi\phi} = \rho^2 - K^2$ vanishes on $r^3 = a^2\rho$, the surface Astesiano and his coauthors' (3.13) writes as $z = \pm\sqrt{a^{4/3}\rho^{2/3} - \rho^2}$ and Bratek, Jałocha and Kutschera write as $\rho = a\sin\alpha\sqrt{\sin\alpha}$, $z = a\cos\alpha\sqrt{\sin\alpha}$.
It is a torus that meets the equatorial plane in the circle $\rho = a$ and closes on the centre.
Inside it the circles of constant $t$, $\rho$ and $z$ are timelike.
Since $g(-\partial_\phi, \partial_t) = -K < 0$, the future directed sense round them is toward $-\phi$.

## Step 5. The spherical chart

With $\rho = r\sin\theta$ and $z = r\cos\theta$,

$$ds^2 = -\left(dt - \frac{a^2\sin^2\theta}{r}\,d\phi\right)^2 + r^2\sin^2\theta\,d\phi^2 + e^{a^4\sin^2\theta(9\sin^2\theta - 8)/8r^4}\left(dr^2 + r^2d\theta^2\right),$$

the form Bratek, Jałocha and Kutschera expand in multipoles and Ilyas and his coauthors use for their discs; Astesiano and his coauthors' (3.8) is the same with the latitude $\alpha = \pi/2 - \theta$.
`bonnor_rotating_dust_pullback` holds it to being the cylindrical chart pulled back in every slot.
Far out $g_{t\phi} \to a^2\sin^2\theta/r$, Kerr's frame dragging term $2GJ\sin^2\theta/c^3r$, with no mass term beside it.

## Step 6. The exponent

Both charts carry one exponential, $e^\mu$, and every value is printed with a single power of it, written as a rational multiple of the exponent the line element spells, by `one_exponential`.
The exponent stands over $8r^8$ with $r^2 = \rho^2 + z^2$, and `expand` multiplies the 8 into the sum, so the checker's `_exponent_terms` takes a number out of a denominator's sum before it names the generator; `_tools/test_roots.py` holds that.

## Step 7. What the diagrams draw

Every drawing is in units of $a$.
The axis $\rho = 0$ is fixed by the rotations and the plane $z = 0$ by the reflection $z \to -z$, so both are totally geodesic.
On the axis $K = \mu = 0$ and the metric is $-c^2dt^2 + dz^2$: its rays run at 45°, and the half axis is Minkowski's triangle in the conformal diagram, with the singularity a timelike line.
In the plane $z = 0$ a ray with no angular momentum has $c\,dt/d\rho = \pm e^{a^4/16\rho^4}\sqrt{1 - a^4/\rho^4}$ and turns as $d\phi = -a^2c\,dt/(\rho^3 - a^4/\rho)$; the spacetime diagram and the conformal diagram draw these with the circles of $\phi$ divided out, from $\rho = a$ outward, since inside it the circles are timelike and the quotient has no light cone.
The circles $\rho = a/2$ and $3a/2$ of that plane are drawn through time as van Stockum's cylinders are: the metric on each is $-c^2dt^2 + (2a^2/\rho)\,c\,dt\,d\phi + (\rho^2 - a^4/\rho^2)\,d\phi^2$ and its null lines are $c\,dt = (a^2/\rho \pm \rho)\,d\phi$.
The figure in three dimensions is the slice $z = 0$ with cones round the circles $a/2$, $a$ and $3a/2$, tipping toward $-\phi$.
The embedding diagram is the moment $t = 0$ of the plane $z = 0$ outside $\rho = a$.
The circle at $\rho$ has radius $R = \sqrt{\rho^2 - a^4/\rho^2}$, and with $x = a^4/\rho^4$, $(dR/d\rho)^2 = (1 + x)^2/(1 - x)$ exceeds $g_{\rho\rho} = e^{x/8}$ for every $0 < x < 1$, so no surface of revolution in flat space carries any part of the slice and the whole of it is drawn in Minkowski space, leaving the axis along the light cone and flattening far out, where $(dR/d\rho)^2 - g_{\rho\rho} \to 23a^4/8\rho^4$.
