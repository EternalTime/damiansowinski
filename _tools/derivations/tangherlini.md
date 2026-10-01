# The Schwarzschild-Tangherlini black hole

Tangherlini's metric in $D$ spacetime dimensions is $ds^2 = -f\,c^2dt^2 + dr^2/f + r^2d\Omega_{D-2}^2$ with $f = 1 - (r_h/r)^{D-3}$.
A chart needs its coordinates named one by one, so the collection publishes $D = 5$ in three charts and $D = 6$ in one, and the prose carries the general $D$.

## Step 1: the parameter

Emparan and Reall write $f = 1 - \mu/r^{D-3}$ with $\mu = 16\pi GM/((D-2)\Omega_{D-2})$ at $c = 1$.
$\mu$ carries the dimension $L^{D-3}$, which changes with $D$, so the charts quote the horizon radius $r_h = \mu^{1/(D-3)}$, a length in every dimension.
With $\Omega_3 = 2\pi^2$ and $\Omega_4 = 8\pi^2/3$, $r_h^2 = 8GM/(3\pi c^2)$ in five dimensions and $r_h^3 = 3GM/(2\pi c^2)$ in six.

## Step 2: the angles

The unit sphere of each dimension is built on the one below: $d\Omega_3^2 = d\psi^2 + \sin^2\psi\,(d\theta^2 + \sin^2\theta\,d\phi^2)$ and $d\Omega_4^2 = d\chi^2 + \sin^2\chi\,d\Omega_3^2$.
So $\theta$ and $\phi$ mean what they mean in Schwarzschild's chart, and every plane drawn holds the new angles at $\pi/2$.

## Step 3: the tortoise coordinate

In five dimensions $1/f = r^2/(r^2 - r_h^2) = 1 + \tfrac{r_h}{2}\left(\tfrac{1}{r - r_h} - \tfrac{1}{r + r_h}\right)$, so $r_* = r + \tfrac{1}{2}r_h\ln|(r - r_h)/(r + r_h)|$, which vanishes at $r = 0$.
The Eddington-Finkelstein charts are $u = ct - r_*$ and $v = ct + r_*$, and $ds^2 = -f\,dv^2 + 2\,dv\,dr + r^2d\Omega_3^2$.
In six, $1/f = 1 + 1/(r^3 - r_h^3)$ and, at $r_h = 1$, $r_* = r + \tfrac{1}{3}\ln|r - 1| - \tfrac{1}{6}\ln(r^2 + r + 1) - \tfrac{1}{\sqrt 3}\left(\arctan\tfrac{2r + 1}{\sqrt 3} - \tfrac{\pi}{6}\right)$, which vanishes at $r = 0$ as well.

## Step 4: the curvature

Both charts are vacuum, so the Ricci tensor vanishes and the Weyl tensor is the Riemann tensor.
The Kretschmann scalar is $K = (D-1)(D-2)^2(D-3)\,r_h^{2D-6}/r^{2D-2}$: $12r_s^2/r^6$ at $D = 4$, $72r_h^4/r^8$ at $D = 5$ and $240r_h^6/r^{10}$ at $D = 6$, each checked against sympy as it is written.
The printer factors $r^n - r_h^n$, and every value is rewritten to keep it whole, read back and compared.

## Step 5: the diagrams

The surface gravity is $\kappa = c^2f'(r_h)/2 = (D-3)c^2/2r_h$, and the Kruskal coordinates $U = -e^{-\kappa u}$, $V = e^{\kappa v}$ give $UV = -e^{2\kappa r_*}$ outside.
Since $r_*(0) = 0$, $UV = 1$ at the singularity in both dimensions, and with $p, q = \arctan U, \arctan V$ it lies on $T = \pm\pi/2$: Kruskal and Szekeres's square, as in four dimensions.
In five dimensions $f$ has the two real roots $\pm r_h$ and the conformal generator's `Tower` takes both; in six the other two roots are complex and `TangherliniSix` carries the arctangent.

The slice of constant $t$ on the plane of $r$ and $\phi$ has $dz/dr = 1/\sqrt{(r/r_h)^{D-3} - 1}$.
At $D = 4$ that is Flamm's paraboloid, at $D = 5$ the catenoid $z = r_h\,\mathrm{arcosh}(r/r_h)$, and at $D = 6$ $z = z_\infty - 2\sqrt{r_h^3/r}\,{}_2F_1(\tfrac12, \tfrac16; \tfrac76; r_h^3/r^3)$ with $z_\infty = \tfrac13B(\tfrac16, \tfrac12)\,r_h = 2.4286\,r_h$.
The integrand falls as $r^{-(D-3)/2}$, so the height is finite for every $D \ge 6$ and unbounded for $D = 4$ and $5$.
