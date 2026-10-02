# Ori's time machine: the charts and what the diagrams draw

Amos Ori's vacuum core (2005) is published in `ori_time_machine.json` in three charts.
This note records where each comes from, the units, and the numbers the diagrams carry.
`print_charts.py` writes the mathematics and checks every statement below before it writes.

## Step 1. The core and its units

Ori's equation (1) is $ds^2 = dx^2 + dy^2 - 2\,dz\,dT + [f(x, y, z) - T]\,dz^2$, with $z$ periodic in $L$ and $f$ periodic in $z$.
The term $dz\,dT$ has to be an area and $f - T$ a pure number times nothing else, so with $z$ a pure number $T$ and $f$ are areas, as Misner's $T$ is in `misner.json`.
Then $a$ in $f = a(x^2 - y^2)/2$ and $e$ in $t = T - a(x^2 - y^2)/2 + e(x^2 + y^2)$ are pure numbers, which Ori's condition $e > (2e + a)^2$ needs.
No coordinate is a time, so no component carries a factor of $c$.

`ori_core_check` holds the chart `vacuum_core` to $\det g = -1$, to his equation (3), $R_{izjz} = -\tfrac{1}{2}f_{,ij}$ with every other component zero up to the symmetries, and to his footnote 1: the one Ricci component is $R_{zz} = -\tfrac{1}{2}(f_{,xx} + f_{,yy})$.
The chart leaves $f$ free in every tensor, so no component assumes his equation (2), and the core is a vacuum exactly where $f$ is harmonic in $x$ and $y$.

## Step 2. The foliation

The chart `foliation` is his equation (5), the core with $f = a(x^2 - y^2)/2$ in the time $t$.
`ori_pullback` checks it slot by slot to be the first chart pulled back along $T = t + a(x^2 - y^2)/2 - e(x^2 + y^2)$, to be a vacuum with $\det g = -1$, and to have his $g^{tt} = t + (2e - a)^2x^2 + (2e + a)^2y^2 - e(x^2 + y^2)$.
With $e > (2e + a)^2$ that $g^{tt}$ is negative for every $t < 0$, so those surfaces are spacelike, and at $t = 0$ it vanishes only on $x = y = 0$.
The circle $x = y = t = 0$ has $g_{zz} = 0$ and $\Gamma^t{}_{zz} = \Gamma^x{}_{zz} = \Gamma^y{}_{zz} = 0$ there, which is checked: it is the closed null geodesic $N$.
The inequality has solutions only for $a < 1/8$, since $4e^2 + (4a - 1)e + a^2 < 0$ needs $(4a - 1)^2 > 16a^2$.
Every diagram is drawn at $a = 1/16$ and $e = 1/8$, where $(2e + a)^2 = 25/256 < 1/8$.

## Step 3. The Brinkmann chart

Ori states that the class is locally a vacuum plane fronted wave and gives no map.
With $u = -2e^{-z/2}$ and $v = Te^{z/2}$ (Euler's $e$ here), $-2\,du\,dv = -2\,dz\,dT - T\,dz^2$ and $f\,dz^2 = 4f\,du^2/u^2$, so the core is the pp-wave $ds^2 = (4f/u^2)\,du^2 - 2\,du\,dv + dx^2 + dy^2$ on $u < 0$.
For Ori's example the profile is $2a(x^2 - y^2)/u^2$, a plane wave whose amplitude falls as $1/u^2$.
A circuit of $z$ is $u \to ue^{-L/2}$, $v \to ve^{L/2}$, a boost, and that plane wave is carried into itself by it, which is why the quotient exists.
The chart `brinkmann` is that form, checked by `ori_pullback` to be the first chart pulled back along $T = -uv/2$, $z = -2\ln(-u/2)$.
The circles of constant $T$, $x$ and $y$ are the hyperbolae $uv = -2T$, timelike where $uv < -a(x^2 - y^2)$.
On $x = y = 0$ the ray $v = 0$ is $N$, and $u$ is an affine parameter along it, so each circuit is shorter by $e^{-L/2}$ and $N$ is incomplete to the future, as Ori says.

## Step 4. The reader and the letter e

Ori's parameter is called $e$, and `verify_metrics.Reader` read every `e^` as the exponential.
It now does so only in a system that does not declare `e` as a name of its own, so the foliation's $4e^2$ is the square of the parameter, and the Brinkmann chart, which has no such parameter, writes its exponentials as `\exp`.
`_tools/test_reader_names.py` holds both readings.

## Step 5. The spacetime diagrams

The plane $x = y = 0$ is totally geodesic for $f = a(x^2 - y^2)/2$, since $\Gamma^x{}_{zz} = -ax/2$ and $\Gamma^y{}_{zz} = ay/2$ vanish on it, and its metric is Misner's, $-2\,dz\,dT - T\,dz^2$.
Each of the first two charts draws that cylinder and the cylinder at $x = 4\,\ell$, $y = 0$, at $L = 2\pi$.
On a cylinder of fixed $x$ and $y$ one null family keeps $z$ and the other keeps $(T - f)e^{z/2}$, or $(t - e(x^2 + y^2))e^{z/2}$, which `--verify` checks; at $x = 4\,\ell$ the circle is null at $T = f = \ell^2/2$ and at $t = ex^2 = 2\,\ell^2$.
Off the central circle the tipping family is null and not geodesic, turned toward larger $x$ by $\Gamma^x{}_{zz}$, and the captions say so.
The Brinkmann chart draws the plane of $u$ and $v$ through the central circle against $(v - u)/2$ and $(u + v)/2$, with the edges $u = -2$ and $u = -2e^{-\pi}$ of one copy.

## Step 6. The conformal diagram

`conformal.ori_time_machine` draws the plane $x = y = 0$ in the flat plane that covers it, with $p = \arctan u$ and $q = \arctan 2v$, one view per chart, at $L = 2$ as Misner's views are.
The coordinates cover the half $p < 0$; the null line $q = 0$ is $N$ in the quotient.
Each view carries a restriction, since off that plane the covering space is the plane wave of Step 3 and no picture of the whole is faithful.

## Step 7. The embedding diagram

On a surface of constant $t < 0$ the slice $y = 0$ has the metric $dx^2 + 2(2e - a)x\,dx\,dz + (ex^2 - t)\,dz^2$.
With $z = \phi - \tfrac{2e - a}{2e}\ln(ex^2 - t)$ it becomes $[1 - (2e - a)^2x^2/(ex^2 - t)]\,dx^2 + (ex^2 - t)\,d\phi^2$, a surface of revolution of radius $\sqrt{ex^2 - t}$ at $L = 2\pi$.
It lies in flat space where $(d\rho/dx)^2 \le g_{xx}$, which is $e^2x^2 + (2e - a)^2x^2 \le ex^2 - t$, true for every $x$ and every $t < 0$ when $e \ge e^2 + (2e - a)^2$: at $a = 1/16$, $e = 1/8$ that is $1/8 \ge 13/256$.
The movie runs from $t = -2\,\ell^2$ to $-0.1\,\ell^2$, the throat of radius $\sqrt{-t}$ closing toward $N$.
Its moments are the lines $t = t_k$ on the foliation's cylinders, $T = t_k$ on the central circle and $T = t_k - 3\ell^2/2$ at $x = 4\,\ell$ in the first chart, and the hyperbolae $uv = -2t_k$ on the Brinkmann plane, which `slices.checks()` holds to the pullbacks of Steps 2 and 3.
