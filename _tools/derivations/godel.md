# The Gödel universe in Gödel's cylindrical coordinates

The Gödel universe carries a second chart beside its Cartesian one: the cylindrical coordinates Gödel wrote in 1949, in which the circles about one world line of the dust are coordinate lines, so that the closed timelike curves are circles of constant $t$, $r$ and $z$.
Its mathematics is written by `_tools/derivations/print_charts.py --metric godel`, and `verify_metrics.py --system godel/cylindrical` checks it in seconds.

## Step 1. The transformation

With $t_x$, $x$, $y$ and $z_x$ the Cartesian coordinates,

$$e^x = \cosh 2r + \cos\phi\sinh 2r,\qquad ye^x = \sqrt{2}\sin\phi\sinh 2r,\qquad t_x = 2t + \sqrt{2}\left(2\arctan\left(e^{-2r}\tan\tfrac{1}{2}\phi\right) - \phi\right),\qquad z_x = 2z,$$

which is Gödel's, with his time scaled to the Cartesian chart's.
Pulled back through it, the Cartesian line element

$$\frac{1}{2\omega^2}\left(-dt_x^2 + dx^2 - \tfrac{1}{2}e^{2x}dy^2 + dz_x^2 - 2e^x\,dt_x\,dy\right)$$

becomes

$$\frac{2}{\omega^2}\left(-dt^2 + dr^2 + \left(\sinh^2 r - \sinh^4 r\right)d\phi^2 - 2\sqrt{2}\sinh^2 r\,dt\,d\phi + dz^2\right).$$

`godel_pullback` in `print_charts.py` computes $J^TgJ$, with $g$ the published Cartesian metric and $J$ the Jacobian of the map, and simplifies its difference from the second form to zero in every slot, symbolically, in about seven seconds, before the chart is written; it names the slot that misses if one does.
The two charts' Kretschmann scalars are both $12\omega^4$ and their Ricci scalars both $-2\omega^2$, as scalars of one geometry must be.

## Step 2. Printing

The checker's `Geometry` hands back every value of this chart in exponentials of $r$.
`chart_printer.hyperbolic` writes each in $S = \sinh r$ and $C = \cosh r$, reduces it by $C^2 = 1 + S^2$, clears its denominator of $C$ with the conjugate, factors it, and writes every $1 + S^2$ as $C^2$ and $(S - 1)(S + 1)$ as $S^2 - 1$.
Every printed value is read back through the checker's `Reader` before it is written, as for every chart the script prints.

## Step 3. The critical radius

$g_{\phi\phi} = 2(1 - \sinh^2 r)\sinh^2 r/\omega^2$ vanishes at $\sinh r = 1$, $r = \ln(1 + \sqrt{2})$.
Inside it the circles of constant $t$, $r$ and $z$ are spacelike, on it null, and beyond it timelike.
The spacetime diagrams find it as the zero of the published $g_{\phi\phi}$ rather than from this formula.
