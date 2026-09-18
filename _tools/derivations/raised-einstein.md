# Raising both indices on the Einstein tensor of the two rotating dust solutions

This is the working behind the `einstein_tensor.uu` blocks of `MFS/assets/data/metrics/godel.json` and `MFS/assets/data/metrics/stockum_dust.json`.

`_tools/derivations/audit-2026-09-18.md` groups these two entries together under the heading "the doubly raised Einstein tensor of the rotating dust entries", eight disagreements in all, four from each.
Both entries agreed with sympy on $G_{\mu\nu}$ and on $G^\mu{}_\nu$ and disagreed on $G^{\mu\nu}$, and in both the published $uu$ block was not symmetric.
A symmetric tensor stays symmetric under raising with a symmetric metric, so the asymmetry alone was proof that the published values were not the raised $ll$ values, whatever else they were.

Nothing is fixed here by adjusting what was published until it becomes symmetric.
Each $uu$ block is recomputed from that entry's own $G_{\mu\nu}$ and that entry's own $g^{\mu\nu}$, by

$$G^{\mu\nu} = g^{\mu\alpha}g^{\nu\beta}G_{\alpha\beta},$$

and the result is then checked three ways: against the entry's own $G^\mu{}_\nu$, against the invariant form the source of the field has, and against the entry's own Ricci scalar through the trace.
The $ll$ and $ul$ blocks were already right and are untouched.

All eight wrong values turn out to be one mistake made twice.
Write $[G]$ for the matrix of $G^\mu{}_\nu$ with $\mu$ the row and $\nu$ the column, and $[g^{-1}]$ for the matrix of $g^{\mu\nu}$.
The raising the $uu$ block is supposed to be, $G^{\mu\nu} = G^\mu{}_\alpha g^{\alpha\nu}$, is the matrix product

$$[G]\,[g^{-1}],$$

and every one of the eight wrong values is reproduced exactly, to the last term, by the other order,

$$[g^{-1}]\,[G].$$

That product is the inverse metric applied to the index that was already up, leaving the index that was down where it was.
It is worth noticing that it cannot be written in index notation at all: $g^{\mu\alpha}G^\alpha{}_\nu$ repeats $\alpha$ upstairs twice and so is not a contraction, which is exactly the error the notation exists to forbid and the reason it is stated here as a matrix product instead.
Step 4 shows why the slip is hard to catch once made: the two orders have the same trace, so the one cheap test a reader would apply passes on both.

---

## Step 1. Conventions

The signature is $(-,+,+,+)$ and the Ricci tensor is the standard contraction $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$, settled for the collection on 2026-09-18.
The Einstein tensor is $G_{\mu\nu} = R_{\mu\nu} - \tfrac{1}{2}Rg_{\mu\nu}$.

Neither entry needs a factor of $c$ carried anywhere, because neither chart has a coordinate that carries dimensions of time.
Godel's $e^x$ forces $x$ dimensionless, and the line element then forces the other three coordinates dimensionless with it, so the whole length scale of that solution sits in $1/\omega$.
Van Stockum's $t$ is already a length, declared as such in the checker's `DIMENSIONS` table, so the chart coordinate $x^0$ is that $t$ itself and not $ct$.
This is why nothing below picks up a $c$ on raising a time index, and it is worth saying out loud because raising a time index is exactly where a factor of $c$ would appear if the chart had one.

Raising an index does move the dimension.
With every chart coordinate a length except Godel's dimensionless four and van Stockum's $\phi$, an upper index contributes $[x^\mu]/L$ and a lower one $L/[x^\mu]$, on top of the $1/L^2$ the Einstein tensor carries when every coordinate is a length.
That arithmetic is what Step 3 uses to catch the two van Stockum components the dimensional pass flagged.

## Step 2. Godel

### The metric and its inverse

The entry publishes the line element

$$ds^2 = -\frac{1}{2\omega^2}dt^2 + \frac{1}{2\omega^2}dx^2 - \frac{e^x}{\omega^2}\,dt\,dy - \frac{e^{2x}}{4\omega^2}dy^2 + \frac{1}{2\omega^2}dz^2,$$

so in the order $(t,x,y,z)$ the metric and the published inverse are

$$g_{\mu\nu} = \frac{1}{2\omega^2}\begin{pmatrix} -1 & 0 & -e^x & 0 \\ 0 & 1 & 0 & 0 \\ -e^x & 0 & -\tfrac{1}{2}e^{2x} & 0 \\ 0 & 0 & 0 & 1\end{pmatrix}, \qquad g^{\mu\nu} = 2\omega^2\begin{pmatrix} 1 & 0 & -2e^{-x} & 0 \\ 0 & 1 & 0 & 0 \\ -2e^{-x} & 0 & 2e^{-2x} & 0 \\ 0 & 0 & 0 & 1\end{pmatrix}.$$

Only the $(t,y)$ block is nontrivial, and there

$$\det\begin{pmatrix} -1 & -e^x \\ -e^x & -\tfrac{1}{2}e^{2x}\end{pmatrix} = \frac{e^{2x}}{2} - e^{2x} = -\frac{e^{2x}}{2},$$

so the inverse of that block is $\dfrac{2\omega^2}{-e^{2x}/2}\begin{pmatrix} -\tfrac{1}{2}e^{2x} & e^x \\ e^x & -1\end{pmatrix} = 2\omega^2\begin{pmatrix} 1 & -2e^{-x} \\ -2e^{-x} & 2e^{-2x}\end{pmatrix}$, which is what the entry publishes.
The published inverse does invert the published metric, so the raising below is done with the entry's own matrix and not with a repaired one.

### The lowered Einstein tensor

The entry publishes $R_{tt} = 1$, $R_{ty} = R_{yt} = e^x$, $R_{yy} = e^{2x}$ and $R = -2\omega^2$, and from $G_{\mu\nu} = R_{\mu\nu} - \tfrac{1}{2}Rg_{\mu\nu}$,

$$G_{tt} = 1 - \tfrac{1}{2}(-2\omega^2)\!\left(-\tfrac{1}{2\omega^2}\right) = 1 - \tfrac{1}{2} = \tfrac{1}{2},$$
$$G_{ty} = e^x - \tfrac{1}{2}(-2\omega^2)\!\left(-\tfrac{e^x}{2\omega^2}\right) = e^x - \tfrac{e^x}{2} = \tfrac{e^x}{2},$$
$$G_{xx} = 0 - \tfrac{1}{2}(-2\omega^2)\!\left(\tfrac{1}{2\omega^2}\right) = \tfrac{1}{2}, \qquad G_{zz} = \tfrac{1}{2},$$
$$G_{yy} = e^{2x} - \tfrac{1}{2}(-2\omega^2)\!\left(-\tfrac{e^{2x}}{4\omega^2}\right) = e^{2x} - \tfrac{e^{2x}}{4} = \tfrac{3e^{2x}}{4},$$

which is the published $ll$ block exactly.
That block is the starting point and it does not move.

### Raising one index

Take it in two passes, so that the intermediate can be checked against what the entry already publishes.
With $G^\mu{}_\nu = g^{\mu\alpha}G_{\alpha\nu}$,

$$G^t{}_t = g^{tt}G_{tt} + g^{ty}G_{yt} = 2\omega^2\!\left(\tfrac{1}{2}\right) + \left(-4\omega^2e^{-x}\right)\!\left(\tfrac{e^x}{2}\right) = \omega^2 - 2\omega^2 = -\omega^2,$$
$$G^t{}_y = g^{tt}G_{ty} + g^{ty}G_{yy} = 2\omega^2\!\left(\tfrac{e^x}{2}\right) + \left(-4\omega^2e^{-x}\right)\!\left(\tfrac{3e^{2x}}{4}\right) = \omega^2e^x - 3\omega^2e^x = -2\omega^2e^x,$$
$$G^y{}_t = g^{yt}G_{tt} + g^{yy}G_{yt} = \left(-4\omega^2e^{-x}\right)\!\left(\tfrac{1}{2}\right) + 4\omega^2e^{-2x}\!\left(\tfrac{e^x}{2}\right) = -2\omega^2e^{-x} + 2\omega^2e^{-x} = 0,$$
$$G^y{}_y = g^{yt}G_{ty} + g^{yy}G_{yy} = \left(-4\omega^2e^{-x}\right)\!\left(\tfrac{e^x}{2}\right) + 4\omega^2e^{-2x}\!\left(\tfrac{3e^{2x}}{4}\right) = -2\omega^2 + 3\omega^2 = \omega^2,$$
$$G^x{}_x = g^{xx}G_{xx} = 2\omega^2\!\left(\tfrac{1}{2}\right) = \omega^2, \qquad G^z{}_z = \omega^2.$$

The entry's published $ul$ block is $G^t{}_t = -\omega^2$, $G^t{}_y = -2e^x\omega^2$, $G^x{}_x = G^y{}_y = G^z{}_z = \omega^2$, with $G^y{}_t$ absent and therefore zero.
Every one of the six agrees, which is the first confirmation that the $ll$ block and the inverse metric being used are the right pair.
Note that $G^\mu{}_\nu$ is not symmetric and is not supposed to be; the two index slots on it are of different kinds.

### Raising the second index

Now $G^{\mu\nu} = G^\mu{}_\alpha g^{\alpha\nu}$:

$$G^{tt} = G^t{}_t g^{tt} + G^t{}_y g^{yt} = \left(-\omega^2\right)\!\left(2\omega^2\right) + \left(-2\omega^2e^x\right)\!\left(-4\omega^2e^{-x}\right) = -2\omega^4 + 8\omega^4 = 6\omega^4,$$
$$G^{ty} = G^t{}_t g^{ty} + G^t{}_y g^{yy} = \left(-\omega^2\right)\!\left(-4\omega^2e^{-x}\right) + \left(-2\omega^2e^x\right)\!\left(4\omega^2e^{-2x}\right) = 4\omega^4e^{-x} - 8\omega^4e^{-x} = -4\omega^4e^{-x},$$
$$G^{yt} = G^y{}_t g^{tt} + G^y{}_y g^{yt} = 0 + \left(\omega^2\right)\!\left(-4\omega^2e^{-x}\right) = -4\omega^4e^{-x},$$
$$G^{yy} = G^y{}_t g^{ty} + G^y{}_y g^{yy} = 0 + \left(\omega^2\right)\!\left(4\omega^2e^{-2x}\right) = 4\omega^4e^{-2x},$$
$$G^{xx} = G^x{}_x g^{xx} = \left(\omega^2\right)\!\left(2\omega^2\right) = 2\omega^4, \qquad G^{zz} = 2\omega^4.$$

$G^{ty} = G^{yt} = -4\omega^4e^{-x}$, so the block is symmetric.
The two routes to the off-diagonal component are genuinely different arithmetic, $4\omega^4e^{-x} - 8\omega^4e^{-x}$ against $0 - 4\omega^4e^{-x}$, and they meet, which is the second confirmation.

### What was published before

| component | published | correct |
| --- | --- | --- |
| $G^{tt}$ | $-2\omega^4$ | $6\omega^4$ |
| $G^{ty}$ | $-8\omega^4\cosh x$ | $-4\omega^4e^{-x}$ |
| $G^{yt}$ | $4\omega^4e^{-x}$ | $-4\omega^4e^{-x}$ |
| $G^{yy}$ | $4(2 + e^{-2x})\omega^4$ | $4\omega^4e^{-2x}$ |
| $G^{xx}$ | $2\omega^4$ | $2\omega^4$ |
| $G^{zz}$ | $2\omega^4$ | $2\omega^4$ |

The two diagonal components in the $(x,z)$ plane were already right, because there both the inverse metric and $G^\mu{}_\nu$ are diagonal and the two orders of the product agree.
All four of the components that touch the $(t,y)$ block were wrong, and all four are entries of $[g^{-1}][G]$:

$$g^{tt}G^t{}_t + g^{ty}G^y{}_t = 2\omega^2\left(-\omega^2\right) + 0 = -2\omega^4,$$
$$g^{tt}G^t{}_y + g^{ty}G^y{}_y = 2\omega^2\left(-2\omega^2e^x\right) + \left(-4\omega^2e^{-x}\right)\omega^2 = -4\omega^4e^x - 4\omega^4e^{-x} = -8\omega^4\cosh x,$$
$$g^{yt}G^t{}_t + g^{yy}G^y{}_t = \left(-4\omega^2e^{-x}\right)\left(-\omega^2\right) + 0 = 4\omega^4e^{-x},$$
$$g^{yt}G^t{}_y + g^{yy}G^y{}_y = \left(-4\omega^2e^{-x}\right)\left(-2\omega^2e^x\right) + 4\omega^2e^{-2x}\omega^2 = 8\omega^4 + 4\omega^4e^{-2x} = 4\left(2 + e^{-2x}\right)\omega^4,$$

which are the four published values term for term.
The $\cosh$ is the signature of the slip: it is $-4\omega^4e^x$ picked up from the wrong index beside the $-4\omega^4e^{-x}$ that belongs there.

### The invariant form

Godel's source is dust with a cosmological constant, and the check that the recomputed block is right and not merely symmetric is that it has the form the source has.
The dust is at rest in this chart, $u^\mu = \sqrt{2}\,\omega\,\delta^\mu{}_t$, which is correctly normalised since

$$g_{\mu\nu}u^\mu u^\nu = -\frac{1}{2\omega^2}\left(\sqrt{2}\,\omega\right)^2 = -1,$$

and lowering gives $u_\mu = \left(-\dfrac{1}{\sqrt{2}\,\omega},\,0,\,-\dfrac{e^x}{\sqrt{2}\,\omega},\,0\right)$.
Writing $G_{\mu\nu} = A\,u_\mu u_\nu + B\,g_{\mu\nu}$ and reading off the $xx$ component gives $B = \omega^2$, and the $tt$ component then gives $A - B = \omega^2$, so $A = 2\omega^2$.
Those two constants were fixed by two components; the other three are then predictions, and

$$2\omega^2u_tu_y + \omega^2g_{ty} = 2\omega^2\frac{e^x}{2\omega^2} - \frac{e^x}{2} = \frac{e^x}{2}, \qquad 2\omega^2u_yu_y + \omega^2g_{yy} = e^{2x} - \frac{e^{2x}}{4} = \frac{3e^{2x}}{4},$$

with $G_{zz} = \omega^2g_{zz} = \tfrac{1}{2}$, all three of which match.
So

$$G_{\mu\nu} = 2\omega^2u_\mu u_\nu + \omega^2 g_{\mu\nu},$$

which is Einstein's equation with $\kappa\rho = 2\omega^2$ and $\Lambda = -\omega^2$, the values Godel's solution is known by.

Raising is now trivial, because $u_\mu$ raises to $u^\mu$ and $g_{\mu\nu}$ raises to $g^{\mu\nu}$ by definition:

$$G^{\mu\nu} = 2\omega^2u^\mu u^\nu + \omega^2 g^{\mu\nu}.$$

Since $u^\mu$ has only a $t$ component, $2\omega^2u^\mu u^\nu$ contributes $2\omega^2\left(\sqrt{2}\omega\right)^2 = 4\omega^4$ to $G^{tt}$ and nothing to anything else, and $\omega^2g^{\mu\nu}$ is the published inverse metric times $\omega^2$.
That gives

$$G^{tt} = 4\omega^4 + 2\omega^4 = 6\omega^4, \quad G^{ty} = G^{yt} = -4\omega^4e^{-x}, \quad G^{yy} = 4\omega^4e^{-2x}, \quad G^{xx} = G^{zz} = 2\omega^4,$$

the same six numbers, reached without raising anything component by component.
Symmetry here is not a check that passed but a property the expression cannot fail to have, since $u^\mu u^\nu$ and $g^{\mu\nu}$ are both symmetric.
This is the third confirmation.

### Dimensions

Godel's four coordinates are all dimensionless and $\omega$ carries $1/L$.
The Einstein tensor carries $1/L^2$ with every coordinate a length, an upper index contributes $[x^\mu]/L = 1/L$ here, so $G^{\mu\nu}$ carries $1/L^4$ in every slot and $\omega^4$ is exactly that.
Every published value, old and new, is a pure number times $\omega^4$, which is why the dimensional pass never flagged Godel: a dimensional check cannot see an error made in the coefficients alone.

## Step 3. Lanczos-van Stockum

### The metric and its inverse

$$ds^2 = -dt^2 - \frac{2r^2}{R}\,dt\,d\phi + e^{-r^2/R^2}dr^2 + r^2\!\left(1 - \frac{r^2}{R^2}\right)d\phi^2 + e^{-r^2/R^2}dz^2,$$

so in the order $(t,r,\phi,z)$ the nonzero components are $g_{tt} = -1$, $g_{t\phi} = g_{\phi t} = -r^2/R$, $g_{rr} = g_{zz} = e^{-r^2/R^2}$ and $g_{\phi\phi} = r^2\left(1 - r^2/R^2\right)$.
The $(t,\phi)$ block has determinant

$$(-1)\,r^2\!\left(1 - \frac{r^2}{R^2}\right) - \frac{r^4}{R^2} = -r^2 + \frac{r^4}{R^2} - \frac{r^4}{R^2} = -r^2,$$

which is independent of $R$, and inverting it gives

$$g^{tt} = \frac{r^2(1-r^2/R^2)}{-r^2} = -1 + \frac{r^2}{R^2}, \qquad g^{t\phi} = g^{\phi t} = \frac{r^2/R}{-r^2} = -\frac{1}{R}, \qquad g^{\phi\phi} = \frac{-1}{-r^2} = \frac{1}{r^2},$$

with $g^{rr} = g^{zz} = e^{r^2/R^2}$.
That is the published inverse metric, so again the raising uses the entry's own matrix.

### The lowered Einstein tensor

The entry publishes $R = 4e^{r^2/R^2}/R^2$ and, abbreviating $E = e^{r^2/R^2}$,

$$R_{tt} = \frac{2E}{R^2}, \quad R_{t\phi} = R_{\phi t} = \frac{2Er^2}{R^3}, \quad R_{rr} = R_{zz} = \frac{2}{R^2}, \quad R_{\phi\phi} = \frac{2Er^2(r^2+R^2)}{R^4}.$$

Then $G_{\mu\nu} = R_{\mu\nu} - \tfrac{1}{2}Rg_{\mu\nu}$ gives

$$G_{tt} = \frac{2E}{R^2} - \frac{2E}{R^2}(-1) = \frac{4E}{R^2}, \qquad G_{rr} = \frac{2}{R^2} - \frac{2E}{R^2}e^{-r^2/R^2} = \frac{2}{R^2} - \frac{2}{R^2} = 0,$$
$$G_{t\phi} = \frac{2Er^2}{R^3} - \frac{2E}{R^2}\!\left(-\frac{r^2}{R}\right) = \frac{4Er^2}{R^3}, \qquad G_{zz} = 0,$$
$$G_{\phi\phi} = \frac{2Er^2(r^2+R^2)}{R^4} - \frac{2E}{R^2}r^2\!\left(1 - \frac{r^2}{R^2}\right) = \frac{2Er^4}{R^4} + \frac{2Er^2}{R^2} - \frac{2Er^2}{R^2} + \frac{2Er^4}{R^4} = \frac{4Er^4}{R^4},$$

which is the published $ll$ block: four nonzero components, all in the $(t,\phi)$ block, and the $rr$ and $zz$ components cancelling exactly.
That cancellation is the statement that this dust has no pressure.

### Raising

$$G^t{}_t = g^{tt}G_{tt} + g^{t\phi}G_{\phi t} = \left(-1 + \frac{r^2}{R^2}\right)\frac{4E}{R^2} - \frac{1}{R}\cdot\frac{4Er^2}{R^3} = -\frac{4E}{R^2} + \frac{4Er^2}{R^4} - \frac{4Er^2}{R^4} = -\frac{4E}{R^2},$$
$$G^t{}_\phi = g^{tt}G_{t\phi} + g^{t\phi}G_{\phi\phi} = \left(-1 + \frac{r^2}{R^2}\right)\frac{4Er^2}{R^3} - \frac{1}{R}\cdot\frac{4Er^4}{R^4} = -\frac{4Er^2}{R^3} + \frac{4Er^4}{R^5} - \frac{4Er^4}{R^5} = -\frac{4Er^2}{R^3},$$
$$G^\phi{}_t = g^{\phi t}G_{tt} + g^{\phi\phi}G_{\phi t} = -\frac{1}{R}\cdot\frac{4E}{R^2} + \frac{1}{r^2}\cdot\frac{4Er^2}{R^3} = -\frac{4E}{R^3} + \frac{4E}{R^3} = 0,$$
$$G^\phi{}_\phi = g^{\phi t}G_{t\phi} + g^{\phi\phi}G_{\phi\phi} = -\frac{1}{R}\cdot\frac{4Er^2}{R^3} + \frac{1}{r^2}\cdot\frac{4Er^4}{R^4} = -\frac{4Er^2}{R^4} + \frac{4Er^2}{R^4} = 0.$$

The entry publishes exactly $G^t{}_t = -4E/R^2$ and $G^t{}_\phi = -4Er^2/R^3$ and nothing else, so the $ul$ block agrees.
Both of the components that vanish do so by a cancellation that the factor $r^2$ in $g^{\phi\phi}G_{\phi t}$ makes exact, and that cancellation is what drives the second raising to zero.

Raising the second index,

$$G^{tt} = G^t{}_tg^{tt} + G^t{}_\phi g^{\phi t} = -\frac{4E}{R^2}\left(-1 + \frac{r^2}{R^2}\right) - \frac{4Er^2}{R^3}\left(-\frac{1}{R}\right) = \frac{4E}{R^2} - \frac{4Er^2}{R^4} + \frac{4Er^2}{R^4} = \frac{4E}{R^2},$$
$$G^{t\phi} = G^t{}_tg^{t\phi} + G^t{}_\phi g^{\phi\phi} = -\frac{4E}{R^2}\left(-\frac{1}{R}\right) - \frac{4Er^2}{R^3}\cdot\frac{1}{r^2} = \frac{4E}{R^3} - \frac{4E}{R^3} = 0,$$
$$G^{\phi t} = G^\phi{}_tg^{tt} + G^\phi{}_\phi g^{\phi t} = 0, \qquad G^{\phi\phi} = G^\phi{}_tg^{t\phi} + G^\phi{}_\phi g^{\phi\phi} = 0.$$

So the whole doubly raised Einstein tensor of the van Stockum dust is one component,

$$G^{tt} = \frac{4e^{r^2/R^2}}{R^2},$$

and it is symmetric because there is nothing off the diagonal left to be asymmetric.

### The invariant form

That collapse is not an accident of the algebra, and the reason is the same one that organised Godel.
This dust is at rest in this chart: $u^\mu = \delta^\mu{}_t$ is normalised, since $g_{tt} = -1$, and lowering it gives $u_\mu = \left(-1,\,0,\,-r^2/R,\,0\right)$.
Then

$$\frac{4E}{R^2}u_\mu u_\nu$$

has $tt$ component $4E/R^2$, $t\phi$ component $4Er^2/R^3$, $\phi\phi$ component $4Er^4/R^4$ and nothing else, which is the published $ll$ block term for term.
So

$$G_{\mu\nu} = \frac{4e^{r^2/R^2}}{R^2}\,u_\mu u_\nu, \qquad \kappa\rho = \frac{4e^{r^2/R^2}}{R^2},$$

pure dust with no cosmological term, and raising both indices gives

$$G^{\mu\nu} = \frac{4e^{r^2/R^2}}{R^2}\,u^\mu u^\nu,$$

which has a $tt$ component and nothing else because $u^\mu$ has a $t$ component and nothing else.
The four published $uu$ components were three too many: the energy density of a comoving dust sits entirely in $G^{tt}$ when the indices are up, and the off-diagonal terms that appear with the indices down are there only because $u_\phi = -r^2/R$ does not vanish.
The omitted components are required to vanish by the checker, so dropping them from the file is itself checked rather than assumed.

### What was published before

| component | published | correct |
| --- | --- | --- |
| $G^{tt}$ | $-\dfrac{4E(r-R)(r+R)}{R^4}$ | $\dfrac{4E}{R^2}$ |
| $G^{t\phi}$ | $-\dfrac{4Er^2(r-R)(r+R)}{R^5}$ | $0$ |
| $G^{\phi t}$ | $\dfrac{4E}{R^3}$ | $0$ |
| $G^{\phi\phi}$ | $\dfrac{4Er^2}{R^4}$ | $0$ |

These four are entries of $[g^{-1}][G]$ as well, and here the slip is especially visible because $[G]$ has only its $t$ row nonzero, so the product collapses to $g^{\mu t}G^t{}_\nu$:

$$g^{tt}G^t{}_t = \left(-1 + \frac{r^2}{R^2}\right)\left(-\frac{4E}{R^2}\right) = -\frac{4E(r-R)(r+R)}{R^4}, \qquad g^{tt}G^t{}_\phi = \left(-1 + \frac{r^2}{R^2}\right)\left(-\frac{4Er^2}{R^3}\right) = -\frac{4Er^2(r-R)(r+R)}{R^5},$$
$$g^{\phi t}G^t{}_t = \left(-\frac{1}{R}\right)\left(-\frac{4E}{R^2}\right) = \frac{4E}{R^3}, \qquad g^{\phi t}G^t{}_\phi = \left(-\frac{1}{R}\right)\left(-\frac{4Er^2}{R^3}\right) = \frac{4Er^2}{R^4},$$

all four of which are the published values exactly.
The whole block is the column $g^{\mu t}$ times the row $G^t{}_\nu$, which is why it is a rank one matrix that is not symmetric: $g^{\mu t}$ and $G^t{}_\nu$ are not proportional.

### Dimensions

Here $t$, $r$, $z$ and $R$ are lengths and $\phi$ is dimensionless.
The Einstein tensor carries $1/L^2$ before the index positions are taken into account, an upper $t$ index contributes $[t]/L = 1$, and an upper $\phi$ index contributes $[\phi]/L = 1/L$.
So

$$[G^{tt}] = \frac{1}{L^2}, \qquad [G^{t\phi}] = \frac{1}{L^3}, \qquad [G^{\phi\phi}] = \frac{1}{L^4}.$$

The published $G^{t\phi} = -4Er^2(r-R)(r+R)/R^5$ expands to two terms, $-4Er^4/R^5$ and $4Er^2/R^3$, each carrying $1/L$ rather than $1/L^3$, and the published $G^{\phi\phi} = 4Er^2/R^4$ carries $1/L^2$ rather than $1/L^4$.
That is the three flagged terms standing for the two flagged components in the audit, and both components are now zero, which carries any dimension asked of it.

Which two got flagged and which two did not follows from the diagnosis above.
In $g^{\mu t}G^t{}_\nu$ the first index is raised twice and the second is never raised at all.
The first index comes out right anyway, because $[g^{\mu t}] = \dfrac{[x^\mu]}{L}\cdot\dfrac{[t]}{L} = \dfrac{[x^\mu]}{L}$ once $[t] = L$ is used, which is precisely what an upper $\mu$ is supposed to contribute.
The second index is wrong by $L^2$ whenever raising it would have changed anything, since a lower $\nu$ carries $L/[x^\nu]$ where an upper one carries $[x^\nu]/L$.
For $\nu = t$ those are both $1$ and nothing shows; for $\nu = \phi$ they differ by $L^2$, which is the gap the pass reports in both flagged components.
So the dimensional pass saw the half of this mistake that fell on a dimensionless coordinate and was blind to the other half, and Godel, whose coordinates are all dimensionless, it could not see at all.
The published $G^{tt}$ carried $1/L^2$ correctly and was still the wrong number, which is the audit's point that the dimensional pass is a filter and not a superset.

## Step 4. Why the error survived this long

The cheapest test there is on an Einstein tensor with both indices up is to contract it back down and compare against the Ricci scalar, since $G^\mu{}_\mu = R - 2R = -R$.
Both published blocks pass it.
For Godel,

$$g_{\mu\nu}G^{\mu\nu}_{\text{published}} = 2\omega^2 = -R,$$

and for van Stockum,

$$g_{\mu\nu}G^{\mu\nu}_{\text{published}} = -\frac{4e^{r^2/R^2}}{R^2} = -R,$$

and the recomputed blocks give the same two numbers.

This is not luck in either case.
Since $g_{\mu\nu}$ is symmetric, contracting any matrix $[M]$ of upper-upper components back down is $g_{\mu\nu}M^{\mu\nu} = \operatorname{tr}\left([g]\,[M]\right)$, and then

$$\operatorname{tr}\left([g]\,[G]\,[g^{-1}]\right) = \operatorname{tr}\left([G]\right) \qquad\text{and}\qquad \operatorname{tr}\left([g]\,[g^{-1}]\,[G]\right) = \operatorname{tr}\left([G]\right),$$

the first by cyclic invariance of the trace and the second because the two metrics cancel where they stand.
Both orders return $\operatorname{tr}[G] = G^\mu{}_\mu$, which is $-R$, and the $ul$ block they are both built from was correct.
The identity needs nothing of this spacetime, so raising the wrong index can never be caught by the trace in any entry in the collection.

The test that does separate the two orders is the one the audit used.
$[G][g^{-1}]$ is symmetric, because it is $[g^{-1}][G_{ll}][g^{-1}]$ with $[G_{ll}]$ and $[g^{-1}]$ both symmetric, while $[g^{-1}][G]$ is a product of two matrices that do not commute, and neither published block was symmetric.

## Step 5. What changed in the files

In `godel.json`, `cartesian.einstein_tensor.uu`: $G^{tt}$ from $-2\omega^4$ to $6\omega^4$, $G^{ty}$ from $-8\omega^4\cosh x$ to $-4e^{-x}\omega^4$, $G^{yt}$ from $4e^{-x}\omega^4$ to $-4e^{-x}\omega^4$, and $G^{yy}$ from $4(2 + e^{-2x})\omega^4$ to $4e^{-2x}\omega^4$.
$G^{xx}$ and $G^{zz}$ were already $2\omega^4$ and are untouched.

In `stockum_dust.json`, `cylindrical.einstein_tensor.uu`: $G^{tt}$ becomes $\dfrac{4e^{r^2/R^2}}{R^2}$, and the entries for $G^{t\phi}$, $G^{\phi t}$ and $G^{\phi\phi}$ are removed, because they are zero and the block lists nonzero components.

Nothing else in either file moved.
The `ll` and `ul` variants of both Einstein tensors, the Ricci tensors, the Ricci scalars, the Christoffel symbols, the Riemann tensors, the Weyl tensors corrected earlier the same day, and the geodesics Godel gained the same day are all as they were.
