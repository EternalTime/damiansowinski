# The Lentz soliton

The Lentz soliton takes the line element Lentz wrote, in Cartesian coordinates, and each piece of it traces to his paper.

The components were computed and written by `_tools/derivations/print_charts.py`, and `print_charts.py` takes about eleven minutes over it, most of them grouping the Weyl tensor, and `_tools/derivations/verify_metrics.py --system lentz/cartesian` recomputes everything from the line element in about fifty seconds.

---

## Step 1. Lentz's form

Lentz (Classical and Quantum Gravity 38, 075015, 2021; arXiv 2006.07125) writes every spacetime he considers in the 3+1 form

$$ds^2 = -\left(N^2 - N^iN_i\right)dt^2 - 2N_i\,dx^i\,dt + h_{ij}\,dx^i\,dx^j,$$

in units with $G = c = 1$, and then fixes three things.
The lapse is one, $N = 1$.
The slices are flat and the coordinates Cartesian, $h_{ij} = \delta_{ij}$.
The shift is the gradient of a potential, $N_i = \partial_i\phi$, which is the first line of his section on constructing positive energy solutions and the step that defines the class.

With those three the line element is

$$ds^2 = -c^2dt^2 + \left(dx - \partial_x\phi\,c\,dt\right)^2 + \left(dy - \partial_y\phi\,c\,dt\right)^2 + \left(dz - \partial_z\phi\,c\,dt\right)^2,$$

with $c$ restored and $\phi$ a length so that each $N_i$ is a pure number.
It is Natário's Cartesian flow chart with the flow $\vec V = c\,\nabla\phi$, and the two can be read against each other.

The rest of Lentz's construction is a choice of $\phi$, and none of it is assumed.
His solitons move rigidly along $z$, $\phi = \phi(x, y, z - z_s(t))$, which in this chart is $\partial_t\phi = -v_s\,\partial_z\phi$.
He ties the second derivatives together with the wave equation $\partial_x^2\phi + \partial_y^2\phi - (2/v_h^2)\partial_z^2\phi = \rho_h$ on each slice and builds $\phi$ from rhomboid cells of the source $\rho_h$ through its Green's function; that $\phi$ exists only as a numerical integral, which is why the line element carries the class and not a member of it.
The potential is left a function of all four coordinates, as Alcubierre's shape function is left free, so the components also hold while a soliton is built, accelerated or taken apart.

## Step 2. Lentz's three statements, as components

Every identity in this step was checked in sympy for an arbitrary $\phi(t, x, y, z)$.

The Eulerian observers have $n_\mu = (-1, 0, 0, 0)$ and $n^\mu = (1, \partial_x\phi, \partial_y\phi, \partial_z\phi)$, and $G^{\mu\nu}n_\mu n_\nu = G^{tt}$.
For any potential

$$G^{tt} = \sigma_2(H) = H_{xx}H_{yy} - H_{xy}^2 + H_{xx}H_{zz} - H_{xz}^2 + H_{yy}H_{zz} - H_{yz}^2,\qquad H_{ij} = \partial_i\partial_j\phi,$$

the sum of the principal minors of the Hessian.
With $K_{ij} = -H_{ij}$ it is $\frac{1}{2}\left(K^2 - K^i{}_jK^j{}_i\right)$, Lentz's $8\pi E$.

The momentum density is $-n_\mu G^\mu{}_i = G^t{}_i$, and it vanishes identically: the momentum constraint is $\partial_jK^j{}_i - \partial_iK$, which is $-\partial_j\partial_j\partial_i\phi + \partial_i\partial_j\partial_j\phi = 0$ for a gradient.
That is Lentz's $J_i = 0$, and the mixed Einstein tensor has no $G^t{}_x$, $G^t{}_y$ or $G^t{}_z$.

## Step 3. The Kretschmann scalar is a sum of squares

Split the Riemann tensor in the frame of the Eulerian observers, $n$ and $\partial_x$, $\partial_y$, $\partial_z$, which is orthonormal here.
The spatial components are fixed by the Gauss equation on a flat slice, $R_{ijkl} = K_{ik}K_{jl} - K_{il}K_{jk}$, so its nonzero components are the cofactors $C_{ab}$ of the Hessian.
The components mixing $n$ with three spatial indices are fixed by the Codazzi equation, $\partial_jK_{ik} - \partial_kK_{ij}$, and vanishes because third derivatives of $\phi$ commute.
The tidal field is

$$E_{ij} = -n^an^cR_{aicj} = \partial_t H_{ij} + N^k\partial_kH_{ij} + H_{ik}H_{kj},$$

the rate of change of the Hessian along the observers' flow plus its square.
With the mixed components gone, every term of $R_{abcd}R^{abcd}$ with an odd number of $n$ indices is absent, and

$$K = 4\sum_{a,b}C_{ab}^2 + 4\sum_{i,j}E_{ij}^2,$$

written with the off diagonal terms doubled.
It is a sum of squares, so it vanishes only where the space is flat, and `print_charts.py` checks it against sympy's full contraction before writing it.

## Step 4. How the values are written

Every value is a polynomial in the derivatives of $\phi$, with no denominator, because the inverse metric of a unit lapse flat slice metric is polynomial in the shift.
`chart_printer.py` collects each value by the first derivatives $\partial_x\phi$, $\partial_y\phi$ and $\partial_z\phi$, which are the components of the flow, and factors what multiplies them.

## Step 5. Dimensions

`DIMENSIONS` declares $t$ a time, $x$, $y$ and $z$ lengths, and $\phi$ a length.
$\partial_t$ is along $x^0 = ct$, so every derivative of $\phi$ of order $k$ carries $L^{1-k}$, and every term balances.
