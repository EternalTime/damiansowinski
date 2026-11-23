# The Roman ring of wormholes

Matt Visser, "Traversable wormholes: The Roman ring", Phys. Rev. D 55, 5212 (1997), gr-qc/9702043, with the single wormhole of Visser, Phys. Rev. D 47, 554 (1993), hep-th/9202090, and the two-wormhole Roman configuration of Visser, Phys. Rev. D 49, 3963 (1994), gr-qc/9311026.

## Step 1. The chart

`rest` is Minkowski space in the Cartesian coordinates $(t, x, y, z)$ of the rest frame of the mouths, $ds^2 = -c^2dt^2 + dx^2 + dy^2 + dz^2$.
Visser starts from "$N$ identical wormholes in otherwise flat Minkowski spacetime" with all mouths at rest with respect to each other (1997, p. 2), and models each as two timelike world lines identified, a spatial jump $L$ apart with a time shift $T$ (1997, Fig. 1).
In the rest frame of the mouths his single wormhole of 1993 is the identification $(t, 0, 0, 0) \equiv (t + T, 0, 0, \ell)$, eq. (3.1), where his $\ell$ is the distance between the two mouths, written $L$ here so that $\ell$ keeps the ring's meaning.
His wormholes of 1994 in general positions are $(t, \mathbf{x}^{\rm in}_k) \equiv (t + T_k, \mathbf{x}^{\rm out}_k)$, eq. (124).
The ring's crossings run from an exit mouth to the next entrance mouth and gain the time $T$, and the jump through each wormhole gives it back, so that once round, eq. (2), the closed geodesic has $s^2 = N^2(\ell^2 - c^2T^2)$: the entrance mouth at $t + T$ is the exit mouth at $t$, which the domain $\mathcal{A}_k(t + T) \sim \mathcal{B}_k(t)$ states for the spheres of radius $R$ about the two.
Visser works with $c = 1$; every relation on the page writes $cT$ where his is $T$, and $\hbar c$ where his is $\hbar$.
The domains state $cT \ll L$, $\ell \ll L$ and $R \ll \ell$, the regime of his estimates (1997, p. 2 and eq. (7)), and $\ell > cT$, where the closed geodesic round the ring is spacelike.
Visser gives no radius or angle for the polygon, only the lengths and its symmetry, so the chart states none.

Outside the mouths every Christoffel symbol and every curvature component vanishes, and the geodesics are straight lines.
On the identified spheres the Riemann tensor is a delta function of the jump in the second fundamental form, Visser 1994 eq. (91), with its traces giving a thin layer of exotic matter, Sec. III and eq. (83).

`print_charts.roman_ring_check` holds the chart to Minkowski's metric, and checks Visser's formulas where they meet.
The polynomials of Visser 1994 eq. (65), $P_N(z) = \sum_n \binom{N + n + 1}{2n + 1} z^n$, are $U_N(1 + z/2)$ for $N = 0$ to $8$, Chebyshev polynomials of the second kind.
The argument printed in eq. (66), $(z + 4)/2$, fails already at $N = 1$, where $P_1 = 2 + z$; with $1 + z/2$ the throat determinant $N^2P_{N-1}(2\ell/R)^{-2}$ of eqs. (115) and (116) is the $[N/U_{N-1}(1 + \ell/R)]^2$ of Visser 1997 eq. (6), checked at three values of $\ell/R$ for $N = 2$ to $6$, and it falls as $N^2(R/2\ell)^{2(N - 1)}$, eq. (7), when $\ell$ is long beside $R$.
The length of a closed curve once round, $N$ crossings taking times that add to $NT$, is stationary where they are equal, the condition Visser 1994 eq. (130) puts on the curve through two wormholes, and there $s^2 = N^2(\ell^2 - c^2T^2)$.

## Step 2. The tags

The chart is Minkowski's, so the facts computed from it would make the ring flat, a vacuum and conformally flat; `OVERRULED` in `metric_tags.py` sets all three false, since the mouths hold the delta function curvature and the exotic matter of Step 1, and the ring is tagged locally flat by hand.
The domain of the identification names $t$, so the chart alone does not read as stationary; a translation of $t$ carries each identified pair of spheres to another, so $\partial_t$ is a timelike Killing vector, which Visser 1993 Sec. IV uses for a wormhole whose mouths are at rest, and `OVERRULED` sets stationary true.
The identification shifts $t$ by $T$ at every wormhole, so the slices of constant $t$ do not close up round the ring, and the ring is not tagged static.

## Step 3. The figure in three dimensions

The figure is the slice $z = 0$ of $t$, $x$ and $y$, $t$ up, of four wormholes, the number Visser's Fig. 2 draws, in units of the crossing $\ell$ with $L = 4\ell$, $cT = 0.8\ell$ and $R = 0.15\ell$.
The mouths stand at the corners of a polygon of eight sides inscribed in a circle, alternately $L$ and $\ell$ long and turned into itself by a quarter turn, `projections.roman_ring_mouths`; the half angles $A$ and $B$ that the two kinds of side subtend satisfy $A + B = \pi/N$ and $L/\sin A = \ell/\sin B$.
The closed geodesic leaves each exit mouth at $t = 0$, reaches the next entrance at $t = T$, and is there the event at $t = 0$ on the next exit mouth.
`roman_ring_loop` checks every crossing against the published metric to have the interval $\ell^2 - c^2T^2$, positive at the drawn values, the whole $s^2 = N^2(\ell^2 - c^2T^2)$, the polygon's sides, and each jump to fall by $T$ over the distance $L$.
The future light cones stand where the crossings start, so each crossing shows just outside the cone at its start, at $1.25\,c$.

## Step 4. What is not drawn

No conformal diagram is drawn: every plane of $t$ and one direction is a plane of Minkowski space, and the identifications join it at each mouth to another place a time $T$ away, so the causal structure of the ring lies on a closed curve round a polygon, which no plane of two coordinates holds and the figure of Step 3 draws.
No embedding diagram is drawn: a slice of constant $t$ is the plane with the balls about the mouths removed, and each entrance sphere is joined to its exit sphere on another slice, a time $T$ away, so no surface of one moment closes over a mouth.
No curvature figure is drawn: every curvature scalar vanishes away from the mouths, and on them the curvature is a delta function with no profile to plot.
