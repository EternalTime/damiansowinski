# Black hole fireworks: the charts, the junction and the drawing

Hal M. Haggard and Carlo Rovelli, Phys. Rev. D 92, 104020 (2015), arXiv:1407.0989, read in its second arXiv version, which is the latest.
Carlos Barceló, Raúl Carballo-Rubio and Luis J. Garay, Int. J. Mod. Phys. D 23, 1442022 (2014), arXiv:1407.1391, and with Gil Jannes, Class. Quantum Grav. 32, 035012 (2015), arXiv:1409.1501.
Tommaso De Lorenzo and Alejandro Perez, Phys. Rev. D 93, 124018 (2016), arXiv:1512.04566.
Marios Christodoulou, Carlo Rovelli, Simone Speziale and Ilya Vilensky, Phys. Rev. D 94, 084035 (2016), arXiv:1605.05268.
Units: r_s = 2m = 1 and c = 1 wherever a number is drawn.

## Step 1: the charts and their sources

Every chart is a piece of Minkowski's or Schwarzschild's spacetime, and `black_to_white_hole_check` in `print_charts.py` holds each to the published chart of `minkowski` or `schwarzschild`, pulled back, and to a vanishing Ricci tensor.

| chart | source | what it covers |
| --- | --- | --- |
| `interior` | Haggard and Rovelli (24), De Lorenzo and Perez (9) | inside the shell, falling in (v < 0) and going out (u > 0) |
| `kruskal` | Haggard and Rovelli (3), (4), (26), with 32 m^3 = 4 r_s^3 | region II, and its time reverse (U, V) -> (-V, -U) |
| `schwarzschild` | Haggard and Rovelli (34) to (38); Christodoulou et al., Fig. 4 | each flap outside r_s |
| `painleve_gullstrand_ingoing` | Barceló, Carballo-Rubio and Garay (2), (3); with Jannes (1), (5) | the flap before the bounce |
| `painleve_gullstrand_outgoing` | the same with v_E = -v_I | the flap after the bounce |
| `lemaitre` | Christodoulou et al. (8), (9), r_S^3 = (9m/2)(r - t)^2 | the flap before the bounce |

The quantum region III has no chart.
Haggard and Rovelli's ansatz for it, their (31) and (32), F = (32 m^3/r_q) e^(r_q/2m) with r_q = (v_q - u_q)/2, does not join region II continuously: on their boundary the radius r_q differs from Kruskal's r, at Delta 0.73 against 7/6 in units of r_s, and r_q is a pure number where r is a length.
Every drawing leaves the region empty and says so.

Lemaître's values are printed as rational functions of r and rho - c tau, with each square of rho - c tau written 4r^3/(9 r_s): `black_to_white_hole_lemaitre` writes the derivatives of r by d_0 r = -2r/(3(rho - x^0)) and d_rho r = 2r/(3(rho - x^0)), which hold no root.
The checker holds r without declared rates and compares each value with r written out, since a root of rho - c tau, a sum that leads with its negative term in sympy's order, does not cancel against the same root written as sqrt(r_s/r).

## Step 2: the junction across the shell

Inside, ds^2 = -du dv + r^2 dOmega^2 with r = (v - u)/2, and the shell falls in along v = 0, where r = -u/2.
Outside, (1 - r) e^r = UV, and the shell is V = V0.
An outgoing ray crosses the shell with one radius on both sides, so (1 + u/2) e^(-u/2) = U V0, and

    U(u) = (1 + u/2) e^(-u/2) / V0.

Haggard and Rovelli's (28) to (30) print the exponent as e^(+u_I/4m); with r_I = (v_I - u_I)/2 on v_I = 0 their own (4) gives (1 + u_I/4m) e^(-u_I/4m), the form above at 4m = 2 r_s = 2.
U rises with u on u < 0, from -infinity to 1/V0 at the bounce, so every outgoing ray of region II continues one ray of the interior.
The apparent horizon U = 0 continues u = -2 r_s, which reaches the centre at ct = -r_s.

## Step 3: E, Delta and the surface of time symmetry

Haggard and Rovelli put E at (u_I, v_I) = (-2 eps, 0), at the radius eps, and Delta at r = 2m + delta on u + v = 0, with delta = m/3, so r(Delta) = 7/6 r_s and V(Delta) = sqrt((7/6 - 1) e^(7/6)) = 0.7316.
The shell crosses U + V = 0 where (1 - r) e^r = -V0^2, which must lie inside Delta, Christodoulou et al.'s (7), so V0 < 0.7316.
Schwarzschild's t outside r_s is U = -sqrt(r - 1) e^((r - t)/2), V = sqrt(r - 1) e^((r + t)/2), so the surface of time symmetry is t = 0.

## Step 4: the geodesic from Delta to E

The radial geodesics of -(4/r) e^(-r) dU dV obey U'' = -(d_U ln F) U'^2 and V'' = -(d_V ln F) V'^2, with d ln F/dr = -1/r - 1 and d_U r = -V e^(-r)/r.
Shot from Delta toward smaller V, a geodesic meets the shell V = V0 at a U that first rises with the angle of departure and then falls.
The largest U reached on the shell, and so the smallest radius at which E can stand, is

| V0 | largest U on the shell | radius there |
| --- | --- | --- |
| 0.7 | 1.419 | 0.112 r_s |
| 0.6 | 1.548 | 0.338 r_s |
| 0.5 | 1.654 | 0.497 r_s |
| 0.3 | 1.840 | 0.736 r_s |

so the smaller V0 is, the further out E must be for a spacelike geodesic to reach it.
Haggard and Rovelli's V0 = exp(-k m/2 l_P) is tiny and their eps Planckian, which this table does not reach.
The drawings take V0 = 3/5 and eps = r_s/2, where two geodesics reach E: one runs in to r = 0.012 r_s and back out, and the other, drawn, keeps its radius falling from 7/6 to 1/2.
Its tangent leaves Delta at dV/dU = -0.029, almost along the ingoing light ray through Delta.

## Step 5: the conformal diagram

Inside the shell p = arctan u and q = arctan v.
Outside it p = arctan u(U), continuous across the shell by Step 2, and q = q(V) with q(V0) = 0 and q(V) = -p(-V) for V >= V(Delta), so that U + V = 0 outside Delta is T = p + q = 0.
Between the shell and the ray through Delta, q = H(V) (V - V0)/(V(Delta) - V0), where H(V) is the value of q that would put the geodesic's point on V on T = 0; then q < H, the geodesic lies below T = 0, and E is at p = -pi/4, q = 0.
The time reverse is (p, q) -> (-q, -p).
The geodesic leaves Delta nearly along the ingoing ray, rising in Kruskal's T, and its mirror image leaves it nearly along the outgoing ray, falling, so in Kruskal's tangent plane at Delta the two flaps overlap and there is no room for the quantum region between them: q has a corner on the ray through Delta, where its slope drops, and the curves of constant r turn there.
Haggard and Rovelli's figure 4 opens the flaps the same way, and Christodoulou et al. call Delta the point where the map between spacetime and the Kruskal geometry bifurcates.
