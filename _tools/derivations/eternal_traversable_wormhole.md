# The eternal traversable wormhole of Maldacena and Qi

Anti-de Sitter space of two dimensions in Jackiw-Teitelboim gravity, with its two boundaries coupled to each other and a dilaton that does not depend on the global time.
This note records the three charts, where each comes from, the field equations each is held to, and what each diagram draws.
The source writes $L = 1$; the metric file keeps $L$, and every coordinate of every chart is a number, so no factor of $c$ appears.
The three charts are written by `_tools/derivations/print_charts.py --metric eternal_traversable_wormhole`, and `verify_metrics.py --system eternal_traversable_wormhole/<chart>` checks each in a second.

The source was read in its arXiv copy, 1804.00491v3, the TeX source and the PDF, whose equation numbers are the ones cited here.
Gao, Jafferis and Wall, 1608.05687; Maldacena, Stanford and Yang, 1606.01857 and 1704.05333; Jensen, 1605.06098; Engelsöy, Mertens and Verlinde, 1606.03438; Maldacena and Stanford, 1604.07818; Maldacena, Milekhin and Popov, 1807.04726; and Maldacena, Michelson and Strominger, hep-th/9812073, were read for what the History says of them, and Sachdev and Ye, cond-mat/9212030, and Galloway, Schleich, Witt and Woolgar, gr-qc/9902061, for their abstracts and author lines.
Maldacena and Qi's paper has no journal on INSPIRE's record, and is cited as the arXiv preprint.

## Step 1. The field equations

The action is Maldacena and Qi's (2.5), $\frac{1}{2}\int\sqrt{-g}\,\phi(R + 2) + \dots + S_{\rm matter}$ with the topological term $\phi_0$ dropped, here with $2/L^2$ for $2$.
Varying $\phi$ gives $R = -2/L^2$, so every chart is locally anti-de Sitter space of two dimensions, with $K = R^2$ and no Einstein or Weyl tensor.
Varying the metric gives $\nabla_\mu\nabla_\nu\phi - g_{\mu\nu}\nabla^2\phi + g_{\mu\nu}\phi/L^2 = -T_{\mu\nu}$.
Its null component in the global chart is their (2.10), $-\partial_+(\sin^2\sigma\,\partial_+\phi) = T_{++}\sin^2\sigma$, since $\nabla_+\nabla_+\phi = \sin^{-2}\sigma\,\partial_+(\sin^2\sigma\,\partial_+\phi)$ for $ds^2 = -dx^+dx^-/\sin^2\sigma$; that fixes the sign of $T$.

## Step 2. The dilaton and its source

Appendix C of Maldacena and Qi works out a free fermion on the strip with the boundary conditions the coupling makes, their (C.177) to (C.181).
The part of the stress tensor that is not proportional to the metric is their (C.180), $\frac{1}{2}(T_{tt} + T_{\sigma\sigma}) = -\epsilon(1 - \epsilon)/4\pi$ for each fermion, with $T_{t\sigma} = 0$.
For $N$ fermions the dilaton is their (C.181), $\phi = N\frac{\epsilon(1 - \epsilon)}{4\pi}\left[(\pi/2 - \sigma)\cot\sigma + 1\right] + \frac{c}{24\pi}$.
The constant $c/24\pi$ comes from the part of the stress tensor proportional to the metric, and they note that a shift of the dilaton, with the topological term, removes it, so the metric file drops it and writes $\phi = \phi_h\left(1 + (\pi/2 - \sigma)\cot\sigma\right)$ with $\phi_h = N\epsilon(1 - \epsilon)/4\pi$.

Then $\nabla_\mu\nabla_\nu\phi - g_{\mu\nu}\nabla^2\phi + g_{\mu\nu}\phi/L^2 = \phi_h\,\mathrm{diag}(1, 1)$ in the global chart, exactly, which `eternal_traversable_wormhole_check` and the test file both compute.
So $T_{tt} = T_{\sigma\sigma} = -\phi_h$ and $T_{t\sigma} = 0$, traceless and constant in the global coordinates, which is their (C.180) times $N$ with $\phi_h = N\epsilon(1 - \epsilon)/4\pi$.
With $\sin^2\sigma\,\phi' = -\phi_h\left(\sin\sigma\cos\sigma + \pi/2 - \sigma\right)$, the derivative $(\sin^2\sigma\,\phi')' = 2\phi_h\sin^2\sigma$, so $T_{++} = -\phi_h/2$, their $-N\epsilon(1 - \epsilon)/8\pi$.

The dilaton is least at $\sigma = \pi/2$, where it is $\phi_h$, and near the boundary $\phi \to \pi\phi_h/2\sigma$.
With the boundary conditions (2.6), $ds = du/\epsilon$ and $\phi = \phi_r/\epsilon$ on the boundary, the boundary sits where $L\,dt/\sigma_b = du/\epsilon$, $\sigma_b = \epsilon\,dt/du$, which Maldacena and Qi write as $\sigma = \epsilon t'$ below their (4.36), and $\pi\phi_h/2\sigma_b = \phi_r/\epsilon$ gives $dt/du = \pi\phi_h/2\phi_r$ at $L = 1$.

In the embedding $-Y_{-1}^2 - Y_0^2 + Y_1^2 = -L^2$, with $Y_{-1} = L\cos t/\sin\sigma$, $Y_0 = L\sin t/\sin\sigma$ and $Y_1 = L\cot\sigma$ in the global chart, the dilaton is $\phi_h(1 + x\arctan x)$ with $x = Y_1/L$, since $\arctan(\cot\sigma) = \pi/2 - \sigma$ on $0 < \sigma < \pi$.
It depends on $Y_1$ alone, so the one isometry it keeps is the rotation of $Y_{-1}$ and $Y_0$, the shift of $t$.

## Step 3. The charts

`global`, $ds^2 = L^2(-dt^2 + d\sigma^2)/\sin^2\sigma$ on $0 < \sigma < \pi$, Maldacena and Qi's (2.1) with $L$ written out.

`poincare`, $ds^2 = L^2(-dt_P^2 + dz^2)/z^2$, the Poincaré chart of their (2.2), which prints $-dt_P + dz^2$ for $-dt_P^2 + dz^2$.
The map is $t_P \pm z = \tan((t \pm \sigma)/2)$, which at the boundary $\sigma = 0$ is their (2.4), $t_P = \tan(t_r/2)$.
Then $x = \cot\sigma = (1 + t_P^2 - z^2)/2z$, the boundary $\sigma = 0$ is $z = 0$, the middle $\sigma = \pi/2$ is $z^2 - t_P^2 = 1$, and the patch is $|t \pm \sigma| < \pi$.

`thermal`, $ds^2 = L^2(-\sinh^2\rho\,dt_R^2 + d\rho^2)$, the Rindler or thermal chart of their (2.2).
The map is $Y_{-1} = L\cosh\rho$, $Y_0 = L\sinh\rho\sinh t_R$ and $Y_1 = L\sinh\rho\cosh t_R$, so $\tanh t_R = \sin t/\cos\sigma$ and $\cosh\rho = \cos t/\sin\sigma$ on the wedge $\cos\sigma > |\sin t|$, and $x = \sinh\rho\cosh t_R$.
Its edge $\rho = 0$ is the event $t = 0$, $\sigma = \pi/2$ and, as $t_R \to \pm\infty$, the two light rays through it, which run from one boundary to the other.

`eternal_traversable_wormhole_check` in `print_charts.py` holds every chart to $R = -2/L^2$, $K = R^2$ and $G_{\mu\nu} = 0$, the global dilaton to $T_{\mu\nu} = -\phi_h\delta_{\mu\nu}$ exactly, every chart's source to a vanishing trace, and the Poincaré and thermal charts to being the global chart carried along their maps, metric, dilaton and stress tensor, at a point of rational coordinates to forty digits; sympy does not clear the thermal chart's trace, a sum of arctangents of $\sinh\rho\cosh t_R$, so that one is taken at two such points.

## Step 4. The diagrams

Every drawing is at $L = 1$; none depends on $\phi_h$.

The spacetime diagrams are the three planes of the charts, `null_rays.py` rows under `ETW`: the global strip $0 < \sigma < \pi$, $|t| < \pi/2$, with the middle $\sigma = \pi/2$ marked and the two rays through it at $t = 0$, which reach the boundaries at $t = \mp\pi/2$ and $\pm\pi/2$; the Poincaré plane out to $z = 3$, with the hyperbola $z^2 - t_P^2 = 1$ and the same two rays through $t_P = 0$, $z = 1$; and the thermal plane out to $\rho = 3$.
`--verify` holds every ray to $t \pm \sigma$, $t_P \pm z$ and $t_R \pm \ln\tanh(\rho/2)$.

The conformal diagram is the global chart itself, $T = t$ and $X = \sigma - \pi/2$, Maldacena and Qi's Figure 2(a), over $|T| < 1.1\pi$ with the strip faded above and below: the boundaries are the vertical sides, the middle is $X = 0$, and the rays through the middle at $t = 0$ run corner to corner of the square $|T| < \pi/2$.
The Poincaré patch is the triangle $|t \pm \sigma| < \pi$ beside $X = -\pi/2$, their Figure 1(c), carried in by $t \pm \sigma = 2\arctan(t_P \pm z)$, and the thermal wedge the triangle $|T| < -X$, their Figure 1(b), carried in through the embedding; every map is checked against its chart's published metric, and one event is checked to land on one point through all three.

The embedding diagram is the Euclidean section with the global time made periodic, the low temperature phase of their section 4.5.1, thermal anti-de Sitter space, whose circle of Euclidean time never closes, their Figure 7(a).
It is drawn at the period $\beta' = 2\pi$, with $t = -i\theta$: $L^2(d\sigma^2 + d\theta^2)/\sin^2\sigma$, a cylinder whose circle has the radius $L/\sin\sigma$, least, $L$, at the middle.
A surface of revolution in flat space needs $g_{\sigma\sigma} \ge (d\rho/d\sigma)^2$, $|\cot\sigma| \le 1$, so it is drawn from $\sigma = \pi/4$ to $3\pi/4$, where the radius is $\sqrt{2}\,L$ and the profile runs flat, and `stops` checks that no surface of revolution in flat space carries it beyond.
In $u = \cot\sigma$ the height is $z = \int_0^{-u}\sqrt{(1 - w^2)/(1 + w^2)}\,dw$, an elliptic integral, checked by quadrature.
Its meridians $\theta = 0$ and $\pi$ are the moment $t = 0$, which the other diagrams mark: across the strip, at $z = \tan(\sigma/2)$ on the Poincaré plane, and in the thermal chart from $\rho = 0$ out to $\cosh\rho = 1/\sin\sigma$ on the half $\sigma \le \pi/2$.

No curvature figure is drawn, since both scalars are the constants $R = -2/L^2$ and $K = 4/L^4$.
