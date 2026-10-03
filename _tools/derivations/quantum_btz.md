# The quantum BTZ black hole

The quantum BTZ black hole is the metric induced on a brane at $x = 0$ of the anti-de Sitter C-metric, read by braneworld holography as the black hole of three dimensions with the exact backreaction of a strongly coupled conformal field theory.
Its five charts are written by `_tools/derivations/print_charts.py --metric quantum_btz`, and `verify_metrics.py --system quantum_btz/<chart>` checks each in seconds.
Equation numbers are those of Emparan, Frassino and Way, arXiv:2007.15999v4, the text of JHEP 11 (2020) 137, counted from its source.

## Step 1. The charts and their sources

- `static`, their (2.37): $ds^2 = -H\,c^2dt^2 + dr^2/H + r^2d\phi^2$ with $H = r^2/\ell_3^2 - M - \ell F/r$ and $\phi$ of period $2\pi$.
  Their $\bar t$, $\bar r$ and $\bar\phi$ are written $t$, $r$ and $\phi$, and their $8\mathcal{G}_3M$ is written $M$, the reading the BTZ entry gives its $M = 8Gm/c^2$.
  Emparan, Fabbri and Kaloper (hep-th/0206155) first wrote this metric, as their (4.3), with $r_1(M)$ for $\ell F$.
- `eddington_finkelstein_ingoing` and `eddington_finkelstein_outgoing`, the advanced and retarded times $v, u = ct \pm r_*$ built on the tortoise coordinate $dr_*/dr = 1/H$, with $r_*(0) = 0$, which the conformal diagram's tower takes too.
- `brane`, their (2.35): the metric (2.1) induces at $x = 0$, $H = r^2/\ell_3^2 + \kappa - \mu\ell/r$ from their (2.2), with $\phi$ of period $2\pi\Delta$, $\Delta = 2x_1/(3 - \kappa x_1^2)$, their (2.34).
- `rotating`, their (3.10): the metric their rotating C-metric (3.1) induces at $x = 0$, with $H$ from their (3.2), $r^2/\ell_3^2 + \kappa - \mu\ell/r + a^2/r^2$, and the shift $a/r^2$.
  Its points are identified along $\phi \to \phi + 2\pi\Delta$ at fixed $ct + \tilde a\ell_3\phi$, the orbits of their Killing vector (3.8), with $\Delta$ from (3.9) and $\tilde a$ from (3.7), and the sentence after (3.9); the chart's convention states it.
  The canonical form (3.15), in $\bar t$, $\bar r$ and $\bar\phi$ with $r$ a function of $\bar r$, is this chart under the linear map (3.11) and the change of radius (3.13); its tensors hold the radical of (3.13) in every value, and it is not printed.

## Step 2. What every chart is held to

`quantum_btz_check` in `print_charts.py` holds each chart to the following before it is written.

The Ricci scalar is $-6/\ell_3^2$ in every chart, so the stress tensor of the quantum fields is traceless.
With $G^a{}_b - \delta^a{}_b/\ell_3^2 = 8\pi G_3\langle T^a{}_b\rangle_0$, their (2.44), the static chart gives $(\ell F/2r^3)\,\mathrm{diag}(1, 1, -2)$, their (2.46); the brane chart the same with $\mu$ for $F$; and the rotating chart their (3.19), the same three with $\mu$ and $\langle T^\phi{}_t\rangle_0 = 3\ell\mu a/(16\pi G_3 r^5)$.
Each Eddington-Finkelstein chart is the static one pulled back along $c\,dt = dv - dr/H$ or $du + dr/H$.
The brane chart is the static one under their rescaling (2.36), $t = \Delta\bar t$, $r = \bar r/\Delta$, $\phi = \Delta\bar\phi$, with $M = -\kappa\Delta^2$ and $F = \mu\Delta^3$, their (2.38) with (2.40) and (2.41).
The rotating chart at $a = 0$ is the brane chart.

The Kretschmann scalar is $12/\ell_3^4 + 6\ell^2F^2/r^6$ in the static charts and $12/\ell_3^4 + 6\mu^2\ell^2/r^6$ in the brane and rotating charts, where $a$ drops out of it.

## Step 3. The mass and the strength of the stress tensor

With $k = -\kappa x_1^2$, their (2.38) and (2.41) are $M = 4k/(3 + k)^2$ and $F = 8(1 + k)/(3 + k)^3$, which the parameter descriptions state.
$k = -1$ is global anti-de Sitter space of three dimensions, $M = -1$ and $F = 0$; $-1 < k < 0$ is branch $1a$, the conical singularities dressed with horizons; $0 < k < 3$ is branch $1b$ and $k > 3$ branch 2, which meet at the largest mass $M = 1/3$, their (2.50) in these units.
`_tools/test_quantum_btz.py` holds those, the temperature (2.68) against the published brane chart, and the first law $dM = T\,dS_{\rm gen}$ of (2.82) from (2.63), (2.68) and (2.70).

## Step 4. What is drawn

- The static hole at $k = 1$, $M = F = 1/4$, on branch $1b$, with $\ell = 15\ell_3/16$, so that $H = (r - 3/4)(r^2 + 3r/4 + 5/16)/r$ in units of $\ell_3$, one horizon at $r_+ = 3\ell_3/4$ against the classical hole's $\ell_3/2$, and a complex pair, which is the shape of Schwarzschild-anti-de Sitter's $f$.
  Its surface gravity is $23/(24\ell_3)$ and its tortoise coordinate tends to $R = 0.679\,\ell_3$, so its conformal diagram is `AdSHoleTower`'s with $B = e^{2\kappa R} = 3.67$.
- The dressed cone in the brane chart at $\kappa = +1$, $\mu = 6$, $\ell = \ell_3/3$: $x_1 = 1/2$, $\Delta = 4/11$ and $M = -16/121$; $H = r^2 + 1 - 2/r$, which is Schwarzschild-anti-de Sitter's $f$ at $r_s = 2$ and $L = 1$, with its horizon at $\ell_3$.
- The rotating hole at $\kappa = -1$, $\mu = 19/8$, $\ell = 3\ell_3/19$, $a^2 = 3\ell_3^2/8$, with $x_1 = 1$ and $\tilde a = a$: $r^2H = (r - 1)(r - 1/2)(r^2 + 3r/2 + 3/4)$, horizons at $\ell_3$ and $\ell_3/2$, the tower of Reissner and Nordström's hole in anti-de Sitter space, which `CarterAdSTower` and `HyperbolicHoleDrawing` draw.
  The ergosurface, $g_{tt} = 0$, is the root $1.151$ of $r^3 - r - 3/8$.
- The embedding diagram is the moment $t = 0$ of the static hole, in flat space out to $H = 1$ at the root $1.202$ of $64r^3 - 80r - 15$ and in Minkowski space beyond, as the BTZ hole's is.
  The dressed cone's circles have circumference $2\pi\Delta r$, and the rotating hole's moment is drawn at no other example, so neither has a surface of its own.
