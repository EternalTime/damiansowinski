# The maximally supersymmetric plane wave

The plane wave of type IIB supergravity in ten dimensions that keeps all 32 supersymmetries, held up by a null 5-form field.
The three charts are written by `_tools/derivations/print_charts.py --metric maximally_supersymmetric_plane_wave`, which took about 3 seconds for all three on 2 October 2026, and `verify_metrics.py --system maximally_supersymmetric_plane_wave/<chart>` checks each in seconds.

## Step 1. The papers and the source of each chart

Every paper of the candidate entry was verified on Crossref and arXiv before anything was built, and read in its TeX source on arXiv: Blau, Figueroa-O'Farrill, Hull and Papadopoulos (hep-th/0110242 and hep-th/0201081), Berenstein, Maldacena and Nastase (hep-th/0202021), Metsaev (hep-th/0112044), Berenstein and Nastase (hep-th/0205048), Marolf and Ross (hep-th/0208197), and Hubeny and Rangamani (hep-th/0210234).
Cahen and Wallach (1970) was read in full from the Bulletin of the American Mathematical Society, which is open.
Kowalski-Glikman (1984) and Penrose (1976) were verified on Crossref and INSPIRE but have no open copy, and the History says of them only what their titles and Blau and his coauthors say.
Names are from the author lines; Ruslan R. Metsaev's first name is INSPIRE's, Nolan Wallach's and Michel Cahen's Crossref's records of their other papers.

- `brinkmann`: Blau et al.'s solution of hep-th/0110242 with $\lambda = \mu/2$, which is Metsaev's line element and Marolf and Ross's, $ds^2 = -2\,dx^+dx^- - \mu^2x^2(dx^+)^2 + dx^2$, with $x^+ = cu$ and $x^- = v$.
- `rosen`: the Penrose limit of anti-de Sitter space times a 5-sphere in Rosen coordinates, hep-th/0201081 at $\rho = 1$, $R^{-2}\bar g = du\,dv + \sin^2(u/2)\,ds^2(\mathbb{E}^8)$, rescaled to $-2c\,du\,dv + \sin^2(\mu cu)\,dy^2$; the map onto the Brinkmann chart is $x_i = y_i\sin(\mu cu)$, $v_B = v + \frac{\mu}{2}y^2\sin(\mu cu)\cos(\mu cu)$, their change to Brinkmann coordinates.
- `conformally_flat`: Berenstein and Nastase's map onto Minkowski space, $\mu cU = \tan(\mu cu)$, $x_i = X_i\cos(\mu cu)$, $v = V - \mu^2cU X^2/2(1 + \mu^2c^2U^2)$, so that $ds^2 = (-2c\,dU\,dV + dX^2)/(1 + \mu^2c^2U^2)$ on $|\mu cu| < \pi/2$.

The five-form is Metsaev's, $F_{u1234} = F_{u5678} = 2\mu$ in the chart $x^0 = cu$, in his normalization $R_{MN} = F_{MPQRS}F_N{}^{PQRS}/24$, which puts the one component $R_{uu} = 8\mu^2$ of the Ricci tensor on the square of $F$.

## Step 2. What the script checks

`msw_check` holds the Brinkmann chart to $R_{MN} = F_{MPQRS}F_N{}^{PQRS}/24$ in every slot, $F$ to being equal to its Hodge dual up to the sign the orientation sets, and the Weyl tensor to vanishing, which is Blau et al.'s conformal flatness for a matrix $A_{ij}$ proportional to the identity.
The Rosen and conformally flat charts are the Brinkmann chart pulled back, $J^TgJ$ slot by slot through the maps above.
The Ricci scalar and the Kretschmann scalar vanish in every chart.

## Step 3. The Einstein static universe

Berenstein and Nastase's chain of maps carries the wave on to the Einstein static universe, $\mathbb{R}\times S^9$, and its last line reads
$$ds^2 = \frac{1}{4|e^{i\psi} - \cos\alpha\,e^{i\beta}|^2}\left(-d\psi^2 + d\alpha^2 + \cos^2\alpha\,d\beta^2 + \sin^2\alpha\,d\Omega_7^2\right),$$
which Marolf and Ross copy.
The factor $1/4$ is a slip: the line above it, $(1 + \tilde u^2)(1 + \tilde v^2)/4(1 + u^2)$ times the metric of the Einstein static universe, is $1/|e^{i\psi} - \cos\alpha\,e^{i\beta}|^2$ times it, and at $\psi = 0$, $\alpha = \pi/2$, which is the event $u = 0$ of their Minkowski space, where its conformal factor is $1$, the last line gives $1/4$.
Their chain of maps was followed numerically from the Einstein static universe back to the Brinkmann chart, and the pulled back metric is $1/|e^{i\psi} - \cos\alpha\,e^{i\beta}|^2$ times the Einstein static universe's to 40 digits at every point tried; `_tools/test_maximally_supersymmetric_plane_wave.py` holds that.
Written as a chart of ten coordinates, its curvature printed to 17 MB, since the factor's $\cos(\psi - \beta)$ is expanded through every component, so it is not a published chart and the conformal diagram draws it.

On the axis $x = 0$ of the Brinkmann chart, $\alpha = 0$, the map is
$$\mu cu = \frac{\psi + \beta - \pi}{2}, \qquad \mu v = -\frac{1}{2}\cot\frac{\psi - \beta}{2},$$
so the axis plane fills the cylinder of $\psi$ and $\beta$ once, but for the null line $\psi - \beta \in 2\pi\mathbb{Z}$, the conformal boundary, which is $v = \pm\infty$.
The Rosen chart's axis covers the band $\pi < \psi + \beta < 3\pi$ and the conformally flat chart's the band $0 < \psi + \beta < 2\pi$.

## Step 4. The drawings

Everything is drawn at $\mu = 1$, in units of $1/\mu$.
The spacetime diagrams are figures of light rays on the surface where every direction across the wave but the first is zero, seen from the side with the null coordinate left out, in each of the three charts (`msw_rays` in `projections.py`).
Each ray is traced with the published Christoffel symbols and checked null and against its closed form: $x_1 = s\sin(\mu cu)$ and $a\cos(\mu cu)$ in the Brinkmann chart, Marolf and Ross's null geodesics; $y_1 = s$ and $a\cot(\mu cu)$ in the Rosen chart; $X_1 = s\,cU$ and $a$ in the conformally flat chart.
The conformal diagram is Marolf and Ross's first figure, the cylinder of $\psi$ and $\beta$ at $\alpha = 0$ with its edges one line, with each chart's axis plane drawn on it by the map of Step 3.
The embedding diagram is a ring of free particles at rest at $u = 0$ in the plane of $x_1$ and $x_2$, which falls through the axis together at $\mu cu = \pi/2$ and comes back turned through half a circle at $\mu cu = \pi$, as the pp-wave's ring does; each wave front of constant $u$ is flat.
