# Bowen-York initial data

## A maximal spatial slice

The initial data are a Riemannian metric $\gamma_{ij}$ and an extrinsic curvature $K_{ij}$ on a spatial hypersurface, with $K = \gamma^{ij}K_{ij} = 0$.
In geometrized units $G = c = 1$, the conformal decomposition is $\gamma_{ij} = \psi^4\delta_{ij}$ and $K_{ij} = \psi^{-2}\widehat A_{ij}$, with $\psi > 0$.
Brandt and Brügmann, Phys. Rev. Lett. 78, 3606 (1997), equations (1)-(3), distinguish these physical and conformal tensors.
All curvature of $\gamma_{ij}$ is intrinsic curvature in three dimensions.

## The momentum constraint

For $r = |\boldsymbol{x}|$ and $n^i = x^i/r$, the conformal tensor is

$$\widehat A^{ij} = \frac{3}{2r^2}\left(P^i n^j + P^j n^i - (\delta^{ij} - n^i n^j)P^k n_k\right) + \frac{3}{r^3}\left(\epsilon^{ikl}S_k n_l n^j + \epsilon^{jkl}S_k n_l n^i\right).$$

Here $P^i$ and $S^i$ are constant linear and angular momentum parameters, and all indices on $\widehat A^{ij}$ are raised and lowered with the Euclidean metric $\delta_{ij}$.
Brandt and Brügmann's equation (5) states this tensor, and equations (2) and (7) state its divergence constraint and its superposition over punctures.
The Hamiltonian equation is

$$\nabla^2\psi = -\frac{1}{8}\psi^{-7}\widehat A_{ij}\widehat A^{ij},$$

where $\nabla^2$ is the Euclidean Laplacian.

## Punctures

Brandt and Brügmann's equations (8)-(10) set $\psi = \sum_A m_A/(2r_A) + u$, with $u\to1$ at spatial infinity and $u$ of class $C^2$ at each puncture.
Here $r_A = |\boldsymbol{x} - \boldsymbol{x}_A|$, $\boldsymbol{x}_A$ is the puncture position, and $m_A>0$ is its bare mass parameter.
Each puncture is another asymptotically flat end, as their equations (11)-(14) show.
Their individual momentum and spin parameters agree with individual ADM quantities in the limit of large separation, as the text following equation (7) specifies.
