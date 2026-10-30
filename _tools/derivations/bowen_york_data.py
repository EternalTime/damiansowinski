"""The axisymmetric single-spin throat, in units of its coordinate radius a.

With t=a/r and mu=cos(theta), psi=1+t+u satisfies
t^2 psi_tt + d_mu((1-mu^2) psi_mu) = -(9/4) j^2 t^4 (1-mu^2) psi^-7,
where j=S/a^2. At t=0 psi=1; at t=1 psi_t=psi/2.
Even Legendre modes make the axes regular and impose equatorial reflection.
Chebyshev collocation includes spatial infinity without an outer cutoff.
"""
import functools

import numpy as np
from numpy.polynomial import chebyshev as cheb, legendre as leg
from scipy.linalg import lu_factor, lu_solve


class SpinThroat:
    """Solve the full angular Hamiltonian equation by fixed-point iteration."""

    def __init__(self, spin=1.0, nodes=32, modes=8):
        self.spin, self.nodes, self.modes = float(spin), nodes, modes
        k = np.arange(nodes + 1)
        x = np.cos(np.pi * k / nodes)
        self.t = (x + 1) / 2
        c = np.where((k == 0) | (k == nodes), 2.0, 1.0) * (-1.0)**k
        delta = x[:, None] - x[None, :]
        D = np.outer(c, 1/c) / (delta + np.eye(nodes + 1))
        D -= np.diag(D.sum(axis=1))
        D *= 2
        D2 = D @ D
        mu, w = leg.leggauss(4 * modes + 8)
        P = np.array([leg.legval(mu, [0]*(2*l) + [1]) for l in range(modes)])
        projection = (4*np.arange(modes)[:, None] + 1) * P * w / 2
        factors = []
        for l in 2*np.arange(modes):
            A = self.t[:, None]**2 * D2 - l*(l+1)*np.eye(nodes + 1)
            A[0] = D[0]
            A[0, 0] -= 0.5
            A[-1] = 0
            A[-1, -1] = 1
            factors.append(lu_factor(A))
        u = np.zeros((modes, nodes + 1))
        for iteration in range(100):
            psi = 1 + self.t[:, None] + u.T @ P
            source = -2.25 * spin**2 * self.t[:, None]**4 * (1-mu**2) / psi**7
            b = projection @ source.T
            b[:, [0, -1]] = 0
            new = np.array([lu_solve(factor, rhs) for factor, rhs in zip(factors, b)])
            miss = np.max(np.abs(new-u))
            u = new
            if miss < 2e-14:
                break
        else:
            raise RuntimeError("the spinning-throat Hamiltonian iteration did not converge")
        self.iterations = iteration + 1
        self.coeff = np.array([cheb.chebfit(x, row, nodes) for row in u])
        # Coefficients below roundoff contribute only differentiation noise.
        self.coeff[np.abs(self.coeff) < 2e-15] = 0

    def field(self, r, mu, dr=0, dm=0):
        """psi and its first two derivatives in r and mu, outside the throat."""
        r, mu = np.broadcast_arrays(np.asarray(r, dtype=float), np.asarray(mu, dtype=float))
        t = 1/r
        out = np.zeros_like(r)
        for l, c in enumerate(self.coeff):
            angular = leg.legder([0]*(2*l)+[1], dm)
            if dr == 0:
                radial = cheb.chebval(2*t-1, c)
            elif dr == 1:
                radial = -2*t*t * cheb.chebval(2*t-1, cheb.chebder(c))
            elif dr == 2:
                radial = (4*t**3 * cheb.chebval(2*t-1, cheb.chebder(c))
                          + 4*t**4 * cheb.chebval(2*t-1, cheb.chebder(c, 2)))
            else:
                raise ValueError(dr)
            out += radial * leg.legval(mu, angular)
        if dm == 0:
            out += (1+t, -t*t, 2*t**3)[dr]
        return out

    def residual(self, r, mu):
        p = self.field(r, mu)
        lap = (self.field(r, mu, 2) + 2*self.field(r, mu, 1)/r
               + ((1-mu*mu)*self.field(r, mu, dm=2) - 2*mu*self.field(r, mu, dm=1))/r**2)
        source = 2.25*self.spin**2*(1-mu*mu)/r**6/p**7
        return lap + source

    def equator(self, r):
        return self.field(r, 0), self.field(r, 0, 1)

    def formula(self):
        """A finite polynomial in 1/r and cos(theta), using the existing function schema.

        Horner forms keep the numerical declaration short and its evaluation stable.
        The polynomial is the numerical approximation, not an exact vacuum solution.
        """
        def horner(c, variable):
            keep = np.flatnonzero(np.abs(c) > 2e-15)
            if not len(keep):
                return "0"
            c = c[:keep[-1]+1]
            out = f"({c[-1]:.17g})"
            for v in c[-2::-1]:
                out = f"(({v:.17g})+({variable})*{out})"
            return out
        terms = ["1+1/r"]
        for l, c in enumerate(self.coeff):
            keep = np.flatnonzero(np.abs(c) > 2e-15)
            if not len(keep):
                continue
            radial = cheb.cheb2poly(c[:keep[-1]+1])
            angular = leg.leg2poly([0]*(2*l)+[1])
            terms.append(horner(radial, "2/r-1") + "*" + horner(angular, "cos(theta)"))
        return "+".join(terms)


@functools.lru_cache(maxsize=None)
def solution(spin=1.0, nodes=32, modes=8):
    return SpinThroat(spin, nodes, modes)
