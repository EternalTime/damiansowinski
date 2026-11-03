"""The four-parameter vacuum field in prolate and spheroidal coordinates.

The functions are those of Kurmanov et al., arXiv:2503.18028, equations (1)-(11).
The regular form of alpha extends the expression in that paper to zero rotation.
"""
import mpmath as mp

METRIC = 'quevedo_mashhoon'
SYSTEMS = ('prolate_spheroidal', 'spherical')


def potentials(x, y, mass, spin, deformation, delta):
    """Return f, omega, exp(2 gamma) on the real exterior chart."""
    sigma = mp.sqrt(mass * mass - spin * spin)
    alpha = -spin / (mass + sigma)
    p = (3 * y * y - 1) / 2
    Q0 = mp.log((x + 1) / (x - 1)) / 2
    Q1 = x * Q0 - 1
    Q2 = (3 * x * x - 1) * Q0 / 2 - 3 * x / 2
    Q2p = 3 * x * Q0 - (3 * x * x - 1) / (2 * (x * x - 1)) - mp.mpf(3) / 2
    def d(sign):
        return (mp.log((x + sign * y) ** 2 / (x * x - 1)) / 2
                + mp.mpf(3) / 2 * (1 - y * y - sign * x * y)
                + mp.mpf(3) / 4 * (x * (1 - y * y) - sign * y * (x * x - 1))
                * mp.log((x - 1) / (x + 1)))
    lp = alpha * (x*x-1)**(1-delta) * (x+y)**(2*delta-2) * mp.exp(2*deformation*delta*d(1))
    lm = alpha * (x*x-1)**(1-delta) * (x-y)**(2*delta-2) * mp.exp(2*deformation*delta*d(-1))
    def ab(sign):
        power = (x + sign)**(delta - 1)
        return (power * (x*(1-lp*lm) + sign*(1+lp*lm)),
                power * (y*(lp+lm) - sign*(lp-lm)))
    ap, bp = ab(1)
    am, bm = ab(-1)
    R, L = ap*am + bp*bm, ap*ap + bp*bp
    H = (x+1)**(delta-1) * (x*(1-y*y)*(lp+lm)*ap + y*(x*x-1)*(1-lp*lm)*bp)
    chi = ((1+deformation)**2 * mp.log((x*x-1)/(x*x-y*y))/2
           + 2*deformation*(1-p)*Q1
           + deformation**2*(1-p)*((1+p)*(Q1*Q1-Q2*Q2)
               + (x*x-1)*(2*Q2*Q2-3*x*Q1*Q2+3*Q0*Q2-Q2p)/2))
    e = mp.exp(2*deformation*delta*p*Q2)
    return R/L/e, -2*spin - 2*sigma*H/R*e, ((1+mass/sigma)**2/4 * R/(x*x-1)**delta
                                                    * mp.exp(2*delta*delta*chi))
