"""Drawings of the cylindrical wave geometry at p=1, alpha=2 and ell=1."""

import numpy as np

METRIC = 'string_cylindrical_waves'
PARAMS = {'p': '1', 'alpha': '2', 'ell': '1'}


def register_null(nr):
    charts = [
        ('prolate', r'\tau', r'\sigma', (.8, 4.2, -2., 2.), 'tau'),
        ('hyperbolic', r'\psi', r'\theta', (-.15, 2.85, -1.5, 1.5), 'psi'),
    ]
    for system, time, radius, box, time_name in charts:
        nr.DIAGRAMS.append(nr.Diagram(METRIC, system, 'radial', f'${time}$ and ${radius}$',
                                     (time, radius), box, f'${radius}$', f'${time}$',
                                     PARAMS, {'phi': '0', 'z': '0'}, tau=time_name))
        nr.CAPTIONS[(METRIC, system, 'radial')] = [
            f'The plane of ${time}$ and ${radius}$ ($\\phi=0$, $z=0$), with radial null geodesics.']
    nr.CLOSED_FORMS[(METRIC, 'prolate', 'radial')] = (
        lambda t, r: np.arcsinh(t)+np.arccosh(np.maximum(r, 1.)),
        lambda t, r: np.arcsinh(t)-np.arccosh(np.maximum(r, 1.)), lambda t, r: r > 1.001)
    nr.CLOSED_FORMS[(METRIC, 'hyperbolic', 'radial')] = (
        lambda t, r: t+r, lambda t, r: t-r, lambda t, r: r > .001)
