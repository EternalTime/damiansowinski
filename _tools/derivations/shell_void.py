"""Accreting dust shell in an Einstein-de Sitter exterior, G=c=1.

Sakai-Haines (2-8)-(2-10), (2-15), (2-18), (2-20), (2-24)-(2-28).
The angular junction fixes density; the temporal and normal junctions fix motion.
The Newtonian late-time power law is not used as an evolution equation.
"""
from functools import lru_cache

import numpy as np
from scipy.integrate import solve_ivp

T0, R0, V0, END = 1.0, 0.15, 0.0, 4.0


def scale(t):
    return np.asarray(t) ** (2 / 3)


def quantities(t, radius, speed):
    hubble = 2 / (3 * t)
    rho = 1 / (6 * np.pi * t * t)
    gamma = 1 / np.sqrt(1 - speed * speed)
    x = hubble * radius
    clock = np.sqrt(1 + 2 * x * speed + x * x)
    sigma = hubble * hubble * radius / (4 * np.pi * gamma * (clock + 1 + x * speed))
    return hubble, rho, gamma, clock, sigma


def rhs(t, y):
    radius, speed = y[:2]
    hubble, rho, gamma, clock, sigma = quantities(t, radius, speed)
    acceleration = (-gamma * speed * hubble + 2 * np.pi * sigma
                    - rho * gamma * gamma * speed * speed / sigma) / gamma ** 3
    return np.array([hubble * radius + speed, acceleration, clock])


@lru_cache(maxsize=1)
def trajectory():
    solution = solve_ivp(rhs, (T0, END), [R0, V0, 0.0], method="DOP853",
                         rtol=2e-12, atol=2e-14, dense_output=True, max_step=0.025)
    if not solution.success:
        raise RuntimeError(solution.message)
    return solution


def state(t):
    if np.any(np.asarray(t) < T0) or np.any(np.asarray(t) > END):
        raise ValueError("shell trajectory is defined only for 1 <= t <= 4")
    return trajectory().sol(t)


def radius_at_inner_time(inner):
    from scipy.optimize import brentq
    end = state(END)[2]
    values = np.asarray(inner)
    def one(value):
        if value < 0 or value > end:
            return np.nan
        t = brentq(lambda t: state(t)[2] - value, T0, END, xtol=1e-13)
        return state(t)[0]
    return np.vectorize(one, otypes=[float])(values)
