"""The static magnetic exterior in Eiroa and Simeone (36)-(39)."""
import sympy as sp

METRIC = 'superconducting_string'
SYSTEM = 'cylindrical'
DIMENSIONS = {(METRIC, SYSTEM): {'t':'L', 'r':'L', r'\phi':'1', 'z':'L',
    'm':'1', 'k':'1', 'W':'1', r'\ell':'L', 'F':'1', 'H':'1', 'P':'1', 'r_c':'L'}}
HELD = {(METRIC, SYSTEM): ('F', 'H', 'P')}
RATES = {(METRIC, SYSTEM): {'F': {'r':'mH/r'}, 'H': {'r':'mF/r'}, 'P': {'r':'2m^2P/r'}}}


def chart_spec():
    """The normalized cylindrical exterior."""
    raise NotImplementedError
