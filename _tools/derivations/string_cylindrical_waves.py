"""Xanthopoulos's vacuum solution as printed by Gurtug and Sakalli, equations (1)-(3)."""

METRIC = 'string_cylindrical_waves'
SYSTEMS = ('prolate', 'cylindrical', 'hyperbolic')


def register_dimensions(vm):
    pass


def spec(system):
    raise NotImplementedError(system)


def register(pc):
    pc.CHARTS[METRIC] = [lambda system=system: spec(system) for system in SYSTEMS]
