"""The binary icosahedral group I* acting on the unit 3-sphere, for the Poincare dodecahedral universe.

A point of the 3-sphere is the unit quaternion w + x i + y j + z ij, which Aurich, Lustig and
Steiner's (7) writes as (cos chi, sin chi sin theta cos phi, sin chi sin theta sin phi,
sin chi cos theta). The group is the closure of their two generators, gamma_1 = j and
gamma_2 = sigma/2 + i/(2 sigma) + j/2 with sigma the golden ratio, acting by left multiplication,
their (18). `turned()` conjugates it by a rotation of the 3-sphere that fixes the point 1, so that
the nearest image of 1 along +i, the centre of a face of the cell about 1, is cos(pi/5) + sin(pi/5) i:
the face then lies in the direction theta = pi/2, phi = 0 of (7). Pure Python, so the tests can use it.
"""
import math

SIGMA = (1 + math.sqrt(5)) / 2
ONE = (1.0, 0.0, 0.0, 0.0)


def mul(p, q):
    """The quaternion product in the basis (1, i, j, ij), with i^2 = j^2 = -1 and ij = -ji."""
    a1, b1, c1, d1 = p
    a2, b2, c2, d2 = q
    return (a1 * a2 - b1 * b2 - c1 * c2 - d1 * d2, a1 * b2 + b1 * a2 + c1 * d2 - d1 * c2,
            a1 * c2 - b1 * d2 + c1 * a2 + d1 * b2, a1 * d2 + b1 * c2 - c1 * b2 + d1 * a2)


def conj(q):
    return (q[0], -q[1], -q[2], -q[3])


def key(q):
    return tuple(round(x, 9) + 0.0 for x in q)


def distance(p, q):
    """The distance on the unit 3-sphere, cos d = p.q."""
    return math.acos(max(-1.0, min(1.0, sum(a * b for a, b in zip(p, q)))))


def group():
    """The 120 elements, the closure of Aurich, Lustig and Steiner's generators, 1 first."""
    generators = ((0.0, 0.0, 1.0, 0.0), (SIGMA / 2, 1 / (2 * SIGMA), 0.5, 0.0))
    found = {key(ONE): ONE}
    frontier = [ONE]
    while frontier:
        fresh = []
        for q in frontier:
            for g in generators:
                r = mul(g, q)
                if key(r) not in found:
                    found[key(r)] = r
                    fresh.append(r)
        frontier = fresh
    return list(found.values())


def neighbours(elements):
    """The elements that carry 1 the shortest distance, pi/5: one for each face of the cell."""
    return [g for g in elements if abs(distance(ONE, g) - math.pi / 5) < 1e-9]


def turning(g, target=(1.0, 0.0, 0.0)):
    """The unit quaternion u with u g u^-1 = cos(pi/5) + sin(pi/5) target, for g of real part cos(pi/5):
    the rotation about the axis g's axis cross target, through the angle between them."""
    axis = [x / math.sin(math.pi / 5) for x in g[1:]]
    cross = (axis[1] * target[2] - axis[2] * target[1], axis[2] * target[0] - axis[0] * target[2],
             axis[0] * target[1] - axis[1] * target[0])
    half = math.acos(max(-1.0, min(1.0, sum(a * b for a, b in zip(axis, target))))) / 2
    size = math.sqrt(sum(x * x for x in cross))
    if size < 1e-12:
        return ONE if sum(a * b for a, b in zip(axis, target)) > 0 else (0.0, 0.0, 0.0, 1.0)
    return (math.cos(half), *(math.sin(half) * x / size for x in cross))


def turned():
    """The group conjugated so that one face of the cell about 1 has its centre along +i."""
    elements = group()
    u = turning(neighbours(elements)[0])
    return [mul(mul(u, g), conj(u)) for g in elements]
