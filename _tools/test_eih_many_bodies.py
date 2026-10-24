"""Einstein, Infeld and Hoffmann's field of many bodies to post-Newtonian order:
.venv.noindex/bin/python -m unittest discover -s _tools

What its texts state, held on the published files: the potentials of point masses on any paths
solve the field equations the parameters state; for two bodies the harmonic chart is Blanchet, Faye
and Ponsot's metric (7.2) and, on a circular orbit, Alvi's in corotating coordinates; one body at
rest is the published post-Newtonian metric at general relativity's values; the equations of motion
(17.2) of the paper of 1938 are the two-body acceleration of Will's review, close a circular orbit
at the angular velocity the drawings declare, and turn an eccentric one by 6 pi m/p a revolution;
and the numbers the History and the captions give. It also holds the rule by which the checker
counts a potential and its time derivatives as post-Newtonian orders. The tests need sympy, numpy
and scipy, and are skipped where they are absent.
"""
import importlib.util
import json
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "MFS" / "assets" / "data"
HAS_SYMPY = all(importlib.util.find_spec(name) is not None for name in ("sympy", "numpy", "scipy"))
sys.path.insert(0, str(Path(__file__).resolve().parent / "derivations"))

POTENTIALS = ("U", "psi", "chi", "V_1", "V_2", "V_3")


def flat(tensor):
    if isinstance(tensor, list):
        for item in tensor:
            yield from flat(item)
    else:
        yield tensor


def bodies(sp, t, places, velocities, accelerations, masses):
    """The potentials of point masses on the paths x_a + v_a t + a_a t^2/2, as the parameters of the
    charts define them, with c = 1 and each mass a length: U, psi, chi and the three components of V
    as expressions in t and the symbols x, y and z."""
    x, y, z = sp.symbols("x y z", real=True)
    paths = [[p + v * t + a * t ** 2 / 2 for p, v, a in zip(*body)] for body in zip(places, velocities, accelerations)]
    speeds = [[sp.diff(c, t) for c in path] for path in paths]
    distance = [sp.sqrt(sum((q - c) ** 2 for q, c in zip((x, y, z), path))) for path in paths]
    U = sum(m / r for m, r in zip(masses, distance))
    chi = sum(m * r for m, r in zip(masses, distance))
    V = [sum(m * v[i] / r for m, v, r in zip(masses, speeds, distance)) for i in range(3)]
    psi = 0
    for a, (m, r) in enumerate(zip(masses, distance)):
        others = sum(masses[b] / sp.sqrt(sum((p - q) ** 2 for p, q in zip(paths[a], paths[b])))
                     for b in range(len(masses)) if b != a)
        psi += m / r * (sp.Rational(3, 2) * sum(v ** 2 for v in speeds[a]) - others)
    return (x, y, z), {"U": U, "psi": psi, "chi": chi, "V_1": V[0], "V_2": V[1], "V_3": V[2]}


@unittest.skipUnless(HAS_SYMPY, "sympy, numpy and scipy are not installed")
class Published(unittest.TestCase):
    """The published charts, read through the checker's Reader, with c = 1."""

    @classmethod
    def setUpClass(cls):
        import sympy
        import verify_metrics
        cls.sp, cls.vm = sympy, verify_metrics
        cls.metric = json.loads((DATA / "metrics" / "eih_many_bodies.json").read_text(encoding="utf-8"))

    def components(self, system, field="metric_components"):
        entry = next(c for c in self.metric["coordinates"] if c["id"] == system)
        reader = self.vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
        return reader, {tuple(e["indices"]): self.sp.sympify(reader(e["value"])) for e in entry[field]}

    def with_bodies(self, system, potentials, space):
        """The published metric of a chart with the potentials of declared bodies put in, as
        expressions in the chart's t and the symbols of `space`."""
        sp = self.sp
        reader, g = self.components(system)
        t = reader.symbol["t"]
        at = dict(zip([reader.symbol[k] for k in "xyz"], space))
        out = {}
        for index, value in g.items():
            for name in POTENTIALS:
                function = reader.parameters[name]
                value = value.replace(function.func, sp.Lambda(tuple(function.args), potentials[name].subs(
                    dict(zip(space, function.args[1:])))))
            out[index] = value.doit().subs(reader.c, 1).subs(at)
        return t, out

    def two_bodies(self):
        """Two bodies at rational places with rational velocities, each falling toward the other as
        Newton has it, and an event off both."""
        sp = self.sp
        R = sp.Rational
        m = (R(3, 100), R(1, 50))
        places = ((R(2), R(1, 3), R(-1, 2)), (R(-1), R(3, 4), R(1)))
        velocities = ((R(1, 10), R(-1, 20), R(1, 30)), (R(-1, 15), R(1, 25), R(-1, 40)))
        apart = [p - q for p, q in zip(*places)]
        d = sp.sqrt(sum(c ** 2 for c in apart))
        accelerations = (tuple(-m[1] * c / d ** 3 for c in apart), tuple(m[0] * c / d ** 3 for c in apart))
        return m, places, velocities, accelerations, d, (R(1, 2), R(-3, 2), R(2))

    def test_the_potentials_of_point_masses_solve_the_field_equations_the_parameters_state(self):
        """On any paths: U, psi and each component of V solve Laplace's equation off the bodies, the
        Laplacian of chi is 2U, and d_t U + div V = 0, which is what makes every Ricci component
        kept vanish between the bodies."""
        sp = self.sp
        t = sp.Symbol("t", real=True)
        m, places, velocities, accelerations, _, event = self.two_bodies()
        (x, y, z), f = bodies(sp, t, places, velocities, accelerations, m)
        at = {x: event[0], y: event[1], z: event[2], t: sp.Rational(1, 7)}
        laplacian = lambda e: sum(sp.diff(e, q, 2) for q in (x, y, z))
        zero = lambda e: self.assertLess(abs(float(e.subs(at))), 1e-13)
        for name in ("U", "psi", "V_1", "V_2", "V_3"):
            zero(laplacian(f[name]))
        zero(laplacian(f["chi"]) - 2 * f["U"])
        zero(sp.diff(f["U"], t) + sp.diff(f["V_1"], x) + sp.diff(f["V_2"], y) + sp.diff(f["V_3"], z))

    def test_two_bodies_in_the_harmonic_chart_are_blanchet_faye_and_ponsots_metric(self):
        """Their (7.2) through the first post-Newtonian order, with G = c = 1: g_00 + 1 is
        2 m_1/r_1 + (m_1/r_1)(4 v_1^2 - (n_1 v_1)^2) - 2 m_1^2/r_1^2 + m_1 m_2 (-2/(r_1 r_2)
        - r_1/(2 r_12^3) + r_1^2/(2 r_2 r_12^3) - 5/(2 r_2 r_12)) and the same with the bodies
        exchanged, g_0i is -4 m_1 v_1^i/r_1 - 4 m_2 v_2^i/r_2, and g_ij is (1 + 2 m_1/r_1 + 2 m_2/r_2)
        delta_ij. It holds where each body falls toward the other as Newton has it, which is where
        they wrote the accelerations out."""
        sp = self.sp
        m, places, velocities, accelerations, d, event = self.two_bodies()
        space = sp.symbols("x y z", real=True)
        t0 = sp.Symbol("t", real=True)
        _, f = bodies(sp, t0, places, velocities, accelerations, m)
        t, g = self.with_bodies("harmonic", {k: v.subs(t0, sp.Symbol("t", real=True)) for k, v in f.items()}, space)
        at = {**dict(zip(space, event)), t: 0, t0: 0}
        r = [sp.sqrt(sum((q - p) ** 2 for q, p in zip(event, place))) for place in places]
        n = [[(q - p) / r[a] for q, p in zip(event, places[a])] for a in range(2)]
        nv = [sum(a * b for a, b in zip(n[a], velocities[a])) for a in range(2)]
        v2 = [sum(c ** 2 for c in velocities[a]) for a in range(2)]
        theirs = 0
        for a, b in ((0, 1), (1, 0)):
            theirs += (2 * m[a] / r[a] + m[a] / r[a] * (4 * v2[a] - nv[a] ** 2) - 2 * m[a] ** 2 / r[a] ** 2
                       + m[a] * m[b] * (-2 / (r[a] * r[b]) - r[a] / (2 * d ** 3) + r[a] ** 2 / (2 * r[b] * d ** 3)
                                        - sp.Rational(5, 2) / (r[b] * d)))
        self.assertAlmostEqual(float((g[("t", "t")] + 1).subs(at)), float(theirs), places=13)
        for i, k in enumerate("xyz"):
            shift = -4 * sum(m[a] * velocities[a][i] / r[a] for a in range(2))
            self.assertAlmostEqual(float(g[("t", k)].subs(at)), float(shift), places=13)
            self.assertAlmostEqual(float(g[(k, k)].subs(at)), float(1 + 2 * m[0] / r[0] + 2 * m[1] / r[1]), places=13)

    def test_on_a_circular_orbit_the_harmonic_chart_is_alvis_metric_in_corotating_coordinates(self):
        """Alvi's near zone metric, his equation for ds^2 in corotating coordinates, for two bodies a
        distance b apart turning at omega = sqrt(m/b^3), with mu = m_1 m_2/m and eps = sqrt(m/b): the
        published harmonic chart with the potentials of that orbit, carried to the turning coordinates
        at t = 0, where the two charts' axes agree."""
        sp = self.sp
        R = sp.Rational
        m1, m2, b = R(1, 40), R(1, 60), R(3)
        m = m1 + m2
        w = sp.sqrt(m / b ** 3)
        mu, eps2 = m1 * m2 / m, m / b
        t0 = sp.Symbol("t", real=True)
        x, y, z = space = sp.symbols("x y z", real=True)
        a1, a2 = m2 * b / m, m1 * b / m
        paths = ([a1 * sp.cos(w * t0), a1 * sp.sin(w * t0), 0], [-a2 * sp.cos(w * t0), -a2 * sp.sin(w * t0), 0])
        dist = [sp.sqrt((x - p[0]) ** 2 + (y - p[1]) ** 2 + z ** 2) for p in paths]
        speed = [[sp.diff(c, t0) for c in p] for p in paths]
        f = {"U": m1 / dist[0] + m2 / dist[1], "chi": m1 * dist[0] + m2 * dist[1],
             "psi": (m1 / dist[0]) * (R(3, 2) * (w * a1) ** 2 - m2 / b) + (m2 / dist[1]) * (R(3, 2) * (w * a2) ** 2 - m1 / b)}
        for i in range(3):
            f[f"V_{i + 1}"] = m1 * speed[0][i] / dist[0] + m2 * speed[1][i] / dist[1]
        t, g = self.with_bodies("harmonic", f, space)
        event = {x: R(7, 5), y: R(-9, 10), z: R(4, 5), t: 0, t0: 0}
        ours = sp.Matrix(4, 4, lambda i, j: g.get(("txyz"[i], "txyz"[j]), 0)).subs(event)
        # dx' = dx - w y dt and dy' = dy + w x dt at t = 0, with the inertial coordinates primed.
        turn = sp.eye(4)
        turn[1, 0], turn[2, 0] = -w * event[y], w * event[x]
        ours = turn.T * ours * turn
        X, Y, Z = event[x], event[y], event[z]
        r1, r2 = sp.sqrt((X - a1) ** 2 + Y ** 2 + Z ** 2), sp.sqrt((X + a2) ** 2 + Y ** 2 + Z ** 2)
        S = 1 + 2 * m1 / r1 + 2 * m2 / r2
        alvi = sp.zeros(4, 4)
        alvi[0, 0] = (-1 + 2 * m1 / r1 + 2 * m2 / r2 - 2 * (m1 / r1 + m2 / r2) ** 2 + 3 * mu / b * (m2 / r1 + m1 / r2)
                      - mu / b * (m2 / r1 ** 3 + m1 / r2 ** 3) * Y ** 2 - 2 * mu * eps2 * (1 / r1 + 1 / r2)
                      - 7 * mu * eps2 * (1 / r1 - 1 / r2) * X / b + w ** 2 * S * (X ** 2 + Y ** 2))
        alvi[0, 1] = alvi[1, 0] = -w * S * Y
        alvi[0, 2] = alvi[2, 0] = w * S * X - 4 * mu * sp.sqrt(eps2) * (1 / r1 - 1 / r2)
        alvi[1, 1] = alvi[2, 2] = alvi[3, 3] = S
        for i in range(4):
            for j in range(4):
                self.assertAlmostEqual(float(ours[i, j]), float(alvi[i, j]), places=12, msg=(i, j))

    def test_one_body_at_rest_is_the_published_post_newtonian_metric_at_general_relativitys_values(self):
        """With U = m/r, chi = m r and nothing moving, both charts are the Cartesian chart of the
        parametrised post-Newtonian metric at beta = gamma = 1, which is Schwarzschild's in isotropic
        coordinates far from the horizon."""
        sp = self.sp
        m = sp.Symbol("m", positive=True)
        x, y, z = space = sp.symbols("x y z", real=True)
        r = sp.sqrt(x ** 2 + y ** 2 + z ** 2)
        f = {"U": m / r, "psi": sp.Integer(0), "chi": m * r, "V_1": sp.Integer(0), "V_2": sp.Integer(0), "V_3": sp.Integer(0)}
        ppn = json.loads((DATA / "metrics" / "ppn_metric.json").read_text(encoding="utf-8"))
        entry = next(c for c in ppn["coordinates"] if c["id"] == "cartesian")
        theirs = self.vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
        same = {theirs.symbol[k]: s for k, s in zip("xyz", space)}
        same.update({theirs.parameters["m"]: m, theirs.parameters["beta"]: 1, theirs.parameters["gamma"]: 1})
        published = {tuple(e["indices"]): sp.sympify(theirs(e["value"])).subs(theirs.held).doit().subs(same)
                     for e in entry["metric_components"]}
        for system in ("harmonic", "standard"):
            _, g = self.with_bodies(system, f, space)
            for i in "txyz":
                for j in "txyz":
                    self.assertEqual(sp.simplify(g.get((i, j), 0) - published.get((i, j), 0)), 0, (system, i, j))

    def test_the_inverse_printed_whole_is_the_inverse_of_the_metric(self):
        sp = self.sp
        for system in ("harmonic", "standard"):
            reader, g = self.components(system)
            _, gi = self.components(system, "inverse_metric_components")
            matrix = lambda d: sp.Matrix(4, 4, lambda i, j: d.get(("txyz"[i], "txyz"[j]), 0))
            product = (matrix(g) * matrix(gi)).applyfunc(sp.simplify)
            self.assertEqual(product, sp.eye(4), system)

    def test_the_kretschmann_scalar_of_one_body_is_schwarzschilds(self):
        """K is 8 d_i d_j U d_i d_j U + 4 (Laplacian of U)^2, and for U = m/r that is 48 m^2/r^6."""
        sp = self.sp
        entry = next(c for c in self.metric["coordinates"] if c["id"] == "harmonic")
        reader = self.vm.Reader(entry["coords"], [p["symbol"] for p in entry["parameters"]], ())
        K = sp.sympify(reader(entry["kretschmann"].partition("=")[2]))
        t, x, y, z = (reader.symbol[k] for k in "txyz")
        U = reader.parameters["U"]
        m = sp.Symbol("m", positive=True)
        r = sp.sqrt(x ** 2 + y ** 2 + z ** 2)
        one = K.replace(U.func, sp.Lambda((t, x, y, z), m / r)).doit()
        self.assertEqual(sp.simplify(one - 48 * m ** 2 / r ** 6), 0)


@unittest.skipUnless(HAS_SYMPY, "sympy, numpy and scipy are not installed")
class Motion(unittest.TestCase):
    """The equations of motion of Einstein, Infeld and Hoffmann, their (17.2), with G = c = 1."""

    @classmethod
    def setUpClass(cls):
        import numpy
        import sympy as sp
        cls.np, cls.sp = numpy, sp
        eta = sp.symbols("eta1:4", real=True)
        zeta = sp.symbols("zeta1:4", real=True)
        de = sp.symbols("deta1:4", real=True)
        dz = sp.symbols("dzeta1:4", real=True)
        m1, m2 = sp.symbols("m1 m2", positive=True)
        r = sp.sqrt(sum((a - b) ** 2 for a, b in zip(eta, zeta)))
        dot = lambda a, b: sum(p * q for p, q in zip(a, b))
        bracket = dot(de, de) + sp.Rational(3, 2) * dot(dz, dz) - 4 * dot(de, dz) - 4 * m2 / r - 5 * m1 / r
        first = []
        for k in range(3):
            value = m2 * sp.diff(1 / r, eta[k]) + m2 * (
                bracket * sp.diff(1 / r, eta[k])
                + sum((4 * de[s] * (dz[k] - de[k]) + 3 * de[k] * dz[s] - 4 * dz[s] * dz[k]) * sp.diff(1 / r, eta[s])
                      for s in range(3))
                + sp.Rational(1, 2) * sum(sp.diff(r, eta[s], eta[q], eta[k]) * dz[s] * dz[q]
                                          for s in range(3) for q in range(3)))
            first.append(value)
        # The other body's equation is the same with the two exchanged.
        swap = {**dict(zip(eta, zeta)), **dict(zip(zeta, eta)), **dict(zip(de, dz)), **dict(zip(dz, de)), m1: m2, m2: m1}
        second = [e.subs(swap, simultaneous=True) for e in first]
        cls.arguments = (*eta, *zeta, *de, *dz, m1, m2)
        cls.acceleration = staticmethod(sp.lambdify(cls.arguments, [first, second], "numpy"))

    def relative(self, x, v, m1, m2):
        """The acceleration of the first body less that of the second, in the frame of their centre
        of mass, for the separation x and the relative velocity v."""
        np = self.np
        m = m1 + m2
        a = self.acceleration(*(m2 / m * x), *(-m1 / m * x), *(m2 / m * v), *(-m1 / m * v), m1, m2)
        return np.array(a[0]) - np.array(a[1])

    def test_the_two_body_acceleration_is_the_one_of_wills_review(self):
        """Will's (2014) relative acceleration to first post-Newtonian order,
        (m/r^2)(-n + ((4 + 2 eta) m/r - (1 + 3 eta) v^2 + (3/2) eta rdot^2) n + (4 - 2 eta) rdot v)."""
        np = self.np
        m1, m2 = 0.03, 0.011
        m, eta = m1 + m2, m1 * m2 / (m1 + m2) ** 2
        x, v = np.array([3.0, -1.5, 0.7]), np.array([0.04, 0.09, -0.05])
        r = np.linalg.norm(x)
        n = x / r
        rdot = n @ v
        will = m / r ** 2 * (-n + ((4 + 2 * eta) * m / r - (1 + 3 * eta) * (v @ v) + 1.5 * eta * rdot ** 2) * n
                             + (4 - 2 * eta) * rdot * v)
        self.assertLess(np.abs(self.relative(x, v, m1, m2) - will).max(), 1e-15)

    def test_a_circular_orbit_turns_at_the_rate_the_drawings_declare(self):
        """A circle of separation d closes where Omega^2 d = (m/d^2)(1 + (1 + 3 eta) Omega^2 d^2 -
        (4 + 2 eta) m/d), which to first order in m/d is Omega^2 = (m/d^3)(1 - (3 - eta) m/d). For two
        bodies of mass 1 a distance 20 apart that is 29/160000, each moving at sqrt(29)/40, and psi is
        (3 v^2/2 - 1/20) U = -(73/3200) U, the declared input of every drawing."""
        import null_rays as nr
        sp = self.sp
        m, d, eta, W = sp.symbols("m d eta W", positive=True)       # W is Omega^2
        closes = sp.solve(sp.Eq(W * d, m / d ** 2 * (1 + (1 + 3 * eta) * W * d ** 2 - (4 + 2 * eta) * m / d)), W)[0]
        eps = sp.Symbol("eps", positive=True)
        first = sp.series(closes.subs(m, eps * m), eps, 0, 3).removeO().subs(eps, 1)
        self.assertEqual(sp.simplify(first - m / d ** 3 * (1 - (3 - eta) * m / d)), 0)
        rate = (m / d ** 3 * (1 - (3 - eta) * m / d)).subs({m: 2, d: 20, eta: sp.Rational(1, 4)})
        self.assertEqual(rate, sp.Rational(29, 160000))
        v2 = rate * 100
        self.assertEqual(sp.Rational(3, 2) * v2 - sp.Rational(1, 20), -sp.Rational(73, 3200))
        t, x, y, z = sp.symbols("t x y z", real=True)
        names = {"t": t, "x": x, "y": y, "z": z}
        ours = {k: sp.sympify(v, locals=names) for k, v in nr.EIH_BINARY.items()}
        w = sp.sqrt(rate)
        paths = ((10 * sp.cos(w * t), 10 * sp.sin(w * t), 0), (-10 * sp.cos(w * t), -10 * sp.sin(w * t), 0))
        dist = [sp.sqrt((x - p[0]) ** 2 + (y - p[1]) ** 2 + (z - p[2]) ** 2) for p in paths]
        theirs = {"U": 1 / dist[0] + 1 / dist[1], "chi": dist[0] + dist[1],
                  "psi": -sp.Rational(73, 3200) * (1 / dist[0] + 1 / dist[1]), "V_3": sp.Integer(0)}
        for i in range(2):
            theirs[f"V_{i + 1}"] = sum(sp.diff(p[i], t) / r for p, r in zip(paths, dist))
        at = {t: sp.Rational(37, 3), x: sp.Rational(3, 2), y: sp.Rational(-7, 2), z: sp.Rational(5, 3)}
        for name in POTENTIALS:
            self.assertAlmostEqual(float(ours[name].subs(at)), float(theirs[name].subs(at)), places=13, msg=name)

    def test_an_eccentric_orbit_turns_by_six_pi_m_over_p_each_revolution(self):
        """Robertson's result, the Einstein-Robertson formula: the pericentre of two bodies of
        comparable mass advances by 6 pi m/p a revolution, with m the sum of the masses and p the
        semilatus rectum. Integrated from (17.2) for masses 0.6 and 0.4 on an orbit of p = 400 and
        e = 0.4, between two passages of the pericentre."""
        from scipy.integrate import solve_ivp
        np = self.np
        m1, m2 = 0.6, 0.4
        p, e = 400.0, 0.4
        start = p / (1 + e)
        speed = math.sqrt((m1 + m2) / p) * (1 + e)

        def rhs(_, s):
            return np.concatenate([s[3:], self.relative(s[:3], s[3:], m1, m2)])

        inward = lambda _, s: float(np.dot(s[:3], s[3:]))
        inward.direction = 1.0                  # the separation stops falling: a pericentre
        out = solve_ivp(rhs, (0.0, 1.0e5), np.array([start, 0.0, 0.0, 0.0, speed, 0.0]), events=inward, rtol=1e-11, atol=1e-13)
        passages = [s for t, s in zip(out.t_events[0], out.y_events[0]) if t > 1.0]
        advance = math.atan2(passages[0][1], passages[0][0]) % (2 * math.pi)
        self.assertAlmostEqual(advance / (6 * math.pi * (m1 + m2) / p), 1.0, delta=0.03)

    def test_the_numbers_of_the_history_and_the_captions(self):
        # The pulsar's orbit against Mercury's, and its period in hours.
        self.assertEqual(round(4.226598 * 3600 * 100 / 42.98, -3), 35000)
        self.assertEqual(round(0.322997448911 * 24, 2), 7.75)
        # On the axis of the declared binary, midway between the bodies, U = 1/5.
        u = 2 / 10
        self.assertEqual(round(math.sqrt((1 - 2 * u + 2 * u * u + 2 * 73 / 3200 * u) / (1 + 2 * u)), 2), 0.70)
        # The near zone ends at c/(2 Omega), Alvi's outer limit.
        self.assertEqual(round(1 / (2 * math.sqrt(29 / 160000)), 1), 37.1)


@unittest.skipUnless(HAS_SYMPY, "sympy, numpy and scipy are not installed")
class Orders(unittest.TestCase):
    """The post-Newtonian orders the checker keeps where the small quantities are potentials."""

    @classmethod
    def setUpClass(cls):
        import sympy
        import verify_metrics
        cls.sp, cls.vm = sympy, verify_metrics

    def test_a_potential_and_each_of_its_time_derivatives_carry_their_orders(self):
        sp, vm = self.sp, self.vm
        reader = vm.Reader(["t", "x"], ["U = U(t,x)", "m"], (), kept=({"U": 2, "m": 2}, 2, "t"))
        t, x = reader.symbol["t"], reader.symbol["x"]
        U, m = reader.parameters["U"], reader.parameters["m"]
        D = sp.Derivative
        value = U + U ** 2 + U ** 3 + D(U, x) + D(U, t) + D(U, t, t) * U + D(U, (t, 3)) + m * D(U, t)
        self.assertEqual(sp.expand(reader.through(value, 4) - (U + U ** 2 + D(U, x) + D(U, t))), 0)
        self.assertEqual(sp.expand(reader.through(value, 5) - (U + U ** 2 + D(U, x) + D(U, t) + D(U, (t, 3))
                                                              + m * D(U, t))), 0)
        self.assertEqual(sp.expand(reader.through(value, 2) - (U + D(U, x))), 0)
        # A system that names no slow time counts constants alone.
        with self.assertRaises(vm.LatexError):
            vm.Reader(["t", "x"], ["U = U(t,x)"], (), kept=({"U": 2}, 2))

    def test_terms_of_the_next_order_in_the_metric_move_nothing_that_is_kept(self):
        """With an arbitrary function of the four coordinates put in each of the ten components at
        the first order the post-Newtonian metric leaves open, sixth in g_tt, fifth in g_ti and
        fourth in g_ij, every tensor the checker hands over for the harmonic chart is what it was,
        time dependent potentials and all."""
        import print_charts
        sp, vm = self.sp, self.vm
        spec = print_charts.eih_many_bodies("harmonic")
        coords = spec["system"]["coords"]
        weights = dict(vm.ORDERS[("eih_many_bodies", "harmonic")][0])

        def geometry(open_terms):
            names = [f"F_{i}{j} = F_{i}{j}(t,x,y,z)" for i in range(4) for j in range(i, 4)] if open_terms else []
            extra = {f"F_{i}{j}": 6 - (i != 0) - (j != 0) for i in range(4) for j in range(i, 4)} if open_terms else {}
            reader = vm.Reader(coords, spec["system"]["parameters"] + names, (), kept=({**weights, **extra}, 2, "t"))
            g = vm.metric_from_line_element(reader, spec["chart_line_element"], coords)
            symbols = [reader.symbol[name] for name in coords]
            for i in range(4):
                for j in range(i, 4):
                    if open_terms:
                        g[i, j] += reader.parameters[f"F_{i}{j}"]
                        g[j, i] = g[i, j]
            return vm.geometry_of(g, symbols, 10 ** 6, reader)

        fixed, opened = geometry(False), geometry(True)
        for name in ("christoffel_ull", "riemann_llll", "ricci_ll", "einstein_ll", "weyl_llll", "ricci_scalar"):
            ours, theirs = getattr(fixed, name)(), getattr(opened, name)()
            self.assertTrue(all(vm.norm(a - b) == 0 for a, b in zip(flat(ours), flat(theirs))), name)


class Texts(unittest.TestCase):
    """What the files say, without sympy."""

    @classmethod
    def setUpClass(cls):
        cls.metric = json.loads((DATA / "metrics" / "eih_many_bodies.json").read_text(encoding="utf-8"))

    def test_both_charts_keep_the_same_potentials(self):
        charts = {c["id"]: c for c in self.metric["coordinates"]}
        self.assertEqual(list(charts), ["harmonic", "standard"])
        symbols = [[p["symbol"] for p in c["parameters"]] for c in charts.values()]
        self.assertEqual(symbols[0], symbols[1])
        self.assertIn("\\partial_t^2\\chi", charts["harmonic"]["line_element"])
        self.assertNotIn("\\partial_t^2\\chi", charts["standard"]["line_element"])
        self.assertIn("\\partial_t\\partial_x\\chi", charts["standard"]["line_element"])

    def test_the_drawings_hold_the_declared_binary(self):
        diagrams = json.loads((DATA / "diagrams" / "eih_many_bodies.json").read_text(encoding="utf-8"))
        embedding = json.loads((DATA / "embedding" / "eih_many_bodies.json").read_text(encoding="utf-8"))
        text = json.dumps(diagrams) + json.dumps(embedding)
        self.assertIn("d = 20\\\\,m", text)
        self.assertEqual(sorted(diagrams["systems"]), ["harmonic", "standard"])
        self.assertEqual([v["id"] for v in embedding["views"]], ["midplane"])


if __name__ == "__main__":
    unittest.main()
