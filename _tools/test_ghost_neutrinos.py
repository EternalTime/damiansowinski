"""The ghost neutrinos of Davis and Ray, held to the stored charts and drawings.

Pin Taub's vacuum in every chart, the neutrino field that solves the massless Dirac
equation with no energy or momentum when its current runs along the axis, the flux
whose three factors cancel, the length L that a rescaling removes, the Kasner
continuation through the singularity, and the curvature at the two ends of the chart.
Run only this module with the collection's Python environment.
"""
import importlib.util
import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / 'MFS' / 'assets' / 'data'
DERIVATIONS = ROOT / '_tools' / 'derivations'
HAS_SYMPY = importlib.util.find_spec('sympy') is not None
sys.path.insert(0, str(DERIVATIONS))
SYSTEMS = ('taub', 'areal', 'null')


def data(kind, name='ghost_neutrinos'):
    return json.loads((DATA / kind / (name + '.json')).read_text())


def chart(system):
    """The stored chart read back: its reader and its metric with c = 1."""
    import sympy as sp
    import verify_metrics as vm
    entry = next(c for c in data('metrics')['coordinates'] if c['id'] == system)
    reader = vm.Reader(entry['coords'], [p['symbol'] for p in entry['parameters']], ())
    g = sp.zeros(4)
    for item in entry['metric_components']:
        i, j = (entry['coords'].index(c) for c in item['indices'])
        g[i, j] = reader(item['value']).subs(reader.c, 1)
    return entry, reader, g


def taub_z(system, reader):
    """Taub's proper distance z in the chart's coordinates, from the chart's convention."""
    import sympy as sp
    L = reader.parameters['L']
    if system == 'taub':
        return reader.symbol['z']
    if system == 'areal':
        r = reader.symbol['r']
    else:
        r = L * sp.sqrt(2 * (reader.symbol['\\mu'] + reader.symbol['\\nu']))
    return sp.Rational(2, 3) * L * (r / L) ** sp.Rational(3, 2)


@unittest.skipUnless(HAS_SYMPY, 'sympy is not installed')
class TaubVacuum(unittest.TestCase):
    """The geometry around the neutrinos is empty space, in every chart."""

    POINTS = {'taub': [(0.3, -0.7, 1.1, 0.4), (-0.8, 2.0, 0.3, 7.5)],
              'areal': [(0.3, -0.7, 1.1, 0.4), (-0.8, 2.0, 0.3, 7.5)],
              'null': [(0.2, 0.5, 0.3, -0.4), (3.0, -1.2, -0.5, 0.9)]}

    def test_every_chart_is_ricci_flat_at_points(self):
        import verify_metrics as vm
        for system in SYSTEMS:
            entry, reader, g = chart(system)
            geo = vm.Geometry(g, [reader.symbol[c] for c in entry['coords']], 60)
            ricci = geo.ricci_ll()
            coords = [reader.symbol[c] for c in entry['coords']]
            for values in self.POINTS[system]:
                at = dict(zip(coords, values))
                at[reader.parameters['L']] = 1.3
                for i in range(4):
                    for j in range(4):
                        self.assertLess(abs(complex(vm._at(ricci, (i, j)).subs(at).evalf(30))), 1e-20,
                                        f'{system}: R_{i}{j} at {values}')

    def test_stored_kretschmann_is_64_over_27_z_to_the_fourth_in_every_chart(self):
        for system in SYSTEMS:
            entry, reader, _ = chart(system)
            stored = reader(entry['kretschmann'].removeprefix('K = '))
            z = taub_z(system, reader)
            coords = [reader.symbol[c] for c in entry['coords']]
            for values in self.POINTS[system]:
                at = dict(zip(coords, values))
                at[reader.parameters['L']] = 1.3
                want = complex((64 / (27 * z ** 4)).subs(at).evalf(30))
                self.assertAlmostEqual(complex(stored.subs(at).evalf(30)) / want, 1, places=12, msg=system)

    def test_curvature_diverges_at_the_singularity_and_dies_away_far_from_it(self):
        entry, reader, _ = chart('taub')
        stored = reader(entry['kretschmann'].removeprefix('K = '))
        z = reader.symbol['z']
        import sympy as sp
        self.assertEqual(sp.limit(stored, z, 0, '+'), sp.oo)
        self.assertEqual(sp.limit(stored, z, sp.oo), 0)

    def test_the_length_L_is_removed_by_rescaling_t_x_and_y(self):
        import sympy as sp
        _, reader, g = chart('taub')
        L, z = reader.parameters['L'], reader.symbol['z']
        # dt -> L^(-1/3) dt and dx, dy -> L^(2/3) dx, dy carry the metric at L to the metric at L = 1.
        scale = sp.diag(L ** sp.Rational(-1, 3), L ** sp.Rational(2, 3), L ** sp.Rational(2, 3), 1)
        for Lv, zv in ((0.4, 0.7), (2.5, 3.1)):
            at = {L: Lv, z: zv}
            rescaled = (scale * g * scale).subs(at).evalf(30)
            unit = g.subs({L: 1, z: zv}).evalf(30)
            self.assertLess(max(abs(complex(v)) for v in rescaled - unit), 1e-25)

    def test_continuing_r_through_zero_gives_kasner_with_exponents_two_thirds_two_thirds_minus_one_third(self):
        import sympy as sp
        _, reader, g = chart('areal')
        L, r = reader.parameters['L'], reader.symbol['r']
        s = sp.Symbol('s', positive=True)  # s = -r > 0 beyond the singularity
        g = g.subs({r: -s, L: sp.Symbol('L', positive=True)})
        self.assertTrue(bool(g[3, 3] < 0) and bool(g[0, 0] > 0))  # r is now the time
        tau = sp.integrate(sp.sqrt(-g[3, 3]), (s, 0, s))  # proper time from the singularity
        exponents = []
        for i in range(3):
            # g_ii as a power of tau: d ln g_ii / d ln tau, which is 2p.
            exponents.append(sp.simplify(sp.diff(sp.log(g[i, i]), s) / sp.diff(sp.log(tau), s)) / 2)
        self.assertEqual(exponents, [sp.Rational(-1, 3), sp.Rational(2, 3), sp.Rational(2, 3)])
        self.assertEqual(sum(exponents), 1)
        self.assertEqual(sum(p ** 2 for p in exponents), 1)


@unittest.skipUnless(HAS_SYMPY, 'sympy is not installed')
class GhostNeutrino(unittest.TestCase):
    """The neutrino field in Taub's chart, with the tetrad diag(cA, C, C, 1) read off the stored metric."""

    @classmethod
    def setUpClass(cls):
        import sympy as sp
        cls.sp = sp
        entry, reader, g = chart('taub')
        # The chart's z and L are positive, which the reader does not declare; real t, x and y.
        L, z = sp.symbols('L z', positive=True)
        t, x, y = sp.symbols('t x y', real=True)
        swap = {reader.parameters['L']: L, reader.symbol['z']: z, reader.symbol['t']: t,
                reader.symbol['x']: x, reader.symbol['y']: y}
        g = g.subs(swap)
        X = [swap[reader.symbol[c]] for c in entry['coords']]
        cls.reader, cls.g, cls.X, cls.L, cls.z = reader, g, X, L, z
        # Orthonormal coframe from the stored diagonal metric, in the chart order (t, x, y, z).
        E = sp.diag(*[sp.sqrt(abs(g[i, i])) for i in range(4)]).subs({L: L, z: z})
        E = E.applyfunc(lambda e: sp.powsimp(sp.powdenest(e, force=True), force=True))
        cls.Einv = E.inv()
        cls.E = E
        eta = sp.diag(-1, 1, 1, 1)
        I2, Z2 = sp.eye(2), sp.zeros(2)
        pauli = [sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -sp.I], [sp.I, 0]]), sp.Matrix([[1, 0], [0, -1]])]

        def block(a, b, c, d):
            return sp.BlockMatrix([[a, b], [c, d]]).as_explicit()
        # Chiral gammas of signature (+,-,-,-), times i for (-,+,+,+); x, y, z are the frame's 1, 2, 3.
        cls.gt = [block(Z2, I2, I2, Z2)] + [block(Z2, s, -s, Z2) for s in pauli]
        gam = [sp.I * m for m in cls.gt]
        ginv = sp.diag(*[1 / g[i, i] for i in range(4)])
        Gam = [[[sum(ginv[m, s] * (sp.diff(g[s, n], X[r]) + sp.diff(g[s, r], X[n]) - sp.diff(g[n, r], X[s]))
                     for s in range(4)) / 2 for r in range(4)] for n in range(4)] for m in range(4)]
        lower = [sum((eta[a, b] * gam[b] for b in range(4)), sp.zeros(4)) for a in range(4)]
        spin = []
        for mu in range(4):
            S = sp.zeros(4)
            for a in range(4):
                for b in range(4):
                    w = sum(E[a, nu] * (sp.diff(cls.Einv[nu, b], X[mu])
                                        + sum(Gam[nu][mu][lam] * cls.Einv[lam, b] for lam in range(4)))
                            for nu in range(4))
                    S += sp.simplify(w) * lower[a] * gam[b] / 4
            spin.append(S)
        cls.spin = spin
        cls.gup = [sum((cls.Einv[mu, a] * gam[a] for a in range(4)), sp.zeros(4)) for mu in range(4)]
        cls.gdown = [sum((g[mu, nu] * cls.gup[nu] for nu in range(4)), sp.zeros(4)) for mu in range(4)]

    def field(self, chi):
        sp = self.sp
        return (self.L / self.z) ** sp.Rational(1, 2) * sp.Matrix(chi)

    def dirac(self, psi):
        sp = self.sp
        return sum((self.gup[m] * (sp.diff(psi, self.X[m]) + self.spin[m] * psi) for m in range(4)), sp.zeros(4, 1))

    def bar(self, psi):
        return psi.H * self.gt[0]

    def stress(self, psi):
        """The symmetric energy-momentum of the field, up to its constant factor."""
        sp = self.sp
        X = self.X

        def nab(m):
            return sp.diff(psi, X[m]) + self.spin[m] * psi

        def nabbar(m):
            return sp.diff(self.bar(psi), X[m]) - self.bar(psi) * self.spin[m]
        return sp.Matrix(4, 4, lambda m, n: sp.simplify(
            (self.bar(psi) * self.gdown[m] * nab(n) - nabbar(n) * self.gdown[m] * psi
             + self.bar(psi) * self.gdown[n] * nab(m) - nabbar(m) * self.gdown[n] * psi)[0]))

    def current(self, psi):
        return [self.sp.simplify((self.bar(psi) * self.sp.I * self.gup[m] * psi)[0]) for m in range(4)]

    def test_the_field_solves_the_massless_dirac_equation_for_every_constant_spinor(self):
        sp = self.sp
        chi = sp.symbols('a0:4')
        residual = self.dirac(self.field(chi))
        self.assertEqual(sp.simplify(residual), sp.zeros(4, 1))

    def test_energy_momentum_vanishes_when_the_current_runs_along_z(self):
        sp = self.sp
        for chi in ((1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 2, 0), (0, 0, 0, sp.I)):
            psi = self.field(chi)
            self.assertEqual(self.stress(psi), sp.zeros(4), chi)
            j = self.current(psi)
            self.assertEqual(j[1:3], [0, 0], chi)
            self.assertNotEqual(j[3], 0, chi)

    def test_energy_momentum_does_not_vanish_when_the_current_leans_off_the_axis(self):
        sp = self.sp
        psi = self.field((1, 1, 0, 0))  # left handed, spin along x
        j = self.current(psi)
        self.assertNotEqual(j[1], 0)
        self.assertNotEqual(self.stress(psi), sp.zeros(4))

    def test_the_current_is_null_and_its_flux_is_the_same_on_every_plane(self):
        sp = self.sp
        g, z = self.g, self.z
        for chi in ((1, 0, 0, 0), (0, 1, 0, 0)):
            j = self.current(self.field(chi))
            self.assertEqual(sp.simplify(sum(g[m, m] * j[m] ** 2 for m in range(4))), 0)
            root = sp.sqrt(-g[0, 0] * g[1, 1] * g[2, 2] * g[3, 3])
            flux = sp.simplify(sp.powsimp(root * j[3], force=True))
            self.assertEqual(sp.diff(flux, z), 0, chi)
            self.assertNotEqual(flux, 0, chi)

    def test_the_three_factors_of_the_flux_cancel(self):
        sp = self.sp
        z, L = self.z, self.L
        density = sp.powsimp(self.current(self.field((1, 0, 0, 0)))[3], force=True)
        area = self.g[1, 1]  # area of a patch of plane, per dx dy
        clock = sp.sqrt(-self.g[0, 0])  # rate of a static clock, per dt
        exponents = [sp.simplify(sp.diff(sp.log(f), z) * z) for f in (density, area, clock)]
        self.assertEqual(exponents, [-1, sp.Rational(4, 3), sp.Rational(-1, 3)])
        self.assertEqual(sum(exponents), 0)
        self.assertEqual(sp.simplify(sp.powsimp(density * area * clock, force=True)).has(z), False)
        self.assertTrue(L in sp.sympify(density).free_symbols)


if __name__ == '__main__':
    unittest.main()
