import unittest

from sympy import Symbol, exp

from pyProximation import Interpolation, Measure, OrthSystem, RationalAprox, RationalApprox


class MeasureTests(unittest.TestCase):
    def test_discrete_integral_accepts_tuple_argument_styles(self):
        measure = Measure({(0, 0): 1, (1, 2): 3})

        self.assertEqual(measure.integral(lambda point: point[0] + point[1]), 9)
        self.assertEqual(measure.integral(lambda x, y: x + y), 9)

    def test_discrete_sampling_rejects_negative_weights(self):
        measure = Measure({0: 1, 1: -1})

        with self.assertRaisesRegex(ValueError, 'negative weights'):
            measure.sample(2)

    def test_continuous_integral_and_norm_with_constant_weight(self):
        measure = Measure([(0, 1)], 2)

        self.assertAlmostEqual(measure.integral(lambda x: x), 1.0)
        self.assertAlmostEqual(measure.norm(2, lambda x: x), (2.0 / 3.0) ** 0.5)


class OrthSystemTests(unittest.TestCase):
    def test_dependent_basis_raises_clear_error(self):
        x = Symbol('x')
        system = OrthSystem([x], [(-1, 1)])
        system.Basis([1, 1])

        with self.assertRaisesRegex(ValueError, 'linearly dependent'):
            system.FormBasis()

    def test_fourier_basis_rejects_degenerate_interval(self):
        x = Symbol('x')
        system = OrthSystem([x], [(0, 0)])

        with self.assertRaisesRegex(ValueError, 'non-degenerate intervals'):
            system.FourierBasis(2)


class InterpolationTests(unittest.TestCase):
    def test_dimension_mismatch_raises_clear_error(self):
        x = Symbol('x')
        interpolation = Interpolation([x])

        with self.assertRaisesRegex(ValueError, '1-dimensional'):
            interpolation.Interpolate([(0, 1)], [1])


class RationalTests(unittest.TestCase):
    def test_public_alias_matches_legacy_class(self):
        x = Symbol('x')
        system = OrthSystem([x], [(-1, 1)])
        system.Basis(system.PolyBasis(4))
        system.FormBasis()

        legacy = RationalAprox(system).RatLSQ(2, 2, exp(x))
        alias = RationalApprox(system).RatLSQ(2, 2, exp(x))

        self.assertEqual(legacy, alias)

    def test_requires_formed_basis(self):
        x = Symbol('x')
        system = OrthSystem([x], [(-1, 1)])

        with self.assertRaisesRegex(ValueError, r'Call `FormBasis\(\)`'):
            RationalAprox(system).RatLSQ(1, 1, exp(x))


if __name__ == '__main__':
    unittest.main()