import unittest

from pqc_values import PQC_Values

class PQCValuesTest(unittest.TestCase):

    def test_properties(self):
        expected_value = .1
        stray = .5
        default_min_allowed = expected_value * (1 - stray)
        default_max_allowed = expected_value * (1 + stray)

        values = PQC_Values(expected_value=.1)

        self.assertEqual(values.values, [])
        self.assertEqual(values.expected_value, expected_value)
        self.assertEqual(values.stray, stray)
        self.assertEqual(values.min_allowed, default_min_allowed)
        self.assertEqual(values.max_allowed, default_max_allowed)

        values.min_allowed = 42.
        self.assertEqual(values.min_allowed, 42.)
        self.assertEqual(values.max_allowed, default_max_allowed)

        values.max_allowed = 43.
        self.assertEqual(values.min_allowed, 42.)
        self.assertEqual(values.max_allowed, 43.)

        values.min_allowed = None
        self.assertEqual(values.min_allowed, default_min_allowed)
        self.assertEqual(values.max_allowed, 43.)

        values.max_allowed = None
        self.assertEqual(values.min_allowed, default_min_allowed)
        self.assertEqual(values.max_allowed, default_max_allowed)

        stray = .25
        default_min_allowed = expected_value * (1 - stray)
        default_max_allowed = expected_value * (1 + stray)

        values.stray = stray
        self.assertEqual(values.expected_value, expected_value)
        self.assertEqual(values.stray, stray)
        self.assertEqual(values.min_allowed, default_min_allowed)
        self.assertEqual(values.max_allowed, default_max_allowed)
