import unittest

from validation.jkr_tools.plotting import sweep_series


class SweepSeriesTests(unittest.TestCase):
    def valid_row(self, velocity):
        return {
            "impact_velocity": velocity,
            "rebound_velocity": 0.5 * velocity,
            "effective_restitution": 0.5,
            "captured": 0.0,
            "maximum_overlap": velocity * 1.0e-9,
            "contact_duration": 1.0e-7,
        }

    def test_series_are_sorted_by_impact_velocity(self):
        series = sweep_series((self.valid_row(2.0), self.valid_row(1.0)))
        self.assertEqual(series["impact_velocity"], (1.0, 2.0))
        self.assertEqual(series["rebound_velocity"], (0.5, 1.0))

    def test_missing_field_is_rejected(self):
        row = self.valid_row(1.0)
        del row["captured"]
        with self.assertRaisesRegex(ValueError, "missing fields"):
            sweep_series((row,))

    def test_duplicate_impact_velocity_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "unique"):
            sweep_series((self.valid_row(1.0), self.valid_row(1.0)))


if __name__ == "__main__":
    unittest.main()
