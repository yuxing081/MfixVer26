"""Unit tests for the in-memory validation comparison primitives."""

from __future__ import annotations

import unittest

from validation.jkr_tools.comparison import (
    Tolerance,
    align_by_contact_pair,
    align_by_global_particle_id,
    canonical_contact_pair,
    compare_contact_records,
    compare_field_sequences,
    compare_fields,
    compare_particle_records,
    make_tolerance,
)


class ComparisonTest(unittest.TestCase):
    def test_absolute_and_relative_tolerance_and_summary(self) -> None:
        result = compare_field_sequences(
            [1.0, 2.0, 100.0], [1.0001, 2.001, 100.5], "position_x",
            make_tolerance(1.0e-3, 1.0e-2),
            [0.0, 0.1, 0.2], [0.0, 0.1, 0.2],
        )
        self.assertTrue(result.passed)
        self.assertIsNone(result.first_difference)
        self.assertEqual(result.field_summaries[0].time, 0.2)
        self.assertAlmostEqual(result.field_summaries[0].maximum_absolute_error, 0.5)

    def test_first_difference_contains_required_context(self) -> None:
        result = compare_field_sequences(
            [1.0, 2.0], [1.0, 2.5], "velocity_x", Tolerance(0.1, 0.0),
            [0.0, 0.25], [0.0, 0.25],
        )
        self.assertFalse(result.passed)
        difference = result.first_difference
        assert difference is not None
        self.assertEqual(difference.field, "velocity_x")
        self.assertEqual(difference.index, 1)
        self.assertEqual(difference.time, 0.25)
        self.assertEqual(difference.left, 2.0)
        self.assertEqual(difference.right, 2.5)
        self.assertEqual(difference.absolute_error, 0.5)
        self.assertEqual(difference.tolerance, 0.1)
        self.assertEqual(result.field_summaries[0].index, 1)

    def test_fields_report_first_failure_and_each_maximum(self) -> None:
        result = compare_fields(
            {"x": (0.0, 1.0), "y": (0.0, 2.0)},
            {"x": (0.0, 1.2), "y": (0.4, 2.0)},
            {"x": (0.05, 0.0), "y": (0.1, 0.0)},
            (0.0, 1.0), (0.0, 1.0),
        )
        self.assertFalse(result.passed)
        self.assertEqual(result.first_difference.field, "y")
        self.assertEqual(result.first_difference.index, 0)
        self.assertAlmostEqual(
            next(summary.maximum_absolute_error for summary in result.field_summaries if summary.field == "x"),
            0.2,
        )
        self.assertAlmostEqual(
            next(summary.maximum_absolute_error for summary in result.field_summaries if summary.field == "y"),
            0.4,
        )

    def test_time_mismatch_is_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "time mismatch"):
            compare_field_sequences([0.0], [0.0], "force", (0.0, 0.0), [1.0], [1.1])

    def test_particle_alignment_rejects_duplicates_and_missing_ids(self) -> None:
        left = [{"global_particle_id": 2}, {"global_particle_id": 1}]
        right = [{"global_particle_id": 1}, {"global_particle_id": 2}]
        aligned = align_by_global_particle_id(left, right)
        self.assertEqual([pair[1]["global_particle_id"] for pair in aligned], [2, 1])
        with self.assertRaisesRegex(ValueError, "duplicate"):
            align_by_global_particle_id(left, right + [{"global_particle_id": 2}])
        with self.assertRaisesRegex(ValueError, "particle IDs differ"):
            align_by_global_particle_id(left, [{"global_particle_id": 1}])

    def test_particle_comparison_uses_global_id_not_row_order(self) -> None:
        result = compare_particle_records(
            [{"global_particle_id": 10, "position_x": 1.0}, {"global_particle_id": 20, "position_x": 2.0}],
            [{"global_particle_id": 20, "position_x": 2.0}, {"global_particle_id": 10, "position_x": 1.0}],
            {"position_x": (0.0, 0.0)}, time=3.0,
        )
        self.assertTrue(result.passed)

    def test_canonical_contact_pair_and_alignment_are_order_independent(self) -> None:
        self.assertEqual(canonical_contact_pair(9, 3), (3, 9))
        left = [{"global_particle_id_a": 9, "global_particle_id_b": 3, "normal_force": 2.0}]
        right = [{"global_particle_id_a": 3, "global_particle_id_b": 9, "normal_force": 2.0}]
        self.assertEqual(len(align_by_contact_pair(left, right)), 1)
        result = compare_contact_records(left, right, {"normal_force": (0.0, 0.0)})
        self.assertTrue(result.passed)

    def test_contact_alignment_rejects_duplicate_canonical_pair(self) -> None:
        records = [
            {"global_particle_id_a": 1, "global_particle_id_b": 2},
            {"global_particle_id_a": 2, "global_particle_id_b": 1},
        ]
        with self.assertRaisesRegex(ValueError, "duplicate contact pair"):
            align_by_contact_pair(records, records[:1])

    def test_invalid_tolerances_are_rejected(self) -> None:
        with self.assertRaises(ValueError):
            make_tolerance(-1.0, 0.0)
        with self.assertRaises(ValueError):
            make_tolerance(0.0, float("inf"))


if __name__ == "__main__":
    unittest.main()
