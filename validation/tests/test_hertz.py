"""Synthetic tests for Hertz trajectory utilities."""
from __future__ import annotations
import unittest
from validation.jkr_tools.hertz import compare_trajectories, extract_contact_metrics, trajectory_from_records, frame_to_records
from validation.jkr_tools.io import VtpArray, VtpData

def _frame(time, rows):
    return {"time": time, "records": tuple(rows)}

def _row(pid, x, velocity, diameter=2.0):
    return {"Particle_ID": pid, "Coordinates": (x, 0.0, 0.0), "Velocity": (velocity, 0.0, 0.0), "Diameter": diameter, "Angular_velocity": (0.0, 0.0, 0.0), "Collision_force": (0.0, 0.0, 0.0)}

class HertzTrajectoryTests(unittest.TestCase):
    def test_compare_aligns_by_global_id(self):
        left = trajectory_from_records((_frame(0.0, (_row(10, 0.0, 1.0), _row(20, 2.0, -1.0))),))
        right = trajectory_from_records((_frame(0.0, (_row(20, 2.0, -1.0), _row(10, 0.0, 1.0))),))
        result = compare_trajectories(left, right)
        self.assertEqual(set(result), {"Coordinates", "Velocity", "Angular_velocity", "Collision_force"})
        self.assertTrue(all(item.passed for item in result.values()))

    def test_duplicate_missing_ids_and_time_mismatch(self):
        with self.assertRaisesRegex(ValueError, "duplicate"):
            trajectory_from_records((_frame(0.0, (_row(1, 0, 0), _row(1, 1, 0))),))
        left = trajectory_from_records((_frame(0.0, (_row(1, 0, 0),)),))
        with self.assertRaisesRegex(ValueError, "particle IDs differ"):
            compare_trajectories(left, trajectory_from_records((_frame(0.0, (_row(2, 0, 0),)),)))
        with self.assertRaisesRegex(ValueError, "time mismatch"):
            compare_trajectories(left, trajectory_from_records((_frame(1.0, (_row(1, 0, 0),)),)))

    def test_identical_terminal_frame_is_normalized(self):
        frame = _frame(0.0, (_row(1, 0, 0),))
        self.assertEqual(len(trajectory_from_records((frame, frame))), 1)

    def test_pair_metrics(self):
        frames = trajectory_from_records((
            _frame(0.0, (_row(10, 0.0, 1.0), _row(20, 2.2, -1.0))),
            _frame(1.0, (_row(10, 0.4, 0.5), _row(20, 1.6, -0.5))),
            _frame(2.0, (_row(10, 0.9, -0.5), _row(20, 1.1, 0.5))),
            _frame(3.0, (_row(10, 0.0, -0.5), _row(20, 3.0, 0.5))),
        ))
        metrics = extract_contact_metrics(frames)
        self.assertAlmostEqual(metrics["maximum_overlap"], 1.8)
        self.assertEqual(metrics["contact_duration"], 1.0)
        self.assertEqual(metrics["impact_velocity"], 2.0)
        self.assertEqual(metrics["rebound_normal_velocity"], 1.0)
        self.assertEqual(metrics["effective_restitution"], 0.5)

    def test_synthetic_vtp_decoding(self):
        arrays = {
            "Particle_ID": VtpArray(1, (20.0, 10.0), "Float64"),
            "coordinates": VtpArray(3, (2.0, 0.0, 0.0, 0.0, 0.0, 0.0), "Float64"),
            "Velocity": VtpArray(3, (-1.0, 0.0, 0.0, 1.0, 0.0, 0.0), "Float64"),
        }
        records = frame_to_records(VtpData(0.0, arrays, "LittleEndian", "UInt32"))
        self.assertEqual([row["Particle_ID"] for row in records["records"]], [20, 10])
        self.assertEqual(records["records"][0]["Coordinates"], (2.0, 0.0, 0.0))

if __name__ == "__main__":
    unittest.main()
