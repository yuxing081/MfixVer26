"""Unit tests for DMP trajectory identity and contact-set helpers."""

import unittest

from validation.jkr_tools.comparison import canonical_contact_pair
from validation.jkr_tools.dmp import (
    compare_contact_sets,
    contact_sets,
    initial_physical_id_mapping,
    rank_transitions,
    remap_to_serial_ids,
)


def frame(time, records):
    return {"time": time, "records": tuple(records)}


def particle(
    particle_id,
    coordinates,
    velocity=(0.0, 0.0, 0.0),
    diameter=2.0,
    rank=0,
    force=(0.0, 0.0, 0.0),
):
    return {
        "Particle_ID": particle_id,
        "Coordinates": coordinates,
        "Velocity": velocity,
        "Diameter": diameter,
        "Particle_Rank": rank,
        "Collision_force": force,
    }


class DmpIdentityTests(unittest.TestCase):
    def test_initial_mapping_uses_physical_state_and_remap_preserves_serial_ids(self):
        serial = [
            frame(0.0, (
                particle(1, (0.0, 0.0, 0.0), velocity=(1.0, 0.0, 0.0)),
                particle(2, (4.0, 0.0, 0.0), velocity=(-1.0, 0.0, 0.0)),
            ))
        ]
        dmp = [
            frame(0.0, (
                particle(20, (4.0, 0.0, 0.0), velocity=(-1.0, 0.0, 0.0)),
                particle(10, (0.0, 0.0, 0.0), velocity=(1.0, 0.0, 0.0)),
            )),
            frame(1.0, (
                particle(20, (3.0, 0.0, 0.0), velocity=(-1.0, 0.0, 0.0)),
                particle(10, (1.0, 0.0, 0.0), velocity=(1.0, 0.0, 0.0)),
            )),
        ]
        mapping = initial_physical_id_mapping(serial, dmp, 1.0e-12, 1.0e-12)
        self.assertEqual(mapping, {1: 10, 2: 20})
        remapped = remap_to_serial_ids(dmp, mapping)
        self.assertEqual(
            tuple(record["Particle_ID"] for record in remapped[0]["records"]),
            (2, 1),
        )

    def test_canonical_pair_is_order_independent(self):
        self.assertEqual(canonical_contact_pair(20, 10), (10, 20))


class DmpTransitionAndContactTests(unittest.TestCase):
    def test_rank_transitions_report_only_owner_changes(self):
        trajectory = (
            frame(0.0, (
                particle(1, (0.0, 0.0, 0.0), rank=0),
                particle(2, (4.0, 0.0, 0.0), rank=0),
            )),
            frame(1.0, (
                particle(1, (1.0, 0.0, 0.0), rank=1),
                particle(2, (3.0, 0.0, 0.0), rank=0),
            )),
            frame(2.0, (
                particle(1, (2.0, 0.0, 0.0), rank=1),
                particle(2, (2.0, 0.0, 0.0), rank=1),
            )),
        )
        self.assertEqual(
            rank_transitions(trajectory),
            {
                1: ({"time": 0.0, "rank": 0}, {"time": 1.0, "rank": 1}),
                2: ({"time": 0.0, "rank": 0}, {"time": 2.0, "rank": 1}),
            },
        )

    def test_contact_set_comparison_reports_first_difference(self):
        serial = (
            {"time": 0.0, "contacts": ()},
            {"time": 1.0, "contacts": ((1, 2),)},
            {"time": 2.0, "contacts": ((1, 2),)},
        )
        dmp = (
            {"time": 0.0, "contacts": ()},
            {"time": 1.0, "contacts": ()},
            {"time": 2.0, "contacts": ((1, 2),)},
        )
        result = compare_contact_sets(serial, dmp)
        self.assertFalse(result["passed"])
        self.assertEqual(result["first_difference"]["time"], 1.0)
        self.assertEqual(result["first_difference"]["serial_contacts"], ((1, 2),))
        self.assertEqual(result["first_active_time"], 1.0)

    def test_pair_contact_inference_uses_geometry_and_force(self):
        trajectory = (
            frame(0.0, (
                particle(1, (0.0, 0.0, 0.0)),
                particle(2, (2.1, 0.0, 0.0)),
            )),
            frame(1.0, (
                particle(1, (0.0, 0.0, 0.0), force=(1.0, 0.0, 0.0)),
                particle(2, (2.1, 0.0, 0.0), force=(-1.0, 0.0, 0.0)),
            )),
        )
        contacts = contact_sets(
            trajectory,
            wall_contact=False,
            force_tolerance=0.1,
            overlap_tolerance=1.0e-12,
        )
        self.assertEqual(contacts[0]["contacts"], ())
        self.assertEqual(contacts[1]["contacts"], ((1, 2),))


if __name__ == "__main__":
    unittest.main()
