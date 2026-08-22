"""Tests for independent normal JKR trajectory comparison."""

from __future__ import annotations

import unittest
from pathlib import Path
from unittest.mock import patch

from validation.jkr_tools.config import default_material_config
from validation.jkr_tools.jkr import analyze_normal_jkr_trajectory
from validation.jkr_tools.reference import (
    contact_parameters,
    elastic_normal_force,
    normal_damping_coefficient,
    stable_contact_radius,
)


def _parameters():
    material = default_material_config()
    return contact_parameters(
        material.particle_radius,
        material.particle_density,
        material.surface_energy,
        material.young_modulus,
        material.poisson_ratio,
        material.normal_restitution,
        None,
        None,
        material.wall_surface_energy,
        material.wall_young_modulus,
        material.wall_poisson_ratio,
    )


def _frame(time, overlap, speed, force):
    material = default_material_config()
    return {
        "time": time,
        "records": ({
            "Particle_ID": 7,
            "Coordinates": (10.0e-6, material.particle_radius - overlap, 10.0e-6),
            "Velocity": (0.0, speed, 0.0),
            "Collision_force": (0.0, force, 0.0),
        },),
    }


def _force(overlap, speed, parameters):
    radius = stable_contact_radius(overlap, parameters)
    return (
        elastic_normal_force(radius, parameters)
        - normal_damping_coefficient(radius, parameters) * speed
    )


class JkrTrajectoryTests(unittest.TestCase):
    def test_capture_uses_contact_history_through_negative_overlap(self):
        parameters = _parameters()
        overlaps = (-1.0e-7, 1.0e-9, -0.5 * parameters.critical_separation)
        speeds = (-0.1, -0.02, 0.0)
        frames = (
            _frame(0.0, overlaps[0], speeds[0], 0.0),
            _frame(
                1.0, overlaps[1], speeds[1],
                _force(overlaps[1], speeds[1], parameters),
            ),
            _frame(
                2.0, overlaps[2], speeds[2],
                _force(overlaps[2], speeds[2], parameters),
            ),
        )
        with patch(
            "validation.jkr_tools.jkr.read_vtp_trajectory",
            return_value=frames,
        ):
            result = analyze_normal_jkr_trajectory(
                (Path("synthetic.vtp"),),
                "wall",
                parameters,
                default_material_config().particle_radius,
                1.0e-15,
                1.0e-10,
                1.0e-8,
            )
        self.assertTrue(result["passed"])
        self.assertTrue(result["captured"])
        self.assertTrue(result["final_contact_active"])
        self.assertAlmostEqual(result["minimum_overlap"], -1.0e-7)

    def test_rebound_requires_break_beyond_critical_separation(self):
        parameters = _parameters()
        positive_overlap = 1.0e-9
        speed = -0.02
        frames = (
            _frame(0.0, -1.0e-7, -0.1, 0.0),
            _frame(
                1.0, positive_overlap, speed,
                _force(positive_overlap, speed, parameters),
            ),
            _frame(2.0, -1.01 * parameters.critical_separation, 0.03, 0.0),
            _frame(3.0, -2.0 * parameters.critical_separation, 0.03, 0.0),
        )
        with patch(
            "validation.jkr_tools.jkr.read_vtp_trajectory",
            return_value=frames,
        ):
            result = analyze_normal_jkr_trajectory(
                (Path("synthetic.vtp"),),
                "wall",
                parameters,
                default_material_config().particle_radius,
                1.0e-15,
                1.0e-10,
                1.0e-8,
            )
        self.assertTrue(result["passed"])
        self.assertFalse(result["captured"])
        self.assertFalse(result["final_contact_active"])
        self.assertAlmostEqual(result["rebound_velocity"], 0.03)

    def test_force_difference_reports_first_failure(self):
        parameters = _parameters()
        frames = (
            _frame(0.0, 1.0e-9, -0.02, 0.0),
            _frame(1.0, 2.0e-9, 0.01, 0.0),
        )
        with patch(
            "validation.jkr_tools.jkr.read_vtp_trajectory",
            return_value=frames,
        ):
            result = analyze_normal_jkr_trajectory(
                (Path("synthetic.vtp"),),
                "wall",
                parameters,
                default_material_config().particle_radius,
                0.0,
                0.0,
                0.0,
            )
        self.assertFalse(result["passed"])
        self.assertEqual(result["force_comparison"]["first_difference"]["time"], 0.0)


if __name__ == "__main__":
    unittest.main()
