"""Tests for the independent e=1 trajectory and contact figures."""

import importlib.util
import tempfile
import unittest
from pathlib import Path

from validation.jkr_tools.config import elastic_jkr_material_config
from validation.jkr_tools.plotting import plot_contact_validation
from validation.jkr_tools.reference import (
    adhesion_only_critical_speed,
    simulate_normal_contact,
)
from validation.run_validation import normal_contact_parameters


class ElasticReferenceTrajectoryTests(unittest.TestCase):
    def test_adhesion_only_threshold_brackets_capture_and_rebound(self):
        material = elastic_jkr_material_config()
        for wall_contact in (True, False):
            parameters = normal_contact_parameters(material, wall_contact)
            critical_speed = adhesion_only_critical_speed(parameters)
            with self.subTest(contact="wall" if wall_contact else "pair"):
                capture = simulate_normal_contact(
                    (0.0, 2.0e-6),
                    -0.2e-6,
                    0.95 * critical_speed,
                    parameters,
                    2.0e-11,
                )
                rebound = simulate_normal_contact(
                    (0.0, 2.0e-6),
                    -0.2e-6,
                    1.05 * critical_speed,
                    parameters,
                    2.0e-11,
                )
                self.assertTrue(capture[-1]["active"])
                self.assertFalse(rebound[-1]["active"])
                self.assertIsNone(rebound[-1]["contact_radius"])


@unittest.skipUnless(
    importlib.util.find_spec("matplotlib"), "matplotlib is not installed"
)
class ContactPlotTests(unittest.TestCase):
    def test_four_standard_contact_plots_are_written(self):
        mfix_rows = (
            {
                "time": 0.0,
                "normal_speed": -0.1,
                "angular_speed_1": 0.0,
                "angular_speed_2": None,
                "overlap": -1.0e-7,
                "observed_normal_force": 0.0,
                "contact_radius": None,
                "active": False,
            },
            {
                "time": 1.0e-7,
                "normal_speed": 0.0,
                "angular_speed_1": 0.0,
                "angular_speed_2": None,
                "overlap": 1.0e-9,
                "observed_normal_force": 1.0e-8,
                "contact_radius": 1.0e-7,
                "active": True,
            },
        )
        reference_rows = (
            {
                "time": 0.0,
                "normal_speed": -0.1,
                "overlap": -1.0e-7,
                "normal_force": 0.0,
                "contact_radius": None,
                "active": False,
            },
            {
                "time": 1.0e-7,
                "normal_speed": 0.0,
                "overlap": 1.0e-9,
                "normal_force": 1.0e-8,
                "contact_radius": 1.0e-7,
                "active": True,
            },
        )
        with tempfile.TemporaryDirectory() as temporary:
            paths = plot_contact_validation(
                mfix_rows,
                reference_rows,
                Path(temporary),
                "contact",
            )
            self.assertEqual(len(paths), 4)
            self.assertTrue(all(path.stat().st_size > 0 for path in paths))


if __name__ == "__main__":
    unittest.main()
