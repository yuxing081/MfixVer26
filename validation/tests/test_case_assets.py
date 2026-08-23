"""Tests for explicit material profiles and committed case assets."""

import tempfile
import unittest
from pathlib import Path

from validation.jkr_tools.case_generation import (
    canonical_case_sync_errors,
    copy_canonical_case_asset,
    pair_overlaps,
    wall_gap,
)
from validation.jkr_tools.config import (
    damped_frictional_material_config,
    default_numerical_config,
    elastic_jkr_material_config,
    material_config,
)
from validation.run_validation import (
    CASES_ROOT,
    canonical_contact_cases,
    normal_jkr_cases,
)


class MaterialProfileTests(unittest.TestCase):
    def test_profiles_split_normal_restitution_explicitly(self):
        self.assertEqual(
            elastic_jkr_material_config().normal_restitution, 1.0
        )
        self.assertEqual(
            damped_frictional_material_config().normal_restitution, 0.6
        )
        with self.assertRaisesRegex(ValueError, "Unknown material profile"):
            material_config("implicit_default")

    def test_elastic_and_damped_normal_case_names_and_factors(self):
        elastic = normal_jkr_cases(elastic_jkr_material_config())
        damped = normal_jkr_cases(damped_frictional_material_config())
        self.assertEqual(
            tuple(case.name for case in elastic),
            (
                "wall_normal_capture",
                "wall_normal_rebound",
                "pair_normal_capture",
                "pair_normal_rebound",
            ),
        )
        self.assertEqual(
            tuple(case.critical_speed_factor for case in elastic),
            (0.95, 1.05, 0.95, 1.05),
        )
        self.assertEqual(
            tuple(case.name for case in damped),
            (
                "wall_normal_damped_capture",
                "wall_normal_damped_rebound",
                "pair_normal_damped_capture",
                "pair_normal_damped_rebound",
            ),
        )
        self.assertEqual(
            tuple(case.critical_speed_factor for case in damped),
            (0.95, 3.0, 0.95, 3.0),
        )

    def test_heterogeneous_wall_profile_is_elastic_with_unlike_wall_energy(self):
        material = material_config("elastic_jkr_wall_heterogeneous")
        self.assertEqual(material.normal_restitution, 1.0)
        self.assertEqual(material.tangential_restitution, 1.0)
        self.assertEqual(material.surface_energy, 0.015)
        self.assertEqual(material.wall_surface_energy, 0.030)


class CanonicalCaseAssetTests(unittest.TestCase):
    def test_all_committed_assets_match_generator_byte_for_byte(self):
        numerical = default_numerical_config()
        cases = canonical_contact_cases()
        self.assertEqual(len(cases), 29)
        for case in cases:
            with self.subTest(case=case.name):
                self.assertEqual(
                    canonical_case_sync_errors(
                        CASES_ROOT,
                        case,
                        material_config(case.material_profile),
                        numerical,
                    ),
                    (),
                )

    def test_runtime_copy_preserves_both_asset_hash_inputs(self):
        case = canonical_contact_cases()[0]
        material = material_config(case.material_profile)
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / case.name
            copy_canonical_case_asset(
                CASES_ROOT,
                output,
                case,
                material,
                default_numerical_config(),
            )
            source = CASES_ROOT / case.category / case.name
            for name in ("mfix.dat", "particle_input.dat"):
                self.assertEqual(
                    (source / name).read_bytes(),
                    (output / name).read_bytes(),
                )

    def test_standard_approach_cases_use_five_percent_radius_gap(self):
        numerical = default_numerical_config()
        for case in canonical_contact_cases():
            if case.category in {"restart", "dmp"}:
                continue
            material = material_config(case.material_profile)
            expected_gap = (
                numerical.normal_initial_gap_fraction
                * material.particle_radius
            )
            with self.subTest(case=case.name):
                if case.wall_contact:
                    self.assertAlmostEqual(
                        wall_gap(case.particles[0]), expected_gap
                    )
                else:
                    overlaps = pair_overlaps(case.particles)
                    self.assertEqual(len(overlaps), 1)
                    self.assertAlmostEqual(overlaps[0], -expected_gap)


if __name__ == "__main__":
    unittest.main()
