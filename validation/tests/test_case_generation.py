import tempfile
import unittest
from pathlib import Path

from validation.jkr_tools.case_generation import (
    case_contact_time,
    keyword_count,
    mfix_dat_text,
    pair_overlaps,
    particle_data_rows,
    particle_input_text,
    validate_generated_case,
    wall_gap,
    write_case_files,
)
from validation.jkr_tools.config import (
    default_material_config,
    dmp_case_configs,
    elastic_jkr_material_config,
    default_numerical_config,
    elastic_jkr_material_config,
    heterogeneous_surface_energy_case_configs,
    hertz_case_configs,
    jkr_normal_case_configs,
    material_config,
    restart_case_config,
)


class RestartCaseGenerationTests(unittest.TestCase):
    def setUp(self):
        self.material = default_material_config()
        self.numerical = default_numerical_config()

    def generated_text(self, case_name):
        case = restart_case_config(case_name, self.material)
        mfix_text = mfix_dat_text(
            case=case,
            material=self.material,
            numerical=self.numerical,
            run_name="TEST_CASE",
            run_type="new",
            stop_time=case.stop_time,
            partition=self.numerical.serial_partition,
        )
        particle_text = particle_input_text(case.particles)
        return case, mfix_text, particle_text

    def test_restart_cases_generate_valid_inputs(self):
        for case_name in ("restart_pair_pfr", "restart_wall_oblique"):
            with self.subTest(case_name=case_name):
                case, mfix_text, particle_text = self.generated_text(case_name)
                errors = validate_generated_case(
                    mfix_text=mfix_text,
                    particle_text=particle_text,
                    case=case,
                    material=self.material,
                    numerical=self.numerical,
                )
                self.assertEqual(errors, ())

    def test_generated_header_is_mfix_26(self):
        _, mfix_text, _ = self.generated_text("restart_wall_oblique")
        self.assertIn("MFIX 26.1.2 JKR validation input", mfix_text)
        self.assertNotIn("MFIX 23.4.1", mfix_text)

    def test_required_keywords_are_unique(self):
        _, mfix_text, _ = self.generated_text("restart_pair_pfr")
        for keyword in (
            "run_name",
            "run_type",
            "dt",
            "particles",
            "use_jkr",
        ):
            self.assertEqual(keyword_count(mfix_text, keyword), 1)

    def test_particle_count_and_user_columns_match(self):
        case, _, particle_text = self.generated_text("restart_wall_oblique")
        rows = particle_data_rows(particle_text)
        self.assertEqual(len(rows), 1)
        self.assertEqual(len(rows[0]), 9 + case.des_usr_var_size)

    def test_pair_initial_overlap_is_intentional(self):
        case = restart_case_config("restart_pair_pfr", self.material)
        self.assertAlmostEqual(
            max(pair_overlaps(case.particles)),
            case.expected_initial_overlap,
            delta=1.0e-15,
        )

    def test_pair_contact_time_preserves_current_case_scale(self):
        case = restart_case_config("restart_pair_pfr", self.material)
        self.assertAlmostEqual(
            case_contact_time(case, self.material),
            6.136935272935e-8,
            delta=1.0e-20,
        )

    def test_duplicate_required_keyword_is_rejected(self):
        case, mfix_text, particle_text = self.generated_text(
            "restart_pair_pfr"
        )
        errors = validate_generated_case(
            mfix_text=mfix_text + "\nuse_jkr = .False.\n",
            particle_text=particle_text,
            case=case,
            material=self.material,
            numerical=self.numerical,
        )
        self.assertTrue(any("use_jkr occurs 2 times" in error for error in errors))

    def test_write_case_files_refuses_existing_directory(self):
        case = restart_case_config("restart_pair_pfr", self.material)
        with tempfile.TemporaryDirectory() as temporary_directory:
            output_directory = Path(temporary_directory) / "generated"
            write_case_files(
                output_directory=output_directory,
                case=case,
                material=self.material,
                numerical=self.numerical,
                run_name="TEST_PAIR",
                run_type="new",
                stop_time=case.stop_time,
                partition=self.numerical.serial_partition,
            )
            self.assertTrue((output_directory / "mfix.dat").is_file())
            with self.assertRaises(FileExistsError):
                write_case_files(
                    output_directory=output_directory,
                    case=case,
                    material=self.material,
                    numerical=self.numerical,
                    run_name="TEST_PAIR",
                    run_type="new",
                    stop_time=case.stop_time,
                    partition=self.numerical.serial_partition,
                )


class NativeHertzCaseGenerationTests(unittest.TestCase):
    def setUp(self):
        self.material = default_material_config()
        self.numerical = default_numerical_config()

    def generated_text(self, case):
        mfix_text = mfix_dat_text(
            case=case,
            material=self.material,
            numerical=self.numerical,
            run_name=case.name.upper(),
            run_type="new",
            stop_time=case.stop_time,
            partition=self.numerical.serial_partition,
            include_jkr_keywords=False,
        )
        return mfix_text, particle_input_text(case.particles)

    def test_official_compatible_inputs_omit_all_jkr_extensions(self):
        extension_keywords = (
            "use_jkr",
            "des_surface_energy",
            "des_surface_energy_wall",
            "jkr_alpha_safety",
            "jkr_thetar_crit",
            "wall_jkr_thetar_crit",
        )
        for case in hertz_case_configs(
            self.material, self.numerical.normal_initial_gap_fraction
        ):
            with self.subTest(case_name=case.name):
                mfix_text, particle_text = self.generated_text(case)
                errors = validate_generated_case(
                    mfix_text,
                    particle_text,
                    case,
                    self.material,
                    self.numerical,
                    include_jkr_keywords=False,
                )
                self.assertEqual(errors, ())
                for keyword in extension_keywords:
                    self.assertEqual(keyword_count(mfix_text, keyword), 0)

    def test_native_hertz_cases_enable_vtp_and_have_valid_geometry(self):
        for case in hertz_case_configs(
            self.material, self.numerical.normal_initial_gap_fraction
        ):
            with self.subTest(case_name=case.name):
                mfix_text, _ = self.generated_text(case)
                self.assertFalse(case.use_jkr)
                self.assertIn("write_vtk_files = .True.", mfix_text)
                if case.wall_contact:
                    self.assertGreater(case.expected_initial_overlap, 0.0)
                else:
                    self.assertLess(case.expected_initial_overlap, 0.0)

    def test_jkr_case_cannot_omit_extension_keywords(self):
        case = restart_case_config("restart_pair_pfr", self.material)
        mfix_text = mfix_dat_text(
            case,
            self.material,
            self.numerical,
            "INVALID_NATIVE",
            "new",
            case.stop_time,
            self.numerical.serial_partition,
            include_jkr_keywords=False,
        )
        errors = validate_generated_case(
            mfix_text,
            particle_input_text(case.particles),
            case,
            self.material,
            self.numerical,
            include_jkr_keywords=False,
        )
        self.assertIn("cannot omit JKR keywords for a JKR-enabled case", errors)


class JkrNormalCaseGenerationTests(unittest.TestCase):
    def setUp(self):
        self.material = default_material_config()
        self.numerical = default_numerical_config()
        self.wall_critical_speed = 2.0
        self.pair_critical_speed = 4.0
        self.cases = jkr_normal_case_configs(
            self.material,
            self.wall_critical_speed,
            self.pair_critical_speed,
            self.numerical.normal_initial_gap_fraction,
        )

    def test_case_names_and_speed_factors_are_explicit(self):
        expected = (
            ("wall_normal_damped_capture", 0.95 * self.wall_critical_speed),
            ("wall_normal_damped_rebound", 3.0 * self.wall_critical_speed),
            ("pair_normal_damped_capture", 0.95 * self.pair_critical_speed),
            ("pair_normal_damped_rebound", 3.0 * self.pair_critical_speed),
        )
        self.assertEqual(
            tuple(case.name for case in self.cases),
            tuple(item[0] for item in expected),
        )
        for case, (_, expected_speed) in zip(self.cases, expected):
            with self.subTest(case_name=case.name):
                self.assertTrue(case.use_jkr)
                self.assertAlmostEqual(case.contact_time_speed, expected_speed)

    def test_jkr_inputs_validate_and_contain_required_extensions(self):
        for case in self.cases:
            with self.subTest(case_name=case.name):
                mfix_text = mfix_dat_text(
                    case=case,
                    material=self.material,
                    numerical=self.numerical,
                    run_name=case.name.upper(),
                    run_type="new",
                    stop_time=case.stop_time,
                    partition=self.numerical.serial_partition,
                )
                particle_text = particle_input_text(case.particles)
                self.assertEqual(
                    validate_generated_case(
                        mfix_text=mfix_text,
                        particle_text=particle_text,
                        case=case,
                        material=self.material,
                        numerical=self.numerical,
                    ),
                    (),
                )
                self.assertEqual(keyword_count(mfix_text, "use_jkr"), 1)
                self.assertIn("use_jkr = .True.", mfix_text)

    def test_wall_and_pair_start_with_gap_not_unintended_overlap(self):
        expected_gap = (
            self.numerical.normal_initial_gap_fraction
            * self.material.particle_radius
        )
        for case in self.cases:
            with self.subTest(case_name=case.name):
                if case.wall_contact:
                    self.assertEqual(len(case.particles), 1)
                    self.assertAlmostEqual(
                        wall_gap(case.particles[0]), expected_gap
                    )
                    self.assertGreater(wall_gap(case.particles[0]), 0.0)
                else:
                    overlaps = pair_overlaps(case.particles)
                    self.assertEqual(len(overlaps), 1)
                    self.assertAlmostEqual(overlaps[0], -expected_gap)
                    self.assertLess(overlaps[0], 0.0)

    def test_non_positive_critical_speeds_are_rejected(self):
        for wall_speed, pair_speed in (
            (0.0, self.pair_critical_speed),
            (-1.0, self.pair_critical_speed),
            (self.wall_critical_speed, 0.0),
            (self.wall_critical_speed, -1.0),
        ):
            with self.subTest(wall_speed=wall_speed, pair_speed=pair_speed):
                with self.assertRaisesRegex(
                    ValueError, "critical speeds must be positive"
                ):
                    jkr_normal_case_configs(
                        self.material,
                        wall_speed,
                        pair_speed,
                        self.numerical.normal_initial_gap_fraction,
                    )
class HeterogeneousSurfaceEnergyCaseGenerationTests(unittest.TestCase):
    def setUp(self):
        self.material = elastic_jkr_material_config()
        self.numerical = default_numerical_config()
        self.cases = heterogeneous_surface_energy_case_configs(
            self.material,
            wall_impact_speed=1.0,
            pair_impact_speed=1.0,
            initial_gap_fraction=self.numerical.normal_initial_gap_fraction,
        )
        self.by_name = {case.name: case for case in self.cases}

    def test_case_names_and_elastic_profile_are_explicit(self):
        self.assertEqual(
            tuple(case.name for case in self.cases),
            (
                "pair_gamma_homogeneous",
                "pair_gamma_low_high",
                "pair_gamma_high_low",
                "pair_gamma_id_swap",
                "pair_gamma_fallback",
                "pair_gamma_invalid_nonpositive",
                "wall_gamma_heterogeneous",
            ),
        )
        for case in self.cases:
            with self.subTest(case_name=case.name):
                material = material_config(case.material_profile)
                self.assertEqual(material.normal_restitution, 1.0)
                expected_profile = (
                    "elastic_jkr_wall_heterogeneous"
                    if case.wall_contact else "elastic_jkr"
                )
                self.assertEqual(case.material_profile, expected_profile)
                self.assertTrue(case.use_jkr)

    def test_all_approach_cases_use_five_percent_radius_gap(self):
        expected_gap = (
            self.numerical.normal_initial_gap_fraction
            * self.material.particle_radius
        )
        for case in self.cases:
            with self.subTest(case_name=case.name):
                if case.wall_contact:
                    self.assertAlmostEqual(
                        wall_gap(case.particles[0]), expected_gap
                    )
                else:
                    self.assertAlmostEqual(
                        max(pair_overlaps(case.particles)), -expected_gap
                    )

    def test_low_high_geometric_mean_and_high_low_symmetry(self):
        low_high = self.by_name["pair_gamma_low_high"]
        high_low = self.by_name["pair_gamma_high_low"]
        low_high_gamma = tuple(
            particle.surface_energy for particle in low_high.particles
        )
        high_low_gamma = tuple(
            particle.surface_energy for particle in high_low.particles
        )
        self.assertEqual(low_high_gamma, (0.0075, 0.030))
        self.assertEqual(high_low_gamma, (0.030, 0.0075))
        self.assertAlmostEqual(
            (low_high_gamma[0] * low_high_gamma[1]) ** 0.5,
            0.015,
        )
        self.assertEqual(sorted(low_high_gamma), sorted(high_low_gamma))

    def test_id_swap_only_reverses_row_identity(self):
        low_high = self.by_name["pair_gamma_low_high"].particles
        id_swap = self.by_name["pair_gamma_id_swap"].particles
        physical_state = lambda particles: sorted(
            (particle.position, particle.velocity, particle.surface_energy)
            for particle in particles
        )
        self.assertEqual(physical_state(low_high), physical_state(id_swap))
        self.assertEqual(id_swap, tuple(reversed(low_high)))

    def test_fallback_omits_user_scalar_and_invalid_values_are_preserved(self):
        fallback = self.by_name["pair_gamma_fallback"]
        self.assertEqual(fallback.des_usr_var_size, 0)
        fallback_text = particle_input_text(fallback.particles)
        self.assertNotIn("| usr_", fallback_text)
        self.assertTrue(
            all(len(row) == 9 for row in particle_data_rows(fallback_text))
        )

        invalid = self.by_name["pair_gamma_invalid_nonpositive"]
        invalid_rows = particle_data_rows(
            particle_input_text(invalid.particles)
        )
        self.assertEqual(tuple(row[9] for row in invalid_rows), (0.0, -0.001))

    def test_wall_effective_surface_energy_uses_wall_profile(self):
        case = self.by_name["wall_gamma_heterogeneous"]
        material = material_config(case.material_profile)
        self.assertAlmostEqual(material.surface_energy, 0.015)
        self.assertAlmostEqual(material.wall_surface_energy, 0.030)
        self.assertAlmostEqual(
            (0.0075 * material.wall_surface_energy) ** 0.5,
            0.015,
        )
class DmpCaseGenerationTests(unittest.TestCase):
    def setUp(self):
        self.elastic = elastic_jkr_material_config()
        self.damped = default_material_config()
        self.numerical = default_numerical_config()
        self.cases = dmp_case_configs(self.elastic, self.damped)

    def test_fixed_two_rank_case_registry_and_partitions(self):
        self.assertEqual(
            tuple(case.name for case in self.cases),
            (
                "pair_cross_jkr",
                "pair_new_contact_jkr",
                "pair_active_contact_migration",
                "wall_cross_jkr",
                "cloud_collision_jkr",
            ),
        )
        self.assertEqual(self.numerical.dmp_ranks, 2)
        self.assertEqual(self.numerical.dmp_partition, (2, 1, 1))
        for case in self.cases:
            with self.subTest(case_name=case.name):
                self.assertEqual(case.category, "dmp")
                self.assertTrue(case.use_jkr)
                self.assertEqual(case.restart_time, None)
                expected_profile = (
                    "elastic_jkr"
                    if case.name in {"pair_cross_jkr", "pair_new_contact_jkr"}
                    else "damped_frictional"
                )
                self.assertEqual(case.material_profile, expected_profile)

    def test_serial_and_dmp_inputs_share_particle_bytes_and_only_partition_differs(self):
        def without_partition(text):
            return tuple(
                line for line in text.splitlines()
                if not line.strip().startswith(("nodesi =", "nodesj =", "nodesk ="))
            )

        for case in self.cases:
            material = (
                self.elastic
                if case.material_profile == "elastic_jkr"
                else self.damped
            )
            with self.subTest(case_name=case.name):
                serial_mfix = mfix_dat_text(
                    case,
                    material,
                    self.numerical,
                    "DMP_CASE",
                    "new",
                    case.stop_time,
                    self.numerical.serial_partition,
                )
                dmp_mfix = mfix_dat_text(
                    case,
                    material,
                    self.numerical,
                    "DMP_CASE",
                    "new",
                    case.stop_time,
                    self.numerical.dmp_partition,
                )
                self.assertEqual(
                    without_partition(serial_mfix),
                    without_partition(dmp_mfix),
                )
                self.assertNotEqual(serial_mfix, dmp_mfix)
                serial_particle = particle_input_text(case.particles)
                dmp_particle = particle_input_text(case.particles)


if __name__ == "__main__":
    unittest.main()
