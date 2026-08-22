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
    write_case_files,
)
from validation.jkr_tools.config import (
    default_material_config,
    default_numerical_config,
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


if __name__ == "__main__":
    unittest.main()
