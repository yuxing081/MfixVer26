"""Tests for the replaceable, category-based result snapshot."""

import tempfile
import unittest
from pathlib import Path

from validation.jkr_tools.config import (
    default_material_config,
    restart_case_config,
)
from validation.run_validation import (
    case_result_root,
    clear_results_directory,
    hertz_validation_cases,
    serial_jkr_validation_cases,
    write_run_info,
)


class ResultPathTests(unittest.TestCase):
    def test_contact_cases_are_grouped_by_physical_subcategory(self):
        root = Path("/validation/results")
        cases = {
            case.name: case
            for case in (
                *hertz_validation_cases(),
                *serial_jkr_validation_cases(),
            )
        }
        self.assertEqual(
            case_result_root(root, cases["pair_normal_capture"]),
            root / "contact_physics" / "pair_normal" / "pair_normal_capture",
        )
        self.assertEqual(
            case_result_root(root, cases["wall_normal_rebound"]),
            root / "contact_physics" / "wall_normal" / "wall_normal_rebound",
        )
        self.assertEqual(
            case_result_root(root, cases["wall_hertz_oblique"]),
            root / "contact_physics" / "wall_oblique" / "wall_hertz_oblique",
        )

    def test_restart_cases_are_grouped_by_history_kind(self):
        root = Path("/validation/results")
        material = default_material_config()
        pair = restart_case_config("restart_pair_pfr", material)
        wall = restart_case_config("restart_wall_oblique", material)
        self.assertEqual(
            case_result_root(root, pair),
            root / "restart" / "pair_history" / pair.name,
        )
        self.assertEqual(
            case_result_root(root, wall),
            root / "restart" / "wall_history" / wall.name,
        )

    def test_all_current_case_roots_are_unique(self):
        root = Path("/validation/results")
        material = default_material_config()
        cases = (
            *hertz_validation_cases(),
            *serial_jkr_validation_cases(),
            restart_case_config("restart_pair_pfr", material),
            restart_case_config("restart_wall_oblique", material),
        )
        paths = tuple(case_result_root(root, case) for case in cases)
        self.assertEqual(len(paths), len(set(paths)))


class ResultSnapshotTests(unittest.TestCase):
    def test_clear_preserves_gitkeep_and_removes_previous_snapshot(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / ".gitkeep").write_text("", encoding="utf-8")
            old_case = root / "contact_physics" / "pair_normal" / "old"
            old_case.mkdir(parents=True)
            (old_case / "solver.log").write_text("old", encoding="utf-8")
            (root / "RUN_INFO.txt").write_text("old", encoding="utf-8")

            clear_results_directory(root)

            self.assertEqual(tuple(path.name for path in root.iterdir()), (".gitkeep",))

    def test_run_info_records_version_commit_solver_and_process_count(self):
        metadata = {
            "mfix_version": "MFIX 26.1.2",
            "jkr_git": {
                "commit": "abc123",
                "branch": "codex/jkr-validation",
                "dirty": False,
            },
            "solver": {
                "path": "/solver",
                "sha256": "deadbeef",
            },
            "execution": {
                "mode": "serial",
                "mpi_processes": 1,
                "partition": (1, 1, 1),
            },
        }
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "RUN_INFO.txt"
            write_run_info(path, "all", metadata)
            text = path.read_text(encoding="utf-8")
        self.assertIn("mfix_version=MFIX 26.1.2", text)
        self.assertIn("git_commit=abc123", text)
        self.assertIn("solver_sha256=deadbeef", text)
        self.assertIn("mpi_processes=1", text)


if __name__ == "__main__":
    unittest.main()
