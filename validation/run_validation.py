#!/usr/bin/env python3
"""Run reproducible MFIX 26.1.2 JKR validation suites."""

from __future__ import annotations

import argparse
import csv
import dataclasses
import datetime
import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

from validation.jkr_tools.case_generation import write_case_files
from validation.jkr_tools.config import (
    default_material_config,
    default_numerical_config,
    restart_case_configs,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SOLVER = PROJECT_ROOT / "code" / "build" / "mfixsolver"
RESULTS_ROOT = PROJECT_ROOT / "validation" / "results" / "runs"


def command_text(command: list[str]) -> str:
    completed = subprocess.run(
        command,
        cwd=PROJECT_ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    return (completed.stdout + completed.stderr).strip()


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def git_metadata() -> dict[str, Any]:
    status = command_text(["git", "status", "--porcelain=v1"])
    return {
        "commit": command_text(["git", "rev-parse", "HEAD"]),
        "branch": command_text(["git", "branch", "--show-current"]),
        "dirty": bool(status),
        "status_porcelain": status.splitlines(),
    }


def compiler_metadata(solver: Path) -> dict[str, Any]:
    cache_path = solver.parent / "CMakeCache.txt"
    compiler_path = None
    if cache_path.is_file():
        prefix = "CMAKE_Fortran_COMPILER:FILEPATH="
        for line in cache_path.read_text(encoding="utf-8").splitlines():
            if line.startswith(prefix):
                compiler_path = Path(line.removeprefix(prefix))
                break
    if compiler_path is None or not compiler_path.is_file():
        return {
            "cmake_cache": str(cache_path),
            "fortran_compiler": None,
            "version": None,
        }
    return {
        "cmake_cache": str(cache_path),
        "fortran_compiler": str(compiler_path),
        "version": command_text([str(compiler_path), "--version"]).splitlines()[:1],
    }


def environment_metadata(solver: Path) -> dict[str, Any]:
    material = default_material_config()
    numerical = default_numerical_config()
    return {
        "mfix_version": command_text([str(solver), "--version"]),
        "jkr_git": git_metadata(),
        "solver": {
            "path": str(solver),
            "sha256": file_sha256(solver),
        },
        "execution": {
            "mode": "serial",
            "mpi_processes": 1,
            "partition": numerical.serial_partition,
        },
        "compiler": compiler_metadata(solver),
        "python": sys.version,
        "material": dataclasses.asdict(material),
        "numerical": dataclasses.asdict(numerical),
    }


def timestamped_run_directory() -> Path:
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime(
        "%Y%m%dT%H%M%S.%fZ"
    )
    directory = RESULTS_ROOT / f"{stamp}-{os.getpid()}"
    directory.mkdir(parents=True, exist_ok=False)
    return directory


def run_solver(solver: Path, case_directory: Path) -> tuple[bool, str]:
    log_path = case_directory / "solver.log"
    with log_path.open("w", encoding="utf-8") as log:
        completed = subprocess.run(
            [str(solver)],
            cwd=case_directory,
            check=False,
            stdout=log,
            stderr=subprocess.STDOUT,
            text=True,
        )
    return completed.returncode == 0, f"solver exit code {completed.returncode}"


def restart_artifacts(stage_directory: Path, run_name: str) -> tuple[Path, ...]:
    candidates = [
        stage_directory / f"{run_name}.RES",
        stage_directory / f"{run_name}.TXT",
        stage_directory / f"{run_name}_DES.RES",
    ]
    candidates.extend(sorted(stage_directory.glob(f"{run_name}.SP*")))
    return tuple(path for path in candidates if path.is_file())


def copy_restart_artifacts(
    stage_directory: Path,
    branch_directory: Path,
    run_name: str,
) -> None:
    artifacts = restart_artifacts(stage_directory, run_name)
    if not any(path.name == f"{run_name}.RES" for path in artifacts):
        raise FileNotFoundError(f"missing restart RES file for {run_name}")
    for source in artifacts:
        shutil.copy2(source, branch_directory / source.name)


def remove_restart_2_history(branch_directory: Path, run_name: str) -> None:
    for path in branch_directory.glob(f"{run_name}.SP*"):
        path.unlink()
    text_path = branch_directory / f"{run_name}.TXT"
    if text_path.exists():
        text_path.unlink()


def compare_particle_output(
    continuous_directory: Path,
    restart_directory: Path,
) -> dict[str, Any]:
    continuous = continuous_directory / "particle_output.csv"
    restarted = restart_directory / "particle_output.csv"
    if not continuous.is_file() or not restarted.is_file():
        return {
            "matches_continuous": False,
            "reason": "particle_output.csv is missing",
        }
    continuous_hash = file_sha256(continuous)
    restarted_hash = file_sha256(restarted)
    return {
        "matches_continuous": continuous_hash == restarted_hash,
        "comparison": "byte-identical double-precision particle output",
        "continuous_sha256": continuous_hash,
        "restart_sha256": restarted_hash,
    }


def write_case_metadata(
    path: Path,
    metadata: dict[str, Any],
    case: Any,
) -> None:
    payload = {
        **metadata,
        "case": dataclasses.asdict(case),
        "recorded_at_utc": datetime.datetime.now(
            datetime.timezone.utc
        ).isoformat(),
    }
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True),
        encoding="utf-8",
    )


def run_restart_case(
    solver: Path,
    case: Any,
    run_directory: Path,
    metadata: dict[str, Any],
) -> dict[str, Any]:
    material = default_material_config()
    numerical = default_numerical_config()
    case_root = run_directory / case.name
    case_root.mkdir()
    write_case_metadata(case_root / "metadata.json", metadata, case)
    continuous_name = f"R26_{case.name.upper()}_CONT"
    split_name = f"R26_{case.name.upper()}_SPLIT"
    branches = {
        "continuous": ("new", case.stop_time, continuous_name),
        "stage1": ("new", case.restart_time, split_name),
        "restart1": ("restart_1", case.stop_time, split_name),
        "restart2": ("restart_2", case.stop_time, split_name),
    }
    for branch, (run_type, stop_time, run_name) in branches.items():
        write_case_files(
            output_directory=case_root / branch,
            case=case,
            material=material,
            numerical=numerical,
            run_name=run_name,
            run_type=run_type,
            stop_time=stop_time,
            partition=numerical.serial_partition,
        )
    statuses = {}
    for branch in ("continuous", "stage1"):
        statuses[branch] = run_solver(solver, case_root / branch)
        if not statuses[branch][0]:
            return {
                "case": case.name,
                "passed": False,
                "reason": statuses[branch][1],
                "directory": str(case_root),
            }
    for branch in ("restart1", "restart2"):
        copy_restart_artifacts(
            case_root / "stage1", case_root / branch, split_name
        )
    remove_restart_2_history(case_root / "restart2", split_name)
    for branch in ("restart1", "restart2"):
        statuses[branch] = run_solver(solver, case_root / branch)
    restart1 = compare_particle_output(
        case_root / "continuous", case_root / "restart1"
    )
    restart2 = compare_particle_output(
        case_root / "continuous", case_root / "restart2"
    )
    return {
        "case": case.name,
        "passed": (
            statuses["restart1"][0]
            and restart1["matches_continuous"]
        ),
        "directory": str(case_root),
        "solver_status": {
            branch: {"passed": value[0], "message": value[1]}
            for branch, value in statuses.items()
        },
        "restart_1": restart1,
        "restart_2_sensitivity": {
            **restart2,
            "expected_difference_observed": not restart2["matches_continuous"],
        },
    }


def run_unit_suite(run_directory: Path) -> dict[str, Any]:
    log_path = run_directory / "unit-tests.log"
    with log_path.open("w", encoding="utf-8") as log:
        completed = subprocess.run(
            [
                sys.executable,
                "-m",
                "unittest",
                "discover",
                "-s",
                "validation/tests",
                "-v",
            ],
            cwd=PROJECT_ROOT,
            check=False,
            stdout=log,
            stderr=subprocess.STDOUT,
            text=True,
        )
    return {
        "case": "unit",
        "passed": completed.returncode == 0,
        "directory": str(log_path),
        "exit_code": completed.returncode,
    }


def write_summary(run_directory: Path, results: list[dict[str, Any]]) -> None:
    payload = {
        "passed": all(result["passed"] for result in results),
        "results": results,
    }
    (run_directory / "summary.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    with (run_directory / "summary.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(
            stream, fieldnames=("case", "passed", "directory")
        )
        writer.writeheader()
        for result in results:
            writer.writerow({
                "case": result["case"],
                "passed": result["passed"],
                "directory": result["directory"],
            })


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("suite", choices=("unit", "restart", "all"))
    parser.add_argument(
        "--solver",
        type=Path,
        default=DEFAULT_SOLVER,
        help="MFIX 26.1.2 serial solver",
    )
    return parser.parse_args()


def main() -> int:
    arguments = parse_arguments()
    solver = arguments.solver.expanduser().resolve()
    if not solver.is_file() or not os.access(solver, os.X_OK):
        raise FileNotFoundError(f"solver is not executable: {solver}")
    run_directory = timestamped_run_directory()
    metadata = environment_metadata(solver)
    (run_directory / "metadata.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    results = []
    if arguments.suite in {"unit", "all"}:
        results.append(run_unit_suite(run_directory))
    if arguments.suite in {"restart", "all"}:
        for case in restart_case_configs(default_material_config()):
            try:
                results.append(
                    run_restart_case(solver, case, run_directory, metadata)
                )
            except Exception as error:
                results.append({
                    "case": case.name,
                    "passed": False,
                    "directory": str(run_directory / case.name),
                    "reason": f"{type(error).__name__}: {error}",
                })
    write_summary(run_directory, results)
    print(run_directory)
    for result in results:
        print(f"{result['case']}: {'PASS' if result['passed'] else 'FAIL'}")
    return 0 if all(result["passed"] for result in results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
