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

from validation.jkr_tools.case_generation import (
    case_contact_time,
    copy_canonical_case_asset,
    write_case_files,
)
from validation.jkr_tools.config import (
    CaseConfig,
    MaterialConfig,
    comparison_tolerances,
    damped_frictional_material_config,
    damped_jkr_normal_case_configs,
    default_material_config,
    default_numerical_config,
    elastic_hertz_case_configs,
    elastic_jkr_material_config,
    elastic_jkr_normal_case_configs,
    hertz_case_configs,
    material_config,
    restart_case_configs,
)
from validation.jkr_tools.hertz import (
    compare_case_trajectories,
    extract_contact_metrics,
)
from validation.jkr_tools.jkr import analyze_normal_jkr_trajectory
from validation.jkr_tools.reference import (
    ContactParameters,
    adhesion_only_critical_speed,
    contact_parameters,
    jkr_key_points,
    simulate_normal_contact,
)
from validation.jkr_tools.plotting import plot_contact_validation



PROJECT_ROOT = Path(__file__).resolve().parents[1]
CASES_ROOT = PROJECT_ROOT / "validation" / "cases"
DEFAULT_SOLVER = PROJECT_ROOT / "code" / "build" / "mfixsolver"
DEFAULT_NATIVE_WRAPPER = Path("/home/yuxing/miniforge3/envs/mfix-26.1.2/bin/mfixsolver")
DEFAULT_NATIVE_MODULE = Path("/home/yuxing/miniforge3/envs/mfix-26.1.2/bin/mfixsolver.so")
RESULTS_ROOT = PROJECT_ROOT / "validation" / "results"


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


def clear_results_directory(results_root: Path) -> None:
    """Remove the replaceable result snapshot while preserving .gitkeep."""
    results_root.mkdir(parents=True, exist_ok=True)
    for path in results_root.iterdir():
        if path.name == ".gitkeep":
            continue
        if path.is_dir():
            shutil.rmtree(path)
        else:
            path.unlink()


def contact_physics_subcategory(case: CaseConfig) -> str:
    participant = "wall" if case.wall_contact else "pair"
    has_tangential_velocity = any(
        abs(particle.velocity[0]) > 0.0
        or abs(particle.velocity[2]) > 0.0
        for particle in case.particles
    )
    collision = "oblique" if has_tangential_velocity else "normal"
    return f"{participant}_{collision}"


def case_result_root(results_root: Path, case: CaseConfig) -> Path:
    if case.category == "restart":
        history = "wall_history" if case.wall_contact else "pair_history"
        return results_root / "restart" / history / case.name
    return (
        results_root
        / "contact_physics"
        / contact_physics_subcategory(case)
        / case.name
    )


def write_run_info(
    path: Path,
    suite: str,
    metadata: dict[str, Any],
) -> None:
    git = metadata["jkr_git"]
    solver = metadata["solver"]
    execution = metadata["execution"]
    lines = (
        "MFIX 26.1.2 JKR validation result snapshot",
        f"recorded_at_utc={datetime.datetime.now(datetime.timezone.utc).isoformat()}",
        f"suite={suite}",
        f"mfix_version={metadata['mfix_version']}",
        f"git_commit={git['commit']}",
        f"git_branch={git['branch']}",
        f"git_dirty={git['dirty']}",
        f"solver_path={solver['path']}",
        f"solver_sha256={solver['sha256']}",
        f"execution_mode={execution['mode']}",
        f"mpi_processes={execution['mpi_processes']}",
        f"partition={execution['partition']}",
        "",
        "This directory is a replaceable snapshot. Re-run the committed",
        "validation code and canonical cases to reproduce it.",
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_solver(command: tuple[str, ...], case_directory: Path) -> tuple[bool, str]:
    log_path = case_directory / "solver.log"
    with log_path.open("w", encoding="utf-8") as log:
        completed = subprocess.run(
            command,
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


def _file_hashes(directory: Path, pattern: str) -> dict[str, str]:
    return {
        path.name: file_sha256(path)
        for path in sorted(directory.glob(pattern))
        if path.is_file()
    }


def compare_hertz_file_evidence(
    official_directory: Path,
    migrated_directory: Path,
) -> dict[str, Any]:
    official_particle = official_directory / "particle_output.csv"
    migrated_particle = migrated_directory / "particle_output.csv"
    if not official_particle.is_file() or not migrated_particle.is_file():
        return {"passed": False, "reason": "particle_output.csv is missing"}
    official_vtp = _file_hashes(official_directory, "PARTICLES_*.vtp")
    migrated_vtp = _file_hashes(migrated_directory, "PARTICLES_*.vtp")
    if not official_vtp or not migrated_vtp:
        return {"passed": False, "reason": "particle VTP trajectory is missing"}
    particle_hashes = {
        "official": file_sha256(official_particle),
        "migrated": file_sha256(migrated_particle),
    }
    shared_names = sorted(set(official_vtp) & set(migrated_vtp))
    differing_vtp = tuple(
        name for name in shared_names
        if official_vtp[name] != migrated_vtp[name]
    )
    extra_official = tuple(sorted(set(official_vtp) - set(migrated_vtp)))
    extra_migrated = tuple(sorted(set(migrated_vtp) - set(official_vtp)))
    extra_names = extra_official + extra_migrated
    last_shared_name = shared_names[-1] if shared_names else None
    last_shared_hash = official_vtp[last_shared_name] if last_shared_name else None
    trailing_duplicates = bool(shared_names) and all(
        name > last_shared_name
        and (official_vtp | migrated_vtp)[name] == last_shared_hash
        for name in extra_names
    )
    return {
        "passed": (
            particle_hashes["official"] == particle_hashes["migrated"]
            and not differing_vtp
            and (not extra_names or trailing_duplicates)
            and bool(shared_names)
        ),
        "particle_output": {
            "comparison": "byte-identical double-precision final output",
            "sha256": particle_hashes,
            "matches": particle_hashes["official"] == particle_hashes["migrated"],
        },
        "vtp": {
            "comparison": "byte-identical Float32 trajectory evidence",
            "file_count": {
                "official": len(official_vtp),
                "migrated": len(migrated_vtp),
            },
            "differing_files": differing_vtp,
            "extra_official": extra_official,
            "extra_migrated": extra_migrated,
            "extra_files_are_duplicate_final_frames": trailing_duplicates,
        },
    }


def run_hertz_case(
    solver: Path,
    native_wrapper: Path,
    native_module: Path,
    case: Any,
    case_root: Path,
    metadata: dict[str, Any],
) -> dict[str, Any]:
    material = material_config(case.material_profile)
    numerical = default_numerical_config()
    official_directory = case_root / "official"
    migrated_directory = case_root / "migrated"
    for directory in (official_directory, migrated_directory):
        copy_canonical_case_asset(
            CASES_ROOT,
            directory,
            case,
            material,
            numerical,
        )
    write_case_metadata(case_root / "metadata.json", {
        **metadata,
        "native_solver": {
            "wrapper": str(native_wrapper),
            "wrapper_sha256": file_sha256(native_wrapper),
            "module": str(native_module),
            "module_sha256": file_sha256(native_module),
        },
        "material_profile": case.material_profile,
        "material": dataclasses.asdict(material),
        "input_contract": "identical stock-compatible input without JKR extensions",
    }, case)
    statuses = {
        "official": run_solver(
            (
                str(native_wrapper),
                "-f",
                "mfix.dat",
                "-m",
                str(native_module),
            ),
            official_directory,
        ),
        "migrated": run_solver((str(solver),), migrated_directory),
    }
    if not all(status[0] for status in statuses.values()):
        return {
            "case": case.name,
            "passed": False,
            "directory": str(case_root),
            "solver_status": statuses,
        }
    comparison = compare_hertz_file_evidence(
        official_directory, migrated_directory
    )
    configured_tolerances = comparison_tolerances()
    trajectory_comparison = compare_case_trajectories(
        official_directory,
        migrated_directory,
        pattern="PARTICLES_*.vtp",
        tolerances={
            "Coordinates": configured_tolerances["vtp_position"],
            "Velocity": configured_tolerances["vtp_velocity"],
            "Angular_velocity": configured_tolerances["vtp_angular_velocity"],
            "Collision_force": configured_tolerances["vtp_force"],
        },
    )
    official_paths = tuple(sorted(official_directory.glob("PARTICLES_*.vtp")))
    migrated_paths = tuple(sorted(migrated_directory.glob("PARTICLES_*.vtp")))
    metric_arguments = {
        "contact_kind": "wall" if case.wall_contact else "pair",
        "particle_radius": material.particle_radius,
    }
    official_metrics = extract_contact_metrics(
        official_paths, **metric_arguments
    )
    migrated_metrics = extract_contact_metrics(
        migrated_paths, **metric_arguments
    )
    metrics = {
        "official": official_metrics,
        "migrated": migrated_metrics,
        "matches": official_metrics == migrated_metrics,
    }
    elastic_tolerance = 5.0e-3
    elastic_behavior_passed = (
        case.material_profile != "elastic_jkr"
        or (
            not official_metrics["captured"]
            and abs(official_metrics["effective_restitution"] - 1.0)
            <= elastic_tolerance
        )
    )
    return {
        "case": case.name,
        "passed": (
            comparison["passed"]
            and trajectory_comparison["passed"]
            and elastic_behavior_passed
            and metrics["matches"]
        ),
        "directory": str(case_root),
        "solver_status": statuses,
        "trajectory_comparison": trajectory_comparison,
        "contact_metrics": metrics,
        "comparison": comparison,
        "elastic_behavior": {
            "passed": elastic_behavior_passed,
            "effective_restitution_tolerance": elastic_tolerance,
        },
    }


def normal_contact_parameters(
    material: MaterialConfig,
    wall_contact: bool,
) -> ContactParameters:
    return contact_parameters(
        material.particle_radius,
        material.particle_density,
        material.surface_energy,
        material.young_modulus,
        material.poisson_ratio,
        material.normal_restitution,
        None if wall_contact else material.particle_radius,
        None if wall_contact else material.particle_density,
        material.wall_surface_energy if wall_contact else material.surface_energy,
        material.wall_young_modulus if wall_contact else material.young_modulus,
        material.wall_poisson_ratio if wall_contact else material.poisson_ratio,
    )


def critical_speeds(material: MaterialConfig) -> tuple[float, float]:
    return (
        adhesion_only_critical_speed(normal_contact_parameters(material, True)),
        adhesion_only_critical_speed(normal_contact_parameters(material, False)),
    )


def normal_jkr_cases(
    material: MaterialConfig,
) -> tuple[CaseConfig, ...]:
    wall_speed, pair_speed = critical_speeds(material)
    initial_gap_fraction = (
        default_numerical_config().normal_initial_gap_fraction
    )
    if material.normal_restitution == 1.0:
        return elastic_jkr_normal_case_configs(
            material, wall_speed, pair_speed, initial_gap_fraction
        )
    return damped_jkr_normal_case_configs(
        material, wall_speed, pair_speed, initial_gap_fraction
    )


def hertz_validation_cases() -> tuple[CaseConfig, ...]:
    damped = damped_frictional_material_config()
    elastic = elastic_jkr_material_config()
    wall_speed, pair_speed = critical_speeds(elastic)
    initial_gap_fraction = (
        default_numerical_config().normal_initial_gap_fraction
    )
    return (
        *hertz_case_configs(damped, initial_gap_fraction),
        *elastic_hertz_case_configs(
            elastic, wall_speed, pair_speed, initial_gap_fraction
        ),
    )


def serial_jkr_validation_cases() -> tuple[CaseConfig, ...]:
    return (
        *normal_jkr_cases(elastic_jkr_material_config()),
        *normal_jkr_cases(damped_frictional_material_config()),
    )


def canonical_contact_cases() -> tuple[CaseConfig, ...]:
    return (
        *hertz_validation_cases(),
        *serial_jkr_validation_cases(),
        *restart_case_configs(damped_frictional_material_config()),
    )


def write_rows_csv(
    path: Path,
    rows: tuple[dict[str, Any], ...],
) -> None:
    if not rows:
        raise ValueError("cannot write an empty evidence CSV")
    fieldnames = tuple(dict.fromkeys(
        key for row in rows for key in row
    ))
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def key_points_have_expected_signs(
    parameters: ContactParameters,
) -> bool:
    points = {point.name: point for point in jkr_key_points(parameters)}
    return (
        points["critical_separation"].overlap < 0.0
        and points["critical_separation"].normal_force < 0.0
        and points["maximum_tension"].overlap < 0.0
        and points["maximum_tension"].normal_force
        == -parameters.pull_off_force_magnitude
        and abs(points["zero_overlap"].overlap) <= 1.0e-22
        and points["zero_overlap"].normal_force < 0.0
        and points["zero_force_equilibrium"].overlap > 0.0
        and abs(points["zero_force_equilibrium"].normal_force) <= 1.0e-22
    )


def run_normal_jkr_case(
    solver: Path,
    case: CaseConfig,
    case_directory: Path,
    metadata: dict[str, Any],
) -> dict[str, Any]:
    material = material_config(case.material_profile)
    numerical = default_numerical_config()
    parameters = normal_contact_parameters(material, case.wall_contact)
    critical_speed = adhesion_only_critical_speed(parameters)
    expected_behavior = "capture" if case.name.endswith("capture") else "rebound"
    copy_canonical_case_asset(
        CASES_ROOT,
        case_directory,
        case,
        material,
        numerical,
    )
    force_absolute_tolerance = 1.0e-9
    force_relative_tolerance = 2.0e-2
    rebound_velocity_tolerance = 1.0e-6
    key_points = tuple(jkr_key_points(parameters))
    reference_metadata = {
        "contact_parameters": dataclasses.asdict(parameters),
        "key_points": [dataclasses.asdict(point) for point in key_points],
        "key_point_signs_passed": key_points_have_expected_signs(parameters),
        "adhesion_only_critical_speed": critical_speed,
        "impact_speed_factor": case.critical_speed_factor,
        "expected_behavior": expected_behavior,
    }
    tolerance_metadata = {
        "vtp_force_absolute": force_absolute_tolerance,
        "vtp_force_relative": force_relative_tolerance,
        "rebound_velocity": rebound_velocity_tolerance,
    }
    base_metadata = {
        **metadata,
        "material_profile": case.material_profile,
        "material": dataclasses.asdict(material),
        "reference": reference_metadata,
        "tolerances": tolerance_metadata,
    }
    write_case_metadata(
        case_directory / "metadata.json", base_metadata, case
    )
    solver_status = run_solver((str(solver),), case_directory)
    if not solver_status[0]:
        return {
            "case": case.name,
            "passed": False,
            "directory": str(case_directory),
            "solver_status": solver_status,
        }
    paths = tuple(sorted(case_directory.glob("PARTICLES_*.vtp")))
    analysis = analyze_normal_jkr_trajectory(
        paths,
        "wall" if case.wall_contact else "pair",
        parameters,
        material.particle_radius,
        force_absolute_tolerance,
        force_relative_tolerance,
        rebound_velocity_tolerance,
    )
    mfix_rows = tuple(analysis.pop("samples"))
    sample_times = tuple(float(row["time"]) for row in mfix_rows)
    initial_overlap = (
        -case.expected_initial_overlap
        if case.wall_contact
        else case.expected_initial_overlap
    )
    reference_rows = simulate_normal_contact(
        sample_times,
        initial_overlap,
        case.contact_time_speed,
        parameters,
        case_contact_time(case, material) / 2000.0,
    )
    mfix_csv = case_directory / "mfix_contact.csv"
    reference_csv = case_directory / "python_reference_contact.csv"
    write_rows_csv(mfix_csv, mfix_rows)
    write_rows_csv(reference_csv, reference_rows)
    plot_paths = plot_contact_validation(
        mfix_rows,
        reference_rows,
        case_directory,
        case.name,
    )
    if expected_behavior == "capture":
        behavior_passed = (
            analysis["captured"]
            and analysis["final_contact_active"]
            and analysis["final_overlap"] >= -parameters.critical_separation
        )
    else:
        behavior_passed = (
            not analysis["captured"]
            and not analysis["final_contact_active"]
            and analysis["rebound_velocity"] > rebound_velocity_tolerance
            and analysis["final_overlap"] < -parameters.critical_separation
        )
    passed = (
        analysis["passed"]
        and behavior_passed
        and reference_metadata["key_point_signs_passed"]
    )
    evidence = {
        "mfix_csv": str(mfix_csv),
        "python_reference_csv": str(reference_csv),
        "plots": [str(path) for path in plot_paths],
        "force_comparison": analysis["force_comparison"],
    }
    write_case_metadata(
        case_directory / "metadata.json",
        {
            **base_metadata,
            "comparison_evidence": evidence,
            "trajectory_summary": analysis,
            "behavior_passed": behavior_passed,
        },
        case,
    )
    return {
        "case": case.name,
        "passed": passed,
        "directory": str(case_directory),
        "solver_status": solver_status,
        "expected_behavior": expected_behavior,
        "behavior_passed": behavior_passed,
        "trajectory": analysis,
        "evidence": evidence,
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
    case_root: Path,
    metadata: dict[str, Any],
) -> dict[str, Any]:
    material = default_material_config()
    numerical = default_numerical_config()
    case_root.mkdir(parents=True)
    write_case_metadata(case_root / "metadata.json", metadata, case)
    split_name = f"R26_{case.name.upper()}_SPLIT"
    copy_canonical_case_asset(
        CASES_ROOT,
        case_root / "continuous",
        case,
        material,
        numerical,
    )
    branches = {
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
        statuses[branch] = run_solver((str(solver),), case_root / branch)
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
        statuses[branch] = run_solver((str(solver),), case_root / branch)
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


def run_unit_suite(unit_directory: Path) -> dict[str, Any]:
    unit_directory.mkdir(parents=True, exist_ok=False)
    log_path = unit_directory / "unit-tests.log"
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


def write_summary(results_root: Path, results: list[dict[str, Any]]) -> None:
    payload = {
        "passed": all(result["passed"] for result in results),
        "results": results,
    }
    (results_root / "summary.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    with (results_root / "summary.csv").open(
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
    parser.add_argument(
        "suite", choices=("unit", "hertz", "serial", "restart", "all")
    )
    parser.add_argument(
        "--solver",
        type=Path,
        default=DEFAULT_SOLVER,
        help="MFIX 26.1.2 serial solver",
    )
    parser.add_argument(
        "--native-wrapper",
        type=Path,
        default=DEFAULT_NATIVE_WRAPPER,
        help="official MFIX 26.1.2 pymfix wrapper",
    )
    parser.add_argument(
        "--native-module",
        type=Path,
        default=DEFAULT_NATIVE_MODULE,
        help="official unmodified MFIX 26.1.2 solver module",
    )
    return parser.parse_args()


def main() -> int:
    arguments = parse_arguments()
    solver = arguments.solver.expanduser().resolve()
    if not solver.is_file() or not os.access(solver, os.X_OK):
        raise FileNotFoundError(f"solver is not executable: {solver}")
    native_wrapper = arguments.native_wrapper.expanduser().resolve()
    native_module = arguments.native_module.expanduser().resolve()
    if arguments.suite in {"hertz", "all"}:
        if not native_wrapper.is_file() or not os.access(native_wrapper, os.X_OK):
            raise FileNotFoundError(
                f"native solver wrapper is not executable: {native_wrapper}"
            )
        if not native_module.is_file():
            raise FileNotFoundError(
                f"native solver module does not exist: {native_module}"
            )
    clear_results_directory(RESULTS_ROOT)
    metadata = environment_metadata(solver)
    (RESULTS_ROOT / "metadata.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    write_run_info(RESULTS_ROOT / "RUN_INFO.txt", arguments.suite, metadata)
    results = []
    if arguments.suite in {"unit", "all"}:
        results.append(run_unit_suite(RESULTS_ROOT / "unit"))
    if arguments.suite in {"hertz", "all"}:
        for case in hertz_validation_cases():
            case_root = case_result_root(RESULTS_ROOT, case)
            try:
                results.append(run_hertz_case(
                    solver,
                    native_wrapper,
                    native_module,
                    case,
                    case_root,
                    metadata,
                ))
            except Exception as error:
                results.append({
                    "case": case.name,
                    "passed": False,
                    "directory": str(case_root),
                    "reason": f"{type(error).__name__}: {error}",
                })
    if arguments.suite in {"serial", "all"}:
        for case in serial_jkr_validation_cases():
            case_root = case_result_root(RESULTS_ROOT, case)
            serial_directory = case_root / "serial"
            try:
                results.append(
                    run_normal_jkr_case(
                        solver, case, serial_directory, metadata
                    )
                )
            except Exception as error:
                results.append({
                    "case": case.name,
                    "passed": False,
                    "directory": str(serial_directory),
                    "reason": f"{type(error).__name__}: {error}",
                })
    if arguments.suite in {"restart", "all"}:
        for case in restart_case_configs(default_material_config()):
            case_root = case_result_root(RESULTS_ROOT, case)
            try:
                results.append(
                    run_restart_case(solver, case, case_root, metadata)
                )
            except Exception as error:
                results.append({
                    "case": case.name,
                    "passed": False,
                    "directory": str(case_root),
                    "reason": f"{type(error).__name__}: {error}",
                })
    write_summary(RESULTS_ROOT, results)
    print(RESULTS_ROOT)
    for result in results:
        print(f"{result['case']}: {'PASS' if result['passed'] else 'FAIL'}")
    return 0 if all(result["passed"] for result in results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
