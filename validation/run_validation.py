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
    dmp_case_configs,
    default_numerical_config,
    elastic_hertz_case_configs,
    elastic_jkr_material_config,
    elastic_jkr_normal_case_configs,
    hertz_case_configs,
    heterogeneous_surface_energy_case_configs,
    material_config,
    restart_case_configs,
)
from validation.jkr_tools.hertz import (
    compare_case_trajectories,
    extract_contact_metrics,
)
from validation.jkr_tools.dmp import audit_serial_dmp_trajectories
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
DEFAULT_DMP_SOLVER = PROJECT_ROOT / "code" / "build_dmp" / "mfixsolver_dmp"
DEFAULT_NATIVE_WRAPPER = Path("/home/yuxing/miniforge3/envs/mfix-26.1.2/bin/mfixsolver")
DEFAULT_NATIVE_MODULE = Path("/home/yuxing/miniforge3/envs/mfix-26.1.2/bin/mfixsolver.so")
DEFAULT_MPIEXEC = Path("/usr/bin/mpiexec")
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


def environment_metadata(
    solver: Path,
    execution_mode: str = "serial",
    mpi_processes: int = 1,
    partition: tuple[int, int, int] | None = None,
) -> dict[str, Any]:
    material = default_material_config()
    numerical = default_numerical_config()
    actual_partition = numerical.serial_partition if partition is None else partition
    return {
        "mfix_version": command_text([str(solver), "--version"]),
        "jkr_git": git_metadata(),
        "solver": {
            "path": str(solver),
            "sha256": file_sha256(solver),
        },
        "execution": {
            "mode": execution_mode,
            "mpi_processes": mpi_processes,
            "partition": actual_partition,
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
    if case.category == "dmp":
        subcategories = {
            "pair_cross_jkr": "pair_cross",
            "pair_new_contact_jkr": "pair_new_contact",
            "pair_active_contact_migration": "active_contact_migration",
            "wall_cross_jkr": "wall_cross",
            "cloud_collision_jkr": "cloud",
        }
        return results_root / "dmp" / subcategories[case.name] / case.name
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


def resolved_particle_surface_energies(
    case: CaseConfig,
    material: MaterialConfig,
) -> tuple[float, ...]:
    """Apply the production positive-override/non-positive-fallback contract."""
    return tuple(
        material.surface_energy
        if particle.surface_energy is None or particle.surface_energy <= 0.0
        else particle.surface_energy
        for particle in case.particles
    )


def case_contact_parameters(
    case: CaseConfig,
    material: MaterialConfig,
) -> ContactParameters:
    """Build explicit reference parameters from the case's physical surfaces."""
    resolved = resolved_particle_surface_energies(case, material)
    first = case.particles[0]
    second = None if case.wall_contact else case.particles[1]
    return contact_parameters(
        0.5 * first.diameter,
        first.density,
        resolved[0],
        material.young_modulus,
        material.poisson_ratio,
        material.normal_restitution,
        None if second is None else 0.5 * second.diameter,
        None if second is None else second.density,
        material.wall_surface_energy if second is None else resolved[1],
        material.wall_young_modulus if second is None else material.young_modulus,
        material.wall_poisson_ratio if second is None else material.poisson_ratio,
    )


def surface_energy_contract(
    case: CaseConfig,
    material: MaterialConfig,
    parameters: ContactParameters,
) -> dict[str, Any]:
    """Return auditable raw, resolved, and combined surface-energy evidence."""
    raw = tuple(particle.surface_energy for particle in case.particles)
    resolved = resolved_particle_surface_energies(case, material)
    contract = {
        "particle_values_raw": raw,
        "particle_values_resolved": resolved,
        "phase_default": material.surface_energy,
        "fallback_rule": "missing or non-positive particle value uses phase default",
        "gamma_i": resolved[0],
        "gamma_effective": parameters.effective_surface_energy,
        "work_of_adhesion": parameters.work_of_adhesion,
    }
    if case.wall_contact:
        contract["gamma_wall"] = material.wall_surface_energy
        contract["combination_rule"] = (
            "sqrt(gamma_particle * gamma_wall)"
        )
    else:
        contract["gamma_j"] = resolved[1]
        contract["combination_rule"] = "sqrt(gamma_i * gamma_j)"
    return contract


def heterogeneous_surface_energy_cases() -> tuple[CaseConfig, ...]:
    """Return the fixed analytical-speed heterogeneous serial registry."""
    material = elastic_jkr_material_config()
    wall_speed, pair_speed = critical_speeds(material)
    return heterogeneous_surface_energy_case_configs(
        material,
        0.95 * wall_speed,
        0.95 * pair_speed,
        default_numerical_config().normal_initial_gap_fraction,
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
        *heterogeneous_surface_energy_cases(),
    )



def dmp_validation_cases() -> tuple[CaseConfig, ...]:
    return dmp_case_configs(
        elastic_jkr_material_config(),
        damped_frictional_material_config(),
    )


def canonical_contact_cases() -> tuple[CaseConfig, ...]:
    return (
        *hertz_validation_cases(),
        *serial_jkr_validation_cases(),
        *dmp_validation_cases(),
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
    parameters = case_contact_parameters(case, material)
    critical_speed = adhesion_only_critical_speed(parameters)
    expected_behavior = (
        "capture"
        if case.critical_speed_factor is not None
        and case.critical_speed_factor < 1.0
        else "rebound"
    )
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
        "surface_energy_contract": surface_energy_contract(case, material, parameters),
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



def configure_partition(
    mfix_path: Path,
    partition: tuple[int, int, int],
) -> None:
    """Change only the three MPI partition keywords in a generated input."""
    replacements = {
        "nodesi": partition[0],
        "nodesj": partition[1],
        "nodesk": partition[2],
    }
    output = []
    seen = set()
    for line in mfix_path.read_text(encoding="utf-8").splitlines():
        keyword = line.partition("=")[0].strip().casefold()
        if keyword in replacements:
            output.append(f"{keyword} = {replacements[keyword]}")
            seen.add(keyword)
        else:
            output.append(line)
    if seen != set(replacements):
        raise ValueError(f"partition keywords missing from {mfix_path}: {seen}")
    mfix_path.write_text("\n".join(output) + "\n", encoding="utf-8")


def partition_independent_mfix_sha256(path: Path) -> str:
    """Hash an mfix.dat after normalizing only nodesi/nodesj/nodesk."""
    lines = []
    for line in path.read_text(encoding="utf-8").splitlines():
        keyword = line.partition("=")[0].strip().casefold()
        lines.append(f"{keyword} = <partition>" if keyword in {
            "nodesi", "nodesj", "nodesk"
        } else line)
    return hashlib.sha256(("\n".join(lines) + "\n").encode("utf-8")).hexdigest()


def dmp_case_contract(case: CaseConfig, audit: Mapping[str, Any]) -> dict[str, Any]:
    """Evaluate case-specific natural contact/migration prerequisites."""
    contact = audit["contact_set_comparison"]
    first_contact = contact["first_active_time"]
    histories = audit["rank_transitions"]
    transition_times = tuple(
        float(history[1]["time"])
        for history in histories.values()
        if len(history) > 1
    )
    initial_ranks = {
        int(history[0]["rank"]) for history in histories.values() if history
    }
    passed = contact["active_frame_count"] > 0
    requirements = {
        "contact_naturally_formed": contact["active_frame_count"] > 0,
        "initial_ranks": tuple(sorted(initial_ranks)),
        "rank_transition_times": transition_times,
        "first_contact_time": first_contact,
    }
    if case.name == "pair_cross_jkr":
        requirements["particles_initially_on_distinct_ranks"] = len(initial_ranks) == 2
        passed = passed and requirements["particles_initially_on_distinct_ranks"]
    elif case.name == "pair_new_contact_jkr":
        formed_after_cross = bool(transition_times) and first_contact is not None and (
            min(transition_times) <= float(first_contact)
        )
        requirements["first_contact_after_rank_cross"] = formed_after_cross
        passed = passed and formed_after_cross
    elif case.name == "pair_active_contact_migration":
        active_before_migration = (
            first_contact is not None
            and any(time > float(first_contact) for time in transition_times)
        )
        histories_nonzero = audit["maximums"]["serial_angular_speed"] > 0.0
        requirements["active_contact_before_owner_migration"] = active_before_migration
        requirements["nonzero_rotation_from_natural_oblique_contact"] = histories_nonzero
        passed = passed and active_before_migration and histories_nonzero
    elif case.name == "wall_cross_jkr":
        wall_history_crossed = bool(transition_times) and first_contact is not None
        requirements["active_wall_contact_crosses_rank"] = wall_history_crossed
        passed = passed and wall_history_crossed
    return {"passed": passed, **requirements}


def run_dmp_case(
    serial_solver: Path,
    dmp_solver: Path,
    mpiexec: Path,
    case: CaseConfig,
    case_root: Path,
    metadata: dict[str, Any],
) -> dict[str, Any]:
    """Run byte-synchronized serial/2-rank inputs and compare by physical IDs."""
    material = material_config(case.material_profile)
    numerical = default_numerical_config()
    serial_directory = case_root / "serial"
    dmp_directory = case_root / "dmp"
    copy_canonical_case_asset(
        CASES_ROOT, serial_directory, case, material, numerical
    )
    copy_canonical_case_asset(
        CASES_ROOT, dmp_directory, case, material, numerical
    )
    configure_partition(dmp_directory / "mfix.dat", numerical.dmp_partition)
    input_contract = {
        "particle_input_byte_identical": (
            file_sha256(serial_directory / "particle_input.dat")
            == file_sha256(dmp_directory / "particle_input.dat")
        ),
        "mfix_physics_byte_identical_except_partition": (
            partition_independent_mfix_sha256(serial_directory / "mfix.dat")
            == partition_independent_mfix_sha256(dmp_directory / "mfix.dat")
        ),
        "serial_partition": numerical.serial_partition,
        "dmp_partition": numerical.dmp_partition,
    }
    base_metadata = {
        **metadata,
        "material_profile": case.material_profile,
        "material": dataclasses.asdict(material),
        "serial_solver": {
            "path": str(serial_solver),
            "sha256": file_sha256(serial_solver),
        },
        "dmp_solver": {
            "path": str(dmp_solver),
            "sha256": file_sha256(dmp_solver),
        },
        "executions": {
            "serial": {
                "execution_mode": "serial",
                "mpi_processes": 1,
                "partition": numerical.serial_partition,
            },
            "dmp": {
                "execution_mode": "dmp",
                "mpi_processes": numerical.dmp_ranks,
                "partition": numerical.dmp_partition,
            },
        },
        "input_contract": input_contract,
    }
    write_case_metadata(case_root / "metadata.json", base_metadata, case)
    serial_status = run_solver((str(serial_solver),), serial_directory)
    dmp_status = run_solver(
        (
            str(mpiexec),
            "-np",
            str(numerical.dmp_ranks),
            str(dmp_solver),
        ),
        dmp_directory,
    )
    statuses_passed = serial_status[0] and dmp_status[0]
    if not statuses_passed:
        return {
            "case": case.name,
            "passed": False,
            "directory": str(case_root),
            "solver_status": {"serial": serial_status, "dmp": dmp_status},
            "input_contract": input_contract,
        }
    tolerances = comparison_tolerances()
    trajectory_tolerances = {
        "Coordinates": tolerances["vtp_position"],
        "Velocity": tolerances["dmp_vtp_velocity"],
        "Angular_velocity": tolerances["dmp_vtp_angular_velocity"],
        "Collision_force": tolerances["vtp_force"],
    }
    strict_reporting_tolerances = {
        "Coordinates": tolerances["vtp_position"],
        "Velocity": tolerances["vtp_velocity"],
        "Angular_velocity": tolerances["vtp_angular_velocity"],
        "Collision_force": tolerances["vtp_force"],
    }
    serial_paths = tuple(sorted(serial_directory.glob("PARTICLES_*.vtp")))
    dmp_paths = tuple(sorted(dmp_directory.glob("PARTICLES_*.vtp")))
    audit = audit_serial_dmp_trajectories(
        serial_paths,
        dmp_paths,
        trajectory_tolerances,
        case.wall_contact,
    )
    strict_audit = audit_serial_dmp_trajectories(
        serial_paths,
        dmp_paths,
        strict_reporting_tolerances,
        case.wall_contact,
    )
    strict_differences = tuple(
        {
            **comparison["first_difference"],
            "comparison_field": field,
        }
        for field, comparison in strict_audit["trajectory_comparison"].items()
        if comparison["first_difference"] is not None
    )
    audit["first_physical_difference"] = min(
        strict_differences,
        key=lambda difference: difference["time"],
        default=None,
    )
    audit["strict_reporting_tolerances"] = strict_reporting_tolerances
    contract = dmp_case_contract(case, audit)
    inputs_match = (
        input_contract["particle_input_byte_identical"]
        and input_contract["mfix_physics_byte_identical_except_partition"]
    )
    passed = inputs_match and audit["passed"] and contract["passed"]
    evidence_path = case_root / "comparison.json"
    evidence_path.write_text(
        json.dumps(
            {
                "passed": passed,
                "input_contract": input_contract,
                "audit": audit,
                "case_contract": contract,
            },
            indent=2,
            sort_keys=True,
        ),
        encoding="utf-8",
    )
    write_case_metadata(
        case_root / "metadata.json",
        {
            **base_metadata,
            "tolerances": trajectory_tolerances,
            "audit": audit,
            "case_contract": contract,
            "comparison_path": str(evidence_path),
        },
        case,
    )
    return {
        "case": case.name,
        "passed": passed,
        "directory": str(case_root),
        "solver_status": {"serial": serial_status, "dmp": dmp_status},
        "input_contract": input_contract,
        "audit": audit,
        "case_contract": contract,
        "evidence": str(evidence_path),
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
        "suite", choices=("unit", "hertz", "serial", "dmp", "restart", "all")
    )
    parser.add_argument(
        "--solver",
        type=Path,
        default=DEFAULT_SOLVER,
        help="MFIX 26.1.2 serial solver",
    )
    parser.add_argument(
        "--dmp-solver",
        type=Path,
        default=DEFAULT_DMP_SOLVER,
        help="MFIX 26.1.2 DMP solver",
    )
    parser.add_argument(
        "--mpiexec",
        type=Path,
        default=DEFAULT_MPIEXEC,
        help="MPI launcher for the fixed 2-rank cases",
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
    dmp_solver = arguments.dmp_solver.expanduser().resolve()
    mpiexec = arguments.mpiexec.expanduser().resolve()
    if arguments.suite in {"dmp", "all"}:
        if not dmp_solver.is_file() or not os.access(dmp_solver, os.X_OK):
            raise FileNotFoundError(
                f"DMP solver is not executable: {dmp_solver}"
            )
        if not mpiexec.is_file() or not os.access(mpiexec, os.X_OK):
            raise FileNotFoundError(f"MPI launcher is not executable: {mpiexec}")

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
    numerical = default_numerical_config()
    uses_dmp = arguments.suite in {"dmp", "all"}
    execution_mode = (
        "dmp"
        if arguments.suite == "dmp"
        else "mixed_serial_dmp"
        if uses_dmp
        else "serial"
    )
    metadata = environment_metadata(
        solver,
        execution_mode=execution_mode,
        mpi_processes=numerical.dmp_ranks if uses_dmp else 1,
        partition=(
            numerical.dmp_partition
            if uses_dmp
            else numerical.serial_partition
        ),
    )
    if uses_dmp:
        metadata["dmp_solver"] = {
            "path": str(dmp_solver),
            "sha256": file_sha256(dmp_solver),
        }
        metadata["mpiexec"] = str(mpiexec)
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
    if arguments.suite in {"dmp", "all"}:
        prerequisites_passed = True
        for case in dmp_validation_cases():
            case_root = case_result_root(RESULTS_ROOT, case)
            if case.name == "cloud_collision_jkr" and not prerequisites_passed:
                results.append({
                    "case": case.name,
                    "passed": False,
                    "directory": str(case_root),
                    "reason": (
                        "not run: one or more focused DMP prerequisites failed"
                    ),
                })
                continue
            try:
                result = run_dmp_case(
                    solver,
                    dmp_solver,
                    mpiexec,
                    case,
                    case_root,
                    metadata,
                )
            except Exception as error:
                result = {
                    "case": case.name,
                    "passed": False,
                    "directory": str(case_root),
                    "reason": f"{type(error).__name__}: {error}",
                }
            results.append(result)
            if case.name != "cloud_collision_jkr":
                prerequisites_passed = prerequisites_passed and result["passed"]
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
