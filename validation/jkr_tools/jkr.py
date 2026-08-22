"""Independent normal JKR trajectory checks for single-contact cases."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any, Sequence

from .hertz import read_vtp_trajectory
from .reference import (
    ContactParameters,
    contact_is_active,
    elastic_normal_force,
    normal_damping_coefficient,
    stable_contact_radius,
)


def _dot(first: Sequence[float], second: Sequence[float]) -> float:
    return sum(float(left) * float(right) for left, right in zip(first, second))


def _vector(record: dict[str, Any], field: str) -> tuple[float, float, float]:
    values = tuple(float(value) for value in record[field])
    if len(values) != 3 or any(not math.isfinite(value) for value in values):
        raise ValueError(f"{field} must contain three finite values")
    return values


def _normal_sample(
    records: dict[int, dict[str, Any]],
    contact_kind: str,
    particle_ids: tuple[int, ...],
    particle_radius: float,
) -> tuple[float, float, float]:
    selected = tuple(records[particle_id] for particle_id in particle_ids)
    positions = tuple(_vector(record, "Coordinates") for record in selected)
    velocities = tuple(_vector(record, "Velocity") for record in selected)
    forces = tuple(_vector(record, "Collision_force") for record in selected)
    if contact_kind == "wall":
        wall_normal = (0.0, 1.0, 0.0)
        overlap = particle_radius - positions[0][1]
        normal_speed = _dot(velocities[0], wall_normal)
        normal_force = _dot(forces[0], wall_normal)
        return overlap, normal_speed, normal_force
    displacement = tuple(
        second - first for first, second in zip(positions[0], positions[1])
    )
    distance = math.sqrt(_dot(displacement, displacement))
    if distance == 0.0:
        raise ValueError("pair particle centers coincide")
    direction = tuple(value / distance for value in displacement)
    relative_velocity = tuple(
        second - first for first, second in zip(velocities[0], velocities[1])
    )
    overlap = 2.0 * particle_radius - distance
    normal_speed = _dot(relative_velocity, direction)
    normal_force = -_dot(forces[0], direction)
    return overlap, normal_speed, normal_force


def analyze_normal_jkr_trajectory(
    paths: Sequence[Path],
    contact_kind: str,
    parameters: ContactParameters,
    particle_radius: float,
    force_absolute_tolerance: float,
    force_relative_tolerance: float,
    rebound_velocity_tolerance: float,
) -> dict[str, Any]:
    """Compare a normal single-contact trajectory with the JKR force law."""
    if contact_kind not in {"wall", "pair"}:
        raise ValueError("contact_kind must be 'wall' or 'pair'")
    if min(
        force_absolute_tolerance,
        force_relative_tolerance,
        rebound_velocity_tolerance,
    ) < 0.0:
        raise ValueError("comparison tolerances must be non-negative")
    trajectory = read_vtp_trajectory(paths)
    ids = tuple(sorted(int(row["Particle_ID"]) for row in trajectory[0]["records"]))
    required_count = 1 if contact_kind == "wall" else 2
    if len(ids) != required_count:
        raise ValueError(f"{contact_kind} analysis requires {required_count} particles")

    samples = []
    active = False
    first_active_time = None
    last_active_time = None
    maximum_absolute_error = 0.0
    maximum_relative_error = 0.0
    maximum_error_time = None
    first_difference = None
    rebound_velocity = 0.0
    for frame in trajectory:
        records = {
            int(row["Particle_ID"]): row for row in frame["records"]
        }
        overlap, normal_speed, observed_force = _normal_sample(
            records, contact_kind, ids, particle_radius
        )
        was_active = active
        active = contact_is_active(
            overlap, active, parameters.critical_separation
        )
        time = float(frame["time"])
        expected_force = None
        absolute_error = None
        relative_error = None
        passed = True
        if active:
            contact_radius = stable_contact_radius(overlap, parameters)
            expected_force = (
                elastic_normal_force(contact_radius, parameters)
                - normal_damping_coefficient(contact_radius, parameters)
                * normal_speed
            )
            absolute_error = abs(observed_force - expected_force)
            relative_error = absolute_error / max(
                abs(expected_force), force_absolute_tolerance
            )
            passed = (
                absolute_error <= force_absolute_tolerance
                or relative_error <= force_relative_tolerance
            )
            first_active_time = time if first_active_time is None else first_active_time
            last_active_time = time
            if absolute_error > maximum_absolute_error:
                maximum_absolute_error = absolute_error
                maximum_relative_error = relative_error
                maximum_error_time = time
            if not passed and first_difference is None:
                first_difference = {
                    "time": time,
                    "overlap": overlap,
                    "normal_speed": normal_speed,
                    "observed_force": observed_force,
                    "expected_force": expected_force,
                    "absolute_error": absolute_error,
                    "relative_error": relative_error,
                }
        elif was_active:
            rebound_velocity = max(rebound_velocity, max(0.0, normal_speed))
        elif first_active_time is not None:
            rebound_velocity = max(rebound_velocity, max(0.0, normal_speed))
        samples.append({
            "time": time,
            "overlap": overlap,
            "normal_speed": normal_speed,
            "observed_normal_force": observed_force,
            "expected_normal_force": expected_force,
            "active": active,
            "passed": passed,
        })

    impact_velocity = max(0.0, -samples[0]["normal_speed"])
    contact_duration = (
        0.0
        if first_active_time is None or last_active_time is None
        else last_active_time - first_active_time
    )
    return {
        "passed": first_active_time is not None and first_difference is None,
        "contact_kind": contact_kind,
        "particle_ids": ids,
        "sample_count": len(samples),
        "active_sample_count": sum(sample["active"] for sample in samples),
        "impact_velocity": impact_velocity,
        "rebound_velocity": rebound_velocity,
        "effective_restitution": (
            rebound_velocity / impact_velocity if impact_velocity > 0.0 else 0.0
        ),
        "captured": active and rebound_velocity <= rebound_velocity_tolerance,
        "final_contact_active": active,
        "contact_duration": contact_duration,
        "maximum_overlap": max(sample["overlap"] for sample in samples),
        "minimum_overlap": min(sample["overlap"] for sample in samples),
        "final_overlap": samples[-1]["overlap"],
        "final_normal_speed": samples[-1]["normal_speed"],
        "maximum_normal_force": max(
            sample["observed_normal_force"] for sample in samples
        ),
        "minimum_normal_force": min(
            sample["observed_normal_force"] for sample in samples
        ),
        "force_comparison": {
            "absolute_tolerance": force_absolute_tolerance,
            "relative_tolerance": force_relative_tolerance,
            "maximum_absolute_error": maximum_absolute_error,
            "relative_error_at_maximum_absolute_error": maximum_relative_error,
            "maximum_error_time": maximum_error_time,
            "first_difference": first_difference,
        },
    }
