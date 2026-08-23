"""Pure functions for serial/2-rank MFIX particle trajectory audits."""

from __future__ import annotations

import math
from collections import defaultdict
from pathlib import Path
from typing import Any, Mapping, Sequence

from .hertz import (
    compare_trajectories,
    frame_to_records,
    trajectory_from_records,
)
from .io import read_vtp


def read_partitioned_vtp_trajectory(
    paths: Sequence[Path],
) -> tuple[dict[str, Any], ...]:
    """Merge owned-particle VTP pieces by physical output time."""
    if not paths:
        raise ValueError("partitioned trajectory contains no VTP files")
    records_by_time: dict[float, list[dict[str, Any]]] = defaultdict(list)
    ids_by_time: dict[float, set[int]] = defaultdict(set)
    records_by_time_id: dict[tuple[float, int], dict[str, Any]] = {}
    for path in sorted(Path(item) for item in paths):
        frame = frame_to_records(read_vtp(path))
        time = float(frame["time"])
        for source in frame["records"]:
            record = dict(source)
            particle_id = int(record["Particle_ID"])
            if particle_id in ids_by_time[time]:
                if records_by_time_id[(time, particle_id)] == record:
                    continue
                raise ValueError(
                    "conflicting duplicate Particle_ID "
                    f"{particle_id} across VTP files at t={time}"
                )
            ids_by_time[time].add(particle_id)
            records_by_time[time].append(record)
            records_by_time_id[(time, particle_id)] = record
    frames = tuple(
        {"time": time, "records": tuple(records_by_time[time])}
        for time in sorted(records_by_time)
    )
    return trajectory_from_records(frames)


def _components(record: Mapping[str, Any], field: str) -> tuple[float, ...]:
    value = record[field]
    if isinstance(value, (tuple, list)):
        return tuple(float(item) for item in value)
    return (float(value),)


def _within(
    left: Mapping[str, Any],
    right: Mapping[str, Any],
    field: str,
    tolerance: float,
) -> bool:
    first = _components(left, field)
    second = _components(right, field)
    return len(first) == len(second) and all(
        abs(a - b) <= tolerance for a, b in zip(first, second)
    )


def initial_physical_id_mapping(
    serial: Sequence[Mapping[str, Any]],
    dmp: Sequence[Mapping[str, Any]],
    position_tolerance: float,
    velocity_tolerance: float,
) -> dict[int, int]:
    """Map serial IDs to DMP IDs using the initial physical particle state."""
    if min(position_tolerance, velocity_tolerance) < 0.0:
        raise ValueError("mapping tolerances must be non-negative")
    serial_records = tuple(serial[0]["records"])
    dmp_records = tuple(dmp[0]["records"])
    if len(serial_records) != len(dmp_records):
        raise ValueError("serial and DMP initial particle counts differ")
    available = {int(record["Particle_ID"]): record for record in dmp_records}
    mapping: dict[int, int] = {}
    for left in serial_records:
        candidates = []
        for right_id, right in available.items():
            same_surface = (
                "Surface_Energy" not in left
                or "Surface_Energy" not in right
                or _within(left, right, "Surface_Energy", 1.0e-12)
            )
            if (
                same_surface
                and _within(left, right, "Coordinates", position_tolerance)
                and _within(left, right, "Velocity", velocity_tolerance)
                and _within(left, right, "Diameter", position_tolerance)
            ):
                candidates.append(right_id)
        if len(candidates) != 1:
            raise ValueError(
                "initial physical state does not define a unique ID mapping: "
                f"serial_id={left['Particle_ID']}, candidates={candidates}"
            )
        serial_id = int(left["Particle_ID"])
        dmp_id = candidates[0]
        mapping[serial_id] = dmp_id
        del available[dmp_id]
    return mapping


def remap_to_serial_ids(
    dmp: Sequence[Mapping[str, Any]],
    serial_to_dmp: Mapping[int, int],
) -> tuple[dict[str, Any], ...]:
    """Return DMP frames keyed by the corresponding serial global IDs."""
    dmp_to_serial = {dmp_id: serial_id for serial_id, dmp_id in serial_to_dmp.items()}
    if len(dmp_to_serial) != len(serial_to_dmp):
        raise ValueError("DMP ID mapping is not one-to-one")
    frames = []
    for frame in dmp:
        records = []
        for source in frame["records"]:
            dmp_id = int(source["Particle_ID"])
            if dmp_id not in dmp_to_serial:
                raise ValueError(f"unmapped DMP Particle_ID {dmp_id}")
            record = dict(source)
            record["Particle_ID"] = dmp_to_serial[dmp_id]
            records.append(record)
        frames.append({"time": float(frame["time"]), "records": tuple(records)})
    return trajectory_from_records(tuple(frames))


def rank_transitions(
    trajectory: Sequence[Mapping[str, Any]],
) -> dict[int, tuple[dict[str, float | int], ...]]:
    """Return rank changes by remapped global particle ID."""
    histories: dict[int, list[dict[str, float | int]]] = defaultdict(list)
    previous: dict[int, int] = {}
    for frame in trajectory:
        time = float(frame["time"])
        for record in frame["records"]:
            particle_id = int(record["Particle_ID"])
            rank = int(round(float(record["Particle_Rank"])))
            if particle_id not in previous or previous[particle_id] != rank:
                histories[particle_id].append({"time": time, "rank": rank})
                previous[particle_id] = rank
    return {particle_id: tuple(items) for particle_id, items in histories.items()}


def _norm(values: Sequence[float]) -> float:
    return math.sqrt(sum(float(value) ** 2 for value in values))


def contact_sets(
    trajectory: Sequence[Mapping[str, Any]],
    wall_contact: bool,
    force_tolerance: float,
    overlap_tolerance: float,
) -> tuple[dict[str, Any], ...]:
    """Infer stable pair/wall contact sets from geometry and collision force."""
    rows = []
    for frame in trajectory:
        records = {
            int(record["Particle_ID"]): record for record in frame["records"]
        }
        contacts: set[tuple[int, ...]] = set()
        ids = sorted(records)
        if wall_contact:
            for particle_id in ids:
                record = records[particle_id]
                radius = 0.5 * float(record["Diameter"])
                gap = float(record["Coordinates"][1]) - radius
                force = _norm(record["Collision_force"])
                if gap <= overlap_tolerance or force > force_tolerance:
                    contacts.add((particle_id,))
        else:
            for left_index, left_id in enumerate(ids):
                for right_id in ids[left_index + 1 :]:
                    left = records[left_id]
                    right = records[right_id]
                    distance = math.dist(
                        tuple(float(value) for value in left["Coordinates"]),
                        tuple(float(value) for value in right["Coordinates"]),
                    )
                    overlap = (
                        0.5 * float(left["Diameter"])
                        + 0.5 * float(right["Diameter"])
                        - distance
                    )
                    pair_force = min(
                        _norm(left["Collision_force"]),
                        _norm(right["Collision_force"]),
                    )
                    if overlap >= -overlap_tolerance or (
                        len(ids) == 2 and pair_force > force_tolerance
                    ):
                        contacts.add((left_id, right_id))
        rows.append({
            "time": float(frame["time"]),
            "contacts": tuple(sorted(contacts)),
        })
    return tuple(rows)


def compare_contact_sets(
    serial: Sequence[Mapping[str, Any]],
    dmp: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    """Report the first contact-set difference without hiding later execution."""
    if len(serial) != len(dmp):
        return {
            "passed": False,
            "first_difference": {
                "reason": "frame_count",
                "serial": len(serial),
                "dmp": len(dmp),
            },
        }
    first_difference = None
    for left, right in zip(serial, dmp):
        if left["time"] != right["time"] or left["contacts"] != right["contacts"]:
            first_difference = {
                "time": left["time"],
                "serial_time": left["time"],
                "dmp_time": right["time"],
                "serial_contacts": left["contacts"],
                "dmp_contacts": right["contacts"],
            }
            break
    return {
        "passed": first_difference is None,
        "first_difference": first_difference,
        "active_frame_count": sum(bool(row["contacts"]) for row in serial),
        "first_active_time": next(
            (row["time"] for row in serial if row["contacts"]),
            None,
        ),
    }


def comparison_result_dict(result: Any) -> dict[str, Any]:
    """Convert the comparison NamedTuple tree into JSON-ready mappings."""
    return {
        "passed": result.passed,
        "first_difference": (
            None
            if result.first_difference is None
            else result.first_difference._asdict()
        ),
        "field_summaries": [
            summary._asdict() for summary in result.field_summaries
        ],
        "compared_values": result.compared_values,
    }


def maximum_vector_norm(
    trajectory: Sequence[Mapping[str, Any]],
    field: str,
) -> float:
    """Return the largest vector magnitude in a particle trajectory."""
    return max(
        _norm(_components(record, field))
        for frame in trajectory
        for record in frame["records"]
    )


def audit_serial_dmp_trajectories(
    serial_paths: Sequence[Path],
    dmp_paths: Sequence[Path],
    tolerances: Mapping[str, tuple[float, float]],
    wall_contact: bool,
) -> dict[str, Any]:
    """Compare merged trajectories after an initial-state global-ID remap."""
    serial = read_partitioned_vtp_trajectory(serial_paths)
    dmp_raw = read_partitioned_vtp_trajectory(dmp_paths)
    mapping = initial_physical_id_mapping(
        serial,
        dmp_raw,
        position_tolerance=tolerances["Coordinates"][0],
        velocity_tolerance=tolerances["Velocity"][0],
    )
    dmp = remap_to_serial_ids(dmp_raw, mapping)
    comparisons = compare_trajectories(
        serial,
        dmp,
        tolerances=tolerances,
        time_tolerance=1.0e-15,
    )
    serial_contacts = contact_sets(
        serial,
        wall_contact,
        force_tolerance=tolerances["Collision_force"][0],
        overlap_tolerance=tolerances["Coordinates"][0],
    )
    dmp_contacts = contact_sets(
        dmp,
        wall_contact,
        force_tolerance=tolerances["Collision_force"][0],
        overlap_tolerance=tolerances["Coordinates"][0],
    )
    contact_comparison = compare_contact_sets(serial_contacts, dmp_contacts)
    comparison_evidence = {
        field: comparison_result_dict(result)
        for field, result in comparisons.items()
    }
    return {
        "passed": (
            all(result["passed"] for result in comparison_evidence.values())
            and contact_comparison["passed"]
        ),
        "serial_to_dmp_global_id": mapping,
        "trajectory_comparison": comparison_evidence,
        "contact_set_comparison": contact_comparison,
        "rank_transitions": rank_transitions(dmp),
        "maximums": {
            "serial_angular_speed": maximum_vector_norm(
                serial, "Angular_velocity"
            ),
            "dmp_angular_speed": maximum_vector_norm(dmp, "Angular_velocity"),
            "serial_collision_force": maximum_vector_norm(
                serial, "Collision_force"
            ),
            "dmp_collision_force": maximum_vector_norm(dmp, "Collision_force"),
        },
        "serial_frame_count": len(serial),
        "dmp_frame_count": len(dmp),
    }
