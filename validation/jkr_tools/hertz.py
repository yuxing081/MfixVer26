"""Pure-Python Hertz VTP trajectory comparison and contact metrics."""
from __future__ import annotations
import math
from pathlib import Path
from typing import Any, Mapping, Sequence
from .comparison import ComparisonResult, compare_fields
from .io import VtpArray, VtpData, read_vtp

COMPARABLE_FIELDS = ("Coordinates", "Velocity", "Angular_velocity", "Collision_force")
ALIASES = {"Coordinates": ("coordinates", "Coordinates", "Points", "Position"), "Velocity": ("Velocity",), "Angular_velocity": ("Angular_velocity",), "Collision_force": ("Collision_force",), "Diameter": ("Diameter",)}

def _id(value: Any, index: int) -> int:
    number = float(value)
    if not math.isfinite(number) or not number.is_integer():
        raise ValueError(f"Particle_ID at row {index} is not a finite integer: {value!r}")
    return int(number)

def _rows(array: VtpArray, field: str) -> tuple[tuple[float, ...], ...]:
    if array.components < 1 or len(array.values) % array.components:
        raise ValueError(f"VTP field {field!r} has an invalid component layout")
    return tuple(tuple(float(x) for x in array.values[i:i + array.components]) for i in range(0, len(array.values), array.components))

def frame_to_records(frame: VtpData) -> dict[str, Any]:
    if frame.time is None or not math.isfinite(float(frame.time)):
        raise ValueError("VTP frame has no finite Time metadata")
    ids = frame.arrays.get("Particle_ID")
    if ids is None:
        raise ValueError("VTP frame is missing required Particle_ID")
    if ids.components != 1:
        raise ValueError("Particle_ID must be a scalar VTP array")
    particle_ids = tuple(_id(x, i) for i, x in enumerate(ids.values))
    if len(set(particle_ids)) != len(particle_ids):
        raise ValueError("VTP frame contains duplicate Particle_ID values")
    decoded: dict[str, tuple[tuple[float, ...], ...]] = {}
    for field in (*COMPARABLE_FIELDS, "Diameter"):
        array = next((frame.arrays[name] for name in ALIASES[field] if name in frame.arrays), None)
        if array is None:
            continue
        rows = _rows(array, field)
        if len(rows) != len(particle_ids):
            raise ValueError(f"VTP field {field!r} row count does not match Particle_ID count")
        decoded[field] = rows
    records = []
    for i, particle_id in enumerate(particle_ids):
        record: dict[str, Any] = {"Particle_ID": particle_id}
        for field, rows in decoded.items():
            values = rows[i]
            record[field] = values[0] if field == "Diameter" and len(values) == 1 else values
        records.append(record)
    return {"time": float(frame.time), "records": tuple(records)}

def _record_map(frame: Mapping[str, Any], frame_index: int) -> dict[int, Mapping[str, Any]]:
    result: dict[int, Mapping[str, Any]] = {}
    for row_index, record in enumerate(frame["records"]):
        if "Particle_ID" not in record:
            raise ValueError(f"trajectory frame {frame_index} row {row_index} is missing Particle_ID")
        particle_id = _id(record["Particle_ID"], row_index)
        if particle_id in result:
            raise ValueError(f"trajectory frame {frame_index} contains duplicate Particle_ID {particle_id}")
        result[particle_id] = record
    return result

def _validate(frames: Sequence[Mapping[str, Any]]) -> None:
    expected: set[int] | None = None
    previous: float | None = None
    for index, frame in enumerate(frames):
        time = float(frame["time"])
        if not math.isfinite(time) or (previous is not None and time <= previous):
            raise ValueError("trajectory frame times must be finite and strictly increasing")
        previous = time
        ids = set(_record_map(frame, index))
        if expected is None:
            expected = ids
        elif ids != expected:
            raise ValueError(f"particle IDs differ between trajectory frames at index {index}: missing={sorted(expected - ids)!r}, extra={sorted(ids - expected)!r}")
def _drop_duplicate_terminal_frames(frames: Sequence[Mapping[str, Any]]) -> tuple[dict[str, Any], ...]:
    result = list(frames)
    while len(result) > 1 and float(result[-1]["time"]) == float(result[-2]["time"]):
        if _record_map(result[-1], len(result) - 1) != _record_map(result[-2], len(result) - 2):
            raise ValueError("duplicate trajectory time has different Particle_ID records")
        result.pop()
    return tuple(result)

def trajectory_from_vtp_data(frames: Sequence[VtpData]) -> tuple[dict[str, Any], ...]:
    if not frames:
        raise ValueError("trajectory contains no VTP frames")
    result = _drop_duplicate_terminal_frames(tuple(frame_to_records(frame) for frame in frames))
    _validate(result)
    return result

def read_vtp_trajectory(paths: Sequence[Path]) -> tuple[dict[str, Any], ...]:
    if not paths:
        raise ValueError("trajectory contains no VTP paths")
    return trajectory_from_vtp_data(tuple(read_vtp(Path(path)) for path in paths))

def trajectory_from_records(frames: Sequence[Mapping[str, Any]]) -> tuple[dict[str, Any], ...]:
    if not frames:
        raise ValueError("trajectory contains no frames")
    result = []
    for index, frame in enumerate(frames):
        if "time" not in frame or "records" not in frame:
            raise ValueError(f"trajectory frame {index} needs time and records")
        records = []
        for row, source in enumerate(frame["records"]):
            if "Particle_ID" not in source:
                raise ValueError(f"trajectory frame {index} row {row} is missing Particle_ID")
            record = dict(source)
            record["Particle_ID"] = _id(source["Particle_ID"], row)
            records.append(record)
        result.append({"time": float(frame["time"]), "records": tuple(records)})
    output = _drop_duplicate_terminal_frames(tuple(result))
    _validate(output)
    return output

def _coerce(frames: Sequence[Path] | Sequence[VtpData] | Sequence[Mapping[str, Any]]) -> tuple[dict[str, Any], ...]:
    if not frames:
        raise ValueError("trajectory contains no frames")
    first = frames[0]
    if isinstance(first, (str, Path)):
        return read_vtp_trajectory(tuple(Path(x) for x in frames))
    if isinstance(first, VtpData):
        return trajectory_from_vtp_data(tuple(frames))
    return trajectory_from_records(tuple(frames))

def _components(record: Mapping[str, Any], field: str) -> tuple[float, ...]:
    if field not in record:
        raise KeyError(field)
    value = record[field]
    values = tuple(float(x) for x in value) if isinstance(value, (tuple, list)) else (float(value),)
    if not values or any(not math.isfinite(x) for x in values):
        raise ValueError(f"field {field!r} contains non-finite or empty values")
    return values

def compare_trajectories(left_frames: Sequence[Path] | Sequence[VtpData] | Sequence[Mapping[str, Any]], right_frames: Sequence[Path] | Sequence[VtpData] | Sequence[Mapping[str, Any]], tolerances: Mapping[str, Any] | None = None, time_tolerance: float = 0.0) -> dict[str, ComparisonResult]:
    tolerance = float(time_tolerance)
    if not math.isfinite(tolerance) or tolerance < 0:
        raise ValueError("time_tolerance must be finite and non-negative")
    left, right = _coerce(left_frames), _coerce(right_frames)
    if len(left) != len(right):
        raise ValueError(f"trajectory frame count mismatch: {len(left)} != {len(right)}")
    left_values: dict[str, list[float]] = {}
    right_values: dict[str, list[float]] = {}
    times_by_key: dict[str, list[float]] = {}
    counts: dict[str, int] = {}
    for index, (lframe, rframe) in enumerate(zip(left, right)):
        lt, rt = float(lframe["time"]), float(rframe["time"])
        if abs(lt - rt) > tolerance:
            raise ValueError(f"time mismatch at frame {index}: {lt!r} != {rt!r}")
        lm, rm = _record_map(lframe, index), _record_map(rframe, index)
        if set(lm) != set(rm):
            raise ValueError(f"particle IDs differ at frame {index}: missing={sorted(set(lm) - set(rm))!r}, extra={sorted(set(rm) - set(lm))!r}")
        for particle_id in sorted(lm):
            lf, rf = set(lm[particle_id]) & set(COMPARABLE_FIELDS), set(rm[particle_id]) & set(COMPARABLE_FIELDS)
            if lf != rf:
                raise ValueError(f"trajectory fields differ for Particle_ID {particle_id} at frame {index}")
            for field in sorted(lf):
                lv, rv = _components(lm[particle_id], field), _components(rm[particle_id], field)
                if len(lv) != len(rv):
                    raise ValueError(f"component count mismatch for {field!r}, Particle_ID {particle_id}")
                if field in counts and counts[field] != len(lv):
                    raise ValueError(f"component count for {field!r} changes across frames")
                counts[field] = len(lv)
                for component, (a, b) in enumerate(zip(lv, rv)):
                    key = f"{field}[{component}]"
                    left_values.setdefault(key, []).append(a)
                    right_values.setdefault(key, []).append(b)
                    times_by_key.setdefault(key, []).append(lt)
    output = {}
    for field in COMPARABLE_FIELDS:
        keys = tuple(key for key in left_values if key.startswith(f"{field}["))
        if not keys:
            continue
        field_tolerances = {}
        for key in keys:
            value = None if tolerances is None else tolerances.get(field, tolerances.get(key))
            field_tolerances[key] = (0.0, 0.0) if value is None else value
        output[field] = compare_fields({key: tuple(left_values[key]) for key in keys}, {key: tuple(right_values[key]) for key in keys}, field_tolerances, times_by_key[keys[0]], times_by_key[keys[0]])
    return output

def compare_vtp_trajectories(left_paths: Sequence[Path], right_paths: Sequence[Path], tolerances: Mapping[str, Any] | None = None, time_tolerance: float = 0.0) -> dict[str, ComparisonResult]:
    return compare_trajectories(left_paths, right_paths, tolerances, time_tolerance)

def _dot(a: Sequence[float], b: Sequence[float]) -> float:
    if len(a) != 3 or len(b) != 3:
        raise ValueError("contact metric vectors must have three components")
    return sum(float(x) * float(y) for x, y in zip(a, b))

def _unit(vector: Sequence[float]) -> tuple[float, float, float]:
    if len(vector) != 3:
        raise ValueError("wall_normal must have three components")
    norm = math.sqrt(_dot(vector, vector))
    if not math.isfinite(norm) or norm == 0:
        raise ValueError("wall_normal must be a finite, non-zero three-vector")
    return tuple(float(x) / norm for x in vector)  # type: ignore[return-value]

def extract_contact_metrics(frames: Sequence[Path] | Sequence[VtpData] | Sequence[Mapping[str, Any]], contact_kind: str = "pair", particle_ids: Sequence[int] | None = None, particle_radius: float | None = None, wall_position: float = 0.0, wall_normal: Sequence[float] = (0.0, 1.0, 0.0), overlap_tolerance: float = 0.0, force_tolerance: float = 0.0, rebound_velocity_tolerance: float = 0.0) -> dict[str, float | bool]:
    trajectory = _coerce(frames)
    kind = contact_kind.casefold()
    if kind not in {"pair", "wall"}:
        raise ValueError("contact_kind must be 'pair' or 'wall'")
    if overlap_tolerance < 0 or force_tolerance < 0 or rebound_velocity_tolerance < 0:
        raise ValueError("contact tolerances must be non-negative")
    ids = sorted(_record_map(trajectory[0], 0))
    selected = tuple(ids if particle_ids is None else particle_ids)
    if (kind == "pair" and len(selected) != 2) or (kind == "wall" and len(selected) != 1):
        raise ValueError(f"{kind} contact metrics require exactly {'two' if kind == 'pair' else 'one'} particle IDs")
    selected = tuple(int(x) for x in selected)
    normal = _unit(wall_normal)
    overlaps, speeds, active = [], [], []
    for index, frame in enumerate(trajectory):
        records = _record_map(frame, index)
        try:
            chosen = tuple(records[pid] for pid in selected)
        except KeyError as error:
            raise ValueError(f"Particle_ID {error.args[0]!r} is missing from frame {index}") from error
        positions = tuple(_components(record, "Coordinates") for record in chosen)
        velocities = tuple(_components(record, "Velocity") for record in chosen)
        if kind == "pair":
            displacement = tuple(b - a for a, b in zip(positions[0], positions[1]))
            distance = math.sqrt(_dot(displacement, displacement))
            if distance == 0:
                raise ValueError(f"pair particle centers coincide at frame {index}")
            radii = [float(record["Diameter"]) / 2 for record in chosen if "Diameter" in record]
            if len(radii) != 2:
                if particle_radius is None:
                    raise ValueError("pair metrics need Diameter fields or particle_radius")
                radii = [float(particle_radius)] * 2
            direction = tuple(x / distance for x in displacement)
            overlap = sum(radii) - distance
            speed = _dot(tuple(a - b for a, b in zip(velocities[1], velocities[0])), direction)
        else:
            radius = float(particle_radius) if particle_radius is not None else None
            if radius is None and "Diameter" in chosen[0]:
                radius = float(chosen[0]["Diameter"]) / 2
            if radius is None or not math.isfinite(radius) or radius <= 0:
                raise ValueError("wall metrics need a positive particle_radius or Diameter field")
            overlap = radius - (_dot(positions[0], normal) - wall_position)
            speed = _dot(velocities[0], normal)
        force = _components(chosen[0], "Collision_force") if "Collision_force" in chosen[0] else ()
        force_active = bool(force) and math.sqrt(_dot(force, force)) > force_tolerance
        overlaps.append(overlap)
        speeds.append(speed)
        active.append(overlap > overlap_tolerance or force_active)
    active_indices = tuple(i for i, value in enumerate(active) if value)
    if active_indices:
        first, last = active_indices[0], active_indices[-1]
        duration = float(trajectory[last]["time"]) - float(trajectory[first]["time"])
        impact = max(0.0, -speeds[0])
        rebound = max((max(0.0, speed) for speed in speeds[last + 1:]), default=0.0)
    else:
        duration, impact = 0.0, max(0.0, -speeds[0])
        rebound = max((max(0.0, speed) for speed in speeds[1:]), default=0.0)
    return {"maximum_overlap": max(overlaps), "minimum_overlap": min(overlaps), "final_overlap": overlaps[-1], "final_normal_velocity": speeds[-1], "contact_active_at_end": active[-1], "contact_duration": duration, "impact_velocity": impact, "rebound_normal_velocity": rebound, "rebound_velocity": rebound, "effective_restitution": rebound / impact if impact > 0 else 0.0, "captured": rebound <= rebound_velocity_tolerance}

def compare_case_trajectories(left_directory: Path, right_directory: Path, pattern: str = "*.vtp", tolerances: Mapping[str, Any] | None = None, time_tolerance: float = 0.0) -> dict[str, Any]:
    left = tuple(sorted(Path(left_directory).glob(pattern)))
    right = tuple(sorted(Path(right_directory).glob(pattern)))
    result = compare_vtp_trajectories(left, right, tolerances, time_tolerance)
    return {"passed": all(item.passed for item in result.values()), "fields": {field: {"passed": item.passed, "compared_values": item.compared_values, "maximum_absolute_error": max((s.maximum_absolute_error for s in item.field_summaries), default=0.0)} for field, item in result.items()}}
