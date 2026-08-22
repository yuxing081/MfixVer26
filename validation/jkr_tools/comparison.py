"""Deterministic, in-memory comparisons for JKR validation data.

The functions in this module deliberately do not know how MFIX output is
stored. Readers should turn output into sequences of mappings and pass those
values here. This keeps alignment and pass/fail rules independent of file
parsing.
"""

from __future__ import annotations

import math
from typing import Any, Iterable, Mapping, NamedTuple, Sequence


class Tolerance(NamedTuple):
    """Absolute and relative tolerance for one scalar field."""

    absolute: float
    relative: float


class FirstDifference(NamedTuple):
    """The first scalar value that exceeds its allowed tolerance."""

    field: str
    index: Any
    time: Any
    left: float
    right: float
    absolute_error: float
    tolerance: float


class FieldSummary(NamedTuple):
    """Largest absolute error observed for one field."""

    field: str
    maximum_absolute_error: float
    time: Any
    index: Any
    left: float
    right: float
    tolerance: float


class ComparisonResult(NamedTuple):
    """Immutable result returned by comparison functions."""

    passed: bool
    first_difference: FirstDifference | None
    field_summaries: tuple[FieldSummary, ...]
    compared_values: int


def make_tolerance(absolute: float, relative: float) -> Tolerance:
    """Create a validated tolerance pair."""

    absolute_value = float(absolute)
    relative_value = float(relative)
    if not math.isfinite(absolute_value) or not math.isfinite(relative_value):
        raise ValueError("tolerances must be finite")
    if absolute_value < 0.0 or relative_value < 0.0:
        raise ValueError("tolerances must be non-negative")
    return Tolerance(absolute_value, relative_value)


def _normalise_tolerance(value: Tolerance | Sequence[float]) -> Tolerance:
    if isinstance(value, Tolerance):
        return make_tolerance(value.absolute, value.relative)
    if len(value) != 2:
        raise ValueError("a tolerance must contain absolute and relative values")
    return make_tolerance(value[0], value[1])


def _absolute_error(left: float, right: float) -> float:
    left_value = float(left)
    right_value = float(right)
    if math.isnan(left_value) or math.isnan(right_value):
        return math.inf
    if left_value == right_value:
        return 0.0
    if math.isinf(left_value) or math.isinf(right_value):
        return math.inf
    return abs(left_value - right_value)


def allowed_error(left: float, right: float, tolerance: Tolerance) -> float:
    """Return absolute + relative * max(|left|, |right|)."""

    checked_tolerance = _normalise_tolerance(tolerance)
    left_value = float(left)
    right_value = float(right)
    return checked_tolerance.absolute + checked_tolerance.relative * max(
        abs(left_value), abs(right_value)
    )


def _time_values(
    count: int,
    left_times: Sequence[Any] | None,
    right_times: Sequence[Any] | None,
) -> tuple[Any, ...]:
    if (left_times is None) != (right_times is None):
        raise ValueError("left_times and right_times must both be supplied or omitted")
    if left_times is None or right_times is None:
        return tuple(range(count))
    if len(left_times) != count or len(right_times) != count:
        raise ValueError("time sequence length does not match value sequence length")
    for index, (left_time, right_time) in enumerate(zip(left_times, right_times)):
        if left_time != right_time:
            raise ValueError(
                f"time mismatch at index {index}: {left_time!r} != {right_time!r}"
            )
    return tuple(left_times)


def compare_field_sequences(
    left_values: Sequence[float],
    right_values: Sequence[float],
    field: str,
    tolerance: Tolerance | Sequence[float],
    left_times: Sequence[Any] | None = None,
    right_times: Sequence[Any] | None = None,
) -> ComparisonResult:
    """Compare two scalar sequences and report first and maximum differences.

    If explicit times are omitted, the sequence index is used as the reported
    time. Output readers should pass actual times.
    """

    if len(left_values) != len(right_values):
        raise ValueError(
            f"sequence length mismatch: {len(left_values)} != {len(right_values)}"
        )
    checked_tolerance = _normalise_tolerance(tolerance)
    times = _time_values(len(left_values), left_times, right_times)
    first_difference: FirstDifference | None = None
    maximum_error = -1.0
    maximum_index: int | None = None
    maximum_time: Any = None
    maximum_left = 0.0
    maximum_right = 0.0
    maximum_allowed = 0.0
    for index, (left, right, time) in enumerate(zip(left_values, right_values, times)):
        error = _absolute_error(left, right)
        allowed = allowed_error(left, right, checked_tolerance)
        if error > maximum_error:
            maximum_error = error
            maximum_index = index
            maximum_time = time
            maximum_left = float(left)
            maximum_right = float(right)
            maximum_allowed = allowed
        if error > allowed and first_difference is None:
            first_difference = FirstDifference(
                field=field,
                index=index,
                time=time,
                left=float(left),
                right=float(right),
                absolute_error=error,
                tolerance=allowed,
            )
    if maximum_index is None:
        summary = FieldSummary(field, 0.0, None, None, 0.0, 0.0, checked_tolerance.absolute)
    else:
        summary = FieldSummary(
            field, maximum_error, maximum_time, maximum_index,
            maximum_left, maximum_right, maximum_allowed,
        )
    return ComparisonResult(
        first_difference is None, first_difference, (summary,), len(left_values)
    )


def compare_fields(
    left_fields: Mapping[str, Sequence[float]],
    right_fields: Mapping[str, Sequence[float]],
    tolerances: Mapping[str, Tolerance | Sequence[float]],
    left_times: Sequence[Any] | None = None,
    right_times: Sequence[Any] | None = None,
) -> ComparisonResult:
    """Compare several named scalar sequences using field-specific tolerances."""

    fields = tuple(tolerances)
    if set(left_fields) != set(fields) or set(right_fields) != set(fields):
        raise ValueError(
            "field names differ: expected "
            f"{sorted(fields)!r}, left={sorted(left_fields)!r}, "
            f"right={sorted(right_fields)!r}"
        )
    results = tuple(
        compare_field_sequences(
            left_fields[field], right_fields[field], field, tolerances[field],
            left_times, right_times,
        )
        for field in fields
    )
    first_candidates = tuple(
        (result.first_difference, order)
        for order, result in enumerate(results)
        if result.first_difference is not None
    )
    first_difference = None
    if first_candidates:
        first_difference, _ = min(
            first_candidates, key=lambda item: (item[0].index, item[1])
        )
    return ComparisonResult(
        first_difference is None,
        first_difference,
        tuple(summary for result in results for summary in result.field_summaries),
        sum(result.compared_values for result in results),
    )


def _require_mapping_field(record: Mapping[str, Any], field: str) -> Any:
    if field not in record:
        raise ValueError(f"record is missing required field {field!r}")
    return record[field]


def _unique_id_map(
    records: Sequence[Mapping[str, Any]], id_field: str
) -> dict[Any, Mapping[str, Any]]:
    records_by_id: dict[Any, Mapping[str, Any]] = {}
    for index, record in enumerate(records):
        particle_id = _require_mapping_field(record, id_field)
        if particle_id in records_by_id:
            raise ValueError(
                f"duplicate {id_field} {particle_id!r} at record index {index}"
            )
        records_by_id[particle_id] = record
    return records_by_id


def align_by_global_particle_id(
    left_records: Sequence[Mapping[str, Any]],
    right_records: Sequence[Mapping[str, Any]],
    id_field: str = "global_particle_id",
) -> tuple[tuple[Mapping[str, Any], Mapping[str, Any]], ...]:
    """Align particle records by global ID and reject duplicate or missing IDs."""

    left_by_id = _unique_id_map(left_records, id_field)
    right_by_id = _unique_id_map(right_records, id_field)
    left_ids = set(left_by_id)
    right_ids = set(right_by_id)
    if left_ids != right_ids:
        missing = sorted(left_ids - right_ids)
        extra = sorted(right_ids - left_ids)
        raise ValueError(f"particle IDs differ: missing={missing!r}, extra={extra!r}")
    return tuple(
        (record, right_by_id[particle_id])
        for particle_id, record in left_by_id.items()
    )


def canonical_contact_pair(first_id: Any, second_id: Any) -> tuple[Any, Any]:
    """Return the order-independent physical pair ID."""

    if first_id == second_id:
        raise ValueError(f"a contact pair cannot contain the same ID twice: {first_id!r}")
    try:
        return (first_id, second_id) if first_id < second_id else (second_id, first_id)
    except TypeError as error:
        raise ValueError("contact IDs must be mutually orderable") from error


def align_by_contact_pair(
    left_records: Sequence[Mapping[str, Any]],
    right_records: Sequence[Mapping[str, Any]],
    first_id_field: str = "global_particle_id_a",
    second_id_field: str = "global_particle_id_b",
) -> tuple[tuple[Mapping[str, Any], Mapping[str, Any]], ...]:
    """Align contact records by canonical (min_gid, max_gid) pair ID."""

    left_by_pair: dict[tuple[Any, Any], Mapping[str, Any]] = {}
    right_by_pair: dict[tuple[Any, Any], Mapping[str, Any]] = {}
    for side, records, destination in (
        ("left", left_records, left_by_pair),
        ("right", right_records, right_by_pair),
    ):
        for index, record in enumerate(records):
            pair = canonical_contact_pair(
                _require_mapping_field(record, first_id_field),
                _require_mapping_field(record, second_id_field),
            )
            if pair in destination:
                raise ValueError(
                    f"duplicate contact pair {pair!r} in {side} records at record index {index}"
                )
            destination[pair] = record
    left_pairs = set(left_by_pair)
    right_pairs = set(right_by_pair)
    if left_pairs != right_pairs:
        missing = sorted(left_pairs - right_pairs)
        extra = sorted(right_pairs - left_pairs)
        raise ValueError(f"contact pairs differ: missing={missing!r}, extra={extra!r}")
    return tuple(
        (record, right_by_pair[pair]) for pair, record in left_by_pair.items()
    )


def _record_fields(
    aligned_records: Iterable[tuple[Mapping[str, Any], Mapping[str, Any]]],
    fields: Sequence[str],
) -> tuple[dict[str, tuple[float, ...]], dict[str, tuple[float, ...]]]:
    left_values = {field: [] for field in fields}
    right_values = {field: [] for field in fields}
    for left_record, right_record in aligned_records:
        for field in fields:
            left_values[field].append(float(_require_mapping_field(left_record, field)))
            right_values[field].append(float(_require_mapping_field(right_record, field)))
    return (
        {field: tuple(values) for field, values in left_values.items()},
        {field: tuple(values) for field, values in right_values.items()},
    )


def compare_aligned_records(
    aligned_records: Sequence[tuple[Mapping[str, Any], Mapping[str, Any]]],
    tolerances: Mapping[str, Tolerance | Sequence[float]],
    time: Any = None,
) -> ComparisonResult:
    """Compare already aligned records at one physical time."""

    fields = tuple(tolerances)
    left_fields, right_fields = _record_fields(aligned_records, fields)
    times = tuple(time for _ in aligned_records)
    return compare_fields(left_fields, right_fields, tolerances, times, times)


def compare_particle_records(
    left_records: Sequence[Mapping[str, Any]],
    right_records: Sequence[Mapping[str, Any]],
    tolerances: Mapping[str, Tolerance | Sequence[float]],
    id_field: str = "global_particle_id",
    time: Any = None,
) -> ComparisonResult:
    """Align and compare one particle frame by global particle ID."""

    aligned_records = align_by_global_particle_id(left_records, right_records, id_field)
    return compare_aligned_records(aligned_records, tolerances, time)


def compare_contact_records(
    left_records: Sequence[Mapping[str, Any]],
    right_records: Sequence[Mapping[str, Any]],
    tolerances: Mapping[str, Tolerance | Sequence[float]],
    first_id_field: str = "global_particle_id_a",
    second_id_field: str = "global_particle_id_b",
    time: Any = None,
) -> ComparisonResult:
    """Align and compare one contact frame by canonical physical pair ID."""

    aligned_records = align_by_contact_pair(
        left_records, right_records, first_id_field, second_id_field
    )
    return compare_aligned_records(aligned_records, tolerances, time)

