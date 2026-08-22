"""Functional readers for MFIX validation output.

VTP precision metadata is preserved: MFIX particle VTP is normally Float32,
and Python float storage does not make it a double-precision diagnostic.
"""

from __future__ import annotations

import csv
import math
import re
import struct
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence


class VtpFormatError(ValueError):
    """Raised for inconsistent appended-binary VTP files."""


@dataclass(frozen=True)
class VtpArray:
    components: int
    values: tuple[float, ...]
    value_type: str


@dataclass(frozen=True)
class VtpData:
    time: float | None
    arrays: Mapping[str, VtpArray]
    byte_order: str
    header_type: str


@dataclass(frozen=True)
class _ArrayHeader:
    name: str
    components: int
    offset: int
    value_type: str


def _require_file(path: Path) -> Path:
    if not isinstance(path, Path):
        raise TypeError(f"path must be pathlib.Path, got {type(path).__name__}")
    if not path.is_file():
        raise FileNotFoundError(f"input file does not exist: {path}")
    return path


def _clean_field(value: str) -> str:
    return value.strip().strip('"').strip()


def _parse_numeric_row(values: Sequence[str], path: Path, line_number: int) -> tuple[float, ...]:
    numbers: list[float] = []
    for column, value in enumerate(values, 1):
        cleaned = value.strip()
        try:
            number = float(cleaned)
        except ValueError as error:
            raise ValueError(
                f"invalid numeric value {cleaned!r} in {path} at line "
                f"{line_number}, column {column}"
            ) from error
        if not cleaned or not math.isfinite(number):
            raise ValueError(
                f"non-finite or empty numeric value in {path} at line "
                f"{line_number}, column {column}"
            )
        numbers.append(number)
    return tuple(numbers)


def read_numeric_csv(path: Path) -> tuple[dict[str, float], ...]:
    """Read a comma-separated table with numeric fields."""

    input_path = _require_file(path)
    with input_path.open(newline="", encoding="utf-8") as stream:
        reader = csv.reader(stream)
        try:
            names = tuple(_clean_field(name) for name in next(reader))
        except StopIteration as error:
            raise ValueError(f"numeric CSV is empty: {input_path}") from error
        if not names or any(not name for name in names):
            raise ValueError(f"numeric CSV has an empty column name: {input_path}")
        if len(set(names)) != len(names):
            raise ValueError(f"numeric CSV has duplicate column names: {input_path}")
        rows: list[dict[str, float]] = []
        for line_number, values in enumerate(reader, 2):
            if not values or all(not value.strip() for value in values):
                continue
            if len(values) != len(names):
                raise ValueError(
                    f"wrong column count in {input_path} at line {line_number}: "
                    f"expected {len(names)}, got {len(values)}"
                )
            rows.append(dict(zip(names, _parse_numeric_row(values, input_path, line_number))))
    if not rows:
        raise ValueError(f"numeric CSV has no data rows: {input_path}")
    return tuple(rows)


def read_monitor(path: Path) -> tuple[dict[str, float], ...]:
    """Read an MFIX monitor CSV, skipping comments and preamble lines."""

    input_path = _require_file(path)
    lines = input_path.read_text(encoding="utf-8").splitlines()
    header_index: int | None = None
    names: tuple[str, ...] | None = None
    for index, line in enumerate(lines):
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        fields = tuple(_clean_field(value) for value in next(csv.reader([line])))
        if any(field.casefold() == "time" for field in fields):
            header_index, names = index, tuple(field.lower() for field in fields)
            break
    if header_index is None or names is None:
        raise ValueError(f"MFIX monitor header containing Time not found in {input_path}")
    if any(not name for name in names) or len(set(names)) != len(names):
        raise ValueError(f"MFIX monitor has empty or duplicate column names: {input_path}")
    rows: list[dict[str, float]] = []
    for line_number, line in enumerate(lines[header_index + 1:], header_index + 2):
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        values = tuple(next(csv.reader([line])))
        if len(values) != len(names):
            raise ValueError(
                f"wrong monitor column count in {input_path} at line {line_number}: "
                f"expected {len(names)}, got {len(values)}"
            )
        rows.append(dict(zip(names, _parse_numeric_row(values, input_path, line_number))))
    if not rows:
        raise ValueError(f"No monitor data rows found in {input_path}")
    return tuple(rows)


def interpolate(rows: Sequence[Mapping[str, float]], time: float, field: str) -> float:
    """Linearly interpolate field at time with endpoint clamping."""

    if not rows:
        raise ValueError("cannot interpolate an empty row sequence")
    if not math.isfinite(time):
        raise ValueError(f"interpolation time must be finite, got {time!r}")
    if field not in rows[0]:
        raise KeyError(f"field {field!r} is not present in interpolation rows")
    raw_times = [row.get("time") for row in rows]
    if any(value is None for value in raw_times):
        raise KeyError("interpolation rows must contain a 'time' field")
    times = [float(value) for value in raw_times if value is not None]
    if any(not math.isfinite(value) for value in times):
        raise ValueError("interpolation row times must be finite")
    if any(first >= second for first, second in zip(times, times[1:])):
        raise ValueError("interpolation row times must be strictly increasing")
    if time <= times[0]:
        return float(rows[0][field])
    if time >= times[-1]:
        return float(rows[-1][field])
    upper = next(index for index, value in enumerate(times) if value > time)
    lower = upper - 1
    fraction = (time - times[lower]) / (times[upper] - times[lower])
    return float(rows[lower][field]) + fraction * (float(rows[upper][field]) - float(rows[lower][field]))


_ATTRIBUTE_PATTERN = re.compile(rb"""([A-Za-z_:][\w:.-]*)\s*=\s*(["'])(.*?)\2""")
_DATA_ARRAY_PATTERN = re.compile(rb"<DataArray\b([^>]*)/?>", re.IGNORECASE)


def _attributes(raw_attributes: bytes) -> dict[str, str]:
    return {
        name.decode("ascii").lower(): value.decode("utf-8")
        for name, _quote, value in _ATTRIBUTE_PATTERN.findall(raw_attributes)
    }


def _vtp_metadata(header: bytes) -> tuple[str, str]:
    match = re.search(rb"<VTKFile\b([^>]*)>", header, re.IGNORECASE)
    if match is None:
        raise VtpFormatError("VTP header does not contain a VTKFile element")
    attributes = _attributes(match.group(1))
    byte_order = attributes.get("byte_order")
    header_type = attributes.get("header_type")
    if byte_order not in {"BigEndian", "LittleEndian"}:
        raise VtpFormatError(f"VTP byte_order must be BigEndian or LittleEndian; got {byte_order!r}")
    if header_type not in {"UInt32", "UInt64"}:
        raise VtpFormatError(f"VTP header_type must be UInt32 or UInt64; got {header_type!r}")
    return byte_order, header_type


def _array_headers_with_types(header: bytes) -> tuple[_ArrayHeader, ...]:
    parsed: list[_ArrayHeader] = []
    for match in _DATA_ARRAY_PATTERN.finditer(header):
        attributes = _attributes(match.group(1))
        if attributes.get("format", "").lower() != "appended":
            continue
        name = attributes.get("name")
        value_type = attributes.get("type")
        offset = attributes.get("offset")
        if not name or value_type not in {"Float32", "Float64"} or offset is None:
            raise VtpFormatError(
                "appended VTP DataArray requires Name, type Float32/Float64, and offset"
            )
        try:
            numeric_offset = int(offset)
            components = int(attributes.get("numberofcomponents", "1"))
        except ValueError as error:
            raise VtpFormatError(f"invalid appended VTP attributes for {name!r}") from error
        if numeric_offset < 0 or components < 1:
            raise VtpFormatError(f"invalid offset/components for appended VTP array {name!r}")
        parsed.append(_ArrayHeader(name, components, numeric_offset, value_type))
    if not parsed:
        raise VtpFormatError("VTP header contains no appended-binary DataArray elements")
    names = [item.name for item in parsed]
    if len(set(names)) != len(names):
        raise VtpFormatError("VTP contains duplicate appended DataArray names")
    return tuple(parsed)


def array_headers(header: bytes) -> tuple[tuple[str, int, int], ...]:
    """Return the historical name, components, offset view."""

    return tuple((item.name, item.components, item.offset)
                 for item in _array_headers_with_types(header))


def array_headers_with_types(header: bytes) -> tuple[tuple[str, int, int, str], ...]:
    """Return appended VTP headers including each declared value type."""

    return tuple((item.name, item.components, item.offset, item.value_type)
                 for item in _array_headers_with_types(header))


def _read_vtp_bytes(content: bytes) -> VtpData:
    marker = content.find(b"<AppendedData")
    if marker < 0:
        raise VtpFormatError("VTP file does not contain an AppendedData element")
    header = content[:marker]
    byte_order, header_type = _vtp_metadata(header)
    headers = _array_headers_with_types(header)
    underscore = content.find(b"_", marker)
    if underscore < 0:
        raise VtpFormatError("VTP AppendedData has no '_' payload marker")
    endian = ">" if byte_order == "BigEndian" else "<"
    length_format = "I" if header_type == "UInt32" else "Q"
    length_size = struct.calcsize(length_format)
    value_formats = {"Float32": "f", "Float64": "d"}
    time_match = re.search(rb"<!--\s*Time\s*=\s*([^\s]+)", header, re.IGNORECASE)
    time = None if time_match is None else float(time_match.group(1))
    arrays: dict[str, VtpArray] = {}
    appended_start = underscore + 1
    for item in headers:
        block_start = appended_start + item.offset
        block_end = block_start + length_size
        if block_start < appended_start or block_end > len(content):
            raise VtpFormatError(f"VTP array {item.name!r} length header is outside the file")
        byte_count = struct.unpack(endian + length_format, content[block_start:block_end])[0]
        data_start, data_end = block_end, block_end + byte_count
        if data_end > len(content):
            raise VtpFormatError(f"VTP array {item.name!r} payload is truncated")
        value_size = struct.calcsize(value_formats[item.value_type])
        if byte_count % value_size:
            raise VtpFormatError(f"VTP array {item.name!r} has an invalid byte count")
        value_count = byte_count // value_size
        if value_count % item.components:
            raise VtpFormatError(
                f"VTP array {item.name!r} value count is not divisible by "
                f"NumberOfComponents={item.components}"
            )
        values = struct.unpack(
            endian + f"{value_count}{value_formats[item.value_type]}",
            content[data_start:data_end],
        )
        arrays[item.name] = VtpArray(item.components, tuple(values), item.value_type)
    return VtpData(time, arrays, byte_order, header_type)


def read_vtp(path: Path) -> VtpData:
    """Read one appended-binary VTP file with precision metadata."""

    return _read_vtp_bytes(_require_file(path).read_bytes())


def read_arrays(content: bytes) -> tuple[float | None, dict[str, tuple[int, tuple[float, ...]]]]:
    """Read VTP bytes using the historical component/value view."""

    parsed = _read_vtp_bytes(content)
    return parsed.time, {
        name: (array.components, array.values) for name, array in parsed.arrays.items()
    }


def read_file(path: Path) -> tuple[float | None, dict[str, tuple[int, tuple[float, ...]]]]:
    """Compatibility wrapper for callers that do not need VTP metadata."""

    parsed = read_vtp(path)
    return parsed.time, {
        name: (array.components, array.values) for name, array in parsed.arrays.items()
    }


def values_only(
    arrays: Mapping[str, VtpArray | tuple[int, tuple[float, ...]]],
) -> dict[str, tuple[float, ...]]:
    """Drop component and precision metadata from a parsed array mapping."""

    return {
        name: array.values if isinstance(array, VtpArray) else array[1]
        for name, array in arrays.items()
    }

