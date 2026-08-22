"""Unit tests for validation.jkr_tools.io."""

from __future__ import annotations

import struct
import tempfile
import unittest
from pathlib import Path

from validation.jkr_tools.io import (
    VtpFormatError,
    array_headers,
    array_headers_with_types,
    interpolate,
    read_arrays,
    read_monitor,
    read_numeric_csv,
    read_vtp,
    values_only,
)


def _write(path: Path, content: str | bytes) -> None:
    if isinstance(content, bytes):
        path.write_bytes(content)
    else:
        path.write_text(content, encoding="utf-8")


def _minimal_vtp(byte_order: str = "BigEndian", header_type: str = "UInt32") -> bytes:
    endian = ">" if byte_order == "BigEndian" else "<"
    length_format = "I" if header_type == "UInt32" else "Q"
    values = struct.pack(endian + "2f", 1.25, -2.5)
    block = struct.pack(endian + length_format, len(values)) + values
    header = (
        f'<VTKFile type="PolyData" version="0.1" byte_order="{byte_order}" '
        f'header_type="{header_type}">'
        '<PolyData><Piece NumberOfPoints="1"><PointData>'
        '<DataArray type="Float32" Name="force" NumberOfComponents="2" '
        'format="appended" offset="0"/>'
        '</PointData></Piece></PolyData>'
        '<AppendedData encoding="raw">_'
    ).encode("ascii")
    return header + block + b"</AppendedData></VTKFile>"


class IoTests(unittest.TestCase):
    def test_numeric_csv(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "data.csv"
            _write(path, "time,force\n0,1.5\n2,-3\n")
            self.assertEqual(
                read_numeric_csv(path),
                ({"time": 0.0, "force": 1.5}, {"time": 2.0, "force": -3.0}),
            )

    def test_monitor_csv_and_errors(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "monitor.csv"
            _write(path, '# Run type: NEW\n"Time","VEL_X"\n0,1\n1,2\n')
            self.assertEqual(read_monitor(path)[1]["vel_x"], 2.0)
            _write(path, "# no header\n")
            with self.assertRaisesRegex(ValueError, "header.*Time"):
                read_monitor(path)

    def test_interpolate_and_validation(self) -> None:
        rows = ({"time": 0.0, "value": 0.0}, {"time": 2.0, "value": 4.0})
        self.assertEqual(interpolate(rows, -1.0, "value"), 0.0)
        self.assertEqual(interpolate(rows, 1.0, "value"), 2.0)
        self.assertEqual(interpolate(rows, 3.0, "value"), 4.0)
        with self.assertRaisesRegex(ValueError, "strictly increasing"):
            interpolate(
                ({"time": 1.0, "value": 1.0}, {"time": 1.0, "value": 2.0}),
                1.0,
                "value",
            )

    def test_minimal_vtp_and_header_metadata(self) -> None:
        content = _minimal_vtp()
        header = content.split(b"<AppendedData")[0]
        self.assertEqual(array_headers(header), (("force", 2, 0),))
        self.assertEqual(
            array_headers_with_types(header),
            (("force", 2, 0, "Float32"),),
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "sample.vtp"
            _write(path, content)
            parsed = read_vtp(path)
        self.assertEqual(parsed.byte_order, "BigEndian")
        self.assertEqual(parsed.header_type, "UInt32")
        self.assertIsNone(parsed.time)
        self.assertEqual(parsed.arrays["force"].value_type, "Float32")
        self.assertEqual(parsed.arrays["force"].values, (1.25, -2.5))
        time, arrays = read_arrays(content)
        self.assertIsNone(time)
        self.assertEqual(values_only(arrays), {"force": (1.25, -2.5)})

    def test_vtp_rejects_bad_metadata(self) -> None:
        bad = _minimal_vtp().replace(
            b'byte_order="BigEndian"', b'byte_order="Native"'
        )
        with self.assertRaisesRegex(VtpFormatError, "byte_order"):
            read_arrays(bad)


if __name__ == "__main__":
    unittest.main()

