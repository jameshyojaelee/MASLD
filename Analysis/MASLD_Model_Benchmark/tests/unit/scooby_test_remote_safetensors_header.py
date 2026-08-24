from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import struct
import unittest


ROOT = Path(__file__).parents[2]
MODULE = ROOT / "scripts" / "scooby_remote_safetensors_header.py"
SPEC = importlib.util.spec_from_file_location("scooby_safetensors_tested", MODULE)
assert SPEC is not None and SPEC.loader is not None
header = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(header)


def _safetensors_bytes() -> bytes:
    payload = {
        "decoder.bias": {
            "dtype": "F32",
            "shape": [2],
            "data_offsets": [0, 8],
        },
        "decoder.weight": {
            "dtype": "F32",
            "shape": [2, 3],
            "data_offsets": [8, 32],
        },
        "__metadata__": {"format": "pt"},
    }
    raw = json.dumps(payload, separators=(",", ":")).encode("utf-8")
    return struct.pack("<Q", len(raw)) + raw + bytes(range(32))


class RemoteSafetensorsHeaderTests(unittest.TestCase):
    def test_bounded_header_inspection_does_not_read_payload(self) -> None:
        reader = header.BytesRangeReader(_safetensors_bytes())
        receipt = header.inspect_reader(reader, expected_parameters=8)
        self.assertEqual(receipt["tensor_count"], 2)
        self.assertEqual(receipt["parameter_count"], 8)
        self.assertEqual(reader.request_count, 2)
        self.assertEqual(
            reader.bytes_downloaded, 8 + receipt["header_size_bytes"]
        )
        self.assertFalse(receipt["checkpoint_payload_downloaded"])

    def test_parameter_count_tamper_fails_closed(self) -> None:
        with self.assertRaisesRegex(
            header.SafetensorsHeaderError, "parameter count"
        ):
            header.inspect_reader(
                header.BytesRangeReader(_safetensors_bytes()),
                expected_parameters=7,
            )

    def test_noncontiguous_spans_fail_closed(self) -> None:
        payload = {
            "tensor": {
                "dtype": "F32",
                "shape": [1],
                "data_offsets": [1, 5],
            }
        }
        raw = json.dumps(payload).encode("utf-8")
        binary = struct.pack("<Q", len(raw)) + raw + b"12345"
        with self.assertRaisesRegex(
            header.SafetensorsHeaderError, "not contiguous"
        ):
            header.inspect_reader(header.BytesRangeReader(binary), 1)


if __name__ == "__main__":
    unittest.main()
