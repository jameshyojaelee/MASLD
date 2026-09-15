from __future__ import annotations

import unittest

from scripts.download_gse105127_rna_bundle import GSE105127DownloadError
from scripts.download_gse105127_rna_bundle_resumable import validate_response_headers


class _Headers(dict):
    pass


class _Response:
    def __init__(self, status: int, **headers: str) -> None:
        self.status = status
        self.headers = _Headers(headers)


class GSE105127ResumableDownloadTests(unittest.TestCase):
    def test_exact_initial_and_resume_headers(self) -> None:
        validate_response_headers(_Response(200, **{"Content-Length": "100"}), offset=0, expected_bytes=100)
        validate_response_headers(
            _Response(206, **{"Content-Length": "60", "Content-Range": "bytes 40-99/100"}),
            offset=40,
            expected_bytes=100,
        )

    def test_resume_rejects_ignored_range(self) -> None:
        with self.assertRaisesRegex(GSE105127DownloadError, "not an exact byte range"):
            validate_response_headers(_Response(200, **{"Content-Length": "100"}), offset=40, expected_bytes=100)


if __name__ == "__main__":
    unittest.main()
