from __future__ import annotations

import unittest
from unittest.mock import patch

from masld_cl.scvi_adapter import assert_runtime_versions


class TestScviRuntime(unittest.TestCase):
    def test_exact_runtime_passes(self):
        config = {"software": {"python": "3.10", "scvi_tools": "1.3.3", "torch": "2.3.1"}}
        with patch("masld_cl.scvi_adapter.platform.python_version", return_value="3.10.99"), patch(
            "masld_cl.scvi_adapter.scvi.__version__", "1.3.3"
        ), patch("masld_cl.scvi_adapter.torch.__version__", "2.3.1+cu121"):
            assert_runtime_versions(config)

    def test_torch_version_mismatch_fails_closed(self):
        config = {"software": {"python": "3.10", "scvi_tools": "1.3.3", "torch": "2.3.1"}}
        with patch("masld_cl.scvi_adapter.platform.python_version", return_value="3.10.99"), patch(
            "masld_cl.scvi_adapter.scvi.__version__", "1.3.3"
        ), patch("masld_cl.scvi_adapter.torch.__version__", "2.4.0"):
            with self.assertRaises(RuntimeError):
                assert_runtime_versions(config)


if __name__ == "__main__":
    unittest.main()
