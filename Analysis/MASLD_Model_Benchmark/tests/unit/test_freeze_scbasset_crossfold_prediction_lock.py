"""Tests for the scBasset outcome-free prediction-lock CLI."""

from __future__ import annotations

from pathlib import Path
from unittest import mock
import unittest

from scripts import freeze_scbasset_crossfold_prediction_lock as lock_builder


class ScBassetPredictionLockTests(unittest.TestCase):
    def test_cli_maps_repeated_input_roots(self) -> None:
        with mock.patch(
            "sys.argv",
            [
                "freeze_scbasset_crossfold_prediction_lock.py",
                "--chain-root",
                "/chain",
                "--input-root",
                "/fold1",
                "--input-root",
                "/fold2",
                "--output",
                "/lock",
            ],
        ):
            with mock.patch.object(
                lock_builder,
                "build_lock",
                return_value={"status": "pass"},
            ) as build:
                with mock.patch("builtins.print"):
                    lock_builder.main()
        build.assert_called_once_with(
            chain_root=Path("/chain"),
            input_roots=[Path("/fold1"), Path("/fold2")],
            output=Path("/lock"),
        )


if __name__ == "__main__":
    unittest.main()
