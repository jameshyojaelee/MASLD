"""Unit checks for the fail-closed scBasset campaign audit."""

from __future__ import annotations

from pathlib import Path
import unittest

from scripts import audit_scbasset_all5_campaign as audit


ROOT = Path(__file__).parents[2]
WRAPPER = ROOT / "slurm/audit_scbasset_all5_campaign.sbatch"


class ScBassetCampaignAuditTests(unittest.TestCase):
    def test_terminal_disposition_blocks_gpu_after_all_lineages_lose(self) -> None:
        aggregate = {
            "outer_folds": [0, 1, 2, 3, 4],
            "n_donors": 39,
            "mixed_seed_meta_aggregate": True,
            "fixed_seed_five_fold_cv": False,
            "development_gate": {
                "overall_threshold_passed": False,
                "improved_lineages": 0,
            },
            "lineage_relative_deviance_reduction": {
                lineage: -0.1
                for lineage in (
                    "cholangiocyte",
                    "fibroblast",
                    "hepatocyte",
                    "macrophage",
                    "t_cell",
                )
            },
        }
        result = audit.terminal_disposition(aggregate)
        self.assertFalse(result["development_gate_passed"])
        self.assertFalse(result["gpu_production_authorized"])

    def test_terminal_disposition_fails_closed_on_drift(self) -> None:
        with self.assertRaises(audit.ScBassetCampaignAuditError):
            audit.terminal_disposition({"outer_folds": [0, 1, 2, 3, 4]})

    def test_wrapper_is_cpu_only_and_does_not_submit(self) -> None:
        text = WRAPPER.read_text(encoding="utf-8")
        header = "\n".join(
            line for line in text.splitlines() if line.startswith("#SBATCH")
        )
        self.assertIn("--partition=cpu", header)
        self.assertIn("--account=nslab", header)
        self.assertIn("--qos=nslab", header)
        self.assertIn("--cpus-per-task=1", header)
        self.assertIn("--mem=16G", header)
        self.assertNotIn("--gres", header)
        self.assertNotIn("--array", header)
        self.assertFalse(
            any(line.lstrip().startswith("sbatch ") for line in text.splitlines())
        )
        self.assertIn(
            "79ecbd2d473c84ebb75e5e5476dcb4059ac27fe1d061263c733e9d35dc1f9a3f",
            text,
        )


if __name__ == "__main__":
    unittest.main()
