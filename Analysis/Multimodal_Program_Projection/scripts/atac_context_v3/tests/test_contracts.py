#!/usr/bin/env python3

from __future__ import annotations

import math
import importlib.util
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPT_DIR))

from atac_context_v3_lib import (  # noqa: E402
    ContractError,
    bh_adjust,
    candidate_root,
    classify_da,
    classify_genetic,
    write_tsv,
)
from genetics_context import (  # noqa: E402
    aggregate_lineage_mass,
    normalize_against_reference,
    reverse_complement,
)


class ContractTests(unittest.TestCase):
    def test_bh(self) -> None:
        observed = bh_adjust([0.01, 0.04, 0.03, 0.002])
        expected = [0.02, 0.04, 0.04, 0.008]
        for left, right in zip(observed, expected):
            self.assertTrue(math.isclose(left, right, rel_tol=0, abs_tol=1e-12))

    def test_da_states(self) -> None:
        self.assertEqual(classify_da(0.01, 0.02, 1, 2), "supported")
        self.assertEqual(classify_da(0.01, 0.02, 1, -2), "discordant")
        self.assertEqual(classify_da(0.01, 0.20, 1, 2), "source_dependent")
        self.assertEqual(classify_da(0.20, 0.20, 1, 2), "indeterminate")

    def test_genetic_states(self) -> None:
        self.assertEqual(classify_genetic(0.94, 0.8, 0.8, 0.8), "untestable")
        self.assertEqual(classify_genetic(0.99, 0.6, 0.7, 0.7), "replicated_accessible")
        self.assertEqual(classify_genetic(0.99, 0.2, 0.6, 0.2), "source_dependent")
        self.assertEqual(classify_genetic(0.99, 0.1, 0.2, 0.1), "partial")
        self.assertEqual(classify_genetic(0.99, 0.0, 0.0, 0.0), "indeterminate")

    def test_fixture_path_and_no_overwrite(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = candidate_root(Path(tmp) / "fixture", fixture_mode=True)
            self.assertEqual(root.name, "fixture")
            path = root / "x.tsv"
            write_tsv(path, ("a",), [{"a": 1}])
            with self.assertRaises(ContractError):
                write_tsv(path, ("a",), [{"a": 2}])

    def test_snv_complement_and_indel_normalization(self) -> None:
        reference = {"chr1": "AACCGGTTAA"}

        def fetch(chrom: str, start: int, end: int) -> str:
            return reference[chrom][start:end]

        snv = normalize_against_reference("chr1", 3, "G", "T", fetch)
        self.assertIsNotNone(snv)
        self.assertEqual((snv.position_1based, snv.ref, snv.alt), (3, "C", "A"))
        insertion = normalize_against_reference("chr1", 2, "A", "AAC", fetch)
        self.assertIsNotNone(insertion)
        self.assertEqual(insertion.ref, "A")
        self.assertEqual(reverse_complement("ACG"), "CGT")

    def test_multi_lineage_mass_is_not_forced_to_one_lineage(self) -> None:
        variants = [
            {"posterior": 0.6, "mapped": True, "chrom": "chr1", "position_1based": 10},
            {"posterior": 0.4, "mapped": True, "chrom": "chr1", "position_1based": 20},
        ]
        first = aggregate_lineage_mass(
            variants,
            lambda chrom, pos: pos == 10,
            lambda chrom, pos: pos == 10,
        )
        second = aggregate_lineage_mass(
            variants,
            lambda chrom, pos: pos in {10, 20},
            lambda chrom, pos: pos in {10, 20},
        )
        self.assertEqual(first["state"], "replicated_accessible")
        self.assertEqual(second["state"], "replicated_accessible")

    def test_genetics_gate_fails_closed_without_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "candidate"
            root.mkdir()
            script = SCRIPT_DIR / "07_prepare_genetic_replay.py"
            result = subprocess.run(
                [sys.executable, str(script), "--candidate-root", str(root), "--fixture-mode"],
                text=True,
                capture_output=True,
                check=True,
            )
            gate = root / "genetics/BLOCKED_UPSTREAM_RELEASE.tsv"
            self.assertTrue(gate.is_file())
            self.assertIn("blocked", result.stdout.lower())
            self.assertFalse((root / "genetics/replay_plan.tsv").exists())

    def test_replay_runner_patch_is_isolated(self) -> None:
        module_path = SCRIPT_DIR / "07_prepare_genetic_replay.py"
        spec = importlib.util.spec_from_file_location("prepare_replay", module_path)
        module = importlib.util.module_from_spec(spec)
        assert spec.loader is not None
        spec.loader.exec_module(module)
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "source.R"
            destination = Path(tmp) / "runner" / "copy.R"
            source.write_text(
                'OUT_DIR <- file.path(FM_DIR, paste0("results/susie_coloc", OUT_SUFFIX), gwas_name)\n'
                "egenes <- unique(eqtl_all$ENSG)\n"
                'susie_method <- "susie"\n',
                encoding="utf-8",
            )
            original = source.read_bytes()
            module.build_runner(source, destination)
            self.assertEqual(source.read_bytes(), original)
            patched = destination.read_text(encoding="utf-8")
            self.assertIn("ATAC_V3_REPLAY_GENE_FILTER", patched)
            self.assertIn("variant_posteriors.tsv.gz", patched)


if __name__ == "__main__":
    unittest.main()
