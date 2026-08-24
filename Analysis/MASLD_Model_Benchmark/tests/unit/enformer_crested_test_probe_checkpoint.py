from __future__ import annotations

import importlib.util
from pathlib import Path
import tempfile
import unittest


MODULE = Path(__file__).parents[2] / "scripts" / "enformer_crested_probe_checkpoint.py"
SPEC = importlib.util.spec_from_file_location("enformer_crested_probe_tested", MODULE)
assert SPEC is not None and SPEC.loader is not None
probe = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(probe)


class EnformerCREstedProbeTests(unittest.TestCase):
    def test_safe_member_rejects_traversal(self) -> None:
        with self.assertRaisesRegex(probe.CheckpointProbeError, "unsafe"):
            probe._safe_member("../model.keras")

    def test_one_hot_and_reverse_complement(self) -> None:
        encoded = probe._one_hot("ACGTN")
        self.assertEqual(encoded.shape, (5, 4))
        self.assertEqual(encoded.sum(axis=1).tolist(), [1.0, 1.0, 1.0, 1.0, 0.0])
        self.assertEqual("AACGT".translate(probe.COMPLEMENT)[::-1], "ACGTT")

    def test_enformer_fixture_header_contract(self) -> None:
        fixture_id = "trainfold2_ccre_032022"
        self.assertEqual(f"{fixture_id}|enformer|REF", "trainfold2_ccre_032022|enformer|REF")

    def test_fasta_parser_rejects_sequence_before_header(self) -> None:
        import gzip

        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "bad.fa.gz"
            with gzip.open(path, "wt", encoding="ascii") as handle:
                handle.write("ACGT\n")
            with self.assertRaisesRegex(probe.CheckpointProbeError, "precedes"):
                probe._read_fasta(path)

    def test_custom_activation_source_fails_closed(self) -> None:
        probe._validate_gelu_source(probe.GELU_ENF_DEFINITION)
        with self.assertRaisesRegex(probe.CheckpointProbeError, "gelu_enf"):
            probe._validate_gelu_source("def gelu_enf(x): return x")


if __name__ == "__main__":
    unittest.main()
