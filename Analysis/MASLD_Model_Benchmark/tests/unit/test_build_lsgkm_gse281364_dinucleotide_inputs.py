from __future__ import annotations

import copy
import gzip
from pathlib import Path
import tempfile
import unittest

from scripts.audit_lsgkm_gse281364_dinucleotide_inputs import (
    DinucleotideInputAuditError,
    audit_pair,
    read_fasta_gzip,
)
from scripts.build_lsgkm_gse281364_dinucleotide_inputs import (
    DinucleotideMaterializationError,
    PeakWindow,
    load_interval_index,
    overlaps,
    select_rectangle,
    validate_config,
    write_fasta_gzip,
)
from scripts.freeze_lsgkm_gse281364_dinucleotide_null_readiness import (
    SEEDS,
    canonical_sequence,
    dinucleotide_counts,
)


ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "config/lsgkm_gse281364_dinucleotide_input_materialization.json"


class LSGKMDinucleotideInputTests(unittest.TestCase):
    def setUp(self) -> None:
        import json

        self.config = json.loads(CONFIG.read_text(encoding="utf-8"))

    def test_config_allows_materialization_and_blocks_every_outcome_or_fit_action(self) -> None:
        validate_config(self.config)
        firewall = self.config["action_firewall"]
        self.assertTrue(firewall["production_input_materialization_authorized"])
        self.assertFalse(firewall["outcome_access_authorized"])
        self.assertFalse(firewall["reporter_count_access_authorized"])
        self.assertFalse(firewall["prediction_value_access_authorized"])
        self.assertFalse(firewall["production_training_authorized"])
        self.assertFalse(firewall["production_prediction_authorized"])

    def test_exact_pair_floor_cannot_be_relaxed(self) -> None:
        changed = copy.deepcopy(self.config)
        changed["terminal_gate"]["exact_pairs_every_split_seed"] = 9999
        with self.assertRaises(DinucleotideMaterializationError):
            validate_config(changed)

    def test_emitted_fasta_parser_and_pair_audit_use_sequence_content(self) -> None:
        positive = ("ACGTTGCAAGTCGATCGGATCCGATGCTAGCTAGGCTA" * 8)[:300]
        candidates = [
            PeakWindow("peak_1", "chr1", 100, 400, positive),
        ]
        selected, negatives, rejected = select_rectangle(
            candidates=candidates,
            scoring_canonical=set(),
            design_id=self.config["design_id"],
            split_id="donor0_genomic0",
            selected_count=1,
            max_attempts=256,
        )
        self.assertEqual(len(selected), 1)
        self.assertEqual(rejected, 0)
        self.assertEqual(set(negatives), set(SEEDS))
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "negative.fa.gz"
            write_fasta_gzip(path, [("negative", negatives[1103][0])])
            emitted = read_fasta_gzip(path)
        self.assertEqual(emitted, [("negative", negatives[1103][0])])
        audit_pair(positive, emitted[0][1])
        self.assertEqual(dinucleotide_counts(positive), dinucleotide_counts(emitted[0][1]))

    def test_emitted_pair_audit_rejects_a_monomer_preserving_dinucleotide_mutation(self) -> None:
        positive = "AC" * 150
        negative = "CA" * 150
        self.assertEqual(sorted(positive), sorted(negative))
        with self.assertRaises(DinucleotideInputAuditError):
            audit_pair(positive, negative)

    def test_selection_rejects_scoring_allele_collision(self) -> None:
        sequence = ("ACGTTGCAAGTCGATCGGATCCGATGCTAGCTAGGCTA" * 8)[:300]
        candidate = PeakWindow("peak_1", "chr1", 100, 400, sequence)
        with self.assertRaises(DinucleotideMaterializationError):
            select_rectangle(
                candidates=[candidate],
                scoring_canonical={canonical_sequence(sequence)},
                design_id=self.config["design_id"],
                split_id="donor0_genomic0",
                selected_count=1,
                max_attempts=4,
            )

    def test_blacklist_overlap_uses_half_open_coordinates(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "blacklist.bed.gz"
            with gzip.open(path, "wt", encoding="utf-8") as handle:
                handle.write("chr1\t100\t200\n")
                handle.write("chr1\t180\t250\n")
            index = load_interval_index(path)
        self.assertFalse(overlaps(index, "chr1", 0, 100))
        self.assertTrue(overlaps(index, "chr1", 99, 101))
        self.assertTrue(overlaps(index, "chr1", 249, 251))
        self.assertFalse(overlaps(index, "chr1", 250, 300))


if __name__ == "__main__":
    unittest.main()
