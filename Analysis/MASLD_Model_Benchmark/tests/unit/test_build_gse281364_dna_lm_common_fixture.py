from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import unittest

from scripts.build_gse281364_dna_lm_common_fixture import (
    CommonFixtureError,
    reverse_complement,
    validate_config,
    validate_pair,
)


ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "config/gse281364_dna_lm_common_fixture.json"


def sequence_digest(sequence: str) -> str:
    return sha256(sequence.encode("ascii")).hexdigest()


class GSE281364DnaLmCommonFixtureTests(unittest.TestCase):
    def setUp(self) -> None:
        ref = "A" * 2_048 + "C" + "G" * 2_047
        alt = ref[:2_048] + "T" + ref[2_049:]
        self.sequences = {
            "REF": ref,
            "ALT": alt,
            "REF_RC": reverse_complement(ref),
            "ALT_RC": reverse_complement(alt),
        }
        self.row = {
            "ref": "C",
            "alt": "T",
            "reference_sequence_sha256": sequence_digest(ref),
            "alternative_sequence_sha256": sequence_digest(alt),
            "reference_reverse_complement_sha256": sequence_digest(
                self.sequences["REF_RC"]
            ),
            "alternative_reverse_complement_sha256": sequence_digest(
                self.sequences["ALT_RC"]
            ),
        }

    def test_ref_alt_and_reverse_complement_pair_geometry(self) -> None:
        validate_pair(
            self.row,
            self.sequences,
            input_length=4_096,
            forward_index=2_048,
            reverse_index=2_047,
        )

    def test_tampered_reverse_complement_is_rejected(self) -> None:
        tampered = dict(self.sequences)
        tampered["ALT_RC"] = "A" + tampered["ALT_RC"][1:]
        with self.assertRaisesRegex(CommonFixtureError, "reverse-complement"):
            validate_pair(
                self.row,
                tampered,
                input_length=4_096,
                forward_index=2_048,
                reverse_index=2_047,
            )

    def test_configuration_freezes_models_folds_and_identical_heads(self) -> None:
        config = json.loads(CONFIG.read_text(encoding="utf-8"))
        lane = validate_config(config)
        self.assertEqual(
            lane["models"], ["dnabert2", "nucleotide_transformer", "hyenadna"]
        )
        self.assertEqual(lane["identical_downstream_heads"], ["linear", "two_layer"])
        self.assertIn("outer_locus_sequence_group", lane["head_training"])


if __name__ == "__main__":
    unittest.main()
