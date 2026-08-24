from __future__ import annotations

import csv
from hashlib import sha256
import json
from pathlib import Path
import tempfile
import unittest

from scripts import borzoi_enformer_build_fixture as fixture


ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "config" / "borzoi_enformer_admission.json"


def write_fasta(path: Path) -> Path:
    records = {
        "chrA": ("ACGT" * 150000),
        "chrB": ("TGCA" * 150000),
        "chrC": ("GATC" * 150000),
    }
    index_rows = []
    offset = 0
    with path.open("wb") as handle:
        for name, sequence in records.items():
            header = f">{name}\n".encode("ascii")
            handle.write(header)
            offset += len(header)
            sequence_offset = offset
            for start in range(0, len(sequence), 60):
                line = (sequence[start : start + 60] + "\n").encode("ascii")
                handle.write(line)
                offset += len(line)
            index_rows.append(f"{name}\t{len(sequence)}\t{sequence_offset}\t60\t61\n")
    fai = Path(f"{path}.fai")
    fai.write_text("".join(index_rows), encoding="ascii")
    return fai


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


class BorzoiEnformerFixtureTests(unittest.TestCase):
    def test_fixture_is_deterministic_shared_and_outcome_free(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fasta = root / "fixture.fa"
            fai = write_fasta(fasta)
            config = json.loads(CONFIG.read_text(encoding="utf-8"))
            config["shared_fixture"]["reference_fasta_sha256"] = digest(fasta)
            config_path = root / "config.json"
            config_path.write_text(json.dumps(config), encoding="utf-8")
            windows = root / "windows.tsv"
            fields = [
                "contig",
                "output_start",
                "output_end",
                "window_id",
                "genomic_fold",
                "selection_hash",
            ]
            with windows.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
                writer.writeheader()
                for fold, contig in ((2, "chrA"), (3, "chrB"), (4, "chrC")):
                    writer.writerow(
                        {
                            "contig": contig,
                            "output_start": 299500,
                            "output_end": 300500,
                            "window_id": f"window{fold}",
                            "genomic_fold": fold,
                            "selection_hash": str(fold) * 64,
                        }
                    )
            first = root / "first"
            second = root / "second"
            first_summary = fixture.build_fixture(
                config_path=config_path,
                fasta_path=fasta,
                fai_path=fai,
                windows_path=windows,
                output=first,
            )
            second_summary = fixture.build_fixture(
                config_path=config_path,
                fasta_path=fasta,
                fai_path=fai,
                windows_path=windows,
                output=second,
            )
            self.assertEqual(first_summary, second_summary)
            self.assertEqual(first_summary["anchors"], 3)
            self.assertFalse(first_summary["observed_outcomes_loaded"])
            self.assertFalse(first_summary["checkpoint_bytes_loaded"])
            for first_path in sorted(first.iterdir()):
                self.assertEqual(first_path.read_bytes(), (second / first_path.name).read_bytes())
            with (first / "sequence_manifest.tsv").open(
                encoding="utf-8", newline=""
            ) as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 6)
            self.assertEqual({int(row["genomic_fold"]) for row in rows}, {2, 3, 4})
            self.assertEqual(
                {row["model_id"] for row in rows},
                {"borzoi_ensemble", "enformer"},
            )
            self.assertTrue(all(row["allele_effect_sign"] == "ALT_minus_REF" for row in rows))
            self.assertTrue(all(row["ref"] != row["alt"] for row in rows))

    def test_outcome_and_checkpoint_firewalls_fail_closed(self) -> None:
        raw = json.loads(CONFIG.read_text(encoding="utf-8"))
        for section, field, value in (
            ("shared_fixture", "observed_atac_loaded", True),
            ("shared_fixture", "sealed_labels_loaded", True),
            ("admission_disposition", "checkpoint_deserialization_allowed", True),
        ):
            mutated = json.loads(json.dumps(raw))
            mutated[section][field] = value
            with tempfile.TemporaryDirectory() as temporary:
                path = Path(temporary) / "config.json"
                path.write_text(json.dumps(mutated), encoding="utf-8")
                with self.assertRaises(fixture.AdmissionFixtureError):
                    fixture._load_contract(path)

    def test_enformer_native_manifest_has_no_strand_pair_contract(self) -> None:
        raw = json.loads(CONFIG.read_text(encoding="utf-8"))
        enformer = raw["models"]["enformer"]
        self.assertNotIn("strand_pair", enformer["target_manifest"]["fields"])
        self.assertIn("target_axis_identity", enformer["reverse_complement_restore"])
        self.assertIn("separately_named_adaptation", enformer["native_inference"])


if __name__ == "__main__":
    unittest.main()
