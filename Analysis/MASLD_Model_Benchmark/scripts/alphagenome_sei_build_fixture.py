#!/usr/bin/env python3
"""Build outcome-free AlphaGenome and Sei sequence/variant fixtures."""

from __future__ import annotations

import argparse
import csv
import gzip
from hashlib import sha256
import io
import json
from pathlib import Path
from typing import Iterable, Mapping, Sequence


ALT_CYCLE = {"A": "C", "C": "G", "G": "T", "T": "A"}
COMPLEMENT = str.maketrans("ACGT", "TGCA")
MODEL_LENGTHS = {"alphagenome_all_folds": 1_048_576, "sei": 4_096}
MANIFEST_FIELDS = (
    "fixture_id",
    "model_id",
    "genomic_fold",
    "contig",
    "anchor0",
    "source_window_id",
    "source_selection_hash",
    "input_start",
    "input_end",
    "input_length",
    "variant_pos_1based",
    "variant_index0",
    "ref",
    "alt",
    "reference_sequence_sha256",
    "alternative_sequence_sha256",
    "reverse_complement_reference_sha256",
    "reverse_complement_alternative_sha256",
    "reference_one_hot_uint8_sha256",
    "alternative_one_hot_uint8_sha256",
    "allele_effect_sign",
)


class AdmissionFixtureError(ValueError):
    """Raised when a sequence fixture does not meet its frozen requirements."""


class IndexedFasta:
    """Minimal read-only FAI-backed FASTA accessor with no dependency."""

    def __init__(self, fasta_path: Path, fai_path: Path):
        self.path = fasta_path
        self.index: dict[str, tuple[int, int, int, int]] = {}
        for line_number, line in enumerate(
            fai_path.read_text(encoding="utf-8").splitlines(), start=1
        ):
            fields = line.split("\t")
            if len(fields) < 5:
                raise AdmissionFixtureError(f"FAI row {line_number} differs")
            name = fields[0]
            if name in self.index:
                raise AdmissionFixtureError("FAI contains duplicate contigs")
            try:
                values = tuple(int(value) for value in fields[1:5])
            except ValueError as exc:
                raise AdmissionFixtureError("FAI numeric field differs") from exc
            length, offset, line_bases, line_width = values
            if min(values) < 1 or line_width < line_bases:
                raise AdmissionFixtureError("FAI geometry differs")
            self.index[name] = (length, offset, line_bases, line_width)
        if not self.index:
            raise AdmissionFixtureError("FAI is empty")

    def fetch(self, contig: str, start: int, end: int) -> str:
        if contig not in self.index:
            raise AdmissionFixtureError(f"contig is absent from FAI: {contig}")
        length, offset, line_bases, line_width = self.index[contig]
        if start < 0 or end <= start or end > length:
            raise AdmissionFixtureError("FASTA interval is outside the contig")
        pieces: list[bytes] = []
        position = start
        with self.path.open("rb") as handle:
            while position < end:
                line_offset = position % line_bases
                take = min(end - position, line_bases - line_offset)
                byte_offset = offset + (position // line_bases) * line_width + line_offset
                handle.seek(byte_offset)
                block = handle.read(take)
                if len(block) != take:
                    raise AdmissionFixtureError("FASTA ended inside an indexed interval")
                pieces.append(block)
                position += take
        try:
            return b"".join(pieces).decode("ascii").upper()
        except UnicodeDecodeError as exc:
            raise AdmissionFixtureError("FASTA sequence is not ASCII") from exc

    def contig_length(self, contig: str) -> int:
        if contig not in self.index:
            raise AdmissionFixtureError(f"contig is absent from FAI: {contig}")
        return self.index[contig][0]


def _digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def _sequence_digest(sequence: str) -> str:
    return sha256(sequence.encode("ascii")).hexdigest()


def _one_hot_digest(sequence: str) -> str:
    channels = {
        "A": b"\x01\x00\x00\x00",
        "C": b"\x00\x01\x00\x00",
        "G": b"\x00\x00\x01\x00",
        "T": b"\x00\x00\x00\x01",
        "N": b"\x00\x00\x00\x00",
    }
    value = sha256()
    try:
        for base in sequence:
            value.update(channels[base])
    except KeyError as exc:
        raise AdmissionFixtureError(f"unsupported sequence base: {exc.args[0]}") from exc
    return value.hexdigest()


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def _write_fasta_gz(path: Path, records: Iterable[tuple[str, str]]) -> None:
    with path.open("xb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="ascii", newline="\n") as handle:
                for name, sequence in records:
                    handle.write(f">{name}\n")
                    for start in range(0, len(sequence), 80):
                        handle.write(sequence[start : start + 80] + "\n")


def _load_contract(path: Path) -> Mapping[str, object]:
    raw = json.loads(path.read_text(encoding="utf-8"))
    if raw.get("schema_version") != "masld-bench-alphagenome-sei-admission-v1":
        raise AdmissionFixtureError("admission config schema differs")
    shared = raw.get("shared_fixture")
    disposition = raw.get("admission_disposition")
    models = raw.get("models")
    if not isinstance(shared, dict) or not isinstance(disposition, dict) or not isinstance(models, dict):
        raise AdmissionFixtureError("admission config fields differ")
    if (
        shared.get("eligible_genomic_folds") != [2, 3, 4]
        or shared.get("anchors_per_fold") != 1
        or shared.get("anchor_source_is_outcome_independent") is not True
        or shared.get("model_checkpoints_loaded") is not False
        or any(
            shared.get(field) is not False
            for field in (
                "observed_atac_loaded",
                "rna_loaded",
                "variant_or_regulatory_outcomes_loaded",
                "sealed_labels_loaded",
            )
        )
        or disposition.get("observed_outcomes_allowed") is not False
        or disposition.get("alphagenome_gpu_forward_allowed") is not False
        or disposition.get("sei_gpu_forward_allowed") is not False
    ):
        raise AdmissionFixtureError("outcome or checkpoint firewall differs")
    if set(models) != set(MODEL_LENGTHS) or any(
        models[model_id].get("input_length_bp") != length
        for model_id, length in MODEL_LENGTHS.items()
    ):
        raise AdmissionFixtureError("model sequence geometry differs")
    return raw


def _read_windows(path: Path) -> list[dict[str, str]]:
    required = {
        "contig",
        "output_start",
        "output_end",
        "window_id",
        "genomic_fold",
        "selection_hash",
    }
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not required.issubset(reader.fieldnames or []):
            raise AdmissionFixtureError("anchor manifest fields differ")
        rows = [dict(row) for row in reader]
    if not rows:
        raise AdmissionFixtureError("anchor manifest is empty")
    return rows


def _select_anchors(
    rows: Sequence[Mapping[str, str]], reference: IndexedFasta
) -> list[dict[str, object]]:
    selected: dict[int, dict[str, object]] = {}
    half = max(MODEL_LENGTHS.values()) // 2
    for row in rows:
        try:
            fold = int(row["genomic_fold"])
            output_start = int(row["output_start"])
            output_end = int(row["output_end"])
        except (KeyError, ValueError) as exc:
            raise AdmissionFixtureError("anchor manifest coordinates differ") from exc
        if fold not in (2, 3, 4) or fold in selected:
            continue
        contig = row["contig"]
        anchor0 = (output_start + output_end) // 2
        start, end = anchor0 - half, anchor0 + half
        if start < 0 or end > reference.contig_length(contig):
            continue
        sequence = reference.fetch(contig, start, end)
        if len(sequence) == 2 * half and set(sequence) <= set("ACGT"):
            selected[fold] = {
                "genomic_fold": fold,
                "contig": contig,
                "anchor0": anchor0,
                "source_window_id": row["window_id"],
                "source_selection_hash": row["selection_hash"],
            }
        if len(selected) == 3:
            break
    if set(selected) != {2, 3, 4}:
        raise AdmissionFixtureError("could not select one complete ACGT anchor per fold")
    return [selected[fold] for fold in (2, 3, 4)]


def build_fixture(
    *,
    config_path: Path,
    fasta_path: Path,
    fai_path: Path,
    windows_path: Path,
    output: Path,
) -> dict[str, object]:
    if output.exists():
        raise AdmissionFixtureError("fixture output already exists")
    for path in (config_path, fasta_path, fai_path, windows_path):
        if path.is_symlink() or not path.is_file():
            raise AdmissionFixtureError("fixture input is not a regular file")
    contract = _load_contract(config_path)
    shared = contract["shared_fixture"]
    if _digest(fasta_path) != shared["reference_fasta_sha256"]:
        raise AdmissionFixtureError("reference FASTA checksum differs")
    reference = IndexedFasta(fasta_path, fai_path)
    anchors = _select_anchors(_read_windows(windows_path), reference)
    output.mkdir(mode=0o750)

    manifest_rows: list[dict[str, object]] = []
    fasta_records: dict[str, list[tuple[str, str]]] = {
        model_id: [] for model_id in MODEL_LENGTHS
    }
    for anchor in anchors:
        fixture_id = f"trainfold{anchor['genomic_fold']}_{anchor['source_window_id']}"
        for model_id, input_length in MODEL_LENGTHS.items():
            half = input_length // 2
            start = int(anchor["anchor0"]) - half
            end = start + input_length
            reference_sequence = reference.fetch(str(anchor["contig"]), start, end)
            if len(reference_sequence) != input_length or set(reference_sequence) - set("ACGT"):
                raise AdmissionFixtureError("selected sequence geometry or alphabet differs")
            variant_index = int(anchor["anchor0"]) - start
            if variant_index != half:
                raise AdmissionFixtureError("variant is not at the exact sequence center")
            ref = reference_sequence[variant_index]
            alt = ALT_CYCLE[ref]
            alternative_sequence = (
                reference_sequence[:variant_index]
                + alt
                + reference_sequence[variant_index + 1 :]
            )
            if sum(a != b for a, b in zip(reference_sequence, alternative_sequence)) != 1:
                raise AdmissionFixtureError("alternative allele sequence differs")
            rc_ref = _reverse_complement(reference_sequence)
            rc_alt = _reverse_complement(alternative_sequence)
            for allele, sequence in (
                ("REF", reference_sequence),
                ("ALT", alternative_sequence),
                ("REF_RC", rc_ref),
                ("ALT_RC", rc_alt),
            ):
                fasta_records[model_id].append(
                    (f"{fixture_id}|{model_id}|{allele}", sequence)
                )
            manifest_rows.append(
                {
                    "fixture_id": fixture_id,
                    "model_id": model_id,
                    "genomic_fold": anchor["genomic_fold"],
                    "contig": anchor["contig"],
                    "anchor0": anchor["anchor0"],
                    "source_window_id": anchor["source_window_id"],
                    "source_selection_hash": anchor["source_selection_hash"],
                    "input_start": start,
                    "input_end": end,
                    "input_length": input_length,
                    "variant_pos_1based": int(anchor["anchor0"]) + 1,
                    "variant_index0": variant_index,
                    "ref": ref,
                    "alt": alt,
                    "reference_sequence_sha256": _sequence_digest(reference_sequence),
                    "alternative_sequence_sha256": _sequence_digest(alternative_sequence),
                    "reverse_complement_reference_sha256": _sequence_digest(rc_ref),
                    "reverse_complement_alternative_sha256": _sequence_digest(rc_alt),
                    "reference_one_hot_uint8_sha256": _one_hot_digest(reference_sequence),
                    "alternative_one_hot_uint8_sha256": _one_hot_digest(alternative_sequence),
                    "allele_effect_sign": "ALT_minus_REF",
                }
            )

    for model_id, records in fasta_records.items():
        _write_fasta_gz(output / f"{model_id}.alleles.fa.gz", records)
    with (output / "sequence_manifest.tsv").open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest_rows)

    firewall = {
        "schema_version": "masld-bench-alphagenome-sei-outcome-firewall-v1",
        "allowed_inputs": [
            "GRCh38.p14 sequence FASTA and FAI",
            "outcome-independent fixed cCRE anchor manifest",
            "frozen admission configuration",
        ],
        "checkpoint_bytes_loaded": False,
        "observed_atac_loaded": False,
        "rna_loaded": False,
        "variant_or_regulatory_outcomes_loaded": False,
        "sealed_labels_loaded": False,
        "fixture_genomic_roles": ["train"],
        "eligible_genomic_folds": [2, 3, 4],
        "biological_outcomes_calculated": False,
    }
    (output / "outcome_firewall.json").write_text(
        json.dumps(firewall, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    summary = {
        "schema_version": "masld-bench-alphagenome-sei-sequence-fixture-v1",
        "status": "pass",
        "reference_build": "GRCh38.p14",
        "reference_fasta_sha256": shared["reference_fasta_sha256"],
        "anchor_manifest_sha256": _digest(windows_path),
        "anchors": len(anchors),
        "models": list(MODEL_LENGTHS),
        "sequence_rows": len(manifest_rows),
        "fasta_records_per_model": 12,
        "eligible_genomic_folds": [2, 3, 4],
        "allele_effect_sign": "ALT_minus_REF",
        "observed_outcomes_loaded": False,
        "checkpoint_bytes_loaded": False,
        "runtime_forward_executed": False,
        "model_specific_output_contracts": {
            "alphagenome_all_folds": {
                "input_length_bp": 1_048_576,
                "native_human_output_slots": 6_126,
                "native_task": "fixed_multimodal_genomic_tracks",
                "four_lineage_signed_RNA_roster_complete": False,
            },
            "sei": {
                "input_length_bp": 4_096,
                "forward_shape": [21_907],
                "sequence_class_shape": [40],
                "native_task": "center_position_binary_regulatory_probability",
                "reverse_complement_rule": contract["models"]["sei"][
                    "reverse_complement_rule"
                ],
            },
        },
    }
    (output / "summary.json").write_text(
        json.dumps(summary, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--fai", type=Path, required=True)
    parser.add_argument("--windows", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    summary = build_fixture(
        config_path=args.config,
        fasta_path=args.fasta,
        fai_path=args.fai,
        windows_path=args.windows,
        output=args.output,
    )
    print(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()
