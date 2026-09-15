#!/usr/bin/env python3
"""Quantify one participant-contained GSE105127 RNA bundle with RSEM/STAR."""

from __future__ import annotations

import argparse
import csv
import gzip
from hashlib import sha256
import json
import math
from pathlib import Path
import shutil
import subprocess

from masld_bench.artifacts import freeze_tree, publish_directory_noreplace, verify_frozen_tree


class GSE105127QuantificationError(RuntimeError):
    """Raised when a raw RNA quantification differs."""


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_rows(path: Path, bundle_id: int) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        forbidden = {"phenotype", "label", "outcome", "disease", "fibrosis", "nas", "sex", "age", "bmi", "outer_fold"}
        if forbidden & set(reader.fieldnames or ()):
            raise GSE105127QuantificationError("RNA plan violates the label firewall")
        rows = [row for row in reader if int(row["bundle_id"]) == bundle_id]
    if not rows:
        raise GSE105127QuantificationError("RNA quantification bundle is empty")
    groups = {row["participant_group_id"] for row in rows}
    if any(sum(member["participant_group_id"] == group for member in rows) != 3 for group in groups):
        raise GSE105127QuantificationError("participant zones are split inside RNA quantification bundle")
    return rows


def audit_gene_results(path: Path) -> dict[str, object]:
    features = 0
    fractional_expected = 0
    expected_sum = 0.0
    tpm_sum = 0.0
    seen: set[str] = set()
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"gene_id", "transcript_id(s)", "length", "effective_length", "expected_count", "TPM", "FPKM"}
        if not required <= set(reader.fieldnames or ()):
            raise GSE105127QuantificationError("RSEM gene-results schema differs")
        for row in reader:
            gene = row["gene_id"]
            if not gene or gene in seen:
                raise GSE105127QuantificationError("RSEM gene axis is empty or duplicated")
            seen.add(gene)
            try:
                values = [float(row[key]) for key in ("length", "effective_length", "expected_count", "TPM", "FPKM")]
            except ValueError as error:
                raise GSE105127QuantificationError("RSEM measurement is nonnumeric") from error
            if not all(math.isfinite(value) and value >= 0 for value in values):
                raise GSE105127QuantificationError("RSEM measurement is invalid")
            expected = values[2]
            expected_sum += expected
            tpm_sum += values[3]
            fractional_expected += not expected.is_integer()
            features += 1
    # RSEM text output can round each gene TPM independently. The aggregate
    # tolerance is therefore tied to half a unit in the final printed decimal,
    # not an arbitrarily exact one-million check. Integer-only expected counts
    # in one library are eligible; the cohort consolidation separately checks
    # that continuous values survive anywhere they are present.
    tpm_tolerance = 0.005 * features + 1.0
    if features == 0 or abs(tpm_sum - 1_000_000.0) > tpm_tolerance:
        raise GSE105127QuantificationError("RSEM feature/TPM census differs")
    return {
        "genes": features,
        "fractional_expected_count_rows": fractional_expected,
        "expected_count_sum": expected_sum,
        "tpm_sum": tpm_sum,
    }


def quantify_one(
    row: dict[str, str], raw_root: Path, reference_root: Path,
    output: Path, attempt: Path, threads: int, seed: int,
) -> dict[str, object]:
    final = output / "participants" / row["row_id"]
    if (final / "COMPLETE").is_file():
        verify_frozen_tree(final)
        receipt = json.loads((final / "receipt.json").read_text(encoding="utf-8"))
        expected = {
            "row_id": row["row_id"],
            "participant_group_id": row["participant_group_id"],
            "zone": row["zone"],
            "seed": seed,
        }
        if any(receipt.get(key) != value for key, value in expected.items()):
            raise GSE105127QuantificationError("replayed RNA quantification differs from plan")
        return receipt
    if final.exists():
        raise GSE105127QuantificationError("incomplete final RNA quantification exists")
    raw = raw_root / "participants" / row["row_id"]
    verify_frozen_tree(raw)
    raw_receipt = json.loads((raw / "receipt.json").read_text(encoding="utf-8"))
    fastq = raw / "reads.fastq.gz"
    if (
        raw_receipt.get("status") != "complete_verified"
        or raw_receipt.get("row_id") != row["row_id"]
        or raw_receipt.get("md5") != row["fastq_md5"].lower()
        or raw_receipt.get("labels_accessed") is not False
        or raw_receipt.get("fit_or_score_performed") is not False
        or not fastq.is_file()
    ):
        raise GSE105127QuantificationError("verified raw RNA member differs")
    reference_receipt = json.loads((reference_root / "receipt.json").read_text(encoding="utf-8"))
    verify_frozen_tree(reference_root)
    if (
        reference_receipt.get("status") != "passed"
        or reference_receipt.get("assembly") != "GRCh38.p14"
        or reference_receipt.get("annotation") != "GENCODE_v49"
        or reference_receipt.get("labels_accessed") is not False
        or reference_receipt.get("fit_or_score_performed") is not False
    ):
        raise GSE105127QuantificationError("RSEM reference contract differs")
    if shutil.which("rsem-calculate-expression") is None or shutil.which("STAR") is None:
        raise GSE105127QuantificationError("RSEM or STAR is unavailable")
    stage = attempt / f"{row['row_id']}.staging"
    stage.mkdir(parents=True, exist_ok=False)
    prefix = stage / "rsem"
    command = [
        "rsem-calculate-expression", "--star", "--star-path",
        str(Path(shutil.which("STAR") or "").parent), "--star-gzipped-read-file",
        "--num-threads", str(threads), "--no-bam-output", "--strandedness", "none",
        "--seed", str(seed),
        str(fastq), str(reference_root / "reference/gse105127_gencode_v49"), str(prefix),
    ]
    with (stage / "command.stdout.txt").open("x", encoding="utf-8") as stdout, (stage / "command.stderr.txt").open("x", encoding="utf-8") as stderr:
        result = subprocess.run(command, stdout=stdout, stderr=stderr, text=True, check=False)
    if result.returncode != 0:
        raise GSE105127QuantificationError("RSEM quantification failed")
    genes = stage / "rsem.genes.results"
    metrics = audit_gene_results(genes)
    compressed = stage / "rsem.genes.results.gz"
    with genes.open("rb") as source, compressed.open("xb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as target:
            shutil.copyfileobj(source, target, length=8 * 1024 * 1024)
    genes.unlink()
    for path in stage.glob("rsem.*"):
        if path != compressed:
            if path.is_dir():
                shutil.rmtree(path)
            else:
                path.unlink()
    receipt = {
        "schema_version": "masld-bench-gse105127-rna-rsem-v1",
        "status": "passed_raw_scale",
        "row_id": row["row_id"],
        "participant_group_id": row["participant_group_id"],
        "zone": row["zone"],
        "pairing_topology": "adjacent_section",
        "input_layout": "single_end_76bp",
        "strandedness": "none",
        "reference": "GRCh38.p14_GENCODE_v49",
        "measurements": ["continuous_expected_count", "TPM", "effective_length", "length"],
        "rounding_applied": False,
        "rsem_native_tpm_retained": True,
        "post_quantification_normalization_or_feature_filter_applied": False,
        "labels_accessed": False,
        "fit_or_score_performed": False,
        "seed": seed,
        "raw_artifacts_sha256": sha256_file(raw / "ARTIFACTS.json"),
        "reference_artifacts_sha256": sha256_file(reference_root / "ARTIFACTS.json"),
        "gene_results_sha256": sha256_file(compressed),
        **metrics,
    }
    (stage / "receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    freeze_tree(stage, {"artifact_class": "gse105127_rna_rsem", "row_id": row["row_id"], "status": "passed"})
    verify_frozen_tree(stage)
    final.parent.mkdir(parents=True, exist_ok=True)
    publish_directory_noreplace(stage, final)
    verify_frozen_tree(final)
    return receipt


def run_bundle(
    *, plan_root: Path, raw_root: Path, reference_root: Path, bundle_id: int,
    output: Path, attempt: Path, threads: int, seed: int,
) -> dict[str, object]:
    verify_frozen_tree(plan_root)
    rows = read_rows(plan_root / "rna_rows.tsv", bundle_id)
    output.mkdir(parents=True, exist_ok=True)
    attempt.mkdir(parents=True, exist_ok=True)
    receipts = [quantify_one(row, raw_root, reference_root, output, attempt, threads, seed) for row in rows]
    genes = {int(receipt["genes"]) for receipt in receipts}
    if len(genes) != 1:
        raise GSE105127QuantificationError("RSEM gene census differs across rows")
    result = {
        "schema_version": "masld-bench-gse105127-rna-quant-bundle-v1",
        "status": "complete",
        "bundle_id": bundle_id,
        "participants": len({row["participant_group_id"] for row in rows}),
        "participant_zone_rows": len(rows),
        "genes": next(iter(genes)),
        "labels_accessed": False,
        "fit_or_score_performed": False,
    }
    bundle = output / f"bundle_{bundle_id:02d}"
    if (bundle / "COMPLETE").is_file():
        verify_frozen_tree(bundle)
        if json.loads((bundle / "receipt.json").read_text(encoding="utf-8")) != result:
            raise GSE105127QuantificationError("replayed RNA quantification bundle differs")
    elif bundle.exists():
        raise GSE105127QuantificationError("incomplete final RNA quantification bundle exists")
    else:
        bundle_stage = attempt / f"bundle_{bundle_id:02d}.staging"
        bundle_stage.mkdir(exist_ok=False)
        (bundle_stage / "receipt.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        freeze_tree(bundle_stage, {"artifact_class": "gse105127_rna_quant_bundle", "bundle_id": bundle_id, "status": "passed"})
        publish_directory_noreplace(bundle_stage, bundle)
        verify_frozen_tree(bundle)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-root", type=Path, required=True)
    parser.add_argument("--raw-root", type=Path, required=True)
    parser.add_argument("--reference-root", type=Path, required=True)
    parser.add_argument("--bundle-id", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--attempt", type=Path, required=True)
    parser.add_argument("--threads", type=int, required=True)
    parser.add_argument("--seed", type=int, default=105127)
    args = parser.parse_args()
    print(json.dumps(run_bundle(
        plan_root=args.plan_root, raw_root=args.raw_root,
        reference_root=args.reference_root, bundle_id=args.bundle_id,
        output=args.output, attempt=args.attempt, threads=args.threads, seed=args.seed,
    ), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
