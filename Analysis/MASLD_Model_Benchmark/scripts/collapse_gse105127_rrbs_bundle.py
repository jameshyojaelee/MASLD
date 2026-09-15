#!/usr/bin/env python3
"""Collapse BisSNP strand rows for one participant-contained RRBS bundle."""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path
import subprocess

from masld_bench.artifacts import freeze_tree, publish_directory_noreplace, verify_frozen_tree


AWK = r'''
BEGIN { OFS="\t"; print "contig_rank","source_chrom","source_start0","source_end0","methylated_reads","coverage","strand_state","auxiliary_state" }
function emit( state, aux) {
    if (!have) return
    if (n == 1) {
        state = plus_seen ? "plus_only" : "minus_only"
        aux = only_aux == 0 ? "unpaired_aux_zero" : "unpaired_aux_nonzero"
        single++
    } else if (n == 2 && plus_seen && minus_seen) {
        state = "both_strands"
        if (plus_aux == minus_cov && minus_aux == plus_cov) { aux="reciprocal_partner_coverage"; reciprocal++ }
        else { aux="auxiliary_inconsistent_ignored"; aux_inconsistent++ }
        paired++
    } else { hard_error++ ; return }
    print current_rank,current_chrom,current_start,current_start+2,meth_sum,cov_sum,state,aux
    cpgs++; total_meth+=meth_sum; total_cov+=cov_sum
}
function reset_group(rank, chrom, start) {
    current_rank=rank; current_chrom=chrom; current_start=start
    n=0; plus_seen=0; minus_seen=0; meth_sum=0; cov_sum=0
    plus_aux=0; minus_aux=0; plus_cov=0; minus_cov=0; only_aux=0; have=1
}
/^track / { track++; next }
{
    source_rows++
    if (NF != 11) { hard_error++; next }
    chrom=$1
    if (chrom == "chrX") rank=23
    else if (chrom == "chrY") rank=24
    else if (chrom == "chrM") rank=25
    else if (chrom ~ /^chr[0-9]+$/) { rank=substr(chrom,4)+0; if (rank < 1 || rank > 22) rank=0 }
    else rank=0
    start=$2+0; finish=$3+0; strand=$6; cov=$5+0; pct=$4+0; aux=$11+0
    canonical = strand == "+" ? start : strand == "-" ? start-1 : -1
    implied=pct*cov/100.0; meth=int(implied+0.5); residual=implied-meth; if(residual<0) residual=-residual
    tolerance=cov*0.00005+0.00001
    if (rank == 0 || start < 0 || finish != start+1 || canonical < 0 || cov < 1 || pct < 0 || pct > 100 || residual > tolerance) { hard_error++; next }
    if (have && (rank < current_rank || (rank == current_rank && canonical < current_start))) { hard_error++; next }
    if (!have || rank != current_rank || canonical != current_start) { emit(); reset_group(rank,chrom,canonical) }
    n++; meth_sum+=meth; cov_sum+=cov; only_aux=aux
    if (strand == "+") { if (plus_seen) hard_error++; plus_seen=1; plus_cov=cov; plus_aux=aux }
    else if (strand == "-") { if (minus_seen) hard_error++; minus_seen=1; minus_cov=cov; minus_aux=aux }
    else hard_error++
    detail_nonzero += ($10+0 > 0)
}
END {
    emit()
    print "source_rows\t" source_rows > metrics
    print "track_rows\t" track >> metrics
    print "canonical_cpgs\t" cpgs >> metrics
    print "paired_cpgs\t" paired >> metrics
    print "single_strand_cpgs\t" single >> metrics
    print "reciprocal_aux_pairs\t" reciprocal >> metrics
    print "auxiliary_inconsistent_pairs\t" aux_inconsistent >> metrics
    print "detail_nonzero_rows\t" detail_nonzero >> metrics
    print "methylated_reads\t" total_meth >> metrics
    print "coverage\t" total_cov >> metrics
    print "hard_errors\t" hard_error >> metrics
    if (track != 1 || hard_error != 0) exit 42
}
'''


class GSE105127CollapseError(RuntimeError):
    """Raised when BisSNP strand collapse differs."""


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
            raise GSE105127CollapseError("RRBS plan violates the label firewall")
        rows = [row for row in reader if int(row["bundle_id"]) == bundle_id]
    if not rows:
        raise GSE105127CollapseError("RRBS bundle is empty")
    groups = {row["participant_group_id"] for row in rows}
    if any(sum(item["participant_group_id"] == group for item in rows) != 3 for group in groups):
        raise GSE105127CollapseError("participant zones are split inside RRBS bundle")
    return rows


def parse_metrics(path: Path) -> dict[str, int]:
    observed: dict[str, int] = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            key, value = line.rstrip("\n").split("\t", 1)
            observed[key] = int(value or 0)
    return observed


def collapse_one(row: dict[str, str], output: Path, attempt: Path) -> dict[str, object]:
    final = output / "participants" / row["row_id"]
    if (final / "COMPLETE").is_file():
        verify_frozen_tree(final)
        receipt = json.loads((final / "receipt.json").read_text(encoding="utf-8"))
        expected = {
            "row_id": row["row_id"],
            "participant_group_id": row["participant_group_id"],
            "zone": row["zone"],
            "source_bed_sha256": row["source_bed_sha256"],
        }
        if any(receipt.get(key) != value for key, value in expected.items()):
            raise GSE105127CollapseError("replayed RRBS member differs from plan")
        return receipt
    if final.exists():
        raise GSE105127CollapseError("incomplete final RRBS member exists")
    source = Path(row["source_bed_path"])
    if (
        not source.is_file()
        or source.stat().st_size != int(row["source_bed_bytes"])
        or sha256_file(source) != row["source_bed_sha256"]
    ):
        raise GSE105127CollapseError("frozen RRBS source identity differs")
    stage = attempt / f"{row['row_id']}.staging"
    stage.mkdir(parents=True, exist_ok=False)
    uncompressed = stage / "cpg_counts.hg19.tsv"
    compressed = stage / "cpg_counts.hg19.tsv.gz"
    metrics_path = stage / "metrics.tsv"
    decompressor = subprocess.Popen(["gzip", "-dc", str(source)], stdout=subprocess.PIPE)
    if decompressor.stdout is None:
        raise GSE105127CollapseError("RRBS gzip pipe is unavailable")
    with uncompressed.open("xb") as handle:
        audit = subprocess.run(
            ["awk", "-v", f"metrics={metrics_path}", AWK],
            stdin=decompressor.stdout, stdout=handle, stderr=subprocess.PIPE,
            text=False, check=False,
        )
    decompressor.stdout.close()
    gzip_status = decompressor.wait()
    if gzip_status != 0 or audit.returncode != 0:
        raise GSE105127CollapseError(
            f"RRBS collapse failed: gzip={gzip_status} awk={audit.returncode} {audit.stderr.decode(errors='replace')}"
        )
    metrics = parse_metrics(metrics_path)
    if (
        metrics.get("source_rows") != int(row["source_cytosine_rows"])
        or metrics.get("track_rows") != 1
        or metrics.get("hard_errors") != 0
        or metrics.get("canonical_cpgs", 0) <= 0
        or metrics.get("paired_cpgs", 0) + metrics.get("single_strand_cpgs", 0)
        != metrics.get("canonical_cpgs")
    ):
        raise GSE105127CollapseError("RRBS collapse census differs")
    with uncompressed.open("rb") as source_handle, compressed.open("xb") as target_handle:
        gzip_run = subprocess.run(["gzip", "-n", "-c"], stdin=source_handle, stdout=target_handle, check=False)
    if gzip_run.returncode != 0 or subprocess.run(["gzip", "-t", str(compressed)], check=False).returncode != 0:
        raise GSE105127CollapseError("collapsed RRBS gzip failed")
    uncompressed.unlink()
    receipt = {
        "schema_version": "masld-bench-gse105127-collapsed-rrbs-v1",
        "status": "passed",
        "row_id": row["row_id"],
        "participant_group_id": row["participant_group_id"],
        "zone": row["zone"],
        "pairing_topology": "adjacent_section",
        "source_bed_sha256": row["source_bed_sha256"],
        "collapsed_sha256": sha256_file(compressed),
        "coordinate_build": "1000_Genomes_GRCh37",
        "coordinate_convention": "zero_based_half_open_two_base_CpG",
        "measurement": "integer_methylated_reads_and_coverage",
        "auxiliary_columns_used_for_measurement": False,
        "missing_as_zero": False,
        "labels_accessed": False,
        "fit_or_score_performed": False,
        **metrics,
    }
    (stage / "receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    freeze_tree(stage, {"artifact_class": "gse105127_collapsed_rrbs", "row_id": row["row_id"], "status": "passed"})
    verify_frozen_tree(stage)
    final.parent.mkdir(parents=True, exist_ok=True)
    publish_directory_noreplace(stage, final)
    verify_frozen_tree(final)
    return receipt


def run_bundle(*, plan_root: Path, bundle_id: int, output: Path, attempt: Path) -> dict[str, object]:
    verify_frozen_tree(plan_root)
    rows = read_rows(plan_root / "rrbs_rows.tsv", bundle_id)
    output.mkdir(parents=True, exist_ok=True)
    attempt.mkdir(parents=True, exist_ok=True)
    receipts = [collapse_one(row, output, attempt) for row in rows]
    result = {
        "schema_version": "masld-bench-gse105127-rrbs-collapse-bundle-v1",
        "status": "complete",
        "bundle_id": bundle_id,
        "participants": len({row["participant_group_id"] for row in rows}),
        "participant_zone_rows": len(rows),
        "source_rows": sum(int(receipt["source_rows"]) for receipt in receipts),
        "canonical_cpgs_across_samples": sum(int(receipt["canonical_cpgs"]) for receipt in receipts),
        "auxiliary_inconsistent_pairs": sum(int(receipt["auxiliary_inconsistent_pairs"]) for receipt in receipts),
        "labels_accessed": False,
        "fit_or_score_performed": False,
    }
    bundle = output / f"bundle_{bundle_id:02d}"
    if (bundle / "COMPLETE").is_file():
        verify_frozen_tree(bundle)
        if json.loads((bundle / "receipt.json").read_text(encoding="utf-8")) != result:
            raise GSE105127CollapseError("replayed RRBS bundle receipt differs")
    elif bundle.exists():
        raise GSE105127CollapseError("incomplete final RRBS bundle receipt exists")
    else:
        bundle_stage = attempt / f"bundle_{bundle_id:02d}.staging"
        bundle_stage.mkdir(exist_ok=False)
        (bundle_stage / "receipt.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        freeze_tree(bundle_stage, {"artifact_class": "gse105127_rrbs_collapse_bundle", "bundle_id": bundle_id, "status": "passed"})
        publish_directory_noreplace(bundle_stage, bundle)
        verify_frozen_tree(bundle)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-root", type=Path, required=True)
    parser.add_argument("--bundle-id", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--attempt", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(run_bundle(plan_root=args.plan_root, bundle_id=args.bundle_id, output=args.output, attempt=args.attempt), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
