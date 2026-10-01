#!/usr/bin/env python3
"""Validate the archived P4 saturation payloads, not just their directory presence.

Every P4 number so far carries the qualifier `present_not_integrity_validated`:
the reduction confirmed that a region directory and its six scorer files exist,
which is not the same as confirming that each file opens, carries the quantiles
layer the features are built from, and holds the track count its own request.json
recorded. A silently truncated or partially written HDF5 would have been counted
as present.

Per region this opens all six archived scorer files and checks:
  1. request.json parses and names the interval the directory name encodes;
  2. every scorer in the request is present as a readable .h5ad;
  3. the quantiles layer exists, because the features refuse a file without it;
  4. X and quantiles share a shape, and the variant axis equals 3 x the number of
     substitutable reference bases (three alternate bases per ACGT position), the
     axis the position statistic reduces over; a reference N carries no
     substitution, so the count is read from the reference rather than assumed to
     be 3 x interval bp;
  5. the archived track count equals request.json's `tracks_archived`;
  6. a deterministic slice of the quantiles layer is finite and inside [-1, 1]
     (AVI_SCORE is [0, 1]), which catches truncation that still opens cleanly.

No feature, statistic or family outcome is computed or read here. This is a
payload check over an existing deposit; it retrieves nothing and writes nothing
into the archive.
"""
import argparse
import csv
import gzip
import json
from multiprocessing import get_context
from pathlib import Path
import time

SCORERS = ("ATAC", "DNASE", "CHIP_HISTONE", "CHIP_TF", "CAGE", "AVI_SCORE")
SIGNED = {"ATAC", "DNASE", "CHIP_HISTONE", "CHIP_TF", "CAGE"}
SLICE_STRIDE = 97  # coprime with 3, so the slice crosses all three alternate bases
FASTA = "/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa"
_GENOME = None


def _reference(chrom, start0, end):
    """One FASTA handle per worker process; 18,570 reopens would dominate the run."""
    global _GENOME
    if _GENOME is None:
        import pysam
        _GENOME = pysam.FastaFile(FASTA)
    return _GENOME.fetch(str(chrom), int(start0), int(end)).upper()


def check_region(task):
    archive, universe, region_key, chrom, start0, end = task
    d = Path(archive) / universe / region_key.replace(":", "_")
    row = {"region_key": region_key, "universe": universe, "state": "ok", "detail": "",
           "scorers_checked": 0, "n_variants": 0, "reference_gap_bases": 0}
    try:
        import anndata
        import numpy as np
        req_path = d / "request.json"
        if not req_path.exists():
            row.update(state="missing_request", detail="request.json absent")
            return row
        request = json.loads(req_path.read_text())
        if request.get("interval") != region_key:
            row.update(state="interval_mismatch", detail=f"request={request.get('interval')}")
            return row
        counts = request.get("track_counts") or {}
        # Three alternate bases per SUBSTITUTABLE position. A reference N carries no
        # substitution, so 3 x interval_bp is the wrong expectation wherever the
        # assembly has a gap base; chr10:58267943 is one, and the two 500-bp windows
        # covering it archive 1,497 rows, not 1,500. The expectation is read from the
        # reference, so a genuinely truncated payload is still caught.
        _seq = _reference(chrom, start0, end)
        n_substitutable = sum(1 for _b in _seq if _b in "ACGT")
        expected_variants = 3 * n_substitutable
        row["reference_gap_bases"] = (int(end) - int(start0)) - n_substitutable
        for scorer in request.get("requested_scorers", SCORERS):
            path = d / f"{scorer}.h5ad"
            if not path.exists():
                row.update(state="missing_scorer", detail=scorer)
                return row
            ad = anndata.read_h5ad(path)
            if "quantiles" not in ad.layers:
                row.update(state="no_quantiles_layer", detail=scorer)
                return row
            x = np.asarray(ad.X)
            q = np.asarray(ad.layers["quantiles"])
            if x.shape != q.shape:
                row.update(state="shape_mismatch", detail=f"{scorer} X{x.shape} q{q.shape}")
                return row
            if x.shape[0] != expected_variants:
                row.update(state="variant_axis_mismatch",
                           detail=f"{scorer} {x.shape[0]} != 3x{n_substitutable} substitutable")
                return row
            archived = (counts.get(scorer) or {}).get("tracks_archived")
            if archived is not None and int(archived) != int(x.shape[1]):
                row.update(state="track_count_mismatch",
                           detail=f"{scorer} {x.shape[1]} != recorded {archived}")
                return row
            sample = q[::SLICE_STRIDE]
            if not np.isfinite(sample).all():
                row.update(state="nonfinite_quantiles", detail=scorer)
                return row
            low = 0.0 if scorer not in SIGNED else -1.0
            if sample.size and (float(sample.min()) < low - 1e-6 or float(sample.max()) > 1.0 + 1e-6):
                row.update(state="quantiles_out_of_range",
                           detail=f"{scorer} [{float(sample.min()):.4f},{float(sample.max()):.4f}]")
                return row
            row["scorers_checked"] += 1
            row["n_variants"] = int(x.shape[0])
            del ad, x, q
    except Exception as exc:  # a payload that cannot be opened is the thing being detected
        row.update(state="exception", detail=f"{type(exc).__name__}: {exc}"[:300])
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True,
                        help="Saturation run root holding tables/saturation_universe.tsv.gz and raw/saturation")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--limit", type=int, default=0,
                        help="Check only the first N regions; for a smoke run before the full pass.")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)

    with gzip.open(args.run / "tables/saturation_universe.tsv.gz", "rt") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    seen, tasks = set(), []
    archive = str(args.run / "raw/saturation")
    for r in rows:
        key = (r["universe"], r["region_key"])
        if key in seen:
            continue
        seen.add(key)
        tasks.append((archive, r["universe"], r["region_key"], r["chrom"], r["start0"], r["end"]))

    if args.limit:
        tasks = tasks[:args.limit]
    started = time.time()
    failures, states = [], {}
    with gzip.open(args.out / "archive_integrity.tsv.gz", "wt", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["region_key", "universe", "state", "detail",
                                                "scorers_checked", "n_variants", "reference_gap_bases"],
                                delimiter="\t", lineterminator="\n")
        writer.writeheader()
        with get_context("fork").Pool(args.workers) as pool:
            for i, row in enumerate(pool.imap_unordered(check_region, tasks, chunksize=8), 1):
                writer.writerow(row)
                states[row["state"]] = states.get(row["state"], 0) + 1
                if row["state"] != "ok":
                    failures.append(row)
                if i % 2000 == 0:
                    rate = i / max(time.time() - started, 1e-9)
                    print(f"checked {i}/{len(tasks)} bad={len(failures)} "
                          f"{rate*3600:.0f}/h eta={(len(tasks)-i)/max(rate,1e-9)/3600:.2f}h", flush=True)

    with open(args.out / "integrity_failures.tsv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["region_key", "universe", "state", "detail",
                                                "scorers_checked", "n_variants", "reference_gap_bases"],
                                delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(failures)

    summary = {"status": "archive_integrity_validated" if not failures else "archive_integrity_failures",
               "run": str(args.run), "regions": len(tasks), "ok": states.get("ok", 0),
               "failures": len(failures), "states": states, "workers": args.workers,
               "wall_seconds": time.time() - started,
               "checks": ["request_parses_and_matches_interval", "all_requested_scorers_readable",
                          "quantiles_layer_present", "X_and_quantiles_same_shape",
                          "variant_axis_equals_3x_substitutable_bases", "track_count_matches_request",
                          "sampled_quantiles_finite_and_in_range"],
               "slice_stride": SLICE_STRIDE,
               "limit": "a strided sample of the quantiles layer is range-checked, not every element",
               "no_saturation_feature_or_family_outcome_read": True}
    (args.out / "integrity_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)
    if failures:
        raise SystemExit(4)


if __name__ == "__main__":
    main()
