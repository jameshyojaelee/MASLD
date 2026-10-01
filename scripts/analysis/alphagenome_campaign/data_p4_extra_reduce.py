#!/usr/bin/env python3
"""Reduce a P4 saturation deposit other than the original universe.

`data_p4_full.py` is pinned to the original run and its 103,390 registered
regions. The F5 null background (U3_f5_background) and the genetic-overlap family
(U2c_genetic_overlap) are separate deposits with their own universes, so the same
reduction is pointed at them here. The feature definitions, track grouping,
life-stage separation and column order are the originals, imported unchanged, so
a background feature and a U1 feature remain the same quantity.

Like the original, this refuses to reduce a partial deposit: every region in the
universe must be present before any feature is written, because a family tested
on whichever regions happened to arrive first is tested on an arrival order, not
on a population.

No retrieval, raw copy or substitution rescore occurs. No family outcome is read.
"""
import argparse
import csv
import gzip
import json
from multiprocessing import get_context
from pathlib import Path
import time

import data_p4 as p4
from data_p4_full import feature_columns, process


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True,
                        help="Saturation run root: tables/saturation_universe.tsv.gz and raw/saturation")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--expect", type=int, required=True,
                        help="Region count this deposit must hold before any feature is written")
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)

    # census() and the region paths read module-level constants; point them at this deposit.
    p4.SAT = args.run.resolve()
    p4.OLD.group_columns = p4.group_columns

    groups = p4.track_registry(args.out)
    rows, counts, n_source = p4.census(args.out)
    missing = [r for r in rows if r["archive_state"] != "present_not_integrity_validated"]
    if len(rows) != args.expect or missing:
        reason = {"status": "dependency_incomplete", "run": str(args.run),
                  "expected_regions": args.expect, "observed_regions": len(rows),
                  "missing_or_incomplete": len(missing), "archive_presence": counts,
                  "full_reduction_run": False}
        (args.out / "INCOMPLETE.json").write_text(json.dumps(reason, indent=2) + "\n")
        print(json.dumps(reason), flush=True)
        raise SystemExit(3)

    columns = feature_columns()
    tfcols = ["region_key", "group", "rank", "track_name", "transcription_factor",
              "mean_rel_effect", "mean_max_abs_quantile"]
    tasks = [(str(p4.SAT / "raw/saturation"), r["region_key"], r["universe"], r["chrom"],
              int(r["start0"]), int(r["end"])) for r in rows]
    started = time.time()
    failures = []
    with gzip.open(args.out / "saturation_region_features.tsv.gz", "wt", newline="") as fh, \
         gzip.open(args.out / "saturation_region_tf_top5.tsv.gz", "wt", newline="") as th:
        fw = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", lineterminator="\n")
        tw = csv.DictWriter(th, fieldnames=tfcols, delimiter="\t", lineterminator="\n")
        fw.writeheader(); tw.writeheader()
        with get_context("fork").Pool(args.workers) as pool:
            for i, (row, top, state) in enumerate(pool.imap(process, tasks, chunksize=8), 1):
                row["archive_state"] = state
                u1 = row["universe"] == "U1_h3k27ac"
                row["counted_interval_start0_inferred"] = int(row["start0"]) - int(u1)
                row["initial_counted_interval_bp_not_scored"] = int(u1)
                fw.writerow(row); tw.writerows(top)
                if state != "ok":
                    failures.append({"region_key": row["region_key"], "universe": row["universe"],
                                     "state": state})
                if i % 1000 == 0:
                    print(f"reduced {i}/{len(tasks)} failures={len(failures)} "
                          f"elapsed={time.time()-started:.1f}s", flush=True)
    p4.write_rows(args.out / "reduction_failures.tsv", failures, ["region_key", "universe", "state"])
    summary = {"status": "reduction_failed" if failures else "full_reduction_complete",
               "run": str(args.run), "regions": len(tasks), "failures": len(failures),
               "workers": args.workers, "wall_seconds": time.time() - started,
               "source_universe_rows": n_source, "archive_presence": counts,
               "feature_definitions": "imported unchanged from the original P4 reduction",
               "group_counts": {f"{a}:{b}": n for (a, b), n in groups.items()},
               "family_outcomes_read": False}
    (args.out / "reduction_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)
    if failures:
        raise SystemExit(4)


if __name__ == "__main__":
    main()
