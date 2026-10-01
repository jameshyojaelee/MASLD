#!/usr/bin/env python3
"""Complete P4 reduction only after every registered archive member is present.

No retrieval, raw copy, substitution rescore, or training feature export occurs.
Any reduction failure keeps scientific family tests closed. U1 archive coordinates
retain their historical one-base omission relative to inferred counted intervals.
"""
import argparse
import csv
import gzip
import json
from multiprocessing import get_context
from pathlib import Path
import time

import data_p4 as p4


def process(task):
    return p4.OLD.process_region(task)


def feature_columns():
    names=set()
    for scorer in p4.OLD.SCORERS:
        with (p4.META / f"{scorer}.track_metadata.tsv").open() as handle:
            rows=list(csv.DictReader(handle,delimiter="\t"))
        var={key:[row[key] for row in rows] for key in rows[0]} if rows else {}
        names.update(p4.group_columns(scorer,var))
    names=sorted(names)
    identity=["region_key","universe","chrom","start0","end","n_positions","archive_state","counted_interval_start0_inferred","initial_counted_interval_bp_not_scored"]
    return identity+[f"{g}_n_tracks" for g in names]+[f"{g}_{sc}_{s}" for g in names for sc in p4.OLD.SCALES for s in p4.OLD.STATS]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",type=Path,required=True)
    parser.add_argument("--workers",type=int,default=4)
    args=parser.parse_args()
    args.out.mkdir(parents=True,exist_ok=False)
    groups=p4.track_registry(args.out)
    rows,counts,n_source=p4.census(args.out)
    missing=[r for r in rows if r["archive_state"]!="present_not_integrity_validated"]
    if len(rows)!=103390 or missing:
        reason={"status":"dependency_incomplete","expected_regions":103390,"observed_regions":len(rows),"missing_or_incomplete":len(missing),"archive_presence":counts,"full_reduction_run":False}
        (args.out / "INCOMPLETE.json").write_text(json.dumps(reason,indent=2)+"\n")
        print(json.dumps(reason),flush=True)
        raise SystemExit(3)
    p4.OLD.group_columns=p4.group_columns
    columns=feature_columns()
    tfcols=["region_key","group","rank","track_name","transcription_factor","mean_rel_effect","mean_max_abs_quantile"]
    tasks=[(str(p4.SAT / "raw/saturation"),r["region_key"],r["universe"],r["chrom"],int(r["start0"]),int(r["end"])) for r in rows]
    started=time.time()
    failures=[]
    with gzip.open(args.out / "saturation_region_features.tsv.gz","wt",newline="") as fh, gzip.open(args.out / "saturation_region_tf_top5.tsv.gz","wt",newline="") as th:
        fw=csv.DictWriter(fh,fieldnames=columns,delimiter="\t",lineterminator="\n")
        tw=csv.DictWriter(th,fieldnames=tfcols,delimiter="\t",lineterminator="\n")
        fw.writeheader();tw.writeheader()
        with get_context("fork").Pool(args.workers) as pool:
            for i,(row,top,state) in enumerate(pool.imap(process,tasks,chunksize=8),1):
                row["archive_state"]=state
                u1=row["universe"]=="U1_h3k27ac"
                row["counted_interval_start0_inferred"]=int(row["start0"])-int(u1)
                row["initial_counted_interval_bp_not_scored"]=int(u1)
                fw.writerow(row);tw.writerows(top)
                if state!="ok":
                    failures.append({"region_key":row["region_key"],"universe":row["universe"],"state":state})
                if i%1000==0:
                    print(f"reduced {i}/103390 failures={len(failures)} elapsed={time.time()-started:.1f}s",flush=True)
    p4.write_rows(args.out / "reduction_failures.tsv",failures,["region_key","universe","state"])
    summary={"status":"reduction_failed" if failures else "full_reduction_complete","regions":len(tasks),"failures":len(failures),"workers":args.workers,"wall_seconds":time.time()-started,"source_universe_rows":n_source,"archive_presence":counts,"U1_coordinate_limit":"archived[a,b);inferred_counted_interval[a-1,b);initial_base_not_scored","hosted_features_training_allowed":False,"group_counts":{f"{a}:{b}":n for (a,b),n in groups.items()},"scientific_F1_F5_complete":False}
    (args.out / "reduction_summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    if failures:
        raise SystemExit(4)


if __name__=="__main__":
    main()
