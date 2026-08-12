#!/usr/bin/env python3
"""Real-data smoke test for pooled plus two-donor peak reproducibility."""

from __future__ import annotations

import argparse
import importlib.util
import sys
import tempfile
from pathlib import Path

import snapatac2 as snap


SCRIPT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPT_DIR))
MODULE_PATH = SCRIPT_DIR / "02b_call_consensus_peaks.py"
SPEC = importlib.util.spec_from_file_location("atac_v3_peak_recovery", MODULE_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"cannot load {MODULE_PATH}")
PEAKS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PEAKS)


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cohort", default="GSE244832")
    parser.add_argument("--lineage", default="macrophage")
    parser.add_argument("--n-jobs", type=int, default=4)
    return parser.parse_args()


def main() -> None:
    args = arguments()
    with tempfile.TemporaryDirectory(prefix="atac_v3_peak_support_") as directory:
        temp = Path(directory)
        dataset, adatas = PEAKS.build_dataset(args.cohort, temp / "dataset.h5ads")
        try:
            groups = PEAKS.eligible_donor_groups(dataset)[args.lineage]
            pooled = snap.tl.macs3(
                dataset,
                groupby="lineage_v3",
                selections={args.lineage},
                qvalue=PEAKS.POOLED_QVALUE,
                blacklist=PEAKS.BLACKLIST,
                shift=-100,
                extsize=200,
                tempdir=temp,
                inplace=False,
                n_jobs=args.n_jobs,
            )
            donor = snap.tl.macs3(
                dataset,
                groupby="donor_lineage_v3",
                selections=set(groups),
                qvalue=PEAKS.DONOR_QVALUE,
                blacklist=PEAKS.BLACKLIST,
                shift=-100,
                extsize=200,
                tempdir=temp,
                inplace=False,
                n_jobs=args.n_jobs,
            )
            frame = PEAKS.normalize_peak_frame(pooled[args.lineage])
            support = [0] * frame.height
            for group in groups:
                flags = PEAKS.overlap_flags(frame, PEAKS.normalize_peak_frame(donor[group]))
                support = [value + int(flag) for value, flag in zip(support, flags)]
            retained = sum(value >= PEAKS.MIN_SUPPORTING_DONORS for value in support)
            if frame.height == 0 or retained == 0:
                raise RuntimeError(
                    f"invalid peak-support smoke result: pooled={frame.height}, retained={retained}"
                )
            print(
                f"Peak-support smoke passed: cohort={args.cohort} lineage={args.lineage} "
                f"eligible_donors={len(groups)} pooled={frame.height} retained={retained}"
            )
        finally:
            PEAKS.close_dataset(dataset, adatas)


if __name__ == "__main__":
    main()
