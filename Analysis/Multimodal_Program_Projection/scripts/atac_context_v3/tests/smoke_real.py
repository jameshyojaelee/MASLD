#!/usr/bin/env python3
"""Bounded real-data smoke test for SnapATAC2 dataset, MACS3, and counting APIs."""

from __future__ import annotations

import argparse
import tempfile
from pathlib import Path

import snapatac2 as snap

import sys

SCRIPT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPT_DIR))

from snap_helpers import build_dataset, close_dataset  # noqa: E402


OLD_PEAKS = (
    Path(__file__).resolve().parents[5]
    / "Analysis/ATAC/Human_Multiome/results/label_transfer/cell_type_peak_sets_v2/Hepatocytes_peaks.bed"
)
BLACKLIST = Path("/gpfs/commons/home/jameslee/reference_genome/blacklists/hg38-blacklist.v2.bed")


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-jobs", type=int, default=4)
    return parser.parse_args()


def main() -> None:
    args = arguments()
    with tempfile.TemporaryDirectory(prefix="atac_v3_smoke_") as directory:
        temp = Path(directory)
        with OLD_PEAKS.open("r", encoding="utf-8") as source, (temp / "peaks50.bed").open(
            "x", encoding="utf-8"
        ) as destination:
            for _, line in zip(range(50), source):
                destination.write(line)
        for cohort in ("GSE244832", "GSE281367"):
            dataset, adatas = build_dataset(cohort, temp / f"{cohort}.h5ads")
            try:
                if dataset.n_obs <= 0 or len(dataset.obs["lineage_v3"]) != dataset.n_obs:
                    raise RuntimeError(f"invalid assembled dataset: {cohort}")
                matrix = snap.pp.make_peak_matrix(
                    dataset,
                    peak_file=temp / "peaks50.bed",
                    counting_strategy="fragment",
                    file=temp / f"{cohort}.matrix.h5ad",
                    inplace=False,
                )
                try:
                    if matrix.n_obs != dataset.n_obs or matrix.n_vars != 50:
                        raise RuntimeError(f"real-data peak matrix shape mismatch: {cohort}")
                    if (
                        len(matrix.obs["lineage_v3"]) != matrix.n_obs
                        or len(matrix.obs["donor_id"]) != matrix.n_obs
                    ):
                        raise RuntimeError(f"peak matrix lost annotations: {cohort}")
                finally:
                    matrix.close()
                if cohort == "GSE244832":
                    peaks = snap.tl.macs3(
                        dataset,
                        groupby="lineage_v3",
                        replicate="donor_id",
                        selections={"hepatocyte"},
                        qvalue=0.05,
                        replicate_qvalue=0.05,
                        blacklist=BLACKLIST,
                        tempdir=temp,
                        inplace=False,
                        n_jobs=args.n_jobs,
                    )
                    if peaks is None or "hepatocyte" not in peaks or len(peaks["hepatocyte"]) == 0:
                        raise RuntimeError("bounded donor-replicated MACS3 smoke test returned no peaks")
            finally:
                close_dataset(dataset, adatas)
    print("Real-data SnapATAC2/MACS3/exact-count smoke test passed")


if __name__ == "__main__":
    main()
