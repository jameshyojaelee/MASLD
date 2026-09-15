#!/usr/bin/env python3
"""D1 task 1: per-donor, per-lineage barcode lists for the two snATAC cohorts.

Label source. The assignment names `Analysis/ATAC/Human_External/snapatac2/per_donor/*.h5ad` for
GSE281367. That directory holds only Z01..Z04 and all four read as 0 x 0 under snapatac2 2.8.0; the
producing log stops during Z04, so the run was interrupted and no cell-type label is stored there.
The joint co-embedding export carries the labels for all 30 donors of both cohorts, and its
GSE244832 half reproduces `cell_donor_condition_map.tsv.gz` (the label source the assignment names
for that cohort) exactly. The deviation is recorded in d1_prespec.json, written before any count.

Barcode reconstruction. The joint index is `<barcode>` with one `-<k>` suffix per concat pass.
GSE281367 cells carry a 16-base cellranger barcode, so the CB tag is the first 16 bases plus `-1`.
GSE244832 cells carry a 22-base combinatorial barcode, which is the read-name prefix in that BAM.

Writes:
  barcodes/<donor>__<lineage>.txt   one barcode per line
  tables/lineage_cell_counts.tsv    donor x lineage cell counts
"""
from __future__ import annotations

import pathlib
import re
import sys

import pandas as pd

PROJECT = pathlib.Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
JOINT = PROJECT / "Analysis/ATAC/Human_External/results/coembed/umap_export_joint.tsv.gz"
MULTIOME_MAP = PROJECT / "Analysis/ATAC/Human_Multiome/results/snapatac2/cell_donor_condition_map.tsv.gz"

PRIMARY = ["Hepatocyte", "Endothelial", "Stellate_Cell", "Macrophage", "Kupffer_Cell",
           "Cholangiocyte", "LSEC", "Plasma_Cell", "NK_T_Cell", "B_Cell"]
POOLS = {"Myeloid": ["Macrophage", "Kupffer_Cell"]}
LEAD = re.compile(r"^([ACGT]+)")


def cb_from_index(raw: str, cohort: str) -> str:
    """Strip anndata make-unique suffixes and return the barcode as it appears in the BAM."""
    m = LEAD.match(raw)
    if m is None:
        raise ValueError(f"index {raw!r} does not start with a barcode")
    lead = m.group(1)
    if cohort == "GSE281367":
        if len(lead) < 16:
            raise ValueError(f"GSE281367 index {raw!r} shorter than 16 bases")
        return lead[:16] + "-1"
    if len(lead) < 22:
        raise ValueError(f"GSE244832 index {raw!r} shorter than 22 bases")
    return lead[:22]


def load_labels() -> pd.DataFrame:
    d = pd.read_csv(JOINT, sep="\t").rename(columns={"Unnamed: 0": "joint_index"})
    d["barcode"] = [cb_from_index(r, c) for r, c in zip(d["joint_index"], d["cohort"])]
    d["donor"] = d["donor_id"].str.split("_", n=1).str[1]
    dup = d.duplicated(["donor", "barcode"]).sum()
    if dup:
        raise SystemExit(f"{dup} duplicated (donor, barcode) pairs after reconstruction")
    # cross-check against the cohort-native label file the assignment names for GSE244832
    m = pd.read_csv(MULTIOME_MAP, sep="\t")
    native = m.groupby(["donor_id", "cell_type"]).size().rename("n_native")
    joint = (d[d["cohort"] == "GSE244832"].groupby(["donor", "cell_type"]).size()
             .rename("n_joint").rename_axis(["donor_id", "cell_type"]))
    chk = pd.concat([native, joint], axis=1).fillna(0)
    mismatch = int((chk["n_native"] != chk["n_joint"]).sum())
    print(f"[d1_01] GSE244832 cross-check: {len(chk)} donor x cell_type cells, {mismatch} mismatched")
    if mismatch:
        raise SystemExit("joint export does not reproduce the GSE244832 native label file")
    return d


def main() -> None:
    out = pathlib.Path(sys.argv[1]).resolve()
    bcdir = out / "barcodes"
    bcdir.mkdir(parents=True, exist_ok=True)
    (out / "tables").mkdir(parents=True, exist_ok=True)
    d = load_labels()

    rows = []
    for (donor, cohort), g in d.groupby(["donor", "cohort"]):
        streams = {lin: g.loc[g["cell_type"] == lin, "barcode"] for lin in PRIMARY}
        for pool, members in POOLS.items():
            streams[pool] = g.loc[g["cell_type"].isin(members), "barcode"]
        streams["ALL_LABELLED"] = g["barcode"]
        for lin, s in streams.items():
            p = bcdir / f"{donor}__{lin}.txt"
            p.write_text("\n".join(sorted(s)) + ("\n" if len(s) else ""))
            rows.append({"cohort": cohort, "donor": donor, "lineage": lin, "n_cells": int(len(s))})
    t = pd.DataFrame(rows).sort_values(["cohort", "donor", "lineage"])
    t.to_csv(out / "tables/lineage_cell_counts.tsv", sep="\t", index=False)
    piv = t.pivot_table(index="lineage", columns="cohort", values="n_cells", aggfunc="sum")
    print(piv.to_string())
    print(f"[d1_01] wrote {len(rows)} barcode lists for {t['donor'].nunique()} donors -> {bcdir}")


if __name__ == "__main__":
    main()
