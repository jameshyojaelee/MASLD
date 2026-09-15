#!/usr/bin/env python3
"""Step 41 (P4): fix the saturation region universe BEFORE any saturation score is read.

Universe U1 (primary): the 96,460 H3K27ac regions of the released chromatin predictor (keys as shipped; the
coordinate convention is inferred, immaterial at region width).
Universe U2 (snATAC consensus peaks, 500 bp, hg38 0-based half-open), restricted to peaks that carry a measured
or genetic reason to be scored — the full 676,658-peak manifest is out of budget:
  U2a  hepatocyte or stellate peaks with a disease-change evidence state other than `indeterminate`
       (supported / discordant / source_dependent) in `da/da_peak_results.tsv.gz`;
  U2b  hepatocyte/stellate peaks significant in either cohort (BH q < 0.05) — includes U2a;
  U2c  peaks of any lineage overlapping a queried Track 0 variant (variant_peak_overlaps.tsv.gz, when present);
  U2d  peaks whose promoter_genes field names a gene of one of the 117 read-only programs.
Each peak keeps its membership flags so families stay separable downstream.

Outputs (tables/): saturation_universe.tsv.gz, saturation_universe_summary.json
"""

from __future__ import annotations

import json
from collections import defaultdict

import numpy as np
import pandas as pd

import lib_atlas as la

ROOT = la.out_root()
TABLES = ROOT / "tables"
BENCH = la.PROJECT / "Analysis/MASLD_Model_Benchmark"
RELEASE = BENCH / "release/masld-liver-chromatin-state-v1.2/weights/chromatin_state_v1_2.npz"
CTX = la.PROJECT / "Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-candidate-2026-08-11-r1"
PROGRAMS = la.PROJECT / "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot/program_membership_v2.tsv"
TRACK0 = la.track0_root() / "tables/variant_peak_overlaps.tsv.gz"


def main() -> None:
    rel = np.load(RELEASE, allow_pickle=True)
    rows = []
    for k in rel["prof_region_key"]:
        chrom, rng_ = str(k).split(":"); a, b = rng_.split("-")
        rows.append({"universe": "U1_h3k27ac", "region_key": str(k), "chrom": chrom, "start0": int(a), "end": int(b), "lineage": "bulk_liver_h3k27ac",
                     "u2a_da_state": "", "u2b_sig_either": False, "u2c_genetic_overlap": False, "u2d_program_linked": False})
    peaks = pd.read_csv(CTX / "consensus_peak_manifest.tsv", sep="\t", dtype=str)
    peaks = peaks[peaks["blacklist_overlap"] != "TRUE"]
    peaks["pid"] = peaks["lineage"] + "|" + peaks["chrom"] + ":" + peaks["start0"] + "-" + peaks["end"]
    da = pd.read_csv(CTX / "da/da_peak_results.tsv.gz", sep="\t", keep_default_na=False, dtype=str)
    da["pid"] = da["lineage"] + "|" + da["peak_coordinate"]
    q1 = pd.to_numeric(da["qvalue_gse244832"], errors="coerce"); q2 = pd.to_numeric(da["qvalue_gse281367"], errors="coerce")
    da["sig_either"] = (q1 < 0.05) | (q2 < 0.05)
    da_state = dict(zip(da["pid"], da["evidence_state"])); sig = dict(zip(da["pid"], da["sig_either"]))
    prog_genes = set(pd.read_csv(PROGRAMS, sep="\t", dtype=str)["canonical_gene"].dropna())
    prom = dict(zip(da["pid"], da["promoter_genes"]))
    genetic = set()
    if TRACK0.exists():
        vp = pd.read_csv(TRACK0, sep="\t", dtype=str)
        genetic = set(vp["lineage"] + "|" + vp["peak"])
    else:
        la.log("Track 0 variant_peak_overlaps not present yet: U2c left empty (rerun after the tables job to add it)")
    n = defaultdict(int)
    for r in peaks.itertuples(index=False):
        pid = r.pid
        st = da_state.get(pid, "")
        u2a = st not in ("", "indeterminate") and r.lineage in ("hepatocyte", "stellate")
        u2b = bool(sig.get(pid, False)) and r.lineage in ("hepatocyte", "stellate")
        u2c = pid in genetic
        u2d = any(g in prog_genes for g in str(prom.get(pid, "")).replace(";", ",").split(",") if g)
        if not (u2a or u2b or u2c or u2d):
            continue
        rows.append({"universe": "U2_snatac", "region_key": f"{r.chrom}:{r.start0}-{r.end}", "chrom": r.chrom, "start0": int(r.start0), "end": int(r.end), "lineage": r.lineage,
                     "u2a_da_state": st if u2a else "", "u2b_sig_either": u2b, "u2c_genetic_overlap": u2c, "u2d_program_linked": u2d})
        for k, v in (("u2a", u2a), ("u2b", u2b), ("u2c", u2c), ("u2d", u2d)):
            n[k] += int(v)
    df = pd.DataFrame(rows).drop_duplicates(["universe", "region_key", "lineage"])
    df["width"] = df["end"] - df["start0"]
    df.to_csv(TABLES / "saturation_universe.tsv.gz", sep="\t", index=False)
    summary = {"U1_h3k27ac": int((df.universe == "U1_h3k27ac").sum()), "U2_snatac_peaks": int((df.universe == "U2_snatac").sum()), "U2_flags": dict(n),
               "U2_unique_coordinates": int(df[df.universe == "U2_snatac"]["region_key"].nunique()), "total_bp": int(df.drop_duplicates(["region_key"])["width"].sum()),
               "fixed_before": "any saturation score was read (2026-09-09)", "rule": __doc__.strip().split("\n")[2:14]}
    json.dump(summary, (TABLES / "saturation_universe_summary.json").open("w"), indent=1)
    la.log(f"saturation universe: {summary}")


if __name__ == "__main__":
    main()
