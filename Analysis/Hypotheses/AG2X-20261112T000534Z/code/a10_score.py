#!/usr/bin/env python3
"""Stage A scoring (PRESPEC_AG2X.md sections 2, 4, 5, 6).

1. Build the deduplicated variant list for every Stage A endpoint (A1, A2, A2-hep, A3a, A3b/A4 tops at
   PIP >= 0.9 and A4 controls), assert REF against the GRCh38 FASTA, write inputs/stageA_variants.tsv.
2. AlphaGenome 2 (ALL_FOLDS), conditions bulk liver UBERON:0002107 and liver hepatocyte
   UBERON:0002107_CL:0000182, scorers ATAC, DNASE, SC_ATAC, RNA_SEQ, SC_GENE_EXPRESSION, 1-Mb interval.
   Raw tidy tables written whole to raw/v2_<condition>/.
3. AlphaGenome version 1 hosted, RNA_SEQ only, on the A3b/A4 variants; only rows of tracks with ontology
   UBERON:0002107, UBERON:0001114, UBERON:0001115 or CL:0000182 and of the A3b/A4 target genes are kept
   (raw/v1_rna_liver/). A 3-variant unfiltered pilot (raw/v1_pilot/) checks the metadata columns first.
Hosted-call ceiling 30,000 for Stage 0 + A, counted from all chunk logs plus the 220 Stage 0 calls.
"""
from __future__ import annotations

import json
import sys

import pandas as pd

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import ag2x_lib as L  # noqa: E402

PROG = L.REPO / "GWAS/finemapping/results/alphagenome_program"
C1 = PROG / "c1-endpoint2-v2-20260915T084559Z/tables"
RAW = L.PKG / "raw"
IN = L.PKG / "inputs"
CAP = 30_000
STAGE0_CALLS = 220
SCORERS = ["ATAC", "DNASE", "SC_ATAC", "RNA_SEQ", "SC_GENE_EXPRESSION"]
CONDITIONS = {"bulk_liver": "UBERON:0002107", "hepatocyte": "UBERON:0002107_CL:0000182"}
V1_ONT = {"UBERON:0002107", "UBERON:0001114", "UBERON:0001115", "CL:0000182"}


def calls_on_disk() -> int:
    n = 0
    for f in RAW.glob("*/chunk_log.jsonl"):
        n += sum(json.loads(x)["n_variants"] for x in open(f))
    return n


def variant_table():
    rows = []
    a1 = pd.read_csv(C1 / "endpoint1_matched_tierA4.tsv.gz", sep="\t", usecols=["chr_x", "pos_hg38", "ref", "alt"])
    rows += [(f"{c}:{p}:{r}:{a}", "A1") for c, p, r, a in zip(a1.chr_x, a1.pos_hg38, a1.ref, a1.alt)]
    rows += [(u, "A2") for u in pd.read_csv(C1 / "endpoint2_matched.tsv", sep="\t", usecols=["uid"]).uid]
    rows += [(u, "A2hep") for u in pd.read_csv(L.PKG / "tables/a2_hepatocyte_lineage_v1.tsv", sep="\t", usecols=["uid"]).uid]
    rows += [(u, "A3a") for u in pd.read_csv(L.PKG / "tables/a3a_gate_set_v1.tsv", sep="\t", usecols=["variant_uid"]).variant_uid]
    tops = pd.read_csv(IN / "a3b_a4_tops.tsv", sep="\t")
    tops = tops[tops.pip >= 0.9]
    rows += [(u, "A34top") for u in tops.uid]
    pairs = pd.read_csv(IN / "a4_pairs.tsv", sep="\t")
    rows += [(u, "A4ctrl") for u in pairs.ctrl_uid]
    v = pd.DataFrame(rows, columns=["uid", "endpoint"]).drop_duplicates()
    wide = v.assign(x=1).pivot_table(index="uid", columns="endpoint", values="x", fill_value=0).reset_index()
    genes = set(tops.ensg) | set(pairs.ensg)
    return wide.sort_values("uid").reset_index(drop=True), genes


def main():
    counter = {}
    wide, genes = variant_table()
    uids = wide.uid.tolist()
    bad = L.assert_ref(uids)
    wide["ref_ok"] = ~wide.uid.isin([b[0] for b in bad])
    wide.to_csv(IN / "stageA_variants.tsv", sep="\t", index=False)
    uids = wide.loc[wide.ref_ok, "uid"].tolist()
    eq = wide.loc[wide.ref_ok & ((wide.get("A34top", 0) == 1) | (wide.get("A4ctrl", 0) == 1)), "uid"].tolist()
    used = STAGE0_CALLS + calls_on_disk()
    L.log(f"variants {len(wide)}; REF mismatches {len(bad)}; eQTL variants {len(eq)}; target genes {len(genes)}; "
          f"hosted calls already on disk incl. Stage 0: {used}")
    planned_new = 2 * len(uids) + len(eq) + 3
    L.log(f"planned total {STAGE0_CALLS + planned_new} (cap {CAP})")
    if STAGE0_CALLS + planned_new > CAP:
        raise SystemExit("planned calls exceed cap; stopping before any call")

    v2 = L.create_v2()
    for name, cond in CONDITIONS.items():
        L.score_chunked(v2, "v2", uids, cond, SCORERS, RAW / f"v2_{name}", chunk=50, max_workers=4,
                        counter=counter, cap=CAP, used_before=used)
        L.log(f"v2 {name} done; counter {counter}")

    v1 = L.create_v1()
    pil = L.score_chunked(v1, "v1", eq[:3], None, ["RNA_SEQ"], RAW / "v1_pilot", chunk=3, max_workers=1,
                          counter=counter, cap=CAP, used_before=used)
    cols = pd.read_parquet(pil[0]).columns.tolist()
    L.log(f"v1 tidy columns: {cols}")
    if "ontology_curie" not in cols or "gene_id" not in cols:
        raise SystemExit("v1 tidy table lacks ontology_curie or gene_id; v1 rescoring stopped")

    def keep(t):
        g = t.gene_id.astype(str).str.split(".").str[0]
        return t[t.ontology_curie.astype(str).isin(V1_ONT) & g.isin(genes)]

    L.score_chunked(v1, "v1", eq, None, ["RNA_SEQ"], RAW / "v1_rna_liver", chunk=25, max_workers=4,
                    counter=counter, keep=keep, cap=CAP, used_before=used)
    total = STAGE0_CALLS + calls_on_disk()
    summary = {"utc": L.utc(), "package_version": L.package_version(), "variants": len(wide),
               "ref_mismatches": len(bad), "eqtl_variants": len(eq), "this_job_counter": counter,
               "hosted_calls_completed_incl_stage0": total,
               "hosted_calls_incl_resends_this_job": used + counter.get("attempted_calls", 0) + counter.get("resent_calls", 0)}
    json.dump(summary, open(RAW / "scoring_summary.json", "w"), indent=1)
    L.log(json.dumps(summary))


if __name__ == "__main__":
    main()
