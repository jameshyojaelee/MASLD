#!/usr/bin/env python
"""
66_pops_magma_prep.py  --  PoPS Axis-2 (which-gene) MAGMA prep
==============================================================
ADDITIVE (new src/66). Does NOT edit any existing 46d/78/27a/60-65 script.

Builds the two inputs a PoPS-compatible MAGMA gene analysis needs on the
PRIMARY MASLD GWAS (MVP_NAFLD_EUR, EUR, N=432,709, 29,342 cases -- the
best-powered NAFLD case-control in the finemapping registry):

  1. An ENSGID-keyed MAGMA gene-location file, derived from the PoPS gene
     annotation (tools/pops/example/data/utils/gene_annot_jun10.txt, GRCh37).
     REQUIRED because PoPS feature rows are keyed by ENSGID, whereas the repo's
     existing magma_annot.genes.annot is Entrez-keyed (NCBI37.3) -> incompatible.

  2. A MAGMA --pval input (SNP=rsID, P=pval, N) for MVP_NAFLD_EUR, mapping
     chr:pos (hg19) -> rsID via the 1000G-EUR lookup already built by the repo
     (tools/magma/work/chr_pos_to_rsid.tsv), exactly as scripts/figures/
     magma_gwas_prep.py does for the other GWAS.

Outputs:
  tools/magma/work/pops/gene_annot_jun10.ensgid.gene.loc
  tools/magma/work/pops/MVP_NAFLD_EUR.magma_input.tsv
"""
from __future__ import annotations
import os, sys, time
from pathlib import Path
import pandas as pd

BASE = Path(os.environ.get("MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"))

# --- primary MASLD GWAS ------------------------------------------------------
GWAS_ID   = "MVP_NAFLD_EUR"
GWAS_N    = 432709          # N_tot from GWAS/finemapping/config/gwas_registry.tsv
SUMSTATS  = BASE / "GWAS/finemapping/data/sumstats/MVP_NAFLD_EUR_reformatted_hg19.tsv"

LOOKUP    = BASE / "tools/magma/work/chr_pos_to_rsid.tsv"
GENE_ANNOT= BASE / "tools/pops/example/data/utils/gene_annot_jun10.txt"
OUT_DIR   = BASE / "tools/magma/work/pops"
OUT_DIR.mkdir(parents=True, exist_ok=True)
GENE_LOC  = OUT_DIR / "gene_annot_jun10.ensgid.gene.loc"
MAGMA_IN  = OUT_DIR / f"{GWAS_ID}.magma_input.tsv"


def log(m): print(f"[66_pops_magma_prep] {time.strftime('%H:%M:%S')}  {m}", flush=True)


def build_gene_loc():
    """ENSGID-keyed MAGMA .gene.loc : GENE CHR START STOP  (0kb window at annotate)."""
    ga = pd.read_csv(GENE_ANNOT, sep="\t", dtype=str)
    # columns: ENSGID NAME CHR START END TSS
    ga = ga.dropna(subset=["ENSGID", "CHR", "START", "END"])
    loc = ga[["ENSGID", "CHR", "START", "END"]].copy()
    loc.to_csv(GENE_LOC, sep="\t", header=False, index=False)
    log(f"wrote ENSGID gene.loc: {GENE_LOC}  ({len(loc):,} genes)")


def build_magma_input():
    log(f"loading chr:pos->rsID lookup {LOOKUP}")
    lut = pd.read_csv(LOOKUP, sep="\t", dtype=str)
    lut = dict(zip(lut["chr_pos"], lut["rsid"]))
    log(f"  lookup entries: {len(lut):,}")

    log(f"reading sumstats {SUMSTATS.name}")
    ss = pd.read_csv(SUMSTATS, sep="\t",
                     usecols=["chromosome", "position", "pval"],
                     dtype={"chromosome": str, "position": str, "pval": float})
    n_in = len(ss)
    ss["chr_pos"] = ss["chromosome"].str.removeprefix("chr") + ":" + ss["position"]
    ss["SNP"] = ss["chr_pos"].map(lut)
    n_match = int(ss["SNP"].notna().sum())
    log(f"  rows in: {n_in:,}   joined to rsID: {n_match:,} ({100*n_match/n_in:.1f}%)")

    out = ss.dropna(subset=["SNP", "pval"]).copy()
    out = out[(out["pval"] > 0) & (out["pval"] <= 1)]
    out["N"] = GWAS_N
    out[["SNP", "pval", "N"]].rename(columns={"pval": "P"}).to_csv(
        MAGMA_IN, sep="\t", index=False)
    log(f"wrote MAGMA input: {MAGMA_IN}  ({len(out):,} SNPs, N={GWAS_N})")


if __name__ == "__main__":
    for p in (SUMSTATS, LOOKUP, GENE_ANNOT):
        if not p.exists():
            log(f"FATAL missing input: {p}"); sys.exit(1)
    build_gene_loc()
    build_magma_input()
    log("done.")
