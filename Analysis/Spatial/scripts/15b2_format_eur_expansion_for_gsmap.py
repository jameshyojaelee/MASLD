#!/usr/bin/env python3
"""
15b2_format_eur_expansion_for_gsmap.py — Format 11 additional EUR GWAS for gsMap.

Ancestry-safe expansion of the gsMap trait panel (7 -> 18). All added traits are
EUR (matched to the EUR-only gsMap LD resource bundle); BBJ/EAS is deliberately
NOT added here (needs a genome-wide EAS LD-score/baseline build — deferred).

Source = finemapping preprocessed/reformatted sumstats
(`GWAS/finemapping/data/sumstats/*_{preprocessed,reformatted_hg19}.tsv`), which are
hg19 and carry chr:pos:allele1:allele2:beta:se:pval but NO rsID. gsMap keys on rsID
(HapMap3), so we annotate chr:pos(hg19) -> rsID using the gsMap bundle's OWN
1000G EUR Phase3 bim (`data/gsmap_resource/LD_Reference_Panel/1000G_EUR_Phase3_plink`).
This guarantees every emitted rsID is in gsMap's SNP universe.

gsMap spatial LDSC is chi^2-based (chisq = Z**2), so allele orientation / Z sign is
immaterial to the enrichment; we only need correct rsID, a valid Z magnitude, and N.

Writes {trait}.sumstats.gz (SNP,A1,A2,Z,N) and REWRITES gwas_config.yaml to include
all 18 traits (original backed up). Re-running 15e then computes spatial_ldsc +
cauchy for the 11 new traits only (per-sample steps resume from .done).

Run on a COMPUTE node (reads ~0.3-1.2 GB TSVs, builds a ~10M-SNP bim map).
SLURM: --partition=cpu --cpus=4 --mem=96G --time=48:00:00
"""
import gzip
import pathlib
import sys

import numpy as np
import pandas as pd

PROJECT_ROOT = pathlib.Path(
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
SUMSTATS_DIR = PROJECT_ROOT / "GWAS/finemapping/data/sumstats"
BIM_DIR = (PROJECT_ROOT
           / "data/gsmap_resource/LD_Reference_Panel/1000G_EUR_Phase3_plink")
OUT_DIR = PROJECT_ROOT / "Analysis/Spatial/data/gsmap_gwas"
CONFIG = OUT_DIR / "gwas_config.yaml"

# ── New EUR trait registry: gsMap_key -> (source filename, N_total) ───────────
# N from GWAS/finemapping/config/gwas_registry.tsv (N_tot). Binary traits use
# total N (matches the existing 7-trait convention in 15b: FinnGen/Ghodsian).
NEW_TRAITS = {
    "mvp_nafld":          ("MVP_NAFLD_EUR_reformatted_hg19.tsv",                    432709),
    "mvp_alt":            ("MVP_ALT_EUR_reformatted_hg19.tsv",                      406803),
    "mvp_ast":            ("MVP_AST_EUR_reformatted_hg19.tsv",                      389192),
    "decode_nafld":       ("2023_36280732_NAFLD_deCode_EUR_reformatted_hg19.tsv",  359187),
    "ukbb2023_nafld":     ("2023_36280732_NAFLD_UKBB_EUR_reformatted_hg19.tsv",    397020),
    "intermtn_nafld":     ("2023_36280732_NAFLD_Intermountain_EUR_reformatted_hg19.tsv", 32689),
    "anstee2020_nafld":   ("2020_32298765_NAFLD_EUR_preprocessed.tsv",              9491),
    "nafld_2019":         ("2019_31311600_NAFLD_EUR_preprocessed.tsv",              8434),
    "pdff_2021a":         ("2021_34128465_PDFF_EUR_preprocessed.tsv",              36116),
    "pdff_2021b":         ("2021_34957434_PDFF_EUR_preprocessed.tsv",              32858),
    "pdff_2022":          ("2022_36402844_PDFF_EUR_preprocessed.tsv",              44867),
}

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
AMBIG = {"AT", "TA", "CG", "GC"}


def build_bim_map():
    """chr:pos(hg19) -> rsID via 1000G EUR Phase3 bim. Returns DataFrame keyed
    on (chrom, pos) with columns rsid, b1, b2 (uppercase alleles)."""
    frames = []
    for chrom in range(1, 23):
        bim = BIM_DIR / f"1000G.EUR.QC.{chrom}.bim"
        if not bim.exists():
            print(f"  WARNING: missing {bim.name}")
            continue
        df = pd.read_csv(bim, sep="\t", header=None,
                         usecols=[0, 1, 3, 4, 5],
                         names=["chrom", "rsid", "pos", "b1", "b2"],
                         dtype={"chrom": "int32", "rsid": str, "pos": "int64",
                                "b1": str, "b2": str})
        df = df[df["rsid"].str.startswith("rs", na=False)]
        df["b1"] = df["b1"].str.upper()
        df["b2"] = df["b2"].str.upper()
        frames.append(df)
    bim_map = pd.concat(frames, ignore_index=True)
    # Drop positions with >1 rsID (ambiguous join target)
    dup = bim_map.duplicated(subset=["chrom", "pos"], keep=False)
    n_dup = int(dup.sum())
    bim_map = bim_map[~dup]
    print(f"  bim map: {len(bim_map):,} unique-position SNPs "
          f"(dropped {n_dup:,} multi-rsID positions)")
    return bim_map.set_index(["chrom", "pos"])


def allele_set_match(a1, a2, b1, b2):
    """True if {a1,a2} matches {b1,b2} on the same or complementary strand."""
    g = {a1, a2}
    if g == {b1, b2}:
        return True
    gc = {COMP.get(a1, "N"), COMP.get(a2, "N")}
    return gc == {b1, b2}


def format_one(trait, filename, N, bim_map):
    src = SUMSTATS_DIR / filename
    if not src.exists():
        print(f"  ERROR: {src} not found — skipping {trait}")
        return None
    print(f"\n  [{trait}] {filename}  (N={N:,})")

    df = pd.read_csv(src, sep="\t", low_memory=False,
                     usecols=["chromosome", "position", "allele1", "allele2",
                              "beta", "se"])
    n_raw = len(df)
    df = df.rename(columns={"chromosome": "chrom", "position": "pos",
                            "allele1": "A1", "allele2": "A2"})
    df["chrom"] = pd.to_numeric(df["chrom"], errors="coerce")
    df = df.dropna(subset=["chrom", "pos"])
    df["chrom"] = df["chrom"].astype("int32")
    df["pos"] = df["pos"].astype("int64")
    df["A1"] = df["A1"].str.upper()
    df["A2"] = df["A2"].str.upper()

    # Join rsID by (chrom, pos)
    df = df.join(bim_map, on=["chrom", "pos"], how="inner")
    n_pos = len(df)

    # Allele-set validation (drops indels / mismatches / wrong multi-allelic)
    ok = [allele_set_match(a1, a2, b1, b2) for a1, a2, b1, b2
          in zip(df["A1"], df["A2"], df["b1"], df["b2"])]
    df = df[pd.Series(ok, index=df.index)]
    n_allele = len(df)

    # Z from beta/se
    df["beta"] = pd.to_numeric(df["beta"], errors="coerce")
    df["se"] = pd.to_numeric(df["se"], errors="coerce")
    df["Z"] = df["beta"] / df["se"]

    out = df[["rsid", "A1", "A2", "Z"]].rename(columns={"rsid": "SNP"})
    out["N"] = N

    # QC filters (mirror 15b)
    out = out.dropna(subset=["SNP", "A1", "A2", "Z", "N"])
    out = out[np.isfinite(out["Z"])]
    pair = out["A1"] + out["A2"]
    out = out[~pair.isin(AMBIG)]          # strand-ambiguous
    out = out.drop_duplicates(subset="SNP", keep="first")
    n_final = len(out)

    chisq = out["Z"] ** 2
    print(f"    raw={n_raw:,} -> pos-join={n_pos:,} -> allele-ok={n_allele:,} "
          f"-> final={n_final:,}  ({100*n_final/max(n_raw,1):.1f}% kept)")
    print(f"    mean chi2={chisq.mean():.3f}  lambdaGC={chisq.median()/0.4549:.3f}  "
          f"N GW-sig={(chisq>29).sum():,}")
    if chisq.mean() < 1.0:
        print("    WARNING: mean chi2 < 1.0 — check formatting")

    out_path = OUT_DIR / f"{trait}.sumstats.gz"
    out.to_csv(out_path, sep="\t", index=False, float_format="%.4f",
               compression="gzip")
    print(f"    saved {out_path.name} ({out_path.stat().st_size/1e6:.1f} MB)")
    return str(out_path)


def rewrite_config(new_paths):
    import yaml
    # Load existing 7-trait config, back it up, merge
    existing = {}
    if CONFIG.exists():
        with open(CONFIG) as f:
            existing = yaml.safe_load(f) or {}
        bak = CONFIG.with_suffix(".yaml.bak_7trait")
        if not bak.exists():
            with open(bak, "w") as f:
                yaml.dump(existing, f, default_flow_style=False, sort_keys=False)
            print(f"\n  backed up original config -> {bak.name}")
    merged = dict(existing)
    merged.update(new_paths)
    with open(CONFIG, "w") as f:
        yaml.dump(merged, f, default_flow_style=False, sort_keys=False)
    print(f"  wrote {CONFIG.name}: {len(merged)} traits "
          f"({len(existing)} existing + {len(new_paths)} new)")
    print(f"  traits: {list(merged.keys())}")


def main():
    print("=" * 70)
    print("  15b2: Format 11 EUR GWAS for gsMap (ancestry-safe expansion)")
    print("=" * 70)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("\n  Building chr:pos(hg19) -> rsID map from 1000G EUR bim ...")
    bim_map = build_bim_map()

    new_paths = {}
    for trait, (fn, N) in NEW_TRAITS.items():
        p = format_one(trait, fn, N, bim_map)
        if p:
            new_paths[trait] = p

    print(f"\n  Formatted {len(new_paths)}/{len(NEW_TRAITS)} traits")
    rewrite_config(new_paths)
    print("\n  15b2: Complete")


if __name__ == "__main__":
    main()
