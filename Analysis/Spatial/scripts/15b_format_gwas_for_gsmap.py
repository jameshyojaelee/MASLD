#!/usr/bin/env python3
"""
15b_format_gwas_for_gsmap.py — Format GWAS summary statistics for gsMap.

Converts 7 harmonized GWAS files into gsMap-compatible .sumstats.gz format
(columns: SNP, A1, A2, Z, N) and generates a YAML config listing all traits.

gsMap resource bundle uses hg19 with rsID matching (build-agnostic for SNP
matching), so hg38 harmonized files work as long as rsIDs are present.

SLURM: --partition=cpu --cpus=4 --mem=64G --time=2:00:00
"""

import pathlib
import sys
import gzip

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from spatial_utils import PROJECT_ROOT, load_config, print_header, print_step

# ── GWAS registry ───────────────────────────────────────────────────────────
# Each entry: (file_path_relative, sample_size, gwas_type, column_mapping)
# column_mapping: dict mapping source columns -> gsMap columns
# gsMap requires: SNP (rsID), A1 (effect), A2 (other), Z, N

GWAS_REGISTRY = {
    "ghodsian_nafld": {
        "file": "GWAS/MR_Data/Ghodsian_2021_NAFLD_harmonised.tsv.gz",
        "N": 778614,
        "type": "binary",
        "columns": {
            "snp": "hm_rsid",
            "a1": "hm_effect_allele",
            "a2": "hm_other_allele",
            "beta": "hm_beta",
            "se": "standard_error",
            "pval": "p_value",
        },
    },
    "ukbb_alt": {
        "file": "GWAS/MR_Data/GCST90019492_UKBB_ALT_harmonised.tsv.gz",
        "N": 343850,
        "type": "quant",
        "columns": {
            "snp": "rsid",
            "a1": "effect_allele",
            "a2": "other_allele",
            "beta": "beta",
            "se": "standard_error",
            "pval": "p_value",
        },
    },
    "ukbb_ast": {
        "file": "GWAS/MR_Data/GCST90019497_UKBB_AST_harmonised.tsv.gz",
        "N": 343850,
        "type": "quant",
        "columns": {
            "snp": "rsid",
            "a1": "effect_allele",
            "a2": "other_allele",
            "beta": "beta",
            "se": "standard_error",
            "pval": "p_value",
        },
    },
    "ukbb_ggt": {
        "file": "GWAS/MR_Data/GCST90019507_UKBB_GGT_harmonised.tsv.gz",
        "N": 343850,
        "type": "quant",
        "columns": {
            "snp": "rsid",
            "a1": "effect_allele",
            "a2": "other_allele",
            "beta": "beta",
            "se": "standard_error",
            "pval": "p_value",
        },
    },
    "pdff": {
        "file": "GWAS/MR_Data/GCST90267352_PDFF_Pazoki2022.tsv.gz",
        "N": None,  # per-SNP N in column 'n'
        "type": "quant",
        "columns": {
            "snp": "variant_id",  # rsIDs already
            "a1": "effect_allele",
            "a2": "other_allele",
            "z": "z_value",
            "n": "n",
        },
    },
    "finngen_nafld": {
        "file": "GWAS/MR_Data/FinnGen/finngen_R12_NAFLD.gz",
        "N": 438857,  # 4,614 cases + 434,243 controls
        "type": "binary",
        "columns": {
            "snp": "rsids",
            "a1": "alt",  # FinnGen: alt = effect allele
            "a2": "ref",  # FinnGen: ref = other allele
            "beta": "beta",
            "se": "sebeta",
            "pval": "pval",
        },
    },
    "finngen_nash": {
        "file": "GWAS/MR_Data/FinnGen/finngen_R12_CHIRHEP_NAS.gz",
        "N": 436066,  # 1,823 cases + 434,243 controls
        "type": "binary",
        "columns": {
            "snp": "rsids",
            "a1": "alt",
            "a2": "ref",
            "beta": "beta",
            "se": "sebeta",
            "pval": "pval",
        },
    },
}


def format_one_gwas(trait_name, trait_info, output_dir):
    """Format a single GWAS file to gsMap sumstats format."""
    gwas_path = PROJECT_ROOT / trait_info["file"]
    if not gwas_path.exists():
        print(f"  WARNING: {gwas_path} not found, skipping {trait_name}")
        return None

    print(f"\n  Loading {trait_name}: {gwas_path.name}")

    # Detect separator: FinnGen uses tab, others use tab
    # FinnGen files start with #chrom header
    df = pd.read_csv(gwas_path, sep="\t", comment=None, low_memory=False)
    # Handle FinnGen header with # prefix
    if df.columns[0].startswith("#"):
        df.columns = [c.lstrip("#") for c in df.columns]

    print(f"    Raw: {len(df):,} SNPs, columns: {list(df.columns)[:8]}...")

    cols = trait_info["columns"]
    out = pd.DataFrame()

    # SNP (rsID)
    out["SNP"] = df[cols["snp"]]

    # Alleles
    out["A1"] = df[cols["a1"]].str.upper()
    out["A2"] = df[cols["a2"]].str.upper()

    # Z-score: use directly if available, else compute from beta/SE or beta/P
    if "z" in cols:
        out["Z"] = df[cols["z"]].astype(float)
    elif "beta" in cols and "se" in cols:
        beta = df[cols["beta"]].astype(float)
        se = df[cols["se"]].astype(float)
        out["Z"] = beta / se
    elif "beta" in cols and "pval" in cols:
        from scipy.stats import chi2
        beta = df[cols["beta"]].astype(float)
        pval = df[cols["pval"]].astype(float)
        out["Z"] = np.sqrt(chi2.isf(pval, 1)) * np.where(beta < 0, -1, 1)
    else:
        print(f"    ERROR: Cannot compute Z for {trait_name}")
        return None

    # Sample size
    if "n" in cols:
        out["N"] = df[cols["n"]].astype(float).astype(int)
    elif trait_info["N"] is not None:
        out["N"] = trait_info["N"]
    else:
        print(f"    ERROR: No sample size for {trait_name}")
        return None

    # ── QC filters ──────────────────────────────────────────────────────
    n_before = len(out)

    # Drop missing
    out = out.dropna(subset=["SNP", "A1", "A2", "Z", "N"])

    # Remove non-rsID SNPs
    out = out[out["SNP"].str.startswith("rs", na=False)]

    # Remove strand-ambiguous SNPs (A/T, C/G)
    valid_pairs = {"AC", "AG", "CA", "CT", "GA", "GT", "TC", "TG"}
    allele_pair = out["A1"] + out["A2"]
    out = out[allele_pair.isin(valid_pairs)]

    # Remove duplicated rsIDs
    out = out.drop_duplicates(subset="SNP", keep="first")

    # Remove infinite or extreme Z
    out = out[np.isfinite(out["Z"])]

    n_after = len(out)
    print(f"    QC: {n_before:,} -> {n_after:,} SNPs "
          f"(removed {n_before - n_after:,})")

    # Summary stats
    chisq = out["Z"] ** 2
    print(f"    Mean chi2 = {chisq.mean():.3f}, "
          f"Lambda GC = {chisq.median() / 0.4549:.3f}, "
          f"N GW-sig = {(chisq > 29).sum():,}")

    # Save
    out_path = output_dir / f"{trait_name}.sumstats.gz"
    out.to_csv(out_path, sep="\t", index=False, float_format="%.4f",
               compression="gzip")
    size_mb = out_path.stat().st_size / 1e6
    print(f"    Saved: {out_path} ({size_mb:.1f} MB)")

    return str(out_path)


def write_gwas_config(trait_paths, output_dir):
    """Write gwas_config.yaml for gsMap sumstats_config_file."""
    import yaml

    config_path = output_dir / "gwas_config.yaml"
    with open(config_path, "w") as f:
        yaml.dump(trait_paths, f, default_flow_style=False, sort_keys=False)
    print(f"\n  GWAS config saved: {config_path}")
    print(f"  Traits: {list(trait_paths.keys())}")
    return config_path


def main():
    print_header("15b: Format GWAS for gsMap")

    config = load_config()
    output_dir = PROJECT_ROOT / "Analysis/Spatial/data/gsmap_gwas"
    output_dir.mkdir(parents=True, exist_ok=True)

    trait_paths = {}
    for i, (trait_name, trait_info) in enumerate(GWAS_REGISTRY.items(), 1):
        print_step(f"Formatting {trait_name}", i, len(GWAS_REGISTRY))
        out_path = format_one_gwas(trait_name, trait_info, output_dir)
        if out_path is not None:
            trait_paths[trait_name] = out_path

    # Write YAML config
    write_gwas_config(trait_paths, output_dir)

    print(f"\n  Summary: {len(trait_paths)}/{len(GWAS_REGISTRY)} GWAS formatted")

    print_header("15b: Complete")


if __name__ == "__main__":
    main()
