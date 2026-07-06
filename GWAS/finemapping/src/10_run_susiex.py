#!/usr/bin/env python3
"""
10_run_susiex.py
Run SuSiEX joint multi-ancestry fine-mapping for one locus.

Generalized to N-ancestry (default EUR+EAS, opt-in EUR+EAS+AFR+SAS via
ANCESTRIES env var). Per-ancestry GWAS is derived from the locus row's
trait_pair (ALT/AST/GGT):
  EUR  → UKBB_<trait>            (n=343,850)
  EAS  → BBJ_<trait>             (n=160,000)
  AFR  → PanUKBB_AFR_<trait>     (n=6,636)
  SAS  → PanUKBB_CSA_<trait>     (n=8,876)

If a locus does not have ≥ MIN_SNPS in a given ancestry's window after
sumstats×LD intersection, that arm is silently dropped and SuSiEX is
invoked on the remaining arms.

Usage:
  ANCESTRIES=EUR,EAS,AFR,SAS python 10_run_susiex.py --locus-row N

LD reference layout (LD_PANEL env var picks the panel; default '1kg'):
  data/ld_ref/{1kg_eur,1kg_eas,1kg_afr,1kg_sas,polyfun_eur,topld_*}/
    approx_LD_blocks.txt     — tab-separated: chr start stop
    chr{N}/                   — block-level dirs (EUR / panel-specific)
    chr{N}_<anc>.{bed,bim,fam}  — chr-level PLINK (1kg_<anc> for non-EUR)
"""

import argparse
import os
import subprocess
import sys
import tempfile
import shutil

import pandas as pd

FM_DIR = os.environ.get(
    "FM_DIR",
    os.path.join(
        os.environ.get(
            "MASLD_PROJECT_ROOT",
            "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
        ),
        "GWAS/finemapping",
    ),
)

SUSIEX_BIN = os.path.join(FM_DIR, "bin/SuSiEx")
PLINK_BIN  = os.path.join(FM_DIR, "bin/plink")

# ---------------------------------------------------------------------------
# Per-ancestry registry (study name template + N + chr-level PLINK basename)
# ---------------------------------------------------------------------------
ANCESTRY_REGISTRY = {
    "EUR": {
        "gwas_template": "UKBB_{trait}",
        "n_gwas":        343850,
        "ref_kind":      "block",   # block-level .bed/.bim/.fam (per BP block)
    },
    "EAS": {
        "gwas_template": "BBJ_{trait}",
        "n_gwas":        160000,
        "ref_kind":      "chr",     # chr-level PLINK extract (chr{N}_eas.*)
        "chr_basename":  "chr{chrom}_eas",
    },
    "AFR": {
        "gwas_template": "PanUKBB_AFR_{trait}",
        "n_gwas":        6636,
        "ref_kind":      "chr",
        "chr_basename":  "chr{chrom}_afr",
    },
    "SAS": {
        "gwas_template": "PanUKBB_CSA_{trait}",
        "n_gwas":        8876,
        "ref_kind":      "chr",
        "chr_basename":  "chr{chrom}_sas",
    },
    # AMR: MVP-only arm (no enzyme UKBB/BBJ/PanUKBB template). gwas_template and
    # n_gwas are omitted deliberately — the MVP cohort path resolves study name
    # and N from the registry (see _build_sumstats_path / _lookup_n_gwas).
    "AMR": {
        "ref_kind":      "chr",
        "chr_basename":  "chr{chrom}_amr",
    },
}

# ---------------------------------------------------------------------------
# GWAS registry lookup (config/gwas_registry.tsv) — used for the MVP cohort,
# whose per-stratum study names and sample sizes are NOT a fixed per-ancestry
# constant (they vary by trait x ancestry).
# ---------------------------------------------------------------------------
_REGISTRY_CACHE = None


def _load_registry():
    global _REGISTRY_CACHE
    if _REGISTRY_CACHE is not None:
        return _REGISTRY_CACHE
    reg_path = os.path.join(FM_DIR, "config/gwas_registry.tsv")
    reg = {}
    with open(reg_path) as f:
        header = None
        for line in f:
            if line.lstrip().startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if header is None:
                header = parts
                continue
            row = dict(zip(header, parts))
            reg[row["study_name"]] = row
    _REGISTRY_CACHE = reg
    return reg


# Sumstats paths (relative to FM_DIR) for each (ancestry, trait). For the legacy
# enzyme cohort the study is derived from the per-ancestry gwas_template; for the
# MVP cohort it is MVP_<trait>_<ANC> resolved against the registry (so both the
# path and N follow the registry, not a hardcoded template).
def _build_sumstats_path(ancestry, trait, cohort="legacy"):
    if cohort == "MVP":
        study = f"MVP_{trait}_{ancestry}"
        reg = _load_registry()
        if study in reg and reg[study].get("sumstats_path"):
            return study, os.path.join(FM_DIR, reg[study]["sumstats_path"])
        # Not in registry (e.g. MVP_Cirrhosis_EAS does not exist) — signal missing.
        return study, None
    template = ANCESTRY_REGISTRY[ancestry].get("gwas_template")
    if template is None:
        return f"{ancestry}_{trait}", None
    study = template.format(trait=trait)
    rel = f"data/sumstats/{study}_reformatted_hg19.tsv"
    return study, os.path.join(FM_DIR, rel)


def _lookup_n_gwas(ancestry, trait, study, cohort="legacy"):
    """Per-arm GWAS sample size. MVP -> registry N_tot (varies per stratum);
    legacy enzyme -> fixed per-ancestry constant in ANCESTRY_REGISTRY."""
    if cohort == "MVP":
        reg = _load_registry()
        if study in reg and reg[study].get("N_tot"):
            return int(float(reg[study]["N_tot"]))
        raise KeyError(f"No N_tot in registry for MVP stratum {study}")
    return ANCESTRY_REGISTRY[ancestry]["n_gwas"]


def get_active_ancestries():
    """ANCESTRIES env var: comma-separated, default 'EUR,EAS' for backward compat."""
    raw = os.environ.get("ANCESTRIES", "EUR,EAS")
    ancs = [a.strip().upper() for a in raw.split(",") if a.strip()]
    for a in ancs:
        if a not in ANCESTRY_REGISTRY:
            raise ValueError(f"ANCESTRIES contains unsupported '{a}'; valid: {list(ANCESTRY_REGISTRY)}")
    return ancs


# ---------------------------------------------------------------------------
# LD panel dispatch (mirrors get_ld_base_dir() in finemapping_functions.R)
# ---------------------------------------------------------------------------
def _resolve_ld_dir(ancestry):
    override_env = {"EUR": "UKBB_LD_DIR", "EAS": "EAS_LD_DIR",
                    "AFR": "AFR_LD_DIR", "AMR": "AMR_LD_DIR",
                    "SAS": "SAS_LD_DIR"}[ancestry]
    if os.environ.get(override_env):
        return os.environ[override_env]
    panel = os.environ.get("LD_PANEL", "1kg").lower()
    if panel not in ("1kg", "topld", "ukbb", "polyfun"):
        raise ValueError(f"LD_PANEL='{panel}' not recognized")
    if panel == "polyfun" and ancestry == "EUR":
        return os.path.join(FM_DIR, "data/ld_ref", "polyfun_eur")
    if panel == "polyfun":  # non-EUR falls back to 1kg
        return os.path.join(FM_DIR, "data/ld_ref", f"1kg_{ancestry.lower()}")
    if panel == "ukbb" and ancestry == "EUR":
        return os.path.join(FM_DIR, "data/ld_ref", "ukbb_eur")
    prefix = panel if panel != "ukbb" else "1kg"
    return os.path.join(FM_DIR, "data/ld_ref", f"{prefix}_{ancestry.lower()}")


# ---------------------------------------------------------------------------
# LD block utilities
# ---------------------------------------------------------------------------

def load_ld_blocks(ld_dir):
    blocks_file = os.path.join(ld_dir, "approx_LD_blocks.txt")
    df = pd.read_csv(blocks_file, sep=r"\s+", engine="python")
    df.columns = [c.strip() for c in df.columns]
    df["chr"]   = df["chr"].astype(int)
    df["start"] = df["start"].astype(int)
    df["stop"]  = df["stop"].astype(int)
    return df


def find_ld_block_prefix(ld_dir, chrom, window_start, window_end):
    """List of block-level path prefixes overlapping the window."""
    blocks = load_ld_blocks(ld_dir)
    overlapping = blocks[
        (blocks["chr"] == chrom)
        & (blocks["start"] < window_end)
        & (blocks["stop"]  > window_start)
    ]
    if overlapping.empty:
        raise ValueError(f"No LD blocks for chr{chrom}:{window_start}-{window_end} in {ld_dir}")
    prefixes = []
    for _, row in overlapping.iterrows():
        block_name = f"{row['start']}.{row['stop']}"
        prefix = os.path.join(ld_dir, f"chr{chrom}", block_name, block_name)
        if not os.path.exists(prefix + ".bim"):
            raise FileNotFoundError(f"Block .bim missing: {prefix}.bim")
        if not os.path.exists(prefix + ".ld"):
            raise FileNotFoundError(f"Block .ld missing: {prefix}.ld")
        prefixes.append(prefix)
    return prefixes


# ---------------------------------------------------------------------------
# Sumstats utilities
# ---------------------------------------------------------------------------

def detect_format(df):
    cols = set(df.columns)
    if "pos_hg19" in cols:
        return {"CHR": "chromosome", "POS": "pos_hg19",
                "A1": "effect_allele", "A2": "other_allele",
                "BETA": "beta", "SE": "standard_error", "PVAL": "p_value"}
    elif "position" in cols and "allele1" in cols:
        return {"CHR": "chromosome", "POS": "position",
                "A1": "allele1", "A2": "allele2",
                "BETA": "beta", "SE": "se", "PVAL": "pval"}
    raise ValueError(f"Unrecognised sumstats format. Columns: {list(df.columns)}")


def load_sumstats_window(study_name, sumstats_path, chrom, window_start, window_end):
    print(f"  Loading {study_name} from {os.path.basename(sumstats_path)} ...")
    df = pd.read_csv(sumstats_path, sep="\t", low_memory=False)
    mapping = detect_format(df)
    df = df.rename(columns={v: k for k, v in mapping.items()})
    df["CHR"] = pd.to_numeric(df["CHR"], errors="coerce")
    df["POS"] = pd.to_numeric(df["POS"], errors="coerce")
    df["BETA"] = pd.to_numeric(df["BETA"], errors="coerce")
    df["SE"] = pd.to_numeric(df["SE"], errors="coerce")
    df = df.dropna(subset=["CHR", "POS", "BETA", "SE", "PVAL"])
    df = df[df["SE"] > 0]
    df["CHR"] = df["CHR"].astype(int)
    df["POS"] = df["POS"].astype(int)
    win = df[(df["CHR"] == chrom) & (df["POS"] >= window_start) & (df["POS"] <= window_end)].copy()
    if len(win) == 0:
        return None
    win["SNP"] = (win["CHR"].astype(str) + ":" + win["POS"].astype(str) + ":"
                  + win["A1"].astype(str) + ":" + win["A2"].astype(str))
    win["stat"] = win["BETA"] / win["SE"]
    out = win[["CHR", "SNP", "POS", "A1", "A2", "BETA", "SE", "stat", "PVAL"]].copy()
    out.columns = ["chr", "snp", "bp", "A1", "A2", "BETA", "SE", "stat", "PVAL"]
    out = out.drop_duplicates(subset=["snp"])
    print(f"    {len(out)} variants in window chr{chrom}:{window_start}-{window_end}")
    return out


# ---------------------------------------------------------------------------
# Reference panel preparation
# ---------------------------------------------------------------------------

CHR_COL, SNP_COL, BP_COL = 1, 2, 3
A1_COL, A2_COL = 4, 5
EFF_COL, SE_COL = 6, 7
PVAL_COL = 9


def rename_bim_to_chrpos(bim_path):
    with open(bim_path) as f:
        lines = f.readlines()
    seen = {}
    new_lines = []
    for line in lines:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 6:
            new_lines.append(line); continue
        chrom_b, snp_id, cm, pos, a1, a2 = parts[:6]
        snp_id = f"{chrom_b}:{pos}:{a1}:{a2}"
        if snp_id in seen:
            seen[snp_id] += 1
            snp_id = f"{snp_id}_{seen[snp_id]}"
        else:
            seen[snp_id] = 1
        new_lines.append("\t".join([chrom_b, snp_id, cm, pos, a1, a2]) + "\n")
    with open(bim_path, "w") as f:
        f.writelines(new_lines)


def prepare_block_ref(ld_dir, chrom, win_start, win_end, tmpdir, locus_id, anc):
    """EUR-style: copy/merge per-block .bed/.bim/.fam into tmpdir; reuse adjacent .ld."""
    blocks = find_ld_block_prefix(ld_dir, chrom, win_start, win_end)
    print(f"  {anc} LD blocks: {len(blocks)}")
    ref_prefix = os.path.join(tmpdir, f"{locus_id}_{anc.lower()}_ref")
    ld_prefix  = os.path.join(tmpdir, f"{locus_id}_{anc.lower()}_ld")

    if len(blocks) == 1:
        for ext in [".bed", ".bim", ".fam"]:
            shutil.copy2(blocks[0] + ext, ref_prefix + ext)
    else:
        merge_list_file = os.path.join(tmpdir, f"{anc.lower()}_merge_list.txt")
        with open(merge_list_file, "w") as f:
            for b in blocks[1:]:
                f.write(f"{b}\n")
        cmd = [PLINK_BIN, "--bfile", blocks[0], "--merge-list", merge_list_file,
               "--make-bed", "--out", ref_prefix, "--allow-extra-chr"]
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            print(f"  PLINK merge warning ({anc}): {r.stderr[-300:]}")
            for ext in [".bed", ".bim", ".fam"]:
                shutil.copy2(blocks[0] + ext, ref_prefix + ext)
    rename_bim_to_chrpos(ref_prefix + ".bim")
    return ref_prefix, ld_prefix


def prepare_chr_ref(ld_dir, chr_basename, chrom, win_start, win_end, tmpdir, locus_id, anc):
    """Non-EUR-style: extract window from chr-level PLINK; SuSiEX recomputes LD on-the-fly."""
    chr_prefix = os.path.join(ld_dir, chr_basename.format(chrom=chrom))
    for ext in [".bed", ".bim", ".fam"]:
        if not os.path.exists(chr_prefix + ext):
            raise FileNotFoundError(f"{anc} chr-level PLINK missing: {chr_prefix}{ext}")
    ref_prefix = os.path.join(tmpdir, f"{locus_id}_{anc.lower()}_ref")
    ld_prefix  = os.path.join(tmpdir, f"{locus_id}_{anc.lower()}_ld")
    cmd = [PLINK_BIN, "--bfile", chr_prefix,
           "--chr", str(chrom), "--from-bp", str(win_start), "--to-bp", str(win_end),
           "--make-bed", "--out", ref_prefix, "--allow-extra-chr"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"PLINK extract failed for {anc}:\n{r.stderr[-500:]}")
    rename_bim_to_chrpos(ref_prefix + ".bim")
    return ref_prefix, ld_prefix


def prepare_ancestry_ref(anc, chrom, win_start, win_end, tmpdir, locus_id):
    ld_dir = _resolve_ld_dir(anc)
    cfg = ANCESTRY_REGISTRY[anc]
    if cfg["ref_kind"] == "block":
        return prepare_block_ref(ld_dir, chrom, win_start, win_end, tmpdir, locus_id, anc)
    return prepare_chr_ref(ld_dir, cfg["chr_basename"], chrom, win_start, win_end,
                           tmpdir, locus_id, anc)


# ---------------------------------------------------------------------------
# Main per-locus runner
# ---------------------------------------------------------------------------

MIN_SNPS = 50  # minimum variants in window to keep an ancestry arm


def run_susiex_for_locus(locus, output_dir, ancestries, threads=4):
    locus_id   = locus["locus_id"]
    chrom      = int(locus["chr"])
    win_start  = int(locus["window_start"])
    win_end    = int(locus["window_end"])
    trait_pair = locus["trait_pair"]
    # cohort: "MVP" (within-MVP multi-ancestry) vs legacy enzyme UKBB<->BBJ.
    cohort     = str(locus["cohort"]) if "cohort" in locus.index else "legacy"

    print(f"\n{'='*60}")
    print(f"Locus: {locus_id}  Region: chr{chrom}:{win_start}-{win_end}")
    print(f"Trait: {trait_pair}  Cohort: {cohort}  Ancestries requested: {','.join(ancestries)}")
    os.makedirs(output_dir, exist_ok=True)

    summary_file = os.path.join(output_dir, f"{locus_id}.summary")
    if os.path.exists(summary_file):
        print(f"  Already completed: {summary_file} — skipping")
        return

    # 1. Load sumstats per ancestry; drop arms with no variants
    ss_per_anc = {}
    for anc in ancestries:
        study, path = _build_sumstats_path(anc, trait_pair, cohort)
        if path is None or not os.path.exists(path):
            print(f"  [{anc}] sumstats missing ({study}); dropping arm")
            continue
        ss = load_sumstats_window(study, path, chrom, win_start, win_end)
        if ss is None or len(ss) < MIN_SNPS:
            n = 0 if ss is None else len(ss)
            print(f"  [{anc}] only {n} variants in window (<{MIN_SNPS}); dropping arm")
            continue
        ss_per_anc[anc] = (study, ss)

    if len(ss_per_anc) < 2:
        print(f"  Fewer than 2 viable ancestry arms ({list(ss_per_anc)}) — skipping")
        return

    active = list(ss_per_anc.keys())

    with tempfile.TemporaryDirectory(prefix=f"susiex_{locus_id}_") as tmpdir:
        # 2. Prepare ref/ld panel per active ancestry
        ref_prefixes = []
        ld_prefixes  = []
        for anc in active:
            try:
                rp, lp = prepare_ancestry_ref(anc, chrom, win_start, win_end, tmpdir, locus_id)
            except (FileNotFoundError, ValueError, RuntimeError) as e:
                print(f"  [{anc}] LD prep failed: {e} — dropping arm")
                ss_per_anc.pop(anc); continue
            ref_prefixes.append(rp)
            ld_prefixes.append(lp)

        active = list(ss_per_anc.keys())  # may have shrunk
        if len(active) < 2:
            print(f"  Fewer than 2 viable arms after LD prep — skipping")
            return

        # 3. Write per-ancestry sumstats temp files (preserving ancestry order)
        ss_files = []
        n_gwas_list = []
        study_list = []
        for anc, rp, lp in zip(active, ref_prefixes, ld_prefixes):
            study, ss = ss_per_anc[anc]
            ss_file = os.path.join(tmpdir, f"{locus_id}_{anc}.txt")
            ss.to_csv(ss_file, sep="\t", index=False)
            ss_files.append(ss_file)
            n_gwas_list.append(_lookup_n_gwas(anc, trait_pair, study, cohort))
            study_list.append(study)

        print(f"  Active arms: {','.join(f'{a}={s}' for a, s in zip(active, study_list))}")

        # 4. Build SuSiEX command (comma-separated lists across ancestries)
        n_anc = len(active)
        cmd = [
            SUSIEX_BIN,
            "--sst_file=" + ",".join(ss_files),
            "--n_gwas="   + ",".join(str(n) for n in n_gwas_list),
            "--ref_file=" + ",".join(ref_prefixes),
            "--ld_file="  + ",".join(ld_prefixes),
            f"--out_dir={output_dir}",
            f"--out_name={locus_id}",
            f"--chr={chrom}",
            f"--bp={win_start},{win_end}",
            "--chr_col=" + ",".join([str(CHR_COL)]  * n_anc),
            "--snp_col=" + ",".join([str(SNP_COL)]  * n_anc),
            "--bp_col="  + ",".join([str(BP_COL)]   * n_anc),
            "--a1_col="  + ",".join([str(A1_COL)]   * n_anc),
            "--a2_col="  + ",".join([str(A2_COL)]   * n_anc),
            "--eff_col=" + ",".join([str(EFF_COL)]  * n_anc),
            "--se_col="  + ",".join([str(SE_COL)]   * n_anc),
            "--pval_col=" + ",".join([str(PVAL_COL)] * n_anc),
            f"--plink={PLINK_BIN}",
            "--keep-ambig=True",
            "--mult-step=True",
            f"--threads={threads}",
            "--level=0.95",
            "--min_purity=0.5",
            "--pval_thresh=1e-5",
            "--maf=0.005",
        ]
        print(f"  Running SuSiEX with {n_anc} ancestries: {','.join(active)}")
        r = subprocess.run(cmd, capture_output=True, text=True, cwd=FM_DIR)
        if r.returncode != 0:
            print(f"  SuSiEX STDERR:\n{r.stderr[-2000:]}", file=sys.stderr)
            raise RuntimeError(f"SuSiEX failed for {locus_id} (exit {r.returncode})")
        if r.stdout:
            print(r.stdout[-1000:])
    print(f"  Done. Results in {output_dir}/{locus_id}.*")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--locus-row", type=int, required=True)
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()

    ancestries = get_active_ancestries()

    # SHARED_LOCI_FILE overrides the default enzyme loci table (used by the MVP
    # runner to point at results/susiex_mvp/shared_loci.csv).
    shared_loci_path = os.environ.get(
        "SHARED_LOCI_FILE",
        os.path.join(FM_DIR, "results/susiex/shared_loci.csv"))
    if not os.path.exists(shared_loci_path):
        print(f"ERROR: shared_loci.csv not found at {shared_loci_path}", file=sys.stderr)
        sys.exit(1)
    shared = pd.read_csv(shared_loci_path)
    if args.locus_row < 1 or args.locus_row > len(shared):
        print(f"ERROR: locus-row {args.locus_row} out of range (1..{len(shared)})", file=sys.stderr)
        sys.exit(1)

    locus = shared.iloc[args.locus_row - 1]
    suffix = os.environ.get("SUSIEX_RESULTS_SUFFIX", "")
    output_dir = os.path.join(FM_DIR, f"results/susiex{suffix}", locus["trait_pair"])
    os.makedirs(output_dir, exist_ok=True)

    try:
        run_susiex_for_locus(locus, output_dir, ancestries, threads=args.threads)
    except RuntimeError as e:
        print(f"WARNING: {e}", file=sys.stderr)
        print(f"Skipping {locus['locus_id']}")
        sys.exit(0)


if __name__ == "__main__":
    main()
