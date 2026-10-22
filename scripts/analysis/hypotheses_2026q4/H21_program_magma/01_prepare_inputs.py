#!/usr/bin/env python
"""H21 step 1: build MAGMA inputs. No outcome is read here.

Writes to <out>/inputs:
  gene_covar.tsv            Entrez, hep_spec, mean_expr (donor-level pseudobulk)
  program_sets.txt          one row per program: uid then Entrez IDs
  program_set_manifest.tsv  mapping counts per program
  random_sets.txt           matched random sets, seeded
  all_sets.txt              programs + random sets
  mhc_genes.txt             Entrez IDs in chr6:25-34 Mb (GRCh37)
  gwas/<trait>.pval.tsv     SNP, P, N
  gwas_prep_audit.tsv       row counts per trait
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

SEED = 20261007
N_RANDOM_PER_PROGRAM = 10
ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
GENE_LOC = ROOT / "tools/magma/ref/NCBI37.3.gene.loc"
BIM = ROOT / "tools/magma/ref/g1000_eur.bim"
SCDIR = ROOT / "Analysis/SingleCell/results_gpu_v2/disease_signatures/celltype_fstage_pseudobulk"
HOT = ROOT / "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot"
REG = HOT / "program_registry_v2.tsv"
MEM = HOT / "program_membership_v2.tsv"
PAIRING = {ds: ROOT / f"data/{ds}/metadata/donor_pairing.csv"
           for ds in ("GSE244832", "GSE185477", "GSE202379", "GSE136103")}
LINEAGES = ["Hepatocytes", "Fibroblasts", "Cholangiocytes", "Macrophages",
            "T_cells", "Endothelial_cells"]
MR = ROOT / "GWAS/MR_Data"


def log(msg: str) -> None:
    print(f"[H21 prep] {msg}", flush=True)


def load_gene_loc() -> pd.DataFrame:
    g = pd.read_csv(GENE_LOC, sep="\t", header=None,
                    names=["entrez", "chr", "start", "end", "strand", "symbol"],
                    dtype={"entrez": str, "chr": str})
    g["length"] = g["end"] - g["start"] + 1
    return g


def build_covariates(gl: pd.DataFrame) -> pd.DataFrame:
    """Donor-level mean log2(CPM+1) per lineage, same recipe as the adopted
    label-versus-content analysis (plc-20260907T105813Z/label_vs_content.py)."""
    srr2donor = {}
    for ds, fp in PAIRING.items():
        dp = pd.read_csv(fp, dtype=str)
        for _, r in dp.iterrows():
            for s in str(r["rna_srrs"]).split(";"):
                if s and s != "nan":
                    srr2donor[s.strip()] = f"{ds}_{r['donor_id']}"
    prof, ndon = {}, {}
    for ln in LINEAGES:
        m = pd.read_csv(SCDIR / f"pseudobulk_{ln}.tsv.gz", sep="\t", index_col=0)
        m = m.groupby(level=0).sum()
        cols = pd.Series(m.columns).map(srr2donor).fillna(pd.Series(m.columns)).values
        m = m.T.groupby(cols).sum().T
        cpm = m.div(m.sum(axis=0).replace(0, np.nan), axis=1) * 1e6
        prof[ln] = np.log2(cpm + 1).mean(axis=1)
        ndon[ln] = int(m.shape[1])
        log(f"lineage {ln}: {ndon[ln]} donor-level profiles")
    P = pd.DataFrame(prof).dropna(how="any")
    sd = P.std(axis=1)
    z = P.sub(P.mean(axis=1), axis=0).div(sd.replace(0, np.nan), axis=0)
    cov = pd.DataFrame({"symbol": P.index,
                        "hep_spec": z["Hepatocytes"].fillna(0.0).values,
                        "mean_expr": P.mean(axis=1).values})
    sym2entrez = gl.drop_duplicates("symbol").set_index("symbol")["entrez"]
    cov["entrez"] = cov["symbol"].map(sym2entrez)
    n_all = len(cov)
    cov = cov.dropna(subset=["entrez"]).drop_duplicates("entrez")
    log(f"covariates: {n_all} atlas genes, {len(cov)} mapped to NCBI37.3 Entrez")
    return cov[["entrez", "symbol", "hep_spec", "mean_expr"]], ndon


def build_program_sets(gl: pd.DataFrame, universe: set[str]):
    reg = pd.read_csv(REG, sep="\t", dtype=str)
    mem = pd.read_csv(MEM, sep="\t", dtype=str)
    assert reg["program_uid"].nunique() == 117, reg["program_uid"].nunique()
    symbols = set(gl["symbol"])
    sym2entrez = gl.drop_duplicates("symbol").set_index("symbol")["entrez"].to_dict()

    def pick(row):
        for c in ("mapped_symbol", "canonical_gene", "source_gene"):
            v = row[c]
            if isinstance(v, str) and v in symbols:
                return sym2entrez[v]
        return None

    mem["entrez"] = mem.apply(pick, axis=1)
    sets, rows = {}, []
    for uid, sub in mem.groupby("program_uid", sort=False):
        ids = sorted(set(sub["entrez"].dropna()))
        in_uni = [i for i in ids if i in universe]
        sets[uid] = in_uni
        r = reg.loc[reg["program_uid"] == uid].iloc[0]
        rows.append({"program_uid": uid, "cell_type": r["cell_type"], "module": r["module"],
                     "module_name": r["module_name"], "n_source_genes": len(sub),
                     "n_mapped_entrez": len(ids), "n_in_candidate_universe": len(in_uni)})
    man = pd.DataFrame(rows)
    assert len(man) == 117
    return sets, man


def matched_random_sets(sets: dict, cov: pd.DataFrame, gl: pd.DataFrame):
    rng = np.random.default_rng(SEED)
    u = cov.merge(gl[["entrez", "length"]], on="entrez", how="left")
    u["q_spec"] = pd.qcut(u["hep_spec"].rank(method="first"), 5, labels=False)
    u["q_len"] = pd.qcut(u["length"].rank(method="first"), 5, labels=False)
    u["stratum"] = u["q_spec"] * 5 + u["q_len"]
    by_stratum = {s: g["entrez"].to_numpy() for s, g in u.groupby("stratum")}
    stratum_of = u.set_index("entrez")["stratum"].to_dict()
    out = {}
    for uid in sorted(sets):
        strata = pd.Series([stratum_of[g] for g in sets[uid]]).value_counts()
        for k in range(N_RANDOM_PER_PROGRAM):
            draw = []
            for s, n in strata.sort_index().items():
                pool = by_stratum[s]
                draw.extend(rng.choice(pool, size=min(n, len(pool)), replace=False))
            out[f"RND__{uid}__{k}"] = sorted(draw)
    return out


def write_sets(path: Path, sets: dict) -> None:
    with open(path, "w") as fh:
        for name, ids in sets.items():
            fh.write(name + " " + " ".join(ids) + "\n")


def finish_gwas(df: pd.DataFrame, trait: str, outdir: Path, n_in: int) -> dict:
    df = df.copy()
    df["P"] = pd.to_numeric(df["P"], errors="coerce")
    df["N"] = pd.to_numeric(df["N"], errors="coerce")
    df = df.dropna(subset=["SNP", "P", "N"])
    df = df[df["SNP"].astype(str).str.startswith("rs")]
    df = df[(df["P"] >= 0) & (df["P"] <= 1)].copy()
    df.loc[df["P"] < 1e-300, "P"] = 1e-300
    df = df.drop_duplicates("SNP", keep=False)
    df["N"] = df["N"].round().astype(int)
    fp = outdir / f"{trait}.pval.tsv"
    df[["SNP", "P", "N"]].to_csv(fp, sep="\t", index=False)
    log(f"{trait}: {n_in} rows in, {len(df)} written, N median {int(df['N'].median())}")
    return {"trait": trait, "rows_in": n_in, "rows_written": len(df),
            "N_min": int(df["N"].min()), "N_median": int(df["N"].median()),
            "N_max": int(df["N"].max()), "file": str(fp)}


def load_bim() -> pd.DataFrame:
    bim = pd.read_csv(BIM, sep="\t", header=None, usecols=[0, 1, 3, 4, 5],
                      names=["chr", "SNP", "pos", "a1", "a2"],
                      dtype={"chr": str, "SNP": str, "pos": np.int64, "a1": str, "a2": str})
    bim["alleles"] = np.where(bim["a1"] < bim["a2"], bim["a1"] + "/" + bim["a2"],
                              bim["a2"] + "/" + bim["a1"])
    return bim[["chr", "pos", "alleles", "SNP"]]


def join_position(df: pd.DataFrame, bim: pd.DataFrame) -> pd.DataFrame:
    """df has chr, pos, x1, x2. Match on GRCh37 position and unordered allele pair."""
    df["alleles"] = np.where(df["x1"] < df["x2"], df["x1"] + "/" + df["x2"],
                             df["x2"] + "/" + df["x1"])
    return df.merge(bim, on=["chr", "pos", "alleles"], how="inner")


def prep_gwas(outdir: Path) -> pd.DataFrame:
    outdir.mkdir(parents=True, exist_ok=True)
    audit = []

    d = pd.read_csv(MR / "GCST90267352_PDFF_Pazoki2022.tsv.gz", sep="\t",
                    usecols=["variant_id", "p_value", "n"])
    audit.append(finish_gwas(d.rename(columns={"variant_id": "SNP", "p_value": "P", "n": "N"}),
                             "pdff", outdir, len(d)))

    d = pd.read_csv(MR / "Ghouse_Cirrhosis/GCST90319877_harmonised_hg38.tsv.gz", sep="\t",
                    usecols=["rsid", "p_value", "samplesize"])
    audit.append(finish_gwas(d.rename(columns={"rsid": "SNP", "p_value": "P", "samplesize": "N"}),
                             "cirrhosis", outdir, len(d)))

    d = pd.read_csv(MR / "Ghodsian_2021_NAFLD_harmonised.tsv.gz", sep="\t",
                    usecols=["hm_rsid", "p_value"])
    d["N"] = 8434 + 770180
    audit.append(finish_gwas(d.rename(columns={"hm_rsid": "SNP", "p_value": "P"}),
                             "nafld", outdir, len(d)))

    d = pd.read_csv(MR / "Yengo_2018_BMI_GWAS.txt.gz", sep=r"\s+", usecols=["SNP", "P", "N"])
    audit.append(finish_gwas(d, "bmi", outdir, len(d)))
    del d

    bim = load_bim()
    log(f"reference bim: {len(bim)} SNPs")

    parts, n_in = [], 0
    usecols = ["chr", "pos", "ref", "alt", "neglog10_pval_EUR", "low_confidence_EUR"]
    for ch in pd.read_csv(MR / "PanUKBB/biomarkers-30620-ALT-both_sexes-irnt.tsv.bgz",
                          sep="\t", compression="gzip", usecols=usecols,
                          dtype={"chr": str, "ref": str, "alt": str, "low_confidence_EUR": str},
                          chunksize=2_000_000):
        n_in += len(ch)
        ch = ch[ch["neglog10_pval_EUR"].notna()]
        ch = ch[ch["low_confidence_EUR"].str.lower() != "true"]
        ch = ch.rename(columns={"ref": "x1", "alt": "x2"})
        m = join_position(ch, bim)
        m["P"] = np.power(10.0, -m["neglog10_pval_EUR"].clip(upper=300))
        parts.append(m[["SNP", "P"]])
    d = pd.concat(parts, ignore_index=True)
    d["N"] = 400822
    audit.append(finish_gwas(d, "alt", outdir, n_in))

    d = pd.read_csv(ROOT / "GWAS/finemapping/data/sumstats/MVP_Cirrhosis_EUR_reformatted_hg19.tsv",
                    sep="\t", usecols=["chromosome", "position", "allele1", "allele2", "pval"],
                    dtype={"chromosome": str, "allele1": str, "allele2": str})
    n_in = len(d)
    d = d.rename(columns={"chromosome": "chr", "position": "pos", "allele1": "x1",
                          "allele2": "x2", "pval": "P"})
    d["chr"] = d["chr"].str.replace("chr", "", regex=False)
    d["x1"] = d["x1"].str.upper()
    d["x2"] = d["x2"].str.upper()
    m = join_position(d, bim)
    m["N"] = 447231
    audit.append(finish_gwas(m[["SNP", "P", "N"]], "cirrhosis_mvp", outdir, n_in))
    return pd.DataFrame(audit)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--skip-gwas", action="store_true")
    a = ap.parse_args()
    inp = Path(a.out) / "inputs"
    inp.mkdir(parents=True, exist_ok=True)
    (inp / "SEED.txt").write_text(f"{SEED}\n")

    gl = load_gene_loc()
    mhc = gl[(gl["chr"] == "6") & (gl["end"] >= 25_000_000) & (gl["start"] <= 34_000_000)]
    mhc["entrez"].to_csv(inp / "mhc_genes.txt", index=False, header=False)
    log(f"MHC genes excluded: {len(mhc)}")

    cov, ndon = build_covariates(gl)
    cov[["entrez", "hep_spec", "mean_expr"]].rename(columns={"entrez": "GENE"}).to_csv(
        inp / "gene_covar.tsv", sep="\t", index=False, float_format="%.6g")
    cov.to_csv(inp / "gene_covar_with_symbols.tsv", sep="\t", index=False)
    pd.Series(ndon, name="n_donor_profiles").to_csv(inp / "donors_per_lineage.tsv", sep="\t")

    cand = cov[~cov["entrez"].isin(set(mhc["entrez"]))]
    universe = set(cand["entrez"])
    log(f"candidate universe (covariates, non-MHC): {len(universe)} genes")

    sets, man = build_program_sets(gl, universe)
    man.to_csv(inp / "program_set_manifest.tsv", sep="\t", index=False)
    log(f"programs: {len(sets)}; mapped size min/median/max "
        f"{man.n_in_candidate_universe.min()}/{int(man.n_in_candidate_universe.median())}/"
        f"{man.n_in_candidate_universe.max()}; below 10 genes: "
        f"{int((man.n_in_candidate_universe < 10).sum())}")
    rnd = matched_random_sets(sets, cand, gl)
    write_sets(inp / "program_sets.txt", sets)
    write_sets(inp / "random_sets.txt", rnd)
    write_sets(inp / "all_sets.txt", {**sets, **rnd})
    log(f"random sets: {len(rnd)} (seed {SEED})")

    if not a.skip_gwas:
        prep_gwas(inp / "gwas").to_csv(inp / "gwas_prep_audit.tsv", sep="\t", index=False)
    log("done")


if __name__ == "__main__":
    sys.exit(main())
