# DIG1 job 73 (remediated 2026-11-18, N05 + N14):
# - Reads COLOC from EQTL-DEDUP-PSD-20261118/results/susie_coloc_all_gwas_refit_20261118.csv
# - Fixes the AST-only orientation bug (N14):
#   Do NOT blindly orient every SNP by sign(z_MVP_ALT_EUR) before locus grouping.
#   Instead:
#   1. Merge C with raw (allele-aligned to MVP_ALT_EUR a1/a2, NOT sign-flipped yet) W_raw, group genes into 500 kb loci.
#   2. Per locus, define defining_enzyme:
#      "MVP_ALT_EUR" if any colocalized gene in the locus has gwas_name == "MVP_ALT_EUR" AND max(|z_MVP_ALT_EUR|) >= 5.45
#      (or if no gene in the locus has an MVP_AST_EUR colocalization passing h4 > 0.8),
#      else "MVP_AST_EUR" (with fallback to "MVP_ALT_EUR" only if all z_MVP_AST_EUR are missing).
#   3. Pick the locus lead SNP as the row in that locus with maximum |z_defining|.
#   4. Orient ALL z_* traits for that locus to the enzyme-raising allele of defining_enzyme: sgn = +1 if z_defining >= 0 else -1.
#      Record defining_enzyme, z_defining_enzyme, enzyme_raising_allele, and alt_raising_allele (= enzyme_raising_allele for compatibility).
#   5. Apply the class rule using z_defining = z_defining_enzyme.
# Usage: python -I 73_release_vs_disease_loci_20261118.py <new_out_dir>
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

out = sys.argv[1]
assert not os.path.exists(out), f"{out} already exists; refusing to overwrite"
os.makedirs(out)

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
RUN = f"{ROOT}/Analysis/Hypotheses/DIG1-20261024T113951Z"
SS = f"{ROOT}/GWAS/finemapping/data/sumstats"
COLOC = f"{ROOT}/Analysis/Hypotheses/EQTL-DEDUP-PSD-20261118/results/susie_coloc_all_gwas_refit_20261118.csv"
DEG = f"{ROOT}/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv"
TRAITS = [
    "MVP_ALT_EUR", "MVP_AST_EUR", "MVP_NAFLD_EUR", "MVP_Cirrhosis_EUR", "MVP_ChronLiver_EUR",
    "MVP_Albumin_EUR", "MVP_Platelet_EUR", "UKBB_ALT", "UKBB_AST", "UKBB_GGT",
    "FinnGen_NAFLD", "FinnGen_NASH", "Pazoki_PDFF",
]
BENCH = {
    "PNPLA3_rs738409": (22, 44324727),
    "TM6SF2_rs58542926": (19, 19379549),
    "GCKR_rs1260326": (2, 27730940),
}

# 1. universe
C_raw = pd.read_csv(
    COLOC,
    usecols=["gwas_name", "gene", "ensembl", "PP.H4.abf", "PP.H4.susie", "top_snp"],
    low_memory=False,
)
C_raw = C_raw[C_raw.gwas_name.isin(["MVP_ALT_EUR", "MVP_AST_EUR"])].copy()
C_raw["PP.H4.abf"] = pd.to_numeric(C_raw["PP.H4.abf"], errors="coerce")
C_raw["PP.H4.susie"] = pd.to_numeric(C_raw["PP.H4.susie"], errors="coerce")
C_raw["h4"] = C_raw[["PP.H4.abf", "PP.H4.susie"]].max(axis=1, skipna=True)
C_pass = C_raw[C_raw.h4 > 0.8].copy()

gene_gwas = C_pass.groupby("gene")["gwas_name"].agg(set).to_dict()
C = C_pass.copy()
C["pref"] = (C.gwas_name != "MVP_ALT_EUR").astype(int)
C = C.sort_values(["gene", "pref", "h4"], ascending=[True, True, False]).drop_duplicates("gene")
C["has_alt_coloc"] = C["gene"].map(lambda g: "MVP_ALT_EUR" in gene_gwas[g])
C["has_ast_coloc"] = C["gene"].map(lambda g: "MVP_AST_EUR" in gene_gwas[g])

sp = C.top_snp.astype(str).str.replace("chr", "").str.split(":", expand=True)
C["chr"] = pd.to_numeric(sp[0], errors="coerce")
C["pos"] = pd.to_numeric(sp[1].str.split("_").str[0], errors="coerce")
print("coloc genes:", len(C), "unparsed top_snp:", int(C.chr.isna().sum()))
C = C.dropna(subset=["chr", "pos"]).astype({"chr": int, "pos": int})
want = set(zip(C.chr, C.pos)) | set(BENCH.values())


# 2. extract raw z (allele-aligned to MVP_ALT_EUR a1/a2, NOT sign-flipped yet)
def extract(t):
    keep = []
    for ch in pd.read_csv(
        f"{SS}/{t}_reformatted_hg19.tsv",
        sep="\t",
        chunksize=2_000_000,
        usecols=["chromosome", "position", "allele1", "allele2", "beta", "se"],
        dtype={"chromosome": str},
    ):
        ch["chromosome"] = pd.to_numeric(ch.chromosome.str.replace("chr", ""), errors="coerce")
        ch = ch.dropna(subset=["chromosome"])
        m = [(int(a), int(b)) in want for a, b in zip(ch.chromosome, ch.position)]
        keep.append(ch[m])
    d = pd.concat(keep, ignore_index=True).rename(columns={"chromosome": "chr", "position": "pos"})
    d["chr"] = d.chr.astype(int)
    d["z"] = d.beta / d.se
    d["allele1"] = d.allele1.str.upper()
    d["allele2"] = d.allele2.str.upper()
    return d.drop_duplicates(["chr", "pos", "allele1", "allele2"])


Z = {t: extract(t) for t in TRAITS}
ref = Z["MVP_ALT_EUR"].drop_duplicates(["chr", "pos"], keep=False)[["chr", "pos", "allele1", "allele2", "z"]]
ref = ref.rename(columns={"allele1": "a1", "allele2": "a2", "z": "z_MVP_ALT_EUR"})
W_raw = ref.copy()
for t in TRAITS[1:]:
    d = Z[t].merge(ref[["chr", "pos", "a1", "a2"]], on=["chr", "pos"])
    same = (d.allele1 == d.a1) & (d.allele2 == d.a2)
    swap = (d.allele1 == d.a2) & (d.allele2 == d.a1)
    d = d[same | swap].copy()
    d.loc[swap[same | swap].values, "z"] *= -1
    d = d.drop_duplicates(["chr", "pos"], keep=False)
    W_raw = W_raw.merge(d[["chr", "pos", "z"]].rename(columns={"z": f"z_{t}"}), on=["chr", "pos"], how="left")

zc = [c for c in W_raw.columns if c.startswith("z_")]

# 3. merge C with raw W_raw and group genes into 500 kb loci BEFORE orientation
G_raw = C.merge(W_raw, on=["chr", "pos"], how="left")
print("coloc genes with MVP ALT z at lead:", int(G_raw.z_MVP_ALT_EUR.notna().sum()), "of", len(G_raw))
G_raw = G_raw.dropna(subset=["z_MVP_ALT_EUR"]).sort_values(["chr", "pos"]).reset_index(drop=True)

lid, prev, k = [], (None, None), -1
for ch, ps in zip(G_raw.chr, G_raw.pos):
    if ch != prev[0] or ps - prev[1] > 500_000:
        k += 1
    lid.append(k)
    prev = (ch, ps)
G_raw["locus"] = lid

hs = pd.read_csv(f"{RUN}/out_hep_share/gene_hep_share.tsv.gz", sep="\t")
hs = hs.groupby("gene").agg(mn=("median_hep_share", "min"), k=("cohort", "nunique"))
HEP = set(hs[(hs.mn >= 0.85) & (hs.k == 3)].index)
STATE = set(open(f"{RUN}/out_state_magma/state.sets").readline().split()[1:])
deg = pd.read_csv(DEG)
deg["ens"] = deg.gene.str.split(".").str[0]
UP = set(deg[(deg.padj < 0.05) & (deg.logFC > 0)].ens)
G_raw["state_up"] = G_raw.ensembl.str.split(".").str[0].isin(STATE)
G_raw["hep_made"] = G_raw.gene.isin(HEP)
G_raw["bulk_up"] = G_raw.ensembl.str.split(".").str[0].isin(UP)

# Per-locus defining_enzyme, lead SNP selection, and orientation to enzyme-raising allele
locus_rows = []
gene_rows = []
for loc_id, grp in G_raw.groupby("locus", sort=True):
    grp = grp.copy()
    has_alt = bool((grp["gwas_name"] == "MVP_ALT_EUR").any())
    has_ast = bool(grp["has_ast_coloc"].any())
    max_abs_alt = grp["z_MVP_ALT_EUR"].abs().max()
    max_abs_ast = grp["z_MVP_AST_EUR"].abs().max()

    if has_alt and pd.notna(max_abs_alt) and max_abs_alt >= 5.45:
        def_enz = "MVP_ALT_EUR"
    elif not has_ast:
        def_enz = "MVP_ALT_EUR"
    elif pd.notna(max_abs_ast) and (not has_alt or max_abs_ast >= 5.45):
        def_enz = "MVP_AST_EUR"
    elif pd.notna(max_abs_ast) and (pd.isna(max_abs_alt) or max_abs_ast > max_abs_alt):
        def_enz = "MVP_AST_EUR"
    else:
        def_enz = "MVP_ALT_EUR"

    zcol = f"z_{def_enz}"
    if grp[zcol].isna().all():
        def_enz = "MVP_ALT_EUR"
        zcol = "z_MVP_ALT_EUR"

    lead_idx = grp[zcol].abs().idxmax()
    lead_row = grp.loc[lead_idx].copy()
    z_def_raw = float(lead_row[zcol])
    sgn = 1.0 if z_def_raw >= 0 else -1.0

    for c in zc:
        lead_row[c] = lead_row[c] * sgn
    enz_allele = lead_row["a1"] if sgn > 0 else lead_row["a2"]

    rec = {
        "locus": int(loc_id),
        "chr": int(lead_row["chr"]),
        "pos": int(lead_row["pos"]),
        "a1": lead_row["a1"],
        "a2": lead_row["a2"],
        "defining_enzyme": def_enz,
        "z_defining_enzyme": float(lead_row[zcol]),
        "enzyme_raising_allele": enz_allele,
        "alt_raising_allele": enz_allele,
    }
    for c in zc:
        rec[c] = lead_row[c]
    rec["genes"] = ",".join(sorted(set(grp["gene"])))
    rec["n_genes"] = int(grp["gene"].nunique())
    rec["state_up"] = bool(grp["state_up"].any())
    rec["hep_made"] = bool(grp["hep_made"].any())
    rec["bulk_up"] = bool(grp["bulk_up"].any())
    locus_rows.append(rec)

    # Orient each gene row in coloc_genes_lead_z.tsv to the enzyme-raising allele of def_enz at that gene's SNP
    for _, gr in grp.iterrows():
        gr = gr.copy()
        zv = gr[zcol] if pd.notna(gr[zcol]) else gr["z_MVP_ALT_EUR"]
        gsgn = 1.0 if float(zv) >= 0 else -1.0
        for c in zc:
            gr[c] = gr[c] * gsgn
        g_allele = gr["a1"] if gsgn > 0 else gr["a2"]
        gr["defining_enzyme"] = def_enz
        gr["z_defining_enzyme"] = float(gr[zcol]) if pd.notna(gr[zcol]) else np.nan
        gr["enzyme_raising_allele"] = g_allele
        gr["alt_raising_allele"] = g_allele
        gene_rows.append(gr)

G = pd.DataFrame(gene_rows)
G.to_csv(f"{out}/coloc_genes_lead_z.tsv", sep="\t", index=False)

L = pd.DataFrame(locus_rows)

# Benchmarks oriented to their ALT-raising allele
B = W_raw[[(a, b) in set(BENCH.values()) for a, b in zip(W_raw.chr, W_raw.pos)]].copy()
bsgn = np.where(B["z_MVP_ALT_EUR"] >= 0, 1.0, -1.0)
for c in zc:
    B[c] = B[c] * bsgn
B["alt_raising_allele"] = np.where(bsgn > 0, B.a1, B.a2)
B["genes"] = [next(k for k, v in BENCH.items() if v == (a, b)) for a, b in zip(B.chr, B.pos)]
print("benchmarks:\n", B[["genes"] + zc].round(1).to_string(index=False))

# 4. classes using z_defining_enzyme
kk = float(
    np.mean(
        [
            B.loc[B.genes == g, "z_MVP_NAFLD_EUR"].iloc[0] / B.loc[B.genes == g, "z_MVP_ALT_EUR"].iloc[0]
            for g in ["PNPLA3_rs738409", "TM6SF2_rs58542926"]
        ]
    )
)
zd, zn = L.z_defining_enzyme, L.z_MVP_NAFLD_EUR
L["cls"] = np.select(
    [zn >= 3.29, zn <= -3.29, (zd >= 5.45) & (zn.abs() < 1.96) & (kk * zd - zn >= 3)],
    ["disease", "opposite", "release"],
    "intermediate",
)
L.loc[zn.isna(), "cls"] = "no_disease_z"
L["ratio_nafld_enzyme"] = zn / zd
L["ratio_nafld_alt"] = zn / L.z_MVP_ALT_EUR
L.to_csv(f"{out}/loci.tsv", sep="\t", index=False)

print(f"k (benchmark z_NAFLD/z_ALT) = {kk:.3f}; loci = {len(L)}")
print("defining_enzyme counts:")
print(L.defining_enzyme.value_counts().to_string())
print("\nclass counts:")
print(L.cls.value_counts().to_string())

# Explicit diagnostics for AST-defined loci and key named loci
pd.set_option("display.width", 250)
pd.set_option("display.max_colwidth", 60)
diag_cols = [
    "locus", "chr", "pos", "genes", "defining_enzyme", "enzyme_raising_allele",
    "z_defining_enzyme", "z_MVP_ALT_EUR", "z_MVP_AST_EUR", "z_MVP_NAFLD_EUR",
    "z_MVP_Cirrhosis_EUR", "z_Pazoki_PDFF", "z_UKBB_ALT", "z_UKBB_AST", "z_UKBB_GGT", "cls",
]
print("\n== AST-defined loci ==")
print(L[L.defining_enzyme == "MVP_AST_EUR"][diag_cols].round(2).to_string(index=False))

print("\n== Named locus diagnostics (ALDH16A1, TNFSF10, PANX1, GCKR, TM6SF2/MAU2, PPP1R3B) ==")
mask = (
    L.genes.str.contains("ALDH16A1|RASIP1|TNFSF10|PANX1|GCKR|MAU2|PPP1R3B", regex=True)
    | ((L.chr == 19) & (L.pos == 19379549))
)
print(L[mask][diag_cols].round(2).to_string(index=False))

# 5. tests
rows = []
for a in ["state_up", "hep_made", "bulk_up"]:
    r = L[L.cls == "release"][a]
    d = L[L.cls == "disease"][a]
    tab = [[int(r.sum()), int((~r).sum())], [int(d.sum()), int((~d).sum())]]
    o, p = fisher_exact(tab)
    rows.append(
        dict(
            annotation=a,
            release_yes=tab[0][0],
            release_n=len(r),
            disease_yes=tab[1][0],
            disease_n=len(d),
            odds_ratio=o,
            p=p,
        )
    )
T = pd.DataFrame(rows).sort_values("p")
m = len(T)
T["p_holm"] = np.minimum(1, np.maximum.accumulate((m - np.arange(m)) * T.p.values))
T.to_csv(f"{out}/class_annotation_tests.tsv", sep="\t", index=False)
print("\n== Class annotation tests ==")
print(T.to_string(index=False))

D = L.groupby("cls")[zc].median().round(2).T
D.to_csv(f"{out}/class_median_z.tsv", sep="\t")
print("\n== Class median z ==")
print(D.to_string())

show = [
    "locus", "genes", "defining_enzyme", "z_defining_enzyme", "z_MVP_ALT_EUR", "z_MVP_AST_EUR",
    "z_MVP_NAFLD_EUR", "z_MVP_Cirrhosis_EUR", "z_Pazoki_PDFF", "z_UKBB_GGT",
    "state_up", "hep_made", "bulk_up",
]
for c in ["release", "opposite", "disease"]:
    print(f"\n== {c} ({int((L.cls == c).sum())})")
    print(L[L.cls == c].sort_values("z_defining_enzyme", ascending=False)[show].round(1).to_string(index=False))

os.system(f"/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/pip freeze > {out}/pip_freeze.txt")
print("DIG1 JOB 73 20261118 DONE")
