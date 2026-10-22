"""AN6 / PLAN A8: Gene Catalog sidecar flag for the 233 calls that change status when four RNA-genotype
ancestry PCs enter the bulk disease-vs-control model (H79 T1, M1). Reads H79 outputs and Catalog r9d; edits neither.
No random draws. Usage: python build_h79_pc_changed_sidecar.py <out_dir>"""
import sys, hashlib, subprocess
from pathlib import Path
import numpy as np, pandas as pd
OUT = Path(sys.argv[1]); OUT.mkdir(parents=True, exist_ok=True)
BASE = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
H79 = BASE / "Analysis/Hypotheses/H79-20261015T125125Z/results/run-20261015T132031Z"
CAT = BASE / "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/catalog_r9d-20261016T105550Z/build_out/catalog_r9d_gene_table.parquet"
def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()

g = pd.read_csv(H79 / "D1_gene_table.tsv.gz", sep="\t", float_precision="round_trip")
F1 = pd.read_csv(H79 / "primary_F1.tsv", sep="\t")
t1 = F1[F1.test.str.startswith("T1")].iloc[0]
assert len(g) == 24115 and int(g.adopted_call.sum()) == 1854
ch = g[g.M1_status.isin(["lost", "gained"])].copy()
vc = ch.M1_status.value_counts().to_dict()
assert len(ch) == int(t1.observed_C) == 233 and vc == {"lost": int(t1.lost), "gained": int(t1.gained)}, vc
assert (ch.M1_status == "lost").eq(ch.adopted_call.astype(bool)).all()  # lost = adopted calls; gained = non-calls

ch["ensembl_id"] = ch.gene.str.split(".").str[0]
lfc_ok_a = ch.adopted_logFC.abs() > 0.5; q_ok_a = ch.adopted_padj < 0.05
lfc_ok_m = ch.M1_logFC.abs() > 0.5; q_ok_m = ch.M1_padj < 0.05
ch["change_reason"] = np.select(
    [(ch.M1_status == "lost") & ~lfc_ok_m & q_ok_m, (ch.M1_status == "lost") & lfc_ok_m & ~q_ok_m, (ch.M1_status == "lost") & ~lfc_ok_m & ~q_ok_m,
     (ch.M1_status == "gained") & ~lfc_ok_a & q_ok_a, (ch.M1_status == "gained") & lfc_ok_a & ~q_ok_a, (ch.M1_status == "gained") & ~lfc_ok_a & ~q_ok_a],
    ["lost_abs_log2fc_le_0.5", "lost_padj_only", "lost_both", "gained_abs_log2fc_crossed_0.5", "gained_padj_only", "gained_both"], "other")
ch["sign_flip"] = np.sign(ch.adopted_logFC) != np.sign(ch.M1_logFC)

W = pd.read_parquet(CAT)
keep = ["ensembl_id", "symbol", "primary_evidence_class", "class_display_r9", "r8_rna_call_state", "rna_canonical_call",
        "rna_ancestry_sensitive_h79", "rna_ancestry_change_h79"]
side = ch.merge(W[keep].rename(columns={"symbol": "catalog_symbol"}), on="ensembl_id", how="left", indicator=True)
side["in_catalog_r9d"] = side._merge == "both"
side = side.drop(columns="_merge")
# a catalog flag already set must match the H79 rule column
m = side.in_catalog_r9d & side.rna_ancestry_sensitive_h79.notna()
assert (side.loc[m, "rna_ancestry_sensitive_h79"].astype(bool) == side.loc[m, "D1_rule_M1"].astype(bool)).all()

side["rna_ancestry_pc_changed_h79"] = side.M1_status  # proposed sidecar flag: lost | gained (233 genes)
cols = ["ensembl_id", "gene", "symbol", "catalog_symbol", "in_catalog_r9d", "primary_evidence_class", "class_display_r9", "r8_rna_call_state",
        "rna_ancestry_pc_changed_h79", "change_reason", "sign_flip", "adopted_call", "adopted_logFC", "adopted_padj", "M1_logFC", "M1_padj",
        "ME_status", "ME_logFC", "ME_padj", "M1_PC1_coef", "M1_PC1_t", "n_perm_changed", "n_perm", "D1_rule_M1", "rna_ancestry_sensitive_h79", "AveExpr"]
side = side[cols].sort_values(["rna_ancestry_pc_changed_h79", "change_reason", "symbol"]).reset_index(drop=True)
side.to_csv(OUT / "sidecar_h79_pc_changed_calls.tsv", sep="\t", index=False)

summ = []
def add(k, v, note=""): summ.append({"quantity": k, "value": v, "note": note})
add("changed_calls_T1", len(side)); add("lost", int((side.rna_ancestry_pc_changed_h79 == "lost").sum())); add("gained", int((side.rna_ancestry_pc_changed_h79 == "gained").sum()))
for k, v in side.change_reason.value_counts().items(): add(f"reason_{k}", int(v))
add("sign_flips", int(side.sign_flip.sum()))
add("meets_64_gene_D1_rule_M1", int(side.D1_rule_M1.astype(bool).sum()), "the existing Catalog r9d flag rna_ancestry_sensitive_h79")
add("in_catalog_r9d", int(side.in_catalog_r9d.sum()))
for k, v in side.loc[side.rna_ancestry_pc_changed_h79 == "lost", "primary_evidence_class"].fillna("NA").value_counts().items(): add(f"lost_by_r8_class_{k}", int(v))
for k, v in side.loc[side.rna_ancestry_pc_changed_h79 == "gained", "primary_evidence_class"].fillna("NA").value_counts().items(): add(f"gained_by_r8_class_{k}", int(v))
for k, v in side.ME_status.value_counts().items(): add(f"EUR_only_status_{k}", int(v))
add("n_perm_changed_median", float(side.n_perm_changed.median()), "of 200 within-cell permuted-PC fits")
# Holm resolution check from the producer's F1 table
for _, r in F1.iterrows():
    B = int(r.B); k = int(r.n_ref_ge_obs)
    add(f"{r.test.split()[0]}_B", B); add(f"{r.test.split()[0]}_n_ref_ge_obs", k)
    add(f"{r.test.split()[0]}_p_recomputed", (1 + k) / (B + 1), f"producer {r.p_empirical_one_sided}")
    add(f"{r.test.split()[0]}_p_floor", 1 / (B + 1), f"producer {r.p_floor}")
pv = F1.p_empirical_one_sided.values; o = np.argsort(pv); h = np.empty(2); h[o[0]] = min(1, 2 * pv[o[0]]); h[o[1]] = min(1, max(h[o[0]], pv[o[1]]))
add("holm_T1_recomputed", h[0], f"producer {F1.p_holm.values[0]}"); add("holm_T2_recomputed", h[1], f"producer {F1.p_holm.values[1]}")
fl = F1.p_floor.values; of = np.argsort(fl); hf = np.empty(2); hf[of[0]] = 2 * fl[of[0]]; hf[of[1]] = max(hf[of[0]], fl[of[1]])
add("holm_floor_T1_and_T2", float(max(hf)), "smallest attainable Holm p for both tests with B 200 and 100 (0 reference draws at or above obs)")
pd.DataFrame(summ).to_csv(OUT / "sidecar_summary.tsv", sep="\t", index=False)
pd.DataFrame({"file": [str(H79 / "D1_gene_table.tsv.gz"), str(H79 / "primary_F1.tsv"), str(CAT)],
              "sha256": [sha(H79 / "D1_gene_table.tsv.gz"), sha(H79 / "primary_F1.tsv"), sha(CAT)]}).to_csv(OUT / "inputs_sha256.tsv", sep="\t", index=False)
pd.DataFrame([{"column": "rna_ancestry_pc_changed_h79", "type": "lost | gained | empty",
  "definition": "Bulk call status changed when PC1-PC4 (RNA-called genotypes) were added to the adopted 774-unit disease-vs-control fit (H79 T1, M1). lost = adopted call no longer called; gained = non-call that became a call. Descriptive; no per-gene test; changes no class.",
  "source": str(H79 / "D1_gene_table.tsv.gz")}]).to_csv(OUT / "column_dictionary_addition.tsv", sep="\t", index=False)
(OUT / "pip_freeze.txt").write_text(subprocess.run([sys.executable, "-m", "pip", "freeze"], capture_output=True, text=True).stdout)
print(pd.DataFrame(summ).to_string())
