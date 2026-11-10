#!/usr/bin/env python3
"""Code review U7: Catalog r9d H-family fields against independent recomputation (deterministic parts only).
H81: decile = AveExpr decile of the primary fit over the 24,115-gene universe (METHODS "Spike-in calibration");
     22,261 = tested genes without an adopted call; count in deciles 1-2; detection-limit lookup per decile.
H79: ancestry-sensitive count and its split by r8 class; missing RNA never coded False.
H87: coding PIP >= 0.5 genes by class (join only; coding annotation not re-derived).
Usage: u7_catalog_fields.py <out_dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd

OUT = Path(sys.argv[1]); OUT.mkdir(parents=True, exist_ok=True)
ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
CAT = ROOT / "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/catalog_r9d-20261016T105550Z/build_out"
b = pd.read_csv(ROOT / "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv")
b["id"] = b.gene.str.split(".").str[0]
b["call"] = (b.padj < 0.05) & (b.logFC.abs() > 0.5)
cols = ["ensembl_id", "primary_evidence_class", "r8_rna_testability_state", "rna_canonical_call", "rna_ancestry_sensitive_h79",
        "rna_ancestry_change_h79", "rna_expression_decile_h81", "rna_detection_limit_log2fc_h81", "rna_below_detection_limit_h81",
        "coding_cs_max_pip", "coding_cs_state", "gen_susie_state_t12", "r8_genetics_testability_state"]
d = pd.read_csv(CAT / "catalog_r9d_gene_table.tsv.gz", sep="\t", usecols=cols, low_memory=False)
mine = pd.read_csv(OUT.parent.parent / "u6/out/u6_independent_classes.tsv.gz", sep="\t")
L = []
def say(s): print(s); L.append(str(s))

# H81 deciles
b = b.sort_values(["AveExpr", "id"]).reset_index(drop=True)
b["dec"] = pd.qcut(np.arange(len(b)), 10, labels=False) + 1
b["dec_qcut_value"] = pd.qcut(b.AveExpr, 10, labels=False) + 1
say(f"bulk {len(b)} calls {int(b.call.sum())} non-called {int((~b.call).sum())}")
nc = b[~b.call]
say(f"non-called in deciles 1-2 (rank deciles): {int((nc.dec <= 2).sum())} = d1 {int((nc.dec == 1).sum())} + d2 {int((nc.dec == 2).sum())}; (value deciles) {int((nc.dec_qcut_value <= 2).sum())}")
say(f"decile AveExpr ranges: " + "; ".join(f"{k}: {g.AveExpr.min():.2f} to {g.AveExpr.max():.2f}" for k, g in b.groupby("dec")))
m = d.merge(b[["id", "dec", "call"]], left_on="ensembl_id", right_on="id", how="left")
has = m.rna_expression_decile_h81.notna()
say(f"Catalog genes with H81 decile {int(has.sum())}; agree with my decile {int((m.loc[has, 'rna_expression_decile_h81'] == m.loc[has, 'dec']).sum())}")
say(f"Catalog H81 decile present on called genes: {int((has & m.call.fillna(False)).sum())}; absent on non-called bulk genes in Catalog: {int((~has & (m.call == False)).sum())}")
say(f"Catalog non-called with decile 1-2: {int((has & (m.rna_expression_decile_h81 <= 2)).sum())}; 4,323 minus this = sidecar count needed {4323 - int((has & (m.rna_expression_decile_h81 <= 2)).sum())}")
say(f"bulk non-called genes without a Catalog record: {int((~nc.id.isin(d.ensembl_id)).sum())} (printed 113); of them in deciles 1-2: {int(((~nc.id.isin(d.ensembl_id)) & (nc.dec <= 2)).sum())} (printed 73)")
say("rna_below_detection_limit_h81 values: " + str(m.rna_below_detection_limit_h81.value_counts(dropna=False).to_dict()))
say("below-limit by decile: " + str(m[has].groupby("rna_expression_decile_h81").rna_below_detection_limit_h81.agg(lambda s: s.value_counts(dropna=False).to_dict()).to_dict()))
say("detection limit by decile: " + str(m[has].groupby("rna_expression_decile_h81").rna_detection_limit_log2fc_h81.agg(lambda s: sorted(s.dropna().unique().tolist())).to_dict()))
say("H81 fields on RNA-untestable Catalog genes: " + str(m[m.r8_rna_testability_state != "testable"].rna_expression_decile_h81.notna().sum()) +
    " of " + str(int((m.r8_rna_testability_state != "testable").sum())) + "; testability states " + str(m.r8_rna_testability_state.value_counts(dropna=False).to_dict()))

# H79
a = d.rna_ancestry_sensitive_h79
say("H79 ancestry_sensitive values: " + str(a.value_counts(dropna=False).to_dict()))
say("H79 sensitive by class x canonical call: " + str(d[a == True].groupby(["primary_evidence_class", "rna_canonical_call"]).size().to_dict()))
say("H79 sensitive field on RNA-untestable genes (should be NA): " + str(d[d.r8_rna_testability_state != "testable"].rna_ancestry_sensitive_h79.value_counts(dropna=False).to_dict()))

# H87
pip = pd.to_numeric(d.coding_cs_max_pip, errors="coerce")
g = d[pip >= 0.5].merge(mine, left_on="ensembl_id", right_on=mine.columns[0], how="left")
say(f"H87 coding PIP>=0.5 genes in Catalog: {len(g)}; by Catalog class {g.primary_evidence_class.value_counts().to_dict()}")
say(f"H87 PIP>=0.1: {int((pip >= 0.1).sum())}; any PIP: {int(pip.notna().sum())}")
(OUT / "u7_catalog_fields_log.txt").write_text("\n".join(L) + "\n")
