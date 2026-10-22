#!/usr/bin/env python3
"""Independent check of r9d: r9c columns identical; new columns recounted from the H79/H81 files."""
import sys
from pathlib import Path
import pandas as pd
OUT = Path(sys.argv[1]); R4 = OUT.parents[1]
R9C = R4 / "catalog_r9c_candidate-20261015T191932Z/build_out"; HYP = R4.parents[1]
a = pd.read_parquet(R9C / "catalog_r9c_candidate_gene_table.parquet"); b = pd.read_parquet(OUT / "catalog_r9d_gene_table.parquet")
res = []
def ck(n, ok, v=""): res.append((n, bool(ok), str(v)))
ck("rows", len(a) == len(b) == 28233, len(b))
diff = [c for c in a.columns if not a[c].equals(b[c])]
ck("r9c columns identical", not diff, diff)
f = pd.read_csv(HYP / "H79-20261015T125125Z/results/run-20261015T132031Z/D1_flagged_genes.tsv", sep="\t")
s = set(f.loc[f.D1_rule_M1.astype(str).str.upper() == "TRUE", "gene"].str.split(".").str[0])
ck("H79 flagged set", set(b.ensembl_id[b.rna_ancestry_sensitive_h79.fillna(False)]) == s, len(s))
ck("H79 lost/gained", b.rna_ancestry_change_h79.value_counts().drop("", errors="ignore").to_dict() == {"lost": 40, "gained": 24})
ck("H79 flag NA without bulk test", (b.rna_ancestry_sensitive_h79.isna() == b.rna_log2fc.isna()).all())
d = pd.read_csv(HYP / "H81-20261015T125233Z/run-20261015T141242Z/summary/D5_detection_limit_lookup_not_supported_genes.tsv.gz", sep="\t")
d["e"] = d.gene.str.split(".").str[0]
mm = b.merge(d, left_on="ensembl_id", right_on="e")
q = pd.read_csv(OUT / "sidecar_ids_without_catalog_record.tsv", sep="\t", keep_default_na=False, dtype=str)
q0 = pd.read_csv(R9C / "sidecar_ids_without_catalog_record.tsv", sep="\t", keep_default_na=False, dtype=str)
ck("sidecar r9c columns identical", q[q0.columns].equals(q0))
qm = q.assign(e=q.ensembl_id.str.split(".").str[0]).merge(d, on="e")
ck("H81 genes matched (table + sidecar)", len(mm) + len(qm) == len(d) and len(set(mm.e) | set(qm.e)) == len(d), f"{len(mm)}+{len(qm)}/{len(d)}")
ck("sidecar H81 decile equal", (qm.rna_expression_decile_h81.astype(int) == qm.decile).all())
ck("sidecar H81 limit equal", (qm.rna_detection_limit_log2fc_h81 == qm.detection_limit_log2FC.astype(str)).all())
ck("H81 decile equal", (mm.rna_expression_decile_h81.astype(int) == mm.decile).all())
ck("H81 limit equal", (mm.rna_detection_limit_log2fc_h81 == mm.detection_limit_log2FC).all())
nb = int(b.rna_below_detection_limit_h81.sum()); nq = int((q.rna_below_detection_limit_h81 == "True").sum())
ck("H81 below-DL count 4,323 (table + sidecar)", nb + nq == 4323 == int(d.decile.isin([1, 2]).sum()), f"{nb}+{nq}")
ck("H81 no supported call has a limit", not (b.rna_canonical_call.fillna(False) & b.rna_detection_limit_log2fc_h81.notna()).any())
cd = pd.read_csv(OUT / "column_dictionary.tsv", sep="\t", keep_default_na=False)
ck("dictionary covers every column", set(cd[cd.table == "catalog_r9d_gene_table"].column) >= set(b.columns), set(b.columns) - set(cd.column))
ck("no candidate release_state left", not (cd.release_state == "candidate").any())
r = pd.DataFrame(res, columns=["check", "pass", "value"]); r.to_csv(OUT / "check.tsv", sep="\t", index=False); print(r.to_string()); sys.exit(0 if r["pass"].all() else 1)
