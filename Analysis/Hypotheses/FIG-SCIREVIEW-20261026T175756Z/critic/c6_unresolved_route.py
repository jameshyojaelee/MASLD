"""Critic check: of the 13,127 'unresolved' Catalog genes routed to 'measure the single most decisive
missing layer' (and counted in the abstract's 94.5% 'obtain the missing evidence layer'), how many
already have both a finite tier-1/2 SuSiE-COLOC result and a bulk RNA test, i.e. no missing layer of the
two the Catalog calls on?"""
import pandas as pd
t = pd.read_csv("fig6/verify_work/v_gene_table.tsv", sep="\t")
print(len(t), t.primary_evidence_class.value_counts().to_dict())
u = t[t.primary_evidence_class == "unresolved"]
fin = u.best_susie.notna()
print("unresolved:", len(u), "| finite SuSiE:", int(fin.sum()), "| no SuSiE:", int((~fin).sum()))
print("unresolved with finite SuSiE, max PP.H3 > 0.8:", int((fin & (u.max_h3 > 0.8)).sum()))
print("rna_call among unresolved:", u.rna_call.value_counts().to_dict())
both = fin & u.rna_call.eq("indeterminate")
print("unresolved with finite SuSiE AND an RNA test (both layers measured, neither supported):", int(both.sum()))
print("share of 28,233:", round(100 * both.sum() / len(t), 2), "%; 26,678 minus these:", 26678 - int(both.sum()), "=", round(100 * (26678 - both.sum()) / len(t), 1), "%")
