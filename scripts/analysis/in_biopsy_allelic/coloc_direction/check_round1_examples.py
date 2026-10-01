"""Re-derive the round-1 reviewer examples and counts from a direction run (read-only)."""
import sys
import numpy as np
import pandas as pd

pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40); pd.set_option("display.max_colwidth", 70)
d = sys.argv[1]
t = pd.read_csv(f"{d}/direction_by_trait.tsv", sep="\t", low_memory=False)
p = pd.read_csv(f"{d}/direction_signal_pairs.tsv", sep="\t", low_memory=False)
g = pd.read_csv(f"{d}/tag_ld_gene_summary.tsv", sep="\t", low_memory=False)

print("== trait rows by release and state")
t["state"] = np.select([t.direction_sign > 0, t.direction_sign < 0], ["increases", "decreases"], "not_directional")
print(t.groupby(["release", "state"]).size().unstack())
print("genes:", t.groupby("release").ensembl.nunique().to_dict(), "union", t.ensembl.nunique())
print("genes with >=1 directional trait row:", t[t.direction_sign != 0].groupby("release").ensembl.nunique().to_dict())

dirn = t[t.direction_sign != 0]
print("\n== directional rows: marginal p (recomputed from columns)")
for rel, x in dirn.groupby("release"):
    print(rel, "n", len(x), "gwas_p>=1e-3", int((x.gwas_p >= 1e-3).sum()), "eqtl_p>=1e-3", int((x.eqtl_p_hit2 >= 1e-3).sum()),
          "both<1e-5", int(((x.gwas_p < 1e-5) & (x.eqtl_p_hit2 < 1e-5)).sum()),
          "excluded span", int((x.hit1_excluded_span.notna() | x.hit2_excluded_span.notna()).sum()),
          "component_opp", int((x.n_pairs_component_blocked_opposing_call > 0).sum()))

print("\n== reviewer examples")
ex = [("AFF4", "ALT"), ("LEAP2", "ALT"), ("UQCRQ", "ALT"), ("PNPLA6", "NAFLD"), ("MAP3K3", "AST"), ("SCGN", "AST"),
      ("EVI2A", "AST"), ("RFC2", "ALT"), ("HLA-DPA1", "AST"), ("HLA-DPB1", "AST"), ("MAU2", "ALT"), ("MAU2", "AST"),
      ("MAU2", "NAFLD"), ("CHID1", "AST"), ("SLC39A8", "AST"), ("C2orf16", "PDFF"), ("LRRC14", "ALT"),
      ("LIPC", "GGT"), ("CYFIP2", "AST"), ("CWF19L1", "ALT"), ("CWF19L1", "AST"), ("HMGN4", "AST"), ("FADS1", "ALT")]
cols = ["release", "reported_study", "gwas_p", "eqtl_p_hit2", "direction_sign", "not_directional_reason"]
flag = [c for c in t.columns if c.startswith("trait_tag_r2_ge_0.8__")]
for gene, tr in ex:
    x = t[(t.gene == gene) & (t.trait == tr)]
    for r in x.itertuples():
        print(f"{gene:9s} {tr:5s} {r.release[:9]:9s} {str(r.reported_study):24s} gp={r.gwas_p:.3g} ep={r.eqtl_p_hit2:.3g} "
              f"sign={int(r.direction_sign):+d} reason={r.not_directional_reason} | "
              + " ".join(f"{c.split('__')[1]}={getattr(r, c) if hasattr(r, c) else x.loc[r.Index, c]}" for c in flag))

print("\n== frequency-block examples (pairs)")
for gene, hit in [("FADS1", "11:61581438"), ("ELOC", "8:74866886"), ("LY96", "8:74866886"), ("MSL2", "3:135910309"),
                  ("SLC33A1", "3:155420125")]:
    x = p[(p.gene == gene) & (p.hit1 == hit)]
    for r in x.itertuples():
        print(gene, r.release[:9], r.study, r.same_lead, int(r.direction_sign), r.all_block_reasons)

print("\n== gene-level tag flags by list")
for c in [c for c in g.columns if c.startswith("gene_tag_r2_ge_0.8__")]:
    print(c.split("__")[1], g.groupby("release")[c].value_counts().unstack(fill_value=0).to_dict("index"))
print("FTSJ3", g[g.gene == "FTSJ3"][[c for c in g.columns if c.startswith("gene_tag_r2_ge_0.8__")]].to_dict("records"))
print(g.ld_reference_eur.unique())

print("\n== directional rows per trait-level tag flag (t3prime_v2_kept)")
c = "trait_tag_r2_ge_0.8__t3prime_v2_kept"
print(t.groupby(["release", c]).apply(lambda x: f"{int((x.direction_sign != 0).sum())}/{len(x)}").unstack(0))
