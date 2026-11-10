"""RV8 post hoc (written after the M2 check failed): does the pair's own GWAS credible set (same study, containing the
anchor) include a coding member of another gene? Not prespecified. Usage: python -I rv8_posthoc_cs_membership.py <rv8_out_dir>"""
import sys
from pathlib import Path

import pandas as pd

OUT = Path(sys.argv[1])
MEM = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/results/candidates/noncoding-dna-precoloc-2026-08-12/run_19783321/credible_set_members_annotated.tsv")
mem = pd.read_csv(MEM, sep="\t", low_memory=False, usecols=["credible_set_uid", "study", "chromosome", "position", "pip_sum_gate_passed"])
mem = mem[mem.pip_sum_gate_passed.astype(str).str.upper().eq("TRUE")]
mem["chr"] = mem.chromosome.astype(str).str.replace("chr", "").astype(int)
cod = pd.read_csv(OUT / "coding_members_mapped.tsv", sep="\t", low_memory=False)
pf = pd.read_csv(OUT / "pair_flags.tsv", sep="\t")
pf = pf[pf.analysis == "primary"].copy()
pf["chr"] = pf.anchor.str.split(":").str[0].astype(int)
pf["pos"] = pf.anchor.str.split(":").str[1].astype(int)
cs = mem.merge(pf[["pair", "gwas", "chr", "pos"]].rename(columns={"gwas": "study"}), left_on=["study", "chr", "position"],
               right_on=["study", "chr", "pos"])[["pair", "credible_set_uid"]].drop_duplicates()
x = cs.merge(cod[["credible_set_uid", "chr", "position", "normalized_susie_pip", "coding_genes", "coding_gene_ids", "consequence_category"]],
             on="credible_set_uid")
x = x.merge(pf[["pair", "gene", "ensembl"]], on="pair")
x["other"] = [str(e) not in str(i).split(";") and str(i) != "nan" for e, i in zip(x.ensembl, x.coding_gene_ids)]
rows = []
for r in pf.itertuples():
    y = x[x.pair == r.pair]
    o = y[y.other].sort_values("normalized_susie_pip", ascending=False)
    rows.append(dict(pair=r.pair, gene=r.gene, gwas=r.gwas, tier=r.tier, anchor=r.anchor, anchor_cs_found=r.pair in set(cs.pair),
                     cs_other_coding=";".join(f"{g}@{c}:{p} pip {q:.3f}" for g, c, p, q in zip(o.coding_genes, o.chr, o.position, o.normalized_susie_pip)),
                     cs_other_coding_max_pip=o.normalized_susie_pip.max() if len(o) else 0.0,
                     cs_own_coding=bool((~y.other).any()), prespec_cls=r.cls))
t = pd.DataFrame(rows)
t.to_csv(OUT / "posthoc_cs_membership_pairs.tsv", sep="\t", index=False)
g = t.sort_values("cs_other_coding_max_pip", ascending=False).drop_duplicates("gene")
g["direct"] = g.gene.isin(t[t.tier == 1].gene)
g.to_csv(OUT / "posthoc_cs_membership_genes.tsv", sep="\t", index=False)
with open(OUT / "posthoc_cs_summary.txt", "w") as fh:
    fh.write(f"Pairs whose anchor was found in a same-study gated credible set: {t.anchor_cs_found.sum()} of {len(t)}\n")
    for thr in (0.1, 0.0):
        s = g[g.cs_other_coding_max_pip > thr] if thr == 0 else g[g.cs_other_coding_max_pip >= thr]
        fh.write(f"Genes whose credible set holds another gene's coding member (PIP {'>' if thr == 0 else '>='} {thr}): {len(s)} of 436 ({int(s.direct.sum())} of 22 direct)\n")
    fh.write("\n" + g[g.cs_other_coding_max_pip > 0][["gene", "direct", "gwas", "anchor", "cs_other_coding", "prespec_cls"]].to_string(index=False) + "\n\nM2 check\n")
    for gn in ["SAMM50", "MAU2", "SPATA31H1", "SUPT7L", "IFI27L2", "SERPINA11", "CNOT1", "SPP1"]:
        y = g[g.gene == gn]
        fh.write(f"  {gn}: " + ("not among the 436" if y.empty else (y.cs_other_coding.iloc[0] or "none") + f" (anchor CS found: {y.anchor_cs_found.iloc[0]})") + "\n")
print("done")
