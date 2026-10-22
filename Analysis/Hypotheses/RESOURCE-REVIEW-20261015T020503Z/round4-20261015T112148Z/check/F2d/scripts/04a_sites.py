"""Write unique gated member sites (GRCh37) for the R re-annotation."""
import pandas as pd
ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
OUT = ROOT + "/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/check/F2d/out"
co = pd.read_csv(ROOT + "/figures/candidates/fig2-coordinate-rerun-20261013T213349Z_22218816/multiplicity/corrected_locus_coordinates.tsv", sep="\t")
m = pd.read_csv(ROOT + "/GWAS/finemapping/results/candidates/noncoding-dna-precoloc-2026-08-12/run_19783321/credible_set_members_annotated.tsv", sep="\t", low_memory=False)
print("member rows", len(m))
g = m[m.credible_set_uid.isin(co.credible_set_uid)]
print("gated member rows", len(g), "gated CS", g.credible_set_uid.nunique())
s = g[["chromosome", "position"]].drop_duplicates().rename(columns={"chromosome": "chrom", "position": "pos"})
print("unique sites", len(s))
s.to_csv(OUT + "/gated_member_sites.tsv", sep="\t", index=False)
