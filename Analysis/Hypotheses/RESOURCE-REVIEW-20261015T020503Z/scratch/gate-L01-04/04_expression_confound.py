"""Gate L01-04: are COLOC genes more highly expressed in hepatocytes than their window neighbours?
Expression = raw-UMI hepatocyte pseudobulk (39 donors) summed / 46,286 nuclei. Windows: TSS within
500 kb of the locus best peak centre; loci = the 34 testable loci from job 3051."""
import numpy as np, pandas as pd
from scipy import stats
ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
S = f"{ROOT}/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/gate-L01-04"
J = f"{ROOT}/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/jobs/3051"
D = f"{ROOT}/Analysis/MASLD_Model_Benchmark/executions/gse296875-crossmodal-recurrence-20260830T193000Z/donor"
rg = pd.read_csv(f"{D}/rna_genes.tsv", sep="\t")
x = np.load(f"{D}/donor_rna_hepatocyte.npy")
print("rna_genes rows", len(rg), "matrix", x.shape, "columns", list(rg.columns))
mean_umi = x.sum(0) / 46286 if x.shape[1] == len(rg) else x.sum(1) / 46286
genes = pd.read_csv(f"{S}/job_gene_table.tsv", sep="\t")
assert len(genes) == len(mean_umi)
genes["mean_umi"] = mean_umi
t = pd.read_csv(f"{J}/locus_testability.tsv", sep="\t")
t = t[t.peak_passes_5pct & (t.coloc_genes_expressed_5pct > 0) & (t.noncoloc_expressed_5pct > 0)]
inp = pd.read_csv(f"{S}/job_locus_inputs.tsv", sep="\t").set_index("locus")
rows = []
for _, L in t.iterrows():
    chrom, s, e = L.peak_id.split("-"); c = (int(s) + int(e)) // 2
    w = genes[(genes.chrom == chrom) & genes.tss.notna() & ((genes.tss - c).abs() <= 500_000) & (genes.mean_umi >= 0.05)]
    cs = set(inp.loc[L.locus, "coloc_ensembl"].split(";"))
    for _, g in w.iterrows():
        rows.append((L.locus, g.ensembl, int(g.ensembl in cs), g.mean_umi))
r = pd.DataFrame(rows, columns=["locus", "ensembl", "coloc", "mean_umi"])
a, b = r[r.coloc == 1].mean_umi, r[r.coloc == 0].mean_umi
print(f"genes with mean UMI/nucleus >= 0.05 in the 34 windows: COLOC {len(a)}, other {len(b)}")
print(f"median mean UMI per hepatocyte nucleus: COLOC {a.median():.3f}, other {b.median():.3f}; ratio {a.median()/b.median():.2f}")
# within-locus comparison: COLOC gene's rank of expression within its window
r["pct"] = r.groupby("locus").mean_umi.rank(pct=True)
print("mean within-locus expression percentile of COLOC genes:", round(r.loc[r.coloc == 1, "pct"].mean(), 3), "(0.5 expected if unrelated)")
print("Mann-Whitney (COLOC > other), descriptive:", stats.mannwhitneyu(a, b, alternative="greater"))
