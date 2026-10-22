"""H66 audit part 2: fixed CosMx shared-probe expansion table checked against HGNC and the Catalog;
pseudobulk row mapping under the audit-patched rule; OR genes with a pseudobulk row. Names only."""
import re, json
import pandas as pd

R = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
OUT = f"{R}/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/audit/H66"
h = pd.read_csv(f"{R}/Analysis/MASLD_Model_Benchmark/executions/variant_program_substrate-20260830T000000Z/reference/hgnc_complete_set.txt", sep="\t", dtype=str, low_memory=False)
h = h[(h.status == "Approved")]
hens = h.dropna(subset=["ensembl_gene_id"])
cur = dict(zip(hens.symbol, hens.ensembl_gene_id))
cur_any = set(h.symbol)
prev = {}
for e, v in zip(hens.ensembl_gene_id, hens.prev_symbol):
    if isinstance(v, str):
        for s in v.split("|"):
            prev.setdefault(s.strip(), set()).add(e)
cat = pd.read_parquet(f"{R}/RNA-seq/results/evidence_passports/candidates/gene-catalog-corrected-r8-20261009T142027Z/passport_gene_index.parquet", columns=["ensembl_id", "symbol", "primary_evidence_class"])
cens = set(cat.ensembl_id)
hgnc_ens = set(hens.ensembl_gene_id)
cat_nohgnc_sym = dict(zip(cat[~cat.ensembl_id.isin(hgnc_ens)].symbol, cat[~cat.ensembl_id.isin(hgnc_ens)].ensembl_id))

EXP = {
    "CCL3/L1/L3": ["CCL3", "CCL3L1", "CCL3L3"], "CCL4/L1/L2": ["CCL4", "CCL4L1", "CCL4L2"],
    "CXCL1/2/3": ["CXCL1", "CXCL2", "CXCL3"], "EIF5A/L1": ["EIF5A", "EIF5AL1"],
    "FCGR3A/B": ["FCGR3A", "FCGR3B"], "HBA1/2": ["HBA1", "HBA2"], "HCAR2/3": ["HCAR2", "HCAR3"],
    "HLA-DQB1/2": ["HLA-DQB1", "HLA-DQB2"], "HSPA1A/B": ["HSPA1A", "HSPA1B"], "IFNA1/13": ["IFNA1", "IFNA13"],
    "IFNL2/3": ["IFNL2", "IFNL3"], "KRT6A/B/C": ["KRT6A", "KRT6B", "KRT6C"], "MAP1LC3B/2": ["MAP1LC3B", "MAP1LC3B2"],
    "MZT2A/B": ["MZT2A", "MZT2B"], "PF4/V1": ["PF4", "PF4V1"], "SAA1/2": ["SAA1", "SAA2"], "TNXA/B": ["TNXA", "TNXB"],
    "TPSAB1/B2": ["TPSAB1", "TPSB2"], "XCL1/2": ["XCL1", "XCL2"],
}
rows = []
for p, gs in EXP.items():
    for g in gs:
        e = cur.get(g)
        rows.append(dict(probe=p, gene=g, hgnc_current_with_ens=e is not None, hgnc_current_any=g in cur_any,
                         ensembl=e, in_catalog=(e in cens) if e else (g in set(cat.symbol))))
t = pd.DataFrame(rows)
t.to_csv(f"{OUT}/cosmx_shared_probe_expansion_audit.tsv", sep="\t", index=False)
print(t.to_string())
print("expanded genes", len(t), "in catalog", int(t.in_catalog.sum()))

# pseudobulk rows under the patched rule
pbn = pd.read_csv(f"{R}/Analysis/SingleCell/results_gpu_v2/pseudobulk/Hepatocytes_pseudobulk.csv", usecols=[0]).iloc[:, 0].astype(str)
same = True
for L in ["B_cells", "Cholangiocytes", "Endothelial_cells", "Fibroblasts", "Macrophages", "Plasma_cells", "Resident_NK", "T_cells"]:
    o = pd.read_csv(f"{R}/Analysis/SingleCell/results_gpu_v2/pseudobulk/{L}_pseudobulk.csv", usecols=[0]).iloc[:, 0].astype(str)
    same &= o.tolist() == pbn.tolist()
print("row names identical across 9 lineages", same)
route = {}
for s in pbn:
    m = re.match(r"^(ENSG\d+)(\.\d+)?$", s)
    if m:
        e = m.group(1); r = "ensembl"
    elif s in cur:
        e = cur[s]; r = "hgnc_current"
    elif len(prev.get(s, set())) == 1:
        e = next(iter(prev[s])); r = "hgnc_prev_unique"
    elif s in cat_nohgnc_sym:
        e = cat_nohgnc_sym[s]; r = "catalog_symbol_no_hgnc"
    else:
        e = None; r = "unmapped"
    route[s] = (e, r)
rt = pd.DataFrame([(s, e, r) for s, (e, r) in route.items()], columns=["row", "ens", "route"])
rt["in_catalog"] = rt.ens.isin(cens)
print(rt.groupby(["route", "in_catalog"]).size())
mc = rt[rt.in_catalog]
dup = mc.ens.value_counts()
print("catalog genes with >=1 row", mc.ens.nunique(), "genes with >1 row", int((dup > 1).sum()))
print(mc[mc.ens.isin(dup[dup > 1].index)].sort_values("ens").head(30).to_string())
ORpat = r"^OR\d+[A-Z]+\d+[A-Z]?\d*$"
orr = cat[cat.symbol.str.match(ORpat, na=False)]
print("OR genes", len(orr), "with a pseudobulk row", int(orr.ensembl_id.isin(set(mc.ens)).sum()))
ka = ["KRT19", "EPCAM", "LUM", "DCN", "COL1A1", "ALB", "PECAM1", "CD68", "CD3E", "MS4A1", "JCHAIN"]
kae = cat[cat.symbol.isin(ka)].set_index("symbol").ensembl_id
print("sc known-answer genes with a pseudobulk row", {k: bool(v in set(mc.ens)) for k, v in kae.items()})
rt.to_csv(f"{OUT}/pseudobulk_row_map_audit.tsv.gz", sep="\t", index=False)
