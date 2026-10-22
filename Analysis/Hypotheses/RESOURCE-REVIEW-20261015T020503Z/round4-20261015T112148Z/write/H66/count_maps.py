"""H66 plan: identifier-mapping coverage and set sizes only. No detection-by-class or class contrast is computed."""
import gzip, re
import pandas as pd
R = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
h = pd.read_csv(f"{R}/Analysis/MASLD_Model_Benchmark/executions/variant_program_substrate-20260830T000000Z/reference/hgnc_complete_set.txt", sep="\t", dtype=str, low_memory=False)
print("HGNC rows", len(h), "status", h.status.value_counts().to_dict(), "max date_modified", h.date_modified.dropna().max())
h = h[h.status == "Approved"]
acc2ens = {}
for e, u in zip(h.ensembl_gene_id, h.uniprot_ids):
    if isinstance(u, str) and isinstance(e, str):
        for a in u.split("|"):
            acc2ens.setdefault(a, set()).add(e)
sym2ens = {s: {e} for s, e in zip(h.symbol, h.ensembl_gene_id) if isinstance(e, str)}
prev = {}
for col in ["prev_symbol", "alias_symbol"]:
    for e, v in zip(h.ensembl_gene_id, h[col]):
        if isinstance(v, str) and isinstance(e, str):
            for s in v.split("|"):
                prev.setdefault((col, s), set()).add(e)
cat = pd.read_parquet(f"{R}/RNA-seq/results/evidence_passports/candidates/gene-catalog-corrected-r8-20261009T142027Z/passport_gene_index.parquet", columns=["ensembl_id", "symbol", "primary_evidence_class"])
print("catalog", len(cat), "unique ens", cat.ensembl_id.nunique(), "unique symbol", cat.symbol.nunique())
print(cat.primary_evidence_class.astype(str).value_counts().to_dict())
cens = set(cat.ensembl_id)
print("catalog ens in HGNC", len(cens & set(h.ensembl_gene_id.dropna())))
for tag, f in [("liver", "liver_protein_quant.txt"), ("plasma", "plasma_protein_quant.txt")]:
    t = pd.read_csv(f"{R}/data/PXD051911/{f}", sep="\t", usecols=["ProteinAccessions", "Genes"], dtype=str)
    n_acc_map = n_sym_only = n_none = 0
    ens_all = set()
    for acc, g in zip(t.ProteinAccessions, t.Genes):
        accs = [re.sub(r"-\d+$", "", a) for a in str(acc).split(";")]
        es = set().union(*[acc2ens.get(a, set()) for a in accs])
        if es:
            n_acc_map += 1
        else:
            syms = [s for s in str(g).split(";") if s and s != "nan"]
            es = set().union(*[sym2ens.get(s, set()) for s in syms])
            if es: n_sym_only += 1
            else: n_none += 1
        ens_all |= es
    print(tag, "groups", len(t), "accession-mapped", n_acc_map, "symbol-only", n_sym_only, "unmapped", n_none, "ensembl genes", len(ens_all), "in catalog", len(ens_all & cens))
cols = gzip.open(f"{R}/data/external/govaere2026_natgenetics/cosmx/GSM9351704_Leuven_1_exprMat_file.csv.gz", "rt").readline().strip().split(",")
tg = [c for c in cols if c not in ("fov", "cell_ID") and not re.match(r"^(SystemControl\d+|Negative\d+|NegativeAdd)$", c)]
multi = [c for c in tg if "/" in c]
single = [c for c in tg if "/" not in c]
cur = sum(c in sym2ens for c in single)
viaprev = [c for c in single if c not in sym2ens and len(prev.get(("prev_symbol", c), set())) == 1]
none = [c for c in single if c not in sym2ens and c not in viaprev]
print("cosmx targets", len(tg), "multi-gene probes", len(multi), "single current symbol", cur, "via unique prev_symbol", len(viaprev), viaprev, "other", len(none), none)
ka = ["ALB","APOB","CYP3A4","PNPLA3","HSD17B13","SERPINA1","TM6SF2","MTARC1","APOA1","TF","HP","FGA","CRP","KRT19","EPCAM","LUM","DCN","COL1A1","PECAM1","CD68","CD3E","MS4A1","JCHAIN"]
print(cat[cat.symbol.isin(ka)].assign(c=lambda d: d.primary_evidence_class.astype(str)).set_index("symbol")["c"].to_dict())
orr = cat[cat.symbol.str.match(r"^OR\d+[A-Z]+\d+[A-Z]?\d*$", na=False)]
print("olfactory receptor genes in catalog", len(orr), orr.primary_evidence_class.astype(str).value_counts().to_dict())
