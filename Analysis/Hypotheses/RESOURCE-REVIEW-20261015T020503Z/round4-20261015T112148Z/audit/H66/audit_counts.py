"""H66 audit: independent recount of every n stated in PRESPEC_H66.md.
Identifier coverage, sample structure and AveExpr pools only. No detection-by-class is computed."""
import gzip, re, json
import numpy as np, pandas as pd

R = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
OUT = f"{R}/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/audit/H66"
res = {}

# HGNC
h = pd.read_csv(f"{R}/Analysis/MASLD_Model_Benchmark/executions/variant_program_substrate-20260830T000000Z/reference/hgnc_complete_set.txt", sep="\t", dtype=str, low_memory=False)
res["hgnc_rows"] = len(h); res["hgnc_status"] = h.status.value_counts().to_dict(); res["hgnc_max_date_modified"] = h.date_modified.dropna().max()
h = h[h.status == "Approved"]
hens = h.dropna(subset=["ensembl_gene_id"])
acc2ens = {}
for e, u in zip(hens.ensembl_gene_id, hens.uniprot_ids):
    if isinstance(u, str):
        for a in u.split("|"):
            acc2ens.setdefault(a.strip(), set()).add(e)
cur = {}
for s, e in zip(hens.symbol, hens.ensembl_gene_id):
    cur.setdefault(s, set()).add(e)
prev = {}
for e, v in zip(hens.ensembl_gene_id, hens.prev_symbol):
    if isinstance(v, str):
        for s in v.split("|"):
            prev.setdefault(s.strip(), set()).add(e)

def sym_to_ens(s):
    if s in cur:
        return cur[s], "current"
    p = prev.get(s, set())
    if len(p) == 1:
        return p, "prev_unique"
    if len(p) > 1:
        return set(), "prev_ambiguous"
    return set(), "none"

# Catalog
cat = pd.read_parquet(f"{R}/RNA-seq/results/evidence_passports/candidates/gene-catalog-corrected-r8-20261009T142027Z/passport_gene_index.parquet", columns=["ensembl_id", "symbol", "primary_evidence_class"])
cat["cls"] = cat.primary_evidence_class.astype(str)
res["catalog_n"] = len(cat); res["catalog_unique_ens"] = int(cat.ensembl_id.nunique()); res["catalog_unique_sym"] = int(cat.symbol.nunique())
res["catalog_classes"] = cat.cls.value_counts().to_dict()
cens = set(cat.ensembl_id)
res["catalog_ens_in_hgnc"] = len(cens & set(hens.ensembl_gene_id))
res["catalog_ens_not_in_hgnc"] = len(cens - set(hens.ensembl_gene_id))
sym2cat = dict(zip(cat.symbol, cat.ensembl_id))

# Proteomes: route counts incl. the prespec's full rule (current, then unique prev)
meta = pd.read_csv(f"{R}/data/PXD051911/meta_data.txt", sep="\t", dtype=str)
res["meta_rows"] = len(meta); res["meta_sample_group"] = meta.sample_group.value_counts(dropna=False).to_dict()
for tag, f, mcol in [("liver", "liver_protein_quant.txt", "liver_proteomics_filename"), ("plasma", "plasma_protein_quant.txt", "plasma_proteomics_filename")]:
    t = pd.read_csv(f"{R}/data/PXD051911/{f}", sep="\t", dtype=str)
    samp = list(t.columns[3:])
    r = {"groups": len(t), "first3": list(t.columns[:3]), "sample_cols": len(samp), "unique_sample_cols": len(set(samp))}
    n_acc = n_sym_cur = n_sym_prev = n_none = 0; ens_all = set(); multi = 0
    for acc, g in zip(t.ProteinAccessions, t.Genes):
        accs = [re.sub(r"-\d+$", "", a.strip()) for a in str(acc).split(";")]
        es = set().union(*[acc2ens.get(a, set()) for a in accs])
        if es:
            n_acc += 1
        else:
            syms = [s.strip() for s in str(g).split(";") if s.strip() and s != "nan"]
            es_cur = set().union(*[cur.get(s, set()) for s in syms]) if syms else set()
            if es_cur:
                n_sym_cur += 1; es = es_cur
            else:
                es_prev = set()
                for s in syms:
                    p = prev.get(s, set())
                    if len(p) == 1: es_prev |= p
                if es_prev:
                    n_sym_prev += 1; es = es_prev
                else:
                    n_none += 1
        if len(es) > 1: multi += 1
        ens_all |= es
    r.update(accession=n_acc, symbol_current=n_sym_cur, symbol_prev_unique=n_sym_prev, unmapped=n_none,
             groups_multi_gene=multi, ens_genes=len(ens_all), ens_in_catalog=len(ens_all & cens))
    if tag == "liver":
        mtarc = t[t.ProteinAccessions.str.contains("Q5VT66|Q969Z3", regex=True)][["ProteinAccessions", "Genes"]]
        r["marc_rows"] = mtarc.values.tolist()
        r["Q5VT66_hgnc"] = sorted(acc2ens.get("Q5VT66", set()))
        m = meta[meta[mcol].notna() & (meta[mcol] != "NA")]
        r["meta_files"] = len(m); r["meta_files_in_quant"] = int(m[mcol].isin(samp).sum())
        r["quant_cols_in_meta"] = int(pd.Series(samp).isin(m[mcol]).sum())
        r["participants"] = int(m.patient_name.nunique())
        r["sample_group"] = m.sample_group.value_counts().to_dict()
        r["saf"] = m.saf_diagnosis.value_counts(dropna=False).to_dict()
    else:
        b = meta[meta.plasma_proteomics_filename.notna() & (meta.plasma_proteomics_filename != "NA")]
        fu = meta[meta.plasma_follow_up_proteomics_filename.notna() & (meta.plasma_follow_up_proteomics_filename != "NA")]
        r["baseline_files"] = len(b); r["baseline_in_quant"] = int(b.plasma_proteomics_filename.isin(samp).sum())
        r["baseline_participants"] = int(b.patient_name.nunique())
        r["baseline_sample_group"] = b.sample_group.value_counts().to_dict()
        r["followup_files"] = len(fu); r["followup_in_quant"] = int(fu.plasma_follow_up_proteomics_filename.isin(samp).sum())
        r["followup_participants"] = int(fu.patient_name.nunique())
        r["followup_participants_in_baseline"] = int(fu.patient_name.isin(b.patient_name).sum())
        r["quant_cols_unassigned"] = int((~pd.Series(samp).isin(set(b.plasma_proteomics_filename) | set(fu.plasma_follow_up_proteomics_filename))).sum())
        r["followup_same_row_as_baseline"] = int((fu.plasma_proteomics_filename.notna() & (fu.plasma_proteomics_filename != "NA")).sum())
    res[tag] = r

# CosMx
cols = gzip.open(f"{R}/data/external/govaere2026_natgenetics/cosmx/GSM9351704_Leuven_1_exprMat_file.csv.gz", "rt").readline().strip().split(",")
sysc = [c for c in cols if re.match(r"^SystemControl\d+$", c)]
neg = [c for c in cols if re.match(r"^(Negative\d+|NegativeAdd)$", c)]
tg = [c for c in cols if c not in ("fov", "cell_ID") and c not in sysc and c not in neg]
shared = [c for c in tg if "/" in c]
single = [c for c in tg if "/" not in c]
routes = {c: sym_to_ens(c)[1] for c in single}
res["cosmx"] = {"header_cols": len(cols), "systemcontrol": len(sysc), "negative": len(neg), "negative_names": neg,
                "targets": len(tg), "shared_probes": shared, "n_shared": len(shared),
                "single_routes": pd.Series(routes).value_counts().to_dict(),
                "single_prev": [c for c, v in routes.items() if v == "prev_unique"],
                "single_other": [c for c, v in routes.items() if v not in ("current", "prev_unique")],
                "duplicate_targets": len(tg) - len(set(tg))}

# Known-answer genes and OR regex
ka = ["ALB","APOB","CYP3A4","PNPLA3","HSD17B13","SERPINA1","TM6SF2","MTARC1","APOA1","TF","HP","FGA","CRP","KRT19","EPCAM","LUM","DCN","COL1A1","PECAM1","CD68","CD3E","MS4A1","JCHAIN"]
kad = cat[cat.symbol.isin(ka)].set_index("symbol").cls
res["known_answer_in_catalog"] = int(len(kad)); res["known_answer_classes"] = kad.value_counts().to_dict()
sc_checks = [("KRT19","Cholangiocytes"),("EPCAM","Cholangiocytes"),("LUM","Fibroblasts"),("DCN","Fibroblasts"),("COL1A1","Fibroblasts"),("ALB","Hepatocytes"),("PECAM1","Endothelial_cells"),("CD68","Macrophages"),("CD3E","T_cells"),("MS4A1","B_cells"),("JCHAIN","Plasma_cells")]
res["sc_known_answer_pairs_listed"] = len(sc_checks)
for name, pat in [("prespec_text", r"^OR\d+[A-Z]+\d+"), ("count_maps_script", r"^OR\d+[A-Z]+\d+[A-Z]?\d*$")]:
    o = cat[cat.symbol.str.match(pat, na=False)]
    res[f"OR_{name}"] = {"n": len(o), "classes": o.cls.value_counts().to_dict()}
o1 = set(cat[cat.symbol.str.match(r"^OR\d+[A-Z]+\d+", na=False)].symbol); o2 = set(cat[cat.symbol.str.match(r"^OR\d+[A-Z]+\d+[A-Z]?\d*$", na=False)].symbol)
res["OR_text_minus_script"] = sorted(o1 - o2)[:40]

# AveExpr, biotype, NC-2 pools
gen = pd.read_csv(f"{R}/data/gencode_v49_gene_metadata.tsv.gz", sep="\t")
deg = pd.read_csv(f"{R}/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv")
deg["ens"] = deg.gene.str.split(".").str[0]
res["deg_rows"] = len(deg); res["deg_unique_ens"] = int(deg.ens.nunique())
c2 = cat.merge(gen[["ensembl_base", "gene_biotype"]], left_on="ensembl_id", right_on="ensembl_base", how="left").merge(deg[["ens", "AveExpr"]], left_on="ensembl_id", right_on="ens", how="left")
res["catalog_rows_after_merge"] = len(c2)
a = c2[c2.cls.isin(["genetically_anchored", "established_state_associated"])]
res["secondary_n"] = len(a); res["secondary_aveexpr_missing"] = int(a.AveExpr.isna().sum())
res["secondary_pc"] = a[a.gene_biotype == "protein_coding"].groupby("cls").size().to_dict()
u = c2[(c2.cls == "unresolved") & (c2.gene_biotype == "protein_coding") & c2.AveExpr.notna()]
res["A3_unresolved_pc_with_aveexpr"] = len(u)
qs = np.quantile(a.AveExpr, np.linspace(0, 1, 11))
dec = pd.qcut(a.AveExpr, 10, labels=False)
need = a.assign(dec=dec).groupby("dec").size()
pool = c2[c2.cls.isin(["unresolved", "untested"]) & (c2.gene_biotype == "protein_coding") & c2.AveExpr.notna()]
res["nc2_pool_total"] = len(pool)
pool_in = pool[(pool.AveExpr >= qs[0]) & (pool.AveExpr <= qs[-1])]
res["nc2_pool_inside_range"] = len(pool_in); res["nc2_pool_outside_range"] = len(pool) - len(pool_in)
pdec = pd.cut(pool_in.AveExpr, bins=qs, labels=False, include_lowest=True)
res["nc2_need_per_decile"] = need.to_dict()
res["nc2_pool_per_decile"] = pd.Series(pdec).value_counts().sort_index().to_dict()

# Pseudobulk row identifiers (names only)
pbn = pd.read_csv(f"{R}/Analysis/SingleCell/results_gpu_v2/pseudobulk/Hepatocytes_pseudobulk.csv", usecols=[0]).iloc[:, 0].astype(str)
isens = pbn.str.match(r"^ENSG\d+(\.\d+)?$")
ens_rows = set(pbn[isens].str.split(".").str[0])
sym_rows = pbn[~isens]
rt = {"ens_rows": int(isens.sum()), "symbol_rows": int((~isens).sum())}
mapped = set(ens_rows & cens)
rt["ens_rows_in_catalog"] = len(mapped)
r_cur = r_prev = r_catsym = r_none = 0
for s in sym_rows:
    e, how = sym_to_ens(s)
    e = e & cens
    if e and how == "current": r_cur += 1; mapped |= e
    elif e and how == "prev_unique": r_prev += 1; mapped |= e
    elif s in sym2cat: r_catsym += 1; mapped.add(sym2cat[s])
    else: r_none += 1
rt.update(symbol_rows_current_in_catalog=r_cur, symbol_rows_prev_in_catalog=r_prev, symbol_rows_catalog_symbol_only=r_catsym, symbol_rows_unmapped=r_none, catalog_genes_with_row=len(mapped))
# what the prespec text rule alone gives (symbols via HGNC current then unique prev; ENSG rows ignored)
m_text = set()
for s in sym_rows:
    e, how = sym_to_ens(s)
    m_text |= (e & cens)
rt["catalog_genes_with_row_text_rule"] = len(m_text)
res["pseudobulk_rows"] = rt
json.dump(res, open(f"{OUT}/audit_counts.json", "w"), indent=1, default=str)
print(json.dumps(res, indent=1, default=str))
