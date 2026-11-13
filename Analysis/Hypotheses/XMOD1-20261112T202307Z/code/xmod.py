#!/usr/bin/env python3
"""XMOD1: genes significant in more than one independent modality, with network and pathway context.

Independent modalities (several analyses inside one modality count once):
  G  genetics: SuSiE-COLOC PP.H4 > 0.5 (post eQTL refit), coding fine-mapped PIP >= 0.5, or caQTL coloc
  B  bulk liver RNA: canonical DEG, SteatoSITE stage replication, continuum, EISA, GTEx postmortem fat
  S  snRNA lineage pseudobulk DE (any lineage, BH q < 0.05 within lineage) or within-hepatocyte stage (CONF1 C3)
  P  liver proteomics DEP (PXD051911, BH q < 0.05)
  V  (annotation only, not counted) bulk RNA genotype signature, PNPLA3 dosage / HSD17B13 carrier
Hotspot program membership and NMF top loadings are not gene-level tests; they are used as context only.

Network: STRING v12 (combined score >= 700) among multi-modality genes; edge count vs 10,000 random sets of the
same size drawn from genes testable in >= 2 modalities (seed 20261112).
Pathways: hypergeometric over-representation (Hallmark, Reactome, KEGG, GO BP), BH within each collection,
universe = genes testable in >= the same number of modalities.
Usage: xmod.py <out_dir>
"""
import os, sys, subprocess
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

R = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"; os.chdir(R)
OUT = sys.argv[1]; os.makedirs(OUT, exist_ok=False)
rng = np.random.default_rng(20261112)
SIG = "RNA-seq/results/evidence_passports/candidates/modality-significant-genes-20261112T200941Z/tables/significant_genes_long.tsv.gz"
CAT = "RNA-seq/results/evidence_passports/candidates/gene-catalog-conf1-eqtlfix-20261110T130857Z/catalog_annotated_gene_table_conf1_eqtlfix.tsv.gz"
DEG = "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv"
LIN = "Analysis/Hypotheses/REVIEW-FIX-20261107T020001Z/singlecell/results/binary"
PROT = "Analysis/Multimodal_Program_Projection/results/proteomics/protein_de_adjusted.tsv"
HEPC = "Analysis/Hypotheses/CONF1-20261109T221846Z/out_c3b_ambient/gene_calls.tsv.gz"
MEM = "RNA-seq/results/evidence_passports/releases/gene-catalog-eligibility-20261106/passport_program_membership.parquet"
PSS = "Analysis/Hypotheses/SS-20261025T135237Z/out_s3d-20261026T233706Z/FSS5_all_programs.tsv"
PW = "figures/candidates/continuum-774-panels-20261011T152106Z/molecular_systems_atlas/source_tables/outcome_blind_gene_membership.tsv.gz"
SYS = "figures/candidates/continuum-774-panels-20261011T152106Z/molecular_systems_atlas/overlays/canonical_atlas_coordinates.tsv"
REG = "data/reference_signatures/scenic_hep_regulons.tsv"
SINFO = "data/string_ppi/9606.protein.info.v12.0.txt.gz"; SLINK = "data/string_ppi/9606.protein.links.v12.0.txt.gz"

cat = pd.read_csv(CAT, sep="\t", low_memory=False).set_index("ensembl_id")
sym2ens = cat.symbol.dropna().reset_index().drop_duplicates("symbol", keep=False).set_index("symbol").ensembl_id
def to_ens(x):
    x = pd.Series(x).astype(str)
    return x.str.split(".").str[0].where(x.str.startswith("ENSG")).fillna(x.map(sym2ens))

# ---- significant calls per modality ----
L = pd.read_csv(SIG, sep="\t", low_memory=False)
gmap = {"coloc_susie_eqtl": "G", "coding_finemap": "G", "caqtl_coloc": "G",
        "bulk_deg_canonical": "B", "steatosite_stage_replicated": "B", "continuum_associated": "B",
        "eisa_exon_intron": "B", "gtex_postmortem_fat": "B",
        "snRNA_within_hepatocyte_stage": "S", "liver_protein_DEP": "P",
        "genotype_PNPLA3_G_dosage": "V", "genotype_HSD17B13_TA_carrier": "V"}
L["mod"] = L.analysis.map(gmap)
calls = L.dropna(subset=["mod"])[["ensembl_id", "mod", "analysis", "effect", "q"]]
lin_rows, lin_tested = [], set()
for f in sorted(os.listdir(LIN)):
    if not f.endswith("_de.csv"): continue
    d = pd.read_csv(f"{LIN}/{f}"); d["ensembl_id"] = to_ens(d.gene).values
    lin_tested |= set(d.ensembl_id.dropna())
    s = d[(d.padj < 0.05) & d.ensembl_id.notna()]
    lin_rows.append(pd.DataFrame({"ensembl_id": s.ensembl_id, "mod": "S", "analysis": "snRNA_" + f[:-7] + "_DE",
                                  "effect": s.logFC, "q": s.padj}))
calls = pd.concat([calls] + lin_rows, ignore_index=True)
calls = calls[calls.ensembl_id.astype(str).str.startswith("ENSG")]

# ---- testable universe per modality ----
deg = pd.read_csv(DEG); deg_ids = set(deg.gene.str.split(".").str[0])
gen_t = set(cat.index[cat.gen_susie_state_t12.isin(["positive", "tested_le_0.5"]) |
                      (cat.coding_cs_state != "no_coding_cs_member") |
                      cat.caqtl_region_state.isin(["colocalized accessible region", "colocalized, lead LD not confirmed",
                                                   "different variant (unsupported)", "indeterminate"])])
pr = pd.read_csv(PROT, sep="\t"); pr["ensembl_id"] = to_ens(pr.gene).values
prot_t = set(pr.ensembl_id.dropna())
sc_t = lin_tested | set(pd.read_csv(HEPC, sep="\t").ensembl_id)
tested = pd.DataFrame({"G": pd.Series(1, index=list(gen_t)), "B": pd.Series(1, index=list(deg_ids)),
                       "S": pd.Series(1, index=list(sc_t)), "P": pd.Series(1, index=list(prot_t))}).fillna(0)
tested["n_testable"] = tested.sum(axis=1)

# ---- per-gene modality summary ----
core = calls[calls["mod"] != "V"]
def sgn(x):
    x = x.dropna(); return "" if x.empty else ("up" if (x > 0).all() else "down" if (x < 0).all() else "mixed")
gene = core.groupby("ensembl_id").agg(mods=("mod", lambda m: "".join(sorted(set(m)))),
                                       analyses=("analysis", lambda a: ";".join(sorted(set(a)))))
gene["n_mod"] = gene.mods.str.len()
for m in "BSP":
    gene[f"dir_{m}"] = core[core["mod"] == m].groupby("ensembl_id").effect.apply(sgn)
# bulk direction from the canonical fit (log2FC sign, BH q) regardless of call, for concordance
deg["ensembl_id"] = deg.gene.str.split(".").str[0]
dd = deg.set_index("ensembl_id")
gene["bulk_log2FC"] = dd.logFC.reindex(gene.index); gene["bulk_q"] = dd.padj.reindex(gene.index)
prd = pr.dropna(subset=["ensembl_id"]).drop_duplicates("ensembl_id").set_index("ensembl_id")
gene["protein_log2FC"] = prd.logFC.reindex(gene.index); gene["protein_q"] = prd.padj.reindex(gene.index)
gene["V_genotype"] = calls[calls["mod"] == "V"].groupby("ensembl_id").analysis.agg(lambda a: ";".join(sorted(set(a))))
gene["directions_concordant"] = gene[["dir_B", "dir_S", "dir_P"]].apply(
    lambda r: (lambda v: "" if len(v) < 2 else ("yes" if len(set(v)) == 1 else "no"))([x for x in r if x in ("up", "down")]), axis=1)
keep = ["symbol", "primary_evidence_class", "gen_coloc_best_susie_pp4", "gen_coloc_best_susie_study",
        "coding_cs_top_variant", "coding_cs_top_consequence", "caqtl_gene_link", "disease_liability_candidate",
        "liability_z_anchor", "liability_z_finngen_nafld", "liability_z_pazoki_pdff", "liability_z_cirrhosis_matched",
        "sc_best_lineage", "zon_mrna_call_yakubovsky2026", "state_partition", "steatosite_stage_num_coef"]
gene = gene.join(cat[[k for k in keep if k in cat]], how="left")
gene = gene.join(tested.n_testable, how="left")

# ---- context: Hotspot programs, pathways, systems, regulons ----
mem = pd.read_parquet(MEM)
pss = pd.read_csv(PSS, sep="\t"); pss["feature_id"] = pss.feature_id.str.strip('"')
mem = mem.merge(pss[["feature_id", "cell_type", "module_name", "beta", "q_bh"]], left_on="program_uid", right_on="feature_id", how="left")
mem["lab"] = mem.cell_type + ":" + mem.module_name + np.where(mem.q_bh < 0.05, np.where(mem.beta > 0, " [stage up]", " [stage down]"), "")
gene["hotspot_programs"] = mem.groupby("ensembl_id").lab.agg(lambda x: " | ".join(sorted(set(x.dropna()))))
pw = pd.read_csv(PW, sep="\t", usecols=["gene_id_base", "feature_id", "collection"])
sysmap = pd.read_csv(SYS, sep="\t").set_index("feature_id").community_id
pw["system"] = pw.feature_id.map(sysmap)
gs = pw.dropna(subset=["system"]).groupby(["gene_id_base", "system"]).size().reset_index(name="n")
gene["continuum_systems"] = gs.sort_values("n", ascending=False).groupby("gene_id_base")[["system", "n"]].apply(
    lambda d: ";".join(f"{s}({n})" for s, n in zip(d.system[:3], d.n[:3])))
reg = pd.read_csv(REG, sep="\t"); reg["ensembl_id"] = to_ens(reg.gene).values
gene["scenic_hep_regulons"] = reg.groupby("ensembl_id").program.agg(lambda x: ";".join(sorted(set(x))))

# ---- STRING network among multi-modality genes ----
info = pd.read_csv(SINFO, sep="\t").rename(columns={"#string_protein_id": "sid"})
sid2ens = info.set_index("sid").preferred_name.map(sym2ens).dropna()
links = pd.read_csv(SLINK, sep=" ")
links = links[links.combined_score >= 700]
links["a"] = links.protein1.map(sid2ens); links["b"] = links.protein2.map(sid2ens)
links = links.dropna(subset=["a", "b"]); links = links[links.a < links.b][["a", "b", "combined_score"]].drop_duplicates(["a", "b"])
multi = set(gene.index[gene.n_mod >= 2]); multi3 = set(gene.index[gene.n_mod >= 3])
in_str = set(links.a) | set(links.b)
node_ix = {g: i for i, g in enumerate(sorted(set(links.a) | set(links.b)))}
ea = links.a.map(node_ix).to_numpy(); eb = links.b.map(node_ix).to_numpy()
def n_edges_idx(idx):
    mask = np.zeros(len(node_ix), bool); mask[idx] = True
    return int((mask[ea] & mask[eb]).sum())
pool = np.array([node_ix[g] for g in sorted(set(tested.index[tested.n_testable >= 2]) & in_str)])
net_res = []
for name, S in [("ge2_modalities", multi), ("ge3_modalities", multi3)]:
    S = S & in_str; obs = n_edges_idx([node_ix[g] for g in S])
    null = np.array([n_edges_idx(rng.choice(pool, len(S), replace=False)) for _ in range(10000)])
    net_res.append(dict(set=name, n_genes_in_string=len(S), edges=obs, null_mean=null.mean(),
                        enrichment=obs / null.mean(), p_perm=(1 + (null >= obs).sum()) / (1 + len(null)), n_perm=len(null)))
pd.DataFrame(net_res).to_csv(f"{OUT}/string_network_enrichment.tsv", sep="\t", index=False)
sub = links[links.a.isin(multi) & links.b.isin(multi)].copy()
deg_in = pd.concat([sub.a, sub.b]).value_counts()
gene["string_degree_in_multimodal"] = deg_in.reindex(gene.index).fillna(0).astype(int)
nb = pd.concat([sub.rename(columns={"a": "g", "b": "n"}), sub.rename(columns={"b": "g", "a": "n"})])
nb["nsym"] = nb.n.map(cat.symbol)
gene["string_neighbors_multimodal"] = nb.groupby("g").nsym.agg(lambda x: ";".join(sorted(x.dropna())))
sub.assign(sym_a=sub.a.map(cat.symbol), sym_b=sub.b.map(cat.symbol),
           mods_a=sub.a.map(gene.mods), mods_b=sub.b.map(gene.mods)).to_csv(f"{OUT}/string_edges_multimodal.tsv", sep="\t", index=False)

# ---- pathway over-representation ----
pwg = pw[pw.collection.isin(["hallmark", "reactome", "kegg", "go_bp"])]
enr = []
for name, S, k in [("ge2_modalities", multi, 2), ("ge3_modalities", multi3, 3)]:
    U = set(tested.index[tested.n_testable >= k]); S = S & U
    sets = pwg[pwg.gene_id_base.isin(U)].groupby(["collection", "feature_id"]).gene_id_base.agg(set)
    for (col, fid), members in sets.items():
        if not (10 <= len(members) <= 500): continue
        hits = members & S
        if len(hits) < 3: continue
        p = hypergeom.sf(len(hits) - 1, len(U), len(members), len(S))
        enr.append(dict(set=name, collection=col, pathway=fid, size=len(members), hits=len(hits),
                        expected=len(members) * len(S) / len(U), p=p,
                        genes=";".join(sorted(cat.symbol.reindex(list(hits)).dropna()))))
enr = pd.DataFrame(enr)
enr["q_bh"] = enr.groupby(["set", "collection"]).p.transform(lambda p: multipletests(p, method="fdr_bh")[1])
enr["fold"] = enr.hits / enr.expected
enr.sort_values(["set", "q_bh"]).to_csv(f"{OUT}/pathway_enrichment.tsv", sep="\t", index=False)

gene = gene.sort_values(["n_mod", "string_degree_in_multimodal"], ascending=[False, False])
gene.reset_index().rename(columns={"index": "ensembl_id"}).to_csv(f"{OUT}/multimodal_genes.tsv", sep="\t", index=False)
calls.to_csv(f"{OUT}/calls_used.tsv.gz", sep="\t", index=False)
tested.to_csv(f"{OUT}/testable_universe.tsv.gz", sep="\t")
with open(f"{OUT}/pip_freeze.txt", "w") as fh:
    fh.write(subprocess.run([sys.executable, "-m", "pip", "freeze"], capture_output=True, text=True).stdout)

pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 80)
print(gene.mods.value_counts().to_string())
print(pd.DataFrame(net_res).to_string(index=False))
print(enr[enr.q_bh < 0.05].sort_values("p").groupby("set").head(15)[["set", "collection", "pathway", "size", "hits", "fold", "q_bh"]].to_string(index=False))
