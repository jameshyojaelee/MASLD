"""Independent blind re-derivation of Fig 5E genetic context from release files."""
import pandas as pd, gzip, hashlib
from pathlib import Path
P = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
OUT = P / "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/check/F5ef/out"
PAIRS = P / "Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-coloc-r2-20260924T053244Z/genetics/context/genetic_lineage_context_primary_pairs.tsv"
MEM = P / "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-stage-corrected-candidate-2026-08-13-v3/hotspot/program_membership_v2.tsv"
REG = P / "figures/main/fig4_singlecell_programs/source_tables/source_diagnosis_correction_20261014/program_registry_v2.tsv"
for f in (PAIRS, MEM, REG):
    print(hashlib.sha256(f.read_bytes()).hexdigest(), f)
pairs = pd.read_csv(PAIRS, sep="\t")
print("rows", len(pairs), "pairs", pairs[["gwas_name","ensembl"]].drop_duplicates().shape[0], "genes", pairs.ensembl.nunique())
print(pairs.groupby("trait_class")[["gwas_name","ensembl"]].apply(lambda d: d.drop_duplicates().shape[0]))
print(pairs.evidence_state.value_counts())
# check state definition vs masses
pairs["max_src"] = pairs[["gse244832_accessible_mass","gse281367_accessible_mass"]].max(axis=1)
pairs["min_src"] = pairs[["gse244832_accessible_mass","gse281367_accessible_mass"]].min(axis=1)
print(pd.crosstab(pairs.evidence_state, pairs.max_src >= 0.5))
print(pd.crosstab(pairs.evidence_state, pairs.min_src >= 0.5))
mem = pd.read_csv(MEM, sep="\t")
mem = mem[mem.mapped_symbol_status != "unmapped_or_ambiguous"]
print("mapped memberships", len(mem))
reg = pd.read_csv(REG, sep="\t", usecols=["cell_type","module","program_uid","module_name","membership_sha256"])
print("registry", len(reg), reg.program_uid.nunique())
chk = mem[["program_uid","membership_sha256"]].drop_duplicates().merge(reg[["program_uid","membership_sha256"]], on="program_uid", suffixes=("_mem","_reg"))
print("sha mismatches", (chk.membership_sha256_mem != chk.membership_sha256_reg).sum(), "programs", len(chk))
# GENCODE v49 symbol for each pair ensembl
g = pd.read_csv(P / "data/gencode_v49_gene_metadata.tsv.gz", sep="\t")
print(g.columns.tolist()[:12])
LIN = {"hepatocytes":"hepatocyte","cholangiocytes":"cholangiocyte","macrophages":"macrophage","fibroblasts":"stellate","tcells":"t_nk"}
mem["lineage"] = mem.cell_type.map(LIN)
open_states = {"replicated_accessible","source_dependent"}
pairs["is_open"] = pairs.evidence_state.isin(open_states)
# symbol-mismatch check: pair gene symbols not in GENCODE v49 symbol space
idcol = [c for c in g.columns if "id" in c.lower()][0]
symcol = [c for c in g.columns if "name" in c.lower() or "symbol" in c.lower()][0]
g["ens"] = g[idcol].astype(str).str.replace(r"\.\d+$","",regex=True)
e2s = dict(zip(g.ens, g[symcol]))
pairs["v49_symbol"] = pairs.ensembl.map(e2s)
mm = pairs[["gene","ensembl","v49_symbol"]].drop_duplicates()
print("pair genes whose v49 symbol differs from pair symbol:", (mm.gene != mm.v49_symbol).sum())
print(mm[mm.gene != mm.v49_symbol].to_string())
res = []
for key, sym in (("pair_symbol","gene"), ("v49_symbol","v49_symbol")):
    j = mem.merge(pairs, left_on=["mapped_symbol","lineage"], right_on=[sym,"lineage"], how="inner")
    for uid, d in mem.groupby("program_uid"):
        jj = j[j.program_uid == uid]
        op = jj[jj.is_open]
        og = sorted(op.mapped_symbol.unique())
        direct_open = sorted(op[op.trait_class=="direct_MASLD"].mapped_symbol.unique())
        direct_any = sorted(jj[jj.trait_class=="direct_MASLD"].mapped_symbol.unique())
        res.append(dict(join=key, program_uid=uid, cell_type=d.cell_type.iloc[0], module=d.module.iloc[0],
            n_coloc_genes=jj.mapped_symbol.nunique(), n_coloc_pairs=jj[["gwas_name","ensembl"]].drop_duplicates().shape[0],
            n_open_genes=len(og), n_open_direct_genes=len(direct_open), open_genes=";".join(og),
            direct_open_genes=";".join(direct_open), direct_any_genes=";".join(direct_any)))
r = pd.DataFrame(res).merge(reg[["program_uid","module_name"]], on="program_uid")
r.to_csv(OUT / "c01_genetic_context_117_blind.tsv", sep="\t", index=False)
for key in ("pair_symbol","v49_symbol"):
    rr = r[r["join"]==key]
    print(key, "programs:", len(rr), ">=1 open gene:", (rr.n_open_genes>0).sum(), ">=1 coloc gene:", (rr.n_coloc_genes>0).sum(),
          "with open direct:", (rr.n_open_direct_genes>0).sum(), "direct genes:", sorted(set(";".join(rr.direct_open_genes).split(";"))-{""}))
rr = r[r["join"]=="pair_symbol"]
print(rr[rr.n_open_genes>0][["cell_type","module","module_name","n_open_genes","n_open_direct_genes","open_genes","direct_any_genes"]].to_string())
print(rr[rr.program_uid.isin(["hotspot_hepatocytes_f05c535ae5bbc0b9","hotspot_hepatocytes_48f39dd4d817a10e"])].to_string())
