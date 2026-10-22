# H62 audit: metadata, annotation and gene-list checks only. No GSE66676 or GSE185051 expression value is read or summarised.
import gzip, collections, re, pandas as pd
R="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
deg=pd.read_csv(R+"/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv")
sig=deg[(deg.padj<0.05)&(deg.logFC.abs()>0.5)].copy(); sig["dir"]=["up" if x>0 else "down" for x in sig.logFC]
ens=sig[sig.symbol.astype(str).str.startswith("ENSG")]
print("ENSG-only symbol calls by dir",ens.dir.value_counts().to_dict())
for g in ["IGFBP2","IL32","RPS4Y1"]:
    print(g, deg.loc[deg.symbol==g,["gene","logFC","P.Value","padj"]].to_string(index=False,header=False))
print("min nonzero padj in table", deg.padj[deg.padj>0].min(), "n padj==0", (deg.padj==0).sum())
# GENCODE v49 names for ENSG-only calls
gm=pd.read_csv(R+"/data/gencode_v49_gene_metadata.tsv.gz",sep="\t")
gm["id"]=gm.gene_id.str.replace(r"\..*","",regex=True)
ens=ens.assign(id=ens.gene.str.replace(r"\..*","",regex=True)).merge(gm[["id","gene_name","gene_biotype"]],on="id",how="left")
ens["named"]=~ens.gene_name.fillna("ENSG").str.startswith("ENSG")
print("ENSG-only calls with a non-ENSG GENCODE v49 gene_name, by dir", ens.groupby("dir").named.sum().to_dict())
print("biotype of ENSG-only calls", ens.gene_biotype.value_counts().head(6).to_dict())
# GPL6244 annotation
ann=pd.read_csv(R+"/data/public_r4/GSE66676/GPL6244.annot.gz",sep="\t",skiprows=27,dtype=str,quoting=3)
ann=ann[ann["ID"].str.match(r"^\d+$",na=False)]
single=ann.dropna(subset=["Gene symbol"]); single=single[~single["Gene symbol"].str.contains("///")]
mp=single.groupby("Gene symbol")["ID"].apply(list)
named=ens[ens.named]
print("ENSG-only calls whose GENCODE name is a single-symbol GPL6244 symbol, by dir", named[named.gene_name.isin(mp.index)].groupby("dir").size().to_dict())
for g in ["RPS4Y1","DDX3Y","EIF1AY","UTY","KDM5D","XIST","IGFBP2","IL32","NLGN4Y","USP9Y","ZFY","TXLNGY"]:
    print("GPL6244 single-symbol", g, mp.get(g))
print("rows containing XIST anywhere", ann[ann["Gene symbol"].fillna("").str.contains(r"\bXIST\b")][["ID","Gene symbol"]].values.tolist())
print("rows containing IL32 or IGFBP2 anywhere", ann[ann["Gene symbol"].fillna("").str.contains(r"\b(IL32|IGFBP2)\b")][["ID","Gene symbol"]].values.tolist())
# probe IDs present in the GSE66676 matrix table (first field only)
ids=[]; inside=False
with gzip.open(R+"/data/public_r4/GSE66676/GSE66676_series_matrix.txt.gz","rt") as f:
    for line in f:
        if line.startswith("!series_matrix_table_begin"): inside=True; next(f); continue
        if line.startswith("!series_matrix_table_end"): break
        if inside: ids.append(line.split("\t",1)[0].strip('"'))
        elif line.startswith(("!Sample_data_processing","!Series_overall_design","!Sample_source_name","!Sample_characteristics_ch1\t\"sample type","!Series_summary","!Series_pubmed","!Sample_platform","!Sample_label_ch1","!Sample_extract_protocol")):
            vals=line.rstrip("\n").split("\t"); print(vals[0], sorted(collections.Counter(vals[1:]).items())[:3] if vals[0].startswith("!Sample") else vals[1:])
ids=set(ids); print("matrix probe rows",len(ids))
sigsym={d:set(sig[sig.dir==d].symbol) for d in ["up","down"]}
for d in ["up","down"]:
    hit=[s for s in sigsym[d] if s in mp.index]; inmat=[s for s in hit if any(p in ids for p in mp[s])]
    allp=[p for s in hit for p in mp[s]]
    print(d,"mapped by annotation",len(hit),"with >=1 probe in matrix",len(inmat),"probes",len(allp),"probes in matrix",sum(p in ids for p in allp))
print("single-symbol probes in matrix", sum(p in ids for l in mp for p in l), "of", sum(len(l) for l in mp))
for p in ["8176375","8176624","8176719","8177137","8177232","8048205","7992828"]: print("probe",p,"in matrix",p in ids)
# GSE185051 counts: gene column only
g=pd.read_csv(R+"/data/public_r4/GSE185051/GSE185051_Gene_count.csv.gz",usecols=[0]).iloc[:,0]
print("GSE185051 rows",len(g),"duplicated symbols",sorted(g[g.duplicated(keep=False)].unique().tolist()))
with gzip.open(R+"/data/public_r4/GSE185051/GSE185051_Gene_count.csv.gz","rt") as f: cols=f.readline().strip().split(",")
print("GSE185051 count columns", cols)
with gzip.open(R+"/data/public_r4/GSE185051/GSE185051_series_matrix.txt.gz","rt") as f:
    for line in f:
        if line.startswith(("!Sample_data_processing","!Sample_title","!Series_overall_design","!Sample_source_name")): print(line[:700].rstrip())
# --- part 2: matrix-restricted mapping detail
print("=== part 2")
for d in ["up","down"]:
    hit=[s for s in sigsym[d] if s in mp.index]
    inm={s:[p for p in mp[s] if p in ids] for s in hit}
    inm={s:v for s,v in inm.items() if v}
    print(d,"symbols with probes in matrix",len(inm),"probes per symbol",dict(sorted(collections.Counter(len(v) for v in inm.values()).items())),"probes",sum(len(v) for v in inm.values()))
    print(d,"annotation-mapped symbols with no matrix probe",sorted(set(hit)-set(inm)))
    print(d,"symbols losing some probes",{s:(len(mp[s]),len(inm.get(s,[]))) for s in hit if 0<len(inm.get(s,[]))<len(mp[s])})
print("GSE185051 has IL32, IGFBP2, XIST, RPS4Y1:", [x in set(g) for x in ["IL32","IGFBP2","XIST","RPS4Y1"]])
print("GSE185051 date-like symbols",sorted([x for x in g if re.match(r"^\d+-(Mar|Sep|Dec|Feb|Apr|Mai|Nov|Oct)$",str(x))])[:40])
# --- part 3: NC-1 pool size, multi-probe symbols, GSE66676 sample titles (header only)
print("=== part 3")
allsig=set(sig.symbol)
pool=[s for s in mp.index if s not in allsig and any(p in ids for p in mp[s])]
print("NC-1 pool: single-symbol GPL6244 symbols with >=1 matrix probe, not a signature symbol:", len(pool))
for d in ["up","down"]:
    print(d,"symbols with >2 matrix probes",{s:len([p for p in mp[s] if p in ids]) for s in sigsym[d] if s in mp.index and len([p for p in mp[s] if p in ids])>2})
with gzip.open(R+"/data/public_r4/GSE66676/GSE66676_series_matrix.txt.gz","rt") as f:
    for line in f:
        if line.startswith("!Sample_title"):
            t=[x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]; print("GSE66676 titles",len(t),"unique",len(set(t)),"first",t[:4])
        if line.startswith("!Sample_description"): print("description", sorted(collections.Counter(line.rstrip("\n").split("\t")[1:]).items())[:3])
        if line.startswith("!series_matrix_table_begin"): break
