# H62 prespec: signature-to-platform mapping counts from annotation and gene lists only; no expression values read.
import pandas as pd, gzip, collections
R="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
deg=pd.read_csv(R+"/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv")
print("deg cols",list(deg.columns)); 
sig=deg[(deg.padj<0.05)&(deg.logFC.abs()>0.5)].copy(); sig["dir"]=["up" if x>0 else "down" for x in sig.logFC]
print("sig",sig.dir.value_counts().to_dict(),"symbol NA",sig.symbol.isna().sum(),"symbol startswith ENSG",sig.symbol.astype(str).str.startswith("ENSG").sum(),"dup symbols",sig.symbol.duplicated().sum())
for g in ["IGFBP2","IL32","RPS4Y1"]:
    print(g, deg.loc[deg.symbol==g,["logFC","padj"]].to_dict("records"))
ann=pd.read_csv(R+"/data/public_r4/GSE66676/GPL6244.annot.gz",sep="\t",skiprows=27,dtype=str,quoting=3).dropna(subset=["Gene symbol"])
ann=ann[ann["ID"].str.match(r"^\d+$")]
single=ann[~ann["Gene symbol"].str.contains("///")]
mp=single.groupby("Gene symbol")["ID"].apply(list)
print("annot rows with symbol",len(ann),"single-symbol",len(single),"unique symbols",mp.size,"multi-symbol rows",len(ann)-len(single))
for d in ["up","down"]:
    s=sig[sig.dir==d]; syms=set(s.symbol.dropna())
    hit=[x for x in syms if x in mp.index]; nprobe=collections.Counter(len(mp[x]) for x in hit)
    print(d,"symbols",len(syms),"mapped",len(hit),"probes per mapped symbol",dict(nprobe))
    print(d,"biotype of mapped vs unmapped" if "gene_biotype" in s.columns else "")
cnt=pd.read_csv(R+"/data/public_r4/GSE185051/GSE185051_Gene_count.csv.gz",usecols=["Gene"])
gs=set(cnt.Gene)
for d in ["up","down"]:
    syms=set(sig[sig.dir==d].symbol.dropna()); print("GSE185051",d,len(syms & gs))
