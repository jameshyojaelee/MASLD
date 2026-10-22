import pandas as pd, json, re
pd.set_option("display.width",250); pd.set_option("display.max_columns",80); pd.set_option("display.max_colwidth",90); pd.set_option("display.max_rows",200)
C="RNA-seq/results/evidence_passports/candidates/gene-catalog-corrected-r8-20261009T142027Z/"
a=pd.read_csv(C+"passport_assay_status.tsv",sep="\t")
print("assay_status",a.shape); print(list(a.columns))
for c in a.columns:
    if a[c].nunique()<=12: print(c,":",a[c].value_counts(dropna=False).to_dict())
show=[c for c in a.columns if c in("input_id","entity_id","entity_label","result_grain","status","call_state","testability_state","estimate","effect_unit","q_value","p_value","n_biological_units","n_technical_units","dataset_id","assay")]
print(a[show].to_string())
print("STALE in limitation:",a.limitation.str.contains("STALE").sum() if "limitation" in a else "n/a")
print(a.limitation.value_counts().to_string() if "limitation" in a else "")
print(a.allowed_wording.value_counts().to_string() if "allowed_wording" in a else "")
pi=pd.read_parquet(C+"passport_program_index.parquet"); pm=pd.read_parquet(C+"passport_program_membership.parquet"); pc=pd.read_parquet(C+"passport_program_context.parquet")
print("\nprogram_index",pi.shape,"unique",pi.iloc[:,1].nunique(),list(pi.columns)); print("membership",pm.shape,list(pm.columns)); print("context",pc.shape,list(pc.columns))
pid=[c for c in pm.columns if "program" in c][0]; print("programs in membership",pm[pid].nunique(), "unique gene-program pairs",len(pm.drop_duplicates([pid]+[c for c in pm.columns if c in("ensembl_id","passport_id")][:1])))
for c in pc.columns:
    if pc[c].nunique()<=14 and c not in("analysis_release_id",): print(c,":",pc[c].value_counts(dropna=False).to_dict())
for c in pi.columns:
    if pi[c].nunique()<=8: print("PI",c,":",pi[c].value_counts(dropna=False).to_dict())
for f in ["passport_summary.parquet","passport_evidence.parquet","passport_program_context.parquet","passport_program_index.parquet","passport_program_membership.parquet"]:
    d=pd.read_parquet(C+"portal_export/"+f); hit=[c for c in d.columns if re.search("myojin",c,re.I)]
    cells=sum(d[c].astype(str).str.contains("yojin",case=False).sum() for c in d.columns if d[c].dtype==object or str(d[c].dtype).startswith("string"))
    print("EXPORT",f,d.shape,"myojin cols",hit,"cells mentioning myojin",cells)
s=pd.read_parquet(C+"portal_export/passport_summary.parquet"); print(list(s.columns))
dc=pd.read_csv(C+"passport_domain_coverage.tsv",sep="\t"); print(dc.to_string())
txt=open(C+"portal_candidate/index.html").read()
print("index.html STALE count",txt.count("STALE"), "34_of_447",txt.count("34_of_447"), "228 genetic_only",txt.count("228 genetic_only"))
for f in ["portal_export/passport_source_graph.json","portal_export/passport_domain_coverage.json","portal_export/passport_dictionary.json"]:
    t=open(C+f).read(); print(f,"STALE",t.count("STALE"),"34_of_447",t.count("34_of_447"),"228 genetic_only",t.count("228 genetic_only"),"myojin",len(re.findall("myojin",t,re.I)))
m=re.search(r'assayGrainRecords":\s*\[(\{.*?\})',txt)
print(m.group(1)[:1500] if m else "no match")
