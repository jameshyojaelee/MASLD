# H62 prespec: metadata-only counts. Reads series-matrix header rows only; no expression values.
import gzip, re, statistics as st, collections
D="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/data/public_r4"
def meta(path):
    rows=collections.defaultdict(list); n=None
    with gzip.open(path,"rt") as f:
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                hdr=next(f).rstrip("\n").split("\t"); break
            if line.startswith("!Sample_"):
                p=line.rstrip("\n").split("\t"); key=p[0]; vals=[v.strip('"') for v in p[1:]]
                if key=="!Sample_characteristics_ch1":
                    k=vals[0].split(":")[0]; rows["ch:"+k]=[v.split(": ",1)[1] if ": " in v else v for v in vals]
                else: rows[key]=vals
    return rows,[h.strip('"') for h in hdr]
r,h=meta(f"{D}/GSE66676/GSE66676_series_matrix.txt.gz")
print("GSE66676 samples",len(r["!Sample_geo_accession"]),"matrix columns",len(h)-1, "order equal", h[1:]==r["!Sample_geo_accession"])
hist=r["ch:histology"]; sex=r["ch:gender"]; age=[float(x) for x in r["ch:age (yrs)"]]; bmi=[float(x) for x in r["ch:bmi"]]
print("characteristic keys",[k for k in r if k.startswith("ch:")])
c=collections.Counter(hist); print("histology",dict(c))
for g in sorted(c):
    ix=[i for i,x in enumerate(hist) if x==g]
    print(g,"n",len(ix),"F/M",sum(sex[i]=="female" for i in ix),"/",sum(sex[i]=="male" for i in ix),
      "age med %.2f range %.2f-%.2f"%(st.median([age[i] for i in ix]),min(age[i] for i in ix),max(age[i] for i in ix)),
      "bmi med %.2f range %.2f-%.2f"%(st.median([bmi[i] for i in ix]),min(bmi[i] for i in ix),max(bmi[i] for i in ix)))
dis=[x!="No NAFLD" for x in hist]
print("binary NoNAFLD",dis.count(False),"NAFLD",dis.count(True),"F",sex.count("female"),"M",sex.count("male"))
print("missing age",sum(x!=x for x in age),"missing bmi",sum(x!=x for x in bmi))
r2,h2=meta(f"{D}/GSE185051/GSE185051_series_matrix.txt.gz")
print("GSE185051 samples",len(r2["!Sample_geo_accession"]),"keys",[k for k in r2 if k.startswith("ch:")])
for k in r2:
    if k.startswith("ch:") and k not in ("ch:age","ch:nas"): print(k,dict(collections.Counter(r2[k])))
print("nas",dict(sorted(collections.Counter(r2["ch:nas"]).items())))
a2=[float(x) for x in r2["ch:age"]]; print("age range",min(a2),max(a2))
t=r2["!Sample_title"]; print("titles",len(t),"unique",len(set(t)))
with gzip.open(f"{D}/GSE185051/GSE185051_Gene_count.csv.gz","rt") as f: cols=f.readline().strip().split(",")[1:]
lb=set("L"+x.split("-")[1] for x in t if x.startswith("LB-")); print("count cols",len(cols),"LB titles mapped to Lxx in counts",len(lb & set(cols)),"titles not LB",[x for x in t if not x.startswith("LB-")])
ch=[i for i,x in enumerate(r2["ch:disease"]) if x=="NAFLD disease"]
for k in ("ch:nas","ch:fibrosis_stage","ch:nash_status","ch:gender"):
    print("children",k,dict(sorted(collections.Counter(r2[k][i] for i in ch).items())))
print("children age range",min(a2[i] for i in ch),max(a2[i] for i in ch))
ct=[i for i,x in enumerate(r2["ch:disease"]) if x!="NAFLD disease"]
print("controls",[(t[i],r2["ch:gender"][i],r2["ch:age"][i],r2["ch:nas"][i],r2["ch:fibrosis_stage"][i]) for i in ct])
