import csv, collections, math, sys
from scipy.stats import fisher_exact
FM="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping"
R=FM+"/results/susie_coloc_r2_20260924T053244Z"
lab={r["study_name"]:r["tier_label"] for r in csv.DictReader(open(FM+"/config/gwas_trait_tier.tsv"),delimiter="\t")}
def num(x):
    try:
        v=float(x); return v if v==v else None
    except: return None
rows=collections.Counter(); susie_tested=collections.Counter(); susie_pos=collections.Counter(); untest=collections.Counter()
abf_pos=collections.Counter(); ncs=collections.Counter(); methods=collections.defaultdict(collections.Counter)
studies=collections.defaultdict(set)
gmax=collections.defaultdict(lambda: {"direct_MASLD":None,"liver_enzyme":None})
ncs_tested=collections.Counter()
for r in csv.DictReader(open(R+"/susie_coloc_all_gwas.csv")):
    c=lab.get(r["gwas_name"])
    if c not in ("direct_MASLD","liver_enzyme"): continue
    g=r["ensembl"]
    if not g.startswith("ENSG"): continue
    rows[c]+=1; studies[c].add(r["gwas_name"])
    methods[c][r["method"]]+=1
    s=num(r["PP.H4.susie"]); a=num(r["PP.H4.abf"])
    if a is not None and a>0.5: abf_pos[c]+=1
    if r["method"]=="susie_untestable_insufficient_shared_posterior": untest[c]+=1
    if s is not None:
        susie_tested[c]+=1
        n=num(r["n_cs_pairs"]); 
        if n: ncs_tested[c]+=n
        if s>0.5: susie_pos[c]+=1
        m=gmax[g][c]
        if m is None or s>m: gmax[g][c]=s
for c in ("direct_MASLD","liver_enzyme"):
    print(c,"studies",len(studies[c]),"rows",rows[c],"susie_tested",susie_tested[c],"susie_pos",susie_pos[c],
          "rate",round(susie_pos[c]/susie_tested[c],4) if susie_tested[c] else None,"untestable",untest[c],
          "sum_n_cs_pairs_tested",ncs_tested[c],"abf_pos",abf_pos[c],"abf_rate",round(abf_pos[c]/rows[c],4))
    print("  methods",dict(methods[c]))
d,e="direct_MASLD","liver_enzyme"
t=[[susie_pos[d],susie_tested[d]-susie_pos[d]],[susie_pos[e],susie_tested[e]-susie_pos[e]]]
print("fisher gene-study susie",t,fisher_exact(t))
t=[[abf_pos[d],rows[d]-abf_pos[d]],[abf_pos[e],rows[e]-abf_pos[e]]]
print("fisher gene-study abf",t,fisher_exact(t))
t=[[30,141],[810,3034]]
print("fisher pairs (claim)",t,fisher_exact(t))
dpos={g for g,v in gmax.items() if v[d] is not None and v[d]>0.5}
epos={g for g,v in gmax.items() if v[e] is not None and v[e]>0.5}
anypos=dpos|epos
print("genes direct>0.5",len(dpos),"enzyme>0.5",len(epos),"both",len(dpos&epos),"union",len(anypos))
drv=collections.Counter()
for g in anypos:
    v=gmax[g]; dv=v[d] if v[d] is not None else -1; ev=v[e] if v[e] is not None else -1
    drv["direct" if dv>ev else "enzyme"]+=1
print("driver by class max",dict(drv))
both=sorted(dpos&epos)
print("direct>0.5 genes whose enzyme max is higher:",sum(1 for g in dpos if (gmax[g][e] or -1)>gmax[g][d]))
print("direct>0.5 genes whose direct max is higher or equal:",sum(1 for g in dpos if (gmax[g][e] or -1)<=gmax[g][d]))
# cross-check with tier12 file driver labels
t12={r["ensembl"]:r for r in csv.DictReader(open(R+"/gene_level_coloc_tier12.csv"))}
cc=collections.Counter()
for g,r in t12.items():
    s=num(r["coloc_best_susie_pp4"])
    if s is not None and s>0.5: cc[lab.get(r["susie_driving_gwas"],"")]+=1
print("tier12 file susie>0.5 by susie_driving_gwas class",dict(cc))
sym={}
for r in csv.DictReader(open(R+"/gene_level_coloc_tier12.csv")): pass
print("direct>0.5 genes with enzyme higher:")
for g in sorted(dpos):
    print(" ",g,round(gmax[g][d],4),None if gmax[g][e] is None else round(gmax[g][e],4), lab.get(t12[g]["susie_driving_gwas"],""))
