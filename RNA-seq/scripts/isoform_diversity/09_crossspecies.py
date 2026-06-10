#!/usr/bin/env python3
"""
Reproducible cross-species isoform-switch concordance (review P0/P1).
- STRICT 1:1 orthologs as PRIMARY (mutual-unique in tier-H map); many-to-many = sensitivity.
- Per-diet hypergeometric + BH correction; union test as the single primary.
- Honest framing: writes all stats to cross_species_concordance_stats.tsv.
Outputs: cross_species_dtu_concordance.tsv + cross_species_concordance_stats.tsv
"""
import csv, gzip, collections, math
PROJ = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
RES  = f"{PROJ}/RNA-seq/results/isoform_diversity"
strip = lambda x: (x or "").split(".")[0]

def dtu(f, col="gene_screen_padj"):
    tested, hits = set(), set()
    for r in csv.DictReader(open(f), delimiter="\t"):
        g = strip(r["gene_id"]); tested.add(g)
        v = r.get(col)
        if v not in (None,"","NA") and float(v) < 0.05: hits.add(g)
    return tested, hits

DIETS = ["FPC","HFD"]
mouse = {D: dtu(f"{RES}/mouse/dtu_mouse_group_diet-{D}.tsv") for D in DIETS}
m_tested = set().union(*[mouse[D][0] for D in DIETS])
# use the WELL-CALIBRATED human stratum (total-RNA, null sigma 1.39) not the
# collapse-flagged pooled set (review P0) as the human DTU reference
h_tested, h_dtu = dtu(f"{RES}/human/dtu_human_group_total-RNA.tsv")

# ortholog maps from tier-H pairs
pairs = []
with gzip.open(f"{PROJ}/data/external/orthologs/master_ortholog_table.tsv.gz","rt") as fh:
    r = csv.DictReader(fh, delimiter="\t"); cols = r.fieldnames
    hc = [c for c in cols if c.lower().startswith("tier_h")]
    for row in r:
        m,h = strip(row["mouse_ensembl"]), strip(row["human_ensembl"])
        if m and h and any(str(row.get(c,"")).upper() in ("TRUE","T","1","YES") for c in hc):
            pairs.append((m,h))
m2h_all = collections.defaultdict(set); h2m_all = collections.defaultdict(set)
for m,h in pairs: m2h_all[m].add(h); h2m_all[h].add(m)
# strict 1:1 = mutual unique
m2h_11 = {m:next(iter(hs)) for m,hs in m2h_all.items()
          if len(hs)==1 and len(h2m_all[next(iter(hs))])==1}

def hyper_sf(k,M,K,n):
    if M==0 or n==0: return float("nan")
    return sum(math.comb(K,i)*math.comb(M-K,n-i) for i in range(k,min(K,n)+1))/math.comb(M,n)

def run(m2h, label, sym=None):
    # universe: mouse genes tested in mouse AND with a human ortholog tested in human
    def orth(g):
        v = m2h.get(g);
        return ({v} if isinstance(v,str) else set(v)) if v else set()
    Muniv = [g for g in m_tested if any(h in h_tested for h in orth(g))]
    M = len(Muniv)
    K = sum(1 for g in Muniv if any(h in h_dtu for h in orth(g)))
    sw = set().union(*[mouse[D][1] for D in DIETS])
    draw = [g for g in sw if g in set(Muniv)]
    n = len(draw); k = sum(1 for g in draw if any(h in h_dtu for h in orth(g)))
    p = hyper_sf(k,M,K,n); exp = n*K/M if M else float("nan")
    rows=[(label,"union",M,K,n,k,round(exp,2),round(k/exp,2) if exp else float("nan"),round(p,4),"")]
    # per-diet with BH
    pds=[]
    for D in DIETS:
        dd=[g for g in mouse[D][1] if g in set(Muniv)]
        nk=sum(1 for g in dd if any(h in h_dtu for h in orth(g)))
        pp=hyper_sf(nk,M,K,len(dd)); pds.append((D,len(dd),nk,pp))
    pv=[x[3] for x in pds]; order=sorted(range(len(pv)),key=lambda i:pv[i])
    bh=[0]*len(pv); m_=len(pv)
    for rank,i in enumerate(order,1): bh[i]=min(1,pv[i]*m_/rank)
    for j,(D,nd,nk,pp) in enumerate(pds):
        rows.append((label,D,M,K,nd,nk,round(len(draw)*0,2),"",round(pp,4),f"BH={round(bh[j],4)}"))
    return rows, draw, M, K

rows_11, draw11, M11, K11 = run(m2h_11, "strict_1to1")
rows_all, drawall, Mall, Kall = run(m2h_all, "many_to_many")

# concordance table (1:1 primary): genes switching in mouse with a human-DTU 1:1 ortholog
sym={}
with gzip.open(f"{PROJ}/data/external/orthologs/master_ortholog_table.tsv.gz","rt") as fh:
    for row in csv.DictReader(fh, delimiter="\t"): sym[strip(row["mouse_ensembl"])]=row.get("mouse_symbol","")
sw = set().union(*[mouse[D][1] for D in DIETS])
with open(f"{RES}/cross_species_dtu_concordance.tsv","w") as out:
    w=csv.writer(out,delimiter="\t"); w.writerow(["mouse_ensg","mouse_symbol","human_ensg","diets","ortholog_class"])
    for g in sorted(sw):
        h11 = m2h_11.get(g)
        hall = [h for h in m2h_all.get(g,()) if h in h_dtu]
        if h11 and h11 in h_dtu:
            w.writerow([g, sym.get(g,""), h11, ",".join(D for D in DIETS if g in mouse[D][1]), "1to1"])
        elif hall:
            w.writerow([g, sym.get(g,""), ";".join(hall), ",".join(D for D in DIETS if g in mouse[D][1]), "many"])
with open(f"{RES}/cross_species_concordance_stats.tsv","w") as out:
    w=csv.writer(out,delimiter="\t"); w.writerow(["mapping","test","M_universe","K_bg","n_draw","k_obs","expected","fold","hyper_p","note"])
    for r in rows_11+rows_all: w.writerow(r)
print("=== PRIMARY (strict 1:1) ==="); [print("  ",*r) for r in rows_11]
print("=== SENSITIVITY (many-to-many) ==="); [print("  ",*r) for r in rows_all]
print(f"wrote {RES}/cross_species_dtu_concordance.tsv + cross_species_concordance_stats.tsv")
