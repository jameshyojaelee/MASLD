#!/usr/bin/env python
"""
spotcheck_100_guides.py — SECOND independent manual-QC pass over the Cas13 guide
table. Now a 500-guide pass (was 100) that COVERS ALL GROUND: stratified so every
structural corner of the library is exercised, drawn from genes NOT in the first
spot check.

Same verification logic as verify_guide_tracks.py (the authoritative checks):
  * guide_seq == reverse-complement of target_seq         (internal consistency)
  * target_seq actually occurs in a real annotated mouse transcript of the claimed
    gene, located by SEQUENCE SEARCH against an independent gffread transcriptome
    (NOT sfriedman's stored `position`)                    (targets the right gene)
  * the guide lands in an EXON of the mature mRNA; region (CDS/UTR/lncRNA) is
    recomputed from the GENCODE vM38 GTF and compared to the pool's region label
  * the guide maps to a genomic coordinate on the gene's locus
  * inter-guide spacing in mature-mRNA nt (should be >= 23, the non-overlap rule)
  * TIGER + Cas13Design scores echoed per guide

Sampling (whole-gene, fixed seed, EXCL = 40 already-checked genes removed):
  1. Original-100 prefix preserved (seed=42 accumulation to >=100) so this run is a
     strict superset of the first 100-guide pass.
  2. Stratified top-up guaranteeing >=1 gene per under-represented category:
     every tier (incl. rare mouse_confirmed/positive_control), both biotypes,
     every region label (incl. 5'UTR|CDS / 3'UTR|CDS), mouse_untestable
     (COLOC-exempt), direction_conflict, ortholog_ambiguous, lncRNA.
  3. Random-fill whole genes until >= 500 guides.
  Whole-gene sampling keeps the spacing check and gene-model tracks meaningful.

Outputs:
  Cas13_Library_Design/data/guide_qc_spotcheck2.tsv               (per-guide verdicts)
  Cas13_Library_Design/figures/guide_qc/spotcheck2_gene_tracks.pdf (paginated tracks)
"""
from __future__ import annotations
import csv, gzip, os, random
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Rectangle, Patch
from matplotlib.lines import Line2D

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
GCSV = ROOT / "Cas13_Library_Design/data/guides/cas13_library_guides_vM38.csv"
TXFA = ROOT / "Cas13_Library_Design/data/guides/cache/nt_screen/transcriptome.fa"
GTF  = Path("/gpfs/commons/home/jameslee/reference_genome/refdata-cellranger-GRCm39-vM38/genes/genes.gtf.gz")
OUTFIG = ROOT / "Cas13_Library_Design/figures/guide_qc"
OUTTSV = ROOT / "Cas13_Library_Design/data/guide_qc_spotcheck2.tsv"

# 40 genes already manually checked in the first pass — excluded here.
EXCL = set("""Plin2 Pnpla3 Dgat2 Tm6sf2 Gckr Lipa Lrrk2 Ptprb Slc39a8 Fads1 Rc3h2 Fads3
Hnf1aos1 Carmn 2010310C07Rik 1700055D18Rik Neat1 Gm12610 Snhg15 Mir155hg
Ccl2 Akr1b10 Abcb1a Filip1l Dyrk2 Naa10 Dohh Gm68527 Aprt Cebpa Plp2 Zfp738
Slc27a1 Gls Adora1 Sel1l3 Piezo1 Atg2a Lpl Cers4""".split())

SEED = 42
BASE_GUIDES = 100          # original spot-check-2 prefix (preserved as a subset)
TARGET_GUIDES = 500        # expanded total (guide-count mode, default)
GENES_PER_PAGE = 8
GUIDE_LEN = 23

# --- gene-count mode (companion 1000-gene coordinate table for the plotgardener
#     coverage figure). When SPOTCHECK_TARGET_GENES>0 the sampler fills to that many
#     GENES (not guides), writes the coordinate TSV to SPOTCHECK_OUT_TSV, and (by
#     default) skips the heavy per-gene matplotlib PDF. The default 500-guide QC
#     behaviour is unchanged.
TARGET_GENES = int(os.environ.get("SPOTCHECK_TARGET_GENES", "0"))
OUT_TSV = Path(os.environ["SPOTCHECK_OUT_TSV"]) if os.environ.get("SPOTCHECK_OUT_TSV") else OUTTSV
DRAW_MPL = os.environ.get("SPOTCHECK_DRAW_MPL", "0" if TARGET_GENES else "1") == "1"

# Stratified coverage: (label, per-guide predicate, n_genes to guarantee).
# n_genes >= 999 means "take all eligible genes" (used for the tiny rare strata).
STRATA = [
    ("tier:mouse_confirmed",  lambda r: r["tier"] == "mouse_confirmed",        999),
    ("tier:positive_control", lambda r: r["tier"] == "positive_control",         8),
    ("tier:coloc",            lambda r: r["tier"] == "coloc",                     8),
    ("mouse_untestable",      lambda r: r["mouse_untestable"] == "True",          6),
    ("direction_conflict",    lambda r: r["direction_conflict"] == "True",        5),
    ("ortholog_ambiguous",    lambda r: r["ortholog_ambiguous"] == "True",        6),
    ("region:5'UTR",          lambda r: r["region"] == "5'UTR",                   4),
    ("region:3'UTR",          lambda r: r["region"] == "3'UTR",                   4),
    ("region:5'UTR|CDS",      lambda r: r["region"] == "5'UTR|CDS",               4),
    ("region:3'UTR|CDS",      lambda r: r["region"] == "3'UTR|CDS",               4),
    ("biotype:lncRNA",        lambda r: r["biotype"] == "lncRNA",                20),
]
CDS_C, UTR_C, LNC_C, GUIDE_C = "#3B6FB6", "#BBBBBB", "#8E6FC7", "#C0143C"

def rc(s): return s.translate(str.maketrans("ACGT","TGCA"))[::-1]

def load_all_guides():
    per = {}
    with open(GCSV) as f:
        for x in csv.DictReader(f):
            per.setdefault(x["gene_symbol_mouse"], []).append(x)
    return per

def base_genes(per):
    """Original spot-check-2 prefix: shuffle non-excluded genes (seed) -> accumulate to >=100."""
    cand = sorted(g for g in per if g not in EXCL)   # sort -> deterministic order pre-shuffle
    rng = random.Random(SEED)
    rng.shuffle(cand)
    chosen, n = [], 0
    for g in cand:
        chosen.append(g); n += len(per[g])
        if n >= BASE_GUIDES:
            break
    return chosen

def sample_genes(per):
    """Stratified whole-gene sample to >=500 guides; original-100 prefix kept as a subset."""
    rng = random.Random(SEED + 1)        # +1 so strata/fill shuffles are independent of the base prefix
    order = base_genes(per)              # 1) original-100 prefix
    chosen = set(order)
    ncum = sum(len(per[g]) for g in order)
    src = {g: "base100" for g in order}

    # 2) stratified top-up: guarantee each category has genes in the sample
    for label, pred, k in STRATA:
        elig = sorted(g for g in per
                      if g not in EXCL and g not in chosen
                      and any(pred(r) for r in per[g]))
        rng.shuffle(elig)
        take = elig if k >= 999 else elig[:k]
        for g in take:
            chosen.add(g); order.append(g); ncum += len(per[g]); src[g] = label

    # 3) random-fill whole genes until the target is met. In gene-count mode
    #    (TARGET_GENES>0) the stop condition is a number of GENES; otherwise it is
    #    a cumulative number of GUIDES (the default 500-guide QC pass).
    rest = sorted(g for g in per if g not in EXCL and g not in chosen)
    rng.shuffle(rest)
    for g in rest:
        if TARGET_GENES > 0:
            if len(chosen) >= TARGET_GENES:
                break
        elif ncum >= TARGET_GUIDES:
            break
        chosen.add(g); order.append(g); ncum += len(per[g]); src[g] = "random_fill"

    return order, ncum, src

def tx_ids_of(g): return [s.split(":")[0] for s in g["tx_id_pos"].split("|") if ":" in s]

def load_seqs(ids):
    want=set(ids); seqs={}; cur=None; buf=[]
    with open(TXFA) as f:
        for line in f:
            if line.startswith(">"):
                if cur in want: seqs[cur]="".join(buf)
                cur=line[1:].strip().split()[0]; buf=[]
            elif cur in want: buf.append(line.strip())
    if cur in want: seqs[cur]="".join(buf)
    return seqs

def parse_gtf(tx_ids):
    want=set(tx_ids); tx={}
    with gzip.open(GTF,"rt") as f:
        for line in f:
            if line.startswith("#"): continue
            c=line.rstrip("\n").split("\t")
            if len(c)<9 or c[2] not in ("exon","CDS"): continue
            i=c[8].find('transcript_id "')
            if i<0: continue
            t=c[8][i+15:c[8].find('"',i+15)]
            if t not in want: continue
            d=tx.setdefault(t,{"exon":[],"CDS":[],"strand":c[6],"chrom":c[0]})
            d[c[2]].append((int(c[3]),int(c[4])))
    return tx

def t2g(tx, tpos):
    ex=sorted(tx["exon"], key=lambda e:e[0], reverse=(tx["strand"]=="-"))
    cum=0
    for s,e in ex:
        L=e-s+1
        if tpos<=cum+L:
            off=tpos-cum-1
            return (s+off) if tx["strand"]=="+" else (e-off)
        cum+=L
    return None

def merge(ivs):
    ivs=sorted(ivs); out=[]
    for s,e in ivs:
        if out and s<=out[-1][1]+1: out[-1]=(out[-1][0],max(out[-1][1],e))
        else: out.append((s,e))
    return out

def union_model(txs):
    ex=[]; cds=[]; strand=None; chrom=None
    for t in txs.values():
        ex+=t["exon"]; cds+=t["CDS"]; strand=strand or t["strand"]; chrom=chrom or t.get("chrom")
    ex=merge(ex); cds=merge(cds)
    cds_g=(cds[0][0],cds[-1][1]) if cds else None
    span=(ex[0][0],ex[-1][1]) if ex else (0,1)
    return {"exons":ex,"cds_g":cds_g,"strand":strand,"chrom":chrom,"span":span}

def tx_region(p, ref_tx):
    ex=sorted(ref_tx["exon"], key=lambda e:e[0], reverse=(ref_tx["strand"]=="-"))
    if not ref_tx["CDS"]: return "lncRNA"
    cmin=min(s for s,_ in ref_tx["CDS"]); cmax=max(e for _,e in ref_tx["CDS"])
    def g2t(G):
        c=0
        for s,e in ex:
            L=e-s+1
            if s<=G<=e: return c+((G-s) if ref_tx["strand"]=="+" else (e-G))+1
            c+=L
    a,b=(g2t(cmin),g2t(cmax)) if ref_tx["strand"]=="+" else (g2t(cmax),g2t(cmin))
    a,b=min(a,b),max(a,b)
    if p+GUIDE_LEN-1 < a: return "5'UTR"
    if p > b: return "3'UTR"
    if p < a <= p+GUIDE_LEN-1 or p <= b < p+GUIDE_LEN-1: return "CDS(boundary)"
    return "CDS"

def analyse(per, genes):
    allids=set()
    for g in genes:
        for gu in per.get(g,[]): allids.update(tx_ids_of(gu))
    seqs=load_seqs(allids); mdls=parse_gtf(allids)
    res={}
    for g in genes:
        gl=per.get(g,[]); bt=gl[0]["biotype"]
        cand=set();  [cand.update(tx_ids_of(gu)) for gu in gl]
        cov=lambda t: sum(1 for gu in gl if gu["target_seq"] and gu["target_seq"] in seqs.get(t,""))
        ref=max(cand, key=lambda t:(cov(t), len(seqs.get(t,""))), default=None)
        refseq=seqs.get(ref,"")
        um=union_model({t:mdls[t] for t in cand if t in mdls})
        guides=[]
        for gu in sorted(gl, key=lambda x:x["guide_id"]):
            tgt=gu["target_seq"]; rc_ok=(rc(gu["guide_seq"])==tgt)
            host,hpos=(ref, refseq.find(tgt)+1) if (tgt and tgt in refseq) else (None,None)
            if host is None:
                for t in tx_ids_of(gu):
                    if tgt and tgt in seqs.get(t,""): host,hpos=t, seqs[t].find(tgt)+1; break
            gpos = t2g(mdls[host], hpos) if (host in mdls and hpos) else None
            ppos = refseq.find(tgt)+1 if (tgt and tgt in refseq) else None
            reg = tx_region(ppos, mdls[ref]) if (ppos and ref in mdls) else gu["region"]
            guides.append({"id":gu["guide_id"],"region_data":gu["region"],"region_calc":reg,
                           "rc_ok":rc_ok,"targets":host is not None,"on_ref":host==ref,
                           "ppos":ppos,"gpos":gpos,"tiger":gu.get("tiger_score",""),
                           "cas13":gu.get("cas13_score",""),
                           "tier":gu.get("tier",""),"is_pos":gu.get("is_positive_control",""),
                           "untestable":gu.get("mouse_untestable",""),
                           "single_iso":gu.get("single_isoform",""),
                           "has_coloc":gu.get("has_coloc",""),
                           "ortho_amb":gu.get("ortholog_ambiguous",""),
                           "dir_conf":gu.get("direction_conflict","")})
        onref=sorted([x for x in guides if x["ppos"]], key=lambda z:z["ppos"])
        spacing=[onref[i+1]["ppos"]-onref[i]["ppos"] for i in range(len(onref)-1)]
        res[g]={"biotype":bt,"ref":ref,"reflen":len(refseq),"um":um,
                "guides":guides,"onref":onref,"spacing":spacing}
    return res

def plot_gene(ax, g, d):
    um=d["um"]; ex=um["exons"]; bt=d["biotype"]
    if not ex:
        ax.text(0.5,0.5,f"{g}: no model", transform=ax.transAxes); ax.axis("off"); return
    gmin,gmax=um["span"]; W=max(gmax-gmin,1); strand=um["strand"]; y0=0.0
    cds_g=um["cds_g"]
    def X(p): return p-gmin
    ax.plot([0,W],[y0,y0], color="#9aa0a6", lw=0.8, zorder=1)
    for xx in np.linspace(0, W, 14)[1:-1]:
        dx = -W*0.006 if strand=="+" else W*0.006
        ax.plot([xx+dx, xx],[y0+0.03, y0], color="#bbb", lw=0.5, zorder=1)
        ax.plot([xx+dx, xx],[y0-0.03, y0], color="#bbb", lw=0.5, zorder=1)
    for s,e in ex:
        xs,xe=X(s),X(e)
        if bt=="lncRNA" or not cds_g:
            ax.add_patch(Rectangle((xs,y0-0.07), xe-xs, 0.14, color=LNC_C, alpha=.6, zorder=3))
        else:
            a,b=cds_g
            cs,ce=max(s,a),min(e,b)
            if s<min(e,a):  ax.add_patch(Rectangle((xs,y0-0.05), X(min(e,a))-xs, 0.10, color=UTR_C, zorder=2))
            if cs<ce:       ax.add_patch(Rectangle((X(cs),y0-0.11), X(ce)-X(cs), 0.22, color=CDS_C, alpha=.65, zorder=3))
            if max(s,b)<e:  ax.add_patch(Rectangle((X(max(s,b)),y0-0.05), xe-X(max(s,b)), 0.10, color=UTR_C, zorder=2))
    yg=0.58
    gl=[x for x in d["guides"] if x["gpos"] is not None]
    for x in gl:
        gx=X(x["gpos"])
        ax.plot([gx,gx],[y0+0.12, yg-0.05], color=GUIDE_C, lw=0.7, zorder=4)
        ax.add_patch(Rectangle((gx-W*0.004, yg-0.05), W*0.008, 0.12, color=GUIDE_C, zorder=5))
        ax.text(gx, yg+0.10, x["id"].split("_")[-1], ha="center", va="bottom",
                fontsize=6, fontweight="bold")
        ax.text(gx, yg+0.27, x["region_calc"], ha="center", va="bottom", fontsize=4.6, color="#555")
        def _f(v):
            try: return f"{float(v):.2f}"
            except (ValueError, TypeError): return "NA"
        ax.text(gx, yg+0.39, f"TIGER:{_f(x['tiger'])} / Cas13Design:{_f(x['cas13'])}",
                ha="center", va="bottom", fontsize=3.8, color=GUIDE_C)
    for i in range(len(d["onref"])-1):
        a,b=d["onref"][i], d["onref"][i+1]
        if a["gpos"] is None or b["gpos"] is None: continue
        x1,x2=X(a["gpos"]), X(b["gpos"])
        ax.annotate("", xy=(x2,yg-0.14), xytext=(x1,yg-0.14),
                    arrowprops=dict(arrowstyle="<->", color="#999", lw=0.5))
        ax.text((x1+x2)/2, yg-0.22, f"{b['ppos']-a['ppos']} nt (mRNA)",
                ha="center", va="top", fontsize=4.2, color="#777")
    cds_txt = "no CDS (lncRNA)" if not cds_g else f"CDS {cds_g[0]:,}-{cds_g[1]:,}"
    arrow = "5'->3'  +strand" if strand=="+" else "3'->5'  -strand"
    ax.set_title(f"{g} ({bt}) · {um['chrom']}:{gmin:,}-{gmax:,} ({arrow}) · genomic span {W:,} nt · "
                 f"{len([x for x in d['guides']])} guides", fontsize=8, loc="left")
    ax.set_xlim(-W*0.03, W*1.03); ax.set_ylim(-0.95, 1.10)
    ax.set_yticks([]); ax.set_xlabel("genomic position in locus (nt)", fontsize=7, labelpad=2)
    for s in ("top","left","right"): ax.spines[s].set_visible(False)
    ax.tick_params(labelsize=6)
    leg=[Patch(facecolor=CDS_C, alpha=.65, label="CDS exon"),
         Patch(facecolor=UTR_C, label="UTR exon"),
         Patch(facecolor=LNC_C, alpha=.6, label="lncRNA exon"),
         Line2D([0],[0], color="#9aa0a6", lw=0.9, label="intron"),
         Patch(facecolor=GUIDE_C, label="guide (23 nt)")]
    ax.legend(handles=leg, loc="upper center", bbox_to_anchor=(0.5,-0.30), ncol=5,
              fontsize=5.5, frameon=False, handlelength=1.1, columnspacing=1.0, handletextpad=0.4)

def main():
    OUTFIG.mkdir(parents=True, exist_ok=True)
    per = load_all_guides()
    genes, n_guides, src = sample_genes(per)
    from collections import Counter
    _tgt = f">={TARGET_GENES} genes" if TARGET_GENES else f">={TARGET_GUIDES} guides"
    print(f"[sample] seed={SEED}  genes={len(genes)}  guides={n_guides}  (target {_tgt})")
    print(f"[sample] gene source: {dict(Counter(src.values()))}")
    print(f"[sample] TSV -> {OUT_TSV}  (matplotlib PDF: {'on' if DRAW_MPL else 'skipped'})")
    res = analyse(per, genes)

    # --- per-guide TSV (now self-documenting for coverage) ---
    with open(OUT_TSV,"w") as out:
        out.write("gene\tbiotype\ttier\tref_transcript\tguide\tregion_pool\tregion_recomputed\t"
                  "targets_gene\trc_matches_target\ton_reference_tx\tgenomic_pos\tmrna_pos\ttiger\tcas13design\t"
                  "is_positive_control\tmouse_untestable\thas_coloc\tortholog_ambiguous\t"
                  "direction_conflict\tsingle_isoform\tgene_source\n")
        for g in genes:
            d=res[g]
            for x in d["guides"]:
                out.write(f"{g}\t{d['biotype']}\t{x['tier']}\t{d['ref']}\t{x['id']}\t{x['region_data']}\t"
                          f"{x['region_calc']}\t{x['targets']}\t{x['rc_ok']}\t{x['on_ref']}\t"
                          f"{x['gpos']}\t{x['ppos']}\t{x['tiger']}\t{x['cas13']}\t{x['is_pos']}\t{x['untestable']}\t"
                          f"{x['has_coloc']}\t{x['ortho_amb']}\t{x['dir_conf']}\t{x['single_iso']}\t{src[g]}\n")

    # --- aggregate pass/fail summary ---
    allg=[x for g in genes for x in res[g]["guides"]]
    n=len(allg)
    rc_ok=sum(x["rc_ok"] for x in allg)
    targ=sum(x["targets"] for x in allg)
    onref=sum(x["on_ref"] for x in allg)
    reg_match=sum(1 for x in allg
                  if x["region_calc"]==x["region_data"]
                  or x["region_calc"].startswith("CDS") and x["region_data"]=="CDS")
    spacings=[s for g in genes for s in res[g]["spacing"]]
    min_space=min(spacings) if spacings else None
    n_overlap=sum(1 for s in spacings if s<GUIDE_LEN)
    off_ref=[x for x in allg if x["targets"] and not x["on_ref"]]
    print("\n================ SPOT-CHECK #2 SUMMARY (500-guide expansion) ================")
    print(f"genes / guides checked    : {len(genes)} genes / {n} guides")
    print(f"RC(guide)==target_seq     : {rc_ok}/{n}")
    print(f"target_seq found in a tx  : {targ}/{n}")
    print(f"  ...on the ref transcript: {onref}/{n}  (off-ref = sister isoform of same gene)")
    print(f"region pool==recomputed   : {reg_match}/{n}  (CDS(boundary) counted as CDS)")
    print(f"min inter-guide spacing   : {min_space} nt   overlaps(<{GUIDE_LEN}nt): {n_overlap}")
    if targ < n:
        print(f"  !! {n-targ} guide(s) did NOT find target_seq in any transcript -- INSPECT")
    print(f"off-ref biotypes          : {dict(Counter(x['region_data'] for x in off_ref))} "
          f"(expect lncRNA-dominated)")
    print("============================================================================")

    # --- coverage matrix: prove all ground is covered ---
    def cov(name, key):
        c=Counter(x[key] for x in allg)
        print(f"  {name:24s}: {dict(c)}")
    print("\nCOVERAGE (guides per category in the 500-set):")
    print("  biotype                 :", dict(Counter(res[g]['biotype'] for g in genes for _ in res[g]['guides'])))
    cov("tier",          "tier")
    cov("region(pool)",  "region_data")
    cov("is_positive_control","is_pos")
    cov("mouse_untestable",  "untestable")
    cov("has_coloc",         "has_coloc")
    cov("ortholog_ambiguous","ortho_amb")
    cov("direction_conflict","dir_conf")
    cov("single_isoform",    "single_iso")

    print("\nper-gene:")
    for g in genes:
        d=res[g]
        print(f"  {g:18s} {d['biotype']:14s} src={src[g]:18s} ref={str(d['ref']):22s} "
              f"guides={len(d['guides'])} targeting={sum(1 for x in d['guides'] if x['targets'])} "
              f"spacing(mRNA)={d['spacing']}")

    # --- paginated gene tracks (matplotlib QC pass; skipped in gene-count mode,
    #     where the plotgardener coverage figure renders the 1000-gene tracks) ---
    if DRAW_MPL:
        pdf_path=OUTFIG/"spotcheck2_gene_tracks.pdf"
        with PdfPages(pdf_path) as pdf:
            for i in range(0, len(genes), GENES_PER_PAGE):
                chunk=genes[i:i+GENES_PER_PAGE]
                fig,axs=plt.subplots(len(chunk),1, figsize=(9.0, 2.1*len(chunk)), squeeze=False)
                for ax,g in zip(axs[:,0], chunk): plot_gene(ax,g,res[g])
                fig.tight_layout(h_pad=2.6)
                pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)
        print(f"figure -> {pdf_path}")
    print(f"\nTSV    -> {OUT_TSV}")

if __name__=="__main__":
    main()
