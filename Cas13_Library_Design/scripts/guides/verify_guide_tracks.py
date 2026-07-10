#!/usr/bin/env python
"""
verify_guide_tracks.py — manual-QC the Cas13 guides for a handful of genes and draw
a GENOMIC gene-model track per gene (overview only):
  * confirm each guide TRULY targets the gene (target_seq present in a transcript;
    guide_seq == reverse-complement of target_seq),
  * show WHERE each guide lands on the real gene architecture — exons as BOXES
    (CDS tall/blue, UTR short/gray; lncRNA purple), introns as connecting LINES with
    strand chevrons — so exon-vs-intron is explicit,
  * report the SPACING between a gene's guides in mature-mRNA nt,
  * TIGER + Cas13Design score per guide.

Targets are located by SEQUENCE SEARCH (the stored `position` uses sfriedman's
transcript build, offset a few nt from this gffread transcriptome). The gene model
is the UNION of exons across the gene's transcripts, so every guide lands in an exon.

Output: Cas13_Library_Design/figures/guide_qc/overview_all_genes.pdf
        Cas13_Library_Design/data/guide_qc_verification.tsv
"""
from __future__ import annotations
import csv, gzip
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
from matplotlib.lines import Line2D

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
GCSV = ROOT / "Cas13_Library_Design/data/guides/cas13_library_guides_vM38.csv"
TXFA = ROOT / "Cas13_Library_Design/data/guides/cache/nt_screen/transcriptome.fa"
GTF  = Path("/gpfs/commons/home/jameslee/reference_genome/refdata-cellranger-GRCm39-vM38/genes/genes.gtf.gz")
OUTFIG = ROOT / "Cas13_Library_Design/figures/guide_qc"
OUTTSV = ROOT / "Cas13_Library_Design/data/guide_qc_verification.tsv"
GENES = ["Dgat2","Lpl","Lrrk2","Rc3h2","Abcb1a","Adora1","Plp2","Dohh","Snhg15","2010310C07Rik"]
GUIDE_LEN = 23
CDS_C, UTR_C, LNC_C, GUIDE_C = "#3B6FB6", "#BBBBBB", "#8E6FC7", "#C0143C"

def rc(s): return s.translate(str.maketrans("ACGT","TGCA"))[::-1]

def load_guides():
    per = {}
    with open(GCSV) as f:
        for x in csv.DictReader(f):
            if x["gene_symbol_mouse"] in GENES:
                per.setdefault(x["gene_symbol_mouse"], []).append(x)
    return per

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
    """transcript position (1-based, 5'->3') -> genomic coordinate."""
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
    """region from transcript position p on the reference (CDS/5'UTR/3'UTR/lncRNA)."""
    ex=sorted(ref_tx["exon"], key=lambda e:e[0], reverse=(ref_tx["strand"]=="-"))
    if not ref_tx["CDS"]: return "lncRNA"
    # CDS in tx coords
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

def analyse(per):
    allids=set()
    for g in GENES:
        for gu in per.get(g,[]): allids.update(tx_ids_of(gu))
    seqs=load_seqs(allids); mdls=parse_gtf(allids)
    res={}
    for g in GENES:
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
                           "cas13":gu.get("cas13_score","")})
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
    def X(p): return p-gmin                       # gene-relative genomic coord (nt)
    # intron line + strand chevrons across the whole locus
    ax.plot([0,W],[y0,y0], color="#9aa0a6", lw=0.8, zorder=1)
    for xx in np.linspace(0, W, 14)[1:-1]:
        dx = -W*0.006 if strand=="+" else W*0.006
        ax.plot([xx+dx, xx],[y0+0.03, y0], color="#bbb", lw=0.5, zorder=1)
        ax.plot([xx+dx, xx],[y0-0.03, y0], color="#bbb", lw=0.5, zorder=1)
    # exons
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
    # guides (genomic positions)
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
    # mature-mRNA spacing between consecutive guides (by transcript order)
    for i in range(len(d["onref"])-1):
        a,b=d["onref"][i], d["onref"][i+1]
        if a["gpos"] is None or b["gpos"] is None: continue
        x1,x2=X(a["gpos"]), X(b["gpos"])
        ax.annotate("", xy=(x2,yg-0.14), xytext=(x1,yg-0.14),
                    arrowprops=dict(arrowstyle="<->", color="#999", lw=0.5))
        ax.text((x1+x2)/2, yg-0.22, f"{b['ppos']-a['ppos']} nt (mRNA)",
                ha="center", va="top", fontsize=4.2, color="#777")
    cds_txt = "no CDS (lncRNA)" if not cds_g else f"CDS {cds_g[0]:,}-{cds_g[1]:,}"
    arrow = "5'→3'  +strand" if strand=="+" else "3'→5'  −strand"
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
    per=load_guides(); res=analyse(per)
    fig,axs=plt.subplots(len(GENES),1, figsize=(9.0, 2.1*len(GENES)))
    for ax,g in zip(axs,GENES): plot_gene(ax,g,res[g])
    fig.tight_layout(h_pad=2.6)
    fig.savefig(OUTFIG/"overview_all_genes.pdf", bbox_inches="tight"); plt.close(fig)
    with open(OUTTSV,"w") as out:
        out.write("gene\tbiotype\tref_transcript\tguide\tregion_pool\tregion_recomputed\t"
                  "targets_gene\trc_matches_target\ton_reference_tx\tgenomic_pos\ttiger\tcas13design\n")
        for g in GENES:
            d=res[g]
            for x in d["guides"]:
                out.write(f"{g}\t{d['biotype']}\t{d['ref']}\t{x['id']}\t{x['region_data']}\t"
                          f"{x['region_calc']}\t{x['targets']}\t{x['rc_ok']}\t{x['on_ref']}\t"
                          f"{x['gpos']}\t{x['tiger']}\t{x['cas13']}\n")
    for g in GENES:
        d=res[g]
        print(f"{g:16s} {d['biotype']:14s} ref={d['ref']:22s} guides={len(d['guides'])} "
              f"targeting={sum(1 for x in d['guides'] if x['targets'])} spacing(mRNA)={d['spacing']}")
    print(f"\nfigure -> {OUTFIG}/overview_all_genes.pdf")

if __name__=="__main__":
    main()
