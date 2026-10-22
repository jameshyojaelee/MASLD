#!/usr/bin/env python3
"""Render the missing Figure 6 structural/entry/access panels as honest candidates."""
from __future__ import annotations
import argparse
import csv
import json
import os
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

plt.rcParams.update({"font.size":6,"axes.labelsize":6,"xtick.labelsize":6,"ytick.labelsize":6,
                     "pdf.fonttype":42,"ps.fonttype":42,"font.family":os.environ.get("MASLD_FIG_FONT_FAMILY","DejaVu Sans")})
GRAY="#9E9E9E";BLUE="#1565C0";PINK="#C9265E";TEAL="#00695C"


def canvas(width,height):
    fig,ax=plt.subplots(figsize=(width,height));ax.set(xlim=(0,1),ylim=(0,1));ax.axis("off")
    fig.subplots_adjust(left=.015,right=.985,top=.98,bottom=.02)
    return fig,ax


def box(ax,x,y,w,h,label,color):
    ax.add_patch(Rectangle((x,y),w,h,facecolor="white",edgecolor=color,linewidth=.6))
    ax.text(x+w/2,y+h/2,label,ha="center",va="center",fontsize=6)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--sidecar",type=Path,required=True)
    # Display options; the defaults reproduce the earlier panels.
    p.add_argument("--hide-status-lines",action="store_true",help="omit the candidate/local-review status lines on 6C and 6E")
    p.add_argument("--class-header",default="Class in source bundle",help="6C header of the class column")
    args=p.parse_args();root=args.sidecar;panels=root/"panels";panels.mkdir(exist_ok=False)
    status=json.loads((root/"sidecar_completion.json").read_text())
    with (root/"worked_entries.tsv").open() as f:entries=list(csv.DictReader(f,delimiter="\t"))
    qualification="Provisional inputs; bulk identity and counting pending" if status["provisional_inputs"] else "Candidate; scientific review pending"
    # KEY MESSAGE: gene-level calls and assay/program context remain separate before object-specific routing.
    fig,ax=canvas(4.4,1.75)
    box(ax,.01,.58,.27,.24,"Inherited liver-trait\nregulatory evidence",BLUE)
    box(ax,.01,.18,.27,.24,"Established disease-state\nRNA evidence",PINK)
    box(ax,.36,.34,.25,.32,"Stable gene identity\nSeparate assay calls\nTestability / provenance",GRAY)
    box(ax,.70,.58,.29,.24,"Molecular object\nGene / RNA / DNA element",TEAL)
    box(ax,.70,.18,.29,.24,"Unresolved mechanism\nDiscriminating experiment",GRAY)
    for start,end in [((.28,.70),(.36,.55)),((.28,.30),(.36,.45)),((.61,.52),(.70,.70)),((.845,.58),(.845,.42))]:
        ax.annotate("",xy=end,xytext=start,arrowprops={"arrowstyle":"->","lw":.6,"color":GRAY})
    ax.text(.5,.04,"Protein, chromatin and spatial context retain their assay or program unit",ha="center",fontsize=6)
    fig.savefig(panels/"fig6a_catalog_structure.pdf");plt.close(fig)
    # KEY MESSAGE: fixed worked entries expose separate calls and routes; their present classes are provisional.
    fig,ax=canvas(5.1,1.75)
    headers=["Gene","Genetics","RNA",args.class_header,"Gene route / coding-object route"]
    x=[.01,.15,.28,.41,.70]
    for xx,label in zip(x,headers):ax.text(xx,.91,label,ha="left",fontsize=6)
    class_labels={"genetically_anchored":"Genetically anchored","concordant":"Concordant",
                  "unresolved":"Unresolved","untested":"Untested"}
    legacy={"EXP_ALLELE_AWARE_V1":"Allele-aware editing","EXP_MULTI_LAYER_V1":"Allele edit + gene perturbation",
            "EXP_UNRESOLVED_V1":"Measure decisive missing layer","EXP_DIRECT_ASSAY_V1":"Acquire missing decisive assay"}
    obj={"EXP_PROTEIN_STATE_CONTEXT_V2":"Coding-state perturbation","EXP_MEASURE_MISSING_ASSAY_V2":"Measure missing assay"}
    for i,r in enumerate(sorted(entries,key=lambda r:r["symbol"])):
        y=.73-i*.16
        vals=[r["symbol"],r["genetics_call"],r["rna_call"],class_labels.get(r["primary_evidence_class"],r["primary_evidence_class"]),
              legacy.get(r["legacy_next_experiment_rule_id"],r["legacy_next_experiment_rule_id"])+"\n"+obj.get(r["recommended_experiment_rule_id"],r["recommended_experiment_rule_id"])]
        for j,(xx,value) in enumerate(zip(x,vals)):
            ax.text(xx,y,value,ha="left",va="center",fontsize=6,fontstyle="italic" if j==0 else "normal")
    if not args.hide_status_lines:ax.text(.5,.025,qualification,ha="center",fontsize=6,color="#555555")
    fig.savefig(panels/"fig6c_worked_entries.pdf");plt.close(fig)
    # KEY MESSAGE: entry details, object routes and source/download records refer to the same local candidate.
    fig,ax=canvas(4.1,1.7)
    box(ax,.01,.53,.31,.32,"View worked gene entries\nSeparate evidence calls",BLUE)
    box(ax,.36,.53,.29,.32,"Inspect object\nand uncertainty\nGene and object routes",TEAL)
    box(ax,.70,.53,.29,.32,"Read source references\nAnnotation / assay states",GRAY)
    for a,b in [(.32,.36),(.65,.70)]:ax.annotate("",xy=(b,.69),xytext=(a,.69),arrowprops={"arrowstyle":"->","lw":.6,"color":GRAY})
    box(ax,.17,.12,.66,.22,"Download gene calls + molecular-object records\nSame selected inputs; no combined score or rank",GRAY)
    if not args.hide_status_lines:
        ax.text(.5,.94,"Local review interface; unpublished",ha="center",fontsize=6)
        ax.text(.5,.015,qualification,ha="center",fontsize=6,color="#555555")
    fig.savefig(panels/"fig6e_portal_source_access.pdf");plt.close(fig)
    captions=[{"panel":"6A","file":"fig6a_catalog_structure.pdf","caption":"Gene calls remain distinct from program/assay context; molecular objects have explicit unresolved mechanisms and routes."},
              {"panel":"6C","file":"fig6c_worked_entries.pdf","caption":"The four existing review genes use the source bundle classes. Gene-level v1 routes and molecular-object v2 routes remain separate. "+qualification+"."},
              {"panel":"6E","file":"fig6e_portal_source_access.pdf","caption":"Local candidate access and downloads; no portal publication or final manual acceptance is implied. "+qualification+"."}]
    with (root/"panel_captions.tsv").open("x") as f:
        w=csv.DictWriter(f,fieldnames=["panel","file","caption"],delimiter="\t",lineterminator="\n");w.writeheader();w.writerows(captions)
    print(json.dumps({"candidate_not_adopted":True,"panels":[r["file"] for r in captions],"provisional_inputs":status["provisional_inputs"]}))


if __name__=="__main__":main()
