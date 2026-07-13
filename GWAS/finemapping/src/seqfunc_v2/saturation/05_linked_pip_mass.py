#!/usr/bin/env python3
"""Calculate source-separated ABC and SCENIC+ linked uniform-PIP masses."""
from __future__ import annotations

import argparse,csv,gzip
from collections import defaultdict
from pathlib import Path

from common import OUT,ROOT,atomic_json

def opener(path): return gzip.open(path,"rt",newline="") if str(path).endswith(".gz") else path.open(newline="")

def overlap_links(anchors, path, build, source, chrom_col, start_col, end_col, gene_col, detail_col=None, score_col=None, detail_prefix=""):
    query=defaultdict(list)
    for a in anchors:
        if build=="GRCh37":
            bits=a["representative_hg19"].split(":"); chrom="chr"+bits[0].removeprefix("chr"); pos=int(bits[1])
        else: chrom=a["chr_hg38"];pos=int(a["pos_hg38"])
        query[chrom].append((pos,a))
    links=[]
    with opener(path) as h:
        for row in csv.DictReader(h,delimiter="," if str(path).endswith(".csv") else "\t"):
            chrom=row[chrom_col]
            if chrom not in query: continue
            start=int(float(row[start_col]))+1; end=int(float(row[end_col]))
            for pos,a in query[chrom]:
                if start<=pos<=end:
                    detail=row.get(detail_col,"") if detail_col else ""
                    links.append({"source":source,"source_detail":detail_prefix+detail,
                      "source_score":row.get(score_col,"") if score_col else "","coordinate_build":build,
                      "anchor_rank":a["anchor_rank"],"variant_id_hg38":a["variant_id_hg38"],"target_gene":row[gene_col],
                      "uniform_pip":a["uniform_max_pip"]})
    return links

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--anchors",type=Path,default=OUT/"anchors/anchor_manifest.tsv")
    ap.add_argument("--out-dir",type=Path,default=OUT/"integration");args=ap.parse_args();args.out_dir.mkdir(parents=True,exist_ok=True)
    with args.anchors.open() as h: anchors=list(csv.DictReader(h,delimiter="\t"))
    support_path=OUT/"anchors/anchor_pip_support.tsv"
    if not support_path.is_file(): raise SystemExit(f"missing locus-resolved PIP support: {support_path}")
    with support_path.open() as h: support=list(csv.DictReader(h,delimiter="\t"))
    selected={a["variant_id_hg38"] for a in anchors}
    if not support or any(x["variant_id_hg38"] not in selected for x in support):
        raise SystemExit("anchor PIP support is empty or contains a nonselected anchor")
    abc=ROOT/"data/external/abc_liver/abc_liver_hepg2_enhancers.tsv.gz"
    if not abc.is_file(): raise SystemExit(f"missing source ABC substrate: {abc}")
    links=overlap_links(anchors,abc,"GRCh37","ABC","chr","start","end","TargetGene","CellType","ABC.Score")
    scenic_dir=ROOT/"Analysis/ATAC/Human_Multiome/scenic_plus"
    scenic_files=sorted(scenic_dir.glob("*enhancer_gene_links.csv"))
    if not scenic_files: raise SystemExit(f"no SCENIC+ enhancer-gene files in {scenic_dir}")
    for path in scenic_files:
        links += overlap_links(anchors,path,"GRCh38","SCENIC+","enhancer_chr","enhancer_start","enhancer_end","target_gene","tf_name","correlation_rna_atac",path.stem+"|")
    # Deduplicate enhancer/source repetitions at variant-gene level before mass.
    dedup={}
    for x in links:
        key=(x["source"],x["source_detail"],x["variant_id_hg38"],x["target_gene"])
        dedup[key]=x
    links=sorted(dedup.values(),key=lambda x:(x["source"],x["target_gene"],int(x["anchor_rank"])))
    link_fields=["source","source_detail","source_score","coordinate_build","anchor_rank","variant_id_hg38","target_gene","uniform_pip"]
    with (args.out_dir/"variant_gene_links.tsv").open("w",newline="") as h:
        w=csv.DictWriter(h,fieldnames=link_fields,delimiter="\t");w.writeheader();w.writerows(links)
    # Link annotations are physical-variant properties.  Re-expand them to the
    # original study+locus fine-mapping units before aggregating PIPs.  The sum
    # within a SuSiE locus is an expected linked causal-variant count, not the
    # probability that at least one linked variant is causal.  Values from
    # different studies/loci are never added together.
    annotated=defaultdict(set)
    for x in links: annotated[x["variant_id_hg38"]].add((x["source"],x["target_gene"]))
    unit_mass=defaultdict(dict)
    unit_meta={}
    for x in support:
        vid=x["variant_id_hg38"]
        for source,gene in annotated.get(vid,set()):
            key=(source,x["study_name"],x["locus_id"],gene)
            unit_mass[key][vid]=float(x["pip"])
            unit_meta[key]=x
    rows=[]
    for (source,study,locus,gene),vals in sorted(unit_mass.items()):
        meta=unit_meta[(source,study,locus,gene)]
        rows.append({"source":source,"study_name":study,"locus_id":locus,
          "ancestry":meta.get("ancestry",""),"target_gene":gene,
          "linked_pip_sum":sum(vals.values()),"n_linked_anchors":len(vals),
          "estimand":"expected linked causal-variant count within this SuSiE study-locus",
          "probability_at_least_one_available":"FALSE",
          "missing_semantics":"absent means unannotated, not unlinked"})
    unit_path=args.out_dir/"source_separated_locus_linked_pip.tsv"
    fields=["source","study_name","locus_id","ancestry","target_gene","linked_pip_sum",
      "n_linked_anchors","estimand","probability_at_least_one_available","missing_semantics"]
    with unit_path.open("w",newline="") as h:
        w=csv.DictWriter(h,fieldnames=fields,delimiter="\t");w.writeheader();w.writerows(rows)

    # A cross-study view uses a maximum plus support counts.  It is explicitly
    # descriptive and cannot be read as a pooled posterior probability.
    grouped=defaultdict(list)
    for row in rows: grouped[(row["source"],row["target_gene"])].append(row)
    summaries=[]
    for (source,gene),vals in sorted(grouped.items()):
        summaries.append({"source":source,"target_gene":gene,
          "max_study_locus_linked_pip_sum":max(float(x["linked_pip_sum"]) for x in vals),
          "n_study_locus_units":len(vals),"n_studies":len({x["study_name"] for x in vals}),
          "n_loci":len({(x["study_name"],x["locus_id"]) for x in vals}),
          "cross_study_sum_prohibited":"TRUE","calibrated_probability":"FALSE",
          "interpretation":"descriptive cross-study support summary"})
    summary_path=args.out_dir/"source_separated_cross_study_summary.tsv"
    summary_fields=["source","target_gene","max_study_locus_linked_pip_sum","n_study_locus_units",
      "n_studies","n_loci","cross_study_sum_prohibited","calibrated_probability","interpretation"]
    with summary_path.open("w",newline="") as h:
        w=csv.DictWriter(h,fieldnames=summary_fields,delimiter="\t");w.writeheader();w.writerows(summaries)
    covered=defaultdict(set)
    for x in links: covered[x["source"]].add(x["variant_id_hg38"])
    coverage=[{"source":s,"anchors_total":len(anchors),"anchors_annotated":len(covered[s]),
      "anchors_unannotated":len(anchors)-len(covered[s]),"missing_means":"unannotated_not_unlinked"} for s in ("ABC","SCENIC+")]
    with (args.out_dir/"source_coverage.tsv").open("w",newline="") as h:
        w=csv.DictWriter(h,fieldnames=list(coverage[0]),delimiter="\t");w.writeheader();w.writerows(coverage)
    atomic_json({"status":"PASS","anchors":len(anchors),"variant_gene_links":len(links),
      "locus_linked_pip_rows":len(rows),"cross_study_summary_rows":len(summaries),
      "pip":"primary uniform-prior SuSiE PIP retained per study and locus",
      "aggregation":"sum unique linked-anchor PIPs only within source+study+locus+target gene; estimand is expected linked causal-variant count",
      "probability_at_least_one_available":False,
      "cross_study_summary":"maximum study-locus linked-PIP sum plus support counts; descriptive, not a pooled posterior",
      "sources_kept_separate":True,"cross_source_sum_prohibited":True,"cross_study_sum_prohibited":True,
      "missing_semantics":"missing means unannotated, never evidence of no link",
      "evidence_class":"predicted_link","apply_only_firewall":True},args.out_dir/"integration_contract.json")
    print(f"PASS: {len(links)} source-separated links, {len(rows)} locus-resolved linked-PIP rows")

if __name__=="__main__":main()
