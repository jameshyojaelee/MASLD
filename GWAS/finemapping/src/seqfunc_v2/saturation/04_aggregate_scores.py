#!/usr/bin/env python3
"""Audit and aggregate held-out ChromBPNet saturation scores.

Only descriptive ranks and empirical percentiles are produced. Scorer p-values
are intentionally not propagated: the 30,000 shuffled variants are not a
genome-wide inferential null and therefore do not support formal FDR claims.
"""
from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import math
from collections import defaultdict
from pathlib import Path

from common import OUT, atomic_json, read_json, sha256


def open_text(path, mode="rt"):
    return gzip.open(path, mode, newline="") if str(path).endswith(".gz") else path.open(mode, newline="")


def fnum(row, names, required=True):
    for name in names:
        if name in row and row[name] not in ("", "NA", "nan"):
            value = float(row[name])
            if not math.isfinite(value): raise ValueError(f"nonfinite {name}")
            return value
    if required: raise KeyError(f"none of score columns present: {names}")
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=OUT)
    args = ap.parse_args(); seq=args.root/"sequences"; score_dir=args.root/"scores"; out=args.root/"aggregate"
    out.mkdir(parents=True,exist_ok=True)
    gate=read_json(args.root/"gates/upstream_gate.json")
    if not gate.get("saturation_submission_allowed",False): raise SystemExit("upstream gate BLOCKED")
    required_contracts = [args.root/"anchors/anchor_contract.json", seq/"sequence_contract.json",
                          args.root/"integration/integration_contract.json"]
    for contract_path in required_contracts:
        if not contract_path.is_file() or read_json(contract_path).get("status") != "PASS":
            raise SystemExit(f"missing or failed preparation/integration contract: {contract_path}")
    with (seq/"scoring_shards.tsv").open() as h: shards=list(csv.DictReader(h,delimiter="\t"))
    score_by_id={}; shard_qc=[]
    for shard in shards:
        sid=shard["shard_id"]; expected_path=Path(shard["variant_table"])
        result=score_dir/f"{sid}.variant_scores.tsv"; provenance=score_dir/f"{sid}.provenance.json"
        for path in (expected_path,result,provenance):
            if not path.is_file() or path.stat().st_size==0: raise SystemExit(f"missing shard artifact: {path}")
        prov=read_json(provenance)
        if prov.get("status")!="PASS" or prov.get("fold_id")!=shard["fold_id"] or not prov.get("apply_only_firewall"):
            raise SystemExit(f"invalid provenance: {provenance}")
        with expected_path.open() as h: expected={r["variant_id"]:r for r in csv.DictReader(h,delimiter="\t")}
        with result.open() as h: observed=list(csv.DictReader(h,delimiter="\t"))
        if len(observed)!=len(expected) or len({r["variant_id"] for r in observed})!=len(observed):
            raise SystemExit(f"row/identity failure: {sid}")
        for row in observed:
            vid=row["variant_id"]
            if vid not in expected: raise SystemExit(f"unexpected variant {vid} in {sid}")
            e=expected[vid]
            chrom=row.get("chr",row.get("chrom")); pos=row.get("pos",row.get("pos_hg38"))
            a1=row.get("allele1",row.get("ref")); a2=row.get("allele2",row.get("alt"))
            if (chrom,str(pos),a1,a2)!=(e["chrom"],e["pos_hg38"],e["ref"],e["alt"]):
                raise SystemExit(f"coordinate/allele failure: {vid}")
            count_signed=fnum(row,["logfc","count_logfc","count_effect"])
            count_abs=fnum(row,["abs_logfc","abs_count_logfc"],required=False)
            if count_abs is None: count_abs=abs(count_signed)
            if abs(count_abs-abs(count_signed))>1e-5: raise SystemExit(f"absolute count inconsistency: {vid}")
            profile=fnum(row,["jsd","profile_jsd"])
            profile_signed=fnum(row,["profile_logfc","signed_profile_effect"],required=False)
            score_by_id[vid]={"count_effect_signed":count_signed,"count_effect_abs":count_abs,
                              "profile_jsd":profile,"profile_effect_abs":profile,"profile_effect_signed":profile_signed,
                              "model_fold_id":shard["fold_id"],"score_shard_id":sid}
        shard_qc.append({"shard_id":sid,"fold_id":shard["fold_id"],"n_expected":len(expected),
                         "n_scored":len(observed),"status":"PASS","score_sha256":sha256(result)})
    if len(score_by_id)!=sum(int(x["n_variants"]) for x in shards): raise SystemExit("global unique score count failure")

    metric_names=["count_effect_signed","count_effect_abs","profile_jsd"]
    global_sorted={m:sorted(x[m] for x in score_by_id.values()) for m in metric_names}
    by_target=defaultdict(list)
    map_path=seq/"target_variant_map.tsv.gz"
    with gzip.open(map_path,"rt") as h:
        for row in csv.DictReader(h,delimiter="\t"):
            if row["genomic_variant_id"] not in score_by_id: raise SystemExit(f"unscored map variant: {row['genomic_variant_id']}")
            by_target[row["target_id"]].append(row)
    detail_fields=["target_perturbation_id","genomic_variant_id","target_id","anchor_rank","anchor_variant_id_hg38",
      "chrom","pos_hg38","ref","alt","offset_from_anchor","is_observed_anchor_alt","model_fold_id",
      "uniform_max_pip","count_effect_signed","count_effect_abs","profile_jsd","profile_effect_abs","profile_effect_signed",
      "rank_target_count_signed","rank_target_count_abs","rank_target_profile_jsd","empirical_percentile_global_count_signed","empirical_percentile_global_count_abs",
      "empirical_percentile_global_profile_jsd","evidence_class","formal_fdr_available"]
    summaries=[]; detail_path=out/"saturation_scores.tsv.gz"; tmp=out/"saturation_scores.tsv.tmp.gz"
    with gzip.open(tmp,"wt",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=detail_fields,delimiter="\t");writer.writeheader()
        for target in sorted(by_target,key=lambda x:int(by_target[x][0]["anchor_rank"])):
            rows=by_target[target]
            if len(rows)!=1503: raise SystemExit(f"target {target} has {len(rows)} rather than 1503 substitutions")
            vals=[]
            for row in rows:
                z=dict(row); z.update(score_by_id[row["genomic_variant_id"]]); vals.append(z)
            order_abs={id(z):rank for rank,z in enumerate(sorted(vals,key=lambda x:(-x["count_effect_abs"],x["genomic_variant_id"])),1)}
            order_signed={id(z):rank for rank,z in enumerate(sorted(vals,key=lambda x:(-x["count_effect_signed"],x["genomic_variant_id"])),1)}
            order_jsd={id(z):rank for rank,z in enumerate(sorted(vals,key=lambda x:(-x["profile_jsd"],x["genomic_variant_id"])),1)}
            for z in vals:
                z["rank_target_count_signed"]=order_signed[id(z)]; z["rank_target_count_abs"]=order_abs[id(z)]; z["rank_target_profile_jsd"]=order_jsd[id(z)]
                z["empirical_percentile_global_count_signed"]=bisect.bisect_right(global_sorted["count_effect_signed"],z["count_effect_signed"])/len(score_by_id)
                z["empirical_percentile_global_count_abs"]=bisect.bisect_right(global_sorted["count_effect_abs"],z["count_effect_abs"])/len(score_by_id)
                z["empirical_percentile_global_profile_jsd"]=bisect.bisect_right(global_sorted["profile_jsd"],z["profile_jsd"])/len(score_by_id)
                z["evidence_class"]="predicted"; z["formal_fdr_available"]="FALSE"; writer.writerow({k:z.get(k,"") for k in detail_fields})
            anchor=[z for z in vals if z["is_observed_anchor_alt"].upper()=="TRUE"]
            if len(anchor)!=1: raise SystemExit(f"observed anchor allele is not unique for {target}")
            top_count=max(vals,key=lambda x:x["count_effect_abs"]); top_profile=max(vals,key=lambda x:x["profile_jsd"]); a=anchor[0]
            summaries.append({"anchor_rank":a["anchor_rank"],"target_id":target,"anchor_variant_id_hg38":a["anchor_variant_id_hg38"],
              "uniform_max_pip":a["uniform_max_pip"],"anchor_count_effect_signed":a["count_effect_signed"],
              "anchor_count_effect_abs":a["count_effect_abs"],"anchor_profile_jsd":a["profile_jsd"],"anchor_profile_effect_abs":a["profile_effect_abs"],
              "anchor_rank_target_count_signed":order_signed[id(a)],"anchor_rank_target_count_abs":order_abs[id(a)],"anchor_rank_target_profile_jsd":order_jsd[id(a)],
              "top_count_variant":top_count["genomic_variant_id"],"top_count_effect_abs":top_count["count_effect_abs"],
              "top_profile_variant":top_profile["genomic_variant_id"],"top_profile_jsd":top_profile["profile_jsd"],
              "formal_fdr_available":"FALSE","interpretation":"mechanistic_nomination_not_causal_evidence"})
    tmp.replace(detail_path)
    for name,rows in (("shard_qc.tsv",shard_qc),("saturation_anchor_summary.tsv",summaries)):
        with (out/name).open("w",newline="") as h:
            w=csv.DictWriter(h,fieldnames=list(rows[0]),delimiter="\t");w.writeheader();w.writerows(rows)
    atomic_json({"status":"PASS","n_shards":len(shards),"n_unique_genomic_substitutions":len(score_by_id),
      "n_target_substitution_rows":sum(map(len,by_target.values())),"n_anchors":len(by_target),
      "metrics":{"count_signed":"log fold change, ALT relative to REF","count_absolute":"absolute log fold change",
                 "profile_absolute":"Jensen-Shannon distance (nonnegative divergence)","profile_signed":"reported only if emitted natively by scorer; never imputed from count direction"},
      "uncertainty":"descriptive empirical percentile against the scored saturation universe",
      "formal_p_values_reported":False,"formal_fdr_reported":False,"reason_no_fdr":"30000 shuffled variants are not a genome-wide inferential null",
      "evidence_class":"predicted","apply_only_firewall":True,"canonical_outputs_mutated":False},out/"aggregate_contract.json")
    atomic_json({"status":"PASS","complete":True,"upstream_unlocked":True,
      "n_anchors":len(by_target),"n_unique_genomic_substitutions":len(score_by_id),
      "n_target_substitution_rows":sum(map(len,by_target.values())),
      "formal_fdr_reported":False,"source_link_scores_summed":False,
      "interpretation":"mechanistic_nomination_not_causal_evidence",
      "apply_only_firewall":True,"canonical_outputs_mutated":False},args.root/"gate_verdict.json")
    print(f"PASS: {len(score_by_id)} unique scores expanded to {sum(map(len,by_target.values()))} target rows")


if __name__=="__main__": main()
