#!/usr/bin/env python3
"""Build a deduplicated genomic scoring set plus target-to-variant map.

Overlapping 501-bp windows share genomic SNVs. ChromBPNet centers on the SNV,
so those SNVs are evaluated once and expanded back to target windows afterward.
"""
import csv, glob, gzip, json, os
from pathlib import Path

ROOT=Path(os.environ.get("MASLD_PROJECT_ROOT","/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"))
SEQ=ROOT/"GWAS/finemapping/results/seqfunc/haplotype_saturation/sequences"
SHARD_SIZE=int(os.environ.get("GENOMIC_SHARD_SIZE","25000"))
sources=sorted(glob.glob(str(SEQ/"saturation_variants.shard[0-9][0-9][0-9].tsv")))
maps=[]; genomic={}
for srcname in sources:
    with open(srcname) as fi:
        for x in csv.DictReader(fi,delimiter="\t"):
            old=x["variant_id"]
            genomic_id=x.get("genomic_variant_id") or (old.split("__",1)[1] if "__" in old else old)
            perturbation_id=x["target_id"]+"__"+genomic_id
            x["variant_id"]=genomic_id; x["genomic_variant_id"]=genomic_id
            x["target_perturbation_id"]=perturbation_id; maps.append(x)
            key=(x["chrom"],x["pos_hg38"],x["ref"],x["alt"])
            if genomic_id in genomic and genomic[genomic_id]!=key:
                raise RuntimeError("inconsistent genomic identity: "+genomic_id)
            genomic[genomic_id]=key
if len(maps)!=751500 or len({x["target_perturbation_id"] for x in maps})!=751500:
    raise SystemExit("target-map identity gate failed")
if len(genomic)!=671187: raise SystemExit(f"expected 671187 unique genomic SNVs; saw {len(genomic)}")

# Immutable full target map; original per-target shards remain provenance inputs.
fields=list(maps[0]);
for col in ("genomic_variant_id","target_perturbation_id"):
    if col not in fields: fields.insert(fields.index("variant_id")+1,col)
tmp=SEQ/"target_variant_map.tsv.gz.tmp"
with gzip.open(tmp,"wt",newline="") as fo:
    w=csv.DictWriter(fo,fieldnames=fields,delimiter="\t"); w.writeheader(); w.writerows(maps)
os.replace(tmp,SEQ/"target_variant_map.tsv.gz")

# Deterministic genomic ordering and 27 scoring shards.
ordered=sorted(genomic.items(),key=lambda z:(int(z[1][0].removeprefix("chr")) if z[1][0].removeprefix("chr").isdigit() else 99,z[1][1],z[1][2],z[1][3]))
for old in glob.glob(str(SEQ/"genomic_variants.shard*")): os.remove(old)
for start in range(0,len(ordered),SHARD_SIZE):
    part=ordered[start:start+SHARD_SIZE]; si=start//SHARD_SIZE
    tab=SEQ/f"genomic_variants.shard{si:03d}.tsv"; cbp=SEQ/f"genomic_variants.shard{si:03d}.chrombpnet.tsv"
    with tab.open("w",newline="") as fo, cbp.open("w") as fc:
        w=csv.writer(fo,delimiter="\t"); w.writerow(("chrom","pos_hg38","ref","alt","variant_id"))
        for vid,(chrom,pos,ref,alt) in part:
            w.writerow((chrom,pos,ref,alt,vid)); fc.write("\t".join((chrom,pos,ref,alt,vid))+"\n")

contract_path=SEQ/"sequence_contract.json"; contract=json.loads(contract_path.read_text())
for stale in ("unique_perturbation_ids","perturbation_identity","overlapping_window_contexts_are_distinct"):
    contract.pop(stale,None)
contract.update({"target_variant_rows":len(maps),"unique_target_perturbation_ids":751500,
 "unique_genomic_snvs_to_score":len(genomic),"genomic_scoring_shards":27,
 "scoring_identity":"genomic_variant_id","expansion_map":"target_variant_map.tsv.gz",
 "overlapping_genomic_snvs_scored_once":True})
tmp=contract_path.with_suffix(".json.tmp"); tmp.write_text(json.dumps(contract,indent=2)+"\n"); os.replace(tmp,contract_path)
print(f"PASS target rows={len(maps)} unique genomic SNVs={len(genomic)} shards=27")
