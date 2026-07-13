#!/usr/bin/env python3
"""Generate reference-validated, model-ready saturation SNVs.

This is an apply-only experimental-design artifact. It creates no learned
features, never calls AlphaGenome, and cannot update fine-mapping or convergence.
The output uses hg38 REF alleles after liftover, with a hard FASTA assertion.
"""
import argparse, csv, gzip, hashlib, json, os
from pathlib import Path

BASES = "ACGT"

def open_text(path, mode="rt"):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--n-targets", type=int, default=500)
    ap.add_argument("--shard-size", type=int, default=25000)
    ap.add_argument("--sequence-window", type=int, default=501)
    a=ap.parse_args()
    if a.sequence_window % 2 != 1: raise SystemExit("--sequence-window must be odd")
    import pysam
    out=Path(a.out_dir); out.mkdir(parents=True,exist_ok=True)
    fa=pysam.FastaFile(a.fasta)
    with open(a.manifest) as f: rows=list(csv.DictReader(f,delimiter="\t"))[:a.n_targets]
    emitted=[]; rejected=[]; fasta_records=[]
    half=a.sequence_window//2
    for rank,r in enumerate(rows,1):
        chrom=r["chr"]; pos=int(r["pos_hg38"]); ref=r["ref_hg38"].upper()
        observed=fa.fetch(chrom,pos-1,pos).upper()
        if len(ref)!=1 or ref not in BASES or observed != ref:
            rejected.append({"rank":rank,"chr":chrom,"pos_hg38":pos,"expected_ref":ref,
                             "fasta_ref":observed,"reason":"hg38_reference_mismatch"}); continue
        start=pos-1-half; end=start+a.sequence_window
        if start < 0 or end > fa.get_reference_length(chrom):
            rejected.append({"rank":rank,"chr":chrom,"pos_hg38":pos,"expected_ref":ref,
                             "fasta_ref":observed,"reason":"window_outside_contig"}); continue
        seq=fa.fetch(chrom,start,end).upper()
        if seq[half] != ref: raise RuntimeError("internal center/reference error")
        target=f"sat{rank:04d}_{chrom}_{pos}_{ref}"
        fasta_records.append((target+"_REF",seq))
        for offset,base in enumerate(seq):
            if base not in BASES: continue
            mut_pos=start+offset+1
            for alt in BASES:
                if alt==base: continue
                genomic_vid=f"{chrom}_{mut_pos}_{base}_{alt}"
                # Perturbation identity must include the target window: candidate
                # windows may overlap, while each target-context score is distinct.
                vid=f"{target}__{genomic_vid}"
                emitted.append({"chrom":chrom,"pos_hg38":mut_pos,"ref":base,"alt":alt,
                    "variant_id":vid,"genomic_variant_id":genomic_vid,
                    "target_id":target,"target_rank":rank,
                    "offset_from_anchor":mut_pos-pos,"anchor_pos_hg38":pos,
                    "primary_max_pip":r["primary_max_pip"],"primary_n_studies":r["primary_n_studies"],
                    "source_variant_id_hg19":r["variant_id_hg19"],
                    "is_credible_allele":str(mut_pos==pos and alt==r["alt_hg38"]).upper(),
                    "apply_only":"TRUE","alphagenome_training_prohibited":"TRUE"})
                mut=seq[:offset]+alt+seq[offset+1:]
                fasta_records.append((target+f"_{mut_pos}_{alt}",mut))
    fields=list(emitted[0]) if emitted else []
    for si in range(0,len(emitted),a.shard_size):
        part=emitted[si:si+a.shard_size]; p=out/f"saturation_variants.shard{si//a.shard_size:03d}.tsv"
        with open(p,"w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=fields,delimiter="\t"); w.writeheader(); w.writerows(part)
        # Headerless schema consumed by the vendored Kundaje variant scorer.
        with open(out/f"saturation_variants.shard{si//a.shard_size:03d}.chrombpnet.tsv","w") as f:
            for x in part:
                f.write("\t".join((x["chrom"],str(x["pos_hg38"]),x["ref"],x["alt"],x["variant_id"]))+"\n")
    with gzip.open(out/"saturation_sequences.fa.gz","wt") as f:
        for name,seq in fasta_records: f.write(f">{name}\n{seq}\n")
    rejfields=["rank","chr","pos_hg38","expected_ref","fasta_ref","reason"]
    with open(out/"reference_rejections.tsv","w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=rejfields,delimiter="\t"); w.writeheader(); w.writerows(rejected)
    contract={"status":"prepared","coordinate_system":"GRCh38","targets_requested":len(rows),
      "targets_passing_reference":len(rows)-len(rejected),"targets_rejected":len(rejected),
      "snvs_emitted":len(emitted),"sequence_window":a.sequence_window,
      "scorer_schema":"headerless chrom,pos1,ref,alt,variant_id; one file per shard",
      "apply_only_firewall":True,"alphagenome_mode":"zero_shot_only_no_training_or_calibration",
      "observed_haplotype_claim_allowed":False,"manifest_sha256":hashlib.sha256(Path(a.manifest).read_bytes()).hexdigest()}
    (out/"sequence_contract.json").write_text(json.dumps(contract,indent=2)+"\n")
    if rejected: raise SystemExit(f"hard reference gate failed for {len(rejected)} targets; outputs quarantined")
    print(json.dumps(contract,indent=2))
if __name__=="__main__": main()
