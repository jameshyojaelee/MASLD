#!/usr/bin/env python3
"""Trace common reporter feature differences without running model inference."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
import pysam

PROJ=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(PROJ/"scripts/analysis/alphagenome_program"))
import i1_extract_mpra_embeddings as old


def main(args):
    rows=old.read_manifest();seqs=old.read_fixture_fasta()
    with args.manifest.open() as stream:
        new={r["element_id"]:r for r in csv.DictReader(stream,delimiter="\t")}
    fasta=pysam.FastaFile(old.FASTA_PATH)
    digest_old=hashlib.sha256();digest_new=hashlib.sha256();failures=[]
    for row in rows:
        n=new[row["element_id"]]
        pos0=int(n["variant_pos1"])-1;start0=pos0-1024
        reference=fasta.fetch(n["contig"],start0,start0+2048).upper()
        alternate=reference[:1024]+n["genomic_alt"]+reference[1025:]
        records,_=old.build_windows(row,seqs,fasta,2048)
        current=dict(zip(old.ALLELE_ORDER,(reference,alternate,old.revcomp(reference),old.revcomp(alternate))))
        reasons=[]
        if start0 != int(row["input_start0"])+1024 or pos0 != int(row["variant_pos0"]):
            reasons.append("input_start_or_variant_coordinate")
        if n["genomic_ref"] != row["ref"] or n["genomic_alt"] != row["alt"]:
            reasons.append("allele_identity")
        for allele in old.ALLELE_ORDER:
            sequence,index=records[allele]
            digest_old.update((row["element_id"]+allele+sequence).encode())
            digest_new.update((row["element_id"]+allele+current[allele]).encode())
            if current[allele] != sequence:reasons.append(allele+"_sequence")
            if index//128 != (7 if allele.endswith("_RC") else 8):reasons.append(allele+"_pool_bin")
        if reasons:failures.append({"element_id":row["element_id"],"mismatches":reasons})
    result={"common_reporters":len(rows),"four_allele_sequences":4*len(rows),"identity_mismatches":failures,
        "archived_sequence_sha256":digest_old.hexdigest(),"new_sequence_sha256":digest_new.hexdigest(),
        "variant_index_forward":1024,"variant_index_reverse":1023,"pool_bin_forward":8,"pool_bin_reverse":7,
        "checkpoint_path_both":str(old.CHECKPOINT),"trunk_loader_both":"i1_extract_mpra_embeddings.load_trunk",
        "state_mode":"same_inference_default_frozen_state",
        "one_hot_dtype":"installed_DNAOneHotEncoder_default_float32; old_numpy_device_put_and_new_jnp_asarray_float32",
        "input_dtype_explanation":"not_supported_both_float32",
        "feature_difference_cause":"unresolved; numerical_kernel_or_compilation_variation_is_an_inference_not_established",
        "scientific_interpretation":"full_source_B1_changes_population_and_numerical_representation; not_a_clean_population_only_comparison"}
    args.out.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    main(parser.parse_args())
