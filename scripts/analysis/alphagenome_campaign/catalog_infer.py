#!/usr/bin/env python3
"""Local variant/phased-SNV inference with explicit molecular scope and source phase.

Run on an allocated L40S using alphagenome_local. Only static liver ATAC/DNase
tracks are supported here. No donor context, disease risk or validated molecular
interaction is inferred from a constructed haplotype. Indel haplotypes require
an independently validated coordinate/readout mapping and are refused.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import sys

from catalog_query import variant_id

ROOT = Path(__file__).resolve().parents[3]


def parse_request(variants, phase, phase_source, length):
    if length not in (2048, 16384, 1048576):
        raise ValueError("Supported input lengths are 2048, 16384 and 1048576")
    ids = [variant_id(v) for v in variants]
    if not ids or len(set(ids)) != len(ids):
        raise ValueError("Provide distinct explicit variant identities")
    if len(ids) > 1 and (phase != "same_haplotype" or not phase_source.strip()):
        raise ValueError("unknown_phase: phased inference requires a declared same-haplotype source")
    parsed = []
    for key in ids:
        build, chrom, pos, ref, alt = key.split(":")
        if len(ref) != 1 or len(alt) != 1:
            raise ValueError("indel_haplotype_readout_mapping_not_validated; use the separately recorded indel route")
        parsed.append((chrom, int(pos), ref, alt))
    parsed.sort()
    if len({p[0] for p in parsed}) != 1 or len({p[1] for p in parsed}) != len(parsed):
        raise ValueError("Local haplotypes require one chromosome and nonconflicting positions")
    center0 = (parsed[0][1] + parsed[-1][1]) // 2 - 1
    start0 = center0 - length // 2
    if start0 < 0 or parsed[0][1]-1 < start0+250 or parsed[-1][1]-1 >= start0+length-250:
        raise ValueError("Variants/readout do not fit the supported window")
    return ids, parsed, start0, start0 + length


def construct(fetch, parsed, start0, end0):
    chrom = parsed[0][0]
    reference = fetch(chrom, start0, end0).upper()
    if len(reference) != end0-start0 or set(reference) - set("ACGT"):
        raise ValueError("incomplete_or_ambiguous_reference_window")
    alternate = list(reference)
    for _, pos1, ref, alt in parsed:
        offset = pos1-1-start0
        if not 0 <= offset < len(reference) or reference[offset] != ref:
            raise ValueError(f"reference_mismatch:{chrom}:{pos1}:{ref}")
        alternate[offset] = alt
    return reference, "".join(alternate)


def native_readout(output, lo, hi):
    import numpy as np
    records = {}
    for assay, attr in (("ATAC", "atac"), ("DNase", "dnase")):
        track = getattr(output, attr)
        values = np.asarray(track.values, dtype=np.float32)
        if values.ndim != 2 or values.shape[0] < hi or not np.isfinite(values[lo:hi]).all():
            raise ValueError("Missing or nonfinite native track values")
        scores = np.log2(1 + np.sum(values[lo:hi], axis=0, dtype=np.float32))
        metadata = track.metadata.to_dict(orient="records")
        if len(metadata) != len(scores) or not len(scores):
            raise ValueError("Missing or mismatched native track identities")
        for i, (score, record) in enumerate(zip(scores, metadata)):
            name = str(record.get("name", i))
            if (assay, name) in records:
                raise ValueError("Duplicate native track identity")
            records[(assay, name)] = {"value": float(score), "track_metadata": {
                str(k): str(v) for k, v in record.items()}, "assay": assay, "track": name}
    return records


def predict(variants, phase, phase_source, length, out, check_repeat=False):
    ids, parsed, start0, end0 = parse_request(variants, phase, phase_source, length)
    if out.exists():
        raise FileExistsError(out)
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Scientific inference requires an allocated compute node")
    sys.path.insert(0, str(ROOT / "scripts/analysis/alphagenome_program"))
    import h1_runtime_common as H
    import i1_common as C
    import pysam
    import numpy as np
    from alphagenome.data import genome
    from alphagenome.models import dna_client
    if C.gpu_card().get("card_tag") != "l40s":
        raise RuntimeError("This recorded inference route requires an L40S")
    with pysam.FastaFile(H.FASTA_PATH) as fasta:
        ref, alt = construct(fasta.fetch, parsed, start0, end0)
    model, device, restore_seconds = H.load_model()
    interval = genome.Interval(chromosome=parsed[0][0], start=start0, end=end0)
    # One fixed source-coordinate readout covers the entire SNV haplotype.
    lo = parsed[0][1]-1-start0-250
    hi = parsed[-1][1]-start0+250
    def call(sequence):
        return native_readout(model.predict_sequence(sequence=sequence, interval=interval,
            ontology_terms=H.LIVER_TERMS,
            requested_outputs=[dna_client.OutputType.ATAC, dna_client.OutputType.DNASE]), lo, hi)
    r, a = call(ref), call(alt)
    if r.keys() != a.keys():
        raise ValueError("REF/ALT track identities disagree")
    repeated_max = None
    if check_repeat:
        repeat = call(ref)
        repeated_max = max(abs(repeat[k]["value"]-r[k]["value"]) for k in r)
        np.testing.assert_allclose([repeat[k]["value"] for k in r], [r[k]["value"] for k in r], atol=1e-5, rtol=1e-5)
    rows = []
    for key in sorted(r):
        effect = float(np.float32(a[key]["value"])-np.float32(r[key]["value"]))
        if not math.isfinite(effect):
            raise ValueError("Nonfinite molecular effect")
        rows.append({"assay":key[0], "track":key[1], "ref_value":r[key]["value"], "alt_value":a[key]["value"],
            "effect":effect, "unit":"ALT_minus_REF_log2_1_plus_predicted_track_sum", "kind":"predicted",
            "track_metadata":r[key]["track_metadata"], "uncertainty":None, "uncertainty_reason":"not_calibrated_for_this_request",
            "source_model":str(H.CHECKPOINT), "context":"packaged_liver_panel; life_stage_retained_per_track",
            "unsupported_reason":"donor_or_disease_condition_effect_not_evaluated"})
    receipt = {"variants":ids, "phase":phase if len(ids)>1 else "not_applicable_single_variant",
        "phase_source":phase_source, "phase_verification":"caller_declared; not_inferred_from_population_LD",
        "sequence_background":"GRCh38_reference_with_declared_SNV_edits", "input_length":length,
        "window":{"chrom":parsed[0][0],"start0":start0,"end0":end0},
        "readout":{"start0":start0+lo,"end0":start0+hi,"coordinate_convention":"zero_based_half_open"},
        "reference_sha256":hashlib.sha256(ref.encode()).hexdigest(), "alternate_sha256":hashlib.sha256(alt.encode()).hexdigest(),
        "effects":rows, "release_state":"candidate_not_adopted", "biological_interaction_claim":False,
        "individual_MASLD_risk":"not_supported", "terms":"AlphaGenome_noncommercial; https://huggingface.co/google/alphagenome-all-folds",
        "calibration":"none; source_native_model_scale", "repeat_max_abs_track_change":repeated_max,
        "memory":H.memory_stats(device), "restore_seconds":restore_seconds,
        "job_id":os.environ["SLURM_JOB_ID"]}
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("x") as handle:
        json.dump(receipt, handle, indent=2, allow_nan=False)
        handle.write("\n")
    print(json.dumps({"out":str(out),"tracks":len(rows),"variants":len(ids),"risk_score":False}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("variants", nargs="+")
    parser.add_argument("--phase", choices=("unknown","same_haplotype"), default="unknown")
    parser.add_argument("--phase-source", default="")
    parser.add_argument("--length", type=int, default=2048)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--check-repeat", action="store_true")
    args = parser.parse_args()
    predict(args.variants,args.phase,args.phase_source,args.length,args.out,args.check_repeat)


if __name__ == "__main__":
    main()
