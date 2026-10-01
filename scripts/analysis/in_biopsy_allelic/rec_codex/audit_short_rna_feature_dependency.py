"""Audit frozen raw profile dependence on the 40 k31-unobservable features.

Only the allowed 99-person development RNA and frozen model are opened.
Removing counts is a developmental sensitivity calculation, not a proposed
receiving input recipe. No H3 values, receiving values or prediction refits.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import platform
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[4]
BENCH = ROOT / "Analysis/MASLD_Model_Benchmark"
FIX = BENCH / "executions/model-data-064-21079902/fixture/molecular"
REC = BENCH / "executions/codex-rec-20260929T142434Z"
RELEASE = BENCH / "release/masld-liver-chromatin-state-v1.2"
SHORT = REC / "ensembl98_salmon_metadata_21998488/modeled_genes_no_unclipped_k31_transcript.tsv"
WEIGHTS = RELEASE / "weights/chromatin_state_v1_2.npz"
SCORER = RELEASE / "score.py"
GUARDS = {
    FIX / "rna_values.npy": "812d5fce6e7cc88f6939ae493db7b7e7470057c21cf3b219abc39f587ff64ddc",
    FIX / "rna_feature_axis.tsv": "baf983ad18a893b7e5db04364818073facc68b47ae25f63070373481e3b217dd",
    FIX / "participant_axis.tsv": "7e5be1a16df505e5daad451398e8b1a4b3cbbc48a32cd1381113e1ae4b7347bc",
    SHORT: "4f348047d98367f9e041716dbe9a7a8952adcf54bba78c0dc97f55477151f7bf",
    WEIGHTS: "72f9eaba359f7c4c61cac2e00da2d5539b5d3ec868dc33b1bacc2ef0105694ec",
    SCORER: "0c9bff99757dcc934529117209397bb7dc9ef60d2293cba1a5b0eb48fadecc39",
}
BLOCK = 2048


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def table(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write(path, rows):
    with path.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def logcpm(counts):
    totals = counts.sum(1)
    require(np.isfinite(totals).all() and (totals > 0).all(), "Invalid modeled-axis count totals")
    return np.log2(1 + counts / totals[:, None] * 1e6)


def local_delta(delta, indices, coefficients, scale, lo, hi):
    out = np.zeros((len(delta), hi - lo))
    for k in range(indices.shape[1]):
        ids = indices[lo:hi, k]
        hit = ids >= 0
        if hit.any():
            out[:, hit] += delta[:, ids[hit]] * (coefficients[lo:hi, k][hit] / scale[ids[hit]])
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    require(os.environ.get("SLURM_JOB_ID"), "Read developmental values only inside compute job")
    require(not args.output.exists(), "Refusing output overwrite")
    hashes = {str(p): sha(p) for p in GUARDS}
    require(all(hashes[str(p)] == expected for p, expected in GUARDS.items()), "Frozen audit input changed")
    participants, features, short = table(FIX / "participant_axis.tsv"), table(FIX / "rna_feature_axis.tsv"), table(SHORT)
    require(len(participants) == len({p["participant_id"] for p in participants}) == 99, "Participant roster changed")
    require(len(features) == 42163 and len(short) == 40, "Feature/short census changed")
    short_idx = np.array([int(r["rna_feature_index"]) for r in short], dtype=int)
    require(len(set(short_idx)) == 40, "Repeated short feature")
    for row, idx in zip(short, short_idx):
        require(0 <= idx < 42163 and features[idx]["stable_gene_id"] == row["stable_gene_id"]
                and int(row["eligible_unclipped_k31"]) == 0, "Short feature identity mismatch")
    counts = np.load(FIX / "rna_values.npy", allow_pickle=False)
    require(counts.shape == (99, 42163) and counts.dtype == np.float64, "RNA fixture shape/dtype changed")
    require(np.isfinite(counts).all() and (counts >= 0).all(), "Invalid development counts")
    weights = np.load(WEIGHTS, allow_pickle=True)
    axis = np.asarray(weights["gene_axis"]).astype(str)
    require(np.array_equal(axis, [r["stable_gene_id"] for r in features]), "Release RNA axis differs")
    require(int(weights["prof_n_train"]) == 99, "Profile training n changed")
    # The raw RRR plus Ocis arm is affine in log2 CPM only if the state
    # coefficients are identically zero. Refuse, rather than drop a route.
    state_coef = np.asarray(weights["prof_Ocis_coef_glob"], dtype=float)
    require(np.isfinite(state_coef).all() and np.count_nonzero(state_coef) == 0,
            "Ocis has a state/global route: current affine derivation is incomplete")
    rrr_scale = np.asarray(weights["prof_rrr_xs"], dtype=float)
    cis_scale = np.asarray(weights["prof_x_sd"], dtype=float)
    vtdv = np.asarray(weights["prof_rrr_VtDv"], dtype=float)
    wvr = np.asarray(weights["prof_rrr_WVr"], dtype=float)
    vr = np.asarray(weights["prof_rrr_Vr"], dtype=float)
    idx = np.asarray(weights["prof_Ocis_cis_idx"])
    cis = np.asarray(weights["prof_Ocis_coef_cis"], dtype=float)
    region_axis = np.asarray(weights["prof_region_key"]).astype(str)
    require(rrr_scale.shape == cis_scale.shape == (42163,) and (rrr_scale > 0).all()
            and (cis_scale > 0).all(), "Invalid frozen feature scales")
    require(vtdv.shape[0] == 42163 and vtdv.shape[1] == wvr.shape[0]
            and wvr.shape[1] == vr.shape[1] == 20 and vr.shape[0] == len(region_axis) == 96460,
            "Frozen RRR dimensions changed")
    require(idx.shape == cis.shape and idx.shape[0] == 96460 and np.issubdtype(idx.dtype, np.integer)
            and ((idx >= -1) & (idx < 42163)).all(), "Invalid frozen local feature axis")
    for value in (rrr_scale, cis_scale, vtdv, wvr, vr, cis):
        require(np.isfinite(value).all(), "Nonfinite frozen model coefficient")
    latent_coef = (vtdv @ wvr) / rrr_scale[:, None]
    x = logcpm(counts)
    # Direct log-expression change keeps the original denominator fixed.
    direct = np.zeros_like(x)
    direct[:, short_idx] = -x[:, short_idx]
    reduced_counts = counts.copy()
    reduced_counts[:, short_idx] = 0
    total_delta = logcpm(reduced_counts) - x
    denominator_delta = total_delta - direct
    del reduced_counts
    latent_delta = {"direct": direct @ latent_coef,
                    "denominator": denominator_delta @ latent_coef,
                    "total": total_delta @ latent_coef}
    cases = ("direct", "denominator", "total")
    deltas = dict(direct=direct, denominator=denominator_delta, total=total_delta)
    metrics = {name: {arm: dict(sum_sq=np.zeros(99), max_abs=np.zeros(99))
                      for arm in ("rrr", "rrr_offset_cis")} for name in cases}
    stats = {name: dict(nonzero=np.zeros(40, dtype=int), max_abs=np.zeros(40), sum_sq=np.zeros(40))
             for name in ("rrr", "cis", "rrr_offset_cis")}
    split_max_error = 0.0
    for lo in range(0, 96460, BLOCK):
        hi = min(lo + BLOCK, 96460)
        global_coef = latent_coef[short_idx] @ vr[lo:hi].T
        local_coef = np.zeros_like(global_coef)
        for position, gene in enumerate(short_idx):
            # Sum every matching local slot; do not assume indices unique.
            local_coef[position] = np.sum(np.where(idx[lo:hi] == gene, cis[lo:hi], 0), axis=1) / cis_scale[gene]
        for name, coef in (("rrr", global_coef), ("cis", local_coef), ("rrr_offset_cis", global_coef + local_coef)):
            stats[name]["nonzero"] += np.count_nonzero(coef, axis=1)
            stats[name]["max_abs"] = np.maximum(stats[name]["max_abs"], np.abs(coef).max(1))
            stats[name]["sum_sq"] += np.square(coef).sum(1)
        block = {}
        for name in cases:
            glob = latent_delta[name] @ vr[lo:hi].T
            local = local_delta(deltas[name], idx, cis, cis_scale, lo, hi)
            block[name] = {"rrr": glob, "rrr_offset_cis": glob + local}
            for arm, values in block[name].items():
                metrics[name][arm]["sum_sq"] += np.square(values).sum(1)
                metrics[name][arm]["max_abs"] = np.maximum(metrics[name][arm]["max_abs"], np.abs(values).max(1))
        for arm in ("rrr", "rrr_offset_cis"):
            split_max_error = max(split_max_error, float(np.abs(block["total"][arm]
                                  - block["direct"][arm] - block["denominator"][arm]).max()))
    require(split_max_error <= 1e-8, "Direct plus denominator split disagrees")
    # Independently compare the frozen public scorer to the affine derivation.
    # Add one synthetic count to every short feature in two development copies,
    # so this check is nonvacuous even when the original forty counts are zero.
    spec = importlib.util.spec_from_file_location("frozen_raw_profile_scorer", SCORER)
    scorer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(scorer)
    probe = counts[:2].copy()
    probe[:, short_idx] += 1
    probe_delta = logcpm(probe) - x[:2]
    probe_latent = probe_delta @ latent_coef
    scorer_errors = {}
    for arm in ("rrr", "rrr_offset_cis"):
        old = scorer.profile(counts[:2].T, axis, str(WEIGHTS), form=arm, transport="raw", min_coverage=1.0)["profile"]
        new = scorer.profile(probe.T, axis, str(WEIGHTS), form=arm, transport="raw", min_coverage=1.0)["profile"]
        pred_delta = probe_latent @ vr.T
        if arm == "rrr_offset_cis":
            pred_delta += local_delta(probe_delta, idx, cis, cis_scale, 0, 96460)
        scorer_errors[arm] = float(np.abs(new - old - pred_delta).max())
        require(scorer_errors[arm] <= 1e-8, "Effective coefficient differs from actual frozen scorer")
    rows = []
    for j, (metadata, gene) in enumerate(zip(short, short_idx)):
        values = counts[:, gene]
        row = dict(stable_gene_id=metadata["stable_gene_id"], rna_feature_index=int(gene),
                   training_nonzero_participants=int(np.count_nonzero(values)), training_count_min=float(values.min()),
                   training_count_max=float(values.max()), training_count_sum=float(values.sum()),
                   training_logcpm_min=float(x[:, gene].min()), training_logcpm_max=float(x[:, gene].max()),
                   rrr_latent_row_exactly_zero=bool(np.count_nonzero(latent_coef[gene]) == 0),
                   local_slots=int(np.count_nonzero(idx == gene)))
        for arm in stats:
            row[arm + "_nonzero_effective_region_coefficients"] = int(stats[arm]["nonzero"][j])
            row[arm + "_coefficient_max_abs_per_log2cpm"] = float(stats[arm]["max_abs"][j])
            row[arm + "_coefficient_rms_per_log2cpm"] = float(np.sqrt(stats[arm]["sum_sq"][j] / 96460))
        rows.append(row)
    donor_rows = []
    mass = counts[:, short_idx].sum(1) / counts.sum(1)
    for i, person in enumerate(participants):
        row = dict(participant_id=person["participant_id"], short_feature_count_fraction=float(mass[i]))
        for name in cases:
            for arm in metrics[name]:
                row[name + "_" + arm + "_prediction_delta_rms"] = float(np.sqrt(metrics[name][arm]["sum_sq"][i] / 96460))
                row[name + "_" + arm + "_prediction_delta_max_abs"] = float(metrics[name][arm]["max_abs"][i])
        # Conditional bound, only for the defined direct perturbation:
        # ||sum_g dx_g B_g|| / sqrt(R) <= sum_g |dx_g| ||B_g|| / sqrt(R).
        for arm in ("rrr", "rrr_offset_cis"):
            row["direct_" + arm + "_rms_triangle_upper_bound"] = float(np.abs(direct[i, short_idx])
                    @ np.sqrt(stats[arm]["sum_sq"] / 96460))
        donor_rows.append(row)
    args.output.mkdir(parents=True, exist_ok=False)
    write(args.output / "short_feature_dependence.tsv", rows)
    write(args.output / "developmental_perturbation.tsv", donor_rows)
    summary = dict(schema_version="short-rna-frozen-profile-dependency-v1", biological_n=99, n_features=42163,
                   n_short_features=40, n_regions=96460, receiving_values_read=False, h3_values_read=False,
                   prediction_model_retrained=False, receiving_input_recipe_changed=False, primary_admission_changed=False,
                   development_counts_all_zero_features=sum(r["training_nonzero_participants"] == 0 for r in rows),
                   coefficient_exact_zero_features={arm: int(np.count_nonzero(stats[arm]["nonzero"] == 0)) for arm in stats},
                   short_feature_count_fraction_max=float(mass.max()),
                   checks=dict(ocis_state_coefficients_nonzero=0, direct_denominator_split_max_difference=split_max_error,
                               actual_frozen_scorer_probe_max_difference=scorer_errors),
                   interpretation="Developmental count zeros do not prove biological absence or receiving zeros. Effective coefficients are conditional sensitivities in raw log2 modeled-axis CPM. Perturbations change predictions, not measured H3 error. No receiving bounds follow without justified receiving count/expression bounds; denominator changes affect all modeled genes.",
                   perturbation="Set only forty developmental features to zero, preserve every axis row; separate fixed-denominator direct and changed-denominator terms. This is diagnostic only, not an admitted quantification or zero-fill route.",
                   conditional_rms_bound="For fixed denominator and bounds abs(delta_log2cpm_g)<=b_g, prediction RMS<=sum_g b_g*coefficient_rms_g. No receiving b_g are established.",
                   stochastic_operations=False, input_sha256=hashes, script_sha256=sha(Path(__file__)),
                   launcher_sha256=sha(Path(__file__).with_name("run_audit_short_rna_feature_dependency.sbatch")),
                   environment=dict(python=sys.version, numpy=np.__version__, platform=platform.platform(), slurm_job_id=os.environ["SLURM_JOB_ID"]),
                   output_sha256={p.name: sha(p) for p in args.output.iterdir() if p.is_file()})
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({k: summary[k] for k in ("development_counts_all_zero_features", "coefficient_exact_zero_features", "checks")}, sort_keys=True))


if __name__ == "__main__":
    main()
