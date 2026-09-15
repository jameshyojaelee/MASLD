#!/usr/bin/env python3
"""Recheck existing Figure 7 predictions and retain source-specific limitations.

No training, model selection, or canonical output writes. Input/output hashes
make the reproduced population and the precise scorer version recoverable.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys

import numpy as np
import pandas as pd
import scipy
from scipy.stats import rankdata, spearmanr


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    root, out = args.root.resolve(), args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    bench = root / "Analysis/MASLD_Model_Benchmark"
    inputs, checks, summaries, memberships = set(), [], [], []

    def track(path):
        path = Path(path).resolve(strict=True)
        inputs.add(path)
        return path

    def table(path):
        return pd.read_csv(track(path), sep="\t")

    def js(path):
        return json.loads(track(path).read_text())

    def check(name, passed, detail):
        checks.append(dict(check=name, passed=bool(passed), detail=str(detail)))
        if not passed:
            print(f"CHECK FAILED: {name}: {detail}", flush=True)

    def require_unique(df, key, name):
        ok = df[key].notna().all() and not df[key].duplicated().any()
        check(name, ok, len(df))
        if not ok:
            raise ValueError(f"Ambiguous identity: {name}")

    def metric(scores, y):
        return float(spearmanr(scores, y).statistic)

    def ceiling(y):
        return float(np.corrcoef(np.sort(rankdata(y)), np.arange(1, len(y) + 1))[0, 1])

    def record_members(cohort, analysis, ids, y):
        memberships.extend(dict(cohort=cohort, analysis=analysis, participant_id=str(i), outcome=float(v)) for i, v in zip(ids, y))

    # Reference score: use the frozen released scorer, with explicit input axes.
    release = bench / "release/masld-severity-v1.2"
    ood = bench / "executions/full-axis-rebuild-and-ood-transfer-20260901T122932Z"
    expected = js(release / "expected_ood.json")
    counts = np.load(track(ood / "gse268273_model_axis_expected_counts.npy"))
    genes = table(ood / "gse268273_model_axis_genes.tsv")["stable_gene_id"].astype(str)
    people = table(ood / "gse268273_participants.tsv")
    require_unique(people, "row_id", "reference_count_axis_unique")
    require_unique(pd.DataFrame({"gene": genes}), "gene", "reference_gene_axis_unique")
    assert counts.shape == (len(people), len(genes))
    np.savez_compressed(out / "reference_counts.npz", counts=counts.T,
                        genes=genes.to_numpy(str), samples=people.row_id.to_numpy(str))
    for name in ("score.py", "kleiner.py", "MODEL_CARD.md"):
        track(release / name)
    for name in ("masld_severity_v1_weights.npz", "activity_head_v1.npz", "shift_diagnostic_v1.npz"):
        track(release / "weights" / name)
    cmd = [sys.executable, str(release / "score.py"), "--counts", str(out / "reference_counts.npz"),
           "--out", str(out / "reference_scores.tsv"), "--json", str(out / "reference_score_report.json"),
           "--weights", str(release / "weights/masld_severity_v1_weights.npz")]
    pr = subprocess.run(cmd, capture_output=True, text=True, check=False)
    (out / "scorer_stdout.txt").write_text(pr.stdout)
    (out / "scorer_stderr.txt").write_text(pr.stderr)
    (out / "scorer_command.json").write_text(json.dumps(cmd, indent=2) + "\n")
    if pr.returncode:
        raise RuntimeError(f"Released scorer refused input: {pr.stderr}")
    s = pd.read_csv(out / "reference_scores.tsv", sep="\t")
    ph = table(bench / "executions/gse268273-evaluator-phenotype-20260830T205701Z/gse268273_participant_phenotype.tsv")
    require_unique(ph, "row_id", "reference_phenotype_unique")
    require_unique(s, "sample_id", "reference_prediction_unique")
    m = s.merge(ph[["row_id", "fibrosis_stage"]], left_on="sample_id", right_on="row_id", validate="one_to_one")
    assert len(m) == len(s) == expected["n_participants"]
    rho = metric(m.latent_severity, m.fibrosis_stage)
    check("reference_spearman_reproduced", abs(rho - expected["ood_spearman"]) <= expected["tolerance"], rho)
    record_members("GSE268273", "released_severity", m.row_id, m.fibrosis_stage)
    summaries.append(dict(task="released_severity", cohort="GSE268273", n=len(m), metric="Spearman", value=rho,
      baseline_value=0.609434, verification="Frozen scorer and participant join reproduced; baseline not refitted",
      training_unit="521 sample records; donor identities unresolved in model card", baseline="In-cohort ridge, as recorded",
      exposure="External to this severity fit; 109/109 also occur in BulkFormer pretraining, a different model",
      uncertainty="Existing model card; no new retraining uncertainty estimate", source=str(release / "MODEL_CARD.md")))

    ext = bench / "executions/external-validation-expansion-20260901T190000Z"
    old = js(ext / "results.json")
    track(ext / "external_validation.py")
    for cohort in ("GSE213621", "GSE276114"):
        s = table(ext / f"{cohort}_scores.tsv")
        require_unique(s, "sample_id", f"{cohort}_scores_unique")
        if cohort == "GSE213621":
            # The deposited table has a repeated leading run field. pandas
            # represents that extra field as the index, matching the source reader.
            ph = table(root / "Analysis/Deconvolution/bulk/GSE213621/GSE213621_metadata.tsv")
            ph = ph.rename(columns={"sample_id": "sample", "fibrotic_stage": "stage_raw"})
            ph["participant_id"] = ph["sample"].astype(str)
            mappings = {"primary_3level_no_control": {"F0F1": 0, "F2": 1, "F3F4": 2},
                        "sensitivity_4level_with_control": {"Control": 0, "F0F1": 1, "F2": 2, "F3F4": 3}}
        else:
            ph = table(bench / "executions/etiology-transfer-fibrosis-distinctness-20260901T135543Z/gse276114_design.tsv")
            ph = ph.rename(columns={"matrix_column_name": "sample", "disease_group": "stage_raw"})
            ph["participant_id"] = ph.geo_accession.astype(str)
            mappings = {"primary_3level": {"F0-2": 0, "F3": 1, "F4": 2}}
        require_unique(ph, "sample", f"{cohort}_phenotype_join_key_unique")
        m = s.merge(ph[["sample", "stage_raw", "participant_id"]], left_on="sample_id", right_on="sample", validate="one_to_one")
        require_unique(m, "participant_id", f"{cohort}_evaluation_unit_unique")
        assert len(m) == len(s)
        for arm, mapping in mappings.items():
            y = m.stage_raw.astype(str).map(mapping)
            ok = y.notna()
            yy = y[ok].to_numpy(float)
            rho, ceil = metric(m.loc[ok, "latent_severity"], yy), ceiling(yy)
            rr = old["results"][f"{cohort}__{arm}"]
            check(f"{cohort}_{arm}_n", len(yy) == rr["n_participants"], len(yy))
            check(f"{cohort}_{arm}_rho", abs(rho - rr["transfer_spearman"]) < 2e-6, rho)
            check(f"{cohort}_{arm}_ceiling_ratio", abs(rho / ceil - rr["spearman_over_ceiling"]) < 2e-6, rho / ceil)
            record_members(cohort, arm, m.loc[ok, "participant_id"], yy)
            summaries.append(dict(task=arm, cohort=cohort, n=len(yy), metric="Spearman", value=rho, tie_ceiling=ceil,
              value_over_ceiling=rho / ceil, baseline_value=rr["native_in_cohort_spearman_ridge"],
              paired_delta_recorded=rr["vs_native_paired_delta"], paired_ci_recorded=json.dumps(rr["vs_native_paired_delta_ci95"]),
              verification="Transfer score and participant join reproduced; native individual predictions not saved by source script",
              baseline="Transductive in-cohort ridge: full-cohort expression scaling; outer-training HVGs/PCA; inner preprocessing not refitted",
              training_unit="Severity fit: 521 samples, unresolved donor identity; baseline evaluation: declared one sample per participant",
              exposure="External to severity fit; GSE213621 is also a Resource bulk cohort; no paper-wide untouched claim",
              uncertainty="Source: 10000 stage-stratified paired participant bootstrap draws, seed 213621; conditional on predictions; not independently rerun",
              source=str(ext / "results.json")))

    # Verify biological pairing and fold membership for the existing chromatin task.
    fixture = bench / "executions/model-data-064-21079902/fixture"
    part = table(fixture / "molecular/participant_axis.tsv")
    folds = table(fixture / "folds/participant_outer_folds.tsv")
    require_unique(part, "participant_id", "chromatin_participant_axis_unique")
    require_unique(folds, "participant_id", "chromatin_fold_assignment_unique")
    paired = part.merge(folds, on="participant_id", validate="one_to_one")
    check("chromatin_pairing_and_folds", len(paired) == 99 and paired.outer_fold.nunique() == 5, len(paired))
    paired.to_csv(out / "chromatin_participant_folds.tsv", sep="\t", index=False)
    croot = bench / "executions/chromatin-state-thesis-20260906T124724Z"
    cr = js(croot / "out/a2_results.json")
    track(croot / "run_a2.py")
    for name in ("F2_forward_RNA_to_H3K27ac", "R1_reverse_H3K27ac_to_RNA"):
        rr = cr["arms"][name]
        summaries.append(dict(task=name, cohort="GSE267145", n=99, metric="log2_CPM_squared_error_skill", value=rr["skill"],
          verification="Pairing and fold identifiers rechecked; aggregate result read; per-participant losses not saved in this output",
          baseline="Training-fold target mean; rank-4 reduced-rank regression", training_unit="Participants, five fixed folds",
          uncertainty=json.dumps(rr["ci95"]) + "; source participant bootstrap, not independently reproduced",
          exposure="Development-exposed single cohort; age and BMI absent; rank 4 is a floor, not optimal rank", source=str(croot / "out/a2_results.json")))
    # Preserve a complete index of other Figure 7 claims without inventing validation.
    for task, source, limit in [
        ("repeat_biopsy_ICC", release / "MODEL_CARD.md", "58 paired participants; source card updated 2026-09-08: only 3 pairs have both biopsies within shipped acceptance bounds; ICC requires explicit out-of-contract scoring; possible true change and shared technical effects; not recalculated here"),
        ("activity_head", release / "MODEL_CARD.md", "Separate endpoint and evaluation populations; does not establish independent biological axes"),
        ("gene_absence_sensitivity", bench / "executions/rerun-floor-20260905T172141Z/SESSION_FINDINGS.md", "Model-specific coverage experiment; not recomputed"),
        ("TabICLv2_vs_ridge", bench / "MODELS.md", "Task-specific recorded comparisons; no general superiority claim; not recomputed"),
        ("RNA_H3K27ac_late_fusion", croot / "out/b_results.json", "Paired bootstrap interval includes potentially meaningful gain; not evidence of equivalence"),
        ("DNA_signed_effect_direction", bench / "MODEL_RESULTS_SUMMARY.md", "Variant/locus task, not patient n; source-reported allele baseline; not recomputed"),
        ("refusal_bounds", release / "MODEL_CARD.md", "Evaluation samples informed some acceptance boundaries; separate score evaluation from boundary calibration"),
    ]:
        track(source)
        summaries.append(dict(task=task, verification="Existing record inspected; no independent numerical reproduction", source=str(source), exposure=limit))
    track(bench / "executions/rerun-floor-20260905T172141Z/EXPOSURE_MEASURED.md")
    pd.DataFrame(summaries).to_csv(out / "model_evaluation_summary.tsv", sep="\t", index=False, na_rep="NA")
    pd.DataFrame(memberships).to_csv(out / "evaluation_membership.tsv", sep="\t", index=False)
    pd.DataFrame(checks).to_csv(out / "model_checks.tsv", sep="\t", index=False)
    hashes = []
    for p in sorted(inputs):
        h = hashlib.sha256()
        with p.open("rb") as fh:
            for block in iter(lambda: fh.read(1 << 20), b""):
                h.update(block)
        hashes.append(dict(source=str(p), sha256=h.hexdigest()))
    pd.DataFrame(hashes).to_csv(out / "input_hashes.tsv", sep="\t", index=False)
    (out / "environment.json").write_text(json.dumps(dict(python=sys.version, executable=sys.executable,
        numpy=np.__version__, pandas=pd.__version__, scipy=scipy.__version__, hostname=platform.node(),
        slurm_job_id=os.getenv("SLURM_JOB_ID"), new_bootstrap_draws=0,
        reason="Existing uncertainty retained; absent paired prediction arrays prevent reproducing paired intervals without refitting"), indent=2) + "\n")
    print(f"Model checks: {sum(c['passed'] for c in checks)}/{len(checks)}; summary rows: {len(summaries)}", flush=True)
    if not all(c["passed"] for c in checks):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
