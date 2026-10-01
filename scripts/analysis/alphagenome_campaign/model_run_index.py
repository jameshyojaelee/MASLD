#!/usr/bin/env python3
"""Write a model-work checkpoint from actual Slurm accounting and result files."""
import argparse
import csv
from datetime import datetime,timezone
import io
import json
from pathlib import Path
import re
import subprocess

RUNS={
    21771783:("full_native_1Mb","native","Full existing32322-lead comparison; existing significance-selected/deduplicated population. Full attempt accounting required before performance."),
    21772093:("scalar_100step_pilot","pilot_21772093","100steps and invariants passed; initial native-function drift included padded channels and is nonfinite, unusable; preserved."),
    21772542:("first_sweep_probe","sweep_21772542","100step probe passed; broader12arm sweep deferred by conservative per-example evaluation projection. No12arm result from this job."),
    21772543:("fixed_feature_extraction","frozen","Both2048/16384 lengths and variant/symmetric/whole measured-target pools, explicit exclusions."),
    21772545:("reporter_wrong_environment","lx2_21772545","Failed before fits: sklearn absent from alphagenome_local; no scientific result."),
    21772726:("bounded_adaptation_sweep","sweep_21772726","Eight adapters,two partial q/v,two fixed-head controls; all12 completed5000steps, one seed/validation fold, no formal finalist claim."),
    21772727:("reporter_historical1033","lx2_21772727","Completed1033 constructs,239source blocks,4 experimental replicate pairs percell; all fitted point errors exceed zero baseline."),
    21772729:("fixed_12_recipes","frozen/comparisons","Completed12recipes with0failedfolds on28600commonrows/2628one-Mb bins; fixed single-seed fivefold development."),
    21772919:("final_matched_uncertainty","final_21772919","Dependency on full native accounting; paired one-Mb-bin intervals plus chromosome-resampling sensitivity; exact specialist target intersection."),
    21772928:("full_reporter_target_features","reporter_targets_21772928","4358/4359 two-kb and4349/4359 sixteen-kb reporters;2765whole-target pairs; zero cross-fold overlapping16kbwindows."),
    21773114:("transfer_with_null_selection","transfer_21773114","Both methods selected zero in all5inner selections; degenerate identical predictions. Original precision-exclusion labels do not establish an information limit."),
    21773115:("reporter_full_population","lx2_full_21773115","4358constructs/239blocks/4replicatepairs; all fitted point errors exceed zero. Common1033 sequences match archive but numerical features differ; not a population-only comparison."),
    21773124:("bounded_sweep_summary","sweep_summary_21773124","14arms including allele and zero controls,6834validationvariants/524one-Mb bins; no fold0 evaluation, no native comparison while population incomplete."),
    21773128:("transfer_learned_diagnostic","transfer_r2_21773128","Development amendment retains null selection and compares learned-only candidates from same grid/folds/strata. All9 estimable CI include zero;12plannedstrata BH family includes3untestable."),
    21773146:("zero_adapter_output_invariant","zero_outputs_21773146","Output-only32sequence check of actual compiled fixed/zero-adapter routes, repeatfloor, perallele and scalar effects; no training/outcomes.")}


def main(args):
    if args.out.exists():raise ValueError("Use a fresh checkpoint filename; preserve prior run indexes")
    raw=subprocess.check_output(["sacct","-j",",".join(map(str,RUNS)),"--format=JobID,State,ElapsedRaw,AllocCPUS,ReqTRES,MaxRSS","-P"],text=True)
    accounting=list(csv.DictReader(io.StringIO(raw),delimiter="|"));by_id={r["JobID"]:r for r in accounting}
    records=[];gpu_hours=0.;cpu_hours=0.
    for job,(name,path,limitation) in RUNS.items():
        row=by_id.get(str(job),{})
        seconds=int(row.get("ElapsedRaw",0));cpus=int(row.get("AllocCPUS",0))
        match=re.search(r"(?:^|,)gres/gpu=(\d+)",row.get("ReqTRES",""));gpus=int(match.group(1)) if match else 0
        cpu=seconds*cpus/3600;gpu=seconds*gpus/3600;cpu_hours+=cpu;gpu_hours+=gpu
        entry={"job_id":job,"assignment":name,"slurm_state":row.get("State","unavailable"),"elapsed_seconds_at_snapshot":seconds,
            "allocated_cpu_core_hours_at_snapshot":cpu,"allocated_gpu_hours_at_snapshot":gpu,
            "max_RSS_accounting_raw":by_id.get(str(job)+".batch",{}).get("MaxRSS") or None,
            "output_directory":str(args.model/path),"scientific_disposition_and_limits":limitation}
        for filename in ("complete.json","completion.json","analysis.json","disposition.json"):
            candidate=args.model/path/filename
            if candidate.exists():entry[filename]=json.loads(candidate.read_text())
        records.append(entry)
    output={"recorded_UTC":datetime.now(timezone.utc).isoformat(),"scope":"model_worker_jobs_only; coordinator_budget_ledger_is_campaign_authority",
        "status_source":"live_sacct_at_recorded_time; running_and_pending_values_are_not_terminal_costs",
        "allocated_gpu_hours_at_snapshot":gpu_hours,"allocated_CPU_core_hours_at_snapshot":cpu_hours,"runs":records,
        "unchanged_boundaries":["protected_outcomes_closed","fixed117_unchanged","genetic_membership_unchanged","Cas13_separate","no_weight_publication_or_final_adoption"],
        "unrun_or_unresolved":["five_seed_nested_finalists","full_backbone_adaptation","Gaussian_SE_vs_MSE_training_comparison",
            "native_function_accuracy_with_measured_labels; drift_only_measured","numeric_cause_of_archived_vs_new_reporter_feature_difference",
            "compiled_zero_adapter_output_agreement_pending_not_run_at_this_checkpoint",
            "full_native_results_and_final_uncertainty_until_dependency_finishes"],
        "common_reporter_feature_difference":{"shared_variants":1033,"sequence_identity_mismatches":0,
            "allele_feature_relative_L2":0.004002501785940633,"signed_effect_feature_relative_L2":0.030797127484605582,
            "signed_effect_feature_global_correlation":0.9995256607634009,
            "cause":"unresolved; same_checkpoint_path_input_sequences_pooling_bins_float32_inputs_inference_helper_and_L40S",
            "interpretation":"not_a_clean_population_only_reporter_comparison"},
        "source_effect_inference":"source_QTL_beta_is_participant_based_association; variant_counts_are_not_biological_donor_n; per_variant_effective_n_not_established",
        "source_intervals":"one_Mb_bins_do_not_establish_LD_or_window_independence; chromosome_sensitivity_is_separate",
        "checks":"scientific invariant receipts in inspection/ and per-pilot pretraining_checks/feasibility; output-route invariant receipt kept in its own job directory"}
    args.out.write_text(json.dumps(output,indent=2)+"\n")
    print(json.dumps({"runs":len(records),"gpu_hours_at_snapshot":gpu_hours,"CPU_core_hours_at_snapshot":cpu_hours,"out":str(args.out)}))


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    main(parser.parse_args())
