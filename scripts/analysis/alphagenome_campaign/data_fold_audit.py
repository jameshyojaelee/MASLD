#!/usr/bin/env python3
"""Resolve operative reporter folds and tabulate fixed physical-overlap strata.

No transfer predictions or effects are evaluated. Do not use the historical
outer_fold in endpoints-21771964 reporter files: join this companion by element_id.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[3]
BENCH=ROOT / "Analysis/MASLD_Model_Benchmark/executions"


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--endpoints",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    args.out.mkdir(parents=True,exist_ok=False)
    historical=pd.read_csv(BENCH / "gse281364-outcome-blind-splits-21068885/split/elements.tsv",sep="\t")
    historical=historical.rename(columns={"outer_fold":"source_fold_historical"})
    groups=pd.read_csv(BENCH / "gse281364-borzoi-native-fixture-21083008/fixture/source_group_map.tsv",sep="\t")
    groups=groups.rename(columns={"outer_fold":"outer_fold","borzoi_long_range_group_id":"long_range_block_id"})
    mapping=historical[["element_id","outer_locus_sequence_group_id","source_fold_historical"]].merge(groups[["source_locus_group_id","outer_fold","long_range_block_id"]],left_on="outer_locus_sequence_group_id",right_on="source_locus_group_id",validate="many_to_one")
    if mapping.groupby("long_range_block_id").outer_fold.nunique().max()!=1:
        raise ValueError("one long-range sequence block spans multiple folds")
    authority=pd.read_csv(BENCH / "model-check-220-21088696/contract/row_universe.tsv",sep="\t")
    authority["outer_fold"]=authority.outer_fold.str.removeprefix("fold-").astype(int)
    authority=authority[["element_id","outer_fold","long_range_block_id"]].drop_duplicates()
    if authority.element_id.duplicated().any() or len(authority)!=1033:
        raise ValueError("row-authority fold/block is not unique per1033 elements")
    matched=mapping.merge(authority,on="element_id",suffixes=("","_authority"),validate="one_to_one")
    if len(matched)!=1033 or not ((matched.outer_fold==matched.outer_fold_authority)&(matched.long_range_block_id==matched.long_range_block_id_authority)).all():
        raise ValueError("Borzoi group folds disagree with operative comparator authority")
    mapping["fold_source"]="gse281364-borzoi-native-fixture-21083008/fixture/source_group_map.tsv;1033_comparator_rows_verified_against_model-check-220-21088696"
    mapping["all_assays_conditions_same_split"]=True
    mapping.to_csv(args.out / "reporter_fold_companion.tsv",sep="\t",index=False)
    nominal=pd.read_csv(args.endpoints / "reporter_caqtl_all_nominal_pairs.tsv.gz",sep="\t")
    nominal=nominal.rename(columns={"outer_fold":"source_fold_historical"}).drop(columns=["borzoi_long_range_group_id"])
    nominal=nominal.merge(mapping[["element_id","outer_fold","long_range_block_id"]],on="element_id",validate="many_to_one")
    valid=nominal.eligible.astype(str).str.lower().eq("true")
    strata={"all_allele_compatible":valid,"variant_inside_peak":valid & nominal.physical_variant_peak_overlap.astype(str).str.lower().eq("true"),"construct_overlaps_peak":valid & nominal.physical_construct_peak_overlap_bp.gt(0),"center_within_500bp":valid & nominal.distance_from_peakCenter.abs().le(500)}
    census=[]
    for name,mask in strata.items():
        for fold in range(5):
            d=nominal.loc[mask & nominal.outer_fold.eq(fold)]
            census.append({"stratum":name,"outer_fold":fold,"variant_peak_pairs":len(d),"reporter_elements":d.element_id.nunique(),"independent_sequence_blocks":d.long_range_block_id.nunique(),"peak_targets":d.peak.nunique(),"context":"bulk_liver_caqtl_vs_reporter_cell_line;assay_context_mismatch_retained","orientation":"native_Currin_beta_oriented_to_reporter_genomic_ALT"})
    pd.DataFrame(census).to_csv(args.out / "transfer_strata_by_operative_fold.tsv",sep="\t",index=False)
    summary={"elements":len(mapping),"long_range_blocks":mapping.long_range_block_id.nunique(),"authority_rows_verified":len(matched),"historical_fold_disagreements_on_1033":int((matched.source_fold_historical!=matched.outer_fold).sum()),"block_fold_conflicts":0,"transfer_fit_run":False,"precision_assessment":"block_counts_only; native_margin_and_attainable_CI_require_training-only_pilot_before_protected_evaluation","minimum_for_computable_fold_uncertainty":"at_least_two_blocks_per_fold_is_necessary_not_sufficient; no_numeric_power_claim"}
    (args.out / "summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary,indent=2),flush=True)


if __name__=="__main__":
    main()
