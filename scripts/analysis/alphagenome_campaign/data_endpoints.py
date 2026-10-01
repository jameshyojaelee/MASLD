#!/usr/bin/env python3
"""Build assay-native development endpoint tables without choosing targets by outcome.

Run on a compute node. Existing releases and protected evaluation data are read-only.
The wider Currin population is all allele-compatible variant/peak tests at reporter
positions from the existing complete chromosome scans, not a minimum-p target per SNP.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
RES = ROOT / "GWAS/finemapping/results"
BENCH = ROOT / "Analysis/MASLD_Model_Benchmark"
CURRIN = ROOT / "GWAS/finemapping/data/seqfunc_external/currin2025_caqtl_v1"
E1 = RES / "alphagenome_program/e1-bridge-20260914T194145Z"
VAL = BENCH / "executions/gse281364-validated-reconstruction-21066470/validated"
SPLIT = BENCH / "executions/gse281364-outcome-blind-splits-21068885/split"
ASE = RES / "alphagenome_atlas/p3-ase-20260909T192926Z/tables"
FASTA = Path("/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa")


def read(path, **kwargs):
    return pd.read_csv(path, sep="\t", **kwargs)


def write(frame, path):
    if path.exists():
        raise FileExistsError(path)
    frame.to_csv(path, sep="\t", index=False)


def identify(frame, variant_column):
    parts = frame[variant_column].str.split(":", expand=True)
    frame = frame.copy()
    frame["variant_id"] = frame[variant_column]
    frame["assembly"] = "GRCh38"
    frame["chrom"] = parts[0]
    frame["position1"] = parts[1].astype(int)
    frame["ref"] = parts[2].str.upper()
    frame["alt"] = parts[3].str.upper()
    frame["block_1mb"] = frame.chrom + ":" + ((frame.position1 - 1) // 1_000_000).astype(str)
    frame["strand_ambiguous"] = [set((r, a)) in ({"A", "T"}, {"C", "G"}) for r, a in zip(frame.ref, frame.alt)]
    frame["orientation"] = "ALT_minus_REF"
    return frame


def currin_tables(out, fold_map, summary):
    src = CURRIN / "liver_significant_caQTL_leadVariants_1kb_analysis_with_populationAlleleFrequencies.bed.gz"
    lead = identify(read(src), "lead_variant_ID")
    lead["effect"] = lead.beta
    lead["effect_unit"] = "source_normalized_accessibility_slope_per_ALT_dosage"
    lead["source"] = "Currin_2025"
    lead["source_path"] = str(src.relative_to(ROOT))
    lead["endpoint"] = "marginal_caqtl_association"
    lead["inferential_unit"] = "locus; source association estimated across participants"
    lead["target_id"] = lead.peak_ID
    lead["outer_fold"] = lead.chrom.map(fold_map)
    lead["population"] = "source_significant_peak_leads; retain_multiple_targets_per_variant"
    lead["uncertainty_state"] = "SE_not_in_lead_file"
    lead["effect_se"] = np.nan
    nominal_columns = ["tag", "variant_id", "target_id", "distance", "p", "beta_nominal", "varbeta", "imputation", "maf"]
    matched = []
    for chrom in range(1, 23):
        scan = read(E1 / f"scan/chr{chrom}.tsv", header=None, names=nominal_columns)
        matched.append(scan.loc[scan.tag.eq("L"), ["variant_id", "target_id", "beta_nominal", "varbeta"]])
    nominal = pd.concat(matched, ignore_index=True)
    if nominal.duplicated(["variant_id", "target_id"]).any():
        raise ValueError("duplicate source lead association")
    lead = lead.merge(nominal, on=["variant_id", "target_id"], how="left", validate="many_to_one")
    for col in ("beta_nominal", "varbeta"):
        lead[col] = pd.to_numeric(lead[col], errors="coerce")
    overlap = lead.beta_nominal.notna()
    if not np.allclose(lead.loc[overlap, "beta"], lead.loc[overlap, "beta_nominal"], rtol=1e-5, atol=1e-7):
        raise ValueError("lead and nominal source slopes disagree")
    lead["effect_se"] = np.sqrt(lead.varbeta.where(lead.varbeta >= 0))
    lead.loc[lead.effect_se.notna(), "uncertainty_state"] = "sqrt_source_nominal_varbeta_matched_by_variant_and_peak"
    lead["eligible_scalar_snv"] = (lead.ref.str.len() == 1) & (lead.alt.str.len() == 1) & lead.chrom.isin(fold_map)
    lead["exclusion_reason"] = np.where(lead.eligible_scalar_snv, "", "non_autosomal_or_non_SNV")
    lead["source_orientation_consistent"] = ((lead.beta > 0) & (lead.EA == lead.alt)) | ((lead.beta < 0) & (lead.EA == lead.ref))
    if not lead.source_orientation_consistent.all():
        raise ValueError("Currin EA is no longer consistent with accessibility-increasing allele")
    import pysam
    genome = pysam.FastaFile(str(FASTA))
    lead["reference_match"] = [genome.fetch(c, p-1, p-1+len(r)).upper() == r for c,p,r in zip(lead.chrom, lead.position1, lead.ref)]
    lead.loc[~lead.reference_match, "eligible_scalar_snv"] = False
    lead.loc[~lead.reference_match, "exclusion_reason"] = "reference_mismatch"
    write(lead, out / "caqtl_source_leads.tsv.gz")
    summary["caqtl_source_leads"] = {"rows": len(lead), "variants": lead.variant_id.nunique(), "eligible_snv_rows": int(lead.eligible_scalar_snv.sum()), "reference_mismatches": int((~lead.reference_match).sum()), "weak_effect_filter": False}
    return genome


def reporter_tables(out, summary):
    elements = read(SPLIT / "elements.tsv")
    elements = elements.rename(columns={"outer_fold": "source_fold_historical"})
    block = read(BENCH / "executions/gse281364-borzoi-native-fixture-21083008/fixture/source_group_map.tsv")
    elements = elements.merge(block[["source_locus_group_id", "borzoi_long_range_group_id", "outer_fold"]], left_on="outer_locus_sequence_group_id", right_on="source_locus_group_id", validate="many_to_one")
    authority = read(BENCH / "executions/model-check-220-21088696/contract/row_universe.tsv")
    authority["outer_fold"] = authority.outer_fold.str.removeprefix("fold-").astype(int)
    authority = authority[["element_id", "outer_fold", "long_range_block_id"]].drop_duplicates()
    common = elements.merge(authority, on="element_id", suffixes=("", "_authority"), validate="one_to_one")
    if len(common) != 1033 or not ((common.outer_fold == common.outer_fold_authority) & common.borzoi_long_range_group_id.eq(common.long_range_block_id)).all():
        raise ValueError("source-group folds disagree with operative comparator folds")
    if elements.groupby("borzoi_long_range_group_id").outer_fold.nunique().max() != 1:
        raise ValueError("one reporter long-range block crosses folds")
    elements["fold_source"] = "Borzoi_long_range_group_mapping;1033_comparator_elements_verified"
    raw = read(VAL / "replicate_outcomes.tsv.gz", keep_default_na=False)
    raw = raw.loc[raw.element_id.isin(elements.element_id)].copy()
    keys = ["element_id", "context_id", "experimental_replicate"]
    if raw.duplicated(keys + ["allele"]).any():
        raise ValueError("duplicate construct replicate measurement")
    ref = raw.loc[raw.allele == "ref"].drop(columns="allele")
    alt = raw.loc[raw.allele == "alt"].drop(columns="allele")
    paired = ref.merge(alt, on=keys, how="outer", suffixes=("_ref", "_alt"), validate="one_to_one")
    paired = paired.merge(elements, on="element_id", validate="many_to_one")
    paired["eligible"] = paired.assay_state_ref.eq("observed") & paired.assay_state_alt.eq("observed")
    paired["exclusion_reason"] = np.where(paired.eligible, "", "missing_or_source_QC_failed_construct")
    for allele in ("ref", "alt"):
        rna = pd.to_numeric(paired[f"RNA_{allele}"], errors="coerce")
        dna = pd.to_numeric(paired[f"DNA_{allele}"], errors="coerce")
        paired[f"activity_{allele}"] = np.log2((rna + .5) / (dna + .5))
    paired["effect"] = (paired.activity_alt - paired.activity_ref).where(paired.eligible)
    paired["effect_unit"] = "log2_RNA_over_DNA_ALT_minus_REF; pseudocount_0.5"
    paired["orientation"] = "genomic_ALT_minus_REF"
    paired["inferential_unit"] = "experimental_replicate_within_construct; loci_clustered"
    paired["uncertainty_state"] = "retain_four_replicates; no_per_replicate_SE"
    paired["source"] = "GSE281364"
    paired["source_path"] = str((VAL / "replicate_outcomes.tsv.gz").relative_to(ROOT))
    write(paired, out / "reporter_allele_replicates.tsv.gz")
    interactions = []
    for cell, treated in (("LX2", "LX2_TGFb"), ("HepG2", "HepG2_PAOA")):
        keys2 = ["element_id", "experimental_replicate"]
        t = paired.loc[paired.context_id.eq(treated)]
        c = paired.loc[paired.context_id.eq(f"{cell}_control")]
        x = t.merge(c[keys2 + ["effect", "eligible", "sample_id_ref", "sample_id_alt"]], on=keys2, how="outer", suffixes=("_treated", "_control"), validate="one_to_one")
        x["effect"] = x.effect_treated - x.effect_control
        x["eligible"] = x.eligible_treated.fillna(False) & x.eligible_control.fillna(False)
        x["exclusion_reason"] = np.where(x.eligible, "", "missing_or_failed_condition_allele_contrast")
        x["effect_unit"] = "difference_of_log2_RNA_over_DNA_allele_contrasts"
        x["cell_line"] = cell
        x["contrast"] = f"{treated}_minus_{cell}_control"
        x["pairing"] = "same_replicate_index; sample_IDs_retained; shared_control_covariance_not_assumed_zero"
        x["measurement_limit"] = "principal_interaction_comparison" if cell == "LX2" else "measurement_limited_comparison"
        interactions.append(x)
    inter = pd.concat(interactions, ignore_index=True)
    write(inter, out / "reporter_treatment_replicates.tsv.gz")
    summary["reporter"] = {"elements": len(elements), "allele_replicates": len(paired), "eligible_allele_replicates": int(paired.eligible.sum()), "interaction_replicates": len(inter), "independent_patients": 0}
    return elements, paired


def reporter_endogenous(out, elements, paired, genome, summary):
    cols = ["tag", "source_variant", "peak", "distance_from_peakCenter", "pvalue", "beta", "varbeta", "imputation_R2", "MAF"]
    chunks = []
    for chrom in range(1, 23):
        d = read(E1 / f"scan/chr{chrom}.tsv", header=None, names=cols)
        chunks.append(d.loc[d.tag.eq("M")].copy())
    d = pd.concat(chunks, ignore_index=True)
    d = identify(d, "source_variant")
    d = d.merge(elements, left_on=["chrom", "position1"], right_on=["contig", "variant_pos1"], how="left", validate="many_to_many")
    direct = d.ref.eq(d.genomic_ref) & d.alt.eq(d.genomic_alt)
    swapped = d.alt.eq(d.genomic_ref) & d.ref.eq(d.genomic_alt)
    d["orientation_multiplier"] = np.where(direct, 1, np.where(swapped, -1, 0))
    d["eligible"] = direct | swapped
    d["exclusion_reason"] = np.where(d.eligible, "", "alleles_do_not_match; complement_not_inferred")
    for col in ("beta", "varbeta", "distance_from_peakCenter", "pvalue"):
        d[col] = pd.to_numeric(d[col], errors="coerce")
    d["effect"] = (d.beta * d.orientation_multiplier).where(d.eligible)
    d["effect_se"] = np.sqrt(d.varbeta.where(d.varbeta >= 0))
    d["effect_unit"] = "source_normalized_accessibility_slope_per_reporter_ALT_dosage"
    d["source"] = "Currin_2025_nominal_at_GSE281364_variant_positions"
    d["source_path"] = str((E1 / "scan").relative_to(ROOT))
    d["eligibility_definition"] = "all_nominal_tests_at_reporter_positions; no_pvalue_or_effect_filter; no_target_selection"
    peaks = read(CURRIN / "supplementalData1_liver_ATAC_peaks.bed.gz").rename(columns={"#chr":"peak_chrom", "start":"peak_start0", "end":"peak_end0"})
    # Verify the peak-ID join against an independent file carrying explicit hg38 columns.
    check = read(CURRIN / "liver_significant_caQTL_leadVariants_1kb_analysis_with_populationAlleleFrequencies.bed.gz")
    check = check.merge(peaks, left_on="peak_ID", right_on="peakID", validate="many_to_one")
    if not ((check.peak_start_hg38 == check.peak_start0) & (check.peak_stop_hg38 == check.peak_end0)).all():
        raise ValueError("deposited peak coordinates disagree with explicit hg38 lead coordinates")
    d = d.merge(peaks[["peakID", "peak_chrom", "peak_start0", "peak_end0", "roadmap_liverTissue_chromatinState"]], left_on="peak", right_on="peakID", how="left", validate="many_to_one")
    known = d.peak_start0.notna() & d.peak_end0.notna() & d.chrom.eq(d.peak_chrom)
    d["physical_variant_peak_overlap"] = ((d.position1-1 >= d.peak_start0) & (d.position1-1 < d.peak_end0)).where(known)
    d["physical_construct_peak_overlap_bp"] = np.maximum(0, np.minimum(d.reference_match_end0, d.peak_end0) - np.maximum(d.reference_match_start0, d.peak_start0)).where(known)
    d["peak_coordinate_state"] = np.where(known, "source_contract_GRCh38_BED0_half_open", "missing_peak_coordinate")
    distance = d.distance_from_peakCenter.abs()
    d["distance_stratum"] = pd.cut(distance, [-1,500,1000,5000,10000,100000,np.inf], labels=["0-500", "501-1000", "1001-5000", "5001-10000", "10001-100000", ">100000"]).astype(str)
    d["target_link_state"] = "not_supplied_in_Currin_variant_peak_association"
    # Activity is source-native and frozen before new transfer evaluation. Its use as a
    # model feature or selected stratum must still be trained/selected inside folds.
    activity = paired.loc[paired.eligible & paired.context_id.eq("LX2_control")].groupby("element_id")[["activity_ref", "activity_alt"]].mean().mean(axis=1)
    d["LX2_control_mean_log2_activity"] = d.element_id.map(activity)
    d["activity_stratum"] = np.where(d.LX2_control_mean_log2_activity.isna(), "missing", np.where(d.LX2_control_mean_log2_activity > 0, "RNA_over_DNA_gt1", "RNA_over_DNA_le1"))
    write(d, out / "reporter_caqtl_all_nominal_pairs.tsv.gz")
    eligible = d.loc[d.eligible]
    census = []
    for name, mask in (("all_allele_compatible", np.ones(len(eligible), dtype=bool)), ("variant_physically_inside_peak", eligible.physical_variant_peak_overlap.eq(True)), ("construct_overlaps_peak", eligible.physical_construct_peak_overlap_bp.gt(0)), ("nearest_center_within_500bp", eligible.distance_from_peakCenter.abs().le(500))):
        sub = eligible.loc[mask]
        census.append({"stratum":name,"pairs":len(sub),"variants":sub.element_id.nunique(),"blocks":sub.borzoi_long_range_group_id.nunique(),"peaks":sub.peak.nunique()})
    write(pd.DataFrame(census), out / "reporter_physical_overlap_census.tsv")
    summary["reporter_endogenous"] = {"all_position_match_rows":len(d),"allele_compatible_rows":len(eligible),"missing_peak_coordinates":int((~known).sum()),"census":census,"transfer_fitted":False}


def allelic_tables(out, summary):
    counts = []
    for cohort, fname in (("GSE281367","allelic_donor_counts.tsv.gz"),("GSE244832","gse244832_allelic_donor_counts.tsv.gz")):
        src = ASE / fname
        d = identify(read(src), "uid")
        d["source"] = cohort
        d["source_path"] = str(src.relative_to(ROOT))
        d["donor_id"] = cohort + ":" + d.donor.astype(str)
        d["allele_informative_reads"] = d.n_ref + d.n_alt
        d["eligible"] = d.is_het.astype(str).str.lower().eq("true") & d.allele_informative_reads.gt(0)
        d["exclusion_reason"] = np.where(d.eligible,"","source_read_based_heterozygosity_or_no_allele_reads")
        d["effect_unit"] = "ALT_reads_conditional_on_ALT_plus_REF_reads"
        d["inferential_unit"] = "donor_crossed_with_locus"
        d["genotype_state"] = "read_based_heterozygote; independent_genotype_not_joined"
        d["phase_state"] = "unknown"
        d["mapping_bias_state"] = "source_counts; use_separate_remapping_sensitivity_for_claims"
        counts.append(d)
    d = pd.concat(counts, ignore_index=True)
    write(d, out / "allelic_accessibility_donor_counts.tsv.gz")
    summary["allelic"] = {"rows":len(d),"eligible_rows":int(d.eligible.sum()),"donors":d.donor_id.nunique(),"sites":d.variant_id.nunique(),"phase_known":False}


def splice_tables(out, summary):
    src = ROOT / "data/external/gtex_v8_liver_sqtl/GTEx_Analysis_v8_sQTL/Liver.v8.sqtl_signifpairs.txt.gz"
    d = read(src)
    # Entire source-native event family retained. No model score defines event eligibility.
    d["source_variant_id"] = d.variant_id
    d["variant_id"] = d.variant_id.str.replace("_b38", "", regex=False).str.replace("_", ":", regex=False)
    d = identify(d,"variant_id")
    e = d.phenotype_id.str.split(":",expand=True)
    d["event_chrom"],d["intron_start_source"],d["intron_end_source"],d["cluster_id"],d["gene_id_versioned"] = [e[i] for i in range(5)]
    d["source"] = "GTEx_v8_Liver_sQTL"
    d["source_path"] = str(src.relative_to(ROOT))
    d["effect"] = d.slope
    d["effect_se"] = d.slope_se
    d["effect_unit"] = "source_normalized_intron_excision_phenotype_per_ALT_dosage; not_delta_PSI"
    d["strand"] = "unresolved_in_signifpairs"
    d["event_identity_state"] = "exact_source_cluster_and_boundaries; strand_and_full_competing_junction_set_required"
    d["new_event_evaluation_eligible"] = False
    d["exclusion_reason"] = "requires_full_source_cluster_membership_strand_and_REF_ALT_junction_union"
    d["pretrained_exposure"] = "documented_GTEx_v8_training_overlap"
    d["population"] = "source_significant_pairs; not_all_tested_events"
    write(d,out / "splicing_source_events.tsv.gz")
    summary["splicing"] = {"pairs":len(d),"clusters":d.cluster_id.nunique(),"new_event_comparison_run":False,"missing":"strand; full cluster denominator; raw REF/ALT junction outputs"}


def donor_table(out, summary):
    src = BENCH / "executions/gse267145-authoritative-join-21064930/participant_join.tsv"
    d = read(src)
    if len(d) != 99 or d.participant_id.nunique() != 99:
        raise ValueError("GSE267145 participant contract changed")
    d["rna_measurement_state"] = "nonnegative_fractional_expression_estimates"
    d["h3k27ac_measurement_state"] = "deposited_integer_region_counts"
    d["fragment_or_read_count_convention"] = "unresolved_by_inspected_metadata"
    d["coordinate_semantics"] = "counted_interval; inferred_1_based_inclusive"
    d["coordinate_uncertainty_bp"] = 1
    d["sequence_use"] = "permitted_by_2026-09-15_owner_decision"
    d["called_peak_base_correspondence"] = "unresolved"
    d["genotype_state"] = "not_joined"
    d["source_path"] = str(src.relative_to(ROOT))
    d["redistribution"] = "local_internal_analysis_only; no_source_data_rehosting"
    write(d,out / "donor_chromatin_participants.tsv")
    summary["donor_chromatin"] = {"participants":99,"rna_features_source":43285,"rna_features_admitted":42163,"retired_rna_features_masked":1122,"H3K27ac_regions":96460,"new_matrix_copy":False,"primary_error":"logCPM; calibration_fit_inside_training_donors_and_regions"}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out",required=True,type=Path)
    args=p.parse_args()
    args.out.mkdir(parents=True,exist_ok=False)
    f=read(RES / "seqfunc/adult_liver_chrombpnet/hepatocyte_5fold_v2/fold_manifest.tsv")
    folds={c:int(r.fold) for r in f.itertuples() for c in r.test_chromosomes.split(",")}
    summary={"protected_outcomes_read":False,"source_data_release":False,"selection":"no_new_pvalue_or_effect_filters","statistical_tests":"not_run; endpoint_construction_only"}
    genome=currin_tables(args.out,folds,summary)
    elements,paired=reporter_tables(args.out,summary)
    reporter_endogenous(args.out,elements,paired,genome,summary)
    allelic_tables(args.out,summary)
    splice_tables(args.out,summary)
    donor_table(args.out,summary)
    genome.close()
    (args.out / "summary.json").write_text(json.dumps(summary,indent=2,default=int)+"\n")
    print(json.dumps(summary,indent=2,default=int),flush=True)


if __name__ == "__main__":
    main()
