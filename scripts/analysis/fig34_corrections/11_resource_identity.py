"""Resource-only identity screens and source-qualified analysis rosters.

Run on an allocated compute node. A1 and Model A inputs are read-only. The old
other-allele fraction and expression-sex markers are retained as diagnostics;
neither is a mixture diagnosis or an exclusion rule here. Genetic identity
screens cannot distinguish duplicate samples from monozygotic twins.
"""
import argparse
import hashlib
import itertools
import json
import os
import platform
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
A1 = ROOT / "Analysis/MASLD_Model_Benchmark/executions/in-biopsy-allelic-a1-full-20260924T025858Z/identity"
QC = ROOT / "Analysis/MASLD_Model_Benchmark/executions/in-biopsy-allelic-obs-20260926T155651Z/qc/v2-20260929T143503Z"
VCF = ROOT / "RNA-seq/Human/Patient_Cohorts/results/rna_genotypes_restricted/a1-full-20260924T025858Z/vcf/identity_panel.merged.vcf.gz"
METADATA = ROOT / "RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv"
HISTOLOGY = ROOT / "figures/candidates/fig34-scientific-corrections-20260929-v1/continuum/participant_histology.tsv"
BULK = ROOT / "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20260923T231939Z/primary_source_controls/model_design.tsv"
SOURCE130 = ROOT / "RNA-seq/Human/Patient_Cohorts/metadata/source/GSE130970/GSE130970_source_diagnosis.tsv"
SOURCE193 = ROOT / "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/workstreams/BULK-STAGE/PREFLIGHT-v1/audits/gse193066_participant_crosswalk.tsv"
JOIN267 = ROOT / "Analysis/MASLD_Model_Benchmark/executions/gse267145-authoritative-join-21064930/participant_join.tsv"
ENDPOINT = ROOT / "figures/candidates/bulk-release-corrected-controls-20260925T175215Z/source_tables/fig3d_endpoint_pca_coordinates.tsv"
STAGE = ROOT / "figures/candidates/bulk-release-corrected-controls-20260925T175215Z/analysis/stage_extensions/stage_extension_sample_manifest.tsv"
MIN_SITES, MIN_SHARED, MIN_HALF_SHARED = 1000, 500, 250
SAME_KINSHIP = 0.354


def read(path):
    return pd.read_csv(path, sep="\t", keep_default_na=False, na_values=["", "NA"])


def yes(values):
    return values.astype(str).str.lower().isin(["true", "1"])


def write(frame, path):
    if path.exists():
        raise FileExistsError(path)
    frame.to_csv(path, sep="\t", index=False, na_rep="NA")


def key(a, b):
    return tuple(sorted((str(a), str(b))))


def components(samples, edges):
    parent = {s: s for s in samples}

    def find(s):
        while parent[s] != s:
            parent[s] = parent[parent[s]]
            s = parent[s]
        return s

    for a, b in edges:
        aa, bb = find(a), find(b)
        parent[max(aa, bb)] = min(aa, bb)
    return {s: find(s) for s in samples}


def command(args, log):
    with log.open("a") as handle:
        handle.write("COMMAND " + repr(list(map(str, args))) + "\n")
        handle.flush()
        subprocess.run(list(map(str, args)), check=True, stdout=handle, stderr=handle)


def fingerprint(vcf, annotation, work, out, threads, minimum_depth=10):
    """Apply existing v2 site exclusions; repeat KING on deterministic site halves.

    Thresholds are conservative coverage requirements, not a calibrated claim
    that every passing library is pure RNA or every cluster is one person.
    Genotypes and temporary PLINK files stay in an owner-only temporary folder.
    """
    plink = os.environ.get("PLINK2", "plink")
    log = out / "filtered_fingerprint_commands.log"
    clean = annotation[yes(annotation.error_site) & annotation.chrom.str.fullmatch(r"chr([1-9]|1[0-9]|2[0-2])")].copy()
    clean["half"] = [int(hashlib.sha256(f"{c}:{int(p)}".encode()).hexdigest()[:2], 16) % 2
                     for c, p in zip(clean.chrom, clean.pos)]
    write(clean[["chrom", "pos", "half"]], out / "filtered_identity_site_selection.tsv")
    pairs, counts = {}, []
    for label in ["full", "half0", "half1"]:
        sites = clean if label == "full" else clean[clean.half.eq(int(label[-1]))]
        regions = work / f"{label}.regions.tsv"
        sites[["chrom", "pos"]].to_csv(regions, sep="\t", header=False, index=False)
        subset = work / f"{label}.vcf.gz"
        # Existing calls already use unique STAR MAPQ 255 and base Q >=20.
        # Preserve the A1 genotype depth rule in the primary screen. The DP20
        # used for homozygous-error QC is a separate strict-depth sensitivity.
        first = subprocess.Popen(["bcftools", "view", "-R", str(regions), "-Ou", str(vcf)], stdout=subprocess.PIPE)
        with log.open("a") as handle:
            subprocess.run(["bcftools", "+setGT", "-", "-Oz", "-o", str(subset), "--", "-t", "q", "-n", ".",
                            "-i", f"FMT/DP<{minimum_depth}"],
                           stdin=first.stdout, stdout=handle, stderr=handle, check=True)
        first.stdout.close()
        if first.wait() != 0:
            raise RuntimeError("bcftools view failed")
        prefix = work / label
        command([plink, "--vcf", subset, "--vcf-half-call", "m", "--max-alleles", "2", "--chr", "1-22",
                 "--set-all-var-ids", "@:#", "--make-pgen", "--threads", threads, "--out", prefix], log)
        command([plink, "--pfile", prefix, "--missing", "sample-only", "--out", prefix], log)
        miss = pd.read_csv(str(prefix) + ".smiss", sep=r"\s+")
        miss["called"] = miss.OBS_CT - miss.MISSING_CT
        counts.append(miss[["#IID", "called"]].rename(columns={"#IID": "sample_id", "called": label + "_called"}))
        # No sample pruning: keep low-coverage records explicitly unqualified.
        command([plink, "--pfile", prefix, "--make-king-table", "counts", "--king-table-filter", "0.0884",
                 "--threads", threads, "--out", prefix], log)
        table = pd.read_csv(str(prefix) + ".kin0", sep=r"\s+")
        table["pair_key"] = ["||".join(key(a, b)) for a, b in zip(table["#IID1"], table.IID2)]
        pairs[label] = table
        write(table, out / f"{label}_filtered_related_pairs.tsv")
    count = counts[0]
    for c in counts[1:]:
        count = count.merge(c, on="sample_id", validate="one_to_one")
    write(count, out / "filtered_callable_site_counts.tsv")
    return pairs, count


def qualify_pairs(original, filtered, counts, half_min_shared=MIN_HALF_SHARED):
    old = {key(a, b): r for a, b, r in original[["#IID1", "IID2", "relation"]].itertuples(index=False, name=None)}
    tables = {label: {key(row["#IID1"], row["IID2"]): row for _, row in table.iterrows()}
              for label, table in filtered.items()}
    candidate_keys = {pair for pair, relation in old.items() if relation == "same_individual"}
    callable_sites = counts.set_index("sample_id")
    candidate_keys |= {pair for pair, row in tables["full"].items() if row.KINSHIP >= SAME_KINSHIP and row.NSNP >= MIN_SHARED and
                      all(s in callable_sites.index and callable_sites.loc[s, "full_called"] >= MIN_SITES for s in pair)}
    rows, stable, unresolved = [], set(), set()
    for a, b in sorted(candidate_keys):
        row = {"sample_a": a, "sample_b": b, "original_relation": old.get((a, b), "not_original_same_identity"),
               "unit_evidence": "same_genetic_identity_screen_not_verified_person"}
        passes = True
        for label, minimum in [("full", MIN_SHARED), ("half0", half_min_shared), ("half1", half_min_shared)]:
            result = tables[label].get((a, b))
            row[label + "_shared_sites"] = int(result.NSNP) if result is not None else 0
            row[label + "_kinship"] = float(result.KINSHIP) if result is not None else np.nan
            row[label + "_ibs0"] = int(result.IBS0) if result is not None else np.nan
            passes &= result is not None and result.NSNP >= minimum and result.KINSHIP >= SAME_KINSHIP
        passes &= all(s in callable_sites.index and callable_sites.loc[s, "full_called"] >= MIN_SITES for s in [a, b])
        row["stable_filtered_same_identity"] = bool(passes)
        row["failure_reason"] = "" if passes else "insufficient_or_discordant_filtered_identity_evidence"
        rows.append(row)
        (stable if passes else unresolved).add((a, b))
    return pd.DataFrame(rows), stable, unresolved


def planted_control(out, threads):
    """Check the real VCF/PLINK route on known genotypes without using outcomes.

    This checks implementation and coverage abstention. It does not establish
    RNA contamination specificity or supply independent donor ground truth.
    """
    control = out / "controls"
    control.mkdir()
    rng = np.random.default_rng(20261005)
    n = 4000
    genotype = rng.binomial(2, .35, (n, 8))
    matrix = np.column_stack([genotype, genotype[:, 0], genotype[:, 1]])
    names = [f"planted_distinct_{i}" for i in range(8)] + ["planted_repeat_0", "planted_low_coverage_1"]
    annotation = pd.DataFrame({"chrom": "chr1", "pos": np.arange(10001, 10001 + n), "error_site": True})
    with tempfile.TemporaryDirectory(prefix="resource_identity_control_", dir=os.environ.get("TMPDIR", "/tmp")) as tmp:
        work = Path(tmp)
        vcf = work / "control.vcf"
        with vcf.open("w") as handle:
            handle.write("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=248956422>\n")
            handle.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
            handle.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n')
            handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(names) + "\n")
            for i in range(n):
                values = [f"{['0/0', '0/1', '1/1'][g]}:{30 if j < 9 or i < 200 else 5}" for j, g in enumerate(matrix[i])]
                handle.write(f"chr1\t{10001+i}\t.\tA\tC\t.\tPASS\t.\tGT:DP\t" + "\t".join(values) + "\n")
        command(["bcftools", "view", "-Oz", "-o", work / "control.vcf.gz", vcf], control / "control_commands.log")
        command(["bcftools", "index", "-t", work / "control.vcf.gz"], control / "control_commands.log")
        filtered, count = fingerprint(work / "control.vcf.gz", annotation, work, control, threads)
    original = pd.DataFrame({"#IID1": [names[0]], "IID2": [names[8]], "relation": ["same_individual"]})
    checks, stable, unresolved = qualify_pairs(original, filtered, count)
    expected = {key(names[0], names[8])}
    if stable != expected or unresolved or count.set_index("sample_id").loc[names[9], "full_called"] != 200:
        raise AssertionError("Planted repeat/distinct/low-coverage fingerprint checks failed")
    write(checks, control / "planted_identity_pair_checks.tsv")
    (control / "known_unit_check.json").write_text(json.dumps({"passed": True, "seed": 20261005, "known_distinct_genotypes": 8,
        "known_identical_pair_recovered": True, "low_coverage_sample_abstained": True,
        "limit": "implementation_check_not_RNA_mixture_truth_or_real_person_validation"}, indent=2) + "\n")


def source_roster(identity, histology, source130, source193, metadata):
    """Preserve full source categories; no coarse Model A S0/S1/S2 codes."""
    m = identity.rename(columns={"run": "sample_id", "cohort": "dataset", "individual_id": "a1_component"}).copy()
    meta = metadata.rename(columns={"sex": "reported_sex"})
    cols = ["sample_id", "condition", "group_binary", "reported_sex", "age", "fibrosis_stage", "nas_score", "diagnosis_harmonized"]
    m = m.merge(meta[cols], on="sample_id", how="left", validate="one_to_one")
    m["source_fibrosis_category"] = m.fibrosis_stage.map(lambda v: f"F{int(v)}" if pd.notna(v) else "Not recorded")
    grouped = m.dataset.eq("GSE213621")
    m.loc[grouped, "source_fibrosis_category"] = m.loc[grouped, "condition"].map({
        "Control": "Control", "Fibrosis_F0F1": "F0–F1", "Fibrosis_F2": "F2", "Fibrosis_F3F4": "F3–F4"})
    if m.loc[grouped, "source_fibrosis_category"].isna().any():
        raise ValueError("Unknown GSE213621 source fibrosis category")
    m["fibrosis_stage_exact"] = m.fibrosis_stage
    m.loc[grouped & ~m.source_fibrosis_category.eq("F2"), "fibrosis_stage_exact"] = np.nan
    m["stage_resolution"] = np.where(grouped, "grouped", np.where(m.fibrosis_stage_exact.notna(), "exact", "not_recorded"))
    # Candidate histology supplies corrected source controls and assay/model sex.
    h = histology.set_index("sample_id")
    for column in ["group_binary", "diagnosis_harmonized", "inferred_sex", "sex_source"]:
        old = m[column] if column in m else pd.Series(np.nan, index=m.index)
        m[column] = m.sample_id.map(h[column]).combine_first(old)
    s130 = source130.set_index("run")
    in130 = m.dataset.eq("GSE130970")
    m.loc[in130, "group_binary"] = np.where(m.loc[in130, "sample_id"].map(s130.source_control_status).eq("Control"), "Control", "Disease")
    for target, source in [("nas_score", "nafld_activity_score"), ("fibrosis_stage", "fibrosis_stage"), ("age", "age_at_biopsy")]:
        m.loc[in130, target] = m.loc[in130, "sample_id"].map(s130[source])
    m.loc[in130, "reported_sex"] = m.loc[in130, "sample_id"].map(s130.sex_geo).map({"female": "F", "male": "M"})
    m.loc[in130, "source_fibrosis_category"] = m.loc[in130, "fibrosis_stage"].map(lambda v: f"F{int(v)}" if pd.notna(v) else "Not recorded")
    m.loc[in130, "fibrosis_stage_exact"] = m.loc[in130, "fibrosis_stage"]
    m["source_control"] = m.group_binary.eq("Control")
    m["source_person_id"] = pd.Series(None, index=m.index, dtype=object)
    m["first_biopsy"] = "unknown"
    m["first_biopsy_source"] = "not_recorded"
    c193 = source193.set_index("run_id")
    in193 = m.dataset.eq("GSE193066")
    m.loc[in193, "source_person_id"] = m.loc[in193, "sample_id"].map(c193.participant_token)
    m.loc[in193, "first_biopsy"] = m.loc[in193, "sample_id"].map(c193.is_first_biopsy).map({True: "True", False: "False", "TRUE": "True", "FALSE": "False"})
    m.loc[in193, "first_biopsy_source"] = "source_participant_crosswalk"
    m["old_identity_marker_flag"] = m.expr_sex.eq("both") | m.other_allele_fraction_at_hom.gt(.02)
    m["old_marker_interpretation"] = "diagnostic_only_not_mixture_truth_or_exclusion"
    return m


def source_conflicts(frame):
    columns = ["dataset", "group_binary", "source_fibrosis_category", "nas_score", "reported_sex", "diagnosis_harmonized"]
    # Documented repeated biopsies may change phenotype. Their source person
    # crosswalk still has to pass the separate cohort-wide genotype checks.
    if frame.dataset.eq("GSE193066").all() and frame.source_person_id.notna().all() and frame.source_person_id.nunique() == 1:
        columns = ["dataset", "reported_sex"]
    return [c for c in columns if frame[c].dropna().astype(str).nunique() > 1]


def select_context(roster, samples, name):
    """Choose within the actual input roster, independently of effect estimates."""
    if len(samples) != len(set(samples)):
        raise ValueError(f"Duplicate samples in {name} input")
    c = roster[roster.sample_id.isin(samples)].copy()
    if set(c.sample_id) != set(samples):
        raise ValueError(f"Unjoined {name} samples: {sorted(set(samples)-set(c.sample_id))}")
    c["primary_selected"] = False
    c["unit_choice"] = "excluded"
    eligible = c[c.identity_qualified].copy()
    eligible["first_rank"] = eligible.first_biopsy.map({"True": 2, "unknown": 1, "False": 0}).fillna(1)
    for _, group in eligible.groupby("filtered_component", sort=True):
        ranked = group.sort_values(["first_rank", "full_called", "sample_id"], ascending=[False, False, True])
        chosen = ranked.index[0]
        choice = "single_available_library"
        if len(ranked) > 1:
            a, b = ranked.iloc[0], ranked.iloc[1]
            choice = "documented_first_biopsy" if a.first_rank != b.first_rank else ("callable_site_coverage" if a.full_called != b.full_called else "run_id")
        c.loc[chosen, "primary_selected"] = True
        c.loc[chosen, "unit_choice"] = choice
        extra = ranked.index[1:]
        c.loc[extra, "exclusion_reason"] = "concordant_repeated_genetic_identity_library"
        c.loc[extra, "unit_choice"] = "not_selected_outcome_free_representative"
    c["context"] = name
    c["analysis_unit_id"] = "ResourceGeneticCluster::" + c.filtered_component
    c["unit_evidence"] = "source_labels_plus_stable_filtered_genetic_identity_cluster_not_verified_person"
    primary = c[c.primary_selected][["sample_id", "analysis_unit_id", "unit_evidence"]]
    if primary.sample_id.duplicated().any() or primary.analysis_unit_id.duplicated().any():
        raise ValueError("Primary unit duplication")
    return c, primary


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate-root", type=Path, required=True)
    parser.add_argument("--identity-directory", default="identity_resource")
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--genotype-depth", choices=[10, 20], type=int, default=10)
    parser.add_argument("--half-min-shared", choices=[250, 500], type=int, default=250)
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Identity reanalysis requires an allocated compute node")
    if args.threads < 1 or args.threads > 4:
        raise ValueError("Resource identity checks allow one to four threads")
    if Path(args.identity_directory).name != args.identity_directory:
        raise ValueError("Identity output directory must be a single relative name")
    out = args.candidate_root / args.identity_directory
    out.mkdir(parents=True, exist_ok=False)
    inputs = [A1 / "sample_to_individual.tsv", A1 / "kinship_related_pairs.tsv", QC / "identity_site_annotation.tsv",
              VCF, METADATA, HISTOLOGY, BULK, SOURCE130, SOURCE193, JOIN267, ENDPOINT, STAGE, Path(__file__)]
    hashes = [{"path": str(p), "sha256": hashlib.sha256(p.read_bytes()).hexdigest()} for p in inputs]
    write(pd.DataFrame(hashes), out / "inputs.tsv")
    planted_control(out, args.threads)
    with tempfile.TemporaryDirectory(prefix="resource_identity_", dir=os.environ.get("TMPDIR", "/tmp")) as tmp:
        filtered, count = fingerprint(VCF, read(QC / "identity_site_annotation.tsv"), Path(tmp), out, args.threads, args.genotype_depth)
    pair_checks, stable, unresolved = qualify_pairs(read(A1 / "kinship_related_pairs.tsv"), filtered, count, args.half_min_shared)
    write(pair_checks, out / "filtered_identity_pair_checks.tsv")
    m = source_roster(read(A1 / "sample_to_individual.tsv"), read(HISTOLOGY), read(SOURCE130), read(SOURCE193),
                      pd.read_csv(METADATA, keep_default_na=False, na_values=[""]))
    m = m.merge(count, on="sample_id", how="left", validate="one_to_one")
    m["filtered_component"] = m.sample_id.map(components(m.sample_id.tolist(), stable))
    # Quarantine the union of old and new components when an old duplicate call
    # cannot be reproduced; do not create independence by dropping a noisy edge.
    all_candidate = stable | unresolved
    m["screen_component"] = m.sample_id.map(components(m.sample_id.tolist(), all_candidate))
    bad_components = set(m.loc[m.sample_id.isin({s for p in unresolved for s in p}), "screen_component"])
    for _, group in m.groupby("filtered_component"):
        if any(key(a, b) not in stable for a, b in itertools.combinations(group.sample_id, 2)):
            bad_components.update(group.screen_component)
    m["exclusion_reason"] = np.where(m.full_called.lt(MIN_SITES) | m.full_called.isna(), "insufficient_filtered_identity_sites", "")
    m.loc[m.screen_component.isin(bad_components), "exclusion_reason"] = "unresolved_filtered_identity_component"
    conflict_rows = []
    for component, group in m.groupby("screen_component", sort=True):
        conflicts = source_conflicts(group)
        if conflicts:
            m.loc[group.index, "exclusion_reason"] = "conflicting_full_source_labels:" + ",".join(conflicts)
            conflict_rows.append({"screen_component": component, "n_libraries": len(group), "conflicting_fields": ",".join(conflicts),
                                  "sample_ids": ",".join(sorted(group.sample_id))})
    write(pd.DataFrame(conflict_rows, columns=["screen_component", "n_libraries", "conflicting_fields", "sample_ids"]), out / "source_label_conflicts.tsv")
    # A source-defined repeat-biopsy cohort must corroborate its person crosswalk
    # in every available comparison before supplying person-level replication.
    g193 = m[m.dataset.eq("GSE193066")]
    source_pairs = {key(a, b) for _, group in g193.dropna(subset=["source_person_id"]).groupby("source_person_id")
                    for a, b in itertools.combinations(sorted(group.sample_id), 2)}
    observed_cross = {p for p in stable if all(s in set(g193.sample_id) for s in p) and
                      m.set_index("sample_id").loc[p[0], "source_person_id"] != m.set_index("sample_id").loc[p[1], "source_person_id"]}
    source193_qualified = bool(source_pairs and source_pairs.issubset(stable) and not observed_cross and
                               g193.source_person_id.notna().all() and g193.exclusion_reason.eq("").all())
    if not source193_qualified:
        m.loc[m.dataset.eq("GSE193066"), "exclusion_reason"] = "GSE193066_source_person_crosswalk_not_qualified"
    m["identity_qualified"] = m.exclusion_reason.eq("")
    m["genetic_person_verified"] = False
    m["identity_evidence"] = np.where(m.identity_qualified, "stable_filtered_genetic_identity_screen", "unresolved_or_source_conflicting_identity")
    write(m, out / "resource_identity_roster.tsv")
    contexts = {"bulk": read(BULK).sample_id.tolist(), "continuum": read(HISTOLOGY).sample_id.tolist(),
                "endpoint": read(ENDPOINT).sample_id.tolist(), "stage": read(STAGE).sample_id.drop_duplicates().tolist(),
                "all_screened": m.sample_id.tolist()}
    census = []
    for name, samples in contexts.items():
        roster, primary = select_context(m, samples, name)
        write(roster, out / f"{name}_roster.tsv")
        write(primary, out / f"{name}_primary_sample_ids.tsv")
        sensitivity = roster.copy()
        sensitivity["evidence_state"] = "historical_all_library_roster_conditional_sensitivity_not_independent_units"
        write(sensitivity, out / f"{name}_conditional_sensitivity_roster.tsv")
        for cohort, group in roster.groupby("dataset"):
            selected = group[group.primary_selected]
            census.append({"context": name, "dataset": cohort, "source_libraries": len(group), "primary_clusters": len(selected),
                           "primary_control_clusters": int(selected.source_control.sum()), "primary_disease_clusters": int((~selected.source_control).sum()),
                           "verified_person_count": "not_established", "excluded_libraries": int((~group.primary_selected).sum())})
    write(pd.DataFrame(census), out / "analysis_unit_census.tsv")
    join267 = read(JOIN267)
    if join267.participant_id.duplicated().any():
        raise ValueError("Authoritative GSE267145 participant join has duplicates")
    join267["unit_evidence"] = "authoritative_source_participant_pairing_not_screened_in_A1"
    join267["genetic_identity_screen_available"] = False
    write(join267, out / "GSE267145_source_participant_units.tsv")
    decisions = {
        "schema_version": "resource-identity-candidate-v1", "state": "candidate_not_adopted",
        "thresholds": {"genotype_depth": args.genotype_depth, "full_called_sites": MIN_SITES, "full_pair_shared_sites": MIN_SHARED,
                       "each_half_pair_shared_sites": args.half_min_shared, "same_genetic_identity_king": SAME_KINSHIP},
        "role": "primary_identity_screen" if (args.genotype_depth, args.half_min_shared) == (10, 250) else "strict_depth_identity_sensitivity",
        "threshold_interpretation": "conservative_coverage_and_stability_screen_not_proof_of_pure_RNA_or_person_identity",
        "old_0_02_and_sex_markers": "diagnostic_only_not_exclusion_or_confirmed_mixture",
        "phenotype_conflicts": "exclude_entire_component_no_label_adjudication_or_effect_based_representative_choice",
        "representative_order": ["source_documented_first_biopsy", "filtered_callable_site_coverage", "run_id"],
        "GSE193066_source_person_crosswalk_qualified": source193_qualified,
        "GSE193066_source_repeat_pairs": len(source_pairs), "GSE193066_source_pairs_passing_filtered_screen": len(source_pairs & stable),
        "GSE193066_filtered_pairs_crossing_source_persons": len(observed_cross),
        "GSE193066_replication_ICC_and_derived_ceilings_primary": "eligible_for_separate_refit" if source193_qualified else "omit_unqualified_cohort",
        "environment": {"python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__,
                        "slurm_job": os.environ["SLURM_JOB_ID"], "threads": args.threads},
        "primary_bulk_clusters": int(sum(c["primary_clusters"] for c in census if c["context"] == "bulk")),
        "limitations": ["KING cannot distinguish monozygotic twins from sample duplicates", "site filtering does not establish contamination truth",
                        "passing singleton screens do not establish universal person independence", "source-conflicting records are not adjudicated"]}
    (out / "analysis_definition.json").write_text(json.dumps(decisions, indent=2) + "\n")
    (out / "IDENTITY_CANDIDATE_COMPLETE").write_text("Filtered identity screen and source-qualified cluster rosters complete; not adopted.\n")
    print(json.dumps(decisions, indent=2))


if __name__ == "__main__":
    main()
