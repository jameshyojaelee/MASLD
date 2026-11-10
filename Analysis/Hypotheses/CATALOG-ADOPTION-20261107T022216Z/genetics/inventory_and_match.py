#!/usr/bin/env python3
"""Bounded saved-signal inventory; no fitting and no genetic class changes."""
import csv
import gzip
import hashlib
import itertools
import json
import math
import os
import platform
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
FM = ROOT / "GWAS/finemapping"
R2 = FM / "results/susie_coloc_r2_20260924T053244Z"
REPLAY = ROOT / "Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-candidate-2026-08-11-r1/genetics"
RV9 = ROOT / "Analysis/Hypotheses/RV9-20261102T151012Z/out-20261102T151110Z"
OUT = Path(sys.argv[1])


def read(path, delimiter="\t"):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def write(name, rows):
    assert rows, name
    with (OUT / name).open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


sources = [HERE / "SPEC.md", Path(__file__), R2 / "susie_coloc_signal_pairs.tsv.gz",
           R2 / "gene_level_coloc_tier12.csv", REPLAY / "prepared/replay_plan.tsv",
           RV9 / "rows.tsv", RV9 / "gene_states.tsv", FM / "config/gwas_registry.tsv",
           FM / "config/gwas_trait_tier.tsv", FM / "config/gwas_af_sources.tsv",
           FM / "src/06_susie_coloc.R", FM / "src/coloc_signal_eligibility.R"]
initial = {str(path): sha(path) for path in sources}
OUT.mkdir(parents=True, exist_ok=False)
if sha(HERE / "SPEC.md") != (HERE / "SPEC.sha256").read_text().split()[0]:
    raise RuntimeError("specification changed after hashing")

genes = read(RV9 / "gene_states.tsv")
assert len(genes) == 436 and len({row["ensembl"] for row in genes}) == 436
native = read(R2 / "gene_level_coloc_tier12.csv", ",")
assert {row["ensembl"] for row in native if row["coloc_best_susie_pp4"] and float(row["coloc_best_susie_pp4"]) > .5} == {row["ensembl"] for row in genes}
signals = defaultdict(list)
for row in read(R2 / "susie_coloc_signal_pairs.tsv.gz"):
    if row["eligible"] == "TRUE" and row["PP.H4"] not in {"", "NA"} and float(row["PP.H4"]) > .5:
        signals[row["gwas"], row["ensembl"]].append(row)
plan = {(row["gwas_name"], row["ensembl"]): row for row in read(REPLAY / "prepared/replay_plan.tsv")}
registry = {row["study_name"]: row for row in read(FM / "config/gwas_registry.tsv")}
tiers = {row["study_name"]: row for row in read(FM / "config/gwas_trait_tier.tsv")}
af = {row["study_name"]: row for row in read(FM / "config/gwas_af_sources.tsv")}

inventory = []
for study, row in registry.items():
    if int(tiers[study]["tier"]) > 2:
        continue
    original = FM / row["sumstats_path"]
    augmented = Path(af.get(study, {}).get("af_sumstats_path", "/not_available"))
    selected = augmented if augmented.is_file() else original
    inventory.append(dict(study=study, trait=tiers[study]["trait"], ancestry=row["ancestry"],
        genome_build="GRCh37/hg19 in production harmonization", registry_ld_panel=row["ld_panel"],
        production_ld_dispatch="polyfun_eur for EUR; recorded Gram panels for non-EUR; registry label alone is not execution provenance",
        original_sumstats_path=str(original), selected_sumstats_path=str(selected),
        selected_exists=selected.is_file(), selected_bytes=selected.stat().st_size if selected.is_file() else "",
        selected_mtime_ns=selected.stat().st_mtime_ns if selected.is_file() else "",
        n_total=row["N_tot"], n_cases=row["N_cases"], source_matched_in_sample_ld="not_established_available"))
write("study_input_inventory.tsv", inventory)


def export_vectors(study, ensembl, selected):
    """Read only the relevant historical export, retaining exact allele pairs."""
    p = plan.get((study, ensembl))
    if p is None:
        return {}, "historical_export_not_planned", []
    folder = REPLAY / "replay_execution" / f"batch_{p['batch_id']}" / "exports" / study / f"chr{p['chr']}" / ensembl
    pairs_file, variants_file = folder / "signal_pairs.tsv", folder / "variant_posteriors.tsv.gz"
    if not pairs_file.is_file() or not variants_file.is_file():
        return {}, "historical_export_missing", []
    pairs = read(pairs_file)
    index = {}
    for row in selected:
        matches = [pair for pair in pairs if pair["idx1"] == row["idx1"] and pair["idx2"] == row["idx2"]
                   and pair["gwas_name"] == study and pair["ensembl"] == ensembl
                   and abs(float(pair["PP.H4.abf"]) - float(row["PP.H4"])) <= 1e-4]
        if len(matches) == 1:
            index[matches[0]["signal_pair_index"]] = (row["idx1"], row["idx2"])
    vectors = {key: {} for key in index.values()}
    invalid = set()
    with gzip.open(variants_file, "rt", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = index.get(row["signal_pair_index"])
            if key is None:
                continue
            if row["gwas_name"] != study or row["ensembl"] != ensembl:
                invalid.add(key); continue
            try:
                probability = float(row["SNP.PP.H4"])
            except ValueError:
                invalid.add(key); continue
            if not math.isfinite(probability) or probability < 0:
                invalid.add(key); continue
            alleles = tuple(sorted([row["allele1"].upper(), row["allele2"].upper()]))
            variant = (row["chr"], row["hg19_position"], *alleles)
            if variant in vectors[key]:
                invalid.add(key)
            vectors[key][variant] = probability
    for key, vector in list(vectors.items()):
        if key in invalid or not vector or abs(sum(vector.values()) - 1) > 1e-5:
            del vectors[key]
    provenance = [dict(study=study, ensembl=ensembl, path=str(path), bytes=path.stat().st_size, sha256=sha(path))
                  for path in [pairs_file, variants_file]]
    return vectors, "historical_export_validated" if vectors else "historical_export_no_compatible_valid_vectors", provenance


def credible_set(vector):
    result, cumulative = set(), 0.0
    for key in sorted(vector, key=lambda key: (-vector[key], key)):
        result.add(key)
        cumulative += vector[key]
        if cumulative >= .95:
            break
    return result


comparisons = defaultdict(list)
for row in read(RV9 / "rows.tsv"):
    comparisons[tuple([row["ensembl"]] + sorted([row["discovery"], row["replication"]]))].append(row)
pair_output, signal_output, export_provenance = [], [], []
eqtl_inventory = {}
for number, ((ensembl, first, second), rows) in enumerate(sorted(comparisons.items()), 1):
    trait = tiers[first]["trait"]
    assert trait == tiers[second]["trait"]
    s1, s2 = signals[first, ensembl], signals[second, ensembl]
    recurrent = bool(s1 and s2)
    same_eqtl = any(a["idx2"] == b["idx2"] for a, b in itertools.product(s1, s2))
    vectors1, state1, prov1 = export_vectors(first, ensembl, s1) if recurrent else ({}, "not_needed_nonrecurrent_gene", [])
    vectors2, state2, prov2 = export_vectors(second, ensembl, s2) if recurrent else ({}, "not_needed_nonrecurrent_gene", [])
    export_provenance.extend(prov1 + prov2)
    for a, b in itertools.product(s1, s2):
        v1, v2 = vectors1.get((a["idx1"], a["idx2"])), vectors2.get((b["idx1"], b["idx2"]))
        overlap, intersects = "", ""
        if v1 is not None and v2 is not None:
            overlap = sum(min(probability, v2.get(key, 0.0)) for key, probability in v1.items())
            intersects = bool(credible_set(v1) & credible_set(v2))
        signal_output.append(dict(ensembl=ensembl, gene=rows[0]["gene"], trait=trait,
            trait_group="direct_disease" if trait == "NAFLD" else "liver_enzyme", study1=first, study2=second,
            gwas_index1=a["idx1"], gwas_index2=b["idx1"], eqtl_index1=a["idx2"], eqtl_index2=b["idx2"],
            same_recorded_eqtl_signal=a["idx2"] == b["idx2"], pp4_1=a["PP.H4"], pp4_2=b["PP.H4"],
            gwas_anchor1=a["hit1"], gwas_anchor2=b["hit1"], eqtl_anchor1=a["hit2"], eqtl_anchor2=b["hit2"],
            shared_variant_posterior_overlap=overlap, shared_variant_95pct_sets_intersect=intersects,
            gwAS_only_credible_sets_available=False,
            terminal_mechanism_state="unresolved_gwas_only_credible_sets_and_adequate_LD_not_established",
            historical_export_state1=state1, historical_export_state2=state2))
    selected = s1 or s2
    chrom = selected[0]["chr"] if selected else ""
    eqtl = FM / "results/eqtl_susie_polyfun" / f"chr{chrom}" / f"{ensembl}_susie.rds"
    if ensembl not in eqtl_inventory:
        eqtl_inventory[ensembl] = dict(ensembl=ensembl, gene=rows[0]["gene"], path=str(eqtl),
            exists=eqtl.is_file(), bytes=eqtl.stat().st_size if eqtl.is_file() else "",
            current_file_sha256=sha(eqtl) if eqtl.is_file() else "",
            historical_file_hash_verified=False,
            limitation="Current file hash recorded; historical fit identity rests on recorded common production source and signal index, not a newly verified archived file digest")
    pair_output.append(dict(ensembl=ensembl, gene=rows[0]["gene"], trait=trait,
        trait_group="direct_disease" if trait == "NAFLD" else "liver_enzyme", study1=first, study2=second,
        original_directional_rows=len(rows), gene_level_recurrence=recurrent,
        supported_signals1=len(s1), supported_signals2=len(s2), same_recorded_eqtl_signal=same_eqtl,
        export_state1=state1, export_state2=state2,
        terminal_state="same_eqtl_index_gwas_signal_unresolved" if same_eqtl else
            "different_supported_eqtl_indices_mechanism_unresolved" if recurrent else
            "partner_signal_unavailable_or_not_supported_not_a_signal_identity_test"))
    if number % 25 == 0:
        print(f"matched saved comparison {number}/{len(comparisons)}", flush=True)

write("unordered_gene_study_comparisons.tsv", pair_output)
write("supported_signal_pair_comparisons.tsv", signal_output)
write("eqtl_fit_inventory.tsv", list(eqtl_inventory.values()))
write("historical_export_sources.tsv", list({row["path"]: row for row in export_provenance}.values()))
write("all_436_genes_reporting_universe.tsv", genes)
assert len({(row["ensembl"], row["study1"], row["study2"]) for row in pair_output}) == len(pair_output)
assert all(row["study1"] < row["study2"] for row in pair_output)
assert {str(path): sha(path) for path in sources} == initial
summary = dict(status="completed_saved_signal_qualification_only", no_new_fit=True,
    genes_preserved=436, unordered_gene_study_comparisons=len(pair_output), supported_signal_comparisons=len(signal_output),
    trait_counts=dict(Counter(row["trait"] for row in pair_output)),
    gene_recurrence_comparisons=sum(row["gene_level_recurrence"] for row in pair_output),
    same_recorded_eqtl_signal_comparisons=sum(row["same_recorded_eqtl_signal"] for row in pair_output),
    terminal_states=dict(Counter(row["terminal_state"] for row in pair_output)),
    no_final_signal_mechanism_replication_claim=True, source_matched_MVP_LD="not_established_available",
    input_sha256=initial, python=platform.python_version(), platform=platform.platform(), slurm_job_id=os.getenv("SLURM_JOB_ID"),
    limitation="Posterior comparisons use historical unsigned shared-variant H4 vectors validated against retained r2 signal pairs; they are not GWAS-only credible sets or independent mechanistic validation.")
(OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps({key: value for key, value in summary.items() if key != "input_sha256"}, indent=2))
