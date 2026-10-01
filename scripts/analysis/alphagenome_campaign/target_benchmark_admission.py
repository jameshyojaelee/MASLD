#!/usr/bin/env python3
"""Admit public measured element--gene labels; no fitting or nucleotide inference.

The manifest pins public author blobs. All deposited rows are retained, including
low-power, non-distal and repeated records. Eligibility never favors significance.
Run on a compute node. Do not reinterpret this cell-line assay as liver validation.
"""
import argparse
import collections
import csv
import datetime
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import sys
import urllib.request


MISSING = {"", "NA", "NaN", "nan", "None", "."}


def present(value):
    return value is not None and str(value) not in MISSING


def first(row, *keys):
    return next((row[k] for k in keys if present(row.get(k))), "")


def number(value):
    try:
        x = float(value)
        return x if math.isfinite(x) else None
    except (ValueError, TypeError):
        return None


def boolean(value):
    if str(value).upper() in {"TRUE", "T", "1"}:
        return True
    if str(value).upper() in {"FALSE", "F", "0"}:
        return False
    return None


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_rows(path):
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def get_asset(asset, destination):
    """Validate exact pinned bytes and Git blob identity; never overwrite."""
    path = destination / asset["filename"]
    expected = asset["bytes"]
    h = hashlib.sha1(f"blob {expected}\0".encode())
    sha = hashlib.sha256()
    count = 0
    request = urllib.request.Request(asset["url"], headers={"User-Agent": "MASLD-target-admission/1.0"})
    with urllib.request.urlopen(request, timeout=60) as response, path.open("xb") as out:
        for block in iter(lambda: response.read(1024 * 1024), b""):
            count += len(block)
            if count > expected:
                raise ValueError(f"asset exceeded deposited byte size: {path.name}")
            h.update(block)
            sha.update(block)
            out.write(block)
    if count != expected or h.hexdigest() != asset["git_blob_sha1"]:
        raise ValueError(f"pinned asset verification failed: {path.name}")
    return {"filename": path.name, "url": asset["url"], "commit": asset["commit"],
            "bytes": count, "git_blob_sha1": h.hexdigest(), "sha256": sha.hexdigest(),
            "retrieved_utc": datetime.datetime.now(datetime.timezone.utc).isoformat()}


def normalize(row, asset, row_number):
    dc = asset["format"] == "dc_tap"
    historical = asset["role"] == "historical_selected_comparator"
    chrom = first(row, "targeting_chr_hg38") if dc else first(row, "chrom")
    start = first(row, "targeting_start_hg38") if dc else first(row, "chromStart")
    end = first(row, "targeting_end_hg38") if dc else first(row, "chromEnd")
    symbol = first(row, "gene_symbol", "measuredGeneSymbol")
    gene_id = first(row, "gene_id", "measuredEnsemblID", "measuredGeneEnsemblId")
    cell = first(row, "cell_type", "CellType") or asset.get("cell_type", "")
    gene_chrom = first(row, "chrTSS_hg38", "chrTSS")
    tss_start = first(row, "startTSS_hg38", "startTSS")
    tss_end = first(row, "endTSS_hg38", "endTSS")
    strand = first(row, "strandGene", "gene_strand")
    distance = first(row, "distance_to_gencode_gene_TSS", "distToTSS", "distanceToTSS")
    significance = first(row, "significant", "Significant")
    native_effect = first(row, "pct_change_effect_size") if dc else first(row, "EffectSize")
    units = "percent_change_expression" if dc else asset["effect_units"]
    power = first(row, "power_at_effect_size_15", "PowerAtEffectSize15")
    native_valid = first(row, "DistalElement_Gene") if dc else first(row, "ValidConnection")
    reasons = []
    s, e, d = number(start), number(end), number(distance)
    if not chrom or s is None or e is None or s != int(s) or e != int(e) or s < 0 or e <= s:
        reasons.append("invalid_or_missing_native_element_coordinates")
    if not gene_id:
        reasons.append("missing_exact_ensembl_gene_id")
    if not symbol:
        reasons.append("missing_source_gene_symbol")
    if not cell:
        reasons.append("missing_assay_context")
    if gene_chrom != chrom:
        reasons.append("missing_or_non_cis_gene_chromosome")
    if number(tss_start) is None or number(tss_end) is None:
        reasons.append("missing_native_TSS_coordinates")
    if boolean(native_valid) is not True:
        reasons.append("source_not_distal_valid_or_validity_unknown")
    if d is None or not 1000 < abs(d) <= 1_000_000:
        reasons.append("outside_or_missing_prespecified_1kb_to_1Mb_distance")
    # Source roles are determined from processing paths, never from this row's outcome.
    if asset["role"] != "prefilter_source":
        reasons.append("source_population_already_filtered_or_unverified")
    geometric = not reasons
    p = number(power)
    powered = geometric and p is not None and .8 <= p <= 1
    sig, effect = boolean(significance), number(native_effect)
    if sig is None or effect is None:
        label = "unknown_measurement"
    elif sig and effect < 0:
        label = "measured_expression_decrease"
    elif sig:
        label = "measured_other_direction"
    elif p is None or not .8 <= p <= 1:
        label = "nonsignificant_insufficient_or_unknown_power"
    else:
        label = "power_qualified_nonsignificant"
    # Uniform power restriction applies to positives too. The public benchmark's
    # significance-dependent power rescue is deliberately not applied.
    analysis = powered and sig is not None and effect is not None
    decrease_label = int(sig and effect < 0) if sig is not None and effect is not None else ""
    any_change_label = int(sig) if sig is not None and effect is not None else ""
    key_gene = gene_id or ("source_symbol_only:" + symbol)
    tuple_key = "|".join(("GRCh38", chrom, start, end, key_gene, cell))
    native_hash = hashlib.sha256(json.dumps(row, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    return {
        "source_asset": asset["id"], "source_row": row_number, "source_row_sha256": native_hash,
        "source_commit": asset["commit"], "source_family": asset["family"],
        "source_dataset_native": first(row, "Dataset", "Reference") or asset["family"],
        "source_role": asset["role"], "assay": asset["assay"], "cell_type": cell,
        "context_status": "cultured_cell_line_not_primary_adult_liver", "genome_build": "GRCh38",
        "element_chrom_source": chrom, "element_start_source": start, "element_end_source": end,
        "element_coordinate_convention": asset["element_coordinates"],
        "element_id_source": first(row, "intended_target_name_hg38", "PerturbationTargetID", "name"),
        "gene_id_source": gene_id, "gene_symbol_source": symbol,
        "gene_chrom_source": gene_chrom, "gene_TSS_start_source": tss_start,
        "gene_TSS_end_source": tss_end, "gene_strand_source": strand or "unknown",
        "TSS_coordinate_convention": asset["tss_coordinates"],
        "distance_to_TSS_source": distance, "distance_units": "bp_source_convention",
        "effect_native": native_effect, "effect_units": units,
        "effect_direction": "perturbed_element_relative_to_source_control",
        "effect_SE_native": first(row, "standard_error_pct_change") if dc else "",
        "effect_CI95_low_native": first(row, "lower_CI_95_pct_change") if dc else first(row, "EffectSize95ConfidenceIntervalLow"),
        "effect_CI95_high_native": first(row, "upper_CI95_pct_change", "upper_CI_95_pct_change") if dc else first(row, "EffectSize95ConfidenceIntervalHigh"),
        "p_adjusted_source": first(row, "sceptre_adj_p_value", "pValueAdjusted"),
        "significance_source": significance, "source_significance_definition": asset["significance"],
        "power_15percent_source": power, "power_definition": asset["power_definition"],
        "distal_validity_source": native_valid, "measurement_state": label,
        "source_random_distal_design": first(row, "Random_DistalElement_Gene"),
        "source_positive_control_distal_design": first(row, "Positive_Control_DistalElement_Gene"),
        "source_design_type": first(row, "design_file_type", "target_type"),
        "candidate_geometry_eligible": geometric, "candidate_exclusion_reasons": ";".join(reasons),
        "uniform_power_eligible": powered, "binary_endpoint_eligible": analysis,
        "expression_decrease_binary_label": decrease_label,
        "any_expression_change_binary_label": any_change_label,
        "decrease_label_zero_interpretation": "not_a_supported_expression_loss_link;significant_increase_is_separate_measured_state_not_no_effect",
        "DHS_RPM_source": first(row, "DHS.RPM"), "DHS_percentile_source": first(row, "DHS.percentile"),
        "activity_scope": "historical_outcome_selected_annotation_only" if historical else "direct_source_field_if_present",
        "source_training_benchmark": "ENCODE_rE2G_scE2G" if asset["id"] == "historical_training" else "",
        "model_exposure_status": asset["exposure"], "nucleotide_REF": "", "nucleotide_ALT": "",
        "nucleotide_effect_supported": False, "element_gene_context_key_source": tuple_key,
        "duplicate_handling": "all_source_rows_retained_no_significance_priority_no_independent_replication_claim",
    }


def invariants():
    asset = dict(id="fixture", format="encode", role="prefilter_source", family="fixture",
                 cell_type="K562", assay="CRISPRi", effect_units="log2_fold_change_expression",
                 commit="fixture", element_coordinates="BED0", tss_coordinates="source",
                 significance="fixture", power_definition="fixture", exposure="fixture")
    r = dict(chrom="chr1", chromStart="100", chromEnd="200", chrTSS="chr1", startTSS="5000",
             endTSS="5001", measuredGeneSymbol="TEST", measuredEnsemblID="ENSG_fixture",
             EffectSize="-0.2", Significant="TRUE", PowerAtEffectSize15="0.79",
             ValidConnection="TRUE", distToTSS="4850")
    a = normalize(r, asset, 1)
    assert a["candidate_geometry_eligible"] and not a["uniform_power_eligible"]
    r["Significant"] = "FALSE"
    b = normalize(r, asset, 1)
    assert a["candidate_geometry_eligible"] == b["candidate_geometry_eligible"]
    assert a["uniform_power_eligible"] == b["uniform_power_eligible"]
    assert not b["binary_endpoint_eligible"]
    r["PowerAtEffectSize15"] = "NA"
    assert normalize(r, asset, 1)["measurement_state"] == "nonsignificant_insufficient_or_unknown_power"
    r["PowerAtEffectSize15"] = ".8"
    assert normalize(r, asset, 1)["binary_endpoint_eligible"]
    r["Significant"] = "TRUE"
    negative_effect = normalize(r, asset, 1)
    r["EffectSize"] = ".2"
    positive_effect = normalize(r, asset, 1)
    assert negative_effect["binary_endpoint_eligible"] == positive_effect["binary_endpoint_eligible"] is True
    assert negative_effect["candidate_geometry_eligible"] == positive_effect["candidate_geometry_eligible"]
    assert negative_effect["expression_decrease_binary_label"] == 1
    assert positive_effect["expression_decrease_binary_label"] == 0
    assert negative_effect["any_expression_change_binary_label"] == positive_effect["any_expression_change_binary_label"] == 1
    r["measuredEnsemblID"] = ""
    assert not normalize(r, asset, 1)["candidate_geometry_eligible"]
    assert normalize(r, asset, 1)["gene_strand_source"] == "unknown"
    assert not normalize(r, asset, 1)["nucleotide_effect_supported"]
    return {"outcome_independent_geometry_and_uniform_power": "passed",
            "sign_reversal_preserves_endpoint_population": "passed",
            "missing_power_not_negative": "passed", "missing_gene_id_not_inferred": "passed",
            "unknown_strand_and_no_nucleotide_claim": "passed"}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Scientific admission checks must run on a compute allocation")
    args.out.mkdir(parents=True, exist_ok=False)
    (args.out / "source").mkdir()
    manifest = json.loads(args.manifest.read_text())
    if sum(a["bytes"] for a in manifest["assets"]) > 25_000_000:
        raise ValueError("fixed acquisition ceiling exceeded")
    receipts = []
    for asset in manifest["assets"]:
        try:
            receipts.append(get_asset(asset, args.out / "source"))
        except Exception as error:
            partial = args.out / "source" / asset["filename"]
            receipts.append({"filename": asset["filename"], "url": asset["url"],
                             "commit": asset["commit"], "status": "failed_retained_if_created",
                             "error": repr(error), "bytes": partial.stat().st_size if partial.exists() else 0,
                             "sha256": digest(partial) if partial.exists() else None})
            (args.out / "download_receipts.json").write_text(json.dumps(receipts, indent=2) + "\n")
            raise
        (args.out / "download_receipts.json").write_text(json.dumps(receipts, indent=2) + "\n")
    counts = collections.Counter()
    state_hashes = collections.defaultdict(set)
    historical = collections.defaultdict(set)
    summaries = []
    # Pass one records multiplicity and historical exposure using exact source tuples.
    # This neither chooses a preferred measurement nor merges nearby elements.
    for asset in manifest["assets"]:
        for i, raw in enumerate(read_rows(args.out / "source" / asset["filename"]), 1):
            r = normalize(raw, asset, i)
            key = (asset["id"], r["element_gene_context_key_source"])
            counts[key] += 1
            state_hashes[key].add(r["source_row_sha256"])
            if asset["role"] == "historical_selected_comparator":
                historical[r["element_gene_context_key_source"]].add(asset["id"])
    for asset in manifest["assets"]:
        counter = collections.Counter()
        genes, elements, cells, chromosomes, headers = set(), set(), set(), set(), []
        outfile = args.out / (asset["id"] + ".normalized.tsv.gz")
        with gzip.open(outfile, "wt", newline="") as handle:
            writer = None
            for i, raw in enumerate(read_rows(args.out / "source" / asset["filename"]), 1):
                if not headers:
                    headers = list(raw)
                r = normalize(raw, asset, i)
                key = (asset["id"], r["element_gene_context_key_source"])
                r["same_tuple_source_row_multiplicity"] = counts[key]
                r["same_tuple_distinct_full_row_states"] = len(state_hashes[key])
                r["exact_source_tuple_historical_exposure"] = ";".join(sorted(historical.get(key[1], [])))
                r["exposure_join_limit"] = "exact_native_tuple_only_nearby_or_resized_element_and_gene_alias_overlap_unresolved"
                if writer is None:
                    writer = csv.DictWriter(handle, fieldnames=list(r), delimiter="\t")
                    writer.writeheader()
                writer.writerow(r)
                counter["rows"] += 1
                for flag in ["candidate_geometry_eligible", "uniform_power_eligible", "binary_endpoint_eligible"]:
                    counter[flag] += int(r[flag])
                counter["state:" + r["measurement_state"]] += 1
                counter["missing_power15"] += int(not present(r["power_15percent_source"]))
                counter["missing_gene_id"] += int(not r["gene_id_source"])
                counter["missing_activity_DHS_RPM"] += int(not present(r["DHS_RPM_source"]))
                counter["missing_native_SE"] += int(not present(r["effect_SE_native"]))
                counter["rows_with_repeated_native_tuple"] += int(counts[key] > 1)
                genes.add(r["gene_id_source"] or "symbol_only:" + r["gene_symbol_source"])
                elements.add((r["element_chrom_source"], r["element_start_source"], r["element_end_source"]))
                cells.add(r["cell_type"])
                chromosomes.add(r["element_chrom_source"])
        summaries.append(dict(asset=asset["id"], counts=dict(counter), genes=len(genes), elements=len(elements),
                              cell_types=sorted(cells), chromosomes=sorted(chromosomes), native_columns=headers,
                              normalized_sha256=digest(outfile)))
    receipt = {"status": "admission_inventory_complete_no_model_fit", "job_id": os.environ["SLURM_JOB_ID"],
               "scientific_invariants": invariants(), "assets": summaries,
               "manifest_sha256": digest(args.manifest), "script_sha256": digest(Path(__file__)),
               "python": sys.version, "platform": platform.platform(),
               "candidate_rule": "source-prefilter distal valid cis pairs with exact gene ID and native 1kb<distance<=1Mb; uniform >=80% power at 15% applied regardless of significance",
               "replication": "biological replication remains source experiment; repeated rows, cells and candidate pairs are not independent donors",
               "benchmark_status": "not_run; exact feature coverage, source/version exposure and interval-coordinate review precede fitting",
               "limits": ["No primary adult-liver target validation", "No nucleotide alleles or nucleotide-effect labels",
                          "No conversion or reconstruction of missing effect scales, SEs, phase or power",
                          "Historical combined data are outcome-selected, never the new primary candidate population",
                          "Exact tuple exposure is a lower bound; source-family and overlapping/resized intervals must remain grouped",
                          "A count of power-qualified pairs does not establish achievable precision or independent biological sample size"]}
    (args.out / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({"status": receipt["status"], "source_bytes": sum(x["bytes"] for x in receipts),
                      "total_rows": sum(x["counts"].get("rows", 0) for x in summaries)}, indent=2))


if __name__ == "__main__":
    main()
