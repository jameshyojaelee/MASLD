#!/usr/bin/env python3
"""Build a label-separated, no-fit GSE274114 etiology fixture and audit."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import gzip
from hashlib import sha256
import io
import json
import math
from pathlib import Path
import re


class GSE274114ActivationError(RuntimeError):
    """Raised when the source-native activation requirement differs."""


EXPECTED_GROUPS = {"CTRL": 9, "ENEG": 11, "NASH": 10, "ENEG_NASH": 9}
GROUP_SEMANTICS = {
    "CTRL": ("healthy_liver_donor_control", "not_observed", "not_observed"),
    "ENEG": ("hbeag_negative_chronic_hbv_only", "not_observed", "observed"),
    "NASH": ("biopsy_proven_mash_only", "observed", "not_observed"),
    "ENEG_NASH": (
        "biopsy_proven_mash_plus_hbeag_negative_chronic_hbv",
        "observed",
        "observed",
    ),
}
GROUP_INSTRUMENT = {
    "CTRL": "Illumina HiSeq 4000",
    "ENEG": "Illumina HiSeq 4000",
    "NASH": "Illumina NovaSeq 6000",
    "ENEG_NASH": "Illumina NovaSeq 6000",
}
WITHIN_INSTRUMENT_CONTRAST = {
    "CTRL": "hiseq_CTRL_vs_ENEG",
    "ENEG": "hiseq_CTRL_vs_ENEG",
    "NASH": "novaseq_NASH_vs_ENEG_NASH",
    "ENEG_NASH": "novaseq_NASH_vs_ENEG_NASH",
}
ATTRIBUTE = re.compile(r'(?P<key>[A-Za-z_]+) "(?P<value>[^"]*)";')
ENSEMBL_GENE = re.compile(r"^(ENSG[0-9]+)(?:\.[0-9]+)?(?P<par>_PAR_Y)?$")


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_gene_id(value: str) -> tuple[str, bool]:
    match = ENSEMBL_GENE.fullmatch(value)
    if match is None:
        raise GSE274114ActivationError(f"non-Ensembl source feature: {value}")
    return match.group(1), match.group("par") is not None


def source_characteristics(row: dict[str, str]) -> dict[str, str]:
    result: dict[str, str] = {}
    for value in json.loads(row["characteristics_json"]):
        if ": " not in value:
            raise GSE274114ActivationError("source characteristic is malformed")
        key, item = value.split(": ", 1)
        if key in result:
            raise GSE274114ActivationError("source characteristic is duplicated")
        result[key] = item
    return result


def relation_accessions(row: dict[str, str]) -> tuple[str, str]:
    values = json.loads(row["relations_json"])
    biosamples = [
        match.group(1)
        for value in values
        if (match := re.search(r"/biosample/(SAMN[0-9]+)/?$", value))
    ]
    experiments = [
        match.group(1)
        for value in values
        if (match := re.search(r"[?&]term=(SRX[0-9]+)$", value))
    ]
    if len(biosamples) != 1 or len(experiments) != 1:
        raise GSE274114ActivationError("BioSample/SRA experiment relation differs")
    return biosamples[0], experiments[0]


def row_id(participant_id: str) -> str:
    return "g274_" + sha256(
        b"masld-bench-gse274114-row-v1\0" + participant_id.encode("utf-8")
    ).hexdigest()[:20]


def assign_folds(labels: list[dict[str, str]]) -> dict[str, int]:
    totals = [0] * 5
    group_totals: dict[str, list[int]] = defaultdict(lambda: [0] * 5)
    result: dict[str, int] = {}
    for group in sorted(EXPECTED_GROUPS):
        identifiers = sorted(row["row_id"] for row in labels if row["source_group"] == group)
        for identifier in identifiers:
            fold = min(
                range(5),
                key=lambda value: (group_totals[group][value], totals[value], value),
            )
            result[identifier] = fold
            group_totals[group][fold] += 1
            totals[fold] += 1
    return result


def parse_gtf(path: Path) -> dict[str, dict[str, str]]:
    genes: dict[str, dict[str, str]] = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="strict") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            attrs = {match["key"]: match["value"] for match in ATTRIBUTE.finditer(fields[8])}
            gene_id = attrs.get("gene_id", "")
            stable, is_par = stable_gene_id(gene_id)
            if is_par or stable in genes:
                raise GSE274114ActivationError("GENCODE v49 stable gene axis differs")
            genes[stable] = {
                "gencode_v49_gene_id": gene_id,
                "gencode_v49_gene_name": attrs.get("gene_name", ""),
                "gencode_v49_gene_type": attrs.get("gene_type", ""),
            }
    if not genes:
        raise GSE274114ActivationError("GENCODE v49 gene axis is empty")
    return genes


def build_crosswalk(features: list[str], genes: dict[str, dict[str, str]]) -> list[dict[str, str]]:
    provisional: list[dict[str, str]] = []
    for feature in features:
        stable, is_par = stable_gene_id(feature)
        target = genes.get(stable)
        state = (
            "unmapped_stable_id"
            if target is None
            else "par_y_alias_to_v49_canonical"
            if is_par
            else "stable_id_exact"
        )
        provisional.append(
            {
                "source_feature_id": feature,
                "source_stable_locus": stable + ("_PAR_Y" if is_par else ""),
                "mapping_state": state,
                "gencode_v49_stable_id": stable if target else "",
                "gencode_v49_gene_id": target["gencode_v49_gene_id"] if target else "",
                "gencode_v49_gene_name": target["gencode_v49_gene_name"] if target else "",
                "gencode_v49_gene_type": target["gencode_v49_gene_type"] if target else "",
            }
        )
    multiplicity = Counter(
        row["gencode_v49_stable_id"]
        for row in provisional
        if row["gencode_v49_stable_id"]
    )
    for row in provisional:
        size = multiplicity.get(row["gencode_v49_stable_id"], 0)
        row["target_source_row_count"] = str(size)
        row["aggregation_rule"] = (
            "not_applicable_unmapped"
            if size == 0
            else "pass_through"
            if size == 1
            else "sum_NumReads_and_TPM_before_normalization;TPM_weighted_EffectiveLength"
        )
    return provisional


def aggregate_quant_records(
    records: list[dict[str, str]], crosswalk: list[dict[str, str]]
) -> tuple[list[dict[str, float | str]], dict[str, float | int]]:
    """Map one Salmon gene file to v49 without normalization or count rounding."""
    if len(records) != len(crosswalk):
        raise GSE274114ActivationError("Salmon/crosswalk row count differs")
    by_source = {row["source_feature_id"]: row for row in crosswalk}
    if len(by_source) != len(crosswalk):
        raise GSE274114ActivationError("source feature axis is not unique")

    target_order: list[str] = []
    accumulators: dict[str, dict[str, float | int | str]] = {}
    mapped_numreads = 0.0
    mapped_tpm = 0.0
    unmapped_numreads = 0.0
    unmapped_tpm = 0.0
    for record in records:
        source_id = record["Name"]
        mapping = by_source.get(source_id)
        if mapping is None:
            raise GSE274114ActivationError("Salmon feature is outside crosswalk")
        length = float(record["Length"])
        effective_length = float(record["EffectiveLength"])
        tpm = float(record["TPM"])
        numreads = float(record["NumReads"])
        if not all(
            math.isfinite(value) and value >= 0
            for value in (length, effective_length, tpm, numreads)
        ):
            raise GSE274114ActivationError("Salmon quant value is invalid")
        target_id = mapping["gencode_v49_gene_id"]
        if not target_id:
            unmapped_numreads += numreads
            unmapped_tpm += tpm
            continue
        stable_id = mapping["gencode_v49_stable_id"]
        if stable_id not in accumulators:
            target_order.append(stable_id)
            accumulators[stable_id] = {
                "Name": target_id,
                "Length_tpm_numerator": 0.0,
                "EffectiveLength_tpm_numerator": 0.0,
                "Length_unweighted_sum": 0.0,
                "EffectiveLength_unweighted_sum": 0.0,
                "TPM": 0.0,
                "NumReads": 0.0,
                "source_rows": 0,
            }
        accumulator = accumulators[stable_id]
        if accumulator["Name"] != target_id:
            raise GSE274114ActivationError("v49 target version differs within stable gene")
        accumulator["Length_tpm_numerator"] += tpm * length
        accumulator["EffectiveLength_tpm_numerator"] += tpm * effective_length
        accumulator["Length_unweighted_sum"] += length
        accumulator["EffectiveLength_unweighted_sum"] += effective_length
        accumulator["TPM"] += tpm
        accumulator["NumReads"] += numreads
        accumulator["source_rows"] += 1
        mapped_numreads += numreads
        mapped_tpm += tpm

    output: list[dict[str, float | str]] = []
    zero_tpm_multirow_targets = 0
    for stable_id in target_order:
        accumulator = accumulators[stable_id]
        tpm = float(accumulator["TPM"])
        source_rows = int(accumulator["source_rows"])
        if tpm > 0:
            length = float(accumulator["Length_tpm_numerator"]) / tpm
            effective_length = (
                float(accumulator["EffectiveLength_tpm_numerator"]) / tpm
            )
        else:
            if float(accumulator["NumReads"]) != 0:
                raise GSE274114ActivationError(
                    "zero aggregate TPM has nonzero aggregate estimated count"
                )
            length = float(accumulator["Length_unweighted_sum"]) / source_rows
            effective_length = (
                float(accumulator["EffectiveLength_unweighted_sum"]) / source_rows
            )
            if source_rows > 1:
                zero_tpm_multirow_targets += 1
        output.append(
            {
                "Name": str(accumulator["Name"]),
                "Length": length,
                "EffectiveLength": effective_length,
                "TPM": tpm,
                "NumReads": float(accumulator["NumReads"]),
            }
        )
    expected_targets = len(
        {
            row["gencode_v49_stable_id"]
            for row in crosswalk
            if row["gencode_v49_stable_id"]
        }
    )
    if len(output) != expected_targets:
        raise GSE274114ActivationError("aggregated v49 feature count differs")
    if not math.isclose(
        sum(float(row["NumReads"]) for row in output),
        mapped_numreads,
        rel_tol=1e-12,
        abs_tol=1e-8,
    ) or not math.isclose(
        sum(float(row["TPM"]) for row in output),
        mapped_tpm,
        rel_tol=1e-12,
        abs_tol=1e-8,
    ):
        raise GSE274114ActivationError("mapped Salmon abundance was not conserved")
    return output, {
        "source_rows": len(records),
        "mapped_source_rows": len(records) - sum(
            1 for row in crosswalk if not row["gencode_v49_gene_id"]
        ),
        "output_v49_genes": len(output),
        "mapped_numreads": mapped_numreads,
        "output_numreads": sum(float(row["NumReads"]) for row in output),
        "unmapped_numreads_masked": unmapped_numreads,
        "mapped_tpm": mapped_tpm,
        "output_tpm": sum(float(row["TPM"]) for row in output),
        "unmapped_tpm_masked": unmapped_tpm,
        "zero_tpm_multirow_targets": zero_tpm_multirow_targets,
    }


def write_salmon_quant(path: Path, rows: list[dict[str, float | str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("xb") as raw_handle:
        with gzip.GzipFile(fileobj=raw_handle, mode="wb", filename="", mtime=0) as zipped:
            with io.TextIOWrapper(zipped, encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=("Name", "Length", "EffectiveLength", "TPM", "NumReads"),
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writeheader()
                for row in rows:
                    writer.writerow(
                        {
                            key: row[key]
                            if key == "Name"
                            else format(float(row[key]), ".17g")
                            for key in writer.fieldnames
                        }
                    )


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise GSE274114ActivationError(f"cannot write empty TSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=tuple(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def audit_exact_identifier_overlap(
    topology: list[dict[str, object]],
    run_join: list[dict[str, object]],
    comparison_paths: list[Path],
) -> list[dict[str, object]]:
    target = {
        "study_accession": {"GSE274114", "PRJNA1144975", "SRP524555"},
        "sample_accession": {str(row["sample_accession"]) for row in topology},
        "biosample_accession": {str(row["biosample_accession"]) for row in topology},
        "sra_experiment": {str(row["sra_experiment"]) for row in topology},
        "sra_run": {str(row["sra_run"]) for row in run_join},
        "participant_local_alias": {str(row["participant_id"]) for row in topology},
    }
    comparison_text = ""
    participant_aliases: set[str] = set()
    for configured in comparison_paths:
        candidates = [configured] if configured.is_file() else sorted(configured.rglob("*"))
        for path in candidates:
            if not path.is_file() or path.is_symlink() or path.stat().st_size > 100_000_000:
                continue
            if path.suffix.lower() not in {".toml", ".json", ".tsv", ".csv", ".txt"}:
                continue
            try:
                text = path.read_text(encoding="utf-8", errors="strict")
            except UnicodeDecodeError:
                continue
            comparison_text += "\n" + text
            if path.suffix.lower() in {".tsv", ".csv"}:
                delimiter = "\t" if path.suffix.lower() == ".tsv" else ","
                with path.open(encoding="utf-8", newline="") as handle:
                    reader = csv.DictReader(handle, delimiter=delimiter)
                    for row in reader:
                        for key in ("participant_id", "donor_id", "person_id", "source_participant_id"):
                            if row.get(key):
                                participant_aliases.add(row[key])
    results: list[dict[str, object]] = []
    for namespace, values in target.items():
        if namespace == "participant_local_alias":
            matches = sorted(values & participant_aliases)
            disposition = (
                "exact_local_alias_match_detected"
                if matches
                else "no_exact_local_alias_match_identity_still_join_unresolved"
            )
        else:
            matches = sorted(value for value in values if value in comparison_text)
            disposition = "exact_overlap_detected" if matches else "no_exact_overlap"
        results.append(
            {
                "identifier_namespace": namespace,
                "target_identifiers": len(values),
                "exact_overlap_count": len(matches),
                "matched_identifiers": ";".join(matches) if matches else "none",
                "disposition": disposition,
            }
        )
    if any(row["exact_overlap_count"] for row in results):
        raise GSE274114ActivationError("exact identifier overlap with another cohort family detected")
    return results


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--source-artifacts-sha256", required=True)
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--matrix-artifacts-sha256", required=True)
    parser.add_argument("--prior-activation", type=Path, required=True)
    parser.add_argument("--prior-activation-artifacts-sha256", required=True)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--gtf-sha256", required=True)
    parser.add_argument("--runinfo", type=Path, required=True)
    parser.add_argument("--runinfo-sha256", required=True)
    parser.add_argument("--article-pdf", type=Path, required=True)
    parser.add_argument("--article-pdf-sha256", required=True)
    parser.add_argument("--supplement-pdf", type=Path, required=True)
    parser.add_argument("--supplement-pdf-sha256", required=True)
    parser.add_argument("--article-text", type=Path, required=True)
    parser.add_argument("--comparison-path", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    locked = (
        (args.source / "ARTIFACTS.json", args.source_artifacts_sha256),
        (args.matrix / "ARTIFACTS.json", args.matrix_artifacts_sha256),
        (args.prior_activation / "ARTIFACTS.json", args.prior_activation_artifacts_sha256),
        (args.gtf, args.gtf_sha256),
        (args.runinfo, args.runinfo_sha256),
        (args.article_pdf, args.article_pdf_sha256),
        (args.supplement_pdf, args.supplement_pdf_sha256),
    )
    for path, expected in locked:
        if sha256_file(path) != expected:
            raise GSE274114ActivationError(f"input SHA-256 differs: {path}")

    prior_activation_audit_path = args.prior_activation / "activation_audit.json"
    prior_activation_audit = json.loads(
        prior_activation_audit_path.read_text(encoding="utf-8")
    )
    if prior_activation_audit.get("model_training_activated") is not False:
        raise GSE274114ActivationError(
            "frozen prior activation does not prove training remained blocked"
        )

    article = args.article_text.read_text(encoding="utf-8", errors="strict")
    for phrase in (
        "Table 1. Patient characteristics.",
        "All biopsies were re-scored by an",
        "HiSeq platform",
        "NovaSeq 6000",
        "GSE-274114",
    ):
        if phrase not in article:
            raise GSE274114ActivationError(f"source article evidence absent: {phrase}")

    sample_path = args.source / "samples" / "GSE274114.tsv"
    with sample_path.open(encoding="utf-8", newline="") as handle:
        source_rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(source_rows) != 39:
        raise GSE274114ActivationError("GEO sample count differs")

    topology: list[dict[str, object]] = []
    labels: list[dict[str, str]] = []
    sample_by_experiment: dict[str, dict[str, str]] = {}
    sample_by_gsm: dict[str, dict[str, str]] = {}
    for source_row in source_rows:
        chars = source_characteristics(source_row)
        group = chars["group"]
        if group not in EXPECTED_GROUPS:
            raise GSE274114ActivationError("source group differs")
        biosample, experiment = relation_accessions(source_row)
        participant = chars["sample id"]
        identifier = row_id(participant)
        topology_row = {
            "row_id": identifier,
            "participant_id": participant,
            "sample_accession": source_row["accession"],
            "biosample_accession": biosample,
            "sra_experiment": experiment,
            "biological_unit": "participant",
            "source_biopsy_count": 1,
            "source_material": chars["material"],
            "pairing_topology": "same_study_unpaired",
        }
        topology.append(topology_row)
        semantics, mash_state, hbv_state = GROUP_SEMANTICS[group]
        labels.append(
            {
                "row_id": identifier,
                "source_group": group,
                "source_group_semantics": semantics,
                "mash_state": mash_state,
                "hbv_state": hbv_state,
                "within_instrument_contrast": WITHIN_INSTRUMENT_CONTRAST[group],
                "expected_instrument": GROUP_INSTRUMENT[group],
            }
        )
        sample_by_experiment[experiment] = {**topology_row, "source_group": group}
        sample_by_gsm[source_row["accession"]] = {**topology_row, "source_group": group}

    for field in ("row_id", "participant_id", "sample_accession", "biosample_accession", "sra_experiment"):
        if len({str(row[field]) for row in topology}) != 39:
            raise GSE274114ActivationError(f"participant topology is not one-to-one: {field}")
    if Counter(row["source_group"] for row in labels) != Counter(EXPECTED_GROUPS):
        raise GSE274114ActivationError("source group counts differ")

    with args.runinfo.open(encoding="utf-8", newline="") as handle:
        run_rows = list(csv.DictReader(handle))
    if len(run_rows) != 71:
        raise GSE274114ActivationError("SRA technical run count differs")
    run_join: list[dict[str, object]] = []
    instruments_by_group: dict[str, set[str]] = defaultdict(set)
    for run in run_rows:
        source = sample_by_experiment.get(run["Experiment"])
        if source is None:
            raise GSE274114ActivationError("SRA experiment is outside GEO join")
        if (
            run["BioProject"] != "PRJNA1144975"
            or run["SRAStudy"] != "SRP524555"
            or run["BioSample"] != source["biosample_accession"]
            or run["LibraryName"] != source["sample_accession"]
            or run["LibraryStrategy"] != "RNA-Seq"
            or run["LibraryLayout"] != "PAIRED"
        ):
            raise GSE274114ActivationError("SRA run topology differs")
        instruments_by_group[source["source_group"]].add(run["Model"])
        run_join.append(
            {
                "row_id": source["row_id"],
                "sra_run": run["Run"],
                "sra_experiment": run["Experiment"],
                "sra_sample": run["Sample"],
                "biosample_accession": run["BioSample"],
                "instrument": run["Model"],
                "read_layout": run["LibraryLayout"],
                "spots": run["spots"],
                "bases": run["bases"],
                "technical_run_is_biological_replicate": False,
            }
        )
    if len({row["sra_run"] for row in run_join}) != 71:
        raise GSE274114ActivationError("SRA run axis is not unique")
    run_counts = Counter(str(row["row_id"]) for row in run_join)
    if set(run_counts) != {str(row["row_id"]) for row in topology}:
        raise GSE274114ActivationError("one or more participants lack an SRA run")
    for row in topology:
        row["technical_sra_run_count"] = run_counts[str(row["row_id"])]
    observed_instruments = {
        group: next(iter(values)) if len(values) == 1 else "MULTIPLE"
        for group, values in instruments_by_group.items()
    }
    if observed_instruments != GROUP_INSTRUMENT:
        raise GSE274114ActivationError("group/instrument confounding structure differs")
    overlap_rows = audit_exact_identifier_overlap(
        topology, run_join, args.comparison_path
    )

    quant_paths = sorted((args.matrix / "raw" / "source").glob("*_quant.genes.sf.gz"))
    if len(quant_paths) != 39:
        raise GSE274114ActivationError("Salmon quant file count differs")
    features: list[str] | None = None
    joined_quant_paths: list[tuple[Path, dict[str, str]]] = []
    for quant in quant_paths:
        gsm = quant.name.split("_", 1)[0]
        source = sample_by_gsm.get(gsm)
        if source is None:
            raise GSE274114ActivationError("Salmon filename does not join to GEO")
        current_features: list[str] = []
        fractional_counts = 0
        with gzip.open(quant, "rt", encoding="utf-8", errors="strict") as handle:
            if handle.readline().rstrip("\n").split("\t") != [
                "Name", "Length", "EffectiveLength", "TPM", "NumReads"
            ]:
                raise GSE274114ActivationError("Salmon quant header differs")
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 5:
                    raise GSE274114ActivationError("Salmon quant width differs")
                values = [float(value) for value in fields[1:]]
                if not all(math.isfinite(value) and value >= 0 for value in values):
                    raise GSE274114ActivationError("Salmon quant value is invalid")
                if not values[3].is_integer():
                    fractional_counts += 1
                current_features.append(fields[0])
        if len(current_features) != 61_598 or fractional_counts == 0:
            raise GSE274114ActivationError("Salmon feature/count semantics differ")
        if features is None:
            features = current_features
        elif current_features != features:
            raise GSE274114ActivationError("Salmon feature order differs")
        joined_quant_paths.append((quant, source))
    assert features is not None

    genes = parse_gtf(args.gtf)
    crosswalk = build_crosswalk(features, genes)
    mapping_counts = Counter(row["mapping_state"] for row in crosswalk)
    target_counts = Counter(
        row["gencode_v49_stable_id"]
        for row in crosswalk
        if row["gencode_v49_stable_id"]
    )
    duplicate_targets = {key: value for key, value in target_counts.items() if value > 1}
    if (
        mapping_counts
        != Counter(
            {
                "stable_id_exact": 60_324,
                "par_y_alias_to_v49_canonical": 44,
                "unmapped_stable_id": 1_230,
            }
        )
        or len(target_counts) != 60_324
        or Counter(duplicate_targets.values()) != Counter({2: 44})
    ):
        raise GSE274114ActivationError("GENCODE v40-to-v49 mapping disposition differs")

    model_inputs: list[dict[str, object]] = []
    aggregation_rows: list[dict[str, object]] = []
    quant_dir = args.output / "model_inputs" / "quant_v49"
    quant_dir.mkdir(parents=True, exist_ok=False)
    for quant, source in joined_quant_paths:
        with gzip.open(quant, "rt", encoding="utf-8", errors="strict", newline="") as handle:
            records = list(csv.DictReader(handle, delimiter="\t"))
        aggregated, summary = aggregate_quant_records(records, crosswalk)
        destination = quant_dir / f"{source['row_id']}.quant.genes.sf.gz"
        write_salmon_quant(destination, aggregated)
        model_inputs.append(
            {
                "row_id": source["row_id"],
                "quant_path": destination.relative_to(args.output).as_posix(),
                "quant_sha256": sha256_file(destination),
                "fields": "Name;Length;EffectiveLength;TPM;NumReads",
                "gene_axis": "GENCODE_v49_60324_unique_targets",
                "source_annotation": "GENCODE_v40",
                "target_annotation": "GENCODE_v49",
                "par_y_x_aggregation_applied": True,
                "unmapped_source_rows_masked": 1_230,
                "normalization_or_fit_applied": False,
                "labels_present": False,
            }
        )
        aggregation_rows.append(
            {
                "row_id": source["row_id"],
                "source_quant_sha256": sha256_file(quant),
                "output_quant_sha256": sha256_file(destination),
                **summary,
                "normalization_run": False,
                "model_fit_or_scoring_run": False,
            }
        )

    folds = assign_folds(labels)
    fold_rows = [
        {
            "row_id": row["row_id"],
            "outer_fold": folds[row["row_id"]],
            "source_group": row["source_group"],
            "fold_assignment_used_labels": True,
            "available_to_model_input": False,
        }
        for row in labels
    ]

    aggregate_metadata = [
        {"source_group": "CTRL", "participants": 9, "age_median_iqr": "52 (40-56)", "male": "4 (44.4%)", "fibrosis": "F0/F1 5; not determined 4", "steatosis": "not_reported", "ballooning": "not_reported", "lobular_inflammation": "not_reported", "BMI_median_iqr": "not_reported"},
        {"source_group": "ENEG", "participants": 11, "age_median_iqr": "36 (27-41)", "male": "8 (72.7%)", "fibrosis": "F0/F1 7; F2 4", "steatosis": "less than 5% by source text", "ballooning": "not_reported", "lobular_inflammation": "not_reported", "BMI_median_iqr": "25.3 (23.4-28.6)"},
        {"source_group": "NASH", "participants": 10, "age_median_iqr": "46 (36-51)", "male": "5 (50%)", "fibrosis": "F0/F1 10", "steatosis": "5-33% 4; 34-66% 3; >66% 3", "ballooning": "10 (100%)", "lobular_inflammation": "grade 1 8; grade 2 2", "BMI_median_iqr": "28.4 (27.4-31.2)"},
        {"source_group": "ENEG_NASH", "participants": 9, "age_median_iqr": "44 (36-47)", "male": "7 (77.8%)", "fibrosis": "F0/F1 9", "steatosis": "5-33% 2; 34-66% 7", "ballooning": "9 (100%)", "lobular_inflammation": "grade 1 5; grade 2 4", "BMI_median_iqr": "28.4 (26.5-29.7)"},
    ]
    for row in aggregate_metadata:
        row["availability"] = "group_aggregate_only"
        row["participant_input_allowed"] = False

    args.output.mkdir(parents=True, exist_ok=True)
    write_tsv(args.output / "authorities" / "participant_join.tsv", topology)
    write_tsv(args.output / "authorities" / "run_join.tsv", run_join)
    write_tsv(args.output / "authorities" / "quant_aggregation.tsv", aggregation_rows)
    write_tsv(args.output / "authorities" / "exact_identifier_overlap.tsv", overlap_rows)
    write_tsv(args.output / "model_inputs" / "quant_manifest.tsv", model_inputs)
    write_tsv(args.output / "model_inputs" / "gene_crosswalk.tsv", crosswalk)
    write_tsv(args.output / "evaluator_only" / "source_labels.tsv", labels)
    write_tsv(args.output / "evaluator_only" / "participant_folds.tsv", fold_rows)
    write_tsv(args.output / "descriptive_only" / "group_aggregate_metadata.tsv", aggregate_metadata)

    rights = {
        "schema_version": "masld-bench-gse274114-rights-v1",
        "access_tier": "public",
        "article_license": "CC-BY-4.0",
        "article_doi": "10.1016/j.jhep.2024.10.032",
        "article_pdf_sha256": args.article_pdf_sha256,
        "supplement_pdf_sha256": args.supplement_pdf_sha256,
        "deposited_data_explicit_license_detected": False,
        "new_DUA_or_controlled_access_required": False,
        "internal_nonclinical_research_use": True,
        "redistribution": "prohibited_pending_source_specific_review",
        "derivative_weight_release": "requires_model_specific_terms_and_review",
    }
    topology_authority = {
        "schema_version": "masld-bench-gse274114-topology-v1",
        "participants": 39,
        "biopsies": 39,
        "geo_samples": 39,
        "biosamples": 39,
        "sra_experiments": 39,
        "sra_runs": 71,
        "technical_runs_are_biological_replicates": False,
        "pairing_topology": "same_study_unpaired",
        "participant_identity_status": "authoritative_one_patient_one_biopsy_one_experiment",
        "technical_run_count_range": [min(run_counts.values()), max(run_counts.values())],
    }
    labels_authority = {
        "schema_version": "masld-bench-gse274114-labels-v1",
        "source_group_counts": EXPECTED_GROUPS,
        "labels_location": "evaluator_only/source_labels.tsv",
        "group_aggregate_metadata_location": "descriptive_only/group_aggregate_metadata.tsv",
        "participant_level_age_sex_ethnicity_histology_BMI_metabolic_metadata": "unavailable",
        "HBV_only_is_MASLD_negative_control": False,
        "four_class_performance_allowed": False,
        "MASH_vs_non_MASH_performance_allowed": False,
        "supported_contrasts": ["hiseq_CTRL_vs_ENEG", "novaseq_NASH_vs_ENEG_NASH"],
    }
    qc = {
        "schema_version": "masld-bench-gse274114-qc-v1",
        "all_source_quant_files_checked": True,
        "quant_files": 39,
        "features_per_file": 61_598,
        "header": ["Name", "Length", "EffectiveLength", "TPM", "NumReads"],
        "NumReads_semantics": "Salmon_estimated_counts_continuous_no_rounding",
        "TPM_semantics": "relative_abundance_not_count_likelihood",
        "EffectiveLength_semantics": "required_for_count_scale_offset_reconstruction",
        "normalization_run": False,
        "model_fit_or_scoring_run": False,
        "output_v49_quant_files": 39,
        "output_v49_genes_per_file": 60_324,
        "par_y_x_aggregation_applied_before_normalization": True,
        "unmapped_source_rows_masked": 1_230,
    }
    reference = {
        "schema_version": "masld-bench-gse274114-reference-v1",
        "source_genome": "GRCh38",
        "source_annotation": "GENCODE_v40",
        "target_genome": "GRCh38.p14",
        "target_annotation": "GENCODE_v49",
        "source_features": 61_598,
        "mapped_source_rows": 60_368,
        "unique_v49_target_genes": 60_324,
        "par_y_x_two_row_target_genes": 44,
        "unmapped_source_rows": 1_230,
        "duplicate_rule": "sum_NumReads_and_TPM_before_normalization;TPM_weighted_EffectiveLength",
        "zero_TPM_duplicate_fallback": "unweighted_mean_Length_and_EffectiveLength_after_asserting_zero_aggregate_NumReads",
        "aggregation_conservation_evidence": "authorities/quant_aggregation.tsv",
        "rounding_allowed": False,
    }
    metadata = {
        "schema_version": "masld-bench-gse274114-metadata-v1",
        "source": "article_Table_1",
        "level": "source_group_aggregate_only",
        "participant_level_metadata_available": False,
        "participant_input_allowed": False,
        "descriptive_table": "descriptive_only/group_aggregate_metadata.tsv",
        "fields_reviewed": [
            "age",
            "sex",
            "ethnicity",
            "fibrosis",
            "steatosis",
            "ballooning",
            "lobular_inflammation",
            "BMI",
            "diabetes",
            "smoking",
            "alcohol",
            "ALT",
            "AST",
            "HBV_markers",
        ],
        "rule": "Never assign group summaries to participants, impute them as individual covariates, or use them as model inputs.",
    }
    exposure = {
        "schema_version": "masld-bench-gse274114-exposure-v1",
        "project_exposure": "downstream_demo",
        "prior_source_artifacts": [
            args.source_artifacts_sha256,
            args.matrix_artifacts_sha256,
            args.prior_activation_artifacts_sha256,
        ],
        "sealed_eligible": False,
        "exact_accession_bioproject_biosample_experiment_run_overlap_other_active_families": 0,
        "participant_alias_overlap_other_families": "no_exact_local_alias_match_but_join_unresolved_deidentified_local_ids",
        "exact_overlap_evidence": "authorities/exact_identifier_overlap.tsv",
        "expression_fingerprint_overlap": "derivable_not_processed",
        "previous_model_use": "not_detected_in_frozen_prior_activation_artifacts",
        "previous_model_use_evidence": {
            "path": "gse274114-activation-21066117/activation_audit.json",
            "sha256": sha256_file(prior_activation_audit_path),
            "model_training_activated": False,
            "prior_status": prior_activation_audit.get("status"),
        },
        "audit_limit": "Exact identifier comparison and frozen prior activation review do not replace a permissioned cross-cohort expression or genotype fingerprint.",
    }
    instrument = {
        "schema_version": "masld-bench-gse274114-instrument-confounding-v1",
        "group_instrument": GROUP_INSTRUMENT,
        "perfect_MASH_superclass_instrument_confounding": True,
        "four_class_performance_allowed": False,
        "MASH_vs_non_MASH_performance_allowed": False,
        "OOD_accuracy_using_MASH_vs_healthy_allowed": False,
        "supported_within_instrument_contrasts": [
            {"instrument": "Illumina HiSeq 4000", "groups": ["CTRL", "ENEG"]},
            {"instrument": "Illumina NovaSeq 6000", "groups": ["NASH", "ENEG_NASH"]},
        ],
    }
    write_json(args.output / "authorities" / "rights.json", rights)
    write_json(args.output / "authorities" / "topology.json", topology_authority)
    write_json(args.output / "authorities" / "labels.json", labels_authority)
    write_json(args.output / "authorities" / "qc.json", qc)
    write_json(args.output / "authorities" / "reference.json", reference)
    write_json(args.output / "authorities" / "metadata.json", metadata)
    write_json(args.output / "authorities" / "exposure_audit.json", exposure)
    write_json(args.output / "authorities" / "instrument_confounding.json", instrument)
    write_json(
        args.output / "activation_audit.json",
        {
            "schema_version": "masld-bench-gse274114-no-fit-activation-v1",
            "status": "pass_no_fit_fixture_two_within_instrument_contrasts_only",
            "participants": 39,
            "source_group_counts": EXPECTED_GROUPS,
            "sra_runs": 71,
            "model_input_labels_separated": True,
            "group_aggregate_metadata_descriptive_only": True,
            "four_class_performance_allowed": False,
            "MASH_vs_non_MASH_performance_allowed": False,
            "HBV_only_is_MASLD_negative_control": False,
            "normalization_run": False,
            "model_fit_or_scoring_run": False,
            "output_v49_genes_per_participant": 60_324,
            "par_y_x_aggregation_applied": True,
            "unmapped_source_rows_masked": 1_230,
        },
    )
    print(json.dumps({"status": "pass", "participants": 39, "runs": 71}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
