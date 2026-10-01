#!/usr/bin/env python3
"""Independent bounded BAM recount and excluded-template geometry diagnosis."""
import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil

import pysam

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / "GWAS/finemapping/results/alphagenome_campaign/week1-20260915"
PRODUCER = BASE / "data/read-phase-21773997"
FASTA = Path("/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa")
FIXED = [("GSE281367", donor, "chr8:9325848:A:G", "chr8:9326086:A:G") for donor in ("Z01", "Z04", "Z05", "Z08")]
FIXED += [("GSE244832", "D15", "chr19:45326124:G:A", "chr19:45326536:G:C"),
          ("GSE244832", "D17", "chr22:43928847:C:G", "chr22:43928850:C:T")]


def rows(path, delimiter="\t"):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(2 ** 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def emit(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def parse(uid):
    chrom, position, ref, alt = uid.split(":")
    assert len(ref) == len(alt) == 1 and ref != alt and not set(ref + alt) - set("ACGT")
    return chrom, int(position) - 1, ref, alt


def allele_calls(read, variants):
    """Query-position to reference-position mapping, distinct from aligned-pair producer."""
    reference_positions = read.get_reference_positions(full_length=True)
    target_to_query = {position: query for query, position in enumerate(reference_positions) if position is not None}
    calls = [None, None]
    for site, (_, position, ref, alt) in enumerate(variants):
        query = target_to_query.get(position)
        if query is None or read.query_qualities is None or read.query_qualities[query] < 20:
            continue
        base = read.query_sequence[query].upper()
        calls[site] = "0" if base == ref else "1" if base == alt else "other"
    return calls


def basic_exclusions(read):
    excluded = [name for name in ("is_unmapped", "is_secondary", "is_supplementary", "is_qcfail", "is_duplicate") if getattr(read, name)]
    if not read.is_paired or not read.is_proper_pair or read.mate_is_unmapped:
        excluded.append("not_proper_mapped_pair")
    if read.mapping_quality < 30:
        excluded.append("low_MAPQ")
    return excluded


def barcode(read, source):
    if source == "GSE281367":
        value = read.get_tag("CB") if read.has_tag("CB") else ""
        return value if re.fullmatch(r"[ACGTN]+(?:-[0-9]+)?", value) else None
    value = read.query_name.split(":", 1)[0] if ":" in read.query_name else ""
    return value if re.fullmatch("[ACGTN]+", value) else None


def template_key(read):
    return (read.get_tag("RG") if read.has_tag("RG") else "", read.query_name)


def join_calls(records):
    states = []
    for site in range(2):
        observed = {calls[site] for calls in records if calls[site] is not None}
        states.append(next(iter(observed)) if len(observed) == 1 else "conflict" if observed else None)
    return states


def collapse_template_calls(templates):
    """A two-site template must exist before duplicate collapse can count phase."""
    states = join_calls([template["state"] for template in templates])
    complete = {"".join(template["state"]) for template in templates if all(allele in ("0", "1") for allele in template["state"])}
    conflict = "conflict" in states or len(complete) > 1
    haplotype = next(iter(complete)) if len(complete) == 1 and not conflict else None
    return states, haplotype, any(template["same_read"] for template in templates), conflict


def geometry_diagnostic(path, invalid, variants, source):
    """No relaxation of the census; examine at most one mate per excluded template."""
    grouped = defaultdict(list)
    for read in invalid:
        grouped[template_key(read)].append(read)
    if len(grouped) > 100:
        raise ValueError("Excluded-template diagnostic exceeds fixed bound")
    counts, differences, cigars = Counter(), Counter(), Counter()
    with pysam.AlignmentFile(str(path), "rb") as mate_bam:
        for key, excluded in grouped.items():
            representative = excluded[0]
            counts["excluded_templates"] += 1
            for read in excluded:
                left = min(read.reference_start, read.next_reference_start)
                differences[str(read.reference_end - (left + abs(read.template_length)))] += 1
                cigars[read.cigarstring] += 1
            try:
                mate = mate_bam.mate(representative)
            except (ValueError, OSError, StopIteration):
                counts["mate_lookup_unavailable"] += 1
                continue
            counts["indexed_mate_lookups_succeeded"] += 1
            if template_key(mate) != key or mate.is_read1 == representative.is_read1:
                raise ValueError("Indexed mate identity differs")
            if basic_exclusions(mate) or barcode(representative, source) is None or barcode(mate, source) != barcode(representative, source):
                counts["mate_fails_original_quality_or_barcode_filter"] += 1
                continue
            pair = [representative, mate]
            callset = join_calls([allele_calls(read, variants) for read in pair])
            if all(allele in ("0", "1") for allele in callset):
                counts["excluded_templates_callable_at_both_sites_before_geometry_rule"] += 1
            elif any(allele in ("0", "1") for allele in callset):
                counts["excluded_templates_callable_at_only_one_site"] += 1
            else:
                counts["excluded_templates_without_callable_target"] += 1
            observed_span = max(read.reference_end for read in pair) - min(read.reference_start for read in pair)
            counts["observed_alignment_span_exceeds_abs_TLEN"] += observed_span > abs(representative.template_length)
            counts["observed_alignment_span_matches_abs_TLEN"] += observed_span == abs(representative.template_length)
    return {"excluded_alignment_records": len(invalid), "counts": dict(counts),
            "alignment_end_minus_declared_template_end_histogram": dict(differences), "CIGAR_histogram": dict(cigars),
            "scope": "Same quality and allele filters; diagnostic templates only, not replacement molecule counts or phase calls. TLEN disagreement does not by itself establish Cell Ranger semantics.",
            "raw_read_names_barcodes_sequences_exported": False}


def recount(item, deposited, fasta):
    variants = [parse(item["variant1"]), parse(item["variant2"])]
    chrom, low, high = variants[0][0], variants[0][1], variants[1][1]
    path = Path(item["bam"])
    assert path.stat().st_size == deposited["bam_bytes"] and path.stat().st_mtime_ns == deposited["bam_mtime_ns"]
    assert sha256(Path(str(path) + ".bai")) == deposited["index_sha256"]
    for c, position, ref, alt in variants:
        assert c == chrom and fasta.fetch(c, position, position + 1).upper() == ref
    accepted = defaultdict(list)
    filters, invalid = Counter(), []
    with pysam.AlignmentFile(str(path), "rb") as bam:
        assert bam.get_reference_length(chrom) == fasta.get_reference_length(chrom)
        sample_ids = {entry["SM"] for entry in bam.header.to_dict().get("RG", []) if entry.get("SM")}
        expected_ids = {item[key] for key in ("donor", "source_srr", "source_gsm", "source_atac_sample")}
        assert not sample_ids - expected_ids
        for read in bam.fetch(chrom, low, high + 1):
            filters["alignments_fetched"] += 1
            if filters["alignments_fetched"] > 1000000:
                raise ValueError("Bounded interval limit exceeded")
            excluded = basic_exclusions(read)
            filters.update("excluded_" + reason for reason in excluded)
            if excluded:
                continue
            if read.next_reference_id != read.reference_id or read.template_length == 0:
                filters["excluded_invalid_template_geometry"] += 1
                invalid.append(read)
                continue
            cell = barcode(read, item["source"])
            if cell is None:
                filters["excluded_missing_or_unrecognized_barcode"] += 1
                continue
            left = min(read.reference_start, read.next_reference_start)
            right = left + abs(read.template_length)
            if left < 0 or right <= left or read.reference_end > right:
                filters["excluded_invalid_template_geometry"] += 1
                invalid.append(read)
                continue
            calls = allele_calls(read, variants)
            accepted[template_key(read)].append({"geometry": (cell, chrom, left, right), "calls": calls,
                "arm": 1 if read.is_read1 else 2 if read.is_read2 else 0})
            filters["alignments_after_filters"] += 1
    molecules = defaultdict(list)
    for records in accepted.values():
        geometries = {record["geometry"] for record in records}
        arms = [record["arm"] for record in records]
        if len(geometries) != 1 or len(set(arms)) != len(arms):
            filters["excluded_template_identity_or_geometry_conflict"] += 1
            continue
        molecules[next(iter(geometries))].append({"state": join_calls([record["calls"] for record in records]),
            "same_read": any(all(call in ("0", "1") for call in record["calls"]) for record in records)})
    result, families = Counter(), Counter()
    for geometry, templates in molecules.items():
        families[str(len(templates))] += 1
        states, haplotype, same_read, conflict = collapse_template_calls(templates)
        result["unique_barcode_fragment_geometries"] += 1
        if geometry[2] <= low and geometry[3] > high:
            result["geometrically_spans_both_sites_among_fetched_templates"] += 1
        if any(state in ("0", "1") for state in states):
            result["molecules_with_any_callable_target_allele"] += 1
        for site, state in enumerate(states, 1):
            result[f"site_{site}_{state}"] += 1
        if conflict:
            result["molecules_with_call_conflict"] += 1
        if haplotype:
            result["haplotype_" + haplotype] += 1
            result["molecules_directly_cocovering_both"] += 1
            result["same_read_cocoverage" if same_read else "mate_combined_cocoverage"] += 1
    for key in set(result) | set(deposited["counts"]):
        assert result[key] == deposited["counts"].get(key, 0), (item["donor"], key, result[key], deposited["counts"].get(key, 0))
    for key, value in filters.items():
        assert value == deposited["read_filters"].get(key, 0), (item["donor"], key, value)
    assert dict(families) == deposited["templates_per_fragment_geometry"]
    numerator = result["molecules_directly_cocovering_both"]
    for key, denominator in (("phase_information_fraction_any_callable", result["molecules_with_any_callable_target_allele"]),
        ("direct_cocoverage_fraction_geometric_span_among_fetched", result["geometrically_spans_both_sites_among_fetched_templates"])):
        assert deposited[key] == (numerator / denominator if denominator else None)
    answer = {"source": item["source"], "donor": item["donor"], "counts": dict(result),
        "filter_counts_checked": dict(filters), "templates_per_geometry": dict(families), "BAM_sample_identity_state": deposited["BAM_sample_identity_state"]}
    if item["donor"] == "Z04":
        answer["excluded_geometry_diagnostic"] = geometry_diagnostic(path, invalid, variants, item["source"])
    return answer


def self_checks():
    read = pysam.AlignedSegment()
    read.query_sequence = "AG"; read.query_qualities = [40, 40]
    read.reference_start = 100; read.cigartuples = [(0, 2)]
    sites = [("chr1", 100, "A", "T"), ("chr1", 101, "C", "G")]
    assert allele_calls(read, sites) == ["0", "1"]
    read.flag = 16
    assert allele_calls(read, sites) == ["0", "1"]
    for operation in (2, 3):
        read.cigartuples = [(0, 1), (operation, 1), (0, 1)]
        assert allele_calls(read, sites) == ["0", None]
    read.cigartuples = [(0, 2)]; read.query_qualities = [19, 40]
    assert allele_calls(read, sites) == [None, "1"]
    def template(a, b, same=False):
        return {"state": [a, b], "same_read": same}
    assert collapse_template_calls([template("0", None), template(None, "1")])[1] is None
    assert collapse_template_calls([template("0", "1"), template("1", "1")])[1] is None
    assert collapse_template_calls([template("1", "1"), template("1", "1")])[1] == "11"
    return ["independent_query_reference_mapping", "reverse_orientation_unchanged", "deletion_skip_and_base_quality_excluded",
            "separate_templates_cannot_supply_phase", "conflicting_duplicates_excluded", "identical_duplicates_count_once"]


def main(args):
    if not os.environ.get("SLURM_JOB_ID"):
        raise SystemExit("Compute allocation required")
    args.out.mkdir(parents=True, exist_ok=False, mode=0o700)
    shutil.copy2(__file__, args.out / "executed_review_phase_counts.py")
    tests = self_checks()
    receipt = json.loads((PRODUCER / "receipt.json").read_text())
    assert receipt["status"] == "targeted_ATAC_phase_information_census_complete" and receipt["BAMs"] == 6
    for path, checksum in receipt["source_sha256"].items():
        if path.endswith("data_read_phase.py"):
            path = str(PRODUCER / "executed_data_read_phase.py")
        assert sha256(Path(path)) == checksum
    manifest = rows(PRODUCER / "fixed_candidate_manifest.tsv")
    assert [(row["source"], row["donor"], row["variant1"], row["variant2"]) for row in manifest] == FIXED
    source_maps = {source: {row["donor_id"]: row for row in rows(ROOT / f"data/{source}/metadata/donor_pairing.csv", ",")} for source in {x[0] for x in FIXED}}
    historical_path = next(Path(path) for path in receipt["source_sha256"] if path.endswith("f1_read_backed_phase_candidates.tsv"))
    historical = {(row["dataset"], row["site_pair"]): row for row in rows(historical_path)}
    ase_root = ROOT / "GWAS/finemapping/results/alphagenome_atlas/p3-ase-20260909T192926Z/tables"
    source_alleles = {}
    for source, filename in (("GSE281367", "allelic_donor_counts.tsv.gz"), ("GSE244832", "gse244832_allelic_donor_counts.tsv.gz")):
        source_alleles[source] = {(row["donor"], row["uid"]) for row in rows(ase_root / filename)}
    results = []
    with pysam.FastaFile(str(FASTA)) as fasta:
        for item in manifest:
            prior = historical[(item["source"], item["variant1"] + " | " + item["variant2"])]
            assert item["donor"] in prior["donors_het_at_both"].split(",")
            assert all((item["donor"], item[key]) in source_alleles[item["source"]] for key in ("variant1", "variant2"))
            source = source_maps[item["source"]][item["donor"]]
            assert source["atac_srr"] == item["source_srr"] and source["atac_gsm"] == item["source_gsm"]
            deposited = json.loads((PRODUCER / (item["source"] + "_" + item["donor"] + ".json")).read_text())
            assert all(deposited["identity"][key] == item[key] for key in ("source", "donor", "variant1", "variant2", "bam", "source_srr", "source_gsm"))
            results.append(recount(item, deposited, fasta))
    emit(args.out / "checks.json", {"status": "independent_six_BAM_phase_recount_passed", "job_id": os.environ["SLURM_JOB_ID"],
        "checks": tests, "results": results, "source_receipt_sha256": sha256(PRODUCER / "receipt.json"),
        "base_call_implementation": "get_reference_positions(full_length=True); no_producer_helper_import",
        "source_identity_scope": "Exact saved source metadata plus BAM/index identity; absent SM remains path/metadata assignment, not independent genotype truth",
        "denominator_scope": "Same fixed fetched intervals; unsequenced inserts with both mates outside the interval are not counted",
        "scientific_scope": "Read information availability only; no genotype, interaction, mapping-bias correction, disease or GNMT claim",
        "new_genotype_or_raw_downloads": False, "individual_read_names_barcodes_or_sequences_exported": False,
        "pysam_version": pysam.__version__})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    try:
        main(args)
    except Exception as exc:
        if args.out.is_dir():
            emit(args.out / "checks.json", {"status": "independent_phase_recount_failed", "error": repr(exc), "successful_result_claim": False})
        raise
