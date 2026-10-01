"""Assemble complete native featureCounts single-read rows into seventeen donors.

No receiving-count CLI is enabled. A separately authorized caller must verify
RNA measurement and frozen prediction prerequisites before invoking the callable.
The separate checker uses synthetic counts only; no alignment or target fitting.
"""
from dataclasses import dataclass
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import shlex

import numpy as np

ROOT = Path(__file__).resolve().parents[4]
REC = ROOT / "Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z"
ROSTER = REC / "paired_h3k27ac_raw_roster_21997742"
ANNOTATION = REC / "fixed_h3_saf_21998535"
SAF = ANNOTATION / "fixed_h3_regions.saf"
DONORS = ("B1", "B7", "B8", "B9", "B10", "B12", "B15", "B19", "B21", "B22", "B24", "B26", "B36", "B38", "B41", "B46", "B47")
REGIONS = 96460
EXACT_MAX = 2**53
INPUT_BYTE_CAP = 8_000_000_000
PUBLIC_FIELDS = ("run_accession", "donor", "GSM", "experiment_accession", "source_biosample")
GUARDS = {
    ROSTER / "summary.json": "6f835d7114d15f3b44ff056321a3dc2e59466c1c537ff32a5a9ed334374c5c33",
    ROSTER / "raw_runs.tsv": "69ede048c0975d464c71703193fb279bd77c9e018bf40fcbcdcd6c1ea4fe483b",
    ROSTER / "raw_files.tsv": "f91c00fb37387020a624a02a388ff00d440efc7fa5bdac5ef41da2a8b1ae7e78",
    ROSTER / "source_samples.tsv": "6f11fe9457bcecf470f9a2b6e80d97adc4e09e7643e208495df60ed8a506cbc1",
    SAF: "3c0ae95792c0c67b022bddf5a702b5ac84ba04040be24030edc9b22029e3d598",
    ANNOTATION / "region_coordinate_crosswalk.tsv": "0e5342853edc26da7bde478906d1b92841ad771ffc902bd447bbc9e5d91facc6",
    ANNOTATION / "summary.json": "feb1a674d911725ebf344438597881f06a6c4f2e08bbacbd2383cec487130237",
}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def tsv_rows(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


@dataclass(frozen=True)
class CountTable:
    """Explicit binding for native output columns; fields reuse public roster names.

    columns maps exact featureCounts BAM column strings to PUBLIC_FIELDS dicts.
    count_command is the actual featureCounts argv; prefilter_commands maps each
    column to the declared samtools-view argv that produced that column's BAM.
    These declarations do not prove BAM contents or actual execution themselves.
    """
    path: Path
    columns: dict
    count_command: tuple
    prefilter_commands: dict


@dataclass(frozen=True)
class Contract:
    donors: tuple
    regions: tuple
    runs: dict
    saf_path: Path


def load_fixed_contract():
    require(os.environ.get("SLURM_JOB_ID"), "Metadata loading/assembly requires compute")
    for path, expected in GUARDS.items():
        require(sha(path) == expected, "Frozen roster/annotation changed: " + str(path))
    summary = json.loads((ROSTER / "summary.json").read_text())
    require(tuple(summary["selection"]["paired_donors"]) == DONORS, "Fixed donor order changed")
    samples = {(r["GSM"], r["experiment_accession"]): r for r in tsv_rows(ROSTER / "source_samples.tsv")}
    h3 = [r for r in tsv_rows(ROSTER / "raw_runs.tsv") if r["assay"] == "H3K27ac"]
    files = [r for r in tsv_rows(ROSTER / "raw_files.tsv") if r["assay"] == "H3K27ac"]
    require(len(h3) == len(files) == 202, "Complete H3 run/file roster differs")
    run_map, file_map = {}, {}
    for row in files:
        run = row["run_accession"]
        require(run not in file_map and row["file_number"] == "1", "Repeated or non-single H3 file")
        file_map[run] = row
    for row in h3:
        run = row["run_accession"]
        require(run not in run_map and re.fullmatch(r"SRR[0-9]+", run), "Repeated/invalid H3 run")
        require(row["donor"] in DONORS and row["library_layout"] == "SINGLE"
                and row["library_strategy"] == "ChIP-Seq" and row["fastq_file_count"] == "1", "Wrong source assay/layout")
        sample = samples.get((row["GSM"], row["experiment_accession"]))
        require(sample is not None and sample["assay"] == "H3K27ac", "H3 sample/experiment join absent")
        require(all(sample[k] == row[k] for k in PUBLIC_FIELDS if k != "run_accession"), "Source sample/run donor assignment differs")
        require(run in file_map and all(file_map[run][k] == row[k] for k in PUBLIC_FIELDS), "Source file/run assignment differs")
        require(row["source_biosample"] == row["sample_accession"], "Native biosample join differs")
        run_map[run] = {k: row[k] for k in PUBLIC_FIELDS}
    require(set(run_map) == set(file_map) and {r["donor"] for r in run_map.values()} == set(DONORS), "Incomplete fixed donor/run population")
    saf, crosswalk = tsv_rows(SAF), tsv_rows(ANNOTATION / "region_coordinate_crosswalk.tsv")
    require(len(saf) == len(crosswalk) == REGIONS, "Complete fixed region population changed")
    regions = []
    for i, (row, bridge) in enumerate(zip(saf, crosswalk)):
        start, end = int(row["Start"]), int(row["End"])
        require(row["Strand"] == "+" and 1 <= start <= end, "Invalid fixed SAF geometry")
        require(int(bridge["h3k27ac_feature_index"]) == i and bridge["opaque_source_feature_key"] == row["GeneID"]
                and bridge["alignment_contig"] == row["Chr"] and int(bridge["start1"]) == start
                and int(bridge["end1"]) == end and int(bridge["interval_width_bp"]) == end - start + 1,
                "SAF/fixed source region order differs")
        regions.append((row["GeneID"], row["Chr"], str(start), str(end), "+", str(end - start + 1)))
    require(len({r[0] for r in regions}) == REGIONS, "Repeated region identity")
    return Contract(DONORS, tuple(regions), run_map, SAF)


def validate_count_command(command, columns, saf_path, count_path):
    require(command and Path(command[0]).name == "featureCounts", "Declare actual featureCounts argv")
    options, inputs = {}, []
    i = 1
    while i < len(command):
        token = command[i]
        if token == "--primary":
            require(token not in options, "Repeated counting option")
            options[token] = True
            i += 1
        elif token in ("-T", "-a", "-F", "-o", "-s", "-Q"):
            require(token not in options and i + 1 < len(command), "Repeated/missing counting argument")
            options[token] = command[i + 1]
            i += 2
        else:
            require(not token.startswith("-"), "Unsupported counting flag; no -O/-M/fraction/paired/dedup/overlap changes")
            inputs.append(token)
            i += 1
    require(options.get("-F") == "SAF" and options.get("-s") == "0" and options.get("--primary") is True,
            "Require explicit unstranded primary SAF counting")
    require("-a" in options and Path(options["-a"]).resolve() == saf_path.resolve(), "Wrong fixed SAF path")
    require("-o" in options and Path(options["-o"]).resolve() == Path(count_path).resolve(), "Counting output/path differs")
    require(options.get("-Q", "10") == "10", "MAPQ cutoff differs from fixed10 recipe")
    if "-T" in options:
        require(re.fullmatch(r"[1-9][0-9]*", options["-T"]) is not None, "Invalid counting threads")
    require(tuple(inputs) == tuple(columns) and len(set(inputs)) == len(inputs), "Command input BAMs/columns differ")


def validate_prefilter(command, output_bam, run):
    require(len(command) >= 3 and Path(command[0]).name == "samtools" and command[1] == "view", "Declare samtools-view prefilter argv")
    options, positional = {}, []
    i = 2
    while i < len(command):
        token = command[i]
        if token == "-b":
            require(token not in options, "Repeated filter option")
            options[token] = True
            i += 1
        elif token in ("-q", "-F", "-o", "-@"):
            require(token not in options and i + 1 < len(command), "Repeated/missing filter argument")
            options[token] = command[i + 1]
            i += 2
        else:
            require(not token.startswith("-"), "Unsupported primary-read filtering option")
            positional.append(token)
            i += 1
    require(options.get("-b") and options.get("-q") == "10" and "-F" in options
            and int(options["-F"], 0) == 0x904, "Require MAPQ10 and unmapped/secondary/supplementary exclusion only")
    require(options.get("-o") == output_bam and len(positional) == 1, "Prefilter output/one source BAM differs")
    require(re.findall(r"SRR[0-9]+", Path(positional[0]).name) == [run], "Prefilter source BAM/run differs")


def _assemble_tables(contract, tables):
    """Private core used with a small synthetic region axis by the checker."""
    require(tables, "No native count tables")
    paths, observed_runs = set(), set()
    for spec in tables:
        path = Path(spec.path).resolve()
        require(path not in paths and path.is_file(), "Repeated/absent count file")
        paths.add(path)
        require(spec.columns and set(spec.prefilter_commands) == set(spec.columns), "Column/filter bindings differ")
        for column, binding in spec.columns.items():
            require(set(binding) == set(PUBLIC_FIELDS), "Use exact public run/donor/sample fields")
            run = binding["run_accession"]
            require(run in contract.runs and binding == contract.runs[run], "Unknown or misassigned source H3 run")
            require(run not in observed_runs, "Repeated H3 run across native columns")
            require(re.findall(r"SRR[0-9]+", Path(column).name) == [run], "Native BAM column/run accession differs")
            observed_runs.add(run)
            validate_prefilter(spec.prefilter_commands[column], column, run)
        validate_count_command(spec.count_command, spec.columns, contract.saf_path, spec.path)
    require(observed_runs == set(contract.runs), "Missing source H3 runs; no donor/run omission")
    require(sum(p.stat().st_size for p in paths) <= INPUT_BYTE_CAP, "Count-table byte population exceeds8GB bound")
    counts = np.zeros((len(contract.donors), len(contract.regions)), dtype=np.uint64)
    donor_positions = {donor: i for i, donor in enumerate(contract.donors)}
    totals = [0] * len(contract.donors)
    evidence = []
    for spec in tables:
        digest, count_rows = hashlib.sha256(), 0
        with Path(spec.path).open("rb") as handle:
            first = handle.readline(1_000_001)
            digest.update(first)
            require(len(first) <= 1_000_000 and first.startswith(b"# Program:featureCounts v2.1.1;")
                    and b"Command:" in first, "Require native featureCounts2.1.1 command header")
            emitted_command = tuple(shlex.split(first.decode().split("Command:", 1)[1].strip()))
            require(emitted_command == tuple(spec.count_command), "Emitted and declared count command differ")
            header = handle.readline(1_000_001)
            digest.update(header)
            fields = header.decode().rstrip("\r\n").split("\t")
            require(fields[:6] == ["Geneid", "Chr", "Start", "End", "Strand", "Length"]
                    and tuple(fields[6:]) == tuple(spec.columns), "Native count columns differ or duplicate")
            positions = [donor_positions[spec.columns[column]["donor"]] for column in fields[6:]]
            while True:
                line = handle.readline(1_000_001)
                if not line:
                    break
                digest.update(line)
                require(len(line) <= 1_000_000 and count_rows < len(contract.regions), "Extra/unbounded count row")
                row = line.decode().rstrip("\r\n").split("\t")
                require(len(row) == len(fields) and tuple(row[:6]) == contract.regions[count_rows],
                        "Missing/repeated/misordered region or altered Chr/Start/End/Strand/Length")
                for donor, text in zip(positions, row[6:]):
                    require(re.fullmatch(r"[0-9]+", text) is not None and len(text) <= 16, "Count is not a native nonnegative integer")
                    value = int(text)
                    combined = int(counts[donor, count_rows]) + value
                    totals[donor] += value
                    require(value <= EXACT_MAX and combined <= EXACT_MAX and totals[donor] <= EXACT_MAX,
                            "Integer/count-total exceeds exact float64 downstream range")
                    counts[donor, count_rows] = combined
                count_rows += 1
        require(count_rows == len(contract.regions), "Missing native region rows; no zero fill")
        evidence.append(dict(path=str(Path(spec.path).resolve()), sha256=digest.hexdigest(),
                             rows=count_rows, native_column_runs=[spec.columns[c]["run_accession"] for c in spec.columns],
                             count_command=list(spec.count_command)))
    require(all(total > 0 for total in totals), "Every fixed donor needs a positive complete-region count total")
    require(counts.shape == (len(contract.donors), len(contract.regions)), "Aggregate axis changed")
    return dict(counts=counts, donor_ids=contract.donors, region_keys=tuple(r[0] for r in contract.regions),
                evidence=dict(input_tables=evidence, run_count=len(observed_runs), biological_n=len(contract.donors),
                              count_unit="primary nonsupplementary single reads; MAPQ10; unstranded; no -O/-M/fraction/dedup",
                              zero_rows_retained_only_from_complete_native_tables=True, missing_rows_created=False,
                              source_run_grouping_not_private_identity_verification=True,
                              bam_contents_or_command_execution_verified_by_assembler=False,
                              normalization_or_target_fitting=False, biological_validation=False))


def assemble_fixed_h3_counts(tables):
    """Authorized production callable; caller verifies RNA/prediction prerequisites.

    Reads actual native counts only when invoked on compute. Source assignments
    and declarations must also be proven by the caller's alignment/count receipts.
    No native input row is invented; all202 runs and96460 rows are mandatory.
    """
    require(os.environ.get("SLURM_JOB_ID"), "Native count assembly requires compute")
    contract = load_fixed_contract()
    result = _assemble_tables(contract, tables)
    for path, expected in GUARDS.items():
        require(sha(path) == expected, "Frozen public metadata changed during assembly")
    result["evidence"]["public_metadata_sha256"] = {str(p):v for p,v in GUARDS.items()}
    return result


if __name__ == "__main__":
    raise SystemExit("No receiving-count CLI enabled. Use the synthetic checker; an authorized caller must verify RNA/frozen-prediction prerequisites before count assembly.")
