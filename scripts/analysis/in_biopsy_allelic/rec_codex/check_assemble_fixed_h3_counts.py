"""Synthetic-only assembler checks; public run metadata, no receiving counts."""
from dataclasses import replace
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys

import numpy as np

from assemble_fixed_h3_counts import (ANNOTATION, Contract, CountTable, GUARDS, REC,
                                      _assemble_tables, load_fixed_contract, require, sha)


def actual_tool_smoke(fixed, output):
    """Native output from planted SAM records, never study reads or alignment."""
    binaries = Path("/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin")
    samtools, featurecounts = binaries / "samtools", binaries / "featureCounts"
    inventory = {}
    for name, binary, version_arg, expected in [("samtools", samtools, "--version", "samtools 1.22.1"),
                                               ("featureCounts", featurecounts, "-v", "featureCounts v2.1.1")]:
        version = subprocess.run([str(binary), version_arg], capture_output=True, text=True, check=True, timeout=60)
        require(expected in version.stdout + version.stderr, "Actual-tool smoke version changed: " + name)
        inventory[name] = dict(path=str(binary), sha256=sha(binary), version_stdout=version.stdout,
                               version_stderr=version.stderr, version_argv=[str(binary), version_arg])
    run = next(iter(fixed.runs))
    binding = dict(fixed.runs[run])
    regions = fixed.regions[:3]
    first = regions[0]
    start, end = int(first[2]), int(first[3])
    require(end - start + 1 >= 20, "Synthetic20M read exceeds first fixed interval")
    require(all(r[1] != first[1] or int(r[3]) < start or int(r[2]) > start + 19 for r in regions[1:]),
            "Planted first-region read also overlaps a different toy interval")
    directory = output / "actual_tool_smoke"
    directory.mkdir()
    saf = directory / "toy_fixed_regions.saf"
    with saf.open("x") as handle:
        handle.write("GeneID\tChr\tStart\tEnd\tStrand\n")
        for row in regions:
            handle.write("\t".join(row[:5]) + "\n")
    source = directory / (run + ".aligned.sam")
    # Both primary records count; the duplicate flag is retained by this recipe.
    records = [("accepted_primary",0,20), ("accepted_duplicate_flag",1024,20),
               ("excluded_mapq9",0,9), ("excluded_secondary",256,20), ("excluded_supplementary",2048,20)]
    with source.open("x") as handle:
        handle.write("@HD\tVN:1.6\tSO:unsorted\n")
        for chromosome in dict.fromkeys(row[1] for row in regions):
            length = max(int(row[3]) for row in regions if row[1] == chromosome) + 1000
            handle.write(f"@SQ\tSN:{chromosome}\tLN:{length}\n")
        for name,flag,mapq in records:
            handle.write("\t".join(map(str,(name,flag,first[1],start,mapq,"20M","*",0,0,"A"*20,"I"*20))) + "\n")
    bam = directory / (run + ".primary.bam")
    filter_command = (str(samtools),"view","-b","-q","10","-F","2308","-o",str(bam),str(source))
    native = directory / "native_featurecounts.tsv"
    count_command = (str(featurecounts),"-T","1","-F","SAF","-a",str(saf),"-o",str(native),
                     "-s","0","-Q","10","--primary",str(bam))
    commands = []
    for name, command in [("samtools_filter",filter_command), ("featureCounts",count_command)]:
        completed = subprocess.run(command, capture_output=True, text=True, timeout=60)
        (directory / (name + ".stdout")).write_text(completed.stdout)
        (directory / (name + ".stderr")).write_text(completed.stderr)
        completed.check_returncode()
        commands.append(list(command))
    inspect_command = [str(samtools),"view",str(bam)]
    inspected = subprocess.run(inspect_command, capture_output=True, text=True, check=True, timeout=60)
    observed = [(line.split("\t")[0],int(line.split("\t")[1])) for line in inspected.stdout.splitlines()]
    require(observed == [("accepted_primary",0),("accepted_duplicate_flag",1024)], "Actual MAPQ/secondary/supplementary filter differs")
    spec = CountTable(native,{str(bam):binding},count_command,{str(bam):filter_command})
    toy = Contract((binding["donor"],),regions,{run:binding},saf)
    result = _assemble_tables(toy,[spec])
    require(np.array_equal(result["counts"],np.array([[2,0,0]],dtype=np.uint64)),
            "Native output differs from hand-derived two primary single-read counts and two zero rows")
    with native.open() as handle:
        native_command_header = handle.readline().rstrip("\n")
    receipt = dict(synthetic_only=True, no_genome_alignment_performed=True, run_id_used_as_public_schema_label=run,
                   planted_records=[dict(name=name,flag=flag,mapq=mapq) for name,flag,mapq in records],
                   expected_counts=[[2,0,0]], actual_counts=result["counts"].tolist(), actual_primary_records=observed,
                   native_command_header=native_command_header, actual_argv=commands,
                   samtools_inspect_argv=inspect_command, tool_inventory=inventory,
                   files_sha256={p.name:sha(p) for p in directory.iterdir() if p.is_file()},
                   receiving_counts_reads_or_alignments_opened=False, biological_validation=False)
    (directory / "receipt.json").write_text(json.dumps(receipt,indent=2) + "\n")
    return receipt


def main():
    require(os.environ.get("SLURM_JOB_ID"), "Synthetic assembly checks require compute")
    fixed = load_fixed_contract()
    # Three real coordinate metadata rows; synthetic numbers on all real202 run IDs.
    toy = Contract(fixed.donors, fixed.regions[:3], fixed.runs, fixed.saf_path)
    output = REC / ("fixed_h3_assembly_checks_" + os.environ["SLURM_JOB_ID"])
    output.mkdir(exist_ok=False)
    tool_smoke = actual_tool_smoke(fixed,output)
    runs = tuple(toy.runs)
    numbers = {run: [i + 1, 2 * (i + 1), 0] for i, run in enumerate(runs)}

    def emit(name, selected=runs, values=None, rows=None, command_extra=()):
        path = output / (name + ".counts.tsv")
        columns = {"/synthetic/" + run + ".primary.bam": dict(toy.runs[run]) for run in selected}
        prefilters = {column: ("samtools", "view", "-b", "-q", "10", "-F", "2308", "-o", column,
                              "/synthetic/" + binding["run_accession"] + ".aligned.bam") for column,binding in columns.items()}
        command = ("featureCounts", "-T", "1", "-F", "SAF", "-a", str(toy.saf_path), "-o", str(path),
                   "-s", "0", "-Q", "10", "--primary") + tuple(command_extra) + tuple(columns)
        values = numbers if values is None else values
        rows = toy.regions if rows is None else rows
        with path.open("x") as handle:
            handle.write("# Program:featureCounts v2.1.1; Command:" + shlex.join(command) + "\n")
            handle.write("\t".join(("Geneid", "Chr", "Start", "End", "Strand", "Length") + tuple(columns)) + "\n")
            for i, row in enumerate(rows):
                handle.write("\t".join(tuple(row) + tuple(str(values[run][min(i,2)]) for run in selected)) + "\n")
        return CountTable(path, columns, command, prefilters)

    all_table = emit("valid_combined")
    result = _assemble_tables(toy, [all_table])
    expected = np.zeros((17, 3), dtype=np.uint64)
    # Independent Python-integer donor sums, not a copy of the streaming parser.
    for j, donor in enumerate(toy.donors):
        weight = sum(i + 1 for i,run in enumerate(runs) if toy.runs[run]["donor"] == donor)
        expected[j] = [weight, 2 * weight, 0]
    require(np.array_equal(result["counts"], expected), "Hand-derived exact donor sums disagree")
    require(result["donor_ids"] == fixed.donors and result["region_keys"] == tuple(r[0] for r in toy.regions)
            and np.all(result["counts"][:,2] == 0), "Native zeros or fixed axes changed")
    split = [emit("valid_part1", runs[:101]), emit("valid_part2", runs[101:])]
    partition = _assemble_tables(toy, split[::-1])
    require(np.array_equal(partition["counts"], expected), "Technical file partition/order changes exact donor sums")
    # The public callable's fixed metadata load must still retain the full axis.
    require(len(fixed.regions) == 96460 and len(fixed.runs) == 202 and len(fixed.donors) == 17,
            "Synthetic core replaced the production metadata population")
    cases = []
    cases.append(("missing_run", [emit("missing_run", runs[:-1])]))
    cases.append(("duplicate_run", [all_table, emit("duplicate_full_run_set")]))
    cases.append(("repeated_file", [all_table, all_table]))
    first_column = next(iter(all_table.columns))
    for key,value in [("run_accession", "SRR0000000"), ("donor", "B7"), ("GSM", "GSM0000000"),
                      ("experiment_accession", "SRX0000000"), ("source_biosample", "SAMN0000000")]:
        columns = {c: dict(binding) for c,binding in all_table.columns.items()}
        columns[first_column][key] = value
        cases.append(("wrong_" + key, [replace(all_table, columns=columns)]))
    first_run = all_table.columns[first_column]["run_accession"]
    edge = {r:list(v) for r,v in numbers.items()}
    old_total = sum(sum(values) for run,values in numbers.items() if toy.runs[run]["donor"] == "B1")
    edge[first_run][0] = 2**53 - (old_total - numbers[first_run][0])
    edge_result = _assemble_tables(toy, [emit("valid_exact_total_boundary", values=edge)])
    require(sum(map(int, edge_result["counts"][0])) == 2**53
            and int(edge_result["counts"][0,0]) == sum(values[0] for run,values in edge.items() if toy.runs[run]["donor"] == "B1"),
            "Valid exact-integer boundary rounded or overflowed")
    other_run = runs[-1]
    changed_columns = {c: dict(binding) for c,binding in all_table.columns.items()}
    changed_columns[first_column] = dict(toy.runs[other_run])
    cases.append(("column_run_mismatch", [replace(all_table, columns=changed_columns)]))
    for flag in ("-O", "-M", "--fraction", "-p", "--countReadPairs", "--ignoreDup", "--largestOverlap"):
        cases.append(("forbidden_" + flag, [emit("forbidden_" + flag.replace("-", ""), command_extra=(flag,))]))
    for label, bad_rows in [("missing_region", toy.regions[:-1]), ("extra_region", toy.regions + (toy.regions[-1],)),
                            ("duplicate_region", (toy.regions[0], toy.regions[0], toy.regions[2])),
                            ("reordered_region", toy.regions[::-1])]:
        cases.append((label, [emit(label, rows=bad_rows)]))
    for index,value in [(1,"2"), (2,str(int(toy.regions[0][2])+1)), (3,str(int(toy.regions[0][3])+1)), (4,"-"), (5,"1")]:
        row = list(toy.regions[0]); row[index] = value
        cases.append(("wrong_geometry_" + str(index), [emit("wrong_geometry_" + str(index), rows=(tuple(row),) + toy.regions[1:])]))
    for label,value in [("fractional_count","0.5"), ("negative_count","-1"), ("exponent_count","1e3"),
                        ("nonfinite_count","NaN"), ("inexact_count",2**53+1)]:
        values = {r:list(v) for r,v in numbers.items()}; values[first_run][0] = value
        cases.append((label, [emit(label, values=values)]))
    values = {r:list(v) for r,v in numbers.items()}; values[first_run][0] = 2**53
    cases.append(("exact_total_overflow", [emit("exact_total_overflow", values=values)]))
    values = {r:([0,0,0] if toy.runs[r]["donor"] == "B1" else list(v)) for r,v in numbers.items()}
    cases.append(("zero_donor_total", [emit("zero_donor_total", values=values)]))
    prefilters = dict(all_table.prefilter_commands)
    for label,index,value in [("supplementary_not_filtered",6,"260"), ("wrong_mapq",4,"0"),
                              ("wrong_filter_run",9,"/synthetic/SRR0000000.aligned.bam")]:
        command = list(all_table.prefilter_commands[first_column]); command[index] = value
        changed = dict(prefilters); changed[first_column] = tuple(command)
        cases.append((label, [replace(all_table, prefilter_commands=changed)]))
    # A truthful but changed native command is refused independently of count values.
    cmd = list(all_table.count_command); cmd[cmd.index("-s")+1] = "1"
    cases.append(("stranded_command", [replace(all_table, count_command=tuple(cmd))]))
    cmd = list(all_table.count_command); cmd[cmd.index("-T")+1] = "2"
    cases.append(("native_command_disagreement", [replace(all_table, count_command=tuple(cmd))]))
    wrong_version = emit("native_version_disagreement")
    wrong_version.path.write_text(wrong_version.path.read_text().replace("v2.1.1;", "v2.1.0;", 1))
    cases.append(("native_version_disagreement", [wrong_version]))
    refusals = []
    for label,specs in cases:
        try:
            _assemble_tables(toy, specs)
        except ValueError as error:
            refusals.append(dict(case=label, refusal=str(error)))
        else:
            raise AssertionError("Invalid native input was accepted: " + label)
    for path,expected_hash in GUARDS.items():
        require(sha(path) == expected_hash, "Public metadata changed during checks")
    require(sum(p.stat().st_size for p in output.rglob("*") if p.is_file()) <= 50_000_000,
            "Synthetic output byte cap exceeded")
    summary = dict(synthetic_only=True, synthetic_regions=3, production_regions=len(fixed.regions),
                   source_h3_runs=len(fixed.runs), fixed_donors=list(fixed.donors),
                   hand_derived_donor_sums_agree=True, technical_partition_and_file_order_invariant=True,
                   exact_integer_total_boundary_retained=True,
                   native_all_zero_region_retained=True, invalid_input_refusals=refusals,
                   actual_tool_format_and_single_read_smoke=tool_smoke,
                   script_sha256=sha(Path(__file__)), assembler_sha256=sha(Path(__file__).with_name("assemble_fixed_h3_counts.py")),
                   launcher_sha256=sha(Path(__file__).with_name("run_check_assemble_fixed_h3_counts.sbatch")),
                   public_metadata_sha256={str(p):v for p,v in GUARDS.items()},
                   receiving_counts_or_reads_opened=False, targets_or_predictions_evaluated=False,
                   biological_validation=False, stochastic_operations=False,
                   python=sys.version, numpy=np.__version__, slurm_job_id=os.environ["SLURM_JOB_ID"])
    (output / "executed_checker.py").write_bytes(Path(__file__).read_bytes())
    (output / "executed_assembler.py").write_bytes(Path(__file__).with_name("assemble_fixed_h3_counts.py").read_bytes())
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(dict(synthetic_only=True, hand_derived_donor_sums_agree=True,
                         technical_partition_invariant=True, invalid_refusals=len(refusals))))


if __name__ == "__main__":
    main()
