"""Build a fixed, approximate Ensembl98 RNA index; no receiving measurements.

Complete native transcript estimates must later be verified before gene
aggregation. Index presence alone is not measurement or predictive validation.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

from export_training_h3_target_transform import require, sha

ROOT = Path(__file__).resolve().parents[4]
REC = ROOT / "Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z"
TX_DIR = REC / "ensembl98_transcript_reconstruction_21998470"
GENOME_DIR = REC / "ensembl98_primary_genome_21998416"
TX = TX_DIR / "ensembl98_gtf_complete_transcripts.fa"
GENOME = GENOME_DIR / "Homo_sapiens.GRCh38.dna.primary_assembly.fa"
FAI = Path(str(GENOME) + ".fai")
SALMON = Path("/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/salmon")
INPUTS = {
    TX_DIR / "summary.json": "1e7d48b28c1178523df3bbc455a6b9cf2c7786124878bf4743566aeb9469e144",
    TX_DIR / "transcript_gene_axis.tsv": "0ac3dd2c4cc9b6923afea3041f5d4e94b67b9ded0f33acc4c028e6f3845eb8d2",
    TX: "cf3fac2b075b649c5bf8c80ef2ec860f068f63d92357560fbcb3f1c3d12d38e2",
    GENOME_DIR / "summary.json": "9f93a8291cf676bbf1bd1bc552d2097b3e0f9850a05f86431e1d93df65392f6b",
    GENOME: "78777b0886e8dfa5e14e4957fbbaa53736fcbaa5668d59e09b6b7945fca93d8c",
    FAI: "57f6ed6f8b07437708fff09814489cfda4b976652afabe2c9b1d54463bcba0ad",
}
DISK_CAP = 300_000_000_000
GENTROME_CAP = 4_000_000_000
THREADS = 8


def fasta_headers(path):
    ids = []
    with path.open("rb") as handle:
        for line in handle:
            if line.startswith(b">"):
                ids.append(line[1:].split()[0].decode("ascii"))
    require(len(ids) == len(set(ids)), "Duplicate FASTA identifiers: " + str(path))
    return ids


def output_bytes(path):
    total = 0
    for p in path.rglob("*"):
        try:
            if p.is_file():
                total += p.stat().st_size
        except FileNotFoundError:
            # Index construction may remove a temporary file during inspection.
            continue
    return total


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    require(os.environ.get("SLURM_JOB_ID"), "Index construction requires compute")
    require(not args.output.exists(), "Refusing output overwrite")
    require(int(os.environ["SLURM_CPUS_PER_TASK"]) >= THREADS, "Insufficient requested CPUs")
    for path, expected in INPUTS.items():
        require(sha(path) == expected, "Frozen reference changed: " + str(path))
    require(json.loads((TX_DIR / "summary.json").read_text())["reference_identity_qualified"],
            "Reference reconstruction is not qualified")
    version = subprocess.run([str(SALMON), "--version"], capture_output=True, text=True, check=True)
    require(version.stdout.strip() == "salmon 1.10.3", "Installed Salmon version changed")
    transcripts = fasta_headers(TX)
    decoys = fasta_headers(GENOME)
    require(len(transcripts) == 227368 and len(decoys) == 194, "Reference header population changed")
    require(not set(transcripts).intersection(decoys), "Transcript/decoy ID collision")
    with (TX_DIR / "transcript_gene_axis.tsv").open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    require(set(transcripts) == {r["transcript_id"] for r in rows}, "Native transcript metadata join differs")
    fai_ids = [line.split("\t")[0] for line in FAI.read_text().splitlines()]
    require(decoys == fai_ids, "Genome sequence and declared contig order differ")
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / "decoys.txt").write_text("\n".join(decoys) + "\n")
    gentrome = args.output / "gentrome.fa"
    combined, written = hashlib.sha256(), 0
    with gentrome.open("xb") as dest:
        for source in (TX, GENOME):
            copied_sha = hashlib.sha256()
            with source.open("rb") as src:
                last = b""
                for chunk in iter(lambda: src.read(1 << 20), b""):
                    written += len(chunk)
                    require(written <= GENTROME_CAP, "Gentrome cap exceeded; partial preserved")
                    dest.write(chunk)
                    combined.update(chunk)
                    copied_sha.update(chunk)
                    last = chunk[-1:]
                require(last == b"\n", "FASTA source lacks terminal newline")
            require(copied_sha.hexdigest() == INPUTS[source], "Reference changed during copying")
    require(written == 3524663123, "Unexpected full reference FASTA length")
    index = args.output / "index"
    command = ["/usr/bin/time", "-v", "-o", str(args.output / "resource_usage.txt"),
               str(SALMON), "index", "-t", str(gentrome), "-i", str(index),
               "-d", str(args.output / "decoys.txt"), "-k", "31", "-p", str(THREADS),
               "--keepDuplicates", "--no-clip", "--keepFixedFasta"]
    protocol = dict(transcripts=227368, decoys=194, k=31, index_type="puff",
                    transcript_order="native complete reconstructed transcripts before all verified genomic decoys",
                    native_gene_estimate_contract="all native transcript estimate rows required, including declared structural zeros; exact versioned tx-to-native-gene mapping before modeled42163-axis extraction",
                    keep_duplicates=True, poly_a_clipping=False,
                    reason_for_no_clip="preserve native annotation-derived terminal sequence; declared departure from Salmon default, not private pipeline reconstruction",
                    identical_sequence_individual_gene_identifiability_asserted=False,
                    receiving_quantification_plan=dict(library_type="A, automatically inferred from receiving RNA mappings before mate arguments", bias_flags=[],
                                                     gene_aggregation="tximport countsFromAbundance=no, txOut=TRUE then exact summarizeToGene; no version/bar stripping or silently dropped unmapped transcripts"),
                    input_sha256={str(p): h for p, h in INPUTS.items()},
                    gentrome_bytes=written, gentrome_sha256=combined.hexdigest(),
                    saved_output_cap_bytes=DISK_CAP, cap_monitor_interval_seconds=30,
                    command=command, salmon_version=version.stdout.strip(),
                    source_sha256=sha(Path(__file__)),
                    launcher_sha256=sha(Path(__file__).with_name("run_codex_build_native_rna_index.sbatch")),
                    python=sys.version, slurm_job_id=os.environ["SLURM_JOB_ID"],
                    source_pipeline_equivalence=False, receiving_data_read=False,
                    model_fitted=False, biological_observability_asserted=False)
    (args.output / "protocol.json").write_text(json.dumps(protocol, indent=2) + "\n")
    started = time.monotonic()
    with (args.output / "salmon_index.log").open("x") as log:
        child = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
        try:
            while child.poll() is None:
                size = output_bytes(args.output)
                print(json.dumps(dict(index_running=True, output_bytes=size, elapsed_seconds=time.monotonic()-started)), flush=True)
                require(size <= DISK_CAP, "Own index output cap exceeded; partials preserved")
                try:
                    child.wait(timeout=30)
                except subprocess.TimeoutExpired:
                    pass
        except BaseException:
            os.killpg(child.pid, signal.SIGTERM)
            try:
                child.wait(timeout=20)
            except subprocess.TimeoutExpired:
                os.killpg(child.pid, signal.SIGKILL)
                child.wait()
            raise
    result = dict(index_returncode=child.returncode, elapsed_seconds=time.monotonic()-started,
                  output_bytes=output_bytes(args.output), protocol_sha256=sha(args.output / "protocol.json"),
                  index_built=False, native_transcript_quantification_not_run=True,
                  receiving_data_read=False, predictive_validation=False)
    (args.output / "summary.json").write_text(json.dumps(result, indent=2) + "\n")
    require(child.returncode == 0, "Salmon index failed; logs and partials preserved")
    require(result["output_bytes"] <= DISK_CAP, "Completed index exceeds declared allocation")
    require((index / "versionInfo.json").is_file() and (index / "info.json").is_file(),
            "Salmon success lacks expected index metadata")
    result["index_metadata"] = {name: json.loads((index / name).read_text())
                                for name in ("versionInfo.json", "info.json")}
    result["index_metadata_sha256"] = {name: sha(index / name)
                                       for name in ("versionInfo.json", "info.json")}
    result["index_built"] = True
    (args.output / "summary.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
