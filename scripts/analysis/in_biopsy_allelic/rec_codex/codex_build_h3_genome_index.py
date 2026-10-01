"""Index only the fixed public genome; no receiving chromatin measurements.

Reconstruct indexed references and compare canonical sequence hashes against
the exact input genome, with Bowtie2's declared ambiguous-base-to-N handling.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys

from codex_build_native_rna_index import output_bytes
from export_training_h3_target_transform import require, sha

ROOT = Path(__file__).resolve().parents[4]
REC = ROOT / "Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z"
REFERENCE = REC / "ensembl98_primary_genome_21998416"
FASTA = REFERENCE / "Homo_sapiens.GRCh38.dna.primary_assembly.fa"
FAI = Path(str(FASTA) + ".fai")
GUARDS = {
    FASTA: "78777b0886e8dfa5e14e4957fbbaa53736fcbaa5668d59e09b6b7945fca93d8c",
    FAI: "57f6ed6f8b07437708fff09814489cfda4b976652afabe2c9b1d54463bcba0ad",
    REFERENCE / "summary.json": "9f93a8291cf676bbf1bd1bc552d2097b3e0f9850a05f86431e1d93df65392f6b",
}
CAP = 50_000_000_000
THREADS = 8
SEED = 20260930
TRANSLATION = bytes(i if i in b"ACGT" else ord("N") for i in range(256))


def canonical_sequences(handle):
    records = {}
    name, digest, length, ambiguous, total = None, None, 0, 0, 0

    def finish():
        require(name is not None and length > 0 and name not in records,
                "Empty or repeated reference sequence")
        records[name] = dict(length=length, canonical_sha256=digest.hexdigest(),
                             ambiguous_bases=ambiguous)

    for line in handle:
        require(len(line) <= 1_000_000, "Unexpected unbounded FASTA line")
        if line.startswith(b">"):
            if name is not None:
                finish()
            name = line[1:].split()[0].decode("ascii")
            digest, length, ambiguous = hashlib.sha256(), 0, 0
        else:
            seq = b"".join(line.split()).upper()
            if not seq:
                continue
            require(name is not None, "Sequence before reference name")
            length += len(seq)
            total += len(seq)
            require(total <= 3_300_000_000, "Reference sequence population exceeds expected bounds")
            ambiguous += len(seq) - sum(seq.count(base) for base in (b"A", b"C", b"G", b"T"))
            digest.update(seq.translate(TRANSLATION))
    require(name is not None, "Empty reference FASTA")
    finish()
    return records


def run_index(command, output):
    with (output / "bowtie2_build.log").open("x") as log:
        child = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
        try:
            while child.poll() is None:
                size = output_bytes(output)
                require(size <= CAP, "Own index output cap exceeded; partial preserved")
                print(json.dumps(dict(index_running=True, saved_output_bytes=size)), flush=True)
                try:
                    child.wait(timeout=30)
                except subprocess.TimeoutExpired:
                    pass
            require(child.returncode == 0, "Genome index build failed; partial preserved")
        except BaseException:
            if child.poll() is None:
                os.killpg(child.pid, signal.SIGTERM)
                try:
                    child.wait(timeout=20)
                except subprocess.TimeoutExpired:
                    os.killpg(child.pid, signal.SIGKILL)
                    child.wait()
            raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    require(os.environ.get("SLURM_JOB_ID"), "Public genome indexing requires compute")
    require(int(os.environ["SLURM_CPUS_PER_TASK"]) >= THREADS, "Insufficient CPUs")
    require(args.output.resolve().parent == REC.resolve() and not args.output.exists(),
            "Require a new output directory under owned execution root")
    for path, expected in GUARDS.items():
        require(sha(path) == expected, "Frozen genome changed: " + str(path))
    builder, inspector = shutil.which("bowtie2-build"), shutil.which("bowtie2-inspect")
    require(builder and inspector, "Load Bowtie2/2.5.4-linux-x86_64")
    version = subprocess.run([builder, "--version"], capture_output=True, text=True, check=True)
    require(" version 2.5.4" in version.stdout, "Installed Bowtie2 version differs")
    expected_lengths = {p[0]: int(p[1]) for p in (line.split("\t") for line in FAI.read_text().splitlines())}
    require(len(expected_lengths) == 194 and sum(expected_lengths.values()) == 3099750718,
            "Fixed genome contig population differs")
    with FASTA.open("rb") as handle:
        expected_sequences = canonical_sequences(handle)
    require({name: row["length"] for name, row in expected_sequences.items()} == expected_lengths,
            "Genome FASTA and fixed FAI differ")
    args.output.mkdir(exist_ok=False)
    prefix = args.output / "GRCh38_ensembl98_primary"
    command = ["/usr/bin/time", "-v", "-o", str(args.output / "resource_usage.txt"),
               builder, "--threads", str(THREADS), "--seed", str(SEED), str(FASTA), str(prefix)]
    receipt = dict(command=command, bowtie2_version=version.stdout, module="Bowtie2/2.5.4-linux-x86_64",
                   build_seed=SEED, contigs=194, reference_bases=3099750718,
                   input_sha256={str(path): value for path, value in GUARDS.items()},
                   tool_paths=dict(builder=builder, inspector=inspector),
                   wrappers_sha256=dict(builder=sha(Path(builder)), inspector=sha(Path(inspector))),
                   expected_canonical_sequences=expected_sequences,
                   saved_output_cap_bytes=CAP, cap_monitor_seconds=30,
                   source_sha256=sha(Path(__file__)), python=sys.version,
                   launcher_sha256=sha(Path(__file__).with_name("run_codex_build_h3_genome_index.sbatch")),
                   slurm_job_id=os.environ["SLURM_JOB_ID"], receiving_data_read=False,
                   alignment_or_counting_run=False, predictive_validation=False,
                   interpretation="Matching public reference only; no reconstruction of source CUT&RUN or ChIP measurement")
    (args.output / "protocol.json").write_text(json.dumps(receipt, indent=2) + "\n")
    run_index(command, args.output)
    suffixes = (".1", ".2", ".3", ".4", ".rev.1", ".rev.2")
    families = [[Path(str(prefix) + suffix + ext) for suffix in suffixes] for ext in (".bt2", ".bt2l")]
    complete = [family for family in families if all(path.is_file() and path.stat().st_size > 0 for path in family)]
    require(len(complete) == 1 and not any(path.exists() for family in families if family != complete[0] for path in family),
            "Expected one complete index format and no mixed index files")
    info = subprocess.run([inspector, "-s", str(prefix)], capture_output=True, text=True, check=True)
    (args.output / "bowtie2_inspect_summary.txt").write_text(info.stdout)
    contigs = {}
    for line in info.stdout.splitlines():
        fields = line.split("\t")
        if fields[0].startswith("Sequence-"):
            require(len(fields) == 3 and fields[1].split(), "Malformed index sequence summary")
            name = fields[1].split()[0]
            require(name not in contigs, "Repeated canonical index reference name")
            contigs[name] = int(fields[2])
    require(contigs == expected_lengths, "Index name/length metadata differs from full genome")
    with (args.output / "bowtie2_inspect_reference.err").open("x") as log:
        child = subprocess.Popen([inspector, "-a", "100000", str(prefix)], stdout=subprocess.PIPE, stderr=log)
        try:
            reconstructed = canonical_sequences(child.stdout)
        except BaseException:
            child.terminate()
            child.wait()
            raise
        finally:
            child.stdout.close()
        require(child.wait() == 0, "Indexed genome reconstruction failed")
    require(reconstructed == expected_sequences, "Indexed canonical genome sequence differs from declared input")
    for path, expected in GUARDS.items():
        require(sha(path) == expected, "Frozen genome changed during indexing")
    require(output_bytes(args.output) <= CAP, "Completed genome index exceeds own cap")
    result = dict(index_built=True, reconstructed_canonical_sequence_agrees=True,
                  contigs=194, bases=3099750718, index_prefix=str(prefix),
                  index_files={str(path): dict(bytes=path.stat().st_size, sha256=sha(path)) for path in complete[0]},
                  protocol_sha256=sha(args.output / "protocol.json"),
                  saved_output_bytes=output_bytes(args.output), receiving_data_read=False,
                  alignment_or_counting_run=False, predictive_validation=False)
    (args.output / "summary.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result), flush=True)


if __name__ == "__main__":
    main()
