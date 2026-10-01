"""Public reference only: declared Ensembl98 full-primary-genome-decoy k31 index."""
import argparse
from collections import defaultdict
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import resource
import signal
import struct
import subprocess
import sys
import time
import urllib.request

ROOT = Path(__file__).resolve().parents[4]
REC = ROOT / "Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z"
TX = REC / "ensembl98_transcript_reconstruction_21998470"
GENOME = REC / "ensembl98_primary_genome_21998416"
FA = TX / "ensembl98_gtf_complete_transcripts.fa"
DNA = GENOME / "Homo_sapiens.GRCh38.dna.primary_assembly.fa"
AXIS = ROOT / "Analysis/MASLD_Model_Benchmark/executions/model-data-064-21079902/fixture/molecular/rna_feature_axis.tsv"
SALMON = Path("/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/salmon")
GUARDS = {
    TX / "summary.json": "1e7d48b28c1178523df3bbc455a6b9cf2c7786124878bf4743566aeb9469e144",
    TX / "transcript_gene_axis.tsv": "0ac3dd2c4cc9b6923afea3041f5d4e94b67b9ded0f33acc4c028e6f3845eb8d2",
    GENOME / "summary.json": "9f93a8291cf676bbf1bd1bc552d2097b3e0f9850a05f86431e1d93df65392f6b",
    GENOME / "acquisition_receipt.json": "72c4190cc3fbd521b9ae81326f50c7f3a517f2251a6d0844c9704ebad74c17a7",
    Path(str(DNA) + ".fai"): "57f6ed6f8b07437708fff09814489cfda4b976652afabe2c9b1d54463bcba0ad",
    AXIS: "baf983ad18a893b7e5db04364818073facc68b47ae25f63070373481e3b217dd",
    SALMON: "9d45524e437fa95de9ccb06e00a2570487ee398895d3a6218883e240ffbf6b6f",
}
SEQUENCE_GUARDS = {
    FA: "cf3fac2b075b649c5bf8c80ef2ec860f068f63d92357560fbcb3f1c3d12d38e2",
    DNA: "78777b0886e8dfa5e14e4957fbbaa53736fcbaa5668d59e09b6b7945fca93d8c",
}
OFFICIAL_SOURCES = {
    "salmon/v1.10.3/scripts/fetchPufferfish.sh": "1b61e06c996bd3c0086ad2094f07ea0b34e81b6f85a667be2715dd965ae8c5c3",
    "pufferfish/salmon-v1.10.3/src/FixFasta.cpp": "8b7d9195f519bf167e769e13cef6d58c1fd4ccd2bf774aba1332816f31ee836b",
    "pufferfish/salmon-v1.10.3/src/PufferfishIndexer.cpp": "963c8e8736b0f503ea35ad6376f3b604bf79d5250fa9c322029eee2da7f6ffa6",
    "pufferfish/salmon-v1.10.3/src/PufferfishBinaryGFAReader.cpp": "d5c59321f3588c50a87960b58283b8f4b719097171d6798707f4f751f0c266f3",
    "pufferfish/salmon-v1.10.3/include/cereal/types/string.hpp": "473ff82531599b7df54c368a884ec1b388262cdb59b68e3325bfd19a64a72e7f",
    "pufferfish/salmon-v1.10.3/include/cereal/types/vector.hpp": "2d4de56deff7c5f4c17e5a95239fdca389e95a53f524b5da82e366181239601d",
}
K = 31
BYTE_CAP = 100_000_000_000
STOP_BYTES = 90_000_000_000
OWN_OUTPUT = None


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def native_gene(value):
    return re.sub(r"\.\d+$", "", re.sub(r"_PAR_Y$", "", value.upper()))


def rows(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def table(path, records):
    require(bool(records), "Refuse empty qualification table")
    with path.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(records)


def fasta_metadata(path):
    """Stream record digests; at most one FASTA line is resident (fixed genome lines are large)."""
    result = {}
    name = None
    digest = None
    length = non_acgt = 0
    with path.open("rb") as handle:
        for line in handle:
            if line.startswith(b">"):
                if name is not None:
                    result[name] = dict(length=length, sha256=digest.hexdigest(), non_acgt=non_acgt)
                name = line[1:].split()[0].decode("ascii")
                require(name not in result, "Duplicate FASTA name: " + name)
                digest = hashlib.sha256()
                length = non_acgt = 0
            else:
                require(name is not None, "Sequence before FASTA header")
                sequence = line.strip()
                require(not any(x in sequence for x in (b" ", b"\t")), "Whitespace inside sequence")
                digest.update(sequence)
                length += len(sequence)
                non_acgt += len(sequence.translate(None, b"ACGTacgt"))
        if name is not None:
            result[name] = dict(length=length, sha256=digest.hexdigest(), non_acgt=non_acgt)
    require(result and all(r["length"] > 0 for r in result.values()), "Empty FASTA record")
    return result


def tree_bytes(path):
    total = 0
    for parent, _, files in os.walk(path):
        for name in files:
            try:
                total += (Path(parent) / name).stat().st_size
            except FileNotFoundError:
                pass  # Salmon removes its own temporary files during construction.
    return total


def index_reference(command, output):
    """Keep all products; stop the owned process group at the 90GB sampled threshold."""
    peak = 0
    with (output / "salmon.stdout").open("xb") as out, (output / "salmon.stderr").open("xb") as err:
        process = subprocess.Popen(command, stdout=out, stderr=err, start_new_session=True)
        try:
            while process.poll() is None:
                peak = max(peak, tree_bytes(output))
                if peak >= STOP_BYTES:
                    raise RuntimeError("Owned reference/index directory reached 90GB stop threshold")
                time.sleep(1)
            require(process.returncode == 0, "Salmon index failed; exit " + str(process.returncode))
        finally:
            if process.poll() is None:
                os.killpg(process.pid, signal.SIGTERM)
                try:
                    process.wait(timeout=15)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGKILL)
                    process.wait()
    peak = max(peak, tree_bytes(output))
    require(peak <= BYTE_CAP, "Reference/index byte allowance exceeded")
    return peak


def cereal_names_lengths(index):
    """Pinned Cereal native-endian vector<string> prefix and vector<uint32_t> lengths."""
    require(sys.byteorder == "little", "Binary metadata reader qualified only for little endian")

    def u64(handle):
        data = handle.read(8)
        require(len(data) == 8, "Truncated Cereal size")
        return struct.unpack("<Q", data)[0]

    with (index / "ctable.bin").open("rb") as handle:
        count = u64(handle)
        require(227368 <= count <= 227562, "Unexpected retained reference-name count")
        names = []
        for _ in range(count):
            size = u64(handle)
            require(0 < size <= 1024, "Invalid Cereal reference-name length")
            name = handle.read(size)
            require(len(name) == size, "Truncated reference name")
            names.append(name.decode("ascii"))
    require(len(set(names)) == len(names), "Repeated retained reference name")
    with (index / "reflengths.bin").open("rb") as handle:
        require(u64(handle) == len(names), "Name/length metadata count mismatch")
        data = handle.read(4 * len(names))
        require(len(data) == 4 * len(names) and handle.read(1) == b"", "Invalid reference lengths")
    return dict(zip(names, struct.unpack("<" + "I" * len(names), data)))


def main():
    global OWN_OUTPUT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    require(os.environ.get("SLURM_JOB_ID"), "Public FASTA processing/indexing must run on compute")
    require(args.output.resolve().parent == REC.resolve(), "Output outside owned execution root")
    require(int(os.environ.get("SLURM_CPUS_PER_TASK", "0")) >= 4, "Four allocated CPUs required")
    for path, expected in GUARDS.items():
        require(sha(path) == expected, "Frozen metadata or executable changed: " + str(path))
    require(json.loads((TX / "summary.json").read_text())["reference_identity_qualified"], "Unqualified reconstructed reference")
    txs = rows(TX / "transcript_gene_axis.tsv")
    modeled = rows(AXIS)
    require(len(txs) == 227368 and len(modeled) == 42163, "Frozen axes changed")
    tx_by_name = {row["transcript_id"]: row for row in txs}
    require(len(tx_by_name) == len(txs), "Duplicate transcript IDs")
    contigs = {}
    with Path(str(DNA) + ".fai").open() as handle:
        for line in handle:
            fields = line.split("\t")
            require(fields[0] not in contigs, "Duplicate genome contig")
            contigs[fields[0]] = int(fields[1])
    require(len(contigs) == 194 and not (set(contigs) & set(tx_by_name)), "Contig/transcript name collision or incomplete decoys")
    args.output.mkdir(exist_ok=False)
    OWN_OUTPUT = args.output
    (args.output / "executed_source.py").write_bytes(Path(__file__).read_bytes())
    metadata = args.output / "source_metadata"
    metadata.mkdir()
    for path, expected in OFFICIAL_SOURCES.items():
        url = "https://raw.githubusercontent.com/COMBINE-lab/" + path
        with urllib.request.urlopen(url, timeout=60) as response:
            data = response.read(1_000_001)
        require(len(data) <= 1_000_000 and hashlib.sha256(data).hexdigest() == expected, "Pinned official source changed: " + url)
        (metadata / Path(path).name).write_bytes(data)
    version = subprocess.run([str(SALMON), "--version"], capture_output=True, text=True, check=True)
    help_text = subprocess.run([str(SALMON), "index", "--help"], capture_output=True, text=True, check=True)
    require(version.stdout.strip() == "salmon 1.10.3", "Unexpected Salmon version")
    for flag in ("--keepDuplicates", "--keepFixedFasta", "--no-clip", "--decoys", "--tmpdir"):
        require(flag in help_text.stdout + help_text.stderr, "Installed CLI missing " + flag)
    (metadata / "salmon_version.txt").write_text(version.stdout + version.stderr)
    (metadata / "salmon_index_help.txt").write_text(help_text.stdout + help_text.stderr)
    for path, expected in SEQUENCE_GUARDS.items():
        require(sha(path) == expected, "Public FASTA changed: " + str(path))
    tx_records, dna_records = fasta_metadata(FA), fasta_metadata(DNA)
    require(set(tx_records) == set(tx_by_name), "Transcript FASTA/axis names differ")
    for name, row in tx_by_name.items():
        require(tx_records[name]["length"] == int(row["length"]) and tx_records[name]["sha256"] == row["sequence_sha256"], "Transcript sequence/axis mismatch: " + name)
    require(set(dna_records) == set(contigs), "FAI/genome contig names differ")
    require(all(dna_records[name]["length"] == length for name, length in contigs.items()), "FAI/genome lengths differ")
    gentrome = args.output / "gentrome.fa"
    digest = hashlib.sha256()
    with gentrome.open("xb") as target:
        for path in (FA, DNA):
            last = None
            with path.open("rb") as source:
                for chunk in iter(lambda: source.read(1 << 20), b""):
                    target.write(chunk)
                    digest.update(chunk)
                    last = chunk[-1:]
            require(last == b"\n", "Input FASTA lacks final newline")
    decoys = args.output / "decoys.txt"
    decoys.write_text("\n".join(contigs) + "\n")
    index = args.output / "index"
    command = [str(SALMON), "index", "-t", str(gentrome), "-d", str(decoys), "-i", str(index),
               "--keepDuplicates", "--keepFixedFasta", "--no-clip", "-k", "31", "-p", "4",
               "--tmpdir", str(args.output / "salmon_temporary")]
    receipt = dict(command=command, job_id=os.environ["SLURM_JOB_ID"], python=sys.version,
                   guards={str(p):v for p,v in {**GUARDS, **SEQUENCE_GUARDS}.items()},
                   official_source_sha256=OFFICIAL_SOURCES, gentrome_sha256=digest.hexdigest(),
                   decoys_sha256=sha(decoys), requested_decoys=len(contigs),
                   script_sha256=sha(Path(__file__)), byte_allowance=BYTE_CAP,
                   launcher_sha256=sha(Path(__file__).with_name("run_qualify_ensembl98_salmon_index.sbatch")),
                   sampled_stop_threshold_bytes=STOP_BYTES, sampling_seconds=1,
                   byte_limit_is_sampled_not_filesystem_quota=True,
                   receiving_sequence_or_count_values_read=False, alignment_or_quantification_run=False,
                   private_pisces_equivalence_asserted=False, complete_observed_input_admission=False)
    (args.output / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    resource.setrlimit(resource.RLIMIT_FSIZE, (BYTE_CAP, BYTE_CAP))
    peak = index_reference(command, args.output)
    fixed = fasta_metadata(index / "ref_k31_fixed.fa")
    retained = cereal_names_lengths(index)
    info = json.loads((index / "info.json").read_text())
    require(info["k"] == K and info["keep_duplicates"] is True, "Index flags disagree with recipe")
    require(set(retained) == set(tx_records) | set(dna_records), "Retained names omit requested transcript/decoy")
    require(info["num_decoys"] == 194, "Salmon did not retain all 194 distinct decoys")
    expected_fixed = {name for name, row in {**tx_records, **dna_records}.items() if row["length"] > K}
    require(set(fixed) == expected_fixed, "Processed FASTA differs from declared unclipped k31 eligibility")
    for name, row in {**tx_records, **dna_records}.items():
        require(retained[name] == row["length"], "Retained reference length changed: " + name)
        if name in fixed:
            require(fixed[name]["length"] == row["length"], "Processed reference length changed: " + name)
            if row["non_acgt"] == 0:
                require(fixed[name]["sha256"] == row["sha256"], "ACGT-only reference changed: " + name)
    groups, processed_groups, by_gene = defaultdict(list), defaultdict(list), defaultdict(list)
    for name, row in tx_by_name.items():
        by_gene[native_gene(row["raw_gene_id"])].append(name)
        groups[(tx_records[name]["length"], tx_records[name]["sha256"])].append(name)
        if name in fixed:
            processed_groups[(fixed[name]["length"], fixed[name]["sha256"])].append(name)
    cross_gene = set()
    for names in groups.values():
        if len({tx_by_name[n]["raw_gene_id"] for n in names}) > 1:
            cross_gene.update(names)
    processed_cross_gene = set()
    for names in processed_groups.values():
        if len({tx_by_name[n]["raw_gene_id"] for n in names}) > 1:
            processed_cross_gene.update(names)
    qualification = []
    for row in modeled:
        names = by_gene[native_gene(row["stable_gene_id"])]
        require(names, "Modeled gene lacks annotated transcript")
        qualification.append(dict(**row, native_transcripts=len(names), retained_named_transcripts=sum(n in retained for n in names),
                                  searchable_fixed_transcripts=sum(n in fixed for n in names), structural_short_transcripts=sum(tx_records[n]["length"] <= K for n in names),
                                  cross_gene_identical_sequence_transcripts=sum(n in cross_gene for n in names),
                                  cross_gene_processed_identical_sequence_transcripts=sum(n in processed_cross_gene for n in names),
                                  non_acgt_input_transcripts=sum(tx_records[n]["non_acgt"] > 0 for n in names),
                                  biological_observability_established=False, observed_zero_created=False))
    table(args.output / "modeled_gene_reference_qualification.tsv", qualification)
    table(args.output / "short_reference_presence.tsv", [dict(transcript_id=n, raw_gene_id=tx_by_name[n]["raw_gene_id"], length=r["length"],
          retained_name=n in retained, searchable_fixed_fasta=n in fixed, biological_zero_asserted=False) for n,r in tx_records.items() if r["length"] <= K])
    inventory = [dict(path=str(p.relative_to(args.output)), bytes=p.stat().st_size, sha256=sha(p))
                 for p in sorted(args.output.rglob("*")) if p.is_file()]
    table(args.output / "file_inventory.tsv", inventory)
    summary = dict(reference_index_identity_qualified=True, transcripts=len(tx_records), requested_decoys=194,
                   retained_named_references=len(retained), searchable_references=len(fixed),
                   short_named_transcripts=sum(r["length"] <= K for r in tx_records.values()),
                   modeled_genes=len(qualification), modeled_genes_without_searchable_transcript=sum(r["searchable_fixed_transcripts"] == 0 for r in qualification),
                   modeled_genes_with_cross_gene_identical_sequence=sum(r["cross_gene_identical_sequence_transcripts"] > 0 for r in qualification),
                   modeled_genes_with_cross_gene_processed_identical_sequence=sum(r["cross_gene_processed_identical_sequence_transcripts"] > 0 for r in qualification),
                   sampled_peak_saved_temporary_bytes=peak, final_directory_bytes=tree_bytes(args.output),
                   observed_input_admission=False, biological_observability_established=False,
                   private_pisces_equivalence_asserted=False, receiving_measurement_unresolved=True,
                   receipt_sha256=sha(args.output / "receipt.json"), inventory_sha256=sha(args.output / "file_inventory.tsv"))
    summary_text = json.dumps(summary, indent=2) + "\n"
    require(tree_bytes(args.output) + len(summary_text.encode()) <= BYTE_CAP,
            "Final reference/index products exceed byte allowance")
    (args.output / "summary.json").write_text(summary_text)
    print(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        # Only this invocation's freshly created directory may receive failure metadata.
        if OWN_OUTPUT is not None and not (OWN_OUTPUT / "failure.json").exists():
            (OWN_OUTPUT / "failure.json").write_text(json.dumps(dict(
                reference_index_identity_qualified=False, error_type=type(error).__name__, error=str(error),
                receiving_sequence_or_count_values_read=False, observed_input_admission=False), indent=2) + "\n")
        raise
