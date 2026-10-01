"""Fixed B1 RNA-only feasibility calculation, before receiving H3 measurements.

Native Salmon expected counts are an approximate assay input. No model is
fitted, no missing gene row is supplied, and no predictive loss is calculated.
"""
import argparse
import csv
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import signal
import struct
import subprocess
import sys
from urllib.parse import urlsplit
from urllib.request import urlopen

import numpy as np

from codex_build_native_rna_index import output_bytes
from export_training_h3_target_transform import FIX, read_axis, require, sha

ROOT = Path(__file__).resolve().parents[4]
REC = ROOT / "Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z"
ROSTER = REC / "paired_h3k27ac_raw_roster_21997742/raw_files.tsv"
ROSTER_SHA = "f91c00fb37387020a624a02a388ff00d440efc7fa5bdac5ef41da2a8b1ae7e78"
INDEX_DIR = REC / "codex_native_rna_index_21998573"
MAPPING = REC / "ensembl98_transcript_reconstruction_21998470/transcript_gene_axis.tsv"
RELEASE = ROOT / "Analysis/MASLD_Model_Benchmark/release/masld-liver-chromatin-state-v1.2"
WEIGHTS = RELEASE / "weights/chromatin_state_v1_2.npz"
SCORER = RELEASE / "score.py"
ENV = Path("/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin")
ADAPTER = Path(__file__).with_name("codex_native_tximport_gene_counts.R")
DONOR = "B1"
THREADS = 8
CAP = 50_000_000_000
FROZEN = {
    ROSTER: ROSTER_SHA,
    MAPPING: "0ac3dd2c4cc9b6923afea3041f5d4e94b67b9ded0f33acc4c028e6f3845eb8d2",
    WEIGHTS: "72f9eaba359f7c4c61cac2e00da2d5539b5d3ec868dc33b1bacc2ef0105694ec",
    SCORER: "0c9bff99757dcc934529117209397bb7dc9ef60d2293cba1a5b0eb48fadecc39",
    FIX / "h3k27ac_feature_axis.tsv": "cbe35aeb188e531283a1bcd636ba14de4cfbf215dfc544b99e3ebb97f9831986",
    ADAPTER: "8aec882b16a9202348f0232c9f7d4365dd0c237386b2f98cea6cfb0219e86582",
    ENV / "salmon": "9d45524e437fa95de9ccb06e00a2570487ee398895d3a6218883e240ffbf6b6f",
}


def qualify_index_names_lengths():
    """Bounded Pufferfish1.10.3 Cereal metadata read before receiving reads.

    Short names may be reinserted into native metadata after search indexing;
    their presence is not biological observability or a quant.sf guarantee.
    """
    require(sys.byteorder == "little", "Index reader requires little endian")

    def u64(handle):
        data = handle.read(8)
        require(len(data) == 8, "Truncated native index size")
        return struct.unpack("<Q", data)[0]

    names = []
    with (INDEX_DIR / "index/ctable.bin").open("rb") as handle:
        require(u64(handle) == 227562, "Native index reference-name population differs")
        for _ in range(227562):
            size = u64(handle)
            require(0 < size <= 1024, "Invalid native index reference-name length")
            data = handle.read(size)
            require(len(data) == size, "Truncated native index name")
            names.append(data.decode("ascii"))
    require(len(set(names)) == len(names), "Duplicate native index name")
    with (INDEX_DIR / "index/reflengths.bin").open("rb") as handle:
        require(u64(handle) == len(names), "Native index name/length populations differ")
        data = handle.read(4 * len(names))
        require(len(data) == 4 * len(names) and handle.read(1) == b"", "Truncated/extra native lengths")
    actual = dict(zip(names, struct.unpack("<" + "I" * len(names), data)))
    with MAPPING.open(newline="") as handle:
        expected = {row["transcript_id"]: int(row["length"])
                    for row in csv.DictReader(handle, delimiter="\t")}
    require(len(expected) == 227368, "Native mapping population differs")
    genome_fai = REC / "ensembl98_primary_genome_21998416/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai"
    require(sha(genome_fai) == "57f6ed6f8b07437708fff09814489cfda4b976652afabe2c9b1d54463bcba0ad",
            "Fixed decoy length metadata changed")
    decoys = {row.split("\t")[0]: int(row.split("\t")[1])
              for row in genome_fai.read_text().splitlines()}
    require(len(decoys) == 194 and not set(expected).intersection(decoys), "Native decoy population differs")
    require(actual == {**expected, **decoys}, "Native index omits, adds or changes reference identity/length")
    return dict(retained_native_transcripts=227368, retained_decoys=194,
                retained_short_transcripts=sum(length <= 31 for length in expected.values()),
                searchable_or_biological_observability_asserted=False)


def run_checked(command, log, output):
    """Bound saved files; preserve partial outputs on a failure."""
    with log.open("x") as handle:
        child = subprocess.Popen(command, stdout=handle, stderr=subprocess.STDOUT,
                                 start_new_session=True)
        try:
            while child.poll() is None:
                require(output_bytes(output) <= CAP, "Own RNA output cap exceeded")
                try:
                    child.wait(timeout=30)
                except subprocess.TimeoutExpired:
                    pass
            require(child.returncode == 0, "Subprocess failed: " + str(log))
            require(output_bytes(output) <= CAP, "Completed RNA output cap exceeded")
        except BaseException:
            if child.poll() is None:
                os.killpg(child.pid, signal.SIGTERM)
                try:
                    child.wait(timeout=20)
                except subprocess.TimeoutExpired:
                    os.killpg(child.pid, signal.SIGKILL)
                    child.wait()
            raise


def acquire(row, destination):
    url = "https://" + row["fastq_url"]
    parsed = urlsplit(url)
    require(parsed.hostname == "ftp.sra.ebi.ac.uk" and parsed.path.startswith("/vol1/fastq/")
            and not parsed.query and not parsed.fragment, "Unexpected source file URL")
    expected = int(row["fastq_bytes"])
    md5, sha256, size = hashlib.md5(), hashlib.sha256(), 0
    with urlopen(url, timeout=120) as response, destination.open("xb") as dest:
        final = urlsplit(response.geturl())
        require(final.scheme == "https" and final.hostname == parsed.hostname
                and final.path == parsed.path, "Unexpected sequencing-file redirect")
        length = response.headers.get("Content-Length")
        require(length is None or int(length) == expected, "Source Content-Length changed")
        for chunk in iter(lambda: response.read(1 << 20), b""):
            size += len(chunk)
            require(size <= expected, "Sequencing file exceeds frozen source size")
            dest.write(chunk)
            md5.update(chunk)
            sha256.update(chunk)
    require(size == expected and md5.hexdigest() == row["fastq_md5"],
            "Source sequencing file size/MD5 mismatch; partial preserved")
    # Complete decompression checks the gzip CRC/trailer without saved expansion.
    expanded = 0
    with gzip.open(destination, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            expanded += len(chunk)
            require(expanded <= 50 * expected, "Unexpected decompression expansion")
    return dict(**row, file=str(destination), sha256=sha256.hexdigest(),
                verified_bytes=size, gzip_crc_verified=True, expanded_bytes=expanded)


def predict_native(counts_file, output):
    spec = importlib.util.spec_from_file_location("fixed_chromatin_scorer", SCORER)
    scorer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(scorer)
    counts, genes, samples = scorer.read_counts(str(counts_file))
    require(tuple(samples) == (DONOR,), "Native donor column changed")
    collapsed, native_ids, duplicates = scorer.collapse_duplicates(counts, scorer.normalise_ids(genes))
    with np.load(WEIGHTS, allow_pickle=True) as weights:
        axis = np.asarray(weights["gene_axis"]).astype(str)
        baseline = np.asarray(weights["prof_rrr_ym"], float).copy()
    require(len(axis) == len(set(axis)) == 42163, "Frozen modeled RNA axis changed")
    pos = {gene: j for j, gene in enumerate(native_ids)}
    missing = [gene for gene in axis if gene not in pos]
    require(not missing, "Native quantifier lacks modeled gene rows; no padding permitted")
    modeled = collapsed[[pos[gene] for gene in axis]]
    require(modeled.shape == (42163, 1) and np.isfinite(modeled).all()
            and (modeled >= 0).all() and modeled.sum() > 0, "Invalid complete modeled estimates")
    with (output / "modeled_geneCounts.tsv").open("x") as handle:
        handle.write("gene_id\t" + DONOR + "\n")
        for gene, value in zip(axis, modeled[:, 0]):
            handle.write(f"{gene}\t{value:.17g}\n")
    expected_regions = tuple(row["opaque_source_feature_key"]
                             for row in read_axis(FIX / "h3k27ac_feature_axis.tsv"))
    require(len(expected_regions) == 96460 and baseline.shape == (96460,)
            and np.isfinite(baseline).all(), "Frozen region/baseline axis changed")
    results = {}
    for name, form in (("global", "rrr"), ("local", "rrr_offset_cis")):
        result = scorer.profile(modeled, axis, weights_path=str(WEIGHTS),
                                form=form, transport="raw", min_coverage=1.0)
        require(result["axis_coverage"] == 1 and result["n_matched"] == 42163
                and result["profile"].shape == (1, 96460)
                and tuple(result["region_key"]) == expected_regions
                and np.isfinite(result["profile"]).all(), "Frozen prediction input/output differs")
        path = output / (name + "_prediction.npz")
        np.savez_compressed(path, profile=result["profile"], region_key=result["region_key"],
                            sample_id=samples, form=form, transport="raw",
                            training_mean=baseline)
        results[name] = dict(path=str(path), sha256=sha(path), form=form, transport="raw")
    return dict(predictions=results, native_genes=len(genes), normalized_native_genes=len(native_ids),
                duplicate_normalized_rows_summed=duplicates, modeled_genes=42163,
                modeled_zero_estimates=int((modeled == 0).sum()),
                modeled_count_sum=float(modeled.sum()), native_count_sum=float(counts.sum()),
                denominator="all42163 modeled native gene estimated counts, unchanged",
                rows_fabricated=False, biological_absence_asserted=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    require(os.environ.get("SLURM_JOB_ID"), "Receiving RNA work requires compute")
    require(int(os.environ["SLURM_CPUS_PER_TASK"]) >= THREADS, "Insufficient allocated CPUs")
    require(not args.output.exists(), "Refusing output overwrite")
    for path, digest in FROZEN.items():
        require(sha(path) == digest, "Fixed input changed: " + str(path))
    summary = json.loads((INDEX_DIR / "summary.json").read_text())
    protocol = json.loads((INDEX_DIR / "protocol.json").read_text())
    require(summary["index_built"] and summary["index_returncode"] == 0
            and sha(INDEX_DIR / "protocol.json") == summary["protocol_sha256"], "RNA index incomplete")
    require(protocol["transcripts"] == 227368 and protocol["decoys"] == 194
            and protocol["k"] == 31 and protocol["keep_duplicates"]
            and not protocol["poly_a_clipping"], "RNA reference recipe changed")
    for name, digest in summary["index_metadata_sha256"].items():
        require(sha(INDEX_DIR / "index" / name) == digest, "Index metadata changed")
    native_index = qualify_index_names_lengths()
    with ROSTER.open(newline="") as handle:
        rows = [row for row in csv.DictReader(handle, delimiter="\t")
                if row["donor"] == DONOR and row["assay"] == "RNA"]
    rows.sort(key=lambda row: (row["run_accession"], int(row["file_number"])))
    runs = sorted({row["run_accession"] for row in rows})
    require(runs == [f"SRR870243{i}" for i in range(1, 9)] and len(rows) == 16
            and sum(int(row["fastq_bytes"]) for row in rows) == 6261964285,
            "Fixed first-donor RNA roster changed")
    for run in runs:
        pair = [row for row in rows if row["run_accession"] == run]
        require([row["file_number"] for row in pair] == ["1", "2"]
                and len({row["source_biosample"] for row in pair}) == 1,
                "RNA mate join differs")
    version = subprocess.run([str(ENV / "salmon"), "--version"],
                             capture_output=True, text=True, check=True).stdout.strip()
    require(version == "salmon 1.10.3", "Salmon version changed")
    args.output.mkdir(parents=True, exist_ok=False)
    raw = args.output / "raw"
    raw.mkdir()
    paths = {row["fastq_url"]: raw / Path(row["fastq_url"]).name for row in rows}
    require(len(set(paths.values())) == len(rows), "RNA basename collision")
    command = [str(ENV / "salmon"), "quant", "-i", str(INDEX_DIR / "index"), "-l", "A",
               "-1", *[str(paths[row["fastq_url"]]) for row in rows if row["file_number"] == "1"],
               "-2", *[str(paths[row["fastq_url"]]) for row in rows if row["file_number"] == "2"],
               "-p", str(THREADS), "-o", str(args.output / "quant")]
    receipt = dict(donor=DONOR, biological_n=1, selection="first donor in fixed17 roster, not selected by measured RNA/H3",
                   technical_runs=runs, expected_download_bytes=6261964285,
                   raw_cap_and_saved_output_cap_bytes=CAP, quant_command=command,
                   bias_flags=[], quantification_algorithm="installed Salmon1.10.3 defaults, native paired-file pooling, library A before mates",
                   quantification_stochastic_seed="No seed flag in installed quant help; no inferential resampling requested; multithread numerical determinism not asserted",
                   native_input_contract="complete227368 native transcript rows; exact native tximport sums then full42163 modeled axis; structural zeros retained without biological-absence claims",
                   private_PISCES_equivalence=False, receiving_H3_read=False,
                   model_fitted=False, target_or_prediction_recalibrated=False,
                   primary_loss_evaluated=False, independently_validated=False,
                   input_sha256={str(path): digest for path, digest in FROZEN.items()},
                   index_summary_sha256=sha(INDEX_DIR / "summary.json"),
                   native_index_metadata_check=native_index,
                   source_sha256=sha(Path(__file__)),
                   adapter_sha256=sha(ADAPTER),
                   salmon_version=version, python=sys.version, numpy=np.__version__,
                   slurm_job_id=os.environ["SLURM_JOB_ID"])
    (args.output / "protocol.json").write_text(json.dumps(receipt, indent=2) + "\n")
    acquired = []
    for row in rows:
        acquired.append(acquire(row, paths[row["fastq_url"]]))
        require(output_bytes(args.output) <= CAP, "Own RNA saved-output cap exceeded")
        (args.output / "acquisition_receipts.json").write_text(json.dumps(acquired, indent=2) + "\n")
        print(json.dumps(dict(acquired_files=len(acquired), expected_files=16)), flush=True)
    run_checked(command, args.output / "salmon_quant.log", args.output)
    counts = args.output / "geneCounts.tsv"
    adapter_command = [str(ENV / "Rscript"), str(ADAPTER),
                       "--quant", str(args.output / "quant/quant.sf"), "--tx2gene", str(MAPPING),
                       "--sample", DONOR, "--out", str(counts), "--report", str(args.output / "geneCounts.report.json")]
    run_checked(adapter_command, args.output / "tximport.log", args.output)
    prediction = predict_native(counts, args.output)
    for path, digest in FROZEN.items():
        require(sha(path) == digest, "Fixed input changed during RNA calculation")
    require(output_bytes(args.output) <= CAP, "RNA calculation exceeds saved-output cap")
    result = dict(**prediction, donor=DONOR, biological_n=1, acquired_files=16,
                  saved_output_bytes=output_bytes(args.output),
                  protocol_sha256=sha(args.output / "protocol.json"),
                  native_quantification_completed=True, frozen_predictions_completed=True,
                  receiving_H3_read=False, predictive_accuracy_evaluated=False,
                  model_fitted=False, independently_validated=False)
    (args.output / "summary.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result), flush=True)


if __name__ == "__main__":
    main()
