"""Unclipped k31 length and sequence-hash duplicate census; metadata only."""
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[4]
REC = ROOT / "Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z"
SOURCE = REC / "ensembl98_transcript_reconstruction_21998470"
AXIS = ROOT / "Analysis/MASLD_Model_Benchmark/executions/model-data-064-21079902/fixture/molecular/rna_feature_axis.tsv"
GUARDS = {SOURCE/"summary.json": "1e7d48b28c1178523df3bbc455a6b9cf2c7786124878bf4743566aeb9469e144",
          SOURCE/"transcript_gene_axis.tsv": "0ac3dd2c4cc9b6923afea3041f5d4e94b67b9ded0f33acc4c028e6f3845eb8d2",
          AXIS: "baf983ad18a893b7e5db04364818073facc68b47ae25f63070373481e3b217dd"}
K = 31


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for b in iter(lambda: handle.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def normalize(x):
    return re.sub(r"\.\d+$", "", re.sub(r"_PAR_Y$", "", x.upper()))


def write(path, records, fields):
    with path.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t"); writer.writeheader(); writer.writerows(records)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    require(os.environ.get("SLURM_JOB_ID"), "Metadata census must run on compute")
    for p, expected in GUARDS.items():
        require(sha(p) == expected, "Frozen reference metadata changed: "+str(p))
    require(json.loads((SOURCE/"summary.json").read_text())["reference_identity_qualified"], "Reference reconstruction not qualified")
    transcripts, groups, by_gene = [], defaultdict(list), defaultdict(list)
    ids = set()
    with (SOURCE/"transcript_gene_axis.tsv").open(newline="") as handle:
        for r in csv.DictReader(handle, delimiter="\t"):
            tx = r["transcript_id"]
            require(tx not in ids, "Duplicate versioned transcript ID")
            ids.add(tx)
            compact = dict(transcript_id=tx, gene_id=r["gene_id"], raw_gene_id=r["raw_gene_id"],
                           raw_transcript_id=r["raw_transcript_id"], length=int(r["length"]), sequence_sha256=r["sequence_sha256"])
            require(compact["length"] > 0 and re.fullmatch(r"[0-9a-f]{64}", compact["sequence_sha256"]), "Invalid transcript metadata")
            transcripts.append(compact)
            groups[(compact["sequence_sha256"], compact["length"])].append(compact)
            by_gene[normalize(r["raw_gene_id"])].append(compact)
    require(len(transcripts) == 227368, "Complete transcript count changed")
    with AXIS.open(newline="") as handle:
        modeled = list(csv.DictReader(handle, delimiter="\t"))
    require(len(modeled) == 42163, "Modeled gene count changed")
    dup_rows, cross_genes, within_genes = [], set(), set()
    for (digest, length), members in sorted(groups.items()):
        if len(members) <= 1:
            continue
        native_genes = {r["raw_gene_id"] for r in members}
        classification = "cross_native_genes" if len(native_genes) > 1 else "within_native_gene"
        (cross_genes if len(native_genes) > 1 else within_genes).update(normalize(g) for g in native_genes)
        dup_rows.append(dict(sequence_sha256=digest, length=length, eligible_unclipped_k31=length>K,
                             classification=classification, transcripts=len(members), native_genes=len(native_genes),
                             transcript_ids=";".join(r["transcript_id"] for r in members),
                             versioned_gene_ids=";".join(sorted({r["gene_id"] for r in members})),
                             raw_gene_ids=";".join(sorted(native_genes))))
    coverage = []
    for r in modeled:
        gene = normalize(r["stable_gene_id"]); source_gene = normalize(r["matrix_gene_id"])
        txs = by_gene.get(gene, [])
        coverage.append(dict(**r, normalized_modeled_gene_id=gene, normalized_source_gene_id=source_gene,
                             transcripts=len(txs), eligible_unclipped_k31=sum(t["length"]>K for t in txs),
                             short_transcripts=sum(t["length"]<=K for t in txs),
                             has_cross_gene_identical_sequence=gene in cross_genes,
                             has_within_gene_identical_sequence=gene in within_genes))
    args.output.mkdir(parents=True, exist_ok=False)
    write(args.output/"modeled_gene_k31_coverage.tsv", coverage, list(coverage[0]))
    write(args.output/"modeled_genes_no_unclipped_k31_transcript.tsv", [r for r in coverage if not r["eligible_unclipped_k31"]], list(coverage[0]))
    write(args.output/"modeled_genes_with_short_or_duplicate_transcripts.tsv", [r for r in coverage if r["short_transcripts"] or r["has_cross_gene_identical_sequence"] or r["has_within_gene_identical_sequence"]], list(coverage[0]))
    dup_fields = ["sequence_sha256", "length", "eligible_unclipped_k31", "classification", "transcripts", "native_genes", "transcript_ids", "versioned_gene_ids", "raw_gene_ids"]
    write(args.output/"sequence_identical_alias_groups.tsv", dup_rows, dup_fields)
    write(args.output/"cross_gene_sequence_identical_alias_groups.tsv", [r for r in dup_rows if r["classification"]=="cross_native_genes"], dup_fields)
    write(args.output/"short_transcripts_k31.tsv", [r for r in transcripts if r["length"]<=K], list(transcripts[0]))
    write(args.output/"tx2gene_versioned.tsv", [dict(TXNAME=r["transcript_id"], GENEID=r["gene_id"]) for r in transcripts], ["TXNAME", "GENEID"])
    write(args.output/"tx2gene_native.tsv", [dict(TXNAME=r["transcript_id"], GENEID=r["raw_gene_id"]) for r in transcripts], ["TXNAME", "GENEID"])
    write(args.output/"transcript_gene_version_mapping.tsv", transcripts, list(transcripts[0]))
    summary = dict(k=K, eligibility_rule="unmodified transcript length >31; no clipping/non-ACGT cleaning simulated",
                   transcripts=len(transcripts), short_transcripts=sum(r["length"]<=K for r in transcripts),
                   modeled_genes=len(coverage), modeled_genes_unclipped_k31_covered=sum(r["eligible_unclipped_k31"]>0 for r in coverage),
                   modeled_genes_with_any_short_transcript=sum(r["short_transcripts"]>0 for r in coverage),
                   duplicate_groups=dict(Counter(r["classification"] for r in dup_rows)),
                   duplicate_transcripts=sum(r["transcripts"] for r in dup_rows),
                   modeled_genes_with_cross_gene_duplicate=sum(r["has_cross_gene_identical_sequence"] for r in coverage),
                   gene_counting_axis_selected=False, transcript_removal_or_zero_fill=False,
                   poly_a_clipping_assessed=False, post_salmon_cleanup_coverage_established=False,
                   sequence_or_count_values_read=False, salmon_index_built=False,
                   input_sha256={str(p):v for p,v in GUARDS.items()}, script_sha256=sha(Path(__file__)),
                   launcher_sha256=sha(Path(__file__).with_name("run_diagnose_salmon_reference_metadata.sbatch")),
                   python=sys.version, slurm_job_id=os.environ["SLURM_JOB_ID"],
                   output_sha256={p.name:sha(p) for p in args.output.iterdir() if p.is_file()})
    (args.output/"summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    print(json.dumps({k:summary[k] for k in ("short_transcripts", "modeled_genes_unclipped_k31_covered", "duplicate_groups", "modeled_genes_with_cross_gene_duplicate")}, sort_keys=True))


if __name__ == "__main__":
    main()
